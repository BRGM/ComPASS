#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import sys
import copy
from collections import defaultdict

import numpy as np

import ComPASS
import ComPASS.mpi as mpi
from ComPASS.utils.units import *

if len(sys.argv) > 1:
    nb_random_wells = int(sys.argv[1])
else:
    nb_random_wells = 10  # nb_random_wells producers + nb_random_wells injectors

# fmt: off
np.random.seed(12345)  # set the seed to always have the same well pattern
pres = 15.0 * MPa      # initial reservoir pressure
Tres = degC2K(60)      # temperature of the reservoir
k_reservoir = 1e-12    # reservoir permeability in m^2
K_reservoir = 2        # bulk thermal conductivity in W/m/K
phi_reservoir = 0.15   # reservoir porosity
gravity = 0
rw = 0.1               # well radius
Qw = 100.0             # well imposed flowrate
# fmt: on

Lx, Ly, Lz = 1000.0, 1000.0, 100.0
nx, ny, nz = 40, 40, 1  # discretization

simulation = ComPASS.load_eos("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(gravity)

grid = ComPASS.Grid(shape=(nx, ny, nz), extent=(Lx, Ly, Lz),)


def make_producer(x, y, Qw):
    well = simulation.create_vertical_well((x, y), rw)
    well.operate_on_flowrate = Qw, -np.inf
    well.produce()
    return well


def make_injector(x, y, Qw):
    well = simulation.create_vertical_well((x, y), rw)
    well.operate_on_flowrate = -Qw, np.inf
    well.inject(Tres)
    return well


def set_wells(n):
    def make_wells(id_offset, factory):
        xy = np.random.random(2 * n)
        xy.shape = -1, 2
        # We do not want well to be on the grid boundaries
        # so that they do not interfere with dirichlet boundary conditions
        assert nx > 2
        dx = Lx / nx
        xy[:, 0] = xy[:, 0] * (Lx - 2 * dx) + dx
        assert ny > 2
        dy = Ly / ny
        xy[:, 1] = xy[:, 1] * (Ly - 2 * dy) + dy
        wells = []
        for wid, pos in enumerate(xy):
            wells.append(factory(pos[0], pos[1], Qw))
            wells[-1].id = id_offset + wid
        return wells

    return make_wells(0, make_producer) + make_wells(n, make_injector)


simulation.init(
    mesh=grid,
    set_dirichlet_nodes=simulation.vertical_boundaries(grid),
    wells=lambda: set_wells(nb_random_wells),
    cell_porosity=phi_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=K_reservoir,
    display_well_ids=True,
)


X0 = simulation.build_state(simulation.Context.liquid, p=pres, T=Tres)
simulation.all_states().set(X0)
simulation.dirichlet_node_states().set(X0)

# We must define well connections between source wells (producers) and target wells (injectors)
doublets = np.reshape(np.arange(2 * nb_random_wells), (-1, 2), order="F")
all_wells = set(doublets.ravel())
# make all doublet wells available on master proc
proc_requests = [(mpi.master_proc_rank, all_wells)]
simulation.add_well_connections(well_pairs=doublets, proc_requests=proc_requests)

# A fuction giving the temperature delta imposed by heat network at time t
def network_deltaT(t):
    return 20 * (1 + np.sin(t * (2 * np.pi) / year))


def chain_wells(tick):
    for source, target in doublets:
        data = simulation.get_well_data(target)
        if data is not None:
            wellhead = simulation.well_connections[source]
            assert wellhead.molar_flowrate.shape == (1,)  # a single component
            molar_flowrate = wellhead.molar_flowrate[0]
            assert (
                molar_flowrate >= 0
            ), f"source well {source} with flowrate {molar_flowrate} should be a producer"
            data.imposed_flowrate = -molar_flowrate
            data.injection_temperature = max(
                degC2K(20), wellhead.temperature - network_deltaT(tick.time)
            )


iteration_callbacks = [chain_wells]

if mpi.is_on_master_proc:
    # This will be kept in memory for the whole simulation
    # you may need to flush data to a file at some point
    wellhead_states = defaultdict(list)

    def collect_wellhead_states(tick):
        for well in all_wells:
            wellhead_states[well].append(
                [tick.time, *simulation.well_connections[well]]
            )

    iteration_callbacks.append(collect_wellhead_states)  # only on master proc...


simulation.standard_loop(
    initial_timestep=day,
    output_period=10 * day,
    final_time=year,
    iteration_callbacks=iteration_callbacks,
)

for wk in range(2 * nb_random_wells):
    perforations = simulation.get_well_perforations_state(wk)
    print(perforations)
    if perforations is not None:
        print([perf.pressure for perf in perforations])
    print(simulation.get_wellhead(wk))

if mpi.is_on_master_proc:
    for well, transient in wellhead_states.items():
        # Will produce a text file with 5 columns
        # time (in s), mass_florate, energy_flowrate, pressure (Pa), temperature (K)
        np.savetxt(f"well_{well:04d}.txt", np.asarray(transient), fmt="%.5g")
