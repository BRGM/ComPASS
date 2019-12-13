#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

from collections import namedtuple
import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop, TimeStepManager
import ComPASS.mpi as mpi


p_right = 1.0 * bar  # initial reservoir pressure
T_right = degC2K(20.0)  # initial reservoir temperature
k_matrix = 1e-12  # column permeability in m^2 (low permeability -> bigger time steps)
K_matrix = 2  # bulk thermal conductivity in W/m/K
phi_matrix = 0.15  # column porosity
mass_flux = 1e-1


Lx, Ly, Lz = 100.0, 1.0, 1.0  # column dimensions
nx, ny, nz = 100, 1, 1  # discretization

# Numbers to be exactly on faces between cells
probe_locations = [(0.021 * Lx, 0, 0), (0.561 * Lx, 0, 0), (0.972 * Lx, 0, 0)]

simulation = ComPASS.load_eos("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(0)  # no gravity

grid = ComPASS.Grid(
    shape=(nx, ny, nz), extent=(Lx, Ly, Lz), origin=(0, -0.5 * Ly, -0.5 * Lz)
)


def left_nodes():
    return simulation.global_vertices()[:, 0] <= 0


def right_nodes():
    return simulation.global_vertices()[:, 0] >= Lx


simulation.init(
    mesh=grid,
    cell_permeability=k_matrix,
    cell_porosity=phi_matrix,
    # set_pressure_dirichlet_nodes = right_nodes,
    # set_temperature_dirichlet_nodes = lambda: left_nodes() | right_nodes(),
    set_dirichlet_nodes=right_nodes,
    cell_thermal_conductivity=K_matrix,
)


def set_initial_states(states):
    states.context[:] = 2
    states.p[:] = p_right
    states.T[:] = T_right
    states.S[:] = [0, 1]
    states.C[:] = 1.0


for states in [
    simulation.dirichlet_node_states(),
    simulation.node_states(),
    simulation.cell_states(),
]:
    set_initial_states(states)


def set_boundary_flux():
    Neumann = ComPASS.NeumannBC()
    Neumann.molar_flux[:] = mass_flux  # one component
    Neumann.heat_flux = mass_flux * simulation.liquid_molar_enthalpy(p_right, T_right)
    face_centers = simulation.face_centers()
    simulation.set_Neumann_faces(face_centers[:, 0] <= 0, Neumann)


set_boundary_flux()


Probe = namedtuple("Probe", ["id", "location", "cell"])


def locate_probes():
    centers = simulation.cell_centers()
    #print('centers on proc', mpi.proc_rank, ':', centers)
    dx = Lx / nx
    probes = []
    nb_cells_own = simulation.nb_cells_own()[mpi.proc_rank]
    for pk, P in enumerate(probe_locations):
        d = np.linalg.norm(centers - np.asarray(P), axis=1)
        where = np.nonzero(d <= 0.5 * dx)[0]
        if where.shape[0] > 0:
            assert where.shape[0] == 1, "setting probe in two locations"
            cell = where[0]
            if cell < nb_cells_own:  # This is important not to store ghost
                probes.append(Probe(pk, P, cell))
    return probes


probes = locate_probes()
probe_cells = np.array([probe.cell for probe in probes])
probe_pressure = []
print(len(probes), "probes on proc", mpi.proc_rank, probe_cells)

def store_data(n, t):
    if len(probe_cells)>0:
        probe_pressure.append((t, simulation.cell_states().p[probe_cells]))
    else:
        probe_pressure.append((t, []))

final_time = 1.2#3600
output_period = 0.1 * final_time
standard_loop(
    simulation,
    final_time=final_time,
    time_step_manager=TimeStepManager(1, output_period),
    output_period=output_period,
    iteration_callbacks=[store_data],
)

# Here each proc could dump the probes that is responsible for
for k, probe in enumerate(probes):
    data = np.array([[data[0], data[1][k]] for data in probe_pressure])
    np.savetxt(f"probe_{probe.id:03d}_pressure", data)

# We can also collect data to write a single file
# We send all data to master proc
comm = mpi.communicator()
data = comm.gather(
    [[probe.id for probe in probes], np.array([data[1] for data in probe_pressure])],
    root=mpi.master_proc_rank,
)
if mpi.is_on_master_proc:
    t = np.array([data[0] for data in probe_pressure])
    # Expand data one item per prob
    probe_data = []
    for proc_data in data:
        for k, pk in enumerate(proc_data[0]):
            probe_data.append((pk, proc_data[1][:, k]))
    probe_data.sort(key=lambda data: data[0])  # sort by probe id
    table = np.vstack([t, np.vstack([data[1] for data in probe_data])])
    np.savetxt("full_data", np.transpose(table))
