#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

import ComPASS
from ComPASS.utils.units import *

wid = 0  # well id
ds = 40  # horizontal cell size
H = 100  # height of the well
nv = 10  # number of vertical layers
hns = 5  # half the number of cells along basis side
rw = 0.1  # well radius
ptop = (
    10 * MPa
)  # pressure at the top of the reservoir, 10*MPa one phase. 1*MPa for two phases
Ttop = degC2K(130)  # temperature at the top of the reservoir
vgradT = 0 / km  # degrees per km - constant to see injection effect
geotherm = lambda zeta: Ttop + vgradT * (H - zeta)
gravity = 9.81
omega_reservoir = 0.15  # reservoir porosity
k_reservoir = (
    1e-12  # reservoir permeability in m^2 - low value to limit reservoir convection
)
K_reservoir = 2  # bulk thermal conductivity in W/m/K
Qw = 300.0


ComPASS.set_output_directory_and_logfile(__file__)

simulation = ComPASS.load_eos("water2ph")
simulation.set_gravity(gravity)

slim = hns * ds
ns = 2 * hns
grid = ComPASS.Grid(
    shape=(ns, ns, nv), extent=(ns * ds, ns * ds, H), origin=(-slim, -slim, 0),
)


def hydrostatic_pressure(zbottom, ztop, nz):
    assert zbottom < ztop
    nbsteps = 100
    z = np.linspace(zbottom, ztop, nz)[::-1]  # from top to bottom
    rho = simulation.liquid_molar_density
    p = ptop
    pressures = [p]
    for zbot, ztop in zip(z[1:], z[:-1]):
        zeta = ztop
        dz = (ztop - zbot) / nbsteps
        for _ in range(nbsteps):
            p += gravity * rho(p, geotherm(zeta)) * dz
            zeta -= dz
        pressures.append(p)
    pressures = np.asarray(pressures)
    return lambda zeta: np.interp(zeta, z[::-1], pressures[::-1])


def make_producer():
    well = simulation.create_vertical_well((0, 0), rw)
    well.id = wid
    well.operate_on_flowrate = Qw, -np.inf
    well.produce()
    return [well]


# This could be read from a table
# it can be a numpy array..
# it must have two columns: time and flowrate
# well will be closed if flowrate == 0
history = [
    (0, 100),
    (3 * day, 0),  # well is closed
    (5 * day, 300),
]
well_operations = simulation.well_production_history(wid, history)


def select_dirichlet_nodes():
    vertices = simulation.global_vertices()
    x, y = vertices[:, 0], vertices[:, 1]
    return (x <= -slim) | (x >= slim) | (y <= -slim) | (y >= slim)


simulation.init(
    mesh=grid,
    set_dirichlet_nodes=select_dirichlet_nodes,
    wells=lambda: make_producer(),
    cell_porosity=omega_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=K_reservoir,
)


# -- Set initial state and boundary conditions
hp = hydrostatic_pressure(0, H, 2 * nv + 1)
initial_state = simulation.build_state(simulation.Context.liquid, p=ptop, T=Ttop)
simulation.all_states().set(initial_state)
dirichlet = simulation.dirichlet_node_states()
dirichlet.set(initial_state)  # will init all variables: context, states...


def set_pT_distribution(states, z):
    states.p[:] = hp(z)
    states.T[:] = geotherm(z)


set_pT_distribution(simulation.node_states(), simulation.vertices()[:, 2])
set_pT_distribution(simulation.cell_states(), simulation.compute_cell_centers()[:, 2])
set_pT_distribution(dirichlet, simulation.vertices()[:, 2])

simulation.standard_loop(
    initial_timestep=hour,
    final_time=10 * day,
    output_period=0.5 * day,
    events=well_operations,
)
