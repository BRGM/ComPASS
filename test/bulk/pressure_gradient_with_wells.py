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
from ComPASS.timeloops import standard_loop, TimeStepManager

np.random.seed(12345) # Set the seed to have always the same pattern
nb_random_wells = 2

# T_right = degC2K( 20. )         # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
k_reservoir = 1E-12           # column permeability in m^2 (low permeability -> bigger time steps)
K_reservoir = 2                   # bulk thermal conductivity in W/m/K
phi_reservoir = 0.15          # column porosity
p_origin = 15. * MPa              # initial reservoir pressure
gradp = 1 * bar / 500 # m
gravity = 9.81
ztop_reservoir = 1500
Ttop = degC2K(60)  # temperature at the top of the reservoir
vgradT = 0 / km  # degrees per km - constant to see injection effect
geotherm = lambda zeta: Ttop + vgradT * (ztop_reservoir - zeta)
rw = 0.1
Tinjection = degC2K(60)
epsilon = 0.1

Lx, Ly, Lz = 1000., 500., 10.                  # column dimensions
nx, ny, nz = 100, 50, 1  # discretization

pressure_gradient = lambda x, y: p_origin + gradp * x

simulation = ComPASS.load_eos('water2ph')
ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(gravity)
simulation.info.activate_direct_solver = True

def hydrostatic_pressure(zbottom, ztop, nz, nbsteps=100):
    assert zbottom < ztop
    z = np.linspace(zbottom, ztop, nz)[::-1]  # from top to bottom
    rho = simulation.liquid_molar_density
    p = 0
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


grid = ComPASS.Grid(
    shape = (nx, ny, nz),
    extent = (Lx, Ly, Lz),
)

def create_well(x, y):
    return simulation.create_vertical_well((x, y), rw)


def make_producer(x, y, Qw):
    well = create_well(x, y)
    well.operate_on_flowrate = Qw, 0.0 * MPa
    well.produce()
    return well


def make_injector(x, y, Qw):
    well = create_well(x, y)
    well.operate_on_flowrate = Qw, max(pressure_gradient(0,0), pressure_gradient(Lx, Ly)) + 100.0 * MPa
    well.inject(Tinjection)
    return well


def select_dirichlet_nodes():
    vertices = simulation.global_vertices()
    x, y = vertices[:, 0], vertices[:, 1]
    return (x <= epsilon) | (x >= Lx - epsilon) | (y <= epsilon) | (y >= Ly - epsilon)

def set_wells(n):
    toss = np.random.random(3*n)
    welltype = np.asarray(np.round(toss[:n], 0), dtype=np.int)
    # print(welltype)
    xy = np.reshape(toss[n:], (n, 2))
    xy[:, 0]*= Lx
    xy[:, 1]*= Ly
    wells = []
    for wid, (wt, pos) in enumerate(zip(welltype, xy)):
        # print(wt, pos)
        if wt==0:
            wells.append(make_producer(pos[0], pos[1], 0))
        else:
            wells.append(make_injector(pos[0], pos[1], 0))
        wells[-1].id = wid
    return wells


simulation.init(
    mesh=grid,
    set_dirichlet_nodes=select_dirichlet_nodes,
    wells=lambda: set_wells(nb_random_wells),
    cell_porosity=phi_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=K_reservoir,
    display_well_ids = True,
)


# -- Set initial state and boundary conditions
hp = hydrostatic_pressure(0, Lz, 10*nz, 10)
initial_state = simulation.build_state(simulation.Context.liquid, p=p_origin, T=Ttop)
simulation.all_states().set(initial_state)
dirichlet = simulation.dirichlet_node_states()
dirichlet.set(initial_state)  # will init all variables: context, states...

def set_pT_distribution(states, xyz):
    x, y, z = [xyz[:, j] for j in range(3)]
    states.p[:] = pressure_gradient(x, y) + hp(z)
    states.T[:] = geotherm(z)
set_pT_distribution(simulation.node_states(), simulation.vertices())
set_pT_distribution(simulation.cell_states(), simulation.compute_cell_centers())
set_pT_distribution(dirichlet, simulation.vertices())

# Close all wells
for wid in range(nb_random_wells):
    simulation.close_well(wid)

final_time = 500 * year
output_period = 0.1 * final_time
standard_loop(
    simulation,
    final_time = final_time,
    time_step_manager = TimeStepManager(1 * year, output_period),
    output_period = output_period,
)
