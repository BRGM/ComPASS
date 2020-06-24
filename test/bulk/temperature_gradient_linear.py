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


rhof = 1e3  # specific mass in kg/m^3
cpf = 4200  # specific heat in J/kg/K
rhofcpf = rhof * cpf  # volumetric heat capacity
muf = 1e-3  # viscosity Pa.s
p_right = 1.0 * bar  # initial reservoir pressure
T_left = degC2K(30.0)  # temperature of left heat flux
T_right = degC2K(
    20.0
)  # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
k_matrix = 1e-12  # column permeability in m^2 (low permeability -> bigger time steps)
phi_matrix = 0.15  # column porosity
K_matrix = 2  # bulk cell thermal conductivity W/m/K
mass_flux = 1e-1

Lx, Ly, Lz = 100.0, 10.0, 10.0  # column dimensions
nx, ny, nz = 100, 1, 1  # discretization

simulation = ComPASS.load_eos("linear_water")
fluid_properties = simulation.get_fluid_properties()
fluid_properties.specific_mass = rhof
fluid_properties.volumetric_heat_capacity = rhofcpf
fluid_properties.dynamic_viscosity = muf

ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(0)  # no gravity

grid = ComPASS.Grid(
    shape=(nx, ny, nz), extent=(Lx, Ly, Lz), origin=(0, -0.5 * Ly, -0.5 * Lz),
)


def left_nodes():
    return simulation.global_vertices()[:, 0] <= 0


def right_nodes():
    return simulation.global_vertices()[:, 0] >= Lx


def both_ends():
    return left_nodes() | right_nodes()


simulation.init(
    mesh=grid,
    cell_permeability=k_matrix,
    cell_porosity=phi_matrix,
    cell_thermal_conductivity=K_matrix,
    set_dirichlet_nodes=both_ends,
)

# initial values
X0 = simulation.build_state(p=p_right, T=T_right)
simulation.all_states().set(X0)

# dirichlet boundary conditions
simulation.dirichlet_node_states().set(X0)
states = simulation.dirichlet_node_states()
xd = simulation.vertices()[:, 0]
states.T[xd <= 0] = T_left

final_time = 100 * year
output_period = 0.1 * final_time
standard_loop(
    simulation,
    final_time=final_time,
    time_step_manager=TimeStepManager(1e-5, output_period),
    output_period=output_period,
)

if ComPASS.mpi.communicator().size == 1:
    assert ComPASS.mpi.is_on_master_proc
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING - matplotlib was not found - no graphics will be generated")
        plt = None
    else:
        x = simulation.cell_centers()[:, 0]
        states = simulation.cell_states()
        plt.clf()
        plt.plot(x, K2degC(states.T))
        plt.xlabel("x in meters")
        plt.ylabel("temperature")
        plt.savefig(ComPASS.to_output_directory("temperature"))
