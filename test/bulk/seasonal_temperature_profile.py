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
from ComPASS.timeloops import standard_loop
from ComPASS.linalg.factory import linear_solver
from ComPASS.newton import Newton
from ComPASS.timestep_management import FixedTimeStep
from ComPASS.properties.enthalpies import build_pure_phase_enthalpy


p0 = 0.0 * bar  # dummy pressure no gravity
Tmean = 10.0  # average surface temperature
deltaT = 10.0  # seasonal amplitude
bottom_heat_flux = 0.08  # W/m2
k_matrix = 1e-18  # permeability - not relevant here but cannot be 0
phi_matrix = 0.15  # porosity - not really relevant either as long as
# fluid and rock volumetric heat capacities are the same
K_matrix = 2.0  # bulk thermal conductivity in W/m/K
rhor = 2200.0  # rock density kg/m^3
cpr = 800.0  # rock specific heat capacity J/K/kg

H = 100.0  # column height
nx, ny, nz = 1, 1, 200  # discretization
dz = H / nz

simulation = ComPASS.load_physics("linear_water")
simulation.set_rock_volumetric_heat_capacity(rhor * cpr)
simulation.set_molar_enthalpy_functions(
    build_pure_phase_enthalpy(volumetric_heat_capacity=rhor * cpr),
)
simulation.set_gravity(0)
ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(dz, dz, H),
    origin=(-0.5 * dz, -0.5 * dz, -H),
)


def top_nodes():
    return simulation.global_vertices()[:, 2] >= 0


def set_node_flags():
    simulation.global_nodeflags()[:] = np.asarray(top_nodes(), dtype=np.int)


simulation.init(
    mesh=grid,
    cell_permeability=k_matrix,
    cell_porosity=phi_matrix,
    cell_thermal_conductivity=K_matrix,
    set_dirichlet_nodes=top_nodes,
    set_global_flags=set_node_flags,
)

X0 = simulation.build_state(p=p0, T=Tmean)
simulation.all_states().set(X0)
simulation.dirichlet_node_states().set(X0)
geotherm = lambda z: Tmean - (bottom_heat_flux / K_matrix) * z
simulation.all_states().T[:] = geotherm(simulation.all_positions()[:, 2])
simulation.dirichlet_node_states().T[:] = geotherm(simulation.vertices()[:, 2])


def set_boundary_heat_flux():
    Neumann = ComPASS.NeumannBC()
    Neumann.heat_flux = bottom_heat_flux
    face_centers = simulation.face_centers()
    simulation.set_Neumann_faces(face_centers[:, 2] <= -H, Neumann)


set_boundary_heat_flux()

# locate dirichlet nodes - not mandatory
# we could have identified different regions using nodeflags
dirichlet_nodes = np.nonzero(simulation.nodeflags())[0]
dirichlet_T = simulation.dirichlet_node_states().T


def change_surface_temperature(tick):
    dirichlet_T[dirichlet_nodes] = Tmean + deltaT * np.sin(
        tick.time * (2 * np.pi / year)
    )


# Construct the linear solver and newton objects outside the time loop
# to set their parameters. Here direct solving is activated
lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-5, 8, lsolver)

final_time = 5 * year
output_period = year / 12
standard_loop(
    simulation,
    newton=newton,
    final_time=final_time,
    time_step_manager=FixedTimeStep(5.0 * day),
    output_period=output_period,
    # iteration callbacks are function of the form f(n, t)
    # where n is the iteration number and t is time
    # they are called at the end of each successful iteration
    # you can put as many of them
    iteration_callbacks=[change_surface_temperature],
)
