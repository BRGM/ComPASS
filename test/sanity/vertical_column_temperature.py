# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import sys
import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.physics.utils import constant_physical_property
from ComPASS.physics.densities import build_pure_phase_volumetric_mass_density
from ComPASS.timeloops import standard_loop

from ComPASS.exceptions import SanityCheckFailure


p0 = 1.0 * bar  # initial reservoir pressure
T0 = degC2K(
    20.0
)  # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
bottom_heat_flux = 0.1  # W/m2
k_matrix = 1e-18  # column permeability in m^2 (low permeability -> bigger time steps)
phi_matrix = 0.15  # column porosity
K_matrix = 2.0  # bulk thermal conductivity in W/m/K

H = 3000.0  # column height
nx, ny, nz = 1, 1, 300  # discretization

rhow = 1e3
rhocp = 2000 * 800  # volumetric heat capacity
muf = 1e-3
g = 10.0

simulation = ComPASS.load_eos("linear_water")
ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(g)
fluid_properties = simulation.get_fluid_properties()
simulation.set_molar_density_functions(
    build_pure_phase_volumetric_mass_density(specific_mass=rhow),
)
simulation.set_viscosity_functions(constant_physical_property(muf))
simulation.set_rock_volumetric_heat_capacity(rhocp)

grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(10.0, 10.0, H),
    origin=(-5, -5, -H),
)


def top_nodes():
    return simulation.global_vertices()[:, 2] >= 0


simulation.init(
    mesh=grid,
    cell_permeability=k_matrix,
    cell_porosity=phi_matrix,
    cell_thermal_conductivity=K_matrix,
    set_dirichlet_nodes=top_nodes,
)


def set_initial_states(states):
    states.context[:] = 1
    states.p[:] = p0
    states.T[:] = T0
    states.S[:] = 1
    states.C[:] = 1.0


for states in [
    simulation.dirichlet_node_states(),
    simulation.node_states(),
    simulation.cell_states(),
]:
    set_initial_states(states)


def set_boundary_heat_flux():
    Neumann = ComPASS.NeumannBC()
    Neumann.heat_flux = bottom_heat_flux
    face_centers = simulation.face_centers()
    simulation.set_Neumann_faces(face_centers[:, 2] <= -H, Neumann)


set_boundary_heat_flux()

final_time = 2e6 * year
# ComPASS.set_maximum_timestep(0.1*final_time)
current_time = standard_loop(
    simulation,
    initial_timestep=1 * day,
    final_time=final_time,
    nitermax=140,
)

T = simulation.cell_states().T
zc = simulation.cell_centers()[:, 2]
# print(np.linalg.norm(T0-zc*bottom_heat_flux/K_matrix-T, np.inf)<1e-2)
# print(ComPASS.get_current_time() < )
# print('Final time:', current_time, current_time/year, current_time >= final_time)
# print('Final error:', np.linalg.norm(T0-zc*bottom_heat_flux/K_matrix-T, np.inf))
if current_time < final_time:
    raise SanityCheckFailure("final time was not reached")

error_Linf = np.linalg.norm(T0 - zc * bottom_heat_flux / K_matrix - T, np.inf)
if error_Linf > 1e-2:
    raise SanityCheckFailure(
        f"mismatch with reference solution ||error||inf = {error_Linf}"
    )
