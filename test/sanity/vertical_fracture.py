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

p0 = 1.0 * bar  # initial reservoir pressure
T0 = degC2K(
    20.0
)  # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
qmass = 1e-1  #
Tbottom = degC2K(50.0)  # temperature influx
k_matrix = 1e-15  # domain permeability in m^2
phi_matrix = 0.15  # domain porosity
k_fracture = 1e-12  # fracture permeability in m^2
phi_fracture = 0.3  # fracture porosity
thermal_cond = 2.0  # bulk thermal conductivity in W/m/K

H = 1000.0  # domain height
nH = 50  # discretization
nx, ny, nz = 2, 1, nH
Lx, Ly, Lz = 2 * H / nH, H / nH, H

simulation = ComPASS.load_physics("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)

# thermodynamic functions are only available once the physics is loaded
pbottom = simulation.get_gravity() * H * 1000.0
hbottom = simulation.liquid_molar_enthalpy(pbottom, Tbottom)

grid = ComPASS.Grid(
    shape=(nx, ny, nz), extent=(Lx, Ly, Lz), origin=(-0.5 * Lx, -0.5 * Ly, -H)
)


def top_nodes():
    return simulation.global_vertices()[:, 2] >= 0


def select_fractures():
    centers = simulation.compute_global_face_centers()
    xc = centers[:, 0]
    zc = centers[:, 2]
    return xc == 0


simulation.init(
    mesh=grid,
    cell_permeability=k_matrix,
    cell_porosity=phi_matrix,
    cell_thermal_conductivity=thermal_cond,
    fracture_faces=select_fractures,
    fracture_permeability=k_fracture,
    fracture_porosity=phi_fracture,
    fracture_thermal_conductivity=thermal_cond,
    set_dirichlet_nodes=top_nodes,
)


def set_initial_states(states):
    states.context[:] = simulation.Context.liquid
    states.p[:] = p0
    states.T[:] = T0
    states.S[:, simulation.phase_index(simulation.Phase.gas)] = 0
    states.S[:, simulation.phase_index(simulation.Phase.liquid)] = 1
    states.C[:] = 1


for states in [
    simulation.dirichlet_node_states(),
    simulation.node_states(),
    simulation.cell_states(),
    simulation.fracture_states(),
]:
    set_initial_states(states)


def set_boundary_fluxes():
    Neumann = ComPASS.NeumannBC()
    Neumann.molar_flux[:] = qmass
    Neumann.heat_flux = qmass * hbottom
    face_centers = simulation.face_centers()
    bottom_fracture_edges = simulation.find_fracture_edges(face_centers[:, 2] <= -H)
    simulation.set_Neumann_fracture_edges(bottom_fracture_edges, Neumann)


set_boundary_fluxes()

# Construct the linear solver and newton objects outside the time loop
# to set their parameters. Here direct solving is activated
lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-5, 8, lsolver)

final_time = 50 * year
output_period = 0.1 * final_time
standard_loop(
    simulation,
    newton=newton,
    initial_timestep=1 * hour,
    final_time=final_time,
    output_period=output_period,
)
