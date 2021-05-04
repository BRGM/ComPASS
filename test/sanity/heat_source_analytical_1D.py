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
from ComPASS.timeloops import standard_loop

import MeshTools as MT
import MeshTools.vtkwriters as vtkw

L = 2 * np.pi
n = 50  # discretization
axis = 0


def u(pts):
    return np.sin((2 * np.pi / L) * pts[:, axis])


simulation = ComPASS.load_eos("linear_water")
ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(0)
fluid_properties = simulation.get_fluid_properties()
fluid_properties.specific_mass = 1.0
fluid_properties.dynamic_viscosity = 1.0
fluid_properties.volumetric_heat_capacity = 1.0
simulation.set_rock_volumetric_heat_capacity(1.0)
fluid_properties.thermal_expansivity = 0


shape = [1,] * 3
shape[axis] = n
extent = [1,] * 3
extent[axis] = L
origin = [-0.5,] * 3
origin[axis] = -0.5 * L
grid_info = {"shape": shape, "extent": extent, "origin": origin}
# grid = MT.grid3D(
grid = ComPASS.Grid(**grid_info)


def dump_solution(filename="solution.vtu"):
    mesh = MT.grid3D(**grid_info)
    MT.to_vtu(
        mesh,
        ComPASS.to_output_directory(filename),
        pointdata={"u": u(mesh.vertices_array())},
    )


dump_solution()


def set_dirichlet_nodes():
    x = simulation.global_vertices()[:, axis]
    return (x <= -0.5 * L) | (x >= 0.5 * L)


def cell_heat_source():
    centers = simulation.compute_global_cell_centers()
    return (2 * np.pi / L) ** 2 * u(centers)


simulation.init(
    mesh=grid,
    cell_permeability=1.0,
    cell_porosity=0.5,
    cell_thermal_conductivity=1.0,
    # we cannot call directly the functions below because the mesh must be created first
    set_dirichlet_nodes=set_dirichlet_nodes,
    cell_heat_source=cell_heat_source,
)

X0 = simulation.build_state(p=0, T=0)
simulation.all_states().set(X0)

dirichlet = simulation.dirichlet_node_states()
dirichlet.set(X0)
vertices = simulation.vertices()
dirichlet.T[:] = u(vertices)

final_time = 1e2
standard_loop(
    simulation,
    initial_timestep=1.0,
    final_time=final_time,
    output_period=0.1 * final_time,
    # nitermax=1,
)

usol = simulation.node_states().T
assert np.allclose(usol, u(vertices), atol=1e-3)
