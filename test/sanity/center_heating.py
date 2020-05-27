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
from ComPASS.simulation_context import SimulationContext
import MeshTools as MT
import MeshTools.vtkwriters as vtkw

L = 1  # 2*np.pi
nx, ny, nz = (10,) * 3  # discretization


def u(pts):
    x, y, z = [pts[:, j] for j in range(3)]
    return 0 * np.cos((2 * np.pi / L) * x) * np.sin((2 * np.pi / L) * y)


ComPASS.load_eos("linear_water")
ComPASS.set_output_directory_and_logfile(__file__)
ComPASS.set_gravity(0)
fluid_properties = ComPASS.get_fluid_properties()
fluid_properties.specific_mass = 1.0
fluid_properties.dynamic_viscosity = 1.0
fluid_properties.volumetric_heat_capacity = 1.0
ComPASS.set_rock_volumetric_heat_capacity(1.0)
fluid_properties.thermal_expansivity = 0

grid_info = {
    "shape": (nx, ny, nz),
    "extent": (L, L, L),
    "origin": (-0.5 * L, -0.5 * L, -0.5 * L),
}
# grid = MT.grid3D(
grid = ComPASS.Grid(**grid_info)

# def dump_solution(filename='solution.vtu'):
# mesh = MT.grid3D(**grid_info)
# MT.to_vtu(
# mesh, ComPASS.to_output_directory(filename),
# pointdata={'u': u(mesh.vertices_array())},
# )
# dump_solution()


def cell_heat_source():
    centers = ComPASS.compute_global_cell_centers()
    res = -((2 * np.pi / L) ** 2) * u(centers)
    res[:] = 0
    res[np.linalg.norm(centers, axis=1) < 0.2 * L] = 1
    return res


def set_dirichlet_nodes():
    vertices = ComPASS.global_vertices()
    boundary = ComPASS.get_global_boundary_vertices()
    vtkw.write_vtu(
        vtkw.points_as_vtu_doc(
            vertices, pointdata={"boundary": np.array(boundary, dtype=np.int)},
        ),
        ComPASS.to_output_directory("boundary_vertices"),
    )
    return ComPASS.get_global_boundary_vertices()


ComPASS.init(
    mesh=grid,
    cell_permeability=1.0,
    cell_porosity=0.5,
    cell_thermal_conductivity=1.0,
    # we cannot call directly the functions below because the mesh must be created first
    set_dirichlet_nodes=set_dirichlet_nodes,
    cell_heat_source=cell_heat_source,
)


def set_initial_states(states):
    states.context[:] = 1
    states.p[:] = 0
    states.T[:] = 0
    states.S[:] = 1
    states.C[:] = 1.0


for states in [ComPASS.node_states(), ComPASS.cell_states()]:
    set_initial_states(states)


def set_dirichlet_states():
    states = ComPASS.dirichlet_node_states()
    set_initial_states(states)
    states.T[:] = u(ComPASS.vertices())


set_dirichlet_states()

context = SimulationContext()
context.abort_on_ksp_failure = False
context.dump_system_on_ksp_failure = False
context.abort_on_newton_failure = False

final_time = 10.0
standard_loop(
    fixed_timestep=1.0,
    final_time=final_time,
    output_period=0.1 * final_time,
    # nitermax=1,
    context=context,
)


def dump_solution(filename="solution.vtu"):
    mesh = MT.grid3D(**grid_info)
    MT.to_vtu(
        mesh,
        ComPASS.to_output_directory(filename),
        pointdata={"u": u(ComPASS.vertices())},
        celldata={
            "u": u(ComPASS.compute_cell_centers()),
            "uvol": u(ComPASS.compute_cell_centers())
            * ((2 * np.pi) ** 3 / ComPASS.number_of_cells()),
        },
    )


dump_solution("toto")

for s in [
    "cellthermalsource",
    "nodethermalsource",
    "porovolfouriercell",
    "porovolfouriernode",
]:
    a = getattr(ComPASS, s)()
    print(s, a.min(), a.mean(), a.max())

vertices = ComPASS.vertices()
usol = ComPASS.node_states().T
assert np.allclose(usol, u(vertices))
