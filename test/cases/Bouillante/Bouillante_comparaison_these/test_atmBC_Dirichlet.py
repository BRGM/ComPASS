#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#
# polygone with prisms (only two nodes are not Dirichlet,
# one is homogeneus Neumann, one is atm BC)
# atm BC at the top, Dirichlet each sides
# gravity = 0

import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.linalg.factory import linear_solver
from ComPASS.newton import Newton
from ComPASS.timestep_management import TimeStepManager
import MeshTools as MT
import gmsh_reader

omega_reservoir = 0.35  # reservoir porosity
k_reservoir = 1e-12  # reservoir permeability in m^2, 1D = 10^-12 m^
cell_thermal_cond = 3.0  # reservoir thermal conductivity : no thermal diffusion
Pporous = 1.5 * bar  # porous Pressure (used also to init the freeflow nodes)
Tporous = degC2K(30)  # porous Temperature (used also to init the freeflow nodes)
Ttop = degC2K(25)
CpRoche = 2.0e6
gravity = 0.0

simulation = ComPASS.load_physics("diphasic")
simulation.set_gravity(gravity)
simulation.set_rock_volumetric_heat_capacity(CpRoche)
ComPASS.set_output_directory_and_logfile(__file__)

filename = "polygone.msh"  # 4 cells
nodes, elements = gmsh_reader.retrieve_mesh_elements(filename)
# grep 3d element
cells = [
    elt for elt, _ in elements if type(elt) in (MT.Tetrahedron, MT.Wedge, MT.Hexahedron)
]
mesh = MT.HybridMesh.create(nodes, cells)

eps = 0.1


def select_dirichlet_nodes():
    vertices = simulation.global_vertices()
    d0 = vertices[:, 0] ** 2 + vertices[:, 1] ** 2
    return d0 > 30.0 - eps


simulation.init(
    mesh=mesh,
    set_dirichlet_nodes=select_dirichlet_nodes,
    cell_porosity=omega_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=cell_thermal_cond,
)

fc = simulation.compute_face_centers()
ff_faces = fc[:, -1] >= 30 - eps
simulation.set_freeflow_faces(ff_faces)
is_ff = simulation.get_freeflow_nodes()  # array of bool of size n_nodes

X0 = simulation.build_state(simulation.Context.liquid, p=Pporous, T=Tporous, Cal=0.0)
X_top = simulation.build_state(
    simulation.Context.diphasic_FF_liq_outflow, p=Pporous, T=Ttop, Cal=0.0
)

simulation.all_states().set(X0)
simulation.dirichlet_node_states().set(X0)
simulation.node_states().set(is_ff, X_top)

tsmger = TimeStepManager(
    initial_timestep=0.2 * year,
    minimum_timestep=1,
    maximum_timestep=10.0 * year,
    increase_factor=1.3,
    decrease_factor=0.6,
)

final_time = 40.0 * year

# Construct the linear solver and newton objects outside the time loop
# to set their parameters. Here direct solving is activated
lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-7, 15, lsolver)

simulation.standard_loop(
    final_time=final_time,
    newton=newton,
    time_step_manager=tsmger,
    # output_every=1,
)

# simulation.postprocess()
