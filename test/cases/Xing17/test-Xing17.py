# -*- coding: utf-8 -*-
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
from ComPASS.timestep_management import TimeStepManager
import MeshTools as MT

# ComPASS.load_eos('water2ph')

H = 3 * km
k_fracture = 1e-12  # fracture permeability in m^2
k_matrix = 1e-16  # reservoir permeability in m^2
phi_fracture = 0.5  # fracture porosity
phi_matrix = 0.25  # reservoir porosity
lambda_reservoir = 2  # reservoir thermal conductivity ???
lambda_fracture = 2  # fracture thermal conductivity ???
ptop = 1 * bar
pbot = 30.0 * MPa
Ttop = degC2K(25.0)
Tbot_matrix = degC2K(200.0)
Tbot_fracture = degC2K(325.0)
thfrac = 1.0  # fracture thickness in meters
qmass = 1.0  # surface mass flux (kg/m2)

# no flux first
# bottom flux condion

ComPASS.load_eos("water2ph")
ComPASS.set_fracture_thickness(thfrac)

basename = "47K/input.2"

# extract information from TetGen file
nodefile = basename + ".node"
elementfile = basename + ".ele"
facefile = basename + ".face"
vertices = np.loadtxt(nodefile, skiprows=1, usecols=(1, 2, 3))
tets = np.loadtxt(elementfile, skiprows=1, usecols=(1, 2, 3, 4), dtype=int)
tets -= 1  # start indexing at 0
tets = MT.idarray(tets)
faces = np.loadtxt(facefile, skiprows=1, usecols=(1, 2, 3), dtype=int)
faces -= 1  # start indexing at 0
faces = MT.idarray(faces)

# filter out box boundaries (unit box)
not_boundary = np.ones(faces.shape[0], dtype=bool)
for iaxis in range(3):
    coordi = vertices[:, iaxis][faces]
    for value in (0, 1):
        on_boundary = np.all(coordi == value, axis=1)
        not_boundary[on_boundary] = False
faces = faces[not_boundary]

# scale domain
vertices *= H
mesh = MT.TetMesh.make(vertices, tets)
mesh_faces = mesh.connectivity.faces
# should be available through C++
fracture_faces = [mesh_faces.id(MT.Triangle(triangle)) for triangle in faces]
fracture_faces = np.array(fracture_faces)


def dirichlet_nodes():
    vertices = ComPASS.global_vertices().view(np.double).reshape((-1, 3))
    on_top = vertices[:, 2] >= H
    on_bottom = vertices[:, 2] <= 0
    fracture_nodes = np.rec.array(ComPASS.global_node_info(), copy=False).frac == ord(
        "y"
    )
    return on_top | (on_bottom & (np.logical_not(fracture_nodes)))


ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    mesh=mesh,
    cell_permeability=k_matrix,
    cell_porosity=phi_matrix,
    cell_thermal_conductivity=lambda_reservoir,
    fracture_faces=lambda: fracture_faces,
    fracture_permeability=k_fracture,
    fracture_porosity=phi_fracture,
    fracture_thermal_conductivity=lambda_fracture,
    set_dirichlet_nodes=dirichlet_nodes,
)


def set_states(states, z):
    states.context[:] = 2
    states.p[:] = ptop - ((pbot - ptop) / H) * (z - H)
    states.T[:] = Ttop - ((Tbot_matrix - Ttop) / H) * (z - H)
    states.S[:] = [0, 1]
    states.C[:] = 1.0


set_states(ComPASS.dirichlet_node_states(), ComPASS.vertices()[:, 2])
set_states(ComPASS.node_states(), ComPASS.vertices()[:, 2])
set_states(ComPASS.cell_states(), ComPASS.compute_cell_centers()[:, 2])
set_states(ComPASS.fracture_states(), ComPASS.compute_fracture_centers()[:, 2])


def set_fracture_dirichlet_bottom_temperature():
    bottom_nodes = ComPASS.vertices()[:, 2] == 0
    fracture_nodes = np.rec.array(ComPASS.node_info(), copy=False).frac == ord("y")
    ComPASS.dirichlet_node_states().T[bottom_nodes & fracture_nodes] = Tbot_fracture
    ComPASS.node_states().T[bottom_nodes & fracture_nodes] = Tbot_fracture


set_fracture_dirichlet_bottom_temperature()


def set_boundary_fluxes():
    Neumann = ComPASS.NeumannBC()
    Neumann.molar_flux[:] = qmass
    Neumann.heat_flux = qmass * ComPASS.liquid_molar_enthalpy(pbot, Tbot_fracture)
    face_centers = np.rec.array(ComPASS.face_centers())
    bottom_fracture_edges = ComPASS.find_fracture_edges(
        ComPASS.compute_face_centers()[:, 2] <= 0
    )
    ComPASS.set_Neumann_fracture_edges(bottom_fracture_edges, Neumann)


set_boundary_fluxes()

final_time = 200 * year
output_period = 0.1 * final_time
standard_loop(
    final_time=final_time,
    output_period=output_period,
    time_step_manager=TimeStepManager(
        initial_timestep=1 * minute, maximum_timestep=output_period
    ),
)
