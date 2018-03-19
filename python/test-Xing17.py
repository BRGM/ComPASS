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
import MeshTools as MT

# ComPASS.load_eos('water2ph')

H = 3 * km
k_fracture = 1E-12             # fracture permeability in m^2
k_matrix = 1E-16               # reservoir permeability in m^2
phi_fracture = 0.5             # fracture porosity
phi_matrix = 0.25              # reservoir porosity
ptop = 1 * bar
pbot = 30. * MPa
Ttop = degC2K( 25. )
Tbot_matrix = degC2K( 200. )
Tbot_fracture = degC2K( 325. )

# no flux first
# bottom flux condion

ComPASS.load_eos('water2ph')

basename = '47K/input.2'

# extract information from TetGen file
nodefile = basename + '.node'
elementfile = basename + '.ele'
facefile = basename + '.face'
vertices = np.loadtxt(nodefile, skiprows=1, usecols=(1, 2, 3))
tets = np.loadtxt(elementfile, skiprows=1, usecols=(1, 2, 3, 4), dtype=int)
tets-= 1 # start indexing at 0
tets = MT.idarray(tets)
faces = np.loadtxt(facefile, skiprows=1, usecols=(1, 2, 3), dtype=int)
faces-= 1 # start indexing at 0
faces = MT.idarray(faces)

# filter out box boundaries (unit box)
not_boundary = np.ones(faces.shape[0], dtype=bool)
for iaxis in range(3):
    coordi = vertices[:, iaxis][faces]
    for value in (0, 1):
        on_boundary = np.all(coordi==value, axis=1)
        not_boundary[on_boundary] = False
faces = faces[not_boundary]
# z = vertices[:, 2][faces]
# faces = faces[np.all(z>0.1, axis=1)]

# scale domain
vertices*= H
mesh = MT.TetMesh.make(vertices, tets)
mesh_faces = mesh.connectivity.faces
# should be available through C++
fracture_faces = [mesh_faces.id(MT.Triangle(triangle)) for triangle in faces]
fracture_faces = np.array(fracture_faces)

def pressure_dirichlet_boundaries():
    vertices = np.rec.array(ComPASS.global_vertices())
    on_top = (vertices.z >= H)
    on_bottom = (vertices.z <= 0)
    fracture_nodes = np.rec.array(ComPASS.global_node_info(), copy=False).frac == ord('y')
    return on_top | (on_bottom & fracture_nodes)
    #return on_top | on_bottom
    #return on_top

def temperature_dirichlet_boundaries():
    vertices = np.rec.array(ComPASS.global_vertices())
    on_top = (vertices.z >= H)
    on_bottom = (vertices.z <= 0)
    return on_top | on_bottom

ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    mesh = mesh,
    cell_permeability = k_matrix,
    cell_porosity = phi_matrix,
    fracture_faces = lambda: fracture_faces,
    fracture_permeability = k_fracture,
    fracture_porosity = phi_fracture,
    set_pressure_dirichlet_nodes = pressure_dirichlet_boundaries,
    set_temperature_dirichlet_nodes = temperature_dirichlet_boundaries,
    #set_dirichlet_nodes = temperature_dirichlet_boundaries,
)

def set_states(states, z):
    states.context[:] = 2
    states.p[:] = ptop - ((pbot - ptop) / H) * (z - H)
    states.T[:] = Ttop - ((Tbot_matrix - Ttop) / H) * (z - H)
    states.S[:] = [0, 1]
    states.C[:] = 1.
set_states(ComPASS.dirichlet_node_states(), np.rec.array(ComPASS.vertices()).z)
set_states(ComPASS.node_states(), np.rec.array(ComPASS.vertices()).z)
set_states(ComPASS.cell_states(), ComPASS.compute_cell_centers()[:,2])
#def set_fracture_states(states, z):
#    states.context[:] = 2
#    states.p[:] = ptop - ((pbot - ptop) / H) * (z - H)
#    states.T[:] = Ttop - ((Tbot_fracture - Ttop) / H) * (z - H)
#    states.S[:] = [0, 1]
#    states.C[:] = 1.
#set_fracture_states(ComPASS.fracture_states(), ComPASS.compute_fracture_centers()[:,2])
set_states(ComPASS.fracture_states(), ComPASS.compute_fracture_centers()[:,2])
def set_fracture_dirichlet_bottom_temperature():
    bottom_nodes = np.rec.array(ComPASS.vertices()).z == 0
    fracture_nodes = np.rec.array(ComPASS.node_info(), copy=False).frac == ord('y')
    ComPASS.dirichlet_node_states().T[bottom_nodes & fracture_nodes] = Tbot_fracture
    ComPASS.node_states().T[bottom_nodes & fracture_nodes] = Tbot_fracture
set_fracture_dirichlet_bottom_temperature()

#shooter = ComPASS.timeloops.Snapshooter()
#shooter.shoot(0)

final_time = 200 * year
output_period = 10 * year
ComPASS.set_maximum_timestep(output_period)
standard_loop(initial_timestep= 30 * day, final_time = final_time, output_period = output_period)
