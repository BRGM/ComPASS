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
from ComPASS.utils.wells import create_vertical_well
from ComPASS.RawMesh import RawMesh
import MeshTools as MT
from MeshTools.io import petrel

pres = 20. * MPa                  # initial reservoir pressure
Tres = degC2K( 70. )              # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
Tinjection = degC2K( 30. )        # injection temperature - convert Celsius to Kelvin degrees
Qm = 300. * ton / hour            # production flowrate
omega_reservoir = 0.15            # reservoir porosity
k_reservoir = 1E-12               # reservoir permeability in m^2
K_reservoir = 2                   # bulk thermal conductivity in W/m/K
g = 9.81

ComPASS.load_eos('water2ph')
ComPASS.set_gravity(g)

def load_eclipse_grid(filename):
    import os
    grdecl, kwargs = petrel.read_file_argument(filename)
    name = os.path.basename(grdecl.split('.')[0])
    pgrid = petrel.PetrelGrid.build_from_files(grdecl, **kwargs)
    hexa, vertices, cell_faces, face_nodes = pgrid.process()
    perm = pgrid.get_perm()
#    mesh = MT.HybridMesh.Mesh()
#    vertices = mesh.vertices
#    for P in vertices:
#        vertices.append(MT.Point(P))
#    cellnodes = mesh.connectivity.cells.nodes
#    for cell in hexa:
#        cellnodes.append(MT.Hexahedron(cell))
#    mesh.connectivity.update_from_cellnodes()
    cell_nodes = [
        np.unique(np.hstack(face_nodes[faces]))
        for faces in cell_faces
    ]
    mesh = RawMesh(
        vertices=vertices,
        cell_faces=list(cell_faces),
        face_nodes=list(face_nodes),
        cell_nodes=cell_nodes,
    )
    return mesh, perm

mesh, perm = load_eclipse_grid('sample.grdecl')

def build_permeability_tensors():
    nbcells = len(perm[0])
    k = np.zeros((nbcells, 3, 3), dtype=np.double)
    for i in range(3): # kx, ky, kz
        k[:, i, i] = np.asarray(perm[i], dtype=np.double)
    k*= 1E-15 # to mdarcy
    return k

def select_dirichlet():
    xyz = ComPASS.global_vertices()
    x, y, z = [xyz[:, i] for i in range(3)]
    return (x==x.min()) | (y==y.min()) | (x==x.max()) | (y==y.max())

def make_wells():
    xyz = ComPASS.global_vertices()
    x, y, z = [xyz[:, i] for i in range(3)]
    Lx, Ly, Lz = (x.max() - x.min(), y.max() - y.min(), z.max() - z.min())
    interwell_distance = 0.5 * Ly
    Cx, Cy = (x.min() + 0.5 * Lx, y.min() + 0.5 * Ly)
    producer = create_vertical_well((Cx, Cy - 0.5 * interwell_distance))
    producer.operate_on_flowrate = Qm , 1. * bar
    producer.produce()
    injector = create_vertical_well((Cx, Cy + 0.5 * interwell_distance))
    injector.operate_on_flowrate = Qm, pres + 100. * MPa
    injector.inject(Tinjection)
    return (producer, injector)

ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    mesh = mesh,
    set_dirichlet_nodes = select_dirichlet,
    wells = make_wells,
    cell_porosity = omega_reservoir,
    cell_permeability = build_permeability_tensors,
    cell_thermal_conductivity = K_reservoir,
)

zref = 3000
pref = 1 * bar # zref = 3000m
rhoref = 1000. # kg / m3

def set_states(states, xyz):
    z = xyz[:, 2]
    states.context[:] = 2
    states.p[:] = pref + rhoref * g * (zref - z)
    states.T[:] = Tres
    states.S[:] = [0, 1]
    states.C[:] = 1.
for states, xyz in [
    (ComPASS.dirichlet_node_states(), ComPASS.vertices()),
    (ComPASS.node_states(), ComPASS.vertices()),
    (ComPASS.cell_states(), ComPASS.compute_cell_centers()),
]:
    set_states(states, xyz)

standard_loop(
    initial_timestep = 1 * day, final_time = 30 * year,
    output_period = year,
)
