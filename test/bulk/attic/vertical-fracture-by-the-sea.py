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

import ComPASS.GridTools as GT

MT = ComPASS.ComPASS.MeshTools

# import vtkwriters as vtkw

patm = 1e5
Tsea, Tatm = 5.0, 30.0
Tgrad = 0.03  # 30°C / km depth
Lsea = 50.0  # transition length for seafloor temperature
fracture_permeability = 1e-12
matrix_permeability = 1e-15

L, P, H = 3.0e3, 1.0e3, 3.0e3
slope = 200 / L
# nb of cells along Ox and Oz
n = 20

with_fractures = True
with_bottom_dirichlet_conditions = True

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gridshape = (n, 2 * max(1, int(0.5 * P / (L / n))), n)
# gridshape = (n, 1, n)
gridextent = (L, P, H)
# vertices, cells = GT.grid2tets(gridshape, gridextent)
vertices, hexs = GT.grid2hexs(gridshape, gridextent)
topo_nodes = np.nonzero(vertices[:, 2] == H)[0]
bottom_nodes = np.nonzero(vertices[:, 2] == 0)[0]
vertices[:, 0] -= 0.5 * L
vertices[:, 1] -= 0.5 * P
vertices[:, 2] += -H + (vertices[:, 2] / H) * slope * vertices[:, 0]

# tets = np.asarray(tets, dtype=MT.idtype())
# mesh = MT.tetmesh(vertices, tets)
hexs = np.asarray(hexs, dtype=MT.idtype())
mesh = MT.hexmesh(vertices, hexs)
# del tets
del hexs


def compute_depth(pts):
    res = slope * pts[:, 0]
    res -= pts[:, 2]
    res[res < 0] = 0
    return res


face_centers = mesh.all_face_centers()
fault_faces_id = np.nonzero(face_centers[:, 1] == 0)[0]


def set_node_flags():
    flags = ComPASS.global_nodeflags()
    flags[:] = 0
    flags[topo_nodes] = 1
    flags[bottom_nodes] = 2


def select_dirichlet_nodes():
    dirichlet = np.zeros(mesh.nb_vertices(), dtype=bool)
    dirichlet[topo_nodes] = True
    dirichlet[bottom_nodes] = with_bottom_dirichlet_conditions
    return dirichlet


assert ComPASS.gravity() > 0


def select_fractures():
    where = np.zeros(mesh.nb_faces(), dtype=bool)
    where[fault_faces_id] = with_fractures
    return where


# This can be used to set different fractures id (and select properties accordingly)
def set_fracture_flags():
    flags = ComPASS.global_faceflags()
    flags[:] = 0
    flags[fault_faces_id] = 1


def set_global_flags():
    set_node_flags()
    set_fracture_flags()


def set_boundary_conditions():
    vertices = ComPASS.vertices().view(np.double).reshape((-1, 3))
    toponodes = np.nonzero(ComPASS.nodeflags() == 1)[0]
    bottnodes = np.nonzero(ComPASS.nodeflags() == 2)[0]
    # we test where z of topo nodes is < 0
    topoelevation = vertices[toponodes, -1]
    dirichlet = ComPASS.dirichlet_node_states()
    dirichlet.p[:] = patm
    # considering a sea water density of 1030 kg/m3
    p = patm - topoelevation * 1.03e3 * ComPASS.gravity()
    p[p < patm] = patm
    dirichlet.p[toponodes] = p
    botdepth = compute_depth(vertices[bottnodes])
    dirichlet.p[bottnodes] = patm + botdepth * 1.0e3 * ComPASS.gravity()
    # considering 5 Celsius degrees on the sea floor with gradual transition
    dirichlet.T[:] = degC2K(Tatm)
    Ttop = Tatm + (Tatm - Tsea) * topoelevation / Lsea
    Ttop[Ttop < Tsea] = Tsea
    Ttop[Ttop > Tatm] = Tatm
    dirichlet.T[toponodes] = degC2K(Ttop)
    # Tbot = Tatm + (Tatm-Tsea) * vertices[bottnodes, -1]/Lsea
    # Tbot[Tbot<Tsea] = Tsea
    # Tbot[Tbot>Tatm] = Tatm
    # Tbot += Tgrad * botdepth
    dirichlet.T[bottnodes] = degC2K(Tatm + Tgrad * botdepth)
    dirichlet.context[:] = 2
    dirichlet.S[:] = [0, 1]
    dirichlet.C[:] = 1.0


def set_initial_values():
    def set_states(states, pts):
        assert states.shape[0] == pts.shape[0]
        depth = compute_depth(pts)
        states.context[:] = 2
        assert np.all(depth >= 0)
        states.p[:] = patm + depth * 1.0e3 * ComPASS.gravity()
        states.T[:] = degC2K(Tatm + Tgrad * depth)
        states.S[:] = [0, 1]
        states.C[:] = 1.0

    xyz = ComPASS.vertices().view(np.double).reshape((-1, 3))
    set_states(ComPASS.node_states(), xyz)
    xyz = ComPASS.compute_face_centers()[ComPASS.frac_face_id() - 1]  # Fotran indexing
    set_states(ComPASS.fracture_states(), xyz)
    # states = ComPASS.fracture_states()
    # vtkw.write_vtu(
    #    vtkw.vtu_doc(ComPASS.vertices().view(np.double).reshape((-1, 3)),
    #                 # Fotran indexing
    #                 np.array([np.array(a, copy=False) - 1 for a in ComPASS.get_connectivity().NodebyFace])[ComPASS.frac_face_id()-1],
    #                 celltypes=np.tile(vtkw.vtk_celltype['quad'], fracvert.shape[0]),
    #                 celldata={'pressure': np.ascontiguousarray(states.p),
    #                           'temperature': np.ascontiguousarray(states.T),
    #                           'centerdepth': np.ascontiguousarray(compute_depth(xyz)),
    #                           'elevation': np.ascontiguousarray(xyz[:,2])}),
    #    'frac_state_%03d.vtu' % ComPASS.mpi.proc_rank
    # )
    xyz = ComPASS.compute_cell_centers()
    set_states(ComPASS.cell_states(), xyz)
    # states = ComPASS.cell_states()
    # vtkw.write_vtu(
    #    vtkw.vtu_doc(ComPASS.vertices().view(np.double).reshape((-1, 3)),
    #                 # Fotran indexing
    #                 np.array([np.array(a, copy=False) - 1 for a in ComPASS.get_connectivity().NodebyCell]),
    #                 celldata={'pressure': np.ascontiguousarray(states.p),
    #                           'temperature': np.ascontiguousarray(states.T),
    #                           'centerdepth': np.ascontiguousarray(compute_depth(xyz))}),
    #    'cell_state_%03d.vtu' % ComPASS.mpi.proc_rank
    # )


ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    mesh=mesh,
    set_dirichlet_nodes=select_dirichlet_nodes,
    fracture_faces=select_fractures,
    set_global_flags=set_global_flags,
    cells_permeability=lambda: matrix_permeability,
    fractures_permeability=lambda: fracture_permeability,
)

set_initial_values()
set_boundary_conditions()

# no longer useful for computation
del topo_nodes
del bottom_nodes
del vertices
del mesh
del fault_faces_id

ComPASS.set_maximum_timestep(1e2 * year)
standard_loop(
    final_time=1e4 * year, output_period=1e3 * year, initial_timestep=1 * year
)
