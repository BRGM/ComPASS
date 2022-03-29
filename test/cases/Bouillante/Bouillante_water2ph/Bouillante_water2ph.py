#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import ComPASS
import numpy as np
import MeshTools as MT
import MeshTools.CGALWrappers as CGAL
import sys
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop


omega_reservoir = 0.10  # reservoir porosity  .1
fracture_porosity = 0.4  # fracture porosity   .4
k_reservoir = 1e-16  # reservoir permeability in m^2
k_frac = 1e-11  # fracture permeability   1E-11
cell_thermal_cond = 2.0  # reservoir thermal conductivity
fracture_thermal_cond = 2.0  # fracture thermal conductivity
p0 = 0.1 * MPa  # top Pressure
T0 = degC2K(20)  # top Temperature
Tbot_frac = degC2K(250)  # input flux Temperature
gravity = 10.0
qbin = 1.0  # 3e-2
bottom_heat_flux = 0.1
geotherm = 0.02  # bottom_heat_flux / cell_thermal_cond # gradient Temperature/m
# je triche avec geoth pour avoir convergence ???

# pb avec identification edge bottom fracture, je triche
Zdepth = -10000.0

x_source, y_source, radius = (
    634.0e3,
    1781.9e3,
    1.0e3,
)  # coordinate and radius of the source center, almost in one fracture
x_top_diphasic, y_top_diphasic, top_diphasic_radius = (
    640.0e3,
    1794.0e3,
    4.0e3,
)  # coordinate and radius of the mountain (xmax, ymax, radius) top diphasic : not possible because only water


bottom_flag = 1
sea_flag = 30
top_gas_flag = 31
top_diphasic_flag = 32

basename = "Bouillante_Test02"

ComPASS.load_eos("water2ph")
ComPASS.set_gravity(gravity)
ComPASS.set_output_directory_and_logfile(__file__)

topography_vertices = np.load(basename + "_topography.npy")
dtm = CGAL.triangulate_points(topography_vertices)

if ComPASS.mpi.is_on_master_proc:

    # load mesh using CGal
    c3t3 = CGAL.C3t3(basename + ".c3t3.binary.cgal")
    mesh = MT.TetMesh.make(
        c3t3.vertices, MT.idarray(c3t3.cells)
    )  # meshtool objects readable by Compass

    def retrieve_named_indices(filename):
        d = {}
        with open(filename) as f:
            for line in f:
                l = line.split()
                if len(l) > 1:
                    d[l[1]] = int(l[0])
        return d

    box_boundaries = retrieve_named_indices(basename + "_box.txt")
    faults = retrieve_named_indices(basename + "_faults.txt")

    facet_tags = c3t3.facet_tags
    # Legacy facet encoding for box boundaries and faults in GeoModeller :-(
    boundary_facets = {
        name: np.where(facet_tags[:, 0] == -val)[0]
        for name, val in box_boundaries.items()
    }
    faults_facets = {
        key: np.where(
            ((facet_tags[:, 0] // 1000) % 10 == 1)
            & ((facet_tags[:, 1] // 1000) % 10 == 2)
            & (facet_tags[:, 0] % 1000 == val)
        )[0]
        for key, val in faults.items()
    }
    toponodes = MT.idarray(c3t3.facets[boundary_facets["zmax"]])
    leftnodes = MT.idarray(c3t3.facets[boundary_facets["xmin"]])
    rightnodes = MT.idarray(c3t3.facets[boundary_facets["xmax"]])
    frontnodes = MT.idarray(c3t3.facets[boundary_facets["ymin"]])
    backnodes = MT.idarray(c3t3.facets[boundary_facets["ymax"]])
    bottomnodes = MT.idarray(c3t3.facets[boundary_facets["zmin"]])

    def P_Dirichlet_node():
        result = np.zeros(mesh.nb_vertices, dtype=bool)
        result[toponodes] = True
        #    result[leftnodes] = True
        #    result[rightnodes] = True
        #    result[frontnodes] = True
        #    result[backnodes] = True
        # fracn = (np.rec.array(ComPASS.global_node_info(), copy=False).frac == ord('y'))
        # for nn in bottomnodes:
        # edge_center = np.array([ComPASS.vertices[nn].mean(axis=0)])
        # for n in nn:
        # if fracn[n]:
        # edge_centers = np.array([vertices[edge].mean(axis=0) for edge in bottom_fracture_edges])
        # result[n] = ((edge_centers[:,0]-x_source)**2+(edge_centers[:,1]-y_source)**2)>radius**2)
        ### result[n] = (((c3t3.vertices[n,0]-x_source)**2+(c3t3.vertices[n,1]-y_source)**2)>radius**2)
        # else:
        # result[n] = True
        return result

    def fracture_faces():
        sys.stdout.write(
            "selecting fractures faces for:" + str(faults_facets.keys()) + "\n"
        )
        result = np.zeros(mesh.nb_faces, dtype=bool)
        for fs in faults_facets.values():
            for f in fs:
                result[mesh.face_id(MT.Triangle(MT.idarray(c3t3.facets[f])))] = True
        return result

    def set_global_flags():
        nodeflags = ComPASS.global_nodeflags()
        nodeflags[:] = 0
        nodeflags[bottomnodes] = bottom_flag
        for (
            nodes
        ) in (
            toponodes
        ):  # inutile: fait 3 fois le test sur chaque noeud car boucle sur les top mailles...
            for nn in nodes:
                xyz_topo = c3t3.vertices[nn]
                if xyz_topo[2] < 0.0:
                    nodeflags[nn] = sea_flag
                elif (xyz_topo[0] - x_top_diphasic) ** 2 + (
                    xyz_topo[1] - y_top_diphasic
                ) ** 2 < top_diphasic_radius**2:
                    nodeflags[nn] = top_diphasic_flag
                else:
                    nodeflags[nn] = top_gas_flag
        faceflags = ComPASS.global_faceflags()
        faceflags[:] = 0
        faceflags[boundary_facets["zmin"]] = bottom_flag


if not ComPASS.mpi.is_on_master_proc:
    mesh = P_Dirichlet_node = fracture_faces = set_global_flags = None

ComPASS.init(
    mesh=mesh,
    fracture_faces=fracture_faces,
    set_dirichlet_nodes=P_Dirichlet_node,
    cell_porosity=omega_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=cell_thermal_cond,
    fracture_porosity=fracture_porosity,
    fracture_permeability=k_frac,
    fracture_thermal_conductivity=fracture_thermal_cond,
    set_global_flags=set_global_flags,
)

sys.stdout.write("Maillage distribue" + "\n")


def to_array(xyz):
    array = np.zeros((xyz.shape[0], 3))
    for i, elt in enumerate(xyz):
        array[i] = np.fromiter(elt, dtype=np.float)
    return array


def lininterp(depths, top, gradient):
    return top + (gradient) * (depths)


def inside_heat_source(pts):
    return (pts[:, 0] - x_source) ** 2 + (pts[:, 1] - y_source) ** 2 < radius**2


def set_fracture_state(state, depths):
    state.context[:] = 2
    state.p[:] = lininterp(depths, p0, 800 * gravity)
    state.T[:] = lininterp(depths, T0, geotherm)
    state.S[:] = [0, 1]
    state.C[:] = 1.0


def set_Dirichlet_state(state, pts):  # top nodes
    node_flags = ComPASS.nodeflags()
    state.context[:] = 2
    # sea
    state.p[node_flags == sea_flag] = lininterp(
        pts[node_flags == sea_flag, 2], p0, -800 * gravity
    )
    Ttemp = []
    Tinterp = lininterp(pts[node_flags == sea_flag, 2], T0, 0.22)
    for Te in Tinterp:
        Ttemp.append(max(Te, degC2K(5)))
    state.T[node_flags == sea_flag] = Ttemp
    state.S[node_flags == sea_flag] = [0, 1]
    state.C[node_flags == sea_flag] = 1.0
    # monophasic gas
    state.p[node_flags == top_gas_flag] = p0
    state.T[node_flags == top_gas_flag] = T0
    state.S[node_flags == top_gas_flag] = [0, 1]  # [1, 0] il n'y a pas encore de gaz
    state.C[node_flags == top_gas_flag] = 1.0
    # diphasic
    state.p[node_flags == top_diphasic_flag] = p0
    state.T[node_flags == top_diphasic_flag] = T0
    state.S[node_flags == top_diphasic_flag] = [
        0,
        1,
    ]  # [.72, .28] il n'y a pas encore de gaz
    state.C[node_flags == top_diphasic_flag] = 1.0


def set_states(state, depths):
    # context 1:Ig; 2:Il; 3:Ig+Il
    state.context[:] = 2
    state.p[:] = lininterp(depths, p0, 800 * gravity)
    state.T[:] = lininterp(depths, T0, geotherm)
    state.S[:] = [0, 1]
    state.C[:] = 1.0


def set_variable_initial_bc_values():
    vertices_depths = dtm.depths(ComPASS.vertices())
    set_states(ComPASS.node_states(), vertices_depths)
    set_states(ComPASS.cell_states(), dtm.depths(ComPASS.compute_cell_centers()))
    set_Dirichlet_state(ComPASS.dirichlet_node_states(), ComPASS.vertices())
    set_fracture_state(
        ComPASS.fracture_states(), dtm.depths(ComPASS.compute_fracture_centers())
    )


# no molar flux but energy flux at the bottom everywhere
# in some part of the fracture, molar flux and energy flux corresponding to Tbot_frac.
def set_variable_boundary_heat_flux():
    Neumann = ComPASS.NeumannBC()
    Neumann.molar_flux[:] = qbin
    Neumann.heat_flux = qbin * ComPASS.liquid_molar_enthalpy(110 * MPa, Tbot_frac)
    face_centers = np.rec.array(ComPASS.face_centers())
    bottom_faces = face_centers[:, 2] <= Zdepth
    bottom_fracture_edges = ComPASS.find_fracture_edges(bottom_faces)
    vertices = ComPASS.vertices()
    edge_centers = np.array(
        [vertices[edge].mean(axis=0) for edge in bottom_fracture_edges]
    )
    neumann_edges = bottom_fracture_edges[inside_heat_source(edge_centers)]
    ComPASS.set_Neumann_fracture_edges(neumann_edges, Neumann)
    Neumann = ComPASS.NeumannBC()
    Neumann.molar_flux[:] = 0
    Neumann.heat_flux = bottom_heat_flux
    ComPASS.set_Neumann_faces(bottom_faces, Neumann)


sys.stdout.write("set initial and BC" + "\n")
set_variable_initial_bc_values()
sys.stdout.write("set Neumann BC" + "\n")
set_variable_boundary_heat_flux()

sys.stdout.flush()

init_dt = 0.15 * hour
final_time = 1000.0 * year
output_period = 0.1 * final_time
ComPASS.set_maximum_timestep(0.7 * year)


standard_loop(
    initial_timestep=init_dt,
    final_time=final_time,
    output_every=20,
    # , output_period = output_period, specific_outputs=[1. * day], output_every=20,
)
