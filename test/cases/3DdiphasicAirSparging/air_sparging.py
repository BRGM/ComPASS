# -*- coding: utf-8 -*-
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#
# Air injection into groundwater (air sparging)
# in a homogeneous axially symmetric porous medium.
# Modeled using a two-phase flow approach over 3D mesh.
#
# axially symmetric domain :
#                 Dirichlet p_atm
#          |****************************** z=1.5m gas only
#          |        Unsaturated          *
#          |           zone              *
#          |-----------------------------* z=0m, liquid water below
#          |                             *
#          |                             *
# No-flow  |                             *
#          |                             *
#          |                             * Dirichlet (hydrostatic pressure)
#          ->                            *
#          -> Injection                  *
#          -> Faces                      *
#          ->                            *
#          |                             *
#          |           No-flow           *
#          |****************************** z=-6m


import numpy as np
import MeshTools as MT
import os
from mpi4py import MPI
import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop, TimeStepManager
import vtkwriters as vtkw
from ComPASS.utils.wells import create_vertical_well
import GROUPS as groups
import ComPASS.mpi as mpi
from ComPASS.simulation_context import SimulationContext


##############################################################
### Simulation parameters & Physics
##############################################################
ComPASS.load_eos("immiscible2ph")
output_dir = "output-air-sparging"  # used for salome
ComPASS.set_output_directory_and_logfile(__file__)

p_atm = 1 * bar  # top pressure
# initial reservoir temperature - convert Celsius degrees to Kelvin degrees
Tall = degC2K(30)

qmass_inj = 5.0 / 3600.0  # 5m^3/h
pure_phase_molar_fraction = [[1.0, 0.0], [0.0, 1.0]]
liq_molar_fraction_inj = pure_phase_molar_fraction
Sg_top = 1.0
# # vDarcy = 4E-1 # m/s                         # qmass = velocity[m/s]*Section[m2]*Massevol
# # qmass_vdh = (40/3600)*997
# # qmass = (1/3600)*997 #                     # m3/h*density : qmass[kg/s]
# # qmass =  (0.81/3600)*997#                  # m3/h*density : qmass[kg/s] => 13.5 L/min => 0.81 m3/h
# # Qm = (rhof * 1E-3) / minute                # production/injection mass flowrate
# Surface_inj=2*np.pi*radius_well*aperture   # Surface of injection
# # Edgelength=2*np.pi*radius_well              # Surface of injection
# # qmass_inj = qmass/Surface_inj              # on surface (fracture opening) normalized flow rate
# Qinj = 13.5 # [L min-1]
# qmass_inj = 0.5 * ((((Qinj/1000)/60)*997)/Surface_inj)  # mass flux [kg s-1 m-2] if faces => (Qinj m-3 s-1 * 997 kg m-3) / (2*pi*0.08 m * aperture m)
# # qmass_inj = (((Qinj/1000)/60)*997)/Edgelength     # if edges => [kg s-1 m-1]
# qmass_ext = qmass_inj                      # extraction similar to injection

k_matrix = 5.3e-11  # domain permeability in m^2
phi_matrix = 0.39  # domain porosity
cell_thermal_cond = 0.0
gravity = ComPASS.get_gravity()
# rho_fluid_density = simulation.liquid_volumetric_mass_density((p_atm), Tall)


Topflag = 3
Bottomflag = 4
Externflag = 5
Injflag = 101
liquid_context = ComPASS.Context.liquid
diphasic_context = ComPASS.Context.diphasic
gas_context = ComPASS.Context.gas
Topz = None


if mpi.is_on_master_proc:

    ##############################################################
    ### Salome import / mesh creation
    ##############################################################

    nodes_file = "NODES.txt"
    tetras_file = "TETRAS.txt"

    vertices = np.loadtxt(nodes_file, usecols=(1, 2, 3))
    tetras = np.loadtxt(tetras_file, dtype=np.ulonglong, usecols=(1, 2, 3, 4))

    # consistency checks
    assert np.all(
        np.arange(1, vertices.shape[0] + 1)
        == np.loadtxt(nodes_file, dtype=np.ulonglong, usecols=0)
    )
    assert np.all(
        np.arange(1, tetras.shape[0] + 1)
        == np.loadtxt(tetras_file, dtype=np.ulonglong, usecols=0)
    )

    cells = MT.idarray(tetras - 1)
    mesh = MT.TetMesh.make(vertices, cells)  # Salome indexing starts at 1

    nodes2cell = {tuple(cell): ck for ck, cell in enumerate(cells)}

    # Flags limits
    TopNodes = MT.idarray(groups.Top) - 1  # Salome indexing starts at 1
    BottomNodes = MT.idarray(groups.Bottom) - 1
    ExternNodes = MT.idarray(groups.Extern) - 1
    InjectionNodes = MT.idarray(groups.Inj_zone_Nodes) - 1
    InjectionFaces = MT.idarray(groups.Inj_zone_Faces) - 1

    # Subdomains /volumes
    # MatrixDomain = MT.idarray(groups.Matrix)- 1

    # Find ids of cells (tetras)
    # MatrixDomain_compass = np.array([
    # nodes2cell[tuple(cell)] for cell in MatrixDomain
    # ]
    # )

    def set_global_flags():
        # nodes
        nodeflags = ComPASS.global_nodeflags()
        nodeflags[:] = 0
        nodeflags[TopNodes] = Topflag
        nodeflags[BottomNodes] = Bottomflag
        nodeflags[ExternNodes] = Externflag
        nodeflags[InjectionNodes] = Injflag
        # faces
        faceflags = ComPASS.global_faceflags()
        faceflags[InjectionFaces] = Injflag

    # for visualization only
    def flag_nodes(a):
        flag = np.zeros(mesh.nb_vertices)
        flag[a] = 666
        return flag

    # output a VTU mesh file with point propreties
    MT.to_vtu(
        mesh,
        output_dir + "-nodesflags",
        pointdata={
            "Top": flag_nodes(TopNodes),
            "Bottom": flag_nodes(BottomNodes),
            "Extern": flag_nodes(ExternNodes),
            "Inj": flag_nodes(InjectionNodes),
        },
    )

    # Limits
    vertices = mesh.vertices_array()
    x = vertices[:, 0]
    y = vertices[:, 1]
    z = vertices[:, -1]
    xmax, xmin = vertices[:, 0].max(), vertices[:, 0].min()
    ymax, ymin = vertices[:, 1].max(), vertices[:, 1].min()
    zmax, zmin = vertices[:, -1].max(), vertices[:, -1].min()
    Topz = zmax
    print("xmax=", xmax, "xmin=", xmin)
    print("ymax=", ymax, "ymin=", ymin)
    print("zmax=", zmax, "zmin=", zmin)

    def select_dirichlet_nodes():
        where = np.zeros(mesh.nb_vertices, dtype=np.bool)
        where[TopNodes] = True
        # where[ExternNodes] = True
        return where


if not mpi.is_on_master_proc:
    mesh = (
        select_dirichlet_nodes
    ) = k_matrix = phi_matrix = cell_thermal_cond = set_global_flags = None

Topz = mpi.communicator().bcast(Topz, root=0)

ComPASS.init(
    mesh=mesh,
    cell_permeability=k_matrix,
    cell_porosity=phi_matrix,
    cell_thermal_conductivity=cell_thermal_cond,
    set_dirichlet_nodes=select_dirichlet_nodes,
    set_global_flags=set_global_flags,
)


##############################################################
### Boundary and initial conditions
##############################################################
def lininterp(depths, top, gradient):
    assert np.all(depths) >= 0, "depths is not positif"
    return top + (gradient) * (depths)


def set_states(state, z, mask=True):
    lower = z <= 0
    upper = np.logical_not(lower)
    lower &= mask
    upper &= mask
    # liquid below 0
    state.context[lower] = liquid_context
    state.p[lower] = lininterp(
        -z[lower],
        p_atm,
        gravity
        * simulation.liquid_volumetric_mass_density(
            p_atm, Tall, pure_phase_molar_fraction
        ),
    )
    state.T[lower] = Tall
    state.S[lower] = [0, 1]
    state.C[lower] = pure_phase_molar_fraction
    # diphasic above 0,
    state.context[upper] = diphasic_context
    state.p[upper] = p_atm
    state.T[upper] = Tall
    state.S[upper, 0] = Sg_top * z[upper] / Topz  # gas
    state.S[upper, 1] = 1 - state.S[upper, 0]  # liquid
    state.C[upper] = pure_phase_molar_fraction


def set_Dirichlet_state(state, z):
    node_flags = ComPASS.nodeflags()
    # top nodes
    state.context[node_flags == Topflag] = gas_context
    state.p[node_flags == Topflag] = p_atm
    state.T[node_flags == Topflag] = Tall
    state.S[node_flags == Topflag] = [Sg_top, 1 - Sg_top]
    state.C[node_flags == Topflag] = pure_phase_molar_fraction
    # Extern nodes
    set_states(state, z, node_flags == Externflag)


def set_Neumann_fluxes():
    # find faces wich are InjectionFaces
    face_flags = ComPASS.faceflags()
    Neumann_inj_faces = np.zeros(len(face_flags), dtype=bool)
    Neumann_inj_faces[face_flags == Injflag] = True
    Density_inj = simulation.liquid_volumetric_mass_density(
        p_atm, Tall, liq_molar_fraction_inj
    )
    Pinj = lininterp(4.5, p_atm, gravity * Density_inj)
    # init Neumann ComPASS object
    Neumann = ComPASS.NeumannBC()
    Neumann.molar_flux[:] = liq_molar_fraction_inj[:] * qmass_inj
    Neumann.heat_flux = qmass_inj * (
        ComPASS.liquid_molar_enthalpy(Pinj, Tall, liq_molar_fraction_inj)
    )
    ComPASS.set_Neumann_faces(Neumann_inj_faces, Neumann)


def set_BC_and_initial_values():
    set_states(ComPASS.node_states(), ComPASS.vertices()[:, 2])
    set_states(ComPASS.cell_states(), ComPASS.compute_cell_centers()[:, 2])
    # set_Neumann_fluxes()
    set_Dirichlet_state(ComPASS.dirichlet_node_states(), ComPASS.vertices()[:, 2])


mpi.master_print("set initial and BC")
set_BC_and_initial_values()

# output a VTU mesh file with point propreties
MT.to_vtu(
    mesh,
    output_dir + "-initialization",
    pointdata={
        "P": np.ascontiguousarray(ComPASS.node_states().p),
        "T": np.ascontiguousarray(ComPASS.node_states().T),
        "Sg": np.ascontiguousarray(ComPASS.node_states().S[:, 0]),
        "Cag": np.ascontiguousarray(ComPASS.node_states().C[:, 0, 0]),
        "Cwl": np.ascontiguousarray(ComPASS.node_states().C[:, -1, -1]),
    },
    celldata={
        "P": np.ascontiguousarray(ComPASS.cell_states().p),
        "T": np.ascontiguousarray(ComPASS.cell_states().T),
        "Sg": np.ascontiguousarray(ComPASS.cell_states().S[:, 0]),
        "Cag": np.ascontiguousarray(ComPASS.cell_states().C[:, 0, 0]),
        "Cwl": np.ascontiguousarray(ComPASS.cell_states().C[:, -1, -1]),
    },
)


# ######################################
# ## Plot & Save results at each points OP1, OP2, OP3
# ######################################
# OP = [8.925, 0., 1.375]
# OPcoord = vertices[np.argmin(np.linalg.norm(vertices-OP, axis=1))]
# print('OP', OPcoord[0], OPcoord[1], OPcoord[2])

# def graph_frac(n, t):
#     #########################
#     # Data in fracture
#     #########################
#     vertices = np.rec.array(ComPASS.vertices())
#     states = ComPASS.node_states()
#     whereOP = (vertices[:,0] == OPcoord[0]) & (vertices[:,1] == OPcoord[1]) & (vertices[:,2] == OPcoord[2])
#     # whereOP2 = (vertices[:,0] == OP2[0]) & (vertices[:,1] == OP2[1]) & (vertices[:,2] == OP2[2])
#     # whereOP3 = (vertices[:,0] == OP3[0]) & (vertices[:,1] == OP3[1]) & (vertices[:,2] == OP3[2])

#     T_OP = K2degC(states.T[whereOP])
#     # T_OP2 = K2degC(states.T[whereOP2])
#     # T_OP3 = K2degC(states.T[whereOP3])
#     # p_OP1 = K2degC(states.p[whereOP1])
#     # p_OP2 = K2degC(states.p[whereOP2])
#     # p_OP3 = K2degC(states.p[whereOP3])

#     Tp_t = t, T_OP
#     Tp_t_list = np.hstack(np.hstack(Tp_t)).tolist()

#     save_Tp_OP123.append(Tp_t_list)
#     np.savetxt("save_Tp_OP123.txt", save_Tp_OP123)
#     # print('time==', t)
#     ##########################
#     ## Plot Results
#     ##########################
#     #plt.clf()
#     # cmap = plt.get_cmap("tab10")
#     # plt.figure(1)
#     # plt.plot(xcoord, T, color='red', label='Simulated '+ str(round(t,0)) +'s' )
#     # plt.legend()


##############################################################
### set linear solver properties
##############################################################
# save_Tp_OP123 = []

context = SimulationContext()
context.abort_on_ksp_failure = False
context.dump_system_on_ksp_failure = False
context.abort_on_newton_failure = False

timestep = TimeStepManager(
    initial_timestep=1,
    minimum_timestep=1e-7,
    maximum_timestep=10.0 * year,
    increase_factor=1.2,
    decrease_factor=0.2,
)

final_time = 100 * hour
output_period = final_time * 0.1

end_of_simu = standard_loop(
    final_time=final_time,
    output_period=output_period,
    context=context,
    time_step_manager=timestep,
    # iteration_callbacks = [graph_frac]
    # output_callbacks=[graph_frac]
)
print("time after the time loop", end_of_simu / year, "years")


# set_boundary_fluxes2new()

# final_time2 = Time1 + 23 * hour #2E4 #150 * year # * day
# initial_timestep2 = Time1 # Time1 +
# output_period2 = final_time2 * 1/23 #3
# Time2 = standard_loop(
#     final_time = final_time2,
#     output_period = output_period2,
#     time_step_manager = TimeStepManager(initial_timestep2, output_period2),
#     # output_callbacks=[graph_frac]
# )

# print('Time2', Time2)

# ### Printing info
# # p = ComPASS.cell_states().p
# # print('pressure:', p.min(), p.mean(), p.max())
# # T = ComPASS.cell_states().T
# # print('temperature:', T.min(), T.mean(), T.max())

# # states = ComPASS.node_states()
# # wherecenter = (vertices[:,0] == center[0]) & (vertices[:,1] == center[1]) & (vertices[:,2] == center[2])
# # Tcenter = states.T[wherecenter]
# # Pcenter = states.p[wherecenter]
# # Rhocenter = simulation.liquid_volumetric_mass_density(Pcenter, Tcenter)
# # Viscocenter = ComPASS.liquid_dynamic_viscosity(Pcenter, Tcenter)
# # print('Tcenter (K):', Tcenter, 'Pcenter (Pa):', Pcenter, 'Rhocenter:', Rhocenter, 'Viscocenter:', Viscocenter)


# # if ComPASS.mpi.communicator().size==1:
#     # assert ComPASS.mpi.is_on_master_proc
#     # try:
#         # import matplotlib
#         # matplotlib.use('Agg')
#         # import matplotlib.pyplot as plt
#     # else:
#         # x = ComPASS.cell_centers()[:, 0]
#         # states = ComPASS.cell_states()
#         # plt.clf()
#         # plt.plot(x, K2degC(states.T))
#         # plt.xlabel('x in meters')
#         # plt.ylabel('temperature')
#         # plt.savefig(ComPASS.to_output_directory('temperature'))
