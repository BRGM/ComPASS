# -*- coding: utf-8 -*-
#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

# load the diphasic physics, set some initialisations values,
# load the mesh, identify global flags and Dirichlet BC,
# set the rocktypes, Pc, kr functions
from small_cigeo_anisotropic_reload import *
import math

ComPASS.set_output_directory_and_logfile(__file__)

# identify an other domain for this simu
ZFC_nodes_id = mesh_data["ZFC_nodes_id"]
ZFC_flag = 9


def set_global_flags_with_ZFC_nodes():
    set_global_flags()
    nodeflags = simulation.global_nodeflags()
    nodeflags[ZFC_nodes_id] = ZFC_flag


simulation.init(
    mesh=mesh,
    set_dirichlet_nodes=select_dirichlet_nodes,
    cell_porosity=omega,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=K_reservoir,
    set_global_flags=set_global_flags_with_ZFC_nodes,
    set_global_rocktype=set_global_rocktype,
)

# reload state from initialization, some modifications are needed after
# Init done with atm nodes on a smaller mesh
# need a mapping
snapshot_directory = "output-small_cigeo_anisotropic_atm_gal_init"
reload_states(snapshot_directory)

# the known states have been reloaded
# necessary to change the context of the nodes which were FF
# the other values of the states are kept
reset_freeflow_context()

# necessary to init the states of the nodes and cells
# inside the circulation gallery and bentonite
# (were eliminated in the init state and replaced by atmospheric BC)
init_gallery_states()

# export_initial_states()


# identify ZFC_faces for Neumann
def face_identify(flag):
    node_flags = simulation.nodeflags()
    domain_nodes_id = np.where(node_flags == flag)[0]
    nodes_by_faces = simulation.get_connectivity().NodebyFace
    domain_faces = []
    for face_id, nodes in enumerate(nodes_by_faces):
        nodes_by_face = np.array(nodes, copy=False) - 1
        domain_face = np.isin(nodes_by_face, domain_nodes_id).all()
        if domain_face:
            1 / 0
            domain_faces.append(face_id)
    return domain_faces, domain_nodes_id


ZFC_faces, ZFC_nodes_id = face_identify(ZFC_flag)
# ZFC_faces is empty because of our mesh
# this function was tested and approved with another domain
breakpoint()


# def convert_colis_m2(length_colis,width_colis,heat_flux):
#     heat_flux = heat_flux / (length_colis * width_colis)  # in w/m^2
#     return heat_flux

# called at the end of each time iteration
# no Neumann flux at the first iteration
def Neumann_change_molarflux(tick):
    length_colis = 2.1  # in m
    width_colis = 1  # in m

    entiere_time = math.ceil(tick.time / year) - 1
    time = np.concatenate(
        (
            np.arange(0, 10, 1),
            np.arange(10, 100, 5),
            np.arange(100, 300, 10),
            np.arange(300, 1000, 100),
            np.arange(1e3, 1e4, 1e3),
            np.arange(1e4, 1e5, 1e4),
            np.arange(1e5, 1.1e6, 1e5),
        )
    )
    heat_flux = np.arange(2100, 0.119, -25.3)
    # inspired from andra values
    # the values should bre readen from excel file
    Neumann_heat_flux = heat_flux / (length_colis * width_colis)
    index_of_heat_flux = np.flatnonzero(entiere_time == time)
    Neumann_heat_flux = heat_flux[index_of_heat_flux]
    molar_flux = [0, 0]
    Neumann = ComPASS.NeumannBC(molar_flux, Neumann_heat_flux)
    simulation.set_Neumann_faces(ZFC_faces, Neumann)


simulation.alignment = AlignmentMethod.manual
# simulation.alignment = AlignmentMethod.inverse_diagonal # it is worse
lsolver = linear_solver(simulation, tolerance=1e-8, direct=False)
newton = Newton(simulation, 1e-6, 10, lsolver)
tsmger = TimeStepManager(
    initial_timestep=0.01 * year,
    increase_factor=1.2,
    decrease_factor=0.5,
)

run_loop = lambda initial_time=0, final_time=1000 * year, output_period=1e4 * year, output_every=None, nitermax=None: simulation.standard_loop(
    initial_time=initial_time,
    final_time=final_time,
    newton=newton,
    time_step_manager=tsmger,
    output_period=output_period,
    output_every=output_every,
    nitermax=nitermax,
    iteration_callbacks=[Neumann_change_molarflux],
)


current_time = run_loop(
    final_time=30 * year,
    output_period=5 * year,
    nitermax=10,
)

tsmger.increase_factor = 1.3
newton.maximum_number_of_iterations = 6
current_time = run_loop(
    initial_time=current_time,
    final_time=30 * year,
    output_period=5 * year,
)

breakpoint()
tsmger.maximum = 5.1 * year
current_time = run_loop(
    initial_time=current_time,
    final_time=500 * year,
    output_period=50 * year,
)
breakpoint()
tsmger.maximum = 700 * year
current_time = run_loop(
    initial_time=current_time,
    output_period=100 * year,
)
