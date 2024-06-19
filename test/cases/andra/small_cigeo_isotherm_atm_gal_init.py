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
import ComPASS.io.mesh as io

from ComPASS.timestep_management import TimeStepManager
from ComPASS.linalg.factory import linear_solver
from ComPASS.simulation import AlignmentMethod
from ComPASS.newton import Newton
import MeshTools as MT

ComPASS.set_output_directory_and_logfile(__file__)

""" ComPASS simulation"""
simulation = ComPASS.load_physics("diphasic")
gravity = 9.81
simulation.set_gravity(gravity)

# ---------------------------------------------
# Read mesh and identify the domains
# Load and create MeshTool mesh from binary file
# (binary file uses keywords "vertices", "hexas"...)
filename = "small_mesh_no_gal.npz"
mesh_data = np.load(filename)
nodes = mesh_data["vertices"]
mesh = MT.HexMesh.make(nodes, mesh_data["hexas"])
# extract the domain id from each cell, it will be
# a tag to identify the different rocktypes
cells_tag = mesh_data["domain"]
hexas = mesh_data["hexas"]
# Identify the nodes to set the atm BC
gallery_nodes_id = mesh_data["gallery_nodes_id"]

Pbottom = 6.0 * MPa
Tinit = degC2K(20)
Tgal = degC2K(20)
Tbento = degC2K(20)
rho = simulation.liquid_volumetric_mass_density
# identify cells with different properties
# bentonite_cells_selection = cells_tag == 5
EDZ_cells_selection = cells_tag == 6
porous_cells_selection = cells_tag == 7
topcover_cells_selection = cells_tag == 8
# bottomcover_cells_selection = cells_tag==9
# specific properties for specific cells_tag
# reservoir permeability in m^2
k_reservoir = 1e-20 * np.ones(len(cells_tag))
k_reservoir[EDZ_cells_selection] = 5e-18

# bulk thermal conductivity in W/m/K
K_reservoir = 2 * np.ones(len(cells_tag))

# reservoir porosity
omega = 0.15 * np.ones(len(cells_tag))
# omega[bentonite_cells_selection] = 0.35

y_min = min(nodes[:, 1])
y_max = max(nodes[:, 1])
z_min = min(nodes[:, -1])
z_max = max(nodes[:, -1])

EDZ_flag = 6
# bentonite_flag = 5
gallery_flag = 4

simulation.set_vanGenuchten_capillary_pressure()  # contains ox, cc, ...
from data.van_genuchten_kr import kr_functions

simulation.set_kr_functions(kr_functions)


def set_global_flags():
    def nodes_extraction(cellmask):
        selected_nodes_id = np.unique(hexas[cellmask])
        selected_nodes = np.zeros(len(nodes), dtype=bool)
        selected_nodes[selected_nodes_id] = True
        return selected_nodes

    nodeflags = simulation.global_nodeflags()
    EDZ_nodes = nodes_extraction(EDZ_cells_selection)
    # bentonite_nodes = nodes_extraction(bentonite_cells_selection)

    nodeflags[EDZ_nodes] = EDZ_flag
    # nodeflags[bentonite_nodes] = bentonite_flag
    nodeflags[gallery_nodes_id] = gallery_flag

    cellflags = simulation.global_cellflags()
    cellflags[EDZ_cells_selection] = EDZ_flag
    # cellflags[bentonite_cells_selection] = bentonite_flag


def set_global_rocktype():
    cellrocktype = simulation.global_cell_rocktypes()
    cellrocktype[:] = np.stack((cells_tag, cells_tag), axis=-1)


def select_dirichlet_nodes():
    on_the_bottom = lambda z: z == z_min
    z = simulation.global_vertices()[:, 2]
    return on_the_bottom(z)


def select_gallery_faces(node_flags):
    gallery_nodes_id = np.where(node_flags == gallery_flag)[0]
    nodes_by_faces = simulation.get_connectivity().NodebyFace
    gallery_faces = []
    for face_id, nodes in enumerate(nodes_by_faces):
        nodes_by_face = np.array(nodes, copy=False) - 1
        gallery_face = np.isin(nodes_by_face, gallery_nodes_id).all()
        if gallery_face:
            gallery_faces.append(face_id)
    return gallery_faces


def init_IC_and_BC(gallery_faces):
    dirichlet = simulation.dirichlet_node_states()
    all_states = simulation.all_states()
    node_states = simulation.node_states()

    X0 = simulation.build_state(simulation.Context.liquid, p=Pbottom, T=Tinit)
    dirichlet.set(X0)
    all_states.set(X0)
    # hydrostatic pressure everywhere
    liquid_phase = simulation.phase_index(simulation.Phase.liquid)
    all_states.p[:] = Pbottom - gravity * rho(X0.p, X0.T, X0.C[liquid_phase]) * (
        simulation.all_positions()[:, 2] - z_min
    )

    # init the freeflow state
    simulation.set_freeflow_faces(gallery_faces)
    is_ff = simulation.get_freeflow_nodes()  # mask
    assert (
        np.where(is_ff)[0] == gallery_nodes_id
    ).all(), "is_ff != gallery_nodes_id (can only be when Dirichlet nodes aside)"

    # With previous execution of the file, I know an estimated value for Sg
    pff = 1 * bar
    Xff = simulation.build_state(
        simulation.Context.diphasic_FF_no_liq_outflow,
        p=pff,
        T=Tinit,
        Sg=1.0e-3,
    )
    node_states.set(is_ff, Xff)


def export_initial_states():
    node_states = simulation.node_states()
    cell_states = simulation.cell_states()
    node_flags = simulation.nodeflags()
    cell_flags = simulation.cellflags()
    petrophysics = simulation.petrophysics()

    pointdata = {
        "dirichlet pressure": simulation.pressure_dirichlet_values(),
        "dirichlet temperature": K2degC(simulation.temperature_dirichlet_values()),
        "initial pressure": node_states.p,
        "initial temperature": K2degC(node_states.T),
        "initial gas saturation": node_states.S[:, 0],
        "node flags": node_flags,
    }
    celldata = {
        "initial pressure": cell_states.p,
        "initial gas saturation": cell_states.S[:, 0],
        "initial temperature": K2degC(cell_states.T),
        "phi": petrophysics.cell_porosity,
        "cells flag": cell_flags,
    }
    save_name = simulation.runtime.output_directory + "/initial_states_cigeo_simplified"
    io.write_mesh(simulation, save_name, pointdata=pointdata, celldata=celldata)


if __name__ == "__main__":
    simulation.init(
        mesh=mesh,
        set_dirichlet_nodes=select_dirichlet_nodes,
        cell_porosity=omega,
        cell_permeability=k_reservoir,
        cell_thermal_conductivity=K_reservoir,
        set_global_flags=set_global_flags,
        set_global_rocktype=set_global_rocktype,
    )

    node_flags = simulation.nodeflags()
    # identify the atmospheric faces
    gallery_faces = select_gallery_faces(node_flags)

    # init
    init_IC_and_BC(gallery_faces)

    # export_initial_states()

    simulation.alignment = AlignmentMethod.manual
    lsolver = linear_solver(simulation, tolerance=1e-8, direct=False)
    newton = Newton(simulation, 1e-6, 50, lsolver)
    tsmger = TimeStepManager(
        initial_timestep=5000,
        increase_factor=1.2,
        decrease_factor=0.5,
    )

    run_loop = lambda initial_time=0, final_time=1e6 * year, output_period=1e4 * year, output_every=None, nitermax=None: simulation.standard_loop(
        initial_time=initial_time,
        final_time=final_time,
        newton=newton,
        time_step_manager=tsmger,
        output_period=output_period,
        output_every=output_every,
        nitermax=nitermax,
    )

    current_time = run_loop(
        final_time=year,
        output_every=1,
        nitermax=4,
    )

    tsmger.increase_factor = 1.4
    newton.maximum_number_of_iterations = 8
    current_time = run_loop(
        initial_time=current_time,
        final_time=year,
        output_every=2,
    )

    # simulation.postprocess()
