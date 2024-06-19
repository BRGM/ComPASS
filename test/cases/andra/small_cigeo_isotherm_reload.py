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
filename = "small_mesh_with_gal.npz"
mesh_data = np.load(filename)
nodes = mesh_data["vertices"]
mesh = MT.HexMesh.make(nodes, mesh_data["hexas"])
# extract the domain id from each cell, it will be
# a tag to identify the different rocktypes
cells_tag = mesh_data["domain"]
hexas = mesh_data["hexas"]

# identify cells with different properties
circulation_cells_selection = cells_tag == 4
bentonite_cells_selection = cells_tag == 5
EDZ_cells_selection = cells_tag == 6
porous_cells_selection = cells_tag == 7
topcover_cells_selection = cells_tag == 8
# bottomcover_cells_selection = cells_tag==9
# specific properties for specific cells_tag
# reservoir permeability in m^2
k_reservoir = 1e-20 * np.ones(len(cells_tag))
k_reservoir[circulation_cells_selection] = 5e-17
k_reservoir[EDZ_cells_selection] = 5e-18

# bulk thermal conductivity in W/m/K
K_reservoir = 2 * np.ones(len(cells_tag))

# reservoir porosity
omega = 0.15 * np.ones(len(cells_tag))
omega[circulation_cells_selection] = 0.4
omega[bentonite_cells_selection] = 0.35

y_min = min(nodes[:, 1])
y_max = max(nodes[:, 1])
z_min = min(nodes[:, -1])
z_max = max(nodes[:, -1])

EDZ_flag = 6
bentonite_flag = 5
circulation_flag = 4

simulation.set_vanGenuchten_capillary_pressure()  # contains ox, cc, ...
from data.van_genuchten_kr import kr_functions

simulation.set_kr_functions(kr_functions)


def set_global_flags():
    def nodes_extraction(cellmask):
        selected_nodes_id = np.unique(hexas[cellmask])
        selected_nodes = np.zeros(len(nodes), dtype=bool)
        selected_nodes[selected_nodes_id] = True
        return selected_nodes

    cellflags = simulation.global_cellflags()
    cellflags[EDZ_cells_selection] = EDZ_flag
    cellflags[bentonite_cells_selection] = bentonite_flag
    cellflags[circulation_cells_selection] = circulation_flag

    nodeflags = simulation.global_nodeflags()
    circulation_nodes = nodes_extraction(circulation_cells_selection)
    nodeflags[circulation_nodes] = circulation_flag
    bentonite_nodes = nodes_extraction(bentonite_cells_selection)
    nodeflags[bentonite_nodes] = bentonite_flag


def set_global_rocktype():
    cellrocktype = simulation.global_cell_rocktypes()
    cellrocktype[:] = np.stack((cells_tag, cells_tag), axis=-1)


def select_dirichlet_nodes():
    on_the_bottom = lambda z: z == z_min
    z = simulation.global_vertices()[:, 2]
    return on_the_bottom(z)


def reload_states(snapshot_directory):
    mapping = {}
    mapping["node"] = mesh_data["node_mapping"]
    mapping["cell"] = mesh_data["cell_mapping"]
    current_time = simulation.reload_snapshot(
        snapshot_directory,
        mapping=mapping,
        reset_dirichlet=True,  # done by default
    )
    return current_time


def reset_freeflow_context():
    node_states = simulation.node_states()
    FF_diphasic = simulation.Context.diphasic_FF_no_liq_outflow.value
    FF_gas = simulation.Context.gas_FF_no_liq_outflow.value
    FF_liquid = simulation.Context.diphasic_FF_liq_outflow.value
    node_states.context[
        node_states.context[:] == FF_diphasic
    ] = simulation.Context.diphasic
    node_states.context[node_states.context[:] == FF_liquid] = simulation.Context.liquid
    node_states.context[node_states.context[:] == FF_gas] = simulation.Context.gas


def init_gallery_states():
    Tgal = degC2K(20)
    Xgal0 = simulation.build_state(
        simulation.Context.diphasic,
        p=1 * bar,
        T=Tgal,
        Sg=0.3,
    )

    node_states = simulation.node_states()
    cell_states = simulation.cell_states()
    node_flags = simulation.nodeflags()
    cell_flags = simulation.cellflags()

    node_states.set(node_flags == bentonite_flag, Xgal0)
    cell_states.set(cell_flags == bentonite_flag, Xgal0)
    node_states.set(node_flags == circulation_flag, Xgal0)
    cell_states.set(cell_flags == circulation_flag, Xgal0)


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
    save_name = simulation.runtime.output_directory + "/initial_states_small_cigeo"
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

    # reload state from initialization, some modifications are needed after
    # Init done with atm nodes on a smaller mesh
    # need a mapping
    snapshot_directory = "output-small_cigeo_isotherm_atm_gal_init"
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

    simulation.alignment = AlignmentMethod.manual
    lsolver = linear_solver(simulation, tolerance=1e-8, direct=False)
    newton = Newton(simulation, 1e-6, 50, lsolver)
    tsmger = TimeStepManager(
        initial_timestep=25 * day,
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
        final_time=30 * year,
        nitermax=1,
    )

    newton.maximum_number_of_iterations = 8
    current_time = run_loop(
        initial_time=current_time,
        final_time=30 * year,
        output_every=2,
        nitermax=10,
    )

    tsmger.increase_factor = 1.3
    newton.maximum_number_of_iterations = 6
    current_time = run_loop(
        initial_time=current_time,
        final_time=30 * year,
        output_period=5 * year,
    )

    # breakpoint()
    # tsmger.maximum = 5.1 * year
    current_time = run_loop(
        initial_time=current_time,
        final_time=3000 * year,  # 5000 * year
        output_period=1e2 * year,
    )

    # breakpoint()
    # tsmger.maximum = 1e4 * year
    # current_time = run_loop(
    #     initial_time = current_time,
    #     output_period = 1e5 * year,
    # )

    # simulation.postprocess()
