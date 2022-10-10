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
from ComPASS.newton import Newton
from ComPASS.timeloops import standard_loop
import MeshTools as MT

import gmsh_reader

ComPASS.set_output_directory_and_logfile(__file__)

""" ComPASS simulation"""
simulation = ComPASS.load_physics("diphasic")
gravity = 9.81
simulation.set_gravity(gravity)
bulk_thermal_conductivity = 2.0  # W . K^-1 . m^-1
rock_internal_energy = 2.0e6  # volumetric heat capacity J/K/m^3
simulation.set_rock_volumetric_heat_capacity(rock_internal_energy)


""" mesh : box of size 0 < X < 20 ; -20 < Y < 20 ; -10 < Z < 10 """
""" two rocktypes : Callova Oxfordian clay (ox) and concrete layer (cc) """
filename = "./small_disposal_gallery.msh"
nodes, elements = gmsh_reader.retrieve_mesh_elements(filename)

# grep 3d elements, they are the cells of the mesh
cells = []
cells_tag = []
for elt, tag in elements:
    if type(elt) in (MT.Tetrahedron, MT.Wedge, MT.Hexahedron):
        cells.append(elt)
        cells_tag.append(tag[0])
cells_tag = np.array(cells_tag)

mesh = MT.HybridMesh.create(nodes, cells)

""" Dirichlet nodes at the bottom, top and gallery """
r = 2.625
y_min = -20
y_max = 20
z_min = -10
z_max = 10

eps = 1e-6

bottom_node = [abs(z - z_min) <= eps for x, y, z in nodes]
bottom_flag = 3
top_node = [abs(z - z_max) <= eps for x, y, z in nodes]
top_flag = 2
gallery_node = [y * y + z * z <= r * r + eps for x, y, z in nodes]
gallery_flag = 1


def set_physical_flags():
    nodeflags = simulation.global_nodeflags()
    nodeflags[:] = 0

    nodeflags[bottom_node] = bottom_flag
    nodeflags[top_node] = top_flag
    nodeflags[gallery_node] = gallery_flag

    cellflags = simulation.global_cellflags()
    cellflags[:] = cells_tag


def select_dirichlet_nodes():
    nodeflags = simulation.global_nodeflags()
    return nodeflags != 0


""" The domain is composed of two rocktypes of index 1 (ox) or 2 (cc) """


def select_global_rocktype():
    cellflags = simulation.global_cellflags()
    rocktype = 1 + ((cellflags - 1) % 2)

    cellrocktype = simulation.global_cell_rocktypes().reshape((-1, 2))
    cellrocktype[:] = np.stack((rocktype, rocktype), axis=-1)


def set_permeabilities():  # domain permeability in m^2
    cellrocktypes = simulation.global_cell_rocktypes()
    k_matrices = np.zeros(np.size(cells_tag))
    k_matrices[cellrocktypes[:, 0] == 1] = 5e-20  # ox
    k_matrices[cellrocktypes[:, 0] == 2] = 1e-18  # cc
    return k_matrices


def set_porosity():  # domain porosity
    cellrocktypes = simulation.global_cell_rocktypes()
    phi_matrices = np.zeros(np.size(cells_tag))
    phi_matrices[cellrocktypes[:, 0] == 1] = 0.15  # ox
    phi_matrices[cellrocktypes[:, 0] == 2] = 0.30  # cc
    return phi_matrices


simulation.init(
    mesh=mesh,
    cell_thermal_conductivity=bulk_thermal_conductivity,
    set_global_flags=set_physical_flags,
    set_global_rocktype=select_global_rocktype,
    set_dirichlet_nodes=select_dirichlet_nodes,
    cell_porosity=set_porosity,
    cell_permeability=set_permeabilities,
)

""" Set petrophyics functions after initialization """
simulation.set_vanGenuchten_capillary_pressure()  # contains ox and cc
from data.van_genuchten_kr import kr_functions

simulation.set_kr_functions(kr_functions)


""" Initial and Dirichlet values """


def hydrostatic_pressure(Ptop, T, C, dz):
    rho = simulation.liquid_volumetric_mass_density(Ptop, T, C)
    return Ptop + gravity * rho * dz


def set_diphasic_equilibrium(Pg, T, Hur):
    Psat = simulation.Psat(T)
    Cwg = Hur * Psat / Pg
    Cag = 1.0 - Cwg

    T1 = 293.0
    T2 = 353.0
    H1 = 6.0e9
    H2 = 10.0e9
    Ha = H1 + (H2 - H1) * (T - T1) / (T2 - T1)

    Cal = Cag * Pg / Ha
    Cwl = 1.0 - Cal

    RZetal = 8.314 * 1000.0 / 0.018
    # fw^l = Psat * exp(- (Psat-pl)/(zeta^l R T))
    # fw^l * Cwl = fw^g * Cwg = Pg * Cwg
    Pl = Psat + np.log(Hur / Cwl) * RZetal * T
    Pc = Pg - Pl
    # rocktype = 2 along the gallery
    n = 1.54
    m = 1.0 - 1.0 / n
    Pr = 2.0e6  # Pa
    Slr = 0.01
    Sgr = 0.0
    Sl = Slr + (1.0 - Slr - Sgr) * (1.0 + (Pc / Pr) ** n) ** (-m)

    return 1.0 - Sl


# in the porous media
Cal = 0.0
Cl = [Cal, 1.0 - Cal]
PlTop = 4e6
TPor = 303  # K
Pc = 0.0  # Sg = 0, there is no entry pressure
# in the galery "gallery_flag"
PgGal = 1.0e5
TGal = 303.0
HurGal = 0.5

# Dirichlet nodes need to be init at equilibrum, thus use build_state
def set_Dirichlet_state(state):
    node_flags = simulation.nodeflags()

    # in the galery value, get Sg from Hur for build_state
    SgGal = set_diphasic_equilibrium(PgGal, TGal, HurGal)
    # in build_state p = Pref = Pg in diphasic physics
    XGal = simulation.build_state(
        simulation.Context.diphasic,
        p=PgGal,
        T=TGal,
        Sg=SgGal,
        rocktype=2,
    )
    state.set(node_flags == gallery_flag, XGal)
    # in porous media : top
    # in build_state p = Pref = Pg in diphasic physics
    # compute Pc and then Pg = Pc + Pl
    PgTop = Pc + PlTop
    XPorTop = simulation.build_state(simulation.Context.liquid, p=PgTop, T=TPor)
    state.set(node_flags == top_flag, XPorTop)
    # bottom
    PlBot = hydrostatic_pressure(PlTop, TPor, Cl, z_max - z_min)
    PgBot = Pc + PlBot
    XPorBot = simulation.build_state(simulation.Context.liquid, p=PgBot, T=TPor)
    state.set(node_flags == bottom_flag, XPorBot)


# set_states not called with Dirichlet nodes, beacause they need to be at equilibrium,
# only T, Pref (ie Pg), Sl and Cal are necessary with liquid context
def set_states(states, dz):
    states.context[:] = simulation.Context.liquid
    states.p[:] = hydrostatic_pressure(PlTop, TPor, Cl, dz)
    states.T[:] = TPor
    states.S[:] = [0, 1]
    states.C[:] = [[0.8, 0.2], Cl]  # [[Cag, Cwg], [Cal, Cwl]]


def set_initial_bc_values():
    set_states(simulation.node_states(), z_max - simulation.vertices()[:, 2])
    set_states(simulation.cell_states(), z_max - simulation.cell_centers()[:, 2])
    set_Dirichlet_state(simulation.dirichlet_node_states())


set_initial_bc_values()


def export_initial_states():
    node_states = simulation.node_states()
    cell_states = simulation.cell_states()
    petrophysics = simulation.petrophysics()

    pointdata = {
        "dirichlet pressure": simulation.pressure_dirichlet_values(),
        "dirichlet temperature": K2degC(simulation.temperature_dirichlet_values()),
        "initial pressure": node_states.p,
        "initial temperature": K2degC(node_states.T),
        "initial gas saturation": node_states.S[:, 0],
    }
    celldata = {
        "initial pressure": cell_states.p,
        "initial gas saturation": cell_states.S[:, 0],
        "initial temperature": K2degC(cell_states.T),
        "phi": petrophysics.cell_porosity,
    }
    save_name = simulation.runtime.output_directory + "/initial_states_andra"
    io.write_mesh(simulation, save_name, pointdata=pointdata, celldata=celldata)


# export_initial_states()

final_time = 1e4 * year
output_period = 0.05 * final_time
intermediate_time = 5 * year

timestep = TimeStepManager(
    initial_timestep=100.0,
    minimum_timestep=1e-1,
    maximum_timestep=output_period,
    increase_factor=1.2,
    decrease_factor=0.2,
)

# linear solver tolerance must be smaller than the newton tolerance
lsolver = linear_solver(simulation, tolerance=1e-8, legacy=False, from_options=True)
newton = Newton(simulation, 1e-7, 14, lsolver)

# call twice the standard_loop to have more outputs at the beggining
# use output_every
current_time = standard_loop(
    simulation,
    final_time=intermediate_time,
    newton=newton,
    output_every=10,
    nitermax=209,
    time_step_manager=timestep,
)

timestep.current_step = current_time / 30.0
# use output_period
standard_loop(
    simulation,
    initial_time=current_time,
    final_time=final_time,
    newton=newton,
    output_period=output_period,
    time_step_manager=timestep,
)

# simulation.post_process()
