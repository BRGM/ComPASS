# -*- coding: utf-8 -*-
# to create the submission mesh
# python make_structured_mesh.py --variant B -nx 420 -ny 140
# mpirun -n 128 SPE11_b.py

import numpy as np
from math import sqrt
import ComPASS
import ComPASS.io.mesh as io
from ComPASS.utils.units import *
from ComPASS.timeloops import TimeStepManager
from ComPASS.utils.various import tensor_coordinates
from ComPASS.newton import Newton
from ComPASS.linalg.factory import linear_solver
from ComPASS.dumps_spe11 import DumperSPE11
from pathlib import Path
import MeshTools as MT
import ComPASS.mpi as mpi
from geometries.msh2compass import convert_mesh
import shutil


ComPASS.set_output_directory_and_logfile(__file__)

""" Parameters"""
year: float = 31536.0e3  # seconds
p0: float = 211.74 * bar  # 221.74434 * bar # Pressure approximation at the top
T0: float = degC2K(70.0)  # Temperature at the bottom (zmin = 0)
pwell: float = 300 * bar  # Pressure at well1
Twell: float = degC2K(10.0)  # Temperature at the well1
pure_phase_molar_fraction = [[1.0, 0.0], [0.0, 1.0]]  # [[CO2_g, W_g], [CO2_l, W_l]]
bottom_boundary: float = 0.0
top_boundary: float = 1200.0
xmin_boundary: float = 0.0
xmax_boundary: float = 8400.0

# Wellnodes coordinates: well1(2700, 1000, 300) and well2(2700, 4000, 300)
# x1 = 2700., y1 = 50. , z1 = 300 # y1=100/2 for Spe11_b, with well_depth = 100

# The SPE11b case is in 2D, however ComPASS is only in 3D
# (the mesh contains one cell in the y direction).
# Then the well is a horizontal line instead of a point.
# The molar flux is defined over the perimeter of the well,
# thus we mutilpy by 2*pi*r*well_depth to get the molar_flux.
# Also, the given flux (0.035) is a mass flux, we transform it into a molar flux using M_CO2=44.01 g/mol.
# *** WARNING ***:
# IF CHANGING THE MESH, change well_depth,
# m (with nx = 300, ny = 100, change if you change the mesh !!!)
# well_depth: float = 15
# m (with nx = 450, ny = 150, change if you change the mesh !!!)
# well_depth: float = 2 / 3
# m (with nx = 360, ny = 120, change if you change the mesh !!!)
# well_depth: float = 12.5
# m (with nx = 420, ny = 140, change if you change the mesh !!!)
well_depth: float = 10.7
well_radius: float = 0.15
# wells coordinates : well1, well2
wells_coord = np.array([[2700, well_depth / 2, 300], [5100, well_depth / 2, 700]])
molar_flux = 0.035 / 44.01e-3 * 2 * 3.14 * well_radius * well_depth

k_Facies1: float = 1e-16  # domain permeability in m^2
k_Facies2: float = 1e-13  # domain permeability in m^2
k_Facies3: float = 2e-13  # domain permeability in m^2
k_Facies4: float = 5e-13  # domain permeability in m^2
k_Facies5: float = 1e-12  # domain permeability in m^2
k_Facies6: float = 2e-12  # domain permeability in m^2
k_Facies7: float = 1e-20  # domain permeability in m^2 should be 0
aniso = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 0.1]])
k_Facies = np.array(
    [
        k_Facies1 * aniso,
        k_Facies2 * aniso,
        k_Facies3 * aniso,
        k_Facies4 * aniso,
        k_Facies5 * aniso,
        k_Facies6 * aniso,
        k_Facies7 * aniso,
    ]
)

phi_Facies1: float = 0.1  # porosity
phi_Facies2: float = 0.2  # porosity
phi_Facies3: float = 0.2  # porosity
phi_Facies4: float = 0.2  # porosity
phi_Facies5: float = 0.25  # porosity
phi_Facies6: float = 0.35  # porosity
phi_Facies7: float = 0.0001  # 0.0001 porosity should be 0
phi_Facies = np.array(
    [
        phi_Facies1,
        phi_Facies2,
        phi_Facies3,
        phi_Facies4,
        phi_Facies5,
        phi_Facies6,
        phi_Facies7,
    ]
)

kth_Facies1: float = 1.9  # rock heat conductivity in W/m/K
kth_Facies2: float = 1.25  # rock heat conductivity in W/m/K
kth_Facies3: float = 1.25  # rock heat conductivity in W/m/K
kth_Facies4: float = 1.25  # rock heat conductivity in W/m/K
kth_Facies5: float = 0.92  # rock heat conductivity in W/m/K
kth_Facies6: float = 0.26  # rock heat conductivity in W/m/K
kth_Facies7: float = 2.0  # rock heat conductivity in W/m/K
kth_Facies = np.array(
    [
        kth_Facies1,
        kth_Facies2,
        kth_Facies3,
        kth_Facies4,
        kth_Facies5,
        kth_Facies6,
        kth_Facies7,
    ]
)


heat_capacity_rock: float = 850.0  # [J.kg-1.K-1] spe11b 0.85 kJ.kg-1.K-1
rho_rock_density: float = 2500.0  # [kg.m-3]

gravity: float = 9.80665

# Tgeo(x) = 70 - 0.025z
Tgrad = -0.025  # 25°C/km # Geothermal gradient

""" MESH """
grid, rocktype = convert_mesh(
    "./geometries/spe11b_structured.msh",
    well_depth,
    "cells_rkt.vtu",
)
# identify the nodes in Dirichlet rocktype cells
if mpi.is_on_master_proc:
    dir_rkt = (rocktype >= 2) & (rocktype <= 5)
    dir_nodes_mask = np.zeros(len(grid.vertices), dtype=bool)
    dir_nodes_mask[np.unique(grid.cell_nodes[dir_rkt])] = True

Rightflag = 1
Leftflag = 2
Topflag = 3
Bottomflag = 4
Well1flag = 5
Well2flag = 6
BottomRightflag = 7
DirichletVerticalflag = 8

""" ComPASS simulation """
simulation = ComPASS.load_physics("diphasicCO2")

""" kr and Capillary functions """
# set pc for rocktypes from 1 to 6
simulation.set_extendedBrooksCorey_pc_SPE11b()
# Brooks Corey kr
from data.brooks_corey_kr import kr_functions_SPE11b as kr_functions

simulation.set_kr_functions(kr_functions)


""" global parameters """
liquid_context = simulation.Context.liquid
gas_context = simulation.Context.gas
diphasic_context = simulation.Context.diphasic

simulation.set_gravity(gravity)
simulation.set_rock_volumetric_heat_capacity(rho_rock_density * heat_capacity_rock)


def set_global_flags():
    # nodes
    nodeflags = simulation.global_nodeflags()
    x = simulation.global_vertices()[:, 0]
    nodeflags[dir_nodes_mask & (x <= xmin_boundary)] = DirichletVerticalflag
    nodeflags[dir_nodes_mask & (x >= xmax_boundary)] = DirichletVerticalflag

    cellflags = simulation.global_cellflags()
    # get the global index of the cell containing each well
    centers = simulation.compute_global_cell_centers()
    Well1_cell = np.argmin(np.linalg.norm(centers - wells_coord[0], axis=1))
    Well2_cell = np.argmin(np.linalg.norm(centers - wells_coord[1], axis=1))
    cellflags[Well1_cell] = Well1flag
    cellflags[Well2_cell] = Well2flag


def set_global_rocktype():
    cellrocktype = simulation.global_cell_rocktypes()
    # duplicate because rocktype of Darcy and Fourier
    cellrocktype[:, 0] = rocktype
    cellrocktype[:, 1] = rocktype
    cellrocktype[rocktype == 7, 0] = 0  # no pc
    cellrocktype[rocktype == 7, 1] = 0  # no pc


""" Functions to define permeabilities: Units (Kx=Ky et Kz=0.1*Kx) """


def subdomains_permeability():
    k_cells = k_Facies[rocktype - 1]
    return k_cells


def subdomains_porosity():
    phi_cells = phi_Facies[rocktype - 1]
    return phi_cells


def subdomains_conductivity():
    kth_cells = kth_Facies[rocktype - 1]
    return kth_cells


""" Functions to define Dirichlet boundary condition """


def P_dirichlet_nodes():
    x = simulation.global_vertices()[:, 0]
    z = simulation.global_vertices()[:, 2]
    # BottomRight_nodes
    return np.logical_and(z <= bottom_boundary, x >= xmax_boundary)


def Temp_dirichlet_nodes():
    z = simulation.global_vertices()[:, 2]
    return np.logical_or(z <= bottom_boundary, z >= top_boundary)


dumper = DumperSPE11(simulation)

simulation.init(
    mesh=grid,
    cell_permeability=subdomains_permeability,
    cell_porosity=subdomains_porosity,
    cell_thermal_conductivity=subdomains_conductivity,
    set_temperature_dirichlet_nodes=Temp_dirichlet_nodes,
    set_pressure_dirichlet_nodes=P_dirichlet_nodes,
    set_global_flags=set_global_flags,
    # init the rocktypes (for the pc)
    set_global_rocktype=set_global_rocktype,
)

if mpi.is_on_master_proc:
    # move cells_rkt.vtu into the output_dir, otherwise not copied
    # back from computation node in Leto
    shutil.move("cells_rkt.vtu", ComPASS.to_output_directory("cells_rkt.vtu"))


""" Initial conditions """


def lininterp(z, top, gradient):
    return top + (gradient) * (z)


def set_initial_states():
    X0 = simulation.build_state(liquid_context, p=300 * bar, T=T0)

    def set_states(states, z):
        # init with liquid context at 300 bar and T0
        states.set(X0)
        # linear gradient for the pressure and temperature
        states.p[:] = lininterp(
            300 - z,
            300 * bar,
            gravity * 1005.0,
        )
        states.T[:] = lininterp(
            z,
            T0,
            Tgrad,
        )

    z = simulation.all_positions()[:, 2]
    set_states(simulation.all_states(), z)
    # # the Dirichlet node states could be copied from the actual states values
    # simulation.reset_dirichlet_nodes(dirichlet_nodes)
    z = simulation.vertices()[:, 2]
    set_states(simulation.dirichlet_node_states(), z)


set_initial_states()


def export_states(filename):
    """Function to export values as vtufile"""
    cellrocktype = simulation.cell_rocktypes()
    node_states = simulation.node_states()
    cell_states = simulation.cell_states()
    rho_l = simulation.cpp_liquid_molar_density
    rho_co2 = simulation.cpp_gas_molar_density
    petrophysics = simulation.petrophysics()

    pointdata = {
        "dirichlet pressure": simulation.pressure_dirichlet_values(),
        "dirichlet temperature": K2degC(simulation.temperature_dirichlet_values()),
        "pressure": node_states.p,
        "temperature": K2degC(node_states.T),
        # "context": node_states.C,
        # "saturation": node_states.S[:,0],
        "liquid molar density": rho_l(node_states.p, node_states.T, node_states.S),
        "liquid mass density": rho_l(node_states.p, node_states.T, node_states.S)
        * (44.01e-3 * node_states.C[:, 1, 0] + 18.0e-3 * node_states.C[:, 1, 1]),
        "co2_molar_density": rho_co2(node_states.p, node_states.T, node_states.S),
        "co2_mass_density": rho_co2(node_states.p, node_states.T, node_states.S)
        * 44.01e-3,
        # "Psat for T reservoir": simulation.Psat(node_states.T),
        # "Tsat for p reservoir": K2degC(simulation.Tsat(node_states.p)),
    }
    celldata = {
        "rocktype": cellrocktype,
        "pressure": cell_states.p,
        # "saturation": cell_states.S,
        "temperature": K2degC(cell_states.T),
        "phi": petrophysics.cell_porosity,
        "co2_mass_density": rho_co2(cell_states.p, cell_states.T, cell_states.S)
        * 44.01e-3,
        "liquid mass density": rho_l(cell_states.p, cell_states.T, cell_states.S)
        * (44.01e-3 * cell_states.C[:, 1, 0] + 18.0e-3 * cell_states.C[:, 1, 1]),
    }
    celldata.update(
        tensor_coordinates(petrophysics.cell_permeability, "k", diagonal_only=True)
    )
    celldata.update(
        tensor_coordinates(
            petrophysics.cell_thermal_conductivity, "k_th", diagonal_only=True
        )
    )
    io.write_mesh(simulation, filename, pointdata=pointdata, celldata=celldata)


export_states(ComPASS.to_output_directory("states_0yr"))


""" Solver parameters """
tsmger = TimeStepManager(
    initial_timestep=30 * year,
    minimum_timestep=1,  # s
    increase_factor=10,
    decrease_factor=0.5,
)

""" Simulation """
final_time0 = 1e3 * year
simu_time = simulation.standard_loop(
    final_time=final_time0,
    time_step_manager=tsmger,
    nb_output=2,
    output_after_loop=True,
    output_before_start=True,
    dumper=dumper,
)

# simulation.postprocess(time_unit="second")


""" Set Dirichlet nodes at xmin and xmax boundaries with initial states """
nodeflags = simulation.nodeflags()
dirichlet_nodes = np.zeros(len(nodeflags), dtype=bool)
dirichlet_nodes[nodeflags == DirichletVerticalflag] = True

simulation.reset_dirichlet_nodes(dirichlet_nodes)


""" Export values as vtufile"""
export_states(ComPASS.to_output_directory("states_1e3yr"))


# set CO2 and heat sources for well 1
def set_molar_heat_sources(well_flag):
    cellflags = simulation.cellflags()
    simulation.cell_molar_sources_vol()[cellflags == well_flag, 0] = molar_flux
    gas_mol_enth = simulation.cpp_gas_molar_enthalpy(pwell, Twell, [1, 0])
    simulation.cellthermalsource()[cellflags == well_flag] = molar_flux * gas_mol_enth


set_molar_heat_sources(Well1flag)


""" Update Solver parameters """
lsolver = linear_solver(simulation, direct=False)
newton = Newton(simulation, 1e-5, 50, lsolver, maxit_before_decrease=4)
tsmger.current_step = 0.01 * year
# tsmger.maximum_timestep = 0.05 * year
tsmger.increase_factor = 1.1
tsmger.decrease_factor = 0.7

final_time1 = 1025 * year
simu_time = simulation.standard_loop(
    initial_time=simu_time,
    final_time=final_time1,
    time_step_manager=tsmger,
    output_period=0.5 * year,
    newton=newton,
    dumper=dumper,
)
# simulation.postprocess(time_unit="second")
# export also the densities...
export_states(ComPASS.to_output_directory("states_1025yr"))

# reactivate the CO2 and heat sources for well 2
set_molar_heat_sources(Well2flag)


tsmger.current_step = 0.01 * year

final_time2 = 1050 * year
simu_time = simulation.standard_loop(
    initial_time=simu_time,
    final_time=final_time2,
    time_step_manager=tsmger,
    output_period=1 * year,
    newton=newton,
    dumper=dumper,
)
# simulation.postprocess(time_unit="second")
# export also the densities...
export_states(ComPASS.to_output_directory("states_1050yr"))

# reset the molar and heat sources to end without sources
simulation.all_molar_sources_vol()[:] = 0.0
simulation.all_thermal_sources()[:] = 0.0

final_time3 = 2000 * year
simu_time = simulation.standard_loop(
    initial_time=simu_time,
    final_time=final_time3,
    time_step_manager=tsmger,
    output_period=10 * year,
    newton=newton,
    dumper=dumper,
    # nitermax=20,
)
# simulation.postprocess(time_unit="second")
# export also the densities...
export_states(ComPASS.to_output_directory("states_2000yr"))
