# -*- coding: utf-8 -*-

import numpy as np

import ComPASS
import ComPASS.mpi as mpi
from ComPASS.utils.units import *
from ComPASS.timeloops import TimeStepManager
from ComPASS.linalg.petsc_linear_solver import *
from ComPASS.utils.salome import SalomeWrapper

ComPASS.set_output_directory_and_logfile(__file__)

""" Parameters"""
p0: float = 1 * bar
T0: float = degC2K(20.0)

k_matrix: float = 1e-19  # domain permeability in m^2
phi_matrix: float = 0.1  # 0.2  # domain porosity
phi_fracture: float = 0.2  # 0.5 ## fracture porosity
lambda_th_water: float = 0.65  # Thermal conductivity of liquid phase kg.m.s^-3.K^-1
lambda_th_rock: float = 2.0  # Thermal conductivity of solid phase kg.m.s^-3.K^-1
heat_capacity_rock: float = 1000.0  # [J.kg-1.K-1] spe_heat_capacity_rock
rho_rock_density: float = 2600.0  # [kg.m-3]

gravity: float = 9.80665

fracture_thickness: float = 50.0
zmax = 0.0
zmin = -1000.0

Tgrad = 0.1  # 30°C/km à 65°C/km # Geothermal gradient


""" ComPASS simulation """
simulation = ComPASS.load_physics("water2ph")
# simulation.info.activate_direct_solver = True

sw = SalomeWrapper(simulation)

""" global parameters """
simulation.set_gravity(gravity)
simulation.set_fracture_thickness(fracture_thickness)
simulation.set_rock_volumetric_heat_capacity(rho_rock_density * heat_capacity_rock)
bulk_thermal_conductivity = (
    lambda_th_rock * (1 - phi_matrix) + lambda_th_water * phi_matrix
)

if mpi.is_on_master_proc:
    vertices = sw.mesh.vertices_array()
    zmax, zmin = vertices[:, -1].max(), vertices[:, -1].min()

""" Functions to define permeabilities: Units and Faults """


def subdomains_permeability(display=False):
    k_cells = np.full(sw.info.mesh.nb_cells, k_matrix)
    k_cells[sw.info.unit_socle.cells] = 1e-19
    k_cells[sw.info.unit_middle.cells] = 1e-12
    k_cells[sw.info.unit_top.cells] = 1e-16
    if display:
        MT.to_vtu(mesh, "k_cells.vtu", celldata={"perm": k_cells})
    return k_cells


def subdomains_thconductivity(display=False):
    th_cells = np.full(sw.info.mesh.nb_cells, 3.0)
    th_cells[sw.info.unit_middle.cells] = 2.75
    th_cells[sw.info.unit_top.cells] = 1.0
    if display:
        MT.to_vtu(mesh, "th_cells.vtu", celldata={"th_cond": th_cells})
    return th_cells


def fracture_sets_permeability(display=False):
    k_frac = np.zeros(sw.mesh.nb_faces)
    # k_frac[:] = 2e-14
    k_frac[sw.info.Fault_socle.faces] = 1e-12
    k_frac[sw.info.Fault_middle.faces] = 5e-13
    k_frac[sw.info.Fault_top.faces] = 1e-15
    return k_frac[sw.info.Fault.faces]


""" Dirichlet Temperature at the top and Bottom if no neumman at the bottom boundary """


def topbottom_dirichlet_nodes():
    where = np.zeros(sw.mesh.nb_vertices, dtype=bool)
    where[sw.info.Top.nodes] = True
    where[sw.info.Bottom.nodes] = True
    return where


""" Dirichlet Pressure at the top if no neumman at the bottom boundary """


def top_dirichlet_nodes():
    where = np.zeros(sw.mesh.nb_vertices, dtype=bool)
    where[sw.info.Top.nodes] = True
    return where


simulation.init(
    mesh=sw.mesh,
    cell_permeability=subdomains_permeability,
    cell_porosity=phi_matrix,
    cell_thermal_conductivity=lambda_th_rock,  # subdomains_thconductivity,
    fracture_thermal_conductivity=lambda_th_water,
    fracture_faces=lambda: sw.info.Fault.faces,
    fracture_permeability=fracture_sets_permeability,
    fracture_porosity=phi_fracture,
    set_temperature_dirichlet_nodes=topbottom_dirichlet_nodes,
    set_pressure_dirichlet_nodes=top_dirichlet_nodes,
    # set_dirichlet_nodes=set_dirichlet_nodes, # to activate if neumann
    set_global_flags=sw.flags_setter,
)

sw.rebuild_locally()


""" Set global initial conditions: P0, TO = f(z)"""
""" Function: linear gradient useful for Pre and Temp"""


def lininterp(z, top, gradient):
    return top + (gradient) * (z)


def set_initial_states():
    gravity = simulation.get_gravity()

    def set_states(states, z):
        states.context[:] = simulation.Context.liquid
        states.p[:] = lininterp(
            zmax - z,
            p0,
            gravity * simulation.liquid_molar_density(p0, T0),
        )
        states.T[:] = lininterp(zmax - z, T0, Tgrad)
        states.S[:] = [0, 1]
        states.C[:] = 1

    z = simulation.vertices()[:, 2]
    set_states(simulation.node_states(), z)
    set_states(simulation.dirichlet_node_states(), z)
    z = simulation.compute_cell_centers()[:, 2]
    set_states(simulation.cell_states(), z)
    z = simulation.compute_fracture_centers()[:, 2]
    set_states(simulation.fracture_states(), z)


set_initial_states()

""" Simulation """

nitermax = 50  # We know it must be 51 iterations
final_time = 1e3 * year

tsmger = TimeStepManager(
    initial_timestep=1 * day,  # 10
    increase_factor=2,
    decrease_factor=1 / 1.5,
    karma_threshold=10,
    karma_reset=False,
    maximum_timestep=0.1 * final_time,
)

t = simulation.standard_loop(
    final_time=final_time,
    time_step_manager=tsmger,
    no_output=True,
    nitermax=nitermax,
)

assert t == final_time
