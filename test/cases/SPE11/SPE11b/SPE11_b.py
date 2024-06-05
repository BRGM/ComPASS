# -*- coding: utf-8 -*-

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
import GROUPS as groups

# from ComPASS.utils.salome import SalomeWrapper
# import ComPASS.mpi as mpi
# from ComPASS.linalg.petsc_linear_solver import *

ComPASS.set_output_directory_and_logfile(__file__)

""" Parameters"""
year: float = 31536.0e3  # seconds
p0: float = 211.74 * bar  # 221.74434 * bar # Pressure approximation at the top
T0: float = degC2K(70.0)  # Temperature at the bottom (zmin = 0)
pwell: float = 300 * bar  # Pressure at well1
Twell: float = degC2K(10.0)  # Temperature at the well1
pure_phase_molar_fraction = [[1.0, 0.0], [0.0, 1.0]]  # [[CO2_g, W_g], [CO2_l, W_l]]

# Wellnodes coordinates: well1(2700, 1000, 300) and well2(2700, 4000, 300)
# x1 = 2700., y1 = 50. , z1 = 300 # y1=100/2 for Spe11_b, with well_depth = 100

# The SPE11b case is in 2D, however ComPASS is only in 3D
# (the mesh contains one cell in the y direction).
# Then the well is a horizontal line instead of a point.
# The molar flux is defined over the perimeter of the well,
# thus we mutilpy by 2*pi*r*well_depth to get the molar_flux.
# We need to transform it into a volumetric flux, it is done
# by dividing by the volume of the cell containing the well.
# Also, the given flux (0.035) is a mass flux, we transform it into a molar flux using M_CO2=44.01 g/mol.
# *** WARNING ***:
# IF CHANGING THE MESH, change well_depth (up to now 100m),
well_depth: float = 100  # m (change if you change the mesh !!!)
well_radius: float = 0.15
# wells coordinates : well1, well2
wells_coord = np.array([[2700, 50, 300], [5100, 50, 700]])
# measure of the cells containing the wells m^3 (change if you change the mesh !!!)
# the correspondance between this values and the mesh is tested in check_source_vol
wells_vol = np.array([85216.3204, 137235.7283])
molar_flux = 0.035 / 44.01e-3 * 2 * 3.14 * well_radius * well_depth
well1_mol_flux = molar_flux / wells_vol[0]
well2_mol_flux = molar_flux / wells_vol[1]

k_Facies1: float = 1e-16  # domain permeability in m^2
k_Facies2: float = 1e-13  # domain permeability in m^2
k_Facies3: float = 2e-13  # domain permeability in m^2
k_Facies4: float = 5e-13  # domain permeability in m^2
k_Facies5: float = 1e-12  # domain permeability in m^2
k_Facies6: float = 2e-12  # domain permeability in m^2
k_Facies7: float = 1e-20  # domain permeability in m^2 should be 0

phi_Facies1: float = 0.1  # porosity
phi_Facies2: float = 0.2  # porosity
phi_Facies3: float = 0.2  # porosity
phi_Facies4: float = 0.2  # porosity
phi_Facies5: float = 0.25  # porosity
phi_Facies6: float = 0.35  # porosity
phi_Facies7: float = 0.0001  # 0.0001 porosity should be 0

kth_Facies1: float = 1.9  # rock heat conductivity in W/m/K
kth_Facies2: float = 1.25  # rock heat conductivity in W/m/K
kth_Facies3: float = 1.25  # rock heat conductivity in W/m/K
kth_Facies4: float = 1.25  # rock heat conductivity in W/m/K
kth_Facies5: float = 0.92  # rock heat conductivity in W/m/K
kth_Facies6: float = 0.26  # rock heat conductivity in W/m/K
kth_Facies7: float = 2.0  # rock heat conductivity in W/m/K


heat_capacity_rock: float = 850.0  # [J.kg-1.K-1] spe11b 0.85 kJ.kg-1.K-1
rho_rock_density: float = 2500.0  # [kg.m-3]

gravity: float = 9.80665

# Tgeo(x) = 70 - 0.025z
Tgrad = -0.025  # 25°C/km # Geothermal gradient

""" MESH """
# sw = SalomeWrapper(simulation)
""" Mesh import from Salome and export as paraview file """


def extract(filename, cols=None, dtype=int, sep=" "):
    filepath = Path(filename)
    assert filepath.exists()
    res = []
    if cols is None:
        append_line = lambda l: res.append([dtype(s) for s in l])
    elif type(cols) is slice:
        append_line = lambda l: res.append([dtype(s) for s in l[cols]])
    else:
        append_line = lambda l: res.append([dtype(l[k]) for k in cols])
    with filepath.open() as f:
        for line in f:
            l = line.strip().split(sep)
            append_line(l)
    return res


# *** WARNING ***:
# IF CHANGING THE MESH, change well_depth (up to now 100m),
vertices = np.array(extract("NODES.txt", dtype=float, cols=slice(1, None)))
cells = extract("PRISMS.txt", cols=slice(1, None))

# Zero based indexing
cells = [MT.idarray(cell) - 1 for cell in cells]
elements = [MT.Wedge(cell) for cell in cells]

mesh = MT.WedgeMesh.create(vertices, elements)

# MT.to_vtu(mesh, "my_wedge_mesh")

# Limits
vertices = mesh.vertices_array()
x = vertices[:, 0]
y = vertices[:, 1]
z = vertices[:, -1]
xmax, xmin = vertices[:, 0].max(), vertices[:, 0].min()
ymax, ymin = vertices[:, 1].max(), vertices[:, 1].min()
zmax, zmin = vertices[:, -1].max(), vertices[:, -1].min()
print("xmax=", xmax, "xmin=", xmin)
print("ymax=", ymax, "ymin=", ymin)
print("zmax=", zmax, "zmin=", zmin)

nodes2cell = {tuple(cell): ck for ck, cell in enumerate(cells)}

""" Groupes of volumes = Facies """
# Facies1_elements = [MT.idarray(elem) - 1 for elem in groups.Facies1]
# Facies1_Do = np.array([
#    nodes2cell[tuple(cell)] for cell in Facies1_elements
# ])
Facies1_Layer = np.array(
    [nodes2cell[tuple(MT.idarray(elem) - 1)] for elem in groups.Facies1]
)
Facies2_Layer = np.array(
    [nodes2cell[tuple(MT.idarray(elem) - 1)] for elem in groups.Facies2]
)
Facies3_Layer = np.array(
    [nodes2cell[tuple(MT.idarray(elem) - 1)] for elem in groups.Facies3]
)
Facies4_Layer = np.array(
    [nodes2cell[tuple(MT.idarray(elem) - 1)] for elem in groups.Facies4]
)
Facies5_Layer = np.array(
    [nodes2cell[tuple(MT.idarray(elem) - 1)] for elem in groups.Facies5]
)
Facies6_Layer = np.array(
    [nodes2cell[tuple(MT.idarray(elem) - 1)] for elem in groups.Facies6]
)
Facies7_Layer = np.array(
    [nodes2cell[tuple(MT.idarray(elem) - 1)] for elem in groups.Facies7]
)

Well1_cell = np.array(
    [nodes2cell[tuple(MT.idarray(elem) - 1)] for elem in groups.Well1]
)
Well2_cell = np.array(
    [nodes2cell[tuple(MT.idarray(elem) - 1)] for elem in groups.Well2]
)

Well1_node = MT.idarray(groups.Well1_node) - 1
# Well1_nodes = MT.idarray(groups.Well1_nodes) - 1 # 2 nodes
Left_nodes = MT.idarray(groups.Left_Boundary_nodes) - 1
Right_nodes = MT.idarray(groups.Right_Boundary_nodes) - 1
Top_nodes = MT.idarray(groups.Top_Boundary_nodes) - 1
Bottom_nodes = MT.idarray(groups.Bottom_Boundary_nodes) - 1
BottomRight_nodes = np.intersect1d(Right_nodes, Bottom_nodes)

Rightflag = 1
Leftflag = 2
Topflag = 3
Bottomflag = 4
Well1flag = 5
Well2flag = 6
BottomRightflag = 7

""" ComPASS simulation """
simulation = ComPASS.load_physics("diphasicCO2")
# simulation.info.activate_direct_solver = True

""" kr and Capillary functions """
# simulation.set_liquid_capillary_pressure("Beaude2018")
# set pc for rocktypes from 1 to 6
simulation.set_extendedBrooksCorey_pc_SPE11b()
# kr = S^2
# from data.kr import kr_functions
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
    nodeflags[Top_nodes] = Topflag
    nodeflags[Bottom_nodes] = Bottomflag
    nodeflags[BottomRight_nodes] = BottomRightflag
    # nodeflags[Well1_nodes] = Well1flag
    # nodeflags[Well2_nodes] = Well2flag

    cellflags = simulation.global_cellflags()
    cellflags[Well1_cell] = Well1flag
    cellflags[Well2_cell] = Well2flag


def set_global_rocktype():
    cellrocktype = simulation.global_cell_rocktypes()
    cellrocktype[Facies1_Layer] = 1
    cellrocktype[Facies2_Layer] = 2
    cellrocktype[Facies3_Layer] = 3
    cellrocktype[Facies4_Layer] = 4
    cellrocktype[Facies5_Layer] = 5
    cellrocktype[Facies6_Layer] = 6
    cellrocktype[Facies7_Layer] = 0  # no pc


""" Functions to define permeabilities: Units (Kx=Ky et Kz=0.1*Kx) """


def subdomains_permeability(display=False):
    k_cells = np.zeros((mesh.nb_cells, 3, 3))
    k_cells[Facies1_Layer] = k_Facies1 * np.array([[1, 0, 0], [0, 1, 0], [0, 0, 0.1]])
    k_cells[Facies2_Layer] = k_Facies2 * np.array([[1, 0, 0], [0, 1, 0], [0, 0, 0.1]])
    k_cells[Facies3_Layer] = k_Facies3 * np.array([[1, 0, 0], [0, 1, 0], [0, 0, 0.1]])
    k_cells[Facies4_Layer] = k_Facies4 * np.array([[1, 0, 0], [0, 1, 0], [0, 0, 0.1]])
    k_cells[Facies5_Layer] = k_Facies5 * np.array([[1, 0, 0], [0, 1, 0], [0, 0, 0.1]])
    k_cells[Facies6_Layer] = k_Facies6 * np.array([[1, 0, 0], [0, 1, 0], [0, 0, 0.1]])
    k_cells[Facies7_Layer] = k_Facies7 * np.array([[1, 0, 0], [0, 1, 0], [0, 0, 0.1]])
    if display:
        MT.to_vtu(mesh, "k_cells.vtu", celldata={"perm": k_cells})
    return k_cells


def subdomains_porosity(display=True):
    phi_cells = np.zeros(mesh.nb_cells)
    phi_cells[Facies1_Layer] = phi_Facies1
    phi_cells[Facies2_Layer] = phi_Facies2
    phi_cells[Facies3_Layer] = phi_Facies3
    phi_cells[Facies4_Layer] = phi_Facies4
    phi_cells[Facies5_Layer] = phi_Facies5
    phi_cells[Facies6_Layer] = phi_Facies6
    phi_cells[Facies7_Layer] = phi_Facies7
    if display:
        MT.to_vtu(mesh, "phi_cells.vtu", celldata={"phi": phi_cells})
    return phi_cells


def subdomains_conductivity(display=False):
    kth_cells = np.zeros(mesh.nb_cells)
    kth_cells[Facies1_Layer] = kth_Facies1
    kth_cells[Facies2_Layer] = kth_Facies2
    kth_cells[Facies3_Layer] = kth_Facies3
    kth_cells[Facies4_Layer] = kth_Facies4
    kth_cells[Facies5_Layer] = kth_Facies5
    kth_cells[Facies6_Layer] = kth_Facies6
    kth_cells[Facies7_Layer] = kth_Facies7
    if display:
        MT.to_vtu(mesh, "kth_cells.vtu", celldata={"kth": kth_cells})
    return kth_cells


""" Functions to define Dirichlet boundary condition """


def P_dirichlet_nodes():
    where = np.zeros(mesh.nb_vertices, dtype=bool)
    # where[Well1_node] = True
    # where[Top_nodes] = True
    # where[Bottom_nodes] = True
    where[BottomRight_nodes] = True
    return where


def Temp_dirichlet_nodes():
    where = np.zeros(mesh.nb_vertices, dtype=bool)
    # where[Well1_node] = True
    where[Top_nodes] = True
    where[Bottom_nodes] = True
    # where[BottomRight_nodes] = True
    return where


# # We must always create wells to have the same mesh parts
# def create_wells():
#     def _well_from_nodes(nodes, well_id):
#         z = simulation.global_vertices()[nodes, 2]
#         well = simulation.create_single_branch_well(
#             nodes[np.argsort(z)[::-1]], well_radius
#         )
#         well.id = well_id
#         return well

#     producer_id : int = 3
#     producer = _well_from_nodes(Well1_nodes, producer_id)
#     producer.operate_on_pressure = 300 * bar, 1e1*ton / hour
#     producer.produce()

#     return [producer]


def reset_well_sources(w_coord, eps=1e-13):
    # get the objects positions where the volumetric flux is non zero
    non_zero_sources = np.where(np.abs(simulation.all_thermal_sources()) > eps)[0]
    source_pos = simulation.all_positions()[non_zero_sources]
    # get the objects close to the well
    close_to_well = non_zero_sources[np.all(np.abs(source_pos - w_coord) < 60, axis=1)]
    assert len(close_to_well) <= 7
    # reset the volumetric sources
    simulation.all_thermal_sources()[close_to_well] = 0.0
    simulation.all_molar_sources_vol()[close_to_well] = 0.0


def check_source_vol(nb_wells, w_coord, w_cell_measures, w_depth, eps=1e-13):
    def triangle_surface(all_nodes):
        # Heron's Formula
        l0 = sqrt(np.sum((all_nodes[1] - all_nodes[0]) ** 2))
        l1 = sqrt(np.sum((all_nodes[2] - all_nodes[0]) ** 2))
        l2 = sqrt(np.sum((all_nodes[2] - all_nodes[1]) ** 2))
        p = (l0 + l1 + l2) / 2.0
        return sqrt(p * (p - l0) * (p - l1) * (p - l2))

    mol_source = np.where(np.abs(simulation.all_molar_sources_vol()) > eps)[0]
    therm_source_bool = np.abs(simulation.all_thermal_sources()) > eps
    therm_source = np.where(therm_source_bool)[0]
    assert np.all(mol_source == therm_source)
    # 3D tetrahedral mesh (tetra extruded),
    # the volumetric source is splitted over the cell and
    # the cell nodes, depending on the rocktypes
    assert len(mol_source) <= 7 * nb_wells

    # check with the objects positions
    source_pos = simulation.all_positions()[therm_source_bool]
    # check the well depth is coherent with the given one
    if source_pos.size :
        assert abs(np.max(source_pos[:, 1]) - np.min(source_pos[:, 1]) - w_depth) < eps

    # check that 7 objects (cell and 6 nodes) max are closed to each well,
    # and all objects from total_source_i are closed to any well
    close_to_well = np.array([], int)
    for i in range(nb_wells):
        close_to_i = therm_source[np.all(np.abs(source_pos - w_coord[i]) < 60, axis=1)]
        close_to_well = np.append(close_to_well, close_to_i)
        assert len(close_to_i) <= 7
        # check the volumes are coherent with the one given
        # get one base triangle of the tetra
        well_nodes = simulation.all_positions()[close_to_i]
        tri_nodes = well_nodes[np.abs(well_nodes[:, 1]) < eps]
        # if some nodes are not in the same rocktype, the distribution does
        # not include these nodes. Then it is not possible anymore to find
        # the cell node using the non-zero values of the distribution
        if len(tri_nodes) == 3:
            surf = triangle_surface(tri_nodes)
            assert (
                abs(w_cell_measures[i] - surf * w_depth) / abs(w_cell_measures[i])
                < 1.0e-6
            )

    # check that all objects are closed to any well
    close_to_well.sort()
    assert np.all(close_to_well == mol_source)


def cell_molar_source():
    # centers = simulation.compute_global_cell_centers()
    wells = np.zeros((mesh.nb_cells, 2))
    wells[Well1_cell, 0] = well1_mol_flux
    wells[Well2_cell, 0] = well2_mol_flux  # CO2 flux in well 2
    return wells


# def cell_heat_source():
#     wells = np.zeros(mesh.nb_cells)
#     # energy inflow is approximated using pwell, Twell
#     gas_mol_enth = simulation.cpp_gas_molar_enthalpy(pwell, Twell, [1, 0])
#     wells[Well1_cell] = well1_mol_flux * gas_mol_enth
#     wells[Well2_cell] = well2_mol_flux * gas_mol_enth
#     return wells


dumper = DumperSPE11(simulation)

simulation.init(
    mesh=mesh,
    # it seems that the well has very little impact on the results,
    # but slower convergence
    # wells = create_wells,
    cell_permeability=subdomains_permeability,
    cell_porosity=subdomains_porosity,
    cell_thermal_conductivity=subdomains_conductivity,
    set_temperature_dirichlet_nodes=Temp_dirichlet_nodes,
    set_pressure_dirichlet_nodes=P_dirichlet_nodes,
    # init molar and energy sources (energy is always necessary with molar
    # otherwise the mater enters at 0K), the energy is set later
    # will be set to 0 for the initial phase then reactivated
    cell_molar_sources=cell_molar_source,
    # carefull, do not use cell_heat_source in init because the distribution
    # depends on the rocktype then the molar and heat distribution could be different
    # cell_heat_source=cell_heat_source,
    set_global_flags=set_global_flags,
    # init the rocktypes (for the pc)
    set_global_rocktype=set_global_rocktype,
)

""" Initial conditions """
# compute the thermal source with all_molar_sources_vol and molar enthalpy
# carefull, do not use cell_heat_source in init because the distribution
# depends on the rocktype then the molar and heat distribution could be different
gas_mol_enth = simulation.cpp_gas_molar_enthalpy(pwell, Twell, [1, 0])
simulation.all_thermal_sources()[:] = (
    simulation.all_molar_sources_vol()[:, 0] * gas_mol_enth
)
# check some coherencies
check_source_vol(2, wells_coord, wells_vol, well_depth)
# save values then set to 0 the molar and energy sources
# copy is important, are pointers
final_molar_sources = simulation.all_molar_sources_vol().copy()
final_heat_source = simulation.all_thermal_sources().copy()
# reset the molar and heat sources to begin without sources
simulation.all_molar_sources_vol()[:] = 0.0
simulation.all_thermal_sources()[:] = 0.0


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


export_states("states_0yr")


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

simulation.postprocess()

""" Export values as vtufile"""
export_states("states_1e3yr")


""" Set Dirichlet nodes at top and bottom boundaries with initial states """
nodeflags = simulation.nodeflags()
dirichlet_nodes = np.zeros(len(nodeflags), dtype=bool)
dirichlet_nodes[nodeflags == Topflag] = True
dirichlet_nodes[nodeflags == Bottomflag] = True

simulation.reset_dirichlet_nodes(dirichlet_nodes)


# reactivate the CO2 and heat sources for both wells
simulation.all_molar_sources_vol()[:] = final_molar_sources
simulation.all_thermal_sources()[:] = final_heat_source
# deactivate Well2 (open later), need to reset the well 2 cell and
# its neighbouring nodes (VAG distribution of the volumetric flux)
reset_well_sources(wells_coord[1])

""" Update Solver parameters """
lsolver = linear_solver(simulation, direct=False)
newton = Newton(simulation, 1e-5, 50, lsolver, maxit_before_decrease=20)
tsmger.current_step = 0.01 * year
# tsmger.maximum_timestep = 0.05 * year
tsmger.increase_factor = 1.2
tsmger.decrease_factor = 0.3

final_time1 = 1025 * year
simu_time = simulation.standard_loop(
    initial_time=simu_time,
    final_time=final_time1,
    time_step_manager=tsmger,
    output_period=year,
    newton=newton,
    dumper=dumper,
)
simulation.postprocess()
# export also the densities...
export_states("states_1025yr")

# reactivate the CO2 and heat sources for the 2 wells
simulation.all_molar_sources_vol()[:] = final_molar_sources
simulation.all_thermal_sources()[:] = final_heat_source
tsmger.current_step = 0.01 * year

final_time2 = 1050 * year
simu_time = simulation.standard_loop(
    initial_time=simu_time,
    final_time=final_time2,
    time_step_manager=tsmger,
    output_period=year,
    newton=newton,
    dumper=dumper,
)
simulation.postprocess()
# export also the densities...
export_states("states_1050yr")

# reset the molar and heat sources to end without sources
simulation.all_molar_sources_vol()[:] = 0.0
simulation.all_thermal_sources()[:] = 0.0

# final_time3 = 2000 * year
# simu_time = simulation.standard_loop(
#     initial_time=simu_time,
#     final_time=final_time3,
#     time_step_manager=tsmger,
#     output_period=50 * year,
#     newton=newton,
#     dumper=dumper,
#     # nitermax=20,
# )
# simulation.postprocess()
# # export also the densities...
# export_states("states_2000yr")
