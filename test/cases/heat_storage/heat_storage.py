# -*- coding: utf-8 -*-

import numpy as np

from MeshTools import TetWedgeMesh, Wedge, idarray
from MeshTools.utils import axis_extrusion

import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import TimeStepManager, Event
from ComPASS.newton import Newton
from ComPASS.linalg.factory import linear_solver
import ComPASS.mpi as mpi

from utils import hydrostatic_pressure

# -- Parameters ---------------------------------------------------------
# fmt: off
with_wells = True            # flag to test simulation with or without wells
dt_ref = 30                  # timestep in s used after closing or opening a well
layer_thickness = 1          # m
nb_layers = 9                #
reservoir_layers = [3, 4, 5] # index of reservoir layers : starts indexing at 0
radius_well = 0.1            # well radius [m]
ptop = 160 * bar             # initial reservoir pressure
Ttop = degC2K(65)            # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
k_reservoir = 3E-12          # [m²]
k_overburden = 1E-18         # [m²]
phi = 0.15                   # reservoir porosity
rho_rock = 2680.             # [kg.m-3]
thermal_conductivity = 2.5   # Thermal conductivity of solid phase kg.m.s^-3.K^-1 or W.m^-1.K^-1
heat_capacity_rock = 833.    # [J.kg-1.K-1] spe_heat_capacity_rock
gravity = 9.80665
vgradT = 0 / km              # thermal gradient degrees per km - constant to see injection effect
# fmt: on
# ------------------------------------------------------------------------

ztop = nb_layers * layer_thickness
geotherm = lambda z: Ttop + vgradT * (ztop - z)

# well ids: two wells 1 and 2 each can be producer (p) or injector (i)
idW1i, idW1p, idW2i, idW2p = range(4)

simulation = ComPASS.load_eos("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(gravity)
# Waiting function to modify rock_volumetric heat capacity as function of rock
simulation.set_rock_volumetric_heat_capacity(rho_rock * heat_capacity_rock)

vertices = np.loadtxt("doublet.nodes")
triangles = np.loadtxt("doublet.triangles")
vertices = np.column_stack([vertices, np.tile(0, vertices.shape[0])])  # to 3D
vertices, wedges = axis_extrusion(vertices, triangles, [layer_thickness] * nb_layers)
elements = [Wedge(idarray(wedge)) for wedge in wedges]
mesh = TetWedgeMesh.create(vertices, elements)


def select_dirichlet_nodes():
    xyz = simulation.global_vertices()
    x, y = [xyz[:, j] for j in range(2)]
    return (x <= x.min()) | (x >= x.max()) | (y >= y.max())


def create_well():
    make_well = lambda xy: simulation.create_vertical_well(xy, radius_well)

    def create_reversible_well(xy, idi, idp):
        I = make_well(xy)
        I.id = idi
        I.operate_on_flowrate = 0, np.inf  # mass flow rate, pressure limit (maximum)
        I.inject(degC2K(30))
        P = make_well(xy)
        P.id = idp
        P.produce()
        P.operate_on_flowrate = 0, 1 * bar  # pressure limit to avoid gaz phase
        return [I, P]

    return create_reversible_well((-600, 0), idW1i, idW1p) + create_reversible_well(
        (600, 0), idW2i, idW2p
    )


def permeability():
    zc = simulation.compute_global_cell_centers()[:, 2]
    k = np.full_like(zc, k_overburden)
    for ik in reservoir_layers:
        reservoir_layer = (
            np.fabs(zc - (ik + 0.5) * layer_thickness) < 0.5 * layer_thickness
        )
        k[reservoir_layer] = k_reservoir
    return k


simulation.init(
    mesh=mesh,
    wells=create_well,
    cell_permeability=permeability,
    cell_porosity=phi,
    cell_thermal_conductivity=thermal_conductivity,
    set_dirichlet_nodes=select_dirichlet_nodes,
)

# -- Set initial state and boundary conditions
initial_state = simulation.build_state(simulation.Context.liquid, p=ptop, T=Ttop)
simulation.all_states().set(initial_state)
dirichlet = simulation.dirichlet_node_states()
dirichlet.set(initial_state)  # will init all variables: context, states...

rho = lambda p, z: simulation.liquid_volumetric_mass_density(p, geotherm(z))
hp = hydrostatic_pressure(0, ztop, ptop, rho, 2 * nb_layers, gravity=gravity)


def set_pT_distribution(states, z):
    states.p[:] = hp(z)
    states.T[:] = geotherm(z)


set_pT_distribution(simulation.node_states(), simulation.vertices()[:, 2])
set_pT_distribution(simulation.cell_states(), simulation.compute_cell_centers()[:, 2])
set_pT_distribution(dirichlet, simulation.vertices()[:, 2])

# We set an *eager* timestep strategy to quickly increase timestep after closing/opening events
tsmger = TimeStepManager(
    initial_timestep=dt_ref,
    maximum_timestep=day,
    increase_factor=2.0,
    decrease_factor=0.5,
)

# Build well histories

well_operations = []
# wells_history in colums: t(week), Q1, Q2, T1, T2
histories = np.loadtxt("wells_history.txt", usecols=(0, 1, 2, 3, 4))
t = histories[:, 0] * week

if with_wells:

    rho_ref = simulation.liquid_volumetric_mass_density(ptop, Ttop)

    def mass_flowrate(Qv):  # m3/h to kg/s /2. => Halfdomain
        return 0.5 * (rho_ref / hour) * Qv

    Q1, Q2 = [mass_flowrate(histories[:, col]) for col in (1, 2)]

    def convert_temperature(TdegC):
        return np.where(TdegC > 0, degC2K(TdegC), 0)

    T1, T2 = [convert_temperature(histories[:, col]) for col in (3, 4)]

    # p producer, i injector
    def dispatch_flowrates(Q):
        return np.where(Q < 0, -Q, 0), np.where(Q > 0, -Q, 0)

    Q1p, Q1i = dispatch_flowrates(Q1)
    Q2p, Q2i = dispatch_flowrates(Q2)

    def filter_production(t, Q):
        result = []
        for ti, Qi in zip(t, Q):
            if len(result) == 0 or result[-1][1] != Qi:
                result.append((ti, Qi))
        return np.array(result)

    def filter_injection(t, Q, T):
        result = []
        for ti, Qi, Ti in zip(t, Q, T):
            if len(result) == 0 or result[-1][1] != Qi or result[-1][2] != Ti:
                result.append((ti, Qi, Ti))
        return np.array(result)

    fQ1p = filter_production(t, Q1p)
    well_operations.extend(simulation.well_production_history(idW1p, fQ1p))
    fQ2p = filter_production(t, Q2p)
    well_operations.extend(simulation.well_production_history(idW2p, fQ2p))
    fQ1i = filter_injection(t, Q1i, T1)
    well_operations.extend(simulation.well_injection_history(idW1i, fQ1i))
    fQ2i = filter_injection(t, Q2i, T2)
    well_operations.extend(simulation.well_injection_history(idW2i, fQ2i))

    def timestep_setter(dt):
        def setter(tick):
            if tick.latest_timestep is not None and tick.latest_timestep > dt:
                mpi.master_print(
                    f"Changing timestep from {tick.latest_timestep} to {dt}"
                )
                tsmger.current_step = dt

        return setter

    # We force a small time steps at closing/opening times
    well_on_off = set()
    for fQ in [fQ1p, fQ2p, fQ1i, fQ2i]:
        well_on_off.update(fQ[fQ[:, 1] == 0, 0])  # closing
        well_on_off.update(fQ[1:][fQ[:-1, 1] == 0, 0])  # opening
    well_operations.extend([Event(t, [timestep_setter(dt_ref)]) for t in well_on_off])

for wid in [
    idW1p,
    idW1i,
    idW2p,
    idW2i,
]:
    simulation.close_well(wid)

if mpi.is_on_master_proc:
    from time import process_time

    now = process_time()

# Construct the linear solver and newton objects outside the time loop
# to set their parameters. Here direct solving is activated
lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-5, 8, lsolver)

simulation.standard_loop(
    newton=newton,
    time_step_manager=tsmger,
    final_time=t[-1],
    output_period=30 * day,
    events=well_operations,
)

if mpi.is_on_master_proc:
    print(f"Elapsed time: {process_time()-now}")
