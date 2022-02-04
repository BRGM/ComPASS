#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#
# Cartesian grid, box of 800m in depth
# with 1 vertical fracture at Lx/2.
# Homogeneous Neumann BC at both sides and at the bottom
# atm BC at the top.
# Imposed flux at the bottom of the fracture.

import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop, TimeStepManager
from ComPASS.newton import Newton
from ComPASS.linalg.factory import linear_solver
import ComPASS.io.mesh as io
from ComPASS.utils.grid import on_zmax

pure_phase_molar_fraction = [[1.0, 0.0], [0.0, 1.0]]
p0 = 1.0 * bar  # initial reservoir pressure
T0 = degC2K(
    20.0
)  # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
qmass = 1e-1  #
Tbottom = degC2K(350.0)  # temperature influx
Tatm = 300.0
k_matrix = 1e-15  # domain permeability in m^2
phi_matrix = 0.15  # domain porosity
k_fracture = 1e-12  # fracture permeability in m^2
phi_fracture = 0.3  # fracture porosity
thermal_cond = 2.0  # bulk thermal conductivity in W/m/K
CpRoche = 2.0e6

H = 800.0  # domain height
ztop = 0.0
nH = 40  # discretization
nx, ny, nz = int(1.2 * nH), 1, nH
Lx, Ly, Lz = 1.2 * H, 0.1 * H, H


simulation = ComPASS.load_eos("diphasic")
simulation.set_liquid_capillary_pressure("Beaude2018")
simulation.set_atm_temperature(Tatm)
ptop = simulation.get_atm_pressure()
simulation.set_atm_rain_flux(0.0)
simulation.set_rock_volumetric_heat_capacity(CpRoche)
ComPASS.set_output_directory_and_logfile(__file__)

# thermodynamic functions are only available once the eos is loaded
pbottom = simulation.get_gravity() * H * 1000.0
hbottom = simulation.liquid_molar_enthalpy(pbottom, Tbottom, pure_phase_molar_fraction)

grid = ComPASS.Grid(
    shape=(nx, ny, nz), extent=(Lx, Ly, Lz), origin=(-0.5 * Lx, -0.5 * Ly, ztop - Lz)
)


def select_fractures():
    centers = simulation.compute_global_face_centers()
    xc = centers[:, 0]
    # zc = centers[:, -1]
    return xc == 0  # & (zc > -0.5 * H)


simulation.init(
    mesh=grid,
    cell_permeability=k_matrix,
    cell_porosity=phi_matrix,
    cell_thermal_conductivity=thermal_cond,
    fracture_faces=select_fractures,
    fracture_permeability=k_fracture,
    fracture_porosity=phi_fracture,
    fracture_thermal_conductivity=thermal_cond,
)


def set_iso_pressure():
    gravity = simulation.get_gravity()

    def lininterp(z, top, gradient):
        return top + (gradient) * (z)

    def set_states(states, z):
        states.p[:] = lininterp(
            ztop - z,
            p0,
            gravity * 1000.0,
        )

    z = simulation.vertices()[:, 2]
    set_states(simulation.node_states(), z)
    z = simulation.compute_cell_centers()[:, 2]
    set_states(simulation.cell_states(), z)
    z = simulation.compute_fracture_centers()[:, 2]
    set_states(simulation.fracture_states(), z)


X0 = simulation.build_state(simulation.Context.liquid, p=p0, T=T0)
simulation.all_states().set(X0)
set_iso_pressure()  # correct the pressure to init with isostatic pressure


def set_boundary_fluxes():
    Neumann = ComPASS.NeumannBC()
    Neumann.molar_flux[:] = qmass
    Neumann.heat_flux = qmass * hbottom
    face_centers = simulation.face_centers()
    bottom_fracture_edges = simulation.find_fracture_edges(face_centers[:, -1] <= -H)
    simulation.set_Neumann_fracture_edges(bottom_fracture_edges, Neumann)


set_boundary_fluxes()

# select freeflow faces
fc = simulation.compute_face_centers()
simulation.set_freeflow_faces(on_zmax(grid)(fc))
# read freeflow nodes
is_ff = simulation.get_freeflow_nodes()  # array of bool of size n_nodes

X_top = simulation.build_state(
    simulation.Context.diphasic_FF_liq_outflow,
    p=p0,
    T=T0,
)
simulation.node_states().set(is_ff, X_top)


def export_initial_states():
    node_states = simulation.node_states()
    cell_states = simulation.cell_states()
    petrophysics = simulation.petrophysics()

    pointdata = {
        # "dirichlet pressure": simulation.pressure_dirichlet_values(),
        # "dirichlet temperature": K2degC(simulation.temperature_dirichlet_values()),
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
    io.write_mesh(
        simulation,
        f"{simulation.runtime.output_directory}/initial_states",
        pointdata=pointdata,
        celldata=celldata,
    )


# export_initial_states()

final_time = 100 * year
output_period = 0.05 * final_time
timestep = TimeStepManager(
    initial_timestep=1.0 * day,
    minimum_timestep=1e-3,
    maximum_timestep=50 * year,
    increase_factor=1.4,
    decrease_factor=0.2,
)

# Construct the linear solver and newton objects outside the time loop
# to set their parameters.
lsolver = linear_solver(simulation, tolerance=1e-7, legacy=False, direct=False)
# linear solver tolerance < newton tolerance
newton = Newton(simulation, 1e-6, 13, lsolver)

standard_loop(
    simulation,
    final_time=final_time,
    time_step_manager=timestep,
    output_period=output_period,
    newton=newton,
)
