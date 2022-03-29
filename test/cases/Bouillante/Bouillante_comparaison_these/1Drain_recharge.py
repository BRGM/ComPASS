#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#
# Cartesian grid, isolated box of 1000m in depth
# Homogeneous Neumann BC at both sides and bottom
# atm BC at the top
# start with no rain : desaturation at the top
# after 10^4 years : impose lots of rain
# to obtain liquid outflow
# gravity = 9.81

import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.linalg.factory import linear_solver
from ComPASS.newton import Newton
from ComPASS.timestep_management import TimeStepManager
from ComPASS.utils.grid import on_zmax
from ComPASS.postprocess import postprocess

Lz = 200.0
nz = 50
nx = ny = 1
dz = Lz / nz
Lx = Ly = nx * dz
Ox, Oy, Oz = 0.0, 0.0, 0.0
Topz = Oz + Lz

omega_reservoir = 0.35  # reservoir porosity
k_reservoir = 1e-12  # reservoir permeability in m^2, 1D = 10^-12 m^2
cell_thermal_cond = 3.0  # reservoir thermal conductivity : no thermal diffusion
Pporous = 1 * bar  # porous Pressure
Tatm = 300.0
Train = 300.0
Tporous = 350.0  # porous Temperature
Ttop = 320  # init atm BC temperature between Tporous and Tatm
CpRoche = 2.0e6
gravity = 9.81
pure_phase_molar_fraction = [[1.0, 0.0], [0.0, 1.0]]

simulation = ComPASS.load_eos("diphasic")
simulation.set_gravity(gravity)
simulation.set_liquid_capillary_pressure("Beaude2018")
simulation.set_atm_temperature(Tatm)
simulation.set_rain_temperature(Train)
simulation.set_rock_volumetric_heat_capacity(CpRoche)
ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(Lx, Ly, Lz),
    origin=(Ox, Oy, Oz),
)

simulation.init(
    mesh=grid,
    cell_porosity=omega_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=cell_thermal_cond,
)

fc = simulation.compute_face_centers()
simulation.set_freeflow_faces(on_zmax(grid)(fc))
is_ff = simulation.get_freeflow_nodes()  # array of bool of size n_nodes


def lininterp(z, top, gradient):
    return top + (gradient) * (z)


def set_iso_states():
    gravity = simulation.get_gravity()

    def set_states(states, z):
        states.context[:] = simulation.Context.liquid
        states.p[:] = lininterp(
            Topz - z,
            Pporous,
            gravity * 999.0,
        )  # Sg = 0 and no entry pressure if capillary pressure, so pl = pg
        states.T[:] = Tporous
        states.S[:] = [0, 1]
        states.C[:] = pure_phase_molar_fraction

    z = simulation.vertices()[:, 2]
    set_states(simulation.node_states(), z)
    z = simulation.compute_cell_centers()[:, 2]
    set_states(simulation.cell_states(), z)
    z = simulation.compute_fracture_centers()[:, 2]
    set_states(simulation.fracture_states(), z)


set_iso_states()
X_top = simulation.build_state(
    simulation.Context.diphasic_FF_liq_outflow,
    p=Pporous,
    T=Ttop,
    Cal=0.01,
)
simulation.node_states().set(is_ff, X_top)

tsmger = TimeStepManager(
    initial_timestep=10 * day,
    minimum_timestep=1,
    maximum_timestep=1.0e3 * year,
    increase_factor=1.2,
    decrease_factor=0.5,
)

initial_time = 0.0
final_time = 1.0e4 * year

# Construct the linear solver and newton objects outside the time loop
# to set their parameters. Here direct solving is activated
lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-7, 10, lsolver)

run_loop = lambda initial_time, final_time, no_output=True: simulation.standard_loop(
    reset_iteration_counter=True,
    initial_time=initial_time,
    final_time=final_time,
    newton=newton,
    time_step_manager=tsmger,
    # output_period=output_period,
    no_output=no_output,
)

# to init the domain without rain flux
simulation.set_atm_rain_flux(0.0)
# output_period = 0.1 * final_time
run_loop(initial_time, final_time, no_output=True)

# simulation.reload_snapshot(simulation.runtime.output_directory + "_init/")

simulation.set_atm_rain_flux(-3.2e-2)  # mol/m^2/s
initial_time = final_time
final_time += 20 * year
tsmger.current_step = 10.0 * day
# output_period = (final_time - initial_time) / 100.0
run_loop(initial_time, final_time, no_output=True)

simulation.postprocess()
