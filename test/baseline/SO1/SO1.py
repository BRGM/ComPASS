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
from ComPASS.timeloops import standard_loop
from ComPASS.timestep_management import TimeStepManager
from ComPASS.linalg.factory import linear_solver
from ComPASS.newton import Newton

p0 = 1.0 * bar
T0_degC = 5.0
T0 = degC2K(T0_degC)
p_outlet = p0
T_injection_degC = 33.0
T_injection = degC2K(T_injection_degC)
permeability = 1e-12  # m^2
porosity = 0.2
flowrate = 1.0e-4  # m^3/s
final_time = 1e8  # s

rhow = 1000  # kg/m3
Cw = 4200.0  # J/K/kg
Kw = 0.6  # W/m/K
rhor = 2600.0  # kg/m3
Cr = 800.0  # J/K/kg
Kr = 2.0  # W/m/K

Keq = (1 - porosity) * Kr + porosity * Kw

# grid specifications
extent = Lx, Ly, Lz = 1000, 10, 10
dx, dy, dz = 0.5, 10, 10
shape = nx, ny, nz = int(Lx / dx), int(Ly / dy), int(Lz / dz)
origin = (0, -0.5 * Ly, -0.5 * Lz)

simulation = ComPASS.load_physics("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_rock_volumetric_heat_capacity(rhor * Cr)

grid = ComPASS.Grid(shape=shape, extent=extent, origin=origin)


def outlet_nodes():
    return simulation.global_vertices()[:, 0] >= Lx


simulation.init(
    mesh=grid,
    cell_permeability=permeability,
    cell_porosity=porosity,
    cell_thermal_conductivity=Keq,
    # FIXME: should be an outflow condition cf. issue https://gitlab.inria.fr/charms/ComPASS/issues/129
    set_dirichlet_nodes=outlet_nodes,
)

initial_state = simulation.build_state(simulation.Context.liquid, p=p0, T=T0)
simulation.all_states().set(initial_state)
simulation.dirichlet_node_states().set(initial_state)


def set_boundary_flux():
    Neumann = ComPASS.NeumannBC()
    specific_massflux = (
        flowrate * simulation.liquid_molar_density(p0, T_injection) / (Ly * Lz)
    )
    Neumann.molar_flux[:] = specific_massflux
    # energy inflow is approximated using p0
    Neumann.heat_flux = specific_massflux * simulation.liquid_molar_enthalpy(
        p0, T_injection
    )
    # print('flux',  flowrate * simulation.liquid_molar_density(p0, T_injection) * ComPASS.liquid_molar_enthalpy(p0, T_injection))
    simulation.set_Neumann_faces(simulation.face_centers()[:, 0] <= 0, Neumann)


set_boundary_flux()

output_period = 0.1 * final_time


# we will store in this list the cell temperature at different times
T = []


def store_T(iteration, time):
    # the copy is important here not to have only a view of the latest array
    # K2degC will generate one
    T.append((time, K2degC(simulation.cell_states().T)))


# https://charms.gitlabpages.inria.fr/ComPASS/python_reference/ComPASS.html#ComPASS.timestep_management.TimeStepManager
ts = TimeStepManager(initial_timestep=1 * hour, maximum_timestep=0.2 * output_period)

# Construct the linear solver and newton objects outside the time loop
# to set their parameters. Here direct solving is activated
lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-5, 8, lsolver)

standard_loop(
    simulation,
    newton=newton,
    final_time=final_time,
    output_period=output_period,
    output_callbacks=[store_T],
    time_step_manager=ts,
)

# save table with collected temperatures
with open("SO1-T.csv", "w") as f:
    cell_centers = simulation.cell_centers()
    print(";", " ; ".join([str(xi) for xi in cell_centers[:, 0]]), file=f)
    for output in T:
        print(output[0], ";", " ; ".join([str(theta) for theta in output[1]]), file=f)

# build analytical solution
from SO1analytical import build_solution

solution = build_solution(
    Qinj=flowrate,
    Syz=Ly * Lz,
    porosity=0.2,
    permeability=permeability / 1e12,  # in darcies
    rhow=rhow,
    Cw=Cw,
    Kw=Kw,
    rhor=rhor,
    Cr=Cr,
    Kr=Kr,
    Tres=T0_degC,
    Tinj=T_injection_degC,
)
# extract final result from computation
last_time, last_T = T[-1]

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("No figure is drawn as matpplotlib is not installed.")
    raise
else:
    plt.clf()
    plt.plot(
        cell_centers[:, 0],
        solution(last_time, cell_centers[:, 0]),
        "-k",
        label="analytical",
    )
    # plot 1 out of 50 points
    plt.plot(cell_centers[:, 0][::50], last_T[::50], "xk", label="ComPASS")
    plt.legend()
    plt.xlabel("distance (m)")
    plt.ylabel("temperature (deg C)")
    plt.grid(True)
    plt.savefig("SO1-comparison.png")
    plt.clf()
    for output in T:
        plt.plot(cell_centers[:, 0], output[1], label="t=%.1f y" % (output[0] / year))
    plt.legend()
    plt.xlabel("distance (m)")
    plt.ylabel("temperature (deg C)")
    plt.grid(True)
    plt.savefig("SO1-fronts.png")
