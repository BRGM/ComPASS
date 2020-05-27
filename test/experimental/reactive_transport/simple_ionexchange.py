#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

#
# The example in this file is from the paper
# AJ Valocchi, RL Street, PV Roberts Transport of ion-exchanging solutes in groundwater:
# Chromatographic theory and field simulation Water Resources Research 17 (5),
# 1517-1527, 1981
# See Fig. 1
#

import ComPASS
from ComPASS.utils.units import *
from ComPASS.timestep_management import FixedTimeStep, TimeStepManager
import ComPASS.timestep as timestep
from ComPASS.simulation_context import SimulationContext

from scipy.optimize import newton_krylov
import numpy as np

import sys
from copy import copy


Lx = 30
nx = 100

grid = ComPASS.Grid(shape=(nx, 1, 1), extent=(Lx, Lx / nx, Lx / nx))

muf = 1  # dynamic viscosity of water
pleft, pright = 0.369 * Lx, 0  # Darcy velocity given 0.369 m / h
k_reservoir = 1  # m^2
omega_reservoir = 0.2
vdarcy = (k_reservoir / muf) * (pleft - pright) / Lx
print("Darcy velocity: ", vdarcy, "m / h")
volcell = (Lx / nx) ** 3

rhow = 1
b = 1
rhocp = 1  # g / l
K_reservoir = 0.069  # m^2 / h
rhos = 2  # g /l ?
Keq = 1

# All in mmol / l, last in mmol / g
ctot_init = 100
ctot_inj = 10
c1_init = 20
c1_inj = 9.5
cTbar = 6

on_the_left = lambda x: x <= grid.origin[0]
on_the_right = lambda x: x >= grid.origin[0] + grid.extent[0]


def select_dirichlet_nodes():
    x = simulation.global_vertices()[:, 0]
    return on_the_left(x) | on_the_right(x)


# concentrations (=temperature !) treated separately, as there are 2 species
def set_boundary_conditions():
    def set_states(states, x):
        left = on_the_left(x)
        states.p[left] = pleft
        right = on_the_right(x)
        states.p[right] = pright
        both = left | right
        states.context[both] = 1
        states.S[both] = 1
        states.C[both] = 1.0

    set_states(simulation.dirichlet_node_states(), simulation.vertices()[:, 0])


def set_states_inj(states, x, name):
    left = on_the_left(x)
    right = on_the_right(x)
    if name == "cT":
        states.T[left] = ctot_inj
        states.T[right] = ctot_init
    elif name == "c1":
        states.T[left] = c1_inj
        states.T[right] = c1_init
    else:
        raise SystemExit("wrong name %s % (name)")


def set_injection(name):
    set_states_inj(
        simulation.dirichlet_node_states(), simulation.vertices()[:, 0], name
    )


# concentrations (=temperature !) treated separately, as there are 2 species
def set_initial_values():
    def set_states(states, x):
        states.context[:] = 1
        states.p[:] = pleft + (pright - pleft) * (x - grid.origin[0]) / Lx
        states.S[:] = 1
        states.C[:] = 1.0

    set_states(simulation.all_states(), simulation.all_positions()[:, 0])


simulation = ComPASS.load_eos("linear_water")
ComPASS.set_output_directory_and_logfile(__file__)

fluid_properties = simulation.get_fluid_properties()
fluid_properties.specific_mass = rhow
fluid_properties.volumetric_heat_capacity = b
fluid_properties.dynamic_viscosity = muf
simulation.set_rock_volumetric_heat_capacity(rhocp)

simulation.init(
    mesh=grid,
    set_dirichlet_nodes=select_dirichlet_nodes,
    cell_porosity=omega_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=K_reservoir,
)

set_initial_values()
set_boundary_conditions()

#%%-- The four following functions are here just to hide
#     the fact that we use temperatures as concentrations


def retrieve_concentrations():
    # copy needed
    return np.copy(simulation.all_states().T)


def set_concentrations(C):
    simulation.all_states().T[:] = C


def set_source_term(S):
    # WARNING: only porous volume at dof not including rock volume
    porous_volume = simulation.all_Fourier_porous_volumes()
    thermal_sources = simulation.all_thermal_sources()
    thermal_sources[:] = -porous_volume / omega_reservoir * S


def clear_source_term():
    simulation.all_thermal_sources()[:] = 0


#%%---------------------------------------------------------


def make_one_timestep(t, dt, cTold, c1old):

    ts_manager = FixedTimeStep(dt)
    # ts_manager = TimeStepManager(initial_timestep=720,)
    newton = ComPASS.newton.Newton(
        simulation, 1e-5, 20, ComPASS.newton.LegacyLinearSolver(1e-6, 150)
    )  # simulation.default_Newton() #
    context = SimulationContext()
    # do cT first (linear)
    set_concentrations(cTold)
    set_injection("cT")
    clear_source_term()

    timestep.make_one_timestep(newton, ts_manager.steps(), simulation_context=context)

    cTnew = retrieve_concentrations()

    print("---------- cT --> c1 ------------------")
    # Next do c1, solve non-linear system with Newton Krylov
    def freac(c1, cT):
        return Keq * rhos * cTbar * (c1 / (cT + (Keq - 1) * c1))

    fc1old = freac(c1old, cTold)

    def fchim(cprev):
        set_concentrations(c1old)
        set_injection("c1")
        set_source_term((freac(cprev, cTnew) - fc1old) / dt)

        timestep.make_one_timestep(
            newton, ts_manager.steps(), simulation_context=context
        )
        return retrieve_concentrations()

    cinit = c1old
    c1new = newton_krylov(lambda c: c - fchim(c), cinit, method="lgmres", verbose=1)

    return cTnew, c1new


def plot_1D_concentrations(t, cT, c1):
    # Cell concentrations are last in the cT / c1 vectors
    xc = simulation.compute_cell_centers()[:, 0]
    xn = simulation.vertices()[:, 0]
    import ComPASS.utils.mpl_backends as mpl_backends

    plt = mpl_backends.import_pyplot(False)
    if plt:
        plt.clf()
        plt.subplot(211)
        plt.plot(xc, cT[-nbCells:])
        plt.xticks(np.arange(0, Lx + 0.01, step=5))
        plt.title(f"t={t:.2f}")
        plt.ylabel("Total conc")
        plt.subplot(212)
        plt.plot(xc, c1[-nbCells:])
        plt.xticks(np.arange(0, Lx + 0.01, step=5))
        plt.xlabel("x in meters")
        plt.ylabel("C1 conc")
        plt.draw()
        plt.pause(0.1)

        plt.savefig(ComPASS.to_output_directory(f"conc_{t:.2f}.png"), format="png")


t = 0
final_time = 55  # hrs
dt = 0.2
check = True  # show solution at t=25 and t=55 (cf figure 1 in above Ref.)

nbNodes = simulation.global_number_of_nodes()
nbCells = simulation.global_number_of_cells()
nbDofs = nbNodes + nbCells

cTold = np.tile(ctot_init, nbDofs)
c1old = np.tile(c1_init, nbDofs)


while t < final_time:
    print("===== Doing time ", t)

    cTnew, c1new = make_one_timestep(t, dt, cTold, c1old)
    plot_1D_concentrations(t, cTnew, c1new)
    t = t + dt
    cTold = cTnew
    c1old = c1new
    if check:
        if np.abs(t - 25) < dt / 2:
            sol25 = np.vstack(([cTold[-nbCells:], c1old[-nbCells:]]))
        if np.abs(t - 55) < dt / 2:
            sol55 = np.vstack(([cTold[-nbCells:], c1old[-nbCells:]]))

if check:
    import ComPASS.utils.mpl_backends as mpl_backends

    plt = mpl_backends.import_pyplot(False)
    if plt:
        plt.ion()
        xc = simulation.compute_cell_centers()[:, 0]
        plt.figure(2)
        plt.subplot(211)
        plt.ylim(0, 110)
        plt.plot(xc, sol25[0, :], label="T=25")
        plt.plot(xc, sol55[0, :], label="T=55")
        plt.ylabel("Total conc")
        plt.legend()
        plt.subplot(212)
        plt.ylim(0, 22)
        plt.plot(xc, sol25[1, :], label="T=25")
        plt.plot(xc, sol55[1, :], label="T=55")
        plt.xlabel("x in meters")
        plt.ylabel("C1 conc")
        plt.legend()
        plt.draw()
        plt.savefig(ComPASS.to_output_directory(f"concs_25_55.png"), format="png")
