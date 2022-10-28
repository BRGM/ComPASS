#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

import ComPASS
from ComPASS.utils.grid import on_zmax, on_vertical_boundaries
from ComPASS.timestep_management import TimeStepManager
from ComPASS.utils.units import *
import ComPASS.mpi as mpi
from mpi4py import MPI
import os
from ComPASS._kernel import get_kernel
from ComPASS.linalg.solver import IterativeSolverSettings
from ComPASS.linalg.petsc_linear_solver import PetscDirectSolver, PetscIterativeSolver
from ComPASS.timestep_management import TimeStepManager, FixedTimeStep
from ComPASS.linalg.factory import linear_solver
from ComPASS.newton import Newton
from ComPASS.newton_coupled import NewtonCoupled
from ComPASS.utils.grid import on_vertical_boundaries
from ComPASS.simulation import AlignmentMethod

# fmt: off
# A vertical well in the middle of a regular grid
# with dh * dh squared cells, H total thickness and, nv vertical layers
L = 1000               # half-width of the reservoir
dh = 100               # horizontal cell size
H = 200                # thickness of the reservoir / height of the well
nv = 10                 # number of vertical layers
rw = 0.1               # well radius (m)
ptop = 4 * MPa
Teps = 1               # temperature variation temperature at the top of the reservoir is set to Tsat(ptop) - Teps
gravity = 9.81
omega = 0.15           # reservoir porosity
kh = 5e-14             # reservoir horizontal permeability in m^2
kvh = 0.001            # ratio of vertical permeability to horizontal permeability - low value to prevent natural convection
K = 2                  # bulk thermal conductivity in W/m/K
Qm = 200 * ton / hour
wid = 0                # well id - could be any number
ip_monitor = True
sleep = False
nprocs = MPI.COMM_WORLD.Get_size()
# fmt: on
################################################################################
class NewtonParameters:
    def __init__(self):
        #########################################################################
        # Newton and solver parameters
        # self.newtonMaxIter = 40
        self.newtonMaxIter = 400
        self.TotalMaxIter = 1000
        self.crit_resrel = 1.0e-5
        self.crit_dxmax = 1.0e-10
        #########################################################################


##################################################################################

simulation = ComPASS.load_physics("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(gravity)

nh = 2 * int(L / dh)
assert nh > 0
grid = ComPASS.Grid(
    shape=(nh, nh, nv),
    extent=(2 * L, 2 * L, H),
    origin=(-L, -L, 0),
)


def make_producer():
    well = simulation.create_vertical_well((0, 0), rw, multi_segmented=True)
    well.id = wid
    well.operate_on_flowrate = Qm, 1.0 * bar
    #    #well.operate_on_pressure =  4.0e5  , Qm
    well.produce()
    return [well]


def permeability():
    k = kh * np.eye(3, dtype="d")
    k[2, 2] *= kvh
    return k


simulation.init(
    mesh=grid,
    wells=make_producer,
    cell_porosity=omega,
    cell_permeability=permeability,
    cell_thermal_conductivity=K,
    # well_model="two_phases",
)

Ttop = simulation.Tsat(ptop) - Teps
print(f"temperature at reservoir top: {K2degC(Ttop):.1f} °C")
print(f"   delta pressure to boiling: {(simulation.Psat(Ttop) - ptop)/bar:.1f} bar")

## ------------------------------
## First step - equilibrium state
## ------------------------------

X0 = simulation.build_state(simulation.Context.liquid, p=ptop, T=Ttop)
# approximate rho
rho = simulation.liquid_molar_density(ptop, Ttop)
ztop = grid.origin[2] + grid.extent[2]


def set_states(states, z):
    states.set(X0)
    states.p[:] = ptop + rho * gravity * (ztop - z)


set_states(simulation.node_states(), simulation.vertices()[:, 2])
set_states(simulation.cell_states(), simulation.cell_centers()[:, 2])
set_states(simulation.fracture_states(), simulation.compute_fracture_centers()[:, 2])
simulation.reset_dirichlet_nodes(on_zmax(grid))
assert np.all(simulation.dirichlet_node_states().p == simulation.node_states().p)
kernel = get_kernel()
kernel.mswells_copy_states_from_reservoir()
simulation.close_well(wid)
simulation.alignment = (
    AlignmentMethod.manual
)  # simulations that have mswells need manual alignment to be set


tsmger = TimeStepManager(
    initial_timestep=1 * year,
    increase_factor=2.0,
    decrease_factor=0.2,
)

# We have a closed/open mswell, thus we need a DirectSolver and  need to specify not to use a LegacyLSolver
lsolver = linear_solver(simulation, legacy=False, direct=True)
newton = Newton(simulation, 1e-5, maxit=8, lsolver=lsolver)
#
simulation.standard_loop(
    final_time=10 * year,  # more than enough to reach pressure equilibrium
    time_step_manager=tsmger,
    newton=newton,
    no_output=True,
)


## ------------------------
## Second step - production
## ------------------------

# Dirichlet nodes will be locked to their equilibrium values
simulation.reset_dirichlet_nodes(on_vertical_boundaries(grid))

simulation.open_well(wid)
simulation.set_well_property(wid, imposed_flowrate=Qm)
# Init
kernel = get_kernel()
#########################################################
# **Computes the coupling between reservoir and mswells**#
kernel.IncCVMSWells_set_compute_coupling(True)
#########################################################
kernel.mswells_copy_states_from_reservoir()
kernel.mswells_init_edge_data()
kernel.IncPrimSecdMSWells_compute()
kernel.LoisThermoHydroMSWells_compute()
# kernel.mswells_init_leaf_data()
# kernel.VSHydroMSWells_init()  # Computes initial pressure, don't call this when reinit from file


##factor_t0 = 0.0005
factor_t = 0.5
timestep = TimeStepManager(
    initial_timestep=0.5,
    minimum_timestep=0.01,
    maximum_timestep=factor_t * day,
    increase_factor=1.2,
    decrease_factor=0.2,
)


# timestep = FixedTimeStep(step=0.5)

output_period = day
# final_time = 0.5*day
# final_time = 30 * day
final_time = 10

if mpi.is_on_master_proc:
    print("\n\n\tRunning second part...\n\n", flush=True)

# We have a closed/open mswell, thus we need a DirectSolver and  need to specify not to use a LegacyLSolver
lsolver = linear_solver(simulation, legacy=False, direct=True)
mynewtonparam = NewtonParameters()

newton = NewtonCoupled(
    simulation,
    tol=mynewtonparam.crit_resrel,
    dxmax_tol=mynewtonparam.crit_dxmax,
    maxit=mynewtonparam.newtonMaxIter,
    lsolver=lsolver,
)


simulation.standard_loop(
    initial_time=0,
    newton=newton,
    time_step_manager=timestep,
    output_period=output_period,
    final_time=final_time,
)


# assert simulation.all_states().S[:, 0].max() > 0.1
