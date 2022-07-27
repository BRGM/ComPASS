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
import ComPASS.dump_wells as dw
import ComPASS.mpi as mpi
from mpi4py import MPI
import os
from ComPASS._kernel import get_kernel

# from ComPASS.linalg.petsc_linear_solver_mswell import PetscLinearSystemMSWell
# from ComPASS.linalg.solver import IterativeSolverSettings
# from ComPASS.linalg.petsc_linear_solver import PetscDirectSolver, PetscIterativeSolver
# from ComPASS.distributed_system_mswells import DistributedSystemMSWells
# from ComPASS.ghosts.synchronizer import Synchronizer

################################################################################################
# A vertical well with 2x2 grid basis and nv horizontal layers over H thickness
ds = 100  # horizontal cell size
H = 1000  # height of the well
nv = 10  #  number of vertical layers
hns = 1  # half the number of cells along basis side
rw = 0.05  # well radius
ptop = 5.0e5  # pressure at the top of the reservoir, 10*MPa one phase. 1*MPa for two phases
Ttop = 350  # , Kelvin, temperature at the top of the reservoir
vgradT = 170 / km  # degrees per km - to reach 300 degC at the bottom
geotherm = lambda zeta: Ttop
gravity = 9.81
epsilon = 0.0001
omega_reservoir = 0.15  # reservoir porosity
k_reservoir = (
    1e-16  # reservoir permeability in m^2 - low value to limit reservoir convection
)
K_reservoir = 2  # bulk thermal conductivity in W/m/K
QwOut = 15.0  # Maximum imposed flow rate for producer
PwOut = 5.0e5  # Minmum  imposed pressure for producer
dt = 0.5  # initial timestep
# dt = 10.0  # initial timestep
verbose = True
water2phase = True
ip_monitor = False
sleep = False
nprocs = MPI.COMM_WORLD.Get_size()
dummy_swell_id = 0
mswell_id = dummy_swell_id + 1  # well id - could be any number
################################################################################
class NewtonParameters:
    def __init__(self, dt):
        #########################################################################
        # Newton and solver parameters
        self.t0 = 0.0
        self.tf = 2000.0  # 3000.0  # 3000.0
        self.dtmax = 40.0
        if water2phase:
            self.newtonMaxIter = 40
        else:
            self.newtonMaxIter = 100
        self.TotalMaxIter = 2000
        self.crit_resrel = 1.0e-9
        self.crit_dxmax = 1.0e-10
        #########################################################################
        self.dt = dt


########################################################################
## Postprocessing
if mpi.is_on_master_proc:
    f_reswell = open("reswell_c.val", "wt")

if water2phase:
    simulation = ComPASS.load_eos("water2ph")
else:
    simulation = ComPASS.load_eos("immiscible2ph")

ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(gravity)
########################################################################


def hydrostatic_pressure(zbottom, ztop, nz):
    assert zbottom < ztop
    nbsteps = 1
    z = np.linspace(zbottom, ztop, nz)[::-1]  # from top to bottom
    rho = simulation.liquid_molar_density
    C = np.array([1])
    p = ptop
    pressures = [p]
    for zbot, ztop in zip(z[1:], z[:-1]):
        zeta = ztop
        dz = (ztop - zbot) / nbsteps
        for _ in range(nbsteps):
            if water2phase:
                p += gravity * rho(ptop, Ttop) * dz  # water2phase interface
            else:
                p += gravity * rho(ptop, Ttop, C) * dz  # immisible2ph interface
            zeta -= dz
        pressures.append(p)
    pressures = np.asarray(pressures)
    return lambda zeta: np.interp(zeta, z[::-1], pressures[::-1])


slim = hns * ds
ns = 2 * hns
grid = ComPASS.Grid(
    shape=(ns, ns, nv),
    extent=(ns * ds, ns * ds, H),
    origin=(-slim, -slim, 0),
)


def create_well(multi_flag):
    return simulation.create_vertical_well(
        (0, 0), rw, multi_segmented=multi_flag, zmin=0
    )


def make_producers():
    dummy_swell = create_well(False)
    dummy_swell.id = dummy_swell_id
    dummy_swell.operate_on_flowrate = QwOut, PwOut
    dummy_swell.produce()

    mswell = create_well(True)
    mswell.id = mswell_id

    if ip_monitor:
        mswell.operate_on_flowrate = QwOut, PwOut  # water2phase, IP-Monitor
    else:
        mswell.operate_on_pressure = (
            PwOut,
            100000,
        )  # immiscible, water2phase no monitoring

    mswell.produce()
    return [dummy_swell, mswell]


def select_dirichlet_nodes():
    vertices = simulation.global_vertices()
    x, y = vertices[:, 0], vertices[:, 1]
    return (x <= -slim) | (x >= slim) | (y <= -slim) | (y >= slim)


pid = os.getpid()
print("MPI rank:", mpi.proc_rank)
print("Process id:", pid)
if sleep:
    print("Sleeping...")
    kernel = get_kernel()
    kernel.init_wait_for_debug(True)


simulation.init(
    mesh=grid,
    wells=make_producers,
    set_dirichlet_nodes=select_dirichlet_nodes,
    cell_porosity=omega_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=K_reservoir,
)

# -- Set initial state and boundary conditions
hp = hydrostatic_pressure(0, H, 2 * nv + 1)
initial_state = simulation.build_state(simulation.Context.liquid, p=ptop, T=Ttop)
simulation.all_states().set(initial_state)
# simulation.mswell_node_states().set(initial_state)
dirichlet = simulation.dirichlet_node_states()
dirichlet.set(initial_state)  # will init all variables: context, states...


def set_pT_distribution(states, z):
    states.p[:] = hp(z)
    states.T[:] = Ttop  # simulation.Tsat(states.p[:])# for liquid_context
    # states.T[:] = simulation.Tsat(states.p[:])  # for gas_liquid_context
    states.accumulation[:] = 0


set_pT_distribution(simulation.node_states(), simulation.vertices()[:, 2])
set_pT_distribution(simulation.cell_states(), simulation.compute_cell_centers()[:, 2])
set_pT_distribution(dirichlet, simulation.vertices()[:, 2])

##########################################################################################
class MSWellsNewton:
    def __init__(self, newton_parameters):
        self.kernel = get_kernel()
        self.residuals = self.kernel.ResidualsMSWells()
        # self.increments = self.kernel.NewtonIncrements()
        # self.increments.init()
        # self.distro_system = DistributedSystemMSWells(self.kernel)
        self.linear_system = None
        self.newton_parameters = newton_parameters
        self.dxmax = 1.0

    def increment(self):
        # ghosts_synchronizer = Synchronizer(self.distro_system, True)
        # x = self.linear_system.x
        # ghosts_synchronizer.synchronize(x)
        # ghosts_synchronizer.retrieve_solution(self.increments)
        # self.kernel.IncPrimSecdMSWells_PrimToSecd(self.increments)
        # relaxation = self.kernel.Newton_compute_relaxation_mswells(self.increments)
        ## relaxation = 1.0
        # self.dxmax = self.kernel.Newton_get_last_max_inc_mswells()
        # self.kernel.IncCV_NewtonIncrement(self.increments, relaxation)
        if verbose:
            mpi.master_print("   Relaxation factor/dxmax:", relaxation, self.dxmax)
        # self.kernel.NN_flash_all_control_volumes()

    def init_iteration(self):
        self.kernel.IncPrimSecdMSWells_compute()
        self.kernel.LoisThermoHydroMSWells_compute()  # it must be compute again after VSHydroMSWells_init
        self.kernel.VSHydroMSWells_compute()

    def solve(self):

        # Perform Additional Initializations
        self.kernel.mswells_copy_states_from_reservoir()
        #######################################################################
        # For reinit (also modify ResiduMSWells_allocate)
        # also do not call self.kernel.VSHydroMSWells_init() below
        # self.kernel.IncCVMSWells_init_from_file ()
        #######################################################################
        self.kernel.mswells_init_edge_data()
        self.kernel.IncPrimSecdMSWells_compute()
        self.kernel.LoisThermoHydroMSWells_compute()
        self.kernel.mswells_init_leaf_data()
        self.kernel.VSHydroMSWells_init()  # Computes initial pressure, don't call this when reinit from file

        # time iteration
        t = self.newton_parameters.t0
        dt = self.newton_parameters.dt
        total_iters = 0
        while t < self.newton_parameters.tf:

            self.kernel.IncCV_SaveIncPreviousTimeStep()
            self.init_iteration()
            self.kernel.ResiduMSWells_reset_history()
            self.kernel.ResiduMSWells_compute(dt)
            self.kernel.JacobianMSWells_ComputeJacSm(dt, True)  # Perform alignment
            exit()
            ##res0 = 0
            ##resrel = 1.0
            ##mpi.master_print("------- Time iteration at t,dt:", t, dt, "-------")
            ##t = t + dt
            ##for iteration in range(self.newton_parameters.newtonMaxIter):
            ##    ## if nprocs == 1 and total_iters == 297:
            ##    ##            print("MPI rank:", mpi.proc_rank)
            ##    ##            print("Process id:", pid)
            ##    ##            self.kernel.init_wait_for_debug(True)

            ##    pcsp = np.copy(simulation.mswell_node_states().p)
            ##    pcsT = np.copy(simulation.mswell_node_states().T)
            ##    pcsS = np.copy(simulation.mswell_node_states().S)
            ##    total_iters = total_iters + 1

            ##    self.linear_system = PetscLinearSystemMSWell(self.kernel)
            ##    self.linear_system.set_from_jacobian()
            ##    # Save Jacobian, Residual data
            ##    # self.kernel.JacobianMSWells_print_LA_info_to_file(dt,total_iters-1)
            ##    #######################################################################
            ##    # Direct Solver
            ##    linear_solver = PetscDirectSolver(self.linear_system)
            ##    # settings = IterativeSolverSettings(
            ##    #    "gmres",
            ##    #    1e-6,
            ##    #    150,
            ##    #    30,
            ##    # )
            ##    # linear_solver=  PetscIterativeSolver(self.linear_system,settings)
            ##    #######################################################################
            ##    linear_solver.solve()
            ##    # Save Petsc lsystem data
            ##    # self.linear_system.dump_ascii(".")
            ##    # exit()
            ##    if verbose:
            ##        mpi.master_print("   Newton % 3d " % (iteration + 1))
            ##    self.increment()
            ##    res_norm = np.linalg.norm(self.residuals.own_mswell_nodes, 1, axis=0)
            ##    global_res_norm = np.zeros(res_norm.shape[0], dtype=np.double)
            ##    MPI.COMM_WORLD.Allreduce(res_norm, global_res_norm, MPI.SUM)
            ##    if iteration == 0:
            ##        res0 = np.sum(global_res_norm)
            ##    else:
            ##        resrel = np.sum(global_res_norm) / res0
            ##    if verbose:
            ##        mpi.master_print(
            ##            "           residuals",
            ##            global_res_norm,
            ##            "res0",
            ##            res0,
            ##            "rel",
            ##            resrel,
            ##        )

            ##    if mpi.is_on_master_proc:
            ##        f_reswell.write(
            ##            str(total_iters)
            ##            + str("\t")
            ##            + str(np.sum(global_res_norm))
            ##            + str("\t")
            ##            + str(resrel)
            ##            + str("\n")
            ##        )
            ##    self.init_iteration()
            ##    self.kernel.ResiduMSWells_compute(dt)
            ##    self.kernel.JacobianMSWells_ComputeJacSm(dt, True)  # Perform alignment
            ##    mpi.synchronize()
            ##    if (
            ##        resrel < self.newton_parameters.crit_resrel
            ##        or self.dxmax < self.newton_parameters.crit_dxmax
            ##    ):
            ##        break

            ##if (iteration == (self.newton_parameters.newtonMaxIter - 1)) and (
            ##    resrel > self.newton_parameters.crit_resrel
            ##):
            ##    mpi.master_print("--------------------------------------------")
            ##    mpi.master_print("------ Warning:  Newton not converged!!-------")
            ##    t = t - dt
            ##    dt = 0.5 * dt
            ##    self.kernel.IncCV_LoadIncPreviousTimeStep()
            ##    if total_iters > (self.newton_parameters.TotalMaxIter - 1):
            ##        exit()
            ##else:
            ##    mpi.master_print(
            ##        "--- ** Succed, t,dt, nllew_iter, tnew_iter, resrel, :",
            ##        t - dt,
            ##        dt,
            ##        iteration + 1,
            ##        total_iters,
            ##        resrel,
            ##    )
            ##    if verbose:
            ##        allreduce = mpi.communicator().allreduce
            ##        mpi.master_print(
            ##            "max p variation",
            ##            allreduce(
            ##                np.fabs(simulation.mswell_node_states().p - pcsp).max(),
            ##                mpi.MPI.MAX,
            ##            ),
            ##        )
            ##        mpi.master_print(
            ##            "max T variation",
            ##            allreduce(
            ##                np.fabs(simulation.mswell_node_states().T - pcsT).max(),
            ##                mpi.MPI.MAX,
            ##            ),
            ##        )
            ##        mpi.master_print(
            ##            "max S variation",
            ##            allreduce(
            ##                np.fabs(simulation.mswell_node_states().S - pcsS).max(),
            ##                mpi.MPI.MAX,
            ##            ),
            ##        )
            ##    dt = min(1.1 * dt, self.newton_parameters.dtmax)
            ##    self.kernel.IncCVMSWells_print_info_to_file()
            ##    self.kernel.JacobianMSWells_print_IP_info_to_file(t)
            ##if mpi.is_on_master_proc:
            ##    f_reswell.write(str("\n"))


######################################################
## Run Newton ##
mynewton_par = NewtonParameters(dt)
mynewton = MSWellsNewton(mynewton_par)
mynewton.solve()
