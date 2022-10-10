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
from mswell_only_newton import MSWellsNewtonParam
from mswell_only_newton import MSWellsNewton

################################################################################################
# A vertical well with 4x4 grid basis and nv horizontal layers over H thickness
ds = 100  # horizontal cell size
H = 1000  # height of the well
nv = 10  #  number of vertical layers
hns = 2  # half the number of cells along basis side
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
mswid = 0  # mswell id - could be any number
################################################################################################
# Important Note:
# It is crucial to be verified that:
# i) mswells/LeafMSWells.F90:LeafMSWells_allocate()
# ii)mswells/LeafMSWells.F90:LeafMSWells_init()
# are configured for the two-mswells and the  type of diphasic model
#################################################################################################


########################################################################
if water2phase:
    simulation = ComPASS.load_physics("water2ph")
else:
    simulation = ComPASS.load_physics("immiscible2ph")

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


def create_well(center):
    return simulation.create_vertical_well(center, rw, multi_segmented=True, zmin=0)


def make_producers():
    wid = mswid
    centers = []
    centers.append((-100, 0))
    centers.append((100, 0))
    mswells = []
    nbwells = 2

    for i in range(nbwells):
        mswells.append(create_well(centers[i]))
        if ip_monitor:
            mswells[i].operate_on_flowrate = QwOut, PwOut  # water2phase, IP-Monitor
        else:
            mswells[i].operate_on_pressure = (
                PwOut,
                100000,
            )  # immiscible, water2phase no monitoring

        mswells[i].produce()
        mswells[i].id = wid
        wid = wid + 1
    return mswells


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
    set_dirichlet_nodes=select_dirichlet_nodes,
    wells=make_producers,
    cell_porosity=omega_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=K_reservoir,
)


# -- Set initial state and boundary conditions
hp = hydrostatic_pressure(0, H, 2 * nv + 1)
initial_state = simulation.build_state(simulation.Context.liquid, p=ptop, T=Ttop)
simulation.all_states().set(initial_state)
dirichlet = simulation.dirichlet_node_states()
dirichlet.set(initial_state)  # will init all variables: context, states...


def set_pT_distribution(states, z):
    states.p[:] = hp(z)
    states.T[:] = Ttop  # simulation.Tsat(states.p[:])# for liquid_context
    # states.T[:] = simulation.Tsat(states.p[:])  # for gas_liquid_context


set_pT_distribution(simulation.node_states(), simulation.vertices()[:, 2])
set_pT_distribution(simulation.cell_states(), simulation.compute_cell_centers()[:, 2])
set_pT_distribution(dirichlet, simulation.vertices()[:, 2])


# simulation.close_well(mswid)#close one well
################################################################################
# Newton and solver parameters
t0 = 0.0
tf = 1000  # 2000.0  # 3000.0  # 3000.0
dtmax = 40.0
if water2phase:
    newtonMaxIter = 40
else:
    newtonMaxIter = 100
TotalMaxIter = 2000
crit_resrel = 1.0e-9
crit_dxmax = 1.0e-10
######################################################
mynewton_par = MSWellsNewtonParam(
    t0=t0,
    tf=tf,
    dt=dt,
    dtmax=dtmax,
    newtonMaxIter=newtonMaxIter,
    TotalMaxIter=TotalMaxIter,
    verbose=verbose,
)
mynewton = MSWellsNewton(mynewton_par, simulation=simulation)
## Run Newton ##
mynewton.solve()
