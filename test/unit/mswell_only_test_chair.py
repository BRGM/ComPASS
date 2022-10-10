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
from ComPASS.utils.various import tensor_coordinates
import ComPASS.io.mesh as io
import ComPASS.dump_wells as dw
import ComPASS.mpi as mpi
from mpi4py import MPI
import os
from mswell_only_newton import MSWellsNewtonParam
from mswell_only_newton import MSWellsNewton


################################################################################################
# A vertical well with 2x2 grid basis and nv horizontal layers over H thickness
H = 1000  # height of the well
#################################################
# Coarse
nv = 10  # number of vertical layers
nh = 20  # number of horizontal layers
#################################################
# Fine I
# nv = 100  # number of vertical layers
# nh = 200  # number of x-horizontal layers
################################################
rw = 0.05  # well radius
ptop = 500000  # pressure at the top of the reservoir, 10*MPa one phase. 1*MPa for two phases
Ttop = 350  # , Kelvin, temperature at the top of the reservoir
vgradT = 170 / km  # degrees per km - to reach 300 degC at the bottom
geotherm = lambda zeta: Ttop
gravity = 9.81
injection = False
Tinjection = degC2K(30.0)  # injection temperature - convert Celsius to Kelvin degrees
epsilon = 0.0001
omega_reservoir = 0.15  # reservoir porosity
k_reservoir = (
    1e-16  # reservoir permeability in m^2 - low value to limit reservoir convection
)
K_reservoir = 2  # bulk thermal conductivity in W/m/K
QwOut = 15.0  # Maximum imposed flow rate for producer
PwOut = 5.0e5  # Minmum  imposed pressure for producer
dt = 0.5  # initial timestep
# dt = 0.25  # initial timestep
verbose = True
water2phase = True
ip_monitor = True
sleep = False
nprocs = MPI.COMM_WORLD.Get_size()
######################################################################
# Mesh
Lx, Ly, Lz = 1000.0, 1000.0, 1000.0
Ox, Oy, Oz = -Lx, -Ly, 0
nx, ny, nz = nh, 2, nv

dx = Lx / nx
lhb = 250.0  # length of the horizontal branch
assert abs(lhb % dx) < 1e-5
nhb = int(lhb / dx)  # length of the horizontal branch in cells

################################################################################################
# Important Note:
# It is crucial to be verified that:
# i) mswells/LeafMSWells.F90:LeafMSWells_allocate()
# ii)mswells/LeafMSWells.F90:LeafMSWells_init()
# are configured for the chair mswell and the  type of diphasic model
#################################################################################################

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


slim = H
grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(2 * Lx, 2 * Ly, Lz),
    origin=(Ox, Oy, Oz),
)


def create_well():
    vertices = simulation.global_vertices()
    v = vertices.reshape((nx + 1, ny + 1, nz + 1, 3), order="F")
    indices = np.reshape(
        np.arange(vertices.shape[0]), (nx + 1, ny + 1, nz + 1), order="F"
    )
    # print(v[:,0,0]) # along Ox
    # print(v[1,:,0]) # along Oy starting from v[1,0,0]
    # print(indices[1,:,0]) # indices
    # print(vertices[indices[1,:,0]])
    ic, jc, kc = (nx + 1) // 2, (ny + 1) // 2, (nz + 1) // 2  # center indices
    # print("Vertical branch", "-" * 20)
    vertical_branch = indices[ic, jc, :]
    # print(vertices[vertical_branch])
    # print("Horizontal branch", "-" * 20)
    horizontal_branch = indices[ic : ic + nhb + 1, jc, kc]
    # print(vertices[horizontal_branch])
    # print("Small vertical branch", "-" * 20)
    small_vertical_branch = indices[ic + nhb, jc, : kc + 1]
    # print(vertices[small_vertical_branch])
    # exit()
    vertical_branch.shape = -1, 1
    horizontal_branch.shape = -1, 1
    small_vertical_branch.shape = -1, 1
    # well segments must be oriented from wellhead downwards
    segments = np.vstack(
        [
            # nodes are from bottom to top
            np.hstack([vertical_branch[1:], vertical_branch[:-1]]),
            # order of nodes is ok
            np.hstack([horizontal_branch[:-1], horizontal_branch[1:]]),
            # nodes are from bottom to top
            np.hstack([small_vertical_branch[1:], small_vertical_branch[:-1]]),
        ]
    )
    return simulation.create_well_from_segments(
        segments, well_radius=rw, multi_segmented=True
    )


def make_producer():
    wid = 0  # well id - could be any number
    mswell = create_well()

    if ip_monitor:
        mswell.operate_on_flowrate = QwOut, PwOut  # water2phase, IP-Monitor
    else:
        mswell.operate_on_pressure = (
            PwOut,
            100000,
        )  # immiscible, water2phase no monitoring

    mswell.produce()
    mswell.id = wid
    return [mswell]


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
    wells=make_producer,
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
##########################################################################################
# Newton and solver parameters
t0 = 0.0
tf = 2000  # 2000.0  # 3000.0  # 3000.0
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
