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
import os
from ComPASS._kernel import get_kernel

# fmt: off
# A vertical well in the middle of a regular grid
# with dh * dh squared cells, H total thickness and, nv vertical layers
L = 1000               # half-width of the reservoir
dh = 200               # horizontal cell size
H = 200                # thickness of the reservoir / height of the well
nv = 10                 # number of vertical layers
rw = 0.1               # well radius (m)
ptop = 4 * MPa
Teps = 1               # temperature variation temperature at the top of the reservoir is set to Tsat(ptop) - Teps
gravity = 9.81
omega = 0.15           # reservoir porosity
k = 5e-14              # reservoir horizontal permeability in m^2
K = 2                  # bulk thermal conductivity in W/m/K
Qm = 200 * ton / hour
sleep = False
# fmt: on

simulation = ComPASS.load_eos("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(gravity)

nh = 2 * int(L / dh)
assert nh > 0
grid = ComPASS.Grid(
    shape=(nh, nh, nv),
    extent=(2 * L, 2 * L, H),
    origin=(-L, -L, 0),
)


def make_wells():
    wid = 0  # well id - could be any number
    pwell = simulation.create_vertical_well((0, 0), rw, multi_segmented=False)
    pwell.id = wid
    pwell.operate_on_flowrate = Qm, 1.0 * bar
    pwell.produce()

    ##wid = wid + 1
    ##iwell = simulation.create_vertical_well((0, 0), rw, multi_segmented=False)
    ##iwell.id = wid
    ##iwell.operate_on_flowrate = Qm, 1.0 * bar
    ##iwell.inject(1.0)

    wid = wid + 1
    mswell = simulation.create_vertical_well((0, 0), rw, multi_segmented=True)
    mswell.id = wid
    mswell.operate_on_flowrate = Qm, 1.0 * bar
    mswell.produce()
    # return [pwell,iwell, mswell]
    return [pwell, mswell]
    # return [ pwell, iwell]


if sleep:
    pid = os.getpid()
    print("MPI rank:", mpi.proc_rank)
    print("Process id:", pid)
    print("Sleeping...")
    kernel = get_kernel()
    kernel.init_wait_for_debug(True)


simulation.init(
    mesh=grid,
    wells=make_wells,
    cell_porosity=omega,
    cell_permeability=k,
    cell_thermal_conductivity=K,
    well_model="two_phases",
)
