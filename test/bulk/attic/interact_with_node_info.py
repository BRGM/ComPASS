#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

import ComPASS
import doublet_utils
from ComPASS.utils.units import *

ComPASS.load_eos("water2ph")

omega_reservoir = 0.15  # reservoir porosity
k_reservoir = 1e-12  # reservoir permeability in m^2
K_reservoir = 2  # bulk thermal conductivity in W/m/K

Lx, Ly, Lz = 3000.0, 2000.0, 100.0
Ox, Oy, Oz = -1500.0, -1000.0, -1600.0
nx, ny, nz = 30, 20, 10

ComPASS.set_gravity(0)

grid = ComPASS.Grid(shape=(nx, ny, nz), extent=(Lx, Ly, Lz), origin=(Ox, Oy, Oz),)

ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    mesh=grid,
    cell_porosity=omega_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=K_reservoir,
)


def array_description(a):
    return {val: np.count_nonzero(a == val) for val in np.unique(a)}


def output_node_info(info):
    # This a way to access field as attibutes
    info = np.rec.array(info)
    print("Nb of nodes:", info.proc.shape[0])
    print(array_description(info.proc.view("c")))
    print(array_description(info.frac.view("c")))
    print(array_description(info.pressure.view("c")))
    print(array_description(info.temperature.view("c")))


@ComPASS.mpi.on_master_proc
def output_global_node_info():
    print("Proc info (own/ghost) is not relevant at global scale.")
    info = output_node_info(ComPASS.global_node_info())


def output_local_node_info():
    print("Node info for proc:", ComPASS.mpi.proc_rank)
    info = output_node_info(ComPASS.node_info())


# The following line will crash as the mesh is distributed at the end of ComPASS.init
# consequently the global node info array is deallocated
# output_global_node_info()

output_local_node_info()
