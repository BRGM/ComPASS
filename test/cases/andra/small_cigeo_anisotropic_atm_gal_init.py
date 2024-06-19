# -*- coding: utf-8 -*-
#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

# load the diphasic physics, set some initialisationsvalues,
# load the mesh, identify global flags and Dirichlet BC,
# set the rocktypes, Pc, kr functions
from small_cigeo_isotherm_atm_gal_init import *

ComPASS.set_output_directory_and_logfile(__file__)


# modify the bulk thermal conductivity in W/m/K
# to set anisotropic values
K_reservoir = np.zeros((len(cells_tag), 3, 3), dtype=float)
# K0 = np.zeros((3, 3), dtype=float)
# K0[0, 0] = 2.0  # horizontal
# K0[1, 1] = 2.0  # horizontal
# K0[2, 2] = 1.3  # vertical
K1 = np.zeros((3, 3), dtype=float)
K1[0, 0] = 2.1  # horizontal
K1[1, 1] = 2.1  # horizontal
K1[2, 2] = 1.4  # vertical
K_reservoir[:] = K1

simulation.init(
    mesh=mesh,
    set_dirichlet_nodes=select_dirichlet_nodes,
    cell_porosity=omega,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=K_reservoir,  # anisotropic value
    set_global_flags=set_global_flags,
    set_global_rocktype=set_global_rocktype,
)

node_flags = simulation.nodeflags()
# identify the atmospheric faces
gallery_faces = select_gallery_faces(node_flags)

# init
init_IC_and_BC(gallery_faces)

# export_initial_states()

simulation.alignment = AlignmentMethod.manual
lsolver = linear_solver(simulation, tolerance=1e-8, direct=False)
newton = Newton(simulation, 1e-6, 50, lsolver)
tsmger = TimeStepManager(
    initial_timestep=4500,
    increase_factor=1.2,
    decrease_factor=0.5,
)

run_loop = lambda initial_time=0, final_time=1e6 * year, output_period=1e4 * year, output_every=None, nitermax=None: simulation.standard_loop(
    initial_time=initial_time,
    final_time=final_time,
    newton=newton,
    time_step_manager=tsmger,
    output_period=output_period,
    output_every=output_every,
    nitermax=nitermax,
)

current_time = run_loop(
    final_time=year,
    output_every=1,
    nitermax=4,
)

tsmger.increase_factor = 1.4
newton.maximum_number_of_iterations = 8
current_time = run_loop(
    initial_time=current_time,
    final_time=year,
    output_every=2,
)

# simulation.postprocess()
