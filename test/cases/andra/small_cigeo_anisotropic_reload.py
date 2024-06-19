# -*- coding: utf-8 -*-
#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

# load the diphasic physics, set some initialisations values,
# load the mesh, identify global flags and Dirichlet BC,
# set the rocktypes, Pc, kr functions
from small_cigeo_isotherm_reload import *

ComPASS.set_output_directory_and_logfile(__file__)


# modify the bulk thermal conductivity in W/m/K
# to set anisotropic values
K_reservoir = np.zeros((len(cells_tag), 3, 3), dtype=np.float)
K0 = np.zeros((3, 3), dtype=np.float)
K0[0, 0] = 2.0  # horizontal
K0[1, 1] = 2.0  # horizontal
K0[2, 2] = 1.3  # vertical
K1 = np.zeros((3, 3), dtype=np.float)
K1[0, 0] = 2.1  # horizontal
K1[1, 1] = 2.1  # horizontal
K1[2, 2] = 1.4  # vertical
K_reservoir[circulation_cells_selection] = K0
K_reservoir[bentonite_cells_selection] = K0
K_reservoir[EDZ_cells_selection] = K1
K_reservoir[porous_cells_selection] = K1
K_reservoir[topcover_cells_selection] = K1


if __name__ == "__main__":

    simulation.init(
        mesh=mesh,
        set_dirichlet_nodes=select_dirichlet_nodes,
        cell_porosity=omega,
        cell_permeability=k_reservoir,
        cell_thermal_conductivity=K_reservoir,
        set_global_flags=set_global_flags,
        set_global_rocktype=set_global_rocktype,
    )

    # reload state from initialization, some modifications are needed after
    # Init done with atm nodes on a smaller mesh
    # need a mapping
    snapshot_directory = "output-small_cigeo_anisotropic_atm_gal_init"
    reload_states(snapshot_directory)

    # the known states have been reloaded
    # necessary to change the context of the nodes which were FF
    # the other values of the states are kept
    reset_freeflow_context()

    # necessary to init the states of the nodes and cells
    # inside the circulation gallery and bentonite
    # (were eliminated in the init state and replaced by atmospheric BC)
    init_gallery_states()

    # export_initial_states()

    simulation.alignment = AlignmentMethod.manual
    # simulation.alignment = AlignmentMethod.inverse_diagonal # it is worse
    lsolver = linear_solver(simulation, tolerance=1e-8, direct=False)
    newton = Newton(simulation, 1e-6, 15, lsolver)
    tsmger = TimeStepManager(
        initial_timestep=0.01 * year,
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
        final_time=30 * year,
        output_period=5 * year,
        nitermax=3,
    )

    tsmger.increase_factor = 1.3
    newton.maximum_number_of_iterations = 8
    current_time = run_loop(
        initial_time=current_time,
        final_time=30 * year,
        output_period=5 * year,
    )

    tsmger.maximum = 5.1 * year
    current_time = run_loop(
        initial_time=current_time,
        final_time=5000 * year,
        output_period=1e2 * year,
    )
    tsmger.maximum = 1e4 * year
    current_time = run_loop(
        initial_time=current_time,
        output_period=1e5 * year,
    )

    # simulation.postprocess()
