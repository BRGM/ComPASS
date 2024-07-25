import sys
import os
from pathlib import Path

import numpy as np

import click
import os
import matplotlib.pyplot as plt
import yaml
from yaml_utilities import *

minute = 60  # s
hour = 60 * minute
day = 24 * hour
year = 365.25 * day


# def plot_advanced_simulation_log(yaml_directory, plot_directory):
#     """
#     A function which plots an overview of the simulation behaviour:
#     for each newton algorithm attempt, plots the dt value that is attempted and
#     the number of iterations performed with distinct colors for converged,
#     diverged and linear failure attempts
#     """
#     # load the timeloop_log.yaml file from yaml_directory
#     advanced_conv_data = safe_load_advanced(yaml_directory)

#     # interpret the yaml files and creates arrays with time, newton
#     # and linear solver information

#     newton_success = []
#     newton_linear_failure = []
#     newton_failure = []

#     dt = get_simulation_data_as_list("dt", simulation_log.newton_attempts)
#     timestep_iteration = get_simulation_data_as_list_count(
#         simulation_log.newton_attempts
#     )
#     get_newton_satuts_iterations_list(
#         simulation_log.newton_attempts,
#         newton_success,
#         newton_linear_failure,
#         newton_failure,
#     )
#     fig, axe_number_newton_iteration = plt.subplots(figsize=(15, 6))
#     axe_number_newton_iteration.set_xlabel("Time step attempts")
#     axe_number_newton_iteration.set_ylabel("Number of Newton iterations", color="b")
#     axe_number_newton_iteration.bar(
#         timestep_iteration,
#         newton_success,
#         width=0.6,
#         color="b",
#         label="Newton convergence",
#     )
#     axe_number_newton_iteration.bar(
#         timestep_iteration,
#         newton_linear_failure,
#         width=0.6,
#         color="r",
#         label="Linear failure",
#     )
#     axe_number_newton_iteration.bar(
#         timestep_iteration,
#         newton_failure,
#         width=0.6,
#         color="tab:orange",
#         label="Newton divergence",
#     )
#     axe_number_newton_iteration.tick_params(axis="y", labelcolor="b")

#     axe_timestep = axe_number_newton_iteration.twinx()

#     axe_timestep.set_ylabel("Timestep (s)", color="g")
#     axe_timestep.set_yscale("log", base=10)
#     axe_timestep.plot(timestep_iteration, dt, color="g", linewidth=4)
#     axe_timestep.tick_params(axis="y", labelcolor="g")

#     fig.tight_layout()
#     axe_number_newton_iteration.legend(loc="upper left", ncol=1)
#     figfile = plot_directory / "advanced_log_plot"
#     plt.savefig(figfile)
#     print(f"Figure {figfile} has been created")


def plot_timeloop_analysis(yaml_directory, plot_directory):
    plot_directory = Path(plot_directory)
    if not plot_directory.is_dir():
        raise IOError(f"Could not find {str(plot_directory)} directory.")

    # load the timeloop_log.yaml file from yaml_directory
    conv_data = safe_load(yaml_directory)

    # interpret the yaml file and creates arrays with time, newton
    # and linear solver information
    conv_log = ConvergenceLog(conv_data)

    # number of Newton algo attempts for the whole timeloop
    total_nb_newton_attempts = np.sum([len(n) for n in conv_log.newton_iterations])
    # total number of Newton iterations for each timestep
    nb_newton_it = np.asarray([n.sum() for n in conv_log.newton_iterations])
    # number of Newton iterations in convergent Newton algo for each timestep
    nb_it_in_success_newton = np.asarray([n[-1] for n in conv_log.newton_iterations])

    # number of linear solver iterations in convergent Newton algo for each timestep
    nb_ls_it_in_success_newton = np.asarray(
        [
            ls[-1].sum() for ls in conv_log.linear_solver
        ]  # last Newton attempt has converged
    )
    # mean of linear solver iterations in convergent Newton algo for each timestep
    mean_ls_it_in_success_newton = np.asarray(
        [
            ls[-1].mean() for ls in conv_log.linear_solver
        ]  # last Newton attempt has converged
    )
    # number of linear solver iterations in all Newton iterations for each timestep
    nb_ls_it_in_newton = []
    for ls_in_newton_attempts in conv_log.linear_solver:
        nb_ls_it_in_newton.append(
            np.sum([ls.sum() for ls in ls_in_newton_attempts])  # 2D array
        )
    nb_ls_it_in_newton = np.asarray(nb_ls_it_in_newton)

    n_time_steps = len(conv_log.time_step.keys())
    # sum over all timesteps
    nb_failed_newton = total_nb_newton_attempts - n_time_steps
    total_nb_newton_it = nb_newton_it.sum()
    total_nb_it_in_success_newton = nb_it_in_success_newton.sum()
    total_nb_ls_it_in_newton = nb_ls_it_in_newton.sum()
    total_nb_ls_it_in_success_newton = nb_ls_it_in_success_newton.sum()

    print(f"Number of time step iterations: {n_time_steps}")

    duration = conv_log.success_dt.sum()
    min_dt = conv_log.success_dt.min()
    max_dt = conv_log.success_dt.max()

    # the computation time for each timestep
    timeloop_comp_time = conv_log.comp_time.sum()

    print(f"Physical duration: {duration:6g} s = {duration/year:6g} y")
    # the timeloop computation time is the sum of the time of each time step
    timeloop_comp_time = np.sum(conv_log.comp_time)
    print(
        f"Timeloop computation time: {timeloop_comp_time:6g} s",
        f"= {timeloop_comp_time/minute:6g} min",
        f"= {timeloop_comp_time/hour:6g} h",
    )

    print(
        f"Timestep interval: {min_dt:6g} s = {min_dt/year:6g} y"
        f" < {max_dt:6g} s = {max_dt/year:6g} y"
    )
    print(
        f"Newton algorithms: {nb_failed_newton} Newton attempts failed out of {total_nb_newton_attempts} "
        f"({100*(1.0 - nb_failed_newton/total_nb_newton_attempts):.0f}% succeeded)."
    )
    if total_nb_newton_it > 0:
        print(
            f"    It means {total_nb_newton_it - total_nb_it_in_success_newton} useless Newton iterations out of {total_nb_newton_it} "
            f"({100*(total_nb_it_in_success_newton/total_nb_newton_it):.0f}% useful).\n"
            f"    There is an average of {nb_it_in_success_newton.mean():.0f} Newton iterations (when the Newton algorithm converged)."
        )
        if total_nb_ls_it_in_newton > 0:
            print(
                f"Linear solver: {total_nb_ls_it_in_newton - total_nb_ls_it_in_success_newton} useless iterations out of {total_nb_ls_it_in_newton} "
                f"({100*(total_nb_ls_it_in_success_newton/total_nb_ls_it_in_newton):.0f}% useful).\n"
            )

    plt.clf()
    plt.plot(range(n_time_steps), conv_log.success_dt / year)
    plt.xlabel("iteration")
    plt.ylabel("timestep (y)")
    plt.yscale("log")
    figfile = plot_directory / "timesteps-iterations"
    plt.savefig(figfile)
    print(f"Figure {figfile} has been created")

    plt.clf()
    plt.plot(conv_log.time / year, conv_log.success_dt / year)
    plt.xlabel("time (y)")
    plt.ylabel("timestep (y)")
    plt.yscale("log")
    figfile = plot_directory / "timesteps-time"
    plt.savefig(figfile)
    print(f"Figure {figfile} has been created")

    plt.clf()
    plt.plot(conv_log.time / year, 100 * (1 - nb_it_in_success_newton / nb_newton_it))
    plt.xlabel("time (y)")
    plt.ylabel("percentage of useless Newton iterations")
    plt.ylim(0, 100)
    figfile = plot_directory / "Newton-useless-iterations"
    plt.savefig(figfile)
    print(f"Figure {figfile} has been created")

    plt.clf()
    fig, ax1 = plt.subplots()
    color = "tab:red"
    ax1.set_xlabel("time (y)")
    ax1.plot(conv_log.time / year, nb_it_in_success_newton, color=color)
    ax1.set_ylabel("number of Newton iterations at success attempts", color=color)
    ax1.tick_params(axis="y", labelcolor=color)
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = "tab:blue"
    ax2.plot(conv_log.time / year, mean_ls_it_in_success_newton, color=color)
    ax2.set_ylabel("mean of Linear solver iterations at success Newton", color=color)
    ax2.tick_params(axis="y", labelcolor=color)
    figfile = plot_directory / "Newton-linear-solver-cv"
    plt.savefig(figfile)
    print(f"Figure {figfile} has been created")

    plt.clf()
    plt.plot(conv_log.time / year, nb_ls_it_in_success_newton)
    plt.xlabel("time (y)")
    plt.ylabel("linear solver iterations in convergent Newton algo")
    figfile = plot_directory / "LinearSolver-iterations-cv-newton"
    plt.savefig(figfile)
    print(f"Figure {figfile} has been created")

    if (nb_ls_it_in_newton > 0).all():
        plt.clf()
        plt.plot(
            conv_log.time / year,
            100 * (1 - nb_ls_it_in_success_newton / nb_ls_it_in_newton),
        )
        plt.xlabel("time (y)")
        plt.ylabel("percentage of useless Linear solver iterations")
        plt.ylim(0, 100)
        figfile = plot_directory / "LinearSolver-useless-iterations"
        plt.savefig(figfile)
        print(f"Figure {figfile} has been created")

    # Plotting the computation time for each timestep
    plt.clf()
    plt.plot(conv_log.time / year, conv_log.comp_time)
    plt.xlabel("time (y)")
    plt.ylabel("computation time per timestep (s)")
    figfile = plot_directory / "computation-time-per-timestep"
    plt.savefig(figfile)
    print(f"Figure {figfile} has been created")

    return conv_log


def plot_timestep(directories, directory_conv_log):
    plt.clf()
    for i, directory in enumerate(directories):
        path = "/mnt/d/devel/ComPASS/test/cases/andra/comparaison"
        if not os.path.isdir(path):
            os.makedirs(path)
        os.chdir(path)
        plt.plot(
            directory_conv_log[i].time / year,
            directory_conv_log[i].success_dt / year,
            label="{}".format(os.path.basename(os.path.normpath(directory))),
        )
    plt.xlabel("time (y)")
    plt.ylabel("timestep (y)")
    plt.yscale("log")
    plt.xscale("log")
    plt.legend()
    plt.legend(fontsize=8)


def advanced_timeloop_analysis(conv_log, timestep_directory, plot_directory):
    raise "advanced_timeloop_analysis not implemented yet"
    for it in conv_log.time_step.keys():
        yaml_file = timestep_directory / f"time_step_{it}_log.yaml"
        with open(f"{yaml_file}", "r") as file:
            it_data = yaml.safe_load(file)
            breakpoint()


# @click.command(name="timeloop-analysis")
# @click.option(
#     "-p",
#     "--procs",
#     "collect_procs_id",
#     is_flag=True,
#     default=False,
#     help="ouput paraview/mesh.pvtu file with cell distribution (proc variable)",
# )
# @click.argument("directories", nargs=-1)
def timeloop_analysis_command(
    directories,
    advanced=False,
):
    analysis_all_dir = []
    for directory in directories:

        print("\n------------------------------------------------------")
        print(f"Processing info from {str(directory)} directory.\n")
        analysis_all_dir.append(plot_timeloop_analysis(directory, directory))

        if advanced:
            timestep_directory = Path(directory + "/time_step_log")
            if not timestep_directory.is_dir():
                raise IOError(f"Could not find {str(timestep_directory)} directory.")

            print("\n------------------------------------------------------")
            print(f"Processing info from {str(timestep_directory)} directory.\n")
            advanced_timeloop_analysis(
                analysis_all_dir[-1], timestep_directory, directory
            )


if __name__ == "__main__":
    """
    Analyse the information stored in the yaml file timeloop_log.yaml
    The yaml file timeloop_log.yaml is created with the option callbacks.timeloop_log,
    it contains a dictionary with all the time steps,
    for each time step it contains convergence informations.
    Give as argument the path to the output directory of the script you want to analyse.

    Example of timeloop_log.yaml with one time step:
    {
        'time_step 1': {
            'time': 0.0,
            'newton_iterations_per_attempt': [10, 9],
            'lsolver_iterations_per_newton': [[8, 28, 34, 92, 64, 92, 92, 92, 92, 92],
                [8, 28, 34, 61, 90, 90, 86, 69, 72]],
            'success_dt': 75.0,
            'computing_time': 8.746352195739746
        }
    }

    Remark :
    More detailed informations in files time_step_log/time_step_{tick.iteration}_log.yaml,
    are analysed in this file using the "advanced" option.
    """
    timeloop_analysis_command(sys.argv[1:], advanced=False)
