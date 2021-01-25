import sys
from pathlib import Path

import numpy as np

import click
import matplotlib.pyplot as plt


def successfull_iterations(data):
    nb_failed = data[:, 3]
    res = 1 + np.nonzero(nb_failed[:-1] == nb_failed[1:])[0]
    assert res[-1] + 1 == data.shape[0]
    # is the first iteration ok ?
    if nb_failed[0] == 0:
        return np.hstack([[0], res])
    return res


minute = 60
hour = 60 * minute
year = 365.25 * 24 * hour


def timeloop_analysis(directory):
    directory = Path(directory)
    if not directory.is_dir():
        raise IOError(f"Could not find {str(directory)} directory.")
    datapath = directory / "timeloop"
    if not datapath.exists():
        raise IOError(f"Could not find {str(datapath)} file.")
    data = np.loadtxt(datapath)

    print(f"Number of iterations: {int(data[-1, 0]-data[0, 0]) + 1}")
    dt = data[:, 2]
    duration = data[-1, 1] - data[0, 1] + dt[0]  # add first timestep
    print(f"Physical duration: {duration:6g} s = {duration/year:6g} y")
    print(
        f"Total simulation time: {data[-1, 5]:6g} s",
        f"= {data[-1, 5]/minute:6g} min",
        f"= {data[-1, 5]/hour:6g} h",
    )

    print(
        f"Timestep: {dt.min():6g} s = {dt.min()/year:6g} y"
        f" < {dt.max():6g} s = {dt.max()/year:6g} y"
    )
    nb_good = int(data[-1, 3])
    nb_failed = int(data[-1, 4])
    print(
        f"Newon iterations: {nb_failed} failed out of {nb_good + nb_failed} "
        f"({nb_good} = {100*(nb_good/(nb_good + nb_failed)):.0f}% succeeded)"
    )

    plt.clf()
    plt.plot(data[:, 0], data[:, 2] / year)
    plt.xlabel("iteration")
    plt.ylabel("timestep (y)")
    plt.yscale("log")
    plt.savefig(directory / "timesteps-iterations")

    plt.clf()
    plt.plot(data[:, 1] / year, data[:, 2] / year)
    plt.xlabel("time (y)")
    plt.ylabel("timestep (y)")
    plt.yscale("log")
    plt.savefig(directory / "timesteps-time")

    plt.clf()
    plt.plot(data[:, 1] / year, 100 * (data[:, 4] / (data[:, 3] + data[:, 4])))
    plt.xlabel("time (y)")
    plt.ylabel("percentage of useless Newton iterations")
    plt.ylim(0, 100)
    plt.savefig(directory / "useless-iterations")


@click.command(name="timeloop-analysis")
# @click.option(
#     "-p",
#     "--procs",
#     "collect_procs_id",
#     is_flag=True,
#     default=False,
#     help="ouput paraview/mesh.pvtu file with cell distribution (proc variable)",
# )
@click.argument("directories", nargs=-1)
def timeloop_analysis_command(directories,):
    for directory in directories:
        timeloop_analysis(directory)


if __name__ == "__main__":
    timeloop_analysis_command(sys.argv[1:])
