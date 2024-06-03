# spe11_performance_time_series.csv
import sys
import numpy as np

# change it !
sys.path.append("../../../utilities")
from yaml_utilities import ConvergenceLog, safe_load


class SPE11_perf:
    def __init__(self, log):
        self.mass = "n/a"
        self.nres = "n/a"
        self.tlinsol = "n/a"
        # dof is not clear: is it the size of the linear system ?
        # or the number of sites ?
        # cells and secondary unknowns are eliminated
        # we decide that dof = size of the linear system (cst in time)
        # todo: get the value of NbNode !
        NbNode = 10952
        self.dof = 3 * NbNode
        # careful time must begin with 0 at injection starting point
        year: float = 31536.0e3  # seconds
        self.time = log.time - 1000 * year
        # average dt
        self.tstep = log.success_dt
        # nb of time step failure
        self.fsteps = [newt_it.size - 1 for newt_it in log.newton_iterations]
        # total number of non linear iterations (converged and failure) for each timestep
        self.nliter = np.asarray([n.sum() for n in log.newton_iterations])
        # total number of linear iterations (converged and failure) for each timestep
        nb_ls_it_in_newton = []
        for ls_in_newton_attempts in log.linear_solver:
            nb_ls_it_in_newton.append(
                np.sum([ls.sum() for ls in ls_in_newton_attempts])  # 2D array
            )
        self.liniter = np.asarray(nb_ls_it_in_newton)

        # iteration computation time
        self.runtime = log.comp_time


def format(*args):
    output = ""
    for a in args:
        output += f"{a:1.3e}, "
    return output


def write_legend(file):
    print(
        "# t [s], tstep [s], fsteps [-], mass [kg], dof [-], nliter [-], nres [-], liniter [-], runtime [s], tlinsol [s]",
        file=file,
    )


def write_csv(
    filename,
    perf,
):
    with open(filename, "w") as f:
        write_legend(f)
        print(
            format(perf.time[0], perf.tstep[0], perf.fsteps[0]),
            f"{perf.mass}, ",
            file=f,
            end="",
        )
        print(format(perf.dof, perf.nliter[0]), f"{perf.nres}, ", file=f, end="")
        print(format(perf.liniter[0], perf.runtime[0]), f"{perf.tlinsol}", file=f)


def process_yaml_info(yaml_directory):
    print("\n------------------------------------------------------")
    print(f"Processing info from {str(yaml_directory)} directory\n")
    # load the timeloop_log.yaml file from yaml_directory
    conv_data = safe_load(yaml_directory)
    # interpret the yaml file and creates arrays with time, newton
    # and linear solver information
    log = ConvergenceLog(conv_data)
    # transform into SPE11 format
    return SPE11_perf(log)


if __name__ == "__main__":

    perf = process_yaml_info(sys.argv[1])
    write_csv(
        "spe11b_performance_time_series.csv",
        perf,
    )
