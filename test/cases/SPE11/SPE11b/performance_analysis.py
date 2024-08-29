# the argument is the output DIRECTORY (no file)
# python performance_analysis.py .../output-.../
# it creates spe11_performance_time_series.csv in the same directory
import sys
import numpy as np
from pathlib import Path
import yaml


YEAR = int(3.1536e7)  # seconds


class ConvergenceLog:
    """
    A data class to store the contents of a yaml timeloop log file.
    It flattens the nested dictionaries loaded by
    yaml.safe_load("timeloop_log.yaml") into several arrays that contain
    the time step (all info), the computation time, the timestep which succeed,
    the Newton iterations (the last element contains the Newton iterations of the
    successful attempt), the linear solver iterations for each Newton iteration
    for every Newton attempts.

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
    """

    def __init__(self, dictionary):
        # dict with key the index/numero of each time step,
        # value the info (dict) of each time step
        self.time_step = {}
        # computation time for each timestep
        self.comp_time = []
        # starting time of each timestep
        self.time = []
        # time step of the successful Newton attempt
        self.success_dt = []
        # number of Newton iterations of all attempts for each timestep
        self.newton_iterations = []
        # number of linear solver iterations of the Newton iterations
        # of all attempts for each timestep
        self.linear_solver = []
        self.get_yaml_data(dictionary)

    def get_yaml_data(self, dictionary):
        for (key, ts) in dictionary.items():
            # key is str 'time step 1'
            # dict with time step index and time_step summary (also a dict)
            self.time_step[int(key[10:])] = ts
            self.comp_time.append(ts["computation_time"])
            self.time.append(ts["time"])
            self.success_dt.append(ts["success_dt"])
            self.newton_iterations.append(
                np.asarray(ts["newton_iterations_per_attempt"])
            )
            linear_solver = []
            for ls in ts["lsolver_iterations_per_newton"]:
                linear_solver.append(np.asarray(ls))
            self.linear_solver.append(linear_solver)
        self.comp_time = np.asarray(self.comp_time)
        self.time = np.asarray(self.time)
        self.success_dt = np.asarray(self.success_dt)


def safe_load(directory):
    dir = Path(directory)
    if not dir.is_dir():
        raise IOError(f"Could not find {str(dir)} directory.")
    yaml_file = "timeloop_log.yaml"
    datapath = dir / yaml_file
    if not datapath.exists():
        raise IOError(f"Could not find {str(datapath)} file.")

    with open(datapath, "r") as f:
        return yaml.safe_load(f)


class SPE11_perf:
    def __init__(self, log):
        self.mass = "n/a"
        self.nres = "n/a"
        self.tlinsol = "n/a"
        # dof is not clear: is it the size of the linear system ?
        # or the number of sites ?
        # cells and secondary unknowns are eliminated
        # we decide that dof = size of the linear system (cst in time)
        # it is misleading in our case because 2D transformed
        # into a 3D mesh so we solve twice the pb
        # todo: get the value of NbNode !
        # with nx = 300, ny = 100
        # NbNode = 301 * 101 * 2
        # with nx = 360, ny = 120
        NbNode = 421 * 141 * 2
        self.dof = 3 * NbNode
        # careful time must begin with 0 at injection starting point
        self.time = log.time - 1000 * YEAR
        # average dt
        self.tstep = log.success_dt
        # nb of time step failure
        self.fsteps = np.array([newt_it.size - 1 for newt_it in log.newton_iterations])
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
        # average value over the last 5 years
        # no average for t = 0 (start of the injection)
        # check close to t = 0 (in sec)
        assert abs(perf.time[0]) < 10
        print(
            format(0, perf.tstep[0], perf.fsteps[0]),
            f"{perf.mass}, ",
            file=f,
            end="",
        )
        print(format(perf.dof, perf.nliter[0]), f"{perf.nres}, ", file=f, end="")
        print(format(perf.liniter[0], perf.runtime[0]), f"{perf.tlinsol}", file=f)

        # loop over all report time (report_t = 0 already done)
        dtb = [*range(5 * YEAR, 1001 * YEAR, 5 * YEAR)]
        first_l = 1
        last_l = 0
        for report_t in dtb:
            while last_l < len(perf.time) and perf.time[last_l] <= report_t:
                last_l += 1
            # print all values (average, sum, ... between first_l and last_l) in file
            print(
                format(
                    report_t,
                    perf.tstep[first_l:last_l].mean(),
                    perf.fsteps[first_l:last_l].mean(),
                ),
                f"{perf.mass}, ",
                file=f,
                end="",
            )
            print(
                format(perf.dof, perf.nliter[first_l:last_l].sum()),
                f"{perf.nres}, ",
                file=f,
                end="",
            )
            print(
                format(
                    perf.liniter[first_l:last_l].sum(),
                    perf.runtime[first_l:last_l].sum(),
                ),
                f"{perf.tlinsol}",
                file=f,
            )

            first_l = last_l
            # end of available data
            if first_l == len(perf.time):
                break

            assert perf.time[first_l] > report_t


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

    data_output_dir = sys.argv[1]
    perf = process_yaml_info(data_output_dir)
    write_csv(
        data_output_dir + "spe11b_performance_time_series.csv",
        perf,
    )
