import sys
import os
from pathlib import Path

import numpy as np

import click
import matplotlib.pyplot as plt
import yaml
from ComPASS.utils.units import *
import statistics as st


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

    Remark :
    More detailed informations in files time_step_log/time_step_{tick.iteration}_log.yaml,
    which are loaded in the TimeStepLog class.
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
            # breakpoint()
        self.comp_time = np.asarray(self.comp_time)
        self.time = np.asarray(self.time)
        self.success_dt = np.asarray(self.success_dt)


class TimeStepLog:
    """
    A data class to store the contents of a yaml
    time_step_log/time_step_{tick.iteration}_log file.
    It flattens the nested dictionaries loaded by
    yaml.safe_load() into several arrays that contain ???



    Example of time_step_log/time_step_{tick.iteration}_log.yaml
    when the first newton attempt converged:
    {
        'attempt 1': {
            'dt': 3600,
            'newton 1': {
                'linear_iterations': 2,
                'linear_solver_status': 'CONVERGED_RTOL',
                'cpu_time': 0.08084297180175781,
                'residual': 0.0034475365594200794
            },
            'newton 2': {
                'linear_iterations': 2,
                'linear_solver_status': 'CONVERGED_RTOL',
                'cpu_time': 0.04762387275695801,
                'residual': 0.0003103529511358692
            },
            'newton 3': {
                'linear_iterations': 2,
                'linear_solver_status': 'CONVERGED_RTOL',
                'cpu_time': 0.04584383964538574,
                'residual': 1.8224308015252549e-07
            },
            'status': 'success'
        }
    }
    """

    def __init__(self, dictionary):
        self.get_yaml_data(dictionary)

    def get_yaml_data(self, dictionary):
        for (key, ts) in dictionary.items():
            self.time_step.append(int(key[10:]))
            self.time_step_summary.append(ts)
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


def find(key, dictionary):
    for k, v in dictionary.items():
        if k == key:
            yield v
        elif isinstance(v, dict):
            for result in find(key, v):
                yield result
        elif isinstance(v, list):
            for d in v:
                for result in find(key, d):
                    yield result


def get_simulation_data_as_list_count(dict):
    list = []
    for value in range(len(dict)):
        list.append(value + 1)
    return list


def get_simulation_data_as_list(name, dict):
    list = []
    for value in dict:
        for data in value:
            if data == name:
                list.append(value[data])
    return list


def get_newton_satuts_iterations_list(
    dict, newton_success, newton_linear_failure, newton_failure
):
    for value in range(len(dict)):
        newton_success.append(0)
        newton_linear_failure.append(0)
        newton_failure.append(0)
    count_id = 0
    for value in dict:
        count_newton = 0
        for data in value:
            if data != "status" and data != "dt":
                count_newton += 1
        for data in value:
            if data == "status":
                if value[data] == "success":
                    newton_success[count_id] = count_newton
                elif value[data] == "newton_failure":
                    newton_failure[count_id] = count_newton
                else:
                    newton_linear_failure[count_id] = count_newton
        count_id += 1


def yaml_file(simulation_file_py):
    assert simulation_file_py.is_file(), f"{simulation_file_py} does not exist"
    output_path = simulation_file_py.parent / f"output-{simulation_file_py.stem}"
    assert output_path.is_dir(), f"{output_path} not found"
    return output_path


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


def safe_load_advanced(directory):
    dir = Path(directory)
    if not dir.is_dir():
        raise IOError(f"Could not find {str(dir)} directory.")
    # yaml_file = "timeloop_log.yaml"
    # datapath = dir / yaml_file
    # if not datapath.exists():
    #     raise IOError(f"Could not find {str(datapath)} file.")

    # with open(datapath, "r") as f:
    #     return yaml.safe_load(f)


def analysed_timeloops_log(output_path):
    log = safe_load(output_path)

    analysed_log = {}
    for it_dict in log.values():
        for it_key, it_value in it_dict.items():
            analysed_log.setdefault(it_key, [])
            analysed_log[it_key].append(it_value)
            if it_key == "newton_iterations_per_attempt":
                analysed_log.setdefault("sum_newton_it", [])
                analysed_log["sum_newton_it"].append(sum(it_value))
                analysed_log.setdefault("nb_it_in_success_newton", [])
                analysed_log["nb_it_in_success_newton"].append(it_value[-1])

            if it_key == "lsolver_iterations_per_newton":
                analysed_log.setdefault("mean_lsolver_iterations_per_newton", [])
                analysed_log["mean_lsolver_iterations_per_newton"].append(
                    st.mean(it_value[0])
                )
    return analysed_log


def plot_keys(dict, x, y, output_dir):
    assert x in dict, "First argument nor found in dictionary"
    assert y in dict, "Second argument nor found in dictionary"
    if x == "newton_iterations_per_attempt" or y == "newton_iterations_per_attempt":
        print("""ERROR: USE OF 'newton_iterations_per_attempt' NOT ALLOWED""")
        return
    plt.clf()
    plt.plot(dict[x], dict[y])
    plt.xlabel(x)
    plt.ylabel(y)
    figfile = output_dir / f"{str(y)}-functions-of-{str(x)}"
    plt.savefig(figfile)
    print(f"Figure {figfile} has been created")


def plot_read_data(directories, x, y):
    for directory in directories:
        plot_keys(
            analysed_timeloops_log(yaml_file(directory)),
            x,
            y,
            directory,  # directory where to save figures
        )
