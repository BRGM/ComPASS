import numpy as np
import matplotlib.pyplot as plt
import yaml
import sys


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


class SimulationLog:
    """
    A data class to store the contents of a yaml simulation log file.
    It flattens the nested dictionaries loaded by
    yaml.safe_load("simulation_log.yaml") into three lists that contain
    the subdictionaries for time steps, newton loop attempts
    and newton iterations respectively
    """

    def __init__(self, dictionary):
        self.time_steps = []
        self.newton_attempts = []
        self.newton_iterations = []
        self.get_yaml_data(dictionary)

    def get_yaml_data(self, dictionary):
        for time_steps in dictionary:
            for timestep_iteration in find(time_steps, dictionary):
                self.time_steps.append(timestep_iteration)
        for timestep_dictionary in self.time_steps:
            for attempts in timestep_dictionary:
                if attempts != "time" and attempts != "success_dt":
                    for attempt in find(attempts, timestep_dictionary):
                        self.newton_attempts.append(attempt)
        for attempt_dictionary in self.newton_attempts:
            for newtons in attempt_dictionary:
                if newtons != "dt" and newtons != "status":
                    for newton in find(newtons, attempt_dictionary):
                        self.newton_iterations.append(newton)


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


def plot_simulation_log(simulation_log, png_filename):
    """
    A function which plots the general overview of the simulation behaviour:
    For each newton solve attempt, plots the dt value that is attempted and
    the number of iterations performed with distinctive colors for converged,
    diverged and linear failure attempts
    """
    newton_success = []
    newton_linear_failure = []
    newton_failure = []

    dt = get_simulation_data_as_list("dt", simulation_log.newton_attempts)
    timestep_iteration = get_simulation_data_as_list_count(
        simulation_log.newton_attempts
    )
    get_newton_satuts_iterations_list(
        simulation_log.newton_attempts,
        newton_success,
        newton_linear_failure,
        newton_failure,
    )
    fig, axe_number_newton_iteration = plt.subplots(figsize=(15, 6))
    axe_number_newton_iteration.set_xlabel("Time step attempts")
    axe_number_newton_iteration.set_ylabel("Number of Newton iterations", color="b")
    axe_number_newton_iteration.bar(
        timestep_iteration,
        newton_success,
        width=0.6,
        color="b",
        label="Newton convergence",
    )
    axe_number_newton_iteration.bar(
        timestep_iteration,
        newton_linear_failure,
        width=0.6,
        color="r",
        label="Linear failure",
    )
    axe_number_newton_iteration.bar(
        timestep_iteration,
        newton_failure,
        width=0.6,
        color="tab:orange",
        label="Newton divergence",
    )
    axe_number_newton_iteration.tick_params(axis="y", labelcolor="b")

    axe_timestep = axe_number_newton_iteration.twinx()

    axe_timestep.set_ylabel("Timestep (s)", color="g")
    axe_timestep.set_yscale("log", base=10)
    axe_timestep.plot(timestep_iteration, dt, color="g", linewidth=4)
    axe_timestep.tick_params(axis="y", labelcolor="g")

    fig.tight_layout()
    axe_number_newton_iteration.legend(loc="upper left", ncol=1)
    plt.savefig(png_filename)


if __name__ == "__main__":
    # Type in your yaml simulation log file as the first command line argument
    yaml_file = sys.argv[1]

    with open(f"{yaml_file}", "r") as f:
        data = yaml.safe_load(f)

    simulation_log = SimulationLog(data)
    plot_simulation_log(simulation_log, "simulation_log.png")
