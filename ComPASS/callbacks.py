from .options import compass_config
from .runtime import to_output_directory
from .linalg.exceptions import explain_reason
import os
from . import mpi
import time
import yaml


class Flag:
    def __init__(self):
        self.is_on = False

    def __call__(self, *args, **kwargs):
        pass


class AlwaysOnFlag(Flag):
    def __init__(self):
        super().__init__()
        self.is_on = True


class TimestepFlag(Flag):
    """
    A flag that will turn on at the first timestep where tick.time > self.t
    and turn off for the rest of the timeloop
    """

    def __init__(self, t):
        super().__init__()
        self.has_been_triggered = False
        self.t = t

    def __call__(self, tick):

        if not self.has_been_triggered:
            if tick.time > self.t:
                if not self.is_on:
                    self.is_on = True
                    self.t = tick.time
                elif self.is_on:
                    self.is_on = False
                    self.has_been_triggered = True


class NewtonDumper:
    def __init__(self, dump_functions, basedir, flag=None):

        self.dump_functions = dump_functions
        self.basedir = to_output_directory(basedir)
        self.flag = flag or AlwaysOnFlag()

    def __call__(self, newton_tick):

        time_tick = newton_tick.timeloop_tick
        self.flag(time_tick)

        if self.flag.is_on:
            tit = time_tick.iteration
            time = time_tick.time
            dt = newton_tick.current_dt
            nit = newton_tick.iteration
            dump_dirname = (
                f"{self.basedir}/it={tit}_t={time:.4e}/dt={dt:.4e}/newton_it={nit}"
            )

            try:
                os.makedirs(dump_dirname)
            except FileExistsError:
                pass

            for function in self.dump_functions:
                function(basename=dump_dirname)


class InterruptTrigger:
    def __init__(self, message=None):
        self.message = message

    def __call__(self, *args, **kwargs):
        mpi.master_print(f"\n ComPASS - Abortion requested {self.message}\n")
        mpi.abort()


class TimestepInterruptTrigger(InterruptTrigger):
    def __init__(self, timeflag, message=None):
        self.flag = timeflag
        self.message = message or f" at time t > {self.flag.t:.4e} s\n"

    def __call__(self, tick):
        # Update flag
        self.flag(tick)
        if self.flag.is_on:
            super().__call__()


class SimulationLogCallback:
    def __init__(self, filename, newton):
        self.dict = {}
        self.last_point_in_time = time.time()
        self.attempts = [{}]
        self.filename = filename
        self.newton = newton
        with open(to_output_directory("simulation_log.yaml"), "w") as f:
            pass  # Clearing the file if it already exists

    def timeloop_callback(self, tick):
        self.attempts[-1]["status"] = "success"
        timestep_dict = {"time": tick.time - tick.latest_timestep}
        for i, attempt_dict in enumerate(self.attempts):
            timestep_dict.update({f"attempt {i+1}": attempt_dict})
        timestep_dict.update({"success_dt": tick.latest_timestep})
        self.dict.update({f"time_step {tick.iteration}": timestep_dict})
        if mpi.is_on_master_proc:
            with open(to_output_directory("simulation_log.yaml"), "a") as f:
                yaml.safe_dump(
                    {f"time_step {tick.iteration}": timestep_dict},
                    f,
                    default_flow_style=False,
                    sort_keys=False,
                )
        self.attempts = [{}]

    def newton_iteration_callback(self, newton_tick):
        self.attempts[-1]["dt"] = newton_tick.current_dt
        self.attempts[-1][f"newton {newton_tick.iteration}"] = {
            "linear_iterations": self.newton.lsolver.nit,
            "linear_solver_status": explain_reason(self.newton.lsolver.ksp_reason),
            "cpu_time": time.time() - self.last_point_in_time,
            "residuals": f"{self.newton.convergence_scheme.pv_norms()[:]} {self.newton.relative_residuals[-1]}",
        }
        self.last_point_in_time = time.time()

    def linear_failure_callback(self, newton_tick):
        self.attempts[-1]["dt"] = newton_tick.current_dt
        self.attempts[-1][f"newton {newton_tick.iteration}"] = {
            "linear_iterations": self.newton.lsolver.nit,
            "linear_solver_status": explain_reason(self.newton.lsolver.ksp_reason),
            "cpu_time": time.time() - self.last_point_in_time,
        }
        self.attempts[-1]["status"] = "linear_failure"
        self.last_point_in_time = time.time()
        self.attempts.append({})

    def newton_failure_callback(self, newton_tick):
        self.newton_iteration_callback(newton_tick)
        self.attempts[-1]["status"] = "newton_failure"
        self.attempts[-1]["dt"] = newton_tick.current_dt
        self.attempts.append({})


def get_callbacks_from_options(
    newton, tick0,
):

    timestep_callbacks = []
    newton_iteration_callbacks = []
    newton_failure_callbacks = []
    linear_failure_callbacks = []
    linear_system = newton.lsolver.linear_system

    if compass_config.get("callbacks.dump_system_on_linear_failure"):
        dump_trigger = NewtonDumper(
            (linear_system.dump_binary, newton.lsolver.write_history),
            "linear_systems/",
        )
        linear_failure_callbacks.extend([dump_trigger])
    if compass_config.get("callbacks.abort_on_linear_failure"):
        linear_failure_callbacks.append(InterruptTrigger("on linear failure"))
    if compass_config.get("callbacks.abort_on_newton_failure"):
        newton_failure_callbacks.append(InterruptTrigger("on Newton failure"))

    if compass_config.get("callbacks.linear_system_dump"):
        t_dump_raw = compass_config["callbacks.linear_system_dump"]
        t_dump_list = t_dump_raw.split(",")
        for t_dump in t_dump_list:
            dump_trigger = NewtonDumper(
                (linear_system.dump_ascii, newton.lsolver.write_history),
                "linear_systems/",
                flag=TimestepFlag(float(t_dump)),
            )
            newton_iteration_callbacks.extend([dump_trigger])

    if compass_config.get("callbacks.linear_system_binary_dump"):
        t_dump_raw = compass_config["callbacks.linear_system_binary_dump"]
        t_dump_list = t_dump_raw.split(",")
        for t_dump in t_dump_list:
            dump_trigger = NewtonDumper(
                (linear_system.dump_binary, newton.lsolver.write_history),
                "linear_systems/",
                flag=TimestepFlag(float(t_dump)),
            )
            newton_iteration_callbacks.extend([dump_trigger])

    if compass_config.get("callbacks.simulation_log"):
        simulation_log_filename = compass_config["callbacks.simulation_log"]
        simulation_log_callback = SimulationLogCallback(simulation_log_filename, newton)
        timestep_callbacks.append(simulation_log_callback.timeloop_callback)
        newton.iteration_callbacks += (
            simulation_log_callback.newton_iteration_callback,
        )
        newton.lsolver.failure_callbacks += (
            simulation_log_callback.linear_failure_callback,
        )
        newton.failure_callbacks += (simulation_log_callback.newton_failure_callback,)

    if compass_config.get("callbacks.abort") is not None:
        t_kill = compass_config["callbacks.abort"]
        t_kill = float(t_kill)
        kill_flag = TimestepFlag(t_kill)
        kill_trigger = TimestepInterruptTrigger(kill_flag)
        timestep_callbacks.append(kill_trigger)

    newton.iteration_callbacks += tuple(newton_iteration_callbacks)
    newton.failure_callbacks += tuple(newton_failure_callbacks)
    newton.lsolver.failure_callbacks += tuple(linear_failure_callbacks)

    return tuple(timestep_callbacks)
