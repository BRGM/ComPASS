from . import options
from .runtime import to_output_directory
import os
from . import mpi


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
    def __init__(self, t):
        super().__init__()
        self.has_been_triggered = False
        self.t = t

    def __call__(self, tick):

        if not self.has_been_triggered:
            if tick.time >= self.t:
                if not self.is_on:
                    self.is_on = True
                elif self.is_on:
                    self.is_on = False
                    self.has_been_triggered = True


class TimestepDirpathMaker:
    def __init__(self, basename, flag):

        self.base_dirname = to_output_directory(basename)
        self.flag = flag
        self.current_dirpath = None
        try:
            os.mkdir(self.base_dirname)
        except FileExistsError:
            pass

    def __call__(self, tick):

        self.flag(tick)
        if self.flag.is_on:
            self.current_dirpath = (
                self.base_dirname + f"/it={tick.iteration}_t={tick.time:.3e}"
            )


class FileDumper:
    def __init__(self, dirpath_maker, dump_function, suffix=None):

        self.dirpath_maker = dirpath_maker
        self.dump_function = dump_function
        self.suffix = suffix or "/"

    def __call__(self, dt):

        if self.dirpath_maker.flag.is_on:
            timestep_dirname = self.dirpath_maker.current_dirpath + f"_dt={dt:.3e}"
            dump_dirname = timestep_dirname + f"{self.suffix}"
            try:
                os.mkdir(timestep_dirname)
            except FileExistsError:
                pass
            try:
                os.mkdir(dump_dirname)
            except FileExistsError:
                pass

            self.dump_function(basename=dump_dirname)


class NewtonIterationFileDumper(FileDumper):
    def __call__(self, dt, it):

        self.suffix = f"/Newton_{it}/"
        super().__call__(dt)


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


class NewtonLogCallback:
    def __init__(self, filename, newton):
        self.filename = to_output_directory(filename)
        self.newton = newton
        with open(self.filename, "w+") as f:
            f.write(
                "Iteration number, Time, Successful timestep, Newton iterations, Linear iterations\n\n"
            )

    def __call__(self, tick):
        def writer(tick):
            with open(self.filename, "a") as f:
                # FIXME : Opening file at every timestep isn't performance-wise ideal
                f.write(
                    f"{tick.iteration},   {tick.time:.13e}, {tick.latest_timestep:.13e}, {len(self.newton.lsolver_iterations)}, {self.newton.lsolver_iterations}\n"
                )

        mpi.on_master_proc(writer)(tick)


def get_callbacks_from_options(newton, tick0):

    timestep_callbacks = []
    newton_iteration_callbacks = []
    newton_failure_callbacks = []
    linear_failure_callbacks = []
    linear_system = newton.lsolver.linear_system

    if options.database["dump_system_on_linear_failure"]:
        dirpath_maker = TimestepDirpathMaker("linear_systems", AlwaysOnFlag())
        timestep_callbacks.append(dirpath_maker)
        sys_dumper = FileDumper(
            dirpath_maker, linear_system.dump_binary, suffix="/Linear failure/"
        )
        log_dumper = FileDumper(
            dirpath_maker, newton.lsolver.write_history, suffix="/Linear failure/"
        )
        linear_failure_callbacks.extend([sys_dumper, log_dumper])
    if options.database["abort_on_linear_failure"]:
        linear_failure_callbacks.append(InterruptTrigger("on linear failure"))
    if options.database["abort_on_newton_failure"]:
        newton_failure_callbacks.append(InterruptTrigger("on Newton failure"))

    t_dump_raw = options.database["dump_ls"]
    if t_dump_raw is not None:
        t_dump_list = t_dump_raw.split(",")
        for t_dump in t_dump_list:
            t_dump = float(t_dump)
            dump_flag = TimestepFlag(t_dump)
            dirpath_maker = TimestepDirpathMaker(
                basename="linear_systems", flag=dump_flag
            )
            dump_trigger = NewtonIterationFileDumper(
                dirpath_maker, linear_system.dump_ascii
            )
            solver_log_trigger = NewtonIterationFileDumper(
                dirpath_maker, newton.lsolver.write_history
            )
            timestep_callbacks.append(dirpath_maker)
            newton_iteration_callbacks.extend([dump_trigger, solver_log_trigger])

    t_dump_b_raw = options.database["dump_ls_binary"]
    if t_dump_b_raw is not None:
        t_dump_b_list = t_dump_b_raw.split(",")
        for t_dump_b in t_dump_b_list:
            t_dump_b = float(t_dump_b)
            dump_flag_b = TimestepFlag(t_dump_b)
            dirpath_maker = TimestepDirpathMaker(
                basename="linear_systems", flag=dump_flag_b
            )
            dump_trigger = NewtonIterationFileDumper(
                dirpath_maker, linear_system.dump_binary
            )
            solver_log_trigger = NewtonIterationFileDumper(
                dirpath_maker, newton.lsolver.write_history
            )
            timestep_callbacks.append(dirpath_maker)
            newton_iteration_callbacks.extend([dump_trigger, solver_log_trigger])

    newton_log_filename = options.database["newton_log"]
    if newton_log_filename is not None:
        newton_log_callback = NewtonLogCallback(newton_log_filename, newton)
        timestep_callbacks.append(newton_log_callback)

    t_kill = options.database["kill"]
    if t_kill is not None:
        t_kill = float(t_kill)
        kill_flag = TimestepFlag(t_kill)
        kill_trigger = TimestepInterruptTrigger(kill_flag)
        timestep_callbacks.append(kill_trigger)

    newton.iteration_callbacks += tuple(newton_iteration_callbacks)
    newton.failure_callbacks += tuple(newton_failure_callbacks)
    newton.lsolver.failure_callbacks += tuple(linear_failure_callbacks)
    return tuple(timestep_callbacks)
