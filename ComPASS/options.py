import sys
from . import mpi
from inept import Config


class LinalgConfig(Config):
    """
    Set of options for ComPASS linear solvers
    """

    with _.options.on:
        view: bool
        with _.switch.on:
            with _.switch.on as legacy:
                with _.group.on as iterative:
                    with _.group.on as gmres:
                        restart: int = 30
                    tolerance: float = 1e-6
                    maxit: int = 150
                    activate_cpramg: bool = True
                direct: bool
            with _.switch as new:
                with _.group.on as iterative:
                    with _.switch:
                        with _.group.on as gmres:
                            restart: int = 30
                        bcgs: bool
                    tolerance: float = 1e-6
                    maxit: int = 150
                    with _.switch as pc:
                        with _.switch.on as cpramg:
                            hypre: bool = True
                            gamg: bool
                        bjacobi: bool
                        none: bool
                direct: bool


class CallbackConfig(Config):
    """
    Options for triggering callbacks
    """

    with _.options:
        abort_on_linear_failure: bool
        abort_on_newton_failure: bool
        dump_system_on_linear_failure: bool
        abort: float  # Format : --callbacks.abort <time>
        simulation_log: bool  # Format --callbacks.simulation_log True
        linear_system_dump: str  # Format : --callbacks.linear_system_dump t1,t2,t3
        linear_system_binary_dump: str  # Format : --callbacks.linear_system_binary_dump t1,t2,t3


class ComPASSConfig(Config):

    with _.options:
        lsolver = LinalgConfig.root
        callbacks = CallbackConfig.root


compass_config = ComPASSConfig()
compass_config.load_cli()


def get(name, default=None):
    _, *args = sys.argv
    if name not in args:
        return default
    idx = args.index(name) + 1
    assert idx < len(args), "no value provided after option"
    return args[idx]


def get_bool(name, default=False):
    _, *args = sys.argv
    if name not in args:
        return default
    else:
        return True
