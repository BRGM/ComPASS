import sys
from . import mpi


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


class Database(dict):
    """
    A dictionary structure which holds the command line options
    """

    def __init__(self):

        # Failure options
        self["abort_on_linear_failure"] = get_bool("--abort_on_linear_failure", False)
        self["dump_system_on_linear_failure"] = get_bool(
            "--dump_system_on_linear_failure", False
        )
        self["abort_on_newton_failure"] = get_bool("--abort_on_newton_failure", False)
        # Linear solver parameters
        self["linear_solver_version"] = get("--linear_solver_version", None)
        self["direct_linear_solver"] = get_bool("--direct_linear_solver", False)
        self["cpr_amg_type"] = get("--cpr_amg_type", None)
        self["disable_cpramg"] = get_bool("--disable_cpramg", False)
        self["linear_solver_view"] = get_bool("--linear_solver_view", False)
        # Dump options
        self["dump_ls"] = get("--dump_ls", None)
        self["dump_ls_binary"] = get("--dump_ls_binary", None)
        # Miscellaneous
        self["newton_log"] = get("--newton_log", None)
        self["kill"] = get("--kill", None)


database = Database()
