from importlib import import_module

import numpy as np

# FIXME: is this a pybind11 bug?
# This is used as global variable
# If not use the reference counter of the closure
# returned by _convert_pc_to_phase_pressure_function
# goes mad at the end of the program execution
holder = None


def _convert_pc_to_phase_pressure_function(pc, dpcdS):
    pc = np.vectorize(pc)
    dpcdS = np.vectorize(dpcdS)

    def phase_pressure_function(states, rocktypes, dpadS):
        pa = states.pa
        nb_phases = pa.shape[1]
        assert nb_phases == 2
        p = states.p  # reference pressure
        pa[:, 0] = p
        Sg = states.S[:, 0]
        pa[:, 1] = p - pc(Sg)
        dpadS[:, 1] = dpcdS(Sg)  # pc(Sl) with Sg = 1 - Sl

    return phase_pressure_function


def set_liquid_capillary_pressure(simulation, model):
    # FIXME: cf. holder definition above
    global holder
    if isinstance(model, str):
        module = None
        try:  # locally
            module = import_module(model)
        except ModuleNotFoundError:
            try:  # in ComPASS models
                module = import_module(
                    f".petrophysics.models.{model}", package="ComPASS"
                )
            except ModuleNotFoundError:
                raise RuntimeError(
                    f"Could not find the petrophysics model {model} neither locally nor in ComPASS libraries."
                )
        try:
            f, df = module.Pc, module.dPcdS
        except:
            raise RuntimeError(
                f"Could not extact cappillary pressure functions from model {model} located at: {model.__file__}."
            )
    else:
        try:
            f, df = model
        except:
            raise RuntimeError(
                "Could not extact cappillary pressure functions from model."
            )
    holder = _convert_pc_to_phase_pressure_function(f, df)
    simulation.set_phase_pressure_functions(holder)
