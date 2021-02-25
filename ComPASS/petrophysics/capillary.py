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


# def Beaude():
#     Pc_cst = 2.0e5
#     Sg0 = 1.0 - 1.0e-2
#     Sl0 = 1.0 - Sg0
#     A = -Pc_cst * np.log(Sl0) - (Pc_cst / Sl0) * Sg0

#     def f(Sg):
#         if Sg < Sg0:
#             return -Pc_cst * np.log(1.0 - Sg)
#         return Pc_cst * Sg / Sl0 + A

#     def df(Sg):
#         if Sg < S0:
#             return Pc_cst / (1.0 - Sg)
#         return Pc_cst / Sl0

#     return f, df


def set_liquid_capillary_pressure(simulation, f, df):
    # FIXME: cf. holder definition above
    global holder
    holder = _convert_pc_to_phase_pressure_function(f, df)
    simulation.set_phase_pressure_functions(holder)
