import numpy as np


def _convert_to_phase_pressure_function(f, df):
    f = np.vetorize(f)
    df = np.vetorize(df)

    def phase_pressure_function(states, rocktypes, pa, dpadS):
        nb_phases = pa.shape[1]
        assert np_phases == 2
        p = states.p  # reference pressure
        pa[:, 0] = p
        Sg = states.S[:, 0]
        pa[:, 1] = p + f(Sg)
        dpadS[:, 1] = -df(Sg)  # Sg = 1 - Sl

    return phase_pressure_function


def Beaude():
    Pc_cst = 2.0e5
    Sg0 = 1.0 - 1.0e-2
    Sl0 = 1.0 - Sg0
    A = -Pc_cst * np.log(1.0 - Sg0) - Pc_cst / (1.0 - Sg0) * Sg0

    def f(Sg):
        if Sg < Sg0:
            return -Pc_cst * log(1.0 - Sg)
        return Pc_cst * Sg / Sl0 + A

    def df(Sg):
        if Sg < S0:
            return Pc_cst / (1.0 - Sg)
        return Pc_cst / Sl0

    return f, df


def set_liquid_capillary_pressure(simulation, f, df):
    simulation.set_phase_pressure_functions(_convert_to_phase_pressure_function(f, df))
