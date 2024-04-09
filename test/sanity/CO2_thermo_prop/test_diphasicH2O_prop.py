# %%
import numpy as np
from ComPASS.utils.dbgtools import check_derivatives, check_partial_derivatives
from ComPASS.utils.units import *
import matplotlib.pyplot as plt
from thermo import *

# Load EoS
import ComPASS.physics.diphasicCO2 as s


def calculate_density(T, P):
    eos = SRK(Tc=647.096, Pc=220.640 * bar, omega=0.3443, T=T, P=P * bar)
    return eos.rho_l  # , eos.drho_dP_l, eos.drho_dT_l


def calculate_enthalpy(T, P):
    eos = SRK(Tc=647.096, Pc=220.640 * bar, omega=0.3443, T=T, P=P * bar)
    return eos.H_dep_l  # , eos.dH_dep_dP_l, eos.dH_dep_dT_l


def plot_results():
    pp = np.linspace(100, 150, 5)
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 10))

    T = 300  # Example of temperature

    for ax, property_name in zip(axes.flat, ["Density", "Enthalpy"]):
        ax.set_ylabel(property_name)
        ax.set_xlabel("Pressure (Bar)")

    density_compass = np.vectorize(
        lambda T, x: s.cpp_liquid_molar_density(x * bar, T, np.array([0, 1]))
    )
    print(density_compass(T, pp))
    density_thermo = np.vectorize(lambda T, x: calculate_density(T, x))

    ax = axes[0]
    ax.plot(pp, density_compass(T, pp), "-b", label="ComPASS")
    ax.plot(pp, density_thermo(T, pp), ":r", label="Thermo")
    ax.legend()

    # Enthalpy comparison
    energy_compass = np.vectorize(
        lambda T, x: s.cpp_liquid_molar_enthalpy(x * bar, T, np.array([0, 1])) * 0.001
    )
    # energy_thermo = np.vectorize(lambda T, x: calculate_enthalpy(T, x))
    print(energy_compass(T, pp))

    ax = axes[1]
    ax.plot(pp, energy_compass(T, pp), "-b", label="ComPASS")
    # ax.plot(pp, energy_thermo(T, pp), ":r", label="Thermo")
    ax.legend()

    plt.tight_layout()
    plt.show()


plot_results()
# %%


def calculate_fugacity_coeff(T, P):
    eos = SRK(Tc=304.13, Pc=73.77 * bar, omega=0.22394, T=T, P=P * bar)
    return eos.fugacity_g / bar


def calculate_density(T, P):
    eos = SRK(Tc=304.13, Pc=73.77 * bar, omega=0.22394, T=T, P=P * bar)
    return eos.rho_g


def calculate_enthalpy(T, P):
    eos = SRK(Tc=304.13, Pc=73.77 * bar, omega=0.22394, T=T, P=P * bar)
    return eos.H_dep_g


def plot_results():
    pp = np.linspace(10, 50, 5)
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 10))

    T = 300  # Example of temperature

    for ax, property_name in zip(
        axes.flat, ["Fugacity of gas", "Density pf gas", "Enthalpy of gas"]
    ):
        ax.set_ylabel(property_name)
        ax.set_xlabel("Pressure (Bar)")

    # Fugacity comparison
    fg_compass = np.vectorize(
        lambda T, x: s.cpp_gas_CO2_fugacity(x * bar, T, np.array([1, 0]))
    )  # fugacity in bar?
    fg_thermo = np.vectorize(lambda T, x: calculate_fugacity_coeff(T, x))

    ax = axes[0, 0]
    ax.plot(pp, fg_compass(T, pp), "-b", label="ComPASS")
    ax.plot(pp, fg_thermo(T, pp), ":r", label="Thermo")
    ax.legend()

    # Density comparison
    density_compass = np.vectorize(
        lambda T, x: s.cpp_gas_molar_density(x * bar, T, np.array([1, 0]))
    )
    density_thermo = np.vectorize(lambda T, x: calculate_density(T, x))

    ax = axes[0, 1]
    ax.plot(pp, density_compass(T, pp), "-b", label="ComPASS")
    ax.plot(pp, density_thermo(T, pp), ":r", label="Thermo")
    ax.legend()

    # Enthalpy comparison
    energy_compass = np.vectorize(
        lambda T, x: s.cpp_gas_molar_enthalpy(x * bar, T, np.array([1, 0]))
    )
    energy_thermo = np.vectorize(lambda T, x: calculate_enthalpy(T, x))

    ax = axes[1, 0]
    ax.plot(pp, energy_compass(T, pp), "-b", label="ComPASS")
    ax.plot(pp, energy_thermo(T, pp), ":r", label="Thermo")
    ax.legend()

    # Enthalpy comparison
    # enthalpy_compass = np.vectorize(lambda T, x: calculate_enthalpy(T, x)[0])
    # enthalpy_thermo = np.vectorize(lambda T, x: calculate_enthalpy(T, x)[1])
    #
    # ax = axes[1, 1]
    # ax.plot(pp, enthalpy_compass(T, pp), "-b", label="ComPASS")
    # ax.plot(pp, enthalpy_thermo(T, pp), ":r", label="Thermo")
    # ax.legend()

    plt.tight_layout()
    plt.show()


plot_results()

# %%
