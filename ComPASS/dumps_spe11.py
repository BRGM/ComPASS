import numpy as np

from . import mpi
from . import dump_wells as dw
from .dumps import Dumper
from .physics import diphasicCO2


CO2_molar_mass = 44.01e-3
H2O_molar_mass = 18e-3


def mass_density_fraction(molar_density, molar_fraction, molar_mass):
    mass_densities = molar_density[:, None] * molar_fraction * molar_mass
    mass_density = mass_densities.sum(axis=1)
    mass_fraction = mass_densities / mass_density[:, None]
    return mass_density, mass_fraction


def CO2_solubility_ratio(P, T, Xg, Xl):
    return diphasicCO2.cpp_liquid_CO2_fugacity(
        P, T, Xl
    ) / diphasicCO2.cpp_gas_CO2_fugacity(P, T, Xg)


class DumperSPE11(Dumper):
    def dump_states(self, tag="", dump_fluxes=True):
        assert self.simulation_running
        result = {}
        simulation = self.simulation
        states_locations = simulation.states_locations()
        phases = simulation.phases()
        if len(phases) == 1:
            phases = ["fluid"]
        components = simulation.components()
        for location, states in states_locations:
            result[f"{location} context"] = states.context
            result[f"{location} pressure"] = states.p
            result[f"{location} temperature"] = states.T
            result[
                f"{location} total specific enthalpy"
            ] = simulation.total_specific_enthalpy(states)
            if len(phases) > 1:
                for phk, phase in enumerate(phases):
                    result[f"{location} {phase} saturation"] = states.S[:, phk]
            if len(components) > 1:
                for ci, comp in enumerate(components):
                    for phk, phase in enumerate(phases):
                        name = f"{location} {comp} fraction in {phase}"
                        result[name] = states.C[:, phk, ci]
            #
            result[f"{location} CO2 solubility ratio"] = CO2_solubility_ratio(
                states.p, states.T, states.C[:, 0], states.C[:, 1]
            )
            molar_density_fct = [
                simulation.cpp_gas_molar_density,
                simulation.cpp_liquid_molar_density,
            ]
            for phk, phase in enumerate(phases):
                molar_fraction = states.C[:, phk]
                molar_density = molar_density_fct[phk](
                    states.p, states.T, molar_fraction
                )
                mass_density, mass_fraction = mass_density_fraction(
                    molar_density,
                    molar_fraction,
                    [CO2_molar_mass, H2O_molar_mass],
                )
                result[f"{location} {phase} mass density"] = mass_density
                for ci, comp in enumerate(components):
                    result[
                        f"{location} {comp} mass fraction in {phase}"
                    ] = mass_fraction[:, ci]
            #
        if dump_fluxes:
            mass_fluxes_locations = simulation.mass_fluxes_locations()
            for location, fluxes in mass_fluxes_locations:
                result[f"{location} total mass flux"] = fluxes.sum(axis=1)
                if len(components) > 1:
                    for ci, comp in enumerate(components):
                        result[f"{location} {comp} mass flux"] = fluxes[:, ci, :]
            enthalpy_fluxes_locations = simulation.enthalpy_fluxes_locations()
            for location, fluxes in enthalpy_fluxes_locations:
                result[f"{location} enthalpy flux"] = fluxes
            for mf, hf in zip(mass_fluxes_locations, enthalpy_fluxes_locations):
                location, mfluxes = mf
                _, hfluxes = hf
                assert location == _
                h = np.linalg.norm(hfluxes, axis=1) / np.linalg.norm(
                    mfluxes.sum(axis=1), axis=1
                )
                result[f"{location} flowing enthalpy"] = h
        np.savez(self.states_filename(mpi.proc_rank, tag), **result)
        dw.dump_all_wells(simulation, self.to_wells_directory(), tag)
