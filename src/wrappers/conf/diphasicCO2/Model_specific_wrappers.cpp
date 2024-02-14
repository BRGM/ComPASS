#include <pybind11/numpy.h>
// Used to expose Phase_vector
// FIXME: Use a specific binding fro Phase_vector
// (which are just contiguous buffers)
#include <pybind11/stl.h>

#include <array>

#include "DefModel.h"
#include "Equilibriums.h"
#include "LoisThermoHydro.h"
#include "Model_wrappers.h"
#include "StateObjects.h"
#include "Thermodynamics.h"
#include "enum_to_rank.h"

static_assert(X::Model::np == 2, "Wrong number of phases.");
static_assert(X::Model::nc == 2, "Wrong number of components.");

// FIXME: assuming liquid phase is the latest phase
// FIXME: cpp starts at 0 where fortran at 1
static_assert(enum_to_rank(Phase::gas) == 0, "Wrong gas phase rank.");
static_assert(to_underlying(Phase::gas) == 1, "Wrong gas phase id.");
static_assert(enum_to_rank(Phase::liquid) == 1, "Wrong liquid phase rank.");
static_assert(to_underlying(Phase::liquid) == 2, "Wrong liquid phase id.");

static_assert(std::is_same_v<X::Model::Real, double>, "Wrong real type");
using Component_vector = X::Model::Component_vector;
using Phase_vector = X::Model::Phase_vector;

extern "C" {
void DiphasicFlash_enforce_consistent_molar_fractions(X &);
}

void init_model() {}

void finalize_model() {}

// nothing to add, all the physical properties are defined in fortran
void add_physical_properties_wrappers(py::module &module) {}

inline double phase_molar_density(const Phase &phase, double p, double T,
                                  py::array_t<double, py::array::c_style> &C) {
   return FluidThermodynamics_molar_density(to_underlying(phase), p, T,
                                            C.data());
}

inline double cpp_liquid_molar_density(
    double p, double T, py::array_t<double, py::array::c_style> &C) {
   return phase_molar_density(Phase::liquid, p, T, C);
}

inline double cpp_gas_molar_density(
    double p, double T, py::array_t<double, py::array::c_style> &C) {
   return phase_molar_density(Phase::gas, p, T, C);
}

inline double phase_molar_enthalpy(const Phase &phase, double p, double T,
                                   py::array_t<double, py::array::c_style> &C) {
   return FluidThermodynamics_molar_enthalpy(to_underlying(phase), p, T,
                                             C.data());
}

inline double cpp_liquid_molar_enthalpy(
    double p, double T, py::array_t<double, py::array::c_style> &C) {
   return phase_molar_enthalpy(Phase::liquid, p, T, C);
}

inline double cpp_gas_molar_enthalpy(
    double p, double T, py::array_t<double, py::array::c_style> &C) {
   return phase_molar_enthalpy(Phase::gas, p, T, C);
}

inline double phase_dynamic_viscosity(
    const Phase &phase, double p, double T,
    py::array_t<double, py::array::c_style> &C) {
   return FluidThermodynamics_dynamic_viscosity(to_underlying(phase), p, T,
                                                C.data());
}

inline double cpp_gas_dynamic_viscosity(
    double p, double T, py::array_t<double, py::array::c_style> &C) {
   return phase_dynamic_viscosity(Phase::gas, p, T, C);
}

inline double cpp_liquid_dynamic_viscosity(
    double p, double T, py::array_t<double, py::array::c_style> &C) {
   return phase_dynamic_viscosity(Phase::liquid, p, T, C);
}

inline double cpp_volumetric_mass_density(
    const Phase &phase, double p, double T,
    py::array_t<double, py::array::c_style> &C) {
   return FluidThermodynamics_volumetric_mass_density(to_underlying(phase), p,
                                                      T, C.data());
}

inline double Psat(double T) {
   double result, dPsatdT;
   FluidThermodynamics_Psat(T, result, dPsatdT);
   return result;
}

inline double Tsat(double p) {
   double result, dTsatdp;
   FluidThermodynamics_Tsat(p, result, dTsatdp);
   return result;
}

void add_specific_model_wrappers(py::module &module) {
   module.def("cpp_liquid_molar_enthalpy",
              py::vectorize(cpp_liquid_molar_enthalpy));
   module.def("cpp_gas_molar_enthalpy", py::vectorize(cpp_gas_molar_enthalpy));
   module.def("cpp_liquid_molar_density",
              py::vectorize(cpp_liquid_molar_density));
   module.def("cpp_gas_molar_density", py::vectorize(cpp_gas_molar_density));
   module.def("cpp_liquid_dynamic_viscosity",
              py::vectorize(cpp_liquid_dynamic_viscosity));
   module.def("cpp_gas_dynamic_viscosity",
              py::vectorize(cpp_gas_dynamic_viscosity));
   module.def("Psat", py::vectorize(Psat));
   module.def("Tsat", py::vectorize(Tsat));

   py::enum_<Component>(module, "Component")
       .value("air", Component::air)
       .value("water", Component::water);

   py::enum_<Context>(module, "Context")
       .value("gas", Context::gas)
       .value("liquid", Context::liquid)
       .value("diphasic", Context::diphasic);

   py::enum_<Phase>(module, "Phase")
       .value("gas", Phase::gas)
       .value("liquid", Phase::liquid);

   module.def(
       "build_state",
       [](py::object context, py::object p, py::object T, py::object Sg,
          py::object Cag, py::object Cal, py::object rocktype) {
          constexpr auto gas = enum_to_rank(Phase::gas);
          constexpr auto liquid = enum_to_rank(Phase::liquid);
          constexpr auto air = enum_to_rank(Component::air);
          constexpr auto water = enum_to_rank(Component::water);

          auto set_gas_state = [&](X &state) {
             if (!Sg.is_none())
                throw std::runtime_error(
                    "You dont need to provide saturation for gas state.");
             if (!Cal.is_none())
                throw std::runtime_error(
                    "You dont need to provide liquid molar fractions for gas "
                    "state.");
             state.p = p.cast<double>();
             state.T = T.cast<double>();
             state.S.fill(0);
             state.S[gas] = 1;
             state.C[gas][air] = Cag.is_none() ? 1. : Cag.cast<double>();
             state.C[liquid][air] = 0;
             DiphasicFlash_enforce_consistent_molar_fractions(state);
          };

          auto set_liquid_state = [&](X &state) {
             if (!Sg.is_none())
                throw std::runtime_error(
                    "You dont need to provide saturation for liquid state.");
             if (!Cag.is_none())
                throw std::runtime_error(
                    "You dont need to provide gas molar fractions for liquid "
                    "state.");
             state.p = p.cast<double>();
             state.T = T.cast<double>();
             state.S.fill(0);
             state.S[liquid] = 1;
             state.C[gas][air] = 1;
             state.C[liquid][air] = Cal.is_none() ? 0 : Cal.cast<double>();
             DiphasicFlash_enforce_consistent_molar_fractions(state);
          };

          auto set_diphasic_state = [&](X &state) {
             if (!(Cag.is_none() && Cal.is_none()))
                throw std::runtime_error(
                    "Don't prescribe molar fractions for diphasic contexts.");
             state.p = p.cast<double>();
             state.T = T.cast<double>();
             const double S = Sg.cast<double>();
             state.S[gas] = S;
             state.S[liquid] = 1. - S;
             if (rocktype.is_none()) {
                update_phase_pressures(state);
             } else {
                update_phase_pressures(state, rocktype.cast<int>());
             }
             auto [Cga, Cla] =
                 diphasic_equilibrium<Component, Phase>(state.pa, state.T);
             state.C[gas][air] = Cga;
             state.C[gas][water] = 1. - Cga;
             state.C[liquid][air] = Cla;
             state.C[liquid][water] = 1. - Cla;
             DiphasicFlash_enforce_consistent_molar_fractions(state);
          };

          Context context_value = context.cast<Context>();
          X result;
          result.context = static_cast<int>(context_value);
          switch (context_value) {
             case Context::gas:
                set_gas_state(result);
                break;
             case Context::liquid:
                set_liquid_state(result);
                break;
             case Context::diphasic:
                set_diphasic_state(result);
                break;
             default:
                throw std::runtime_error("Requested context does not exist!");
          }
          return result;
       },
       py::arg("context").none(false), py::arg("p") = py::none{},
       py::arg("T") = py::none{}, py::arg("Sg") = py::none{},
       py::arg("Cag") = py::none{}, py::arg("Cal") = py::none{},
       py::arg("rocktype") = py::none{},
       R"doc(
Construct a state given a specific context and physical parameters.

Parameters
----------

:param context: context (i.e. liquid, gas or diphasic)
:param p: pressure
:param T: temperature
:param Sg: gaz phase saturation
:param Cag: gaz phase air molar fraction
:param Cal: liquid phase air molar fraction
:param rocktype: rocktype of the site

)doc");
}
