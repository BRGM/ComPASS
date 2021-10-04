#include <pybind11/numpy.h>
// Used to expose Phase_vector
// FIXME: Use a specific binding fro Phase_vector (which are just contiguous
// buffers)
#include <pybind11/stl.h>

#include <array>

#include "DefModel.h"
#include "Equilibriums.h"
#include "LoisThermoHydro.h"
#include "Model_wrappers.h"
#include "StateObjects.h"
#include "Thermodynamics.h"
#include "enum_to_rank.h"

constexpr int NC = ComPASS_NUMBER_OF_COMPONENTS;
constexpr int NP = ComPASS_NUMBER_OF_PHASES;
static_assert(NP == X::Model::np, "Wrong number of phases.");
static_assert(NC == X::Model::nc, "Wrong number of components.");
static_assert(NP == 2, "Wrong number of phases.");
static_assert(NC == 2, "Wrong number of components.");
static_assert(ComPASS_NUMBER_OF_CONTEXTS == 3, "Wrong number of contexts.");

using Phase_vector = X::Model::Phase_vector;
using Real = X::Model::Real;

// FIXME: assuming liquid phase is the latest phase
constexpr int GAS_PHASE = 0;
constexpr int LIQUID_PHASE = 1;

void init_model() {}

void finalize_model() {}

template <int PHASE>
inline double phase_molar_density(double p, double T,
                                  py::array_t<double, py::array::c_style> &C) {
   double xsi, dxsidp, dxsidT;
   double dxsidC[NC] = {0};
   FluidThermodynamics_molar_density(PHASE + 1, p, T, C.data(), xsi, dxsidp,
                                     dxsidT, dxsidC);
   return xsi;
}

inline double liquid_molar_density(double p, double T,
                                   py::array_t<double, py::array::c_style> &C) {
   return phase_molar_density<LIQUID_PHASE>(p, T, C);
}

inline double gas_molar_density(double p, double T,
                                py::array_t<double, py::array::c_style> &C) {
   return phase_molar_density<GAS_PHASE>(p, T, C);
}

template <int PHASE>
inline double phase_molar_enthalpy(double p, double T,
                                   py::array_t<double, py::array::c_style> &C) {
   double h, dhdp, dhdT;
   double dhdC[NC] = {0};
   FluidThermodynamics_molar_enthalpy(PHASE + 1, p, T, C.data(), h, dhdp, dhdT,
                                      dhdC);
   return h;
}

inline double liquid_molar_enthalpy(
    double p, double T, py::array_t<double, py::array::c_style> &C) {
   return phase_molar_enthalpy<LIQUID_PHASE>(p, T, C);
}

inline double gas_molar_enthalpy(double p, double T,
                                 py::array_t<double, py::array::c_style> &C) {
   return phase_molar_enthalpy<GAS_PHASE>(p, T, C);
}

template <int PHASE>
inline double phase_dynamic_viscosity(
    double p, double T, py::array_t<double, py::array::c_style> &C) {
   double mu, dmudp, dmudT;
   double dmudC[NC] = {0};
   FluidThermodynamics_dynamic_viscosity(PHASE + 1, p, T, C.data(), mu, dmudp,
                                         dmudT, dmudC);
   return mu;
}

inline double gas_dynamic_viscosity(
    double p, double T, py::array_t<double, py::array::c_style> &C) {
   return phase_dynamic_viscosity<GAS_PHASE>(p, T, C);
}

inline double liquid_dynamic_viscosity(
    double p, double T, py::array_t<double, py::array::c_style> &C) {
   return phase_dynamic_viscosity<LIQUID_PHASE>(p, T, C);
}

inline double Psat(double T) {
   double result;
   double dPsatdT;
   FluidThermodynamics_Psat(T, result, dPsatdT);
   return result;
}

inline double Tsat(double p) {
   double result;
   double dTsatdp;
   FluidThermodynamics_Tsat(p, result, dTsatdp);
   return result;
}

inline double specific_mass(const Phase &phase, const X &x) {
   double rho, drhodp, drhodT;
   double drhodC[NC];
   FluidThermodynamics_specific_mass(enum_to_rank(phase) + 1, x.p, x.T,
                                     x.C[enum_to_rank(phase)].data(), rho,
                                     drhodp, drhodT, drhodC);
   return rho;
}

void add_specific_model_wrappers(py::module &module) {
   module.def("liquid_molar_enthalpy", py::vectorize(liquid_molar_enthalpy));
   module.def("gas_molar_enthalpy", py::vectorize(gas_molar_enthalpy));
   module.def("liquid_molar_density", py::vectorize(liquid_molar_density));
   module.def("gas_molar_density", py::vectorize(gas_molar_density));
   module.def("liquid_dynamic_viscosity",
              py::vectorize(liquid_dynamic_viscosity));
   module.def("gas_dynamic_viscosity", py::vectorize(gas_dynamic_viscosity));
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
          py::object Cag, py::object Cal) {
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
             update_phase_pressures(state);
             auto [Cga, Cla] = diphasic_equilibrium(state.pa, state.T);
             state.C[gas][air] = Cga;
             state.C[gas][water] = 1 - Cga;
             state.C[liquid][air] = Cla;
             state.C[liquid][water] = 1 - Cla;
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

)doc");

   module.def("specific_mass", &specific_mass);

   module.def(
       "diphasic_equilibrium",
       [](py::tuple pa, const double T, const double atol,
          const std::size_t maxiter) {
          return diphasic_equilibrium({pa[0].cast<Real>(), pa[1].cast<Real>()},
                                      T, atol, maxiter);
       },
       py::arg("pa"), py::arg("T"), py::arg("atol") = 1.e-8,
       py::arg("maxiter") = 1000);
}
