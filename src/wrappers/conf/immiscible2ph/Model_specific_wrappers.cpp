#include <pybind11/numpy.h>

#include "DefModel.h"
#include "LoisThermoHydro.h"
#include "Model_wrappers.h"
#include "StateObjects.h"
#include "Thermodynamics.h"
#include "enum_to_rank.h"

constexpr int NC = ComPASS_NUMBER_OF_COMPONENTS;
constexpr int NP = ComPASS_NUMBER_OF_PHASES;
static_assert(NP == 2, "Wrong numpber of phases.");
static_assert(NC == 2, "Wrong numpber of components.");
static_assert(ComPASS_NUMBER_OF_CONTEXTS == 3, "Wrong number of contexts.");

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
   // PHASE + 1 : convert Python to Fortran convention
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
   // PHASE + 1 : convert Python to Fortran convention
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
   // PHASE + 1 : convert Python to Fortran convention
   return FluidThermodynamics_dynamic_viscosity(PHASE + 1, p, T, C.data());
}

inline double cpp_gas_dynamic_viscosity(
    double p, double T, py::array_t<double, py::array::c_style> &C) {
   return phase_dynamic_viscosity<GAS_PHASE>(p, T, C);
}

inline double cpp_liquid_dynamic_viscosity(
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

void add_specific_model_wrappers(py::module &module) {
   module.def("liquid_molar_enthalpy", py::vectorize(liquid_molar_enthalpy));
   module.def("gas_molar_enthalpy", py::vectorize(gas_molar_enthalpy));
   module.def("liquid_molar_density", py::vectorize(liquid_molar_density));
   module.def("gas_molar_density", py::vectorize(gas_molar_density));
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
       [](py::object context, py::object p, py::object T, py::object Sg) {
          constexpr auto gas = enum_to_rank(Phase::gas);
          constexpr auto liquid = enum_to_rank(Phase::liquid);
          constexpr auto air = enum_to_rank(Component::air);
          constexpr auto water = enum_to_rank(Component::water);

          auto set_gas_state = [&](X &state) {
             if (p.is_none())
                throw std::runtime_error(
                    "You dont need to provide pressure for gas state.");
             if (T.is_none())
                throw std::runtime_error(
                    "You dont need to provide temperature for gas state.");
             if (!Sg.is_none())
                throw std::runtime_error(
                    "You dont need to provide saturation for gas state.");
             state.p = p.cast<double>();
             state.T = T.cast<double>();
             state.S.fill(0);
             state.S[gas] = 1;
          };

          auto set_liquid_state = [&](X &state) {
             if (p.is_none())
                throw std::runtime_error(
                    "You dont need to provide pressure for liquid state.");
             if (T.is_none())
                throw std::runtime_error(
                    "You dont need to provide temperature for liquid state.");
             if (!Sg.is_none())
                throw std::runtime_error(
                    "You dont need to provide saturation for liquid state.");
             state.p = p.cast<double>();
             state.T = T.cast<double>();
             state.S.fill(0);
             state.S[liquid] = 1;
          };

          auto set_diphasic_state = [&](X &state) {
             if (p.is_none())
                throw std::runtime_error(
                    "You dont need to provide pressure for diphasic state.");
             if (T.is_none())
                throw std::runtime_error(
                    "You dont need to provide temperature for diphasic state.");
             if (Sg.is_none())
                throw std::runtime_error(
                    "You need to provide gas saturation for diphasic states.");
             state.p = p.cast<double>();
             state.T = T.cast<double>();
             const double S = Sg.cast<double>();
             state.S[gas] = S;
             state.S[liquid] = 1. - S;
          };

          auto set_molar_fractions = [&](X &state) {
             state.C[gas].fill(0);
             state.C[liquid].fill(0);
             state.C[gas][air] = 1;
             state.C[liquid][water] = 1;
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
          set_molar_fractions(result);
          update_phase_pressures(result);
          return result;
       },
       py::arg("context").none(false), py::arg("p") = py::none{},
       py::arg("T") = py::none{}, py::arg("Sg") = py::none{},
       R"doc(
Construct a state given a specific context and physical parameters.

Parameters
----------

:param context: context (i.e. liquid, gas or diphasic)
:param p: pressure
:param T: temperature
:param Sg: gaz phase saturation

)doc");
}
