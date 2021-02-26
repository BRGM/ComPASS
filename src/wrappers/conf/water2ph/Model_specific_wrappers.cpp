#include <pybind11/numpy.h>

#include <array>

#include "../common/enum_to_rank.h"
#include "../common/fortran_thermodynamics.h"
#include "Model_wrappers.h"
#include "StateObjects.h"

constexpr int NC = ComPASS_NUMBER_OF_COMPONENTS;
constexpr int NP = ComPASS_NUMBER_OF_PHASES;
static_assert(NP == 2, "Wrong numpber of phases.");
static_assert(NC == 1, "Wrong numpber of components.");
static_assert(ComPASS_NUMBER_OF_CONTEXTS == 3, "Wrong number of contexts.");
// FIXME: assuming liquid phase is the latest phase
constexpr int GAS_PHASE = 0;  // FIXME: cpp starts at 0 where fortran at 1
constexpr int LIQUID_PHASE = 1;

enum struct Component { water = ComPASS_SINGLE_COMPONENT };

enum struct Phase { gas = ComPASS_GAS_PHASE, liquid = ComPASS_LIQUID_PHASE };

enum struct Context {
   gas = ComPASS_GAS_CONTEXT,
   liquid = ComPASS_LIQUID_CONTEXT,
   diphasic = ComPASS_DIPHASIC_CONTEXT
};

void init_model() {}

void finalize_model() {}

template <int PHASE>
inline double phase_molar_density(double p, double T) {
   double xsi, dxsidp, dxsidT;
   double C[NC] = {1};
   double dxsidC[NC] = {0};
   FluidThermodynamics_molar_density(PHASE + 1, p, T, C, xsi, dxsidp, dxsidT,
                                     dxsidC);
   return xsi;
}

inline double gas_molar_density(double p, double T) {
   return phase_molar_density<GAS_PHASE>(p, T);
}

inline double liquid_molar_density(double p, double T) {
   return phase_molar_density<LIQUID_PHASE>(p, T);
}

template <int PHASE>
inline double phase_molar_enthalpy(double p, double T) {
   double h, dhdp, dhdT;
   double C[NC] = {1};
   double dhdC[NC] = {0};
   FluidThermodynamics_molar_enthalpy(PHASE + 1, p, T, C, h, dhdp, dhdT, dhdC);
   return h;
}

inline double gas_molar_enthalpy(double p, double T) {
   return phase_molar_enthalpy<GAS_PHASE>(p, T);
}

inline double liquid_molar_enthalpy(double p, double T) {
   return phase_molar_enthalpy<LIQUID_PHASE>(p, T);
}

template <int PHASE>
inline double phase_dynamic_viscosity(double p, double T) {
   double mu, dmudp, dmudT;
   double C[NC] = {1};
   double dmudC[NC] = {0};
   FluidThermodynamics_dynamic_viscosity(PHASE + 1, p, T, C, mu, dmudp, dmudT,
                                         dmudC);
   return mu;
}

inline double gas_dynamic_viscosity(double p, double T) {
   return phase_dynamic_viscosity<GAS_PHASE>(p, T);
}

inline double liquid_dynamic_viscosity(double p, double T) {
   return phase_dynamic_viscosity<LIQUID_PHASE>(p, T);
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
   module.def("Psat", py::vectorize(Psat));
   module.def("Tsat", py::vectorize(Tsat));
   module.def("gas_molar_density", py::vectorize(gas_molar_density));
   module.def("gas_molar_enthalpy", py::vectorize(gas_molar_enthalpy));
   module.def("gas_dynamic_viscosity", py::vectorize(gas_dynamic_viscosity));
   module.def("liquid_molar_density", py::vectorize(liquid_molar_density));
   module.def("liquid_molar_enthalpy", py::vectorize(liquid_molar_enthalpy));
   module.def("liquid_dynamic_viscosity",
              py::vectorize(liquid_dynamic_viscosity));

   py::enum_<Component>(module, "Component").value("water", Component::water);

   py::enum_<Context>(module, "Context")
       .value("gas", Context::gas)
       .value("liquid", Context::liquid)
       .value("diphasic", Context::diphasic);

   py::enum_<Phase>(module, "Phase")
       .value("gas", Phase::gas)
       .value("liquid", Phase::liquid);

   use_context_locker(module);

   module.def(
       "build_state",
       [](py::object context, py::object p, py::object T, py::object Sg) {
          auto find_context = [&]() {
             if (!Sg.is_none() || p.is_none() || T.is_none())
                throw std::runtime_error(
                    "You must set a monophasic states providing p and T.");
             double Tsat, foo;  // dummy value
             FluidThermodynamics_Tsat(p.cast<double>(), Tsat, foo);
             if (T.cast<double>() >= Tsat) return Context::gas;
             return Context::liquid;
          };
          auto set_monophasic_state = [&](X &state) {
             if (!Sg.is_none() || p.is_none() || T.is_none())
                throw std::runtime_error(
                    "You must set a monophasic states providing p and T.");
             state.p = p.cast<double>();
             state.T = T.cast<double>();
             state.S.fill(0);
             for (auto &&Ck : state.C) Ck.fill(1);
          };
          auto set_diphasic_state = [&](X &state) {
             bool ok = !(p.is_none() && T.is_none());
             ok &= p.is_none() || T.is_none();
             ok &= !Sg.is_none();
             if (!ok)
                throw std::runtime_error(
                    "You must set diphasic states providing p or T, and Sg.");
             double dummy;
             if (T.is_none()) {
                state.p = p.cast<double>();
                FluidThermodynamics_Tsat(state.p, state.T, dummy);
             } else {
                state.T = T.cast<double>();
                FluidThermodynamics_Psat(state.T, state.p, dummy);
             }
             const double S = Sg.cast<double>();
             state.S[enum_to_rank(Phase::gas)] = S;
             state.S[enum_to_rank(Phase::liquid)] = 1. - S;
             for (auto &&Ck : state.C) Ck.fill(1);
          };
          Context context_value =
              context.is_none() ? find_context() : context.cast<Context>();
          X result;
          result.context = static_cast<int>(context_value);
          switch (context_value) {
             case Context::gas:
                set_monophasic_state(result);
                result.S[enum_to_rank(Phase::gas)] = 1;
                break;
             case Context::liquid:
                set_monophasic_state(result);
                result.S[enum_to_rank(Phase::liquid)] = 1;
                break;
             case Context::diphasic:
                set_diphasic_state(result);
                break;
             default:
                throw std::runtime_error("Requested context does not exist!");
          }
          // As we have a single component concentrations in both phases are
          // always 1
          for (auto &&Cphi : result.C) {
             Cphi.fill(1);
          }
          return result;
       },
       py::arg("context") = py::none{}, py::arg("p") = py::none{},
       py::arg("T") = py::none{}, py::arg("Sg") = py::none{},
       R"doc(
Construct a state given a specific context and physical parameters.

Monophasic states (i.e. liquid and gas) must be set-up using pressure and temperature.
Diphasic state must be set-up using either pressure or temperature and gaz phase saturation.
The saturation conditions will be use to compute missing parameters (pressure or temperature).

Parameters
----------

:param context: context (i.e. liquid, gas or diphasic)
:param p: pressure
:param T: temperature
:param Sg: gaz phase saturation

)doc");
}
