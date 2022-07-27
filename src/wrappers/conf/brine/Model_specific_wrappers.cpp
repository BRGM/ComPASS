#include <pybind11/numpy.h>

#include <array>
#include <string>

#include "DefModel.h"
#include "Model_wrappers.h"
#include "StateObjects.h"
#include "Thermodynamics.h"
#include "enum_to_rank.h"

constexpr int NC = ComPASS_NUMBER_OF_COMPONENTS;
constexpr int NP = ComPASS_NUMBER_OF_PHASES;
static_assert(NP == 1, "Wrong numpber of phases.");
static_assert(NC == 2, "Wrong numpber of components.");
static_assert(ComPASS_NUMBER_OF_CONTEXTS == 1, "Wrong number of contexts.");
// FIXME: assuming liquid phase is the latest phase
constexpr int LIQUID_PHASE = 0;

using fluid_property = decltype(FluidThermodynamics_molar_density);

void init_model() {
   py::print();
   py::print();
   py::print();
   py::print("                   ---------- WARNING ----------");
   py::print();
   py::print(
       "This physics is EXPERIMENTAL and correlations have not been checked!");
   py::print();
   py::print("              DO NOT use it for production purposes!");
   py::print();
   py::print("                   ---------- WARNING ----------");
   py::print();
   py::print();
   py::print();
}

void finalize_model() {}

inline X build_state(double p, double T, double Cs) {
   X result;
   result.context = static_cast<X::Model::Context>(Context::single_context);
   result.p = p;
   result.T = T;
   result.S.fill(1);
   if (Cs < 0 || Cs > 1)
      throw std::runtime_error("Invalid salt molar fraction.");
   result.C[0][enum_to_rank(Component::salt)] = Cs;
   result.C[0][enum_to_rank(Component::water)] = 1 - Cs;
   return result;
}

inline double compute_property(fluid_property f, double p, double T,
                               double Cs) {
   double prop, dpropdp, dpropdT;
   auto state = build_state(p, T, Cs);
   X::Model::Phase_component_matrix dpropdC;
   X::Model::Phase_vector dpropdS;
   f(ComPASS_SINGLE_PHASE, state.p, state.T, state.C.data()->data(), prop,
     dpropdp, dpropdT, dpropdC.data()->data());
   return prop;
}

inline double molar_density(double p, double T, double Cs) {
   return compute_property(FluidThermodynamics_molar_density, p, T, Cs);
}

inline double molar_enthalpy(double p, double T, double Cs) {
   return compute_property(FluidThermodynamics_molar_enthalpy, p, T, Cs);
}

inline double cpp_dynamic_viscosity(double p, double T, double Cs) {
   auto state = build_state(p, T, Cs);
   return FluidThermodynamics_dynamic_viscosity(
       ComPASS_SINGLE_PHASE, state.p, state.T, state.C.data()->data());
}

void add_specific_model_wrappers(py::module &module) {
   auto help_add_parameters = [](const std::string &s) {
      return (s + R"doc(
        Parameters
        ----------

        :param p: pressure
        :param T: temperature
        :param Csalt: salt molar fraction

        )doc")
          .c_str();
   };

   module.def("molar_density", py::vectorize(molar_density),
              help_add_parameters(R"doc(
        Brine molar density

        )doc"));
   module.def("cpp_dynamic_viscosity", py::vectorize(cpp_dynamic_viscosity),
              help_add_parameters(R"doc(
        Brine dynamic visosity

        )doc"));
   module.def("molar_enthalpy", py::vectorize(molar_enthalpy),
              help_add_parameters(R"doc(
        Brine molar enthalpy

        )doc"));

   py::enum_<Component>(module, "Component")
       .value("water", Component::water)
       .value("salt", Component::salt);

   py::enum_<Context>(module, "Context")
       .value("single_context", Context::single_context);

   py::enum_<Phase>(module, "Phase").value("single_phase", Phase::single_phase);

   module.def(
       "build_state",
       [](py::object p, py::object T, py::object Csalt) {
          const double Cs = Csalt.is_none() ? 0 : Csalt.cast<double>();
          return build_state(p.cast<double>(), T.cast<double>(), Cs);
       },
       py::arg("p") = py::none{}, py::arg("T") = py::none{},
       py::arg("Csalt") = py::none{}, help_add_parameters(R"doc(
Construct a state given a specific context and physical parameters.

Monophasic states (i.e. liquid and gas) must be set-up using pressure and temperature.

)doc"));
}
