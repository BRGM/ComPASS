#include <pybind11/numpy.h>

#include "DefModel.h"
#include "Model_wrappers.h"
#include "StateObjects.h"
#include "Thermodynamics.h"

void init_model() {}

void finalize_model() {}

template <typename Function>
inline double call_physical_function(Function function, double p, double T) {
   constexpr int NC = ComPASS_NUMBER_OF_COMPONENTS;
   constexpr int NP = ComPASS_NUMBER_OF_PHASES;
   static_assert(NC == 1, "we assume there is only one component");
   static_assert(NP == 1, "we assume there is only one phase");
   double C[NC] = {1};
   return function(1, p, T, C);
}

inline double cpp_molar_density(double p, double T) {
   return call_physical_function(FluidThermodynamics_molar_density, p, T);
}

inline double cpp_molar_enthalpy(double p, double T) {
   return call_physical_function(FluidThermodynamics_molar_enthalpy, p, T);
}

inline double cpp_dynamic_viscosity(double p, double T) {
   return call_physical_function(FluidThermodynamics_dynamic_viscosity, p, T);
}

void add_specific_model_wrappers(py::module &module) {
   // provided mainly for consistency... we are filling arrays with constant
   // scalars...
   module.def("cpp_molar_density", py::vectorize(cpp_molar_density));
   module.def("cpp_molar_enthalpy", py::vectorize(cpp_molar_enthalpy));
   module.def("cpp_dynamic_viscosity", py::vectorize(cpp_dynamic_viscosity));

   py::enum_<Component>(module, "Component")
       .value("single_component", Component::single_component);

   py::enum_<Context>(module, "Context")
       .value("single_context", Context::single_context);

   py::enum_<Phase>(module, "Phase").value("single_phase", Phase::single_phase);

   module.def(
       "build_state",
       [](py::object p, py::object T) {
          if (p.is_none() || T.is_none())
             throw std::runtime_error("You must provide both p and T.");
          X result;
          result.context = static_cast<int>(Context::single_context);
          result.p = p.cast<double>();
          result.T = T.cast<double>();
          result.S.fill(1);    // single phase
          result.C[0][0] = 1;  // single component
          return result;
       },
       py::arg("p") = py::none{}, py::arg("T") = py::none{},
       R"doc(
Construct a state pressure and temperature.

Parameters
----------

:param p: pressure
:param T: temperature

)doc");
}
