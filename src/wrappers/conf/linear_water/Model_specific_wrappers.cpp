#include <pybind11/numpy.h>

#include "DefModel.h"
#include "Model_wrappers.h"
#include "StateObjects.h"
#include "Thermodynamics.h"

struct Fluid_properties {
   double specific_mass;
   double compressibility;
   double thermal_expansivity;
   double volumetric_heat_capacity;
   double reference_pressure;
   double reference_temperature;
};

// Fortran functions
extern "C" {
Fluid_properties *get_fluid_properties();
}

void init_model() {
   auto properties = get_fluid_properties();
   properties->specific_mass = 1;
   properties->compressibility = 0;
   properties->thermal_expansivity = 0;
   properties->volumetric_heat_capacity = 1;
   properties->reference_pressure = 1E5;        // 1 bar
   properties->reference_temperature = 293.15;  // 20°C
}

void finalize_model() {}

template <typename Function>
inline double call_physical_subroutine(Function function, double p, double T) {
   constexpr int NC = ComPASS_NUMBER_OF_COMPONENTS;
   constexpr int NP = ComPASS_NUMBER_OF_PHASES;
   static_assert(NC == 1, "we assume there is only one component");
   static_assert(NP == 1, "we assume there is only one phase");
   double f, dfdp, dfdT;
   double C[NC] = {1};
   double dfdC[NC] = {0};
   function(1, p, T, C, f, dfdp, dfdT, dfdC);
   return f;
}

template <typename Function>
inline double call_physical_function(Function function, double p, double T) {
   constexpr int NC = ComPASS_NUMBER_OF_COMPONENTS;
   constexpr int NP = ComPASS_NUMBER_OF_PHASES;
   static_assert(NC == 1, "we assume there is only one component");
   static_assert(NP == 1, "we assume there is only one phase");
   double C[NC] = {1};
   return function(1, p, T, C);
}

inline double molar_density(double p, double T) {
   return call_physical_subroutine(FluidThermodynamics_molar_density, p, T);
}

inline double molar_enthalpy(double p, double T) {
   return call_physical_subroutine(FluidThermodynamics_molar_enthalpy, p, T);
}

inline double cpp_dynamic_viscosity(double p, double T) {
   return call_physical_function(FluidThermodynamics_dynamic_viscosity, p, T);
}

void add_specific_model_wrappers(py::module &module) {
   py::class_<Fluid_properties>(module, "FluidProperties")
       .def_readwrite("specific_mass", &Fluid_properties::specific_mass)
       .def_readwrite("compressibility", &Fluid_properties::compressibility)
       .def_readwrite("thermal_expansivity",
                      &Fluid_properties::thermal_expansivity)
       .def_readwrite("volumetric_heat_capacity",
                      &Fluid_properties::volumetric_heat_capacity)
       .def_readwrite("reference_pressure",
                      &Fluid_properties::reference_pressure)
       .def_readwrite("reference_temperature",
                      &Fluid_properties::reference_temperature);

   module.def("get_fluid_properties", &get_fluid_properties,
              py::return_value_policy::reference);

   // provided mainly for consistency... we are filling arrays with constant
   // scalars...
   module.def("molar_density", py::vectorize(molar_density));
   module.def("molar_enthalpy", py::vectorize(molar_enthalpy));
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
