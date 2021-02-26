#include <pybind11/numpy.h>

#include "DefModel.h"
#include "Model_wrappers.h"
#include "Thermodynamics.h"

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
}
