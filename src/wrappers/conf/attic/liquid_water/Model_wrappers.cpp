#include "Model_wrappers.h"

#include <pybind11/numpy.h>

#include <array>

constexpr int NC = ComPASS_NUMBER_OF_COMPONENTS;
constexpr int NP = ComPASS_NUMBER_OF_PHASES;
static_assert(NP == 1, "Wrong number of phases.");
constexpr int SINGLE_PHASE = 0;

// Fortran functions
extern "C" {
void FluidThermodynamics_molar_density(int, double, double, const double *,
                                       double &, double &, double &, double *);
void FluidThermodynamics_molar_enthalpy(int, double, double, const double *,
                                        double &, double &, double &, double *);
void FluidThermodynamics_dynamic_viscosity(int, double, double, const double *,
                                           double &, double &, double &,
                                           double *);
}

void init_model() {}

void finalize_model() {}

template <int PHASE>
inline double phase_molar_density(double p, double T) {
   double xsi, dxsidp, dxsidT;
   double C[NC] = {1};
   double dxsidC[NC] = {0};
   FluidThermodynamics_molar_density(2, p, T, C, xsi, dxsidp, dxsidT, dxsidC);
   return xsi;
}

inline double molar_density(double p, double T) {
   return phase_molar_density<SINGLE_PHASE>(p, T);
}

template <int PHASE>
inline double phase_molar_enthalpy(double p, double T) {
   double h, dhdp, dhdT;
   double C[NC] = {1};
   double dhdC[NC] = {0};
   FluidThermodynamics_molar_enthalpy(2, p, T, C, h, dhdp, dhdT, dhdC);
   return h;
}

inline double molar_enthalpy(double p, double T) {
   return phase_molar_enthalpy<SINGLE_PHASE>(p, T);
}

template <int PHASE>
inline double phase_dynamic_viscosity(double p, double T) {
   double mu, dmudp, dmudT;
   double C[NC] = {1};
   double dmudC[NC] = {0};
   FluidThermodynamics_dynamic_viscosity(2, p, T, C, mu, dmudp, dmudT, dmudC);
   return mu;
}

inline double dynamic_viscosity(double p, double T) {
   return phase_dynamic_viscosity<SINGLE_PHASE>(p, T);
}

void add_model_wrappers(py::module &module) {
   module.def("init_model", &init_model);
   module.def("finalize_model", &finalize_model);
   module.def("molar_density", py::vectorize(molar_density));
   module.def("molar_enthalpy", py::vectorize(molar_enthalpy));
   module.def("dynamic_viscosity", py::vectorize(dynamic_viscosity));
}
