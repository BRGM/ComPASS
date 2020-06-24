//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <pybind11/numpy.h>

#include <array>

#include "Model_wrappers.h"

constexpr int NC = ComPASS_NUMBER_OF_COMPONENTS;
constexpr int NP = ComPASS_NUMBER_OF_PHASES;
// FIXME: assumin liquid phase is the latest phase
constexpr int LIQUID_PHASE = NP - 1;

// Fortran functions
extern "C" {
void FluidThermodynamics_molar_density(int, double, double, const double *,
                                       const double *, double &, double &,
                                       double &, double *, double *);
void FluidThermodynamics_molar_enthalpy(int, double, double, const double *,
                                        const double *, double &, double &,
                                        double &, double *, double *);
void FluidThermodynamics_dynamic_viscosity(int, double, double, const double *,
                                           const double *, double &, double &,
                                           double &, double *, double *);
void FluidThermodynamics_Psat(double, double &, double &);
void FluidThermodynamics_Tsat(double, double &, double &);
void check_array_interop(const double *, double *);
}

inline auto liquid_saturation() {
   auto S = std::array<double, NP>{};  // zero initialization
   S[LIQUID_PHASE] = 1;
   assert(std::accumulate(begin(S), end(S), 0) == 1);
   return S;
}

inline double liquid_molar_density(double p, double T) {
   double xsi, dxsidp, dxsidT;
   double C[NC] = {1};
   auto S = liquid_saturation();
   double dxsidC[NC] = {0};
   auto dxsidS = std::array<double, NP>{};  // zero initialization
   FluidThermodynamics_molar_density(2, p, T, C, S.data(), xsi, dxsidp, dxsidT,
                                     dxsidC, dxsidS.data());
   return xsi;
}

inline double liquid_molar_enthalpy(double p, double T) {
   double h, dhdp, dhdT;
   double C[NC] = {1};
   auto S = liquid_saturation();
   double dhdC[NC] = {0};
   auto dhdS = std::array<double, NP>{};  // zero initialization
   FluidThermodynamics_molar_enthalpy(2, p, T, C, S.data(), h, dhdp, dhdT, dhdC,
                                      dhdS.data());
   return h;
}

inline double liquid_dynamic_viscosity(double p, double T) {
   double mu, dmudp, dmudT;
   double C[NC] = {1};
   auto S = liquid_saturation();
   double dmudC[NC] = {0};
   auto dmudS = std::array<double, NP>{};  // zero initialization
   FluidThermodynamics_dynamic_viscosity(2, p, T, C, S.data(), mu, dmudp, dmudT,
                                         dmudC, dmudS.data());
   return mu;
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

void add_Model_wrappers(py::module &module) {
   module.def("Psat", py::vectorize(Psat));
   module.def("Tsat", py::vectorize(Tsat));
   module.def("liquid_molar_density", py::vectorize(liquid_molar_density));
   module.def("liquid_molar_enthalpy", py::vectorize(liquid_molar_enthalpy));
   module.def("liquid_dynamic_viscosity",
              py::vectorize(liquid_dynamic_viscosity));
}
