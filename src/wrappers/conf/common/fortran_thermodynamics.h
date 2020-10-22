#pragma once

extern "C" {
void FluidThermodynamics_molar_density(const int &, const double &,
                                       const double &, const double *, double &,
                                       double &, double &, double *);
void FluidThermodynamics_molar_enthalpy(const int &, const double &,
                                        const double &, const double *,
                                        double &, double &, double &, double *);
void FluidThermodynamics_dynamic_viscosity(const int &, const double &,
                                           const double &, const double *,
                                           double &, double &, double &,
                                           double *);
void FluidThermodynamics_Psat(const double &, double &, double &);
void FluidThermodynamics_Tsat(const double &, double &, double &);
}
