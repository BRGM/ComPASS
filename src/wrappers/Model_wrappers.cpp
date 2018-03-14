//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both the terms
// of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
// and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include "Model_wrappers.h"

#include <array>

#include <pybind11/numpy.h>

constexpr int NC = ComPASS_NUMBER_OF_COMPONENTS;
constexpr int NP = ComPASS_NUMBER_OF_PHASES;

// Fortran functions
extern "C"
{
    void FluidThermodynamics_molar_density(int, double, double, const double *, const double *, double&, double&, double&, double *, double *);
    void FluidThermodynamics_molar_enthalpy(int, double, double, const double *, const double *, double&, double&, double&, double *, double *);
    void FluidThermodynamics_Psat(double, double&, double&);
    void FluidThermodynamics_Tsat(double, double&, double&);
    void check_array_interop(const double*, double*);
}

inline double liquid_molar_density(double p, double T)
{
    double xsi, dxsidp, dxsidT;
    double C[NC] = { 1 };
    double S[NP] = { 0, 1 };
    double dxsidC[NC] = { 0 };
    double dxsidS[NP] = { 0, 0 };
    FluidThermodynamics_molar_density(2, p, T, C, S, xsi, dxsidp, dxsidT, dxsidC, dxsidS);
    return xsi;
}

inline double liquid_molar_enthalpy(double p, double T)
{
    double h, dhdp, dhdT;
    double C[NC] = { 1 };
    double S[NP] = { 0, 1 };
    double dhdC[NC] = { 0 };
    double dhdS[NP] = { 0, 0 };
    FluidThermodynamics_molar_enthalpy(2, p, T, C, S, h, dhdp, dhdT, dhdC, dhdS);
    return h;
}

inline double Psat(double T)
{
    double result;
    double dPsatdT;
    FluidThermodynamics_Psat(T, result, dPsatdT);
    return result;
}


inline double Tsat(double p)
{
    double result;
    double dTsatdp;
    FluidThermodynamics_Tsat(p, result, dTsatdp);
    return result;
}

void add_Model_wrappers(py::module& module)
{

    module.def("Psat", py::vectorize(Psat));
    module.def("Tsat", py::vectorize(Tsat));
    module.def("liquid_molar_density", py::vectorize(liquid_molar_density));
    module.def("liquid_molar_enthalpy", py::vectorize(liquid_molar_enthalpy));

}
