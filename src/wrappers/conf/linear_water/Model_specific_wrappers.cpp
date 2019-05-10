#include "Model_wrappers.h"

#include <pybind11/numpy.h>

struct Fluid_properties
{
    double specific_mass;
    double compressibility;
    double thermal_expansivity;
    double volumetric_heat_capacity;
    double dynamic_viscosity;
    double reference_pressure;
    double reference_temperature;
};

// Fortran functions
extern "C"
{
    Fluid_properties * get_fluid_properties();
    void FluidThermodynamics_molar_density(int, double, double, const double *, const double *, double&, double&, double&, double *, double *);
    void FluidThermodynamics_molar_enthalpy(int, double, double, const double *, const double *, double&, double&, double&, double *, double *);
    void FluidThermodynamics_dynamic_viscosity(int, double, double, const double *, const double *, double&, double&, double&, double *, double *);
}

void init_model() 
{
    auto properties = get_fluid_properties();
    properties->specific_mass = 1;
    properties->compressibility = 0;
    properties->thermal_expansivity = 0;
    properties->volumetric_heat_capacity = 1;
    properties->dynamic_viscosity = 1;
    properties->reference_pressure = 1E5 ; // 1 bar
    properties->reference_temperature = 293.15 ; // 20°C
}

void finalize_model() {}

template <typename Function>
inline double call_physical_function(Function function, double p, double T)
{
    constexpr int NC = ComPASS_NUMBER_OF_COMPONENTS;
    constexpr int NP = ComPASS_NUMBER_OF_PHASES;
    static_assert(NC == 1, "we assume there is only one component");
    static_assert(NP == 1, "we assume there is only one phase");
    double f, dfdp, dfdT;
    double C[NC] = { 1 };
    double S[NP] = { 1 };
    double dfdC[NC] = { 0 };
    double dfdS[NP] = { 0 };
    function(1, f, T, C, S, f, dfdp, dfdT, dfdC, dfdS);
    return f;
}

inline double molar_density(double p, double T)
{
    return call_physical_function(FluidThermodynamics_molar_density, p, T);
}

inline double molar_enthalpy(double p, double T)
{
    return call_physical_function(FluidThermodynamics_molar_enthalpy, p, T);
}

inline double dynamic_viscosity(double p, double T)
{
    return call_physical_function(FluidThermodynamics_dynamic_viscosity, p, T);
}

void add_specific_model_wrappers(py::module& module)
{

    py::class_<Fluid_properties>(module, "FluidProperties")
        .def_readwrite("specific_mass", &Fluid_properties::specific_mass)
        .def_readwrite("compressibility", &Fluid_properties::compressibility)
        .def_readwrite("thermal_expansivity", &Fluid_properties::thermal_expansivity)
        .def_readwrite("volumetric_heat_capacity", &Fluid_properties::volumetric_heat_capacity)
        .def_readwrite("dynamic_viscosity", &Fluid_properties::dynamic_viscosity)
        .def_readwrite("reference_pressure", &Fluid_properties::reference_pressure)
        .def_readwrite("reference_temperature", &Fluid_properties::reference_temperature)
        ;

    module.def("get_fluid_properties", &get_fluid_properties, py::return_value_policy::reference);

    // provided mainly for consistency... we are filling arrays with constant scalars...
    module.def("molar_density", py::vectorize(molar_density));
    module.def("molar_enthalpy", py::vectorize(molar_enthalpy));
    module.def("dynamic_viscosity", py::vectorize(dynamic_viscosity));

}
