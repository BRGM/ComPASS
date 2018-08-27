#include "Model_wrappers.h"

struct Fluid_properties
{
    double specific_mass;
    double compressibility;
    double thermal_expansivity;
    double specific_enthalpy;
    double viscosity;
};

// Fortran functions
extern "C"
{
    Fluid_properties * get_fluid_properties();
}

void init_model() 
{
    auto properties = get_fluid_properties();
    properties->specific_mass = 1;
    properties->compressibility = 0;
    properties->thermal_expansivity = 0;
    properties->specific_enthalpy = 1;
    properties->viscosity = 1;
}

void finalize_model() {}

void add_model_wrappers(py::module& module)
{

    module.def("init_model", &init_model);
    module.def("finalize_model", &finalize_model);

    py::class_<Fluid_properties>(module, "FluidProperties")
        .def_readwrite("specific_mass", &Fluid_properties::specific_mass)
        .def_readwrite("compressibility", &Fluid_properties::compressibility)
        .def_readwrite("thermal_expansivity", &Fluid_properties::thermal_expansivity)
        .def_readwrite("specific_enthalpy", &Fluid_properties::specific_enthalpy)
        .def_readwrite("viscosity", &Fluid_properties::viscosity)
        ;

    module.def("get_fluid_properties", &get_fluid_properties, py::return_value_policy::reference);

}
