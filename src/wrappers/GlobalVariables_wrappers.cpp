//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both the terms
// of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
// and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

// FIXME: should be defined with common types/defines
//        in ComPASS namespace
typedef int ComPASS_Fortran_size_type;

// Fortran functions
extern "C"
{
    double get_current_time();
    void set_current_time(double);
    double get_delta_t();
    double get_final_time();
    double get_initial_timestep();
    double get_maximum_timestep();
    double get_gravity();
    void set_gravity(double);
    double get_rock_volumetric_heat_capacity();
    void set_rock_volumetric_heat_capacity(double);
    double get_fracture_thickness();
    void set_fracture_thickness(double);
    void set_final_time(double);
    void set_initial_timestep(double);
    void set_maximum_timestep(double);
    // CHECKME: We rely on pybind11 access to Py_ssize_t 
    //          which is a signed integral type such that
    //          sizeof(Py_ssize_t) == sizeof(size_t).
    //          See pybind11 code and PEP 353 for details.
    void pass_and_dump_dim_array(const double *, ComPASS_Fortran_size_type, const ComPASS_Fortran_size_type*);
    void increment_first_column(double *, ComPASS_Fortran_size_type, const ComPASS_Fortran_size_type*);
}

#include "GlobalVariables_wrappers.h"
#include <pybind11/numpy.h>

void add_global_variables_wrappers(py::module& module)
{

    // module.def("get_current_time", &get_current_time);
    // module.def("set_current_time", &set_current_time);
    // module.def("get_timestep", &get_delta_t);
    // module.def("get_final_time", &get_final_time);
    // module.def("set_final_time", &set_final_time);
    // module.def("get_initial_timestep", &get_initial_timestep);
    // module.def("set_initial_timestep", &set_initial_timestep);
    // module.def("get_maximum_timestep", &get_maximum_timestep);
    // module.def("set_maximum_timestep", &set_maximum_timestep);
    module.def("get_gravity", &get_gravity);
    module.def("set_gravity", &set_gravity);
    module.def("get_rock_volumetric_heat_capacity", &get_rock_volumetric_heat_capacity);
    module.def("set_rock_volumetric_heat_capacity", &set_rock_volumetric_heat_capacity);
    module.def("get_fracture_thickness", &get_fracture_thickness);
    module.def("set_fracture_thickness", &set_fracture_thickness);
    module.def("has_energy_transfer_enabled", []() {
#ifdef _THERMIQUE_
        return true;
#else
        return false;
#endif
    });
    
/**
 CHECKME: we could pass py::array_t::shape() which is a pointer
          to the underlying shape yet we would be wrapping through an array of 
          C unsigned longs which have no equivalent in Fortran
          and we should have a compilation time check that 
          sizeof(Py_ssize_t) == sizeof(size_t)
          == the kind of Fortran integer used here
          c_long_long => KIND=8
          furthermore we need to reverse the C array shape to access colums
          instead of rows in Fortran and vice versa
*/

    auto fortran_shape = [](py::array_t<double, py::array::c_style>& a) {
        std::vector<ComPASS_Fortran_size_type> res;
        res.reserve(a.ndim());
        for(int i=a.ndim(); i>0; --i) {
            res.push_back( a.shape(i-1) );
        }
        return res;
    };
    
    module.def("dump_array_in_fortran",
        [fortran_shape](py::array_t<double, py::array::c_style>& a) {
            auto fshape = fortran_shape(a);
            pass_and_dump_dim_array(a.data(), a.ndim(), fshape.data());
        }
    );

    module.def("increment_first_column_in_fortran",
        [fortran_shape](py::array_t<double, py::array::c_style>& a) {
            auto fshape = fortran_shape(a);
            increment_first_column(a.mutable_data(), a.ndim(), fshape.data());
        }
    );

}
