//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both the terms
// of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
// and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

extern "C"
{
    // CHECKME: We rely on pybind11 access to Py_ssize_t 
    //          which is a signed integral type such that
    //          sizeof(Py_ssize_t) == sizeof(size_t).
    //          See pybind11 code and PEP 353 for details.
    void pass_and_dump_dim_array(const double *, ComPASS_Fortran_size_type, const ComPASS_Fortran_size_type*);
    void increment_first_column(double *, ComPASS_Fortran_size_type, const ComPASS_Fortran_size_type*);
}

#include <pybind11/numpy.h>

// Not used to be tested...

void add_wrappers(py::module& module)
{

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
