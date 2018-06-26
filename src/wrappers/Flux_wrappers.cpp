//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both the terms
// of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
// and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include "Flux_wrappers.h"
#include <pybind11/numpy.h>

struct cpp_narray_wrapper
{
    const std::size_t * shape;
    const double * data;
};

// Fortran functions
extern "C"
{
    void retrieve_mass_fluxes(const cpp_narray_wrapper&, const cpp_narray_wrapper&);
    std::size_t number_of_nodes();
    std::size_t number_of_cells();
    std::size_t number_of_fractures();
}

void add_flux_wrappers(py::module& module)
{

    constexpr std::size_t NC = ComPASS_NUMBER_OF_COMPONENTS;
    constexpr std::size_t NP = ComPASS_NUMBER_OF_PHASES;

    module.def("retrieve_mass_fluxes", [] {
        const auto Fcs = std::vector<std::size_t>{ 3, NC, number_of_cells() };
        auto Fc = py::array_t<double, py::array::c_style>{ Fcs };
        const auto Ffs = std::vector<std::size_t>{ 3, NC, number_of_fractures() };
        auto Ff = py::array_t<double, py::array::c_style>{ Ffs };
        retrieve_mass_fluxes(
            {Fcs.data(), reinterpret_cast<const double *>(Fc.request().ptr) },
            {Ffs.data(), reinterpret_cast<const double *>(Ff.request().ptr) }
        );
        return py::make_tuple(Fc, Ff);
    });

}
