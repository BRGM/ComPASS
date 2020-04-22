//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both the terms
// of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
// and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include "VAGFrac_wrappers.h"


// Fortran functions
extern "C"
{
    void VAGFrac_TransDarcy();
    void VAGFrac_TransFourier();
    void VAGFrac_VolsDarcy(double, double);
    void VAGFrac_VolsFourier(double, double);
}


void add_VAGFrac_wrappers(py::module& module)
{

    module.def("VAGFrac_TransDarcy", &VAGFrac_TransDarcy);
    module.def("VAGFrac_TransFourier", &VAGFrac_TransFourier);
    module.def("VAGFrac_VolsDarcy", &VAGFrac_VolsDarcy);
    module.def("VAGFrac_VolsFourier", &VAGFrac_VolsFourier);

}
