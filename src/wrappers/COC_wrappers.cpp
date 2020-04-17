//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both the terms
// of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
// and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include "COC_wrappers.h"
#include "COC_exposers.h"

void add_coc_wrappers(py::module& module)
{

	auto exposed = expose_coc<int, true>(module, "");
	auto& container_class = std::get<2>(exposed);

	container_class.def_buffer([](COC_container &cocc) -> py::buffer_info {
		return py::buffer_info(
			std::begin(cocc),                               /* Pointer to buffer */
			sizeof(int),                          /* Size of one scalar */
			py::format_descriptor<int>::format(), /* Python struct-style format descriptor */
			1,                                      /* Number of dimensions */
			{ cocc.length() },                 /* Buffer dimensions */
			{ sizeof(int) }
		);
	})
	;

}
