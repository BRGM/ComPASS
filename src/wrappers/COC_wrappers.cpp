//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both the terms
// of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
// and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include "COC.h"

#include "COC_wrappers.h"
#include <pybind11/numpy.h>

void add_coc_wrappers(py::module& module)
{

	py::class_<COC_container>(module, "COCcontainer", py::buffer_protocol())
		.def_buffer([](COC_container &cocc) -> py::buffer_info {
		return py::buffer_info(
			std::begin(cocc),                               /* Pointer to buffer */
			sizeof(int),                          /* Size of one scalar */
			py::format_descriptor<int>::format(), /* Python struct-style format descriptor */
			1,                                      /* Number of dimensions */
			{ cocc.length() },                 /* Buffer dimensions */
			{ sizeof(int) }
		);
	})
		.def("__len__", [](const COC_container &cocc) { return cocc.length(); });

	py::class_<COC_iterator>(module, "COCiterator");

	py::class_<COC>(module, "COC")
		.def("__iter__", [](COC& coc) {
		return py::make_iterator(coc.begin(), coc.end());
	},
			py::keep_alive<0, 1>() /* Keep COC alive while iterator is used */
		)
		.def("__getitem__", [](COC& coc, int i) -> COC_container { return coc[i]; })
		.def("__len__", [](const COC& coc) { return coc.number_of_containers(); })
		.def("offsets", [](py::object object) {
		auto coc = object.cast<COC&>();
		return py::array_t<typename COC::offset_type, py::array::c_style>{
			{ coc.number_of_containers() + 1 }, coc.offset_data(), object
		};
	}, py::keep_alive<0,1>())		
		.def("contiguous_content", [](py::object object) {
		auto coc = object.cast<COC&>();
		return py::array_t<typename COC::value_type, py::array::c_style>{
			{ coc.size() }, coc.content_data(), object
		};
	}, py::keep_alive<0, 1>())
		;

}
