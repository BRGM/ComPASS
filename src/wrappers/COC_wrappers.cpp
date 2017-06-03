#include "COC.h"

#include "COC_wrappers.h"

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
	;

}