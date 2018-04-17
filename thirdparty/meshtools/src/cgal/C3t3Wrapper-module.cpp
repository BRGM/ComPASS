#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "C3t3Wrapper.h"
#include "C3t3Wrapper-module.h"

namespace py = pybind11;

template <typename Buffer>
auto to_array(const Buffer& buffer)
{
    typedef typename Buffer::value_type Scalar;
    return py::array_t<Scalar, py::array::c_style> {
        buffer.shape, buffer.stride, buffer.data
    };
}

void add_c3t3_wrapper(py::module& module)
{

    module.doc() = "pybind11 homemade CGAL C3t3 interface";

	py::class_<C3t3Wrapper>(module, "C3t3")
        .def(py::init([](const std::string& filename, bool binary)
            {
                auto wrapper = std::make_unique<C3t3Wrapper>();
                if (binary) wrapper->binary_reload(filename);
                else wrapper->reload(filename);
                return wrapper;
            }
        ), py::arg("filename"), py::arg("binary") = true)
        .def("reload", &C3t3Wrapper::reload)
                .def("binary_reload", &C3t3Wrapper::binary_reload)
                // default return value policy for def_property is keep_alive<0,1>
                .def_property_readonly("vertices", [](C3t3Wrapper& self) {
                return to_array(self.vertices_buffer());
            })
                .def_property_readonly("cells", [](C3t3Wrapper& self) {
                return to_array(self.cells_buffer());
            })
                .def_property_readonly("domain_tags", [](C3t3Wrapper& self) {
                return to_array(self.domain_tags_buffer());
            })
                .def_property_readonly("facets", [](C3t3Wrapper& self) {
                return to_array(self.facets_buffer());
            })
                .def_property_readonly("facet_tags", [](C3t3Wrapper& self) {
                return to_array(self.facet_tags_buffer());
            })
        ;

}
