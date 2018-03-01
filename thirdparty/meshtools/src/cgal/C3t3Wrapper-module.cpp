#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "C3t3Wrapper.h"

namespace py = pybind11;

template <typename Buffer>
auto to_array(const Buffer& buffer)
{
    typedef typename Buffer::value_type Scalar;
    return py::array_t<Scalar, py::array::c_style> {
        buffer.shape, buffer.stride, buffer.data
    };
}

PYBIND11_MODULE(C3t3Wrapper, module)
{

    module.doc() = "pybind11 homemade CGAL C3t3 interface";

	py::class_<C3t3Wrapper>(module, "C3t3")
        .def(py::init<>())
        .def(py::init([](const std::string& filename)
            {
                auto wrapper = std::make_unique<C3t3Wrapper>();
                wrapper->reload(filename);
                return wrapper;
            }
        ))
        .def("reload", &C3t3Wrapper::reload)
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
