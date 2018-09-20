#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "C3t3Wrapper.h"
#include "C3t3Wrapper-module.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include "mesh-pyutils.h"
#include "mesh_implicit_domains.h"
#include "implicit_functions.h"

typedef CGAL::Surface_mesh<typename CGAL::Epick::Point_3> Mesh;

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

    // quick and dirty, should be elsewhere
    py::class_<Mesh>(module, "SMesh")
        .def("as_arrays", [](const Mesh& self) {
        return mesh_as_arrays(self);
    });

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
                .def("output_boundaries_to_off", [](C3t3Wrapper& self, typename C3t3Wrapper::Domain_tag subdomain, py::str filename) {
                std::ofstream output(filename);
                self.get().output_boundary_to_off(output, subdomain);
                output.close();
            })
                .def("collect_boundaries", [](C3t3Wrapper& self, typename C3t3Wrapper::Domain_tag subdomain, bool normals_point_outside) {
                Mesh mesh;
                copy_boundary_to_surface_mesh(self.get(), subdomain, mesh, normals_point_outside);
                return mesh;
            }, py::arg("subdomain"), py::arg("normals_point_outside") = true)
                .def("collect_all_boundaries", [](C3t3Wrapper& self, bool normals_point_outside) {
                return collect_all_boundaries_as_surface_meshes<Mesh>(self.get(), normals_point_outside);
            }, py::arg("normals_point_outside") = true)
                ;

            module.def("mesh_implicit_domains_boundaries",
                &mesh_implicit_domains_boundaries);

            module.def("torus", &torus_function);
            module.def("sphere", &sphere_function<3>);

}
