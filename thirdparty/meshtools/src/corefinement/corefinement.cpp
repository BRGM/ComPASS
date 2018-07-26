#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

typedef CGAL::Epick Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

auto mesh_as_arrays(const Mesh& mesh)
{
    typedef std::size_t New_index;
    std::vector<New_index> reindex;
    reindex.resize(mesh.num_vertices());
    auto nv = New_index{ 0 };
    for (auto&& v : mesh.vertices()) {
        reindex[v] = nv++;
    }
    auto vertices = py::array_t<double, py::array::c_style>{ 
        { static_cast<std::size_t>(nv), static_cast<std::size_t>(Point::Ambient_dimension::value) }
    };
    {
        auto p = reinterpret_cast<Point*>(vertices.request().ptr);
        static_assert(sizeof(Point) == Point::Ambient_dimension::value * sizeof(double), "Inconsistent sizes in memory!");
        for (auto&& v : mesh.vertices()) {
            *(p++) = mesh.point(v);
        }
    }
    auto triangles = py::array_t<New_index, py::array::c_style>{ 
        { static_cast<std::size_t>(mesh.number_of_faces()), static_cast<std::size_t>(3) }
    };
    {
        auto p = reinterpret_cast<New_index*>(triangles.request().ptr);
        for (auto&& f : mesh.faces()) {
            assert(mesh.degree(f) == 3);
            for (auto&& v : CGAL::vertices_around_face(mesh.halfedge(f), mesh)) {
                *(p++) = reindex[v];
            }
        }
    }
    return py::make_tuple(vertices, triangles);
}

auto test()
{
    Mesh tm1;
    Mesh::Vertex_index u = tm1.add_vertex(Point(0, 1, 0));
    Mesh::Vertex_index v = tm1.add_vertex(Point(0, 0, 0.2));
    Mesh::Vertex_index w = tm1.add_vertex(Point(1, 0, 0));
    tm1.add_face(u, v, w);
    Mesh tm2;
    u = tm2.add_vertex(Point(  0, 0.4, -1));
    v = tm2.add_vertex(Point(  1, 0.5, -1));
    w = tm2.add_vertex(Point(0.5, 0.6,  1));
    tm2.add_face(u, v, w);
    CGAL::Polygon_mesh_processing::corefine(tm1, tm2, true);
    std::ofstream os("test.off");
    CGAL::write_off(os, tm1);
    os.close();
    return py::make_tuple(mesh_as_arrays(tm1), mesh_as_arrays(tm2));
}

PYBIND11_MODULE(Corefinement, module)
{

    module.doc() = "pybind11 with CGAL surface corefinement (quick and dirty!!!)";

    module.def("test", &test);

}
