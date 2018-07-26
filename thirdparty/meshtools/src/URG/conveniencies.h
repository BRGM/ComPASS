#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "common.h"

namespace py = pybind11;

template <typename T>
struct Array_view
{
    typedef std::remove_cv_t<T> value_type;
    T * begin_;
    T * end_;
    template <typename TT>
    Array_view(py::array_t<TT, py::array::c_style>& a) :
        begin_{ reinterpret_cast<T*>(a.data()) },
        end_{ begin_ } {
        static_assert(sizeof(T) % sizeof(TT) == 0, "inconsistent size in memory");
        assert(a.size() % (sizeof(T) / sizeof(TT)) == 0);
        end_ += a.size() / (sizeof(T) / sizeof(TT));
    }
    auto begin() { return begin_; }
    auto end() { return end_; }
    auto begin() const { return begin_; }
    auto end() const { return end_; }
};

template <typename Surface_mesh>
auto as_numpy_arrays(const Surface_mesh& mesh) -> py::tuple
{

    typedef typename Surface_mesh::Point Surface_point;
    typedef typename CGAL::Kernel_traits<Surface_point>::Kernel Surface_kernel;
    typedef CGAL::Epick::Point_3 Point;
    auto to_epick = CGAL::Cartesian_converter<Surface_kernel, CGAL::Epick>{};

    typedef MeshTools::ElementId New_index;
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
            *(p++) = to_epick(mesh.point(v));
        }
    }
    auto triangles = py::array_t<New_index, py::array::c_style>{
        { static_cast<std::size_t>(mesh.number_of_faces()), static_cast<std::size_t>(3) }
    };
    {
        auto p = reinterpret_cast<New_index*>(triangles.request().ptr);
        for (auto&& f : mesh.faces()) {
            //std::cerr << "Dumping: " << f << "(" << mesh.degree(f) << ")" << std::endl;
            assert(mesh.degree(f) == 3);
            for (auto&& v : CGAL::vertices_around_face(mesh.halfedge(f), mesh)) {
                *(p++) = reindex[v];
            }
        }
    }
    return py::make_tuple(vertices, triangles);
}

