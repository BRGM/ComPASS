#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include "TSurfBlob.h"

typedef TSurfBlob<typename CGAL::Epick::Point_3> Blob;
typedef typename Blob::Point Point;
typedef typename Blob::TSurf Mesh;

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "mesh-pyutils.h"
#include "collect_consrained_edges_as_curves.h"

namespace py = pybind11;

auto corefine_surfaces(Blob& blob, Mesh& S1, Mesh& S2)
{

    typedef Mesh::Vertex_index Vertex_index;
    typedef Mesh::Edge_index Edge_index;

    auto add_constraint_map = [](Mesh& mesh) {
        auto result = mesh.add_property_map<Edge_index, bool>("e:constrained", false);
        assert(result.second);
        return result.first;
    };
    auto constraints1 = add_constraint_map(S1);
    auto constraints2 = add_constraint_map(S2);
    namespace parameters = CGAL::Polygon_mesh_processing::parameters;
    CGAL::Polygon_mesh_processing::corefine(
        S1, S2,
        parameters::edge_is_constrained_map(constraints1),
        parameters::edge_is_constrained_map(constraints2),
        true
    );
    auto curves1 = collect_consrained_edges_as_curves(S1, constraints1);
    //auto curves2 = collect_curves(S2, constraints2);
    S1.remove_property_map(constraints1);
    S2.remove_property_map(constraints2);
/**
    // ------------------------------------------------------------------
    std::cout << "curve 1 ------------------------" << std::endl;
    for (auto&& curve : curves1) {
        std::cout << "-- curve" << std::endl;
        for (auto&& v : curve) {
            std::cout << "  " << S1.point(v) << std::endl;
        }
    }
    std::cout << "curve 2 ------------------------" << std::endl;
    for (auto&& curve : curves2) {
        std::cout << "-- curve" << std::endl;
        for (auto&& v : curve) {
            std::cout << "  " << S2.point(v) << std::endl;
        }
    }
/**/
    //associate_curves(S1, curves1, S2, curves2);
    for (auto&& curve : curves1) {
        auto p = blob.new_intersection(&S1, &S2);
        for (auto&& v1 : curve) {
            assert(S1.is_valid(v1));
            p->emplace_back(S1.point(v1));
        }
    }
    std::cout << "Number of blob curves: " << blob.curves_factory.already_created() << std::endl;
}

auto test()
{
    typedef typename CGAL::Epick::Point_3 Point;
    Mesh tm1;
    Mesh::Vertex_index u = tm1.add_vertex(Point(0, 1, 0));
    Mesh::Vertex_index v = tm1.add_vertex(Point(0, 0, 0.2));
    Mesh::Vertex_index w = tm1.add_vertex(Point(1, 0, 0));
    tm1.add_face(u, v, w);
    Mesh tm2;
    u = tm2.add_vertex(Point(0, 0.4, -1));
    v = tm2.add_vertex(Point(1, 0.5, -1));
    w = tm2.add_vertex(Point(0.5, 0.6, 1));
    Mesh::Vertex_index x = tm2.add_vertex(Point(-0.5, 0.8, 1));
    tm2.add_face(u, v, w);
    tm2.add_face(u, w, x); // beware of face orientation - otherwise face is not added
    Mesh tm3;
    u = tm3.add_vertex(Point(-1.5, 1.1, 1));
    v = tm3.add_vertex(Point(0.8, 1.2, -1.3));
    w = tm3.add_vertex(Point(0.7, -1., 1.5));
    tm3.add_face(u, v, w);

    Blob blob;
    corefine_surfaces(blob, tm1, tm2);
    corefine_surfaces(blob, tm1, tm3);
    corefine_surfaces(blob, tm2, tm3);

    //std::ofstream os("test.off");
    //CGAL::write_off(os, tm1);
    //os.close();
    return py::make_tuple(mesh_as_arrays(tm1), mesh_as_arrays(tm2), mesh_as_arrays(tm3));

}


template <typename Functor>
py::object apply_functor_on_array(Functor F,
    py::array_t<double, py::array::c_style> a)
{
    typedef Point Arg;
    typedef double DType;
    constexpr std::size_t dim = 3;
    static_assert(sizeof(Arg) == dim * sizeof(DType), "Inconsistent sizes in memory.");
    assert(a.size() % dim == 0);
    const std::size_t n = a.size() / dim;
    // we do NOT want to alter the input array shape so we work on  
    auto p = reinterpret_cast<const Arg *>(a.request().ptr);
    if (n == 1) {
        return py::float_{ F(*p) }; // returns a scalar wrapped as py::object
    }
    auto result = py::array_t<DType, py::array::c_style>(n);
    auto res = result.template mutable_unchecked<1>();
    for (std::size_t i = 0; i < n; ++i) {
        *(res.mutable_data(i)) = F(*p);
        ++p;
    }
    return result; // returns an array wrapped as py::object
}



inline auto f1(const Point& P)
{
    return P.z() - (P.x()*P.x() + P.y()*P.y());
}

inline auto f2(const Point& P)
{
    return P.z() + (P.x()*P.x() + P.y()*P.y());
}

PYBIND11_MODULE(Corefinement, module)
{

    module.doc() = "pybind11 with CGAL surface corefinement (quick and dirty!!!)";

    module.def("test", &test);
    module.def("f1", [](py::array_t<double, py::array::c_style> a) {
        return apply_functor_on_array(f1, a);
    });
    module.def("f2", [](py::array_t<double, py::array::c_style> a) {
        return apply_functor_on_array(f2, a);
    });

}
