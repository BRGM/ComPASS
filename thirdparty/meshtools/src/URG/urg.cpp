#include <array>
#include <vector>

#include <CGAL/Kernel_traits.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Convex_hull_3.h>
//#include <CGAL/Surface_mesh.h>
//#include <CGAL/Polygon_mesh_processing/corefinement.h>
//#include <CGAL/barycenter.h>
//#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>


//FIXME: prefer standard compliant optional
#include <boost/optional.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

//#include "common.h"

#include "On_side.h"
#include "Boundary_hull.h"
#include "Existence_domain.h"
#include "Indexed_plane.h"
#include "conveniencies.h"

// FIXME: std::size_t -> pointer like !
template <typename Handle_type, std::size_t n>
struct Handle {
    Handle_type handle;
    std::array<std::size_t, n> indices;
    template <typename ...Ts>
    Handle(const Handle_type& h, const Ts&... is) :
        handle{ h },
        indices{ { is... } } {
        static_assert(sizeof...(Ts) == n, "wrong number of indices");
        std::sort(begin(indices), end(indices));
    }
};

typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<CGAL::Epeck, CGAL::Default, Itag> CDT;

template <typename Kernel, typename Constraint_range, typename Surface_mesh>
void triangulate(const Boundary_hull<Kernel>& hull, const typename Kernel::Plane_3& plane, const Constraint_range& constraints, Surface_mesh& mesh)
{
    //hull.intersect(plane, mesh);
    auto cdt = CDT{};
    //for (auto&& e : mesh.edges()) {
    //    if(mesh.is_border(e)) {
    //        cdt.insert_constraint(
    //            plane.to_2d(mesh.point(mesh.vertex(e, 0))),
    //            plane.to_2d(mesh.point(mesh.vertex(e, 1)))
    //        );
    //    }
    //}
    for (auto&& S : constraints) {
        cdt.insert_constraint(
            plane.to_2d(S.source()),
            plane.to_2d(S.target())
        );
    }
    std::map<CDT::Vertex_handle, typename Surface_mesh::Vertex_index> vmap;
    mesh.clear();
    for (auto v = cdt.finite_vertices_begin(); v != cdt.finite_vertices_end(); ++v) {
        vmap[v] = mesh.add_vertex(plane.to_3d(v->point()));
    }
    for (auto f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); ++f) {
        mesh.add_face(
            vmap[f->vertex(0)],
            vmap[f->vertex(1)],
            vmap[f->vertex(2)]
        );
    }
}

template <typename Kernel>
struct URGMesh
{
    typedef Boundary_hull<Kernel> Hull;
    typedef typename Hull::Surface_mesh Surface_mesh;
    std::vector<Surface_mesh> surface_meshes;
    template <typename Definitions>
    URGMesh(const Hull& hull, const Definitions& definitions)
    {
        typedef typename Kernel::Point_3 Point;
        typedef typename Kernel::Vector_3 Vector;
        typedef typename Kernel::Line_3 Line;
        typedef typename Kernel::Segment_3 Segment;
        typedef typename Kernel::Plane_3 Plane;
        typedef Indexed_plane<Plane> Indexed_plane;
        typedef Existence_domain<Indexed_plane> Domain;

        std::vector<Point> origin;
        for (auto&& def : definitions) {
            origin.emplace_back(def.first);
        }
        origin.shrink_to_fit();
        std::vector<Plane> plane;
        plane.reserve(origin.size());
        for (auto&& def : definitions) {
            plane.emplace_back(def.first, def.second);
        }

        const auto n = origin.size();
        std::vector<Domain> domain;
        std::cerr << "build domains" << std::endl;
        for (std::size_t j = 0; j < n; ++j) {
            domain.emplace_back();
            auto& Dj = domain.back();
            const auto& Oj = origin[j];
            for (std::size_t i = 0; i < j; ++i) {
                if (domain[i](Oj)) {
                    Dj.add_side(Indexed_plane{ plane[i], i }, Oj);
                }
            }
            //std::cerr << "plane " << j
            //          << " has " << Dj.sides.size()
            //          << " constraints" << std::endl;
        }

        std::cerr << "build plane intersections" << std::endl;

        using Segment_handle = Handle<std::size_t, 2>;
        std::vector< std::vector<Segment_handle> > segment_handles;
        segment_handles.resize(n);
        std::vector<Segment> segments;
        for (std::size_t i = 0; i < n; ++i) {
            const auto& Pi = plane[i];
            const auto& Di = domain[i];
            for (auto&& side : Di.sides) {
                const auto j = side.sided_object.index;
                assert(j < i);
                const auto& Pj = plane[j];
                const auto oL = CGAL::intersection(Pi, Pj);
                if (oL) {
                    if (const Line *pL = boost::get<Line>(&(*oL))) {
                        //std::cerr << "      intersection line " << i << " " << j << std::endl;
                        const auto oS = hull.intersect(*pL);
                        if (oS) {
                            if (const Segment *pS = boost::get<Segment>(&(*oS))) {
                                //std::cerr << "      intersection segment " << i << " " << j << std::endl;
                                if (const auto oSij = Di.clip(*pS)) {
                                    //std::cerr << "      clipped segment " << i << " " << j << std::endl;
                                    auto hSij = Segment_handle{ segments.size(), j, i };
                                    segments.emplace_back(*oSij);
                                    segment_handles[i].push_back(hSij);
                                    segment_handles[j].push_back(hSij);
                                }
                            }
                            else {
                                assert(false);
                            }
                        }
                    }
                }
            }
        }
        std::cerr << "   found " << segments.size() << " segments" << std::endl;

        /*

        std::cerr << "build corners" << std::endl;

        using Corner_handle = Handle<std::size_t, 3>;
        std::vector< std::vector<Corner_handle> > corner_handles;
        corner_handles.resize(n);
        std::vector<Point> corners;
        for (std::size_t i = 0; i < n; ++i) {
        const auto& hS = segment_handles[i];
        const auto nhS = hS.size();
        for (std::size_t hj = 0; hj < nhS; ++hj) {
        for (std::size_t hk = hj+1; hk < nhS; ++hk) {
        assert(hS[hj].indices[0] == i || hS[hj].indices[1] == i);
        assert(hS[hk].indices[0] == i || hS[hk].indices[1] == i);
        const auto& Sj = segments[hS[hj].handle];
        const auto& Sk = segments[hS[hk].handle];
        const auto j = hS[hj].indices[0] == i ? hS[hj].indices[1] : hS[hj].indices[0];
        const auto k = hS[hk].indices[0] == i ? hS[hk].indices[1] : hS[hk].indices[0];
        if (i < j && j < k) {
        auto I = CGAL::intersection(Sj, Sk);
        if (I) {
        if (const Point *P = boost::get<Point>(&*I)) {
        auto hCijk = Corner_handle{ corners.size(), i, j, k };
        corners.emplace_back(*P);
        corner_handles[i].push_back(hCijk);
        corner_handles[j].push_back(hCijk);
        corner_handles[k].push_back(hCijk);
        }
        else {
        assert(false);
        }
        }
        }
        }
        }
        }
        std::cerr << "   found " << corners.size() << " corners" << std::endl;

        */

        surface_meshes.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            const auto& Pi = plane[i];
            auto& mesh = surface_meshes[i];
            hull.intersect(Pi, mesh);
            auto cdt = CDT{};
            for (auto&& e : mesh.edges()) {
                if(mesh.is_border(e)) {
                    auto border_edge = domain[i].clip(Segment{ mesh.point(mesh.vertex(e, 0)), mesh.point(mesh.vertex(e, 1)) });
                    if (border_edge) {
                        //std::cerr << "clipped border edge " << i << std::endl;
                        cdt.insert_constraint(
                            Pi.to_2d(border_edge->source()),
                            Pi.to_2d(border_edge->target())
                        );
                    }
                }
            }
            for (auto&& Sh : segment_handles[i]) {
                const auto& S = segments[Sh.handle];
                cdt.insert_constraint(
                    Pi.to_2d(S.source()),
                    Pi.to_2d(S.target())
                );
            }
            mesh.clear();
            std::map<CDT::Vertex_handle, typename Surface_mesh::Vertex_index> vmap;
            for (auto v = cdt.finite_vertices_begin(); v != cdt.finite_vertices_end(); ++v) {
                vmap[v] = mesh.add_vertex(Pi.to_3d(v->point()));
            }
            for (auto f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); ++f) {
                mesh.add_face(
                    vmap[f->vertex(0)],
                    vmap[f->vertex(1)],
                    vmap[f->vertex(2)]
                );
            }
        }
    }

};

typedef CGAL::Epick Kernel;
typedef typename Kernel::Point_3 Point;
typedef typename Kernel::Vector_3 Vector;
typedef typename Kernel::Plane_3 Plane;
typedef URGMesh<CGAL::Epeck> Mesh;
typedef typename Mesh::Hull Hull;
typedef typename Mesh::Surface_mesh Surface_mesh;

namespace py = pybind11;

auto make_mesh(const Hull& hull, py::list& l)
{
    typedef typename Hull::Kernel Hull_kernel;
    auto to_hull_kernel = CGAL::Cartesian_converter<Kernel, Hull_kernel>{};
    std::vector < std::pair < typename Hull_kernel::Point_3, typename Hull_kernel::Vector_3> > definitions;
    for (auto&& p : l) {
        auto t = p.cast<py::tuple>();
        definitions.emplace_back(
            to_hull_kernel(t[0].cast<Point>()),
            to_hull_kernel(t[1].cast<Vector>())
        );
    }
    return std::make_unique<Mesh>(hull, definitions);
}

PYBIND11_MODULE(URG, module)
{

    module.doc() = "plane intersection for URG model";

    py::class_<Point>(module, "Point")
        .def(py::init<double, double, double>())
        ;

    py::class_<Vector>(module, "Vector")
        .def(py::init<double, double, double>())
        ;

    py::class_<Plane>(module, "Plane")
        .def(py::init<Point, Vector>())
        ;

    py::class_<Surface_mesh>(module, "Surface_mesh")
        .def("as_arrays", &as_numpy_arrays<Surface_mesh>)
        ;

    py::class_<Hull>(module, "Hull")
        .def(py::init([](py::list l) {
        std::vector<Plane> planes;
        planes.reserve(py::len(l));
        for (auto& P : l) {
            planes.emplace_back(P.cast<Plane>());
        }
        return std::make_unique<Hull>(cbegin(planes), cend(planes));
    }))
        .def("boundaries", [](const Hull& self) {
        return py::make_iterator(begin(self.boundaries), end(self.boundaries));
    }, py::keep_alive<0, 1>())
        ;

    py::class_<Mesh>(module, "Mesh")
        .def(py::init(&make_mesh))
        .def("surfaces", [](const Mesh& self) {
        return py::make_iterator(begin(self.surface_meshes), end(self.surface_meshes));
    }, py::keep_alive<0, 1>())
        ;


}
