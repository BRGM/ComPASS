#include <array>
#include <vector>

#include <CGAL/Kernel_traits.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Delaunay_triangulation_2.h>
//#include <CGAL/Delaunay_mesher_2.h>
//#include <CGAL/Delaunay_mesh_face_base_2.h>
//#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

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

constexpr auto edge_constraints_map_name = "e:is_constrained";

template <typename Surface_mesh>
auto add_edge_constraints_map(Surface_mesh& mesh) {
    typedef typename Surface_mesh::Edge_index Edge_index;
    auto pmap = mesh.template add_property_map<
        Edge_index, bool
    >(edge_constraints_map_name, false);
    assert(pmap.second);
    return pmap.first;
};

template <typename Surface_mesh>
auto get_edge_constraints_map(Surface_mesh& mesh) {
    typedef typename Surface_mesh::Edge_index Edge_index;
    auto pmap = mesh.template property_map<
        Edge_index, bool
    >(edge_constraints_map_name);
    assert(pmap.second);
    return pmap.first;
};

template <typename Triangulation2, typename Surface_mesh, typename Converter>
auto triangulation_to_surface_mesh_vmap(const Triangulation2& triangulation, Surface_mesh& mesh, const Converter& convert)
{
    std::map<typename Triangulation2::Vertex_handle, typename Surface_mesh::Vertex_index> vmap;
    for (auto v = triangulation.finite_vertices_begin(); v != triangulation.finite_vertices_end(); ++v) {
        vmap[v] = mesh.add_vertex(convert(v->point()));
    }
    return vmap;
}

template <typename Triangulation2, typename VMap, typename Surface_mesh>
void triangulation_to_surface_mesh_add_constraints(const Triangulation2& triangulation, const VMap& vmap, Surface_mesh& mesh)
{
    auto is_constrained_edge = get_edge_constraints_map(mesh);
    for (auto ce = triangulation.constrained_edges_begin(); ce != triangulation.constrained_edges_end(); ++ce) {
        const auto v1 = vmap.at(ce->first->vertex((ce->second + 1) % 3));
        const auto v2 = vmap.at(ce->first->vertex((ce->second + 2) % 3));
        assert(mesh.is_valid(v1));
        auto h = mesh.null_halfedge();
        for (auto&& h2 : CGAL::halfedges_around_source(v1, mesh)) {
            if (mesh.target(h2) == v2) {
                h = h2;
                break;
            }
        }
        assert(h != mesh.null_halfedge());
        is_constrained_edge[mesh.edge(h)] = true;
    }
}

template <typename Triangulation2, typename VMap, typename Surface_mesh>
void triangulation_to_surface_mesh_add_faces(const Triangulation2& triangulation, const VMap& vmap, Surface_mesh& mesh)
{
    for (auto f = triangulation.finite_faces_begin(); f != triangulation.finite_faces_end(); ++f) {
        mesh.add_face(
            vmap.at(f->vertex(0)),
            vmap.at(f->vertex(1)),
            vmap.at(f->vertex(2))
        );
    }
}

template <typename Triangulation2, typename Surface_mesh, typename Converter>
void triangulation_to_surface_mesh(const Triangulation2& triangulation, Surface_mesh& mesh, const Converter& convert)
{
    mesh.clear();
    auto vmap = triangulation_to_surface_mesh_vmap(triangulation, mesh, convert);
    triangulation_to_surface_mesh_add_faces(triangulation, vmap, mesh);
}

template <typename Triangulation2, typename Surface_mesh, typename Converter>
void triangulation_to_surface_mesh_with_constraints(const Triangulation2& triangulation, Surface_mesh& mesh, const Converter& convert)
{
    mesh.clear();
    auto vmap = triangulation_to_surface_mesh_vmap(triangulation, mesh, convert);
    triangulation_to_surface_mesh_add_faces(triangulation, vmap, mesh);
    triangulation_to_surface_mesh_add_constraints(triangulation, vmap, mesh);
}

struct Default_converter
{
    template <typename T>
    auto operator()(const T& t) const noexcept {
        return t;
    }
};

template <typename Triangulation2, typename Surface_mesh>
void triangulation_to_surface_mesh(const Triangulation2& triangulation, Surface_mesh& mesh)
{
    triangulation_to_surface_mesh(triangulation, mesh, Default_converter{});
}

template <typename Triangulation2, typename Surface_mesh>
void triangulation_to_surface_mesh_with_constraints(const Triangulation2& triangulation, Surface_mesh& mesh)
{
    triangulation_to_surface_mesh_with_constraints(triangulation, mesh, Default_converter{});
}

template <typename PointRange, typename Surface_mesh>
void elevation_surface(const PointRange& points, Surface_mesh& mesh)
{

    typedef typename PointRange::value_type Point;
    typedef typename CGAL::Kernel_traits<Point>::Kernel Point_kernel;
    typedef CGAL::Projection_traits_xy_3<Point_kernel>  Geometric_traits;
    typedef CGAL::Delaunay_triangulation_2<Geometric_traits> Delaunay;
    typedef typename Surface_mesh::Point Surface_point;
    typedef typename CGAL::Kernel_traits<Surface_point>::Kernel Surface_kernel;
    
    auto dtm = Delaunay{ std::begin(points), std::end(points) };
    auto to_surface_kernel = CGAL::Cartesian_converter<Point_kernel, Surface_kernel>{};
    triangulation_to_surface_mesh(dtm, mesh, to_surface_kernel);

}

template <typename Kernel>
struct URGMesh
{
    typedef Boundary_hull<Kernel> Hull;
    typedef CGAL::Epick Surface_kernel;
    typedef typename Surface_kernel::Point_3 Surface_point;
    typedef CGAL::Surface_mesh<Surface_point> Surface_mesh;
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

        typedef Surface_kernel Tk; // Triangulation_kernel
        //typedef CGAL::Triangulation_vertex_base_2<Tk> Vb;
        //typedef CGAL::Delaunay_mesh_face_base_2<Tk> Fb;
        //typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
        typedef CGAL::Exact_predicates_tag Itag;
        //typedef CGAL::Constrained_Delaunay_triangulation_2<Tk, Tds, Itag> CDT;
        typedef CGAL::Constrained_Delaunay_triangulation_2<Tk, CGAL::Default, Itag> CDT;
        //typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Mesh_criteria;
        //typedef CGAL::Delaunay_mesher_2<CDT, Mesh_criteria> Mesher;
        auto to_Tk = CGAL::Cartesian_converter<Kernel, Tk>{};

        surface_meshes.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            const auto& Pi = plane[i];
            typedef typename Hull::Surface_mesh Hull_surface;
            auto hull_mesh = Hull_surface{};
            hull.intersect(Pi, hull_mesh);
            auto cdt = CDT{};
            for (auto&& e : hull_mesh.edges()) {
                if(hull_mesh.is_border(e)) {
                    auto border_edge = domain[i].clip(
                        Segment{ 
                        hull_mesh.point(hull_mesh.vertex(e, 0)), 
                        hull_mesh.point(hull_mesh.vertex(e, 1)) 
                    });
                    if (border_edge) {
                        //std::cerr << "clipped border edge " << i << std::endl;
                        cdt.insert_constraint(
                            to_Tk(Pi.to_2d(border_edge->source())),
                            to_Tk(Pi.to_2d(border_edge->target()))
                        );
                    }
                }
            }
            for (auto&& Sh : segment_handles[i]) {
                const auto& S = segments[Sh.handle];
                cdt.insert_constraint(
                    to_Tk(Pi.to_2d(S.source())),
                    to_Tk(Pi.to_2d(S.target()))
                );
            }
            //auto mesher = Mesher{ cdt };
            //mesher.set_criteria(Mesh_criteria{ 0.125, 500. });
            //mesher.refine_mesh();
            auto& mesh = surface_meshes[i];
            add_edge_constraints_map(mesh);
            auto aPi = to_Tk(Pi);
            triangulation_to_surface_mesh_with_constraints(cdt, mesh, [&aPi](const auto& P) { return aPi.to_3d(P); });
        }
    }

    void isotropic_remeshing(const double& target_edge_length) {
        namespace PMP = CGAL::Polygon_mesh_processing;
        for (auto&& mesh : surface_meshes) {
            //typedef typename Surface_mesh::Edge_index Edge_index;
            auto constraints = get_edge_constraints_map(mesh);
            //auto result = mesh.add_property_map<Edge_index, bool>("e:tmp", true);
            //assert(result.second);
            //auto tmp = result.first;
            //std::vector<Edge_index> constrained_edges;
            //for (auto&& e : mesh.edges()) {
            //    if (constraints[e] || mesh.is_border(e)) {
            //        constrained_edges.emplace_back(e);
            //    }
            //    else {
            //        tmp[e] = false;
            //    }
            //}
            //CGAL::Polygon_mesh_processing::split_long_edges(
            //    constrained_edges, target_edge_length, mesh
            //);
            //for (auto&& e : mesh.edges()) {
            //    if (tmp[e]) {
            //        constraints[e] = true;
            //    }
            //}
            PMP::isotropic_remeshing(
                mesh.faces(), target_edge_length, mesh,
                PMP::parameters::edge_is_constrained_map(constraints)
                //.protect_constraints(true)
                //.relax_constraints(false)
            );
            //mesh.remove_property_map(tmp);
        }
    }

    void add_surface(const Surface_mesh& new_mesh) {
        auto mesh = Surface_mesh{ new_mesh };
        auto emap_new = get_edge_constraints_map(mesh);
        for (auto&& S: surface_meshes) {
            auto emap = get_edge_constraints_map(S);
            namespace parameters = CGAL::Polygon_mesh_processing::parameters;
            CGAL::Polygon_mesh_processing::corefine(
                S, mesh,
                parameters::edge_is_constrained_map(emap),
                parameters::edge_is_constrained_map(emap_new)
            );
        }
        surface_meshes.emplace_back(mesh);
    }

};

typedef CGAL::Epick Kernel;
typedef typename Kernel::Point_3 Point;
typedef typename Kernel::Vector_3 Vector;
typedef typename Kernel::Plane_3 Plane;
typedef URGMesh<CGAL::Epeck> Mesh;
typedef typename Mesh::Hull Hull;
typedef typename Hull::Surface_mesh Hull_surface;
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

    py::class_<Hull_surface>(module, "Hull_surface")
        .def("as_arrays", (decltype(&as_numpy_arrays<Hull_surface>))&as_numpy_arrays<Hull_surface>) // decltype is due to gcc bug
        ;

    py::class_<Surface_mesh>(module, "Surface_mesh")
        .def("as_arrays", (decltype(&as_numpy_arrays<Surface_mesh>))&as_numpy_arrays<Surface_mesh>) // decltype is due to gcc bug
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
        .def("isotropic_remeshing", &Mesh::isotropic_remeshing)
        .def("add_surface", &Mesh::add_surface)
        ;

    module.def("elevation_surface", [](py::array_t<double, py::array::c_style>& a) {
        auto mesh = std::make_unique<Surface_mesh>();
        add_edge_constraints_map(*mesh);
        elevation_surface(Array_view<const Point>{a}, *mesh);
        return mesh;
    });

}
