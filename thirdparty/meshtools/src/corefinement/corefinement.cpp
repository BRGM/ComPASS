#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include "UIdFactory.h"

typedef CGAL::Epick Kernel;
typedef Kernel::Point_3 Point;

typedef std::size_t Id_type;

template <typename Point_type>
using Surface_type = CGAL::Surface_mesh<Point_type>;

template <typename Point_type>
struct Point_id
{
    typedef Surface_type<Point_type> Surface;
    Surface * surface;
    Id_type on_surface_id;
};

template <typename Point_type>
using Curve_node = std::pair<Point_id<Point_type>, Point_id<Point_type>>;

template <typename Point_type>
struct Curve : std::list<Curve_node<Point_type>>
{};

template <typename Point_type>
struct On_curve_constraint
{
    Curve<Point_type> * curve;
    typename Curve<Point_type>::iterator position;
};

struct Point_on_surface : Point
{
    Point_id<Point> id;
    On_curve_constraint<Point> constraints;
    Point_on_surface() :
        Point{}, id{}, constraints{}
    {}
    Point_on_surface(double x, double y, double z) :
        Point{ x, y, z },
        id{}, constraints{}
    {}
    Point_on_surface(const Point& P) :
        Point{ P },
        id{}, constraints{}
    {}
    Point_on_surface(Point&& P) :
        Point{ std::forward<Point>(P) },
        id{}, constraints{}
    {}
};


//Surface_point id
//(surface, id_on_surface) = spi
//constraints 

//on_curve -> curve_node -> list<pair<spi, spi>)

// Curve : pair<Suface,Surface>, list<pair<pid,pid>> on_curve_constraint pair<curve, list_iterator>

//on_corner_constraint -> 3 curves


//corner -> ((S1, p1), (S2, p2), (S3, p3)) (Sk can be boundary of Sl)
// over_constrained corner 

//curve_patch = 1 or 2 points if only one point is corner (you don't want to move it)
//connected_point_id, (other_surface, id_on_other_surface)
//
//on deletion -> remove from connected points / test no longer connected

// split / check connected points arouns


// CHECKME: ? Encapsulate id_factory so that we have id{} instead of id{ id_factory.make_uid() }
//Id{
//    static Factory * factory; // = 0
//}

//struct Point_with_index: Point
//{
//    typedef UIdFactory<Id_type> Id_factory;
//    typedef typename Id_factory::Id_type Id;
//    static  Id_factory id_factory;
//    //typedef int Degree;
//    //typedef int Surface_index;
//    // CHECKME: ? Encapsulate id_factory so that we have id{} instead of id{ id_factory.make_uid() }
//    Id id;
//    //struct Constraints
//    //{
//    //    Degree degree; // 0 initialisation
//    //    std::array<Surface_index, 3> surfaces;
//    //    auto begin() { return surfaces.data(); }
//    //    auto end() { return surfaces.data() + degree; }
//    //    bool has_constraint(Surface_index si) {
//    //        return std::find(begin(), end(), si) != end();
//    //    }
//    //    void add_constraint(Surface_index si) {
//    //        assert(!has_constraint(si));
//    //        assert(degree < surfaces.size());
//    //        surfaces[degree] = si;
//    //        ++degree;
//    //        if (degree > 1) {
//    //            std::sort(begin(), end());
//    //        }
//    //    }
//    //};
//    //Constraints constraints;
//    Point_with_index() :
//        Point{},
//        id{ id_factory.make() }
//    {}
//    Point_with_index(double x, double y, double z) :
//        Point{ x, y, z },
//        id{ id_factory.make() }
//    {}
//    Point_with_index(const Point& P) :
//        Point{ P },
//        id{ id_factory.make() }
//    {}
//    Point_with_index(Point&& P) :
//        Point{ std::forward<Point>(P) },
//        id{ id_factory.make() }
//    {}
//};
//
//Point_with_index::Id_factory Point_with_index::id_factory;

typedef Surface_type<Point_on_surface> Mesh;

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

template <typename Mesh, typename Constraint_map>
auto collect_curves(const Mesh& mesh, const Constraint_map& constraints)
{

    typedef typename Mesh::Vertex_index Vertex_index;
    typedef typename Mesh::Edge_index Edge_index;

    typedef std::list<Vertex_index> Curve;
    std::vector<Curve> curves;

    auto collect_constrained_edges = [](const Mesh& mesh, const auto& constraints) {
        std::vector<Edge_index> result;
        for (auto&& edge : mesh.edges()) {
            if (constraints[edge]) {
                //std::cout << "ce:" << static_cast<Point>(mesh.point(mesh.vertex(edge, 0)))
                //    << "->" << static_cast<Point>(mesh.point(mesh.vertex(edge, 1))) << std::endl;
                result.emplace_back(edge);
            }
        }
        return result;
        //std::vector<std::pair<Vertex_index, Vertex_index>> result;
        //for (auto&& edge : mesh.edges()) {
        //    if (constraints[edge]) {
        //        result.emplace_back(
        //            mesh.vertex(edge, 0),
        //            mesh.vertex(edge, 1)
        //        );
        //    }
        //}
        //return result;
    };

    auto constrained_edges = collect_constrained_edges(mesh, constraints);

    auto edge_in_curve = std::map<Edge_index, bool>{};
    for (auto&& edge : constrained_edges) {
        edge_in_curve.emplace(edge, false);
    }
    auto find_free_edge = [&edge_in_curve]() {
        typedef boost::optional<Edge_index> Result;
        auto free_edge = std::find_if(
            begin(edge_in_curve), end(edge_in_curve),
            // find the first pair <key, value> with value = false, i.e. edge is not in curve
            [](auto&& pair) { return !pair.second; }
        );
        if (free_edge == end(edge_in_curve)) return Result{};
        return Result{ free_edge->first };
    };
    auto collect_vertices_until_corner = [&mesh, &edge_in_curve](Vertex_index start, auto out) {
        auto find_single_exit = [&edge_in_curve, &mesh](Vertex_index v) {
            std::cerr << "looking for exit at " << mesh.point(v) << " "
                      << std::count_if(begin(edge_in_curve), end(edge_in_curve), [](auto p) { return p.second; })
                      << " edges in curves" << std::endl;
            typedef std::pair<Edge_index, Vertex_index> Exit;
            typedef boost::optional<Exit> Result;
            auto result = Result{};
            for (auto&& h : CGAL::halfedges_around_source(v, mesh)) {
                auto p = edge_in_curve.find(mesh.edge(h));
                if (p != end(edge_in_curve)) {
                    if (!p->second) { // edge is not on curve
                        if (result) { // already one exit found
                            return Result{};
                        }
                        else {
                            result = Exit{ p->first, mesh.target(h) };
                        }
                    }
                }
            }
            return result;
        };
        auto next_vertex_found = find_single_exit(start);
        while (next_vertex_found) {
            edge_in_curve[next_vertex_found->first] = true;
            start = next_vertex_found->second;
            *out = start;
            ++out;
            next_vertex_found = find_single_exit(start);
        }
    };
    for (
        auto edge_left = find_free_edge();
        edge_left;
        edge_left = find_free_edge()
        ) {
        auto starting_edge = *edge_left;
        std::cerr << "#### new curve" << std::endl;
        curves.emplace_back();
        auto& curve = curves.back();
        edge_in_curve[starting_edge] = true;
        curve.emplace_back(mesh.vertex(starting_edge, 0));
        curve.emplace_back(mesh.vertex(starting_edge, 1));
        std::cerr << "#### new curve backward" << std::endl;
        collect_vertices_until_corner(curve.back(), std::back_inserter(curve));
        std::cerr << "#### new curve forward" << std::endl;
        collect_vertices_until_corner(curve.front(), std::front_inserter(curve));
    }
    assert(std::all_of(begin(edge_in_curve), end(edge_in_curve), [](auto&& p) {return p.second; }));
    return curves;

}


template <typename Mesh, typename Curves>
auto associate_curves(const Mesh& mesh1, const Curves& curves1, const Mesh& mesh2, const Curves& curves2)
{

    typedef typename Curves::value_type Curve;
    typedef typename Curve::value_type Vertex_index;
    static_assert(std::is_same<typename Mesh::Vertex_index, Vertex_index>::value, "Inconsistent types");

    assert(curves1.size() == curves2.size());
    std::vector<const Curve *> already_associated;
    already_associated.reserve(curves2.size());
    std::map<Vertex_index, Vertex_index> v1tov2, v2tov1;

    for (auto& curve1 : curves1) {
        bool found_corresponding_curve = false;
        for (auto& curve2 : curves2) {
            auto p = &curve2; 
            if (std::find(begin(already_associated), end(already_associated), p) == end(already_associated)) {
                if (curve1.size() == curve2.size()) {
                    bool all_points_associated = true;
                    for (auto&& v1 : curve1) {
                        if (v1tov2.count(v1) == 0) {
                            auto P1 = Point{ mesh1.point(v1) };
                            auto pv2 = std::find_if(begin(curve2), end(curve2), [&mesh2, &P1](const auto& v2) { return P1 == mesh2.point(v2); });
                            if (pv2 == end(curve2)) {
                                all_points_associated = false;
                                break;
                            }
                            else {
                                v1tov2[v1] = *pv2;
                            }
                        }
                    }
                    if (all_points_associated) {
                        found_corresponding_curve = true;
                        already_associated.emplace_back(&curve2);
                        break;
                    }
                }
            }
        }
        assert(found_corresponding_curve);
    }

    for (auto&& pair : v1tov2) {
        v2tov1[pair.second] = pair.first;
    }

#ifndef NDEBUG
    auto check_all_points_have_correspondances = [](const auto& curves, const auto& vmap) {
        for (auto&& curve : curves) {
            for (auto&& v : curve) {
                if (vmap.count(v) == 0) return false;
            }
        }
        return true;
    };
    assert(check_all_points_have_correspondances(curves1, v1tov2));
    assert(check_all_points_have_correspondances(curves2, v2tov1));
#endif

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
    w = tm2.add_vertex(Point(0.5, 0.6, 1));
    Mesh::Vertex_index x = tm2.add_vertex(Point(-0.5, 0.8, 1));
    tm2.add_face(u, v, w);
    tm2.add_face(u, w, x); // beware of face orientation - otherwise face is not added

    typedef Mesh::Vertex_index Vertex_index;
    typedef Mesh::Edge_index Edge_index;
    auto add_constraint_map = [](Mesh& mesh) {
        auto result = mesh.add_property_map<Edge_index, bool>("e:constrained", false);
        assert(result.second);
        return result.first;
    };
    auto constraints1 = add_constraint_map(tm1);
    auto constraints2 = add_constraint_map(tm2);
    namespace parameters = CGAL::Polygon_mesh_processing::parameters;
    CGAL::Polygon_mesh_processing::corefine(
        tm1, tm2,
        parameters::edge_is_constrained_map(constraints1),
        parameters::edge_is_constrained_map(constraints2),
        true
    );


    auto curves1 = collect_curves(tm1, constraints1);
    auto curves2 = collect_curves(tm2, constraints2);
    std::cout << "curve 1 ------------------------" << std::endl;
    for (auto&& curve : curves1) {
        std::cout << "-- curve" << std::endl;
        for (auto&& v : curve) {
            std::cout << "  " << tm1.point(v) << std::endl;
        }
    }
    std::cout << "curve 2 ------------------------" << std::endl;
    for (auto&& curve : curves2) {
        std::cout << "-- curve" << std::endl;
        for (auto&& v : curve) {
            std::cout << "  " << tm2.point(v) << std::endl;
        }
    }
    
    associate_curves(tm1, curves1, tm2, curves2);

    //std::ofstream os("test.off");
    //CGAL::write_off(os, tm1);
    //os.close();
    return py::make_tuple(mesh_as_arrays(tm1), mesh_as_arrays(tm2));

}

PYBIND11_MODULE(Corefinement, module)
{

    module.doc() = "pybind11 with CGAL surface corefinement (quick and dirty!!!)";

    module.def("test", &test);

}
