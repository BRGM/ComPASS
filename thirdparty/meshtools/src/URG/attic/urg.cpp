#include <array>
#include <vector>

#include <CGAL/Kernel_traits.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Convex_hull_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/barycenter.h>
#include <CGAL/boost/graph/Euler_operations.h>

//FIXME: prefer standard compliant optional
#include <boost/optional.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "common.h"
#include "conveniencies.h"

template <typename Point, typename Index>
struct Indexed_point
{
    Point point;
    Index index;
};

template <typename Segment, typename Index>
struct Indexed_segment
{
    Segment segment;
    Index index;
};

template <typename Surface_index>
struct Surface_vertex_index
{
    // There would be much less lines of code direcly replacing the variant by a std::set<Surface_index>
    //struct Not_indexed {};
    //typedef boost::variant <
    //    Not_indexed,
    //    Surface_index,
    //    std::array<Surface_index, 2>,
    //    std::array<Surface_index, 3>,
    //    std::set<Surface_index>
    //> Index_type;
    //Index_type index;
    //struct Adder : boost::static_visitor<Index_type> {
    //    Surface_index added;
    //    Adder(Surface_index si) :
    //        boost::static_visitor<Index_type>{},
    //        added{ si }
    //    {}
    //    result_type operator()(const Not_indexed&) const {
    //        return added;
    //    }
    //    result_type operator()(const Surface_index& si) const {
    //        if (added == si) return si;
    //        if (si < added) return std::array<Surface_index, 2>{ {si, added } };
    //        return std::array<Surface_index, 2>{ {added, si } };
    //    }
    //    result_type operator()(const std::array<Surface_index, 2>& sia) const {
    //        assert(sia[0] < sia[1]);
    //        if (added == sia[0] || added == sia[1]) return sia;
    //        if (added < sia[0]) return std::array<Surface_index, 3>{ { added, sia[0], sia[1] } };
    //        if (added < sia[1]) return std::array<Surface_index, 3>{ { sia[0], added, sia[1] } };
    //        return std::array<Surface_index, 3>{ { sia[0], sia[1], added } };
    //    }
    //    result_type operator()(const std::array<Surface_index, 3>& sia) const {
    //        if (added == sia[0] || added == sia[1] || added == sia[2]) return sia;
    //        auto s = std::set<Surface_index>{ begin(sia), end(sia) };
    //        s.insert(added);
    //        return s;
    //    }
    //    result_type operator()(std::set<Surface_index>& s) const {
    //        s.insert(added);
    //        return s;
    //    }
    //};
    //struct Degree : boost::static_visitor<std::size_t> {
    //    constexpr result_type operator()(const Not_indexed&) const noexcept {
    //        return 0;
    //    }
    //    constexpr result_type operator()(const Surface_index&) const noexcept {
    //        return 1;
    //    }
    //    constexpr result_type operator()(const std::array<Surface_index, 2>&) const noexcept {
    //        return 2;
    //    }
    //    constexpr result_type operator()(const std::array<Surface_index, 3>&) const noexcept {
    //        return 3;
    //    }
    //    result_type operator()(std::set<Surface_index>& s) const noexcept {
    //        return s.size();
    //    }
    //};
    //Surface_vertex_index& add(Surface_index si) {
    //    index = boost::apply_visitor(Adder{ si }, index);
    //    return *this;
    //}
    //std::size_t degree() const {
    //    return boost::apply_visitor(Degree{}, index);
    //}
    typedef std::set<Surface_index> Index_type;
    Index_type index;
    Surface_vertex_index() = default;
    Surface_vertex_index(const Surface_vertex_index&) = default;
    Surface_vertex_index& operator=(const Surface_vertex_index&) = default;
    Surface_vertex_index(Surface_index si) :
        index{ si }
    {}
    Surface_vertex_index& add(Surface_index si) {
        index.insert(si);
        return *this;
    }
    Surface_vertex_index& remove(Surface_index si) {
        index.erase(si);
        return *this;
    }
    auto degree() const {
        return index.size();
    }
};

template <typename Surface_index>
bool operator==(const Surface_vertex_index<Surface_index>& svi1, const Surface_vertex_index<Surface_index>& svi2)
{
    const auto& s1 = svi1.index;
    const auto& s2 = svi2.index;
    if (s1.size() != s2.size()) return false;
    std::vector<Surface_index> intersection;
    intersection.reserve(s1.size());
    std::set_intersection(begin(s1), end(s1), begin(s2), end(s2), std::back_inserter(intersection));
    return intersection.size() == s1.size();
}

//template <typename Surface_index, typename Surface_mesh>
//auto add_vertex_degree_map(Surface_mesh& mesh)
//{
//    typedef typename Surface_mesh::Vertex_index Vertex_index;
//    typedef Surface_vertex_index<Surface_index> Vertex_degree;
//    auto pmap = mesh.add_property_map<Vertex_index, Vertex_degree>("vertex_degree", Vertex_degree{});
//    assert(pmap.second);
//    return pmap.first;
//}
//
//template <typename Surface_mesh>
//auto add_is_contrained_edge_map(Surface_mesh& mesh)
//{
//    typedef typename Surface_mesh::Edge_index Edge_index;
//    auto pmap = mesh.add_property_map<Edge_index, bool>("is_constrained_edge", false);
//    assert(pmap.second);
//    return pmap.first;
//}

//template <typename Surface_index, typename Surface_mesh>
//void add_properties(Surface_mesh& mesh)
//{
//    add_vertex_degree_map<Surface_index>(mesh);
//    add_is_contrained_edge_map(mesh);
//}

template <typename Point>
using MyMesh = CGAL::Surface_mesh<Point>;

struct Split_recorder
{
    typedef MyMesh<CGAL::Epeck::Point_3> Mesh_type;
    typedef typename Mesh_type::Halfedge_index Halfedge_index;
    typedef std::vector<std::pair<Halfedge_index, Halfedge_index>> Records;
    typedef std::map<Mesh_type*, Records> Records_map;
    static Records_map records_map;
    static void reset(Mesh_type& M1, Mesh_type& M2)
    {
        records_map.clear();
        records_map.insert({ &M1, Records{} });
        records_map.insert({ &M2, Records{} });
    }
    static void record_split(Mesh_type& mesh, Halfedge_index hold, Halfedge_index hnew)
    { 
        assert(mesh.next(hnew) == hold);
        assert(records_map.count(&mesh) == 1);
        records_map[&mesh].emplace_back(hold, hnew);
    }
    auto associate_records()
    {
        assert(records_map.size() == 2);

    }
};



Split_recorder::Records_map Split_recorder::records_map;

namespace CGAL {
    namespace Euler {
        template <typename Point>
        typename boost::graph_traits<MyMesh<Point>>::halfedge_descriptor
            split_edge(typename boost::graph_traits<MyMesh<Point>>::halfedge_descriptor h, MyMesh<Point>& g)
        {
            auto hnew = opposite(CGAL::Euler::split_vertex(prev(h, g), opposite(h, g), g), g);
            Split_recorder::record_split(g, h, hnew);
            return hnew;
        }
    }
}

template <typename Point_type, typename Surface_index_type>
struct Surface_factory
{
    typedef Point_type Point;
    //typedef CGAL::Surface_mesh<Point> Surface_mesh;
    typedef MyMesh<Point> Surface_mesh;
    typedef Surface_index_type Surface_index;
    struct Surface
    {
        typedef Surface_mesh Surface_mesh;
        typedef Surface_index_type Surface_index;
        struct Border_constraint {
            typedef Surface_index value_type;
            value_type id;
            constexpr bool operator<(const Border_constraint& other) const noexcept { return id < other.id; }
            constexpr bool operator==(const Border_constraint& other) const noexcept { return id == other.id; }
            static constexpr value_type default_value() noexcept { return std::numeric_limits<Surface_index>::max(); }
            Border_constraint() :
                id{ default_value() }
            {}
            Border_constraint(const Border_constraint&) = default;
            Border_constraint& operator=(const Border_constraint&) = default;
            Border_constraint(value_type s) :
                id{ s }
            {}
            bool unspecified() const noexcept { return id == default_value(); }
        };
        struct Surface_constraint {
            typedef Surface * value_type;
            //typedef Surface_index value_type;
            value_type source;
            constexpr bool operator<(const Surface_constraint& other) const noexcept { return source < other.source; }
            constexpr bool operator==(const Surface_constraint& other) const noexcept { return source == other.source; }
            //static constexpr value_type default_value() noexcept { return std::numeric_limits<Surface_index>::max(); }
            static constexpr value_type default_value() noexcept { return nullptr; }
            Surface_constraint() :
                source{ default_value() }
            {}
            Surface_constraint(const Surface_constraint&) = default;
            Surface_constraint& operator=(const Surface_constraint&) = default;
            Surface_constraint(value_type s) :
                source{ s }
            {}
            bool unspecified() const noexcept { return source == default_value(); }
        };
        typedef boost::variant<Border_constraint, Surface_constraint> Constraint;
        struct Unspecified : boost::static_visitor<bool> {
            template <typename T>
            result_type operator()(const T& t) const {
                return t.unspecified();
            }
        };
        static bool unspecified(const Constraint& variant) {
            return boost::apply_visitor(Unspecified{}, variant);
        }
        struct Human_readable : boost::static_visitor<std::string> {
            result_type operator()(const Border_constraint& BC) const {
                assert(!Surface::unspecified(BC));
                return std::string("B") + std::to_string(BC.id);
            }
            result_type operator()(const Surface_constraint& SC) const {
                assert(!Surface::unspecified(SC));
                return std::string("S") + std::to_string(SC.source->index);
            }
        };
        static std::string human_readable(const Constraint& variant) {
            return boost::apply_visitor(Human_readable{}, variant);
        }
        typedef typename Surface_mesh::Vertex_index Vertex_index;
        typedef typename Surface_mesh::Edge_index Edge_index;
        //typedef Surface_vertex_index<Surface_index> Vertex_info;
        //constexpr static const auto vertex_info_map_name = "v:info";
        constexpr static const auto corners_map_name = "v:is_corner";
        constexpr static const auto edge_constraints_map_name = "e:is_constrained";
        constexpr static const auto edge_constraint_source_name = "e:constraint_source";
        Surface_mesh mesh;
        Surface_index index;
        friend struct Surface_factory<Surface_mesh, Surface_index>;
        typedef std::set<Edge_index> Constrained_edges;
        std::map<Constraint, Constrained_edges> external_constraints_map;
        std::map<Vertex_index, std::set<Constraint>> corner_constraints;
        Surface() = delete;
        template <typename Location, typename value_type>
        void create_named_property_map(const char* name, value_type default_value = value_type{}) {
            auto pmap = mesh.add_property_map<
                Location, value_type
            >(name, default_value);
            assert(pmap.second);
        }
        template <typename Location, typename value_type>
        auto retrieve_named_property_map(const char* name) {
            auto pmap = mesh.property_map<
                Location, value_type
            >(name);
            assert(pmap.second);
            return pmap.first;
        }
        Surface(Surface_index si) :
            mesh{},
            index{ si },
            external_constraints_map{}
        {
            //create_named_property_map<Vertex_index, Vertex_info>(vertex_info_map_name, { index });
            create_named_property_map<Vertex_index, bool>(corners_map_name, false);
            create_named_property_map<Edge_index, bool>(edge_constraints_map_name, false);
            create_named_property_map<Edge_index, Constraint>(edge_constraint_source_name);
        }
        //auto vertex_info_map() {
        //    return retrieve_named_property_map<Vertex_index, Vertex_info>(vertex_info_map_name);
        //}
        auto corners_map() {
            return retrieve_named_property_map<Vertex_index, bool>(corners_map_name);
        }
        auto edge_constraints_map() {
            return retrieve_named_property_map<Edge_index, bool>(edge_constraints_map_name);
        }
        auto edge_constraint_source_map() {
            return retrieve_named_property_map<Edge_index, Constraint>(edge_constraint_source_name);
        }
        template <typename Location, typename value_type>
        auto property_map(value_type t = value_type{}) {
            auto pair = mesh.add_property_map<Location, value_type>("", t);
            assert(pair.second);
            return pair.first;
        }
        template <typename Property_map>
        void remove_property_map(Property_map& p) {
            mesh.remove_property_map(p);
        }
        Constrained_edges& external_constraints(Constraint c) {
            assert(!Surface::unspecified(c));
            auto p = external_constraints_map.find(c);
            if (p == external_constraints_map.end()) {
                auto pair = external_constraints_map.insert({ c, Constrained_edges{} });
                assert(pair.second);
                p = pair.first;
            }
            return p->second;
        }
        void register_external_constraint(Constraint c, Edge_index e) {
            external_constraints(c).insert(e);
        }
        void register_external_constraints(Constraint c, Constrained_edges&& edges) {
            external_constraints(c).merge(std::forward<Constrained_edges>(edges));
        }
        bool is_constrained_by(Surface& other) const noexcept {
            return external_constraints_map.count(&other) > 0;
        }
        std::size_t number_of_constraints(Surface& other) const noexcept {
            auto p = external_constraints_map.find(&other);
            if (p == external_constraints_map.end()) {
                return 0;
            }
            return p->second.size();
        }
        void mark_all_vertices_as_corners() {
            auto is_corner = corners_map();
            for (auto&& v : mesh.vertices()) {
                is_corner[v] = true;
            }
        }
        void register_corner(Vertex_index v, Constraint C1, Constraint C2) {
            assert(!mesh.is_removed(v));
            auto is_corner = corners_map();
            assert(is_corner[v] || corner_constraints.count(v) == 0);
            if (!is_corner[v]) {
                is_corner[v] = true;
                corner_constraints.emplace(v, {});
            }
            auto& s = corner_constraints[v];
            s.insert(Constraint{ this });
            s.insert(C1);
            s.insert(C2);
        }
        void relax_corner(Vertex_index v, Constraint C) {
            assert(!mesh.is_removed(v));
            assert(!C.source == this);
            auto is_corner = corners_map();
            assert(is_corner[v] && corner_constraints.count(v) == 1 && corner_constraints[v].count(C) == 1);
            auto& s = corner_constraints[v];
            s.erase(C);
            if (s.size() < 3) {
                corner_constraints.erase(C);
                is_corner[v] = false;
            }
        }
        static void register_corner(
            Vertex_index v1, Constraint C1,
             Vertex_index v2, Constraint C2,
             Vertex_index v3, Constraint C3
        ) {
            C1.source->register_corner(v1, C2, C3);
            C2.source->register_corner(v2, C1, C3);
            C3.source->register_corner(v3, C1, C2);
        }
        //struct Constraint_remover
        //{

        //    void operator()(Border_contraint&) {
        //        assert(false);
        //    }
        //    void operator()(Surface_contraint& SC) {
        //        assert(false);
        //    }
        //};

        //static void remove_constraint(Constraint& c)
        //{
        //    typedef typename Surface::Constraint Constraint;
        //    auto is_constrained = surface.edge_constraints_map();
        //    auto constraint_source = surface.edge_constraint_source_map();
        //    auto& constraints = surface.external_constraints(c);
        //    for (auto&& e : constraints) {
        //        is_constrained[e] = false;
        //        constraint_source[e] = Constraint{};
        //    }
        //    surface.external_constraints_map.erase(c);
        //}
    };
    Surface_index new_index;
    Surface_factory() :
        new_index{ 0 }
    {}
    Surface_factory(const Surface_factory&) = delete;
    Surface_factory& operator=(const Surface_factory&) = delete;
    auto make() {
        return Surface{ new_index++ };
    }
    auto make_unique_pointer() {
        return std::make_unique<Surface>( new_index++ );
    }
};

template <typename Kernel_type>
struct Boundary_hull
{

    typedef Kernel_type Kernel;
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Segment_3 Segment;
    typedef typename Kernel::Plane_3 Plane;
    typedef std::vector<Plane> Plane_collection;
    typedef typename Plane_collection::size_type Plane_index;
    typedef std::array<Plane_index, 3> Corner_index;
    typedef Indexed_point<Point, Corner_index> Corner;
    typedef std::array<Plane_index, 2> Edge_index;
    typedef Indexed_segment<Segment, Edge_index> Edge;
    typedef std::vector<Edge> Edge_collection;
    typedef CGAL::Surface_mesh<Point> Surface_mesh;
    typedef std::vector<Surface_mesh> Mesh_collection;

    Plane_collection planes;
    Edge_collection edges;
    std::vector<std::vector<Corner>> corners;
    Mesh_collection boundaries;

    inline bool is_inside(const Point& P) const noexcept {
        return std::all_of(
            begin(planes), end(planes),
            [P](const auto& plane) {
                return !plane.has_on_positive_side(P);
            }
        );
    }

    template <typename Plane_iterator>
    Boundary_hull(Plane_iterator first, Plane_iterator last):
        planes{},
        edges{},
        corners{},
        boundaries{}
    {
        typedef typename std::iterator_traits<Plane_iterator>::value_type Source_plane;
        typedef typename CGAL::Kernel_traits<Source_plane>::Kernel Source_kernel;
        typedef CGAL::Cartesian_converter<Source_kernel, Kernel> Kernel_converter;
        std::transform(first, last, std::back_inserter(planes), Kernel_converter{});
        const auto n = planes.size();
        std::vector<std::vector<Corner>> corners(n);
        std::map<Edge_index, std::vector<Corner>> edge_map;
        for (Plane_index i = 0; i < n; ++i) {
            for (Plane_index j = i + 1; j < n; ++j) {
                for (Plane_index k = j + 1; k < n; ++k) {
                    auto intersection = CGAL::intersection(planes[i], planes[j], planes[k]);
                    if (intersection) {
                        if (const Point * p = boost::get<Point>(&*intersection)) {
                            if (is_inside(*p)) {
                                const auto P = Corner{ *p,{ i, j, k } };
                                corners[i].push_back(P);
                                corners[j].push_back(P);
                                corners[k].push_back(P);
                                edge_map[{i, j}].push_back(P);
                                edge_map[{i, k}].push_back(P);
                                edge_map[{j, k}].push_back(P);
                            }
                        }
                    }
                }
            }
        }
        edges.reserve(edge_map.size());
        for (const auto& edge : edge_map) {
            const auto& index = edge.first;
            const auto& points = edge.second;
            std::cerr << "Edge (" << index[0] << "," << index[1] << ") with" << points.size() << " points" << std::endl;
            assert(points.size() == 2);
            // FIXME: could be emplace_back here but does not work on MSVC 2015
            edges.push_back(Edge{ Segment{ points.front().point, points.back().point }, index });
        }
        edges.shrink_to_fit();
        std::cerr << "Collected " << edges.size() << " edges" << std::endl;
        assert(boundaries.empty());
        boundaries.reserve(n);
        for (Plane_index i = 0; i < n; ++i) {
            const auto& plane_corners = corners[i];
            std::cerr << "Size corners" << plane_corners.size() << std::endl;
            std::vector<Point> points;
            points.reserve(plane_corners.size());
            for (const auto& corner : plane_corners) points.push_back(corner.point);
            boundaries.emplace_back();
            auto& boundary = boundaries.back();
            CGAL::convex_hull_3(cbegin(points), cend(points), boundary);
            std::cout << "The convex hull contains " << boundary.number_of_vertices() << " vertices and " << boundary.number_of_faces() << " faces" << std::endl;
            for (const auto& face : boundary.faces()) {
                std::cerr << "Face with " << boundary.degree(face) << " points" << std::endl;
            }
        }
        boundaries.shrink_to_fit();
    }

    void intersect(const Plane& plane, Surface_mesh& mesh) const
    {
        std::vector<Point> points;
        for (const auto& edge : edges) {
            auto intersection = CGAL::intersection(plane, edge.segment);
            if (intersection) {
                if (const Point * p = boost::get<Point>(&*intersection)) {
                    if (is_inside(*p)) {
                        points.push_back(*p);
                    }
                }
                else {
                    if (const Segment * S = boost::get<Segment>(&*intersection)) {
                        // here *S and edge.segment are the same
                        points.push_back(S->source());
                        points.push_back(S->target());
                    }
                }
            }
        }
        std::cerr << "Hull intersection with " << points.size() << " points." << std::endl;
        CGAL::convex_hull_3(cbegin(points), cend(points), mesh);
    }

};

namespace py = pybind11;

template <typename Surface>
auto constrained_midpoints(Surface& surface)
{

    typedef typename Surface::Surface_mesh Surface_mesh;
    typedef typename Surface_mesh::Point Surface_point;
    typedef typename CGAL::Kernel_traits<Surface_point>::Kernel Surface_kernel;
    typedef CGAL::Epick::Point_3 Point;
    auto to_epick = CGAL::Cartesian_converter<Surface_kernel, CGAL::Epick>{};

    auto& mesh = surface.mesh;
    auto is_constrained_edge = surface.edge_constraints_map();
    std::vector<Point> midpoints;
    for (auto&& e : mesh.edges()) {
        if (is_constrained_edge[e]) {
            midpoints.emplace_back(to_epick(CGAL::midpoint(mesh.point(mesh.vertex(e, 0)), mesh.point(mesh.vertex(e, 1)))));
        }
    }
    auto result = py::array_t<double, py::array::c_style>{
            { midpoints.size(), static_cast<std::size_t>(Point::Ambient_dimension::value) }
    };
    auto p = reinterpret_cast<Point*>(result.request().ptr);
    static_assert(sizeof(Point) == Point::Ambient_dimension::value * sizeof(double), "Inconsistent sizes in memory!");
    for (auto&& M : midpoints) {
        *(p++) = M;
    }
    return result;
}


template <typename Plane_type>
struct On_plane_side
{
    typedef Plane_type Plane;
    typedef typename CGAL::Kernel_traits<Plane>::Kernel Kernel;
    typedef typename Kernel::Point_3 Point;
protected:
    Plane plane;
    CGAL::Oriented_side side;
public:
    On_plane_side(const Plane& P, const Point& A) :
        plane{ P },
        side{ plane.oriented_side(A) }
    {}
    On_plane_side(Plane&& P, const Point& A) :
        plane{ std::forward<Plane>(P) },
        side{ plane.oriented_side(A) }
    {}
    bool operator()(const Point& P) const noexcept {
        return plane.oriented_side(P) == side;
    }
};

template <typename Surface_mesh>
inline auto face_barycenter(const typename Surface_mesh::Halfedge_index h, const Surface_mesh& mesh)
{

    typedef typename Surface_mesh::Point Surface_point;
    typedef typename CGAL::Kernel_traits<Surface_point>::Kernel Surface_kernel;
    // special type to be used with CGAL:barycenter
    typedef std::pair<Surface_point, typename Surface_kernel::FT> Face_point;

    assert(!mesh.is_border(h));
    std::vector<Face_point> face_points;
    for (auto&& v : CGAL::vertices_around_face(h, mesh)) {
        face_points.emplace_back(mesh.point(v), 1);
    }
    assert(!face_points.empty());
    return CGAL::barycenter(cbegin(face_points), cend(face_points));

}

//template <typename Surface>
//auto update_vertex_properties(Surface& surface)
//{
//    typedef typename Surface::Vertex_info Vertex_info;
//    auto& mesh = surface.mesh;
//    auto vertex_info = surface.vertex_info_map();
//    for (auto&& v : mesh.vertices()) {
//        if (mesh.is_removed(v)) {
//            auto& vi = vertex_info[v];
//            assert(vi.degree() == 0 || vi == Vertex_info{ surface.index });
//            if (!vi.degree() == 0) {
//                vi = Vertex_info{};
//            }
//        }
//    }
//}

//template <typename Surface>
//auto update_properties(Surface& surface)
//{
//    typedef typename Surface::Surface_index Surface_index;
//    auto& mesh = surface.mesh;
//    //auto vertex_info = surface.vertex_info_map();
//    auto is_constrained_edge = surface.edge_constraints_map();
//    auto edge_constraint_source = surface.edge_constraint_source_map();
//    std::set<Surface_index> relaxed_edge_constraints;
//    for (auto&& e : mesh.edges()) {
//        if (mesh.is_removed(e)) {
//            if (is_constrained_edge[e]) {
//                const auto constraint_source = edge_constraint_source[e];
//                relaxed_edge_constraints.insert(constraint_source);
//                //for (int i = 0; i < 2; ++i) {
//                //    vertex_info[mesh.vertex(e, i)].remove(constraint_source);
//                //}
//                is_constrained_edge[e] = false;
//                edge_constraint_source[e] = std::numeric_limits<Surface_index>::max();
//            }
//            assert(!is_constrained_edge[e]);
//            assert(edge_constraint_source[e] == std::numeric_limits<Surface_index>::max());
//        }
//    }
//    //update_vertex_properties(surface);
//    assert(
//        std::all_of(begin(relaxed_edge_constraints), end(relaxed_edge_constraints),
//            [&surface](auto si) {
//        const auto& edges = surface.external_constraints(si);
//        return std::all_of(begin(edges), end(edges),
//            [&surface](auto e) {
//            return surface.mesh.is_removed(e) || surface.mesh.is_border(e);
//        });
//    })
//    );
//    std::cerr << "Relaxed constraints on " << surface.index << " : ";
//    std::copy(begin(relaxed_edge_constraints), end(relaxed_edge_constraints), std::ostream_iterator<Surface_index>(std::cerr, " "));
//    std::cerr << std::endl;
//    return relaxed_edge_constraints;
//}

// TODO improve by locating connected component with point
template <typename Surface, typename Mask>
auto select_connected_components(Surface& surface, const Mask& is_kept)
{

    typedef typename Surface::Surface_index Surface_index;
    typedef typename Surface::Constraint Constraint;
    typedef typename Surface::Surface_mesh Surface_mesh;
    typedef typename Surface_mesh::Face_index Face_index;
    typedef typename Surface_mesh::faces_size_type Component_id;

    namespace parameters = CGAL::Polygon_mesh_processing::parameters;
    auto components = surface.property_map<Face_index, Component_id>();
    auto& mesh = surface.mesh;
    auto is_constrained_edge = surface.edge_constraints_map();
    auto constraint_source = surface.edge_constraint_source_map();
    CGAL::Polygon_mesh_processing::connected_components(
        mesh, components,
        parameters::edge_is_constrained_map(is_constrained_edge)
    );
    std::map<Component_id, bool> is_kept_component;
    for (auto&& f : mesh.faces()) {
        const auto C = components[f];
        if (is_kept_component.count(C) == 0) {
            is_kept_component[C] = is_kept(face_barycenter(mesh.halfedge(f), mesh));
        }
        assert(is_kept_component[C] == is_kept(face_barycenter(mesh.halfedge(f), mesh)));
    }
    std::set<Constraint> relaxed_edge_constraints;
    std::cerr << "Checking constraint: ";
    for (auto&& p : surface.external_constraints_map) {
        const auto c = p.first;
        std::cerr << Surface::human_readable(c) << " (" << p.second.size() << ") ";
        if (!p.second.empty()) {
            const auto e_test = *(p.second.begin());
            bool relax = true;
            for (int i = 0; i < 2; ++i) {
                const auto hi = mesh.halfedge(e_test, i);
                if (!mesh.is_border(hi)) {
                    relax = relax && !(is_kept_component[components[mesh.face(hi)]]);
                }
            }
            if (relax) {
                relaxed_edge_constraints.insert(c);
                for (auto&& e : p.second) {
                    assert(is_constrained_edge[e]);
                    assert(constraint_source[e]==c);
                    is_constrained_edge[e] = false;
                    constraint_source[e] = Constraint{};
                }
            }
            if (relax) std::cerr << "R ";
            else std::cerr << "U ";
        }
    }
    std::cerr << std::endl;
    //for (auto&& c : relaxed_edge_constraints) {
    //    assert(!Surface::unspecified(c));
    //    remove_constraint(c, Surface_constraint{ &surface });
    //    surface.external_constraints_map.erase(c);
    //}
    std::cerr << "Keeping on " << surface.index << ": ";
    for (auto&& p : is_kept_component) {
        if (p.second) {
            std::cerr << p.first << " ";
        }
    }
    std::cerr << "leaving: ";
    for (auto&& p : is_kept_component) {
        if (!p.second) {
            std::cerr << p.first << " ";
        }
    }
    std::cerr << std::endl;
    std::cerr << "Relaxed constraints on " << surface.index << " : ";
    for (auto&& c : relaxed_edge_constraints) std::cerr << Surface::human_readable(c) << " ";
    std::cerr << std::endl;
    std::vector<Component_id> kept_components;
    for (auto&& p : is_kept_component) {
        if (p.second) {
            kept_components.emplace_back(p.first);
        }
    }
    CGAL::Polygon_mesh_processing::keep_connected_components(mesh, kept_components, components);
    surface.remove_property_map(components);
    return relaxed_edge_constraints;
}


// TODO: use connected component ?
template <typename Surface, typename Mask>
inline auto select_faces(Surface& surface, const Mask& is_kept)
{
    typedef typename Surface::Vertex_index Vertex_index;
    //typedef typename Surface::Vertex_info Vertex_info;
    typedef typename Surface::Edge_index Edge_index;
    //std::set<Vertex_index> possibly_removed_vertices;
    //std::set<Edge_index> possibly_removed_edges;
    auto& mesh = surface.mesh;
    typedef typename Surface::Surface_mesh::Point Surface_point;
    typedef typename CGAL::Kernel_traits<Surface_point>::Kernel Surface_kernel;
    auto to_epick = CGAL::Cartesian_converter<Surface_kernel, CGAL::Epick>{};
    std::ofstream removedfaces{ "removedfaces.csv" };
    for (auto&& f : mesh.faces()) {
        auto h = mesh.halfedge(f);
        if (!is_kept(face_barycenter(h, mesh))) {
            removedfaces << to_epick(face_barycenter(h, mesh)) << std::endl;
            //for (auto&& v : CGAL::vertices_around_face(h, mesh)) {
            //    possibly_removed_vertices.insert(v);
            //}
            //for (auto&& e : CGAL::edges_around_face(h, mesh)) {
            //    possibly_removed_edges.insert(e);
            //}
            CGAL::Euler::remove_face(h, mesh);
        }
    }
    //auto info = surface.vertex_info_map();
    //// TODO: only edges must be collected
    ////for (auto&& v : possibly_removed_vertices) {
    ////    if (mesh.is_removed(v)) {
    ////        info[v] = Vertex_info{};
    ////    }
    ////}
    //auto is_constrained = surface.edge_constraints_map();
    //auto constraint_source = surface.edge_constraint_source_map();
    //std::set<Surface_index> relaxed_constraints;
    //for (auto&& e : possibly_removed_edges) {
    //    if (mesh.is_removed(e)) {
    //        //std::cerr << "removing " << e;
    //        //if (is_constrained[e])
    //        //    std::cerr << " with constraint " << constraint_source[e];
    //        //std::cerr << std::endl;
    //        if (is_constrained[e]) {
    //            is_constrained[e] = false;
    //            relaxed_constraints.insert(constraint_source[e]);
    //            for (int i = 0; i < 2; ++i) {
    //                const auto v = mesh.vertex(e, i);
    //                if (mesh.is_removed(v)) {
    //                    info[v] = Vertex_info{};
    //                }
    //        else {
    //                    info[v].remove(constraint_source[e]);
    //                }
    //            }
    //            constraint_source[e] = std::numeric_limits<Surface_index>::max();
    //        }
    //    }
    //}
    //for (auto&& si : relaxed_constraints) {
    //    bool good = true;
    //    for (auto&& e : surface.external_constraints(si)) {
    //        if (!(mesh.is_removed(e) || mesh.is_border(e))) {
    //            std::cerr << "!!!!!!!!!!!!!!!!!!!! " << e << " on " << surface.index << " is constrained by " << si << "(" << mesh.is_removed(e) << mesh.is_border(e) << ")" << std::endl;
    //            good = false;
    //        }
    //    }
    //    assert(good);
    //}
    //assert(
    //    std::all_of(begin(relaxed_constraints), end(relaxed_constraints),
    //        [&surface](auto si) {
    //    const auto& edges = surface.external_constraints(si);
    //     return std::all_of(begin(edges), end(edges),
    //         [&surface](auto e) {
    //         return surface.mesh.is_removed(e) || surface.mesh.is_border(e);
    //     });
    //        })
    //);
    //return update_properties(surface);
}

//typedef Boundary_hull<CGAL::Epeck> Hull;
//typedef typename Hull::Surface_mesh Surface_mesh;
//struct Mask
//{
//    template <typename T>
//    bool operator()(const T&) const { return true; }
//};

template <typename Surface_mesh>
void simplify_faces(Surface_mesh& mesh)
{

    auto h = *(begin(mesh.halfedges()));
    if (mesh.is_border(h)) h = mesh.opposite(h);
    typedef boost::optional<decltype(h)> Flag;
    auto loop_entry = Flag{};
    while (!loop_entry || h != *loop_entry) {
        if (mesh.is_border(mesh.opposite(h))) {
            if (!loop_entry) loop_entry = Flag{ h };
            h = mesh.next(h);
        }
        else {
            loop_entry = Flag{};
            h = CGAL::Euler::join_face(h, mesh);
        }
    }
    auto start = h;
    typedef typename Surface_mesh::Point Surface_point;
    typedef typename CGAL::Kernel_traits<Surface_point>::Kernel Surface_kernel;
    auto to_epick = CGAL::Cartesian_converter<Surface_kernel, CGAL::Epick>{};
    std::cerr << "Points (" << mesh.degree(mesh.face(h)) << "): ";
    for (h = mesh.next(h); h != start; h = mesh.next(h)) {
        std::cerr << to_epick(mesh.point(mesh.target(h))) << " | ";
    }
    std::cerr << std::endl;
    assert(!mesh.is_border(h) && mesh.is_border(mesh.opposite(h)));
    const auto P = face_barycenter(h, mesh);
    h = CGAL::Euler::add_center_vertex(h, mesh);
    mesh.point(mesh.target(h)) = P;

}

template <typename Surface>
void remove_center_vertices(Surface& surface)
{

    typedef typename Surface::Surface_mesh Surface_mesh;
    typedef typename Surface_mesh::Vertex_index Vertex_index;

    auto& mesh = surface.mesh;
    auto is_constrained_edge = surface.edge_constraints_map();
    auto edge_constraint_source = surface.edge_constraint_source_map();
    auto& external_constraints = surface.external_constraints_map;
    std::set<Vertex_index> center_vertices;
    for (auto&& v : mesh.vertices()) {
        if (center_vertices.find(v) == center_vertices.end()) {
            bool is_center = true;
            for (auto& h : CGAL::halfedges_around_source(v, mesh)) {
                const auto e = mesh.edge(h);
                if (mesh.is_border(e) || is_constrained_edge[e]) {
                    is_center = false;
                    break;
                }
            }
            if (is_center) {
                assert(
                    std::all_of(
                        begin(CGAL::halfedges_around_source(v, mesh)),
                        end(CGAL::halfedges_around_source(v, mesh)),
                        [&mesh, &is_constrained_edge, &edge_constraint_source](auto h) {
                    const auto e = mesh.edge(h);
                    return (!is_constrained_edge[e]) && Surface::unspecified(edge_constraint_source[e]);
                        }
                        )
                );
                center_vertices.insert(v);
            }
        }
    }
    for (auto&& v : center_vertices) {
        CGAL::Euler::remove_center_vertex(mesh.halfedge(v), mesh);
    }
    //update_vertex_properties(surface);

}


template <typename Surface>
void remove_inner_edges(Surface& surface)
{

    typedef typename Surface::Surface_mesh Surface_mesh;
    typedef typename Surface::Surface_index Surface_index;
    typedef typename Surface::Constraint Constraint;
    typedef typename Surface_mesh::Edge_index Edge_index;

    auto& mesh = surface.mesh;
    auto is_constrained_edge = surface.edge_constraints_map();
    auto edge_constraint_source = surface.edge_constraint_source_map();
    auto& external_constraints = surface.external_constraints_map;
    std::vector<Edge_index> inner_edges;
    for (auto&& e : mesh.edges()) {
        if (!(mesh.is_border(e) || is_constrained_edge[e])) {
            inner_edges.emplace_back(e);
        }
    }
    for (auto&& e : inner_edges) {
        if (is_constrained_edge[e]) {
            external_constraints[edge_constraint_source[e]].erase(e);
            is_constrained_edge[e] = false;
            edge_constraint_source[e] = Constraint{};
        }
        assert(!is_constrained_edge[e]);
        assert(Surface::unspecified(edge_constraint_source[e]));
        CGAL::Euler::join_face(mesh.halfedge(e), mesh);
    }

}

template <typename Surface>
void simplify_face_borders(Surface& surface)
{

    typedef typename Surface::Surface_mesh Surface_mesh;
    typedef typename Surface::Constraint Constraint;
    typedef typename Surface::Surface_index Surface_index;
    typedef typename Surface_mesh::Halfedge_index Halfedge_index;
    //typedef std::vector<Vertex_index> Chain;
    
    //std::vector<Chain> chains;
    auto& mesh = surface.mesh;
    //auto vertex_info = surface.vertex_info_map();
    auto is_corner = surface.corners_map();
    auto is_constrained_edge = surface.edge_constraints_map();
    auto edge_constraint_source = surface.edge_constraint_source_map();
    auto& external_constraints = surface.external_constraints_map;
    //chains.emplace_back();
    //std::vector<Halfedge_index> removed_vertices;
    std::multiset<Constraint> constraints_around_vertex;
    for (auto&& v : mesh.vertices()) {
        if (is_corner[v]) continue;
        int borders_around_v = 0;
        Halfedge_index removed;
        for (auto&& h : CGAL::halfedges_around_source(v, mesh)) {
            const auto e = mesh.edge(h);
            if (mesh.is_border(e)) {
                ++borders_around_v;
                removed = h;
            }
            if (is_constrained_edge[e]) {
                constraints_around_vertex.insert(edge_constraint_source[e]);
                removed = h;
            }
        }
        // assertion failure means a corner has been missed or a center vertex has not been removed 
        if (
            !(borders_around_v == 2 && constraints_around_vertex.empty()) &&
            !(borders_around_v == 2 && constraints_around_vertex.size() == 2 && constraints_around_vertex.count(*(constraints_around_vertex.begin())) == 2) &&
            !(borders_around_v == 0 && constraints_around_vertex.size() == 2 && constraints_around_vertex.count(*(constraints_around_vertex.begin())) == 2)
            ) {
            typedef typename Surface::Surface_mesh Surface_mesh;
            typedef typename Surface_mesh::Point Surface_point;
            typedef typename CGAL::Kernel_traits<Surface_point>::Kernel Surface_kernel;
            typedef CGAL::Epick::Point_3 Point;
            auto to_epick = CGAL::Cartesian_converter<Surface_kernel, CGAL::Epick>{};
            std::cerr << "Corner pb: " << to_epick(mesh.point(v)) << " " << borders_around_v << " " << constraints_around_vertex.size() << std::endl;
        }
        assert(
            (borders_around_v == 2 && constraints_around_vertex.empty()) ||
            (borders_around_v == 2 && constraints_around_vertex.size() == 2 && constraints_around_vertex.count(*(constraints_around_vertex.begin())) == 2) ||
            (borders_around_v == 0 && constraints_around_vertex.size() == 2 && constraints_around_vertex.count(*(constraints_around_vertex.begin())) == 2)
        );
        const auto e = mesh.edge(removed);
        if (is_constrained_edge[e]) {
            external_constraints[edge_constraint_source[e]].erase(e);
        }
        is_constrained_edge[e] = false;
        edge_constraint_source[e] = Constraint{};
        CGAL::Euler::join_vertex(removed, mesh);
        constraints_around_vertex.clear();
    }

}

template <typename Surface>
void triangulate_faces(Surface& surface)
{

    typedef typename Surface::Surface_mesh Surface_mesh;
    typedef typename Surface::Surface_index Surface_index;
    //typedef typename Surface::Vertex_info Vertex_info;
    typedef typename Surface_mesh::Face_index Face_index;

    auto& mesh = surface.mesh;
    //auto vertex_info = surface.vertex_info_map();
    auto is_constrained_edge = surface.edge_constraints_map();
    auto edge_constraint_source = surface.edge_constraint_source_map();
    std::vector<Face_index> triangulation_candidates;
    for (auto&& f : mesh.faces()) {
        if (mesh.degree(f) > 3) {
            triangulation_candidates.emplace_back(f);
        }
    }
    for (auto&& f : triangulation_candidates) {
        auto h = mesh.halfedge(f);
        const auto P = face_barycenter(h, mesh);
        h = CGAL::Euler::add_center_vertex(h, mesh);
        const auto v = mesh.target(h);
        //vertex_info[v] = Vertex_info{ surface.index };
        assert(
            std::all_of(
                begin(CGAL::halfedges_around_source(v, mesh)),
                end(CGAL::halfedges_around_source(v, mesh)),
                [&mesh, &is_constrained_edge, &edge_constraint_source](auto h) {
            const auto e = mesh.edge(h);
            return !is_constrained_edge[e] && Surface::unspecified(edge_constraint_source[e]);
        }
            )
        );
        //for (auto&& h2 : CGAL::halfedges_around_source(v, mesh)) {
        //    const auto e = mesh.edge(h2);
        //    is_constrained_edge[e] = false;
        //    edge_constraint_source[e] = std::numeric_limits<Surface_index>::max();
        //}
        mesh.point(v) = P;
    }

}
template <typename Surface>
void simplify_connected(Surface& surface)
{

    remove_center_vertices(surface);
    remove_inner_edges(surface);
    simplify_face_borders(surface);
    triangulate_faces(surface);

}

template <typename SurfaceMesh, typename Location>
struct Property_accessor;

template <typename Point>
struct Property_accessor<CGAL::Surface_mesh<Point>, typename CGAL::Surface_mesh<Point>::Vertex_index>
{
    auto number_of_active_locations(const CGAL::Surface_mesh<Point>& mesh) const noexcept {
        return mesh.number_of_vertices();
    }
    auto active_locations(const CGAL::Surface_mesh<Point>& mesh) const noexcept {
        return mesh.vertices();
    }
};

template <typename Point>
struct Property_accessor<CGAL::Surface_mesh<Point>, typename CGAL::Surface_mesh<Point>::Face_index>
{
    auto number_of_active_locations(const CGAL::Surface_mesh<Point>& mesh) const noexcept {
        return mesh.number_of_faces();
    }
    auto active_locations(const CGAL::Surface_mesh<Point>& mesh) const noexcept {
        return mesh.faces();
    }
};

template <typename SurfaceMesh, typename PropertyMap>
auto property_accessor()
{
    return Property_accessor<SurfaceMesh, typename PropertyMap::key_type>{};
}

template <typename Mesh, typename PropertyMap>
auto propertymap_as_array(const Mesh& mesh, const PropertyMap& pmap)
{
    auto accessor = property_accessor<Mesh, PropertyMap>();
    const auto n = accessor.number_of_active_locations(mesh);
    typedef typename PropertyMap::value_type Value_type;
    auto result = py::array_t< Value_type, py::array::c_style >{ n };
    auto p = reinterpret_cast<Value_type*>(result.request().ptr);
    for (auto&& e : accessor.active_locations(mesh)) {
        *p = pmap[e];
        ++p;
    }
    return result;
}

template <typename Surface>
auto corefine(Surface& S1, Surface& S2) {

    typedef typename Surface::Surface_mesh Surface_mesh;
    typedef Surface_mesh::Vertex_index Vertex_index;
    typedef Surface_mesh::Edge_index Edge_index;
    auto constraints1 = S1.property_map<Edge_index, bool>(false);
    auto constraints2 = S2.property_map<Edge_index, bool>(false);
    namespace parameters = CGAL::Polygon_mesh_processing::parameters;
    Split_recorder::reset(S1.mesh, S2.mesh);
    CGAL::Polygon_mesh_processing::corefine(
        S1.mesh, S2.mesh,
        parameters::edge_is_constrained_map(constraints1),
        parameters::edge_is_constrained_map(constraints2)
    );
    std::cerr << "Recorded: ";
    for (auto&& m : Split_recorder::records_map) {
        std::cerr << m.second.size() << " ";
    }
    std::cerr << std::endl;
    auto apply_constraints = [](Surface& S, const auto& edge_is_constrained, Surface& other) {
        typedef typename Surface::Constraint Constraint;
        auto& mesh = S.mesh;
        auto has_constraint = S.edge_constraints_map();
        auto constraint_source = S.edge_constraint_source_map();
        const auto other_as_constraint = Constraint{ &other };
        auto is_corner = S.corners_map();
        auto& previously_constrained = Split_recorder::records_map[&mesh];
        for (auto&& p : previously_constrained) {
            assert(!mesh.is_removed(p.first));
            const auto eold = mesh.edge(p.first);
            if (has_constraint[eold]) {
                assert(!(constraint_source[eold] == other_as_constraint));
                const auto enew = mesh.edge(p.second);
                has_constraint[enew] = true;
                assert(!Surface::unspecified(constraint_source[eold]));
                constraint_source[enew] = constraint_source[eold];
                assert(!mesh.is_removed(mesh.source(p.first)));
                assert(!mesh.is_removed(mesh.target(p.second)));
                assert(mesh.source(p.first) == mesh.target(p.second));
                is_corner[mesh.target(p.second)] = true;
            }
        }
        auto& external_constraints = S.external_constraints(other_as_constraint);
        assert(external_constraints.empty());
        auto constrained_vertices = std::set<Vertex_index>{};
        for (auto&& e : mesh.edges()) {
            if (edge_is_constrained[e]) {
                if (has_constraint[e]) {
                    std::cerr << "Inconsistency on surface " << S.index
                        << " corefined with " << other.index
                        << " edge " << e << " already constrained by "
                        << Surface::human_readable(constraint_source[e]) << std::endl;
                    typedef typename Surface::Surface_mesh Surface_mesh;
                    typedef typename Surface_mesh::Point Surface_point;
                    typedef typename CGAL::Kernel_traits<Surface_point>::Kernel Surface_kernel;
                    typedef CGAL::Epick::Point_3 Point;
                    auto to_epick = CGAL::Cartesian_converter<Surface_kernel, CGAL::Epick>{};
                    auto f = std::ofstream("suspect_edge.csv");
                    for (int i = 0; i < 2; ++i) {
                        f << to_epick(mesh.point(mesh.vertex(e, i))) << std::endl;
                    }
                }
                assert(!has_constraint[e]);
                has_constraint[e] = true;
                constraint_source[e] = other_as_constraint;
                external_constraints.insert(e);
                for (int i = 0; i < 2; ++i) {
                    constrained_vertices.insert(mesh.vertex(e, i));
                }
            }
        }
        if (external_constraints.empty()) {
            S.external_constraints_map.erase(other_as_constraint);
        }
        std::multiset<Constraint> vertex_constraints;
        for (auto&& v : constrained_vertices) {
            if (!is_corner[v]) {
                for (auto&& h : CGAL::halfedges_around_source(v, mesh)) {
                    const auto e = mesh.edge(h);
                    if (has_constraint[e]) {
                        vertex_constraints.insert(constraint_source[e]);
                    }
                }
                if (!(vertex_constraints.size() == 2 && vertex_constraints.count(*(vertex_constraints.begin())) == 2)) {
                    is_corner[v] = true;
                }
                vertex_constraints.clear();
            }
        }
    };
    apply_constraints(S1, constraints1, S2);
    apply_constraints(S2, constraints2, S1);
    const auto n = std::count(std::begin(constraints1), std::end(constraints1), true);
    assert(std::count(std::begin(constraints2), std::end(constraints2), true) == n);
    S1.remove_property_map(constraints1);
    S2.remove_property_map(constraints2);
    return n;
}

template <typename Surface>
auto connected_components_as_array(Surface& surface)
{

    typedef typename Surface::Surface_mesh Surface_mesh;
    typedef typename Surface_mesh::Face_index Face_index;
    typedef typename Surface_mesh::faces_size_type Component_id;

    namespace parameters = CGAL::Polygon_mesh_processing::parameters;
    auto components = surface.property_map<Face_index, Component_id>();
    auto& mesh = surface.mesh;
    CGAL::Polygon_mesh_processing::connected_components(
        mesh, components,
        parameters::edge_is_constrained_map(surface.edge_constraints_map())
    );
    auto result = propertymap_as_array(mesh, components);
    surface.remove_property_map(components);
    return result;

}

typedef CGAL::Epick Kernel;
typedef typename Kernel::Point_3 Point;
typedef typename Kernel::Vector_3 Vector;
typedef typename Kernel::Plane_3 Plane;
typedef std::size_t Surface_index;
typedef Boundary_hull<CGAL::Epeck> Hull;
typedef typename Hull::Surface_mesh Hull_boundary;
typedef On_plane_side<typename Hull::Kernel::Plane_3> On_side;
typedef Surface_factory<typename Hull::Point, Surface_index> Factory;
typedef typename Factory::Surface Surface;

Factory surface_factory;

Factory& get_surface_factory()
{
    return surface_factory;
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

    py::class_<On_side>(module, "On_side")
        .def(py::init([](const Plane& plane, const Point& P) {
        auto to_hull_kernel = CGAL::Cartesian_converter<Kernel, typename Hull::Kernel>{};
        return std::make_unique<On_side>(
            to_hull_kernel(plane),
            to_hull_kernel(P)
            );
    }))
        .def("__call__", [](const On_side& self, const Point& P) {
        return self(CGAL::Cartesian_converter<Kernel, typename Hull::Kernel>{}(P));
        })
        ;

        py::class_<Hull_boundary>(module, "Hull_boundary")
            .def("as_arrays", &as_numpy_arrays<Hull_boundary>)
            ;
        
        py::class_<Surface>(module, "Surface")
            .def_readonly("index", &Surface::index)
            .def("as_arrays", [](Surface& self) {
            return as_numpy_arrays(self.mesh);
        })
            .def("simplify_connected_components", [](Surface& self) {
            simplify_connected(self);
        })
        //    .def("remove_constraint", [](Surface& self, Surface_index si) {
        //    self.remove_constraint(si);
        //})
        //    .def("keep", [](Surface& self, const On_side& mask) {
        //    auto relaxed = select_faces(self, mask);
        //    auto l = py::list{};
        //    for (auto&& si : relaxed) {
        //        l.append(si);
        //    }
        //    return l;
        //})
            .def("keep_connected", [](Surface& self, const On_side& mask) {
            auto relaxed = select_connected_components(self, mask);
            auto l = py::list{};
            for (auto&& c : relaxed) {
                l.append(Surface::human_readable(c));
            }
            return l;
        })
            .def("connected_components", [](Surface& self) {
            return connected_components_as_array(self);
        })
            .def("is_constrained_by", &Surface::is_constrained_by)
            .def("constraints", [](Surface& self) {
            auto result = py::list{};
            for (auto&& p : self.external_constraints_map) {
                result.append(py::make_tuple(Surface::human_readable(p.first), p.second.size()));
            }
            return result;
        })
            .def("number_of_constraints", &Surface::number_of_constraints)
            .def("constrained_midpoints", &constrained_midpoints<Surface>)
            .def("mark_all_vertices_as_corners", &Surface::mark_all_vertices_as_corners)
                ;

        py::class_<Factory>(module, "Surface_factory")
            .def("__call__", [](Factory& self) {
            return self.make_unique_pointer();
        });

        module.def("surface_factory", &get_surface_factory);

        py::class_<Hull>(module, "Hull")
        .def(py::init([](py::list plane_list) {
        std::vector<Plane> planes;
        planes.reserve(py::len(plane_list));
        for (auto& plane : plane_list) {
            planes.emplace_back(plane.cast<Plane>());
        }
        return std::make_unique<Hull>(cbegin(planes), cend(planes));
    }))
        .def("boundaries", [](const Hull& self) {
        return py::make_iterator(begin(self.boundaries), end(self.boundaries));
    }, py::keep_alive<0, 1>())
        .def("intersect", [](const Hull& self, const Plane& plane) {
        auto surface = surface_factory.make_unique_pointer();
        self.intersect(
            CGAL::Cartesian_converter<Kernel, typename Hull::Kernel>{}(plane),
            surface->mesh
        );
        return surface;
    })
            ;

    module.def("corefine", &corefine<Surface>)
        ;

    module.def("simplify", [](Hull_boundary& B) {
        simplify_faces(B);
    });

    module.def("simplify", [](Surface& S) {
        simplify_faces(S.mesh);
    });
    
    //module.def("simplify_connected_components", [](Surface& S) {
    //    simplify_connected_components(S);
    //});

}
