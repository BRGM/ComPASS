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

//FIXME: prefer standard compliant optional
#include <boost/optional.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "common.h"

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
    Surface_vertex_index& add(Surface_index si) {
        index.insert(si);
        return *this;
    }
    auto degree() const {
        return index.size();
    }
};

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

template <typename Surface_mesh_type, typename Surface_index_type>
struct Surface_factory
{
    typedef Surface_mesh_type Surface_mesh;
    typedef Surface_index_type Surface_index;
    struct Surface
    {
        typedef Surface_mesh_type Surface_mesh;
        typedef Surface_index_type Surface_index;
        typedef Surface_vertex_index<Surface_index> Vertex_degree;
        constexpr static const auto vertex_degree_map_name = "vertex_degree";
        constexpr static const auto edge_constraint_map_name = "is_constrained_edge";
        Surface_mesh mesh;
        Surface_index index;
        friend struct Surface_factory<Surface_mesh, Surface_index>;
        Surface() = delete;
        Surface(Surface_index si) :
            mesh{},
            index{ si }
        {
                auto vpmap = mesh.add_property_map<
                    typename Surface_mesh::Vertex_index, Vertex_degree
                >(vertex_degree_map_name, Vertex_degree{});
                assert(vpmap.second);
                auto epmap = mesh.add_property_map<
                    typename Surface_mesh::Edge_index, bool
                >(edge_constraint_map_name, false);
                assert(epmap.second);
        }
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

template <typename Surface_mesh>
auto as_numpy_arrays(const Surface_mesh& mesh)
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
            assert(mesh.degree(f) == 3);
            for (auto&& v : CGAL::vertices_around_face(mesh.halfedge(f), mesh)) {
                *(p++) = reindex[v];
            }
        }
    }
    return py::make_tuple(vertices, triangles);
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

template <typename Surface_mesh, typename Mask>
inline void select_faces(Surface_mesh& mesh, const Mask& is_kept)
{

    for (auto&& face : mesh.faces()) {
        auto h = mesh.halfedge(face);
        if (!is_kept(face_barycenter(h, mesh))) {
            CGAL::Euler::remove_face(h, mesh);
        }
    }
    mesh.collect_garbage();

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

typedef CGAL::Epick Kernel;
typedef typename Kernel::Point_3 Point;
typedef typename Kernel::Vector_3 Vector;
typedef typename Kernel::Plane_3 Plane;
typedef std::size_t Surface_index;
typedef Boundary_hull<CGAL::Epeck> Hull;
typedef typename Hull::Surface_mesh Hull_boundary;
typedef On_plane_side<typename Hull::Kernel::Plane_3> On_side;
typedef Surface_factory<Hull_boundary, Surface_index> Factory;
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
        .def("keep", [](Surface& self, const On_side& mask) {
            select_faces(self.mesh, mask);
        })
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

    module.def("corefine", [](Surface& S1, Surface& S2) {
        typedef typename Surface::Surface_mesh Surface_mesh;
        typedef Surface_mesh::Edge_index Edge_index;
        auto& mesh1 = S1.mesh;
        auto& mesh2 = S2.mesh;
        const auto map_name = std::string{ "intersections_constraint" };
        auto tmp_map1 = mesh1.add_property_map<Edge_index, bool>(map_name, false);
        assert(tmp_map1.second);
        auto tmp_map2 = mesh2.add_property_map<Edge_index, bool>(map_name, false);
        assert(tmp_map2.second);
        namespace parameters = CGAL::Polygon_mesh_processing::parameters;
        CGAL::Polygon_mesh_processing::corefine(mesh1, mesh2,
            parameters::edge_is_constrained_map(tmp_map1.first),
            parameters::edge_is_constrained_map(tmp_map2.first)
        );
        std::cerr << "Map sizes: " << std::count(std::begin(tmp_map1.first), std::end(tmp_map1.first), true) << std::endl;
        std::cerr << "Map sizes: " << std::count(std::begin(tmp_map2.first), std::end(tmp_map2.first), true) << std::endl;
    });

    module.def("simplify", [](Hull_boundary& B) {
        simplify_faces(B);
    });

    module.def("simplify", [](Surface& S) {
        simplify_faces(S.mesh);
    });

}
