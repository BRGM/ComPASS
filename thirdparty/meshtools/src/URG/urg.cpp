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
        boundaries.resize(n);
        for (Plane_index i = 0; i < n; ++i) {
            const auto& plane_corners = corners[i];
            std::cerr << "Size corners" << plane_corners.size() << std::endl;
            std::vector<Point> points;
            points.reserve(plane_corners.size());
            for (const auto& corner : plane_corners) points.push_back(corner.point);
            auto& boundary = boundaries[i];
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

//typedef Boundary_hull<CGAL::Epeck> Hull;
//typedef typename Hull::Surface_mesh Surface_mesh;
//struct Mask
//{
//    template <typename T>
//    bool operator()(const T&) const { return true; }
//};

template <typename Surface_mesh, typename Mask>
void select_faces(Surface_mesh& mesh, const Mask& is_kept)
{

    typedef typename Surface_mesh::Point Surface_point;
    typedef typename CGAL::Kernel_traits<Surface_point>::Kernel Surface_kernel;
    typedef typename Surface_kernel::Triangle_3 Surface_triangle;

    // special type to be used with CGAL:barycenter
    typedef std::pair<Surface_point, typename Surface_kernel::FT> Face_point;
    std::vector<Face_point> face_points;
    for (auto&& face : mesh.faces()) {
        face_points.clear();
        for (auto&& v : CGAL::vertices_around_face(mesh.halfedge(face), mesh)) {
            face_points.emplace_back(mesh.point(v), 1);
        }
        if (!is_kept(CGAL::barycenter(cbegin(face_points), cend(face_points)))) {
            mesh.remove_face(face);
        }
    }
    mesh.collect_garbage();

}


PYBIND11_MODULE(URG, module)
{

    module.doc() = "plane intersection for URG model";

    typedef CGAL::Epick Kernel;
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::Plane_3 Plane;

    py::class_<Point>(module, "Point")
        .def(py::init<double, double, double>())
        ;

    py::class_<Vector>(module, "Vector")
        .def(py::init<double, double, double>())
        ;

    py::class_<Plane>(module, "Plane")
        .def(py::init<Point, Vector>())
        ;

    typedef Boundary_hull<CGAL::Epeck> Hull;
    typedef typename Hull::Surface_mesh Surface_mesh;
    typedef On_plane_side<typename Hull::Kernel::Plane_3> On_side;

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

    py::class_<Surface_mesh>(module, "Surface")
        .def("as_arrays", &as_numpy_arrays<Surface_mesh>)
        .def("keep", [](Surface_mesh& self, const On_side& mask) {
        select_faces(self, mask);
    })
        ;

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
        auto mesh = std::make_unique<Surface_mesh>();
        self.intersect(
            CGAL::Cartesian_converter<Kernel, typename Hull::Kernel>{}(plane),
            *mesh
        );
        return mesh;
    })
            ;

    module.def("corefine", [](Surface_mesh& mesh1, Surface_mesh& mesh2) {
        CGAL::Polygon_mesh_processing::corefine(mesh1, mesh2);
    });

}
