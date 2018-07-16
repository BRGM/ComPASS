#include <array>
#include <vector>

#include <CGAL/Kernel_traits.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Convex_hull_3.h>
#include <CGAL/Surface_mesh.h>

//FIXME: prefer standard compliant optional
#include <boost/optional.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

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

template <typename Kernel>
struct Boundary_hull
{

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
        Edge_collection edges;
        edges.reserve(edge_map.size());
        for (const auto& edge : edge_map) {
            const auto& index = edge.first;
            const auto& points = edge.second;
            std::cerr << "Edge (" << index[0] << "," << index[1] << ") with" << points.size() << " points" << std::endl;
            assert(points.size() == 2);
            // FIXME: could be emplace_back here but does not work on MSVC 2015
            edges.push_back(Edge{ Segment{ points.front().point, points.back().point }, index });
        }
        std::cerr << "Collected " << edges.size() << " edges" << std::endl;
        boundaries.resize(n);
        for (Plane_index i = 0; i < n; ++i) {
            const auto& plane_corners = corners[i];
            std::cerr << "Size corners" << plane_corners.size() << std::endl;
            std::vector<Point> points;
            points.reserve(plane_corners.size());
            for (const auto& corner : plane_corners) points.push_back(corner.point);
            auto& boundary = boundaries[i];
            CGAL::convex_hull_3(begin(points), end(points), boundary);
            std::cout << "The convex hull contains " << boundary.number_of_vertices() << " vertices and " << boundary.number_of_faces() << " faces" << std::endl;
            for (const auto& face : boundary.faces()) {
                std::cerr << "Face with " << boundary.degree(face) << " points" << std::endl;
            }
        }
    }
};

void test_cube()
{
    typedef CGAL::Epick::Point_3 Point;
    typedef CGAL::Epick::Vector_3 Vector;
    typedef CGAL::Epick::Plane_3 Plane;
    // cube
    auto planes = std::vector<Plane>{};
    planes.emplace_back(Point(-1, 0, 0), Vector(-1, 0, 0));
    planes.emplace_back(Point(1, 0, 0), Vector(1, 0, 0));
    planes.emplace_back(Point(0, -1, 0), Vector(0, -1, 0));
    planes.emplace_back(Point(0, 1, 0), Vector(0, 1, 0));
    planes.emplace_back(Point(0, 0, -1), Vector(0, 0, -1));
    planes.emplace_back(Point(0, 0, 1), Vector(0, 0, 1));
    Boundary_hull<CGAL::Epeck> hull{ cbegin(planes), cend(planes) };
    std::cerr << "Done" << std::endl;
}


PYBIND11_MODULE(URG, module)
{

    module.doc() = "plane intersection for URG model";
    
    module.def("test_cube", &test_cube);
    //module.def("split_and_mesh", &split_and_mesh);
    //module.def("foo", []() {
    //    boost::optional<int> test;
    //    test = false;
    //    if (test) py::print("ya bon");
    //});
}
