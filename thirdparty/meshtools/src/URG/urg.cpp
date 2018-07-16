#include <array>
#include <vector>

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

template <typename Plane_iterator>
auto test(Plane_iterator first, Plane_iterator last)
{
    typedef CGAL::Epeck Kernel;
    typedef CGAL::Cartesian_converter<CGAL::Epick, Kernel> IK_to_EK;
    typedef CGAL::Cartesian_converter<Kernel, CGAL::Epick> EK_to_IK;
    typedef Kernel::Plane_3 Plane;
    typedef std::vector<Plane> Plane_collection;
    typedef Plane_collection::size_type Plane_index;
    auto planes = std::vector<Plane>{};
    auto to_exact = IK_to_EK{};
    std::transform(first, last, std::back_inserter(planes), to_exact);
    typedef Kernel::Point_3 Point;
    typedef std::array<Plane_index, 3> Corner_index;
    typedef Indexed_point<Point, Corner_index> Corner;
    const auto n = planes.size();
    std::vector<std::vector<Corner>> corners(n);
    typedef std::array<Plane_index, 2> Edge_index;
    std::map<Edge_index, std::vector<Corner>> edge_map;
    for (Plane_index i = 0; i < n; ++i) {
        for (Plane_index j = i+1; j < n; ++j) {
            for (Plane_index k = j+1; k < n; ++k) {
                auto intersection = CGAL::intersection(planes[i], planes[j], planes[k]);
                if (intersection) {
                    if (const Point * p = boost::get<Point>(&*intersection)) {
                        std::cerr << "Found " << EK_to_IK{}(*p) << std::endl;
                        if (std::all_of(begin(planes), end(planes),
                            [p](const auto& plane) { return !plane.has_on_positive_side(*p); })) {
                            std::cerr << "Added one" << EK_to_IK{}(*p) << std::endl;
                            const auto P = Corner{ *p, {i, j, k} };
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
    typedef Kernel::Segment_3 Segment;
    typedef Indexed_segment<Segment, Edge_index> Edge;
    typedef std::vector<Edge> Edge_collection;
    Edge_collection edges;
    edges.reserve(edge_map.size());
    for (const auto& edge : edge_map) {
        const auto& index = edge.first;
        const auto& points = edge.second;
        std::cerr << "Edge(" << index[0] << "," << index[1] << ") with" << points.size() << "points" << std::endl;
        assert(points.size() == 2);
        // FIXME: could be emplace_back here but does not work on MSVC 2015
        edges.push_back(Edge{ Segment{ points.front().point, points.back().point }, index });
    }
    typedef CGAL::Surface_mesh<Point> Surface_mesh;
    for (const auto& plane_corners : corners) {
        Surface_mesh sm;
        std::cerr << "Size corners" << plane_corners.size() << std::endl;
        std::vector<Point> points;
        points.reserve(plane_corners.size());
        for (const auto& corner : plane_corners) points.push_back(corner.point);
        CGAL::convex_hull_3(begin(points), end(points), sm);
        std::cout << "The convex hull contains " << sm.number_of_vertices() << " vertices and " << sm.number_of_faces() << "faces" << std::endl;
        for (const auto& face : sm.faces()) {
            std::cerr << "Face with " << sm.degree(face) << " points" << std::endl;
        }
    }
    std::cerr << "Done" << std::endl;
}

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
    test(cbegin(planes), cend(planes));
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
