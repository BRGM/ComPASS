#include <array>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Convex_hull_3.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Cartesian_converter<CGAL::Epick, CGAL::Epeck> IK_to_EK;
typedef CGAL::Cartesian_converter<CGAL::Epeck, CGAL::Epick> EK_to_IK;

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

template <typename Plane_iterator>
auto test(Plane_iterator first, Plane_iterator last)
{
    typedef CGAL::Epeck::Plane_3 Plane;
    auto planes = std::vector<Plane>{};
    auto to_exact = IK_to_EK{};
    std::transform(first, last, std::back_inserter(planes), to_exact);
    typedef CGAL::Epeck::Point_3 Point;
    const auto n = planes.size();
    std::vector<std::vector<Point>> corners(n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i+1; j < n; ++j) {
            for (std::size_t k = j+1; k < n; ++k) {
                auto intersection = CGAL::intersection(planes[i], planes[j], planes[k]);
                if (intersection) {
                    if (const Point * p = boost::get<Point>(&*intersection)) {
                        std::cerr << "Found " << EK_to_IK{}(*p) << std::endl;
                        if (std::all_of(begin(planes), end(planes),
                            [p](const auto& plane) { return !plane.has_on_positive_side(*p); })) {
                            std::cerr << "Added one" << EK_to_IK{}(*p) << std::endl;
                            corners[i].push_back(*p);
                            corners[j].push_back(*p);
                            corners[k].push_back(*p);
                        }
                    }
                }
            }
        }
    }
    typedef CGAL::Surface_mesh<Point> Surface_mesh;
    for (auto&& points : corners) {
        Surface_mesh sm;
        std::cerr << "Size corners" << points.size() << std::endl;
        CGAL::convex_hull_3(begin(points), end(points), sm);
        std::cout << "The convex hull contains " << sm.number_of_vertices() << " vertices and " << sm.number_of_faces() << "faces" << std::endl;
        for (auto&& face : sm.faces()) {

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
