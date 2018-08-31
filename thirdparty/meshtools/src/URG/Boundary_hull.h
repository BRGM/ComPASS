#pragma once

#include <array>
#include <vector>

#include <CGAL/Kernel_traits.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Surface_mesh.h>

//FIXME: prefer standard compliant optional
#include <boost/optional.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "Indexed_elements.h"

template <typename Kernel_type>
struct Boundary_hull
{

    typedef Kernel_type Kernel;
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Plane_3 Plane;
    typedef typename Kernel::Line_3 Line;
    typedef typename Kernel::Segment_3 Segment;
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
            //std::cerr << "Edge (" << index[0] << "," << index[1] << ") with" << points.size() << " points" << std::endl;
            assert(points.size() == 2);
            // FIXME: could be emplace_back here but does not work on MSVC 2015
            edges.push_back(Edge{ Segment{ points.front().point, points.back().point }, index });
        }
        edges.shrink_to_fit();
        //std::cerr << "Collected " << edges.size() << " hull edges" << std::endl;
        assert(boundaries.empty());
        boundaries.reserve(n);
        for (Plane_index i = 0; i < n; ++i) {
            const auto& plane_corners = corners[i];
            //std::cerr << "Size corners" << plane_corners.size() << std::endl;
            std::vector<Point> points;
            points.reserve(plane_corners.size());
            for (const auto& corner : plane_corners) points.push_back(corner.point);
            boundaries.emplace_back();
            auto& boundary = boundaries.back();
            CGAL::convex_hull_3(cbegin(points), cend(points), boundary);
            //std::cout << "The convex hull contains " << boundary.number_of_vertices() << " vertices and " << boundary.number_of_faces() << " faces" << std::endl;
            //for (const auto& face : boundary.faces()) {
            //    std::cerr << "Face with " << boundary.degree(face) << " points" << std::endl;
            //}
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
        //std::cerr << "Hull intersection with " << points.size() << " points." << std::endl;
        CGAL::convex_hull_3(cbegin(points), cend(points), mesh);
    }

    auto intersect(const Line& l) const
    {
        typedef boost::optional<boost::variant<Point, Segment>> Result_type;
        std::vector<Point> points;
        points.reserve(2);
        for (const auto& P : planes) {
            auto intersection = CGAL::intersection(P, l);
            if (intersection) {
                if (const Point * p = boost::get<Point>(&*intersection)) {
                    if (is_inside(*p)) {
                        points.push_back(*p);
                    }
                }
            }
        }
        if (points.size() == 2) {
            if (CGAL::squared_distance(points.front(), points.back()) == 0)
                return Result_type{ points.front() };
            return Result_type{ Segment{ points.front(), points.back() } };
        }
        if (points.empty()) 
            return Result_type{};
        std::vector<Point> filtered_points;
        for (auto&& p : points) {
            bool add = true;
            for (auto&& fp : filtered_points) {
                if (CGAL::squared_distance(p, fp) == 0) {
                    add = false;
                    break;
                }
            }
            if (add) {
                filtered_points.push_back(p);
            }
        }
        assert(!filtered_points.empty());
        assert(filtered_points.size()<=2);
        if(filtered_points.size()==1)
            return Result_type{ points.front() };
        return Result_type{ Segment{ points.front(), points.back() } };
    }

};
