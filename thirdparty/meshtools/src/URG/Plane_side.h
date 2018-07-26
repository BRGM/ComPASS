#pragma once

#include <vector>
#include <CGAL/Kernel_traits.h>

template <typename Plane_type>
struct Indexed_plane
{
    typedef Plane_type Plane;
    typedef typename CGAL::Kernel_traits<Plane>::Kernel Kernel;
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Segment_3 Segment;

    const Plane * plane;
    Indexed_plane() = delete;
    Indexed_plane(const Indexed_plane&) = default;
    Indexed_plane(Indexed_plane&&) = default;
    Indexed_plane& operator=(const Indexed_plane&) = default;
    auto oriented_side(const Point& P) const noexcept {
        assert(plane != nullptr);
        return plane->oriented_side(P);
    }
    auto split(const Segment& S) const noexcept {
        assert(plane != nullptr);
        assert(S.squared_length() > 0);
        std::vector<Segment> result;
        auto I = CGAL::intersection(*plane, S);
        if (I) {
            if (const Point * P = boost::get<Point>(&*I)) {
                if (CGAL::squared_distance(S.source(), *P)>0)
                    result.emplace_back(S.source(), *P);
                if (CGAL::squared_distance(*P, S.target())>0)
                    result.emplace_back(*P, S.target());
            }
            else {
                assert(boost::get<Segment>(&*I) != nullptr);
                result.emplace_back(S);
            }
        }
        else {
            result.emplace_back(S);
        }
        return result;
    }
};

