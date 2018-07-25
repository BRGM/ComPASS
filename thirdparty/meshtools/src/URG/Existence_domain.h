#pragma once

#include "On_side.h"

template <typename Sided_object_type>
struct Existence_domain
{
    typedef Sided_object_type Sided_object;
    typedef On_side<Sided_object> On_side_type;
    typedef typename Sided_object::Kernel Kernel;
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Segment_3 Segment;

    std::vector<On_side_type> sides;
    void add_side(const Sided_object& S, const Point& A) {
        sides.emplace_back(On_side_type{ S, A });
    }
    void add_side(Sided_object&& S, const Point& A) {
        sides.emplace_back(On_side_type{ std::forward<Sided_object>(S), A });
    }
    bool operator()(const Point& P) const noexcept {
        return std::all_of(
            begin(sides), end(sides),
            [&P](auto side) { return side(P); }
        );
    }
    // FIXME: could return several pieces
    auto clip(const Segment& S) const noexcept {
        typedef boost::optional<Segment> Result;
        auto result = Result{S};
        for (auto&& side : sides) {
            result = side.clip(*result);
            if (!result) break;
        }
        return result;
    }
};
