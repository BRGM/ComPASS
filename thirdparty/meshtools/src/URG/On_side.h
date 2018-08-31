#pragma once

#include <CGAL/enum.h>
#include <CGAL/Kernel_traits.h>

template <typename Sided_oject_type>
struct On_side
{
    typedef Sided_oject_type Sided_oject;
    typedef typename Sided_oject::Kernel Kernel;
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Segment_3 Segment;

    Sided_oject sided_object;
    CGAL::Oriented_side side;
    On_side(const Sided_oject& S, const Point& A) :
        sided_object{ S },
        side{ S.oriented_side(A) }
    {}
    On_side(Sided_oject&& S, const Point& A) :
        sided_object{ std::forward<Sided_oject>(S) },
        side{ S.oriented_side(A) }
    {}
    bool operator()(const Point& P) const noexcept {
        // negation to include boundary
        return !(sided_object.oriented_side(P) == CGAL::opposite(side));
    }
    // FIXME: could return several pieces
    auto clip(const Segment& S) const noexcept {
        typedef boost::optional<Segment> Result;
        const auto pieces = sided_object.split(S);
        assert(pieces.size() > 0);
        for (auto&& Si : pieces) {
            if (this->operator()(CGAL::midpoint(Si.source(), Si.target()))) {
                return Result{ Si };
            }
        }
        return Result{};
    }
};