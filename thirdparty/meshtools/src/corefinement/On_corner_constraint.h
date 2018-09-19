#pragma once

#include "Curves_link.h"

namespace TSurfBlobTraits
{

    template <typename Curve_type>
    struct On_corner_constraint
    {
        typedef Curve_type Curve;
        typedef Curves_link<Curve> Link;
        std::array<Link, 3> links;
        bool has_curve(const Curve * p) const noexcept {
            for (auto&& link : links) {
                if (link.first == p || link.second == p) return true;
            }
            return false;
        }
        bool is_valid() const noexcept {
            return std::all_of(begin(links), end(links),
                [](auto l) { return l.is_valid(); });
        }
        bool has_weak_link() const {
            assert(links[0].is_valid);
            assert(
                links[1].is_weak() ||
                (links[1].is_valid() && links[2].is_weak())
            );
            return links[1].is_weak() || links[2].is_weak();
        }
    };

} // namespace TSurfBlobTraits
