#pragma once

#include "Curves_link.h"

namespace TSurfBlobTraits
{

    template <typename Curve_type>
    struct On_junction_constraint
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
            if (!links[0].is_valid()) return false;
            if (links[0].is_weak()) return false;
            if (!links[1].is_valid()) return false;
            if (!links[1].is_weak()) {
                if (!links[2].is_valid()) return false;
                if (links[2].is_weak()) return false;
            }
            return true;
        }
        bool has_weak_link() const {
            assert(is_valid());
            return links[1].is_weak();
        }
    };

} // namespace TSurfBlobTraits
