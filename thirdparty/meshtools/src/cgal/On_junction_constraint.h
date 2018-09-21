#pragma once

#include "Curves_junction.h"

namespace TSurfBlobTraits
{

    template <typename Curve_type>
    struct On_junction_constraint
    {
        typedef Curve_type Curve;
        typedef Curves_junction<Curve> Junction;
        std::vector<Junction> junctions;
        bool has_incoming_curve(const Curve * p) const noexcept {
            for (auto&& junction : junctions) {
                if (junction.first == p) return true;
            }
            return false;
        }
        bool has_outgoing_curve(const Curve * p) const noexcept {
            for (auto&& junction : junctions) {
                if (junction.second == p) return true;
            }
            return false;
        }
        bool is_valid() const noexcept {
            if (junctions.size() < 2) return false;
            if (junctions.size() == 2) {
                return junctions[0].is_weak() && junctions[1].is_weak();
            }
            for (auto&& junction : junctions) {
                if (junction.is_weak()) return false;
            }
            return true;
        }
        bool has_weak_junction() const {
            assert(is_valid());
            return junctions[1].is_weak();
        }
    };

} // namespace TSurfBlobTraits
