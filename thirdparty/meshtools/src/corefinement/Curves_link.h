#pragma once

#include "On_curve_constraint.h"

namespace TSurfBlobTraits
{

    template <typename Curve_type>
    struct Curves_link : std::pair<On_curve_constraint<Curve_type>, On_curve_constraint<Curve_type>>
    {
        typedef Curve_type Curve;
        //typedef On_curve_constraint<Curve> Constraint;
        //typedef typename Constraint_type::Curve Curve;
        bool is_weak() const {
            assert(first.curve != nullptr);
            assert(second.curve != nullptr);
            return first.curve == second.curve;
        }
        bool is_valid() const {
            assert(first.is_valid());
            assert(second.is_valid());
            assert(first.curve != nullptr);
            assert(second.curve != nullptr);
            if (first.curve == second.curve) { // weak edge to be cut later
                return next(first.position) == second.position;
            }
            if (!first.is_on_extremity()) return false;
            if (!second.is_on_extremity()) return false;
            if (!((*first.position) == (*second.position))) return false;
            return true;
        }
    };

} // namespace TSurfBlobTraits
