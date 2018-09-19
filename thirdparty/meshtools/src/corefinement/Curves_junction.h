#pragma once

#include "On_curve_constraint.h"

namespace TSurfBlobTraits
{

    template <typename Curve_type>
    struct Curves_junction : std::pair<On_curve_constraint<Curve_type>, On_curve_constraint<Curve_type>>
    {
        typedef Curve_type Curve;
        typedef On_curve_constraint<Curve> On_curve;
        //typedef typename Constraint_type::Curve Curve;
        Curves_junction() = delete;
        Curves_junction(const On_curve& c1, const On_curve& c2) :
            std::pair<On_curve, On_curve>{ c1, c2 }
        {}
        Curves_junction(On_curve&& c1, On_curve&& c2) :
            std::pair<On_curve, On_curve>{
                std::forward<On_curve>(c1),
                std::forward<On_curve>(c2)
            } 
        {}
        bool is_weak() const noexcept {
            assert(is_valid());
            return this->first.position->is_weak_corner();
        }
        bool is_valid() const noexcept {
            assert(this->first.curve);
            assert(this->second.curve);
            assert(this->first.curve!=this->second.curve);
            if (!this->first.is_on_extremity()) return false;
            if (!this->second.is_on_extremity()) return false;
            if (!this->first.curve->is_valid()) return false;
            if (!this->second.curve->is_valid()) return false;
            if (!((*this->first.position) == (*this->second.position))) return false;
            if (this->first.position->is_weak_corner()) {
               if(!this->second.position->is_weak_corner()) return false;
               return true;
            }
            if (!(this->first.position->point() == this->second.position->point())) return false;
            return true;
        }
    };

} // namespace TSurfBlobTraits
