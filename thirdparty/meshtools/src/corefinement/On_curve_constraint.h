#pragma once

namespace TSurfBlobTraits
{

    template <typename Curve_type>
    struct On_curve_constraint
    {
        typedef Curve_type Curve;
        typedef typename Curve::Position Position;
        typedef typename Curve::Node Node;
        Curve * curve;
        Position position;
        On_curve_constraint() = delete;
        On_curve_constraint(Curve * p, Position pos) :
            curve{ p },
            position{ pos }
        {}
        bool is_valid() const noexcept {
            assert(curve);
            if (!curve->is_valid()) return false;
            if (!position->is_valid()) return false;
            for (auto p = curve->begin(); p != curve->end(); ++p) {
                if (p == position) return true;
            }
            return false;
        }
        bool is_on_extremity() const noexcept {
            return position == curve->begin() || next(position) == curve->end();
        }
        bool operator==(const On_curve_constraint& other) const noexcept {
            return curve == other.curve && position == other.position;
        }
    };

} // namespace TSurfBlobTraits
