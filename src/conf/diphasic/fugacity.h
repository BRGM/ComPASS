#pragma once

#include "Thermodynamics.h"
#include "enum_to_rank.h"

// Fugacity function of Cpha = air fraction in phase ph
template <typename Components, typename Phases>
inline auto fugacity(const Components cp, const Phases ph, const double &p,
                     const double &T, const double &Cpha) {
   static_assert(std::is_enum_v<Components>);
   static_assert(std::is_enum_v<Phases>);
   assert(cp == Components::air || cp == Components::water);
   assert(ph == Phases::gas || ph == Phases::liquid);
   double f, _, Cph[2], dfdC[2];
   Cph[enum_to_rank(Components::air)] = Cpha;
   Cph[enum_to_rank(Components::water)] = 1 - Cpha;
   FluidThermodynamics_fugacity(to_underlying(cp), to_underlying(ph), p, T, Cph,
                                f, _, _, dfdC);
   return std::make_tuple(f, dfdC[enum_to_rank(Components::air)] -
                                 dfdC[enum_to_rank(Components::water)]);
}
