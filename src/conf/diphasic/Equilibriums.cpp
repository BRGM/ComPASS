#include "Equilibriums.h"

#include <cmath>
#include <stdexcept>

#include "DefModel.h"
#include "Thermodynamics.h"
#include "enum_to_rank.h"

// Fugacity function of Cpha = air fraction in phase ph
template <Component cp, Phase ph>
inline auto fugacity(const Real &p, const Real &T, const Real &Cpha) {
   Real f, _, Cph[2], dfdC[2];
   Cph[enum_to_rank(Component::air)] = Cpha;
   Cph[enum_to_rank(Component::water)] = 1 - Cpha;
   FluidThermodynamics_fugacity(static_cast<int>(cp), static_cast<int>(ph), p,
                                T, Cph, f, _, _, dfdC);
   return std::make_tuple(
       f, dfdC[enum_to_rank(Phase::gas)] - dfdC[enum_to_rank(Phase::liquid)]);
}

Phase_vector diphasic_equilibrium(const Phase_vector &pa, const Real &T,
                                  const double atol, std::size_t maxiter) {
   const Real pg = pa[enum_to_rank(Phase::gas)];
   const Real pl = pa[enum_to_rank(Phase::liquid)];
   Real Cga = 1, Cla = 0;  // Newton start: no water in gas no air in liquid
   Real f, df, Jag, Jal, Jwg, Jwl, det, Ra, Rw;

   for (; maxiter != 0; --maxiter) {
      std::tie(f, df) = fugacity<Component::air, Phase::gas>(pg, T, Cga);
      Ra = f * Cga;
      Jag = df * Cga + f;
      std::tie(f, df) = fugacity<Component::air, Phase::liquid>(pl, T, Cla);
      Ra -= f * Cla;
      Jal = -df * Cla - f;
      std::tie(f, df) = fugacity<Component::water, Phase::gas>(pg, T, Cga);
      Rw = f * (1 - Cga);
      Jwg = df * (1 - Cga) - f;
      std::tie(f, df) = fugacity<Component::water, Phase::liquid>(pl, T, Cla);
      Rw -= f * (1 - Cla);
      Jwl = -df * (1 - Cla) + f;
      if (fabs(Ra) + fabs(Rw) < atol) return Phase_vector{{Cga, Cla}};
      // Newton: X -> X - J^-1 R
      det = Jag * Jwl - Jwg * Jal;
      assert(det != 0);
      Cga -= (Jwl * Ra - Jal * Rw) / det;
      Cla -= (-Jwg * Ra + Jag * Rw) / det;
   }

   throw std::runtime_error(
       "Maximum number of iterations in diphasic_equilibrium exceeded.");
}
