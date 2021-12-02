#pragma once

#include "fugacity.h"

template <typename Components, typename Phases, typename PhaseVector>
PhaseVector diphasic_equilibrium(const PhaseVector& pa, const double& T,
                                 const double atol = 1e-7,
                                 std::size_t maxiter = 1000) {
   static_assert(std::is_enum_v<Components>);
   static_assert(std::is_enum_v<Phases>);

   const double pg = pa[enum_to_rank(Phases::gas)];
   const double pl = pa[enum_to_rank(Phases::liquid)];
   double Cga = 1, Cla = 0;  // Newton start: no water in gas no air in liquid
   double f, df, Jag, Jal, Jwg, Jwl, det, Ra, Rw;

   for (; maxiter != 0; --maxiter) {
      std::tie(f, df) = fugacity(Components::air, Phases::gas, pg, T, Cga);
      Ra = f * Cga;
      Jag = df * Cga + f;
      std::tie(f, df) = fugacity(Components::air, Phases::liquid, pl, T, Cla);
      Ra -= f * Cla;
      Jal = -df * Cla - f;
      std::tie(f, df) = fugacity(Components::water, Phases::gas, pg, T, Cga);
      Rw = f * (1 - Cga);
      Jwg = df * (1 - Cga) - f;
      std::tie(f, df) = fugacity(Components::water, Phases::liquid, pl, T, Cla);
      Rw -= f * (1 - Cla);
      Jwl = -df * (1 - Cla) + f;
      if (fabs(Ra) + fabs(Rw) < atol) return PhaseVector{{Cga, Cla}};
      // Newton: X -> X - J^-1 R
      det = Jag * Jwl - Jwg * Jal;
      assert(det != 0);
      Cga -= (Jwl * Ra - Jal * Rw) / det;
      Cla -= (-Jwg * Ra + Jag * Rw) / det;
   }

   throw std::runtime_error(
       "Maximum number of iterations in diphasic_equilibrium exceeded.");
}
