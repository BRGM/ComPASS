#pragma once

#include "fugacity.h"

// there is no gas water, the equilibrium writes for the CO2 component:
// fug(CO2, gas) = fug(CO2, liq)
template <typename Components, typename Phases, typename PhaseVector>
double CO2_liquid_molar_fraction(const PhaseVector& pa, const double& T,
                                 const double atol = 1e-7,
                                 std::size_t maxiter = 1000) {
   static_assert(std::is_enum_v<Components>);
   static_assert(std::is_enum_v<Phases>);

   const double pg = pa[enum_to_rank(Phases::gas)];
   const double pl = pa[enum_to_rank(Phases::liquid)];
   double Cga = 1.e0;  // there is no gas water
   double Cla = 0.e0;  // Newton start: no CO2 in liquid
   double fg, dfg, f, df, Jag, Ra;

   // constant during the loop, because Cga is cst
   std::tie(fg, dfg) = fugacity(Components::air, Phases::gas, pg, T, Cga);
   for (; maxiter != 0; --maxiter) {
      std::tie(f, df) = fugacity(Components::air, Phases::liquid, pl, T, Cla);
      //   residual = fug(CO2, gas) - fug(CO2, liq)
      Ra = fg - f;
      //   derivative wrt CO2
      Jag = dfg - df;
      if (fabs(Ra) < atol) return Cla;
      // Newton: X -> X - J^-1 R
      assert(Jag != 0);
      Cla -= Ra / Jag;
   }

   throw std::runtime_error(
       "Maximum number of iterations in CO2_liquid_molar_fraction exceeded.");
}
