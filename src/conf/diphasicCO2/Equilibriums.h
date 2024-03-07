#pragma once

#include "Thermodynamics.h"
#include "enum_to_rank.h"

// there is no gas water, the equilibrium writes for the CO2 component:
// fug(CO2, gas) = fug(CO2, liq),
// fug(CO2, gas) is known and fug(CO2, liq) = coeff_CO2_liquid_fug * Cl(CO2)
template <typename Components, typename Phases, typename PhaseVector>
double CO2_liquid_molar_fraction(const PhaseVector& pa, const double& T) {
   static_assert(std::is_enum_v<Components>);
   static_assert(std::is_enum_v<Phases>);

   const double pg = pa[enum_to_rank(Phases::gas)];
   const double pl = pa[enum_to_rank(Phases::liquid)];
   double Cg[2], fg, fl, _, dfdC[2];

   Cg[enum_to_rank(Components::CO2)] = 1.e0;    // there is no gas water
   Cg[enum_to_rank(Components::water)] = 0.e0;  // there is no gas water

   // fug(CO2, gas)
   FluidThermodynamics_fugacity(to_underlying(Components::CO2),
                                to_underlying(Phases::gas), pg, T, Cg, fg, _, _,
                                dfdC);

   // Fugacity is coeff_CO2_liquid_fug * Cl(CO2)
   // coeff_CO2_liquid_fug does not depend on C
   FluidThermodynamics_coeff_CO2_liquid_fug(pl, T, fl, _, _);

   // return Cl(CO2) s.t. fug(CO2, gas) = coeff_CO2_liquid_fug * Cl(CO2)
   return fg / fl;
}
