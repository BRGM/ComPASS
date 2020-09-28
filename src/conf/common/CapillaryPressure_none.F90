! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module CapillaryPressure

   use, intrinsic :: iso_c_binding, only: c_double, c_int
   use DefModel, only: NbPhase, IndThermique

   implicit none

   public :: f_PressionCapillaire

contains

   ! P(iph) = Pref + f_PressionCapillaire(iph)
   !< rt is the rocktype identifier
   !< iph is the phase identifier : GAS_PHASE or LIQUID_PHASE
   !< S is all the saturations
   pure subroutine f_PressionCapillaire(rt, iph, S, f, DSf)

      integer(c_int), intent(in) :: rt(IndThermique + 1)
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: S(NbPhase)
      real(c_double), intent(out) :: f, DSf(NbPhase)

      f = 0.d0
      dSf(:) = 0.d0

   end subroutine f_PressionCapillaire

end module CapillaryPressure
