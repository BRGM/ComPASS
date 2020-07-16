!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module RelativePermeabilities

   use, intrinsic :: iso_c_binding, only: c_double, c_int
   use DefModel, only: NbPhase, IndThermique

   implicit none

   public :: f_PermRel

contains

   ! Permeabilites
   !< rt is the rocktype identifier
   !< iph is an identifier for each phase, here only one phase: LIQUID_PHASE
   !< S is all the saturations
   pure subroutine f_PermRel(rt, iph, S, kr, dkrdS)

      integer(c_int), intent(in) :: rt(IndThermique + 1)
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: S(NbPhase)
      real(c_double), intent(out) :: kr
      real(c_double), intent(out) :: dkrdS(NbPhase)

      kr = 1.d0
      dkrdS = 0.d0

   end subroutine f_PermRel

end module RelativePermeabilities
