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

   ! Permeabilites = S**2
   !< rt is the rocktype identifier
   !< iph is the phase identifier : GAS_PHASE or LIQUID_PHASE
   !< S is all the saturations
   ! No interaction between phases
#ifdef NDEBUG
   pure &
#endif
      subroutine f_PermRel(rt, iph, S, kr, dkrdS)

      integer(c_int), intent(in) :: rt(IndThermique + 1)
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: S(NbPhase)
      real(c_double), intent(out) :: kr
      real(c_double), intent(out) :: dkrdS(NbPhase)

      kr = S(iph)
      dkrdS = 0.d0
      dkrdS(iph) = 1.d0

#ifndef NDEBUG
      if (iph /= GAS_PHASE .and. iph /= LIQUID_PHASE) &
         call CommonMPI_abort('unknow phase in f_PermRel')
      if (any(S < 0.d0) .or. any(S > 1.d0)) &
         call CommonMPI_abort('Unconsistent saturations in f_PermRel')
#endif

   end subroutine f_PermRel

end module RelativePermeabilities
