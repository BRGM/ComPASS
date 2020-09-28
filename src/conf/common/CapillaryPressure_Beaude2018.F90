! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module CapillaryPressure

   use, intrinsic :: iso_c_binding, only: c_double, c_int
   use DefModel, only: &
      NbPhase, IndThermique, &
      GAS_PHASE, LIQUID_PHASE

   implicit none

   public :: f_PressionCapillaire

contains

   ! P(iph) = Pref + f_PressionCapillaire(iph)
   !< rt is the rocktype identifier
   !< iph is the phase identifier : GAS_PHASE or LIQUID_PHASE
   !< S is all the saturations
   ! FIXME: IF f_PressionCapillaire DEPENDS ON THE ROCKTYPE,
   ! MODIFY f_EnergieInterne AND f_DensiteMolaire
   pure subroutine f_PressionCapillaire(rt, iph, S, f, DSf)

      integer(c_int), intent(in) :: rt(IndThermique + 1)
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: S(NbPhase)
      real(c_double), intent(out) :: f, DSf(NbPhase)

      real(c_double) :: Pc_cst, Sg0, A
      real(c_double) :: Sg

      dSf = 0.d0

      if (iph == GAS_PHASE) then
         f = 0.d0
      else if (iph == LIQUID_PHASE) then
         Sg = S(GAS_PHASE)
         Pc_cst = 2.d5
         Sg0 = 1.d0 - 1.d-2
         A = -Pc_cst*log(1.d0 - Sg0) - Pc_cst/(1.d0 - Sg0)*Sg0
         if (Sg < Sg0) then
            f = -Pc_cst*log(1.d0 - Sg)
            dSf(iph) = Pc_cst/(1.d0 - Sg) ! wrt Sg
         else
            f = Pc_cst*Sg/(1.d0 - Sg0) + A
            dSf(iph) = Pc_cst/(1.d0 - Sg0)  ! wrt Sg
         endif
         ! FIXME: f_PressionCapillaire(LIQUID_PHASE) = - Pc
         f = -f   ! because P(LIQUID_PHASE) = Pref(=Pg) + f_PressionCapillaire(LIQUID_PHASE)
         ! NO modification of sign of dSf because f = -f and dSf(iph) = -dSf(iph) wrt Sl = 1 - Sg
      endif

   end subroutine f_PressionCapillaire

end module CapillaryPressure
