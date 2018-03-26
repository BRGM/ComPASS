!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

! Model: 2 phase 1 comp, thermal well

!> \brief Define the flash to determine the phases
!! which are actualy present in each cell, and
!! the mode of the well (flowrate or pressure).
module DefFlash

   use IncCV
   use Physics
   use VAGFrac ! to have rocktypes

   implicit none

   public :: &
      DefFlash_Flash ! Flash after each Newton iteration

   private :: &
      DefFlash_Flash_cv

contains

   !> \brief Main surboutine, after each Newton iteration
   !! execute the flash to determine the phases
   !! which are actualy present, and
   !! the mode of the well (flowrate or pressure).
   subroutine DefFlash_Flash

      integer :: k

      do k = 1, NbNodeLocal_Ncpus(commRank + 1)
         call DefFlash_Flash_cv(IncNode(k), NodeRocktypeLocal(:, k), PoroVolDarcyNode(k))
      end do

      do k = 1, NbFracLocal_Ncpus(commRank + 1)
         call DefFlash_Flash_cv(IncFrac(k), FracRocktypeLocal(:, k), PoroVolDarcyFrac(k))
      end do

      do k = 1, NbCellLocal_Ncpus(commRank + 1)
         call DefFlash_Flash_cv(IncCell(k), CellRocktypeLocal(:, k), PoroVolDarcyCell(k))
      end do

      ! choose between linear or non-linear update of the Newton unknown Pw
      ! The next subroutines also compute the mode of the wells ('pressure' or 'flowrate')
      call DefFlashWells_NewtonFlashLinWells

   end subroutine DefFlash_Flash

   !> \brief Determine the phases
   !! which are actualy present.
   !!
   !! Applied to IncNode, IncFrac and IncCell.
   !! \param[in]      porovol   porous Volume ?????
   !! \param[inout]   inc       Unknown (IncNode, IncFrac or IncCell)
   subroutine DefFlash_Flash_cv(inc, rt, porovol)

      type(Type_IncCV), intent(inout) :: inc
      INTEGER, INTENT(IN) :: rt(IndThermique + 1)
      double precision, intent(in) :: porovol ! porovol

      integer :: i, iph, j, icp, m, mph, ic
      double precision :: DensiteMolaire(NbComp), acc1, acc2, &
         dPf, dTf, dCf(NbComp), dSf(NbPhase)

      double precision :: T, Ha
      double precision :: RZetal
      double precision :: Psat, dTSat
      double precision :: S(NbPhase), Pc, DSPc(NbPhase)
      double precision :: Cal, Cel
      double precision :: PgCag, PgCeg
      double precision :: Pg
      double precision :: Cag
      double precision :: Slrk

      ic = inc%ic
      Pg = inc%Pression
      T = inc%Temperature
      S = inc%Saturation

      IF (rt(1) == 1) THEN
         Slrk = 0.4d0
      ELSEIF (rt(1) == 2) THEN
         Slrk = 0.01d0
      ELSE
         PRINT *, 'error'
         STOP
      ENDIF

      iph = 2
      CALL f_PressionCapillaire(rt, iph, S, Pc, DSPc)

      !   write(*,*)' S Pg Pl ',ic,S,Pg,Pg+Pc

      IF (ic == 2) THEN
         CALL air_henry(T, Ha)
         PgCag = inc%Comp(1, PHASE_WATER)*Ha

         RZetal = 8.314d0*1000.d0/0.018d0
         CALL DefModel_Psat(T, Psat, dTSat)

         iph = 2
         CALL f_PressionCapillaire(rt, iph, S, Pc, DSPc)

         PgCeg = inc%Comp(2, PHASE_WATER)*Psat*DEXP(Pc/(T*RZetal))

         IF (PgCag + PgCeg > Pg) THEN

!        write(*,*)' apparition gas ', Pg, T

            inc%ic = 3
            IF (Pg < PgCeg) THEN
               inc%Pression = PgCeg
            ENDIF
            inc%Saturation(PHASE_GAS) = 0
            inc%Saturation(PHASE_WATER) = 1
            inc%Comp(1, PHASE_GAS) = 0.d0
            inc%Comp(2, PHASE_GAS) = 1.d0

         ENDIF

      ELSE IF (ic == 3) THEN

         IF (S(PHASE_GAS) < 0.d0) THEN

!        write(*,*)' disp du gaz ', Pg, T

            inc%ic = 2
            inc%Saturation(PHASE_GAS) = 0
            inc%Saturation(PHASE_WATER) = 1

         ELSE IF (S(PHASE_WATER) < Slrk) THEN

            write (*, *) ' slrk ', Pg, T

            inc%Saturation(PHASE_GAS) = 1.d0 - 1.d-12-Slrk
            inc%Saturation(PHASE_WATER) = Slrk + 1.d-12
         ENDIF

         Cag = MIN(MAX(inc%Comp(1, PHASE_GAS), 0.d0), 1.d0)
         inc%Comp(1, PHASE_GAS) = Cag
         inc%Comp(2, PHASE_GAS) = 1.d0 - Cag

         Cal = MIN(MAX(inc%Comp(1, PHASE_WATER), 0.d0), 1.d0)
         inc%Comp(1, PHASE_WATER) = Cal
         inc%Comp(2, PHASE_WATER) = 1.d0 - Cal

      ELSE
         PRINT *, "Error in Flash: no such context"
         PRINT *, "only gas in porous medium"
         STOP
      END IF

   end subroutine DefFlash_Flash_cv

end module DefFlash
