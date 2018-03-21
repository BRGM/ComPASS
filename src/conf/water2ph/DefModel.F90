!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

! Model: 2 phase 1 comp thermal, MCP=(1,1)

! 1: Gas
! 2: Water

module DefModel

   use CommonType

   implicit none

   ! ! ****** Model ****** ! !

   integer, parameter :: &
      NbComp = ComPASS_NUMBER_OF_COMPONENTS, &
      NbPhase = ComPASS_NUMBER_OF_PHASES

   integer, parameter :: &
      NbContexte = 2**NbPhase - 1

   ! MCP
   integer, parameter, dimension(NbComp, NbPhase) :: &
      MCP = transpose(reshape( &
                      (/1, 1/), (/NbPhase, NbComp/)))

   ! Phase PHASE_WATER is Liquid; PHASE_GAS is gas
   integer, parameter :: PHASE_GAS = 1
   integer, parameter :: PHASE_WATER = 2

   ! CpRoche
   double precision, parameter :: CpRoche = 800.d0*2000.d0 !< ???

   ! thickness of frac
   double precision, parameter :: Thickness = 1.d0 !< Thickness of the fractures

   ! Thermique
#ifdef _THERMIQUE_
   integer, parameter :: IndThermique = 1
#else
   integer, parameter :: IndThermique = 0
#endif

   ! ! ****** Constants derived from model (do not edit) ****** ! !

   ! Nombre Max d'eq d'equilibre
   !               d'eq de fermeture thermodynamique
   !               d'inc P (T) C
   !               d'inc P (T) C primaires
   integer, parameter :: &
      NbEqEquilibreMax = NbComp*(NbPhase - 1), & !< Max number of balance equations
      NbEqFermetureMax = NbPhase + NbEqEquilibreMax, & !< Max number of closure laws
      NbIncPTCMax = 1 + IndThermique + sum(MCP), &
      NbIncPTCSecondMax = NbEqFermetureMax, &
      NbIncPTCSMax = NbIncPTCMax + NbPhase, &
      NbIncPTCSPrimMax = NbComp + IndThermique, &
      NbCompThermique = NbComp + IndThermique

   ! ! ****** How to choose primary variables ****** ! !

   ! Served in module LoisthermoHydro.F90

   ! pschoice=1: manually
   !     it is necessary to give PTCS Prim and PTC Secd for each context: psprim
   ! pschoice=2: Glouton method
   !     the matrix psprim and pssecd are defined formally for compile
   ! pschoise=3: Gauss method
   !     the matrix psprim and pssecd are defined formally for compile

   integer, parameter :: pschoice = 1

   integer, parameter, dimension(NbIncPTCSPrimMax, NbContexte) :: &
      psprim = reshape((/ &
                       1, 2, & ! ic=1
                       1, 2, & ! ic=2
                       1, 5 & ! ic=3
                       /), (/NbIncPTCSPrimMax, NbContexte/))

   integer, parameter, dimension(NbIncPTCSecondMax, NbContexte) :: &
      pssecd = reshape((/ &
                       3, 4, 0, & ! ic=1
                       3, 4, 0, & ! ic=2
                       2, 3, 4 & ! ic=3
                       /), (/NbIncPTCSecondMax, NbContexte/))

   ! ! ****** Alignment method ****** ! !

   ! Served in module Jacobian.F90

   ! aligmethod=1, manually
   !     it is necessary to give a three-dimension matrix: aligmat
   !     aligmat(:,:,ic) is the alignment matrix for context ic
   !     the index order of aligmat(:,:,ic) is (col,row), it allows us
   !     to define aligmat(:,:,ic) without considering that the matrix
   !     in Fortran is column-major
   !
   ! aligmethod=2, inverse diagnal
   !     it is necessary to define aligmat formally for compile

   integer, parameter :: aligmethod = 2

   double precision, parameter, &
      dimension(NbCompThermique, NbCompThermique, NbContexte) :: &
      aligmat = reshape((/ &
                        1.d0, 0.d0, & ! ic=1
                        0.d0, 1.d0, &
                        !
                        1.d0, 0.d0, & ! ic=2
                        0.d0, 1.d0, &
                        !
                        1.d0, 0.d0, & ! ic=3
                        0.d0, 1.d0 &
                        /), (/NbCompThermique, NbCompThermique, NbContexte/))

contains

   !> \brief User set permeability
   !!
   !! \param[in] NbCellG,NbFracG Global number of cell and fracture face
   !! \param[in] IdCellG it is possible to set different permeability by rocktype
   !! \param[in,out] PermCellG Permeability tensor for each cell
   !! \param[in,out] PermFracG Permeability constant for each fracture face
   subroutine DefModel_SetPerm( &
      NbCellG, CellRocktypeG, NbFracG, FracRocktypeG, &
      PermCellG, PermFracG)

      integer, intent(in) :: NbCellG
      integer, dimension(:, :), intent(in) :: CellRocktypeG
      integer, intent(in) :: NbFracG
      integer, dimension(:, :), intent(in) :: FracRocktypeG
      ! ouptuts:
      double precision, dimension(:, :, :), allocatable, intent(inout) :: &
         PermCellG
      double precision, dimension(:), allocatable, intent(inout) :: &
         PermFracG

      integer :: i

      ! allocate and set values
      allocate (PermCellG(3, 3, NbCellG))
      do i = 1, NbCellG
         PermCellG(:, :, i) = 0.d0
         PermCellG(1, 1, i) = 1.d-14
         PermCellG(2, 2, i) = 1.d-14
         PermCellG(3, 3, i) = 1.d-14
      end do

      allocate (PermFracG(NbFracG))
      PermFracG(:) = 1.d-11

   end subroutine DefModel_SetPerm

   subroutine DefModel_SetPorosite( &
      NbCellG, CellRocktypeG, NbFracG, FracRocktypeG, &
      PorositeCell, PorositeFrac)

      integer, intent(in) :: NbCellG
      integer, dimension(:, :), intent(in) :: CellRocktypeG
      integer, intent(in) :: NbFracG
      integer, dimension(:, :), intent(in) :: FracRocktypeG
      ! ouptuts:
      double precision, dimension(:), allocatable, intent(inout) :: &
         PorositeCell
      double precision, dimension(:), allocatable, intent(inout) :: &
         PorositeFrac

      allocate (PorositeCell(NbCellG))
      allocate (PorositeFrac(NbFracG))

      PorositeCell(:) = 1.d-1
      PorositeFrac(:) = 4.d-1
   end subroutine DefModel_SetPorosite

#ifdef _THERMIQUE_

   subroutine DefModel_SetCondThermique( &
      NbCellG, CellRocktypeG, NbFracG, FracRocktypeG, &
      CondThermalCellG, CondThermalFracG)

      integer, intent(in) :: NbCellG
      integer, dimension(:, :), intent(in) :: CellRocktypeG
      integer, intent(in) :: NbFracG
      integer, dimension(:, :), intent(in) :: FracRocktypeG

      ! ouptuts:
      double precision, dimension(:, :, :), allocatable, intent(inout) :: &
         CondThermalCellG
      double precision, dimension(:), allocatable, intent(inout) :: &
         CondThermalFracG

      integer :: i

      allocate (CondThermalCellG(3, 3, NbCellG))
      do i = 1, NbCellG
         CondThermalCellG(:, :, i) = 0.d0
         CondThermalCellG(1, 1, i) = 2.d0
         CondThermalCellG(2, 2, i) = 2.d0
         CondThermalCellG(3, 3, i) = 2.d0
      end do

      allocate (CondThermalFracG(NbFracG))
      CondThermalFracG(:) = 2.d0

   end subroutine DefModel_SetCondThermique

   SUBROUTINE DefModel_SetThermalSource( &
      NbCell, &
      CellThermalSourceType, &
      NbFrac, &
      FracThermalSourceType, &
      CellThermalSource, &
      FracThermalSource)

      INTEGER, INTENT(IN) :: NbCell
      INTEGER, DIMENSION(:), INTENT(IN) :: CellThermalSourceType
      INTEGER, INTENT(IN) :: NbFrac
      INTEGER, DIMENSION(:), INTENT(IN) :: FracThermalSourceType
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: CellThermalSource
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: FracThermalSource

      ALLOCATE (CellThermalSource(NbCell))
      CellThermalSource = 0.d0

      ALLOCATE (FracThermalSource(NbFrac))
      FracThermalSource = 0.d0
   END SUBROUTINE DefModel_SetThermalSource

#endif

end module DefModel
