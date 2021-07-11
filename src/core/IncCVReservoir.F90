!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module IncCVReservoir

   use iso_c_binding, only: c_int, c_double, c_size_t

   use CommonMPI, only: ComPASS_COMM_WORLD, commRank

   use DefModel, only: &
      NbPhase, NbComp, NbContexte, NbEqEquilibreMax, NbIncPTCMax, &
      NbCompThermique, NbIncTotalMax, &
      NumPhasePresente_ctx, NbPhasePresente_ctx, &
      IndThermique, MCP, NbIncTotalPrimMax

   use MeshSchema, only: &
      NbCellOwn_Ncpus, NbFracOwn_Ncpus, NbNodeOwn_Ncpus, &
      NbNodeLocal_Ncpus, NbFracLocal_Ncpus, NbCellLocal_Ncpus, &
      SubArrayInfo, MeshSchema_subarrays_info

   use NumbyContext, only: &
      NbIncPTC_ctx, NbIncTotal_ctx, &
      NumIncPTC2NumIncComp_comp_ctx, NumIncPTC2NumIncComp_phase_ctx, &
      NumCompCtilde_ctx, NbCompCtilde_ctx

   use SchemeParameters, only: &
      NewtonIncreObj_C, NewtonIncreObj_P, NewtonIncreObj_S, NewtonIncreObj_T, &
      eps

   use IncCVReservoirTypes, only: TYPE_IncCVReservoir
   use Thermodynamics, only: f_DensiteMassique

   implicit none

   ! Inc for current time step: cell, fracture faces and nodes
   type(TYPE_IncCVReservoir), allocatable, dimension(:), target, public :: IncAll
   type(TYPE_IncCVReservoir), dimension(:), pointer, public :: &
      IncCell, & !< Cell unknowns for current time step
      IncFrac, & !< Fracture Face unknowns for current time step
      IncNode !< Node unknowns for current time step

   ! Inc for previous time step: current time step - 1
   TYPE(TYPE_IncCVReservoir), allocatable, dimension(:), public :: &
      IncCellPreviousTimeStep, & !< Cell unknowns for previous time step
      IncFracPreviousTimeStep, & !< Fracture Face unknowns for previous time step
      IncNodePreviousTimeStep !< Node unknowns for previous time step

   real(c_double), allocatable, dimension(:, :, :), target, public :: &
      divPhasePressureAll
   real(c_double), dimension(:, :, :), pointer, public :: &
      divPhasePressureCell, divPhasePressureFrac, divPhasePressureNode
   real(c_double), allocatable, dimension(:, :), target, public :: &
      dPhasePressuredSAll
   real(c_double), dimension(:, :), pointer, public :: &
      dPhasePressuredSCell, dPhasePressuredSFrac, dPhasePressuredSNode

   public :: &
      IncCVReservoir_allocate, &
      IncCVReservoir_NewtonIncrement, &
      IncCVReservoir_LoadIncPreviousTimeStep, &
      IncCVReservoir_SaveIncPreviousTimeStep, &
      IncCVReservoir_free, &
      IncCVReservoir_compute_density

   private :: &
      IncCVReservoir_NewtonIncrement_reservoir

contains

#ifndef NDEBUG

   subroutine dump_incv_info() &
      bind(C, name="dump_incv_info")
      integer(c_size_t) :: k, n
      type(SubArrayInfo) :: info

      call MeshSchema_subarrays_info(info)

      n = info%nb%nodes + info%nb%fractures + info%nb%cells
      write (*, *) info%nb%nodes, "nodes"
      write (*, *) info%nb%fractures, "fractures"
      write (*, *) info%nb%cells, "cells"
      write (*, *) "-->", n, "dofs"
      write (*, *) "%% - IncAll"
      do k = 1, n
         write (*, *) "context", k, ":", IncAll(k)%ic
      end do
      write (*, *) "%% - IncNode"
      do k = 1, size(IncNode)
         write (*, *) "context node", k, ":", IncNode(k)%ic
      end do
      write (*, *) "%% - IncCell"
      do k = 1, size(IncCell)
         write (*, *) "context cell", k, ":", IncCell(k)%ic
      end do
      write (*, *) "%% - IncFrac"
      do k = 1, size(IncFrac)
         write (*, *) "context fracture", k, ":", IncFrac(k)%ic
      end do

   end subroutine dump_incv_info

#endif

   !> \brief Allocate unknowns vectors
   subroutine IncCVReservoir_allocate

      type(SubArrayInfo) :: info
      integer(c_size_t) :: k, n
      integer(c_size_t) :: begin_node, end_node
      integer(c_size_t) :: begin_frac, end_frac
      integer(c_size_t) :: begin_cell, end_cell

      call MeshSchema_subarrays_info(info)
      n = info%nb%nodes + info%nb%fractures + info%nb%cells
      begin_node = info%offset%nodes
      end_node = info%offset%nodes - 1 + info%nb%nodes
      begin_frac = info%offset%fractures
      end_frac = info%offset%fractures - 1 + info%nb%fractures
      begin_cell = info%offset%cells
      end_cell = info%offset%cells - 1 + info%nb%cells

      allocate (IncAll(n))
      IncNode => IncAll(begin_node:end_node)
      IncFrac => IncAll(begin_frac:end_frac)
      IncCell => IncAll(begin_cell:end_cell)

      allocate (divPhasePressureAll(NbIncTotalPrimMax, NbPhase, n))
      divPhasePressureNode => divPhasePressureAll(:, :, begin_node:end_node)
      divPhasePressureFrac => divPhasePressureAll(:, :, begin_frac:end_frac)
      divPhasePressureCell => divPhasePressureAll(:, :, begin_cell:end_cell)
      divPhasePressureAll = 0.d0
      do k = 1, n
         divPhasePressureAll(1, :, k) = 1.d0 ! CHECKME: reference pressure - first primary unknown
      end do

      allocate (dPhasePressuredSAll(NbPhase, n))
      dPhasePressuredSNode => dPhasePressuredSAll(:, begin_node:end_node)
      dPhasePressuredSFrac => dPhasePressuredSAll(:, begin_frac:end_frac)
      dPhasePressuredSCell => dPhasePressuredSAll(:, begin_cell:end_cell)
      dPhasePressuredSAll = 0.d0

      allocate (IncCellPreviousTimeStep(NbCellLocal_Ncpus(commRank + 1)))
      allocate (IncFracPreviousTimeStep(NbFracLocal_Ncpus(commRank + 1)))
      allocate (IncNodePreviousTimeStep(NbNodeLocal_Ncpus(commRank + 1)))

   end subroutine IncCVReservoir_allocate

   !> \brief Deallocate unknowns vectors
   subroutine IncCVReservoir_free

      nullify (IncNode, IncFrac, IncCell)
      deallocate (IncAll)

      nullify (divPhasePressureNode)
      nullify (divPhasePressureFrac)
      nullify (divPhasePressureCell)
      deallocate (divPhasePressureAll)

      nullify (dPhasePressuredSNode)
      nullify (dPhasePressuredSFrac)
      nullify (dPhasePressuredSCell)
      deallocate (dPhasePressuredSAll)

      deallocate (IncCellPreviousTimeStep)
      deallocate (IncFracPreviousTimeStep)
      deallocate (IncNodePreviousTimeStep)

   end subroutine IncCVReservoir_free

   !> \brief Loops over nodes, fracs and cells to increment the unknowns
   subroutine IncCVReservoir_NewtonIncrement( &
      NewtonIncreNode, NewtonIncreFrac, NewtonIncreCell, &
      relax)

      double precision, dimension(:, :), intent(in) :: &
         NewtonIncreNode, &
         NewtonIncreFrac, &
         NewtonIncreCell

      double precision, intent(in) :: relax

      integer :: k

      ! nodes
      do k = 1, NbNodeLocal_Ncpus(commRank + 1)
         call IncCVReservoir_NewtonIncrement_reservoir(IncNode(k), NewtonIncreNode(:, k), relax)
      end do

      ! fracture faces
      do k = 1, NbFracLocal_Ncpus(commRank + 1)
         call IncCVReservoir_NewtonIncrement_reservoir(IncFrac(k), NewtonIncreFrac(:, k), relax)
      end do

      ! cells
      do k = 1, NbCellLocal_Ncpus(commRank + 1)
         call IncCVReservoir_NewtonIncrement_reservoir(IncCell(k), NewtonIncreCell(:, k), relax)
         !     write(*,*) ' increment cell ',k,NewtonIncreCell(:,k)
      end do

   end subroutine IncCVReservoir_NewtonIncrement

   !> \brief Realize Newton increment of each control volume           <br>
   !! NOMBERING IS FIXED:                                              <br>
   !!  Pressure=1,                                                     <br>
   !!  Temperature=2,                                                   <br>
   !!  Molar fraction of present component (without Ctilde),           <br>
   !!  Saturation of present phase,                                     <br>
   !!  Molar fraction only present in absent phase: Ctilde (put into inc%AccVol(icp) ???)
   subroutine IncCVReservoir_NewtonIncrement_reservoir(inc, incre, relax)

      type(TYPE_IncCVReservoir), intent(inout) :: inc
      double precision, intent(in) :: incre(NbIncTotalMax), relax

      integer :: i, icp, iph
      integer :: NbIncPTC, NbIncTotal, NbPhasePresente
      integer :: ic

      ic = inc%ic
      NbIncPTC = NbIncPTC_ctx(ic)
      NbIncTotal = NbIncTotal_ctx(ic)
      NbPhasePresente = NbPhasePresente_ctx(ic)

      ! increment Pressure
      inc%Pression = inc%Pression + relax*incre(1)

      !    write(*,*)' increment P ',relax,incre(1)

#ifdef _THERMIQUE_

      ! increment Temperature
      inc%Temperature = inc%Temperature + relax*incre(2)
#endif

      ! increment comp
      do i = 2 + IndThermique, NbIncPTC

         icp = NumIncPTC2NumIncComp_comp_ctx(i, ic)
         iph = NumIncPTC2NumIncComp_phase_ctx(i, ic)
         inc%Comp(icp, iph) = inc%Comp(icp, iph) + relax*incre(i)
      enddo

      ! increment saturation
      do i = 1, NbPhasePresente
         iph = NumPhasePresente_ctx(i, ic)

         inc%Saturation(iph) = inc%Saturation(iph) + relax*incre(iph + NbIncPTC)
      end do

#ifdef _WIP_FREEFLOW_STRUCTURES_
      ! increment freeflow molar flowrate
      if (ic >= 2**NbPhase) then ! FIXME: loop over freeflow dof only, avoid reservoir node
         do i = 1, NbPhasePresente
            iph = NumPhasePresente_ctx(i, ic)

            inc%FreeFlow_flowrate(iph) = inc%FreeFlow_flowrate(iph) + relax*incre(iph + NbIncPTC + NbPhasePresente)
         enddo
      endif
#endif

      ! Contribution to Ctilde (to test appearance of phase)
      do i = 1, NbCompCtilde_ctx(ic)
         icp = NumCompCtilde_ctx(i, ic)
         inc%AccVol(icp) = incre(NbIncTotal + i)
      end do

   end subroutine IncCVReservoir_NewtonIncrement_reservoir

   !> \brief Save current status if it is necessary to start again current time iteration.
   !!
   !! Copy IncObj to IncObjPreviousTimeStep
   subroutine IncCVReservoir_SaveIncPreviousTimeStep

      integer :: k

      ! save current status
      do k = 1, NbNodeLocal_Ncpus(commRank + 1)
         IncNodePreviousTimeStep(k) = IncNode(k)
      end do
      do k = 1, NbFracLocal_Ncpus(commRank + 1)
         IncFracPreviousTimeStep(k) = IncFrac(k)
      end do
      do k = 1, NbCellLocal_Ncpus(commRank + 1)
         IncCellPreviousTimeStep(k) = IncCell(k)
      end do

   end subroutine IncCVReservoir_SaveIncPreviousTimeStep

   !> \brief Load previous status to start again current time iteration.
   !!
   !! Copy IncObjPreviousTimeStep to IncObj
   subroutine IncCVReservoir_LoadIncPreviousTimeStep

      integer :: k

      do k = 1, NbNodeLocal_Ncpus(commRank + 1)
         IncNode(k) = IncNodePreviousTimeStep(k)
      end do
      do k = 1, NbFracLocal_Ncpus(commRank + 1)
         IncFrac(k) = IncFracPreviousTimeStep(k)
      end do
      do k = 1, NbCellLocal_Ncpus(commRank + 1)
         IncCell(k) = IncCellPreviousTimeStep(k)
      end do

   end subroutine IncCVReservoir_LoadIncPreviousTimeStep

   !> \brief Transform IncCVReservoir format to a vector,
   !! necessary for visualization.
   !!
   !! Transform Inc into output vector datavisu
   !! The structure of the vector is                                                 <br>
   !!   (Pressure, Temperature,                                                      <br>
   !!        Comp(1), ... , Comp(n),                                                 <br>
   !!            Saturation(1), ..., Saturation(n))
   SUBROUTINE IncCVReservoir_ToVec_cv(NbIncOwn, Inc, datavisu)

      INTEGER, INTENT(IN) :: NbIncOwn
      TYPE(Type_IncCVReservoir), DIMENSION(:), INTENT(IN) :: Inc
      DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: datavisu

      integer :: i_data, j_data
      integer :: i, j

      ! Pressure
      i_data = 1
      j_data = NbIncOwn

      datavisu(i_data:j_data) = Inc(1:NbIncOwn)%Pression

      ! Temperature
#ifdef _THERMIQUE_
      i_data = j_data + 1
      j_data = j_data + NbIncOwn
      datavisu(i_data:j_data) = Inc(1:NbIncOwn)%Temperature
#endif

      ! Comp
      DO j = 1, NbPhase
         DO i = 1, NbComp

            IF (MCP(i, j) == 1) THEN
               i_data = j_data + 1
               j_data = j_data + NbIncOwn
               datavisu(i_data:j_data) = Inc(1:NbIncOwn)%Comp(i, j)
            ENDIF
         ENDDO
      ENDDO

      ! Saturation
      DO i = 1, NbPhase
         i_data = j_data + 1
         j_data = j_data + NbIncOwn
         datavisu(i_data:j_data) = Inc(1:NbIncOwn)%Saturation(i)
      ENDDO
   ENDSUBROUTINE IncCVReservoir_ToVec_cv

   SUBROUTINE IncCVReservoir_ToVec( &
      datavisucell, &
      datavisufrac, &
      datavisunode, &
      datavisuwellinj, &
      datavisuwellprod)

      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: datavisucell
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: datavisufrac
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: datavisunode
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: datavisuwellinj
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: datavisuwellprod

      INTEGER :: NbCellOwn, NbFracOwn, NbNodeOwn

      NbCellOwn = NbCellOwn_Ncpus(commRank + 1)
      NbFracOwn = NbFracOwn_Ncpus(commRank + 1)
      NbNodeOwn = NbNodeOwn_Ncpus(commRank + 1)

      CALL IncCVReservoir_ToVec_cv(NbCellOwn, IncCell, datavisucell)
      CALL IncCVReservoir_ToVec_cv(NbFracOwn, IncFrac, datavisufrac)
      CALL IncCVReservoir_ToVec_cv(NbNodeOwn, IncNode, datavisunode)

      datavisuwellinj = 1.d0 ! not implemented
      datavisuwellprod = 1.d0 ! not implemented
   ENDSUBROUTINE IncCVReservoir_ToVec

   function IncCVReservoir_compute_density(inc) result(rho)
      type(TYPE_IncCVReservoir), intent(in) :: inc
      real(c_double) :: rho
      integer :: m, mph
      real(c_double) :: rhoph, drhodp, drhodT, drhodC(NbComp)

      rho = 0.d0
      do m = 1, NbPhasePresente_ctx(inc%ic)
         mph = NumPhasePresente_ctx(m, inc%ic)
         call f_DensiteMassique(mph, inc%phase_pressure(mph), inc%Temperature, &
                                inc%Comp(:, mph), rhoph, drhodp, drhodT, drhodC)
         rho = rho + rhoph*inc%Saturation(mph)
      end do

   end function IncCVReservoir_compute_density

end module IncCVReservoir
