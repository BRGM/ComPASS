!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module DefWellFlash

   use IncCVWells, only: PerfoWellProd, NodebyWellProdLocal
   use IncCVReservoir, only: IncNode
   use Thermodynamics, only: f_DensiteMassique
   use DefModel, only: NbPhase, NbComp, LIQUID_PHASE
   use WellState, only: &
      WellState_FlowrateWellProd, WellState_solve_for_temperature, &
      summolarFluxProd, sumnrjFluxProd

   implicit none

   public :: DefWellFlash_TimeFlashWellProd

contains

   !! \brief Execute the flash to determine T and the molar fractions
   !! to update PerfoWellProd(s)%Temperature and PerfoWellProd(s)%Density.
   !! This Flash is performed for a monophasic multicomponent fluid.
   !!
   !! Loop over the nodes s from head to tail to
   !! to determine which phases are present, then computes the temperature
   !! and the mean density.
   !! The well pressure and the pressure drop are taken from
   !! the previous Newton iteration and the previous time step respectively.
   subroutine DefWellFlash_TimeFlashWellProd() &
      bind(C, name="DefWellFlash_TimeFlashWellProd")

      double precision :: T, ResT, Pws, Ci(NbComp), sumni, E
      double precision :: rhoMean
      double precision :: dPf, dTf, dCf(NbComp) ! not used for now, empty passed to f_Enthalpie
      integer :: k, s
      logical :: converged

      call WellState_FlowrateWellProd

      do k = 1, NodebyWellProdLocal%Nb

         ! looping from head to queue
         do s = NodebyWellProdLocal%Pt(k + 1), NodebyWellProdLocal%Pt(k) + 1, -1

            Pws = PerfoWellProd(s)%Pression

            ! Newton method to compute T: R = E-Enthalpie*n = 0
            E = sumnrjFluxProd(s) ! energy

            sumni = sum(summolarFluxProd(:, s))

            ! initialize newton with reservoir temperature
            T = IncNode(NodebyWellProdLocal%Num(s))%Temperature
            call WellState_solve_for_temperature(LIQUID_PHASE, E, Pws, T, sumni, converged, ResT)
            if (.not. converged) then
               print *, "Warning: Newton in DefWellFlash_TimeFlashWellProd has not converged"
               print *, "Residue is", abs(ResT), "Well_idx= ", k, "node_idx= ", s
            end if

            ! update PhysPerfoWell
            PerfoWellProd(s)%Temperature = T
            call f_DensiteMassique(LIQUID_PHASE, Pws, T, Ci, rhoMean, dPf, dTf, dCf)
            PerfoWellProd(s)%Density = rhoMean

         end do ! end of well k
      end do

   end subroutine DefWellFlash_TimeFlashWellProd

end module DefWellFlash
