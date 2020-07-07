!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module DefWellFlash

   use CommonMPI, only: CommonMPI_abort
   use IncCVWells, only: PerfoWellProd, NodebyWellProdLocal
   use IncCVReservoir, only: IncNode
   use Thermodynamics, only: &
      f_DensiteMassique, f_DensiteMolaire, f_Enthalpie, FluidThermodynamics_Tsat
   use DefModel, only: NbPhase, NbComp, LIQUID_PHASE, GAS_PHASE
   use MeshSchema, only: DataWellProdLocal
   use WellState, only: &
      WellState_FlowrateWellProd, WellState_solve_for_temperature, &
      summolarFluxProd, sumnrjFluxProd

   implicit none

   public :: DefWellFlash_TimeFlashWellProd

contains

   !! \brief Execute the flash to determine T and the molar fractions
   !! to update PerfoWellProd(s)%Temperature and PerfoWellProd(s)%Density.
   !! This Flash is performed for a diphasique monocomponent fluid (H20).
   !!
   !! Loop over the nodes s from head to tail to
   !! to determine which phases are present, then computes the temperature
   !! and the mean density.
   !! The well pressure and the pressure drop are taken from
   !! the previous Newton iteration and the previous time step respectively.
   subroutine DefWellFlash_TimeFlashWellProd

      double precision :: Temp, Pws, liq_molarfrac, ResT
      double precision :: Hgas, Hliq, rhogas, rholiq, avg_molar_dens
      double precision :: sumni, E
      ! not used, empty passed to f_Enthalpie
      double precision :: dPf, dTf, dP_Tsat, sat(NbPhase), molarFrac(NbComp), dCf(NbComp), dSf(NbPhase)
      integer :: nWell, s, ID_PHASE ! ID_PHASE=(-1 if diphasique, GAS_PHASE if gas, LIQUID_PHASE if liq)
      logical :: converged

#ifndef _THERMIQUE_
      call CommonMPI_abort("DefWellFlash_TimeFlashWellProd: thermal transfer must be activated")
#endif

      ! compute flowrate of well nWell (fill summolarFluxProd and sumnrjFluxProd)
      call WellState_FlowrateWellProd

      ! loop over production well
      do nWell = 1, NodebyWellProdLocal%Nb

         ! looping from head to queue
         do s = NodebyWellProdLocal%Pt(nWell + 1), NodebyWellProdLocal%Pt(nWell) + 1, -1

            Pws = PerfoWellProd(s)%Pression
            E = sumnrjFluxProd(s) ! energy
            sumni = sum(summolarFluxProd(:, s))
            if (sumni < 1.0D-12) then
               cycle !keep everything as the previous  timestep
            end if

            if (DataWellProdLocal(nWell)%IndWell == 'c') then ! well is closed,perform ....

               Temp = IncNode(NodebyWellProdLocal%Num(s))%Temperature
               call WellState_solve_for_temperature(LIQUID_PHASE, E, Pws, Temp, sumni, converged, ResT)
               if (.not. converged) then
                  print *, "Warning: Newton in DefWellFlash_TimeFlashWellProd has not converged"
                  print *, "Residue is", abs(ResT), "Well_idx= ", nWell, "node_idx= ", s
               end if

               ! update PhysPerfoWell
               PerfoWellProd(s)%Temperature = Temp
               PerfoWellProd(s)%Saturation(GAS_PHASE) = 0.d0
               PerfoWellProd(s)%Saturation(LIQUID_PHASE) = 1.d0
               Sat(:) = 0.d0
               Sat(LIQUID_PHASE) = 1.d0
               call f_DensiteMassique(LIQUID_PHASE, Pws, Temp, molarFrac, Sat, rholiq, &
                                      dPf, dTf, dCf, dSf)
               PerfoWellProd(s)%Density = rholiq

            else

               ! suppose that the two phases are present at the node
               ID_PHASE = -1
               ! set temperature to saturation temperature at Pws
               call FluidThermodynamics_Tsat(Pws, Temp, dP_Tsat)
               ! thus compute liq_molarfrac thanks to the energy, and the enthalpies
               ! molarFrac is not used in the computation of the enthalpies
               call f_Enthalpie(GAS_PHASE, Pws, Temp, molarFrac, sat, Hgas, dPf, dTf, dCf, dSf)
               call f_Enthalpie(LIQUID_PHASE, Pws, Temp, molarFrac, sat, Hliq, dPf, dTf, dCf, dSf)
               !! and compute liq_molarfrac
               liq_molarfrac = (Hgas - E/sumni)/(Hgas - Hliq)

               if (liq_molarfrac < 0.d0) then ! the hypothesis that the two phases are present is wrong: only gas
                  liq_molarfrac = 0.d0
                  ID_PHASE = GAS_PHASE
                  PerfoWellProd(s)%Saturation(GAS_PHASE) = 1.d0
                  PerfoWellProd(s)%Saturation(LIQUID_PHASE) = 0.d0
               else if (liq_molarfrac > 1.d0) then ! the hypothesis that the two phases are present is wrong: only liquid
                  liq_molarfrac = 1.d0
                  ID_PHASE = LIQUID_PHASE
                  PerfoWellProd(s)%Saturation(GAS_PHASE) = 0.d0
                  PerfoWellProd(s)%Saturation(LIQUID_PHASE) = 1.d0
               endif

               if (ID_PHASE > 0) then
                  !solve temperature
                  ! initialize newton with reservoir temperature
                  Temp = IncNode(NodebyWellProdLocal%Num(s))%Temperature

                  call WellState_solve_for_temperature(ID_PHASE, E, Pws, Temp, sumni, converged, ResT)
                  if (.not. converged) then
                     print *, "Warning: Newton in DefWellFlash_TimeFlashWellProd has not converged"
                     print *, "Residue is", abs(ResT), "Well_idx= ", nWell, "node_idx= ", s
                  end if
               end if

               ! we deduce the mean density
               ! molarFrac is not used in the computation of the massique densities
               call f_DensiteMassique(GAS_PHASE, Pws, Temp, molarFrac, sat, rhogas, dPf, dTf, dCf, dSf)
               call f_DensiteMassique(LIQUID_PHASE, Pws, Temp, molarFrac, sat, rholiq, dPf, dTf, dCf, dSf)
               PerfoWellProd(s)%Density = liq_molarfrac*rholiq + (1.d0 - liq_molarfrac)*rhogas

               !Compute Saturations
               call f_DensiteMolaire(GAS_PHASE, Pws, Temp, molarFrac, sat, rhogas, dPf, dTf, dCf, dSf)
               call f_DensiteMolaire(LIQUID_PHASE, Pws, Temp, molarFrac, sat, rholiq, dPf, dTf, dCf, dSf)
               avg_molar_dens = liq_molarfrac*rholiq + (1.d0 - liq_molarfrac)*rhogas

               PerfoWellProd(s)%Saturation(GAS_PHASE) = (1.d0 - liq_molarfrac)*rhogas/avg_molar_dens
               PerfoWellProd(s)%Saturation(LIQUID_PHASE) = liq_molarfrac*rholiq/avg_molar_dens

               ! fill PhysPerfoWell%T
               PerfoWellProd(s)%Temperature = Temp

            end if !well is cloed

         enddo ! node s

      enddo ! nWell

   end subroutine DefWellFlash_TimeFlashWellProd

end module DefWellFlash
