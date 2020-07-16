!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

#ifndef _THERMIQUE_
#error DefWellFlash: thermal transfer must be activated
#endif

module DefWellFlash

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

   subroutine DefWellFlash_copy_reservoir_states(wk)
      integer, intent(in) :: wk

      double precision :: Tws, Pws, Sg, Sl
      double precision :: rhog, rhol, rho
      ! not used, empty passed to f_Enthalpie
      double precision :: dPf, dTf, dP_Tsat, sat(NbPhase), molarFrac(NbComp), dCf(NbComp), dSf(NbPhase)
      integer :: s, sr

      do s = NodebyWellProdLocal%Pt(wk) + 1, NodebyWellProdLocal%Pt(wk)
         sr = NodebyWellProdLocal%Num(s) ! reservoir node index
         Pws = IncNode(sr)%Pression
         Tws = IncNode(sr)%Temperature
         Sg = IncNode(sr)%Saturation(GAS_PHASE)
         Sl = IncNode(sr)%Saturation(LIQUID_PHASE)
         PerfoWellProd(s)%Pression = Pws
         PerfoWellProd(s)%Temperature = Tws
         PerfoWellProd(s)%Saturation = IncNode(sr)%Saturation
         rho = 0.d0
         if (Sg > 0) then
            call f_DensiteMassique(GAS_PHASE, Pws, Tws, molarFrac, sat, rhog, dPf, dTf, dCf, dSf)
            rho = rho + Sg*rhog
         end if
         if (Sl > 0) then
            call f_DensiteMassique(LIQUID_PHASE, Pws, Tws, molarFrac, sat, rhol, dPf, dTf, dCf, dSf)
            rho = rho + Sl*rhol
         end if
         PerfoWellProd(s)%Density = rho
      enddo ! loop on perforations s of well wk

   end subroutine DefWellFlash_copy_reservoir_states

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

      double precision :: Tws, Pws, xg, Sg, Sl, ResT
      double precision :: Hgas, Hliq, rhog, rhol, rho
      double precision :: sumni, E
      ! not used, empty passed to f_Enthalpie
      double precision :: dPf, dTf, dP_Tsat, sat(NbPhase), molarFrac(NbComp), dCf(NbComp), dSf(NbPhase)
      integer :: wk, s, sr, ID_PHASE ! ID_PHASE=(-1 if diphasique, GAS_PHASE if gas, LIQUID_PHASE if liq)
      logical :: converged

      ! compute flowrate of well wk (fill summolarFluxProd and sumnrjFluxProd)
      call WellState_FlowrateWellProd

      ! loop over production well
      do wk = 1, NodebyWellProdLocal%Nb

         if (DataWellProdLocal(wk)%IndWell == 'c') then ! well is closed

            call DefWellFlash_copy_reservoir_states(wk)

         else ! well is open

            ! looping from head to queue
            do s = NodebyWellProdLocal%Pt(wk + 1), NodebyWellProdLocal%Pt(wk) + 1, -1

               Pws = PerfoWellProd(s)%Pression
               E = sumnrjFluxProd(s) ! energy
               sumni = sum(summolarFluxProd(:, s))
               if (sumni < 1.0D-12) then
                  cycle !keep everything as the previous  timestep
               end if

               ! suppose that the two phases are present at the node
               ID_PHASE = -1
               ! set temperature to saturation temperature at Pws
               call FluidThermodynamics_Tsat(Pws, Tws, dP_Tsat)
               ! thus compute liq_molarfrac thanks to the energy, and the enthalpies
               ! molarFrac is not used in the computation of the enthalpies
               call f_Enthalpie(GAS_PHASE, Pws, Tws, molarFrac, sat, Hgas, dPf, dTf, dCf, dSf)
               call f_Enthalpie(LIQUID_PHASE, Pws, Tws, molarFrac, sat, Hliq, dPf, dTf, dCf, dSf)

               xg = (E/sumni - Hliq)/(Hgas - Hliq)

               if (xg <= 0.d0) then ! the hypothesis that the two phases are present is wrong: only liquid
                  xg = 0.d0
                  ID_PHASE = LIQUID_PHASE
                  PerfoWellProd(s)%Saturation(GAS_PHASE) = 0.d0
                  PerfoWellProd(s)%Saturation(LIQUID_PHASE) = 1.d0
               else if (xg >= 1.d0) then ! the hypothesis that the two phases are present is wrong: only gas
                  xg = 1.d0
                  ID_PHASE = GAS_PHASE
                  PerfoWellProd(s)%Saturation(GAS_PHASE) = 1.d0
                  PerfoWellProd(s)%Saturation(LIQUID_PHASE) = 0.d0
               endif

               if (ID_PHASE > 0) then ! perforation is monophasic
                  ! compute temperature: initialize newton with reservoir temperature
                  Tws = IncNode(NodebyWellProdLocal%Num(s))%Temperature
                  call WellState_solve_for_temperature(ID_PHASE, E, Pws, Tws, sumni, converged, ResT)
                  if (.not. converged) then
                     print *, "Warning: Newton in DefWellFlash_TimeFlashWellProd has not converged"
                     print *, "Residue is", abs(ResT), "Well_idx= ", wk, "node_idx= ", s
                  end if
                  PerfoWellProd(s)%Temperature = Tws
                  call f_DensiteMassique(LIQUID_PHASE, Pws, Tws, molarFrac, sat, PerfoWellProd(s)%Density, dPf, dTf, dCf, dSf)
               else ! perforation is diphasic
                  PerfoWellProd(s)%Temperature = Tws
                  ! molarFrac is not used in the computation of the massique densities
                  call f_DensiteMassique(GAS_PHASE, Pws, Tws, molarFrac, sat, rhog, dPf, dTf, dCf, dSf)
                  call f_DensiteMassique(LIQUID_PHASE, Pws, Tws, molarFrac, sat, rhol, dPf, dTf, dCf, dSf)
                  Sg = (xg/rhog)/((xg/rhog) + ((1.d0 - xg)/rhol))
                  Sl = 1.d0 - Sg
                  PerfoWellProd(s)%Density = Sl*rhol + Sg*rhog
                  PerfoWellProd(s)%Saturation(GAS_PHASE) = Sg
                  PerfoWellProd(s)%Saturation(LIQUID_PHASE) = Sl
               end if ! perforation is monophasic / diphasic

            enddo ! loop on perforations s of well wk

         end if ! well is closed / open

      enddo ! loop on wells wk

   end subroutine DefWellFlash_TimeFlashWellProd

end module DefWellFlash
