   !> \brief Non linear update of Pw and
   !! determine the mode of the injection well
   !! (flowrate or pressure). The injection
   !! well is monophasic liquid
   !!
   !! As long as the pressure is less or egal to
   !! the pressure max, the flowrate of the well
   !! is imposed. If the pressure is too high,
   !! the flowrate is no more fixed and the pressure
   !! is set as Pressure max.
   subroutine DefFlashWells_NewtonFlashNonLinWellInj

      integer :: nWell
      double precision :: Flowrate_head

      do nWell = 1, NodebyWellInjLocal%Nb

         if (DataWellInjLocal(nWell)%IndWell == 'c') cycle ! well is closed

         if (DataWellInjLocal(nWell)%IndWell == 'f') then ! flowrate mode

            ! non linear update of the unknown pressure in well
            call DefFlashWells_NonLinPressureUpdateWellInj(nWell)

            if (IncPressionWellInj(nWell) > DataWellInjLocal(nWell)%PressionMax) then
               DataWellInjLocal(nWell)%IndWell = 'p' ! change to pressure mode
               IncPressionWellInj(nWell) = DataWellInjLocal(nWell)%PressionMax ! Pw = PwMax
            endif

         else if (DataWellInjLocal(nWell)%IndWell == 'p') then ! pressure mode

            if (IncPressionWellInj(nWell) > DataWellInjLocal(nWell)%PressionMax) then
               IncPressionWellInj(nWell) = DataWellInjLocal(nWell)%PressionMax ! With Newton inc, Pw>Pmax, then change it
            endif

            ! compute the new flowrate at the head node of well nWell
            call LoisThermoHydro_divP_wellinj(nWell) ! update thermo Laws of nodes in well num_Well
            call DefFlashWells_PressureToFlowrateWellInj(nWell, Flowrate_head)

            if (Flowrate_head < DataWellInjLocal(nWell)%ImposedFlowrate) then ! inj well then DataWellInjLocal(nWell)%flowrate < 0
               DataWellInjLocal(nWell)%IndWell = 'f' ! change to flowrate mode
               ! non linear update of the unknown pressure in well
               call DefFlashWells_NonLinPressureUpdateWellInj(nWell)
            endif

         else
            print *, "Error in Flash Well: no such context"
         end if

      enddo ! well

   end subroutine DefFlashWells_NewtonFlashNonLinWellInj

   !> \brief Non linear update of Pw and
   !! determine the mode of the projection well
   !! (flowrate or pressure).
   !!
   !! As long as the pressure is less or egal to
   !! the pressure max, the flowrate of the well
   !! is imposed. If the pressure is too high,
   !! the flowrate is no more fixed and the pressure
   !! is set as Pressure max.
   subroutine DefFlashWells_NewtonFlashNonLinWellProd

      integer :: nWell
      double precision :: Flowrate_head

      do nWell = 1, NodebyWellProdLocal%Nb

         if (DataWellProdLocal(nWell)%IndWell == 'c') cycle ! well is closed

         if (DataWellProdLocal(nWell)%IndWell == 'f') then ! flowrate mode

            ! non linear update of the unknown pressure in well
            call DefFlashWells_NonLinPressureUpdateWellProd(nWell)
            if (IncPressionWellProd(nWell) < DataWellProdLocal(nWell)%PressionMin) then
               DataWellProdLocal(nWell)%IndWell = 'p' ! change to pressure mode
               IncPressionWellProd(nWell) = DataWellProdLocal(nWell)%PressionMin ! Pw = PwMin
            endif

         else if (DataWellProdLocal(nWell)%IndWell == 'p') then ! pressure mode

            if (IncPressionWellProd(nWell) < DataWellProdLocal(nWell)%PressionMin) then
               IncPressionWellProd(nWell) = DataWellProdLocal(nWell)%PressionMin ! With Newton inc, Pw<Pmin, then change it
            endif
            ! compute the new flowrate at the head node of well nWell
            call DefFlashWells_PressureToFlowrateWellProd(nWell, Flowrate_head)
            if (Flowrate_head > DataWellProdLocal(nWell)%ImposedFlowrate) then ! Prod well then DataWellProdLocal(nWell)%ImposedFlowrate > 0
               DataWellProdLocal(nWell)%IndWell = 'f' ! change to flowrate mode
               ! non linear update of the unknown pressure in well
               call DefFlashWells_NonLinPressureUpdateWellProd(nWell)
            endif

         else
            print *, "Error in Flash Well: no such context"
         end if

      enddo ! well

   end subroutine DefFlashWells_NewtonFlashNonLinWellProd

   !> \brief Nonlinear update the unknown of injection well (Pressure)
   !! by computing Pw such as the well is at equilibrium wit the others
   !! unknows (Qmol_w, matrix,...)
   subroutine DefFlashWells_NonLinPressureUpdateWellInj(nWell)
      integer, intent(in) :: nWell ! numero of the injection well

      double precision :: Pws, Tws, Sw(NbPhase), Cw(NbComp)
      double precision :: Viscosity, DensiteMolaire, PermRel
      double precision :: dPf, dTf, dCf(NbComp), dSf(NbPhase) ! not used for now, empty passed to f_DensiteMolaire
      double precision :: WIDws, SumMob, SumMobR
      integer :: s, j, nums
      integer :: rt(IndThermique + 1)

      Mob(:, :) = 0.d0
      R(:, :) = 0.d0
      Flow(:) = 0.d0

      do s = NodebyWellInjLocal%Pt(nWell) + 1, NodebyWellInjLocal%Pt(nWell + 1)
         WIDws = NodeDatabyWellInjLocal%Val(s)%WID
         nums = NodebyWellInjLocal%Num(s)
         Pws = IncPressionWellInj(nWell) + PerfoWellInj(s)%PressureDrop ! P_{w,s} = P_w^{n} + PressureDrop_{w,s}^{n-1}
         Tws = PerfoWellInj(s)%Temperature ! T_{w,s}
         R(LIQUID_PHASE, s) = IncNode(nums)%Pression - PerfoWellInj(s)%PressureDrop ! R_s = P_s^{n} - PressureDrop_{w,s}^{n-1}
         Sw(:) = 0.d0
         Sw(LIQUID_PHASE) = 1.d0 ! Monophasic liquid
         Cw = DataWellInjLocal(nWell)%CompTotal
         rt = NodeRocktype(:, s)

         ! LIQUID_PHASE (monophasic in injection well)
         ! Molar density
         call f_DensiteMolaire(LIQUID_PHASE, Pws, Tws, Cw, Sw, &
                               DensiteMolaire, dPf, dTf, dCf, dSf)
         ! viscosity
         call f_Viscosite(LIQUID_PHASE, Pws, Tws, Cw, Sw, &
                          Viscosity, dPf, dTf, dCf, dSf)
         ! Permrel
         call f_PermRel(rt, LIQUID_PHASE, Sw, PermRel, dSf)

         Mob(LIQUID_PHASE, s) = PermRel*DensiteMolaire/Viscosity*WIDws
         ! initialization of RSorted before calling QuickSortCSR
         RSortedInj%Val(s) = R(LIQUID_PHASE, s)
         RSortedInj%Num(s) = s
      enddo ! node s

      write (*, *) 'before sort RSortedInj', RSortedInj%Val(NodebyWellInjLocal%Pt(nWell) + 1:NodebyWellInjLocal%Pt(nWell + 1))

      ! sort CSR R in increasing order of R (carreful, Ri can be egual to Ri+1)
      call QuickSortCSR(RSortedInj, RSortedInj%Pt(nWell) + 1, RSortedInj%Pt(nWell + 1), 'i')
      ! sort CSR Mob in same order than RSortedInj
      do s = RSortedInj%Pt(nWell) + 1, RSortedInj%Pt(nWell + 1)
         MobSortedInj(s) = Mob(LIQUID_PHASE, RSortedInj%Num(s))
      enddo

      ! compute Flow(i) = sum_{j=1}^{i-1} MobSortedInj(j) * (RSortedInj%Val(i) - RSortedInj%Val(j))
      ! Flow(i)<=Flow(i+1) due to the order of MobSortedInj and Rsorted
      do s = NodebyWellInjLocal%Pt(nWell) + 1, NodebyWellInjLocal%Pt(nWell + 1)
         do j = NodebyWellInjLocal%Pt(nWell) + 1, s - 1
            Flow(s) = Flow(s) + MobSortedInj(j)*(RSortedInj%Val(s) - RSortedInj%Val(j))
         enddo
      enddo

      write (*, *) 'RSortedInj', RSortedInj%Val(NodebyWellInjLocal%Pt(nWell) + 1:NodebyWellInjLocal%Pt(nWell + 1))
      write (*, *) 'Flow', Flow(NodebyWellInjLocal%Pt(nWell) + 1:NodebyWellInjLocal%Pt(nWell + 1))

      ! if Flow(n) < -Qmol_w then Pw > RSortedInj(n)
      ! -Qmol_w because flowrate is negatif (injection well)
      if (Flow(RSortedInj%Pt(nWell + 1)) <= -DataWellInjLocal(nWell)%ImposedFlowrate) then
         j = RSortedInj%Pt(nWell + 1)
         write (*, *) 'Flow(n) ,-qmol', Flow(j), -DataWellInjLocal(nWell)%ImposedFlowrate
      else
         ! find j such that Flow(j) <= -Qmol_w <= Flow(j+1)
         ! then r(i) <= Pw <= r(i+1)
         j = RSortedInj%Pt(nWell) + 1
         do while (Flow(j) <= -DataWellInjLocal(nWell)%ImposedFlowrate)
            j = j + 1
         enddo
         j = j - 1
         write (*, *) 'Flow(j) ,-qmol,Flow(j+1)', Flow(j), -DataWellInjLocal(nWell)%ImposedFlowrate, Flow(j + 1)
      endif

      SumMob = 0.d0
      SumMobR = 0.d0
      do s = NodebyWellInjLocal%Pt(nWell) + 1, j
         SumMob = SumMob + MobSortedInj(s)
         SumMobR = SumMobR + MobSortedInj(s)*RSortedInj%Val(s)
      enddo
      write (*, *) 'before update IncPressionWellInj(nWell)', IncPressionWellInj(nWell)
      IncPressionWellInj(nWell) = (-DataWellInjLocal(nWell)%ImposedFlowrate + SumMobR)/SumMob
      write (*, *) 'after update IncPressionWellInj(nWell)', IncPressionWellInj(nWell)

   end subroutine DefFlashWells_NonLinPressureUpdateWellInj

   !> \brief Nonlinear update the unknown of production well (Pressure)
   !! by computing Pw such as the well is at equilibrium wit the others
   !! unknows (Qmol_w, matrix,...)
   subroutine DefFlashWells_NonLinPressureUpdateWellProd(nWell)
      integer, intent(in) :: nWell ! numero of the production well

      double precision :: Ps, Ts, Sat(NbPhase), C(NbComp, NbPhase)
      double precision :: Viscosity, DensiteMolaire, PermRel
      double precision :: dPf, dTf, dCf(NbComp), dSf(NbPhase) ! not used for now, empty passed to f_DensiteMolaire
      double precision :: WIDws, SumMob, SumMobR
      integer :: s, j, nums, m, mph, comptn
      integer :: rt(IndThermique + 1)

      Mob(:, :) = 0.d0
      R(:, :) = 0.d0
      Flow(:) = 0.d0

      comptn = 0

      do s = NodebyWellProdLocal%Pt(nWell) + 1, NodebyWellProdLocal%Pt(nWell + 1)
         WIDws = NodeDatabyWellProdLocal%Val(s)%WID
         nums = NodebyWellProdLocal%Num(s)
         Ts = IncNode(nums)%Temperature ! Ts: Temperature in matrix
         Sat(:) = IncNode(nums)%Saturation(:) ! Sat in matrix
         rt = NodeRocktype(:, s)
         Ps = IncNode(nums)%Pression ! Ps: Reference Pressure in matrix

         ! loop over alpha in Q_s
         do m = 1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
            comptn = comptn + 1 ! comptn = sum(s) sum(alpha in Q_s) 1

            mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)
            R(mph, s) = Ps - PerfoWellProd(s)%PressureDrop ! R_s,alpha = P_s,alpha^{n} - PressureDrop_{w,s}^{n-1}, does not depend on mph
            C(:, mph) = IncNode(nums)%Comp(:, mph) ! Comp in matrix

            ! Molar density
            call f_DensiteMolaire(mph, Ps, Ts, C(:, mph), Sat, &
                                  DensiteMolaire, dPf, dTf, dCf, dSf)
            ! viscosity
            call f_Viscosite(mph, Ps, Ts, C(:, mph), Sat, &
                             Viscosity, dPf, dTf, dCf, dSf)
            ! Permrel
            call f_PermRel(rt, mph, Sat, PermRel, dSf)

            Mob(mph, s) = PermRel*DensiteMolaire/Viscosity*WIDws

            ! initialization of RSorted and MobSortedTempProd before calling QuickSortCSR
            RSortedProd%Num(RSortedProd%Pt(nWell) + comptn) = NodebyWellProdLocal%Pt(nWell) + comptn
            RSortedProd%Val(RSortedProd%Pt(nWell) + comptn) = R(mph, s)
            MobSortedTempProd(RSortedProd%Pt(nWell) + comptn) = Mob(mph, s)
         enddo

      enddo ! node s

      write (*, *) 'before sort RSortedProd good size', RSortedProd%Val(RSortedProd%Pt(nWell) + 1:RSortedProd%Pt(nWell) + comptn)
      write (*, *) 'before sort RSortedProd longer size', RSortedProd%Val(RSortedProd%Pt(nWell) + 1:RSortedProd%Pt(nWell + 1))

      ! sort CSR R in increasing order of R_s^alpha (carreful, Ri can be egual to Ri+1)
      call QuickSortCSR(RSortedProd, RSortedProd%Pt(nWell) + 1, RSortedProd%Pt(nWell) + comptn, 'i')
      ! sort CSR Mob in same order than RSortedProd
      do s = RSortedProd%Pt(nWell) + 1, RSortedProd%Pt(nWell) + comptn
         MobSortedProd(s) = MobSortedTempProd(RSortedProd%Num(s))
      enddo

      ! compute Flow(i) = sum_{j=i+1}^{n} MobSortedProd(j) * (RSortedProd%Val(j) - RSortedProd%Val(i))
      ! Flow(i+1)<=Flow(i) due to the order of MobSortedProd and Rsorted
      do s = RSortedProd%Pt(nWell) + 1, RSortedProd%Pt(nWell) + comptn
         do j = s + 1, RSortedProd%Pt(nWell) + comptn
            Flow(s) = Flow(s) + MobSortedProd(j)*(RSortedProd%Val(j) - RSortedProd%Val(s))
         enddo
      enddo

      write (*, *) 'RSortedProd', RSortedProd%Val(RSortedProd%Pt(nWell) + 1:RSortedProd%Pt(nWell) + comptn)
      write (*, *) 'Flow', Flow(RSortedProd%Pt(nWell) + 1:RSortedProd%Pt(nWell) + comptn)

      ! if Qmol_w > Flow(1)   then    Pw < RSortedProd(1)
      if (Flow(RSortedProd%Pt(nWell) + 1) <= DataWellProdLocal(nWell)%ImposedFlowrate) then
         j = RSortedProd%Pt(nWell) + 1

         write (*, *) 'Flow(1) ,qmol', Flow(j), DataWellProdLocal(nWell)%ImposedFlowrate

      else
         ! search for the index j such that Flow(j+1)<= Qmol_w <= Flow(j)
         ! then r(i) <= Pw <= r(i+1)
         j = RSortedProd%Pt(nWell) + comptn
         do while (Flow(j) <= DataWellProdLocal(nWell)%ImposedFlowrate)
            j = j - 1
         enddo
         j = j + 1
         write (*, *) 'Flow(j+1) ,qmol,Flow(j)', Flow(j + 1), DataWellProdLocal(nWell)%ImposedFlowrate, Flow(j)
      endif

      SumMob = 0.d0
      SumMobR = 0.d0
      do s = j, RSortedProd%Pt(nWell) + comptn
         SumMob = SumMob + MobSortedProd(s)
         SumMobR = SumMobR + MobSortedProd(s)*RSortedProd%Val(s)
      enddo
      write (*, *) 'before update IncPressionWellProd(nWell)', IncPressionWellProd(nWell)
      IncPressionWellProd(nWell) = (-DataWellProdLocal(nWell)%ImposedFlowrate + SumMobR)/SumMob
      write (*, *) 'after update IncPressionWellProd(nWell)', IncPressionWellProd(nWell)

   end subroutine DefFlashWells_NonLinPressureUpdateWellProd
