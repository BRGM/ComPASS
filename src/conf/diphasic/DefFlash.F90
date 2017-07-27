! Model: 2 phase 1 comp, thermal well

!> \brief Define the flash to determine the phases
!! which are actualy present in each cell, and 
!! the mode of the well (flowrate or pressure).
module DefFlash

  use DefModel
  use IncCV
  use VAGFrac
  use LoisThermoHydro

  implicit none

  double precision, parameter, private :: Tol = 1.0d-5  !< Newton convergence precision
  integer, parameter, private:: Maxiter = 1000 !< Max number of iteration for newton
  integer, parameter, private:: Nslice = 100   !< Number of discretization for the approximation of the integral
  integer, private :: fdFl  !< debug output

  type(CSRdble) :: ZSortedInj   !< CSR vector storing the nodes of each injection well in Z coordinate order
  type(CSRdble) :: RSortedInj   !< CSR vector storing the R of each node of injection well in R order (R_s = P_s - PressureDrop)

  double precision, allocatable, dimension(:) :: &
       MobSortedInj !< like CSR vector storing the Mob of each node of injection well in R order (R_s = P_s - PressureDrop)

  type(CSRdble) :: RSortedProd   !< CSR vector storing the R of each node of production well in R order (R_s^alpha = P_s^alpha - PressureDrop)

  double precision, allocatable, dimension(:) :: &
       MobSortedTempProd,& !< like CSR vector storing the Mob of each phase of each node
       MobSortedProd       !  of production well in R order (R_s^alpha = P_s^alpha - PressureDrop)

  double precision, allocatable, dimension(:,:) :: R, Mob  ! CSR vectors used to compute the non linear update of the pressure
  double precision, allocatable, dimension(:) :: Flow      ! CSR vectors used to compute the non linear update of the pressure


  ! molar/energy fluxes for production well(s)
  double precision, allocatable, dimension(:,:) :: summolarFluxProd !< Molar flux for production well: summolarFluxProd(component,node_well)
  double precision, allocatable, dimension(:) :: sumnrjFluxProd     !< Energy flux for production well: sumnrjFluxProd(node_well)

  ! molar fluxes for injection well(s), head node
  double precision, allocatable, dimension(:) :: headmolarFluxInj !< Molar flux for injection well: headmolarFluxProd(node_well)

  
  public :: &  
       DefFlash_allocate,  &  ! Allocation, initialization and deallocation (in NN.F90)
       DefFlash_Flash,     &  ! Flash after each Newton iteration
       DefFlash_TimeFlash, &  ! Flash after each time step
       DefFlash_free                          

  private :: &
       DefFlash_Flash_cv

  private :: &
       DefFlash_NewtonFlashLinWellInj,  &    ! Flash after each Newton iteratio
       DefFlash_NewtonFlashLinWellProd, &    ! Flash after each Newton iteration
       DefFlash_NewtonFlashNonLinWellInj, &  ! Flash after each Newton iteration
       DefFlash_NewtonFlashNonLinWellProd, & ! Flash after each Newton iteration
       DefFlash_NonLinPressureUpdateWellInj, &
       DefFlash_NonLinPressureUpdateWellProd

  private :: &
       DefFlash_TimeFlashSinglePhaseWellProd, & ! Flash after time step to compute T and rho
       DefFlash_TimeFlashTwoPhasesProd          ! Flash after time step to compute T and rho

  private :: &
       DefFlash_FlowrateWellProd, &
       DefFlash_PressureToFlowrateWellProd, &
       DefFlash_PressureToFlowrateWellInj,  &
       QuickSortCSR, &                              
       DefFlash_SortHeights_and_Init
  ! DefFlash_PressureDropInj ! Pressure drop calculation allong the injection well

contains


  !> \brief Main surboutine, after each Newton iteration
  !! execute the flash to determine the phases
  !! which are actualy present, and
  !! the mode of the well (flowrate or pressure).
  subroutine DefFlash_Flash

    integer :: k
    double precision :: Psat, dTsat, Tsat
    integer :: rocktype

    rocktype = 1
    do k=1, NbNodeLocal_Ncpus(commRank+1)
       call DefFlash_Flash_cv(rocktype, IncNode(k), PoroVolDarcyNode(k))
    end do

    rocktype = 2
    do k=1, NbFracLocal_Ncpus(commRank+1)
       call DefFlash_Flash_cv(rocktype, IncFrac(k), PoroVolDarcyFrac(k))
    end do

    do k=1, NbCellLocal_Ncpus(commRank+1)
       call DefFlash_Flash_cv(rocktype, IncCell(k), PoroVolDarcyCell(k))
    end do

    ! choose between linear or non-linear update of the Newton unknown Pw
    ! The next subroutines also compute the mode of the well ('pressure' or 'flowrate')
    call DefFlash_NewtonFlashLinWellInj
    call DefFlash_NewtonFlashLinWellProd

  end subroutine DefFlash_Flash

  !> \brief Main surboutine, after each time iteration
  subroutine DefFlash_TimeFlash

    integer :: num_Well
    
    call DefFlash_TimeFlashSinglePhaseWellProd

    ! compute 
    do num_Well=1, NbWellInjLocal_Ncpus(commRank+1)
       call DefFlash_PressureToFlowrateWellInj(num_Well, headmolarFluxInj(num_Well))
       ! print*, "head ", headmolarFluxInj(num_Well)
    end do

  end subroutine DefFlash_TimeFlash


  !> Allocate global vectors used only in this file
  subroutine DefFlash_allocate
    integer :: k, Nnz
    character(len=300) :: fn, pid

    ! put debugging magic here, open a file descriptor per processor
    ! for more readability
#define _DEBUG_LVL1_
#if defined _DEBUG_ && defined _DEBUG_LVL1_
    if (.false.) then
       fdFl = 6 ! stdout: 6, scratch (temprary discarded file): 12
    else
       fdFl = 12
       write(pid, *) getpid()
       write(fn, '(a,a,a)') 'proc.', trim(adjustl(pid)), '.log'
       open(unit=fdFl, file=trim(fn), action='WRITE')
    end if
#else
    fdFl = 12
    open(unit=fdFl, status='SCRATCH')
#endif

    ! allocate flowrate
    allocate(headmolarFluxInj(NbWellInjLocal_Ncpus(commRank+1)))
    
    Nnz = NodebyWellProdLocal%Pt(NodebyWellProdLocal%Nb+1)
    allocate(summolarFluxProd(NbComp, Nnz))
    allocate(sumnrjFluxProd(Nnz))
    
    ZSortedInj%Nb = NodebyWellInjLocal%Nb
    allocate(ZSortedInj%Pt(ZSortedInj%Nb+1))
    ZSortedInj%Pt = NodebyWellInjLocal%Pt
    allocate(ZSortedInj%Num(ZSortedInj%Pt(ZSortedInj%Nb+1)))
    allocate(ZSortedInj%Val(ZSortedInj%Pt(ZSortedInj%Nb+1)))

    RSortedInj%Nb = NodebyWellInjLocal%Nb
    allocate(RSortedInj%Pt(RSortedInj%Nb+1))
    RSortedInj%Pt = NodebyWellInjLocal%Pt
    allocate(RSortedInj%Num(RSortedInj%Pt(RSortedInj%Nb+1)))
    allocate(RSortedInj%Val(RSortedInj%Pt(RSortedInj%Nb+1)))

    ! prod well : multi-phase
    ! RSortedProd will contain the vector r_s^alpha sorted with respect to s AND alpha
    ! then size max is NbPhase * NodebyWellProdLocal
    RSortedProd%Nb = NodebyWellProdLocal%Nb
    allocate(RSortedProd%Pt(RSortedProd%Nb+1))
    RSortedProd%Pt = NbPhase * NodebyWellProdLocal%Pt
    allocate(RSortedProd%Num(RSortedProd%Pt(RSortedProd%Nb+1)))
    allocate(RSortedProd%Val(RSortedProd%Pt(RSortedProd%Nb+1)))

    allocate(MobSortedInj(RSortedInj%Pt(RSortedInj%Nb+1)))
    allocate(MobSortedTempProd(RSortedProd%Pt(RSortedProd%Nb+1)))
    allocate(MobSortedProd(RSortedProd%Pt(RSortedProd%Nb+1)))

    allocate(Mob(NbPhase, NodebyWellInjLocal%Pt(NodebyWellInjLocal%Nb+1)))
    allocate(R(NbPhase, NodebyWellInjLocal%Pt(NodebyWellInjLocal%Nb+1)))
    allocate(Flow(NbPhase * NodebyWellInjLocal%Pt(NodebyWellInjLocal%Nb+1)))

    ! init sort height
    call DefFlash_SortHeights_and_Init

  end subroutine DefFlash_allocate


  subroutine DefFlash_free()

    ! call CommonType_deallocCSRdble(ZSortedInj)
    ! close(unit=12)
    ! call CommonType_deallocCSRdble(ZSortedInj)
    ! call CommonType_deallocCSRdble(RSortedInj)
    ! call CommonType_deallocCSRdble(RSortedProd)
    deallocate(MobSortedInj)
    deallocate(MobSortedTempProd)
    deallocate(MobSortedProd)
    deallocate(Mob, R, Flow)
    deallocate(summolarFluxProd, sumnrjFluxProd)
    deallocate(headmolarFluxInj)

  end subroutine DefFlash_free

  !> \brief Performs some initialization, and sorts the injector well nodes
  subroutine DefFlash_SortHeights_and_Init
    integer :: k, s

    ! init these module global variables
    do k=1, NbWellInjLocal_Ncpus(commRank+1)
       do s=ZSortedInj%Pt(k)+1, ZSortedInj%Pt(k+1)
          ZSortedInj%Num(s) = s
          ZSortedInj%Val(s) = XNodeLocal(3, NodebyWellInjLocal%Num(s))
       end do
    end do

    !! sorting the well nodes according to the heights
    do k=1, NbWellInjLocal_Ncpus(commRank+1)
       call QuickSortCSR(ZSortedInj, ZSortedInj%Pt(k)+1, ZSortedInj%Pt(k+1), 'd')
    end do

    write(fdFl, *) '{'
    write(fdFl, *) ZSortedInj%Num
    write(fdFl, *) ZSortedInj%Val
    write(fdFl, *) '}{'
    do s=1, ZSortedInj%Pt(ZSortedInj%Nb+1)
       write(fdFl, *)  'ptL', ZSortedInj%Num(s), 'numL', NodebyWellInjLocal%Num(ZSortedInj%Num(s)), & 
            ZSortedInj%Val(s)
    end do
    write(fdFl, *) '}'

  end subroutine DefFlash_SortHeights_and_Init


  !> \brief Determine the phases
  !! which are actualy present.
  !!
  !! Applied to IncNode, IncFrac and IncCell.
  !! \param[in]      porovol   porous Volume ?????
  !! \param[inout]   inc       Unknown (IncNode, IncFrac or IncCell)
  subroutine DefFlash_Flash_cv(rocktype, inc, porovol)

    INTEGER, INTENT(IN) :: rocktype
    type(Type_IncCV), intent(inout) :: inc
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

    Slrk = 0.4d0     

    iph = 2
    CALL f_PressionCapillaire(rocktype,iph,S,Pc,DSPc)

 !   write(*,*)' S Pg Pl ',ic,S,Pg,Pg+Pc
    
    IF(ic == 2)THEN
      CALL air_henry(T,Ha)
      PgCag = inc%Comp(1,PHASE_WATER) * Ha

      RZetal = 8.314d0 * 1000.d0 / 0.018d0
      CALL DefModel_Psat(T, Psat, dTSat)

      iph = 2
      CALL f_PressionCapillaire(rocktype,iph,S,Pc,DSPc)

      PgCeg = inc%Comp(2,PHASE_WATER) * Psat * DEXP(Pc/(T*RZetal))

      IF(PgCag + PgCeg > Pg)THEN

        write(*,*)' apparition gas ', Pg, T
        
        inc%ic = 3
        inc%Saturation(PHASE_GAS) = 0.d0
        inc%Saturation(PHASE_WATER) = 1.d0

        if (T<=273.d0) then             
           inc%Temperature = 273.d0 
        endif
        CALL DefModel_Psat(T, Psat, dTSat)

        if (Pg<=1.d-3) then             
           inc%Pression = Psat
        endif        

        inc%Comp(1,PHASE_GAS) = 0.d0 
        inc%Comp(2,PHASE_GAS) = 1.d0
      ENDIF

      
   ELSE IF(ic == 3)THEN
      
      IF(S(PHASE_GAS) < 0.d0)THEN
         
        write(*,*)' disp du gaz ', Pg, T
         
        inc%ic = 2
        inc%Saturation(PHASE_GAS) = 0
        inc%Saturation(PHASE_WATER) = 1.d0
        
      ELSE IF(S(PHASE_WATER) < Slrk)THEN
        inc%Saturation(PHASE_GAS) = 1.d0 - 1.d-12 - Slrk
        inc%Saturation(PHASE_WATER) = Slrk+1.d-12
     ENDIF
     
      Cag = MIN(MAX(inc%Comp(1,PHASE_GAS),0.d0),1.d0)
      inc%Comp(1,PHASE_GAS) = Cag
      inc%Comp(2,PHASE_GAS) = 1.d0 - Cag

      Cal = MIN(MAX(inc%Comp(1,PHASE_WATER),0.d0),1.d0)
      inc%Comp(1,PHASE_WATER) = Cal
      inc%Comp(2,PHASE_WATER) = 1.d0 - Cal

        if (T<=273.d0) then             
           inc%Temperature = 273.d0 
        endif
        CALL DefModel_Psat(T, Psat, dTSat)


        if (Pg<=1.d-3) then             
           inc%Pression = Psat
        endif     
          

    ELSE
      PRINT*, "Error in Flash: no such context"
      PRINT*, "only gas in porous medium"
      STOP
    END IF

  end subroutine DefFlash_Flash_cv

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
  subroutine DefFlash_NewtonFlashNonLinWellInj

    integer :: nWell
    double precision :: Flowrate_head

    do nWell=1, NodebyWellInjLocal%Nb

       if(DataWellInjLocal(nWell)%IndWell == 'f') then ! flowrate mode

          ! non linear update of the unknown pressure in well
          call DefFlash_NonLinPressureUpdateWellInj(nWell)

          if(IncPressionWellInj(nWell) > DataWellInjLocal(nWell)%PressionMax) then
             DataWellInjLocal(nWell)%IndWell = 'p'  ! change to pressure mode
             IncPressionWellInj(nWell) = DataWellInjLocal(nWell)%PressionMax  ! Pw = PwMax
          endif

       else if(DataWellInjLocal(nWell)%IndWell == 'p') then ! pressure mode

          if(IncPressionWellInj(nWell) > DataWellInjLocal(nWell)%PressionMax) then
             IncPressionWellInj(nWell) = DataWellInjLocal(nWell)%PressionMax  ! With Newton inc, Pw>Pmax, then change it
          endif

          ! compute the new flowrate at the head node of well nWell
          call LoisThermoHydro_divP_wellinj(nWell) ! update thermo Laws of nodes in well num_Well
          call DefFlash_PressureToFlowrateWellInj(nWell, Flowrate_head)

          if (Flowrate_head < DataWellInjLocal(nWell)%FlowrateImposed) then  ! inj well then DataWellInjLocal(nWell)%flowrate < 0
             DataWellInjLocal(nWell)%IndWell = 'f'  ! change to flowrate mode
             ! non linear update of the unknown pressure in well
             call DefFlash_NonLinPressureUpdateWellInj(nWell)
          endif

       else
          print*, "Error in Flash Well: no such context"
       end if

    enddo ! well

  end subroutine DefFlash_NewtonFlashNonLinWellInj

  !> \brief Nonlinear update the unknown of injection well (Pressure)
  !! by computing Pw such as the well is at equilibrium wit the others
  !! unknows (Qmol_w, matrix,...)
  subroutine DefFlash_NonLinPressureUpdateWellInj(nWell)
    integer, intent(in) :: nWell ! numero of the injection well

    double precision :: Pws, Tws, Sw(NbPhase), Cw(NbComp)
    double precision :: Viscosity, DensiteMolaire, PermRel
    double precision :: dPf, dTf, dCf(NbComp), dSf(NbPhase) ! not used for now, empty passed to f_DensiteMolaire
    double precision :: WIDws, SumMob, SumMobR
    integer :: s, j, nums

    Mob(:,:) = 0.d0
    R(:,:) = 0.d0
    Flow(:) = 0.d0

    do s=NodebyWellInjLocal%Pt(nWell)+1, NodebyWellInjLocal%Pt(nWell+1)
       WIDws = NodeDatabyWellInjLocal%Val(s)%WID
       nums = NodebyWellInjLocal%Num(s)
       Pws = IncPressionWellInj(nWell) + PerfoWellInj(s)%PressureDrop   ! P_{w,s} = P_w^{n} + PressureDrop_{w,s}^{n-1}
       Tws = PerfoWellInj(s)%Temperature                                ! T_{w,s}
       R(PHASE_WATER,s) = IncNode(nums)%Pression - PerfoWellInj(s)%PressureDrop       ! R_s = P_s^{n} - PressureDrop_{w,s}^{n-1}
       Sw(PHASE_WATER) = 1.d0                                           ! Monophasic liquid
!       Sw(PHASE_GAS) = 0.d0
       Cw = DataWellInjLocal(nWell)%CompTotal

       ! PHASE_WATER (monophasic in injection well)
       ! Molar density
       call f_DensiteMolaire(PHASE_WATER, Pws, Tws, Cw, Sw, &
            DensiteMolaire, dPf, dTf, dCf, dSf)
       ! viscosity
       call f_Viscosite(PHASE_WATER, Pws, Tws, Cw, Sw, &
            Viscosity, dPf, dTf, dCf, dSf)
       ! Permrel
       call f_PermRel(PHASE_WATER, Sw, PermRel, dSf)

       Mob(PHASE_WATER, s) = PermRel * DensiteMolaire / Viscosity * WIDws
       ! initialization of RSorted before calling QuickSortCSR
       RSortedInj%Val(s) = R(PHASE_WATER,s)
       RSortedInj%Num(s) = s
    enddo ! node s

    write(*,*)'before sort RSortedInj',RSortedInj%Val(NodebyWellInjLocal%Pt(nWell)+1: NodebyWellInjLocal%Pt(nWell+1))


    ! sort CSR R in increasing order of R (carreful, Ri can be egual to Ri+1)
    call QuickSortCSR(RSortedInj, RSortedInj%Pt(nWell)+1, RSortedInj%Pt(nWell+1), 'i')
    ! sort CSR Mob in same order than RSortedInj
    do s=RSortedInj%Pt(nWell)+1, RSortedInj%Pt(nWell+1)
       MobSortedInj(s) = Mob(PHASE_WATER, RSortedInj%Num(s))
    enddo

    ! compute Flow(i) = sum_{j=1}^{i-1} MobSortedInj(j) * (RSortedInj%Val(i) - RSortedInj%Val(j))
    ! Flow(i)<=Flow(i+1) due to the order of MobSortedInj and Rsorted
    do s=NodebyWellInjLocal%Pt(nWell)+1, NodebyWellInjLocal%Pt(nWell+1)
       do j=NodebyWellInjLocal%Pt(nWell)+1,s-1
          Flow(s) = Flow(s) + MobSortedInj(j) * (RSortedInj%Val(s)-RSortedInj%Val(j))
       enddo
    enddo


    write(*,*)'RSortedInj',RSortedInj%Val(NodebyWellInjLocal%Pt(nWell)+1: NodebyWellInjLocal%Pt(nWell+1))
    write(*,*)'Flow',Flow(NodebyWellInjLocal%Pt(nWell)+1: NodebyWellInjLocal%Pt(nWell+1))


    ! if Flow(n) < -Qmol_w then Pw > RSortedInj(n)
    ! -Qmol_w because flowrate is negatif (injection well)
    if(Flow(RSortedInj%Pt(nWell+1)) <= -DataWellInjLocal(nWell)%FlowrateImposed) then
       j=RSortedInj%Pt(nWell+1)
       write(*,*)'Flow(n) ,-qmol',Flow(j), -DataWellInjLocal(nWell)%FlowrateImposed
    else
       ! find j such that Flow(j) <= -Qmol_w <= Flow(j+1)
       ! then r(i) <= Pw <= r(i+1)
       j=RSortedInj%Pt(nWell)+1
       do while(Flow(j) <= -DataWellInjLocal(nWell)%FlowrateImposed)
          j=j+1
       enddo
       j=j-1
       write(*,*)'Flow(j) ,-qmol,Flow(j+1)',Flow(j), -DataWellInjLocal(nWell)%FlowrateImposed, Flow(j+1)
    endif

    SumMob = 0.d0
    SumMobR = 0.d0
    do s=NodebyWellInjLocal%Pt(nWell)+1,j
       SumMob = SumMob + MobSortedInj(s)
       SumMobR = SumMobR + MobSortedInj(s) * RSortedInj%Val(s)
    enddo
    write(*,*)'before update IncPressionWellInj(nWell)',IncPressionWellInj(nWell)
    IncPressionWellInj(nWell) = (-DataWellInjLocal(nWell)%FlowrateImposed + SumMobR)/SumMob
    write(*,*)'after update IncPressionWellInj(nWell)',IncPressionWellInj(nWell)


  end subroutine DefFlash_NonLinPressureUpdateWellInj

  !> \brief Non linear update of Pw and
  !! determine the mode of the projection well
  !! (flowrate or pressure).
  !!
  !! As long as the pressure is less or egal to 
  !! the pressure max, the flowrate of the well 
  !! is imposed. If the pressure is too high, 
  !! the flowrate is no more fixed and the pressure 
  !! is set as Pressure max.
  subroutine DefFlash_NewtonFlashNonLinWellProd

    integer :: nWell, s, j, nums
    double precision :: Flowrate_head

    do nWell=1, NodebyWellProdLocal%Nb

       if(DataWellProdLocal(nWell)%IndWell == 'f') then ! flowrate mode

          ! non linear update of the unknown pressure in well
          call DefFlash_NonLinPressureUpdateWellProd(nWell)
          if(IncPressionWellProd(nWell) < DataWellProdLocal(nWell)%PressionMin) then
             DataWellProdLocal(nWell)%IndWell = 'p'  ! change to pressure mode
             IncPressionWellProd(nWell) = DataWellProdLocal(nWell)%PressionMin  ! Pw = PwMin
          endif

       else if(DataWellProdLocal(nWell)%IndWell == 'p') then ! pressure mode

          if(IncPressionWellProd(nWell) < DataWellProdLocal(nWell)%PressionMin) then
             IncPressionWellProd(nWell) = DataWellProdLocal(nWell)%PressionMin  ! With Newton inc, Pw<Pmin, then change it
          endif
          ! compute the new flowrate at the head node of well nWell
          call DefFlash_PressureToFlowrateWellProd(nWell, Flowrate_head)
          if (Flowrate_head > DataWellProdLocal(nWell)%FlowrateImposed) then  ! Prod well then DataWellProdLocal(nWell)%FlowrateImposed > 0
             DataWellProdLocal(nWell)%IndWell = 'f'  ! change to flowrate mode
             ! non linear update of the unknown pressure in well
             call DefFlash_NonLinPressureUpdateWellProd(nWell)
          endif

       else
          print*, "Error in Flash Well: no such context"
       end if

    enddo ! well

  end subroutine DefFlash_NewtonFlashNonLinWellProd

  !> \brief Nonlinear update the unknown of production well (Pressure)
  !! by computing Pw such as the well is at equilibrium wit the others
  !! unknows (Qmol_w, matrix,...)
  subroutine DefFlash_NonLinPressureUpdateWellProd(nWell)
    integer, intent(in) :: nWell ! numero of the production well

    double precision :: Ps(NbPhase), Ts, Sat(NbPhase), C(NbComp,NbPhase)
    double precision :: Viscosity, DensiteMolaire, PermRel
    double precision :: dPf, dTf, dCf(NbComp), dSf(NbPhase) ! not used for now, empty passed to f_DensiteMolaire
    double precision :: WIDws, SumMob, SumMobR
    integer :: s, j, nums, m, mph, comptn

    Mob(:,:) = 0.d0
    R(:,:) = 0.d0
    Flow(:) = 0.d0

    comptn = 0

    do s=NodebyWellProdLocal%Pt(nWell)+1, NodebyWellProdLocal%Pt(nWell+1)
       WIDws = NodeDatabyWellProdLocal%Val(s)%WID
       nums = NodebyWellProdLocal%Num(s)
       Ts = IncNode(nums)%Temperature      ! Ts: Temperature in matrix
       Sat(:) = IncNode(nums)%Saturation(:)  ! Sat in matrix

       ! loop over alpha in Q_s
       do m=1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
          comptn = comptn + 1    ! comptn = sum(s) sum(alpha in Q_s) 1

          mph = NumPhasePresente_ctx(m,IncNode(nums)%ic)
          Ps(mph) = IncNode(nums)%Pression      ! Ps: Pressure in matrix (no capillary pressure)    ???
          R(mph,s) = Ps(mph) - PerfoWellProd(s)%PressureDrop     ! R_s,alpha = P_s,alpha^{n} - PressureDrop_{w,s}^{n-1}
          C(:,mph) = IncNode(nums)%Comp(:,mph)    ! Comp in matrix

          ! Molar density
          call f_DensiteMolaire(mph, Ps(mph), Ts, C(:,mph), Sat, &
               DensiteMolaire, dPf, dTf, dCf, dSf)
          ! viscosity
          call f_Viscosite(mph, Ps(mph), Ts, C(:,mph), Sat, &
               Viscosity, dPf, dTf, dCf, dSf)
          ! Permrel
          call f_PermRel(mph, Sat, PermRel, dSf)

          Mob(mph,s) = PermRel * DensiteMolaire / Viscosity * WIDws

          ! initialization of RSorted and MobSortedTempProd before calling QuickSortCSR
          RSortedProd%Num(RSortedProd%Pt(nWell)+comptn) = NodebyWellProdLocal%Pt(nWell)+comptn
          RSortedProd%Val(RSortedProd%Pt(nWell)+comptn) = R(mph,s)
          MobSortedTempProd(RSortedProd%Pt(nWell)+comptn) = Mob(mph,s)
       enddo

    enddo ! node s

    write(*,*)'before sort RSortedProd good size',RSortedProd%Val(RSortedProd%Pt(nWell)+1: RSortedProd%Pt(nWell)+comptn)
    write(*,*)'before sort RSortedProd longer size',RSortedProd%Val(RSortedProd%Pt(nWell)+1: RSortedProd%Pt(nWell+1))


    ! sort CSR R in increasing order of R_s^alpha (carreful, Ri can be egual to Ri+1)
    call QuickSortCSR(RSortedProd, RSortedProd%Pt(nWell)+1, RSortedProd%Pt(nWell)+comptn, 'i')
    ! sort CSR Mob in same order than RSortedProd
    do s=RSortedProd%Pt(nWell)+1, RSortedProd%Pt(nWell)+comptn
       MobSortedProd(s) = MobSortedTempProd(RSortedProd%Num(s))
    enddo

    ! compute Flow(i) = sum_{j=i+1}^{n} MobSortedProd(j) * (RSortedProd%Val(j) - RSortedProd%Val(i))
    ! Flow(i+1)<=Flow(i) due to the order of MobSortedProd and Rsorted
    do s=RSortedProd%Pt(nWell)+1, RSortedProd%Pt(nWell)+comptn
       do j=s+1,RSortedProd%Pt(nWell)+comptn
          Flow(s) = Flow(s) + MobSortedProd(j) * (RSortedProd%Val(j)-RSortedProd%Val(s))
       enddo
    enddo


    write(*,*)'RSortedProd',RSortedProd%Val(RSortedProd%Pt(nWell)+1: RSortedProd%Pt(nWell)+comptn)
    write(*,*)'Flow',Flow(RSortedProd%Pt(nWell)+1: RSortedProd%Pt(nWell)+comptn)


    ! if Qmol_w > Flow(1)   then    Pw < RSortedProd(1)
    if(Flow(RSortedProd%Pt(nWell)+1) <= DataWellProdLocal(nWell)%FlowrateImposed) then
       j=RSortedProd%Pt(nWell)+1

       write(*,*)'Flow(1) ,qmol',Flow(j), DataWellProdLocal(nWell)%FlowrateImposed

    else
       ! search for the index j such that Flow(j+1)<= Qmol_w <= Flow(j)
       ! then r(i) <= Pw <= r(i+1)
       j=RSortedProd%Pt(nWell)+comptn
       do while(Flow(j) <= DataWellProdLocal(nWell)%FlowrateImposed)
          j=j-1
       enddo
       j=j+1
       write(*,*)'Flow(j+1) ,qmol,Flow(j)',Flow(j+1), DataWellProdLocal(nWell)%FlowrateImposed, Flow(j)
    endif

    SumMob = 0.d0
    SumMobR = 0.d0
    do s=j,RSortedProd%Pt(nWell)+comptn
       SumMob = SumMob + MobSortedProd(s)
       SumMobR = SumMobR + MobSortedProd(s) * RSortedProd%Val(s)
    enddo
    write(*,*)'before update IncPressionWellProd(nWell)',IncPressionWellProd(nWell)
    IncPressionWellProd(nWell) = (-DataWellProdLocal(nWell)%FlowrateImposed + SumMobR)/SumMob
    write(*,*)'after update IncPressionWellProd(nWell)',IncPressionWellProd(nWell)


  end subroutine DefFlash_NonLinPressureUpdateWellProd

  !> \brief Determine the mode of the injection well
  !! (flowrate or pressure).
  !!
  !! As long as the pressure is less or egal to 
  !! the pressure max, the flowrate of the well 
  !! is imposed. If the pressure is too high, 
  !! the flowrate is no more fixed and the pressure 
  !! is set as Pressure max.
  subroutine DefFlash_NewtonFlashLinWellInj

    double precision :: Flowrate_head
    integer :: num_Well

    ! print*, NbWellInjLocal_Ncpus

    do num_Well=1, NbWellInjLocal_Ncpus(commRank+1)

       if(DataWellInjLocal(num_Well)%IndWell == 'f') then ! flowrate mode

          if(IncPressionWellInj(num_Well) > DataWellInjLocal(num_Well)%PressionMax) then
             
             DataWellInjLocal(num_Well)%IndWell = 'p'  ! change to pressure mode
             IncPressionWellInj(num_Well) = DataWellInjLocal(num_Well)%PressionMax  ! Pw = PwMax
          endif

       else if(DataWellInjLocal(num_Well)%IndWell == 'p') then ! pressure mode

          if(IncPressionWellInj(num_Well) > DataWellInjLocal(num_Well)%PressionMax) then
             IncPressionWellInj(num_Well) = DataWellInjLocal(num_Well)%PressionMax  ! With Newton inc, Pw>Pmax, then change it
          endif

          ! compute the new flowrate at the head node of well num_Well
          call LoisThermoHydro_divP_wellinj(num_Well) ! update thermo Laws of nodes in well num_Well
          call DefFlash_PressureToFlowrateWellInj(num_Well, Flowrate_head)

          if (abs(Flowrate_head) > abs(DataWellInjLocal(num_Well)%FlowrateImposed)) then  ! inj well then DataWellInjLocal(num_Well)%flowrate < 0
             DataWellInjLocal(num_Well)%IndWell = 'f'  ! change to flowrate mode
          endif
       else

          print*, "Error in Newton Flash Injection Well", num_Well, "  no such index well: ", DataWellInjLocal(num_Well)%IndWell
       end if

    end do ! well

  end subroutine DefFlash_NewtonFlashLinWellInj

  !> \brief Determine the mode of the production well
  !! (flowrate or pressure).
  !!
  !! As long as the pressure is greater or egal to
  !! the pressure min, the flowrate of the well
  !! is imposed. If the pressure is too low,
  !! the flowrate is no more fixed and the pressure 
  !! is set as Pressure min.
  subroutine DefFlash_NewtonFlashLinWellProd

    double precision :: Flowrate_head
    integer :: num_Well

    do num_Well=1, NbWellProdLocal_Ncpus(commRank+1)

       if (DataWellProdLocal(num_Well)%IndWell == 'f') then ! flowrate mode
          
          if (IncPressionWellProd(num_Well) < DataWellProdLocal(num_Well)%PressionMin) then

             DataWellProdLocal(num_Well)%IndWell = 'p'  ! change to pressure mode
             IncPressionWellProd(num_Well) = DataWellProdLocal(num_Well)%PressionMin  ! Pw = PwMin
          endif

       else if(DataWellProdLocal(num_Well)%IndWell == 'p') then ! pressure mode
          
          if(IncPressionWellProd(num_Well) < DataWellProdLocal(num_Well)%PressionMin) then
             IncPressionWellProd(num_Well) = DataWellProdLocal(num_Well)%PressionMin  ! With Newton inc, Pw<Pmin, then change it
          endif

          ! compute the new flowrate at the head node of the well num_Well
          call DefFlash_PressureToFlowrateWellProd(num_Well, Flowrate_head)

          if (abs(Flowrate_head) > abs(DataWellProdLocal(num_Well)%FlowrateImposed)) then  ! Prod well then DataWellProdLocal(num_Well)%flowrate > 0
             DataWellProdLocal(num_Well)%IndWell = 'f'  ! change to flowrate mode
          endif
       else

          print*, "Error in Newton Flash Production Well", num_Well, "  no such index well: ", DataWellProdLocal(num_Well)%IndWell
       end if

    enddo ! well

  end subroutine DefFlash_NewtonFlashLinWellProd

  !> \brief Compute the flowrate at the each node of the production well num_Well:
  !! Fill the CSR vectors summolarFluxProd and sumnrjFluxProd.
  !! Use PerfoWellProd\%Pression and IncNode\%Pression.
  !!
  !! \param[in]      num_Well             Numero of the production well
  subroutine DefFlash_FlowrateWellProd(num_Well)

    integer, intent(in) :: num_Well

    integer :: s, nums, icp, m, mph, sparent
    double precision :: Pws, Ps, WIDws
    double precision:: Flux_ks(NbComp), FluxT_ks

    ! nodes of well k
    do s=NodebyWellProdLocal%Pt(num_Well)+1, NodebyWellProdLocal%Pt(num_Well+1)
       nums = NodebyWellProdLocal%Num(s)

       sparent = NodeDatabyWellProdLocal%Val(s)%PtParent

       Pws = PerfoWellProd(s)%Pression ! P_{w,s}
       Ps = IncNode(nums)%Pression     ! P_s
       WIDws = NodeDatabyWellProdLocal%Val(s)%WID

       Flux_ks(:) = 0.d0
       FluxT_ks = 0.d0

       do m=1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
          mph = NumPhasePresente_ctx(m,IncNode(nums)%ic)

          !
          if ((Ps - Pws) > 0.d0) then
             do icp=1, NbComp
                if (MCP(icp,mph) == 1) then ! \cap P_i
                   Flux_ks(icp) = Flux_ks(icp) + DensitemolaireKrViscoCompNode(icp,m,nums) * WIDws * (Ps - Pws)
                end if
             end do

#ifdef _THERMIQUE_
             FluxT_ks = FluxT_ks + DensitemolaireKrViscoEnthalpieNode(m,nums) * WIDws * (Ps - Pws)
#endif
          end if
       end do

       summolarFluxProd(:,s) = summolarFluxProd(:,s) + Flux_ks(:)

#ifdef _THERMIQUE_
       sumnrjFluxProd(s) = sumnrjFluxProd(s) + FluxT_ks
#endif
       
       if(sparent /= -1) then ! head node if sparent = -1
          summolarFluxProd(:,sparent) = summolarFluxProd(:,sparent) + summolarFluxProd(:,s)
          sumnrjFluxProd(sparent) = sumnrjFluxProd(sparent) + sumnrjFluxProd(s)
       end if
    end do

    ! compute qw, i.e. sum_{i} q_{w,s,i} where s is head
    ! Flowrate_head = 0.d0
    ! do icp=1, NbComp
    !    Flowrate_head = Flowrate_head + summolarFluxProd(icp,NodebyWellProdLocal%Pt(num_Well+1))
    ! end do

  end subroutine DefFlash_FlowrateWellProd

  !> \brief Compute the flowrate at the head node of the production well num_Well
  !! with the new value of the unknows
  !!
  !! Use IncPressionWellProd, PerfoWellProd\%PressureDrop and IncNode\%Pression.
  !!
  !! \param[in]      num_Well             Numero of the production well
  !! \param[out]     Qw                   Flowrate at the head node
  subroutine DefFlash_PressureToFlowrateWellProd(num_Well, Qw)

    integer, intent(in) :: num_Well
    double precision, intent(out) :: Qw

    integer :: s, k, nums, icp, m, mph, &
         rocktypeinc
    double precision :: Pws, Ps, WIDws
    double precision:: Flux_ks(NbComp)

    Qw = 0.d0

    ! update thermo Laws of all nodes
    call LoisThermoHydro_divPrim_nodes

    ! nodes of well num_Well
    do s=NodebyWellProdLocal%Pt(num_Well)+1, NodebyWellProdLocal%Pt(num_Well+1)
       nums = NodebyWellProdLocal%Num(s)

       Pws = IncPressionWellProd(num_Well) + PerfoWellProd(s)%PressureDrop ! P_{w,s}       
       Ps = IncNode(nums)%Pression     ! P_s
       WIDws = NodeDatabyWellProdLocal%Val(s)%WID

       Flux_ks(:) = 0.d0

       do m=1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
          mph = NumPhasePresente_ctx(m,IncNode(nums)%ic)

          !
          if ((Ps - Pws) > 0.d0) then
             do icp=1, NbComp
                if (MCP(icp,mph) == 1) then ! \cap P_i
                   Flux_ks(icp) = Flux_ks(icp) + DensitemolaireKrViscoCompNode(icp,m,nums) * WIDws * (Ps - Pws)
                end if
             end do
          end if
       end do

       do icp=1, NbComp
          Qw = Qw + Flux_ks(icp)
       end do

    end do

  end subroutine DefFlash_PressureToFlowrateWellProd

  !> \brief Compute the flowrate at the head node of the injection well num_Well.
  !! Use PerfoWellInj, IncPressionWellInj, WID and IncNode\%Pression.
  !!
  !! \param[in]      num_Well             Numero of the injection well
  !! \param[out]     Qw    Flowrate at the head node of the injection well
  subroutine DefFlash_PressureToFlowrateWellInj(num_Well, Qw)

    integer, intent(in) :: num_Well
    double precision, intent(inout) :: Qw

    integer :: s, nums, icp, m, mph
    double precision :: Pws, Ps, WIDws
    double precision:: Flux_ks(NbComp)

    Qw = 0.d0

    ! nodes of well num_Well
    do s=NodebyWellInjLocal%Pt(num_Well)+1, NodebyWellInjLocal%Pt(num_Well+1)
       nums = NodebyWellInjLocal%Num(s)

       Pws = IncPressionWellInj(num_Well) + PerfoWellInj(s)%PressureDrop ! P_{w,s}       
       Ps = IncNode(nums)%Pression     ! P_s

       WIDws = NodeDatabyWellInjLocal%Val(s)%WID

       Flux_ks(:) = 0.d0

       do m=1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
          mph = NumPhasePresente_ctx(m,IncNode(nums)%ic)

          if((Ps - Pws) < 0.d0) then ! < 0
             do icp=1, NbComp
                Flux_ks(icp) = DensitemolaireKrViscoCompWellInj(icp,s) * WIDws * (Ps - Pws)
             end do
          end if
       end do

       do icp=1, NbComp
          Qw = Qw + Flux_ks(icp)
       end do
    end do

  end subroutine DefFlash_PressureToFlowrateWellInj

  !! \brief Execute the flash to determine T and the molar fractions
  !! to update PerfoWellProd(s)%Temperature and PerfoWellProd(s)%Density.
  !! This Flash is performed for a monophasique multicomponent fluid.
  !!
  !! Loop over the nodes s from head to tail to
  !! to determine wich phases are present, the temperature and the mean density.
  !! The pressure of the following node depends on the mean density, this is why 
  !! the loop is done from head to tail (mean density is updated before being used).  
  subroutine DefFlash_TimeFlashSinglePhaseWellProd

    double precision :: T, RT, Pws, Ci(NbComp), sumci, E, zp, zs
    double precision :: H, rhoMean, dP_Tsat, Pdrop
    double precision :: dPf, dTf, Sat(NbPhase), dCf(NbComp), dSf(NbPhase) ! not used for now, empty passed to f_Enthalpie
    integer :: k, s, icp, sparent, Nnz, i
    logical :: converged

!    Sat(PHASE_GAS) = 0.d0
    Sat(PHASE_WATER) = 1.d0

    summolarFluxProd(:,:) = 0.d0
    sumnrjFluxProd(:) = 0.d0

    do k=1, NodebyWellProdLocal%Nb

       ! compute flowrate of well k (fill summolarFluxProd and sumnrjFluxProd)
       call DefFlash_FlowrateWellProd(k)

       ! looping from head to queue
       do s=NodebyWellProdLocal%Pt(k+1), NodebyWellProdLocal%Pt(k)+1, -1

          ! compute P_{w,s}
          if(s==NodebyWellProdLocal%Pt(k+1)) then ! head node, P = Pw

             Pws = IncPressionWellProd(k)
             PerfoWellProd(s)%Pression = Pws
             PerfoWellProd(s)%PressureDrop = 0.d0

          else ! Pws = P_{w,parent} + \Delta P_{w,parent}

             zs = XNodeLocal(3,NodebyWellProdLocal%Num(s)) ! z-cordinate of node s          
             zp = XNodeLocal(3,NodeDatabyWellProdLocal%Val(s)%Parent) ! z-cordinate of parent of s

             sparent = NodeDatabyWellProdLocal%Val(s)%PtParent ! parent pointer

             ! as the loop is done from head to queue, %Density is updated before being used
             Pdrop = PerfoWellProd(sparent)%Density * Gravite * (zp - zs)             
             Pws = PerfoWellProd(sparent)%Pression + Pdrop ! Pws

             PerfoWellProd(s)%Pression = Pws
             PerfoWellProd(s)%PressureDrop = PerfoWellProd(sparent)%PressureDrop + Pdrop
          end if

          ! Newton method to compute T: R = E-Enthalpie*n = 0
          E = sumnrjFluxProd(s) ! energy

          do icp=1, NbComp
             sumci = summolarFluxProd(icp,s) ! sum_i {n_i}
          end do
          
          ! initialize newton with Tsat
          call DefModel_Tsat(PerfoWellProd(s)%Pression, T, dP_Tsat)

          converged = .false.
          
          do i=1, Maxiter

             call f_Enthalpie(PHASE_WATER, Pws, T, Ci(:), Sat(:), &
                  H, dPf, dTf, dCf, dSf)

             RT = E - H * sumci ! residu

             if (abs(RT) < Tol) then
                converged = .true.                
                exit
             else
                T = T + RT / (dTf * sumci)
             end if
          end do

          ! check if the Newton has converged
          if(converged .eqv. .false.) then
             print*, "Warning: Newton in DefFlash_TimeFlashSinglePhaseProd has not converged"
             print*, "Residu is", abs(RT), "k= ", k, "s= ", s
          end if

          ! update PhysPerfoWell
          PerfoWellProd(s)%Temperature = T
          call f_DensiteMassique(PHASE_WATER, Pws, T, Ci, Sat, rhoMean, &
               dPf, dTf, dCf, dSf)
          PerfoWellProd(s)%Density = rhoMean

       end do ! end of well k
    end do

  end subroutine DefFlash_TimeFlashSinglePhaseWellProd

  !! \brief Execute the flash to determine T and the molar fractions
  !! to update PerfoWellProd(s)%Temperature and PerfoWellProd(s)%Density.
  !! This Flash is performed for a diphasique monocomponent fluid.
  !!
  !! Loop over the nodes s from head to tail to compute the thermodynamical flash 
  !! to determine wich phases are present, the temperature and the mean density.
  !! The pressure of the following node depends on the mean density, this is why 
  !! the loop is done from head to tail (mean density is updated before being used).
  subroutine DefFlash_TimeFlashTwoPhasesProd

    double precision :: Temp, Pws, Tsat, dP_Tsat, liq_molarfrac, zs, zp
    double precision :: totalMolarFlux, Hgas, Hliq, Hphase, rhogas, rholiq
    double precision :: sumci, E, Res, Pdrop
    ! not used, empty passed to f_Enthalpie
    double precision :: dPf, dTf, sat(NbPhase), molarFrac(NbComp), dCf(NbComp), dSf(NbPhase)
    integer :: Nnz, nWell, s, sparent, icp, i, ID_PHASE ! ID_PHASE=(-1 if diphasique, PHASE_GAS if gas, PHASE_WATER if liq)

    summolarFluxProd(:,:) = 0.d0
    sumnrjFluxProd(:) = 0.d0

    ! loop over production well
    do nWell=1, NodebyWellProdLocal%Nb

       ! compute flowrate of well nWell (fill summolarFluxProd and sumnrjFluxProd)
       call DefFlash_FlowrateWellProd(nWell)

       ! looping from head to queue (to update Pws with the value of rho computed in current flash)
       do s=NodebyWellProdLocal%Pt(nWell+1), NodebyWellProdLocal%Pt(nWell)+1, -1

          ! compute P_{w,s}
          if(s==NodebyWellProdLocal%Pt(nWell+1)) then ! head node, P = Pw

             Pws = IncPressionWellProd(NWell)
             PerfoWellProd(s)%Pression = Pws
             PerfoWellProd(s)%PressureDrop = 0.d0

          else ! Pws = P_{w,parent} + \Delta P_{w,parent}

             zs = XNodeLocal(3,NodebyWellProdLocal%Num(s)) ! z-cordinate of node s          
             zp = XNodeLocal(3,NodeDatabyWellProdLocal%Val(s)%Parent) ! z-cordinate of parent of s

             sparent = NodeDatabyWellProdLocal%Val(s)%PtParent ! parent pointer

             ! as the loop is done from head to queue, %Density is updated before being used
             Pdrop = PerfoWellProd(sparent)%Density * Gravite * (zp - zs)             
             Pws = PerfoWellProd(sparent)%Pression + Pdrop ! Pws

             PerfoWellProd(s)%Pression = Pws
             PerfoWellProd(s)%PressureDrop = PerfoWellProd(sparent)%PressureDrop + Pdrop
          end if

          ! Flash
          ! we sum the molar flux at each node
          totalMolarFlux = 0.d0
          do icp=1, NbComp
             totalMolarFlux = totalMolarFlux + summolarFluxProd(icp, s) ! total number of moles at node s
          end do

          ! suppose that the two phases are present at the node
          ID_PHASE = -1
          ! then T=Tsat(P)
#ifdef _THERMIQUE_
          call DefModel_Tsat(Pws, Tsat, dP_Tsat)
#endif
          Temp = Tsat

          ! thus compute liq_molarfrac thanks to the energy, and the enthalpies
          ! molarFrac is not used in the computation of the enthalpies
#ifdef _THERMIQUE_
!          call f_Enthalpie(PHASE_GAS, Pws, Temp, molarFrac, sat, Hgas, dPf, dTf, dCf, dSf)
          call f_Enthalpie(PHASE_WATER, Pws, Temp, molarFrac, sat, Hliq, dPf, dTf, dCf, dSf)
#endif
          ! and compute liq_molarfrac
          liq_molarfrac = (Hgas - sumnrjFluxProd(s)/totalMolarFlux) / (Hgas-Hliq)

          if(liq_molarfrac<0.d0) then  ! the hypothesis that the two phases are present is wrong: only gas
             liq_molarfrac = 0.d0
             ID_PHASE = PHASE_GAS
          else if(liq_molarfrac>1.d0) then  ! the hypothesis that the two phases are present is wrong: only liquid
             liq_molarfrac = 1.d0
             ID_PHASE = PHASE_WATER
          endif

          if(ID_PHASE > 0) then

             ! Newton algo on Energy=n*Enthalpie(P,T) to determine T.
             ! Is initialized with Temp = Tsat
             E = sumnrjFluxProd(s) ! energy

             do icp=1, NbComp
                sumci = summolarFluxProd(icp,s) ! sum_i {n_i}
             end do

             do i=1, Maxiter
#ifdef _THERMIQUE_
                call f_Enthalpie(ID_PHASE, Pws, Temp, molarFrac, sat, &    ! molarFrac is not used in Enthalpie
                     Hphase, dPf, dTf, dCf, dSf)
#endif

                Res = E - Hphase * sumci ! residu
                if (abs(Res) < Tol) then
                   exit
                else
                   Temp = Temp - Res / (dTf*sumci)
                end if
             end do

             ! check if the Newton converged
             if(i==Maxiter) then
                print*, "Warning: Newton in DefFlash_TimeFlashTwoPhasesProd has not converged"
             end if

          endif

          ! we deduce the mean density
          ! molarFrac is not used in the computation of the massique densities
          call f_DensiteMassique(PHASE_GAS,Pws,Temp,molarFrac,sat,rhogas,dPf,dTf,dCf,dSf)
          call f_DensiteMassique(PHASE_WATER,Pws,Temp,molarFrac,sat,rholiq,dPf,dTf,dCf,dSf)
          PerfoWellProd(s)%Density = liq_molarfrac * rholiq + (1.d0-liq_molarfrac) * rhogas

          ! fill PhysPerfoWell%T
          PerfoWellProd(s)%Temperature = Temp
       enddo ! node s

    enddo ! nWell

  end subroutine DefFlash_TimeFlashTwoPhasesProd

  !> \brief Sorting the heights contained in mycsr%Value, and update 
  !! the corresponding node values (indexes) stored in mycsr%Num
  !! mycsr%Pt is constructed on NodebyWellInjLocal
  recursive subroutine QuickSortCSR(myCSR, left, right, mode)
    use commontype
    implicit none
    type(CSRdble), intent(inout) :: myCSR 
    integer, intent(in) :: left, right
    character(len=1), intent(in) :: mode
    double precision :: x, tmp_val
    integer :: i, j, tmp_num

    ! write(fdFl, *) '--> quicksortCSR called with left: ', left, ' right: ', right
    ! we are sorting CSR%Num(f), f=CSR%Pt(k)+1,CSR%Pt(k+1), for well k

    x = myCSR%Val( (left + right) / 2 ) ! pivot point
    ! write(fdFl, *) 'pivot is ', x
    i = left; j = right ! left and right boundaries
    do
       if (mode == 'd') then
          do while (x < myCSR%Val(i))
             i = i + 1
          end do
          do while (x > myCSR%Val(j))
             j = j - 1
          end do
       elseif (mode == 'i') then 
          do while (x > myCSR%Val(i))
             i = i + 1
          end do
          do while (x < myCSR%Val(j))
             j = j - 1
          end do
       else
          write(*,*) 'WARNING: sorting mode unknown !'
          return
       end if
       if (i >= j) then 
          exit
       end if

       ! swapping the %Val field and %Num fields
       tmp_val = myCSR%Val(i); tmp_num = myCSR%Num(i)
       myCSR%Val(i) = myCSR%Val(j); myCSR%Num(i) = myCSR%Num(j)
       myCSR%Val(j) = tmp_val; myCSR%Num(j) = tmp_num
       i = i + 1; j = j - 1
    end do
    if (left < i-1) then 
       call QuickSortCSR(myCSR, left, i-1, mode)
    end if
    if (j+1 < right) then
       call QuickSortCSR(myCSR, j+1, right, mode)
    end if
  end subroutine QuickSortCSR

  !> \brief Update PressureDrop and Pression of all the nodes of a given well
  ! we start by the head node and update the new value using rho*gravity*delta_z increments
  ! LIMITATIONS: the head node of a given well must be at the top height
  ! pressure update loop:  p^{n+1} = p^{n} + rho(p^{n}) * g * (z^{n+1} - z^{n})
  ! (z axis)
  ! ^
  ! | ----- (PerfoWellInj(pts1)%Pression)
  ! |    ^
  ! |    .
  ! |    .
  ! |    delta(P)
  ! |    .
  ! |    .
  ! |    .
  ! |    v
  ! | ----- (PerfoWellInj(pts2)%Pression)

  ! subroutine DefFlash_PressureDropInj

  !   double precision :: Sat(NbPhase), dCf(NbComp), dSf(NbPhase)
  !   double precision :: Rhotmp, dz, Ptmp, ztmp, dPf, dTf
  !   integer :: i, s, k, nums1, nums2, pts1, pts2

  !   do k=1, NbWellInjLocal_Ncpus(commRank+1)

  !      ! initialize the head pressure drop to zero, and the head pressure 
  !     ! with the Unknown value for the Well
  !     PerfoWellInj(ZSortedInj%Num(ZSortedInj%Pt(k)+1))%PressureDrop = 0.d0
  !     PerfoWellInj(ZSortedInj%Num(ZSortedInj%Pt(k)+1))%Pression = IncPressionWellInj(k)

  !     do s=ZSortedInj%Pt(k)+1, ZSortedInj%Pt(k+1) - 1
  !       write(fdFl, *) '.. .. .. .. .. .. .. .. .. .. .. .. .. .. ..'
  !       ! node pointer
  !       pts1 = ZsortedInj%Num(s)
  !       pts2 = ZsortedInj%Num(s+1)
  !       ! local node number
  !       nums1 = NodebyWellInjLocal%Num(pts1)
  !       nums2 = NodebyWellInjLocal%Num(pts2)

  !       Sat(PHASE_GAS) = 0.d0
  !       Sat(PHASE_WATER) = 1.d0

  !       Ptmp = PerfoWellInj(pts1)%Pression
  !       ztmp = XNodeLocal(3, nums1) ! just to check we are OK

  !       dz = (XNodeLocal(3, nums2) - XNodeLocal(3, nums1)) / dble(Nslice)
  !       ! pressure update loop:  p^{n+1} = p^{n} + rho(p^{n}) * g * (z^{n+1} - z^{n})
  !       do i=1, Nslice
  !         ztmp = ztmp + dz
  !         call f_DensiteMolaire(PHASE_WATER, Ptmp, PerfoWellInj(pts1)%Temperature, &
  !             DataWellInjLocal(k)%CompTotal, Sat, Rhotmp, dPf, dTf, dCf, dSf)
  !         ! gravity points downwards, heights points upwards, hence the negative sign
  !         Ptmp = Ptmp - Gravite * Rhotmp * dz
  !       end do

  !       PerfoWellInj(pts2)%Pression = Ptmp
  !       PerfoWellInj(pts2)%PressureDrop = PerfoWellInj(pts2)%Pression - IncPressionWellInj(k)

  !       ! debugging
  !       write(fdFl, '(4(a,i2),a)') 'nodes (', pts1, 'ptL,', nums1, &
  !           'numL) -> (', pts2, 'ptL,', nums2, 'numL)'

  !       write(fdFl, '(a,f6.2,a,f6.2,f6.2,a,e14.6,a,e14.6)'), &
  !           'integrate z1: ', XNodeLocal(3, nums1), ' -> z2: ', XNodeLocal(3, nums2), &
  !           ztmp, ' P1: ', PerfoWellInj(pts1)%Pression, ' P2: ', Ptmp
  !       write(fdFl, '(a,e14.6)') 'pressure drop', PerfoWellInj(pts2)%PressureDrop 
  !     end do
  !   end do
  ! end subroutine DefFlash_PressureDropInj

end module DefFlash
