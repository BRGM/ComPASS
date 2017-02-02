module Residu

  use DefModel
  use IncCV
  use LoisThermoHydro
  use Flux

  implicit none

  ! Residu cell, fracture face, node
  double precision, allocatable, dimension(:,:), protected :: &
       ResiduCell,    &
       ResiduFrac,    &
       ResiduNode

  ! Residu injection well and production well
  double precision, allocatable, dimension(:), protected :: &
       ResiduWellInj, &
       ResiduWellProd

  ! AccVol of time step n-1
  double precision, allocatable, dimension(:,:), private :: &
       AccVolCell_1, &
       AccVolFrac_1, &
       AccVolNode_1

  public :: &
       Residu_allocate,    &
       Residu_free,        &
       Residu_compute,     &
       Residu_RelativeNorm

  private :: &
       Residu_AccVol, &
       Residu_Conservations

contains

  ! compute Residu
  subroutine Residu_compute(Delta_t, IterNewton)

    double precision, intent(in) :: Delta_t
    integer, intent(in) :: IterNewton

    integer :: i, k

    ! init Residu as zero
    ResiduCell(:,:) = 0.d0
    ResiduFrac(:,:) = 0.d0
    ResiduNode(:,:) = 0.d0

    ResiduWellInj(:) = 0.d0
    ResiduWellProd(:) = 0.d0

    ! term AccVol
    call Residu_AccVol

    ! Newton is initialized with the unknown at time step n-1
    if(IterNewton==1) then

       do k=1, NbCellLocal_Ncpus(commRank+1) ! cell
          do i=1, NbCompThermique
             AccVolCell_1(i,k) = IncCell(k)%AccVol(i)
          end do
       end do

       do k=1, NbFracLocal_Ncpus(commRank+1) ! frac
          do i=1, NbCompThermique
             AccVolFrac_1(i,k) = IncFrac(k)%AccVol(i)
          end do
       end do

       do k=1, NbNodeLocal_Ncpus(commRank+1) ! node
          do i=1, NbCompThermique
             AccVolNode_1(i,k) = IncNode(k)%AccVol(i)
          end do
       end do
    end if

    ! Residu from terms of accumulation
    do k=1, NbCellLocal_Ncpus(commRank+1) ! cell
       do i=1, NbCompThermique
          ResiduCell(i,k) = &
               ( IncCell(k)%AccVol(i)-AccVolCell_1(i,k) ) / Delta_t
       end do
    end do

    do k=1, NbFracLocal_Ncpus(commRank+1) ! Frac
       do i=1, NbCompThermique
          ResiduFrac(i,k) = &
               ( IncFrac(k)%AccVol(i)-AccVolFrac_1(i,k) ) / Delta_t
       end do
    end do

    do k=1, NbNodeLocal_Ncpus(commRank+1) ! node 
       do i=1, NbCompThermique
          ResiduNode(i,k) = &
               ( IncNode(k)%AccVol(i)-AccVolNode_1(i,k) ) / Delta_t
       end do
    end do

    ! Residu for dir node will be reset as 0 in the end of this subroutine
    ! We do not consider much about dir nodes during computing the residu such that
    ! the code is more clear

    ! Residu from conservation composants
    call Residu_Conservations

    ! Residu for dir node should be reset as 0 for computing norm of Residu
#ifdef _THERMIQUE_
    do k=1, NbNodeLocal_Ncpus(commRank+1) ! node

       if(IdNodeLocal(k)%P=="d") then
          ResiduNode(1,k) = 0.d0
          ResiduNode(3:NbCompThermique,k) = 0.d0
       end if

       if(IdNodeLocal(k)%T=="d") then
          ResiduNode(2,k) = 0.d0
       end if
    end do
#else
    do k=1, NbNodeLocal_Ncpus(commRank+1) 
       if(IdNodeLocal(k)%P=="d") then
          ResiduNode(1:NbCompThermique,k) = 0.d0
       end if
    end do
#endif

    ! if(commRank==1) then
    !    open(unit=12,file="log.txt",status="unknown")

    !    do k=1, 3, 2
    !       do i=1, NbCompThermique
    !          write(12,'(ES22.12)') ResiduNode(i,k)
    !          write(*,'(ES22.12)') ResiduNode(i,k)
    !       end do
    !       write(12,*) ""
    !       write(*,*) ""
    !    end do
    !    close(12)
    ! end if

  end subroutine Residu_compute


  ! allocate
  subroutine Residu_allocate

    ! For node/frac/well, both of own and ghost are computed,
    ! however, only the own are correct and used in second member of Jacobian
    allocate( ResiduCell(NbCompThermique, NbCellLocal_Ncpus(commRank+1)))
    allocate( ResiduFrac(NbCompThermique, NbFracLocal_Ncpus(commRank+1)))
    allocate( ResiduNode(NbCompThermique, NbNodeLocal_Ncpus(commRank+1)))

    allocate( ResiduWellInj(NbWellInjLocal_Ncpus(commRank+1)))
    allocate( ResiduWellProd(NbWellProdLocal_Ncpus(commRank+1)))

    allocate( AccVolCell_1(NbCompThermique, NbCellLocal_Ncpus(commRank+1)))
    allocate( AccVolFrac_1(NbCompThermique, NbFracLocal_Ncpus(commRank+1)))
    allocate( AccVolNode_1(NbCompThermique, NbNodeLocal_Ncpus(commRank+1)))

  end subroutine Residu_allocate


  subroutine Residu_AccVol

    integer :: k, m, mph, i, ic, icp
    integer :: tmpAccVol(NbComp)

    ! Cell    
    do k=1, NbCellLocal_Ncpus(commRank+1)

       ic = IncCell(k)%ic

       ! if icp is not in Ctilde, we compute its values,
       !    firstly we need AccVol(:) = 0
       ! But, if icp is in Ctilde, it is an unknown independent
       !     we must keep its values
       do i=1, NbCompCtilde_ctx(ic)
          icp = NumCompCtilde_ctx(i,ic)
          tmpAccVol(i) = IncCell(k)%AccVol(icp)
       end do

       IncCell(k)%AccVol(1:NbCompThermique) = 0.d0

       do i=1, NbCompCtilde_ctx(ic)
          icp = NumCompCtilde_ctx(i,ic)
          IncCell(k)%AccVol(icp) = tmpAccVol(i)
       end do

       do m=1, NbPhasePresente_ctx(ic) ! Q_k
          mph = NumPhasePresente_ctx(m,ic)

          do icp=1, NbComp
             if(MCP(icp,mph)==1) then ! Q_k \cap P_i

                IncCell(k)%AccVol(icp) = IncCell(k)%AccVol(icp) &
                     + PoroVolDarcyCell(k) * DensitemolaireSatCompCell(icp,m,k)
             end if
          end do

#ifdef _THERMIQUE_

          IncCell(k)%AccVol(NbComp+1) = IncCell(k)%AccVol(NbComp+1) &
               + PoroVolFourierCell(k) * DensitemolaireEnergieInterneSatCell(m,k)
#endif          
       end do ! end of phase m

#ifdef _THERMIQUE_

       IncCell(k)%AccVol(NbComp+1) = IncCell(k)%AccVol(NbComp+1) &
            + Poro_1volFourierCell(k) * CpRoche * IncCell(k)%Temperature
#endif

    end do ! end of cell k    

    ! Frac    
    do k=1, NbFracLocal_Ncpus(commRank+1)

       ic = IncFrac(k)%ic

       do i=1, NbCompCtilde_ctx(ic)
          icp = NumCompCtilde_ctx(i,ic)
          tmpAccVol(i) = IncFrac(k)%AccVol(icp)
       end do

       IncFrac(k)%AccVol(1:NbCompThermique) = 0.d0

       do i=1, NbCompCtilde_ctx(ic)
          icp = NumCompCtilde_ctx(i,ic)
          IncFrac(k)%AccVol(icp) = tmpAccVol(i)
       end do

       do m=1, NbPhasePresente_ctx(ic) ! Q_k
          mph = NumPhasePresente_ctx(m,ic)

          do icp=1, NbComp
             if(MCP(icp,mph)==1) then ! Q_k \cap P_i

                IncFrac(k)%AccVol(icp) = IncFrac(k)%AccVol(icp) &
                     + PoroVolDarcyFrac(k) * DensitemolaireSatCompFrac(icp,m,k)
             end if
          end do

#ifdef _THERMIQUE_

          IncFrac(k)%AccVol(NbComp+1) = IncFrac(k)%AccVol(NbComp+1) &
               + PoroVolFourierFrac(k) * DensitemolaireEnergieInterneSatFrac(m,k)
#endif
       end do ! end of phase m

#ifdef _THERMIQUE_

       IncFrac(k)%AccVol(NbComp+1) = IncFrac(k)%AccVol(NbComp+1) &
            + Poro_1volFourierFrac(k) * CpRoche * IncFrac(k)%Temperature
#endif
    end do ! end of frac k


    ! Node    
    do k=1, NbNodeLocal_Ncpus(commRank+1)

       ic = IncNode(k)%ic

       do i=1, NbCompCtilde_ctx(ic)
          icp = NumCompCtilde_ctx(i,ic)
          tmpAccVol(i) = IncNode(k)%AccVol(icp)
       end do

       IncNode(k)%AccVol(1:NbCompThermique) = 0.d0

       do i=1, NbCompCtilde_ctx(ic)
          icp = NumCompCtilde_ctx(i,ic)
          IncNode(k)%AccVol(icp) = tmpAccVol(i)
       end do

       do m=1, NbPhasePresente_ctx(IncNode(k)%ic) ! Q_k
          mph = NumPhasePresente_ctx(m,IncNode(k)%ic)

          do icp=1, NbComp
             if(MCP(icp,mph)==1) then ! Q_k \cap P_i

                IncNode(k)%AccVol(icp) = IncNode(k)%AccVol(icp) &
                     + PoroVolDarcyNode(k) * DensitemolaireSatCompNode(icp,m,k)
             end if
          end do

#ifdef _THERMIQUE_

          IncNode(k)%AccVol(NbComp+1) = IncNode(k)%AccVol(NbComp+1) &
               + PoroVolFourierNode(k) * DensitemolaireEnergieInterneSatNode(m,k)
#endif
       end do ! end of phase m

#ifdef _THERMIQUE_

       IncNode(k)%AccVol(NbComp+1) = IncNode(k)%AccVol(NbComp+1) &
            + Poro_1volFourierNode(k) * CpRoche * IncNode(k)%Temperature
#endif
    end do ! end of node k

  end subroutine Residu_AccVol


  ! compute residu conservation components and energy
  subroutine Residu_Conservations

    integer :: k, s, fs, fk, nums, m, mph, icp
    integer :: NbNodeCell, NbFracCell, NbNodeFrac, nbNodeWell

    double precision :: Flux_ks(NbComp), FluxT_ks
    double precision :: Pws, Tws, Ps, Ts, WIDws, WIFws, qw, Ps_Pws

    do k=1, NbCellLocal_Ncpus(commRank+1)

       ! number of nodes/fracs in cell k
       NbNodeCell = NodebyCellLocal%Pt(k+1) - NodebyCellLocal%Pt(k)
       NbFracCell = FracbyCellLocal%Pt(k+1) - FracbyCellLocal%Pt(k)

       ! s is node in cell k
       do s=1, NbNodeCell
          nums = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+s)

          Flux_ks(:) = 0.d0
          FluxT_ks = 0.d0

          ! alpha in Q_{K_{k,s}^{alpha}} =
          !      ( {alpha | alpha in Q_k} \cap {alpha | K_{k,s}^{alpha}=k})
          ! \cup ( {alpha | alpha in Q_s} \cap {alpha | K_{k,s}^{alpha}=s})
          do m=1, NbPhasePresente_ctx( IncCell(k)%ic) ! Q_k
             mph = NumPhasePresente_ctx(m, IncCell(k)%ic) 

             if( FluxDarcyKI(mph,s,k)>=0.d0) then ! K_{k,s}^{alpha}=k

                do icp=1, NbComp
                   if(MCP(icp,mph)==1) then ! \cap P_i

                      Flux_ks(icp) = Flux_ks(icp) &
                           + DensitemolaireKrViscoCompCell(icp,m,k) * FluxDarcyKI(mph,s,k)

                      ! if(commRank==1 .and. k==1) then
                      !    print*, s, mph, FluxDarcyKI(mph,s,k)*DensitemolaireKrViscoCompCell(icp,m,k)
                      ! end if

                   end if
                end do

#ifdef _THERMIQUE_

                FluxT_ks = FluxT_ks + DensitemolaireKrViscoEnthalpieCell(m,k) &
                     * FluxDarcyKI(mph,s,k)

                ! if(commRank==1 .and. k==1) then
                !    print*, s, mph, FluxDarcyKI(mph,s,k)* DensitemolaireKrViscoEnthalpieCell(m,k)
                ! end if
#endif

             end if
          end do

          do m=1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
             mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)

             if(FluxDarcyKI(mph,s,k)<0.d0) then ! K_{k,s}^{alpha}=s             

                do icp=1, NbComp
                   if(MCP(icp,mph)==1) then ! \cap P_i

                      Flux_ks(icp) = Flux_ks(icp) &
                           + DensitemolaireKrViscoCompNode(icp,m,nums) * FluxDarcyKI(mph,s,k)

                      ! if(commRank==1 .and. k==15813) then
                      !    print*, k,s, nums, Flux_ks
                      !    ! print*, XNodeLocal(:,NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+s))
                      ! end if

                      ! if(commRank==1 .and. k==1) then
                      !    print*, s, mph, FluxDarcyKI(mph,s,k)* DensitemolaireKrViscoCompNode(icp,m,nums)
                      ! end if

                   end if
                end do

#ifdef _THERMIQUE_

                FluxT_ks = FluxT_ks + DensitemolaireKrViscoEnthalpieNode(m,nums) &
                     * FluxDarcyKI(mph,s,k)
#endif

             end if
          end do

          ! if(commRank==1 .and. k==15813) then
          !    print*, k,s, Flux_ks, &
          !         XNodeLocal(:,NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+s))
          ! end if

          ResiduCell(1:NbComp,k) = ResiduCell(1:NbComp,k) + Flux_ks(:)
          ResiduNode(1:NbComp,nums) = ResiduNode(1:NbComp,nums) - Flux_ks(:) 
          
#ifdef _THERMIQUE_

          ResiduCell(NbComp+1,k) = ResiduCell(NbComp+1,k) &
               + FluxT_ks &
               + FluxFourierKI(s,k)

          ResiduNode(NbComp+1,nums) = ResiduNode(NbComp+1,nums) &
               - FluxT_ks &
               - FluxFourierKI(s,k)
#endif

       end do ! end of node s in cell k

       ! s is frac in cell k
       do s=1, NbFracCell

          fs = FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+s) 
          nums = FaceToFracLocal(fs) ! fs is face number, nums is frac number

          Flux_ks(:) = 0.d0
          FluxT_ks = 0.d0

          ! alpha in Q_{K_{k,s}^{alpha}} =
          !      ( {alpha | alpha in Q_k} \cap {alpha | K_{k,s}^{alpha}=k})
          ! \cup ( {alpha | alpha in Q_s} \cap {alpha | K_{k,s}^{alpha}=s})
          do m=1, NbPhasePresente_ctx(IncCell(k)%ic) ! Q_k
             mph = NumPhasePresente_ctx(m, IncCell(k)%ic) 

             if(FluxDarcyKI(mph,s+NbNodeCell,k)>=0.d0) then ! K_{k,s}^{alpha}=k

                do icp=1, NbComp
                   if(MCP(icp,mph)==1) then ! \cap P_i

                      Flux_ks(icp) = Flux_ks(icp) &
                           + DensitemolaireKrViscoCompCell(icp,m,k) * FluxDarcyKI(mph,s+NbNodeCell,k)

                      ! if(commRank==1 .and. k==1) then
                      !    print*, s, mph, DensitemolaireKrViscoCompCell(icp,m,k)*FluxDarcyKI(mph,s+NbNodeCell,k)
                      ! end if

                   end if
                end do

#ifdef _THERMIQUE_

                FluxT_ks = FluxT_ks + DensitemolaireKrViscoEnthalpieCell(m,k) &
                     * FluxDarcyKI(mph,s+NbNodeCell,k)
#endif

             end if
          end do

          do m=1, NbPhasePresente_ctx(IncFrac(nums)%ic) ! Q_s
             mph = NumPhasePresente_ctx(m, IncFrac(nums)%ic)

             if(FluxDarcyKI(mph,s+NbNodeCell,k)<0.d0) then ! K_{k,s}^{alpha}=s             

                do icp=1, NbComp
                   if(MCP(icp,mph)==1) then ! \cap P_i

                      Flux_ks(icp) = Flux_ks(icp) &
                           + DensitemolaireKrViscoCompFrac(icp,m,nums) * FluxDarcyKI(mph,s+NbNodeCell,k)

                      ! if(commRank==1 .and. k==1) then
                      !    print*, s, mph, DensitemolaireKrViscoCompFrac(icp,m,nums)* FluxDarcyKI(mph,s+NbNodeCell,k)
                      ! end if

                   end if
                end do

#ifdef _THERMIQUE_

                FluxT_ks = FluxT_ks + DensitemolaireKrViscoEnthalpieFrac(m,nums) &
                     * FluxDarcyKI(mph,s+NbNodeCell,k)
#endif

             end if
          end do

          ResiduCell(1:NbComp,k) = ResiduCell(1:NbComp,k) + Flux_ks(:)
          ResiduFrac(1:NbComp,nums) = ResiduFrac(1:NbComp,nums) - Flux_ks(:) 

#ifdef _THERMIQUE_

          ResiduCell(NbComp+1,k) = ResiduCell(NbComp+1,k) &
               + FluxT_ks &
               + FluxFourierKI(s+NbNodeCell,k)

          ResiduFrac(NbComp+1,nums) = ResiduFrac(NbComp+1,nums) &
               - FluxT_ks &
               - FluxFourierKI(s+NbNodeCell,k)
#endif

       end do ! end of frac s in cell k

    end do ! end of cell k

    do k=1, NbFracLocal_Ncpus(commRank+1)

       fk = FracToFaceLocal(k) ! fk is face number

       ! number of nodes in a frac
       NbNodeFrac = NodebyFaceLocal%Pt(fk+1) - NodebyFaceLocal%Pt(fk)

       do s=1, NbNodeFrac
          nums = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk)+s)

          Flux_ks(:) = 0.d0
          FluxT_ks = 0.d0

          ! alpha in Q_{K_{k,s}^{alpha}} =
          !      ( {alpha | alpha in Q_k} \cap {alpha | K_{k,s}^{alpha}=k})
          ! \cup ( {alpha | alpha in Q_s} \cap {alpha | K_{k,s}^{alpha}=s})
          do m=1, NbPhasePresente_ctx(IncFrac(k)%ic) ! Q_k
             mph = NumPhasePresente_ctx(m, IncFrac(k)%ic) 

             if(FluxDarcyFI(mph,s,k)>=0.d0) then ! K_{k,s}^{alpha}=k

                do icp=1, NbComp
                   if(MCP(icp,mph)==1) then ! \cap P_i

                      Flux_ks(icp) = Flux_ks(icp) &
                           + DensitemolaireKrViscoCompFrac(icp,m,k) * FluxDarcyFI(mph,s,k)

                      ! if(commRank==1 .and. k==1) then
                      !    print*, k, s, mph, DensitemolaireKrViscoCompFrac(icp,m,k), FluxDarcyFI(mph,s,k)
                      ! end if

                   end if
                end do

#ifdef _THERMIQUE_

                FluxT_ks = FluxT_ks + DensitemolaireKrViscoEnthalpieFrac(m,k) &
                     * FluxDarcyFI(mph,s,k)
#endif                
             end if
          end do

          do m=1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
             mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)

             if(FluxDarcyFI(mph,s,k)<0.d0) then ! K_{k,s}^{alpha}=s             

                do icp=1, NbComp
                   if(MCP(icp,mph)==1) then ! \cap P_i

                      Flux_ks(icp) = Flux_ks(icp) &
                           + DensitemolaireKrViscoCompNode(icp,m,nums) * FluxDarcyFI(mph,s,k)

                      ! if(commRank==1 .and. k==1) then
                      !    print*, k, s, DensitemolaireKrViscoCompNode(icp,m,nums), FluxDarcyFI(mph,s,k)
                      ! end if

                   end if
                end do

#ifdef _THERMIQUE_

                FluxT_ks = FluxT_ks + DensitemolaireKrViscoEnthalpieNode(m,nums) &
                     * FluxDarcyFI(mph,s,k)
#endif
             end if
          end do

          ResiduFrac(1:NbComp,k) = ResiduFrac(1:NbComp,k) + Flux_ks(:)
          ResiduNode(1:NbComp,nums) = ResiduNode(1:NbComp,nums) - Flux_ks(:) 

#ifdef _THERMIQUE_

          ResiduFrac(NbComp+1,k) = ResiduFrac(NbComp+1,k) &
               + FluxT_ks &
               + FluxFourierFI(s,k)

          ResiduNode(NbComp+1,nums) = ResiduNode(NbComp+1,nums) &
               - FluxT_ks &
               - FluxFourierFI(s,k)
#endif

       end do ! end of node s in frac k

    end do ! end of frac k

    ! Injection well
    ! It is possible that one ghost well contains some own nodes,
    ! thus, when assembly of flux, we need to take into account
    ! not only own well, but also ghost well
    do k=1, NbWellInjLocal_Ncpus(commRank+1) ! injection well k

       qw = 0.d0

       ! nodes of well k
       do s=NodebyWellInjLocal%Pt(k)+1, NodebyWellInjLocal%Pt(k+1)
          nums = NodebyWellInjLocal%Num(s)

          Pws = PerfoWellInj(s)%Pression    ! P_{w,s}
          Tws = PerfoWellInj(s)%Temperature ! T_{w,s}
          Ps = IncNode(nums)%Pression       ! P_s
          Ts = IncNode(nums)%Temperature    ! T_s

          WIDws = NodeDatabyWellInjLocal%Val(s)%WID
          WIFws = NodeDatabyWellInjLocal%Val(s)%WIF

          Flux_ks(:) = 0.d0
          FluxT_ks = 0.d0

          ! Flux q_{w,s,i} and q_{w,s,e}
          Ps_Pws = Ps - Pws

          if(Ps_Pws < 0.d0) then ! < 0

             do icp=1, NbComp
                Flux_ks(icp) = DensitemolaireKrViscoCompWellInj(icp,s) * WIDws * Ps_Pws
             end do
#ifdef _THERMIQUE_
             FluxT_ks = DensitemolaireKrViscoEnthalpieWellInj(s) * WIDws * Ps_Pws
#endif
          end if

          ! node equation 
          ResiduNode(1:NbComp,nums) = ResiduNode(1:NbComp,nums) + Flux_ks(:)

#ifdef _THERMIQUE_
          ResiduNode(NbComp+1,nums) = ResiduNode(NbComp+1,nums) + FluxT_ks ! + WIFws * (Ts-Tws)
#endif
          ! qw
          do icp=1, NbComp
             qw = qw + Flux_ks(icp)
          end do
       end do

       ! inj well equation
       if( DataWellInjLocal(k)%IndWell == 'p') then                 
          ResiduWellInj(k) = DataWellInjLocal(k)%PressionMax - IncPressionWellInj(k)
       else if( DataWellInjLocal(k)%IndWell == 'f') then
          ResiduWellInj(k) = qw - DataWellInjLocal(k)%FlowrateImposed
       else
          print*, "indwell inj should be p or f"
       end if
          
    end do ! end of node s in injection well k

    ! if(commRank==0) then
    !    print*, "qw inj", qw
    ! end if
    
    ! Production well
    do k=1, NbWellProdLocal_Ncpus(commRank+1)
       
       qw = 0.d0

       ! nodes of well k
       do s=NodebyWellProdLocal%Pt(k)+1, NodebyWellProdLocal%Pt(k+1)
          nums = NodebyWellProdLocal%Num(s)

          Pws = PerfoWellProd(s)%Pression ! P_{w,s}
          Ps = IncNode(nums)%Pression     ! P_s
          Ps_Pws = Ps - Pws
          WIDws = NodeDatabyWellProdLocal%Val(s)%WID

          ! print*, "ps_pw: ", Ps_Pws
          
          Flux_ks(:) = 0.d0
          FluxT_ks = 0.d0

          do m=1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
             mph = NumPhasePresente_ctx(m,IncNode(nums)%ic)

             if(Ps_Pws > 0.d0) then

                do icp=1, NbComp
                   if(MCP(icp,mph)==1) then ! \cap P_i
                      Flux_ks(icp) = Flux_ks(icp) + DensitemolaireKrViscoCompNode(icp,m,nums) * WIDws * Ps_Pws
                   end if
                end do

#ifdef _THERMIQUE_
                FluxT_ks = FluxT_ks + DensitemolaireKrViscoEnthalpieNode(m,nums) * WIDws * Ps_Pws 
#endif
             end if
          end do

          ResiduNode(1:NbComp,nums) = ResiduNode(1:NbComp,nums) + Flux_ks(:)

#ifdef _THERMIQUE_
          ResiduNode(NbComp+1,nums) = ResiduNode(NbComp+1,nums) + FluxT_ks
#endif
          ! qw
          do icp=1, NbComp
             qw = qw + Flux_ks(icp)
          end do
       end do

       ! prod well equation
       if( DataWellProdLocal(k)%IndWell == 'p' ) then
          ResiduWellProd(k) = IncPressionWellProd(k) - DataWellProdLocal(k)%PressionMin
       else if ( DataWellProdLocal(k)%IndWell == 'f') then
          ResiduWellProd(k) = DataWellProdLocal(k)%FlowrateImposed - qw
       else
          print*, "indwell prod should be p or f"
       end if

       ! if(commRank==25) then
       !    print*, "qw prod", qw
       ! end if
       
    end do

    ! if(commRank==0) then
    !    print*, ResiduWellInj(:)
    !    print*, ResiduWellProd(:)
    ! end if

  end subroutine Residu_Conservations


  subroutine Residu_RelativeNorm(Iter, Delta_t, &
       ResNormRel, ResConvInit, ResClosInit)

    integer, intent(in) :: Iter ! Newton iteration: Iter
    double precision, intent(in) :: Delta_t
    double precision, intent(out) :: ResNormRel

    ! Residu relative norm from first Newton iteration
    ! used for all Newton iterations
    double precision, intent(inout) :: &
         ResConvInit(NbCompThermique), &
         ResClosInit

    integer :: k, i, ic, start, Ierr

    double precision :: &
         ResConvInitLocal(NbCompThermique), ResClosInitLocal

    double precision :: &
         ResConvRefLocal(NbCompThermique), ResClosRef

    double precision:: &
         ResConvLocal(NbCompThermique), ResClosLocal,  &
         ResNormRelLocal

    ! residu conservation: ResConv
    ResConvLocal(:) = 0.d0

    do k=1, NbNodeOwn_Ncpus(commRank+1) ! node
       ResConvLocal(:) = ResConvLocal(:) + abs(ResiduNode(:,k))
    enddo

    do k=1, NbFracOwn_Ncpus(commRank+1) ! frac
       ResConvLocal(:) = ResConvLocal(:) + abs(ResiduFrac(:,k))
    enddo

    do k=1, NbCellOwn_Ncpus(commRank+1) ! cell
       ResConvLocal(:) = ResConvLocal(:) + abs(ResiduCell(:,k))
    enddo

    do k=1, NbWellInjOwn_Ncpus(commRank+1) ! well inj
       ResConvLocal(1) = ResConvLocal(1) + abs(ResiduWellInj(k))
    enddo

    do k=1, NbWellProdOwn_Ncpus(commRank+1) ! well prod
       ResConvLocal(1) = ResConvLocal(1) + abs(ResiduWellProd(k))
    enddo
    
    ! residu closing equations: ResClos
    ResClosLocal = 0.d0

    do k=1, NbNodeOwn_Ncpus(commRank+1) ! node
       ic = IncNode(k)%ic 

       start = NbPhasePresente_ctx(ic)       
       do i=1, NbEqEquilibre_ctx(ic)
          ResClosLocal = ResClosLocal + abs(SmFNode(start+i,k))
       enddo
    enddo

    do k=1, NbFracOwn_Ncpus(commRank+1) ! frac
       ic = IncFrac(k)%ic

       start = NbPhasePresente_ctx(ic)       
       do i=1, NbEqEquilibre_ctx(ic)
          ResClosLocal = ResClosLocal + abs(SmFFrac(start+i,k))
       enddo
    enddo

    do k=1, NbCellOwn_Ncpus(commRank+1) ! cell
       ic = IncCell(k)%ic

       start = NbPhasePresente_ctx(ic)
       do i=1, NbEqEquilibre_ctx(ic)
          ResClosLocal = ResClosLocal + abs(SmFCell(start+i,k))
       enddo
    enddo

    ! If the first Newton iteration(Iter==1),
    ! compute references: ResConvInit and ResClosInit
    if(Iter==1) then

       ! compute a reference for relative residues
       ResConvRefLocal(:) = 1.d0
       ResClosInitLocal = 1.d0

       do k=1, NbNodeOwn_Ncpus(commRank+1)
          ResConvRefLocal(:) = ResConvRefLocal(:) + &
               abs(IncNode(k)%AccVol(:))/Delta_t/1000.d0
       end do
       do k=1, NbFracOwn_Ncpus(commRank+1)
          ResConvRefLocal(:) = ResConvRefLocal(:) + &
               abs(IncFrac(k)%AccVol(:))/Delta_t/1000.d0
       end do
       do k=1, NbCellOwn_Ncpus(commRank+1)
          ResConvRefLocal(:) = ResConvRefLocal(:) + &
               abs(IncCell(k)%AccVol(:))/Delta_t/1000.d0
       end do

       do i=1, NbCompThermique
          ResConvInitLocal(i) = max(ResConvLocal(i), ResConvRefLocal(i))
       end do
       ResClosInitLocal = max(ResClosLocal, ResClosInitLocal)

       ! if(commRank==0) then
          ! print*, ResConvRefLocal(:)
          ! print*, ResConvLocal(:)
          ! print*, ResConvInitLocal(:)
       ! end if
       
       ! max for all mpi
       call MPI_Allreduce(ResConvInitLocal, ResConvInit, NbCompThermique, &
            MPI_DOUBLE, MPI_MAX, ComPASS_COMM_WORLD, Ierr)

       call MPI_Allreduce(ResClosInitLocal, ResClosInit, 1, &
            MPI_DOUBLE, MPI_MAX, ComPASS_COMM_WORLD, Ierr)
    end if

    ! Relative residu
    ResNormRelLocal = ResClosLocal/ResClosInit
    do i=1, NbCompThermique
       ResNormRelLocal = max(ResConvLocal(i)/ResConvInit(i), ResNormRelLocal) 
    enddo

    ! max for all mpi
    call MPI_Allreduce(ResNormRelLocal, ResNormRel, 1, &
         MPI_DOUBLE, MPI_MAX, ComPASS_COMM_WORLD, Ierr)

  end subroutine Residu_RelativeNorm


  ! free
  subroutine Residu_free

    deallocate( ResiduCell)
    deallocate( ResiduFrac)
    deallocate( ResiduNode)

    deallocate( ResiduWellInj)
    deallocate( ResiduWellProd)

    deallocate( AccVolCell_1)
    deallocate( AccVolFrac_1)
    deallocate( AccVolNode_1)

  end subroutine Residu_free

end module Residu



! ! residu conservation energy
! subroutine Residu_ConservationsEnergy

!   integer :: k, s, fs, fk, nums, m, mph, icp
!   integer :: NbNodeCell, NbFracCell, NbNodeFrac

!   double precision :: Flux_ks


!   do k=1, NbCellLocal_Ncpus(commRank+1)

!      ! number of nodes/fracs in cell k
!      NbNodeCell = NodebyCellLocal%Pt(k+1) - NodebyCellLocal%Pt(k)
!      NbFracCell = FracbyCellLocal%Pt(k+1) - FracbyCellLocal%Pt(k)

!      ! s is node in cell k
!      do s=1, NbNodeCell
!         nums = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+s)

!         Flux_ks = 0.d0

!    ! alpha in Q_{K_{k,s}^{alpha}} =
!    !      ( {alpha | alpha in Q_k} \cap {alpha | K_{k,s}^{alpha}=k})
!    ! \cup ( {alpha | alpha in Q_s} \cap {alpha | K_{k,s}^{alpha}=s})
!    do m=1, NbPhasePresente_ctx(IncCell(k)%ic) ! Q_k
!       mph = NumPhasePresente_ctx(m, IncCell(k)%ic) 

!       if(FluxDarcyKI(mph,s,k)>=0.d0) then ! K_{k,s}^{alpha}=k

!          Flux_ks = Flux_ks + DensitemolaireKrViscoEnthalpieCell(m,k) &
!               * FluxDarcyKI(mph,s,k)
!       end if
!    end do

!    do m=1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
!       mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)

!       if(FluxDarcyKI(mph,s,k)<0.d0) then ! K_{k,s}^{alpha}=s             

!          Flux_ks = Flux_ks + DensitemolaireKrViscoEnthalpieNode(m,nums) &
!               * FluxDarcyKI(mph,s,k)
!       end if
!    end do

!    ResiduCell(NbComp+1,k) = ResiduCell(NbComp+1,k) &
!         + Flux_ks &
!         + FluxFourierKI(s,k)

!    ResiduNode(NbComp+1,nums) = ResiduNode(NbComp+1,nums) &
!         - Flux_ks &
!         - FluxFourierKI(s,k)

! end do ! end of node s in cell k


!    ! s is frac in cell k
!    do s=1, NbFracCell

!       fs = FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+s) 
!       nums = FaceToFracLocal(fs) ! fs is face number, nums is frac number

!       Flux_ks = 0.d0

!       ! alpha in Q_{K_{k,s}^{alpha}} =
!       !      ( {alpha | alpha in Q_k} \cap {alpha | K_{k,s}^{alpha}=k})
!       ! \cup ( {alpha | alpha in Q_s} \cap {alpha | K_{k,s}^{alpha}=s})
!       do m=1, NbPhasePresente_ctx(IncCell(k)%ic) ! Q_k
!          mph = NumPhasePresente_ctx(m, IncCell(k)%ic) 

!          if(FluxDarcyKI(mph,s+NbNodeCell,k)>=0.d0) then ! K_{k,s}^{alpha}=k

!             Flux_ks = Flux_ks + DensitemolaireKrViscoEnthalpieCell(m,k) &
!                  * FluxDarcyKI(mph,s+NbNodeCell,k)
!          end if
!       end do

!       do m=1, NbPhasePresente_ctx(IncFrac(nums)%ic) ! Q_s
!          mph = NumPhasePresente_ctx(m, IncFrac(nums)%ic)

!          if(FluxDarcyKI(mph,s+NbNodeCell,k)<0.d0) then ! K_{k,s}^{alpha}=s             
!             Flux_ks = Flux_ks + DensitemolaireKrViscoEnthalpieFrac(m,nums) &
!                  * FluxDarcyKI(mph,s+NbNodeCell,k)
!          end if
!       end do

!       ResiduCell(NbComp+1,k) = ResiduCell(NbComp+1,k) &
!            + Flux_ks &
!            + FluxFourierKI(s+NbNodeCell,k)

!       ResiduFrac(NbComp+1,nums) = ResiduFrac(NbComp+1,nums) &
!            - Flux_ks &
!            - FluxFourierKI(s+NbNodeCell,k)

!    end do ! end of frac s in cell k

! end do ! end of cell k


!   do k=1, NbFracLocal_Ncpus(commRank+1)

!      fk = FracToFaceLocal(k) ! fk is face number

!      ! number of nodes in a frac
!      NbNodeFrac = NodebyFaceLocal%Pt(fk+1) - NodebyFaceLocal%Pt(fk)

!      do s=1, NbNodeFrac
!         nums = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk)+s)

!         Flux_ks = 0.d0

!         ! alpha in Q_{K_{k,s}^{alpha}} =
!         !      ( {alpha | alpha in Q_k} \cap {alpha | K_{k,s}^{alpha}=k})
!         ! \cup ( {alpha | alpha in Q_s} \cap {alpha | K_{k,s}^{alpha}=s})
!         do m=1, NbPhasePresente_ctx(IncFrac(k)%ic) ! Q_k
!            mph = NumPhasePresente_ctx(m, IncFrac(k)%ic) 

!            if(FluxDarcyFI(mph,s,k)>=0.d0) then ! K_{k,s}^{alpha}=k

!               Flux_ks = Flux_ks + DensitemolaireKrViscoEnthalpieFrac(m,k) &
!                    * FluxDarcyFI(mph,s,k)
!            end if
!         end do

!         do m=1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
!            mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)

!            if(FluxDarcyFI(mph,s,k)<0.d0) then ! K_{k,s}^{alpha}=s             

!               Flux_ks = Flux_ks + DensitemolaireKrViscoEnthalpieNode(m,nums) &
!                    * FluxDarcyFI(mph,s,k)
!            end if
!         end do

!         ResiduFrac(NbComp+1,k) = ResiduFrac(NbComp+1,k) &
!              + Flux_ks &
!              + FluxFourierFI(s,k)

!         ResiduNode(NbComp+1,nums) = ResiduNode(NbComp+1,nums) &
!              - Flux_ks &
!              - FluxFourierFI(s,k)

!      end do ! end of node s in frac k

!   end do ! end of frac k

! end subroutine Residu_ConservationsEnergy