! Model: 2 phase 2 comp thermal, MCP=(1)
!        fractures z=0

! 1: GAS
! 2: WATER

! 1: Air
! 2: H2O

! *** Part of module IncCV.F90 *** !

! Set dir bc values: IncNodeDirBC
subroutine IncCV_SetDirBCValue

  integer :: i
  double precision :: sizeZ, sol

  double precision :: PgGal
  double precision :: PcGal
  double precision :: TGal
  double precision :: HurGal
  double precision :: SlGal
  double precision :: PsatGal, dT_PSatGal
  double precision :: ceg, cag
  double precision :: cel, cal
  double precision :: Rgp
    
  Rgp = 8.314d0

  PgGal = 1.d5
  TGal = 303.d0
  HurGal = 0.5d0

  CALL DefModel_Psat(TGal, PsatGal, dT_PSatGal)

  ceg = PsatGal * HurGal / PgGal
  cag = 1.d0 - ceg

  cal = cag * PgGal / HurGal
  cel = 1.d0 - cal

!  CALL f_Fugacity(PHASE_GAS,1,P,TGal,C,S,fa,DPf,DTf,DCf)
!  CALL f_Fugacity(PHASE_WATER,1,P,TGal,C,S,fe,DPf,DTf,DCf)

  CALL f_Sl(PcGal,SlGal)

  do i=1, NbNodeLocal_Ncpus(commRank+1)

    IF( IdNodeLocal(i)%P == "d") THEN

        if( abs(XNodeLocal(3,i)-Mesh_zmin)<eps ) then

          IncNodeDirBC(i)%ic = 3
          IncNodeDirBC(i)%Pression = PgGal
          IncNodeDirBC(i)%Saturation(PHASE_GAS) = 1.d0 - SlGal
          IncNodeDirBC(i)%Saturation(PHASE_WATER) = SlGal
          IncNodeDirBC(i)%Temperature = TGal
          IncNodeDirBC(i)%Comp(1,PHASE_GAS) = cag
          IncNodeDirBC(i)%Comp(2,PHASE_GAS) = ceg
          IncNodeDirBC(i)%Comp(1,PHASE_WATER) = cal
          IncNodeDirBC(i)%Comp(2,PHASE_WATER) = cel
          IncNodeDirBC(i)%AccVol(:) = 0.d0

        else if( abs(XNodeLocal(3,i)-Mesh_zmax)<eps ) then

          IncNodeDirBC(i)%ic = 2
          IncNodeDirBC(i)%Pression = 40.0d+5
          IncNodeDirBC(i)%Saturation(PHASE_GAS) = 0.d0
          IncNodeDirBC(i)%Saturation(PHASE_WATER) = 1.d0
          IncNodeDirBC(i)%Temperature = 303.d0

          ! # NbComp, NbPhase
          IncNodeDirBC(i)%Comp(1,PHASE_GAS) = 1.d0
          IncNodeDirBC(i)%Comp(2,PHASE_GAS) = 0.d0
          IncNodeDirBC(i)%Comp(1,PHASE_WATER) = 0.d0
          IncNodeDirBC(i)%Comp(2,PHASE_WATER) = 1.d0
          IncNodeDirBC(i)%AccVol(:) = 0.d0

          !IncNodeDirBC(i)%Temperature = 30.d0 + 273.d0   
        end if
    ENDIF
  end do

end subroutine IncCV_SetDirBCValue


subroutine IncCV_SetInitialValue

  integer :: i, s, m, nums, mph, fi, sparent
  double precision :: sizeZ
  double precision :: Tsat, dP_Tsat
  double precision :: Rhotmp, Pws, zp, zs, Pdrop
  double precision :: dPf, dTf, dCf(NbComp), dSf(NbPhase)

  sizeZ = meshSize_zmax - meshSize_zmin

  ! Node
  do i=1, NbNodeLocal_Ncpus(commRank+1)

    IncNode(i)%ic = 1
    IncNode(i)%Pression = 2.d5
    IncNode(i)%Saturation(1) = 1.d0
    IncNode(i)%Temperature = 20.d0 + 273.d0   
    IncNode(i)%Comp(1,1) = 1.d0
    IncNode(i)%AccVol(:) = 0.d0

  end do

  ! Cell
  do i=1, NbCellLocal_Ncpus(commRank+1)

    IncCell(i)%ic = 1
    IncCell(i)%Pression = 2.d5
    IncCell(i)%Saturation(1) = 1.d0
    IncCell(i)%Temperature = 20.d0 + 273.d0   
    IncCell(i)%Comp(1,1) = 1.d0
    IncCell(i)%AccVol(:) = 0.d0

  end do

end subroutine IncCV_SetInitialValue
