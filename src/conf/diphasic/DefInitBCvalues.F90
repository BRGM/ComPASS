! Model: 2 phase 2 comp thermal, MCP=(1)
!        fractures z=0

! 1: GAS
! 2: WATER

! 1: Air
! 2: H2O

! *** Part of module IncCV.F90 *** !

! Set dir bc values: IncNodeDirBC
subroutine IncCV_SetDirBCValue

  integer :: i, icPor, icGal, rocktype
  double precision :: sizeZ, sol

  double precision :: PgPor, PlPor, PcPor, PPor
  double precision :: TPor
  double precision :: SgPor, SlPor, SPor(NbPhase)
  double precision :: DSf(NbPhase)
  double precision :: PcGal, PgGal, PGal
  double precision :: TGal
  double precision :: HurGal
  double precision :: PsatGal, dT_PSatGal
  double precision :: CegGal, CagGal
  double precision :: CelGal, CalGal
  double precision :: RZetal
  double precision :: SlGal, SgGal
  double precision :: Ha

  rocktype = 1

  PlPor = 40.0d+5
  icPor = 2
  SgPor = 0.d0
  SlPor = 1.d0
  TPor = 303.d0

  SPor = (/ SgPor, SlPor /)
  CALL f_PressionCapillaire(rocktype,2,SPor,PcPor,DSf)
  PPor = PlPor - PcPor
    
  PgGal = 2.d5
  TGal = 303.d0
  HurGal = 0.5d0

  PGal = PgGal

  CALL DefModel_Psat(TGal, PsatGal, dT_PSatGal)

  CegGal = MIN(MAX(PsatGal*HurGal/PgGal,0.d0),1.d0)
  CagGal = 1.d0 - CegGal

  RZetal = 8.314d0 * 1000.d0 / 0.018d0
  PcGal = LOG(CegGal*PgGal/PsatGal) * RZetal * TGal

  IF(PcGal < 0)THEN
    icGal=2
    SlGal = 1.d0
    SgGal = 0.d0
  ELSE
    icGal=3
    CALL f_Sl(PcGal,SlGal)
    SgGal = 1 - SlGal

    CALL air_henry(TGal,Ha)
    CalGal = MIN(MAX(CagGal*PgGal/Ha,0.d0),1.d0)
    CelGal = 1.d0 - CalGal
  ENDIF

  do i=1, NbNodeLocal_Ncpus(commRank+1)

    IF( IdNodeLocal(i)%P == "d") THEN

        if( abs(XNodeLocal(3,i)-Mesh_zmin)<eps ) then

          IncNodeDirBC(i)%ic = icGal
          IncNodeDirBC(i)%Pression = PGal
          IncNodeDirBC(i)%Saturation(PHASE_GAS) = SgGal
          IncNodeDirBC(i)%Saturation(PHASE_WATER) = SlGal
          IncNodeDirBC(i)%Temperature = TGal
          IncNodeDirBC(i)%Comp(1,PHASE_GAS) = CagGal
          IncNodeDirBC(i)%Comp(2,PHASE_GAS) = CegGal
          IncNodeDirBC(i)%Comp(1,PHASE_WATER) = CalGal
          IncNodeDirBC(i)%Comp(2,PHASE_WATER) = CelGal
          IncNodeDirBC(i)%AccVol(:) = 0.d0

        else if( abs(XNodeLocal(3,i)-Mesh_zmax)<eps ) then

          IncNodeDirBC(i)%ic = icPor
          IncNodeDirBC(i)%Pression = PPor
          IncNodeDirBC(i)%Temperature = TPor
          IncNodeDirBC(i)%Saturation(PHASE_GAS) = SgPor
          IncNodeDirBC(i)%Saturation(PHASE_WATER) = SlPor

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

  integer :: i, icPor, rocktype
  double precision :: sizeZ, sol

  double precision :: PgPor, PlPor, PcPor, PPor
  double precision :: TPor
  double precision :: SgPor, SlPor, SPor(NbPhase)
  double precision :: DSf(NbPhase)

  rocktype = 1

  PlPor = 40.0d+5
  icPor = 2
  SgPor = 0.d0
  SlPor = 1.d0
  TPor = 303.d0

  SPor = (/ SgPor, SlPor /)
  CALL f_PressionCapillaire(rocktype,2,SPor,PcPor,DSf)
  PPor = PlPor - PcPor

  ! Node
  do i=1, NbNodeLocal_Ncpus(commRank+1)

    IncNode(i)%ic = icPor
    IncNode(i)%Pression = PPor
    IncNode(i)%Saturation(PHASE_GAS) = SgPor
    IncNode(i)%Saturation(PHASE_WATER) = SlPor
    IncNode(i)%Temperature = TPor
    IncNode(i)%Comp(1,PHASE_GAS) = 1.d0
    IncNode(i)%Comp(2,PHASE_GAS) = 0.d0
    IncNode(i)%Comp(1,PHASE_WATER) = 0.d0
    IncNode(i)%Comp(2,PHASE_WATER) = 1.d0
    IncNode(i)%AccVol(:) = 0.d0

  end do

  ! Cell
  do i=1, NbCellLocal_Ncpus(commRank+1)

    IncCell(i)%ic = icPor
    IncCell(i)%Pression = PPor
    IncCell(i)%Saturation(PHASE_GAS) = SgPor
    IncCell(i)%Saturation(PHASE_WATER) = SlPor
    IncCell(i)%Temperature = TPor
    IncCell(i)%Comp(1,PHASE_GAS) = 1.d0
    IncCell(i)%Comp(2,PHASE_GAS) = 0.d0
    IncCell(i)%Comp(1,PHASE_WATER) = 0.d0
    IncCell(i)%Comp(2,PHASE_WATER) = 1.d0
    IncCell(i)%AccVol(:) = 0.d0

  end do

end subroutine IncCV_SetInitialValue
