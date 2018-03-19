!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

! Model: 2 phase 2 comp thermal, MCP=(1)
!        fractures z=0

! 1: GAS
! 2: WATER

! 1: Air
! 2: H2O

! *** Part of module IncCV.F90 *** !

! Set dir bc values: IncNodeDirBC
subroutine IncCV_SetDirBCValue

  integer :: i, rt(IndThermique+1)
  double precision :: sizeZ, sol

  integer :: icPor
  integer :: rtPor(IndThermique+1)
  double precision :: PlPor, PlTop, PlRef, zRef
  double precision :: SPor(NbPhase)
  double precision :: PcPor, PgPor, PPor
  double precision :: TPor
  double precision :: HurPor
  double precision :: PsatPor, dT_PSatPor
  double precision :: CegPor, CagPor
  double precision :: CelPor, CalPor
  double precision :: SlPor, SgPor
  double precision :: CPor(NbComp)
  double precision :: rho, dPf, dTf, dCf(NbComp)
  integer :: icGal
  integer :: rtGal(IndThermique+1)
  double precision :: PcGal, PgGal, PGal
  double precision :: TGal
  double precision :: HurGal
  double precision :: PsatGal, dT_PSatGal,Psatt
  double precision :: CegGal, CagGal
  double precision :: CelGal, CalGal
  double precision :: SlGal, SgGal
  double precision :: SGal(NbPhase)
  double precision :: RZetal
  double precision :: DSf(NbPhase)
  double precision :: Ha
  double precision :: Pc

  icPor = 2

  PlTop = 4.0d+6
  SgPor = 0.d0
  SlPor = 1.d0
  TPor = 303.d0

  SPor = (/ SgPor, SlPor /)

  CagPor = 0.d0
  CegPor = 1.d0 - CagPor

  CalPor = 0.d0
  CelPor = 1.d0

  CPor = (/ CalPor, CelPor /)

  zRef = meshSizeS_zmin
  PlRef = PlTop
  CALL f_DensiteMassique(PHASE_WATER,PlRef,TPor,CPor,SPor,rho,dPf,dTf,dCf,dSf)

!   !!!!!!!!!!!!!!!!!!!!!!!
!
!   icPor=3
!
!   PgPor = 1.d5
!   TPor = 303.d0
!   HurPor = 0.9d0
!
!   PPor = PgPor
!
!   CALL DefModel_Psat(TPor, PsatPor, dT_PSatPor)
!
!   CegPor = PsatPor*HurPor/PgPor
!   CagPor = 1.d0 - CegPor
!
!   CALL air_henry(TPor,Ha)
!   CalPor = CagPor*PgPor/Ha
!   CelPor = 1.d0 - CalPor
!
!   RZetal = 8.314d0 * 1000.d0 / 0.018d0
!   PcPor = DLOG(CelPor/HurPor) * RZetal * TPor
!
!   CALL f_Sl(PcPor,SlPor)
!   SgPor = 1 - SlPor
!
!  !!!!!!!!!!!!!!!!!!!!!!!

  icGal=3

!  PgGal = 1.d5
!  TGal = 303.d0
!  SlGal = 0.5d0
!  SgGal = 1 - SlGal
!  PGal = PgGal
!  SGal = (/ SgGal, SlGal /)
!  CALL f_PressionCapillaire(rt,2,SGal,PcGal,DSf)
!  CALL DefModel_Psat(TGal, PsatGal, dT_PSatGal)
!  RZetal = 8.314d0 * 1000.d0 / 0.018d0
!  Psatt =  PsatGal*dexp(PcGal/RZetal / TGal)
!  CALL air_henry(TGal,Ha)
!  CelGal = (Ha-PgGal)/(Ha-Psatt)
!  CalGal = 1.d0 - CelGal
!  CagGal = Ha*CalGal/PgGal
!  CegGal = 1.d0 - CagGal

  PgGal = 1.d5
  TGal = 303.d0
  HurGal = 0.5d0
  PGal = PgGal

  CALL DefModel_Psat(TGal, PsatGal, dT_PSatGal)

  CegGal = PsatGal*HurGal/PgGal
  CagGal = 1.d0 - CegGal

  CALL air_henry(TGal,Ha)
  CalGal = CagGal*PgGal/Ha
  CelGal = 1.d0 - CalGal

  RZetal = 8.314d0 * 1000.d0 / 0.018d0
  PcGal = DLOG(CelGal/HurGal) * RZetal * TGal

  do i=1, NbNodeLocal_Ncpus(commRank+1)

    IF( IdNodeLocal(i)%P == "d") THEN

        if( NodeFlagsLocal(i) == 1 ) then

          rtGal = NodeRocktypeLocal(:,i)
          CALL f_Sl(rtGal,PcGal,SlGal)
          SgGal = 1 - SlGal

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

        else if( NodeFlagsLocal(i) == 2 )then

          PlPor = liquid_pressure(zRef, PlRef, rho, gravity, XNodeLocal(3,i))

          rtPor = NodeRocktypeLocal(:,i)
          CALL f_PressionCapillaire(rtPor,2,SPor,PcPor,DSf)

          PPor = PlPor - PcPor
          PgPor = PPor

        IncNodeDirBC(i)%ic = icPor
        IncNodeDirBC(i)%Pression = PPor
        IncNodeDirBC(i)%Saturation(PHASE_GAS) = SgPor
        IncNodeDirBC(i)%Saturation(PHASE_WATER) = SlPor
        IncNodeDirBC(i)%Temperature = TPor
        IncNodeDirBC(i)%Comp(1,PHASE_GAS) = CagPor
        IncNodeDirBC(i)%Comp(2,PHASE_GAS) = CegPor
        IncNodeDirBC(i)%Comp(1,PHASE_WATER) = CalPor
        IncNodeDirBC(i)%Comp(2,PHASE_WATER) = CelPor
          IncNodeDirBC(i)%AccVol(:) = 0.d0

        ELSE
          PRINT*, 'ERROR STOP'
          STOP
        end if
    ENDIF
  end do

end subroutine IncCV_SetDirBCValue


subroutine IncCV_SetInitialValue

  integer :: i, rt(IndThermique+1)
  double precision :: sizeZ, sol

  integer :: icPor
  double precision :: PlPor, PlTop, PlRef, zRef
  double precision :: PgPor, PcPor, PPor
  double precision :: TPor
  double precision :: HurPor
  double precision :: PsatPor, dT_PSatPor
  double precision :: CegPor, CagPor
  double precision :: CelPor, CalPor
  double precision :: RZetal
  double precision :: SlPor, SgPor, SPor(NbPhase)
  double precision :: Ha
  double precision :: DSf(NbPhase)
  double precision :: CPor(NbComp)
  double precision :: rho, dPf, dTf, dCf(NbComp)

  PlPor = 40.0d+5
  icPor = 2
  SgPor = 0.d0
  SlPor = 1.d0
  TPor = 303.d0

  SPor = (/ SgPor, SlPor /)

  CalPor = 0.d0
  CelPor = 1.d0 - CalPor

  !!!!!!!!!!!!!!!!!!!

  icPor = 2

  PlTop = 4.0d+6
  SgPor = 0.d0
  SlPor = 1.d0
  TPor = 303.d0

  SPor = (/ SgPor, SlPor /)

  CagPor = 0.d0
  CegPor = 1.d0 - CagPor

  CalPor = 0.d0
  CelPor = 1.d0

  CPor = (/ CalPor, CelPor /)

  zRef = meshSizeS_zmin
  PlRef = PlTop
  CALL f_DensiteMassique(PHASE_WATER,PlRef,TPor,CPor,SPor,rho,dPf,dTf,dCf,dSf)

!  !!!!!!!!!!!!!!!!!!!
!
!  icPor=3
!
!  PgPor = 1.d5
!  TPor = 303.d0
!  HurPor = 0.9d0
!
!  PPor = PgPor
!
!  CALL DefModel_Psat(TPor, PsatPor, dT_PSatPor)
!
!  CegPor = PsatPor*HurPor/PgPor
!  CagPor = 1.d0 - CegPor
!
!  CALL air_henry(TPor,Ha)
!  CalPor = CagPor*PgPor/Ha
!  CelPor = 1.d0 - CalPor
!
!  RZetal = 8.314d0 * 1000.d0 / 0.018d0
!  PcPor = DLOG(CelPor/HurPor) * RZetal * TPor
!
!  CALL f_Sl(PcPor,SlPor)
!  SgPor = 1 - SlPor

  PlPor = PlTop

  ! Node
  do i=1, NbNodeLocal_Ncpus(commRank+1)

    rt = NodeRocktypeLocal(:,i)
    CALL f_PressionCapillaire(rt,PHASE_WATER,SPor,PcPor,DSf)

    PPor = PlPor - PcPor
    PgPor = PPor

    IncNode(i)%ic = icPor
    IncNode(i)%Pression = PPor
    IncNode(i)%Saturation(PHASE_GAS) = SgPor
    IncNode(i)%Saturation(PHASE_WATER) = SlPor
    IncNode(i)%Temperature = TPor
    IncNode(i)%Comp(1,PHASE_GAS) = CagPor
    IncNode(i)%Comp(2,PHASE_GAS) = CegPor
    IncNode(i)%Comp(1,PHASE_WATER) = CalPor
    IncNode(i)%Comp(2,PHASE_WATER) = CelPor
    IncNode(i)%AccVol(:) = 0.d0

  end do

  ! Cell
  do i=1, NbCellLocal_Ncpus(commRank+1)

    rt = CellRocktypeLocal(:,i)
    CALL f_PressionCapillaire(rt,PHASE_WATER,SPor,PcPor,DSf)

    PPor = PlPor - PcPor
    PgPor = PPor


    IncCell(i)%ic = icPor
    IncCell(i)%Pression = PPor
    IncCell(i)%Saturation(PHASE_GAS) = SgPor
    IncCell(i)%Saturation(PHASE_WATER) = SlPor
    IncCell(i)%Temperature = TPor
    IncCell(i)%Comp(1,PHASE_GAS) = CagPor
    IncCell(i)%Comp(2,PHASE_GAS) = CegPor
    IncCell(i)%Comp(1,PHASE_WATER) = CalPor
    IncCell(i)%Comp(2,PHASE_WATER) = CelPor
    IncCell(i)%AccVol(:) = 0.d0

  end do

end subroutine IncCV_SetInitialValue
