!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

! Model: 2 phase 1 comp thermal, MCP=(1,1)
!        fractures z=0

! 1: Gas
! 2: Water

! *** Part of module IncCV.F90 *** !

! Set dir bc values: IncNodeDirBC
subroutine IncCV_SetDirBCValue

  integer :: i
  double precision :: sizeZ, sol

  do i=1, NbNodeLocal_Ncpus(commRank+1)

     IncNodeDirBC(i)%ic = 2

     IncNodeDirBC(i)%Pression = 2.d7

     IncNodeDirBC(i)%Temperature = 140.d0 + 273.d0
     IncNodeDirBC(i)%Saturation(PHASE_WATER) = 1.d0
     IncNodeDirBC(i)%Saturation(PHASE_GAS) = 0.d0

     IncNodeDirBC(i)%Comp(1,1) = 1.d0
     IncNodeDirBC(i)%Comp(1,2) = 1.d0

     IncNodeDirBC(i)%AccVol(:) = 0.d0
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


     IncNode(i)%ic = 2

     IncNode(i)%Pression = 2.d7
     IncNode(i)%Temperature = 140.d0 + 273.d0

     IncNode(i)%Saturation(PHASE_WATER) = 1.d0
     IncNode(i)%Saturation(PHASE_GAS) = 1.d0

     IncNode(i)%Comp(1,1) = 1.d0
     IncNode(i)%Comp(1,2) = 1.d0

     IncNode(i)%AccVol(:) = 0.d0

  end do

  ! Cell
  do i=1, NbCellLocal_Ncpus(commRank+1)

     IncCell(i)%Pression = 2.d7
     IncCell(i)%Temperature = 140.d0 + 273.d0

     IncCell(i)%ic = 2
     IncCell(i)%Saturation(PHASE_WATER) = 1.d0
     IncCell(i)%Saturation(PHASE_GAS) = 0.d0

     IncCell(i)%Comp(1,1) = 1.d0
     IncCell(i)%Comp(1,2) = 1.d0

     IncCell(i)%AccVol(:) = 0.d0
  end do

  ! Frac
  do i=1, NbFracLocal_Ncpus(commRank+1)

     IncFrac(i)%Pression = 2.d7
     IncFrac(i)%Temperature = 140.d0 + 273.d0

     IncFrac(i)%ic = 2
     IncFrac(i)%Saturation(PHASE_WATER) = 1.d0
     IncFrac(i)%Saturation(PHASE_GAS) = 0.d0

     IncFrac(i)%Comp(1,1) = 1.d0
     IncFrac(i)%Comp(1,2) = 1.d0

     IncFrac(i)%AccVol(:) = 0.d0
  end do

  ! cas test monophasique liquide non isotherme avec puits injecteur et producteur et une faille


  ! Init production well
  do i=1, NbWellInjLocal_Ncpus(commRank+1)

     ! Pressure
     IncPressionWellInj(i) = 2.d7

     ! Temperature
     ! for the injector, the temperature is fixed all along the well (Tw)
     do s=NodebyWellInjLocal%Pt(i)+1, NodebyWellInjLocal%Pt(i+1)
        PerfoWellInj(s)%Temperature = DataWellInjLocal(i)%Temperature
     end do
  end do

  ! Compute PerfoWellInj%Pressure, %PressureDrop with Pw
  call IncCVWells_PressureDropWellInj

  do i=1, NbWellProdLocal_Ncpus(commRank+1)

     ! Pressure
     IncPressionWellProd(i) = 2.d7

     ! Temperature
     do s=NodebyWellProdLocal%Pt(i+1), NodebyWellProdLocal%Pt(i)+1, -1
        nums = NodebyWellProdLocal%Num(s)
        PerfoWellProd(s)%Temperature = IncNode(nums)%Temperature
     end do
  end do

  ! Compute PerfoWellProd%Pressure, %PressureDrop, %Density with Pw
  call IncCVWells_PressureDropWellProd

end subroutine IncCV_SetInitialValue
