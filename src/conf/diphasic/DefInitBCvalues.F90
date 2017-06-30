! Model: 1 phase 1 comp thermal, MCP=(1)
!        fractures z=0

! 1: Water

! *** Part of module IncCV.F90 *** !

! Set dir bc values: IncNodeDirBC
subroutine IncCV_SetDirBCValue

  integer :: i
  double precision :: sizeZ, sol

  do i=1, NbNodeLocal_Ncpus(commRank+1)

    IF( IdNodeLocal(i)%P == "d") THEN

        if( abs(XNodeLocal(3,i)-Mesh_zmin)<eps ) then

          IncNodeDirBC(i)%ic = 1
          IncNodeDirBC(i)%Pression = 2.d5
          IncNodeDirBC(i)%Saturation(1) = 1.d0
          IncNodeDirBC(i)%Temperature = 20.d0 + 273.d0   
          IncNodeDirBC(i)%Comp(1,1) = 1.d0
          IncNodeDirBC(i)%AccVol(:) = 0.d0

        else if( abs(XNodeLocal(3,i)-Mesh_zmax)<eps ) then

          IncNodeDirBC(i)%ic = 1
          IncNodeDirBC(i)%Pression = 2.d7
          IncNodeDirBC(i)%Saturation(1) = 1.d0
          IncNodeDirBC(i)%Temperature = 30.d0 + 273.d0   
          IncNodeDirBC(i)%Comp(1,1) = 1.d0
          IncNodeDirBC(i)%AccVol(:) = 0.d0

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
