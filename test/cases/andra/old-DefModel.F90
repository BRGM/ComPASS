  !> \brief User set permeability
  !!
  !! \param[in] NbCellG,NbFracG Global number of cell and fracture face
  !! \param[in] IdCellG it is possible to set different permeability by rocktype
  !! \param[in,out] PermCellG Permeability tensor for each cell
  !! \param[in,out] PermFracG Permeability constant for each fracture face
  subroutine DefModel_SetPerm( &
      NbCellG, CellRocktypeG, NbFracG, FracRocktypeG, &
      PermCellG, PermFracG)

    integer, intent(in) :: NbCellG
    integer, dimension(:,:), intent(in) :: CellRocktypeG
    integer, intent(in) :: NbFracG
    integer, dimension(:,:), intent(in) :: FracRocktypeG
    ! ouptuts:
    double precision, dimension(:,:,:), allocatable, intent(inout) :: &
      PermCellG
    double precision, dimension(:), allocatable, intent(inout) :: &
      PermFracG

    integer :: i

    ! allocate and set values
    allocate(PermCellG(3,3,NbCellG))
    do i=1, NbCellG
      PermCellG(:,:,i) = 0.d0

     IF(CellRocktypeG(1,i) == 1) THEN
        PermCellG(1,1,i) = 5.d-20
        PermCellG(2,2,i) = 5.d-20
        PermCellG(3,3,i) = 5.d-20
     ELSEIF(CellRocktypeG(1,i) == 2)THEN
       PermCellG(1,1,i) = 1.d-18
       PermCellG(2,2,i) = 1.d-18
       PermCellG(3,3,i) = 1.d-18
     ELSE
       PRINT*, 'error DefModel_SetPerm, unknow rocktype'
       PRINT*, i, CellRocktypeG(1,i)
       STOP
     ENDIF
    end do

    allocate(PermFracG(NbFracG))
    PermFracG(:) = 1.d-11
  end subroutine DefModel_SetPerm


  subroutine DefModel_SetPorosite( &
      NbCellG, CellRocktypeG, NbFracG, FracRocktypeG, &
      PorositeCell, PorositeFrac)

    integer, intent(in) :: NbCellG
    integer, dimension(:,:), intent(in) :: CellRocktypeG
    integer, intent(in) :: NbFracG
    integer, dimension(:,:), intent(in) :: FracRocktypeG
    ! ouptuts:
    double precision, dimension(:), allocatable, intent(inout) :: &
      PorositeCell
    double precision, dimension(:), allocatable, intent(inout) :: &
      PorositeFrac

    integer :: i

    allocate(PorositeCell(NbCellG))
    do i=1, NbCellG
      IF(CellRocktypeG(1,i) == 1) THEN
        PorositeCell(i) = 0.15d0
      ELSEIF(CellRocktypeG(1,i) == 2)THEN
        PorositeCell(i) = 0.3d0
      ELSE
        PRINT*, 'error in DefModel_SetPorosite, unknow CellRocktype'
        PRINT*, i, CellRocktypeG(1,i)
        STOP
      ENDIF
    end do

    allocate(PorositeFrac(NbFracG))
    do i=1, NbFracG
      IF(FracRocktypeG(1,i) == 1) THEN
        PorositeFrac(i) = 0.15d0
      ELSEIF(FracRocktypeG(1,i) == 2)THEN
        PorositeFrac(i) = 0.3d0
      ELSE
        PRINT*, 'error in DefModel_SetPorosite, unknow FracRocktype'
        PRINT*, i, FracRocktypeG(1,i)
        STOP
      ENDIF
    ENDDO
  end subroutine DefModel_SetPorosite


#ifdef _THERMIQUE_

  subroutine DefModel_SetCondThermique( &
      NbCellG, CellRocktypeG, NbFracG, FracRocktypeG, &
      CondThermalCellG, CondThermalFracG)

    integer, intent(in) :: NbCellG
    integer, dimension(:,:), intent(in) :: CellRocktypeG
    integer, intent(in) :: NbFracG
    integer, dimension(:,:), intent(in) :: FracRocktypeG

    ! ouptuts:
    double precision, dimension(:,:,:), allocatable, intent(inout) :: &
      CondThermalCellG
    double precision, dimension(:), allocatable, intent(inout) :: &
      CondThermalFracG

    integer :: i

    allocate(CondThermalCellG(3,3,NbCellG))
    do i=1, NbCellG
      CondThermalCellG(:,:,i) = 0.d0
      IF(CellRocktypeG(2,i) == 1) THEN
        CondThermalCellG(1,1,i) = 2.d0
        CondThermalCellG(2,2,i) = 2.d0
        CondThermalCellG(3,3,i) = 2.d0
      ELSEIF(CellRocktypeG(2,i) == 2)THEN
        CondThermalCellG(1,1,i) = 2.d0
        CondThermalCellG(2,2,i) = 2.d0
        CondThermalCellG(3,3,i) = 2.d0
      ELSE
        PRINT*, 'error in DefModel_SetCondThermique, unknow FracRocktype'
        PRINT*, i, FracRocktypeG(2,i)
        STOP
      ENDIF
    end do

    allocate(CondThermalFracG(NbFracG))
    CondThermalFracG(:) = 2.d0

  end subroutine DefModel_SetCondThermique


  !SUBROUTINE DefModel_SetThermalSource( &
  !    NbCell, &
  !    CellThermalSourceType, &
  !    NbFrac, &
  !    FracThermalSourceType, &
  !    CellThermalSource, &
  !    FracThermalSource)
  !
  !  INTEGER, INTENT(IN) :: NbCell
  !  INTEGER, DIMENSION(:), INTENT(IN) :: CellThermalSourceType
  !  INTEGER, INTENT(IN) :: NbFrac
  !  INTEGER, DIMENSION(:), INTENT(IN) :: FracThermalSourceType
  !  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: CellThermalSource
  !  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: FracThermalSource
  !
  !  ALLOCATE(CellThermalSource(NbCell))
  !  CellThermalSource = MERGE(30.d0, 0.d0, CellThermalSourceType == 1)
  !
  !  ALLOCATE(FracThermalSource(NbFrac))
  !  FracThermalSource = 0.d0
  !END SUBROUTINE DefModel_SetThermalSource


