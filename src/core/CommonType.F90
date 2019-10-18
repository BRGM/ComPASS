!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module CommonType

  ! This for arrays that are interfaced with python/C++
  use iso_c_binding, only: &
     c_int, c_int8_t, c_char, c_double

  use CommonMPI, only: &
     CommonMPI_abort

  implicit none

  !> Type to be used to pass model configuration when no performance is needed
  type ModelConfiguration
    integer(c_int) :: nb_phases, nb_components, nb_contexts
    integer(c_int) :: IndThermique, NbEqEquilibreMax, NbIncPTCMax
    integer(c_int), allocatable :: NbPhasePresente_ctx(:)
    integer(c_int), allocatable :: NumPhasePresente_ctx(:,:)
    integer(c_int), allocatable :: MCP(:,:)
    integer(c_int), allocatable :: pssecd(:,:)
  end type ModelConfiguration

  !> Array 1d integer
  type ARRAY1Int
    integer(c_int), allocatable, dimension(:) :: Val
  end type ARRAY1Int

  !> Array 1d integer
  type ARRAY1Int8
    integer(c_int8_t), allocatable, dimension(:) :: Val
  end type ARRAY1Int8

  !> Array 1d double precision
  type ARRAY1dble
    double precision, allocatable, dimension(:) :: Val
  end type ARRAY1dble

  !> Array 2d integer
  type ARRAY2Int
    integer(c_int), allocatable, dimension(:,:) :: Array2d
  end type ARRAY2Int

  !> Array 2d double precision
  type ARRAY2dble
    real(c_double), allocatable, dimension(:,:) :: Array2d
  end type ARRAY2dble

  !> Array 3d double precision
  type ARRAY3dble
    double precision, allocatable, dimension(:,:,:) :: Array3d
  end type ARRAY3dble

  !> Standard type CSR with Pt, Num, and Val (1d integer)
    ! kind=c_int is because csr are to be interfaced with python/C++
  type CSR
    integer(kind=c_int) :: Nb
    integer(kind=c_int), allocatable, dimension(:) :: Pt
    integer(kind=c_int), allocatable, dimension(:) :: Num
    integer(kind=c_int), allocatable, dimension(:) :: Val
  end type CSR

  type COC
    integer(kind=c_int)                            :: Nb
    integer(kind=c_int), allocatable, dimension(:) :: Pt
    integer(kind=c_int), allocatable, dimension(:) :: Num
  end type COC

  type FractureInfo
    integer(kind=c_int)                            :: face
    integer(kind=c_int)                            :: fracture
  end type FractureInfo

  type FractureInfoCOC
    integer(kind=c_int)                            :: Nb
    integer(kind=c_int), allocatable, dimension(:) :: Pt
    type(FractureInfo), allocatable, dimension(:)  :: Num
  end type FractureInfoCOC

  !> Standar type CSR with Pt, Num, and Val (1d double precision)
  type CSRdble
    integer :: Nb
    integer, allocatable, dimension(:) :: Pt
    integer, allocatable, dimension(:) :: Num
    double precision, allocatable, dimension(:) :: Val
  end type CSRdble

  !> Standar type CSR with Pt, Num, and Val (3d double precision)
  type CSRArray2dble
    integer :: Nb
    integer, allocatable, dimension(:) :: Pt
    integer, allocatable, dimension(:) :: Num
    double precision, allocatable, dimension(:,:,:) :: val
  end type CSRArray2dble

  !> Store data of Node about own/ghost; fractures and boundary conditions
  ! FIXME: P stands for "mass balance equationS"
  ! FIXME: T stands for "energy balance equation"
  type Type_IdNode
     sequence
     character(c_char) :: Proc  !< "o"/"g": own/ghost
     character(c_char) :: Frac  !< "y"/"n": node in fracture/not in fracture
     character(c_char) :: P     !< "d"/"n"/"i": dirichlet/newmann/interior for the Pressure
     character(c_char) :: T     !< "d"/"n"/"i"/"o": dirichlet/newmann/interior/outflow for the Temperature
  end type Type_IdNode

  !> Array 1d type_IdNode
  type Array1IdNode
    type(Type_IdNode), allocatable, dimension(:) :: val
  end type Array1IdNode

  ! tmp constant values
  integer, parameter, public :: &
       VALSIZE_ZERO = 0, & !< Identifier used in communication, val is not used in CSR
       VALSIZE_NB   = 1, & !< Identifier used in communication, val is used in CSR and size(%Val)=%Nb
       VALSIZE_NNZ  = 2    !< Identifier used in communication, val is used in CSR and size(%Val)=Nnz=size(%Num)

  public :: &
       CommonType_csrcopy, &            ! copy CSR
       CommonType_deallocCSR, &         ! free CSR
       CommonType_printCSR, &           ! print CSR (\%Val is integer)
       CommonType_printCSRdble          ! print CSR (\%Val is double)

  !> to allow = between two TYPE_IdNode
  interface assignment(=)
     module procedure assign_IdNode_equal
  end interface assignment(=)

    interface first_of_all
      module procedure CSR_first_of_all
    end interface first_of_all

    interface last_of_all
      module procedure CSR_last_of_all
    end interface last_of_all

    interface associate
      module procedure CSR_associate_row
    end interface associate

    interface row_view
      module procedure CSR_row_view
    end interface row_view

contains

  !> \brief Copy CSR1 to CSR2
  subroutine CommonType_csrcopy(CSR1, CSR2, valsize)

    type(CSR), intent(in)  :: CSR1
    type(CSR), intent(out) :: CSR2
    integer,   intent(in)  :: valsize

    integer :: Nnz, i

    CSR2%Nb = CSR1%Nb

    allocate(CSR2%Pt(CSR1%Nb+1))
    CSR2%Pt(:) = CSR1%Pt(:)

    Nnz = CSR1%Pt(CSR1%Nb+1)

    allocate( CSR2%Num(Nnz))
    do i=1, Nnz
      CSR2%Num(i) = CSR1%Num(i)
    end do

    if(valsize==VALSIZE_NB) then
      allocate( CSR2%Val( CSR1%Nb))
      CSR2%Val(:) = CSR1%Val(:)
    else if (valsize==VALSIZE_NNZ) then
      allocate( CSR2%Val(Nnz))
      do i=1, Nnz
        CSR2%Val(i) = CSR1%Val(i)
      end do
    end if

  end subroutine CommonType_csrcopy

  !> \brief Deallocate CSR (\%Pt, \%Num and \%Val)
  subroutine CommonType_deallocCSR(CSR1)

    type(CSR), intent(inout) :: CSR1

    deallocate(CSR1%Pt)

    if(allocated(CSR1%Num)) then
      deallocate(CSR1%Num)
    end if

    if (allocated(CSR1%Val)) then
      deallocate(CSR1%Val)
    end if

  end subroutine CommonType_deallocCSR

  subroutine CommonType_deallocCOC(COC1)

    type(COC), intent(inout) :: COC1

    deallocate(COC1%Pt)

    if(allocated(COC1%Num)) then
      deallocate(COC1%Num)
    end if

  end subroutine CommonType_deallocCOC

  subroutine CommonType_deallocFracInfoCOC(COC)

    type(FractureInfoCOC), intent(inout) :: COC

    if(allocated(COC%Pt)) then
      deallocate(COC%Pt)
    end if
    if(allocated(COC%Num)) then
      deallocate(COC%Num)
    end if

  end subroutine CommonType_deallocFracInfoCOC

  !> \brief Deallocate CSRdble (\%Pt, \%Num and \%Val)
  subroutine CommonType_deallocCSRdble(CSR1)

    type(CSRdble), intent(inout) :: CSR1

    deallocate(CSR1%Pt,CSR1%Num)
    if (allocated(CSR1%Val)) then
      deallocate(CSR1%Val)
    end if

  end subroutine CommonType_deallocCSRdble

  !> Print CSR
  subroutine CommonType_printCSR(CSR1)

    type(CSR), intent(in) :: CSR1
    integer :: i, j

    write(*,'(A4,I5)') "Nb:  ", CSR1%Nb
    do i=1,CSR1%Nb
      write(*,"(A4,I5,A8,I3)") "Row  ", i, " : size=", CSR1%Pt(i+1)-CSR1%Pt(i)

      do j=CSR1%Pt(i)+1,CSR1%Pt(i+1)
        write(*,"(I5)",advance="no") CSR1%Num(j)
      end do
      print*,""
    end do

  end subroutine CommonType_printCSR


  !> Print CSRdble
  subroutine CommonType_printCSRdble(CSR1)

    type(CSRdble), intent(inout) :: CSR1
    integer :: i, j

    write(*,'(A4,I3)') "Nb:  ", CSR1%Nb
    do i=1,CSR1%Nb
      write(*,"(A4,I3,A8,I2)") "Row  ", i, " : size=", CSR1%Pt(i+1)-CSR1%Pt(i)

      do j=CSR1%Pt(i)+1,CSR1%Pt(i+1)
        write(*,"(I3)",advance="no") CSR1%Num(j)
      end do
      print*,""
    end do

    print*, ""
    print*, ""

    do i=1,CSR1%Nb
      write(*,"(A4,I3,A8,I2)") "Row  ", i, " : size=", CSR1%Pt(i+1)-CSR1%Pt(i)

      do j=CSR1%Pt(i)+1,CSR1%Pt(i+1)
        ! write(*,"(D2.6)",advance="no"), CSR1%Val(j)
        print*, CSR1%Val(j)
      end do
      print*,""
    end do

  end subroutine CommonType_printCSRdble

  !> Print CSRdble
  subroutine CommonType_printCSRArray2dble(CSR1)

    type(CSRArray2dble), intent(inout) :: CSR1
    integer :: i, j

    write(*,'(A4,I3)') "Nb:  ", CSR1%Nb
    do i=1,CSR1%Nb
      write(*,"(A4,I3,A8,I2)") "Row  ", i, " : size=", CSR1%Pt(i+1)-CSR1%Pt(i)

      do j=CSR1%Pt(i)+1,CSR1%Pt(i+1)
        write(*,"(I3)",advance="no") CSR1%Num(j)
      end do
      print*,""
    end do

    print*, ""
    print*, ""

    do i=1,CSR1%Nb
      write(*,"(A4,I3,A8,I2)") "Row  ", i, " : size=", CSR1%Pt(i+1)-CSR1%Pt(i)

      do j=CSR1%Pt(i)+1,CSR1%Pt(i+1)
        ! write(*,"(D2.6)",advance="no"), CSR1%Val(j)
        print*, CSR1%Val(j, :, :)
      end do
      print*,""
    end do

  end subroutine CommonType_printCSRArray2dble


  !> \brief Define operator = between two TYPE_IdNode:  x2 = x1
  subroutine assign_IdNode_equal(x2, x1)

    TYPE(TYPE_IdNode), intent(in) :: x1
    TYPE(TYPE_IdNode), intent(out) :: x2

    x2%Proc = x1%Proc
    x2%Frac = x1%Frac
    x2%P = x1%P
    x2%T = x1%T

  end subroutine assign_IdNode_equal


  function CSR_first_of_all(obj) result(j)
    type(CSR), intent(in) :: obj
    integer(c_int) :: j

#ifndef NDEBUG
    if(obj%Pt(1)/=0) &
        call CommonMPI_abort('inconsistent CSR')
#endif

    j = obj%Pt(1) + 1
  
  end function CSR_first_of_all

  function CSR_last_of_all(obj) result(j)
    type(CSR), intent(in) :: obj
    integer(c_int) :: j

#ifndef NDEBUG
    if(obj%Pt(1)/=0) &
        call CommonMPI_abort('inconsistent CSR')
#endif

    j = obj%Pt(obj%Nb+1)
  
  end function CSR_last_of_all

  subroutine CSR_associate_row(obj, i, p)
    type(CSR), target, intent(in) :: obj
    integer(c_int), intent(in) :: i
    integer(c_int), pointer, dimension(:), intent(out) :: p
  
#ifndef NDEBUG
    if(i<1.or.i>obj%Nb) &
        call CommonMPI_abort('wrong CSR index')
    if(size(obj%Pt)/=obj%Nb+1) &
        call CommonMPI_abort('inconsistent CSR')
#endif
    
    p=>obj%Num(obj%Pt(i)+1:obj%Pt(i+1))
  
  end subroutine CSR_associate_row

  function CSR_row_view(obj, i) result(p)
    type(CSR), target, intent(in) :: obj
    integer(c_int), intent(in) :: i
    integer(c_int), pointer, dimension(:) :: p

    call CSR_associate_row(obj, i, p)
  
  end function CSR_row_view
  

end module CommonType
