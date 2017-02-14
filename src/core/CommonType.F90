module CommonType

  implicit none

  !> Array 1d integer
  type ARRAY1Int
    integer, allocatable, dimension(:) :: Val
  end type ARRAY1Int

  !> Array 1d double precision
  type ARRAY1dble
    double precision, allocatable, dimension(:) :: Val
  end type ARRAY1dble
  
  !> Array 2d double precision
  type ARRAY2dble
    double precision, allocatable, dimension(:,:) :: Array2d
  end type ARRAY2dble

  !> Array 3d double precision
  type ARRAY3dble
    double precision, allocatable, dimension(:,:,:) :: Array3d
  end type ARRAY3dble

  !> Standar type CSR with Pt, Num, and Val (1d integer)
  type CSR
    integer :: Nb
    integer, allocatable, dimension(:) :: Pt
    integer, allocatable, dimension(:) :: Num
    integer, allocatable, dimension(:) :: Val
  end type CSR

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
     character :: Proc  !< "o"/"g": own/ghost
     character :: Frac  !< "y"/"n": node in fracture/not in fracture
     character :: P     !< "d"/"n"/"i": dirichlet/newmann/interior for the Pressure
     character :: T     !< "d"/"n"/"i": dirichlet/newmann/interior for the Temperature
  end type Type_IdNode

  !> Array 1d type_IdNode
  type Array1IdNode
    type(Type_IdNode), allocatable, dimension(:) :: val
  end type Array1IdNode

  ! tmp constant values
  integer, parameter, private :: &
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

  
  !> \brief Define operator = between two TYPE_IdNode:  x2 = x1
  subroutine assign_IdNode_equal(x2, x1)

    TYPE(TYPE_IdNode), intent(in) :: x1
    TYPE(TYPE_IdNode), intent(out) :: x2

    x2%Proc = x1%Proc
    x2%Frac = x1%Frac
    x2%P = x1%P
    x2%T = x1%T

  end subroutine assign_IdNode_equal

end module CommonType
