module VisuVTK

  use PathUtilities

  use CommonType
  use CommonMPI
  use MeshSchema

  use iso_c_binding

  implicit none

  ! ptr to class VisuVTK_Time
  type(c_ptr), private :: visuptr

  ! output dir 
  character(len=200), private :: OutputDir

  ! Two types of Mesh
  !  1. Cartesien/Tetra/Hexahedron: in vtk, not need to define faces
  !  2. General mesh: cell is defined by FacebyCell, face is defined by NodebyFace

  ! macro to distinguish mesh types for visu
  integer, parameter :: &
      MESH_GEN = 0,    & ! general mesh:          vtk_polyhedron
      MESH_CAR = 1,    & ! cartesispen mesh:      vtk_voxel,       vtk_quad
      MESH_TET = 2,    & ! tet mesh:              vtk_tetra,       vtk_triangle
      MESH_HEX = 3,    & ! hex mesh (hexahedron): vtk_hexahedron,  vtk_quad
      MESH_WEDGE = 4     ! wedge mesh           : vtk_wedge,       vtk_quad, vtk_triangle


  ! the times which are saved for visu
  double precision, dimension(:), allocatable, private :: VisuTimes  

  ! Nb of times which are saved for visu
  integer, private :: NbVisuTimes = 0

  interface

    ! 1. interface for visu time
    function visuvtk_time_initcxx(meshtypecpp, outputdircpp, &
        commRankcpp, commSizecpp, &
        NbCompcpp, NbPhasecpp, MCPcpp, IndThermiquecpp, &
        !
        NbCellOwncpp, NbFaceOwncpp, NbNodeLocalcpp, &
        NbWellInjOwncpp, NbWellProdOwncpp, &
        !
        NbFracOwncpp, FracToFaceLocalcpp, &
        !
        NbNodeOwncpp, &
        !
        NodebyCellLocal_Nbcpp, NodebyCellLocal_Ptcpp, NodebyCellLocal_Numcpp, &
        FacebycellLocal_Nbcpp, FacebycellLocal_Ptcpp, FacebycellLocal_Numcpp, &
        NodebyfaceLocal_Nbcpp, NodebyfaceLocal_Ptcpp, NodebyfaceLocal_Numcpp, &
                                !
        NbEdgebyWellInjcpp, NbEdgebyWellProdcpp, &
        NumNodebyEdgebyWellInjcpp, NumNodebyEdgebyWellProdcpp, &
                                !
        XNodeLocalcpp, XCellLocalcpp) &

        result(this) bind(C, name="visuvtk_time_initcxx_")

      use iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_null_char

      type(c_ptr) :: this

      ! mesh type
      integer(c_int), value :: meshtypecpp

      ! output dir
      character(c_char) :: outputdircpp(*)

      ! commRank and commSize
      integer(c_int), value :: commRankcpp
      integer(c_int), value :: commSizecpp

      ! NbPhase, NbComp, MCP, IndThermique
      integer(c_int), value :: NbCompcpp
      integer(c_int), value :: NbPhasecpp
      integer(c_int) :: MCPcpp(*)
      integer(c_int), value :: IndThermiquecpp

      ! Nb
      integer(c_int), value :: NbCellOwncpp
      integer(c_int), value :: NbFaceOwncpp
      integer(c_int), value :: NbNodeLocalcpp
      integer(c_int), value :: NbWellInjOwncpp  ! nb of inj well own
      integer(c_int), value :: NbWellProdOwncpp ! nb of prod well own 

      ! Nb of edges of each well
      integer(c_int) :: NbEdgebyWellInjcpp(*)
      integer(c_int) :: NbEdgebyWellProdcpp(*)

      ! Nodes connected to edges in wells
      integer(c_int) :: NumNodebyEdgebyWellInjcpp(*)
      integer(c_int) :: NumNodebyEdgebyWellProdcpp(*)

      ! NodebyCell
      integer (c_int), value :: NodebyCellLocal_Nbcpp
      integer (c_int) :: NodebyCellLocal_Ptcpp(*)       
      integer (c_int) :: NodebyCellLocal_Numcpp(*)

      ! FacebyCell
      integer (c_int), value :: FacebyCellLocal_Nbcpp
      integer (c_int) :: FacebyCellLocal_Ptcpp(*)       
      integer (c_int) :: FacebyCellLocal_Numcpp(*)

      ! NodebyFace
      integer (c_int), value :: NodebyFaceLocal_Nbcpp
      integer (c_int) :: NodebyFaceLocal_Ptcpp(*)       
      integer (c_int) :: NodebyFaceLocal_Numcpp(*)

      ! XNodeLocal
      real (c_double) :: XNodeLocalcpp(*)

      ! XCellLocal
      real (c_double) :: XCellLocalcpp(*)

      ! Nb of frac
      integer (c_int), value :: NbFracOwncpp

      ! Nb of node
      integer (c_int), value :: NbNodeOwncpp

      ! FractoFaceLocal
      integer (c_int) :: FracToFaceLocalcpp(*)

    end function visuvtk_time_initcxx


    SUBROUTINE visuvtk_time_writedatacxx( &
        this, &
        NbVisuTimescpp, &
        datacellcpp, &
        datafraccpp, &
        datanodecpp, &
        datawellinjcpp, &
        datawellprodcpp) &
        BIND(C, name="visuvtk_time_writedatacxx_")

      USE ISO_C_BINDING, ONLY: c_ptr, c_int, c_double

      TYPE(c_ptr), VALUE :: this

      INTEGER(c_int), VALUE :: NbVisuTimescpp

      REAL(c_double) :: datacellcpp(*)
      REAL(c_double) :: datafraccpp(*)
      REAL(c_double) :: datanodecpp(*)
      REAL(c_double) :: datawellinjcpp(*)
      REAL(c_double) :: datawellprodcpp(*)

    ENDSUBROUTINE visuvtk_time_writedatacxx


    subroutine visuvtk_time_freecxx(this) &
        bind(C, name="visuvtk_time_freecxx_")

      use iso_c_binding, only: c_ptr, c_int

      type(c_ptr), value :: this

    end subroutine visuvtk_time_freecxx


    subroutine visuvtk_pvdwritercxx(dirname, Nbvisu, Timesvisu) &
        bind(C, name="visuvtk_pvdwritercxx_")

      use iso_c_binding       

      character(c_char) :: dirname(*)
      integer (c_int), value :: Nbvisu
      real (c_double) :: Timesvisu(*)

    end subroutine visuvtk_pvdwritercxx


    ! ! interface for visu without time
    ! subroutine VisuVTK_Visucxx( &
    !      meshtype,  &                                                       ! meshtype
    !      commRank, commSize, &                                              ! mpi info
    !      NbCellOwn,   NbFaceOwn,   NbNodeOwn,   &                           ! nb of ...
    !      NbCellLocal, NbFaceLocal, NbNodeLocal, &                           
    !      NodebyCellLocal_Nb, NodebyCellLocal_Pt, NodebyCellLocal_Num, &     !
    !      FacebyCellLocal_Nb, FacebyCellLocal_Pt, FacebyCellLocal_Num, &     ! connectivites
    !      NodebyFaceLocal_Nb, NodebyFaceLocal_Pt, NodebyFaceLocal_Num, &     !
    !      XNodeLocal, XCellLocal, &                                          ! X of cell and frac
    !      NbFracOwn, FracToFaceLocal, &                                      ! info of frac
    !      datacell, datafrac)                                                ! data to visu, cell and frac

    !   use iso_c_binding

    !   ! meshtype
    !   integer (c_int), VALUE :: meshtype

    !   ! mpi
    !   integer (c_int), VALUE :: commRank
    !   integer (c_int), VALUE :: commSize

    !   ! Nb of ...
    !   integer (c_int), VALUE :: NbCellLocal
    !   integer (c_int), VALUE :: NbFaceLocal
    !   integer (c_int), VALUE :: NbNodeLocal
    !   integer (c_int), VALUE :: NbCellOwn
    !   integer (c_int), VALUE :: NbFaceOwn
    !   integer (c_int), VALUE :: NbNodeOwn

    !   integer (c_int), VALUE :: NbFracOwn

    !   ! NodebyCell
    !   integer (c_int), VALUE :: NodebyCellLocal_Nb
    !   integer (c_int) :: NodebyCellLocal_Pt(*)       
    !   integer (c_int) :: NodebyCellLocal_Num(*)

    !   ! FacebyCell
    !   integer (c_int), VALUE :: FacebyCellLocal_Nb
    !   integer (c_int) :: FacebyCellLocal_Pt(*)       
    !   integer (c_int) :: FacebyCellLocal_Num(*)

    !   ! NodebyFace
    !   integer (c_int), VALUE :: NodebyFaceLocal_Nb
    !   integer (c_int) :: NodebyFaceLocal_Pt(*)       
    !   integer (c_int) :: NodebyFaceLocal_Num(*)

    !   ! XNodeLocal
    !   real (c_double) :: XNodeLocal(*)

    !   ! XCellLocal
    !   real (c_double) :: XCellLocal(*)

    !   ! FractoFaceLocal
    !   integer (c_int) :: FracToFaceLocal(*)

    !   ! data to visu
    !   real (c_double) :: datacell(*)
    !   real (c_double) :: datafrac(*)

    ! end subroutine VisuVTK_Visucxx

  end interface

contains

  ! new class VisuVTK_Time
  ! structures used for every time step
  subroutine VisuVTK_VisuTime_Init( &
      meshtype, dirname, &
      Tf, output_frequency, &
      NbComp, NbPhase, MCP, IndThermique)

    integer, intent(in) :: meshtype

    character(*), intent(in) :: dirname

    double precision, intent(in) :: &
        Tf, &      ! final time
        output_frequency ! output frequency

    integer, intent(in) :: &
        NbComp, NbPhase, MCP(NbComp, NbPhase), &
        IndThermique

    OutputDir = dirname

    ! init visu structure
    visuptr = visuvtk_time_initcxx(meshtype, trim(OutputDir)//C_NULL_CHAR, &
        commRank, commSize, &
        NbComp, NbPhase, MCP, IndThermique, &
        !
        NbCellOwn_Ncpus(commRank+1), NbFaceOwn_Ncpus(commRank+1), NbNodeLocal_Ncpus(commRank+1), &
        NbWellInjOwn_Ncpus(commRank+1), NbWellProdOwn_Ncpus(commRank+1), &
        NbFracOwn_Ncpus(commRank+1), FracToFaceLocal, &
        NbNodeOwn_Ncpus(commRank+1), &
        !
        NodebyCellLocal%Nb, NodebyCellLocal%Pt, NodebyCellLocal%Num, &
        FacebyCellLocal%Nb, FacebycellLocal%Pt, FacebycellLocal%Num, &
        NodebyFaceLocal%Nb, NodebyFaceLocal%Pt, NodebyFaceLocal%Num, &
        !
        NbEdgebyWellInjLocal, NbEdgebyWellProdLocal, &
        NumNodebyEdgebyWellInjLocal, NumNodebyEdgebyWellProdLocal, &
        !
        XNodeLocal, XCellLocal)

    ! allocate VisuTimes 
    ! FIXME: The array is copied and expanded each time a visualaisation ouput is made
    allocate(VisuTimes(NbVisuTimes))

  end subroutine VisuVTK_VisuTime_Init


  ! write data for time step
  SUBROUTINE VisuVTK_VisuTime_writedata( &
      t, &
      datacell, &
      datafrac, &
      datanode, &
      datawellinj, &
      datawellprod)

    DOUBLE PRECISION, INTENT(IN) :: t
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: datacell
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: datafrac
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: datanode
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: datawellinj
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: datawellprod

    character(len=200) :: output_path
    double precision, dimension(size(VisuTimes)) :: TmpVisuTimes  

    ! FIXME: The allocation/deallocation is done at each output
    TmpVisuTimes = VisuTimes
    deallocate(VisuTimes)
    
    NbVisuTimes = NbVisuTimes + 1
    
    allocate(VisuTimes(NbVisuTimes))
    VisuTimes(1:NbVisuTimes-1) = TmpVisuTimes
 
    VisuTimes(NbVisuTimes) = t ! all timesteps for .pvd are stored here

    write(output_path, '(A,I0)')  trim(OutputDir) // "/time_", NbVisuTimes-1
    call make_directory(output_path)

    ! write data 
    call visuvtk_time_writedatacxx( &
      visuptr, &
      NbVisuTimes, &
      datacell, &
      datafrac, &
      datanode, &
      datawellinj, &
      datawellprod)

  end subroutine VisuVTK_VisuTime_writedata


  ! free
  subroutine VisuVTK_VisuTime_free

    call visuvtk_time_freecxx(visuptr)

    deallocate( VisuTimes)

  end subroutine VisuVTK_VisuTime_free


  ! write pvd file 
  subroutine VisuVTK_VisuTime_pvdwriter

    if(commRank==0) then
      call VisuVTK_pvdwritercxx(trim(OutputDir)//C_NULL_CHAR, &
          NbVisuTimes, VisuTimes)
    end if

  end subroutine VisuVTK_VisuTime_pvdwriter

  ! ! visu without time steps
  ! subroutine VisuVTK_Visu(meshtype, datacell, datafrac)

  !   integer, intent(in) :: meshtype ! type of mesh car/tet/gen 

  !   double precision, dimension(:), allocatable, intent(in) :: &
  !        datacell, &  ! data of cell to visu
  !        datafrac     ! data of frac to visu

  !   call visuvtk_visucxx( &
  !        meshtype, &
  !        commRank, commSize, &
  !        NbCellOwn_Ncpus(commRank+1),   NbFaceOwn_Ncpus(commRank+1),   NbNodeOwn_Ncpus(commRank+1), &
  !        NbCellLocal_Ncpus(commRank+1), NbFaceLocal_Ncpus(commRank+1), NbNodeLocal_Ncpus(commRank+1), &
  !        NodebyCellLocal%Nb,NodebyCellLocal%Pt,NodebyCellLocal%Num, &
  !        FacebyCellLocal%Nb,FacebyCellLocal%Pt,FacebyCellLocal%Num, &
  !        NodebyFaceLocal%Nb,NodebyFaceLocal%Pt,NodebyFaceLocal%Num, &
  !        XNodeLocal, XCellLocal,                                    &
  !        NbFracOwn_Ncpus(commRank+1), FracToFaceLocal,              &
  !        datacell, datafrac)

  ! end subroutine VisuVTK_Visu

end module VisuVTK
