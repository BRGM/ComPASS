

MODULE Mesh

  USE File
  USE CompassMesh
  USE GmshMesh

  IMPLICIT NONE

  INTEGER, PARAMETER :: GMSH_FILE_FORMAT = 1
  INTEGER, PARAMETER :: COMPASS_FILE_FORMAT = 2


CONTAINS


  SUBROUTINE Mesh_select_mesh_format(mesh_file_name, mesh_format)

    CHARACTER(LEN=*), INTENT(IN) :: mesh_file_name
    INTEGER, INTENT(OUT) :: mesh_format

    CHARACTER(LEN=256) :: line

    CALL File_open(mesh_file_name)
    CALL File_read_line(line)
    CALL File_close()

    IF(TRIM(line) == "$MeshFormat")THEN
      mesh_format = GMSH_FILE_FORMAT
    ELSE IF(TRIM(line) == "Cartesian mesh")THEN
      mesh_format = COMPASS_FILE_FORMAT
    ELSE 
      PRINT *,"Error in Mesh_select_mesh_format, format not supported:"
      PRINT *, TRIM(mesh_file_name)
      STOP
    ENDIF
  END SUBROUTINE Mesh_select_mesh_format


  SUBROUTINE Mesh_read_mesh_file(mesh_file_name, mesh_format)

    CHARACTER(LEN=*), INTENT(IN) :: mesh_file_name
    INTEGER, INTENT(IN) :: mesh_format

    IF(mesh_format == GMSH_FILE_FORMAT)THEN ! TODO
      PRINT *,"Error in Mesh_read_mesh_file, format not supported:"
      PRINT *, TRIM(mesh_file_name), mesh_format
      STOP
    ELSE IF(mesh_format == COMPASS_FILE_FORMAT)THEN
      CALL CompassMesh_read_mesh(mesh_file_name)
    ELSE 
      PRINT *,"Error in Mesh_read_mesh_file, format not supported:"
      PRINT *, TRIM(mesh_file_name), mesh_format
      STOP
    ENDIF
  END SUBROUTINE Mesh_read_mesh_file


END MODULE Mesh
