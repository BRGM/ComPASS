

MODULE CompassMesh

  USE GlobalMesh

  IMPLICIT NONE

CONTAINS

  SUBROUTINE CompassMesh_read_mesh(mesh_file_name)

    CHARACTER(LEN=*), INTENT(IN) :: mesh_file_name

    CALL GlobalMesh_Make_read_file(mesh_file_name)
  END SUBROUTINE CompassMesh_read_mesh

end module CompassMesh
