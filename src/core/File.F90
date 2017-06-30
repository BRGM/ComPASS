MODULE File

  IMPLICIT NONE

CONTAINS

  SUBROUTINE File_open(file_name)
    
    CHARACTER(LEN=*), INTENT(IN) :: file_name

    LOGICAL :: is_open
    INTEGER :: ios

    INQUIRE(UNIT=333, OPENED=is_open)
    IF(is_open)THEN
       PRINT *, "Error in File_open, a file is already open:"
       PRINT *, TRIM(file_name)
       STOP
    ENDIF

    OPEN(333, FILE=file_name, IOSTAT=ios)

    IF (ios /= 0) THEN 
       PRINT *, "Error in File_open, file not found or busy:"
       PRINT *, TRIM(file_name)
       STOP
    ENDIF
  END SUBROUTINE File_open


  SUBROUTINE File_close()

    LOGICAL :: is_open

    INQUIRE(UNIT=333, OPENED=is_open)
    IF(.NOT. is_open)THEN
       PRINT *,"Error in File_close, no opened file"
       STOP
    ENDIF

    CLOSE(333)
  END SUBROUTINE File_close


  SUBROUTINE File_is_EOF(is_EOF)

    LOGICAL, INTENT(OUT) :: is_EOF

    LOGICAL :: is_open
    INTEGER :: ios

    INQUIRE(UNIT=333, OPENED=is_open, IOSTAT=ios)

    IF(.NOT. is_open)THEN
       PRINT *,"Error in File_is_EOF, no opened file"
       STOP
    ENDIF

    is_EOF = ios < 0
  END SUBROUTINE File_is_EOF


  SUBROUTINE File_read_int_array(n, m, A)

    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: m
    INTEGER, INTENT(OUT) :: A(n,m)

    LOGICAL :: is_open
    INTEGER :: ios

    INQUIRE(UNIT=333, OPENED=is_open)
    IF(.NOT. is_open)THEN
       PRINT *,"Error in File_read_int_array, no opened file"
       PRINT *, n, m
       STOP
    ENDIF

    READ(333,*,IOSTAT=ios) A

    IF (ios /= 0) THEN 
       PRINT *,"Error in File_read_int_array, wrong matrix size or invalid character:"
       PRINT *, n, m
       STOP
    ENDIF
  END SUBROUTINE File_read_int_array


  SUBROUTINE File_read_double_array(n, m, A)

    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: m
    DOUBLE PRECISION, INTENT(OUT) :: A(n,m)

    LOGICAL :: is_open
    INTEGER :: ios

    INQUIRE(UNIT=333, OPENED=is_open)
    IF(.NOT. is_open)THEN
       PRINT *,"Error in File_read_double_array, no opened file"
       STOP
    ENDIF

    READ(333,*,IOSTAT=ios) A

    IF (ios /= 0) THEN 
       PRINT *,"Error in File_read_double_array, wrong matrix size or invalid character:"
       PRINT *, n, m
       STOP
    ENDIF
  END SUBROUTINE File_read_double_array


  SUBROUTINE File_read_line(line)

    CHARACTER(LEN=256), INTENT(OUT) :: line

    LOGICAL :: is_open
    INTEGER :: ios

    INQUIRE(UNIT=333, OPENED=is_open)
    IF(.NOT. is_open)THEN
       PRINT *,"Error in File_read_line, no opened file"
       STOP
    ENDIF

    READ(333,'(A)',IOSTAT=ios) line
  END SUBROUTINE File_read_line

END MODULE File


