!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

    program ComPASS

    use CommonMPI
    use NN

    implicit none

    logical :: EnterMainLoop = .True.
    integer :: TimeIter = 0
    integer :: argc
    !FIXME: 250 (is this always enough for long pathnames ?)
    character (250), target :: &
        meshfile, logfile, outputdir, &
        TimeIterStr

    ! Input mesh name
    argc = IARGC()

    call GetArg(1, meshfile)
    meshfile = trim(meshfile)

    call GetArg(2, logfile)
    logfile = trim(logfile)

    call GetArg(3, outputdir)
    outputdir = trim(outputdir)

    ! This initialize MPI first consequently PetscInitialize will not call it
    ! and will adapt a symmetric behavior letting us call MPI_finaliaze at the end
    ! of our program
    call MPI_Init(Ierr)

    call NN_init(meshfile, logfile, outputdir)

    if(argc==4) then
        call GetArg(4, TimeIterStr)
        read (TimeIterStr,'(I10)') TimeIter
        if(TimeIter<0) then
            EnterMainLoop = .False.
            TimeIter = 0
        end if
    else
        TimeIter = 0
    end if

    if(EnterMainLoop) then
        call NN_main(TimeIter, OutputDir)
    end if

    call NN_finalize()

    ! MPI_Init symmetric call (cf. supra)
    call MPI_Finalize(Ierr)

    end program ComPASS
