    program ComPASS

    use NN

    logical :: EnterMainLoop = .True.
    integer :: TimeIter = 0
    integer :: argc
    !FIXME: 250 (is this always enough for long pathnames ?)
    character (250) :: &
        meshfile, logfile, outputdir, &
        TimeIterStr

    ! Input mesh name
    argc = IARGC()

    call GetArg(1, meshfile)
    call GetArg(2, logfile)
    call GetArg(3, outputdir)

    call NN_init(trim(meshfile), trim(logfile), trim(outputdir))

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

    end program ComPASS
