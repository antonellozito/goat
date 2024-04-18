!======================================================================!
!                                                                      !
!                            DOCUMENTATION                             !
!                                                                      !
!======================================================================!
! Description
!============
! Small program to translate the GOAToptions.dat file from Matlab to 
! Fortran format or from Fortran to Matlab. It uses the input file 
! parser implementation of GOAT to do the reading of the key-value 
! pairs, and then, if the pair should be replaced, replaces it with the
! Matlab/Fortran variant. Values are not adjusted. A list is maintained
! that saves all the keywords that should be replaced. The program 
! expects the following command line input (in this order):
! 
!   - filepath to the GOAToptions.dat file
!   - filepath to where the output should be written
!   - string that is either 'matlab' or 'fortran', which corresponds to 
!   - the format of the input file. So if set to 'matlab', the output 
!       will be in 'fortran' format and vice versa. 

! Notes
!======
! Note 1: due to the use of the key value pair extraction routine, only active
! (not outcommented) lines will be read and replaced! Furthermore, any
! comments present in changed lines will be removed (only the key value
! pair is written ...)
!
! Note 2: some options are not simply translated by replacing the 
! keywords. For example, if formerly a file was specified using separate
! filename and directory inputs and is now done using a single file path,
! this requires special attention. These keys should be manually adjusted!
!
! Note 3: some options may or may not be available in one or the other 
! module. No operations will be done on these options as they will be 
! ignored anyway. Note that this may lead to different behavior between 
! Matlab and fortran version (anyway it is never guaranteed that 
! exactly the same behavior will be observed)! 
!
! Note 4: the mapping between options may not be bijective. 

program TranslateGOAToptions

    ! Modules
    !========
    use mod_inputfileparser

    ! The usual
    implicit none

    ! Declare variables
    !==================
    ! Arguments
    type(StringUDT), allocatable       :: args(:)
    character(:), allocatable   :: inputfilepath, outputfilepath, &
        inputformat

    ! Auxiliary
    type(StringUDT), allocatable        ::  matlabkeys(:), &
        fortrankeys(:), inputkeys(:), outputkeys(:)
    character(:), allocatable   :: outputformat, &
        mgo, fgo, tempkey, tempvalue, thisline, &
        writeline, dirkey, filekey, &
        dirval, fileval
    character(len=4192), allocatable :: temp(:)
    character(len=4192)              :: temparg

    logical                     :: reachedeof, haspair, isfound

    integer                     :: nargin, inputfid, outputfid, &
        openstatus
    
    ! Loop
    integer(I8)                 :: i, k

    
    ! Initialize
    !===========
    ! File IDs
    inputfid = 10
    outputfid = 20

    ! Print header
    print *, '!=======================================================!'
    print *, '!             GOAT INPUT FILE TRANSLATOR                !'
    print *, '!=======================================================!'

    ! Get number of inputs
    nargin = command_argument_count()

    ! Check
    if (nargin /= 3) then 
        print *, 'Error: three input arguments expected ' // &
            '(input filepath, output filepath, input format)'
        print *, 'input format should be "matlab" or "fortran"'
        stop
    end if 

    ! Get inputs
    allocate(args(3))
    do i = 1, nargin 
        call get_command_argument(i, temparg)
        args(i)%val = trim(temparg)
    end do

    ! Unpack
    inputfilepath   = trim(args(1)%val)
    outputfilepath  = trim(args(2)%val)
    inputformat     = trim(args(3)%val)

    
    
    ! Check
    if ((inputformat == 'matlab') .or. (inputformat == 'Matlab') .or. &
        (inputformat == 'MATLAB')) then 
            inputformat = 'matlab'
            outputformat = 'fortran'
    elseif ((inputformat == 'fortran') .or. (inputformat == 'Fortran') .or. &
        (inputformat == 'FORTRAN')) then 
            inputformat = 'fortran'
            outputformat = 'matlab'
    else
        print *, 'Error: unknown input format: ' // inputformat 
        print *, 'input format should be "matlab" or "fortran"'
        stop
    end if
    if (inputfilepath == outputfilepath) then 
        print *, 'Input and output filepath should be different, exiting...'
        stop 
    end if 
    
    ! Print inputs
    print *, '!-------'
    print *, '! Inputs'
    print *, '!-------'
    print *, 'Reading from file: ' // inputfilepath
    print *, 'Writing to file: ' // outputfilepath 
    print *, 'Input format: ' // inputformat 
    print *, 'Output format: ' // outputformat 
    
    ! Set key mapping
    !================
    ! Set some useful prefixes etc
    mgo = 'GOAToptions.' ! goatoptions prefix
    fgo = 'goat.'
    

    ! Set mapping
    temp = [character(len=4192) :: &
        mgo // 'mfloaddir', &
        mgo // 'mfloadfile', &
        'gd.grid.gridloaddir', &
        'gd.grid.gridloadfile', &
        'gd.grid.vesselloaddir', &
        'gd.grid.vesselloadfile', &
        'gd.mf.loaddir', &
        'gd.mf.loadfile', &
        mgo // 'readtype', &
        mgo // 'readfile',&
        mgo // 'debug', &
        mgo // 'meth', &
        mgo // 'writedir', &
        mgo // 'strfile', &
        mgo // 'GDinputfilepath', &
        mgo // 'write_final', &
        mgo // 'write_traduitb2us', &
        mgo // 'write_b2agdat', &
        mgo // 'write_Xpointdata', &
        mgo // 'write_OMPdata', &
        mgo // 'vesselmode', &
        mgo // 'slab', &
        mgo // 'artificial_slab', &
        mgo // 'OMP_r', &
        mgo // 'OMP_z', &
        mgo // 'IMP_r', &
        mgo // 'IMP_z', &
        mgo // 'GDtoGA.facelabelsubfrom', &
        mgo // 'GDtoGA.facelabelsubto', &
        mgo // 'GDtoGA.facelabelmappingGG', &
        mgo // 'GDtoGA.facelabelmappingGD', &
        mgo // 'GGtoGA.facelabelsubfrom', &
        mgo // 'GGtoGA.facelabelsubto', &
        mgo // 'GGtoGA.facelabelmappingGG', &
        mgo // 'GGtoGA.facelabelmappingGA', &
        mgo // 'data.RBtor', &
        'gd.design.cfv.par.PRPB.lengthparam', &
        'gd.design.cfv.par.PRPB.biassep', &
        'gd.design.cfv.par.PRPB.biasinf', &
        'gd.design.cfv.par.PRPB.weightlengthparam', &
        'gd.design.cfv.par.PRPB.weightsep', &
        'gd.design.ec.par.fixedfluxvalues.coreflux', &
        'gd.design.ec.par.fixedfluxvalues.outerflux', &
        'gd.grid.vesselshapemeth', &
        'gd.grid.vesselresx', &
        'gd.grid.vesselresy', &
        'gd.grid.vesseloffsetfracx', &
        'gd.grid.vesseloffsetfracy', &
        'gd.grid.vesselC', &
        'gd.grid.vesselM', &
        'gd.grid.TP', &
        'gd.grid.TPind', &  
        'gd.grid.exclude', &
        'gd.grid.vesselloadmeth', &
        'gd.grid.refinevessel' , &
        'gd.grid.maxvesseldist', &
        'gd.mf.readmeth', &
        'gd.mf.interpC', &
        'gd.mf.interpM' , &
        'gd.mf.interpmeth', &
        'gd.num.opt.itmax' , &
        'gd.num.opt.tol' , &
        'gd.num.KKTdriver.docfrelax' , &
        'gd.num.KKTdriver.rxf', &
        'gd.num.KKTdriver.rxfdesign', &
        'gd.num.KKTdriver.minrxf', &
        'gd.num.KKTdriver.rxfdec', &
        'gd.num.KKTsolveroptions.ls.dolinesearch', &
        'gd.num.KKTsolveroptions.ls.itmax', &
        'gd.num.KKTdriver.doeqconrelax',& 
        'gd.num.KKTdriver.rxfeqcon',& 
        'gd.num.KKTdriver.minrxfeqcon', &
        'gd.num.KKTdriver.rxfeqcondec', &
        'gd.num.KKTdriver.doineqconrelax', &
        'gd.num.KKTdriver.rxfineqcon'       , &
        'gd.num.KKTdriver.minrxfineqcon'    , &
        'gd.num.KKTdriver.rxfineqcondec'    , &
        'gd.num.KKTdriver.dovisualization' , &
        'gd.num.KKTdriver.makevideo'        , &
        'gd.num.KKTdriver.stabilize'        , &
        'gd.num.KKTdriver.doconhess'        , &
        'gd.num.KKTsolver'                  &
    ]
        
    allocate(matlabkeys(size(temp, 1)))
    do i = 1, size(temp, 1)
        matlabkeys(i)%val = trim(temp(i))
    end do

    temp = [character(len=4192) :: &
        ' ', &
        ' ', &
        ' ', &
        ' ', &
        ' ', &
        ' ', &
        ' ', &
        ' ', &
        fgo // 'grid.readmeth', &
        fgo // 'gridfilepath', &
        fgo // 'debug', &
        fgo // 'meth', &
        fgo // 'writefilepath', &
        fgo // 'structurefilepath', &
        fgo // 'GDinputfilepath', &
        fgo // 'write_final', &
        fgo // 'write_traduitb2us', &
        fgo // 'write_b2agdat', &
        fgo // 'write_Xpointdata', &
        fgo // 'write_OMPdata', &
        fgo // 'vesselmode', &
        fgo // 'slab', &
        fgo // 'artificial_slab', &
        fgo // 'OMP_r', &
        fgo // 'OMP_z', &
        fgo // 'IMP_r', &
        fgo // 'IMP_z', &
        fgo // 'GGtoGD.facelabelsubfrom', &
        fgo // 'GGtoGD.facelabelsubto', &
        fgo // 'GGtoGD.facelabelmappingGG', &
        fgo // 'GGtoGD.facelabelmappingGD', &
        fgo // 'GGtoGA.facelabelsubfrom', &
        fgo // 'GGtoGA.facelabelsubto', &
        fgo // 'GGtoGA.facelabelmappingGG', &
        fgo // 'GGtoGA.facelabelmappingGA', &
        fgo // 'data.RBtor', &
        'gd.design.cfv.par.PRPB.biasdecaylength', &
        'gd.design.cfv.par.PRPB.biasatsep', &
        'gd.design.cfv.par.PRPB.biasatinf', &
        'gd.design.cfv.par.PRPB.weightdecaylength', &
        'gd.design.cfv.par.PRPB.weightatsep', &
        'gd.design.ec.par.fixedfluxvalues.corefluxval', &
        'gd.design.ec.par.fixedfluxvalues.outerfluxval', &
        'goat.vessel.shapemeth', &
        'goat.vessel.resx', &
        'goat.vessel.resy', &
        'goat.vessel.offsetfracx', &
        'goat.vessel.offsetfracy', &
        'goat.vessel.interpC', &
        'goat.vessel.interpM', &
        'goat.vessel.TP' , &
        'goat.vessel.TPind' , &
        'goat.vessel.exclude', &
        'goat.vessel.readmeth', &
        'goat.vessel.refinevessel' , &
        'goat.vessel.maxvesseldist', &   
        'goat.mf.readmeth', &
        'goat.mf.interpC', &
        'goat.mf.interpM', &
        'goat.mf.interpmeth', &
        'gd.opt.num.itmax', &
        'gd.opt.num.tol', &
        'gd.opt.num.docfrelax', &
        'gd.opt.num.rxf', &
        'gd.opt.num.rxfdesign', &
        'gd.opt.num.rxfmin'                 , &
        'gd.opt.num.rxfdec'                 , &
        'gd.opt.num.ls.dolinesearch'        , &
        'gd.opt.num.ls.itmax'               , &
        'gd.num.KKTdriver.doeqconrelax'     , &
        'gd.num.KKTdriver.rxfeqcon'         , &
        'gd.num.KKTdriver.minrxfeqcon'      , &
        'gd.num.KKTdriver.rxfeqcondec'      , &
        'gd.num.KKTdriver.doineqconrelax'   , &
        'gd.num.KKTdriver.rxfineqcon'       , &
        'gd.num.KKTdriver.minrxfineqcon'    , &
        'gd.num.KKTdriver.rxfineqcondec'    , &
        'gd.num.KKTdriver.dovisualization'  , &
        'gd.num.KKTdriver.makevideo'        , &
        'gd.num.KKTdriver.stabilize'        , &
        'gd.num.KKTdriver.doconhess'        , &
        'gd.num.KKTsolver'                  &
        ]

    allocate(fortrankeys(size(temp, 1)))
    do i = 1, size(temp, 1)
        fortrankeys(i)%val = trim(temp(i))
    end do

    ! Check which should be input 
    if (inputformat == 'matlab') then 
        inputkeys = matlabkeys 
        outputkeys = fortrankeys 
    else
        inputkeys = fortrankeys 
        outputkeys = matlabkeys
    end if 

    ! Print mapping
    !==============
    print *, '!--------'
    print *, '! Mapping'
    print *, '!--------'
    print *, 'matlab keys   <=>     fortran keys'
    do i = 1, size(matlabkeys, 1)
        print *, matlabkeys(i)%val, '    <=>    ', fortrankeys(i)%val
    end do

    ! Read and map
    !=============
    ! Open files
    ! Open the input file, check if it exists
    open(unit=inputfid, file=inputfilepath, status='old', &
        iostat=openstatus)

    if (openstatus > 0) then 
        ! Something wrong when reading file - stop
        print *, 'Could not open input file, exiting...'
        stop
    elseif (openstatus < 0) then 
        ! File appears to be empty
        print *, 'File appears to be empty, exiting...'
    end if

    ! Open the output file
    open(unit=outputfid, file=outputfilepath, status='unknown', &
        iostat=openstatus)

    if (openstatus /= 0) then 
        ! Something wrong when reading file - stop
        print *, 'Could not open output file, exiting...'
        stop
    end if

    do while ( (.not. reachedeof) )

        ! Read single line
        call ReadSingleLine(inputfid, thisline, reachedeof)

        ! Set initial write line
        writeline = thisline 

        ! Check if it contains a key-value pair
        call GetKeyValuePair(thisline, tempkey, tempvalue, haspair)
        
        if (haspair) then 
            ! Check if the key matches any of the keys to be changed, if 
            ! so: replace and exit
            k = 1
            isfound = .false.
            do while ((k <= size(inputkeys, 1)) .and. (.not. isfound))
                ! Check
                if (inputkeys(k)%val == tempkey) then 
                    ! Found
                    isfound = .true. 

                    ! Write line
                    writeline = delimiter // outputkeys(k)%val // delimiter // &
                        '    ' // delimiter // tempvalue // delimiter 
                end if 

                ! Update counter
                k = k + 1

            end do
        end if 

        ! Write (not if empty)
        if (haspair) then 
            if (outputkeys(k-1)%val(1:1) /= ' ') then 
                write(outputfid, '(a)') writeline
            end if 
        else
            write(outputfid, '(a)') writeline
        end if 

    end do

    ! Check special cases for matlab input
    !=====================================
    if (inputformat == 'matlab') then 

        ! Filepath magnetic field
        dirkey = 'GOAToptions.mfloaddir'
        filekey = 'GOAToptions.mfloadfile'
        call GetValueWithKey(inputfid, dirkey, dirval, isfound)
        if (isfound) then 
            call GetValueWithKey(inputfid, filekey, fileval, isfound)
            if (.not. isfound) then 
                ! Weird, throw warning and don't continue
                print *, 'Magnetic field loading path translation: only found dir, not file. Not replacing...'
            else
                ! Replace
                writeline = delimiter // 'goat.magneticfieldfilepath' // &
                    delimiter // '    ' // delimiter // dirval // filesepchar // &
                    fileval // delimiter
                write(outputfid, '(a)') writeline
            end if
        end if

        ! Filepath grid from grid options
        dirkey = 'gd.grid.gridloaddir'
        filekey = 'gd.grid.gridloadfile'
        call GetValueWithKey(inputfid, dirkey, dirval, isfound)
        if (isfound) then 
            call GetValueWithKey(inputfid, filekey, fileval, isfound)
            if (.not. isfound) then 
                ! Weird, throw warning and don't continue
                print *, 'Grid loading path translation: only found dir, not file. Not replacing...'
            else
                ! Replace
                writeline = delimiter // 'goat.gridfilepath' // &
                    delimiter // '    ' // delimiter // dirval // filesepchar // &
                    fileval // delimiter
                write(outputfid, '(a)') writeline
            end if
        end if

        ! Filepath vessel from grid options
        dirkey = 'gd.grid.vesselloaddir'
        filekey = 'gd.grid.vesselloadfile'
        call GetValueWithKey(inputfid, dirkey, dirval, isfound)
        if (isfound) then 
            call GetValueWithKey(inputfid, filekey, fileval, isfound)
            if (.not. isfound) then 
                ! Weird, throw warning and don't continue
                print *, 'Vessel loading path translation: only found dir, not file. Not replacing...'
            else
                ! Replace
                writeline = delimiter // 'goat.structurefilepath' // &
                    delimiter // '    ' // delimiter // dirval // filesepchar // &
                    fileval // delimiter
                write(outputfid, '(a)') writeline
            end if
        end if

        ! Filepath magnetic field from grid
        dirkey = 'gd.mf.loaddir'
        filekey = 'gd.mf.loadfile'
        call GetValueWithKey(inputfid, dirkey, dirval, isfound)
        if (isfound) then 
            call GetValueWithKey(inputfid, filekey, fileval, isfound)
            if (.not. isfound) then 
                ! Weird, throw warning and don't continue
                print *, 'Magnetic field loading path translation: only found dir, not file. Not replacing...'
            else
                ! Replace
                writeline = delimiter // 'goat.magneticfieldfilepath' // &
                    delimiter // '    ' // delimiter // dirval // filesepchar // &
                    fileval // delimiter
                write(outputfid, '(a)') writeline
            end if
        end if

    end if

    ! Check special cases for fortran input
    if (inputformat == 'fortran') then 

        ! Filepath magnetic field
        dirkey = 'goat.magneticfieldfilepath'
        call GetValueWithKey(inputfid, dirkey, dirval, isfound)
        if (isfound) then 
            ! Write dir by default as '.', file should be the rest
            writeline = delimiter // 'GOAToptions.mfloaddir' // &
                delimiter // '    ' // delimiter // '.' // delimiter
            write(outputfid, '(a)') writeline

            ! Check format
            if (dirval(1:1) == filesepchar) then 
                ! Write without the file separator
                writeline = delimiter // 'GOAToptions.mfloadfile' // &
                    delimiter // '    ' // delimiter // dirval(2:len(dirval)) // delimiter
                write(outputfid, '(a)') writeline
            elseif (dirval(1:1) == '.') then 
                ! Write without the file separator
                writeline = delimiter // 'GOAToptions.mfloadfile' // &
                    delimiter // '    ' // delimiter // dirval(3:len(dirval)) // delimiter
                write(outputfid, '(a)') writeline
            else
                ! Write 
                writeline = delimiter // 'GOAToptions.mfloadfile' // &
                    delimiter // '    ' // delimiter // dirval // delimiter
                write(outputfid, '(a)') writeline
            end if 
        end if

        ! Filepath grid from grid options
        dirkey = 'goat.gridfilepath'
        call GetValueWithKey(inputfid, dirkey, dirval, isfound)
        if (isfound) then 
            ! Write dir by default as '.', file should be the rest
            writeline = delimiter // 'gd.grid.gridloaddir' // &
                delimiter // '    ' // delimiter // '.' // delimiter
            write(outputfid, '(a)') writeline

            ! Check format
            if (dirval(1:1) == filesepchar) then 
                ! Write without the file separator
                writeline = delimiter // 'gd.grid.gridloadfile' // &
                    delimiter // '    ' // delimiter // dirval(2:len(dirval)) // delimiter
                write(outputfid, '(a)') writeline
            elseif (dirval(1:1) == '.') then 
                ! Write without the file separator
                writeline = delimiter // 'gd.grid.gridloadfile' // &
                    delimiter // '    ' // delimiter // dirval(3:len(dirval)) // delimiter
                write(outputfid, '(a)') writeline
            else
                ! Write 
                writeline = delimiter // 'gd.grid.gridloadfile' // &
                    delimiter // '    ' // delimiter // dirval // delimiter
                write(outputfid, '(a)') writeline
            end if 
        end if

        ! Filepath vessel from grid options
        dirkey = 'goat.structurefilepath'
        call GetValueWithKey(inputfid, dirkey, dirval, isfound)
        if (isfound) then 
            ! Write dir by default as '.', file should be the rest
            writeline = delimiter // 'gd.grid.vesselloaddir' // &
                delimiter // '    ' // delimiter // '.' // delimiter
            write(outputfid, '(a)') writeline

            ! Check format
            if (dirval(1:1) == filesepchar) then 
                ! Write without the file separator
                writeline = delimiter // 'gd.grid.vesselloadfile' // &
                    delimiter // '    ' // delimiter // dirval(2:len(dirval)) // delimiter
                write(outputfid, '(a)') writeline
            elseif (dirval(1:1) == '.') then 
                ! Write without the file separator
                writeline = delimiter // 'gd.grid.vesselloadfile' // &
                    delimiter // '    ' // delimiter // dirval(3:len(dirval)) // delimiter
                write(outputfid, '(a)') writeline
            else
                ! Write 
                writeline = delimiter // 'gd.grid.vesselloadfile' // &
                    delimiter // '    ' // delimiter // dirval // delimiter
                write(outputfid, '(a)') writeline
            end if 
        end if

        ! Filepath magnetic field from grid
        dirkey = 'goat.magneticfieldfilepath'
        call GetValueWithKey(inputfid, dirkey, dirval, isfound)
        if (isfound) then 
            ! Write dir by default as '.', file should be the rest
            writeline = delimiter // 'gd.mf.loaddir' // &
                delimiter // '    ' // delimiter // '.' // delimiter
            write(outputfid, '(a)') writeline

            ! Check format
            if (dirval(1:1) == filesepchar) then 
                ! Write without the file separator
                writeline = delimiter // 'gd.mf.loadfile' // &
                    delimiter // '    ' // delimiter // dirval(2:len(dirval)) // delimiter
                write(outputfid, '(a)') writeline
            elseif (dirval(1:1) == '.') then 
                ! Write without the file separator
                writeline = delimiter // 'gd.mf.loadfile' // &
                    delimiter // '    ' // delimiter // dirval(3:len(dirval)) // delimiter
                write(outputfid, '(a)') writeline
            else
                ! Write 
                writeline = delimiter // 'gd.mf.loadfile' // &
                    delimiter // '    ' // delimiter // dirval // delimiter
                write(outputfid, '(a)') writeline
            end if 
        end if

    end if 

    ! Close files
    close(inputfid)
    close(outputfid)


end program