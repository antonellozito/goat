subroutine ReadTraduitUS(grid, filepath)

    ! Description
    !============
    ! Read in the necessary grid data from a traduit.out.b2us file.


    ! Initialize
    !===========
    ! Declare modules
    use gdmod_types 
    use mod_readwrite

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Main variables
    type(GridUDT)               :: grid

    integer                     :: filespec
    integer(I8)                 :: idum(0:9), idum2(1)    
    integer(I8)                 :: nc, nf, nv, nfsFc, nftCv, nftFc
    character(*)               :: filepath

    logical                     :: reachedeof, readTopologicalData 

    character(:), allocatable   :: chardummy   ! dummy array
    integer(I8), allocatable    :: cdummy(:,:) ! dummy array
    integer(I8), allocatable    :: fdummy(:,:) ! dummy array
    integer(I8), allocatable    :: vdummy(:,:) ! dummy array
    integer(I8), allocatable    :: ftdummy(:) ! dummy array

    real(R8), allocatable       :: cdummyr(:,:) ! dummy array
    real(R8), allocatable       :: fdummyr(:,:) ! dummy array
    real(R8), allocatable       :: vdummyr(:,:) ! dummy array
    real(R8), allocatable       :: fsdummyr(:) ! dummy array
    real(R8), allocatable       :: facelistdummy(:) ! dummy array
    integer(I8), allocatable    :: ftCvdum(:), ftFcdum(:), fsFcdum(:)

    real(R8), allocatable       :: vdata(:, :), cdata(:, :), &
        fsdata(:, :)
        
    integer(I8), allocatable    :: vlist(:), clist(:), &
        flist(:), ftlist(:), fslist(:), cdatai1(:, :), &
        fdatai(:, :), ftdatai(:, :), fsdatai(:, :), cdatai2(:, :)

    ! Loop
    integer(I8)                 :: i

    ! Data
    data filespec /60/
    
    ! Read grid dimensions & allocate
    !================================
    ! Open the file
    print *, 'reading grid in traduit format from file: ' // filepath
    open(unit = filespec, file = filepath)
    rewind(filespec)

    ! First, read the header with the version
    call ReadSingleLine(filespec, chardummy, reachedeof)
    if (reachedeof) then 
        call gdErrorHandler('ReadTraduitUS: reached EOF prematurely')
    end if 

    ! Check the version to determine what to read in 
    readTopologicalData = .false. 
    if (chardummy(8:17) >= '03.002.001') then 
        ! Topological data should be present
        readTopologicalData = .true.
    end if 
    ! call cfverr(filespec,b2fgmtryversion)

    ! Primary array dimensions
    call cfruin (filespec,6,idum,'nCi,nFc,nVx,nCg,nFs,nFt')
    nc = idum(0) ! note: only reading in actual cells, no guard cells
    nf = idum(1)
    nv = idum(2)

    ! Add to grid
    grid%cell%ntot         = nc
    grid%face%ntot         = nf
    grid%vert%ntot          = nv
    grid%data%fluxdata%nFs  = idum(4)
    grid%data%fluxdata%nFt  = idum(5)

    ! Secondary array dimensions
    call cfruin (filespec,5,idum,'nCmxVx,nCmxFc,nFmxCv,nVmxCv,nVmxFc')

    ! Add to grid
    grid%cell%nvert = idum(0)
    grid%cell%nface = idum(1)
    grid%vert%ncell  = idum(3)
    grid%vert%nface  = idum(4)

    ! Initialize topological data
    if (allocated(grid%data%isprimaryxp)) deallocate(grid%data%isprimaryxp)
    if (allocated(grid%data%xpointID)) deallocate(grid%data%xpointID)
    if (allocated(grid%data%opointID)) deallocate(grid%data%opointID)
    if (allocated(grid%data%spointID)) deallocate(grid%data%spointID)
    if (allocated(grid%data%divFcP)) deallocate(grid%data%divFcP)
    if (allocated(grid%data%divFc)) deallocate(grid%data%divFc)
    if (allocated(grid%data%spointdivID)) deallocate(grid%data%spointdivID)
    if (allocated(grid%data%tpointdivID)) deallocate(grid%data%tpointdivID)
    
    ! Read topological data
    if (readTopologicalData) then 
        ! Read topological mesh flag
        call cfruin(filespec, 1, idum, 'topoflag')
        grid%data%topoflag = idum(0)

        ! Read number of topological points
        call cfruin(filespec, 6, idum, 'nX,nO,nS,nT,Div,nDivFc')
        grid%data%nxp = idum(0)
        grid%data%nop = idum(1)
        grid%data%nsp = idum(2)
        grid%data%ntp = idum(3)
        grid%data%ndiv = idum(4)
        grid%data%ndivFc = idum(5)

        ! Allocate
        allocate(grid%data%xpointID(idum(0)), grid%data%opointID(idum(1)), &
            grid%data%spointID(idum(2)), grid%data%isprimaryxp(idum(0)), &
            grid%data%divFcP(grid%data%ndiv, 2), grid%data%divFc(grid%data%ndivFc), &
            grid%data%spointdivID(grid%data%nsp), grid%data%tpointdivID(grid%data%ntp))

        ! Read X-point data
        call ReadSingleLine(filespec, chardummy, reachedeof) ! header
        do i = 1, grid%data%nxp  
            ! Read 
            read(filespec, *) grid%data%xpointID(i), grid%data%isprimaryxp, &
                idum(0)
        end do

        ! Read O-point data
        call ReadSingleLine(filespec, chardummy, reachedeof) ! header
        do i = 1, grid%data%nop  
            ! Read 
            read(filespec, *) grid%data%opointID(i)
        end do

        ! Read strike point data
        call ReadSingleLine(filespec, chardummy, reachedeof) ! header
        do i = 1, grid%data%nsp  
            ! Read 
            read(filespec, *) grid%data%spointID(i), grid%data%spointxpID(i), &
                grid%data%spointdivID(i), idum(0)
        end do

        ! Read tangency point data
        call ReadSingleLine(filespec, chardummy, reachedeof) ! header
        do i = 1, grid%data%ntp  
            ! Read 
            read(filespec, *) grid%data%tpointID(i), grid%data%tpointdivID(i), &
                idum(0)
        end do

        ! Read divertor data
        call ReadSingleLine(filespec, chardummy, reachedeof) ! header
        do i = 1, grid%data%ndiv  
            ! Read 
            read(filespec, *) idum(0), grid%data%divFcP(i, 1), &
                grid%data%divFcP(i, 2)
        end do

        ! Read divertor face data
        call cfruin(filespec, grid%data%ndivFc, grid%data%divFc, 'divFc')

    else 
        ! Initialize to zero
        grid%data%topoflag = 0
        grid%data%nxp = 0
        grid%data%nop = 0
        grid%data%nsp = 0
        grid%data%ntp = 0
        grid%data%ndiv = 0
        grid%data%ndivFc = 0

        ! Allocate
        allocate(grid%data%xpointID(0), grid%data%opointID(0), &
            grid%data%spointID(0), grid%data%isprimaryxp(0), &
            grid%data%divFcP(0, 2), grid%data%divFc(0), &
            grid%data%spointdivID(0), grid%data%tpointdivID(0))

    end if 

    ! Allocate grid
    call AllocateGrid(grid)

    ! Read grid interconnection data
    !===============================
    ! Initialize
    !-----------
    ! Allocate dumies
    allocate(cdummy(nc,4))
    allocate(fdummy(nf,4))
    allocate(vdummy(nv,4))
    allocate(cdummyr(nc,4))
    allocate(fdummyr(nf,4))
    allocate(vdummyr(nv,4))
    allocate(fsdummyr(grid%data%fluxdata%nFs))
    allocate(facelistdummy(grid%cell%nface))
    allocate(ftdummy(grid%data%fluxdata%nFt))

    ! Read data for structured grid remapping (to be deleted in future)
    call cfruin (filespec,1,idum2,'isClassicalGrid') 
    grid%data%sglegacy%isClassicalGrid = int(idum2(1), I4) ! cast to I4 type
    if (grid%data%sglegacy%isClassicalGrid == 1) then 
        call cfruin (filespec,3,idum,'nx,ny,nncut') ! this seems to be wrongly formatted for now - to be checked in the future
        grid%data%sglegacy%nx = idum(0)
        grid%data%sglegacy%ny = idum(1)
    end if

    ! Attempt to read in OMP/IMP data 
    !--------------------------------
    ! OMP
    call ReadUntilFound(filespec, 'OMPr', reachedeof)
    if (reachedeof) then 
        ! Not found, issue warning and rewind
        print *, 'ReadTraduitUS: could not find field "OMPr", setting to zero...'
        rewind(filespec)
    else 
        ! Found, read in coordinates (assumed two)
        backspace(filespec)
        call cfrure(filespec, 2, grid%data%OMPr, 'OMPr')
    end if 
    call ReadUntilFound(filespec, 'OMPz', reachedeof)
    if (reachedeof) then 
        ! Not found, issue warning and rewind
        print *, 'ReadTraduitUS: could not find field "OMPz", setting to zero...'
        rewind(filespec)
    else 
        ! Found, read in coordinates (assumed two)
        backspace(filespec)
        call cfrure(filespec, 2, grid%data%OMPz, 'OMPz')
    end if 

    ! IMP
    call ReadUntilFound(filespec, 'IMPr', reachedeof)
    if (reachedeof) then 
        ! Not found, issue warning and rewind
        print *, 'ReadTraduitUS: could not find field "IMPr", setting to zero...'
        rewind(filespec)
    else 
        ! Found, read in coordinates (assumed two)
        backspace(filespec)
        call cfrure(filespec, 2, grid%data%IMPr, 'IMPr')
    end if 
    call ReadUntilFound(filespec, 'IMPz', reachedeof)
    if (reachedeof) then 
        ! Not found, issue warning and rewind
        print *, 'ReadTraduitUS: could not find field "IMPz", setting to zero...'
        rewind(filespec)
    else 
        ! Found, read in coordinates (assumed two)
        backspace(filespec)
        call cfrure(filespec, 2, grid%data%IMPz, 'IMPz')
    end if 

    ! Rewind and reset to expected location
    rewind(filespec)
    call ReadUntilFound(filespec, 'vxFfbz', reachedeof)
    if (reachedeof) then 
        ! Something wrong, exit
        call gdErrorHandler('ReadTraduitUS: could not reset file, check file content')
    else
        backspace(filespec)
    end if

    ! Vertex data
    !------------
    ! Here we simply loop ourselves, as cfruin etc
    ! do not seem to be capable to load in this data (?)
    allocate(vdata(nv, 6), vlist(nv))
    
    ! Skip the header
    call ReadSingleLine(filespec, chardummy, reachedeof)
    if (reachedeof) then 
        call gdErrorHandler('ReadTraduitUS: reached EOF prematurely')
    end if 

    ! Read in data
    do i = 1, nv 
        ! Read 
        read(filespec, *) vlist(i), vdata(i,:) ! apparently this works fine...
    end do

    ! Add to grid
    grid%vert%x(vlist)      = vdata(:, 1)
    grid%vert%y(vlist)      = vdata(:, 2)
    grid%vert%psi(vlist)    = vdata(:, 3)
    grid%vert%bx(vlist)     = vdata(:, 4)
    grid%vert%by(vlist)     = vdata(:, 5)
    grid%vert%ffbz(vlist)   = vdata(:, 6)

    ! Cell data
    !----------
    allocate(cdatai1(nc, 2), cdata(nc, 5), cdatai2(nc, 3), clist(nc))

    ! Skip the header
    call ReadSingleLine(filespec, chardummy, reachedeof)
    if (reachedeof) then 
        call gdErrorHandler('ReadTraduitUS: reached EOF prematurely')
    end if 

    ! Read in data
    do i = 1, nc 
        ! Read 
        read(filespec, *) clist(i), cdatai1(i,:), cdata(i, :), cdatai2(i, :) ! apparently this works fine...
    end do 

    ! Add to grid
    grid%cell%vertP(clist, 1)          = cdatai1(:, 1)
    grid%cell%vertP(clist, 2)          = cdatai1(:, 2)
    grid%cell%reg(clist)               = cdatai2(:, 2)
    grid%cell%ft(clist)                = cdatai2(:, 3)
    grid%cell%cflags(clist)            = cdatai2(:, 1)
    grid%cell%x(clist)                 = cdata(:, 1)
    grid%cell%y(clist)                 = cdata(:, 2)
    grid%cell%psi(clist)               = cdata(:, 3)
    grid%cell%bp(clist)                = cdata(:, 4)
    grid%cell%bt(clist)                = cdata(:, 5)

    ! Cell vertices and faces
    call cfruin (filespec, grid%cell%nvert, grid%cell%vert,  'cvVx')
    call cfruin (filespec, grid%cell%nface, grid%cell%face,  'cvFc')

    ! Faces
    !------
    allocate(fdatai(nf, 5), flist(nf))

    ! Skip the header
    call ReadSingleLine(filespec, chardummy, reachedeof)
    if (reachedeof) then 
        call gdErrorHandler('ReadTraduitUS: reached EOF prematurely')
    end if 

    ! Read in data
    do i = 1, nf 
        ! Read 
        read(filespec, *) flist(i), fdatai(i,:) ! apparently this works fine...
    end do 

    ! Add to grid
    grid%face%vert(flist, 1)            = fdatai(:, 1)
    grid%face%vert(flist, 2)            = fdatai(:, 2)
    grid%face%label(flist)              = fdatai(:, 3)
    grid%face%reg(flist)                = fdatai(:, 4)
    grid%face%aligned                   = fdatai(:, 5) 

    ! Flux tubes
    !-----------
    allocate(ftdatai(grid%data%fluxdata%nFt, 5), &
        ftlist(grid%data%fluxdata%nFt))

    ! Skip the header
    call ReadSingleLine(filespec, chardummy, reachedeof)
    if (reachedeof) then 
        call gdErrorHandler('ReadTraduitUS: reached EOF prematurely')
    end if 

    ! Read in data
    do i = 1, grid%data%fluxdata%nFt 
        ! Read 
        read(filespec, *) ftlist(i), ftdatai(i,:) ! apparently this works fine...
    end do

    ! Compute the number of ftCv and ftFc that are actually there
    nftCv = ftdatai(grid%data%fluxdata%nFt, 2) + &
        ftdatai(grid%data%fluxdata%nFt, 1)-1
    nftFc = ftdatai(grid%data%fluxdata%nFt, 4) + &
        ftdatai(grid%data%fluxdata%nFt, 3)-1
    
    ! Read ftCv, ftFc
    allocate(ftCvdum(nftCv), ftFcdum(nftFc))
    call cfruin (filespec, nftCv, ftCvdum, 'ftCv')
    call cfruin (filespec, nftFc, ftFcdum, 'ftFc')
    

    ! Add to grid
    grid%data%fluxdata%fluxtubecells(1:nftCv) = ftCvdum
    grid%data%fluxdata%fluxtubefaces(1:nftFc) = ftFcdum
    grid%data%fluxdata%fluxtubecellsP(ftlist, 1) = ftdatai(:, 1)
    grid%data%fluxdata%fluxtubecellsP(ftlist, 2) = ftdatai(:, 2)
    grid%data%fluxdata%fluxtubefacesP(ftlist, 1) = ftdatai(:, 3)
    grid%data%fluxdata%fluxtubefacesP(ftlist, 2) = ftdatai(:, 4)
    grid%data%fluxdata%fluxtuberegID(ftlist) = ftdatai(:, 5)

    ! Flux surfaces
    !--------------
    allocate(fsdatai(grid%data%fluxdata%nFs, 2), &
        fsdata(grid%data%fluxdata%nFs, 1), &
        fslist(grid%data%fluxdata%nFs)) 

    ! Skip the header
    call ReadSingleLine(filespec, chardummy, reachedeof)
    if (reachedeof) then 
        call gdErrorHandler('ReadTraduitUS: reached EOF prematurely')
    end if 

    ! Read in data
    do i = 1, grid%data%fluxdata%nFs 
        ! Read 
        read(filespec, *) fslist(i), fsdatai(i,:), fsdata(i, :) ! apparently this works fine...
    end do

    ! Compute the number of fsFc
    nfsFc = fsdatai(grid%data%fluxdata%nFs, 2) + &
        fsdatai(grid%data%fluxdata%nFs, 1)-1

    ! Read fsFc
    allocate(fsFcdum(nfsFc))
    call cfruin (filespec, nfsFc, fsFcdum ,  'fsFc')

    ! Add to grid
    grid%data%fluxdata%fluxsurfacefaces(1:nfsFc) = fsFcdum
    grid%data%fluxdata%fluxsurfacefacesP(fslist, 1) = fsdatai(:, 1)
    grid%data%fluxdata%fluxsurfacefacesP(fslist, 2) = fsdatai(:, 2)
    grid%data%fluxdata%fluxsurfacepsi(fslist)       = fsdata(:, 1)
    ! fsdata(:, 1) contains psi value of flux surfaces, to be read in in the future?
    
    ! Housekeeping
    !=============
    deallocate(vdata, cdata, fsdata, vlist, flist, clist, ftlist, &
        fslist, cdatai1, cdatai2, fdatai, ftdatai, fsdatai, ftCvdum, &
        ftFcdum, fsFcdum)
    close(filespec)

end subroutine