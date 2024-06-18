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

    logical                     :: reachedeof

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

    ! First, read the header with the versionv - just ignore that alltogether
    call ReadSingleLine(filespec, chardummy, reachedeof)
    if (reachedeof) then 
        call gdErrorHandler('ReadTraduitUS: reached EOF prematurely')
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
    grid%vert%x(vlist) = vdata(:, 1)
    grid%vert%y(vlist) = vdata(:, 2)

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
    grid%data%regions%fluxtuberegID(ftlist) = ftdatai(:, 5)

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