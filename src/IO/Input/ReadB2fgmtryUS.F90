subroutine ReadB2fgmtryUS(grid, filepath)

    ! Description
    !============
    ! Read in the necessary grid data from a b2fgmtry_us file for the 
    ! grid deformation module. This routine can serve as wrapper in the 
    ! future to call the dedicated reading and writing routines 
    ! present in b2mod_geo for example. 
    
    ! Notes
    !======
    ! Note 1: typically, there are still spurious x-point vertices and 
    ! guard cells in the b2fgmtry format. These are removed here. To 
    ! this end, we first read in the grid in a classical way. Then, 
    ! we check which vertices and cells have to be removed and 
    ! construct a grid without guard cells etc. It is assumed that all
    ! guard cells are added at the end of the cell list

    ! Note 2: contrary to the traduit.out.b2us file, we cannot assume
    ! that the cell vertices are ordened either clockwise or 
    ! counterclockwise. This necessitates an additional sorting loop
    ! to comply with the assumptions made in GOAT. (we could at some
    ! point include this sorting in the grid interconnection computation
    ! routine, but this would lead to unnecessary overhead for traduit
    ! files)

    ! Initialize
    !===========
    ! Declare modules
    use gdmod_types 
    use mod_polygon

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Main variables
    type(GridUDT)               :: grid
    character(*), intent(in)    :: filepath

    ! Auxiliary
    type(GridUDT)               :: grido
    integer(I8)                 :: idum(0:9), idum2(1), filespec   
    character(10)               :: b2fgmtryversion
    integer(I8)                 :: nc,nf,nv ! total number of cells, faces, vertices

    character(120)              :: chardummy   ! dummy array
    integer(I8), allocatable    :: cdummy(:,:), cdummy2(:) ! dummy array
    integer(I8), allocatable    :: fdummy(:,:), fdummy2(:, :) ! dummy array
    integer(I8), allocatable    :: vdummy(:,:) ! dummy array
    integer(I8), allocatable    :: ftdummy(:) ! dummy array

    real(R8), allocatable       :: fcQalf(:, :) ! to reconstruct aligned faces
    real(R8), allocatable       :: cdummyr(:,:) ! dummy array
    real(R8), allocatable       :: fdummyr(:,:) ! dummy array
    real(R8), allocatable       :: vdummyr(:,:) ! dummy array
    real(R8), allocatable       :: fsdummyr(:) ! dummy array
    real(R8), allocatable       :: facedummy(:) ! dummy array
    integer(I8), allocatable       :: n2dummy(:) ! dummy array
    integer(I8), allocatable       :: nxdummy(:) ! dummy array

    integer(I8)                 :: n2,nx  ! legacy structured data

    integer(I8)                                 :: ngv, ngc 
    integer(I8), allocatable, dimension(:)      :: vertmap, cellmap, &
        sortindex, tcf, tv, tftc
    integer(I8), allocatable, dimension(:, :)   :: tempvertfaceP, &
        tempvertcellP, tempcellfaceP, tempcellvertP, tcfv

    logical, allocatable, dimension(:)      :: isnoghostvert, keepvertface, &
        keepvertcell, isnoguardcell, keepcellface, keepcellvert, &
        ispolygonstart 

    ! Loop
    integer(I8)                 :: i, j, k
    
    ! Read (inlucding guard cells)
    !=============================
    ! Open the file
    print *, 'reading grid in b2gmtry_us format from file: ' // filepath
    open(unit = filespec, file = filepath)

    ! First, read the header with the version
    call cfverr(filespec,b2fgmtryversion)

    ! Primary array dimensions
    call cfruin (filespec,7,idum,'nCi,nCg,nCv,nFc,nVx,nFs,nFt')
    ngc = int(idum(1), I8)
    nc = int(idum(2), I8) ! make sure to cast to correct type
    nf = int(idum(3), I8)
    nv = int(idum(4), I8)

    ! Add to grid
    grid%cell%ntot         = nc
    grid%face%ntot         = nf
    grid%vert%ntot          = nv
    grid%data%fluxdata%nFs  = idum(5)
    grid%data%fluxdata%nFt  = idum(6)

    ! Secondary array dimensions
    call cfruin (filespec,5,idum,'nCmxVx,nCmxFc,nVmxCv,nVmxFc,nCmxNv')

    ! Add to grid
    grid%cell%nvert = idum(0)
    grid%cell%nface = idum(1)
    grid%vert%ncell  = idum(2)
    grid%vert%nface  = idum(3)

    ! Allocate grid
    call AllocateGrid(grid)

    ! Check
    if (.not. allocated(grid%vert%face)) then 
        allocate(grid%vert%face(grid%vert%nface))
    end if 
    if (.not. allocated(grid%vert%cell)) then 
        allocate(grid%vert%cell(grid%vert%ncell))
    end if 

    ! Allocate dumies
    allocate(cdummy(nc,4), cdummy2(nc))
    allocate(fdummy(nf,4))
    allocate(vdummy(nv,4))
    allocate(cdummyr(nc,4))
    allocate(fdummyr(nf,4))
    allocate(vdummyr(nv,4))
    allocate(fsdummyr(grid%data%fluxdata%nFs))
    allocate(facedummy(grid%cell%nface))
    allocate(ftdummy(grid%data%fluxdata%nFt))
    allocate(fcQalf(nf, 2))

    ! Read data for structured grid remapping (to be deleted in future)
    call cfruin (filespec,1,idum2,'isClassicalGrid')
    grid%data%sglegacy%isClassicalGrid = int(idum2(1), I4)
    call cfruin (filespec,3,idum,'nx,ny,nncut')
    grid%data%sglegacy%nx       = idum(0)
    grid%data%sglegacy%ny       = idum(1)
    grid%data%sglegacy%nncut = idum(2)
    
    ! Read data that is not used
    call cfruch (filespec,120,chardummy,'label')
    call cfruin (filespec,1,idum, 'isymm')

    ! Add grid mapping data
    allocate(fdummy2(nf, 2))
    call cfruin (filespec, nc*2,    grid%cell%faceP, 'cvFcP')
    call cfruin (filespec, grid%cell%nface, grid%cell%face,  'cvFc')
    call cfruin (filespec, nf*2,    fdummy2,  'fcCv')
    call cfruin (filespec, nf*2,    grid%face%vert,  'fcVx')
    call cfruin (filespec, nc*2,    grid%cell%vertP, 'cvVxP')
    call cfruin (filespec, grid%cell%nvert, grid%cell%vert,  'cvVx')
    call cfruin (filespec, nv*2,     grid%vert%faceP, 'vxFcP')
    call cfruin (filespec, grid%vert%nface,  grid%vert%face,  'vxFc')
    call cfruin (filespec, nv*2,     grid%vert%cellP, 'vxCvP')
    call cfruin (filespec, grid%vert%ncell,  grid%vert%cell,  'vxCv')
    call cfruin (filespec, grid%data%fluxdata%nFt*2,   grid%data%fluxdata%fluxtubecellsP, 'ftCvP')
    call cfruin (filespec, nc,     grid%data%fluxdata%fluxtubecells,  'ftCv')
    call cfruin (filespec, grid%data%fluxdata%nFt*2,   grid%data%fluxdata%fluxtubefacesP, 'ftFcP')
    call cfruin (filespec, nf,     grid%data%fluxdata%fluxtubefaces,  'ftFc')
    call cfruin (filespec, nc,     grid%cell%ft,  'cvFt')   
    call cfruin (filespec, grid%data%fluxdata%nFs*2,   grid%data%fluxdata%fluxsurfacefacesP, 'fsFcP')
    call cfruin (filespec, nf,     grid%data%fluxdata%fluxsurfacefaces,  'fsFc')
    call cfruin (filespec, nf,     grid%face%reg, 'fcReg')
    call cfruin (filespec, nc,     grid%cell%reg, 'cvReg')
    call cfruin (filespec, grid%data%fluxdata%nFt,     grid%data%regions%fluxtuberegID, 'ftReg')
    call cfrure (filespec, grid%cell%nface,  facedummy,'intcellP') ! not used
    call cfrure (filespec, grid%cell%nface,  facedummy,'intcellR') ! not used
    deallocate(fdummy2)

    ! Add additional data from underlying structured grid (to be removed)
    n2 = 0
    nx = 0
    if (grid%data%sglegacy%isClassicalGrid == 1) n2 = (grid%data%sglegacy%nx+2)*(grid%data%sglegacy%ny+2)
    if (grid%data%sglegacy%isClassicalGrid == 1) nx =  grid%data%sglegacy%nx

    ! For now, put it into dummies
    allocate(n2dummy(n2))
    allocate(nxdummy(nx))
    call cfruin (filespec, n2, n2dummy,  'imapCv')
    call cfruin (filespec, n2, n2dummy, 'imapFcx')
    call cfruin (filespec, n2, n2dummy, 'imapFcy')
    call cfruin (filespec, n2, n2dummy,  'imapVx')
    call cfruin (filespec, nx, nxdummy, 'icornVx')

    ! Read some labels, ignore
    call cfruin (filespec, nf, grid%face%label,'fcLbl')
    call cfruin (filespec, nc, cdummy,'cvLbl')
    call cfruin (filespec, grid%data%fluxdata%nFt, ftdummy,'ftLbl')

    ! Add geometry data
    ! cell data - ignore all
    call cfrure (filespec, nc*4, cdummyr(:,1:4),   'cvBb')
    call cfrure (filespec, nc*3, cdummyr(:,1:3),   'cvEb')
    call cfrure (filespec, nc,   cdummyr(:,1),    'cvX')
    call cfrure (filespec, nc,   cdummyr(:,1),    'cvY')
    call cfrure (filespec, nc,   cdummyr(:,1),   'cvSz')
    call cfrure (filespec, nc,   cdummyr(:,1),   'cvHz')
    call cfrure (filespec, nc,   cdummyr(:,1),   'cvHx') !WG temp!
    call cfrure (filespec, nc*2, cdummyr(:,1:2), 'cvQgam')
    call cfrure (filespec, nc,   cdummyr(:,1),  'cvVol')

    ! face quantities - ignore all except fcQalf
    call cfrure (filespec, nf*4, fdummyr(:,1:4),   'fcBb')
    call cfrure (filespec, nf,   fdummyr(:,1),    'fcS')
    call cfrure (filespec, nf*2, fdummyr(:,1:2),   'fcHc')
    call cfrure (filespec, nf,   fdummyr(:,1),   'fcHt')
    call cfrure (filespec, nf*2, fdummyr(:,1:2), 'fcQgam')
    call cfrure (filespec, nf*2, fcQalf, 'fcQalf')
    call cfrure (filespec, nf*2, fdummyr(:,1:2), 'fcQbet')
    call cfrure (filespec, nf,   fdummyr(:,1),  'fcPbs')

    ! Determine aligned faces as those faces with cos(alpha) = 0
    ! But do it with a very small tolerance to avoid numerical garbage
    where (fcQalf(:, 1) < 1e-10) 
        grid%face%aligned = 1
    elsewhere
        grid%face%aligned = 0
    end where


    ! vertex quantities - only keep coordinates (and ffbz)
    call cfrure (filespec, nv*4, vdummyr(:,1:4),   'vxBb')
    call cfrure (filespec, nv,   grid%vert%x,    'vxX')
    call cfrure (filespec, nv,   grid%vert%y,    'vxY')
    call cfrure (filespec, nv,   grid%vert%ffbz, 'vxFfbz')
    call cfrure (filespec, nv,   vdummyr(:,1), 'vxFpsi')

    ! flux surface quantities
    call cfrure (filespec, nc,   cdummyr(:,1), 'cvConn')
    call cfrure (filespec, grid%data%fluxdata%nFs,   &
        grid%data%fluxdata%fluxsurfacepsi, 'fsPsi')

    ! Eliminate 'ghost' vertices
    !===========================
    ! Allocate & initialize
    allocate(isnoghostvert(grid%vert%ntot))
    isnoghostvert(:) = .false.

    ! Loop
    do j = 1, 2
        do i = 1, grid%face%ntot
            isnoghostvert(grid%face%vert(i,j)) = .true.
        enddo
    enddo

    ! Count number of ghost vertices
    ngv = count(.not. isnoghostvert) ! total number of ghost vertices

    ! Construct vertex mapping for vertex arrays
    allocate(vertmap(grid%vert%ntot))
    vertmap = 0
    k = 0
    do i = 1, grid%vert%ntot
        ! Skip if ghostvert
        if (isnoghostvert(i)) then 

            ! Update k
            k = k + 1

            ! Set mapping
            vertmap(i) = k 
        end if
    end do

    ! Store original grid
    grido = grid 

    ! Rebuild vertex fields
    grid%vert%ntot = grid%vert%ntot - ngv
    grid%vert%x = pack(grid%vert%x, isnoghostvert)
    grid%vert%y = pack(grid%vert%y, isnoghostvert)
    grid%vert%fieldlineID = pack(grid%vert%fieldlineID, isnoghostvert)
    grid%vert%bx = pack(grid%vert%bx, isnoghostvert)
    grid%vert%by = pack(grid%vert%by, isnoghostvert)
    grid%vert%BV = pack(grid%vert%BV, isnoghostvert)
    grid%vert%ffbz = pack(grid%vert%ffbz, isnoghostvert)
    grid%vert%psi = pack(grid%vert%psi, isnoghostvert)
    if (allocated(grid%vert%neigP)) then 
        deallocate(grid%vert%neigP)
        allocate(grid%vert%neigP(grid%vert%ntot, 2))
    end if 

    ! Rebuild face fields
    grid%face%vert(:, 1) = vertmap(grid%face%vert(:, 1))
    grid%face%vert(:, 2) = vertmap(grid%face%vert(:, 2))

    ! Rebuild cell fields
    grid%cell%vert = vertmap(grid%cell%vert)

    ! Sanity checks
    if (any(grid%cell%vert == 0)) then 
        ! Zero vertices in cells detected, shouldn't happen
        call gdErrorHandler('ReadB2fgmtryUS: cells with zero vertex ID ' // & 
            'detected after removing ghost vertices. Check input grid')
    end if 
    if (any(any(grid%face%vert == 0, 1))) then 
        ! Zero vertices in faces detected, shouldn't happen
        call gdErrorHandler('ReadB2fgmtryUS: faces with zero vertex ID' // & 
            'detected after removing ghost vertices. Check input grid')
    end if 

    ! Rebuild vertex interconnection data
    allocate(keepvertface(grid%vert%nface), keepvertcell(grid%vert%ncell), &
        tempvertfaceP(grid%vert%ntot, 2), tempvertcellP(grid%vert%ntot, 2))
    keepvertface = .true. 
    keepvertcell = .true. 
    where (grid%vert%face == 0) keepvertface = .false. ! eliminate zero values 
    where (grid%vert%cell == 0) keepvertcell = .false.
    do i = 1, grid%vert%ntot + ngv
        if (.not. isnoghostvert(i)) then 
            ! Set to false
            keepvertface(grid%vert%faceP(i, 1):grid%vert%faceP(i, 1)+grid%vert%faceP(i, 2)-1) = .false.
            keepvertcell(grid%vert%cellP(i, 1):grid%vert%cellP(i, 1)+grid%vert%cellP(i, 2)-1) = .false.
        end if
    end do 
    tempvertfaceP(:, 1) = pack(grid%vert%faceP(:, 1), isnoghostvert)
    tempvertfaceP(:, 2) = pack(grid%vert%faceP(:, 2), isnoghostvert)
    tempvertcellP(:, 1) = pack(grid%vert%cellP(:, 1), isnoghostvert)
    tempvertcellP(:, 2) = pack(grid%vert%cellP(:, 2), isnoghostvert)
    grid%vert%faceP = tempvertfaceP 
    grid%vert%faceP(1, 1)  = 1
    grid%vert%cellP = tempvertcellP 
    grid%vert%cellP(1, 1) = 1
    do i = 2, grid%vert%ntot
        grid%vert%faceP(i, 1) = grid%vert%faceP(i-1, 1) + grid%vert%faceP(i-1, 2)
        grid%vert%cellP(i, 1) = grid%vert%cellP(i-1, 1) + grid%vert%cellP(i-1, 2)
    end do 
    grid%vert%face = pack(grid%vert%face, keepvertface)
    grid%vert%cell = pack(grid%vert%cell, keepvertcell)
    grid%vert%nface = size(grid%vert%face)
    grid%vert%ncell = size(grid%vert%cell)

    ! Eliminate guard cells
    !======================
    ! Allocate & initialize
    allocate(isnoguardcell(grid%cell%ntot))
    isnoguardcell(:) = .false.
    isnoguardcell(1:grid%cell%ntot-ngc) = .true.

    ! Construct cell mapping for cell arrays
    allocate(cellmap(grid%cell%ntot))
    cellmap = 0
    k = 0
    do i = 1, grid%cell%ntot
        ! Skip if ghostvert
        if (isnoguardcell(i)) then 

            ! Update k
            k = k + 1

            ! Set mapping
            cellmap(i) = k 
        end if
    end do

    ! Rebuild vertex fields
    grid%vert%cell = pack(grid%vert%cell, isnoguardcell(grid%vert%cell))
    grid%vert%cell = cellmap(grid%vert%cell)

    ! Rebuild face fields (nothing to be done)

    ! Rebuild cell fields
    grid%cell%ntot = grid%cell%ntot - ngc 
    grid%cell%ft = pack(grid%cell%ft, isnoguardcell)
    grid%cell%reg = pack(grid%cell%reg, isnoguardcell)
    
    ! Rebuild flux tube data
    grid%data%fluxdata%fluxtubecells = cellmap(grid%data%fluxdata%fluxtubecells)
    do i = 1, grid%data%fluxdata%nft 
        tftc = GetFTCell(grid%data%fluxdata, i)
        grid%data%fluxdata%fluxtubecellsP(i, 2) = grid%data%fluxdata%fluxtubecellsP(i, 2) - count(tftc == 0)
    end do 
    grid%data%fluxdata%fluxtubecells = &
        pack(grid%data%fluxdata%fluxtubecells, grid%data%fluxdata%fluxtubecells /= 0)
    grid%data%fluxdata%fluxtubecellsP(1, 1) = 1
    grid%data%fluxdata%fluxtubefacesP(1, 1) = 1 ! hedge for junk here from input
    do i = 2, grid%data%fluxdata%nft 
        grid%data%fluxdata%fluxtubecellsP(i, 1) = &
            grid%data%fluxdata%fluxtubecellsP(i-1, 1) + grid%data%fluxdata%fluxtubecellsP(i-1, 2)
        grid%data%fluxdata%fluxtubefacesP(i, 1) = &
            grid%data%fluxdata%fluxtubefacesP(i-1, 1) + grid%data%fluxdata%fluxtubefacesP(i-1, 2)
    end do 


    ! Rebuild cell interconnection data
    allocate(keepcellface(grid%cell%nface), keepcellvert(grid%cell%nvert), &
        tempcellfaceP(grid%cell%ntot, 2), tempcellvertP(grid%cell%ntot, 2))
    keepcellface = .true. 
    keepcellvert = .true.
    where (grid%cell%face == 0) keepcellface = .false.
    where (grid%cell%vert == 0) keepcellvert = .false. 
    do i = 1, grid%cell%ntot + ngc
        if (.not. isnoguardcell(i)) then 
            ! Set to false
            keepcellface(grid%cell%faceP(i, 1):grid%cell%faceP(i, 1)+grid%cell%faceP(i, 2)-1) = .false.
            keepcellvert(grid%cell%vertP(i, 1):grid%cell%vertP(i, 1)+grid%cell%vertP(i, 2)-1) = .false.
        end if
    end do 
    tempcellfaceP(:, 1) = pack(grid%cell%faceP(:, 1), isnoguardcell)
    tempcellfaceP(:, 2) = pack(grid%cell%faceP(:, 2), isnoguardcell)
    tempcellvertP(:, 1) = pack(grid%cell%vertP(:, 1), isnoguardcell)
    tempcellvertP(:, 2) = pack(grid%cell%vertP(:, 2), isnoguardcell)
    grid%cell%faceP = tempcellfaceP 
    grid%cell%faceP(1, 1)  = 1
    grid%cell%vertP = tempcellvertP 
    grid%cell%vertP(1, 1) = 1
    do i = 2, grid%cell%ntot
        grid%cell%faceP(i, 1) = grid%cell%faceP(i-1, 1) + grid%cell%faceP(i-1, 2)
        grid%cell%vertP(i, 1) = grid%cell%vertP(i-1, 1) + grid%cell%vertP(i-1, 2)
    end do 
    grid%cell%face = pack(grid%cell%face, keepcellface)
    grid%cell%vert = pack(grid%cell%vert, keepcellvert)
    grid%cell%nface = size(grid%cell%face)
    grid%cell%nvert = size(grid%cell%vert)

    ! Reorden cell vertices
    !======================
    ! Since we know the cell faces from the format, we can do this 
    ! easily by using the SortPolygonEdges and ExtractPolygonVertices 
    ! routines
    do i = 1, grid%cell%ntot 
        ! Get cell faces
        tcf = GetCellFace(grid%cell, i)

        ! Get vertices
        tcfv = grid%face%vert(tcf, :)

        ! Sort
        allocate(sortindex(size(tcf)), ispolygonstart(size(tcf)))
        call SortPolygonEdges(tcfv, size(tcf), sortindex, ispolygonstart)
        tcfv = tcfv(sortindex, :)

        ! Check
        if ( (count(ispolygonstart) /= 1)) then 
            call gdErrorHandler('ReadB2fgmtry: multiple polygons detected ' // &
                'for single cell, not supported. Check grid input')
        end if 

        ! Extract vertices
        allocate(tv(size(tcf)+1))
        call ExtractPolygonVertices(tcfv, size(tcf), tv)

        ! Check
        if (tv(1) /= tv(size(tcf)+1)) then 
            call gdErrorHandler('ReadB2fgmtry: cell vertices do not form ' // &
                'closed polygon, not supported. Check grid input')
        end if 

        ! Add to grid vertices
        grid%cell%vert(grid%cell%vertP(i, 1):grid%cell%vertP(i, 1)+grid%cell%vertP(i, 2)-1) = &
            tv(1:size(tcf))

        ! Housekeeping
        deallocate(tv, sortindex, ispolygonstart)

    end do

    ! Housekeeping
    !=============
    ! Close file
    close(filespec)

end subroutine