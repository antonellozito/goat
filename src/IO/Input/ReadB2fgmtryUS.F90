subroutine ReadB2fgmtryUS(filespec,grid)

    ! Description
    !============
    ! Read in the necessary grid data from a b2fgmtry_us file for the 
    ! grid deformation module. This routine can serve as wrapper in the 
    ! future to call the dedicated reading and writing routines 
    ! present in b2mod_geo for example. 

    ! Initialize
    !===========
    ! Declare modules
    use gdmod_types 

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Main variables
    type(GridUDT)               :: grid

    integer, intent(in)         :: filespec ! file specifier
    integer(I8)                 :: idum(0:9), idum2(1)    
    character(10)               :: b2fgmtryversion
    integer(I8)                 :: nc,nf,nv ! total number of cells, faces, vertices

    character(120)              :: chardummy   ! dummy array
    integer(I8), allocatable    :: cdummy(:,:) ! dummy array
    integer(I8), allocatable    :: fdummy(:,:), fdummy2(:, :) ! dummy array
    integer(I8), allocatable    :: vdummy(:,:) ! dummy array
    integer(I8), allocatable    :: ftdummy(:) ! dummy array

    real(R8), allocatable       :: cdummyr(:,:) ! dummy array
    real(R8), allocatable       :: fdummyr(:,:) ! dummy array
    real(R8), allocatable       :: vdummyr(:,:) ! dummy array
    real(R8), allocatable       :: fsdummyr(:) ! dummy array
    real(R8), allocatable       :: facedummy(:) ! dummy array
    integer(I8), allocatable       :: n2dummy(:) ! dummy array
    integer(I8), allocatable       :: nxdummy(:) ! dummy array

    integer(I8)                 :: n2,nx  ! legacy structured data
    
    ! Read 
    !=====
    ! First, read the header with the version
    call cfverr(filespec,b2fgmtryversion)

    ! Primary array dimensions
    call cfruin (filespec,7,idum,'nCi,nCg,nCv,nFc,nVx,nFs,nFt')
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

    ! Allocate dumies
    allocate(cdummy(nc,4))
    allocate(fdummy(nf,4))
    allocate(vdummy(nv,4))
    allocate(cdummyr(nc,4))
    allocate(fdummyr(nf,4))
    allocate(vdummyr(nv,4))
    allocate(fsdummyr(grid%data%fluxdata%nFs))
    allocate(facedummy(grid%cell%nface))
    allocate(ftdummy(grid%data%fluxdata%nFt))

    ! Read data for structured grid remapping (to be deleted in future)
    call cfruin (filespec,1,idum2,'isClassicalGrid')
    grid%data%sglegacy%isClassicalGrid = int(idum2(1), I4)
    call cfruin (filespec,3,idum,'nx,ny,nncut')
    grid%data%sglegacy%nx = idum(0)
    grid%data%sglegacy%ny = idum(1)
    
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
    call cfruin (filespec, nc,     grid%data%fluxdata%cellfluxtubeID,  'cvFt')   
    call cfruin (filespec, grid%data%fluxdata%nFs*2,   grid%data%fluxdata%fluxsurfacefacesP, 'fsFcP')
    call cfruin (filespec, nf,     grid%data%fluxdata%fluxsurfacefaces,  'fsFc')
    call cfruin (filespec, nf,     grid%data%regions%faceregID, 'fcReg')
    call cfruin (filespec, nc,     grid%data%regions%cellregID, 'cvReg')
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
    call cfruin (filespec, nf, grid%data%regions%facelabel,'fcLbl')
    call cfruin (filespec, nc, grid%data%regions%celllabel,'cvLbl')
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

    ! face quantities - ignore all
    call cfrure (filespec, nf*4, fdummyr(:,1:4),   'fcBb')
    call cfrure (filespec, nf,   fdummyr(:,1),    'fcS')
    call cfrure (filespec, nf*2, fdummyr(:,1:2),   'fcHc')
    call cfrure (filespec, nf,   fdummyr(:,1),   'fcHt')
    call cfrure (filespec, nf*2, fdummyr(:,1:2), 'fcQgam')
    call cfrure (filespec, nf*2, fdummyr(:,1:2), 'fcQalf')
    call cfrure (filespec, nf*2, fdummyr(:,1:2), 'fcQbet')
    call cfrure (filespec, nf,   fdummyr(:,1),  'fcPbs')

        ! vertex quantities - only keep coordinates
    call cfrure (filespec, nv*4, vdummyr(:,1:4),   'vxBb')
    call cfrure (filespec, nv,   grid%vert%x,    'vxX')
    call cfrure (filespec, nv,   grid%vert%y,    'vxY')
    call cfrure (filespec, nv,   vdummyr(:,1), 'vxFfbz')
    call cfrure (filespec, nv,   vdummyr(:,1), 'vxFpsi')

    ! flux surface quantities
    call cfrure (filespec, nc,   cdummyr(:,1), 'cvConn')
    call cfrure (filespec, grid%data%fluxdata%nFs,   fsdummyr,  'fsPsi')

end subroutine