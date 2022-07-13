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
    type(GridUDT)           :: grid

    integer, intent(in)     :: filespec ! file specifier
    integer(I8)             :: idum(0:9)    
    character(10)           :: b2fgmtryversion

    character(120)          :: chardummy   ! dummy array
    integer, allocatable    :: cdummy(:,:) ! dummy array
    integer, allocatable    :: fdummy(:,:) ! dummy array
    integer, allocatable    :: vdummy(:,:) ! dummy array

    ! Read 
    !=====
    ! First, read the header with the version
    call cfverr(filespec,b2fgmtryversion)

    ! Primray array dimensions
    call cfruin (filespec,7,idum,'nCi,nCg,nCv,nFc,nVx,nFs,nFt')

    ! Other things we don't need
    !call cfruin (ninp(1),5,idum,'nCmxVx,nCmxFc,nVmxCv,nVmxFc,nCmxNv')
    !call cfruin (ninp(1),1,idum,'isClassicalGrid')
    !call cfruin (ninp(1),3,idum,'nx,ny,nncut')

    ! Add to grid
    grid%cells%ntot = idum(2)
    grid%faces%ntot = idum(3)
    grid%vert%ntot  = idum(4)

    ! Secondary array dimensions
    call cfruin (filespec,5,idum,'nCmxVx,nCmxFc,nVmxCv,nVmxFc,nCmxNv')

    ! Add to grid
    grid%cells%nvertlist = idum(0)
    grid%cells%nfacelist = idum(1)
    grid%vert%ncelllist  = idum(2)
    grid%vert%nfacelist  = idum(3)

    ! Allocate grid
    call AllocateGrid(grid)

    ! Allocate dumies
    allocate(cdummy(grid%cells%ntot,2))
    allocate(fdummy(grid%faces%ntot,2))
    allocate(vdummy(grid%vert%ntot,2))

    ! Continue reading (data not required)
    call cfruin (filespec,1,idum,'isClassicalGrid') 
    call cfruin (filespec,3,idum,'nx,ny,nncut')
    call cfruch (filespec,120,chardummy,'label')
    call cfruin (filespec,1,idum, 'isymm')

    ! Add cell data
    call cfruin (filespec, grid%cells%ntot*2,    grid%cells%faceP, 'cvFcP')
    call cfruin (filespec, grid%cells%nfacelist, grid%cells%facelist,  'cvFc')
    call cfruin (filespec, grid%faces%ntot*2,    grid%faces%neig,  'fcCv')
    call cfruin (filespec, grid%faces%ntot*2,    grid%faces%vert,  'fcVx')
    call cfruin (filespec, grid%cells%ntot*2,    grid%cells%vertP, 'cvVxP')
    call cfruin (filespec, grid%cells%nvertlist, grid%cells%vertlist,  'cvVx')
    call cfruin (filespec, grid%vert%ntot*2,     grid%vert%faceP, 'vxFcP')
    call cfruin (filespec, grid%vert%nfacelist,  grid%vert%facelist,  'vxFc')
    call cfruin (filespec, grid%vert%ntot*2,     grid%vert%cellP, 'vxCvP')
    call cfruin (filespec, grid%vert%ncelllist,  grid%vert%celllist,  'vxCv')
    !call cfruin (filespec, grid%flu%nFt*2,   mpg%ftCvP, 'ftCvP')
    !call cfruin (filespec, grid%cells%ntot,     mpg%ftCv,  'ftCv')
    !call cfruin (filespec, mpg%nFt*2,   mpg%ftFcP, 'ftFcP')
    !call cfruin (filespec, grid%faces%ntot,     mpg%ftFc,  'ftFc')
    !call cfruin (filespec, grid%cells%ntot,     mpg%cvFt,  'cvFt')   
    !call cfruin (filespec, mpg%nFs*2,   mpg%fsFcP, 'fsFcP')
    !call cfruin (filespec, grid%faces%ntot,     mpg%fsFc,  'fsFc')
    !call cfruin (filespec, grid%faces%ntot,     mpg%fcReg, 'fcReg')
    !call cfruin (filespec, grid%cells%ntot,     mpg%cvReg, 'cvReg')
    !call cfruin (filespec, mpg%nFt,     mpg%ftReg, 'ftReg')
    !call cfrure (filespec, mpg%nCmxFc,  mpg%intcellP,'intcellP')
    !call cfrure (filespec, mpg%nCmxFc,  mpg%intcellR,'intcellR')


    !call cfrure (filespec, grid%vert%ntot,   grid%vert%x,    'vxX')
    !call cfrure (filespec, grid%vert%ntot,   grid%vert%y,    'vxY')

    

    ! Add to grid
    !============

    


end subroutine