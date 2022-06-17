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
    integer                 :: idum(0:9)    
    character(10)           :: b2fgmtryversion

    ! Read 
    !=====
    ! First, read the header with the version
    call cfverr(filespec,b2fgmtryversion)

    ! Array dimensions
    call cfruin (filespec,7,idum,'nCi,nCg,nCv,nFc,nVx,nFs,nFt')

    ! Add to grid
    grid%cells%ntot = idum(0) + idum(1)

    ! Add to grid
    !============

    


end subroutine