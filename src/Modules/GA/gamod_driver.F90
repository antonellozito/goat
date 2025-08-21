!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains all routines related to the driver for the grid adaptation module, such as initialization and postprocessing routines.

module gamod_driver

    ! Initialize
    !============
    ! Load modules
    use goatmod_types
    use goatmod_userinput
    use gamod_utility
    use gamod_types

    ! The usual
    implicit none
    save
    public     
    
    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!  
    
    contains

    subroutine GridAdaptor(grid,environment,magneticField,options)

        ! Description
        !============
        ! Overarching driver for grid adaptation, one lever lower than GADriver. ! Adapts topology of unstructured grid with the goal to improve grid  quality based on grid metric and user inputs

        ! Initialize
        !===========
        ! Declare modules  
        
        ! Declare variables
        !==================
        ! Arguments
        type(GAGridUDT), intent(inout)              :: grid
        type(EnvironmentUDT), intent(in)            :: environment
        type(MagneticFieldUDT), intent(in)          :: magneticField     
        
        type(GAoptionsUDT), intent(in)              :: options

        ! Initialize grid adaptation
        !===========================
        call GAinit(grid,options,environment,magneticField)
        

        ! Driver Selection
        !=================
        select case (options%meth)

        case ('simple')

            ! Regular grid adaption
            call GAInternalDriver(grid,options,environment,magneticField)

        case ('aposteriori')

            ! Grid adaptation based on simulation information
            ! call GAapostDriver

        case default

            ! Call error handler
            call gdErrorHandler('GridAdaptor: unknown driver option')
        
        end select

        ! Postprocessing
        !===============
        ! call PostProcessGA

    end subroutine

    subroutine GAinit(grid,options,environment,magneticField)

        ! Description
        !============

        ! Declare variables
        !==================
        ! Arguments
        type(GAGridUDT), intent(inout)          :: grid 
        type(GAoptionsUDT), intent(in)          :: options
        type(EnvironmentUDT), intent(in)        :: environment
        type(MagneticFieldUDT), intent(in)      :: magneticField

        ! Variables
        integer(I8) :: i
        integer(I8), allocatable, dimension(:) :: tv, cvLookUp
        logical :: cells(grid%cell%ntot), is_ordered(grid%cell%ntot), &
            use_nsep, use_sepID, start
        character(:), allocatable :: base_func


        ! Initialize
        !===========
 
        ! Set gaoptions
        call options%Set()

        associate(&
            c  => grid%cell, &
            v  => grid%vert &
            )
        ! Recompute cell centers
        do i = 1, c%ntot
            ! Get cell vertices
            tv = GetCellVertGA(c, i)

            ! Compute coordinates
            call c%x%SetSingleElement(i,sum(v%x%GetMultipleElements(tv))/real(size(tv), kind=R8)) 
            call c%y%SetSingleElement(i,sum(v%y%GetMultipleElements(tv))/real(size(tv), kind=R8)) 

        end do


        ! Check order of vertices  (see GetGeo_usCouples.m)
        cells = .true.
        call grid%CheckVertOrder(is_ordered, cells) 

        if (.not. all(is_ordered)) then

            ! Do ordening
            call grid%ReorderCellConn(is_ordered)

        end if 

        ! Get fsVx from fsFc
        call grid%GetFsVxFromFsFc()

        ! Determine Xpoints and separatrices
        cvLookUp = GetCvLookUp(c)
        use_nsep = .false.
        call grid%GiveXpoints(use_nsep,cvLookUp)
        use_sepID = .false.
        start = .true.
        call grid%GiveSeparatrices(use_nsep, use_sepID, start, cvLookUp)

        ! Identify aligned faces
        call grid%IdentifyAlignedFaces(options,magneticField)

        ! Set up the distance functions
        if (options%dist_function) then

            ! Define base function to compute the distance function
            base_func = options%base_func !'exp(-dist/d)'

            ! Free distance function based on user-defined input
            grid%fun%d_char_type = options%d_char_type
            grid%fun%dist_type   = options%dist_type
            grid%fun%d_rescale   = options%d_rescale
            call grid%fun%ComputeDistanceFunction(grid,options,base_func)

            ! Distance function for high poloidal flux next to the separatrix
            grid%fun_r%d_char_type = options%d_char_type
            grid%fun_r%dist_type   = 'pol_flux_est'
            grid%fun_r%d_rescale   = options%d_rescale
            call grid%fun_r%ComputeDistanceFunction(grid,options,base_func)

            ! Distance function for wall proximity
            grid%fun_wall%d_char_type = options%d_char_type
            grid%fun_wall%dist_type = options%dist_type_wall
            grid%fun_wall%d_rescale = options%d_rescale_wall
            call grid%fun_wall%ComputeDistanceFunction(grid,options,base_func)

        end if 


        ! Check the connectivity
        if (options%debug) then
            call grid%CheckUnstructuredGrid(.false.)
        end if

        ! correct face labels

        ! Remove some connectivity field

        ! Identify farSOL cells

        end associate        


    end subroutine

    subroutine GAInternalDriver(grid,options,environment,magneticField)

        ! Description
        !============
        ! Internal driver for the grid adaptation where all real adaptation take place such as removal of small triangles, stacked triangles

        ! Declare variables
        !==================
        ! Argument
        type(GAGridUDT), intent(inout)       :: grid
        type(GAoptionsUDT), intent(in)       :: options
        type(EnvironmentUDT), intent(in)     :: environment
        type(MagneticFieldUDT), intent(in)   :: magneticField

        ! Calculate quality metric

        ! Remove Small triangles

        ! Remove flux tubes with only two triangles

        ! Stacked to cutcell

        ! Splitting non-alinged quads

        ! Splitting trapezoids in concave shaved-off flux tube

        ! Splitting  and merging

        ! Stacked triangles


        ! Remove sticking out triangles

        ! Remove boundary flux tubes with only two triangles

        ! Remove stickout quad

        ! Boundary layer grid




    end subroutine






end module