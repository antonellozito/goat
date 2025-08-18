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
        type(GridUDT), intent(inout)                :: grid
        type(EnvironmentUDT), intent(in)            :: environment
        type(MagneticFieldUDT), intent(in)          :: magneticField     
        
        type(GAoptionsUDT), intent(in)              :: options

        ! Initialize grid adaptation
        !===========================
        call GAinit(grid,options,environment)
        

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
        ! call postprocessGA

    end subroutine

    subroutine GAinit(grid,options,environment)

        ! Description
        !============

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(inout)            :: grid 
        type(GAoptionsUDT), intent(in)          :: options
        type(EnvironmentUDT), intent(in)        :: environment

        ! Variables
        integer(I8) :: i
        integer(I8), allocatable, dimension(:) :: tv
        logical :: cells(grid%cell%ntot), is_ordered(grid%cell%ntot)


        ! Initialize
        !===========
 
        ! Set gaoptions
        call options%Set()    ! Also do other options, or make then a subtype???

        associate(&
            c  => grid%cell, &
            v  => grid%vert &
            )
        ! Recompute cell centers
        do i = 1, c%ntot
            ! Get cell vertices
            tv = GetCellVert(c, i)

            ! Compute coordinates
            c%x(i) = sum(v%x(tv))/real(size(tv), kind=R8)
            c%y(i) = sum(v%y(tv))/real(size(tv), kind=R8)
        end do
        end associate

        ! Check order of vertices  (see GetGeo_usCouples.m)
        cells = .true.
        call CheckVertOrder(grid, is_ordered, cells) 

        if (.not. all(is_ordered)) then
            ! Do ordening
            call ReorderCellConn(grid, is_ordered)
        end if 

        ! Get fsVx from fsFc

        ! Determine Xpoints and separatrices

        ! Identify aligned faces


        ! Check the connectivity

        ! correct face labels

        ! Remove some connectivity field

        ! Identify farSOL cells
        


    end subroutine

    subroutine GAInternalDriver(grid,options,environment,magneticField)

        ! Description
        !============
        ! Internal driver for the grid adaptation where all real adaptation take place such as removal of small triangles, stacked triangles

        ! Declare variables
        !==================
        ! Argument
        type(GridUDT), intent(inout)         :: grid
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