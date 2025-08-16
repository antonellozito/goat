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
    use gamod_userinput

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
        type(GridUDT)               :: grid
        type(EnvironmentUDT)        :: environment
        type(MagneticFieldUDT)      :: magneticField
        type(GoatoptionsUDT)        :: options        
        
        type(GAoptionsUDT)          :: gaoptions

        ! Initialize grid adaptation
        !===========================
        ! TODO call GAinit
        call gaoptions%Set()    ! Also do other options???

        ! Driver Selection
        !=================
        select case (gaoptions%meth)

        case ('simple')

            ! Regular grid adaption
            ! call GAInternalDriver

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






end module