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
        
        type(GAoptionsUDT)          :: options

        ! Initialize grid adaptation
        !===========================
        ! TODO call GAinit
        

        ! Driver Selection
        !=================
        select case (options%meth)

        case ('simple')

            ! Regular grid adaption
            ! call GAInternalDriver(grid,gaoptions,environment,magneticField)

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
        class(GridUDT)          :: grid 
        class(GAoptionsUDT)     :: options
        class(EnvironmentUDT)   :: environment

        ! Variables
        integer(I8) :: icv, s, nv, vc(1:100)
        real(R8) :: cvX(grid%cell%ntot), cvY(grid%cell%ntot)


        ! Set gaoptions
        call options%Set()    ! Also do other options, or make then a subtype???


        ! Recompute cell centers
        do icv = 1, grid%cell%ntot
            s = grid%cell%vertP(icv,1)
            nv = grid%cell%vertP(icv,2)
            vc(1:nv) = grid%cell%vert(s:s+nv-1)
            cvX(icv) = sum(grid%vert%x(vc(1:nv)))/nv
            cvY(icv) = sum(grid%vert%y(vc(1:nv)))/nv
        end do

        grid%cell%x = cvX
        grid%cell%y = cvY

        ! Check order of vertices  (see GetGeo_usCouples.m)
        ! call CheckVertOrder, 

        ! Get fsVx from fsFc

        ! Determine Xpoints and separatrices

        ! Identify aligned faces


        ! Check the connectivity

        ! correct face labels

        ! Remove some connectivity field

        ! Identify farSOL cells
        


    end subroutine

    subroutine GAInternalDriver()

        ! Description
        !============
    end subroutine






end module