!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains all the necessary routines to manipulate the 
! design variables (except for the actual optimization routines). The 
! only implementation provided here is an abstract type and abstract 
! interfaces to the generic procedures that are called by the optimizer.
! The actual implementation of these routines should be made in a 
! user-defined module, where the main object is inherited. For example:
!
! type, extends(DesignVariablesUDT) :: DesignVariablesCoordinatesUDT
! ! <define some additional fields if desired here>
! contains
! ! <point towards the correct routines here>
! InitializeDesign => InitializeDesignCoordinates
! ... < other generic functions come here >
! end type
! 
! The required arguments etc can be found in the 'interface' section, 
! where the generic interface of each of these procedures is defined. 

! Notes
!======
! Note 1: this module is specifically built for grid optimization. Most
! routines therefore expect at least the input arguments 'grid', 
! 'magneticField', and 'environment'. 

! Note 2: although the code structure (should) allow(s) for any type of 
! design variable, currently only the 'coordinates' type is implemented.

module optmod_designvariables
    
    ! Initialize
    !============
    ! Load modules
    use optmod_state

    ! The usual
    implicit none
    save
    public 

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    ! Abstract types
    !===============
    type, abstract :: DesignVariablesUDT

        ! Description
        !============
        integer(kind=8)                         :: nphi
        real(kind=8), allocatable               :: phi(:)
        
    contains

        procedure(InitializeDesign), deferred :: InitializeDesign

    end type

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    ! Abstract interfaces
    !====================
    abstract interface 

        ! Design initialization
        subroutine InitializeDesign(designvariables, state)
            ! This routine should initialize the design, based on the 
            ! data given in the state structure. The user has to provide
            ! the initialization (with possible choices between 
            ! different design variables) and, if desired, unpacking 
            ! routines for the state parameters. 
            import:: DesignVariablesUDT, StateUDT
            class(DesignVariablesUDT)   :: designvariables 
            class(StateUDT)             :: state

        end subroutine

    end interface 

    contains

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

end module