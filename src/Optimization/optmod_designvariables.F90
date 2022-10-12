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
    use mod_precision

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
    type :: DesignVariablesUDT

        ! Description
        !============
        integer(I8)                         :: nphi
        real(R8), allocatable               :: phi(:)
        
    contains

        ! Housekeeping procedures
        procedure :: Allocate           => AllocateDesign
        procedure :: Deallocate         => DeallocateDesign

    end type

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    ! Abstract interfaces
    !====================

    contains

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    ! Allocation
    subroutine AllocateDesign(designvariables)

        ! Description
        !============
        ! Allocate the design variables

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariablesUDT)       :: designvariables

        ! Allocate
        !=========
        allocate(designvariables%phi(designvariables%nphi))

    end subroutine

    ! Deallocation
    subroutine DeallocateDesign(designvariables)

        ! Description
        !============
        ! Allocate the design variables

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariablesUDT)       :: designvariables

        ! Allocate
        !=========
        deallocate(designvariables%phi)

    end subroutine

end module