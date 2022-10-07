!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains the generic 'state' class, which is a sort of 
! wrapper for all the possible state arguments that may be needed when
! evaluating the cost function, updating the design, ... 

! Usage
!======
! In order to use this module, the user should define his own state 
! module, where the type definition inherits from this type. For 
! example:
!
! type, extends(StateUDT) :: MyStateUDT
!   real(kind=8), allocatable :: myfield
! end type

module optmod_state
    
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

    ! Base types
    type, abstract :: StateUDT

        ! Description
        !============
        ! Defines the basic optimization problem: it has a set of 
        ! design variables, constraints, and a cost function. 

    end type



    ! Abstract types
    !===============

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    contains

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

end module