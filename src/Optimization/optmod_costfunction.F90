!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains the generic costfunction class, which contains
! the generic cost function parameters and evaluation routines. 

module optmod_costfunction
    
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
    ! General cost function type
    type, abstract :: CostfunctionUDT

        ! Description
        !============
        ! Defines the basic optimization problem: it has a set of 
        ! design variables, constraints, and a cost function. 

    end type

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