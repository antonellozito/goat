!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module serves to define the required 'state' variables to solve 
! the optimization problem for the grid optimization problem. Here, the 
! state consists of the grid, magneticField, and environment types. 
! These types are set in the gdmod_types module.  

! Notes
!======

module gdmod_state
    
    ! Initialize
    !============
    ! Load modules
    use gdmod_types
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

    ! Derived types
    !==============
    type, extends(StateUDT) :: StateGDUDT

        ! Description
        !============
        ! State definition for the grid deformation 

        ! Fields
        type(GridUDT)               :: grid 
        type(MagneticFieldUDT)      :: magneticField
        type(EnvironmentOptionsUDT) :: environment

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