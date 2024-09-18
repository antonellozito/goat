!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module provides some basic utilities when interfacing C in 
! Fortran


module CUtilities

    ! Modules
    !========
    use, intrinsic :: iso_c_binding

    implicit none 
    private 
    public :: FreeCptr

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    !==================================================================!
    !                                                                  !
    !                             INTERFACES                           !
    !                                                                  !
    !==================================================================!
    
    interface 

        ! C pointer memory deallocation 
        subroutine FreeCptr(cptr) &
            bind(c, name='free')
            use, intrinsic :: iso_c_binding, only: c_ptr
            type(c_ptr) :: cptr 
        end subroutine 

    end interface 

end module
