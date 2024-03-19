!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! General overarching module for definitions of constants

module mod_constants

    use mod_precision

    implicit none
    public 

    ! Pi
    real(R8), parameter :: myone_R8 = 1
    real(R8), parameter :: pi = 4*atan(myone_R8)

end module