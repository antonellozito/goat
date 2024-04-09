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
    use, intrinsic :: ieee_arithmetic

    implicit none
    public 

    ! Pi
    real(R8), parameter :: myone_R8 = 1
    real(R8), parameter :: pi_R8 = 4*atan(myone_R8)

    contains 

    ! Function that returns NaN
    function nanval_R8() result(c)

        ! Description
        !============
        ! Returns NaN
        real(R8)   :: c
        
        c = ieee_value(c, ieee_quiet_nan)


    end function

    ! Function that returns positive infinity
    function posinfval_R8() result(c)

        ! Description
        !============
        ! Returns NaN
        real(R8)   :: c
        
        c = ieee_value(c, ieee_positive_inf)


    end function

end module