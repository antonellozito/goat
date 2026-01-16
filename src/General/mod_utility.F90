!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains some general utility routines (e.g. wall clock
! timer)

module mod_utility

    ! Initialize
    !============
    ! Load modules
    use mod_precision

    ! The usual
    implicit none
    save
    public 

    contains 

    ! Wall time computation
    subroutine wall_time(t)

        ! Description
        !============
        ! Compute wall clock time by calling the system clock 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(out)   :: t 

        ! Auxiliary
        integer(I16)            :: CPUcount, CPUrate 
        
        ! Compute
        !========
        call system_clock(CPUcount, CPUrate)
        t = real(CPUcount, kind=R8)/real(CPUrate, kind = R8)
        return 

    end 

end module