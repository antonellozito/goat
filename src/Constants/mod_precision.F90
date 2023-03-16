!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! General overarching module for precision definition. 

module mod_precision

    ! Initialize
    !============
    ! The usual
    implicit none
    save
    public 

    ! Precision types
    integer, parameter       :: R4      = selected_real_kind(6)
    integer, parameter       :: R8      = selected_real_kind(14)
    integer, parameter       :: R16     = selected_real_kind(33)

    integer, parameter       :: C32     = 32
    integer, parameter       :: C64     = 64
    integer, parameter       :: C128    = 128

    integer, parameter       :: I4      = selected_int_kind(4)
    integer, parameter       :: I8      = selected_int_kind(8)
    integer, parameter       :: I16     = selected_int_kind(16)

end module