!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! General overarching module for definition of characters with special
! meaning. Only to be used for string manipulation. 

module mod_specialchars

    ! Initialize
    !============
    ! The usual
    implicit none
    save
    public 

    character, parameter   :: commentchar  = '#'
    character, parameter   :: delimiter    = "'"
    character, parameter   :: matstart     = '['
    character, parameter   :: matend       = ']'
    character, parameter   :: decpoint     = '.'
    character, parameter   :: rowsep       = ','
    character, parameter   :: rowdel       = ';'
    character, parameter   :: veccon       = ':'
    character, parameter   :: repeatchar   = '*'

end module