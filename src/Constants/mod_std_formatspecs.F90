! Description
!============
! General overarching module for definition of standard format 
! specifications. Should make it easier to consistently read/write
! output. 

module mod_std_formatspecs

    ! Initialize
    !============
    ! The usual
    implicit none
    save
    public 

    character(*), parameter :: spacefm      = '1X'
    character(*), parameter :: Ifm          = 'I11'
    character(*), parameter :: Rfm          = 'F15.8'


end module