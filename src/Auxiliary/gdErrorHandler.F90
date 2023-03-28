subroutine gdErrorHandler(msg)
    ! Description
    !============
    ! Wrapper for error handling in the grid deformation module. 
    ! Currently only prints out a message but does not stop the program
    ! execution. 

    ! Initialize
    !===========
    ! The usual
    implicit none

    ! Declare variables
    character(*)     :: msg

    ! Error handling
    !===============
    ! Print
    print *, 'gd error: ', msg
    stop 

end subroutine