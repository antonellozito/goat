!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! Very simple module for error and warning handling. 

module mod_errorhandler

    ! Initialize
    !===========
    use mod_precision

    ! The usual
    implicit none
    save
    private
    public :: ErrorStack, ErrorHandlerUDT, gdErrorHandler 

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!
    
    ! Basic error type
    type :: ErrorMessageUDT 

        ! Description
        !============
        ! Basic type that contains the information of an error. 
        character(:), allocatable       :: msg ! error message

    end type

    ! Error handler type
    type :: ErrorHandlerUDT 
        
        ! Description
        !============
        ! General error handler object. Contains a list of the thrown 
        ! errors. 
        integer(I8)                         :: n = 0 ! current number of errors
        integer(I8)                         :: t = 0 ! current number of trackers
        integer(I8), allocatable            :: ID(:) ! error IDs
        logical, allocatable                :: trackstack(:) ! tracking stack
        type(ErrorMessageUDT), allocatable  :: msgstack(:)
        

    contains 

        ! Add error
        procedure :: Add        => AddError

        ! Start exception tracking
        procedure :: StartTrack

        ! Error checker
        procedure :: ErrorState

        ! End exception tracking
        procedure :: EndTrack 

        ! Check storage
        procedure :: CheckStorage

        ! Printer
        procedure :: Print

        ! Final
        final :: FinalizeErrorHandler

    end type

    ! Define error handlers
    !======================
    type(ErrorHandlerUDT) :: ErrorStack

    contains 

    !==================================================================!
    !                                                                  !
    !                             ROUTINES                             !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                              GENERAL                             !
    !------------------------------------------------------------------!

    subroutine gdErrorHandler(msg, IDin, severityin)

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
        character(*)                        :: msg
        integer, intent(in), optional       :: IDin, severityin

        integer                             :: ID, severity

        ! Check inputs
        !=============
        if (present(IDin)) then 
            ID = IDin 
        else
            ID = 0 ! default: no ID
        end if 
        if (present(severityin)) then 
            severity = severityin 
        else 
            severity = 1 ! default: stop program
        end if 
    
        ! Error handling
        !===============
        ! Add to stack
        call ErrorStack%Add(msg, ID)

        ! Print
        print *, 'gd error: ', msg

        ! Check for stop
        if (severity > 0) then 
            stop 
        end if 
    
    end subroutine

    !------------------------------------------------------------------!
    !                           ERROR HANDLER                          !
    !------------------------------------------------------------------!

    ! Add error
    subroutine AddError(errorhandler, msg, ID)

        ! Description
        !============
        ! Add error to the error handler

        ! Declare variables
        !==================
        ! Arguments
        class(ErrorHandlerUDT)      :: errorhandler 
        character(*), intent(in)    :: msg 
        integer, intent(in)         :: ID 

        ! Add
        !====
        ! Associate
        associate(n         => errorhandler%n,      &
            t               => errorhandler%t)

        ! Update warning counter
        n = n + 1

        ! Check field sizes etc
        call errorhandler%CheckStorage()

        ! Add the error
        errorhandler%msgstack(n)%msg    = msg 
        errorhandler%ID(n)              = ID 

        ! Check if we're tracking
        if (t > 0) then 
            ! Error encountered, set to true
            errorhandler%trackstack(t) = .true. 
        end if 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Check storage
    subroutine CheckStorage(errorhandler)

        ! Description
        !============
        ! Check if sufficient storage is available to store any new 
        ! messages.

        ! Declare variables
        !==================
        ! Arguments
        class(ErrorHandlerUDT)              :: errorhandler 

        ! Auxiliary
        type(ErrorMessageUDT), allocatable  :: tempmsgstack(:)
        integer, allocatable                :: tempID(:)
        logical, allocatable                :: temptrackstack(:)

        ! Check allocation
        !=================
        if (.not. allocated(errorhandler%ID)) then 
            ! Allocate to modest initial storage
            allocate(errorhandler%ID(10))
        end if 
        if (.not. allocated(errorhandler%msgstack)) then 
            allocate(errorhandler%msgstack(10))
        end if 
        if (.not. allocated(errorhandler%trackstack)) then 
            allocate(errorhandler%trackstack(10))
        end if 

        ! Check size
        !===========
        if (size(errorhandler%ID) <= errorhandler%n) then 
            ! Extend 
            tempID = errorhandler%ID 
            errorhandler%ID = [tempID, tempID]
        end if 
        if (size(errorhandler%msgstack) <= errorhandler%n) then 
            tempmsgstack = errorhandler%msgstack
            deallocate(errorhandler%msgstack)
            allocate(errorhandler%msgstack(2*errorhandler%n))
            errorhandler%msgstack(1:size(tempmsgstack)) = tempmsgstack
        end if 
        if (size(errorhandler%trackstack) <= errorhandler%t) then 
            temptrackstack = errorhandler%trackstack 
            errorhandler%trackstack = [temptrackstack, temptrackstack .and. (.not. temptrackstack)]
        end if 



    end subroutine

    ! Start exception tracking
    subroutine StartTrack(errorhandler)

        ! Description
        !============
        ! Start tracking whether an exception is executed until the 
        ! next EndTrack function is called 

        ! Declare variables
        !==================
        ! Arguments
        class(ErrorHandlerUDT)      :: errorhandler 

        ! Initialize
        !===========
        ! Associate
        associate(t         => errorhandler%t)

        ! Update the tracking counter
        t = t + 1

        ! Check storage
        call errorhandler%CheckStorage()

        ! Add the tracker
        errorhandler%trackstack(t) = .false. ! no error upon entry

        ! Housekeeping
        !=============
        end associate



    end subroutine

    ! Check if error has occured during current tracking
    function ErrorState(errorhandler) result(out)

        ! Description
        !============
        ! Simple function to check if an error has occured while 
        ! tracking errors with StartTrack. If we are tracking, the result
        ! is >= 0, with > 0 indicating an active error state and == 0 an
        ! inactive error state. If we are not tracking, the result is 
        ! < 0. 

        ! Declare variables
        !==================
        ! Arguments
        class(ErrorHandlerUDT)      :: errorhandler 
        integer                     :: out 

        ! Initialize
        !===========
        out = -1 

        ! Check
        !======
        ! Are we tracking?
        if (errorhandler%t <= 0) then 
            ! Not tracking, return
            return 
        end if 

        ! Check the latest track
        if (errorhandler%trackstack(errorhandler%t)) then 
            ! Error is active
            out = 1
        else
            out = 0
        end if 

    end function

    ! End exception tracking
    subroutine EndTrack(errorhandler)

        ! Description
        !============
        ! End tracking whether an exception is executed 

        ! Declare variables
        !==================
        ! Arguments
        class(ErrorHandlerUDT)      :: errorhandler 

        ! Initialize
        !===========
        ! Associate
        associate(t         => errorhandler%t)

        ! Remove the tracker
        errorhandler%trackstack(t) = .false. ! no error upon exit

        ! Update the tracking counter
        t = t - 1

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Printer
    subroutine Print(errorhandler)

        ! Description
        !============
        ! Print the current errors

        ! Declare variables
        !==================
        ! Arguments
        class(ErrorHandlerUDT)       :: errorhandler 

        ! Auxiliary
        
        ! Loop
        integer(I8)                 :: i 
        
        ! Print messages
        !===============
        if (errorhandler%n == 0) then 
            print *, 'Errorhandler: no recorded errors'
        else
            print *, 'Errorhandler: ', errorhandler%n, ' recorded errors.'
            print *, 'Errors (ID, message):'
            do i = 1, errorhandler%n 
                print *, 'ID: ', errorhandler%ID(i), ', message: ', errorhandler%msgstack(i)%msg 
            end do
        end if 
        if (errorhandler%t == 0) then 
            print *, 'Errorhandler: not tracking errors'
        else
            print *, 'Errorhandler: currently ', errorhandler%t, ' trackers active'
            print *, 'Errorhandler: currently ', count(errorhandler%trackstack(1:errorhandler%t)), ' triggered trackers'
        end if 

    end subroutine

    ! Finalizer
    subroutine FinalizeErrorHandler(errorhandler)

        ! Description
        !============
        ! Finalization routine of the error handler, which should be
        ! called when the object is destroyed (so normally only at
        ! exit of the program). Here, we print out the number of errors
        ! etc 

        ! Declare variables
        !==================
        ! Arguments
        type(ErrorHandlerUDT)       :: errorhandler 

        ! Auxiliary
        
        ! Loop
        integer(I8)                 :: i 
        
        ! Print messages
        !===============
        if (errorhandler%n == 0) then 
            print *, 'Errorhandler exited with no recorded errors'
        else
            print *, 'Errorhandler exited with ', errorhandler%n, ' recorded errors.'
            print *, 'Errors (ID, message):'
            do i = 1, errorhandler%n 
                print *, 'ID: ', errorhandler%ID(i), ', message: ', errorhandler%msgstack(i)%msg 
            end do
        end if 
        if (errorhandler%t == 0) then 
            print *, 'Errorhandler: exited normal with no active trackers anymore'
        else
            print *, 'WARNING: errorhandler is still tracking, this indicates a missing EndTrack somewhere in the code!'
            print *, 'Errorhandler: currently ', errorhandler%t, ' trackers active'
            print *, 'Errorhandler: currently ', count(errorhandler%trackstack(1:errorhandler%t)), ' triggered trackers'
        end if 

    end subroutine 

end module