!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! Module that contains (relatively) general purpose read-write routines.

module mod_readwrite

    ! Initialize
    !============
    ! Modules
    use mod_precision

    ! The usual
    implicit none
    save
    public 

    contains 

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    ! Read in single line
    subroutine ReadSingleLine(fid, thisline, reachedeof)

        ! Description
        !============
        ! This routine reads the next line of the file with unit 'fid', 
        ! which should already be opened. This means that the next line 
        ! is read in. The line is returned as a character array in 
        ! 'line' (which should be parsed as a variable length character
        ! array, which is not yet allocated). The routine will first 
        ! compute how long the line is, then allocate the character 
        ! array, and read in the line contents. 

        ! Declare modules
        !================
        use iso_fortran_env, only : iostat_eor

        ! Declare variables
        !==================
        ! Arguments
        integer, intent(in)             :: fid 
        character(len=:), allocatable   :: thisline
        logical, intent(out)            :: reachedeof

        ! Auxiliary
        integer                     :: readstatus, chunksize, totalsize
        integer, parameter          :: buffersize = 10
        character(buffersize)       :: buffer

        ! Initialize
        !===========
        ! Initialize total size of line
        totalsize = 0
        reachedeof = .false.

        ! Check if line is already allocated
        if (allocated(thisline)) then
            deallocate(thisline)
        end if

        ! Read to get size of record
        do while (.true.) 

            ! Read next chunk
            read (fid, '(A)', advance='no', size=chunksize, iostat=readstatus) buffer 

            ! Add length to total size
            totalsize = totalsize + chunksize

            ! Check EOF/EOR
            reachedeof = is_iostat_end(readstatus)
            if (reachedeof .or. (readstatus .eq. iostat_eor)) then 
                ! Exit the loop
                exit 
            end if

        end do 

        ! Allocate the line
        allocate(character(totalsize) :: thisline)

        ! Go back one line 
        backspace(fid) 

        ! Read
        read (fid, '(A)',  iostat=readstatus) thisline 
        print *, thisline 

    end subroutine

end module