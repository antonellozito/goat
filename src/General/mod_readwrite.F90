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
    use mod_errorhandler
    use mod_specialchars 

    ! The usual
    implicit none
    save
    public 



    !==================================================================!
    !                                                                  !
    !                          INTERFACES                              !
    !                                                                  !
    !==================================================================!

    ! Read array
    interface ReadArray 
        module procedure ReadArrayI8, ReadArrayR8, ReadNamedArrayI8, &
            ReadNamedArrayR8, ReadNamedArray2DI8, ReadNamedArray2DR8, &
            ReadNamedArray3DR8, ReadArray3DR8
    end interface

    ! Write array
    interface WriteArray
        module procedure WriteArrayI8, WriteArrayR8
    end interface

    contains 

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                       Basic reading                              !
    !------------------------------------------------------------------!

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

    end subroutine

    ! Read until substring is found in line (position is next line)
    subroutine ReadUntilFound(fid, substring, reachedeof)

        ! Description
        !============
        ! This routine keeps reading next lines until either the 
        ! a line contains the desired substring or until the end of file
        ! has been reached. Note that the current line being read in
        ! may depend on previous read statements and is not reset here!

        ! Declare variables
        !==================
        ! Arguments
        integer, intent(in)             :: fid 
        character(*), intent(in)        :: substring
        logical, intent(out)            :: reachedeof

        ! Auxiliary
        character(:), allocatable       :: thisline
        
        ! Loop
        logical                         :: isnotfound 

        ! Initialize
        !===========
        isnotfound = .true. 
        reachedeof = .false. 

        ! Loop
        !=====
        do while (isnotfound .and. (.not. reachedeof))

            ! Read next line
            call ReadSingleLine(fid, thisline, reachedeof)

            ! Check 
            if (reachedeof) then 
                ! Will exit, do nothing
            end if 
            if (index(thisline, substring) /= 0) then 
                ! Substring found, exit
                isnotfound = .false. 
            end if

        end do


    end subroutine

    ! Read array of reals with predefined size
    subroutine ReadArrayR8(fid, n, out)

        ! Description
        !============
        ! Read an array with predefined size n (array goes from 1 to n)
        ! from a file that has already been opened and is located at 
        ! the correct position to commence reading. 

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)         :: n, fid 
        real(R8), intent(inout)         :: out(1:n)

        ! Loop
        integer(I8)                     :: i 

        ! Read
        read(fid, *) (out(i), i = 1, n)

    end subroutine 

    subroutine ReadArray2DR8(fid, n, out)

        ! Description
        !============
        ! Read a 2D array from a file. 

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)         :: n, fid 
        real(R8), intent(inout)         :: out(:, :)

        ! Auxiliary
        real(R8), allocatable, dimension(:)     :: temp 

        ! Loop
        integer(I8)                     :: i 

        ! Read
        !=====
        ! First read in 1D array
        allocate(temp(size(out, 1)*size(out, 2)))
        temp = 0
        read(fid, *) (temp(i), i = 1, n)

        ! Reshape
        out = reshape(temp, [size(out, 1), size(out, 2)])

    end subroutine

    subroutine ReadArray3DR8(fid, n, out)

        ! Description
        !============
        ! Read a 2D array from a file. 

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)         :: n, fid 
        real(R8), intent(inout)         :: out(:, :, :)

        ! Auxiliary
        real(R8), allocatable, dimension(:)     :: temp 

        ! Loop
        integer(I8)                     :: i 

        ! Read
        !=====
        ! First read in 1D array
        allocate(temp(size(out, 1)*size(out, 2)))
        temp = 0
        read(fid, *) (temp, i = 1, n)

        ! Reshape
        out = reshape(temp, [size(out, 1), size(out, 2), size(out, 3)])

    end subroutine

    ! Read array of reals with header and predefined size
    subroutine ReadNamedArrayR8(fid, n, out, name)

        ! Description
        !============
        ! Read an array that has a header line that starts with 'name'.
        ! We read from the current position in fid.

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)         :: n, fid 
        real(R8), intent(inout)         :: out(1:n)
        character(*), intent(in)        :: name 

        ! Auxiliary
        logical                         :: reachedeof 

        ! Read
        !=====
        ! Find header
        call ReadUntilFound(fid, name, reachedeof)
        if (reachedeof) then 
            call gdErrorHandler('ReadNamedArrayR8: could not find array ' // &
                'with (part of) header: "' // name // '" in the current file ' // &
                'starting from the given position')
        end if 

        ! Read
        call ReadArrayR8(fid, n, out)

    end subroutine 

    ! Read 2D array of reals with header and predefined size
    subroutine ReadNamedArray2DR8(fid, n, out, name)

        ! Description
        !============
        ! Read an array that has a header line that starts with 'name'.
        ! We read from the current position in fid.

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)         :: n, fid 
        real(R8), intent(inout)         :: out(:, :)
        character(*), intent(in)        :: name 

        ! Auxiliary
        logical                         :: reachedeof 

        ! Read
        !=====
        ! Find header
        call ReadUntilFound(fid, name, reachedeof)
        if (reachedeof) then 
            call gdErrorHandler('ReadNamedArray2DR8: could not find array ' // &
                'with (part of) header: "' // name //'" in the current file ' // &
                'starting from the given position')
        end if 

        ! Read
        call ReadArray2DR8(fid, n, out)

    end subroutine

    ! Read 3D array of reals with header and predefined size
    subroutine ReadNamedArray3DR8(fid, n, out, name)

        ! Description
        !============
        ! Read an array that has a header line that starts with 'name'.
        ! We read from the current position in fid.

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)         :: n, fid 
        real(R8), intent(inout)         :: out(:, :, :)
        character(*), intent(in)        :: name 

        ! Auxiliary
        logical                         :: reachedeof 

        ! Read
        !=====
        ! Find header
        call ReadUntilFound(fid, name, reachedeof)
        if (reachedeof) then 
            call gdErrorHandler('ReadNamedArray3DR8: could not find array ' // &
                'with (part of) header: "' // name //'" in the current file ' // &
                'starting from the given position')
        end if 

        ! Read
        call ReadArray3DR8(fid, n, out)

    end subroutine

    ! Read array of integers with predefined size
    subroutine ReadArrayI8(fid, n, out)

        ! Description
        !============
        ! Read an array with predefined size n (array goes from 1 to n)
        ! from a file that has already been opened and is located at 
        ! the correct position to commence reading. 

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)         :: n, fid 
        integer(I8), intent(inout)      :: out(1:n)

        ! Loop
        integer(I8)                     :: i 

        ! Read
        read(fid, *) (out(i), i = 1, n)

    end subroutine 

    subroutine ReadArray2DI8(fid, n, out)

        ! Description
        !============
        ! Read a 2D array from a file. 

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)         :: n, fid 
        integer(I8), intent(inout)      :: out(:, :)

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: temp 

        ! Loop
        integer(I8)                     :: i 

        ! Read
        !=====
        ! First read in 1D array
        allocate(temp(size(out, 1)*size(out, 2)))
        temp = 0
        read(fid, *) (temp(i), i = 1, n)

        ! Reshape
        out = reshape(temp, [size(out, 1), size(out, 2)])

    end subroutine

    ! Read array of integers with header and predefined size
    subroutine ReadNamedArrayI8(fid, n, out, name)

        ! Description
        !============
        ! Read an array that has a header line that starts with 'name'.
        ! We read from the current position in fid.

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)         :: n, fid 
        integer(I8), intent(inout)      :: out(1:n)
        character(*), intent(in)        :: name 

        ! Auxiliary
        logical                         :: reachedeof 

        ! Read
        !=====
        ! Find header
        call ReadUntilFound(fid, name, reachedeof)
        if (reachedeof) then 
            call gdErrorHandler('ReadNamedArrayI8: could not find array ' // &
                'with (part of) header: "' // name //'" in the current file ' // &
                'starting from the given position')
        end if 

        ! Read
        call ReadArrayI8(fid, n, out)

    end subroutine 

    ! Read 2D array of integers with header and predefined size
    subroutine ReadNamedArray2DI8(fid, n, out, name)

        ! Description
        !============
        ! Read an array that has a header line that starts with 'name'.
        ! We read from the current position in fid.

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)         :: n, fid 
        integer(I8), intent(inout)      :: out(:, :)
        character(*), intent(in)        :: name 

        ! Auxiliary
        logical                         :: reachedeof 

        ! Read
        !=====
        ! Find header
        call ReadUntilFound(fid, name, reachedeof)
        if (reachedeof) then 
            call gdErrorHandler('ReadNamedArray2DI8: could not find array ' // &
                'with (part of) header: "' // name //'" in the current file ' // &
                'starting from the given position')
        end if 

        ! Read
        call ReadArray2DI8(fid, n, out)

    end subroutine

    !------------------------------------------------------------------!
    !                           Writing                                !
    !------------------------------------------------------------------!
    
    subroutine WriteArrayI8(a, filename)

        ! Description
        !============
        ! Write an integer array in a file

        ! Declare variables
        !==================
        ! Modules 
        use mod_specialchars, only : filesepchar

        ! Arguments
        integer(I8), intent(in)  :: a(:)
        character(*), intent(in) :: filename 

        ! Auxiliary
        integer :: fu     
        integer(I8) :: i
        character(:), allocatable :: dir

        ! Construct writing directory
        dir = './output' // filesepchar // filename // '.dat'

        ! Open file
        open (action='write', file=trim(dir), newunit=fu, &
             status='unknown')

        ! Size data
        write (fu, *) 'Elements'
        write (fu, *) size(a)

        ! Array
        write (fu, *) 'ID val(ID)'
        do i = 1, size(a)
            write(fu, *) i, a(i)
        end do

        close(fu)

    end subroutine

    subroutine WriteArrayR8(a, filename)

        ! Description
        !============
        ! Write an integer array in a file

        ! Declare variables
        !==================
        ! Modules  
        use mod_specialchars, only : filesepchar

        ! Arguments
        real(R8), allocatable :: a(:)
        character(*), intent(in) :: filename 

        ! Auxiliary
        integer :: fu     
        integer(I8) :: i
        character(:), allocatable :: dir

        ! Construct writing directory
        dir = './output' // filesepchar // filename // '.dat'

        ! Open file
        open (action='write', file=trim(dir), newunit=fu, &
             status='unknown')

        ! Size data
        write (fu, *) 'Elements'
        write (fu, *) size(a)

        ! Array
        write (fu, *) 'ID val(ID)'
        do i = 1, size(a)
            write(fu, *) i, a(i)
        end do

        close(fu)

    end subroutine

end module