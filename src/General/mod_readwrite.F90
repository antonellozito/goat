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
    use mod_specialchars

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

    !------------------------------------------------------------------!
    !                       Basic reading & writing                    !
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

    !------------------------------------------------------------------!
    !                         String manipulation                      !
    !------------------------------------------------------------------!

    ! Scalar integer extraction
    subroutine ExtractIntegerFromString0D(stringval, intval, islegal)

        ! Description
        !============
        ! This routine extracts an integer from a string. Some things 
        ! are checked first: 
        ! - square brackets must be removed first before extraction, if
        !   present
        ! - It must be checked whether the integer is actually an 
        !   integer 
        !
        ! The logical islegal that is returned indicates whether a 
        ! correct value could be extracted (false if this is not the 
        ! case).

        ! Declare variables
        !==================
        ! Arguments
        character(:), allocatable, intent(in)   :: stringval
        integer(I8), intent(out)                :: intval 

        ! Auxiliary
        integer(I8)                             :: stringlength, &
            readstatus
        logical                                 :: islegal ! My lord, is that legal?
        character(:), allocatable               :: tempc 

        logical, allocatable                    :: check(:), templ(:)

        ! Loop
        integer(I8)                             :: i, k

        ! Initialize
        !===========
        ! Length of the string
        stringlength = len(stringval) 

        ! Allocate
        allocate(check(stringlength), templ(stringlength))
        allocate(character(stringlength) :: tempc)

        ! Initialize
        islegal     = .false. 
        check(:)    = .false. 

        ! Checks
        !=======
        ! Is there a decimal point? 
        call CompareStringWithCharacter(stringval, decpoint, templ)
        if (any(templ)) then 
            return 
        end if

        ! Is there a row separator? 
        call CompareStringWithCharacter(stringval, rowsep, templ) 
        if (any(templ)) then 
            return 
        end if

        ! Is there a column separator? 
        call CompareStringWithCharacter(stringval, rowdel, templ) 
        if (any(templ)) then 
            return 
        end if

        ! Remove any square brackets if present
        call CompareStringWithCharacter(stringval, matstart, templ)
        check = check .or. templ 
        call CompareStringWithCharacter(stringval, matend, templ)
        check = check .or. templ 
        check = .not. check 

        ! Trim the string
        do i = 1, stringlength 
            if (check(i)) then 
                ! Add
                tempc(i:i) = stringval(i:i)
            else
                ! Replace by whitespace
                tempc(i:i) = ' ' 
            end if 
        end do 

        ! Read
        !=====
        ! Read
        read (tempc, *, iostat=readstatus) intval

        ! Check if read succeeded
        if (readstatus == 0) then 
            islegal = .true. ! I'll make it legal
        end if

    end subroutine

    ! Array integer extraction
    subroutine ExtractIntegerFromString1D(stringval, intval, islegal)

        ! Description
        !============
        ! This routine extracts an integer array from a string. 
        ! Some things are checked first: 
        ! - square brackets must be removed first before extraction, if
        !   present
        ! - It must be checked whether the integer is actually an 
        !   integer 
        ! - The amount of comma's are counted as this indicates the 
        !   number of elements of the vector.
        !
        ! The logical islegal that is returned indicates whether a 
        ! correct value could be extracted (false if this is not the 
        ! case).

        ! Declare variables
        !==================
        ! Arguments
        character(:), allocatable, intent(in)   :: stringval
        integer(I8), allocatable, intent(out)   :: intval(:) 

        ! Auxiliary
        integer(I8)                             :: stringlength, &
            readstatus, nval
        logical                                 :: islegal ! My lord, is that legal?
        character(:), allocatable               :: tempc 

        logical, allocatable                    :: check(:), templ(:)

        ! Loop
        integer(I8)                             :: i, k
        

        ! Initialize
        !===========
        ! Length of the string
        stringlength = len(stringval) 

        ! Allocate
        allocate(check(stringlength), templ(stringlength))
        allocate(character(stringlength) :: tempc)

        ! Initialize
        islegal     = .false. 
        check(:)    = .false. 

        ! Checks
        !=======
        ! Is there a decimal point? 
        call CompareStringWithCharacter(stringval, decpoint, templ)
        if (any(templ)) then 
            return 
        end if

        ! Is there a column separator? 
        call CompareStringWithCharacter(stringval, rowdel, templ) 
        if (any(templ)) then 
            return 
        end if

        ! Remove any square brackets if present
        call CompareStringWithCharacter(stringval, matstart, templ)
        check = check .or. templ 
        call CompareStringWithCharacter(stringval, matend, templ)
        check = check .or. templ 
        
        ! Count the amount of row delimiters to obtain size
        call CompareStringWithCharacter(stringval, rowsep, templ)
        nval = count(templ)+1

        ! Eliminate row delimiters 
        check = check .or. templ

        ! Trim the string
        check = .not. check 
        do i = 1, stringlength 
            if (check(i)) then 
                ! Add
                tempc(i:i) = stringval(i:i)
            else
                ! Replace by whitespace
                tempc(i:i) = ' ' 
            end if 
        end do 

        ! Read
        !=====
        ! Allocate
        allocate(intval(nval))

        ! Read
        read (tempc, *, iostat=readstatus) intval

        ! Check if read succeeded
        if (readstatus == 0) then 
            islegal = .true. ! I'll make it legal
        end if

    end subroutine

    ! Matrix integer extraction
    subroutine ExtractIntegerFromString2D(stringval, intval, islegal)

        ! Description
        !============
        ! This routine extracts an integer matrix from a string. 
        ! Some things are checked first: 
        ! - square brackets must be removed first before extraction, if
        !   present
        ! - It must be checked whether the integer is actually an 
        !   integer 
        ! - The amount of row and column separators are counted as this 
        !   yields the number of rows and columns
        !
        ! The logical islegal that is returned indicates whether a 
        ! correct value could be extracted (false if this is not the 
        ! case).

        ! Declare variables
        !==================
        ! Arguments
        character(:), allocatable, intent(in)   :: stringval
        integer(I8), allocatable, intent(out)   :: intval(:, :) 

        ! Auxiliary
        integer(I8)                             :: stringlength, &
            readstatus, nval, nrow, ncol, remdivcheck
        logical                                 :: islegal ! My lord, is that legal?
        character(:), allocatable               :: tempc 

        logical, allocatable                    :: check(:), templ(:)

        integer, allocatable                    :: tempi(:)

        ! Loop
        integer(I8)                             :: i, j, k
        

        ! Initialize
        !===========
        ! Length of the string
        stringlength = len(stringval) 

        ! Allocate
        allocate(check(stringlength), templ(stringlength))
        allocate(character(stringlength) :: tempc)

        ! Initialize
        islegal     = .false. 
        check(:)    = .false. 

        ! Checks
        !=======
        ! Is there a decimal point? 
        call CompareStringWithCharacter(stringval, decpoint, templ)
        if (any(templ)) then 
            return 
        end if

        ! Remove any square brackets if present
        call CompareStringWithCharacter(stringval, matstart, templ)
        check = check .or. templ 
        call CompareStringWithCharacter(stringval, matend, templ)
        check = check .or. templ 
        

        ! Count the amount of row delimiters to obtain number of rows
        call CompareStringWithCharacter(stringval, rowdel, templ)
        nrow = count(templ)+1 ! always number + 1
        check = check .or. templ ! remove delimiters

        ! Count the amount of row separators to obtain number of columns
        call CompareStringWithCharacter(stringval, rowsep, templ)
        ncol = count(templ)/nrow + 1
        check = check .or. templ ! remove separators

        ! Check remainder of division
        remdivcheck = mod(count(templ), nrow)
        if (remdivcheck .ne. 0) then 
            ! Something wrong with input, return
            return 
        end if

        ! Trim the string
        check = .not. check 
        do i = 1, stringlength 
            if (check(i)) then 
                ! Add
                tempc(i:i) = stringval(i:i)
            else
                ! Replace by whitespace
                tempc(i:i) = ' ' 
            end if 
        end do 

        ! Read
        !=====
        ! Allocate
        allocate(intval(nrow, ncol))
        allocate(tempi(nrow*ncol))

        ! Read
        read (tempc, *, iostat=readstatus) tempi

        ! Check if read succeeded
        if (readstatus == 0) then 
            ! Set logical
            islegal = .true. ! I'll make it legal

            ! Reformat
            k = 1
            do i = 1, nrow 
                do j = 1, ncol 
                    intval(i, j) = tempi(k)
                    k = k + 1
                end do 
            end do
        end if

    end subroutine

    ! Scalar integer extraction
    subroutine ExtractRealFromString0D(stringval, realval, islegal)

        ! Description
        !============
        ! This routine extracts an integer from a string. Some things 
        ! are checked first: 
        ! - square brackets must be removed first before extraction, if
        !   present
        ! - It must be checked whether the integer is actually an 
        !   integer 
        !
        ! The logical islegal that is returned indicates whether a 
        ! correct value could be extracted (false if this is not the 
        ! case).

        ! Declare variables
        !==================
        ! Arguments
        character(:), allocatable, intent(in)   :: stringval
        real(R8), intent(out)                   :: realval 

        ! Auxiliary
        integer(I8)                             :: stringlength, &
            readstatus
        logical                                 :: islegal ! My lord, is that legal?
        character(:), allocatable               :: tempc 

        logical, allocatable                    :: check(:), templ(:)

        ! Loop
        integer(I8)                             :: i, k

        ! Initialize
        !===========
        ! Length of the string
        stringlength = len(stringval) 

        ! Allocate
        allocate(check(stringlength), templ(stringlength))
        allocate(character(stringlength) :: tempc)

        ! Initialize
        islegal     = .false. 
        check(:)    = .false. 

        ! Checks
        !=======
        ! Is there a row separator? 
        call CompareStringWithCharacter(stringval, rowsep, templ) 
        if (any(templ)) then 
            return 
        end if

        ! Is there a column separator? 
        call CompareStringWithCharacter(stringval, rowdel, templ) 
        if (any(templ)) then 
            return 
        end if

        ! Remove any square brackets if present
        call CompareStringWithCharacter(stringval, matstart, templ)
        check = check .or. templ 
        call CompareStringWithCharacter(stringval, matend, templ)
        check = check .or. templ 
        check = .not. check 

        ! Trim the string
        do i = 1, stringlength 
            if (check(i)) then 
                ! Add
                tempc(i:i) = stringval(i:i)
            else
                ! Replace by whitespace
                tempc(i:i) = ' ' 
            end if 
        end do 

        ! Read
        !=====
        ! Read
        read (tempc, *, iostat=readstatus) realval

        ! Check if read succeeded
        if (readstatus == 0) then 
            islegal = .true. ! I'll make it legal
        end if

    end subroutine

    ! Array integer extraction
    subroutine ExtractRealFromString1D(stringval, realval, islegal)

        ! Description
        !============
        ! This routine extracts an integer array from a string. 
        ! Some things are checked first: 
        ! - square brackets must be removed first before extraction, if
        !   present
        ! - It must be checked whether the integer is actually an 
        !   integer 
        ! - The amount of comma's are counted as this indicates the 
        !   number of elements of the vector.
        !
        ! The logical islegal that is returned indicates whether a 
        ! correct value could be extracted (false if this is not the 
        ! case).

        ! Declare variables
        !==================
        ! Arguments
        character(:), allocatable, intent(in)   :: stringval
        real(R8), allocatable, intent(out)      :: realval(:) 

        ! Auxiliary
        integer(I8)                             :: stringlength, &
            readstatus, nval
        logical                                 :: islegal ! My lord, is that legal?
        character(:), allocatable               :: tempc 

        logical, allocatable                    :: check(:), templ(:)

        ! Loop
        integer(I8)                             :: i, k
        

        ! Initialize
        !===========
        ! Length of the string
        stringlength = len(stringval) 

        ! Allocate
        allocate(check(stringlength), templ(stringlength))
        allocate(character(stringlength) :: tempc)

        ! Initialize
        islegal     = .false. 
        check(:)    = .false. 

        ! Checks
        !=======
        ! Is there a column separator? 
        call CompareStringWithCharacter(stringval, rowdel, templ) 
        if (any(templ)) then 
            return 
        end if

        ! Remove any square brackets if present
        call CompareStringWithCharacter(stringval, matstart, templ)
        check = check .or. templ 
        call CompareStringWithCharacter(stringval, matend, templ)
        check = check .or. templ 
        
        ! Count the amount of row delimiters to obtain size
        call CompareStringWithCharacter(stringval, rowsep, templ)
        nval = count(templ)+1

        ! Eliminate row delimiters 
        check = check .or. templ

        ! Trim the string
        check = .not. check 
        do i = 1, stringlength 
            if (check(i)) then 
                ! Add
                tempc(i:i) = stringval(i:i)
            else
                ! Replace by whitespace
                tempc(i:i) = ' ' 
            end if 
        end do 

        ! Read
        !=====
        ! Allocate
        allocate(realval(nval))

        ! Read
        read (tempc, *, iostat=readstatus) realval

        ! Check if read succeeded
        if (readstatus == 0) then 
            islegal = .true. ! I'll make it legal
        end if

    end subroutine

    ! Matrix integer extraction
    subroutine ExtractRealFromString2D(stringval, intval, islegal)

        ! Description
        !============
        ! This routine extracts a real matrix from a string. 
        ! Some things are checked first: 
        ! - square brackets must be removed first before extraction, if
        !   present
        ! - The amount of row and column separators are counted as this 
        !   yields the number of rows and columns
        !
        ! The logical islegal that is returned indicates whether a 
        ! correct value could be extracted (false if this is not the 
        ! case).

        ! Declare variables
        !==================
        ! Arguments
        character(:), allocatable, intent(in)   :: stringval
        real(R8), allocatable, intent(out)      :: intval(:, :) 

        ! Auxiliary
        integer(I8)                             :: stringlength, &
            readstatus, nval, nrow, ncol, remdivcheck
        logical                                 :: islegal ! My lord, is that legal?
        character(:), allocatable               :: tempc 

        logical, allocatable                    :: check(:), templ(:)

        real(R8), allocatable                   :: tempr(:)

        ! Loop
        integer(I8)                             :: i, j, k
        

        ! Initialize
        !===========
        ! Length of the string
        stringlength = len(stringval) 

        ! Allocate
        allocate(check(stringlength), templ(stringlength))
        allocate(character(stringlength) :: tempc)

        ! Initialize
        islegal     = .false. 
        check(:)    = .false. 

        ! Checks
        !=======
        ! Remove any square brackets if present
        call CompareStringWithCharacter(stringval, matstart, templ)
        check = check .or. templ 
        call CompareStringWithCharacter(stringval, matend, templ)
        check = check .or. templ 
        
        ! Count the amount of row delimiters to obtain number of rows
        call CompareStringWithCharacter(stringval, rowdel, templ)
        nrow = count(templ)+1 ! always number + 1
        check = check .or. templ ! remove delimiters

        ! Count the amount of row separators to obtain number of columns
        call CompareStringWithCharacter(stringval, rowsep, templ)
        ncol = count(templ)/nrow + 1
        check = check .or. templ ! remove separators

        ! Check remainder of division
        remdivcheck = mod(count(templ), nrow)
        if (remdivcheck .ne. 0) then 
            ! Something wrong with input, return
            return 
        end if

        ! Trim the string
        check = .not. check 
        do i = 1, stringlength 
            if (check(i)) then 
                ! Add
                tempc(i:i) = stringval(i:i)
            else
                ! Replace by whitespace
                tempc(i:i) = ' ' 
            end if 
        end do 

        ! Read
        !=====
        ! Allocate
        allocate(intval(nrow, ncol))
        allocate(tempr(nrow*ncol))

        ! Read
        read (tempc, *, iostat=readstatus) tempr

        ! Check if read succeeded
        if (readstatus == 0) then 
            ! Set logical
            islegal = .true. ! I'll make it legal

            ! Reformat
            k = 1
            do i = 1, nrow 
                do j = 1, ncol 
                    intval(i, j) = tempr(k)
                    k = k + 1
                end do 
            end do
        end if

    end subroutine

    ! Auxiliary routines
    subroutine CompareStringWithCharacter(s, c, l) 

        ! Description
        !============
        ! This routine compares the characters  of a string s to the 
        ! character given in c. The logical l is of length len(s) and 
        ! indicates which characters are the same. 

        ! Declare variables
        !==================
        ! Arguments
        character(:), allocatable, intent(in)   :: s 
        character, intent(in)                   :: c 
        logical, allocatable, intent(inout)     :: l(:)

        ! Loop
        integer(I8)                             :: i 

        ! Initialize
        !===========
        ! Allocate
        if (allocated(l)) then  
            deallocate(l) 
            allocate(l(len(s)))
        end if

        ! Initialize
        l(:) = .false.

        ! Compare
        !========
        do i = 1, len(s) 
            if (s(i:i) == c) then 
                l(i) = .true. 
            end if
        end do

    end subroutine

end module