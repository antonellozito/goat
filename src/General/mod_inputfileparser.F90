!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains the definition of the types and routines used to
! process user input for goat. The goat input file is eventually a 
! set of key-value pairs. The intention is to read in options or a small
! amount of data (so no grid data or others, there are dedicated 
! routines for that). 

! The following conventions are adopted:
! 
! - The key comes before the value
! - The key is delimited by single quotation marks, as is the value. 
!   This means that at least four quotation marks have to be present in 
!   order to interpret a line. If more are present, it is assumed that 
!   the key is between the first two quotations marks and the value is 
!   between the last two. Other quotation marks are assumed to be typos
!   and are therefore ignored. Notice that keys and values therefore 
!   should not contain any (single) quotation marks! If strings have to 
!   read in, simply specify in the reading option that the data should 
!   be formatted to a string. 
! - Lines starting with '#' are interpreted as comments and are ignored
! - The following types can be read in: 
!       * scalar integers (kind: I8)
!       * arrays of integers (one- or two-dimensional) (kind: I8)
!       * scalar reals (kind: R8)
!       * arrays of reals (one- or two-dimensional) (kind: R8)
!       * arrays of characters (one-dimensional, case-sensitive)
!       * scalar logicals
!   Here, the (arrays of) integers and scalars can be delimited by 
!   square brackets (not required). Values in arrays should be 
!   separated by commas and semicolons. The commas are used to separate 
!   row-wise entries, the semicolons indicate the end of rows. For one-
!   dimensional arrays, only commas should be present. Note that the
!   precision of reals and integers is fixed. This may be extended in 
!   the future by defining subroutines that deal with different 
!   precision specifications. 
!

! Notes
!======
! Note 1: if a header is present, it will be ignored as long as it 
! does not contain any quotation marks. 

! Note 2: it is important to remark that array and matrix dimensions are 
! derived from the presence of commas and semicolons. This makes 
! incrementing some options more easy, but is not suited for large data!
! This implementation may be extended for different input format, e.g. 
! where the sizes of the array are parsed beforehand and where one 
! makes assumptions on the data format used. However, other routines
! may be more useful for this kind of data. 

! Note 3: though checks are performed on the data to see if everything 
! is read in correctly, there is no guarantee that all possible errors
! will be caught. Visual checks by the user can be done based on the 
! printed output. Matrices will be printed row per row. 

! Note 4: concerning 2D matrices, the following assumptions are made:
! - The number of row delimiters + 1 is equal to the number of columns.
!   This means that semicolons should *not* appear at the end of the 
!   matrix! this is not checked for explicitly
! - The same holds for row separators: the amount of elements in one
!   row is assumed to be equal to the number of row separators + 1.

! Note 5: the definition of delimiters etc and which characters are 
! used for those is given in the general module mod_specialchars. If 
! desired, one can change the definition there.

! Note 6: arrays may be specified using the array constructor defined
! by the 'veccon' symbol in mod_specialchars. For example, if veccon is 
! the colon, then '1:10' would specify an array going from 1 to 10. Note 
! that it is assumed that the increment is always one! To repeat a 
! number, one can use the 'repeatchar' symbol from mod_specialchars. If 
! this character is '*', then '5*10' would repeat the number 5 ten 
! times. 

module mod_inputfileparser

    ! Initialize
    !============
    use mod_precision
    use mod_specialchars 
    use mod_readwrite 

    contains 

    !------------------------------------------------------------------!
    !                           File processing                        !
    !------------------------------------------------------------------!

    ! Extract scalar integer
    subroutine ExtractOptionValueInteger0D(fid, key, val)

        ! Description
        !============
        ! Main driver to extract scalar integer value from a formatted
        ! .dat input file. It is assumed that the file has been opened 
        ! and exists (its unit is then given by fid). 

        ! Declare variables
        !==================
        ! Arguments
        integer, intent(in)                     :: fid
        character(:), allocatable, intent(in)   :: key 
        integer(I8)                             :: val 

        ! Auxiliary
        logical                                 :: islegal, isfound
        integer(I8)                             :: tempi 
        character(:), allocatable               :: temp

        ! Initialize
        !===========
        islegal = .false. 

        ! Search
        !=======
        ! Get the value belonging to the key in character array format
        call GetValueWithKey(fid, key, temp, isfound)

        ! Check if it is found, otherwise exit
        if (isfound) then 

            ! Attribute
            call ExtractIntegerFromString0D(temp, tempi, islegal)

            ! Check
            if (islegal) then 
                ! Attribute
                val = tempi

                ! Print
                print *, key, ' = ', val 
            else 
                ! Print
                print *, key, ' = ', val, '(default, illegal value in file)'
            end if 

        else 

            ! Print  
            print *, key, ' = ', val, ' (default, could not find in file)'

        end if

    end subroutine

    ! Extract array integer
    subroutine ExtractOptionValueInteger1D(fid, key, val)

        ! Description
        !============
        ! Main driver to extract array integer values from a formatted
        ! .dat input file. It is assumed that the file has been opened 
        ! and exists (its unit is then given by fid). 

        ! Declare variables
        !==================
        ! Arguments
        integer, intent(in)                     :: fid
        character(:), allocatable, intent(in)   :: key 
        integer(I8), allocatable                :: val(:) 

        ! Auxiliary
        logical                                 :: islegal, isfound
        integer(I8), allocatable                :: tempi(:) 
        character(:), allocatable               :: temp

        ! Initialize
        !===========
        ! Set logicals
        islegal = .false.

        ! Search
        !=======
        ! Get the value belonging to the key in character array format
        call GetValueWithKey(fid, key, temp, isfound)

        ! Check if it is found, otherwise exit
        if (isfound) then 

            ! Attribute
            call ExtractIntegerFromString1D(temp, tempi, islegal)

            ! Check
            if (islegal) then 

                ! Deallocate if already allocated - can't know the size
                if (allocated(val)) then 
                    deallocate(val) 
                end if 

                ! Attribute
                allocate(val(size(tempi)))
                val = tempi

                ! Print
                print *, key, ' = ', val 
            else 
                ! Print
                print *, key, ' = ', val, '(default, illegal value in file)'
            end if 

        else 

            ! Print  
            print *, key, ' = ', val, ' (default, could not find in file)'

        end if

    end subroutine

    ! Extract matrix integer
    subroutine ExtractOptionValueInteger2D(fid, key, val)

        ! Description
        !============
        ! Main driver to extract matrix integer values from a formatted
        ! .dat input file. It is assumed that the file has been opened 
        ! and exists (its unit is then given by fid). 

        ! Declare variables
        !==================
        ! Arguments
        integer, intent(in)                     :: fid
        character(:), allocatable, intent(in)   :: key 
        integer(I8), allocatable                :: val(:, :) 

        ! Auxiliary
        logical                                 :: islegal, isfound
        integer(I8), allocatable                :: tempi(:, :) 
        character(:), allocatable               :: temp

        ! Loop
        integer(I8)                             :: i

        ! Initialize
        !===========
        ! Set logicals
        islegal = .false.

        ! Search
        !=======
        ! Get the value belonging to the key in character array format
        call GetValueWithKey(fid, key, temp, isfound)

        ! Check if it is found, otherwise exit
        if (isfound) then 

            ! Attribute
            call ExtractIntegerFromString2D(temp, tempi, islegal)

            ! Check
            if (islegal) then 

                ! Deallocate if already allocated - can't know the size
                if (allocated(val)) then 
                    deallocate(val) 
                end if 

                ! Attribute
                allocate(val(size(tempi, 1), size(tempi, 2)))
                val = tempi

                ! Print
                print *, key, ' = ', val(1, :)
                do i = 2, size(val, 1) ! print row per row
                    print *, val(i, :)
                end do 
            else 
                ! Print
                print *, key, ' = ', val, '(default, illegal value in file)'
            end if 

        else 

            ! Print  
            print *, key, ' = ', val, ' (default, could not find in file)'

        end if

    end subroutine

    ! Extract scalar real
    subroutine ExtractOptionValueReal0D(fid, key, val)

        ! Description
        !============
        ! Main driver to extract scalar real value from a formatted
        ! .dat input file. It is assumed that the file has been opened 
        ! and exists (its unit is then given by fid). 

        ! Declare variables
        !==================
        ! Arguments
        integer, intent(in)                     :: fid
        character(:), allocatable, intent(in)   :: key 
        real(R8)                                :: val 

        ! Auxiliary
        logical                                 :: islegal, isfound
        real(R8)                                :: tempr
        character(:), allocatable               :: temp

        ! Initialize
        !===========
        islegal = .false. 

        ! Search
        !=======
        ! Get the value belonging to the key in character array format
        call GetValueWithKey(fid, key, temp, isfound)

        ! Check if it is found, otherwise exit
        if (isfound) then 

            ! Attribute
            call ExtractRealFromString0D(temp, tempr, islegal)

            ! Check
            if (islegal) then 
                ! Attribute
                val = tempr

                ! Print
                print *, key, ' = ', val 
            else 
                ! Print
                print *, key, ' = ', val, '(default, illegal value in file)'
            end if 

        else 

            ! Print  
            print *, key, ' = ', val, ' (default, could not find in file)'

        end if

    end subroutine

    ! Extract array real
    subroutine ExtractOptionValueReal1D(fid, key, val)

        ! Description
        !============
        ! Main driver to extract array real values from a formatted
        ! .dat input file. It is assumed that the file has been opened 
        ! and exists (its unit is then given by fid). 

        ! Declare variables
        !==================
        ! Arguments
        integer, intent(in)                     :: fid
        character(:), allocatable, intent(in)   :: key 
        real(R8), allocatable                   :: val(:) 

        ! Auxiliary
        logical                                 :: islegal, isfound
        real(R8), allocatable                   :: tempr(:) 
        character(:), allocatable               :: temp

        ! Initialize
        !===========
        ! Set logicals
        islegal = .false.

        ! Search
        !=======
        ! Get the value belonging to the key in character array format
        call GetValueWithKey(fid, key, temp, isfound)

        ! Check if it is found, otherwise exit
        if (isfound) then 

            ! Attribute
            call ExtractRealFromString1D(temp, tempr, islegal)

            ! Check
            if (islegal) then 

                ! Deallocate if already allocated - can't know the size
                if (allocated(val)) then 
                    deallocate(val) 
                end if 

                ! Attribute
                allocate(val(size(tempr)))
                val = tempr

                ! Print
                print *, key, ' = ', val 
            else 
                ! Print
                print *, key, ' = ', val, '(default, illegal value in file)'
            end if 

        else 

            ! Print  
            print *, key, ' = ', val, ' (default, could not find in file)'

        end if

    end subroutine

    ! Extract matrix real
    subroutine ExtractOptionValueReal2D(fid, key, val)

        ! Description
        !============
        ! Main driver to extract matrix integer values from a formatted
        ! .dat input file. It is assumed that the file has been opened 
        ! and exists (its unit is then given by fid). 

        ! Declare variables
        !==================
        ! Arguments
        integer, intent(in)                     :: fid
        character(:), allocatable, intent(in)   :: key 
        real(R8), allocatable                   :: val(:, :) 

        ! Auxiliary
        logical                                 :: islegal, isfound
        real(R8), allocatable                   :: tempr(:, :) 
        character(:), allocatable               :: temp

        ! Loop
        integer(I8)                             :: i

        ! Initialize
        !===========
        ! Set logicals
        islegal = .false.

        ! Search
        !=======
        ! Get the value belonging to the key in character array format
        call GetValueWithKey(fid, key, temp, isfound)

        ! Check if it is found, otherwise exit
        if (isfound) then 

            ! Attribute
            call ExtractRealFromString2D(temp, tempr, islegal)

            ! Check
            if (islegal) then 

                ! Deallocate if already allocated - can't know the size
                if (allocated(val)) then 
                    deallocate(val) 
                end if 

                ! Attribute
                allocate(val(size(tempr, 1), size(tempr, 2)))
                val = tempr

                ! Print
                print *, key, ' = ', val(1, :)
                do i = 2, size(val, 1) ! print row per row
                    print *, val(i, :)
                end do 
            else 
                ! Print
                print *, key, ' = ', val, '(default, illegal value in file)'
            end if 

        else 

            ! Print  
            print *, key, ' = ', val, ' (default, could not find in file)'

        end if

    end subroutine

    ! Extract character array
    subroutine ExtractOptionValueCharacter(fid, key, val)

        ! Description
        !============
        ! See ExtractOptionValueInteger0D, but now for character array.

        ! Declare variables
        !==================
        ! Arguments
        integer, intent(in)                     :: fid
        character(:), allocatable, intent(in)   :: key 
        character(:), allocatable               :: val 

        ! Auxiliary
        logical                                 :: isfound
        character(:), allocatable               :: temp

        ! Search
        !=======
        ! Get the value belonging to the key in character array format
        call GetValueWithKey(fid, key, temp, isfound)

        ! Check if it is found, otherwise exit
        if (isfound) then 

            ! Attribute
            val = temp 

            ! Print
            print *, key, ' = ', val

        else 

            ! Print  
            print *, key, ' = ', val, ' (default, could not find in file)'

        end if

    end subroutine

    ! Extract scalar logical
    subroutine ExtractOptionValueLogical0D(fid, key, val)

        ! Description
        !============
        ! Main driver to extract scalar logical value from a formatted
        ! .dat input file. It is assumed that the file has been opened 
        ! and exists (its unit is then given by fid). 

        ! Declare variables
        !==================
        ! Arguments
        integer, intent(in)                     :: fid
        character(:), allocatable, intent(in)   :: key 
        integer(I8)                             :: valint 

        ! Auxiliary
        logical                                 :: islegal, isfound, &
            val
        integer(I8)                             :: tempi 
        character(:), allocatable               :: temp

        ! Initialize
        !===========
        islegal = .false. 

        ! Search
        !=======
        ! Get the value belonging to the key in character array format
        call GetValueWithKey(fid, key, temp, isfound)

        ! Check if it is found, otherwise exit
        if (isfound) then 

            ! Attribute
            call ExtractIntegerFromString0D(temp, tempi, islegal)

            ! Check
            if (islegal) then 
                ! Attribute
                valint = tempi
                if (valint .lt. 0) then 
                    ! Value is true
                    val = .true.
                else 
                    val = .false.
                end if

                ! Print
                print *, key, ' = ', val
            else 
                ! Print
                print *, key, ' = ', val, '(default, illegal value in file)'
            end if 

        else 

            ! Print  
            print *, key, ' = ', val, ' (default, could not find in file)'

        end if

    end subroutine

    ! Value getter (character array format)
    subroutine GetValueWithKey(fid, key, val, isfound)

        ! Description
        !============
        ! Get the value of a key (if it exists). Until the key is found
        ! or the EOF is reached, a line (record) from the file is read 
        ! in and subsequently processed. 

        ! Declare variables
        !==================
        ! Arguments
        integer, intent(in)                     :: fid
        character(:), allocatable, intent(in)   :: key 
        character(:), allocatable               :: val
        logical                                 :: isfound  
        
        ! Auxiliary
        logical                                 :: reachedeof, haspair
        character(:), allocatable               :: tempkey, tempvalue, &
            thisline

        ! Initialize
        !===========
        ! Rewind to make sure we start at the beginning
        rewind(fid) 

        ! Set logicals
        isfound     = .false. 
        reachedeof  = .false. 
        haspair     = .false. 

        ! Loop
        !=====
        do while ( (.not. reachedeof) .and. (.not. isfound) ) 

            ! Read line
            call ReadSingleLine(fid, thisline, reachedeof)

            ! Check if it contains a key-value pair
            call GetKeyValuePair(thisline, tempkey, tempvalue, haspair)
            if (haspair) then 
                ! Check if the key matches the given key
                if (key == tempkey) then 
                    ! Found the correct key, exit
                    isfound = .true. 

                    ! Set value
                    val = tempvalue
                end if
            end if
        end do

    end subroutine

    ! Key-value pair getter
    subroutine GetKeyValuePair(curline, key, val, haspair)

        ! Description
        !============
        ! This routine extracts from a single line (i.e. character 
        ! array) the key and value if they are present. To do so, the 
        ! rules specified in the beginning of this module are followed. 
        ! To allow for easy adaptation, the delimiter and comment
        ! characters are stored in variables (i.e. indirection is used)

        ! Declare variables
        !==================
        ! Arguments
        character(:), allocatable, intent(in)   :: curline 
        character(:), allocatable, intent(out)  :: key, val 
        logical                                 :: haspair 

        ! Auxiliary
        integer(I8)                             :: ndelim

        integer(I8), allocatable                :: delimloc(:)

        logical, allocatable                    :: isdelim(:)

        ! Loop
        integer(I8)                             :: k

        ! Initialize
        !===========
        ! Set logicals
        haspair = .false. 

        ! Check where the delimiters occur
        allocate(isdelim(len(curline))) 

        ! Count the number of quotation marks
        do k = 1, len(curline)
            isdelim(k) = curline(k:k) == delimiter
        end do
        ndelim = count(isdelim)

        ! Check string
        !=============
        if ( ((curline(1:1) == commentchar)) .or. (ndelim < 4) ) then ! Comment or bogus line
            ! No key-value pair present
        else

            ! Extract the location of the delimiters
            allocate(delimloc(ndelim))
            delimloc = pack( [(k, k = 1, len(curline))], isdelim)

            ! Get the key-value pair
            key = curline(delimloc(1)+1:delimloc(2)-1)
            val = curline(delimloc(ndelim-1)+1:delimloc(ndelim)-1)

            ! Set logicals
            haspair = .true.

            deallocate(delimloc)

        end if 

        ! Housekeeping
        !=============
        deallocate(isdelim)

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

        ! Is there a vector constructor?
        call CompareStringWithCharacter(stringval, veccon, templ) 
        if (any(templ)) then 
            return 
        end if

        ! Is there a repeater character?
        call CompareStringWithCharacter(stringval, repeatchar, templ) 
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
            readstatus, nval, ntempval, tval, tempi 
        integer(I8), allocatable                :: index(:), &
            tempintval(:)
        logical                                 :: islegal ! My lord, is that legal?
        character(:), allocatable               :: tempc 

        logical, allocatable                    :: check(:), templ(:), &
            isveccon(:), isarraystart(:), isrepchar(:), isrepstart(:)

        ! Loop
        integer(I8)                             :: i, j, k
        

        ! Initialize
        !===========
        ! Length of the string
        stringlength = len(stringval) 

        ! Allocate
        allocate(check(stringlength), templ(stringlength), &
            isveccon(stringlength), index(stringlength), & 
            isrepchar(stringlength))
        allocate(character(stringlength) :: tempc)

        ! Initialize
        islegal     = .false. 
        check(:)    = .false. 
        isveccon(:) = .false. 
        isrepchar(:)    = .false.
        index(:)        = 0
        

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

        ! Check for any vector constructor symbols
        call CompareStringWithCharacter(stringval, veccon, isveccon)
        nval = nval + count(isveccon) ! +1 already accounted for 

        ! Eliminate vector constructors
        check = check .or. isveccon

        ! Check for any repeater symbols
        call CompareStringWithCharacter(stringval, repeatchar, isrepchar)
        nval = nval + count(isrepchar) ! +1 already accounted for

        ! Eliminate repeater symbols
        check = check .or. isrepchar

        ! Trim the string
        check = .not. check 
        k = 0
        do i = 1, stringlength 
            if (check(i)) then 
                ! Add
                tempc(i:i) = stringval(i:i)

                ! Checks for index
                if (i > 1) then 
                    ! Check if previous character was also a number. 
                    ! Otherwise, increase k
                    if (.not. check(i-1)) then 
                        k = k + 1
                    end if
                else 
                    ! First index, don't check previous
                    k = k + 1
                end if
                index(i) = k

            else
                ! Replace by whitespace
                tempc(i:i) = ' ' 
            end if 
        end do 

        ! Check if only whitespace remains. In that case, an empty array
        ! was specified
        call CompareStringWithCharacter(tempc, ' ', templ)
        if (all(templ)) then 
            nval = 0
        end if

        ! Read
        !=====
        ! Allocate
        allocate(tempintval(nval))

        ! Read (only if nonempty)
        if (nval == 0) then 
            islegal = .true.
            return 
        end if
        read (tempc, *, iostat=readstatus) tempintval

        ! Check if read succeeded
        if (readstatus .ne. 0) then 
            return 
        end if

        ! Check if we need to execute vector constructors or repeaters
        if (any(isveccon, 1) .or. any(isrepchar, 1)) then 
            ! First, we extract which number are start/end of array
            allocate(isarraystart(nval), isrepstart(nval))
            
            ! Initialize
            isarraystart(:) = .false. 
            isrepstart(:) = .false.

            ! Loop over isveccon to determine which one(s) are start
            do i = 2, stringlength-1 ! skip first and last entry - would/should be illegal
                if (isveccon(i)) then 
                    ! Check previous and next character (no whitespace allowed!)
                    if ((index(i-1) == 0) .or. (index(i+1) ==0)) then
                        ! Not allowed, return
                        return 
                    end if 

                    ! Set previous to true
                    isarraystart(index(i-1)) = .true.
                elseif (isrepchar(i)) then 
                    ! Check previous and next character (no whitespace allowed!)
                    if ((index(i-1) == 0) .or. (index(i+1) ==0)) then
                        ! Not allowed, return
                        return 
                    end if 

                    ! Set previous to true
                    isrepstart(index(i-1)) = .true.
                end if

            end do

            ! Determine the actual number of values for intval
            ntempval = nval
            nval = 0
            do i = 1, ntempval 
                if (isarraystart(i)) then 
                    ! Compute number of elements and add
                    nval = nval + tempintval(i+1)-tempintval(i)+1
                elseif (isrepstart(i)) then 
                    ! Compute number of elements and add
                    nval = nval + tempintval(i+1)
                elseif (isrepstart(i-1) .or. isarraystart(i-1)) then 
                    ! Skip addition, already added before
                else
                    ! Normal entry, +1
                    nval = nval + 1
                end if
            end do 

            ! Construct the final value vector
            k = 0
            allocate(intval(nval))
            do i = 1, ntempval 
                if (isarraystart(i)) then 
                    ! Compute number of elements and add
                    tempi = tempintval(i+1)-tempintval(i)+1
                    intval(k+1:k+tempi) = [(j, j = tempintval(i), tempintval(i+1))]
                    k = k + tempi
                elseif (isrepstart(i)) then 
                    ! Compute number of elements and add
                    tempi = tempintval(i+1)
                    intval(k+1:k+tempi) = tempintval(i)
                    k = k + tempi
                elseif (isrepstart(i-1) .or. isarraystart(i-1)) then 
                    ! Skip addition, already added before
                else
                    ! Normal entry
                    intval(k+1) = tempintval(i)
                    k = k + 1
                end if
            end do 
        else 
            ! Set intval equal to tempintval
            intval = tempintval 
        end if

        ! If we managed to get here, reading etc was successful
        islegal = .true.

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
        ! - Repeater characters and vector constructors are not yet 
        ! supported
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

        ! Is there a vector constructor?
        call CompareStringWithCharacter(stringval, veccon, templ) 
        if (any(templ)) then 
            return 
        end if

        ! Is there a repeater character?
        call CompareStringWithCharacter(stringval, repeatchar, templ) 
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

        ! Check if only whitespace remains. In that case, an empty array
        ! was specified
        call CompareStringWithCharacter(tempc, ' ', templ)
        if (all(templ)) then 
            nrow = 0
            ncol = 0
        end if

        ! Check remainder of division
        call CompareStringWithCharacter(stringval, rowsep, templ)
        if (nrow > 0) then 
            remdivcheck = mod(count(templ), nrow)
            if (remdivcheck .ne. 0) then 
                ! Something wrong with input, return
                return 
            end if
        end if

        ! Read
        !=====
        ! Allocate
        allocate(intval(nrow, ncol))
        allocate(tempi(nrow*ncol))

        ! Read (only if nonempty)
        if ( (nrow == 0) .or. (ncol == 0) ) then 
            islegal = .true.
            return 
        end if
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

        ! Is there a vector constructor?
        call CompareStringWithCharacter(stringval, veccon, templ) 
        if (any(templ)) then 
            return 
        end if

        ! Is there a repeater character?
        call CompareStringWithCharacter(stringval, repeatchar, templ) 
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
        ! - Repeater characters and vector constructors are not yet 
        ! supported
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

        ! Is there a vector constructor?
        call CompareStringWithCharacter(stringval, veccon, templ) 
        if (any(templ)) then 
            return 
        end if

        ! Is there a repeater character?
        call CompareStringWithCharacter(stringval, repeatchar, templ) 
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

        ! Check if only whitespace remains. In that case, an empty array
        ! was specified
        call CompareStringWithCharacter(tempc, ' ', templ)
        if (all(templ)) then 
            nval = 0
        end if

        ! Read
        !=====
        ! Allocate
        allocate(realval(nval))

        ! Read (only if nonempty)
        if (nval == 0) then 
            islegal = .true.
            return 
        end if
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
        ! - Repeater characters and vector constructors are not yet 
        ! supported
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
        ! Is there a vector constructor?
        call CompareStringWithCharacter(stringval, veccon, templ) 
        if (any(templ)) then 
            return 
        end if

        ! Is there a repeater character?
        call CompareStringWithCharacter(stringval, repeatchar, templ) 
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

        ! Check if only whitespace remains. In that case, an empty array
        ! was specified
        call CompareStringWithCharacter(tempc, ' ', templ)
        if (all(templ)) then 
            nrow = 0
            ncol = 0
        end if

        ! Check remainder of division
        call CompareStringWithCharacter(stringval, rowsep, templ)
        if (nrow > 0) then 
            remdivcheck = mod(count(templ), nrow)
            if (remdivcheck .ne. 0) then 
                ! Something wrong with input, return
                return 
            end if
        end if

        ! Read
        !=====
        ! Allocate
        allocate(intval(nrow, ncol))
        allocate(tempr(nrow*ncol))

        ! Read (only if nonempty)
        if ( (nrow == 0) .or. (ncol == 0) ) then 
            islegal = .true.
            return 
        end if
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
        else
        print *, stringval
        print *, check
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