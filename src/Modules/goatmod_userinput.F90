!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains the definition of the types and routines used to
! process user input for goat. The goat input file is eventually a 
! set of key-value pairs. The following conventions are adopted:
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
!   Here, the arrays of integers and scalars should be delimited by 
!   square brackets. Scalar values may also be delimited with square 
!   brackets, but this is not required. Values in arrays should be 
!   separated by commas and semicolons. The commas are used to delimit 
!   row-wise entries, the semicolons indicate the columns. For one-
!   dimensional arrays, only commas should be present. Note that the
!   precision of reals and integers is fixed. This may be extended in 
!   the future by defining subroutines that deal with different 
!   precision specifications. 
!
! Note that if a header is present, it will be ignored as long as it 
! does not contain any quotation marks. 

module goatmod_userinput

    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_readwrite

    ! The usual
    implicit none
    save
    public 

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    ! Abstract option type
    type, abstract :: OptionsUDT  

        ! General abstract type for options. Only states which routines 
        ! should be provided (set)

    contains 

        procedure(SetOptionsINT), deferred      :: Set 
        procedure(SetDefaultsINT), deferred     :: SetDefaults

    end type 

    ! General goat options
    type, extends(OptionsUDT) :: GoatoptionsUDT

        ! Structure containing the options for goat. The following 
        ! fields are present:
        ! - driver:     specify the driver to be used in the goat 
        !               program (see Goat.F90 for the options)
        ! - 
        character(:), allocatable   :: driver ! driver to be taken for goat
        integer                     :: itmax 
        character(:), allocatable   :: filepath ! file path to options file

    contains

        ! Routines to manipulate the options
        procedure   :: Set              => SetGoatOptions
        procedure   :: SetDefaults      => SetDefaultGoatOptions
        

    end type

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    ! Abstract interfaces
    abstract interface 

        ! Options
        subroutine SetOptionsINT(options, filepath)

            import :: OptionsUDT 
            class(OptionsUDT)           :: options
            character(*), intent(in)    :: filepath

        end subroutine

        subroutine SetDefaultsINT(options)

            import :: OptionsUDT 
            class(OptionsUDT) :: options 

        end subroutine

    end interface

    contains 

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                            Option setters                        !
    !------------------------------------------------------------------!

    ! Goat options routines
    subroutine SetGoatOptions(options, filepath)

        ! Declare variables
        !==================
        ! Arguments
        class(GoatoptionsUDT)       :: options
        character(*), intent(in)    :: filepath
        
        ! Set default values
        !===================
        options%filepath = filepath
        call SetDefaultGoatOptions(options)

        ! User-specified options
        !=======================
        call ReadGoatOptions(options)

    end subroutine

    subroutine SetDefaultGoatOptions(options)

        ! Description
        !============
        ! This routine sets the default goat options. These can be later 
        ! overridden by reading in the user-specified options. 

        ! Declare variables
        !==================
        ! Arguments
        class(GoatoptionsUDT)       :: options 

        ! Set default options
        !====================
        options%driver = 'GD' ! default driver: grid deformation
        options%itmax = 10

    end subroutine

    ! Goat options reader
    subroutine ReadGoatOptions(options)

        ! Description
        !============
        ! This routine reads in the goat options from a file of which 
        ! the full path should be given in options%filepath. The default
        ! options should have already been set at this point, as this 
        ! routine will only overwrite options that are present in the 
        ! user-specified options file. If no options file is present, 
        ! nothing is read in and a message will be shown. 

        ! Notes
        !======
        ! Note 1: this routine starts with opening the file and checking
        ! whether it exists. Afterwards, for each user-adjustable option 
        ! we read through the entire file to 

        ! Declare variables
        !==================
        ! Arguments
        type(GoatoptionsUDT)            :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: thisline, field
        integer, parameter              :: fid = 10 
        logical                         :: reachedeof

        ! Initialize
        !===========
        ! Variables
        reachedeof = .false. 

        ! Open the file, check if it exists
        open(unit=fid, file=options%filepath, status='old', &
            iostat=openstatus)

        if (openstatus > 0) then 
            ! Something wrong when reading file - continue with default
            ! values
            print *, 'ReadGoatOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadGoatOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Driver
        field = 'GOAToptions.driver'
        call ExtractOptionValueCharacter(fid, field, options%driver)

        ! Max number of iterations
        field = 'GOAToptions.itmax'
        call ExtractOptionValueInteger0D(fid, field, options%itmax)

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)


    end subroutine

    !------------------------------------------------------------------!
    !                           File processing                        !
    !------------------------------------------------------------------!

    ! Option extractor implementations
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
            call ExtractIntegerFromString(temp, tempi, islegal)

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


        ! Delimiters etc
        character, parameter                 :: commentchar  = '#'
        character, parameter                 :: delimiter    = "'"

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

    ! Extraction routines
    subroutine ExtractIntegerFromString(stringval, intval, islegal)

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

        ! Delimiters
        character, parameter                    :: matstart = '['
        character, parameter                    :: matend   = ']'
        character, parameter                    :: decpoint = '.'
        character, parameter                    :: rowsep   = ','
        character, parameter                    :: colsep   = ';'

        ! Initialize
        !===========
        ! Length of the string
        stringlength = len(stringval) 

        ! Allocate
        allocate(check(stringlength), templ(stringlength))

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
        call CompareStringWithCharacter(stringval, colsep, templ) 
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
        k = 1
        allocate(character(count(check)) :: tempc) 
        do i = 1, stringlength 
            if (check(i)) then 
                tempc(k:k) = stringval(i:i)
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

    


end module goatmod_userinput