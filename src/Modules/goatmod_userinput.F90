!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!



module goatmod_userinput

    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_readwrite
    use mod_inputfileparser

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
        ! - facelabelsubfrom:   labels in the original grid file to 
        !                       replace by other files
        ! - facelabelsubto:     replacement labels

        character(:), allocatable   :: driver ! driver to be taken for goat
        integer                     :: itmax 
        character(:), allocatable   :: filepath ! file path to options file

        ! Mappings
        integer(I8), allocatable    :: facelabelsubfrom(:) ! labels

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
        ! General options
        options%driver = 'GD' ! default driver: grid deformation
        options%itmax = 10

        ! Mappings
        allocate(options%facelabelsubfrom(0)) ! assume no mappings

    end subroutine

    !------------------------------------------------------------------!
    !                            Option readers                        !
    !------------------------------------------------------------------!

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

        ! Mappings
        field = 'GOAToptions.GDtoGA.facelabelsubfrom'
        call ExtractOptionValueInteger1D(fid, field, &
            options%facelabelsubfrom)
        

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)


    end subroutine

end module goatmod_userinput