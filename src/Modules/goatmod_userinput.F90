!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains the definition of the types and routines used to
! process user input for goat. 

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
        character(C32)              :: driver ! driver to be taken for goat
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

    ! Goat options routines
    subroutine SetGoatOptions(options, filepath)

        ! Declare variables
        !==================
        ! Arguments
        class(GoatoptionsUDT)       :: options
        character(*), intent(in)    :: filepath
        character(len=:), allocatable      :: thisline

        ! Auxiliary
        integer, parameter          :: fid = 10 
        logical                     :: reachedeof

        ! Initialize
        !===========
        reachedeof = .false. 

        ! Set default values
        !===================
        options%filepath = filepath
        call SetDefaultGoatOptions(options)

        ! User-specified options
        !=======================
        open(unit=fid, file=options%filepath)
        do while (.not. reachedeof)
            call ReadSingleLine(fid, thisline, reachedeof)
        end do
        close(unit=fid)

        

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

    end subroutine

end module goatmod_userinput