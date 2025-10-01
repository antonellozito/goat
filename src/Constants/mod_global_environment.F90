!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module defines variables depending on compiler options etc. 
! This cannot avoid all compiler directives in the code, but at least 
! some things (e.g. saving/loading behavior dependent on whether other
! modules/codes are/are not used) can be accounted for from here. 

module mod_global_environment

    implicit none 
    save 
    public 

    ! SOLPS-related switches
    !=======================
#ifdef SOLPS
    logical, parameter          :: solps = .true. 
#else 
    logical, parameter          :: solps = .false.
#endif 
    character(*), parameter     :: solps_inputfilepath = 'GOAToptions.dat', &
        solps_gridfilepath              = 'traduit.out.b2us', &
        solps_structurefilepath         = 'structure.dat', &
        solps_magneticfieldfilepath     = 'rzpsi.dat', &
        solps_statefilepath             = 'b2fplasmf', &
        solps_writefilepath             = 'traduit_goat.out.b2us', & 
        solps_outputfilepath            = './output/'


    contains 

    ! Pre-run routine for solps
    subroutine SolpsPreamble(filename)

        ! Description
        !============
        ! Do solps related preamble, such as writing a .prt file for 
        ! a specific driver routine (filename should be given without
        ! .prt extension!). Also, we generate some directories here
        ! using the command line. 

        ! Declare variables
        !==================
        ! Arguments
        character(*), intent(in)        :: filename 

        ! Auxiliary
        integer, parameter              :: fid = 10
        integer                         :: iostat 

        ! Construct .prt file
        !====================
        ! Make the goat.prt file to check 
        open(unit=fid, file=filename//'.prt', status='new', iostat=iostat)

        ! If this was not possible, call error and exit
        close(unit=fid)
        if (iostat /= 0) then 
            ! Print and stop
            print *, 'SolpsPreamble: goat.prt already exists, exiting...'
            stop
        end if 

        ! Construct the output directory (one level above current directory
        ! to comply with run/baserun assumption)
        call execute_command_line('mkdir ./output')
        
    end subroutine


end module
