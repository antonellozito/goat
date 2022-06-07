!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains all the user-defined input parameters and related
! routines for the grid deformation module. 

module gdmod_userinput

    ! Initialize
    !============
    ! Load modules
    use gdmod_types

    ! The usual
    implicit none
    save
    private 

    ! All functions
    public SetRunfileOptions

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!
    
    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    contains 

    subroutine SetRunfileOptions(options)
        
        ! Description
        !============
        ! Set the run options. The following fields have to be set:
        ! - runtype: define how to run the grid deformation. 'optimize'
        ! will use the optimization drivers (which can be further 
        ! specified using 'meth') - other options are not yet 
        ! implemented. 
        ! - gridtype: type of grid to be considered. Can only be 
        ! 'plasma' for now, though in the future e.g. EIRENE grids could
        ! also be considered.
        ! - meth: method to be used for the specific runtype. For 
        ! 'optimize', only 'KKT' is currently available. Here, the full
        ! KKT conditions of the grid deformation optimization problem 
        ! are set up and solved.
        ! - export: export the results in a certain format to a 
        ! specified file (can be true or false). See also the
        ! SetExportOptions subroutine in this file. 

        ! Declaration
        type (RunfileOptions), intent(inout)    :: options

        ! Default options
        options%runtype     = 'optimize'
        options%gridtype    = 'plasma'
        options%meth        = 'KKT'
        options%export      = .true.  

    end subroutine


end module