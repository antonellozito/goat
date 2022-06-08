!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains all the type definitions used in the grid 
! deformation module. See below for an explanation of each type. 

module gdmod_types

    ! Initialize
    !============
    ! The usual
    implicit none
    save
    private 

    ! All types
    public RunfileOptionsUDT
    public ExportOptionsUDT
    public GridUDT
    public DesignParamsUDT 
    

    ! All functions
    public SetDefaultRunfileOptions

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    ! Option types
    !=============
    ! Each option type has a setting routine called 
    ! SetDefault<optionname>, see subroutines after the 'contain' 
    ! statement. 

    ! Options for the main runfile type
    type RunfileOptionsUDT
        character(32)       :: runtype ! type of run: 'optimize' or 'test'
        character(32)       :: gridtype ! type of grid: 'plasma' 
        character(32)       :: meth ! method for grid deformation: 'KKT'
        logical             :: export ! do export? 
    end type

    ! Options for exporting
    type ExportOptionsUDT
        character(32)       :: gridformat ! format of grid to export
    end type


    ! Grid types
    !===========
    ! Main grid structure
    type GridUDT
        ! Coordinates
        real, allocatable   :: x, y ! grid coordinates

    end type

    ! Optimization
    !=============
    ! Design parameter structure
    type DesignParamsUDT
        ! Design variables
        real, allocatable   :: phi 

    end type

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    contains 

    subroutine SetDefaultRunfileOptions(options)
        ! Description
        !============
        ! Set the default runfile options

        ! Declaration
        type (RunfileOptionsUDT), intent(inout)    :: options

        ! Default options
        options%runtype     = 'optimize'
        options%gridtype    = 'plasma'
        options%meth        = 'KKT'
        options%export      = .true.  

    end subroutine


end module