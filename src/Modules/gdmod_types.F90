!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains all the type definitions used in the grid 
! deformation module. See the option setting routines in gdmod_userinput
! for an explanation of the different fields and options. 

module gdmod_types

    ! Initialize
    !============
    ! The usual
    implicit none
    save
    private 

    ! All public types
    public RunfileOptionsUDT
    public ExportOptionsUDT
    public DesignOptionsUDT
    public GridOptionsUDT
    public NumOptionsUDT
    public MagneticFieldOptionsUDT
    public GridUDT
    public DesignParamsUDT 
    
    ! All public functions
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

    ! Options for design variables
    type DesignVariableOptionsUDT
        character(32)       :: type ! type of design variables
    end type

    ! Options for the cost function
    type CostFunctionOptionsUDT
        character(32)       :: type ! the cost function format
    end type

    ! Options for the constraints
    type ConstraintOptionsUDT
        ! Fields for equality constraints
        logical             :: fluxfunction ! impose constraints on flux
        logical             :: xpoints ! impose x-point location
        logical             :: edgelengths ! impose edge length cons
        logical             :: orthogonality 

        ! Fields for inequality constraints
        logical             :: linefolding ! prevent flux line folding
    end type

    ! Options for design optimization
    type DesignOptionsUDT
        type(CostFunctionOptionsUDT)       :: costfunction
        type(DesignVariableOptionsUDT)     :: variables
        type(ConstraintOptionsUDT)         :: constraints
    end type

    ! Options for the grid
    type GridOptionsUDT
        logical             :: vesselrefine
        real(kind=8)        :: vesselmaxdist
    end type  

    ! Options for numerics
    type NumOptionsUDT
        integer(kind=8)     :: itmax 
    end type

    ! Options for magnetic field
    type MagneticFieldOptionsUDT
        character(32)       :: readmeth
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