!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains all the type definitions used in the grid 
! deformation module. See the option setting routines in gdmod_userinput
! for an explanation of the different fields and options. Additionally,
! precision types are set here that can (should) be used in all 
! subroutines. 

! Note: default accessibility is set to public, since all these types 
! and corresponding functions are free to be used. 

module gdmod_types

    ! Initialize
    !============
    ! The usual
    implicit none
    save
    public 

    ! Precision types
    integer, parameter       :: R4      = selected_real_kind(6)
    integer, parameter       :: R8      = selected_real_kind(14)
    integer, parameter       :: R16     = selected_real_kind(33)

    integer, parameter       :: C32     = 32
    integer, parameter       :: C64     = 64
    integer, parameter       :: C132    = 132

    integer, parameter       :: I4      = selected_int_kind(4)
    integer, parameter       :: I8      = selected_int_kind(8)
    integer, parameter       :: I16     = selected_int_kind(16)

    ! All private types
    
    ! All private functions

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
        character(C32)          :: runtype ! type of run: 'optimize' or 'test'
        character(C32)          :: gridtype ! type of grid: 'plasma' 
        character(C32)          :: meth ! method for grid deformation: 'KKT'
        logical                 :: export ! do export? 
    end type

    ! Options for exporting
    type ExportOptionsUDT
        character(C32)       :: gridformat ! format of grid to export
    end type

    ! Options for design variables
    type DesignVariableOptionsUDT
        character(C32)       :: type ! type of design variables
    end type

    ! Options for the cost function
    type CostFunctionOptionsUDT
        character(C32)       :: type ! the cost function format
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
        character(C32)          :: inputtype
        logical                 :: vesselrefine
        real(R8)                :: vesselmaxdist
    end type  

    ! Options for numerics
    type NumOptionsUDT
        integer(R8)     :: itmax 
    end type

    ! Options for magnetic field
    type MagneticFieldOptionsUDT
        character(C32)       :: readmeth
    end type


    ! Grid types
    !===========
    ! Vertex structure
    type VertexUDT
        ! Coordinates
        real(R8), allocatable               :: x(:),y(:) 

        ! Logicals and indices
        logical, allocatable                :: BV
        integer(I8), allocatable            :: fieldlineID 
        integer(I8)                         :: ntot

    end type

    ! Face structure
    type FaceUDT
        ! Logicals and indices
        integer(I8), allocatable            :: vert
        logical, allocatable                :: V 
    end type

    ! Cell structure
    type CellUDT
        ! Logicals and indices
        integer(I8), allocatable            :: vert
        integer(I8)                         :: ntot
    end type

    ! Main grid structure
    type GridUDT
        ! Vertices
        type(VertexUDT)                     :: vert

        ! Faces
        type(FaceUDT)                       :: faces

        ! Cells
        type(CellUDT)                       :: cells
    end type

    ! Optimization
    !=============
    ! Design parameter structure
    type DesignParamsUDT
        ! Design variables
        real, allocatable   :: phi(:) 

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

    ! Grid routines
    !==============
    ! Allocation
    subroutine AllocateGrid(grid)
        ! Description
        !============
        ! Allocate the fields in the vertex, faces, and cell structures
        ! of the grid. 

        ! The usual
        implicit none

        ! Declare variables
        type(GridUDT)       :: grid

        ! Return if already allocated
        if (allocated(grid%vert%x)) return

        ! Allocate vertex data
        allocate(grid%vert%x(grid%vert%ntot))     


    end subroutine

    ! Deallocation


end module