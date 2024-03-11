!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains all the necessary routines to manipulate the 
! design variables (except for the actual optimization routines). The
! design variable UDT should be set in gdmod_types, user input in 
! gdmod_userinput. 

! Notes
!======
! Note 1: this module is specifically built for grid optimization. Most
! routines therefore expect at least the input arguments 'grid', 
! 'magneticField', and 'environment'. 

! Note 2: although the code structure (should) allow(s) for any type of 
! design variable, currently only the 'coordinates' type is implemented.

module gdmod_designvariables
    
    ! Initialize
    !============
    ! Load modules
    use gdmod_types
    use gdmod_userinput 
    use optmod_designvariables

    ! The usual
    implicit none
    save
    public 

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    ! Abstract types
    !===============
    ! Main design variable type for the grid deformation. Others should 
    ! be derived from this type. 
    type, abstract, extends(DesignVariablesUDT) :: DesignVariablesGDUDT

        character(:), allocatable                :: type 

    contains 

        ! Design initialization
        procedure(InitializeINT), deferred      :: Initialize

        ! Design update
        procedure(UpdateDesignINT), deferred      :: UpdateDesign

    end type

    ! Derived types
    !==============
    ! Should all extend the base type DesignVariablesUDT, defined in 
    ! optmod_designvariables.

    ! Design variables for coordinates
    type, extends(DesignVariablesGDUDT) :: DesignVariablesCoordinatesUDT 

        ! Coordinate indices
        integer(I8), allocatable        :: xind(:), yind(:) 

    contains 
    
        ! Design initialization
        procedure :: Initialize     => InitializeDesignCoordinates

        ! Design update
        procedure :: UpdateDesign   => UpdateDesignCoordinates

        ! Housekeeping
        procedure :: AllocateDesignCoordinates
        procedure :: DeallocateDesignCoordinates
        final :: DestroyDesignCoordinates

    end type

    ! Design variables for flux values
    type, extends(DesignVariablesGDUDT) :: DesignVariablesFluxValuesUDT

        ! Flux value indices
        integer(I8), allocatable        :: psiind(:)

    contains

        ! Design initialization
        procedure :: Initialize     => InitializeDesignFluxValues

        ! Design update
        procedure :: UpdateDesign   => UpdateDesignFluxValues

        ! Housekeeping
        procedure :: AllocateDesignFluxValues
        procedure :: DeallocateDesignFluxValues
        final :: DestroyDesignFluxValues

    end type

    ! Design variables, combined coordinates and flux
    type, extends(DesignVariablesGDUDT) :: DesignVariablesCoordinatesFluxUDT

        ! Coordinate & flux indices
        integer(I8), allocatable        :: xind(:), yind(:), psiind(:) 

    contains

        ! Design initialization
        procedure :: Initialize     => InitializeDesignCoordinatesFlux

        ! Design update
        procedure :: UpdateDesign   => UpdateDesignCoordinatesFlux

        ! Housekeeping
        procedure :: DeallocateDesignCoordinatesFlux
        final :: DestroyDesignCoordinatesFlux
    
    end type

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    ! Abstract interfaces
    !====================
    abstract interface

        ! Design initialization, specific for our purposes
        subroutine InitializeINT(designvariables, grid, magneticField, &
            environment)

            ! Description
            !============
            ! Interface for the design initialization routine that each
            ! derived type should have. Since this is specific for the
            ! grid deformation routines, the grid, magnetic field and
            ! environment structures can/have to be passed. 

            ! Import
            import :: DesignVariablesGDUDT, GridUDT, MagneticFieldUDT, & 
                EnvironmentUDT

            ! Declare
            class(DesignVariablesGDUDT)         :: designvariables 
            type(GridUDT), intent(in)           :: grid
            type(MagneticFieldUDT), intent(in)  :: magneticField 
            type(EnvironmentUDT) , intent(in)   :: environment

        end subroutine

        ! Design update
        subroutine UpdateDesignINT(designvariables, grid, &
            magneticfield, environment)

            ! Description
            !============
            ! Interface for the design update routine that each
            ! derived type should have. Since this is specific for the
            ! grid deformation routines, the grid, magnetic field and
            ! environment structures can/have to be passed. 
            ! This routine should update the grid, magnetic field, and
            ! environment structure according to the values given in 
            ! designvariables%phi

            ! Import
            import :: DesignVariablesGDUDT, GridUDT, MagneticFieldUDT, & 
                EnvironmentUDT

            ! Declare
            class(DesignVariablesGDUDT)             :: designvariables 
            type(GridUDT), intent(inout)            :: grid
            type(MagneticFieldUDT), intent(inout)   :: magneticField 
            type(EnvironmentUDT) , intent(inout)    :: environment

        end subroutine

    end interface

    contains

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                              Coordinates                         !
    !------------------------------------------------------------------!

    ! Design initialization
    subroutine InitializeDesignCoordinates(designvariables, grid, &
        magneticField, environment)

        ! Description
        !============
        ! Initialize the design variables, which are in this case only 
        ! the x and y coordinates. 

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariablesCoordinatesUDT)        :: designvariables
        type(gridUDT), intent(in)                   :: grid 
        type(MagneticFieldUDT), intent(in)          :: magneticField 
        type(EnvironmentUDT), intent(in)            :: environment
    
        ! Loop variables
        integer(I8)                                 :: i

        ! Auxiliary variables 

        ! Data

        ! Initialize
        !===========
        ! Set the type
        designvariables%type = 'coordinates' 
        
        ! Set the number of design variables
        designvariables%nphi = 2*grid%vert%ntot 

        ! Allocate (specific routine for this design variable type)
        call designvariables%AllocateDesignCoordinates(grid)

        ! Set the design variables
        designvariables%phi(1:grid%vert%ntot) = grid%vert%x 
        designvariables%phi(grid%vert%ntot+1:designvariables%nphi) = &
            grid%vert%y 

        ! Set other fields
        designvariables%xind = (/(i, i = 1, grid%vert%ntot)/)
        designvariables%yind = designvariables%xind + grid%vert%ntot

    end subroutine

    ! Design update
    subroutine UpdateDesignCoordinates(designvariables, grid, &
        magneticField, environment)

        ! Description
        !============
        ! Update the design coordinates according to the phi values. 

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariablesCoordinatesUDT)        :: designvariables
        type(gridUDT), intent(inout)                :: grid 
        type(MagneticFieldUDT), intent(inout)       :: magneticField 
        type(EnvironmentUDT), intent(inout)         :: environment
    
        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Update
        !=======
        ! Grid coordinates
        grid%vert%x = designvariables%phi(designvariables%xind)
        grid%vert%y = designvariables%phi(designvariables%yind)

    end subroutine

    ! Housekeeping
    subroutine AllocateDesignCoordinates(designvariables, grid)

        ! Description
        !============
        ! Allocate. Assumed that nphi is given. 

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariablesCoordinatesUDT)        :: designvariables 
        type(GridUDT)                               :: grid

        ! Allocate
        !=========
        ! Call the parent allocation routine
        call designvariables%Allocate()

        ! Allocate the remaining fields, specific for this type
        allocate(designvariables%xind(grid%vert%ntot))
        allocate(designvariables%yind(grid%vert%ntot))

    end subroutine

    subroutine DeallocateDesignCoordinates(designvariables)

        ! Description
        !============
        ! Deallocate

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariablesCoordinatesUDT)        :: designvariables 

        ! Allocate
        !=========
        ! Call the parent allocation routine
        call designvariables%Deallocate()

        ! Allocate the remaining fields, specific for this type
        deallocate(designvariables%xind)
        deallocate(designvariables%yind)

    end subroutine

    subroutine DestroyDesignCoordinates(designvariables)

        ! Description
        !============
        ! Destructor
        
        ! Declare variables
        !==================
        ! Arguments
        type(DesignVariablesCoordinatesUDT) :: designvariables

        ! Destroy
        !========
        call designvariables%Deallocate()
        
    end subroutine

    !------------------------------------------------------------------!
    !                              Flux values                         !
    !------------------------------------------------------------------!

    ! Design initialization
    subroutine InitializeDesignFluxValues(designvariables, grid, &
        magneticField, environment)

        ! Description
        !============
        ! Initialize the design variables, which are in this case only 
        ! the x and y coordinates. 

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariablesFluxValuesUDT)        :: designvariables
        type(gridUDT), intent(in)                   :: grid 
        type(MagneticFieldUDT), intent(in)          :: magneticField 
        type(EnvironmentUDT), intent(in)            :: environment
    
        ! Loop variables
        integer(I8)                                 :: i

        ! Auxiliary variables 

        ! Data

        ! Initialize
        !===========
        ! Set the type
        designvariables%type = 'fluxvalues' 

        ! Further initialization can only be done after initializing
        ! the flux value constraints - see implementation under 
        ! FinalizeInitialization routine of GD optimization problem type

    end subroutine

    ! Design update
    subroutine UpdateDesignFluxValues(designvariables, grid, &
        magneticField, environment)

        ! Description
        !============
        ! Void routine - no grid, magnetic field or environment 
        ! parameters to be updated based on design variables

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariablesFluxValuesUDT)        :: designvariables
        type(gridUDT), intent(inout)                :: grid 
        type(MagneticFieldUDT), intent(inout)       :: magneticField 
        type(EnvironmentUDT), intent(inout)         :: environment

    end subroutine

    ! Housekeeping
    subroutine AllocateDesignFluxValues(designvariables)

        ! Description
        !============
        ! Allocate. Assumed that nphi is given. 

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariablesFluxValuesUDT)        :: designvariables 

        ! Allocate
        !=========
        ! Call the parent allocation routine
        call designvariables%Allocate()

        ! Allocate others
        allocate(designvariables%psiind(designvariables%nphi))

    end subroutine

    subroutine DeallocateDesignFluxValues(designvariables)

        ! Description
        !============
        ! Deallocate

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariablesFluxValuesUDT)        :: designvariables 

        ! Allocate
        !=========
        ! Call the parent allocation routine
        call designvariables%Deallocate()

        ! Allocate the remaining fields, specific for this type
        deallocate(designvariables%psiind)

    end subroutine

    subroutine DestroyDesignFluxValues(designvariables)

        ! Description
        !============
        ! Destructor
        
        ! Declare variables
        !==================
        ! Arguments
        type(DesignVariablesFluxValuesUDT) :: designvariables

        ! Destroy
        !========
        call designvariables%Deallocate()
        
    end subroutine

    !------------------------------------------------------------------!
    !                       Coordinates & Flux values                  !
    !------------------------------------------------------------------!

    ! Design initialization
    subroutine InitializeDesignCoordinatesFlux(designvariables, grid, &
        magneticField, environment)

        ! Description
        !============
        ! Initialize the design variables, which are in this case only 
        ! the x and y coordinates. 

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariablesCoordinatesFluxUDT)    :: designvariables
        type(gridUDT), intent(in)                   :: grid 
        type(MagneticFieldUDT), intent(in)          :: magneticField 
        type(EnvironmentUDT), intent(in)            :: environment
    
        ! Loop variables
        integer(I8)                                 :: i

        ! Auxiliary variables 

        ! Data

        ! Initialize
        !===========
        ! Set the type
        designvariables%type = 'coordinates_desiredflux' 

        ! All other initialization is deferred until after the problem
        ! initialization, since the amount of design variables etc 
        ! depends on the number of flux surface constraints
        

    end subroutine

    ! Design update
    subroutine UpdateDesignCoordinatesFlux(designvariables, grid, &
        magneticField, environment)

        ! Description
        !============
        ! Update the design coordinates according to the phi values. 

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariablesCoordinatesFluxUDT)    :: designvariables
        type(gridUDT), intent(inout)                :: grid 
        type(MagneticFieldUDT), intent(inout)       :: magneticField 
        type(EnvironmentUDT), intent(inout)         :: environment
    
        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Update
        !=======
        ! Grid coordinates
        grid%vert%x = designvariables%phi(designvariables%xind)
        grid%vert%y = designvariables%phi(designvariables%yind)
        
        ! Flux values should be updated in constraints...

    end subroutine

    subroutine DeallocateDesignCoordinatesFlux(designvariables)

        ! Description
        !============
        ! Deallocate

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariablesCoordinatesFluxUDT)        :: designvariables 

        ! Allocate
        !=========
        ! Call the parent allocation routine
        call designvariables%Deallocate()

        ! Allocate the remaining fields, specific for this type
        deallocate(designvariables%xind)
        deallocate(designvariables%yind)
        deallocate(designvariables%psiind)

    end subroutine

    subroutine DestroyDesignCoordinatesFlux(designvariables)

        ! Description
        !============
        ! Destructor
        
        ! Declare variables
        !==================
        ! Arguments
        type(DesignVariablesCoordinatesFluxUDT) :: designvariables

        ! Destroy
        !========
        call designvariables%Deallocate()
        
    end subroutine


end module