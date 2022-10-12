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

    contains 

        ! Design initialization
        procedure(InitializeINT), deferred  :: Initialize

        ! Design update
        !procedure(UpdateINT), deferred      :: Update

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

        ! Housekeeping
        procedure :: AllocateDesignCoordinates
        procedure :: DeallocateDesignCoordinates
        final :: DestroyDesignCoordinates

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
        ! Set the number of design variables
        designvariables%nphi = 2*grid%vert%ntot 

        ! Allocate (specific routine for this design variable type)
        call designvariables%AllocateDesignCoordinates(grid)

        print *, designvariables%nphi

        ! Set the design variables
        designvariables%phi(1:grid%vert%ntot) = grid%vert%x 
        designvariables%phi(grid%vert%ntot+1:designvariables%nphi) = &
            grid%vert%y 

        ! Set other fields
        designvariables%xind = (/(i, i = 1, grid%vert%ntot)/)
        designvariables%yind = designvariables%xind + grid%vert%ntot

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

end module