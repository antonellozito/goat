!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains the constraint classes specific for the grid 
! deformation. It relies on user input defined in the gdmod_userinput
! module, and on the gdmod_user types. 

! The constraints are structured as follows:
! - All derived constraint types inherit from the 'mother' type 
!   'GenericConstraintsGDUDT', which contains the field 'ncon', the 
!   number of constraints, and the initialization and evaluation 
!   routines. 
! - The overarching constraints structure contains an object 'eqcon' and 
!   'ineqcon', which are objects that contain the specific equality and
!   inequality constraints, respectively. Both objects have general
!   initialization, evaluation, and destruction routines that should be
!   used in the optimizer. 
! - The equality constraints contain different fields (e.g. fluxfunction
!   ) with a logical (e.g. dofluxfunction) that indicates whether the
!   constraint should be considered. 
! - Each specific constraints (e.g. fluxfunction) has its own evaluation
!   , initialization, and destruction routines (type-bound). 

module gdmod_constraints
    
    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use optmod_constraints
    use gdmod_types

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
    ! Generic constraint type
    type, abstract, extends(ConstraintsUDT) :: GenericConstraintsGDUDT

        ! Description
        !============
        ! Generic type for grid deformation constraints. Inherits from 
        ! the mother constraint type defined in gdmod_constraints.

        ! Fields:
        integer(I8)                 :: ncon = 0 ! number of constraints

    contains 

        ! Initialization
        procedure(InitializeConstraintsINT), deferred :: Initialize 

    end type

    ! Specific constraint types
    !==========================
    ! Flux function constraints
    type, extends(GenericConstraintsGDUDT) :: FluxfunctionConstraintsUDT

        ! Description
        !============
        ! Flux function constraints

        ! Fields:
        



    contains

        ! Initialization
        procedure :: Initialize     => InitializeFluxfunctionConstraints
        ! procedure :: Allocate       => AllocateFluxfunctionConstraints

        ! Destructor
        ! procedure :: DestroyFluxfunctionConstraints

    end type

    ! Boundary function constraints
    type, extends(ConstraintsUDT) :: BoundaryFunctionConstraintsUDT

    end type

    ! Overarching types
    !==================
    ! Equality constraints
    type, extends(ConstraintsUDT) :: EqConGDUDT

        ! Description
        !============
        ! This type contains all the different constraints as different
        ! derived types. For each type of constraint, a different
        ! type is defined. 

        ! Total number of constraints 
        integer(I8)                             :: neqcon = 0

        ! Flux function constraint
        logical                                 :: dofluxfunction = .false.
        type(FluxfunctionConstraintsUDT)        :: fluxfunction 
        

    contains

        ! Procedure to initialize constraints
        procedure :: Initialize         => InitializeEqCon

        ! Procedure to evaluate constraints
        !procedure :: Evaluate           => EvaluateEqCon

    end type 

    ! Inequality constraints
    type, extends(ConstraintsUDT) :: IneqConGDUDT

        ! Total number of inequality constraints
        integer(I8)             :: nineqcon = 0

    end type

    ! All constraints for the grid deformation
    type, extends(ConstraintsUDT) :: ConstraintsGDUDT

        ! Description
        !============
        ! Defines the basic optimization problem: it has a set of 
        ! design variables, constraints, and a cost function. 

        ! Fields: 

        ! Equality constraints 
        type(EqConGDUDT)        :: eqcon

        ! Inequality constraints 
        type(IneqConGDUDT)      :: ineqcon

    contains

        ! Initialization
        procedure :: Initialize                 => InitializeConstraints

        ! Number of constraints getter
        procedure :: GetConstraintsDimensions  

        ! Evaluation
        ! procedure(EvaluateConstraintsINT), deferred :: Evaluate
        
        ! Housekeeping

    end type

    ! Monitor
    !========
    type ConstraintsMonitorUDT

        ! Description
        !============
        ! This object can be used to keep track of the amount of 
        ! constraints that are imposed per vertex. As such, it can be
        ! checked whether the problem will be overly constrained and to 
        ! prevent constraints to be imposed on certain vertices if 
        ! a maximum number is reached. For now, the equality constraints
        ! per vertex are counted in eqvcc, the inequality in ineqvcc. 
        ! The maximum number of equality and inequality constraints per 
        ! vertex are set in maxeqvcc and maxineqvcc. The methods are 
        ! the following:

        ! - Initialize: initialization routine that takes the grid, 
        !   magnetic field, and environment as input (in case of 
        !   extension of this object in the future)
        ! - CheckConstraints: will check if the problem is already 
        !   overly constrained, or may become overly constrained if 
        !   inequality constraints become active. 

        ! Fields:
        integer(I8), allocatable           :: eqvcc(:), ineqvcc(:)
        integer(I8)                        :: maxeqvcc, maxineqvcc 
        
    contains 

        ! Routines
        procedure :: Initialize         => InitializeMonitor
        ! procedure :: CheckConstraints   => CheckConstraintsMonitor

        ! Finalization
        final :: DestroyMonitor

    end type

    

    ! Derived types
    !==============
    

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    ! Abstract interfaces
    !====================
    abstract interface

        ! Constraint initialization
        subroutine InitializeConstraintsINT(constraints, grid, &
            magneticField, environment, monitor)

            ! Description
            !============
            ! This routine serves as general initialization routine for
            ! a generic grid deformation constraint (that inherits from
            ! the GenericConstraintUDT type)

            ! Import
            import :: GenericConstraintsGDUDT, GridUDT, &
                MagneticFieldUDT, EnvironmentUDT, ConstraintsMonitorUDT

            ! Declare
            class(GenericConstraintsGDUDT)      :: constraints 
            type(GridUDT)                       :: grid 
            type(MagneticFieldUDT)              :: magneticField 
            type(EnvironmentUDT)                :: environment 
            type(ConstraintsMonitorUDT)         :: monitor

        end subroutine

    end interface

    contains

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                           CONSTRAINTS MONITOR                    !
    !------------------------------------------------------------------!
    ! Initialization
    subroutine InitializeMonitor(monitor, grid, magneticField, &
        environment)

        ! Description
        !============
        ! Initializes the constraints monitor structure. It is assumed
        ! that the grid, magnetic field and environment are properly 
        ! allocated and initialized. 

        ! Declare variables
        !==================
        ! Arguments
        class(ConstraintsMonitorUDT)        :: monitor 
        type(GridUDT)                       :: grid 
        type(MagneticFieldUDT)              :: magneticField 
        type(EnvironmentUDT)                :: environment

        ! Loop variables

        ! Auxiliary 

        ! Initialize
        !===========
        ! Set the maximal number of constraints
        monitor%maxeqvcc = 2 ! for most cases this is fine
        monitor%maxineqvcc = 1000 ! a stupid large number - can impose any number

        ! Allocate
        allocate(monitor%eqvcc(grid%vert%ntot))
        allocate(monitor%ineqvcc(grid%vert%ntot))

        ! Initialize
        monitor%eqvcc(:)    = 0
        monitor%ineqvcc(:)  = 0

    end subroutine

    ! Destructor
    subroutine DestroyMonitor(monitor)

        ! Description
        !============
        ! Finalization for the constraints monitor object

        ! Declare variables
        !==================
        ! Arguments
        type(ConstraintsMonitorUDT)         :: monitor 

        ! Destroy
        !========
        deallocate(monitor%eqvcc, monitor%ineqvcc)

    end subroutine

    !------------------------------------------------------------------!
    !                           GENERAL CONSTRAINTS                    !
    !------------------------------------------------------------------!
    ! Initialization
    subroutine InitializeConstraints(constraints, grid, magneticField, &
        environment, constraintoptions)

        ! Description
        !============
        ! Routine that initializes the equality and inequality 
        ! constraints, using the initialization routines of those 
        ! objects. 

        ! Declare variables
        !==================
        ! Arguments
        class(ConstraintsGDUDT)     :: constraints 
        type(GridUDT)               :: grid 
        type(MagneticFieldUDT)      :: magneticField 
        type(EnvironmentUDT)        :: environment 
        type(ConstraintOptionsUDT)  :: constraintoptions
        type(ConstraintsMonitorUDT) :: monitor

        ! Loop variables

        ! Auxiliary variables

        ! Initialize monitor
        !===================
        call monitor%Initialize(grid, magneticField, environment)

        ! Initialize constraints
        !=======================
        ! Equality constraints
        call constraints%eqcon%Initialize(grid, magneticField, &
            environment, constraintoptions, monitor)

        ! Inequality constraints (to do)

    end subroutine

    ! Dimension getter
    subroutine GetConstraintsDimensions(constraints, neqcon, nineqcon)

        ! Description
        !============
        ! Return the current dimensions of the equality and inequality
        ! constraints. Can be used for initialization of other 
        ! quantities at higher levels. 

        ! Declare variables
        !==================
        ! Arguments
        class(ConstraintsGDUDT)         :: constraints 
        integer(I8), intent(out)        :: neqcon, nineqcon

        ! Extract dimensions
        !===================
        ! Stored in eqcon, ineqcon
        neqcon = constraints%eqcon%neqcon 
        nineqcon = constraints%ineqcon%nineqcon

    end subroutine

    !------------------------------------------------------------------!
    !                           EQUALITY CONSTRAINTS                   !
    !------------------------------------------------------------------!
    ! Initialization
    subroutine InitializeEqCon(constraints, grid, magneticField, &
        environment, constraintoptions, monitor)

        ! Description
        !============
        ! Routine that initializes the desired constraints, based on the
        ! user options d

        ! Declare variables
        !==================
        ! Arguments
        class(EqConGDUDT)           :: constraints 
        type(GridUDT)               :: grid 
        type(MagneticFieldUDT)      :: magneticField 
        type(EnvironmentUDT)        :: environment 
        type(ConstraintOptionsUDT)  :: constraintoptions
        type(ConstraintsMonitorUDT) :: monitor

        ! Loop variables

        ! Auxiliary variables

        ! Initialize constraints
        !=======================
        ! Flux function
        if (constraintoptions%fluxfunction) then 
            ! Set the logical
            constraints%dofluxfunction = .true.

            ! Initialize
            call constraints%fluxfunction%Initialize(grid, &
                magneticField, environment, monitor)

        else
            ! Set to false, don't initialize
            constraints%dofluxfunction = .false.

        end if



    end subroutine

    !------------------------------------------------------------------!
    !                           FLUX FUNCTION                          !
    !------------------------------------------------------------------!
    ! Initialize
    subroutine InitializeFluxfunctionConstraints(constraints, grid, &
        magneticField, environment, monitor)

        ! Description
        !============
        ! Initialize the required fields of the flux function 
        ! constraints. These constraints impose for each vertex that 
        ! lies on a flux surface (i.e. has a flux surface ID which is
        ! non-zero) by fixing its flux value. The flux values can be
        ! determined in different ways. Here, we compute the initial
        ! psi values by averaging the current psi values of the vertices
        ! that lie on a flux surface. This hedges a bit for 
        ! discretization errors originating from the re-interpretation 
        ! of the magnetic field data as a bicubic spline interpolant 
        ! (here) instead of a linear interpolant in the grid generator. 

        ! Notes
        !======
        ! Note 1: the last argument of this function is a derived type 
        ! used to monitor the constraints. It can be used to make sure
        ! that the problem is - at least not by the equality constraints 
        ! - is overly constrained. It is UP TO THE DEVELOPER to use this
        ! monitor properly!
        
        ! Declare variables
        !==================
        ! Arguments 
        class(FluxfunctionConstraintsUDT)       :: constraints 
        type(GridUDT)                           :: grid 
        type(MagneticFieldUDT)                  :: magneticField 
        type(EnvironmentUDT)                    :: environment 
        type(ConstraintsMonitorUDT)             :: monitor

        ! Loop variables

        ! Auxiliary variables

        ! Data


        ! Initialize
        !===========
        constraints%ncon = 10

        ! Determine flux values to impose
        !================================

        ! Determine constraints
        

    end subroutine

end module