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
    use mod_sparseinterface
    use optmod_constraints
    use gdmod_types
    use gdmod_designvariables
    use gdmod_utility_optimization
    use gdmod_plots
    use PolygonShapeFunction
    use BicubicSplineInterpolant
    

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

        ! Evaluation
        procedure(EvaluateConstraintsINT), deferred :: Evaluate

    end type

    ! Specific constraint types
    !==========================
    ! Flux function constraints
    type, extends(GenericConstraintsGDUDT) :: FluxfunctionConstraintsUDT

        ! Description
        !============
        ! Flux function constraints. Fixes the flux function values 
        ! for a set of desired nodes. The following fields are added:
        ! - vert:       the vertices to consider (1 constraint per 
        !               entry of this array, so ncon-by-1 dimension)
        ! - PsiD:       Desired psi value of each vertex (ncon-by-1)
        ! - ncon:       (inherited) number of constraints 

        ! No other routines than the standard initialization, evaluation
        ! and destruction routines are implemented nor needed. 

        ! Fields: 
        integer(I8), allocatable        :: vert(:)
        real(R8), allocatable           :: PsiD(:)

    contains

        ! Initialization
        procedure :: Initialize     => InitializeFluxfunctionConstraints

        ! Evaluation
        procedure :: Evaluate       => EvaluateFluxfunctionConstraints

        ! Destructor
        final :: DestroyFluxfunctionConstraints

    end type

    ! Boundary function constraints
    type, extends(GenericConstraintsGDUDT) :: BoundaryFunctionConstraintsUDT 

        ! Description
        !============
        ! Boundary function constraints. The boundary function is given
        ! by the 'psf' polygon (here, only one polygon can be specified
        ! though this may be extended in the future). The vertices on
        ! which these conditions have to be imposed are given in the 
        ! 'vertID' array (number of vertices given in nvertID)

        ! Fields: 
        type(PolygonShapeFunctionUDT)       :: psf 
        integer(I8), allocatable            :: vert(:)
        integer(I8)                         :: nvert 

    contains

        ! Initialization
        procedure :: Initialize => InitializeBoundaryFunctionConstraints
        
        ! Evaluation
        procedure :: Evaluate   => EvaluateBoundaryFunctionConstraints
    
    end type

    ! X-point constraints
    type, extends(GenericConstraintsGDUDT) :: XPointConstraintsUDT 

        ! Description
        !============
        ! X-point constraints. Constrains the location of the x-point by 
        ! fixing the initial x-point coordinates (or, perhaps in the
        ! future, look at where the gradient of the flux function 
        ! vanishes). The following fields are added:
        !
        ! - xpind:      the x-point vertex IDs 
        ! - nxpind:     the total number of x-points    
        ! - locx/y:     x and y coordinate locations of the x-points
    
        ! Note 1: the initial x-point location and x-point vertices
        ! are determined by the initial grid. The x-point indices are 
        ! determined using the DetermineXPoints routine in 
        ! gdmod_utility_optimization. 

        ! Fields: 
        integer(I8), allocatable            :: xpind(:)
        integer(I8)                         :: nxpind 
        real(R8), allocatable               :: locx(:)
        real(R8), allocatable               :: locy(:)

    contains

        ! Initialization
        procedure :: Initialize => InitializeXPointConstraints
        
        ! Evaluation
        procedure :: Evaluate   => EvaluateXPointConstraints
    
    end type

    ! Edge lengths constraints
    type, extends(GenericConstraintsGDUDT) :: EdgeLengthsConstraintsUDT 

        ! Description
        !============
        ! Edge length constraints. Constrain an arbitrary set of edges
        ! to have a certain edge length (in meters) - see also 
        ! InitializeEdgeLengthsConstraints for additional predefined
        ! edge sets. The following fields are added:
        !
        ! - nedges:     (scalar) total number of edges to constrain
        ! - edgevert:   (nedges-by-2) edge vertex IDs 
        ! - d           (nedges-by-1) desired length for each edge   

        ! Fields: 
        integer(I8), allocatable            :: edgevert(:, :)
        integer(I8)                         :: nedges
        real(R8), allocatable               :: d(:)

    contains

        ! Initialization
        procedure :: Initialize => InitializeEdgeLengthsConstraints
        
        ! Evaluation
        procedure :: Evaluate   => EvaluateEdgeLengthsConstraints
    
    end type

    ! Orthogonality constraints
    type, extends(GenericConstraintsGDUDT) :: OrthogonalityConstraintsUDT 

        ! Description
        !============
        ! Orthogonality constraints. These constraints fix certain edges
        ! to be orthogonal to the magnetic field. See the initialization
        ! routine (InitializeOrthogonalityConstraints) to see the 
        ! different options that are available to determine which 
        ! edges to constrain. Typically, edges near the core and SOL
        ! are constrained. 
        ! The following fields are added:
        !
        ! - nedges:     (scalar) total number of edges to constrain
        ! - edgevert:   (nedges-by-2) edge vertex IDs 

        ! Fields: 
        integer(I8), allocatable            :: edgevert(:, :)
        integer(I8)                         :: nedges

    contains

        ! Initialization
        procedure :: Initialize => InitializeOrthogonalityConstraints
        
        ! Evaluation
        procedure :: Evaluate   => EvaluateOrthogonalityConstraints
    
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
        integer(I8)                         :: neqcon = 0

        ! Constraint switches
        logical                             :: dofluxfunction = .false.
        logical                             :: doboundaryfunction = .false.
        logical                             :: doxpoints = .false.
        logical                             :: doedgelengths = .false.
        logical                             :: doorthogonality = .false.

        type(FluxfunctionConstraintsUDT)    :: fluxfunction 
        type(BoundaryFunctionConstraintsUDT):: boundaryfunction
        type(XPointConstraintsUDT)          :: xpoints
        type(EdgeLengthsConstraintsUDT)     :: edgelengths
        type(OrthogonalityConstraintsUDT)   :: orthogonality

    contains

        ! Procedure to initialize constraints
        procedure :: Initialize         => InitializeEqCon

        ! Procedure to evaluate constraints
        procedure :: Evaluate           => EvaluateEqCon

    end type 

    ! Inequality constraints
    type, extends(ConstraintsUDT) :: IneqConGDUDT

        ! Total number of inequality constraints
        integer(I8)             :: nineqcon = 0

    contains

        ! Constraints initialization
        procedure :: Initialize         => InitializeIneqCon

        ! Constraints evaluation
        procedure :: Evaluate           => EvaluateIneqCon

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
        procedure :: Initialize         => InitializeConstraints

        ! Number of constraints getter
        procedure :: GetConstraintsDimensions  

        ! Evaluation
        ! procedure :: Evaluate           => EvaluateEqualityConstraintsGD
        
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

        ! Constraint evaluation
        subroutine EvaluateConstraintsINT(constraints, G, gradG, & 
            hessG, grid, magneticField, environment, &
            dogradient, dohessian, designvariables, lambda)

            ! Description
            !============
            ! This reoutine serves as a general evaluation routine for 
            ! a generic grid deformation constraint. 

            ! Import
            import :: GenericConstraintsGDUDT, MySparseUDT, GridUDT, &
                R8, MagneticFieldUDT, EnvironmentUDT, & 
                DesignVariablesGDUDT
            
            ! Declare
            class(GenericConstraintsGDUDT)  :: constraints 
            real(R8), allocatable           :: G(:), lambda(:)
            type(MySparseUDT)               :: hessG, gradG 
            type(GridUDT)                   :: grid 
            type(MagneticFieldUDT)          :: magneticField 
            type(EnvironmentUDT)            :: environment 
            logical                         :: dogradient, dohessian
            class(DesignVariablesGDUDT)     :: designvariables

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

        ! Inequality constraints
        call constraints%ineqcon%Initialize(grid, magneticField, &
            environment, constraintoptions, monitor)

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
        constraints%neqcon = 0

        ! X-points
        if (constraintoptions%xpoints == 1) then 

            ! Set the logical
            constraints%doxpoints = .true.

            ! Initialize
            call constraints%xpoints%Initialize(grid, &
                magneticField, environment, monitor)

            ! Add constraints number
            constraints%neqcon = constraints%neqcon + &
                constraints%xpoints%ncon 

            ! Print
            print *, 'number of x-point function constraints: ', &
                constraints%xpoints%ncon


        else
            ! Set to false, don't initialize
            constraints%doxpoints = .false.

        end if

        ! Boundary function
        if (constraintoptions%boundaryfunctions == 1) then 

            ! Set the logical
            constraints%doboundaryfunction = .true.

            ! Initialize
            call constraints%boundaryfunction%Initialize(grid, &
                magneticField, environment, monitor)

            ! Add constraints number
            constraints%neqcon = constraints%neqcon + &
                constraints%boundaryfunction%ncon 

            ! Print
            print *, 'number of boundary function constraints: ', &
                constraints%boundaryfunction%ncon


        else
            ! Set to false, don't initialize
            constraints%doboundaryfunction = .false.

        end if

        ! Flux function
        if (constraintoptions%fluxfunction == 1) then 
            ! Set the logical
            constraints%dofluxfunction = .true.

            ! Initialize
            call constraints%fluxfunction%Initialize(grid, &
                magneticField, environment, monitor)

            ! Add constraints number
            constraints%neqcon = constraints%neqcon + &
                constraints%fluxfunction%ncon 

            ! Print
            print *, 'number of flux function constraints: ', &
                constraints%fluxfunction%ncon

        else
            ! Set to false, don't initialize
            constraints%dofluxfunction = .false.

        end if

        ! Edge lengths
        if (constraintoptions%edgelengths == 1) then 
            ! Set the logical
            constraints%doedgelengths = .true.

            ! Initialize
            call constraints%edgelengths%Initialize(grid, &
                magneticField, environment, monitor)

            ! Add constraints number
            constraints%neqcon = constraints%neqcon + &
                constraints%edgelengths%ncon 

            ! Print
            print *, 'number of edge lengths constraints: ', &
                constraints%edgelengths%ncon

        else
            ! Set to false, don't initialize
            constraints%doedgelengths = .false.

        end if

        ! Orthogonality
        if (constraintoptions%orthogonality == 1) then 
            ! Set the logical
            constraints%doorthogonality = .true.

            ! Initialize
            call constraints%orthogonality%Initialize(grid, &
                magneticField, environment, monitor)

            ! Add constraints number
            constraints%neqcon = constraints%neqcon + &
                constraints%orthogonality%ncon 

            ! Print
            print *, 'number of orthogonality constraints: ', &
                constraints%orthogonality%ncon

        else
            ! Set to false, don't initialize
            constraints%doorthogonality = .false.

        end if

        

    end subroutine

    ! Constraint evaluation
    subroutine EvaluateEqCon(constraints, G, gradG, hessG, &
        grid, magneticField, environment, dogradient, dohessian, & 
        designvariables, lambda)

        ! Description
        !============
        ! This routine evaluates the constraints G and the corresponding
        ! gradient and hessian. To do so, every type of constraint is 
        ! checked whether it is imposed, and the contributions are 
        ! added by calling the evaluation routine of each constraint. 

        ! Declare variables
        !==================
        ! Arguments
        class(EqConGDUDT)               :: constraints
        real(R8), intent(inout)         :: G(:)
        real(R8), intent(in)            :: lambda(:)
        type(MySparseUDT)               :: gradG, hessG 
        type(GridUDT)                   :: grid
        type(MagneticFieldUDT)          :: magneticField 
        type(EnvironmentUDT)            :: environment
        logical                         :: dogradient, dohessian 
        class(DesignVariablesGDUDT)     :: designvariables 

        ! Loop
        integer(I8)                     :: ic, ivg, ivh, k
        integer(I8), allocatable        :: conindex(:)

        ! Auxiliary
        real(R8), allocatable           :: G_flux(:), lambda_flux(:)
        type(MySparseUDT)               :: gradG_flux, hessG_flux

        real(R8), allocatable           :: G_bnd(:), lambda_bnd(:)
        type(MySparseUDT)               :: gradG_bnd, hessG_bnd

        real(R8), allocatable           :: G_xp(:), lambda_xp(:)
        type(MySparseUDT)               :: gradG_xp, hessG_xp

        real(R8), allocatable           :: G_el(:), lambda_el(:)
        type(MySparseUDT)               :: gradG_el, hessG_el

        real(R8), allocatable           :: G_orth(:), lambda_orth(:)
        type(MySparseUDT)               :: gradG_orth, hessG_orth

        ! Initialize
        !===========
        ! Set the constraint counter
        ic = 0

        ! Constraint values
        !==================

        ! X-point constraints
        !--------------------
        if (constraints%doxpoints) then 
            ! Construct the constraint index
            allocate(conindex(constraints%xpoints%ncon))
            conindex = [(k, k = ic+1, ic+constraints%xpoints%ncon)]

            ! Allocate & initialize
            allocate(lambda_xp(constraints%xpoints%ncon))
            lambda_xp = lambda(conindex)

            ! Call the evaluation routine
            call constraints%xpoints%Evaluate(G_xp, &
                gradG_xp, hessG_xp, &
                grid, magneticField, environment, dogradient, &
                dohessian, designvariables, &
                lambda_xp)

            ! Assign
            G(conindex) = G_xp

            ! Update the gradient column indices
            if (dogradient) then 
                ! For easier concatenation later on
                gradG_xp%col = gradG_xp%col + ic

            end if

            ! Update the constraint counter
            ic = ic + constraints%xpoints%ncon

            ! Housekeeping
            deallocate(conindex, lambda_xp) 
            if (allocated(G_xp)) then
                deallocate(G_xp)
            end if

        end if

        ! Boundary function constraints
        !------------------------------
        if (constraints%doboundaryfunction) then 
            ! Construct the constraint index
            allocate(conindex(constraints%boundaryfunction%ncon))
            conindex = [(k, k = ic+1, ic+constraints%boundaryfunction%ncon)]

            ! Allocate & initialize
            allocate(lambda_bnd(constraints%boundaryfunction%ncon))
            lambda_bnd = lambda(conindex)

            ! Call the evaluation routine
            call constraints%boundaryfunction%Evaluate(G_bnd, &
                gradG_bnd, hessG_bnd, &
                grid, magneticField, environment, dogradient, &
                dohessian, designvariables, &
                lambda_bnd)

            ! Assign
            G(conindex) = G_bnd

            ! Update the gradient column indices
            if (dogradient) then 
                ! For easier concatenation later on
                gradG_bnd%col = gradG_bnd%col + ic

            end if

            ! Update the constraint counter
            ic = ic + constraints%boundaryfunction%ncon

            ! Housekeeping
            deallocate(conindex, lambda_bnd) 
            if (allocated(G_bnd)) then
                deallocate(G_bnd)
            end if

        end if

        ! Flux function constraints
        !--------------------------
        if (constraints%dofluxfunction) then 
            ! Construct the constraint index
            allocate(conindex(constraints%fluxfunction%ncon))
            conindex = [(k, k = ic+1, ic+constraints%fluxfunction%ncon)]

            ! Allocate & initialize
            allocate(lambda_flux(constraints%fluxfunction%ncon))
            lambda_flux = lambda(conindex)

            ! Call the evaluation routine
            call constraints%fluxfunction%Evaluate(G_flux, &
                gradG_flux, hessG_flux, &
                grid, magneticField, environment, dogradient, &
                dohessian, designvariables, &
                lambda_flux)

            ! Assign
            G(conindex) = G_flux

            ! Update the gradient column indices
            if (dogradient) then 
                ! For easier concatenation later on
                gradG_flux%col = gradG_flux%col + ic

            end if

            ! Update the constraint counter
            ic = ic + constraints%fluxfunction%ncon

            ! Housekeeping
            deallocate(conindex, lambda_flux) 
            if (allocated(G_flux)) then
                deallocate(G_flux)
            end if

        end if

        ! Edge lengths constraints
        !-------------------------
        if (constraints%doedgelengths) then 
            ! Construct the constraint index
            allocate(conindex(constraints%edgelengths%ncon))
            conindex = [(k, k = ic+1, ic+constraints%edgelengths%ncon)]

            ! Allocate & initialize
            allocate(lambda_el(constraints%edgelengths%ncon))
            lambda_el = lambda(conindex)

            ! Call the evaluation routine
            call constraints%edgelengths%Evaluate(G_el, &
                gradG_el, hessG_el, &
                grid, magneticField, environment, dogradient, &
                dohessian, designvariables, &
                lambda_el)

            ! Assign
            G(conindex) = G_el

            ! Update the gradient column indices
            if (dogradient) then 
                ! For easier concatenation later on
                gradG_el%col = gradG_el%col + ic

            end if

            ! Update the constraint counter
            ic = ic + constraints%edgelengths%ncon

            ! Housekeeping
            deallocate(conindex, lambda_el) 
            if (allocated(G_el)) then
                deallocate(G_el)
            end if

        end if

        ! Orthogonality constraints
        !--------------------------
        if (constraints%doorthogonality) then 
            ! Construct the constraint index
            allocate(conindex(constraints%orthogonality%ncon))
            conindex = [(k, k = ic+1, ic+constraints%orthogonality%ncon)]

            ! Allocate & initialize
            allocate(lambda_orth(constraints%orthogonality%ncon))
            lambda_orth = lambda(conindex)

            ! Call the evaluation routine
            call constraints%orthogonality%Evaluate(G_orth, &
                gradG_orth, hessG_orth, &
                grid, magneticField, environment, dogradient, &
                dohessian, designvariables, &
                lambda_orth)

            ! Assign
            G(conindex) = G_orth

            ! Update the gradient column indices
            if (dogradient) then 
                ! For easier concatenation later on
                gradG_orth%col = gradG_orth%col + ic

            end if

            ! Update the constraint counter
            ic = ic + constraints%orthogonality%ncon

            ! Housekeeping
            deallocate(conindex, lambda_orth) 
            if (allocated(G_orth)) then
                deallocate(G_orth)
            end if

        end if


        ! Concatenate gradient
        !=====================
        if (dogradient) then 

            ! Determine sizes
            !----------------
            ! Size of the gradient
            gradG%ncol = constraints%neqcon 
            gradG%nrow = designvariables%nphi

            ! Allocate
            if (.not. allocated(gradG%val)) then 
                ! Number of values (to be determined)
                gradG%nval = 0

                ! Add values of each constraint, if used
                if (constraints%doxpoints) then 
                    gradG%nval = gradG%nval + gradG_xp%nval 
                end if
                if (constraints%doboundaryfunction) then 
                    gradG%nval = gradG%nval + gradG_bnd%nval 
                end if
                if (constraints%dofluxfunction) then 
                    gradG%nval = gradG%nval + gradG_flux%nval  
                end if 
                if (constraints%doedgelengths) then 
                    gradG%nval = gradG%nval + gradG_el%nval  
                end if 
                if (constraints%doorthogonality) then 
                    gradG%nval = gradG%nval + gradG_orth%nval  
                end if 

                ! Allocate
                call gradG%Allocate()
            end if

            ! Add contributions
            !------------------
            ! Initialize counter
            ivg = 0

            ! x-points
            if (constraints%doxpoints) then 
                ! Associate 
                associate(&
                    nc      => constraints%xpoints%ncon, &
                    nval    => gradG_xp%nval)

                ! Add values
                gradG%row(ivg+1:ivg+nval) = gradG_xp%row 
                gradG%col(ivg+1:ivg+nval) = gradG_xp%col
                gradG%val(ivg+1:ivg+nval) = gradG_xp%val

                ! Update counter
                ivg = ivg + nval 

                ! End associate
                end associate

            end if

            ! Boundary function
            if (constraints%doboundaryfunction) then 
                ! Associate 
                associate(&
                    nc      => constraints%boundaryfunction%ncon, &
                    nval    => gradG_bnd%nval)

                ! Add values
                gradG%row(ivg+1:ivg+nval) = gradG_bnd%row 
                gradG%col(ivg+1:ivg+nval) = gradG_bnd%col
                gradG%val(ivg+1:ivg+nval) = gradG_bnd%val

                ! Update counter
                ivg = ivg + nval 

                ! End associate
                end associate

            end if
            
            ! Flux function
            if (constraints%dofluxfunction) then 
                ! Associate 
                associate(&
                    nc      => constraints%fluxfunction%ncon, &
                    nval    => gradG_flux%nval)

                ! Add values
                gradG%row(ivg+1:ivg+nval) = gradG_flux%row 
                gradG%col(ivg+1:ivg+nval) = gradG_flux%col
                gradG%val(ivg+1:ivg+nval) = gradG_flux%val

                ! Update counter
                ivg = ivg + nval 

                ! End associate
                end associate

            end if

            ! Edge lengths
            if (constraints%doedgelengths) then 
                ! Associate 
                associate(&
                    nc      => constraints%edgelengths%ncon, &
                    nval    => gradG_el%nval)

                ! Add values
                gradG%row(ivg+1:ivg+nval) = gradG_el%row 
                gradG%col(ivg+1:ivg+nval) = gradG_el%col
                gradG%val(ivg+1:ivg+nval) = gradG_el%val

                ! Update counter
                ivg = ivg + nval 

                ! End associate
                end associate

            end if

            ! Orthogonality
            if (constraints%doorthogonality) then 
                ! Associate 
                associate(&
                    nc      => constraints%orthogonality%ncon, &
                    nval    => gradG_orth%nval)

                ! Add values
                gradG%row(ivg+1:ivg+nval) = gradG_orth%row 
                gradG%col(ivg+1:ivg+nval) = gradG_orth%col
                gradG%val(ivg+1:ivg+nval) = gradG_orth%val

                ! Update counter
                ivg = ivg + nval 

                ! End associate
                end associate

            end if

        end if

        ! Concatenate the hessian
        !========================
        if (dohessian) then 

            ! Determine sizes
            !----------------
            ! Size of the hessian
            hessG%ncol = designvariables%nphi
            hessG%nrow = designvariables%nphi

            ! Allocate
            if (.not. allocated(hessG%val)) then 
                ! Number of values (to be determined)
                hessG%nval = 0

                ! Add values of each constraint, if used
                if (constraints%doxpoints) then 
                    hessG%nval = hessG%nval + hessG_xp%nval  
                end if 
                if (constraints%doboundaryfunction) then 
                    hessG%nval = hessG%nval + hessG_bnd%nval  
                end if 
                if (constraints%dofluxfunction) then 
                    hessG%nval = hessG%nval + hessG_flux%nval  
                end if 
                if (constraints%doedgelengths) then 
                    hessG%nval = hessG%nval + hessG_el%nval  
                end if 
                if (constraints%doorthogonality) then 
                    hessG%nval = hessG%nval + hessG_orth%nval  
                end if 
                
                ! Allocate
                call hessG%Allocate()

            end if

            ! Add contributions
            !------------------
            ! Initialize counter
            ivh = 0

            ! Boundary function
            if (constraints%doxpoints) then 
                ! Associate 
                associate(&
                    nc      => constraints%xpoints%ncon, &
                    nval    => hessG_xp%nval)

                ! Add values
                hessG%row(ivh+1:ivh+nval) = hessG_xp%row 
                hessG%col(ivh+1:ivh+nval) = hessG_xp%col
                hessG%val(ivh+1:ivh+nval) = hessG_xp%val

                ! Update counter
                ivh = ivh + nval 

                ! End associate
                end associate

            end if

            ! Boundary function
            if (constraints%doboundaryfunction) then 
                ! Associate 
                associate(&
                    nc      => constraints%boundaryfunction%ncon, &
                    nval    => hessG_bnd%nval)

                ! Add values
                hessG%row(ivh+1:ivh+nval) = hessG_bnd%row 
                hessG%col(ivh+1:ivh+nval) = hessG_bnd%col
                hessG%val(ivh+1:ivh+nval) = hessG_bnd%val

                ! Update counter
                ivh = ivh + nval 

                ! End associate
                end associate

            end if
            
            ! Flux function
            if (constraints%dofluxfunction) then 
                ! Associate 
                associate(&
                    nc      => constraints%fluxfunction%ncon, &
                    nval    => hessG_flux%nval)

                ! Add values
                hessG%row(ivh+1:ivh+nval) = hessG_flux%row 
                hessG%col(ivh+1:ivh+nval) = hessG_flux%col
                hessG%val(ivh+1:ivh+nval) = hessG_flux%val

                ! Update counter
                ivh = ivh + nval 

                ! End associate
                end associate

            end if

            ! Edge lengths
            if (constraints%doedgelengths) then 
                ! Associate 
                associate(&
                    nc      => constraints%edgelengths%ncon, &
                    nval    => hessG_el%nval)

                ! Add values
                hessG%row(ivh+1:ivh+nval) = hessG_el%row 
                hessG%col(ivh+1:ivh+nval) = hessG_el%col
                hessG%val(ivh+1:ivh+nval) = hessG_el%val

                ! Update counter
                ivh = ivh + nval 

                ! End associate
                end associate

            end if

            ! Orthogonality
            if (constraints%doorthogonality) then 
                ! Associate 
                associate(&
                    nc      => constraints%orthogonality%ncon, &
                    nval    => hessG_orth%nval)

                ! Add values
                hessG%row(ivh+1:ivh+nval) = hessG_orth%row 
                hessG%col(ivh+1:ivh+nval) = hessG_orth%col
                hessG%val(ivh+1:ivh+nval) = hessG_orth%val

                ! Update counter
                ivh = ivh + nval 

                ! End associate
                end associate

            end if

        end if


    end subroutine


    !------------------------------------------------------------------!
    !                          INEQUALITY CONSTRAINTS                  !
    !------------------------------------------------------------------!
    ! Initialization
    subroutine InitializeIneqCon(constraints, grid, magneticField, &
        environment, constraintoptions, monitor)

        ! Description
        !============
        ! Routine that initializes the desired constraints, based on the
        ! user options d

        ! Declare variables
        !==================
        ! Arguments
        class(IneqConGDUDT)         :: constraints 
        type(GridUDT)               :: grid 
        type(MagneticFieldUDT)      :: magneticField 
        type(EnvironmentUDT)        :: environment 
        type(ConstraintOptionsUDT)  :: constraintoptions
        type(ConstraintsMonitorUDT) :: monitor

        ! Loop variables

        ! Auxiliary variables

        ! Initialize
        !===========

        ! Initialize constraints
        !=======================
        constraints%nineqcon = 0

    end subroutine

    ! Constraint evaluation
    subroutine EvaluateIneqCon(constraints, G, gradG, hessG, &
        grid, magneticField, environment, dogradient, dohessian, & 
        designvariables, lambda)

        ! Description
        !============
        ! This routine evaluates the constraints G and the corresponding
        ! gradient and hessian. To do so, every type of constraint is 
        ! checked whether it is imposed, and the contributions are 
        ! added by calling the evaluation routine of each constraint. 

        ! Declare variables
        !==================
        ! Arguments
        class(IneqConGDUDT)             :: constraints
        real(R8), intent(inout)         :: G(:)
        real(R8), intent(in)            :: lambda(:)
        type(MySparseUDT)               :: gradG, hessG 
        type(GridUDT)                   :: grid
        type(MagneticFieldUDT)          :: magneticField 
        type(EnvironmentUDT)            :: environment
        logical                         :: dogradient, dohessian 
        class(DesignVariablesGDUDT)     :: designvariables 

        ! Loop
        integer(I8)                     :: ic, ivg, ivh, k
        integer(I8), allocatable        :: conindex(:)

        ! Auxiliary

        ! Initialize
        !===========
        ! Set the constraint counter
        ic = 0

        ! Initialize
        G(:) = 0
    
        ! Concatenate gradient
        !=====================
        if (dogradient) then 

            ! Determine sizes
            !----------------
            ! Size of the gradient
            gradG%ncol = constraints%nineqcon 
            gradG%nrow = designvariables%nphi

            ! Allocate
            if (.not. allocated(gradG%val)) then 
                ! Number of values (to be determined)
                gradG%nval = 0

                ! Add values of each constraint, if used

                ! Allocate
                call gradG%Allocate()
            end if

            ! Add contributions
            !------------------
            ! Initialize counter
            ivg = 0
            
        end if

        ! Concatenate the hessian
        !========================
        if (dohessian) then 

            ! Determine sizes
            !----------------
            ! Size of the hessian
            hessG%ncol = designvariables%nphi
            hessG%nrow = designvariables%nphi

            ! Allocate
            if (.not. allocated(hessG%val)) then 
                ! Number of values (to be determined)
                hessG%nval = 0

                ! Add values of each constraint, if used

                ! Allocate
                call hessG%Allocate()
            end if

            ! Add contributions
            !------------------
            ! Initialize counter
            ivh = 0

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

        ! Note 2: the constraint options are passed to this function in 
        ! order to determine whether e.g. boundary nodes should be 
        ! considered for the constraints. 

        ! Note 3: for vessel mode grids (grids reaching up to the 
        ! vessel wall), some vertices may not belong to a flux surface.
        ! To fix them anyway on the vessel wall, the logical 
        ! 'fixfarvesselflux' can be set to 'true'. This is the 
        ! recommended default value. 

        ! Note 4: for target mode grids, the flux value at the 'corners'
        ! is typically also determined by averaging over the field line.
        ! However, this gives issues in some cases where the targets are
        ! nearly flux-aligned. To avoid this, set the logical 
        ! 'fixfluxalignedtargets' to true. This is the recommended 
        ! default value. 

        ! Note 5: it is assumed that the grid contains the flux surface
        ! data and that the flux surfaces are numbered from 1 to nFs

        ! Initialize
        !===========
        ! Modules
        use BicubicSplineInterpolant
        use gdmod_plots
        
        ! Declare variables
        !==================
        ! Arguments 
        class(FluxfunctionConstraintsUDT)       :: constraints 
        type(GridUDT)                           :: grid 
        type(MagneticFieldUDT)                  :: magneticField 
        type(EnvironmentUDT)                    :: environment 
        type(ConstraintsMonitorUDT)             :: monitor

        ! Loop variables
        integer(I8)                 :: i

        ! Auxiliary variables
        integer(I8), allocatable    :: vert_tmp(:), vertID(:) 
        real(R8), allocatable       :: PsiD_tmp(:) 
        logical, allocatable        :: delind(:), mask(:)
        real(R8)                    :: tpsi

        ! Data
        logical                     :: fixfluxalignedtargets = .true. 
        logical                     :: fixfarvesselflux = .true. 

        ! Initialize
        !===========
        ! Number of constraints
        constraints%ncon = 0

        ! Allocate temporary arrays
        allocate(vert_tmp(grid%vert%ntot))
        allocate(PsiD_tmp(grid%vert%ntot))

        ! Allocate auxiliary arrays
        allocate(mask(grid%vert%ntot))
        allocate(delind(grid%vert%ntot))
        allocate(vertID(grid%vert%ntot))

        ! Initialize
        PsiD_tmp(:) = 0
        mask(:) = .false.
        delind(:) = .false.
        vert_tmp = [(i, i = 1, grid%vert%ntot)]

        ! Associate
        associate(&
            vert    => grid%vert,       x       => grid%vert%x,     &
            y       => grid%vert%y,     cc      => monitor%eqvcc,   &
            maxcc   => monitor%maxeqvcc)
        
        ! Determine flux values to impose
        !================================
        ! Evaluate flux at all nodes
        call EvaluateBicubicSplineInterpolant(x, y, PsiD_tmp, &
            magneticField%interp, '0', '0')

        !call Plot2DUnstructuredField(PsiD_tmp, grid, 'v', '-p')
        !call PlotFluxSurfaces(grid, '-p')

        ! Loop over all the flux surfaces to compute desired flux
        do i = 1, grid%data%fluxdata%nFs
            ! Get all vertices with this ID
            mask(:) = (vert%fieldlineID == i)

            ! Get the flux values
            if (any(pack(vert%BV,mask))) then 
                ! Average only over boundary vertices
                tpsi = sum(pack(PsiD_tmp, (mask .and. vert%BV))) & 
                    /count((mask .and. vert%BV))
            else
                ! Average over all vertices
                tpsi = sum(pack(PsiD_tmp, mask))/count(mask)
            end if

            ! Adjust PsiD
            where(mask) PsiD_tmp = tpsi 

        end do

        ! Compensate for flux aligned targets?
        if (fixfluxalignedtargets) then
            ! This is still to do
            print *, 'fix for flux aligned targets is not yet available'
        end if

        ! Fix flux values of vessel boundaries? 
        if (fixfarvesselflux) then 
            ! This is still to do
            print *, 'fix for vessel boundaries is not yet available'
        end if

        ! Set the deletion vector
        where(cc >= maxcc) delind = .true. ! don't constrain

        ! Update monitor
        !===============
        where (.not. delind) cc = cc + 1

        ! Allocate and assign
        !====================
        ! Allocate
        constraints%ncon = count( .not. delind)
        allocate(constraints%vert(constraints%ncon))
        allocate(constraints%PsiD(constraints%ncon))

        ! Assign
        !=======
        constraints%vert = pack(vert_tmp, (.not. delind))
        constraints%PsiD = pack(PsiD_tmp, (.not. delind))

        ! Housekeeping
        !=============
        ! End associate
        end associate

        ! Deallocate temporary arrays
        deallocate(vert_tmp)
        deallocate(PsiD_tmp)

        ! Deallocate auxiliary arrays
        deallocate(mask)
        deallocate(delind)
        deallocate(vertID)
        
    end subroutine

    ! Evaluation
    subroutine EvaluateFluxfunctionConstraints(constraints, G, gradG, & 
        hessG, grid, magneticField, environment, dogradient, &
        dohessian, designvariables, lambda)

        ! Description
        !============
        ! Evaluate the flux function constraints imposed on the 
        ! vertices. For each  vertex considered (see InitDesign), the 
        ! flux function is imposed mathematically as:
        ! 
        !       G_i = Psi(x_i,y_i) - Psi_D,
        !
        ! where x_i and y_i are the i-th vertex's coordinates, Psi is 
        ! the underlying flux function, characterized by a bicubic 
        ! spline interpolant, and Psi_D is a vector containing the
        ! desired flux function values. 

        ! The Hessian of each ith constraint is:
        !
        !       Hjk,i = d^2F/dx^2 if j == k == i
        !       Hjk,i = d^2F/dy^2 if j == k == i+numel(x)
        !       Hjk,i = d^2F/dxdy if j == k+numel(x)
        !       Otherwise zero
        !
        ! Therefore, the multiplication Hjk,i lambda_i is equal to:
        !
        !       Hjk,i lambda_i = d^2F/dx^2 lambda_i (idem for y)

        ! Notes
        !======
        ! Note 1: not the true hessian of the constraint vector is 
        ! returned, but the hessian-vector multiplication with the 
        ! vector lambda, which should be of suitable size. 

        ! Note 2: the row and column indices for the constraints are 
        ! local, meaning that in no way other constraints are accounted
        ! for in positioning the elements in the matrix. This should be
        ! done in an overarching routine. 

        ! Note 3: at first, we compute the linearization of the 
        ! constraints, meaning that we actually compute the Jacobian.
        ! Afterwards, we switch the row and column indices (i.e. 
        ! transpose) to obtain the gradient. 

        ! Initialize
        !===========
        ! Modules
        use BicubicSplineInterpolant
        
        ! Declare variables
        !==================
        ! Arguments 
        class(FluxfunctionConstraintsUDT)   :: constraints 
        real(R8), allocatable               :: G(:) 
        real(R8), allocatable               :: lambda(:)
        type(MySparseUDT)                   :: hessG, gradG, jacG 
        type(GridUDT)                       :: grid 
        type(MagneticFieldUDT)              :: magneticField 
        type(EnvironmentUDT)                :: environment 
        logical                             :: dogradient, dohessian
        class(DesignVariablesGDUDT)         :: designvariables                

        ! Loop variables
        integer(I8)                         :: ic, ivg, ivh, k
        integer(I8), allocatable            :: valindex(:), conindex(:)

        ! Auxiliary variables
        real(R8), allocatable               :: psival(:), dpsidx(:), &
            dpsidy(:), valxx(:), valxy(:), valyy(:)
        integer(I8)                         :: ntv

        ! Data

        ! Initialize
        !===========
        ! Checks
        if ( (.not. allocated(lambda)) .and. dohessian) then
            ! Throw error
            call gdErrorHandler('When evaluating the hessian vector' &
                // ' multiplication, lambda must be given')
        end if

        if (size(lambda) .ne. constraints%ncon) then
            ! Lambda should have the same size as the constraints
            call gdErrorHandler('Lambda should have the same size ' &
                // 'as the constraint vector')
        end if

        ! Counters
        ic = 0 ! constraint counter (local)
        ivg = 0 ! value index for gradient
        ivh = 0 ! value index for hessian

        ! Associate
        associate(&
            nc      => constraints%ncon,        &
            psiD    => constraints%PsiD,        &
            tv      => constraints%vert,        &
            Psifun  => magneticField%interp,    & 
            x       => grid%vert%x,             & 
            y       => grid%vert%y              & 
            )

        ! Constraint value
        !=================
        ! Initialize
        ntv = size(tv)

        ! Allocate
        allocate(psival(nc))
        allocate(G(nc))
        psival(:) = 0

        ! Evaluate
        call EvaluateBicubicSplineInterpolant(x(tv), y(tv), &
            psival, Psifun, '0', '0')
        G(:) = psival - psiD

        ! Constraint gradient
        !====================
        if (dogradient) then 
            ! Initialize
            jacG%nrow = nc 
            jacG%ncol = designvariables%nphi

            ! Check design variables
            select case(designvariables%type)

            case ('coordinates')

                ! Order in jacobian: first x, then y. Has as many 
                ! non-zero elements as there are design variables. 

                ! Allocate
                jacG%nval = 2*ntv
                call jacG%Allocate() 
                allocate(dpsidx(designvariables%nphi))
                allocate(dpsidy(designvariables%nphi))
                allocate(conindex(ntv))
                allocate(valindex(ntv))

                ! Compute the derivative values
                call EvaluateBicubicSplineInterpolant(&
                    x(tv), y(tv), dpsidx, Psifun, '1', '0')
                call EvaluateBicubicSplineInterpolant(&
                    x(tv), y(tv), dpsidy, Psifun, '0', '1')

                ! x-contribution
                !---------------
                ! Build indices
                conindex = [(k, k = ic+1, ic+ntv)]
                valindex = [(k, k = ivg+1, ivg+ntv)]

                ! Add values
                jacG%row(valindex) = conindex  
                jacG%col(valindex) = tv
                jacG%val(valindex) = dpsidx 

                ! y-contribution
                !---------------
                ! Build indices
                ivg = ivg + ntv
                valindex = valindex + ntv

                ! Add values
                jacG%row(valindex) = conindex 
                jacG%col(valindex) = tv + grid%vert%ntot 
                jacG%val(valindex) = dpsidy 

                ! Build gradient
                gradG%nrow = jacG%ncol 
                gradG%ncol = jacG%nrow 
                gradG%nval = jacG%nval 
                
                call gradG%Allocate()
                gradG%row = jacG%col 
                gradG%col = jacG%row
                gradG%val = jacG%val

                ! Housekeeping
                call jacG%Deallocate()

            case default

                ! Unknown, throw error
                call gdErrorHandler('Gradient not implemented for ' &
                    // 'this type of design variable')

            end select

        end if

        ! Constraint hessian
        !===================
        if (dohessian) then 

            ! Initialize
            hessG%nrow = designvariables%nphi 
            hessG%ncol = designvariables%nphi 

            ! Check design variables
            select case(designvariables%type)

            case ('coordinates')
            
                ! Allocate
                hessG%nval = 4*ntv
                if (.not. allocated(valindex)) then
                    allocate(valindex(ntv))
                end if
                if (.not. allocated(conindex)) then 
                    allocate(conindex(ntv))
                end if 
                if (.not. allocated(hessG%val)) then
                    call hessG%Allocate()
                end if
                allocate(valxx(ntv))
                allocate(valxy(ntv))
                allocate(valyy(ntv))

                ! Initialize
                valxx(:) = 0
                valyy(:) = 0
                valxy(:) = 0

                ! Compute contributions
                call EvaluateBicubicSplineInterpolant(&
                    x(tv), y(tv), valxx, Psifun, '2', '0')
                call EvaluateBicubicSplineInterpolant(&
                    x(tv), y(tv), valyy, Psifun, '0', '2')
                call EvaluateBicubicSplineInterpolant(&
                    x(tv), y(tv), valxy, Psifun, '1', '1') ! 

                ! xx-contribution
                !----------------
                k = 1
                ! Build indices
                valindex = [(k, k = ivh+1, ivh+ntv)] 

                ! Add values
                hessG%row(valindex) = tv 
                hessG%col(valindex) = tv 
                hessG%val(valindex) = valxx*lambda ! element-wise mult. 
                
                ! xy-contribution
                !----------------
                ! Build indices
                ivh = ivh + ntv
                valindex = valindex + ntv

                ! Add values
                hessG%row(valindex) = tv
                hessG%col(valindex) = tv + grid%vert%ntot 
                hessG%val(valindex) = valxy*lambda ! element-wise mult. 

                ! yx-contribution
                !----------------
                ! symmetric with xy
                ! Build indices
                ivh = ivh + ntv
                valindex = valindex + ntv 

                ! Add values
                hessG%row(valindex) = tv + grid%vert%ntot 
                hessG%col(valindex) = tv 
                hessG%val(valindex) = valxy*lambda ! element-wise mult. 

                ! yy-contribution
                !----------------
                ! Build indices
                ivh = ivh + ntv
                valindex = valindex + ntv

                ! Add values
                hessG%row(valindex) = tv + grid%vert%ntot 
                hessG%col(valindex) = tv + grid%vert%ntot 
                hessG%val(valindex) = valyy*lambda ! element-wise mult. 

            case default

                ! Unknown, throw error
                call gdErrorHandler('Gradient not implemented for ' &
                    // 'this type of design variable')

            end select

        end if

        ! Housekeeping
        !=============
        ! End associate
        end associate

        ! Deallocate
        deallocate(psival)

        if (dogradient) then 
            deallocate(dpsidx, dpsidy, valindex, conindex)
        end if 

        if (dohessian) then 
            if (allocated(valindex)) then 
                deallocate(valindex, conindex)
            end if 
            deallocate(valxx, valxy, valyy)
        end if

    end subroutine

    ! Destructor
    subroutine DestroyFluxfunctionConstraints(constraints)

        ! Description
        !============
        ! Destructor of the flux function constraints

        ! Declare variables
        !==================
        ! Arguments
        type(FluxfunctionConstraintsUDT)        :: constraints 

        ! Destroy
        !========
        if (allocated(constraints%vert)) then
            deallocate(constraints%vert)
        end if
        if (allocated(constraints%PsiD)) then
            deallocate(constraints%PsiD)
        end if

    end subroutine

    !------------------------------------------------------------------!
    !                         BOUNDARY FUNCTION                        !
    !------------------------------------------------------------------!
    ! Initialize
    subroutine InitializeBoundaryFunctionConstraints(constraints, &
        grid, magneticField, environment, monitor)

        ! Description
        !============
        ! Initialize the boundary function constraints. Here, the 
        ! function equals zero when a point lies precisely on the 
        ! boundary. Though the 'boundary' on which the nodes should lie
        ! can be defined in many different ways, we chose here to base 
        ! it on the vessel bounding polygon, which basically 
        ! encompasses the plasma chamber. 

        ! Only the vertices lying on the target plate are considered
        ! here to be constrained, though one could, if for some reason
        ! desired, also extend this towards other boundaries/boundary
        ! descriptions and other nodes. 

        ! Declare variables
        !==================
        ! Arguments
        class(BoundaryFunctionConstraintsUDT)       :: constraints
        type(GridUDT)                               :: grid 
        type(MagneticFieldUDT)                      :: magneticField 
        type(EnvironmentUDT)                        :: environment 
        type(ConstraintsMonitorUDT)                 :: monitor

        ! Auxiliary
        real(R8), allocatable                       :: xe(:, :), &
            ye(:, :)
        logical                                     :: debugplots
        integer(I8)                                 :: ic, nv 
        integer(I8), allocatable                    :: tv(:)
        logical, allocatable                        :: mask(:)

        ! Loop
        integer(I8)                                 :: i 

        ! Data
        data debugplots /.false./

        ! Initialize
        !===========
        ! Associate
        associate(&
            psf         => constraints%psf,     &
            vessel      => environment%vessel  &
            )

        ! Construct boundary
        !===================
        ! Get vessel edges
        call vessel%polygonset%GetEdges(xe, ye)

        ! Construct polygon description
        call psf%Initialize(xe, ye, size(xe, 1))

        ! Visualize if desired
        if (debugplots) then 
            call psf%Visualize()
        end if 

        ! Set the constraints
        !====================
        ! Compute number of constraints
        constraints%ncon = 0
        constraints%nvert = 0        
        do i = 1, size(grid%bnd)

            ! Check if target plate - hard coded here... 
            if (any(grid%bnd(i)%ID == [1, 2])) then 

                ! Get the current vertices
                allocate(tv(grid%bnd(i)%nvert))
                tv(:) = grid%bnd(i)%vert

                ! Construct the mask
                allocate(mask(grid%bnd(i)%nvert))
                mask(:) = .true.
                
                ! Check the monitor
                where (monitor%eqvcc(tv) .ge. monitor%maxeqvcc) mask = .false.
                nv = count(mask)

                ! Add these nodes
                constraints%ncon = constraints%ncon + nv
                constraints%nvert = constraints%nvert + nv

                ! Deallocate
                deallocate(tv, mask)

            end if

        end do

        ! Allocate
        allocate(constraints%vert(constraints%nvert))

        ! Add vertices
        ic = 0
        do i = 1, size(grid%bnd)
            if (any(grid%bnd(i)%ID == [1, 2])) then 

                ! Get the current vertices
                allocate(tv(grid%bnd(i)%nvert))
                tv(:) = grid%bnd(i)%vert

                ! Construct the mask
                allocate(mask(grid%bnd(i)%nvert))
                mask(:) = .true.
                
                ! Check the monitor
                where (monitor%eqvcc(tv) .ge. monitor%maxeqvcc) mask = .false.
                nv = count(mask)

                ! Add these nodes
                constraints%vert(ic+1:nv) = pack(tv, mask)
                
                ! Update counter
                ic = ic + nv

                ! Update monitor
                monitor%eqvcc(pack(tv, mask)) = &
                    monitor%eqvcc(pack(tv, mask)) + 1

                ! Deallocate
                deallocate(tv, mask)

            end if

        end do

        ! End associate
        end associate




    end subroutine

    ! Evaluation
    subroutine EvaluateBoundaryFunctionConstraints(constraints, G, gradG, & 
        hessG, grid, magneticField, environment, dogradient, &
        dohessian, designvariables, lambda)

        ! Description
        !============
        ! Evaluate the boundary function constraints based on the
        ! boundary function F. For each point, this is given as:
        ! 
        !       G_i = F(x_i,y_i) = 0
        !
        ! where x_i and y_i are the i-th vertex's coordinates, F is 
        ! the underlying boundary shape function which is only zero
        ! on the boundary itself. F is given by a boundary shape 
        ! function type in the constraints (currently it has to be
        ! the same for each node, though this may be extended in 
        ! the future by extending the polygon structure and vertex
        ! IDs)

        ! The Hessian of each ith constraint is:
        !
        !       Hjk,i = d^2F/dx^2 if j == k == i
        !       Hjk,i = d^2F/dy^2 if j == k == i+numel(x)
        !       Hjk,i = d^2F/dxdy if j == k+numel(x)
        !       Otherwise zero
        !
        ! Therefore, the multiplication Hjk,i lambda_i is equal to:
        !
        !       Hjk,i lambda_i = d^2F/dx^2 lambda_i (idem for y)

        ! Notes
        !======
        ! Note 1: not the true hessian of the constraint vector is 
        ! returned, but the hessian-vector multiplication with the 
        ! vector lambda, which should be of suitable size. 

        ! Note 2: the row and column indices for the constraints are 
        ! local, meaning that in no way other constraints are accounted
        ! for in positioning the elements in the matrix. This should be
        ! done in an overarching routine. 

        ! Note 3: at first, we compute the linearization of the 
        ! constraints, meaning that we actually compute the Jacobian.
        ! Afterwards, we switch the row and column indices (i.e. 
        ! transpose) to obtain the gradient. 
        
        ! Declare variables
        !==================
        ! Arguments 
        class(BoundaryFunctionConstraintsUDT)   :: constraints 
        real(R8), allocatable               :: G(:) 
        real(R8), allocatable               :: lambda(:)
        type(MySparseUDT)                   :: hessG, gradG, jacG 
        type(GridUDT)                       :: grid 
        type(MagneticFieldUDT)              :: magneticField 
        type(EnvironmentUDT)                :: environment 
        logical                             :: dogradient, dohessian
        class(DesignVariablesGDUDT)         :: designvariables 

        ! Loop variables
        integer(I8)                         :: ic, ivg, ivh, k
        integer(I8), allocatable            :: valindex(:), conindex(:)

        ! Auxiliary variables
        real(R8), allocatable               :: dpsfdx(:), dpsfdy(:), &
            valxx(:), valxy(:), valyy(:)
        integer(I8)                         :: ntv
        
        ! Initialize
        !===========
        ! Checks
        if ( (.not. allocated(lambda)) .and. dohessian) then
            ! Throw error
            call gdErrorHandler('When evaluating the hessian vector' &
                // ' multiplication, lambda must be given')
        end if

        if (size(lambda) .ne. constraints%ncon) then
            ! Lambda should have the same size as the constraints
            call gdErrorHandler('Lambda should have the same size ' &
                // 'as the constraint vector')
        end if

        ! Counters
        ic = 0 ! constraint counter (local)
        ivg = 0 ! value index for gradient
        ivh = 0 ! value index for hessian

        ! Associate
        associate(&
            nc      => constraints%ncon,        &
            psf     => constraints%psf,         &
            tv      => constraints%vert,        &
            ntv     => constraints%nvert,       &
            x       => grid%vert%x,             & 
            y       => grid%vert%y              & 
            )

        ! Constraint value
        !=================
        ! Allocate
        if (.not. allocated(G)) then 
            allocate(G(nc))
        else
            if (size(G) .ne. nc) then 

                ! Print a warning and reallocate
                print *, 'EvaluateBoundaryFunctionConstraints: ' &
                    // 'Wrong dimension of G, reallocating'
                
                ! Deallocate and reallocate
                deallocate(G)
                allocate(G(nc))

            end if
        end if

        ! Evaluate
        call psf%Evaluate(x(tv), y(tv), G, '0', '0') 

        ! Constraint gradient
        !====================
        if (dogradient) then 
            ! Initialize
            jacG%nrow = nc 
            jacG%ncol = designvariables%nphi

            ! Check design variables
            select case(designvariables%type)

            case ('coordinates')

                ! Order in jacobian: first x, then y. Has as many 
                ! non-zero elements as there are design variables. 

                ! Allocate
                jacG%nval = 2*ntv
                call jacG%Allocate() 
                allocate(dpsfdx(ntv))
                allocate(dpsfdy(ntv))
                allocate(conindex(ntv))
                allocate(valindex(ntv))

                ! Compute the derivative values
                call psf%Evaluate(x(tv), y(tv), dpsfdx, '1', '0')
                call psf%Evaluate(x(tv), y(tv), dpsfdy, '0', '1')

                ! x-contribution
                !---------------
                ! Build indices
                conindex = [(k, k = ic+1, ic+ntv)]
                valindex = [(k, k = ivg+1, ivg+ntv)]

                ! Add values
                jacG%row(valindex) = conindex  
                jacG%col(valindex) = tv
                jacG%val(valindex) = dpsfdx 

                ! y-contribution
                !---------------
                ! Build indices
                ivg = ivg + ntv
                valindex = valindex + ntv

                ! Add values
                jacG%row(valindex) = conindex 
                jacG%col(valindex) = tv + grid%vert%ntot 
                jacG%val(valindex) = dpsfdy 

                ! Build gradient
                gradG%nrow = jacG%ncol 
                gradG%ncol = jacG%nrow 
                gradG%nval = jacG%nval 
                
                call gradG%Allocate()
                gradG%row = jacG%col 
                gradG%col = jacG%row
                gradG%val = jacG%val

                ! Housekeeping
                call jacG%Deallocate()

            case default

                ! Unknown, throw error
                call gdErrorHandler('Gradient not implemented for ' &
                    // 'this type of design variable')

            end select

        end if

        ! Constraint hessian
        !===================
        if (dohessian) then 

            ! Initialize
            hessG%nrow = designvariables%nphi 
            hessG%ncol = designvariables%nphi 

            ! Check design variables
            select case(designvariables%type)

            case ('coordinates')
            
                ! Allocate
                hessG%nval = 4*ntv
                if (.not. allocated(valindex)) then
                    allocate(valindex(ntv))
                end if
                if (.not. allocated(conindex)) then 
                    allocate(conindex(ntv))
                end if 
                if (.not. allocated(hessG%val)) then
                    call hessG%Allocate()
                end if
                allocate(valxx(ntv))
                allocate(valxy(ntv))
                allocate(valyy(ntv))

                ! Initialize
                valxx(:) = 0
                valyy(:) = 0
                valxy(:) = 0

                ! Compute contributions
                call psf%Evaluate(x(tv), y(tv), valxx, '2', '0')
                call psf%Evaluate(x(tv), y(tv), valyy, '0', '2')
                call psf%Evaluate(x(tv), y(tv), valxy, '1', '1')

                ! xx-contribution
                !----------------
                k = 1
                ! Build indices
                valindex = [(k, k = ivh+1, ivh+ntv)] 

                ! Add values
                hessG%row(valindex) = tv 
                hessG%col(valindex) = tv 
                hessG%val(valindex) = valxx*lambda ! element-wise mult. 
                
                ! xy-contribution
                !----------------
                ! Build indices
                ivh = ivh + ntv
                valindex = valindex + ntv

                ! Add values
                hessG%row(valindex) = tv
                hessG%col(valindex) = tv + grid%vert%ntot 
                hessG%val(valindex) = valxy*lambda ! element-wise mult. 

                ! yx-contribution
                !----------------
                ! symmetric with xy
                ! Build indices
                ivh = ivh + ntv
                valindex = valindex + ntv 

                ! Add values
                hessG%row(valindex) = tv + grid%vert%ntot 
                hessG%col(valindex) = tv 
                hessG%val(valindex) = valxy*lambda ! element-wise mult. 

                ! yy-contribution
                !----------------
                ! Build indices
                ivh = ivh + ntv
                valindex = valindex + ntv

                ! Add values
                hessG%row(valindex) = tv + grid%vert%ntot 
                hessG%col(valindex) = tv + grid%vert%ntot 
                hessG%val(valindex) = valyy*lambda ! element-wise mult. 

            case default

                ! Unknown, throw error
                call gdErrorHandler('Gradient not implemented for ' &
                    // 'this type of design variable')

            end select

        end if
        
        ! Housekeeping
        !=============
        ! End associate
        end associate

        ! Deallocate

        if (dogradient) then 
            deallocate(dpsfdx, dpsfdy, valindex, conindex)
        end if 

        if (dohessian) then 
            if (allocated(valindex)) then 
                deallocate(valindex, conindex)
            end if 
            deallocate(valxx, valxy, valyy)
        end if

    end subroutine

    ! Destructor

    !------------------------------------------------------------------!
    !                              X-POINTS                            !
    !------------------------------------------------------------------!

    ! Initialize
    subroutine InitializeXPointConstraints(constraints, &
        grid, magneticField, environment, monitor)

        ! Description
        !============
        ! Initialize the x-point constraints. Here, the x-point is 
        ! constrained to its initial location (so it is based on the
        ! initial grid, not yet on the magnetic field flux gradient).

        ! Declare variables
        !==================
        ! Arguments
        class(XPointConstraintsUDT)                 :: constraints
        type(GridUDT)                               :: grid 
        type(MagneticFieldUDT)                      :: magneticField 
        type(EnvironmentUDT)                        :: environment 
        type(ConstraintsMonitorUDT)                 :: monitor

        ! Auxiliary
        integer(I8)                                 :: nxpind
        integer(I8), allocatable                    :: xpind(:), order(:)

        logical, allocatable                        :: mask(:)

        ! Data

        ! Initialize
        !===========

        ! Get x-points
        !=============
        ! Get
        call DetermineXPoints(xpind, nxpind, order, grid)

        ! Allocate
        allocate(mask(nxpind))

        ! Check which x-points are not yet constrained
        mask(:) = .false. 
        where (monitor%eqvcc(xpind) <= monitor%maxeqvcc-2) mask = .true.

        ! Allocate and assign
        constraints%nxpind = count(mask)
        allocate(constraints%xpind(count(mask)))
        constraints%xpind = pack(xpind, mask)

        ! Get current coordinates, assign
        constraints%locx = grid%vert%x(constraints%xpind)
        constraints%locy = grid%vert%y(constraints%xpind)

        ! Update constraint quantities
        !=============================
        ! Constraints
        constraints%ncon = 2*constraints%nxpind

        ! Monitor
        monitor%eqvcc(constraints%xpind) = &
            monitor%eqvcc(constraints%xpind) + 2

    end subroutine

    ! Evaluation
    subroutine EvaluateXPointConstraints(constraints, G, gradG, & 
        hessG, grid, magneticField, environment, dogradient, &
        dohessian, designvariables, lambda)

        ! Description
        !============
        ! Evaluate the X-point constraints
        ! 
        !       G(2*i-1) = (x_i - x_i0)**2 + (x_i - x_i0) = 0
        !       G(2*i)   = (y_i - y_i0)**2 + (y_i - y_i0) = 0
        !
        ! where x_i and y_i are the i-th x-point coordinates, x_i0, y_i0
        ! are the desired (and constant) x-point coordinates. 
        
        ! The gradient is computed as follows:
        !
        !       J(2*i-1, xind(i))      = 2*(x_i - x_i0) + 1
        !       J(2*i, xind(i) + nv)   = 2*(y_i - y_i0) + 1
        !
        ! where nv is the number of vertices in the grid, and xind the 
        ! vertex ID vector of all x-points considered. 

        ! The Hessian of each ith constraint is:
        !
        !       H(2*i-1, xind(i), xind(i)) = 2
        !       H(2*i, xind(i) + nv, xind(i) + nv) = 2
        !       Otherwise zero
        !
        ! Therefore, the multiplication Hjk,i lambda_i is equal to:
        !
        !       Hjk,i lambda_i = 2 * lambda_i 

        ! Notes
        !======
        ! Note 1: not the true hessian of the constraint vector is 
        ! returned, but the hessian-vector multiplication with the 
        ! vector lambda, which should be of suitable size. 

        ! Note 2: the row and column indices for the constraints are 
        ! local, meaning that in no way other constraints are accounted
        ! for in positioning the elements in the matrix. This should be
        ! done in an overarching routine. 

        ! Note 3: at first, we compute the linearization of the 
        ! constraints, meaning that we actually compute the Jacobian.
        ! Afterwards, we switch the row and column indices (i.e. 
        ! transpose) to obtain the gradient. 
        
        ! Declare variables
        !==================
        ! Arguments 
        class(XPointConstraintsUDT)         :: constraints 
        real(R8), allocatable               :: G(:) 
        real(R8), allocatable               :: lambda(:)
        type(MySparseUDT)                   :: hessG, gradG, jacG 
        type(GridUDT)                       :: grid 
        type(MagneticFieldUDT)              :: magneticField 
        type(EnvironmentUDT)                :: environment 
        logical                             :: dogradient, dohessian
        class(DesignVariablesGDUDT)         :: designvariables 

        ! Loop variables
        integer(I8)                         :: i, ic, ivg, ivh, k
        integer(I8), allocatable            :: valindex(:), conindex(:)

        ! Auxiliary variables
        real(R8), allocatable               :: dpsfdx(:), dpsfdy(:), &
            valxx(:), valxy(:), valyy(:)
        integer(I8)                         :: ntv
        
        ! Initialize
        !===========
        ! Checks
        if ( (.not. allocated(lambda)) .and. dohessian) then
            ! Throw error
            call gdErrorHandler('When evaluating the hessian vector' &
                // ' multiplication, lambda must be given')
        end if

        if (size(lambda) .ne. constraints%ncon) then
            ! Lambda should have the same size as the constraints
            call gdErrorHandler('Lambda should have the same size ' &
                // 'as the constraint vector')
        end if

        ! Counters
        ic = 0 ! constraint counter (local)
        ivg = 0 ! value index for gradient
        ivh = 0 ! value index for hessian

        ! Associate
        associate(&
            nc      => constraints%ncon,        &
            tv      => constraints%xpind,       &
            ntv     => constraints%nxpind,      &
            locx    => constraints%locx,        &
            locy    => constraints%locy,        &
            x       => grid%vert%x,             & 
            y       => grid%vert%y              & 
            )

        ! Constraint value
        !=================
        ! Allocate
        if (.not. allocated(G)) then 
            allocate(G(nc))
        else
            if (size(G) .ne. nc) then 

                ! Print a warning and reallocate
                print *, 'EvaluateBoundaryFunctionConstraints: ' &
                    // 'Wrong dimension of G, reallocating'
                
                ! Deallocate and reallocate
                deallocate(G)
                allocate(G(nc))

            end if
        end if

        ! Evaluate
        do i = 1, ntv
            G(2*i-1) = ( x(tv(i)) - locx(i) )**2 + ( x(tv(i)) - locx(i) )
            G(2*i)   = ( y(tv(i)) - locy(i) )**2 + ( y(tv(i)) - locy(i) )
        end do

        ! Constraint gradient
        !====================
        if (dogradient) then 
            ! Initialize
            jacG%nrow = nc 
            jacG%ncol = designvariables%nphi

            ! Check design variables
            select case(designvariables%type)

            case ('coordinates')

                ! Order in jacobian: first x, then y. 

                ! Allocate
                jacG%nval = 2*ntv
                call jacG%Allocate() 
                allocate(conindex(ntv))
                allocate(valindex(ntv))

                ! x-contribution
                !---------------
                ! Build indices
                conindex = [(k, k = ic+1, ic+2*ntv-1, 2)]
                valindex = [(k, k = ivg+1, ivg+ntv)]

                ! Add values
                jacG%row(valindex) = conindex  
                jacG%col(valindex) = tv
                jacG%val(valindex) = 2*( x(tv) - locx ) + 1

                ! y-contribution
                !---------------
                ! Build indices
                ivg = ivg + ntv
                conindex = [(k, k = ic+2, ic+2*ntv, 2)]
                valindex = valindex + ntv

                ! Add values
                jacG%row(valindex) = conindex 
                jacG%col(valindex) = tv + grid%vert%ntot 
                jacG%val(valindex) = 2*( y(tv) - locy ) + 1 

                ! Build gradient
                gradG%nrow = jacG%ncol 
                gradG%ncol = jacG%nrow 
                gradG%nval = jacG%nval 
                
                call gradG%Allocate()
                gradG%row = jacG%col 
                gradG%col = jacG%row
                gradG%val = jacG%val

                ! Housekeeping
                call jacG%Deallocate()

            case default

                ! Unknown, throw error
                call gdErrorHandler('Gradient not implemented for ' &
                    // 'this type of design variable')

            end select

        end if

        ! Constraint hessian
        !===================
        if (dohessian) then 

            ! Initialize
            hessG%nrow = designvariables%nphi 
            hessG%ncol = designvariables%nphi 

            ! Check design variables
            select case(designvariables%type)

            case ('coordinates')
            
                ! Allocate
                hessG%nval = 4*ntv
                if (.not. allocated(valindex)) then
                    allocate(valindex(ntv))
                end if
                if (.not. allocated(conindex)) then 
                    allocate(conindex(ntv))
                end if 
                if (.not. allocated(hessG%val)) then
                    call hessG%Allocate()
                end if
                allocate(valxx(ntv))
                allocate(valxy(ntv))
                allocate(valyy(ntv))

                ! Compute contributions
                valxx(:) = 2
                valyy(:) = 2
                valxy(:) = 0

                ! xx-contribution
                !----------------
                k = 1
                ! Build indices
                conindex = [(k, k = ic+1, ic+2*ntv-1, 2)]
                valindex = [(k, k = ivh+1, ivh+ntv)] 

                ! Add values
                hessG%row(valindex) = tv 
                hessG%col(valindex) = tv 
                hessG%val(valindex) = valxx*lambda(conindex) ! element-wise mult. 
                
                ! xy-contribution
                !----------------
                ! Build indices
                ivh = ivh + ntv
                valindex = valindex + ntv

                ! Add values
                hessG%row(valindex) = tv
                hessG%col(valindex) = tv + grid%vert%ntot 
                hessG%val(valindex) = 0 ! element-wise mult. 

                ! yx-contribution
                !----------------
                ! symmetric with xy
                ! Build indices
                ivh = ivh + ntv
                valindex = valindex + ntv 

                ! Add values
                hessG%row(valindex) = tv + grid%vert%ntot 
                hessG%col(valindex) = tv 
                hessG%val(valindex) = 0 ! element-wise mult. 

                ! yy-contribution
                !----------------
                ! Build indices
                ivh = ivh + ntv
                conindex = [(k, k = ic+2, ic+2*ntv, 2)]
                valindex = valindex + ntv

                ! Add values
                hessG%row(valindex) = tv + grid%vert%ntot 
                hessG%col(valindex) = tv + grid%vert%ntot 
                hessG%val(valindex) = valyy*lambda(conindex) ! element-wise mult. 

            case default

                ! Unknown, throw error
                call gdErrorHandler('Gradient not implemented for ' &
                    // 'this type of design variable')

            end select

        end if
        
        ! Housekeeping
        !=============
        ! End associate
        end associate

        ! Deallocate

        if (dogradient) then 
            deallocate(valindex, conindex)
        end if 

        if (dohessian) then 
            if (allocated(valindex)) then 
                deallocate(valindex, conindex)
            end if 
            deallocate(valxx, valxy, valyy)
        end if

    end subroutine

    !------------------------------------------------------------------!
    !                           EDGE LENGTHS                           !
    !------------------------------------------------------------------!

    ! Initialize
    subroutine InitializeEdgeLengthsConstraints(constraints, &
        grid, magneticField, environment, monitor)

        ! Description
        !============
        ! Initialize the (desired) edge length constraints. There are 
        ! two main parameters here: 
        !
        ! - the edge (or face in 2D) indices to be constrained
        ! - the desired length (in [m]) of these edges. 
        !
        ! For the former, some preset options are implemented through 
        ! routines in the gdmod_utility_optimization module, which can 
        ! be called here. These currently include:
        !
        ! - edges that have one boundary vertex that lies on the vessel
        !   (these are basically the edges that determine the width of 
        !   the boundary cell near the target plates or other vessel 
        !   segments)
        ! - edges that have an x-point as vertex. 
        !
        ! For the latter (i.e. the length), it is in principle possible 
        ! to apply different lengths for all edges, though the options 
        ! here allow only to specify a uniform edge length automatically
        ! which is sufficient for most purposes. Otherwise, this should
        ! be implemented/given manually. 

        ! Notes
        !======
        ! Note 1: in determining to which vertex to assing this 
        ! constraint (for the monitor), priority is given to the 
        ! boundary vertices, if it is possible to assign it there. 

        ! Declare variables
        !==================
        ! Arguments
        class(EdgeLengthsConstraintsUDT)            :: constraints
        type(GridUDT)                               :: grid 
        type(MagneticFieldUDT)                      :: magneticField 
        type(EnvironmentUDT)                        :: environment 
        type(ConstraintsMonitorUDT)                 :: monitor

        ! Auxiliary
        real(R8)                    :: edgedistvessel, edgedistxpoint

        real(R8), allocatable       :: dvesseledges(:), dxpointedges(:)

        integer(I8)                 :: nvesseledges, nxpointedges, cc, &
            ev(1:2)

        integer(I8), allocatable    :: vesseledges(:, :), &
            xpointedges(:, :), tempvesseledges(:, :), &
            tempxpointedges(:, :)

        logical                     :: dovesseledges, doxpointedges, &
            doTP, doWG
        
        logical, allocatable        :: dovesseledgecon(:), &
            doxpointedgecon(:)
        
        ! Loop
        integer(I8)                 :: i, j

        ! Data

        ! Initialize
        !===========
        ! Associate
        associate(&
            vert        => grid%vert,           &
            vcc         => monitor%eqvcc,       &
            maxvcc      => monitor%maxeqvcc)

        ! Do vessel edge lengths?
        dovesseledges   = .true.
        edgedistvessel  = 5e-3 ! desired edge length in [m]
        doTP            = .true. ! do target plates? 
        doWG            = .true. ! do wide grid boundaries?

        ! Do x-point edge lengths?
        doxpointedges   = .false. 
        edgedistxpoint  = 5e-3 ! desired edge length in [m]

        ! Determine edges
        !================
        ! Edges near vessel 
        !------------------
        nvesseledges = 0
        if (dovesseledges) then 
            call DetermineFluxAlignedVesselEdges(nvesseledges, &
                tempvesseledges, grid, doTP, doWG)

            ! Check for these edges whether they can be constrained
            allocate(dovesseledgecon(nvesseledges))
            dovesseledgecon(:) = .true.
            where ( (vcc(tempvesseledges(:, 1)) >= maxvcc) .and. &
                (vcc(tempvesseledges(:, 2)) >= maxvcc)) &
                dovesseledgecon = .false.

            ! Recompute edges
            nvesseledges = count(dovesseledgecon)
            allocate(vesseledges(nvesseledges, 2))
            vesseledges(:, 1) = pack(tempvesseledges(:, 1), &
                dovesseledgecon)
            vesseledges(:, 2) = pack(tempvesseledges(:, 2), &
                dovesseledgecon)

            ! Set lengths
            allocate(dvesseledges(nvesseledges))
            dvesseledges(:) = edgedistvessel
        end if 

        ! Edges near x-point(s)
        !----------------------
        nxpointedges = 0
        if (doxpointedges) then 
            call DetermineFluxAlignedXPointEdges(nxpointedges, &
                tempxpointedges, grid)

            ! Check for these edges whether they can be constrained
            allocate(doxpointedgecon(nxpointedges))
            doxpointedgecon(:) = .true.
            where ( (vcc(tempxpointedges(:, 1)) >= maxvcc) .and. &
                (vcc(tempxpointedges(:, 2)) >= maxvcc)) &
                doxpointedgecon = .false.

            ! Recompute edges
            nxpointedges = count(doxpointedgecon)
            allocate(xpointedges(nxpointedges, 2))
            xpointedges(:, 1) = pack(tempxpointedges(:, 1), &
                doxpointedgecon)
            xpointedges(:, 2) = pack(tempxpointedges(:, 2), &
                doxpointedgecon)

            ! Set lengths
            allocate(dxpointedges(nxpointedges))
            dxpointedges(:) = edgedistxpoint
        end if

        ! Update constraint quantities
        !=============================
        ! Constraints
        constraints%ncon = nvesseledges + nxpointedges
        constraints%nedges = nvesseledges + nxpointedges 
        allocate(constraints%edgevert(constraints%nedges, 2))
        allocate(constraints%d(constraints%nedges))

        cc = 0 ! constraint counter
        if (dovesseledges) then 
            do j = 1, 2
                constraints%edgevert(cc+1:cc+nvesseledges, j) = &
                    vesseledges(:, j) 
            end do
            constraints%d(cc+1:cc+nvesseledges) = &
                dvesseledges
            cc = cc + nvesseledges 
        end if
        if (doxpointedges) then 
            do j = 1, 2
                constraints%edgevert(cc+1:cc+nxpointedges, j) = &
                    xpointedges(:, j)
            end do
            constraints%d(cc+1:cc+nxpointedges) = &
                dxpointedges
            cc = cc + nxpointedges
        end if

        ! Monitor
        do i = 1, constraints%nedges 
            ! Unpack for ease
            ev = constraints%edgevert(i, 1:2)

            ! Update counter
            if (vert%BV(ev(1)) .and. (vcc(ev(1)) < maxvcc)) then 
                ! Assign to boundary vertex
                vcc(ev(1)) = vcc(ev(1)) + 1
            elseif (vert%BV(ev(2)) .and. (vcc(ev(2)) < maxvcc)) then 
                ! Assign to boundary vertex
                vcc(ev(2)) = vcc(ev(2)) + 1
            elseif (vcc(ev(1)) <= maxvcc) then 
                ! Assign to first vertex, no boundary vertex
                vcc(ev(1)) = vcc(ev(1)) + 1
            elseif (vcc(ev(2)) <= maxvcc) then 
                ! Assign to first vertex, no boundary vertex
                vcc(ev(2)) = vcc(ev(2)) + 1
            else
                ! Something wrong - indicates that this edge shouldn't 
                ! have been added, though this should've been catched 
                ! before
                call gdErrorHandler('InitializeEdgelengthsConstraints:' &
                    // 'edge with vertex IDs as shown above should not' &
                    // ' have been included in the constraints')
            end if
        end do

        ! Housekeeping
        !=============
        ! Deallocate
        if (dovesseledges) then 
            deallocate(vesseledges, tempvesseledges, dovesseledgecon, &
            dvesseledges)
        end if 
        if (doxpointedges) then 
            deallocate(xpointedges, tempxpointedges, doxpointedgecon, &
            dxpointedges)
        end if

        ! Deassociate
        end associate

    end subroutine

    ! Evaluation
    subroutine EvaluateEdgeLengthsConstraints(constraints, G, gradG, & 
        hessG, grid, magneticField, environment, dogradient, &
        dohessian, designvariables, lambda)

        ! Description
        !============
        ! The edge length constraint for the i-th edge (and therefore 
        ! in this case the i-th constraint) can be mathematically 
        ! formulated as: 
        !
        !       G(i) = 0.5 * (l_i ** 2  - d0 ** 2)
        !
        ! This form ensures that at least one part of the gradient 
        ! w.r.t. the coordinate values is non-zero (for non-zero 
        ! desired distance, which should always be the case!), while 
        ! keeping the constraint formulation as a simple quadratic 
        ! function in terms of the vertex coordinates. Here, l_i is 
        ! defined as the L2 norm of the edge vector, where the latter is
        ! formed by the vector between its vertices.
        !
        ! The Jacobian w.r.t. the vertex coordinates is then found as
        !
        !       J(i, j) = xv(:, 1) - xv(:, 2), if j = edgevert(i, 1)
        !       J(i, j) = -(xv(:, 1) - xv(:, 2)), if j = edgevert(i, 2)
        !       J(i, j+nv) = yv(:, 1) - yv(:, 2), if j = edgevert(i, 1)
        !       J(i, j+nv) = -(yv(:, 1) - yv(:, 2)), if j = edgevert(i, 2)
        !
        ! where i is the constraint index, j is the vertex index of the 
        ! i-th edge (first or second, as specified above), and nv is the
        ! number of grid vertices.
        !
        ! The Hessian w.r.t. the grid vertex coordinates is found as
        !
        !       H(i, j, k) = 1 (j = k)
        !       H(i, j, k) = -1 (j != k)
        !       H(i, j+nv, k+nv) = 1 (j = k)
        !       H(i, j+nv, k+nv) = -1 (j != k)
        !       otherwise zero

        ! Notes
        !======
        ! Note 1: not the true hessian of the constraint vector is 
        ! returned, but the hessian-vector multiplication with the 
        ! vector lambda, which should be of suitable size. 

        ! Note 2: the row and column indices for the constraints are 
        ! local, meaning that in no way other constraints are accounted
        ! for in positioning the elements in the matrix. This should be
        ! done in an overarching routine. 

        ! Note 3: at first, we compute the linearization of the 
        ! constraints, meaning that we actually compute the Jacobian.
        ! Afterwards, we switch the row and column indices (i.e. 
        ! transpose) to obtain the gradient. 
        
        ! Declare variables
        !==================
        ! Arguments 
        class(EdgeLengthsConstraintsUDT)    :: constraints 
        real(R8), allocatable               :: G(:) 
        real(R8), allocatable               :: lambda(:)
        type(MySparseUDT)                   :: hessG, gradG, jacG 
        type(GridUDT)                       :: grid 
        type(MagneticFieldUDT)              :: magneticField 
        type(EnvironmentUDT)                :: environment 
        logical                             :: dogradient, dohessian
        class(DesignVariablesGDUDT)         :: designvariables 

        ! Loop variables
        integer(I8)                         :: i, ic, ivg, ivh, k
        integer(I8), allocatable            :: valindex(:), conindex(:)

        ! Auxiliary variables
        real(R8), allocatable               :: valxx(:), valxy(:), &
            valyy(:), xv1(:), xv2(:), yv1(:), yv2(:), dist(:)
        integer(I8)                         :: ntv
        
        ! Initialize
        !===========
        ! Checks
        if ( (.not. allocated(lambda)) .and. dohessian) then
            ! Throw error
            call gdErrorHandler('When evaluating the hessian vector' &
                // ' multiplication, lambda must be given')
        end if

        if (size(lambda) .ne. constraints%ncon) then
            ! Lambda should have the same size as the constraints
            call gdErrorHandler('Lambda should have the same size ' &
                // 'as the constraint vector')
        end if

        ! Counters
        ic = 0 ! constraint counter (local)
        ivg = 0 ! value index for gradient
        ivh = 0 ! value index for hessian

        ! Associate
        associate(&
            nc      => constraints%ncon,        &
            ev      => constraints%edgevert,    &
            d0      => constraints%d,           &
            nv      => grid%vert%ntot,          &
            x       => grid%vert%x,             & 
            y       => grid%vert%y              & 
            )

        ! Constraint value
        !=================
        ! Allocate
        if (.not. allocated(G)) then 
            allocate(G(nc))
        else
            if (size(G) .ne. nc) then 

                ! Print a warning and reallocate
                print *, 'EvaluateEdgelengthsConstraints: ' &
                    // 'Wrong dimension of G, reallocating'
                
                ! Deallocate and reallocate
                deallocate(G)
                allocate(G(nc))

            end if
        end if

        ! Compute lengths
        allocate(xv1(nc), xv2(nc), yv1(nc), yv2(nc), dist(nc))
        xv1 = x(ev(:, 1))
        yv1 = y(ev(:, 1))
        xv2 = x(ev(:, 2))
        yv2 = y(ev(:, 2))

        ! Evaluate
        dist = (xv1 - xv2)**2 + (yv1 - yv2)**2
        G = 0.5*(dist - d0**2)

        ! Constraint gradient
        !====================
        if (dogradient) then 
            ! Initialize
            jacG%nrow = nc 
            jacG%ncol = designvariables%nphi

            ! Check design variables
            select case(designvariables%type)

            case ('coordinates')

                ! Order in jacobian: first x, then y. 

                ! Allocate
                jacG%nval = 4*nc
                call jacG%Allocate() 
                allocate(conindex(nc))
                allocate(valindex(nc))

                ! Build constraint indices
                conindex = [(k, k = ic+1, ic+nc)]

                ! x-contribution
                !---------------
                ! Build indices for xv1
                valindex = [(k, k = ivg+1, ivg+nc)]

                ! Add values
                jacG%row(valindex) = conindex  
                jacG%col(valindex) = ev(:, 1)
                jacG%val(valindex) = xv1 - xv2

                ! Update counters
                ivg = ivg + nc

                ! Build indices for xv2
                valindex = [(k, k = ivg+1, ivg+nc)]

                ! Add values
                jacG%row(valindex) = conindex  
                jacG%col(valindex) = ev(:, 2)
                jacG%val(valindex) = -(xv1 - xv2)

                ! Update counters
                ivg = ivg + nc

                ! y-contribution
                !---------------
                ! Build indices for yv1
                valindex = [(k, k = ivg+1, ivg+nc)]

                ! Add values
                jacG%row(valindex) = conindex 
                jacG%col(valindex) = ev(:, 1) + grid%vert%ntot 
                jacG%val(valindex) = yv1 - yv2

                ! Update counters
                ivg = ivg + nc

                ! Build indices for yv2
                valindex = [(k, k = ivg+1, ivg+nc)]

                ! Add values
                jacG%row(valindex) = conindex 
                jacG%col(valindex) = ev(:, 2) + grid%vert%ntot 
                jacG%val(valindex) = -(yv1 - yv2)

                ! Update counters
                ivg = ivg + nc

                ! Build gradient
                !---------------
                gradG%nrow = jacG%ncol 
                gradG%ncol = jacG%nrow 
                gradG%nval = jacG%nval 
                
                call gradG%Allocate()
                gradG%row = jacG%col 
                gradG%col = jacG%row
                gradG%val = jacG%val

                ! Housekeeping
                call jacG%Deallocate()

            case default

                ! Unknown, throw error
                call gdErrorHandler('Gradient not implemented for ' &
                    // 'this type of design variable')

            end select

        end if

        ! Constraint hessian
        !===================
        if (dohessian) then 

            ! Initialize
            ic = 0
            hessG%nrow = designvariables%nphi 
            hessG%ncol = designvariables%nphi 

            ! Check design variables
            select case(designvariables%type)

            case ('coordinates')
            
                ! Allocate
                hessG%nval = 8*nc 
                if (.not. allocated(valindex)) then
                    allocate(valindex(nc))
                end if
                if (.not. allocated(conindex)) then 
                    allocate(conindex(nc))
                end if 
                if (.not. allocated(hessG%val)) then
                    call hessG%Allocate()
                end if
                allocate(valxx(nc))
                allocate(valxy(nc))
                allocate(valyy(nc))

                ! Compute contributions
                valxx(:) = 1
                valyy(:) = 1
                valxy(:) = 0

                ! Build constraint indices
                conindex = [(k, k = ic+1, ic+nc)]

                ! xx-contribution
                !----------------
                ! Build indices
                valindex = [(k, k = ivh+1, ivh+nc)] 

                ! Add values
                hessG%row(valindex) = ev(:, 1) 
                hessG%col(valindex) = ev(:, 1) 
                hessG%val(valindex) = valxx*lambda(conindex) ! x1x1

                ! Update counters
                ivh = ivh + nc

                ! Build indices
                valindex = [(k, k = ivh+1, ivh+nc)] 

                ! Add values
                hessG%row(valindex) = ev(:, 1) 
                hessG%col(valindex) = ev(:, 2) 
                hessG%val(valindex) = -valxx*lambda(conindex) ! x1x2

                ! Update counters
                ivh = ivh + nc

                ! Build indices
                valindex = [(k, k = ivh+1, ivh+nc)] 

                ! Add values
                hessG%row(valindex) = ev(:, 2) 
                hessG%col(valindex) = ev(:, 1) 
                hessG%val(valindex) = -valxx*lambda(conindex) ! x2x1

                ! Update counters
                ivh = ivh + nc

                ! Build indices
                valindex = [(k, k = ivh+1, ivh+nc)] 

                ! Add values
                hessG%row(valindex) = ev(:, 2) 
                hessG%col(valindex) = ev(:, 2) 
                hessG%val(valindex) = valxx*lambda(conindex) ! x2x2

                ! Update counters
                ivh = ivh + nc
                
                ! xy-contribution
                !----------------
                ! no contributions

                ! yx-contribution
                !----------------
                ! no contributions

                ! yy-contribution
                !----------------
                ! Build indices
                valindex = [(k, k = ivh+1, ivh+nc)] 

                ! Add values
                hessG%row(valindex) = ev(:, 1) + nv
                hessG%col(valindex) = ev(:, 1) + nv 
                hessG%val(valindex) = valxx*lambda(conindex) ! y1y1

                ! Update counters
                ivh = ivh + nc

                ! Build indices
                valindex = [(k, k = ivh+1, ivh+nc)] 

                ! Add values
                hessG%row(valindex) = ev(:, 1) + nv 
                hessG%col(valindex) = ev(:, 2) + nv 
                hessG%val(valindex) = -valxx*lambda(conindex) ! y1y2

                ! Update counters
                ivh = ivh + nc

                ! Build indices
                valindex = [(k, k = ivh+1, ivh+nc)] 

                ! Add values
                hessG%row(valindex) = ev(:, 2) + nv 
                hessG%col(valindex) = ev(:, 1) + nv 
                hessG%val(valindex) = -valxx*lambda(conindex) ! y2y1

                ! Update counters
                ivh = ivh + nc

                ! Build indices
                valindex = [(k, k = ivh+1, ivh+nc)] 

                ! Add values
                hessG%row(valindex) = ev(:, 2) + nv 
                hessG%col(valindex) = ev(:, 2) + nv 
                hessG%val(valindex) = valxx*lambda(conindex) ! y2y2

                ! Update counters
                ivh = ivh + nc

            case default

                ! Unknown, throw error
                call gdErrorHandler('Gradient not implemented for ' &
                    // 'this type of design variable')

            end select

        end if
        
        ! Housekeeping
        !=============
        ! End associate
        end associate

        ! Deallocate

        if (dogradient) then 
            deallocate(valindex, conindex)
        end if 

        if (dohessian) then 
            if (allocated(valindex)) then 
                deallocate(valindex, conindex)
            end if 
            deallocate(valxx, valxy, valyy)
        end if

    end subroutine

    !------------------------------------------------------------------!
    !                           ORTHOGONALITY                          !
    !------------------------------------------------------------------!

    ! Initialize
    subroutine InitializeOrthogonalityConstraints(constraints, &
        grid, magneticField, environment, monitor)

        ! Description
        !============
        ! Initialize the edge length constraints on the desired edges.
        ! As there are many ways to determine the edges, this is 
        ! cast into a function located in the gdmod_utility_optimization
        ! module, i.e. DetermineEdgesOrthogonalityConstraints. Different
        ! preset options are present:
        ! 
        ! - simple boxes that determine, starting from the original
        !   grid, which edges to be constrained (or not) by checking if
        !   edges are within the box. 
        ! - flux value based: edges with a flux value between certain 
        !   limits are constrained (useful for e.g. the core)
        ! - initial orthogonality based: numerically compute the 
        !   deviation from orthogonality of the original edges and 
        !   decide based on that whether to include the edges. 
        !
        ! See the routine for more options and details. 

        ! Notes
        !======

        ! Declare variables
        !==================
        ! Arguments
        class(OrthogonalityConstraintsUDT)          :: constraints
        type(GridUDT)                               :: grid 
        type(MagneticFieldUDT)                      :: magneticField 
        type(EnvironmentUDT)                        :: environment 
        type(ConstraintsMonitorUDT)                 :: monitor

        ! Auxiliary
        integer(I8)                 :: nincludebox, nexcludebox, tv, &
            startindex, endindex, tID, nbID(1:2), vpc, tbv, tnbv, tf
        real(R8)                    :: tx, ty, tn, bx, by, bn, nb, dotprod, &
            epsortho, epsperp
        logical                     :: checkperp, isfaceperp, &
            doincludebox, doexcludebox, debugplots 
        
        integer(I8), allocatable    :: cvertlist(:), temp(:), &
            northcon(:), maxnorthcon(:), tvn(:), vpairs(:, :)

        real(R8), allocatable       :: includeboxlbx(:), &
            includeboxubx(:), includeboxlby(:), includeboxuby(:), &
            excludeboxubx(:), excludeboxlby(:), excludeboxuby(:), &
            excludeboxlbx(:), Btx(:), Bty(:), xf(:), yf(:)

        logical, allocatable        :: cvert(:), boxcheck(:), &
            movetoback(:), movetofront(:), ismarked(:), isconstrained(:), &
            isperp(:), cID(:)
        
        ! Loop
        integer(I8)                 :: i, j, k

        ! Data

        ! Initialize
        !===========
        ! Associate
        associate(&
            interp      => magneticField%interp,    &
            vert        => grid%vert,               &
            vcc         => monitor%eqvcc,           &
            maxvcc      => monitor%maxeqvcc,        &
            x           => grid%vert%x,             &
            y           => grid%vert%y)

        ! Debug plots?
        debugplots = .false.

        ! Boxes for edges to be included?
        doincludebox    = .true.
        nincludebox     = 1 ! number of boxes

        ! Boxes for edges to be excluded
        doexcludebox    = .false. 
        nexcludebox     = 1 ! number of boxes

        ! Tolerances (to be moved to input)
        epsortho    = 0.1
        epsperp     = 0.1 

        ! Check perpendicularity? (to be moved to input)
        checkperp = .true.

        ! Allocate
        allocate(&
            includeboxlbx(nincludebox), includeboxlby(nincludebox), &
            includeboxubx(nincludebox), includeboxuby(nincludebox), &
            excludeboxlbx(nexcludebox), excludeboxlby(nexcludebox), &
            excludeboxubx(nexcludebox), excludeboxuby(nexcludebox))

        ! Set the values of the boxes (to be read in in the future)
        includeboxlbx(1:nincludebox) = [-1]
        includeboxubx(1:nincludebox) = [1]
        includeboxlby(1:nincludebox) = [-0.2]
        includeboxuby(1:nincludebox) = [1]

        excludeboxlbx(1:nexcludebox) = [-1]
        excludeboxubx(1:nexcludebox) = [1]
        excludeboxlby(1:nexcludebox) = [-1]
        excludeboxuby(1:nexcludebox) = [1]

        ! Determine edge vertices
        !========================
        ! Compute magnetic field vector components Btx, Bty
        allocate(Btx(vert%ntot), Bty(vert%ntot))
        call EvaluateBicubicSplineInterpolant(vert%x, vert%y, Btx, interp, &
            '0', '1')
        call EvaluateBicubicSplineInterpolant(vert%x, vert%y, Bty, interp, &
            '1', '0')
        Btx = -Btx ! adjust sign 

        ! Determine which nodes to consider
        allocate(cvert(vert%ntot))
        allocate(boxcheck(vert%ntot))
        allocate(northcon(vert%ntot), maxnorthcon(vert%ntot))

        cvert(:) = .false.
        northcon(:) = 0
        maxnorthcon(:) = 2
        where (vert%BV) maxnorthcon = 1

        if (doincludebox) then 
            do i = 1, nincludebox ! include points in the box
                boxcheck(:) = .false.
                boxcheck = (vert%x >= includeboxlbx(i)) .and. &
                    (vert%x <= includeboxubx(i)) .and. &
                    (vert%y >= includeboxlby(i)) .and. &
                    (vert%y <= includeboxuby(i))
                cvert = cvert .or. boxcheck 
            end do
        end if 

        if (doexcludebox) then 
            do i = 1, nincludebox ! exclude points outside the box
                boxcheck(:) = .false.
                boxcheck = (vert%x >= includeboxlbx(i)) .and. &
                    (vert%x <= includeboxubx(i)) .and. &
                    (vert%y >= includeboxlby(i)) .and. &
                    (vert%y <= includeboxuby(i))
                where (boxcheck) cvert = .false.
            end do
        end if

        ! Prioritize inner vertices (typically yields better results, 
        ! but may not be a general approach)
        allocate(cvertlist(count(cvert)))
        allocate(ismarked(vert%ntot))
        allocate(movetoback(count(cvert)))
        allocate(movetofront(count(cvert)))

        ismarked(:) = .false. 
        where (cvert) ismarked = .true. 
        movetoback(:) = .false. 
        movetofront(:) = .false. 
        cvertlist = pack( [(k, k = 1, vert%ntot)], cvert) ! node indices

        do i = 1, count(cvert)
            ! Get the current vertex 
            tv = cvertlist(i)

            ! Get the neighbours
            allocate(tvn(vert%neigP(tv, 2)))
            tvn = vert%neiglist(vert%neigP(tv, 1):vert%neigP(tv, 1)+vert%neigP(tv, 2)-1)

            ! Check if any are not marked
            if (any(.not. ismarked(tvn))) then
                movetoback(i) = .true.
            end if

            ! Check if it is an x-point
            if (.not. movetoback(i)) then
                if (vert%neigP(tv, 2) > 4) then ! crude check for x-point 
                    ! Prioritize, move to front 
                    movetofront(i) = .true. 
                end if
            end if

            ! Deallocate 
            deallocate(tvn)
        end do

        ! Rebuild cvertlist
        allocate(temp(size(cvertlist, 1)))
        temp = cvertlist
        startindex = 1
        endindex = count(movetofront)
        cvertlist(startindex:endindex) = pack(temp, movetofront)
        startindex = startindex + endindex 
        endindex = endindex + count((.not. movetoback) .and. (.not. movetofront))
        cvertlist(startindex:endindex) = pack(temp, (.not. movetoback) .and. (.not. movetofront))
        startindex = endindex + 1
        endindex = endindex + count(movetoback)
        cvertlist(startindex:endindex) = pack(temp, movetoback)

        deallocate(temp)
        deallocate(movetoback, movetofront)
        
        ! Initialize   
        allocate(isconstrained(grid%faces%ntot))
        isconstrained(:) = .false.

        ! Loop 
        vpc = 0 ! face counter
        allocate(vpairs(grid%faces%ntot, 2))! allocate too big, trim later

        do i = 1, size(cvertlist, 1)
            ! Get the current vertex
            tv = cvertlist(i)

            ! Get the neighbours of this vertex
            allocate(tvn(vert%neigP(tv, 2)))
            tvn = vert%neiglist(vert%neigP(tv, 1):vert%neigP(tv, 1)+vert%neigP(tv, 2)-1)

            ! Get the fieldline ID 
            tID = vert%fieldlineID(tv)

            ! Check if zero 
            if (tID .eq. 0) then ! this is a vertex without fieldline 
                ! Check which faces are aligned 
                allocate(isperp(size(tvn, 1)))

                ! Loop over all vertex neighbours 
                do j = 1, size(tvn, 1)
                    ! Get normalized face vector 
                    tx = x(tvn(j)) - x(tv)
                    ty = y(tvn(j)) - y(tv)
                    tn = (tx**2 + tx**2)**0.5
                    tx = tx/tn 
                    ty = ty/tn
                    
                    ! Get normalized magnetic field vector 
                    bx = (Btx(tv) + Btx(tvn(j)))*0.5
                    by = (Bty(tv) + Bty(tvn(j)))*0.5
                    bn = (bx**2 + by**2)**0.5
                    bx = bx/bn ! normalize 
                    by = by/bn

                    ! Compute dot product 
                    dotprod = tx*bx + ty*by 

                    ! Check 
                    if (abs(dotprod) < epsortho) then 
                        isperp(j) = .true.
                    endif 

                end do

                ! Retain perpendicular faces only if 2 found
                if (count(isperp) == 2) then
                    allocate(temp(size(tvn, 1)))
                    temp = tvn 
                    deallocate(tvn)
                    allocate(tvn(count(isperp)))
                    tvn = pack(temp, isperp)

                    deallocate(temp)
                else 
                    deallocate(tvn)
                    allocate(tvn(0))
                end if

                ! Deallocate 
                deallocate(isperp)
            else
                ! Check which vertices have the same ID 
                allocate(cID(size(tvn, 1)))
                cID = tID .eq. vert%fieldlineID(tvn)

                ! Extract vertices that have NOT the same ID
                allocate(temp(size(tvn, 1)))
                temp = tvn
                deallocate(tvn)
                allocate(tvn(count(.not. cID)))
                tvn = pack(temp, .not. cID)
                deallocate(temp)

                ! If two nodes remain, and they have different field
                ! line IDs, then constrain. Otherwise don't (probably 
                ! stacked triangles)
                if (size(tvn, 1) == 2) then 
                    nbID = vert%fieldlineID(tvn)
                    if (nbID(1) == nbID(2)) then 
                        deallocate(tvn)
                        allocate(tvn(0))
                    end if 
                else
                    deallocate(tvn)
                    allocate(tvn(0))
                end if

                ! Deallocate 
                deallocate(cID)
            end if

            ! Constrain each pair 
            do j = 1, size(tvn, 1)
                ! Get the face 
                call MapVertexPairToFace(tv, tvn(j), grid%faces%vert, &
                    grid%faces%ntot, tf)

                ! Check initial perpendicularity
                isfaceperp = .true.
                if (checkperp) then 
                    ! Get normalized face vector 
                    tx = x(tvn(j)) - x(tv)
                    ty = y(tvn(j)) - y(tv)
                    tn = (tx**2 + tx**2)**0.5
                    tx = tx/tn 
                    ty = ty/tn
                    
                    ! Get normalized magnetic field vector 
                    bx = (Btx(tv) + Btx(tvn(j)))*0.5
                    by = (Bty(tv) + Bty(tvn(j)))*0.5
                    bn = (bx**2 + by**2)**0.5
                    bx = bx/bn ! normalize 
                    by = by/bn

                    ! Compute dot product 
                    dotprod = tx*bx + ty*by 

                    ! Check 
                    if (abs(dotprod) > epsperp) then
                        ! Not perpendicular, don't consider this face 
                        isfaceperp = .false.
                    end if 
                end if

                ! Check if there is a boundary vertex in this edge which 
                ! already has been constrained. 
                if ( (.not. isconstrained(tf)) & ! should not be constrained
                    .and. ( (vcc(tvn(j)) < maxvcc) .or. (vcc(tv) < maxvcc) ) & ! a vertex has less than 2 constraints already imposed
                    .and. ( (.not. vert%BV(tvn(j))) .or. (.not. vert%BV(tv))) & ! at least one is an internal vertex
                    .and. ( isfaceperp ) & ! the face is initially almost orthogonal 
                    .and. ( ( northcon(tvn(j)) < maxnorthcon(tvn(j)) ) .and. ( northcon(tv) < maxnorthcon(tv) ) ) &
                    ) then

                    ! Add
                    vpc = vpc + 1
                    vpairs(vpc, :) = [tv, tvn(j)]

                    ! Update counter - first attribute to BV if possible
                    if (vert%BV(tv) .or. vert%BV(tvn(j)) ) then 
                        if (vert%BV(tv)) then 
                            tbv = tv 
                            tnbv = tvn(j)
                        else
                            tbv = tvn(j) 
                            tnbv = tv 
                        end if

                        ! Check if we can add the constraint there
                        if ( ( northcon(tbv) < maxnorthcon(tbv) ) &
                            .and. vcc(tbv) < maxvcc ) then 
                            vcc(tbv) = vcc(tbv) + 1 
                        else
                            vcc(tnbv) = vcc(tnbv) + 1 
                        end if
                    else
                        if ( vcc(tv) < maxvcc ) then 
                            vcc(tv) = vcc(tv) + 1
                        else
                            vcc(tvn(j)) = vcc(tvn(j)) +  1
                        end if 
                    end if

                    ! Update counters 
                    northcon(tvn(j)) = northcon(tvn(j)) + 1
                    northcon(tv) = northcon(tv) + 1
                    isconstrained(tf) = .true. 

                end if 
                    
            end do

            ! Deallocate 
            deallocate(tvn)
        end do

        ! Update constraint quantities
        !=============================
        ! Constraints
        constraints%ncon = vpc 
        constraints%nedges = vpc
        allocate(constraints%edgevert(constraints%nedges, 2))
        constraints%edgevert = vpairs(1:vpc, :)

        ! Visualize
        if (debugplots) then 

            ! Plot the faces that are constrained

            ! Compute face centers
            allocate(xf(vpc), yf(vpc))
            xf = 0.5*(x(constraints%edgevert(:, 1)) + x(constraints%edgevert(:, 2)))
            yf = 0.5*(y(constraints%edgevert(:, 1)) + y(constraints%edgevert(:, 2)))
            call PlotGridWithPoints(grid, xf, yf, '-p')
            deallocate(xf, yf)

        end if
        

        ! Housekeeping
        !=============
        ! Deallocate
        deallocate(includeboxlbx, includeboxlby, includeboxubx, &
            includeboxuby, excludeboxlbx, excludeboxlby, excludeboxubx, &
            excludeboxuby, Btx, Bty, cvert, boxcheck, vpairs)

        ! Deassociate
        end associate

    end subroutine

    ! Evaluation
    subroutine EvaluateOrthogonalityConstraints(constraints, G, gradG, & 
        hessG, grid, magneticField, environment, dogradient, &
        dohessian, designvariables, lambda)

        ! Description
        !============
        ! The orthogonality constraints are evaluated per vertex pair, 
        ! where one tries to satisfy the following condition:
        !
        !       gx*dx + gy*dy = 0.
        !
        ! This ensures orthogonality of the face w.r.t. the magnetic 
        ! field, given by its components gx and gy. dx and dy are the 
        ! face tangents (they can be normalized, but this is not 
        ! strictly necessary, same for the magnetic field. This does, 
        ! however, influence the absolute value of the constraint). 


        ! Notes
        !======
        
        ! Declare variables
        !==================
        ! Arguments 
        class(OrthogonalityConstraintsUDT)  :: constraints 
        real(R8), allocatable               :: G(:) 
        real(R8), allocatable               :: lambda(:)
        type(MySparseUDT)                   :: hessG, gradG, jacG 
        type(GridUDT)                       :: grid 
        type(MagneticFieldUDT)              :: magneticField 
        type(EnvironmentUDT)                :: environment 
        logical                             :: dogradient, dohessian
        class(DesignVariablesGDUDT)         :: designvariables 

        ! Loop variables
        integer(I8)                         :: i, ic, ivg, ivh, k
        integer(I8), allocatable            :: valindex(:), conindex(:), &
            row(:), col(:)

        ! Auxiliary variables
        real(R8), allocatable               :: valxx(:), valxy(:), valyx(:), &
            valyy(:), xv(:, :), yv(:, :), xf(:), yf(:), gxf(:), gyf(:), &
            dx(:), dy(:), gxxf(:), gxyf(:), gyxf(:), gyyf(:), gxxxf(:), &
            gxyxf(:), gyxxf(:), gyyxf(:), gxxyf(:), gxyyf(:), gyxyf(:), &
            gyyyf(:)
        
        ! Initialize
        !===========
        ! Checks
        if ( (.not. allocated(lambda)) .and. dohessian) then
            ! Throw error
            call gdErrorHandler('When evaluating the hessian vector' &
                // ' multiplication, lambda must be given')
        end if

        if (size(lambda) .ne. constraints%ncon) then
            ! Lambda should have the same size as the constraints
            call gdErrorHandler('Lambda should have the same size ' &
                // 'as the constraint vector')
        end if

        ! Counters
        ic = 0 ! constraint counter (local)
        ivg = 0 ! value index for gradient
        ivh = 0 ! value index for hessian

        ! Associate
        associate(&
            interp  => magneticField%interp,    &
            nc      => constraints%ncon,        &
            ev      => constraints%edgevert,    &
            vert    => grid%vert,               &
            nv      => grid%vert%ntot,          &
            x       => grid%vert%x,             & 
            y       => grid%vert%y              & 
            )

        ! Constraint value
        !=================
        ! Allocate
        if (.not. allocated(G)) then 
            allocate(G(nc))
        else
            if (size(G) .ne. nc) then 

                ! Print a warning and reallocate
                print *, 'EvaluateEdgelengthsConstraints: ' &
                    // 'Wrong dimension of G, reallocating'
                
                ! Deallocate and reallocate
                deallocate(G)
                allocate(G(nc))

            end if
        end if

        ! Precompute
        allocate(xv(nc, 2), yv(nc, 2), xf(nc), yf(nc), gxf(nc), gyf(nc), &
            dx(nc), dy(nc))

        xv(:, 1) = x(ev(:, 1))
        xv(:, 2) = x(ev(:, 2)) 
        yv(:, 1) = y(ev(:, 1))
        yv(:, 2) = y(ev(:, 2))
        dx = xv(:, 2) - xv(:, 1)
        dy = yv(:, 2) - yv(:, 1)
        xf = 0.5*sum(xv, 2)
        yf = 0.5*sum(yv, 2) 
        call EvaluateBicubicSplineInterpolant(xf, yf, gxf, &
            interp, '0', '1')
        call EvaluateBicubicSplineInterpolant(xf, yf, gyf, &
            interp, '1', '0')
        gxf = -gxf ! take correct sign

        ! Evaluate
        G = gxf*dx + gyf*dy

        ! Constraint gradient
        !====================
        if (dogradient) then 
            ! Initialize
            jacG%nrow = nc 
            jacG%ncol = designvariables%nphi

            ! Check design variables
            select case(designvariables%type)

            case ('coordinates')

                ! Order in jacobian: first x, then y. 

                ! Allocate
                jacG%nval = 4*nc
                call jacG%Allocate() 
                allocate(conindex(nc))
                allocate(valindex(nc))

                ! Precompute
                allocate(gxxf(nc), gxyf(nc), gyxf(nc), gyyf(nc))
                call EvaluateBicubicSplineInterpolant(xf, yf, gxxf, &
                    interp, '1', '1')
                call EvaluateBicubicSplineInterpolant(xf, yf, gxyf, &
                    interp, '0', '2')
                call EvaluateBicubicSplineInterpolant(xf, yf, gyxf, &
                    interp, '1', '1')
                call EvaluateBicubicSplineInterpolant(xf, yf, gyyf, &
                    interp, '2', '0')
                gxxf = -gxxf ! correct sign
                gxyf = -gxyf ! correct sign

                ! Build constraint indices
                conindex = [(k, k = ic+1, ic+nc)]

                ! x-contribution
                !---------------
                ! Build indices for xv1
                valindex = [(k, k = ivg+1, ivg+nc)]

                ! Add values
                jacG%row(valindex) = conindex  
                jacG%col(valindex) = ev(:, 1)
                jacG%val(valindex) = (-gxf + 0.5*dx*gxxf + 0.5*dy*gyxf)

                ! Update counters
                ivg = ivg + nc

                ! Build indices for xv2
                valindex = [(k, k = ivg+1, ivg+nc)]

                ! Add values
                jacG%row(valindex) = conindex  
                jacG%col(valindex) = ev(:, 2)
                jacG%val(valindex) = (gxf + 0.5*dx*gxxf + 0.5*dy*gyxf)

                ! Update counters
                ivg = ivg + nc

                ! y-contribution
                !---------------
                ! Build indices for yv1
                valindex = [(k, k = ivg+1, ivg+nc)]

                ! Add values
                jacG%row(valindex) = conindex 
                jacG%col(valindex) = ev(:, 1) + grid%vert%ntot 
                jacG%val(valindex) = (-gyf + 0.5*dy*gyyf + 0.5*dx*gxyf)

                ! Update counters
                ivg = ivg + nc

                ! Build indices for yv2
                valindex = [(k, k = ivg+1, ivg+nc)]

                ! Add values
                jacG%row(valindex) = conindex 
                jacG%col(valindex) = ev(:, 2) + grid%vert%ntot 
                jacG%val(valindex) = (gyf + 0.5*dy*gyyf + 0.5*dx*gxyf)

                ! Update counters
                ivg = ivg + nc

                ! Build gradient
                !---------------
                gradG%nrow = jacG%ncol 
                gradG%ncol = jacG%nrow 
                gradG%nval = jacG%nval 
                
                call gradG%Allocate()
                gradG%row = jacG%col 
                gradG%col = jacG%row
                gradG%val = jacG%val

                ! Housekeeping
                call jacG%Deallocate()

            case default

                ! Unknown, throw error
                call gdErrorHandler('Gradient not implemented for ' &
                    // 'this type of design variable')

            end select

        end if

        ! Constraint hessian
        !===================
        if (dohessian) then 

            ! Initialize
            ic = 0
            hessG%nrow = designvariables%nphi 
            hessG%ncol = designvariables%nphi 

            ! Check design variables
            select case(designvariables%type)

            case ('coordinates')
            
                ! Initialize
                !===========
                ! Allocate
                hessG%nval = 16*nc 
                if (.not. allocated(valindex)) then
                    allocate(valindex(nc))
                end if
                if (.not. allocated(conindex)) then 
                    allocate(conindex(nc))
                end if 
                if (.not. allocated(hessG%val)) then
                    call hessG%Allocate()
                end if
                allocate(valxx(4*nc))
                allocate(valxy(4*nc))
                allocate(valyx(4*nc))
                allocate(valyy(4*nc))
                allocate(row(4*nc), col(4*nc))

                ! Precompute
                allocate(gxxxf(nc), gxyxf(nc), gyxxf(nc), gyyxf(nc), &
                    gxxyf(nc), gxyyf(nc), gyxyf(nc), gyyyf(nc))

                call EvaluateBicubicSplineInterpolant(xf, yf, gxxxf, &
                    interp, '2', '1')
                call EvaluateBicubicSplineInterpolant(xf, yf, gxyxf, &
                    interp, '1', '2')
                call EvaluateBicubicSplineInterpolant(xf, yf, gyxxf, &
                    interp, '2', '1')
                call EvaluateBicubicSplineInterpolant(xf, yf, gyyxf, &
                    interp, '3', '0') 
                call EvaluateBicubicSplineInterpolant(xf, yf, gxxyf, &
                    interp, '1', '2')
                call EvaluateBicubicSplineInterpolant(xf, yf, gxyyf, &
                    interp, '0', '3')
                call EvaluateBicubicSplineInterpolant(xf, yf, gyxyf, &
                    interp, '1', '2')
                call EvaluateBicubicSplineInterpolant(xf, yf, gyyyf, &
                    interp, '2', '1')
                gxxxf = -gxxxf 
                gxyxf = -gxyxf 
                gxxyf = -gxxyf 
                gxyyf = -gxyyf 

                ! Compute contributions
                !======================

                ! Build constraint indices
                conindex = [(k, k = ic+1, ic+nc)]

                ! v1, v1
                valindex = [(k, k = ivh+1, ivh+nc)] 
                row(valindex) = ev(:, 1)
                col(valindex) = ev(:, 1)
                valxx(valindex) = (0.25*dx*gxxxf - 1.0*gxxf + 0.25*dy*gyxxf)*lambda ! x1x1
                valxy(valindex) = (0.25*dx*gxyxf - 0.5*gyxf - 0.5*gxyf + 0.25*dy*gyyxf)*lambda ! x1y1
                valyx(valindex) = (0.25*dx*gxyxf - 0.5*gyxf - 0.5*gxyf + 0.25*dy*gyyxf)*lambda ! y1x1
                valyy(valindex) = (0.25*dx*gxyyf - 1.0*gyyf + 0.25*dy*gyyyf)*lambda ! y1y1
                ivh = ivh + nc 

                ! v1, v2
                valindex = [(k, k = ivh+1, ivh+nc)] 
                row(valindex) = ev(:, 1)
                col(valindex) = ev(:, 2)
                valxx(valindex) = (0.25*dx*gxxxf + 0.25*dy*gyxxf)*lambda ! x1x2
                valxy(valindex) = (0.5*gyxf - 0.5*gxyf + 0.25*dx*gxyxf + 0.25*dy*gyyxf)*lambda ! x1y2
                valyx(valindex) = (0.5*gxyf - 0.5*gyxf + 0.25*dx*gxyxf + 0.25*dy*gyyxf)*lambda ! y1x2
                valyy(valindex) = (0.25*dx*gxyyf + 0.25*dy*gyyyf)*lambda ! y1y2
                ivh = ivh + nc 

                ! v2, v1
                valindex = [(k, k = ivh+1, ivh+nc)] 
                row(valindex) = ev(:, 2)
                col(valindex) = ev(:, 1)
                valxx(valindex) = (0.25*dx*gxxxf + 0.25*dy*gyxxf)*lambda ! x2x1
                valxy(valindex) = (0.5*gxyf - 0.5*gyxf + 0.25*dx*gxyxf + 0.25*dy*gyyxf)*lambda ! x2y1
                valyx(valindex) = (0.5*gyxf - 0.5*gxyf + 0.25*dx*gxyxf + 0.25*dy*gyyxf)*lambda ! y2x1
                valyy(valindex) = (0.25*dx*gxyyf + 0.25*dy*gyyyf)*lambda ! y1y2
                ivh = ivh + nc 

                ! v2, v2
                valindex = [(k, k = ivh+1, ivh+nc)] 
                row(valindex) = ev(:, 2)
                col(valindex) = ev(:, 2)
                valxx(valindex) = (1.0*gxxf + 0.25*dx*gxxxf + 0.25*dy*gyxxf)*lambda ! x2x2
                valxy(valindex) = (0.5*gxyf + 0.5*gyxf + 0.25*dx*gxyxf + 0.25*dy*gyyxf)*lambda ! x2y2
                valyx(valindex) = (0.5*gxyf + 0.5*gyxf + 0.25*dx*gxyxf + 0.25*dy*gyyxf)*lambda ! y2x2
                valyy(valindex) = (1.0*gyyf + 0.25*dx*gxyyf + 0.25*dy*gyyyf)*lambda ! y2y2
                ivh = ivh + nc 

                ! Build full hessian
                !===================
                hessG%row = [row, row, row+vert%ntot, row+vert%ntot]
                hessG%col = [col, col+vert%ntot, col, col+vert%ntot]
                hessG%val = [valxx, valxy, valyx, valyy]

            case default

                ! Unknown, throw error
                call gdErrorHandler('Gradient not implemented for ' &
                    // 'this type of design variable')

            end select

        end if
        
        ! Housekeeping
        !=============
        ! End associate
        end associate

        ! Deallocate
        deallocate(gxf, gyf, xv, yv, xf, yf, dx, dy)

        if (dogradient) then 
            deallocate(valindex, conindex)
            deallocate(gxxf, gxyf, gyxf, gyyf)
        end if 

        if (dohessian) then 
            if (allocated(valindex)) then 
                deallocate(valindex, conindex)
            end if 
            deallocate(valxx, valxy, valyx, valyy)
            deallocate(gxxxf, gxyxf, gyxxf, gyyxf, gxxyf, gxyyf, &
                gyxyf, gyyyf)  
        end if

    end subroutine


end module
