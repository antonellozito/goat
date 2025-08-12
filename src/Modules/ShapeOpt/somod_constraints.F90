!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains the constraint classes specific for shape 
! optimization

! The constraints are structured as follows:
! - All derived constraint types inherit from the 'mother' type 
!   'GenericConstraintsSOUDT', which contains the field 'ncon', the 
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

module somod_constraints
    
    ! Initialize
    !============
    ! Load modules
    use optmod_constraints
    use optmod_optimizationengine
    use somod_designvariables
    use somod_userinput
    use gdmod_constraints ! to be cleaned up?
    use gdmod_optimizationengine
    use, intrinsic :: ieee_arithmetic
    use PolygonLevelsetFunction2D

    ! The usual
    implicit none
    save
    public 

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    ! Auxiliary types
    !================

    ! Abstract types
    !===============
    ! Generic constraint type
    type, abstract, extends(ConstraintsUDT) :: GenericConstraintsSOUDT

        ! Description
        !============
        ! Generic type for grid deformation constraints. Inherits from 
        ! the mother constraint type defined in gdmod_constraints.

        ! Fields:
        integer(I8)                 :: ncon = 0 ! number of constraints

    contains 

        ! Initialization
        procedure(InitializeConstraintsSOINT), deferred :: Initialize 

        ! Evaluation
        procedure(EvaluateConstraintsSOINT), deferred :: Evaluate

    end type

    ! Specific constraint types
    !==========================
    ! Fix vertices
    type, extends(GenericConstraintsSOUDT) :: FixedVesselPointsConstraintsUDT

        ! Description
        !============
        ! This constraint fixes certain vessel coordinate points to 
        ! a specific location. Only possible if design variables contain
        ! the vessel coordinates of course. The following fields are 
        ! added:
        ! 
        !   xv, yv:         allocatable arrays with coordinates of fixed
        !                   vertices
        !   vertIDs:        vessel vertex indices (local) of the 
        !                   constrained vertices

        ! Fields
        real(R8), allocatable           :: xv(:), yv(:) ! coordinates
        integer(I8), allocatable        :: vertIDs(:) ! indices 

    contains

        ! Initialization
        procedure :: Initialize     => InitializeFixedVesselPointsConstraints

        ! Evaluation
        procedure :: Evaluate       => EvaluateFixedVesselPointsConstraints

        ! Output
        procedure :: WriteData      => WriteDataFixedVesselPointsConstraints

    end type

    ! GOAT 
    type, extends(GenericConstraintsSOUDT) :: GoatConstraintsUDT 

        ! Description
        !============
        ! Constraint type to impose the goat optimization problem as 
        ! constraint to the shape optimization problem. We assume that 
        ! the goat problem, which has all the necessary evaluation 
        ! routines etc, is present in the shape optimization problem 
        ! structure and can be passed to the initialization and 
        ! evaluation routines. Therefore, no additional fields are 
        ! needed here. 

    contains 

        ! Initialization
        procedure :: Initialize     => InitializeGoatConstraints

        ! Evaluation
        procedure :: Evaluate       => EvaluateGoatConstraints

    end type

    ! Flux constraints
    type, extends(GenericConstraintsSOUDT) :: FixedVesselFluxConstraintsUDT

        ! Description
        !============
        ! Simple constraints to fix the flux values of the vessel to 
        ! the initial flux value. It is assumed no x-points or other
        ! special points are present. 

        integer(I8)                     :: nID
        integer(I8), allocatable        :: ID(:)
        real(R8), allocatable           :: PsiD(:)

    contains
    
        ! Initialization
        procedure :: Initialize     => InitializeFixedVesselFluxConstraints

        ! Evaluation
        procedure :: Evaluate       => EvaluateFixedVesselFluxConstraints

        ! Write
        procedure :: WriteData      => WriteDataFixedVesselFluxConstraints

    end type

    ! Distance constraints
    type, extends(GenericConstraintsSOUDT) :: VesselDistanceConstraintsUDT
        
        ! Description
        !============
        ! Constraints that fix the maximal vessel Distance based on a 
        ! closed polygon description. Should only be used as an 
        ! inequality constraint. The distance is based on the vessel
        ! levelset function value and should be interpreted (roughly) 
        ! as the perpendicular signed distance to the vessel. As always,
        ! the constraints can be imposed on vessel polygon vertices 
        ! only, and the vertex IDs can be given either directly or based
        ! on the original vessel structure IDs. 

        ! Note that it is possible to set upper and lower bounds by 
        ! having two instances of this constraint and constraining the 
        ! maximal distance to be zero. Furthermore, by defining
        ! two different polygons, even the shape can be free (even 
        ! disjoint polygons are allowed, as long as they are closed)

        ! Fields:
        ! - nID:    number of vertices
        ! - ID:     vessel vertex indices to be considered 
        ! - d:      maximal/minimal signed distance w.r.t. polygon 
        ! - mode:   either 'max' or 'min', to determine if distance is 
        !           maximal or minimal
        ! - ps:     polygon set describing the (set of) closed polygons
        ! - plf:    polygon levelset function based on ps for 
        !           evaluation

        integer(I8)                     :: nID
        integer(I8), allocatable        :: ID(:)
        real(R8)                        :: d 
        character(:), allocatable       :: mode 
        type(PolygonSetUDT)             :: ps 
        class(PolygonLevelsetFunction2DUDT), allocatable    :: plf

    contains

        ! Initialization
        procedure :: Initialize     => InitializeVesselDistanceConstraints

        ! Evaluation
        procedure :: Evaluate       => EvaluateVesselDistanceConstraints

        ! Write
        procedure :: WriteData      => WriteDataVesselDistanceConstraints

    end type

    ! Vessel angle constraints 
    type, extends(GenericConstraintsSOUDT) :: VesselAngleDifferenceConstraintsUDT
        
        ! Description
        !============
        ! Constraints that fix how much the angle between vessel faces
        ! can deviate from the initial angle. The upper and lower bounds
        ! on this deviation are based on user input and on the initial 
        ! angle (i.e. the minimal and maximal angle are always capped by
        ! pi-alpha and -pi+alpha, where alpha > 0 is a minimal angle that
        ! needs to be guaranteed). The constraint on the upper and lower
        ! bound are therefore formulated as:
        !
        !       theta < theta0 + dtheta_upper
        !       theta > theta0 - dtheta_lower 
        ! 
        ! where dtheta_upper, dtheta_lower > 0. These angles are:
        !
        !       dtheta_upper = dtheta_user if theta0+dtheta > pi-alpha
        !       dtheta_upper = pi - alpha otherwise
        !       (analogous for dtheta_lower)
        !
        ! In principle it's possible to determine for each possible 
        ! vertex a different maximal and minimal angle, and that the
        ! maximal and minimal angle are different. The latter, however,
        ! makes not much sense since it is a priori unknown how the 
        ! vessel polygon is actually sorted, and since this constraint
        ! should anyway reflect a minimal feature size, which shouldn't
        ! depend on orientation. The only reason we distinguish here 
        ! is that one of the two may be effectively bounded by the 
        ! minimal feature size while the other is not. 

        ! Fields:
        ! - nvpairs:    number of vertex pairs
        ! - vpairs:     vessel vertex pairs to be considered 
        ! - theta0:     initial angle
        ! - dtheta:     angle differences (upper and lower)

        integer(I8)                     :: nvpairs
        integer(I8), allocatable        :: vpairs(:, :)
        real(R8), allocatable           :: theta0(:), dtheta_ub(:), &
            dtheta_lb(:)

    contains

        ! Initialization
        procedure :: Initialize     => InitializeVesselAngleDifferenceConstraints

        ! Evaluation
        procedure :: Evaluate       => EvaluateVesselAngleDifferenceConstraints

        ! Write
        procedure :: WriteData      => WriteDataVesselAngleDifferenceConstraints

    end type


    ! Overarching types
    !==================
    ! Equality constraints
    type, extends(ConstraintsUDT) :: EqConSOUDT

        ! Description
        !============
        ! This type contains all the different constraints as different
        ! derived types. For each type of constraint, a different
        ! type is defined. 

        ! Total number of constraints 
        integer(I8)                         :: neqcon = 0

        ! Constraint switches
        logical                             :: dofixedvesselpoints = .false.
        logical                             :: dofixedvesselflux = .false.
        logical                             :: dogoat = .true.

        type(FixedVesselPointsConstraintsUDT)   :: fixedvesselpoints 
        type(FixedVesselFluxConstraintsUDT)     :: fixedvesselflux
        type(GoatConstraintsUDT)                :: goat

    contains

        ! Procedure to initialize constraints
        procedure :: Initialize         => InitializeEqConSO

        ! Procedure to evaluate constraints
        procedure :: Evaluate           => EvaluateEqConSO

    end type 

    ! Inequality constraints
    type, extends(ConstraintsUDT) :: IneqConSOUDT

        ! Total number of inequality constraints
        integer(I8)             :: nineqcon = 0

        ! Constraint switches
        logical                 :: dovesselupperbound = .false. 
        logical                 :: dovessellowerbound = .false. 
        logical                 :: dovesselangledifference = .false.

        type(VesselDistanceConstraintsUDT)  :: vesselupperbound 
        type(VesselDistanceConstraintsUDT)  :: vessellowerbound
        type(VesselAngleDifferenceConstraintsUDT)   :: vesselangledifference

    contains

        ! Constraints initialization
        procedure :: Initialize         => InitializeIneqConSO

        ! Constraints evaluation
        procedure :: Evaluate           => EvaluateIneqConSO

    end type

    ! All constraints for the grid deformation
    type, extends(ConstraintsUDT) :: ConstraintsSOUDT

        ! Description
        !============
        ! Basic constraint holder for shape optimization

        ! Fields: 

        ! Equality constraints 
        type(EqConSOUDT)        :: eqcon

        ! Inequality constraints 
        type(IneqConSOUDT)      :: ineqcon

    contains

        ! Initialization
        procedure :: Initialize         => InitializeConstraintsSO

        ! Number of constraints getter
        procedure :: GetConstraintsDimensions  => GetConstraintsDimensionsSO

        ! Evaluation
        ! procedure :: Evaluate           => EvaluateEqualityConstraintsGD
        
        ! Housekeeping

    end type

    ! Monitor
    !========
    type ConstraintsMonitorSOUDT

        ! Description
        !============
        ! This object can be used to keep track of the amount of 
        ! constraints that are imposed per design variable (or other 
        ! quantity of interest). Not always an airtight way of 
        ! preventing infeasibility, but can help to detect it 
        ! beforehand... 

        ! Fields:
        integer(I8), allocatable           :: eqvcc(:), ineqvcc(:)
        integer(I8), allocatable           :: maxeqvcc(:), maxineqvcc(:) 
        
    contains 

        ! Routines
        procedure :: Initialize         => InitializeMonitorSO
        ! procedure :: CheckConstraints   => CheckConstraintsMonitor

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
        subroutine InitializeConstraintsSOINT(constraints, goat, &
            monitor, designvariables, options)

            ! Description
            !============
            ! This routine serves as general initialization routine for
            ! a generic grid deformation constraint (that inherits from
            ! the GenericConstraintUDT type)

            ! Import
            import :: GenericConstraintsSOUDT, OptimizationProblemGDUDT, &
                ConstraintsMonitorSOUDT, ConstraintOptionsSOUDT, &
                DesignVariablesSOUDT

            ! Declare
            class(GenericConstraintsSOUDT)      :: constraints 
            type(OptimizationProblemGDUDT)      :: goat
            type(ConstraintsMonitorSOUDT)         :: monitor
            type(ConstraintOptionsSOUDT)        :: options
            class(DesignVariablesSOUDT)         :: designvariables

        end subroutine

        ! Constraint evaluation
        subroutine EvaluateConstraintsSOINT(constraints, G, gradG, & 
            hessG, goat, dogradient, dohessian, designvariables, lambda)

            ! Description
            !============
            ! This reoutine serves as a general evaluation routine for 
            ! a generic grid deformation constraint. 

            ! Import
            import :: GenericConstraintsSOUDT, MySparseUDT, &
                R8, OptimizationProblemGDUDT, DesignVariablesSOUDT
            
            ! Declare
            class(GenericConstraintsSOUDT)  :: constraints 
            real(R8), allocatable           :: G(:), lambda(:)
            type(MySparseUDT)               :: hessG, gradG 
            type(OptimizationProblemGDUDT)  :: goat
            logical                         :: dogradient, dohessian
            class(DesignVariablesSOUDT)     :: designvariables

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
    subroutine InitializeMonitorSO(monitor, goat, designvariables)

        ! Description
        !============
        ! Initializes the constraints monitor structure. It is assumed
        ! that the grid, magnetic field and environment are properly 
        ! allocated and initialized. 

        ! Declare variables
        !==================
        ! Arguments
        class(ConstraintsMonitorSOUDT)      :: monitor 
        type(OptimizationProblemGDUDT)      :: goat 
        class(DesignVariablesSOUDT)         :: designvariables

        ! Loop variables

        ! Auxiliary 
        integer(I8)                             :: nv
        real(R8), allocatable, dimension(:)     :: xv, yv

        ! Initialize
        !===========
        ! Get vessel coordinates
        call goat%environment%vessel%GetVesselCoordinates(xv, yv)
        nv = size(xv)

        ! Allocate
        select type (designvariables)

        type is (DesignVariablesVesselCoordinatesUDT)

            allocate(monitor%eqvcc(nv), monitor%maxeqvcc(nv))
            allocate(monitor%ineqvcc(nv), monitor%maxineqvcc(nv))

        type is (DesignVariablesVesselCoordinatesGoatUDT)

            allocate(monitor%eqvcc(nv), monitor%maxeqvcc(nv))
            allocate(monitor%ineqvcc(nv), monitor%maxineqvcc(nv))

        class default 

            call gdErrorHandler('InitializeMonitorSO: unknown design variable type')
            
        end select

        ! Initialize
        monitor%eqvcc(:)        = 0
        monitor%ineqvcc(:)      = 0
        monitor%maxeqvcc(:)     = 2 
        monitor%maxineqvcc(:)   = 1000 ! a stupid large number - can impose any number

    end subroutine

    !------------------------------------------------------------------!
    !                           GENERAL CONSTRAINTS                    !
    !------------------------------------------------------------------!
    ! Initialization
    subroutine InitializeConstraintsSO(constraints, goat, designvariables, options)

        ! Description
        !============
        ! Routine that initializes the equality and inequality 
        ! constraints, using the initialization routines of those 
        ! objects. 

        ! Declare variables
        !==================
        ! Arguments
        class(ConstraintsSOUDT)         :: constraints 
        type(OptimizationProblemGDUDT)  :: goat
        type(ConstraintOptionsSOUDT)    :: options
        type(ConstraintsMonitorSOUDT)   :: monitor
        class(DesignVariablesSOUDT)     :: designvariables

        ! Loop variables

        ! Auxiliary variables

        ! Initialize monitor
        !===================
        call monitor%Initialize(goat, designvariables)

        ! Initialize constraints
        !=======================
        ! Equality constraints
        call constraints%eqcon%Initialize(goat, options, &
            designvariables, monitor)

        ! Inequality constraints
        call constraints%ineqcon%Initialize(goat, options, &
            designvariables, monitor)

    end subroutine

    ! Dimension getter
    subroutine GetConstraintsDimensionsSO(constraints, neqcon, nineqcon)

        ! Description
        !============
        ! Return the current dimensions of the equality and inequality
        ! constraints. Can be used for initialization of other 
        ! quantities at higher levels. 

        ! Declare variables
        !==================
        ! Arguments
        class(ConstraintsSOUDT)         :: constraints 
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
    subroutine InitializeEqConSO(constraints, goat, constraintoptions, &
        designvariables, monitor)

        ! Description
        !============
        ! Routine that initializes the desired constraints. 

        ! Declare variables
        !==================
        ! Arguments
        class(EqConSOUDT)               :: constraints 
        type(OptimizationProblemGDUDT)  :: goat
        type(ConstraintOptionsSOUDT)    :: constraintoptions
        type(ConstraintsMonitorSOUDT)   :: monitor
        class(DesignVariablesSOUDT)     :: designvariables

        ! Loop variables

        ! Auxiliary variables
        
        ! Initialize constraints
        !=======================
        constraints%neqcon = 0

        ! Fixed vessel points
        if (constraintoptions%fixedvesselpoints == 1) then 

            ! Set the logical
            constraints%dofixedvesselpoints = .true. 

            ! Initialize
            call constraints%fixedvesselpoints%Initialize(goat, &
                monitor, designvariables, constraintoptions)

            ! Add constraints number
            constraints%neqcon = constraints%neqcon + &
                constraints%fixedvesselpoints%ncon 

            ! Print
            print *, 'number of fixed vessel points constraints: ', &
                constraints%fixedvesselpoints%ncon

        else

            ! Set to false, don't initialize
            constraints%dofixedvesselpoints = .false.

        end if

        ! Fixed vessel flux
        if (constraintoptions%fixedvesselflux == 1) then 

            ! Set the logical
            constraints%dofixedvesselflux = .true. 

            ! Initialize
            call constraints%fixedvesselflux%Initialize(goat, &
                monitor, designvariables, constraintoptions)

            ! Add constraints number
            constraints%neqcon = constraints%neqcon + &
                constraints%fixedvesselflux%ncon 

            ! Print
            print *, 'number of fixed vessel points constraints: ', &
                constraints%fixedvesselflux%ncon

        else

            ! Set to false, don't initialize
            constraints%dofixedvesselflux = .false.

        end if 

        ! Goat
        if (constraintoptions%goat == 1) then 

            ! Set the logical
            constraints%dogoat = .true. 

            ! Initialize
            call constraints%goat%Initialize(goat, &
                monitor, designvariables, constraintoptions)

            ! Add constraints number
            constraints%neqcon = constraints%neqcon + &
                constraints%goat%ncon 

            ! Print
            print *, 'number of goat constraints: ', &
                constraints%goat%ncon

        else

            ! Set to false, don't initialize
            constraints%dogoat = .false.

        end if 


    end subroutine

    ! Constraint evaluation
    subroutine EvaluateEqConSO(constraints, G, gradG, hessG, &
        goat, dogradient, dohessian, designvariables, lambda)

        ! Description
        !============
        ! This routine evaluates the constraints G and the corresponding
        ! gradient and hessian. To do so, every type of constraint is 
        ! checked whether it is imposed, and the contributions are 
        ! added by calling the evaluation routine of each constraint. 

        ! Declare variables
        !==================
        ! Arguments
        class(EqConSOUDT)               :: constraints
        real(R8), intent(inout)         :: G(:)
        real(R8), intent(in)            :: lambda(:)
        type(MySparseUDT)               :: gradG, hessG 
        type(OptimizationProblemGDUDT)  :: goat
        logical                         :: dogradient, dohessian 
        class(DesignVariablesSOUDT)     :: designvariables 

        ! Loop
        integer(I8)                     :: ic, k
        integer(I8), allocatable        :: conindex(:)

        ! Auxiliary
        real(R8), allocatable           :: G_fvp(:), lambda_fvp(:)
        type(MySparseUDT)               :: gradG_fvp, hessG_fvp

        real(R8), allocatable           :: G_fvf(:), lambda_fvf(:)
        type(MySparseUDT)               :: gradG_fvf, hessG_fvf

        real(R8), allocatable           :: G_goat(:), lambda_goat(:)
        type(MySparseUDT)               :: gradG_goat, hessG_goat

        ! Initialize
        !===========
        ! Set the constraint counter
        ic = 0

        ! Constraint values
        !==================
        gradG = SpZeros(designvariables%nphi, 0) ! will concatenate
        hessG = SpZeros(designvariables%nphi, designvariables%nphi) ! will add

        ! Fixed vessel points constraints
        !--------------------------------
        if (constraints%dofixedvesselpoints) then 
            ! Construct the constraint index
            conindex = [(k, k = ic+1, ic+constraints%fixedvesselpoints%ncon)]

            ! Allocate & initialize
            lambda_fvp = lambda(conindex)

            ! Call the evaluation routine
            call constraints%fixedvesselpoints%Evaluate(G_fvp, &
                gradG_fvp, hessG_fvp, goat, dogradient, &
                dohessian, designvariables, lambda_fvp)

            ! Assign
            G(conindex) = G_fvp
            if (dogradient) then 
                gradG = gradG%Concatenate(gradG_fvp, 2)
            end if
            if (dohessian) then 
                hessG = hessG + hessG_fvp
            end if 

            ! Update the constraint counter
            ic = ic + constraints%fixedvesselpoints%ncon

        end if

        ! Fixed vessel flux constraints
        !--------------------------------
        if (constraints%dofixedvesselflux) then 
            ! Construct the constraint index
            conindex = [(k, k = ic+1, ic+constraints%fixedvesselflux%ncon)]

            ! Allocate & initialize
            lambda_fvf = lambda(conindex)

            ! Call the evaluation routine
            call constraints%fixedvesselflux%Evaluate(G_fvf, &
                gradG_fvf, hessG_fvf, goat, dogradient, &
                dohessian, designvariables, lambda_fvf)

            ! Assign
            G(conindex) = G_fvf
            if (dogradient) then 
                gradG = gradG%Concatenate(gradG_fvf, 2)
            end if
            if (dohessian) then 
                hessG = hessG + hessG_fvf
            end if 

            ! Update the constraint counter
            ic = ic + constraints%fixedvesselflux%ncon

        end if

        ! GOAT constraints
        !-----------------
        if (constraints%dogoat) then 
            ! Construct the constraint index
            conindex = [(k, k = ic+1, ic+constraints%goat%ncon)]

            ! Allocate & initialize
            lambda_goat = lambda(conindex)

            ! Call the evaluation routine
            call constraints%goat%Evaluate(G_goat, &
                gradG_goat, hessG_goat, goat, dogradient, &
                dohessian, designvariables, lambda_goat)

            ! Assign
            G(conindex) = G_goat
            if (dogradient) then 
                gradG = gradG%Concatenate(gradG_goat, 2)
            end if
            if (dohessian) then 
                hessG = hessG + hessG_goat
            end if 

            ! Update the constraint counter
            ic = ic + constraints%goat%ncon
        end if 

    end subroutine

    !------------------------------------------------------------------!
    !                          INEQUALITY CONSTRAINTS                  !
    !------------------------------------------------------------------!
    ! Initialization
    subroutine InitializeIneqConSO(constraints, goat, constraintoptions, &
        designvariables, monitor)

        ! Description
        !============
        ! Routine that initializes the desired constraints

        ! Declare variables
        !==================
        ! Arguments
        class(IneqConSOUDT)             :: constraints 
        type(OptimizationProblemGDUDT)  :: goat
        type(ConstraintOptionsSOUDT)    :: constraintoptions
        type(ConstraintsMonitorSOUDT)   :: monitor
        class(DesignVariablesSOUDT)     :: designvariables

        ! Loop variables

        ! Auxiliary variables

        ! Initialize
        !===========

        ! Initialize constraints
        !=======================
        constraints%nineqcon = 0

        ! Vessel upper bound
        if (constraintoptions%vesselupperbound == 1) then 

            ! Set the logical
            constraints%dovesselupperbound = .true. 

            ! Initialize
            constraintoptions%vdoptions = constraintoptions%vduboptions 
            call constraints%vesselupperbound%Initialize(goat, &
                monitor, designvariables, constraintoptions)

            ! Add constraints number
            constraints%nineqcon = constraints%nineqcon + &
                constraints%vesselupperbound%ncon 

            ! Print
            print *, 'number of vessel upperbound constraints: ', &
                constraints%vesselupperbound%ncon

        else

            ! Set to false, don't initialize
            constraints%dovesselupperbound = .false.

        end if

        ! Vessel lower bound
        if (constraintoptions%vessellowerbound == 1) then 

            ! Set the logical
            constraints%dovessellowerbound = .true. 

            ! Initialize
            constraintoptions%vdoptions = constraintoptions%vdlboptions
            call constraints%vessellowerbound%Initialize(goat, &
                monitor, designvariables, constraintoptions)

            ! Add constraints number
            constraints%nineqcon = constraints%nineqcon + &
                constraints%vessellowerbound%ncon 

            ! Print
            print *, 'number of vessel upperbound constraints: ', &
                constraints%vessellowerbound%ncon

        else

            ! Set to false, don't initialize
            constraints%dovessellowerbound = .false.

        end if

        ! Vessel angle difference
        if (constraintoptions%vesselangledifference == 1) then 

            ! Set the logical
            constraints%dovesselangledifference = .true. 

            ! Initialize
            call constraints%vesselangledifference%Initialize(goat, &
                monitor, designvariables, constraintoptions)

            ! Add constraints number
            constraints%nineqcon = constraints%nineqcon + &
                constraints%vesselangledifference%ncon 

            ! Print
            print *, 'number of vessel angle difference constraints: ', &
                constraints%vesselangledifference%ncon

        else

            ! Set to false, don't initialize
            constraints%dovesselangledifference = .false.

        end if


    end subroutine

    ! Constraint evaluation
    subroutine EvaluateIneqConSO(constraints, G, gradG, hessG, &
        goat, dogradient, dohessian, designvariables, lambda)

        ! Description
        !============
        ! This routine evaluates the constraints G and the corresponding
        ! gradient and hessian. To do so, every type of constraint is 
        ! checked whether it is imposed, and the contributions are 
        ! added by calling the evaluation routine of each constraint. 

        ! Declare variables
        !==================
        ! Arguments
        class(IneqConSOUDT)             :: constraints
        real(R8), intent(inout)         :: G(:)
        real(R8), intent(in)            :: lambda(:)
        type(MySparseUDT)               :: gradG, hessG 
        type(OptimizationProblemGDUDT)  :: goat
        logical                         :: dogradient, dohessian 
        class(DesignVariablesSOUDT)     :: designvariables 

        ! Loop
        integer(I8)                     :: ic, k

        ! Auxiliary
        integer(I8), allocatable        :: conindex(:)
        real(R8), allocatable           :: G_vub(:), lambda_vub(:)
        type(MySparseUDT)               :: gradG_vub, hessG_vub

        real(R8), allocatable           :: G_vlb(:), lambda_vlb(:)
        type(MySparseUDT)               :: gradG_vlb, hessG_vlb

        real(R8), allocatable           :: G_vad(:), lambda_vad(:)
        type(MySparseUDT)               :: gradG_vad, hessG_vad

        ! Initialize
        !===========
        ! Set the constraint counter
        ic = 0

        ! Constraint values
        !==================
        gradG = SpZeros(designvariables%nphi, 0) ! will concatenate
        hessG = SpZeros(designvariables%nphi, designvariables%nphi) ! will add

        ! Vessel upper bound
        !-------------------
        if (constraints%dovesselupperbound) then 
            ! Construct the constraint index
            conindex = [(k, k = ic+1, ic+constraints%vesselupperbound%ncon)]

            ! Allocate & initialize
            lambda_vub = lambda(conindex)

            ! Call the evaluation routine
            call constraints%vesselupperbound%Evaluate(G_vub, &
                gradG_vub, hessG_vub, goat, dogradient, &
                dohessian, designvariables, lambda_vub)

            ! Assign
            G(conindex) = G_vub
            if (dogradient) then 
                gradG = gradG%Concatenate(gradG_vub, 2)
            end if
            if (dohessian) then 
                hessG = hessG + hessG_vub
            end if 

            ! Update the constraint counter
            ic = ic + constraints%vesselupperbound%ncon

        end if

        ! Vessel lower bound
        !-------------------
        if (constraints%dovessellowerbound) then 
            ! Construct the constraint index
            conindex = [(k, k = ic+1, ic+constraints%vessellowerbound%ncon)]

            ! Allocate & initialize
            lambda_vlb = lambda(conindex)

            ! Call the evaluation routine
            call constraints%vessellowerbound%Evaluate(G_vlb, &
                gradG_vlb, hessG_vlb, goat, dogradient, &
                dohessian, designvariables, lambda_vlb)

            ! Assign
            G(conindex) = G_vlb
            if (dogradient) then 
                gradG = gradG%Concatenate(gradG_vlb, 2)
            end if
            if (dohessian) then 
                hessG = hessG + hessG_vlb
            end if 

            ! Update the constraint counter
            ic = ic + constraints%vessellowerbound%ncon

        end if

        ! Vessel angle difference
        !------------------------
        if (constraints%dovesselangledifference) then 
            ! Construct the constraint index
            conindex = [(k, k = ic+1, ic+constraints%vesselangledifference%ncon)]

            ! Allocate & initialize
            lambda_vad = lambda(conindex)

            ! Call the evaluation routine
            call constraints%vesselangledifference%Evaluate(G_vad, &
                gradG_vad, hessG_vad, goat, dogradient, &
                dohessian, designvariables, lambda_vad)

            ! Assign
            G(conindex) = G_vad
            if (dogradient) then 
                gradG = gradG%Concatenate(gradG_vad, 2)
            end if
            if (dohessian) then 
                hessG = hessG + hessG_vad
            end if 

            ! Update the constraint counter
            ic = ic + constraints%vesselangledifference%ncon

        end if
        
    end subroutine

    !------------------------------------------------------------------!
    !                        FIXED VESSEL POINTS                       !
    !------------------------------------------------------------------!
    ! Initialize
    subroutine InitializeFixedVesselPointsConstraints(constraints, goat, &
        monitor, designvariables, options)

        ! Description
        !============
        ! Initialize the required fields of the fixed vessel point 
        ! constraints. Vertices that should be fixed can be defined
        ! either through the structure ID or through the unique vertex
        ! ID. For the structures, if any of the vertex labels is equal, 
        ! the vertex is constrained. In case specific vertices should be
        ! excluded from the constraint, this can be done by specifying 
        ! the unique vertex ID. It is possible that certain vertex IDs
        ! won't be available anymore, and if any of these IDs are 
        ! specified, they are ignored (a message will be printed).

        ! Notes
        !======
        ! Note 1: Currently, the implementation constrains the vertices
        ! to the initial position, though in principle this position
        ! can be an input parameter. Might be extended in the future. 
        
        ! Initialize
        !===========
        implicit none
        
        ! Declare variables
        !==================
        ! Arguments 
        class(FixedVesselPointsConstraintsUDT)  :: constraints 
        type(OptimizationProblemGDUDT)          :: goat 
        type(ConstraintsMonitorSOUDT)           :: monitor
        type(ConstraintOptionsSOUDT)            :: options 
        class(DesignVariablesSOUDT)             :: designvariables

        ! Auxiliary variables
        logical, allocatable                    :: isconstrained(:), &
            isvertextreated(:)
        integer(I8), allocatable                :: labels(:, :)

        real(R8), allocatable                   :: xv(:), yv(:)

        ! Loop variables
        integer(I8)                 :: i, k

        ! Initialize
        !===========
        ! Number of constraints
        constraints%ncon = 0

        ! Check design variable type
        select case (designvariables%type)

        case ('vesselcoordinates', 'vesselcoordinates_goat')

            ! All good

        case default

            ! All bad
            call gdErrorHandler('InitializeFixedVesselPointsConstraints: ' // &
                'design variable type not implemented')

        end select

        ! Associate
        associate(&
            ps      => goat%environment%vessel%polygonset,          &
            cc      => monitor%eqvcc,                               &
            maxcc   => monitor%maxeqvcc,                            &
            opt     => options%fvpoptions                           &
            )

        ! Get vessel coordinates and labels (should be same as design variables)
        call ps%GetLabels(labels)
        call ps%GetVertices(xv, yv)

        ! Initialize logical indicating if vertex is constrained
        allocate(isconstrained(size(labels, 1)))
        isconstrained = .false. 

        ! Initialize logical indicating if desired constrained vertex
        ! was present and constrained
        allocate(isvertextreated(size(opt%vertIDs)))
        isvertextreated = .false. 

        ! Determine constrained vertices
        !===============================
        ! Constrain per vessel structure (label 1 and 2)
        do i = 1, size(opt%structureIDs)
            ! Unpack ID
            associate(tID       => opt%structureIDs(i))

            ! Check vertices
            where ( (labels(:, 1) == tID) .or. (labels(:, 2) == tID) ) &
                isconstrained = .true. 

            ! Housekeeping
            end associate
        end do

        ! Constrain per vertex ID
        do i = 1, size(opt%vertIDs)
            ! Unpack ID
            associate(tID       => opt%vertIDs(i))

            ! Check vertices
            where( (labels(:, 3) == tID)) isconstrained = .true. 

            ! Check if found
            if (any(labels(:, 3) == tID)) then 
                isvertextreated(i) = .true. 
            end if

            ! Housekeeping
            end associate
        end do

        ! Check if we can constrain, set to false if not the case
        where (cc + 2 > maxcc) isconstrained = .false. 

        ! Add to constraints
        !===================
        ! Vertices, but locally indexed (i.e. not vertex ID, but ID 
        ! according to polygon structure)
        constraints%vertIDs = pack([(k, k = 1, size(labels, 1))], isconstrained)

        ! Get current coordinates
        constraints%xv = pack(xv, isconstrained)
        constraints%yv = pack(yv, isconstrained)

        ! set number of constraints
        constraints%ncon = 2*size(constraints%vertIDs)

        ! Update
        cc(constraints%vertIDs) = cc(constraints%vertIDs) + 2

        ! Output
        !=======
        ! Display warning message if necessary
        if (any(.not. isvertextreated)) then 
            ! Display message
            print *, 'Vertices with following IDs were not present and ' // &
                'are not constrained: ', pack(opt%vertIDs, .not. isvertextreated)
        end if

        ! Print data
        call constraints%WriteData()
        
        ! Housekeeping
        !=============
        end associate
        
    end subroutine

    ! Evaluation
    subroutine EvaluateFixedVesselPointsConstraints(constraints, G, gradG, & 
        hessG, goat, dogradient, dohessian, designvariables, lambda)

        ! Description
        !============
        ! Evaluate the fixed vessel points constraints. These are
        ! imposed by simply imposing the coordinates:
        !
        !   Gx      = xv - xv0,
        !   Gy      = yv - yv0.
        ! 
        ! Gradients and hessians w.r.t. xv, yv are simply one and zero,
        ! resp.  

        ! Initialize
        !===========
        ! Modules
        
        ! Declare variables
        !==================
        ! Arguments 
        class(FixedVesselPointsConstraintsUDT)  :: constraints 
        real(R8), allocatable                   :: G(:) 
        real(R8), allocatable                   :: lambda(:)
        type(MySparseUDT)                       :: hessG, gradG
        type(OptimizationProblemGDUDT)          :: goat 
        logical                                 :: dogradient, dohessian
        class(DesignVariablesSOUDT)             :: designvariables   
        
        ! Auxiliary
        integer(I8)                             :: ntv 
        integer(I8), allocatable, dimension(:)  :: row, col
        real(R8), allocatable, dimension(:)     :: xv, yv, val
        type(MySparseUDT)                       :: jacG

        ! Loop 
        integer(I8)                         :: ic, ivg, ivh, k

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
        if (allocated(G)) then 
            if (size(G, 1) .ne. constraints%ncon) then 
                deallocate(G)
                allocate(G(constraints%ncon))
            end if 
        else
            allocate(G(constraints%ncon))
        end if 

        ! Counters
        ic = 0 ! constraint counter (local)
        ivg = 0 ! value index for gradient
        ivh = 0 ! value index for hessian

        ! Associate
        associate(&
            ps          => goat%environment%vessel%polygonset,    &
            tv          => constraints%vertIDs,           &
            xv0         => constraints%xv,                &
            yv0         => constraints%yv                 &
            )

        ! Get coordinates
        call ps%GetVertices(xv, yv)

        ! Get number of constrained vertices
        ntv = size(tv)

        ! Evaluate
        !=========
        ! x-coordinate
        G(ic+1:ic+ntv) = xv(tv) - xv0 
        ic = ic + ntv 

        ! y-coordinate
        G(ic+1:ic+ntv) = yv(tv) - yv0 
        ic = ic + ntv 

        ! Derivatives
        !============
        ! Initialize
        jacG%nrow = constraints%ncon 
        jacG%ncol = designvariables%nphi 
        
        hessG%nrow = designvariables%nphi 
        hessG%ncol = designvariables%nphi 

        ! Check design variable type
        select case(designvariables%type)

        case ('vesselcoordinates', 'vesselcoordinates_goat') ! same since vessel coordinates come first in design variable vector

            ! Gradient
            !---------
            if (dogradient) then 
                ! Set contributions
                row = [(k, k = 1, constraints%ncon)]
                col = [tv, tv+size(xv)]
                val = spread(1.0, 1, 2*size(tv))
                jacG = ConstructMySparse(row, col, val, constraints%ncon, &
                    designvariables%nphi)
                
                ! Transpose
                gradG = jacG%Transpose()
            end if 

            ! Hessian
            !--------
            if (dohessian) then 
                ! Simply zero
                hessG = SpZeros(designvariables%nphi, designvariables%nphi)
            end if 

        case default

            call gdErrorHandler('EvaluateFixedVesselPointsConstraints: ' // &
                'derivatives not implemented for design variable type: ' // &
                designvariables%type)

        end select
        
        ! Housekeeping
        !=============
        ! End associate
        end associate

    end subroutine 

    ! Data output
    subroutine WriteDataFixedVesselPointsConstraints(constraints)

        ! Description
        !============
        ! Write vessel nodes in the following format:
        ! ID, x, y 
        ! Different files are written for special vertices, fixed
        ! vertices, flux surface vertices, and tangency points

        ! The usual
        implicit none

        ! Declare variables
        !==================
        ! Arguments
        class(FixedVesselPointsConstraintsUDT)          :: constraints

        ! Auxiliary
        character(:), allocatable                       :: filepath

        ! Initialize
        !===========
        ! Set the correct directories
        allocate(character(len('so_con_fvp_vertices')) :: filepath)
        filepath = 'so_con_fvp_vertices'

        ! Associate
        associate(&
            x               => constraints%xv,     &
            y               => constraints%yv,     &
            ID              => constraints%vertIDs  &
            )

        ! Write
        call WriteVertexData(ID, x, y, filepath)

        ! Housekeeping
        !=============
        end associate

    end subroutine

    !------------------------------------------------------------------!
    !                         FIXED VESSEL FLUX                        !
    !------------------------------------------------------------------!
    ! Initialize
    subroutine InitializeFixedVesselFluxConstraints(constraints, goat, &
        monitor, designvariables, options)

        ! Description
        !============
        ! Initialize the fixed vessel flux constraints. 

        ! Notes
        !======
        
        ! Initialize
        !===========
        implicit none
        
        ! Declare variables
        !==================
        ! Arguments 
        class(FixedVesselFluxConstraintsUDT)    :: constraints 
        type(OptimizationProblemGDUDT)          :: goat 
        type(ConstraintsMonitorSOUDT)           :: monitor
        type(ConstraintOptionsSOUDT)            :: options 
        class(DesignVariablesSOUDT)             :: designvariables

        ! Auxiliary variables
        logical, allocatable                    :: isconstrained(:), &
            isvertextreated(:)
        integer(I8), allocatable                :: labels(:, :)

        real(R8), allocatable                   :: xv(:), yv(:)

        ! Loop variables
        integer(I8)                 :: i, k

        ! Initialize
        !===========
        ! Number of constraints
        constraints%ncon = 0

        ! Check design variable type
        select case (designvariables%type)

        case ('vesselcoordinates', 'vesselcoordinates_goat')

            ! All good

        case default

            ! All bad
            call gdErrorHandler('InitializeFixedVesselFluxConstraints: ' // &
                'design variable type not implemented')

        end select

        ! Associate
        associate(&
            ps      => goat%environment%vessel%polygonset,          &
            cc      => monitor%eqvcc,                               &
            maxcc   => monitor%maxeqvcc,                            &
            opt     => options%fvfoptions                           &
            )

        ! Get vessel coordinates and labels (should be same as design variables)
        call ps%GetLabels(labels)
        call ps%GetVertices(xv, yv)

        ! Initialize logical indicating if vertex is constrained
        allocate(isconstrained(size(labels, 1)))
        isconstrained = .false. 

        ! Initialize logical indicating if desired constrained vertex
        ! was present and constrained
        allocate(isvertextreated(size(opt%vertIDs)))
        isvertextreated = .false. 

        ! Determine constrained vertices
        !===============================
        ! Constrain per vessel structure (label 1 and 2)
        do i = 1, size(opt%structureIDs)
            ! Unpack ID
            associate(tID       => opt%structureIDs(i))

            ! Check vertices
            where ( (labels(:, 1) == tID) .or. (labels(:, 2) == tID) ) &
                isconstrained = .true. 

            ! Housekeeping
            end associate
        end do

        ! Constrain per vertex ID
        do i = 1, size(opt%vertIDs)
            ! Unpack ID
            associate(tID       => opt%vertIDs(i))

            ! Check vertices
            where( (labels(:, 3) == tID)) isconstrained = .true. 

            ! Check if found
            if (any(labels(:, 3) == tID)) then 
                isvertextreated(i) = .true. 
            end if

            ! Housekeeping
            end associate
        end do

        ! Check if we can constrain, set to false if not the case
        where (cc + 1> maxcc) isconstrained = .false. 

        ! Add to constraints
        !===================
        ! Vertices, but locally indexed (i.e. not vertex ID, but ID 
        ! according to polygon structure)
        constraints%ID = pack([(k, k = 1, size(labels, 1))], isconstrained)
        constraints%nID = size(constraints%ID)

        ! Get current psi values
        if (allocated(constraints%PsiD)) then 
            deallocate(constraints%PsiD)
        end if 
        allocate(constraints%PsiD(size(constraints%ID)))
        call goat%magneticField%interp%Evaluate(&
            pack(xv, isconstrained), pack(yv, isconstrained), 0, 0, &
            constraints%PsiD)

        ! set number of constraints
        constraints%ncon = size(constraints%ID)

        ! Update
        cc(constraints%ID) = cc(constraints%ID) + 1

        ! Output
        !=======
        ! Display warning message if necessary
        if (any(.not. isvertextreated)) then 
            ! Display message
            print *, 'Vertices with following IDs were not present and ' // &
                'are not constrained: ', pack(opt%vertIDs, .not. isvertextreated)
        end if

        ! Print data
        call constraints%WriteData(goat)
        
        ! Housekeeping
        !=============
        end associate
        
    end subroutine

    ! Evaluation
    subroutine EvaluateFixedVesselFluxConstraints(constraints, G, gradG, & 
        hessG, goat, dogradient, dohessian, designvariables, lambda)

        ! Description
        !============
        ! Evaluate the fixed vessel flux constraints
        ! 
        !       Psi(xv, yv) - PsiD = 0
        ! 
        ! Gradients and hessians w.r.t. xv, yv are simply the 
        ! derivatives of the Psi function

        ! Initialize
        !===========
        ! Modules
        
        ! Declare variables
        !==================
        ! Arguments 
        class(FixedVesselFluxConstraintsUDT)    :: constraints 
        real(R8), allocatable                   :: G(:) 
        real(R8), allocatable                   :: lambda(:)
        type(MySparseUDT)                       :: hessG, gradG
        type(OptimizationProblemGDUDT)          :: goat 
        logical                                 :: dogradient, dohessian
        class(DesignVariablesSOUDT)             :: designvariables   
        
        ! Auxiliary
        integer(I8)                             :: nv, ntv
        integer(I8), allocatable, dimension(:)  :: row, col
        real(R8), allocatable, dimension(:)     :: xv, yv, val, psi, &
            dpsidx, dpsidy, d2psidx2, d2psidy2, d2psidxdy
        type(MySparseUDT)                       :: jacG

        ! Loop 
        integer(I8)                         :: ic, ivg, ivh, k

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
        if (allocated(G)) then 
            if (size(G, 1) .ne. constraints%ncon) then 
                deallocate(G)
                allocate(G(constraints%ncon))
            end if 
        else
            allocate(G(constraints%ncon))
        end if 

        ! Counters
        ic = 0 ! constraint counter (local)
        ivg = 0 ! value index for gradient
        ivh = 0 ! value index for hessian

        ! Associate
        associate(&
            mfinterp    => goat%magneticField%interp,     &
            tv          => constraints%ID,                &
            PsiD        => constraints%PsiD,              &
            ps          => goat%environment%vessel%polygonset   &
            )

        ! Get coordinates
        call ps%GetVertices(xv, yv)

        ! Get number of vessel vertices
        nv = size(xv)

        ! Take only coordinates that are considered
        xv = xv(tv)
        yv = yv(tv)

        ! Get number of constrained vertices
        ntv = constraints%nID 

        ! Evaluate
        !=========
        ! Psi values at vessel coordinates
        allocate(Psi(ntv))
        call mfinterp%Evaluate(xv, yv, 0, 0, Psi)

        ! Constraint
        G(ic+1:ic+ntv) = Psi - PsiD 
        
        ! Derivatives
        !============
        ! Initialize
        jacG%nrow = constraints%ncon 
        jacG%ncol = designvariables%nphi 
        
        hessG%nrow = designvariables%nphi 
        hessG%ncol = designvariables%nphi 

        ! Check design variable type
        select case(designvariables%type)

        case ('vesselcoordinates', 'vesselcoordinates_goat') ! same since vessel coordinates come first in design variable vector

            ! Gradient
            !---------
            if (dogradient) then
                ! Evaluate derivatives
                allocate(dpsidx(ntv), dpsidy(ntv))
                call mfinterp%Evaluate(xv, yv, 1, 0, dpsidx)
                call mfinterp%Evaluate(xv, yv, 0, 1, dpsidy) 
                 
                ! All contributions
                row = [[(k, k = 1, constraints%ncon)], [(k, k = 1, constraints%ncon)]]
                col = [tv, tv+nv]
                val = [dpsidx, dpsidy]
                jacG = ConstructMySparse(row, col, val, constraints%ncon, &
                    designvariables%nphi)
                                
                ! Transpose
                gradG = jacG%Transpose()
            end if 

            ! Hessian
            !--------
            if (dohessian) then
                
                ! Evaluate derivatives
                allocate(d2psidx2(ntv), d2psidy2(ntv), d2psidxdy(ntv))
                call mfinterp%Evaluate(xv, yv, 2, 0, d2psidx2)
                call mfinterp%Evaluate(xv, yv, 1, 1, d2psidxdy)
                call mfinterp%Evaluate(xv, yv, 0, 2, d2psidy2)

                ! All contributions
                row = [tv, tv, tv+nv, tv+nv]
                col = [tv, tv+nv, tv, tv+nv]
                val = [d2psidx2*lambda, d2psidxdy*lambda, &
                    d2psidxdy*lambda, d2psidy2*lambda]

                ! Simply zero
                hessG = ConstructMySparse(row, col, val, &
                    designvariables%nphi, designvariables%nphi)
            end if 

        case default

            call gdErrorHandler('EvaluateFixedVesselPointsConstraints: ' // &
                'derivatives not implemented for design variable type: ' // &
                designvariables%type)

        end select
        
        ! Housekeeping
        !=============
        ! End associate
        end associate

    end subroutine 

    ! Data output
    subroutine WriteDataFixedVesselFluxConstraints(constraints, goat)

        ! Description
        !============
        ! Write vessel nodes in the following format:
        ! ID, x, y 
        ! Different files are written for special vertices, fixed
        ! vertices, flux surface vertices, and tangency points

        ! The usual
        implicit none

        ! Declare variables
        !==================
        ! Arguments
        class(FixedVesselFluxConstraintsUDT)            :: constraints
        type(OptimizationProblemGDUDT)                  :: goat

        ! Auxiliary
        real(R8), allocatable, dimension(:)             :: x, y
        character(:), allocatable                       :: filepath

        ! Initialize
        !===========
        ! Set the correct directories
        allocate(character(len('so_con_fvf_vertices')) :: filepath)
        filepath = 'so_con_fvf_vertices'

        ! Associate
        associate(&
            ps              => goat%environment%vessel%polygonset,     &
            ID              => constraints%ID  &
            )

        ! Get vessel vertices
        call ps%GetVertices(x, y)

        ! Write
        call WriteVertexData(ID, x(ID), y(ID), filepath)

        ! Housekeeping
        !=============
        end associate

    end subroutine

    !------------------------------------------------------------------!
    !                           VESSEL DISTANCE                        !
    !------------------------------------------------------------------!
    ! Initialize
    subroutine InitializeVesselDistanceConstraints(constraints, goat, &
        monitor, designvariables, options)

        ! Description
        !============
        ! Initialize the vessel Distance constraints. 

        ! Notes
        !======
        
        ! Initialize
        !===========
        implicit none
        
        ! Declare variables
        !==================
        ! Arguments 
        class(VesselDistanceConstraintsUDT) :: constraints 
        type(OptimizationProblemGDUDT)          :: goat 
        type(ConstraintsMonitorSOUDT)           :: monitor
        type(ConstraintOptionsSOUDT)            :: options 
        class(DesignVariablesSOUDT)             :: designvariables

        ! Auxiliary variables
        class(PLF2DOptionsUDT), allocatable     :: plfoptions
        logical, allocatable                    :: isconstrained(:), &
            isvertextreated(:)
        integer(I8), allocatable                :: labels(:, :)

        real(R8), allocatable                   :: xv(:), yv(:)

        ! Loop variables
        integer(I8)                 :: i, k

        ! Initialize
        !===========
        ! Number of constraints
        constraints%ncon = 0

        ! Check design variable type
        select case (designvariables%type)

        case ('vesselcoordinates', 'vesselcoordinates_goat')

            ! All good

        case default

            ! All bad
            call gdErrorHandler('InitializeFixedVesselFluxConstraints: ' // &
                'design variable type not implemented')

        end select

        ! Associate
        associate(&
            ps      => goat%environment%vessel%polygonset,          &
            cc      => monitor%eqvcc,                               &
            maxcc   => monitor%maxeqvcc,                            &
            opt     => options%vdoptions                            &
            )

        ! Get vessel coordinates and labels (should be same as design variables)
        call ps%GetLabels(labels)
        call ps%GetVertices(xv, yv)

        ! Initialize logical indicating if vertex is constrained
        allocate(isconstrained(size(labels, 1)))
        isconstrained = .false. 

        ! Initialize logical indicating if desired constrained vertex
        ! was present and constrained
        allocate(isvertextreated(size(opt%vertIDs)))
        isvertextreated = .false. 

        ! Determine constrained vertices
        !===============================
        ! Constrain per vessel structure (label 1 and 2)
        do i = 1, size(opt%structureIDs)
            ! Unpack ID
            associate(tID       => opt%structureIDs(i))

            ! Check vertices
            where ( (labels(:, 1) == tID) .or. (labels(:, 2) == tID) ) &
                isconstrained = .true. 

            ! Housekeeping
            end associate
        end do

        ! Constrain per vertex ID
        do i = 1, size(opt%vertIDs)
            ! Unpack ID
            associate(tID       => opt%vertIDs(i))

            ! Check vertices
            where( (labels(:, 3) == tID)) isconstrained = .true. 

            ! Check if found
            if (any(labels(:, 3) == tID)) then 
                isvertextreated(i) = .true. 
            end if

            ! Housekeeping
            end associate
        end do

        ! Check if we can constrain, set to false if not the case
        where (cc + 1> maxcc) isconstrained = .false. 

        ! Add to constraints
        !===================
        ! Vertices, but locally indexed (i.e. not vertex ID, but ID 
        ! according to polygon structure)
        constraints%ID = pack([(k, k = 1, size(labels, 1))], isconstrained)

        ! set number of constraints
        constraints%ncon = size(constraints%ID)
        constraints%nID = constraints%ncon

        ! Update
        cc(constraints%ID) = cc(constraints%ID) + 1

        ! Construct constraint polygon levelset
        !======================================
        ! Construct polygon set
        call constraints%ps%Construct(opt%xp, opt%yp)

        ! Set options for levelset
        select case (opt%plftype)

        case ('closedpolygon_exact')

            ! Allocate
            allocate(PLF2DClosedExactOptionsUDT::plfoptions)

        case ('closedpolygon_smoothapproximation')

            ! Allocate
            allocate(PLF2DClosedApproximationOptionsUDT::plfoptions)

        case default 

            ! Throw error
            call gdErrorHandler('InitializeVesselDistanceConstraints: ' // & 
                'polygon levelset function: ' // opt%plftype // & 
                ' not implemented')

        end select

        select type (plfoptions)

        type is (PLF2DClosedExactOptionsUDT)

            ! Set options (none to set right now)

        type is (PLF2DClosedApproximationOptionsUDT)
            ! Set options
            plfoptions%resx     = opt%resx
            plfoptions%resy     = opt%resy
            plfoptions%C        = opt%C
            plfoptions%M        = opt%M
            plfoptions%offsetx  = opt%offsetx
            plfoptions%offsety  = opt%offsety
            plfoptions%meth     = opt%meth

        class default

            ! Throw error - should actually already be thrown above
            call gdErrorHandler('InitializeVesselDistanceConstraints: ' // & 
            'polygon levelset function: ' // opt%plftype // & 
                ' not implemented')

        end select

        ! Construct levelset
        call InitializePolygonLevelsetFunction2D(constraints%plf, &
            constraints%ps, plfoptions)

        ! Output
        !=======
        ! Display warning message if necessary
        if (any(.not. isvertextreated)) then 
            ! Display message
            print *, 'Vertices with following IDs were not present and ' // &
                'are not constrained: ', pack(opt%vertIDs, .not. isvertextreated)
        end if

        ! Print data
        call constraints%WriteData(goat)
        
        ! Housekeeping
        !=============
        end associate
        
    end subroutine

    ! Evaluation
    subroutine EvaluateVesselDistanceConstraints(constraints, G, gradG, & 
        hessG, goat, dogradient, dohessian, designvariables, lambda)

        ! Description
        !============
        ! Evaluate the vessel Distance constraints
        ! 
        !       plf(xv, yv) - d = 0
        ! 
        ! Gradients and hessians w.r.t. xv, yv are simply the 
        ! derivatives of the plf function. It is assumed that the plf
        ! gives some measure of the distance (only important if d is 
        ! non-zero)

        ! Initialize
        !===========
        ! Modules
        
        ! Declare variables
        !==================
        ! Arguments 
        class(VesselDistanceConstraintsUDT)     :: constraints 
        real(R8), allocatable                   :: G(:) 
        real(R8), allocatable                   :: lambda(:)
        type(MySparseUDT)                       :: hessG, gradG
        type(OptimizationProblemGDUDT)          :: goat 
        logical                                 :: dogradient, dohessian
        class(DesignVariablesSOUDT)             :: designvariables   
        
        ! Auxiliary
        integer(I8)                             :: nv, ntv
        integer(I8), allocatable, dimension(:)  :: row, col
        real(R8), allocatable, dimension(:)     :: xv, yv, val, V, &
            dVdx, dVdy, d2Vdx2, d2Vdy2, d2Vdxdy
        type(MySparseUDT)                       :: jacG

        ! Loop 
        integer(I8)                         :: ic, ivg, ivh, k

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
        if (allocated(G)) then 
            if (size(G, 1) .ne. constraints%ncon) then 
                deallocate(G)
                allocate(G(constraints%ncon))
            end if 
        else
            allocate(G(constraints%ncon))
        end if 

        ! Counters
        ic = 0 ! constraint counter (local)
        ivg = 0 ! value index for gradient
        ivh = 0 ! value index for hessian

        ! Associate
        associate(&
            plf         => constraints%plf,                 &
            tv          => constraints%ID,                  &
            d           => constraints%d,                   &
            ps          => goat%environment%vessel%polygonset   &
            )

        ! Get coordinates
        call ps%GetVertices(xv, yv)

        ! Get number of vessel vertices
        nv = size(xv)

        ! Take only coordinates that are considered
        xv = xv(tv)
        yv = yv(tv)

        ! Get number of constrained vertices
        ntv = constraints%nID 

        ! Evaluate
        !=========
        ! Values at vessel coordinates
        allocate(V(ntv))
        call plf%Evaluate(xv, yv, 0, 0, V)

        ! Constraint
        G(ic+1:ic+ntv) = V - d
        
        ! Derivatives
        !============
        ! Initialize
        jacG%nrow = constraints%ncon 
        jacG%ncol = designvariables%nphi 
        
        hessG%nrow = designvariables%nphi 
        hessG%ncol = designvariables%nphi 

        ! Check design variable type
        select case(designvariables%type)

        case ('vesselcoordinates', 'vesselcoordinates_goat') ! same since vessel coordinates come first in design variable vector

            ! Gradient
            !---------
            if (dogradient) then
                ! Evaluate derivatives
                allocate(dVdx(ntv), dVdy(ntv))
                call plf%Evaluate(xv, yv, 1, 0, dVdx)
                call plf%Evaluate(xv, yv, 0, 1, dVdy) 
                 
                ! All contributions
                row = [[(k, k = 1, constraints%ncon)], [(k, k = 1, constraints%ncon)]]
                col = [tv, tv+nv]
                val = [dVdx, dVdy]
                jacG = ConstructMySparse(row, col, val, constraints%ncon, &
                    designvariables%nphi)
                                
                ! Transpose
                gradG = jacG%Transpose()
            end if 

            ! Hessian
            !--------
            if (dohessian) then
                
                ! Evaluate derivatives
                allocate(d2Vdx2(ntv), d2Vdy2(ntv), d2Vdxdy(ntv))
                call plf%Evaluate(xv, yv, 2, 0, d2Vdx2)
                call plf%Evaluate(xv, yv, 1, 1, d2Vdxdy)
                call plf%Evaluate(xv, yv, 0, 2, d2Vdy2)

                ! All contributions
                row = [tv, tv, tv+nv, tv+nv]
                col = [tv, tv+nv, tv, tv+nv]
                val = [d2Vdx2*lambda, d2Vdxdy*lambda, &
                    d2Vdxdy*lambda, d2Vdy2*lambda]

                ! Simply zero
                hessG = ConstructMySparse(row, col, val, &
                    designvariables%nphi, designvariables%nphi)
            end if 

        case default

            call gdErrorHandler('EvaluateFixedVesselPointsConstraints: ' // &
                'derivatives not implemented for design variable type: ' // &
                designvariables%type)

        end select
        
        ! Housekeeping
        !=============
        ! End associate
        end associate

    end subroutine 

    ! Data output
    subroutine WriteDataVesselDistanceConstraints(constraints, goat)

        ! Description
        !============
        ! Write vessel nodes in the following format:
        ! ID, x, y 
        ! Different files are written for special vertices, fixed
        ! vertices, flux surface vertices, and tangency points

        ! The usual
        implicit none

        ! Declare variables
        !==================
        ! Arguments
        class(VesselDistanceConstraintsUDT)             :: constraints
        type(OptimizationProblemGDUDT)                  :: goat

        ! Auxiliary
        real(R8), allocatable, dimension(:)             :: x, y
        character(:), allocatable                       :: filepath

        ! Initialize
        !===========
        ! Set the correct directories
        allocate(character(len('so_con_vd_vertices')) :: filepath)
        filepath = 'so_con_vd_vertices'

        ! Associate
        associate(&
            ps              => goat%environment%vessel%polygonset,     &
            ID              => constraints%ID  &
            )

        ! Get vessel vertices
        call ps%GetVertices(x, y)

        ! Write
        call WriteVertexData(ID, x(ID), y(ID), filepath)

        ! Housekeeping
        !=============
        end associate

    end subroutine

    !------------------------------------------------------------------!
    !                     VESSEL FACE ANGLE DIFFERENCE                 !
    !------------------------------------------------------------------!
    ! Initialize
    subroutine InitializeVesselAngleDifferenceConstraints(constraints, goat, &
        monitor, designvariables, options)

        ! Description
        !============
        ! Initialize the vessel angle difference constraints. 

        ! Notes
        !======
        
        ! Initialize
        !===========
        implicit none
        
        ! Declare variables
        !==================
        ! Arguments 
        class(VesselAngleDifferenceConstraintsUDT)  :: constraints 
        type(OptimizationProblemGDUDT)          :: goat 
        type(ConstraintsMonitorSOUDT)           :: monitor
        type(ConstraintOptionsSOUDT)            :: options 
        class(DesignVariablesSOUDT)             :: designvariables

        ! Auxiliary variables
        logical, allocatable, dimension(:)      :: isconstrained, &
            isvertextreated, excludepairs

        integer(I8)                             :: tID, sID1, sID2, &
            nvpairs, tloc
        integer(I8), allocatable                :: labels(:, :), &
            vpairs(:, :)

        real(R8)                                :: thisangle
        real(R8), allocatable, dimension(:)     :: xv, yv, dthetas, &
            x1, x2, x3, y1, y2, y3, dx1, dx2, dy1, dy2, cp, dp, theta0

        ! Loop variables
        integer(I8)                 :: i

        ! Initialize
        !===========
        ! Associate
        associate(&
            ps      => goat%environment%vessel%polygonset,          &
            cc      => monitor%eqvcc,                               &
            maxcc   => monitor%maxeqvcc,                            &
            opt     => options%vadoptions                           &
            )

        ! Number of constraints
        constraints%ncon = 0

        ! Check design variable type
        select case (designvariables%type)

        case ('vesselcoordinates', 'vesselcoordinates_goat')

            ! All good

        case default

            ! All bad
            call gdErrorHandler('InitializeVesselAngleDifferenceConstraints: ' // &
                'design variable type not implemented')

        end select

        ! Check options
        if (size(opt%structureIDs) /= size(opt%dthetas)) then 
            call gdErrorHandler('InitializeVesselAngleDifferenceConstraints: ' // & 
                'incompatible sizes of structureIDs and dthetas in options')
        end if 
        if (size(opt%vertIDs) /= size(opt%dthetav)) then 
            call gdErrorHandler('InitializeVesselAngleDifferenceConstraints: ' // & 
                'incompatible sizes of vertIDs and dthetav in options')
        end if 

        ! Get vessel coordinates and labels (should be same as design variables)
        call ps%GetLabels(labels)
        call ps%GetVertices(xv, yv)

        ! Initialize logical indicating if vertex is constrained
        allocate(isconstrained(size(labels, 1)))
        isconstrained = .false. 

        ! Initialize logical indicating if desired constrained vertex
        ! was present and constrained
        allocate(isvertextreated(size(opt%vertIDs)))
        isvertextreated = .false. 

        ! Determine constrained vertices
        !===============================
        ! Also determine desired angles
        ! Constrain per vessel structure (label 1 and 2)
        do i = 1, size(opt%structureIDs)
            ! Unpack ID
            associate(tID       => opt%structureIDs(i))

            ! Check vertices
            where ( (labels(:, 1) == tID) .or. (labels(:, 2) == tID) ) &
                isconstrained = .true. 

            ! Housekeeping
            end associate
        end do

        ! Constrain per vertex ID
        do i = 1, size(opt%vertIDs)
            ! Unpack ID
            associate(tID       => opt%vertIDs(i))

            ! Check vertices
            where( (labels(:, 3) == tID)) isconstrained = .true. 

            ! Check if found
            if (any(labels(:, 3) == tID)) then 
                isvertextreated(i) = .true. 
            end if

            ! Housekeeping
            end associate
        end do

        ! Check if we can constrain, set to false if not the case
        where (cc + 1> maxcc ) isconstrained = .false. 

        ! Determine initial vertex pairs
        call goat%environment%vessel%GetVesselVertexPairs(vpairs, &
            opt%structureIDs, opt%vertIDs)
        nvpairs = size(vpairs, 1)

        ! Determine angles
        !=================
        ! Initialize
        allocate(constraints%dtheta_ub(nvpairs), constraints%dtheta_lb(nvpairs))
        constraints%dtheta_ub = posinfval_R8()
        constraints%dtheta_lb = -posinfval_R8()

        ! Exclude pairs of non-constrained vertices
        allocate(excludepairs(nvpairs))
        excludepairs = .false. 
        where (.not. isconstrained(vpairs(:, 2))) excludepairs = .true. 

        ! Determine initial angle
        x1 = xv(vpairs(:, 1))
        x2 = xv(vpairs(:, 2))
        x3 = xv(vpairs(:, 3))

        y1 = yv(vpairs(:, 1))
        y2 = yv(vpairs(:, 2))
        y3 = yv(vpairs(:, 3))

        dx1 = x2 - x1
        dx2 = x3 - x2 
        dy1 = y2 - y1
        dy2 = y3 - y2 

        dp = dx1*dx2 + dy1*dy2 
        cp = dx1*dy2 - dx2*dy1 

        theta0 = atan2(cp, dp) ! in radians

        ! Exclude non-treated pairs
        theta0 = pack(theta0, .not. excludepairs)
        nvpairs = size(theta0)

        ! Add vertex pairs and initial angle
        constraints%theta0 = theta0 
        allocate(constraints%vpairs(nvpairs, 3))
        do i = 1, 3
            constraints%vpairs(:, i) = pack(vpairs(:, i), .not. excludepairs)
        end do
        constraints%nvpairs = nvpairs 
        constraints%ncon = 2*nvpairs ! upper and lower bound

        ! Determine angle difference for each structure
        allocate(dthetas(max(maxval(labels(:, 1)), maxval(labels(:, 2)))))
        dthetas = posinfval_R8()
        do i = 1, size(opt%structureIDs)
            dthetas(opt%structureIDs(i)) = opt%dthetas(i)
        end do 

        ! Determine maximal angle for each vertex
        do i = 1, size(theta0)
            
            ! Get structure IDs
            tID = vpairs(i, 2)
            sID1 = labels(tID, 1)
            sID2 = labels(tID, 2)

            ! Get minimal angle according to structures
            if (sID2 /= 0) then 
                thisangle = min(dthetas(sID1), dthetas(sID2))
            else
                thisangle = dthetas(sID1)
            end if 
            tloc = findloc(opt%vertIDs, tID, 1)
            if (tloc /= 0) then 
                thisangle = min(thisangle, opt%dthetav(tloc)) ! also in radians
            end if 

            ! Check
            if (theta0(i) + thisangle <= (Pi_R8 - opt%alpha)) then 
                constraints%dtheta_ub(i) = Pi_R8 - opt%alpha 
            else
                constraints%dtheta_ub(i) = thisangle 
            end if 
            if (theta0(i) - thisangle <= (-Pi_R8 - opt%alpha)) then 
                constraints%dtheta_lb(i) = -Pi_R8 + opt%alpha 
            else
                constraints%dtheta_lb(i) = thisangle
            end if 

        end do 
        
        ! Output
        !=======
        ! Display warning message if necessary
        if (any(.not. isvertextreated)) then 
            ! Display message
            print *, 'InitializeVesselAngleDifferenceConstraints: ' // & 
                'Vertices with following IDs were not present and ' // &
                'are not constrained: ', pack(opt%vertIDs, .not. isvertextreated)
        end if

        ! Print data
        call constraints%WriteData(xv, yv)
        
        ! Housekeeping
        !=============
        end associate
        
    end subroutine

    ! Evaluation
    subroutine EvaluateVesselAngleDifferenceConstraints(constraints, G, gradG, & 
        hessG, goat, dogradient, dohessian, designvariables, lambda)

        ! Description
        !============
        ! Evaluate the vessel angle difference constraints
        ! 
        !       theta <= theta0 + dtheta_upper
        !       -theta <= -theta0 + dtheta_lower 
        ! 

        ! Initialize
        !===========
        ! Modules
        
        ! Declare variables
        !==================
        ! Arguments 
        class(VesselAngleDifferenceConstraintsUDT)     :: constraints 
        real(R8), allocatable                   :: G(:) 
        real(R8), allocatable                   :: lambda(:)
        type(MySparseUDT)                       :: hessG, gradG
        type(OptimizationProblemGDUDT)          :: goat 
        logical                                 :: dogradient, dohessian
        class(DesignVariablesSOUDT)             :: designvariables   
        
        ! Auxiliary
        integer(I8)                             :: nv
        integer(I8), allocatable, dimension(:)  :: row, col
        real(R8), allocatable, dimension(:)     :: xv, yv, &
            x1, x2, x3, y1, y2, y3, dx1, dx2, dy1, dy2, cp, &
            dp, theta, rat, val, valxx, valxy, valyx, valyy,lambdaub, &
            lambdalb
        type(MySparseUDT)                       :: jacG

        ! Loop 
        integer(I8)                         :: ic, ivg, ivh, k

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
        if (allocated(G)) then 
            if (size(G, 1) .ne. constraints%ncon) then 
                deallocate(G)
                allocate(G(constraints%ncon))
            end if 
        else
            allocate(G(constraints%ncon))
        end if 

        ! Counters
        ic = 0 ! constraint counter (local)
        ivg = 0 ! value index for gradient
        ivh = 0 ! value index for hessian

        ! Associate
        associate(&
            vpairs      => constraints%vpairs,              &
            nvpairs     => constraints%nvpairs,             &
            dthetaub    => constraints%dtheta_ub,           &
            dthetalb    => constraints%dtheta_lb,           &
            theta0      => constraints%theta0,              &
            ps          => goat%environment%vessel%polygonset   &
            )

        ! Get coordinates
        call ps%GetVertices(xv, yv)
        nv = size(xv)

        ! Evaluate
        !=========
        ! Evaluate angle
        x1 = xv(vpairs(:, 1))
        x2 = xv(vpairs(:, 2))
        x3 = xv(vpairs(:, 3))

        y1 = yv(vpairs(:, 1))
        y2 = yv(vpairs(:, 2))
        y3 = yv(vpairs(:, 3))

        dx1 = x2 - x1
        dx2 = x3 - x2 
        dy1 = y2 - y1
        dy2 = y3 - y2 

        dp = dx1*dx2 + dy1*dy2 
        cp = dx1*dy2 - dx2*dy1 

        rat = cp/dp

        theta = atan2(cp, dp) ! in radians

        ! Upper bound
        G(ic+1:ic+nvpairs) = theta - theta0 - dthetaub 
        ic = ic + nvpairs

        ! Lower bound
        G(ic+1:ic+nvpairs) = -(theta - theta0) - dthetalb 
        ic = ic + nvpairs
        
        ! Derivatives
        !============
        ! Initialize
        jacG%nrow = constraints%ncon 
        jacG%ncol = designvariables%nphi 
        
        hessG%nrow = designvariables%nphi 
        hessG%ncol = designvariables%nphi 

        ! Check design variable type
        select case(designvariables%type)

        case ('vesselcoordinates', 'vesselcoordinates_goat') ! same since vessel coordinates come first in design variable vector

            ! Gradient
            !---------
            if (dogradient) then
                ! Reset counter
                ivg = 0

                ! Allocate
                !---------
                ! Only for upper bound - lower bound is the same, but
                ! then with minus sign
                allocate(row(nvpairs*6), col(nvpairs*6), val(nvpairs*6))

                ! Upper bound
                !------------
                ! x1
                row(ivg+1:ivg+nvpairs) = [(k, k = 1, nvpairs)]
                col(ivg+1:ivg+nvpairs) = vpairs(:, 1)
                val(ivg+1:ivg+nvpairs) = -(dy2/dp - (dx2*rat)/dp)/(rat**2 + 1) !x1
                ivg = ivg + nvpairs

                ! x2
                row(ivg+1:ivg+nvpairs) = [(k, k = 1, nvpairs)]
                col(ivg+1:ivg+nvpairs) = vpairs(:, 2) 
                val(ivg+1:ivg+nvpairs) = -((y1 - y3)/dp &
                    - (rat*(dx1 - dx2))/dp)/(rat**2 + 1) !x2
                ivg = ivg + nvpairs

                ! x3
                row(ivg+1:ivg+nvpairs) = [(k, k = 1, nvpairs)]
                col(ivg+1:ivg+nvpairs) = vpairs(:, 3) 
                val(ivg+1:ivg+nvpairs) = -(dy1/dp + (dx1*rat)/dp)/(rat**2 + 1) !x3
                ivg = ivg + nvpairs

                ! y1
                row(ivg+1:ivg+nvpairs) = [(k, k = 1, nvpairs)] 
                col(ivg+1:ivg+nvpairs) = vpairs(:, 1) + nv
                val(ivg+1:ivg+nvpairs) = (dx2/dp + (dy2*rat)/dp)/(rat**2 + 1) !y1
                ivg = ivg + nvpairs

                ! y2
                row(ivg+1:ivg+nvpairs) = [(k, k = 1, nvpairs)] 
                col(ivg+1:ivg+nvpairs) = vpairs(:, 2) + nv
                val(ivg+1:ivg+nvpairs) = ((x1 - x3)/dp &
                    + (rat*(dy1 - dy2))/dp)/(rat**2 + 1) !y2
                ivg = ivg + nvpairs

                ! y3
                row(ivg+1:ivg+nvpairs) = [(k, k = 1, nvpairs)] 
                col(ivg+1:ivg+nvpairs) = vpairs(:, 3) + nv
                val(ivg+1:ivg+nvpairs) = (dx1/dp - (dy1*rat)/dp)/(rat**2 + 1) !y3
                ivg = ivg + nvpairs

                ! Lower bound
                !------------
                ! Not done explicitly here

                ! Concatenate
                !------------
                ! Construct Jacobian
                jacG = ConstructMySparse([row, row + nvpairs], [col, col], &
                    [val, -val], constraints%ncon, designvariables%nphi)
                                
                ! Transpose
                gradG = jacG%Transpose()
            end if 

            ! Hessian
            !--------
            if (dohessian) then

                ! Initialize counter
                ivh = 0

                ! Allocate
                if (allocated(row)) then ! assume all allocated
                    deallocate(row, col)
                end if
                allocate(row(nvpairs*9), col(nvpairs*9), valxx(nvpairs*9), &
                    valxy(nvpairs*9), valyx(nvpairs*9), valyy(nvpairs*9))

                ! Upper bound
                !------------
                ! Note: we don't account yet for multiplication with
                ! lambda, this comes afterwards

                ! v1v1
                row(ivh+1:ivh+nvpairs) = vpairs(:, 1)
                col(ivh+1:ivh+nvpairs) = vpairs(:, 1)
                valxx(ivh+1:ivh+nvpairs) = ((2*dx2**2*rat)/dp**2 &
                    - (2*dx2*dy2)/dp**2)/(rat**2 + 1) + (((2*dx2*rat**2)/dp &
                    - (2*dy2*rat)/dp)*(dy2/dp - (dx2*rat)/dp))/(rat**2 + 1)**2 !x1x1
                valxy(ivh+1:ivh+nvpairs) = (dx2**2/dp**2 - dy2**2/dp**2 &
                    + (2*dx2*dy2*rat)/dp**2)/(rat**2 + 1) + (((2*dy2*rat**2)/dp &
                    + (2*dx2*rat)/dp)*(dy2/dp - (dx2*rat)/dp))/(rat**2 + 1)**2 !x1y1
                valyx(ivh+1:ivh+nvpairs) = valxy(ivh+1:ivh+nvpairs)
                valyy(ivh+1:ivh+nvpairs) = ((2*dy2**2*rat)/dp**2 &
                    + (2*dx2*dy2)/dp**2)/(rat**2 + 1) - (((2*dy2*rat**2)/dp &
                    + (2*dx2*rat)/dp)*(dx2/dp + (dy2*rat)/dp))/(rat**2 + 1)**2 !y1y1
                ivh = ivh + nvpairs

                ! v1v2
                row(ivh+1:ivh+nvpairs) = vpairs(:, 1)
                col(ivh+1:ivh+nvpairs) = vpairs(:, 2)
                valxx(ivh+1:ivh+nvpairs) = - (rat/dp + (dy2*(dx1 - dx2))/dp**2 &
                    + (dx2*(y1 - y3))/dp**2 - (2*dx2*rat*(dx1 - dx2))/dp**2)/(rat**2 + 1) &
                    - (((2*rat*(y1 - y3))/dp - (2*rat**2*(dx1 - dx2))/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 !x1x2
                valxy(ivh+1:ivh+nvpairs) = (1/dp - (dy2*(dy1 - dy2))/dp**2 &
                    + (dx2*(x1 - x3))/dp**2 + (2*dx2*rat*(dy1 - dy2))/dp**2)/(rat**2 + 1) &
                    + (((2*rat*(x1 - x3))/dp + (2*rat**2*(dy1 - dy2))/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 !x1y2
                valyx(ivh+1:ivh+nvpairs) = (((2*rat*(y1 - y3))/dp &
                    - (2*rat**2*(dx1 - dx2))/dp)*(dx2/dp + (dy2*rat)/dp))/(rat**2 + 1)**2 &
                    - (1/dp - (dx2*(dx1 - dx2))/dp**2 + (dy2*(y1 - y3))/dp**2 &
                    - (2*dy2*rat*(dx1 - dx2))/dp**2)/(rat**2 + 1) !y1x2
                valyy(ivh+1:ivh+nvpairs) = ((dx2*(dy1 - dy2))/dp**2 - rat/dp &
                    + (dy2*(x1 - x3))/dp**2 + (2*dy2*rat*(dy1 - dy2))/dp**2)/(rat**2 + 1) &
                    - (((2*rat*(x1 - x3))/dp + (2*rat**2*(dy1 - dy2))/dp)*(dx2/dp &
                    + (dy2*rat)/dp))/(rat**2 + 1)**2 !y1y2
                ivh = ivh + nvpairs

                ! v2v1
                row(ivh+1:ivh+nvpairs) = vpairs(:, 2)
                col(ivh+1:ivh+nvpairs) = vpairs(:, 1)
                valxx(ivh+1:ivh+nvpairs) = valxx(ivh-nvpairs+1:ivh)
                valxy(ivh+1:ivh+nvpairs) = valyx(ivh-nvpairs+1:ivh)
                valyx(ivh+1:ivh+nvpairs) = valxy(ivh-nvpairs+1:ivh)
                valyy(ivh+1:ivh+nvpairs) = valyy(ivh-nvpairs+1:ivh)
                ivh = ivh + nvpairs

                ! v1v3
                row(ivh+1:ivh+nvpairs) = vpairs(:, 1)
                col(ivh+1:ivh+nvpairs) = vpairs(:, 3)
                valxx(ivh+1:ivh+nvpairs) = (rat/dp + (dx1*dy2)/dp**2 &
                    - (dx2*dy1)/dp**2 - (2*dx1*dx2*rat)/dp**2)/(rat**2 + 1) &
                    - (((2*dx1*rat**2)/dp + (2*dy1*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 !x1x3
                valxy(ivh+1:ivh+nvpairs) = - (1/dp - (dx1*dx2)/dp**2 &
                    - (dy1*dy2)/dp**2 + (2*dx2*dy1*rat)/dp**2)/(rat**2 + 1) &
                    - (((2*dy1*rat**2)/dp - (2*dx1*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 !x1y3
                valyx(ivh+1:ivh+nvpairs) = (((2*dx1*rat**2)/dp &
                    + (2*dy1*rat)/dp)*(dx2/dp + (dy2*rat)/dp))/(rat**2 + 1)**2 &
                    - ((dx1*dx2)/dp**2 - 1/dp + (dy1*dy2)/dp**2 &
                    + (2*dx1*dy2*rat)/dp**2)/(rat**2 + 1) !y1x3
                valyy(ivh+1:ivh+nvpairs) = (rat/dp + (dx1*dy2)/dp**2 &
                    - (dx2*dy1)/dp**2 - (2*dy1*dy2*rat)/dp**2)/(rat**2 + 1) &
                    + (((2*dy1*rat**2)/dp - (2*dx1*rat)/dp)*(dx2/dp &
                    + (dy2*rat)/dp))/(rat**2 + 1)**2 !y1y3
                ivh = ivh + nvpairs

                ! v3v1
                row(ivh+1:ivh+nvpairs) = vpairs(:, 3)
                col(ivh+1:ivh+nvpairs) = vpairs(:, 1)
                valxx(ivh+1:ivh+nvpairs) = valxx(ivh-nvpairs+1:ivh)
                valxy(ivh+1:ivh+nvpairs) = valyx(ivh-nvpairs+1:ivh)
                valyx(ivh+1:ivh+nvpairs) = valxy(ivh-nvpairs+1:ivh)
                valyy(ivh+1:ivh+nvpairs) = valyy(ivh-nvpairs+1:ivh)
                ivh = ivh + nvpairs

                ! v2v2
                row(ivh+1:ivh+nvpairs) = vpairs(:, 2)
                col(ivh+1:ivh+nvpairs) = vpairs(:, 2)
                valxx(ivh+1:ivh+nvpairs) = ((2*rat)/dp &
                    - (2*(dx1 - dx2)*(y1 - y3))/dp**2 + (2*rat*(dx1 &
                    - dx2)**2)/dp**2)/(rat**2 + 1) - (((2*rat*(y1 - y3))/dp &
                    - (2*rat**2*(dx1 - dx2))/dp)*((y1 - y3)/dp &
                    - (rat*(dx1 - dx2))/dp))/(rat**2 + 1)**2 !x2x2
                valxy(ivh+1:ivh+nvpairs) = (((dx1 - dx2)*(x1 - x3))/dp**2 &
                    - ((dy1 - dy2)*(y1 - y3))/dp**2 &
                    + (2*rat*(dx1 - dx2)*(dy1 - dy2))/dp**2)/(rat**2 + 1) &
                    + (((2*rat*(x1 - x3))/dp + (2*rat**2*(dy1 - dy2))/dp)*((y1 - y3)/dp &
                    - (rat*(dx1 - dx2))/dp))/(rat**2 + 1)**2 !x2y2
                valyx(ivh+1:ivh+nvpairs) = valxy(ivh+1:ivh+nvpairs)
                valyy(ivh+1:ivh+nvpairs) = ((2*rat)/dp + (2*(dy1 - dy2)*(x1 - x3))/dp**2 &
                    + (2*rat*(dy1 - dy2)**2)/dp**2)/(rat**2 + 1) &
                    - (((2*rat*(x1 - x3))/dp + (2*rat**2*(dy1 - dy2))/dp)*((x1 - x3)/dp &
                    + (rat*(dy1 - dy2))/dp))/(rat**2 + 1)**2 !y2y2
                ivh = ivh + nvpairs

                ! v2v3
                row(ivh+1:ivh+nvpairs) = vpairs(:, 2)
                col(ivh+1:ivh+nvpairs) = vpairs(:, 3)
                valxx(ivh+1:ivh+nvpairs) = - (rat/dp + (dy1*(dx1 - dx2))/dp**2 &
                    - (dx1*(y1 - y3))/dp**2 + (2*dx1*rat*(dx1 - dx2))/dp**2)/(rat**2 + 1) &
                    - (((y1 - y3)/dp - (rat*(dx1 - dx2))/dp)*((2*dx1*rat**2)/dp &
                    + (2*dy1*rat)/dp))/(rat**2 + 1)**2 !x2x3
                valxy(ivh+1:ivh+nvpairs) = (1/dp + (dx1*(dx1 - dx2))/dp**2 &
                    + (dy1*(y1 - y3))/dp**2 - (2*dy1*rat*(dx1 - dx2))/dp**2)/(rat**2 + 1) &
                    - (((y1 - y3)/dp - (rat*(dx1 - dx2))/dp)*((2*dy1*rat**2)/dp &
                    - (2*dx1*rat)/dp))/(rat**2 + 1)**2 !x2y3
                valyx(ivh+1:ivh+nvpairs) = (((x1 - x3)/dp &
                    + (rat*(dy1 - dy2))/dp)*((2*dx1*rat**2)/dp &
                    + (2*dy1*rat)/dp))/(rat**2 + 1)**2 - (1/dp + (dy1*(dy1 - dy2))/dp**2 &
                    + (dx1*(x1 - x3))/dp**2 + (2*dx1*rat*(dy1 - dy2))/dp**2)/(rat**2 + 1) !y2x
                valyy(ivh+1:ivh+nvpairs) = (((x1 - x3)/dp &
                    + (rat*(dy1 - dy2))/dp)*((2*dy1*rat**2)/dp &
                    - (2*dx1*rat)/dp))/(rat**2 + 1)**2 - (rat/dp - (dx1*(dy1 - dy2))/dp**2 &
                    + (dy1*(x1 - x3))/dp**2 + (2*dy1*rat*(dy1 - dy2))/dp**2)/(rat**2 + 1) !y2y3
                ivh = ivh + nvpairs

                ! v3v2
                row(ivh+1:ivh+nvpairs) = vpairs(:, 3)
                col(ivh+1:ivh+nvpairs) = vpairs(:, 2)
                valxx(ivh+1:ivh+nvpairs) = valxx(ivh-nvpairs+1:ivh)
                valxy(ivh+1:ivh+nvpairs) = valyx(ivh-nvpairs+1:ivh)
                valyx(ivh+1:ivh+nvpairs) = valxy(ivh-nvpairs+1:ivh)
                valyy(ivh+1:ivh+nvpairs) = valyy(ivh-nvpairs+1:ivh)
                ivh = ivh + nvpairs

                ! v3v3
                row(ivh+1:ivh+nvpairs) = vpairs(:, 3)
                col(ivh+1:ivh+nvpairs) = vpairs(:, 3)
                valxx(ivh+1:ivh+nvpairs) = ((2*dx1**2*rat)/dp**2 &
                    + (2*dx1*dy1)/dp**2)/(rat**2 + 1) - (((2*dx1*rat**2)/dp &
                    + (2*dy1*rat)/dp)*(dy1/dp + (dx1*rat)/dp))/(rat**2 + 1)**2 !x3x3
                valxy(ivh+1:ivh+nvpairs) = (dy1**2/dp**2 - dx1**2/dp**2 &
                    + (2*dx1*dy1*rat)/dp**2)/(rat**2 + 1) - (((2*dy1*rat**2)/dp &
                    - (2*dx1*rat)/dp)*(dy1/dp + (dx1*rat)/dp))/(rat**2 + 1)**2 !x3y3
                valyx(ivh+1:ivh+nvpairs) = valxy(ivh+1:ivh+nvpairs)
                valyy(ivh+1:ivh+nvpairs) = ((2*dy1**2*rat)/dp**2 &
                    - (2*dx1*dy1)/dp**2)/(rat**2 + 1) + (((2*dy1*rat**2)/dp &
                    - (2*dx1*rat)/dp)*(dx1/dp - (dy1*rat)/dp))/(rat**2 + 1)**2 !y3y3
                ivh = ivh + nvpairs

                ! Lower bound
                !------------
                ! Not done explicitly here

                ! Concatenate and multiply
                !-------------------------
                ! Need to deal with extra dimension introduced by spread
                lambdaub = reshape(spread(lambda(1:nvpairs), 2, 36), [36*nvpairs])
                lambdalb = reshape(spread(lambda(nvpairs+1:2*nvpairs), 2, 36), [36*nvpairs])
                hessG = ConstructMySparse(&
                    [row, row, row+nv, row+nv, row, row, row+nv, row+nv], &
                    [col, col+nv, col, col+nv, col, col+nv, col, col+nv], &
                    [[valxx, valxy, valyx, valyy]*lambdaub, &
                    -[valxx, valxy, valyx, valyy]*lambdalb], &
                    designvariables%nphi, designvariables%nphi) 
                
            end if 

        case default

            call gdErrorHandler('EvaluateFixedVesselPointsConstraints: ' // &
                'derivatives not implemented for design variable type: ' // &
                designvariables%type)

        end select
        
        ! Housekeeping
        !=============
        ! End associate
        end associate

    end subroutine 

    ! Data output
    subroutine WriteDataVesselAngleDifferenceConstraints(constraints, xv, yv)

        ! Description
        !============
        ! This routine writes out the cost function data (vertex pairs)
        ! for post-processing

        ! Declare variables
        !==================
        ! Arguments
        class(VesselAngleDifferenceConstraintsUDT)  :: constraints 
        real(R8), intent(in)                :: xv(:), yv(:)

        ! Auxiliary
        integer(I8)                     :: ncol, nrow 

        integer(I8), allocatable        :: IDn(:, :) 
        real(R8), allocatable           :: xn(:, :), yn(:, :)
        character(:), allocatable       :: filename 

        ! Loop
        integer(I8)                     :: j 

        ! Initialize
        !===========
        ! Set filename
        filename = 'so_con_vad_vpairs'

        ! Allocate
        nrow = size(constraints%vpairs, 1)
        ncol = size(constraints%vpairs, 2)
        allocate(IDn(nrow, ncol), xn(nrow, ncol), yn(nrow, ncol))

        ! Unpack
        associate(vpairs      => constraints%vpairs)

        ! Loop
        xn = 0
        yn = 0
        do j = 1, nrow 
            IDn(j, :) = vpairs(j, :) 
            xn(j, :) = xv(vpairs(j, :))
            yn(j, :) = yv(vpairs(j, :))
        end do

        ! Call writer
        !============
        call WriteVertexPairData(IDn, xn, yn, filename)

        ! Housekeeping
        !=============
        end associate
        deallocate(IDn, xn, yn)


    end subroutine

    !------------------------------------------------------------------!
    !                               GOAT                               !
    !------------------------------------------------------------------!

    ! Initialize
    subroutine InitializeGoatConstraints(constraints, goat, &
        monitor, designvariables, options)

        ! Description
        !============
        ! Basically a dummy routine since all parameters of the goat
        ! optimization should have already been set beforehand. The only
        ! thing we need to compute is the total amount of constraints. 
        ! This amount is equal to the amount of design variables and 
        ! constraints in the goat problem (we take all inequality 
        ! constraints since we work with the ncp function)
        
        ! Initialize
        !===========
        implicit none
        
        ! Declare variables
        !==================
        ! Arguments 
        class(GoatConstraintsUDT)               :: constraints 
        type(OptimizationProblemGDUDT)          :: goat 
        type(ConstraintsMonitorSOUDT)           :: monitor
        type(ConstraintOptionsSOUDT)            :: options 
        class(DesignVariablesSOUDT)             :: designvariables

        ! Determine number of constraints
        !================================
        ! Initialize
        constraints%ncon = 0

        ! Design contribution
        constraints%ncon = constraints%ncon + goat%designvariables%nphi 

        ! Equality constraints contribution
        constraints%ncon = constraints%ncon + goat%constraints%eqcon%neqcon 

        ! Inequality constraints contribution
        constraints%ncon = constraints%ncon + goat%constraints%ineqcon%nineqcon 
        
    end subroutine

    ! Evaluation
    subroutine EvaluateGoatConstraints(constraints, G, gradG, & 
        hessG, goat, dogradient, dohessian, designvariables, lambda)

        ! Description
        !============
        ! Evaluate the goat constraints. This means that we evaluate
        ! the value but also the gradient w.r.t the design variables.
        ! The hessian is (very) hard to compute analytically, since this
        ! involves computing third order derivatives of the constraints 
        ! in the goat optimization problem. Therefore, if the hessian
        ! (vector) product needs to be evaluated, it is returned as a 
        ! zero matrix of size nphi-by-nphi 

        ! Note: the gradients and hessian returned by goat are only
        ! w.r.t. vertex coordinates! We need to compute other 
        ! contributions separately and concatenate the final result
        ! ourselves. 

        ! Initialize
        !===========
        ! Modules
        
        ! Declare variables
        !==================
        ! Arguments 
        class(GoatConstraintsUDT)               :: constraints 
        real(R8), allocatable                   :: G(:) 
        real(R8), allocatable                   :: lambda(:)
        type(MySparseUDT)                       :: hessG, gradG
        type(OptimizationProblemGDUDT)          :: goat 
        logical                                 :: dogradient, dohessian
        class(DesignVariablesSOUDT)             :: designvariables   
        
        ! Auxiliary
        logical                                 :: dogradientg, &
            dohessiang
        logical, allocatable, dimension(:)      :: A, I 

        real(R8)                                :: Jg, L
        real(R8), allocatable, dimension(:)     :: xv, yv, &
            gradJg, Gg, Hg, gradL, rhs, ncp

        type(MySparseUDT)                       :: jacG, &
            hessJg, gradGg, hessGg, gradHg, hessHg, gradncpphi, &
            gradncpmu, hessL, lhs
        type(NumNCPUDT)                         :: num
        type(OptimizationSolverKKTUDT)          :: kktsolver

        ! Loop 

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
        if (allocated(G)) then 
            if (size(G, 1) .ne. constraints%ncon) then 
                deallocate(G)
                allocate(G(constraints%ncon))
            end if 
        else
            allocate(G(constraints%ncon))
        end if 

        ! Set numerics
        num%ncpfun = 'max'
        num%alpha = 1

        ! Check derivatives etc
        dogradientg = .true. ! always to be evaluated
        dohessiang = .false.
        if (dogradient) then 
            ! Also evaluate hessian
            dohessiang = .true. 
        end if 

        ! Initialize & allocate
        associate(nphi         => goat%designvariables%nphi, &
            neqcon              => goat%constraints%eqcon%neqcon, &
            nineqcon            => goat%constraints%ineqcon%nineqcon)
        allocate(gradJg(nphi), Gg(neqcon), Hg(nineqcon), ncp(nineqcon), &
            A(nineqcon), I(nineqcon), gradL(nphi+neqcon+nineqcon), &
            rhs(nphi+neqcon+nineqcon))
        hessJg = SpZeros(goat%designvariables%nphi, goat%designvariables%nphi)
        hessGg = hessJg 
        hessHg = hessJG 
        gradGg = SpZeros(goat%designvariables%nphi, neqcon)
        gradHg = SpZeros(goat%designvariables%nphi, nineqcon)
        

        ! Evaluate
        !=========
        ! Cost function
        call goat%EvaluateCostFunction(Jg, gradJg, hessJg, dogradientg, &
            dohessiang)

        ! Equality constraints
        call goat%EvaluateEqualityConstraints(Gg, gradGg, hessGg, dogradientg, &
            dohessiang, goat%lambda)

        ! Inequality constraints
        call goat%EvaluateInequalityConstraints(Hg, gradHg, hessHg, dogradientg, &
            dohessiang, goat%mu)

        ! Nonlinear complementarity function
        call EvaluateNCPfunction(ncp, A, I, gradncpphi, gradncpmu, &
            Hg, gradHg, goat%mu, num, dogradientg)

        ! Evaluate lagrangian
        call goat%EvaluateLagrangian(L, gradL, hessL, Jg, gradJg, &
            hessJg, Gg, gradGg, hessGg, goat%lambda, Hg, gradHg, hessHg, goat%mu, A, &
            dogradient, dohessian)

        ! Build goat hessian & rhs
        call kktsolver%SetupCorrectionEquation(lhs, rhs, &
            gradJg, hessJg, Gg, gradGg, hessGg, goat%lambda, Hg, gradHg, hessHg, goat%mu, A, &
            ncp, gradncpphi, gradncpmu)

        ! Relax lhs? Only if hessian is evaluated! Will yield wrong 
        ! constraint gradients, but that shouldn't be an issue 
        !if (dohessiang) then 
        !    kktsolver%numKKT%rxf = 1-5
        !    call kktsolver%RelaxKKTSystem(lhs, nphi, neqcon, nineqcon)
        !end if 
        
        ! Set constraint value
        G = -rhs ! need to compensate for minus sign in SetUpCorrectionEquation

        if (any(ieee_is_nan(rhs))) then 
            print *, 'NaN in goat constraints'
        end if 
        if (any(.not. ieee_is_finite(rhs))) then 
            print *, 'Inf in goat constraints'
        end if 
            
        ! Construct linearization
        !========================
        if (dogradient) then

            select case (designvariables%type)

            case ('vesselcoordinates')

                ! Only derivatives w.r.t. vessel coordinates - probably not
                ! used? 
                call goat%environment%vessel%polygonset%GetVertices(xv, yv)
                call goat%EvaluateJacobian('vesselcoordinates', [xv, yv], jacG)

                ! Transpose
                gradG = jacG%Transpose()

            case ('vesselcoordinates_goat')

                ! Derivatives w.r.t. vessel coordinates first, then goat
                ! linearization
                call goat%environment%vessel%polygonset%GetVertices(xv, yv)
                call goat%EvaluateJacobian('vesselcoordinates', [xv, yv], jacG)

                ! Concatenate with goat hessian
                jacG = jacG%Concatenate(lhs, 2)

                ! Transpose
                gradG = jacG%Transpose()

            case default 

                ! Throw error
                call gdErrorHandler('EvaluateGoatConstraints: design variable' // &
                    ' not implemented: ' // designvariables%type)

            end select

        else 

            ! Simply zero
            gradG = SpZeros(designvariables%nphi, constraints%ncon)

        end if 

        ! Check
        if (any(ieee_is_nan(gradG%val))) then 
            print *, 'NaN in goat constraints linearization'
        end if 
        if (any(.not. ieee_is_finite(gradG%val))) then 
            print *, 'Inf in goat constraints linearization'
        end if 

        ! Hessian
        !========
        ! Simply set to zero
        hessG = SpZeros(designvariables%nphi, designvariables%nphi)

        ! Housekeeping
        !=============
        end associate
        
    end subroutine


    
        

end module
