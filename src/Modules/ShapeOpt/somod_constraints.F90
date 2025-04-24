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

        type(VesselDistanceConstraintsUDT)  :: vesselupperbound 
        type(VesselDistanceConstraintsUDT)  :: vessellowerbound

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
        where (cc > maxcc + 2) isconstrained = .false. 

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
        where (cc > maxcc + 1) isconstrained = .false. 

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
        where (cc > maxcc + 1) isconstrained = .false. 

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
