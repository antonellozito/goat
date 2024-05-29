!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains the implementation of the optimization engine for
! shape optimization (without SOLPS for now). 

! Notes
!======
! Note 1: Descriptions of what the deferred procedures should do are 
! provided in the abstract interface. 

module somod_optimizationengine
    
    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use optmod_optimizationengine
    use somod_designvariables
    use somod_costfunction
    use somod_constraints
    use gdmod_optimizationengine
    use gdmod_designvariables
    use gdmod_costfunction
    use gdmod_constraints
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

    ! Optimization problem
    !=====================
    type, extends(OptimizationProblemUDT) :: OptimizationProblemSOUDT

        ! Description
        !============
        ! Overwrite the initial design routines with our own
        ! implemented and derived types. 

        ! The classical derived types
        class(DesignVariablesSOUDT), allocatable :: designvariables
        class(CostfunctionSOUDT), allocatable    :: costfunction
        type(ConstraintsSOUDT)                   :: constraints

        ! Additional fields that are needed
        type(OptimizationProblemGDUDT)      :: goat 
        type(DesignOptionsSOUDT)            :: designoptions 
        type(CGStructureUDT), allocatable   :: congroups(:) 
        type(DOFGStructureUDT), allocatable :: dofgroups(:)
        
    contains

        ! Overwrite 
        !==========
        ! Problem initialization
        procedure :: Initialize      => InitializeOptimizationProblemSO

        ! Dimension query
        procedure :: GetProblemDimensions   => GetProblemDimensionsSO 

        ! Design variables query
        procedure :: GetProblemDesignVariables => &
            GetProblemDesignVariablesSO

        ! Design update
        procedure :: UpdateDesign           => UpdateDesignSO

        ! Problem update
        procedure :: UpdateProblem          => UpdateProblemSO

        ! Cost function evaluation
        procedure :: EvaluateCostFunction   => EvaluateCostFunctionSO

        ! Equality constraints evaluation
        procedure :: EvaluateEqualityConstraints    &
                        => EvaluateEqualityConstraintsSO

        ! Inequality constraints evaluation
        procedure :: EvaluateInequalityConstraints &
                        => EvaluateInequalityConstraintsSO

        ! Additional routines
        !====================
        ! Initialization finalizer to account for cross-design/cfv/con
        ! initialization requirements
        procedure :: FinalizeInitialization => FinalizeInitializationSO  

    end type 

    ! Optimization engine
    !====================
    type, extends(OptimizationEngineUDT) :: OptimizationEngineSOUDT 


    contains 

        ! Procedure for initialization of engine
        procedure :: SetupOptimizationDriver => SetupOptimizationDriverSO

    end type

    ! Abstract types
    !===============

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    contains

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                       OPTIMIZATION ENGINE                        !
    !------------------------------------------------------------------!
    ! Optimization engine initialization
    subroutine SetupOptimizationDriverSO(optimizationdriver) 

        ! Description
        !============
        ! This routine should be used to initialize the solver and
        ! problem to the desired type. Further initialization is done
        ! through the problem and solver initialization routines. 

        ! The usual
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationEngineSOUDT)      :: optimizationdriver 

        ! Auxiliary
        type(OptimizationProblemSOUDT)      :: thisproblem
        type(OptimizationSolverKKTUDT)      :: thissolver

        ! Data

        ! Initialize
        !===========
        ! Allocate the problem type
        allocate(optimizationdriver%problem, source=thisproblem)

        ! Allocate the solver
        allocate(optimizationdriver%solver, source=thissolver)

    end subroutine

    !------------------------------------------------------------------!
    !                       OPTIMIZATION PROBLEM                       !
    !------------------------------------------------------------------!
    ! Dimension query
    subroutine GetProblemDimensionsSO(problem, nphi, neq, nineq)

        ! Description
        !============
        ! Return the problem dimensions nphi, neq, nineq (number of 
        ! design variables, equality contraints, and inequality 
        ! constraints, resp. ). It is assumed that the problem is
        ! already properly initialized. 

        ! nphi is obtained from designvariables%nphi
        ! TO DO: change neq, nineq 

        ! Initialize
        !===========
        ! Declare modules

        ! The usual
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemSOUDT)     :: problem 
        integer(I8), intent(out)            :: nphi, neq, nineq

        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Output
        !=======
        ! Design variables
        nphi = problem%designvariables%nphi

        ! Constraints - call routine
        call problem%constraints%GetConstraintsDimensions(neq, nineq)

    end subroutine

    ! Design variables query
    subroutine GetProblemDesignVariablesSO(problem, phi)

        ! Description
        !============
        ! This routine returns the current design variable vector phi

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemSOUDT)         :: problem
        real(R8), allocatable, intent(out)      :: phi(:)

        ! Get!
        !=====
        phi = problem%designvariables%phi

    end subroutine

    ! Optimization problem initialization
    subroutine InitializeOptimizationProblemSO(problem) 

        ! Description
        !============
        ! This routine further initializes the design variables, cost
        ! function, and constraints. It is assumed that the auxiliary
        ! problem structures have been initialized beforehand. 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemSOUDT)                 :: problem

        ! Auxiliary
        type(OptimizationEngineGDUDT)      :: goatengine  
        type(ShapeOptimizationOptionsUDT)  :: SOoptions

        ! Data

        ! Filepaths
        !==========
        ! Shape optimization problem filepaths
        SOoptions%inputfilepath = problem%inputfilepath

        ! Design options
        problem%designoptions%inputfilepath = problem%inputfilepath

        ! Design options
        !===============
        ! Set shape optimization options
        call SOoptions%Set()

        ! Set the design options
        call problem%designoptions%Set()

        ! GOAT
        !=====
        ! Initialize the GOAT optimization driver
        call GDinitialize(SOoptions%goatfilepath, goatengine)

        ! Assign problem
        associate(goat      => goatengine%problem)
        select type(goat)

        type is (OptimizationProblemGDUDT)

            problem%goat = goat

        class default

            call gdErrorHandler('InitializeOptimizationProblemSO: ' // & 
                'goat optimization problem type not supported')

        end select
        end associate

        ! Design variables
        !=================
        ! Allocate the design variables, depending on the type.
        select case (trim(problem%designoptions%variables%type))

        case ('vesselcoordinates')

            ! Vessel coordinates
            allocate(DesignVariablesVesselCoordinatesUDT::problem%designvariables)

        case default

            ! Throw error
            call gdErrorHandler('Unknown design variable type')

        end select

        ! Initialize the design variables
        call problem%designvariables%Initialize(problem%goat)
        
        ! Cost function
        !==============
        ! Allocate the cost function, depending on the type
        select case (trim(problem%designoptions%costfunction%type))

        case ('PLF')

            ! Allocate 
            allocate(CostFunctionPLFUDT::problem%costfunction)

        case default
            
            ! Throw error
            call gdErrorHandler('Unknown cost function type')

        end select

        ! Initialize the cost function
        call problem%costfunction%Initialize(problem%goat, &
            problem%designoptions%costfunction)

        ! Constraints
        !============
        ! Given the (many) possible options for the constraints, the 
        ! constraints are set in its own initialization. 
        call problem%constraints%Initialize(problem%goat, & 
            problem%designvariables, problem%designoptions%constraints)

        ! 
        !=================
        ! Initialize design variables further for constraint/cfv 
        ! dependent fields
        call problem%FinalizeInitialization()


    end subroutine

    ! Finalize the problem initialization
    subroutine FinalizeInitializationSO(problem)

        ! Description
        !============
        ! This routine further initializes the optimization problem 
        ! after specific routines for initializing the design variables,
        ! constraints, and cost function have been called. As such, 
        ! interdependencies between these objects can be accounted for.
        ! For example, the amount of design variables for the 
        ! 'fluxvalue' (desired psi) type of design variables depends 
        ! on the number of flux function constraints. 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemSOUDT)     :: problem 

        ! Auxiliary

        ! Loop

        ! Initialize
        !===========
        ! Associate

        ! Update constraints
        !===================
        ! Construct constraint groups
        ! call problem%ConstructInequalityConstraintGroups()

        ! Housekeeping
        !=============

    end subroutine

    ! Problem update
    subroutine UpdateProblemSO(problem)

        ! Description
        !============
        ! Update the problem according to the design variables. Here,
        ! the grid, magnetic field, and environment should be updated,
        ! depending on the type of design variable. 

        ! Initialize
        !===========
        ! Declare modules

        ! The usual
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemSOUDT)     :: problem 

        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Initialize
        !===========
        ! Associate to use select type
        associate(&
            designvariables     => problem%designvariables, &
            constraints         => problem%constraints,     &
            goat                => problem%goat             &
            )

        ! Update design
        !==============
        ! Simply call the update routine from the design variables
        call designvariables%UpdateDesign(goat)

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Design update
    subroutine UpdateDesignSO(problem, dx)

        ! Description
        !============
        ! Update the design variables according to the update dx. 

        ! Initialize
        !===========
        ! Declare modules
        use gdmod_plots

        ! The usual
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemSOUDT)     :: problem 
        real(R8), intent(in)                :: dx(:)
    
        ! Loop variables

        ! Auxiliary variables 
        character(:), allocatable           :: vertpath, cellpath

        ! Data

        ! Update design
        !==============
        ! Simply update designvariables%phi
        problem%designvariables%phi = problem%designvariables%phi + dx 

        ! Debug plots?
        !=============
        ! Call grid vertex writing routine
        allocate(character(len('vertices_iterate')) :: vertpath)
        allocate(character(len('cells_iterate')) :: cellpath)
        vertpath = 'vertices_iterate'
        cellpath = 'cells_iterate'
        call WriteGridVertices(problem%goat%grid, vertpath) 
        call WriteGridCells(problem%goat%grid, cellpath)

    end subroutine

    ! Cost function evaluation
    subroutine EvaluateCostFunctionSO(problem, J, gradJ, hessJ, &
        dogradient, dohessian)

        ! Description
        !============
        ! Evaluate the cost function, based on the current state of the
        ! optimization problem. Besides the scalar value itself, also
        ! the gradient and hessian are evaluated, if needed. This can 
        ! be set using the dogradient and dohessian logicals. 

        ! Initialize
        !===========
        ! Declare modules

        ! The usual
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemSOUDT)     :: problem 
        real(R8)                            :: J 
        real(R8), allocatable               :: gradJ(:)
        type(MySparseUDT)                   :: hessJ
        logical                             :: dogradient, dohessian                            

        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Compute the cost function
        !==========================
        ! Simply call the cost function computation routine
        call problem%costfunction%Evaluate(J, gradJ, hessJ, problem%goat, &
            dogradient, dohessian, problem%designvariables)

    end subroutine

    ! Equality constraints evaluation
    subroutine EvaluateEqualityConstraintsSO(problem, G, gradG, hessG, &
        dogradient, dohessian, lambda)

        ! Description
        !============
        ! Evaluate the constraints, based on the current state of the
        ! problem. Note that the 'hessian' that is returned, is in fact
        ! the hessian-vector multiplication between the hessian and the
        ! vector lambda. Lambda should have been initialized with the 
        ! proper size. 

        ! Initialize
        !===========
        ! Declare modules

        ! The usual
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemSOUDT)     :: problem 
        real(R8), allocatable               :: G(:), lambda(:)
        type(MySparseUDT)                   :: gradG, hessG
        logical                             :: dogradient, dohessian                            

        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Compute the constraints
        !========================
        ! Simply call the constraint computation routine
        call problem%constraints%eqcon%Evaluate(G, gradG, hessG, &
            problem%goat, dogradient, &
            dohessian, problem%designvariables, lambda)

    end subroutine

    ! Inequality constraints evaluation
    subroutine EvaluateInequalityConstraintsSO(problem, H, gradH, &
        hessH, dogradient, dohessian, mu)

        ! Description
        !============
        ! Evaluate the constraints, based on the current state of the
        ! problem. Note that the 'hessian' that is returned, is in fact
        ! the hessian-vector multiplication between the hessian and the
        ! vector mu. mu should have been initialized with the 
        ! proper size. 

        ! Initialize
        !===========
        ! Declare modules

        ! The usual
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemSOUDT)     :: problem 
        real(R8), allocatable               :: H(:), mu(:)
        type(MySparseUDT)                   :: gradH, hessH
        logical                             :: dogradient, dohessian                            

        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Compute the constraints
        !========================
        ! Simply call the constraint computation routine
        call problem%constraints%ineqcon%Evaluate(H, gradH, hessH, &
            problem%goat, &
            dogradient, dohessian, problem%designvariables, mu)

        ! Deactivate constraints
        !=======================
        ! Only if necessary - determined in the routine itself
        ! Note: should we add switch here? And adjust gradient/Hessian
        ! as well? At least gradient should be tackled by activeness
        ! of constraint...
        ! call problem%DeactivateInequalityConstraints(H)

    end subroutine


    


end module