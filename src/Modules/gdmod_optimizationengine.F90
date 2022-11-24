!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains the implementation of the optimization engine for
! the grid deformation. It inherits from the base generic optimization
! structures in the optimization modules. 

! Notes
!======
! Note 1: Descriptions of what the deferred procedures should do are 
! provided in the abstract interface. 

module gdmod_optimizationengine
    
    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use optmod_optimizationengine
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
    type, extends(OptimizationProblemUDT) :: OptimizationProblemGDUDT

        ! Description
        !============
        ! Overwrite the initial design routines with our own
        ! implemented and derived types. 

        ! The classical derived types
        class(DesignVariablesGDUDT), allocatable :: designvariables
        class(CostfunctionGDUDT), allocatable    :: costfunction
        type(ConstraintsGDUDT)                   :: constraints

        ! Additional fields that are needed
        type(GridUDT)                       :: grid 
        type(MagneticFieldUDT)              :: magneticField
        type(EnvironmentUDT)                :: environment
        
    contains

        ! Overwrite 
        !==========
        ! Problem initialization
        procedure :: Initialize      => InitializeOptimizationProblemGD

        ! Dimension query
        procedure :: GetProblemDimensions   => GetProblemDimensionsGD 

        ! Design variables query
        procedure :: GetProblemDesignVariables => &
            GetProblemDesignVariablesGD

        ! Design update
        procedure :: UpdateDesign           => UpdateDesignGD

        ! Problem update
        procedure :: UpdateProblem          => UpdateProblemGD

        ! Cost function evaluation
        procedure :: EvaluateCostFunction   => EvaluateCostFunctionGD

        ! Equality constraints evaluation
        procedure :: EvaluateEqualityConstraints    &
                        => EvaluateEqualityConstraintsGD

        ! Inequality constraints evaluation
        procedure :: EvaluateInequalityConstraints &
                        => EvaluateInequalityConstraintsGD

        ! Additional routines
        !====================

    end type 

    ! Optimization engine
    !====================
    type, extends(OptimizationEngineUDT) :: OptimizationEngineGDUDT 


    contains 

        ! Procedure for initialization of engine
        procedure :: SetupOptimizationDriver => SetupOptimizationDriverGD

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
    subroutine SetupOptimizationDriverGD(optimizationdriver) 

        ! Description
        !============
        ! This is only a dummy routine for setting up the optimization 
        ! engine and should be overwritten. Here, the design options etc
        ! are called to initialize the design. 

        ! Initialize
        !===========
        ! Declare modules
        !use gdmod_types
        !use gdmod_userinput

        ! The usual
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationEngineGDUDT)      :: optimizationdriver 

        ! Loop variables
        type(DesignVariablesCoordinatesUDT) :: thisdesign
        type(OptimizationProblemGDUDT)      :: thisproblem

        ! Auxiliary variables 
        type(DesignOptionsUDT)              :: designoptions
        real(R8) :: dummy

        ! Data

        ! Design variables
        !=================
        ! Set the design options
        call SetDesignOptions(designoptions)
        allocate(optimizationdriver%problem, source=thisproblem)
        !allocate(optimizationdriver%problem%designvariables, source=thisdesign)
        !allocate(optimizationdriver%state, source = thisstate)
        !call optimizationdriver%problem%designvariables%Initialize(dummy)
        !call thisdesign%Initialize(dummy)

    end subroutine

    !------------------------------------------------------------------!
    !                       OPTIMIZATION PROBLEM                       !
    !------------------------------------------------------------------!
    ! Dimension query
    subroutine GetProblemDimensionsGD(problem, nphi, neq, nineq)

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
        class(OptimizationProblemGDUDT)     :: problem 
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
    subroutine GetProblemDesignVariablesGD(problem, phi)

        ! Description
        !============
        ! This routine returns the current design variable vector phi

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemGDUDT)         :: problem
        real(R8), allocatable, intent(out)      :: phi(:)

        ! Get!
        !=====
        phi = problem%designvariables%phi

    end subroutine

    ! Optimization problem initialization
    subroutine InitializeOptimizationProblemGD(problem) 

        ! Description
        !============
        ! This routine further initializes the design variables, cost
        ! function, and constraints. It is assumed that the grid, 
        ! magnetic field, and environment structures are properly 
        ! assigned. The allocatable components of the problem, namely 
        ! the cost function and design variables.

        ! Initialize
        !===========
        ! Declare modules
        use gdmod_types
        use gdmod_userinput

        ! The usual
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemGDUDT)      :: problem

        ! Loop variables

        ! Auxiliary variables 
        type(DesignOptionsUDT)              :: designoptions

        ! Data

        ! Design options
        !===============
        ! Set the design options
        call SetDesignOptions(designoptions)

        ! Design variables
        !=================
        ! Allocate the design variables, depending on the type.
        select case (trim(designoptions%variables%type))

        case ('coordinates')

            ! Only coordinates
            
            allocate(DesignVariablesCoordinatesUDT::problem%designvariables)

        case default

            ! Throw error
            call gdErrorHandler('Unknown design variable type')

        end select

        ! Initialize the design variables
        call problem%designvariables%Initialize(problem%grid, &
            problem%magneticField, problem%environment)
        
        ! Cost function
        !==============
        ! Allocate the cost function, depending on the type
        select case (trim(designoptions%costfunction%type))

        case ('lengthratio')

            ! Allocate
            allocate(CostFunctionLRUDT::problem%costfunction)

        case ('lengthratio2')

            ! Allocate
            allocate(CostFunctionLRUDT2::problem%costfunction)

        case ('FAD')

            ! Allocate
            allocate(CostfunctionFADUDT::problem%costfunction)

        case default
            
            call gdErrorHandler('Unknown cost function type')

        end select

        ! Initialize the cost function
        call problem%costfunction%Initialize(problem%grid, &
            problem%magneticField, problem%environment)

        ! Constraints
        !============
        ! Given the (many) possible options for the constraints, the 
        ! constraints are set in its own initialization. 
        call problem%constraints%Initialize(problem%grid, &
            problem%magneticField, problem%environment, & 
            designoptions%constraints)


    end subroutine

    ! Problem update
    subroutine UpdateProblemGD(problem)

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
        class(OptimizationProblemGDUDT)     :: problem 

        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Update design
        !==============
        ! Simply call the update routine from the design variables
        call problem%designvariables%UpdateDesign(problem%grid, &
            problem%magneticField, problem%environment)

    end subroutine

    ! Design update
    subroutine UpdateDesignGD(problem, dx)

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
        class(OptimizationProblemGDUDT)     :: problem 
        real(R8), intent(in)                :: dx(:)
    
        ! Loop variables

        ! Auxiliary variables 
        logical                             :: dodebugplots

        ! Data

        ! Do NOT set this to true unless you
        ! know what you're doing! prints out a figure, which may be in 
        ! each optimization iteration!
        
        data dodebugplots /.false./ 

        ! Update design
        !==============
        ! Simply update designvariables%phi
        problem%designvariables%phi = problem%designvariables%phi + dx 

        ! Debug plots?
        !=============
        if (dodebugplots) then
            ! Call grid visualization
            call PlotGridCells(problem%grid, '-p')
        end if

    end subroutine

    ! Cost function evaluation
    subroutine EvaluateCostFunctionGD(problem, J, gradJ, hessJ, &
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
        class(OptimizationProblemGDUDT)     :: problem 
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
        call problem%costfunction%Evaluate(J, gradJ, hessJ, problem%grid, &
            problem%magneticField, problem%environment, dogradient, &
            dohessian, problem%designvariables)

    end subroutine

    ! Equality constraints evaluation
    subroutine EvaluateEqualityConstraintsGD(problem, G, gradG, hessG, &
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
        class(OptimizationProblemGDUDT)     :: problem 
        real(R8), allocatable               :: G(:), lambda(:)
        type(MySparseUDT)                   :: gradG, hessG
        logical                             :: dogradient, dohessian                            

        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Compute the constraints
        !========================
        ! Simply call the constraint computation routine
        call problem%constraints%eqcon%Evaluate(G, gradG, hessG, problem%grid, &
            problem%magneticField, problem%environment, dogradient, &
            dohessian, problem%designvariables, lambda)

    end subroutine

    ! Inequality constraints evaluation
    subroutine EvaluateInequalityConstraintsGD(problem, H, gradH, &
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
        class(OptimizationProblemGDUDT)     :: problem 
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
            problem%grid, problem%magneticField, problem%environment, &
            dogradient, dohessian, problem%designvariables, mu)

    end subroutine


end module