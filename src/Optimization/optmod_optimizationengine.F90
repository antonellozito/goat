!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains the basic optimization engine class. 
! Instantiation (i.e. creating an object) yields an object with the 
! 'problem' and 'solver' derived types, where the solver has a method
! 'solve' that can be applied to the 'problem' type to yield the 
! solution of the optimization problem at the end. To abstract the 
! optimization engine from any external data structures, routines etc
! needed to evaluate the cost function etc., the 'state' class is made
! available (see optmod_state.F90), that can contain the necessary 
! structures, data, ... for evaluation of any other subroutines. Though
! this results in some programming overhead, it ensures perfect 
! encapsulation of the optimization routines and modularity of the code. 

! Notes
!======
! Note 1: Descriptions of what the deferred procedures should do are 
! provided in the abstract interface. 

module optmod_optimizationengine
    
    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_constants
    use mod_utility
    use mod_sparseinterface
    use mod_diagnostics
    use mod_linearsolverinterface
    use mod_plotter
    use mod_errorhandler
    use optmod_designvariables
    use optmod_costfunction
    use optmod_constraints
    use optmod_numerics
    use optmod_monitor

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
    type, abstract :: OptimizationProblemUDT

        ! Description
        !============
        ! Defines the basic optimization problem. All that is actually
        ! required, are deferred procedures that update the design, 
        ! compute the cost function, constraints (and their derivatives)
        ! and a routine to update other structures of the problem, if 
        ! necessary. It is highly adviseable to use the provided 
        ! design variable, cost function, and constraint types in optmod, 
        ! as these form the basis of any decent optimization problem. 
        ! However, this is not a strict requirement to be able to use 
        ! this module. The only field that is added, is a monitor 
        ! structure to keep track of the progress of the algorithm. 
        ! (to be moved to the engine in the future)

        ! Note: additionally, we store the lagrange multipliers of
        ! equality and inequality constraints, and the active set
        character(:), allocatable                   :: inputfilepath
        type(OptimizationMonitorUDT)                :: monitor

        real(R8), allocatable                       :: lambda(:), mu(:)
        logical, allocatable                        :: A(:)

    contains 

        ! Dummy setup routine - to be replaced by the user
        ! General initialization routine 
        ! Design initialization

        ! Get problem dimensions
        procedure(GetProblemDimensionsINT), deferred :: & 
            GetProblemDimensions

        ! Get problem design variables
        procedure(GetProblemDesignVariablesINT), deferred :: & 
            GetProblemDesignVariables

        ! Design initialization
        procedure(InitializeINT), deferred :: &
            Initialize 

        ! Design updates
        procedure(UpdateDesignINT), deferred :: &
            UpdateDesign     

        ! Problem updates
        procedure(UpdateProblemINT), deferred :: &
            UpdateProblem
        
        ! Cost function evaluation
        procedure(EvaluateCostFunctionINT), deferred:: &
            EvaluateCostFunction

        ! Equality constraints evaluation
        procedure(EvaluateEqualityConstraintsINT), deferred :: &
            EvaluateEqualityConstraints

        ! Inequality constraints evaluation
        procedure(EvaluateInequalityConstraintsINT), deferred :: &
            EvaluateInequalityConstraints

        ! KKT system relaxation (not always used)
        procedure(RelaxProblemKKTSystemINT), deferred  :: &
            RelaxProblemKKTSystem

        ! Data writing per optimization iteration
        procedure(WriteIterationDataOptimizationProblemINT), deferred :: &
            WriteIterationData 

        ! Merit function evaluation
        procedure :: EvaluateMeritFunction
        procedure :: EvaluateMeritFunctionL1

        ! Lagrangian evaluation
        procedure :: EvaluateLagrangian

    end type

    ! Optimization solver
    !====================
    ! General abstract solver type
    type, abstract :: OptimizationSolverUDT

        ! Defines the optimization solver and its numerics. It should 
        ! have a generic 'solve' method which acts on the optimization
        ! problem
        character(:), allocatable                   :: inputfilepath
        character(:), allocatable                   :: inputfileprefix

        type(NumUDT)                                :: num

    contains 

        ! Initialization
        procedure(InitializeOptimizationSolverINT), deferred :: &
            Initialize

        ! Solution procedure
        procedure(SolveOptimizationProblemINT), deferred    :: &
            SolveOptimizationProblem 

    end type

    ! KKT solver
    type, extends(OptimizationSolverUDT) :: OptimizationSolverKKTUDT

        ! Description
        !============
        ! KKT system based optimization solver. Contains additional
        ! numerics fields and auxiliary routines
        type(NumKKTUDT)         :: numKKT
        type(numLSUDT)          :: numLS 
        type(NumNCPUDT)         :: numNCP

    contains 

        ! Initialization
        procedure ::    Initialize                  => InitializeKKTSolver

        ! Solution procedure using KKT solver
        procedure ::    SolveOptimizationProblem    => SolveOptimizationProblemKKT

        ! Convergence checking
        procedure :: CheckConvergenceKKT 

        ! Setup of correction equation
        procedure :: SetupCorrectionEquation

        ! Relaxation of KKT system
        procedure :: RelaxKKTSystem

    end type

    ! Quasi-Newton solver 
    type, extends(OptimizationSolverUDT) :: OptimizationSolverQNUDT

        ! Description
        !============  
        ! Quasi-Newton solver that uses approximate Hessian to compute
        ! the step direction. 

        ! Numerics
        type(NumQNUDT)          :: numQN
        type(numLSUDT)          :: numLS
        
        ! Hessian approximation
        !class(HessianApprximationUDT)   :: hessian

    contains

        ! Initialization
        procedure ::    Initialize                  => InitializeQNSolver

        ! Solution procedure using QN solver
        procedure ::    SolveOptimizationProblem    => SolveOptimizationProblemQN

        ! Convergence checking
        procedure :: CheckConvergenceQN 

        ! Hessian 

    end type

    ! Optimization engine
    !====================
    type OptimizationEngineUDT 

        ! Description
        !============
        ! The engine simply contains the problem and solver structures
        ! (and in the future the monitor). Also contains a 
        ! character array through which an input file can be specified
        ! for reading in numerics.
        character(:), allocatable                   :: inputfilepath
        character(:), allocatable                   :: inputfileprefix
        class(OptimizationProblemUDT), allocatable  :: problem
        class(OptimizationSolverUDT), allocatable   :: solver

    contains

        ! Main driver to solve a problem
        procedure :: Driver                 => OptimizationEngineDriver

    end type

    ! Diagnostics
    !============
    ! Cost function for diagnostics
    type, extends(DiagnosticsFunctionUDT) :: DFCostfunctionUDT 

        ! Description
        !============
        ! Function type for evaluating and checking the cost function.
        class(OptimizationProblemUDT), allocatable      :: problem 

    contains 

        ! Evaluation procedure
        procedure :: Evaluate       => EvaluateDFCostfunction 

        ! Dimension getter
        procedure :: GetDimensions  => GetProblemDimensions

        ! Argument getter
        procedure :: GetArguments   => GetProblemArguments

    end type

    ! Lagrangian for diagnostics
    type, extends(DiagnosticsFunctionUDT) :: DFLagrangianUDT 

        ! Description
        !============
        ! Function type for evaluating and checking the cost function.
        class(OptimizationProblemUDT), allocatable      :: problem 
        class(OptimizationSolverUDT), allocatable       :: solver

        ! Lagrange multipliers
        real(R8), allocatable, dimension(:)             :: lambda, mu 


    contains 

        ! Evaluation procedure
        procedure :: Evaluate       => EvaluateDFLagrangian 

        ! Dimension getter
        procedure :: GetDimensions  => GetProblemDimensionsLagrangian

        ! Argument getter
        procedure :: GetArguments   => GetProblemArgumentsLagrangian

    end type

    ! Eqcon for diagnostics
    type, extends(DiagnosticsFunctionUDT) :: DFEqconUDT 

        ! Description
        !============
        ! Function type for evaluating and checking the cost function.
        class(OptimizationProblemUDT), allocatable      :: problem 

        ! Lagrange multipliers
        real(R8), allocatable, dimension(:)             :: lambda
        
        ! Equality constraint ID
        integer(I8)                                     :: eqID

        ! Indicator to multiply with lambda
        logical                                         :: multlambda


    contains 

        ! Evaluation procedure
        procedure :: Evaluate       => EvaluateDFEqcon

        ! Dimension getter
        procedure :: GetDimensions  => GetProblemDimensionsEqcon

        ! Argument getter
        procedure :: GetArguments   => GetProblemArgumentsEqcon

    end type

    ! Ineqcon for diagnostics
    type, extends(DiagnosticsFunctionUDT) :: DFIneqconUDT 

        ! Description
        !============
        ! Function type for evaluating and checking the cost function.
        class(OptimizationProblemUDT), allocatable      :: problem 

        ! Lagrange multipliers
        real(R8), allocatable, dimension(:)             :: mu
        
        ! Inequality constraint ID
        integer(I8)                                     :: ineqID

        ! Indicator to multiply with mu
        logical                                         :: multmu


    contains 

        ! Evaluation procedure
        procedure :: Evaluate       => EvaluateDFIneqcon

        ! Dimension getter
        procedure :: GetDimensions  => GetProblemDimensionsIneqcon

        ! Argument getter
        procedure :: GetArguments   => GetProblemArgumentsIneqcon

    end type


    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    ! Optimization problem
    !=====================
    abstract interface

        ! Problem initialization
        subroutine InitializeINT(problem)
            
            ! Description
            !============
            ! This should be a general subroutine that initializes the 
            ! design, given only the problem as input. Possibly, this 
            ! problem can contain different substructures that allow 
            ! the design to be initialized. 

            ! Define interface
            !=================
            import :: OptimizationProblemUDT 
            class(OptimizationProblemUDT) :: problem

        end subroutine

        ! Get problem dimensions
        subroutine GetProblemDimensionsINT(problem, nphi, neq, nineq)

            ! Description
            !============
            ! This routine should return the general problem dimensions.
            ! Here, nphi is the number of design variables, neq the 
            ! number of equality constraints, and nineq the number of
            ! inequality constraints.
            
            ! Import
            import :: OptimizationProblemUDT, I8 

            ! Declare
            class(OptimizationProblemUDT)       :: problem 
            integer(I8), intent(out)            :: nphi,  neq, nineq

        end subroutine

        ! Get problem design variables
        subroutine GetProblemDesignVariablesINT(problem, phi)

            ! Description
            !============
            ! This routine should return the problem's design variables
            ! as a real array 'phi'

            ! Import
            import :: OptimizationProblemUDT, R8 

            ! Declare
            class(OptimizationProblemUDT)       :: problem 
            real(R8), allocatable, intent(out)  :: phi(:)

        end subroutine

        ! Design update
        subroutine UpdateDesignINT(problem, dx)

            ! Description
            !============
            ! This routine should update only the design variable vector
            ! with the increment dx. This increment should have the 
            ! correct dimensions

            ! Import
            import :: OptimizationProblemUDT, R8 

            ! Declare
            class(OptimizationProblemUDT)       :: problem 
            real(R8), intent(in)                :: dx(:)

        end subroutine

        ! Problem update
        subroutine UpdateProblemINT(problem)

            ! Description
            !============
            ! This routine should update all problem structures after
            ! the design step has been computed. 

            ! Import
            import :: OptimizationProblemUDT

            ! Declare
            class(OptimizationProblemUDT)       :: problem 

        end subroutine

        ! Cost function evaluation
        subroutine EvaluateCostFunctionINT(problem, J, gradJ, hessJ, &
            dogradient, dohessian)

            ! Description
            !============
            ! This routine returns the value of the cost function, J, 
            ! and (if dogradient is true) the gradient gradJ, and 
            ! (if dohessian is true) the hessian of J. In case the 
            ! the gradient or hessian shouldn't be computed, the values
            ! of gradJ and hessJ are garbage. Note that the hessian is 
            ! stored in a derived type with row, col and val vectors.
            ! This should be used later on to construct the actual 
            ! hessian in any suitable format. For this end, an external
            ! module 'mod_sparseinterface' is used. 
            
            ! Import
            import :: OptimizationProblemUDT, MySparseUDT, R8

            ! Declare
            class(OptimizationProblemUDT)       :: problem
            real(R8)                            :: J
            real(R8), allocatable               :: gradJ(:)
            type(MySparseUDT)                   :: hessJ 
            logical                             :: dogradient, dohessian

        end subroutine

        ! Equality constraints evaluation
        subroutine EvaluateEqualityConstraintsINT(problem, G, gradG, &
            hessG, dogradient, dohessian, lambda)

            ! Description
            !============
            ! This routine returns the equality constraints, G, the
            ! gradient, gradG (in mysparse format), and the hessian-
            ! vector multiplication with the vector lambda hessG (also 
            ! in mysparse format).

            ! Import
            import :: OptimizationProblemUDT, MySparseUDT, R8

            ! Declare
            class(OptimizationProblemUDT)       :: problem 
            real(R8), allocatable               :: G(:), lambda(:)
            type(MySparseUDT)                   :: gradG, hessG 
            logical                             :: dogradient, dohessian

        end subroutine

        ! Inequality constraints evaluation
        subroutine EvaluateInequalityConstraintsINT(problem, H, gradH, &
            hessH, dogradient, dohessian, mu)

            ! Description
            !============
            ! This routine returns the inequality constraints, H, the
            ! gradient, gradH (in mysparse format), and the hessian-
            ! vector multiplication with the vector mu hessG (also 
            ! in mysparse format).

            ! Import
            import :: OptimizationProblemUDT, MySparseUDT, R8

            ! Declare
            class(OptimizationProblemUDT)       :: problem 
            real(R8), allocatable               :: H(:), mu(:)
            type(MySparseUDT)                   :: gradH, hessH 
            logical                             :: dogradient, dohessian

        end subroutine
    
        ! KKT system relaxation
        subroutine RelaxProblemKKTSystemINT(problem, KKT)
            import :: MySparseUDT, OptimizationProblemUDT
            class(OptimizationProblemUDT)   :: problem 
            type(MySparseUDT)               :: KKT
        end subroutine

        ! Data output
        subroutine WriteIterationDataOptimizationProblemINT(problem, itopt)
            import :: OptimizationProblemUDT, I8
            class(OptimizationProblemUDT)   :: problem 
            integer(I8)                     :: itopt
        end subroutine

        ! Solver initialization
        subroutine InitializeOptimizationSolverINT(solver)
            import :: OptimizationSolverUDT
            class(OptimizationSolverUDT)    :: solver 
        end subroutine

        ! Solver driver
        subroutine SolveOptimizationProblemINT(solver, problem)
            import :: OptimizationSolverUDT, OptimizationProblemUDT 
            class(OptimizationSolverUDT)    :: solver 
            class(OptimizationProblemUDt)   :: problem 
        end subroutine

    end interface

    contains

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                       OPTIMIZATION ENGINE                        !
    !------------------------------------------------------------------!

    ! Main driver
    subroutine OptimizationEngineDriver(optimizationengine)

        ! Description
        !============
        ! This is the main generic driver of the optimization problem, 
        ! which normally should not change anymore, as all user specific
        ! data etc should be set through different routines. This driver
        ! sets up and solves an optimization problem of a specific type
        ! by calling the setup routines of the problem and solver, and 
        ! then calling the main drivers of these routines. 

        ! Notes
        !======
        ! Note 1: as currently only a KKT solver is implemented, these
        ! routines are called directly. Should be encapsulated in the 
        ! future. 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationEngineUDT) :: optimizationengine

        ! Initialize
        !===========
        !call optimizationengine%SetupOptimizationDriver()
        
        ! Set up the problem
        optimizationengine%problem%inputfilepath = optimizationengine%inputfilepath
        call optimizationengine%problem%Initialize()

        ! Set up the solver - now KKT
        print *, 'reading optimization engine input from file: ' // optimizationengine%inputfilepath
        optimizationengine%solver%inputfilepath = optimizationengine%inputfilepath
        optimizationengine%solver%inputfileprefix = optimizationengine%inputfileprefix
        call optimizationengine%solver%Initialize()

        ! Solve
        !======
        ! Solve the optimization problem by calling the KKT solver
        call optimizationengine%solver%SolveOptimizationProblem( &
            optimizationengine%problem)

    end subroutine

    !------------------------------------------------------------------!
    !                       OPTIMIZATION PROBLEM                       !
    !------------------------------------------------------------------!

    ! Merit function wrapper
    recursive subroutine EvaluateMeritFunction(problem, f, DJf, dx, lambda, mu, &
        doderiv, meritfunction, num)

        ! Description
        !============
        ! Wrapper for merit function evaluation. It is assumed that all
        ! data is given at the current (not updated) iterate, and that 
        ! the next iterate is achieved by updating the design, lambda, 
        ! and mu with dx. 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemUDT)   :: problem 
        real(R8)                        :: f, DJf
        real(R8)                        :: dx(:)
        real(R8), allocatable           :: lambda(:), mu(:)
        logical                         :: doderiv 
        character(*), intent(in)        :: meritfunction 
        type(numLSUDT)                  :: num

        ! Check which merit function to evaluate
        select case (meritfunction) 

        case ('l1')

            call problem%EvaluateMeritFunctionL1(f, DJf, dx, lambda, &
                mu, doderiv, num)

        case default 

            call gdErrorHandler('Unknown merit function type')

        end select

        ! Check if the result is NaN or inf - in that case, return inf
        if ( (.not. ieee_is_finite(f)) .or. (ieee_is_nan(f)) ) then 
            ! Set to infinity
            f = posinfval_R8()
        end if 

    end subroutine

    ! L1 merit function
    recursive subroutine EvaluateMeritFunctionL1(problem, f, DJf, dx, lambda, mu, &
        doderiv, num)

        ! Description
        !============
        ! Evaluate the L1 merit function. This typically suffers from
        ! Maratos effect in line searches, so this should be countered 
        ! (e.g. with the second order correction). We assume that all
        ! fields are already correctly updated (including lambda and mu)

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemUDT)   :: problem 
        real(R8)                        :: f, DJf
        real(R8)                        :: dx(:)
        real(R8), allocatable           :: lambda(:), mu(:), dl(:), dm(:)
        logical                         :: doderiv 
        type(numLSUDT)                  :: num

        ! Auxiliary
        logical                         :: dogradient, dohessian 
        integer(I8)                     :: nphi, neq, nineq
        real(R8)                        :: penfac, delta, sf, &
            maxvalmu, maxvallambda, gnorm, hnorm

        ! Cost function 
        real(R8)                    :: J 
        real(R8), allocatable       :: gradJ(:)
        type(MySparseUDT)           :: hessJ 

        ! Equality constraints 
        real(R8), allocatable       :: G(:)
        type(MySparseUDT)           :: gradG, hessG  

        ! Inequality constraints 
        logical, allocatable        :: A(:)
        real(R8), allocatable       :: H(:)
        type(MySparseUDT)           :: gradH, hessH

        ! Initialize
        !===========
        ! Get problem dimensions
        call problem%GetProblemDimensions(nphi, neq, nineq)

        ! Set scaling factor
        sf = 1e-4

        ! Extract lambda and mu updates
        allocate(dl(neq), dm(nineq))
        dl = dx(nphi+1:nphi+neq)
        dm = dx(nphi+neq+1:nphi+neq+nineq)

        ! Cost function 
        allocate(gradJ(nphi))
        J = 0
        gradJ(:) = 0
        hessJ%nrow = nphi 
        hessJ%ncol = nphi 

        ! Equality constraint 
        allocate(G(neq))
        G(:) = 0
        gradG%nrow = nphi 
        gradG%ncol = neq 
        hessG%nrow = nphi 
        hessG%ncol = nphi

        ! Inequality constraint 
        allocate(H(nineq), A(nineq))
        H(:) = 0
        A(:) = .false. 
        gradH%nrow = nphi 
        gradH%ncol = nineq 
        hessH%nrow = nphi 
        hessH%ncol = nphi

        ! Compute delta
        !==============
        ! Check for maxima
        maxvalmu        = 0
        maxvallambda    = 0
        if (nineq > 0) then 
            maxvalmu = maxval(abs(mu + dm)) ! just in case the updates lead to negative mu
        end if 
        if (neq > 0) then 
            maxvallambda = maxval(abs(lambda + dl)) 
        end if 

        ! Set constant
        delta = 1e-4*(max(maxval(lambda + dl), maxval(mu + dm)))
        if (delta <= 0) then 
            ! Set default
            delta = 1
        end if

        ! Compute cost function and constraints
        !======================================
        ! Check what we need to evaluate
        dogradient  = .true.
        dohessian   = .false.
        if (.not. doderiv) then 
            ! No gradients required
            dogradient = .false.
        end if 

        ! Cost function
        call problem%EvaluateCostFunction(J, gradJ, hessJ, dogradient, &
            dohessian)

        ! Equality constraints
        call problem%EvaluateEqualityConstraints(G, gradG, hessG, &
            dogradient, dohessian, lambda)

        ! Inequality constraints
        call problem%EvaluateInequalityConstraints(H, gradH, hessH, &
            dogradient, dohessian, mu)

        ! Determine which inequality constraints are active and should
        ! contribute
        A = H > 0

        ! Compute merit function
        !=======================
        ! Compute penalty factor
        penfac = max(maxvallambda, maxvalmu) + delta
        
        ! Compute L1 norm of constraints
        gnorm = 0
        hnorm = 0
        if (neq > 0) then 
            gnorm = sum(abs(G))
        end if 
        if (nineq > 0) then 
            hnorm = sum(abs(pack(H, A)))
        end if 

        ! Compute values
        f = J + penfac*(gnorm + hnorm)
        if (doderiv) then 
            if ((neq > 0) .or. (nineq > 0)) then 
                DJf = sum(gradJ*dx(1:nphi)) - penfac*(gnorm + hnorm)
            else 
                ! Unconstrained value
                DJf = sum(gradJ*dx(1:nphi))
            end if 
        end if 

        ! Housekeeping
        !=============

    end subroutine
    
    !------------------------------------------------------------------!
    !                     KKT OPTIMIZATION SOLVER                      !
    !------------------------------------------------------------------!

    ! KKT solver initialization
    subroutine InitializeKKTSolver(solver)

        ! Description
        !============
        ! Initialize the necessary fields and structures for the KKT 
        ! solver. This is basically only the numerics - the optimization
        ! problem itself should already be initialized beforehand!

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationSolverKKTUDT)             :: solver    

        ! Data

        ! Initialize
        !===========
        ! Numerics
        solver%num%inputfilepath    = solver%inputfilepath
        solver%num%fieldprefix      = solver%inputfileprefix
        solver%numKKT%inputfilepath = solver%inputfilepath
        solver%numKKT%fieldprefix   = solver%inputfileprefix 
        solver%numLS%inputfilepath  = solver%inputfilepath
        solver%numLS%fieldprefix    = solver%inputfileprefix
        solver%numNCP%inputfilepath = solver%inputfilepath
        solver%numNCP%fieldprefix   = solver%inputfileprefix
        call solver%num%InitializeNumParams()
        call solver%numKKT%InitializeNumParams() 
        call solver%numLS%InitializeNumParams()
        call solver%numNCP%InitializeNumParams()

    end subroutine

    ! KKT solver
    recursive subroutine SolveOptimizationProblemKKT(solver, problem)

        ! Description
        !============
        ! KKT solver for the optimization problem defined by the generic
        ! 'problem'. It is assumed that the optimization problem is
        ! properly initialized. The solver is initialized here. 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationSolverKKTUDT)             :: solver    
        class(OptimizationProblemUDT)               :: problem
    
        ! Loop variables
        integer(I8)                 :: itopt, maxit, k
        logical                     :: converged
        
        ! Auxiliary variables 
        real(R8)                    :: convnorm
        logical                     :: dogradient, dohessian 

        ! Cost function 
        real(R8)                    :: J 
        real(R8), allocatable       :: gradJ(:)
        type(MySparseUDT)           :: hessJ 

        ! Equality constraints 
        real(R8), allocatable       :: G(:), lambda(:) 
        type(MySparseUDT)           :: gradG, hessG  

        ! Inequality constraints 
        real(R8), allocatable       :: H(:), mu(:) 
        type(MySparseUDT)           :: gradH, hessH  

        ! nonlinear complementarity function
        real(R8), allocatable       :: ncp(:) 
        type(MySparseUDT)           :: gradncpphi, gradncpmu
        logical, allocatable        :: A(:), I(:) 

        ! Lagrangian 
        real(R8)                    :: L 
        real(R8), allocatable       :: gradL(:)
        type(MySparseUDT)           :: hessL 

        ! Solver & updates
        double precision, allocatable :: dx(:), dxl(:)
        real(R8), allocatable       :: rhs(:)
        real(R8)                    :: alphals
        integer(I8)                 :: flag, flagls
        integer(I8), allocatable    :: phiind(:), eqconind(:), &
            ineqconind(:), activeineqconind(:), inactiveineqconind(:)
        type(MySparseUDT)           :: hessLJ, lhs, hessLC

        ! FD checkers
        type(FDcheckerUDT)          :: FDcfv, FDeqcon, FDineqcon

        ! Diagnostics
        integer                     :: errstat

        ! Timing
        real(R8)                    :: t_it_s, t_it_e, &
            t_eval_s, t_eval_e, t_linsolve_s, t_linsolve_e

        ! Data

        ! Temporary variables (to be deleted in the future)
        integer(I8)                 :: nphi, neq, nineq

        ! Initialize & unpack
        !====================
        ! Logicals
        dogradient  = .true. 
        dohessian   = .true. 

        ! Initialize the monitor - only temporary here
        call problem%GetProblemDimensions(nphi, neq, nineq)
        call problem%monitor%Initialize(solver%num%maxit, nphi, neq,&
            nineq, solver%num%tol)

        ! Initialize FD checkers
        call FDcfv%Initialize(solver%numKKT%checkcfvvars, &
            solver%numKKT%FDsteps, solver%numKKT%checkoutputfile // '_cfv')
        call FDeqcon%Initialize(solver%numKKT%checkeqconvars, &
            solver%numKKT%FDsteps, solver%numKKT%checkoutputfile // '_eqcon')
        call FDineqcon%Initialize(solver%numKKT%checkineqconvars, &
            solver%numKKT%FDsteps, solver%numKKT%checkoutputfile // '_ineqcon')

        ! Cost function 
        allocate(gradJ(nphi))
        J = 0
        gradJ(:) = 0
        hessJ%nrow = nphi 
        hessJ%ncol = nphi 

        ! Equality constraint 
        allocate(G(neq), lambda(neq))
        G(:) = 0
        lambda(:) = 0.0
        gradG%nrow = nphi 
        gradG%ncol = neq 
        hessG%nrow = nphi 
        hessG%ncol = nphi

        ! Inequality constraint 
        allocate(H(nineq), mu(nineq))
        H(:) = 0
        mu(:) = 0
        gradH%nrow = nphi 
        gradH%ncol = nineq 
        hessH%nrow = nphi 
        hessH%ncol = nphi

        ! NCP function
        allocate(ncp(nineq), A(nineq), I(nineq))
        gradncpphi%nrow = nphi 
        gradncpphi%ncol = nineq
        gradncpmu%nrow = nineq 
        gradncpmu%ncol = nineq
        A(:) = .false.
        I(:) = .not. A

        ! Lagrangian 
        allocate(gradL(nphi + neq + nineq))
        hessL%nrow = nphi + neq + nineq
        hessL%ncol = nphi + neq + nineq
        L = 0
        gradL(:) = 0

        ! Solver
        allocate(dx(hessL%nrow))
        allocate(rhs(hessl%nrow))
        dx = 0
        rhs = 0
        allocate(phiind(nphi), eqconind(neq), ineqconind(nineq))
        phiind = [(k, k = 1, nphi)]
        eqconind = [(k, k = nphi+1, nphi+neq)]
        ineqconind = [(k, k = nphi+neq+1, nphi+neq+nineq)]

        ! Diagnostics
        !checkgradients  = .false. ! check gradients in each iteration?
        !checkhessians   = .false. ! check hessians in each iteration?

        ! Initialize counter(s)
        itopt = 1
        problem%monitor%itopt = itopt
        maxit = solver%num%maxit

        ! Unpack 
        associate(&
            rxf         => solver%numKKT%rxf, &
            rxfdesign   => solver%numKKT%rxfdesign, &
            rxfdec      => solver%numKKT%rxfdec, &
            rxfmin      => solver%numKKT%rxfmin, &
            verbosity   => solver%numKKT%verbosity, &
            num         => solver%numKKT)

        ! Main loop
        !==========
        ! Set convergence
        converged = .false. 

        ! Print solver header
        if (verbosity > 0) then
            ! Print out the header
            call problem%monitor%PrintHeader()

        end if

        ! Loop
        do while ( (.not. converged) .and. (itopt <= maxit))

            ! Start timing
            call wall_time(t_it_s)
            call wall_time(t_eval_s)

            ! Update the monitor
            problem%monitor%itopt = itopt
            problem%monitor%rxf = rxf 

            ! Start tracking errors
            call ErrorStack%StartTrack()

            ! Update the optimization problem 
            call problem%UpdateProblem()

            ! Check gradients & hessians if needed
            if (num%checkcfvgradient .or. num%checkcfvhessian) then 
                call CheckCostFunctionLinearization(problem, FDcfv, &
                    num%checkcfvgradient, num%checkcfvhessian)
            end if 
            if (num%checkeqcongradient .or. num%checkeqconhessian) then 
                call CheckEqconLinearization(problem, FDeqcon, lambda, &
                    num%checkeqcongradient, num%checkeqconhessian, &
                    num%checkeqconeqs)
            end if 
            if (num%checkineqcongradient .or. num%checkineqconhessian) then 
                call CheckIneqconLinearization(problem, FDineqcon, mu, &
                    num%checkineqcongradient, num%checkineqconhessian, &
                    num%checkineqconeqs)
            end if 

            !if (checkgradients .or. checkhessians) then 
                !call CheckLagrangianLinearization(problem, solver, lambda, & 
                !    mu, checkgradients, checkhessians)
            !end if

            ! Evaluate cost function
            call problem%EvaluateCostFunction(J, gradJ, & 
                hessJ, dogradient, dohessian)

            ! Evaluate the equality constraints
            call problem%EvaluateEqualityConstraints(G, gradG, &
                hessG, dogradient, dohessian, lambda)

            ! Evaluate the inequality constraints
            call problem%EvaluateInequalityConstraints(H, gradH, &
                hessH, dogradient, dohessian, mu)

            ! Evaluate the nonlinear complementarity function 
            call EvaluateNCPfunction(ncp, A, I, gradncpphi, gradncpmu, &
                H, gradH, mu, solver%numNCP, dogradient)

            ! Evaluate the Lagrangian
            call problem%EvaluateLagrangian(L, gradL, hessL, &
                J, gradJ, hessJ, &
                G, gradG, hessG, lambda, &
                H, gradH, hessH, mu, A, &
                dogradient, dohessian)

            ! Check if an error was encountered
            errstat = ErrorStack%ErrorState()
            call ErrorStack%EndTrack()
            if (errstat > 0) then 
                ! Call error, exit the loop
                call gdErrorHandler('SolveOptimizationProblemKKT: could ' // &
                    'not evaluate problem, exiting', severityin=0)
                exit 
            end if 

            ! Check if NaNs are encountered in residual
            if (any(isnan(gradL))) then 
                ! Call error, exit the loop
                call gdErrorHandler('SolveOptimizationProblemKKT: NaNs ' // &
                    'detected when evaluating the problem, exiting', severityin=0)
                exit 
            end if 

            ! Check convergence
            call solver%CheckConvergenceKKT(gradL, converged, convnorm)

            ! Timers
            call wall_time(t_eval_e)

            ! Solve 
            if (.not. converged) then 

                ! Set up correction equation
                call solver%SetupCorrectionEquation(lhs, rhs, &
                    gradJ, hessJ, &
                    G, gradG, hessG, lambda, &
                    H, gradH, hessH, mu, A, &
                    ncp, gradncpphi, gradncpmu)

                ! Relax
                if (solver%numKKT%useproblemrelaxation) then 
                    ! Call problem-specific KKT relaxation routine
                    call problem%RelaxProblemKKTSystem(lhs)
                else
                    ! Use default KKT solver routine
                    call solver%RelaxKKTSystem(lhs, nphi, neq, nineq)
                end if 

                !print *, 'hessian size: ', hessL%nrow, hessL%ncol
                !call SpyPlot(lhs%row, lhs%col, lhs%nval, '-p')
                
                ! Call the sparse solver
                call wall_time(t_linsolve_s)
                call SolveSparseLinearSystemDI(lhs, rhs, dx, flag)
                call wall_time(t_linsolve_e)

            else
                ! Don't solve again, already converged. Exit 
                dx(:) = 0
            end if

            ! Do linesearch?
            alphals = 1
            if (solver%numLS%dolinesearch) then 
                ! Compute the step length for the line search, don't 
                ! apply relaxation using rxfdesign. Note: also the 
                ! Lagrange multipliers may change!
                if (flag == 0) then 

                    call ComputeStepLengthLS(problem, solver%numLS, dx, lambda, mu, alphals, flagls) ! dx is changed during linesearch
                else 
                    ! Something wrong during linear solver, try with relaxation
                    flagls = 1
                end if 

                ! Check the linesearch output
                if (flagls == 0) then 
                    ! All good

                elseif (flagls == 1) then 
                    ! Non-descent direction, print message and skip remainder of iterate
                    if (flag /= 0) then 
                        print *, 'step direction computation not succeeded, ' // &
                            'reattempting with damped Hessian'

                    else

                        print *, 'non-descent direction, skipping update ' // &
                            'and reattempt with damped Hessian'

                    end if 

                    ! Set step to zero
                    dx(:) = 0
                    alphals = 0

                    ! Add relaxation
                    if (rxf > 0) then 
                        rxf = 2*rxf 
                    else
                        ! Apparently no relaxation, add
                        print *, 'No relaxation detected, adding relaxation'
                        rxf = 1
                        if (rxfdec > 0) then 
                        else 
                            rxfdec = 0.9
                        end if 
                    end if 
                end if 

                ! Update lagrange multipliers using least-squares approach
                ! for active constraints
                if ( (flag == 0) .and. (flagls == 0)) then 
                    allocate(activeineqconind(count(A)), inactiveineqconind(count(I)))
                    activeineqconind = pack(ineqconind, A)
                    inactiveineqconind = pack(ineqconind, I)
                    hessLJ = lhs%DeleteColumns([eqconind, ineqconind])
                    hessLJ = hessLJ%DeleteRows([eqconind, ineqconind])
                    hessLC = lhs%DeleteColumns([phiind, inactiveineqconind])
                    hessLC = hessLC%DeleteRows([eqconind, ineqconind])
                    allocate(dxl(neq + count(A)))
                    !call SolveSparseLinearSystemDI((gradG%Transpose()*gradG), &
                    !    MatrixVectorProduct(gradG%Transpose(), (gradL(1:nphi) + MatrixVectorProduct(hessLJ, dx(1:nphi)))), &
                    !    dxl2, flag)
                    !dxl2 = -dxl2
                    call SolveSparseLinearSystemDI((hessLC%Transpose()*hessLC), &
                        MatrixVectorProduct(hessLC%Transpose(), &
                        (rhs(phiind) - MatrixVectorProduct(hessLJ, dx(phiind)))), &
                        dxl, flag)

                    !print *, maxval(abs(dx(nphi+1:nphi+neq+nineq) - dxl))
                    
                    dx(nphi+1:nphi+neq) = dxl(1:neq)
                    dx(activeineqconind) = dxl(neq+1:neq+count(A))
                    deallocate(dxl, activeineqconind, inactiveineqconind)
                end if 
                
            else
                ! Check convergence of solver
                if (flag == 0) then 
                    ! Directly update design without linesearch, apply the
                    ! relaxation using rxfdesign
                    dx(1:nphi) = rxfdesign*dx(1:nphi)
                else
                    ! Set step to zero
                    dx(:) = 0
                    alphals = 0

                    ! Add relaxation
                    if (rxf > 0) then 
                        rxf = 2*rxf 
                    else
                        ! Apparently no relaxation, add
                        print *, 'No relaxation detected, adding relaxation'
                        rxf = 1
                        if (rxfdec > 0) then 
                        else 
                            rxfdec = 0.9
                        end if 
                    end if 
                end if 
            end if

            ! Update the design
            call problem%UpdateDesign(dx(1:nphi))

            ! Update lagrange multipliers
            lambda  = lambda + dx(nphi+1:nphi+neq)
            mu      = mu + dx(nphi+neq+1:nphi+neq+nineq)

            ! Set mu of non-active constraints to zero
            where (.not. A) mu = 0

            ! Update problem multipliers
            problem%lambda  = lambda 
            problem%mu      = mu 
            problem%A       = A
            
            ! Timers
            call wall_time(t_it_e)

            ! Update the monitor again
            problem%monitor%J(itopt)        = J
            problem%monitor%L(itopt)        = L
            problem%monitor%G(itopt)        = maxval(abs(G))
            problem%monitor%H(itopt)        = maxval(H)
            problem%monitor%alpha(itopt)    = alphals
            problem%monitor%convnorm(itopt) = convnorm
            problem%monitor%evaltime        = t_eval_e - t_eval_s
            problem%monitor%ittime          = t_it_e - t_it_s 
            problem%monitor%linsolvetime    = t_linsolve_e - t_linsolve_s
            problem%monitor%maxdphi = maxval(dx(1:nphi))

            ! Print the current iterate
            if (verbosity > 0) then 
                ! Print out the iterate
                call problem%monitor%PrintIterate()

            end if

            ! Write out the problem data
            if (verbosity > 1) then 
                call problem%WriteIterationData(itopt)
            end if 

            ! Update the iteration counter
            itopt = itopt+1

            ! Update the hessian relaxation factor
            rxf = rxf*rxfdec
            rxf = max(rxf, rxfmin)

        end do

        ! Housekeeping
        deallocate(G, H, gradJ, lambda, mu)
        end associate 
        
    end subroutine

    ! KKT system relaxation
    subroutine RelaxKKTSystem(solver, hessL, nphi, neq, nineq)

        ! Description
        !============
        ! This routine applies a relaxation procedure on the hessian of
        ! the KKT system, which should be given by hessL (in mysparse
        ! format). The relaxation is based on the absolute value of the
        ! sum of the columns of the hessian. 

        ! We relax here only the design variable part of the hessian, 
        ! i.e. the d2L/dphi2 part. Here, a diagonal matrix is added that
        ! is equal to rxf * sum(abs(hessL), 2) on the diagonal. The 
        ! relaxation factor rxf should be given in the numerics of the
        ! solver object. 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationSolverKKTUDT)     :: solver 
        type(MySparseUDT)                   :: hessL, hessL_copy
        integer(I8)                         :: nphi, neq, nineq

        ! Auxiliary
        type(MySparseUDT)                   :: hessL_temp, hessrelax
        real(R8), allocatable               :: diag(:)

        ! Loop
        integer(I8)                         :: k


        ! Construct relaxator
        !====================
        ! Initialize & allocate
        hessL_copy = hessL
        hessrelax%nrow = nphi + neq + nineq 
        hessrelax%ncol = hessrelax%nrow 
        hessrelax%nval = nphi 
        call hessrelax%Allocate()
        allocate(diag(hessrelax%nrow))
        
        ! Set the diagonal elements
        hessL_copy%val = abs(hessL_copy%val)
        call hessL_copy%SumColumnwiseFull(diag)

        ! Set to one if diagonal is zero
        where (diag == 0.0) diag= 1.0

        ! Set the diagonal elements
        hessrelax%val = solver%numKKT%rxf*diag(1:nphi)
        hessrelax%row = [(k, k=1, nphi)]
        hessrelax%col = hessrelax%row

        ! Relax
        !======
        ! Simply add up
        hessL_temp = hessrelax + hessL

        ! And assign output
        hessL = hessL_temp

        ! Deallocate
        deallocate(diag)




    end subroutine

    ! Convergence checker
    subroutine CheckConvergenceKKT(solver, gradient, converged, infnorm)
        
        ! Description
        !============
        ! Check if the optimization problem is converged by checking
        ! the infinity norm of the Lagrange gradient. Note that it is
        ! assumed that any projection on box constraints, or accounting
        ! for other constraints, has been done already. 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationSolverKKTUDT)     :: solver 
        real(R8), intent(in)                :: gradient(:)
        logical, intent(inout)              :: converged 
        real(R8), intent(out)               :: infnorm

        ! Check convergence
        !==================
        ! Compute infinity norm
        infnorm = maxval(abs(gradient))

        ! Compare
        converged = infnorm < solver%num%tol

    end subroutine

    ! Lagrangian 
    subroutine EvaluateLagrangian(problem, L, gradL, hessL, J, gradJ, &
        hessJ, G, gradG, hessG, lambda, H, gradH, hessH, mu, A, &
        dogradient, dohessian)

        ! Description
        !============
        ! Evaluate the lagrangian, its gradient and hessian based on
        ! the cost function, constraints, and their multipliers. The 
        ! lagrangian is evaluated as follows:
        !
        !   L = J + sum(lambda(i) * G(i)) + sum(mu(j) * H(j), if A(j))
        !
        ! where A denotes the active inequality constraints. The 
        ! Jacobian is then as follows (only computed if dogradient
        ! is true, transpose is returned):
        !
        !   dLdphi      = dJdphi + sum(lambda(i) * dG(i)dphi) 
        !               + sum(mu(j), dH(j)dphi, A(j))
        !   dLdlambda   = transp(G)
        !   dLdmu       = transp(H, A) ! only active constraints
        !
        ! And the hessian is equal to (only computed if dohessian is
        ! true):
        !
        !   dLdphi2         = dJdphi2 + sum(lambda(i) * dG(i)dphi2) 
        !                   + sum(mu(j), dH(j)dphi2, A(j))
        !   dLdphidlambda   = dGdphi
        !   dLdphidmu       = dHdphi ! only for active constraints 
        !   dLdlambdadphi   = transp(dGdphi)
        !   dLdlambda2      = 0 
        !   dLdlambdadmu    = 0  
        !   dLdmdphi        = transp(dHdphi)
        !   dLdmdlambda     = 0
        !   dLdm2           = 0
        !
        ! It is assumed that the hessian vector product lambda * dGdphi2
        ! and mu * dHdphi2 is given, the latter one accounting for the
        ! activeness of the constraints. 

        ! Notes
        !======
        ! Note 1: it is assumed that the hessian and gradient of the
        ! inequality constrains are correctly adjusted for activeness
        ! of the constraints! To be sure, mu is copied into a local 
        ! variable and set to zero where A is false. For the Lagrangian
        ! gradient, also the values where A is false are set to zero. 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemUDT)       :: problem 

        real(R8), intent(inout)             :: L, gradL(:) 
        type(MySparseUDT), intent(inout)    :: hessL 

        real(R8), intent(in)                :: J, G(:), H(:), mu(:), &
                                            lambda(:), gradJ(:)
        type(MySparseUDT), intent(in)       :: hessJ, gradG, hessG, &
                                            gradH, hessH
        logical, intent(in)                 :: A(:), dogradient, &
                                            dohessian 

        ! Loop variables
        integer(I8)                         :: k 

        ! Auxiliary variables
        real(R8), allocatable               :: tmu(:), tH(:), val(:)
        integer(I8), allocatable            :: row(:), col(:)
                                        
        ! Compute Lagrangian
        !===================
        ! To be sure: set mu(.not. A) equal to zero
        allocate(tmu(size(mu)), tH(size(H)))
        tmu = mu 
        tH = H 
        where (.not. A) 
            tmu = 0
            tH = 0
        end where
       
        L = J + sum(lambda * G) + sum(tmu * H)

        ! Compute gradient
        !=================
        if (dogradient) then 

            ! Cost function contribution
            gradL(1:size(gradJ)) = gradJ(:) + &
                gradG%MatrixVectorProduct(lambda) + &
                gradH%MatrixVectorProduct(tmu)

            ! Equality constraints contribution
            gradL(size(gradJ)+1:size(G)+size(gradJ)) = G(:)

            ! Inequality constraints contribution
            gradL(size(gradJ)+size(G)+1:size(gradL)) = tH(:)

        end if 

        ! Compute hessian
        !================
        if (dohessian) then 

            ! Initialize & allocate
            k = 0
            hessL%nval  = hessJ%nval + hessG%nval + hessH%nval &
                        + 2*gradG%nval + 2*gradH%nval   

            allocate(row(hessL%nval), col(hessL%nval), val(hessL%nval))

            ! dLdpsi2
            !--------
            ! Cost function contribution
            val(k+1:k+hessJ%nval) = hessJ%val
            row(k+1:k+hessJ%nval) = hessJ%row 
            col(k+1:k+hessJ%nval) = hessJ%col 
            k = k + hessJ%nval 

            ! Equality constraints contribution
            val(k+1:k+hessG%nval) = hessG%val
            row(k+1:k+hessG%nval) = hessG%row 
            col(k+1:k+hessG%nval) = hessG%col 
            k = k + hessG%nval 

            ! Inequality constraints contribution
            val(k+1:k+hessH%nval) = hessH%val
            row(k+1:k+hessH%nval) = hessH%row 
            col(k+1:k+hessH%nval) = hessH%col 
            k = k + hessH%nval 

            ! dLdlambdadphi
            !--------------
            val(k+1:k+gradG%nval) = gradG%val 
            row(k+1:k+gradG%nval) = gradG%row
            col(k+1:k+gradG%nval) = gradG%col + hessJ%nrow
            k = k + gradG%nval 

            ! dLdphidlambda 
            !--------------
            val(k+1:k+gradG%nval) = gradG%val 
            row(k+1:k+gradG%nval) = gradG%col + hessJ%nrow
            col(k+1:k+gradG%nval) = gradG%row  
            k = k + gradG%nval 

            ! dLdmudphi
            !----------
            val(k+1:k+gradH%nval) = gradH%val 
            row(k+1:k+gradH%nval) = gradH%row
            col(k+1:k+gradH%nval) = gradH%col + hessJ%nrow + gradG%ncol
            k = k + gradH%nval 

            ! dLdphidmu
            !----------
            val(k+1:k+gradH%nval) = gradH%val 
            row(k+1:k+gradH%nval) = gradH%col + hessJ%nrow + gradG%ncol
            col(k+1:k+gradH%nval) = gradH%row
            k = k + gradH%nval 

            ! Construct 
            !==========
            hessL = ConstructMySparse(row, col, val, hessL%nrow, hessL%ncol)

        end if

    end subroutine

    ! Nonlinear complementarity function 
    subroutine EvaluateNCPfunction(ncp, A, I, gradncpphi, gradncpmu, H, &
        gradH, mu, num, dogradient)

        ! Description
        !============
        ! NonlinearComplementarityFunction returns the nonlinear complementarity
        ! (problem)
        ! function and its derivatives with respect to the consistency and internal
        ! variables evaluated at the current state. The nonlinear complementarity
        ! function is used to enforce the complementarity condition, mu >= 0,
        ! in combination with the inequality constraints, h <= 0, i.e. mu*h = 0. For
        ! numerical reasons, this is relaxed to the following ncp:
        !
        !       ncp = max(alpha*h + mu,0) - mu = 0,
        !
        ! where mu is the lagrange multiplier, and alpha an arbitrarily yet
        ! strictly positive constant. Clearly, this function is not continuous, and
        ! its linearization depends on the outcome of the max operator. Therefore,
        ! the active and inactive sets are determined (active if mu + alpha*h >
        ! 0).

        ! Input
        !======
        ! - h:          inequality constraint values at the current design point
        !               (nh-by-1)
        ! - mu:         current estimated value of the lagrange multipliers
        !               (nh-by-1)
        ! - alpha:      strictly positive constant (scalar)
        ! - gradh:      (required if more than 3 output arguments) gradient of the
        !               inequality constraint values at the current design point
        !               (nh-by-nphi)

        ! Output
        !=======
        ! - ncp:        the residuals of the ncp function, i.e.
        !               max(alpha*f + lambda,0) - lambda evaluated at the current
        !               iterate (nh-by-1)
        ! - gradncp:    structure containing linearization of the ncp function with
        !               respect to the design and lagrange multipliers lambda and
        !               mu (the second one is always the zero matrix).
        ! - A,I:        the active (A) and inactive (I) sets.

        ! Notes
        !======
        ! Note 1: different implementations of the ncp function are available, yet
        ! the 'max' one is recommended. Should be set in SetNumParams, otherwise
        ! default 'max' is set.

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)                :: H(:), mu(:)
        logical                             :: dogradient
        logical, intent(inout)              :: A(size(H, 1)), &
            I(size(H, 1))
        type(MySparseUDT), intent(in)       :: gradH 
        type(MySparseUDT), intent(inout)    :: gradncpphi, gradncpmu 
        real(R8), intent(out)               :: ncp(size(H, 1))
        type(NumNCPUDT), intent(in)         :: num

        ! Auxiliary
        real(R8)                            :: temp(size(H, 1), 2)

        integer(I8)                         :: nh

        ! Loop
        integer(I8)                         :: k

        ! Initialize
        !===========
        nh = size(H, 1)

        ! Ncp function value
        !===================
        select case (num%ncpfun)
            
        case ('max')
            ! Based on maximal value
            temp(:, 1) = num%alpha*H + mu
            temp(:, 2) = 0  
            ncp = maxval(temp, 2) - mu
            
            ! Active and inactive sets
            A = (num%alpha*H + mu) >= 0
            I = .not. A
            
        case ('FB')

            ! Based on non-smooth Fischer-Burmeister function
            ncp = -(sqrt(H**2 + mu**2) - mu + H);
            
            ! Active and inactive sets
            A = .true. ! not required
            I = .false.
            
        case ('FBsmooth')

            ! Based on smoothed Fischer-Burmeister function
            ncp = -(sqrt((H)**2 + mu**2 + 2*num%alpha) - mu + (H));
            
            ! Active and inactive sets
            A = .true. ! not required
            I = .false.
            
        case default
            
            call gdErrorHandler('NonlinearComplementarityFunction: unknown ncp function')
            
        end select


        ! Ncp function derivative
        !========================
        if (dogradient) then 

            ! Initialize
            gradncpphi = gradH 
            gradncpmu%nrow = nh
            gradncpmu%ncol = nh
            gradncpmu%nval = nh ! diagonal matrix basically
            call gradncpmu%Allocate()
            gradncpmu%row = [(k, k = 1, nh)]
            gradncpmu%col = [(k, k = 1, nh)]
            
            select case (num%ncpfun)
                
            case ('max')

                ! Adjust values
                where (A(gradncpphi%col)) 
                    gradncpphi%val = num%alpha*gradH%val 
                elsewhere 
                    gradncpphi%val = 0
                end where
                where (A) 
                    gradncpmu%val = 0
                elsewhere 
                    gradncpmu%val = -1
                end where
                
            case ('FB')
                
                call gdErrorHandler('derivative not yet implemented')
                
            case ('FBsmooth')
                
                call gdErrorHandler('derivative not yet implemented')

            case default
                
                ! Error thrown above
            end select

        end if

    end subroutine

    ! Step length computation
    recursive subroutine ComputeStepLengthLS(problem, numLS, dx, lambda, mu, alpha, flag)

        ! Description
        !============
        ! Compute the step length to take, given the step direction dx.
        ! It is assumed that all necessary data etc can be derived from
        ! the problem definition (e.g. dimensions of design variables).
        ! The step length computation is typically done using a 
        ! linesearch approach.

        ! The line search method should be defined in numLS%type
        ! and can be 'backtracking', 'wolfe', 'backtracking_soc'. The 
        ! last one also applies a second-order correction. Note that 
        ! no constraints are explicitly accounted for in this line 
        ! search routine! Use an appropriate merit function for this
        ! (see also the EvaluateMeritFunction subroutine)

        ! Modules
        !========
        use ieee_arithmetic

        ! Declare variables
        !==================
        ! Arguments
        type(numLSUDT)                      :: numLS  
        class(OptimizationProblemUDT)       :: problem 
        real(R8), allocatable               :: dx(:), lambda(:), mu(:), dphi(:)
        real(R8)                            :: alpha
        integer(I8)                         :: flag 
        
        ! Auxiliary
        logical                             :: conv, doderiv 
        integer(I8)                         :: nphi, neq, nineq, na
        real(R8)                            :: f0, DJf0, fk, DJfk, &
            alpha_bot, alpha_top

        real(R8), allocatable               :: x0(:), x(:)
        real                                :: inf 

        ! Loop
        integer(I8)                         :: itls

        ! Variables for second order correction
        real(R8), allocatable               :: G(:), H(:), ck(:), &
            wkt(:), wk(:)
        type(MySparseUDT)                   :: gradG, gradH, hessG, &
            hessH, Ak, LSA  
        logical                             :: dogradient, dohessian 
        logical, allocatable                :: A(:)
        integer(I8)                         :: flag2
        integer                             :: errstat

        ! Initialize
        !===========
        ! Set flag 
        flag = 0

        ! Set convergence parameters
        conv = .false. 
        itls = 0

        ! Set infinity
        inf = ieee_value(inf, ieee_positive_inf)

        ! Store initial design point
        call problem%GetProblemDesignVariables(x0)

        ! Get problem parameters
        call problem%GetProblemDimensions(nphi, neq, nineq)
        
        ! Unpack design update
        allocate(dphi(nphi))
        dphi = dx(1:nphi)

        ! Initialize second order correction variables
        allocate(G(neq), H(nineq))
        gradG%nrow = nphi 
        gradG%ncol = neq 
        hessG%nrow = nphi 
        hessG%ncol = nphi
        gradH%nrow = nphi 
        gradH%ncol = nineq
        hessH%nrow = nphi 
        hessH%ncol = nphi
        dogradient  = .true. 
        dohessian   = .false.

        ! Associate
        associate(&
            maxit           => numLS%maxit, &
            c1              => numLS%c1, &
            c2              => numLS%c2, &
            dec             => numLS%dec, &
            inc             => numLS%inc, &
            meritfunction   => numLS%meritfunction)

        ! Descent check
        !==============
        ! Evaluate merit function and directional derivative
        doderiv = .true. 
        call problem%EvaluateMeritFunction(f0, DJf0, dx, lambda, mu, &
            doderiv, meritfunction, numLS)

        ! If no descent, exit with flag 1
        if (DJf0 >= 0) then 
            flag = 1
            return 
        end if 

        ! Step length
        !============
        select case (numLS%type)

        case ('backtracking')

            ! Don't compute any derivatives
            doderiv = .false.

            ! Loop
            do while ( (.not. conv) .and. (itls <= maxit) )
            
                ! Update current iterate
                x = x0 + alpha*dphi

                ! Start tracking for possible problems
                call ErrorStack%StartTrack()

                ! Update the design
                call problem%UpdateDesign(alpha*dphi)

                ! Update the problem
                call problem%UpdateProblem()

                ! Check error status
                errstat = ErrorStack%ErrorState()
                call ErrorStack%EndTrack()
                if (errstat > 0) then 
                    ! Error encountered, set value to infinity - don't
                    ! bother trying to compute the merit function
                    fk = inf
                else
                    ! Calculate new cost function value
                    call problem%EvaluateMeritFunction(fk, DJfk, dx, lambda, &
                        mu, doderiv, meritfunction, numLS)
                end if

                
                ! Check Armijo condition
                if (fk < f0 + c1*alpha*DJf0) then 
                    
                    ! Sufficient decrease, terminate
                    conv = .true.
                    
                else
                    
                    ! Decrease alpha
                    alpha = dec*alpha
                    
                end if 
                
                ! Update counter
                itls = itls + 1

                ! De-update the design
                call problem%UpdateDesign(x0-x)

            end do
            
            ! Checks
            if (itls-1 == maxit) then 
                ! Print message, set flag
                if (numLS%verbosity > 0) then 
                    print *, 'linesearch did not converge'
                end if 
            end if 
            

        case ('wolfe')

            ! Compute derivatives
            doderiv = .true. 

            ! Bounds for alpha (hard coded here)
            alpha_bot = 0
            alpha_top = inf 

            ! Loop
            do while ( (.not. conv) .and. (itls <= maxit) )
            
                ! Update current iterate
                x = x0 + alpha*dphi

                ! Start tracking for possible problems
                call ErrorStack%StartTrack()

                ! Update the design
                call problem%UpdateDesign(alpha*dphi)

                ! Update the problem
                call problem%UpdateProblem()

                ! Check error status
                errstat = ErrorStack%ErrorState()
                call ErrorStack%EndTrack()
                if (errstat > 0) then 
                    ! Error encountered, set value to infinity - don't
                    ! bother trying to compute the merit function
                    fk = inf
                else
                    ! Calculate new cost function value
                    call problem%EvaluateMeritFunction(fk, DJfk, dx, lambda, &
                        mu, doderiv, meritfunction, numLS)
                end if
                

                
                ! Check Armijo & Wolfe conditions
                if (fk > f0 + c1*alpha*DJf0) then 
                    
                    ! Decrease alpha
                    alpha_top   = alpha 
                    alpha       = dec*(alpha_top + alpha_bot)
                    alpha = dec*alpha

                elseif (DJfk < c2*DJf0) then 

                    ! Increase alpha
                    alpha_bot = alpha 
                    if (itls < maxit) then 
                        if (alpha_top == inf) then 
                            alpha = inc*alpha_bot 
                        else 
                            alpha = dec*(alpha_top + alpha_bot) 
                        end if 
                    end if 
                    
                else 

                    ! Sufficient decrease
                    conv = .true. 
                    
                end if 
                
                ! Update counter
                itls = itls + 1

                ! De-update the design
                call problem%UpdateDesign(x0-x)

            end do
            
            ! Checks
            if (itls-1 == maxit) then 
                ! Print message, set flag
                if (numLS%verbosity > 0) then 
                    print *, 'linesearch did not converge'
                end if 
            end if 


        case ('backtracking_soc')

            ! Don't compute any derivatives
            doderiv = .false.

            ! Loop
            do while ( (.not. conv) .and. (itls <= maxit) )
            
                ! Update current iterate
                x = x0 + alpha*dphi

                ! Start tracking for possible problems
                call ErrorStack%StartTrack()

                ! Update the design
                call problem%UpdateDesign(alpha*dphi)

                ! Update the problem
                call problem%UpdateProblem()

                ! Check error status
                errstat = ErrorStack%ErrorState()
                call ErrorStack%EndTrack()
                if (errstat > 0) then 
                    ! Error encountered, set value to infinity - don't
                    ! bother trying to compute the merit function
                    fk = inf
                else
                    ! Calculate new cost function value
                    call problem%EvaluateMeritFunction(fk, DJfk, dx, lambda, &
                        mu, doderiv, meritfunction, numLS)
                end if
                
                ! Check Armijo condition
                if (fk < f0 + c1*alpha*DJf0) then 
                    
                    ! Sufficient decrease, terminate
                    conv = .true.
                    
                elseif (errstat > 0) then 

                    ! Error when updating problem, don't even bother 
                    ! trying a second order correction

                    ! Decrease alpha
                    alpha = dec*alpha

                else
                    
                    ! Try if we can get there by applying a second order
                    ! correction

                    ! Note: in Nocedal, this is only done if alpha = 1, but
                    ! this seems to work better if we do it at each attempt. 

                    ! Compute constraints & linearization
                    call problem%EvaluateEqualityConstraints(G, gradG, &
                        hessG, dogradient, dohessian, lambda)
                    call problem%EvaluateInequalityConstraints(H, gradH, &
                        hessH, dogradient, dohessian, mu)

                    ! Compute active set
                    A = H > 0
                    na = count(A)

                    ! Construct problem 
                    allocate(ck(neq + na))
                    ck = [G, pack(H, A)]
                    Ak = gradG%Concatenate(gradH%DeleteColumns(.not. A), 2)
                    LSA = Ak%Transpose()*Ak

                    ! Check if constraints are bounded, otherwise skip
                    if (all(ieee_is_finite(ck))) then
                        ! Compute correction step
                        call SolveSparseLinearSystemDI(LSA, ck, wkt, flag2)

                        ! Check if it converged, otherwise skip update
                        if (flag2 /= 0) then 
                            print *, 'linesearch backtracking_soc: could not converge problem, skipping soc update'
                            wkt(:) = 0
                        end if 

                        ! Compute correction
                        wk = MatrixVectorProduct(Ak, -wkt)

                        ! Recompute cost function at step x0 + alpha*d + wk
                        ! Note: the problem is already updated to x + alpha*d!
                        x = x0 + alpha*dphi + wk ! update x to ensure proper downdate later
                        call problem%UpdateDesign(wk)
                        call problem%UpdateProblem()
                        call problem%EvaluateMeritFunction(fk, DJfk, dx, &
                            lambda, mu, doderiv, meritfunction, numLS)
                    end if 

                    ! Housekeeping
                    deallocate(ck)

                    ! Check the Armijo condition again
                    if (fk < f0 + c1*alpha*DJf0) then 

                        ! Sufficient decrease, terminate
                        conv = .true.

                    else

                        ! Decrease alpha
                        alpha = dec*alpha

                    end if 
                    
                end if 
                
                ! Update counter
                itls = itls + 1

                ! De-update the design
                call problem%UpdateDesign(x0-x)

            end do
            
            ! Checks
            if (itls-1 == maxit) then 
                ! Print message, set flag
                if (numLS%verbosity > 0) then 
                    print *, 'linesearch did not converge'
                end if 
            end if 

        case default

            ! Unknown option, throw error
            call gdErrorHandler('ComputeStepLength: unknown line search option: ' &
                 // numLS%type)

        end select 

        ! Apply step length to dphi
        dx(1:nphi) = dx(1:nphi)*alpha

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Correction equation setup
    subroutine SetupCorrectionEquation(solver, lhs, rhs, &
        gradJ, hessJ, G, gradG, hessG, lambda, H, gradH, hessH, mu, A, &
        ncp, gradncpphi, gradncpmu)

        ! Description
        !============
        ! This routine sets up the correction equation where lhs is a
        ! sparse matrix and rhs the residual vector. The system of 
        ! equations is assumed to be the KKT system where the 
        ! inequality constraints are replaced by the ncp function. 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationSolverKKTUDT)     :: solver 

        real(R8), intent(inout)             :: rhs(:) 
        type(MySparseUDT), intent(inout)    :: lhs 

        real(R8), intent(in)                :: G(:), H(:), mu(:), &
                                            lambda(:), gradJ(:), ncp(:)
        type(MySparseUDT), intent(in)       :: hessJ, gradG, hessG, &
                                            gradH, hessH, gradncpphi, &
                                            gradncpmu
        logical, intent(in)                 :: A(:)

        ! Loop variables
        integer(I8)                         :: k 

        ! Auxiliary variables
        real(R8), allocatable               :: tmu(:), val(:)
        integer(I8), allocatable            :: row(:), col(:)

        ! Set up rhs
        !===========
        if (.not. allocated(tmu)) then 
            allocate(tmu, source=mu)
        end if
        tmu = mu
        where (.not. A) tmu = 0
        rhs = -[gradJ + gradG%MatrixVectorProduct(lambda) + &
            gradH%MatrixVectorProduct(tmu), G, ncp]

        ! Set up lhs
        !===========
        ! Dimensions
        lhs%nrow = size(rhs, 1)
        lhs%ncol = size(rhs, 1)
        lhs%nval = hessJ%nval + hessG%nval + hessH%nval + 2*gradG%nval &
            + 2*gradncpphi%nval + gradncpmu%nval
        allocate(row(lhs%nval), col(lhs%nval), val(lhs%nval))

        ! dLdpsi2
        !--------
        ! Initialize
        k = 0

        ! Cost function contribution
        val(k+1:k+hessJ%nval) = hessJ%val
        row(k+1:k+hessJ%nval) = hessJ%row 
        col(k+1:k+hessJ%nval) = hessJ%col 
        k = k + hessJ%nval 

        ! Equality constraints contribution
        val(k+1:k+hessG%nval) = hessG%val
        row(k+1:k+hessG%nval) = hessG%row 
        col(k+1:k+hessG%nval) = hessG%col 
        k = k + hessG%nval 

        ! Inequality constraints contribution
        val(k+1:k+hessH%nval) = hessH%val
        row(k+1:k+hessH%nval) = hessH%row 
        col(k+1:k+hessH%nval) = hessH%col 
        k = k + hessH%nval 

        ! dLdlambdadphi
        !--------------
        val(k+1:k+gradG%nval) = gradG%val 
        row(k+1:k+gradG%nval) = gradG%row 
        col(k+1:k+gradG%nval) = gradG%col + hessJ%nrow
        k = k + gradG%nval 

        ! dLdphidlambda 
        !--------------
        val(k+1:k+gradG%nval) = gradG%val 
        row(k+1:k+gradG%nval) = gradG%col + hessJ%nrow
        col(k+1:k+gradG%nval) = gradG%row  
        k = k + gradG%nval 

        ! dLdmudphi
        !----------
        val(k+1:k+gradncpphi%nval) = gradncpphi%val 
        row(k+1:k+gradncpphi%nval) = gradncpphi%row
        col(k+1:k+gradncpphi%nval) = gradncpphi%col + hessJ%nrow + gradG%ncol
        k = k + gradncpphi%nval 

        ! dLdphidmu
        !----------
        val(k+1:k+gradncpphi%nval) = gradncpphi%val 
        row(k+1:k+gradncpphi%nval) = gradncpphi%col + hessJ%nrow + gradG%ncol
        col(k+1:k+gradncpphi%nval) = gradncpphi%row
        k = k + gradncpphi%nval 

        ! dLdmudmu
        !----------
        val(k+1:k+gradncpmu%nval) = gradncpmu%val 
        row(k+1:k+gradncpmu%nval) = gradncpmu%col + hessJ%nrow + gradG%ncol
        col(k+1:k+gradncpmu%nval) = gradncpmu%row + hessJ%nrow + gradG%ncol
        k = k + gradncpmu%nval 

        ! Construct sparse
        !=================
        lhs = ConstructMySparse(row, col, val, lhs%nrow, lhs%ncol)

    end subroutine

    !------------------------------------------------------------------!
    !                      QN OPTIMIZATION SOLVER                      !
    !------------------------------------------------------------------!

    ! QN solver initialization
    subroutine InitializeQNSolver(solver)

        ! Description
        !============
        ! Initialize the necessary fields and structures for the QN 
        ! solver. This is basically only the numerics - the optimization
        ! problem itself should already be initialized beforehand!

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationSolverQNUDT)             :: solver    

        ! Data

        ! Initialize
        !===========
        ! Numerics
        solver%num%inputfilepath    = solver%inputfilepath
        solver%num%fieldprefix      = solver%inputfileprefix
        solver%numLS%inputfilepath  = solver%inputfilepath
        solver%numLS%fieldprefix    = solver%inputfileprefix
        solver%numQN%inputfilepath  = solver%inputfilepath
        solver%numQN%fieldprefix    = solver%inputfileprefix
        call solver%num%InitializeNumParams()
        call solver%numLS%InitializeNumParams()
        call solver%numQN%InitializeNumParams()

    end subroutine

    ! QN solver
    subroutine SolveOptimizationProblemQN(solver, problem)

        ! Description
        !============
        ! QN solver for the optimization problem defined by the generic
        ! 'problem'. It is assumed that the optimization problem is
        ! properly initialized. The solver solves the optimization 
        ! problem by computing the solution to quadratic subproblems 
        ! where the hessian is estimated during the optimization. A 
        ! linesearch procedure is used to determine the step size. 

        ! Note that, due to the implementation of the optimization 
        ! problem structure, a hessian output argument is expected from
        ! cost function and constraint routines. However, this may be 
        ! left unitialized/unused, as it will not be used by the solver.

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationSolverQNUDT)              :: solver    
        class(OptimizationProblemUDT)               :: problem

    end subroutine

    ! QN convergence checker
    subroutine CheckConvergenceQN(solver, gradient, converged, infnorm)

        ! Description
        !============
        ! Convergence checker for QN solver. Convergence is based on 
        ! the infinity norm of the provided gradient (which is assumed
        ! to be the gradient of the Lagrangian)

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationSolverQNUDT)      :: solver 
        real(R8), intent(in)                :: gradient(:)
        logical, intent(inout)              :: converged 
        real(R8), intent(out)               :: infnorm

        ! Check convergence
        !==================
        ! Compute infinity norm
        infnorm = maxval(abs(gradient))

        ! Compare
        converged = infnorm < solver%num%tol

    end subroutine

 

    !==================================================================!
    !                                                                  !
    !                            DIAGNOSTICS                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                           COST FUNCTION                          !
    !------------------------------------------------------------------!
    ! Stand-alone driver to check the cost function
    subroutine CheckCostFunctionLinearization(problem, FDchecker, checkgradient, &
        checkhessian)

        ! Description
        !============
        ! This routine compares the gradient computed by finite
        ! differences with the actual gradient computation by 
        ! using the mod_diagnostics module. The errors are printed to
        ! the terminal. 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemUDT)       :: problem 
        logical, intent(in)                 :: checkgradient, &
            checkhessian

        ! Auxiliary
        type(FDcheckerUDT)                  :: FDchecker 

        integer(I8)                         :: nphi, neq, nineq 

        ! Initialize
        !===========
        ! Get the problem dimensions
        call problem%GetProblemDimensions(nphi, neq, nineq)

        ! Sanity checks
        if (any(FDchecker%vars > nphi)) then
            ! Throw error
            call gdErrorHandler('Design indices exceed the number &
                & of design variables!')
        end if

        ! Initialize checker 
        if (.not. allocated(FDchecker%fun)) then 
            allocate(DFCostfunctionUDT::FDchecker%fun)
        end if

        ! Associate
        associate(&
            fun         => FDchecker%fun)
        
        ! Initialize checker functino
        select type(fun)

        type is (DFCostfunctionUDT)

            ! Set the problem
            fun%problem = problem

        class default

            ! Throw error
            call gdErrorHandler('CheckCostFunctionLinearization: unexpected cost function type')

        end select

        ! End associate
        end associate

        ! Compute errors
        !===============
        if (checkgradient) then 
            call FDchecker%CheckGradient()
        end if
        if (checkhessian) then 
            call FDchecker%CheckHessian()
        end if


    end subroutine

    ! Evaluation
    subroutine EvaluateDFCostfunction(fun, x, f, df, d2f, dogradient, &
        dohessian)

        ! Description
        !============
        ! Evaluate the cost function by first updating the design and
        ! then computing the cost function value. 

        ! Declare variables
        !==================
        ! Arguments
        class(DFCostfunctionUDT)        :: fun 
        real(R8), allocatable           :: x(:), df(:)
        real(R8)                        :: f
        type(MySparseUDT)               :: d2f 
        logical                         :: dogradient, dohessian

        ! Auxiliary
        integer(I8)                     :: nx, neq, nineq
        real(R8), allocatable           :: xref(:)

        ! Update design
        !==============
        ! Get dimensions
        call fun%problem%GetProblemDimensions(nx, neq, nineq)

        ! Allocate
        allocate(xref(nx))

        ! Get current design variables
        call fun%problem%GetProblemDesignVariables(xref)

        ! Update with dx 
        call fun%problem%UpdateDesign(x - xref)

        ! Update the optimization problem 
        call fun%problem%UpdateProblem()

        ! Evaluate cost function
        call fun%problem%EvaluateCostFunction(f, df, & 
            d2f, dogradient, dohessian)

    end subroutine

    ! Dimension getter
    subroutine GetProblemDimensions(fun, dimx)

        ! Description
        !============
        ! Get cost function dimensions

        ! Declare variables
        !==================
        ! Arguments
        class(DFCostfunctionUDT)        :: fun 
        integer(I8)                     :: dimx, neq, nineq

        ! Get dimensions
        !===============
        call fun%problem%GetProblemDimensions(dimx, neq, nineq)

    end subroutine

    ! Arguments getter
    subroutine GetProblemArguments(fun, x)

        ! Description
        !============
        ! Get cost function arguments

        ! Declare variables
        !==================
        ! Arguments
        class(DFCostfunctionUDT)        :: fun 
        real(R8), allocatable           :: x(:)

        ! Get dimensions
        !===============
        call fun%problem%GetProblemDesignVariables(x)

    end subroutine

    !------------------------------------------------------------------!
    !                           LAGRANGIAN                             !
    !------------------------------------------------------------------!
    ! Stand-alone driver 
    subroutine CheckLagrangianLinearization(problem, solver, lambda, mu, &
        checkgradient, checkhessian)

        ! Description
        !============
        ! This routine compares the gradient computed by finite
        ! differences with the actual gradient computation by 
        ! using the mod_diagnostics module. The errors are printed to
        ! the terminal. 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemUDT)       :: problem 
        class(OptimizationSolverUDT)        :: solver
        real(R8), allocatable               :: lambda(:), mu(:)
        logical                             :: checkgradient, &
            checkhessian

        ! Auxiliary
        type(FDcheckerUDT)                  :: FDchecker 

        integer(I8)                         :: nvars, nphi, neq, nineq 
        integer(I8), allocatable            :: vars(:)

        ! Initialize
        !===========
        ! Get the problem dimensions
        call problem%GetProblemDimensions(nphi, neq, nineq)

        ! Set the design variables to check
        nvars = 5
        allocate(vars(nvars))
        vars = [1, 1+860, 3, 3+860, 5] ! some random variables for now

        ! Sanity checks
        if (any(vars > nphi)) then
            ! Throw error
            call gdErrorHandler('Design indices exceed the number &
                & of design variables!')
        end if

        ! Initialize checker 
        call FDchecker%Initialize(vars, real([1e-2, 1e-4, 1e-6, 1e-8], kind=R8), 'fd_check_lagrange')
        allocate(DFLagrangianUDT::FDchecker%fun)

        ! Associate
        associate(&
            fun         => FDchecker%fun)
        
        ! Initialize checker functino
        select type(fun)

        type is (DFLagrangianUDT)

            ! Set the problem
            fun%problem = problem
            fun%solver  = solver
            allocate(fun%lambda(size(lambda)))
            allocate(fun%mu(size(mu)))
            fun%lambda = lambda 
            fun%mu = mu

        end select

        ! End associate
        end associate

        ! Compute errors
        !===============
        if (checkgradient) then 
            call FDchecker%CheckGradient()
        end if
        if (checkhessian) then 
            call FDchecker%CheckHessian()
        end if

        ! Deallocate
        deallocate(vars)


    end subroutine

    ! Evaluation
    subroutine EvaluateDFLagrangian(fun, x, f, df, d2f, dogradient, &
        dohessian)

        ! Description
        !============
        ! Evaluate the cost function by first updating the design and
        ! then computing the cost function value. 

        ! Declare variables
        !==================
        ! Arguments
        class(DFLagrangianUDT)          :: fun 
        real(R8), allocatable           :: x(:), df(:)
        real(R8)                        :: f
        type(MySparseUDT)               :: d2f 
        logical                         :: dogradient, dohessian

        ! Auxiliary
        integer(I8)                     :: nphi, neq, nineq
        real(R8), allocatable           :: phiref(:)

        ! Cost function 
        real(R8)                    :: J 
        real(R8), allocatable       :: gradJ(:)
        type(MySparseUDT)           :: hessJ 

        ! Equality constraints 
        real(R8), allocatable       :: G(:), lambda(:) 
        type(MySparseUDT)           :: gradG, hessG  

        ! Inequality constraints 
        real(R8), allocatable       :: H(:), mu(:) 
        type(MySparseUDT)           :: gradH, hessH  

        ! nonlinear complementarity function
        real(R8), allocatable       :: ncp(:) 
        type(MySparseUDT)           :: gradncp
        logical, allocatable        :: A(:), I(:) 

        ! Lagrangian 
        real(R8), allocatable       :: gradL(:)
        type(MySparseUDT)           :: hessL 

        ! Initialize
        !===========
        ! Dimensions
        call fun%problem%GetProblemDimensions(nphi, neq, nineq)

        ! Cost function 
        allocate(gradJ(nphi))
        J = 0
        gradJ(:) = 0
        hessJ%nrow = nphi 
        hessJ%ncol = nphi 

        ! Equality constraint 
        allocate(G(neq), lambda(neq))
        G(:) = 0
        lambda(:) = fun%lambda
        gradG%nrow = nphi 
        gradG%ncol = neq 
        hessG%nrow = nphi 
        hessG%ncol = nphi

        ! Inequality constraint 
        allocate(H(nineq), mu(nineq))
        H(:) = 0
        mu(:) = fun%mu
        gradH%nrow = nphi 
        gradH%ncol = nineq 
        hessH%nrow = nphi 
        hessH%ncol = nphi

        ! NCP function
        allocate(ncp(nineq), A(nineq), I(nineq))
        gradncp%nrow = nphi 
        gradncp%ncol = nineq
        A(:) = .false.
        I(:) = .not. A

        ! Lagrangian 
        allocate(gradL(nphi + neq + nineq))
        hessL%nrow = nphi + neq + nineq
        hessL%ncol = nphi + neq + nineq


        ! Update design
        !==============
        ! Allocate
        allocate(phiref(nphi))

        ! Get current design variables
        call fun%problem%GetProblemDesignVariables(phiref)

        ! Update with dx - x contains lambda and mu as well!
        call fun%problem%UpdateDesign(x(1:nphi) - phiref)
        fun%lambda  = x(nphi+1:nphi+neq)
        fun%mu      = x(nphi+neq+1:nphi+neq+nineq)

        ! Update the optimization problem 
        call fun%problem%UpdateProblem()

        ! Evaluate cost function
        call fun%problem%EvaluateCostFunction(J, gradJ, & 
            hessJ, dogradient, dohessian)

        ! Evaluate the equality constraints
        call fun%problem%EvaluateEqualityConstraints(G, gradG, &
            hessG, dogradient, dohessian, fun%lambda)

        ! Evaluate the inequality constraints
        call fun%problem%EvaluateInequalityConstraints(H, gradH, &
            hessH, dogradient, dohessian, fun%mu)

        ! Evaluate the nonlinear complementarity function 
        !call solver%EvaluateNCPfunction(ncp, A, I, gradncp, &
        !    H, gradH, mu)

        ! Evaluate the Lagrangian
        call fun%problem%EvaluateLagrangian(f, df, d2f, &
            J, gradJ, hessJ, &
            G, gradG, hessG, fun%lambda, &
            H, gradH, hessH, fun%mu, A, &
            dogradient, dohessian)

    end subroutine

    ! Dimension getter
    subroutine GetProblemDimensionsLagrangian(fun, dimx)

        ! Description
        !============
        ! Get lagrangian dimensions

        ! Declare variables
        !==================
        ! Arguments
        class(DFLagrangianUDT)          :: fun 
        integer(I8)                     :: dimx, neq, nineq, nphi

        ! Get dimensions
        !===============
        call fun%problem%GetProblemDimensions(nphi, neq, nineq)
        dimx = nphi + neq + nineq

    end subroutine

    ! Arguments getter
    subroutine GetProblemArgumentsLagrangian(fun, x)

        ! Description
        !============
        ! Get lagrangian arguments

        ! Declare variables
        !==================
        ! Arguments
        class(DFLagrangianUDT)          :: fun 
        real(R8), allocatable           :: x(:)

        ! Auxiliary
        real(R8), allocatable           :: phi(:)
        integer(I8)                     :: nphi, neq, nineq

        ! Get dimensions
        call fun%problem%GetProblemDimensions(nphi, neq, nineq)

        ! Allocate
        allocate(phi(nphi))

        ! Get design variables
        call fun%problem%GetProblemDesignVariables(phi)

        ! Set arguments
        x(1:nphi)                       = phi
        x(nphi+1:nphi+neq)              = fun%lambda 
        x(nphi+neq+1:nphi+neq+nineq)    = fun%mu

    end subroutine

    !------------------------------------------------------------------!
    !                       EQUALITY CONSTRAINTS                       !
    !------------------------------------------------------------------!
    ! Stand-alone driver 
    subroutine CheckEqconLinearization(problem, FDchecker, lambda, checkgradient, &
        checkhessian, eqID)

        ! Description
        !============
        ! This routine compares the gradient computed by finite
        ! differences with the actual gradient computation by 
        ! using the mod_diagnostics module. The errors are printed to
        ! the terminal. Here, since there are often multiple,
        ! constraints, we also need to specify which constraints to 
        ! check and to loop over. Note that, since we only have a 
        ! general routine to evaluate the constraints as abstraction is 
        ! made, the computational cost of doing this scales with the 
        ! computational cost of computing *all* the constraints, i.e. 
        ! we compute the values for all constraints and afterwards 
        ! extract the ones we need. 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemUDT)       :: problem 
        logical                             :: checkgradient, &
            checkhessian

        ! Auxiliary
        type(FDcheckerUDT)                  :: FDchecker 

        integer(I8)                         :: neq, neqID, nphi, nineq  
        integer(I8), intent(in)             :: eqID(:)

        real(R8), allocatable               :: lambda(:)

        ! Loop
        integer(I8)                         :: i
        

        ! Initialize
        !===========
        ! Get the problem dimensions
        call problem%GetProblemDimensions(nphi, neq, nineq)
        neqID =  size(eqID)


        ! Sanity checks
        if (any(eqID > neq)) then
            ! Throw error
            call gdErrorHandler('Constraint indices exceed the number &
                & of constraints!')
        end if
        if (any(FDchecker%vars > nphi)) then
            ! Throw error
            call gdErrorHandler('Design indices exceed the number &
                & of design variables!')
        end if


        ! Initialize checker 
        if (.not. allocated(FDchecker%fun)) then 
            allocate(DFEqconUDT::FDchecker%fun)
        end if 

        ! Associate
        associate(&
            fun         => FDchecker%fun)
        
        ! Initialize checker functino
        select type(fun)

        type is (DFEqconUDT)

            ! Set the problem
            fun%problem = problem
            allocate(fun%lambda(neq))
            fun%lambda = lambda

            ! Compute errors
            !===============
            ! Loop over all constraints
            do i = 1, neqID
                ! Set the correct ID
                fun%eqID = eqID(i)

                ! Print
                print *, '============================================'
                print *, 'Checking equality constraint number: ', eqID(i)
                print *, '============================================'

                ! Check
                if (checkgradient) then 
                    ! Don't do lambda multiplication
                    fun%multlambda = .false.
                    call FDchecker%CheckGradient()
                end if
                if (checkhessian) then 
                    ! Do lambda multiplication
                    fun%multlambda = .true.
                    call FDchecker%CheckHessian()
                end if
            end do

        class default 

            ! Throw error
            call gdErrorHandler('CheckEqonLinearization: unexpected function type')

        end select

        ! End associate
        end associate

    end subroutine

    ! Evaluation
    subroutine EvaluateDFEqcon(fun, x, f, df, d2f, dogradient, &
        dohessian)

        ! Description
        !============
        ! Evaluate the equality cons by first updating the design and
        ! then computing the constraint values. Afterwards, the correct
        ! constraint is unpacked. 

        ! Notes
        !======
        ! When checking the hessian-vector product, one must make sure 
        ! that only the contributions of the current constraint that is
        ! to be checked are taken into account. To this end, all lambda
        ! values, except the one of the current constraint, are set to
        ! zero. 

        ! Declare variables
        !==================
        ! Arguments
        class(DFEQconUDT)               :: fun 
        real(R8), allocatable           :: x(:), df(:)
        real(R8)                        :: f
        type(MySparseUDT)               :: d2f 
        logical                         :: dogradient, dohessian

        ! Auxiliary
        integer(I8)                     :: nphi, neq, nineq
        real(R8), allocatable           :: phiref(:)

        ! Equality constraints 
        real(R8), allocatable       :: G(:), lambda(:)
        type(MySparseUDT)           :: gradG, hessG  
        

        ! Initialize
        !===========
        ! Dimensions
        call fun%problem%GetProblemDimensions(nphi, neq, nineq)

        ! Equality constraint 
        allocate(G(neq), lambda(neq))
        G(:) = 0
        lambda(:) = 0 ! Important! 
        lambda(fun%eqID) = fun%lambda(fun%eqID) ! only keep the current lambda
        gradG%nrow = nphi 
        gradG%ncol = neq 
        hessG%nrow = nphi 
        hessG%ncol = nphi

        ! Update design
        !==============
        ! Allocate
        allocate(phiref(nphi))

        ! Get current design variables
        call fun%problem%GetProblemDesignVariables(phiref)

        ! Update with dx 
        call fun%problem%UpdateDesign(x - phiref)

        ! Update the optimization problem 
        call fun%problem%UpdateProblem()

        ! Evaluate the equality constraints
        call fun%problem%EvaluateEqualityConstraints(G, gradG, &
            d2f, dogradient, dohessian, lambda)

        ! Extract the correct constraint
        f = G(fun%eqID)
        if (dogradient) then
            call gradG%ExtractColumnFull(df, fun%eqID)

            ! Check if we need to multiply with lambda - e.g. when 
            ! computing FD for hessian
            if (fun%multlambda) then
                df = df*lambda(fun%eqID)
            end if

        end if

    end subroutine

    ! Dimension getter
    subroutine GetProblemDimensionsEqcon(fun, dimx)

        ! Description
        !============
        ! Get lagrangian dimensions

        ! Declare variables
        !==================
        ! Arguments
        class(DFEQconUDT)               :: fun 
        integer(I8)                     :: dimx, neq, nineq, nphi

        ! Get dimensions
        !===============
        call fun%problem%GetProblemDimensions(nphi, neq, nineq)
        dimx = nphi

    end subroutine

    ! Arguments getter
    subroutine GetProblemArgumentsEqcon(fun, x)

        ! Description
        !============
        ! Get lagrangian arguments

        ! Declare variables
        !==================
        ! Arguments
        class(DFEQconUDT)               :: fun 
        real(R8), allocatable           :: x(:)

        ! Auxiliary
        real(R8), allocatable           :: phi(:)
        integer(I8)                     :: nphi, neq, nineq

        ! Get dimensions
        call fun%problem%GetProblemDimensions(nphi, neq, nineq)

        ! Allocate
        allocate(phi(nphi))

        ! Get design variables
        call fun%problem%GetProblemDesignVariables(phi)

        ! Set arguments
        x = phi

    end subroutine

    !------------------------------------------------------------------!
    !                      INEQUALITY CONSTRAINTS                      !
    !------------------------------------------------------------------!
    ! Stand-alone driver 
    subroutine CheckIneqconLinearization(problem, FDchecker, mu, checkgradient, &
        checkhessian, ineqID)

        ! Description
        !============
        ! This routine compares the gradient computed by finite
        ! differences with the actual gradient computation by 
        ! using the mod_diagnostics module. The errors are printed to
        ! the terminal. Here, since there are often multiple,
        ! constraints, we also need to specify which constraints to 
        ! check and to loop over. Note that, since we only have a 
        ! general routine to evaluate the constraints as abstraction is 
        ! made, the computational cost of doing this scales with the 
        ! computational cost of computing *all* the constraints, i.e. 
        ! we compute the values for all constraints and afterwards 
        ! extract the ones we need. 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemUDT)       :: problem 
        logical                             :: checkgradient, &
            checkhessian

        ! Auxiliary
        type(FDcheckerUDT)                  :: FDchecker 

        integer(I8)                         :: nineq, nineqID, nphi, neq  
        integer(I8), intent(in)             :: ineqID(:)

        real(R8), allocatable               :: mu(:)

        ! Loop
        integer(I8)                         :: i
        

        ! Initialize
        !===========
        ! Get the problem dimensions
        call problem%GetProblemDimensions(nphi, neq, nineq)
        nineqID =  size(ineqID)

        ! Sanity checks
        if (any(ineqID > nineq)) then
            ! Throw error
            call gdErrorHandler('Constraint indices exceed the number &
                & of constraints!')
        end if
        if (any(FDchecker%vars > nphi)) then
            ! Throw error
            call gdErrorHandler('Design indices exceed the number &
                & of design variables!')
        end if


        ! Initialize checker 
        if (.not. allocated(FDchecker%fun)) then
            allocate(DFIneqconUDT::FDchecker%fun)
        end if 

        ! Associate
        associate(&
            fun         => FDchecker%fun)
        
        ! Initialize checker functino
        select type(fun)

        type is (DFIneqconUDT)

            ! Set the problem
            fun%problem = problem
            !allocate(fun%mu(nineq))
            fun%mu = mu

            ! Compute errors
            !===============
            ! Loop over all constraints
            do i = 1, nineqID
                ! Set the correct ID
                fun%ineqID = ineqID(i)

                ! Print
                print *, '============================================'
                print *, 'Checking inequality constraint number: ', ineqID(i)
                print *, '============================================'

                ! Check
                if (checkgradient) then 
                    ! Don't do mu multiplication
                    fun%multmu = .false.
                    call FDchecker%CheckGradient()
                end if
                if (checkhessian) then 
                    ! Do mu multiplication
                    fun%multmu = .true.
                    call FDchecker%CheckHessian()
                end if
            end do

        class default 

            ! Throw error
            call gdErrorHandler('CheckIneqconLinearization: unexpected function type')

        end select

        ! End associate
        end associate

    end subroutine

    ! Evaluation
    subroutine EvaluateDFIneqcon(fun, x, f, df, d2f, dogradient, &
        dohessian)

        ! Description
        !============
        ! Evaluate the equality cons by first updating the design and
        ! then computing the constraint values. Afterwards, the correct
        ! constraint is unpacked. 

        ! Notes
        !======
        ! When checking the hessian-vector product, one must make sure 
        ! that only the contributions of the current constraint that is
        ! to be checked are taken into account. To this end, all lambda
        ! values, except the one of the current constraint, are set to
        ! zero. 

        ! Declare variables
        !==================
        ! Arguments
        class(DFIneqconUDT)             :: fun 
        real(R8), allocatable           :: x(:), df(:)
        real(R8)                        :: f
        type(MySparseUDT)               :: d2f 
        logical                         :: dogradient, dohessian

        ! Auxiliary
        integer(I8)                     :: nphi, neq, nineq
        real(R8), allocatable           :: phiref(:)

        ! Equality constraints 
        real(R8), allocatable       :: H(:), mu(:)
        type(MySparseUDT)           :: gradH, hessH  
        

        ! Initialize
        !===========
        ! Dimensions
        call fun%problem%GetProblemDimensions(nphi, neq, nineq)

        ! Equality constraint 
        allocate(H(nineq), mu(nineq))
        H(:) = 0
        mu(:) = 0 ! Important! 
        mu(fun%ineqID) = fun%mu(fun%ineqID) ! only keep the current mu
        gradH%nrow = nphi 
        gradH%ncol = neq 
        hessH%nrow = nphi 
        hessH%ncol = nphi

        ! Update design
        !==============
        ! Allocate
        allocate(phiref(nphi))

        ! Get current design variables
        call fun%problem%GetProblemDesignVariables(phiref)

        ! Update with dx 
        call fun%problem%UpdateDesign(x - phiref)

        ! Update the optimization problem 
        call fun%problem%UpdateProblem()

        ! Evaluate the equality constraints
        call fun%problem%EvaluateInequalityConstraints(H, gradH, &
            d2f, dogradient, dohessian, mu)

        ! Extract the correct constraint
        f = H(fun%ineqID)
        if (dogradient) then
            call gradH%ExtractColumnFull(df, fun%ineqID)

            ! Check if we need to multiply with lambda - e.g. when 
            ! computing FD for hessian
            if (fun%multmu) then
                df = df*mu(fun%ineqID)
            end if

        end if

    end subroutine

    ! Dimension getter
    subroutine GetProblemDimensionsIneqcon(fun, dimx)

        ! Description
        !============
        ! Get lagrangian dimensions

        ! Declare variables
        !==================
        ! Arguments
        class(DFIneqconUDT)             :: fun 
        integer(I8)                     :: dimx, neq, nineq, nphi

        ! Get dimensions
        !===============
        call fun%problem%GetProblemDimensions(nphi, neq, nineq)
        dimx = nphi

    end subroutine

    ! Arguments getter
    subroutine GetProblemArgumentsIneqcon(fun, x)

        ! Description
        !============
        ! Get lagrangian arguments

        ! Declare variables
        !==================
        ! Arguments
        class(DFIneqconUDT)             :: fun 
        real(R8), allocatable           :: x(:)

        ! Auxiliary
        real(R8), allocatable           :: phi(:)
        integer(I8)                     :: nphi, neq, nineq

        ! Get dimensions
        call fun%problem%GetProblemDimensions(nphi, neq, nineq)

        ! Allocate
        allocate(phi(nphi))

        ! Get design variables
        call fun%problem%GetProblemDesignVariables(phi)

        ! Set arguments
        x = phi

    end subroutine

end module
