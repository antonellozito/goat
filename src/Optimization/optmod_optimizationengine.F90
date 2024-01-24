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
    use mod_sparseinterface
    use mod_diagnostics
    use mod_linearsolverinterface
    use mod_plotter
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
        character(:), allocatable                   :: inputfilepath
        type(OptimizationMonitorUDT)                :: monitor

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

    end type

    ! Optimization solver
    !====================
    type :: OptimizationSolverUDT

        ! Description
        !============
        ! Defines the optimization solver and its numerics. It should 
        ! have a generic 'solve' method which acts on the optimization
        ! problem
        character(:), allocatable                   :: inputfilepath
        type(NumKKTUDT)         :: numKKT

    contains 

        ! Solution procedure using KKT solver
        procedure ::    InitializeKKTSolver     => InitializeKKTSolver
        procedure ::    SolveOptimizationProblemKKT &
                        => SolveOptimizationProblemKKT

        ! Convergence checking
        procedure :: CheckConvergenceKKT 

        ! Lagrangian evaluation
        procedure :: EvaluateLagrangian

        ! Relaxation of KKT system
        procedure :: RelaxKKTSystem

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
        class(OptimizationProblemUDT), allocatable  :: problem
        type(OptimizationSolverUDT)                 :: solver

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
        type(OptimizationSolverUDT)                     :: solver

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
        optimizationengine%solver%inputfilepath = optimizationengine%inputfilepath
        call optimizationengine%solver%InitializeKKTSolver()

        ! Solve
        !======
        ! Solve the optimization problem by calling the KKT solver
        call optimizationengine%solver%SolveOptimizationProblemKKT( &
            optimizationengine%problem)

    end subroutine

    !------------------------------------------------------------------!
    !                       OPTIMIZATION SOLVER                        !
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
        class(OptimizationSolverUDT)                :: solver    

        ! Data

        ! Initialize
        !===========
        ! Numerics
        solver%numKKT%inputfilepath = solver%inputfilepath
        call solver%numKKT%InitializeNumParams() 

    end subroutine

    ! KKT solver
    subroutine SolveOptimizationProblemKKT(solver, problem)

        ! Description
        !============
        ! KKT solver for the optimization problem defined by the generic
        ! 'problem'. It is assumed that the optimization problem is
        ! properly initialized. The solver is initialized here. 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationSolverUDT)                :: solver    
        class(OptimizationProblemUDT)               :: problem
    
        ! Loop variables
        integer(I8)                 :: itopt, maxit, verbosity
        logical                     :: converged
        
        ! Auxiliary variables 
        real(R8)                    :: rxf, rxfdesign, rxfdec, rxfmin, &
            convnorm
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
        type(MySparseUDT)           :: gradncp
        logical, allocatable        :: A(:), I(:) 

        ! Lagrangian 
        real(R8)                    :: L 
        real(R8), allocatable       :: gradL(:)
        type(MySparseUDT)           :: hessL 

        ! Solver & updates
        real(R8), allocatable       :: fullmat(:, :)
        double precision, allocatable :: dx(:)

        ! Diagnostics
        logical                     :: checkgradients, checkhessians 

        ! Data

        ! Temporary variables (to be deleted in the future)
        integer(I8)                 :: nphi, neq, nineq
        real(R8)                    :: opttol

            external dgesv

        ! Initialize & unpack
        !====================
        ! Logicals
        dogradient  = .true. 
        dohessian   = .true. 

        ! Initialize the solver
        call solver%InitializeKKTSolver()

        ! Initialize the monitor - only temporary here
        opttol = 1e-8
        call problem%GetProblemDimensions(nphi, neq, nineq)
        call problem%monitor%Initialize(solver%numKKT%maxit, nphi, neq,&
            nineq, opttol)

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
        gradncp%nrow = nphi 
        gradncp%ncol = nineq
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
        allocate(fullmat(hessL%nrow, hessL%ncol))

        ! Diagnostics
        checkgradients  = .false. ! check gradients in each iteration?
        checkhessians   = .false. ! check hessians in each iteration?

        ! Initialize counter(s)
        itopt = 1
        maxit = solver%numKKT%maxit

        ! Unpack numerical options
        rxf = solver%numKKT%rxf 
        rxfdesign = solver%numKKT%rxfdesign 
        rxfdec = solver%numKKT%rxfdec 
        rxfmin = solver%numKKT%rxfmin 
        verbosity = solver%numKKT%verbosity

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

            ! Update the optimization problem 
            call problem%UpdateProblem()

            ! Check gradients & hessians if needed
            if (checkgradients .or. checkhessians) then 
                call CheckCostFunctionLinearization(problem, &
                    checkgradients, checkhessians)
                call CheckLagrangianLinearization(problem, solver, lambda, & 
                    mu, checkgradients, checkhessians)
                call CheckEqconLinearization(problem, lambda, &
                    checkgradients, checkhessians)
            end if

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
            !call solver%EvaluateNCPfunction(ncp, A, I, gradncp, &
            !    H, gradH, mu)

            ! Evaluate the Lagrangian
            call solver%EvaluateLagrangian(L, gradL, hessL, &
                J, gradJ, hessJ, &
                G, gradG, hessG, lambda, &
                H, gradH, hessH, mu, A, &
                dogradient, dohessian)

            ! Check convergence
            call solver%CheckConvergenceKKT(gradL, converged, convnorm)

            ! Solve 
            if (.not. converged) then 

                ! Relax
                call solver%RelaxKKTSystem(hessL, nphi, neq, nineq)

                !print *, 'hessian size: ', hessL%nrow, hessL%ncol
                !call SpyPlot(hessL%row, hessL%col, hessL%nval, '-p')
                
                ! Call the sparse solver
                call SolveSparseLinearSystemDI(hessL, -gradL, dx)

            else
                ! Don't solve again, already converged. Exit 
                dx(:) = 0
            end if

            ! Update the design
            call problem%UpdateDesign(solver%numKKT%rxfdesign*dx(1:nphi))

            ! Update lambda
            lambda(:) = lambda(:) + dx(nphi+1:nphi+neq)

            ! Update the monitor
            problem%monitor%itopt = itopt

            ! Update the monitor again
            problem%monitor%J(itopt)        = J
            problem%monitor%dJ(:,itopt)     = gradJ
            problem%monitor%dL(:,itopt)     = 0
            problem%monitor%G(:,itopt)      = G
            problem%monitor%H(:,itopt)      = H
            problem%monitor%convnorm(itopt) = convnorm
            problem%monitor%rxf = solver%numKKT%rxf 

            ! Print the current iterate
            if (verbosity > 0) then 
                ! Print out the iterate
                call problem%monitor%PrintIterate()

            end if

            ! Update the iteration counter
            itopt = itopt+1

            ! Update the hessian relaxation factor
            solver%numKKT%rxf = solver%numKKT%rxf*solver%numKKT%rxfdec
            solver%numKKT%rxf = max(solver%numKKT%rxf, solver%numKKT%rxfmin)

        end do

        ! Housekeeping
        deallocate(G, H, gradJ, lambda, mu)


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
        class(OptimizationSolverUDT)        :: solver 
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
        class(OptimizationSolverUDT)        :: solver 
        real(R8), intent(in)                :: gradient(:)
        logical, intent(inout)              :: converged 
        real(R8), intent(out)               :: infnorm

        ! Check convergence
        !==================
        ! Compute infinity norm
        infnorm = maxval(abs(gradient))

        ! Compare
        converged = infnorm < solver%numKKT%tol

    end subroutine

    ! Lagrangian 
    subroutine EvaluateLagrangian(solver, L, gradL, hessL, J, gradJ, &
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

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationSolverUDT)        :: solver 

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
        real(R8), allocatable               :: tmu(:)
                                        
        ! Compute Lagrangian
        !===================
        ! For ease: set mu(.not. A) equal to zero
        allocate(tmu(size(mu)))
        where (.not. A) tmu = 0
        L = J + sum(lambda * G) + sum(tmu * H)

        ! Compute gradient
        !=================
        if (dogradient) then 

            ! Cost function contribution
            gradL(1:size(gradJ)) = gradJ(:) + &
                gradG%MatrixVectorProduct(lambda) + &
                gradH%MatrixVectorProduct(tmu)

            ! Equality constraints contribution
            gradL(size(gradJ)+1:size(G)) = G(:)

            ! Inequality constraints contribution
            gradL(size(gradJ)+size(G)+1:size(gradL)) = H(:)

        end if 

        ! Compute hessian
        !================
        if (dohessian) then 

            ! Initialize & allocate
            k = 0
            hessL%nval  = hessJ%nval + hessG%nval + hessH%nval &
                        + 2*gradG%nval + 2*gradH%nval   

            if (.not. allocated(hessL%val)) then
                call hessL%Allocate()
            end if

            ! dLdpsi2
            !--------
            ! Cost function contribution
            hessL%val(k+1:k+hessJ%nval) = hessJ%val
            hessL%row(k+1:k+hessJ%nval) = hessJ%row 
            hessL%col(k+1:k+hessJ%nval) = hessJ%col 
            k = k + hessJ%nval 

            ! Equality constraints contribution
            hessL%val(k+1:k+hessG%nval) = hessG%val
            hessL%row(k+1:k+hessG%nval) = hessG%row 
            hessL%col(k+1:k+hessG%nval) = hessG%col 
            k = k + hessG%nval 

            ! Inequality constraints contribution
            hessL%val(k+1:k+hessH%nval) = hessH%val
            hessL%row(k+1:k+hessH%nval) = hessH%row 
            hessL%col(k+1:k+hessH%nval) = hessH%col 
            k = k + hessH%nval 

            ! dLdlambdadphi
            !--------------
            hessL%val(k+1:k+gradG%nval) = gradG%val 
            hessL%row(k+1:k+gradG%nval) = gradG%row
            hessL%col(k+1:k+gradG%nval) = gradG%col + hessJ%nrow
            k = k + gradG%nval 

            ! dLdphidlambda 
            !--------------
            hessL%val(k+1:k+gradG%nval) = gradG%val 
            hessL%row(k+1:k+gradG%nval) = gradG%col + hessJ%nrow
            hessL%col(k+1:k+gradG%nval) = gradG%row  
            k = k + gradG%nval 

            ! dLdmudphi
            !----------
            hessL%val(k+1:k+gradH%nval) = gradH%val 
            hessL%row(k+1:k+gradH%nval) = gradH%row
            hessL%col(k+1:k+gradH%nval) = gradH%col + hessJ%nrow + gradG%ncol
            k = k + gradH%nval 

            ! dLdphidmu
            !----------
            hessL%val(k+1:k+gradH%nval) = gradH%val 
            hessL%row(k+1:k+gradH%nval) = gradH%col + hessJ%nrow + gradG%ncol
            hessL%col(k+1:k+gradH%nval) = gradH%row
            k = k + gradH%nval 

        end if

        ! A simple workaround to avoid unused dummy argument warnings 
        ! during compilation
        if (.false.) then 
            print *, solver%inputfilepath
        end if
        
    end subroutine

    ! Nonlinear complementarity function 
    !subroutine EvaluateNCPfunction(ncp, A, I, gradncp, H, gradH, mu)

    !end subroutine

    !==================================================================!
    !                                                                  !
    !                            DIAGNOSTICS                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                           COST FUNCTION                          !
    !------------------------------------------------------------------!
    ! Stand-alone driver to check the cost function
    subroutine CheckCostFunctionLinearization(problem, checkgradient, &
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
        call FDchecker%Initialize(nvars, vars)
        allocate(DFCostfunctionUDT::FDchecker%fun)

        ! Associate
        associate(&
            fun         => FDchecker%fun)
        
        ! Initialize checker functino
        select type(fun)

        type is (DFCostfunctionUDT)

            ! Set the problem
            fun%problem = problem

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
        type(OptimizationSolverUDT)         :: solver
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
        call FDchecker%Initialize(nvars, vars)
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
        call fun%solver%EvaluateLagrangian(f, df, d2f, &
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
    subroutine CheckEqconLinearization(problem, lambda, checkgradient, &
        checkhessian)

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

        integer(I8)                         :: nvars 
        integer(I8), allocatable            :: vars(:)

        integer(I8)                         :: neq, neqID, nphi, nineq  
        integer(I8), allocatable            :: eqID(:)

        real(R8), allocatable               :: lambda(:)

        ! Loop
        integer(I8)                         :: i
        

        ! Initialize
        !===========
        ! Get the problem dimensions
        call problem%GetProblemDimensions(nphi, neq, nineq)

        ! Set the constraints to check
        neqID =  1
        allocate(eqID(neqID))
        eqID = [902] !  some random numbers for now

        ! Set the design variables to check
        nvars = 5
        allocate(vars(nvars))
        vars = [1, 1+860, 860, 3+860, 1720] ! some random variables for now

        ! Sanity checks
        if (any(eqID > neq)) then
            ! Throw error
            call gdErrorHandler('Constraint indices exceed the number &
                & of constraints!')
        end if
        if (any(vars > nphi)) then
            ! Throw error
            call gdErrorHandler('Design indices exceed the number &
                & of design variables!')
        end if


        ! Initialize checker 
        call FDchecker%Initialize(nvars, vars)
        allocate(DFEqconUDT::FDchecker%fun)

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

        end select

        ! End associate
        end associate

        ! Deallocate
        deallocate(vars, eqID)


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

end module