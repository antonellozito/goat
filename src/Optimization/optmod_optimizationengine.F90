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

    end type

    ! Optimization engine
    !====================
    type OptimizationEngineUDT 

        ! Description
        !============
        ! The engine simply contains the problem and solver structures
        ! (and in the future the monitor). 
        class(OptimizationProblemUDT), allocatable  :: problem
        type(OptimizationSolverUDT)                 :: solver

    contains

        ! Main driver to solve a problem
        procedure :: Driver                 => OptimizationEngineDriver

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
        call optimizationengine%problem%Initialize()

        ! Set up the solver - now KKT
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
    
        ! Loop variables
        integer(I8)                 :: itopt, maxit
        
        ! Auxiliary variables 
        real(R8)                    :: rxf, rxfdesign, rxfdec, rxfmin 

        ! Data

        ! Initialize
        !===========
        ! Numerics
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

        ! Data

        ! Temporary variables (to be deleted in the future)
        integer(I8)                 :: nphi, neq, nineq
        real(R8)                    :: opttol

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
        lambda(:) = 0
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

            ! Evaluate cost function
            call problem%EvaluateCostFunction(J, gradJ, & 
                hessJ, dogradient, dohessian)

            ! Evaluate the equality constraints
            call problem%EvaluateEqualityConstraints(G, gradG, &
                hessG, dogradient, dohessian, lambda)

            ! Evaluate the inequality constraints
            call problem%EvaluateInequalityConstraints(H, gradH, &
                hessH, dogradient, dohessian, lambda)

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
            !call solver%CheckConvergenceKKT(gradL, converged, convnorm)

            ! Update the monitor
            problem%monitor%itopt = itopt

            ! Update the monitor again
            problem%monitor%J(itopt)        = J
            problem%monitor%dJ(:,itopt)     = gradJ
            problem%monitor%L(itopt)        = 0
            problem%monitor%dL(:,itopt)     = 0
            problem%monitor%G(:,itopt)      = G
            problem%monitor%H(:,itopt)      = H
            problem%monitor%convnorm(itopt) = convnorm

            ! Print the current iterate
            if (verbosity > 0) then 
                ! Print out the iterate
                call problem%monitor%PrintIterate()

            end if

            ! Update the iteration counter
            itopt = itopt+1

        end do

        ! Housekeeping
        deallocate(G, H, gradJ, lambda, mu)


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
            gradL = gradJ

            ! Equality constraints contribution
            do k = 1, gradG%nval
                gradL(gradG%col(k)) = lambda(gradG%row(k))*gradG%val(k)
            end do 

            ! Inequality constraints contribution
            do k = 1, gradH%nval 
                gradL(gradH%col(k)) = tmu(gradH%row(k))*gradH%val(k)
            end do

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
        


    end subroutine

    ! Nonlinear complementarity function 
    !subroutine EvaluateNCPfunction(ncp, A, I, gradncp, H, gradH, mu)

    !end subroutine



end module