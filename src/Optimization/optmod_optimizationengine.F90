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

        ! Design initialization
        procedure(InitializeINT), deferred :: Initialize
        

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
        procedure :: InitializeKKTSolver => InitializeKKTSolver
        procedure :: SolveOptimizationProblemKKT => SolveOptimizationProblemKKT

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

    end interface

    contains

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

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
        logical                     :: notconverged
        
        ! Auxiliary variables 
        real(R8)                    :: rxf, rxfdesign, rxfdec, rxfmin 

        ! Data

        ! Temporary variables (to be deleted in the future)
        integer(I8)                 :: nphi, neq, nineq
        real(R8)                    :: opttol

        ! Initialize & unpack
        !====================
        ! Initialize the solver
        call solver%InitializeKKTSolver()

        ! Initialize the monitor - only temporary here
        opttol = 1e-8
        nphi = 1
        neq = 1
        nineq = 1
        call problem%monitor%Initialize(solver%numKKT%maxit, nphi, neq,&
            nineq, opttol)


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
        notconverged = .true. 

        ! Print solver header
        if (verbosity > 0) then
            ! Print out the header
            call problem%monitor%PrintHeader()

        end if

        ! Loop
        do while (notconverged .and. (itopt <= maxit))

            ! Update the monitor
            problem%monitor%itopt = itopt

            ! Update the monitor again
            problem%monitor%J(itopt)        = 0
            problem%monitor%dJ(:,itopt)     = 0
            problem%monitor%L(itopt)        = 0
            problem%monitor%dL(:,itopt)     = 0
            problem%monitor%G(:,itopt)      = 0
            problem%monitor%H(:,itopt)      = 0
            problem%monitor%convnorm(itopt) = 0

            ! Print the current iterate
            if (verbosity > 0) then 
                ! Print out the iterate
                call problem%monitor%PrintIterate()

            end if

            ! Update the iteration counter
            itopt = itopt+1

        end do

    end subroutine


end module