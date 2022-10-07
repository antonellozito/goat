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
    use optmod_state
    use optmod_numerics

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
    type OptimizationProblemUDT

        ! Description
        !============
        ! Defines the basic optimization problem: it has a set of 
        ! design variables, constraints, and a cost function. 
        class(DesignVariablesUDT), allocatable      :: designvariables        
        class(CostfunctionUDT), allocatable         :: costfunction 
        class(ConstraintsUDT), allocatable          :: constraints
        class(StateUDT), allocatable                :: state

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
        procedure :: InitializeKKTSolver => InitializeKKTSolverINT
        procedure :: SolveOptimizationProblemKKT => KKTSolver

    end type

    ! Optimization engine
    !====================
    type OptimizationEngineUDT 

        type(OptimizationProblemUDT)                :: problem
        type(OptimizationSolverUDT)                 :: solver

    contains

        procedure :: SetupOptimizationDriver => SetupOptimizationDriverDummy

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

    ! Dummy initialization routine - to be overwritten by the user
    subroutine SetupOptimizationDriverDummy(optimizationdriver) 

        ! Description
        !============
        ! This is only a dummy routine for setting up the optimization 
        ! engine and should be overwritten 

        class(OptimizationEngineUDT) :: optimizationdriver 

        print *, 'please replace this routine with your own setup routines'

    end subroutine

    ! KKT solver initialization
    subroutine InitializeKKTSolverINT(solver)

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
    subroutine KKTSolver(solver, problem)

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
        integer(I8)                 :: itopt, maxit
        
        ! Auxiliary variables 
        real(R8)                    :: rxf, rxfdesign, rxfdec, rxfmin 

        ! Data

        ! Initialize & unpack
        !====================
        ! Initialize the solver
        call solver%InitializeKKTSolver()

        ! Initialize counter(s)
        itopt = 1
        maxit = 1

        ! Unpack numerical options
        associate(maxit => solver%numKKT%maxit)
            print *, maxit
        end associate
        print *, maxit 


    end subroutine

end module