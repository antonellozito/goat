!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains the optimization monitor class, of which the 
! objects can keep track of the progress of the optimization and 
! print out any useful information. The monitor is normally stored in 
! the optimization problem structure. The following fields are being
! tracked 

! Fields
!=======
! (most are R8/I8 precision *-by-maxitopt arrays where maxitopt
! is the maximal number of optimization iterations):
!
! - maxitopt:   The maximum number of iterations   
! - itopt:      The current optimization iterate (scalar)
! - J:          Cost function value (1-by-maxitopt)
! - L:          Lagrangian value (1-by-maxitopt)
! - G:          max. absolute value of the equality constraints (1-by-maxitopt)
! - H:          max. value of the inequality contraints (1-by-maxitopt)
! - convnorm:   Convergence criterion value (1-by-maxitopt)
! - opttol:     Optimization tolerance (scalar)
! - nphi:       Number of design variables
! - neq:        Number of equality constraints
! - nineq:      Number of inequality constraints
! - evaltime:   time it takes to evaluate the problem, hessian, etc (scalar)
! - linsolve:   time it takes to solve the linear system (scalar)
! - ittime:     time it takes for a full iterate (scalar)
! - alpha:      line search step length (if computed, 1-by-maxitopt)

!
! Methods
!========
! Most methods are simply for post-process or monitoring during the 
! optimization process. See the implementation details of each routine
! for additional explanation. 
!
! - PrintHeader:        Print out the header of the optimization
! - PrintIterate:       Print out the current iterate with some data


module optmod_monitor
    
    ! Initialize
    !============
    ! Load modules
    use mod_precision

    ! The usual
    implicit none
    save
    public 

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    ! Monitor
    !========
    type OptimizationMonitorUDT

        ! Description
        !============
        ! See comments at beginning of module
        integer(I8)             :: maxitopt, itopt, nphi, neq, nineq
        real(R8)                :: opttol, rxf, maxdphi
        real(R8), allocatable   :: J(:), L(:), G(:), H(:), convnorm(:), &
            alpha(:) 
        real(R8)                :: evaltime, ittime, linsolvetime

    contains

        ! Printing procedures
        procedure :: PrintHeader        => PrintHeader 
        procedure :: PrintIterate       => PrintIterate

        ! Housekeeping
        procedure :: Initialize         => InitializeMonitor 
        procedure :: Allocate           => AllocateMonitor 
        procedure :: Deallocate         => DeallocateMonitor

    end type

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

    ! Print header
    subroutine PrintHeader(monitor)

        ! Description
        !============
        ! Print out the header at the start of the optimization.

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationMonitorUDT)        :: monitor

        character(:), allocatable            :: myformat 

        ! Print
        !======
        allocate(myformat, source='')
        myformat = repeat('4x, a4, 4x,', 9)
        print *, '!===================================================!'
        print *, '!         Solving the grid deformation problem      !'
        print *, '!===================================================!'
        print *, 'Performing ', monitor%maxitopt, ' iterations'
        print '(' // myformat // '4x, a4, 4x, 4x, a4, 4x, 4x, a4, 4x)' &
            ,   'it  ',  'conv', 'dphi', 'L', 'J', 'Gmax', 'Hmax', 'rxf', &
            'step', 'tol', 'itt', 'eval', 'lst'

    end subroutine

    ! Print iterate
    subroutine PrintIterate(monitor)

        ! Description
        !============
        ! Print out an iteration of the optimization

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationMonitorUDT)        :: monitor 

        character(:), allocatable            :: myformat 

        ! Associate for ease
        associate(&
            itopt => monitor%itopt)

        ! Print
        !======
        allocate(myformat, source='')
        myformat = repeat('es12.4, 4x, ', 12)
        print '(i8, 4x, ' // myformat // '4x)', &
            monitor%itopt, monitor%convnorm(itopt), monitor%maxdphi, &
            monitor%L(itopt), monitor%J(itopt), monitor%G(itopt), &
            monitor%H(itopt), monitor%rxf, monitor%alpha(itopt), &
            monitor%opttol, monitor%ittime, monitor%evaltime, monitor%linsolvetime

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Housekeeping
    subroutine InitializeMonitor(monitor, maxitopt, nphi, neq, nineq, & 
        opttol)

        ! Description
        !============
        ! Monitor initialization (including allocation)

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationMonitorUDT)        :: monitor 
        integer(I8)                         :: maxitopt, nphi, neq, nineq
        real(R8)                            :: opttol 

        ! Initialization
        !===============
        monitor%nphi = nphi
        monitor%neq = neq
        monitor%nineq = nineq
        monitor%maxitopt = maxitopt
        monitor%opttol = opttol 

        ! Allocate
        !=========
        call monitor%Allocate()
        
    end subroutine

    subroutine AllocateMonitor(monitor)

        ! Description
        !============
        ! Monitor allocation. It is assumed that the monitor has been
        ! initialized first. 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationMonitorUDT)        :: monitor 

        ! Allocate
        !=========
        ! Cost function & Lagrangian
        allocate(monitor%J(monitor%maxitopt))
        allocate(monitor%L(monitor%maxitopt))

        ! Constraints
        allocate(monitor%G(monitor%maxitopt))
        allocate(monitor%H(monitor%maxitopt))

        ! Convergence
        allocate(monitor%convnorm(monitor%maxitopt))

        ! Line search
        allocate(monitor%alpha(monitor%maxitopt))

    end subroutine

    subroutine DeallocateMonitor(monitor)

        ! Description
        !============
        ! Monitor deallocation

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationMonitorUDT)        :: monitor 

        ! Allocate
        !=========
        ! Cost function & Lagrangian
        deallocate(monitor%J)
        deallocate(monitor%L)

        ! Constraints
        deallocate(monitor%G)
        deallocate(monitor%H)
        
        ! Line search
        deallocate(monitor%alpha)

        ! Convergence
        deallocate(monitor%convnorm)

    end subroutine

end module