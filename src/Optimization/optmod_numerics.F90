!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains the numerical parameters for the optimization. 

module optmod_numerics
    
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

    ! Abstract types
    !===============
    ! General numerical parameter object
    type :: NumUDT

        ! Description
        !============
        ! Abstract type that contains the general numerical parameters
        ! and routines (if any) for optimization solvers. 
        ! Fields:
        ! - tol:            tolerance to which the solver has to solve
        ! - maxit:          maximum number of iterations
        ! - verbosity:      verbosity level of printing

        real(R8)            :: tol 
        integer(I8)         :: maxit 
        integer(I8)         :: verbosity

    contains 

        procedure :: SetDefaultNumParams => SetDefaultNumParamsINT
        procedure :: InitializeNumParams => InitializeGeneralNumParamsINT

    end type

    ! Derived types
    !==============
    type, extends(NumUDT) :: NumKKTUDT 

        ! Description
        !============
        ! Numerical parameters for the KKT solver. The following fields
        ! are added:
        ! - rxf:    relaxation factor on the cost function hessian 
        !           (sum(abs(hessJ,2))*unity matrix*rxf is added). 
        !           Higher rxf -> more relaxation -> (normally) better 
        !           convergence. Set to zero for no relaxation.
        ! - rxfdesign:  relaxation factor on the design update. Rescales
        !               the design update as 
        !               x_new = x_old + rxfdesign*delta x. 
        !               smaller values -> slower convergence, yet more
        !               stable. Set to one for no relaxation. 
        ! - rxfdec: relaxation reduction factor per iteration, i.e. the
        !           next iteration, rxf_new = rxf_old*rxfdec. Set to one
        !           for constant relaxation factor. Set smaller to one 
        !           to reduce relaxation factor (0.98 works good in most
        !           cases)
        ! - rxfmin: minimal value of relaxation factor

        ! Relaxation factors
        real(R8)            :: rxf
        real(R8)            :: rxfdec 
        real(R8)            :: rxfmin
        real(R8)            :: rxfdesign

    contains 

        procedure :: SetDefaultNumParamsKKT => SetDefaultNumParamsKKTINT
        procedure :: InitializeNumParams => InitializeNumParamsKKTINT

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

    !------------------------------------------------------------------!
    !                           GENERAL OPTIONS                        !
    !------------------------------------------------------------------!
    ! Set the general numerical parameters
    subroutine SetDefaultNumParamsINT(num)

        ! Description
        !============
        ! Set default numerical parameters tol, itmax, verbosity

        ! Declare variables
        !==================
        ! Arguments
        class(NumUDT)            :: num
    
        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Set defaults
        !=============
        num%tol         = 1e-6
        num%maxit       = 10
        num%verbosity   = 1

    end subroutine

    ! Initialize the numerics
    subroutine InitializeGeneralNumParamsINT(num)

        ! Description
        !============
        ! Initialize the general numerical parameters

        ! Declare variables
        !==================
        ! Arguments
        class(NumUDT)            :: num
    
        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Set defaults
        call num%SetDefaultNumParams()

        ! Override with user settings (to be implemented)

    end subroutine

    !------------------------------------------------------------------!
    !                             KKT OPTIONS                          !
    !------------------------------------------------------------------!
    ! Set the KKT numerical parameters
    subroutine SetDefaultNumParamsKKTINT(num)

        ! Description
        !============
        ! Set default numerical parameters tol, itmax, verbosity

        ! Declare variables
        !==================
        ! Arguments
        class(NumKKTUDT)            :: num
    
        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Set defaults
        !=============
        ! General numerics
        call num%SetDefaultNumParams()

        ! Specifics for KKT numerics
        num%rxf             = 1e2
        num%rxfdesign       = 1
        num%rxfdec          = 0.98
        num%rxfmin          = 2e-3

    end subroutine

    ! Initialize the numerics
    subroutine InitializeNumParamsKKTINT(num)

        ! Description
        !============
        ! Initialize the general numerical parameters

        ! Declare variables
        !==================
        ! Arguments
        class(NumKKTUDT)            :: num
    
        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Set defaults 
        call num%SetDefaultNumParamsKKT()

        ! Override with user settings (to be implemented)

    end subroutine

end module