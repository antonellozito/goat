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
    use optmod_optimizationengine
    use optmod_designvariables
    use optmod_costfunction
    use optmod_constraints
    use optmod_state

    ! The usual
    implicit none
    save
    public 

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!


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
        use gdmod_types
        use gdmod_userinput

        ! The usual
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationEngineGDUDT)      :: optimizationdriver 

        ! Loop variables

        ! Auxiliary variables 
        type(DesignOptionsUDT)              :: designoptions

        ! Data

        ! Design variables
        !=================
        ! Set the design options
        call SetDesignOptions(designoptions)


    end subroutine

end module