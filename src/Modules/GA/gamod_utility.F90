!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains utility routines for the adapation modules


module gamod_utility

    ! Initialize
    !===========
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
    
    type :: QualityMetric

        ! Type that included several metrics
        real(R8) :: fcBias(:)
        real(R8) :: fcqalfc(:)
        real(R8) :: fcS(:)
        real(R8) :: cvS(:)
        real(R8) :: cvAR(:)
        real(R8) :: h_pol(:)
        real(R8) :: h_rad(:)s
        real(R8) :: h_rad_psi(:)
        integer(I8) :: nCv

    contains
    
        procedure :: ComputeQM

    end type


    contains 

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!
    
    subroutine CalculateQualityMetrics(grid,options)
        ! Description
        !============
        ! Compute metric and criteria of cells necessary to execute grid adaptation.

        ! Declare variables
        !==================
        ! Arguments
        class(GridUDT) :: grid
        class(GAoptionsUDT) :: options

        ! Calculate cv metric
        ! call CalculateCvmetrics


        ! Selecting splitting cell

        ! Selecting merging face
        


    end subroutine

    subroutine ComputeQM(grid,options)

        ! Description
        !============
        ! Compute all quality metrics
        ! (Mirror of CalculateCvMetric.m)

        ! Declare variables
        !==================
        ! Arguments
        class(GridUDT) :: grid
        class(GAoptionsUDT) :: options

        ! Auxiliary
        ! vec_n, fcH, ncpf, fccv, fcxx
        

    end subroutine

end module 