!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains gradient computation methods on finite volume meshes.

module mod_gradient

    ! Load modules
    use mod_precision
    use mod_errorhandler

    implicit none

    !==================================================================!
    !                                                                  !
    !                            TYPES                                 !
    !                                                                  !
    !==================================================================!
    
    ! Gradient reconstruction type
    type, abstract :: GradientReconstructionUDT

        ! Description
        !============
        ! The GradientReconstruction type contains the following data
        ! - cNv: stencil array 
        ! - cNvP: pointer for stencil array
        ! - type1: can be 'cell', 'face', 'vert' and indicates for which
        !           connectivity type to compute gradients
        ! - type2: can be 'cell', 'face', 'vert' and indicates for which
        !           connectivity type to draw information for computation
        ! - meth: can be 

        character(:), allocatable :: type1
        character(:), allocatable :: type2
        character(:), allocatable :: meth
        integer(I8), allocatable :: cNv(:), cNvP(:,:)

        real(R8), allocatable    :: invA(:,:)

    contains

        ! Setting parameters
        !procedure(SetParametersINT), deferred   :: SetParameters

        ! Constructor
        procedure(ConstructINT), deferred       :: Construct

        ! Evaluate
        procedure(EvaluateINT), deferred        :: Evaluate

    end type

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!
    ! Abstract interface
    abstract interface

        ! Constructor
        subroutine ConstructINT(GR)
            import :: GradientReconstructionUDT
            class(GradientReconstructionUDT) :: GR
        end subroutine

        ! Evaluate
        subroutine EvaluateINT(GR, v)
            import :: GradientReconstructionUDT, R8
            class(GradientReconstructionUDT) :: GR
            real(R8), intent(in) :: v(:)
        end subroutine

    
    end interface

    contains



end module