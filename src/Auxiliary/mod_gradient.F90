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
        !procedure(ConstructINT), deferred       :: Construct

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
        !subroutine ConstructINT(GR)
        !    import :: GradientReconstructionUDT
        !    class(GradientReconstructionUDT) :: GR
        !end subroutine

        ! Evaluate
        subroutine EvaluateINT(GR, v)
            import :: GradientReconstructionUDT, R8
            class(GradientReconstructionUDT) :: GR
            real(R8), intent(in) :: v(:)
        end subroutine

    
    end interface

    contains

    !==================================================================!
    !                                                                  !
    !                             ROUTINES                             !
    !                                                                  !
    !==================================================================!

    subroutine ComputeATA2(dx, dy, ATA)

        ! Description
        !============
        ! Compute the coefficients for Gradient reconstruction for only the gradients, so no interpolated field value.

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in) :: dx(:), dy(:)
        real(R8), allocatable, intent(out) :: ATA(:,:)

        ! Auxiliary
        real(R8) :: c11, c12, c22, E(2,2), det
        real(R8), allocatable :: AT(:,:)

        ! A transpose
        allocate(AT(2,size(dx)))
        AT(1,:) = dx
        AT(2,:) = dy
            
        ! Construct C matrix
        c11 = sum(dx**2)
        c12 = sum(dx*dy)
        c22 = sum(dy**2)

        ! Compute inverse of C
        det = c11*c22 - c12**2
        E(1, 1) = c22
        E(1, 2) = -c12
        E(2, 1) = E(1, 2)
        E(2, 2) = c11
        E = E/det       

        ! Multiply with transpose A
        ATA = E*AT

    end subroutine

    subroutine ComputeATA3(dx, dy, ATA)

        ! Description
        !============
        ! Compute coefficient for gradient reconstruction, with extra interpolation result.

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in) :: dx(:), dy(:)
        real(R8), allocatable, intent(out) :: ATA(:,:)

        ! Auxiliary
        real(R8) :: c11, c12, c13, c22, c23, c33, E(3,3), det1, det2, det3, det
        real(R8), allocatable :: AT(:,:)

        allocate(AT(3, size(dx)))

        AT(1,:) = dx
        AT(2,:) = dy
        AT(3,:) = 1

        c11 = sum(dx**2)
        c12 = sum(dx*dy)
        c13 = sum(dx)
        c22 = sum(dy**2)
        c23 = sum(dy)
        c33 = size(dx)

        det1 = c33*c22 - c23**2
        det2 = c33*c12 - c23*c13
        det3 = c23*c12 - c22*c13
        det = c11*det1 - c12*det2 + c13*det3

        E(1,1) = det1
        E(1,2) = -det2
        E(1,3) = det3
        E(2,1) = E(1,2)
        E(2,2) = c33*c11 - c13**2
        E(2,3) = -(c23*c11 - c12*c13)
        E(3,1) = E(1,3)
        E(3,2) = E(2,3)
        E(3,3) = c22*c11-c12*c12
        E = E/det

        ATA = E*AT

    end subroutine


end module