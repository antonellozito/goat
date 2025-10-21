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
    use mod_linearsolverinterface

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
        ! - meth: can be TODO 

        character(:), allocatable   :: type1
        character(:), allocatable   :: type2
        character(:), allocatable   :: meth
        integer(I8), allocatable    :: cNv(:), cNvP(:,:)

        real(R8), allocatable       :: invA(:,:)
        real(R8), allocatable       :: w(:)
        integer(I8)                 :: deriv


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
        subroutine EvaluateINT(GR, v, deriv_vals)
            import :: GradientReconstructionUDT, R8
            class(GradientReconstructionUDT) :: GR
            real(R8), intent(in) :: v(:)
            real(R8), allocatable, intent(out) :: deriv_vals(:,:)
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
        real(R8), allocatable :: AT(:,:), ATA_dummy(:,:)

        ! A transpose
        allocate(AT(2, size(dx)))
        AT(1,:) = dx
        AT(2,:) = dy

        ! Construct C matrix (transpose(A)*A)
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
        allocate(ATA_dummy(2, size(dx)))
        ATA_dummy = matmul(E, AT)

        ! Do transpose for better memory
        ATA = transpose(ATA_dummy)
     

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
        real(R8), allocatable :: AT(:,:), ATA_dummy(:,:)

        ! Transpose A
        allocate(AT(3, size(dx)))
        AT(1,:) = dx
        AT(2,:) = dy
        AT(3,:) = 1

        ! Construct C matrix sparse (transpose(A)*A)
        c11 = sum(dx**2)
        c12 = sum(dx*dy)
        c13 = sum(dx)
        c22 = sum(dy**2)
        c23 = sum(dy)
        c33 = size(dx)

        ! Compute inverse of C
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

        ! Multiply with transpose A
        allocate(ATA_dummy(3, size(dx)))
        ATA_dummy = matmul(E,AT)

        ! Return transpose for better memory
        ATA = transpose(ATA_dummy)

    end subroutine

    subroutine ComputeATA5(dx, dy, ATA)

        ! Description
        !============
        ! Compute coefficient for second order gradient reconstruction.
        ! The arrangment is: dphidx, dphidy, d2phidx2, d2phidy2, d2phidxdy

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in) :: dx(:), dy(:)
        real(R8), allocatable, intent(out) :: ATA(:,:)

        ! Auxiliary
        integer(I8) :: info
        real(R8) :: C(5, 5), invC(5, 5)
        real(R8), allocatable :: A(:,:), AT(:,:), ATA_dummy(:,:), w(:), &
            temp(:), sol(:)

        ! Compute weigths (w_ij = 1/(dx**2 + dy**2))
        w = 1/(dx**2 + dy**2)

        ! A transpose
        allocate(A(size(dx), 5))
        A(:, 1) = w*dx
        A(:, 2) = w*dy
        A(:, 3) = 0.5_R8*w*dx**2
        A(:, 4) = 0.5_R8*w*dy**2
        A(:, 5) = w*dx*dy

        ! Construct C matrix (transpose(A)*A)
        AT = transpose(A)
        C = matmul(AT, A)

        ! Compute inverse of C (symmetric)
        allocate(temp(size(dx)))
        temp = 0
        sol = temp
        call SolveDenseLinearSystemDI(C, temp, sol, info, invC)

        ! Multiply with transpose A
        allocate(ATA_dummy(5, size(dx)))
        ATA_dummy = matmul(invC, AT)

        ! Do transpose for bettery memory
        ATA = transpose(ATA_dummy)

    end subroutine

    subroutine ComputeATA(dx, dy, deriv, int, ATA)

        ! Description
        !============
        ! Compute coefficient for second order gradient reconstruction.
        ! The arrangment is: dphidx, dphidy, d2phidx2, d2phidy2, d2phidxdy, ...
        ! Deriv is GR%deriv
        ! Int is a logical to indicate also to give the interpolated value

        ! Declare variables
        !==================        
        ! Arguments
        real(R8), intent(in) :: dx(:), dy(:)
        real(R8), allocatable, intent(out) :: ATA(:,:)   
        integer(I8), intent(in)            :: deriv
        logical, intent(in)                 :: int              
        
        ! Auxiliary
        integer(I8) :: n, info
        real(R8) :: det, det1, det2, det3
        real(R8), allocatable :: A(:,:), AT(:,:), ATA_dummy(:,:), w(:), &
            temp(:), sol(:), C(:,:), invC(:,:)
        
        if (deriv == 1 .and. .not. int) then

            ! Number of arrays
            n = 2
            allocate(C(n,n), invC(n,n))

            ! A transpose
            allocate(AT(n, size(dx)))
            AT(1,:) = dx
            AT(2,:) = dy

            ! Construct C matrix (transpose(A)*A)
            C = 0
            C(1,1) = sum(dx**2)
            C(1,2) = sum(dx*dy)
            C(2,2) = sum(dy**2) 

            ! Compute inverse of C
            det = C(1,1)*C(2,2) - C(1,2)**2
            invC(1, 1) =  C(2,2)
            invC(1, 2) = -C(1,2)
            invC(2, 1) = invC(1, 2)
            invC(2, 2) = C(1,1)
            invC = invC/det             

        else if (deriv == 1 .and. int) then

            ! Number of arrays
            n = 3
            allocate(C(n,n))

            ! Transpose A
            allocate(AT(3, size(dx)))
            AT(1,:) = dx
            AT(2,:) = dy
            AT(3,:) = 1

            ! Construct C matrix sparse (transpose(A)*A)
            C = 0
            C(1,1) = sum(dx**2)
            C(1,2) = sum(dx*dy)
            C(1,3) = sum(dx)
            C(2,2) = sum(dy**2)
            C(2,3) = sum(dy)
            C(3,3) = size(dx)

            ! Compute inverse of C
            det1 = C(3,3)*C(2,2) - C(2,3)**2
            det2 = C(3,3)*C(1,2) - C(2,3)*C(1,3)
            det3 = C(2,3)*C(1,2) - C(2,2)*C(1,3)
            det = C(1,1)*det1 - C(1,2)*det2 + C(1,3)*det3

            invC(1,1) = det1
            invC(1,2) = -det2
            invC(1,3) = det3
            invC(2,1) = invC(1,2)
            invC(2,2) = C(3,3)*C(1,1) - C(1,3)**2
            invC(2,3) = -(C(2,3)*C(1,1) - C(1,2)*C(1,3))
            invC(3,1) = invC(1,3)
            invC(3,2) = invC(2,3)
            invC(3,3) = C(2,2)*C(1,1)-C(1,2)**2
            invC = invC/det            

        else if (deriv .gt. 1 .and. .not. int) then

            ! Compute weigths (w_ij = 1/(dx**2 + dy**2))
            w = 1/(dx**2 + dy**2)

            if (deriv == 2) then

                ! Number of arrays
                n = 5

                ! A matrix
                allocate(A(size(dx), n))
                A(:, 1) = w*dx
                A(:, 2) = w*dy
                A(:, 3) = 0.5_R8*w*dx**2
                A(:, 4) = 0.5_R8*w*dy**2
                A(:, 5) = w*dx*dy

            else if (deriv == 3) then

                ! Number of arrays
                n = 9

                ! A matrix
                allocate(A(size(dx), n))
                A(:,1) = w*dx ! dphidx
                A(:, 2) = w*dy ! dphidy
                A(:, 3) = 0.5_R8*w*dx**2 ! dphidx2
                A(:, 4) = 0.5_R8*w*dy**2 ! dphidy2
                A(:, 5) = w*dx*dy   ! dphidxdy
                A(:, 6) = 1.0_R8/6.0_R8 * w*dx**3 ! dphidx3
                A(:, 7) = 1.0_R8/6.0_R8 * w*dy**3 ! dphidx3
                A(:, 8) = 0.5_R8 * w*dx**2*dy ! dphidx2dy
                A(:, 9) = 0.5_R8 * w*dx*dy**2 ! dphidxdy2

            else if (deriv == 4) then
                
                ! Number of arrays
                n = 14

                ! A matrix 
                allocate(A(size(dx),n))
                A(:,1) = w*dx ! dphidx
                A(:, 2) = w*dy ! dphidy
                A(:, 3) = 0.5_R8*w*dx**2 ! dphidx2
                A(:, 4) = 0.5_R8*w*dy**2 ! dphidy2
                A(:, 5) = w*dx*dy   ! dphidxdy
                A(:, 6) = 1.0_R8/6.0_R8 * w*dx**3 ! dphidx3
                A(:, 7) = 1.0_R8/6.0_R8 * w*dy**3 ! dphidx3
                A(:, 8) = 0.5_R8 * w*dx**2*dy ! dphidx2dy
                A(:, 9) = 0.5_R8 * w*dx*dy**2 ! dphidxdy2
                A(:,10) = 1.0_R8/24.0_R8 * w*dx**4  ! dphidx4          
                A(:,11) = 1.0_R8/24.0_R8 * w*dy**4  ! dphidy4  
                A(:,12) = 1.0_R8/6.0_R8 * w*dx**3*dy ! dphidx3dy           
                A(:,13) = 1.0_R8/4.0_R8 * w*dx**2*dy**2 ! dphidx2dy2           
                A(:,14) = 1.0_R8/6.0_R8 * w*dx*dy**3 ! dphidxdy3           

            else if (deriv == 5) then

                ! Number of arrays
                n = 20

                ! A matrix 
                allocate(A(size(dx),n))
                A(:,1) = w*dx ! dphidx
                A(:, 2) = w*dy ! dphidy
                A(:, 3) = 0.5_R8*w*dx**2 ! dphidx2
                A(:, 4) = 0.5_R8*w*dy**2 ! dphidy2
                A(:, 5) = w*dx*dy   ! dphidxdy
                A(:, 6) = 1.0_R8/6.0_R8 * w*dx**3 ! dphidx3
                A(:, 7) = 1.0_R8/6.0_R8 * w*dy**3 ! dphidx3
                A(:, 8) = 0.5_R8 * w*dx**2*dy ! dphidx2dy
                A(:, 9) = 0.5_R8 * w*dx*dy**2 ! dphidxdy2
                A(:,10) = 1.0_R8/24.0_R8 * w*dx**4  ! dphidx4          
                A(:,11) = 1.0_R8/24.0_R8 * w*dy**4  ! dphidy4  
                A(:,12) = 1.0_R8/6.0_R8 * w*dx**3*dy ! dphidx3dy           
                A(:,13) = 1.0_R8/4.0_R8 * w*dx**2*dy**2 ! dphidx2dy2           
                A(:,14) = 1.0_R8/6.0_R8 * w*dx*dy**3 ! dphidxdy3   
                A(:,15) = 1.0_R8/120.0_R8 * w*dx**5  ! dphidx5          
                A(:,16) = 1.0_R8/120.0_R8 * w*dy**5  ! dphidy5
                A(:,17) = 1.0_R8/24.0_R8 * w*dx**4*dy  ! dphidx4dy
                A(:,18) = 1.0_R8/12.0_R8 * w*dx**3*dy**2  ! dphidx3dy2
                A(:,19) = 1.0_R8/12.0_R8 * w*dx**2*dy**3  ! dphidx2dy3
                A(:,20) = 1.0_R8/24.0_R8 * w*dx*dy**4  ! dphidxdy4
                
            else if (deriv == 6) then

                ! Number of arrays
                n = 27

                ! A matrix 
                allocate(A(size(dx),n))
                A(:,1) = w*dx ! dphidx
                A(:, 2) = w*dy ! dphidy
                A(:, 3) = 0.5_R8*w*dx**2 ! dphidx2
                A(:, 4) = 0.5_R8*w*dy**2 ! dphidy2
                A(:, 5) = w*dx*dy   ! dphidxdy
                A(:, 6) = 1.0_R8/6.0_R8 * w*dx**3 ! dphidx3
                A(:, 7) = 1.0_R8/6.0_R8 * w*dy**3 ! dphidx3
                A(:, 8) = 0.5_R8 * w*dx**2*dy ! dphidx2dy
                A(:, 9) = 0.5_R8 * w*dx*dy**2 ! dphidxdy2
                A(:,10) = 1.0_R8/24.0_R8 * w*dx**4  ! dphidx4          
                A(:,11) = 1.0_R8/24.0_R8 * w*dy**4  ! dphidy4  
                A(:,12) = 1.0_R8/6.0_R8 * w*dx**3*dy ! dphidx3dy           
                A(:,13) = 1.0_R8/4.0_R8 * w*dx**2*dy**2 ! dphidx2dy2           
                A(:,14) = 1.0_R8/6.0_R8 * w*dx*dy**3 ! dphidxdy3   
                A(:,15) = 1.0_R8/120.0_R8 * w*dx**5  ! dphidx5          
                A(:,16) = 1.0_R8/120.0_R8 * w*dy**5  ! dphidy5
                A(:,17) = 1.0_R8/24.0_R8 * w*dx**4*dy  ! dphidx4dy
                A(:,18) = 1.0_R8/12.0_R8 * w*dx**3*dy**2  ! dphidx3dy2
                A(:,19) = 1.0_R8/12.0_R8 * w*dx**2*dy**3  ! dphidx2dy3
                A(:,20) = 1.0_R8/24.0_R8 * w*dx*dy**4  ! dphidxdy4
                A(:,21) = 1.0_R8/720.0_R8 * w*dx**6  ! dphidx6
                A(:,22) = 1.0_R8/720.0_R8 * w*dy**6  ! dphidy6
                A(:,23) = 1.0_R8/120.0_R8 * w*dx**5*dy  ! dphidx5dy
                A(:,24) = 1.0_R8/48.0_R8 * w*dx**4*dy**2  ! dphidx4dy2
                A(:,25) = 1.0_R8/36.0_R8 * w*dx**3*dy**3  ! dphidx3dy3
                A(:,26) = 1.0_R8/48.0_R8 * w*dx**2*dy**4  ! dphidx2dy4
                A(:,27) = 1.0_R8/120.0_R8 * w*dx*dy**5  ! dphidxdy5

            else     
                call gdErrorHandler('ComputeATA: not implemented')
            end if 

            ! Construct C matrix (transpose(A)*A)
            allocate(C(n,n), invC(n,n))
            AT = transpose(A)
            C = matmul(AT, A)

            ! Compute inverse of C (symmetric)
            allocate(temp(size(dx)))
            temp = 0
            sol = temp
            call SolveDenseLinearSystemDI(C, temp, sol, info, invC)

        else 

            call gdErrorHandler('ComputeATA: not implemented')

        end if

        ! Multiply with transpose A
        allocate(ATA_dummy(n, size(dx)))
        ATA_dummy = matmul(invC, AT)

        ! Do transpose for better memory
        ATA = transpose(ATA_dummy)
    end 


end module