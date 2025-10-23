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
        ! - meth: can be 'global' to just use the global x-y coordinates system 

        character(:), allocatable   :: type1
        character(:), allocatable   :: type2
        character(:), allocatable   :: meth
        integer(I8), allocatable    :: cNv(:), cNvP(:,:)

        real(R8), allocatable       :: coef(:,:)
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
        integer(I8) :: i, j, k, n, m, info
        real(R8) :: det, det1, det2, det3
        real(R8), allocatable :: A(:,:), AT(:,:), ATA_dummy(:,:), w(:), &
            temp(:), sol(:), C(:,:), invC(:,:), A_test(:,:)
        real(R8), allocatable :: prefact(:)

        ! Initialize
        allocate(prefact(0:deriv))
        do i = 0, deriv
            prefact(i) = 1.0_R8/factorial(i)
        end do
        n = sum((/(i, i = 2, deriv + 1)/)) ! number of elements
        
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

            ! Set weight of first vertex (dx=0, dy=0) to double the max of the array
            w(1) = maxval(w(2:size(dx)))*2

            ! Allocate A matrix
            allocate(A(size(dx), n))       
            A = 0
            A_test = A

            if (deriv == 2) then

                ! A matrix
                A(:, 1) = w*dx ! dphidx
                A(:, 2) = w*dy ! dphidy
                A(:, 3) = prefact(2)*w*dx**2 ! dphidx2
                A(:, 4) = w*dx*dy ! dphidxdy
                A(:, 5) = prefact(2)*w*dy**2 ! dphidy2


            else if (deriv == 3) then

                ! A matrix
                A(:, 1) = w*dx ! dphidx
                A(:, 2) = w*dy ! dphidy
                A(:, 3) = prefact(2)*w*dx**2 ! dphidx2
                A(:, 4) = w*dx*dy   ! dphidxdy
                A(:, 5) = prefact(2)*w*dy**2 ! dphidy2
                A(:, 6) = prefact(3)*w*dx**3 ! dphidx3
                A(:, 7) = prefact(2)*w*dx**2*dy ! dphidx2dy
                A(:, 8) = prefact(2)*w*dx*dy**2 ! dphidxdy2
                A(:, 9) = prefact(3)*w*dy**3 ! dphidy3

            else if (deriv == 4) then
                
                ! A matrix 
                A(:, 1) = w*dx ! dphidx
                A(:, 2) = w*dy ! dphidy
                A(:, 3) = prefact(2)*w*dx**2 ! dphidx2
                A(:, 4) = w*dx*dy   ! dphidxdy
                A(:, 5) = prefact(2)*w*dy**2 ! dphidy2
                A(:, 6) = prefact(3)*w*dx**3 ! dphidx3
                A(:, 7) = prefact(2)*w*dx**2*dy ! dphidx2dy
                A(:, 8) = prefact(2)*w*dx*dy**2 ! dphidxdy2
                A(:, 9) = prefact(3)*w*dy**3 ! dphidy3
                A(:,10) = prefact(4)*w*dx**4  ! dphidx4          
                A(:,11) = prefact(3)* w*dx**3*dy ! dphidx3dy 
                A(:,12) = prefact(2)*prefact(2)* w*dx**2*dy**2 ! dphidx2dy2           
                A(:,13) = prefact(3)*w*dx*dy**3 ! dphidxdy3           
                A(:,14) = prefact(4)*w*dy**4  ! dphidy4        

            else if (deriv == 5) then

                ! A matrix 
                A(:, 1) = w*dx ! dphidx
                A(:, 2) = w*dy ! dphidy
                A(:, 3) = prefact(2)*w*dx**2 ! dphidx2
                A(:, 4) = w*dx*dy   ! dphidxdy
                A(:, 5) = prefact(2)*w*dy**2 ! dphidy2
                A(:, 6) = prefact(3)*w*dx**3 ! dphidx3
                A(:, 7) = prefact(2)*w*dx**2*dy ! dphidx2dy
                A(:, 8) = prefact(2)*w*dx*dy**2 ! dphidxdy2
                A(:, 9) = prefact(3)*w*dy**3 ! dphidy3
                A(:,10) = prefact(4)*w*dx**4  ! dphidx4          
                A(:,11) = prefact(3)* w*dx**3*dy ! dphidx3dy 
                A(:,12) = prefact(2)*prefact(2)* w*dx**2*dy**2 ! dphidx2dy2           
                A(:,13) = prefact(3)*w*dx*dy**3 ! dphidxdy3           
                A(:,14) = prefact(4)*w*dy**4  ! dphidy4   
                A(:,15) = prefact(5)* w*dx**5  ! dphidx5          
                A(:,16) = prefact(4)* w*dx**4*dy  ! dphidx4dy
                A(:,17) = prefact(3)*prefact(2)* w*dx**3*dy**2  ! dphidx3dy2
                A(:,18) = prefact(2)*prefact(3)* w*dx**2*dy**3  ! dphidx2dy3
                A(:,19) = prefact(4)* w*dx*dy**4  ! dphidxdy4
                A(:,20) = prefact(5)* w*dy**5  ! dphidy5
                
            else if (deriv == 6) then

                ! A matrix 
                A(:, 1) = w*dx ! dphidx
                A(:, 2) = w*dy ! dphidy
                A(:, 3) = prefact(2)*w*dx**2 ! dphidx2
                A(:, 4) = w*dx*dy   ! dphidxdy
                A(:, 5) = prefact(2)*w*dy**2 ! dphidy2
                A(:, 6) = prefact(3)*w*dx**3 ! dphidx3
                A(:, 7) = prefact(2)*w*dx**2*dy ! dphidx2dy
                A(:, 8) = prefact(2)*w*dx*dy**2 ! dphidxdy2
                A(:, 9) = prefact(3)*w*dy**3 ! dphidy3
                A(:,10) = prefact(4)*w*dx**4  ! dphidx4          
                A(:,11) = prefact(3)* w*dx**3*dy ! dphidx3dy 
                A(:,12) = prefact(2)*prefact(2)* w*dx**2*dy**2 ! dphidx2dy2           
                A(:,13) = prefact(3)*w*dx*dy**3 ! dphidxdy3           
                A(:,14) = prefact(4)*w*dy**4  ! dphidy4   
                A(:,15) = prefact(5)* w*dx**5  ! dphidx5          
                A(:,16) = prefact(4)* w*dx**4*dy  ! dphidx4dy
                A(:,17) = prefact(3)*prefact(2)* w*dx**3*dy**2  ! dphidx3dy2
                A(:,18) = prefact(2)*prefact(3)* w*dx**2*dy**3  ! dphidx2dy3
                A(:,19) = prefact(4)* w*dx*dy**4  ! dphidxdy4
                A(:,20) = prefact(5)* w*dy**5  ! dphidy5
                A(:,21) = prefact(6)* w*dx**6  ! dphidx6
                A(:,22) = prefact(5)* w*dx**5*dy  ! dphidx5dy
                A(:,23) = prefact(4)*prefact(2)* w*dx**4*dy**2  ! dphidx4dy2
                A(:,24) = prefact(3)*prefact(3)* w*dx**3*dy**3  ! dphidx3dy3
                A(:,25) = prefact(2)*prefact(4)* w*dx**2*dy**4  ! dphidx2dy4
                A(:,26) = prefact(5)* w*dx*dy**5  ! dphidxdy5
                A(:,27) = prefact(6)* w*dy**6  ! dphidy6

            else     
                call gdErrorHandler('ComputeATA: not implemented')
            end if 

            k = 0
            do m = 1, deriv
                do j = 0, m
                    i = m - j
                    k = k + 1
                    A_test(:,k) = prefact(i)*prefact(j)*w*dx**i*dy**j
                end do
            end do

            if (.not. all(A_test == A)) call gdErrorHandler('Mistake')

            ! Construct C matrix (transpose(A)*A)
            allocate(C(n,n), invC(n,n))
            AT = transpose(A)
            C = matmul(AT, A)

            ! Compute inverse of C (symmetric)
            allocate(temp(n))
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
    end subroutine

       ! Factorial computation
    integer(I16) function factorial(n)

        ! Description
        !============
        ! Simple factorial computation routine. Not accurate for n >> 10
        ! likely. 

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)             :: n

        ! Loop
        integer(I8)                         :: i 

        ! Compute
        !========
        ! Initialize
        factorial = 1
        if (n == 0) then 
            ! Exit
            return 
        end if 
        do i = 1, n 
            factorial = factorial * i
        end do 

    end function 


end module