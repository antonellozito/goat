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

    subroutine ComputeATA(dx, dy, deriv, int, ATA)

        ! Description
        !============
        ! Compute coefficient for second order gradient reconstruction.
        ! The arrangment is: dphidx, dphidy, d2phidx2, d2phidxdy, d2phidy2,...
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
        integer(I8) :: i, j, k, n, m, info, start_d
        real(R8), allocatable :: A(:,:), AT(:,:), ATA_dummy(:,:), w(:), &
            temp(:), sol(:), C(:,:), invC(:,:)
        real(R8), allocatable :: prefact(:)

        ! Initialize
        allocate(prefact(0:deriv))
        do i = 0, deriv
            prefact(i) = 1.0_R8/factorial(i)
        end do

        ! Determine number of column is A matrix
        if (.not.int) then
            n = sum((/(i, i = 2, deriv + 1)/)) ! number of elements
            start_d = 1
        else if (int) then
            n = sum((/(i, i = 1, deriv + 1)/)) ! number of elements
            start_d = 0
        end if
        
        if (deriv == 1) then

            ! Weight
            allocate(w(size(dx)))
            w = 1.0_R8

        else if (deriv .gt. 1) then

            ! Compute weigths (w_ij = 1/(dx**2 + dy**2))
            w = 1/(dx**2 + dy**2)

            ! Set weight of first vertex (dx=0, dy=0) to double the max of the array
            w(1) = maxval(w(2:size(dx)))*2

        end if

        ! Allocate A matrix
        allocate(A(size(dx), n))       
        A = 0
        k = 0
        do m = start_d, deriv
            do j = 0, m
                i = m - j
                k = k + 1
                A(:,k) = prefact(i)*prefact(j)*w*dx**i*dy**j
            end do
        end do

        ! Construct C matrix (transpose(A)*A)
        allocate(C(n,n), invC(n,n))
        AT = transpose(A)
        C = matmul(AT, A)

        ! Compute inverse of C (symmetric)
        allocate(temp(n))
        temp = 0
        sol = temp
        call SolveDenseLinearSystemDI(C, temp, sol, info, invC)


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