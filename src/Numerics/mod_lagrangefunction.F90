!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module provides functionality for the use of lagrange basis
! functions or lagrange polynomials. This polynomial is given by 
!
!       L = Sum_i l_i,
!
! where l_i is the i-th Lagrange polynomial. This polynomial is given by:
!
!       l_i(s) = prod( (x - x_m)(x_i - x_m) ) for m = 0..k, m != i
!
! If desired, one can also write the Lagrange polynomial in the form of
! l_i = prod( 1/(x_i - x_m) ) sum_i a_i x^i, which is convenient for
! integration/differentiation. Here, the a_i are computed recursively
! (perhaps there's a closed and elegant formulation but haven't found that
! one yet). If we start with the 'constant' part (m = 0) and move further
! to m = M, then, we can find each set of coefficients of m+1 as function
! of the coefficients of m (note: we assume here that i is skipped, but do
! not reflect this in notation):
!
!       a_j^m+1 = a_j-1^m (-x_m+1) + (a_j^m)
!
! For example:
!
!       m   |   a_0^m   a_1^m-1     ...
!       0   |   1       0     0     ...
!       1   |   1       -x_1    0
!       2   |   1       -(x_1 + x_2)    x_1 x_2
!       3   |   1       -(x_1 + x_2 + x_3)  x_1 x_2 + x_1 x_3 + x_2 x_3     -x_1 x_2 x_3
!
! For the definition of l_i above, this would result for example for l_1
! and M = 3 in:
!
!       m   |   a_0^m   a_1^m-1     ...
!       0   |   1       0     0     ...
!       ----------------------------------
!       2   |   1       -(x_2)
!       3   |   1       -(x_2 + x_3)  x_2 x_3

module mod_lagrangefunctions

    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_errorhandler

    implicit none 
    private 
    public :: ConstructLagrangianBasisFunctions

    !==================================================================!
    !                                                                  !
    !                            TYPES                                 !
    !                                                                  !
    !==================================================================!

    contains 

    !==================================================================!
    !                                                                  !
    !                          ROUTINES                                !
    !                                                                  !
    !==================================================================!

    ! Lagrangian basis function coefficient construction
    subroutine ConstructLagrangianBasisFunctions(order, xi, a, aint)

        ! Description
        !============
        ! Construct coefficients of the Lagrangian basis functions. Here, we
        ! compute the polynomial in the form l_i = a_i x^i, where the coefficients
        ! a_i are computed by the 'recursive' relationship (see also explanation 
        ! above). It is assumed here that the spacing is
        ! distributed over the interval [0, 1] (which should be included in the
        ! xi argument), and that the amount of points equals the order-1.

        ! Note: additionally, we also compute the weights for the integrated
        ! lagrangian (only integrated once) for convenience in the rest of this
        ! routine.

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)                :: xi(:)
        integer(I8), intent(in)             :: order 
        real(R8), allocatable, intent(out)  :: a(:, :), aint(:, :)

        ! Auxiliary
        real(R8)                            :: tp 
        real(R8), allocatable, dimension(:) :: ai, ai_old

        ! Loop
        integer(I8)                         :: i, j, k

        ! Initialize
        !===========
        ! Check
        if ((order+1) /= size(xi)) then 
            call gdErrorHandler('ConstructLagrangianBasisFunctions: ' // & 
                'order and xi dimensions do not match')
        end if 
        if (allocated(a)) then 
            deallocate(a)
        end if 
        if (allocated(aint)) then 
            deallocate(aint)
        end if 
        allocate(a(size(xi), size(xi)), aint(size(xi), size(xi)), &
            ai(size(xi)))
        a = 0
        aint = 0

        ! Construct coefficients
        !=======================
        do i = 0, order
            ! Initialize weights
            ai = 0
            ai(1) = 1
            do j = 0, order
                ! Skip if j == i
                if (j == i) then 
                    cycle
                end if 
                
                ! Build recursively
                ai_old = ai
                do k = 1, order
                    ai(k+1) = ai_old(k)*(-xi(j+1)) + ai_old(k+1)
                end do 
            end do 
            
            ! Build product
            tp = 1
            do j = 0, order
                if (j == i) then 
                    cycle 
                end if 
                tp = tp*1/( xi(i+1) - xi(j+1) )
            end do
            
            ! Add
            a(i+1, :) = tp*ai(order+1:1:-11); ! because x^i will go from i = 0..M, not M .. 0
            aint(i+1, :) = a(i+1, :)/([(k, k = 1, order+1)])
        end do 

    end subroutine 

end module 
