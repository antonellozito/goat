!======================================================================!
!                                                                      !
!                            DOCUMENTATION                             !
!                                                                      !
!======================================================================!

! Short description
!==================
! This is the main driver to run tests for several components/features 
! of the goat module. Different tests are available 

module GOAT_tests 

    implicit none
    save 

    contains 

    !------------------------------------------------------------------!
    !                         General drivers                          !
    !------------------------------------------------------------------!
    ! Subroutine to run all tests
    subroutine RunAllTests()

        ! Interpolant testing
        call TestStructuredInterpolant2D()

    end subroutine

    !------------------------------------------------------------------!
    !                         Output routines                          !
    !------------------------------------------------------------------!
    ! Write starting header
    subroutine DisplayTestStart(testname)

        character(*),  intent(in)       :: testname 

        print *, '%---------------------------------------------------%'
        print *, '%                      TEST START                   %'
        print *, '%---------------------------------------------------%'
        print *, 'Test name: ' // testname 

    end subroutine

    ! Ending header
    subroutine DisplayTestEnd()

        print *, '%---------------------------------------------------%'
        print *, '%                      TEST END                     %'
        print *, '%---------------------------------------------------%'

    end subroutine

    !------------------------------------------------------------------!
    !                         Test routines                            !
    !------------------------------------------------------------------!
    ! Structured interpolant testing
    subroutine TestStructuredInterpolant2D()

        ! Modules
        !========
        use Interpolant

        implicit none
        save 

        ! Declare variables
        !==================
        ! Arguments

        ! Loop variables

        ! Auxiliary
        type(StructuredInterpolant2DUDT)    :: interp

        integer(I8)                         :: nx, ny, nq, nres
        real(R8)                            :: Lx, Ly

        character(:), allocatable           :: meth 
        integer(I8), allocatable            :: reshuffle(:)
        real(R8), allocatable               :: xgv(:), ygv(:), &
            v(:, :), a(:), vq(:, :), xq(:), yq(:), vqan(:, :), &
            dv(:, :), temp(:), relerr(:)

        ! Loop
        integer(I8)                         :: i, j, k 

        ! Initialize
        !===========
        ! Display
        call DisplayTestStart('TestStructuredInterpolant2D')

        ! Build structured grid
        nx = 100 ! number of cells, not vertices!
        ny = 200
        Lx = 10
        Ly = 5 
        allocate(xgv(nx+1), ygv(ny+1), v(nx+1, ny+1))
        xgv = Lx*[(k, k=0, nx)]/nx 
        ygv = Ly*[(k, k=0, ny)]/ny

        ! Set evaluation points (avoid out of bounds)
        nq = 100
        nres = 6 ! number of results 
        allocate(xq(nq), yq(nq), vqan(nq, nres), dv(nq, nres), &
            vq(nq, nres), temp(nq))
        call random_number(xq)
        call random_number(yq)

        !nq = 2
        !nres = 6 ! number of results 
        !allocate(xq(nq), yq(nq), vqan(nq, nres), dv(nq, nres), vq(nq, nres))
        !xq = [0.12, 0.999]
        !yq = [0.16, 0.999]
        
        xq = xq*Lx
        yq = yq*Ly

        ! Linear test case
        !=================
        ! Build interpolant values
        allocate(a(3)) ! coefficients
        a = [0, 1, 3]
        do i = 1, nx+1
            do j = 1, ny+1
                v(i, j) = a(1) + a(2)*xgv(i) + a(3)*ygv(j)
            end do 
        end do

        ! Build interpolant
        meth = 'uniformgrid'
        call interp%SetParameters(meth, 3, 6)
        call interp%Construct(xgv, ygv, v)

        ! Test1: evaluate the interpolant at some query points and 
        ! evaluate the analytic solution
        vqan(:, 1) = a(1) + a(2)*xq + a(3)*yq ! field value
        vqan(:, 2) = a(3) ! dfdy
        vqan(:, 3) = 0 ! d2fdy2
        vqan(:, 4) = a(2) !dfdx
        vqan(:, 5:6) = 0 ! d2fdx2, d2fdxdy

        k = 1
        do i = 0, 2
            do j = 0, 2-i
                call interp%Evaluate(xq, yq, i, j, temp)
                vq(:, k) = temp
                k = k + 1
            end do 
        end do 

        ! Evaluate norm
        dv = vqan - vq
        allocate(relerr(nres), reshuffle(nres))
        do k = 1, nres
            relerr(k) = maxval(abs(dv(:, k)/vqan(:, k)))
        enddo 
        reshuffle = [1, 4, 2, 6, 5, 3]

        ! Print
        print *, 'test case     f     dfdx    dfdy    d2fdx2      d2fdxdy     d2fdy2      '
        print *, 'Linear field test case', relerr(reshuffle)

        ! Housekeeping
        deallocate(a)
        call interp%deallocate()

        ! Quadratic case
        !===============
        ! Build interpolant values
        allocate(a(6)) ! coefficients
        a = [0.0, 1.0, 3.0, 0.01, 0.05, 0.0016]
        do i = 1, nx+1
            do j = 1, ny+1
                v(i, j) = a(1) + a(2)*xgv(i) + a(3)*ygv(j) + &
                    a(4)*xgv(i)*ygv(j) + a(5)*xgv(i)**2 + a(6)*ygv(j)**2
            end do 
        end do

        ! Build interpolant
        meth = 'uniformgrid'
        call interp%SetParameters(meth, 3, 6)
        call interp%Construct(xgv, ygv, v)

        ! Test1: evaluate the interpolant at some query points and 
        ! evaluate the analytic solution
        vqan(:, 1) = a(1) + a(2)*xq + a(3)*yq + a(4)*xq*yq + &
            a(5)*xq**2 + a(6)*yq**2 ! field value
        vqan(:, 2) = a(3) + a(4)*xq + 2*a(6)*yq ! dfdy
        vqan(:, 3) = 2*a(6) ! d2fdy2
        vqan(:, 4) = a(2) + a(4)*yq + 2*a(5)*xq !dfdx
        vqan(:, 5) = a(4) !d2fdxdy
        vqan(:, 6) = 2*a(5) ! d2fdx2, d2fdxdy

        k = 1
        do i = 0, 2
            do j = 0, 2-i
                call interp%Evaluate(xq, yq, i, j, temp)
                vq(:, k) = temp
                k = k + 1
            end do 
        end do 

        ! Evaluate norm
        dv = vqan - vq
        do k = 1, nres
            relerr(k) = maxval(abs(dv(:, k)/vqan(:, k)))
        enddo 

        ! Print
        print *, 'Quadratic field test case', relerr(reshuffle)

        ! Housekeeping
        deallocate(a)
        call interp%deallocate()

        ! Housekeeping
        !=============
        deallocate(xgv, ygv, xq, yq, vqan, dv, temp, v)
        call DisplayTestEnd()



    end subroutine 

    

end module 


program Tests

    ! Description
    !============
    ! Tests for the structured 2D interpolant. We test for a given 
    ! analytical field whether it is correctly represented again by 
    ! the structured polynomial interpolant

    use GOAT_tests 

    call RunAllTests()


end program Tests


                     