!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module serves as an interface between dedicated linear solvers
! (e.g. UMFPACK, BLAS/LAPACK, ...), using our own datatypes. 

module mod_linearsolverinterface
    
    ! Initialize
    !============
    ! Load modules
    use mod_precision

    ! Binding with C
    use, intrinsic :: iso_c_binding 

    ! The usual
    implicit none
    save
    private 

    external umf4def

    ! Public routines
    public TestUMFPACK


    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    ! UMFPACK solvers
    !interface 
    !    subroutine LinSolve_UMFPACK()

    !    end subroutine
    !end interface

    ! UMFPACK auxiliaries (di)
    !=========================
    interface 
        ! Defaults
        subroutine UmfpackDefaultsDI(control) & 
            bind(c, name='umfpack_di_defaults')
            
            ! C module
            use, intrinsic :: iso_c_binding 

            ! Declare
            real(c_double)      :: control(20)

        end subroutine

        ! Symbolic creator
        subroutine UmfpackSymbolicDI(n, m, Ap, Ai, Ax, symbolic, &
            control, info) & 
            bind(c, name='umfpack_di_symbolic')
            
            ! C module
            use, intrinsic :: iso_c_binding 

            ! Declare
            integer(c_int), value           :: n, m
            integer(c_int), intent(in)      :: Ap(*), Ai(*)
            real(c_double)                  :: control(20), info(90)
            real(c_double), intent(in)      :: Ax(*)
            type(c_ptr)                     :: symbolic 

        end subroutine

        ! Symbolic destructor
        subroutine UmfpackFreeSymbolicDI(symbolic) & 
            bind(c, name='umfpack_di_free_symbolic')
            
            ! C module
            use, intrinsic :: iso_c_binding 

            ! Declare
            type(c_ptr)                     :: symbolic 

        end subroutine

        ! Numeric creator
        subroutine UmfpackNumericDI(Ap, Ai, Ax, symbolic, numeric, &
            control, info) &
            bind(c, name='umfpack_di_numeric')

            ! C module
            use, intrinsic :: iso_c_binding 

            ! Declare
            integer(c_int)                  :: Ap(*), Ai(*)
            real(c_double)                  :: control(20), info(90)
            real(c_double)                  :: Ax(*)
            type(c_ptr)                     :: numeric
            type(c_ptr), value              :: symbolic

        end subroutine

        ! Numeric destructor
        subroutine UmfpackFreeNumericDI(numeric) & 
            bind(c, name='umfpack_di_free_numeric')
            
            ! C module
            use, intrinsic :: iso_c_binding 

            ! Declare
            type(c_ptr)                     :: numeric 

        end subroutine

        ! Solver
        subroutine UmfpackSolveDI(sys, Ap, Ai, Ax, x, b, numeric, &
            control, info) &
            bind(c, name='umfpack_di_solve')

            ! C module
            use, intrinsic :: iso_c_binding 

            ! Declare
            integer(c_int)                  :: Ap(*), Ai(*)
            integer(c_int), value           :: sys
            real(c_double)                  :: control(20), info(90)
            real(c_double)                  :: Ax(*), x(*), b(*)
            type(c_ptr), value              :: numeric

        end subroutine

        ! Info reporter
        subroutine UmfpackReportInfoDI(control, info) &
            bind(c, name='umfpack_di_report_info')

            ! C module
            use, intrinsic :: iso_c_binding 

            ! Declare
            real(c_double)                     :: info(90), control(20)
            
        end subroutine


    end interface    

    contains

    !==================================================================!
    !                                                                  !
    !                             ROUTINES                             !
    !                                                                  !
    !==================================================================!

    ! Simple test routine
    subroutine TestUMFPACK()

        ! Description
        !============
        ! Simple test routine to see if umfpack is installed properly. 
        ! Creates a matrix and computes the solution with a fixed rhs
        ! for which the solution is known. This example is highly based
        ! on the example of the UMFPACK manual by Tim Davis. 

        use, intrinsic :: iso_c_binding

        ! Declare variables
        !==================
        integer(c_int)                      :: n, nval, sys
        integer(c_int), allocatable         :: Ap(:), Ai(:)
        real(c_double), allocatable         :: Ax(:), b(:), sol(:)
        real(c_double)                      :: control(20), info(90)
        type(c_ptr)                         :: symbolic, numeric

        ! Allocate & initialize
        !======================
        ! Number of equations
        n = 5

        ! Number of nonzero values
        nval = 12

        ! Set the system type to be solved
        sys = 0 ! standard Ax = b

        ! Allocate
        allocate(Ap(n+1))
        allocate(Ai(nval))
        allocate(Ax(nval))
        allocate(b(n))
        allocate(sol(n))

        ! Construct matrix representation
        Ap = [0, 2, 5, 9, 10, nval] ! column pointer
        Ai = [0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4] ! row index starting from 0
        Ax = [2, 3, 3, -1, 4, 4, -3, 1, 2, 2, 6, 1] ! values 
        b = [8, 45, -3, 3, 19]

        ! Solve
        !======
        ! Solve using the C-routines of UMFPACK. The Double precision, 
        ! integer routines (_di_) are used 

        ! Set defaults
        call UmfpackDefaultsDI(control)

        ! Set symbolic factorization
        call UmfpackSymbolicDI(n, n, Ap, Ai, Ax, symbolic, &
            control, info)
            print *, info(1)

        ! Compute numerical factorization
        call UmfpackNumericDI(Ap, Ai, Ax, symbolic, numeric, &
            control, info)
            print *, info(1)

        ! Compute the solution
        call UmfpackSolveDI(sys, Ap, Ai, Ax, sol, b, numeric, &
            control, info)
            print *, info(1)
        
        ! Print the solution
        print *, sol

        ! Destroy symbolic
        call UmfpackFreeSymbolicDI(symbolic)

        ! Destroy numeric
        call UmfpackFreeNumericDI(numeric)

        ! Deallocate
        deallocate(Ap)
        deallocate(Ai)
        deallocate(Ax)
        deallocate(b)
        deallocate(sol)

    end subroutine

    ! The actual decent solver


end module