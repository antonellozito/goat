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
    use mod_sparseinterface

    ! Binding with C
    use, intrinsic :: iso_c_binding 

    ! The usual
    implicit none
    save
    private 

    external umf4def

    ! Public routines
    public TestUMFPACK !  tester
    public SolveSparseLinearSystemDI ! Sparse system solver
    public SolveDenseLinearSystemDI ! dense system solver

    ! UMFPACK variables - hard coded here...
    integer(c_int) :: umfpack_a = 0
    integer(c_int) :: umfpack_control = 20
    integer(c_int) :: umfpack_info = 90
    ! integer(c_int) :: umfpack_prl = 1 ! printing level

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

    ! UMFPACK double integer routines
    !================================
    ! Note: the final authority of the implementation on the c-side of
    ! these routines lies with the manual provided by Dr. Tim Davis. 
    ! See the github repository 
    ! (https://github.com/DrTimothyAldenDavis/SuiteSparse, and the 
    ! UMFPACK package thereof) and the manual included there. Below,
    ! we provide a brief description in the interface of the expected
    ! arguments for each function. 
    
    interface 
        ! Defaults
        subroutine UmfpackDefaultsDI(control) & 
            bind(c, name='umfpack_di_defaults')

            ! Description
            !============
            ! This routine simply loads in the default values for the 
            ! control variables used. Normally, there is no need to 
            ! change these. If this is the case, refer to the manual. 
            
            ! C module
            use, intrinsic :: iso_c_binding 

            ! Declare
            real(c_double)      :: control(*)

        end subroutine

        ! Symbolic creator
        subroutine UmfpackSymbolicDI(nrow, ncol, Ap, Ai, Ax, symbolic, &
            control, info) & 
            bind(c, name='umfpack_di_symbolic')

            ! Description
            !============
            ! This routine loads in the symbolic factorization, which is
            ! a sort of preprocessing step for the numerical 
            ! factorization. The arguments to be parsed are the number 
            ! of rows and columns of the matrix, the sparse matrix in 
            ! compressed column storage format (Ap is the column 
            ! pointer, Ai the row indices, and Ax the values), and the
            ! control and info arrays. The 'symbolic' argument here 
            ! should be a pointer, which, on exit, points to the 
            ! symbolic object used in the numerical counterpart of this
            ! routine. 
            
            ! C module
            use, intrinsic :: iso_c_binding 

            ! Import
            import umfpack_control, umfpack_info

            ! Declare
            integer(c_int), value           :: nrow, ncol
            integer(c_int), intent(in)      :: Ap(*), Ai(*)
            real(c_double)                  :: &
                control(umfpack_control), info(umfpack_info)
            real(c_double), intent(in)      :: Ax(*)
            type(c_ptr)                     :: symbolic 

        end subroutine

        ! Symbolic destructor
        subroutine UmfpackFreeSymbolicDI(symbolic) & 
            bind(c, name='umfpack_di_free_symbolic')

            ! Description
            !============
            ! This routine destroys the symbolic object
            
            ! C module
            use, intrinsic :: iso_c_binding 

            ! Declare
            type(c_ptr)                     :: symbolic 

        end subroutine

        ! Numeric creator
        subroutine UmfpackNumericDI(Ap, Ai, Ax, symbolic, numeric, &
            control, info) &
            bind(c, name='umfpack_di_numeric')

            ! Description
            !============
            ! Constructs the object to factorize the matrix numerically,
            ! which is used in the solver. The input arguments are 
            ! largely the same as in the symbolic counterpart of this 
            ! routine, though the symbolic argument is now passed by
            ! value, whereas the numeric argument is passed as a 
            ! pointer. 

            ! C module
            use, intrinsic :: iso_c_binding 

            ! Declare
            integer(c_int)                  :: Ap(*), Ai(*)
            real(c_double)                  :: control(*), info(*)
            real(c_double)                  :: Ax(*)
            type(c_ptr)                     :: numeric
            type(c_ptr), value              :: symbolic

        end subroutine

        ! Numeric destructor
        subroutine UmfpackFreeNumericDI(numeric) & 
            bind(c, name='umfpack_di_free_numeric')

            ! Description
            !============
            ! Destroys the numeric object
            
            ! C module
            use, intrinsic :: iso_c_binding 

            ! Declare
            type(c_ptr)                     :: numeric 

        end subroutine

        ! Solver
        subroutine UmfpackSolveDI(sys, Ap, Ai, Ax, x, b, numeric, &
            control, info) &
            bind(c, name='umfpack_di_solve')

            ! Description
            !============
            ! General solver for a system ( Ax = b). Note that the 'A'
            ! in this system may be manipulated by adjusting the 'sys'
            ! argument. To solve a regular system, this should equal 
            ! UMFPACK_A (which is currently 0). This UMFPACK_A is a 
            ! variable that should be supplied by UMFPACK.

            ! C module
            use, intrinsic :: iso_c_binding 

            ! Declare
            integer(c_int)                  :: Ap(*), Ai(*)
            integer(c_int), value           :: sys
            real(c_double)                  :: control(*), info(*)
            real(c_double)                  :: Ax(*), x(*), b(*)
            type(c_ptr), value              :: numeric

        end subroutine

        ! Info reporter - not used
        subroutine UmfpackReportInfoDI(control, info) &
            bind(c, name='umfpack_di_report_info')

            ! C module
            use, intrinsic :: iso_c_binding 

            ! Declare
            real(c_double)                     :: info(*), control(*)
            
        end subroutine

        ! Coordinate to compressed column storage format
        subroutine UmfpackTriplet2ColDI(nrow, ncol, nz, Ti, Tj, Tx, &
            Ap, Ai, Ax, Map) &
            bind(c, name='umfpack_di_triplet_to_col')

            ! C module
            use, intrinsic :: iso_c_binding 

            ! Declare
            integer(c_int), value           :: nrow, ncol, nz
            integer(c_int), intent(in)      :: Ti(*), Tj(*)
            real(c_double), intent(in)      :: Tx(*)
            integer(c_int)                  :: Ap(*), Ai(*)
            real(c_double)                  :: Ax(*)
            integer(c_int)                  :: Map(*)


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
        real(c_double)                      :: &
            control(umfpack_control), info(umfpack_info)
        type(c_ptr)                         :: symbolic, numeric

        ! Allocate & initialize
        !======================
        ! Print
        print *, 'umfpack_a: ', umfpack_a 
        print *, 'umfpack_control: ', umfpack_control
        print *, 'umfpack_info: ', umfpack_info

        ! Number of equations
        n = 5

        ! Number of nonzero values
        nval = 12

        ! Set the system type to be solved
        sys = umfpack_a ! standard Ax = b

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

    ! The sparse solver 
    subroutine SolveSparseLinearSystemDI(A, b, sol, flag)

        ! Description
        !============
        ! This routine solves a sparse linear system by calling the 
        ! UMFPACK solvers. Here, it is assumed that the variables are
        ! in double precision. The input matrix A should be in a 
        ! MySparse type format and will be converted to compressed 
        ! column storage format (CCS) using the UMFPACK routines. The
        ! system is then solved by calling the required UMFPACK solver. 
        ! The rhs should be given in b, the solution will be returned in
        ! sol. Note that the solver here is a direct solver! 

        ! The result whether the solver is successful is parsed through
        ! 'info_solver'. It is zero when successful, and nonzero when
        ! not.

        ! Declare variables
        !==================
        ! Arguments
        type(MySparseUDT)               :: A 
        real(R8), intent(in)            :: b(:)
        real(R8), allocatable           :: sol(:)

        ! Auxiliary 
        integer(c_int), allocatable     :: Ap(:), Ai(:), Map(:)
        integer(I8)                     :: flag
        real(c_double), allocatable     :: Ax(:)

        real(c_double), allocatable         :: sol_c(:)
        real(c_double)                      :: &
            control(umfpack_control), info(umfpack_info)

        type(c_ptr)                         :: symbolic, numeric

        ! Loop

        ! Data

        ! Convert to CCS
        !===============
        ! First, convert the row and column indices to zero based 
        ! indexing
        A%col = A%col - 1
        A%row = A%row - 1

        ! Allocate
        allocate(Ap(A%ncol+1))
        allocate(Ai(A%nval))
        allocate(Ax(A%nval))
        allocate(Map(A%nval))
        allocate(sol_c(A%nrow))

        ! Call converter
        call UmfpackTriplet2ColDI(A%nrow, A%ncol, A%nval, A%row, &
            A%col, A%val, Ap, Ai, Ax, Map)

        ! Solve
        !======
        ! Solve using the C-routines of UMFPACK. The Double precision, 
        ! integer routines (_di_) are used 

        ! Set defaults
        call UmfpackDefaultsDI(control)

        ! Set symbolic factorization
        call UmfpackSymbolicDI(A%nrow, A%ncol, Ap, Ai, Ax, symbolic, &
            control, info)

        ! Compute numerical factorization
        call UmfpackNumericDI(Ap, Ai, Ax, symbolic, numeric, &
            control, info)

        ! Compute the solution
        call UmfpackSolveDI(umfpack_a, Ap, Ai, Ax, sol_c, b, numeric, &
            control, info)
        if (info(1) .ne. 0) then 
            print *, 'Linear solver may not have converged, info: ', info(1)
        end if

        ! Destroy symbolic
        call UmfpackFreeSymbolicDI(symbolic)

        ! Destroy numeric
        call UmfpackFreeNumericDI(numeric)

        ! Post-process
        !=============
        ! Cast solution
        sol = sol_c 

        ! Cast flag
        flag = nint(info(1))

        ! Reconstruct A
        A%col = A%col + 1
        A%row = A%row + 1

        ! Deallocate
        deallocate(Ap, Ai, Ax, Map, sol_c)

    end subroutine

    ! The dense solver
    subroutine SolveDenseLinearSystemDI(A, b, sol, flag)

        ! Description
        !============
        ! Solve (small) dense linear systems. We rely on LAPACK to 
        ! compute the solution

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)    :: A(:, :)
        real(R8), intent(in)    :: b(:)
        real(R8), intent(out)   :: sol(size(b))
        integer(I8)             :: flag 

        ! Auxiliary
        integer(I8), allocatable        :: ipiv(:)
        integer                         :: neq, info
        double precision, allocatable   :: rhs(:), lhs(:, :)

        ! Initialize
        !===========
        ! Check dimensions
        neq = size(b)
        if ( (size(A, 1) /= neq) .or. (size(A, 2) /= neq)) then 
            ! Incompatible sizes
            call gdErrorHandler('SolveDenseLinearSystemDI: incompatible size of A w.r.t. b')
        end if 

        ! Set solver data
        allocate(rhs(neq), lhs(neq, neq))
        rhs = b
        lhs = A
        flag = 0

        ! Call the solver (make sure input is in right format!)
        allocate(ipiv(neq))
        call dgesv(neq, 1, lhs, neq, ipiv, rhs, neq, info)

        ! Check if converged
        if (info .ne. 0) then
            ! Not converged 
            flag = 1
            print *, 'dgesv could not converge, info: ', info 
        end if

        ! Set solution
        sol = rhs

    end subroutine



end module