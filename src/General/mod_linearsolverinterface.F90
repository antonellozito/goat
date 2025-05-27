!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module serves as an interface between dedicated linear solvers
! (e.g. UMFPACK, BLAS/LAPACK, MUMPS, ...), using our own datatypes. 
! It is assumed that at least UMFPACK and BLAS/LAPACK are available. 
! The (d)mumps implementation is only defined if the compiler flag 
! USE_DMUMPS is defined. 

! IMPORTANT: when using MUMPS, the main program must start with an
! MPI_ÏNIT call!

! Note: we provide on the one hand routines that are directly callable
! to solve certain systems of linear equations (sparse and dense). 
! However, these routines are directly based on UMFPACK and LAPACK (resp.).
! A more general approach (yet possibly less elegant) is to initialize a 
! solver object that may or may not persist. In this way, it is easier 
! to optimize the solver (e.g. options can be set by the user, which is
! currently not possible for the callable routines) and to avoid 
! reinitializing data. Best practice is to use the LinearSolverUDT 
! abstract user defined type and initialize it using the public routine
! InitializeLinearSolver with the appropriate linear solver options. This
! will initialize the appropriate derived type and deal with proper 
! initialization/destruction automatically. 

! Note: the solvers provided here are only double precision solvers.

module mod_linearsolverinterface
    
    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_sparseinterface
    use mod_errorhandler

    ! Binding with C
    use, intrinsic :: iso_c_binding 

#if (defined(MUMPS) || defined(USE_MPI))
    use mpi
#endif 
    ! The usual
    implicit none
    save
    private 

#ifdef MUMPS 
    ! Load MUMPS MPI header and structure
    include 'dmumps_struc.h'
#endif 

    

    external umf4def
    integer, external :: ilaenv

    ! Public routines
    !================
    ! Standalone solvers
    public SolveSparseLinearSystemDI ! Sparse system solver (centralized - UMFPACK)
    public SolveDenseLinearSystemDI ! dense system solver
    public ConstructDLinearSolver ! Linear solver object constructor 
#ifdef MUMPS 
    public SolveSparseLinearSystemDIDMUMPS ! sparse system solver (centralized - MUMPS)
#endif

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

    ! Abstract solver type
    type, abstract, public :: DLinearSolverUDT 

        ! Description
        !============
        ! This type should be used as template for other solver types.
        ! Each solver should have a solving method for sparse and dense
        ! linear systems. Each
        ! solution should also return an error status flag that is zero
        ! when the solver was successful, and is otherwise non-zero. 
        ! This structure also contains solver options, such as verbosity
        ! level, reusing factorizations, and so on. This ensures 
        ! abstraction at the user level of the underlying solver 
        ! specific option set (e.g. the umfpack or mumps control 
        ! variables). It is the responsibility of the developer to  
        ! ensure proper behavior using these options for each linear
        ! solver. Some options may of course not be supported/unavailable
        ! for some solvers. 

        ! Options (public)
        !-----------------
        integer(I8)         :: verbosity ! for output printing control

        ! Internal variables (private)
        !-----------------------------
        ! Underlying solver type
        character(:), allocatable   :: solvername 

        ! Matrix and matrix inverse storage (dense and sparse)
        real(R8), allocatable, dimension(:, :), private :: Ad, Adinv
        type(MySparseUDT), private  :: Asp, Aspinv  

    contains 

        ! Initializer (solver specific, only called by constructor)
        procedure(InitializeLinearSystemSolverDINT), deferred ::  &
            Initialize

        ! Linear solver methods
        procedure(SolveSparseLinearSystemDINT), deferred :: SolveSparseLinearSystem 
        procedure(SolveDenseLinearSystemDINT), deferred  :: SolveDenseLinearSystem

        ! Solver finalization routine
        procedure(FinalizeLinearSystemSolverDINT), deferred :: &
            Finalize 

    end type

    ! LAPACK solver (direct)
    type, extends(DLinearSolverUDT), public  :: DLAPACKLinearSolverUDT

        ! Description
        !============
        ! LAPACK dense system solver. Mainly suited for dense 
        ! systems, sparse system solver will throw a warning and then 
        ! convert sparse system into dense matrix format to solve. 

    contains 

        ! Initialization
        procedure :: Initialize                 => InitializeDLAPACKSolver

        ! Solver procedures
        procedure :: SolveSparseLinearSystem    => SolveSparseLinearSystemDLAPACK
        procedure :: SolveDenseLinearSystem     => SolveDenseLinearSystemDLAPACK 

        ! Finalization
        procedure :: Finalize                   => FinalizeDLAPACKSolver
        final :: DestroyDLAPACKSolver

    end type

    ! UMFPACK solver (direct)
    type, extends(DLinearSolverUDT), public  :: DUMFPACKLinearSolverUDT

        ! Description
        !============
        ! Umfpack solver from suitesparse. Mainly suited for sparse 
        ! systems, dense system solver will throw a warning and then 
        ! convert dense system into sparse matrix format to solve. 

    contains 

        ! Initialization
        procedure :: Initialize                 => InitializeDUMFPACKSolver

        ! Solver procedures
        procedure :: SolveSparseLinearSystem    => SolveSparseLinearSystemDUMFPACK
        procedure :: SolveDenseLinearSystem     => SolveDenseLinearSystemDUMFPACK 

        ! Finalization
        procedure :: Finalize                   => FinalizeDUMFPACKSolver
        final :: DestroyDUMFPACKSolver


    end type

#ifdef MUMPS
    ! DMUMPS solver (direct)
    type, extends(DLinearSolverUDT), public  :: DMUMPSLinearSolverUDT

        ! Description
        !============
        ! DMUMPS solver, mainly suited for large sparse systems. 
        ! depending on how the external mumps library is compiled, 
        ! either an MPI, openMP, or fully sequential version is 
        ! available - none of which is visible here, since we only use 
        ! the overarching mumps_type object to call the solver. For more
        ! information, visit the mumps pages (https://mumps-tech.com or 
        ! https://mumps-solver.org)
            
        ! Mumps type
        type(dmumps_struc) :: mumps_par 

    contains 

        ! Initialization
        procedure :: Initialize                 => InitializeDMUMPSSolver

        ! Solver procedures
        procedure :: SolveSparseLinearSystem    => SolveSparseLinearSystemDMUMPS
        procedure :: SolveDenseLinearSystem     => SolveDenseLinearSystemDMUMPS 

        ! Finalization for derived types
        procedure :: Finalize                   => FinalizeDMUMPSSolver
        final :: DestroyDMUMPSSolver

    end type
#endif 

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    ! Solver objects
    !===============
    abstract interface 

        ! Solver initialization 
        subroutine InitializeLinearSystemSolverDINT(ls)
            import :: DLinearSolverUDT
            class(DLinearSolverUDT)             :: ls 
        end subroutine

        ! Sparse system solver
        subroutine SolveSparseLinearSystemDINT(ls, A, b, sol, flag)
            import :: MySparseUDT, R8, I8, DLinearSolverUDT
            class(DLinearSolverUDT)             :: ls 
            type(MySparseUDT), intent(in)       :: A 
            real(R8), dimension(:), intent(in)  :: b
            real(R8), allocatable, dimension(:), intent(out)    :: sol 
            integer(I8), intent(out)            :: flag  
        end subroutine

        ! Dense system solver
        subroutine SolveDenseLinearSystemDINT(ls, A, b, sol, flag)
            import :: R8, I8, DLinearSolverUDT
            class(DLinearSolverUDT)                 :: ls 
            real(R8), dimension(:, :), intent(in)   :: A 
            real(R8), dimension(:), intent(in)      :: b
            real(R8), allocatable, dimension(:), intent(out)    :: sol 
            integer(I8), intent(out)                :: flag  
        end subroutine

        ! Solver finalization 
        subroutine FinalizeLinearSystemSolverDINT(ls)
            import :: DLinearSolverUDT
            class(DLinearSolverUDT)             :: ls 
        end subroutine

    end interface 

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
    
    ! Others
    !=======
    ! Interface for dense system solver for ease
    interface SolveDenseLinearSystemDI
        module procedure SolveDenseLinearSystemDI1D 
        module procedure SolveDenseLinearSystemDI2D
    end interface

    contains

    !==================================================================!
    !                                                                  !
    !                             ROUTINES                             !
    !                                                                  !
    !==================================================================!

    ! General solver
    !===============
    ! Constructor
    function ConstructDLinearSolver(solvername) result(ls)

        ! Description
        !============
        ! This routine constructs a linear solver object with default
        ! settings. These default settings should allow to solve any
        ! type of matrix, but are not optimized for any specific problem.
        ! For this, modify the available public options in the 
        ! linear solver object directly, or call a dedicated 
        ! initializer for a certain type of matrix. 

        ! Declare variables
        !==================
        ! Arguments
        character(*), intent(in)                :: solvername 
        class(DLinearSolverUDT), allocatable    :: ls 

        ! Initialize
        !===========
        ! Check allocation status
        if (allocated(ls)) then
            deallocate(ls)
        end if 

        ! Initialize 
        select case(solvername)

        case ('UMFPACK', 'umfpack', 'Umfpack')

            allocate(DUMFPACKLinearSolverUDT::ls)

#ifdef MUMPS 
        case ('MUMPS', 'mumps', 'Mumps')

            allocate(DMUMPSLinearSolverUDT::ls)

#endif 

        case ('LAPACK', 'lapack', 'Lapack')

            allocate(DLAPACKLinearSolverUDT::ls)

        case default 

            call gdErrorHandler('ConstructDLinearSolver: unknown option: ' // & 
                solvername)

        end select

        ! Set default options and solver type
        ls%solvername   = solvername 
        ls%verbosity    = 1_I8 

        ! Call initializer
        call ls%Initialize()

    end function 

    ! UMFPACK
    !========
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

        ! Check
        !======
        ! Empty matrix?
        if ((a%nrow == 0) .or. (a%ncol == 0)) then
            if (allocated(sol)) then  
                deallocate(sol)
            end if 
            allocate(sol(0))
            flag = 0
            return 
        end if

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

    ! Initialization
    subroutine InitializeDUMFPACKSolver(ls)

        ! Description
        !============
        ! Initialization of umfpack solver - currently nothing to be 
        ! done, but may change in the future. 

        ! Declare variables
        !==================
        ! Arguments
        class(DUMFPACKLinearSolverUDT)       :: ls 

    end subroutine

    ! Finalization
    subroutine FinalizeDUMFPACKSolver(ls)

        ! Description
        !============
        ! Finalization of umfpack solver - currently nothing to be 
        ! done, but may change in the future. 

        ! Declare variables
        !==================
        ! Arguments
        class(DUMFPACKLinearSolverUDT)       :: ls 

    end subroutine

    ! Destruction
    subroutine DestroyDUMFPACKSolver(ls)

        ! Description
        !============
        ! Destruction of umfpack solver. Calls finalization routine

        ! Declare variables
        !==================
        ! Arguments
        type(DUMFPACKLinearSolverUDT)   :: ls 

        ! Finalize
        !=========
        call ls%Finalize()

        ! Print
        if (ls%verbosity > 1) then 
            print *, 'Destroyed umfpack solver'
        end if 

    end subroutine

    ! Sparse solver, object based
    subroutine SolveSparseLinearSystemDUMFPACK(ls, A, b, sol, flag)
        
        ! Description
        !============
        ! For now just a wrapper for SolveSparseLinearSystemDI, but 
        ! may have a more optimized implementation in the future, e.g.
        ! avoiding refactorization in some cases

        ! Declare variables
        !==================
        ! Arguments
        class(DUMFPACKLinearSolverUDT)      :: ls 
        type(MySparseUDT), intent(in)       :: A 
        real(R8), dimension(:), intent(in)  :: b
        real(R8), allocatable, dimension(:), intent(out)    :: sol 
        integer(I8), intent(out)        :: flag 

        ! Auxiliary

        ! Call subroutine
        !================
        call SolveSparseLinearSystemDI(A, b, sol, flag)

    end subroutine

    ! Dense solver, object based
    subroutine SolveDenseLinearSystemDUMFPACK(ls, A, b, sol, flag)

        ! Description
        !============
        ! 'dense' solver for UMFPACK, but actually we just convert the
        ! dense matrix to a sparse matrix format and solve using 
        ! the sparse solver. This is often not desireable, so we print 
        ! a message if this happens. This message can be suppressed by
        ! setting the verbosity level of the solver lower than 1. 

        ! Declare variables
        !==================
        ! Arguments
        class(DUMFPACKLinearSolverUDT)          :: ls 
        real(R8), dimension(:, :), intent(in)   :: A 
        real(R8), dimension(:), intent(in)      :: b
        real(R8), allocatable, dimension(:), intent(out)    :: sol 
        integer(I8), intent(out)        :: flag 

        ! Auxiliary
        type(MySparseUDT)  :: Asparse

        ! Call subroutine
        !================
        ! Print warning
        if (ls%verbosity > 0) then 
            print *, 'SolveDenseLinearSystemDUMFPACK: dense solver of umfpack ' // & 
                'calls sparse solver, may not be desired. Set solver ' // & 
                'verbosity < 1 to suppress this message'
        end if 

        ! Convert
        Asparse = ConstructMySparse(A)

        ! Solve
        call SolveSparseLinearSystemDI(Asparse, b, sol, flag)

    end subroutine

    ! LAPACK
    !=======
    ! The dense solver
    subroutine SolveDenseLinearSystemDI1D(A, b, sol, flag, Ainv)

        ! Description
        !============
        ! Solve (small) dense linear systems. We rely on LAPACK to 
        ! compute the solution. If desired, the inverse is computed

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)    :: A(:, :)
        real(R8), intent(in)    :: b(:)
        real(R8), intent(out)   :: sol(size(b))
        integer(I8)             :: flag 

        ! Optional arguments
        real(R8), intent(out), optional     :: Ainv(size(b), size(b)) ! inv

        ! Auxiliary
        integer(I8), allocatable        :: ipiv(:)
        integer                         :: neq, info, infoinv, nb, lda
        double precision, allocatable   :: rhs(:), lhs(:, :), work(:)

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
        lda = max(neq, 1)
        call dgesv(neq, 1, lhs, lda, ipiv, rhs, lda, info)

        ! Check if converged
        if (info .ne. 0) then
            ! Not converged 
            flag = 1
            print *, 'dgesv could not converge, info: ', info 
        end if

        ! Set solution
        sol = rhs

        ! Compute inverse if desired
        if (present(Ainv)) then 
            ! Compute preliminaries 
            nb = ilaenv( 1, 'DGETRI', ' ', neq, -1, -1, -1 )
            allocate(work(nb*neq))

            ! Compute inverse
            Ainv = lhs
            call dgetri(neq, Ainv, neq, ipiv, work, neq*nb, infoinv)

            ! Check 
            if (infoinv /= 0) then 
                print *, 'dgetri could not determine inverse, info: ', infoinv 
            end if
        end if 

    end subroutine

    subroutine SolveDenseLinearSystemDI2D(A, b, sol, flag, Ainv)

        ! Description
        !============
        ! Solve (small) dense linear systems. We rely on LAPACK to 
        ! compute the solution. If desired, the inverse is computed

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)    :: A(:, :)
        real(R8), intent(in)    :: b(:, :)
        real(R8), intent(out)   :: sol(size(b, 1), size(b, 2))
        integer(I8)             :: flag 

        ! Optional arguments
        real(R8), intent(out), optional     :: Ainv(size(b, 1), size(b, 1)) ! inv

        ! Auxiliary
        integer(I8), allocatable        :: ipiv(:)
        integer                         :: neq, info, infoinv, nb, lda
        double precision, allocatable   :: rhs(:, :), lhs(:, :), work(:)

        ! Initialize
        !===========
        ! Check dimensions
        neq = size(b, 1)
        if ( (size(A, 1) /= neq) .or. (size(A, 2) /= neq)) then 
            ! Incompatible sizes
            call gdErrorHandler('SolveDenseLinearSystemDI: incompatible size of A w.r.t. b')
        end if 

        ! Set solver data
        allocate(rhs(neq, size(b, 2)), lhs(neq, neq))
        rhs = b
        lhs = A
        flag = 0

        ! Call the solver (make sure input is in right format!)
        allocate(ipiv(neq))
        lda = max(neq, 1)
        call dgesv(neq, size(b, 2), lhs, lda, ipiv, rhs, lda, info)

        ! Check if converged
        if (info .ne. 0) then
            ! Not converged 
            flag = 1
            print *, 'dgesv could not converge, info: ', info 
        end if

        ! Set solution
        sol = rhs

        ! Compute inverse if desired
        if (present(Ainv)) then 
            ! Compute preliminaries 
            nb = ilaenv( 1, 'DGETRI', ' ', neq, -1, -1, -1 )
            allocate(work(nb*neq))

            ! Compute inverse
            Ainv = lhs
            call dgetri(neq, Ainv, neq, ipiv, work, neq*nb, infoinv)

            ! Check 
            if (infoinv /= 0) then 
                print *, 'dgetri could not determine inverse, info: ', infoinv 
            end if
        end if 

    end subroutine

    ! Initialization
    subroutine InitializeDLAPACKSolver(ls)

        ! Description
        !============
        ! Initialization of lapack solver - nothing to be done here 

        ! Declare variables
        !==================
        ! Arguments
        class(DLAPACKLinearSolverUDT)       :: ls 

    end subroutine

    ! Finalization
    subroutine FinalizeDLAPACKSolver(ls)

        ! Description
        !============
        ! Finalization of lapack solver - nothing to be done here 

        ! Declare variables
        !==================
        ! Arguments
        class(DLAPACKLinearSolverUDT)       :: ls 

    end subroutine

    ! Destruction
    subroutine DestroyDLAPACKSolver(ls)

        ! Description
        !============
        ! Destruction of LAPACK solver. Calls finalization routine

        ! Declare variables
        !==================
        ! Arguments
        type(DLAPACKLinearSolverUDT)   :: ls 

        ! Finalize
        !=========
        call ls%Finalize()

        ! Print
        if (ls%verbosity > 1) then 
            print *, 'Destroyed lapack solver'
        end if 

    end subroutine

    ! Sparse solver, object based
    subroutine SolveSparseLinearSystemDLAPACK(ls, A, b, sol, flag)

        ! Description
        !============
        ! 'sparse' system solver for LAPACK, but in fact this just 
        ! converts the sparse matrix into a dense matrix and then 
        ! solves it. This is likely unintended behavior, so a warning 
        ! is thrown if verbosity is > 0.

        ! Declare variables
        !==================
        ! Arguments
        class(DLAPACKLinearSolverUDT)       :: ls 
        type(MySparseUDT), intent(in)       :: A 
        real(R8), dimension(:), intent(in)  :: b
        real(R8), allocatable, dimension(:), intent(out)    :: sol 
        integer(I8), intent(out)        :: flag 

        ! Auxiliary
        real(R8), allocatable, dimension(:, :)  :: Afull

        ! Call subroutine
        !================
        ! Print warning
        if (ls%verbosity > 0) then 
            print *, 'SolveSparseLinearSystemDLAPACK: sparse solver of lapack ' // & 
                'calls dense solver, may not be desired. Set solver ' // & 
                'verbosity < 1 to suppress this message'
        end if 

        ! Convert
        call A%ConvertToFull(Afull)

        ! Check
        if (.not. allocated(sol)) then 
            allocate(sol(size(b)))
        elseif (size(sol) /= size(b)) then 
            deallocate(sol)
            allocate(sol(size(b)))
        end if 

        ! Solve
        call SolveDenseLinearSystemDI1D(Afull, b, sol, flag)

    end subroutine 

    ! Dense solver, object based
    subroutine SolveDenseLinearSystemDLAPACK(ls, A, b, sol, flag)

        ! Description
        !============
        ! For now just a wrapper for SolveDenseLinearSystemDI1D, but 
        ! may have a more optimized implementation in the future, e.g.
        ! avoiding refactorization in some cases

        ! Declare variables
        !==================
        ! Arguments
        class(DLAPACKLinearSolverUDT)           :: ls 
        real(R8), dimension(:, :), intent(in)   :: A 
        real(R8), dimension(:), intent(in)      :: b
        real(R8), allocatable, dimension(:), intent(out)    :: sol 
        integer(I8), intent(out)        :: flag 

        ! Call subroutine
        !================
        ! Check
        if (.not. allocated(sol)) then 
            allocate(sol(size(b)))
        elseif (size(sol) /= size(b)) then 
            deallocate(sol)
            allocate(sol(size(b)))
        end if 

        ! Call solver
        call SolveDenseLinearSystemDI1D(A, b, sol, flag)

    end subroutine


#ifdef MUMPS
    ! DMUMPS 
    !=======
    ! Standalone solver
    subroutine SolveSparseLinearSystemDIDMUMPS(A, b, sol, flag)

        ! Declare variables
        !==================
        ! Arguments
        type(MySparseUDT), intent(in)       :: A 
        real(R8), dimension(:), intent(in)  :: b
        real(R8), allocatable, dimension(:), intent(out)    :: sol 
        integer(I8), intent(out)        :: flag 

        ! Auxiliary
        type(dmumps_struc)                  :: mumps_par
        integer                             :: verbosity

        ! Check
        !======
        ! Empty matrix?
        verbosity = 0
        if ((a%nrow == 0) .or. (a%ncol == 0)) then
            if (allocated(sol)) then  
                deallocate(sol)
            end if 
            allocate(sol(0))
            flag = 0
            return 
        end if
        
        ! Check dimensions
        if (A%ncol /= A%nrow) then 
            call gdErrorHandler('SolveSparseLinearSystemDMUMPS: matrix ' // & 
                'is not square, not supported')
        end if 
        if (A%nrow /= size(b)) then 
            call gdErrorHandler('SolveSparseLinearSystemDMUMPS: incompatible ' // & 
                'dimensions between rhs and lhs')
        end if 


        ! Initialize
        !===========
        ! To be done for all processes
        ! Define the communicator
        mumps_par%comm = mpi_comm_world 

        ! Initialize the instance
        mumps_par%job = -1 ! flag for initialization
        mumps_par%sym = 0 ! general unsymmetric matrix
        mumps_par%par = 1 ! let the main process also do work
        call dmumps(mumps_par)

        ! Adjust some default parameters
        mumps_par%icntl(14) = 50
        mumps_par%icntl(4) = 1 ! only error message output
        mumps_par%icntl(16) = 4

        ! Print a message
        if (verbosity > 0) then 
            print *, 'InitializeDMUMPSSolver: solver initialized'
        end if 

        ! Define problem
        !===============
        ! Define problem on the host (process 0)
        if (mumps_par%myid == 0) then 
            ! Initialize
            mumps_par%icntl(14) = 50
            ! Set dimensions
            mumps_par%n = A%nrow  
            mumps_par%nnz = A%nval  

            ! Allocate
            allocate(mumps_par%irn(mumps_par%nnz), mumps_par%jcn(mumps_par%nnz), &
                mumps_par%A(mumps_par%nnz), mumps_par%rhs(mumps_par%n))
            
            ! Assign
            mumps_par%irn    = A%row 
            mumps_par%jcn    = A%col 
            mumps_par%a      = A%val
            mumps_par%rhs    = b 

        end if 

        ! Solve
        !======
        ! Call package for solution
        mumps_par%job = 1 ! Analyze
        call dmumps(mumps_par)
        mumps_par%job = 2 ! solve
        mumps_par%ICNTL(14) = 1000 
        call dmumps(mumps_par)

        mumps_par%job = 3 ! solve
        call dmumps(mumps_par)

        ! Check if successful
        if (mumps_par%infog(1) < 0) then 
            print *, 'SolveSparseLinearSystemDMUMPS: linear solver returned ' // & 
                'with error flags: ', mumps_par%infog(1), mumps_par%infog(2) 
        end if 

        ! Solution has been assembled on the host
        if (mumps_par%myid == 0) then 
            sol = mumps_par%rhs 
        end if 

        ! Deallocate data
        if (mumps_par%myid == 0) then 
            deallocate(mumps_par%irn, mumps_par%jcn, &
                mumps_par%A, mumps_par%rhs)
        end if 

        ! Destroy the instance (if not done already)
        if (mumps_par%job /= -2) then 
            mumps_par%job = -2 
            call dmumps(mumps_par)
        end if 
        
        ! Check if successful
        if (mumps_par%infog(1) < 0) then 
            print *, 'FinalizeDMUMPSSolver: finalization failed and returned with ' // &
                'error flags = ', mumps_par%infog(1), mumps_par%infog(2)
            call gdErrorHandler('FinalizeDMUMPSSolver: could not finalize, exiting...')
        end if 


    end subroutine

    ! Initialization
    subroutine InitializeDMUMPSSolver(ls)

        ! Description
        !============
        ! Initialization of MUMPS solver. Here, only the general
        ! MUMPS data is assigned and default values for solving 
        ! a general system of equations are used. 

        ! Declare variables
        !==================
        ! Arguments
        class(DMUMPSLinearSolverUDT)       :: ls 

        ! Initialize
        !===========
        ! To be done for all processes
        ! Define the communicator
        ls%mumps_par%comm = mpi_comm_world 
        print *, mpi_comm_world

        ! Initialize the instance
        ls%mumps_par%job = -1 ! flag for initialization
        ls%mumps_par%sym = 0 ! general unsymmetric matrix
        ls%mumps_par%par = 1 ! let the main process also do work
        call dmumps(ls%mumps_par)

        ! Adjust some default parameters
        ls%mumps_par%icntl(14) = 35

        ! Print a message
        if (ls%verbosity > 0) then 
            print *, 'InitializeDMUMPSSolver: solver initialized'
        end if 

    end subroutine

    ! Finalization
    subroutine FinalizeDMUMPSSolver(ls)

        ! Description
        !============
        ! Finalize the solver by deallocating any remaining allocatables
        ! specific to the mumps solver and by calling the mumps instance
        ! destructor. 

        ! Declare variables
        !==================
        ! Arguments
        class(DMUMPSLinearSolverUDT)    :: ls 

        ! Deallocate
        !===========
        ! Check for each process if the allocatable arrays irn, jcn, 
        ! a, rhs are allocated. If so, deallocate
        if (associated(ls%mumps_par%irn)) deallocate(ls%mumps_par%irn)
        if (associated(ls%mumps_par%jcn)) deallocate(ls%mumps_par%jcn)
        if (associated(ls%mumps_par%a)) deallocate(ls%mumps_par%a)
        if (associated(ls%mumps_par%rhs)) deallocate(ls%mumps_par%rhs)

        ! Destroy the instance (if not done already)
        if (ls%mumps_par%job /= -2) then 
            ls%mumps_par%job = -2 
            call dmumps(ls%mumps_par)
        end if 
        
        ! Check if successful
        if (ls%mumps_par%infog(1) < 0) then 
            print *, 'FinalizeDMUMPSSolver: finalization failed and returned with ' // &
                'error flags = ', ls%mumps_par%infog(1), ls%mumps_par%infog(2)
            call gdErrorHandler('FinalizeDMUMPSSolver: could not finalize, exiting...')
        end if 

    end subroutine

    ! Destruction
    subroutine DestroyDMUMPSSolver(ls)

        ! Description
        !============
        ! Destruction of MUMPS solver. Calls finalization routine

        ! Declare variables
        !==================
        ! Arguments
        type(DMUMPSLinearSolverUDT)   :: ls 

        ! Finalize
        !=========
        call ls%Finalize()

        ! Print
        if (ls%verbosity > 1) then 
            print *, 'Destroyed mumps solver'
        end if 

    end subroutine

    ! Sparse solver, object based
    subroutine SolveSparseLinearSystemDMUMPS(ls, A, b, sol, flag)

        ! Description
        !============
        ! This routine solves a sparse system of linear equations. For
        ! now, we assume that the input is given by the host process
        ! and hence we also assemble the problem on the host process. 

        ! Note: optimizations are not yet in place, but may be added in 
        ! the future. 

        ! Declare variables
        !==================
        ! Arguments
        class(DMUMPSLinearSolverUDT)        :: ls 
        type(MySparseUDT), intent(in)       :: A 
        real(R8), dimension(:), intent(in)  :: b 
        real(R8), dimension(:), allocatable, intent(out)    :: sol 
        integer(I8), intent(out)            :: flag 

        ! Auxiliary

        ! Check
        !======
        ! Empty matrix?
        if ((a%nrow == 0) .or. (a%ncol == 0)) then
            if (allocated(sol)) then  
                deallocate(sol)
            end if 
            allocate(sol(0))
            flag = 0
            return 
        end if
        
        ! Check dimensions
        if (A%ncol /= A%nrow) then 
            call gdErrorHandler('SolveSparseLinearSystemDMUMPS: matrix ' // & 
                'is not square, not supported')
        end if 
        if (A%nrow /= size(b)) then 
            call gdErrorHandler('SolveSparseLinearSystemDMUMPS: incompatible ' // & 
                'dimensions between rhs and lhs')
        end if 

        ! Initialize
        !===========
        ! Define problem on the host (process 0)
        if (ls%mumps_par%myid == 0) then 
            ! Initialize
            ! Set dimensions
            ls%mumps_par%n = A%nrow  
            ls%mumps_par%nnz = A%nval  


            ! Allocate
            allocate(ls%mumps_par%irn(ls%mumps_par%nnz), ls%mumps_par%jcn(ls%mumps_par%nnz), &
                ls%mumps_par%A(ls%mumps_par%nnz), ls%mumps_par%rhs(ls%mumps_par%n))
            
            ! Assign
            ls%mumps_par%irn    = A%row 
            ls%mumps_par%jcn    = A%col 
            ls%mumps_par%a      = A%val
            ls%mumps_par%rhs    = b 

        end if 

        ! Call package for solution
        ls%mumps_par%job = 1 ! Analyze
        call dmumps(ls%mumps_par)
        ! Set memory larger
        ls%mumps_par%ICNTL(14) = 1000 ! systems are often tough to solve 
        ls%mumps_par%job = 2 ! factorize
        call dmumps(ls%mumps_par)

        ls%mumps_par%job = 3 ! solve
        call dmumps(ls%mumps_par)

        ! Check if successful
        if (ls%mumps_par%infog(1) < 0) then 
            print *, 'SolveSparseLinearSystemDMUMPS: linear solver returned ' // & 
                'with error flags: ', ls%mumps_par%infog(1), ls%mumps_par%infog(2) 
        end if 

        ! Solution has been assembled on the host
        if (ls%mumps_par%myid == 0) then 
            sol = ls%mumps_par%rhs 
        end if 

        ! Deallocate data
        if (ls%mumps_par%myid == 0) then 
            deallocate(ls%mumps_par%irn, ls%mumps_par%jcn, &
                ls%mumps_par%A, ls%mumps_par%rhs)
        end if 

    end subroutine

    ! Dense solver, object based
    subroutine SolveDenseLinearSystemDMUMPS(ls, A, b, sol, flag)

        ! Description
        !============
        ! This routine serves simply as a wrapper to call the sparse
        ! implementation since there is no dense alternative for 
        ! MUMPS. If the verbosity level is > 0, then a message will be
        ! shown, since this is likely undesired behavior. 

        ! Declare variables
        !==================
        ! Arguments
        class(DMUMPSLinearSolverUDT)            :: ls 
        real(R8), dimension(:, :), intent(in)   :: A 
        real(R8), dimension(:), intent(in)      :: b 
        real(R8), dimension(:), allocatable, intent(out)    :: sol 
        integer(I8), intent(out)            :: flag 

        ! Auxiliary
        type(MySparseUDT)  :: Asparse

        ! Call subroutine
        !================
        ! Print warning
        if (ls%verbosity > 0) then 
            print *, 'SolveDenseLinearSystemDMUMPS: dense solver of mumps ' // & 
                'calls sparse solver, may not be desired. Set solver ' // & 
                'verbosity < 1 to suppress this message'
        end if 

        ! Convert
        Asparse = ConstructMySparse(A)

        ! Solve
        call ls%SolveSparseLinearSystem(Asparse, b, sol, flag)



    end subroutine

#endif 

end module