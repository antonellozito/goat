!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains several hessian approximation types that are 
! used in the (S)QP modules. The following formats are present:
! - dense:      classic dense representation, suitable for small 
!               problems. Can easily employ BFGS updating etc. 
! - sparse:     sparse representation, suitable for large problems with
!               many nonzeros in the matrix. Not suitable for BFGS 
!               updating! Useful to have simple representation for 
!               e.g. steepest descent or purely gradient-based 
!               procedures. 

! The types all derive from the abstract hessian approximation type. If 
! no format-specific routines/data is needed, the overarching 
! 'HessianApproximationUDT' type is recommended. 

! TODO: need to add limited memory representation... 

module optmod_hessianapproximation

    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_constants
    use mod_sparseinterface
    use mod_linearsolverinterface

    ! The usual
    implicit none
    save
    private 
    
    ! Set public
    public ConstructHessianApproximation
    public HessianApproximationUDT

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    ! Hessian approximation class (abstract)
    type, abstract :: HessianApproximationUDT

        integer(I8)                 :: nphi
        character(:), allocatable   :: updatemethod 

    contains 

        ! Hessian updating
        procedure(UpdateHessianBFGSINT), deferred :: UpdateHessianBFGS
        generic :: UpdateHessian        => UpdateHessianBFGS

        ! Inverse Hessian vector product
        procedure(InverseHessianVectorProductINT), deferred     :: &
            InverseHessianVectorProduct

    end type

    ! Dense hessian representation
    type, extends(HessianApproximationUDT)  :: DenseHessianApproximationUDT 

        ! Description
        !============
        ! Dense hessian storage format. Values are simply stored in the
        ! 'val' array (nphi-by-nphi). Not memory efficient if the 
        ! Hessian is likely to be sparse.

        ! Values
        real(R8), allocatable       :: val(:, :)

    contains 

        ! Hessian updating
        procedure :: UpdateHessianBFGS  => UpdateDenseHessianBFGS

        ! Inverse Hessian vector product
        procedure :: InverseHessianVectorProduct &   
            => InverseDenseHessianVectorProduct

    end type

    ! Sparse hessian representation
    type, extends(HessianApproximationUDT)  :: SparseHessianApproximationUDT 

        ! Description
        !============
        ! Sparse hessian storage format. Values are stored as sparse
        ! matrix in the 'val' field. Useful when the hessian is either
        ! known exactly, or when updating the matrix preserves the 
        ! sparse matrix structure. Can be used for e.g. (scaled) 
        ! steepest descent methods. 

        ! Values
        type(MySparseUDT)           :: val 

    contains 

        ! Hessian updating
        procedure :: UpdateHessianBFGS  => UpdateSparseHessianBFGS

        ! Inverse Hessian vector product
        procedure :: InverseHessianVectorProduct &   
            => InverseSparseHessianVectorProduct

    end type


    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    ! Constructor interface
    interface ConstructHessianApproximation
        module procedure ConstructSparseHessianApproximation
        module procedure ConstructDenseHessianApproximation
        module procedure ConstructDiagonalHessianApproximation
        module procedure ConstructElementwiseHessianApproximation
    end interface

    ! Abstract interface
    abstract interface 

            ! Hessian updating with BFGS
            subroutine UpdateHessianBFGSINT(hess, phi0, phi1, grad0, grad1)
                import :: R8, HessianApproximationUDT
                class(HessianApproximationUDT)      :: hess
                real(R8), intent(in), dimension(:)  :: phi0, phi1, &
                    grad0, grad1
            end subroutine

            ! Product between inverse hessian and vector
            function InverseHessianVectorProductINT(hess, x) result(y)
                import :: R8, HessianApproximationUDT 
                class(HessianApproximationUDT)      :: hess 
                real(R8), intent(in), dimension(:)  :: x 
                real(R8), allocatable               :: y(:)
            end function

    end interface

    contains

    !==================================================================!
    !                                                                  !
    !                              ROUTINES                            !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                            CONSTRUCTORS                          !
    !------------------------------------------------------------------!

    ! General dense constructor
    function ConstructDenseHessianApproximation(updatemethod, &
        nphi, inithess) result(hess)

        ! Description
        !============
        ! Construct a dense hessian approximation but return as a 
        ! general hessian approximation object.

        ! Declare variables
        !==================
        ! Arguments
        character(*), intent(in)                        :: updatemethod 
        integer(I8), intent(in)                         :: nphi
        real(R8), intent(in)                            :: inithess(1:nphi, 1:nphi)
        class(HessianApproximationUDT), allocatable     :: hess

        ! Auxiliary

        ! Initialize
        !===========
        ! Allocate
        allocate(DenseHessianApproximationUDT::hess)

        ! Set parameters
        !===============
        ! General parameters 
        hess%updatemethod   = updatemethod
        hess%nphi           = nphi

        ! Type-specific parameters
        select type (hess)

        type is (DenseHessianApproximationUDT)

            ! Initialize to initial matrix
            hess%val = inithess

        class default 
        
            ! Nothing to do here, shouldn't be possible
            call gdErrorHandler('ConstructDenseHessianApproximation: implementation bug ')

        end select

    end function

    ! General sparse constructor
    function ConstructSparseHessianApproximation(updatemethod, &
        nphi, inithess) result(hess)

        ! Description
        !============
        ! Construct a dense hessian approximation but return as a 
        ! general hessian approximation object.

        ! Declare variables
        !==================
        ! Arguments
        character(*), intent(in)                        :: updatemethod 
        integer(I8), intent(in)                         :: nphi
        type(MySparseUDT), intent(in)                   :: inithess
        class(HessianApproximationUDT), allocatable     :: hess

        ! Auxiliary

        ! Initialize
        !===========
        ! Allocate
        allocate(SparseHessianApproximationUDT::hess)

        ! Set parameters
        !===============
        ! General parameters 
        hess%updatemethod   = updatemethod
        hess%nphi           = nphi

        ! Check
        if (.not. allocated(inithess%val)) then 
            call gdErrorHandler('ConstructSparseHessianApproximation: ' // &
                'initial hessian not allocated')
        end if 
        if ((inithess%nrow /= nphi) .or. (inithess%ncol /= nphi)) then 
            call gdErrorHandler('ConstructSparseHessianApproximation: ' // &
                'inconsistent dimensions of initial hessian')
        end if 

        ! Type-specific parameters
        select type (hess)

        type is (SparseHessianApproximationUDT)

            ! Set initial hessian
            hess%val = inithess

        class default 
        
            ! Nothing to do here, shouldn't be possible
            call gdErrorHandler('ConstructDenseHessianApproximation: implementation bug ')

        end select

    end function

    ! Diagonal constructor
    function ConstructDiagonalHessianApproximation(updatemethod, &
        nphi, diagind, diagvals, hessiantype) result(hess)

        ! Description
        !============
        ! Construct a diagonal hessian approximation based on the given
        ! diagonal indices (0: main diagonal, offsets go columnwise) and
        ! values for those diagonals. The values ('diagvals') can be either
        ! the same size as the diagonal indices, in which case each 
        ! diagonal has a uniform value, or the same as the number of 
        ! elements to be attributed (for non-uniform, general diagonal
        ! matrices). 'hessiantype' should specify which kind of hessian
        ! (dense, sparse, limited memory, ... ) type should be taken. 

        ! Declare variables
        !==================
        ! Arguments
        character(*), intent(in)            :: updatemethod, hessiantype 
        integer(I8), intent(in)             :: nphi, diagind(:)
        real(R8), intent(in)                :: diagvals(:)
        class(HessianApproximationUDT), allocatable :: hess

        ! Auxiliary
        integer(I8)                         :: nd, nv, nel 
        integer(I8), allocatable            :: row(:), col(:)
        real(R8), allocatable               :: eld(:, :), val(:)
        logical                             :: uniformdiagonals = .false.
        type(MySparseUDT)                   :: diagmat, elsp

        ! Loop
        integer(I8)                         :: i, k, cc

        ! Initialize
        !===========
        ! Hedge for empty hessian & check input type
        select case ('hessiantype')

        case ('dense')

            ! Check nphi
            if (nphi <= 0) then 
                allocate(eld(0, 0))
                eld = 0
                hess = ConstructDenseHessianApproximation(updatemethod, nphi, eld)
                return 
            end if 

        case ('sparse')

            ! Check nphi
            if (nphi <= 0) then 
                elsp = SpZeros(nphi, nphi)
                hess = ConstructSparseHessianApproximation(updatemethod, nphi, elsp)
                return 
            end if 

        case default 

            call gdErrorHandler('ConstructDiagonalHessianApproximation: ' // & 
                'hessian type: "' // hessiantype // '" not implemented')

        end select

        ! Check input format
        nel = 0
        nd = size(diagind)
        nv = size(diagvals)
        if (nd == nv) then 
            ! Uniform values for each diagonal
            uniformdiagonals = .true. 
        else
            ! Non-uniform diagonals, values per element
            uniformdiagonals = .false. 
        end if 

        ! Compute total number of elements
        do i = 1, nd 
            nel = nel + (nphi - abs(diagind(i)))
        end do

        ! Check dimensions
        if (.not. uniformdiagonals .and. (nel /= nv)) then 
            call gdErrorHandler('ConstructDiagonalHessianApproximation: ' // &
                'mismatch between specified diagonal indices and ' // & 
                'number of values. Check if given values correspond to ' // &
                'total number of elements in the hessian or to the ' // &
                'number of specified diagonals')
        end if 

        ! Construct initial sparse matrix
        !================================
        ! Construct row, col indices
        allocate(row(nel), col(nel), val(nel))
        cc = 0
        do i = 1, nd
            if (diagind(i) >= 0) then 
                row(cc+1:cc+(nphi-abs(diagind(i)))) = [(k, k = 1, nphi-abs(diagind(i)))]
                col(cc+1:cc+(nphi-abs(diagind(i)))) = [(k, k = abs(diagind(i)), nphi)]
            else
                col(cc+1:cc+(nphi-abs(diagind(i)))) = [(k, k = 1, nphi-abs(diagind(i)))]
                row(cc+1:cc+(nphi-abs(diagind(i)))) = [(k, k = abs(diagind(i)), nphi)]
            end if 
            cc = cc + (nphi-abs(diagind(i)))
        end do

        ! Construct values
        if (uniformdiagonals) then 
            ! Assign by looping
            cc = 0
            do i = 1, nd 
                val(cc+1:cc+(nphi-abs(diagind(i)))) = diagvals(i)
                cc = cc + (nphi-abs(diagind(i)))
            end do
        else
            val = diagvals
        end if 
        
        ! Construct hessian approximation
        !================================
        hess = ConstructElementwiseHessianApproximation(updatemethod, &
            nphi, row, col, val, hessiantype)

    end function

    ! Elementwise constructor
    function ConstructElementwiseHessianApproximation(updatemethod, &
        nphi, row, col, val, hessiantype) result(hess)
    end function

    !------------------------------------------------------------------!
    !                             OPERATORS                            !
    !------------------------------------------------------------------!

    ! Hessian approximation updating (dense)
    subroutine UpdateDenseHessianBFGS(hess, phi0, phi1,  grad0, grad1)

        ! Description
        !============
        ! Update the hessian using the BFGS formula with Powell's trick 
        ! to avoid non-positive definite hessians. The step is computed 
        ! as sk = phi1 - phi0, and yk = grad1 - grad0. Note that we 
        ! assume that the initial hessian was chosen to be positive 
        ! definite! Otherwise, this will not work. 

        ! Declare variables
        !==================
        ! Arguments
        class(DenseHessianApproximationUDT)         :: hess 
        real(R8), intent(in), dimension(:)          :: phi0, phi1, &
            grad0, grad1 

        ! Auxiliary
        real(R8)                                    :: sbs(1, 1), ys(1, 1), theta(1, 1)
        real(R8), allocatable, dimension(:, :)      :: sk, yk, bs
        real(R8), allocatable, dimension(:, :)      :: Bk, ykt

        ! Compute preliminaries
        !======================
        ! Allocate
        allocate(sk(hess%nphi, 1), yk(hess%nphi, 1), bs(hess%nphi, 1))

        ! Compute differences
        sk(:, 1) = phi1 - phi0 
        yk(:, 1) = grad1 - grad0 

        ! Compute products
        sbs = matmul(matmul(transpose(sk), hess%val), sk)
        ys  =matmul(transpose(yk), sk)
        bs = matmul(Bk, sk)

        ! Compute theta
        if (ys(1, 1) < 0.2*sbs(1, 1)) then 
            theta = (0.2*sbs - ys)/(sbs - ys)
        else
            theta = 0
        end if 

        ! Determine modified yk, ykt
        allocate(ykt(size(hess%val, 1), 1)) ! 2D array for outer product later on
        ykt(:, 1) = yk(:, 1) + theta(1, 1)*(matmul(hess%val, sk(:, 1)) - yk(:, 1))

        ! BFGS Update
        !============
        ! Store Bk
        Bk = hess%val 

        ! Compute new hessian
        hess%val = Bk - matmul(bs, transpose(bs))/sbs(1, 1) &
            + matmul(ykt, transpose(ykt))/ys(1, 1)

        ! Ensure symmetry
        hess%val = 0.5*(hess%val + transpose(hess%val))


    end subroutine

    ! Hessian approximation updating (sparse)
    subroutine UpdateSparseHessianBFGS(hess, phi0, phi1,  grad0, grad1)

        ! Description
        !============
        ! Update the hessian using the BFGS formula with Powell's trick 
        ! to avoid non-positive definite hessians. The step is computed 
        ! as sk = phi1 - phi0, and yk = grad1 - grad0. Note that we 
        ! assume that the initial hessian was chosen to be positive 
        ! definite! Otherwise, this will not work. 

        ! Note: though we provide the implementation here, it is 
        ! absolutely not recommended to be used! The update will likely 
        ! be dense and destroy the sparse structure of the initial
        ! hessian... We throw a warning since this is most likely 
        ! unintended usage

        ! Declare variables
        !==================
        ! Arguments
        class(SparseHessianApproximationUDT)        :: hess 
        real(R8), intent(in), dimension(:)          :: phi0, phi1, &
            grad0, grad1 

        ! Auxiliary
        real(R8)                                    :: sbs(1, 1), ys(1, 1), theta(1, 1)
        real(R8), allocatable, dimension(:, :)         :: sk, yk 
        real(R8), allocatable, dimension(:, :)      :: ykt, bs

        type(MysparseUDT)                           :: Bk, dval

        ! Throw warning
        !==============
        print *, 'UpdateSparseHessianBFGS: updating with BFGS will lead' // &
            'to dense matrix! Consider updating the sparse hessian differently'

        ! Compute preliminaries
        !======================
        ! Allocate
        allocate(sk(hess%nphi, 1), yk(hess%nphi, 1), bs(hess%nphi, 1))

        ! Compute differences
        sk(:, 1) = phi1 - phi0 
        yk(:, 1) = grad1 - grad0 

        ! Compute products
        sbs = matmul(transpose(sk), reshape(hess%val%MatrixVectorProduct(sk(:, 1)), [hess%nphi, 1]))
        ys  = matmul(transpose(yk), sk)

        ! Compute theta
        if (ys(1, 1) < 0.2*sbs(1, 1)) then 
            theta = (0.2*sbs - ys)/(sbs - ys)
        else
            theta = 0
        end if 

        ! Determine modified yk, ykt
        allocate(ykt(hess%nphi, 1)) ! 2D array for outer product later on
        ykt(:, 1) = yk(:, 1) + theta(1, 1)*(hess%val%MatrixVectorProduct(sk(:, 1)) - yk(:, 1))

        ! BFGS Update
        !============
        ! Store Bk
        Bk = hess%val 

        ! Store bs
        bs(:, 1) = hess%val%MatrixVectorProduct(sk(:, 1))

        ! Compute new hessian - Actually dense due to outer products!
        dval = ConstructMySparse(-matmul(bs, transpose(bs))/sbs(1, 1) &
            - matmul(ykt, transpose(ykt))/ys(1, 1))
        hess%val = Bk + dval


    end subroutine

    ! Inverse hessian vector product (dense)
    function InverseDenseHessianVectorProduct(hess, x) result(y)

        ! Description
        !============
        ! Compute y = B^-1 x. Of course we don't compute the inverse
        ! directly, but we call the dense solver routines.

        ! Declare variables
        !==================
        ! Arguments
        class(DenseHessianApproximationUDT)         :: hess 
        real(R8), intent(in), dimension(:)          :: x 
        real(R8), allocatable, dimension(:)         :: y 

        ! Auxiliary
        integer(I8)                                 :: flag

        ! Checks
        !=======
        ! Ensure correct dimensions
        if (.not. (hess%nphi == size(x, 1))) then 
            call gdErrorHandler('InverseDenseHessianVectorProduct:  ' // &
                'right hand side does not have compatible dimensions')
        end if 

        ! Compute solution
        !=================
        ! Initialize
        y = x 

        ! Compute
        call SolveDenseLinearSystemDI(hess%val, x, y, flag)

        ! Check convergence
        if (flag /= 0) then 
            print *, 'InverseDenseHessianVectorProduct: linear solver did not converge'
        end if 

    end function

    ! Inverse hessian vector product (sparse)
    function InverseSparseHessianVectorProduct(hess, x) result(y)

        ! Description
        !============
        ! Compute y = B^-1 x. Of course we don't compute the inverse
        ! directly, but we call the sparse solver routines.

        ! Declare variables
        !==================
        ! Arguments
        class(SparseHessianApproximationUDT)        :: hess 
        real(R8), intent(in), dimension(:)          :: x 
        real(R8), allocatable, dimension(:)         :: y 

        ! Auxiliary
        integer(I8)                                 :: flag

        ! Checks
        !=======
        ! Ensure correct dimensions
        if (.not. (hess%nphi == size(x, 1))) then 
            call gdErrorHandler('InverseDenseHessianVectorProduct:  ' // &
                'right hand side does not have compatible dimensions')
        end if 

        ! Compute solution
        !=================
        ! Initialize
        y = x 

        ! Compute
        call SolveSparseLinearSystemDI(hess%val, x, y, flag)

        ! Check convergence
        if (flag /= 0) then 
            print *, 'InverseSparseHessianVectorProduct: linear solver did not converge'
        end if 

    end function

    
  

end module