!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module serves as an interface between 'our' sparse matrix 
! representation (which only serves as a temporary fix) and any other 
! decent sparse library. Basically, an intermediate derived type is 
! provided that simply stores the row and column indices of the 
! non-zero elements along with their values, and the dimensions nrow and
! ncol of the matrix. This should be sufficient to construct the sparse
! matrix in any desireable format. 

! Notes
!======
! Note 1: tt is likely that more efficient routines for operators, 
! conversion and so on exist. Currently, our implementation relies
! largely on home-brewn functions, which may, however, be interchanged
! or interfaced in the future with more optimized libraries. 

module mod_sparseinterface
    
    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use Clayer

    ! The usual
    implicit none
    save
    public 

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    ! Abstract types
    !===============
    ! Sparse matrix type
    type :: MySparseUDT

        ! Description
        !============
        ! Our own sparse matrix format. 

        ! Fields:
        ! - nrow:           number of rows of the matrix 
        ! - ncol:           number of columns of the matrix 
        ! - nval:           number of non-zero values of the matrix 
        ! - row:            row indices (nval-by-1)
        ! - col:            column indices (nval-by-1)
        ! - val:            values (nval-by-1)

        integer(I8)                 :: nrow = 0, ncol = 0, nval = 0
        integer(I8), allocatable    :: row(:), col(:)
        real(R8), allocatable       :: val(:)

    contains

        ! Constructor
        procedure :: Initialize     => InitializeMySparse

        ! Sparse to full transformation routine
        procedure :: Full           => ConvertToFull

        ! Extraction routines
        procedure :: ExtractColumnFull 
        procedure :: ExtractRowFull

        ! Deletion routines
        procedure :: DeleteRowsLogical 
        procedure :: DeleteRowsInteger 
        generic   :: DeleteRows     => DeleteRowsLogical, DeleteRowsInteger

        procedure :: DeleteColumnsLogical 
        procedure :: DeleteColumnsInteger 
        generic   :: Deletecolumns  => DeleteColumnsLogical, DeleteColumnsInteger

        ! Summation routines
        procedure :: SumColumnwiseFull
        procedure :: SumRowwiseFull

        ! Multiplication (other than overloaded element-wise functions)
        procedure :: MatrixVectorProduct

        ! Transposition
        procedure :: Transpose      => Transp

        ! Concatenation
        procedure :: Concatenate2
        procedure :: ConcatenateN 
        generic   :: Concatenate    => Concatenate2, ConcatenateN

        ! Conversion
        procedure :: ConvertToFull

        ! Housekeeping procedures
        procedure :: Allocate       => AllocateMySparse
        procedure :: Deallocate     => DeallocateMySparse

        ! Destructor
        final :: DestroyMySparse

    end type

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    ! Overload the summation operator
    interface  operator(+)
        module procedure AddSparse
    end interface

    ! Overload the multiplication operator
    interface operator(*)
        module procedure MultiplyWithScalarLeft
        module procedure MultiplyWithScalarRight
        module procedure MultiplyWithVectorRowwiseLeft
        module procedure MultiplyWithVectorRowwiseRight
        module procedure MultiplySparseMatrices
    end interface

    contains

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    ! Constructor
    subroutine InitializeMySparse(mysparse, nrow, ncol, nval)

        ! Description
        !============
        ! Sparse matrix constructor, including allocation

        ! Initialize
        !===========
        ! The usual
        implicit none

        ! Declare variables
        !==================
        ! Arguments
        class(MySparseUDT)           :: mysparse 
        integer(I8), intent(in)     :: nrow, ncol, nval 

        ! Construct
        !==========
        ! Set dimensions
        mysparse%nrow = nrow 
        mysparse%ncol = ncol
        mysparse%nval = nval

        ! Allocate
        call mysparse%Allocate()

    end subroutine

    !------------------------------------------------------------------!
    !                         DATA CONVERSION                          !
    !------------------------------------------------------------------!

    ! Sparse to full conversion
    subroutine ConvertToFull(mysparse, fullmatrix)

        ! Description
        !============
        ! Converts a sparse matrix to a full matrix, which is a 2D array
        ! with the same dimensions as the sparse matrix. 

        ! Initialize
        !===========
        ! The usual
        implicit none

        ! Declare variables
        !==================
        ! Arguments
        class(MySparseUDT)                      :: mysparse 
        real(R8), allocatable, intent(inout)    :: fullmatrix(:, :)

        ! Loop
        integer(I16)                             :: i

        ! Auxiliary 
        real(R8)                                :: tv 
        integer(I8)                             :: tc, tr

        ! Construct
        !==========
        ! Allocate
        if (.not. allocated(fullmatrix)) then 
            allocate(fullmatrix(mysparse%nrow, mysparse%ncol))
        end if

        ! Initialize
        fullmatrix(:, :) = 0
        
        ! Associate
        associate(&
            col     => mysparse%col, &
            row     => mysparse%row, &
            val     => mysparse%val, &
            nval    => mysparse%nval)

        ! Loop over all elements and add
        do i = 1, nval
            ! Unpack
            tc = col(i)
            tr = row(i)
            tv = val(i)

            ! print *, i, tc, tr, tv
            fullmatrix(tr, tc) = tv
            
            
        end do

        ! Deassoicate
        end associate



    end subroutine

    ! Column extraction
    subroutine ExtractColumnFull(mysparse, col, colID)

        ! Description
        !============
        ! Extract the column of a sparse matrix in full form, without
        ! constructing a full version of the sparse matrix. colID should
        ! contain the index of the column. 

        ! Declare variables
        !==================
        ! Arguments
        class(MySparseUDT)              :: mysparse 
        real(R8), allocatable           :: col(:)
        integer(I8)                     :: colID

        ! Loop
        integer(I8)                     :: i

        ! Extract
        !========
        ! Allocate
        if (.not. allocated(col)) then 
            allocate(col(mysparse%nrow))
        end if

        ! Initialize
        col(:) = 0

        ! Set
        do i = 1, mysparse%nval 
            if (mysparse%col(i) == colID) then 
                col(mysparse%row(i)) = col(mysparse%row(i)) + &
                    mysparse%val(i)
            end if
        end do

    end subroutine

    ! Row extraction
    subroutine ExtractRowFull(mysparse, row, rowID)

        ! Description
        !============
        ! Extract the column of a sparse matrix in full form, without
        ! constructing a full version of the sparse matrix. colID should
        ! contain the index of the column. 

        ! Declare variables
        !==================
        ! Arguments
        class(MySparseUDT)              :: mysparse 
        real(R8), allocatable           :: row(:)
        integer(I8)                     :: rowID

        ! Loop
        integer(I8)                     :: i

        ! Extract
        !========
        ! Allocate
        if (.not. allocated(row)) then 
            allocate(row(mysparse%ncol))
        end if 

        ! Initialize
        row(:) = 0

        ! Set
        do i = 1, mysparse%nval 
            if (mysparse%row(i) == rowID) then 
                row(mysparse%col(i)) = row(mysparse%col(i)) + &
                    mysparse%val(i)
            end if
        end do

    end subroutine

    ! Conversion from CSparse to MySparse format
    subroutine ConvertToMySparse(mysparsecs, mysparse)

        ! Description
        !============
        ! Convert sparse matrix in CSparseUDT format to MySparseUDT 
        ! format. It is assumed that the matrix is allocated. Note that
        ! MySparseUDT expects rows and columns to start from one! 

        ! Declare
        !========
        ! Arguments
        type(MySparseUDT)               :: mysparse 
        type(CSparseUDT)       :: mysparsecs 

        ! Auxiliary
        real(c_double), pointer         :: valp(:)
        integer(c_int), pointer         :: rowp(:), colp(:)
        type(MySparseUDT)               :: temp

        ! Convert
        !========
        ! Associate
        associate( &
            nrow        => mysparsecs%nrow, &
            ncol        => mysparsecs%ncol, &
            nval        => mysparsecs%nval, &
            row         => mysparsecs%row, &
            col         => mysparsecs%col, &
            val         => mysparsecs%val &
        )

        ! Initialize & allocate
        temp%nval = nval 
        temp%nrow = nrow 
        temp%ncol = ncol 
        call temp%Allocate()
        
        ! Extract row, col, val data from C pointers
        call c_f_pointer(val, valp, [nval])
        temp%val = valp 
        call c_f_pointer(row, rowp, [nval])
        temp%row = rowp 
        temp%row = temp%row+1
        call c_f_pointer(col, colp, [nval])
        temp%col = colp 
        temp%col = temp%col+1

        ! Housekeeping
        !=============
        end associate
        mysparse = temp

    end subroutine 



    !------------------------------------------------------------------!
    !                            OPERATORS                             !
    !------------------------------------------------------------------!

    ! Add sparse matrices
    function AddSparse(a, b) result(c)

        ! Description
        !============
        ! Overloading of the summation operator for sparse matrices. 
        ! Here, summation is simply carried out by appending the row, 
        ! col, and val indices of both matrices. Checks are performed
        ! on the dimensions provided in nrow, ncol to see if the 
        ! addition is properly defined. 

        ! Declare variables
        !==================
        ! Arguments
        type(MySparseUDT), intent(in)       :: a, b 
        type(MySparseUDT)                   :: c 

        ! Checks
        !=======
        ! Proper dimensions?
        if ( (a%nrow .ne. b%nrow) .or. (a%ncol .ne. b%ncol) ) then 
            call gdErrorHandler('AddSparse: Incompatible matrix' // &
             ' dimensions')
        end if 

        ! Sum
        !====
        ! Initialize 
        c%nrow = a%nrow 
        c%ncol = a%ncol 
        c%nval = a%nval + b%nval 

        ! Allocate
        call c%Allocate()

        ! Sum
        c%row = [a%row, b%row]
        c%col = [a%col, b%col]
        c%val = [a%val, b%val]

    end function

    ! Multiply with scalar
    function MultiplyWithScalarRight(a, b) result(c) 

        ! Description
        !============
        ! Multiply entire matrix with scalar

        ! Declare variables
        !==================
        ! Arguments
        type(MySparseUDT), intent(in)       :: a
        real(R8)                            :: b  
        type(MysparseUDT)                   :: c 

        ! Compute
        !========
        ! Initialize
        c = a

        ! Multiply
        c%val = c%val*b

    end function 

    function MultiplyWithScalarLeft(a, b) result(c) 

        ! Description
        !============
        ! Multiply entire matrix with scalar

        ! Declare variables
        !==================
        ! Arguments
        type(MySparseUDT), intent(in)       :: b
        real(R8)                            :: a 
        type(MysparseUDT)                   :: c 

        ! Compute
        !========
        ! Initialize
        c = b

        ! Multiply
        c%val = c%val*a

    end function 

    ! Multiply element-wise with vector
    function MultiplyWithVectorRowwiseRight(a, b) result(c) 

        ! Description
        !============
        ! Do element-wise multiplication as Aij = Aij*bi

        ! Declare variables
        !==================
        ! Arguments
        type(MySparseUDT), intent(in)       :: a 
        real(R8), allocatable, intent(in)   :: b(:) 
        type(MySparseUDT)                   :: c 

        ! Loop
        integer(I8)                         :: i 

        ! Initialize
        !===========
        ! Checks
        if (size(b) .ne. a%nrow) then 
            ! Throw error
            call gdErrorHandler('SparseMultiplication: incompatible dimensions')
        end if 

        ! Initialize output
        c = a ! Sparsity pattern doesn't change

        ! Loop
        !=====
        do i = 1, a%nval 
            c%val(i) = a%val(i)*b(a%row(i))
        end do

    end function

    function MultiplyWithVectorRowwiseLeft(a, b) result(c) 

        ! Description
        !============
        ! Do element-wise multiplication as Bij = bi*Bij

        ! Declare variables
        !==================
        ! Arguments
        type(MySparseUDT), intent(in)       :: b 
        real(R8), allocatable, intent(in)   :: a(:) 
        type(MySparseUDT)                   :: c 

        ! Loop
        integer(I8)                         :: i 

        ! Initialize
        !===========
        ! Checks
        if (size(a) .ne. b%nrow) then 
            ! Throw error
            call gdErrorHandler('SparseMultiplication: incompatible dimensions')
        end if 

        ! Initialize output
        c = b ! Sparsity pattern doesn't change

        ! Loop
        !=====
        do i = 1, b%nval 
            c%val(i) = b%val(i)*a(b%row(i))
        end do

    end function

    ! Multiply two sparse matrices
    function MultiplySparseMatrices(a, b) result(c) 

        ! Description
        !============
        ! Multiply two sparse matrices by calling the SpMM wrapper 
        ! routine in the C layer that calls the CXSparse routine 
        ! 'cs_multiply'. 

        ! The operation being done is:
        ! 
        !   c_ij = a_ik b_kj,
        !
        ! where summation of k is implied. a and b should have 
        ! compatible dimensions

        ! Modules
        !========
        use Clayer

        ! Declare variables
        !==================
        ! Arguments
        type(MySparseUDT), intent(in)           :: a, b
        type(MySparseUDT)                       :: c 

        ! Auxiliary
        type(CSparseUDT), target                :: csaf, csbf, cscf
        integer(c_int), allocatable, target     :: rowta(:), colta(:)
        real(c_double), allocatable, target     :: valta(:)
        integer(c_int), allocatable, target     :: rowtb(:), coltb(:)
        real(c_double), allocatable, target     :: valtb(:)

        ! Initialize
        !===========
        ! Check if dimensions are compatible
        if (a%ncol .ne. b%nrow) then 
            call gdErrorHandler('MultiplySparseMatrices: dimensions are inconsistent')
        end if 

        ! Determine dimensions of c
        c%nrow = a%nrow
        c%ncol = b%ncol 

        ! Check if one of the matrices is the zero matrix
        if ( (a%nval == 0) .or. (b%nval == 0) ) then 
            ! Return the zero matrix
            c%nval = 0
            call c%Allocate()
            return 
        end if 

        ! Multiply
        !=========
        ! First convert to CS format - note: it appears that setting
        ! the correct c memory location is scope dependent (ore more 
        ! precisely, if we do this in a subroutine, memory goes to shit)
        ! That's why we do it here... 
        rowta = a%row-1 
        colta = a%col-1 
        valta = a%val
        csaf%row = c_loc(rowta)
        csaf%col = c_loc(colta)
        csaf%val = c_loc(valta)
        csaf%nrow = a%nrow 
        csaf%ncol = a%ncol 
        csaf%nval = a%nval

        rowtb = b%row-1 
        coltb = b%col-1 
        valtb = b%val
        csbf%row = c_loc(rowtb)
        csbf%col = c_loc(coltb)
        csbf%val = c_loc(valtb)
        csbf%nrow = b%nrow 
        csbf%ncol = b%ncol 
        csbf%nval = b%nval

        cscf = SpMMF(csaf, csbf)
        
        ! Reconvert
        call ConvertToMySparse(cscf, c)
        
    end function 

    ! Multiply sparse matrices (old)
    function MultiplySparseMatrices_deprecated(a, b) result(c) 

        ! Description
        !============
        ! Multiply two sparse matrices. A very stupid and likely 
        ! inefficient implementation which should by all means by 
        ! replaced by a more powerful alternative in the future 
        ! (e.g. interfacing with CHOLMOD or GraphBLAS of 
        ! the SuiteSparse code suite)

        ! The operation being done is:
        ! 
        !   c_ij = a_ik b_kj,
        !
        ! where summation of k is implied. a and b should have 
        ! compatible dimensions

        ! Declare variables
        !==================
        ! Arguments
        type(MySparseUDT), intent(in)           :: a, b
        type(MySparseUDT)                       :: c 

        ! Auxiliary
        real(R8), allocatable                   :: val(:), tempval(:)
        integer(I8), allocatable                :: row(:), temprow(:), &
            col(:), tempcol(:)
        integer(I8)                             :: grow, nv 
        logical                                 :: doupdate

        ! Loop
        integer(I8)                             :: i, j, k, l

        ! Initialize
        !===========
        ! Check if dimensions are compatible
        if (a%ncol .ne. b%nrow) then 
            call gdErrorHandler('MultiplySparseMatrices: dimensions are inconsistent')
        end if 

        ! Determine dimensions of c
        c%nrow = a%nrow
        c%ncol = b%ncol 

        ! Check if one of the matrices is the zero matrix
        if ( (a%nval == 0) .or. (b%nval == 0) ) then 
            ! Return the zero matrix
            c%nval = 0
            call c%Allocate()
            return 
        end if 
        
        ! Compute
        !========
        ! Set value counter
        nv = 1

        ! Set grow factor
        grow = 2 

        ! Initialize
        allocate(row(a%nval), col(a%nval), val(a%nval))
        val(:) = 0

        ! Loop
        do i = 1, c%nrow 
            do j = 1, c%ncol 
                ! Initialize
                doupdate = .false.

                ! Update value counter
                do k = 1, a%nval 
                    if (a%row(k) == i) then 
                        do l = 1, b%nval 
                            if ( (b%col(l) == j) .and. (a%col(k) == b%row(l)) ) then 
                                val(nv) = val(nv) + a%val(k)*b%val(l)
                                doupdate = .true.
                            end if 
                        end do 
                    end if 
                end do
                
                ! Check if we need to update nv
                if (doupdate) then 
                    ! Set row and column indices of this value
                    row(nv) = i
                    col(nv) = j 

                    ! Check if we need to extend
                    if (nv+1 > size(val, 1)) then 
                        ! Copy and extend
                        allocate(temprow(nv), tempcol(nv), tempval(nv))
                        temprow = row 
                        tempcol = col 
                        tempval = val
                        
                        deallocate(row, col, val)
                        allocate(row(nv*grow), col(nv*grow), val(nv*grow))

                        ! Set value
                        row(1:nv) = temprow(1:nv)
                        col(1:nv) = tempcol(1:nv)
                        val(1:nv) = tempval(1:nv)

                        ! Deallocate
                        deallocate(temprow, tempcol, tempval)
                    end if

                    ! Update counter
                    nv = nv + 1

                end if 
            end do 
        end do 

        ! Downdate
        nv = nv - 1

        ! Build c
        c%nval = nv 
        call c%Allocate()
        c%row = row(1:nv)
        c%col = col(1:nv)
        c%val = val(1:nv)

    end function 

    ! Delete rows
    function DeleteRowsLogical(a, b) result(c)

        ! Description
        !============
        ! Delete the rows of the matrix a that are true in the logical
        ! vector b. b should have adequate dimensions. The matrix's row
        ! dimension will be reduced by COUNT(b).

        ! Declare variables
        !==================
        ! Arguments
        class(MySparseUDT)      :: a 
        type(MySparseUDT)       :: c
        logical, intent(in)     :: b(:) 

        ! Auxiliary
        logical, allocatable    :: reductionvec(:)

        ! Initialize
        !===========
        ! Set output
        c = a 

        ! Check sizes
        if (size(b, 1) .ne. a%nrow) then 
            call gdErrorHandler('DeleteRowsLogical: deletion vector has improper dimensions')
        end if 

        ! Check if allocated
        if (.not. allocated(a%row)) then 
            call gdErrorHandler('DeleteRowsLogical: matrix is not allocated')
        end if 

        ! Delete
        !=======
        ! Adjust dimension
        c%nrow = a%nrow - count(b)
        c%ncol = a%ncol

        ! Remove elements, if any
        if (a%nval > 0) then 
            
            ! Determine elements to retain
            reductionvec = .not. b(a%row)

            ! Construct new matrix
            c%nval  = count(reductionvec)
            c%row   = pack(a%row, reductionvec)
            c%col   = pack(a%col, reductionvec)

        end if 

    end function

    function DeleteRowsInteger(a, b) result(c)

        ! Description
        !============
        ! Same functionality as DeleteRowsLogical, but now the input
        ! vector is an integer array. This routine converts the array
        ! to a logical array (and performs some checks) and calls the
        ! logical deletion routine afterwards. 

        ! Declare variables
        !==================
        ! Arguments
        class(MySparseUDT)      :: a 
        integer(I8), intent(in) :: b(:)
        type(MySparseUDT)       :: c 

        ! Auxiliary
        logical, allocatable    :: bl(:)

        ! Initialize
        !===========
        ! Check b
        if ( (minval(b) < 1) .or. (maxval(b) > a%nrow ) ) then 
            ! Out of bounds, throw error
            call gdErrorHandler('DeleteRowsInteger: some values are out of bounds for deletion')
        end if 

        ! Construct logical
        allocate(bl(a%nrow))
        bl(b) = .true. 

        ! Delete
        !=======
        c = a%DeleteRowsLogical(bl) 


    end function 

    ! Delete columns
    function DeleteColumnsLogical(a, b) result(c)
        
        ! Description
        !============
        ! The same as rows, but now for columns. Here, we simply 
        ! transpose the matrix first, call the row deleter, then 
        ! transpose again. 

        ! Declare variables
        !==================
        ! Arguments
        class(MySparseUDT)      :: a 
        logical, intent(in)     :: b(:)
        type(MySparseUDT)       :: c 

        ! Delete
        !=======
        ! Transpose
        c = a%Transpose()

        ! Delete
        c = c%DeleteRows(b)

        ! Transpose again
        c = c%Transpose()

    end function

    function DeleteColumnsInteger(a, b) result(c)
        
        ! Description
        !============
        ! The same as rows, but now for columns. Here, we simply 
        ! transpose the matrix first, call the row deleter, then 
        ! transpose again. 

        ! Declare variables
        !==================
        ! Arguments
        class(MySparseUDT)      :: a 
        integer(I8), intent(in) :: b(:)
        type(MySparseUDT)       :: c 

        ! Delete
        !=======
        ! Transpose
        c = a%Transpose()

        ! Delete
        c = c%DeleteRows(b)

        ! Transpose again
        c = c%Transpose()

    end function

    ! Sum over columns (returns full)
    subroutine SumColumnwiseFull(mysparse, col)

        ! Description
        !============
        ! Sum the matrix columns into one column, which is returned as
        ! a full vector (usually, this is the most compact format). 
        ! The summation is carried out by looping over all values and
        ! adding on the go. 

        ! Declare variables
        !==================
        ! Arguments
        class(MySparseUDT)          :: mysparse 
        real(R8), allocatable       :: col(:)

        ! Loop variables
        integer(I8)                 :: i

        ! Initialize
        !===========
        ! Check if allocated
        if (.not. allocated(col)) then
            allocate(col(mysparse%nrow))
        else
            ! Check if the dimensions are correct, otherwise reallocate
            if (size(col) .ne. mysparse%nrow) then
                print *, 'incorrect allocation of col, reallocating'
                deallocate(col)
                allocate(col(mysparse%nrow))
            end if
        end if

        ! Associate
        associate(&
            nval    => mysparse%nval,   &
            row     => mysparse%row,    &
            val     => mysparse%val)

        ! Compute
        !========
        ! Initialize
        col(:) = 0

        ! Loop over all values, sum
        do i = 1, nval 
            col(row(i)) = col(row(i)) + val(i)
        end do

        ! End associate
        end associate

    end subroutine

    ! Sum over rows (returns full)
    subroutine SumRowwiseFull(mysparse, row)

        ! Description
        !============
        ! Sum the matrix rows into one row, which is returned as
        ! a full vector (usually, this is the most compact format). 
        ! The summation is carried out by looping over all values and
        ! adding on the go. 

        ! Declare variables
        !==================
        ! Arguments
        class(MySparseUDT)          :: mysparse 
        real(R8), allocatable       :: row(:)

        ! Loop variables
        integer(I8)                 :: i

        ! Initialize
        !===========
        ! Check if allocated
        if (.not. allocated(row)) then
            allocate(row(mysparse%ncol))
        else
            ! Check if the dimensions are correct, otherwise reallocate
            if (size(row) .ne. mysparse%ncol) then
                print *, 'incorrect allocation of row, reallocating'
                deallocate(row)
                allocate(row(mysparse%ncol))
            end if
        end if

        ! Associate
        associate(&
            nval    => mysparse%nval,   &
            col     => mysparse%col,    &
            val     => mysparse%val)

        ! Compute
        !========
        ! Initialize
        row(:) = 0

        ! Loop over all values, sum
        do i = 1, nval 
            row(col(i)) = row(col(i)) + val(i)
        end do

        ! End associate
        end associate

    end subroutine

    ! Transposition
    function Transp(a) result(b)

        ! Description
        !============
        ! Transpose the matrix a

        ! Declare variables
        !==================
        ! Arguments
        class(MySparseUDT), intent(in) :: a 
        type(MySparseUDT)             :: b 

        ! Set
        !====
        b = a 
        b%row = a%col 
        b%col = a%row 
        b%nrow = a%ncol 
        b%ncol = a%nrow

    end function 

    ! Matrix-vector multiplication
    function MatrixVectorProduct(a, b) result(c) 

        ! Description
        !============
        ! Compute the matrix-vector product c = A*b. C is a full vector
        ! of same dimension as b

        ! Declare variables
        !==================
        ! Arguments
        class(MySparseUDT), intent(in)      :: a 
        real(R8), intent(in)                :: b(:) 
        real(R8), allocatable               :: c(:) 

        ! Auxiliary
        integer(I8)                         :: tr, tc

        ! Loop
        integer(I8)                         :: i 

        ! Initialize
        !===========
        ! Checks
        if (size(b) .ne. a%ncol) then 
            print *, size(b), a%ncol
            call gdErrorHandler('MatrixVectorProduct: incompatible dimensions')
        end if 

        ! Allocate
        allocate(c(a%nrow))
        c = 0

        ! Compute
        !========
        do i = 1, a%nval 
            tc = a%col(i)
            tr = a%row(i)
            c(tr) = c(tr) + b(tc)*a%val(i)
        end do 

    end function

    ! Concatenation, 2 sparse matrices
    function Concatenate2(a, b, dim) result(c)

        ! Description
        !============
        ! Concatenate two matrices a and b along dimension dim (can only
        ! be 1 (row) and 2 (column)). Matrix a is always the first 
        ! matrix, so basically this operation yields [a, b] in case of 
        ! column-wise (dim = 2) concatenation, and [a; b] in case of 
        ! rowwise concatenation

        ! Declare variables
        !==================
        ! Arguments
        class(MySparseUDT), intent(in)      :: a, b
        type(MySparseUDT)                   :: c

        ! Auxiliary
        integer(I8), intent(in)             :: dim

        ! Loop
        integer(I8)                         :: i 

        ! Checks
        !=======
        ! Are matrices allocated
        if ( (.not. allocated(a%val)) .or.(.not. allocated(b%val)) ) then 
            call gdErrorHandler('Concatenate2: matrices are not allocated')
        end if 

        ! Is the dimension legal
        if ( (dim > 2) .or. (dim < 1) ) then 
            call gdErrorHandler('Concatenate2: dimension argument must be either 1 or 2')
        end if 

        ! Concatenate
        !============
        if (dim == 1) then 
            ! Concatenate rowwise, i.e. column dimension must be the same
            if ( a%ncol .ne. b%ncol) then 
                call gdErrorHandler('Concatenate2: cannot concatenate ' &
                    // 'over row, column dimensions are inconsistent')
            end if 

            ! Set dimensions
            c%ncol = a%ncol
            c%nrow = a%nrow + b%nrow
            c%nval = a%nval + b%nval 

            ! Allocate
            call c%Allocate()
            
            ! Concatenate
            c%val = [a%val, b%val]
            c%row = [a%row, b%row+a%nrow]
            c%col = [a%col, b%col]

        else
            ! Concatenate columnwise, i.e. row dimension must be the same
            if ( a%nrow .ne. b%nrow) then 
                call gdErrorHandler('Concatenate2: cannot concatenate ' &
                    // 'over column, row dimensions are inconsistent')
            end if 

            ! Set dimensions
            c%nrow = a%nrow
            c%ncol = a%ncol + b%ncol
            c%nval = a%nval + b%nval 

            ! Allocate
            call c%Allocate()
            
            ! Concatenate
            c%val = [a%val, b%val]
            c%row = [a%row, b%row]
            c%col = [a%col, b%col + a%ncol]

        end if 


    end function 

    ! Concatenation, array of N sparse matrices
    function ConcatenateN(b, a, dim) result(c)

        ! Description
        !============
        ! Concatenate N matrices, which should be given in an array
        ! as input. They are concatenated in the same sequence as they
        ! appear in the array. Other functionality is like Concatenate2

        ! Declare variables
        !==================
        ! Arguments
        class(MySparseUDT), intent(in)      :: a(:)
        class(MySparseUDT)                  :: b ! passed object must be scalar
        type(MySparseUDT)                   :: c

        ! Auxiliary
        integer(I8), intent(in)             :: dim
        integer(I8)                         :: n, nv, nr, nc

        ! Loop
        integer(I8)                         :: i 

        ! Checks
        !=======
        ! Get the size
        n = size(a, 1)

        ! Check
        if (n == 0) then 
            ! Empty, return empty matrix
            c%nrow = 0
            c%ncol = 0
            c%nval = 0
            call c%Allocate()
            return 
        end if 

        ! Are matrices allocated
        do i = 1, n
            if ( .not. allocated(a(i)%val) ) then 
                call gdErrorHandler('ConcatenateN: matrices are not allocated')
            end if 
        end do 

        ! Is the dimension legal
        if ( (dim > 2) .or. (dim < 1) ) then 
            call gdErrorHandler('ConcatenateN: dimension argument must be either 1 or 2')
        end if 

        ! Concatenate
        !============
        if (dim == 1) then 
            ! Concatenate rowwise, i.e. column dimension must be the same
            do i = 2, n
                if ( a(1)%ncol .ne. a(i)%ncol) then 
                    call gdErrorHandler('ConcatenateN: cannot concatenate ' &
                        // 'over row, column dimensions are inconsistent')
                end if 
            end do 

            ! Set dimensions
            c%ncol = a(1)%ncol
            c%nrow = 0
            c%nval = 0
            do i = 1, n
                c%nrow = c%nrow + a(i)%nrow
                c%nval = c%nval + a(i)%nval 
            end do

            ! Allocate
            call c%Allocate()
            
            ! Concatenate
            nv = 0
            nr = 0
            do i = 1, n 
                ! Add
                c%val(nv+1:nv+a(i)%nval) = a(i)%val
                c%row(nv+1:nv+a(i)%nval) = a(i)%row + nr
                c%col(nv+1:nv+a(i)%nval) = a(i)%col

                ! Update counters
                nr = nr + a(i)%nrow
                nv = nv + a(i)%nval
            end do

        else
            ! Concatenate columnwise, i.e. row dimension must be the same
            do i = 2, n
                if ( a(1)%nrow .ne. a(i)%nrow) then 
                    call gdErrorHandler('ConcatenateN: cannot concatenate ' &
                        // 'over column, row dimensions are inconsistent')
                end if 
            end do 

            ! Set dimensions
            c%nrow = a(1)%nrow
            c%ncol = 0
            c%nval = 0
            do i = 1, n
                c%ncol = c%ncol + a(i)%ncol
                c%nval = c%nval + a(i)%nval 
            end do

            ! Allocate
            call c%Allocate()
            
            ! Concatenate
            nv = 0
            nc = 0
            do i = 1, n 
                ! Add
                c%val(nv+1:nv+a(i)%nval) = a(i)%val
                c%row(nv+1:nv+a(i)%nval) = a(i)%row 
                c%col(nv+1:nv+a(i)%nval) = a(i)%col + nc

                ! Update counters
                nc = nc + a(i)%ncol
                nv = nv + a(i)%nval
            end do

        end if 


    end function 

    !------------------------------------------------------------------!
    !                           HOUSEKEEPING                           !
    !------------------------------------------------------------------!
    ! Allocation
    subroutine AllocateMySparse(mysparse)

        ! Description
        !============
        ! Allocate. It is assumed that nrow, ncol, nval are set to the
        ! correct values. 

        ! Initialize
        !===========
        ! The usual
        implicit none

        ! Declare variables
        !==================
        ! Arguments
        class(MySparseUDT)       :: mysparse

        ! Allocate
        !=========
        if (.not. allocated(mysparse%col)) then 
            allocate(mysparse%col(mysparse%nval))
            allocate(mysparse%row(mysparse%nval))
            allocate(mysparse%val(mysparse%nval))
        end if 

    end subroutine

    ! Deallocation
    subroutine DeallocateMySparse(mysparse)

        ! Description
        !============
        ! Deallocate.

        ! Initialize
        !===========
        ! The usual
        implicit none

        ! Declare variables
        !==================
        ! Arguments
        class(MySparseUDT)       :: mysparse

        ! Allocate
        !=========
        if (allocated(mysparse%col)) then 
            ! assume rest is also allocated
            deallocate(mysparse%col)
            deallocate(mysparse%row)
            deallocate(mysparse%val)
        end if

    end subroutine

    ! Destructor
    subroutine DestroyMySparse(mysparse)

        ! Description
        !============
        ! Destructor, to be called when the object should be destroyed

        ! Initialize
        !===========
        ! The usual
        implicit none

        ! Declare variables
        !==================
        type(MySparseUDT), intent(in)       :: mysparse

        ! Destroy
        !========
        call mysparse%Deallocate()

    end subroutine

end module