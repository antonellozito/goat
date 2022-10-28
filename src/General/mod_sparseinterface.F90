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

module mod_sparseinterface
    
    ! Initialize
    !============
    ! Load modules
    use mod_precision

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
        allocate(mysparse%col(mysparse%nval))
        allocate(mysparse%row(mysparse%nval))
        allocate(mysparse%val(mysparse%nval))

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