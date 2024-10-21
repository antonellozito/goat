!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains algorithms for sorting (1D) arrays. Only certain
! precisions are available - cast/recast if necessary, or implement
! dedicated routines. Only the general sort interfaces are made
! publicly available. 

! Note 1: the sorting algorithms defined here are not the most performant
! ones. In the future, one might resort to linking to either standard
! fortran sorting routines (if that ever becomes available) or to C
! implementations (there is actually qsort in C, but it does not
! seem to support all the desired functionality, e.g. returning the sort
! index). 

module mod_sort

    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_errorhandler
    
    ! The usual
    implicit none 
    private 
    public :: Sort 

    !==================================================================!
    !                                                                  !
    !                            TYPES                                 !
    !                                                                  !
    !==================================================================!

    ! Currently no additional types to define

    !==================================================================!
    !                                                                  !
    !                          INTERFACES                              !
    !                                                                  !
    !==================================================================!

    ! General sorting
    interface Sort 
        module procedure Sort_I8, Sort_R8
    end interface 

    contains 

    !==================================================================!
    !                                                                  !
    !                           ROUTINES                               !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                           INTEGER SORTING                        !
    !------------------------------------------------------------------!

    ! Integer sorter
    subroutine Sort_I8(a, ind, ascend, algorithm)

        ! Description
        !============
        ! Wrapper for sorting functions. The 'ind' array will yield the 
        ! indices that can be used to sort the array again by doing
        ! a_sorted = a_unsorted(ind)

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(inout)          :: a(:)
        integer(I8), intent(out), optional  :: ind(size(a))
        logical, intent(in), optional       :: ascend 
        character(*), intent(in), optional  :: algorithm 

        ! Auxiliary
        character(:), allocatable           :: meth 
        logical                             :: flip, doindex
        logical, allocatable                :: test(:)

        ! Loop
        integer(I8)                         :: k 

        ! Initialize
        !===========
        ! Check inputs
        if (present(algorithm)) then 
            meth = algorithm 
        else 
            ! Default: quicksort in place
            meth = 'quicksort_inplace'
        end if 
        if (present(ascend)) then 
            flip = .not. ascend 
        else 
            ! Don't flip
            flip = .false.
        end if

        ! Initialize index vector
        if (present(ind)) then 
            doindex = .true. 
        else
            doindex = .false. 
        end if 

        ! Call sorter
        !============
        select case (meth)

        case ('quicksort', 'quicksort_inplace')

            ! Call sorter
            if (present(ind)) then 
                ind = [(k, k = 1, size(a))] 
                call Quicksort_inplace_ind_I8(a, ind)
            else
                call Quicksort_inplace_I8(a)
            end if 

        case default

            ! Throw error
            call gdErrorHandler('Sort_I8: algorithm: ' // meth // ' not yet ' // &
                'implemented for integers')

        end select 

        ! Post-process
        !=============
        ! Do a sanity check - you never know
        test = (a(2:size(a)) - a(1:size(a)-1)) < 0
        if (any(test)) then 
            call gdErrorHandler('Sort_I8: bug detected: array is not sorted')
        end if 

        if (flip) then 
            a = a(size(a):1:-1)
            if (doindex) then 
                ind = ind(size(a):1:-1)
            end if 
        end if 

    end subroutine

    ! Quicksort (in-place, ind)
    recursive subroutine Quicksort_inplace_ind_I8(a, ind)

        ! Description
        !============
        ! Basic in-place quicksort algorithm. When provided, the 'ind'
        ! array is sorted in the same way as the 'a' array. This is 
        ! useful to find the sorting indices and sort any other arrays
        ! later on

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(inout)      :: a(:), ind(:)

        ! Auxiliary
        integer(I8)                     :: left, right, pivot, first, &
            last, t, tind 

        ! Initialize
        !===========
        first = 1
        last = size(a)

        ! Loop
        !=====
        if (size(a) > 1) then 
            ! Determine pivot and set initial bounds
            pivot = a(size(a)/2)
            left = 1
            right = size(a)

            ! Partition by swapping elements that are inverted (i.e. 
            ! one element at the left higher that the pivot, one 
            ! element at the right lower than the pivot)
            do while (left <= right) 
                do while (a(left) < pivot)
                    left = left + 1
                end do
                do while (a(right) > pivot)
                    right = right - 1
                end do
                if (left <= right) then 
                    ! Swap
                    t = a(left)
                    tind = ind(left)
                    a(left) = a(right)
                    ind(left) = ind(right)
                    a(right) = t 
                    ind(right) = tind
                    left = left + 1
                    right = right - 1
                end if  
            end do 

            ! Keep sorting the left and right side until the trivial
            ! case is encountered
            call Quicksort_inplace_ind_I8(a(first:right), ind(first:right))
            call Quicksort_inplace_ind_I8(a(left:last), ind(left:last))

        end if 



    end subroutine

    ! Quicksort (in-place)
    recursive subroutine Quicksort_inplace_I8(a)

        ! Description
        !============
        ! Basic in-place quicksort algorithm. 

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(inout)      :: a(:)

        ! Auxiliary
        integer(I8)                     :: left, right, pivot, first, &
            last, t

        ! Initialize
        !===========
        first = 1
        last = size(a)

        ! Loop
        !=====
        if (size(a) > 1) then 
            ! Determine pivot and set initial bounds
            pivot = a(size(a)/2)
            left = 1
            right = size(a)

            ! Partition by swapping elements that are inverted (i.e. 
            ! one element at the left higher that the pivot, one 
            ! element at the right lower than the pivot)
            do while (left <= right) 
                do while (a(left) < pivot)
                    left = left + 1
                end do
                do while (a(right) > pivot)
                    right = right - 1
                end do
                if (left <= right) then 
                    ! Swap
                    t = a(left)
                    a(left) = a(right)
                    a(right) = t 
                    left = left + 1
                    right = right - 1
                end if  
            end do 

            ! Keep sorting the left and right side until the trivial
            ! case is encountered
            call Quicksort_inplace_I8(a(first:right))
            call Quicksort_inplace_I8(a(left:last))

        end if 



    end subroutine

    !------------------------------------------------------------------!
    !                           REAL SORTING                           !
    !------------------------------------------------------------------!

    ! Real sorter
    subroutine Sort_R8(a, ind, ascend, algorithm)

        ! Description
        !============
        ! Wrapper for sorting functions. The 'ind' array will yield the 
        ! indices that can be used to sort the array again by doing
        ! a_sorted = a_unsorted(ind)

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(inout)             :: a(:)
        integer(I8), intent(out), optional  :: ind(size(a))
        logical, intent(in), optional       :: ascend 
        character(*), intent(in), optional  :: algorithm 

        ! Auxiliary
        character(:), allocatable           :: meth 
        logical                             :: flip, doindex
        logical, allocatable                :: test(:)

        ! Loop
        integer(I8)                         :: k 

        ! Initialize
        !===========
        ! Check inputs
        if (present(algorithm)) then 
            meth = algorithm 
        else 
            ! Default: quicksort in place
            meth = 'quicksort_inplace'
        end if 
        if (present(ascend)) then 
            flip = .not. ascend 
        else 
            ! Don't flip
            flip = .false.
        end if

        ! Initialize index vector
        if (present(ind)) then 
            doindex = .true. 
        else
            doindex = .false. 
        end if 

        ! Call sorter
        !============
        select case (meth)

        case ('quicksort', 'quicksort_inplace')

            ! Call sorter
            if (present(ind)) then 
                ind = [(k, k = 1, size(a))] 
                call Quicksort_inplace_ind_R8(a, ind)
            else
                call Quicksort_inplace_R8(a)
            end if 

        case default

            ! Throw error
            call gdErrorHandler('Sort_R8: algorithm: ' // meth // ' not yet ' // &
                'implemented for integers')

        end select 

        ! Post-process
        !=============
        ! Do a sanity check - you never know
        test = (a(2:size(a)) - a(1:size(a)-1)) < 0
        if (any(test)) then 
            call gdErrorHandler('Sort_R8: bug detected: array is not sorted')
        end if 

        if (flip) then 
            a = a(size(a):1:-1)
            if (doindex) then 
                ind = ind(size(a):1:-1)
            end if 
        end if 

    end subroutine

    ! Quicksort (in-place, ind)
    recursive subroutine Quicksort_inplace_ind_R8(a, ind)

        ! Description
        !============
        ! Basic in-place quicksort algorithm. When provided, the 'ind'
        ! array is sorted in the same way as the 'a' array. This is 
        ! useful to find the sorting indices and sort any other arrays
        ! later on

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(inout)         :: a(:)
        integer(I8), intent(inout)      :: ind(:)

        ! Auxiliary
        real(R8)                        :: pivot, t
        integer(I8)                     :: left, right, first, &
            last, tind 

        ! Initialize
        !===========
        first = 1
        last = size(a)

        ! Loop
        !=====
        if (size(a) > 1) then 
            ! Determine pivot and set initial bounds
            pivot = a(size(a)/2)
            left = 1
            right = size(a)

            ! Partition by swapping elements that are inverted (i.e. 
            ! one element at the left higher that the pivot, one 
            ! element at the right lower than the pivot)
            do while (left <= right) 
                do while (a(left) < pivot)
                    left = left + 1
                end do
                do while (a(right) > pivot)
                    right = right - 1
                end do
                if (left <= right) then 
                    ! Swap
                    t = a(left)
                    tind = ind(left)
                    a(left) = a(right)
                    ind(left) = ind(right)
                    a(right) = t 
                    ind(right) = tind
                    left = left + 1
                    right = right - 1
                end if  
            end do 

            ! Keep sorting the left and right side until the trivial
            ! case is encountered
            call Quicksort_inplace_ind_R8(a(first:right), ind(first:right))
            call Quicksort_inplace_ind_R8(a(left:last), ind(left:last))

        end if 



    end subroutine

    ! Quicksort (in-place)
    recursive subroutine Quicksort_inplace_R8(a)

        ! Description
        !============
        ! Basic in-place quicksort algorithm. 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(inout)         :: a(:)

        ! Auxiliary
        real(R8)                        :: t, pivot
        integer(I8)                     :: left, right, first, &
            last

        ! Initialize
        !===========
        first = 1
        last = size(a)

        ! Loop
        !=====
        if (size(a) > 1) then 
            ! Determine pivot and set initial bounds
            pivot = a(size(a)/2)
            left = 1
            right = size(a)

            ! Partition by swapping elements that are inverted (i.e. 
            ! one element at the left higher that the pivot, one 
            ! element at the right lower than the pivot)
            do while (left <= right) 
                do while (a(left) < pivot)
                    left = left + 1
                end do
                do while (a(right) > pivot)
                    right = right - 1
                end do
                if (left <= right) then 
                    ! Swap
                    t = a(left)
                    a(left) = a(right)
                    a(right) = t 
                    left = left + 1
                    right = right - 1
                end if  
            end do 

            ! Keep sorting the left and right side until the trivial
            ! case is encountered
            call Quicksort_inplace_R8(a(first:right))
            call Quicksort_inplace_R8(a(left:last))

        end if 



    end subroutine

end module