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

! Note 2: some derivatives of sorting algorithms are also provided, e.g.
! the 'unique' function (which is here based on sorting)

module mod_sort

    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_errorhandler
    
    ! The usual
    implicit none 
    private 
    public :: Sort, Unique, Setdiff, CountOccurrence, SearchSortedArray

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

    ! General unique
    interface Unique 
        module procedure Unique_I8, Unique_R8
    end interface

    ! General set difference
    interface Setdiff 
        module procedure Setdiff_I8
    end interface

    ! General element occurrence counter
    interface CountOccurrence
        module procedure CountOccurrence_I8
    end interface

    ! General binary search of a sorted array
    interface SearchSortedArray
        module procedure BinarySearch_I8
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
        if (size(a) > 1) then 
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
        if (size(a) > 1) then 
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

    !------------------------------------------------------------------!
    !                         DERIVED OPERATORS                        !
    !------------------------------------------------------------------!

    ! Integer unique
    subroutine Unique_I8(in, out, ind)

        ! Description
        !============
        ! Returns the unique elements of an integer array. This is done
        ! by first sorting the array, then taking the value difference
        ! and removing elements which are not different. Optionally, 
        ! the index vector is returned that yields out = in(ind)

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)                 :: in(:) 
        integer(I8), allocatable, intent(out)   :: out(:)
        integer(I8), allocatable, intent(out), optional     :: ind(:)

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: diff, in_sorted
        logical, allocatable, dimension(:)      :: keepind

        ! Loop
        integer(I8)                             :: k 

        ! Compute
        !========
        if (allocated(out)) then 
            deallocate(out)
        end if 
        if (size(in) == 0) then 
            out = in 
            return 
        end if 

        in_sorted = in 
        if (present(ind)) then 
            ! Sort
            call Sort(in_sorted, ind=ind)

            ! Compute difference
            diff = in_sorted(2:) - in_sorted(1:size(in_sorted)-1)

            ! Eliminate
            keepind = [.true., diff /= 0_I8]
            allocate(out(count(keepind)))
            out = pack(in_sorted, keepind)
            ind = pack([(k, k = 1, size(in_sorted))], keepind)
        else 
            ! Sort 
            call Sort(in_sorted)

            ! Compute difference
            diff = in_sorted(2:size(in_sorted)) - in_sorted(1:size(in_sorted)-1)

            ! Eliminate
            keepind = [.true., diff /= 0_I8]
            allocate(out(count(keepind)))
            out = pack(in_sorted, keepind)
        end if 

    end subroutine

    ! Real unique
    subroutine Unique_R8(in, out, ind)

        ! Description
        !============
        ! Returns the unique elements of a real array. This is done
        ! by first sorting the array, then taking the value difference
        ! and removing elements which are not different.

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)                                :: in(:) 
        real(R8), allocatable, intent(out)                  :: out(:)
        integer(I8), allocatable, intent(out), optional     :: ind(:)

        ! Auxiliary
        real(R8), allocatable, dimension(:)     :: diff, in_sorted
        logical, allocatable, dimension(:)      :: keepind

        ! Loop
        integer(I8)                             :: k 

        ! Compute
        !========
        if (allocated(out)) then 
            deallocate(out)
        end if 
        if (size(in) == 0) then 
            out = in 
            return 
        end if 

        in_sorted = in 
        if (present(ind)) then 
            ! Sort
            call Sort(in_sorted, ind=ind)

            ! Compute difference
            diff = in_sorted(2:) - in_sorted(1:size(in_sorted)-1)

            ! Eliminate
            keepind = [.true., diff /= 0_I8]
            allocate(out(count(keepind)))
            out = pack(in_sorted, keepind)
            ind = pack([(k, k = 1, size(in_sorted))], keepind)
        else 
            ! Sort 
            call Sort(in_sorted)

            ! Compute difference
            diff = in_sorted(2:) - in_sorted(1:size(in_sorted))

            ! Eliminate
            keepind = [.true., diff /= 0_I8]
            allocate(out(count(keepind)))
            out = pack(in_sorted, keepind)
        end if 

    end subroutine

    ! Integer setdiff
    subroutine Setdiff_I8(a, b, out)

        ! Description
        !============
        ! Compute the set difference between sets a and b (i.e. we return
        ! only the elements of a that are not in b). The result will be 
        ! sorted and unique. if ind is present, the index is returned
        ! that would yield out = a(ind). 

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in), dimension(:)   :: a, b 
        integer(I8), intent(out), dimension(:), allocatable  :: out 

        ! Auxiliary
        integer(I8)                             :: na
        integer(I8), allocatable, dimension(:)  :: c, &
            sortind, c_orig, au, bu
        logical, allocatable, dimension(:)      :: ina, inb, keepind

        ! Loop
        integer(I8)                             :: i 

        ! Hedge for limit cases
        !======================
        if (size(b) == 0) then 
            call Unique(a, out)
            return 
        end if 
        if (size(a) == 0) then 
            out = a 
            return 
        end if 

        ! Sort
        !=====
        ! Take unique of a and b
        call Unique(a, au)
        call Unique(b, bu)

        ! Get size
        na = size(au)

        ! Concatente
        c = [au, bu]
        c_orig = c 
        allocate(sortind(size(c)))
        call Sort(c, ind=sortind, ascend=.true.)
        !call Unique(c, cu, ind=uniqueind)
        !sortind = sortind(uniqueind)
        !c = cu

        ! Eliminate
        !==========
        ! Initialize
        allocate(keepind(size(c)))
        keepind = .true. 
        ina = sortind <= na 
        inb = sortind > na

        ! Remove elements 
        i = 1
        do while (i < size(c))
            if (c(i) == c(i+1)) then 
                keepind(i:i+1) = .false. 
            end if 
            i = i + 1
        end do 

        ! Set output
        allocate(out(count(keepind)))
        out = pack(c, keepind)

    end subroutine

    ! Integer occurrence counter
    subroutine CountOccurrence_I8(in, oc, el, sortindoc, sortindel)

        ! Description
        !============
        ! This routine counts how many times a certain element appears.
        ! To this end, we follow a similar approach as the unique 
        ! function and first sort the array to afterwards count the 
        ! number of sequently same elements (instead of removing them).
        ! The (unique) elements are returned in the 'el' array, and the
        ! number of occurrences 'oc' are returned for this sorted array.
        ! If one wants to map the occurences back to the original 
        ! array 'in' (which includes the doubles), one can use the 'sortindoc'
        ! optional output argument (then oc_in = oc(sortindoc) and in = el(sortindel))

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), dimension(:), intent(in)                           :: in 
        integer(I8), dimension(:), intent(out), allocatable             :: oc, el 
        integer(I8), dimension(:), intent(out), allocatable, optional   :: &
            sortindel, sortindoc

        ! Auxiliary
        integer(I8)                                 :: nel 
        integer(I8), dimension(:), allocatable      :: in_sorted, &
            sortind, diff, indoc, uniqueind
        logical, dimension(:), allocatable          :: keepind

        ! Loop
        integer(I8)                                 :: i, k

        ! Count
        !======
        ! Initialize
        in_sorted = in 
        allocate(sortind(size(in)))

        ! Sort
        call Sort(in_sorted, ind=sortind, ascend=.true.)
        indoc = sortind 

        ! Compute difference
        diff = in_sorted(2:size(in_sorted)) - in_sorted(1:size(in_sorted)-1)

        ! Count
        keepind = [.true., diff /= 0_I8]
        nel = count(keepind)
        allocate(oc(nel), el(nel), uniqueind(nel))
        uniqueind = pack(sortind, keepind)
        oc = 0
        el = 0
        if (nel > 0) then 
            oc(1) = 1
            el(1) = in_sorted(1)
            indoc(1) = 1
        end if 
        k = 1
        do i = 1, size(diff)
            if (diff(i) /= 0_I8) then 
                ! Update element counter and add element
                k = k + 1
                el(k) = in_sorted(i+1)
                oc(k) = oc(k) + 1
            else
                oc(k) = oc(k) + 1
            end if 
            indoc(i+1) = k
        end do 

        ! Output
        if (present(sortindel)) then 
            sortindel = sortind 
        end if 
        if (present(sortindoc)) then 
            sortindoc = indoc ! for initialization only
            sortindoc(sortind) = indoc
        end if 
        
    end subroutine

    ! Binary search for sorted integer arrays
    function BinarySearch_I8(in, el) result(ind)

        ! Description
        !============
        ! Find the index of the element with value 'el' in the input 
        ! array 'in', where 'in' is assumed to be sorted. We do this in 
        ! a binary search way to achieve O(log(N)) behavior instead of 
        ! regular O(N) for unsorted arrays. If the element is not 
        ! present, then 0 is returned as index. 

        ! Declare arguments
        !==================
        ! Arguments
        integer(I8), dimension(:), intent(in)       :: in 
        integer(I8), intent(in)                     :: el 
        integer(I8)                                 :: ind 

        ! Auxiliary
        logical                     :: ascend 
        integer(I8)                 :: nin, low, mid, high, midel

        ! Initialize
        !===========
        ! Set default value
        ind = 0
        nin = size(in)

        ! Hedge for trivial cases
        if (nin == 0) then 
            return 
        end if 
        if (nin == 1) then 
            if (in(1) == el) then 
                ind = 1
            end if 
            return 
        end if 

        ! Check sorting direction
        ascend = .true. 
        if (in(1) > in(nin)) then 
            ascend = .false. 
        end if 

        ! Search
        !=======
        ! Initialize
        low = 1
        high = nin 
        mid = (low + high)/2 ! intended integer division

        ! Loop
        if (ascend) then 
            do while (.true.)
                ! Checks if mid equals el
                midel = in(mid)
                if (in(mid) == el) then
                    ind = mid  
                    exit 
                end if 

                ! If not, check 
                if (el < midel) then 
                    high = mid - 1
                else
                    low = mid + 1
                end if 
                mid = (high + low)/2

                ! Some checks for now
                if (mid < low .or. mid > high) then 
                    call gdErrorHandler('something wrong in binary search')
                end if 

            end do 
        else 
            do while (.true.)
                ! Checks if mid equals el
                midel = in(mid)
                if (in(mid) == el) then
                    ind = mid  
                    exit 
                end if 

                ! If not, check 
                if (el > midel) then 
                    low = mid - 1
                else
                    high = mid + 1
                end if 
                mid = (high + low)/2

                ! Some checks for now
                if (mid < low .or. mid > high) then 
                    call gdErrorHandler('something wrong in binary search')
                end if 

            end do 
        end if 



    end function 



end module