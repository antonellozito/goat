!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains auxiliary methods that are (typically) generally
! useable for any 2D interpolant. 

module Interpolant2D_auxiliaries

    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_errorhandler
    use mod_sort
    use omp_lib

    implicit none

    interface GetIndex 
        module procedure GetLinearIndex
        module procedure GetMatrixIndices
    end interface

    contains 

    ! Factorial computation
    integer(I16) function MyFactorial(n)

        ! Description
        !============
        ! Simple factorial computation routine. Not accurate for n >> 10
        ! likely. 

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)             :: n

        ! Loop
        integer(I8)                         :: i 

        ! Compute
        !========
        ! Initialize
        MyFactorial = 1
        if (n == 0) then 
            ! Exit
            return 
        end if 
        do i = 1, n 
            MyFactorial = MyFactorial * i
        end do 

    end function 

    ! Binner (index retriever), 2D
    subroutine GetLinearIndex(xq, yq, x, y, ind)

        ! Description
        !============
        ! Determine in which 'bin' or cell the query points lie. A point
        ! lies in the k = i + (j-1)*(nx - 1) bin if the following holds:
        !
        !    x(i) <= xq < x(i+1) (for i = 1:nx-1)
        !    y(j) <= yq < y(i+1) (for j = 1:ny-1)
        !
        ! Points that are outside of the domain will get a zero value 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)                    :: xq(:), yq(:), x(:), y(:)
        integer(I8), allocatable, intent(out)   :: ind(:)

        ! Loop variables
        integer(I8)                         :: k, indx, indy

        ! Auxiliary variables 
        integer(I8)                         :: nq, nx, ny

        ! Compute indices
        !================
        ! Check
        if (allocated(ind)) then 
            deallocate(ind)
        end if
        allocate(ind(size(xq, 1)))

        ! Compute sizes
        nq = size(xq, 1)
        nx = size(x, 1)
        ny = size(y, 1)

        ! Loop over all query points
        !$omp parallel do default(none) private(indx, indy, k) &
        !$omp if((.not. omp_in_parallel()) .and. (nq > 100)) schedule(dynamic) & 
        !$omp shared(nq, xq, x, nx, yq, y, ny, ind)
        do k = 1, nq
            ! Get the bin index for x
            call GetBinIndex(xq(k), x, nx, indx)

            ! Get the bin index for y
            call GetBinIndex(yq(k), y, ny, indy)

            ! Compute the total bin index
            ind(k) = indx + (indy - 1)*(nx - 1)
            if ((indx == 0) .or. (indy == 0)) then
                ! Out of bounds, put index to zero
                ind(k) = 0
            end if
            
        end do
        !$omp end parallel do

    end subroutine

    subroutine GetMatrixIndices(xq, yq, x, y, ix, iy)

        ! Description
        !============
        ! Determine in which 'bin' or cell the query points lie. A point
        ! lies in the k = i + (j-1)*(nx - 1) bin if the following holds:
        !
        !    x(i) <= xq < x(i+1) (for i = 1:nx-1)
        !    y(j) <= yq < y(i+1) (for j = 1:ny-1)
        !
        ! Points that are outside of the domain will get a zero value 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)                    :: xq(:), yq(:), x(:), y(:)
        integer(I8), allocatable, intent(out)   :: ix(:), iy(:)

        ! Loop variables
        integer(I8)                         :: k, indx, indy

        ! Auxiliary variables 
        integer(I8)                         :: nq, nx, ny

        ! Compute indices
        !================
        ! Check
        if (allocated(ix)) then 
            deallocate(ix)
        end if
        if (allocated(iy)) then 
            deallocate(iy)
        end if
        allocate(ix(size(xq, 1)), iy(size(xq, 1)))

        ! Compute sizes
        nq = size(xq, 1)
        nx = size(x, 1)
        ny = size(y, 1)

        ! Loop over all query points
        !$omp parallel do default(none) private(indx, indy, k) &
        !$omp if(.not. omp_in_parallel()) schedule(dynamic) & 
        !$omp shared(nq, xq, x, nx, yq, y, ny, ix, iy)
        do k = 1, nq
            ! Get the bin index for x
            call GetBinIndex(xq(k), x, nx, indx)

            ! Get the bin index for y
            call GetBinIndex(yq(k), y, ny, indy)

            ! Add to output
            ix(k) = indx 
            iy(k) = indy
            
        end do
        !$omp end parallel do 

    end subroutine

    ! Get index 1D

    subroutine GetBinIndex(xq, x, nx, ind)

        ! Description
        !============
        ! Get the bin index of a point xq (scalar) in the array x. Not
        ! optimized. It is assumed that x is monotonously increasing.
        ! If the value is outside of x, the index value 0 is returned.
        
        ! A point lies in the ith bin if:
        !
        !    x(i) <= xq < x(i+1) (for i = 1:nx-1)
        !
        ! Points that are outside of the domain will get a zero value 
        ! and are eliminated first 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)            :: xq, x(:)
        integer, intent(out)            :: ind 
        
        ! Loop variables
        
        ! Auxiliary
        integer(I8)                     :: nx

        ! Compute the bin
        !================
        ind = GetBinIndexSortedArray(x, xq)

    end subroutine

    subroutine GetBinIndex_old(xq, x, nx, ind)

        ! Description
        !============
        ! Get the bin index of a point xq (scalar) in the array x. Not
        ! optimized. It is assumed that x is monotonously increasing.
        ! If the value is outside of x, the index value 0 is returned.
        
        ! A point lies in the ith bin if:
        !
        !    x(i) <= xq < x(i+1) (for i = 1:nx-1)
        !
        ! Points that are outside of the domain will get a zero value 
        ! and are eliminated first 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)            :: xq, x(:)
        integer, intent(out)            :: ind 
        
        ! Loop variables
        
        ! Auxiliary
        integer(I8)                     :: nx
        logical                         :: notfound

        ! Compute the bin
        !================
        ! Determine the position in x (not optimized)
        ind = 1
        notfound = .true. 

        ! Check if the value lies inside the range
        if ((xq < (x(1))) .or. (xq > x(size(x, 1)))) then
            ! Outside of bin
            notfound = .false. 
            ind = 0
            return 
        end if

        ! Check if there are NaNs - then return
        if (isnan(xq)) then 
            ind = 0
            return 
        end if

        ! Check the first bin - equality holds for first value
        if ((xq >= x(ind)) .and. (xq <= x(ind+1))) then
            notfound = .false.
            ind = 1
            return 
        else
            ! Not found, increase index
            ind = ind+1
        end if

        do while ((ind < nx))
            if ((xq > x(ind)) .and. (xq <= x(ind+1))) then
                ! Found, exit the loop
                notfound = .false.
                exit
            else 
               ! Update ind
                ind = ind+1
            end if
        end do

        ! Hedge for NaNs etc
        if (notfound) then 
            ! Should be NaN - set out of bounds
            ind = 0
        end if 


    end subroutine

    

end module