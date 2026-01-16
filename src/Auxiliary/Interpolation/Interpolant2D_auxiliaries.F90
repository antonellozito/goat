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
    use mod_triangulation

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

    ! InTriangle routine
    subroutine InTriangle(v1, v2, v3, x, y, xp, yp, in, on)

        ! Description
        !============
        ! This routine checks if a point lies within or on a triangle.
        ! The result is returned in the in and on arrays,
        ! as input the triangle vertices should be given in v1, v2, v3. 
        ! The triangle coordinates, which can be queried as x(v1), should
        ! be given in the x, y arrays, and the query point is given as a
        ! scalar set of coordinates xp, yp. 

        ! Note: the implementation here is very naive and simply loops
        ! over all triangles, assuming there are not many to deal with.
        ! This assumption is only valid within this module, since 
        ! this routine is only used for saddle points which have few
        ! triangles. 

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), dimension(:), intent(in)           :: v1, v2, v3 
        real(R8), dimension(:), intent(in)              :: x, y 
        real(R8), intent(in)                            :: xp, yp 
        logical, allocatable, dimension(:), intent(out) :: in, on

        ! Auxiliary
        real(R8)                        :: x1, y1, dx1, dy1, dx1p, dy1p, &
            x2, y2, dx2, dy2, dx2p, dy2p, x3, y3, dx3, dy3, dx3p, dy3p, & 
            cp(1:3)
        logical                         :: do_parallel

        ! Loop
        integer(I8)                         :: i 

        ! Checks
        !=======
        if ((size(v1) /= size(v2)) .or. (size(v1) /= size(v3))) then 
            call gdErrorHandler('InTriangle: v1, v2, v3 have incompatible sizes')
        end if 
        if (size(x) /= size(y)) then 
            call gdErrorHandler('InTriangle: x, y have incompatible sizes')
        end if 
        if (allocated(in)) then 
            deallocate(in)
        end if 
        if (allocated(on)) then
            deallocate(on)
        end if 

        ! Initialize
        !===========
        allocate(in(size(v1)), on(size(v1)))
        in = .false. 
        on = .false. 

        ! Compute
        !========
        ! Loop
        do_parallel = .not. omp_in_parallel()
        !$omp parallel do default(none) schedule(static) if (do_parallel) &
        !$omp shared(v1, v2, v3, x, y, in, on, xp, yp) &
        !$omp private(i, x1, x2, x3, y1, y2, y3, dx1, dx2, dx3, dy1, dy2, dy3, &
        !$omp dx1p, dx2p, dx3p, dy1p, dy2p, dy3p, cp)
        do i = 1, size(v1)
            ! Get coordinates
            x1 = x(v1(i))
            x2 = x(v2(i))
            x3 = x(v3(i))
            y1 = y(v1(i))
            y2 = y(v2(i))
            y3 = y(v3(i))

            ! Get vectors
            dx1 = x2 - x1 
            dx2 = x3 - x2 
            dx3 = x1 - x3
            dy1 = y2 - y1 
            dy2 = y3 - y2 
            dy3 = y1 - y3

            dx1p = xp - x1 
            dx2p = xp - x2 
            dx3p = xp - x3
            dy1p = yp - y1 
            dy2p = yp - y2 
            dy3p = yp - y3

            ! Compute cross products
            cp(1) = dx1*dy1p - dy1*dx1p 
            cp(2) = dx2*dy2p - dy2*dx2p
            cp(3) = dx3*dy3p - dy3*dx3p

            ! Check signs
            if (all(cp > 0) .or. all(cp < 0)) then 
                in(i) = .true. 
            elseif (all(cp >= 0) .or. all(cp <= 0)) then 
                on(i) = .true.
            end if 

        end do 
        !$omp end parallel do

    end subroutine    

   ! InTriangleGrid routine
    subroutine InTriangleGrid(tria, xp, yp, ctri_array)

        ! Description
        !============
        ! This routine checks if a points lies within or on a triangle.
        ! The result is returned in ctri_arrays, containing the triangle number.
        ! The query point is given as an array set of coordinates xp, yp. 

        ! Declare variables
        !==================
        ! Arguments
        type(TriangulationUDT), intent(in)              :: tria
        real(R8), intent(in)                            :: xp(:), yp(:) 
        integer(I8), allocatable, intent(out)           :: ctri_array(:)

        ! Auxiliary
        real(R8)                        :: dx1p(tria%nc), dy1p(tria%nc), dx2p(tria%nc), &
            dy2p(tria%nc), dx3p(tria%nc), dy3p(tria%nc), cp(tria%nc,1:3), &
            x1(tria%nc), x2(tria%nc), x3(tria%nc), y1(tria%nc), y2(tria%nc), y3(tria%nc), &
            dx1(tria%nc), dx2(tria%nc), dx3(tria%nc), &
            dy1(tria%nc), dy2(tria%nc), dy3(tria%nc)
        logical                         :: do_parallel, in(tria%nc)
        
        ! Loop
        integer(I8)                         :: i, nc

        ! Initialize
        !===========
        nc = tria%nc
        in = .false. 
        allocate(ctri_array(size(xp,1)))
        ctri_array = 0

        ! Compute
        !========
        ! Loop
        do_parallel = .not. omp_in_parallel()

        x1 = tria%x(tria%cvert(:,1)) 
        x2 = tria%x(tria%cvert(:,2))
        x3 = tria%x(tria%cvert(:,3))
        y1 = tria%y(tria%cvert(:,1)) 
        y2 = tria%y(tria%cvert(:,2))
        y3 = tria%y(tria%cvert(:,3))
        dx1 = tria%dx(:,1)
        dx2 = tria%dx(:,2)
        dx3 = tria%dx(:,3)
        dy1 = tria%dy(:,1)
        dy2 = tria%dy(:,2)
        dy3 = tria%dy(:,3)

        !$omp parallel do default(shared) schedule(static) if (do_parallel) &
        !$omp private(i, dx1p, dx2p, dx3p, dy1p, dy2p, dy3p, cp)
        do i = 1, size(xp)

            ! Get vectors
            dx1p = xp(i) - x1
            dx2p = xp(i) - x2
            dx3p = xp(i) - x3
            dy1p = yp(i) - y1
            dy2p = yp(i) - y2
            dy3p = yp(i) - y3   
            
            ! Compute cross products
            cp(:,1) = dx1*dy1p - dy1*dx1p 
            cp(:,2) = dx2*dy2p - dy2*dx2p
            cp(:,3) = dx3*dy3p - dy3*dx3p

            in = all(cp >= 0, 2) .or. all(cp <= 0, 2) 
            ctri_array(i) = findloc(in, .true., 1) 
        end do   
        !$omp end parallel do

    end subroutine    

    ! Barycentric coordinate conversion
    subroutine Cart2Bary(x, y, x1, y1, x2, y2, x3, y3, lambda1, lambda2, lambda3)

        ! Description
        !============
        ! Convert cartesian coordinates (x, y) to barycentric 
        ! coordinates (lambda1, lambda2, lambda3) for the triangle
        ! given by points (1, 2, 3). Lambda1 then corresponds to the 
        ! weight of the first vertex, etc. Conversion is based on the
        ! following coordinate transformation (basic equations plus 
        ! normalization of lambdas):
        !
        ! x = lambda1*x1 + lambda2*x2 + lambda3*x3
        ! y = lambda1*y1 + lambda2*y2 + lambda3*y3
        ! lambda3 = 1 - lambda1 - lambda2
        ! 
        ! Inverting this system yields, with r = (x, y) etc:
        ! 
        ! lambda1 = (r - r3) x (r2 - r3)/(r1 - r3) x (r2 - r3)
        ! lambda2 = (r - r3) x (r3 - r1)/(r1 - r3) x (r2 - r3)
        ! lambda3 = 1 - lambda1 - lambda2

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)            :: x, y, x1, y1, x2, y2, x3, y3
        real(R8), intent(out)           :: lambda1, lambda2, lambda3

        ! Auxiliary
        real(R8)                        :: cp1, cp2, cp3

        ! Compute
        !========
        cp1 = (x - x3)*(y2 - y3) - (x2 - x3)*(y - y3)
        cp2 = (x - x3)*(y3 - y1) - (x3 - x1)*(y - y3)
        cp3 = (x1 - x3)*(y2 - y3) - (x2 - x3)*(y1 - y3)
        lambda1 = cp1/cp3
        lambda2 = cp2/cp3
        lambda3 = 1 - lambda1 - lambda2

    end subroutine

end module