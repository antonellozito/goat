!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module provides all the routines that can operate on a bicubic 
! spline interpolant structure (or object), of which the type is defined
! in gdmod_types. The interpolant can be initialized, constructed, and 
! operated upon using the routines given here. However, for generic 
! programming reasons, the routines given here are private, as the 
! actual operation should go through the interpolation module, which 
! serves as an interface where the correct interpolation routines are 
! specified, depending on the interpolation input. 

! Note that the interpolant is only implemented for (non-uniform) 
! Cartesian grids. 

! Mathematical description
!=========================
! The interpolant is described as  follows in each quadrilateral of 
! the grid:
!
!       f(x,y) = S_i=0^3 S_j=0^3 a_ij x^i y^j = V(i,j)
!       fx(x,y) = S_i=1^3 S_j=0^3 a_ij i x^(i-1) y^j = Vx(i,j)
!       fy(x,y) = S_i=0^3 S_j=1^3 a_ij j x^i y^(j-1) = Vy(i,j)
!       fxy(x,y) = S_i=1^3 S_j=1^3 a_ij ij x^(i-1) y^(j-1) = Vxy(i,j)
!
! Note that fx is the derivative of f w.r.t. x etc. The derivatives of the
! nodal values are computed using finite differences, e.g. using the matlab
! gradient routine. Higher order derivatives can afterwards be computed by
! constructing the analytic differences per cell.

! Algorithm
!==========
! 1) Normalize the coordinates, save normalization operation
! 2) Compute the derivatives of V using finite differences
! 3) Construct the interpolation matrix for a single quadrilateral, invert
! 4) Compute for each quadrilateral the coefficient values by
! multiplication
! 5) Construct the interpolant evaluation functions

module BicubicSplineInterpolant

    ! Initialize
    !============
    ! Load modules
    use gdmod_types
    use gdmod_plots

    ! The usual
    implicit none
    save
    public

    external sgesv
    external findloc
    !private
    
    !public ConstructBicubicSplineInterpolant

    contains 

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    ! Constructor
    !============
    subroutine ConstructBicubicSplineInterpolant(v, x, y, nx, ny, interp)

        ! Description
        !============
        ! Interpolant constructor. v should be a nx-by-ny array 
        ! containing the values on the grid points given by the grid 
        ! vectors x, y, i.e. v(i,j) is the value on x(i), y(j). nx and 
        ! ny should be the number of points in x and y. 

        ! Declare variables
        !==================
        ! Arguments
        type(BicubicSplineInterpolantUDT)   :: interp
        real(R8), intent(in)                :: x(:), y(:), v(:,:)
        integer(I8), intent(in)                :: nx, ny

        ! Loop variables
        integer(I8)                 :: i, j, ind

        ! Auxiliary variables 
        double precision                  :: xg(1:nx, 1:ny), yg(1:nx, 1:ny), & 
            x0(1:nx-1, 1:ny-1), y0(1:nx-1, 1:ny-1), vx(1:nx, 1:ny), &
            vy(1:nx, 1:ny), vxy(1:nx, 1:ny), vyx(1:nx, 1:ny), &
            vxx(1:nx, 1:ny), vyy(1:nx, 1:ny), dx(1:nx-1), dy(1:ny-1), &
            tdx, tdy       

        double precision                    :: c(1:16,1:16)
        double precision, allocatable       :: rhs(:,:)

        ! Linear solver variables (dummies basically)
        integer(I8), allocatable           :: ipiv(:)
        integer                     :: info
        integer                     :: neq


        ! Initialize
        !===========
        neq = 16

        ! First, allocate (routine specified in gdmod_types)
        interp%nx = nx
        interp%ny = ny
        call AllocateBicubicSplineInterpolant(interp)

        ! Construct the 2D grid
        xg = spread(x, 2, ny)
        yg = spread(y, 1, nx)

        ! Construct the reference coordinate for each cell
        x0 = xg(1:nx-1, 1:ny-1)
        y0 = yg(1:nx-1, 1:ny-1)

        ! Construct spacings (for later)
        dx = x(2:nx) - x(1:nx-1)
        dy = y(2:ny) - y(1:ny-1)

        ! Compute data gradients
        !=======================
        ! vx, vy
        call ComputeGradient2DStructured(v, xg, yg, nx, ny, vx, vy)

        ! Cross derivative (take average)
        call ComputeGradient2DStructured(vx, xg, yg, nx, ny, vxx, vxy)
        call ComputeGradient2DStructured(vy, xg, yg, nx, ny, vyx, vyy)
        vxy = (vxy + vyx)*0.5

        !call Plot2DStructuredField(v, x, y, nx, ny, '-p')
        !call Plot2DStructuredField(vx, x, y, nx, ny, '-p')
        !call Plot2DStructuredField(vy, x, y, nx, ny, '-p')
        !call Plot2DStructuredField(vxy, x, y, nx, ny, '-p')

        ! Construct the interpolation matrix
        !===================================
        ! Note: this matrix is determined as [f fx fy fxy]' = C*a, 
        ! where a are the coefficients a_ij (in this case, there are 
        ! 16 coefficients), and C is the linearization matrix w.r.t 
        ! a_ij. In normalized coordinates, there are four points: 
        ! (0,0) (0,1) (1,1) (1,0). Since at these normalized
        ! coordinates the x and y values (including their powers) are 
        ! either 0 or 1, the computation of C simplifies substantially.

        ! Loop
        ind = 1
        c(:,:) = 0
        do i = 0, 3
            do j = 0, 3
                ! f
                c(1,ind) = 0**(i)*0**(j)
                c(2,ind) = 1*0**(j)
                c(3,ind) = 1
                c(4,ind) = 0**(i)*1

                ! fx
                if (i > 0) then
                    c(5,ind) = i*0**(i-1)*0**(j)
                    c(6,ind) = i*1*0**(j)
                    c(7,ind) = i
                    c(8,ind) = i*0**(i-1)*1
                end if

                ! fy
                if (j > 0) then
                    c(9,ind) = j*0**(i)*0**(j-1)
                    c(10,ind) = j*1*0**(j-1)
                    c(11,ind) = j
                    c(12,ind) = j*0**(i)*1
                end if

                ! fxy
                if ((i > 0) .and. (j > 0)) then
                    c(13,ind) = i*j*0**(i-1)*0**(j-1)
                    c(14,ind) = i*j*0**(j-1)
                    c(15,ind) = i*j
                    c(16,ind) = i*j*0**(i-1)*1
                end if

                ! Update index
                ind = ind+1

            end do
        end do 

        ! Compute coefficients
        !=====================
        ! Loop over all quadrilaterals to compute the rhs
        allocate(rhs(16, interp%nc))
        ind = 1
        do j = 1, ny-1
            ! Extract scaling factor
            tdy = dy(j)
            do i = 1, nx-1
                ! Extract scaling factor
                tdx = dx(i)

                ! Add value vector
                rhs(:,ind) = &
                    (/v(i, j), v(i+1, j), v(i+1, j+1), v(i, j+1), &
                    tdx*(/vx(i, j), vx(i+1, j), vx(i+1, j+1), vx(i, j+1)/), &
                    tdy*(/vy(i, j), vy(i+1, j), vy(i+1, j+1), vy(i, j+1)/), &
                    tdx*tdy*(/vxy(i, j), vxy(i+1, j), vxy(i+1, j+1), vxy(i, j+1)/) &
                    /)

                ! Store scaling factors and reference values
                interp%refx(ind) = x0(i, j)
                interp%refy(ind) = y0(i, j)
                interp%refdx(ind) = tdx
                interp%refdy(ind) = tdy

                ! Update index
                ind = ind+1
            end do
        end do

        ! Call the solver (make sure input is in right format!)
        allocate(ipiv(neq))
        call dgesv(neq, interp%nc, c, neq, ipiv, rhs, neq, info)

        ! Check if converged
        if (info .ne. 0) then
            ! Not converged - issue error
            call gdErrorHandler('Could not compute interpolant, check input')
        end if

        ! Add interpolant fields
        !=======================
        interp%x = x
        interp%y = y
        interp%a = transpose(rhs)

        ! Housekeeping
        deallocate(rhs, ipiv)

    end subroutine

    ! Evaluator
    !==========
    subroutine EvaluateBicubicSplineInterpolant(xq, yq, vq, interp, &
        derivx, derivy)

        ! Description
        !============
        ! Evaluator for the interpolant. xq and yq should be the query
        ! points (equal size, column vector format), interp the bicubic
        ! spline interpolant object with all necessary fields 
        ! added (i.e. it has to be constructed), and derivx and derivy 
        ! are characters denoting which derivative should be taken in x
        ! and y direction. vq is the value at the query points (output)

        ! Declare variables
        !==================
        ! Arguments
        type(BicubicSplineInterpolantUDT), intent(in)   :: interp
        real(R8), intent(in)                :: xq(:), yq(:)
        character(1), intent(in)            :: derivx, derivy
        real(R8), intent(out)               :: vq(:)

        ! Loop variables
        integer(I8)                         :: i, j, k

        ! Auxiliary variables 
        integer(I8)                         :: nq 
        integer(I8), allocatable            :: ind(:)
        character(2)                        :: deriv

        real(R8), allocatable               :: aq(:,:), xqn(:), & 
            yqn(:), rep(:,:), sf(:)

        ! Initialize
        !===========
        ! Derive the number of query points
        nq = size(xq, 1)

        ! Concatenate the derivative
        deriv = derivx // derivy

        ! Allocate
        allocate(ind(nq))

        ! Get bins
        !=========
        ! Get the index of which 'bin' the query points belong to (i.e. 
        ! which cell do the points lie in). The first 1:nx bins are 
        ! located at y(1), the second at y(2) etc (so column order)
        call GetIndex(xq, yq, interp%x, interp%y, ind)

        ! Check if any indices are equal to zero, if so: throw error
        if (any(ind == 0)) then
            call gdErrorHandler('EvaluateBicubicSplineInterpolant: points lie outside of the domain, not supported')
        end if

        ! Evaluate
        !=========
        ! Allocate
        allocate(aq(nq, 16))
        allocate(xqn(nq))
        allocate(yqn(nq))
        allocate(rep(nq, 16))
        allocate(sf(nq))

        ! Extract a values
        aq = interp%a(ind,:)

        ! Compute normalized query points
        xqn = (xq - interp%refx(ind))/interp%refdx(ind)
        yqn = (yq - interp%refy(ind))/interp%refdy(ind)

        ! Evaluate
        select case (deriv) 

        case ('00')

            ! Function evaluation

            ! Construct polynomial representation
            k = 1
            do i = 0, 3
                do j = 0, 3
                    rep(:,k) = (xqn**i)*(yqn**j)
                    k = k + 1 
                end do
            end do

        case ('10')

            ! First order derivative w.r.t. x

            ! Construct polynomial representation
            sf = 1/interp%refdx(ind) ! rescale due to normalization
            k = 1
            do i = 0, 3
                do j = 0, 3
                    if (i > 0) then
                        rep(:,k) = sf*(i*xqn**(i-1))*(yqn**j)
                    else
                        rep(:,k) = 0
                    end if
                    k = k+1
                end do
            end do

        case ('01')

            ! First order derivative w.r.t. y

            ! Construct polynomial representation
            sf = 1/interp%refdy(ind) ! rescale due to normalization
            k = 1
            do i = 0, 3
                do j = 0, 3
                    if (j > 0) then
                        rep(:,k) = sf*(xqn**(i))*(j*yqn**(j-1))
                    else
                        rep(:,k) = 0
                    end if
                    k = k+1
                end do
            end do

        case ('11')

            ! Second order derivative w.r.t. x, y

            ! Construct polynomial representation
            sf = 1/(interp%refdy(ind)*interp%refdx(ind)) 
            k = 1
            do i = 0, 3
                do j = 0, 3
                    if ((j > 0) .and. (i > 0)) then
                        rep(:,k) = sf*(i*xqn**(i-1))*(j*yqn**(j-1))
                    else
                        rep(:,k) = 0
                    end if
                    k = k+1
                end do
            end do

        case ('20')

            ! Second order derivative w.r.t. x, x

            ! Construct polynomial representation
            sf = 1/(interp%refdx(ind)**2) 
            k = 1
            do i = 0, 3
                do j = 0, 3
                    if ((i > 1)) then
                        rep(:,k) = sf*(i*(i-1)*xqn**(i-2))*(yqn**(j))
                    else
                        rep(:,k) = 0
                    end if
                    k = k+1
                end do
            end do

        case ('02')

            ! First order derivative w.r.t. y, y

            ! Construct polynomial representation
            sf = 1/(interp%refdy(ind)**2) 
            k = 1
            do i = 0, 3
                do j = 0, 3
                    if ((j > 1)) then
                        rep(:,k) = sf*(xqn**(i))*(j*(j-1)*yqn**(j-2))
                    else
                        rep(:,k) = 0
                    end if
                    k = k+1
                end do
            end do

        case ('21')

            ! Third order derivative w.r.t. x, x, y

            ! Construct polynomial representation
            sf = 1/(interp%refdx(ind)**2*interp%refdy(ind)) 
            k = 1
            do i = 0, 3
                do j = 0, 3
                    if ((i > 1) .and. (j > 0)) then
                        rep(:,k) = sf*(i*(i-1)*xqn**(i-2))*(j*yqn**(j-1))
                    else
                        rep(:,k) = 0
                    end if
                    k = k+1
                end do
            end do

        case ('12')

            ! Third order derivative w.r.t. x, y, y

            ! Construct polynomial representation
            sf = 1/(interp%refdx(ind)*interp%refdy(ind)**2) 
            k = 1
            do i = 0, 3
                do j = 0, 3
                    if ((i > 0) .and. (j > 1)) then
                        rep(:,k) = sf*(i*xqn**(i-1))*(j*(j-1)*yqn**(j-2))
                    else
                        rep(:,k) = 0
                    end if
                    k = k+1
                end do
            end do

        case ('30')

            ! Third order derivative w.r.t. x, x, x

            ! Construct polynomial representation
            sf = 1/(interp%refdx(ind)**3) 
            k = 1
            do i = 0, 3
                do j = 0, 3
                    if ((i > 2)) then
                        rep(:,k) = sf*(i*(i-1)*(i-2)*xqn**(i-3))*(yqn**(j))
                    else
                        rep(:,k) = 0
                    end if
                    k = k+1
                end do
            end do

        case ('03')

            ! Third order derivative w.r.t. y, y, y

            ! Construct polynomial representation
            sf = 1/(interp%refdy(ind)**3) 
            k = 1
            do i = 0, 3
                do j = 0, 3
                    if ((j > 2)) then
                        rep(:,k) = sf*(xqn**(i))*(j*(j-1)*(j-2)*yqn**(j-3))
                    else
                        rep(:,k) = 0
                    end if
                    k = k+1
                end do
            end do

        case default 

            ! Throw error
            call gdErrorHandler('EvaluateBicubicSplineInterpolation: unknown derivative case')

        end select

        ! Evaluate
        vq = sum(aq*rep, 2)

        ! Housekeeping
        deallocate(aq, rep, ind, xqn, yqn, sf)

    end subroutine

    ! Binner (index retriever)
    subroutine GetIndex(xq, yq, x, y, ind)

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
        real(R8), intent(in)                :: xq(:), yq(:), x(:), y(:)
        integer(I8), intent(out)            :: ind(:)

        ! Loop variables
        integer(I8)                         :: k, indx, indy

        ! Auxiliary variables 
        integer(I8)                         :: nq, nx, ny

        ! Compute indices
        !================
        ! Compute sizes
        nq = size(xq, 1)
        nx = size(x, 1)
        ny = size(y, 1)

        ! Loop over all query points
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

    end subroutine

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
        logical                         :: notfound

        ! Compute the bin
        !================
        ! Determine the position in x (not optimized)
        ind = 1
        notfound = .true. 

        ! Check if the value lies inside the range
        if ((xq < (minval(x))) .or. (xq > maxval(x))) then
            ! Outside of bin
            notfound = .false. 
            ind = 0
        end if

        ! Check the first bin - equality holds for first value
        if (notfound) then
            if ((xq >= x(ind)) .and. (xq <= x(ind+1))) then
                notfound = .false.
                ind = 1
            else
                ! Not found, increase index
                ind = ind+1
            end if
        end if

        ! Loop over the remaining bins
        do while ((notfound) .and. (ind < nx))
            if ((xq > x(ind)) .and. (xq <= x(ind+1))) then
                ! Found, exit the loop
                notfound = .false.
            else 
                ! Update ind
                ind = ind+1
            end if
        end do


    end subroutine


    

end module