!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!------------
! PolynomialInterpolatorStructured2D constructs a polynomial
! interpolant starting from the data points V = F(X,Y). V should be in a
! ny-by-nx format, and X (nx-by-1) and Y (ny-by-1) should be monotonically
! increasing coordinate vectors describing a (possibly non-uniform)
! rectangular grid on which V is given. The interpolant is described as
! follows in each quadrilateral of the grid:
!
!       f(x,y) = S_i=0^n S_j=0^n a_ij x^i y^j = V(i,j)
!
! Here, n is related to the desired continuity measure of the interpolant
! as n = 2*C + 1 (C >= 0, C is a natural number).
!
! In order to determine the required coefficients, one needs
! k = 4*n^2 equations, e.g.:
!
!       C       n       k       phi
!       0       1       4       1: phi
!       1       3       16      4: " and dx, dy, dxdy
!       2       5       36      9: " and dx2, dx2dy, dx2dy2, dy2dx, dy2
!       3       7       64      16: " and dx3, dx3dy, dx3dy2 dx3dy3, dy3,
!                                       dy3dx, dy3dx2
!
! These equations are provided by numerically evaluating the function and
! its spatial derivatives at the given grid points. Note that for C = 1,
! all derivatives containing x or y should be given (see column with
! 'phi'), for C = 2 all derivatives containing x^2 or y^2, ... These
! derivatives are computed numerically starting from the given sample
! points in V by constructing a polynomial representation of order M in
! each grid cell, where M >= C. A compact stencil is used to compute the
! required derivatives (for all derivatives, the same number of points is
! used, resulting in higher order approximations for lower order
! derivatives). Note that the given grid should be large enough to support
! the desired order of the approximation! Near boundaries, the stencil is
! shifted to accomodate this.

! Notes
!------
! Note 1: M is preferably an even number, such that the stencil (which is
! size M+1) is symmetric.

! Note 2: the vertex based interpolant to construct the derivatives can be
! done using a 'centered' interpolant (where each interpolant in a point is
! built relative to that coordinate, which is typically better scaled but
! slow) or by the 'fast' method (mathematically same accuracy, but
! numerically less well scaled, which may result in roundoff errors). Note
! that the 'fast' method suppresses any ill-conditioning warnings! 



module StructuredInterpolant2D

    ! Initialize
    !============
    ! Load modules
    use Interpolant2D_auxiliaries
    use Interpolant2D

    implicit none
    save

    ! Structured interpolant type
    !============================
    type, extends(GenericInterpolant2DUDT) :: StructuredInterpolant2DUDT 

        ! Description
        !============
        ! Apart from the fields of the generic interpolant, several 
        ! other fields are defined as well:
        ! - meth: methodology to construct interpolant (only 
        ! 'uniformgrid' currently supported)
        ! - C, M: order of the interpolant and order of the
        ! approximation method to compute derivatives for the 
        ! interpolant construction
        ! - xgv, ygv: grid vectors
        ! - A: interpolation coefficients
        ! - refx, refy, refdx, refdy: reference values used to compute
        ! derivatives etc 
        ! - cellindex: nx-1 by ny-1 array containing the cell indices
        ! - n: number of terms in the interpolant

        character(:), allocatable       :: meth 
        integer(I8)                     :: C, M, n

        real(R8), allocatable           :: xgv(:), ygv(:), A(:, :), &
            refx(:), refy(:), refdx(:), refdy(:)
        integer(I8), allocatable        :: cellindex(:, :)

    contains 

        ! Parameter setter routine
        procedure :: SetParameters

        ! Construct based on structured data
        procedure :: ConstructStructured => ConstructSI2DS 

        ! Construct based on unstructured data
        procedure :: ConstructUnstructured => ConstructSI2DUS

        ! Evaluator
        procedure :: Evaluate   => EvaluateStructuredInterpolant2D

        ! Deallocate
        procedure :: Deallocate => DeallocateStructuredInterpolant2D

    end type 

    contains 

    ! Set parameters
    subroutine SetParameters(interp, meth, C, M)

        ! Description
        !============
        ! Set the parameters of the interpolation routine

        ! Declare variables
        !==================
        ! Arguments
        class(StructuredInterpolant2DUDT)       :: interp
        character(:), allocatable, intent(in)   :: meth
        integer(I8), intent(in)                 :: C, M 

        ! Set
        !====
        interp%meth = meth 
        interp%C    = C 
        interp%M    = M
        
    end subroutine 

    ! Constructor, unstructured
    subroutine ConstructSI2DUS(interp, xg, yg, v)

        ! Declare variables
        !==================
        class(StructuredInterpolant2DUDT)       :: interp 
        real(R8), allocatable                   :: xg(:), yg(:), v(:)

        ! Currently, no implementation yet. Call error handler
        call gdErrorHandler('Unstructured initialization of 2D structured interpolant not yet implemented')

    end subroutine 

    ! Constructor, structured
    subroutine ConstructSI2DS(interp, xg, yg, v)

        ! Description
        !============
        ! This routine is a wrapper for the uniform and non-uniform 
        ! variants of the structured way of constructing the structured
        ! interpolant. 

        ! Declare variables
        !==================
        ! Arguments
        class(StructuredInterpolant2DUDT)       :: interp
        real(R8), allocatable                   :: xg(:), yg(:), v(:, :)

        ! Auxiliary

        ! Loop

        ! Initialize
        !===========
        ! Check which method to follow
        select case(interp%meth)

        case ('uniformgrid')

            call ConstructSI2DSuniform(interp, xg, yg, v)

        case ('nonuniformgrid')

            call ConstructSI2DSnonuniform(interp, xg, yg, v)

        case default 

            ! Unknown case, call error handler
            call gdErrorHandler('ConstructSI2DS: unknown construction method for 2D interpolant')

        end select 

    end subroutine
    
    ! Actual uniform constructor
    subroutine ConstructSI2DSuniform(interp, xgv, ygv, v)
        
        ! Description
        !============
        ! We build the interpolant based on the structured data by doing
        ! the following steps: 
        ! 0) Compute the required derivatives on the vertex nodes.
        ! 1) Normalize the coordinates, save normalization operation
        ! 2) Compute the derivatives of V using finite differences
        ! 3) Construct the interpolation matrix for a single quadrilateral, invert
        ! 4) Compute for each quadrilateral the coefficient values by
        ! multiplication
        ! 5) Construct the interpolant evaluation functions

        ! Here, we exploit the assumption that the grid is uniform in both
        ! directions (though the increments in X and Y may differ) by
        ! observing that in the bulk of the domain, the coefficient matrix
        ! for the vertex interpolant is the same for each bulk vertex. A
        ! vertex is in the bulk if the stencil doesn't have to be adapted.
        ! The boundaries and corners are treated separately. 

        ! Declare variables
        !==================
        ! Arguments
        class(StructuredInterpolant2DUDT)       :: interp 
        real(R8), allocatable                   :: xgv(:), ygv(:), &
            v(:, :)

        ! Auxiliary
        integer(I8)                             :: nx, ny, &
            ixlb, ixub, iylb, iyub, nxr, nyr, vID, termind, &
            nc, vecind 
        real(R8)                                :: dxmean, dymean, &
            reldevx, reldevy, prefac

        integer(I8), allocatable                :: stencil(:), &
            vertind(:, :), xstencil(:), ystencil(:), tvID(:), &
            xrange(:), yrange(:), dxstencil(:), dystencil(:), &
            cellindex(:, :)
        integer(I16), allocatable               :: allprefac(:)
        real(R8), allocatable                   :: xg(:, :), &
            yg(:, :), cij(:, :), A(:, :), dx(:), dy(:), res(:, :), &
            temp(:, :), fvals(:, :, :), refx(:), refy(:), refdx(:), &
            refdy(:), rhs(:, :), vec(:), aij(:, :)
            

        ! Loop
        integer(I8)                             :: i, j, k, l, ix, iy, &
            ii, jj, kk, ll, derivind, derivcount, nterms

        ! Linear solver variables (dummies basically)
        integer(I8), allocatable            :: ipiv(:)
        integer                             :: info
        integer                             :: neq

        ! Parameters
        real(R8), parameter                     :: tol = 1e-8

        ! Initialize
        !===========
        ! Unpack for ease
        associate(&
            M           => interp%M, &
            C           => interp%C &
            )

        ! Construct 2D grid
        nx = size(xgv, 1)
        ny = size(ygv, 1)
        allocate(xg(nx, ny), yg(nx, ny))
        xg = spread(xgv, 2, ny)
        yg = spread(ygv, 1, nx)

        ! Checks
        if ( (nx .ne. size(v, 1)) .or. (ny .ne. size(v, 2)) ) then
            ! Incompatible sizes
            call gdErrorHandler('ConstructSI2DSuniform: incompatible sizes, values should be nx-by-ny array')
        end if 
        if (M < C) then 
            ! Issue warning, set M = C
            print * , 'ConstructSI2DSuniform: M should be larger than or equal to C, setting M = C'
            M = C
        end if 
        if ( (nx <= M+1) .or. (ny <= M+1) ) then 
            ! Not enough points provided, throw error
            call gdErrorHandler('ConstructSI2DSuniform: the given grid is too small, it should be at least M+1-by-M+1')
        end if 

        ! Construct derivatives
        !======================
        ! Initialize coefficient matrix for derivatives (column-wise storage for
        ! nodes, coefficients as: c00, c10, ... cM0, c01, c11, .. cM1, ... cMM) 
        ! - this is different than in the matlab code! first x derivatives, 
        ! then y derivatives...

        ! Initialize stencil and matrices
        nterms = (M + 1)**2
        allocate(cij(nx*ny, nterms), A(nterms, nterms), res(nx*ny, nterms))
        cij(:, :) = 0

        ! Construct the base stencil
        allocate(stencil(M+1), xstencil(M+1), ystencil(M+1))
        if (mod(M, 2) == 0) then 
            ! Even
            do i = 0, M
                stencil(i+1) = i-M/2
            end do
        else 
            ! Uneven
            do i = 0, M
                stencil(i+1) = i-(M+1)/2
            end do
        end if

        ! Construct difference along xgv, ygv
        allocate(dx(nx-1), dy(ny-1))
        dx = xgv(2:nx) - xgv(1:nx-1)
        dy = ygv(2:ny) - ygv(1:ny-1)

        ! Check if the difference is approx. the same everywhere
        dxmean = sum(dx)/(nx-1) 
        dymean = sum(dy)/(ny-1)
        reldevx = maxval(abs( (dx - dxmean)/dxmean ))*100 ! in %
        reldevy = maxval(abs( (dy - dymean)/dymean ))*100

        if (reldevx >= tol) then 
            ! Issue warning, but this may be ok
            print *, 'ConstructSI2Duniform: non-uniformities detected' // &
                ' in X grid vector of approx. ', reldevx, ' %. ' // &
                ' Results may be inaccurate'
            
        end if 
        if (reldevy >= tol) then 
            ! Issue warning, but this may be ok
            print *, 'ConstructSI2Duniform: non-uniformities detected' // &
                ' in Y grid vector of approx. ', reldevx, ' %. ' // &
                ' Results may be inaccurate'
            
        end if 

        ! Set boundary indices
        ixlb = abs(stencil(1)) + 1
        ixub = nx - abs(stencil(M+1)) 
        iylb = abs(stencil(1)) + 1
        iyub = ny - abs(stencil(M+1)) 

        ! Initialize vertex indices
        allocate(vertind(nx, ny))
        vertind = reshape([(i, i = 1, nx*ny)], (/nx, ny/))

        ! Bulk of domain
        !---------------
        ! Loop to compute 'weights' of coefficients
        xstencil = stencil 
        ystencil = stencil 
        k = 1
        do j = 1, M+1
            do i = 1, M+1
                ! Compute coefficient matrix entry
                ! A(k, :) = reshape((spread(xstencil(i)**(/ix, ix=0,M/), 1, M+1)).*(repmat( (ystencil(j).^(0:M)), M+1, 1))', 1, []);
                l = 1
                do iy = 0, M
                    do ix = 0, M
                        A(k, l) = xstencil(i)**ix * ystencil(j)**iy
                        l = l + 1
                    end do 
                end do

                ! Update
                k = k + 1
            end do
        end do

        ! Loop to compute left hand side
        nxr = ixub - ixlb + 1
        nyr = iyub - iylb + 1
        allocate(xrange(nxr), yrange(nyr))
        xrange = [(i, i = ixlb, ixub)]
        yrange = [(i, i = iylb, iyub)]

        do i = 1, nxr
            ! Compute stencil
            ix = xrange(i)
            xstencil = stencil + ix

            do j = 1, nyr 
                ! Compute stencil
                iy = yrange(j)
                ystencil = stencil + iy 

                ! Compute left hand side
                vID = vertind(ix, iy)
                res(vID, :) = reshape(v(xstencil, ystencil), (/nterms/))

            end do 
        end do 

        ! Solve for internal vertices only
        allocate(tvID(nxr*nyr), temp(nterms, nxr*nyr), ipiv(nterms))
        tvID = reshape(vertind(xrange, yrange), (/nxr*nyr/))
        temp = transpose(res(tvID, :))
        call dgesv(nterms, nxr*nyr, A, nterms, ipiv, temp, nterms, info)
        cij(tvID, :) = transpose(temp(:, :))
        deallocate(tvID, temp, ipiv)

        ! Top & bottom boundary
        !----------------------
        ! For each iy we need to construct a new interpolation matrix, but
        ! for ix it remains the same
        deallocate(yrange)
        nyr = iylb + ny-iyub+1
        allocate(yrange(nyr), dystencil(nyr))
        xrange = [(i, i = ixlb, ixub)]
        yrange = [ [(i, i = 1, iylb)], [(i, i = iyub, ny)] ]
        dystencil = [ [(i, i = iylb-1, 0, -1)], -[(i, i = 0, ny-iyub, 1)] ] 
        xstencil = stencil

        ! Loop
        ll = 1
        do j = 1, nyr
            ! Compute stencil
            ystencil = stencil + dystencil(ll) ! adjust for boundary

            ! Compute coefficients
            k = 1
            do jj = 1, M+1
                do ii = 1, M+1
                    
                    ! Compute coefficient matrix entry (c00, c01, ... , c0M,
                    ! c11, ...)
                    l = 1
                    do iy = 0, M
                        do ix = 0, M
                            A(k, l) = xstencil(ii)**ix * ystencil(jj)**iy
                            l = l + 1
                        end do 
                    end do
                    
                    ! Update k
                    k = k + 1; 
                end do
            end do

            ! Loop over bulk vertices to compute weights
            do i = 1, nxr 
                ! Get residuals
                vID = vertind(xrange(i), yrange(j))
                res(vID, :) = reshape(v(xstencil + xrange(i), ystencil + yrange(j)), (/nterms/))
            end do 

            ! Solve for boundary vertices
            allocate(tvID(nxr), temp(nterms, nxr), ipiv(nterms))
            tvID = reshape(vertind(xrange, yrange(j)), (/nxr/))
            temp = transpose(res(tvID, :))
            call dgesv(nterms, nxr, A, nterms, ipiv, temp, nterms, info)
            cij(tvID, :) = transpose(temp(:, :))
            deallocate(tvID, temp, ipiv)

            ! Update ll
            ll = ll + 1
        end do 

        ! Deallocate
        deallocate(dystencil)

        ! Left & right boundary
        !----------------------
        ! For each iy we need to construct a new interpolation matrix, but
        ! for ix it remains the same
        deallocate(xrange, yrange)
        nxr = ixlb + nx-ixub+1
        nyr = iyub - iylb + 1
        allocate(xrange(nxr), yrange(nyr), dxstencil(nxr))
        yrange = [(i, i = iylb, iyub)]
        xrange = [ [(i, i = 1, ixlb)], [(i, i = ixub, nx)] ]
        dxstencil = [ [(i, i = ixlb-1, 0, -1)], -[(i, i = 0, nx-ixub)] ] 
        ystencil = stencil

        ! Loop
        ll = 1
        do i = 1, nxr
            ! Compute stencil
            xstencil = stencil + dxstencil(ll) ! adjust for boundary

            ! Compute coefficients
            k = 1
            do jj = 1, M+1
                do ii = 1, M+1
                    
                    ! Compute coefficient matrix entry (c00, c01, ... , c0M,
                    ! c11, ...)
                    l = 1
                    do iy = 0, M
                        do ix = 0, M
                            A(k, l) = xstencil(ii)**ix * ystencil(jj)**iy
                            l = l + 1
                        end do 
                    end do
                    
                    ! Update k
                    k = k + 1; 
                end do
            end do

            ! Loop over bulk vertices to compute weights
            do j = 1, nyr 
                ! Get residuals
                vID = vertind(xrange(i), yrange(j))
                res(vID, :) = reshape(v(xstencil + xrange(i), ystencil + yrange(j)), (/nterms/))
            end do 

            ! Solve for boundary vertices
            allocate(tvID(nyr), temp(nterms, nyr), ipiv(nterms))
            tvID = reshape(vertind(xrange(i), yrange), (/nyr/))
            temp = transpose(res(tvID, :))
            call dgesv(nterms, nyr, A, nterms, ipiv, temp, nterms, info)
            cij(tvID, :) = transpose(temp(:, :))
            deallocate(tvID, temp, ipiv)

            ! Update ll
            ll = ll + 1

        end do

        ! Housekeeping
        deallocate(dxstencil)

        ! Corners
        !--------
        ! Here, the stencil is different for each 'corner' vertex, and we
        ! have to resort to the general method.
        deallocate(xrange, yrange)
        nxr = ixlb + nx-ixub+1
        nyr = iylb + ny-iyub+1
        allocate(xrange(nxr), yrange(nyr), dxstencil(nxr), dystencil(nyr))
        xrange = [ [(i, i = 1, ixlb)], [(i, i = ixub, nx)] ]
        dxstencil = [ [(i, i = ixlb-1, 0, -1)], -[(i, i = 0, nx-ixub)] ] 
        yrange = [ [(i, i = 1, iylb)], [(i, i = iyub, ny)] ]
        dystencil = [ [(i, i = iylb-1, 0, -1)], -[(i, i = 0, ny-iyub)] ] 

        ! Loop
        ll = 1;
        do i = 1, nxr
            ! Stencils
            xstencil = stencil + dxstencil(ll);
            
            kk = 1;
            do j = 1, nyr
                ! Stencils
                ystencil = stencil + dystencil(kk);

                ! Compute coefficients
                k = 1;
                do jj = 1, M+1
                    do ii = 1, M+1
                        l = 1
                        do iy = 0, M
                            do ix = 0, M
                                A(k, l) = xstencil(ii)**ix * ystencil(jj)**iy
                                l = l + 1
                            end do 
                        end do

                        ! Update k
                        k = k + 1;
                    end do
                end do
                
                ! Get residuals
                vID = vertind(xrange(i), yrange(j));
                res(vID, :) = reshape(v(xstencil + xrange(i), ystencil + yrange(j)), (/nterms/));
                    
                ! Solve for boundary vertices only
                allocate(tvID(1), temp(nterms, 1), ipiv(nterms))
                tvID = vertind(xrange(i), yrange(j))
                temp = transpose(res(tvID, :))
                call dgesv(nterms, 1, A, nterms, ipiv, temp, nterms, info)
                cij(tvID, :) = transpose(temp(:, :))
                deallocate(tvID, temp, ipiv)
                
                ! Update kk
                kk = kk + 1;
                
            end do
            
            ! Update ll 
            ll = ll + 1;
        end do

        ! Compute actual derivatives
        !---------------------------
        ! Correct for scaling
        k = 1
        do j = 0, M
            do i = 0, M 
                cij(:, k) = cij(:, k)/( (dxmean**i) * (dymean**j))
                k = k + 1
            end do 
        end do

        ! Construct derivatives as fields - need to account for 
        ! factorial prefactor!
        allocate(allprefac((M+1)**2), fvals(size(v, 1), size(v, 2), nterms))
        k = 1
        do j = 0, M
            do i  = 0, M 
                allprefac(k) = int(MyFactorial(i)*MyFactorial(j), kind=I16)
                k = k + 1
            end do 
        end do 
        fvals(:, :, :) = 0
        do k = 1, nterms
            fvals(:, :, k) = allprefac(k)*reshape(cij(:, k), (/nx,  ny/))
        end do 

        ! Construct interpolation matrix
        !===============================
        ! Note: this matrix is determined as [f fx fy fxy]' = C*a, where a are the
        ! coefficients a_ij (in this case, there are 16 coefficients), and C is the
        ! linearization matrix w.r.t a_ij. In normalized coordinates, there are
        ! four points: (0,0) (0,1) (1,1) (1,0). Since at these normalized
        ! coordinates the x and y values (including their powers) are either 0 or
        ! 1, the computation of C simplifies substantially.

        ! The derivative of the general polynomial
        !
        !       f = sum_i=0..n, sum_j=0..n a_ij x^i + y^j
        !
        ! is equal to:
        !
        !       d^(nx+ny)f/dx^(nx)dy^(ny) = sum_i=nx..n sum_j=ny..n ...
        !           a_ij i!/(i-nx)! j!/(j-ny)! x^(i-nx) y(j-ny)
        !
        ! where nx is the nx-th partial derivative in the x-direction (ny
        ! analogous). x and y should be substituted by the coordinate points.

        ! It should be noted that at the points where the function/derivative is
        ! evaluated, that x and y are zero in the representation of phi at the
        ! nodes. This means that at the node, f = a00

        ! Initialize
        nterms = (2*(C+1))**2
        deallocate(A)
        allocate(A(nterms, nterms))
        A(:, :) = 0 

        ! Loop over all derivatives
        derivind = 0
        do j = 0, C 
            do i = 0, C
                ! Construct coefficients of the derivative of this 
                ! polynomial
                do l = j, 2*(C+1)-1
                    do k = i, 2*(C+1)-1
                        ! Factorial prefactor
                        prefac = real(MyFactorial(l), kind=R8)/real(MyFactorial(l-j), kind = R8)*&
                            real(MyFactorial(k), kind=R8)/real(MyFactorial(k-i), kind = R8)
                        
                        ! Index of the term
                        termind = l*(2*(C+1)) + k + 1
        
                        ! Contribution of each point
                        A(1+4*derivind, termind) = prefac*0**(k - i)*0**(l - j);
                        A(2+4*derivind, termind) = prefac*1**(k - i)*0**(l - j);
                        A(3+4*derivind, termind) = prefac*1**(k - i)*1**(l - j);
                        A(4+4*derivind, termind) = prefac*0**(k - i)*1**(l - j);
                    end do
                end do

                ! Update counter
                derivind = derivind + 1
            end do
        end do

        ! Compute coefficients
        !=====================
        ! Initialize
        nc = (nx-1)*(ny-1) ! number of cells
        allocate(cellindex(nx-1, ny-1))
        allocate(aij(nc, nterms), refx(nc), refy(nc), refdx(nc), &
            refdy(nc), rhs(nc, nterms), vec((2*(C+1))**2) )
        cellindex = reshape([(i, i = 1, (nx-1)*(ny-1))], (/nx-1, ny-1/))

        ! Loop over all quads
        do i = 1, nx-1
            do j = 1, ny-1
                ! Value vectors at four points
                vec = 0
                do l = 0, C
                    do k = 0, C
                        ! Set indices
                       termind = l*(M + 1) + k
                        vecind = l*(C + 1) + k

                        ! Compute
                        vec(1+vecind*4) = fvals(i, j, termind+1);
                        vec(2+vecind*4) = fvals(i+1, j, termind+1);
                        vec(3+vecind*4) = fvals(i+1, j+1, termind+1);
                        vec(4+vecind*4) = fvals(i, j+1, termind+1);
                    end do 
                end do 

                ! Save residual
                rhs(cellindex(i, j), :) = vec 

                ! Store scaling factors & reference values
                refx(cellindex(i, j))   = xg(i, j)
                refy(cellindex(i, j))   = yg(i, j)
                refdx(cellindex(i, j))  = dxmean 
                refdy(cellindex(i, j))  = dymean
            end do
        end do

        ! Scale 
        derivcount  = 0
        do l = 0, C
            do k = 0, C 
                rhs(:, derivcount+1:derivcount+4) = &
                    rhs(:, derivcount+1:derivcount+4)*(dxmean**k)*(dymean**l)
                derivcount = derivcount + 4
            end do 
        end do 

        ! Solve for coefficients
        allocate(ipiv(nterms), temp(size(rhs, 2), size(rhs, 1)))
        temp = transpose(rhs) 
        call dgesv(nterms, nc, A, nterms, ipiv, temp, nterms, info)
        rhs = transpose(temp)
        aij(reshape(cellindex, (/nc/)), :) = rhs(reshape(cellindex, (/nc/)), :)
        deallocate(ipiv)


        ! Add to interpolant
        !===================
        ! Allocate first
        allocate(interp%xgv(nx), interp%ygv(ny), &
            interp%A(size(aij, 1), size(aij, 2)), interp%refx(nx-1), &
            interp%refy(ny-1), interp%refdx(nx-1), interp%refdy(ny-1), &
            interp%cellindex(nx-1, ny-1))
        
        

        ! Housekeeping
        !=============
        end associate

        ! Attribute
        interp%xgv      = xgv
        interp%ygv      = ygv 
        interp%n        = 2*interp%C+1
        interp%A        = aij 
        interp%refx     = refx
        interp%refy     = refy 
        interp%refdx    = refdx
        interp%refdy    = refdy 
        interp%cellindex    = cellindex 

        deallocate(stencil, vertind, xstencil, ystencil, xrange, &
            yrange, dxstencil, dystencil, allprefac, cellindex, &
            xg, yg, cij, A, dx, dy, res, fvals, refx, refy, &
            refdx, refdy, rhs, vec, aij)
        
        

    end subroutine 

    ! Actual nonuniform constructor
    subroutine ConstructSI2DSnonuniform(interp, xg, yg, v)

        ! Description
        !============

        ! Declare variables
        !==================
        ! Arguments
        class(StructuredInterpolant2DUDT)       :: interp 
        real(R8), allocatable                   :: xg(:), yg(:), v(:, :)

        ! Auxiliary

    end subroutine 

    ! Structured evaluator
    subroutine EvaluateStructuredInterpolant2D(interp, xq, yq, derivx, derivy, vq)

        ! Description
        !============
        ! Evaluate the interpolation operator based on the coefficients per
        ! cell (quadrilateral). First, it is tested which derivative has to
        ! be computed (derivx and derivy give the orders of the derivative,
        ! if both 0, then the function itself is evaluated). Afterwards,
        ! the index of the quadrilateral where to evaluate xq and yq is
        ! determined. The coefficient values are stored as f, fy, fyy, ...,
        ! fx, fxy, fxyy, ... etc
        
        ! The derivative of the general polynomial
        !
        !       f = sum_i=0..n, sum_j=0..n a_ij x^i + y^j
        !
        ! is equal to:
        !
        !       d^(nx+ny)f/dx^(nx)dy^(ny) = sum_i=nx..n sum_j=ny..n ...
        !           1/refx^nx 1/refy^ny a_ij i!/(i-nx)! j!/(j-ny)! x^(i-nx) y(j-ny)
        !
        ! where nx is the nx-th partial derivative in the x-direction (ny
        ! analogous). x and y should be substituted by the normalized
        ! coordinates, refx and refy are the scaling factors with which the
        ! coordinates were scaled to the unity square

        ! Declare variables
        !==================
        ! Arguments
        class(StructuredInterpolant2DUDT)       :: interp 
        real(R8), allocatable, intent(in)       :: xq(:), yq(:)
        real(R8), allocatable, intent(out)      :: vq(:)
        integer(I8), intent(in)                 :: derivx, derivy

        ! Auxiliary
        integer(I8)                             :: nq
        real(R8)                                :: prefac

        integer(I8), allocatable                :: ind(:) 
        real(R8), allocatable                   :: xqn(:), yqn(:), &
            term(:), thisA(:, :)

        ! Loop
        integer(I8)                             :: i, j, indder

        ! Initialize
        !===========
        ! Check inputs
        if ( (.not. allocated(xq)) .or. (.not. allocated(yq)) ) then 
            call gdErrorHandler('EvaluateStructuredInterpolant2D: query points are not allocated')
        end if 
        if (size(xq, 1) .ne. size(yq, 1)) then 
            call gdErrorHandler('EvaluateStructuredInterpolant2D: ' // &
            'query point coordinates xq and yq have dissimilar ' // &
            'dimensions, check input')
        end if 
        if (allocated(vq)) then 
            ! Deallocate just to be sure
            deallocate(vq)
        end if 

        ! Initialize
        nq = size(xq, 1)

        ! Associate
        associate( &
            xv      => interp%xgv, &
            yv      => interp%ygv, &
            n       => interp%n, &
            A       => interp%A, & 
            refx    => interp%refx, &
            refy    => interp%refy, &
            refdx   => interp%refdx, &
            refdy   => interp%refdy, &
            cellindex => interp%cellindex)

        ! Get cell index of query point
        call GetIndex(xq, yq, xv, yv, ind)

        ! Extract
        allocate(thisA(size(ind, 1), size(A, 2)))
        thisA = A(ind, :)

        ! Compute normalized query points
        allocate(xqn(nq), yqn(nq))
        xqn = (xq - refx(ind))/refdx(ind)
        yqn = (yq - refy(ind))/refdy(ind)

        ! Evaluate and sum each term
        allocate(vq(nq), term(nq))
        vq(:) = 0

        do i = derivx, n 
            do j = derivy, n
                ! Derivative index
                indder = j*(n+1) + i + 1

                ! Factorial prefactor
                prefac = real(MyFactorial(i), kind=R8)/&
                    real(MyFactorial(i-derivx), kind = R8)*&
                    real(MyFactorial(j), kind=R8)/real(MyFactorial(j-derivy), kind = R8)
                term = thisA(:, indder)*prefac*xqn**(i - derivx)*yqn**(j - derivy)
                vq = vq + term 
            end do 
        end do

        ! Scale
        vq = vq/( (refdx(ind))**derivx * (refdy(ind))**derivy)

        ! Housekeeping
        !=============
        deallocate(xqn, yqn, term, ind)
        end associate

    end subroutine

    ! Structured deallocator
    subroutine DeallocateStructuredInterpolant2D(interp)

        ! Description
        !============
        ! Deallocate a fully constructued interpolant 
        class(StructuredInterpolant2DUDT)       :: interp 

        if (allocated(interp%cellindex)) then 
            ! Assume all allocated
            deallocate(interp%cellindex, interp%refdx, interp%refdy, &
                interp%refx, interp%refy, interp%xgv, interp%ygv, &
                interp%A)
        end if 
        
    end subroutine
    

end module