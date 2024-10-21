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
    use mod_constants
    use mod_sparseinterface
    use mod_linearsolverinterface
    use omp_lib

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
        ! - precomputedfac: the required factorials precomputed to save
        ! - some time during evaluation

        character(:), allocatable       :: meth 
        integer(I8)                     :: C, M, n

        real(R8), allocatable           :: xgv(:), ygv(:), A(:, :), &
            refx(:), refy(:), refdx(:), refdy(:)
        integer(I8), allocatable        :: cellindex(:, :)
        integer(I16), allocatable       :: precomputedfac(:)

    contains 

        ! Parameter setter routine
        procedure :: SetParameters

        ! Construct based on structured data
        procedure :: ConstructStructured => ConstructSI2DS 

        ! Construct based on unstructured data
        procedure :: ConstructUnstructured => ConstructSI2DUS

        ! Evaluator
        procedure :: Evaluate   => EvaluateStructuredInterpolant2D

        ! Derivatives
        procedure :: EvaluateDiffCoef2Val     => EvaluateCoefficientDerivativesUniformStructured
        procedure :: EvaluateDiffInterp2Coef  => EvaluateInterpolantDerivativesUniformStructured
        procedure :: EvaluateDiffInterp2Val   => EvaluateDerivativesVq2VinitUniformStructured

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
        integer(I16), allocatable               :: allprefac(:), precomputedfac(:)
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
                ! A(k, :) = reshape((spread(xstencil(i)**(/ix, ix=0,M/), 1, M+1)).*(repmat( (ystencil(j).^(0:M)), M+1, 1))', 1, []) 
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
        nyr = iylb + ny-iyub-1
        allocate(yrange(nyr), dystencil(nyr))
        xrange = [(i, i = ixlb, ixub)]
        yrange = [ [(i, i = 1, iylb-1)], [(i, i = iyub+1, ny)] ]
        dystencil = [ [(i, i = iylb-1, 1, -1)], -[(i, i = 1, ny-iyub, 1)] ] 
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
                    k = k + 1  
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
        nxr = ixlb + nx-ixub-1
        nyr = iyub - iylb + 1
        allocate(xrange(nxr), yrange(nyr), dxstencil(nxr))
        yrange = [(i, i = iylb, iyub)]
        xrange = [ [(i, i = 1, ixlb-1)], [(i, i = ixub+1, nx)] ]
        dxstencil = [ [(i, i = ixlb-1, 1, -1)], -[(i, i = 1, nx-ixub)] ] 
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
                    k = k + 1  
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
        nxr = ixlb + nx-ixub-1
        nyr = iylb + ny-iyub-1
        allocate(xrange(nxr), yrange(nyr), dxstencil(nxr), dystencil(nyr))
        xrange = [ [(i, i = 1, ixlb-1)], [(i, i = ixub+1, nx)] ]
        dxstencil = [ [(i, i = ixlb-1, 1, -1)], -[(i, i = 1, nx-ixub)] ] 
        yrange = [ [(i, i = 1, iylb-1)], [(i, i = iyub+1, ny)] ]
        dystencil = [ [(i, i = iylb-1, 1, -1)], -[(i, i = 1, ny-iyub)] ] 

        ! Loop
        ll = 1 
        do i = 1, nxr
            ! Stencils
            xstencil = stencil + dxstencil(ll) 
            
            kk = 1 
            do j = 1, nyr
                ! Stencils
                ystencil = stencil + dystencil(kk) 

                ! Compute coefficients
                k = 1 
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
                        k = k + 1 
                    end do
                end do
                
                ! Get residuals
                vID = vertind(xrange(i), yrange(j)) 
                res(vID, :) = reshape(v(xstencil + xrange(i), ystencil + yrange(j)), (/nterms/)) 
                    
                ! Solve for boundary vertices only
                allocate(tvID(1), temp(nterms, 1), ipiv(nterms))
                tvID = vertind(xrange(i), yrange(j))
                temp = transpose(res(tvID, :))
                call dgesv(nterms, 1, A, nterms, ipiv, temp, nterms, info)
                cij(tvID, :) = transpose(temp(:, :))
                deallocate(tvID, temp, ipiv)
                
                ! Update kk
                kk = kk + 1 
                
            end do
            
            ! Update ll 
            ll = ll + 1 
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
        allocate(allprefac((M+1)**2), fvals(size(v, 1), size(v, 2), nterms), precomputedfac(2*(C+1)))
        precomputedfac(1:2) = 1
        do i = 3, 2*(C+1)
            precomputedfac(i) = precomputedfac(i-1)*int(i-1, kind=I16)
        end do 
        interp%precomputedfac = precomputedfac
        k = 1
        do j = 0, M
            do i  = 0, M 
                allprefac(k) = int(precomputedfac(i+1)*precomputedfac(j+1), kind=I16)
                k = k + 1
            end do 
        end do 
        fvals(:, :, :) = 0
        !$omp parallel do default(shared)
        do k = 1, nterms
            fvals(:, :, k) = allprefac(k)*reshape(cij(:, k), (/nx,  ny/))
        end do 
        !$omp end parallel do 

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
                        prefac = real(precomputedfac(l+1), kind=R8)/real(precomputedfac(l+1-j), kind = R8)*&
                            real(precomputedfac(k+1), kind=R8)/real(precomputedfac(k+1-i), kind = R8)
                        
                        ! Index of the term
                        termind = l*(2*(C+1)) + k + 1
        
                        ! Contribution of each point
                        A(1+4*derivind, termind) = prefac*0**(k - i)*0**(l - j)
                        A(2+4*derivind, termind) = prefac*1**(k - i)*0**(l - j)
                        A(3+4*derivind, termind) = prefac*1**(k - i)*1**(l - j)
                        A(4+4*derivind, termind) = prefac*0**(k - i)*1**(l - j)
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
        !$omp parallel do default(shared) private(termind, vecind, vec)
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
                        vec(1+vecind*4) = fvals(i, j, termind+1)
                        vec(2+vecind*4) = fvals(i+1, j, termind+1)
                        vec(3+vecind*4) = fvals(i+1, j+1, termind+1)
                        vec(4+vecind*4) = fvals(i, j+1, termind+1)
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
        !$omp end parallel do

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
        !allocate(interp%xgv(nx), interp%ygv(ny), &
        !    interp%A(size(aij, 1), size(aij, 2)), interp%refx(nx-1), &
        !    interp%refy(ny-1), interp%refdx(nx-1), interp%refdy(ny-1), &
        !    interp%cellindex(nx-1, ny-1))
        
        

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

        ! Not implemented yet
        call gdErrorHandler('Non-uniform structured 2D interpolant constructor not implemented yet')

    end subroutine 

    ! Structured evaluator
    subroutine EvaluateStructuredInterpolant2D(interp, xq, yq, derivx, &
        derivy, vq)

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
        real(R8), intent(in)                    :: xq(:), yq(:)
        real(R8), intent(out)                   :: vq(:)
        integer(I8), intent(in)                 :: derivx, derivy

        ! Auxiliary
        integer(I8)                             :: nq
        real(R8)                                :: prefac

        integer(I8), allocatable                :: ind(:), ind_orig(:) 
        real(R8), allocatable                   :: xqn(:), yqn(:), &
            term(:), thisA(:, :)

        ! Loop
        integer(I8)                             :: i, j, indder

        ! Initialize
        !===========
        ! Check inputs
        if (size(xq, 1) .ne. size(yq, 1)) then 
            call gdErrorHandler('EvaluateStructuredInterpolant2D: ' // &
            'query point coordinates xq and yq have dissimilar ' // &
            'dimensions, check input')
        end if 
        if (size(vq, 1) .ne. size(xq, 1)) then 
            call gdErrorHandler('EvaluateStructuredInterpolant2D: ' // &
            'query point values vq does not have the same dimensions ' // &
            'as query point coordinates xq, yq, check input')
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
            precfac => interp%precomputedfac, &
            cellindex => interp%cellindex)

        ! Get cell index of query point
        call GetIndex(xq, yq, xv, yv, ind)

        ! Set zero indices temporarily to one (set to NaN later)
        allocate(ind_orig, source=ind)
        ind_orig = ind
        where (ind == 0) ind = 1

        ! Extract
        allocate(thisA(size(ind, 1), size(A, 2)))
        thisA = A(ind, :)

        ! Compute normalized query points
        allocate(xqn(nq), yqn(nq))
        xqn = (xq - refx(ind))/refdx(ind)
        yqn = (yq - refy(ind))/refdy(ind)

        ! Evaluate and sum each term
        allocate(term(nq))
        vq(:) = 0

        !$omp parallel default(private) shared(derivx, derivy, &
        !$omp xqn, yqn, thisA, interp, vq)
        !$omp do
        do i = derivx, n 
            do j = derivy, n
                ! Derivative index
                indder = j*(n+1) + i + 1

                ! Factorial prefactor
                prefac = real(precfac(i+1), kind=R8)/&
                    real(precfac(i+1-derivx), kind = R8)*&
                    real(precfac(j+1), kind=R8)/real(precfac(j+1-derivy), kind = R8)
                term = thisA(:, indder)*prefac*xqn**(i - derivx)*yqn**(j - derivy)
                !$omp critical
                vq = vq + term
                !$omp end critical
            end do 
        end do
        !$omp end do 
        !$omp end parallel

        ! Scale
        vq = vq/( (refdx(ind))**derivx * (refdy(ind))**derivy)

        ! Set zero indices to NaN
        where (ind_orig == 0) vq = nanval_R8()

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

    ! Structured coefficient derivatives
    subroutine EvaluateCoefficientDerivativesUniformStructured(&
        interp, xgv, ygv, v, dadv)

        ! Description
        !============
        ! This routine computes the derivatives of the interpolant 
        ! coefficients a with respect to the input values v that are 
        ! given on a structured mesh. To this end, v is considered to be
        ! reshaped to a rank 1 array in column major order. 

        ! Declare variables
        !==================
        ! Arguments
        class(StructuredInterpolant2DUDT)       :: interp 
        real(R8), intent(in)                    :: xgv(:), ygv(:), &
            v(:, :)
        type(MySparseUDT), intent(out)          :: dadv

        ! Auxiliary
        type(MySparseUDT)                       :: dresdv, dcijdres, &
            dfvalsdcij, drhsdfvals, daijdrhs
        integer(I8)                             :: nx, ny, &
            ixlb, ixub, iylb, iyub, nxr, nyr, vID, termind, &
            nc, vecind 
        real(R8)                                :: dxmean, dymean, &
            reldevx, reldevy, prefac

        integer(I8), allocatable                :: stencil(:), &
            vertind(:, :), xstencil(:), ystencil(:), tvID(:), &
            xrange(:), yrange(:), dxstencil(:), dystencil(:), &
            cellindex(:, :), resID(:, :), cijID(:, :), fvalsID(:, :, :), &
            tcID(:), aijID(:, :), rhsID(:, :)
        integer(I16), allocatable               :: allprefac(:)
        real(R8), allocatable                   :: xg(:, :), &
            yg(:, :), cij(:, :), A(:, :), dx(:), dy(:), res(:, :), &
            temp(:, :), fvals(:, :, :), refx(:), refy(:), refdx(:), &
            refdy(:), rhs(:, :), vec(:), aij(:, :), Ainv(:, :), sol(:, :), &
            dcijdres_aux(:, :, :), multfac(:)
            

        ! Loop
        integer(I8)                             :: i, j, k, l, ix, iy, &
            ii, jj, kk, ll, derivind, derivcount, nterms, cc1, cc2, cc

        ! Linear solver variables (dummies basically)
        integer                             :: info

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
            ! Issue $warning, set M = C
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
        allocate(cij(nx*ny, nterms), A(nterms, nterms), res(nx*ny, nterms), &
            dcijdres_aux(nx*ny, nterms, nterms))
        cij(:, :) = 0
        Ainv = A 

        ! Initialize derivative contributions
        resID = reshape([(k, k = 1, nx*ny*nterms)], [nx*ny, nterms]) ! same as cijID 
        cijID = resID 
        call dresdv%Initialize(nx*ny*nterms, nx*ny, nx*ny*nterms)
        call dcijdres%Initialize(nx*ny*nterms, nx*ny*nterms, nx*ny*nterms*nterms)
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
                ! A(k, :) = reshape((spread(xstencil(i)**(/ix, ix=0,M/), 1, M+1)).*(repmat( (ystencil(j).^(0:M)), M+1, 1))', 1, []) 
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

        ! Initialize counters
        cc1 = 0
        cc2 = 0

        ! Loop
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

                ! Compute linearization
                dresdv%row(cc1+1:cc1+nterms) = resID(vID, :)
                dresdv%col(cc1+1:cc1+nterms) = reshape(vertind(xstencil, ystencil), [nterms])
                dresdv%val(cc1+1:cc1+nterms) = 1
                cc1 = cc1 + nterms

            end do 
        end do 

        ! Solve for internal vertices only
        allocate(tvID(nxr*nyr), temp(nterms, nxr*nyr))
        tvID = reshape(vertind(xrange, yrange), (/nxr*nyr/))
        temp = transpose(res(tvID, :))
        sol = temp
        call SolveDenseLinearSystemDI(A, temp, sol, info, Ainv)
        ! call dgesv(nterms, nxr*nyr, A, nterms, ipiv, temp, nterms, info)
        cij(tvID, :) = transpose(sol)
        
        ! Compute linearization
        cc2 = 0
        do k = 1, size(tvID, 1)
            do j = 1, nterms
                dcijdres%row(cc2+1:cc2+nterms) = cijID(tvID(k), j)
                dcijdres%col(cc2+1:cc2+nterms) = resID(tvID(k), :)
                dcijdres%val(cc2+1:cc2+nterms) = Ainv(j, :)
                cc2 = cc2 + nterms 
            end do
        end do 

        ! Housekeeping
        deallocate(tvID, temp)

        ! Top & bottom boundary
        !----------------------
        ! For each iy we need to construct a new interpolation matrix, but
        ! for ix it remains the same
        deallocate(yrange)
        nyr = iylb + ny-iyub-1
        allocate(yrange(nyr), dystencil(nyr))
        xrange = [(i, i = ixlb, ixub)]
        yrange = [ [(i, i = 1, iylb-1)], [(i, i = iyub+1, ny)] ]
        dystencil = [ [(i, i = iylb-1, 1, -1)], -[(i, i = 1, ny-iyub, 1)] ] 
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
                    k = k + 1  
                end do
            end do

            ! Loop over bulk vertices to compute weights
            do i = 1, nxr 
                ! Get residuals
                vID = vertind(xrange(i), yrange(j))
                res(vID, :) = reshape(v(xstencil + xrange(i), ystencil + yrange(j)), (/nterms/))
                
                ! Compute linearization
                dresdv%row(cc1+1:cc1+nterms) = resID(vID, :)
                dresdv%col(cc1+1:cc1+nterms) = reshape(vertind(xstencil + xrange(i), ystencil + yrange(j)), [nterms])
                dresdv%val(cc1+1:cc1+nterms) = 1
                cc1 = cc1 + nterms
            end do 

            ! Solve for boundary vertices
            allocate(tvID(nxr), temp(nterms, nxr))
            tvID = reshape(vertind(xrange, yrange(j)), (/nxr/))
            temp = transpose(res(tvID, :))
            sol = temp
            call SolveDenseLinearSystemDI(A, temp, sol, info, Ainv)
            ! call dgesv(nterms, nxr, A, nterms, ipiv, temp, nterms, info)
            cij(tvID, :) = transpose(sol(:, :))

            ! Compute linearization
            do k = 1, size(tvID, 1)
                do l = 1, nterms
                    dcijdres%row(cc2+1:cc2+nterms) = cijID(tvID(k), l)
                    dcijdres%col(cc2+1:cc2+nterms) = resID(tvID(k), :)
                    dcijdres%val(cc2+1:cc2+nterms) = Ainv(l, :)
                    cc2 = cc2 + nterms 
                end do
            end do 

            ! Housekeeping
            deallocate(tvID, temp)

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
        nxr = ixlb + nx-ixub-1
        nyr = iyub - iylb + 1
        allocate(xrange(nxr), yrange(nyr), dxstencil(nxr))
        yrange = [(i, i = iylb, iyub)]
        xrange = [ [(i, i = 1, ixlb-1)], [(i, i = ixub+1, nx)] ]
        dxstencil = [ [(i, i = ixlb-1, 1, -1)], -[(i, i = 1, nx-ixub)] ] 
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
                    k = k + 1  
                end do
            end do

            ! Loop over bulk vertices to compute weights
            do j = 1, nyr 
                ! Get residuals
                vID = vertind(xrange(i), yrange(j))
                res(vID, :) = reshape(v(xstencil + xrange(i), ystencil + yrange(j)), (/nterms/))
            
                ! Compute linearization
                dresdv%row(cc1+1:cc1+nterms) = resID(vID, :)
                dresdv%col(cc1+1:cc1+nterms) = reshape(vertind(xstencil + xrange(i), ystencil + yrange(j)), [nterms])
                dresdv%val(cc1+1:cc1+nterms) = 1
                cc1 = cc1 + nterms
            end do 

            ! Solve for boundary vertices
            allocate(tvID(nyr), temp(nterms, nyr))
            tvID = reshape(vertind(xrange(i), yrange), (/nyr/))
            temp = transpose(res(tvID, :))
            sol = temp
            call SolveDenseLinearSystemDI(A, temp, sol, info, Ainv)
            ! call dgesv(nterms, nxr*nyr, A, nterms, ipiv, temp, nterms, info)
            cij(tvID, :) = transpose(sol)
            
            ! Compute linearization
            do k = 1, size(tvID, 1)
                do j = 1, nterms
                    dcijdres%row(cc2+1:cc2+nterms) = cijID(tvID(k), j)
                    dcijdres%col(cc2+1:cc2+nterms) = resID(tvID(k), :)
                    dcijdres%val(cc2+1:cc2+nterms) = Ainv(j, :)
                    cc2 = cc2 + nterms 
                end do
            end do 

            ! Housekeeping
            deallocate(tvID, temp)

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
        nxr = ixlb + nx-ixub-1
        nyr = iylb + ny-iyub-1
        allocate(xrange(nxr), yrange(nyr), dxstencil(nxr), dystencil(nyr))
        xrange = [ [(i, i = 1, ixlb-1)], [(i, i = ixub-1, nx)] ]
        dxstencil = [ [(i, i = ixlb-1, 1, -1)], -[(i, i = 1, nx-ixub)] ] 
        yrange = [ [(i, i = 1, iylb-1)], [(i, i = iyub-1, ny)] ]
        dystencil = [ [(i, i = iylb-1, 1, -1)], -[(i, i = 1, ny-iyub)] ] 

        ! Loop
        ll = 1 
        do i = 1, nxr
            ! Stencils
            xstencil = stencil + dxstencil(ll) 
            
            kk = 1 
            do j = 1, nyr
                ! Stencils
                ystencil = stencil + dystencil(kk) 

                ! Compute coefficients
                k = 1 
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
                        k = k + 1 
                    end do
                end do
                
                ! Get residuals
                vID = vertind(xrange(i), yrange(j)) 
                res(vID, :) = reshape(v(xstencil + xrange(i), ystencil + yrange(j)), (/nterms/)) 
                    
                ! Compute linearization
                dresdv%row(cc1+1:cc1+nterms) = resID(vID, :)
                dresdv%col(cc1+1:cc1+nterms) = reshape(vertind(xstencil + xrange(i), ystencil + yrange(j)), [nterms])
                dresdv%val(cc1+1:cc1+nterms) = 1
                cc1 = cc1 + nterms

                ! Solve for boundary vertices only
                allocate(tvID(1), temp(nterms, 1))
                tvID = vertind(xrange(i), yrange(j))
                temp = transpose(res(tvID, :))
                sol = temp
                call SolveDenseLinearSystemDI(A, temp, sol, info, Ainv)
                ! call dgesv(nterms, nxr*nyr, A, nterms, ipiv, temp, nterms, info)
                cij(tvID, :) = transpose(sol)

                ! Compute linearization
                do k = 1, size(tvID, 1)
                    do l = 1, nterms
                        dcijdres%row(cc2+1:cc2+nterms) = cijID(tvID(k), l)
                        dcijdres%col(cc2+1:cc2+nterms) = resID(tvID(k), :)
                        dcijdres%val(cc2+1:cc2+nterms) = Ainv(l, :)
                        cc2 = cc2 + nterms 
                    end do
                end do 

                ! Housekeeping
                deallocate(tvID, temp)
                
                ! Update kk
                kk = kk + 1 
                
            end do
            
            ! Update ll 
            ll = ll + 1 
        end do

        ! Compute actual derivatives
        !---------------------------
        ! Correct for scaling - note: for the derivatives, we account 
        ! for this later on in dfvalsdcij
        allocate(multfac(nterms))
        k = 1
        do j = 0, M
            do i = 0, M 
                cij(:, k) = cij(:, k)/( (dxmean**i) * (dymean**j))
                multfac(k) = 1.0/( (dxmean**i) * (dymean**j))
                k = k + 1
            end do 
        end do

        ! Construct derivatives as fields - need to account for 
        ! factorial prefactor!
        allocate(allprefac((M+1)**2), fvals(nx, ny, nterms), &
            fvalsID(nx, ny, nterms))

        ! Initialize derivatives
        call dfvalsdcij%Initialize(nx*ny*nterms, nx*ny*nterms, nx*ny*nterms)

        ! Compute
        k = 1
        do j = 0, M
            do i  = 0, M 
                allprefac(k) = int(MyFactorial(i)*MyFactorial(j), kind=I16)
                k = k + 1
            end do 
        end do 
        fvals(:, :, :) = 0
        fvalsID = reshape([(k, k = 1, nx*ny*nterms)], [nx, ny, nterms])
        cc = 0
        do k = 1, nterms
            fvals(:, :, k) = allprefac(k)*reshape(cij(:, k), (/nx,  ny/))

            ! Compute derivative, account for multfac
            dfvalsdcij%row(cc+1:cc+nx*ny) = reshape(fvalsID(:, :, k), [nx*ny])
            dfvalsdcij%col(cc+1:cc+nx*ny) = cijID(:, k)
            dfvalsdcij%val(cc+1:cc+nx*ny) = real(allprefac(k), kind=R8)*multfac(k)
            cc = cc + nx*ny
                
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
                        A(1+4*derivind, termind) = prefac*0**(k - i)*0**(l - j) 
                        A(2+4*derivind, termind) = prefac*1**(k - i)*0**(l - j) 
                        A(3+4*derivind, termind) = prefac*1**(k - i)*1**(l - j) 
                        A(4+4*derivind, termind) = prefac*0**(k - i)*1**(l - j) 
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
        cellindex = reshape([(i, i = 1, nc)], (/nx-1, ny-1/))
        aijID = reshape([(k, k = 1, nc*nterms)], [nc, nterms])
        rhsID = aijID

        ! Initialize derivatives
        call drhsdfvals%Initialize(nc*nterms, size(fvals), nc*nterms)
        call daijdrhs%Initialize(nc*nterms, nc*nterms, nc*nterms*nterms)

        ! Loop over all quads
        cc = 0
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
                        vec(1+vecind*4) = fvals(i, j, termind+1)
                        vec(2+vecind*4) = fvals(i+1, j, termind+1)
                        vec(3+vecind*4) = fvals(i+1, j+1, termind+1)
                        vec(4+vecind*4) = fvals(i, j+1, termind+1)

                        ! Compute linearization, account for later multiplication with dxmean**k, dymean**l
                        drhsdfvals%col(cc+1:cc+4) = [fvalsID(i, j, termind+1), &
                            fvalsID(i+1, j, termind+1), fvalsID(i+1, j+1, termind+1), & 
                            fvalsID(i, j+1, termind+1)]
                        drhsdfvals%row(cc+1:cc+4) = rhsID(cellindex(i, j), 1+vecind*4:4+vecind*4)
                        drhsdfvals%val(cc+1:cc+4) = 1.0*(dxmean**k)*(dymean**l) ! account for later multiplication
                        cc = cc + 4
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
        allocate(temp(size(rhs, 2), size(rhs, 1)))
        temp = transpose(rhs) 
        Ainv = A 
        sol = temp 
        call SolveDenseLinearSystemDI(A, temp, sol, info, Ainv)
        rhs = transpose(sol)
        tcID = reshape(cellindex, [nc])
        aij(tcID, :) = rhs(tcID, :)

        ! Compute linearization
        cc2 = 0
        do k = 1, size(tcID, 1)
            do j = 1, nterms
                daijdrhs%row(cc2+1:cc2+nterms) = aijID(tcID(k), j)
                daijdrhs%col(cc2+1:cc2+nterms) = rhsID(tcID(k), :)
                daijdrhs%val(cc2+1:cc2+nterms) = Ainv(j, :)
                cc2 = cc2 + nterms 
            end do
        end do  

        ! Compute full linearization
        !===========================
        ! Apply chain rule
        dadv = daijdrhs*drhsdfvals
        dadv = dadv*dfvalsdcij 
        dadv = dadv*dcijdres
        dadv = dadv*dresdv 
        !dadv = daijdrhs*drhsdfvals*dfvalsdcij*dcijdres*dresdv

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Structured interpolant derivatives w.r.t. coefficients
    subroutine EvaluateInterpolantDerivativesUniformStructured(&
        interp, xq, yq, derivx, derivy, vq, dvqda)

        ! Description
        !============
        ! This routine evaluates the derivative of the values of the 
        ! interpolant w.r.t. the coefficients of the interpolant. 
        ! Basically, these are the base functions evaluated at the
        ! query coordinates. 

        ! Declare variables
        !==================
        ! Arguments
        class(StructuredInterpolant2DUDT)       :: interp 
        real(R8), intent(in)                    :: xq(:), yq(:)
        real(R8), intent(out)                   :: vq(:)
        integer(I8), intent(in)                 :: derivx, derivy
        type(MySparseUDT)                       :: dvqda

        ! Auxiliary
        integer(I8)                             :: nq
        real(R8)                                :: prefac

        integer(I8), allocatable                :: ind(:), ind_orig(:), &
            aID(:, :)
        real(R8), allocatable                   :: xqn(:), yqn(:), &
            term(:), thisA(:, :), tempval(:)

        ! Loop
        integer(I8)                             :: i, j, k, indder, cc

        ! Initialize
        !===========
        ! Check inputs
        if (size(xq, 1) .ne. size(yq, 1)) then 
            call gdErrorHandler('EvaluateStructuredInterpolant2D: ' // &
            'query point coordinates xq and yq have dissimilar ' // &
            'dimensions, check input')
        end if 
        if (size(vq, 1) .ne. size(xq, 1)) then 
            call gdErrorHandler('EvaluateStructuredInterpolant2D: ' // &
            'query point values vq does not have the same dimensions ' // &
            'as query point coordinates xq, yq, check input')
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

        ! Initialize
        aID = reshape([(k, k = 1, size(A, 1)*size(A, 2))], &
            [size(A, 1), size(A, 2)])
        call dvqda%Initialize(nq, size(aID), size(A, 2)*nq)

        ! Get cell index of query point
        call GetIndex(xq, yq, xv, yv, ind)

        ! Set zero indices temporarily to one (set to NaN later)
        allocate(ind_orig, source=ind)
        ind_orig = ind
        where (ind == 0) ind = 1

        ! Extract
        allocate(thisA(size(ind, 1), size(A, 2)))
        thisA = A(ind, :)

        ! Compute normalized query points
        allocate(xqn(nq), yqn(nq))
        xqn = (xq - refx(ind))/refdx(ind)
        yqn = (yq - refy(ind))/refdy(ind)

        ! Evaluate and sum each term
        allocate(term(nq))
        vq(:) = 0
        cc = 0
        do i = derivx, n 
            do j = derivy, n
                ! Derivative index
                indder = j*(n+1) + i + 1

                ! Factorial prefactor
                prefac = real(MyFactorial(i), kind=R8)/&
                    real(MyFactorial(i-derivx), kind = R8)*&
                    real(MyFactorial(j), kind=R8)/real(MyFactorial(j-derivy), kind = R8)
                term = thisA(:, indder)*prefac*xqn**(i - derivx)*yqn**(j - derivy)
                
                ! Compute value
                vq = vq + term
                
                ! Compute derivatives
                dvqda%col(cc+1:cc+nq) = aID(ind, indder)
                dvqda%row(cc+1:cc+nq) = [(k, k = 1, nq)]
                tempval = prefac*xqn**(i - derivx)*yqn**(j - derivy) & 
                    /( (refdx(ind))**derivx * (refdy(ind))**derivy)
                where (ind_orig == 0) tempval = nanval_R8() ! Set NaNs at out of bounds points
                dvqda%val(cc+1:cc+nq) = tempval
                cc = cc + nq
            end do 
        end do

        ! Scale
        vq = vq/( (refdx(ind))**derivx * (refdy(ind))**derivy)

        ! Set zero indices to NaN
        where (ind_orig == 0) vq = nanval_R8()

        ! Housekeeping
        !=============
        deallocate(xqn, yqn, term, ind)
        end associate

    end subroutine

    ! Derivatives of query values w.r.t. initial values
    subroutine EvaluateDerivativesVq2VinitUniformStructured(&
        interp, xq, yq, derivx, derivy, dvqdv)

        ! Description
        !============
        ! Routine that directly evaluates the derivative (jacobian) of 
        ! the values queried on points xq, yq (vq) w.r.t. the initial
        ! values used to set up the interpolant v (in structured format)
        ! In theory this can be computed by multiplying the derivatives
        ! dvqda*dadv obtained by other subroutines, but dadv is typically
        ! (very) memory-intensive. Here, we exploit the fact that
        ! we can compute the derivative contributions for each point 
        ! separately, since it only depends on the interpolant cell in 
        ! which the point lies. The number of entries in dvqdv should 
        ! therefore be nq*na. Since na is small (< 100) and constant,
        ! the number of non-zeros only scales with nq.  Since:
        !
        !       vq = sum_i sum_j a_ij(b(v)) xq^i yq^j 
        !       a_k = (A^-1 b(v))_(k, :)
        !       b(v) = B^-1 v 
        !
        ! we have that 
        !
        !       dvqdv = sum_i sum_j (A^-1)_(k, :) B^-1 xq^i yq^j  
        !
        ! Since A is the same for each cell, A^-1 can be inverted once. 
        ! B is the same in the bulk of the domain for each cell, 
        ! but varies for most vertices near the edges of the domain. 
        ! For these edges, B^-1 is recomputed from scratch. Given that 
        ! most query points will likely lie within the domain, this 
        ! additional cost should be small. 

        ! The algorithm proceeds as follows:
        ! 0) The necessary 'state' variables are reconstructed (
        ! sampling grid on which v was evaluated, A^-1, B^-1 for bulk
        ! vertices, ...)
        ! 1) For each query point, A^-1*B^-1 is computed and added to the
        ! linearization matrix 

        ! Declare variables
        !==================
        ! Arguments
        class(StructuredInterpolant2DUDT)       :: interp 
        real(R8), intent(in)                    :: xq(:), yq(:)
        integer(I8), intent(in)                 :: derivx, derivy
        type(MySparseUDT), intent(out)          :: dvqdv

        ! Auxiliary
        integer(I8)                             :: nq, nxv, nyv, na, &
            ixlb, ixub, iylb, iyub, nterms, flag, nvals, info, termind, &
            vecind
        integer(I8), allocatable, dimension(:)  :: row, col, stencil, &
            ix, iy, ic, ix_orig, iy_orig, ic_orig, xstencil, ystencil, &
            valindx, valindy, xstencil0, ystencil0
        integer(I8), allocatable, dimension(:, :)   :: tempcol, &
            tempvertind, vertind, tvID
            
        real(R8)                                :: xqn, yqn, dxmean, &
            dymean, tempfac
        real(R8), allocatable, dimension(:)     :: val, &
            eiga, prefac, allfactorials, tempb, sol
        real(R8), allocatable, dimension(:, :)  :: temp, tempval, A, &
            Ainv, B, Binv, Binvbulk

        logical, allocatable                    :: inbulk(:, :)

        ! Loop
        integer(I8)                             :: i, j, k, l, indder, &
            derivind, ii, jj, cc, iq


        ! Initialize
        !===========
        ! Associate
        associate(&
            refx        => interp%refx,     &
            refy        => interp%refy,     &
            refdx       => interp%refdx,    &
            refdy       => interp%refdy,    &
            M           => interp%M,        &
            C           => interp%C,        &
            xgv         => interp%xgv,      &
            ygv         => interp%ygv,      &
            n           => interp%n         &
            )

        ! Compute dimensions
        nq = size(xq)
        nxv = size(interp%xgv)
        nyv = size(interp%ygv)
        na = size(interp%a, 2)
        nvals = (M+2)**2

        ! Initialize linearization
        allocate(row(nq*nvals), col(nq*nvals), val(nq*nvals))
        
        ! Compute preliminaries
        !======================
        ! Determine stencil for derivative approximation in vertices
        allocate(stencil(M+1))
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

        ! Set boundary indices
        ixlb = abs(stencil(1)) + 1
        ixub = nxv - abs(stencil(M+1)) 
        iylb = abs(stencil(1)) + 1
        iyub = nyv - abs(stencil(M+1)) 

        ! Determine some fixed values
        dxmean = sum(xgv(2:nxv) - xgv(1:nxv-1))/(nxv-1)
        dymean = sum(ygv(2:nyv) - ygv(1:nyv-1))/(nyv-1)

        ! Initialize temporary variables
        allocate(tempval(M+2, M+2), tempcol(M+2, M+2), tempvertind(M+2, M+2), &
            temp(M+1, M+1))

        ! Initialize basis function arrays
        allocate(eiga((n+1)**2), allfactorials((n+1)**2))

        ! Determine which points lie in the bulk of the domain
        allocate(inbulk(size(xgv), size(ygv)))
        inbulk = .false. 
        inbulk(ixlb:ixub, iylb:iyub) = .true. 

        ! Determine vertex index
        vertind = reshape([(i, i = 1, nxv*nyv)], [nxv, nyv])

        ! Determine factorial prefactors
        allfactorials = 0
        do i = derivx, n 
            do j = derivy, n 
                ! Compute
                indder = j*(n+1) + i + 1
                allfactorials(indder) = real(MyFactorial(i), kind=R8)/&
                    real(MyFactorial(i-derivx), kind = R8)*&
                    real(MyFactorial(j), kind=R8)/real(MyFactorial(j-derivy), kind = R8)
            end do
        end do
        
        ! Determine A^-1
        !---------------
        nterms = (2*(C+1))**2
        allocate(A(nterms, nterms), tempb(nterms))
        A       = 0 
        Ainv    = A 
        tempb   = 0 ! doesn't matter
        sol     = tempb

        ! Loop over all derivatives
        derivind = 0
        do j = 0, C 
            do i = 0, C
                ! Construct coefficients of the derivative of this 
                ! polynomial
                do l = j, 2*(C+1)-1
                    do k = i, 2*(C+1)-1
                        ! Factorial prefactor
                        tempfac = real(MyFactorial(l), kind=R8)/real(MyFactorial(l-j), kind = R8)*&
                            real(MyFactorial(k), kind=R8)/real(MyFactorial(k-i), kind = R8)
                        
                        ! Index of the term
                        termind = l*(2*(C+1)) + k + 1
        
                        ! Contribution of each point
                        A(1+4*derivind, termind) = tempfac*0**(k - i)*0**(l - j) 
                        A(2+4*derivind, termind) = tempfac*1**(k - i)*0**(l - j) 
                        A(3+4*derivind, termind) = tempfac*1**(k - i)*1**(l - j) 
                        A(4+4*derivind, termind) = tempfac*0**(k - i)*1**(l - j) 
                    end do
                end do

                ! Update counter
                derivind = derivind + 1
            end do
        end do

        ! Compute A^-1 
        call SolveDenseLinearSystemDI(A, tempb, sol, flag, Ainv)
        deallocate(A, tempb, sol)

        ! Determine B^-1 in the bulk of the domain
        !-----------------------------------------
        ! Initialize
        nterms = (M + 1)**2
        allocate(A(nterms, nterms), tempb(nterms), sol(nterms))
        tempb   = 0
        sol     = 0
        Binvbulk = A
        
        ! Loop to compute 'weights' of coefficients
        xstencil = stencil 
        ystencil = stencil 
        k = 1
        do j = 1, M+1
            do i = 1, M+1
                ! Compute coefficient matrix entry
                ! A(k, :) = reshape((spread(xstencil(i)**(/ix, ix=0,M/), 1, M+1)).*(repmat( (ystencil(j).^(0:M)), M+1, 1))', 1, []) 
                l = 1
                do jj = 0, M
                    do ii = 0, M
                        A(k, l) = xstencil(i)**ii * ystencil(j)**jj
                        l = l + 1
                    end do 
                end do

                ! Update
                k = k + 1
            end do
        end do

        ! Invert
        call SolveDenseLinearSystemDI(A, tempb, sol, info, Binvbulk)
        
        ! Compute prefactor
        allocate(prefac((M+1)**2))
        k = 1
        do j = 0, M
            do i = 0, M 
                ! Compute
                prefac(k) = real(int(MyFactorial(i)*MyFactorial(j), kind=I16), kind=R8) & 
                    /(dxmean**i * dymean**j)

                ! Update counter
                k = k + 1
            end do
        end do

        ! Determine cell indices in matrix and linear format
        call GetIndex(xq, yq, xgv, ygv, ix, iy)
        call GetIndex(xq, yq, xgv, ygv, ic) 

        ! Hedge for out of bounds
        ic_orig = ic 
        ix_orig = ix 
        iy_orig = iy 
        where (ic == 0) ic = 1
        where (ix == 0) ix = 1
        where (iy == 0) iy = 1
        

        ! Compute linearization
        !======================
        ! Loop over all query points
        cc = 0
        nterms = (2*(C+1))**2
        allocate(tvID(M+2, M+2))
        do iq = 1, nq

            ! Initialize 
            !-----------
            ! Values
            tempval = 0
            tvID = 0

            ! Normalized coordinates
            xqn = (xq(iq) - refx(ic(iq)))/refdx(ic(iq))
            yqn = (yq(iq) - refy(ic(iq)))/refdy(ic(iq))

            ! Shape function vector 
            eiga = 0
            k = 1
            do i = derivx, n 
                do j = derivy, n
                    indder = j*(n+1) + i + 1
                    eiga(indder) = allfactorials(indder)* &
                        (xqn)**(i - derivx)*yqn**(j - derivy)
                end do
            end do

            ! Contribution of vertex (ix, iy)
            !--------------------------------
            ! Determine B^-1
            if (inbulk(ix(iq), iy(iq))) then 
                ! Just reuse bulk inverse
                Binv = Binvbulk 
                xstencil = stencil + ix(iq)
                ystencil = stencil + iy(iq)
            else
                ! Compute on the fly
                call ComputeCoefficientMatrix(stencil, ix(iq), iy(iq), B, Binv, &
                    M, xstencil, ystencil, nxv, nyv)
                xstencil = xstencil + ix(iq)
                ystencil = ystencil + iy(iq)
            end if 

            ! Save original stencil for comparison
            xstencil0 = xstencil
            ystencil0 = ystencil

            ! Set vertex ID
            tvID(1:M+1, 1:M+1) = vertind(xstencil, ystencil)

            ! Loop 
            temp = 0
            do ii = 1, na
                do l = 0, C 
                    do k = 0, C 
                        termind = l*(M + 1) + k
                        vecind = l*(C + 1) + k
                        temp = temp + reshape(Ainv(ii, 1+vecind*4)*Binv(termind+1, :) &
                            *prefac(termind+1)*eiga(ii)/ &
                            (refdx(ic(iq))**derivx * refdy(ic(iq))**derivy), &
                            [M+1, M+1])*(dxmean**k)*(dymean**l)
                    end do 
                end do
            end do
            !do k = 1, na 
            !    do l = 1, (M + 1)**2
            !        temp = temp + reshape(Ainv(k, l)*Binv(l, :) &
            !            *prefac*eiga(k)/ &
            !            (refdx(ic(iq))**derivx * refdy(ic(iq))**derivy), &
            !            [M+1, M+1]) 
            !    end do
            !end do

            ! Add
            tempval(1:M+1, 1:M+1) = tempval(1:M+1, 1:M+1) + temp

            ! Contribution of vertex (ix+1, iy)
            !--------------------------------
            ! Set value indices
            valindx = [(k, k = 2, M+2)]
            valindy = [(k, k = 1, M+1)]

            ! Determine B^-1
            if (inbulk(ix(iq)+1, iy(iq))) then 
                ! Just reuse bulk inverse
                Binv = Binvbulk 
                xstencil = stencil + ix(iq) + 1
                ystencil = stencil + iy(iq)
            
                
            else
                ! Compute on the fly
                call ComputeCoefficientMatrix(stencil, ix(iq)+1, iy(iq), B, Binv, &
                    M, xstencil, ystencil, nxv, nyv)
                xstencil = xstencil + ix(iq) + 1
                ystencil = ystencil + iy(iq)

            end if 

            ! Check stencil
            if (xstencil(1) == xstencil0(1)) then 
                valindx = valindx - 1
            end if 

            ! Check
            if (.not. all(tvID(valindx(1:M), valindy) == vertind(xstencil(1:M), ystencil))) then 
                ! Stop, this shouldn't happen
                call gdErrorHandler('Implementation bug')
            else
                tvID(valindx, valindy) = vertind(xstencil, ystencil)
            end if

            ! Loop
             
            temp = 0
            do ii = 1, na
                do l = 0, C 
                    do k = 0, C 
                        termind = l*(M + 1) + k
                        vecind = l*(C + 1) + k
                        temp = temp + reshape(Ainv(ii, 2+vecind*4)*Binv(termind+1, :) &
                            *prefac(termind+1)*eiga(ii)/ &
                            (refdx(ic(iq))**derivx * refdy(ic(iq))**derivy), &
                            [M+1, M+1])*(dxmean**k)*(dymean**l)
                    end do 
                end do
            end do
            !do k = 1, na 
            !    do l = 1, (M + 1)**2
            !        temp = temp + reshape(Ainv(k, l)*Binv(l, :) &
            !            *prefac*eiga(k)/ &
            !            (refdx(ic(iq))**derivx * refdy(ic(iq))**derivy), &
            !            [M+1, M+1]) 
            !    end do
            !end do

            ! Add
            tempval(valindx, valindy) = tempval(valindx, valindy) + temp

            ! Contribution of vertex (ix+1, iy+1)
            !--------------------------------
            ! Set value indices
            valindx = [(k, k = 2, M+2)]
            valindy = [(k, k = 2, M+2)]

            ! Determine B^-1
            if (inbulk(ix(iq)+1, iy(iq)+1)) then 
                ! Just reuse bulk inverse
                Binv = Binvbulk 
                xstencil = stencil + ix(iq) + 1
                ystencil = stencil + iy(iq) + 1
            else
                ! Compute on the fly
                call ComputeCoefficientMatrix(stencil, ix(iq)+1, iy(iq)+1, B, Binv, &
                    M, xstencil, ystencil, nxv, nyv)
                xstencil = xstencil + ix(iq) + 1
                ystencil = ystencil + iy(iq) + 1
                
            end if 

            ! Check stencil
            if (xstencil(1) == xstencil0(1)) then 
                valindx = valindx - 1
            end if 
            if (ystencil(1) == ystencil0(1)) then 
                valindy = valindy - 1
            end if 

            ! Check
            if (.not. all(tvID(valindx(1:M), valindy(1:M)) == vertind(xstencil(1:M), ystencil(1:M)))) then 
                ! Stop, this shouldn't happen
                call gdErrorHandler('Implementation bug')
            else
                tvID(valindx, valindy) = vertind(xstencil, ystencil)
            end if

            ! Loop 
            temp = 0
            do ii = 1, na
                do l = 0, C 
                    do k = 0, C 
                        termind = l*(M + 1) + k
                        vecind = l*(C + 1) + k
                        temp = temp + reshape(Ainv(ii, 3+vecind*4)*Binv(termind+1, :) &
                            *prefac(termind+1)*eiga(ii)/ &
                            (refdx(ic(iq))**derivx * refdy(ic(iq))**derivy), &
                            [M+1, M+1])*(dxmean**k)*(dymean**l)
                    end do 
                end do
            end do
            !do k = 1, na 
            !    do l = 1, (M + 1)**2
            !        temp = temp + reshape(Ainv(k, l)*Binv(l, :) &
            !            *prefac*eiga(k)/ &
            !            (refdx(ic(iq))**derivx * refdy(ic(iq))**derivy), &
            !            [M+1, M+1]) 
            !    end do
            !end do


            ! Add
            tempval(valindx, valindy) = tempval(valindx, valindy) + temp

            ! Contribution of vertex (ix, iy+1)
            !--------------------------------
            ! Set value indices
            valindx = [(k, k = 1, M+1)]
            valindy = [(k, k = 2, M+2)]

            ! Determine B^-1
            if (inbulk(ix(iq), iy(iq)+1)) then 
                ! Just reuse bulk inverse
                Binv = Binvbulk 
                xstencil = stencil + ix(iq) 
                ystencil = stencil + iy(iq) + 1
            else
                ! Compute on the fly
                call ComputeCoefficientMatrix(stencil, ix(iq), iy(iq)+1, B, Binv, &
                    M, xstencil, ystencil, nxv, nyv)
                xstencil = xstencil + ix(iq)
                ystencil = ystencil + iy(iq) + 1

            end if 

            ! Check stencil
            if (ystencil(1) == ystencil0(1)) then 
                valindy = valindy - 1
            end if 

            ! Check
            if (.not. all(tvID(valindx, valindy(1:M)) == vertind(xstencil, ystencil(1:M)))) then 
                ! Stop, this shouldn't happen
                call gdErrorHandler('Implementation bug')
            else
                tvID(valindx, valindy) = vertind(xstencil, ystencil)
            end if

            ! Loop 
            temp = 0
            do ii = 1, na
                do l = 0, C 
                    do k = 0, C 
                        termind = l*(M + 1) + k
                        vecind = l*(C + 1) + k
                        temp = temp + reshape(Ainv(ii, 4+vecind*4)*Binv(termind+1, :) &
                            *prefac(termind+1)*eiga(ii)/ &
                            (refdx(ic(iq))**derivx * refdy(ic(iq))**derivy), &
                            [M+1, M+1])*(dxmean**k)*(dymean**l)
                    end do 
                end do
            end do 
            !do k = 1, na 
            !    do l = 1, (M + 1)**2
            !        temp = temp + reshape(Ainv(k, l)*Binv(l, :) &
            !            *prefac*eiga(k)/ &
            !            (refdx(ic(iq))**derivx * refdy(ic(iq))**derivy), &
            !            [M+1, M+1]) 
            !    end do
            !end do


            ! Add
            tempval(valindx, valindy) = tempval(valindx, valindy) + temp

            ! Add 
            !----
            ! Add row, col, val 
            where (tvID == 0) 
                tvID = 1
                tempval = 0
            end where
            row(cc+1:cc+nvals) = iq
            col(cc+1:cc+nvals) = reshape(tvID, [nvals])
            val(cc+1:cc+nvals) = reshape(tempval, [nvals])

            ! Update counter
            cc = cc + nvals

        end do
        
        ! Check
        if (cc /= nq*nvals) then 
            call gdErrorHandler('Implementation bug')
        end if

        ! Construct sparse matrix
        dvqdv = ConstructMySparse(row, col, val, nq, nxv*nyv)

        ! Housekeeping
        !=============
        end associate


    end subroutine

    subroutine ComputeCoefficientMatrix(stencil, ix, iy, B, Binv, M, &
        xstencil, ystencil, nxv, nyv)

        ! Description
        !============
        ! Auxiliary routine to compute coefficient matrix given a certain
        ! stencil (and its inverse)

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)                 :: stencil(:), ix, iy, M, &
            nxv, nyv
        real(R8), intent(out), allocatable      :: B(:, :), Binv(:, :)
           

        ! Auxiliary
        integer(I8)                             :: dx, dy, info, &
            nterms
        integer(I8), intent(out), allocatable   :: xstencil(:), ystencil(:)

        real(R8), allocatable                   :: tempb(:), sol(:)
            

        ! Loop
        integer(I8)                             :: i, j, ii, jj, k, l

        ! Initialize
        !===========
        ! Check allocation status
        if (allocated(B)) then 
            deallocate(B)
        end if
        if (allocated(Binv)) then 
            deallocate(Binv)
        end if 

        ! Allocate & initialize
        nterms = (M + 1)**2
        allocate(B(nterms, nterms), Binv(nterms, nterms), tempb(nterms), &
            sol(nterms))
        B       = 0
        Binv    = 0
        sol     = 0
        tempb   = 0

        ! Determine stencil
        !==================
        xstencil = stencil 
        ystencil = stencil
        if (any( (stencil + ix) <= 0)) then 
            dx = minval(stencil + ix)-1
            xstencil = stencil - dx 
        elseif (any( (stencil + ix) > nxv)) then 
            dx = nxv - maxval(stencil + ix)
            xstencil = stencil + dx 
        end if
        if (any( (stencil + iy) <= 0)) then 
            dy = minval(stencil + iy)-1
            ystencil = stencil - dy 
        elseif (any( (stencil + iy) > nyv)) then 
            dy = nyv - maxval(stencil + iy)
            ystencil = stencil + dy 
        end if

        ! Compute
        !========
        ! Loop to compute 'weights' of coefficients
        k = 1
        do j = 1, M+1
            do i = 1, M+1
                ! Compute coefficient matrix entry
                ! A(k, :) = reshape((spread(xstencil(i)**(/ix, ix=0,M/), 1, M+1)).*(repmat( (ystencil(j).^(0:M)), M+1, 1))', 1, []) 
                l = 1
                do jj = 0, M
                    do ii = 0, M
                        B(k, l) = xstencil(i)**ii * ystencil(j)**jj
                        l = l + 1
                    end do 
                end do

                ! Update
                k = k + 1
            end do
        end do

        ! Invert
        call SolveDenseLinearSystemDI(B, tempb, sol, info, Binv)
        

    end subroutine
    

end module