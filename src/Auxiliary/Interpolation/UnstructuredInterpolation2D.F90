!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!------------


module UnstructuredInterpolant2D

    ! Initialize
    !===========
    ! Load modules
    use mod_triangulation  
    use Interpolant2D
    use Interpolant2D_auxiliaries


    implicit none
    save

    ! Unstructured interpolant type
    !==============================
    type, extends(GenericInterpolant2DUDT) :: UnstructuredInterpolant2DUDT

        ! Description
        !============
        ! Apart from the fields of the generic interpolant, several 
        ! other fields are defined as well:
        ! - meth: methodology to construct interpolant
        ! - C, M: order of the interpolant and order of the approximation
        ! method to compute derivatives for the interpolant construction     
        ! - n: number of terms in the interpolant    
        ! - allowextrapolation: logical to check if we can extrapolate
        ! - triangulation: triangulated mesh to use for interpolation
        ! - base_func: type of used base function, for now only 'polynomial'

        ! - v: the values at the vertex points
        ! - xgv, ygv: grid vectors
        ! - A: interpolation coefficients
        ! - refx, refy, refdx, refdy: reference values used to compute
        ! derivatives etc 
        ! - cellindex: nx-1 by ny-1 array containing the cell indices

        ! - precomputedfac: the required factorials precomputed to save
        ! some time during evaluation 

        character(:), allocatable       :: meth 
        integer(I8)                     :: C, M, n
        logical                         :: allowextrapolation
        type(TriangulationUDT)          :: triangulation
        character(:), allocatable       :: base_func

        type(GradientReconstructionTriaUDT) :: GR

        real(R8), allocatable           :: xgv(:), ygv(:), A(:, :), &
            refx(:), refy(:), refdx(:), refdy(:)
        integer(I8), allocatable        :: cellindex(:, :)
        integer(I16), allocatable       :: precomputedfac(:)  
         
    
    contains

        ! Parameter setter routine
        procedure :: SetParametersUS 

        ! Construct based on structured data
        procedure :: ConstructStructured => ConstructUSI2DS
    
        ! Construct based on unstructured data
        procedure :: ConstructUnstructured => ConstructUSI2DUS

        ! Evaluator
        procedure :: Evaluate   => EvaluateUnstructuredInterpolant2D  
        procedure :: EvaluateWrapper         

    end type

    contains

    ! Set parameters
    subroutine SetParametersUS(interp, meth, C, M, triangulation, base_func)

        ! Description
        !============
        ! Set the parameters of the interpolation routine
        ! - C:    desired continuity of the interpolant 
        ! - M:    order of the interpolant describing the values 
        !               at the grid nodes

        ! Declare variables
        !==================
        ! Arguments
        class(UnstructuredInterpolant2DUDT)     :: interp
        character(:), allocatable, intent(in)   :: meth
        integer(I8), intent(in)                 :: C, M 
        character(*)                            :: base_func
        type(TriangulationUDT), intent(in)      :: triangulation

        ! Set
        !====
        interp%meth = meth 
        interp%C    = C 
        interp%M    = M
        interp%allowextrapolation = .true.
        interp%triangulation = triangulation
        interp%base_func = base_func
        
    end subroutine

    ! Constructor, structured
    subroutine ConstructUSI2DS(interp, xg, yg, v) 

        ! Declare variables
        !==================
        ! Arguments
        class(UnstructuredInterpolant2DUDT)       :: interp
        real(R8), allocatable                     :: xg(:), yg(:), v(:, :)

        ! Currently, no implementation yet
        call gdErrorHandler('Structured initialization of 2D unstructured interpolant not implemented')

    end subroutine    

    ! Constructor, unstructured
    subroutine ConstructUSI2DUS(interp, xg, yg, v)

        ! This routine is a wrapper for the variants of the unstructured
        ! way of constructing the unstructured interpolant.
        ! - xg and yg: coordinates of point where the values are known
        ! - v: field to interpolate

        ! Declare variables
        !==================
        ! Arguments
        class(UnstructuredInterpolant2DUDT)     :: interp
        real(R8), allocatable                   :: xg(:), yg(:), v(:)

        ! Initialize
        !===========
        ! Check which method to follow
        select case(interp%meth)

        case ('barycentric')

            call ConstructUSI2DUSBarycentric(interp, xg, yg, v)

        case ('finite_element')

            call ConstructUSI2DUSFinEelem(interp, xg, yg, v)

        case default

            call gdErrorHandler('ConstructUSI2DUS: unknown construction method for 2D unstructured interpolant')

        end select

    end subroutine

    subroutine ConstructUSI2DUSBarycentric(interp, xg, yg, v)

        ! Description
        !============
        ! We build the interpolant. For the Barycentric it is
        ! just saving the field information in the vertices.

        ! Declare variables
        !==================
        ! Arguments
        class(UnstructuredInterpolant2DUDT)     :: interp 
        real(R8), intent(in)                    :: xg(:), yg(:)
        real(R8), intent(out)                   :: v(:)

        ! Checks
        if (interp%C /= 0) &
            call gdErrorHandler('ConstructUSI2DUSBarycentric: not higher than C0 continuity with ' // &
                'Barycentric interpolation possible')
        if (interp%M .gt. 1) &
            call gdErrorHandler('ConstructUSI2DUSBarycentric: not higher than linear ' // &
                'interpolation possible with Barycentric interpolation')
        if (size(xg, 1) /= size(yg, 1)) &
            call gdErrorHandler('ConstructUSI2DUSBarycentric: size of xg and yg incompatible')
        if (size(xg) /= size(v)) &
            call gdErrorHandler('ConstructUSI2DUSBarycentric: size of xg and v incompatible')
        if (interp%triangulation%nv /= size(xg)) &
            call gdErrorHandler('ConstructUSI2DUSBarycentric: size of xg and triangulation incompatible')

        ! (xg, yg) are the coordinates where the values are known
        ! This should be the same as the vertex coordinates of the triangulation
        ! The values at the vertices are also saved in the interpolant object
        interp%v = v

    
    end subroutine

    subroutine ConstructUSI2DUSFinEelem(interp, xg, yg, v)

        ! Description
        !============
        ! We build the interpolant. The following steps are taken:
        ! 0) Compute the required derivatives on the vertex nodes using gradient reconstruction
        ! 1) Determine required order of the interpolant depended on the wanted continuity
        !   The interpolant is of type phi(x,y) = sum_i^N sum_j^N a_ij x^i y^j
        ! 2) Compute the coefficients of the gradient reconstruction of the correct order


        ! Declare variables
        !==================
        ! Arguments
        class(UnstructuredInterpolant2DUDT)     :: interp 
        real(R8), intent(in)                    :: xg(:), yg(:)
        real(R8), intent(out)                   :: v(:)

        ! Auxiliary
        integer(I8) :: i, ic, order, min_number_of_terms, number_of_terms
        integer(I8), allocatable, dimension(:) :: ar, n_ar, tv
        real(R8), allocatable :: deriv_vals(:,:), A(:,:), xv(:), yv(:), Amask(:,:)
        logical, allocatable :: log(:)
        !type(GradientReconstructionTriaUDT) :: GR


        ! Checks
        if (size(xg, 1) /= size(yg, 1)) &
            call gdErrorHandler('ConstructUSI2DUSFinElem: size of xg and yg incompatible')
        if (size(xg) /= size(v)) &
            call gdErrorHandler('ConstructUSI2DUSFinElem: size of xg and v incompatible')
        if (interp%triangulation%nv /= size(xg)) &
            call gdErrorHandler('ConstructUSI2DUSFinElem: size of xg and triangulation incompatible')
        if (interp%C .gt. 3) &
            call gdErrorHandler('ConstructUSI2DUSFinEelem: higher than 3th order continuity not implemented')

        ! Save field information
        interp%v = v

        ! Determine needed order of derivatives
        ! This depends on the required continuity
        ! F.e. C = 3 => for triangles required 3 * 10 equations
        ! This would correspond with a 5th order interpolant because this has 36 terms
        ! while 4th order only has 25 terms and insufficient.
        print *, 'Constructing unstructured interpolation with C',interp%C, 'continuity on triangulated mesh'

        ! Compute number of terms
        ar = (/(i, i = 1, 10)/)
        n_ar = (ar + 1)**2
        min_number_of_terms = sum(ar(1:interp%C+1)) * 3 ! number of equation

        ! Determine neccesary order of interpolant
        ! (number of term = (N + 1)**2) for interpolant phi = sum_i^N sum_j^N a_ij x^i y^j
        log = [min_number_of_terms .lt. n_ar]
        order = findloc(log, .true., 1)
        number_of_terms = (order + 1)**2
        interp%n = number_of_terms

        ! Continue computation if not yet done before for other field on the same grid
        if (.not. allocated(interp%A)) then

            ! Build b matrix
            ! Compute GR coefficients
            call interp%GR%SetParameters('vert', 'vert', 'global', order)
            call interp%GR%Construct(interp%triangulation)
            call interp%GR%Evaluate(v, deriv_vals)

            ! Build A matrix
            allocate(interp%A(interp%n, interp%n))

            select case (interp%base_func)
            case ('polynomial')

                ! Construct Amask
                call ConstructAmask(order, Amask)

                ! Loop over the cell
                do ic = 1, interp%triangulation%nc

                    ! Get vertices
                    tv = interp%triangulation%cvert(ic, :)

                    ! Get coordinates
                    xv  = interp%triangulation%x(tv)
                    yv  = interp%triangulation%x(tv)

                    xv = [2, 2, 2]
                    yv = [2, 2, 2]
                    
                    ! Construct Afull up to highest supported order
                    allocate(A(interp%n,interp%n))
                    call ConstructA(order, interp%n, xv, yv, Amask, A)

                end do




            case default
                call gdErrorHandler('ConstructUSI2DUSFinEelem: type of base_func not implemented')
            end select


        end if

    end subroutine

    subroutine EvaluateUnstructuredInterpolant2D(interp, xq, yq, derivx, derivy, vq)

        ! Description
        !============
        ! This is a wrapper for the evaluation of the unstructured
        ! 2D interpolant. The two options are 'barycentric' and 
        ! 'finite_element'

        ! Declare variables
        !==================
        ! Arguments
        class(UnstructuredInterpolant2DUDT)     :: interp 
        real(R8), intent(in)                    :: xq(:), yq(:)
        real(R8), intent(out)                   :: vq(:)
        integer(I8), intent(in)                 :: derivx, derivy

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

        select case (interp%meth)

        case ('barycentric')

            call EvaluationUnstructuredInterpolant2DBary(interp, xq, yq, derivx, derivy, vq)           

        case ('finite_element')

            ! TODO

        case default

            call gdErrorHandler('EvaluateStructuredInterpolant2D: methode not implemented')
            
        end select



    end subroutine

    subroutine EvaluationUnstructuredInterpolant2DBary(interp, xq, yq, derivx, derivy, vq)


        ! Declare variables
        !==================
        ! Arguments
        class(UnstructuredInterpolant2DUDT)     :: interp 
        real(R8), intent(in)                    :: xq(:), yq(:)
        real(R8), intent(out)                   :: vq(:)
        integer(I8), intent(in)                 :: derivx, derivy

        ! Auxiliary
        integer(I8) :: i, v1, v2, v3, ctri, ind
        integer(I8), allocatable, dimension(:)  :: vt1, vt2, vt3
        real(R8) :: lambda1, lambda2, lambda3
        real(R8), allocatable, dimension(:)     :: vx, vy, v, dist
        logical, allocatable, dimension(:)      :: in, on

        ! Loop over query points
        vt1 = interp%triangulation%cvert(:,1)
        vt2 = interp%triangulation%cvert(:,2)
        vt3 = interp%triangulation%cvert(:,3)
        vx = interp%triangulation%x
        vy = interp%triangulation%y
        v = interp%v
        do i = 1, size(xq, 1)

            ! Get triangle
            call InTriangle(vt1, vt2, vt3, vx, vy, xq(i), yq(i), in, on)
            ctri = findloc(in .or. on, .true., 1)

            ! Check if triangle was found
            if (ctri /= 0) then 
                                
                ! Get triangle vertices
                v1 = vt1(ctri)
                v2 = vt2(ctri)
                v3 = vt3(ctri)

                ! Interpolate
                !============
                ! Compute barycentric coordinates
                call Cart2Bary(xq(i), yq(i), vx(v1), vy(v1), vx(v2), vy(v2), vx(v3), &
                    vy(v3), lambda1, lambda2, lambda3)

                ! Sanity check
                if (lambda1 < 0_R8 .or. lambda2 < 0_R8 .or. lambda3 < 0_R8) then 
                    print *, 'EvaluateFromSaddlePoint: point should lie in triangle, ' // & 
                        'but negative barycentric coordinates present. May be a ' // &
                        'bug '
                end if 
                
                ! Compute value
                vq(i) = lambda1*v(v1) + lambda2*v(v2) + lambda3*v(v3)

            else if (interp%allowextrapolation) then

                ! No triangle found, some extrapolation needed
                ! Find nearest value point
                dist = sqrt((vx - xq(i))**2 + (vy - yq(i))**2)
                ind = minloc(dist, 1)
                vq(i) = v(ind)

            else 

                call gdErrorHandler('EvaluationUnstructuredInterpolant2DBary: triangle ' // & 
                    'could not be found and no extrapolation allowed, ' // &
                    'check if point actually lies in mesh or enable extrapolation.')

            end if

        end do      
        
    end subroutine

    subroutine EvaluateWrapper(interp, v, xq, yq, derivx, derivy, vq)

        ! Description
        !============
        ! Wrapper for Construction and Evaluation

        ! Declare variables
        !==================
        ! Arguments
        class(UnstructuredInterpolant2DUDT)     :: interp 
        real(R8), intent(in)                    :: v(:)
        real(R8), intent(in)                    :: xq(:), yq(:)
        integer(I8), intent(in)                 :: derivx, derivy
        real(R8), intent(out)                   :: vq(:)

        ! Auxiliary
        real(R8), allocatable, dimension(:) :: vg, xg, yg

        vg = v
        xg = interp%triangulation%x
        yg = interp%triangulation%y
        call interp%ConstructUnstructured(xg, yg, vg)
        call interp%Evaluate(xq, yq, derivx, derivy, vq)

    end subroutine
        
    subroutine DeallocateUnstructuredInterpolant2D(interp)

        ! Description
        !============
        ! Deallocate a fully unstructured interpolation
        class(UnstructuredInterpolant2DUDT) :: interp

        ! TODO
        call gdErrorHandler('DeallocateUnstructuredInterpolant2D: not implemented')

    end subroutine

    subroutine ConstructAmask(order, A)

        ! Description
        !============
        ! Construct Amask matrix to avoid incorrect terms in A matrix

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in) :: order
        integer(I8), intent(out) :: A(36, (order+1)**2)

        ! Auxiliary
        integer(I8) :: i, j, M, k


        ! Initialize
        A = 1

        ! Loop 
        M = order+1
        do j = 1, M
            do i = 1, M
                ! Column index
                k = (j-1)*M + i 
                
                if (j .le. 3) then
                    ! For fourth derivatives of y onward not automatically 0
                    (34:36, k) = 0  ! d4phidy4
                    if (j .le. 2) then
                        ! For third derivatives of y onward not automatically 0
                        A(22:24, k) = 0  ! d3phidy3
                        if (j .le. 1) then
                            ! From second derivatives of y onwards not automatically 0
                            A(13:15, k) = 0  ! d2phidy2
                            A(28:30, k) = 0  ! d3dphidxdy2
                        end if
                    end if
                end if

                if (i .le. 3) then
                    ! For fourth derivatives of y onward not automatically 0
                    A(31:33, k) = 0  ! d4phidx4
                    if (i .le. 2) then
                        ! For third derivatives of x onward not automatically 0
                        A(19:21, k) = 0  ! d3phidx3
                        if (i .le. 1) then
                            ! From second derivatives of x onwards not automatically 0
                            A(10:12, k) = 0  ! d2phidx2
                            A(25:27, k) = 0  ! d3dphidx2dy
                        end if
                    end if
                end if

                !if (j == 1) then
                !    ! From second derivatives of y onwards not automatically 0
                !    A(13:15, k) = 0  ! d2phidy2
                !    A(22:24, k) = 0  ! d3phidy3
                !    A(28:30, k) = 0  ! d3dphidxdy2
                !    A(34:36, k) = 0  ! d4phidy4
                !else if (j == 2) then
                !    ! For third derivatives of y onward not automatically 0
                !    A(22:24, k) = 0  ! d3phidy3
                !    A(34:36, k) = 0  ! d4phidy4
                !else if (j == 3) then
                !    ! For fourth derivatives of y onward not automatically 0
                !    A(34:36, k) = 0  ! d4phidy4  
                !end if

                !if (i == 1) then
                !    ! From second derivatives of x onwards not automatically 0
                !    A(10:12, k) = 0  ! d2phidx2
                !    A(19:21, k) = 0  ! d3phidx3
                !    A(25:27, k) = 0  ! d3dphidx2dy
                !    A(31:33, k) = 0  ! d4phidx4  
                !else if (i == 2) then
                !    ! For third derivatives of x onward not automatically 0
                !    A(19:21, k) = 0  ! d3phidx3
                !    A(31:33, k) = 0  ! d4phidx4 
                !else if (i == 3) then
                !    ! For fourth derivatives of y onward not automatically 0
                !    A(31:33, k) = 0  ! d4phidx4
                !end if
                
            end do
        end do

    end subroutine

    subroutine ConstructA(order, n, xv, yv, Amask, A)

        ! Description
        !============
        ! Construct the whole A matrix up to 5th order

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in) :: order, n
        real(R8), intent(in)    :: xv(3), yv(3)
        integer, intent(in)     :: Amask(36, (order+1)**2)
        real(R8), intent(out)   :: A(n, n) 

        ! Auxiliary
        integer(I8) :: i, j, M, k
        real(R8) :: Afull(36, (order+1)**2)

        ! Initialize
        Afull = 0

        ! Equation phi is constant
        M = order+1
        do j = 1, M
            do i = 1, M

                ! Column index
                k = (j-1)*M + i

                ! Equation phi
                Afull(1:3, k) = Amask(1:3, k)*xv**(i-1) * yv**(j-1)

                ! Equation dphidx 
                Afull(4:6, k) = Amask(4:6, k)*(i-1)*xv**(i-2) * yv**(j-1)

                ! Equation dphidy
                Afull(7:9, k) = Amask(7:9, k)*xv**(i-1) * (j-1)*yv**(j-2)

                ! Equation d2phidx2
                Afull(10:12, k) = Amask(10:12, k)*(i-2)*xv**(i-3) * yv**(j-1)

                ! Equation d2phidy2
                Afull(13:15, k) = Amask(13:15, k)*xv**(i-1) * (j-2)*yv**(j-3)

                ! Equation d2phidxdy
                Afull(16:18, k) = Amask(16:18, k)*(i-1)*xv**(i-2) * (j-1)*yv**(j-2)

                ! Equation d3phidx3
                Afull(19:21, k) = Amask(19:21, k)*(i-3)*xv**(i-4) * yv**(j-1)

                ! Equation d3phidy3
                Afull(22:24, k) = Amask(22:24, k)*xv**(i-1) * (j-3)*yv**(j-4)

                ! Equation d3dphidx2dy
                Afull(25:27, k) = Amask(25:27, k)*(i-2)*xv**(i-3) * (j-1)*yv**(j-2)

                ! Equation d3dphidxdy2
                Afull(28:30, k) = Amask(28:30, k)*(i-1)*xv**(i-2) * (j-2)*yv**(j-3)

                ! Equation d4phidx4
                Afull(31:33, k) = Amask(31:33, k)*(i-4)*xv**(i-5) * yv**(j-1)

                ! Equation d4phidy4
                Afull(34:36, k) = Amask(34:36, k)*xv**(i-1) * (j-4)*yv**(j-5)

                ! Multiply with mask
                Afull(:,k) = Afull(:,k)*Amask(:,k)

            end do
        end do

        call gdErrorHandler('ConstructA: problem with ')

        ! Take the correct slice
        A = Afull(1:n, 1:n)        
 
    end subroutine


end module