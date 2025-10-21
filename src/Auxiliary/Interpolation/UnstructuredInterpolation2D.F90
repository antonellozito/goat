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
        ! - allowextrapolation: logical to check if we can extrapolate
        ! - triangulation: triangulated mesh to use for interpolation
        ! - v: the values at the vertex points
        ! - xgv, ygv: grid vectors
        ! - A: interpolation coefficients
        ! - refx, refy, refdx, refdy: reference values used to compute
        ! derivatives etc 
        ! - cellindex: nx-1 by ny-1 array containing the cell indices
        ! - n: number of terms in the interpolant
        ! - precomputedfac: the required factorials precomputed to save
        ! some time during evaluation 

        character(:), allocatable       :: meth 
        integer(I8)                     :: C, M, n
        logical                         :: allowextrapolation
        type(TriangulationUDT)          :: triangulation

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
    subroutine SetParametersUS(interp, meth, C, M, triangulation)

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
        type(TriangulationUDT), intent(in)      :: triangulation

        ! Set
        !====
        interp%meth = meth 
        interp%C    = C 
        interp%M    = M
        interp%allowextrapolation = .true.
        interp%triangulation = triangulation
        
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
        if (size(interp%triangulation%x) /= size(xg)) &
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
        ! 0) Compute the required derivatives on the vertex nodes.
        ! TODO

        ! Declare variables
        !==================
        ! Arguments
        class(UnstructuredInterpolant2DUDT)     :: interp 
        real(R8), intent(in)                    :: xg(:), yg(:)
        real(R8), intent(out)                   :: v(:)

        ! Auxiliary
        real(R8), allocatable :: deriv_vals(:,:)
        type(GradientReconstructionTriaUDT) :: GR


        ! Checks
        if (size(xg, 1) /= size(yg, 1)) &
            call gdErrorHandler('ConstructUSI2DUSFinElem: size of xg and yg incompatible')
        if (size(xg) /= size(v)) &
            call gdErrorHandler('ConstructUSI2DUSFinElem: size of xg and v incompatible')
        if (size(interp%triangulation%x) /= size(xg)) &
            call gdErrorHandler('ConstructUSI2DUSFinElem: size of xg and triangulation incompatible')

        ! Save field information
        interp%v = v

        ! Test GR results
        call GR%SetParameters('vert', 'vert', 'global', interp%M)
        call GR%Construct(interp%triangulation)
        call GR%Evaluate(v, deriv_vals)


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


end module