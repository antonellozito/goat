!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module provides distribution function types for the grid 
! deformation module. Different types of distribution methods are 
! implemented.

module DistributionFunction

    ! Initialize
    !============
    ! Load modules
    use mod_precision 
    use Interpolant
    use mod_linearsolverinterface
    use mod_polygon

    ! The usual
    implicit none 
    public 
    

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================! 

    ! Abstract type
    !==============
    type, abstract :: DistributionFunctionUDT 

        ! Description
        !============
        ! Provides the main functions any distribution type should
        ! possess. Each distribution type should have an evaluation 
        ! function that takes x, y data as input and outputs an array of
        ! the same size as the input data. 

        ! No specific fields have to be included

    contains

        ! Initialization is not a deferred procedure as this is
        ! input-dependent etc.

        ! Evaluation
        procedure(EvaluateDistributionFunctionINT), deferred :: Evaluate   

        ! Visualization 
        procedure   :: Visualize    => VisualizeDistributionFunction

    end type

    ! Regular distance function
    type, extends(DistributionFunctionUDT)  :: Structured2DDistanceDFUDT 

        ! Description
        !============
        ! Distance function based on a 2D structured interpolant 
        ! function F, which is assumed to yield the distance w.r.t. some 
        ! quantity. Can, for example, be an interpolant of the vessel
        ! boundary function or something alike. The function being 
        ! evaluated is: 
        !
        !   d = a0*exp(-|F(x, y)|/d0) + (1-exp(-|F(x, y)|/d0))*b0,
        !
        ! Where a0 is the desired value at the vessel boundary, b0 is
        ! the asymptotic value for |F(x, y)| -> infinity, and d0 is a decay
        ! length (in [m]). 
        

        ! Fields:
        type(StructuredInterpolant2DUDT)        :: F 
        real(R8)                                :: a0, b0, d0

    contains 

        ! Initialization routine
        procedure :: Initialize     => InitializeStructured2DDistanceDF

        ! Evaluation
        procedure :: Evaluate       => EvaluateStructured2DDistanceDF

    end type

    ! Polygonset and field based, 2D
    type, extends(DistributionFunctionUDT) :: Polygonset2DFieldDistanceDFUDT

        ! Description
        !============
        ! Distance function based on polygon set and 2D (vector) field. The field
        ! should be described using a 2D structured interpolant. The 
        ! polygon set has to be a set of closed, non-intersecting 
        ! (possibly nested) polygons, otherwise the distribution 
        ! function cannot be defined. By setting the 'meth' to 'signed' or
        ! 'unsigned', one can obtain an unsymmetric or symmetric 
        ! distribution based on the sign of the dot product between    
        ! polygon set normal and field normal (field normal is 
        ! determined as [-dFdy, dFdx] where F is the scalar value of 
        ! the field and dFdx, dFdy the gradient)

        ! Fields:
        type(StructuredInterpolant2DUDT)        :: F  ! field
        type(PolygonSetUDT)                     :: PS ! polygon set
        real(R8)                                :: b0, a0, d0

        real(R8), allocatable                   :: xa(:), ya(:), nxa(:), &
            nya(:), d(:), coef(:), dp(:)

        character(:), allocatable               :: meth
        

    contains

        ! Initialization
        procedure :: Initialize     => InitializePolygonset2DFieldDistanceDF

        ! Evaluation
        procedure :: Evaluate       => EvaluatePolygonset2DFieldDistanceDF

    end type

    ! Polygonset and field based, 1D 
    type, extends(DistributionFunctionUDT)  :: Polygonset1DFieldDistanceDFUDT

        ! Description
        !============
        ! Similar distance function as the 2D equivalent above, but now
        ! we construct a distribution for 1D values and use only the
        ! value of the 2D interpolant. Useful for distance to certain 
        ! psi values for example. Do note that it is a scalar field, yet
        ! it is evaluated in two dimensions. Analogous to the 2D 
        ! type, one can use 'meth' to have an asymmetric ('signed') or
        ! symmetric ('unsigned') distribution. 

        ! Note: a0 is now the desired value at the polygonset vertices. 
        ! Furthermore, the polygon does not have to be closed whatsoever
        ! - we simply determine the psi values at the polygon nodes     
        ! and construct attractors to those values. Likewise, b0 is the 
        ! desired distribution value far way from these locations. 

        ! Fields:
        type(StructuredInterpolant2DUDT)        :: F  ! field
        type(PolygonSetUDT)                     :: PS ! polygon set
        real(R8)                                :: b0, a0, d0

        real(R8), allocatable                   :: xa(:), ya(:), &
            coef(:), fval(:)

        character(:), allocatable               :: meth

    contains 

        ! Initialization
        procedure :: Initialize     => InitializePolygonset1DFieldDistanceDF

        ! Evaluation
        procedure :: Evaluate       => EvaluatePolygonset1DFieldDistanceDF

    end type

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================! 

    ! Abstract interfaces
    !====================
    abstract interface 

        ! Evaluation routine
        subroutine EvaluateDistributionFunctionINT(distribution, x, y, v)
            import :: DistributionFunctionUDT, R8
            class(DistributionFunctionUDT) :: distribution
            real(R8), intent(in)        :: x(:), y(:)
            real(R8), intent(out)       :: v(size(x))
        end subroutine

    end interface

    contains 

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                               GENERAL                            !
    !------------------------------------------------------------------!

    ! Visualization
    subroutine VisualizeDistributionFunction(distribution, &
        xrange, yrange, resx, resy, savefilepath)

        ! Description
        !============
        ! Visualize the distribution function by evaluating it on a 
        ! structured 2D mesh. Coordinates and values are written out 
        ! in z = f(x, y) format, which can be used later on for post-
        ! processing

        ! Declare variables
        !==================
        ! Arguments
        class(DistributionFunctionUDT)  :: distribution 
        real(R8), intent(in)            :: xrange(1:2), yrange(1:2)
        integer(I8), intent(in)         :: resx, resy 
        character(*), intent(in)        :: savefilepath 

        ! Auxiliary
        real(R8)                        :: dx, dy 
        real(R8), allocatable           :: xgv(:), ygv(:), xg(:), &
            yg(:), vg(:)

        ! Loop
        integer(I8)                     :: k 

        ! Initialize
        !===========
        ! Construct gridding vectors
        dx = maxval(xrange) - minval(xrange)
        dy = maxval(yrange) - minval(yrange)
        xgv = [(k, k = 0, resx)]*dx/resx + minval(xrange)
        ygv = [(k, k = 0, resy)]*dy/resy + minval(yrange)

        ! Construct grid
        allocate(xg((resx+1)*(resy+1)), yg((resx+1)*(resy+1)), vg((resx+1)*(resy+1)))
        call Construct2DStructuredGrid(xgv, ygv, resx+1, resy+1, xg, yg)

        ! Evaluate
        call distribution%Evaluate(xg, yg, vg)

        ! Write data
        call Write3DCoordinateData(xg, yg, vg, savefilepath)

    end subroutine

    !------------------------------------------------------------------!
    !                         DISTANCE FUNCTION                        !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeStructured2DDistanceDF(distribution, interp, &
        val0, valinf, decaylength)

        ! Description
        !============
        ! Initialization routine. The 'interp' structured interpolant
        ! must be initialized and correctly set up. The argument 'val0'
        ! is the value that the function achieves where the interpolant
        ! is (approximately) zero, the value 'valinf' is achieved at 
        ! regions where the interpolant approaches +-infinity. 
        ! 'decaylength' is a decay length of how fast val0 transitions 
        ! to valinf. 

        ! Declare variables
        !==================
        ! Arguments
        class(Structured2DDistanceDFUDT)                :: distribution 
        type(StructuredInterpolant2DUDT), intent(in)    :: interp 
        real(R8), intent(in)                            :: val0, valinf, &
            decaylength

        ! Set fields
        !===========
        distribution%F  = interp 
        distribution%a0 = val0
        distribution%b0 = valinf 
        distribution%d0 = decaylength

    end subroutine

    ! Evaluation
    subroutine EvaluateStructured2DDistanceDF(distribution, x, y, v)

        ! Description
        !============
        ! Evaluate the distribution function

        ! Declare variables
        !==================
        ! Arguments
        class(Structured2DDistanceDFUDT)    :: distribution 
        real(R8), intent(in)                :: x(:), y(:)
        real(R8), intent(out)               :: v(size(x))

        ! Auxiliary
        real(R8)                            :: d(size(x))

        ! Initialize
        !===========
        ! Check sizes
        if ( (size(v) /= size(x)) .or. (size(x) /= size(y))) then 
            ! Throw error
            call gdErrorHandler('EvaluateStructured2DDistanceDF: incompatible sizes in input')
        end if 

        ! Associate
        associate(&
            a0      => distribution%a0,     & 
            b0      => distribution%b0,     & 
            d0      => distribution%d0,     & 
            F       => distribution%F       & 
        )

        ! Evaluate
        !=========
        ! Distance 
        call F%Evaluate(x, y, 0, 0, d)

        ! Value
        v = a0*exp(-abs(d)/d0) + b0*(1 - exp(-abs(d)/d0))

        ! Housekeeping
        !=============
        end associate


    end subroutine

    !------------------------------------------------------------------!
    !                      POLYGONSET & 2D FIELD                       !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializePolygonset2DFieldDistanceDF(distribution, interp, &
        ps, val0, valinf, decaylength, meth)

        ! Description
        !============
        ! Initialization routine. The 'interp' structured interpolant
        ! must be initialized and correctly set up. The argument 'val0'
        ! is the value that the function achieves where the interpolant
        ! is (approximately) zero, the value 'valinf' is achieved at 
        ! regions where the interpolant approaches +-infinity. 
        ! 'decaylength' is a decay length of how fast val0 transitions 
        ! to valinf. 

        ! Declare variables
        !==================
        ! Arguments
        class(Polygonset2DFieldDistanceDFUDT)           :: distribution 
        type(StructuredInterpolant2DUDT), intent(in)    :: interp 
        type(PolygonSetUDT), intent(in)                 :: ps
        real(R8), intent(in)                            :: val0, valinf, &
            decaylength
        character(*), intent(in)                        :: meth

        ! Auxiliary
        integer(I8)                                     :: flag, na
        integer(I8), allocatable                        :: valindex(:) 

        real(R8)                                        :: tempd
        real(R8), allocatable                           :: xa(:), &
            ya(:), nxa(:), nya(:), b(:), bval(:), A(:, :), &
            tnxa(:), tnya(:), tnna(:), sol(:), dFdx(:), dFdy(:)

        ! Loop
        integer(I8)                                     :: i, j, k, cc

        ! Set fields
        !===========
        ! Data
        distribution%F  = interp 
        distribution%PS = ps
        distribution%b0 = val0
        distribution%a0 = valinf 
        distribution%d0 = decaylength
        distribution%meth = trim(meth)

        ! Associate
        associate(&
            a0      => distribution%a0,     &
            b0      => distribution%b0,     &
            d0      => distribution%d0)

        ! Check
        !======
        ! Check if the method is either 'signed' or 'unsigned'
        if ((meth /= 'signed') .and. (meth /= 'unsigned')) then 
            ! Throw error
            call gdErrorHandler('PolygonsetField2DDistanceDF: ' // &
                'meth should be either "signed" or "unsigned" ')
        end if 

        ! Call orientation routine, if that one fails -> can't make
        ! distribution. Afterwards, all normals should point outward
        ! the domain
        call ps%OrientNestedClosedPolygons(flag)
        if (flag > 0) then 
            ! Throw error
            call gdErrorHandler('PolygonsetField2DDistanceDF: ' //&
                'polygonset is not a set of closed, non-intersecting polygons')
        end if 


        ! Construct attractor function
        !=============================
        ! Determine number of attractor points
        na = 0
        do i = 1, ps%np 
            na = na + ps%polygons(i)%nv
        end do 

        ! Allocate
        allocate(xa(na), ya(na), nxa(na), nya(na), bval(na), &
            b(na), A(na, na))

        ! Determine
        cc = 0
        do i = 1, ps%np 
            ! Get point indices
            valindex = [(k, k = cc+1, cc+ps%polygons(i)%nv)]

            ! Set points
            xa(valindex) = ps%polygons(i)%x
            ya(valindex) = ps%polygons(i)%y

            ! Determine normal in points
            allocate(tnxa(ps%polygons(i)%nv), tnya(ps%polygons(i)%nv))
            do j = 1, 2 
                tnxa = tnxa + ps%polygons(i)%nx/ps%polygons(i)%nn 
                tnya = tnya + ps%polygons(i)%ny/ps%polygons(i)%nn
            end do

            ! Rescale
            tnna = sqrt(tnxa**2 + tnya**2)
            tnxa = tnxa/tnna 
            tnya = tnya/tnna

            ! Add
            nxa(valindex) = tnxa 
            nya(valindex) = tnya

            ! Housekeeping
            deallocate(tnxa, tnya)

            ! Update counter
            cc = cc + ps%polygons(i)%nv
        end do

        ! Construct rhs to compute attractor coefficients
        b = a0 - b0

        ! Compute lhs to compute attractor coefficients
        A = 0
        do j = 1, na
            do i = 1, na
                if (i /= 0) then 
                    tempd = sqrt((xa(i) - xa(j))**2 + (ya(i) - ya(j))**2)
                    A(i, j) = exp(-tempd/d0)
                else 
                    A(i, j) = 1
                end if 
            end do 
        end do

        ! Call solver
        allocate(sol(size(b)))
        call SolveDenseLinearSystemDI(A, b, sol, flag)
        if (flag /= 0) then
            ! Call error
            call gdErrorHandler('InitializePolygonsetField2DDistanceDF: ' // &
                'could not determine attractor function coefficients ' // &
                'due to non-converging linear solver')
        end if 

        ! Evaluate field
        allocate(dFdx(na), dFdy(na))
        call distribution%F%Evaluate(xa, ya, 1, 0, dFdx)
        call distribution%F%Evaluate(xa, ya, 0, 1, dFdy)

        ! Add
        !====
        distribution%xa     = xa 
        distribution%ya     = ya 
        distribution%nxa    = nxa 
        distribution%nya    = nya 
        distribution%d      = spread(d0, 1, size(xa))
        distribution%coef   = sol 
        distribution%dp = - (-dFdy*nxa + dFdx*nya)

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Evaluation
    subroutine EvaluatePolygonset2DFieldDistanceDF(distribution, x, y, v)

        ! Description
        !============
        ! Evaluate the distribution function

        ! Declare variables
        !==================
        ! Arguments
        class(Polygonset2DFieldDistanceDFUDT)   :: distribution 
        real(R8), intent(in)                    :: x(:), y(:)
        real(R8), intent(out)                   :: v(size(x))

        ! Auxiliary
        real(R8), parameter                     :: myone = 1
        real(R8)                                :: d(size(x))

        ! Loop
        integer(I8)                             :: i

        ! Initialize
        !===========
        ! Check sizes
        if ( (size(v) /= size(x)) .or. (size(x) /= size(y))) then 
            ! Throw error
            call gdErrorHandler('EvaluatePolygonsetField2DDistanceDF: incompatible sizes in input')
        end if 

        ! Associate
        associate(&
            a0      => distribution%a0,     & 
            b0      => distribution%coef,   & 
            d0      => distribution%d,      &
            xa      => distribution%xa,     &
            ya      => distribution%ya,     &
            nxa     => distribution%nxa,    &
            nya     => distribution%nya,    &
            coef    => distribution%coef,   &
            dp      => distribution%dp      &
        )

        ! Evaluate
        !=========
        v = 0
        select case (distribution%meth) 

        case ('unsigned')

            do i = 1, size(xa)
                ! Distance 
                d = sqrt( (x - xa(i))**2 + (y - ya(i))**2)

                ! Value
                v = v + b0(i)*exp(-d/d0(i))
            end do

        case ('signed')

            do i = 1, size(xa)
                ! Distance 
                d = sqrt( (x - xa(i))**2 + (y - ya(i))**2)

                ! Value
                v = v + b0(i)*(sign(myone, dp(i)))*exp(-d/d0(i))
            end do

        case default

            call gdErrorHandler('Unknown method')

        end select

        v = a0 - v

        ! Housekeeping
        !=============
        end associate


    end subroutine 

    !------------------------------------------------------------------!
    !                      POLYGONSET & FIELD, 1D                      !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializePolygonset1DFieldDistanceDF(distribution, interp, &
        ps, val0, valinf, decaylength, meth)

        ! Description
        !============
        ! Initialization routine. The 'interp' structured interpolant
        ! must be initialized and correctly set up. The argument 'val0'
        ! is the value that the function achieves where the interpolant
        ! is (approximately) zero, the value 'valinf' is achieved at 
        ! regions where the interpolant approaches +-infinity. 
        ! 'decaylength' is a decay length of how fast val0 transitions 
        ! to valinf. 

        ! Declare variables
        !==================
        ! Arguments
        class(Polygonset1DFieldDistanceDFUDT)           :: distribution 
        type(StructuredInterpolant2DUDT), intent(in)    :: interp 
        type(PolygonSetUDT), intent(in)                 :: ps
        real(R8), intent(in)                            :: val0, valinf, &
            decaylength
        character(*), intent(in)                        :: meth

        ! Auxiliary
        integer(I8)                                     :: flag, na
        integer(I8), allocatable                        :: valindex(:) 

        real(R8), parameter                             :: myone = 1
        real(R8)                                        :: tempd
        real(R8), allocatable                           :: xa(:), &
            ya(:), nxa(:), nya(:), d(:), b(:), bval(:), A(:, :), &
           sol(:), fval(:)

        ! Loop
        integer(I8)                                     :: i, j, k, cc

        ! Set fields
        !===========
        ! Data
        distribution%F  = interp 
        distribution%PS = ps
        distribution%a0 = val0
        distribution%b0 = valinf 
        distribution%d0 = decaylength
        distribution%meth = trim(meth)

        ! Associate
        associate(&
            a0      => distribution%a0,     &
            b0      => distribution%b0,     &
            d0      => distribution%d0)

        ! Check
        !======
        ! Check if the method is either 'signed' or 'unsigned'
        if ((meth /= 'signed') .and. (meth /= 'unsigned')) then 
            ! Throw error
            call gdErrorHandler('PolygonsetField2DDistanceDF: ' // &
                'meth should be either "signed" or "unsigned" ')
        end if 

        ! Construct attractor function
        !=============================
        ! Determine number of attractor points
        na = 0
        do i = 1, ps%np 
            na = na + ps%polygons(i)%nv
        end do 

        ! Allocate
        allocate(xa(na), ya(na), nxa(na), nya(na), d(na), bval(na), &
            b(na), A(na, na))

        ! Determine
        cc = 0
        do i = 1, ps%np 
            ! Get point indices
            valindex = [(k, k = cc+1, cc+ps%polygons(i)%nv)]

            ! Set points
            xa(valindex) = ps%polygons(i)%x
            ya(valindex) = ps%polygons(i)%y

            ! Update counter
            cc = cc + ps%polygons(i)%nv
        end do

        ! Field values
        allocate(fval(na))
        call distribution%F%Evaluate(xa, ya, 0, 0, fval)

        distribution%fval = fval 

        ! Construct rhs to compute attractor coefficients
        b = a0 - b0

        ! Compute lhs to compute attractor coefficients
        A = 0
        do j = 1, na
            do i = 1, na
                if (i /= 0) then 
                    tempd = abs(fval(i) - fval(j))
                    A(i, j) = sign(myone, fval(i)-fval(j))*exp(-tempd/d0)
                else 
                    A(i, j) = 1
                end if 
            end do 
        end do

        ! Call solver
        allocate(sol(size(b)))
        call SolveDenseLinearSystemDI(A, b, sol, flag)
        if (flag /= 0) then
            ! Call error
            call gdErrorHandler('InitializePolygonsetField2DDistanceDF: ' // &
                'could not determine attractor function coefficients ' // &
                'due to non-converging linear solver')
        end if 

        ! Add
        !====
        distribution%xa     = xa 
        distribution%ya     = ya 
        distribution%fval   = fval
        distribution%coef   = sol 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Evaluation
    subroutine EvaluatePolygonset1DFieldDistanceDF(distribution, x, y, v)

        ! Description
        !============
        ! Evaluate the distribution function

        ! Declare variables
        !==================
        ! Arguments
        class(Polygonset1DFieldDistanceDFUDT)   :: distribution 
        real(R8), intent(in)                    :: x(:), y(:)
        real(R8), intent(out)                   :: v(size(x))

        ! Auxiliary
        real(R8), parameter                     :: myone = 1
        real(R8)                                :: fv(size(x)), d(size(x))

        ! Loop
        integer(I8)                             :: i

        ! Initialize
        !===========
        ! Check sizes
        if ( (size(v) /= size(x)) .or. (size(x) /= size(y))) then 
            ! Throw error
            call gdErrorHandler('EvaluatePolygonsetField2DDistanceDF: incompatible sizes in input')
        end if 

        ! Associate
        associate(&
            F       => distribution%F,      &
            a0      => distribution%a0,     & 
            b0      => distribution%coef,   & 
            d0      => distribution%d0,     &
            xa      => distribution%xa,     &
            ya      => distribution%ya,     &
            fval    => distribution%fval    &
        )

        ! Evaluate
        !=========
        ! Evaluate field values in coordinates
        call F%Evaluate(x, y, 0, 0, fv)
        select case (distribution%meth) 

        case ('unsigned')

            do i = 1, size(xa)
                ! Distance 
                d = (fv - fval(i))

                ! Value
                v = v + b0*exp(-abs(d)/d0)
            end do

        case ('signed')

            do i = 1, size(xa)
                ! Distance 
                d = (fv - fval(i))

                ! Value
                v = v + b0*(sign(myone, d))*exp(-abs(d)/d0)
            end do

        case default

            call gdErrorHandler('Unknown method')

        end select

        ! Housekeeping
        !=============
        end associate


    end subroutine 



    

end module