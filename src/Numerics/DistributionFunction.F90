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
    use PolygonLevelsetFunction2D
    use mod_structured2Dgridding
    use omp_lib

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

        ! Derivative evaluation
        procedure(EvaluateDerivativeDistributionFunctionINT), deferred :: &
            EvaluateDerivative

        ! Visualization 
        procedure   :: Visualize    => VisualizeDistributionFunction

    end type

    ! Simple field evaluation based on interpolant
    type, extends(DistributionFunctionUDT) :: Structured2DDFUDT

        ! Description
        !============
        ! Distribution function that serves as a wrapper for a 2D
        ! structured interpolant. May be usefule in some cases.

        ! Fields
        type(StructuredInterpolant2DUDT)    :: F 

    contains 

        ! Evaluation
        procedure :: Evaluate       => EvaluateStructured2DDF

        ! Derivative evaluation
        procedure :: EvaluateDerivative  => EvaluateDerivativeStructured2DDF

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

        ! Derivative evaluation
        procedure :: EvaluateDerivative     => EvaluateDerivativeStructured2DDistanceDF

    end type

    type, extends(DistributionFunctionUDT)  :: StructuredPLF2DDistanceDFUDT 

        ! Description
        !============
        ! Distance function based on a 2D structured interpolant 
        ! function F, which is assumed to yield the distance w.r.t. some 
        ! quantity. Can, for example, be an interpolant of the vessel
        ! boundary function or something alike. The function being 
        ! evaluated is: 
        !
        !   d = a0*dp*exp(-|F(x, y)|/d0) + (1-dp*exp(-|F(x, y)|/d0))*b0,
        !
        ! Where a0 is the desired value at the vessel boundary, b0 is
        ! the asymptotic value for |F(x, y)| -> infinity, and d0 is a decay
        ! length (in [m]). dp is the dot product between the normal of V 
        ! and the tangent of F
        

        ! Fields:
        type(StructuredInterpolant2DUDT)                    :: F
        class(PolygonLevelsetFunction2DUDT), allocatable    :: Vd, Vn
        real(R8)                                :: a0, b0, d0
        character(:), allocatable               :: meth

    contains 

        ! Initialization routine
        procedure :: Initialize     => InitializeStructuredPLF2DDistanceDF

        ! Evaluation
        procedure :: Evaluate       => EvaluateStructuredPLF2DDistanceDF

        ! Derivative evaluation
        procedure :: EvaluateDerivative     => EvaluateDerivativeStructuredPLF2DDistanceDF

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

        ! Derivative evaluation
        procedure :: EvaluateDerivative     => EvaluateDerivativePolygonset2DFieldDistanceDF

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

        ! Derivative evaluation
        procedure :: EvaluateDerivative     => EvaluateDerivativePolygonset1DFieldDistanceDF

    end type

    ! Coordinates, 1D
    type, extends(DistributionFunctionUDT)  :: Coordinates1DFieldDistanceDFUDT
        
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
        real(R8)                                :: b0, a0, d0

        real(R8), allocatable                   :: xa(:), ya(:), &
            coef(:), fval(:)

        character(:), allocatable               :: meth

    contains 

        ! Initialization
        procedure :: Initialize     => InitializeCoordinates1DFieldDistanceDF

        ! Evaluation
        procedure :: Evaluate       => EvaluateCoordinates1DFieldDistanceDF

        ! Derivative evaluation
        procedure :: EvaluateDerivative     => EvaluateDerivativeCoordinates1DFieldDistanceDF

    end type

    ! Coordinates, 2D
    type, extends(DistributionFunctionUDT)  :: Coordinates2DDistanceDFUDT
        
        ! Description
        !============
        ! Simple distribution function based on decaying exponentials
        ! defined in the set points xa, ya (at these points, the value
        ! will be equal to fval). Decay lengths can be set to determine
        ! how fast the value at those points approaches the value at    
        ! infinitely far away from these points. 

        ! Fields:
        real(R8)                                :: b0

        real(R8), allocatable                   :: xa(:), ya(:), &
            coef(:), d0(:), a0(:)

        character(:), allocatable               :: meth

    contains 

        ! Initialization
        procedure :: Initialize     => InitializeCoordinates2DDistanceDF

        ! Evaluation
        procedure :: Evaluate       => EvaluateCoordinates2DDistanceDF

        ! Derivative evaluation
        procedure :: EvaluateDerivative     => EvaluateDerivativeCoordinates2DDistanceDF

    end type

    ! Coordinates and PLF, 2D
    type, extends(DistributionFunctionUDT)  :: CoordinatesPLF2DDistanceDFUDT
        
        ! Description
        !============
        ! Distribution function based on coordinates (like 
        ! Coordiantes2DDistanceDFUDT), but now there's an additional
        ! background distribution based on a polygon levelset function.
        ! The function being evaluated is:
        !
        !   F(x, y) = sum_i (a_i exp(-d(x, y, xi, yi)/d_i)) + b + 
        !               a_plf*exp(-d_plf(x, y)/d_plf)
        !   
        ! Here, the first term originates from the specified points, 
        ! where for each point the distance is taken and divided through
        ! a decay distance d_i. b is the value far away from the plf and
        ! point. The last term is then the contribution of the plf, and,
        ! in absence of any points, a_plf then represents the value at
        ! zero distance. 

        ! Fields:
        real(R8)                                :: b0, a_plf, d_plf
        real(R8), allocatable                   :: xa(:), ya(:), &
            coef(:), d0(:), a0(:)
        class(PolygonLevelsetFunction2DUDT), allocatable    :: plf 

        character(:), allocatable               :: meth

    contains 

        ! Initialization
        procedure :: Initialize     => InitializeCoordinatesPLF2DDistanceDF

        ! Evaluation
        procedure :: Evaluate       => EvaluateCoordinatesPLF2DDistanceDF

        ! Derivative evaluation
        procedure :: EvaluateDerivative  => EvaluateDerivativeCoordinatesPLF2DDistanceDF

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

        ! Derivative evaluation routine
        subroutine EvaluateDerivativeDistributionFunctionINT(distribution, x, y, derivx, derivy, v)
            import :: DistributionFunctionUDT, R8, I8
            class(DistributionFunctionUDT) :: distribution
            real(R8), intent(in)        :: x(:), y(:)
            integer(I8), intent(in)     :: derivx, derivy
            real(R8), intent(out)       :: v(size(x))
        end subroutine

    end interface

    ! Assignment override
    interface assignment(=)
        module procedure AssignDFClass
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
        allocate(xgv(resx+1), ygv(resy+1))
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

    ! Assignment
    subroutine AssignDFClass(a, b)

        class(DistributionFunctionUDT), allocatable, intent(inout)    :: a 
        class(DistributionFunctionUDT), intent(in)                    :: b 

        if (allocated(a)) then 
            deallocate(a)
        end if 

        select type (b)

        class default 

            call gdErrorHandler('Unknown type')

        type is (Structured2DDFUDT)

            allocate(a, source=b)
            select type (a)
            type is (Structured2DDFUDT)
                a = b 
            end select

        type is (Structured2DDistanceDFUDT)

            allocate(a, source=b)
            select type (a)
            type is (Structured2DDistanceDFUDT)
                a = b 
            end select

        type is (StructuredPLF2DDistanceDFUDT)

            allocate(a, source=b)
            select type (a)
            type is (StructuredPLF2DDistanceDFUDT)
                a = b 
            end select

        type is (Polygonset2DFieldDistanceDFUDT)

            allocate(a, source=b)
            select type (a)
            type is (Polygonset2DFieldDistanceDFUDT)
                a = b
            end select

        type is (Polygonset1DFieldDistanceDFUDT)

            allocate(a, source=b)
            select type (a)
            type is (Polygonset1DFieldDistanceDFUDT)
                a = b
            end select

        type is (Coordinates1DFieldDistanceDFUDT)

            allocate(a, source=b)
            select type (a)
            type is (Coordinates1DFieldDistanceDFUDT)
                a = b
            end select

        type is (Coordinates2DDistanceDFUDT)

            allocate(a, source=b)
            select type (a)
            type is (Coordinates2DDistanceDFUDT)
                a = b
            end select

        type is (CoordinatesPLF2DDistanceDFUDT)

            allocate(a, source=b)
            select type (a)
            type is (CoordinatesPLF2DDistanceDFUDT)
                a = b
            end select

        end select 

    
    end subroutine

    !------------------------------------------------------------------!
    !                     2D STRUCTURED INTERPOLANT                    !
    !------------------------------------------------------------------!

    ! Constructor
    function ConstructStructured2DDF(interp) result(distribution)

        ! Description
        !============
        ! Construct the distributor based on the given interpolant

        ! Declare variables
        !==================
        ! Arguments
        class(StructuredInterpolant2DUDT), intent(in)   :: interp 
        class(DistributionFunctionUDT), allocatable     :: distribution 

        ! Initialize
        !===========
        allocate(Structured2DDFUDT::distribution)

        select type(distribution)

        type is (Structured2DDFUDT)

            ! Add interpolant
            distribution%F = interp

        end select

    end function 

    ! Evaluation
    subroutine EvaluateStructured2DDF(distribution, x, y, v)

        ! Description
        !============
        ! Evaluate the distribution function

        ! Declare variables
        !==================
        ! Arguments
        class(Structured2DDFUDT)            :: distribution 
        real(R8), intent(in)                :: x(:), y(:)
        real(R8), intent(out)               :: v(size(x))

        ! Evaluate
        !=========
        ! Just call interpolant evaluator
        call distribution%F%Evaluate(x, y, 0, 0, v)

    end subroutine 

    ! Derivative evaluation
    subroutine EvaluateDerivativeStructured2DDF(distribution, x, y, &
        derivx, derivy, v)

        ! Description
        !============
        ! Evaluate the distribution function derivative

        ! Declare variables
        !==================
        ! Arguments
        class(Structured2DDFUDT)            :: distribution 
        real(R8), intent(in)                :: x(:), y(:)
        real(R8), intent(out)               :: v(size(x))
        integer(I8), intent(in)             :: derivx, derivy 

        ! Evaluate
        !=========
        ! Just call interpolant evaluator
        call distribution%F%Evaluate(x, y, derivx, derivy, v)

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

    ! Derivative evaluation
    subroutine EvaluateDerivativeStructured2DDistanceDF(distribution, x, y, &
        derivx, derivy, v)

        ! Description
        !============
        ! Evaluate the distribution function

        ! Declare variables
        !==================
        ! Arguments
        class(Structured2DDistanceDFUDT)    :: distribution 
        real(R8), intent(in)                :: x(:), y(:)
        real(R8), intent(out)               :: v(size(x))
        integer(I8), intent(in)             :: derivx, derivy

        call gdErrorHandler('EvaluateDerivativeStructured2DDistanceDF: ' // & 
            'derivatives not yet implemented for this distribution type')
    
    end subroutine

    !------------------------------------------------------------------!
    !                       DISTANCE PLF FUNCTION                      !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeStructuredPLF2DDistanceDF(distribution, F, Vd, &
        Vn, val0, valinf, decaylength, meth)

        ! Description
        !============
        ! Initialization routine. The F and V structured interpolants
        ! must be initialized and correctly set up. The argument 'val0'
        ! is the value that the function achieves where the interpolant
        ! is (approximately) zero, the value 'valinf' is achieved at 
        ! regions where the interpolant approaches +-infinity. 
        ! 'decaylength' is a decay length of how fast val0 transitions 
        ! to valinf. 

        ! Declare variables
        !==================
        ! Arguments
        class(StructuredPLF2DDistanceDFUDT)             :: distribution 
        type(StructuredInterpolant2DUDT), intent(in)    :: F
        class(PolygonLevelsetFunction2DUDT), intent(in) :: Vd, Vn 
        real(R8), intent(in)                            :: val0, valinf, &
            decaylength
        character(*), intent(in)                        :: meth

        ! Check
        !======
        ! Check if the method is either 'signed' or 'unsigned'
        if ((meth /= 'signed') .and. (meth /= 'unsigned')) then 
            ! Throw error
            call gdErrorHandler('InitializeStructuredPLF2DDistanceDF: ' // &
                'meth should be either "signed" or "unsigned" ')
        end if 

        ! Set fields
        !===========
        distribution%F  = F
        distribution%Vd = Vd
        distribution%Vn = Vn 
        distribution%a0 = val0
        distribution%b0 = valinf 
        distribution%d0 = decaylength
        distribution%meth = meth 

    end subroutine

    ! Evaluation
    subroutine EvaluateStructuredPLF2DDistanceDF(distribution, x, y, v)

        ! Description
        !============
        ! Evaluate the distribution function

        ! Declare variables
        !==================
        ! Arguments
        class(StructuredPLF2DDistanceDFUDT) :: distribution 
        real(R8), intent(in)                :: x(:), y(:)
        real(R8), intent(out)               :: v(size(x))

        ! Auxiliary
        real(R8)                            :: d(size(x)), &
            dVdx(size(x)), dVdy(size(x)), dFdx(size(x)), dFdy(size(x)), &
            dp(size(x)), myones(size(x))

        ! Initialize
        !===========
        ! set
        myones = 1

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
            Vd      => distribution%Vd,     & 
            Vn      => distribution%Vn,     & 
            F       => distribution%F       & 
        )

        ! Evaluate
        !=========
        ! Distance 
        call Vd%Evaluate(x, y, 0, 0, d)
        
        
        ! Dot product
        select case (distribution%meth)

        case ('signed')

            ! Normal
            call Vn%Evaluate(x, y, 1, 0, dVdx)
            call Vn%Evaluate(x, y, 0, 1, dVdy)

            ! Vector field
            call F%Evaluate(x, y, 1, 0, dFdx)
            call F%Evaluate(x, y, 0, 1, dFdy)

            ! Dot product
            dp = -dVdx*dFdy + dVdy*dFdx 

        case ('unsigned')

            dp = spread(1, 1, size(x, 1))

        end select 

        ! Value
        v = a0*sign(myones, dp)*exp(-abs(d)/d0) + b0*(1 - sign(myones, dp)*exp(-abs(d)/d0))

        ! Housekeeping
        !=============
        end associate


    end subroutine

    ! Derivative evaluation
    subroutine EvaluateDerivativeStructuredPLF2DDistanceDF(distribution, &
        x, y, derivx, derivy, v)

        ! Description
        !============
        ! Evaluate the distribution function

        ! Declare variables
        !==================
        ! Arguments
        class(StructuredPLF2DDistanceDFUDT) :: distribution 
        real(R8), intent(in)                :: x(:), y(:)
        real(R8), intent(out)               :: v(size(x))
        integer(I8), intent(in)             :: derivx, derivy

        call gdErrorHandler('EvaluateDerivativeStructuredPLF2DDistanceDF: ' // & 
            'derivatives not yet implemented for this distribution type')

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
            ! Associate
            associate(  ne  => ps%polygons(i)%ne, &
                        nx  => ps%polygons(i)%nx, &
                        ny  => ps%polygons(i)%ny, &
                        nn  => ps%polygons(i)%nn, &
                        nv  => ps%polygons(i)%nv, & 
                        edges    => ps%polygons(i)%edges)
            ! Get point indices
            valindex = [(k, k = cc+1, cc+nv)]

            ! Set points
            xa(valindex) = ps%polygons(i)%x
            ya(valindex) = ps%polygons(i)%y

            ! Determine normal in points
            allocate(tnxa(nv), tnya(nv))
            tnxa = 0
            tnya = 0
            do j = 1, 2
                do k = 1, ne 
                    tnxa(edges(k, j)) = tnxa(edges(k, j)) + 0.5*nx(k)/nn(k)
                    tnya(edges(k, j)) = tnya(edges(k, j)) + 0.5*ny(k)/nn(k)
                end do 
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

            ! End asoociation
            end associate
        end do

        ! Construct rhs to compute attractor coefficients
        b = a0 - b0

        ! Compute lhs to compute attractor coefficients
        A = 0
        do j = 1, na
            do i = 1, na
                if (i /= j) then 
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
        distribution%dp = -(-dFdy*nxa + dFdx*nya)

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


            !$omp parallel do default(none) schedule(static) if (.not. omp_in_parallel()) &
            !$omp shared(distribution, x, y) &
            !$omp private(i, d) &
            !$omp reduction(+:v)
            do i = 1, size(xa)
                ! Distance 
                d = sqrt( (x - xa(i))**2 + (y - ya(i))**2)

                ! Value
                v = v + b0(i)*exp(-d/d0(i))
            end do
            !$omp end parallel do

        case ('signed')

            !$omp parallel do default(none) schedule(static) if (.not. omp_in_parallel()) &
            !$omp shared(distribution, x, y) &
            !$omp private(i, d) &
            !$omp reduction(+:v)
            do i = 1, size(xa)
                ! Distance 
                d = sqrt( (x - xa(i))**2 + (y - ya(i))**2)

                ! Value
                v = v + b0(i)*(sign(myone, dp(i)))*exp(-d/d0(i))
            end do
            !$omp end parallel do 

        case default

            call gdErrorHandler('Unknown method')

        end select

        v = a0 - v

        ! Housekeeping
        !=============
        end associate


    end subroutine
    
    ! Derivative evaluation
    subroutine EvaluateDerivativePolygonset2DFieldDistanceDF(distribution, &
        x, y, derivx, derivy, v)

        ! Description
        !============
        ! Evaluate the distribution function

        ! Declare variables
        !==================
        ! Arguments
        class(Polygonset2DFieldDistanceDFUDT)   :: distribution 
        real(R8), intent(in)                    :: x(:), y(:)
        real(R8), intent(out)                   :: v(size(x))
        integer(I8), intent(in)                 :: derivx, derivy

        call gdErrorHandler('EvaluateDerivativePolygonset2DFieldDistanceDF: ' // & 
            'derivatives not yet implemented for this distribution type')


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
                if (i /= j) then 
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

            !$omp parallel do default(none) schedule(static) if (.not. omp_in_parallel()) &
            !$omp shared(distribution, fv) &
            !$omp private(i, d) &
            !$omp reduction(+:v)
            do i = 1, size(xa)
                ! Distance 
                d = (fv - fval(i))

                ! Value
                v = v + b0*exp(-abs(d)/d0)
            end do
            !$omp end parallel do 

        case ('signed')

            !$omp parallel do default(none) schedule(static) if (.not. omp_in_parallel()) &
            !$omp shared(distribution, fv) &
            !$omp private(i, d) &
            !$omp reduction(+:v)
            do i = 1, size(xa)
                ! Distance 
                d = (fv - fval(i))

                ! Value
                v = v + b0*(sign(myone, d))*exp(-abs(d)/d0)
            end do
            !$omp end parallel do

        case default

            call gdErrorHandler('Unknown method')

        end select

        ! Housekeeping
        !=============
        end associate


    end subroutine

    ! Derivative evaluation
    subroutine EvaluateDerivativePolygonset1DFieldDistanceDF(distribution, &
        x, y, derivx, derivy, v)

        ! Description
        !============
        ! Evaluate the distribution function

        ! Declare variables
        !==================
        ! Arguments
        class(Polygonset1DFieldDistanceDFUDT)   :: distribution 
        real(R8), intent(in)                    :: x(:), y(:)
        real(R8), intent(out)                   :: v(size(x))
        integer(I8), intent(in)                 :: derivx, derivy

        call gdErrorHandler('EvaluateDerivativePolygonset1DFieldDistanceDF: ' // & 
            'derivatives not yet implemented for this distribution type')


    end subroutine 
    
    !------------------------------------------------------------------!
    !                      COORDINATES & FIELD, 1D                     !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCoordinates1DFieldDistanceDF(distribution, interp, &
        xp, yp, val0, valinf, decaylength, meth)

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
        class(Coordinates1DFieldDistanceDFUDT)           :: distribution 
        type(StructuredInterpolant2DUDT), intent(in)    :: interp 
        real(R8), intent(in)                            :: val0, valinf, &
            decaylength, xp(:), yp(:)
        character(*), intent(in)                        :: meth

        ! Auxiliary
        integer(I8)                                     :: flag, na

        real(R8), parameter                             :: myone = 1
        real(R8)                                        :: tempd
        real(R8), allocatable                           :: d(:), b(:), &
            bval(:), A(:, :), sol(:), fval(:)

        ! Loop
        integer(I8)                                     :: i, j

        ! Set fields
        !===========
        ! Data
        distribution%F  = interp 
        distribution%xa = xp 
        distribution%ya = yp
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
        na = size(xp, 1)

        ! Allocate
        allocate(d(na), bval(na), b(na), A(na, na))

        ! Field values
        allocate(fval(na))
        call distribution%F%Evaluate(xp, yp, 0, 0, fval)

        distribution%fval = fval 

        ! Construct rhs to compute attractor coefficients
        b = a0 - b0

        ! Compute lhs to compute attractor coefficients
        A = 0
        do j = 1, na
            do i = 1, na
                if (i /= j) then 
                    tempd = abs(fval(i) - fval(j))
                    A(i, j) = sign(myone, fval(i)-fval(j))*exp(-tempd/d0)
                else 
                    A(i, j) = 1
                end if 
            end do 
        end do
        where (A <= 1e-10) A = 0

        ! Call solver
        allocate(sol(size(b)))
        call SolveDenseLinearSystemDI(A, b, sol, flag)
        if (flag /= 0) then
            ! Call error
            call gdErrorHandler('InitializeCoordinates1DDistanceDF: ' // &
                'could not determine attractor function coefficients ' // &
                'due to non-converging linear solver')
        end if 

        ! Add
        !====
        distribution%fval   = fval
        distribution%coef   = sol 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Evaluation
    subroutine EvaluateCoordinates1DFieldDistanceDF(distribution, x, y, v)

        ! Description
        !============
        ! Evaluate the distribution function

        ! Declare variables
        !==================
        ! Arguments
        class(Coordinates1DFieldDistanceDFUDT)   :: distribution 
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
            c       => distribution%coef,   & 
            b0      => distribution%b0,     & 
            d0      => distribution%d0,     &
            xa      => distribution%xa,     &
            ya      => distribution%ya,     &
            fval    => distribution%fval    &
        )

        ! Evaluate
        !=========
        ! Evaluate field values in coordinates
        v = 0
        call F%Evaluate(x, y, 0, 0, fv)
        select case (distribution%meth) 

        case ('unsigned')

            !$omp parallel do default(none) schedule(static) if (.not. omp_in_parallel()) &
            !$omp shared(distribution, fv) &
            !$omp private(i, d) &
            !$omp reduction(+:v)
            do i = 1, size(xa)
                ! Distance 
                d = (fv - fval(i))

                ! Value
                v = v + c(i)*exp(-abs(d)/d0)
            end do
            !$omp end parallel do 

        case ('signed')

            !$omp parallel do default(none) schedule(static) if (.not. omp_in_parallel()) &
            !$omp shared(distribution, fv) &
            !$omp private(i, d) &
            !$omp reduction(+:v)
            do i = 1, size(xa)
                ! Distance 
                d = (fv - fval(i))

                ! Value
                v = v + c(i)*(sign(myone, d))*exp(-abs(d)/d0)
            end do
            !$omp end parallel do 

        case default

            call gdErrorHandler('Unknown method')

        end select

        ! Add constant component
        v = v + b0

        ! Housekeeping
        !=============
        end associate


    end subroutine

    ! Derivative evaluation
    subroutine EvaluateDerivativeCoordinates1DFieldDistanceDF(distribution, &
        x, y, derivx, derivy, v)

        ! Description
        !============
        ! Evaluate the distribution function

        ! Declare variables
        !==================
        ! Arguments
        class(Coordinates1DFieldDistanceDFUDT)  :: distribution 
        real(R8), intent(in)                    :: x(:), y(:)
        real(R8), intent(out)                   :: v(size(x))
        integer(I8), intent(in)                 :: derivx, derivy

        call gdErrorHandler('EvaluateDerivativeCoordinates1DFieldDistanceDF: ' // & 
            'derivatives not yet implemented for this distribution type')


    end subroutine 

    !------------------------------------------------------------------!
    !                          COORDINATES, 2D                         !
    !------------------------------------------------------------------!

    ! Constructor
    function ConstructCoordinates2DDistanceDF(xp, yp, val0, valinf, &
        decaylengthp) result(distribution)

        ! Description
        !============
        ! Construct the distributor - simply a wrapper for the 
        ! initialization function

        ! Declare variables
        !==================
        ! Arguments
        class(DistributionFunctionUDT), allocatable     :: distribution 
        real(R8), intent(in)                            :: val0(:), valinf, &
            xp(:), yp(:), decaylengthp(:)

        ! Initialize
        !===========
        allocate(Coordinates2DDistanceDFUDT::distribution)

        select type(distribution)

        type is (Coordinates2DDistanceDFUDT)

            ! Call initializer
            call distribution%Initialize(xp, yp, val0, valinf, decaylengthp)

        end select

    end function 

    ! Initialization
    subroutine InitializeCoordinates2DDistanceDF(distribution, &
        xp, yp, val0, valinf, decaylength)

        ! Description
        !============
        ! Initialization routine. The 'interp' structured interpolant
        ! must be initialized and correctly set up. The argument 'val0'
        ! is the value that the function achieves on the points xp, yp,
        ! the value 'valinf' is achieved at locations very far of these 
        ! points. 'decaylength' is a decay length of how fast val0 transitions 
        ! to valinf. 

        ! Declare variables
        !==================
        ! Arguments
        class(Coordinates2DDistanceDFUDT)               :: distribution 
        real(R8), intent(in)                            :: val0(:), valinf, &
            decaylength(:), xp(:), yp(:)

        ! Auxiliary
        integer(I8)                                     :: flag, na

        real(R8), parameter                             :: myone = 1
        real(R8)                                        :: tempd
        real(R8), allocatable                           :: d(:), b(:), &
            A(:, :), sol(:)

        logical, allocatable, dimension(:)              :: isduplicate

        ! Loop
        integer(I8)                                     :: i, j

        ! Set fields
        !===========
        ! Data
        distribution%xa = xp 
        distribution%ya = yp
        distribution%a0 = val0
        distribution%b0 = valinf 
        distribution%d0 = decaylength

        ! Associate
        associate(&
            a0      => distribution%a0,     &
            b0      => distribution%b0,     &
            d0      => distribution%d0)

        ! Construct attractor function
        !=============================
        ! Determine number of attractor points
        na = size(xp, 1)

        ! Allocate
        allocate(d(na), b(na), A(na, na), isduplicate(na))

        ! Construct rhs to compute attractor coefficients
        b = a0 - b0

        ! Compute lhs to compute attractor coefficients
        isduplicate = .false. 
        A = 0
        !$omp parallel do default(none) schedule(static) collapse(2) &
        !$omp shared(na, xp, yp, isduplicate, A, distribution) &
        !$omp private(tempd) if(.not. omp_in_parallel())
        do j = 1, na
            do i = 1, na
                if (i /= j) then 
                    tempd = sqrt( (xp(i) - xp(j))**2 + (yp(i) - yp(j))**2)
                    if (tempd == 0) then 
                        !$omp critical
                        isduplicate(j) = .true.
                        !$omp end critical
                    end if  
                    A(i, j) = exp(-tempd/d0(j))
                else 
                    A(i, j) = 1.0_R8
                end if 
            end do 
        end do
        !$omp end parallel do

        ! Adjust
        do j = 1, na
            if (isduplicate(j)) then 
                A(j, :) = 0
                A(:, j) = 0
                A(j, j) = 1
                b(j) = 1
            end if 
        end do

        ! Call solver
        allocate(sol(size(b)))
        call SolveDenseLinearSystemDI(A, b, sol, flag)
        if (flag /= 0) then
            ! Call error
            call gdErrorHandler('InitializeCoordinates2DDistanceDF: ' // &
                'could not determine attractor function coefficients ' // &
                'due to non-converging linear solver')
        end if 

        ! Add
        !====
        distribution%coef   = sol 

        ! Housekeeping
        !=============
        end associate

        ! Reset distribution quantities to exclude duplicate points
        distribution%coef = pack(distribution%coef, .not. isduplicate)
        distribution%xa = pack(distribution%xa, .not. isduplicate) 
        distribution%ya = pack(distribution%ya, .not. isduplicate)
        distribution%a0 = pack(distribution%a0, .not. isduplicate)
        distribution%d0 = pack(distribution%d0, .not. isduplicate)

    end subroutine

    ! Evaluation
    subroutine EvaluateCoordinates2DDistanceDF(distribution, x, y, v)

        ! Description
        !============
        ! Evaluate the distribution function

        ! Declare variables
        !==================
        ! Arguments
        class(Coordinates2DDistanceDFUDT)       :: distribution 
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
            c       => distribution%coef,   & 
            b0      => distribution%b0,     & 
            d0      => distribution%d0,     &
            xa      => distribution%xa,     &
            ya      => distribution%ya      &
        )

        ! Evaluate
        !=========
        ! Evaluate field values in coordinates
        v = 0

        !$omp parallel do default(none) schedule(static) if (.not. omp_in_parallel()) &
        !$omp shared(distribution, x, y) &
        !$omp private(i, d) &
        !$omp reduction(+:v)
        do i = 1, size(xa)
            ! Distance 
            d = sqrt((x - xa(i))**2 + (y - ya(i))**2)

            ! Value
            v = v + c(i)*exp(-d/d0(i))
        end do
        !$omp end parallel do 

        ! Add constant component
        v = v + b0

        ! Housekeeping
        !=============
        end associate


    end subroutine

    ! Derivative evaluation
    subroutine EvaluateDerivativeCoordinates2DDistanceDF(distribution, &
        x, y, derivx, derivy, v)

        ! Description
        !============
        ! Evaluate the distribution function

        ! Declare variables
        !==================
        ! Arguments
        class(Coordinates2DDistanceDFUDT)       :: distribution 
        real(R8), intent(in)                    :: x(:), y(:)
        real(R8), intent(out)                   :: v(size(x))
        integer(I8), intent(in)                 :: derivx, derivy 

        ! Auxiliary
        real(R8), parameter                     :: myone = 1
        real(R8)                                :: d(size(x))
        real(R8), allocatable, dimension(:)     :: dddx, dddy 

        ! Loop
        integer(I8)                             :: i

        ! Initialize
        !===========
        ! Check sizes
        if ( (size(v) /= size(x)) .or. (size(x) /= size(y))) then 
            ! Throw error
            call gdErrorHandler('EvaluateDerivativeCoordinates2DDistanceDF: ' // & 
                'incompatible sizes in input')
        end if 

        ! Associate
        associate(&
            a0      => distribution%a0,     & 
            c       => distribution%coef,   & 
            b0      => distribution%b0,     & 
            d0      => distribution%d0,     &
            xa      => distribution%xa,     &
            ya      => distribution%ya      &
        )

        ! Evaluate
        !=========
        ! Evaluate field values in coordinates
        v = 0

        select case (derivx)

        case (0)

            select case (derivy)

            case (0)

                ! Field value
                do i = 1, size(xa)
                    ! Distance 
                    d = sqrt((x - xa(i))**2 + (y - ya(i))**2)

                    ! Value
                    v = v + c(i)*exp(-d/d0(i))
                end do

                ! Add constant component
                v = v + b0

            case (1)

                ! dfdy 
                do i = 1, size(xa)
                    ! Distance 
                    d = sqrt((x - xa(i))**2 + (y - ya(i))**2)
                    dddy = (y - ya(i))/d

                    ! Value
                    v = v - c(i)/d0(i)*exp(-d/d0(i))*dddy
                end do

            case (2)

                ! d2fdy2
                do i = 1, size(xa)
                    ! Distance 
                    d = sqrt((x - xa(i))**2 + (y - ya(i))**2)
                    v = v + c(i)*(-1/d + (y - ya(i))**2/(d**2*d0(i)) + (y - ya(i))**2/d**3)*exp(-d/d0(i))/d0(i)
                end do

            case default

                call gdErrorHandler('EvaluateDerivativeCoordinates2DDistanceDF: ' // & 
                    'derivative not implemented')

            end select

        case (1)

            select case (derivy)

            case (0)

                ! dfdx
                do i = 1, size(xa)
                    ! Distance 
                    d = sqrt((x - xa(i))**2 + (y - ya(i))**2)
                    dddx = (x - xa(i))/d

                    ! Value
                    v = v - c(i)/d0(i)*exp(-d/d0(i))*dddx
                end do

            case (1)

                ! d2fdxdy
                do i = 1, size(xa)
                    ! Distance 
                    d = sqrt((x - xa(i))**2 + (y - ya(i))**2)

                    ! Value
                    v = v + c(i)*(x - xa(i))*(y - ya(i))*(1/(d**2*d0(i)) + d**(-3))*exp(-d/d0(i))/d0(i)
                end do 

            case default

                call gdErrorHandler('EvaluateDerivativeCoordinates2DDistanceDF: ' // & 
                    'derivative not implemented')

            end select

        case (2)

            select case(derivy)

            case (0)

                ! d2fdx2
                do i = 1, size(xa)
                    ! Distance 
                    d = sqrt((x - xa(i))**2 + (y - ya(i))**2)

                    ! Value
                    v = v + c(i)*(-1/d + (x - xa(i))**2/(d**2*d0(i)) + (x - xa(i))**2/d**3)*exp(-d/d0(i))/d0(i)
                end do 

            case default 

                call gdErrorHandler('EvaluateDerivativeCoordinates2DDistanceDF: ' // & 
                    'derivative not implemented')

            end select

        case default 

            call gdErrorHandler('EvaluateDerivativeCoordinates2DDistanceDF: ' // & 
                    'derivative not implemented')

        end select

        ! Housekeeping
        !=============
        end associate


    end subroutine

    !------------------------------------------------------------------!
    !                      COORDINATES AND PLF, 2D                     !
    !------------------------------------------------------------------!

    ! Constructor
    function ConstructCoordinatesPLF2DDistanceDF(plf, valplf, decaylengthplf, &
        xp, yp, val0, valinf, decaylengthp) result(distribution)

        ! Description
        !============
        ! Construct the distributor - simply a wrapper for the 
        ! initialization function

        ! Declare variables
        !==================
        ! Arguments
        class(DistributionFunctionUDT), allocatable     :: distribution 
        real(R8), intent(in)                            :: val0(:), valinf, &
            decaylengthp(:), xp(:), yp(:), valplf, decaylengthplf
        class(PolygonLevelsetFunction2DUDT), intent(in) :: plf 

        ! Initialize
        !===========
        allocate(CoordinatesPLF2DDistanceDFUDT::distribution)

        select type(distribution)

        type is (CoordinatesPLF2DDistanceDFUDT)

            ! Call initializer
            call distribution%Initialize(plf, valplf, decaylengthplf, &
                xp, yp, val0, valinf, decaylengthp)

        end select

    end function 

    ! Initialization
    subroutine InitializeCoordinatesPLF2DDistanceDF(distribution, &
        plf, valplf, decaylengthplf, xp, yp, val0, valinf, decaylengthp)


        ! Declare variables
        !==================
        ! Arguments
        class(CoordinatesPLF2DDistanceDFUDT)            :: distribution 
        real(R8), intent(in)                            :: val0(:), valinf, &
            decaylengthp(:), xp(:), yp(:), valplf, decaylengthplf
        class(PolygonLevelsetFunction2DUDT), intent(in) :: plf 

        ! Auxiliary
        integer(I8)                                     :: flag, na

        real(R8), parameter                             :: myone = 1
        real(R8)                                        :: tempd
        real(R8), allocatable                           :: d(:), b(:), &
            A(:, :), sol(:), vplf(:)

        ! Loop
        integer(I8)                                     :: i, j

        ! Set fields
        !===========
        ! Data
        distribution%xa = xp 
        distribution%ya = yp
        distribution%a0 = val0
        distribution%b0 = valinf 
        distribution%d0 = decaylengthp

        distribution%plf = plf 
        distribution%d_plf = decaylengthplf
        distribution%a_plf = valplf

        ! Associate
        associate(&
            a0      => distribution%a0,     &
            b0      => distribution%b0,     &
            d0      => distribution%d0,     &
            a_plf   => distribution%a_plf,  &
            d_plf   => distribution%d_plf,  &
            plf     => distribution%plf     &
            )

        ! Construct attractor function
        !=============================
        ! Determine number of attractor points
        na = size(xp, 1)

        ! Allocate
        allocate(d(na), b(na), A(na, na))

        ! Evaluate plf at point locations (to subtract later)
        allocate(vplf(na))
        call plf%Evaluate(xp, yp, 0, 0, vplf)

        ! Construct rhs to compute attractor coefficients
        b = a0 - b0 - a_plf*exp(-abs(vplf)/d_plf)

        ! Compute lhs to compute attractor coefficients
        A = 0
        do j = 1, na
            do i = 1, na
                if (i /= j) then 
                    tempd = sqrt( (xp(i) - xp(j))**2 + (yp(i) - yp(j))**2)
                    A(i, j) = exp(-tempd/d0(j))
                else 
                    A(i, j) = 1.0_R8
                end if 
            end do 
        end do

        ! Call solver
        allocate(sol(size(b)))
        call SolveDenseLinearSystemDI(A, b, sol, flag)
        if (flag /= 0) then
            ! Call error
            call gdErrorHandler('InitializeCoordinatesPLF2DDistanceDF: ' // &
                'could not determine attractor function coefficients ' // &
                'due to non-converging linear solver')
        end if 

        ! Add
        !====
        distribution%coef   = sol 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Evaluation
    subroutine EvaluateCoordinatesPLF2DDistanceDF(distribution, x, y, v)

        ! Description
        !============
        ! Evaluate the distribution function

        ! Declare variables
        !==================
        ! Arguments
        class(CoordinatesPLF2DDistanceDFUDT)    :: distribution 
        real(R8), intent(in)                    :: x(:), y(:)
        real(R8), intent(out)                   :: v(size(x))

        ! Auxiliary
        real(R8), parameter                     :: myone = 1
        real(R8)                                :: d(size(x))
        real(R8), allocatable, dimension(:)     :: vplf

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
            a_plf   => distribution%a_plf,  &
            d_plf   => distribution%d_plf,  &
            plf     => distribution%plf,    &
            c       => distribution%coef,   & 
            b0      => distribution%b0,     & 
            d0      => distribution%d0,     &
            xa      => distribution%xa,     &
            ya      => distribution%ya      &
        )

        ! Evaluate
        !=========
        ! Evaluate field values in coordinates
        v = 0

        ! Point contributions
        !$omp parallel do default(none) schedule(static) if (.not. omp_in_parallel()) &
        !$omp shared(distribution, x, y) &
        !$omp private(i, d) &
        !$omp reduction(+:v)
        do i = 1, size(xa)
            ! Distance 
            d = sqrt((x - xa(i))**2 + (y - ya(i))**2)

            ! Value
            v = v + c(i)*exp(-d/d0(i))
        end do
        !$omp end parallel do 

        ! PLF contributions
        allocate(vplf(size(x)))
        call plf%Evaluate(x, y, 0, 0, vplf)
        v = v + a_plf*exp(-abs(vplf)/d_plf)

        ! Add constant component
        v = v + b0

        ! Housekeeping
        !=============
        end associate


    end subroutine

    ! Derivative evaluation
    subroutine EvaluateDerivativeCoordinatesPLF2DDistanceDF(distribution, &
        x, y, derivx, derivy, v)

        ! Description
        !============
        ! Evaluate the distribution function

        ! Declare variables
        !==================
        ! Arguments
        class(CoordinatesPLF2DDistanceDFUDT)    :: distribution 
        real(R8), intent(in)                    :: x(:), y(:)
        real(R8), intent(out)                   :: v(size(x))
        integer(I8), intent(in)                 :: derivx, derivy 

        call gdErrorHandler('EvaluateDerivativeCoordinatesPLF2DDistanceDF: ' // & 
            'derivatives not yet implemented for this distribution type')
        

    end subroutine



    

end module