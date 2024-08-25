!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module provides contouring functionality for 2D fields. Currently
! only data provided on structured 2D Cartesian meshes is supported, 
! though this may be extended in the future to triangular meshes as
! well (these allow much easier local refinement). The contours are 
! typically returned as a contour UDT (though one may implement wrappers
! to return this in more basic formats for interfacing). 

module mod_contour2D

    ! Load modules
    use mod_precision
    use mod_errorhandler
    use mod_dynamicarrays
    use mod_structured2Dgridding

    implicit none
    private 
    public :: ContourUDT, TraceContoursStructured2D, ContourTracerUDT, &
        StructuredContourTracerUDT, ConstructStructuredTracer

    ! Module parameters
    integer, parameter  :: npq          = 2 ! number of padding quads for 2D tracer to determine saddle points
    integer, parameter  :: verbosity    = 0 ! verbosity level (0: only true errors, 1: additional messages)
    real(R8), parameter :: pert         = 1e-13 ! perturbation value 
    real(R8), parameter :: spvalabstol  = 1e-13 ! absolute tolerance in field value to determine if value is equal to saddlepoint value
    real(R8), parameter :: spvalreltol  = 1e-10 ! relative tolerance for ^ 

    !==================================================================!
    !                                                                  !
    !                            TYPES                                 !
    !                                                                  !
    !==================================================================!

    ! Contour type
    type :: ContourUDT 

        ! Description
        !============
        ! Each contour represents a single, simple line on which the 
        ! field value is constant. The contour always forms a simple
        ! polygon, but may be closed upon itself (end point is 
        ! repeated). If the contour is non-simple (e.g. a saddle point
        ! contour), then multiple contour parts will be returned that 
        ! can be stitched back together using the ID information (the 
        ! contour ID indicates from which original seed value the 
        ! contour was computed - so same ID, part of the same total 
        ! contour). Because saddle points play an important role in 
        ! many applications and may be known beforehand, they may 
        ! typically be passed to the routine and accounted for. In the
        ! contour type, we also store the index of the saddle point if 
        ! encountered at the start and/or end. 

        ! The contour type contains the following data:
        ! - x, y   coordinates of the contour (as simple arrays)
        ! - val    field value of the contour
        ! - ID     ID that maps back to the initial trace value (to 
        !          merge contour parts together afterwards if the 
        !          contour is not a simple line)
        ! - isclosed        logical indicating if contour part is closed
        ! - startsaddle     index of the saddle point at which the 
        !                   contour starts (zero if start is no saddle)
        ! - endsaddle       similar as above, but for end point   
    
        real(R8), allocatable   :: x(:), y(:)
        real(R8)                :: val
        integer(I8)             :: ID, startsaddle, endsaddle
        logical                 :: isclosed 

    contains 

    end type

    ! General contour tracer 
    type, abstract :: ContourTracerUDT

        ! Description
        !============
        ! General type that traces contours. Allows abstraction from any
        ! underlying mechanics (structured/unstructured grid etc), except
        ! on the constructor level. All necessary data to trace the 
        ! contours should be saved in the object itself, except for the 
        ! saddle point data - that is already available in this 
        ! type already, since it should be supported in any tracer here.
        ! For easier redefinition of tracing values, routines should be
        ! specified that return the location on which the field value
        ! should be known, and a routine that sets the values again. 
        real(R8), allocatable, dimension(:)     :: xs, ys, vs
        integer(I8), allocatable, dimension(:)  :: order 

    contains 

        ! Tracer
        procedure(TraceContoursValINT), deferred    :: TraceContoursVal 
        procedure(TraceContoursLocINT), deferred    :: TraceContoursLoc
        generic :: TraceContours        => TraceContoursVal, &
            TraceContoursLoc

        ! Coordinate getter
        procedure(GetTracerValueCoordinatesINT), deferred :: GetCoordinates 

        ! Field value setter
        procedure(SetTracerValuesINT), deferred     :: SetValues

    end type

    ! Structured 2D tracer
    type, extends(ContourTracerUDT)     :: StructuredContourTracerUDT

        ! Description
        !============
        ! Tracer for structured 2D grids. Values etc are all saved in 
        ! the structure, the tracer routine is a wrapper for the 
        ! standalone 2D tracer routine. 
        real(R8), allocatable       :: X(:), Y(:), V(:, :), xg(:), yg(:)


    contains 

        ! Tracer
        procedure :: TraceContoursVal  => TraceContoursValStructured2DWrapper
        procedure :: TraceContoursLoc  => TraceContoursLocStructured2DWrapper

        ! Coordinate getter
        procedure :: GetCoordinates => GetCoordinatesStructured2D

        ! Field value setter
        procedure :: SetValues      => SetValuesStructured2D

    end type 

    ! Saddle point structure type (only used in this module)
    type :: sp2DUDT 

        ! Description
        !============
        ! Structure that contains additional saddle point information
        ! used during the contour line tracing to traverse saddle 
        ! points in a correct way. 

        integer(I8)                 :: ID, ixquad, iyquad, order
        integer(I8), allocatable    :: tri(:, :), ixquadtri(:), &
            iyquadtri(:), ixpoints(:), iypoints(:)
        real(R8)                    :: x, y, val 
        real(R8), allocatable       :: valpoints(:)

    contains 

    end type 

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    ! Abstract interfaces
    !====================
    ! Tracer interface
    abstract interface 

        ! Tracer routine
        function TraceContoursValINT(tracer, tracevalues) result(contours)
            
            import :: ContourTracerUDT, R8, ContourUDT 
            class(ContourTracerUDT)         :: tracer 
            real(R8), intent(in)            :: tracevalues(:) 
            type(ContourUDT), allocatable   :: contours(:)

        end function

        function TraceContoursLocINT(tracer, x, y) result(contours)
            
            import :: ContourTracerUDT, R8, ContourUDT 
            class(ContourTracerUDT)         :: tracer 
            real(R8), intent(in)            :: x(:), y(:)
            type(ContourUDT), allocatable   :: contours(:)

        end function

        ! Coordinates getter
        subroutine GetTracerValueCoordinatesINT(tracer, x, y) 
            
            import :: ContourTracerUDT, R8
            class(ContourTracerUDT)             :: tracer 
            real(R8), allocatable, intent(out)  :: x(:), y(:)

        end subroutine

        ! Values setter
        subroutine SetTracerValuesINT(tracer, v) 
            
            import :: ContourTracerUDT, R8
            class(ContourTracerUDT)             :: tracer 
            real(R8), intent(in)                :: v(:)

        end subroutine

    end interface

    ! Normal interfaces
    !==================
    ! Contour addition
    interface AddContours 
        module procedure AddContourArray
        module procedure AddContourScalar
    end interface

    contains 

    !==================================================================!
    !                                                                  !
    !                             ROUTINES                             !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                           CONSTRUCTORS                           !
    !------------------------------------------------------------------!

    ! Contour constructor
    function ConstructContour(x, y, val, ID, ss, es, isclosed) result(contour)

        ! Description
        !============
        ! Contour constructor that requires all contour data to be set
        ! to be passed. Only for single contours

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)        :: x(:), y(:), val 
        integer(I8), intent(in)     :: ID, ss, es
        logical, intent(in)         :: isclosed 
        type(ContourUDT)            :: contour

        ! Construct
        !==========
        contour%x           = x 
        contour%y           = y 
        contour%val         = val 
        contour%ID          = ID 
        contour%startsaddle = ss 
        contour%endsaddle   = es 
        contour%isclosed    = isclosed

    end function

    ! Structured tracer constructor
    function ConstructStructuredTracer(V, X, Y, xs, ys, vs) result(tracer)

        ! Description
        !============
        ! Construct a structured tracer. The result is returned as an 
        ! abstract type though that is allocated in this routine

        ! Declare variables
        !==================
        ! Arguments
        class(ContourTracerUDT), allocatable    :: tracer 
        real(R8), intent(in), dimension(:)      :: X, Y, xs, ys, vs 
        real(R8), intent(in)                    :: V(:, :)

        ! Auxiliary
        integer(I8), allocatable                :: order(:)
        real(R8), allocatable                   :: xg(:), yg(:)


        ! Initialize
        !===========
        ! Allocate
        allocate(StructuredContourTracerUDT::tracer) 
        allocate(order(size(xs)))

        ! Set values
        !===========
        select type (tracer)

        type is (StructuredContourTracerUDT)

            tracer%V = V 
            tracer%X = X 
            tracer%Y = Y 
            tracer%xs = xs 
            tracer%ys = ys 
            tracer%vs = vs 
            tracer%order = order  

            allocate(xg(size(X)*size(Y)), yg(size(X)*size(Y)))
            call Construct2DStructuredGrid(X, Y, size(X), size(Y), xg, yg)
            tracer%xg = xg 
            tracer%yg = yg

        class default 

        end select

    end function 

    !------------------------------------------------------------------!
    !                              TRACERS                             !
    !------------------------------------------------------------------!

    ! 2D structured contour line tracer 
    subroutine TraceContoursStructured2D(V, X, Y, tracevalues, xs, ys, &
        vs, order, contours) 

        ! Description
        !============
        ! Trace the contour lines of a field that is given by its values in 'V' on
        ! a structured (not necessarily uniform) 2D grid, determined by the grid
        ! vectors X and Y. It is assumed that value V(i, j) corresponds to the
        ! value at X(i), Y(j). 'Tracevalues' is a n-by-1 array, containing the to
        ! be tracedfield values.
        ! It is assumed that the field varies linearly between the grid points.

        ! The output consists of a structure with the fields 'V' (trace value,
        ! scalar), 'x', y', (array of coordinates), 'ID', (index of traced value or
        ! point, to be able to merge different parts together afterwards)
        ! Each line represents a simple
        ! polygon that is either closed upon itself (end point is repeated) or
        ! starts and stops at a mesh boundary. If saddle points are present, the
        ! different parts of the contour line are treated as different polygon
        ! pieces, where each end point is either a grid boundary or a saddle point.
        ! The location of the saddle point itself is determined as the center of
        ! the cell where the saddle point is present.

        ! If it is known that there are saddle points present in the field, that
        ! may not be captured in a single quad (this is in fact rather likely),
        ! they can be parsed as optional argument 'saddlepoints'. For each saddle
        ! point, we identify the quad where it lies in and take an additional
        ! layer of quads (actually multiple layers could be possible, it is set as a
        ! hard coded parameter now).

        ! Algorithm
        !==========
        ! When searching contours for specified values F, the following algorithm is
        ! employed:
        !
        ! 1) Check if V >= F
        ! 2) Check each quadrilateral:
        !       if it has the value on all four edges, it is a saddle point
        !       otherwise it is a regular quad
        ! 3) First, trace all closed field lines starting from the boundary of the
        ! domain:
        !   3.1) Determine the next boundary face to start, if none left, go to 4)
        !   A face should always have a counter equal to one or zero, and can only
        !   be considered once per contour (therefore when subtracting/putting
        !   counter to zero, face cannot be considered anymore)
        !   3.2) Determine the current quad and subtract from quad counter. 
        !   If this quad is a saddle point, exit and go to 3.1)
        !   3.3) Find the next face. If it is a boundary face -> subtract, exit,
        !   and go to 3.1. Otherwise, subtract and go to 3.2. 
        ! 4) Normally, only open contour lines remain now:
        !   4.1) Find the next quad: either start with saddle point quad with
        !   counter > 0 or next quad with counter > 0. 
        !   4.2) Find the next face that has not yet been considered and that
        !   contains the value (counter > 0). If no faces are found, we should have
        !   arrived at the original cell (is checked, otherwise error -> would indicate bug in code)
        !   4.3) Determine the quad indices and flag. If it is a saddle point,
        !   subtract, exit, and go to 4.1. Otherwise, go to 4.2. 
        !   

        ! Notes
        !======
        ! Note 1: it is assumed that there is some variation in V over X, Y.
        ! Otherwise, any contour algorithm will return bullshit (e.g when (part of)
        ! the field is constant).

        ! Note 2: if values lie exactly on a node (or multiple nodes), the value in
        ! these nodes is elevated by 1e-13 (a message will be displayed if this
        ! happens).

        ! Note 3: actually, the way of tackling saddle points here can be extended
        ! to an arbitrary high order of saddle points by increasing the number of
        ! quads padded around the quad that contains the saddle point (at the cost
        ! of local resolution of course, but given that fields are typically poorly
        ! discretized anyway close to that point, this is a small price to pay).
        ! Furthermore, the order can be identified by simply counting the amount of
        ! faces that contain the x-point (if n faces contain the value, then the
        ! x-point is of order n/2-1 -> 4 faces, first order, 6 faces, second order
        ! etc etc). If the saddle points are given, this is returned as an extra
        ! output. 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)            :: V(:, :), X(:), Y(:), &
            tracevalues(:)
        real(R8), intent(in)            :: xs(:), ys(:), vs(:)
        type(ContourUDT), allocatable, intent(out)  :: contours(:) 
        integer(I8), intent(inout)      :: order(:) 

        ! Auxiliary
        type(sp2DUDt), allocatable      :: spstruct(:)
        type(ContourUDT), allocatable   :: tempcontours(:)
        logical, allocatable            :: superquadfacexflags(:, :), &
            superquadfaceyflags(:, :)
        integer(I8)                     :: nx, ny 
        integer(I8), allocatable        :: superquadflags(:, :)
        real(R8)                        :: tv
        real(R8), allocatable           :: emptyarrayR8(:)

        ! Loop
        integer(I8)                     :: i, j

        ! Checks
        !=======
        ! Check inputs
        nx = size(X)
        ny = size(Y)
        if ((nx /= size(V, 1)) .or. (ny /= size(V, 2))) then 
            call gdErrorHandler('TraceContourStructured2D: incompatible ' // &
                'dimensions of arguments V, X, Y')
        end if 
        if ((size(order) /= size(xs)) .or. (size(xs) /= size(ys)) .or. &
            (size(ys) /= size(vs))) then 
            call gdErrorHandler('TraceContourStructured2D: incompatible ' // &
            'dimensions of arguments xs, ys, vs, order')
        end if 
        if (size(tracevalues) == 0) then 
            ! Empty trace value array - not allowed
            call gdErrorHandler('TraceContourStructured2D: trace value ' // & 
                'array appears to be empty, not supported')
        end if 

        ! Initialize
        !===========
        ! Empty array
        allocate(emptyarrayR8(0))

        ! Contours
        allocate(contours(0))

        ! Initialize the saddle point structure
        allocate(superquadflags(nx-1, ny-1), superquadfacexflags(nx, ny-1), &
            superquadfaceyflags(nx-1, ny))
        call InitializeSaddlePointStructure2D(spstruct, superquadflags, &
            superquadfacexflags, superquadfaceyflags, V, X, Y, xs, ys, &
            vs)
        do i = 1, size(spstruct)
            order(i) = spstruct(i)%order
        end do 

        ! Main loop
        !==========
        do i = 1, size(tracevalues)
            ! Get current trace value
            tv = tracevalues(i)

            ! Call tracer
            tempcontours = TraceSingleContourStructured2D(V, X, Y, tv, &
                spstruct, superquadflags, superquadfacexflags, &
                superquadfaceyflags, nx, ny)

            ! Set ID
            do j = 1, size(tempcontours)
                tempcontours(j)%ID = i 
            end do

            ! Add contours
            contours = AddContours(contours, tempcontours)

        end do 



    end subroutine

    ! 2D structured contour line tracer that starts from point locations
    subroutine TracecontoursStructured2DPoint(V, X, Y, xt, yt, xs, ys, &
        vs, order, contours)

        ! Description
        !============
        ! Similar routine as TracecontourLineStructured2D, but now we start from a
        ! given set of points x0, y0 instead of trace values. First, we check in
        ! which quad the point is present (if no quad, an empty contour line is
        ! returned - the point lies outside of the domain in that case). Then, we
        ! use bilinear interpolation to determine the contour value to be traced.
        ! The tracing can only start in the current point itself and any saddle
        ! points that are encountered along the way.

        ! The output consists of a structure with the fields 'V' (trace value,
        ! scalar), 'x', y', (array of coordinates), 'ID', (index of traced
        ! point, to be able to merge different parts together afterwards)
        ! Each line represents a simple
        ! polygon that is either closed upon itself (end point is repeated) or
        ! starts and stops at a mesh boundary. If saddle points are present, the
        ! different parts of the contour line are treated as different polygon
        ! pieces, where each end point is either a grid boundary or a saddle point.
        ! The location of the saddle point itself is determined as the center of
        ! the cell where the saddle point is present.

        ! Algorithm
        !==========
        ! When searching contours for specified points, the following algorithm is
        ! employed:
        !
        ! 1) Check if V >= F
        ! 2) Check each quadrilateral:
        !       if it has the value on all four edges, it is a saddle point
        !       otherwise it is a regular quad
        ! 3) Determine from which quad to start and determine the trace value using
        ! bilinear interpolation. Then determine the quad types and which faces
        ! that contain the value like in any contour tracing algorithm. However,
        ! now, we keep separate track of which quads that can be started from (and
        ! how many times). For the initial quad, this value can be two or four.
        ! 4) Trace the contour line:
        !   4.1) Check if there are starting quads left, if yes, take one. IF no,
        !   exit (all contour parts should have been found)
        !   4.2) Determine the first face with the current value that has not yet
        !   been considered, determine the next quad. If no face is found, we
        !   should have ended up in the first quad and the contour line should be
        !   closed (we stop in saddle points before, see below)
        !   4.3) If the next quad is a saddle point, add this saddle point as a
        !   starting quad and stop. Go to 4.1)

        ! Notes
        !======
        ! Note 1: it is assumed that there is some variation in V over X, Y.
        ! Otherwise, any contour algorithm will return bullshit (e.g when (part of)
        ! the field is constant).

        ! Note 2: if values lie exactly on a node (or multiple nodes), the value in
        ! these nodes is elevated by 1e-13 (a message will be displayed if this
        ! happens).

        ! Note 3: the coordinates x0 and y0 will likely not be a part of the traced
        ! contour! It is simply used to determine the starting point of the contour
        ! and to determine which part of the contour one would like.

        ! Note 4: starting from the exact location of an x-point is possible.
        ! However, the interpolated value is likely not equal to the actual value
        ! of the x-point. In order to be able to trace the contour also beyond
        ! other x-point that did originally have exactly the same value (or have
        ! been set to do so), we check whether the starting point is an x-point,
        ! and set the value of all x-points that have exactly the same value as
        ! this x-point equal to the interpolated value. 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)            :: V(:, :), X(:), Y(:), &
            xt(:), yt(:)
        real(R8), intent(in)            :: xs(:), ys(:), vs(:)
        type(ContourUDT), allocatable, intent(out)  :: contours(:) 
        integer(I8), intent(inout)      :: order(:)

        ! Auxiliary
        type(sp2DUDt), allocatable      :: spstruct(:)
        type(ContourUDT), allocatable   :: tempcontours(:)
        logical, allocatable            :: superquadfacexflags(:, :), &
            superquadfaceyflags(:, :)
        integer(I8)                     :: nx, ny, nt
        integer(I8), allocatable        :: superquadflags(:, :)
        real(R8)                        :: txt, tyt

        ! Loop
        integer(I8)                     :: i, j

        ! Initialize
        !===========
        ! Check inputs
        nx = size(X)
        ny = size(Y)
        nt = size(xt)
        if ((nx /= size(V, 1)) .or. (ny /= size(V, 2))) then 
            call gdErrorHandler('TraceContourStructured2DPoint: incompatible ' // &
                'dimensions of arguments V, X, Y')
        end if 
        if ((size(order) /= size(xs)) .or. (size(xs) /= size(ys)) .or. &
            (size(ys) /= size(vs))) then 
            call gdErrorHandler('TraceContourStructured2DPoint: incompatible ' // &
            'dimensions of arguments xs, ys, vs, order')
        end if 
        if (nt == 0) then 
            ! Empty trace value array - not allowed
            call gdErrorHandler('TraceContourStructured2DPoint: trace value ' // & 
                'array appears to be empty, not supported')
        end if 
        if (nt /= size(yt)) then 
            ! Incompatible dimensions of tracing coordinates
            call gdErrorHandler('TraceContourStructured2DPoint: trace coordinates ' // & 
                'xt and yt have incompatible dimensions, check input')
        end if 

        ! Contours
        allocate(contours(0))

        ! Initialize the saddle point structure
        allocate(superquadflags(nx-1, ny-1), superquadfacexflags(nx, ny-1), &
            superquadfaceyflags(nx-1, ny))
        call InitializeSaddlePointStructure2D(spstruct, superquadflags, &
            superquadfacexflags, superquadfaceyflags, V, X, Y, xs, ys, &
            vs)
        do i = 1, size(spstruct)
            order(i) = spstruct(i)%order
        end do 

        ! Main loop
        !==========
        do i = 1, nt
            ! Get current trace value
            txt = xt(i)
            tyt = yt(i)

            ! Call tracer
            tempcontours = TraceSingleContourStructured2DPoint(V, X, Y, txt, tyt, &
                spstruct, superquadflags, superquadfacexflags, &
                superquadfaceyflags, nx, ny)

            ! Set ID
            do j = 1, size(tempcontours)
                tempcontours(j)%ID = i 
            end do

            ! Add contours
            contours = AddContours(contours, tempcontours)

        end do 
        
    end subroutine

    ! Wrapper for structured contour line tracer
    function TraceContoursValStructured2DWrapper(tracer, tracevalues) result(contours)

        ! Description
        !============
        ! Wrapper for the 2D tracer. The order etc is added to the 
        ! tracer routine

        ! Declare variables
        !==================
        ! Arguments
        class(StructuredContourTracerUDT)       :: tracer 
        real(R8), intent(in)                    :: tracevalues(:)
        type(ContourUDT), allocatable           :: contours(:)

        ! Trace
        !======
        call TraceContoursStructured2D(tracer%V, tracer%X, tracer%Y, &
            tracevalues, tracer%xs, tracer%ys, tracer%vs, tracer%order, &
            contours)

    end function 

    function TraceContoursLocStructured2DWrapper(tracer, x, y) result(contours)
        
        ! Description
        !============
        ! Wrapper for the 2D tracer. The order etc is added to the 
        ! tracer routine

        ! Declare variables
        !==================
        ! Arguments
        class(StructuredContourTracerUDT)       :: tracer 
        real(R8), intent(in)                    :: x(:), y(:)
        type(ContourUDT), allocatable           :: contours(:)

        ! Trace
        !======
        call TraceContoursStructured2DPoint(tracer%V, tracer%X, tracer%Y, &
            x, y, tracer%xs, tracer%ys, tracer%vs, tracer%order, &
            contours)

    end function 

    ! 2D structured single contour line tracer 
    function TraceSingleContourStructured2D(V, X, Y, tv, spstruct,  &
        superquadflags, superquadfacexflags, superquadfaceyflags, nx, ny) &
        result(contours) 

        ! Description
        !============
        ! Trace the contour line(s) for a single value. This function
        ! should not be used standalone, but in conjunction with 
        ! TraceContourStructured2D. This is the main workhorse though 
        ! of the contouring algorithm. For the algorithm, see the 
        ! TraceContourStructured2D routine. 

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)         :: nx, ny
        real(R8), intent(in)            :: V(nx, ny), X(nx), Y(ny), tv 
        integer(I8), intent(in)         :: superquadflags(nx-1, ny-1)
        logical, intent(in)             :: superquadfacexflags(nx, ny-1), &
            superquadfaceyflags(nx-1, ny)
        type(sp2DUDT), intent(in)       :: spstruct(:) 
        type(ContourUDT), allocatable   :: contours(:)

        ! Auxiliary
        logical                         :: hasvv(nx, ny), isexactv(nx, ny), &
            xfacec(nx, ny-1), yfacec(nx-1, ny), issuperquad(nx-1, ny-1), &
            haswestface(nx-1, ny-1), haseastface(nx-1, ny-1), &
            hasnorthface(nx-1, ny-1), hassouthface(nx-1, ny-1), isfound, &
            doexit, saddlequads(nx-1,ny-1), issad, addpoint
        integer(I8)                     :: iif(1:2), jjf(1:2), iic, jjc, &
            bndface, nextquad, siic, sjjc
        integer                         :: tf
        integer(I8), allocatable        :: quadflags(:, :), quadc(:, :)
        real(R8)                        :: V1, V2, x1, y1, x2, &
            y2, frac, tx, ty
        real(R8), allocatable           :: emptyarrayR8(:), Vtrace(:, :)
        type(RealDynamicArrayUDT)       :: xc, yc
        type(ContourUDT)                :: thiscontour
        type(sp2DUDT)                   :: thissp

        ! Loop
        integer(I8)                     :: cc

        ! Initialize
        !===========
        ! Empty array
        allocate(emptyarrayR8(0))

        ! Contour counter
        cc = 0

        ! Check nodal values
        Vtrace = V
        hasvv = Vtrace > tv 
        isexactv = Vtrace == tv 

        ! Check if any values found, if not: add empty contour and 
        ! continue to next value
        issad = any(hasvv)
        if (.not. any(hasvv) .and. (.not. any(isexactv))) then 
            ! Empty contour - construct and return
            allocate(contours(1))
            contours(1) = ConstructContour(emptyarrayR8, &
                emptyarrayR8, tv, 0, 0, 0, .false.)

            ! Return
            return
        else
            allocate(contours(0)) 
        end if 

        ! Values found, check if we need to perturb some locations
        if (any(isexactv)) then 
            ! Issue message 
            if (verbosity > 0) then 
                ! Print message
                print *, 'TraceSingleContourStructured2D: some values ' // & 
                    'lie exactly on a vertex - perturbing the value ' // & 
                    'at these vertices'
            end if 

            ! Perturb values & recompute
            where (isexactv) 
                Vtrace = Vtrace + pert
                hasvv = Vtrace > tv
            end where

        end if 

        ! Determine the quadflags
        quadflags = GetQuadFlags(hasvv, superquadflags)

        ! Determine faceflags (true if value is present on face and if 
        ! face is not an internal face flag - can be logical since we 
        ! can only pass each face once in a contour)
        xfacec = ((hasvv(:, 1:ny-1) .and. (.not. hasvv(:, 2:ny))) .or. &
            ((.not. hasvv(:, 1:ny-1)) .and. (hasvv(:, 2:ny)))) .and. &
            (.not. superquadfacexflags)
        yfacec = ((hasvv(1:nx-1, :) .and. (.not. hasvv(2:nx, :))) .or. &
            ((.not. hasvv(1:nx-1, :)) .and. (hasvv(2:nx, :)))) .and. &
            (.not. superquadfaceyflags)
        
        ! Set number of times we may end up in quad 
        quadc = quadflags
        quadc = 0
        issuperquad = superquadflags > 0 
        where ((quadflags == 1) .and. (.not. issuperquad)) quadc = 1
        where ((quadflags == 4) .and. (.not. issuperquad)) quadc = 4
        where (issuperquad) quadc = 0 ! just to be sure but shouldn't be necessary 
        
        ! Only consider quads of superquads that have an external face of
        ! the superquad on which the value is present
        haswestface = issuperquad .and. xfacec(1:nx-1, :) 
        haseastface = issuperquad .and. xfacec(2:nx, :)  
        hasnorthface = issuperquad .and. yfacec(:, 1:ny-1) 
        hassouthface = issuperquad .and. yfacec(:, 2:ny) 
        where (haswestface) quadc = quadc + 1
        where (haseastface) quadc = quadc + 1
        where (hasnorthface) quadc = quadc + 1
        where (hassouthface) quadc = quadc + 1

        ! Open contours
        !==============
        do while (.true.) 

            ! Get boundary face
            !------------------
            ! Check each of the four boundaries, if none found - no open
            ! contours left
            isfound = .false. 
            bndface = findloc(xfacec(1, :), .true., 1) 
            if (bndface /= 0) then 
                xfacec(1, bndface) = .false. 
                isfound = .true. 
                jjf = [bndface, bndface+1]
                iif = [1, 1]
                jjc = bndface 
                iic = 1
            end if 
            if (.not. isfound) then 
                bndface = findloc(xfacec(nx, :), .true., 1)
                if (bndface /= 0) then 
                    xfacec(nx, bndface) = .false.
                    isfound = .true. 
                    jjf = [bndface, bndface+1]
                    iif = [nx, nx]
                    jjc = bndface 
                    iic = nx-1
                end if 
            end if 
            if (.not. isfound) then 
                bndface = findloc(yfacec(:, 1), .true., 1)
                if (bndface /= 0) then 
                    yfacec(bndface, 1) = .false. 
                    isfound = .true. 
                    jjf = [1, 1]
                    iif = [bndface, bndface+1]
                    jjc = 1 
                    iic = bndface
                end if 
            end if   
            if (.not. isfound) then 
                bndface = findloc(yfacec(:, ny), .true., 1)
                if (bndface /= 0) then 
                    yfacec(bndface, ny) = .false.
                    isfound = .true. 
                    jjf = [ny, ny]
                    iif = [bndface, bndface+1]
                    jjc = ny-1 
                    iic = bndface
                end if 
            end if 
            
            ! Check if a boundary face was found, if not exit loop
            if (.not. isfound) then 
                exit 
            end if 

            ! If we got here, we found a boundary face. Start tracing 

            ! Initialize the new contour
            cc = cc + 1
            thiscontour = ConstructContour(emptyarrayR8, emptyarrayR8, &
                tv, 0, 0, 0, .false.) ! ID is set later

            ! Initialize the x, y coordinates
            xc = ConstructRealDynamicArray()
            yc = ConstructRealDynamicArray()


            ! Remove counter from the current quad 
            if (quadc(iic, jjc) > 0) then 
                quadc(iic, jjc) = quadc(iic, jjc) - 1
            else
                ! This shouldn't happen - throw error
                call gdErrorHandler('TraceSingleContourStructured2D: ' // &
                    ' could not subtract counter from quads after identifying ' // &
                    ' boundary face, this is a bug.')
            end if 

            

            ! Compute the point at this face
            V1 = Vtrace(iif(1), jjf(1))
            V2 = Vtrace(iif(2), jjf(2))
            x1 = X(iif(1)); x2 = X(iif(2))
            y1 = Y(jjf(1)); y2 = Y(jjf(2))
            frac = (tv - V1)/(V2 - V1)
            tx = x1 + frac*(x2 - x1)
            ty = y1 + frac*(y2 - y1)

            ! Add point to contour
            call xc%Append(tx)
            call yc%Append(ty)

            ! Loop until we hit the next boundary or a saddle point
            do while (.true.)
                ! Check if the current cell is a saddle point
                if ((quadflags(iic, jjc) == 4) .and. (.not. issuperquad(iic, jjc))) then 
                    ! Compute point at center and add
                    tx = sum(X([iic, iic+1]))/2.0_R8
                    ty = sum(Y([jjc, jjc+1]))/2.0_R8

                    call xc%Append(tx)
                    call yc%Append(ty)

                    ! Exit loop 
                    exit 
                end if 

                ! Check if the current cell belongs to a predefined 
                ! saddle point
                if (issuperquad(iic, jjc)) then 
                    ! Traverse the saddle point region
                    call TraverseSaddlePoint(iic, jjc, iif, jjf, doexit, &
                        hasvv, quadflags, spstruct, tv, xc, yc, thiscontour, &
                        V, X, Y, nx, ny)

                    ! Exit loop?
                    if (doexit) then 
                        exit 
                    end if 

                    ! Otherwise continue and subtract
                    if (quadc(iic, jjc) > 0) then 
                        quadc(iic, jjc) = quadc(iic, jjc) - 1
                    else
                        ! This shouldn't happen - throw error
                        call gdErrorHandler('TraceSingleContourStructured2D: ' // &
                            ' could not subtract counter from quads after traversing ' // &
                            ' saddle point, this is a bug.')
                    end if 
                end if 

                ! Otherwise, check which faces remain of the cell 
                ! (numbering: 1 (north), 2 (east), 3 (south), 4 (west))
                ! Search for next face
                tf = 0
                isfound = .false.
                if ((.not. isfound) .and. (yfacec(iic, jjc+1))) then 
                    ! North face?
                    tf = 1
                    iif = [iic, iic+1]
                    jjf = [jjc+1, jjc+1]
                    yfacec(iic, jjc+1) = .false.
                    isfound = .true.
                end if
                if ((.not. isfound) .and. (yfacec(iic, jjc))) then 
                    ! South face?
                    tf = 3
                    iif = [iic, iic+1] 
                    jjf = [jjc, jjc] 
                    yfacec(iic, jjc) = .false. 
                    isfound = .true. 
                end if
                if ((.not. isfound) .and. (xfacec(iic+1, jjc))) then
                    ! East face?
                    tf = 2
                    iif = [iic+1, iic+1]
                    jjf = [jjc, jjc+1]
                    xfacec(iic+1, jjc) = .false.
                    isfound = .true.
                end if
                if ((.not. isfound) .and. (xfacec(iic, jjc))) then
                    ! West face?
                    tf = 4 
                    iif = [iic, iic] 
                    jjf = [jjc, jjc+1] 
                    xfacec(iic, jjc) = .false.
                    isfound = .true. 
                end if
                
                ! Check 
                if (.not. isfound) then 
                    call gdErrorHandler('TraceSingleContourStructured2D: ' // &
                        'Could not find next face in open contour, ' // & 
                        'something wrong')
                end if
                
                ! Compute point at this face
                V1 = V(iif(1), jjf(1)); V2 = V(iif(2), jjf(2));
                x1 = X(iif(1)); x2 = X(iif(2));
                y1 = Y(jjf(1)); y2 = Y(jjf(2));
                frac = (tv - V1)/(V2 - V1);
                tx = x1 + frac*(x2 - x1);
                ty = y1 + frac*(y2 - y1);
                
                ! Add point to contour
                call xc%Append(tx)
                call yc%Append(ty)
                
                ! Check if the last face was a boundary face. If so -> exit
                if (iif(1) == iif(2)) then 
                    if ((iif(1) == 1) .or. (iif(1) == nx)) then
                        ! Check which boundary face this is
                        xfacec(iif(1), minval(jjf)) = .false.
                        exit
                    end if 
                end if
                if (jjf(1) == jjf(2)) then 
                    if ((jjf(1) == 1) .or. (jjf(1) == ny)) then 
                        yfacec(minval(iif), jjf(1)) = .false.
                        exit
                    end if 
                end if 
                
                ! If we didn't exit, get the next cell
                if (tf == 1) then 
                    ! Go north
                    jjc = jjc + 1 
                elseif (tf == 3) then 
                    ! Go south
                    jjc = jjc - 1 
                elseif (tf == 2) then 
                    ! Go east
                    iic = iic + 1 
                elseif (tf == 4) then
                    ! Go west
                    iic = iic - 1 
                else
                    call gdErrorHandler('Something wrong')
                end if
                
                ! Subtract
                if (quadc(iic, jjc) > 0) then 
                    quadc(iic, jjc) = quadc(iic, jjc) - 1
                else
                    ! This shouldn't happen - throw error
                    call gdErrorHandler('TraceSingleContourStructured2D: ' // &
                        ' could not subtract counter from quads after identifying ' // &
                        ' next internal face, this is a bug.')
                end if 
                 

            end do 

            ! Add the coordinates etc to the contour
            thiscontour%x = xc%Get()
            thiscontour%y = yc%Get()
            thiscontour%isclosed = .false. 
            contours = AddContours(contours, thiscontour)
            
        end do 

        ! Closed contours
        !================
        ! Here, we need to loop over the remaining quads
        saddlequads = (quadflags == 4) .and. (.not. issuperquad)
        do while (.true.)
            
            ! Check if any saddle point quads remain. If so, start with
            ! these
            nextquad = findloc(reshape((quadc > 0) .and. saddlequads, [size(quadc)]), &
                .true., 1)
            issad = .true. ! is it a saddle point? 
            if (nextquad == 0) then 
                ! Check if any superquad remains
                nextquad = findloc(reshape((quadc > 0) .and. issuperquad, [size(quadc)]),  .true., 1)
                issad = .false.
            end if 
                
            if (nextquad == 0) then 
                ! Check if any quad in general remains
                nextquad = findloc(reshape(quadc > 0, [size(quadc)]), .true., 1)
                issad = .false. 
            end if 
            
            if (nextquad == 0) then 
                ! All found, exit
                exit 
            end if 
            
            ! Initialize the new contour
            cc = cc + 1
            thiscontour = ConstructContour(emptyarrayR8, emptyarrayR8, &
                tv, 0, 0, 0, .false.) ! ID is set later

            ! Initialize the x, y coordinates
            xc = ConstructRealDynamicArray()
            yc = ConstructRealDynamicArray()
            
            ! Get indices
            iic = modulo(nextquad, nx-1)
            jjc = nextquad/(nx-1) + 1
            if (iic == 0) then 
                ! Hedge for edge case where modulo becomes zero 
                iic = nx-1 
                jjc = jjc - 1
            end if
            
            ! Save starting indices
            siic = iic
            sjjc = jjc
            
            ! Subtract counter. If it's a saddle point, add the center of 
            ! this quad
            if (issad) then 
                ! Compute point at center and add
                tx = sum(X([iic, iic+1]))/2.0_R8
                ty = sum(Y([jjc, jjc+1]))/2.0_R8
                call xc%Append(tx)
                call yc%Append(ty)
                quadc(iic, jjc) = quadc(iic, jjc) - 1
            end if 
            if (superquadflags(iic, jjc) > 0) then 
                ! Check if the value is exactly the same as the x-point
                ! value - add the point in that case
                quadc(iic, jjc) = quadc(iic, jjc) - 1
                thissp = spstruct(quadflags(iic, jjc))
                addpoint = .false. 
                if (thissp%val - tv < spvalabstol) then 
                    addpoint = .true. 
                elseif ( (thissp%val - tv)/tv < spvalreltol) then 
                    addpoint = .true.
                end if 
                if (addpoint) then 
                    tx = thissp%x
                    ty = thissp%y 
                    call xc%Append(tx)
                    call yc%Append(ty)
                    thiscontour%startsaddle = thissp%ID
                end if 
            end if 
            
            ! Loop
            do while (.true.)
                ! Find a face that contains the value and for which the
                ! neighbour has quadc > 0 - normally, no quads with the value
                ! exactly on the boundary should be present anymore (boundary
                ! quads where the value lies not on the actual boundary may
                ! still be there). Additionally, no turning back is
                ! allowed!
                
                isfound = .false.                
                ! North face?
                if (yfacec(iic, jjc+1) .and. (.not. isfound)) then 
                    isfound = .true.
                    iif = [iic, iic+1]
                    jjf = [jjc+1, jjc+1]
                    yfacec(iic, jjc+1) = .false.
                    jjc = jjc+1
                end if 
                
                ! South face?
                if (yfacec(iic, jjc) .and. (.not. isfound)) then 
                    isfound = .true.
                    iif = [iic, iic+1]
                    jjf = [jjc, jjc]
                    yfacec(iic, jjc) = .false.
                    jjc = jjc-1
                end if 
                
                ! East face?
                if (xfacec(iic+1, jjc) .and. (.not. isfound)) then 
                    isfound = .true.
                    iif = [iic+1, iic+1]
                    jjf = [jjc, jjc+1]
                    xfacec(iic+1, jjc) = .false.
                    iic = iic + 1
                end if 
                
                ! West face?
                if (xfacec(iic, jjc) .and. (.not. isfound)) then 
                    isfound = .true.
                    iif = [iic, iic]
                    jjf = [jjc, jjc+1]
                    xfacec(iic, jjc) = .false.
                    iic = iic - 1
                end if 
                
                ! check if a cell could be found, otherwise exit
                if (.not. isfound) then 
                    ! Normally no other cell left here, if the starting
                    ! quad is encountered again (and it wasn't a saddle
                    ! point), add the first vertex again
                    if (.not. issad) then 
                        ! Do a sanity check: this should be the first cell
                        if ((iic == siic) .and. (jjc == sjjc)) then 
                            
                        else
                            call gdErrorHandler('TraceSingleContourStructured2D: ' // &
                                'Contour should be closed, but did not ' // &
                                'end up at original cell')
                        end if 
                        tx = xc%Get(1)
                        ty = yc%Get(1)
                        call xc%Append(tx)
                        call yc%Append(ty)
                    else
                        call gdErrorHandler('TraceSingleContourStructured2D: ' // &
                            'unknown error, this is likely a bug')
                    end if 
                    
                    ! Check if the contour was closed
                    if ((iic == siic) .and. (jjc == sjjc)) then 
                        thiscontour%isclosed = .true.
                    end if 
                
                    exit 
                end if 
                
                ! Subtract
                if (quadc(iic, jjc) == 0) then 
                    call gdErrorHandler('TraceSingleContourStructured2D: ' // & 
                        'could not subtract counter, this is a bug')
                end if 
                quadc(iic, jjc) = quadc(iic, jjc) - 1 
                
                ! Compute point at this face
                V1 = V(iif(1), jjf(1)) 
                V2 = V(iif(2), jjf(2))
                x1 = X(iif(1)) 
                x2 = X(iif(2))
                y1 = Y(jjf(1)) 
                y2 = Y(jjf(2))
                frac = (tv - V1)/(V2 - V1)
                tx = x1 + frac*(x2 - x1)
                ty = y1 + frac*(y2 - y1)
                
                ! Add point to contour
                call xc%Append(tx)
                call yc%Append(ty)
                
                ! Check if the next quad contains a saddle point -> add and
                ! stop
                if (saddlequads(iic, jjc)) then 
                    ! Compute point at center and add
                    tx = sum(X([iic, iic+1]))/2.0_R8
                    ty = sum(Y([jjc, jjc+1]))/2.0_R8
                    call xc%Append(tx)
                    call yc%Append(ty)
                    
                    ! Exit
                    exit 
                end if 
                
                ! Check if the next quad is part of a superquad -> traverse
                ! through, exit if necessary
                if (superquadflags(iic, jjc) > 0) then 
                    ! Traverse the saddle point
                    call TraverseSaddlePoint(iic, jjc, iif, jjf, doexit, &
                        hasvv, quadflags, spstruct, tv, xc, yc, thiscontour, &
                        V, X, Y, nx, ny)

                    ! Exit?
                    if (doexit) then 
                        exit 
                    end if 
                    
                    ! Subtract
                    if (quadc(iic, jjc) == 0) then 
                        call gdErrorHandler('TraceSingleContourStructured2D: ' // & 
                        'could not subtract counter, this is a bug')
                    end if 
                    quadc(iic, jjc) = quadc(iic, jjc) - 1
                end if 
            end do 

            ! Add contours
            thiscontour%x = xc%Get()
            thiscontour%y = yc%Get()
            contours = AddContours(contours, thiscontour)

        end do



    end function

    ! 2D structured single contour line tracer starting from point
    function TraceSingleContourStructured2DPoint(V, X, Y, x0, y0, spstruct,  &
        superquadflags, superquadfacexflags, superquadfaceyflags, nx, ny) &
        result(contours) 

        ! Description
        !============
        ! Trace the contour line(s) for a single point. This function
        ! should not be used standalone, but in conjunction with 
        ! TraceContourStructured2DPoint. This is the main workhorse though 
        ! of the contouring algorithm. For the algorithm, see the 
        ! TraceContourStructured2DPoint routine. Note that this routine is 
        ! slightly different than the routine that starts from a given 
        ! value, since we can only start/end in one particular quad. 

        ! Declare variables
        !==================
        ! Modules
        use mod_constants, only: posinfval_R8
        ! Arguments
        integer(I8), intent(in)         :: nx, ny
        real(R8), intent(in)            :: V(nx, ny), X(nx), Y(ny), x0, y0 
        integer(I8), intent(in)         :: superquadflags(nx-1, ny-1)
        logical, intent(in)             :: superquadfacexflags(nx, ny-1), &
            superquadfaceyflags(nx-1, ny)
        type(sp2DUDT), intent(in)       :: spstruct(:) 
        type(ContourUDT), allocatable   :: contours(:)

        ! Auxiliary
        logical                         :: hasvv(nx, ny), isexactv(nx, ny), &
            xfacec(nx, ny-1), yfacec(nx-1, ny), issuperquad(nx-1, ny-1), &
            haswestface(nx-1, ny-1), haseastface(nx-1, ny-1), &
            hasnorthface(nx-1, ny-1), hassouthface(nx-1, ny-1), isfound, &
            doexit, saddlequads(nx-1,ny-1), issad, addpoint, issaddlepoint, &
            startquads(nx-1, ny-1), isstartperturbed
        integer(I8)                     :: iif(1:2), jjf(1:2), iic, jjc, &
            nextquad, siic, sjjc, iisq, jjsq, saddlepointID, xloc, yloc
        integer(I8), allocatable        :: quadflags(:, :), quadc(:, :)
        real(R8)                        :: V1, V2, x1, y1, x2, &
            y2, frac, tx, ty, tv1, tv2, tv, x0p, y0p, xpert, ypert
        real(R8), allocatable           :: emptyarrayR8(:), Vtrace(:, :)
        type(RealDynamicArrayUDT)       :: xc, yc
        type(ContourUDT)                :: thiscontour
        type(sp2DUDT)                   :: thissp

        ! Loop
        integer(I8)                     :: i, cc

        ! Initialize
        !===========
        ! Empty array
        allocate(emptyarrayR8(0))

        ! Contour counter
        cc = 0

        ! Hedge for starting indices that lie exactly on a node 
        isstartperturbed = .false.
        x0p = x0 
        y0p = y0
        xloc = findloc(X, x0, 1) 
        yloc = findloc(Y, y0, 1)
        if ( (xloc /= 0) .and. ( yloc /= 0)) then 
            ! Only perturb to compute initial value, then set equal 
            ! again to original starting point
            isstartperturbed = .true.

            ! Compute perturbation 
            x0p = x0 + pert
            y0p = y0 + pert  
        end if 

        ! Find the starting quad indices
        iisq = findloc(x0 > X, .true., 1, back=.true.)
        jjsq = findloc(y0 > Y, .true., 1, back=.true.)

        ! Hedge for out of bounds
        if ((iisq == 0) .or. (jjsq == 0) .or. (iisq == nx) .or. (jjsq == ny)) then 
            ! Out of bounds, make empty contour and return 
            allocate(contours(1))
            contours(1) = ConstructContour(emptyarrayR8, &
                emptyarrayR8, posinfval_R8(), 0, 0, 0, .false.)
            return 
        else 
            allocate(contours(0))
        end if 

        ! Compute value at starting quad by bilinear interpolation
        tv1 = (V(iisq+1, jjsq) - V(iisq, jjsq))/(X(iisq+1) - X(iisq))*(x0p - X(iisq)) + V(iisq, jjsq)
        tv2 = (V(iisq+1, jjsq+1) - V(iisq, jjsq+1))/(X(iisq+1) - X(iisq))*(x0p - X(iisq)) + V(iisq, jjsq+1)
        tv = (tv2 - tv1)/(Y(jjsq+1) - Y(jjsq))*(y0p - Y(jjsq)) + tv1
        
        ! Check if the current point is an x-point. 
        issaddlepoint = .false.
        saddlepointID = 0
        do i = 1, size(spstruct)
            if ((x0 == spstruct(i)%x) .and. (y0 == spstruct(i)%y)) then 
                issaddlepoint = .true.
                saddlepointID = i
            end if 
        end do 
        
        ! Overwrite the trace value
        if (issaddlepoint) then 
            tv = spstruct(saddlepointID)%val
        end if 

        ! Check nodal values
        Vtrace = V
        hasvv = Vtrace > tv 
        isexactv = Vtrace == tv 

        ! Values found, check if we need to perturb some locations
        if (any(isexactv)) then 
            ! Issue message 
            if (verbosity > 0) then 
                ! Print message
                print *, 'TraceSingleContourStructured2DPoint: some values ' // & 
                    'lie exactly on a vertex - perturbing the value ' // & 
                    'at these vertices'
            end if 

            ! Perturb values & recompute
            where (isexactv) 
                Vtrace = Vtrace + pert
                hasvv = Vtrace > tv
            end where

        end if 

        ! Determine the quadflags
        quadflags = GetQuadFlags(hasvv, superquadflags)

        ! Determine faceflags (true if value is present on face and if 
        ! face is not an internal face flag - can be logical since we 
        ! can only pass each face once in a contour)
        xfacec = ((hasvv(:, 1:ny-1) .and. (.not. hasvv(:, 2:ny))) .or. &
            ((.not. hasvv(:, 1:ny-1)) .and. (hasvv(:, 2:ny)))) .and. &
            (.not. superquadfacexflags)
        yfacec = ((hasvv(1:nx-1, :) .and. (.not. hasvv(2:nx, :))) .or. &
            ((.not. hasvv(1:nx-1, :)) .and. (hasvv(2:nx, :)))) .and. &
            (.not. superquadfaceyflags)
        
        ! Set number of times we may end up in quad 
        quadc = quadflags
        startquads = .false. 
        quadc = 0
        issuperquad = superquadflags > 0 
        where ((quadflags == 1) .and. (.not. issuperquad)) quadc = 1
        where ((quadflags == 4) .and. (.not. issuperquad)) quadc = 4
        where (issuperquad) quadc = 0 ! just to be sure but shouldn't be necessary 
        
        ! Only consider quads of superquads that have an external face of
        ! the superquad on which the value is present
        haswestface = issuperquad .and. xfacec(1:nx-1, :) 
        haseastface = issuperquad .and. xfacec(2:nx, :)  
        hasnorthface = issuperquad .and. yfacec(:, 1:ny-1) 
        hassouthface = issuperquad .and. yfacec(:, 2:ny) 
        where (haswestface) quadc = quadc + 1
        where (haseastface) quadc = quadc + 1
        where (hasnorthface) quadc = quadc + 1
        where (hassouthface) quadc = quadc + 1

        ! Determine starting quads
        if (issuperquad(iisq, jjsq)) then 
            ! Set all superquad cells to true with this flag
            where ( (quadc > 0) .and. (superquadflags == superquadflags(iisq, jjsq)) ) startquads = .true.
        else 
            ! Only current cell is start quad
            startquads(iisq, jjsq) = .true. 
        end if 

        ! Starting quads that are not saddle point quads should have two
        ! starts
        saddlequads = (quadflags == 4) .and. (.not. issuperquad)
        where (startquads .and. (.not. saddlequads) .and. (.not. issuperquad)) quadc = 2

        ! All contours
        !=============    
        do while (.true.)
            
            ! Check if any saddle point quads remain. If so, start with
            ! these
            nextquad = findloc(reshape((quadc > 0) .and. saddlequads &
                .and. startquads .and. (.not. issuperquad), [size(quadc)]), &
                .true., 1)
            issad = .true. ! is it a saddle point? 
            if (nextquad == 0) then 
                ! Check if any superquad remains
                nextquad = findloc(reshape((quadc > 0) .and. startquads &
                .and. issuperquad, [size(quadc)]),  .true., 1)
                issad = .false.
            end if 
            if (nextquad == 0) then 
                ! Check if any quad in general remains
                nextquad = findloc(reshape((quadc > 0) .and. (startquads), [size(quadc)]), .true., 1)
                issad = .false. 
            end if 
            if (nextquad == 0) then 
                ! All found, exit
                exit 
            end if 
            
            ! Initialize the new contour
            cc = cc + 1
            thiscontour = ConstructContour(emptyarrayR8, emptyarrayR8, &
                tv, 0, 0, 0, .false.) ! ID is set later

            ! Initialize the x, y coordinates
            xc = ConstructRealDynamicArray()
            yc = ConstructRealDynamicArray()
            
            ! Get indices
            iic = modulo(nextquad, nx-1)
            jjc = nextquad/(nx-1) + 1
            if (iic == 0) then 
                ! Hedge for edge case where modulo becomes zero 
                iic = nx-1 
                jjc = jjc - 1
            end if

            ! Check if we should add the starting point
            if ((iic == iisq) .and. (jjc == jjsq)) then 
                call xc%Append(x0)
                call yc%Append(y0)
            end if 
            
            ! Save starting indices
            siic = iic
            sjjc = jjc
            
            ! Subtract counter. If it's a saddle point, add the center of 
            ! this quad
            if (issad) then 
                ! Compute point at center and add
                tx = sum(X([iic, iic+1]))/2.0_R8
                ty = sum(Y([jjc, jjc+1]))/2.0_R8
                call xc%Append(tx)
                call yc%Append(ty)
                quadc(iic, jjc) = quadc(iic, jjc) - 1
            end if 
            if (superquadflags(iic, jjc) > 0) then 
                ! Check if the value is exactly the same as the x-point
                ! value - add the point in that case
                quadc(iic, jjc) = quadc(iic, jjc) - 1
                thissp = spstruct(quadflags(iic, jjc))
                addpoint = .false. 
                if (thissp%val - tv < spvalabstol) then 
                    addpoint = .true. 
                elseif ( (thissp%val - tv)/tv < spvalreltol) then 
                    addpoint = .true.
                end if 
                if (addpoint) then 
                    tx = thissp%x
                    ty = thissp%y 
                    call xc%Append(tx)
                    call yc%Append(ty)
                    thiscontour%startsaddle = thissp%ID
                end if 
            end if 
            
            ! Loop
            do while (.true.)
                ! Find a face that contains the value and for which the
                ! neighbour has quadc > 0 - normally, no quads with the value
                ! exactly on the boundary should be present anymore (boundary
                ! quads where the value lies not on the actual boundary may
                ! still be there). Additionally, no turning back is
                ! allowed!
                
                isfound = .false.      
                ! North face?
                if (yfacec(iic, jjc+1) .and. (.not. isfound)) then 
                    isfound = .true.
                    iif = [iic, iic+1]
                    jjf = [jjc+1, jjc+1]
                    yfacec(iic, jjc+1) = .false.
                    jjc = jjc+1
                end if 
                
                ! South face?
                if (.not. isfound) then 
                    if (yfacec(iic, jjc)) then 
                        isfound = .true.
                        iif = [iic, iic+1]
                        jjf = [jjc, jjc]
                        yfacec(iic, jjc) = .false.
                        jjc = jjc-1
                    end if 
                end if 
                
                ! East face?
                if (.not. isfound) then 
                    if (xfacec(iic+1, jjc)) then 
                        isfound = .true.
                        iif = [iic+1, iic+1]
                        jjf = [jjc, jjc+1]
                        xfacec(iic+1, jjc) = .false.
                        iic = iic + 1
                    end if 
                end if 
                
                ! West face?
                if (.not. isfound) then 
                    if (xfacec(iic, jjc)) then 
                        isfound = .true.
                        iif = [iic, iic]
                        jjf = [jjc, jjc+1]
                        xfacec(iic, jjc) = .false.
                        iic = iic - 1
                    end if 
                end if 
                
                ! check if a cell could be found, otherwise exit
                if (.not. isfound) then 
                    ! Normally no other cell left here, if the starting
                    ! quad is encountered again (and it wasn't a saddle
                    ! point), add the first vertex again
                    if (.not. issad) then 
                        ! Do a sanity check: this should be the first cell
                        if ((iic == siic) .and. (jjc == sjjc)) then 
                            ! Subtract
                            quadc(iic, jjc) = quadc(iic, jjc) - 1

                            ! Set isclosed
                            thiscontour%isclosed = .true. 
                            
                        else
                            call gdErrorHandler('TraceSingleContourStructured2D: ' // &
                                'Contour should be closed, but did not ' // &
                                'end up at original cell')
                        end if 
                        tx = xc%Get(1)
                        ty = yc%Get(1)
                        call xc%Append(tx)
                        call yc%Append(ty)
                    else
                        call gdErrorHandler('TraceSingleContourStructured2D: ' // &
                            'unknown error, this is likely a bug')
                    end if 
                
                    exit 
                end if 
                
                ! Compute point at this face
                V1 = V(iif(1), jjf(1)) 
                V2 = V(iif(2), jjf(2))
                x1 = X(iif(1)) 
                x2 = X(iif(2))
                y1 = Y(jjf(1)) 
                y2 = Y(jjf(2))
                frac = (tv - V1)/(V2 - V1)
                tx = x1 + frac*(x2 - x1)
                ty = y1 + frac*(y2 - y1)

                ! Check if we hit a boundary face
                if ((iic == nx) .or. (iic == 0) .or. (jjc == ny) .or. (jjc == 0)) then 
                    ! Exit
                    exit
                end if 

                ! Subtract
                quadc(iic, jjc) = quadc(iic, jjc) - 1
                
                ! Add point to contour
                call xc%Append(tx)
                call yc%Append(ty)
                
                ! Check if the next quad contains a saddle point -> add and
                ! stop
                if (saddlequads(iic, jjc)) then 
                    ! Compute point at center and add
                    tx = sum(X([iic, iic+1]))/2.0_R8
                    ty = sum(Y([jjc, jjc+1]))/2.0_R8
                    call xc%Append(tx)
                    call yc%Append(ty)
                    
                    ! Exit
                    exit 
                end if 
                
                ! Check if the next quad is part of a superquad -> traverse
                ! through, exit if necessary
                if (superquadflags(iic, jjc) > 0) then 
                    ! Traverse the saddle point
                    call TraverseSaddlePoint(iic, jjc, iif, jjf, doexit, &
                        hasvv, quadflags, spstruct, tv, xc, yc, thiscontour, &
                        V, X, Y, nx, ny)

                    ! Exit?
                    if (doexit) then 
                        ! Mark this saddle point as potential starting quads
                        where ((superquadflags == superquadflags(iic, jjc)) &
                            .and. (quadc > 0)) startquads = .true.  
                        exit 
                    end if 
                    
                    ! Subtract
                    if (quadc(iic, jjc) == 0) then 
                        call gdErrorHandler('TraceSingleContourStructured2D: ' // & 
                        'could not subtract counter, this is a bug')
                    end if 
                    quadc(iic, jjc) = quadc(iic, jjc) - 1
                end if 
            end do 

            ! Add contours
            thiscontour%x = xc%Get()
            thiscontour%y = yc%Get()
            contours = AddContours(contours, thiscontour)

        end do

    end function

    !------------------------------------------------------------------!
    !                             AUXILIARY                            !
    !------------------------------------------------------------------!

    ! Contour adding routine
    function AddContourArray(contours1, contours2) result(newcontours)

        ! Description
        !============
        ! Simple routine to add contour structs together. Contours2 are 
        ! appended to the end of contours1. Note: may not be the most 
        ! memory-friendly way of handling things, but we assume that
        ! this cost is negligible compared to e.g. computing the contour
        ! itself (which should typically be the case). 

        ! Declare variables
        !==================
        ! Arguments
        type(ContourUDT), intent(in)    :: contours1(:), contours2(:)
        type(ContourUDT), allocatable   :: newcontours(:)

        ! Add
        !====
        ! Allocate
        allocate(newcontours(size(contours1) + size(contours2)))

        ! Attribute
        newcontours(1:size(contours1)) = contours1 
        newcontours(size(contours1)+1:size(contours2)+size(contours1)) = &
            contours2

    end function

    function AddContourScalar(contours1, contours2) result(newcontours)

        ! Description
        !============
        ! Simple routine to add contour structs together. Contours2 are 
        ! appended to the end of contours1. Note: may not be the most 
        ! memory-friendly way of handling things, but we assume that
        ! this cost is negligible compared to e.g. computing the contour
        ! itself (which should typically be the case). 

        ! Declare variables
        !==================
        ! Arguments
        type(ContourUDT), intent(in)    :: contours1(:), contours2
        type(ContourUDT), allocatable   :: newcontours(:)

        ! Add
        !====
        ! Allocate
        allocate(newcontours(size(contours1) + 1))

        ! Attribute
        newcontours(1:size(contours1)) = contours1 
        newcontours(size(contours1)+1) = contours2

    end function 

    ! Saddle point structure initializer for 2D tracer
    subroutine InitializeSaddlePointStructure2D(spstruct, &
        quadflags, facexflags, faceyflags, &
        V, X, Y, xs, ys, vs)

        ! Description
        !============
        ! Build the saddle point structure. Saddle points that lie outside of the 
        ! domain are ignored. For each saddle point, a simple triangular 'mesh' is
        ! constructed, where the first node is always the x-point and where
        ! triangles are sorted in a clockwise fashion, starting from the bottom
        ! left node (ixquad, iyquad). In each direction, an additional number of
        ! cells equal to 'npc' is added, and the vertices that form the boundary of
        ! this 'superquad' are considered for the triangles. Therefore, if npc = 1,
        ! 13 vertices will be present in the mesh (12 for the outermost boundary, 1
        ! for the x-point). For cells close to the boundary, we shift this stencil 
        ! such that it fits and that the same number of cells is taken. 
        ! In this way, the algorithm can proceed but
        ! the results may be inaccurate due to lower local resolution near the
        ! boundary. Only refinement can alleviate this. It may also prove difficult
        ! to determine the exact topology near the boundary, since boundary faces
        ! may not contain the x-point value at some places, leading to improper
        ! order determination. 

        ! Notes 
        !======
        ! Note 1: overlapping saddle point domains are not allowed, adjacent should
        ! be possible. If this poses to be an issue, refine the grid locally or
        ! remove one of the x-points (at the expense of possible inaccuracies).

        ! Note 2: we're now hedging more for wrongly identified saddle
        ! points, which will not be included in the spstruct. This 
        ! identification is based firstly on location (outside of bounds
        ! ) and secondly on the field value. However, we cannot fully
        ! hedge for wrongly passed saddle points, so errors may still 
        ! occur... 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)                :: xs(:), ys(:), vs(:), &
            X(:), Y(:), V(:, :)
        type(sp2DUDT), allocatable, intent(out)     :: spstruct(:)
        integer(I8), intent(out)            :: quadflags(size(X)-1, size(Y)-1)
        logical, intent(out)                :: facexflags(size(X), size(Y)-1), &
            faceyflags(size(X)-1, size(Y))
        
        ! Auxiliary
        integer(I8)                         :: ixquad, iyquad, ntri, &
            nx, ny, m
        integer(I8), allocatable            :: stencilx(:), stencily(:)
        logical, allocatable                :: hasvp(:), hasftri(:), &
            keepind(:)

        ! Loop  
        integer(I8)                         :: i, k, qfc

        ! Initialize
        !===========
        ! Quadflag counter - to account for possible deletion of 
        ! saddle points
        qfc = 0

        ! Set sizes
        nx = size(X)
        ny = size(Y)

        ! Set number of triangles
        ntri = (2*npq+1)*4
        
        ! Initialize flags
        quadflags = 0 
        facexflags = .false. 
        faceyflags = .false.

        ! Initialize deletion index & spstruct
        allocate(keepind(size(xs)), spstruct(size(xs)))
        keepind = .true. 

        ! Set up saddle point structure
        !==============================
        do i = 1, size(xs)

            ! Determine saddle point location
            ixquad = findloc(xs(i) > X, .true., dim=1, back=.true.)
            iyquad = findloc(ys(i) > Y, .true., dim=1, back=.true.)

            ! Check for out-of-bounds
            if ((ixquad == 0) .or. (iyquad == 0)) then 

                ! Print message and skip
                print *, 'InitializeSaddlePointStructure2D: saddle point ' // &
                    'location lies out of bounds, ignoring saddle point ' // &
                    'number: ', i 

                ! Mark for deletion
                keepind(i) = .false. 

                ! Skip
                cycle 
            end if 

            ! Add
            spstruct(i)%ID      = i  
            spstruct(i)%x       = xs(i)
            spstruct(i)%y       = ys(i) 
            spstruct(i)%val     = vs(i)
            spstruct(i)%ixquad  = ixquad 
            spstruct(i)%iyquad  = iyquad 
            
            ! Initialize
            
            allocate(spstruct(i)%ixpoints(ntri+1), &
                spstruct(i)%iypoints(ntri+1), &
                spstruct(i)%valpoints(ntri+1), &
                spstruct(i)%tri(ntri, 3))
            
            ! Triangles
            spstruct(i)%tri(:, 1) = 1  ! first point is saddle point
            spstruct(i)%valpoints(1) = spstruct(i)%val 
            spstruct(i)%tri(:, 2) = [(k, k = 2, ntri+1)] 
            spstruct(i)%tri(:, 3) = [(k, k = 3, ntri+1), 2] 
            
            ! Points
            stencilx = [(k, k = -npq, npq)] + ixquad  ! stencil for cells
            stencily = [(k, k = -npq, npq)] + iyquad 
            if (stencilx(1) < 1) then 
                stencilx = stencilx + 1 - stencilx(1)  
            end if
            if (stencilx(2*npq+1) > nx-1) then 
                stencilx = stencilx - (stencilx(2*npq+1) - (nx - 1) ) 
            end if 
            if (stencily(1) < 1) then 
                stencily = stencily + 1 - stencily(1)  
            end if
            if (stencily(2*npq+1) > ny-1) then 
                stencily = stencily - (stencily(2*npq+1) - (ny - 1) ) 
            end if
            
            spstruct(i)%ixpoints = [0, spread(stencilx(1), 1, 2*npq+1), &
                stencilx, spread(stencilx(2*npq+1)+1, 1, 2*npq+2), &
                stencilx((2*npq+1):2:-1)] 
            spstruct(i)%iypoints = [0, stencily, &
                spread(stencily(2*npq+1)+1, 1, 2*npq+2), &
                stencily((2*npq+1):1:-1), spread(stencily(1), 1, 2*npq)]
            do k = 2, ntri+1
                spstruct(i)%valpoints(k) = V(spstruct(i)%ixpoints(k), spstruct(i)%iypoints(k)) 
            end do
                    
            ! Check
            if (any(any(quadflags(stencilx, stencily) > 0, 1))) then
                print *, 'Saddle point number: ', i  
                call gdErrorHandler('InitializeSaddlePointStructure: ' // &
                 'saddle point domains overlap, not supported. ' // & 
                 'Consider refining the grid or removing saddle points. ')
            end if

            ! Determine x-point order
            hasvp = spstruct(i)%valpoints >= spstruct(i)%val 
            hasftri = (hasvp(spstruct(i)%tri(:, 2)) .and. &
                (.not. hasvp(spstruct(i)%tri(:, 3)))) .or. &
                (hasvp(spstruct(i)%tri(:, 3)) .and. &
                (.not. (hasvp(spstruct(i)%tri(:, 2))))) 
            m = count(hasftri)  ! number of intersections with outer boundary
            
            ! Check
            if (modulo(m, 2) > 0) then 
                ! Something weird going on
                print *, 'InitializeSaddlePointStructure: order of saddle ' // &
                    'point with ID ', i, ' could not be determined, ' // &
                    'returning NaN for this order and not including as '// & 
                    'saddle point for tracing'
                
                ! Mark for deletion
                keepind(i) = .false. 

                ! Skip
                cycle 
            end if
            if (m == 0) then 
                ! This may happen when e.g. an extremum is given as a 
                ! pseudo saddlepoint. Ignore it, but don't issue warning
                !print *, 'InitializeSaddlePointStructure: given saddle ' // & 
                !    'point value is not present on any of the saddle ' // & 
                !    'point domain boundaries - please check input value.' // &
                !    'Ignoring saddle point ', i

                ! Mark for deletion
                keepind(i) = .false. 
                cycle 
            end if 
            
            ! Compute
            spstruct(i)%order = m/2-1   

            ! Update quadflag counter
            qfc = qfc + 1
            
            ! Set quadflags
            quadflags(stencilx, stencily) = qfc
            
            ! Set faceflags - only for inner faces
            facexflags(stencilx(2:size(stencilx)), stencily) = .true. 
            faceyflags(stencilx, stencily(2:size(stencily))) = .true. 
            
            ! Set quad indices for triangles
            spstruct(i)%ixquadtri = [spread(stencilx(1), 1, 2*npq+1), &
                stencilx, spread(stencilx(2*npq+1), 1, 2*npq+1), &
                stencilx(2*npq+1:1:-1)]  
            spstruct(i)%iyquadtri = [stencily, &
                spread(stencily(2*npq+1), 1, 2*npq+1), &
                stencily((2*npq+1):1:-1), spread(stencily(1), 1, 2*npq+1)]  

        end do 

        ! Delete if necessary
        spstruct = pack(spstruct, keepind)


    end subroutine 

    ! Saddle point traversal 
    subroutine TraverseSaddlePoint(iic, jjc, ixv, iyv, doexit, &
        hasvv, quadflags, spstruct, tv, xc, yc, contour, V, X, Y, &
        nx, ny)

        ! Description
        !============
        ! Traverse a saddle point region characterized by the saddle point
        ! index 'isp'. The algorithm is as follows:
        ! - check which nodes are higher/lower than the current value. If
        ! the value is exactly equal to the saddle point value, signal to
        ! exit the loop ('exactly' is rather up to spvalabstol, as defined
        ! in the module above or spvalreltol for relative difference)
        ! - if we didn't exit, continue tracing until we encounter the
        ! first pair of nodes that is fully external to the x-point. Return
        ! the quad indices corresponding to that point as output.

        ! Note: if the saddle point value has been wrongly determined,
        ! it is possible that this routine cannot proceed. 

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(inout)      :: iic, jjc
        integer(I8), intent(in)         :: ixv(1:2), iyv(1:2), quadflags(:, :), &
            nx, ny
        logical, intent(in)             :: hasvv(nx, ny)
        logical, intent(inout)          :: doexit 
        type(sp2DUDT), intent(in)       :: spstruct(:)
        type(ContourUDT), intent(inout) :: contour
        type(RealDynamicArrayUDT), intent(inout)    :: xc, yc
        real(R8), intent(in)            :: tv, V(nx, ny), X(nx), Y(ny)

        ! Auxiliary
        type(sp2DUDT)                   :: thissp
        logical, allocatable            :: hasvsp(:), hasftri(:, :), &
            hasvtri(:, :)
        integer(I8)                     :: isp, thisiic, thisjjc, &
            nodes(1:2), ctri, p2
        real(R8)                        :: tx, ty, V1, V2, x1, x2, y1, &
            y2, frac 

        ! Loop 
        integer(I8)                     :: i

        ! Initialize
        !===========
        ! Get saddlepoint index
        isp = quadflags(iic, jjc)
        
        ! Initialize
        thisiic = iic
        thisjjc = jjc
        doexit = .false.
        
        ! Unpack
        thissp = spstruct(isp)

        ! Check if we need to skip 
        if ( abs(tv - thissp%val) <= spvalabstol) then 
            ! Exit
            doexit = .true.
        else
            ! It should be safe to check for the relative error now, 
            ! otherwise we would've exited due to the first statement 
            ! already
            if (abs((tv - thissp%val)/tv) <= spvalreltol) then 
                doexit = .true.
            end if 
        end if 
        
        ! Check
        if (doexit) then 

            ! Add the saddle point location
            tx = thissp%x
            ty = thissp%y

            ! Add to contour
            call xc%Append(tx)
            call yc%Append(ty)
            contour%endsaddle = thissp%ID
            return
        end if
        
        ! Check values
        allocate(hasvsp(size(thissp%ixpoints)), &
            hasvtri(size(thissp%tri, 1), size(thissp%tri, 2)))
        hasvsp(1) = thissp%val > tv 
        do i = 2, size(thissp%ixpoints)
            hasvsp(i) = hasvv(thissp%ixpoints(i), thissp%iypoints(i))
        end do
        do i = 1, size(thissp%tri, 2)
            hasvtri(:, i) = hasvsp(thissp%tri(:, i))
        end do 
        hasftri = hasvtri 
        hasftri = .false. ! indicator for faces of triangle
        hasftri(:, 1) = (hasvtri(:, 1) .and. (.not. hasvtri(:, 2))) .or.&
            (hasvtri(:, 2) .and. (.not. hasvtri(:, 1)))
        hasftri(:, 2) = (hasvtri(:, 2) .and. (.not. hasvtri(:, 3))) .or.&
            (hasvtri(:, 3) .and. (.not. hasvtri(:, 2)))
        hasftri(:, 3) = (hasvtri(:, 1) .and. (.not. hasvtri(:, 3))) .or.&
            (hasvtri(:, 3) .and. (.not. hasvtri(:, 1)))
        
        ! Trace
        !======
        ! Determine initial triangle
        nodes(1) = findloc( (thissp%ixpoints == ixv(1)) .and. (thissp%iypoints == iyv(1)), .true., 1)
        nodes(2) = findloc( (thissp%ixpoints == ixv(2)) .and. (thissp%iypoints == iyv(2)), .true., 1)
        if ((maxval(nodes)-minval(nodes)) == 1) then 
            ctri = minval(nodes)-1
        else
            ctri = size(thissp%tri, 1)
        end if
        
        ! Traverse through the saddle point region
        hasftri(ctri, 2) = .false.
        do while (.true.)
            ! Consistency check
            if (count(hasftri(ctri, :)) > 1) then 
                call gdErrorHandler('TraverseSaddlePoint: Too many ' // &
                    'values found in saddle point region - check if ' // & 
                    'saddle points were correctly determined')
            elseif (count(hasftri(ctri, :)) < 1) then 
                call gdErrorHandler('TraverseSaddlePoint: Not enough values' // & 
                    'found in saddle point region - check if ' // & 
                    'saddle points were correctly determined')
            end if
            
            ! Find the next face that contains the value
            if (hasftri(ctri, 1)) then 
                ! Compute point at this face
                p2 = thissp%tri(ctri, 2)
                V1 = thissp%val 
                V2 = V(thissp%ixpoints(p2), thissp%iypoints(p2))
                x1 = thissp%x 
                x2 = X(thissp%ixpoints(p2))
                y1 = thissp%y 
                y2 = Y(thissp%iypoints(p2))
                frac = (tv - V1)/(V2 - V1)
                tx = x1 + frac*(x2 - x1)
                ty = y1 + frac*(y2 - y1)
                
                ! Add point to contour
                call xc%Append(tx)
                call yc%Append(ty)
            
                ! Previous triangle, update
                hasftri(ctri, 1) = .false.
                ctri = ctri-1
                if (ctri == 0) then 
                    ctri = size(thissp%tri, 1)
                end if
                
                ! Remove from next triangle
                hasftri(ctri, 3) = .false. 
                
            elseif (hasftri(ctri, 2)) then 
                ! Current triangle, exterior face -> exit (point will
                ! be added in overarching tracing routine later)
                ! hasftri(ctri, 2) = false
                exit
            elseif (hasftri(ctri, 3)) then 
                ! Compute point at this face
                p2 = thissp%tri(ctri, 3)
                V1 = thissp%val 
                V2 = V(thissp%ixpoints(p2), thissp%iypoints(p2))
                x1 = thissp%x 
                x2 = X(thissp%ixpoints(p2))
                y1 = thissp%y 
                y2 = Y(thissp%iypoints(p2))
                frac = (tv - V1)/(V2 - V1)
                tx = x1 + frac*(x2 - x1)
                ty = y1 + frac*(y2 - y1)
                
                ! Add point to contour
                call xc%Append(tx)
                call yc%Append(ty)
                
                ! Next triangle, update
                hasftri(ctri, 3) = .false.
                ctri = ctri+1
                if (ctri > size(thissp%tri, 1)) then 
                    ctri = 1
                end if 
                
                ! Remove from next triangle
                hasftri(ctri, 1) = .false.
            else
                call gdErrorHandler('TraverseSaddlePoint: Could not ' // & 
                    'find next face in saddle point region, this is ' // & 
                    'likely a bug in the code')
            end if
        end do
        
        ! Compute the quad we ended up in
        iic = thissp%ixquadtri(ctri)
        jjc = thissp%iyquadtri(ctri)
        
    end subroutine 

    ! Quad flag determination
    function GetQuadFlags(hasvv, superquadflags) result(quadflags)

        ! Description
        !============
        ! Determine which flag the quad should get:
        ! - 0: no values present
        ! - 1: regular quad with value present on two sides
        ! - 4: quad contains a saddle point up to resolution of the grid
        ! - ?: quad belongs to a superquad of an x-point (number equals
        ! index of saddle point)
        
        ! As input, the logical 2D array hasvv should be passed, which 
        ! should be true at vertex locations where the value is higher 
        ! than the to be traced value (lower will likely also work but
        ! ok)

        ! Declare variables
        !==================
        ! Arguments
        logical, intent(in)         :: hasvv(:, :)
        integer(I8), intent(in)     :: superquadflags(size(hasvv, 1)-1, size(hasvv, 2)-1)
        integer(I8)                 :: nv1, nv2, &
            quadflags(size(hasvv, 1)-1, size(hasvv, 2)-1)

        ! Auxiliary
        logical, dimension(size(hasvv, 1)-1, size(hasvv, 2)-1)  :: &
            issaddle, noval, hasval, c1, c2, c3, c4

        ! Initialize
        !===========
        ! Dimensions
        nv1 = size(hasvv, 1)
        nv2 = size(hasvv, 2)

        ! Initial values
        quadflags = 0 ! default: no value
        c1 = hasvv(1:nv1-1, 1:nv2-1)
        c2 = hasvv(2:nv1, 1:nv2-1)
        c3 = hasvv(2:nv1, 2:nv2) 
        c4 = hasvv(1:nv1-1, 2:nv2)
        
        ! Determine quad types
        issaddle = ( (c1 .and. c3) .and. ((.not. c2) .and. (.not. c4)) ) .or. &
            ( (c2 .and. c4) .and. ((.not. c1) .and. (.not. c3)) ) ! true only at opposite corners of a quad means saddle point
        noval = (c1 .and. c2 .and. c3 .and. c4) .or. &
            ((.not. c1) .and. (.not. c2) .and. (.not. c3) .and. (.not. c4)) ! all true or all false means no values
        hasval = .not. noval .and. .not. issaddle ! otherwise, there is just a value
        
        ! Set quadflags
        where (hasval)      quadflags = 1 
        where (issaddle)    quadflags = 4
        where (superquadflags > 0) quadflags = superquadflags

    end function 

    ! Coordinate/values getters/setters
    subroutine GetCoordinatesStructured2D(tracer, x, y)

        ! Declare variables
        !==================
        ! Arguments
        class(StructuredContourTracerUDT)   :: tracer 
        real(R8), allocatable, intent(out)  :: x(:), y(:)

        ! Get
        !====
        x = tracer%xg 
        y = tracer%yg

    end subroutine

    subroutine SetValuesStructured2D(tracer, v)

        ! Declare variables
        !==================
        ! Arguments
        class(StructuredContourTracerUDT)   :: tracer 
        real(R8), intent(in)                :: v(:)

        ! Set
        !====
        ! Check
        if (size(v) /= size(tracer%xg)) then 
            call gdErrorHandler('SetValuesStructured2D: value dimension '// &
                ' is inconsistent with current number of grid points')
        end if 

        ! Set
        tracer%V = reshape(v, [size(tracer%X), size(tracer%Y)])
        
    end subroutine

    

end module