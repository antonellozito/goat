!==================================================================!
!                                                                  !
!                        DOCUMENTATION                             !
!                                                                  !
!==================================================================!
! This module provides some functionality to trace 2D streamlines for a 
! field discretized on a certain grid. Currently, only 2D structured 
! non-uniform Cartesian grids are supported. Though streamlines are in
! general nothing more than 2D lines, we do return them as a type (array)
! to hedge for possible future extension. 

module mod_streamlinetracing2D

    ! Load modules
    use mod_precision
    use mod_errorhandler
    use mod_dynamicarrays
    use mod_structured2Dgridding

    implicit none
    private 
    public :: ConstructStructuredStreamlineTracer, StreamlineTracerUDT, &
        StructuredStreamlineTracerUDT, StreamlineUDT, &
        TraceStreamlinesStructured2D, CleanStreamlines

    ! Module parameters
    real(R8), parameter :: disttol = 1e-12 ! distance tolerance when cleaning streamlines [m]


    !==================================================================!
    !                                                                  !
    !                            TYPES                                 !
    !                                                                  !
    !==================================================================!

    ! Streamline type
    type :: StreamlineUDT 

        ! Description
        !============
        ! The streamline type currently only stores the streamline 
        ! coordinates x, y
        real(R8), allocatable   :: x(:), y(:)

    contains 

    end type

    ! General streamline tracer
    type, abstract :: StreamlineTracerUDT

        ! Description
        !============
        ! General type that traces streamlines. Allows abstraction from any
        ! underlying mechanics (structured/unstructured grid etc), except
        ! on the constructor level. All necessary data to trace the 
        ! streamlines should be saved in the object itself.
        ! For easier redefinition of tracing values, routines should be
        ! specified that return the location on which the field value
        ! should be known, and a routine that sets the values again. 
        
        ! Fields:
        !   step:       step size (relative) to be taken to trace contours
        !   nsteps:     maximal number of steps
        !   tol:        tolerance on the velocity magnitude to identify stagnation points
        
        real(R8)                :: step, tol 
        integer(I8)             :: nsteps

    contains 

        ! Tracer
        procedure(TraceStreamlinesINT), deferred    :: TraceStreamlines

        ! Coordinate getter
        procedure(GetTracerValueCoordinatesINT), deferred :: GetCoordinates 

        ! Field velocity setter
        procedure(SetTracerVelocityINT), deferred     :: SetVelocity

    end type

    ! Structured 2D tracer
    type, extends(StreamlineTracerUDT)     :: StructuredStreamlineTracerUDT

        ! Description
        !============
        ! Tracer for structured 2D grids. Values etc are all saved in 
        ! the structure, the tracer routine is a wrapper for the 
        ! standalone 2D tracer routine. 
        real(R8), allocatable       :: X(:), Y(:), U(:, :), V(:, :), &
            xg(:), yg(:)

    contains 

        ! Tracer
        procedure :: TraceStreamlines  => TraceStreamlinesStructured2DWrapper

        ! Coordinate getter
        procedure :: GetCoordinates => GetCoordinatesStructured2D

        ! Field value setter
        procedure :: SetVelocity      => SetVelocityStructured2D

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
        function TraceStreamlinesINT(tracer, x0, y0, xb, yb, direction) result(streamlines)
            
            import :: StreamlineTracerUDT, R8, StreamlineUDT, I8
            class(StreamlineTracerUDT)          :: tracer 
            real(R8), intent(in)                :: x0(:), y0(:) 
            real(R8), intent(inout)             :: xb(1:2), yb(1:2)
            integer(I8), intent(in)             :: direction(:)
            type(StreamlineUDT), allocatable    :: streamlines(:)

        end function

        ! Coordinates getter
        subroutine GetTracerValueCoordinatesINT(tracer, x, y) 
            
            import :: StreamlineTracerUDT, R8
            class(StreamlineTracerUDT)              :: tracer 
            real(R8), allocatable, intent(out)      :: x(:), y(:)

        end subroutine

        ! Values setter
        subroutine SetTracerVelocityINT(tracer, U, V) 
            
            import :: StreamlineTracerUDT, R8
            class(StreamlineTracerUDT)             :: tracer 
            real(R8), intent(in)                    :: U(:), V(:)

        end subroutine

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
    function ConstructStreamline(x, y) result(streamline)

        ! Description
        !============
        ! Streamline constructor that requires all streamline data to be set
        ! to be passed. Only for single streamlines

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)        :: x(:), y(:)
        type(StreamlineUDT)         :: streamline

        ! Construct
        !==========
        streamline%x           = x 
        streamline%y           = y 

    end function

    ! Structured tracer constructor
    function ConstructStructuredStreamlineTracer(U, V, X, Y) result(tracer)

        ! Description
        !============
        ! Construct a structured tracer. The result is returned as an 
        ! abstract type though that is allocated in this routine

        ! Declare variables
        !==================
        ! Arguments
        class(StreamlineTracerUDT), allocatable    :: tracer 
        real(R8), intent(in), dimension(:)      :: X, Y
        real(R8), intent(in)                    :: U(:, :), V(:, :)

        ! Auxiliary
        real(R8), allocatable                   :: xg(:), yg(:)


        ! Initialize
        !===========
        ! Allocate
        allocate(StructuredStreamlineTracerUDT::tracer) 

        ! Set values
        !===========
        select type (tracer)

        type is (StructuredStreamlineTracerUDT)

            tracer%U = U
            tracer%V = V 
            tracer%X = X 
            tracer%Y = Y 

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

    ! 2D structured streamline line tracer 
    subroutine TraceStreamlinesStructured2D(U, V, X, Y, x0, y0, &
        direction, xb, yb, s, tol, nmax, streamlines, isoutofbounds)

        ! Description
        !============
        ! This function traces a stream line (or stream lines if x0 and y0 are
        ! arrays) for a given velocity field (V, U) that is discretized on a
        ! structured 2D cartesian (non-uniform) mesh. X and Y should be the
        ! coordinate vectors that construct the 2D cartesian mesh such that V(i, j)
        ! represents the x-component of the velocity at point X(i), Y(j) (same for
        ! U). Stream lines are traced in both directions until either the point
        ! would exceed the min/max values of the coordinate vectors X, Y, or the
        ! given bounds on x, y (xb, yb). The step size has to be given as argument
        ! 's' (a good value is e.g. 0.5)
        
        ! Notes
        !------
        ! Note 1: to avoid scale effects due to the magnitude of the velocity
        ! field, we rescale the field at each point to have unit length. The given
        ! step length is then scaled with respect to the domain extent (we take
        ! min( (X(end)-X(1))/numel(X), (Y(end)-Y(1))/numel(Y)) as the domain 
        ! extent). A step length of 1 would then correspond to the average width of
        ! one cell in the domain in X or Y direction. 
        
        ! Note 2: to deal with stagnation points, we do check the magnitude of the
        ! velocity. If this magnitude is below some predefined tolerance 'tol', the
        ! algorithm will stop. 
        
        ! Note 3: since generally speaking streamlines may become infinite (e.g. a
        ! closed contour), we stop after nmax number of nodes. If this number isn't
        ! given, a default of 10 000 nodes is used. 
        
        ! Algorithm
        !----------
        ! General
        ! 0) Do some sanity checks:
        !   - is the initial point inside the domain? 
        !   - does the given step length make sense? 
        !   - what are the limiting bounds on x, y? 
        ! 1) Scale the velocity field to (local) unity vectors
        ! 2) For each point, trace the 'forward' and 'backward' contour until the
        ! bounds are reached
        ! 3) Concatenate and return
        
        ! Streamline tracing
        ! 1) Check if the current point lies outside of the bounds. If this is the
        ! case, exit. Otherwise, go to 2)
        ! 2) Check in which cell the current point is located. 
        ! 3) Compute the velocity at this point by bilinear interpolation
        ! 4) Update the current position by moving the point a certain length 's'
        ! along the current velocity vector. Go to 1)

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)                :: U(:, :), V(:, :), X(:), &
            Y(:), x0(:), y0(:), s, tol
        real(R8), intent(inout)             :: xb(1:2), yb(1:2)
        integer(I8), intent(in)             :: direction(:), nmax
        type(StreamlineUDT), allocatable, intent(out)   :: streamlines(:)
        logical, allocatable, intent(out)               :: isoutofbounds(:)

        ! Auxiliary
        integer(I8)                             :: nx, ny, nv
        real(R8), allocatable, dimension(:)     :: xfw, yfw, xbw, ybw
        real(R8)                                :: emptyarray(0), step

        ! Loop
        integer(I8)                             :: i 
        
        ! Initialize
        !-----------
        ! Check sizes
        nx = size(X)
        ny = size(Y)
        nv = size(x0)

        ! Checks
        if (nv /= size(y0)) then 
            call gdErrorHandler('TraceStreamlineStructured2D: ' // & 
                'start coordinates for tracing have incompatible dimensions')
        end if
        if ((size(U, 1) /= nx) .or. (size(U, 2) /= ny)) then 
            call gdErrorHandler('TraceStreamlineStructured2D: ' // & 
                'incompatible size of velocity x-component with given grid coordinates')
        end if 
        if (nv /= size(direction)) then 
            call gdErrorHandler('TraceStreamlineStructured2D: ' // & 
                'incompatible size of direction array w.r.t start coordinates')
        end if 
        if ((size(V, 1) /= nx) .or. (size(V, 2) /= ny)) then 
            call gdErrorHandler('TraceStreamlineStructured2D: ' // & 
                'incompatible size of velocity y-component with given grid coordinates')
        end if 

        ! Check limiting bounds
        xb = [max(minval(X), minval(xb)), min(maxval(X), maxval(xb))]
        yb = [max(minval(Y), minval(yb)), min(maxval(Y), maxval(yb))]
        
        ! Rescale step size
        if (abs(X(nx)-X(1)) < abs(Y(ny)-Y(1))) then 
            step = s/(nx*abs(X(nx)-X(1)))
        else
            step = s/(ny*abs(Y(ny)-Y(1)))
        end if 
        
        ! Check points
        isoutofbounds = (x0 < xb(1)) .or. (x0 > xb(2)) .or. (y0 < yb(1)) .or. (y0 > yb(2))
        
        ! Initialize
        allocate(streamlines(nv))
        
        ! Loop and trace
        !---------------
        do i = 1, nv
            ! Check if we should trace
            if (isoutofbounds(i)) then 
                ! Simply initialize empty
                streamlines(i) = ConstructStreamline(emptyarray, emptyarray)
                cycle 
            end if 

            select case (direction(i))

            case (0)

                ! Both directions
                ! Trace forward
                call TraceStreamlineSegment(U, V, X, Y, x0(i), y0(i), xb, yb, s, tol, nmax, xfw, yfw)
                    
                ! Trace backward
                call TraceStreamlineSegment(-U, -V, X, Y, x0(i), y0(i), xb, yb, s, tol, nmax, xbw, ybw)
                
                ! Concatenate
                streamlines(i) = ConstructStreamline([xfw(size(xfw):2:-1), xbw], &
                    [yfw(size(yfw):2:-1), ybw])

            case (1)

                ! Forward
                call TraceStreamlineSegment(U, V, X, Y, x0(i), y0(i), xb, yb, s, tol, nmax, xfw, yfw)
                streamlines(i) = ConstructStreamline(xfw, yfw)

            case (-1)

                ! Backward
                call TraceStreamlineSegment(-U, -V, X, Y, x0(i), y0(i), xb, yb, s, tol, nmax, xbw, ybw)
                streamlines(i) = ConstructStreamline(xbw, ybw)

            case default

                ! Unknown
                call gdErrorHandler('TraceStreamlineStructured2D: ' // & 
                    'unknown tracing direction')
            end select

        end do
        
        
    end subroutine 
        
    ! Wrapper for structured contour line tracer
    function TraceStreamlinesStructured2DWrapper(tracer, x0, y0, &
        xb, yb, direction) result(streamlines)

        ! Description
        !============
        ! Wrapper for the 2D tracer. The order etc is added to the 
        ! tracer routine

        ! Declare variables
        !==================
        ! Arguments
        class(StructuredStreamlineTracerUDT)            :: tracer 
        real(R8), intent(in)                            :: x0(:), y0(:)
        real(R8), intent(inout)                         :: xb(1:2), yb(1:2)
        integer(I8), intent(in)                         :: direction(:)
        type(StreamlineUDT), allocatable                :: streamlines(:)

        ! Auxiliary
        logical, allocatable, dimension(:)              :: isoutofbounds
        
        ! Trace
        !======
        call TraceStreamlinesStructured2D(tracer%U, tracer%V, tracer%X, &
            tracer%Y, x0, y0, direction, xb, yb, tracer%step, tracer%tol, tracer%nsteps, &
            streamlines, isoutofbounds)

    end function 

    !------------------------------------------------------------------!
    !                             AUXILIARY                            !
    !------------------------------------------------------------------!

    ! Streamline tracer in one direction
    subroutine TraceStreamlineSegment(U, V, X, Y, x0, y0, xb, yb, s, tol, nmax, xs, ys)

        ! Description
        !============
        ! Trace a streamline in a specific direction

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)        :: U(:, :), V(:, :), x0, y0, &
            xb(1:2), yb(1:2), s, tol, X(:), Y(:)
        integer(I8), intent(in)     :: nmax
        real(R8), intent(out), allocatable  :: xs(:), ys(:)

        ! Auxiliary
        integer(I8)                         :: initdim, i, j
        real(R8)                            :: dx, dxold, dy, dyold, &
            tv1, tv2, U0, V0, vmag, xt, yt, dotprod

        ! Loop
        integer(I8)                         :: k 
        
        
        ! Initialize
        initdim = ceiling(1_R8/s)
        allocate(xs(initdim), ys(initdim))
        
        ! Point counter
        k = 1
        dx = 0
        dy = 0
        xt = x0
        yt = y0
        do while (.true.)
            ! Check if we need to increment
            if (k > size(xs)) then 
                xs = [xs, spread(0.0_R8, 1, initdim)]
                ys = [ys, spread(0.0_R8, 1, initdim)]
            end if
            
            ! Add current point
            xs(k) = xt
            ys(k) = yt
            k = k + 1
            dxold = dx
            dyold = dy
            
            ! Check if we should exit
            if ((x0 < xb(1)) .or. (x0 > xb(2)) .or. (y0 < yb(1)) .or. (y0 > yb(2))) then 
                exit
            end if 
            
            ! Check in which box we are (if we know we stepped out)
            i = findloc(x0 >= X, .true., 1, back=.true.)
            j = findloc(y0 >= Y, .true., 1, back=.true.)
            
            if (i == size(X)) then 
                i = i - 1
            end if 
            if (j == size(Y)) then 
                j = j - 1
            end if
            
            ! Compute velocity field by bilinear interpolation
            tv1 = (U(i+1, j) - U(i, j))/(X(i+1) - X(i))*(xt - X(i)) + U(i, j)
            tv2 = (U(i+1, j+1) - U(i, j+1))/(X(i+1) - X(i))*(xt - X(i)) + U(i, j+1)
            U0 = (tv2 - tv1)/(Y(j+1) - Y(j))*(yt - Y(j)) + tv1
            
            tv1 = (V(i+1, j) - V(i, j))/(X(i+1) - X(i))*(xt - X(i)) + V(i, j)
            tv2 = (V(i+1, j+1) - V(i, j+1))/(X(i+1) - X(i))*(xt - X(i)) + V(i, j+1)
            V0 = (tv2 - tv1)/(Y(j+1) - Y(j))*(yt - Y(j)) + tv1
            
            ! Compute velocity magnitude
            vmag = sqrt(U0**2 + V0**2);
            
            ! Check for termination
            if (vmag < tol) then 
                exit
            end if
            
            ! Rescale
            U0 = U0/vmag
            V0 = V0/vmag
            
            ! Compute step
            dx = U0*s
            dy = V0*s
            
            ! Update
            xt = xt + dx
            yt = yt + dy
            
            ! Check if we should exit based on new segment
            if (k > 2) then 
                dotprod = dxold*dx + dyold*dy
                if (dotprod < 0) then 
                    exit 
                end if 
            end if 
            if (k >= nmax) then 
                ! Maximum number of nodes reached
                exit 
            end if 
            
        end do 
        
        ! Trim
        xs = xs(1:k-1)
        ys = ys(1:k-1)
            
    end subroutine

    ! Coordinate/values getters/setters
    subroutine GetCoordinatesStructured2D(tracer, x, y)

        ! Declare variables
        !==================
        ! Arguments
        class(StructuredStreamlineTracerUDT)   :: tracer 
        real(R8), allocatable, intent(out)  :: x(:), y(:)

        ! Get
        !====
        x = tracer%xg 
        y = tracer%yg

    end subroutine

    subroutine SetVelocityStructured2D(tracer, U, V)

        ! Declare variables
        !==================
        ! Arguments
        class(StructuredStreamlineTracerUDT)   :: tracer 
        real(R8), intent(in)                :: U(:), V(:)

        ! Set
        !====
        ! Check
        if ((size(V) /= size(tracer%xg)) .or. (size(U) /= size(tracer%xg))) then 
            call gdErrorHandler('SetValuesStructured2D: value dimension '// &
                ' is inconsistent with current number of grid points')
        end if 

        ! Set
        tracer%U = reshape(U, [size(tracer%X), size(tracer%Y)])
        tracer%V = reshape(V, [size(tracer%X), size(tracer%Y)])
        
    end subroutine

    ! Streamline clean-up
    subroutine CleanStreamlines(streamlines)

        ! Description
        !============
        ! This routine cleans up the streamlines, i.e. it removes subsequent
        ! points that are up to disttol coinciding. May be necessary
        ! for later intersection computing etc. 

        ! Declare variables
        !==================
        ! Arguments
        class(StreamlineUDT), intent(inout)        :: streamlines(:)

        ! Auxiliary
        logical, allocatable                    :: delind(:)
        real(R8), allocatable                   :: dx(:), dy(:)

        ! Loop
        integer(I8)                             :: i 
        
        ! Clean
        !======
        do  i = 1, size(streamlines)
            dx = streamlines(i)%x(2:size(streamlines(i)%x)) - &
                streamlines(i)%x(1:size(streamlines(i)%x)-1)
            dy = streamlines(i)%y(2:size(streamlines(i)%y)) - &
                streamlines(i)%y(1:size(streamlines(i)%y)-1)
            allocate(delind(size(dx)))
            delind = (abs(dx) <= disttol) .and. (abs(dy) <= disttol)
            if (any(delind)) then 
                streamlines(i)%x = pack(streamlines(i)%x, [.true., .not. delind])
                streamlines(i)%y = pack(streamlines(i)%y, [.true., .not. delind])
            end if 
            deallocate(delind)
        end do 
    end subroutine 

end module