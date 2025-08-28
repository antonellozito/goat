!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains all 2D interpolant types and implementations. 
! Each interpolant inherits from the generic 2D interpolant type, which 
! defines which routines etc should be present to evaluate. Depends on
! the 'precision' module. 

module Interpolant2D

    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_sparseinterface
    use mod_plotter
    use mod_structured2Dgridding

    !==================================================================!
    !                                                                  !
    !                            TYPES                                 !
    !                                                                  !
    !==================================================================!
    
    ! Generic interpolant type
    !=========================
    type, abstract :: GenericInterpolant2DUDT

        ! Description
        !============
        ! Generic 2D interpolant. Contains the following fields
        ! - xv, yv: coordinates at which the values are given
        ! - v:      the values at the coordinates
        ! Any other required fields are to be added by the specific
        ! interpolants themselves. 

        ! The following routines should be present in any 2D interpolant
        ! - Construct(xv, yv, v): interpolant constructor
        ! - Deconstruct()               : deconstructor
        ! - Evaluate(xq, yq, derivx, derivy) : evaluator, derivx and 
        !   derivy indicate the derivative order in x and y direction
        
        real(R8), allocatable   :: xv(:), yv(:), v(:)

    contains 

        ! Constructors (structured and unstructured)
        procedure :: ConstructStructuredGEN => ConstructFromStructuredData
        procedure :: ConstructUnstructuredGEN => ConstructFromUnstructuredData
        generic :: Construct => ConstructStructuredGEN, ConstructUnstructuredGEN

        procedure(ConstructSINT), deferred      :: ConstructStructured 
        procedure(ConstructUSINT), deferred     :: ConstructUnstructured  

        ! Evaluation
        procedure(EvaluateINT), deferred        :: Evaluate

        ! Visualization 
        procedure :: Visualize          => VisualizeInterpolant2D

    end type

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!
    ! Abstract interface
    abstract interface 
        
        ! Structured constructor
        subroutine ConstructSINT(interp, xg, yg, v)
            import :: GenericInterpolant2DUDT, R8
            class(GenericInterpolant2DUDT)      :: interp 
            real(R8), allocatable               :: xg(:), yg(:), v(:, :)   
        end subroutine  
        
        ! Unstructured constructor
        subroutine ConstructUSINT(interp, xg, yg, v)
            import :: GenericInterpolant2DUDT, R8
            class(GenericInterpolant2DUDT)      :: interp 
            real(R8), allocatable               :: xg(:), yg(:), v(:)   
        end subroutine  

        ! Evaluator
        subroutine EvaluateINT(interp, xq, yq, derivx, derivy, vq)
            import :: GenericInterpolant2DUDT, R8, I8
            class(GenericInterpolant2DUDT)      :: interp 
            real(R8), intent(in)                :: xq(:), yq(:)
            real(R8), intent(out)               :: vq(:)
            integer(I8), intent(in)             :: derivx, derivy
        end subroutine  
        
    end interface        
    
    contains 

    ! Structured constructor routine (wrapper)
    subroutine ConstructFromStructuredData(interp, xg, yg, v) 

        class(GenericInterpolant2DUDT)      :: interp 
        real(R8), allocatable               :: xg(:), yg(:), v(:, :)

        interp%xv = xg
        interp%yv = yg
        interp%v = reshape(v, [size(v)])
        call interp%ConstructStructured(xg, yg, v) 

    end subroutine 

    ! Unstructured constructor routine (wrapper)
    subroutine ConstructFromUnstructuredData(interp, xg, yg, v) 
        
        class(GenericInterpolant2DUDT)      :: interp 
        real(R8), allocatable               :: xg(:), yg(:), v(:)

        interp%xv = xg
        interp%yv = yg
        interp%v = v
        call interp%ConstructUnstructured(xg, yg, v)

    end subroutine 

    ! Visualization
    subroutine VisualizeInterpolant2D(interp, savefilepath, &
        nxin, nyin, xderivin, yderivin)

    ! Description
    !============
    ! This routine provides a crude way of visualizing the 
    ! levelset function by making a 2D patchplot on a structured 
    ! grid. The extent is equal to the coordinate boundaries 
    ! of the polygon set vertices, and the resolution in both 
    ! directions is fixed and constant. 

    ! Declare variables
    !==================
    ! Arguments
    class(GenericInterpolant2DUDT)           :: interp 
    character(*), intent(in)                 :: savefilepath
    integer(I8), intent(in), optional        :: nxin, nyin, xderivin, &
        yderivin 
    
    ! Auxiliary
    integer(I8)                         :: nx, ny, xderiv, yderiv 
    real(R8), allocatable               :: xgv(:), ygv(:), xg(:), &
        yg(:), vg(:)
    real(R8)                            :: xmin, ymin, xmax, ymax, &
        offsetx, offsety, dx, dy, dxgv, dygv

    ! Loop
    integer(I8)                         :: k

    ! Construct a 2D grid
    !====================
    ! Set mesh size
    if (present(nxin)) then 
        nx = nxin 
    else 
        nx = 200
    end if 
    if (present(nyin)) then 
        ny = nyin 
    else
        ny = 200
    end if 
    if (present(xderivin)) then 
        xderiv = xderivin 
    else
        xderiv = 0
    end if 
    if (present(yderivin)) then 
        yderiv = yderivin
    else
        yderiv = 0
    end if 

    ! Allocate
    allocate(xgv(nx), ygv(ny), xg(nx*ny), yg(nx*ny), &
        vg(nx*ny))

    ! Set the gridding vectors
    xmin = minval(interp%xv)
    xmax = maxval(interp%xv)
    ymin = minval(interp%yv)
    ymax = maxval(interp%yv)

    dx = (xmax - xmin) 
    dy = (ymax - ymin)
    
    offsetx = 0.0*dx 
    offsety = 0.0*dy 

    dxgv = (dx + 2*offsetx)/nx 
    dygv = (dy + 2*offsety)/ny 

    xgv(:) = dxgv*[(k, k = 0, nx-1)] - offsetx + xmin
    ygv(:) = dygv*[(k, k = 0, ny-1)] - offsety + ymin
    
    ! Construct
    call Construct2DStructuredGrid(xgv, ygv, nx, ny, xg, yg)

    ! Evaluate
    !=========
    ! Call evaluator
    call interp%Evaluate(xg, yg, xderiv, yderiv, vg)

    ! Write data
    !===========
    call Write3DCoordinateData(xg, yg, vg, savefilepath)

    ! Housekeeping
    !=============
    ! Deallocate
    deallocate(xgv, ygv, xg, yg, vg)

end subroutine
    

end module