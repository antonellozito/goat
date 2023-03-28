!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module provides the PolygonShapeFunctionUDT type, which is an 
! object that provides a function that is zero on the polygon itself
! and strictly positive when not on the polygon. Though non-smooth, it
! yields decent results in combination with the grid deformation 
! framework.  Basically, we create a function of which the contours are
! based on the polygon, and where the zero-contour corresponds to the 
! polygon itself. 

! Mathematical description
!=========================
! The mathematical formulation is based on two key requirements for a 
! point to strictly lie on a polygon edge: 
!
! - The perpendicular distance between the line formed by the edge 
! points and the query point should be zero
! - The lengths of the vectors between the edge points and the query 
! point should equal the length of the vector between the edge points
! 
! The function value can therefore be defined as (for one edge, one 
! sample point):
!
!       f = abs (abs(vxp1) + abs(vxp2) - abs(txp)) 
!           + abs (abs(vyp1) + abs(vyp2) - abs(typ))
!           + abs(dvn)
! 
! Here, vp1, vp2 are the vectors between the query point and the first
! and second edge points, tp is the vector between the two edge points
! (and their components are indicated with 'x' and 'y'). Dvn is the 
! perpendicular length between the query point and the line passing 
! through the edge points, which is computed by projecting the vector
! between the edge center (xf, yf) and query point with the edge normal:
!
!       dvn = vx * nxp + vy * nyp
!
! where vx, vy are the components of the vertex - edge face vector. 

! For the case with multiple edges, we simply loop over all edges and
! compute the minimal value of this function for each of the edges. 
! Though the abs function only results in non-smoothness at the edge,
! the max function leads to non-smoothness at locations away from the 
! edge.  

! Notes 
!======
! Note 1: where f = 0 exactly, the gradient is undefined due to the abs
! function. 


module PolygonShapeFunction 

    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_plotter

    ! The usual
    implicit none
    save
    private 

    ! Types
    public PolygonShapeFunctionUDT

    !==================================================================!
    !                                                                  !
    !                              TYPES                               !
    !                                                                  !
    !==================================================================!

    ! Main type
    !==========
    type PolygonShapeFunctionUDT 

        ! Description
        !============
        ! Main polygon type that contains the necessary functions to
        ! evaluate the polygon and its derivatives. It has to be 
        ! initialized using the polygon coordinates, at which step the
        ! necessary auxiliary data is saved in the object. The following
        ! fields are used:
        !
        ! - nc:         number of polygon edgees
        ! - x, y        nc-by-2 array, polygon egde coordinates ( first
        !               column is first node, second column is second 
        !               node of the polygon) 
        ! - xf, yf      nc-by-1 array with the polygon edge centers
        ! - nxp, nyp    the normalized polygon edge normals
        ! - txp, typ    the NOT normalized polygon tangents
        !
        ! The following routines are provided: 
        !
        ! - Evaluate    Evaluate the function value or its derivatives, 
        !               depending on the passed arguments. See function
        !               implementation for more details. 
        ! - Initialize  Initialize the shape function, based on the 
        !               passed polygon. Can also be used to reconstruct
        !               the representation in case the polygon would
        !               change. 
        ! - Visualize   Visualize the values of the polygon shape 
        !               function

        ! Declare
        integer(I8)                 :: nsp, nc
        real(R8), allocatable       :: x(:, :), y(:, :), nxp(:), &
            nyp(:), txp(:), typ(:), xf(:), yf(:)
        integer(I8), allocatable    :: sp(:)

    
    contains 
        ! Constructor
        procedure :: Initialize     => InitializePolygonShapeFunction

        ! Evaluator
        procedure :: Evaluate       => EvaluatePolygonShapeFunction

        ! Visualization
        procedure :: Visualize      => VisualizePolygonShapeFunction

        ! Destructor
        final :: DestroyPolygonShapeFunction

        ! Allocation
        procedure :: Allocate       => AllocatePolygonShapeFunction

        ! Deallocation
        procedure :: Deallocate     => DeallocatePolygonShapeFunction

    end type

    contains 

    !==================================================================!
    !                                                                  !
    !                            ROUTINES                              !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                       POLYGON SHAPE FUNCTIONS                    !
    !------------------------------------------------------------------!
    ! Constructor
    subroutine InitializePolygonShapeFunction(psf, x, y, nc)

        ! Description
        !============
        ! This function initializes the polygon shape function given:
        ! - the coordinates x, y of the polygon (and the number nc)
        ! - the edges of the polygon ep (and their number nep)
        ! - the vector with polygon starts sp (and their number nsp)
        ! These quantities allow to define the shape function for 
        ! arbitrary, multiple-segment, partially open/closed, polygons.
        ! This should be the most general case. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonShapeFunctionUDT)      :: psf 
        integer(I8), intent(in)             :: nc
        real(R8), dimension(1:nc, 2)        :: x, y 

        ! Auxiliary
        real(R8), allocatable               :: nxp(:), nyp(:), txp(:), &
            typ(:), xf(:), yf(:), nn(:)

        ! Loop

        ! Data

        ! Initialize
        !===========
        ! Set some values
        psf%nc = nc

        ! Allocate
        call psf%Allocate()

        ! Add
        psf%x = x 
        psf%y = y

        ! Compute auxiliary quantities
        !=============================
        ! Allocate
        allocate(nxp(nc), nyp(nc), txp(nc), typ(nc), xf(nc), yf(nc), &
            nn(nc) )
        
        ! Polygon normals - normalized
        nxp = -( y(:, 2) - y(:, 1) ) 
        nyp = x(:, 2) - x(:, 1)
        nn = (nxp**2 + nyp**2)**0.5
        nxp = nxp/nn 
        nyp = nyp/nn 

        ! Polygon tangents - NOT normalized
        txp = x(:, 2) - x(:, 1)
        typ = y(:, 2) - y(:, 1)

        ! Polygon edge centers
        xf = 0.5*(sum(x, 2))
        yf = 0.5*(sum(y, 2))

        ! Add to shape function
        psf%nxp = nxp 
        psf%nyp = nyp 
        psf%txp = txp 
        psf%typ = typ 
        psf%xf = xf
        psf%yf = yf

         ! Deallocate
        deallocate(nxp, nyp, txp, typ, xf, yf)


    end subroutine

    ! Evaluation
    subroutine EvaluatePolygonShapeFunction(psf, xq, yq, vq, &
        derivx, derivy)

        ! Description
        !============
        ! Evaluate the polygon shape function in the query points xq, yq
        ! (value returned is in vq). The 'derivx' and 'derivy' input 
        ! arguments are used to specify the order of the derivative to 
        ! return. 

        ! Note that derivates are returned pointwise as there is no 
        ! dependence of the values on other sample points or polygon
        ! coordinates. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonShapeFunctionUDT)          :: psf 
        real(R8), intent(in)                    :: xq(:), yq(:)
        real(R8), intent(out)                   :: vq(:) 
        character(len=1), intent(in)            :: derivx, derivy 

        ! Auxiliary
        character(len=2)                        :: deriv 
        integer(I8)                             :: nq, nc
        real(R8), allocatable                   :: vnx(:), vny(:), &
            vtx(:, :), vty(:, :), f(:), ind(:), fq(:), vtxq(:, :), &
            vtyq(:, :), dvnq(:), myones(:), vnxq(:), vnyq(:), dvn(:), &
            nxpq(:), nypq(:), txpq(:), typq(:)
        real(R8)                                :: macheps = 0

        ! Loop
        integer(I8)                             :: i, thisloc

        ! Data

        ! Initialize
        !===========
        ! Set macheps
        macheps = 1e-13

        ! Get the number of points to evaluate on, assume proper
        ! dimensions of inputs
        nq = size(xq) 
        nc = psf%nc

        ! Check the derivative, compute
        deriv = derivx // derivy

        ! Allocate
        allocate(vnx(nc), vny(nc), vtx(nc, 2), vty(nc, 2), f(nc), &
            dvn(nc), ind(nc), fq(nq), vtxq(nq, 2), vtyq(nq, 2), &
            dvnq(nq), myones(nq), vnxq(nq), vnyq(nq), nxpq(nq), &
            nypq(nq), txpq(nq), typq(nq))
        
        ! Initialize
        myones(:) = 1

        ! Associate
        associate(&
            x       => psf%x,   &
            y       => psf%y,   &
            xf      => psf%xf,  &
            yf      => psf%yf,  &
            nxp     => psf%nxp, &
            nyp     => psf%nyp, &
            txp     => psf%txp, &
            typ     => psf%typ)

        ! Compute
        !========
        ! Loop over all points
        do i = 1, nq 

            ! Normal distance
            !----------------
            ! Compute vectors between query point and edge center
            vnx(:) = -(xf(:) - xq(i)) 
            vny(:) = -(yf(:) - yq(i))

            ! Project 
            dvn = vnx * nxp + vny * nyp 

            ! Tangential distances
            !---------------------
            ! Compute vectors between query point and vertices
            vtx(:,:) = x(:,:) - xq(i)
            vty(:,:) = y(:,:) - yq(i)
            
            ! Compute 
            !--------
            ! Value
            f(:) = abs( abs(vtx(:,1)) + abs(vtx(:, 2)) - abs(txp) ) &
                + abs( abs(vty(:,1)) + abs(vty(:, 2)) - abs(typ) ) &
                + abs(dvn)

            ! Maximum location
            thisloc = minloc(f, 1)

            ! Store variables for later
            dvnq(i) = dvn(thisloc)
            fq(i) = f(thisloc)
            vnxq(i) = vnx(thisloc)
            vnyq(i) = vny(thisloc)
            vtxq(i, :) = vtx(thisloc, :)
            vtyq(i, :) = vty(thisloc, :)
            txpq(i) = txp(thisloc)
            typq(i) = typ(thisloc)
            nxpq(i) = nxp(thisloc)
            nypq(i) = nyp(thisloc)

        end do

        ! Hedge for machine precision
        where (fq .lt. macheps) fq = 0

        ! Return
        !=======
        ! Check what to compute exactly
        select case(deriv)

        case ('00')

            ! Function evaluation
            vq = fq

        case ('10')

            ! fx
            vq = nxpq*sign(myones, dvnq) - sign(myones, abs(vtxq(:, 1)) &
                - abs(txpq) + abs(vtxq(:, 2))) * &
                (sign(myones, vtxq(:, 1)) + sign(myones, vtxq(:, 2)))

            ! Set the gradient to 1 - can choose this value, since the 
            ! function is exactly zero here anyway. 
            where (fq .eq. 0) vq = 1

        case ('01')

            ! fy
            vq = nypq*sign(myones, dvnq) - sign(myones, abs(vtyq(:, 1)) &
                - abs(typq) + abs(vtyq(:, 2))) * &
                (sign(myones, vtyq(:, 1)) + sign(myones, vtyq(:, 2)))

            ! Set the gradient to 1 - can choose this value, since the 
            ! function is exactly zero here anyway. 
            where (fq .eq. 0) vq = 1

        case ('20')

            ! fxx - simply zero
            vq(:) = 0

        case ('11')

            ! fxy - simply zero
            vq(:) = 0

        case ('02')

            ! fyy - simply zero
            vq(:) = 0

        case default

            ! Throw error
            call gdErrorHandler('EvaluatePolygonShapeFunction: ' & 
                // 'derivative case ' // deriv // ' not implemented')

        end select

        ! End associate
        end associate

        ! Deallocate
        deallocate(vnx, vny, vtx, vty, f, ind, fq, vtxq, vtyq, dvnq)

    end subroutine

    ! Visualization
    subroutine VisualizePolygonShapeFunction(psf)

        ! Description
        !============
        ! Visualize the polygon shape function using the routines 
        ! provided in the mod_plotter routines. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonShapeFunctionUDT)      :: psf 
        
        ! Auxiliary
        integer(I8)                         :: nx, ny
        real(R8), allocatable               :: xgv(:), ygv(:), xg(:), &
            yg(:), vg(:), vsg(:, :)
        real(R8)                            :: xmin, ymin, xmax, ymax, &
            offsetx, offsety, dx, dy, dxgv, dygv

        ! Loop
        integer(I8)                         :: k

        ! Construct a 2D grid
        !====================
        ! Set mesh size
        nx = 100
        ny = 100

        ! Allocate
        allocate(xgv(nx), ygv(ny), xg(nx*ny), yg(nx*ny), vsg(nx, ny), &
            vg(nx*ny))

        ! Set the gridding vectors
        xmin = minval(psf%xf)
        xmax = maxval(psf%xf)
        ymin = minval(psf%yf)
        ymax = maxval(psf%yf)

        dx = (xmax - xmin) 
        dy = (ymax - ymin)
        
        offsetx = 0.05*dx 
        offsety = 0.05*dy 

        dxgv = (dx + 2*offsetx)/nx 
        dygv = (dy + 2*offsety)/ny 

        xgv(:) = dxgv*[(k, k = 0, nx-1)] - offsetx + xmin
        ygv(:) = dygv*[(k, k = 0, ny-1)] - offsety + ymin
        
        ! Construct
        call Construct2DStructuredGrid(xgv, ygv, nx, ny, xg, yg)

        ! Evaluate
        !=========
        ! Call evaluator
        call psf%Evaluate(xg, yg, vg, '0', '0')

        ! Reshape
        vsg = reshape(vg, (/nx, ny/))

        ! Plot
        !=====
        call Plot2DStructuredField(vsg, xgv, ygv, nx, ny, '-p')

        ! Housekeeping
        !=============
        ! Deallocate
        deallocate(xgv, ygv, xg, yg, vsg, vg)


    end subroutine

    ! Allocation
    subroutine AllocatePolygonShapeFunction(psf)

        ! Description
        !============
        ! Allocate

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonShapeFunctionUDT)            :: psf

        ! Allocate
        !=========
        ! Coordinate arrays
        allocate(psf%x(psf%nc, 2))
        allocate(psf%y(psf%nc, 2))

        ! Polygon normals
        allocate(psf%nxp(psf%nc), psf%nyp(psf%nc))
        
        ! Polygon tangents
        allocate(psf%txp(psf%nc), psf%typ(psf%nc))

        ! Polygon edge centers
        allocate(psf%xf(psf%nc), psf%yf(psf%nc))
        
    end subroutine

    ! Deallocate
    subroutine DeallocatePolygonShapeFunction(psf)

        ! Description
        !============
        ! Deallocate

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonShapeFunctionUDT)            :: psf

        ! Deallocate
        !===========
        if (allocated(psf%x)) then ! assume rest allocated as well
            deallocate(psf%x, psf%y, psf%nxp, psf%nyp, psf%txp, psf%typ, &
                psf%xf, psf%yf)
        end if

    end subroutine

    ! Destructor
    subroutine DestroyPolygonShapeFunction(psf)

        ! Description
        !============
        ! Destructor

        ! Declare variables
        !==================
        ! Arguments
        type(PolygonShapeFunctionUDT)       :: psf 

        ! Destroy
        !========
        call psf%Deallocate()

    end subroutine


end module
