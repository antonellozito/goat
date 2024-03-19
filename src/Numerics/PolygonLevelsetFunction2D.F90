!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module provides the (abstract) PolygonLevelsetFunction2D type 
! along with other (non-abstract) derived types. The main object 
! provides (at least) evaluation routines for (x, y) coordinate pairs. 
! Each levelset is based on a polygon set object, where the zero level 
! set corresponds to the (approximation of) the polygon set. Note that 
! most levelset functions here are only defined for closed, 
! non-intersecting (possibly nested) sets of polygons. 

! For mathematical descriptions, see the definition of each extended
! type. Basically, there are three types: 
!
!   1)  the 'general' one, which is valid for any polygon 
!       (self-intersecting, closed/open, ...) but which 
!       is discontinuous. The description is exact. 
!   2)  a smoother (but still non-smooth in several regions) 
!       approximation for nested closed non-intersecting polygons only. 
!       Still an exact representation of the original polygon set
!   3)  a fully smooth yet approximate representation of 2) by using 
!       an interpolation function (hence this module requires the 
!       Interpolant module)

module PolygonLevelsetFunction2D 

    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_constants
    use mod_plotter
    use Interpolant
    use mod_polygon

    ! The usual
    implicit none
    save
    private 

    ! Types
    public PolygonLevelsetFunction2DUDT, PolygonLevelsetFunction2DGeneralUDT, &
        PolygonLevelsetFunction2DClosedExactUDT, &
        PolygonLevelsetFunction2DClosedApproximationUDT
    public PLF2DOptionsUDT, PLF2DGeneralOptionsUDT, &
        PLF2DClosedExactOptionsUDT, PLF2DClosedApproximationOptionsUDT

    ! Routines
    public InitializePolygonLevelsetFunction2D

    !==================================================================!
    !                                                                  !
    !                              TYPES                               !
    !                                                                  !
    !==================================================================! 

    ! Abstract
    !=========
    ! Options type
    type, abstract :: PLF2DOptionsUDT 

        ! Description
        !============
        ! Abstract options type to get initialization routine to work (?)
        character(:), allocatable   :: levelsettype 

    end type

    ! Polygon levelset function type
    type, abstract :: PolygonLevelsetFunction2DUDT 

        ! Description
        !============
        ! Main abstract type

        ! Fields
        type(PolygonSetUDT)     :: ps 

    contains

        ! Initialization
        procedure(InitializePLF2DINT), deferred     :: Initialize

        ! Evaluation
        procedure(EvaluatePLF2DINT), deferred       :: Evaluate 

        ! Visualization
        procedure   :: Visualize        => VisualizePolygonLevelsetFunction2D

    end type

    ! Extended
    !=========
    ! Options for general levelset function
    type, extends(PLF2DOptionsUDT)  :: PLF2DGeneralOptionsUDT

    contains 

        ! Options setter
        procedure :: Set    => SetOptionsGeneral

    end type

    ! Options for closed exact representation
    type, extends(PLF2DOptionsUDT)  :: PLF2DClosedExactOptionsUDT
    contains 

        ! Options setter
        procedure :: Set    => SetOptionsClosedExact

    end type

    ! Options for closed approximate representation
    type, extends(PLF2DOptionsUDT)  :: PLF2DClosedApproximationOptionsUDT

        ! Fields
        ! Options for closed exact representation (none for now)
        type(PLF2DClosedExactOptionsUDT) :: optionsClosedExact

        ! Interpolation options
        real(R8)        :: offsetx, offsety 
        integer(I8)     :: resx, resy, C, M 
        character(:), allocatable   :: meth 

    contains 

        ! Options setter
        procedure :: Set    => SetOptionsClosedApproximation

    end type
    
    ! General levelset function
    type, extends(PolygonLevelsetFunction2DUDT)  :: &
        PolygonLevelsetFunction2DGeneralUDT 

        ! Description
        !============
        ! General levelset type. The following fields are added for 
        ! easier computation later on (computed based on given polygon
        ! set):
        !
        ! - nc:         number of polygon edgees
        ! - x, y        nc-by-2 array, polygon egde coordinates ( first
        !               column is first node, second column is second 
        !               node of the polygon) 
        ! - xf, yf      nc-by-1 array with the polygon edge centers
        ! - nxp, nyp    the normalized polygon edge normals
        ! - txp, typ    the NOT normalized polygon tangents

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

        ! Fields
        type(PLF2DGeneralOptionsUDT)    :: options
        integer(I8)                     :: nc
        real(R8), allocatable           :: xf(:), yf(:), nxp(:), &
            nyp(:), txp(:), typ(:), x(:, :), y(:, :)

    contains 

        ! Initializer
        procedure :: Initialize     => InitializePLF2DGeneral

        ! Evaluate
        procedure :: Evaluate       => EvaluatePLF2DGeneral

    end type

    ! Closed exact representation
    type, extends(PolygonLevelsetFunction2DUDT)  :: &
        PolygonLevelsetFunction2DClosedExactUDT 

        ! Description
        !============
        ! Exact representation of closed, non-intersecting set of 
        ! polygons. Hereto, we subdivide the 2D domain into vertex and
        ! edge regions, which then determines how the distance to the 
        ! polygon is computed (either the Euclidean distance to the 
        ! vertex, or the normal distance to the line formed by a 
        ! certain polygon edge). 

        ! Fields
        type(PLF2DClosedExactOptionsUDT)    :: options
        real(R8), allocatable               :: xp(:), yp(:), xf(:), &
            yf(:), nxp(:), nyp(:), tnp(:), nxpv(:, :), &
            nypv(:, :), crossprod(:), thetav(:)

    contains 

        ! Initialization
        procedure :: Initialize     => InitializePLF2DClosedExact

        ! Evaluate
        procedure :: Evaluate       => EvaluatePLF2DClosedExact

    end type

    ! Closed approximate representation
    type, extends(PolygonLevelsetFunction2DUDT)  :: &
        PolygonLevelsetFunction2DClosedApproximationUDT 

        ! Options
        type(PLF2DClosedApproximationOptionsUDT)    :: options

        ! Interpolant
        type(StructuredInterpolant2DUDT)            :: interp

    contains 

        ! Initialization
        procedure :: Initialize     => InitializePLF2DClosedApproximation

        ! Evaluation
        procedure :: Evaluate       => EvaluatePLF2DClosedApproximation

    end type

    
    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================! 

    ! Abstract
    !=========
    abstract interface 

        ! Initialization
        subroutine InitializePLF2DINT(plf, ps)
            import :: PolygonLevelsetFunction2DUDT, PolygonSetUDT, &
                PLF2DOptionsUDT
            class(PolygonLevelsetFunction2DUDT) :: plf 
            type(PolygonSetUDT), intent(in)     :: ps
        end subroutine

        ! Evaluation
        subroutine EvaluatePLF2DINT(plf, xq, yq, derivx, derivy, vq)
            import :: PolygonLevelsetFunction2DUDT, R8, I8
            class(PolygonLevelsetFunction2DUDT)     :: plf 
            real(R8), intent(in)                    :: xq(:), yq(:)
            real(R8), intent(out)                   :: vq(size(xq))
            integer(I8), intent(in)                 :: derivx, derivy

        end subroutine

    end interface
    
    contains

    !==================================================================!
    !                                                                  !
    !                            ROUTINES                              !
    !                                                                  !
    !==================================================================! 

    !------------------------------------------------------------------!
    !                       ABSTRACT TYPE ROUTINES                     !
    !------------------------------------------------------------------!

    ! Initializer
    subroutine InitializePolygonLevelsetFunction2D(plf, ps, options) 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonLevelsetFunction2DUDT), allocatable        :: plf 
        type(PolygonSetUDT), intent(in)                         :: ps
        class(PLF2DOptionsUDT), intent(in)                      :: options

        ! Initialize
        !===========
        select type (options)

        type is (PLF2DGeneralOptionsUDT)

            ! Allocate
            allocate(PolygonLevelsetFunction2DGeneralUDT::plf)

            ! Initialize 
            select type(plf)

            type is (PolygonLevelsetFunction2DGeneralUDT)

                plf%options = options
                call plf%Initialize(ps)

            end select

        type is (PLF2DClosedExactOptionsUDT) 

            ! Allocate
            allocate(PolygonLevelsetFunction2DClosedExactUDT::plf)

            ! Initialize 
            select type(plf)

            type is (PolygonLevelsetFunction2DClosedExactUDT)

                plf%options = options
                call plf%Initialize(ps)

            end select

        type is (PLF2DClosedApproximationOptionsUDT)

            ! Allocate
            allocate(PolygonLevelsetFunction2DClosedApproximationUDT::plf)

            ! Initialize 
            select type(plf)

            type is (PolygonLevelsetFunction2DClosedApproximationUDT)

                plf%options = options
                call plf%Initialize(ps)

            end select

        class default

        end select 

    end subroutine

    ! Visualization
    subroutine VisualizePolygonLevelsetFunction2D(plf, savefilepath)

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
        class(PolygonLevelsetFunction2DUDT)      :: plf 
        character(*), intent(in)                 :: savefilepath
        
        ! Auxiliary
        integer(I8)                         :: nx, ny
        real(R8), allocatable               :: xgv(:), ygv(:), xg(:), &
            yg(:), vg(:), xe(:, :), ye(:, :)
        real(R8)                            :: xmin, ymin, xmax, ymax, &
            offsetx, offsety, dx, dy, dxgv, dygv

        ! Loop
        integer(I8)                         :: k

        ! Construct a 2D grid
        !====================
        ! Set mesh size
        nx = 200
        ny = 400

        ! Allocate
        allocate(xgv(nx), ygv(ny), xg(nx*ny), yg(nx*ny), &
            vg(nx*ny))

        ! Set the gridding vectors
        call plf%ps%GetEdges(xe, ye)
        xmin = minval(xe)
        xmax = maxval(xe)
        ymin = minval(ye)
        ymax = maxval(ye)

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
        call plf%Evaluate(xg, yg, 0, 0, vg)

        ! Write data
        !===========
        call Write3DCoordinateData(xg, yg, vg, savefilepath)

        ! Housekeeping
        !=============
        ! Deallocate
        deallocate(xgv, ygv, xg, yg, vg)

    end subroutine

    !------------------------------------------------------------------!
    !                           OPTION ROUTINES                        !
    !------------------------------------------------------------------!

    ! Option setter, general
    subroutine SetOptionsGeneral(options)

        ! Description
        !============
        ! Set options for general plf description
        
        ! Declare variables
        !==================
        ! Arguments
        class(PLF2DGeneralOptionsUDT)   :: options 

    end subroutine

    ! Option setter, closed exact
    subroutine SetOptionsClosedExact(options)

        ! Description
        !============
        ! Set options for general plf description
        
        ! Declare variables
        !==================
        ! Arguments
        class(PLF2DClosedExactOptionsUDT)   :: options 

    end subroutine

    ! Options setter, closed approximation
    subroutine SetOptionsClosedApproximation(options, resx, resy, &
        offsetx, offsety, C, M, meth)

        ! Description
        !============
        ! Set options for closed approximation plf description. One 
        ! can set the interpolant resolution (x, y direction), offset
        ! relative to extend of x and y of the polygonset, and the 
        ! continuity parameters C and M of the polygon constructor. 
        
        ! Declare variables
        !==================
        ! Arguments
        class(PLF2DClosedApproximationOptionsUDT)   :: options 
        integer(I8), intent(in)                     :: resx, resy, C, M
        real(R8), intent(in)                        :: offsetx, offsety
        character(*), intent(in)            :: meth 
        
        ! Set options
        !============
        options%resx = resx 
        options%resy = resy 
        options%offsetx = offsetx 
        options%offsety = offsety 
        options%C = C 
        options%M = M 
        options%meth = meth

        call options%optionsClosedExact%Set()

    end subroutine

    !------------------------------------------------------------------!
    !                           PLF GENERAL                            !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializePLF2DGeneral(plf, ps)

        ! Description
        !============
        ! Initialize the general plf object. 
        
        ! Declare variables
        !==================
        ! Arguments
        class(PolygonLevelsetFunction2DGeneralUDT)  :: plf 
        type(PolygonSetUDT), intent(in)             :: ps 

        ! Auxiliary
        integer(I8)                                 :: nc

        real(R8), allocatable                       :: xe(:, :), &
            ye(:, :), nx(:), ny(:), nn(:), tx(:), ty(:), tn(:)

        ! Initialize
        !===========
        ! Add polygonset structure
        plf%ps = ps 

        ! Get edge data
        call ps%GetEdges(xe, ye)
        call ps%GetNormals(nx, ny, nn)
        call ps%GetTangents(tx, ty, tn)
        nc = size(xe, 1)
        
        ! Add to plf
        plf%x   = xe 
        plf%y   = ye
        plf%xf  = 0.5*( xe(:, 1) + xe(:, 2) )
        plf%yf  = 0.5*( ye(:, 1) + ye(:, 2) )
        plf%nxp = nx/nn ! normalize
        plf%nyp = ny/nn
        plf%txp = tx
        plf%typ = ty 
        plf%nc  = size(xe, 1)
        
    end subroutine

    ! Evaluation
    subroutine EvaluatePLF2DGeneral(plf, xq, yq, &
        derivx, derivy, vq)

        ! Description
        !============
        ! Evaluate the polygon levelset function in the query points xq, yq
        ! (value returned is in vq). The 'derivx' and 'derivy' input 
        ! arguments are used to specify the order of the derivative to 
        ! return. 

        ! Note that derivates are returned pointwise as there is no 
        ! dependence of the values on other sample points or polygon
        ! coordinates. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonLevelsetFunction2DGeneralUDT)  :: plf 
        real(R8), intent(in)                        :: xq(:), yq(:)
        real(R8), intent(out)                       :: vq(size(xq)) 
        integer(I8), intent(in)                     :: derivx, derivy 

        ! Auxiliary
        character(:), allocatable               :: deriv, derivxc, &
            derivyc 
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
        nc = plf%nc

        ! Allocate
        allocate(vnx(nc), vny(nc), vtx(nc, 2), vty(nc, 2), f(nc), &
            dvn(nc), ind(nc), fq(nq), vtxq(nq, 2), vtyq(nq, 2), &
            dvnq(nq), myones(nq), vnxq(nq), vnyq(nq), nxpq(nq), &
            nypq(nq), txpq(nq), typq(nq))
        
        ! Initialize
        myones(:) = 1

        ! Associate
        associate(&
            x       => plf%x,   &
            y       => plf%y,   &
            xf      => plf%xf,  &
            yf      => plf%yf,  &
            nxp     => plf%nxp, &
            nyp     => plf%nyp, &
            txp     => plf%txp, &
            typ     => plf%typ)

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
        if ( (derivx > 9) .or. (derivy > 9) ) then 
            ! Throw error, these derivatives should anyway not be necessary...
            call gdErrorHandler('Derivatives are of too high order, cannot compute')
        end if
            
        derivxc = ' ' 
        derivyc = ' ' 
        write(derivxc, '(I0)') derivx 
        write(derivyc, '(I0)') derivy 
        deriv = trim(derivxc) // trim(derivyc)

        select case(trim(deriv))

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

            ! Put to zero - should be higher order derivatve
            vq(:) = 0

        end select

        ! End associate
        end associate

        ! Deallocate
        deallocate(vnx, vny, vtx, vty, f, ind, fq, vtxq, vtyq, dvnq)

    end subroutine

    !------------------------------------------------------------------!
    !                         PLF CLOSED EXACT                         !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializePLF2DClosedExact(plf, ps)

        ! Description
        !============
        ! Initialize the general plf object. 
        
        ! Declare variables
        !==================
        ! Arguments
        class(PolygonLevelsetFunction2DClosedExactUDT)      :: plf 
        type(PolygonSetUDT), intent(in)             :: ps 

        ! Auxiliary
        integer(I8)                                 :: np, flag, ne

        real(R8), allocatable                       :: xe(:, :), &
            ye(:, :), nx(:), ny(:), nn(:), tx(:), ty(:), tn(:), xp(:), &
            yp(:), tempx(:), tempy(:), tempnx(:, :), tempny(:, :), &
            nxpv(:, :), nypv(:, :), theta0(:), crossprod(:)

        ! Loop
        integer(I8)                                 :: i, ce, ce2

        ! Initialize
        !===========
        ! Check polygonset 
        call ps%OrientNestedClosedPolygons(flag)
        if (flag > 0) then 
            ! Polygon set does not comply, throw error
            call gdErrorHandler('InitializePLF2DClosedExact: polygon set should be closed and non-intersecting')
        end if 
        
        ! Add polygonset structure
        plf%ps = ps 

        ! Associate
        associate(&
            pol     => ps%polygons &
        )

        ! Compute additional data
        !========================
        ! Get edge data
        call ps%GetEdges(xe, ye) ! may still change
        call ps%GetNormals(nx, ny, nn) ! may still change
        call ps%GetTangents(tx, ty, tn)
        !call ps%GetPoints(xp, yp) ! may still change

        ! Allocate
        np = (size(xe, 1)+1)*ps%np
        allocate(nxpv(np, 2), nypv(np, 2), xp(np), yp(np))

        ! Compute vertex regions
        ce = 0
        ce2 = 0
        do i = 1, ps%np 
            ! Normally, since each polygon is closed and non-intersecting,
            ! the first and last vertex will be the same and no 
            ! duplicate vertices will exist in vert. 
            
            ! Unpack
            ne = pol(i)%ne 

            ! Coordinates
            tempx = pol(i)%x(pol(i)%vert)
            tempy = pol(i)%y(pol(i)%vert)

            ! Normals
            allocate(tempnx(ne+1, 2), tempny(ne+1, 2))

            tempnx(1:ne, 2) = pol(i)%nx/pol(i)%nn 
            tempnx(1:ne, 1) = [pol(i)%nx(ne)/pol(i)%nn(ne), &
                pol(i)%nx(1:ne-1)/pol(i)%nn(1:ne-1)]
            tempnx(ne+1, 1:2) = [pol(i)%nx(1)/pol(i)%nn(1), pol(i)%nx(ne)/pol(i)%nn(ne)]

            tempny(1:ne, 2) = pol(i)%ny/pol(i)%nn 
            tempny(1:ne, 1) = [pol(i)%ny(ne)/pol(i)%nn(ne), &
                pol(i)%ny(1:ne-1)/pol(i)%nn(1:ne-1)]
            tempny(ne+1, 1:2) = [pol(i)%ny(1)/pol(i)%nn(1), pol(i)%ny(ne)/pol(i)%nn(ne)]


            ! Add
            xp(ce+1:ce+ne+1) = tempx 
            yp(ce+1:ce+ne+1) = tempy 
            nxpv(ce+1:ce+ne+1, :) = tempnx ! normals should be oriented outwards already
            nypv(ce+1:ce+ne+1, :) = tempny 
            nx(ce2+1:ce2+ne) = pol(i)%nx/pol(i)%nn 
            ny(ce2+1:ce2+ne) = pol(i)%ny/pol(i)%nn 
            tn(ce2+1:ce2+ne) = pol(i)%nn 

            ! Update counter
            ce = ce + ne + 1
            ce2 = ce2 + ne

            ! Housekeeping
            deallocate(tempnx, tempny)

            
            
        end do 

        ! Compute cross product
        crossprod = nxpv(:, 2)*nypv(:, 1) - nxpv(:, 1)*nypv(:, 2)

        ! Compute angles
        theta0 = atan2(nxpv(:, 2)*nypv(:, 1) - nxpv(:, 1)*nypv(:, 2), &
            nxpv(:, 1)*nxpv(:, 2) + nypv(:, 1)*nypv(:, 2))
        where (theta0 < 0) theta0 = theta0 + 2*pi         
        
        ! Add to plf
        plf%xp  = xp 
        plf%yp  = yp
        plf%xf  = 0.5*( xe(:, 1) + xe(:, 2) )
        plf%yf  = 0.5*( ye(:, 1) + ye(:, 2) )
        plf%nxp = nx 
        plf%nyp = ny
        plf%tnp = tn 
        plf%nxpv    = nxpv 
        plf%nypv    = nypv
        plf%thetav  = theta0 
        plf%crossprod   = crossprod

        ! Housekeeping
        !=============
        end associate

    end subroutine 

    ! Evaluation
    subroutine EvaluatePLF2DClosedExact(plf, xq, yq, derivx, derivy, vq)

        ! Description
        !============
        ! Evaluate the polygon levelset function in the query points xq, yq
        ! (value returned is in vq). The 'derivx' and 'derivy' input 
        ! arguments are used to specify the order of the derivative to 
        ! return. 

        ! Note that derivates are returned pointwise as there is no 
        ! dependence of the values on other sample points or polygon
        ! coordinates.
        
        use iso_fortran_env
        use, intrinsic :: ieee_arithmetic

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonLevelsetFunction2DClosedExactUDT)  :: plf 
        real(R8), intent(in)                        :: xq(:), yq(:)
        real(R8), intent(out)                       :: vq(size(xq)) 
        integer(I8), intent(in)                     :: derivx, derivy 

        ! Auxiliary
        character(:), allocatable                   :: deriv, derivxc, &
            derivyc

        integer(I8)                             :: nq, np, &
            indmine, indminv, indmin
        integer(I8), allocatable                :: eind(:), vind(:), &
            minind(:)

        real(R8)                                :: inf, signe, signv, &
            fv, fe, xqr, yqr, totsign(1:2)
        real(R8), parameter                     :: myone = 1
        real(R8), allocatable                   :: myones(:), dvn(:), &
            tdistvert(:), tcrossprod(:), vx(:), vy(:), tvn(:), &
            dx(:), dy(:), theta(:), distedge(:), distvert(:), &
            val(:)
        

        logical, allocatable                    :: isinvert(:), &
            onedge(:), te(:), tv(:)

        ! Loop
        integer(I8)                             :: iq

        ! Data
        inf = ieee_value(inf, ieee_positive_inf)

        ! Initialize
        !===========
        ! Associate
        associate(&
            xp      => plf%xp,      &
            yp      => plf%yp,      &
            nxp     => plf%nxp,     &
            nyp     => plf%nyp,     &
            nxpv    => plf%nxpv,    &
            nypv    => plf%nypv,    &
            tnp     => plf%tnp,     &
            theta0  => plf%thetav,  &
            crossprod => plf%crossprod, &
            xf      => plf%xf,      &
            yf      => plf%yf       &
            )

        ! Get sizes
        np = size(xp, 1)
        nq = size(xq, 1)
        
        ! Allocate and initialize
        allocate(eind(nq), vind(nq), minind(nq), tdistvert(nq), &
            tcrossprod(nq), val(nq))

        ! Check what to compute exactly
        if ( (derivx > 9) .or. (derivy > 9) ) then 
            ! Throw error, these derivatives should anyway not be necessary...
            call gdErrorHandler('Derivatives are of too high order, cannot compute')
        end if
        derivxc = ' ' 
        derivyc = ' ' 
        write(derivxc, '(I0)') derivx 
        write(derivyc, '(I0)') derivy 
        deriv = trim(derivxc) // trim(derivyc)

        ! Loop 
        !=====
        allocate(vx(np), vy(np), dvn(np), tvn(np), dx(np), dy(np), &
            theta(np), isinvert(np), onedge(np), distedge(np), &
            distvert(np), myones(np))
        myones = 1
        do iq = 1, nq
            
            ! Unpack
            xqr = xq(iq)
            yqr = yq(iq)

            ! Construct vector between edge centers and points
            vx = (xqr - xf)
            vy = (yqr - yf)
            dvn = (vx*nxp + vy*nyp)
            tvn = vx*nyp - vy*nxp ! distance in tangential direction
            
            ! Check in which regions the query points lie
            ! Vertex regions
            dx = xqr - xp
            dy = yqr - yp
            theta = atan2(dx*nypv(:, 1) - dy*nxpv(:, 1), dx*nxpv(:, 1) + dy*nypv(:, 1))
            where (theta < 0) theta = theta + 2*pi
            isinvert = theta < theta0
            
            ! Edge regions
            onedge = (abs(tvn) <= 0.5*tnp)
            
            ! Compute distances
            distedge = dvn
            distvert = sign(myones, crossprod)*sqrt((xp - xqr)**2 + (yp - yqr)**2)
            where (.not. onedge) distedge = inf 
            where (.not. isinvert) distvert = inf 
            
            ! Compute minimal distance
            indmine = minloc(abs(distedge), 1)
            fe = minval(abs(distedge))
            signe = sign(myone, distedge(indmine))
            eind(iq) = indmine

            indminv = minloc(abs(distvert), 1)
            fv = minval(abs(distvert))
            tdistvert(iq) = distvert(indminv)
            tcrossprod(iq) = crossprod(indminv)

            signv = sign(myone, distvert(indminv))
            vind(iq) = indminv
            totsign = [signe, signv]
            val(iq) = minval([fe, fv])
            indmin = minloc([fe, fv], 1)
            minind(iq) = indmin 
            val(iq) = val(iq)*totsign(indmin)
            
        end do 

        ! Construct function handle depending on type
        select case (trim(deriv))
        
        case ('00')
            
            ! Function value
            vq = val
            
        case ('10')
            
            ! First order derivative w.r.t. xq
            
            ! Compute
            te = minind == 1
            tv = minind == 2
            
            do iq = 1, nq
                if (te(iq)) then 
                    vq(iq) = nxp(eind(iq))
                else
                    vq(iq) = -sign(myone, tdistvert(iq))*sign(myone, tcrossprod(iq))/tdistvert(iq)*(xp(vind(iq)) - xq(iq))
                end if
            end do 
            
        case ('01')
            
            ! First order derivative w.r.t. yq
            
            ! Compute
            te = minind == 1
            tv = minind == 2

            do iq = 1, nq
                if (te(iq)) then 
                    vq(iq) = nyp(eind(iq))
                else
                    vq(iq) = -sign(myone, tdistvert(iq))*sign(myone, tcrossprod(iq))/tdistvert(iq)*(yp(vind(iq)) - yq(iq))
                end if
            end do 

        case ('20')
            
            ! Second order derivative w.r.t. xq, xq
            vq = 0
            te = minind == 1
            tv = minind == 2
            do iq = 1, nq
                if (tv(iq)) then 
                    vq(iq) = -sign(myone, tdistvert(iq))*&
                    sign(myone, tcrossprod(iq))/tdistvert(iq)**3&
                    *(xp(vind(iq)) - xq(iq))**2 + &
                    sign(myone, tdistvert(iq))*sign(myone, tcrossprod(iq))/tdistvert(iq)
                end if
            end do 

        case ('11')
            
            ! Second order derivative w.r.t. xq, yq
            vq = 0
            te = minind == 1;
            tv = minind == 2;
            do iq = 1, nq
                if (tv(iq)) then 
                    vq(iq) = -sign(myone, tdistvert(iq))*&
                    sign(myone, tcrossprod(iq))/tdistvert(iq)**3&
                    *(xp(vind(iq)) - xq(iq))*(yp(vind(iq)) - xq(iq))
                end if
            end do 
            
        case ('02')
            
            vq = 0
            te = minind == 1;
            tv = minind == 2;

            do iq = 1, nq
                if (tv(iq)) then 
                    vq(iq) = -sign(myone, tdistvert(iq))*&
                    sign(myone, tcrossprod(iq))/tdistvert(iq)**3&
                    *(yp(vind(iq)) - yq(iq))**2 + &
                    sign(myone, tdistvert(iq))*sign(myone, tcrossprod(iq))/tdistvert(iq)
                end if
            end do 

            
        case default

            ! Just for now
            vq = 0

        end select

        ! Housekeeping
        !=============
        end associate

    end subroutine

    !------------------------------------------------------------------!
    !                     PLF CLOSED APPROXIMATION                     !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializePLF2DClosedApproximation(plf, ps)

        ! Description
        !============
        ! Initialize the general plf object. 
        
        ! Declare variables
        !==================
        ! Arguments
        class(PolygonLevelsetFunction2DClosedApproximationUDT)  :: plf 
        type(PolygonSetUDT), intent(in)             :: ps 

        ! Auxiliary
        class(PolygonLevelsetFunction2DUDT), allocatable :: plfe

        integer(I8)                         :: nx, ny

        real(R8)                            :: minx, maxx, miny, maxy, &
            dx, dy
        real(R8), allocatable               :: xp(:), yp(:), xg(:), &
            yg(:), xgv(:), ygv(:), vg(:), vg2D(:, :)

        ! Loop
        integer(I8)                         :: k

        ! Initialize
        !===========
        ! Add polygonset structure
        plf%ps = ps 

        ! Get polygon vertices for later
        call ps%GetPoints(xp, yp)

        ! Construct closed exact levelset function
        call InitializePolygonLevelsetFunction2D(plfe, ps, &
            plf%options%optionsClosedExact)

        ! Evaluate exact representation
        !==============================
        ! Set dimensions of vertices (assumed resx, resy are number of 
        ! cells)
        nx = plf%options%resx + 1
        ny = plf%options%resy + 1

        ! Allocate
        allocate(xgv(nx), ygv(ny), xg(nx*ny), yg(nx*ny), vg(nx*ny))

        ! Get polygon extent
        minx = minval(xp) 
        maxx = maxval(xp)
        miny = minval(yp) 
        maxy = maxval(yp)
        dx = maxx - minx 
        dy = maxy - miny

        ! Check how far the grid should extend
        minx = minx - plf%options%offsetx*dx 
        maxx = maxx + plf%options%offsetx*dx 
        miny = miny - plf%options%offsety*dy 
        maxy = maxy + plf%options%offsety*dy 
        
        ! Reset dx, dy
        dx = maxx - minx 
        dy = maxy - miny
        
        ! Gridding vectors
        xgv = [(k, k = 0, plf%options%resx)]*dx/plf%options%resx + minx 
        ygv = [(k, k = 0, plf%options%resy)]*dy/plf%options%resy + miny 
        
        call Construct2DStructuredGrid(xgv, ygv, nx, ny, xg, yg)

        ! Evaluate exact representation
        call plfe%Evaluate(xg, yg, 0, 0, vg)

        ! Reshape
        vg2D = reshape(vg, (/nx, ny/))

        ! Build polygon representation
        !=============================
        ! Construct
        call plf%interp%SetParameters(plf%options%meth, plf%options%C, plf%options%M)
        call plf%interp%Construct(xgv, ygv, vg2D)
        

    end subroutine

    ! Evaluation
    subroutine EvaluatePLF2DClosedApproximation(plf, xq, yq, derivx, &
        derivy, vq)

        ! Description
        !============
        ! Evaluate - simply call interpolant routine

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonLevelsetFunction2DClosedApproximationUDT)  :: plf 
        real(R8), intent(in)            :: xq(:), yq(:)
        integer(I8), intent(in)         :: derivx, derivy
        real(R8), intent(out)           :: vq(size(xq))

        ! Call interpolant
        !=================
        call plf%interp%Evaluate(xq, yq, derivx, derivy, vq)
 
    end subroutine



end module