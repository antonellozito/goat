!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module provides the 'polygonset' type (object), which can be used 
! to handle multiple polygons and polygon-like shapes. Each polygonset 
! contains a single or multiple polygons ('polygon' type)

! The coordinates may contain duplicate points, and the
! polygon(s) may intersect, be closed (this requires at least two duplicate
! points at the start and end), consist of multiple segments (separated by
! NaNs). Polygons can only be constructed by parsing x, y coordinate 
! lists, which may contain NaNs (IEEE standard) to distinguish 
! separate polygons. Due to the nature of the input,
! true branching polygons cannot be generated (which is often
! also not desired). However, by adding duplicate points, the polygon can
! in fact contain multiple times the same point and based on that self
! intersect and become non-simple. Note that a polygon is called 'simple'
! here if each vertex only occurs maximally twice in the polygon (so it can
! still intersect, simply not in a single node). 

! A polygonset object has the following fields
!       
!       np          number of distinct polygons
!       polygons    an np-by-1 structure containing np polygon objects
!       

! A polygon object has the following fields: 
!
!       nv, ne      number of unique polygon vertices (nv) and polygon 
!                   edges (ne)
!       x, y        vertex coordinates of the polygon (nv-by-1)
!       edge        vertex pairs for each edge (ne-by-1, sorted in
!                   arbitrary direction)
!       vert        vertex ordering of the polygon (ne+1-by-1), may contain
!                   duplicates if the polygon is closed/non simple
!       ex, ey      polygon edge centers
!       tx, ty      polygon tangents (not normalized)
!       nx, ny      polygon normals (pointing to interior of polygon, not
!                   normalized)
!       tn, nn      tangent and normal length (for normalization if
!                   desired)
!       
!       isclosed    logical that indicates if the polygon is closed upon
!                   itself (i.e. the first and final edge have a vertex in
!                   common)
!       selfintersecting    logical that indicates if the polygon
!                           intersects with itself. Note: polygons that are
!                           closed but do not intersect anywhere else are
!                           not considered to be self intersecting. 
!       simplepolygon       logical that indicates if the polygon is a
!                           simple polygon, i.e. that each vertex belongs
!                           to maximally 2 edges. 

! The following routines are provided:
!       
!       ConstructPolygonSet(p, x, y) Polygon constructor, x and y should 
!                                   be 1D arrays of equal size, p is the
!                                   polygonset object 

! Notes
!======
! Note 1: a polygon should consist of at least one edge (i.e. 2 
! vertices). Otherwise, an error is thrown during construction. 

! Note 2: NaNs should be used to separate polygons. Therefore, they 
! should not occur in the first and last element of the coordinate 
! vectors (indicating an empty polygon, which is not supported, see 
! note 1). This is checked for and an error will be thrown. 

module mod_polygon

    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
    use mod_plotter

    ! The usual
    implicit none
    save
    public 

    ! Constants
    real(R8), private, parameter        :: disttol = 1e-8 ! tolerance when computing distances
    real(R8), private, parameter        :: macheps = 1e-12 ! 'machine' precision

    !==================================================================!
    !                                                                  !
    !                                TYPES                             !
    !                                                                  !
    !==================================================================!

    ! The polygon type
    type PolygonUDT 

        ! The general polygon type, see documentation of this module for
        ! the fields and procedures. 

        integer(I8)                 :: nv, ne

        logical                     :: isclosed, selfintersecting, &
            issimple

        integer(I8), allocatable    :: edges(:, :), vert(:)

        real(R8), allocatable       :: x(:), y(:), ex(:), ey(:), &
            tx(:), ty(:), tn(:), nx(:), ny(:), nn(:)
        

    contains 

        procedure :: Allocate       => AllocatePolygon
        procedure :: Deallocate     => DeallocatePolygon
        procedure :: Construct      => ConstructPolygon
        procedure :: Initialize     => InitializePolygon 
        procedure :: ComputeMetrics => ComputePolygonMetrics
        procedure :: SetVert        => SetPolygonVertices
        procedure :: RemoveDuplicatePoints
        procedure :: IsClosedPolygon
        procedure :: IsSimplePolygon
        procedure :: Inpolygon

    end type 

    ! The polygonset type
    type PolygonSetUDT 

        ! The general polygonset type, see documentation of this module
        ! for the fields and procedures

        integer(I8)                     :: np 
        type(PolygonUDT), allocatable   :: polygons(:)

    contains 

        procedure :: Construct      => ConstructPolygonSet

    end type 


    !==================================================================!
    !                                                                  !
    !                             INTERFACES                           !
    !                                                                  !
    !==================================================================!

    ! Distance function for different inputs
    interface Distance 
        module procedure Distance0D
        module procedure Distance1D
        ! module procedure Distance2D
    end interface


    contains 

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                         PolygonSet routines                      !
    !------------------------------------------------------------------!

    subroutine ConstructPolygonSet(polygonset, x, y)

        ! Description
        !============
        ! Construct a set of polygons from the given x and y 
        ! coordinates. NaNs should be present to indicate different 
        ! polygon pieces. First, the number of polygons are checked 
        ! by computing the number of NaNs etc. Afterwards, each polygon
        ! is constructed separately, including all metric and logical
        ! data. 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), allocatable, intent(in)   :: x(:), y(:) 
        class(PolygonSetUDT), intent(inout)  :: polygonset 

        ! Auxiliary 
        integer(I8)                     :: nnans, nx, ny
        integer(I8), allocatable        :: nanloc(:), startind(:), &
            endind(:), nvpp(:) 

        real(R8), allocatable           :: tempx(:), tempy(:)

        ! Loop
        integer(I8)                     :: i, k


        ! Initialize
        !===========
        ! Compute and check sizes
        nx = size(x)
        ny = size(y)
        if (nx .ne. ny) then 
            ! Call the error handler
            call PolygonErrorHandler('x and y should have equal number of elements')
        end if 

        ! Count NaNs and check
        nnans = count(isnan(x))
        allocate(nanloc(nnans))
        nanloc = pack([(k, k = 1, nx)], isnan(x))
        if (any(isnan(x) .neqv. isnan(y))) then 
            call PolygonErrorHandler('NaN values do not correspond in x and y coordinates')
        end if
        
        ! Determine polygons
        !===================
        ! Compute number of polygons
        polygonset%np = nnans + 1 

        ! Checks
        if (isnan(x(1)) .or. isnan(x(nx))) then  
            ! NaNs in first or last entry, throw error
            call PolygonErrorHandler('NaN values detected in first ' &
                // 'and/or last element of coordinate array. ' &
                // 'Not supported, check input')
        end if 

        ! Initialize polygons
        allocate(startind(polygonset%np), endind(polygonset%np), &
            nvpp(polygonset%np))
        startind    = [1, nanloc+1]
        endind      = [nanloc-1, nx] 
        nvpp        = endind - startind + 1
        allocate(polygonset%polygons(polygonset%np))

        ! Check if all polygons have more than one vertex 
        if (any(nvpp < 2)) then 
            call PolygonErrorHandler('Polygons with less than two ' &
                // 'vertices found, not supported. Check input.')
        end if 
        
        ! Loop over all polygons
        do i = 1, polygonset%np 

            ! Build polygon
            allocate(tempx(nvpp(i)), tempy(nvpp(i)))
            tempx = x(startind(i):endind(i))
            tempy = y(startind(i):endind(i))
            call polygonset%polygons(i)%Construct(tempx, tempy)
            deallocate(tempx, tempy)

        end do 

        ! Housekeeping
        !=============
        deallocate(startind, endind, nvpp, nanloc)

    end subroutine

    !------------------------------------------------------------------!
    !                          Polygon routines                        !
    !------------------------------------------------------------------!

    ! Constructor
    subroutine ConstructPolygon(polygon, x, y) 

        ! Description
        !============
        ! Construct a single polygon from the given x and y coordinates.
        ! It is assumed that there are no NaNs present anymore in the 
        ! coordinates and that the sizes of x and y are the same. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT)                       :: polygon 
        real(R8), allocatable, intent(in)       :: x(:), y(:) 

        ! Auxiliary

        ! Loop

        ! Construct
        !==========
        ! Set coordinates & initialize vertices etc
        call polygon%Initialize(x, y)

        ! Check for duplicate points, remove
        call polygon%RemoveDuplicatePoints()

        ! Extract the vertices
        call polygon%SetVert()

        ! Check if the polygon is closed
        call polygon%IsClosedPolygon()

        ! Construct the polygon metrics
        call polygon%ComputeMetrics()

        ! Check if the polygon is simple
        call polygon%IsSimplePolygon()

        ! Check if the polygon self intersects
        !call polygon%IsSelfIntersectingPolygon()


    end subroutine

    ! Allocator
    subroutine AllocatePolygon(polygon)

        ! Description
        !============
        ! Allocate the polygon fields for a given number of vertices nv
        ! and edges ne. These numbers must be attributed beforehand. 

        ! Declare variables
        !=================
        ! Arguments
        class(PolygonUDT)        :: polygon


        ! Allocate
        !=========
        associate(nv => polygon%nv, ne => polygon%ne)

        allocate(polygon%x(nv), polygon%y(nv), polygon%vert(ne+1), &
            polygon%edges(ne, 2), polygon%ex(ne), polygon%ey(ne), &
            polygon%tx(ne), polygon%ty(ne), polygon%tn(ne), &
            polygon%nx(ne), polygon%ny(ne), polygon%nn(ne))

        end associate


    end subroutine

    ! Deallocator
    subroutine DeallocatePolygon(polygon)

        ! Description
        !============
        ! Deallocate polygon fields

        ! Declare variables
        !=================
        ! Arguments
        class(PolygonUDT)        :: polygon

        ! Deallocate
        !===========
        deallocate(polygon%x, polygon%y, polygon%vert, &
            polygon%edges, polygon%ex, polygon%ey, &
            polygon%tx, polygon%ty, polygon%tn, &
            polygon%nx, polygon%ny, polygon%nn)


    end subroutine

    ! Initialization
    subroutine InitializePolygon(polygon, x, y)

        ! Description
        !============
        ! Initialize the polygon by allocating the required fields, 
        ! computing the initial polygon dimensions, and adding the 
        ! coordinates. It is assumed that all quantities are compliant
        ! with the polygon definition. 
        
        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT)                   :: polygon 
        real(R8), allocatable, intent(in)   :: x(:), y(:)

        ! Auxiliary

        ! Loop
        integer(I8)                         :: k 

        ! Initialize
        !===========
        ! Sizes
        polygon%nv = size(x)
        polygon%ne = polygon%nv - 1

        ! Allocate
        call polygon%Allocate()

        ! Set coordinates
        polygon%x = x 
        polygon%y = y 

        ! Set edges
        polygon%edges(:, 1) = [(k, k = 1, polygon%nv-1)]
        polygon%edges(:, 2) = [(k, k = 2, polygon%nv)]


    end subroutine

    ! Duplicate point remover
    subroutine RemoveDuplicatePoints(polygon)

        ! Description
        !============
        ! This function checks whether coordinates of a polygon appear multiple
        ! times and replaces them by a single coordinate (and renumbers the
        ! vertices accordingly). If the polygon is destroyed in the process (i.e.
        ! all points are duplicate), an error is thrown.
        ! Note that this is an expensive routine
        ! and should only be called once during construction of the polygon! 

        ! Input
        !------
        ! - polygon:        polygon structure with the x, y, vert fields

        ! Output
        !-------
        ! - polygon:        same structure, but with duplicate points removed and
        !                   empty polygons destroyed.

        ! Algorithm
        !----------
        ! Loop over all points, compute the distance w.r.t. each consecutive point
        ! in the polygon. If the distance is smaller than a small tolerance, the
        ! point is considered to be duplicate. In that case, the next point is
        ! adapted to be the current point, the vertices are changed, and the
        ! algorithm continues on the next points. 

        ! Notes 
        !------
        ! Note 1: the polygon metrics are *not* updated in this routine - this has
        ! to be done afterwards!

        ! Note 2: when duplicate points are removed, the polygon has to 
        ! be rebuilt and dimensions of arrays will change. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT)           :: polygon 

        ! Auxiliary
        real(R8)                    :: d 

        real(R8), allocatable       :: tempx(:), tempy(:)

        integer(I8), allocatable    :: mapping(:), tempedge(:, :), &
            diffindex(:)

        logical, allocatable        :: keepvertex(:), keepedges(:)

        ! Loop
        integer(I8)                 :: i, j, k

        ! Find duplicate points
        !======================
        ! Associate
        associate(&
            x       => polygon%x,       &
            y       => polygon%y,       &
            nv      => polygon%nv,      &
            ne      => polygon%ne,      & 
            edges   => polygon%edges)

        ! Initialize
        allocate(mapping(nv), keepvertex(nv), diffindex(nv)) 
        mapping(:) = [(k, k = 1, nv)]
        keepvertex(:) = .true.
        diffindex(:) = 0
        
        ! Loop over all points but one 
        do i = 1, nv-1 
            ! Compute distance to the next coordinates
            do j = i+1, nv 
                ! Compute the distance
                call Distance(d, x(i), y(i), x(j), y(j))

                ! Check if distance is smaller than distance tolerance
                if (d < disttol) then 
                    ! Mark for removal
                    keepvertex(j) = .false. 

                    ! Update mapping
                    mapping(j) = mapping(i) ! not just i - i can already be mapped to another vertex

                end if
            end do
        end do

        ! Stop associate
        end associate

        ! Remove & rebuild
        !=================
        ! Compute offset for new vertices (renumbering due to deletion)
        where (.not. keepvertex) diffindex = 1
        do i = 2, polygon%nv 
            diffindex(i) = diffindex(i-1) + diffindex(i)
        end do

        ! Remap the vertices in the edges
        polygon%edges(:, 1) = mapping(polygon%edges(:, 1))
        polygon%edges(:, 2) = mapping(polygon%edges(:, 2))
        polygon%edges(:, 1) = polygon%edges(:, 1) - diffindex(polygon%edges(:, 1))
        polygon%edges(:, 2) = polygon%edges(:, 2) - diffindex(polygon%edges(:, 2))


        ! Check if any edges have to be removed
        keepedges = polygon%edges(:, 1) .ne. polygon%edges(:, 2)
        if (.not. all(keepedges)) then 
            call PolygonWarningHandler('RemoveDuplicatePoints: edges ' &
                // 'removed that consist of only duplicate points')
        end if 

        ! Check if the polygon is destroyed in the process
        if (count(keepvertex) < 2) then 
            call PolygonErrorHandler('RemoveDuplicatePoints: removal ' &
                // 'of duplicate points resulted in polygon ' &
                // 'destruction - check input.')
        end if 

        ! Rebuild the polygon
        allocate(tempedge(count(keepedges), 2), &
            tempx(count(keepvertex)), tempy(count(keepvertex)))
        tempedge(:, 1) = pack(polygon%edges(:, 1), keepedges)
        tempedge(:, 2) = pack(polygon%edges(:, 2), keepedges)
        tempx = pack(polygon%x, keepvertex)
        tempy = pack(polygon%y, keepvertex)

        call polygon%Deallocate() 
        polygon%nv = count(keepvertex)
        polygon%ne = count(keepedges)
        call polygon%Allocate() 
        polygon%edges = tempedge 
        polygon%x = tempx 
        polygon%y = tempy

        ! Housekeeping
        !=============
        deallocate(tempedge, tempx, tempy, diffindex, keepvertex, &
            keepedges)


    end subroutine

    ! Vertex setter
    subroutine SetPolygonVertices(polygon)

        ! Description
        !============
        ! This routine extracts the subsequent set of vertices of a 
        ! polygon from its edges. This set of vertices is then stored
        ! in the polygon%vert field. Note that the size of this field
        ! is not necessarily equal to the size of polygon%x! This is 
        ! because multiple points can occur multiple times in the
        ! polygon. 

        ! It is assumed that the polygon has been properly initialized
        ! and that any duplicate points have been removed (will work 
        ! just as well without this last requirement, but then vert will
        ! change after removing duplicate points). 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT)                   :: polygon 

        ! Auxiliary
        integer(I8)                         :: cvind, nv  
        logical                             :: check(1:4)

        ! Loop
        integer(I8)                         :: k 

        ! Initialize
        !===========
        ! Associate
        associate( & 
            ne          => polygon%ne,      &
            edges       => polygon%edges,   &
            vert        => polygon%vert) 
        

        ! Check
        check(:) = .false. 

        ! Construct vert
        !===============
        ! Get the starting edge
        check(1) = edges(1, 1) == edges(2, 1)
        check(2) = edges(1, 1) == edges(2, 2)
        check(3) = edges(1, 2) == edges(2, 1)
        check(4) = edges(1, 2) == edges(2, 2)

        ! Do checks
        if ( (check(1) .and. ( .not. any(check(2:4))) ) .or. &
                (check(2) .and. ( .not. any(check([1, 3, 4]))) ) ) then 
            ! Second vertex is starting vertex
            cvind = 2
        elseif ( (check(4) .and. ( .not. any(check(1:3))) ) .or. &
                (check(3) .and. ( .not. any(check([1, 2, 4]))) ) ) then 
            ! First vertex is starting vertex
            cvind = 1
        else
            ! Something wrong with input - throw error
            call PolygonErrorHandler('ExtractPolygonVertices: ' &
                // 'something wrong with input - please check input ' &
                // 'polygon')
        end if

        ! Loop 
        k = 1
        do while (k < ne) 
            ! Get the current vertex and add
            vert(k) = edges(k, cvind) 

            ! Get the next vertex
            if (cvind == 1) then 
                nv = edges(k, 2)
            else 
                nv = edges(k, 1)
            end if 

            ! Update cvind 
            if (nv == edges(k+1, 1)) then 
                cvind = 1
            elseif (nv == edges(k+1, 2)) then 
                cvind = 2
            else 
                ! Throw error, next edge does not contain the vertex
                call PolygonErrorHandler('ExtractPolygonVertices: ' &
                    // 'could not find next edge - check input polygon')
            end if 

            ! Update counter
            k = k + 1
        end do 
 
        ! Add last vertices
        vert(k) = edges(ne, cvind)
        if (cvind == 1) then 
            vert(k+1) = edges(ne, 2)
        else 
            vert(k+1) = edges(ne, 1)
        end if 

        end associate 

    end subroutine

    ! Is polygon closed
    subroutine IsClosedPolygon(polygon)

        ! Description
        !============
        ! Check if the polygon is closed upon itself. It is assumed that the
        ! polygon's duplicate vertices are removed and that as such simply the
        ! start and end edge can be checked for a duplicate vertex. If one (or more
        ! in the case of a polygon with two edges with different orientation) are
        ! found, the polygon is said to be closed. This is then updated in the
        ! polygon%isclosed logical field.

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT)               :: polygon

        ! Auxiliary
        integer(I8)                     :: sv(1:2), ev(1:2), nv 

        ! Checks
        !=======
        if (size(polygon%edges, 1) > 1) then 
            ! Polygon with multiple edges, need to check start and end

            ! Start and end vertices
            sv = polygon%edges(1, :)
            ev = polygon%edges(polygon%ne, :)

            ! Compare
            nv = count(sv(1) == ev) + count(sv(2) == ev)

            ! Set logical
            polygon%isclosed = nv >= 1

        else 
            ! Polygon with a single edge -> false
            polygon%isclosed = .false. 

        end if 


    end subroutine

    ! Metrics computation
    subroutine ComputePolygonMetrics(polygon)

        ! Description
        !============
        ! Compute the geometric quantities (or metrics) of the polygon. It is
        ! assumed that the polygon's duplicate points have been removed (otherwise
        ! it will work just as well, but NaNs may occur in the tangents/normals). 

        ! Algorithm
        !==========
        ! Most things are quite straightforward computations. To determine the
        ! inward pointing normal, we construct for each edge pair a point projected
        ! along the current normal and determine whether point lies within the 
        ! polygonal region using 'inpolygon'. If this is not the case, the sign of
        ! the normal is reversed. 
        
        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT)                   :: polygon 

        ! Auxiliary
        real(R8)                            :: testdist
        real(R8), allocatable, dimension(:) :: xe1, xe2, ye1, ye2, &
            temptn, tpx, tpy

        logical, allocatable                :: isin(:)

        ! Loop
        integer(I8)                         :: k 

        ! Initialize
        !===========
        ! Distance to construct test points for normal
        testdist = disttol*100 ! larger than tolerance 

        ! Check to be sure if the polygon is closed
        call polygon%IsClosedPolygon()

        ! Compute metrics
        !================
        ! Associate
        associate(&
            x           => polygon%x,       &
            y           => polygon%y,       &
            ne          => polygon%ne,      &
            tx          => polygon%tx,      &
            ty          => polygon%ty,      &
            tn          => polygon%tn,      &
            nx          => polygon%nx,      &
            ny          => polygon%ny,      &
            nn          => polygon%nn,      &
            ex          => polygon%ex,      &
            ey          => polygon%ey,      &
            edges       => polygon%edges,   &
            vert        => polygon%vert)

        ! Polygon edge coordinates
        allocate(xe1(ne), ye1(ne), xe2(ne), ye2(ne), temptn(ne))
        xe1 = x(edges(:, 1)) 
        xe2 = x(edges(:, 2))
        ye1 = y(edges(:, 1)) 
        ye2 = y(edges(:, 2))

        ! Edge centers
        ex = 0.5*(xe1 + xe2)
        ey = 0.5*(ye1 + ye2)

        ! Tangents
        tx = xe2 - xe1
        ty = ye2- ye1
        call Distance(temptn, xe1, ye1, xe2, ye2)
        tn = temptn

        ! Initial normals
        nx = -ty
        ny = tx
        nn = tn 

        ! Determine normal orientation
        !=============================
        ! Construct test points
        allocate(tpx(ne), tpy(ne))
        tpx = ex + testdist*nx/nn 
        tpy = ey + testdist*ny/nn 

        ! See if points lie in polygon
        call polygon%Inpolygon(tpx, tpy, isin)
        where (.not. isin) nx = -nx 
        where (.not. isin) ny = -ny


        end associate

        ! Housekeeping
        !=============
        deallocate(xe1, ye1, xe2, ye2, tpx, tpy, temptn)

    end subroutine

    ! Simple polygon determination
    subroutine IsSimplePolygon(polygon)

        ! Description
        !============
        ! Check if the polygon is a simple polygon, meaning that each vertex only
        ! belongs to maximally two edges. It is assumed that duplicate vertices and
        ! pathological edges have already been removed, and that therefore this
        ! check can be done on a vertex ID basis only. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT)           :: polygon 

        ! Auxiliary
        integer(I8), allocatable    :: counter(:) 

        ! Loop
        integer(I8)                 :: i


        ! Initialize
        !===========
        allocate(counter(size(polygon%x))) 

        ! Check
        !======
        ! Count vertex occurrences 
        counter(:) = 0
        do i = 1, size(polygon%edges, 1) 
            counter(polygon%edges(i, 1)) = counter(polygon%edges(i, 1)) + 1 
            counter(polygon%edges(i, 2)) = counter(polygon%edges(i, 2)) + 1 
        end do 

        ! Set logical
        polygon%issimple = all(counter <= 2)

        ! Housekeeping
        !=============
        deallocate(counter) 



    end subroutine
        

    !------------------------------------------------------------------!
    !                              Auxiliary                           !
    !------------------------------------------------------------------!

    ! Error handler
    subroutine PolygonErrorHandler(msg) 

        ! Description
        !============
        ! Error handler. Will print the message in 'msg' and stop the
        ! execution of the program. 

        ! Declare variables
        !==================
        ! Arguments
        character(*), intent(in)        :: msg 

        ! Handle error
        !=============
        ! Print 
        print *, 'PolygonErrorHandler: ', msg 

        ! Stop program execution
        stop 

    end subroutine

    ! Warning handler
    subroutine PolygonWarningHandler(msg) 

        ! Description
        !============
        ! Warning handler: will print only the message with warning in
        ! front, will not stop program execution.

        ! Declare variables
        !==================
        ! Arguments
        character(*), intent(in)        :: msg 

        ! Handle error
        !=============
        ! Print 
        print *, 'Warning: PolygonWarningHandler: ', msg 

    end subroutine

    !------------------------------------------------------------------!
    !                              Numerics                            !
    !------------------------------------------------------------------!

    ! Euclidean distance between two points
    subroutine Distance0D(d, x1, y1, x2, y2) 

        ! Description
        !============
        ! Compute the euclidean distance between two points

        ! Declare variables
        !==================
        ! Arguments
        real(R8)            :: d, x1, x2, y1, y2 

        ! Compute
        !========
        d = sqrt((x2 - x1)**2 + (y2 - y1)**2)

    end subroutine

    ! Euclidean distance between two sets of points (element-wise)
    subroutine Distance1D(d, x1, y1, x2, y2) 

        ! Description
        !============
        ! Compute the euclidean distance between two points

        ! Declare variables
        !==================
        ! Arguments
        real(R8), allocatable, intent(in)  ::  x1(:), x2(:), y1(:), y2(:) 
        real(R8), allocatable              :: d(:)

        ! Compute
        !========
        d = sqrt((x2 - x1)**2 + (y2 - y1)**2)

    end subroutine

    ! Intersection between two lines with points x11, x12, x21, x22 
    subroutine LineIntersections(x, y, x11, y11, x12, y12, x21, y21, &
        x22, y22)

        ! Description
        !============
        ! Compute the intersection between two lines. Hereto, we solve
        ! the following system for x, y:
        !
        !       (y12 - y11)(x - x11) = (x12 - x11)(y - y11)
        !       (y22 - y21)(x - x21) = (x22 - x21)(y - y21)
        !
        ! The solution can be computed analytically in this case 
        ! (provided some conditions on the differences, which is checked
        ! for with disttol) 

        ! If the intersection does not exist (e.g. parallel lines up to
        ! precision disttol),  x and y are NaN valued. 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)    :: x11, y11, x12, y12, x21, y21, &
            x22, y22 
        real(R8), intent(out)   :: x, y

        ! Auxiliary
        real(R8)                :: det, r1, r2, nan, dx1, dx2, dy1, &
            dy2 

        ! Loop

        ! Initialize
        !===========
        ! Compute
        dx1 = x12 - x11 
        dy1 = y12 - y11 
        dx2 = x22 - x21 
        dy2 = y22 - y21 

        ! Check determinant
        det = -dy1*dx2 + dx1*dy2
        if ( (abs(det) < macheps) ) then 
            x = IEEE_VALUE(nan, IEEE_QUIET_NAN)
            y = x 
            return 
        end if 

        ! Compute intersection
        !=====================
        ! Hedge for small dx when computing slope 
        if (abs(dx1) > macheps) then 
            r1 = dy1/dx1 
            if (abs(dx2) > macheps) then 
                ! Two non-parallel, non-vertical and non-horizontal lines
                r2 = dy2/dx2 
                x = (r1*x11 - r2*x21 -y11 + y21)/(r1 - r2)
                y = r1*(x - x11) + y11
            else 
                ! Second line is vertical line, first one is non-vertical
                x = x21
                y = r1*(x - x11) + y11
            end if 

        else
            ! First line is vertical line, second is non-vertical 
            ! (otherwise, det would have been zero)
            x = x11 
            r2 = dy2/dx2 
            y = r2*(x - x21) + y21
        end if 

    end subroutine

    ! Intersection between segment and line 
    subroutine SegmentLineIntersecions(x, y, x11, y11, x12, y12, x21, y21, &
        x22, y22)

        ! Description
        !============
        ! This routine computes the intersection between a line segment
        ! given by the points (x11, y11) and (x12, y12) and the line
        ! given by the points (x21, y21), (x22, y22). First, it is 
        ! checked whether there is an intersection between the two lines
        ! formed by this point. If there is an intersection, then it is
        ! checked whether this intersection lies on the line segment. 
        ! The latter is done by checking the scalar product between the
        ! vectors formed by the intersection point and the segment's 
        ! points (this should be negative if the point lies on the 
        ! segment). We hedge for machine precision effects by checking
        ! if the point lies within macheps to one of the nodes. 
        
        ! If no intersection is found, NaNs are returned for x and y

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)    :: x11, y11, x12, y12, x21, y21, &
            x22, y22 
        real(R8), intent(out)   :: x, y

        ! Auxiliary
        real(R8)                :: nan, d1, d2, dotprod

        ! Compute intersection
        !=====================
        call LineIntersections(x, y, x11, y11, x12, y12, x21, y21, &
            x22, y22)

        ! Return if no intersection is found - x and y will be NaN 
        ! already
        if (isnan(x)) then 
            return 
        end if 

        ! Check if intersection is on segment
        dotprod = (x - x11)*(x - x12) + (y - y11)*(y - y12)
        call Distance(d1, x, y, x11, y11) 
        call Distance(d2, x, y, x12, y12)
        if ( (d1 < disttol) .or. (d2 < disttol) ) then 
            return 
        end if 
        if (dotprod > 0) then 
            x = IEEE_VALUE(nan, IEEE_QUIET_NAN)
            y = x 
        end if 

    end subroutine

    ! Intersection between two segments
    subroutine SegmentIntersecions(x, y, x11, y11, x12, y12, x21, y21, &
        x22, y22)

        ! Description
        !============
        ! This routine computes the intersection between two line segments
        ! given by the points (x11, y11) and (x12, y12) and 
        ! (x21, y21), (x22, y22). First, it is 
        ! checked whether there is an intersection between the two lines
        ! formed by these points. If there is an intersection, then it is
        ! checked whether this intersection lies on both line segments. 
        ! The latter is done by checking the scalar product between the
        ! vectors formed by the intersection point and the segment's 
        ! points (this should be negative if the point lies on the 
        ! segment). We hedge for machine precision effects by checking
        ! if the point lies within macheps to one of the nodes. 
        
        ! If no intersection is found, NaNs are returned for x and y

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)    :: x11, y11, x12, y12, x21, y21, &
            x22, y22 
        real(R8), intent(out)   :: x, y

        ! Auxiliary
        real(R8)                :: nan, d1, d2, dotprod

        ! Compute intersection
        !=====================
        call LineIntersections(x, y, x11, y11, x12, y12, x21, y21, &
            x22, y22)

        ! Return if no intersection is found - x and y will be NaN 
        ! already
        if (isnan(x)) then 
            return 
        end if 

        ! Check if intersection is on segment 1
        dotprod = (x - x11)*(x - x12) + (y - y11)*(y - y12)
        call Distance(d1, x, y, x11, y11) 
        call Distance(d2, x, y, x12, y12)
        if ( (d1 < disttol) .or. (d2 < disttol) ) then 
            return 
        end if 
        if (dotprod > 0) then 
            x = IEEE_VALUE(nan, IEEE_QUIET_NAN)
            y = x 
            return 
        end if 
    
        ! Check if intersection is on segment 2
        dotprod = (x - x21)*(x - x22) + (y - y21)*(y - y22)
        call Distance(d1, x, y, x21, y21) 
        call Distance(d2, x, y, x22, y22)
        if ( (d1 < disttol) .or. (d2 < disttol) ) then 
            return 
        end if 
        if (dotprod > 0) then 
            x = IEEE_VALUE(nan, IEEE_QUIET_NAN)
            y = x 
            return 
        end if 

    end subroutine

    ! Intersections between a segment and a polygon
    !subroutine SegmentPolygonIntersections(x, y, x1, x2, xpe, ype)

    !end subroutine

    ! Inpolygon routine
    subroutine Inpolygon(polygon, xq, yq, in)

        ! Description
        !============
        ! Compute whether the query points xq, yq are located in the
        ! polygon or not. To compute this, we apply a winding number 
        ! algorithm (often referred to as Sunday's algorithm, but here 
        ! we based ourselves on the paper "A winding number and 
        ! point-in-polygon algorithm" by David G. Alciatore and Rick 
        ! Miranda. Basically, this consists of a ray tracing algorithm,
        ! where a ray is traced along the x-axis from the origin (the 
        ! polygon is shifted by each query point), and it is checked
        ! whether a polygon edge crosses the x-axis at the positive 
        ! side. If it goes from positive to negative, the winding number
        ! decreases by a half, if it goes from negative to positive y, 
        ! the winding number increases by a half. Points that lie 
        ! in the interior of the polygon should have a winding number
        ! different than zero. If they lay outside, it is equal to zero.

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT), intent(in)       :: polygon 
        real(R8), allocatable, intent(in)   :: xq(:), yq(:) 
        logical, allocatable, intent(out)   :: in(:)

        ! Auxiliary 
        integer(I8)                         :: nq, np 
        integer(I8), allocatable            :: w(:), pvert(:) 

        real(R8)                            :: xi
        real(R8), allocatable               :: xp(:), yp(:), xcross(:), &
            tempx(:), tempy(:) 

        ! Loop
        integer(I8)                         :: i, j 

        ! Initialize
        !===========
        ! Associate
        associate(&
            x           => polygon%x,           &
            y           => polygon%y,           &
            vert        => polygon%vert,        &
            isclosed    => polygon%isclosed,    &
            ne          => polygon%ne) 

        ! Number of query points
        nq = size(xq)

        ! Winding number*2 (to be able to store as integer)
        allocate(w(nq), in(nq))
        w(:) =  0
        in(:) = .true.

        ! Shifted polygon points
        if (isclosed) then 
            allocate(xp(ne+1), yp(ne+1), pvert(ne+1))
            pvert = vert
        else
            allocate(xp(ne+2), yp(ne+2), pvert(ne+2))
            pvert = [vert, vert(1)]
        end if
        np = size(xp)
        allocate(xcross(np-1), tempx(np), tempy(np))
        tempx = x(pvert)
        tempy = y(pvert)
        
        ! Loop
        !=====
        do i = 1, nq 
            ! Shift points
            xp = tempx - xq(i)
            yp = tempy - yq(i)

            ! Compute which points cross x-axis
            xcross = yp(1:np-1)*yp(2:np)

            do j = 1, np-1
                if (xcross(j) < 0) then 
                    ! x-coordinate of intersection
                    xi = xp(j) + yp(j)*(xp(j+1) - xp(j))/(yp(j+1) - yp(j))

                    if (xi > 0) then 
                        if (yp(j) <0) then 
                            w(i) = w(i) + 2
                        else 
                            w(i) = w(i) - 2 
                        end if 
                    end if 
                elseif ( (yp(j) == 0) .and. (xp(j+1) > 0) ) then
                    if ( (yp(j+1) > 0) ) then 
                        w(i) = w(i) + 1
                    else 
                        w(i) = w(i) - 1
                    end if 
                elseif ( (yp(j+1) == 0) .and. (xp(j+1) > 0) ) then 
                    if (yp(j) < 0) then 
                        w(i) = w(i) + 1
                    else 
                        w(i) = w(i) - 1
                    end if 
                end if 

            end do

        end do 

        ! Check in/out
        !=============
        where (w == 0) in = .false.

        ! Housekeeping
        !=============
        deallocate(xp, yp, w, xcross)
        end associate 



    end subroutine 

    !------------------------------------------------------------------!
    !                                Tests                             !
    !------------------------------------------------------------------!

    ! Main test routine
    subroutine TestPolygonRoutines()

        ! Description
        !============
        ! Run tests on some polygons

        ! Declare variables
        !==================
        ! Arguments

        ! Auxiliary
        type(PolygonSetUDT)         :: ps 

        integer(I8)                 :: id 
        real(R8), allocatable       :: x(:), y(:) 


        ! Set polygon
        !============
        ! Polygon case
        id = 1

        ! Coordinates
        call ConstructTestPolygonCoordinates(id, x, y)

        ! Test
        !=====
        call ps%Construct(x, y)

        !call Plot2DPolygon(ps%polygons(1)%x, ps%polygons(1)%y, size(ps%polygons(1)%x), '-p')
        !call Quiverplot2D(ps%polygons(1)%ex, ps%polygons(1)%ey, &
        !    ps%polygons(1)%nx, ps%polygons(1)%ny, &
        !    size(ps%polygons(1)%ex), '-p')

        ! Set polygon
        !============
        ! Polygon case
        id = 5

        ! Coordinates
        call ConstructTestPolygonCoordinates(id, x, y)

        ! Test
        !=====
        if (allocated(ps%polygons)) then 
            deallocate(ps%polygons)
        end if
        call ps%Construct(x, y)

    end subroutine

    ! Set of test polygons
    subroutine ConstructTestPolygonCoordinates(id, x, y) 

        ! Description
        !============
        ! Return the test polygon coordinates belonging to index id

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)         :: id 
        real(R8), allocatable           :: x(:), y(:)

        ! Auxiliary
        integer(I8)                     :: np

        ! Construct polygon coordinates
        !==============================
        ! Check
        if (allocated(x)) then 
            deallocate(x) 
        end if 
        if (allocated(y)) then 
            deallocate(y) 
        end if 

        ! Build
        select case (id) 

            case (1) 

                ! Normal polygon
                np = 4
                allocate(x(np), y(np))
                x = [1, 1, 2, 4]
                y = [0, 2, 2, 1]
            
            case (5)

                ! Polygon with some duplicate points
                np = 12
                allocate(x(np), y(np))
                x = [1, 1, 1, 1, 2, 4, 1, 1, 2, 5, 6, 8]
                y = [1, 1, 0, 2, 2, 1, 1, 1, 3, 2, 0, 2]

            case default 

            ! Throw error
            call PolygonErrorHandler('No coordinates found for this id')
            
        end select 




    end subroutine



end module