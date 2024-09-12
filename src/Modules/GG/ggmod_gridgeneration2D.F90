!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module provides functionality to generate 2D grids (potentially in 
! different ways) starting from a topological mesh. The main drivers 
! are made public and can be found under the 'drivers' subroutines tab. 

! The general idea behind the 2D grid generator for aligned,
! unstructured meshes is the following:
!
! - The topological mesh divides the entire region into parts with 
!   distinct flux tubes, which are defined as areas bounded by multiple
!   aligned boundaries (points or faces) and at least two radial faces (
!   in some cases, this may be the same face, e.g. in core regions). 
!   True disc regions with no radial lines are not allowed. 
! - Each flux tube consists of cells, which basically have the same 
!   properties as the flux tube, except that there are exactly two 
!   radial lines. 
! - Each flux tube has a set of field lines along which the grid vertices
!   should be distributed (if we go for non-aligned grids, this is even
!   not a requirement anymore and things become even easier actually). 
!   Each cell that belongs to a certain tube therefore has the same
!   field lines, but every time a different part of that field line. 
! - To grid a cell, we therefore need to distribute the vertices along
!   all the lines (its aligned boundaries and the field lines) and make
!   sure that vertex distributions on the boundaries are propagated to
!   all different cells
! 
! Note that we assume that all cells are present in flux tubes. 
! Below we mention the different algorithms for distributing the 
! vertices along the field lines. Note that the amount of field lines
! and the amount of vertices along a field line can be influenced 
! locally by setting the properties of the poloidal and radial 
! vertex distributors (see also ggmod_Vertexdistribution2D). 
! 
! Forming of cells and vertices can be done in different ways as well,
! though the most important one is probably the 'quads_triangles', since
! it will promote the use of quads wherever it is possible and leads to
! the most orthogonal face with the field. 

! Algorithms
!===========
! Independent
!------------
! This algorithm distributes the vertices for each cell independently. 
! The vertex distribution on the entire cell boundary is assumed to be 
! fixed, which is the only reason why each cell can be dealt with 
! independently. Typically leads to lower quality grids, since 
! distributions may not be propagated decently, and since no effort for
! orthogonality is made. 

! Orthogonal
!-----------
! This algorithm attempts to distribute the vertices as orthogonal as
! possible by tracing lines orthogonal to the field (so-called radial 
! lines). By computing the intersection between radial and poloidal 
! (aligned) lines, we find the required nodes. To do this in a consistent
! way without too many sudden changes near boundaries, we need to 
! propagate these radial lines properly over the cell boundaries, which 
! is non-trivial. To this purpose, we need to carefully choose which 
! cells to grid when. The main rules are as follows:
! - We always distribute vertices from high field to low field values.
!   This defines the direction in which to trace the radial lines, and 
!   defines a direction in which to choose the next cell.
! - A cell can be gridded if:
!       all its high field poloidal boundaries are actual boundary faces
!       i.e. they do not connect to another cell
!   OR
!       all its high field, non-boundary poloidal faces have been gridded
!       previously and therefore their node distribution is known
!   Normally, this should lead to all cells being gridded and to the 
!   radial lines being properly propagated. 
! - Note that it is very much encouraged to do refinement/coarsening
!   in this algorithm, since radial line bunching is often a thing.

module ggmod_gridgeneration2D

    ! Load modules
    use mod_precision
    use mod_constants
    use mod_errorhandler
    use mod_dynamicarrays
    use mod_contour2D
    use mod_polygon
    use mod_sort
    use mod_definitions, only: TMvertexbndID, TMvertexmaxID, &
        TMvertexminID, TMvertexsaddleID, TMvertextp1ID, &
        TMvertextp2ID, TMfacepolID, TMfaceradID, TMfacebndID, &
        TMvertexsplitID, TMfacesepID
    use mod_linearsolverinterface, only: SolveDenseLinearSystemDI
    use goatmod_types, only : magneticFieldUDT, VesselUDT, GridUDT
    use goatmod_userinput, only : GGoptionsUDT
    use ggmod_topology2D
    use ggmod_vertexdistribution2D
    use DistributionFunction
    use mod_plotter
    use mod_utility, only: wall_time
    implicit none
    private 
    public :: GenerateUnstructuredAlignedGrid

    ! Module parameters
    real(R8), parameter, private        :: tprelfieldtol = 1e-10 ! relative field tolerance under which extrema are removed
    real(R8), parameter, private        :: disttol = 1e-12 ! distance tolerance

    !==================================================================!
    !                                                                  !
    !                            TYPES                                 !
    !                                                                  !
    !==================================================================!

    ! Grid data related to topological mesh
    !======================================
    ! Vertex data
    type :: GGTMVertexDataUDT

    end type 
    
    ! Face data
    type :: GGTMFaceDataUDT

        ! Description
        !============
        ! Contains additional face data, including the grid vertices
        ! xv, yv 
        real(R8), allocatable, dimension(:)     :: xv, yv
        integer(I8), allocatable, dimension(:)  :: vert

    end type 

    ! Field line data
    type :: GGTMFieldlineDataUDT

        real(R8), allocatable, dimension(:)         :: xv, yv, xl, yl
        integer(I8), allocatable, dimension(:)      :: vert, facelabels
        integer(I8)                                 :: fsID

    end type

    ! Cell data
    type :: GGTMCellDataUDT

        ! Description
        !============
        ! Contains a lot of additional data for the gridder, including
        ! the field line data, vertices distributed thereon, etc. Note 
        ! that the field line data is stored in an auxiliary structure
        ! 'lines' for each cell. Additionally, the face indices that 
        ! are on the high field side are stored in hffaces, the low field
        ! in lffaces (similar for vertices of the cell). These arrays
        ! are not sorted yet in any particular order (will have to check
        ! face orientation anyway when making the grid)

        type(GGTMFieldlineDataUDT), allocatable     :: lines(:)
        type(GGTMFieldlineDataUDT)                  :: hfline, lfline
        integer(I8)                                 :: srflabel, erflabel, &
            srf, erf, region
        integer(I8), allocatable, dimension(:)      :: srfvert, erfvert, &
            hffaces, lffaces, hfvert, lfvert
        logical                                     :: flipsrf, fliperf

    end type

    ! Tube data
    type :: GGTMTubeDataUDT

        ! Description
        !============
        ! Contains additional flux tube data, such as field value 
        ! data etc 

        real(R8), allocatable, dimension(:)     :: fval
        integer(I8)                             :: distributionface

    end type 


    ! Topological mesh grid generator data
    type :: GGTMDataUDT

        ! Description
        !============
        ! This type contains additional information related to the 
        ! topological mesh to construct a grid. 

        type(GGTMVertexDataUDT), allocatable, dimension(:)  :: vert
        type(GGTMFaceDataUDT), allocatable, dimension(:)    :: face 
        type(GGTMCellDataUDT), allocatable, dimension(:)    :: cell 
        type(GGTMTubeDataUDT), allocatable, dimension(:)    :: tube 

    contains 

        ! Initializer
        procedure :: Initialize         => InitializeGGTMData

    end type 

    ! Temporary grid type
    !====================
    ! Vertices
    type :: GGVertUDT 
        
        ! Description
        !============
        ! Type that contains only the basic vertex structures necessary
        ! to define the grid. Here, these are only the x, y coordinates 
        ! (the IDs are assumed to be equal to the element number of the
        ! coordinate vector). For optimized handling, the types are set 
        ! as allocatable dynamic array classes, such that suitable 
        ! classes can be determined during initialization.

        integer(I8)                                 :: ntot 
        class(RealDynamicArrayUDT), allocatable     :: x, y
        class(IntegerDynamicArrayUDT), allocatable  :: fieldlineID

    contains

        ! Initialization
        procedure :: Initialize         => InitializeGGVert

    end type

    ! Faces
    type :: GGFaceUDT 
        
        ! Description
        !============
        ! Type that contains only the basic face structures necessary
        ! to define the grid. Here, these are:
        ! - the face vertices v1, v2
        ! - the face label  
    
        ! For optimized handling, the types are set 
        ! as allocatable dynamic array classes, such that suitable 
        ! classes can be determined during initialization.

        integer(I8)                                 :: ntot 
        class(IntegerDynamicArrayUDT), allocatable  :: v1, v2, label, &
            region

    contains

        ! Initialization
        procedure :: Initialize         => InitializeGGFace

    end type

    ! Cells
    type :: GGCellUDT 
        
        ! Description
        !============
        ! Type that contains only the basic cell structures necessary
        ! to define the grid. Here, these are:
        ! - the cell vertices 
        ! - the cell vertex pointer
        ! - the cell region
    
        ! For optimized handling, the types are set 
        ! as allocatable dynamic array classes, such that suitable 
        ! classes can be determined during initialization.

        integer(I8)                                 :: ntot 
        class(IntegerDynamicArrayUDT), allocatable  :: vert, vp1, vp2, &
            region

    contains

        ! Initialization
        procedure :: Initialize         => InitializeGGCell

    end type

    ! Grid
    type :: GGGridUDT

        type(GGVertUDT)             :: vert 
        type(GGFaceUDT)             :: face 
        type(GGCellUDT)             :: cell 

    contains 

        ! Initialization
        procedure :: Initialize         => InitializeGGGrid

        ! Vertex addition
        procedure :: AddVert            => AddGGVert

        ! Face addition
        procedure :: AddFace            => AddGGFace 

        ! Cell addition
        procedure :: AddCell            => AddGGCell
    end type 



    contains 

    !==================================================================!
    !                                                                  !
    !                           ROUTINES                               !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                     GRID GENERATION DRIVERS                      !
    !------------------------------------------------------------------!

    ! Unstructured aligned grid generator
    subroutine GenerateUnstructuredAlignedGrid(topomesh, magneticField, &
        vessel, fieldtracer, boundarytracer, options)

        ! Description
        !============
        ! Generate a grid from scratch for a topological mesh given in topomesh. It
        ! is assumed that each cell in the topological mesh has exactly two
        ! 'radial' boundaries (type 1) or one radial and one vessel boundary (type
        ! 3), and an arbitrary number of poloidal boundaries. First, we construct
        ! based on this topology mesh the tubes (list of cell and faces) through
        ! which one bundle of field coordinates should pass. Then, we distribute
        ! vertices on all faces of the topological mesh, where we assume that for
        ! each 'poloidal' face, we have the freedom to distribute however we like
        ! regardless of the other faces, but that the distribution of the radial
        ! faces is determined based on the tubes. Afterwards, we trace the field
        ! lines of each tube. Then, we are ready to distribute for each grid cell
        ! the grid nodes, which should be based on a given distribution (determined
        ! through the options). Finally, we connect all cells and vertices together
        ! and call all other necessary grid initialization routines. 

        ! Notes
        !======
        ! Note 1: all vertices of the topological mesh will become vertices of the
        ! grid

        ! Note 2: vertex distributions of each face are added on the topological
        ! mesh level (facedata)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)          :: topomesh 
        type(MagneticFieldUDT)      :: magneticField
        type(VesselUDT)             :: vessel 
        class(ContourTracerUDT)     :: fieldtracer, boundarytracer 
        type(GGoptionsUDT)          :: options 

        ! Auxiliary
        type(GGTMDataUDT)           :: ggtmdata 
        class(VertexDistributor2DUDT), allocatable      :: &
            poloidalvertexdistributor, radialvertexdistributor
        class(DistributionFunctionUDT), allocatable     :: & 
            magneticFieldDF 
        type(GGGridUDT)             :: grid 

        ! Initialize
        !===========
        ! Required data of topomesh for grid generator
        call ggtmdata%Initialize(topomesh)

        ! Magnetic field distribution function
        magneticFieldDF = ConstructStructured2DDF(magneticField%interp)

        ! Temporary grid structure
        call grid%Initialize('standard')

        ! Set up vertex distribution
        !===========================
        ! Poloidal vertex distributor
        select case (options%vdptype)

        case ('uniform')

            ! Construct uniform distributor with facelength 'options%vdpdfacelength'
            poloidalvertexdistributor = ConstructUniformVertexDistributor(&
                options%vdpdfacelength, options%vdrdfieldwidth)

        case default 

            ! Unknown option
            call gdErrorHandler('GenerateUnstructuredAlignedGrid: ' // & 
                'poloidal vertex distribution option: ' // options%vdptype // &
                ' not implemented')

        end select

        ! Radial vertex distributor
        select case (options%vdrtype)

        case ('uniform')

            ! Construct uniform distributor 
            radialvertexdistributor = ConstructUniformVertexDistributor(&
                options%vdpdfacelength, options%vdrdfieldwidth)

        case Default

            ! Unknown option
            call gdErrorHandler('GenerateUnstructuredAlignedGrid: ' // & 
                'radial vertex distribution option: ' // options%vdrtype // &
                ' not implemented')

        end select

        !pdoptions = options.poloidaldistributor;
        !pdoptions = SetPoloidalDistributorOptions(pdoptions);
        !pdoptions.distribution.distrfield = field;
        !pdoptions.distribution.bnd = bnd;
        !poloidaldistributor = ConstructVertexDistributor(pdoptions);

        ! Extract x-point data
        !isxp = cat(1, topomesh.vertdata.type) == 2;
        !xpdata = cat(1, topomesh.vertdata(isxp).F);

        ! Radial
        !rdoptions = options.radialdistributor;
        !rdoptions = SetRadialDistributorOptions(rdoptions);
        !rdoptions.distribution.distrfield = field;
        !rdoptions.distribution.bnd = bnd;
        !rdoptions.distribution.xpointdata = xpdata;
        !radialdistributor = ConstructFieldDistributor(rdoptions, field);

        ! Distribute vertices on topological faces
        !=========================================
        ! Poloidal faces
        call DistributeVerticesTopologicalMeshFaces(ggtmdata, topomesh, &
            poloidalvertexdistributor, [TMfacepolID, TMfacesepID])

        ! Relevant radial faces of tubes
        call DistributeVerticesTopologicalMeshTubes(ggtmdata, &
            topomesh, radialvertexdistributor, magneticFieldDF, [TMfaceradID, TMfacebndID])

        ! Generate initial grid
        !======================
        ! Construct the required topological cell data
        call AddTopologicalMeshCellGriddingData(ggtmdata, topomesh, &
            fieldtracer, magneticField)

        ! Distribute vertices 
        select case (options%ggmethod)

        case ('independent')

            call DistributeVerticesIndependent(ggtmdata, topomesh, grid, &
                poloidalvertexdistributor)

        case default 

            call gdErrorHandler('GenerateUnstructuredAlignedGrid: ' // & 
                'unknown distribution method: ' // options%ggmethod)

        end select

        ! Write data
        call WriteGGTMData(ggtmdata, 'ggtmdata_after_vertexdistribution')

        ! Extract vertices
        call ConstructGridVertices(ggtmdata, grid, topomesh)

        ! Construct faces and cells
        select case (options%cellconstructionmethod)

        case ('quads_triangles')

            ! Evaluate the magnetic field vector at vertex locations
            call ConstructCellsQuadTria(ggtmdata, grid, magneticField)

        case default

            call gdErrorHandler('GenerateUnstructuredAlignedGrid: ' // & 
                'unknown cell construction method: ' // options%cellconstructionmethod)

        end select

        ! Add interconnections
        ![grid] = ComputeGridInterconnections(grid);

        ! Process grid
        !=============

        ! Visualize
        !----------
        !VisualizeGrid(grid, 1);


    end subroutine

    !------------------------------------------------------------------!
    !                         GRID GENERATION                          !
    !------------------------------------------------------------------!

    ! Independent gridder
    subroutine DistributeVerticesIndependent(ggtmdata, &
        topomesh, grid, vd)

        ! Description
        !============
        ! This routine generates for each topological cell a grid in 
        ! an independent way. It is assumed that grid coordinates are
        ! already present on all topological mesh faces, but that no 
        ! vertex indices have been assigned. This is done as a first
        ! step in this routine. Afterwards, vertices are distributed
        ! per grid cell over the lines defined there. 

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMDataUDT)                          :: ggtmdata 
        class(TopomeshUDT), intent(in)              :: topomesh 
        class(GGGridUDT), intent(inout)             :: grid 
        class(VertexDistributor2DUDT), intent(in)   :: vd

        ! Auxiliary
        integer(I8)                                 :: nv, vertID
        integer(I8), allocatable, dimension(:)      :: tvID, srfvID, &
            erfvID

        ! Loop
        integer(I8)                                 :: i, j, k

        ! Initialize
        !===========
        ! Associate
        associate(&
            vert            => topomesh%vert,   &
            face            => topomesh%face,   &
            facedata        => ggtmdata%face,   &
            cell            => topomesh%cell,   &
            celldata        => ggtmdata%cell    &
            )

        ! Add topological mesh vertices
        !==============================
        ! Only update counter, vertices are added afterwards
        grid%vert%ntot = grid%vert%ntot + vert%ntot

        ! Set vertID
        vertID = grid%vert%ntot

        ! Add face vertices
        !==================
        do i = 1, face%ntot
            ! Compute number of new vertices
            nv = size(facedata(i)%xv) - 2 

            ! Set ID
            tvID = [face%vert(i, 1), (k, k = vertID+1, vertID+nv), face%vert(i, 2)]
            facedata(i)%vert = tvID

            ! Update 
            vertID = vertID + nv
            grid%vert%ntot = grid%vert%ntot + nv
        end do

        ! Add cell vertices
        !==================
        do i = 1, cell%ntot

            ! Get cell starting and ending radial line vertices
            srfvID = facedata(celldata(i)%srf)%vert
            if (celldata(i)%flipsrf) then 
                srfvID = srfvID(size(srfvID):1:-1)
            end if 
            erfvID = facedata(celldata(i)%erf)%vert
            if (celldata(i)%fliperf) then 
                erfvID = erfvID(size(erfvID):1:-1)
            end if 

            ! Distribute lines
            do j = 1, size(celldata(i)%lines)
                ! Distribute over line
                call vd%DistributeOverCurve(celldata(i)%lines(j)%xl, &
                celldata(i)%lines(j)%yl, celldata(i)%lines(j)%xv, &
                celldata(i)%lines(j)%yv, nv)

                ! Set vertex ID
                celldata(i)%lines(j)%vert = [srfvID(j+1), &
                    & (k, k = vertID+1, vertID+nv-2), erfvID(j+1)]

                ! Update total number of vertices
                vertID = vertID+nv-2
                grid%vert%ntot = grid%vert%ntot + nv-2
            end do 

            ! Extract high field line
            call ExtractTMCellAlignedBoundary(celldata(i), 'high', ggtmdata, &
                topomesh, celldata(i)%hfline)

            ! Extract low field line
            call ExtractTMCellAlignedBoundary(celldata(i), 'low', ggtmdata, &
                topomesh, celldata(i)%lfline) 

        end do 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Grid vertex constructor
    subroutine ConstructGridVertices(ggtmdata, grid, topomesh)

        ! Description
        !============
        ! Construct grid vertices by extracting them from the
        ! topological mesh cells. It is assumed that the total
        ! number of grid vertices is known and given in grid%vert%ntot
        ! but that the vertex coordinates have not yet been added. 
        ! Additionally, we do some checks to ensure all vertices have
        ! been accounted for 

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMDataUDT)                      :: ggtmdata 
        class(GGGridUDT), intent(inout)         :: grid 
        class(TopomeshUDT), intent(in)          :: topomesh 

        ! Auxiliary
        real(R8), allocatable, dimension(:)     :: xv, yv 
        integer(I8), allocatable, dimension(:)  :: fieldlineID, tvID
        logical, allocatable, dimension(:)      :: isvertfound

        ! Loop
        integer(I8)                             :: i, j

        ! Initialize
        !===========
        ! Associate
        associate(&
            cell            => topomesh%cell,   &
            celldata        => ggtmdata%cell    &
            )
        ! Allocate
        allocate(xv(grid%vert%ntot), yv(grid%vert%ntot), &
            isvertfound(grid%vert%ntot), fieldlineID(grid%vert%ntot))
        isvertfound = .false. 

        ! Determine vertices
        !===================
        do i = 1, cell%ntot 
            ! Loop over lines
            do j = 1, size(celldata(i)%lines)
                ! Get IDs
                tvID = celldata(i)%lines(j)%vert

                ! Set logicals
                isvertfound(tvID) = .true.

                ! Add coordinates
                xv(tvID) = celldata(i)%lines(j)%xv 
                yv(tvID) = celldata(i)%lines(j)%yv 
                
                ! Set field line ID
                fieldlineID(tvID) = celldata(i)%lines(j)%fsID
            end do

            ! Add high field and low field line
            ! Get IDs
            tvID = celldata(i)%hfline%vert

            ! Set logicals
            isvertfound(tvID) = .true.

            ! Add coordinates
            xv(tvID) = celldata(i)%hfline%xv 
            yv(tvID) = celldata(i)%hfline%yv 
            
            ! Set field line ID
            fieldlineID(tvID) = celldata(i)%hfline%fsID

            ! Get IDs
            tvID = celldata(i)%lfline%vert

            ! Set logicals
            isvertfound(tvID) = .true.

            ! Add coordinates
            xv(tvID) = celldata(i)%lfline%xv 
            yv(tvID) = celldata(i)%lfline%yv 
            
            ! Set field line ID
            fieldlineID(tvID) = celldata(i)%lfline%fsID
        end do 

        ! Checks
        !=======
        ! This is normally impossible
        if (.not. all(isvertfound)) then 
            call gdErrorHandler('ConstructGridVertices: some vertices ' // & 
                'were not found, this is probably a bug')
        end if 

        ! Add
        !====
        call grid%AddVert(xv, yv, fieldlineID)

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Grid cell and face construction, quads & triangles
    subroutine ConstructCellsQuadTria(ggtmdata, grid, magneticField)

        ! Description
        !============
        ! Construct (mostly) quadrilateral cells from the vertex distributions given in lines.
        ! To this end, we assume throughout most of this routine (and possibly
        ! subroutines) that the lines are nearly 'parallel' to each other, and that
        ! the local curvature doesn't vary quickly along the lines. 
        ! To this end, we do the following steps:

        ! 1) Start at the beginning of a pair of subsequent lines
        ! 2) Check which additional face should be added (take one which results in
        ! smallest aligned face length for now). Then, add the additional faces and
        ! cells to the grid using the dedicated routines for that.
        ! 3) Repeat for all line pairs until ended

        ! This is the basic algorithm. In practice, we added the following improved
        ! logic:
        ! - We check if the resulting face would overlap/intersect with one of the
        ! lines by checking vector cross products
        ! - We check if the first face intersects with any of the lines (not in the
        ! initial points). If that's the case, we know that an overlap will be
        ! created, but we try to construct the grid as such that the overlapping
        ! triangles can be removed by the 'RemoveNarrowBoundaryTriangles'
        ! routine. This is only possible if at least a subset of the faces is not
        ! overlapping. For this, we check potential faces of both start nodes in
        ! the two lines with the nodes in the other lines. The one with the least
        ! amount of intersections is taken. 


        ! Notes
        !======
        ! Note 1: we need to check if the last face is the same as the initial one
        ! to hedge for closed cells. 

        ! Note 2: no guarantees can be given on non-overlapping cells. 

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMDataUDT)                      :: ggtmdata
        class(GGGridUDT), intent(inout)         :: grid 
        type(MagneticFieldUDT), intent(in)      :: magneticField

        ! Auxiliary
        real(R8)                                :: dx1, dy1, &
            alpha1, dx2, dy2, alpha2, dx3, dy3, alpha3, bxf(1:3), byf(1:3)
        real(R8), allocatable, dimension(:)     :: xint, yint
        integer(I8)                             :: nf, nc, ncv, mink1, &
            mink2, v1, v3, v2, v4, tff(1:2), indmin
        integer(I8), allocatable, dimension(:)  :: tempfacelabels, &
            tempcellvert, s1, s2, tempfaceregion, tempcellregion
        integer(I8), allocatable, dimension(:, :)   :: tempfacevert, &
            tempcellvertP
        logical                                 :: isintersectingl1, &
            isintersectingl2, doquad, issameface
        type(GGTMFieldlineDataUDT), allocatable :: lines(:)
        type(PolygonUDT)                        :: facep, lp

        ! Loop
        integer(I8)                             :: i, j, k1, k2

        ! Initialize
        !===========
        associate(celldata  => ggtmdata%cell)

        ! Loop
        !=====
        do i = 1, size(celldata)
            ! Concatenate lines for ease
            lines = [celldata(i)%hfline, celldata(i)%lines, &
                celldata(i)%lfline]

            ! Loop over all lines
            do j = 1, size(lines)-1
                ! Unpack
                associate(&
                    l1      => lines(j),    &
                    l2      => lines(j+1),      &
                    n1      => size(lines(j)%xv),   &
                    n2      => size(lines(j+1)%xv), &
                    sff     => celldata(i)%srflabel, &
                    eff     => celldata(i)%erflabel)

                ! Initialize
                k1 = 1
                k2 = 1
                nf = 0
                nc = 0
                ncv = 0
                allocate(tempfacevert(4*(n1+n2), 2), tempfacelabels(4*(n1+n2)), &
                    tempcellvert(3*(n1+n2)), tempcellvertP(4*(n1+n2), 2)) ! overestimations

                ! Check if the first face intersects with one of the lines
                isintersectingl1 = .false. 
                isintersectingl2 = .false.
                mink1 = 0
                mink2 = 0
                call facep%Construct([l1%xv(1), l2%xv(1)], [l1%yv(1), l2%yv(1)])

                ! First line
                if (l1%vert(1) == l1%vert(size(l1%vert))) then 
                    ! Closed 
                    call lp%Construct(l1%xv(2:size(l1%xv)-1), l1%yv(2:size(l1%xv)-1)) ! no need to include 
                    call PolygonIntersections(lp, facep, xint, yint, s1, s2)
                else
                    call lp%Construct(l1%xv(2:), l1%yv(2:))
                    call PolygonIntersections(lp, facep, xint, yint, s1, s2)
                end if
                if (size(xint) > 0) then 
                    isintersectingl1 = .true.
                    mink1 = s1(1)+2
                end if

                ! Second line
                if (l2%vert(1) == l2%vert(size(l2%vert))) then 
                    ! Closed 
                    call lp%Construct(l2%xv(2:size(l2%xv)-1), l2%yv(2:size(l2%xv)-1)) ! no need to include 
                    call PolygonIntersections(lp, facep, xint, yint, s1, s2)
                else
                    call lp%Construct(l1%xv(2:), l1%yv(2:))
                    call PolygonIntersections(lp, facep, xint, yint, s1, s2)
                end if 
                if (size(xint) > 0) then 
                    isintersectingl2 = .true.
                    mink2 = s1(1)+2
                end if
                
                ! Add the first face
                nf = nf + 1
                tempfacevert(nf, :) = [l1%vert(k1), l2%vert(k2)]
                tempfacelabels(nf) = sff

                ! Loop
                do while (.true.)

                    ! Make candidate faces
                    !---------------------
                    ! Face pair 1: vertex k1 and vertex k2+1, vertex k2, k2+1
                    ! (triangle)
                    
                    ! Face pair 2: vertex k1+1 and vertex k2, vertex k1, k1+1
                    ! (triangle)
                    
                    ! Face pair 3: k1+1, k2+1 (quad)

                    ! Determine which face to take
                    !-----------------------------
                    ! First two vertices are always the same
                    v1 = l1%vert(k1)
                    v3 = l2%vert(k2)
                    
                    ! Here, purely based on length of added aligned face
                    doquad = .false.
                    if ((k1 < n1) .and. (k2 < n2)) then 
                                
                        ! Compute angle of face normal with magnetic field
                        dx1 = l2%xv(k2+1) - l1%xv(k1)
                        dy1 = l2%yv(k2+1) - l1%yv(k1)
                        dx2 = l1%xv(k1+1) - l2%xv(k2)
                        dy2 = l1%yv(k1+1) - l2%yv(k2)
                        dx3 = l2%xv(k2+1) - l1%xv(k1+1)
                        dy3 = l2%yv(k2+1) - l1%yv(k1+1)

                        call magneticField%interp%Evaluate(&
                            [l1%xv(k1)+0.5*dx1, l2%xv(k2)+0.5*dx2, l1%xv(k1+1)+0.5*dx3], &
                            [l1%yv(k1)+0.5*dy1, l2%yv(k2)+0.5*dy2, l1%yv(k1+1)+0.5*dy3], &
                            0, 1, bxf)
                        call magneticField%interp%Evaluate(&
                            [l1%xv(k1)+0.5*dx1, l2%xv(k2)+0.5*dx2, l1%xv(k1+1)+0.5*dx3], &
                            [l1%yv(k1)+0.5*dy1, l2%yv(k2)+0.5*dy2, l1%yv(k1+1)+0.5*dy3], &
                            1, 0, byf)
                        bxf = -bxf ! adjust sign

                        alpha1 = abs(atan( (-dy1*byf(1) -dx1*bxf(1))/(-dy1*bxf(1) + dx1*byf(1))))
                        alpha2 = abs(atan( (-dy2*byf(2) -dx2*bxf(2))/(-dy2*bxf(2) + dx2*byf(2))))
                        alpha3 = abs(atan( (-dy3*byf(3) -dx3*bxf(3))/(-dy3*bxf(3) + dx3*byf(3))))
                        


                        ! Check if we still need to hedge for intersecting
                        ! lines
                        if (k1 > mink1) then 
                            isintersectingl1 = .false.
                        end if 
                        if (k2 > mink2) then 
                            isintersectingl2 = .false.
                        end if

                        ! Check which triangles are allowed
                        if (isintersectingl1) then 
                            ! First line is intersected by starting face:
                            ! triangles should be built using second line
                            alpha1 = posinfval_R8()
                            alpha3 = posinfval_R8()
                            
                        elseif (isintersectingl2) then 
                            ! Second line is intersected by starting face:
                            ! triangles should be built using first line
                            alpha2 = posinfval_R8()
                            alpha3 = posinfval_R8()
                        end if
                        
                        ! Find face that makes the smallest angle
                        indmin = minloc([alpha1, alpha2, alpha3], 1)
                        
                        ! Add the face
                        if (indmin == 3) then 
                            ! Add third face, quad
                            doquad = .true.
                            v2 = l2%vert(k2+1)
                            v4 = l1%vert(k1+1)
                            
                            ! Update counter
                            k2 = k2+1
                            k1 = k1+1
                        elseif (indmin == 1) then 
                            ! Add first face, triangle
                            v2 = l2%vert(k2+1)
                            tff = [0, l2%facelabels(k2)]
                            
                            ! Update counter
                            k2 = k2 + 1
                        elseif (indmin == 2 ) then 
                            ! Add second face, triangle
                            v2 = l1%vert(k1+1)
                            tff = [l1%facelabels(k1), 0]
                            
                            ! Update counter
                            k1 = k1 + 1
                        else
                            ! Should never happen
                            call gdErrorHandler('Something wrong')
                        end if 
                                
                        
                    elseif ((k1 == n1) .and. (k2 < n2)) then 
                        ! We have to take the first option, no vertices left in first
                        ! line
                        ! Add first face pair
                        v2 = l2%vert(k2+1)
                        tff = [0, l2%facelabels(k2)]
                        
                        ! Update counter
                        k2 = k2 + 1
                        
                        if (k2 == n2) then 
                            tff(1) = eff
                        end if 
                    elseif ((k1 < n1) .and. (k2 == n2)) then 
                        ! We have to take the second option
                        ! Add second face pair
                        v2 = l1%vert(k1+1)
                        tff = [l1%facelabels(k1), 0]
                        
                        ! Update counter
                        k1 = k1 + 1
                        
                        if (k1 == n1) then  
                            tff(2) = eff
                        end if 
                        
                    else
                        ! This shouldn't happen
                        call gdErrorHandler('Something wrong in quad gridder')
                    end if 
                        
                    ! Add grid faces and cells
                    !-------------------------
                    ! Add face pair
                    if (doquad) then 
                        nf = nf + 1
                        tempfacevert(nf, :) = [v1, v4]
                        tempfacelabels(nf) = l1%facelabels(k1-1)
                        nf = nf + 1
                        tempfacevert(nf, :) = [v3, v2]
                        tempfacelabels(nf) = l2%facelabels(k2-1)
                        nf = nf + 1
                        tempfacevert(nf, :) = [v4, v2]
                        if ((k1 /= n1) .or. (k2 /= n2)) then 
                            tempfacelabels(nf) = 0
                        else
                            tempfacelabels(nf) = eff
                        end if 
                    else
                        nf = nf + 1
                        tempfacevert(nf, :) = [v1, v2]
                        tempfacelabels(nf) = tff(1)
                        nf = nf + 1
                        tempfacevert(nf, :) = [v3, v2]
                        tempfacelabels(nf) = tff(2)
                    end if 

                    ! Add cell
                    nc = nc + 1
                    if (doquad) then 
                        tempcellvert(ncv+1:ncv+4) = [v1, v3, v2, v4]
                        tempcellvertP(nc, :) = [ncv+1, 4]
                        ncv = ncv + 4
                    else
                        tempcellvert(ncv+1:ncv+3) = [v1, v2, v3]
                        tempcellvertP(nc, :) = [ncv+1, 3]
                        ncv = ncv + 3
                    end if 
                    
                    ! Check stop criterium
                    if ((k1 >= n1) .and. (k2 >= n2)) then 
                        ! Check if the last non-aligned face was equal to the first one
                        issameface = (tempfacevert(nf-1, 1) == tempfacevert(1, 1) ) .and. &
                            (tempfacevert(nf-1, 2) == tempfacevert(1, 2))
                        issameface = issameface .or. (tempfacevert(nf-1, 2) == tempfacevert(1, 1) ) .and. &
                            (tempfacevert(nf-1, 1) == tempfacevert(1, 2))
                        if (issameface) then 
                            ! Don't add the last face
                            nf = nf-1
                        end if 
                        exit
                    end if 

                end do

                ! Housekeeping
                end associate
            end do 
    
            
            ! Add to grid
            allocate(tempcellregion(nc), tempfaceregion(nf))
            tempcellregion = celldata(i)%region
            tempfaceregion = celldata(i)%region
            call grid%AddFace(tempfacevert(1:nf, :), tempfacelabels(1:nf), tempfaceregion(1:nf))
            call grid%AddCell(tempcellvert(1:ncv), tempcellvertP(1:nc, :), tempcellregion(1:nc))

            ! Housekeeping
            deallocate(tempcellregion, tempfaceregion)

        end do 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    !------------------------------------------------------------------!
    !                  TOPOMESH GRID DATA ROUTINES                     !
    !------------------------------------------------------------------!

    ! GGTM data initialization 
    subroutine InitializeGGTMData(ggtmdata, topomesh)

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMDataUDT)                  :: ggtmdata 
        class(TopomeshUDT),intent(in)       :: topomesh
        
        ! Initialize
        !===========
        ! Substructures
        allocate(ggtmdata%vert(topomesh%vert%ntot))
        allocate(ggtmdata%face(topomesh%face%ntot))
        allocate(ggtmdata%cell(topomesh%cell%ntot))
        allocate(ggtmdata%tube(topomesh%tube%ntot))

    end subroutine 
    
    ! Vertex distribution over faces
    subroutine DistributeVerticesTopologicalMeshFaces(ggtmdata, topomesh, &
        vd, facetypes)

        ! Description
        !============
        ! Distribute vertices over topological mesh faces of the types
        ! defined in 'facetypes'. The vertices etc are stored in the
        ! ggtmdata structure. The vertex distribution is done based on 
        ! the vertexdistributor that is passed 

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMDataUDT)                      :: ggtmdata 
        class(TopomeshUDT), intent(in)          :: topomesh 
        class(VertexDistributor2DUDT), intent(in)   :: vd 
        integer(I8), intent(in)                 :: facetypes(:)

        ! Auxiliary
        integer(I8)                             :: nv 
        real(R8), allocatable, dimension(:)     :: tx, ty 

        ! Loop
        integer(I8)                             :: i 

        ! Initialize
        !===========
        ! Associate
        associate(&
            vert        => topomesh%vert,   &
            face        => topomesh%face,   &
            facedata    => ggtmdata%face    &
            )

        ! Loop over all faces and distribute
        !===================================
        do i = 1, face%ntot
            if (any(face%type(i) == facetypes)) then 
                ! Distribute
                call vd%DistributeOverCurve(face%x(i)%Get(), face%y(i)%Get(), tx, ty, nv)
                
                ! Adjust start and end to be sure
                tx(1) = vert%x(face%vert(i, 1))
                ty(1) = vert%y(face%vert(i, 1))
                tx(nv) = vert%x(face%vert(i, 2))
                ty(nv) = vert%y(face%vert(i, 2))
                
                ! Add
                facedata(i)%xv = tx
                facedata(i)%yv = ty
                !facedata%faceflags(i) = ConstructIntegerDynamicArray(spread())
                !facedata(i).faceflags = ones(numel(facedata(i).vx)-1, size(face.flags, 2)).*face.flags(i, :);
            end if 
        end do 

        ! Housekeeping
        !=============
        end associate




    end subroutine 

    ! Vertex distribution over tubes
    subroutine DistributeVerticesTopologicalMeshTubes(ggtmdata, topomesh, &
        vd, field, facetypes)

        ! Description
        !============
        ! Distribute field values over topological mesh faces of the types
        ! defined in 'facetypes'. The field values are computed per tube,
        ! first based on an initial distribution on the faces in 
        ! 'facetypes'. This chosen distribution is simply 
        ! taken as the distribution that gives the maximal amount 
        ! of faces, as this is likely the desired one. We store the 
        ! chosen topological face ID for each flux tube (we don't 
        ! propagate the distribution to other faces yet, since we have 
        ! to compute intersections anyway)

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMDataUDT)                      :: ggtmdata 
        class(TopomeshUDT), intent(in)          :: topomesh 
        class(VertexDistributor2DUDT), intent(in)    :: vd 
        class(DistributionFunctionUDT), intent(in)  :: field 
        integer(I8), intent(in)                 :: facetypes(:)

        ! Auxiliary
        integer(I8)                             :: nv, nfl, nflmax, &
            tfmax
        integer(I8), allocatable, dimension(:)  :: tf
        real(R8), allocatable, dimension(:)     :: fc, xc, yc, tx, &
            ty

        ! Loop
        integer(I8)                             :: i, j  

        ! Initialize
        !===========
        ! Associate
        associate(&
            vert        => topomesh%vert,   &
            face        => topomesh%face,   &
            facedata    => ggtmdata%face,   &
            tube        => topomesh%tube,   &
            tubedata    => ggtmdata%tube    &
            )

        ! Distribute over faces
        !======================
        do i = 1, face%ntot 
            if (any(face%type(i) == facetypes)) then 
                ! Unpack
                xc = face%x(i)%Get()
                yc = face%y(i)%Get()
                allocate(fc(face%x(i)%Size()))

                ! Distribute
                call vd%DistributeOverField(xc, yc, field, tx, ty, nv)
                tx(1) = vert%x(face%vert(i, 1))
                ty(1) = vert%y(face%vert(i, 1))
                tx(nv) = vert%x(face%vert(i, 2))
                ty(nv) = vert%y(face%vert(i, 2))

                ! Add data
                facedata(i)%xv = tx
                facedata(i)%yv = ty

                ! Housekeeping
                deallocate(fc)

            end if 
        end do 

        ! Determine tube distribution
        !============================
        do i = 1, tube%ntot 
            ! Get tube faces
            tf = tube%GetFace(i)

            ! Initialize
            nfl = 0
            nflmax = 0

            ! Loop
            do j = 1, size(tf)
                ! Determine number of field lines
                nfl = size(facedata(tf(j))%xv)
                if (nfl > nflmax) then 
                    nflmax = nfl 
                    tfmax = tf(j)
                end if 
            end do 

            ! Add maximal distribution to tube
            tubedata(i)%distributionface = tfmax

        end do 

        ! Overwrite
        !==========

        ! Housekeeping
        !=============
        end associate

    end subroutine 

    ! Cell data
    subroutine AddTopologicalMeshCellGriddingData(ggtmdata, topomesh, &
        fieldtracer, magneticField)

        ! Description
        !============
        ! This routine constructs the required cell data for 
        ! easier grid generation later on. The following data is added:
        ! - phihfaces:      poloidal face IDs with high psi value
        ! - field lines:    lines along the field with x, y coordinates
        !                   along which nodes need to be distributed.
        !                   these lines are sorted from high to low
        !                   field value and start in the starting 
        !                   and ending radial boundary of the cell.

        ! Notes
        !======
        ! Note 1: it is assumed that the tube cell and face data is
        ! properly sorted (some sanity checks are done)

        ! Note 2: in principle, only one intersection of a radial line
        ! with a poloidal line should be found. However, in some cases
        ! it may happen that there are multiple intersections or no 
        ! intersections at all. In the first case, if there are still 
        ! intersections with other radial lines, we simply take one 
        ! of the intersections. This will likely lead to some 
        ! approximation of the actual geometry (in case of boundary 
        ! faces) or approximation of the radial line. Since we do this
        ! consistently the same, it shouldn't yield issues for the 
        ! final grid as long as the resulting geometry is not too extreme.
        ! If no intersection is found with one of the radial lines of
        ! the flux tube, the flux surface is removed from the grid and
        ! not considered anymore. This may be the case if topological 
        ! regions have not been determined accurately (e.g. at 
        ! tangency points) and that therefore some contours do not 
        ! intersect some radial lines. 

        ! Note 3: for ease later on, we already add vertex IDs for the 
        ! radial face intersections with the tube. This ID starts from
        ! the topomesh vertex counter, since normally all topological
        ! mesh vertices are added to the grid. 

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMDataUDT)                      :: ggtmdata 
        class(TopomeshUDT), intent(in)          :: topomesh 
        class(ContourTracerUDT), intent(in)     :: fieldtracer 
        type(MagneticFieldUDT), intent(in)      :: magneticField

        ! Auxiliary
        integer(I8)                             :: tc, srf, erf, inderf, &
            indsrf, tf, cind, nc, ntf, incr, nv, temp, minind, maxind
        integer(I8), allocatable, dimension(:)  :: tubec, tubef, tcf, &
            tcv, tcfv1, tcfv2, hffaces, lffaces, hfvert, lfvert, &
            allIDs, s1, s2, polv, tf1, tf2
        integer(I8), allocatable, dimension(:, :)   :: nint, segrf, &
            segc, vertexID
        real(R8)                                :: hfval, lfval, &
            dhf1, dhf2, dlf1, dlf2, nxc, nyc, nxfv, nyfv
        real(R8), allocatable, dimension(:)     :: tcvfval, tcfv1val, &
            tcfv2val, tx, ty, xl, yl, tfval, sr1, sr2, txint, tyint, &
            nxf, nyf, nnf
        real(R8), allocatable, dimension(:, :)  :: segrrf, segrc, &
            xint, yint
        logical                                 :: isflremoved, &
            isintersectremoved
        logical, allocatable, dimension(:)      :: ishfface, islfface, &
            ishfvert, iscontourfound, keepind
        type(ContourUDT), allocatable           :: tempc(:)
        type(contourUDT)                        :: c1, c2
        type(RealDynamicArrayUDT), allocatable, dimension(:, :)     :: &
            xintda, yintda, segrrfda, segrcda
        type(IntegerDynamicArrayUDT), allocatable, dimension(:, :)     :: &
            segcda, segrfda
        type(PolygonUDT), allocatable           :: polc(:)
        type(PolygonSetUDT)                     :: tempps

        ! Diagnostics
        real(R8)                                :: tstart, tend 

        ! Loop
        integer(I8)                             :: i, j, k, cc

        ! Initialize
        !===========
        ! Associate
        associate(&
            vert        => topomesh%vert,   &
            vertdata    => ggtmdata%vert,   &
            face        => topomesh%face,   &
            facedata    => ggtmdata%face,   &
            cell        => topomesh%cell,   &
            celldata    => ggtmdata%cell,   &
            tube        => topomesh%tube,   &
            tubedata    => ggtmdata%tube    &
            )

        ! Initialize
        nv = topomesh%vert%ntot

        ! Determine cell boundaries
        !==========================
        ! Loop over all tubes
        do i = 1, tube%ntot 
            ! Get the tube cells & faces
            tubec = tube%GetCell(i)
            tubef = tube%GetFace(i)

            if (tube%isclosed(i)) then 
                print *, 'closed tube'
            end if 

            ! Loop over all tube cells
            do j = 1, size(tubec)
                ! Unpack
                tc = tubec(j)

                ! Set cell region
                !----------------
                celldata(tc)%region = tc ! just equal to cell number for now

                ! Extract boundaries
                !-------------------
                ! Get the cell faces & vertices and precompute some values
                tcf = cell%GetFace(tc)
                tcv = cell%GetVert(tc)
                tcfv1 = face%vert(tcf, 1)
                tcfv2 = face%vert(tcf, 2)
                tcvfval = vert%fval(tcv)

                ! Set start radial face
                srf = tubef(j) ! should work 
                celldata(tc)%srf = srf 
                indsrf = findloc(tcf, srf, 1, back=.false.) ! take the first one - hedge for closed cell
                if (indsrf == 0) then 
                    ! This shouldn't be happening
                    call gdErrorHandler('AddTopologicalMeshCellGriddingData: ' // & 
                        'expected cell to have radial faces as defined in ' // & 
                        'flux tube, but face is not present in cell faces. ' // &
                        'Check topological mesh for inconsistencies.')
                end if 

                ! Set end radial face
                erf = tubef(j+1) ! should work 
                celldata(tc)%erf = erf
                inderf = findloc(tcf, erf, 1, back=.true.) ! take the last one - hedge for closed cell
                if (inderf == 0) then 
                    ! This shouldn't be happening
                    call gdErrorHandler('AddTopologicalMeshCellGriddingData: ' // & 
                        'expected cell to have radial faces as defined in ' // & 
                        'flux tube, but face is not present in cell faces. ' // &
                        'Check topological mesh for inconsistencies.')
                end if 

                ! Determine high and low field value
                hfval = maxval(vert%fval(face%vert(srf, :)))
                lfval = minval(vert%fval(face%vert(srf, :)))

                ! Check if, when querying the radial faces, we should 
                ! flip to get them from high to low field
                celldata(tc)%flipsrf = .false.
                celldata(tc)%fliperf = .false.
                if (vert%fval(face%vert(srf, 1)) < vert%fval(face%vert(srf, 2))) then 
                    celldata(tc)%flipsrf = .true.
                end if 
                if (vert%fval(face%vert(erf, 1)) < vert%fval(face%vert(erf, 2))) then 
                    celldata(tc)%fliperf = .true.
                end if 

                ! Determine high and low field poloidal faces (sorted, but not properly oriented)
                tcfv1val = vert%fval(tcfv1)
                tcfv2val = vert%fval(tcfv2)
                allocate(ishfface(size(tcf)), islfface(size(tcf)))
                ishfface = .false. 
                islfface = .false. 
                do k = 1, size(tcf)
                    ! Check if the face is of poloidal/sep type
                    if (any(face%type(tcf(k)) == [TMfacepolID, TMfacesepID])) then 
                        dhf1 = abs(tcfv1val(k) - hfval)
                        dhf2 = abs(tcfv2val(k) - hfval)
                        dlf1 = abs(tcfv1val(k) - lfval)
                        dlf2 = abs(tcfv2val(k) - lfval)
                        if ((dhf1 < dlf1) .and. (dhf2 < dlf2)) then 
                            ! High field face
                            ishfface(k) = .true. 
                        elseif ((dhf1 > dlf1) .and. (dhf1 > dlf1)) then 
                            ! Low field face
                            islfface(k) = .true. 
                        else 
                            ! Undetermined - throw error (should actually
                            ! not happen)
                            call gdErrorHandler('AddTopologicalMeshCellGriddingData: ' // & 
                                'could not determine based on field value ' // & 
                                'if poloidal face is high or low field.')
                        end if 
                    end if 
                end do 

                ! Extract in a sorted way
                minind = min(indsrf, inderf)
                maxind = max(indsrf, inderf)
                if (minind == maxind) then 
                    ! Sanity check failed, this shouldn't happen even 
                    ! if both start and end radial face are the same,
                    ! since they should then appear twice in the cell. 
                    call gdErrorHandler('AddTopologicalMeshCellGriddingData: ' // & 
                        'start and end radial face have the ' // & 
                        'same index, unexpected')
                end if 
                tf1 = [tcf(maxind+1:), tcf(:minind-1)]
                tf2 = tcf(minind+1:maxind-1)

                ! Sanity checks
                if (count(.not. ishfface .and. .not. islfface) /= 2) then 
                    call gdErrorHandler('AddTopologicalMeshCellGriddingData: ' // & 
                        'cell has not exactly two radial faces, this is not yet supported')
                end if 
                if ((any(ishfface(minind+1:maxind-1)) .and. &
                    any(islfface(minind+1:maxind-1))) .or. &
                    (any([ishfface(maxind+1:), ishfface(:minind-1)]) .and. &
                    any([islfface(maxind+1:), islfface(:minind-1)]))) then 
                        call gdErrorHandler('AddTopologicalMeshCellGriddingData: ' // & 
                        'not all extracted faces are strictly high or ' // & 
                        'low field faces, unexpected')
                end if 


                ! Determine high field and low field vertices (unsorted)
                ishfvert = abs(tcvfval - hfval) < abs(tcvfval - lfval)

                ! Extract
                allocate(hffaces(count(ishfface)), lffaces(count(islfface)), &
                    hfvert(count(ishfvert)), lfvert(count(.not.ishfvert)))
                hffaces = pack(tcf, ishfface)
                lffaces = pack(tcf, islfface)
                hfvert = pack(tcv, ishfvert)
                lfvert = pack(tcv, .not. ishfvert)
                
                ! Overwrite to sort 
                if (size(tf1) > 0) then 
                    if (any([ishfface(maxind+1:), ishfface(:minind-1)])) then 
                        hffaces = tf1 
                    else
                        lffaces = tf1 
                    end if 
                end if
                if (size(tf2) > 0) then 
                    if (any(ishfface(minind+1:maxind-1))) then 
                        hffaces = tf2
                    else
                        lffaces = tf2 
                    end if 
                end if

                ! Add
                celldata(tc)%hffaces = hffaces
                celldata(tc)%lffaces = lffaces
                celldata(tc)%hfvert = hfvert
                celldata(tc)%lfvert = lfvert


                ! Housekeeping
                deallocate(hffaces, lffaces, ishfface, islfface, hfvert, &
                    lfvert)
            end do 
        end do 

        ! Trace contours 
        !===============
        ! Loop over all tubes
        do i = 1, tube%ntot 

            ! Initialize
            !-----------
            tubef = tube%GetFace(i)
            tubec = tube%GetCell(i)

            if (tube%isclosed(i)) then 
                print *, 'closed tube'
            end if 

            ! Trace contours
            !---------------
            ! Get initial vertex distribution (skip start and end points)
            tf = tubedata(i)%distributionface
            tx = facedata(tf)%xv
            ty = facedata(tf)%yv

            !  Make sure we trace from high to low field value
            if (allocated(tfval)) then 
                deallocate(tfval)
            end if 
            allocate(tfval(size(tx)))
            call magneticField%interp%Evaluate(tx, ty, 0, 0, tfval)
            if (tfval(1) < tfval(size(tfval))) then 
                tx = tx(size(tx):1:-1)
                ty = ty(size(ty):1:-1)
            end if 

            ! Trace
            tempc = fieldtracer%TraceContours(tx(2:size(tx)-1), ty(2:size(ty)-1))

            ! Clean
            call CleanContours(tempc)

            ! Check if contours make sense and reformat if necessary
            allIDs = tempc%ID
            allocate(iscontourfound(size(tx)-2)) ! normally, IDs go from 1 to number of points
            allocate(keepind(size(tempc)))
            iscontourfound = .false. 
            keepind = .true. 
            if (tube%isclosed(i)) then 
                ! Precompute face normals for each flux surface for 
                ! determining orientation
                nxf = -(facedata(tf)%yv(2:) - facedata(tf)%yv(1:size(facedata(tf)%yv)-1))
                nyf = (facedata(tf)%xv(2:) - facedata(tf)%xv(1:size(facedata(tf)%xv)-1))
                nnf = sqrt(nxf**2 + nyf**2)
                nxf = nxf/nnf
                nyf = nyf/nnf

                if (tfval(1) < tfval(size(tfval))) then 
                    nxf = -nxf(size(nxf):1:-1) 
                    nyf = -nyf(size(nyf):1:-1) 
                end if

                ! Closed contour: only one contour value expected
                do j = 1, size(tempc)
                    ! Check if contour is closed
                    if (.not. tempc(j)%isclosed) then 
                        call gdErrorHandler('AddTopologicalMeshCellGriddingData :' // & 
                            'expected closed contour, but found open contour')
                    end if 
                    if (count(tempc(j)%ID == allIDs) > 1) then 
                        call gdErrorHandler('AddTopologicalMeshCellGriddingData :' // & 
                            'expected single closed contour, but found multiple')
                    end if 

                    ! Mark as found
                    iscontourfound(tempc(j)%ID) = .true. 

                    ! Ensure proper orientation
                    nxc = tempc(j)%x(2) - tempc(j)%x(1)
                    nyc = tempc(j)%y(2) - tempc(j)%y(1)
                    nxfv = 0.5*(nxf(j) + nxf(j+1))
                    nyfv = 0.5*(nyf(j) + nyf(j+1))
                    if ((nxc*nxfv + nyc*nyfv) < 0) then 
                        ! Flip the contour
                        tempc(j)%x = tempc(j)%x(size(tempc(j)%x):1:-1)
                        tempc(j)%y = tempc(j)%y(size(tempc(j)%y):1:-1)
                        temp = tempc(j)%startsaddle
                        tempc(j)%startsaddle = tempc(j)%endsaddle
                        tempc(j)%endsaddle = temp
                    end if 

                end do

            else
                ! Either one closed contour or two open contours
                do j = 1, size(tempc)
                    ! Unpack for ease
                    c1 = tempc(j)

                    if (.not. iscontourfound(c1%ID)) then 

                        if (c1%isclosed) then 
                            ! Check if only one contour exists
                            if (count(c1%ID == allIDs) > 1) then 
                                call gdErrorHandler('AddTopologicalMeshCellGriddingData :' // & 
                                    'expected single closed contour, but found multiple')
                            end if 
                        else
                            ! Open contour, should only have two parts
                            if (count(c1%ID == allIDS) < 2) then 
                                call gdErrorHandler('AddTopologicalMeshCellGriddingData :' // & 
                                    'expected two open contours, but found less')
                            elseif (count(c1%ID == allIDS) > 2) then 
                                call gdErrorHandler('AddTopologicalMeshCellGriddingData :' // & 
                                    'expected two open contours, but found more')
                            end if 

                            ! If we got here, only two contours left. 
                            cind = findloc((allIDs == c1%ID) .and. &
                                (j /= [(k, k = 1, size(allIDs))]), .true., 1)
                            if (cind == 0) then 
                                ! This should be impossible
                                call gdErrorHandler('AddTopologicalMeshCellGriddingData :' // & 
                                    'could not find second contour index, this is a bug')
                            end if
                            c2 = tempc(cind)

                            ! Check if both contours start in the same point
                            if (.not.((c1%x(1) == c2%x(1)) .and. (c1%y(1) == c2%y(1)))) then 
                                call gdErrorHandler('AddTopologicalMeshCellGriddingData :' // & 
                                    'two open contours found, but starting ' // &
                                    'points differ. Cannot merge and proceed')
                            end if 

                            ! Merge the contour
                            c1%x = [c1%x(size(c1%x):2:-1), c2%x]
                            c1%y = [c1%y(size(c1%y):2:-1), c2%y]
                            c1%startsaddle = c1%endsaddle
                            c1%endsaddle = c2%endsaddle
                            c1%isclosed = .false. ! normally, checked later
                            tempc(j) = c1
                            
                            ! Check if the contour is open, should be the case
                            if ((c1%x(1) == c1%x(size(c1%x))) .and. (c1%y(1) == c1%y(size(c1%y)))) then 
                                call gdErrorHandler('AddTopologicalMeshCellGriddingData :' // & 
                                    'two open contours found that form a closed contour, ' // &
                                    'may be a bug in the contouring algorithm')
                            end if

                            ! Mark as found
                            iscontourfound(c1%ID) = .true. 

                            ! Mark for removal
                            keepind(cind) = .false. 
                        end if 
                    end if 
                end do 
            end if 

            ! Check if all contours were found
            if (.not. all(iscontourfound)) then 
                call gdErrorHandler('AddTopologicalMeshCellGriddingData :' // & 
                    'some contours of starting points were not traced')
            end if

            ! Remove contours
            tempc = pack(tempc, keepind)

            ! Housekeeping
            deallocate(keepind, iscontourfound)

            ! Compute intersections 
            !----------------------
            nc = size(tempc)
            ntf = size(tubef)
            allocate(xintda(nc, ntf), yintda(nc, ntf), nint(nc, ntf), &
                segrfda(nc, ntf), segcda(nc, ntf), segrrfda(nc, ntf), &
                segrcda(nc, ntf), polc(nc))
            nint = 0

            ! Convert to polygon
            call wall_time(tstart)
            !$omp parallel do private(j)
            do j = 1, nc
                call polc(j)%Construct(tempc(j)%x, tempc(j)%y)
            end do 
            !$omp end parallel do
            call wall_time(tend)
            !print *, 'time spent in polygon setup:', tend-tstart

            ! Loop over all contours
            call wall_time(tstart)
            !$omp parallel do private(i, j, txint, tyint, s1, s2, sr1, sr2)
            do j = 1, size(tempc)
                ! Compute intersections with each radial face's polygon
                do k = 1, size(tubef)
                    ! Compute intersections
                    call PolygonIntersections(polc(j), face%pol(tubef(k)), &
                        txint, tyint, s1, s2, sr1, sr2)

                    ! Add
                    xintda(j, k)    = ConstructRealDynamicArray(txint)
                    yintda(j, k)    = ConstructRealDynamicArray(tyint)
                    segrfda(j, k)   = ConstructIntegerDynamicArray(s2)
                    segcda(j, k)    = ConstructIntegerDynamicArray(s1)
                    segrrfda(j, k)  = ConstructRealDynamicArray(sr2)
                    segrcda(j, k)   = ConstructRealDynamicArray(sr1)
                    nint(j, k) = size(txint)
                end do 
            end do 
            !$omp end parallel do
            call wall_time(tend)
            !print *, 'time spent in intersections:', tend-tstart
            ! Print
            call tempps%Construct(polc)
            call tempps%WriteData('lines')
            call tempps%Construct(face%pol(tubef))
            call tempps%WriteData('faces')

            !!! we need a different approach for closed ft!
            ! Check intersections & sort
            !---------------------------
            allocate(keepind(nc))
            keepind = .true. 
            isflremoved = .false. 
            isintersectremoved = .false. 
            do j = 1, nc
                ! Check if no intersections with radial lines
                if (any(nint(j, :) == 0)) then 
                    ! Mark for removal
                    keepind(j) = .false. 
                    isflremoved = .true.
                    
                    ! Go to the next line
                    cycle 
                end if 

                ! Check if multiple intersections with radial lines
                if (tube%isclosed(i)) then 
                    ! Two intersections are expected in first and 
                    ! last face, since they are the same face
                    if (any(nint(j, 2:ntf-1) > 1)) then 
                        ! Simply set warning message
                        isintersectremoved = .true. 
                    end if 
                    if ((nint(j, 1) > 2) .or. (nint(j, ntf) > 2)) then 
                        isintersectremoved = .true. 
                    end if 
                else
                    if (any(nint(j, :) > 1)) then 
                        ! Simply set warning message
                        isintersectremoved = .true. 
                    end if 
                end if 
            end do 

            ! Remove field lines
            tempc = pack(tempc, keepind)
            nc = count(keepind)

            ! Issue warnings
            if (isflremoved) then 
                print *, 'AddTopologicalMeshCellGriddingData: ' // & 
                    'field lines were removed since they do not intersect ' // & 
                    'with one or more radial lines for tube: ', i 
            end if 
            if (isintersectremoved) then 
                print *, 'AddTopologicalMeshCellGriddingData: ' // & 
                    'multiple intersections found with some radial lines ' // & 
                    'for tube: ', i, ', taking first intersection'
            end if  

            ! Unpack intersections
            allocate(xint(nc, ntf), yint(nc, ntf), segc(nc, ntf), &
                segrf(nc, ntf), segrc(nc, ntf), segrrf(nc, ntf), &
                vertexID(nc, ntf))
            cc = 0 
            do j = 1, size(nint, 1) 
                ! Skip
                if (keepind(j)) then 
                    ! Update counter
                    cc = cc + 1 

                    ! Loop
                    do k = 1, ntf-1 
                        ! Get first intersection
                        xint(cc, k) = xintda(j, k)%Get(1)
                        yint(cc, k) = yintda(j, k)%Get(1)
                        segc(cc, k) = segcda(j, k)%Get(1)
                        segrf(cc, k) = segrfda(j, k)%Get(1)
                        segrc(cc, k) = segrcda(j, k)%Get(1)
                        segrrf(cc, k) = segrrfda(j, k)%Get(1)

                        ! Set vertex ID
                        nv = nv + 1
                        vertexID(cc, k) = nv
                    end do 
                    
                    ! Hedge for closed flux surfaces
                    if (tube%isclosed(i)) then 
                        ! Intersection with first and last radial line
                        ! should be exactly the same! Only intersection
                        ! coordinate should differ
                        xint(cc, ntf) = xint(cc, 1)
                        yint(cc, ntf) = yint(cc, 1)
                        segc(cc, ntf) = segcda(j, ntf)%Get(segcda(j, ntf)%Size())
                        segrf(cc, ntf) = segrfda(j, ntf)%Get(segrfda(j, ntf)%Size())
                        segrc(cc, ntf) = segrcda(j, ntf)%Get(segrcda(j, ntf)%Size())
                        segrrf(cc, ntf) = segrrfda(j, ntf)%Get(segrrfda(j, ntf)%Size())

                        ! Set vertex ID - should be the same as before
                        nv = nv + 1
                        vertexID(cc, ntf) = vertexID(cc, 1)
                    else
                        ! Get first intersection
                        xint(cc, ntf) = xintda(j, ntf)%Get(1)
                        yint(cc, ntf) = yintda(j, ntf)%Get(1)
                        segc(cc, ntf) = segcda(j, ntf)%Get(1)
                        segrf(cc, ntf) = segrfda(j, ntf)%Get(1)
                        segrc(cc, ntf) = segrcda(j, ntf)%Get(1)
                        segrrf(cc, ntf) = segrrfda(j, ntf)%Get(1)

                        ! Set vertex ID
                        nv = nv + 1
                        vertexID(cc, ntf) = nv
                    end if 
                end if 
            end do 
            deallocate(keepind)

            ! Do sanity checks for closed tube
            if (tube%isclosed(i)) then 
                ! First and last intersection should be at 0 
                ! and ne coordinate by construction
                do j = 1, nc
                    if (segrc(j, 1) /= 0_R8) then 
                        call gdErrorHandler('AddTopologicalMeshCellGriddingData: ' // & 
                            'first intersection of closed contour is not at start of ' // & 
                            'contour for closed flux tube, unexpected')
                    end if
                    if (segrc(j, size(segrc, 2)) /= polc(j)%nv) then 
                        call gdErrorHandler('AddTopologicalMeshCellGriddingData: ' // & 
                            'last intersection of closed contour is not at end of ' // & 
                            'contour for closed flux tube, unexpected')
                    end if 
                end do 
            end if 

            ! Add lines for each cell
            do k = 1, ntf-1
                ! Allocate the amount of lines for this cell
                if (allocated(celldata(tubec(k))%lines)) then
                    deallocate(celldata(tubec(k))%lines)
                end if 
                allocate(celldata(tubec(k))%lines(nc))

                ! Loop over lines
                do j = 1, nc
                    ! Check line order
                    if (segrc(j, k) < segrc(j, k+1)) then 
                        incr = 1
                    else
                        incr = -1
                    end if 

                    ! Note: we ensure no duplicate points by skipping the
                    ! first vertex of the face segment
                    polv = polc(j)%vert(segc(j, k)+2:segc(j, k+1)-1:incr)
                    xl = [xint(j, k), polc(j)%x(polv), xint(j, k+1)]
                    yl = [yint(j, k), polc(j)%y(polv), yint(j, k+1)]

                    ! Add 
                    celldata(tubec(k))%lines(j)%xl = xl
                    celldata(tubec(k))%lines(j)%yl = yl

                end do 

                ! Add vertex IDs of start and end radial face (without
                ! start and end point of cell boundary, only additional
                ! lines)
                celldata(tubec(k))%srfvert = vertexID(:, k)
                celldata(tubec(k))%erfvert = vertexID(:, k+1)

                ! Add labels of radial faces
                celldata(tubec(k))%srflabel = face%label(tubef(k))
                celldata(tubec(k))%erflabel = face%label(tubef(k+1))

            end do 

            ! Housekeeping
            deallocate(xint, yint, segc, segrf, segrc, segrrf, nint, &
                vertexID, polc, xintda, yintda, segcda, segrfda, &
                segrcda, segrrfda,)
        end do 
        
        ! Housekeeping
        !=============
        end associate


    end subroutine 

    !------------------------------------------------------------------!
    !                       GRID DATA ROUTINES                         !
    !------------------------------------------------------------------!

    ! Initializers
    subroutine InitializeGGVert(vert, storagetype)

        ! Declare variables
        !==================
        ! Arguments
        class(GGVertUDT)                    :: vert 
        character(*), intent(in)            :: storagetype 

        ! Initialize
        !===========
        vert%ntot = 0
        select case (storagetype)

        case ('standard')

            vert%x = ConstructRealDynamicArray()
            vert%y = ConstructRealDynamicArray()
            vert%fieldlineID = ConstructIntegerDynamicArray()

        case default 

            call gdErrorHandler('InitializeGGVert: unknown storage type')

        end select

    end subroutine

    subroutine InitializeGGFace(face, storagetype)

        ! Declare variables
        !==================
        ! Arguments
        class(GGFaceUDT)                    :: face 
        character(*), intent(in)            :: storagetype 

        ! Initialize
        !===========
        face%ntot = 0
        select case (storagetype)

        case ('standard')

            face%v1 = ConstructIntegerDynamicArray()
            face%v2 = ConstructIntegerDynamicArray()
            face%label = ConstructIntegerDynamicArray()

        case default 

            call gdErrorHandler('InitializeGGFace: unknown storage type')

        end select

    end subroutine

    subroutine InitializeGGCell(cell, storagetype)

        ! Declare variables
        !==================
        ! Arguments
        class(GGCellUDT)                    :: cell 
        character(*), intent(in)            :: storagetype 

        ! Initialize
        !===========
        cell%ntot = 0
        select case (storagetype)

        case ('standard')

            cell%vert = ConstructIntegerDynamicArray()
            cell%vp1 = ConstructIntegerDynamicArray()
            cell%vp2 = ConstructIntegerDynamicArray()
            cell%region = ConstructIntegerDynamicArray()

        case default 

            call gdErrorHandler('InitializeGGCell: unknown storage type')

        end select

    end subroutine

    subroutine InitializeGGGrid(grid, storagetype)

        ! Declare variables
        !==================
        ! Arguments
        class(GGGridUDT)                    :: grid 
        character(*), intent(in)            :: storagetype 

        ! Initialize
        !===========
        call grid%vert%Initialize(storagetype)
        call grid%face%Initialize(storagetype)
        call grid%cell%Initialize(storagetype)

    end subroutine

    ! Vertex addition
    subroutine AddGGVert(grid, xv, yv, flID)

        ! Description
        !============
        ! Add vertices to the grid (without updating interconnection)

        ! Declare variables
        !==================
        ! Arguments
        class(GGGridUDT)            :: grid 
        real(R8), intent(in)        :: xv(:), yv(:)
        integer(I8), intent(in)     :: flID(:)

        ! Checks
        !=======
        ! Size checks
        if ((size(xv) /= size(yv)) .or. (size(xv) /= size(flID))) then 
            call gdErrorHandler('AddGGVert: incompatible input sizes')
        end if 

        ! Add
        !====
        call grid%vert%x%Append(xv)
        call grid%vert%y%Append(yv)
        call grid%vert%fieldlineID%Append(flID)

    end subroutine

    ! Face addition
    subroutine AddGGFace(grid, facevert, facelabel, faceregion)

        ! Description
        !============
        ! Add faces to the grid (without updating interconnection)

        ! Declare variables
        !==================
        ! Arguments
        class(GGGridUDT)            :: grid 
        integer(I8), intent(in)     :: facevert(:, :), facelabel(:), &
            faceregion(:)

        ! Checks
        !=======
        ! Size checks
        if ((size(facevert, 1) /= size(facelabel)) .or. (size(facelabel) /= size(faceregion))) then 
            call gdErrorHandler('AddGGFace: incompatible input sizes')
        end if 
        if (size(facevert, 2) /= 2) then 
            call gdErrorHandler('AddGGFace: wrong second dimension of face vertices')
        end if 

        ! Add
        !====
        call grid%face%v1%Append(facevert(:, 1))
        call grid%face%v2%Append(facevert(:, 2))
        call grid%face%label%Append(facelabel)
        call grid%face%region%Append(faceregion)

    end subroutine

    ! Cell addition
    subroutine AddGGCell(grid, cellvert, cellvertP, cellregion)

        ! Description
        !============
        ! Add cells to the grid (without updating interconnection)

        ! Declare variables
        !==================
        ! Arguments
        class(GGGridUDT)            :: grid 
        integer(I8), intent(in)     :: cellvertP(:, :), cellvert(:), &
            cellregion(:)

        ! Checks
        !=======
        ! Size checks
        if ((size(cellvertP, 1) /= size(cellvert)) .or. (size(cellvert) /= size(cellregion))) then 
            call gdErrorHandler('AddGGcell: incompatible input sizes')
        end if 
        if (size(cellvertP, 2) /= 2) then 
            call gdErrorHandler('AddGGcell: wrong second dimension of cell vertex pointer')
        end if 

        ! Add
        !====
        call grid%cell%vp1%Append(cellvertP(:, 1))
        call grid%cell%vp2%Append(cellvertP(:, 2))
        call grid%cell%vert%Append(cellvert)
        call grid%cell%region%Append(cellregion)

    end subroutine

    !------------------------------------------------------------------!
    !                           AUXILIARY                              !
    !------------------------------------------------------------------!

    ! Boundary line extractor for cells
    subroutine ExtractTMCellAlignedBoundary(tmcell, loc, ggtmdata, &
        topomesh, line)

        ! Description
        !============
        ! This routine extracts a cell boundary that is oriented in the
        ! correct direction. 'loc' indicates the location ('high', or 
        ! 'low') of the boundary that should be extracted. The 
        ! extraction format is a fieldlinedata type along with 
        ! additional information needed for distributing vertices and
        ! updating information. Note that it is possible that the 
        ! boundary is a single point. In that case, the line only 
        ! holds one coordinate and the additional quantities are 
        ! allocated as empty arrays. 

        ! Note 1: the facedata is assumed to hold vertex distributions 
        ! (xv, yv) for all faces, including vertex IDs. 

        ! Note 2: it is assumed that the face data is already properly
        ! sorted, but that the faces may still need to be flipped. 

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMCellDataUDT)                  :: tmcell 
        character(*), intent(in)                :: loc 
        class(GGTMFieldlineDataUDT)             :: line 
        class(TopomeshUDT), intent(in)          :: topomesh 
        class(GGTMDataUDT)                      :: ggtmdata 

        ! Auxiliary
        real(R8)                                :: nxsrf, nysrf, nxl, &
            nyl
        real(R8), allocatable, dimension(:)     :: tx, ty 
        integer(I8)                             :: startv, endv, &
            thisf, thisfind, nextv
        integer(I8), allocatable, dimension(:)  :: bndf, bndv 
        logical                                 :: doflip, facefound
        logical, allocatable, dimension(:)      :: isnotfound 

        ! Loop
        integer(I8)                             :: i, j 

        ! Initialize
        !===========
        ! Associate
        associate(&
            vert        => topomesh%vert,   &
            face        => topomesh%face,   &
            facedata    => ggtmdata%face    &
            )

        ! Check boundary to extract
        select case (loc)

        case ('high')

            ! Get faces & vert
            bndf = tmcell%hffaces
            bndv = tmcell%hfvert

            ! Get the high field side vertex of start and end radial face
            associate( & 
                srfv => topomesh%face%vert(tmcell%srf, :), &
                erfv => topomesh%face%vert(tmcell%erf, :), &
                srfvx   => facedata(tmcell%srf)%xv,     &
                srfvy   => facedata(tmcell%srf)%yv)

            ! Check field values
            if (vert%fval(srfv(1)) > vert%fval(srfv(2))) then 
                startv = srfv(1)
                nxsrf = -(srfvy(2) - srfvy(1))
                nysrf =  (srfvx(2) - srfvx(1))
            else 
                startv = srfv(2)
                nxsrf =  (srfvy(size(srfvy)) - srfvy(size(srfvy)-1))
                nysrf = -(srfvx(size(srfvx)) - srfvx(size(srfvx)-1))
            end if
            if (vert%fval(erfv(1)) > vert%fval(erfv(2))) then 
                endv = erfv(1)
            else 
                endv = erfv(2)
            end if

            ! Houskeeping
            end associate

        case ('low')

            ! Get faces & vert
            bndf = tmcell%lffaces
            bndv = tmcell%lfvert

            ! Get the low field side vertex of start and end radial face
            associate( & 
                srfv => topomesh%face%vert(tmcell%srf, :), &
                erfv => topomesh%face%vert(tmcell%erf, :), &
                srfvx   => facedata(tmcell%srf)%xv,     &
                srfvy   => facedata(tmcell%srf)%yv)

            ! Check field values
            if (vert%fval(srfv(1)) < vert%fval(srfv(2))) then 
                startv = srfv(1)
                nxsrf =  (srfvy(2) - srfvy(1))
                nysrf = -(srfvx(2) - srfvx(1))
            else 
                startv = srfv(2)
                nxsrf = -(srfvy(size(srfvy)) - srfvy(size(srfvy)-1))
                nysrf =  (srfvx(size(srfvx)) - srfvx(size(srfvx)-1))
            end if
            if (vert%fval(erfv(1)) < vert%fval(erfv(2))) then 
                endv = erfv(1)
            else 
                endv = erfv(2)
            end if

            ! Houskeeping
            end associate

        case default 

            call gdErrorHandler('ExtractTMCellAlignedBoundary: unknown location')

        end select

        ! Check for the trivial case
        if (size(bndf) == 0) then 
            ! Just a sanity check
            if (size(bndv) == 1) then 
                ! Set the line data
                line%xl = vert%x(bndv)
                line%yl = vert%y(bndv)
                line%xv = line%xl
                line%yv = line%yl
                line%vert = bndv
                if (allocated(line%facelabels)) then 
                    deallocate(line%facelabels)
                end if 
                allocate(line%facelabels(0))

            else 
                call gdErrorHandler('ExtractTMCellAlignedBoundary: expected ' // & 
                    'a single boundary vertex but found multiple/none')
            end if
            
            ! Return 
            return 
        end if 

        ! Sort faces
        !===========
        ! Initialize
        allocate(isnotfound(size(bndf)))
        isnotfound = .true. 

        ! Get the starting aligned face (should be first or last face) - 
        ! if it's the last face, flip the boundary order for ease
        thisfind = 1
        if (findloc(face%vert(bndf(1), :), startv, 1) /= 0) then 
            ! No flipping of bndf needed
        elseif (findloc(face%vert(bndf(size(bndf)), :), startv, 1) /= 0) then 
            ! Need to flip bndf
            bndf = bndf(size(bndf):1:-1)
        else
            ! Shouldn't happen
            call gdErrorHandler('ExtractTMCellAlignedBoundary: both start and ' // & 
                'end face do not have the starting vertex, unexpected' )
        end if 
        
        ! Check if we need to flip the orientation
        doflip = .false. 
        if (face%vert(bndf(1), 2) == startv) then 
            doflip = .true. 
        end if  

        ! Set as found
        isnotfound(thisfind) = .false.

        ! Add face data
        thisf   = bndf(thisfind)
        line%xl = face%x(thisf)%Get()
        line%yl = face%y(thisf)%Get()
        line%xv = facedata(thisf)%xv
        line%yv = facedata(thisf)%yv
        line%vert = facedata(thisf)%vert
        line%facelabels = spread(face%label(thisf), 1, size(line%vert)-1)

        if (doflip) then 
            line%xl = line%xl(size(line%xl):1:-1)
            line%yl = line%yl(size(line%yl):1:-1)
            line%xv = line%xv(size(line%xv):1:-1)
            line%yv = line%yv(size(line%yv):1:-1)
            line%vert = line%vert(size(line%vert):1:-1)
            line%facelabels = line%facelabels(size(line%facelabels):1:-1)
        end if 

        ! Get the next vertex
        if (doflip) then 
            nextv = face%vert(thisf, 1)
        else
            nextv = face%vert(thisf, 2)
        end if
        
        ! Loop over remaining faces (should be properly ordened, but we do
        ! some sanity checks)
        do i = 2, size(bndf)
            ! Check
            doflip = .false. 
            if (face%vert(bndf(i), 1) == nextv) then 
                ! all good
                nextv = face%vert(bndf(i), 2)
            elseif (face%vert(bndf(i), 2) == nextv) then 
                ! need to flip
                doflip = .true. 
                nextv = face%vert(bndf(i), 1)
            else 
                ! All bad
                call gdErrorHandler('ExtractTMCellAlignedBoundary: '// & 
                    'could not find next face, check input')
            end if 

            ! Add face data 
            if (doflip) then 
                tx = face%x(bndf(i))%Get()
                ty = face%y(bndf(i))%Get()
                line%xl = [line%xl(1:size(line%xl)-1), tx(size(tx):1:-1)] ! avoid double coordinates
                line%yl = [line%yl(1:size(line%yl)-1), ty(size(ty):1:-1)]
                line%xv = [line%xv(1:size(line%xv)-1), &
                    facedata(bndf(i))%xv(size(facedata(bndf(i))%xv):1:-1)]
                line%yv = [line%yv(1:size(line%yv)-1), &
                    facedata(bndf(i))%yv(size(facedata(bndf(i))%yv):1:-1)]
                line%vert = [line%vert(1:size(line%vert)-1), &
                    facedata(bndf(i))%vert(size(facedata(bndf(i))%vert):1:-1)]
                line%facelabels = [line%facelabels(1:size(line%facelabels)-1), &
                    spread(face%label(bndf(i)), 1, size(line%vert)-1)]
            
            else
                line%xl = [line%xl(1:size(line%xl)-1), face%x(bndf(i))%Get()] ! avoid double coordinates
                line%yl = [line%yl(1:size(line%yl)-1), face%y(bndf(i))%Get()]
                line%xv = [line%xv(1:size(line%xv)-1), facedata(bndf(i))%xv]
                line%yv = [line%yv(1:size(line%yv)-1), facedata(bndf(i))%yv]
                line%vert = [line%vert(1:size(line%vert)-1), facedata(bndf(i))%vert]
                line%facelabels = [line%facelabels(1:size(line%facelabels)-1), &
                    spread(face%label(bndf(i)), 1, size(line%vert)-1)]
            end if 

        end do

        ! Sanity check: last vertex should be end vertex
        if (nextv /= endv) then 
            call gdErrorHandler('ExtractTMCellAlignedBoundary: '// & 
                'all faces found, but last vertex does not equal end ' // & 
                'vertex, unexpected')
        end if 
        
        ! Check if the line closes upon itself. If that is the case, 
        ! then we need to ensure a correct orientation (positive dot 
        ! product with radial line that is sorted from high field to
        ! low field)
        if (startv == endv) then 
            ! Get vector along line
            nxl = line%xl(2) - line%xl(1)
            nyl = line%yl(2) - line%yl(1)

            ! Compute dot product and check
            if ((nxl*nxsrf + nyl*nysrf) < 0) then 
                line%xl = line%xl(size(line%xl):1:-1)
                line%yl = line%yl(size(line%yl):1:-1)
                line%xv = line%xv(size(line%xv):1:-1)
                line%yv = line%yv(size(line%yv):1:-1)
                line%vert = line%vert(size(line%vert):1:-1)
                line%facelabels = line%facelabels(size(line%facelabels):1:-1)
            end if 
            
        end if 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! GGTM data writing 
    subroutine WriteGGTMData(ggtmdata, filename)

        ! Description
        !============
        ! Routine to write ggtmdata in our own format. Also
        ! handy for reading in once constructed. The following data is
        ! written in a column-wise fashion:

        ! 'faces'
        ! <face%ntot, number of printed faces> 
        ! 'face coordinates'
        ! 'face <nf>' <repeated for each face including header>
        ! <vertID, xv, yv> 
        ! 'cells'
        ! <cell%ntot, number of printed> 
        ! 'ID, srf, erf, cell line size>'
        ! <ID, nl>
        ! 'cell lines <nc>' (repeated nc times)
        ! <vertID, x, y> 
        ! 'tubes'
        ! <tubes%ntot> 
        ! 'tube distribution faces'
        ! <tube%distributionface>

        ! This information is best used in conjuction with the topomesh
        ! information to visualize the desired vertex distributions. Note
        ! that the vertices and their coordinates are only implicitly 
        ! present here through the cell lines and face vertex 
        ! distributions

        ! Declare variables
        !==================
        ! Modules 
        use mod_plotter 
        use mod_specialchars, only : filesepchar

        ! Arguments
        class(GGTMDataUDT)                      :: ggtmdata
        character(*), intent(in)                :: filename 

        ! Auxiliary
        integer                                 :: fu, nf, nc, nl
        logical, allocatable, dimension(:)      :: doface, docell, dolines 
        character(:), allocatable               :: dir

        ! Loop
        integer(I8)                             :: i, j, k

        ! Initialize
        !===========
        ! Unpack
        associate(&
            f       => ggtmdata%face,   &
            c       => ggtmdata%cell,   &
            t       => ggtmdata%tube    &
        )

        ! Construct writing directory
        dir = plotdir // filesepchar // filename // '.dat'

        ! Open file
        open (action='write', file=trim(dir), newunit=fu, &
             status='unknown')

        ! Write header
        write(fu, *) 'VERSION3.00.00'

        ! Check data
        !===========
        ! Only the fields that are not added when constructing vertices,
        ! faces or cells. The rest should be always available
        nf = 0
        allocate(doface(size(f)))
        doface = .false. 
        do i = 1, size(f)
            if (allocated(f(i)%xv) .and. allocated(f(i)%yv) .and. &
                allocated(f(i)%vert)) then 
                nf = nf + 1
                doface(i) = .true. 
            end if 
        end do 

        allocate(docell(size(c)), dolines(size(c)))
        docell = .false.
        dolines = .false.
        nc = 0
        do i = 1, size(c)
            if (allocated(c(i)%srfvert)) then 
                ! Assume other allocated as well
                nc = nc + 1
                docell(i) = .true.
                if (allocated(c(i)%lines)) then 
                    ! Lines should be determined as well
                    dolines(i) = .true.
                end if  
            end if 
        end do 

        ! Write face data
        !================
        ! Number of faces
        write (fu, *) 'faces'
        write (fu, *) size(f), nf

        ! Face data
        write (fu, *) 'ID, nc'
        do i = 1, size(f)
            if (doface(i)) then 
                ! Write ID and number of coordinates
                write (fu, *) i, size(f(i)%vert)
            end if 
        end do 
        
        ! Face coordinates
        write (fu, *) 'face coordinates'
        do i = 1, size(f)
            if (doface(i)) then 
                ! Write header
                write (fu, *) 'face ', i 

                ! Write coordinates
                do j = 1, size(f(i)%xv)
                    write(fu, *) f(i)%vert(j), f(i)%xv(j), f(i)%yv(j)
                end do 
            end if 
        end do 

        ! Write cell data
        !================
        ! Number of cells
        write (fu, *) 'cells'
        write (fu, *) size(c), nc

        ! Cell data
        write (fu, *) 'ID, srf, erf, cell line size'
        do i = 1, size(c) 
            if (docell(i)) then 
                if (dolines(i)) then 
                    nl = size(c(i)%lines) + 2 ! account for hfline, lfline
                else
                    nl = 0 !  don't print lines
                end if 
                write (fu, *) i, c(i)%srf, c(i)%erf, nl 
            end if 
        end do 

        ! Cell lines
        write (fu, *) 'cell lines'
        do i = 1, size(c) 
            if (dolines(i)) then 
                ! Write header
                write (fu, *) 'cell ', i 

                ! Write high field line
                write (fu, *) size(c(i)%hfline%vert)
                do j = 1, size(c(i)%hfline%vert)
                    write (fu, *) c(i)%hfline%vert(j), c(i)%hfline%xv(j), c(i)%hfline%yv(j)
                end do 

                ! Write other lines
                do k = 1, size(c(i)%lines)
                    write (fu, *) size(c(i)%lines(k)%vert)
                    do j = 1, size(c(i)%lines(k)%vert)
                        write (fu, *) c(i)%lines(k)%vert(j), c(i)%lines(k)%xv(j), c(i)%lines(k)%yv(j)
                    end do 
                end do 

                ! Write low field line
                write (fu, *) size(c(i)%lfline%vert)
                do j = 1, size(c(i)%lfline%vert)
                    write (fu, *) c(i)%lfline%vert(j), c(i)%lfline%xv(j), c(i)%lfline%yv(j)
                end do 
            end if 
        end do 

        ! Housekeeping
        !=============
        ! Deallocate again

        ! Others
        end associate
        close(fu)

    end subroutine



end module 