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
    use goatmod_types
    use mod_precision
    use mod_constants
    use mod_errorhandler
    use mod_dynamicarrays
    use mod_contour2D
    use mod_streamlinetracing2D
    use interpolant1D, only: Interpolate1D
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
    use omp_lib
    implicit none
    private 
    public :: GenerateUnstructuredAlignedGrid

    ! Module parameters
    real(R8), parameter, private        :: tprelfieldtol = 1e-10 ! relative field tolerance under which extrema are removed
    real(R8), parameter, private        :: disttol = 1e-12 ! distance tolerance
    integer(I8), private                :: verbosity = 1 ! verbosity level

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

    ! Field line data
    type :: GGTMFieldlineDataUDT

        real(R8), allocatable, dimension(:)         :: xv, yv, xl, yl, &
            dll, dllc, dlcv
        integer(I8), allocatable, dimension(:)      :: vert, facelabels
        integer(I8)                                 :: fsID, nv, nl

    contains 

        ! Initialization
        procedure :: Initialize     => InitializeGGTMFieldLineData

        ! Vertex coordinates addition
        procedure :: AddVertexCoordinates 

        ! Vertex ID addition
        procedure :: AddVertexIDs

        ! Data addition
        procedure :: UpdateLineData

        ! GGTM data updating
        procedure :: UpdateGGTMData

    end type
    
    ! Face data
    type :: GGTMFaceDataUDT

        ! Description
        !============
        ! Contains additional face data, including the grid vertices
        ! xv, yv - now stored as a GGTMFieldlineData type
        type(GGTMFieldlineDataUDT)              :: line 

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

        ! Data writing
        procedure :: WriteDAta          => WriteGGGridData
    end type 

    ! GGTM line refinement
    !=====================
    ! Abstract type
    type, abstract :: GGTMLineRefiner2DUDT

        ! Description
        !============
        ! This type has several routines to refine the vertex 
        ! distribution on a GGTM line type. Construction of the 
        ! type is deferred to a dedicated function. 

    contains 

        ! Refine line (single line based)
        procedure(RefineGGTMLineSingleINT), deferred :: RefineLineSingle
        generic :: Refine   => RefineLineSingle

    end type

    ! No refinement (dummy)
    type, extends(GGTMLineRefiner2DUDT)     :: GGTMLineRefinerNoRefUDT

        ! Description
        !============
        ! This refiner doesn't do any refinement and simply returns
        ! the original distribution. No additional data needs to be
        ! stored. 

    contains 

        ! Refine line
        procedure :: RefineLineSingle   => RefineLineSingleNoRef

    end type

    ! Length based refiner 
    type, extends(GGTMLineRefiner2DUDT)     :: GGTMLineRefinerLB2DUDT

        ! Description
        !============
        ! This refiner is based on a minimal and maximal length 
        ! distribution. 
        character(:), allocatable                       :: meth 
        class(DistributionFunctionUDT), allocatable     :: Lmin, Lmax 

    contains 

        ! Refine line
        procedure :: RefineLineSingle   => RefineLineSingleLB 

    end type

    !==================================================================!
    !                                                                  !
    !                          INTERFACES                              !
    !                                                                  !
    !==================================================================!

    ! Abstract 
    !=========
    abstract interface 

        ! Refine GGTM line
        subroutine RefineGGTMLineSingleINT(refiner, line, vertID, &
            keepvert)

            import :: GGTMFieldlineDataUDT, I8, GGTMLineRefiner2DUDT
            class(GGTMLineRefiner2DUDT)                 :: refiner 
            type(GGTMFieldlineDataUDT), intent(inout)   :: line
            integer(I8), intent(inout)                  :: vertID 
            logical, intent(in)                         :: keepvert(:) 

        end subroutine

    end interface

    ! Normal
    !=======
    ! Interface for line refiner constructors
    interface ConstructGGTMLineRefiner 
        module procedure ConstructGGTMLineRefinerLB 
    end interface
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
    subroutine GenerateUnstructuredAlignedGrid(simgrid, topomesh, magneticField, &
        vessel, fieldtracer, boundarytracer, streamlinetracer, options)

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
        type(GridUDT)               :: simgrid

        ! Auxiliary
        type(GGTMDataUDT)           :: ggtmdata 
        class(VertexDistributor2DUDT), allocatable      :: &
            poloidalvertexdistributor, radialvertexdistributor
        class(DistributionFunctionUDT), allocatable     :: & 
            magneticFieldDF 
        class(StreamlineTracerUDT), intent(inout)   :: streamlinetracer
        class(GGTMLineRefiner2DUDT), allocatable    :: GGTMlinerefiner
        type(GGGridUDT)             :: grid 

        ! Initialize
        !===========
        ! Set verbosity
        verbosity  = options%verbosity 

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

        ! Refiner
        GGTMlinerefiner = InitializeGGTMLineRefiner(topomesh, &
            magneticField, vessel, fieldtracer, boundarytracer, &
            poloidalvertexdistributor, radialvertexdistributor, options)

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

        case ('orthogonal')

            call DistributeVerticesOrthogonal(ggtmdata, topomesh, grid, &
                poloidalvertexdistributor, magneticField, streamlinetracer, &
                GGTMlinerefiner, options)

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

        ! Write intermediate file
        call grid%WriteData('grid_after_cellconstruction')

        ! Add interconnections
        ![grid] = ComputeGridInterconnections(grid);

        ! Process grid
        !=============

        ! Extract the necessary gridding data
        !====================================
        ! Extract
        call ExtractSimulationGrid(simgrid, grid, magneticField)



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
        real(R8), allocatable, dimension(:)         :: dlcv

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
            nv = size(facedata(i)%line%xv) - 2 

            ! Set ID
            tvID = [face%vert(i, 1), (k, k = vertID+1, vertID+nv), face%vert(i, 2)]
            call facedata(i)%line%AddVertexIDs(tvID)

            ! Update line data
            call facedata(i)%line%UpdateLineData(topomesh, ggtmdata)

            ! Update 
            vertID = vertID + nv
            grid%vert%ntot = grid%vert%ntot + nv
        end do

        ! Add cell vertices
        !==================
        do i = 1, cell%ntot

            ! Get cell starting and ending radial line vertices
            srfvID = facedata(celldata(i)%srf)%line%vert
            if (celldata(i)%flipsrf) then 
                srfvID = srfvID(size(srfvID):1:-1)
            end if 
            erfvID = facedata(celldata(i)%erf)%line%vert
            if (celldata(i)%fliperf) then 
                erfvID = erfvID(size(erfvID):1:-1)
            end if 

            ! Distribute lines
            do j = 1, size(celldata(i)%lines)
                ! Distribute over line
                call vd%DistributeOverCurve(celldata(i)%lines(j)%xl, &
                    celldata(i)%lines(j)%yl, nv, ldistr=dlcv)
                call celldata(i)%lines(j)%AddVertexCoordinates(dlcv)

                ! Set vertex ID
                call celldata(i)%lines(j)%AddVertexIDs([srfvID(j+1), &
                    (k, k = vertID+1, vertID+nv-2), erfvID(j+1)])

                ! Update line data
                call celldata(i)%lines(j)%UpdateLineData(topomesh, &
                    ggtmdata)

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

    ! Orthogonal gridder
    subroutine DistributeVerticesOrthogonal(ggtmdata, topomesh, &
        grid, vd, magneticField, streamlinetracer, GGTMlinerefiner, options)

        ! Description
        !============
        ! Construct the grid in an orthogonal way by gridding cells in 
        ! a specific order. This order is determined as follows:
        ! - a cell can be gridded if its high field line is fully 
        !   gridded OR if the cell is a boundary cell where the high
        !   field line is the grid boundary
        ! - after each cell is gridded, the low field boundary vertex
        !   distribution is regenerated based on the previous one, and
        !   propagated to the faces. 
        ! Note that we always grid from high field to low field. 

        ! The vertex distribution is therefore determined by 
        ! the initial vertex distribution on the high field faces. 
        ! Refinement/coarsening can be done in principle in different 
        ! ways, but here we simply take a maximal length distribution
        ! and deduce from that a minimal desired length as well. 

        ! Note 1: we exploit the fact that the topological mesh vertices
        ! are added first, so we can easily check if a vertex is a 
        ! topological vertex by checking if ID <= topomesh.vert.ntot

        ! Note 2: we assume that the streamline tracer is based on 
        ! gradient data of the magnetic field, which points from low
        ! to high value (therefore, we need to trace in the backward
        ! direction since we go from high to low)

        ! Note 3: we don't explicitly keep track of all deleted vertices,
        ! as this would be cumbersome. Instead, at the end we check which
        ! vertex IDs are still present and remap those such that the 
        ! numbering goes from 1 to grid%vert%ntot again. 

        ! Modules
        !========
        use mod_definitions, only: TMfacepolID, TMfacesepID, &
            TMvertexbndID, TMvertextp1ID, TMvertextp2ID


        ! Declare variables
        !==================
        ! Arguments
        class(GGTMDataUDT)                          :: ggtmdata 
        class(TopomeshUDT), intent(in)              :: topomesh 
        class(GGGridUDT), intent(inout)             :: grid 
        class(VertexDistributor2DUDT), intent(in)   :: vd
        type(MagneticFieldUDT), intent(in)          :: magneticField
        class(StreamlineTracerUDT), intent(in)      :: streamlinetracer
        class(GGTMLineRefiner2DUDT), intent(in)     :: GGTMlinerefiner
        type(GGoptionsUDT), intent(in)              :: options

        ! Auxiliary
        integer(I8)                                 :: nv, vertID, tc
        integer(I8), allocatable, dimension(:)      :: tvID, srfvID, &
            erfvID, s11, s12, s21, s22, s31, s32, s41, s42, &
            stype, sortind, newID, updatedfaces, allIDs, vertmap
        logical                                     :: addpoint
        logical, allocatable, dimension(:)          :: iscelldone, &
            isfacedone, isstartingcell, isstartingface, istopoverthf, &
            istopovertlf, isvertexdeleted, keepvert
        real(R8)                                    :: xb(1:2), yb(1:2)
        real(R8), allocatable, dimension(:)         :: xt, yt, &
            x1, x2, x3, x4, y1, y2, y3, y4, s11r, s12r, s21r, s22r, &
            s31r, s32r, s41r, s42r, s1r, temps2r, tempx, tempy, newtx, &
            newty, news2r, newdlcv

        type(GGTMFieldlineDataUDT), allocatable     :: tclines(:)
        type(StreamlineUDT), allocatable            :: orthlines(:)

        ! Loop
        integer(I8)                                 :: i, j, k, nnew

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

        ! Initialize
        allocate(iscelldone(cell%ntot), isfacedone(face%ntot), &
            isstartingcell(cell%ntot), isstartingface(face%ntot))
        iscelldone = .false. 
        isfacedone = .false. 
        isstartingcell = .false.
        isstartingface = .false.

        ! Add topological mesh vertices
        !==============================
        ! Keep track of topological mesh vertices
        ! Set vertID
        vertID = vert%ntot

        ! Preprocess
        !===========
        ! Determine which cells and faces can be started from 
        do i = 1, cell%ntot 
            ! Get cell high field faces and vert
            associate(&
                hffaces     => celldata(i)%hffaces,     &
                hfvert      => celldata(i)%hfvert       &
                )

            ! Check if only one vertex -> check vertex type
            if (size(hfvert) == 1) then 
                if (.not. allocated(celldata(i)%hfvert)) then 
                    call gdErrorHandler('DistributeVerticesOrthogonal: ' // & 
                        'hfvert is not yet allocated')
                end if 
                if (.not. all(IsTopomeshVert(hfvert, topomesh))) then 
                    call gdErrorHandler('DistributeVerticesOrthogonal: ' // & 
                        'assumed topomesh vertex is not a topomesh vertex')
                end if 
                if (any(hfvert > size(vert%type))) then 
                    call gdErrorHandler('d')
                end if 
                ! Sanity check
                ! print *, hfvert 
                if (any(vert%type(hfvert(1)) == [TMvertextp1ID, TMvertextp2ID, &
                    TMvertexbndID, TMvertexmaxID, TMvertexminID])) then 
                    isstartingcell(i) = .true. 
                else
                    ! This shouldn't happen, unless we missed some kind of 
                    ! exception case. Throw error
                    call gdErrorHandler('DistributeVerticesOrthogonal: ' // & 
                        'expected single vertex to be boundary vertex, ' // & 
                        'but is not boundary vertex. Check topological mesh')
                end if
                
                ! Skip remainder of the loop
                cycle
            end if 

            ! If we got here, there should be at least one face. 
            do j = 1, size(hffaces)
                if (face%BF(hffaces(j)) .and. any(face%type(hffaces(j)) == [TMfacepolID, TMfacesepID])) then 
                    ! This is a potential starting face
                    isstartingface(hffaces(j)) = .true. 
                end if 
            end do 

            ! Check if we can start from this cell
            if (all(isstartingface(hffaces))) then 
                isstartingcell(i) = .true. 
            end if 

            ! Housekeeping
            end associate 
        end do

        ! Add face vertices
        !==================
        ! Do for all, but only refine starting faces! 
        do i = 1, face%ntot
            ! Compute number of new vertices
            nv = size(facedata(i)%line%xv) - 2 

            ! Set ID
            tvID = [face%vert(i, 1), (k, k = vertID+1, vertID+nv), face%vert(i, 2)]

            ! Update 
            vertID = vertID + nv

            ! Set data
            call facedata(i)%line%AddVertexIDs(tvID)
            call facedata(i)%line%UpdateLineData(topomesh, ggtmdata)

            ! Check if we should refine
            if (isstartingface(i)) then 
                ! Refine
                keepvert = IsTopomeshVert(facedata(i)%line%vert, topomesh)
                call GGTMlinerefiner%Refine(facedata(i)%line, vertID, keepvert)

                ! Update
                call facedata(i)%line%UpdateLineData(topomesh, ggtmdata)
            end if 
        end do

        ! Add cell vertices
        !==================
        do while (.true.)

            ! Preprocess
            !-----------
            ! Get the next cell 
            tc = findloc((.not. iscelldone) .and. isstartingcell, .true., 1)

            ! Sanity check
            if (tc == 0) then 
                ! This shouldn't be happening
                call gdErrorHandler('DistributeVerticesOrthogonal: ' // & 
                    'could not find next cell - this is either a bug ' // & 
                    'or something is wrong in input')
            end if 

            ! Mark cell as being done
            iscelldone(tc) = .true.

            ! Mark cell faces as being done
            isfacedone(celldata(tc)%hffaces) = .true. 
            isfacedone(celldata(tc)%lffaces) = .true. 

            ! Mark low field faces as potential starting faces
            isstartingface(celldata(tc)%lffaces) = .true. 
            isstartingface(celldata(tc)%hffaces) = .true. 

            ! Recompute potential starting cells - do not overwrite 
            ! original starting cells, as cells with only one high field
            ! vertex are not recomputed here!
            do i = 1, cell%ntot 
                ! Skip if already done
                if (iscelldone(i)) then
                    cycle 
                end if 
                
                ! Check
                if (size(celldata(i)%hffaces) /= 0) then 
                    if (all(isstartingface(celldata(i)%hffaces))) then 
                        isstartingcell(i) = .true. 
                    end if 
                end if 
            end do 

            ! Extract high field line
            call ExtractTMCellAlignedBoundary(celldata(tc), 'high', ggtmdata, &
                topomesh, celldata(tc)%hfline)

            ! Extract low field line (will be overwritten later)
            call ExtractTMCellAlignedBoundary(celldata(tc), 'low', ggtmdata, &
                topomesh, celldata(tc)%lfline) 
            
            ! Unpack
            associate(&
                hfline          => celldata(tc)%hfline,     &
                lfline          => celldata(tc)%lfline,     &
                lines           => celldata(tc)%lines       &
            )

            ! Initialize
            istopoverthf = (hfline%vert /=0 ) .and. (hfline%vert <= vert%ntot)
            istopovertlf = (lfline%vert /=0 ) .and. (lfline%vert <= vert%ntot)

            ! Sanity checks
            if (any(hfline%vert == 0)) then
                ! This shouldn't be happening
                call gdErrorHandler('DistributeVerticesOrthogonal: '// & 
                    'high field line as vertex IDs that are zero, which ' // & 
                    'indicates a cell has been taken of which the high ' // &
                    'field line was not yet fully gridded. This is a bug')
            end if 

            ! Get cell starting and ending radial line vertices
            srfvID = facedata(celldata(tc)%srf)%line%vert
            if (celldata(tc)%flipsrf) then 
                srfvID = srfvID(size(srfvID):1:-1)
            end if 
            erfvID = facedata(celldata(tc)%erf)%line%vert
            if (celldata(tc)%fliperf) then 
                erfvID = erfvID(size(erfvID):1:-1)
            end if 

            ! Determine cell vertices
            !------------------------
            ! Concatenate lines for ease
            tclines = [celldata(tc)%hfline, celldata(tc)%lines, celldata(tc)%lfline]

            ! Compute intersections
            do i = 2, size(tclines)
                ! Skip if it is a tangency point
                if (size(tclines(i)%xl) == 1) then 
                    cycle 
                end if 

                ! Trace field lines starting from previous distribution
                xt = tclines(i-1)%xv
                yt = tclines(i-1)%yv 
                xb = [minval([tclines(i)%xl, tclines(i-1)%xl]), &
                    maxval([tclines(i)%xl, tclines(i-1)%xl])]
                yb = [minval([tclines(i)%yl, tclines(i-1)%yl]), &
                    maxval([tclines(i)%yl, tclines(i-1)%yl])]
                orthlines = streamlinetracer%TraceStreamlines(xt, yt, &
                    xb, yb, spread(-1_I8, 1, size(xt))) ! normally, gradient goes from low to high, so need to reverse sign

                ! Initialize potential new coordinates
                allocate(newtx(size(orthlines)), newty(size(orthlines)), &
                    news2r(size(orthlines)), newID(size(orthlines)))
                nnew = 0

                ! Find intersections with all other boundaries
                !$omp parallel default(private) shared(orthlines, vertID, tclines, newtx, nnew, newty, news2r, newID, i) 
                !$omp do 
                do j = 1, size(orthlines)
                    ! Intersections with starting boundary (should only
                    ! intersect in first point)
                    if (tclines(i-1)%nv == 1) then 
                        ! Previous boundary was point - start and end
                        ! should be the same
                        if ((tclines(i-1)%xv(1) == orthlines(j)%x(1)) .and. &
                            (tclines(i-1)%xv(1) == orthlines(j)%x(1))) then 
                            x1 = tclines(i-1)%xv
                            y1 = tclines(i-1)%yv 
                            allocate(s11(size(x1)), s12(size(x1)), &
                                s11r(size(x1)), s12r(size(x1)))
                            s11 = 1_I8
                            s12 = 1_I8
                            s11r = 0_R8
                            s12r = 0_R8
                        else 
                            call gdErrorHandler('DistributeVerticesOrthogonal: ' // & 
                                'starting point of radial line should be ' // &
                                'the same as tangency point')
                        end if 
                    else
                        call SimplePolygonIntersections(orthlines(j)%x, &
                        orthlines(j)%y, tclines(i-1)%xl, tclines(i-1)%yl, &
                        x1, y1, s11, s12, s11r, s12r)
                    end if 
                    
                    ! Intersections with ending boundary
                    call SimplePolygonIntersections(orthlines(j)%x, &
                        orthlines(j)%y, tclines(i)%xl, tclines(i)%yl, &
                        x2, y2, s21, s22, s21r, s22r)
                    
                    ! Intersections with side boundary 1
                    call SimplePolygonIntersections(orthlines(j)%x, &
                        orthlines(j)%y, [tclines(i-1)%xl(1), tclines(i)%xl(1)], &
                        [tclines(i-1)%yl(1), tclines(i)%yl(1)], &
                        x3, y3, s31, s32, s31r, s32r)
                    
                    ! Intersections with side boundary 2
                    call SimplePolygonIntersections(orthlines(j)%x, &
                        orthlines(j)%y, [tclines(i-1)%xl(tclines(i-1)%nl), tclines(i)%xl(tclines(i)%nl)], &
                        [tclines(i-1)%yl(tclines(i-1)%nl), tclines(i)%yl(tclines(i)%nl)], &
                        x4, y4, s41, s42, s41r, s42r)
                    
                    ! Sort
                    stype = [spread(1_I8, 1, size(x1)), &
                        spread(2_I8, 1, size(x2)), spread(3_I8, 1, size(x3)), &
                        spread(4_I8, 1, size(x4))]
                    s1r = [s11r, s21r, s31r, s41r]
                    allocate(sortind(size(s1r)))
                    sortind = 0_I8
                    call Sort(s1r, ind=sortind, ascend=.true.)
                    stype = stype(sortind)
                    
                    ! Checks:
                    ! - The first intersection should be in the point
                    ! itself
                    ! - The second intersection should be with the
                    ! ending boundary
                    addpoint = .true.
                    if (size(stype) < 2) then 
                        ! Only one intersection found - don't add
                        addpoint = .false.
                    elseif (stype(1) /= 1 .or. s11r(1) /= 0) then ! .or. s1(1) /= 0
                        ! We expect that the first point is an
                        ! intersection with the first boundary
                        addpoint = .false.
                    elseif (stype(2) /= 2) then 
                        ! The second point should intersect with the
                        ! second boundary
                        addpoint = .false.
                    end if 
                    
                    ! Add the point if allowed
                    
                    if (addpoint) then 
                        
                        
                        ! Get point
                        tempx = [x1, x2, x3, x4]
                        tempy = [y1, y2, y3, y4]
                        temps2r = [s12r, s22r, s32r, s42r]
                        tempx = tempx(sortind)
                        tempy = tempy(sortind)
                        temps2r = temps2r(sortind)
                        
                        
                        !$omp critical
                        ! Update counter
                        nnew = nnew + 1

                        ! Add
                        newtx(nnew) = tempx(2)
                        newty(nnew) = tempy(2)
                        news2r(nnew) = temps2r(2)
                        newID(nnew) = vertID+1
                        vertID = vertID+1
                        !$omp end critical
                    end if
                    
                    ! Housekeeping
                    deallocate(s11, s12, s11r, s12r, sortind)
                    
                end do 
                !$omp end do
                !$omp end parallel
                
                ! Trim
                newtx = newtx(1:nnew)
                newty = newty(1:nnew)
                news2r = news2r(1:nnew)
                newID = newID(1:nnew)

                ! Add topological mesh vertices for last line (make sure 
                ! to include first and last vertex...)
                allocate(newdlcv(size(news2r)))
                newdlcv = 0_R8
                if (i == size(tclines)) then 
                    call Interpolate1D(news2r, newdlcv, &
                        real([(k, k = 0, tclines(i)%nl-1)], kind=R8), tclines(i)%dllc)
                    newdlcv = [newdlcv, pack(tclines(i)%dlcv, istopovertlf)]
                    newID = [newID, pack(tclines(i)%vert, istopovertlf)]
                else
                    call Interpolate1D(news2r, newdlcv, &
                        real([(k, k = 0, tclines(i)%nl-1)], kind=R8), tclines(i)%dllc)
                    newdlcv = [0.0_R8, newdlcv, tclines(i)%dllc(tclines(i)%nl)]
                    newID = [srfvID(i), newID, erfvID(i)]
                end if 

                ! Sort
                allocate(sortind(size(newdlcv)))
                sortind = 0_I8
                call Sort(newdlcv, ind=sortind, ascend=.true.)
                newID = newID(sortind)
                
                ! Add points
                call tclines(i)%AddVertexCoordinates(newdlcv, ecbased=.false.)
                call tclines(i)%AddVertexIDs(newID)

                ! Refine/coarsen
                keepvert = IsTopomeshVert(tclines(i)%vert, topomesh) ! keep vertices if topomesh vert
                keepvert(1) = .true. ! keep first and last vertex anyway
                keepvert(size(keepvert)) = .true. 
                call GGTMLinerefiner%Refine(tclines(i), vertID, keepvert)

                ! Update line data (face labels, facedata if last line, ...)
                call tclines(i)%UpdateLineData(topomesh, ggtmdata)

                ! Update GGTM data
                call tclines(i)%UpdateGGTMData(topomesh, ggtmdata, updatedfaces)

                ! Housekeeping
                deallocate(newtx, newty, newID, news2r, sortind, newdlcv)

            end do

            ! Extract line data again
            !------------------------
            ! Extract high field line
            call ExtractTMCellAlignedBoundary(celldata(tc), 'high', ggtmdata, &
                topomesh, celldata(tc)%hfline)

            ! Extract low field line (will be overwritten later)
            call ExtractTMCellAlignedBoundary(celldata(tc), 'low', ggtmdata, &
                topomesh, celldata(tc)%lfline) 

            ! Set cell line data
            celldata(tc)%lines = tclines(2:size(tclines)-1)

            ! Housekeeping
            end associate

            ! Check termination 
            !------------------
            if (all(iscelldone)) then 
                exit 
            end if
        end do 
        
        ! Post-process
        !=============
        ! Still need to account for deleted/overwritten vertices
        allocate(isvertexdeleted(vertID))
        isvertexdeleted = .true. 
        do i = 1, cell%ntot 
            ! Check if vertices are present
            isvertexdeleted(celldata(i)%hfline%vert) = .false.
            isvertexdeleted(celldata(i)%lfline%vert) = .false.
            do j = 1, size(celldata(i)%lines) 
                isvertexdeleted(celldata(i)%lines(j)%vert) = .false.
            end do
        end do 

        ! Get all IDs and construct mapping
        allIDs = pack([(k, k = 1, vertID)], .not. isvertexdeleted)
        allocate(vertmap(vertID))
        vertmap = 0_I8
        vertmap(allIDs) = [(k, k = 1, count(.not. isvertexdeleted))]

        ! Check if any topological mesh vertices were deleted (should not happen)
        if (any(isvertexdeleted(1:topomesh%vert%ntot))) then 
            call gdErrorHandler('DistributeVerticesOrthogonal: ' // & 
                'topological mesh vertices were deleted, this is a bug')
        end if 

        ! Loop and adjust IDs
        do i = 1, cell%ntot 
            ! Remap
            celldata(i)%hfline%vert = vertmap(celldata(i)%hfline%vert)
            celldata(i)%lfline%vert = vertmap(celldata(i)%lfline%vert)
            do j = 1, size(celldata(i)%lines) 
                celldata(i)%lines(j)%vert = vertmap(celldata(i)%lines(j)%vert)
            end do
        end do 

        ! Update number of grid vertices
        grid%vert%ntot = count(.not. isvertexdeleted)

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

        ! Note 3: during construction, duplicate faces will exist since
        ! we don't check whether the start/end line are already 
        ! treated. We remove these faces afterwards. 

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

                ! First line
                if (n1 > 2) then ! Hedge for point line
                    if (l1%vert(1) == l1%vert(n1)) then 
                        ! Closed 
                        call SimplePolygonIntersections(&
                            l1%xv(2:n1-1), l1%yv(2:n1-1), [l1%xv(1), l2%xv(1)], &
                            [l1%yv(1), l2%yv(1)], xint, yint, s1, s2)
                    else
                        call SimplePolygonIntersections(&
                            l1%xv(2:), l1%yv(2:), [l1%xv(1), l2%xv(1)], &
                            [l1%yv(1), l2%yv(1)], xint, yint, s1, s2)
                    end if
                    if (size(xint) > 0) then 
                        isintersectingl1 = .true.
                        mink1 = s1(1)+2
                    end if
                end if 

                ! Second line
                if (n2 > 2) then 
                    if (l2%vert(1) == l2%vert(n2)) then 
                        ! Closed 
                        call SimplePolygonIntersections(&
                            l2%xv(2:n2-1), l2%yv(2:n2-1), [l1%xv(1), l2%xv(1)], &
                            [l1%yv(1), l2%yv(1)], xint, yint, s1, s2)
                    else
                        call SimplePolygonIntersections(&
                            l2%xv(2:), l2%yv(2:), [l1%xv(1), l2%xv(1)], &
                            [l1%yv(1), l2%yv(1)], xint, yint, s1, s2)
                    end if 
                    if (size(xint) > 0) then 
                        isintersectingl2 = .true.
                        mink2 = s1(1)+2
                    end if
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

                ! Add to grid
                allocate(tempcellregion(nc), tempfaceregion(nf))
                tempcellregion = celldata(i)%region
                tempfaceregion = celldata(i)%region
                call grid%AddFace(tempfacevert(1:nf, :), tempfacelabels(1:nf), tempfaceregion(1:nf))
                call grid%AddCell(tempcellvert(1:ncv), tempcellvertP(1:nc, 1:2), tempcellregion(1:nc))

                ! Housekeeping
                deallocate(tempfacevert, tempfacelabels, tempcellvert, &
                    tempcellvertP, tempcellregion, tempfaceregion)
                end associate
            end do 

        end do 

        ! Cleanup
        !========
        call RemoveDuplicateGridFaces(grid)

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Removal of duplicate faces in grid
    subroutine RemoveDuplicateGridFaces(grid)

        ! Description
        !============
        ! This routine removes duplicate faces in the grid which may 
        ! originate from the cell and face constructor(s). Here, we 
        ! simply loop over the faces and mark faces for deletion if 
        ! they have already been found. 

        ! Note: this routine should be called directly after grid cell
        ! and face construction, not later
    

        ! Declare variables
        !==================
        ! Arguments
        type(GGGridUDT), intent(inout)          :: grid 

        ! Auxiliary
        integer(I8)                             :: vtemp
        integer(I8), allocatable, dimension(:)  :: v1, v2, fID, &
            sortind, dv1, dv2, tsortind, tv2, tfID
        integer(I8), allocatable                :: fvert(:, :)
        logical, allocatable, dimension(:)      :: delind 
        
        ! Loop
        integer(I8)                             :: i, k, kold

        ! Initialize
        !===========
        ! Unpack
        associate(f       => grid%face)
        
        v1 = f%v1%Get()
        v2 = f%v2%Get()
        allocate(fvert(f%ntot, 2))
        fvert(:, 1) = v1
        fvert(:, 2) = v2

        ! Remove faces
        !=============
        ! Get face indices
        fID = [(k, k = 1, f%ntot)]

        ! Sort along rows
        do i = 1, size(v1)
            if (v1(i) > v2(i)) then 
                vtemp = v1(i)
                v1(i) = v2(i)
                v2(i) = vtemp
            end if 
        end do

        ! Sort faces according to first vertex
        allocate(sortind(size(v1)))
        call Sort(v1, ind=sortind)
        v2 = v2(sortind)
        fID = fID(sortind)

        ! Sort v2 per segment of equal v1
        k = 0
        do while (k < size(v1))
            ! Update kold
            kold = k 

            ! Get segment with all same values of v1
            k = findloc(v1(kold+1:) /= v1(kold+1), .true., 1, back=.false.)

            ! Hedge for end effects
            if (k == 0) then 
                k = size(v1)
            else
                k = k + kold - 1
            end if 

            ! Sort this segment
            tfID = fID(kold+1:k)
            tv2 = v2(kold+1:k)
            allocate(tsortind(size(tfID)))
            call Sort(tv2, ind=tsortind)
            tfID = tfID(tsortind)
            fID(kold+1:k) = tfID 
            v2(kold+1:k) = tv2
            deallocate(tsortind)
        end do 

         ! Check differences
        dv1 = v1(2:size(v1)) - v1(1:size(v1)-1)
        dv2 = v2(2:size(v2)) - v2(1:size(v2)-1)
        
        ! If both are zero, set delind to true
        delind = [.false., (dv1 == 0) .and. (dv2 == 0)]

        ! Remove faces
        fID = pack(fID, delind)
        call RemoveGGFace(grid, fID)

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
        real(R8), allocatable, dimension(:)     :: dlcv 

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
                ! Initialize
                call facedata(i)%line%Initialize(face%x(i)%Get(), face%y(i)%Get(), face%fsID(i))
                
                ! Distribute
                call vd%DistributeOverCurve(face%x(i)%Get(), face%y(i)%Get(), &
                    nv, ldistr=dlcv)
                
                ! Adjust start and end to be sure
                dlcv(1) = 0
                dlcv(size(dlcv)) = facedata(i)%line%dllc(size(facedata(i)%line%dllc))
                
                ! Add vertex coordinates
                call facedata(i)%line%AddVertexCoordinates(dlcv)
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
        real(R8), allocatable, dimension(:)     :: xc, yc, dlcv

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

                ! Initialize
                call facedata(i)%line%Initialize(xc, yc, face%fsID(i))

                ! Distribute
                call vd%DistributeOverField(xc, yc, field, nv,  ldistr=dlcv)
                dlcv(1) = 0
                dlcv(size(dlcv)) = facedata(i)%line%dllc(size(facedata(i)%line%dllc))

                ! Add data
                call facedata(i)%line%AddVertexCoordinates(dlcv)

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
                nfl = size(facedata(tf(j))%line%xv)
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
            indsrf, tf, cind, nc, ntf, incr, nv, temp, minind, maxind, &
            startind, endind, nfs
        integer(I8), allocatable, dimension(:)  :: tubec, tubef, tcf, &
            tcv, tcfv1, tcfv2, hffaces, lffaces, hfvert, lfvert, &
            allIDs, s1, s2, polv, tf1, tf2, fsID
        integer(I8), allocatable, dimension(:, :)   :: nint, segrf, &
            segc, vertexID
        real(R8)                                :: hfval, lfval, &
            dhf1, dhf2, dlf1, dlf2, nxc, nyc, nxfv, nyfv
        real(R8), allocatable, dimension(:)     :: tcvfval, tcfv1val, &
            tcfv2val, tx, ty, xl, yl, tfval, sr1, sr2, txint, tyint, &
            nxf, nyf, nnf, tsegrc
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

                if ((size(hffaces)+1 /= size(hfvert)) .or. (size(lffaces)+1 /= size(lfvert))) then 
                    print *, 'unexpected'
                end if 
                
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
        ! Initialize flux surface counter
        nfs = topomesh%nfs

        ! Loop over all tubes
        do i = 1, tube%ntot 

            ! Initialize
            !-----------
            tubef = tube%GetFace(i)
            tubec = tube%GetCell(i)

            ! Trace contours
            !---------------
            ! Get initial vertex distribution (skip start and end points)
            tf = tubedata(i)%distributionface
            tx = facedata(tf)%line%xv
            ty = facedata(tf)%line%yv

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
                nxf = -(facedata(tf)%line%yv(2:) - facedata(tf)%line%yv(1:size(facedata(tf)%line%yv)-1))
                nyf = (facedata(tf)%line%xv(2:) - facedata(tf)%line%xv(1:size(facedata(tf)%line%xv)-1))
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

            ! Loop over all contours
            call wall_time(tstart)
            !omp parallel do private(i, j, txint, tyint, s1, s2, sr1, sr2)
            do j = 1, size(tempc)
                ! Compute intersections with each radial face's polygon
                do k = 1, size(tubef)
                    ! Compute intersections
                    call SimplePolygonIntersections(tempc(j)%x, tempc(j)%y, &
                        face%x(tubef(k))%Get(), face%y(tubef(k))%Get(), &
                        txint, tyint, s1, s2, sr1, sr2)
                    !call PolygonIntersections(polc(j), face%pol(tubef(k)), &
                    !    txint, tyint, s1, s2, sr1, sr2)

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
            !omp end parallel do
            call wall_time(tend)
            !print *, 'time spent in intersections:', tend-tstart
            ! Print
            !call tempps%Construct(polc)
            !call tempps%WriteData('lines')
            !call tempps%Construct(face%pol(tubef))
            !call tempps%WriteData('faces')

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
                if (verbosity > 0) then 
                    print *, 'AddTopologicalMeshCellGriddingData: ' // & 
                        'field lines were removed since they do not intersect ' // & 
                        'with one or more radial lines for tube: ', i 
                end if 
            end if 
            if (isintersectremoved) then 
                if (verbosity > 0) then 
                    print *, 'AddTopologicalMeshCellGriddingData: ' // & 
                        'multiple intersections found with some radial lines ' // & 
                        'for tube: ', i, ', taking first intersection'
                end if 
            end if  

            ! Unpack intersections
            allocate(xint(nc, ntf), yint(nc, ntf), segc(nc, ntf), &
                segrf(nc, ntf), segrc(nc, ntf), segrrf(nc, ntf), &
                vertexID(nc, ntf), fsID(nc))
            cc = 0 
            do j = 1, size(nint, 1) 
                ! Skip
                if (keepind(j)) then 
                    ! Update counter
                    cc = cc + 1 

                    ! First intersection: should always be the one
                    ! that is at the start of the contour IF it is a 
                    ! closed polygon! Otherwise, we just take the first one...
                    k = 1
                    if (tube%isclosed(i)) then 
                        tsegrc = segrcda(j, k)%Get()
                        startind = findloc(tsegrc, 0.0_R8, 1)
                        if (startind == 0) then 
                            print *, 'AddToplogicalMeshCellGriddingData: ' // & 
                                'contour does not intersect at starting point. ' // & 
                                'contour: ', cc, 'tube: ', i 
                            print *, 'taking first intersection...'
                            startind = 1
                        end if 
                    else
                        startind = 1
                    end if 
                    xint(cc, k) = xintda(j, k)%Get(startind)
                    yint(cc, k) = yintda(j, k)%Get(startind)
                    segc(cc, k) = segcda(j, k)%Get(startind)
                    segrf(cc, k) = segrfda(j, k)%Get(startind)
                    segrc(cc, k) = segrcda(j, k)%Get(startind)
                    segrrf(cc, k) = segrrfda(j, k)%Get(startind)

                    ! Set vertex ID
                    nv = nv + 1
                    vertexID(cc, k) = nv

                    ! Set flux surface ID
                    fsID(cc) = nfs + 1
                    nfs = nfs + 1

                    ! Loop
                    do k = 2, ntf-1 
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

                    ! Last intersection - can be anything with open 
                    ! surface, but not with closed surface
                    endind = 1
                    
                    ! Hedge for closed flux surfaces
                    if (tube%isclosed(i)) then 
                        ! Intersection with first and last radial line
                        ! should be exactly the same! Only intersection
                        ! coordinate should differ
                        tsegrc = segrcda(cc, ntf)%Get()
                        endind = findloc(tsegrc, real(size(tempc(j)%x)-1, kind=R8), 1)
                        if (endind == 0) then 
                            print *, 'AddToplogicalMeshCellGriddingData: ' // & 
                                'closed contour does not intersect at ending point. ' // & 
                                'contour: ', cc, 'tube: ', i 
                            print *, 'taking first intersection...'
                            endind = 1
                        end if 

                        ! Set vertex ID - should be the same as before
                        nv = nv + 1
                        vertexID(cc, ntf) = vertexID(cc, 1)
                    else

                        ! Set vertex ID
                        nv = nv + 1
                        vertexID(cc, ntf) = nv
                    end if 
                    xint(cc, ntf) = xintda(j, ntf)%Get(endind)
                    yint(cc, ntf) = yintda(j, ntf)%Get(endind)
                    segc(cc, ntf) = segcda(j, ntf)%Get(endind)
                    segrf(cc, ntf) = segrfda(j, ntf)%Get(endind)
                    segrc(cc, ntf) = segrcda(j, ntf)%Get(endind)
                    segrrf(cc, ntf) = segrrfda(j, ntf)%Get(endind)

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
                    if (segrc(j, size(segrc, 2)) /= size(tempc(j)%x)-1) then 
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
                    !polv = [(cc, cc = 1, size(tempc(j)%x))]
                    if (incr > 0) then 
                        polv = [(cc, cc = segc(j, k)+2, segc(j, k+1)-1, incr)]
                    else 
                        polv = [(cc, cc = segc(j, k)-2, segc(j, k+1)+1, incr)]
                    end if 
                    !polv = polc(j)%vert(segc(j, k)+2:segc(j, k+1)-1:incr)
                    xl = [xint(j, k), tempc(j)%x(polv), xint(j, k+1)]
                    yl = [yint(j, k), tempc(j)%y(polv), yint(j, k+1)]

                    ! Add 
                    call celldata(tubec(k))%lines(j)%Initialize(xl, yl, fsID(j))

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
                segrcda, segrrfda, fsID)
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
            face%region = ConstructIntegerDynamicArray()

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
        grid%face%ntot = grid%face%ntot + size(facelabel)

    end subroutine

    ! Cell addition
    subroutine AddGGCell(grid, cellvert, cellvertP, cellregion)

        ! Description
        !============
        ! Add cells to the grid (without updating interconnection).
        ! Note: it is assumed that the vertex pointer vp1 is local, 
        ! i.e. that it always starts from one. We adjust this value
        ! to be immediately correct. 

        ! Declare variables
        !==================
        ! Arguments
        class(GGGridUDT)            :: grid 
        integer(I8), intent(in)     :: cellvertP(:, :), cellvert(:), &
            cellregion(:)

        ! Checks
        !=======
        ! Size checks
        if ((size(cellvertP, 1) /= size(cellregion))) then 
            call gdErrorHandler('AddGGcell: incompatible input sizes')
        end if 
        if (size(cellvertP, 2) /= 2) then 
            call gdErrorHandler('AddGGcell: wrong second dimension of cell vertex pointer')
        end if 

        ! Add
        !====
        if (grid%cell%vp1%Size() > 0) then 
            call grid%cell%vp1%Append(cellvertP(:, 1) &
                + grid%cell%vp1%Get(grid%cell%vp1%Size()) & 
                + grid%cell%vp2%Get(grid%cell%vp2%Size()) - 1)
        else
            call grid%cell%vp1%Append(cellvertP(:, 1))
        end if 
        call grid%cell%vp2%Append(cellvertP(:, 2))
        call grid%cell%vert%Append(cellvert)
        call grid%cell%region%Append(cellregion)
        grid%cell%ntot = grid%cell%ntot + size(cellregion)

    end subroutine

    ! Vertex removal
    subroutine RemoveGGVert(grid, delind)

        ! Description
        !============
        ! Remove vertices of the grid (without updating interconnection)

        ! Declare variables
        !==================
        ! Arguments
        class(GGGridUDT)            :: grid 
        integer(I8), intent(in)     :: delind(:)

        ! Remove
        !=======
        call grid%vert%x%Remove(delind)
        call grid%vert%y%Remove(delind)
        call grid%vert%fieldlineID%Remove(delind)
        grid%vert%ntot = grid%vert%x%Size()

    end subroutine

    ! Face removal
    subroutine RemoveGGFace(grid, delind)

        ! Description
        !============
        ! Remove faces to the grid (without updating interconnection)

        ! Declare variables
        !==================
        ! Arguments
        class(GGGridUDT)            :: grid 
        integer(I8), intent(in)     :: delind(:)

        ! Remove
        !=======
        call grid%face%v1%Remove(delind)
        call grid%face%v2%Remove(delind)
        call grid%face%label%Remove(delind)
        call grid%face%region%Remove(delind)
        grid%face%ntot = grid%face%v1%Size()

    end subroutine

    ! Cell removal
    subroutine RemoveGGCell(grid, delind)

        ! Description
        !============
        ! Remove cells to the grid (without updating interconnection).
        ! Note: it is assumed that the vertex pointer vp1 is local, 
        ! i.e. that it always starts from one. We adjust this value
        ! to be immediately correct. 

        ! Declare variables
        !==================
        ! Arguments
        class(GGGridUDT)            :: grid 
        integer(I8), intent(in)     :: delind(:)

        ! Auxiliary
        integer(I8)                 :: lb, ub 
        integer(I8), allocatable    :: vertdel(:)

        ! Loop
        integer(I8)                 :: i, k

        ! Remove
        !=======
        do i = 1, size(delind)
            lb = grid%cell%vp1%Get(delind(i))
            ub = lb  + grid%cell%vp2%Get(delind(i))-1
            vertdel = [(k, k = ub, lb)]
            call grid%cell%vert%Remove(vertdel)
        end do 
        call grid%cell%vp1%Remove(delind)
        call grid%cell%vp2%Remove(delind)
        call grid%cell%region%Remove(delind)
        grid%cell%ntot = grid%cell%vp1%Size()

    end subroutine

    !------------------------------------------------------------------!
    !                      GGTM LINE HANDLING                          !
    !------------------------------------------------------------------!

    ! GGTM line initialization
    subroutine InitializeGGTMFieldLineData(line, xl, yl, fsID)

        ! Description
        !============
        ! This routine initializes the GGTM field line, but only the 
        ! underlying coordinates that remain unchanged during the 
        ! grid generation step. Derivatives such as lengths and 
        ! accumulated length (the length coordinate axis) are computed
        ! here as well. 

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMFieldlineDataUDT)         :: line 
        real(R8), intent(in)                :: xl(:), yl(:)
        integer(I8), intent(in)             :: fsID

        ! Loop
        integer(I8)                         :: i 

        ! Construct
        !==========
        line%xl = xl 
        line%yl = yl 
        line%dll = sqrt((xl(2:) - xl(1:size(xl)-1))**2 &
            + (yl(2:) - yl(1:size(yl)-1))**2)
        line%dllc = xl 
        line%dllc = 0_R8
        do i = 1, size(xl)-1
            line%dllc(i+1) = line%dllc(i) + line%dll(i)
        end do 
        line%nl = size(xl)
        line%nv = 0
        line%fsID = fsID

    end subroutine

    ! GGTM line vertex adder
    subroutine AddVertexCoordinates(line, dlcv, ecbased)

        ! Description
        !============
        ! This routine adds vertices xv, yv on the line by interpolating
        ! the x and y coordinates based on the given length distribution
        ! dlcv. Note that this routinewill overwrite any pre-existing 
        ! vertex distribution, and that the vertex IDs are not yet 
        ! added here. 

        ! Note: we hedge for 'lines' that are literally one vertex long.
        ! In that case, dlcv should equal zero and have size one

        ! Note: if provided, the optional input argument can be used to
        ! give dlcv in terms of edge coordinate (going from 0 to nl-1) if
        ! true. Default is that dlcv is in terms of length (going from 0
        ! to line%dllc(end))

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMFieldlineDataUDT)         :: line 
        real(R8), intent(in)                :: dlcv(:)
        logical, intent(in), optional       :: ecbased 

        ! Auxiliary
        real(R8), allocatable, dimension(:)     :: tdlcv

        ! Loop
        integer(I8)                             :: k 

        ! Check
        !======
        if (present(ecbased)) then 
            if (ecbased) then 
                ! Reformat dlcv to be in terms of length
                tdlcv = dlcv ! just initialization
                call Interpolate1D(dlcv, tdlcv, real([(k, k = 0, line%nl-1)], kind=R8), line%dllc)
                
                ! Check
                where (dlcv == 0_R8) tdlcv = 0
                where (dlcv == int((line%nl-1), kind=R8)) tdlcv = line%dllc(line%nl)
            else
                tdlcv = dlcv 
            end if 
        else
            tdlcv = dlcv
        end if 

        ! Add
        !====
        if (size(tdlcv) == 1) then 
            if (tdlcv(1) == 0.0_R8) then 
                line%dlcv = tdlcv 
                line%xv = line%xl(1:1)
                line%yv = line%yl(1:1)
                line%nv = 1
                return 
            end if 
        end if 

        ! Hedge for out of bounds points
        line%dlcv = tdlcv 
        where (line%dlcv >= line%dllc(line%nl)) line%dlcv = line%dllc(line%nl)     
        
        if (allocated(line%xv)) then 
            deallocate(line%xv, line%yv)
        end if 
        allocate(line%xv(size(tdlcv)), line%yv(size(tdlcv)))
        call Interpolate1D(line%dlcv, line%xv, line%dllc, line%xl)
        call Interpolate1D(line%dlcv, line%yv, line%dllc, line%yl)
        if (any(.not. ieee_is_finite(line%xv))) then 
            print *, 'NaNs detected'
        end if 
        line%nv = size(tdlcv)

    end subroutine

    ! GGTM line vertex ID setter
    subroutine AddVertexIDs(line, vertIDs)

        ! Description
        !============
        ! ID setter (simple wrapper)

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMFieldlineDataUDT)     :: line 
        integer(I8), intent(in)         :: vertIDs(:)

        ! Set
        !====
        line%vert = vertIDs 

    end subroutine

    ! GGTM line data updater
    subroutine UpdateLineData(line, topomesh, ggtmdata)

        ! Description
        !============
        ! This routine updates additional line data required for later
        ! gridding. This includes face labels, vertex field line IDs, 
        ! etc. It is assumed that vertices are already added for this 
        ! line and that vertex IDs are present. It is implicitly assumed
        ! that all vertices with ID <= topomesh%vert%ntot are topological
        ! mesh vertices. 

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMFieldlineDataUDT)         :: line 
        class(TopomeshUDT), intent(in)      :: topomesh 
        class(GGTMDataUDT), intent(in)      :: ggtmdata

        ! Auxiliary
        logical, allocatable, dimension(:)  :: istopovert 
        integer                             :: floc1, floc2, nf
        integer, allocatable, dimension(:)  :: tv, tvind

        ! Loop
        integer(I8)                         :: i, k, fc  

        ! Initialize
        !===========
        ! Associate
        associate(&
            face        => topomesh%face    &
            )
        ! Check if the line is a tangency point - if so, simply initialize 
        ! data and return
        if (size(line%xv) == 1) then 
            ! Sanity check
            if (line%vert(1) > topomesh%vert%ntot) then 
                call gdErrorHandler('UpdateLineData: line with single ' // & 
                    'vertex is expected to be tangency point, but ' // & 
                    'vertex is not a topological mesh vertex')
            end if 

            ! Set facelabels
            if (allocated(line%facelabels)) then 
                deallocate(line%facelabels)
            end if 
            allocate(line%facelabels(0))

            ! Exit
            return 
        end if 

        ! Check if there are topological mesh vertices
        istopovert = IsTopomeshVert(line%vert, topomesh)

        ! Sanity checks
        if (count(istopovert) == 1) then 
            call gdErrorHandler('UpdateLineData: found only a single ' // & 
                'topological mesh vertex in line, this is unexpected ' // &
                '(should be either 0, 2, or more)')
        end if 

        ! Update
        !=======
        ! Initialize
        line%facelabels = line%vert(1:line%nv-1)
        line%facelabels = 0_I8

        ! Check for topoverts
        if (any(istopovert)) then ! Need to check topological face labels
            ! Initialize
            fc = 0

            ! Get topomesh vertices
            allocate(tv(count(istopovert)), tvind(count(istopovert)))
            tv = pack(line%vert, istopovert)
            tvind = pack([(k, k = 1, size(line%vert))], istopovert)
            
            ! Get topomesh faces
            do i = 1, size(tv)-1
                ! Get amount of line faces
                nf = tvind(i+1) - tvind(i)

                ! Get face index
                floc1 = findloc((face%vert(:, 1) == tv(i)) .and. &
                    (face%vert(:, 2) == tv(i+1)), .true., 1)
                floc2 = findloc((face%vert(:, 2) == tv(i)) .and. &
                    (face%vert(:, 1) == tv(i+1)), .true., 1)                

                ! Sanity check
                if ((floc1 /= 0) .and. (floc2 /= 0)) then 
                    print *, 'vertex 1:', tv(i), 'vertex 2:', tv(i+1)
                    call gdErrorHandler('UpdateLineData: multiple ' // &
                        'topological faces found for vertex pair, check input.')
                elseif ((floc1 == 0) .and. (floc2 == 0)) then 
                    print *, 'vertex 1:', tv(i), 'vertex 2:', tv(i+1)
                    call gdErrorHandler('UpdateLineData: no ' // &
                        'topological faces found for vertex pair, check input.')
                elseif (floc1 /= 0) then 
                    ! Set facelabels equal to face index of topomesh face
                    line%facelabels(fc+1:fc+nf) = spread(floc1, 1, nf)
                else
                    ! Set facelabels equal to face index of topomesh face
                    line%facelabels(fc+1:fc+nf) = spread(floc2, 1, nf)
                end if 

                ! Update counter
                fc = fc + nf
            end do 

            ! Housekeeping
            deallocate(tv, tvind)

        end if 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! GGTM data updating
    subroutine UpdateGGTMData(line, topomesh, ggtmdata, adjustedfaces)

        ! Description
        !============
        ! This subroutine updates the GGTM data. More precisely, it 
        ! updates the vertex distribution of any faces that the 
        ! current line may hold. To determine these faces, a similar
        ! approach as in UpdateLineData is taken (i.e. also here we 
        ! assume that if vertID <= topomesh%vert%ntot, the vertex is 
        ! a topological mesh vertex). Additionally, this routine returns
        ! the topological mesh face indices that have been reworked as 
        ! an output argument. 

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMFieldlineDataUDT)             :: line 
        class(TopomeshUDT), intent(in)          :: topomesh 
        class(GGTMDataUDT), intent(inout)       :: ggtmdata 
        integer(I8), allocatable, intent(out)   :: adjustedfaces(:)

        ! Auxiliary
        logical, allocatable, dimension(:)      :: istopovert 
        integer(I8)                             :: floc1, floc2, nf
        integer(I8), allocatable, dimension(:)  :: tv, tvind
        real(R8), allocatable, dimension(:)     :: dlcvf

        ! Loop
        integer(I8)                         :: i, k 

        ! Initialize
        !===========
        ! Associate
        associate(&
            face        => topomesh%face,       &
            vert        => topomesh%vert,       &
            facedata    => ggtmdata%face        &
            )

        ! Initialize
        if (allocated(adjustedfaces)) then 
            deallocate(adjustedfaces)
        end if 

        ! Check if the line is a tangency point - if so, simply initialize 
        ! data and return
        if (size(line%xv) == 1) then 
            ! Sanity check
            if (line%vert(1) > topomesh%vert%ntot) then 
                call gdErrorHandler('UpdateGGTMData: line with single ' // & 
                    'vertex is expected to be tangency point, but ' // & 
                    'vertex is not a topological mesh vertex')
            end if 

            ! Initialize
            allocate(adjustedfaces(0))

            ! Exit
            return 
        end if 

        ! Check if there are topological mesh vertices. If not, return
        istopovert = IsTopomeshVert(line%vert, topomesh)
        if (count(istopovert) == 0) then 
            allocate(adjustedfaces(0))
            return 
        end if 

        ! Sanity checks
        if (count(istopovert) == 1) then 
            call gdErrorHandler('UpdateGGTMData: found only a single ' // & 
                'topological mesh vertex in line, this is unexpected ' // &
                '(should be either 0, 2, or more)')
            return 
        end if 

        ! Update
        !=======
        ! Get topomesh vertices
        allocate(tv(count(istopovert)), tvind(count(istopovert)))
        tv = pack(line%vert, istopovert)
        tvind = pack([(k, k = 1, size(line%vert))], istopovert)

        ! Get topomesh faces
        allocate(adjustedfaces(size(tv)-1))
        do i = 1, size(adjustedfaces)
            ! Get amount of line faces
            nf = tvind(i+1) - tvind(i)

            ! Get face index
            floc1 = findloc((face%vert(:, 1) == tv(i)) .and. &
                (face%vert(:, 2) == tv(i+1)), .true., 1)
            floc2 = findloc((face%vert(:, 2) == tv(i)) .and. &
                (face%vert(:, 1) == tv(i+1)), .true., 1)                

            ! Sanity check
            if ((floc1 /= 0) .and. (floc2 /= 0)) then 
                print *, 'vertex 1:', tv(i), 'vertex 2:', tv(i+1)
                call gdErrorHandler('UpdateGGTMData: multiple ' // &
                    'topological faces found for vertex pair, check input.')
            elseif ((floc1 == 0) .and. (floc2 == 0)) then 
                print *, 'vertex 1:', tv(i), 'vertex 2:', tv(i+1)
                call gdErrorHandler('UpdateGGTMData: no ' // &
                    'topological faces found for vertex pair, check input.')
            elseif (floc1 /= 0) then 
                ! Reset vertices by recomputing the length distribution
                ! for the face
                dlcvf = line%dlcv(tvind(i):tvind(i+1))
                dlcvf = (dlcvf - dlcvf(1))
                dlcvf(1) = 0
                dlcvf(size(dlcvf)) = facedata(floc1)%line%dllc(facedata(floc1)%line%nl)
                call facedata(floc1)%line%AddVertexCoordinates(dlcvf)
                call facedata(floc1)%line%AddVertexIDs(line%vert(tvind(i):tvind(i+1)))
                
                ! Add to adjusted faces
                adjustedfaces(i) = floc1
            else
                ! Reset vertices by recomputing the length distribution
                ! for the face. Need to flip now
                dlcvf = line%dlcv(tvind(i):tvind(i+1))
                dlcvf(size(dlcvf):1:-1) = &
                    facedata(floc2)%line%dllc(facedata(floc2)%line%nl) - (dlcvf - dlcvf(1))
                dlcvf(1) = 0
                dlcvf(size(dlcvf)) = facedata(floc2)%line%dllc(facedata(floc2)%line%nl)
                call facedata(floc2)%line%AddVertexCoordinates(dlcvf)
                call facedata(floc2)%line%AddVertexIDs(line%vert(tvind(i+1):tvind(i):-1))
            
                ! Add to adjusted faces
                adjustedfaces(i) = floc2
            end if 
        end do 

        ! Housekeeping
        deallocate(tv, tvind)

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Refiner initialization
    function InitializeGGTMLineRefiner(topomesh, &
        magneticField, vessel, fieldtracer, boundarytracer, &
        poloidalvertexdistributor, radialvertexdistributor, options)& 
        result(GGTMlinerefiner)

        ! Description
        !============
        ! This routine initializes the line refiner based on the options
        ! given in the GGoptions structure. Since there may be a lot of
        ! possible desired options to do refinement already in the grid
        ! generation stage, we also pass most available data that we 
        ! have.

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT), intent(in)          :: topomesh 
        type(MagneticFieldUDT), intent(in)      :: magneticField
        type(VesselUDT), intent(in)             :: vessel 
        class(ContourTracerUDT), intent(in)     :: fieldtracer, boundarytracer 
        type(GGoptionsUDT), intent(in)          :: options 
        class(VertexDistributor2DUDT), intent(in)      :: &
            poloidalvertexdistributor, radialvertexdistributor
        class(GGTMLineRefiner2DUDT), allocatable    :: GGTMlinerefiner

        ! Auxiliary
        integer(I8), allocatable, dimension(:, :)   :: labels(:, :)
        logical, allocatable, dimension(:)          :: includevesselvertex
        real(R8), allocatable, dimension(:)         :: xp, yp, &
            valpLmin, valpLmax, decaylength, xv, yv, tempLmin, tempLmax, &
            tempdecaylength
        type(Coordinates2DDistanceDFUDT)            :: Lmin, Lmax

        ! Loop
        integer(I8)                                 :: i 

        ! Select refiner
        !===============
        select case (options%refmeth)

        case ('no')

            ! No refinement
            allocate(GGTMLineRefinerNoRefUDT::GGTMlinerefiner)

        case ('lengthbased')

            ! Length-based refinement
            allocate(GGTMLineRefinerLB2DUDT::GGTMlinerefiner)

        case default 

            ! Unknown
            call gdErrorHandler('InitializeGGTMLineRefiner: unknown ' // & 
                'option: ' // options%refmeth)

        end select

        ! Initialize refiner
        !===================
        select type (GGTMlinerefiner)

        type is (GGTMLineRefinerNoRefUDT)

            ! Do nothing

        type is (GGTMLineRefinerLB2DUDT)

            ! Check how to generate min and max distributions - note: 
            ! currently, we only do point-based distribution methods
        
            ! Check which points to include
            allocate(xp(0), yp(0), valpLmin(0), valpLmax(0), decaylength(0))

            ! Include x-point regions?
            if (options%refLBdoxp) then 
                ! Add all x-points
                do i = 1, topomesh%vert%ntot 
                    if (topomesh%vert%type(i) == TMvertexsaddleID) then 
                        xp = [xp, topomesh%vert%x(i)]
                        yp = [yp, topomesh%vert%y(i)]
                        valpLmin = [valpLmin, options%refLBLminxp]
                        valpLmax = [valpLmax, options%refLBLmaxxp]
                        decaylength = [decaylength, options%refLBdecaylengthxp]
                    end if 
                end do 
            end if 

            ! Include vessel vertices (e.g. targets)?
            if (options%refLBdovessel) then 
                ! Get labels & coordinates
                call vessel%polygonset%GetLabels(labels)
                call vessel%polygonset%GetVertices(xv, yv)
                allocate(includevesselvertex(size(xv)), tempLmin(size(xv)), &
                    tempLmax(size(yv)), tempdecaylength(size(xv)))
                includevesselvertex = .false. 
                tempLmin = 0
                tempLmax = 0
                tempdecaylength = 0

                ! Add per vessel structure
                do i = 1, size(options%refLBstructureIDs)
                    ! Unpack ID
                    associate(tID       => options%refLBstructureIDs(i))

                    ! Check vertices
                    where ( (labels(:, 1) == tID) .or. (labels(:, 2) == tID) ) 
                        includevesselvertex = .true. 
                        tempLmin = options%refLBLminstructure(i)
                        tempLmax = options%refLBLmaxstructure(i)
                        tempdecaylength = options%refLBdecaylengthstructure(i)
                    end where 
        
                    ! Housekeeping
                    end associate
                end do 

                ! Add per separate vessel vertex ID
                do i = 1, size(options%refLBvertIDs)
                    ! Unpack ID
                    associate(tID       => options%refLBvertIDs(i))
        
                    ! Check vertices
                    where( (labels(:, 3) == tID)) 
                        includevesselvertex = .true. 
                        tempLmin = options%refLBLminvert(i)
                        tempLmax = options%refLBLmaxvert(i)
                        tempdecaylength = options%refLBdecaylengthvert(i)
                    end where
        
                    ! Housekeeping
                    end associate
                end do

                ! Include
                xp = [xp, pack(xv, includevesselvertex)]
                yp = [yp, pack(yv, includevesselvertex)]
                valpLmin = [valpLmin, pack(tempLmin, includevesselvertex)]
                valpLmax = [valpLmax, pack(tempLmax, includevesselvertex)]
                decaylength = [decaylength, pack(tempdecaylength, includevesselvertex)]

            end if 

            ! Construct distribution functions
            call Lmin%Initialize(xp, yp, valpLmin, options%refLBLmininf, &
                decaylength)
            call Lmax%Initialize(xp, yp, valpLmax, options%refLBLmaxinf, &
                decaylength)

            ! Construct refiner
            GGTMlinerefiner = ConstructGGTMLineRefiner(Lmin, Lmax, 'classic')


        class default 

            call gdErrorHandler('INitializeGGTMLineRefiner: type not implemented')

        end select


    end function 

    ! Length-based refiner constructor
    function ConstructGGTMLineRefinerLB(Lmin, Lmax, meth) result(refiner)

        ! Description
        !============
        ! Constructor for line refinement, based on minimal and maximal
        ! length distributions.

        ! Declare variables
        !==================
        ! Arguments
        type(GGTMLineRefinerLB2DUDT)                :: refiner 
        class(DistributionFunctionUDT), intent(in)  :: Lmin, Lmax
        character(*), intent(in)                    :: meth 
 
        ! Initialize
        !===========
        refiner%Lmin = Lmin 
        refiner%Lmax = Lmax 
        refiner%meth = meth

    end function

    ! Single line refinement, dummy
    subroutine RefineLineSingleNoRef(refiner, line, vertID, keepvert)

        ! Description
        !============
        ! Simply returns the original distribution

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMLineRefinerNoRefUDT)              :: refiner 
        type(GGTMFieldlineDataUDT), intent(inout)   :: line 
        integer(I8), intent(inout)                  :: vertID 
        logical, intent(in)                         :: keepvert(:)

        ! Checks
        !=======
        ! Ensure proper dimensions
        if (size(keepvert) /= size(line%vert)) then 
            call gdErrorHandler('RefineLineSingleLB: keepvert does not '// & 
                'have same number of elements as line%vert, check input')
        end if 

        ! Return
        !=======

    end subroutine

    ! Single line refinement, length based
    subroutine RefineLineSingleLB(refiner, line, vertID, keepvert)

        ! Description
        !============
        ! Refine/coarsen a line based on a maximal and minimal length
        ! distribution. The refinement/coarsening factor is hard coded 
        ! to be two (i.e. either a vertex is added or deleted). In principle
        ! the minimal length distribution should be at least two times
        ! smaller than the maximal length distribution in each evaluation
        ! to avoid infinite refinement. However, here we solve this 
        ! potential issue by marking faces that have been refinement as 
        ! impossible to coarsen, and vice versa. 

        ! Since vertices may be added/deleted, the current vertex index
        ! vertID should be passed. On exit, this will be updated (if 
        ! vertices were added). The deletedvert logical array will be 
        ! equal to or larger than vertID. If an element is true, this 
        ! means that that vertex ID was deleted. On exit, this array will
        ! be up to date from 1:vertID. It should be used, after all 
        ! vertices have been constructed, to remove unused vertices. 

        ! Note: instead of throwing an error, we allow the minimal length
        ! to be larger than the maximal length. Here, we give priority
        ! to the maximal length (i.e. we will refine rather than coarsen)
        ! since this is often more important for simulations etc. 

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMLineRefinerLB2DUDT)               :: refiner 
        type(GGTMFieldlineDataUDT), intent(inout)   :: line 
        integer(I8), intent(inout)                  :: vertID 
        logical, intent(in)                         :: keepvert(:)

        ! Auxiliary
        logical                                     :: ismerged, &
            isnextfacelegal, isprevfacelegal
        logical, allocatable, dimension(:)          :: isreflegal, &
            iscoarselegal, thiskeepvert, newkeepvert, refineface, &
            coarsenface,  newisreflegal, newiscoarselegal, iscoarsenedface
        real(R8), allocatable, dimension(:)         :: dll, &
            Lmaxvert, Lminvert, newdll, newdllc
        integer(I8), allocatable, dimension(:)      :: thisvertID, &
            newvertID

        ! Loop
        integer(I8)                                 :: i, cc 

        ! Initialize
        !===========
        ! Ensure proper dimensions
        if (size(keepvert) /= size(line%vert)) then 
            call gdErrorHandler('RefineLineSingleLB: keepvert does not '// & 
                'have same number of elements as line%vert, check input')
        end if 

        ! Set initial logicals
        allocate(iscoarselegal(line%nv-1))
        iscoarselegal = .true.
        isreflegal = iscoarselegal
        thiskeepvert = keepvert
        thisvertID = line%vert

        ! Initial distribution
        dll = line%dlcv(2:line%nv) - line%dlcv(1:line%nv-1)
        allocate(newdllc(size(dll)+1))
        newdllc = spread(0.0_R8, 1, size(dll)+1)
        do i = 2, size(newdllc)
            newdllc(i) = newdllc(i-1) + dll(i-1)
        end do 

        ! Loop
        !=====
        do while (.true.)

            ! Precompute
            !-----------
            ! Minimal & maximal length @ vertices
            allocate(Lmaxvert(line%nv), Lminvert(line%nv))
            Lmaxvert = 0.0_R8
            Lminvert = Lmaxvert 
            call refiner%Lmax%Evaluate(line%xv, line%yv, Lmaxvert)
            call refiner%Lmin%Evaluate(line%xv, line%yv, Lminvert)

            ! Issue warning if Lminvert >= Lmaxvert 
            if (verbosity > 1) then 
                if (any(Lminvert >= Lmaxvert)) then 
                    print *, 'RefineLineSingleLB: minimal length exceeds ' // & 
                        'maximal length at some vertices'
                end if 
            end if 

            ! Initialize
            allocate(coarsenface(size(dll)), refineface(size(dll)))
            coarsenface = .false.
            refineface = .false.

            ! Determine which faces to refine/coarsen
            where (((dll > Lmaxvert(1:line%nv-1)) .or. (dll > Lmaxvert(2:line%nv))) &
                .and. isreflegal) 
                refineface = .true. 
                iscoarselegal = .false. 
            end where
            where (((dll < Lminvert(1:line%nv-1)) .or. (dll < Lminvert(2:line%nv))) &
                .and. iscoarselegal &
                .and. .not. (thiskeepvert(1:line%nv-1) .and. thiskeepvert(2:line%nv)) &
                .and. .not. ([.false., .not. iscoarselegal(1:line%nv-2)] .and. [iscoarselegal(2:line%nv-1), .false.]))
                coarsenface = .true.
                isreflegal = .false.
            end where
                
            ! Check exit conditions
            if ((.not. any(refineface)) .and. (.not. any(coarsenface))) then 
                exit
            end if 

            ! Refine/coarsen
            !---------------
            select case (refiner%meth)

            case ('classic')

                ! Classic refinement by splitting face, coarsening by
                ! removing vertex (without any other adaptations)

                ! Initialize new length distribution etc - too long, trimmed later
                allocate(newdll(size(dll)+count(refineface)))
                allocate(newkeepvert(size(newdll)+1), newvertID(size(newdll)+1))
                allocate(newiscoarselegal(size(newdll)), newisreflegal(size(newdll)), &
                    iscoarsenedface(size(newdll)))
                newdll = 0.0_R8
                newkeepvert = .false. 
                newvertID = 0_I8
                newiscoarselegal = .true.
                newisreflegal = .true.
                iscoarsenedface = .false.

                ! Loop
                cc = 1 ! from one to size newdll
                i = 1 ! from one to size dll
                newkeepvert(1) = keepvert(1)
                newvertID(cc) = thisvertID(cc) ! always keep first vertex
                do while (i <= size(dll))
                    
                    ! Refine/coarsen face? 
                    if (refineface(i)) then 
                        ! Split face
                        if (cc+1 > size(newdll)) then 
                            call gdErrorHandler('Something wrong')
                        end if 
                        if (iscoarsenedface(cc)) then 
                            cc = cc + 1
                        end if 
                        newdll(cc:cc+1) = dll(i)/2.0_R8 

                        ! Update next vertices
                        newkeepvert(cc+1) = .false.
                        newvertID(cc+1) = vertID+1
                        newkeepvert(cc+2) = thiskeepvert(i+1)
                        newvertID(cc+2) = thisvertID(i+1)

                        ! Coarsening is illegal
                        newiscoarselegal(cc:cc+1) = .false.

                        ! Update counters
                        vertID = vertID + 1
                        i = i + 1
                        cc = cc + 2
                    elseif (coarsenface(i)) then 
                        ! Preliminary checks
                        isnextfacelegal = (i < size(dll)) 
                        if (isnextfacelegal) then 
                            isnextfacelegal = isnextfacelegal .and. (.not. thiskeepvert(i+1))
                        end if 
                        isprevfacelegal = (cc > 1) 
                        if (isprevfacelegal) then 
                            isprevfacelegal = isprevfacelegal &
                                .and. (.not. newkeepvert(cc)) &
                                .and. (.not. newiscoarselegal(cc-1))
                        end if  
                        
                        ! Check which face to merge with
                        ismerged = .false. 
                        if (isnextfacelegal) then 
                            if (coarsenface(i+1)) then ! both faces mergeable
                                ! Merge faces
                                newdll(cc) = dll(i) + dll(i+1)

                                ! Update next vertex
                                newvertID(cc+1) = thisvertID(i+2)
                                newkeepvert(cc+1) = thiskeepvert(i+2)

                                ! Refinement is illegal
                                newisreflegal(cc) = .false. 
                                iscoarsenedface(cc) = .true. 

                                ! Update counters
                                cc = cc + 1 
                                i = i + 2 !skip next face, cause already merged
                                ismerged = .true.
                            end if 
                        end if 

                        ! If not merged, check previous/next face
                        if (.not. ismerged .and. (isnextfacelegal .and. isprevfacelegal)) then
                            ! Need to check both previous and next face
                            ! Next face may be merged
                            if (newdll(cc-1) < dll(i+1) .or. (refineface(i+1))) then 
                                ! Merge previous face
                                newdll(cc-1) = newdll(cc-1)+dll(i)

                                ! Set next vertex 
                                newvertID(cc) = thisvertID(i+1)
                                newkeepvert(cc) = thiskeepvert(i+1)

                                ! Refinement is illegal
                                newisreflegal(cc-1) = .false. 

                                ! Update counter
                                i = i + 1
                            else
                                ! Merge next face
                                newdll(cc) = dll(i) + dll(i+1)

                                ! Set next vertex
                                newvertID(cc+1) = thisvertID(i+2)
                                newkeepvert(cc+1) = thiskeepvert(i+2)

                                ! Refinement is illegal
                                newisreflegal(cc+1) = .false. 

                                ! Update counters
                                cc = cc + 1
                                i = i + 2
                            end if 

                            ! Anyhow merging should succeed
                            ismerged = .true.
                        end if 
                        if (.not. ismerged .and. (isnextfacelegal)) then 
                            if (.not. refineface(i+1)) then ! only merge if next face isn't refined
                                ! Merge next face
                                newdll(cc) = dll(i) + dll(i+1)

                                ! Set new vertex
                                newvertID(cc+1) = thisvertID(i+2)
                                newkeepvert(cc+1) = thiskeepvert(i+2)

                                ! Refinement is illegal
                                newisreflegal(cc) = .false. 

                                ! Update counters
                                cc = cc + 1
                                i = i + 2
                                ismerged = .true.
                            end if 
                        end if 
                        if (.not. ismerged .and. (isprevfacelegal)) then 
                            ! Merge previous face
                            newdll(cc-1) = newdll(cc-1)+dll(i)

                            ! Set next vertex 
                            newvertID(cc) = thisvertID(i+1)
                            newkeepvert(cc) = thiskeepvert(i+1)

                            ! Refinement is illegal
                            newisreflegal(cc-1) = .false. 

                            ! Update counter
                            i = i + 1

                            ! Update statement
                            ismerged = .true.
                        end if 
                        if (.not. ismerged) then 
                            ! Skip for next loop
                            newdll(cc) = dll(i)
                            newvertID(cc+1) = thisvertID(i+1)
                            newkeepvert(cc+1) = thiskeepvert(i+1)
                            newisreflegal(cc) = isreflegal(i)
                            newiscoarselegal(cc) = iscoarselegal(i)
                            if (.not. isnextfacelegal .and. .not. isprevfacelegal) then 
                                newiscoarselegal(cc) = .false. 
                            end if 
                            cc = cc + 1
                            i = i + 1
                        end if
                    else
                        ! Don't do anything, simply copy
                        newdll(cc) = dll(i)
                        newvertID(cc+1) = thisvertID(i+1)
                        newkeepvert(cc+1) = thiskeepvert(i+1)
                        newisreflegal(cc) = isreflegal(i)
                        newiscoarselegal(cc) = iscoarselegal(i)
                        
                        ! Update counters
                        cc = cc + 1
                        i = i + 1
                    end if 
                end do 

                ! Update
                thiskeepvert = newkeepvert(1:cc)
                thisvertID = newvertID(1:cc)
                isreflegal = newisreflegal(1:cc-1) 
                iscoarselegal = newiscoarselegal(1:cc-1)
                dll = newdll(1:cc-1)
                deallocate(newdllc)
                allocate(newdllc(size(dll)+1))
                newdllc = 0_R8
                do i = 2, size(newdllc)
                    newdllc(i) = dll(i-1) + newdllc(i-1)
                end do 
                call line%AddVertexCoordinates(newdllc)

                ! Housekeeping
                deallocate(newvertID, newdll, newkeepvert, newisreflegal, &
                    newiscoarselegal, iscoarsenedface)

            case default

                call gdErrorHandler('RefineSingleLineLB: unknown ' // & 
                    'refiner method: ' // refiner%meth)

            end select

            ! Housekeeping
            deallocate(coarsenface, refineface, Lminvert, Lmaxvert)

        end do

        ! Construct line coordinates
        !===========================
        call line%AddVertexCoordinates(newdllc)
        call line%AddVertexIDs(thisvertID)

    end subroutine

    !------------------------------------------------------------------!
    !                             OUTPUT                               !
    !------------------------------------------------------------------!

    ! Simulation grid extraction
    subroutine ExtractSimulationGrid(simgrid, grid, magneticField)

        ! Description
        !============
        ! This routine extracts the necessary grid data for the 
        ! classical 'gridudt' type used for simulations/deformation/...
        
        ! Notes
        !======
        ! Note 1: face labels etc are not yet translated to logical
        ! values. This should be done in a separate routine
        ! afterwards. 

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(out)              :: simgrid
        type(GGGridUDT), intent(in)             :: grid 
        type(MagneticFieldUDT), intent(in)      :: magneticField

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: tv
        real(R8), allocatable, dimension(:)     :: bpvx

        ! Loop
        integer(I8)                             :: i 

        ! Initialize
        !===========
        ! Associate
        associate(&
            nv      => grid%vert%ntot,  &
            nf      => grid%face%ntot,  &
            nc      => grid%cell%ntot   &
            )
        

        ! Basic data
        !===========
        ! Sizes
        simgrid%vert%ntot = nv 
        simgrid%face%ntot = nf 
        simgrid%cell%ntot = nc 
        simgrid%cell%nvert = grid%cell%vert%Size()
        simgrid%cell%nface = simgrid%cell%nvert ! should be exactly the same
        simgrid%data%fluxdata%nFs = 0 ! changed later, just for initialization now
        simgrid%data%fluxdata%nFt = 0 ! changed later, just for initialization now

        ! Allocate grid
        call AllocateGrid(simgrid)

        ! Vertices
        simgrid%vert%x              = grid%vert%x%Get()
        simgrid%vert%y              = grid%vert%y%Get()
        simgrid%vert%fieldlineID    = grid%vert%fieldlineID%Get()
        simgrid%vert%ntot           = nv 

        ! Faces
        simgrid%face%vert(:, 1) = grid%face%v1%Get()
        simgrid%face%vert(:, 2) = grid%face%v2%Get()
        simgrid%face%label      = grid%face%label%Get()
        simgrid%face%reg        = grid%face%region%Get()
        simgrid%face%ntot       = nf 

        ! Cells
        simgrid%cell%vert           = grid%cell%vert%Get()
        simgrid%cell%vertP(:, 1)    = grid%cell%vp1%Get()
        simgrid%cell%vertP(:, 2)    = grid%cell%vp2%Get()
        simgrid%cell%ntot           = nc 

        associate(&
            c     => simgrid%cell, &
            v   => simgrid%vert)
        do i = 1, nc
            ! Get cell vertices
            tv = GetCellVert(c, i)

            ! Compute coordinates
            c%x(i) = sum(v%x(tv))/real(size(tv), kind=R8)
            c%y(i) = sum(v%y(tv))/real(size(tv), kind=R8)
        end do
        end associate

        ! Magnetic field data
        !====================
        ! Vertices
        call magneticField%interp%Evaluate(simgrid%vert%x, simgrid%vert%y, &
            0, 0, simgrid%vert%psi)
        call magneticField%interp%Evaluate(simgrid%vert%x, simgrid%vert%y, &
            1, 0, simgrid%vert%bx)
        call magneticField%interp%Evaluate(simgrid%vert%x, simgrid%vert%y, &
            0, 1, simgrid%vert%by)
        simgrid%vert%ffbz = magneticField%RBtor*2.0_R8*pi_R8

        ! Faces
        associate(&
            fv      => simgrid%face%vert, &
            vfID    => simgrid%vert%fieldlineID)
        simgrid%face%aligned = 0_I8 
        do i = 1, simgrid%face%ntot 
            if ((vfID(fv(i, 1)) /= 0) .and. (vfID(fv(i, 2)) /= 0)) then
                if (vfID(fv(i, 1)) == vfID(fv(i, 2))) then  
                    simgrid%face%aligned(i) = 1_I8
                end if 
            end if 
        end do 
        end associate 

        ! Cells
        associate(&
            c     => simgrid%cell, &
            v     => simgrid%vert)
        call magneticField%interp%Evaluate(c%x, c%y, 0, 0, c%psi)
        bpvx = sqrt(v%bx**2 + v%by**2)
        do i = 1, nc
            ! Get cell vertices
            tv = GetCellVert(c, i)

            ! Compute poloidal and toroidal field
            c%bp(i) = -sum(bpvx(tv))/(real(size(tv), kind=R8)*2.0_R8*c%x(i)*pi_R8)
            c%bt(i) = sum(v%ffbz(tv))/(real(size(tv), kind=R8)*2.0_R8*c%x(i)*pi_R8)
        end do 
        end associate

        ! Grid interconnection
        !=====================
        ! Classic grid interconnection data
        call ComputeGridInterconnections(simgrid)

        ! Additional grid data
        call ComputeGridData(simgrid, magneticField)

        ! Grid boundary
        !call ComputeGridBoundaries(simgrid)


        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Compute grid data
    subroutine ComputeGridData(simgrid, magneticField)

        ! Description
        !============
        ! Compute additional grid data related to flux surfaces, 
        ! flux tubes, etc. This does not include the grid boundaries
        ! yet -> separate routine. 

        ! Note 1: in principle, this routine only requires simulation grid
        ! data, so perhaps this can be moved to the goatmod_types module?

        ! Note 2: we don't account yet for empty flux surfaces - these 
        ! simply have zero faces/cells etc. Can be deleted later on

        ! Note 3: cut cells etc should be correctly included as long as 
        ! cells only have two boundary faces - otherwise they're not included

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(inout)        :: simgrid 
        type(MagneticFieldUDT), intent(in)  :: magneticField

        ! Auxiliary
        integer(I8)                             :: nfsf, nfsftot, &
            nftftot, nftc, si, ei, np
        integer(I8), allocatable, dimension(:)  :: ffsID, fID, vID, &
            sortind, polind, tfind, ftf, ftc, tcf, IDs, allcellreg, ind, &
            nfsneig, tcv
        integer(I8), allocatable, dimension(:, :)   :: tfnb
        real(R8), allocatable, dimension(:)     :: vpsi, xf, yf, xc, yc, &
            ccx, ccy, bfx, bfy, dp
        logical, allocatable, dimension(:)      :: temp, tf, &
            ispolygonstart, tc

        ! Loop
        integer(I8)                             :: i, j, k, cc, ncell, &
            nface

        ! Initialize
        !===========
        associate(&
            fd      => simgrid%data%fluxdata,   &
            c       => simgrid%cell,            &
            f       => simgrid%face,            &
            v       => simgrid%vert)

        ! Flux surfaces
        !==============
        ! Set initial number of flux surfaces 
        fd%nFs      = maxval(v%fieldlineID)
        
        ! Initialize other fields
        if (allocated(fd%fluxsurfacefacesP)) then 
            deallocate(fd%fluxsurfacefacesP)
        end if 
        allocate(fd%fluxsurfacefacesP(fd%nFs, 2))
        fd%fluxsurfacefacesP = 0_I8
        
        ! Determine flux surface ID of faces
        allocate(ffsID(f%ntot))
        ffsID = 0 
        do i = 1, f%ntot
            ! Check
            if (v%fieldlineID(f%vert(i, 1)) == v%fieldlineID(f%vert(i, 2))) then 
                ffsID(i) = v%fieldlineID(f%vert(i, 1))
            end if 
        end do 

        ! Extract
        nfsftot = count(ffsID /= 0)
        if (allocated(fd%fluxsurfacefaces)) then 
            deallocate(fd%fluxsurfacefaces)
        end if 
        if (allocated(fd%fluxsurfacepsi)) then 
            deallocate(fd%fluxsurfacepsi)
        end if 
        allocate(fd%fluxsurfacefaces(nfsftot), fd%fluxsurfacepsi(fd%nFs))
        fd%fluxsurfacefaces = 0
        fd%fluxsurfacefacesP(1, 1) = 1 
        fID = [(k, k = 1, f%ntot)]
        vID = [(k, k = 1, v%ntot)]
        do i = 1, fd%nFs 
            ! Get faces with this ID
            temp = ffsID == i
            nfsf = count(temp)

            ! Set pointer
            fd%fluxsurfacefacesP(i, 2) = nfsf 
            if (i < fd%nFs) then 
                fd%fluxsurfacefacesP(i+1, 1) = fd%fluxsurfacefacesP(i, 1) &
                    + nfsf
            end if 

            ! Set faces
            fd%fluxsurfacefaces(fd%fluxsurfacefacesP(i, 1):fd%fluxsurfacefacesP(i, 1)+nfsf-1) & 
                = pack(fID, temp)

            ! Set ID
            fd%fluxsurfaceID(i) = i

            ! Compute psi value as mean of vertex psi value (assumed computed)
            allocate(vpsi(count(v%fieldlineID == i)))
            vpsi = pack(v%psi, v%fieldlineID == i)
            fd%fluxsurfacepsi(i) = sum(vpsi)/real(size(vpsi), kind=R8)
            deallocate(vpsi)

        end do 

        ! Flux tubes
        !===========
        ! Get all non-aligned, non-boundary faces
        tf = (.not. f%BF) .and. (f%cellP(:, 2) == 2) .and. (.not. (f%aligned == 1))

        ! Get all single cells with two boundary faces that have the 
        ! same non-zero flux surface ID 
        allocate(tc(c%ntot))
        tc = .false. 
        do i = 1, c%ntot
            ! Get cell faces
            tcf = GetCellFace(c, i)
            tcf = pack(tcf, f%BF(tcf) .and. (.not. (f%aligned(tcf) == 1)))

            ! Get unique vertex ID without zeros
            if (size(tcf) == 2) then 
                call Unique([v%fieldlineID(f%vert(tcf, 1)), &
                    v%fieldlineID(f%vert(tcf, 2))], IDs)
                IDs = pack(IDs, IDs /= 0)
                if (size(IDs) == 2) then 
                    tc(i) = .true.
                end if
            end if
        end do  

        ! Construct the face neighbour
        allocate(tfnb(count(tf), 2))
        tfnb = 0
        cc = 0
        do i = 1, f%ntot 
            if (tf(i)) then 
                cc = cc + 1
                tfnb(cc, 1:f%cellP(i, 2)) = GetFaceCell(f, i)
            end if 
        end do 

        ! Sort the edges (assumed no branching)
        allocate(sortind(size(tfnb, 1)), ispolygonstart(size(tfnb, 1)))
        call SortPolygonEdges(tfnb, count(tf), sortind, ispolygonstart)
        allocate(polind(count(ispolygonstart)))
        polind = pack([(k, k = 1, size(tfnb, 1))], ispolygonstart)
        polind = [polind, size(tfnb, 1)+1]

        ! Sort faces
        allocate(tfind(count(tf)))
        tfind = pack(fID, tf)
        tfind = tfind(sortind)

        ! Extract tubes
        nftftot = 2*count(tf)  + 2*count(tc) ! overestimation
        np = count(ispolygonstart)
        fd%nFt = np + count(tc)
        if (allocated(fd%fluxtubefacesP)) then 
            deallocate(fd%fluxtubefacesP)
        end if 
        if (allocated(fd%fluxtubefaces)) then 
            deallocate(fd%fluxtubefaces)
        end if 
        if (allocated(fd%fluxtubefsIDs)) then 
            deallocate(fd%fluxtubefsIDs)
        end if 
        if (allocated(fd%fluxtuberegID)) then 
            deallocate(fd%fluxtuberegID)
        end if 
        if (allocated(fd%fluxtubecellsP)) then 
            deallocate(fd%fluxtubecellsP)
        end if 
        if (allocated(fd%fluxtubecells)) then 
            deallocate(fd%fluxtubecells)
        end if 
        if (allocated(fd%isclosedft)) then 
            deallocate(fd%isclosedft)
        end if 
        allocate(fd%fluxtubefacesP(fd%nFt, 2), fd%fluxtubefaces(nftftot), &
            fd%fluxtubefsIDs(fd%nFt, 2), fd%fluxtuberegID(fd%nFt), &
            fd%fluxtubecellsP(fd%nFt, 2), fd%fluxtubecells(count(tf)*2+count(tc)), &
            fd%isclosedft(fd%nFt))
        fd%fluxtubefacesP = 0
        if (fd%nFs > 0) then 
            fd%fluxtubefacesP(1, 1) = 1
        end if 
        if (fd%nFt > 0) then 
            fd%fluxtubecellsP(1, 1) = 1
        end if 
        fd%fluxtubefaces = 0
        fd%fluxtubefsIDs = 0
        fd%fluxtuberegID = 0
        fd%isclosedft = .false.

        ncell = 0
        nface = 0
        do i = 1, np
            ! Get index
            si = polind(i)
            ei = polind(i+1)-1

            ! Get cells
            nftc = ei - si + 2
            allocate(ftc(nftc))
            call ExtractPolygonVertices(tfnb(si:ei, :), nftc-1, ftc)

            ! Check for closed tubes
            if (ftc(1) /= ftc(nftc)) then 
                ! Add cells
                fd%fluxtubecells(ncell+1:ncell+nftc) = ftc 
                fd%fluxtubecellsP(i, 2) = nftc 
                if (i < fd%nFt) then 
                    fd%fluxtubecellsP(i+1, 1) = fd%fluxtubecellsP(i, 1) + nftc 
                end if 
                c%ft(ftc) = i 

                ! Update counters
                ncell = ncell + nftc

                ! Get inner faces
                ftf = tfind(si:ei)

                ! Find boundary face(s) of first cell
                tcf = GetCellFace(c, ftc(1))

                ! Remove non-boundary and aligned faces
                tcf = pack(tcf, f%BF(tcf) .and. (.not. (f%aligned(tcf) == 1)))

                ! Check
                if (size(tcf) == 1) then
                    ! Standard case, add face 
                    ftf = [tcf, ftf]
                elseif (size(tcf) < 1) then 
                    ! Unexpected, but not an issue
                    print *, 'ComputeGridData: flux tube: ', i, &
                        ' is not closed but the starting cell (number: ', ftc(1), &
                        ' ) has no boundary faces. Not adding face'
                elseif (size(tcf) > 1) then 
                    ! Unexpected, may be an issue
                    print *, 'ComputeGridData: flux tube: ', i, &
                        ' is not closed but the starting cell (number: ', ftc(1), &
                        ' ) has multiple boundary faces. Not adding face'
                end if 

                ! Find boundary face(s) of last cell
                tcf = GetCellFace(c, ftc(nftc))

                ! Remove non-boundary faces
                tcf = pack(tcf, f%BF(tcf) .and. (.not. (f%aligned(tcf) == 1)))

                ! Check
                if (size(tcf) == 1) then
                    ! Standard case, add face 
                    ftf = [ftf, tcf]
                elseif (size(tcf) < 1) then 
                    ! Unexpected, but not an issue
                    print *, 'ComputeGridData: flux tube: ', i, &
                        ' is not closed but the ending cell (number: ', ftc(1), &
                        ' ) has no boundary faces. Not adding face'
                elseif (size(tcf) > 1) then 
                    ! Unexpected, may be an issue
                    print *, 'ComputeGridData: flux tube: ', i, &
                        ' is not closed but the ending cell (number: ', ftc(1), &
                        ' ) has multiple boundary faces. Not adding face'
                end if 
            else 
                ! Closed tube
                fd%isclosedft(i) = .true. 

                ! Add cells (don't take last one, is duplicate)
                fd%fluxtubecells(ncell+1:ncell+nftc-1) = ftc(1:nftc-1)
                fd%fluxtubecellsP(i, 2) = nftc-1
                if (i < fd%nFt) then 
                    fd%fluxtubecellsP(i+1, 1) = fd%fluxtubecellsP(i, 1) + nftc-1 
                end if 
                c%ft(ftc(1:nftc-1)) = i

                ! Update counters
                ncell = ncell + nftc - 1

                ! Get faces
                ftf = tfind(si:ei)

            end if 

            ! Get vertex IDs (without zero IDs)
            call Unique(v%fieldlineID([f%vert(ftf, 1), f%vert(ftf, 2)]), IDs)

            ! Remove zero IDs
            IDs = pack(IDs, IDs /= 0_I8)

            ! Sanity check
            if (size(IDs) /= 2) then 
                ! This should actually not happen
                print *, 'ComputeGridData: flux tube ', i, ' has ', size(IDs), &
                    ' non-zero IDs, unexpected. Setting flux surface IDs to zero'
                fd%fluxtubefsIDs(i, :) = 0
            else 
                ! Add
                fd%fluxtubefsIDs(i, :) = IDs 
            end if 

            ! Get flux tube region
            allcellreg = c%reg(ftc)

            ! Check
            if (any(allcellreg(1) /= allcellreg)) then 
                ! Print warning
                print *, 'ComputeGridData: multiple cell regions detected ' // & 
                    'for flux tube: ', i, '. Taking first cell region...'

            end if 
            fd%fluxtuberegID(i) = allcellreg(1)

            ! Add faces
            fd%fluxtubefaces(nface+1:nface+size(ftf)) = ftf 
            fd%fluxtubefacesP(i, 2) = size(ftf)
            if (i < fd%nFt) then 
                fd%fluxtubefacesP(i+1, 1) = fd%fluxtubefacesP(i, 1) + size(ftf)
            end if 

            ! Update counters
            nface = nface + size(ftf)

            ! Housekeeping
            deallocate(ftc)

        end do 

        ! Add single cell tubes
        cc = np
        do i = 1, c%ntot
            if (tc(i)) then 
                ! Update counter
                cc = cc + 1

                ! Add cell
                fd%fluxtubecells(ncell+1:ncell+1) = i
                fd%fluxtubecellsP(cc, 2) = 1
                if (cc < fd%nFt) then 
                    fd%fluxtubecellsP(cc+1, 1) = fd%fluxtubecellsP(cc, 1) + 1
                end if 

                ! Get faces
                tcf = GetCellFace(c, i)
                tcf = pack(tcf, f%BF(tcf))
                
                ! Add
                fd%fluxtubefaces(nface+1:nface+2) = tcf 
                fd%fluxtubefacesP(cc, 2) = 2
                if (cc < fd%nFt) then 
                    fd%fluxtubefacesP(cc+1, 1) = fd%fluxtubefacesP(cc, 1) + 2
                end if 

                ! Get flux surface IDs
                call Unique([v%fieldlineID(f%vert(tcf, 1)), &
                    v%fieldlineID(f%vert(tcf, 2))], IDs)
                IDs = pack(IDs, IDs /= 0)  
                fd%fluxtubefsIDS(cc, :) = IDs 

                ! Set cell ft region
                c%ft(i) = cc

                ! Set fluxtube region ID
                fd%fluxtuberegID(cc) = c%reg(i)

                ! Update counters
                nface = nface + 2
                ncell = ncell + 1

            end if  
        end do 

        ! Trim
        fd%fluxtubefaces = fd%fluxtubefaces(1:nface)
        fd%fluxtubecells = fd%fluxtubecells(1:ncell)

        ! Additional data
        !================
        ! Flux tube neighbours & tangency points
        if (allocated(fd%fluxsurfaceneig)) then 
            deallocate(fd%fluxsurfaceneig)
        end if 
        if (allocated(fd%fluxsurfaceneigP)) then 
            deallocate(fd%fluxsurfaceneigP)
        end if 
        allocate(fd%fluxsurfaceneigP(fd%nFs, 2))
        fd%fluxsurfaceneigP = 0

        ! Compute pointer
        fd%fluxsurfaceneigP(1, 1 ) = 1
        do i = 1, fd%nFt
            ! Get flux tube faces
            tcf = GetFTFace(fd, i)

            ! Get unique IDs
            call Unique([v%fieldlineID(f%vert(tcf, 1)), &
                v%fieldlineID(f%vert(tcf, 2))], IDs)
            IDs = pack(IDs, IDs /= 0)

            ! Update counters
            fd%fluxsurfaceneigP(IDs, 2) = fd%fluxsurfaceneigP(IDs, 2) + 1
        end do 
        do i = 2, fd%nFs
            fd%fluxsurfaceneigP(i, 1) = fd%fluxsurfaceneigP(i-1, 1) + &
                fd%fluxsurfaceneigP(i-1, 2)
        end do 

        ! Compute neighbours
        allocate(fd%fluxsurfaceneig(sum(fd%fluxsurfaceneigP(:, 2))))
        allocate(nfsneig(fd%nFs))
        nfsneig = 0
        do i = 1, fd%nFt 
            ! Get flux tube faces
            tcf = GetFTFace(fd, i)

            ! Get unique IDs
            call Unique([v%fieldlineID(f%vert(tcf, 1)), &
                v%fieldlineID(f%vert(tcf, 2))], IDs)
            IDs = pack(IDs, IDs /= 0)

            ! Add
            ind = fd%fluxsurfaceneigP(IDs, 1) + nfsneig(IDs)
            fd%fluxsurfaceneig(ind) = IDs 

            ! Update
            nfsneig(IDs) = nfsneig(IDs) + 1
        end do 

        ! Check orientation
        !==================
        ! Note: previous algorithm should already result in sorted flux
        ! tube cells and faces - only need to check orientation w.r.t. 
        ! magnetic field
        do i = 1, fd%nFt 
            ! Get faces, cells
            ftf = GetFTFace(fd, i)
            ftc = GetFTCell(fd, i)

            ! Compute face and cell center coordinates
            allocate(xf(size(ftf)), yf(size(ftf)), xc(size(ftc)), yc(size(ftc)))
            do j = 1, size(ftf)
                xf(j) = 0.5*sum(v%x(f%vert(ftf(j), :)))
                yf(j) = 0.5*sum(v%y(f%vert(ftf(j), :)))
            end do 
            do j = 1, size(ftc)
                tcv = GetCellVert(c, ftc(j))
                xc(j) = sum(v%x(tcv))/(real(size(tcv), kind=R8))
                yc(j) = sum(v%y(tcv))/(real(size(tcv), kind=R8))
            end do 

            ! Compute magnetic field vector at face centers
            allocate(bfx(size(xf)), bfy(size(xf)))
            call magneticField%interp%Evaluate(xf, yf, 1, 0, bfy)
            call magneticField%interp%Evaluate(xf, yf, 0, 1, bfx)
            bfx = -bfx 

            ! Compute cell connector
            ccx = xc(2:size(xc)) - xc(1:size(xc)-1)
            ccy = yc(2:size(yc)) - yc(1:size(yc)-1)

            ! Compute dot product between inner faces and cell connector
            if (fd%isclosedft(i)) then 
                ! Close the tube
                ccx = [ccx, ccx(1)]
                ccy = [ccy, ccy(1)]
                dp = (ccx*bfx + ccy*bfy)
            else
                dp = (ccx*bfx(2:size(ftf)-1) + ccy*bfy(2:size(ftf)-1))
            end if 

            ! Switch if necessary
            if (sum(dp) < 0) then 
                ftf = ftf(size(ftf):1:-1)
                ftc = ftc(size(ftc):1:-1)
                fd%fluxtubefaces(fd%fluxtubefacesP(i, 1):&
                    fd%fluxtubefacesP(i, 1)+fd%fluxtubefacesP(i, 2)-1) = ftf 
                fd%fluxtubecells(fd%fluxtubecellsP(i, 1):&
                    fd%fluxtubecellsP(i, 1)+fd%fluxtubecellsP(i, 2)-1) = ftc
            end if  

            ! Housekeeping
            deallocate(xf, yf, xc, yc, bfx, bfy)

        end do

        print *, fd%nFt

        ! Housekeeping
        !=============
        end associate


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
        type(GGTMFieldlineDataUDT)             :: line 
        class(TopomeshUDT), intent(in)          :: topomesh 
        class(GGTMDataUDT)                      :: ggtmdata 

        ! Auxiliary
        real(R8)                                :: nxsrf, nysrf, nxl, &
            nyl, thisfl
        real(R8), allocatable, dimension(:)     :: tx, ty, xl, yl, dlcv, &
            tdlcv
        integer(I8)                             :: startv, endv, &
            thisf, thisfind, nextv, fsID
        integer(I8), allocatable, dimension(:)  :: bndf, bndv, tvID, &
            allfsID
        logical                                 :: doflip
        logical, allocatable, dimension(:)      :: isnotfound 

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
                srfvx   => facedata(tmcell%srf)%line%xv,     &
                srfvy   => facedata(tmcell%srf)%line%yv)

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
                srfvx   => facedata(tmcell%srf)%line%xv,     &
                srfvy   => facedata(tmcell%srf)%line%yv)

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
                call line%Initialize(vert%x(bndv), vert%y(bndv), vert%fsID(bndv(1)))
                call line%AddVertexCoordinates([0.0_R8])
                call line%AddVertexIDs(bndv)
                call line%UpdateLineData(topomesh, ggtmdata)

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
        thisfl  = facedata(thisf)%line%dllc(facedata(thisf)%line%nl)
        if (doflip) then 
            ! Reconstruct by flipping line length distribution
            xl = facedata(thisf)%line%xl(facedata(thisf)%line%nl:1:-1)
            yl = facedata(thisf)%line%yl(facedata(thisf)%line%nl:1:-1)
            dlcv = facedata(thisf)%line%dllc(facedata(thisf)%line%nl) - &
                facedata(thisf)%line%dlcv(facedata(thisf)%line%nv:1:-1)
            tvID = facedata(thisf)%line%vert(facedata(thisf)%line%nv:1:-1)
        else
            xl = facedata(thisf)%line%xl 
            yl = facedata(thisf)%line%yl 
            tvID = facedata(thisf)%line%vert 
            dlcv = facedata(thisf)%line%dlcv
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
                tdlcv = facedata(bndf(i))%line%dlcv 
                tdlcv = facedata(bndf(i))%line%dllc(facedata(bndf(i))%line%nl) - tdlcv
                tdlcv = tdlcv + thisfl! Need to add length!
                xl = [xl(1:size(xl)-1), tx(size(tx):1:-1)]
                yl = [yl(1:size(yl)-1), ty(size(ty):1:-1)]
                tvID = [tvID(1:size(tvID)-1), facedata(bndf(i))%line%vert(facedata(bndf(i))%line%nv:1:-1)]
                dlcv = [dlcv(1:size(dlcv)-1), tdlcv(size(tdlcv):1:-1)]
            
            else
                tx = face%x(bndf(i))%Get()
                ty = face%y(bndf(i))%Get()
                tdlcv = facedata(bndf(i))%line%dlcv 
                tdlcv = tdlcv + thisfl! Need to add length!
                xl = [xl(1:size(xl)-1), tx]
                yl = [yl(1:size(yl)-1), ty]
                tvID = [tvID(1:size(tvID)-1), facedata(bndf(i))%line%vert]
                dlcv = [dlcv(1:size(dlcv)-1), tdlcv]
            end if 

            ! Update face length
            thisfl  = thisfl + facedata(bndf(i))%line%dllc(facedata(bndf(i))%line%nl)
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
            nxl = xl(2) - xl(1)
            nyl = yl(2) - yl(1)

            ! Compute dot product and check
            if ((nxl*nxsrf + nyl*nysrf) < 0) then 
                xl = xl(size(xl):1:-1)
                yl = yl(size(yl):1:-1)
                thisfl = sum(sqrt((xl(2:) - xl(1:size(xl)-1))**2 + &
                    (yl(2:) - yl(1:size(yl)-1))**2))
                dlcv = thisfl - dlcv(size(dlcv):1:-1)
                tvID = tvID(size(tvID):1:-1)
            end if 
            
        end if 

        ! Check which flux surface ID we should take, and issue warning
        ! if these are not the same
        allfsID = topomesh%face%fsID(bndf)
        if (any(allfsID(1) /= allfsID)) then 
            print *, 'ExtractTMCellAlignedBoundary: not all faces have the ' // & 
                'same flux surface ID (IDs: ', allfsID, ', taking first one'
        end if 
        fsID = allfsID(1)

        ! Construct line
        call line%Initialize(xl, yl, fsID)
        dlcv(1) = 0 ! ensure start point lies on start
        dlcv(size(dlcv)) = line%dllc(line%nl) ! ensure end point lies on line end
        call line%AddVertexCoordinates(dlcv)
        call line%AddVertexIDs(tvID)
        call line%UpdateLineData(topomesh, ggtmdata)

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Determine if vertex ID is a topomesh vertex ID
    function IsTopomeshVert(vertID, topomesh) result(out)

        ! Description
        !============
        ! Determine if the vertex IDs specified in vertID are topomesh
        ! vertices. This is true if 0 < vertID <= topomesh%vert%ntot

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)                 :: vertID(:)
        type(TopomeshUDT), intent(in)           :: topomesh 
        logical, allocatable                    :: out(:)

        ! Determine
        !==========
        out = (vertID <= topomesh%vert%ntot) .and. &
            vertID > 0

    end function 

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
            if (allocated(f(i)%line%xv) .and. allocated(f(i)%line%yv) .and. &
                allocated(f(i)%line%vert)) then 
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
                write (fu, *) i, size(f(i)%line%vert)
            end if 
        end do 
        
        ! Face coordinates
        write (fu, *) 'face coordinates'
        do i = 1, size(f)
            if (doface(i)) then 
                ! Write header
                write (fu, *) 'face ', i 

                ! Write coordinates
                do j = 1, size(f(i)%line%xv)
                    write(fu, *) f(i)%line%vert(j), f(i)%line%xv(j), f(i)%line%yv(j)
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

    ! Grid data writing
    subroutine WriteGGGridData(grid, filename)

        ! Description
        !============
        ! Write out grid data of intermediate grid structure. Mainly 
        ! for plotting purposes, quite limited in output. 

        ! 'vertices'
        ! <vert%ntot>
        ! 'ID, x, y, fieldlineID'
        ! <ID, x, y, fieldlineID>
        ! 'faces'
        ! <face%ntot> 
        ! 'ID, v1, v2, label, region'
        ! <ID, v1, v2, label, region>
        ! 'cells'
        ! <cell%ntot, cell%nvert> 
        ! 'ID, vp1, vp2, region>'
        ! <ID, vp1, vp2, region>
        ! 'cell vertices'
        ! <cell%vert> 

        ! Declare variables
        !==================
        ! Modules 
        use mod_plotter 
        use mod_specialchars, only : filesepchar

        ! Arguments
        class(GGGridUDT)                        :: grid
        character(*), intent(in)                :: filename 

        ! Auxiliary
        integer                                 :: fu
        real(R8), allocatable, dimension(:)     :: x, y 
        integer(I8), allocatable, dimension(:)  :: fID, v1, v2, region, &   
            label, vc
        character(:), allocatable               :: dir

        ! Loop
        integer(I8)                             :: i

        ! Initialize
        !===========
        ! Unpack
        associate(&
            f       => grid%face,   &
            c       => grid%cell,   &
            v       => grid%vert    &
        )

        ! Construct writing directory
        dir = plotdir // filesepchar // filename // '.dat'

        ! Open file
        open (action='write', file=trim(dir), newunit=fu, &
             status='unknown')

        ! Write header
        write(fu, *) 'VERSION3.00.00'

        ! Write vertex data
        !==================
        ! Unpack
        x = v%x%Get()
        y = v%y%Get()
        fID = v%fieldlineID%Get()

        ! Number of vertices
        write (fu, *) 'vertices'
        write (fu, *) v%ntot 

        ! Vertex data
        write (fu, *) 'ID, x, y, fieldlineID'
        do i = 1, v%ntot 
            write (fu, *) i, x(i), y(i), fID(i)
        end do 

        ! Write face data
        !================
        ! Unpack
        v1 = f%v1%Get()
        v2 = f%v2%Get()
        region = f%region%Get()
        label = f%label%Get()

        ! Number of faces
        write (fu, *) 'faces'
        write (fu, *) f%ntot

        ! Face data
        write (fu, *) 'ID, v1, v2, label, region'
        do i = 1, f%ntot
            write (fu, *) i, v1(i), v2(i), label(i), region(i)
        end do 

        ! Write cell data
        !================
        ! Unpack
        vc = c%vert%Get()
        v1 = c%vp1%Get()
        v2 = c%vp2%Get()
        region = c%region%Get()

        ! Number of cells
        write (fu, *) 'cells'
        write (fu, *) c%ntot, size(vc)

        ! Cell data
        write (fu, *) 'ID, vp1, vp2, region'
        do i = 1, c%ntot
            write (fu, *) i, v1(i), v2(i), region(i)
        end do 

        ! Cell vertices
        write (fu, *) 'cell vertices'
        do i = 1, size(vc)
            write (fu, *) vc(i)
        end do 

        ! Housekeeping
        !=============
        ! Deallocate again

        ! Others
        end associate
        close(fu)

    end subroutine


end module 