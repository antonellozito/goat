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
! - To deal with shitty small features/cell overlap/very large regions
!   with type 1 tangency points (huge amount of stacked triangles 
!   typically), etc etc, we first post-process the initial vertex 
!   distribution, where each vertex lies on a flux surface/line. 
!   Depending on user input, we then construct for each cell, and for
!   each line pair, a flux tube that will be used for gridding later on.
!   The two lines of this pair may be further extended to allow e.g.
!   cut cells, or to insert vertices along an almost aligned vessel 
!   boundary for regions with type 1 tangency points. After this step, 
!   grid vertices will emerge that do not lie on flux surfaces (i.e. the
!   cut cell vertices typically). 
! 
! Note that we assume that all cells are present in flux tubes. 
! Below we mention the different algorithms for distributing the 
! vertices along the field lines. Note that the amount of field lines
! and the amount of vertices along a field line can be influenced 
! locally by setting the properties of the poloidal and radial 
! vertex distributors (see also ggmod_Vertexdistribution2D), as well as 
! the line refiners defined here (this typically leads to less smooth 
! cell transitions, but is often superior to control local refinement, 
! impose boundary layers, etc)
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
    use mod_definitions
    use mod_linearsolverinterface, only: SolveDenseLinearSystemDI
    use goatmod_types, only : magneticFieldUDT, VesselUDT, GridUDT
    use goatmod_userinput, only : GGoptionsUDT
    use ggmod_topology2D
    use ggmod_vertexdistribution2D
    use DistributionFunction
    use mod_plotter
    use mod_utility, only: wall_time
    use omp_lib
    use mod_graph
    implicit none
    private 
    public :: GenerateUnstructuredAlignedGrid, TranslateGridLabels, &
        ComputeTopologicalData, GetGridFaceLabelMappingGD

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

    ! Field line segment data
    type :: GGTMSegmentUDT

        ! Description
        !============
        ! This data type forms the backbone for grid generation as it 
        ! holds the original coordinates of a line segment (can be any
        ! line, but it has to be a simple polygon) and later on the 
        ! coordinates of the vertices and their IDs. Line coordinates are
        ! stored in xl, yl, vertex coordinates in xv, yv. fsID indicates 
        ! the flux surface ID that it corresponds to (zero if none), 
        ! TMfaceID corresponds to the topomesh face ID (zero if none).
        ! sv and ev are the start and end vertices of the grid of the 
        ! line (these should never be removed!) and go from start to 
        ! end of the segment. Segments should be stored in the general
        ! GGTMData structure

        ! Note: a segment can also represent a single vertex (in this 
        ! case 'isvertex' will be true), for which xl, yl still have 
        ! two entries, as do vert etc, but the TMfaceID will correspond
        ! (normally) to the topological mesh vertex index. 

        ! Note: segments can be closed, but they have very limited 
        ! support. Errors will be thrown.

        real(R8), allocatable, dimension(:)         :: xv, yv, xl, yl, &
            dll, dllc, dlcv
        integer(I8), allocatable, dimension(:)      :: vert
        integer(I8)                                 :: nv, nl, sv, ev, &
            fsID, TMfaceID
        logical                                     :: isvertex, isclosed

    contains 

        ! Initialization
        procedure :: Initialize     => InitializeGGTMSegment

        ! Vertex coordinates addition
        procedure :: AddVertices    => AddVerticesGGTMSegment

        ! Cleaning (removal of small edges)
        procedure :: Clean          => CleanGGTMSegment

        ! Segment flipping 
        procedure :: Flip           => FlipGGTMSegment

    end type


    ! Field line data
    type :: GGTMFieldlineDataUDT

        ! Description
        !============
        ! This data type holds information on the field line, such as
        ! the original coordinates xl, yl, the vertex coordinates xv, yv
        ! (and the one-dimensional coordinates dlcv that go from 0 to the
        ! number of edges), the vertex IDs ('vert') and the face labels
        ! of faces formed by vertex pairs formed by subsequent vertices. 
        ! 
        ! This type is now also extended to deal with multiple segments
        ! (though we still assume that it is a simple polygon) that each
        ! can have a flux surface ID (zero if it is not a flux surface), 
        ! and so on. The vertex distribution and initial coordinates are
        ! still defined over the entire field line, the parts can be 
        ! identified by the 'nodes', which contain the index in xl, yl 
        ! coordinates of the different parts. To track facelabels, we 
        ! also store the topological mesh face ID (if the line is based
        ! on that, this is zero if the line is not based on any 
        ! topological mesh face)

        ! For ease of implementation, one can extend the original line
        ! by appending (either at start or end of the original line) 
        ! additional lines. These should, of course, start or end in the
        ! same point. 

        ! It is important to notice that this extension allows things 
        ! such as the cut-cell approach, but that the main idea of a 
        ! field aligned approach (with aligned and non-aligned 
        ! boundaries) remains! 

        real(R8), allocatable, dimension(:)         :: xv, yv, xl, yl, &
            dll, dllc, dlcv
        integer(I8), allocatable, dimension(:)      :: vert, facelabels, &
            fsID, nodes, TMfaceID, nodevert, segID
        logical, allocatable, dimension(:)          :: flipseg, isnodevert
        integer(I8)                                 :: nv, nl, ns 

    contains 

        ! Initialization
        procedure :: Initialize     => InitializeGGTMFieldLineData

        ! Appending
        procedure :: AppendSegment  => AppendGGTMFieldLineSegment

        ! Flipping
        procedure :: Flip           => FlipGGTMFieldLine

        ! Splitting into segments at vertex
        procedure :: SplitAtVertex  => SplitGGTMFieldLineAtVertex
        procedure :: SplitAtVertices => SplitGGTMFieldLineAtVertices
        procedure :: SplitAtNodes   => SplitGGTMFieldLineAtSegmentNodes

        ! Vertex addition
        procedure :: AddVertexCoordinates 

        ! Vertex ID addition
        procedure :: AddVertexIDs

        ! Data addition
        procedure :: UpdateLineData
        procedure :: UpdateSegmentData
        procedure :: UpdateLineGriddingData

        ! Getters
        procedure :: GetSegmentFaceIndices  => GetGGTMFieldLineSegmentFaceIndices
        procedure :: GetSegmentVertIndices  => GetGGTMFieldLineSegmentVertIndices
        procedure :: GetAllSegmentVertIndices   => GetGGTMFieldLineAllSegmentVertIndices

    end type

    ! Field line pair data
    type :: GGTMFieldlinePairDataUDT

        ! Data type that contains data on pairs of field lines, such as 
        ! line-of-sight etc. Also holds (possibly extended) field line
        ! line data for the high field and low field line of the tube.
        ! These lines of course don't have to be at the high or low 
        ! field.  
        type(GGTMFieldlineDataUDT)              :: hfline, lfline 
        real(R8), allocatable, dimension(:)     :: graphxv, graphyv
        integer(I8)                             :: srflabel, erflabel
        integer(I8), allocatable, dimension(:)  :: l1minLOS, l1maxLOS, &
            l2minLOS, l2maxLOS
        logical                                 :: isextendedstart, &
            isextendedend, dograph
        type(IntegerDynamicArrayUDT)    :: hfface, lfface, tubeface, &
            cell 
        type(UGraphUDT)             :: graph

    contains 

        ! Initialization
        procedure :: Initialize         => InitializeGGTMFieldlinePairData

        ! Graph initialization
        procedure :: InitializeGraph    => InitializeGGTMFieldlinePairDataGraph  

        ! Graph visualization 
        procedure :: VisualizeGraph     => VisualizeGGTMFieldlinePairDataGraph

    end type

    ! Field line refinement options
    type :: GGTMFieldlineRefinementOptionsUDT 

        ! Description
        !============
        ! This structure contains all refinement options for refining
        ! a GGTM line with any GGTM line refiner. These options should
        ! be added as a structure to the GGTMCellData structure, such 
        ! that this data transfer can be done smoothly

        ! Data for length-based refinement
        logical                                 :: doBLstart, doBLend, &
            dlBLlengthbased 
        integer(I8)                             :: ncBLstart, ncBLend 
        real(R8), allocatable, dimension(:)     :: dlBLstart, dlBLend

    contains 

        ! Initializer to avoid issues
        procedure :: Initialize         => InitializeGGTMLineRefinementOptions

    end type

    ! Vertex data
    type :: GGTMVertexDataUDT

        ! Description
        !============
        ! This data type holds additional vertex information for the 
        ! grid generation. Right now, it has a 'line' field, which is
        ! only used to represent the vertex as a field 'line' for ease
        ! during the grid generation (to deal with e.g. lines formed
        ! only by tangency points)
        type(GGTMFieldlineDataUDT)              :: line 

    end type 
    
    ! Face data
    type :: GGTMFaceDataUDT

        ! Description
        !============
        ! Contains additional face data, including the grid vertices
        ! xv, yv - now stored as a GGTMFieldlineData type. Also contains
        ! line refinement options (similar to celldata%linerefoptions)
        type(GGTMFieldlineDataUDT)              :: line 
        type(GGTMFieldlineRefinementOptionsUDT) :: linerefoptions

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
        ! are sorted from start to end radial face. Additionally, we 
        ! include any user-defined refinement/mesh size/... options 
        ! here. Also flux tube data is stored for grid
        ! face/cell construction etc after an initial vertex distribution
        ! has been constructed
        type(GGTMFieldlineDataUDT), allocatable     :: lines(:)
        type(GGTMFieldlinePairDataUDT), allocatable :: tubes(:)
        type(GGTMFieldlineRefinementOptionsUDT)     :: linerefoptions
        type(GGTMFieldlineDataUDT)                  :: hfline, lfline
        integer(I8)                                 :: srflabel, erflabel, &
            srf, erf, region
        integer(I8), allocatable, dimension(:)      :: srfvert, erfvert, &
            hffaces, lffaces, hfvert, lfvert
        logical                                     :: flipsrf, fliperf
        character(:), allocatable                   :: legalcellstyle

    end type

    ! Tube data
    type :: GGTMTubeDataUDT

        ! Description
        !============
        ! Contains additional flux tube data, such as field value 
        ! data etc 

        real(R8), allocatable, dimension(:)     :: fval
        integer(I8)                             :: distributionface
        type(GGTMFieldlineRefinementOptionsUDT) :: linerefoptions

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
        type(GGTMSegmentUDT), allocatable, dimension(:)     :: seg
        integer(I8)                                         :: nseg 

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
        ! type is deferred to a dedicated function. Options for line 
        ! refinement that are specific to a certain refiner can be 
        ! set by calling the 'UpdateRefinementOptions' routine (which is
        ! deferred) that takes as input a GGTMCellData structure and the
        ! topological mesh (see also abstract interface of the 
        ! routine). This should be enough information to update the 
        ! desired refinement (otherwise to be extended in the future). 
        ! It is thus assumed that all refinement options are stored in 
        ! the refiner itself

    contains 

        ! Refine line (single line based)
        procedure(RefineGGTMLineSingleINT), deferred :: RefineLineSingle
        generic :: Refine   => RefineLineSingle

        ! Apply vertex distribution of a line to another 
        procedure(ProjectLineVertexDistributionINT), deferred   :: &
            ProjectLineVertexDistribution 

        ! Update refinement options
        procedure(UpdateRefinementOptionsINT), deferred :: UpdateRefinementOptions

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

        ! Update refinement options
        procedure :: UpdateRefinementOptions    => UpdateRefinementOptionsNoRef

        ! Project distribution
        procedure :: ProjectLineVertexDistribution  => &
            ProjectLineVertexDistributionNoRef

    end type

    ! Length based refiner 
    type, extends(GGTMLineRefiner2DUDT)     :: GGTMLineRefinerLB2DUDT

        ! Description
        !============
        ! This refiner is based on a minimal and maximal length 
        ! distribution. Additionally, one can specify a boundary layer 
        ! at each side of the face where a specified length (or set of 
        ! lengths) is specified. 

        ! List of options:
        ! - meth            method of refining/coarsening, currently
        !                   only 'classic' which simply splits up a face
        !                   in two or merges a face with a neighbouring 
        !                   face (boundary layers are not included)
        ! - doBLstart       do starting (or ending if doBLend) boundary
        !                   layer
        ! - dlBLlengthbased switch to determine if it is based on 
        !                   classical (euler) length (true) or on the 
        !                   lengthtype (false)
        ! - ncBLstart       number of boundary layer cells 
        ! - dlBLstart       size of boundary layer cells (in length units 
        !                   specified by lengthtype!)
        ! - lengthtype      type of length to consider when doing 
        !                   refinement. Can be 'euler' (default classic
        !                   length), 'radial' (length, projected
        !                   radially and in absolute value)
        ! - field           field description (magnetic field like) on 
        !                    which field related lengths are based
        ! - linedllc        line length distribution, but in terms of 
        !                   lengthtype
        character(:), allocatable           :: meth, lengthtype 
        logical                             :: doBLstart, &
            doBLend, dlBLlengthbased
        integer(I8)                         :: ncBLstart, ncBLend 
        real(R8), allocatable, dimension(:) :: dlBLstart, dlBLend, &
            linedllc 
        class(DistributionFunctionUDT), allocatable     :: Lmin, Lmax 
        type(MagneticFieldUDT)              :: field

    contains 

        ! Refine line
        procedure :: RefineLineSingle               => RefineLineSingleLB

        ! Update refinement options
        procedure :: UpdateRefinementOptions        => UpdateRefinementOptionsLB

        ! Project distribution
        procedure :: ProjectLineVertexDistribution  => &
            ProjectLineVertexDistributionLB

        ! Auxiliary
        procedure :: GetLineEdgeLength              => GetLineEdgeLengthLB
        procedure :: AddLineVertexCoordinates       => AddLineVertexCoordinatesLB   
        procedure :: InitializeLineData             => InitializeLineDataLB

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

        ! Update refinement options
        subroutine UpdateRefinementOptionsINT(refiner, refoptions, &
            topomesh)

            import :: GGTMLineRefiner2DUDT, GGTMFieldlineRefinementOptionsUDT, &
                TopomeshUDT
            class(GGTMLineRefiner2DUDT)                         :: refiner 
            type(GGTMFieldlineRefinementOptionsUDT), intent(in) :: refoptions 
            type(TopomeshUDT), intent(in)                       :: topomesh

        end subroutine

        ! Project vertex distribution
        subroutine ProjectLineVertexDistributionINT(refiner, linein, &
            lineout, vertID, ggtmdata)
            import :: GGTMFieldlineDataUDT, I8, GGTMLineRefiner2DUDT, &
                GGTMDataUDT
            class(GGTMLineRefiner2DUDT)                 :: refiner 
            type(GGTMFieldlineDataUDT), intent(in)      :: linein 
            type(GGTMFieldlineDataUDT), intent(inout)   :: lineout
            integer(I8), intent(inout)                  :: vertID
            type(GGTMDataUDT), intent(inout)            :: ggtmdata 
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
        real(R8)                    :: valplf, xb(1:2), yb(1:2)
        real(R8), allocatable, dimension(:)         :: xp, yp, dp, valp
        integer(I8)                 :: resx, resy
        integer(I8), allocatable, dimension(:)  :: xpind, spind
        type(GGTMDataUDT)           :: ggtmdata 
        class(VertexDistributor2DUDT), allocatable      :: &
            poloidalvertexdistributor, radialvertexdistributor
        class(DistributionFunctionUDT), allocatable     :: & 
            magneticFieldDF, vdpdensityfunction, vdrdensityfunction
        class(StreamlineTracerUDT), intent(inout)   :: streamlinetracer
        class(GGTMLineRefiner2DUDT), allocatable    :: GGTMlinerefinerpol, &
            GGTMlinerefinerrad
        class(PolygonLevelsetFunction2DUDT), allocatable    :: vdpplf
        type(GGGridUDT)             :: grid 

        ! Loop
        integer(I8)                 :: i

        ! Initialize
        !===========
        ! Set verbosity
        verbosity  = options%verbosity 

        ! Set plotting bounds
        xb = [minval(magneticField%R), maxval(magneticField%R)]
        yb = [minval(magneticField%Z), maxval(magneticField%Z)]
        resx = 100
        resy = 100

        ! Required data of topomesh for grid generator
        call ggtmdata%Initialize(topomesh)

        ! Magnetic field distribution function
        magneticFieldDF = ConstructStructured2DDF(magneticField%interp)

        ! Temporary grid structure
        call grid%Initialize('standard')

        ! Set up vertex distribution 'actors'
        !====================================
        ! Poloidal vertex distributor
        select case (options%vdptype)

        case ('uniform')

            ! Construct uniform distributor with facelength 'options%vdpdfacelength'
            poloidalvertexdistributor = ConstructUniformVertexDistributor(&
                options%vdpdfacelength, options%vdrdfieldwidth)

        case ('densitybased')

            ! Construct distribution function based on vessel and 
            ! additional points. 

            ! Associate for ease
            associate(&
                decaylengthplf  => options%vdpddecaylengthplf,  &
                decaylengthxp   => options%vdpddecaylengthxp,   &
                valxp           => options%vdpddensityatxp,     &
                valinf          => options%vdpddensityatinf     &
            )

            ! Check which levelset function to use
            select case(options%vdpplftype)

            case ('vessel')

                vdpplf = vessel%plfvessel
                valplf = options%vdpddensityatvessel

            case ('target')

                vdpplf = vessel%plftarget
                valplf = options%vdpddensityatvessel

            case ('none')

                ! Just take vessel PLF but set weight to zero
                vdpplf = vessel%plftarget 
                valplf = 0.0_R8
                
            case default

                call gdErrorHandler('GenerateUnstructuredAlignedGrid: ' // & 
                    'unknown polygon levelset function choice: ' // options%vdpplftype)

            end select

            ! Check which points to use
            xp      = options%vdpdx 
            yp      = options%vdpdy 
            valp    = options%vdpdval 
            dp      = options%vdpdd

            if (options%vdpdincludexp) then 
                ! Add all X-points
                do i = 1, topomesh%vert%ntot 
                    if (topomesh%vert%type(i) == TMvertexsaddleID) then 
                        xp = [xp, topomesh%vert%x(i)]
                        yp = [yp, topomesh%vert%y(i)]
                        valp = [valp, options%vdpddensityatxp]
                        dp = [dp, options%vdpddecaylengthxp]
                    end if 
                end do 
            end if 

            ! Construct density function
            vdpdensityfunction = ConstructCoordinatesPLF2DDistanceDF(&
                vdpplf, valplf, decaylengthplf, xp, yp, valp, valinf, dp)

            ! Visualize
            call vdpdensityfunction%Visualize(xb, yb, resx, resy, &
                'gg_vd_poloidaldensityfunction')

            ! Construct density based distribution function
            poloidalvertexdistributor = ConstructDensityBasedVertexDistributor(vdpdensityfunction, 1_I8)

            ! Housekeeping
            !=============
            deallocate(xp, yp, valp, dp)
            end associate


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

        case ('densitybased')

            ! Construct density based distributor
            
            ! Initialize
            associate(&
                decaylengthxp   => options%vdrddecaylengthxp,    &
                valxp           => options%vdrddensityatxp,      &
                valinf          => options%vdrddensityatinf      &
                )

            allocate(xp(0), yp(0), valp(0), dp(0))

            ! Add user-defined points
            xp = [xp, options%vdrdx]
            yp = [yp, options%vdrdy]
            dp = [dp, options%vdrdd]
            valp = [valp, options%vdrdval]
            

            ! Check which points to use
            if (options%vdrdoxp) then 
                ! Include refinement near x-points (and their separatrices)
                do i = 1, topomesh%face%ntot
                    ! Add separatrix points
                    if (topomesh%face%type(i) == TMfacesepID) then 
                        ! Add all points except the end points (added 
                        ! later to avoid duplication)
                        associate(tpol => topomesh%face%pol(i))
                        xp = [xp, tpol%x(tpol%vert(2:size(tpol%vert)-1))]
                        yp = [yp, tpol%y(tpol%vert(2:size(tpol%vert)-1))]
                        valp = [valp, spread(valxp, 1, tpol%ne-1)]
                        dp = [dp, spread(decaylengthxp, 1, tpol%ne-1)]
                        end associate
                    end if

                    ! Add x- and o-points
                    xpind = topomesh%GetXPointIDs()
                    spind = topomesh%GetStrikePointIDs()
                    xp = [xp, topomesh%vert%x([xpind, spind])]
                    yp = [yp, topomesh%vert%y([xpind, spind])]
                    valp = [valp, spread(valxp, 1, size([xpind, spind]))]
                    dp = [dp, spread(decaylengthxp, 1, size([xpind, spind]))]
                end do  
            end if 

            ! Construct density function 
            vdrdensityfunction = ConstructCoordinates2DDistanceDF(&
                xp, yp, valp, valinf, dp)

            ! Visualize
            call vdrdensityfunction%Visualize(xb, yb, resx, resy, &
                'gg_vd_radialdensityfunction')

            ! Construct density based distribution function
            radialvertexdistributor = ConstructDensityBasedVertexDistributor(vdrdensityfunction, 1_I8)

            ! Housekeeping
            end associate
            

        case Default

            ! Unknown option
            call gdErrorHandler('GenerateUnstructuredAlignedGrid: ' // & 
                'radial vertex distribution option: ' // options%vdrtype // &
                ' not implemented')

        end select

        ! Refiner (poloidal)
        GGTMlinerefinerpol = InitializeGGTMLineRefiner(topomesh, &
            magneticField, vessel, fieldtracer, boundarytracer, &
            poloidalvertexdistributor, radialvertexdistributor, options, &
            'poloidal')

        ! Refiner (radial)
        GGTMlinerefinerrad = InitializeGGTMLineRefiner(topomesh, &
            magneticField, vessel, fieldtracer, boundarytracer, &
            poloidalvertexdistributor, radialvertexdistributor, options, &
            'radial')


        ! Add gridding data (and additional refinement options)
        call AddTopologicalMeshGriddingData(ggtmdata, topomesh, &
            fieldtracer, magneticField, options)

        ! Distribute vertices on topological faces
        !=========================================
        ! Initialize the vertex counter (normally, all topomesh vertices 
        ! are added)
        grid%vert%ntot = topomesh%vert%ntot

        ! Poloidal faces
        call DistributeVerticesTopologicalMeshFaces(grid, ggtmdata, topomesh, &
            poloidalvertexdistributor, GGTMlinerefinerpol, TMfacealignedID)

        ! Relevant radial faces of tubes
        call DistributeVerticesTopologicalMeshTubes(grid, ggtmdata, &
            topomesh, radialvertexdistributor, GGTMlinerefinerrad, &
            magneticFieldDF)

        ! Generate initial grid
        !======================
        ! Trace contours
        call TraceTopologicalMeshTubeContours(grid, ggtmdata, topomesh, &
            fieldtracer, magneticField, options)

        ! Generate elemental flux tubes for gridding
        call ConstructTopologicalMeshCellFluxTubes(grid, ggtmdata, topomesh, &
            fieldtracer, magneticField, options)
        
        ! Distribute vertices 
        select case (options%TMcellgriddingorder)

        case ('independent')

            call DistributeVerticesIndependent(ggtmdata, topomesh, grid, &
                poloidalvertexdistributor, magneticField, streamlinetracer, &
                GGTMlinerefinerpol, options)

        case ('sequential')

            call DistributeVerticesSequential(ggtmdata, topomesh, grid, &
                poloidalvertexdistributor, magneticField, streamlinetracer, &
                GGTMlinerefinerpol, options)

        case default 

            call gdErrorHandler('GenerateUnstructuredAlignedGrid: ' // & 
                'unknown topological mesh cell distribution order option: ' & 
                 // options%TMcellgriddingorder)

        end select

        ! Write data
        call WriteGGTMData(ggtmdata, 'ggtmdata_before_pp')

        ! Post-process the distribution
        call PostProcessVertexDistribution(ggtmdata, topomesh, grid, &
            streamlinetracer)

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

        ! Apply grid adaptations
        !=======================
        ! Remove non-convex cells
        call SplitNonConvexCells(ggtmdata, grid) 

        ! 

        ! Extract the necessary gridding data
        !====================================
        ! Extract
        call ExtractSimulationGrid(simgrid, grid, magneticField, topomesh)

        ! Diagnostics
        !============
        call RunGridDiagnostics(simgrid)



        ! Visualize
        !----------
        !VisualizeGrid(grid, 1);


    end subroutine

    !------------------------------------------------------------------!
    !                         GRID GENERATION                          !
    !------------------------------------------------------------------!

    ! Independent gridder (loop for multiple cells)
    subroutine DistributeVerticesIndependent(ggtmdata, topomesh, grid, &
        vd, magneticField, streamlinetracer, GGTMlinerefiner, options)

        ! Description
        !============
        ! This routine generates for each topological cell a grid in 
        ! an independent way. It is assumed that grid coordinates are
        ! already present on all topological mesh faces and that
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
        type(MagneticFieldUDT), intent(in)          :: magneticField
        class(StreamlineTracerUDT), intent(in)      :: streamlinetracer
        class(GGTMLineRefiner2DUDT), intent(in)     :: GGTMlinerefiner
        type(GGoptionsUDT), intent(in)              :: options

        ! Auxiliary
        integer(I8)                                 :: vertID
        integer(I8), allocatable, dimension(:)      :: allIDs, &
            vertmap
        logical, allocatable, dimension(:)          :: isvertexdeleted


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
        ! Set vertID - assumed this was already initialized in a previous
        ! step
        vertID = grid%vert%ntot

        ! Add cell vertices
        !==================
        do i = 1, cell%ntot

            ! Update the refiner
            call GGTMLineRefiner%UpdateRefinementOptions(celldata(i)%linerefoptions, &
                topomesh)

            ! Check which method to use
            select case (options%ggmethod)

            case ('independent')

                ! Simply grid the cell in an independent way - don't 
                ! regrid 
                call DistributeVerticesSingleTMCellIndependent(i, vertID, &
                    .false., .false., ggtmdata, topomesh, vd, GGTMLineRefiner)

            case ('orthogonal')

                ! Grid in an orthogonal way, possibly with refinement,
                ! but not updating original line distribution...
                call DistributeVerticesSingleTMCellOrthogonal(i, vertID, &
                    .true., ggtmdata, topomesh, grid, vd, &
                    magneticField, streamlinetracer, GGTMlinerefiner)

            case default

                ! Unknown
                call gdErrorHandler('DistributeVerticesIndependent: ' // &
                    'unknown vertex distribution method: ' // options%ggmethod)

            end select

        end do 

        ! Post-process
        !=============
        ! Still need to account for deleted/overwritten vertices
        allocate(isvertexdeleted(vertID))
        isvertexdeleted = .true. 
        do i = 1, cell%ntot 
            ! Check if vertices are present & ensure proper updated lines
            do j = 1, size(celldata(i)%tubes)
                call celldata(i)%tubes(j)%hfline%UpdateLineData(ggtmdata)
                call celldata(i)%tubes(j)%lfline%UpdateLineData(ggtmdata)
                isvertexdeleted(celldata(i)%tubes(j)%hfline%vert) = .false.
                isvertexdeleted(celldata(i)%tubes(j)%lfline%vert) = .false.
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

        ! Loop and adjust IDs - first segments, then tubes
        do i = 1, ggtmdata%nseg 
            ggtmdata%seg(i)%vert = vertmap(ggtmdata%seg(i)%vert)
            ggtmdata%seg(i)%sv = vertmap(ggtmdata%seg(i)%sv)
            ggtmdata%seg(i)%ev = vertmap(ggtmdata%seg(i)%ev)
        end do 
        do i = 1, cell%ntot 
            ! Update
            do j = 1, size(celldata(i)%tubes)
                call celldata(i)%tubes(j)%hfline%UpdateLineData(ggtmdata)
                call celldata(i)%tubes(j)%lfline%UpdateLineData(ggtmdata)
            end do 
        end do 

        ! Update number of grid vertices
        grid%vert%ntot = count(.not. isvertexdeleted)

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Sequential gridder
    subroutine DistributeVerticesSequential(ggtmdata, topomesh, &
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
        use mod_definitions, only: TMfacealignedID, &
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
        integer(I8)                                 :: vertID, tc
        integer(I8), allocatable, dimension(:)      :: &
            allIDs, vertmap
        logical, allocatable, dimension(:)          :: iscelldone, &
            isfacedone, isstartingcell, isstartingface, &
            isvertexdeleted

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

        ! Initialize
        allocate(iscelldone(cell%ntot), isfacedone(face%ntot), &
            isstartingcell(cell%ntot), isstartingface(face%ntot))
        iscelldone = .false. 
        isfacedone = .false. 
        isstartingcell = .false.
        isstartingface = .false.

        ! Set vertID
        vertID = grid%vert%ntot

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
                if (face%BF(hffaces(j)) .and. any(face%type(hffaces(j)) == TMfacealignedID)) then 
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

            ! Add cell vertices
            !------------------
            ! Update the refiner
            call GGTMLineRefiner%UpdateRefinementOptions(&
                celldata(tc)%linerefoptions, topomesh)

            ! Check which method to use
            select case (options%ggmethod)

            case ('independent')

                ! Simply grid the cell in an independent way - don't 
                ! regrid 
                call DistributeVerticesSingleTMCellIndependent(tc, vertID, &
                    .false., .true., ggtmdata, topomesh, vd, GGTMLinerefiner)

            case ('orthogonal')

                ! Grid in an orthogonal way, possibly with refinement,
                ! but not updating original line distribution...
                call DistributeVerticesSingleTMCellOrthogonal(tc, vertID, &
                    .true., ggtmdata, topomesh, grid, vd, &
                    magneticField, streamlinetracer, GGTMlinerefiner)

            case default

                ! Unknown
                call gdErrorHandler('DistributeVerticesIndependent: ' // &
                    'unknown vertex distribution method: ' // options%ggmethod)

            end select

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
            ! Check if vertices are present & ensure properly updated lines
            do j = 1, size(celldata(i)%tubes)
                call celldata(i)%tubes(j)%hfline%UpdateLineData(ggtmdata)
                call celldata(i)%tubes(j)%lfline%UpdateLineData(ggtmdata)
                isvertexdeleted(celldata(i)%tubes(j)%hfline%vert) = .false.
                isvertexdeleted(celldata(i)%tubes(j)%lfline%vert) = .false.
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

        ! Loop and adjust IDs - first segments, then tubes
        do i = 1, ggtmdata%nseg 
            ggtmdata%seg(i)%vert = vertmap(ggtmdata%seg(i)%vert)
            ggtmdata%seg(i)%sv = vertmap(ggtmdata%seg(i)%sv)
            ggtmdata%seg(i)%ev = vertmap(ggtmdata%seg(i)%ev)
        end do 
        do i = 1, cell%ntot 
            ! Update
            do j = 1, size(celldata(i)%tubes)
                call celldata(i)%tubes(j)%hfline%UpdateLineData(ggtmdata)
                call celldata(i)%tubes(j)%lfline%UpdateLineData(ggtmdata)
            end do 
        end do 

        ! Update number of grid vertices
        grid%vert%ntot = count(.not. isvertexdeleted)

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Single cell crap gridder
    subroutine DistributeVerticesSingleTMCellIndependent(tc, vertID, &
        dohfline, dolfline, ggtmdata, topomesh, vd, GGTMlinerefiner)

        ! Description
        !============
        ! Distribute the grid cells on a single topological mesh cell in 
        ! an independent way. This means that grid vertices are 
        ! distributed first on all field lines after which cells and 
        ! faces are generated. Although we only grid one cell, we need
        ! to pass the entire ggtmdata structure since we may need
        ! to update data afterwards that is not local to the cell. 
        ! The cell index to be gridded should thus be given as input.
        ! Additionally, one can specifiy whether the high field and 
        ! low field lines of the cell should be regridded (existing
        ! distributions are then overwritten, so handle with care!)

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)                     :: tc
        integer(I8), intent(inout)                  :: vertID 
        type(GGTMDataUDT), intent(inout)            :: ggtmdata
        class(TopomeshUDT), intent(in)              :: topomesh
        logical, intent(in)                         :: dohfline, &
            dolfline
        class(VertexDistributor2DUDT), intent(in)   :: vd
        class(GGTMLineRefiner2DUDT), intent(in)     :: GGTMlinerefiner

        ! Auxiliary
        integer(I8)                                 :: nv 
        integer(I8), allocatable, dimension(:)      :: tvID
        real(R8), allocatable, dimension(:)         :: dlcv
        logical, allocatable, dimension(:)          :: keepvert

        ! Loop
        integer(I8)                                 :: i, j, k 


        ! Initialize
        !===========
        ! Associate
        associate(&
            seg             => ggtmdata%seg,    &
            vert            => topomesh%vert,   &
            face            => topomesh%face,   &
            facedata        => ggtmdata%face,   &
            cell            => topomesh%cell,   &
            celldata        => ggtmdata%cell,   &
            tubes           => ggtmdata%cell(tc)%tubes  &
            )

        ! Determine cell vertices
        !========================
        ! Compute vertices 
        do i = 1, size(tubes)
            ! Unpack
            !-------
            associate(&
                hfline      => tubes(i)%hfline, &
                lfline      => tubes(i)%lfline  &
                )
            
            ! High field line
            !----------------
            ! Distribute if dohfline and first tube
            if (i == 1 .and. dohfline) then 
                ! Distribute over segments
                if (.not. seg(hfline%segID(1))%isvertex) then
                    do j = 1, size(hfline%segID)
                        ! Distribute over segment
                        call vd%DistributeOverCurve(seg(hfline%segID(j))%xl, &
                            seg(hfline%segID(j))%yl, nv, ldistr=dlcv)
                        tvID = [(k, k = vertID+1, vertID+size(dlcv)-2)]
                        call seg(hfline%segID(j))%AddVertices(dlcv(2:size(dlcv)-1), tvID)
                    
                        ! Update vertID
                        vertID = vertID + size(dlcv) - 2
                    end do

                    ! Refresh hfline 
                    call hfline%UpdateLineData(ggtmdata)

                    ! Refine/coarsen
                    keepvert = hfline%isnodevert ! keep vertices on nodes
                    call GGTMLinerefiner%Refine(hfline, vertID, keepvert)

                    ! Update hfline segments
                    call hfline%UpdateSegmentData(ggtmdata)
                end if 
            else
                ! Refresh hfline - should be gridded already
                call hfline%UpdateLineData(ggtmdata)
            end if 

            ! Low field line
            !---------------
            ! Distribute over segments
            if ((i == size(tubes) .and. dolfline) .or. (i < size(tubes))) then
                ! Distribute over segments
                if (.not. seg(lfline%segID(1))%isvertex) then 
                    do j = 1, size(lfline%segID)
                        ! Distribute over segment
                        call vd%DistributeOverCurve(seg(lfline%segID(j))%xl, &
                            seg(lfline%segID(j))%yl, nv, ldistr=dlcv)
                        tvID = [(k, k = vertID+1, vertID+size(dlcv)-2)]
                        call seg(lfline%segID(j))%AddVertices(dlcv(2:size(dlcv)-1), tvID)
                    
                        ! Update vertID
                        vertID = vertID + size(dlcv) - 2
                    end do

                    ! Refresh lfline 
                    call lfline%UpdateLineData(ggtmdata)

                    ! Refine/coarsen
                    if (.not. seg(lfline%segID(1))%isvertex) then 
                        keepvert = lfline%isnodevert ! keep vertices on nodes
                        call GGTMLinerefiner%Refine(lfline, vertID, keepvert)
                    end if 

                    ! Update hfline segments
                    call lfline%UpdateSegmentData(ggtmdata)
                end if 
            else
                ! Refresh lfline - should be gridded already
                call lfline%UpdateLineData(ggtmdata)
            end if 

            ! Housekeeping
            end associate
        end do 

        ! Housekeeping
        end associate

    end subroutine

    ! Single cell orthogonal gridder 
    subroutine DistributeVerticesSingleTMCellOrthogonal(tc, vertID, &
        dolfline, ggtmdata, topomesh, grid, vd, &
        magneticField, streamlinetracer, GGTMlinerefiner)

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
        use mod_definitions, only: TMfacealignedID, &
            TMvertexbndID, TMvertextp1ID, TMvertextp2ID


        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)                     :: tc 
        integer(I8), intent(inout)                  :: vertID
        logical, intent(in)                         :: dolfline
        class(GGTMDataUDT)                          :: ggtmdata 
        class(TopomeshUDT), intent(in)              :: topomesh 
        class(GGGridUDT), intent(inout)             :: grid 
        class(VertexDistributor2DUDT), intent(in)   :: vd
        type(MagneticFieldUDT), intent(in)          :: magneticField
        class(StreamlineTracerUDT), intent(in)      :: streamlinetracer
        class(GGTMLineRefiner2DUDT), intent(in)     :: GGTMlinerefiner

        ! Auxiliary
        integer(I8)                                 :: tf
        integer(I8), allocatable, dimension(:)      :: &
            s11, s12, s21, s22, s31, s32, s41, s42, &
            stype, sortind, newID
        logical                                     :: addpoint
        logical, allocatable, dimension(:)          ::  keepvert, newisnodevert, &
            updateseg
        real(R8)                                    :: xb(1:2), yb(1:2)
        real(R8), allocatable, dimension(:)         :: xt, yt, &
            x1, x2, x3, x4, y1, y2, y3, y4, s11r, s12r, s21r, s22r, &
            s31r, s32r, s41r, s42r, s1r, temps2r, tempx, tempy, newtx, &
            newty, news2r, newdlcv

        type(StreamlineUDT), allocatable            :: orthlines(:)

        ! Loop
        integer(I8)                                 :: i, j, k, nnew

        ! Initialize
        !===========
        ! Associate
        associate(&
            seg             => ggtmdata%seg,    &
            vert            => topomesh%vert,   &
            face            => topomesh%face,   &
            facedata        => ggtmdata%face,   &
            cell            => topomesh%cell,   &
            celldata        => ggtmdata%cell,   &
            tubes           => ggtmdata%cell(tc)%tubes  &
            )

        ! Determine cell vertices
        !========================
        ! Compute intersections;
        do i = 1, size(tubes)
            ! Initialize
            !-----------
            ! Unpack
            associate(&
                hfline      => tubes(i)%hfline,     &
                lfline      => tubes(i)%lfline      &
                )

            ! Check if we need to update the original hfline
            if (i == 1) then 
                
                ! Update the aligned parts of the hf line
                do j = 1, size(celldata(tc)%hffaces)
                    ! Unpack
                    tf = celldata(tc)%hffaces(j)

                    ! Update the face data to propagate previous distribution
                    call facedata(tf)%line%UpdateLineData(ggtmdata) 

                    ! Update the refiner
                    call GGTMLineRefiner%UpdateRefinementOptions(&
                        facedata(tf)%linerefoptions, topomesh)

                    ! Refine
                    keepvert = IsTopomeshVert(facedata(tf)%line%vert, topomesh)
                    call GGTMlinerefiner%Refine(facedata(tf)%line, vertID, keepvert)

                    ! Update segment data
                    call facedata(tf)%line%UpdateSegmentData(ggtmdata) 
                    
                end do 

                ! Update refinement data 
                call GGTMLineRefiner%UpdateRefinementOptions(&
                celldata(tc)%linerefoptions, topomesh)

                ! Refine/coarsen non-aligned parts
                keepvert = hfline%isnodevert
                call GGTMLinerefiner%Refine(hfline, vertID, keepvert)

                ! Update segments, but only those that are not aligned 
                ! (i.e. no fsID)
                updateseg = hfline%fsID == 0
                call hfline%UpdateSegmentData(ggtmdata, updateseg)
                call hfline%UpdateLineData(ggtmdata)

            else

                ! Update line data
                call hfline%UpdateLineData(ggtmdata)

                ! Refine/coarsen 
                keepvert = hfline%isnodevert
                call GGTMLinerefiner%Refine(hfline, vertID, keepvert)

                ! Update segments
                call hfline%UpdateSegmentData(ggtmdata)

            end if 
            
            ! Skip if the lfline is a vertex
            if (seg(lfline%segID(1))%isvertex) then 
                cycle
            end if 

            ! Trace streamlines starting from the hfline
            xt = hfline%xv
            yt = hfline%yv 
            xb = [minval([lfline%xl, hfline%xl]), &
                maxval([lfline%xl, hfline%xl])]
            yb = [minval([lfline%yl, hfline%yl]), &
                maxval([lfline%yl, hfline%yl])]
            orthlines = streamlinetracer%TraceStreamlines(xt, yt, &
                xb, yb, spread(-1_I8, 1, size(xt))) ! normally, gradient goes from low to high, so need to reverse sign

            ! Initialize potential new coordinates
            allocate(newtx(size(orthlines)), newty(size(orthlines)), &
                news2r(size(orthlines)), newID(size(orthlines)))
            nnew = 0

            ! Find intersections with all other boundaries
            !$omp parallel default(private) shared(orthlines, vertID, ggtmdata, &
            !$omp newtx, nnew, newty, news2r, newID, i) 
            !$omp do schedule(dynamic)
            do j = 1, size(orthlines)
                ! Initialize
                addpoint = .true. 

                ! Check for starting and ending point being the same
                if (j == 1) then 
                    if (hfline%vert(1) == lfline%vert(1)) then 
                        addpoint = .false. 
                    end if 
                elseif (j == hfline%nv) then 
                    if (hfline%vert(hfline%nv) == lfline%vert(lfline%nv)) then 
                        addpoint = .false. 
                    end if 
                end if 

                ! Intersections with starting boundary (should only
                ! intersect in first point)
                if (seg(hfline%segID(1))%isvertex) then 
                    ! Previous boundary was point - start and end
                    ! should be the same
                    if ((hfline%xv(1) == orthlines(j)%x(1)) .and. &
                        (hfline%yv(1) == orthlines(j)%y(1))) then 
                        x1 = hfline%xv(1:1)
                        y1 = hfline%yv(1:1) 
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
                    orthlines(j)%y, hfline%xl, hfline%yl, &
                    x1, y1, s11, s12, s11r, s12r)
                end if 
                
                ! Intersections with ending boundary
                call SimplePolygonIntersections(orthlines(j)%x, &
                    orthlines(j)%y, lfline%xl, lfline%yl, &
                    x2, y2, s21, s22, s21r, s22r)
                
                ! Intersections with side boundary 1
                if (hfline%vert(1) /= lfline%vert(1)) then 
                    call SimplePolygonIntersections(orthlines(j)%x, &
                        orthlines(j)%y, [hfline%xl(1), lfline%xl(1)], &
                        [hfline%yl(1), lfline%yl(1)], &
                        x3, y3, s31, s32, s31r, s32r)

                    ! Check for intersections in starting point (may happen
                    ! if j == 1 or j == size(orhtlines))
                    if (j == 1 .or. j == size(orthlines)) then 
                        x3 = pack(x3, s31r /= 0_R8)
                        y3 = pack(x3, s31r /= 0_R8)
                        s32r = pack(x3, s31r /= 0_R8)
                        s31r = pack(x3, s31r /= 0_R8)
                    end if

                else
                    ! Point - no intersections
                    x3 = spread(0_R8, 1, 0)
                    y3 = spread(0_R8, 1, 0)
                    s31r = spread(0_R8, 1, 0)
                    s32r = spread(0_R8, 1, 0)
                end if 
                
                ! Intersections with side boundary 2
                if (hfline%vert(hfline%nv) /= lfline%vert(lfline%nv)) then 
                    call SimplePolygonIntersections(orthlines(j)%x, &
                        orthlines(j)%y, [hfline%xl(hfline%nl), lfline%xl(lfline%nl)], &
                        [hfline%yl(hfline%nl), lfline%yl(lfline%nl)], &
                        x4, y4, s41, s42, s41r, s42r)

                    ! Check for intersections in starting point (may happen
                    ! if j == 1 or j == size(orhtlines))
                    if (j == 1 .or. j == size(orthlines)) then 
                        x4 = pack(x4, s41r /= 0_R8)
                        y4 = pack(x4, s41r /= 0_R8)
                        s42r = pack(x4, s41r /= 0_R8)
                        s41r = pack(x4, s41r /= 0_R8)
                    end if
                else
                    ! Point - no intersections
                    x4 = spread(0_R8, 1, 0)
                    y4 = spread(0_R8, 1, 0)
                    s41r = spread(0_R8, 1, 0)
                    s42r = spread(0_R8, 1, 0)
                end if 
                
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
                ! ending boundary and should not be in its start or end
                if (addpoint) then 
                    if (size(s11r) == 0) then 
                        addpoint = .false.
                    elseif (size(stype) < 2) then 
                        ! Only one intersection found - don't add
                        call Write2DPolygonData(orthlines(j)%x, orthlines(j)%y, 'l1')
                        call Write2DPolygonData(lfline%xl, lfline%yl, 'l2')
                        call Write2DPolygonData(hfline%xl, hfline%yl, 'l3')
                        call Write2DPolygonData(hfline%xv, hfline%yv, 'l4')
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
            ! to include node vertices...)
            allocate(newdlcv(size(news2r)))
            newdlcv = 0_R8
            call Interpolate1D(news2r, newdlcv, &
                real([(k, k = 0, lfline%nl-1)], kind=R8), lfline%dllc)
            newdlcv = [newdlcv, lfline%dllc(lfline%nodes)]
            newID = [newID, pack(lfline%vert, lfline%isnodevert)]
            newisnodevert = [spread(.false., 1, size(news2r)), spread(.true., 1, lfline%ns+1)]

            ! Sort
            allocate(sortind(size(newdlcv)))
            sortind = 0_I8
            call Sort(newdlcv, ind=sortind, ascend=.true.)
            newID = newID(sortind)
            newisnodevert = newisnodevert(sortind)
            
            ! Add points
            call lfline%AddVertexCoordinates(newdlcv, ecbased=.false.)
            call lfline%AddVertexIDs(newID, newisnodevert)
        
            ! Refine/coarsen
            keepvert = lfline%isnodevert
            call GGTMLinerefiner%Refine(lfline, vertID, keepvert)

            ! Update segments
            call lfline%UpdateSegmentData(ggtmdata)

            ! Check if its the final line
            if (i == size(tubes)) then 
                ! Update the aligned parts of the lf line
                do j = 1, size(celldata(tc)%lffaces)
                    ! Unpack
                    tf = celldata(tc)%lffaces(j)

                    ! Update the face data to propagate previous distribution
                    call facedata(tf)%line%UpdateLineData(ggtmdata) 

                    ! Update the refiner
                    call GGTMLineRefiner%UpdateRefinementOptions(&
                        facedata(tf)%linerefoptions, topomesh)

                    ! Refine
                    keepvert = IsTopomeshVert(facedata(tf)%line%vert, topomesh)
                    call GGTMlinerefiner%Refine(facedata(tf)%line, vertID, keepvert)

                    ! Update segment data
                    call facedata(tf)%line%UpdateSegmentData(ggtmdata) 
                    
                end do 

                ! Update line data
                call lfline%UpdateLineData(ggtmdata)

                ! Update refinement data 
                call GGTMLineRefiner%UpdateRefinementOptions(&
                    celldata(tc)%linerefoptions, topomesh)

                ! Refine/coarsen non-aligned parts
                keepvert = lfline%isnodevert
                call GGTMLinerefiner%Refine(lfline, vertID, keepvert)

                ! Update segments, but only those that are not aligned 
                ! (i.e. no fsID)
                updateseg = lfline%fsID == 0
                call lfline%UpdateSegmentData(ggtmdata, updateseg)
                call lfline%UpdateLineData(ggtmdata)

            end if 

            ! Housekeeping
            deallocate(newtx, newty, newID, news2r, sortind, newdlcv)
            end associate
        end do

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Vertex distribution post-processing
    subroutine PostProcessVertexDistribution(ggtmdata, topomesh, grid, &
        streamlinetracer)

        ! Description
        !============
        ! This routine applies any post-processing operation to the 
        ! initial vertex distribution determined beforehand. The 
        ! routine should be called *before* adding the grid vertices, 
        ! since possibly vertices may still be deleted/adjusted 
        ! depending on circumstances. 

        ! Currently, we apply the following checks/add the following 
        ! information:
        ! - tube intersections (1)  we check whether lines of a tube 
        !                           after vertex construction intersect
        !                           with any other tube. We assume that 
        !                           the contour lines themselves do not
        !                           intersect (should be checked in 
        !                           tracing routine upstream), but that 
        !                           intersections arise from the discrete
        !                           vertex distribution and the field 
        !                           curvature. Then, we insert vertices 
        !                           near these intersections to capture 
        !                           the curvature better and hopefully 
        !                           remove the interesections. 
        ! - tube intersections (2)  we check whether the lines of a 
        !                           tube, *after* vertex construction, 
        !                           intersect. These tubes are removed
        !                           to prevent overlapping cells.
        ! - line-of-sight data:     data on LOS is added for vertices 
        !                           to facilitate cell construction in 
        !                           ConstructCellsQuadTria

        ! Declare variables
        !==================
        ! Arguments
        type(GGTMDataUDT), intent(inout)            :: ggtmdata 
        class(TopomeshUDT), intent(in)              :: topomesh
        type(GGGridUDT), intent(inout)              :: grid
        class(StreamlineTracerUDT), intent(in)      :: streamlinetracer
        
        ! Auxiliary
        integer(I8)                             :: nft, nct, t1, t2
        integer(I8), allocatable, dimension(:)  :: tc, s1, s2, allIDs, &
            vertmap, sortind, tracevert, ts1, ts2, &
            newvert
        real(R8)                                :: xb(1:2), yb(1:2), &
            tempr
        real(R8), allocatable, dimension(:)     :: xint, yint, s1r, s2r, &
            txint, tyint, ts1r, ts2r, tempdlcv, newdlcv
        logical                                 :: vertexwasdeleted, &
            foundIntersection
        logical, allocatable, dimension(:)      :: keepind, &
            isvertexdeleted, newisnodevert
        type(StreamlineUDT), allocatable        :: orthlines(:)
        character(6)                            :: tstring

        ! Loop
        integer(I8)                             :: i, j, k, cc, ic, is, &
            ie

        ! Initialize
        !===========
        ! Unpack for ease
        associate(&
            facedata        => ggtmdata%face,   &
            celldata        => ggtmdata%cell,   &
            tube            => topomesh%tube,   &
            cell            => topomesh%cell    &
            )

        ! Check if all necessary data is properly initialized
        do i = 1, size(celldata)
            ! Associate for ease
            associate(tc    => ggtmdata%cell(i))

            ! Check line pair allocation
            if (.not. allocated(tc%tubes)) then 
                ! Should have been allocated alread
                call gdErrorHandler('PostProcessVertexDistribution: ' // & 
                    'cell tubes have not yet been allocated')
            end if 

            ! Housekeeping
            end associate
        end do

        ! Initialize
        vertexwasdeleted = .false. 

        ! Check intersections (1)
        !========================
        ! At this stage, we attempt to remove intersections by adding
        ! grid vertices. The idea is that the contour lines themselves
        ! do not intersect, but that due to the discrete vertex 
        ! distribution and magnetic field curvature, the 'coarsened' 
        ! contours do intersect. By tracing orthogonal lines at vertices
        ! between intersection points, we can find additional vertices 
        ! to be inserted into the other field line to capture this 
        ! curvature. We keep looping over all tubes until no more 
        ! intersections are found. 
        
        ! Loop over all cells 
        !$omp parallel do default(none) private(i, k, xb, yb, xint, yint, &
        !$omp s1, s2, s1r, s2r, ic, tracevert, sortind, is, ie, orthlines, &
        !$omp newdlcv, newvert, newisnodevert, txint, tyint, ts1, ts2, &
        !$omp ts1r, ts2r, tempr, tempdlcv, foundIntersection, cc) &
        !$omp shared(ggtmdata, topomesh, grid, streamlinetracer) schedule(dynamic)
        do i = 1, cell%ntot
            ! Initialize
            allocate(xint(0), yint(0), s1(0), s2(0), s1r(0), s2r(0))

            ! Unpack for ease
            associate(tubes     => celldata(i)%tubes)

            ! Loop over all tubes
            k = 1
            do while (k <= size(tubes))
                ! Unpack
                associate(&
                    hfline      => tubes(k)%hfline,     &
                    lfline      => tubes(k)%lfline      &
                    )

                ! Update lines
                !$omp critical
                call hfline%UpdateLineData(ggtmdata)
                call lfline%UpdateLineData(ggtmdata)
                !$omp end critical

                ! Initialize
                foundIntersection = .false.

                ! Determine trace boxes
                xb = [minval([hfline%xl, lfline%xl]), maxval([hfline%xl, lfline%xl])]
                yb = [minval([hfline%yl, lfline%yl]), maxval([hfline%yl, lfline%yl])]

                ! Sanity check: if the tube line contours intersect, 
                ! skip tube (will likely be deleted downstream)
                call GGTMLineIntersections(ggtmdata, hfline, lfline, xint, yint, &
                    s1, s2, vertbased=.false.)
                if (size(xint) > 0) then
                    ! Print a warning
                    call Write2DPolygonData(hfline%xl, hfline%yl, 'l1')
                    call Write2DPolygonData(lfline%xl, lfline%yl, 'l2')
                    print *, 'cell: ', i, 'tube: ', k 
                    print *, 'PostProcessVertexDistribution: detected ' // & 
                        'self-intersecting tube, which should not occur at ' // & 
                        'this stage, unless tubes were extended that have ' // & 
                        'field lines that intersect multiple times with ' // & 
                        'the vessel boundary, or if contour coarsening has led to ' // & 
                        '(try turning this option off in the input file). ' // & 
                        'These tubes will be deleted ' // &
                        'and may lead to cell overlap and intersections...' 
                    ! Update counter
                    k = k + 1 

                    ! Skip remainder of loop
                    cycle
                end if 
                    
                ! Compute all vertex intersections
                call GGTMLineIntersections(ggtmdata, hfline, lfline, xint, yint, &
                    s1, s2, s1r=s1r, s2r=s2r, vertbased=.true.)

                ! Check if any intersections are present - if not, go 
                ! to the next tube
                if (size(xint) == 0) then
                    ! Update counter 
                    k = k + 1

                    ! Skip remainder of loop
                    cycle
                end if 

                call Write2DPolygonData(hfline%xv, hfline%yv, 'l1')
                call Write2DPolygonData(lfline%xv, lfline%yv, 'l2')

                ! Determine lfline tracing vertices
                allocate(tracevert(0))
                ic = 1 ! initialize intersection counter
                allocate(sortind(size(s1r)))
                call Sort(s1r, ind=sortind, ascend=.true.)
                s2r = s2r(sortind)
                s1 = s1(sortind)
                s2 = s2(sortind)
                deallocate(sortind)
                do while (.true.)
                    ! Check exit conditions
                    if (ic >= size(s1)) then 
                        exit 
                    end if 

                    ! Check if this intersection is in the same face 
                    ! as the next one
                    if (s1(ic) == s1(ic+1)) then 
                        ! Set ic as start index, find end
                        is = ic 
                        ie = findloc(s1, s1(ic), 1, back=.true.)

                        ! Add node for later tracing
                        tracevert = [tracevert, (cc, cc = ceiling(s2r(is))+1, ceiling(s2r(ie)))]

                        ! Update the counter
                        ic = ie + 1
                    else 
                        ic  = ic + 1
                    end if 
                end do

                ! Trace the streamlines from the lfline
                orthlines = streamlinetracer%TraceStreamlines(&
                    lfline%xv(tracevert), lfline%yv(tracevert), &
                    xb, yb, spread(1_I8, 1, size(tracevert)))

                ! Compute intersections with the hfline & build vertices
                newdlcv = hfline%dlcv
                newvert = hfline%vert
                newisnodevert = hfline%isnodevert
                do j = 1, size(orthlines)
                    ! Intersections
                    call SimplePolygonIntersections(orthlines(j)%x, &
                        orthlines(j)%y, hfline%xl, hfline%yl, &
                        txint, tyint, ts1, ts2, ts1r, ts2r)
                    call Write2DPolygonData(orthlines(j)%x, orthlines(j)%y, 'l3')

                    ! Check
                    if (size(txint) > 0) then 
                        ! Take the intersection closest to the start
                        ! of the orthogonal line
                        allocate(sortind(size(ts1r)))
                        call Sort(ts1r, ind=sortind, ascend=.true.)
                        ts2r = ts2r(sortind)
                        tempr = ts2r(1)
                        deallocate(sortind)

                        ! Interpolate to get actual length 
                        call Interpolate1D([tempr], tempdlcv, &
                            real([(cc, cc = 0, hfline%nl-1)], kind=R8), hfline%dllc)

                        ! Update grid vertex counter
                        grid%vert%ntot = grid%vert%ntot + 1

                        ! Add
                        newdlcv = [newdlcv, tempdlcv]
                        newvert = [newvert, grid%vert%ntot]
                        newisnodevert = [newisnodevert, .false.]                        

                        ! Update logicals
                        foundIntersection = .true.

                    else
                        ! A bit weird, but nothing tremendously wrong -
                        ! intersection will be caught downstream and 
                        ! tube will be removed
                    end if 
                end do


                ! Add vertices on the hfline
                allocate(sortind(size(newdlcv)))
                call Sort(newdlcv, ind=sortind, ascend=.true.)
                newvert = newvert(sortind)
                newisnodevert = newisnodevert(sortind)
                deallocate(sortind)

                !$omp critical
                call hfline%AddVertexCoordinates(newdlcv)
                call hfline%AddVertexIDs(newvert, newisnodevert)

                ! Update hfline segment data
                call hfline%UpdateSegmentData(ggtmdata)
                !$omp end critical 

                ! Housekeeping
                deallocate(tracevert)

                ! Determine hfline tracing vertices
                allocate(tracevert(0))
                ic = 1 ! initialize intersection counter
                allocate(sortind(size(s1r)))
                call Sort(s2r, ind=sortind, ascend=.true.)
                s1r = s1r(sortind)
                s1 = s1(sortind)
                s2 = s2(sortind)
                deallocate(sortind)
                do while (.true.)
                    ! Check exit conditions
                    if (ic >= size(s1)) then 
                        exit 
                    end if 

                    ! Check if this intersection is in the same face 
                    ! as the next one
                    if (s2(ic) == s2(ic+1)) then 
                        ! Set ic as start index, find end
                        is = ic 
                        ie = findloc(s2, s2(ic), 1, back=.true.)

                        ! Add node for later tracing
                        tracevert = [tracevert, (cc, cc = ceiling(s2r(is))+1, ceiling(s2r(ie)))]

                        ! Update the counter
                        ic = ie + 1
                    else 
                        ic = ic + 1
                    end if 
                end do

                ! Trace the streamlines from the hfline
                orthlines = streamlinetracer%TraceStreamlines(&
                    hfline%xv(tracevert), hfline%yv(tracevert), &
                    xb, yb, spread(-1_I8, 1, size(tracevert)))

                ! Compute intersections with the hfline & build vertices
                newdlcv = lfline%dlcv
                newvert = lfline%vert
                newisnodevert = lfline%isnodevert
                do j = 1, size(orthlines)
                    ! Intersections
                    call SimplePolygonIntersections(orthlines(j)%x, &
                        orthlines(j)%y, lfline%xl, lfline%yl, &
                        txint, tyint, ts1, ts2, ts1r, ts2r)

                    ! Check
                    if (size(txint) > 0) then 
                        ! Take the intersection closest to the start
                        ! of the orthogonal line
                        allocate(sortind(size(ts2r)))
                        call Sort(ts2r, ind=sortind, ascend=.true.)
                        ts1r = ts1r(sortind)
                        tempr = ts1r(1)
                        deallocate(sortind)

                        ! Interpolate to get actual length 
                        call Interpolate1D([tempr], tempdlcv, &
                            real([(cc, cc = 0, lfline%nl-1)], kind=R8), lfline%dllc)

                        ! Update grid vertex counter
                        grid%vert%ntot = grid%vert%ntot + 1

                        ! Add
                        newdlcv = [newdlcv, tempdlcv]
                        newvert = [newvert, grid%vert%ntot]
                        newisnodevert = [newisnodevert, .false.]

                        ! Update logicals
                        foundIntersection = .true.

                    else
                        ! A bit weird, but nothing tremendously wrong -
                        ! intersection will be caught downstream and 
                        ! tube will be removed
                    end if 
                end do


                ! Add vertices on the lfline
                allocate(sortind(size(newdlcv)))
                call Sort(newdlcv, ind=sortind, ascend=.true.)
                newvert = newvert(sortind)
                newisnodevert = newisnodevert(sortind)
                deallocate(sortind)

                !$omp critical
                call lfline%AddVertexCoordinates(newdlcv)
                call lfline%AddVertexIDs(newvert, newisnodevert)

                ! Update lfline segment data
                call lfline%UpdateSegmentData(ggtmdata)
                !$omp end critical

                ! Housekeeping
                deallocate(tracevert)

                ! Update the tube index
                if (foundIntersection) then 
                    ! Check if we can decrease the counter to recheck
                    ! the previous tube
                    if (k > 1) then 
                        k = k - 1
                    else
                        ! Keep k at one and recheck the tube
                        k = 1
                    end if 
                else
                    ! Increase the counter
                    k = k + 1
                end if 

                ! Housekeeping
                end associate
            end do 

            ! Housekeeping
            deallocate(xint, yint, s1, s2, s1r, s2r)
            end associate
        end do 
        !$omp end parallel do
        
        ! Check intersections (2) 
        !========================
        ! Lines of tube may intersect with the starting or ending
        ! field line or neighbouring tubes
        allocate(xint(0), yint(0), s1(0), s2(0))
        do i = 1, tube%ntot
            ! Initialize
            tc = tube%GetCell(i)
            nft = size(celldata(tc(1))%tubes) ! number of tubes
            allocate(keepind(nft))
            keepind = .true.

            ! Loop over all cells of this tube
            do j = 1, size(tc)
                ! Sanity check
                if (nft /= size(celldata(tc(j))%tubes)) then 
                    print *, 'tube: ', i, 'cell: ', j
                    call gdErrorHandler('ConstructTopologicalMEshCellFluxTubes: ' // & 
                        'cell within same topological mesh tube has different ' // & 
                        'number of flux tubes, not supported')
                end if 

                ! Associate for ease
                associate(ct        => celldata(tc(j))%tubes)

                ! Check intersections of lfline (except if final tube)
                nct = size(ct)
                !$omp parallel do default(none) schedule(dynamic) &
                !$omp private(k, xint, yint, s1, s2) & 
                !$omp shared(nct, ggtmdata, j, keepind)
                do k = 1, nct
                    ! Associate for ease
                    associate(&
                        hfline1     => ct(1)%hfline,    &
                        lfline1     => ct(1)%lfline,    &
                        hflinen     => ct(nct)%hfline,  &
                        lflinen     => ct(nct)%lfline,  &
                        hflinek     => ct(k)%hfline,    &
                        lflinek     => ct(k)%lfline     &
                        )

                    ! Update just to be sure
                    !$omp critical
                    call hfline1%UpdateLineData(ggtmdata)
                    call hflinen%UpdateLineData(ggtmdata)
                    call hflinek%UpdateLineData(ggtmdata)
                    call lfline1%UpdateLineData(ggtmdata)
                    call lflinen%UpdateLineData(ggtmdata)
                    call lflinek%UpdateLineData(ggtmdata)
                    !$omp end critical

                    if (k == 1) then 
                        ! Only need to check the lfline
                        if (keepind(k)) then 
                            ! Compute intersections
                            call GGTMLineIntersections(ggtmdata, lflinek, hfline1, &
                                xint, yint, s1, s2, vertbased=.true.)

                            ! Check
                            if (size(xint) > 0) then
                                !$omp critical
                                keepind(k) = .false.
                                !$omp end critical
                            end if 
                        end if 

                        if (nct > 1) then 
                            ! Check for intersections with last lfline
                            if (keepind(k)) then 
                                ! Use dedicated routine to hedge for end point
                                ! intersections
                                call GGTMLineIntersections(ggtmdata, lflinek, lflinen, &
                                    xint, yint, s1, s2, vertbased=.true.)

                                ! Check
                                if (size(xint) > 0) then
                                    !$omp critical
                                    keepind(k) = .false.
                                    !$omp end critical
                                end if 
                            end if 

                            ! Check for intersections with next lfline
                            if (keepind(k)) then 
                                ! Use dedicated routine to hedge for end point
                                ! intersections
                                call GGTMLineIntersections(ggtmdata, lflinek, ct(k+1)%lfline, &
                                    xint, yint, s1, s2, vertbased=.true.)

                                ! Check
                                if (size(xint) > 0) then
                                    !$omp critical
                                    keepind(k) = .false.
                                    !$omp end critical
                                end if 
                            end if 
                        end if 
                    elseif (k == nct .and. nct > 1) then 
                        ! Only need to check the hfline

                        ! First hfline
                        if (keepind(k)) then 
                            ! Use dedicated routine to hedge for end point
                            ! intersections
                            call GGTMLineIntersections(ggtmdata, hflinek, hfline1, &
                                xint, yint, s1, s2, vertbased=.true.)

                            ! Check
                            if (size(xint) > 0) then
                                !$omp critical
                                keepind(k) = .false.
                                !$omp end critical
                            end if 
                        end if 

                        ! Last lfline
                        if (keepind(k)) then 
                            ! Use dedicated routine to hedge for end point
                            ! intersections
                            call GGTMLineIntersections(ggtmdata, hflinek, lflinen, &
                                xint, yint, s1, s2, vertbased=.true.)

                            ! Check
                            if (size(xint) > 0) then
                                !$omp critical
                                keepind(k) = .false.
                                !$omp end critical
                            end if 
                        end if 

                        ! Previous hfline
                        if (keepind(k)) then 
                            ! Use dedicated routine to hedge for end point
                            ! intersections
                            call GGTMLineIntersections(ggtmdata, hflinek, ct(k-1)%hfline, &
                                xint, yint, s1, s2, vertbased=.true.)

                            ! Check
                            if (size(xint) > 0) then
                                !$omp critical
                                keepind(k) = .false.
                                !$omp end critical
                            end if 
                        end if 

                    else
                        ! Standard case, check both lines

                        ! First hfline, this lfline
                        if (keepind(k)) then 
                            ! Use dedicated routine to hedge for end point
                            ! intersections
                            call GGTMLineIntersections(ggtmdata, lflinek, hfline1, &
                                xint, yint, s1, s2, vertbased=.true.)

                            ! Check
                            if (size(xint) > 0) then
                                !$omp critical
                                keepind(k) = .false.
                                !$omp end critical
                            end if 
                        end if 

                        ! Last lfline, this lfline
                        if (keepind(k)) then 
                            ! Use dedicated routine to hedge for end point
                            ! intersections
                            call GGTMLineIntersections(ggtmdata, lflinek, lflinen, &
                                xint, yint, s1, s2, vertbased=.true.)

                            ! Check
                            if (size(xint) > 0) then
                                !$omp critical
                                keepind(k) = .false.
                                !$omp end critical
                            end if 
                        end if 

                        ! Last lfline, this hfline
                        if (keepind(k)) then 
                            ! Use dedicated routine to hedge for end point
                            ! intersections
                            call GGTMLineIntersections(ggtmdata, hflinek, lflinen, &
                                xint, yint, s1, s2, vertbased=.true.)

                            ! Check
                            if (size(xint) > 0) then
                                !$omp critical
                                keepind(k) = .false.
                                !$omp end critical
                            end if 
                        end if 

                        ! First hfline, this hfline
                        if (keepind(k)) then 
                            ! Use dedicated routine to hedge for end point
                            ! intersections
                            call GGTMLineIntersections(ggtmdata, hflinek, hfline1, &
                                xint, yint, s1, s2, vertbased=.true.)

                            ! Check
                            if (size(xint) > 0) then
                                !$omp critical
                                keepind(k) = .false.
                                !$omp end critical
                            end if 
                        end if
                        
                        ! Previous lfline with this lfline
                        if (keepind(k)) then 
                            ! Use dedicated routine to hedge for end point
                            ! intersections
                            call GGTMLineIntersections(ggtmdata, lflinek, ct(k-1)%lfline, &
                                xint, yint, s1, s2, vertbased=.true.)

                            ! Check
                            if (size(xint) > 0) then
                                !$omp critical
                                keepind(k) = .false.
                                !$omp end critical
                            end if 
                        end if
                        
                        ! Next lfline with this lfline
                        if (keepind(k)) then 
                            ! Use dedicated routine to hedge for end point
                            ! intersections
                            call GGTMLineIntersections(ggtmdata, lflinek, ct(k+1)%lfline, &
                                xint, yint, s1, s2, vertbased=.true.)

                            ! Check
                            if (size(xint) > 0) then
                                !$omp critical
                                keepind(k) = .false.
                                !$omp end critical
                            end if 
                        end if

                        ! Previous hfline with this hfline
                        if (keepind(k)) then 
                            ! Use dedicated routine to hedge for end point
                            ! intersections
                            call GGTMLineIntersections(ggtmdata, hflinek, ct(k-1)%hfline, &
                                xint, yint, s1, s2, vertbased=.true.)

                            ! Check
                            if (size(xint) > 0) then
                                !$omp critical
                                keepind(k) = .false.
                                !$omp end critical
                            end if 
                        end if

                        ! Next hfline with this hfline
                        if (keepind(k)) then 
                            ! Use dedicated routine to hedge for end point
                            ! intersections
                            call GGTMLineIntersections(ggtmdata, hflinek, ct(k+1)%hfline, &
                                xint, yint, s1, s2, vertbased=.true.)

                            ! Check
                            if (size(xint) > 0) then
                                !$omp critical
                                keepind(k) = .false.
                                !$omp end critical
                            end if 
                        end if
                    end if 

                     
                    
                    ! Housekeeping
                    end associate
                end do
                !$omp end parallel do

                ! Housekeeping
                end associate

            end do 

            
            ! Check if any were deleted
            if (any(.not. keepind)) then 
                ! Issue message
                print *, 'Deleting tubes for cells: ', tc, 'since they ' // & 
                    'intersect with the cell boundaries'

                ! Set logical
                vertexwasdeleted = .true. 

                ! Check for exceptional cases
                if (size(keepind) == 1 .and. .not. keepind(1)) then 
                    ! Only one tube that should be deleted -> cell 
                    ! boundaries intersect
                    print *, 'Cells ', tc, 'have only one tube of ' // & 
                        'which the boundaries intersect, cannot delete. ' // & 
                        'Grid will have intersecting cells...'
                    
                    ! Set to true again to prevent deletion
                    keepind = .true. 
                end if 

                ! Sanity checks - tubes to be deleted should always come
                ! in pairs
                if (size(keepind) > 1) then 
                    ! Check first
                    if (.not. keepind(1)) then 
                        if (keepind(2)) then 
                            keepind(2) = .false.
                            !call gdErrorHandler('ConstructTopologicalMeshCellFluxTubes: ' // & 
                            !    'only one tube marked for deletion, unexpected')
                        end if 
                    end if 

                    ! Check last
                    if (.not. keepind(nft)) then 
                        if (keepind(nft-1)) then 
                            keepind(nft-1) = .false.
                            !call gdErrorHandler('ConstructTopologicalMeshCellFluxTubes: ' // & 
                            !    'only one tube marked for deletion, unexpected')
                        end if 
                    end if 

                    ! Check others
                    k = 2
                    do while (k < nft)
                        if (.not. keepind(k)) then 
                            if (keepind(k-1) .and. keepind(k+1)) then 
                                ! Write out the tubes
                                do j = 1, size(tc)
                                    associate(ct        => celldata(tc(j))%tubes)
                                        associate(&
                                        hfline1     => ct(k-1)%hfline,    &
                                        lfline1     => ct(k-1)%lfline,    &
                                        hflinen     => ct(k+1)%hfline,  &
                                        lflinen     => ct(k+1)%lfline,  &
                                        hflinek     => ct(k)%hfline,    &
                                        lflinek     => ct(k)%lfline     &
                                        )
                                    write (tstring, '(a2, i4)') 'l1', tc(j)
                                    call Write2DPolygonData(hfline1%xv, hfline1%yv, tstring)
                                    write (tstring, '(a2, i4)') 'l2', tc(j)
                                    call Write2DPolygonData(lfline1%xv, lfline1%yv, tstring)
                                    write (tstring, '(a2, i4)') 'l3', tc(j)
                                    call Write2DPolygonData(hflinek%xv, hflinek%yv, tstring)
                                    write (tstring, '(a2, i4)') 'l4', tc(j)
                                    call Write2DPolygonData(lflinek%xv, lflinek%yv, tstring)
                                    write (tstring, '(a2, i4)') 'l5', tc(j)
                                    call Write2DPolygonData(hflinen%xv, hflinen%yv, tstring)
                                    write (tstring, '(a2, i4)') 'l6', tc(j)
                                    call Write2DPolygonData(lflinen%xv, lflinen%yv, tstring)
                                    end associate
                                    end associate
                                end do

                                ! Print error
                                call gdErrorHandler('ConstructTopologicalMeshCellFluxTubes: ' // & 
                                'only one tube marked for deletion, unexpected')
                            end if 
                        end if
                        k = k + 1
                    end do
                end if 
                        
                ! Check if we still need to delete
                if (any(.not. keepind)) then 
                    ! Check which tubes to merge
                    t1 = 0
                    t2 = 0
                    do while (t1 < nct) 
                        ! Check
                        if (t2 == nct) then 
                            exit
                        end if 

                        ! Get the tube indices
                        t1 = findloc(keepind(t2+1:size(keepind)), .false., 1) + t2
                        if (t1 == t2) then 
                            ! No more start is found, exit
                            exit 
                        end if 

                        t2 = findloc(keepind(t1+1:size(keepind)), .true., 1) + t1
                        if (t2 == t1) then 
                            ! Should be last one to merge with
                            t2 = nct 
                        else
                            ! Need to subtract 1 since t2 is now a tube to be kept
                            t2 = t2 - 1
                        end if

                        ! Loop over each cell to adjust (first tube is 
                        ! adjusted)
                        do cc = 1, size(tc)
                            ! Take hfline and lfline
                            celldata(tc(cc))%tubes(t1)%lfline = celldata(tc(cc))%tubes(t2)%lfline
                        end do 

                        ! Adjust keepind for first tube (is adapted tube)
                        keepind(t1) = .true. 

                    end do 

                    ! Remove
                    do cc = 1, size(tc)
                        celldata(tc(cc))%tubes = pack(celldata(tc(cc))%tubes, keepind)
                        celldata(tc(cc))%lines = pack(celldata(tc(cc))%lines, keepind(2:size(keepind)))
                    end do 
                end if
            end if 

            ! Housekeeping
            deallocate(keepind)
        end do

        ! Update vertex numbering
        !========================
        ! Initialize
        allocate(isvertexdeleted(grid%vert%ntot))
        isvertexdeleted = .true. 
        do i = 1, cell%ntot 
            ! Check if vertices are present & ensure properly updated lines
            do j = 1, size(celldata(i)%tubes)
                call celldata(i)%tubes(j)%hfline%UpdateLineData(ggtmdata)
                call celldata(i)%tubes(j)%lfline%UpdateLineData(ggtmdata)
                isvertexdeleted(celldata(i)%tubes(j)%hfline%vert) = .false.
                isvertexdeleted(celldata(i)%tubes(j)%lfline%vert) = .false.
            end do
        end do 

        ! Get all IDs and construct mapping
        allIDs = pack([(k, k = 1, grid%vert%ntot)], .not. isvertexdeleted)
        allocate(vertmap(grid%vert%ntot))
        vertmap = 0_I8
        vertmap(allIDs) = [(k, k = 1, count(.not. isvertexdeleted))]

        ! Loop and adjust IDs - first segments, then tubes
        do i = 1, ggtmdata%nseg 
            ggtmdata%seg(i)%vert = vertmap(ggtmdata%seg(i)%vert)
            ggtmdata%seg(i)%sv = vertmap(ggtmdata%seg(i)%sv)
            ggtmdata%seg(i)%ev = vertmap(ggtmdata%seg(i)%ev)
        end do 
        do i = 1, cell%ntot 
            ! Update
            do j = 1, size(celldata(i)%tubes)
                call celldata(i)%tubes(j)%hfline%UpdateLineData(ggtmdata)
                call celldata(i)%tubes(j)%lfline%UpdateLineData(ggtmdata)
            end do 
        end do 

        ! Update number of grid vertices
        grid%vert%ntot = count(.not. isvertexdeleted)

        ! Add LOS
        !========
        ! call DetermineLOSlimits(ggtmdata)

        ! Construct graph
        !$omp parallel default(none) & 
        !$omp private(i, j) shared(ggtmdata)
        !$omp single 
        do i = 1, cell%ntot 
           do j = 1, size(celldata(i)%tubes)
                !$omp task
                call celldata(i)%tubes(j)%InitializeGraph()
                !$omp end task
            end do 
        end do 
        !$omp end single
        !$omp end parallel 
        
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
        integer(I8)                             :: i, j, k, cc

        ! Initialize
        !===========
        ! Associate
        associate(&
            seg             => ggtmdata%seg,    &
            cell            => topomesh%cell,   &
            celldata        => ggtmdata%cell    &
            )
        ! Allocate
        allocate(xv(grid%vert%ntot), yv(grid%vert%ntot), &
            isvertfound(grid%vert%ntot), fieldlineID(grid%vert%ntot))
        isvertfound = .false. 
        fieldlineID = 0 ! initialize

        ! Determine vertices
        !===================
        do i = 1, cell%ntot 
            ! Loop over tubes
            do j = 1, size(celldata(i)%tubes)
                ! Unpack for ease
                associate(tube      => celldata(i)%tubes(j))

                ! Loop over high field line segments
                do k = 1, tube%hfline%ns
                    ! Get IDs
                    tvID = seg(tube%hfline%segID(k))%vert

                    ! Set logicals
                    isvertfound(tvID) = .true.

                    ! Add coordinates
                    xv(tvID) = seg(tube%hfline%segID(k))%xv
                    yv(tvID) = seg(tube%hfline%segID(k))%yv
                    
                    ! Set field line ID (check)
                    do cc = 1, size(tvID)
                        if (fieldlineID(tvID(cc)) == 0) then 
                            fieldlineID(tvID(cc)) =   seg(tube%hfline%segID(k))%fsID
                        end if 
                    end do 
                end do 

                ! Loop over low field line segments
                do k = 1, tube%lfline%ns
                    ! Get IDs
                    tvID = seg(tube%lfline%segID(k))%vert

                    ! Set logicals
                    isvertfound(tvID) = .true.

                    ! Add coordinates
                    xv(tvID) = seg(tube%lfline%segID(k))%xv
                    yv(tvID) = seg(tube%lfline%segID(k))%yv
                    
                    ! Set field line ID
                    do cc = 1, size(tvID)
                        if (fieldlineID(tvID(cc)) == 0) then 
                            fieldlineID(tvID(cc)) =   seg(tube%lfline%segID(k))%fsID
                        end if 
                    end do 
                end do 

                ! Housekeeping
                end associate
            end do
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

        ! Note 4: the graph of the field line pair is used (in the default case) to 
        ! determine if vertex connections can be made. This graph should
        ! be computed beforehand when initializing the tubes of each 
        ! cell. It should only hold connections between hfline and 
        ! lfline vertices (so no lfline-lfline or hfline-hfline 
        ! connections) that do not lead to edges out of the tube or 
        ! edges intersecting the hfline or lfline. These edges should
        ! represent all valid edges. During face formation, we 'remove'
        ! nodes from the graph (not really) and see whether we still 
        ! have a connected graph. If this is not the case, the edge 
        ! cannot be chosen and is set to be illegal. This should lead to
        ! the best possible result, but is not sufficient to prevent
        ! cell overlap if the graph is disconnected from the start. In
        ! that case, cells will always overlap. This may be prevented
        ! in the future by doing checks beforehand and adding/deleting
        ! vertices where necessary. 

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMDataUDT)                      :: ggtmdata
        class(GGGridUDT), intent(inout)         :: grid 
        type(MagneticFieldUDT), intent(in)      :: magneticField

        ! Auxiliary
        real(R8)                                :: dx1, dy1, &
            alpha1, dx2, dy2, alpha2, dx3, dy3, alpha3, bxf(1:3), byf(1:3)
        real(R8), allocatable, dimension(:)     :: dx, dy, &
            dp, bx, by
        integer(I8)                             :: nf, nc, ncv, &
            v1, v3, v2, v4, tff(1:2), indmin, n1, n2, &
            vind, nk
        integer(I8), allocatable, dimension(:)  :: tempfacelabels, &
            tempcellvert, tempfaceregion, tempcellregion, tracingdir, &
            ii, jj
        integer(I8), allocatable, dimension(:, :)   :: tempfacevert, &
            tempcellvertP
        logical                                 :: &
            doquad, islegaltria1, islegaltria2, islegalquad, &
            dolasttriangle, isv2onl1
        logical, allocatable, dimension(:)      :: skipvert
        type(GGTMFieldlineDataUDT)              :: thisline

        ! Loop
        integer(I8)                             :: i, j, k, k1, k2

        ! Initialize
        !===========
        associate(&
            celldata    => ggtmdata%cell,   &
            seg         => ggtmdata%seg     &
            )

        ! Precompute
        allocate(tracingdir(size(celldata)))
        tracingdir = 0_I8
        do i = 1, size(celldata)
            ! Associate for ease
            associate(&
                tubes     => celldata(i)%tubes,     &
                doBLstart     => celldata(i)%linerefoptions%doBLstart,     &
                doBLend       => celldata(i)%linerefoptions%doBLend  ,     &
                ncBLstart     => celldata(i)%linerefoptions%ncBLstart,     &
                ncBLend       => celldata(i)%linerefoptions%ncBLend        &
                )

            ! Check tracing direction based on first line vertices
            if (.not. seg(tubes(1)%hfline%segID(1))%isvertex) then
                thisline = tubes(1)%hfline
            elseif (.not. seg(tubes(1)%lfline%segID(1))%isvertex) then 
                thisline = tubes(1)%lfline
            else 
                call gdErrorHandler('ConstructCellsQuadTria: could not find' // & 
                    'line that is longer than a single point, not supported')
            end if 
            dx = thisline%xv(2:thisline%nv) - thisline%xv(1:thisline%nv-1)
            dy = thisline%yv(2:thisline%nv) - thisline%yv(1:thisline%nv-1)
            allocate(bx(size(dx)), by(size(dy)))
            call magneticField%interp%Evaluate(thisline%xv(1:thisline%nv-1), thisline%yv(1:thisline%nv-1), 1, 0, by)
            call magneticField%interp%Evaluate(thisline%xv(1:thisline%nv-1), thisline%yv(1:thisline%nv-1), 0, 1, bx)
            bx = -bx
            dp = dx*bx + dy*by

            if (sum(dp) >= 0.0_R8) then 
                tracingdir(i) = 1_I8 
            else
                tracingdir(i) = -1_I8 
            end if 
            deallocate(bx, by)
            end associate
        end do 

        ! Pre-collapse double loop over celldata and tubes for 
        ! easier parallellization
        nk = 0
        do i = 1, size(celldata)
            nk = nk + size(celldata(i)%tubes)
        end do 
        allocate(ii(nk), jj(nk))
        k = 0
        do i = 1, size(celldata)
            do j = 1, size(celldata(i)%tubes)
                k = k + 1
                ii(k) = i
                jj(k) = j
            end do 
        end do 

        ! Loop
        !=====
        !$omp parallel default(none) &
        !$omp private(i, j, k, k1, k2, nf, &
        !$omp nc, ncv, n1, n2, tempfacevert, tempfacelabels, skipvert, &
        !$omp dolasttriangle, vind, v1, v2, v3, tempcellvert, tempcellvertP, &
        !$omp islegaltria1, islegaltria2, islegalquad, isv2onl1, doquad, &
        !$omp dx1, dx2, dx3, dy1, dy2, dy3, bxf, byf, alpha1, alpha2, &
        !$omp alpha3, indmin, v4, tff, tempcellregion, tempfaceregion) &
        !$omp shared(ggtmdata, grid, magneticField, tracingdir, ii, jj, nk)
        !$omp do schedule(dynamic)
        
        do k = 1, nk 
            ! Unpack iterates
            i = ii(k)
            j = jj(k)

            ! Associate for ease
            associate(&
                tubes     => celldata(i)%tubes,     &
                doBLstart     => celldata(i)%linerefoptions%doBLstart,     &
                doBLend       => celldata(i)%linerefoptions%doBLend  ,     &
                ncBLstart     => celldata(i)%linerefoptions%ncBLstart,     &
                ncBLend       => celldata(i)%linerefoptions%ncBLend        &
                )

            ! Unpack
            associate(&
                graph   => tubes(j)%graph,      &
                l1      => tubes(j)%hfline,     &
                l2      => tubes(j)%lfline,     &
                sff     => tubes(j)%srflabel,   &
                eff     => tubes(j)%erflabel    &
                )
            ! Initialize
            k1 = 1
            k2 = 1
            nf = 0
            nc = 0
            ncv = 0
            n1 = tubes(j)%hfline%nv 
            n2 = tubes(j)%lfline%nv
            allocate(tempfacevert(4*(n1+n2), 2), tempfacelabels(4*(n1+n2)), &
                tempcellvert(3*(n1+n2)), tempcellvertP(4*(n1+n2), 2)) ! overestimations
            allocate(skipvert(graph%nv))
            skipvert = .false.

            ! Compute face labels
            !$omp critical
            call l1%UpdateLineGriddingData(ggtmdata)
            call l2%UpdateLineGriddingData(ggtmdata)
            !$omp end critical

            ! Hedge for vertex lines
            if (seg(tubes(j)%hfline%segID(1))%isvertex) then
                n1 = 1
            end if 
            if (seg(tubes(j)%lfline%segID(1))%isvertex) then
                n2 = 1
            end if 

            ! Hedge for last vertices being the same
            dolasttriangle = .false. 
            if (l1%vert(n1) == l2%vert(n2)) then 
                ! Last cell should be treated as triangle IF we don't
                ! have a single triangle
                dolasttriangle = .true.

                ! Last vertex does not have to be connected in the 
                ! graph, so set to false
                vind = graph%GetVertexIndex(l1%vert(n1))
                skipvert(vind) = .true. 

                ! Make sure to skip in main loop
                n1 = n1 - 1
                n2 = n2 - 1

                ! Sanity check
                if (n1 == 0 .or. n2 == 0) then 
                    call gdErrorHandler('ConstructCellsQuadsTria: ' // & 
                        'first vertex is the same, yet appears to be ' // & 
                        'tangency point, unexpected')
                end if 
            end if 

            ! Hedge for first vertices being the same (insert triangle)
            if (l1%vert(k1) == l2%vert(k2)) then 
                ! Add the triangle by adding the first two faces (third
                ! face will be added later automatically)

                ! Set vertices
                v1 = l1%vert(k1) ! common vertex
                v2 = l1%vert(k1+1)
                v3 = l2%vert(k2+1)

                ! Set faces
                nf = nf + 1
                call tubes(j)%hfface%Append(nf) ! local index, updated later
                tempfacevert(nf, :) = [v1, v2]
                tempfacelabels(nf) = l1%facelabels(k1)
                nf = nf + 1
                call tubes(j)%lfface%Append(nf)
                tempfacevert(nf, :) = [v1, v3]
                tempfacelabels(nf) = l2%facelabels(k2)
                nf = nf + 1
                call tubes(j)%tubeface%Append(nf)
                tempfacevert(nf, :) = [v2, v3]
                tempfacelabels(nf) = 0

                ! Set cell
                nc = nc + 1
                call tubes(j)%cell%Append(nc) ! local index, updated later
                tempcellvert(ncv+1:ncv+3) = [v1, v2, v3]
                tempcellvertP(nc, :) = [ncv+1, 3]
                ncv = ncv + 3

                ! First vertex can be disconnected from the graph
                vind = graph%GetVertexIndex(l1%vert(k1)) 
                skipvert(vind) = .true. 

                ! Update k1, k2
                k1 = k1 + 1
                k2 = k2 + 1

                if (k1 > n1) then 
                    k1 = n1 
                end if 
                if (k2 > n2) then 
                    k2 = n2
                end if 

                ! Check if we should really add the last triangle
                if (dolasttriangle .and. (l1%nv == 2 .and. l2%nv == 3) &
                    .or. (l1%nv == 3 .and. l2%nv == 2)) then 
                    ! Start and end triangle are the same...
                    dolasttriangle = .false.
                end if 
                    
            
            else
                ! Add the first face 
                nf = nf + 1
                call tubes(j)%tubeface%Append(nf)
                tempfacevert(nf, :) = [l1%vert(k1), l2%vert(k2)]
                tempfacelabels(nf) = sff
            end if 

            ! Loop
            do while (k1 < n1 .or. k2 < n2)

                ! Make candidate faces
                !---------------------
                ! Face pair 1: vertex k1 and vertex k2+1, vertex k2, k2+1
                ! (triangle)
                
                ! Face pair 2: vertex k1+1 and vertex k2, vertex k1, k1+1
                ! (triangle)
                
                ! Face pair 3: k1+1, k2+1 (quad)

                ! Determine which face to take
                !-----------------------------
                ! Check which faces are legal
                call DetermineLegalCellsQuadTria(islegaltria1, &
                    islegaltria2, islegalquad, l1, l2, k1, k2, n1, n2, &
                    tracingdir(i), magneticField, celldata(i)%linerefoptions, &
                    celldata(i)%legalcellstyle, graph, skipvert)

                ! If none are legal, then throw warning for 
                ! overlapping cells and reset
                if (.not. any([islegaltria1, islegaltria2, islegalquad])) then 
                    ! We don't have a fix for this yet...
                    print *, 'ConstructCellsQuadTria: could not ' // & 
                        'find non-overlapping cell. Overlapping ' // &
                        'cells will be present in the grid...'
                    print *, 'cell: ', i, 'line: ', j, 'near vertex ID: ', &
                        l1%vert(k1), 'coordinates: ', l1%xv(k1), l1%yv(k1)
                    call tubes(j)%VisualizeGraph('lpgraph')

                    ! Reset to continue...
                    islegaltria1 = .true. 
                    islegaltria2 = .true. 
                    islegalquad = .true. 

                end if 

                ! First two vertices are always the same
                v1 = l1%vert(k1)
                v3 = l2%vert(k2)
                isv2onl1 = .false.
                
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
                


                    ! Check which triangles are allowed
                    if (.not. islegaltria1) then 
                        alpha1 = posinfval_R8()
                    end if 
                    if (.not. islegaltria2) then 
                        alpha2 = posinfval_R8()
                    end if 
                    if (.not. islegalquad) then 
                        alpha3 = posinfval_R8()
                    end if 
                    
                    ! Find face that makes the smallest angle
                    indmin = minloc([alpha1, alpha2, alpha3], 1)

                    ! Add the face
                    if (indmin == 3) then 
                        ! Add third face, quad
                        doquad = .true.
                        isv2onl1 = .false.
                        v2 = l2%vert(k2+1)
                        v4 = l1%vert(k1+1)
                        
                        ! Disconnect previous vertices from graph
                        vind = graph%GetVertexIndex(l1%vert(k1))
                        skipvert(vind) = .true.
                        vind = graph%GetVertexIndex(l2%vert(k2))
                        skipvert(vind) = .true.

                        ! Update counter
                        k2 = k2+1
                        k1 = k1+1
                    elseif (indmin == 1) then 
                        ! Add first face, triangle
                        v2 = l2%vert(k2+1)
                        isv2onl1 = .false.
                        tff = [0, l2%facelabels(k2)]

                        ! Disconnect previous vertex from graph
                        vind = graph%GetVertexIndex(l2%vert(k2))
                        skipvert(vind) = .true.
                        
                        ! Update counter
                        k2 = k2 + 1
                    elseif (indmin == 2 ) then 
                        ! Add second face, triangle
                        v2 = l1%vert(k1+1)
                        isv2onl1 = .true.
                        tff = [l1%facelabels(k1), 0]

                        ! Disconnect previous vertex from graph
                        vind = graph%GetVertexIndex(l1%vert(k1))
                        skipvert(vind) = .true.

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
                    isv2onl1 = .false.
                    tff = [0, l2%facelabels(k2)]

                    ! Disconnect previous vertex from graph
                    vind = graph%GetVertexIndex(l2%vert(k2))
                    skipvert(vind) = .true.
                    
                    ! Update counter
                    k2 = k2 + 1
                    
                    if (k2 == n2 .and. .not. dolasttriangle) then 
                        tff(1) = eff
                    end if 
                elseif ((k1 < n1) .and. (k2 == n2)) then 
                    ! We have to take the second option
                    ! Add second face pair
                    v2 = l1%vert(k1+1)
                    isv2onl1 = .true.
                    tff = [l1%facelabels(k1), 0]

                    ! Disconnect previous vertex from graph
                    vind = graph%GetVertexIndex(l1%vert(k1))
                    skipvert(vind) = .true.
                    
                    ! Update counter
                    k1 = k1 + 1
                    
                    if (k1 == n1 .and. .not. dolasttriangle) then  
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
                    call tubes(j)%hfface%Append(nf) 
                    tempfacevert(nf, :) = [v1, v4]
                    tempfacelabels(nf) = l1%facelabels(k1-1)
                    nf = nf + 1
                    call tubes(j)%lfface%Append(nf) 
                    tempfacevert(nf, :) = [v3, v2]
                    tempfacelabels(nf) = l2%facelabels(k2-1)
                    nf = nf + 1
                    call tubes(j)%tubeface%Append(nf) 
                    tempfacevert(nf, :) = [v4, v2]
                    if ((k1 /= n1) .or. (k2 /= n2) .or. dolasttriangle) then 
                        tempfacelabels(nf) = 0
                    else
                        tempfacelabels(nf) = eff
                    end if 
                else
                    nf = nf + 1
                    if (isv2onl1) then 
                        call tubes(j)%hfface%Append(nf) 
                    else
                        call tubes(j)%tubeface%Append(nf) 
                    end if 
                    tempfacevert(nf, :) = [v1, v2]
                    tempfacelabels(nf) = tff(1)
                    nf = nf + 1
                    if (isv2onl1) then 
                        call tubes(j)%tubeface%Append(nf) 
                    else
                        call tubes(j)%lfface%Append(nf) 
                    end if
                    tempfacevert(nf, :) = [v3, v2]
                    tempfacelabels(nf) = tff(2)
                end if 

                ! Add cell
                nc = nc + 1
                call tubes(j)%cell%Append(nc)
                if (doquad .and. v1 /= v3) then 
                    tempcellvert(ncv+1:ncv+4) = [v1, v3, v2, v4]
                    tempcellvertP(nc, :) = [ncv+1, 4]
                    ncv = ncv + 4
                elseif (doquad) then 
                    tempcellvert(ncv+1:ncv+3) = [v1, v2, v4]
                    tempcellvertP(nc, :) = [ncv+1, 3]
                    ncv = ncv + 2
                else
                    tempcellvert(ncv+1:ncv+3) = [v1, v2, v3]
                    tempcellvertP(nc, :) = [ncv+1, 3]
                    ncv = ncv + 3
                end if 
                
                ! Check stop criterium
                if ((k1 >= n1) .and. (k2 >= n2)) then 
                    !! Check if the last non-aligned face was equal to the first one
                    !issameface = (tempfacevert(nf-1, 1) == tempfacevert(1, 1) ) .and. &
                    !    (tempfacevert(nf-1, 2) == tempfacevert(1, 2))
                    !issameface = issameface .or. (tempfacevert(nf-1, 2) == tempfacevert(1, 1) ) .and. &
                    !    (tempfacevert(nf-1, 1) == tempfacevert(1, 2))
                    !if (issameface) then 
                    !    ! Don't add the last face
                    !    nf = nf-1
                    !end if 
                    exit
                end if 

                

            end do

            ! Check for last cell
            if (dolasttriangle) then 
                ! Add the triangle by adding the last two faces (third
                ! face already added normally speaking)

                ! Set vertices
                k1 = n1 ! reset to be sure
                k2 = n2
                v1 = l1%vert(k1+1) ! common vertex
                v2 = l1%vert(k1) ! previous vertices
                v3 = l2%vert(k2)

                ! Set faces
                nf = nf + 1
                call tubes(j)%hfface%Append(nf)
                tempfacevert(nf, :) = [v1, v2]
                tempfacelabels(nf) = l1%facelabels(k1)
                nf = nf + 1
                call tubes(j)%lfface%Append(nf)
                tempfacevert(nf, :) = [v1, v3]
                tempfacelabels(nf) = l2%facelabels(k2)

                ! Set cell
                nc = nc + 1
                call tubes(j)%cell%Append(nc)
                tempcellvert(ncv+1:ncv+3) = [v1, v2, v3]
                tempcellvertP(nc, :) = [ncv+1, 3]
                ncv = ncv + 3
            end if 

            ! Add to grid
            allocate(tempcellregion(nc), tempfaceregion(nf))
            tempcellregion = celldata(i)%region
            tempfaceregion = celldata(i)%region
            !$omp critical
            tubes(j)%cell       = tubes(j)%cell + grid%cell%ntot 
            tubes(j)%tubeface   = tubes(j)%tubeface + grid%face%ntot 
            tubes(j)%lfface     = tubes(j)%lfface + grid%face%ntot 
            tubes(j)%hfface     = tubes(j)%hfface + grid%face%ntot
            call grid%AddFace(tempfacevert(1:nf, :), tempfacelabels(1:nf), tempfaceregion(1:nf))
            call grid%AddCell(tempcellvert(1:ncv), tempcellvertP(1:nc, 1:2), tempcellregion(1:nc))
            !$omp end critical

            ! Housekeeping
            deallocate(tempfacevert, tempfacelabels, tempcellvert, &
                tempcellvertP, tempcellregion, tempfaceregion, &
                skipvert)
            end associate

            ! Housekeeping
            end associate
        end do
        !$omp end do
        !$omp end parallel

        ! Cleanup
        !========
        call RemoveDuplicateGridFaces(ggtmdata, grid)

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Auxiliary function for ConstructCellsQuadTria to determine legal
    ! triangles/quads
    subroutine DetermineLegalCellsQuadTria(islegaltria1, islegaltria2, &
        islegalquad, l1, l2, k1, k2, n1, n2, tracingdir, magneticField, &
        refoptions, legalcellstyle, graph, skipvert)

        ! Description
        !============
        ! This is an auxiliary subroutine that should only be used by 
        ! constructors of quads/triangles to determine if a triangle or
        ! quad is legal vor construction. For this, the current two 
        ! gridding lines must be given in l1, l2, and the current node
        ! positions in k1, k2 (resp.). Refinement options shall be 
        ! parsed through the refoptions structure. The routine accounts
        ! for specifications in this refoptions struct, but may override
        ! certain desired features if e.g. overlapping cells would be 
        ! inserted this way. Do note that in no way this routine can 
        ! guarantee non-overlapping cells, although it will catch most 
        ! of them (typically not if at the start of the face the cells
        ! already overlap due to flux surface/bnd shape). 

        ! Additional required inputs are:
        ! - tracingdir:     1 if we go along the magnetic field in 
        !                   increasing k1, k2 direction, -1 otherwise
        ! - magneticField:  magnetic field object to evaluate 
        !                   magnetic field vector (used to deambigue 
        !                   non-convex cases)

        ! Candidate faces are:
        ! Face pair 1: vertex k1 and vertex k2+1, vertex k2, k2+1
        ! (triangle 1)
        
        ! Face pair 2: vertex k1+1 and vertex k2, vertex k1, k1+1
        ! (triangle 2)
        
        ! Face pair 3: k1+1, k2+1 (quad)

        ! Notes
        !======
        ! Note 1: it is assumed that no intersections exist between 
        ! the lines and any starting/ending radial faces. If this is not
        ! the case, output may be unexpected.

        ! Note 2: this routine does *not* hedge for vertices that are 
        ! the same! This should only occur at the start/end of lines and
        ! should be hedged for upstream. 

        ! Note 3: skipvert is adjusted in this routine but should be 
        ! returned unaltered eventually (perhaps making a copy in this
        ! routine would be better at some point to avoid mistakes...)
        

        ! Declare variables
        !==================
        ! Arguments
        type(GGTMFieldlineDataUDT), intent(in)          :: l1, l2 
        integer(I8), intent(in)                         :: k1, k2, tracingdir, &
            n1, n2
        type(GGTMFieldlineRefinementOptionsUDT), intent(in)     :: refoptions
        logical, intent(out)                            :: islegaltria1, &
            islegaltria2, islegalquad
        type(MagneticFieldUDT), intent(in)              :: magneticField 
        class(UGraphUDT), intent(in)                    :: graph
        logical, dimension(:), intent(inout)            :: skipvert 
        character(*), intent(in)                        :: legalcellstyle

        ! Auxiliary
        real(R8)                        :: dx(1:7), dy(1:7), beta(1:4), &
            bx(1:3), by(1:3), dp(1:2), xs, ys
        real(R8), allocatable, dimension(:)     :: xint, yint
        integer(I8)                             :: vindk1, vindk2
        integer(I8), allocatable, dimension(:)  :: sint
        logical                                 :: oldstyle, doBLquad, &
            doBLendcheck, doBLtria1, doBLtria2, isconnected
        logical, allocatable, dimension(:)      :: keepind


        ! Initialize
        !===========
        oldstyle = .false. 
        ! Associate
        associate(&
            doBLstart     => refoptions%doBLstart,     &
            doBLend       => refoptions%doBLend  ,     &
            ncBLstart     => refoptions%ncBLstart,     &
            ncBLend       => refoptions%ncBLend        &
            )

        ! Set initial logical
        doBLquad        = .false. ! true if we should do a boundary layer quad
        doBLendcheck    = .true.
        islegaltria1    = .true.
        islegaltria2    = .true.
        islegalquad     = .true.

        ! Boundary layer checks
        !----------------------
        ! Check if we can/should insert a boundary layer
        ! of quads -> should be given precedence 
        doBLquad = .false.
        doBLtria1 = .false.
        doBLtria2 = .false.
        if (k1 < n1 .and. k2 < n2) then 
            if (doBLstart) then 
                if ((k1 < ncBLstart+1) .and. (k2 < ncBLstart+1)) then 
                    ! Ensure quad 
                    doBLquad = .true. 
                    doBLendcheck = .false. 
                end if 
            end if
            if (doBLend .and. doBLendcheck) then 
                ! Check that we're not marching further than
                ! allowed (not an issue for the initial 
                ! BL since we start there...)
                
                if (k1 > (n1-ncBLend-1) .and. .not. (k2 > (n2-ncBLend-1))) then 
                    ! Need to add triangle from k2:k2+1 if possible
                    doBLtria1 = .true.
                elseif (.not. (k1 > (n1-ncBLend-1)) .and. k2 > (n2-ncBLend-1)) then
                    ! Need to add triangle from k1:k1+1 if possible
                    doBLtria2 = .true.
                elseif ((k1 > (n1-ncBLend-1)) .and. (k2 > (n2-ncBLend-1))) then 
                    ! Ensure quad if possible
                    doBLquad = .true. 
                end if 
            end if
        end if 
        

        ! Checks
        !=======
        select case (legalcellstyle)

        case ('no')

            ! No checks at all

        case ('old')

            ! Old shitty method 

            ! Precompute data
            !----------------
            ! Compute vectors:
            ! 1: k1 -> k2+1
            ! 2: k2 -> k1+1
            ! 3: k1+1 -> k2+1
            ! 4: k1 -> k2
            ! 5: k1 -> k1+1
            ! 6: k2 -> k2+1
            ! 7: normal vector of 4, in tracing direction
            ! b: magnetic field vector in face center of 4
            
            if (k2 < n2) then 
                dx(1) = l2%xv(k2+1) - l1%xv(k1)
                dy(1) = l2%yv(k2+1) - l1%yv(k1)
                dx(6) = l2%xv(k2+1) - l2%xv(k2)
                dy(6) = l2%yv(k2+1) - l2%yv(k2)
            end if 
            if (k1 < n1) then 
                dx(2) = l1%xv(k1+1) - l2%xv(k2)
                dy(2) = l1%yv(k1+1) - l2%yv(k2)
                dx(5) = l1%xv(k1+1) - l1%xv(k1)
                dy(5) = l1%yv(k1+1) - l1%yv(k1)
            end if 
            if ((k1 < n1) .and. (k2 < n2)) then 
                dx(3) = l2%xv(k2+1) - l1%xv(k1+1)
                dy(3) = l2%yv(k2+1) - l1%yv(k1+1)
            end if 
            dx(4) = l2%xv(k2) - l1%xv(k1)
            dy(4) = l2%yv(k2) - l1%yv(k1)
            dx(7) = -dy(4) ! just initial value
            dy(7) = dx(4) ! just initial value

            call magneticField%interp%Evaluate([l1%xv(k1), l2%xv(k2)], &
                [l1%yv(k1), l2%yv(k2)], 1, 0, by(2:3))
            call magneticField%interp%Evaluate([l1%xv(k1), l2%xv(k2)], &
                [l1%yv(k1), l2%yv(k2)], 0, 1, bx(2:3))
            bx = -bx
            bx(1) = 0.5*(bx(2) + bx(3))
            by(1) = 0.5*(by(2) + by(3))

            if ((dx(7)*bx(1) + dy(7)*by(1)) >= 0.0_R8) then 
                ! normal is currently along magnetic field direction
                if (tracingdir < 0_I8) then 
                    ! Need to switch, we're tracing opposite to mf direction
                    dx(7) = -dx(7)
                    dy(7) = -dy(7)
                end if 
            else
                ! normal is currently opposite to magnetic field direction
                if (tracingdir > 0_I8) then 
                    ! Need to switch, we're tracing along mf direction
                    dx(7) = -dx(7)
                    dy(7) = -dy(7)
                end if 
            end if

            ! Compute angles:
            ! beta1: 5 -> 4
            ! beta2: 5 -> 1
            ! beta3: 6 -> -4 (so 4 but opposite direction)
            ! beta4: 6 -> 2
            if ((k1 < n1) .and. (k2 < n2)) then 
                beta(1) = atan2(dx(5)*dy(4) - dy(5)*dx(4), dx(5)*dx(4) + dy(5)*dy(4))
                beta(2) = atan2(dx(5)*dy(1) - dy(5)*dx(1), dx(5)*dx(1) + dy(5)*dy(1))
                beta(3) = atan2(-dx(6)*dy(4) + dy(6)*dx(4), -dx(6)*dx(4) - dy(6)*dy(4))
                beta(4) = atan2(dx(6)*dy(2) - dy(6)*dx(2), dx(6)*dx(2) + dy(6)*dy(2))
            end if 

            ! Compute dot products
            ! dp1: 5, 7
            ! dp2: 6, 7
            if (k1 < n1) then 
                dp(1) = dx(5)*dx(7) + dy(5)*dy(7)
            end if 
            if (k2 < n2) then 
                dp(2) = dx(6)*dx(7) + dy(6)*dy(7)
            end if 

            ! Ensure that angles are computed inside the tube by checking
            ! dot products
            if ((k1 < n1) .and. (k2 < n2)) then 
                if (dp(1) > 0 .and. dp(2) > 0) then 
                    ! Do nothing
                elseif (dp(1) < 0 .and. dp(2) > 0) then 
                    ! We know that l2 is properly oriented. The sign of 
                    ! the angle beta3 should be opposite to the sign of
                    ! beta1 - if not, add/subtract 2*pi
                    if (beta(3) > 0 .and. beta(1) > 0) then 
                        beta(1) = beta(1) - 2*pi_R8
                        if (beta(2) > 0) then 
                            beta(2) = beta(2) - 2*pi_R8
                        end if 
                    elseif (beta(3) < 0 .and. beta(1) < 0) then 
                        beta(1) = beta(1) + 2*pi_R8
                        if (beta(2) > 0) then 
                            beta(2) = beta(2) + 2*pi_R8
                        end if 
                    end if 
                elseif (dp(2) < 0 .and. dp(1) > 0) then 
                    ! We know that l1 is properly oriented. The sign of 
                    ! the angle beta3 should be opposite to the sign of
                    ! beta1 - if not, add/subtract 2*pi
                    if (beta(3) > 0 .and. beta(1) > 0) then 
                        beta(3) = beta(3) - 2*pi_R8
                        if (beta(4) > 0) then 
                            beta(4) = beta(4) - 2*pi_R8
                        end if 
                    elseif (beta(3) < 0 .and. beta(1) < 0) then 
                        beta(3) = beta(3) + 2*pi_R8
                        if (beta(4) > 0) then 
                            beta(4) = beta(4) + 2*pi_R8
                        end if 
                    end if 
                else
                    call Write2DPolygonData(l1%xv, l1%yv, 'l1')
                    call Write2DPolygonData(l2%xv, l2%yv, 'l2')
                    call WriteVertexData([k1, k2], [l1%xv(k1), l2%xv(k2)], &
                        [l1%xv(k1), l2%xv(k2)], 'kdata')
                    ! This shouldn't happen, both dot products negative
                    print *, 'DetermineLegalCellsQuadTria: both ' // & 
                        'dot products are negative, check if tracing direction ' // & 
                        'was correctly computed. Returning...'
                    return
                end if 
            end if 

            ! Checks
            !-------
            ! Are we at the end of one of the lines?
            if ((k1 < n1) .and. (k2 < n2)) then ! no
 
                ! Compute checks
                !---------------
                ! triangle 1: first check is on convexity. If dp2 < 0, then 
                !   the triangle will lie outside of the tube and is 
                !   therefore not legal. If dp2 > 0, then beta2 > beta1 if 
                !   beta2 > 0 to be legal, if beta1 < 0, then beta2 < beta1.
                ! triangle 2: similar to triangle 1, but now dp2 -> dp1, 
                !   beta1 -> beta3, beta2 -> beta4
                ! quad: if either dp1 or dp2 is smaller than zero (both is 
                !   unexpected), we don't allow quads (would be non-convex)
                !   unless explicitly forced by boundary layer and if other
                !   checks pass. If the quad face intersects with the 
                !   previous face (vector 4), then it is illegal to make 
                !   a quad, even with boundary layers. 

                ! Triangle 1
                if (dp(2) <= 0) then 
                    islegaltria1 = .false.
                else
                    if (beta(1) > 0) then 
                        if (beta(2) > beta(1)) then
                            islegaltria1 = .false.
                        end if 
                    else
                        if (beta(2) < beta(1)) then 
                            islegaltria1 = .false.
                        end if 
                    end if  
                end if  

                ! Triangle 2
                if (dp(1) <= 0) then 
                    islegaltria2 = .false.
                else
                    if (beta(3) > 0) then 
                        if (beta(4) > beta(3)) then
                            islegaltria2 = .false.
                        end if 
                    else
                        if (beta(4) < beta(3)) then 
                            islegaltria2 = .false.
                        end if 
                    end if 
                end if  

                ! Quad
                if (dp(1) < 0 .or. dp(2) < 0) then 
                    if (doBLquad) then 
                        ! Check intersection with previous edge
                        call SegmentIntersections(xs, ys, l1%xv(k1), l1%yv(k1), &
                            l2%xv(k2), l2%yv(k2), l1%xv(k1+1), l1%yv(k1+1), &
                            l2%xv(k2+1), l2%yv(k2+1))

                        if (.not. isnan(xs)) then 
                            ! Intersection found, set to false
                            islegalquad = .false. 
                        end if 
                    else
                        ! Default set to false
                        islegalquad = .false. 
                    end if
                    
                    ! Always legal otherwise
                end if    
                
            elseif ((k1 == n1) .and. (k2 < n2)) then 
                ! Only triangle 1 can be formed
                islegaltria2    = .false.
                islegalquad     = .false.

                ! Still need to check if this doesn't lead to overlap etc
                ! - can be used upstream for warnings etc
                if (dp(2) <= 0) then 
                    islegaltria1 = .false.
                end if  

            elseif ((k1 < n1) .and. (k2 == n2)) then 
                ! Only triangle 2 can be formed
                islegaltria1    = .false.
                islegalquad     = .false.

                ! Still need to check if this doesn't lead to overlap etc
                ! - can be used upstream for warnings etc
                if (dp(1) <= 0) then 
                    islegaltria2 = .false.
                end if  
                
            else
                ! This shouldn't happen
                call gdErrorHandler('Something wrong in quad gridder')
            end if 

            ! Compute non-local checks
            !-------------------------
            ! Here, we see if the faces that are formed do not
            ! intersect with other (non-tangent) faces. We only perform 
            ! these checks for faces that are still legal

            ! Check for intersections with remaining part of curve 1
            if (k1 < n1-1) then 
                ! First triangle
                if (k2 < n2) then 
                    if (islegaltria1) then 
                        ! Compute intersections
                        call SegmentSimplePolygonIntersections(l1%xv(k1+1:l1%nv), &
                                l1%yv(k1+1:l1%nv), l1%xv(k1), l1%yv(k1), &
                                l2%xv(k2+1), l2%yv(k2+1), xint, yint, sint)

                        ! Eliminate any intersections with edges that have
                        ! a vertex in common
                        if (size(xint) > 0) then 
                            allocate(keepind(size(sint)))
                            keepind = .true.
                            where (l1%vert(k1+sint) == l1%vert(k1) .or. &
                                l1%vert(k1+sint + 1) == l1%vert(k1))
                                keepind = .false. 
                            end where 
                            xint = pack(xint, keepind)
                            deallocate(keepind)
                        end if 

                        ! Check
                        if (size(xint) > 0) then 
                            ! Eliminate edges with the same vertices
                            islegaltria1 = .false. 
                        end if 
                    end if
                end if 

                ! Second triangle
                if (islegaltria2) then 
                    call SegmentSimplePolygonIntersections(l1%xv(k1+2:l1%nv), &
                        l1%yv(k1+2:l1%nv), l1%xv(k1+1), l1%yv(k1+1), &
                        l2%xv(k2), l2%yv(k2), xint, yint, sint)

                    ! Eliminate any intersections with edges that have
                    ! a vertex in common
                    if (size(xint) > 0) then 
                        allocate(keepind(size(sint)))
                        keepind = .true.
                        where (l1%vert(k1+sint+1) == l1%vert(k1+1) .or. &
                            l1%vert(k1+sint+2) == l1%vert(k1+1))
                            keepind = .false. 
                        end where 
                        xint = pack(xint, keepind)
                        deallocate(keepind)
                    end if 

                    ! Check
                    if (size(xint) > 0) then 
                        islegaltria2 = .false. 
                    end if 
                end if

                ! Quad
                if (k2 < n2) then 
                    ! Check intersections of quad face with other faces
                    if (islegalquad) then 
                        call SegmentSimplePolygonIntersections(l1%xv(k1+2:l1%nv), &
                            l1%yv(k1+2:l1%nv), l1%xv(k1+1), l1%yv(k1+1), &
                            l2%xv(k2+1), l2%yv(k2+1), xint, yint, sint)

                        ! Eliminate any intersections with edges that have
                        ! a vertex in common
                        if (size(xint) > 0) then 
                            allocate(keepind(size(sint)))
                            keepind = .true.
                            where (l1%vert(k1+sint+1) == l1%vert(k1+1) .or. &
                                l1%vert(k1+sint+2) == l1%vert(k1+1))
                                keepind = .false. 
                            end where 
                            xint = pack(xint, keepind)
                            deallocate(keepind)
                        end if 

                        if (size(xint) > 0) then 
                            islegalquad = .false. 
                            islegaltria1 = .false.
                        end if 
                    end if

                    ! Check intersection of next tria 1 with L1
                    if (k2 < n2-1) then 
                        if (islegalquad) then 
                            ! Compute intersections
                            call SegmentSimplePolygonIntersections(l1%xv(k1+2:l1%nv), &
                                    l1%yv(k1+2:l1%nv), l1%xv(k1+1), l1%yv(k1+1), &
                                    l2%xv(k2+2), l2%yv(k2+2), xint, yint, sint)
        
                            ! Eliminate any intersections with edges that have
                            ! a vertex in common
                            if (size(xint) > 0) then 
                                allocate(keepind(size(sint)))
                                keepind = .true.
                                where (l1%vert(k1+sint+1) == l1%vert(k1+1) .or. &
                                    l1%vert(k1+sint + 2) == l1%vert(k1+1))
                                    keepind = .false. 
                                end where 
                                xint = pack(xint, keepind)
                                deallocate(keepind)
                            end if 
        
                            ! Check
                            if (size(xint) > 0) then 
                                islegalquad = .false. 
                            end if 
                        end if
                    end if 
                end if 
            end if 

            ! Check for intersections with remaining part of curve 2
            if (k2 < n2-1) then 
                ! First triangle
                if (islegaltria1) then 
                    call SegmentSimplePolygonIntersections(l2%xv(k2+2:l2%nv), &
                        l2%yv(k2+2:l2%nv), l1%xv(k1), l1%yv(k1), &
                        l2%xv(k2+1), l2%yv(k2+1), xint, yint, sint)

                    ! Eliminate any intersections with edges that have
                    ! a vertex in common
                    if (size(xint) > 0) then 
                        allocate(keepind(size(sint)))
                        keepind = .true.
                        where (l2%vert(k2+sint+1) == l2%vert(k2+1) .or. &
                            l2%vert(k2+sint+2) == l2%vert(k2+1))
                            keepind = .false. 
                        end where 
                        xint = pack(xint, keepind)
                        deallocate(keepind)
                    end if 

                    ! Check
                    if (size(xint) > 0) then 
                        islegaltria1 = .false. 
                    end if 
                end if

                ! Second triangle
                if (k1 < n1) then 
                    if (islegaltria2) then 
                        call SegmentSimplePolygonIntersections(l2%xv(k2+1:l2%nv), &
                            l2%yv(k2+1:l2%nv), l1%xv(k1+1), l1%yv(k1+1), &
                            l2%xv(k2), l2%yv(k2), xint, yint, sint)

                        ! Eliminate any intersections with edges that have
                        ! a vertex in common
                        if (size(xint) > 0) then 
                            allocate(keepind(size(sint)))
                            keepind = .true.
                            where (l2%vert(k2+sint) == l2%vert(k2) .or. &
                                l2%vert(k2+sint+1) == l2%vert(k2))
                                keepind = .false. 
                            end where 
                            xint = pack(xint, keepind)
                            deallocate(keepind)
                        end if 

                        ! Check
                        if (size(xint) > 0) then 
                            islegaltria2 = .false. 
                        end if 
                    end if
                end if 

                ! Quad
                if (k1 < n1) then 

                    ! Check for quad face intersections
                    if (islegalquad) then 
                        call SegmentSimplePolygonIntersections(l2%xv(k2+2:l2%nv), &
                            l2%yv(k2+2:l2%nv), l1%xv(k1+1), l1%yv(k1+1), &
                            l2%xv(k2+1), l2%yv(k2+1), xint, yint, sint)

                        !Eliminate any intersections with edges that have
                        ! a vertex in common
                        if (size(xint) > 0) then 
                            allocate(keepind(size(sint)))
                            keepind = .true.
                            where (l2%vert(k2+sint+1) == l2%vert(k2+1) .or. &
                                l2%vert(k2+sint+2) == l2%vert(k2+1))
                                keepind = .false. 
                            end where 
                            xint = pack(xint, keepind)
                            deallocate(keepind)
                        end if 
        
                        ! Check
                        if (size(xint) > 0) then 
                            islegalquad = .false. 
                            islegaltria2 = .false.
                        end if 
                    end if

                    ! Check if next triangle 2 intersects
                    if (k1 < n1-1) then 
                        if (islegalquad) then 
                            call SegmentSimplePolygonIntersections(l2%xv(k2+2:l2%nv), &
                                l2%yv(k2+2:l2%nv), l1%xv(k1+2), l1%yv(k1+2), &
                                l2%xv(k2+1), l2%yv(k2+1), xint, yint, sint)
        
                            ! Eliminate any intersections with edges that have
                            ! a vertex in common
                            if (size(xint) > 0) then 
                                allocate(keepind(size(sint)))
                                keepind = .true.
                                where (l2%vert(k2+sint+1) == l2%vert(k2+1) .or. &
                                    l2%vert(k2+sint+2) == l2%vert(k2+1))
                                    keepind = .false. 
                                end where 
                                xint = pack(xint, keepind)
                                deallocate(keepind)
                            end if 
        
                            ! Check
                            if (size(xint) > 0) then 
                                islegalquad = .false. 
                            end if 
                        end if
                    end if 
                end if 
            end if 
        
        case default

            ! New style, graph based
            ! Initialize
            vindk1 = graph%GetVertexIndex(l1%vert(k1))
            vindk2 = graph%GetVertexIndex(l2%vert(k2))

            ! Test triangle 1
            !----------------
            if (k2 < n2) then 
                ! Remove vertex at k2 temporarily
                skipvert(vindk2) = .true. 

                ! Check if graph is still connected, if yes, then triangle 
                ! is legal
                isconnected = graph%IsConnected(skipvert)
                if (isconnected) then 
                    islegaltria1 = .true.
                else
                    islegaltria1 = .false.
                end if 

                ! Add vertex again
                skipvert(vindk2) = .false. 
            else
                islegaltria1 = .false.
            end if 

            ! Test triangle 2
            !----------------
            if (k1 < n1) then 
                ! Remove vertex at k1 temporarily
                skipvert(vindk1) = .true. 

                ! Check if graph is still connected, if yes, then triangle 
                ! is legal
                isconnected = graph%IsConnected(skipvert)
                if (isconnected) then 
                    islegaltria2 = .true.
                else
                    islegaltria2 = .false.
                end if 

                ! Add vertex again
                skipvert(vindk1) = .false.
            else 
                islegaltria2 = .false.
            end if 
            
            ! Test quad
            !----------
            ! Only if both triangles are legal we need to test
            if (islegaltria1 .and. islegaltria2) then 
                ! Remove both vertices temporarily
                skipvert(vindk1) = .true.
                skipvert(vindk2) = .true.
                
                ! Check if graph is still connected, if yes, then quad 
                ! is legal
                isconnected = graph%IsConnected(skipvert)
                if (isconnected) then 
                    islegalquad = .true.
                else
                    islegalquad = .false.
                end if 

                ! Add vertices again
                skipvert(vindk1) = .false.
                skipvert(vindk2) = .false.

            else
                islegalquad = .false.
            end if     
        
        end select
        
        ! Final checks
        !=============
        ! Override triangle legality if boundary layer is set to true
        ! and quad is legal (or if first/last vertices have the same ID)
        if (islegalquad .and. doBLquad) then 
            islegaltria1 = .false. 
            islegaltria2 = .false.
        elseif (islegaltria1 .and. doBLtria1) then 
            islegalquad = .false.
            islegaltria2 = .false.
        elseif (islegaltria2 .and. doBLtria2) then 
            islegalquad = .false. 
            islegaltria1 = .false.
        end if 
                
        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Line-of-sight determination for line pairs to determine possible
    ! cell shapes etc
    subroutine DetermineLOSlimits(ggtmdata)

        ! Description
        !============
        ! This routine determines the line-of-sight for each line pair 
        ! of each cell, which can be used during cell generation to 
        ! determine whether a cell can be formed or not, or during 
        ! grid construction to determine if more local refinement is
        ! needed. The line-of-sight limits are determined between 
        ! two lines, l1, l2, as follows:
        !
        ! - minLOS for l1, vertex k1: the vertex index k2, starting from 
        !   the first line, which forms the first face (k1, k2) that does
        !   not intersect either l1 or l2 (i.e. k1, k2-1 intersects). 
        ! - maxLOS for l1, vertex k1: similar to minLOS, but now the 
        !   vertex pair (k1, k2) is the last one, starting from minLOS, 
        !   that does not intersect with l1 or l2 (i.e. k1, k2+1 interects)
        ! - similar definitions hold for l2 of course
        ! - for each of these edges, there's the additional requirement
        !   that it has to be facing to the interior of the domain 
        !   determined by the line pairs. 
        
        ! Notes
        !======
        ! Note 1: it is assumed that all line data is initialized and 
        ! has proper (initial) dimensions. 

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMDataUDT)                      :: ggtmdata 

        ! Auxiliary
        real(R8), allocatable, dimension(:)     :: xint, yint, xint1, &
            xint2
        integer(I8), allocatable, dimension(:)  :: sint
        type(GGTMFieldlineDataUDT)              :: l1, l2

        ! Loop  
        integer(I8)                             :: i, j, k1, k2

        ! Initialize
        !===========
        allocate(xint(0), xint1(0), xint2(0))

        ! Determine LOS
        !==============
        do i = 1, size(ggtmdata%cell)
            ! Associate
            associate(tc        => ggtmdata%cell(i))

            ! Loop over all line pairs
            do j = 1, size(tc%tubes)
                ! Initialize
                tc%tubes(j)%l1maxLOS = spread(0_I8, 1, tc%tubes(j)%hfline%nv)
                tc%tubes(j)%l1minLOS = spread(0_I8, 1, tc%tubes(j)%hfline%nv)
                tc%tubes(j)%l2maxLOS = spread(0_I8, 1, tc%tubes(j)%lfline%nv)
                tc%tubes(j)%l2minLOS = spread(0_I8, 1, tc%tubes(j)%lfline%nv)

                ! Unpack
                l1 = tc%tubes(j)%hfline
                l2 = tc%tubes(j)%lfline

                ! Determine LOS for l1
                do k1 = 1, l1%nv 
                    ! minLOS
                    k2 = 1
                    
                    do while (k2 <= l2%nv)
                        ! Compute segment intersections, l1
                        if (k1 > 1) then 
                            ! Check for intersections in l1(1:k1-1)
                            call SegmentSimplePolygonIntersections(l1%xv(1:k1-1), &
                                l1%yv(1:k1-1), l1%xv(k1), l1%yv(k1), l2%xv(k2), l2%yv(k2), &
                                xint1, yint, sint)
                        end if 
                        if (k1 < l1%nv) then 
                            ! Check for intersections in l1(k1+1:nv)
                            call SegmentSimplePolygonIntersections(l1%xv(k1+1:l1%nv), &
                                l1%yv(k1+1:l1%nv), l1%xv(k1), l1%yv(k1), l2%xv(k2), l2%yv(k2), &
                                xint2, yint, sint)
                        end if 

                        ! Check
                        if (size(xint1) == 0 .and. size(xint2) == 0) then 
                            exit
                        else 
                            ! Update k2
                            k2 = k2 + 1
                        end if
                    end do 
                    

                    ! Set the minimal value
                    tc%tubes(j)%l1minLOS(k1) = k2

                    ! maxLOS
                    k2 = tc%tubes(j)%l1minLOS(k1) ! start from minLOS
                    do while (k2 <= l2%nv)
                        ! Compute segment intersections 
                        if (k1 > 1) then 
                            ! Check for intersections in l1(1:k1-1)
                            call SegmentSimplePolygonIntersections(l1%xv(1:k1-1), &
                                l1%yv(1:k1-1), l1%xv(k1), l1%yv(k1), l2%xv(k2), l2%yv(k2), &
                                xint, yint, sint)
                            
                            if (size(xint) > 0) then 
                                ! Exit the loop
                                exit
                            end if  
                        end if 
                        if (k1 < l1%nv) then 
                            ! Check for intersections in l1(k1+1:nv)
                            call SegmentSimplePolygonIntersections(l1%xv(k1+1:l1%nv), &
                                l1%yv(k1+1:l1%nv), l1%xv(k1), l1%yv(k1), l2%xv(k2), l2%yv(k2), &
                                xint, yint, sint)
                            
                            if (size(xint) > 0) then 
                                ! Exit the loop
                                exit
                            end if  
                        end if 

                        ! If we got here, update k2
                        k2 = k2 + 1
                    end do 

                    ! Set the maximal value
                    tc%tubes(j)%l1maxLOS(k1) = k2

                end do 

                ! Determine LOS for l2
                do k2 = 1, l2%nv 
                    ! minLOS
                    k1 = 1
                    do while (k1 <= l1%nv)
                        ! Compute segment intersections
                        if (k2 > 1) then 
                            ! Check for intersections in l1(1:k1-1)
                            call SegmentSimplePolygonIntersections(l2%xv(1:k2-1), &
                                l2%yv(1:k2-1), l2%xv(k2), l2%yv(k2), l1%xv(k1), l1%yv(k1), &
                                xint1, yint, sint)
                        end if 
                        if (k1 < l1%nv) then 
                            ! Check for intersections in l1(k1+1:nv)
                            call SegmentSimplePolygonIntersections(l2%xv(k2+1:l2%nv), &
                                l2%yv(k2+1:l2%nv), l1%xv(k1), l1%yv(k1), l2%xv(k2), l2%yv(k2), &
                                xint2, yint, sint)
                        end if 

                        ! Check
                        if (size(xint1) == 0 .and. size(xint2) == 0) then 
                            exit
                        else 
                            ! Update k1
                            k1 = k1 + 1
                        end if
                    end do 

                    ! Set the minimal value
                    tc%tubes(j)%l2minLOS(k2) = k1

                    ! maxLOS
                    k1 = tc%tubes(j)%l2minLOS(k2) ! start from minLOS
                    do while (k1 <= l1%nv)
                        ! Compute segment intersections 
                        if (k2 > 1) then 
                            ! Check for intersections in l1(1:k1-1)
                            call SegmentSimplePolygonIntersections(l2%xv(1:k2-1), &
                                l2%yv(1:k2-1), l1%xv(k1), l1%yv(k1), l2%xv(k2), l2%yv(k2), &
                                xint, yint, sint)
                            
                            if (size(xint) > 0) then 
                                ! Exit the loop
                                exit
                            end if  
                        end if 
                        if (k1 < l1%nv) then 
                            ! Check for intersections in l1(k1+1:nv)
                            call SegmentSimplePolygonIntersections(l2%xv(k2+1:l2%nv), &
                                l2%yv(k2+1:l2%nv), l1%xv(k1), l1%yv(k1), l2%xv(k2), l2%yv(k2), &
                                xint, yint, sint)
                            
                            if (size(xint) > 0) then 
                                ! Exit the loop
                                exit
                            end if  
                        end if 

                        ! If we got here, update k1
                        k1 = k1 + 1
                    end do 

                    ! Set the maximal value
                    tc%tubes(j)%l2maxLOS(k2) = k1
                end do 
            end do 


            ! Housekeeping
            end associate
        end do 


    end subroutine

    ! Removal of duplicate faces in grid
    subroutine RemoveDuplicateGridFaces(ggtmdata, grid)

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
        type(GGTMDataUDT), intent(inout)        :: ggtmdata

        ! Auxiliary
        integer(I8)                             :: vtemp, thisfID
        integer(I8), allocatable, dimension(:)  :: v1, v2, fID, &
            sortind, dv1, dv2, tsortind, tv2, tfID, fIDnew
        integer(I8), allocatable                :: fvert(:, :)
        logical, allocatable, dimension(:)      :: delind 
        
        ! Loop
        integer(I8)                             :: i, j, k, kold

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

        ! Construct mapping
        fIDnew = fID
        thisfID = fID(1)
        do i = 1, size(fID)
            if (.not. delind(i)) then 
                thisfID = fID(i)
            end if 
            fIDnew(i) = thisfID
        end do

        ! Remove faces from grid
        fID = pack(fID, delind)
        call RemoveGGFace(grid, fID)

        ! Remap faces in cell tube data 
        do i = 1, size(ggtmdata%cell)
            do j = 1, size(ggtmdata%cell(i)%tubes)
                associate(tube  => ggtmdata%cell(i)%tubes(j))
                call tube%hfface%Set(fIDnew(tube%hfface%Get()))
                call tube%lfface%Set(fIDnew(tube%lfface%Get()))
                call tube%tubeface%Set(fIDnew(tube%tubeface%Get()))
                end associate
            end do 
        end do

        ! Housekeeping
        !=============
        end associate

    end subroutine

    !------------------------------------------------------------------!
    !                         GRID ADAPTATION                          !
    !------------------------------------------------------------------!

    ! Non-convex cell splitting
    subroutine SplitNonConvexCells(ggtmdata, grid)

        ! Description
        !============
        ! This routine removes non-convex cells by splitting the cell
        ! at non-convex parts into triangles. Only 'radial' faces can 
        ! be introduced at this stage, and no additional vertices are
        ! added (faces are added, cells are deleted and added). Tube 
        ! data of the ggtmdata structure is adequately adjusted to 
        ! keep a correct representation of the grid. 

        ! Note: currently, we only support splitting quads, since at 
        ! this stage, no cells with more than 4 vertices are generated.
        ! This is checked for in the beginning. 

        ! Note: though ggtmdata is updated, the cells and faces will very likely
        ! not be in the correct order anymore. To this end, one should
        ! call sorting routines that reconstruct these. 

        ! Declare variables
        !==================
        ! Arguments
        type(GGTMDataUDT), intent(inout)        :: ggtmdata 
        type(GGGridUDT), intent(inout)          :: grid 

        ! Auxiliary
        integer(I8)                             :: ntv, npos, nneg, &
            tvind, nnewc
        integer(I8), allocatable, dimension(:)  :: tv, vertind, &
            splitvertind, vp1, vp2, newcellvert, cellvert, othercellvert, &
            newfacelabel, splitcellind, newcells, tc, cellfacemapping, &
            newcellregion, newfaceregion
        integer(I8), allocatable, dimension(:, :)   :: newfacevert, &
            newcellvertP, cellmapping
        real(R8), allocatable, dimension(:)     :: dx, dy, sintheta, &
            x, y
        logical, allocatable, dimension(:)      :: splitcell, issplitcell

        ! Loop
        integer(I8)                             :: i, j, k, cc

        ! Initialize
        !===========
        ! Associate
        associate(&
            celldata        => ggtmdata%cell,   &
            cell            => grid%cell,       &
            face            => grid%face,       &
            vert            => grid%vert        &
            )

        ! Initialize
        allocate(cellmapping(grid%cell%ntot, 2), cellfacemapping(grid%cell%ntot))
        cellmapping = 0 ! mapping from splitted cell to two new cells
        cellfacemapping = 0
        vp1 = cell%vp1%Get()
        vp2 = cell%vp2%Get()
        cellvert = cell%vert%Get()
        vertind = [2, 3, 4, 1]
        x = vert%x%Get()
        y = vert%y%Get()

        ! Check
        if (any(cell%vp2%Get() > 4)) then 
            call gdErrorHandler('SplitNonConvexCells: cells with more than ' // & 
                '4 vertices detected, not yet supported')
        end if 

        ! Mark cells
        !===========
        ! Initialize
        allocate(splitcell(cell%ntot), splitvertind(cell%ntot))
        splitcell = .false. 
        splitvertind = 0

        ! Determine which cells to split
        do i = 1, cell%ntot
            ! Check amount
            if (vp2(i) < 4) then 
                cycle
            end if 

            ! Get vertices
            tv = GetArrayFromPointer(vp1, vp2, cellvert, i)
            ntv = size(tv)

            ! Sanity check
            if (ntv /= 4) then 
                call gdErrorHandler('SplitNonConvexCells: this is a bug')
            end if 

            ! Close for ease
            tv = [tv, tv(1:2)]

            ! Compute angles 
            dx = x(tv(2:ntv+2)) - x(tv(1:ntv+1))
            dy = y(tv(2:ntv+2)) - y(tv(1:ntv+1))
            sintheta = dx(1:ntv)*dy(2:ntv+1) - dx(2:ntv+1)*dy(1:ntv)

            ! Check sign
            npos = count(sintheta >= 0.0_R8)
            nneg = count(sintheta < 0.0_R8)
            if ((npos == ntv) .or. (nneg == ntv)) then 
                ! All good, skip
                cycle
            elseif ((npos == 1) .or. (nneg == 1)) then 
                ! Mark for splitting
                splitcell(i) = .true. 
                if (npos == 1) then 
                    splitvertind(i) = vertind(findloc(sintheta >= 0.0_R8, .true., 1))
                else
                    splitvertind(i) = vertind(findloc(sintheta < 0.0_R8, .true., 1))
                end if 
            else
                ! Cell is self-intersecting, cannot split. Issue message
                print *, 'SplitNonConvexCells: cell ', i, ' is self-intersecting, ' // & 
                    'can not split. Continuing...'
            end if
        end do

        ! Check if there are any cells to remove. If not, return
        nnewc = count(splitcell)
        if (nnewc == 0) then 
            return 
        end if 

        ! Split cells
        !============
        ! Print message
        print *, 'SplitNonConvexCells: non-convex quads detected, splitting...'

        ! Initialize
        allocate(newfacevert(nnewc, 2), newfacelabel(nnewc)) ! one face per splitted quad
        allocate(newcellvertP(2*nnewc, 2), newcellvert(nnewc*6)) ! two triangles per splitted quad
        allocate(newcellregion(2*nnewc), newfaceregion(nnewc))
        newfacevert = 0
        newcellvert = 0
        newfacelabel = 0
        newcellregion = 0
        newfaceregion = 0
        newcellvertP(:, 1) = [(k, k = 1, 2*nnewc*3-2, 3)]
        newcellvertP(:, 2) = 3

        ! Loop
        cc = 0
        do i = 1, cell%ntot
            if (splitcell(i)) then 
                ! Update counter
                cc = cc + 1

                ! Get vertices
                tv = GetArrayFromPointer(vp1, vp2, cellvert, i)

                ! Check which vertices form the new face 
                tvind = mod(splitvertind(i) + 2, 4) ! because we know only 4 vertices exist...
                if (tvind == 0) then 
                    tvind = 4
                end if 

                ! Add the new face vertices & region
                newfacevert(cc, :) = [tv(splitvertind(i)), tv(tvind)]
                newfaceregion(cc) = cell%region%Get(i)


                ! Add the new cell vertices & region
                call SetDiff(tv, newfacevert(cc, :), othercellvert)
                newcellvert(6*cc-5:6*cc-3) = [newfacevert(cc, :), othercellvert(1)]
                newcellvert(6*cc-2:6*cc) = [newfacevert(cc, :), othercellvert(2)]
                newcellregion(2*cc-1:2*cc) = cell%region%Get(i)

                ! Construct the cell mapping
                cellmapping(i, :) = [2*cc-1, 2*cc]
                cellfacemapping(i) = face%ntot + cc

            end if 
        end do 

        ! Add constant to cellmapping
        cellmapping = cellmapping + cell%ntot

        ! Remap ggmdata
        !==============
        ! If we need to sort afterwards anyway, no need to put in effort
        ! to sort here. Just append/replace... (though nearly sortedness
        ! may likely lead to speed-up, so we do a bit of effort)
        do i = 1, size(celldata)
            do j = 1, size(celldata(i)%tubes)
                ! Associate for ease
                associate(tube  => celldata(i)%tubes(j))

                ! Check tube cells
                issplitcell = splitcell(tube%cell%Get())

                ! Replace if present
                if (any(issplitcell)) then 
                    ! Get indices
                    allocate(splitcellind(count(issplitcell)))
                    splitcellind = pack([(k, k = 1, size(issplitcell))], &
                        issplitcell)

                    ! Loop
                    do k = 1, size(splitcellind)
                        ! Get current tube cells
                        tc = tube%cell%Get()
                        
                        ! Delete the old cell
                        newcells = cellmapping(tc(splitcellind(k)), :)
                        call tube%cell%Set(newcells(1), splitcellind(k))
                        call tube%cell%Insert(newcells(2), splitcellind(k))

                        ! Insert the face at an approximate location
                        call tube%tubeface%Insert(cellfacemapping(tc(splitcellind(k))), &
                            splitcellind(k))

                        ! Update splitcellind to be consistent with extra
                        ! vertex
                        splitcellind = splitcellind + 1

                    end do

                    ! Housekeeping
                    deallocate(splitcellind)
                end if 

                ! Housekeeping
                end associate
            end do 
        end do

        ! Add and delete cells & faces
        !=============================
        ! Faces
        call AddGGFace(grid, newfacevert, newfacelabel, newfaceregion)

        ! Cells
        call AddGGCell(grid, newcellvert, newcellvertP, newcellregion)
        call RemoveGGCell(grid, pack([(k, k = 1, size(splitcell))], splitcell))

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

        ! Auxiliary
        
        ! Loop
        integer(I8)                         :: i 
        
        ! Initialize
        !===========
        ! Substructures (except segments)
        allocate(ggtmdata%vert(topomesh%vert%ntot))
        allocate(ggtmdata%face(topomesh%face%ntot))
        allocate(ggtmdata%cell(topomesh%cell%ntot))
        allocate(ggtmdata%tube(topomesh%tube%ntot))
        ggtmdata%nseg = 0_I8
        allocate(ggtmdata%seg(topomesh%vert%ntot + topomesh%face%ntot)) ! initial size

        ! Construct segments
        !===================
        ! Initialize
        associate(&
                vert        => topomesh%vert,   &
                vertdata    => ggtmdata%vert,   &
                face        => topomesh%face,   &
                facedata    => ggtmdata%face,   &
                nseg        => ggtmdata%nseg    &
                )

        ! Construct vertex 'segments' and 'lines'
        do i = 1, vert%ntot
            ! Update counter
            nseg = nseg + 1

            ! Initialize segment
            call ggtmdata%seg(nseg)%Initialize([vert%x(i), vert%x(i)], &
                [vert%y(i), vert%y(i)], vert%fsID(i), i, i, i)

            ! Initialize line
            call vertdata(i)%line%Initialize(ggtmdata, [nseg])
        end do 

        ! Construct face segments
        do i = 1, face%ntot
            ! Update counter
            nseg = nseg + 1

            ! Initialize segment
            call ggtmdata%seg(nseg)%Initialize(face%x(i)%Get(), &
                face%y(i)%Get(), face%fsID(i), i, face%vert(i, 1), face%vert(i, 2))

            ! Initialize line
            call facedata(i)%line%Initialize(ggtmdata, [nseg])
        end do 

        ! Housekeeping
        !=============
        end associate

    end subroutine 
    
    ! Vertex distribution over faces
    subroutine DistributeVerticesTopologicalMeshFaces(grid, ggtmdata, topomesh, &
        vd, GGTMlinerefiner, facetypes)

        ! Description
        !============
        ! Distribute vertices over topological mesh faces of the types
        ! defined in 'facetypes'. The vertices etc are stored in the
        ! ggtmdata structure. The vertex distribution is done based on 
        ! the vertexdistributor that is passed 

        ! Declare variables
        !==================
        ! Arguments
        class(GGGridUDT), intent(inout)         :: grid
        class(GGTMDataUDT)                      :: ggtmdata 
        class(TopomeshUDT), intent(in)          :: topomesh 
        class(VertexDistributor2DUDT), intent(in)   :: vd 
        class(GGTMLineRefiner2DUDT), intent(in) :: GGTMlinerefiner
        integer(I8), intent(in)                 :: facetypes(:)

        ! Auxiliary
        integer(I8)                             :: nv 
        integer(I8), allocatable, dimension(:)  :: tvID
        logical, allocatable, dimension(:)      :: keepvert
        real(R8), allocatable, dimension(:)     :: dlcv 

        ! Loop
        integer(I8)                             :: i, k

        ! Initialize
        !===========
        ! Associate
        associate(&
            nseg        => ggtmdata%nseg,   &
            vert        => topomesh%vert,   &
            face        => topomesh%face,   &
            facedata    => ggtmdata%face    &
            )
        
        ! Loop over all faces and distribute
        !===================================
        do i = 1, face%ntot
            if (any(face%type(i) == facetypes)) then

                ! Distribute
                call vd%DistributeOverCurve(face%x(i)%Get(), face%y(i)%Get(), &
                    nv, ldistr=dlcv)
                
                ! Adjust start and end to be sure
                dlcv(1) = 0
                dlcv(size(dlcv)) = facedata(i)%line%dllc(size(facedata(i)%line%dllc))
                
                ! Add vertex coordinates
                call facedata(i)%line%AddVertexCoordinates(dlcv)

                ! Add IDs
                tvID = [face%vert(i, 1), (k, k = grid%vert%ntot+1, grid%vert%ntot + size(dlcv)-2), face%vert(i, 2)]
                call facedata(i)%line%AddVertexIDs(tvID, &
                    [.true., spread(.false., 1, size(dlcv)-2), .true.])

                ! Update vertex counter
                grid%vert%ntot = grid%vert%ntot + size(dlcv) - 2

                ! Update the refiner
                call GGTMLineRefiner%UpdateRefinementOptions(&
                    facedata(i)%linerefoptions, topomesh)

                ! Refine
                keepvert = IsTopomeshVert(facedata(i)%line%vert, topomesh)
                call GGTMlinerefiner%Refine(facedata(i)%line, grid%vert%ntot, keepvert)

                ! Update segment data
                call facedata(i)%line%UpdateSegmentData(ggtmdata)
            end if 
        end do 

        ! Housekeeping
        !=============
        end associate

    end subroutine 

    ! Vertex distribution over tubes
    subroutine DistributeVerticesTopologicalMeshTubes(grid, ggtmdata, topomesh, &
        vd, GGTMlinerefiner, field)

        ! Description
        !============
        ! Distribute field values over topological mesh tubes. 
        ! The field values are computed per tube,
        ! first based on an initial distribution on the tube faces.
        ! This chosen distribution is simply 
        ! taken as the distribution that gives the maximal amount 
        ! of vertices, as this is likely the desired one. We store the 
        ! chosen topological face ID for each flux tube (we don't 
        ! propagate the distribution to other faces yet, since we have 
        ! to compute intersections anyway)

        ! Declare variables
        !==================
        ! Arguments
        class(GGGridUDT), intent(inout)         :: grid
        class(GGTMDataUDT)                      :: ggtmdata 
        class(TopomeshUDT), intent(in)          :: topomesh 
        class(VertexDistributor2DUDT), intent(in)    :: vd 
        class(DistributionFunctionUDT), intent(in)  :: field 
        class(GGTMLineRefiner2DUDT), intent(in) :: GGTMlinerefiner

        ! Auxiliary
        integer(I8)                             :: nv, nfl, nflmax, &
            tfmax
        integer(I8), allocatable, dimension(:)  :: tf, tvID, bndv1, &  
            bndv2
        logical, allocatable, dimension(:)      :: keepvert
        real(R8), allocatable, dimension(:)     :: xc, yc, dlcv

        ! Loop
        integer(I8)                             :: i, j, k  

        ! Initialize
        !===========
        ! Associate
        associate(&
            seg         => ggtmdata%seg,    &
            nseg        => ggtmdata%nseg,   &
            vert        => topomesh%vert,   &
            face        => topomesh%face,   &
            facedata    => ggtmdata%face,   &
            tube        => topomesh%tube,   &
            tubedata    => ggtmdata%tube    &
            )
        
        ! Determine tube distribution
        !============================
        do i = 1, tube%ntot 
            ! Get tube faces and boundary vertices (to check later if we
            ! need to project distributions etc)
            tf = tube%GetFace(i)
            bndv1 = tube%GetBndVert(i, 1)
            bndv2 = tube%GetBndVert(i, 2)

            ! Initialize
            nfl = 0
            nflmax = 0 

            ! Loop for the first time to construct initial distribution
            do k = 1, size(tf)
                ! Unpack
                j = tf(k)
                xc = face%x(j)%Get()
                yc = face%y(j)%Get()

                ! Update the refiner
                call GGTMLineRefiner%UpdateRefinementOptions(&
                    facedata(j)%linerefoptions, topomesh)

                ! Distribute
                !call vd%DistributeOverCurve(xc, yc, nv, ldistr=dlcv)
                call vd%DistributeOverField(xc, yc, field, nv,  ldistr=dlcv)
                dlcv(1) = 0
                dlcv(size(dlcv)) = facedata(j)%line%dllc(size(facedata(j)%line%dllc))

                ! Add data
                call facedata(j)%line%AddVertexCoordinates(dlcv)

                ! Add vertex IDs
                tvID = [face%vert(j, 1), (k, k = grid%vert%ntot+1, grid%vert%ntot + size(dlcv) - 2), face%vert(j, 2)]
                call facedata(j)%line%AddVertexIDs(tvID, &
                    [.true., spread(.false., 1, size(dlcv)-2), .true.])

                ! Update
                call facedata(j)%line%UpdateSegmentData(ggtmdata)

                ! Update vertex ID
                grid%vert%ntot = grid%vert%ntot + size(dlcv) - 2

                ! Refine 
                keepvert = IsTopomeshVert(facedata(j)%line%vert, topomesh)
                call GGTMlinerefiner%Refine(facedata(j)%line, grid%vert%ntot, keepvert)

                ! Update
                call facedata(j)%line%UpdateSegmentData(ggtmdata)

            end do 


            ! Loop
            nflmax = size(facedata(tf(1))%line%xv)
            tfmax = tf(1)
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

        ! Housekeeping
        !=============
        end associate

    end subroutine 

    ! Cell data
    subroutine AddTopologicalMeshGriddingData(ggtmdata, topomesh, &
        fieldtracer, magneticField, options)

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
        type(GGOptionsUDT), intent(in)          :: options

        ! Auxiliary
        integer(I8)                             :: tc, srf, erf, inderf, &
            indsrf, minind, maxind, minindloc
        integer(I8), allocatable, dimension(:)  :: tubec, tubef, tcf, &
            tcv, tcfv1, tcfv2, hffaces, lffaces, hfvert, lfvert, &
            tf1, tf2
        real(R8)                                :: hfval, lfval, &
            dhf1, dhf2, dlf1, dlf2
        real(R8), allocatable, dimension(:)     :: tcvfval, tcfv1val, &
            tcfv2val
        logical, allocatable, dimension(:)      :: ishfface, islfface, &
            ishfvert

        ! Diagnostics

        ! Loop
        integer(I8)                             :: i, j, k

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
                    if (any(face%type(tcf(k)) == TMfacealignedID)) then 
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
                minindloc = minloc([indsrf, inderf], 1)
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

                ! Sort to make sure tf1 and tf2 always go from starting to ending radial face
                if (minindloc == 1) then 
                    ! tf2 is correctly sorted, need to flip tf1
                    if (size(tf1) > 0) then 
                        tf1 = tf1(size(tf1):1:-1)
                    end if 
                else
                    ! tf1 is correctly sorted, need to flip tf2
                    if (size(tf2) > 0) then 
                        tf2 = tf2(size(tf2):1:-1)
                    end if 
                end if 

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

    
        ! Add refinement data
        !====================
        call AddTopologicalMeshLineRefinementData(ggtmdata, &
            topomesh, options)
        
        ! Housekeeping
        !=============
        end associate


    end subroutine 

    ! Topological mesh tube contours
    subroutine TraceTopologicalMeshTubeContours(grid, ggtmdata, topomesh, &
        fieldtracer, magneticField, options)

        ! Description
        !============
        ! This routine traces the field contours using the field tracer
        ! for each topological mesh tube. 

        ! Notes
        !======
        ! Note 1: it is assumed that the tube cell and face data is
        ! properly sorted (some sanity checks are done) and that the
        ! initial distributions on the topomesh tube are set.

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

        ! Note 3: it is possible that, due to a difference in contour 
        ! tracer resolution, fieldlines may intersect with other 
        ! field lines or lines present in the topological mesh. These 
        ! lines will be removed afterwards - this may be avoided by 
        ! reconstructing everything from scratch and not loading
        ! anything in. Then, the topomesh resolution should be 
        ! consistent between topomesh and grid generation stage. 

        ! Declare variables
        !==================
        ! Arguments
        class(GGGridUDT), intent(inout)         :: grid
        class(GGTMDataUDT)                      :: ggtmdata 
        class(TopomeshUDT), intent(in)          :: topomesh 
        class(ContourTracerUDT), intent(in)     :: fieldtracer 
        type(MagneticFieldUDT), intent(in)      :: magneticField
        type(GGOptionsUDT), intent(in)          :: options

        ! Auxiliary
        integer(I8)                             :: &
            tf, cind, nc, ntf, incr, nv, temp, &
            startind, endind, nfs, tfloc, tfc, nthf, ntlf
        integer(I8), allocatable, dimension(:)  :: tubec, tubef, &
            allIDs, s1, s2, polv, fsID, sortind, thf, tlf, tvertexID, &
            tf1, tf2, allnbtf, contourind
        integer(I8), allocatable, dimension(:, :)   :: nint, segrf, &
            segc, vertexID, temp2
        real(R8)                                :: &
            nxc, nyc, nxfv, nyfv, txf, tyf, ntxf, tmaxval
        real(R8), allocatable, dimension(:)     :: &
            tx, ty, xl, yl, tfval, sr1, sr2, txint, tyint, &
            nxf, nyf, nnf, tsegrc, tsegrrf, dlcv, newdlcv, newtfval
        real(R8), allocatable, dimension(:, :)  :: segrrf, segrc, &
            xint, yint
        logical                                 :: isflremoved_nointersect, &
            isintersectremoved, issrf, doflip, changesign, &
            isflremoved_multipleintersect, doremoval
        logical, allocatable, dimension(:)      :: &
            iscontourfound, keepind, isdescending
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
            nseg        => ggtmdata%nseg,   &
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
        nv = grid%vert%ntot
        doremoval = .true.

        ! Trace contours 
        !===============
        ! Initialize flux surface counter
        nfs = topomesh%nfs

        ! Loop over all tubes
        do i = 1, tube%ntot 

            ! Initialize
            !-----------
            ! Get tube cells and faces
            tubef = tube%GetFace(i)
            tubec = tube%GetCell(i)
            tf = tubedata(i)%distributionface
            tf1 = tube%GetBndFace(i, 1)
            tf2 = tube%GetBndFace(i, 2)
            allnbtf = [tf1, tf2]

            ! Get cell belonging to this face
            tfloc = findloc(tubef, tf, 1, back=.false.)
            if (tfloc == 0) then 
                ! This should not happen
                print *, 'tube: ', i, 'face: ', tf 
                call gdErrorHandler('AddTopologicalMeshCellGriddingData: ' // &
                    'could not find distribution face in tube, this is a bug' )
            elseif (tfloc == size(tubec)+1) then
                ! Cell is last cell
                tfc = tubec(size(tubec))
            else 
                ! Default case
                tfc = tubec(tfloc)
            end if 

            ! Determine contour starting points
            !----------------------------------
            ! Get initial vertex distribution 
            tx = facedata(tf)%line%xv
            ty = facedata(tf)%line%yv
            dlcv = facedata(tf)%line%dlcv

            ! Make sure we trace from high to low field value and 
            ! ensure that values are monotonous and within the hf and 
            ! lf bound
            if (allocated(tfval)) then 
                deallocate(tfval)
            end if 
            allocate(tfval(size(tx)))
            tfval = fieldtracer%Evaluate(tx, ty)
            if (tfval(1) < tfval(size(tfval))) then 
                tfval = tfval(size(tfval):1:-1)
                dlcv = dlcv(size(dlcv):1:-1)
            end if 

            ! Eliminate values outside of bounds
            allocate(keepind(size(tfval)))
            keepind = .true.
            where ((tfval > tfval(1) ).or. (tfval < tfval(size(tfval))))
                keepind = .false. 
            end where
            allocate(newtfval(count(keepind)),newdlcv(count(keepind)))
            newtfval = pack(tfval, keepind)
            newdlcv = pack(dlcv, keepind)

            ! Eliminate parts that are not (strictly) descending
            allocate(isdescending(size(newtfval)))
            isdescending = .true. 
            tmaxval = newtfval(1)
            do j = 2, size(newtfval) ! this should work since we eliminated out of bounds values
                if (newtfval(j) < tmaxval) then 
                    tmaxval = newtfval(j)
                else
                    isdescending(j) = .false. 
                end if 
            end do 
            if (.not. isdescending(size(isdescending))) then 
                isdescending(size(isdescending)-1) = .false.
                isdescending(size(isdescending)) = .true.
            end if 
            newdlcv = pack(newdlcv, isdescending)

            ! Eliminate lines 
            
            ! Sort from high to low values
            !if (allocated(sortind)) then 
            !    deallocate(sortind)
            !end if 
            !allocate(sortind(size(newtfval)))
            !call Sort(newtfval, ind=sortind, ascend=.false.)
            !newdlcv = newdlcv(sortind)

            ! Redistribute (vertIDs etc are updated later)
            call facedata(tf)%line%AddVertexCoordinates(newdlcv)
            tx = facedata(tf)%line%xv
            ty = facedata(tf)%line%yv

            ! Housekeeping
            deallocate(newdlcv, newtfval, keepind, isdescending)

            ! Trace contours
            !---------------
            ! Trace
            tempc = fieldtracer%TraceContours(tx(2:size(tx)-1), ty(2:size(ty)-1))

            ! Clean
            call CleanContours(tempc)

            ! Coarsen
            if (options%coarsencontours) then 
                call fieldtracer%CoarsenContours(tempc)
            end if 

            ! Concatenate & orient contours
            !------------------------------
            ! Precompute face normals for each flux surface for 
            ! determining orientation
            nxf = -(facedata(tf)%line%yl(2:) - facedata(tf)%line%yl(1:size(facedata(tf)%line%yl)-1))
            nyf = (facedata(tf)%line%xl(2:) - facedata(tf)%line%xl(1:size(facedata(tf)%line%xl)-1))
            nxf = [nxf, nxf(size(nxf))]
            nyf = [nyf, nyf(size(nyf))]
            nnf = sqrt(nxf**2 + nyf**2)
            nxf = nxf/nnf
            nyf = nyf/nnf

            ! Find high field/low field face that has the same vertex
            thf = celldata(tfc)%hffaces
            nthf = size(thf)
            tlf = celldata(tfc)%lffaces 
            ntlf = size(tlf)
            issrf = .false. 
            if (celldata(tfc)%srf == tf) then 
                issrf = .true.
            end if 
            doflip = .false.
            changesign = .false.
            if (nthf > 0) then 
                ! Check with high field face (should be first or last face)
                associate(hf1       => facedata(thf(1))%line, &
                    hf2             => facedata(thf(nthf))%line)

                if (topomesh%face%vert(thf(1), 1) == topomesh%face%vert(tf, 1)) then 
                    ! Take first edge of high field line
                    !print *, 'AddTopologicalMeshCellGriddingData: code part ' // & 
                    !'not yet verified'
                    txf = (hf1%xl(2) - hf1%xl(1))
                    tyf = (hf1%yl(2) - hf1%yl(1))
                elseif (topomesh%face%vert(thf(1), 1) == topomesh%face%vert(tf, 2)) then
                    ! Take first edge of high field line, flip 
                    !print *, 'AddTopologicalMeshCellGriddingData: code part ' // & 
                    !'not yet verified'
                    txf = (hf1%xl(2) - hf1%xl(1))
                    tyf = (hf1%yl(2) - hf1%yl(1))
                    nxf = nxf(size(nxf):1:-1) 
                    nyf = nyf(size(nyf):1:-1) 
                    doflip = .true.
                elseif (topomesh%face%vert(thf(1), 2) == topomesh%face%vert(tf, 1)) then
                    ! Take last edge of high field line
                    !print *, 'AddTopologicalMeshCellGriddingData: code part ' // & 
                    !'not yet verified'
                    txf = (hf1%xl(hf1%nl-1) - hf1%xl(hf1%nl))
                    tyf = (hf1%yl(hf1%nl-1) - hf1%yl(hf1%nl))
                elseif (topomesh%face%vert(thf(1), 2) == topomesh%face%vert(tf, 2)) then
                    ! Take last edge of high field line, flipe
                    !print *, 'AddTopologicalMeshCellGriddingData: code part ' // & 
                    !'not yet verified'
                    txf = (hf1%xl(hf1%nl-1) - hf1%xl(hf1%nl))
                    tyf = (hf1%yl(hf1%nl-1) - hf1%yl(hf1%nl))
                    nxf = nxf(size(nxf):1:-1) 
                    nyf = nyf(size(nyf):1:-1) 
                    doflip = .true.
                elseif (topomesh%face%vert(thf(nthf), 1) == topomesh%face%vert(tf, 1)) then 
                    ! Take first edge of high field line
                    !print *, 'AddTopologicalMeshCellGriddingData: code part ' // & 
                    !'not yet verified'
                    txf = (hf2%xl(2) - hf2%xl(1))
                    tyf = (hf2%yl(2) - hf2%yl(1))
                elseif (topomesh%face%vert(thf(nthf), 1) == topomesh%face%vert(tf, 2)) then
                    ! Take first edge of high field line, flip 
                    !print *, 'AddTopologicalMeshCellGriddingData: code part ' // & 
                    !'not yet verified'
                    txf = (hf2%xl(2) - hf2%xl(1))
                    tyf = (hf2%yl(2) - hf2%yl(1))
                    nxf = nxf(size(nxf):1:-1) 
                    nyf = nyf(size(nyf):1:-1) 
                    doflip = .true.
                elseif (topomesh%face%vert(thf(nthf), 2) == topomesh%face%vert(tf, 1)) then
                    ! Take last edge of high field line
                    !print *, 'AddTopologicalMeshCellGriddingData: code part ' // & 
                    !'not yet verified'
                    txf = (hf2%xl(hf2%nl-1) - hf2%xl(hf2%nl))
                    tyf = (hf2%yl(hf2%nl-1) - hf2%yl(hf2%nl))
                elseif (topomesh%face%vert(thf(nthf), 2) == topomesh%face%vert(tf, 2)) then
                    ! Take last edge of high field line, flip
                    !print *, 'AddTopologicalMeshCellGriddingData: code part ' // & 
                    !'not yet verified'
                    txf = (hf2%xl(hf2%nl-1) - hf2%xl(hf2%nl))
                    tyf = (hf2%yl(hf2%nl-1) - hf2%yl(hf2%nl))
                    nxf = nxf(size(nxf):1:-1) 
                    nyf = nyf(size(nyf):1:-1) 
                    doflip = .true.
                else
                    ! No correspondence between vertices found, throw error
                    call gdErrorHandler('AddTopologicalMeshCellGriddingData: ' // & 
                        'Could not find any corresponding vertices, this is a bug')
                end if 

                ! Check orientation
                !if (issrf) then 
                !    txf = -txf
                !    tyf = -tyf
                !end if 
                ntxf = sqrt(txf**2 + tyf**2)
                txf = txf/ntxf 
                tyf = tyf/ntxf
                if ((nxf(1)*txf + nyf(1)*tyf) < 0.0_R8) then 
                    nxf = -nxf 
                    nyf = -nyf 
                    changesign = .true.
                end if 

                end associate
    
            elseif (size(tlf) > 0) then 
                ! Check with high field face (should be first or last face)
                associate(lf1       => facedata(tlf(1))%line, &
                    lf2             => facedata(tlf(ntlf))%line)
                if (topomesh%face%vert(tlf(1), 1) == topomesh%face%vert(tf, 2)) then 
                    ! Take first edge of high field line
                    !print *, 'AddTopologicalMeshCellGriddingData: code part ' // & 
                    !'not yet verified'
                    txf = (lf1%xl(2) - lf1%xl(1))
                    tyf = (lf1%yl(2) - lf1%yl(1))
                elseif (topomesh%face%vert(tlf(1), 1) == topomesh%face%vert(tf, 1)) then
                    ! Take first edge of high field line, flip 
                    !print *, 'AddTopologicalMeshCellGriddingData: code part ' // & 
                    !'not yet verified'
                    txf = (lf1%xl(2) - lf1%xl(1))
                    tyf = (lf1%yl(2) - lf1%yl(1))
                    nxf = nxf(size(nxf):1:-1) 
                    nyf = nyf(size(nyf):1:-1) 
                    doflip = .true.
                elseif (topomesh%face%vert(tlf(1), 2) == topomesh%face%vert(tf, 2)) then
                    ! Take last edge of high field line
                    !print *, 'AddTopologicalMeshCellGriddingData: code part ' // & 
                    !'not yet verified'
                    txf = (lf1%xl(lf1%nl-1) - lf1%xl(lf1%nl))
                    tyf = (lf1%yl(lf1%nl-1) - lf1%yl(lf1%nl))
                elseif (topomesh%face%vert(tlf(1), 2) == topomesh%face%vert(tf, 1)) then
                    ! Take last edge of high field line, flipe
                    !print *, 'AddTopologicalMeshCellGriddingData: code part ' // & 
                    !'not yet verified'
                    txf = (lf1%xl(lf1%nl-1) - lf1%xl(lf1%nl))
                    tyf = (lf1%yl(lf1%nl-1) - lf1%yl(lf1%nl))
                    nxf = nxf(size(nxf):1:-1) 
                    nyf = nyf(size(nyf):1:-1) 
                    doflip = .true.
                elseif (topomesh%face%vert(tlf(ntlf), 1) == topomesh%face%vert(tf, 2)) then 
                    ! Take first edge of high field line
                    !print *, 'AddTopologicalMeshCellGriddingData: code part ' // & 
                    !'not yet verified'
                    txf = (lf2%xl(2) - lf2%xl(1))
                    tyf = (lf2%yl(2) - lf2%yl(1))
                elseif (topomesh%face%vert(tlf(ntlf), 1) == topomesh%face%vert(tf, 1)) then
                    ! Take first edge of high field line, flip 
                    !print *, 'AddTopologicalMeshCellGriddingData: code part ' // & 
                    !'not yet verified'
                    txf = (lf2%xl(2) - lf2%xl(1))
                    tyf = (lf2%yl(2) - lf2%yl(1))
                    nxf = nxf(size(nxf):1:-1) 
                    nyf = nyf(size(nyf):1:-1) 
                    doflip = .true.
                elseif (topomesh%face%vert(tlf(ntlf), 2) == topomesh%face%vert(tf, 2)) then
                    ! Take last edge of high field line
                    !print *, 'AddTopologicalMeshCellGriddingData: code part ' // & 
                    !'not yet verified'
                    txf = (lf2%xl(lf2%nl-1) - lf2%xl(lf2%nl))
                    tyf = (lf2%yl(lf2%nl-1) - lf2%yl(lf2%nl))
                elseif (topomesh%face%vert(tlf(ntlf), 2) == topomesh%face%vert(tf, 1)) then
                    ! Take last edge of high field line, flip
                    !print *, 'AddTopologicalMeshCellGriddingData: code part ' // & 
                    !'not yet verified'
                    txf = (lf2%xl(lf2%nl-1) - lf2%xl(lf2%nl))
                    tyf = (lf2%yl(lf2%nl-1) - lf2%yl(lf2%nl))
                    nxf = nxf(size(nxf):1:-1) 
                    nyf = nyf(size(nyf):1:-1) 
                    doflip = .true.
                else
                    ! No correspondence between vertices found, throw error
                    call gdErrorHandler('AddTopologicalMeshCellGriddingData: ' // & 
                        'Could not find any corresponding vertices, this is a bug')
                end if 

                ! Check orientation
                !if (issrf) then 
                !    txf = -txf
                !    tyf = -tyf
                !end if 
                ntxf = sqrt(txf**2 + tyf**2)
                txf = txf/ntxf 
                tyf = tyf/ntxf
                if ((nxf(size(nxf))*txf + nyf(size(nyf))*tyf) < 0.0_R8) then 
                    nxf = -nxf 
                    nyf = -nyf 
                    changesign = .true.
                end if 

                end associate
            else 
                ! Cell with only tangency points, not supported
                call gdErrorHandler('AddTopologicalMeshCellGriddingData: ' // &
                    'cell has no actual faces for high and low field boundaries, ' // & 
                    'not supported')
            end if 

            ! Recompute face normals, now with (coarser) vertex values
            if (facedata(tf)%line%dlcv(1) < facedata(tf)%line%dlcv(facedata(tf)%line%nv)) then 
                nxf = -(facedata(tf)%line%yv(2:) - facedata(tf)%line%yv(1:size(facedata(tf)%line%yv)-1))
                nyf = (facedata(tf)%line%xv(2:) - facedata(tf)%line%xv(1:size(facedata(tf)%line%xv)-1))
            else
                nxf = (facedata(tf)%line%yv(2:) - facedata(tf)%line%yv(1:size(facedata(tf)%line%yv)-1))
                nyf = -(facedata(tf)%line%xv(2:) - facedata(tf)%line%xv(1:size(facedata(tf)%line%xv)-1))
            end if 
            nnf = sqrt(nxf**2 + nyf**2)
            nxf = nxf/nnf
            nyf = nyf/nnf
            !nxf = nxf2 
            !nyf = nyf2
            if (doflip) then 
                nxf = nxf(size(nxf):1:-1)
                nyf = nyf(size(nyf):1:-1)
            end if 
            if (changesign) then 
                nxf = -nxf
                nyf = -nyf
            end if 
            
            ! Check if contours make sense and reformat if necessary
            allIDs = tempc%ID
            allocate(iscontourfound(size(tx)-2)) ! normally, IDs go from 1 to number of points
            allocate(keepind(size(tempc)))
            iscontourfound = .false. 
            keepind = .true. 
            if (tube%isclosed(i)) then 
                ! Precompute face normals for each flux surface for 
                ! determining orientation
                !nxf = -(facedata(tf)%line%yv(2:) - facedata(tf)%line%yv(1:size(facedata(tf)%line%yv)-1))
                !nyf = (facedata(tf)%line%xv(2:) - facedata(tf)%line%xv(1:size(facedata(tf)%line%xv)-1))
                !nnf = sqrt(nxf**2 + nyf**2)
                !nxf = nxf/nnf
                !nyf = nyf/nnf

                !if (tfval(1) < tfval(size(tfval))) then 
                !    nxf = -nxf(size(nxf):1:-1) 
                !    nyf = -nyf(size(nyf):1:-1) 
                !end if

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
                    associate(tID       => tempc(j)%ID)
                    nxc = tempc(j)%x(2) - tempc(j)%x(1)
                    nyc = tempc(j)%y(2) - tempc(j)%y(1)
                    nxfv = 0.5*(nxf(tID) + nxf(tID+1)) ! need to work with ID here, since contours may be open etc
                    nyfv = 0.5*(nyf(tID) + nyf(tID+1))
                    if ((nxc*nxfv + nyc*nyfv) < 0) then 
                        ! Flip the contour
                        tempc(j)%x = tempc(j)%x(size(tempc(j)%x):1:-1)
                        tempc(j)%y = tempc(j)%y(size(tempc(j)%y):1:-1)
                        temp = tempc(j)%startsaddle
                        tempc(j)%startsaddle = tempc(j)%endsaddle
                        tempc(j)%endsaddle = temp
                    end if 
                    end associate

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

                            ! Set found tot rue
                            iscontourfound(c1%ID) = .true.

                            ! Ensure proper orientation
                            associate(tID   => tempc(j)%ID)
                            nxc = c1%x(2) - c1%x(1)
                            nyc = c1%y(2) -  c1%y(1)
                            nxfv = 0.5*(nxf(tID) + nxf(tID+1)) ! need to work with ID here, since contours may be open etc
                            nyfv = 0.5*(nyf(tID) + nyf(tID+1))
                            if ((nxc*nxfv + nyc*nyfv) < 0) then 
                                ! Flip the contour
                               c1%x =c1%x(size(tempc(j)%x):1:-1)
                               c1%y =c1%y(size(tempc(j)%y):1:-1)
                                temp =c1%startsaddle
                               c1%startsaddle =c1%endsaddle
                               c1%endsaddle = temp
                            end if 
                            end associate

                            ! Adjust 
                            tempc(j) = c1

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

            ! Check contour-tube intersections
            !---------------------------------
            ! Initialize
            allocate(keepind(size(tempc)))
            keepind = .true.
             

            ! Check if contour lines intersect with the tube boundary 
            ! faces
            !$omp parallel do default(none) schedule(dynamic) collapse(2) &
            !$omp private(j, k, txint, tyint, s1, s2, sr1, sr2) &
            !$omp shared(tempc, allnbtf, grid, topomesh, keepind)
            do j = 1, size(tempc)
                ! Check for intersections with aligned faces
                do k = 1, size(allnbtf)
                    if (keepind(j)) then  ! no need to check if already found intersection
                        txint = spread(0, 1, 0)
                        call SimplePolygonIntersections(tempc(j)%x, tempc(j)%y, &
                            face%x(allnbtf(k))%Get(), face%y(allnbtf(k))%Get(), &
                            txint, tyint, s1, s2, sr1, sr2)
                        if (size(txint) > 0) then 
                            !$omp critical
                            keepind(j) = .false. 
                            !$omp end critical
                        end if 
                    end if 
                end do 
            end do
            !$omp end parallel do

            ! Remove contours
            if (any(.not. keepind)) then 
                ! Issue message
                print *, 'TraceTopologicalMeshTubeContours: some contours ' // & 
                    'intersect with tube boundaries, probably due to ' // & 
                    'different contouring resolution between gridding ' // &
                    'and topomesh construction. Removing these contours...'

                ! Remove
                tempc = pack(tempc, keepind)
            end if 
                

            ! Housekeeping
            deallocate(keepind)

            ! Compute intersections 
            !----------------------
            ! Initialize
            nc = size(tempc)
            ntf = size(tubef)
            allocate(xintda(nc, ntf), yintda(nc, ntf), nint(nc, ntf), &
                segrfda(nc, ntf), segcda(nc, ntf), segrrfda(nc, ntf), &
                segrcda(nc, ntf), polc(nc))
            nint = 0

            ! Loop over all contours
            call wall_time(tstart)
            !$omp parallel do default(none) schedule(dynamic) collapse(2) &
            !$omp private(i, j, txint, tyint, s1, s2, sr1, sr2) &
            !$omp shared(tempc, tubef, xintda, yintda, segrfda, segcda, &
            !$omp segrrfda, segrcda, nint)
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
            !$omp end parallel do
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
            isflremoved_nointersect = .false. 
            isflremoved_multipleintersect = .false.
            isintersectremoved = .false. 
            do j = 1, nc
                ! Check if no intersections with radial lines
                if (any(nint(j, :) == 0)) then 
                    ! Mark for removal
                    keepind(j) = .false. 
                    isflremoved_nointersect = .true.
                    
                    ! Go to the next line
                    cycle 
                end if 

                ! Check if multiple intersections with radial lines
                if (tube%isclosed(i)) then 
                    
                    if (tfloc == 1 .or. tfloc == size(tubef)) then 
                        ! Two intersections are expected in first and 
                        ! last face, since they are the same face (normally)
                        ! Also last face
                        if (any(nint(j, 2:ntf-1) > 1)) then 
                            if (doremoval) then 
                                ! Mark for removal
                                keepind(j) = .false. 
                                isflremoved_multipleintersect = .true.
                            else
                                ! Simply set warning message
                                isintersectremoved = .true. 
                            end if 
                        end if 
                        if ((nint(j, 1) > 2) .or. (nint(j, ntf) > 2)) then 
                            if (doremoval) then 
                                ! Mark for removal
                                keepind(j) = .false. 
                                isflremoved_multipleintersect = .true.
                            else
                                ! Simply set warning message
                                isintersectremoved = .true. 
                            end if
                        end if 
                    else
                        ! Only one face, but should have two intersections
                        if (any(nint(j, [(k, k = 1, tfloc-1), (k, k = tfloc+1, ntf)]) > 1)) then 
                            if (doremoval) then 
                                ! Mark for removal
                                keepind(j) = .false. 
                                isflremoved_multipleintersect = .true.
                            else
                                ! Simply set warning message
                                isintersectremoved = .true. 
                            end if
                        end if 
                        if (nint(j, tfloc) > 2) then 
                            if (doremoval) then 
                                ! Mark for removal
                                keepind(j) = .false. 
                                isflremoved_multipleintersect = .true.
                            else
                                ! Simply set warning message
                                isintersectremoved = .true. 
                            end if
                        end if 
                    end if 
                else
                    if (tempc(j)%isclosed) then 
                        ! Two intersections are expected in the tracing
                        ! face only
                        if (any(nint(j, [(k, k = 1, tfloc-1), (k, k = tfloc+1, ntf)]) > 1)) then 
                            if (doremoval) then 
                                ! Mark for removal
                                keepind(j) = .false. 
                                isflremoved_multipleintersect = .true.
                            else
                                ! Simply set warning message
                                isintersectremoved = .true. 
                            end if
                        end if 
                        if (nint(j, tfloc) > 2) then 
                            if (doremoval) then 
                                ! Mark for removal
                                keepind(j) = .false. 
                                isflremoved_multipleintersect = .true.
                            else
                                ! Simply set warning message
                                isintersectremoved = .true. 
                            end if
                        end if 
                    elseif (any(nint(j, :) > 1)) then 
                        if (doremoval) then 
                            ! Mark for removal
                            keepind(j) = .false. 
                            isflremoved_multipleintersect = .true.
                        else
                            ! Simply set warning message
                            isintersectremoved = .true. 
                        end if
                    end if 
                end if 
            end do 

            ! Remove field lines
            tempc = pack(tempc, keepind)
            nc = count(keepind)

            ! Issue warnings
            if (isflremoved_nointersect) then 
                if (verbosity > 0) then 
                    print *, 'AddTopologicalMeshCellGriddingData: ' // & 
                        'field lines were removed since they do not intersect ' // & 
                        'with one or more radial lines for tube: ', i 
                end if 
            end if 
            if (isflremoved_multipleintersect) then 
                if (verbosity > 0) then 
                    print *, 'AddTopologicalMeshCellGriddingData: ' // & 
                        'field lines were removed since they intersect multiple times ' // & 
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
                vertexID(nc, ntf), fsID(nc), temp2(nc, ntf), &
                contourind(nc))
            !do j = 1, size(nint, 2)
            !    temp2(:, j) = pack(nint(:, j), keepind)
            !end do 
            !nint = temp2
            cc = 0 
            do j = 1, size(nint, 1) 
                ! Skip
                if (keepind(j)) then 
                    ! Update counter
                    cc = cc + 1 
                    contourind(cc) = j

                    ! Intersection in tracing face: should always be the one
                    ! that is at the start of the contour IF it is a 
                    ! closed polygon! Otherwise, we just take the first one...
                    k = tfloc
                    if (tempc(cc)%isclosed) then 
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
                        ! Find intersection that is closest in distance
                        ! to the original tracing starting point (we 
                        ! should actually find that one up to some 
                        ! precision)
                        if (nint(j, k) == 1)  then 
                            ! Trivial case
                            startind = 1
                        else
                            ! Get all coordinates
                            txint = xintda(j, k)%Get()
                            tyint = yintda(j, k)%Get()

                            ! Check which is closest
                            startind = minloc((txint - tx(j+1))**2 + &
                                (tyint - ty(j+1))**2, 1)
                        end if 
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
                    do k = 1, ntf-1 
                        ! Skip the tracing radial face
                        if (k == tfloc) then 
                            cycle 
                        end if 

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
                        if (tfloc == 1 .or. tfloc == size(tubef)) then 
                            ! Intersection with first and last radial line
                            ! should be exactly the same! Only intersection
                            ! coordinate should differ
                            tsegrc = segrcda(cc, ntf)%Get()
                            endind = findloc(tsegrc, real(size(tempc(cc)%x)-1, kind=R8), 1)
                            if (endind == 0) then 
                                print *, 'AddToplogicalMeshCellGriddingData: ' // & 
                                    'closed contour does not intersect at ending point. ' // & 
                                    'contour: ', cc, 'tube: ', i 
                                print *, 'taking first intersection...'
                                endind = 1
                                call Write2DPolygonData(tempc(cc)%x, tempc(cc)%y, 'l1')
                            end if 

                        else

                            ! Just a regular end for last face
                            endind = 1

                        end if 

                        ! First and last face are always the same 
                        ! normally speaking, so vertex should be 
                        ! the same
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
            !if (tube%isclosed(i)) then 
            !    ! First and last intersection should be at 0 
            !    ! and ne coordinate by construction
            !    do j = 1, nc
            !        if (segrc(j, 1) /= 0_R8) then 
            !            call gdErrorHandler('AddTopologicalMeshCellGriddingData: ' // & 
            !                'first intersection of closed contour is not at start of ' // & 
            !                'contour for closed flux tube, unexpected')
            !        end if
            !        if (segrc(j, size(segrc, 2)) /= size(tempc(j)%x)-1) then 
            !            call gdErrorHandler('AddTopologicalMeshCellGriddingData: ' // & 
            !                'last intersection of closed contour is not at end of ' // & 
            !                'contour for closed flux tube, unexpected')
            !        end if 
            !    end do 
            !end if 

            ! Add lines for each cell
            do k = 1, ntf-1
                ! Allocate the amount of lines for this cell
                if (allocated(celldata(tubec(k))%lines)) then
                    deallocate(celldata(tubec(k))%lines)
                end if 
                allocate(celldata(tubec(k))%lines(nc))

                
                ! Loop over lines
                do j = 1, nc
                    ! Hedge for closed tubes of which the starting radial 
                    ! face is not the first or last one (need to adjust 
                    ! segrc on the fly there... First segrc is 0, second is
                    ! size of the line)
                    if (k+1 == tfloc) then 
                        if (tube%isclosed(i) .and. (tfloc /= 1 .or. tfloc /= size(tubef))) then 
                            ! Need to adjust starting position
                            cc = contourind(j)
                            tsegrc = segrcda(cc, k+1)%Get()
                            endind = findloc(tsegrc, real(size(tempc(cc)%x)-1, kind=R8), 1)
                            if (endind == 0) then 
                                print *, 'AddToplogicalMeshCellGriddingData: ' // & 
                                    'closed contour does not intersect at ending point. ' // & 
                                    'contour: ', cc, 'tube: ', i 
                                print *, 'taking first intersection...'
                                endind = 1
                                call Write2DPolygonData(tempc(cc)%x, tempc(cc)%y, 'l1')
                            end if 

                            xint(j, k+1) = xintda(cc, k+1)%Get(endind)
                            yint(j, k+1) = yintda(cc, k+1)%Get(endind)
                            segc(j, k+1) = segcda(cc, k+1)%Get(endind)
                            segrf(j, k+1) = segrfda(cc, k+1)%Get(endind)
                            segrc(j, k+1) = segrcda(cc, k+1)%Get(endind)
                            segrrf(j, k+1) = segrrfda(cc, k+1)%Get(endind)

                        end if 
                    elseif (k == tfloc) then 
                        if (tube%isclosed(i) .and. (tfloc /= 1 .or. tfloc /= size(tubef))) then 
                            ! Need to adjust starting position
                            cc = contourind(j)
                            tsegrc = segrcda(cc, k)%Get()
                            startind = findloc(tsegrc, 0.0_R8, 1)
                            if (startind == 0) then 
                                print *, 'AddToplogicalMeshCellGriddingData: ' // & 
                                    'closed contour does not intersect at starting point. ' // & 
                                    'contour: ', cc, 'tube: ', i 
                                print *, 'taking first intersection...'
                                startind = 1
                                call Write2DPolygonData(tempc(cc)%x, tempc(cc)%y, 'l1')
                            end if 

                            xint(j, k) = xintda(cc, k)%Get(startind)
                            yint(j, k) = yintda(cc, k)%Get(startind)
                            segc(j, k) = segcda(cc, k)%Get(startind)
                            segrf(j, k) = segrfda(cc, k)%Get(startind)
                            segrc(j, k) = segrcda(cc, k)%Get(startind)
                            segrrf(j, k) = segrrfda(cc, k)%Get(startind)

                        end if 
                    end if

                    ! Check line order
                    if (segrc(j, k) < segrc(j, k+1)) then 
                        incr = 1
                    else
                        incr = -1
                    end if 

                    ! Note: we ensure no duplicate points by skipping the
                    ! first vertex of the face segment
                    if (incr > 0) then 
                        polv = [(cc, cc = segc(j, k)+2, segc(j, k+1)-1, incr)]
                    else 
                        polv = [(cc, cc = segc(j, k)-2, segc(j, k+1)+1, incr)]
                    end if 
                    xl = [xint(j, k), tempc(j)%x(polv), xint(j, k+1)]
                    yl = [yint(j, k), tempc(j)%y(polv), yint(j, k+1)]
                    
                    ! Add segment
                    nseg = nseg + 1
                    if (nseg > size(ggtmdata%seg)) then 
                        ! Double the size
                        ggtmdata%seg = [ggtmdata%seg, ggtmdata%seg]
                    end if 
                    if (vertexID(j, k) == vertexID(j, k+1)) then 
                        ! Closed segment, ensure the same to machine precision
                        xl(size(xl)) = xl(1)
                        yl(size(yl)) = yl(1)
                    end if 
                    call ggtmdata%seg(nseg)%Initialize(xl, yl, fsID(j), &
                        0_I8, vertexID(j, k), vertexID(j, k+1))

                    ! Add 
                    call celldata(tubec(k))%lines(j)%Initialize(ggtmdata, [nseg])

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

            ! Adjust distributions of other tube faces
            do j = 1, ntf
                ! Get intersections
                tsegrrf = segrrf(:, j)
                txint = xint(:, j)
                tyint = yint(:, j)
                tvertexID = vertexID(:, j)

                ! Sort
                sortind = [(k, k= 1, size(tsegrrf))]
                call Sort(tsegrrf, ind=sortind)
                txint = txint(sortind)
                tyint = tyint(sortind)
                tvertexID = tvertexID(sortind)

                ! Reconstruct line
                call facedata(tubef(j))%line%AddVertexCoordinates([0.0_R8, &
                    tsegrrf, real(facedata(tubef(j))%line%nl-1, kind=R8)], .true.)

                ! Add IDs
                call facedata(tubef(j))%line%AddVertexIDs(&
                    [face%vert(tubef(j), 1), tvertexID, face%vert(tubef(j), 2)], &
                    [.true., spread(.false., 1, size(tvertexID)), .true.])

                ! Update line data
                call facedata(tubef(j))%line%UpdateSegmentData(ggtmdata)
                
            end do 

            ! Housekeeping
            deallocate(xint, yint, segc, segrf, segrc, segrrf, nint, &
                vertexID, polc, xintda, yintda, segcda, segrfda, &
                segrcda, segrrfda, fsID, temp2, contourind)
        end do 

        ! Housekeeping
        !=============
        ! Trim
        ggtmdata%seg = ggtmdata%seg(1:nseg)

        ! Update
        grid%vert%ntot = nv
        end associate

    end subroutine

    ! Cell tube construction
    subroutine ConstructTopologicalMeshCellFluxTubes(grid, ggtmdata, topomesh, &
        fieldtracer, magneticField, options)

        ! Description
        !============
        ! This routine constructs for each topological mesh cell the 
        ! flux tubes resulting from (a.o.) the contour tracing done
        ! in TraceTopologicalMeshTubeContours. Each flux tube is stored
        ! as a GGTMFieldlineDataPair type, which may hold additional 
        ! information on line of sight etc. The lines in these flux tubes,
        ! contrary to the lines stored in the cell data, can be extended
        ! with vessel parts for example to allow cut cell approaches 
        ! etc. These tubes should be used when generating the grid, 
        ! not the cell field lines (this was the case in the past, though)

        ! Options for checks, cut cell approach, etc can be set in the
        ! general options structure

        ! Declare variables
        !==================
        ! Arguments
        type(GGGridUDT), intent(inout)          :: grid
        type(GGTMDataUDT), intent(inout)       :: ggtmdata
        class(TopomeshUDT), intent(in)          :: topomesh 
        class(ContourTracerUDT), intent(in)     :: fieldtracer 
        type(MagneticFieldUDT), intent(in)      :: magneticField 
        type(GGOptionsUDT), intent(in)          :: options

        ! Auxiliary
        integer(I8)                             :: nt, startsegID, &
            endsegID, vind1, vind2, tsegID
        integer(I8), allocatable, dimension(:)  :: tc, srfvert, erfvert, &
            allsegID, uerfvert, usrfvert
        logical, allocatable, dimension(:)      :: dostart, doend
        real(R8)                                :: dl 

        ! Loop
        integer(I8)                             :: i, j, k

        ! Initialize
        !===========
        ! Unpack for ease
        associate(&
            seg             => ggtmdata%seg,    &
            facedata        => ggtmdata%face,   &
            celldata        => ggtmdata%cell,   &
            tube            => topomesh%tube,   &
            cell            => topomesh%cell    &
            )

        ! Check if we need to split at 

        ! Construct basic tubes
        !======================
        ! Simple tubes without any extension, so purely flux surfaces
        do i = 1, cell%ntot
            ! Initialize
            !-----------
            ! Associate for ease
            associate(tc        => celldata(i))

            ! Sanity checks
            if (.not. allocated(tc%lines)) then 
                print *, 'cell: ', i 
                call gdErrorHandler('ConstructTopologicalMeshCellFluxTubes: ' // & 
                    'cell lines for this cell are not yet allocated')
            end if 

            ! Initialize
            nt = size(tc%lines)+1
            allocate(tc%tubes(nt))

            ! Extract the high/low field line data 
            call ExtractTMCellAlignedBoundary(tc, 'high', ggtmdata, &
                topomesh, tc%hfline)
            call ExtractTMCellAlignedBoundary(tc, 'low', ggtmdata, &
                topomesh, tc%lfline)

            ! Construct tubes
            !----------------
            ! Initialize
            do j = 1, nt 
                call tc%tubes(j)%Initialize()
            end do 

            ! First tube, starts with hfline
            tc%tubes(1)%hfline = tc%hfline
            if (nt > 1) then 
                tc%tubes(1)%lfline = tc%lines(1)
            else
                tc%tubes(1)%lfline = tc%lfline
            end if 

            ! Intermediate tubes
            do j = 2, nt-1
                tc%tubes(j)%hfline = tc%lines(j-1)
                tc%tubes(j)%lfline = tc%lines(j)
            end do 

            ! Final tube
            if (nt > 1) then 
                tc%tubes(nt)%hfline = tc%lines(nt-1)
                tc%tubes(nt)%lfline = tc%lfline 
            end if 

            ! Add labels
            do j = 1, nt
                tc%tubes(j)%srflabel = tc%srflabel
                tc%tubes(j)%erflabel = tc%erflabel
            end do 

            ! Housekeeping
            end associate
        end do 

        ! Extend with vessel parts
        !=========================
        ! Extend the current simple tubes with vessel parts in some 
        ! cases
        do i = 1, cell%ntot
            ! Initialize
            !-----------
            ! Associate for ease
            associate(&
                vertID  => grid%vert%ntot,      &
                tc      => celldata(i),         &
                srfline     => facedata(celldata(i)%srf)%line,   &
                erfline     => facedata(celldata(i)%erf)%line,   &
                tubes   => celldata(i)%tubes    &
                )

            ! Initialize
            nt = size(tubes) 

            ! Tangency point tubes
            !---------------------
            if (options%extendtptubes) then 
                ! Check first tube
                if (ggtmdata%seg(tubes(1)%hfline%segID(1))%isvertex) then 
                    ! Get starting radial line, update its data
                    call facedata(celldata(i)%srf)%line%UpdateLineData(ggtmdata)
                    
                    ! Split the line 
                    call facedata(celldata(i)%srf)%line%SplitAtVertex(&
                        tubes(1)%lfline%vert(1), ggtmdata)
                    
                    ! Get the correct segment for extension
                    startsegID = 0
                    allsegID = facedata(celldata(i)%srf)%line%segID
                    do k = 1, size(allsegID)
                        if (ggtmdata%seg(allsegID(k))%sv == tubes(1)%hfline%vert(1) .or. &
                            ggtmdata%seg(allsegID(k))%ev == tubes(1)%hfline%vert(1)) then 
                            startsegID = allsegID(k)
                            exit 
                        end if 
                    end do 

                    ! Sanity check
                    if (startsegID == 0) then 
                        print *, 'cell: ', i, 'tube: ', 1
                        call gdErrorHandler('ConstructTopologicalMeshCellFluxTubes: ' // & 
                            'could not find start segment for extension')
                    end if 

                    ! Get the ending radial line, update its data
                    call facedata(celldata(i)%erf)%line%UpdateLineData(ggtmdata)

                    ! Split the line
                    call facedata(celldata(i)%erf)%line%SplitAtVertex(&
                        tubes(1)%lfline%vert(tubes(1)%lfline%nv), ggtmdata)

                    ! Get the correct segment for extension
                    endsegID = 0
                    allsegID = facedata(celldata(i)%erf)%line%segID
                    do k = 1, size(allsegID)
                        if (ggtmdata%seg(allsegID(k))%sv == tubes(1)%hfline%vert(1) .or. &
                            ggtmdata%seg(allsegID(k))%ev == tubes(1)%hfline%vert(1)) then 
                            endsegID = allsegID(k)
                            exit 
                        end if 
                    end do 

                    ! Sanity check
                    if (endsegID == 0) then 
                        print *, 'cell: ', i, 'tube: ', 1
                        call gdErrorHandler('ConstructTopologicalMeshCellFluxTubes: ' // & 
                            'could not find end segment for extension')
                    end if 

                    ! Construct line
                    call tubes(1)%hfline%Initialize(ggtmdata, [startsegID, endsegID])

                    ! Check if we need to flip
                    if (ggtmdata%seg(tubes(1)%hfline%segID(1))%TMfaceID == celldata(i)%erf) then 
                        ! Need to flip
                        call tubes(1)%hfline%Flip()
                    elseif (ggtmdata%seg(tubes(1)%hfline%segID(1))%TMfaceID == celldata(i)%srf) then 
                        ! All good
                    else
                        ! All bad
                        call gdErrorHandler('ConstrutTopologicalMeshCellFluxTubes: ' // & 
                            'first segment does not come from starting or ' // & 
                            'ending radial face of cell')
                    end if 

                    ! Set as extended
                    tubes(1)%isextendedstart = .true.
                    tubes(1)%isextendedend = .true.
                end if 

                ! Check last tube
                if (ggtmdata%seg(tubes(nt)%lfline%segID(1))%isvertex) then 
                    ! Get starting radial line, update its data
                    call facedata(celldata(i)%srf)%line%UpdateLineData(ggtmdata)
                    
                    ! Split the line 
                    call facedata(celldata(i)%srf)%line%SplitAtVertex(&
                        tubes(nt)%hfline%vert(1), ggtmdata)
                    
                    ! Get the correct segment for extension
                    startsegID = 0
                    allsegID = facedata(celldata(i)%srf)%line%segID
                    do k = 1, size(allsegID)
                        if (ggtmdata%seg(allsegID(k))%sv == tubes(nt)%lfline%vert(1) .or. &
                            ggtmdata%seg(allsegID(k))%ev == tubes(nt)%lfline%vert(1)) then 
                            startsegID = allsegID(k)
                            exit 
                        end if 
                    end do 

                    ! Sanity check
                    if (startsegID == 0) then 
                        print *, 'cell: ', i, 'tube: ', 1
                        call gdErrorHandler('ConstructTopologicalMeshCellFluxTubes: ' // & 
                            'could not find start segment for extension')
                    end if 

                    ! Get the ending radial line, update its data
                    call facedata(celldata(i)%erf)%line%UpdateLineData(ggtmdata)

                    ! Split the line
                    call facedata(celldata(i)%erf)%line%SplitAtVertex(&
                        tubes(nt)%hfline%vert(tubes(nt)%hfline%nv), ggtmdata)

                    ! Get the correct segment for extension
                    endsegID = 0
                    allsegID = facedata(celldata(i)%erf)%line%segID
                    do k = 1, size(allsegID)
                        if (ggtmdata%seg(allsegID(k))%sv == tubes(nt)%lfline%vert(1) .or. &
                            ggtmdata%seg(allsegID(k))%ev == tubes(nt)%lfline%vert(1)) then 
                            endsegID = allsegID(k)
                            exit 
                        end if 
                    end do 

                    ! Sanity check
                    if (endsegID == 0) then 
                        print *, 'cell: ', i, 'tube: ', 1
                        call gdErrorHandler('ConstructTopologicalMeshCellFluxTubes: ' // & 
                            'could not find end segment for extension')
                    end if 

                    ! Construct line
                    call tubes(nt)%lfline%Initialize(ggtmdata, [startsegID, endsegID])

                    ! Check if we need to flip
                    if (ggtmdata%seg(tubes(nt)%lfline%segID(1))%TMfaceID == celldata(i)%erf) then 
                        ! Need to flip
                        call tubes(nt)%lfline%Flip()
                    elseif (ggtmdata%seg(tubes(nt)%lfline%segID(1))%TMfaceID == celldata(i)%srf) then 
                        ! All good
                    else
                        ! All bad
                        call gdErrorHandler('ConstrutTopologicalMeshCellFluxTubes: ' // & 
                            'first segment does not come from starting or ' // & 
                            'ending radial face of cell')
                    end if 

                    ! Set as extended
                    tubes(nt)%isextendedstart = .true.
                    tubes(nt)%isextendedend = .true.
                end if 
            end if 
                
            ! Other vessel tubes
            !-------------------
            if (options%extendvesseltubes) then
                ! Initialize
                allocate(dostart(nt), doend(nt), srfvert(2*nt), erfvert(2*nt))
                dostart = .false. 
                doend = .false. 
                srfvert = 0
                erfvert = 0

                ! Determine which tubes to treat
                do j = 1, nt
                    ! Associate for ease
                    associate(&
                        thfline         => tubes(j)%hfline, &
                        tlfline         => tubes(j)%lfline  &
                    )

                    ! Start
                    !------
                    ! Check if we can extend start  
                    if (.not. tubes(j)%isextendedstart .and. &
                        topomesh%face%type(tc%srf) == TMfacebndID &
                        .and. .not. (ggtmdata%seg(tubes(j)%hfline%segID(1))%isvertex) & 
                        .and. .not. (ggtmdata%seg(tubes(j)%lfline%segID(1))%isvertex)) then

                        ! Get length of segment
                        vind1 = findloc(srfline%vert, thfline%vert(1), 1)
                        vind2 = findloc(srfline%vert, tlfline%vert(1), 1)
                        if (vind1 == 0 .or. vind2 == 0) then 
                            call gdErrorHandler('ConstructTopologicalMeshCellFluxTubes: ' // & 
                                'could not find vertices in start radial face, unexpected')
                        end if 
                        dl = abs(srfline%dlcv(vind1) - srfline%dlcv(vind2))

                        ! Check if it is too large
                        if (dl > options%evtmaxvessellength) then 
                            dostart(j) = .true.
                            srfvert(2*j-1) = thfline%vert(1)
                            srfvert(2*j) = tlfline%vert(1) 
                        end if 
                    end if 

                    ! Check if we can extend end  
                    if (.not. tubes(j)%isextendedend .and. &
                        topomesh%face%type(tc%erf) == TMfacebndID &
                        .and. .not. (ggtmdata%seg(tubes(j)%hfline%segID(1))%isvertex) & 
                        .and. .not. (ggtmdata%seg(tubes(j)%lfline%segID(1))%isvertex)) then

                        ! Get length of segment
                        vind1 = findloc(erfline%vert, thfline%vert(thfline%nv), 1)
                        vind2 = findloc(erfline%vert, tlfline%vert(tlfline%nv), 1)
                        if (vind1 == 0 .or. vind2 == 0) then 
                            call gdErrorHandler('ConstructTopologicalMeshCellFluxTubes: ' // & 
                                'could not find vertices in end radial face, unexpected')
                        end if 
                        dl = abs(erfline%dlcv(vind1) - erfline%dlcv(vind2))

                        ! Check if it is too large
                        if (dl > options%evtmaxvessellength) then 
                            doend(j) = .true.
                            erfvert(2*j-1) = thfline%vert(thfline%nv)
                            erfvert(2*j) = tlfline%vert(tlfline%nv) 
                        end if 

                    end if 

                    ! Housekeeping
                    end associate
                end do

                ! Split the radial faces at these vertices
                call Unique(pack(erfvert, erfvert /= 0), uerfvert)
                call Unique(pack(srfvert, srfvert /= 0), usrfvert)
                call srfline%SplitAtVertices(usrfvert, ggtmdata)
                call erfline%SplitAtVertices(uerfvert, ggtmdata)

                ! Extend the tubes
                do j = 1, nt 
                    if (dostart(j)) then 
                        ! Get segment ID for extension
                        tsegID = GetSegmentIDFromVertices(srfline, ggtmdata, &
                            srfvert(2*j-1), srfvert(2*j))

                        ! Sanity check
                        if (tsegID == 0) then 
                            call gdErrorHandler('ConstructTopologicalMeshCellFluxTubes: ' // & 
                                'could not find starting segment, unexpected')
                        end if 

                        ! Extend the tube
                        call ExtendTubeWithSegment(i, tubes(j), tsegID, &
                            ggtmdata, .true., vertID, magneticField)

                    end if 
                    if (doend(j)) then 
                        ! Get segment ID for extension
                        tsegID = GetSegmentIDFromVertices(erfline, ggtmdata, &
                            erfvert(2*j-1), erfvert(2*j))

                        ! Sanity check
                        if (tsegID == 0) then 
                            call gdErrorHandler('ConstructTopologicalMeshCellFluxTubes: ' // & 
                                'could not find starting segment, unexpected')
                        end if 

                        ! Extend the tube
                        call ExtendTubeWithSegment(i, tubes(j), tsegID, &
                            ggtmdata, .false., vertID, magneticField)
                    end if 
                end do

                ! Housekeeping
                deallocate(srfvert, erfvert, dostart, doend)
            end if 


            ! Housekeeping
            end associate
        end do 
        
        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Cell data for refining lines 
    subroutine AddTopologicalMeshLineRefinementData(ggtmdata, topomesh, &
        options)

        ! Description
        !============
        ! This routine adds any additional data for line refinement. It 
        ! is assumed that the gridding data is already added, so we can
        ! use the topological mesh to more easily identify refinement 
        ! options etc. This routine writes out at the end a file 
        ! containing all refinement options for each cell. If desired,
        ! the user can adjust that information to influence the grid 
        ! generator by changing the file data and changing the method 
        ! to 'existing' and defining the filepath of the data. Typically 
        ! one would generate an initial set of data (or an original grid) 
        ! and later adjust that one.
        
        ! The refinement options here are propagated to any other 
        ! topological mesh object (faces, tubes, ...).

        ! Declare variables
        !==================
        ! Arguments
        type(GGTMDataUDT), intent(inout)        :: ggtmdata 
        class(TopomeshUDT), intent(in)          :: topomesh
        type(GGoptionsUDT), intent(in)          :: options 

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: targetfaceIDs, &
            vesselfaceIDs, strikepointIDs, tv1, tv2

        ! Loop
        integer(I8)                             :: i, j

        ! Initialize
        !===========
        ! Associate for ease
        associate(&
            celldata        => ggtmdata%cell,   &
            facedata        => ggtmdata%face,   &
            tubedata        => ggtmdata%tube,   &
            tube            => topomesh%tube,   &
            face            => topomesh%face,   &
            vert            => topomesh%vert)

        ! Set basic options
        !==================
        ! Checks to make when determining if cell can be constructed
        do i = 1, size(celldata)
            celldata(i)%legalcellstyle = options%legalcellstyle
        end do     

        ! Refinement options
        !===================
        ! Check if we simply read in the exisintg file
        if (options%readexistingrefdata) then 
            ! Read
            call ReadTopologicalMeshLineRefinementData(ggtmdata, options%refdatafile)

            ! Propagate
            call PropagateTopologicalMeshLineRefinementData(ggtmdata, topomesh)
            return
        end if 

        ! If we got here, this means that we didn't read in an existing
        ! file

        ! Initialize
        do i = 1, size(celldata)
            call celldata(i)%linerefoptions%Initialize()
        end do
        do i = 1, size(tubedata)
            call tubedata(i)%linerefoptions%Initialize()
        end do 
        
        ! Boundary layer, poloidal
        !-------------------------
        ! Vessel faces
        if (options%refBLdovessel) then 
            ! Get all vessel boundaries of the topological mesh
            vesselfaceIDs = topomesh%GetVesselFaceIDs()

            ! Check for each cell
            do i = 1, size(celldata)

                ! Add global data
                celldata(i)%linerefoptions%dlBLlengthbased = options%refdlBLlengthbased

                ! Check starting radial face
                if (any(celldata(i)%srf == vesselfaceIDs)) then 
                    ! Overwrite defaults
                    celldata(i)%linerefoptions%doBLstart = .true. 
                end if 

                ! Add data anyway
                celldata(i)%linerefoptions%ncBLstart = options%refBLncvessel
                celldata(i)%linerefoptions%dlBLstart = options%refBLdlvessel

                ! Check the ending radial face
                if (any(celldata(i)%erf == vesselfaceIDs)) then 
                    ! Overwrite defaults
                    celldata(i)%linerefoptions%doBLend = .true. 
                end if 

                ! Add data anyway
                celldata(i)%linerefoptions%ncBLend = options%refBLncvessel
                celldata(i)%linerefoptions%dlBLend = &
                    options%refBLdlvessel(size(options%refBLdlvessel):1:-1) ! need to flip
            end do 
        end if 

        ! Target faces (will overwrite existing vessel faces as intended)
        if (options%refBLdotarget) then 
            ! Get all target boundaries of the topological mesh
            targetfaceIDs = topomesh%GetTargetFaceIDs()

            ! Check for each cell
            do i = 1, size(celldata)

                ! Add global data
                celldata(i)%linerefoptions%dlBLlengthbased = options%refdlBLlengthbased

                ! Check starting radial face
                if (any(celldata(i)%srf == targetfaceIDs)) then 
                    ! Overwrite defaults
                    celldata(i)%linerefoptions%doBLstart = .true. 
                end if 

                ! Add data anyway
                celldata(i)%linerefoptions%ncBLstart = options%refBLnctarget
                celldata(i)%linerefoptions%dlBLstart = options%refBLdltarget

                ! Check the ending radial face
                if (any(celldata(i)%erf == targetfaceIDs)) then 
                    ! Overwrite defaults
                    celldata(i)%linerefoptions%doBLend = .true. 
                end if 

                ! Add data anyway
                celldata(i)%linerefoptions%ncBLend = options%refBLnctarget
                celldata(i)%linerefoptions%dlBLend = &
                    options%refBLdltarget(size(options%refBLdltarget):1:-1) ! need to flip
            end do 
        end if 

        ! Boundary layer, radial
        !-----------------------
        ! Target faces (will overwrite existing vessel faces as intended)
        if (options%radrefBLdosp) then 
            ! Get all target boundaries of the topological mesh
            targetfaceIDs = topomesh%GetTargetFaceIDs()

            ! Get all strike points of the topological mesh
            strikepointIDs = topomesh%GetStrikePointIDs()

            ! Check for each tube ('start' is the first boundary, 'end'
            ! is the second) 
            do i = 1, tube%ntot
                ! Get first and second tube boundary vertices
                tv1 = tube%GetBndVert(i, 1_I8)
                tv2 = tube%GetBndVert(i, 2_I8)

                ! Check if any are strike point IDs
                do j = 1, size(tv1)
                    if (any(tv1(j) == strikepointIDs) .or. vert%type(tv1(j)) == TMvertexsaddleID) then 
                        ! Do boundary layer at start
                        tubedata(i)%linerefoptions%doBLstart = .true. 
                    end if 
                end do 
                do j = 1, size(tv2)
                    if (any(tv2(j) == strikepointIDs) .or. vert%type(tv2(j)) == TMvertexsaddleID) then 
                        ! Do boundary layer at start
                        tubedata(i)%linerefoptions%doBLend = .true. 
                    end if 
                end do 

                ! Add data anyway
                tubedata(i)%linerefoptions%dlBLlengthbased  = options%radrefdlBLlengthbased
                tubedata(i)%linerefoptions%ncBLstart = options%radrefBLncsp
                tubedata(i)%linerefoptions%dlBLstart = options%radrefBLdlsp
                tubedata(i)%linerefoptions%ncBLend = options%radrefBLncsp
                tubedata(i)%linerefoptions%dlBLend = options%radrefBLdlsp(&
                    size(options%radrefBLdlsp):1:-1)

            end do 
        end if 

        ! Propagate options to faces
        !===========================
        call PropagateTopologicalMeshLineRefinementData(ggtmdata, topomesh)

        ! Write out options
        !==================
        call WriteTopologicalMeshLineRefinementData(ggtmdata, 'refdataTM')

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Topomesh refinement data propagation
    subroutine PropagateTopologicalMeshLineRefinementData(ggtmdata, topomesh)

        ! Description
        !============
        ! This routine propagates the topological mesh data from cells
        ! and tubes to faces (and, if required, other objects of the 
        ! topomesh)

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMDataUDT), intent(inout)       :: ggtmdata 
        class(TopomeshUDT), intent(in)          :: topomesh

        ! Auxiliary
        integer(I8)                             :: nhf, nlf
        integer(I8), allocatable, dimension(:)  :: tv1, tv2, tf

        
        ! Loop
        integer(I8)                             :: i, j

        ! Initialize
        !===========
        ! Associate for ease
        associate(&
            celldata        => ggtmdata%cell,   &
            tubedata        => ggtmdata%tube,   &
            facedata        => ggtmdata%face,   &
            cell            => topomesh%cell,   &
            face            => topomesh%face,   &
            tube            => topomesh%tube,   &
            vert            => topomesh%vert)

        ! Initialize
        do i = 1, size(facedata)
            call facedata(i)%linerefoptions%Initialize()
        end do

        ! Propagate to faces
        !===================
        ! Boundary layer data (cells)
        !----------------------------
        ! Propagate options to aligned faces at start/end of cell.
        ! Preference for BL imposition is imposed implicitly
        do i = 1, size(celldata)
            ! Associate
            associate(&
                refoptions  => celldata(i)%linerefoptions,  & 
                hf      => celldata(i)%hffaces, &
                lf      => celldata(i)%lffaces, &
                srf     => celldata(i)%srf,     &
                erf     => celldata(i)%erf      &
                )

            ! Initialize
            nhf = size(hf)
            nlf = size(lf)

            ! Check boundary layer at cell start
            if (refoptions%doBLstart) then 
                ! Check hf face
                if (nhf > 0) then 
                    ! Add general data
                    facedata(hf(1))%linerefoptions%dlBLlengthbased = refoptions%dlBLlengthbased 

                    ! Check if we need to do specific things
                    if (any(face%vert(hf(1), 1) == face%vert(srf, :))) then 
                        ! Set start
                        facedata(hf(1))%linerefoptions%doBLstart = refoptions%doBLstart
                        facedata(hf(1))%linerefoptions%ncBLstart = refoptions%ncBLstart
                        facedata(hf(1))%linerefoptions%dlBLstart = refoptions%dlBLstart
                    else
                        ! Set end
                        facedata(hf(1))%linerefoptions%doBLend = refoptions%doBLstart
                        facedata(hf(1))%linerefoptions%ncBLend = refoptions%ncBLstart
                        facedata(hf(1))%linerefoptions%dlBLend = &
                            refoptions%dlBLstart(refoptions%ncBLstart:1:-1)
                    end if 
                end if

                ! Check lf face
                if (nlf > 0) then 
                    ! Add general data
                    facedata(lf(1))%linerefoptions%dlBLlengthbased = refoptions%dlBLlengthbased 

                    if (any(face%vert(lf(1), 1) == face%vert(srf, :))) then 
                        ! Set start
                        facedata(lf(1))%linerefoptions%doBLstart = refoptions%doBLstart
                        facedata(lf(1))%linerefoptions%ncBLstart = refoptions%ncBLstart
                        facedata(lf(1))%linerefoptions%dlBLstart = refoptions%dlBLstart
                    else
                        ! Set end
                        facedata(lf(1))%linerefoptions%doBLend = refoptions%doBLstart
                        facedata(lf(1))%linerefoptions%ncBLend = refoptions%ncBLstart
                        facedata(lf(1))%linerefoptions%dlBLend = &
                            refoptions%dlBLstart(refoptions%ncBLstart:1:-1)
                    end if 
                end if 
            end if 

            ! Check boundary layer at cell end
            if (refoptions%doBLend) then 
                ! Check hf face
                if (nhf > 0) then
                    ! Add general data
                    facedata(hf(nhf))%linerefoptions%dlBLlengthbased = refoptions%dlBLlengthbased 

                    if (any(face%vert(hf(nhf), 1) == face%vert(erf, :))) then 
                        ! Set start
                        facedata(hf(nhf))%linerefoptions%doBLstart = refoptions%doBLend
                        facedata(hf(nhf))%linerefoptions%ncBLstart = refoptions%ncBLend
                        facedata(hf(nhf))%linerefoptions%dlBLstart = &
                            refoptions%dlBLend(refoptions%ncBLend:1:-1)
                    else
                        ! Set end
                        facedata(hf(nhf))%linerefoptions%doBLend = refoptions%doBLend
                        facedata(hf(nhf))%linerefoptions%ncBLend = refoptions%ncBLend
                        facedata(hf(nhf))%linerefoptions%dlBLend = refoptions%dlBLend
                    end if 
                end if

                ! Check lf face
                if (nlf > 0) then
                    ! Add general data
                    facedata(lf(nlf))%linerefoptions%dlBLlengthbased = refoptions%dlBLlengthbased 

                    if (any(face%vert(lf(nlf), 1) == face%vert(erf, :))) then 
                        ! Set start
                        facedata(lf(nlf))%linerefoptions%doBLstart = refoptions%doBLend
                        facedata(lf(nlf))%linerefoptions%ncBLstart = refoptions%ncBLend
                        facedata(lf(nlf))%linerefoptions%dlBLstart = &
                            refoptions%dlBLend(refoptions%ncBLend:1:-1)
                    else
                        ! Set end
                        facedata(lf(nlf))%linerefoptions%doBLend = refoptions%doBLend
                        facedata(lf(nlf))%linerefoptions%ncBLend = refoptions%ncBLend
                        facedata(lf(nlf))%linerefoptions%dlBLend = refoptions%dlBLend
                    end if 
                end if 
            end if 

            ! Housekeeping
            end associate
        end do 

        ! Propagate options to hf and lf refoptions

        ! Boundary layer data (tubes)
        !----------------------------
        ! Only need to propagate to radial faces of tubes
        do i = 1, size(tubedata)
            ! Associate for ease
            associate(refoptions    => tubedata(i)%linerefoptions)
            
            ! Get tube faces
            tf = tube%GetFace(i)

            ! Get tube boundary vertices
            tv1 = tube%GetBndVert(i, 1_I8)
            tv2 = tube%GetBndVert(i, 2_I8)

            ! For each face, check if BL should be applied
            if (refoptions%doBLstart) then 
                do j = 1, size(tf)
                    ! Set general data
                    facedata(tf(j))%linerefoptions%dlBLlengthbased = refoptions%dlBLlengthbased

                    ! Check face vertices
                    if (any(face%vert(tf(j), 1) == tv1)) then 
                        facedata(tf(j))%linerefoptions%doBLstart = .true. 
                        facedata(tf(j))%linerefoptions%ncBLstart = refoptions%ncBLstart
                        facedata(tf(j))%linerefoptions%dlBLstart = refoptions%dlBLstart
                    elseif (any(face%vert(tf(j), 2) == tv1)) then
                        facedata(tf(j))%linerefoptions%doBLend = .true. 
                        facedata(tf(j))%linerefoptions%ncBLend = refoptions%ncBLstart
                        facedata(tf(j))%linerefoptions%dlBLend = refoptions%dlBLstart(size(refoptions%dlBLstart):1:-1)
                    else
                        ! This shouldn't happen
                        call gdErrorHandler('PropagateTopologicalMeshLineRefinementData: ' // & 
                            'radial face of tube does not have a vertex in common ' // & 
                            'with tube boundary vertices - this is a bug')
                    end if 
                end do 
            end if 

            ! For each face, check if BL should be applied
            if (refoptions%doBLend) then 
                do j = 1, size(tf)
                    ! Set general data
                    facedata(tf(j))%linerefoptions%dlBLlengthbased = refoptions%dlBLlengthbased

                    ! Check face vertices
                    if (any(face%vert(tf(j), 1) == tv2)) then 
                        facedata(tf(j))%linerefoptions%doBLstart = .true. 
                        facedata(tf(j))%linerefoptions%ncBLstart = refoptions%ncBLend
                        facedata(tf(j))%linerefoptions%dlBLstart = refoptions%dlBLend(size(refoptions%dlBLend):1:-1)
                    elseif (any(face%vert(tf(j), 2) == tv2)) then
                        facedata(tf(j))%linerefoptions%doBLend = .true. 
                        facedata(tf(j))%linerefoptions%ncBLend = refoptions%ncBLend
                        facedata(tf(j))%linerefoptions%dlBLend = refoptions%dlBLend
                    else
                        ! This shouldn't happen
                        call gdErrorHandler('PropagateTopologicalMeshLineRefinementData: ' // & 
                            'radial face of tube does not have a vertex in common ' // & 
                            'with tube boundary vertices - this is a bug')
                    end if 
                end do 
            end if 


            ! Housekeeping
            end associate

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

    subroutine InitializeGGTMLineRefinementOptions(options)

        ! Description
        !============
        ! Option initializer. Initially, every option is set to false.

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMFieldlineRefinementOptionsUDT)    :: options 

        ! Set defaults
        !=============
        ! Boundary layer options
        options%dlBLlengthbased = .false.
        options%doBLstart = .false.
        options%ncBLstart = 0
        options%doBLend = .false. 
        options%ncBLend = 0
        if (allocated(options%dlBLstart)) then 
            deallocate(options%dlBLstart)
        end if 
        if (allocated(options%dlBLend)) then 
            deallocate(options%dlBLend)
        end if 
        allocate(options%dlBLstart(options%ncBLstart), &
            options%dlBLend(options%ncBLend))

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
        integer(I8), allocatable    :: vertdel(:), vp1(:), vp2(:)

        ! Loop
        integer(I8)                 :: i, k

        ! Remove
        !=======
        allocate(vertdel(0))
        do i = 1, size(delind)
            lb = grid%cell%vp1%Get(delind(i))
            ub = lb  + grid%cell%vp2%Get(delind(i))-1
            vertdel = [vertdel, [(k, k = lb, ub)]]
        end do 
        call grid%cell%vert%Remove(vertdel)
        call grid%cell%vp2%Remove(delind)
        call grid%cell%region%Remove(delind)
        grid%cell%ntot = grid%cell%vp2%Size()

        ! Rebuild the pointer
        allocate(vp1(grid%cell%ntot))
        vp1(1) = 1
        vp2 = grid%cell%vp2%Get()
        do i = 2, grid%cell%ntot
            vp1(i) = vp1(i-1) + vp2(i-1)
        end do
        call grid%cell%vp1%Set(vp1)

    end subroutine

    !------------------------------------------------------------------!
    !                      GGTM SEGMENT ROUTINES                       !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeGGTMSegment(segment, xl, yl, fsID, TMfaceID, &
        sv, ev)

        ! Description
        !============
        ! Initialization. We need the original line coordinates, flux
        ! surface ID (can be zero), face ID of the topomesh (can be zero)
        ! and the start and end grid vertex (sv, ev, resp.) for initial
        ! setup. It is assumed that the coordinates of sv, ev correspond
        ! to the first and last coordinate of xl, yl. The vertex 
        ! coordinates are initialized as such. Single-coordinate 
        ! 'segments' are supported to enable 'segments' for special 
        ! topological mesh vertices such as tangency points. In this 
        ! case, the TMfaceID will correspond to the vertex ID in the 
        ! topological mesh. Note that still xl, yl, need to have two 
        ! elements (which should then be the same) and that sv, ev should
        ! have the same value 

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMSegmentUDT)                       :: segment 
        real(R8), intent(in), dimension(:)          :: xl, yl 
        integer(I8), intent(in)                     :: fsID, TMfaceID, &
            sv, ev 

        ! Auxiliary
        
        ! Loop
        integer(I8)                                 :: i 

        ! Initialize
        !===========
        ! Basics
        segment%xl = xl 
        segment%yl = yl 
        segment%fsID = fsID 
        segment%TMfaceID = TMfaceID 
        segment%sv = sv 
        segment%ev = ev 
        segment%nl = size(xl)

        ! Vertex data
        segment%nv = 2_I8
        segment%vert = [sv, ev]
        segment%xv = [xl(1), xl(segment%nl)]
        segment%yv = [yl(1), yl(segment%nl)]
        
        ! Compute
        segment%dll = sqrt((xl(2:) - xl(1:size(xl)-1))**2 &
            + (yl(2:) - yl(1:size(yl)-1))**2)
        segment%dllc = xl 
        segment%dllc = 0.0_R8
        do i = 1, size(xl)-1
            segment%dllc(i+1) = segment%dllc(i) + segment%dll(i)
        end do 
        segment%dlcv = [0.0_R8, segment%dllc(segment%nl)]
        
        ! Check
        segment%isvertex = .false. 
        segment%isclosed = .false.
        if (sv == ev) then 
            if (size(xl) == 2) then 
                if (xl(1) == xl(2) .and. yl(1) == yl(2)) then 
                    segment%isvertex = .true.
                else
                    call gdErrorHandler('InitializeGGTMSegment: ' // & 
                        'start and end vertex are the same, but coordinates are not')
                end if 
            else
                if (xl(1) == xl(size(xl)) .and. (yl(1) == yl(size(yl)))) then 
                    segment%isclosed = .true. 
                else
                    ! Weird
                    call gdErrorHandler('InitializeGGTMSegment: start and ' // & 
                        'end vertex are the same, but coordinates are ' // & 
                        'not. Unexpected')
                end if
            end if 
        end if 

        ! Clean
        call segment%Clean()

        ! Checks
        if (any(segment%dlcv(2:)-segment%dlcv(1:segment%nv-1) < disttol) .and. .not. segment%isvertex) then 
            print *, 'InitializeGGTMSegment: vertices form very small face ' // & 
                '(smaller than predefined disttol)'
        end if 
        if (any(segment%dll < disttol) .and. .not. segment%isvertex) then 
            print *, 'InitializeGGTMSegment: GGTM segment is extremely small ' // & 
                '(smaller than predefined disttol)'
        end if 

    end subroutine

    ! Segment cleaning (removal of almost coinciding points)
    subroutine CleanGGTMSegment(segment)

        ! Description
        !============
        ! This routine removes vertices (except end points) to remove
        ! edges with distance smaller than disttol. To be used 
        ! when constructing a new edge. It is assumed that an initial
        ! length distribution etc is already determined. Note that this
        ! routine only considers the original coordinates xl, yl and that
        ! the vertex distribution is not adjusted!

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMSegmentUDT)               :: segment

        ! Auxiliary
        logical, allocatable, dimension(:)  :: keepvert 

        ! Loop
        integer(I8)                         :: i 

        ! Check
        !======
        ! Hedge for vertex segments
        if (segment%isvertex) then 
            return 
        end if 

        ! Check for almost coinciding vertices
        keepvert = [segment%dll > disttol, .true.]
        if (.not. keepvert(1)) then 
            if (segment%nl > 2) then ! don't do if only two points
                keepvert(2) = .false. 
            end if 
            keepvert(1) = .true. 
        end if 
        if (all(keepvert)) then 
            ! Nothing to do, return
            return
        end if 

        ! Rebuild
        !========
        segment%xl = pack(segment%xl, keepvert)
        segment%yl = pack(segment%yl, keepvert)
        segment%nl = size(segment%xl)

        ! Compute
        segment%dll = sqrt((segment%xl(2:) - segment%xl(1:size(segment%xl)-1))**2 &
            + (segment%yl(2:) - segment%yl(1:size(segment%yl)-1))**2)
        segment%dllc = segment%xl 
        segment%dllc = 0.0_R8
        do i = 1, size(segment%xl)-1
            segment%dllc(i+1) = segment%dllc(i) + segment%dll(i)
        end do 
        segment%dlcv = [0.0_R8, segment%dllc(segment%nl)]

    end subroutine

    ! Vertex coordinate addition
    subroutine AddVerticesGGTMSegment(segment, dlcv, vertID)

        ! Description
        !============
        ! This routine adds vertices xv, yv on the segment by interpolating
        ! the x and y coordinates based on the given length distribution
        ! dlcv. Note that this routine will overwrite any pre-existing 
        ! vertex distribution. One can either provide all vertices (then 
        ! sv and ev should be equal to vertID(1) and vertID(end)) or 
        ! provide the vertices in between the existing vertices (then 
        ! both vertID(1) and vertID(end) should be different)

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMSegmentUDT)               :: segment 
        real(R8), intent(in)                :: dlcv(:)
        integer(I8), intent(in)             :: vertID(:)

        ! Auxiliary

        ! Loop

        ! Sanity checks
        !==============
        if (size(vertID) /= size(dlcv)) then 
            call gdErrorHandler('AddVertexCoordinatesGGTMSegment: ' // & 
                'incompatible input dimensions')
        end if 

        if (size(dlcv) > 1) then
            if (any(abs(dlcv(2:) - dlcv(1:size(dlcv)-1)) < disttol)) then 
                print *, 'AddVerticesGGTMSegment: coinciding vertices'
            end if
        end if

        ! Add
        !====
        ! Check
        if (size(dlcv) == 1) then 
            if (dlcv(1) == 0.0_R8) then 
                ! Check - we should in fact do nothing in this case
                if (segment%sv /= vertID(1) .or. segment%ev /= vertID(1)) then 
                    call gdErrorHandler('AddVerticesGGTMSegment: ' // &
                        'segment consisting of single point should have ' // & 
                        'the same starting and ending vertex')
                end if 
                segment%vert = vertID(1)
                return 
            end if 
        end if 

        ! Hedge for out of bounds points/empty arrays
        if (size(dlcv) == 0) then 
            ! Simply add start and end points
            segment%dlcv = [segment%dlcv(1), segment%dlcv(segment%nv)] 
            segment%vert = [segment%sv, segment%ev]
        elseif (segment%sv /= vertID(1) .and. segment%ev /= vertID(size(vertID))) then 
            ! Only add vertices in between 
            segment%dlcv = [segment%dlcv(1), dlcv, segment%dlcv(segment%nv)] 
            segment%vert = [segment%sv, vertID, segment%ev]
        elseif (segment%sv == vertID(1) .and. segment%ev == vertID(size(vertID))) then 
            ! Add regularly
            segment%dlcv = dlcv
            segment%vert = vertID
        else
            call gdErrorHandler('AddVerticesGGTMSegment: only one start or ' // & 
                'end vertex is the same, unexpected')
        end if 

        segment%nv = size(segment%dlcv)
        where (segment%dlcv >= segment%dllc(segment%nl)) segment%dlcv = segment%dllc(segment%nl)     
        
        if (allocated(segment%xv)) then 
            deallocate(segment%xv, segment%yv)
        end if 
        allocate(segment%xv(segment%nv), segment%yv(segment%nv))
        call Interpolate1D(segment%dlcv, segment%xv, segment%dllc, segment%xl)
        call Interpolate1D(segment%dlcv, segment%yv, segment%dllc, segment%yl)
        if (any(.not. ieee_is_finite(segment%xv))) then 
            print *, 'NaNs detected'
        end if 
        
    end subroutine

    ! Segment flipping routine
    subroutine FlipGGTMSegment(segment)

        ! Description
        !============
        ! This routine reverses the segment direction

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMSegmentUDT)               :: segment 

        ! Auxiliary
        type(GGTMSegmentUDT)                :: temp

        ! Checks
        !=======
        if (segment%isvertex) then 
            ! flipping results in same output
            return
        end if 

        ! Flip
        !=====
        ! Copy 
        temp = segment 

        ! Flip
        segment%xl = temp%xl(temp%nl:1:-1)
        segment%yl = temp%yl(temp%nl:1:-1)
        segment%dll = temp%dll(temp%nl-1:1:-1)
        segment%dllc = temp%dllc(temp%nl) - temp%dllc(temp%nl:1:-1)
        segment%dlcv = temp%dlcv(temp%nv) - temp%dlcv(temp%nv:1:-1)
        segment%xv = temp%xv(temp%nv:1:-1)
        segment%yv = temp%yv(temp%nv:1:-1)
        segment%vert = temp%vert(temp%nv:1:-1)
        segment%sv = temp%ev 
        segment%ev = temp%sv

    end subroutine

    ! Segment getter from vertices
    function GetSegmentIDFromVertices(line, ggtmdata, v1, v2) result(segID)

        ! Description
        !============
        ! Simple function that returns the segment ID (index) from 
        ! ggtmdata that has v1, v2 as either start or end index (if 
        ! multiple, which should in fact not happen, the first 
        ! one encountered is returned). If none is found, segID is zero.

        ! Declare variables
        !==================
        ! Arguments
        type(GGTMFieldlineDataUDT), intent(in)      :: line
        type(GGTMDataUDt), intent(in)               :: ggtmdata
        integer(I8), intent(in)                     :: v1, v2
        integer(I8)                                 :: segID 

        ! Auxiliary

        ! Loop
        integer(I8)                                 :: i 

        ! Determine seg ID
        !=================
        i = 0
        segID = 0
        do while (i < line%ns)
            ! Update counter 
            i = i + 1

            ! Check
            if (any(ggtmdata%seg(line%segID(i))%ev == [v1, v2]) .and. &
                any(ggtmdata%seg(line%segID(i))%sv == [v1, v2])) then 
                segID = line%segID(i)
                exit
            end if 
        end do

    end function

    !------------------------------------------------------------------!
    !                      GGTM LINE HANDLING                          !
    !------------------------------------------------------------------!

    ! GGTM line initialization
    subroutine InitializeGGTMFieldLineData(line, ggtmdata, segID)

        ! Description
        !============
        ! This routine initializes the GGTM field line, based on the 
        ! segments that the field line consists of (these should be 
        ! given in the correct order in segID and should map to the
        ! correct segments in ggtmdata).

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMFieldlineDataUDT)         :: line 
        type(GGTMDataUDT), intent(in)       :: ggtmdata
        integer(I8), intent(in)             :: segID(:)

        ! Auxiliary
        real(R8)                            :: Ltot
        real(R8), allocatable               :: xl(:), yl(:), dlcv(:), &
            xv(:), yv(:)
        integer(I8)                         :: ns
        integer(I8), allocatable            :: fsID(:), nodes(:), &
            TMfaceID(:), vert(:)
        logical, allocatable                :: doflip(:), isnodevert(:)
        type(GGTMSegmentUDT)                :: tseg

        ! Loop
        integer(I8)                         :: i 

        ! Construct
        !==========
        ! Associate
        associate(seg       => ggtmdata%seg)

        ! Initialize
        ns = size(segID)
        allocate(fsID(ns), nodes(ns+1), TMfaceID(ns), doflip(ns), &
            xl(0), yl(0), xv(0), yv(0), dlcv(0), vert(0), isnodevert(0))
        doflip(:) = .false. 

        ! Concatenate segments (need to check in case of multiple segments)
        if (ns > 1) then 
            ! Check that there are no vertex segments, not yet supported
            if (any([seg(segID)%isvertex])) then 
                call gdErrorHandler('InitializeGGTMFieldlineData: ' // & 
                    'field line with multiple segments of which some are ' // & 
                    'vertex segments is not yet supported')
            end if 

            ! Check that there are no closed segments, not yet supported
            if (any([seg(segID)%isclosed])) then 
                call gdErrorHandler('InitializeGGTMFieldlineData: ' // & 
                    'field line with multiple segments of which some are ' // & 
                    'closed segments is not yet supported')
            end if 

            ! Need to check common vertex with next segment
            if (seg(segID(1))%ev == seg(segID(2))%sv .or. &
                seg(segID(1))%ev == seg(segID(2))%ev) then 
                ! End vertex is common vertex and first segment is properly sorted
                doflip(1) = .false.
            elseif (seg(segID(1))%sv == seg(segID(2))%sv .or. &
                seg(segID(1))%sv == seg(segID(2))%ev) then 
                ! Start vertex is common vertex - need to flip first segment
                doflip(1) = .true.
            else ! no common vertex
                call gdErrorHandler('InitializeGGTMFieldLineData: could not ' // &
                    'find common vertex of first segment with next segment')
            end if 
        end if 

        ! Add first segment
        Ltot = 0_R8
        nodes(1) = 1_I8
        isnodevert = [isnodevert, .true.] ! append already first value
        tseg = seg(segID(1))
        if (doflip(1)) then 
            xl = [xl, tseg%xl(tseg%nl:1:-1)]
            yl = [yl, tseg%yl(tseg%nl:1:-1)]
            dlcv = [dlcv, tseg%dlcv(tseg%nv) - tseg%dlcv(tseg%nv:1:-1) + Ltot]
            vert = [vert, tseg%vert(tseg%nv:1:-1)]
        else
            xl = [xl, tseg%xl]
            yl = [yl, tseg%yl]
            dlcv = [dlcv, tseg%dlcv + Ltot]
            vert = [vert, tseg%vert]
        end if 
        Ltot = dlcv(size(dlcv))
        fsID(1) = tseg%fsID
        TMfaceID(1) = tseg%TMfaceID
        nodes(2) = nodes(1) + tseg%nl-1
        isnodevert = [.true., spread(.false., 1, tseg%nv-2), .true.]

        do i = 2, ns
            ! Associate
            associate(tseg    => seg(segID(i)))

            ! Check if we need to flip
            if (seg(segID(i))%sv == seg(segID(i-1))%sv .or. &
                seg(segID(i))%sv == seg(segID(i-1))%ev) then 
                ! Start vertex is common vertex and this segment is properly sorted
                doflip(i) = .false.
            elseif (seg(segID(i))%ev == seg(segID(i-1))%sv .or. &
                seg(segID(i))%ev == seg(segID(i-1))%ev) then 
                ! End vertex is common vertex - need to flip first segment
                doflip(i) = .true.
            else ! no common vertex
                call gdErrorHandler('InitializeGGTMFieldLineData: could not ' // &
                    'find common vertex of current segment with previous segment')
            end if 

            ! Append coordinates etc
            if (doflip(i)) then 
                xl = [xl, tseg%xl(tseg%nl-1:1:-1)]
                yl = [yl, tseg%yl(tseg%nl-1:1:-1)]
                dlcv = [dlcv, tseg%dlcv(tseg%nv) - tseg%dlcv(tseg%nv-1:1:-1) + Ltot]
                vert = [vert, tseg%vert(tseg%nv-1:1:-1)]
            else
                xl = [xl, tseg%xl(2:)]
                yl = [yl, tseg%yl(2:)]
                dlcv = [dlcv, tseg%dlcv(2:) + Ltot]
                vert = [vert, tseg%vert(2:)]
            end if 
            Ltot = dlcv(size(dlcv))
            fsID(i) = tseg%fsID
            TMfaceID(i) = tseg%TMfaceID
            nodes(i+1) = nodes(i) + tseg%nl-1
            isnodevert = [isnodevert, spread(.false., 1, tseg%nv-2), .true.]

            ! Housekeeping
            end associate
        end do 

        ! Compute
        line%xl = xl 
        line%yl = yl 
        line%dll = sqrt((xl(2:) - xl(1:size(xl)-1))**2 &
            + (yl(2:) - yl(1:size(yl)-1))**2)
        line%dllc = xl 
        line%dllc = 0_R8
        do i = 1, size(xl)-1
            line%dllc(i+1) = line%dllc(i) + line%dll(i)
        end do 
        line%dlcv = dlcv

        ! Add
        line%nl         = size(xl)
        line%nv         = size(xv)
        line%ns         = size(fsID)
        line%fsID       = fsID
        line%TMfaceID   = TMfaceID
        line%segID      = segID
        line%flipseg    = doflip
        line%isnodevert = isnodevert
        line%nodes      = nodes

        ! Add vertices
        call line%AddVertexCoordinates(dlcv)
        call line%AddVertexIDs(vert, isnodevert)
                
        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! GGTM segment appending at start/end
    subroutine AppendGGTMFieldLineSegment(line, segmentID, ggtmdata, back)
        
        ! Description
        !============
        ! This routine appends a line segment, given in 'segment', to the
        ! current field line given in 'line'. 'back' should be either 
        ! false (appending at start) or true (appending at end).
        ! Insertion of lines is not (yet) supported. Note that both
        ! lines should already have been initialized. 
        ! 
        ! To merge both segments, we basically do a rebuild of the 
        ! field line from scratch with updated coordinates etc. 

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMFieldlineDataUDT)                 :: line 
        type(GGTMDataUDT), intent(in)               :: ggtmdata
        logical, intent(in)                         :: back 
        integer(I8), intent(in)                     :: segmentID

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: newsegID

        ! Loop


        ! Checks
        !=======
        if (segmentID <= 0) then 
            call gdErrorHandler('AppendGGTMFieldLineSegment: ' // &
                'segment ID must be > 0')
        end if 

        ! Vertex segments are not supported
        if (ggtmdata%seg(segmentID)%isvertex) then 
            call gdErrorHandler('AppendGGTMFieldLineSegment: ' // &
                'segment is vertex segment, appending of these ' // &
                'segments is not yet supported')
        end if 

        ! Though technically possible here, we do not allow extension 
        ! with closed segments (will error in line rebuilding stage)
        if (ggtmdata%seg(segmentID)%isclosed) then 
            call gdErrorHandler('AppendGGTMFieldLineSegment: ' // &
                'segment is closed segment, appending of these ' // &
                'segments is not yet supported')
        end if 

        ! Append segment
        if (back) then 
            newsegID = [line%segID, segmentID]
        else
            newsegID = [segmentID, line%segID]
        end if 

        ! Rebuild line
        call line%Initialize(ggtmdata, newsegID)

    end subroutine

    ! GGTM line flipping
    subroutine FlipGGTMFieldLine(line)

        ! Description
        !============
        ! This subroutine flips the line and all its data

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMFieldlineDataUDT)               :: line 

        ! Auxiliary
        type(GGTMFieldlineDataUDT)                :: temp

        ! Flip
        !=====
        ! Copy 
        temp = line 

        ! Flip
        line%xl         = temp%xl(temp%nl:1:-1)
        line%yl         = temp%yl(temp%nl:1:-1)
        line%dll        = temp%dll(temp%nl-1:1:-1)
        line%dllc       = temp%dllc(temp%nl) - temp%dllc(temp%nl:1:-1)
        line%dlcv       = temp%dlcv(temp%nv) - temp%dlcv(temp%nv:1:-1)
        line%xv         = temp%xv(temp%nv:1:-1)
        line%yv         = temp%yv(temp%nv:1:-1)
        line%vert       = temp%vert(temp%nv:1:-1)
        line%isnodevert = temp%isnodevert(temp%nv:1:-1)
        line%nodes      = temp%nodes(temp%ns+1) - temp%nodes(temp%ns+1:1:-1) + 1
        line%segID      = temp%segID(temp%ns:1:-1)
        line%flipseg    = .not. temp%flipseg(temp%ns:1:-1) 

        ! Face labels may not be allocated/up to date
        if (allocated(line%facelabels)) then 
            if (size(line%facelabels) == temp%nv-1) then 
                line%facelabels = temp%facelabels(temp%nv-1:1:-1)
            end if
        end if
        

    end subroutine

    ! GGTM line splitting
    subroutine SplitGGTMFieldLineAtVertex(line, vID, ggtmdata)

        ! Description
        !============
        ! This routines splits a line at a vertex into (at least) two 
        ! segments. If the vertex ID is a node vertex, nothing is done.
        ! Otherwise, two new segments are constructed and added to the
        ! line. No other operations on the line itself are performed, 
        ! except for adding the segments.  

        ! Note: old segments are not yet deleted nor are other lines
        ! updated! Perhaps need to track this in the future...

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMFieldLineDataUDT)             :: line
        integer(I8), intent(in)                 :: vID 
        type(GGTMDataUDT), intent(inout)        :: ggtmdata

        ! Auxiliary
        integer(I8)                             :: vind, tsegID, segvind(1:2), &
            vindseg, indl1, indl2, tsegIDind, sv, ev
        integer(I8), allocatable, dimension(:)  :: tvertID, newsegID
        real(R8), allocatable, dimension(:)     :: xl, yl, tdlcv, &
            txl, tyl
        type(GGTMSegmentUDT)                    :: tempseg(1:2), tseg
        logical                                 :: isvertonnode 

        ! Loop
        integer(I8)                             :: i 

        ! Checks
        !=======
        ! Initialize
        isvertonnode = .false. 
        sv = line%vert(1)
        ev = line%vert(line%nv)

        ! Hedge for vertex lines
        if (ggtmdata%seg(line%segID(1))%isvertex) then 
            ! Print warning and return
            print *, 'SplitGGTMFieldLineAtVertex: line to be split is ' // & 
                'a vertex, not splitting and continuing...'
            return 
        end if 

        ! Determine vertex index
        vind = findloc(line%vert, vID, 1)
        if (vind == 0) then 
            call gdErrorHandler('SplitGGTMFieldLineAtVertex: line does ' // & 
                'not contain specified vertex, check input')
        end if 

        ! Hedge for vertex node
        if (line%isnodevert(vind)) then 
            ! Simply return
            return 
        end if 

        ! Get the segment that is being split
        tsegID = 0
        do i = 1, line%ns 
            ! Get segment vertex indices
            segvind = line%GetSegmentVertIndices(i)

            ! Check
            if (vind >= segvind(1) .and. vind <= segvind(2)) then 
                tsegIDind = i
                tsegID = line%segID(i)
                exit
            end if 
        end do 

        ! Sanity check - should actually not happen
        if (tsegID == 0) then 
            call gdErrorHandler('SplitGGTMFieldLineAtVertex: could not ' // & 
                'find segment, unexpected')
        end if 

        ! Construct segments
        !===================
        ! Initialize segments
        tseg = ggtmdata%seg(tsegID)
        ggtmdata%nseg = ggtmdata%nseg + 2
        ggtmdata%seg = [ggtmdata%seg, tempseg]

        ! Extract coordinates etc
        xl = tseg%xl
        yl = tseg%yl

        ! Check where the vertex is located in this segment
        vindseg = findloc(tseg%vert, vID, 1)
        if (vindseg == 0) then 
            call gdErrorHandler('SplitGGTMFieldLineAtVertex: segment does ' // & 
                'not contain specified vertex, check input')
        end if 

        ! Check if the vertex lies on a node (approx.)
        indl1 = findloc(abs(tseg%dllc - tseg%dlcv(vindseg)) < disttol, .true., 1)
        if (indl1 /= 0) then 
            isvertonnode = .true. 
            indl2 = indl1
        else
            indl1 = findloc(tseg%dllc < tseg%dlcv(vindseg), .true., 1, back=.true.)
            indl2 = findloc(tseg%dllc > tseg%dlcv(vindseg), .true., 1, back=.false.)
        end if 

        ! Sanity check
        if (indl1 == 0 .or. indl2 == 0) then 
            call gdErrorHandler('SplitGGTMFieldLineAtVertex: segment does ' // & 
                'not contain dlcv specified vertex, unexpected')
        end if 

        ! Construct first segment 
        !------------------------
        ! Extract coordinates
        txl = xl(1:indl1)
        tyl = yl(1:indl1)
        if (.not. isvertonnode) then 
            ! Append the split vertex coordinates
            txl = [txl, line%xv(vind)]
            tyl = [tyl, line%yv(vind)]
        end if 

        ! Get segment data and initialize
        call ggtmdata%seg(ggtmdata%nseg-1)%Initialize(txl, tyl, &
            tseg%fsID, tseg%TMfaceID, line%vert(segvind(1)), line%vert(vind))

        ! Initialize vertices
        tdlcv = line%dlcv(segvind(1)+1:vind-1) - line%dlcv(segvind(1)) ! exclude end vertices
        tvertID = line%vert(segvind(1)+1:vind-1)
        call ggtmdata%seg(ggtmdata%nseg-1)%AddVertices(tdlcv, tvertID)

        ! Construct second segment 
        !-------------------------
        ! Extract coordinates
        txl = xl(indl2:)
        tyl = yl(indl2:)
        if (.not. isvertonnode) then 
            ! Prepend the split vertex coordinates
            txl = [line%xv(vind), txl]
            tyl = [line%yv(vind), tyl]
        end if 

        ! Get segment data and initialize
        call ggtmdata%seg(ggtmdata%nseg)%Initialize(txl, tyl, &
            tseg%fsID, tseg%TMfaceID, line%vert(vind), line%vert(segvind(2)))

        ! Initialize vertices
        tdlcv = line%dlcv(vind+1:segvind(2)-1) - line%dlcv(vind) ! exclude end vertices
        tvertID = line%vert(vind+1:segvind(2)-1)
        call ggtmdata%seg(ggtmdata%nseg)%AddVertices(tdlcv, tvertID)

        ! Reconstruct line
        !=================
        ! Reconstruct based on new segment IDs
        newsegID = [line%segID(1:tsegIDind-1), ggtmdata%nseg-1, ggtmdata%nseg, &
            line%segID(tsegIDind+1:line%ns)]
        call line%Initialize(ggtmdata, newsegID)

        ! Check if we need to flip
        if (line%vert(1) == sv .and. line%vert(line%nv) == ev) then 
            ! All good
        elseif (line%vert(1) == ev .and. line%vert(line%nv) == sv) then 
            ! All good, but flip
            call line%Flip()
        else
            ! All bad
            call gdErrorHandler('splitGGTMFieldLineAtVertex: could not ' // & 
                'find original start and end vertex of line, this ' // & 
                'is a bug')
        end if 

    end subroutine

    ! GGTM line splitting, multiple vertices
    subroutine SplitGGTMFieldLineAtVertices(line, vID, ggtmdata)

        ! Description
        !============
        ! This routines splits a line at a vertex into (at least) two 
        ! segments. If the vertex ID is a node vertex, nothing is done.
        ! Otherwise, two new segments are constructed and added to the
        ! line. No other operations on the line itself are performed, 
        ! except for adding the segments.  

        ! Note: old segments are not yet deleted nor are other lines
        ! updated! Perhaps need to track this in the future...

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMFieldLineDataUDT)             :: line
        integer(I8), intent(in)                 :: vID(:) 
        type(GGTMDataUDT), intent(inout)        :: ggtmdata

        ! Auxiliary
        integer(I8)                             :: segvind(1:2), &
            indlstart, indlend, sv, ev, nv, startv, endv, &
            nsegvID, vindsegstart, vindsegend, segID
        integer(I8), allocatable, dimension(:)  :: tvertID, newsegID, &
            vind, tvID, segvID, linesegID, linevind
        real(R8), allocatable, dimension(:)     :: xl, yl, tdlcv, &
            txl, tyl
        type(GGTMSegmentUDT)                    :: tseg
        type(GGTMSegmentUDT), allocatable       :: tempseg(:)
        logical                                 :: isnodestart, isnodeend
        logical, allocatable, dimension(:)      :: isvertonnode, wasverttreated, &
            isvertonseg

        ! Loop
        integer(I8)                             :: i, j, k

        ! Checks
        !=======
        ! Initialize
        nv = size(vID)
        sv = line%vert(1)
        ev = line%vert(line%nv)

        ! Update segment data to be sure
        call line%UpdateSegmentData(ggtmdata)

        ! Hedge for vertex lines
        if (ggtmdata%seg(line%segID(1))%isvertex) then 
            ! Print warning and return
            print *, 'SplitGGTMFieldLineAtVertices: line to be split is ' // & 
                'a vertex, not splitting and continuing...'
            return 
        end if 

        ! Hedge for no vertices present
        if (nv == 0) then 
            return 
        end if 

        ! Determine vertex index
        allocate(vind(nv))
        do i = 1, nv
            vind(i) = findloc(line%vert, vID(i), 1)
        end do 
        if (any(vind == 0)) then 
            call gdErrorHandler('SplitGGTMFieldLineAtVertices: line does ' // & 
                'not contain specified vertex, check input')
        end if 

        ! Sort
        call Sort(vind, ascend=.true.)

        !! Check for outermost points
        !if (vind(1) /= 1) then
        !    ! Append
        !    vind = [1, vind]
        !    nv = nv + 1
        !end if 
        !if (vind(nv) /= line%nv) then 
        !    vind = [vind, line%nv]
        !    nv = nv + 1
        !end if

        ! Get vertex IDs
        tvID = line%vert(vind)

        ! Further initialization
        isvertonnode = line%isnodevert(vind) 
        wasverttreated = isvertonnode 
        allocate(newsegID(nv-1))
        newsegID = 0
        allocate(linesegID(0))

        ! Construct segments
        !===================
        ! Loop over existing segments
        do i = 1, line%ns

            ! Determine vertices on segment
            !------------------------------
            ! Get segment
            tseg = ggtmdata%seg(line%segID(i))

            ! Flip if necessary to get it in the same direction as the
            ! line
            if (line%flipseg(i)) then 
                call tseg%Flip()
            end if 

            ! Check which vertices we need to split on this segment
            segvind = line%GetSegmentVertIndices(i)
            isvertonseg = (vind >= segvind(1)) .and. (vind <= segvind(2))
            startv = findloc(isvertonseg , .true., 1, back=.false.)
            endv = findloc(isvertonseg, .true., 1, back=.true.)

            ! Check if there is any vertex on this segment
            if (endv == 0 .and. startv == 0) then 
                ! Add the segment to the list of the face
                linesegID = [linesegID, line%segID(i)]

                ! Skip remainder of the loop
                cycle
            elseif ((endv == 0 .and. startv /= 0) .or. (endv == 0 .and. startv /= 0)) then
                ! Sanity check failed
                call gdErrorHandler('SplitGGTMFieldLineAtVertices: ' // & 
                    'only start or end found, this is a bug ')
            end if 

            ! If endv is one more than startv, and both are node vertices, 
            ! then the segment already exists and does not have to be splitted.
            ! Check if both vertices are already node vertices, then skip
            if (line%isnodevert(vind(startv)) .and. line%isnodevert(vind(endv)) .and. &
                (endv - startv == 1)) then 
                ! Add the segment to the list of the face
                linesegID = [linesegID, line%segID(i)]

                ! Skip remainder of the loop
                cycle
            end if 

            ! Get the current vertices to split
            segvID = tvID(startv:endv)
            nsegvID = size(segvID)
            linevind = vind(startv:endv)

            ! Check if the start and end vertices of the segment are
            ! already included, otherwise append
            if (segvID(1) /= tseg%sv) then 
                segvID = [tseg%sv, segvID]
                nsegvID = nsegvID + 1
                linevind = [findloc(line%vert, tseg%sv, 1, back=.false.), linevind]
                if (linevind(1) == 0) then 
                    call gdErrorHandler('SplitGGTMFieldLineAtVertices: ' // & 
                        'could not find segment vertex in line vertices, this is a bug')
                end if 
            end if 
            if (segvID(nsegvID) /= tseg%ev) then 
                segvID = [segvID, tseg%ev]
                nsegvID = nsegvID + 1
                linevind = [linevind, findloc(line%vert, tseg%ev, 1, back=.false.)]
                if (linevind(size(linevind)) == 0) then 
                    call gdErrorHandler('SplitGGTMFieldLineAtVertices: ' // & 
                        'could not find segment vertex in line vertices, this is a bug')
                end if 
            end if 

            ! Split segment
            !--------------
            ! Initialize segment index
            segID = ggtmdata%nseg 

            ! Increase size of ggtmdata
            allocate(tempseg(nsegvID-1))
            ggtmdata%seg = [ggtmdata%seg, tempseg]
            linesegID = [linesegID, [(k, k = ggtmdata%nseg+1, ggtmdata%nseg+size(tempseg))]]
            ggtmdata%nseg = ggtmdata%nseg + size(tempseg)
            deallocate(tempseg)
            

            ! Extract coordinates
            xl = tseg%xl
            yl = tseg%yl

            ! Loop over all new segments to be constructed
            do j = 1, nsegvID-1

                ! Check where the vertex is located in this segment
                vindsegstart = findloc(tseg%vert, segvID(j), 1)
                vindsegend = findloc(tseg%vert, segvID(j+1), 1)

                ! Sanity check
                if (vindsegstart == 0 .or. vindsegend == 0) then 
                    call gdErrorHandler('SplitGGTMFieldLineAtVertices: segment does ' // & 
                        'not contain specified vertex, check input')
                end if 

                ! Check if the vertices lie on a node (approx.)
                indlstart = findloc(abs(tseg%dllc - tseg%dlcv(vindsegstart)) < disttol, .true., 1)
                indlend = findloc(abs(tseg%dllc - tseg%dlcv(vindsegend)) < disttol, .true., 1)
                isnodestart = .false.
                isnodeend = .false.
                if (indlstart /= 0) then 
                    isnodestart = .true. 
                else
                    indlstart = findloc(tseg%dllc > tseg%dlcv(vindsegstart), .true., 1, back=.false.)
                end if 
                if (indlend /= 0) then 
                    isnodeend = .true.
                else
                    indlend = findloc(tseg%dllc < tseg%dlcv(vindsegend), .true., 1, back=.true.)
                end if 

                ! Sanity check
                if (indlstart == 0 .or. indlend == 0) then 
                    call gdErrorHandler('SplitGGTMFieldLineAtVertices: segment does ' // & 
                        'not contain dlcv specified vertex, unexpected')
                end if 

                ! Construct segment coordinates
                txl = xl(indlstart:indlend)
                tyl = yl(indlstart:indlend)
                if (.not. isnodestart) then 
                    txl = [line%xv(linevind(j)), txl]
                    tyl = [line%yv(linevind(j)), tyl]
                end if 
                if (.not. isnodeend) then 
                    txl = [txl, line%xv(linevind(j+1))]
                    tyl = [tyl, line%yv(linevind(j+1))]
                end if 

                ! Initialize segment coordinates
                segID = segID + 1
                call ggtmdata%seg(segID)%Initialize(txl, tyl, &
                    tseg%fsID, tseg%TMfaceID, segvID(j), segvID(j+1))

                ! Initialize segment vertices
                tdlcv = line%dlcv(linevind(j)+1:linevind(j+1)-1) - line%dlcv(linevind(j))
                tvertID = line%vert(linevind(j)+1:linevind(j+1)-1)
                call ggtmdata%seg(segID)%AddVertices(tdlcv, tvertID)

            end do 
        end do

        ! Reconstruct line
        !=================
        ! Reconstruct based on new segment IDs
        call line%Initialize(ggtmdata, linesegID)

        ! Check if we need to flip
        if (line%vert(1) == sv .and. line%vert(line%nv) == ev) then 
            ! All good
        elseif (line%vert(1) == ev .and. line%vert(line%nv) == sv) then 
            ! All good, but flip
            call line%Flip()
        else
            ! All bad
            call gdErrorHandler('splitGGTMFieldLineAtVertex: could not ' // & 
                'find original start and end vertex of line, this ' // & 
                'is a bug')
        end if 

    end subroutine

    ! GGTM line splitting, multiple nodes of a segment
    subroutine SplitGGTMFieldLineAtSegmentNodes(line, segID, nodeind, &
        ggtmdata, vertID, newvertID)

        ! Description
        !============
        ! This routine splits a line at a node position (so xl, yl 
        ! coordinates) of a certain vertex (i.e. the nodeind is relative
        ! index w.r.t. the nodes of the segment in segID). It is assumed
        ! that this nodeID makes sense. Since splitting at a node will 
        ! likely introduce a new vertex in the grid, the current (maximal)
        ! amount of grid nodes should be passed through 'vertID', which 
        ! will increase in case new vertices are introduced. The new 
        ! vertex IDs, in sequence of the given nodeind, will be returned
        ! in 'newvertID' for reference. 

        ! Note: under the hood, this function simply calls the SplitAtVertices
        ! routine after introducing a new vertex at the segID

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMFieldLineDataUDT)             :: line 
        integer(I8), intent(in)                 :: segID 
        integer(I8), intent(in), dimension(:)   :: nodeind 
        type(GGTMDataUDT), intent(inout)        :: ggtmdata
        integer(I8), intent(inout)              :: vertID 
        integer(I8), intent(out), dimension(:), allocatable  :: newvertID

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: newvert, sortind
        real(R8), allocatable, dimension(:)     :: ltot, newdlcv
        logical, allocatable, dimension(:)      :: isnodevert

        ! Loop 
        integer(I8)                             :: i, k 

        ! Insert vertex
        !==============
        ! Check
        if (.not. any(line%segID == segID)) then 
            call gdErrorHandler('SplitGGTMFieldLineAtSegmentNodes: ' // & 
                'segment is not part of the line')
        end if

        ! Determine location in line
        ltot = ggtmdata%seg(segID)%dllc(nodeind)
        i = 1
        do while (line%segID(i) /= segID)
            ! Add segment length
            ltot = ltot + ggtmdata%seg(line%segID(i))%dllc(ggtmdata%seg(line%segID(i))%nl)

            ! Update
            i = i + 1
        end do 

        ! Construct new vertices
        newvertID = [(k, k = vertID+1, vertID+size(nodeind))]
        newdlcv = [line%dlcv, ltot]
        isnodevert = [line%isnodevert, spread(.false., 1, size(nodeind))] ! these should not lie on nodes
        newvert = [line%vert, newvertID]
        vertID = vertID + size(nodeind)

        ! Sort
        allocate(sortind(size(newdlcv)))
        call Sort(newdlcv, ind=sortind)
        newvert = newvert(sortind)

        ! Add
        call line%AddVertexCoordinates(newdlcv)
        call line%AddVertexIDs(newvert, isnodevert)
        call line%UpdateSegmentData(ggtmdata)

        ! Split at vertices
        !==================
        call line%SplitAtVertices(newvertID, ggtmdata)

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
        if (all(tdlcv == 0.0_R8)) then 
            line%dlcv = tdlcv 
            line%xv = spread(line%xl(1), 1, size(tdlcv))
            line%yv = spread(line%yl(1), 1, size(tdlcv))
            line%nv = size(tdlcv)
            return 
        end if 

        ! Hedge for out of bounds points
        line%dlcv = tdlcv 
        if (count(line%dlcv >= line%dllc(line%nl)) > 1) then 
            print *, 'weird'
        end if
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
    subroutine AddVertexIDs(line, vertIDs, isnodevert)

        ! Description
        !============
        ! ID setter (simple wrapper), also sets the logical to check
        ! if a vertex is a node vertex (we assume this is correctly
        ! done... )

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMFieldlineDataUDT)     :: line 
        integer(I8), intent(in)         :: vertIDs(:)
        logical, intent(in)             :: isnodevert(:)

        ! Set
        !====
        if (size(vertIDs) /= size(isnodevert)) then 
            call gdErrorHandler('AddVertexIDs: incompatible input sizes')
        end if 
        line%vert = vertIDs 
        line%isnodevert = isnodevert

    end subroutine

    ! GGTM segment data updater
    subroutine UpdateSegmentData(line, ggtmdata, doseg)

        ! Description
        !============
        ! This routine updates the segment data stored in ggtmdata by
        ! updating the vertices of these segments to correspond to the 
        ! new vertices in a line (since typically we work line-based for
        ! refinement etc, but segment-based when constructing flux tubes
        ! etc since the line-based framework then breaks down)

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMFieldlineDataUDT)                 :: line
        type(GGTMDataUDT), intent(inout)            :: ggtmdata

        ! Auxiliary
        integer(I8)                                 :: tvind(1:2), segID
        integer(I8), allocatable, dimension(:)      :: vertID
        real(R8), allocatable, dimension(:)         :: tdlcv
        logical, allocatable, dimension(:)          :: adjustseg
        logical, allocatable, dimension(:), optional, intent(in)    :: doseg

        ! Loop
        integer(I8)                                 :: i 

        ! Initialize
        !===========
        ! Associate
        associate(&
            seg       => ggtmdata%seg     &
            )

        ! Check
        if (present(doseg)) then 
            adjustseg = doseg 
            if (size(adjustseg) /= line%ns) then 
                call gdErrorHandler('UpdateSegmentData: incompatible ' // &
                    'dimension of optional argument doseg')
            end if 
        else
            allocate(adjustseg(line%ns))
            adjustseg = .true.
        end if 


        ! Update vertices
        !================
        ! Loop over all segments
        do i = 1, size(line%segID)
            ! Get the vertices of this segment 
            tvind = line%GetSegmentVertIndices(i)
            segID = line%segID(i)

            ! Normally, the end vertices should equal the start and 
            ! end vertices of the segment - we also need to check if we 
            ! should flip
            vertID = line%vert(tvind(1):tvind(2))
            tdlcv = line%dlcv(tvind(1):tvind(2))
            tdlcv = tdlcv - line%dllc(line%nodes(i)) ! subtract length of previous segments!
            if (any(ubound(tdlcv) < 1)) then 
                print *, 'weird'
            end if
            tdlcv(1) = 0.0_R8
            tdlcv(size(tdlcv)) = seg(segID)%dllc(seg(segID)%nl)
            if (seg(segID)%sv == vertID(1) .and. seg(segID)%ev == vertID(size(vertID))) then 
                ! Add regularly, nothing to do here
            elseif (seg(segID)%ev == vertID(1) .and. seg(segID)%sv == vertID(size(vertID))) then
                ! Need to flip
                vertID = vertID(size(vertID):1:-1)
                tdlcv = tdlcv(size(tdlcv)) - tdlcv(size(tdlcv):1:-1)
            else
                ! Shouldn't happen
                call gdErrorHandler('UpdateSegmentData: start and end vertex ' // & 
                    'of segment not found')
            end if 

            ! Add vertices 
            if (adjustseg(i)) then 
                call seg(segID)%AddVertices(tdlcv, vertID)
            endif
        end do 
        
        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! GGTM line data updater
    subroutine UpdateLineData(line, ggtmdata)

        ! Description
        !============
        ! This routine updates the line data starting from the segment 
        ! data that this line has. It is assumed that all fields etc are
        ! properly initialized. Gridding data is not yet updated in this
        ! routine! Basically, these are only the vertices since the 
        ! underlying coordinates themselves (and orientation) should not
        ! change anymore. If this would be the case, then one needs to 
        ! fully reinitialize the line from scratch. 

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMFieldLineDataUDT)         :: line 
        type(GGTMDataUDT), intent(in)       :: ggtmdata
        
        ! Auxiliary
        real(R8)                            :: Ltot
        real(R8), allocatable               :: dlcv(:), &
            xv(:), yv(:)
        integer(I8), allocatable            :: vert(:)
        logical, allocatable                :: isnodevert(:)
        type(GGTMSegmentUDT)                :: tseg

        ! Loop
        integer(I8)                         :: i, cc

        ! Construct
        !==========
        ! Associate
        associate(seg       => ggtmdata%seg)

        ! Initialize
        allocate(xv(0), yv(0), dlcv(0), vert(0), isnodevert(0))

        ! Concatenate
        Ltot = 0_R8
        do i = 1, line%ns
            ! Get segment
            tseg = seg(line%segID(i))

            ! Check if we need to flip it
            if (line%flipseg(i)) then 
                call tseg%Flip()
            end if 

            ! Append
            if (i == 1) then 
                xv = [xv, tseg%xv]
                yv = [yv, tseg%yv]
                vert = [vert, tseg%vert]
                isnodevert = [.true., spread(.false., 1, tseg%nv-2), .true.]
                dlcv = tseg%dlcv + Ltot
            else
                xv = [xv, tseg%xv(2:)]
                yv = [yv, tseg%yv(2:)]
                vert = [vert, tseg%vert(2:)]
                isnodevert = [isnodevert, spread(.false., 1, tseg%nv-2), .true.]
                dlcv = [dlcv, tseg%dlcv(2:) + Ltot]
            end if 
            Ltot = dlcv(size(dlcv))
        end do 

        ! Ensure node vertices lie exactly on the node coordinates
        cc = 1
        do i = 1, size(dlcv)
            if (i > size(isnodevert)) then 
                print *, 'weird'
            end if 
            if (isnodevert(i)) then 
                dlcv(i) = line%dllc(line%nodes(cc))
                cc = cc + 1
            end if 
        end do 

        ! Add vertices
        call line%AddVertexCoordinates(dlcv)
        call line%AddVertexIDs(vert, isnodevert)

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! GGTM line gridding data updater
    subroutine UpdateLineGriddingData(line, ggtmdata)

        ! Description
        !============
        ! This routine updates additional line data required for later
        ! gridding. This includes currently only face labels (vertices 
        ! etc are assumed to be present already)

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMFieldlineDataUDT)         :: line 
        class(GGTMDataUDT), intent(in)      :: ggtmdata

        ! Auxiliary
        integer(I8)                             :: sfind(1:2)

        ! Loop
        integer(I8)                         :: i

        ! Update
        !=======
        ! Initialize
        line%facelabels = line%vert(1:line%nv-1)
        line%facelabels = 0_I8

        ! Loop over all segments
        do  i = 1, line%ns
            ! Get the segment face indices
            sfind = line%GetSegmentFaceIndices(i)

            ! Set labels
            line%facelabels(sfind(1):sfind(2)) = line%TMfaceID(i)
        end do

    end subroutine

    ! GGTM field line segment vertex indices getter
    function GetGGTMFieldLineSegmentVertIndices(line, i) result(vertind)

        ! Description
        !============
        ! This function returns the vertex indices (going from 1 to nv-1)
        ! of vertices that lie on the i-th segment. We do not perform 
        ! sanity checks here, it is assumed that i is not out of bounds.
        ! The result is an array of size 2 that holds the start and end
        ! indices for the face indices (i.e. vertind(1):vertind(2) 
        ! gives an array slice equal to the face indices). It is possible
        ! to do this this way since the vertices etc should be contiguous. 

        ! Note: is it assumed that vertices etc have been properly 
        ! initialized and constructed

        ! Note: we use minloc instead of bitwise comparison since small
        ! numerical roundoff errors may occur
        
        ! Declare variables
        !==================
        ! Arguments
        class(GGTMFieldlineDataUDT), intent(in)     :: line 
        integer(I8), intent(in)                     :: i 
        integer(I8)                                 :: vertind(1:2)

        ! Determine face index bounds
        !============================
        vertind = 0 
        vertind(1) = minloc(abs(line%dlcv - line%dllc(line%nodes(i))), 1)
        vertind(2) = minloc(abs(line%dlcv - line%dllc(line%nodes(i+1))), 1)
        if (any(vertind == 0)) then 
            call gdErrorHandler('GetGGTMFieldLineSegmentVertIndices: ' // &
                'could not find start or end of vertex indices, check input')
        elseif (vertind(2) < vertind(1)) then 
            ! Wrongly sorted
            call gdErrorHandler('GetGGTMFieldLineSegmentVertIndices: ' // &
                'vertex index 2 is smaller than vertex index 1, line ' // & 
                'seems not properly sorted')
        end if 

    end function 

    ! GGTM field line all segment vertex indices getter
    function GetGGTMFieldLineAllSegmentVertIndices(line) result(vertind)

        ! Description
        !============
        ! Same as for one segment, but now all (unique and sorted) 
        ! vertex IDs are returned. 
        
        ! Declare variables
        !==================
        ! Arguments
        class(GGTMFieldlineDataUDT), intent(in)     :: line 
        integer(I8), allocatable, dimension(:)      :: vertind(:)

        ! Auxiliary
        integer(I8), allocatable, dimension(:, :)   :: tempvertind(:, :)

        ! Loop
        integer(I8)                                 :: i 

        ! Compute
        !========
        ! Initialize
        allocate(tempvertind(line%ns, 2))

        ! Determine all vertex indices
        do i = 1, line%ns
            tempvertind(i, :) = line%GetSegmentVertIndices(i)
        end do 

        ! Get unique set
        call Unique(reshape(tempvertind, [size(tempvertind)]), vertind)

    end function 

    ! GGTM field line segment face indices getter
    function GetGGTMFieldLineSegmentFaceIndices(line, i) result(faceind)

        ! Description
        !============
        ! This function returns the face indices (going from 1 to nv-1)
        ! of faces that lie on the i-th segment. We do not perform 
        ! sanity checks here, it is assumed that i is not out of bounds.
        ! The result is an array of size 2 that holds the start and end
        ! indices for the face indices (i.e. faceind(1):faceind(2) 
        ! gives an array slice equal to the face indices). It is possible
        ! to do this this way since the faces etc should be contiguous. 

        ! Note: is it assumed that vertices etc have been properly 
        ! initialized and constructed
        
        ! Declare variables
        !==================
        ! Arguments
        class(GGTMFieldlineDataUDT), intent(in)     :: line 
        integer(I8), intent(in)                     :: i 
        integer(I8)                                 :: faceind(1:2)

        ! Determine face index bounds
        !============================
        faceind = line%GetSegmentVertIndices(i)
        faceind(2) = faceind(2) - 1

    end function

    ! GGTM line pair initialization
    subroutine InitializeGGTMFieldlinePairData(linepair)

        ! Description
        !============
        ! Initialize the line pair data 

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMFieldlinePairDataUDT)         :: linepair 

        ! Initialize
        !===========
        ! Check allocation status
        if (allocated(linepair%l1minLOS)) then 
            deallocate(linepair%l1minLOS)
        end if 
        if (allocated(linepair%l2minLOS)) then 
            deallocate(linepair%l2minLOS)
        end if 
        if (allocated(linepair%l1maxLOS)) then 
            deallocate(linepair%l1maxLOS)
        end if 
        if (allocated(linepair%l2maxLOS)) then 
            deallocate(linepair%l2maxLOS)
        end if 

        ! Allocate and initialize
        allocate(linepair%l1minLOS(0), linepair%l1maxLOS(0), &
            linepair%l2minLOS(0), linepair%l2maxLOS(0))
        linepair%l1minLOS = 0_I8
        linepair%l1maxLOS = 0_I8
        linepair%l2minLOS = 0_I8
        linepair%l2maxLOS = 0_I8
        linepair%srflabel = 0_I8
        linepair%erflabel = 0_I8
        linepair%isextendedstart = .false. 
        linepair%isextendedend = .false.
        linepair%hfface     = ConstructIntegerDynamicArray()
        linepair%lfface     = ConstructIntegerDynamicArray()
        linepair%tubeface   = ConstructIntegerDynamicArray()
        linepair%cell       = ConstructIntegerDynamicArray()

    end subroutine

    ! GGTM line pair graph initialization 
    subroutine InitializeGGTMFieldlinePairDataGraph(linepair)

        ! Description
        !============
        ! This routine initializes the graph of the field line pair 
        ! data. This is a separate routine, since this initialization 
        ! is non-trivial and involves some rather costly operations.
        ! It is likely also not necessary in most cases to have this
        ! graph if the grid is 'nicely' behaved (i.e. we get away 
        ! with just relying on the magnetic field and other subroutines
        ! to give a non-overlapping grid). However, in several cases it
        ! may be that cell overlap will occur, and that one may need to 
        ! properly check which vertex can be connected where. This is 
        ! where the graph comes in - it can determine which connections
        ! are possible etc. This is, however, a rather expensive 
        ! operation (not only to construct, but also to evaluate).

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMFieldlinePairDataUDT)         :: linepair 

        ! Auxiliary
        logical                                 :: isintube
        integer(I8)                             :: tvind
        integer(I8), allocatable, dimension(:)  :: ev1, ev2, &
            gvert, hfv, lfv, vp
        real(R8), allocatable, dimension(:)     :: xp, yp, xv, yv, &
            xvert, yvert

        ! Loop
        integer(I8)                             :: i, j, k, svhf, svlf, &
            evhf, evlf

        ! Initialize
        !===========
        ! Unpack
        associate(&
            graph       => linepair%graph,  &
            hfline      => linepair%hfline, &
            lfline      => linepair%lfline  &
            )        

        ! Construct closed polygon coordinates and mapping from hfline
        ! and lfline indices
        xp = hfline%xv
        yp = hfline%yv 
        vp = hfline%vert
        hfv = [(k, k = 1, hfline%nv)]
        lfv = [(k, k = lfline%nv, 1, -1)] + hfline%nv
        if (hfline%vert(hfline%nv) /= lfline%vert(lfline%nv)) then 
            xp = [xp, lfline%xv(lfline%nv:1:-1)]
            yp = [yp, lfline%yv(lfline%nv:1:-1)]
            vp = [vp, lfline%vert(lfline%nv:1:-1)]
        else
            xp = [xp, lfline%xv(lfline%nv-1:1:-1)]
            yp = [yp, lfline%yv(lfline%nv-1:1:-1)]
            vp = [vp, lfline%vert(lfline%nv-1:1:-1)]
            lfv = lfv - 1
        end if 
        if (hfline%vert(1) /= lfline%vert(1)) then 
            xp = [xp, hfline%xv(1)]
            yp = [yp, hfline%yv(1)]
            vp = [vp, hfline%vert(1)]
        end if 

        ! Compute vertices
        !=================
        ! Easy - just take vertices of hfline and lfline 
        gvert = [hfline%vert, lfline%vert]
        xvert = [hfline%xv, lfline%xv]
        yvert = [hfline%yv, lfline%yv]

        ! Compute edges
        !==============
        ! Hard(er) - need to check:
        ! - which edges cross 
        ! - which edges (that don't cross) go out of bounds
        ! Using the dedicated routine from mod_polygon

        ! Initialize
        allocate(ev1(0), ev2(0))

        ! Hedge for the same start/end vertex
        svhf = 1
        svlf = 1
        evhf = hfline%nv
        evlf = lfline%nv
        if (hfline%vert(1) == lfline%vert(1)) then 
            svhf = 2
            svlf = 2
        end if 
        if (hfline%vert(hfline%nv) == lfline%vert(lfline%nv)) then 
            evhf = evhf - 1
            evlf = evlf - 1
        end if 

        ! Loop over hfline vertices
        do i = svhf, evhf
            ! Unpack
            associate(&
                v1 => hfline%vert(i),    &
                x1 => hfline%xv(i),      &
                y1 => hfline%yv(i)       &
                )

            ! Loop over lfline vertices
            do j = svlf, evlf 
                ! Unpack 
                associate(&
                    v2 => lfline%vert(j),    &
                    x2 => lfline%xv(j),      &
                    y2 => lfline%yv(j)       &
                    )

                ! Skip check
                !===========
                ! Add if first or last radial face (if it exists) - 
                ! should always be present 
                if ((i == 1 .and. j == 1) .or. &
                    (i == hfline%nv .and. j == lfline%nv)) then 
                    ev1 = [ev1, v1]
                    ev2 = [ev2, v2]
                    cycle 
                end if 

                ! Edge in polygon check
                !======================
                ! Check if edge starts and ends in the interior 
                isintube = IsEdgeInClosedSimplePolygon(xp, yp, vp, hfv(i), &
                    lfv(j))
                if (.not. isintube) then 
                    cycle
                end if 

                ! If we got here, add the vertices
                ev1 = [ev1, v1]
                ev2 = [ev2, v2]

                ! Housekeeping
                end associate
            end do 
            ! Housekeeping
            end associate
        end do

        ! Construct graph
        !================
        ! Construct
        call graph%Construct(ev1, ev2, gvert)

        ! Get vertex coordinates
        allocate(xv(graph%nv), yv(graph%nv))
        do i = 1, size(gvert)
            tvind = graph%GetVertexIndex(gvert(i))
            xv(tvind) = xvert(i)
            yv(tvind) = yvert(i)
        end do 
        linepair%graphxv = xv 
        linepair%graphyv = yv
        
        ! Write data (to be removed)
        ! call linepair%VisualizeGraph('lpgraph')

        ! Checks? Probably not needed here as this should trigger when 
        ! forming cells/faces etc

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! GGTM line pair graph visualization 
    subroutine VisualizeGGTMFieldlinePairDataGraph(linepair, filename)

        ! Description
        !============
        ! This routine visualizes the field line pair graph by writing
        ! out the data in a .dat file with the name specified in 
        ! 'filename'. This file contains the set of edge coordinates of 
        ! the graph. Basically, we construct a polygonset from the 
        ! edges of the graph and visualize that one.

        ! Declare variables
        !=================
        ! Arguments
        class(GGTMFieldlinePairDataUDT)         :: linepair 
        character(*), intent(in)                :: filename

        ! Auxiliary
        integer(I8), allocatable, dimension(:, :)   :: edges
        type(PolygonSetUDT)     :: tps 

        ! Loop 
        integer(I8)                             :: i 

        ! Construct polygonset
        !=====================
        ! Associate
        associate(graph => linepair%graph)

        ! Construct edges
        allocate(edges(graph%ne, 2))
        do i = 1, graph%ne
            edges(i, 1) = graph%ev1(i)
            edges(i, 2) = graph%ev2(i)
        end do 

        ! Construct polygonset
        call tps%Construct(edges, linepair%graphxv, linepair%graphyv)

        ! Visualize
        !==========
        ! Print edges
        call tps%WriteData(filename // '_edges')

        ! Print vertex coordinates
        call Write2DCoordinateData(linepair%graphxv, linepair%graphyv, filename // '_vert')

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! GGTM segment intersections
    subroutine GGTMSegmentIntersections(l1, l2, xint, yint, s1, s2, s1r, s2r, &
        vertbased)

        ! Description
        !============
        ! More or less a wrapper for intersections between GGTM segments.
        ! However, intersections at end points are ignored if the vertex
        ! IDs are the same there. Intersections are computed based on xl
        ! and yl, not on the vertex coordinates (reason why we can check
        ! vertices, is because by definition the first and last vertex
        ! have to lie on the first and last point of the line). 

        ! Note: this routine is particularly useful when checking
        ! intersections between field lines to determine to remove a tube
        ! or field line. 

        ! Declare variables
        !==================
        ! Arguments
        type(GGTMSegmentUDT), intent(in)      :: l1, l2 
        real(R8), allocatable, dimension(:), intent(out)    :: xint, yint
        integer(I8), allocatable, dimension(:), intent(out) :: s1, s2 
        real(R8), allocatable, dimension(:), intent(out), optional  :: s1r, s2r 
        logical, intent(in), optional                       :: vertbased

        ! Auxiliary
        real(R8), allocatable, dimension(:)         :: s1raux, s2raux
        logical, allocatable, dimension(:)          :: keepx

        ! Compute 
        !========
        ! Intersections
        if (present(vertbased)) then 
            if (vertbased) then
                call SimplePolygonIntersections(l1%xv, l1%yv, l2%xv, l2%yv, &
                    xint, yint, s1, s2, s1raux, s2raux) 
            else
                call SimplePolygonIntersections(l1%xl, l1%yl, l2%xl, l2%yl, &
                    xint, yint, s1, s2, s1raux, s2raux)
            end if
        else
            call SimplePolygonIntersections(l1%xl, l1%yl, l2%xl, l2%yl, &
                xint, yint, s1, s2, s1raux, s2raux)
        end if
        

        ! Process
        if (size(xint) > 0) then
            ! First, hedge for intersections exactly 
            ! in end point ->  simply cut cell tube
            allocate(keepx(size(xint)))
            keepx = .true. 
            if (any(l1%vert(1) == [l2%vert(1), l2%vert(l2%nv)])) then 
                where (xint == l1%xl(1) .and. &
                    yint == l1%yl(1)) keepx = .false. 
            end if 
            if (any(l1%vert(l1%nv) == [l2%vert(1), l2%vert(l2%nv)])) then 
                where (xint == l1%xl(l1%nl) .and. &
                    yint == l1%yl(l1%nl)) keepx = .false. 
            end if 

            ! Remove intersections
            xint = pack(xint, keepx)
            yint = pack(yint, keepx)
            s1 = pack(s1, keepx)
            s2 = pack(s2, keepx)
            s1raux = pack(s1raux, keepx)
            s2raux = pack(s2raux, keepx)
        end if 

        ! Check optional arguments
        if (present(s1r) .and. present(s2r)) then 
            s1r = s1raux 
            s2r = s2raux 
        end if 

    end subroutine

    ! GGTM line intersections
    subroutine GGTMLineIntersections(ggtmdata, l1, l2, xint, yint, s1, s2, s1r, s2r, &
        vertbased)

        ! Description
        !============
        ! More or less a wrapper for intersections between GGTM lines.
        ! However, intersections at segment end points are ignored if the vertex
        ! IDs are the same there. Intersections are computed based on xl
        ! and yl, not on the vertex coordinates (reason why we can check
        ! vertices, is because by definition the first and last vertex
        ! have to lie on the first and last point of the line). 

        ! Note: this routine is particularly useful when checking
        ! intersections between field lines to determine to remove a tube
        ! or field line. 

        ! Declare variables
        !==================
        ! Arguments
        type(GGTMDataUDT), intent(in)               :: ggtmdata
        type(GGTMFieldLineDataUDT), intent(in)      :: l1, l2 
        real(R8), allocatable, dimension(:), intent(out)    :: xint, yint
        integer(I8), allocatable, dimension(:), intent(out) :: s1, s2 
        real(R8), allocatable, dimension(:), intent(out), optional  :: s1r, s2r 
        logical, intent(in), optional                       :: vertbased

        ! Auxiliary
        real(R8)                                    :: Li, Lj
        integer(I8)                                 :: ni, nj
        real(R8), allocatable, dimension(:)         :: s1raux, s2raux, &
            tempx, tempy, temps1r, temps2r
        integer(I8), allocatable, dimension(:)      :: temps1, temps2

        ! Loop
        integer(I8)                                 :: i, j 

        ! Initialize
        !=========
        allocate(s1(0), s2(0), xint(0), yint(0), s1raux(0), &
            s2raux(0))

        ! Compute 
        !========
        ! Intersections
        Li = 0 ! Keep track of accumulated segment 'length'
        ni = 0
        do i = 1, l1%ns 
            ! Keep track of accumulated 'length' 
            Lj = 0
            nj = 0

            ! Unpack
            associate(seg1      => ggtmdata%seg(l1%segID(i)))
            do j = 1, l2%ns
                ! Unpack
                associate(seg2      => ggtmdata%seg(l2%segID(j)))
                
                ! Compute intersections
                call GGTMSegmentIntersections(seg1, seg2, tempx, tempy, &
                    temps1, temps2, temps1r, temps2r, vertbased)

                ! Append
                if (size(tempx) > 0) then 
                    xint = [xint, tempx]
                    yint = [yint, tempy]
                    if (l1%flipseg(i)) then 
                        s1 = [s1, seg1%nv - temps1 - 1 + ni]
                        s1raux = [s1raux, seg1%nv - 1 - temps1r + Li]
                    else
                        s1 = [s1, temps1 + ni]
                        s1raux = [s1raux, temps1r + Li]
                    end if 
                    if (l2%flipseg(j)) then 
                        s2 = [s2, seg2%nv - temps2 - 1 + nj]
                        s2raux = [s2raux, seg2%nv - 1 - temps2r + Lj]
                    else
                        s2 = [s2, temps2 + nj]
                        s2raux = [s2raux, temps2r + Lj]
                    end if 
                end if

                ! Housekeeping
                Lj = Lj + seg2%nv-1
                nj = nj + seg2%nv-1
                end associate
            end do 

            ! Update length
            Li = Li + seg1%nv-1
            ni = ni + seg1%nv-1
            end associate
        end do 
        
        ! Check optional arguments
        if (present(s1r) .and. present(s2r)) then 
            s1r = s1raux 
            s2r = s2raux 
        end if 

    end subroutine

    ! GGTM tube extension
    subroutine ExtendTubeWithSegment(cellID, tube, segID, ggtmdata, start, &
        vertID, magneticField)

        ! Description
        !============
        ! This routine extends an existing tube either at the start or 
        ! end (determined by 'start') with the segment with ID 'segID'. 
        ! Vertex segments are not allowed. Hereto, we need to check
        ! the geometry of the segment w.r.t. the geometry of the tube 
        ! lines, in order to have a decent grid later on. We distinguish 
        ! some particular cases based on the angle between the last line
        ! edge and the first segment edge:
        ! - if the lf angle is larger than 90° and the hf angle is 
        !   smaller (or vice versa), add the segment to the line with 
        !   the largest angle. If both are exactly 90°, pick one (doesn't)
        !   matter, but shouldn't occur actually, since this will lead to 
        !   collapsing of faces and hence bad grids). 
        ! - if both angles are smaller than 90°, recompute the angles 
        !   based on a straight segment and redo (this should result in 
        !   the first use case). 
        ! - if both angles are larger than 90°, both lines get part of 
        !   the segment.
        ! Note that we can do the checks just based on the sign of the 
        ! dot product.

        ! If a line gets a part of the segment, we may have to split the 
        ! segment:
        ! - if only one line gets a part, we split the segment in two. 
        ! - if both lines get a part, it will either be two or three 
        ! - segments. 
        ! Splitting is based on geometry considerations. Currently, the 
        ! splitting point is determined as the location where the sign 
        ! of the scalar product between segment edge tangent vectors and 
        ! magnetic field vector becomes negative. This might 
        ! lead to bad vertex distributions locally and to cell overlap.
        ! Note that we adjust the sign of the magnetic field to point 
        ! along the hfline in the segment direction. 

        ! Note: currently, we don't properly keep track of all the 
        ! segments that may be generated. Normally, this is also not 
        ! necessary anymore since tube extension is the last step 
        ! we can take to modify the lines on which vertices are 
        ! distributed. Also, it is typically not possible to extend 
        ! the tube further than this. 

        ! Declare variables
        !==================
        ! Arguments
        type(GGTMFieldLinePairDataUDT), intent(inout)   :: tube 
        integer(I8), intent(in)                         :: segID, cellID
        integer(I8), intent(inout)                      :: vertID
        type(GGTMDataUDT), intent(inout)                :: ggtmdata 
        logical, intent(in)                             :: start 
        type(MagneticFieldUDT), intent(in)              :: magneticField

        ! Auxiliary
        integer(I8)                             :: newsegID, ind, indhf, &
            indlf, newhfvert, newlfvert
        integer(I8), allocatable, dimension(:)  :: newvertID, segnodes
        real(R8)                                :: hftx, hfty
        real(R8), allocatable, dimension(:)     :: segtx, segty, &
            By, Bx, Bhfx, Bhfy, Byf, Bxf, segdp
        logical                                 :: doflip

        ! Loop
        integer(I8)                             :: k 

        ! Checks
        !=======
        ! Sanity checks
        if (start .and. tube%isextendedstart) then 
            print *, 'ExtendTubeWithSegment: start has already been ' // &
                'extended, returning...'
            return 
        end if 
        if (.not. start .and. tube%isextendedend) then 
            print *, 'ExtendTubeWithSegment: end has already been ' // &
                'extended, returning...'
            return 
        end if 
        if (segID <= 0) then 
            call gdErrorHandler('ExtendTubeWithSegment: segID should be ' // & 
                'greater than 0')
        end if
        
        ! Compute data
        !=============
        ! Get segment
        associate(&
            interp      => magneticField%interp,    &
            srfline     => ggtmdata%face(ggtmdata%cell(cellID)%srf)%line,    &
            erfline     => ggtmdata%face(ggtmdata%cell(cellID)%erf)%line,    & 
            tseg        => ggtmdata%seg(segID),     &
            hfline      => tube%hfline,         &
            lfline      => tube%lfline          &
            )

        ! Evaluate magnetic field at segment nodes
        allocate(By(tseg%nl), Bx(tseg%nl))
        call interp%Evaluate(tseg%xl, tseg%yl, 0, 1, Bx)
        call interp%Evaluate(tseg%xl, tseg%yl, 1, 0, By)
        Bx = -Bx
        Bxf = 0.5*(Bx(2:tseg%nl) + Bx(1:tseg%nl-1))
        Byf = 0.5*(By(2:tseg%nl) + By(1:tseg%nl-1))

        ! Compute dot products and other data - we orient the segment data
        ! from hfline to lfline
        doflip = .false.
        if (start) then 
            ! Evaluate magnetic field at seg
            ! Get segment tangents in hfline starting point
            if (tseg%sv == hfline%vert(1)) then 
                segtx = tseg%xl(2:tseg%nl) - tseg%xl(1:tseg%nl-1)
                segty = tseg%yl(2:tseg%nl) - tseg%yl(1:tseg%nl-1)
            elseif (tseg%ev == hfline%vert(1)) then 
                doflip = .true.
                segtx = -tseg%xl(tseg%nl:2:-1) + tseg%xl(tseg%nl-1:1:-1)
                segty = -tseg%yl(tseg%nl:2:-1) + tseg%yl(tseg%nl-1:1:-1)
                Bx = Bx(size(Bx):1:-1)
                By = By(size(By):1:-1)
                Bxf = Bxf(size(Bxf):1:-1)
                Byf = Byf(size(Byf):1:-1)
            else
                call gdErrorHandler('ExtendTubeWithSegment: segment ' // & 
                    'does not have starting vertex of hfline, unexpected')
            end if 

            ! Compute line tangents
            hftx = hfline%xl(1) - hfline%xl(2)
            hfty = hfline%yl(1) - hfline%yl(2)

            ! Compute magnetic field
            allocate(Bhfx(1), Bhfy(1))
            call interp%Evaluate(hfline%xl(1:1), hfline%yl(1:1), 0, 1, Bhfx)
            call interp%Evaluate(hfline%xl(1:1), hfline%yl(1:1), 1, 0, Bhfy)
            Bhfx = -Bhfx

            ! Compute dot products
            segdp = Bxf*segtx + Byf*segty 
            if ((Bhfx(1)*hftx + Bhfy(1)*hfty) < 0) then 
                segdp = -segdp 
            end if 
        else
            if (tseg%sv == hfline%vert(hfline%nv)) then 
                segtx = tseg%xl(2:tseg%nl) - tseg%xl(1:tseg%nl-1)
                segty = tseg%yl(2:tseg%nl) - tseg%yl(1:tseg%nl-1)
            elseif (tseg%ev == hfline%vert(hfline%nv)) then 
                doflip = .true.
                segtx = -tseg%xl(tseg%nl:2:-1) + tseg%xl(tseg%nl-1:1:-1)
                segty = -tseg%yl(tseg%nl:2:-1) + tseg%yl(tseg%nl-1:1:-1)
                Bx = Bx(size(Bx):1:-1)
                By = By(size(By):1:-1)
                Bxf = Bxf(size(Bxf):1:-1)
                Byf = Byf(size(Byf):1:-1)
            else
                call gdErrorHandler('ExtendTubeWithSegment: segment ' // & 
                    'does not have ending vertex of hfline, unexpected')
            end if 

            ! Compute line tangents
            hftx = hfline%xl(hfline%nl) - hfline%xl(hfline%nl-1)
            hfty = hfline%yl(hfline%nl) - hfline%yl(hfline%nl-1)

            ! Compute magnetic field
            allocate(Bhfx(1), Bhfy(1))
            call interp%Evaluate([hfline%xl(hfline%nl)], [hfline%yl(hfline%nl)], 0, 1, Bhfx)
            call interp%Evaluate([hfline%xl(hfline%nl)], [hfline%yl(hfline%nl)], 1, 0, Bhfy)
            Bhfx = -Bhfx

            ! Compute dot products
            segdp = Bxf*segtx + Byf*segty 
            if ((Bhfx(1)*hftx + Bhfy(1)*hfty) < 0) then 
                segdp = -segdp 
            end if 
        end if 

        ! Construct segment nodes - for when we need to flip the segment
        segnodes = [(k, k = 1, tseg%nl)]
        if (doflip) then 
            segnodes = segnodes(tseg%nl:1:-1)
        end if 

        ! Extend
        !=======
        if (start) then 
            if (segdp(1) >= 0 .and. segdp(size(segdp)) >= 0) then ! Extend the hfline

                ! No extension of lfline possible since edge goes against mf
                ! (segdp remains positive at end)
                
                ! Check if we can take the segment as a whole
                if (all(segdp > 0)) then 
                    ! Just append the segment
                    call hfline%AppendSegment(segID, ggtmdata, .false.)
                else
                    ! Construct a new segment by finding the first location 
                    ! where the sign switches
                    ind = findloc(segdp <= 0, .true., 1, back=.false.)

                    ! Split the line at this node
                    call srfline%SplitAtNodes(segID, [segnodes(ind)], &
                        ggtmdata, vertID, newvertID)

                    ! Get the segment that has the split vertex and the new
                    ! vertex
                    newsegID = GetSegmentIDFromVertices(srfline, ggtmdata, newvertID(1), &
                        hfline%vert(1))
                    
                    ! Append this segment
                    call hfline%AppendSegment(newsegID, ggtmdata, .false.)
                end if 
            
            elseif (segdp(size(segdp)) <= 0 .and. segdp(1) <= 0) then ! Extend the lfline

                ! The hfline cannot be extended (dot product is negative)
                ! but the lfline can (dot product at end is negative) 

                ! Print warning that this code was not yet verified
                !print *, 'ExtendTubeWithSegment: code part not yet verified'

                ! Check if we can take the segment as a whole
                if (all(segdp <= 0)) then 
                    ! Just append the segment
                    call lfline%AppendSegment(segID, ggtmdata, .false.)
                else
                    ! Construct a new segment by finding the first location 
                    ! where the sign switches
                    ind = findloc(segdp >= 0, .true., 1, back=.true.)

                    ! Split the line at this node
                    call srfline%SplitAtNodes(segID, [segnodes(ind)], &
                        ggtmdata, vertID, newvertID)

                    ! Get the segment that has the split vertex and the new
                    ! vertex
                    newsegID = GetSegmentIDFromVertices(srfline, ggtmdata, newvertID(1), &
                        lfline%vert(1))
                    
                    ! Append this segment
                    call lfline%AppendSegment(newsegID, ggtmdata, .false.)
                end if 
            
            elseif (segdp(1) <= 0 .and. segdp(size(segdp)) >= 0) then ! Recompute

                ! Both segments go against the magnetic field - need to make
                ! straight segment between lines

                ! Print warning that this code was not yet verified
                print *, 'ExtendTubeWithSegment: code part not yet verified'

                ! Adjust the segment and update the line
                call tseg%Initialize(tseg%xl([1, tseg%nl]), &
                    tseg%yl([1, tseg%nl]), tseg%fsID, tseg%TMfaceID, &
                    tseg%sv, tseg%ev)
                call srfline%UpdateLineData(ggtmdata)

                ! Recompute the magnetic field
                Bxf = 0.5*(Bx(1:1) + [Bx(size(Bx))])
                Byf = 0.5*(By(1:1) + [By(size(By))])

                ! Recompute the dot product
                ! Get segment tangents in hfline starting point
                if (tseg%sv == hfline%vert(1)) then 
                    segtx = tseg%xl(2:tseg%nl) - tseg%xl(1:tseg%nl-1)
                    segty = tseg%yl(2:tseg%nl) - tseg%yl(1:tseg%nl-1)
                elseif (tseg%ev == hfline%vert(1)) then 
                    doflip = .true.
                    segtx = -tseg%xl(tseg%nl:2:-1) + tseg%xl(tseg%nl-1:1:-1)
                    segty = -tseg%yl(tseg%nl:2:-1) + tseg%yl(tseg%nl-1:1:-1)
                else
                    call gdErrorHandler('ExtendTubeWithSegment: segment ' // & 
                        'does not have starting vertex of hfline, unexpected')
                end if 

                ! Compute dot products
                segdp = Bxf*segtx + Byf*segty 
                if ((Bhfx(1)*hftx + Bhfy(1)*hfty) < 0) then 
                    segdp = -segdp
                end if 

                ! Extend
                if (segdp(1) >= 0) then 
                    call hfline%AppendSegment(segID, ggtmdata, .false.)
                else
                    call lfline%AppendSegment(segID, ggtmdata, .false.)
                end if

            else ! Need to split the line at different points

                ! Print warning that this code was not yet verified
                print *, 'ExtendTubeWithSegment: code part not yet verified'

                ! Construct a new segment by finding the first location 
                ! where the sign switches
                indhf = findloc(segdp <= 0, .true., 1, back=.false.)
                indlf = findloc(segdp >= 0, .true., 1, back=.true.)

                ! Update from edge index to vertex index
                if (indhf /= 0) then  
                    !if (indhf == size(segdp)) then  ! since we trace from the start
                        !indhf = indhf + 1
                    !end if 
                end if 
                if (indlf /= 0) then
                    !if (indlf /= 1) then  ! since we trace from the back,
                        indlf = indlf + 1
                    !end if
                end if 

                ! Check how many segments we get
                if (indhf == 0 .and. indlf == 0) then 
                    ! This shouldn't be happening here, case already 
                    ! hedged for upstream
                    call gdErrorHandler('ExtendTubeWithSegment: could not ' // &
                        'find any split vertices, unexpected')
                elseif (indhf == 0) then 
                    ! Only extension of lfline, only one segment
                    call srfline%SplitAtNodes(segID, [segnodes(indlf)], &
                        ggtmdata, vertID, newvertID)
                    newhfvert = newvertID(1)
                    newlfvert = newvertID(1)
                elseif (indlf == 0) then 
                    ! Only extension of hfline, only one segment
                    call srfline%SplitAtNodes(segID, [segnodes(indhf)], &
                        ggtmdata, vertID, newvertID)
                    newhfvert = newvertID(1)
                    newlfvert = newvertID(1)
                elseif (indhf == indlf) then 
                    ! Just one segment
                    call srfline%SplitAtNodes(segID, [segnodes(indhf)], &
                        ggtmdata, vertID, newvertID)
                    newhfvert = newvertID(1)
                    newlfvert = newvertID(1)
                else
                    ! Two segments
                    call srfline%SplitAtNodes(segID, segnodes([indhf, indlf]), &
                        ggtmdata, vertID, newvertID)
                    newhfvert = newvertID(1)
                    newlfvert = newvertID(2)
                end if

                ! Check which lines to extend
                if (indhf /= 0) then 
                    ! Get the segment that has the split vertex and the new
                    ! vertex for the hfline
                    newsegID = GetSegmentIDFromVertices(srfline, ggtmdata, newhfvert, &
                        hfline%vert(1))
                    
                    ! Append this segment
                    call hfline%AppendSegment(newsegID, ggtmdata, .false.)
                end if 

                if (indlf /= 0) then 
                    ! Get the segment that has the split vertex and the new
                    ! vertex for the lfline
                    newsegID = GetSegmentIDFromVertices(srfline, ggtmdata, newlfvert, &
                        lfline%vert(1))
                    
                    ! Append this segment
                    call lfline%AppendSegment(newsegID, ggtmdata, .false.)
                end if 
            end if   
        else
            if (segdp(1) >= 0 .and. segdp(size(segdp)) >= 0) then ! Extend the hfline

                ! Print warning that this code was not yet verified
                !print *, 'ExtendTubeWithSegment: code part not yet verified'
                
                ! Check if we can take the segment as a whole
                if (all(segdp > 0)) then 
                    ! Just append the segment
                    call hfline%AppendSegment(segID, ggtmdata, .true.)
                else
                    ! Construct a new segment by finding the first location 
                    ! where the sign switches
                    ind = findloc(segdp <= 0, .true., 1, back=.false.)

                    ! Split the line at this node
                    call erfline%SplitAtNodes(segID, [segnodes(ind)], &
                        ggtmdata, vertID, newvertID)

                    ! Get the segment that has the split vertex and the new
                    ! vertex
                    newsegID = GetSegmentIDFromVertices(erfline, ggtmdata, newvertID(1), &
                        hfline%vert(hfline%nv))
                    
                    ! Append this segment
                    call hfline%AppendSegment(newsegID, ggtmdata, .true.)
                end if 
            
            elseif (segdp(size(segdp)) <= 0 .and. segdp(1) <= 0) then ! Extend the lfline

                ! Print warning that this code was not yet verified
                print *, 'ExtendTubeWithSegment: code part not yet verified'
                
                ! Check if we can take the segment as a whole
                if (all(segdp <= 0)) then 
                    ! Just append the segment
                    call lfline%AppendSegment(segID, ggtmdata, .true.)
                else
                    ! Construct a new segment by finding the first location 
                    ! where the sign switches
                    ind = findloc(segdp >= 0, .true., 1, back=.true.)

                    ! Split the line at this node
                    call erfline%SplitAtNodes(segID, [segnodes(ind)], &
                        ggtmdata, vertID, newvertID)

                    ! Get the segment that has the split vertex and the new
                    ! vertex
                    newsegID = GetSegmentIDFromVertices(erfline, ggtmdata, newvertID(1), &
                        lfline%vert(lfline%nv))
                    
                    ! Append this segment
                    call lfline%AppendSegment(newsegID, ggtmdata, .true.)
                end if 
            
            elseif (segdp(1) <= 0 .and. segdp(size(segdp)) >= 0) then ! Recompute

                ! Print warning that this code was not yet verified
                print *, 'ExtendTubeWithSegment: code part not yet verified'

                ! Adjust the segment and update the line
                call tseg%Initialize(tseg%xl([1, tseg%nl]), &
                    tseg%yl([1, tseg%nl]), tseg%fsID, tseg%TMfaceID, &
                    tseg%sv, tseg%ev)
                call erfline%UpdateLineData(ggtmdata)

                ! Recompute the dot product
                ! Get segment tangents in hfline starting point
                if (tseg%sv == hfline%vert(hfline%nv)) then 
                    segtx = tseg%xl(2:tseg%nl) - tseg%xl(1:tseg%nl-1)
                    segty = tseg%yl(2:tseg%nl) - tseg%yl(1:tseg%nl-1)
                elseif (tseg%ev == hfline%vert(hfline%nv)) then 
                    doflip = .true.
                    segtx = -tseg%xl(tseg%nl:2:-1) + tseg%xl(tseg%nl-1:1:-1)
                    segty = -tseg%yl(tseg%nl:2:-1) + tseg%yl(tseg%nl-1:1:-1)
                else
                    call gdErrorHandler('ExtendTubeWithSegment: segment ' // & 
                        'does not have ending vertex of hfline, unexpected')
                end if 

                ! Recompute the magnetic field
                Bxf = 0.5*(Bx(1:1) + [Bx(size(Bx))])
                Byf = 0.5*(By(1:1) + [By(size(By))])

                ! Compute dot products
                segdp = Bxf*segtx + Byf*segty 
                if ((Bhfx(1)*hftx + Bhfy(1)*hfty) < 0) then 
                    segdp = -segdp
                end if 
                    
                ! Extend
                if (segdp(1) >= 0) then 
                    call hfline%AppendSegment(segID, ggtmdata, .true.)
                else
                    call lfline%AppendSegment(segID, ggtmdata, .true.)
                end if

            else ! Need to split the line at different points

                ! Construct a new segment by finding the first location 
                ! where the sign switches
                indhf = findloc(segdp <= 0, .true., 1, back=.false.)
                indlf = findloc(segdp >= 0, .true., 1, back=.true.)

                ! Update from edge index to vertex index
                if (indhf /= 0) then  
                    !if (indhf == size(segdp)) then  ! since we trace from the start
                        !indhf = indhf + 1
                    !end if 
                end if 
                if (indlf /= 0) then
                    !if (indlf /= 1) then  ! since we trace from the back,
                        indlf = indlf + 1
                    !end if
                end if 

                ! Check how many segments we get
                if (indhf == 0 .and. indlf == 0) then 
                    ! This shouldn't be happening here, case already 
                    ! hedged for upstream
                    call gdErrorHandler('ExtendTubeWithSegment: could not ' // &
                        'find any split vertices, unexpected')
                elseif (indhf == 0) then 
                    ! Only extension of lfline, only one segment
                    call erfline%SplitAtNodes(segID, [segnodes(indlf)], &
                        ggtmdata, vertID, newvertID)
                    newhfvert = newvertID(1)
                    newlfvert = newvertID(1)
                elseif (indlf == 0) then 
                    ! Only extension of hfline, only one segment
                    call erfline%SplitAtNodes(segID, [segnodes(indhf)], &
                        ggtmdata, vertID, newvertID)
                    newhfvert = newvertID(1)
                    newlfvert = newvertID(1)
                elseif (indhf == indlf) then 
                    ! Just one segment
                    call erfline%SplitAtNodes(segID, [segnodes(indhf)], &
                        ggtmdata, vertID, newvertID)
                    newhfvert = newvertID(1)
                    newlfvert = newvertID(1)
                else
                    ! Two segments
                    call erfline%SplitAtNodes(segID, segnodes([indhf, indlf]), &
                        ggtmdata, vertID, newvertID)
                    newhfvert = newvertID(1)
                    newlfvert = newvertID(2)
                end if

                ! Check which lines to extend
                if (indhf /= 0) then 
                    ! Get the segment that has the split vertex and the new
                    ! vertex for the hfline
                    newsegID = GetSegmentIDFromVertices(erfline, ggtmdata, newhfvert, &
                        hfline%vert(hfline%nv))
                    
                    ! Append this segment
                    call hfline%AppendSegment(newsegID, ggtmdata, .true.)
                end if 

                if (indlf /= 0) then 
                    ! Get the segment that has the split vertex and the new
                    ! vertex for the lfline
                    newsegID = GetSegmentIDFromVertices(erfline, ggtmdata, newlfvert, &
                        lfline%vert(lfline%nv))
                    
                    ! Append this segment
                    call lfline%AppendSegment(newsegID, ggtmdata, .true.)
                end if 

            end if   
        end if 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    !------------------------------------------------------------------!
    !                      GGTM LINE REFINERS                          !
    !------------------------------------------------------------------!

    ! Refiner initialization
    function InitializeGGTMLineRefiner(topomesh, &
        magneticField, vessel, fieldtracer, boundarytracer, &
        poloidalvertexdistributor, radialvertexdistributor, options, &
        direction) result(GGTMlinerefiner)

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
        character(*), intent(in)                :: direction

        ! Auxiliary
        integer(I8), allocatable, dimension(:)      :: tv
        integer(I8), allocatable, dimension(:, :)   :: labels
        logical, allocatable, dimension(:)          :: includevesselvertex
        real(R8), allocatable, dimension(:)         :: xp, yp, &
            valpLmin, valpLmax, decaylength, xv, yv, tempLmin, tempLmax, &
            tempdecaylength
        type(Coordinates2DDistanceDFUDT)            :: Lmin, Lmax

        ! Loop
        integer(I8)                                 :: i 

        ! Select refiner
        !===============
        select case (direction)

        case ('poloidal')

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

        case ('radial')

            select case (options%radrefmeth)

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

        case default

            call gdErrorHandler('InitializeGGTMLineRefiner: direction ' // & 
                ' "' // direction // '" not implemented')

        end select

        ! Initialize refiner
        !===================
        select type (GGTMlinerefiner)

        type is (GGTMLineRefinerNoRefUDT)

            ! Do nothing

        type is (GGTMLineRefinerLB2DUDT)

            ! Check how to generate min and max distributions - note: 
            ! currently, we only do point-based distribution methods

            select case (direction)

            case ('poloidal')
        
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

                ! Visualize
                call Lmin%Visualize([minval(topomesh%vert%x), maxval(topomesh%vert%x)], &
                    [minval(topomesh%vert%y), maxval(topomesh%vert%y)], &
                    100, 100, 'Lminpolref')
                call Lmax%Visualize([minval(topomesh%vert%x), maxval(topomesh%vert%x)], &
                    [minval(topomesh%vert%y), maxval(topomesh%vert%y)], &
                    100, 100, 'Lmaxpolref')

                ! Construct refiner
                GGTMlinerefiner = ConstructGGTMLineRefiner(Lmin, Lmax, 'classic', &
                options%reflengthtype, magneticField)

            case ('radial') 

                ! Check which points to include
                allocate(xp(0), yp(0), valpLmin(0), valpLmax(0), decaylength(0))

                ! Include x-point/strike point regions?
                if (options%radrefLBdosp) then 
                    ! Add all strike and x-points
                    tv = [topomesh%GetStrikePointIDs(), topomesh%GetXPointIDs()]
                    xp = [xp, topomesh%vert%x(tv)]
                    yp = [yp, topomesh%vert%y(tv)]
                    valpLmin = [valpLmin, spread(options%radrefLBLminsp, 1, size(tv))]
                    valpLmax = [valpLmax, spread(options%radrefLBLmaxsp, 1, size(tv))]
                    decaylength = [decaylength, spread(options%radrefLBdecaylengthsp, 1, size(tv))]

                end if

                ! Construct distribution functions
                call Lmin%Initialize(xp, yp, valpLmin, options%radrefLBLmininf, &
                    decaylength)
                call Lmax%Initialize(xp, yp, valpLmax, options%radrefLBLmaxinf, &
                    decaylength)

                ! Visualize
                call Lmin%Visualize([minval(topomesh%vert%x), maxval(topomesh%vert%x)], &
                    [minval(topomesh%vert%y), maxval(topomesh%vert%y)], &
                    100, 100, 'Lminradref')
                call Lmax%Visualize([minval(topomesh%vert%x), maxval(topomesh%vert%x)], &
                    [minval(topomesh%vert%y), maxval(topomesh%vert%y)], &
                    100, 100, 'Lmaxradref')

                ! Construct refiner
                GGTMlinerefiner = ConstructGGTMLineRefiner(Lmin, Lmax, 'classic', &
                options%radreflengthtype, magneticField)

            case default

                call gdErrorHandler('InitializeGGTMLineRefiner: direction ' // & 
                    ' "' // direction // '" not implemented')

            end select


        class default 

            call gdErrorHandler('INitializeGGTMLineRefiner: type not implemented')

        end select


    end function 

    ! Refiner option updating, dummy
    subroutine UpdateRefinementOptionsNoRef(refiner, refoptions, topomesh)

        ! Nothing to do here, move along
        class(GGTMLineRefinerNoRefUDT)                      :: refiner 
        type(GGTMFieldlineRefinementOptionsUDT), intent(in) :: refoptions 
        type(TopomeshUDT), intent(in)                       :: topomesh

    end subroutine

    ! Projection, dummy
    subroutine ProjectLineVertexDistributionNoRef(refiner, linein, &
        lineout, vertID, ggtmdata)

        ! Description
        !============
        ! Project the current length distribution in the incoming line
        ! onto the length distribution of the outgoing line. For this 
        ! refiner, this simply means interpolating the length 
        ! distribution of the vertices of the incoming line (except for
        ! the nodes to be kept). The original vertex distribution there 
        ! is lost of course. Small faces (smaller than distttol) are 
        ! removed. 

        ! Note: it is assumed that dlcv goes from 0 to the length of the
        ! line. Further, we assume that the vertex distribution of the
        ! incoming line has no very small faces (i.e. they are removed
        ! before this routine) and therefore small faces can only 
        ! originate from/near vertices to be kept in the outgoing line. 

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMLineRefinerNoRefUDT)              :: refiner 
        type(GGTMFieldlineDataUDT), intent(in)      :: linein
        type(GGTMFieldlineDataUDT), intent(inout)   :: lineout
        integer(I8), intent(inout)                  :: vertID 
        type(GGTMDataUDT), intent(inout)            :: ggtmdata

        ! Auxiliary
        real(R8), allocatable, dimension(:)         :: tdlcv, dlv, &
            refinerdllcin, refinerdllcout, refinerdlcvin, refinerdlcvout, &
            refinerdlcvoutinit
        integer(I8), allocatable, dimension(:)      :: tvertID, &
            sortind, segvertind
        logical, allocatable, dimension(:)          :: keepind, &
            tisnodevert

        ! Loop
        integer(I8)                                 :: i, k 

        ! Unpack
        !=======
        associate(&
            dllcin      => linein%dllc,     &
            dlcvin      => linein%dlcv,     &
            dlcvout     => lineout%dlcv,    &
            dllcout     => lineout%dllc     &
            )

        ! Project
        !========
        ! Determine line lengths in correct length measure
        refinerdllcin = linein%dllc
        refinerdllcout = lineout%dllc

        ! Initialize vertices etc
        tvertID = [(k, k = vertID+1, vertID+size(linein%dlcv))]
        tisnodevert = [spread(.false., 1, size(linein%dlcv))]

        ! Update vertex index
        vertID = vertID + size(linein%dlcv)

        ! Interpolate the incoming line length distribution
        call Interpolate1D(linein%dlcv, refinerdlcvin, linein%dllc, refinerdllcin)

        ! Project onto the outgoing line length distribution
        refinerdlcvout = refinerdlcvin/refinerdllcin(size(refinerdllcin))*refinerdllcout(size(refinerdllcout))

        ! Append vertices to be kept
        refinerdlcvoutinit = lineout%dlcv
        segvertind = lineout%GetAllSegmentVertIndices()
        refinerdlcvout = [refinerdlcvout, refinerdlcvoutinit(segvertind)]
        tvertID = [tvertID, lineout%vert(segvertind)]
        tisnodevert = [tisnodevert, lineout%isnodevert(segvertind)]
        
        ! Interpolate the outgoing actual line length
        call Interpolate1D(refinerdlcvout, tdlcv, refinerdllcout, lineout%dllc)
        !tdlcv = refinerdlcvout/refinerdlcvout(size(refinerdlcvout))*&
        !    refinerdllcout(size(refinerdllcout))

        ! Insert vertices to be kept
        allocate(sortind(size(tdlcv)))
        call Sort(tdlcv, ind=sortind, ascend=.true.)
        tvertID     = tvertID(sortind)
        tisnodevert = tisnodevert(sortind)
        deallocate(sortind)

        ! Check for very small distances
        allocate(keepind(size(tdlcv)))
        keepind = .true.
        dlv = abs(tdlcv(2:) - tdlcv(1:size(tdlcv)-1))
        do i = 1, size(dlv)
            if (dlv(i) < disttol) then 
                if (tisnodevert(i) .and. .not. tisnodevert(i+1)) then
                    keepind(i+1) = .false.
                elseif (.not. tisnodevert(i) .and. tisnodevert(i+1)) then 
                    keepind(i) = .false.
                elseif (tvertID(i) == tvertID(i+1)) then 
                    ! Same vertices, can delete one 
                    keepind(i) = .false.
                else
                    ! Can't delete nodes, small faces will be present...
                end if  
            end if 
        end do 

        ! Delete vertices to remove small distances
        tdlcv = pack(tdlcv, keepind)
        tvertID = pack(tvertID, keepind)
        tisnodevert = pack(tisnodevert, keepind)
        tdlcv(1) = 0 ! ensure start and end on line
        tdlcv(size(tdlcv)) = lineout%dllc(lineout%nl)
        tisnodevert(1) = .true.
        tisnodevert(size(tdlcv)) = .true.

        ! Reconstruct line
        !=================
        ! Add coordinates
        call lineout%AddVertexCoordinates(tdlcv)

        ! Add line vertices
        call lineout%AddVertexIDs(tvertID, tisnodevert)

        ! Update segment data
        call lineout%UpdateSegmentData(ggtmdata)

        ! Housekeeping
        !=============
        end associate
        
    end subroutine

    ! Refiner option updating, length-based
    subroutine UpdateRefinementOptionsLB(refiner, refoptions, topomesh)

        ! Description
        !============
        ! Update the length distribution options for the length based 
        ! refiner. 

        ! Note: currently we don't yet update the length distribution
        ! itself, though this may be done for each different cell. 

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMLineRefinerLB2DUDT)                       :: refiner 
        type(GGTMFieldlineRefinementOptionsUDT), intent(in) :: refoptions 
        type(TopomeshUDT), intent(in)                       :: topomesh

        ! Update boundary layer
        !======================
        refiner%dlBLlengthbased     = refoptions%dlBLlengthbased
        refiner%doBLstart   = refoptions%doBLstart 
        refiner%doBLend     = refoptions%doBLend 
        refiner%ncBLstart   = refoptions%ncBLstart 
        refiner%ncBLend     = refoptions%ncBLend 
        refiner%dlBLstart   = refoptions%dlBLstart 
        refiner%dlBLend     = refoptions%dlBLend 

    end subroutine

    ! Length-based refiner constructor
    function ConstructGGTMLineRefinerLB(Lmin, Lmax, meth, lengthtype, &
        magneticField) result(refiner)

        ! Description
        !============
        ! Constructor for line refinement, based on minimal and maximal
        ! length distributions.

        ! Declare variables
        !==================
        ! Arguments
        type(GGTMLineRefinerLB2DUDT)                :: refiner 
        class(DistributionFunctionUDT), intent(in)  :: Lmin, Lmax
        character(*), intent(in)                    :: meth, lengthtype 
        type(MagneticFieldUDT), intent(in)          :: magneticField
 
        ! Initialize
        !===========
        ! General refinement data
        refiner%Lmin        = Lmin 
        refiner%Lmax        = Lmax 
        refiner%meth        = meth
        refiner%lengthtype  = lengthtype 
        refiner%field       = magneticField

        ! Boundary layer (simply set to zero currently)
        refiner%dlBLlengthbased = .false.
        refiner%doBLstart = .false. 
        refiner%doBLend = .false. 
        if (allocated(refiner%dlBLstart)) then 
            deallocate(refiner%dlBLstart)
        end if 
        if (allocated(refiner%dlBLend)) then 
            deallocate(refiner%dlBLend)
        end if 
        allocate(refiner%dlBLstart(0), refiner%dlBLend(0))
        refiner%ncBLstart = 0
        refiner%ncBLend = 0

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

        ! Note: we now also support imposing a boundary layer of cells.
        ! If the boundary layer would exceed the length of the line, then
        ! cells are deleted until the length is matched (it is possible
        ! no boundary layer is then imposed). Also, if nodes should be
        ! kept that are in the range of the boundary layer, then 

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMLineRefinerLB2DUDT)               :: refiner 
        type(GGTMFieldlineDataUDT), intent(inout)   :: line 
        integer(I8), intent(inout)                  :: vertID 
        logical, intent(in)                         :: keepvert(:)

        ! Auxiliary
        logical                                     :: ismerged, &
            isnextfacelegal, isprevfacelegal, dostart
        logical, allocatable, dimension(:)          :: isreflegal, &
            iscoarselegal, thiskeepvert, newkeepvert, refineface, &
            coarsenface,  newisreflegal, newiscoarselegal, iscoarsenedface, &
            keepind, keepvertBLstart, iscoarselegalBLstart, &
            isreflegalBLstart, keepvertBLend, iscoarselegalBLend, &
            isreflegalBLend
        real(R8)                                    :: lline, &
            lBLstart, lBLend
        real(R8), allocatable, dimension(:)         :: dll, &
            Lmaxvert, Lminvert, newdll, newdllc, dlcBLstart, &
            dlcBLend, tempdlcBLstart, tempdlcBLend, dlBLstart, &
            dlBLend, tdlcBLstart, tdlcBLend, dlcBLstarti, dlcBLendi
        integer(I8)                                 :: minind, &
            ncBLstart, ncBLend
        integer(I8), allocatable, dimension(:)      :: thisvertID, &
            newvertID, thisvertIDBLstart, sortind, thisvertIDBLend

        ! Loop
        integer(I8)                                 :: i, k, cc 

        ! Initialize
        !===========
        ! Associate
        associate(&
            lengthbased     => refiner%dlBLlengthbased,    &
            doBLstart       => refiner%doBLstart,   &
            doBLend         => refiner%doBLend,     &
            dlBLstarti      => refiner%dlBLstart,   &
            dlBLendi        => refiner%dlBLend      &
            )

        ! Initialize the refiner
        call refiner%InitializeLineData(line)

        ! Check if distance is (euler) length based or not (if not, 
        ! assumed given in current units).
        if (lengthbased) then 
            ! Need to interpolate
            
            if (doBLstart) then 
                allocate(dlcBLstarti(size(dlBLstarti)+1))
                dlcBLstarti = 0
                do i = 2, size(dlcBLstarti)
                    dlcBLstarti(i) = dlcBLstarti(i-1) + dlBLstarti(i-1)
                end do
                call Interpolate1D(dlcBLstarti, tdlcBLstart, line%dllc, refiner%linedllc)
                dlBLstart = tdlcBLstart(2:size(tdlcBLstart)) - tdlcBLstart(1:size(tdlcBLstart)-1)
                deallocate(dlcBLstarti)
            end if 
            if (doBLend) then 
                allocate(dlcBLendi(size(dlBLendi)+1))
                dlcBLendi = 0
                do i = 2, size(dlcBLendi)
                    dlcBLendi(i) = dlcBLendi(i-1) + dlBLendi(i-1)
                end do
                dlcBLendi = line%dllc(line%nl) - dlcBLendi
                call Interpolate1D(dlcBLendi, tdlcBLend, line%dllc, refiner%linedllc)
                dlBLend = -tdlcBLend(2:size(tdlcBLend)) + tdlcBLend(1:size(tdlcBLend)-1)
                deallocate(dlcBLendi)
            end if 
        else
            ! Just take as is
            dlBLstart = dlBLstarti
            dlBLend = dlBLendi
        end if

        ! Ensure proper dimensions
        if (size(keepvert) /= size(line%vert)) then 
            call gdErrorHandler('RefineLineSingleLB: keepvert does not '// & 
                'have same number of elements as line%vert, check input')
        end if 

        ! Initialize
        ncBLstart = refiner%ncBLstart
        ncBLend = refiner%ncBLend

        ! Set initial logicals
        allocate(iscoarselegal(line%nv-1))
        iscoarselegal = .true.
        isreflegal = iscoarselegal
        thiskeepvert = keepvert
        thisvertID = line%vert

        ! Initial distribution
        newdllc = refiner%GetLineEdgeLength(line)
        dll = newdllc(2:size(newdllc)) - newdllc(1:size(newdllc)-1)
        lline = newdllc(size(newdllc))

        ! Initialize boundary layer distributions
        if (doBLstart) then 
            ! Check if the length allows it
            allocate(dlcBLstart(ncBLstart+1))
            dlcBLstart = 0
            do i = 2, size(dlcBLstart)
                dlcBLstart(i) = dlcBLstart(i-1) + dlBLstart(i-1)
            end do 
            lBLstart = dlcBLstart(ncBLstart+1)
        else
            allocate(dlcBLstart(0))
            lBLstart = 0
        end if 
        if (doBLend) then 
            ! Check if the length allows it
            allocate(dlcBLend(ncBLend+1))
            dlcBLend = lline
            do i = size(dlcBLend)-1, 1, -1 
                dlcBLend(i) = dlcBLend(i+1) - dlBLend(i)
            end do 
            lBLend = lline - dlcBLend(1)
        else
            allocate(dlcBLend(0))
            lBLend = 0
        end if 

        ! Adjust boundary layer distribution if necessary
        dostart = .true. 
        do while ((lBLstart + lBLend > lline) .or. isnan(lBLstart) .or. &
            isnan(lBLend)) ! account for out of bounds interpolation
            if (dostart .and. (lBLstart > 0.0_R8 .or. isnan(lBLstart))) then
                dlcBLstart = dlcBLstart(1:size(dlcBLstart)-1)
                lBLstart = dlcBLstart(size(dlcBLstart))
                dostart = .false.
            elseif (.not. dostart .and. (lBLend > 0.0_R8 .or. isnan(lBLend))) then 
                dlcBLend = dlcBLend(2:size(dlcBLend))
                lBLend = lline - dlcBLend(1)
                dostart = .true. 
            elseif (lBLstart == 0.0_R8 .and. lBLend == 0.0_R8) then 
                exit ! we should in fact not reach this if lline >= 0
            elseif (dostart) then 
                dostart = .false. 
            elseif (.not. dostart) then 
                dostart = .true. 
            end if 
        end do 

        ! Adjust initial distribution, which vertices to keep etc
        ncBLstart = size(dlcBLstart)-1
        ncBLend = size(dlcBLend)-1
        if (size(dlcBLstart) > 0) then 
            ! Delete nodes that are inside boundary layer (unless they 
            ! are kept, then they replace)
            keepind = newdllc > lBLstart

            ! Check if any vertices are kept in the boundary layer part.
            ! If yes, check how to insert them (they replace a BL vert)
            k = 1
            tempdlcBLstart = dlcBLstart
            allocate(keepvertBLstart(ncBLstart+1), iscoarselegalBLstart(ncBLstart), &
                isreflegalBLstart(ncBLstart), thisvertIDBLstart(ncBLstart+1))
            keepvertBLstart         = .true.
            thisvertIDBLstart       = [(k, k = vertID+1, vertID+ncBLstart+2)]
            vertID = vertID + ncBLstart + 2
            iscoarselegalBLstart    = .false.
            isreflegalBLstart       = .false.
            do while (.not. keepind(k))
                if (thiskeepvert(k)) then 
                    ! Need to keep this vertex
                    if (all(tempdlcBLstart == Posinfval_R8())) then 
                        ! Append
                        print *, 'RefineSingleLB: code not verified'
                        dlcBLstart              = [dlcBLstart, newdllc(k)]
                        keepvertBLstart         = [keepvertBLstart, .true.]
                        thisvertIDBLstart       = [thisvertIDBLstart, thisvertID(k)]
                        iscoarselegalBLstart    = [iscoarselegalBLstart, .false.]
                        isreflegalBLstart       = [isreflegalBLstart, .false.]    

                    else
                        ! Check which vertex to replace
                        minind = minloc(abs(tempdlcBLstart - newdllc(k)), 1)
                        tempdlcBLstart(minind) = Posinfval_R8()
                        dlcBLstart(minind) = newdllc(k)
                        thisvertIDBLstart(minind) = thisvertID(k)
                    end if 
                    
                else
                    ! Overwrite this vertex
                    ! (nothing to do)

                end if 
                
                ! Increment counter
                k = k + 1
            end do

            ! Sort
            allocate(sortind(size(dlcBLstart)))
            sortind = [(k, k = 1, size(dlcBLstart))]
            call Sort(dlcBLstart, ind=sortind, ascend=.true.)
            thisvertIDBLstart = thisvertIDBLstart(sortind)
            deallocate(sortind)

            ! Rebuild 
            newdllc = [dlcBLstart, pack(newdllc, keepind)]
            thiskeepvert = [keepvertBLstart, pack(thiskeepvert, keepind)]
            thisvertID = [thisvertIDBLstart, pack(thisvertID, keepind)]
            iscoarselegal = [iscoarselegalBLstart, pack(iscoarselegal, keepind(2:))]
            isreflegal = [isreflegalBLStart, pack(isreflegal, keepind(2:))]
        end if 
        if (size(dlcBLend) > 0) then 
            ! Delete nodes that are inside boundary layer (unless they 
            ! are kept, then they replace)
            keepind = newdllc < (lline - lBLend)

            ! Check if any vertices are kept in the boundary layer part.
            ! If yes, check how to insert them (they replace a BL vert)
            k = size(newdllc)
            tempdlcBLend = dlcBLend
            allocate(keepvertBLend(ncBLend+1), iscoarselegalBLend(ncBLend), &
                isreflegalBLend(ncBLend), thisvertIDBLend(ncBLend+1))
            keepvertBLend         = .true.
            thisvertIDBLend       = [(k, k = vertID+1, vertID+ncBLend+2)]
            vertID = vertID + ncBLend + 2
            iscoarselegalBLend    = .false.
            isreflegalBLend       = .false.
            do while (.not. keepind(k))
                if (thiskeepvert(k)) then 
                    ! Need to keep this vertex
                    if (all(tempdlcBLend == Posinfval_R8())) then 
                        ! Append (will sort later anyway)
                        print *, 'RefineSingleLB: code not verified'
                        dlcBLend              = [dlcBLend, newdllc(k)]
                        keepvertBLend         = [keepvertBLend, .true.]
                        thisvertIDBLend       = [thisvertIDBLend, thisvertID(k)]
                        iscoarselegalBLend    = [iscoarselegalBLend, .false.]
                        isreflegalBLend       = [isreflegalBLend, .false.]    

                    else
                        ! Check which vertex to replace
                        minind = minloc(abs(tempdlcBLend - newdllc(k)), 1)
                        tempdlcBLend(minind) = Posinfval_R8()
                        dlcBLend(minind) = newdllc(k)
                        thisvertIDBLend(minind) = thisvertID(k)
                    end if 
                    
                else
                    ! Overwrite this vertex
                    ! (nothing to do)

                end if 
                
                ! Decrement counter
                k = k - 1
            end do

            ! Sort
            allocate(sortind(size(dlcBLend)))
            sortind = [(k, k = 1, size(dlcBLend))]
            call Sort(dlcBLend, ind=sortind, ascend=.true.)
            thisvertIDBLend = thisvertIDBLend(sortind)
            deallocate(sortind)

            ! Rebuild 
            newdllc         = [pack(newdllc, keepind), dlcBLend]
            thiskeepvert    = [pack(thiskeepvert, keepind), keepvertBLend]
            thisvertID      = [pack(thisvertID, keepind), thisvertIDBLend]
            iscoarselegal   = [pack(iscoarselegal, keepind(1:size(keepind)-1)), iscoarselegalBLend]
            isreflegal      = [pack(isreflegal, keepind(1:size(keepind)-1)), isreflegalBLend]
        end if

        ! Rebuild
        dll = newdllc(2:size(newdllc)) - newdllc(1:size(newdllc)-1)
        call refiner%AddLineVertexCoordinates(line, newdllc)
        ! call line%AddVertexCoordinates(newdllc)
        if (any(thisvertID == 0)) then 
            print *, thisvertID
        end if 

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

            if (size(isreflegal) /= size(dll)) then 
                print *, 'size error'
            end if

            ! Determine which faces to refine/coarsen
            where (((dll > Lmaxvert(1:line%nv-1)) .or. (dll > Lmaxvert(2:line%nv))) &
                .and. isreflegal) 
                refineface = .true. 
                iscoarselegal = .false. 
            end where
            where (((dll < Lminvert(1:line%nv-1)) .or. (dll < Lminvert(2:line%nv))) &
                .and. iscoarselegal &
                .and. .not. (thiskeepvert(1:line%nv-1) .and. thiskeepvert(2:line%nv)) &
                .and. .not. ([.not. iscoarselegal(1:line%nv-2), .false.] .and. [.false., iscoarselegal(2:line%nv-1)]))
                coarsenface = .true.
                isreflegal = .false.
            end where
                
            ! Check exit conditions
            if ((.not. any(refineface)) .and. (.not. any(coarsenface))) then
                if (any(dll < disttol)) then 
                    print *, minloc(dll)
                end if  
                exit
            end if 
            if (any(dll < disttol.and. .not. coarsenface)) then 
                print *, minloc(dll)
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
                                .and. (.not. newkeepvert(cc))
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
                if (any(newdll(1:cc-1) < disttol)) then 
                    print *, minloc(newdll)
                end if 
                dll = newdll(1:cc-1)
                if (any(dll < disttol)) then 
                    print *, minloc(dll)
                end if 
                deallocate(newdllc)
                allocate(newdllc(size(dll)+1))
                newdllc = 0_R8
                do i = 2, size(newdllc)
                    newdllc(i) = dll(i-1) + newdllc(i-1)
                end do 
                call refiner%AddLineVertexCoordinates(line, newdllc)
                ! call line%AddVertexCoordinates(newdllc)

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
        if (any(dll <= disttol)) then 
            print *, 'weird'
        end if 
        call refiner%AddLineVertexCoordinates(line, newdllc)
        ! call line%AddVertexCoordinates(newdllc)
        call line%AddVertexIDs(thisvertID, thiskeepvert)

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Line initialization
    subroutine InitializeLineDataLB(refiner, line)

        ! Description
        !============
        ! Initialize any line-based data for the refiner. To be called
        ! only once at the beginning of each refinement loop. In this 
        ! case, the line length coordinate in terms of the specified 
        ! length measure is computed and stored.

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMLineRefinerLB2DUDT)           :: refiner
        type(GGTMFieldLineDataUDT), intent(in)  :: line 

        ! Auxiliary
        real(R8), allocatable, dimension(:)     :: dx, dy, xf, yf, bx, &
            by, bn, dll, dllc, psi, dpsi

        ! Loop 
        integer(I8)                             :: i 

        ! Compute
        !========
        select case (refiner%lengthtype)

        case ('euler')

            ! Easy
            refiner%linedllc = line%dllc 

        case ('radial')

            ! Initialize
            dx = line%xl(2:line%nl) - line%xl(1:line%nl-1)
            dy = line%yl(2:line%nl) - line%yl(1:line%nl-1) 
            xf = 0.5*(line%xl(2:line%nl) + line%xl(1:line%nl-1))
            yf = 0.5*(line%yl(2:line%nl) + line%yl(1:line%nl-1)) 
            allocate(bx(line%nl-1), by(line%nl-1))
            call refiner%field%interp%Evaluate(xf, yf, 1, 0, bx)
            call refiner%field%interp%Evaluate(xf, yf, 0, 1, by)
            bn = sqrt(bx**2 + by**2)
            bx = bx/bn 
            by = by/bn 

            ! Project and take absolute value
            dll = dx*bx + dy*by
            if (sum(dll)/size(dll) < 0) then 
                where(dll > 0) dll = 0
            else
                where(dll < 0) dll = 0
            end if 
            dll = abs(dll) 
            
            ! Compute accumulative length
            allocate(dllc(line%nl))
            dllc = 0_R8
            do i = 2, line%nl
                dllc(i) = dllc(i-1) + dll(i-1)
            end do 
            refiner%linedllc = dllc

        case ('psi')

            ! Evaluate psi values on line
            allocate(psi(line%nl))
            call refiner%field%interp%Evaluate(line%xl, line%yl, 0, 0, psi)

            ! Take differences along the line
            dpsi = (psi(2:line%nl) - psi(1:line%nl-1))

            ! Check general difference between first and last, remove
            ! parts that are  not ascending/descending
            if (psi(1) > psi(size(psi))) then 
                ! Descending
                where (dpsi > 0) dpsi = 0
                where((psi(2:) > psi(1)) .or. (psi(2:) < psi(size(psi)))) dpsi = 0
            else 
                ! Ascending
                where(dpsi < 0) dpsi = 0
                where((psi(2:) < psi(1)) .or. (psi(2:) > psi(size(psi)))) dpsi = 0
            end if 
            dpsi = abs(dpsi)
                

            ! Compute accumulative length
            allocate(dllc(line%nl))
            dllc = 0_R8
            do i = 2, line%nl 
                dllc(i) = dllc(i-1) + dpsi(i-1)
            end do
            refiner%linedllc = dllc

        case default 
            
            call gdErrorHandler('InitializeLineDataLB: length ' // & 
                'method "' // refiner%lengthtype // '" not implemented')

        end select


    end subroutine

    ! Line length getter, length based
    function GetLineEdgeLengthLB(refiner, line) result(dlcv)

        ! Description
        !============
        ! Wrapper to get line edge lengths. The wrapper here is 
        ! necessary to allow different length definitions and not only
        ! the eulerian length given by the line's dll field. 

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMLineRefinerLB2DUDT)               :: refiner 
        type(GGTMFieldlineDataUDT), intent(in)      :: line 
        real(R8), allocatable, dimension(:)         :: dlcv


        ! Determine length
        !=================
        select case (refiner%lengthtype)

        case ('euler')

            ! Simply return line%dlcv
            dlcv = line%dlcv

        case ('radial', 'psi')

            ! Need to interpolate 
            call Interpolate1D(line%dlcv, dlcv, line%dllc, refiner%linedllc)

        case default 
            
            call gdErrorHandler('RefineLineSingleLB: length ' // & 
                'method "' // refiner%lengthtype // '" not implemented')

        end select

    end function

    ! Line coordinate setter
    subroutine AddLineVertexCoordinatesLB(refiner, line, dlcv)

        ! Description
        !============
        ! Add the line vertex coordinates based on the new length 
        ! coordinates given in dlcv which are in terms of the refiner
        ! length type. 

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMLineRefinerLB2DUDT)               :: refiner 
        type(GGTMFieldlineDataUDT), intent(in)      :: line 
        real(R8), dimension(:), intent(in)          :: dlcv

        ! Auxiliary
        real(R8), allocatable, dimension(:)         :: newdlcv


        ! Determine length
        !=================
        select case (refiner%lengthtype)

        case ('euler')

            ! Simply call line method
            call line%AddVertexCoordinates(dlcv)

        case ('radial', 'psi')

            ! Need to interpolate first
            call Interpolate1D(dlcv, newdlcv, refiner%linedllc, line%dllc)
            newdlcv(1) = 0_R8
            newdlcv(size(newdlcv)) = line%dllc(size(line%dllc))
            
            ! Call line method
            call line%AddVertexCoordinates(newdlcv)
            
        case default 
            
            call gdErrorHandler('RefineLineSingleLB: length ' // & 
                'method "' // refiner%lengthtype // '" not implemented')

        end select


    end subroutine
    
    ! Projection, length-based
    subroutine ProjectLineVertexDistributionLB(refiner, linein, &
        lineout, vertID, ggtmdata)

        ! Description
        !============
        ! Project the current length distribution in the incoming line
        ! onto the length distribution of the outgoing line. For this 
        ! refiner, this simply means interpolating the length 
        ! distribution of the vertices of the incoming line (except for
        ! the nodes to be kept). The original vertex distribution there 
        ! is lost of course. Small faces (smaller than distttol) are 
        ! removed. 

        ! Note: it is assumed that dlcv goes from 0 to the length of the
        ! line. Further, we assume that the vertex distribution of the
        ! incoming line has no very small faces (i.e. they are removed
        ! before this routine) and therefore small faces can only 
        ! originate from/near vertices to be kept in the outgoing line. 

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMLineRefinerLB2DUDT)               :: refiner 
        type(GGTMFieldlineDataUDT), intent(in)      :: linein
        type(GGTMFieldlineDataUDT), intent(inout)   :: lineout
        integer(I8), intent(inout)                  :: vertID 
        type(GGTMDataUDT), intent(inout)            :: ggtmdata

        ! Auxiliary
        real(R8), allocatable, dimension(:)         :: tdlcv, dlv, &
            refinerdllcin, refinerdllcout, refinerdlcvin, refinerdlcvout, &
            refinerdlcvoutinit
        integer(I8), allocatable, dimension(:)      :: tvertID, &
            sortind, segvertind
        logical, allocatable, dimension(:)          :: keepind, &
            tisnodevert

        ! Loop
        integer(I8)                                 :: i, k 

        ! Project
        !========
        ! Determine line lengths in correct length measure
        call refiner%InitializeLineData(linein)
        refinerdllcin = refiner%linedllc 
        call refiner%InitializeLineData(lineout)
        refinerdllcout = refiner%linedllc

        ! Initialize vertices etc
        tvertID = [(k, k = vertID+1, vertID+size(linein%dlcv))]
        tisnodevert = [spread(.false., 1, size(linein%dlcv))]

        ! Update vertex index
        vertID = vertID + size(linein%dlcv)

        ! Interpolate the incoming line length distribution
        call Interpolate1D(linein%dlcv, refinerdlcvin, linein%dllc, refinerdllcin)

        ! Project onto the outgoing line length distribution
        refinerdlcvout = refinerdlcvin/refinerdllcin(size(refinerdllcin))*refinerdllcout(size(refinerdllcout))

        ! Append vertices to be kept
        refinerdlcvoutinit = refiner%GetLineEdgeLength(lineout)
        segvertind = lineout%GetAllSegmentVertIndices()
        refinerdlcvout = [refinerdlcvout, refinerdlcvoutinit(segvertind)]
        tvertID = [tvertID, lineout%vert(segvertind)]
        tisnodevert = [tisnodevert, lineout%isnodevert(segvertind)]
        
        ! Interpolate the outgoing actual line length
        call Interpolate1D(refinerdlcvout, tdlcv, refinerdllcout, lineout%dllc)
        !tdlcv = refinerdlcvout/refinerdlcvout(size(refinerdlcvout))*&
        !    refinerdllcout(size(refinerdllcout))

        ! Insert vertices to be kept
        allocate(sortind(size(tdlcv)))
        call Sort(tdlcv, ind=sortind, ascend=.true.)
        tvertID     = tvertID(sortind)
        tisnodevert = tisnodevert(sortind)
        deallocate(sortind)

        ! Check for very small distances
        allocate(keepind(size(tdlcv)))
        keepind = .true.
        dlv = abs(tdlcv(2:) - tdlcv(1:size(tdlcv)-1))
        do i = 1, size(dlv)
            if (dlv(i) < disttol) then 
                if (tisnodevert(i) .and. .not. tisnodevert(i+1)) then
                    keepind(i+1) = .false.
                elseif (.not. tisnodevert(i) .and. tisnodevert(i+1)) then 
                    keepind(i) = .false.
                elseif (tvertID(i) == tvertID(i+1)) then 
                    ! Same vertices, can delete one 
                    keepind(i) = .false.
                else
                    ! Can't delete nodes, small faces will be present...
                end if  
            end if 
        end do 

        ! Delete vertices to remove small distances
        tdlcv = pack(tdlcv, keepind)
        tvertID = pack(tvertID, keepind)
        tisnodevert = pack(tisnodevert, keepind)
        tdlcv(1) = 0 ! ensure start and end on line
        tdlcv(size(tdlcv)) = lineout%dllc(lineout%nl)
        tisnodevert(1) = .true.
        tisnodevert(size(tdlcv)) = .true.

        ! Reconstruct line
        !=================
        ! Add coordinates
        call lineout%AddVertexCoordinates(tdlcv)

        ! Add line vertices
        call lineout%AddVertexIDs(tvertID, tisnodevert)

        ! Update segment data
        call lineout%UpdateSegmentData(ggtmdata)
        
    end subroutine

    !------------------------------------------------------------------!
    !                             OUTPUT                               !
    !------------------------------------------------------------------!

    ! Simulation grid extraction
    subroutine ExtractSimulationGrid(simgrid, grid, magneticField, topomesh)

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
        class(TopomeshUDT), intent(in)          :: topomesh

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
        ! General
        simgrid%data%sglegacy%isClassicalGrid = 0

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
        simgrid%face%TMfacelabel    = simgrid%face%label ! assumed the same for now
        simgrid%face%reg        = grid%face%region%Get()
        simgrid%face%ntot       = nf 

        ! Cells
        simgrid%cell%vert           = grid%cell%vert%Get()
        simgrid%cell%vertP(:, 1)    = grid%cell%vp1%Get()
        simgrid%cell%vertP(:, 2)    = grid%cell%vp2%Get()
        simgrid%cell%reg            = grid%cell%region%Get()
        simgrid%cell%ntot           = nc 
        simgrid%cell%ngc            = 0

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
            0, 1, simgrid%vert%by)
        call magneticField%interp%Evaluate(simgrid%vert%x, simgrid%vert%y, &
            1, 0, simgrid%vert%bx)
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
        call ComputeTopologicalData(simgrid, topomesh)

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
            ispolygonstart, tc, isbranchingpolygon, hasbndf1, hasbndf2

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
                if (size(IDs) <= 2) then 
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
        allocate(sortind(size(tfnb, 1)), ispolygonstart(size(tfnb, 1)), &
            isbranchingpolygon(size(tfnb, 1)))
        call SortPolygonEdges(tfnb, count(tf), sortind, ispolygonstart, isbranchingpolygon)
        isbranchingpolygon = isbranchingpolygon(sortind)
        allocate(polind(count(ispolygonstart)))
        polind = pack([(k, k = 1, size(tfnb, 1))], ispolygonstart)
        polind = [polind, size(tfnb, 1)+1]
        if (count(isbranchingpolygon) > 0) then 
            ! This shouldn't happen
            print *, count(isbranchingpolygon)
            print *, pack(polind, isbranchingpolygon((polind(1:size(polind)-1))))
            !call gdErrorHandler('ComputeGridData: branching flux ' // & 
            !    'tubes detected, this is a bug')
        end if 
        

        ! Sort faces
        allocate(tfind(count(tf)))
        tfind = pack(fID, tf)
        tfind = tfind(sortind)
        tfnb(:, 1) = tfnb(sortind, 1)
        tfnb(:, 2) = tfnb(sortind, 2)

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
            fd%isclosedft(fd%nFt), hasbndf1(fd%nFt), hasbndf2(fd%nFt))
        hasbndf1 = .true.
        hasbndf2 = .true.
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
                    hasbndf1(i) = .false.
                    print *, 'ComputeGridData: flux tube: ', i, &
                        ' is not closed but the starting cell (number: ', ftc(1), &
                        ' ) has no boundary faces. Not adding face'
                elseif (size(tcf) > 1) then 
                    ! Unexpected, may be an issue
                    hasbndf1(i) = .false.
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
                    hasbndf2(i) = .false.
                    print *, 'ComputeGridData: flux tube: ', i, &
                        ' is not closed but the ending cell (number: ', ftc(1), &
                        ' ) has no boundary faces. Not adding face'
                elseif (size(tcf) > 1) then 
                    ! Unexpected, may be an issue
                    hasbndf2(i) = .false.
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
                if (size(IDs) /= 2) then 
                    fd%fluxtubefsIDs(cc, :) = 0
                else 
                    fd%fluxtubefsIDS(cc, :) = IDs 
                end if 

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

        ! Sanity check: all cells should be in a flux tube
        if (ncell /= simgrid%cell%ntot) then 
            print *, 'ComputeGridData: not all cells ' // & 
                'have been attributed to a flux tube'
        end if 

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
                if (hasbndf1(i) .and. hasbndf2(i)) then 
                    dp = (ccx*bfx(2:size(ftf)-1) + ccy*bfy(2:size(ftf)-1))
                elseif (hasbndf1(i) .and. .not. hasbndf2(i)) then 
                    dp = (ccx*bfx(2:size(ftf)) + ccy*bfy(2:size(ftf)))
                elseif (.not. hasbndf1(i) .and. hasbndf2(i)) then 
                    dp = (ccx*bfx(1:size(ftf)-1) + ccy*bfy(1:size(ftf)-1))
                else
                    dp = (ccx*bfx + ccy*bfy)
                end if 
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


        ! Housekeeping
        !=============
        end associate


    end subroutine

    ! Compute topological data
    subroutine ComputeTopologicalData(simgrid, topomesh)

        ! Description
        !============
        ! This routine adds desired topological data to the grid, such 
        ! as the X point IDs, strike point IDs, ... These are 
        ! added to the simgrid. 

        ! Note: it is currently assumed that all topomesh vertices 
        ! appear in the grid with unaltered numbering. Normally, this 
        ! should be the case (if this ever would not be the case anymore,
        ! then one can identify the vertices by doing distance checks 
        ! probably)

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(inout)            :: simgrid
        type(TopomeshUDT), intent(in)           :: topomesh

        ! Auxiliary
        integer(I8)                             :: tp, tploc, &
            ndivface, ndiv, ps, pe
        integer(I8), allocatable, dimension(:)  :: primaryxp, &
            allpoints, divind, tpbf, tpf, tflabels, divface, sortind, &
            vertdivID
        integer(I8), allocatable, dimension(:, :)   :: divfacevert

        logical, allocatable, dimension(:)      :: isdivface, &
            isbranchingpolygon, ispolygonstart

        ! Loop
        integer(I8)                             :: i, j, k

        ! Compute basic data
        !===================
        ! Determine topological mesh type
        simgrid%data%topoflag = IdentifyTopologicalMeshType(topomesh)

        ! Basic X-, O-, S-, T-point data
        simgrid%data%xpointID = topomesh%GetXPointIDs()
        simgrid%data%opointID = topomesh%GetOPointIDs()
        simgrid%data%spointID = topomesh%GetStrikePointIDs()
        simgrid%data%tpointID = topomesh%GetClosedContourTangencyPointIDs()
        simgrid%data%spointxpID = topomesh%GetStrikePointXPointIDs()
        simgrid%data%nxp = size(simgrid%data%xpointID)
        simgrid%data%nop = size(simgrid%data%opointID)
        simgrid%data%nsp = size(simgrid%data%spointID)
        simgrid%data%ntp = size(simgrid%data%tpointID)        

        ! Compute divertor data
        !======================
        ! A divertor target part is defined as all the (sorted) faces 
        ! with a label near one (or multiple) strike or tangency 
        ! points. Algorithm is simple: first, mark all face labels that
        ! are labels of a divertor plate. Then, extract the faces and
        ! sort them. Each polygon then forms a divertor target. In post-
        ! process, the mapping between strike/tangency points and 
        ! divertor targets can be made. 

        ! Initialize
        allocate(divind(topomesh%vert%ntot))
        divind = 0
        
        ! Concatenate for ease
        allpoints = [simgrid%data%spointID, simgrid%data%tpointID]

        ! Loop
        allocate(isdivface(simgrid%face%ntot))
        isdivface = .false.
        do i = 1, size(allpoints)
            ! Determine the next point
            tp = allpoints(i)

            ! Get the boundary faces of this point
            tpf = GetvertFace(simgrid%vert, tp)
            allocate(tpbf(count(simgrid%face%BF(tpf))))
            tpbf = pack(tpf, simgrid%face%BF(tpf))

            ! Take face labels
            tflabels = simgrid%face%label(tpbf)

            ! Mark faces with these labels
            do j = 1, size(tflabels)
                where (simgrid%face%label == tflabels(j)) isdivface = .true.
            end do 

            ! Housekeeping
            deallocate(tpbf)
        end do

        ! Get all divertor face vertices and sort
        ndivface = count(isdivface)
        allocate(divface(ndivface))
        divface = pack([(k, k = 1, simgrid%face%ntot)], isdivface)
        divfacevert = simgrid%face%vert(divface, :)
        allocate(sortind(ndivface), ispolygonstart(ndivface), &
            isbranchingpolygon(ndivface))
        call SortPolygonEdges(divfacevert, ndivface, sortind, ispolygonstart, &
            isbranchingpolygon)
        divfacevert = divfacevert(sortind, :)
        divface = divface(sortind)

        ! Sanity checks
        if (count(isbranchingpolygon) /= 0) then 
            call gdErrorHandler('ComputeTopologicalData: some divertor ' // & 
                'plates form branching polygons, unexpected')
        end if 

        ! Build list and pointer
        if (allocated(simgrid%data%divFcP)) deallocate(simgrid%data%divFcP)
        simgrid%data%divFc = divface
        simgrid%data%ndivFc = size(divface)
        ndiv = count(ispolygonstart)
        simgrid%data%ndiv = ndiv
        allocate(simgrid%data%divFcP(ndiv, 2))
        allocate(vertdivID(simgrid%vert%ntot))
        vertdivID = 0
        ps = 0
        pe = 0
        do i = 1, ndiv 
            ! Get polygon bounds
            ps = findloc(ispolygonstart(pe+1:), .true., 1, back=.false.)
            ps = ps + pe 
            if (i < ndiv) then 
                pe = findloc(ispolygonstart(ps+1:), .true., 1, back=.false.)
                pe = pe + ps - 1
            else
                pe = size(ispolygonstart)
            end if 

            ! Set pointer
            simgrid%data%divFcP(i, 1) = ps
            simgrid%data%divFcP(i, 2) = pe - ps + 1

            ! Set divertor plate ID for vertices
            vertdivID([divfacevert(ps:pe, 1), divfacevert(ps:pe, 2)]) = i
        end do

        ! Compute additional point data
        !==============================
        ! Additional x-point data
        primaryxp = topomesh%GetPrimaryXPointIDs()
        if (allocated(simgrid%data%isprimaryxp)) deallocate(simgrid%data%isprimaryxp)
        allocate(simgrid%data%isprimaryxp(simgrid%data%nxp))
        do i = 1, simgrid%data%nxp
            if (any(simgrid%data%xpointID(i) == primaryxp)) then 
                simgrid%data%isprimaryxp(i) = 1
            else
                simgrid%data%isprimaryxp(i) = 0 
            end if 
        end do 

        ! Additional s-, t-point data
        simgrid%data%spointdivID = vertdivID(simgrid%data%spointID)
        simgrid%data%tpointdivID = vertdivID(simgrid%data%tpointID) 

    end subroutine

    ! Label translation
    subroutine TranslateGridLabels(simgrid, topomesh, vessel, options, &
        formattype)

        ! Description
        !============
        ! This routine translates the gridding labels from the grid
        ! generator format to a format of choice (this can be e.g. 
        ! the grid deformation or solps format). Definitions of face
        ! labels and cell flags should be given in mod_definitions.F90

        ! Note: information will likely be lost during this step, which 
        ! is exactly the reason why we can't translate back... 

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(inout)    :: simgrid 
        type(VesselUDT), intent(in)     :: vessel
        type(TopomeshUDT), intent(in)   :: topomesh 
        type(GGOptionsUDT), intent(in)  :: options
        character(*), intent(in)        :: formattype ! destination format


        ! Checks
        !=======
        select case (formattype)

        case ('solps')

            ! Call translator
            call TranslateGridLabelsSOLPS(simgrid, topomesh, vessel, options)

        case ('GD')

            ! Call translator
            call TranslateGridLabelsGD(simgrid, topomesh)

        case ('no')

            ! Don't do anything

        case default
            
            call gdErrorHandler('TranslateGridLabels: format type: ' // formattype &
                // ' not implemented')

        end select

    end subroutine

    subroutine TranslateGridLabelsSOLPS(simgrid, topomesh, vessel, options)

        ! Description
        !============
        ! Not all rules are known yet, but in order to at least 
        ! simplify the amount of boundaries we get, we do the following:
        ! - face labels:
        !       Internal boundaries: labels are set to zero
        !       Non-vessel boundaries: concatenated where possible, negative label 
        !       vessel boundaries: labels are given based on original 
        !       structure (positive)
        !       Core boundaries: concatenated, negative label (random)
        ! - cell regions:
        !       Core parts: SOLPScoreregID + SOLPScoreregIDincr
        !       Other parts: all but core part value
        ! - cell flags: 
        !       Internal cells: SOLPSinternalcellID
        !       Boundary cells: SOLPSbndcellID
        ! - face regions: 
        !       zero everywhere, except on targets. There they go from 
        !       1 to the number of targets

        ! Note 1: it is assumed that the initial face labels identify the
        ! topological mesh boundary ID

        ! Declare variables
        !==================
        ! Arguments
        type(VesselUDT), intent(in)                 :: vessel
        type(GridUDT), intent(inout)                :: simgrid 
        type(TopomeshUDT), intent(in)               :: topomesh
        type(GGOptionsUDT), intent(in)              :: options

        ! Auxiliary
        integer(I8)                                 :: ne, &
            TPlabel, fcregID, mySOLPScoreregIDincr
        integer(I8), allocatable, dimension(:)      :: IFlabels, &
            TPlabels, bndlabels, fl_orig, fl_new, Clabels, &
            facelabelmapping, allfID, tfID, &
            sortindex, ind, solpslabels, psind, coreIDs, &
            cellregionmapping, veslabels, WGlabels, OFlabels, tfc, &
            allsepIDs, tfv, tfsepv, allTPlabels, uvesstructlabels, &
            reslabels
        integer(I8), allocatable                    :: edges(:, :), &
            vesstructlabels(:, :), flabels(:, :)
        real(R8), allocatable, dimension(:)         :: xf, yf
        logical, allocatable, dimension(:)          :: &
            ispolygonstart, isbranchingpolygon, islabelfound, &
            keepvert, isvesselface
        ! type(PolygonSetUDT)                         :: tempps

        ! Loop
        integer(I8)                                 :: i, j, k, flc, &
            flcinc, coreIDc, regIDc, flccoreinc, flccore

        ! Initialize
        !===========
        ! Allocate
        allocate(reslabels(0))

        ! Initialize
        isvesselface = simgrid%face%BF ! adjusted later

        ! Face label counter and face label increment
        flc = 0 
        flcinc = -1 ! we set negative face labels
        flccoreinc = SOLPScorefclblIDincr ! we set different core labels
        flccore = SOLPScorefclblID

        ! Set solps temporary labels
        solpslabels = [(k, k = 1, 3)]
        TPlabel = maxval(solpslabels)+1

        ! Get vessel vertex labels (first one, structures)
        call vessel%polygonset%GetLabels(vesstructlabels)
        
        ! Take unique value of first labels
        call Unique(vesstructlabels(:, 1), uvesstructlabels)

        ! Add labels already 
        reslabels = [reslabels, uvesstructlabels]

        ! Map
        !====
        ! Face labels
        !------------
        ! Store original face labels
        fl_orig = simgrid%face%TMfacelabel
        fl_new = fl_orig

        ! Get basic IDs
        IFlabels = topomesh%GetInternalFaceIDs()
        Clabels = topomesh%GetCoreFaceIDs()
        TPlabels = topomesh%GetTargetFaceIDs()

        ! Get derived IDs
        bndlabels = topomesh%GetBoundaryFaceIDs()
        veslabels = topomesh%GetVesselFaceIDs()
        call Setdiff(bndlabels, [veslabels, Clabels], OFlabels)
        call Setdiff(veslabels, TPlabels, WGlabels)

        ! Create (temporary) mapping
        allocate(facelabelmapping(0:maxval(fl_orig))) ! start from zero for ease
        facelabelmapping = 0
        facelabelmapping(IFlabels) = 0
        facelabelmapping(Clabels) = solpslabels(1)
        facelabelmapping(OFlabels) = solpslabels(2)
        facelabelmapping(WGlabels) = solpslabels(3)

        ! Deal with targets to comply to idiot SOLPS conventions
        allsepIDs = topomesh%GetSeparatrixFluxSurfaceIDs()
        allocate(islabelfound(topomesh%face%ntot))
        islabelfound = .false.
        do j = 1, size(TPlabels)
            ! Check if the label was found already, otherwise continue
            if (islabelfound(TPlabels(j))) then 
                cycle 
            end if 

            ! Get vertices of this face
            tfv = topomesh%face%vert(TPlabels(j), :)

            ! Get vertex/vertices with separatrix IDS 
            allocate(keepvert(size(tfv)))
            keepvert = .false.
            do i = 1, size(tfv)
                if (any(topomesh%vert%fsID(tfv(i)) == allsepIDs)) then 
                    keepvert(i) = .true.
                end if
            end do

            ! Sanity check
            if (count(keepvert) == 0) then 
                call gdErrorHandler('Translate GridLabelsSOLPS: target ' // & 
                    'boundary does not have a separatrix vertex, this ' // & 
                    'is likely a bug')
            end if 

            ! Keep only separatrix vertices
            allocate(tfsepv(count(keepvert)))
            tfsepv = pack(tfv, keepvert)

            ! Check with which other target plates this part connects 
            ! through its separatrix vertices
            do i = 1, size(TPlabels)
                ! Skip the same face
                if (i == j) then 
                    cycle
                end if

                ! Get face vertices
                tfv = topomesh%face%vert(TPlabels(i), :)

                ! Check
                if(any(tfv(1) == tfsepv) .or. any(tfv(2) == tfsepv)) then 
                    ! Check if the facelabel mapping is already initialized
                    ! for one of the faces 
                    if (facelabelmapping(TPlabels(i)) == 0 .and. facelabelmapping(TPlabels(j)) == 0) then 
                        ! Both zero, define new label
                        facelabelmapping(TPlabels(i)) = TPlabel
                        facelabelmapping(TPlabels(j)) = TPlabel

                        ! Update
                        solpslabels = [solpslabels, TPlabel]
                        TPlabel = TPlabel + 1
                        
                    elseif (facelabelmapping(TPlabels(i)) /= 0 .and. facelabelmapping(TPlabels(j)) == 0) then 
                        ! Overwrite
                        facelabelmapping(TPlabels(j)) = facelabelmapping(TPlabels(i))
                        
                    elseif (facelabelmapping(TPlabels(j)) /= 0 .and. facelabelmapping(TPlabels(i)) == 0) then 
                        ! Overwrite
                        facelabelmapping(TPlabels(i)) = facelabelmapping(TPlabels(j))

                    else 
                        ! This shouldn't happen, print warning and choose one
                        print *, 'faces: ', i, j
                        print *, 'TranslateGridLabelsSOLPS: both faces ' // & 
                            'already have a target label, this should not ' // & 
                            'happen and may be a bug. We continue and ' // & 
                            'overwrite one of the labels'
                        facelabelmapping(TPlabels(j)) = facelabelmapping(TPlabels(i))
                    end if
                end if 
            end do 

            ! Housekeeping
            deallocate(keepvert, tfsepv)
        end do 

        ! Map
        fl_new = facelabelmapping(fl_orig)
        do i = 1, simgrid%face%ntot
            if (facelabelmapping(simgrid%face%label(i)) == 0) then 
                simgrid%face%label(i) = 0
            end if 
        end do 

        ! Extract the polygon set for each set of boundaries and set 
        ! a different label for each simple polygon
        allfID = [(k, k = 1, simgrid%face%ntot)]

        ! Deal with core boundaries separately (should be 
        ! first solps label)
        ! Get number of faces
        j = 1
        ne = count(fl_new == solpslabels(j))

        ! Extract faces
        allocate(tfID(ne), edges(ne, 2), sortindex(ne), &
            ispolygonstart(ne), isbranchingpolygon(ne)) ! 
        tfID = pack(allfID, fl_new == solpslabels(j))
        edges = simgrid%face%vert(tfID, :)
        call SortPolygonEdges(edges, ne, sortindex, ispolygonstart, &
            isbranchingpolygon)
        allocate(psind(count(ispolygonstart)))

        ! Set labels for each distinct polygon piece
        psind = pack([(k, k = 1, ne)], ispolygonstart)
        psind = [psind, ne+1]
        do i = 1, count(ispolygonstart)

            ! Sanity check
            do while (any(abs(reslabels) == abs(flccore)))  
                ! Print warning
                print *, 'ID: ', flccore
                print *, 'TranslateGridLabelsSOLPS: warning: core boundary ID ' // &
                    'already used for other non-core and predefined boundary'

                ! Update counter
                flccore = flccore + flccoreinc
            end do 

            ! Add to reserved face labels
            reslabels = [reslabels, flccore]

            ! Get indices
            ind = [(k, k = psind(i), psind(i+1)-1)]

            ! Set label
            simgrid%face%label(tfID(sortindex(ind))) = flccore 

            ! Update the face label
            flccore = flccore + flccoreinc

            ! Update vessel faces
            isvesselface(tfID(sortindex(ind))) = .false.

        end do 
        deallocate(tfID, edges, sortindex, ispolygonstart, isbranchingpolygon, psind)

        ! Loop over remaining labels
        do j = 2, size(solpslabels)
            ! Get number of faces
            ne = count(fl_new == solpslabels(j))

            ! Extract faces
            allocate(tfID(ne), edges(ne, 2), sortindex(ne), &
                ispolygonstart(ne), isbranchingpolygon(ne)) ! 
            tfID = pack(allfID, fl_new == solpslabels(j))
            edges = simgrid%face%vert(tfID, :)
            call SortPolygonEdges(edges, ne, sortindex, ispolygonstart, &
                isbranchingpolygon)
            allocate(psind(count(ispolygonstart)))

            ! Set labels for each distinct polygon piece
            psind = pack([(k, k = 1, ne)], ispolygonstart)
            psind = [psind, ne+1]
            do i = 1, count(ispolygonstart)
                ! Update the face label
                flc = flc + flcinc 

                ! Avoid adding any reserved labels
                do while (any(abs(flc) == abs(reslabels)))
                    flc = flc + flcinc
                end do 

                ! Get indices
                ind = [(k, k = psind(i), psind(i+1)-1)]

                ! Set label
                simgrid%face%label(tfID(sortindex(ind))) = flc 

                ! Update vessel faces
                if (j == 2) then 
                    isvesselface(tfID(sortindex(ind))) = .false.
                end if 
            end do 
            deallocate(tfID, edges, sortindex, ispolygonstart, isbranchingpolygon, psind)
        end do 

        ! Overwrite vessel region labels
        !-------------------------------
        if (options%structurebasedlabels) then 
            ! Unpack for ease
            associate(&
                xv      => simgrid%vert%x,    &
                yv      => simgrid%vert%y,    &
                fv      => simgrid%face%vert  &
                )
            
            ! Compute face coordinates
            xf = 0.5*(xv(fv(:, 1)) + xv(fv(:, 2)))
            yf = 0.5*(yv(fv(:, 1)) + yv(fv(:, 2)))

            ! Interpolate
            call vessel%exactplfvessel%EvaluateLabel(xf, yf, flabels)

            ! Extract
            where (isvesselface) simgrid%face%label = abs(flabels(:, 1))
            
            ! Housekeeping
            end associate
        end if 

        ! Cell regions
        !-------------
        ! Get core IDs
        coreIDs = topomesh%GetCoreCellIDs()
        mySOLPScoreregIDincr = SOLPScoreregIDincr

        ! Compute mapping
        allocate(cellregionmapping(0:maxval(simgrid%cell%reg)))
        cellregionmapping = 0
        coreIDc = SOLPScoreregID
        regIDc = 1
        do i = 1, maxval(simgrid%cell%reg)
            if (any(i == coreIDs)) then 
                ! Core region
                cellregionmapping(i) = coreIDc 
                coreIDc = coreIDc + SOLPScoreregIDincr
            else
                ! Check if we should update the region ID
                !if (regIDc == coreIDc) then 
                    if (mySOLPScoreregIDincr /= 0) then 
                        if ((mod(regIDc, mySOLPScoreregIDincr)-solpscoreregID) == 0) then 
                        ! Assumed solpscoreregIDincr larger than one
                        regIDc = regIDc + 1
                        end if 
                    end if
                !end if 
                cellregionmapping(i) = regIDc 
                regIDc = regIDc + 1
            end if 
        end do 
        simgrid%cell%reg = cellregionmapping(simgrid%cell%reg)

        ! Cell flags
        !-----------
        ! Initialize to internal cell
        simgrid%cell%cflags = SOLPSinternalcellID

        ! Set boundary and internal cell flags
        do i = 1, simgrid%face%ntot 
            ! Get neighbour cells
            tfc = GetFaceCell(simgrid%face, i)

            ! Check if boundary face
            if (size(tfc) == 1) then
                simgrid%cell%cflags(tfc) = SOLPSbndcellID
            elseif (size(tfc) == 2) then 
                ! do nothing
            else
                ! Call error
                print *, 'face: ', i, 'cells: ', tfc
                call gdErrorHandler('TranslateGridLabelsSOLPS: face does ' // & 
                    'not have one or two cells, check grid interconnection')
            end if 

        end do

        ! Face regions
        !-------------
        ! Initialize
        simgrid%face%reg = 0
        allTPlabels = solpslabels(4:)
        fcregID = 0

        ! Set target flags
        do i = 1, size(allTPlabels)
            fcregID = fcregID + 1
            where (fl_new == allTPlabels(i)) simgrid%face%reg = fcregID 
        end do

        


    end subroutine

    subroutine TranslateGridLabelsGD(simgrid, topomesh)

        ! Description
        !============
        ! Translate the face labels from the initial grid to be 
        ! compatible with the grid deformation implementation. This 
        ! means that targets, core boundaries, outer boundaries, and 
        ! vessel boundaries have to be identified. It is assumed that 
        ! the input grid has the default labels attributed by the 
        ! grid generator, which point to the topological faces. 

        ! Modules
        !========
        use mod_definitions, only: targetID, coreID, vesselID, &
            outerboundaryID, interiorID

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(inout)        :: simgrid 
        type(TopomeshUDT), intent(in)       :: topomesh 

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: coreIDs, targetIDs, &
            vesselIDs, lastfsIDs, interiorIDs

        ! Loop
        integer(I8)                             :: k 

        ! Initialize
        !===========
        ! Get necessary boundary IDs from topomesh
        coreIDs = topomesh%GetCoreFaceIDs()
        targetIDs = topomesh%GetTargetFaceIDs()
        vesselIDs = topomesh%GetVesselFaceIDs()
        lastfsIDs = topomesh%GetLastFluxSurfaceFaceIDs()
        call SetDiff([(k, k = 1, maxval(simgrid%face%label))], &
            [coreIDs, targetIDs, vesselIDs, lastfsIDs], interiorIDs)

        ! Adjust labels
        !==============
        where (simgrid%face%label == targetIDs)     simgrid%face%label = targetID
        where (simgrid%face%label == coreIDs)       simgrid%face%label = coreID
        where (simgrid%face%label == vesselIDs)     simgrid%face%label = vesselID
        where (simgrid%face%label == lastfsIDs)     simgrid%face%label = outerboundaryID
        where (simgrid%face%label == interiorIDs)   simgrid%face%label = interiorID

    end subroutine

    ! Label mapping getter
    subroutine GetGridFaceLabelMappingGD(simgrid, topomesh, facelabelGG, &
        facelabelGD)

        ! Description
        !============
        ! Routine to get only the mapping without actually changing the
        ! face labels

        ! Modules
        !========
        use mod_definitions, only: targetID, coreID, vesselID, &
            outerboundaryID, interiorID

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(in)               :: simgrid 
        type(TopomeshUDT), intent(in)           :: topomesh 
        integer(I8), allocatable, dimension(:), intent(out)     :: &
            facelabelGG, facelabelGD 

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: coreIDs, targetIDs, &
            vesselIDs, lastfsIDs, interiorIDs

        ! Loop
        integer(I8)                             :: k 

        ! Initialize
        !===========
        ! Get necessary boundary IDs from topomesh
        coreIDs = topomesh%GetCoreFaceIDs()
        targetIDs = topomesh%GetTargetFaceIDs()
        vesselIDs = topomesh%GetVesselFaceIDs()
        lastfsIDs = topomesh%GetLastFluxSurfaceFaceIDs()
        call SetDiff([(k, k = 1, maxval(simgrid%face%label))], &
            [coreIDs, targetIDs, vesselIDs, lastfsIDs], interiorIDs)

        ! Construct mapping
        facelabelGG = [coreIDs, targetIDs, vesselIDs, lastfsIDs, interiorIDs]
        facelabelGD = [spread(coreID, 1, size(coreIDs)), spread(targetID, 1, size(targetIDs)), &
            spread(vesselID, 1, size(vesselIDs)), spread(outerboundaryID, 1, size(lastfsIDs)), &
            spread(interiorID, 1, size(interiorIDs))]
        

    end subroutine

    !------------------------------------------------------------------!
    !                          DIAGNOSTICS                             !
    !------------------------------------------------------------------!

    ! Main diagnostics driver
    subroutine RunGridDiagnostics(simgrid)

        ! Description
        !============
        ! This routine runs some diagnostics on the grid, which print 
        ! out useful information for debugging/checking where the grid
        ! is of poor quality. Currently, these entail:
        ! - face intersections

        ! Declare variables
        !==================
        ! Arguments
        class(GridUDT), intent(in)          :: simgrid 

        ! Diagnostics
        !============
        ! Face intersections
        call CheckFaceIntersections(simgrid, 'gg_faceintersections')

    end subroutine
    
    ! Face intersections
    subroutine CheckFaceIntersections(simgrid, writefile) 

        ! Description
        !============
        ! This routine checks if non-adjacent faces of the grid 
        ! intersect (adjacent faces will always intersect in the common
        ! vertex). Normally, two faces can only intersect at one point
        ! unless they are coinciding. These faces are flagged as well.
        ! The list of face IDs and the intersection coordinates is 
        ! written out in a file.

        ! Declare variables
        !==================
        ! Arguments
        class(GridUDT), intent(in)              :: simgrid 
        character(*), intent(in)                :: writefile

        ! Auxiliary
        real(R8)                                :: xint, yint
        real(R8), allocatable, dimension(:)     :: allxint, &
            allyint 
        integer(I8), allocatable, dimension(:)  :: faceIDs

        ! Loop 
        integer(I8)                             :: i, j 

        ! Initialize
        !===========
        ! Unpack
        associate(&
            fv      => simgrid%face%vert,   &
            xv      => simgrid%vert%x,      &
            yv      => simgrid%vert%y,      &
            nf      => simgrid%face%ntot    &
            )
        
        ! Initialize
        allocate(allxint(0), allyint(0), faceIDs(0))

        ! Compute
        !========
        !$omp parallel do default(private) shared(simgrid, allxint, allyint, faceIDs) schedule(dynamic)
        do i = 1, nf 
            do j = i+1, nf 
                ! Compmute segment-segment intersection
                call SegmentIntersections(xint, yint, xv(fv(i, 1)), &
                    yv(fv(i, 1)), xv(fv(i, 2)), yv(fv(i, 2)), xv(fv(j, 1)), &
                    yv(fv(j, 1)), xv(fv(j, 2)), yv(fv(j, 2)))

                if (.not. isnan(xint)) then 
                    if (any(fv(i, 1) == fv(j, :)) .or. any(fv(i, 2) == fv(j, :))) then 
                        ! Don't do anything - neighbouring face
                    else 
                        !$omp critical 
                        allxint = [allxint, xint, xint]
                        allyint = [allyint, yint, yint]
                        faceIDs = [faceIDs, i, j]
                        !$omp end critical
                    end if 
                end if
            end do 
        end do 
        !$omp end parallel do

        ! Write
        !======
        ! Write to file
        call WriteVertexData(faceIDs, allxint, allyint, writefile)
        
        ! Write to terminal
        if (size(faceIDs) > 0) then 
            print *, 'CheckFaceIntersections: some faces are intersecting, ' // & 
                'IDs and intersection coordinates are written in ' // &
                ' "' // writefile // '".'
        end if 


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
        type(GGTMFieldlineDataUDT)              :: line 
        class(TopomeshUDT), intent(in)          :: topomesh 
        class(GGTMDataUDT)                      :: ggtmdata 

        ! Auxiliary
        integer(I8)                             :: startv, endv
        integer(I8), allocatable, dimension(:)  :: bndf, bndv, allsegIDs
        logical, allocatable, dimension(:)      :: isnotfound 

        ! Loop
        integer(I8)                             :: i

        ! Initialize
        !===========
        ! Associate
        associate(&
            vert        => topomesh%vert,   &
            vertdata    => ggtmdata%vert,   &
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
                erfv => topomesh%face%vert(tmcell%erf, :))

            ! Check field values
            if (vert%fval(srfv(1)) > vert%fval(srfv(2))) then 
                startv = srfv(1)
            else 
                startv = srfv(2)
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
            else 
                startv = srfv(2)
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
                line = vertdata(bndv(1))%line
                !call line%Initialize(vert%x(bndv), vert%y(bndv), [vert%fsID(bndv(1))], &
                !    [1, 1], [0_I8])
                !call line%AddVertexCoordinates([0.0_R8])
                !call line%AddVertexIDs(bndv)
                !call line%UpdateLineData(topomesh, ggtmdata)

            else 
                call gdErrorHandler('ExtractTMCellAlignedBoundary: expected ' // & 
                    'a single boundary vertex but found multiple/none')
            end if
            
            ! Return 
            return 
        end if 

        ! Construct line
        !===============
        ! Initialize
        allocate(isnotfound(size(bndf)))
        isnotfound = .true. 

        ! Get all segments (each face has a line, but that should 
        ! consist of a single segment)
        allocate(allsegIDs(size(bndf)))
        do i = 1, size(bndf)
            if (facedata(bndf(i))%line%ns /= 1) then 
                ! We could support this in the future though
                call gdErrorHandler('ExtractTMCellAlignedBoundary: ' // & 
                    'topomesh face has not exactly one segment, not supported')
            end if 
            allsegIDs(i) = facedata(bndf(i))%line%segID(1)
        end do 

        ! Construct the line
        call line%Initialize(ggtmdata, allsegIDs)

        ! Check orientation
        if (line%vert(1) == startv .and. line%vert(line%nv) == endv) then 
            ! Correctly sorted, nothing to do
        elseif (line%vert(1) == endv .and. line%vert(line%nv) == startv) then
            ! Need to flip
            call line%Flip()
        else
            ! Something wrong
            call gdErrorHandler('ExtractTMCellAlignedBoundary: could not ' // &
                'find start or end vertex after extracting line')
        end if  

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

        ! 'segments'
        ! <nseg, total number of segments>
        ! 'segment coordinates'
        ! 'segment <nseg>' <repeated for each segment>
        ! <vertID, xv, yv>
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
        ! 'cell tubes hfline <nc>' 
        ! <vertID, x, y>
        ! 'cell tubes lfline <nc>'
        ! <vertID, x, y>

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
        logical, allocatable, dimension(:)      :: doface, docell, dolines, &
            dotubes
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

        allocate(docell(size(c)), dolines(size(c)), dotubes(size(c)))
        docell = .false.
        dolines = .false.
        dotubes = .false.
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
                if (allocated(c(i)%tubes)) then 
                    dotubes(i) = .true.
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

        ! Cell tubes
        write (fu, *) 'cell tubes'
        do i = 1, size(c) 
            if (dotubes(i)) then 
                ! Write header
                write (fu, *) 'cell ', i 

                ! Write high field line
                do j = 1, size(c(i)%tubes)
                    write (fu, *) 'hfline tube ', j
                    write (fu, *) c(i)%tubes(j)%hfline%nv
                    do k = 1, c(i)%tubes(j)%hfline%nv
                        write (fu, *) c(i)%tubes(j)%hfline%vert(k), &
                            c(i)%tubes(j)%hfline%xv(k), c(i)%tubes(j)%hfline%yv(k)
                    end do 

                    ! Write low field line
                    write (fu, *) 'lfline tube ', j
                    write (fu, *) c(i)%tubes(j)%lfline%nv
                    do k = 1, c(i)%tubes(j)%lfline%nv
                        write (fu, *) c(i)%tubes(j)%lfline%vert(k), &
                            c(i)%tubes(j)%lfline%xv(k), c(i)%tubes(j)%lfline%yv(k)
                    end do 
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

    ! TM refinement data writing
    subroutine WriteTopologicalMeshLineRefinementData(ggtmdata, &
        filename)

        ! Description
        !============
        ! This routine writes out the refinement data into a file 
        ! specified by 'filename' into the output directory. This file 
        ! may be modified by the user and read in in order to have more
        ! control over the refinement. Structure is as follows:

        ! <header> 
        ! 'celldata'
        ! <cell%ntot>
        ! 'ID, doBLstart, doBLend, ncBLstart, ncBLend, dlBLlengthbased'
        ! <the above for each cell per line>
        ! 'dlBLstart' (in order of celldata, one row per cell)
        ! <dlBLstart>
        ! 'dlBLend' (in order of celldata, one row per cell)
        ! <dlBLend>
        ! 'tubedata'
        ! <cell%ntot>
        ! 'ID, doBLstart, doBLend, ncBLstart, ncBLend'
        ! <the above for each cell per line>
        ! 'dlBLstart' (in order of celldata, one row per cell)
        ! <dlBLstart>
        ! 'dlBLend' (in order of celldata, one row per cell)
        ! <dlBLend>

        ! Declare variables
        !==================
        ! Modules 
        use mod_plotter 
        use mod_specialchars, only : filesepchar

        ! Arguments
        class(GGTMDataUDT), intent(in)          :: ggtmdata
        character(*), intent(in)                :: filename 

        ! Auxiliary
        integer                                 :: fu, lengthbased 
        character(:), allocatable               :: dir

        ! Loop
        integer(I8)                             :: i

        ! Initialize
        !===========
        ! Unpack
        associate(&
            t       => ggtmdata%tube,   & 
            c       => ggtmdata%cell    &
        )

        ! Construct writing directory
        dir = plotdir // filesepchar // filename // '.dat'

        ! Open file
        open (action='write', file=trim(dir), newunit=fu, &
             status='unknown')

        ! Write header
        write(fu, *) 'VERSION3.00.00'

        ! Write scalar cell data
        !=======================
        ! Write header
        write(fu, *) 'celldata'

        ! Write total number of cells
        write(fu, *) size(c)

        ! Write scalar data
        write(fu, *) 'ID, doBLstart, doBLend, ncBLstart, ncBLend, dlBLlengthbased'
        do i = 1, size(c)
            associate(r     => c(i)%linerefoptions)
            if (r%dlBLlengthbased) then 
                lengthbased = 1
            else
                lengthbased = 0
            end if 
            write(fu, *) i, r%doBLstart, r%doBLend, r%ncBLstart, r%ncBLend, lengthbased
            end associate
        end do

        ! write array cell data
        !======================
        ! Write boundary layer data
        !--------------------------
        ! Write dlBLstart
        write (fu, *) 'dlBLstart'
        do i = 1, size(c)
            ! Write header
            write(fu, *) c(i)%linerefoptions%dlBLstart
        end do

        ! Write dlBLstart
        write (fu, *) 'dlBLend'
        do i = 1, size(c)
            ! Write header
            write(fu, *) c(i)%linerefoptions%dlBLend
        end do

        ! Write scalar tube data
        !=======================
        ! Write header
        write(fu, *) 'tubedata'

        ! Write total number of tubes
        write(fu, *) size(t)

        ! Write scalar data
        write(fu, *) 'ID, doBLstart, doBLend, ncBLstart, ncBLend, lengthbased'
        do i = 1, size(t)
            associate(r     => t(i)%linerefoptions)
            if (r%dlBLlengthbased) then 
                lengthbased = 1
            else
                lengthbased = 0
            end if 
            write(fu, *) i, r%doBLstart, r%doBLend, r%ncBLstart, r%ncBLend, lengthbased
            end associate
        end do

        ! write array tube data
        !======================
        ! Write boundary layer data
        !--------------------------
        ! Write dlBLstart
        write (fu, *) 'dlBLstart'
        do i = 1, size(t)
            ! Write header
            write(fu, *) t(i)%linerefoptions%dlBLstart
        end do

        ! Write dlBLstart
        write (fu, *) 'dlBLend'
        do i = 1, size(t)
            ! Write header
            write(fu, *) t(i)%linerefoptions%dlBLend
        end do
        
        ! Housekeeping
        !=============
        end associate
        close(fu)

    end subroutine

    ! TM refinement data reading
    subroutine ReadTopologicalMeshLineRefinementData(ggtmdata, &
        filename)

        ! Description
        !============
        ! This routine reads the refinement data from a file 
        ! specified by 'filename' from the output directory. This file 
        ! may be modified by the user and read in in order to have more
        ! control over the refinement. Structure is as follows:

        ! <header> 
        ! 'celldata'
        ! <cell%ntot>
        ! 'ID, doBLstart, doBLend, ncBLstart, ncBLend'
        ! <the above for each cell per line>
        ! 'dlBLstart' (in order of celldata, one row per cell)
        ! <dlBLstart>
        ! 'dlBLend' (in order of celldata, one row per cell)
        ! <dlBLend>

        ! Declare variables
        !==================
        ! Modules 
        use mod_plotter 
        use mod_specialchars, only : filesepchar
        use mod_readwrite
        use mod_inputfileparser

        ! Arguments
        class(GGTMDataUDT), intent(inout)       :: ggtmdata
        character(*), intent(in)                :: filename 

        ! Auxiliary
        integer                                 :: fu, lengthbased
        integer(I8)                             :: nc, cID, nt
        real(R8), allocatable, dimension(:)     :: temp
        character(:), allocatable               :: thisline
        logical                                 :: reachedeof

        ! Loop
        integer(I8)                             :: i

        ! Initialize
        !===========
        ! Unpack
        associate(&
            t       => ggtmdata%tube,   &
            c       => ggtmdata%cell    &
        )

        ! Construct reading directory

        ! Open file
        open (action='read', file=trim(filename), newunit=fu, &
             status='unknown')

        ! Read cell data
        !===============
        ! Read header
        call ReadUntilFound(fu, 'celldata', reachedeof)
        if (reachedeof) then 
            call gdErrorHandler('ReadTopologicalMeshCellLineRefinementData: ' // & 
                'could not find cell data in file')
        end if 

        ! Read number of cells - just for check, since celldata should 
        ! be allocated already
        read(fu, *) nc 

        ! Check
        if (nc /= size(c)) then 
            call gdErrorHandler('ReadTopologicalMeshCellLineRefinementData: ' // & 
                'cell data does not have the same dimensions as current ' // & 
                'cell data, check input')
        end if 

        
        ! Read scalar cell data
        !----------------------
        ! Read header
        call ReadSingleLine(fu, thisline, reachedeof)

        ! Read scalar data
        do i = 1, size(c)
            associate(r     => c(i)%linerefoptions)
            read(fu, *) cID, r%doBLstart, r%doBLend, r%ncBLstart, r%ncBLend, lengthbased 
            if (lengthbased > 0) then 
                c(i)%linerefoptions%dlBLlengthbased = .true.
            else
                c(i)%linerefoptions%dlBLlengthbased = .false.
            end if 
            end associate
        end do

        ! Read array cell data
        !---------------------
        ! Read dlBLstart
        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, size(c)
            allocate(temp(c(i)%linerefoptions%ncBLstart))
            read(fu, *) temp 
            c(i)%linerefoptions%dlBLstart  = temp 
            deallocate(temp)
        end do

        ! Read dlBLstart
        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, size(c)
            allocate(temp(c(i)%linerefoptions%ncBLend))
            read(fu, *) temp 
            c(i)%linerefoptions%dlBLend = temp 
            deallocate(temp)
        end do

        ! Read tube data
        !===============
        ! Read header
        call ReadUntilFound(fu, 'tubedata', reachedeof)
        if (reachedeof) then 
            call gdErrorHandler('ReadTopologicalMeshCellLineRefinementData: ' // & 
                'could not find tube data in file')
        end if 

        ! Read number of tubes - just for check, since celldata should 
        ! be allocated already
        read(fu, *) nt

        ! Check
        if (nt /= size(t)) then 
            call gdErrorHandler('ReadTopologicalMeshCellLineRefinementData: ' // & 
                'tube data does not have the same dimensions as current ' // & 
                'tube data, check input')
        end if 

        
        ! Read scalar tube data
        !----------------------
        ! Read header
        call ReadSingleLine(fu, thisline, reachedeof)

        ! Read scalar data
        do i = 1, size(t)
            associate(r     => t(i)%linerefoptions)
            read(fu, *) cID, r%doBLstart, r%doBLend, r%ncBLstart, r%ncBLend, lengthbased
            if (lengthbased > 0) then 
                r%dlBLlengthbased = .true.
            else
                r%dlBLlengthbased = .false.
            end if 
            end associate
        end do

        ! Read array tube data
        !---------------------
        ! Read dlBLstart
        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, size(t)
            allocate(temp(t(i)%linerefoptions%ncBLstart))
            read(fu, *) temp 
            t(i)%linerefoptions%dlBLstart  = temp 
            deallocate(temp)
        end do

        ! Read dlBLstart
        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, size(t)
            allocate(temp(t(i)%linerefoptions%ncBLend))
            read(fu, *) temp 
            t(i)%linerefoptions%dlBLend = temp 
            deallocate(temp)
        end do
        
        ! Housekeeping
        !=============
        end associate
        close(fu)

    end subroutine


end module 