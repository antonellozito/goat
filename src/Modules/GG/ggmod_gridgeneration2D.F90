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
        ! xv, yv and the field values at these vertices, fv.
        real(R8), allocatable, dimension(:)     :: xv, yv, fv

    end type 

    ! Field line data
    type :: GGTMFieldlineDataUDT

        real(R8), allocatable, dimension(:)         :: xv, yv, xl, yl
        integer(I8), allocatable, dimension(:)      :: vert

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
        integer(I8)                                 :: srflabel, erflabel, &
            srf, erf 
        integer(I8), allocatable, dimension(:)      :: srfvert, erfvert, &
            hffaces, lffaces, hfvert, lfvert

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

        class(RealDynamicArrayUDT), allocatable     :: x, y

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

        class(IntegerDynamicArrayUDT), allocatable  :: v1, v2, label

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

        class(IntegerDynamicArrayUDT), allocatable  :: v, vp1, vp2, &
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

        ! Initial grid structure
        !vert = struct('x', zeros(0, 1), 'y', zeros(0, 1), 'BV', zeros(0, 1), 'ntot', 0);
        !face = struct('vert', zeros(0, 2) , 'ntot', 0, 'labels', zeros(0, 1), 'region', zeros(0, 1));
        !cell = struct('vert', zeros(0, 1), 'vertP', zeros(0, 2), 'region', zeros(0, 1), 'ntot', 0, 'nvert', 0);
        !fs = struct('ntot', 0);
        !grid = struct('vert', vert, 'face', face, 'cell', cell, 'fs', fs);

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

        case default 

            call gdErrorHandler('GenerateUnstructuredAlignedGrid: ' // & 
                'unknown distribution method: ' // options%ggmethod)

        end select

        ! Distribute the vertices
        !call DistributeVerticesTopologicalMeshCells(ggtmdata, topomesh, &
        !    topomesh, poloidalvertexdistributor, magneticFieldDF, )

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
    !subroutine IndependentGridder(ggtmdata, topomesh, )

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
                nfl = size(facedata(tf(j))%fv)
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
            indsrf, tf, cind, nc, ntf, incr, nv
        integer(I8), allocatable, dimension(:)  :: tubec, tubef, tcf, &
            tcv, tcfv1, tcfv2, hffaces, lffaces, hfvert, lfvert, &
            allIDs, s1, s2
        integer(I8), allocatable, dimension(:, :)   :: nint, segrf, &
            segc, vertexID
        real(R8)                                :: hfval, lfval, &
            dhf1, dhf2, dlf1, dlf2
        real(R8), allocatable, dimension(:)     :: tcvfval, tcfv1val, &
            tcfv2val, tx, ty, xl, yl, tfval, sr1, sr2, txint, tyint
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
                celldata(tubec(k))%srf = srf 
                indsrf = findloc(tcf, srf, 1)
                if (indsrf == 0) then 
                    ! This shouldn't be happening
                    call gdErrorHandler('AddTopologicalMeshCellGriddingData: ' // & 
                        'expected cell to have radial faces as defined in ' // & 
                        'flux tube, but face is not present in cell faces. ' // &
                        'Check topological mesh for inconsistencies.')
                end if 

                ! Set end radial face
                erf = tubef(j+1) ! should work 
                celldata(tubec(k))%erf = erf
                inderf = findloc(tcf, erf, 1)
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

                ! Determine high and low field poloidal faces (unsorted)
                tcfv1val = vert%fval(tcfv1)
                tcfv2val = vert%fval(tcfv2)
                allocate(ishfface(size(tcf)), islfface(size(tcf)))
                ishfface = .false. 
                islfface = .false. 
                do k = 1, size(tcf)
                    ! Check if the face is of poloidal/sep type
                    if (any(face%type(tcf(k)) == [TMfacepolID, TMfacesepID])) then 
                        dhf1 = abs(tcfv1val(tcf(k)) - hfval)
                        dhf2 = abs(tcfv2val(tcf(k)) - hfval)
                        dlf1 = abs(tcfv1val(tcf(k)) - lfval)
                        dlf2 = abs(tcfv2val(tcf(k)) - lfval)
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

                ! Sanity checks
                if (count(.not. ishfface .and. .not. islfface) /= 2) then 
                    call gdErrorHandler('AddTopologicalMEshCellGriddingData: ' // & 
                        'cell has not exactly two radial faces, this is not yet supported')
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

            ! Trace contours
            !---------------
            ! Get initial vertex distribution
            tf = tubedata(i)%distributionface
            tx = facedata(i)%xv
            ty = facedata(i)%yv

            ! Make sure we trace from high to low field value
            allocate(tfval(size(tx)))
            call magneticField%interp%Evaluate(tx, ty, 0, 0, tfval)
            if (tfval(1) < tfval(size(tfval))) then 
                tx = tx(size(tx):1:-1)
                ty = ty(size(ty):1:-1)
            end if 
            deallocate(tfval)

            ! Trace
            tempc = fieldtracer%TraceContours(tx, ty)

            ! Clean
            call CleanContours(tempc)

            ! Check if contours make sense and reformat if necessary
            allIDs = tempc%ID
            allocate(iscontourfound(size(tx))) ! normally, IDs go from 1 to number of points
            allocate(keepind(size(tempc)))
            iscontourfound = .false. 
            keepind = .true. 
            if (tube%isclosed(i)) then 
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
                            
                            ! Check if the contour is open, should be the case
                            if ((c1%x(1) == c1%x(size(c1%x))) .and. (c1%y(1) == c1%y(size(c1%y)))) then 
                                call gdErrorHandler('AddTopologicalMeshCellGriddingData :' // & 
                                    'two open contours found that form a closed contour, ' // &
                                    'may be a bug in the contouring algorithm')
                            end if

                            ! Mark as found
                            iscontourfound(c1%ID) = .true. 

                            ! Mark for removal
                            keepind(j) = .false. 
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
            ! Initialize
            nc = size(tempc)
            ntf = size(tubef)
            allocate(xintda(nc, ntf), yintda(nc, ntf), nint(nc, ntf))
            nint = 0

            ! Convert to polygon
            do j = 1, nc
                call polc(j)%Construct(tempc(j)%x, tempc(j)%y)
            end do 

            ! Loop over all contours
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
                if (any(nint(j, :) > 1)) then 
                    ! Simply set warning message
                    isintersectremoved = .true. 
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
                    do k = 1, ntf 
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
                end if 
            end do 
            deallocate(keepind)

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
                    xl = [xint(j, k), tempc(j)%x(segc(j, k)+2:segc(j, k+1)-1:incr), xint(j, k+1)]
                    yl = [yint(j, k), tempc(j)%y(segc(j, k)+2:segc(j, k+1)-1:incr), yint(j, k+1)]

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
                vertexID)
        end do 
        
        ! Housekeeping
        !=============
        end associate


    end subroutine 

    !------------------------------------------------------------------!
    !                           AUXILIARY                              !
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
        select case (storagetype)

        case ('standard')

            vert%x = ConstructRealDynamicArray()
            vert%y = ConstructRealDynamicArray()

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
        select case (storagetype)

        case ('standard')

            cell%v = ConstructIntegerDynamicArray()
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

        ! Note: the facedata is assumed to hold vertex distributions 
        ! (xv, yv) for all faces. However, since they may be 
        ! overwritten at later stages, vertex indices are not 
        ! propagated in this routine (these may be not initialized
        ! yet...)

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMCellDataUDT)                  :: tmcell 
        character(*), intent(in)                :: loc 
        class(GGTMFieldlineDataUDT)             :: line 
        class(TopomeshUDT), intent(in)          :: topomesh 
        class(GGTMDataUDT)                      :: ggtmdata 

        ! Auxiliary
        real(R8), allocatable, dimension(:)     :: tx, ty 
        integer(I8)                             :: startv, endv, &
            thisf, thisfind, nextv
        integer(I8), allocatable, dimension(:)  :: bndf, bndv 
        logical                                 :: doflip, facefound
        logical, allocatable, dimension(:)      :: isnotfound 

        ! Loop
        integer(I8)                             :: j 

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
                erfv => topomesh%face%vert(tmcell%erf, :))

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
                line%xl = vert%x(bndv)
                line%yl = vert%y(bndv)
                line%xv = line%xl
                line%yv = line%yl
                line%vert = bndv

                ! Initialize rest
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

        ! Get the starting aligned face
        thisfind = findloc(face%vert(bndf, 1), startv, 1)
        doflip = .false. 
        if (thisfind == 0) then 
            thisfind = findloc(face%vert(bndf, 2), startv, 1)
            doflip = .true. 
        end if 
        if (thisfind == 0) then 
            ! Probably something wrong upstream when determining cell data
            call gdErrorHandler('ExtractTMCellAlignedBoundary: could not ' // & 
                'find face with starting vertex, check cell data')
        end if 

        ! Set as found
        isnotfound(thisfind) = .false.

        ! Add face data
        thisf   = bndf(thisfind)
        line%xl = face%x(thisf)%Get()
        line%yl = face%y(thisf)%Get()
        line%xv = facedata(thisf)%xv
        line%yv = facedata(thisf)%yv

        if (doflip) then 
            line%xl = line%xl(size(line%xl):1:-1)
            line%yl = line%yl(size(line%yl):1:-1)
            line%xv = line%xv(size(line%xv):1:-1)
            line%yv = line%yv(size(line%yv):1:-1)
        end if 

        ! Get the next vertex
        if (doflip) then 
            nextv = face%vert(thisf, 1)
        else
            nextv = face%vert(thisf, 2)
        end if
        
        ! Loop over remaining faces
        do while (any(isnotfound))
            
            ! Check if the next vertex is the end vertex
            if (nextv == endv) then 

                ! Sanity check
                if (any(isnotfound)) then 
                    call gdErrorHandler('ExtractTMCellAlignedBoundary: ' // & 
                        'already found end vertex while not all faces ' // & 
                        'were found, this is a bug')
                end if
            end if 

            ! Get the next face & vertex
            j = 1
            doflip = .false. 
            facefound = .true.
            do while (j <= size(bndf))
                if (isnotfound(j)) then 
                    if (face%vert(bndf(j), 1) == nextv) then 
                        ! Set
                        thisfind = j 
                        thisf = bndf(j)
                        nextv = face%vert(thisf, 2)
                        facefound = .true.

                        ! Exit
                        exit 
                    elseif (face%vert(bndf(j), 2) == nextv) then 
                        ! Set
                        thisfind = j 
                        thisf = bndf(j)
                        nextv = face%vert(thisf, 1)
                        facefound = .true. 
                        
                        ! Exit
                        exit 
                    else
                        ! Update j
                        j = j + 1
                    end if 
                end if 
            end do 

            ! Sanity check
            if (.not. facefound) then 
                call gdErrorHandler('ExtractTMCellAlignedBoundary: ' // & 
                    'could not find a next face, this is a bug ')
            end if 

            ! Add face data 
            if (doflip) then 
                tx = face%x(thisf)%Get()
                ty = face%y(thisf)%Get()
                line%xl = [line%xl(1:size(line%xl)-1), tx(size(tx):1:-1)] ! avoid double coordinates
                line%yl = [line%yl(1:size(line%yl)-1), ty(size(ty):1:-1)]
                line%xv = [line%xv(1:size(line%yl)-1), &
                    facedata(thisf)%xv(size(facedata(thisf)%xv):1:-1)]
                line%yv = [line%yv(1:size(line%yl)-1), &
                    facedata(thisf)%yv(size(facedata(thisf)%xv):1:-1)]
            
            else
                line%xl = [line%xl(1:size(line%xl)-1), face%x(thisf)%Get()] ! avoid double coordinates
                line%yl = [line%yl(1:size(line%yl)-1), face%y(thisf)%Get()]
                line%xv = [line%xv(1:size(line%yl)-1), facedata(thisf)%xv]
                line%yv = [line%yv(1:size(line%yl)-1), facedata(thisf)%yv]
            end if 

            ! Set
            isnotfound(thisf) = .false.

        end do 

        ! Housekeeping
        !=============
        end associate

    end subroutine


end module 