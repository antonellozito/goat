!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module provides functionality to determine the 2D topological 
! mesh that is defined by the vessel and magnetic field. Although the 
! topological mesh is also a grid, it does not inherit from the GridUDT
! structure, since this structure typically contains too much or even
! irrelevant data. This does imply some duplication of some basic 
! routines, but that shouldn't be a big issue. The topological mesh can
! be used to generate an unstructured grid later on based on the 
! cells and flux tubes determined here. 

module ggmod_topology2D

    ! Load modules
    use mod_precision
    use mod_errorhandler
    use mod_dynamicarrays
    use mod_contour2D
    use mod_polygon
    use mod_sort
    use mod_definitions
    use mod_linearsolverinterface, only: SolveDenseLinearSystemDI
    use goatmod_types, only : magneticFieldUDT, VesselUDT, &
        ConstructVesselPolygonSet
    use goatmod_userinput, only : TopomeshOptionsUDT
    use Interpolant1D
    use mod_streamlinetracing2D
    use PolygonLevelsetFunction2D
    implicit none
    private 
    public :: TopomeshUDT, ConstructTopologicalMesh, TraceExtrema2D, &
        TraceTangencyPoints2D, ReadTopologicalMesh, WriteTopologicalMesh, &
        IdentifyTopologicalMeshType

    ! Module parameters
    real(R8), parameter, private        :: tprelfieldtol = 1e-14 ! relative field tolerance under which extrema are removed
    real(R8), parameter, private        :: disttol = 1e-12 ! distance tolerance
    real(R8), parameter, private        :: distfrac = 1e-3 ! distance fraction when removing edge vertices

    !==================================================================!
    !                                                                  !
    !                            TYPES                                 !
    !                                                                  !
    !==================================================================!

    ! Topological vertex 
    type :: TopomeshVertUDT 

        ! Description
        !============
        ! Vertex structure for the topological mesh. Contains different 
        ! fields to identify the type of vertex and to navigate the mesh
        ! more easily. 

        integer(I8)                     :: ntot ! total number of vertices
        integer(I8), allocatable        :: ID(:), face(:), cell(:), &
            flags(:), type(:), faceP(:, :), cellP(:, :), fsID(:)
        real(R8), allocatable           :: x(:), y(:), fval(:)
        logical, allocatable            :: BV(:)

    contains 

        ! Initializer
        procedure :: Initialize     => InitializeTopologicalMeshVertex

        ! Getters
        procedure :: GetFace        => GetTMVertFace 
        procedure :: GetFaceNeig    => GetTMVertFaceNeig

    end type 

    ! Topological face
    type :: TopomeshFaceUDT

        ! Description
        !============
        ! Face structure for the topological mesh. Contains different 
        ! fields to identify the type of vertex and to navigate the mesh
        ! more easily. Note that contrary to classical grid faces, these
        ! faces are not a straight line between vessel vertices, but 
        ! contain a (simple) polygon that describes the spatial layout
        ! of the face. 

        integer(I8)                     :: ntot, nfc ! total number of vertices
        integer(I8), allocatable        :: ID(:), vert(:, :), cell(:), &
            label(:), fsID(:), type(:), cellP(:, :)
        type(RealDynamicArrayUDT), allocatable  :: x(:), y(:)
        type(PolygonUDT), allocatable   :: pol(:)
        logical, allocatable            :: BF(:)

    contains 
    
        ! Initializer
        procedure :: Initialize     => InitializeTopologicalMeshFace

        ! Getters
        procedure :: GetCell        => GetTMFaceCell

    end type

    ! Topological cell
    type :: TopomeshCellUDT

        ! Description
        !============
        ! Topological mesh cell

        integer(I8)                     :: ntot, ntotv, ntotf ! total number of vertices
        integer(I8), allocatable        :: ID(:), vert(:), vertP(:, :), &
            face(:), faceP(:, :), flags(:), tube(:)
    contains 

        ! Initializer
        procedure :: Initialize     => InitializeTopologicalMeshCell

        ! Getters
        procedure :: GetVert        => GetTMCellVert
        procedure :: GetFace        => GetTMCellFace 

    end type

    ! Topological tubes
    type :: TopomeshTubeUDT

        ! Description
        !============
        ! Topological mesh tube object. This contains cells that are in 
        ! the same flux tube. The faces that are in this flux tube are 
        ! topological mesh faces that have the tube cells as boundaries.
        ! Note that these tubes can only be constructed as a 
        ! post-processing step after topological mesh construction. The
        ! cell faces and cells present in the 'cell' and 'face' structure
        ! are sorted after construction. Additionally, we store the 
        ! aligned faces and vertices from both sides in separate arrays,
        ! in this case in bndface1 (with pointer bndface1P), bndface2, &
        ! bndvert1, bndvert2. It should be noted that these faces and 
        ! vertices are not sorted in any particular direction (this is 
        ! non-trivial, since tube boundaries may be branching polygons)

        integer(I8)                                 :: ncell, nface, ntot 
        integer(I8), allocatable, dimension(:)      :: cell, face,  &
            bndf1, bndf2, bndv1, bndv2, ftneig1, ftneig2
        integer(I8), allocatable, dimension(:, :)   :: cellP, faceP, &
            bndf1P, bndf2P, bndv1P, bndv2P, ftneig1P, ftneig2P
        logical, allocatable, dimension(:)          :: isclosed 

    contains
    
        ! Initializer
        procedure :: Initialize     => InitializeTopologicalMeshTube

        ! Deallocator
        procedure :: Deallocate     => DeallocateTopologicalMeshTube

        ! Getter
        procedure :: GetFace        => GetTMTubeFace
        procedure :: GetCell        => GetTMTubeCell
        procedure :: GetBndFace     => GetTMTubeBndFace
        procedure :: GetBndVert     => GetTMTubeBndVert
        procedure :: GetNeig        => GetTMTubeNeig

    end type 

    ! General topological mesh type
    type :: TopomeshUDT 

        type(TopomeshVertUDT)   :: vert 
        type(TopomeshFaceUDT)   :: face 
        type(TopomeshCellUDT)   :: cell 
        type(TopomeshTubeUDT)   :: tube 
        integer(I8)             :: nfs 
        type(IntegerDynamicArrayUDT)    :: fsID
        type(RealDynamicArrayUDT)       :: fsfval

    contains 

        ! Initialize
        procedure :: Initialize         => InitializeTopologicalMesh 

        ! Getters
        procedure :: GetInternalFaceIDs, GetBoundaryFaceIDs, &
            GetTargetFaceIDs, GetSeparatrixFaceIDs, GetVesselFaceIDs, &
            GetCoreFaceIDs, GetCoreCellIDs, GetWideGridCellIDs, &
            GetSeparatrixFluxSurfaceIDs, GetStrikePointIDs, GetXPOintIDs, &
            GetOPointIDs, GetPrimaryXPointIDs, GetStrikePointXPointIDs, &
            GetClosedContourTangencyPointIDs, GetLastFluxSurfaceFaceIDs
    end type 

    contains 

    !==================================================================!
    !                                                                  !
    !                           ROUTINES                               !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                  TOPOLOGICAL MESH CONSTRUCTION                   !
    !------------------------------------------------------------------!

    ! Main constructor
    subroutine ConstructTopologicalMesh(vessel, magneticField, options, &
        topomesh, fieldtracer, vesseltracer, streamlinetracer)

        ! Description
        !============
        ! This routine generates the topological mesh based on the given magnetic
        ! field and vessel geometry. 

        ! Algorithm
        !==========
        ! 1) Compute points of interest (minima, maxima, saddle points, tangency
        ! points) 
        ! 2) Create a reduced topological mesh of the extrema to generate the cuts
        ! (should contain vertices and faces, not cells) - this forms the initial
        ! basis of the final topological mesh
        ! 3) Find the vessel boundary, introduce the tangency points and compute
        ! additional contours starting from those tangency points.
        ! 4) Compute all intersections between the current faces, introduce these
        ! as new vertices (and keep track of type). Intersections in end points are
        ! assumed to be known already and not added twice
        ! 5) Construct the topological mesh cells

        ! Modules
        !========
        use mod_structured2Dgridding

        ! Declare variables
        !==================
        ! Arguments
        type(TopomeshUDT)                       :: topomesh
        type(VesselUDT), intent(in)             :: vessel
        type(magneticFieldUDT), intent(in)      :: magneticField 
        type(TopomeshOptionsUDT), intent(in)    :: options
        class(ContourTracerUDT), allocatable, intent(inout)  :: vesseltracer, fieldtracer
        class(StreamlineTracerUDT), intent(in)  :: streamlinetracer

        ! Auxiliary
        type(VesselUDT)                         :: tmvessel

        ! Initialize
        !===========
        ! Copy the vessel structure since this may change due to 
        ! reconstruction of the vessel polygons etc while constructing
        ! the topological mesh (also small geometrical features might 
        ! change, but that should be negligible)
        tmvessel = vessel

        ! Construct basic mesh
        !=====================
        ! Check if we construct from scratch or load in a file
        if (options%readexistingTM) then 
            ! Read
            call ReadTopologicalMesh(topomesh, options%TMfilepath)

            ! Update the field tracer
            call UpdateTracersFromTopomesh(topomesh, fieldtracer, &
                magneticField, tmvessel, options)
        else
            ! Construct from scratch
            call ConstructBasicTopologicalMesh(tmvessel, magneticField, options, &
                topomesh, fieldtracer, vesseltracer, streamlinetracer)
        end if

        ! Apply adaptations
        !==================
        if (options%doadaptations) then 
            call ModifyTopologicalMesh(tmvessel, magneticField, options, &
                topomesh, fieldtracer, vesseltracer, streamlinetracer)
        end if 

        ! Write
        !======
        call WriteTopologicalMesh(topomesh, 'topomesh')

    end subroutine

    ! Basic topological mesh constructor
    subroutine ConstructBasicTopologicalMesh(vessel, magneticField, options, &
        topomesh, fieldtracer, vesseltracer, streamlinetracer)

        ! Description
        !============
        ! This routine generates the topological mesh based on the given magnetic
        ! field and vessel geometry. 

        ! Algorithm
        !==========
        ! 1) Compute points of interest (minima, maxima, saddle points, tangency
        ! points) 
        ! 2) Create a reduced topological mesh of the extrema to generate the cuts
        ! (should contain vertices and faces, not cells) - this forms the initial
        ! basis of the final topological mesh
        ! 3) Find the vessel boundary, introduce the tangency points and compute
        ! additional contours starting from those tangency points.
        ! 4) Compute all intersections between the current faces, introduce these
        ! as new vertices (and keep track of type). Intersections in end points are
        ! assumed to be known already and not added twice
        ! 5) Construct the topological mesh cells

        ! Modules
        !========
        use mod_structured2Dgridding

        ! Declare variables
        !==================
        ! Arguments
        type(TopomeshUDT)                       :: topomesh
        type(VesselUDT), intent(inout)          :: vessel
        type(magneticFieldUDT), intent(in)      :: magneticField 
        type(TopomeshOptionsUDT), intent(in)    :: options
        class(ContourTracerUDT), allocatable, intent(inout)  :: vesseltracer, fieldtracer
        class(StreamlineTracerUDT), intent(in)  :: streamlinetracer

        ! Auxiliary
        real(R8)                                :: dxfracmin, dyfracmin
        real(R8), allocatable, dimension(:)     :: xb, yb, xps, &
            yps, xg, yg, Vf, xgv, ygv
        real(R8), parameter                     :: emptyR8(0)= 0
        real(R8), allocatable, dimension(:)     :: xtp, ytp, Ftp, &
            xe, ye, fe
        integer(I8)                             :: nv, ntp 
        integer(I8), allocatable, dimension(:)  :: typee, IDs
        integer(I8), parameter                  :: emptyI8(0) = 0
        type(VesselUDT)                         :: newvessel

        ! Loop
        integer(I8)                             :: i, j

        ! Initialize
        !===========
        ! Initialize data
        dxfracmin = 1e-4_R8
        dyfracmin = 1e-4_R8

        ! Initialize topomesh 
        call topomesh%Initialize()

        ! Determine domain bounds based on vessel and magnetic field extent
        call vessel%plfvessel%ps%GetVertices(xps, yps)
        xb = [minval([xps, magneticField%interp%xgv]), maxval([xps, magneticField%interp%xgv])]
        yb = [minval([yps, magneticField%interp%ygv]), maxval([yps, magneticField%interp%ygv])]

        ! Construct a 2D structured grid for tracing (may be extended
        ! in the future for different grid types)
        nv = options%vresx*options%vresy
        allocate(xg(nv), yg(nv), xgv(options%vresx), ygv(options%vresy))
        call Construct2DStructuredUniformGrid(xg, yg, xgv, ygv, xb, yb, &
            options%vresx,  options%vresy, 0.0_R8, 0.0_R8)

        ! Compute extrema
        !================
        ! Only for mesh refinement later on
        call TraceExtrema2D(xe, ye, fe, typee, fieldtracer, magneticField, &
            options%fdonewton)

        ! Add extrema & related boundaries
        !=================================
        ! Tangency points & vessel geometry
        call AddTopologicalMeshTangencyPoints2(topomesh, vesseltracer, &
            vessel, magneticField, options)

        ! Construct refined grid based on tangency points and extrema
        ntp = count((topomesh%vert%type == TMvertextp1ID) .or. (topomesh%vert%type == TMvertextp2ID))
        allocate(xtp(ntp), ytp(ntp), Ftp(ntp), IDs(ntp))
        xtp = pack(topomesh%vert%x, &
            (topomesh%vert%type == TMvertextp1ID) .or. (topomesh%vert%type == TMvertextp2ID))
        ytp = pack(topomesh%vert%y, &
            (topomesh%vert%type == TMvertextp1ID) .or. (topomesh%vert%type == TMvertextp2ID))
        Ftp = pack(topomesh%vert%fval, &
            (topomesh%vert%type == TMvertextp1ID) .or. (topomesh%vert%type == TMvertextp2ID))
        IDs = pack(topomesh%vert%ID, &
            (topomesh%vert%type == TMvertextp1ID) .or. (topomesh%vert%type == TMvertextp2ID))
        call ConstructRefined2DStructuredGrid(xg, yg, xgv, ygv, xb, yb, &
            options%vresx, options%vresy, [xtp, xe], [ytp, ye], 5, 5, dxfracmin, dyfracmin)  

        ! Evaluate magnetic field and vessel
        allocate(Vf(size(xg)))
        call magneticField%interp%Evaluate(xg, yg, 0, 0, Vf)
        call magneticField%interp%Evaluate(xtp, ytp, 0, 0, Ftp)

        ! Update the field tracer 
        fieldtracer = ConstructStructuredTracer(&
            reshape(Vf, [size(xgv), size(ygv)]), xgv, ygv, &
            xtp, ytp, Ftp, IDs, fieldtracer%npmin, fieldtracer%npmax, fieldtracer%dl)
        
        !! Visualize by tracing contours
        !resc = 100
        !dv = (maxval(Vf) - minval(Vf))
        !cgv = [(k, k = 0, resc)]*(dv*0.90_R8)/real(resc, kind=R8) + minval(Vf) + dv*0.05
        !contours = fieldtracer%TraceContours(cgv)
        !allocate(pcontours(size(contours)))
        !!$omp parallel do default(shared)
        !do k = 1, size(contours)
        !    call pcontours(k)%Construct(contours(k)%x, contours(k)%y)
        !end do 
        !!$omp end parallel do
        !call tempps%Construct(pcontours)
        !call tempps%WriteData('mfcontours')

        ! Extrema
        call AddTopologicalMeshExtrema(topomesh, fieldtracer, &
            magneticField, options)

        ! Do temporary writing
        call WriteTopologicalMesh(topomesh, 'topomesh_beforecells')

        ! Remove parts that do not lie inside the vessel
        newvessel = vessel
        call TrimTopologicalMesh(topomesh, magneticField, newvessel)

            ! Do temporary writing
        call WriteTopologicalMesh(topomesh, 'topomesh_beforecells')

        ! Process points
        !===============
        ! Set extrema with nearly identical values to be the same value
        do i = 1, topomesh%vert%ntot-1
            if (topomesh%vert%type(i) == TMvertexsaddleID) then 
                do j = i+1, topomesh%vert%ntot 
                    if (topomesh%vert%type(j) == TMvertexsaddleID) then 
                        if (abs(topomesh%vert%fval(i)-topomesh%vert%fval(j)) < options%ffieldtol) then 
                            topomesh%vert%fval(j) = topomesh%vert%fval(i)
                            topomesh%vert%fsID(j) = topomesh%vert%fsID(i)
                        end if 
                    end if 
                end do 
            end if 
        end do 

        ! Compute necessary contours
        !===========================
        ! Add contours and intersections
        call AddTopologicalMeshContours(topomesh, magneticField, newvessel, &
            fieldtracer, options)

        ! Vertex faces (preliminary, for garbage tangency point removal)
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Do temporary writing
        call WriteTopologicalMesh(topomesh, 'topomesh_beforecells')

        ! Remove garbage tangency points
        call RemoveGarbageTangencyPoints(topomesh)

        ! Do temporary writing
        call WriteTopologicalMesh(topomesh, 'topomesh_beforecells')

        ! Eliminate limiter-like configurations
        !======================================
        ! Eliminate
        call RemoveTopologicalMeshLimiterRegions(topomesh, magneticField, &
            fieldtracer, streamlinetracer, options)

        ! Simplify
        call SimplifyTopologicalMeshFaces(topomesh)

        ! Do temporary writing
        call WriteTopologicalMesh(topomesh, 'topomesh_beforecells')

        ! Add necessary data
        !===================
        ! Vertex faces
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Data 
        call AddTopologicalMeshData(topomesh)

        ! Add cells
        call AddTopologicalMeshCells(topomesh)        

        ! Compute interconnection data
        call AddTopologicalMeshInterconnectionData(topomesh)

        ! Write
        !======
        ! Overwrite vessel structure
        vessel = newvessel

        ! Write topological mesh
        call WriteTopologicalMesh(topomesh, 'topomesh_base')

    end subroutine

    ! Main modifier
    subroutine ModifyTopologicalMesh(vessel, magneticField, options, &
        topomesh, fieldtracer, vesseltracer, streamlinetracer)

        ! Description
        !============
        ! This routine contains all possible modifications that may be
        ! applied to a basic topological mesh. Though it is not recommended, one might 
        ! try to apply these modifications to an already modified mesh 
        ! if desired. No guarantees on result though (most of these 
        ! operations are irreversible)

        ! Declare variables
        !==================
        ! Arguments
        type(TopomeshUDT)                       :: topomesh
        type(VesselUDT), intent(inout)          :: vessel
        type(magneticFieldUDT), intent(in)      :: magneticField 
        type(TopomeshOptionsUDT), intent(in)    :: options
        class(ContourTracerUDT), allocatable, intent(inout)  :: vesseltracer, fieldtracer
        class(StreamlineTracerUDT), intent(in)  :: streamlinetracer

        ! Compute additional contours
        !============================
        ! Core boundaries? 
        if (options%addcoreboundaries) then 
            call AddTopologicalMeshCoreBoundaries(topomesh, magneticField, &
                vessel, fieldtracer, options)
        end if 

        ! 'PF' boundaries?
        if (options%addPFboundaries) then 
            call AddTopologicalMeshPFBoundaries(topomesh, magneticField, &
                vessel, fieldtracer, options)
        end if

        ! Simplify
        call SimplifyTopologicalMeshFaces(topomesh)

        ! Do temporary writing
        call WriteTopologicalMesh(topomesh, 'topomesh_beforecells')

        ! Compute additional interconnnection data
        !=========================================
        ! Vertex faces
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Data 
        call AddTopologicalMeshData(topomesh)

        ! Add cells
        call AddTopologicalMeshCells(topomesh)

        ! Recompute interconnection data
        call AddTopologicalMeshInterconnectionData(topomesh)

        ! Do temporary writing
        call WriteTopologicalMesh(topomesh, 'topomesh_beforecells')

        ! Merge tubes?
        if (options%mergetangencypointtubes) then 
            call MergeTopologicalMeshFluxTubes(topomesh, magneticField, &
                vessel, fieldtracer, options)
        end if

        ! Remove parts if desired
        if (options%removecoreregions) then 
            call RemoveTopologicalMeshCoreRegions(topomesh)
        end if 
        if (options%removewidegridregions) then 
            call RemoveTopologicalMeshWideGridRegions(topomesh)
        end if 
        if (options%removenoncoreregions) then 
            call RemoveTopologicalMeshNonCoreRegions(topomesh)
        end if

        

        ! Recompute interconnection data
        call AddTopologicalMeshInterconnectionData(topomesh)

    end subroutine

    ! Extrema addition
    subroutine AddTopologicalMeshExtrema(topomesh, fieldtracer, &
        magneticField, options)

        ! Description
        !============
        ! This routine is largely based on TraceExtrema2DBox, but adds the extrema
        ! and their connection lines to the topology mesh as vertices and faces.
        ! Note that the location of the extrema can be refined using the Newton
        ! solver, but that the connecting faces are simply taken as the gradient
        ! lines and are therefore less accurate. For a description on how the
        ! extrema are computed, see TraceExtrema2DBox. 

        ! IMPORTANT: it is assumed that the boundary faces are already
        ! present in the topological mesh! These are used to determine
        ! which 'radial' faces, resulting from the dfdx, dfdy contours, 
        ! should be retained in the topological mesh. Only those that 
        ! do not intersect the boundary are retained. This means that
        ! only radial faces that start and end in an extremum without
        ! intersecting the vessel boundary will be retained. 

        ! Notes
        !======
        ! Note 1: the contour tracing algorithm used to trace dFdx = 0 and dFdy = 0
        ! contours is assumed to be able to deal with saddle points and should
        ! return only 'simple' polygons. Saddle points may exist even in the dFdx =
        ! 0 and dFdy = 0 fields. Therefore, we've implemented our own tracing
        ! routine in TraceContourLineStructured2D. 

        ! Note 2: when computing the optimum with the Newton solver, it is assumed
        ! that no damping strategy is necessary since the initial point should be
        ! (very) close to the optimum. 

        ! Note 3: we pass the field tracer to reuse the grid already
        ! used before. However, we need to initialize the tracers for
        ! dfdx and dfdy again, since 1) the values change, 2) the 
        ! saddle points of the original field should not be present 
        ! (these are in fact often no saddle points of the derivatives!)

        ! Note 4: only segments of dfdx and dfdy that start or end in a
        ! topological mesh extremum (saddle, min, max) are added to the 
        ! topological mesh

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT), intent(inout)       :: topomesh
        class(ContourTracerUDT), intent(in)     :: fieldtracer 
        type(MagneticFieldUDT), intent(in)      :: magneticField 
        type(TopomeshOptionsUDT), intent(in)    :: options 

        ! Auxiliary
        integer(I8)                             :: ngp, nfxc, nfyc, &
            nx, nvinit
        integer(I8), allocatable, dimension(:)  :: ts1, ts2, tt, typee, &
            vf1, vf2, teid, tsid, sortind, temps1, temps2
        real(R8)                                :: tempx, tempy
        real(R8), allocatable, dimension(:)     :: xg, yg, f, fx, fy, &
            tx, ty, thiseig, tf, tfxx, tfxy, tfyy, xe, ye, fe, tsr1, tsr2, &
            tsrid, tempxint, tempyint
        logical                                 :: conv, donewton, addsegment
        type(RealDynamicArrayUDT)               :: xc, yc, fc
        type(RealDynamicArrayUDT), allocatable  :: xfrda(:), yfrda(:), &
            fxpsrid(:), fypsrid(:) 
        type(IntegerDynamicArrayUDT), allocatable   :: fxpeid(:), &
            fxpsid(:), fypeid(:), fypsid(:)
        type(IntegerDynamicArrayUDT)            :: tc
        type(ContourUDT)                        :: tempc
        type(ContourUDT), allocatable           :: fxc(:), fyc(:)
        type(PolygonUDT)                        :: temppol
        type(PolygonUDT), allocatable           :: fxp(:), fyp(:)
        class(ContourTracerUDT), allocatable    :: fxtracer, fytracer
        type(PolygonSetUDT)                     :: tempps

        ! Loop
        integer(I8)                             :: i, j, k, ec

        ! Initialize
        !===========
        ! Check inputs
        donewton = options%fdonewton
        if (.not. donewton) then 
            print *, 'TraceExtrema2DBox: not applying Newton solver, ' // & 
                'locations of extrema may be slightly inaccurate'
        end if

        ! Extract grid coordinates
        call fieldtracer%GetCoordinates(xg, yg)
        ngp = size(xg)

        ! Initialize dynamic arrays
        xc = ConstructRealDynamicArray()
        yc = ConstructRealDynamicArray()
        fc = ConstructRealDynamicArray()
        tc = ConstructIntegerDynamicArray()

        ! Initialize extrema counter
        ec = 0

        ! Store initial amount of vertices
        nvinit = topomesh%vert%ntot ! to update intersection IDs later on

        ! Trace contours
        !===============
        ! Evaluate derivatives
        allocate(f(ngp), fx(ngp), fy(ngp))
        call magneticField%interp%Evaluate(xg, yg, 0, 0, f)
        call magneticField%interp%Evaluate(xg, yg, 1, 0, fx)
        call magneticField%interp%Evaluate(xg, yg, 0, 1, fy)

        ! Construct new tracers
        fxtracer = fieldtracer
        fytracer = fieldtracer 
        if (allocated(fxtracer%xs)) then 
            deallocate(fxtracer%xs, fxtracer%ys, fxtracer%vs, &
                fxtracer%order, fxtracer%IDs)
        end if 
        allocate(fxtracer%xs(0), fxtracer%ys(0), fxtracer%vs(0), &
            fxtracer%order(0), fxtracer%IDs(0))
        if (allocated(fytracer%xs)) then 
            deallocate(fytracer%xs, fytracer%ys, fytracer%vs, &
            fytracer%order, fytracer%IDs)
        end if 
        allocate(fytracer%xs(0), fytracer%ys(0), fytracer%vs(0), &
            fytracer%order(0), fytracer%IDs(0))
        call fxtracer%SetValues(fx)
        call fytracer%SetValues(fy)
        

        ! Trace 
        fxc = fxtracer%TraceContours([0.0_R8])
        fyc = fytracer%TraceContours([0.0_R8])

        ! Clean
        call CleanContours(fxc)
        call CleanContours(fyc)

        ! Compute all intersections
        !==========================
        ! Convert all to polygons
        nfxc = size(fxc)
        nfyc = size(fyc)
        allocate(fxp(nfxc), fyp(nfyc))

        do i = 1, nfxc 
            call fxp(i)%Construct(fxc(i)%x, fxc(i)%y)
        end do 
        do i = 1, nfyc 
            call fyp(i)%Construct(fyc(i)%x, fyc(i)%y)
        end do 

        ! Write out data
        call tempps%Construct(fxp)
        call tempps%WriteData('extrema_fx_lines')
        call tempps%Construct(fyp)
        call tempps%WriteData('extrema_fy_lines')
        

        ! Initialize intersection trackers
        allocate(fxpeid(nfxc), fxpsid(nfxc), fypeid(nfyc), fypsid(nfyc), &
            fxpsrid(nfxc), fypsrid(nfyc))
        do i = 1, nfxc 
            fxpeid(i) = ConstructIntegerDynamicArray()
            fxpsid(i) = ConstructIntegerDynamicArray()
            fxpsrid(i) = ConstructRealDynamicArray()
        end do 
        do i = 1, nfyc 
            fypeid(i) = ConstructIntegerDynamicArray()
            fypsid(i) = ConstructIntegerDynamicArray() 
            fypsrid(i) = ConstructRealDynamicArray()
        end do 

        ! Compute intersections
        !$omp parallel do default(none) schedule(dynamic) collapse(2) &
        !$omp shared(fxp, fyp, magneticField, fxpeid, fxpsid, fxpsrid, &
        !$omp fypeid, fypsid, fypsrid, xc, yc, fc, tc, ec, nfxc, nfyc) & 
        !$omp private(i, j, tx, ty, ts1, ts2, tsr1, tsr2, nx, tt, k, &
        !$omp tempx, tempy, conv, tf, tfxx, tfxy, tfyy, thiseig)
        do i = 1, nfxc
            ! Only compute intersections with other polygons
            do j = 1, nfyc
                ! Compute intersections with next polygon
                call PolygonIntersections(fxp(i), fyp(j), tx, ty, ts1, ts2, &
                    tsr1, tsr2)
                
                ! Check if found
                nx = size(tx)
            
                ! Add intersections
                allocate(tt(nx))
                do k = 1, nx 
                    ! Refine
                    call TinyNewtonSolver(tempx, tempy, conv, &
                        tx(k), ty(k), magneticField)
                     
                    ! Check 
                    if (conv) then 
                        tx(k) = tempx
                        ty(k) = tempy
                    else
                        print *, 'TraceExtrema2DBox: newton solver did ' // & 
                            'not converge, taking original estimate'
                    end if
                end do 

                ! Compute value at location and second order derivatives
                allocate(tf(nx), tfxx(nx), tfxy(nx), tfyy(nx))
                call magneticField%interp%Evaluate(tx, ty, 0, 0, tf)
                call magneticField%interp%Evaluate(tx, ty, 2, 0, tfxx)
                call magneticField%interp%Evaluate(tx, ty, 1, 1, tfxy)
                call magneticField%interp%Evaluate(tx, ty, 0, 2, tfyy)

                ! Determine order
                do k = 1, nx

                    ! Compute value and check type of extremum
                    thiseig = ComputeEigenvaluesSymmetric2by2Matrix(&
                        tfxx(k), tfyy(k), tfxy(k))
                    if (any(thiseig == 0)) then 
                        print *, 'TraceExtrema2DBox: extrema with zero ' // & 
                            'eigenvalue detected, may be improperly identified'
                    end if 
                    if (all(thiseig > 0)) then 
                        ! Local minimum
                        tt(k) = -1;
                    elseif (all(thiseig < 0)) then 
                        ! Local maximum
                        tt(k) = 1;
                    else
                        ! Saddle point
                        tt(k) = 0;
                    end if

                    ! Update counter
                    !$omp critical
                    ec = ec + 1

                    ! Store intersection data
                    call fxpeid(i)%Append(ec)
                    call fxpsid(i)%Append(ts1(k))
                    call fxpsrid(i)%Append(tsr1(k))
                    call fypeid(j)%Append(ec)
                    call fypsid(j)%Append(ts2(k))
                    call fypsrid(j)%Append(tsr2(k))
                    !$omp end critical
                end do 

                ! Append
                !$omp critical
                call xc%Append(tx)
                call yc%Append(ty)
                call fc%Append(tf)
                call tc%Append(tt)
                !$omp end critical

                ! Housekeeping
                deallocate(tt, tf, tfxx, tfyy, tfxy)
            end do 
        end do
        !$omp end parallel do 
        
        ! Add extrema
        !============
        ! Extract
        xe = xc%Get()
        ye = yc%Get()
        fe = fc%Get()
        typee = tc%Get()

        ! Set vertex type
        do i = 1, size(typee)
            if (typee(i) == -1_I8) then 
                typee(i) = TMvertexminID 
            elseif (typee(i) == 0_I8) then 
                typee(i) = TMvertexsaddleID 
            else
                typee(i) = TMvertexmaxID 
            end if 
        end do 

        ! Loop to add
        do i = 1, ec 
            ! Update ID
            topomesh%nfs = topomesh%nfs+1

            ! Add vertex
            call AddTopologicalMeshVertex(topomesh, xe(i), ye(i), fe(i), &
                typee(i), topomesh%nfs)
            
        end do

        ! Reconstruct faces
        !==================
        ! Loop over fxp
        do i = 1, nfxc 
            ! Unpack
            teid = fxpeid(i)%Get()
            tsid = fxpsid(i)%Get()
            tsrid = fxpsrid(i)%Get()
            allocate(sortind(size(tsrid)))
            call Sort(tsrid, ind=sortind)
            tsid = tsid(sortind)
            teid = teid(sortind)
            deallocate(sortind)

            ! Check
            if (size(teid) == 0) then 
                ! No intersections whatsoever, so no saddle points at 
                ! start and end - continue without adding
                cycle
            end if 
            
            ! Extract the faces from this polygon
            call ExtractTopologicalFacesFromPolygon(fxp(i), teid, &
                tsid, xe, ye, vf1, vf2, xfrda, yfrda)

            ! Insert into topological mesh as contour (need to rebuild contour)
            tempc = fxc(i)
            do j = 1, size(xfrda)
                ! Check
                if ((vf1(j) == 0) .or. (vf2(j) == 0)) then 
                    ! Skip
                    cycle
                end if 

                ! Check if we intersect a boundary face
                call temppol%Construct(xfrda(j)%Get(), yfrda(j)%Get())
                addsegment = .true. 
                do k = 1, topomesh%face%ntot
                    if (topomesh%face%type(k) == TMfacebndID) then 
                        ! Check for intersections
                        call PolygonIntersections(topomesh%face%pol(k), &
                            temppol, tempxint, tempyint, temps1, temps2)
                        if (size(tempxint) /= 0) then 
                            ! Skip and exit
                            addsegment = .false. 
                            exit 
                        end if 
                    end if 
                end do 

                ! Check if we add
                if (addsegment) then 
                    tempc%x = xfrda(j)%Get()
                    tempc%y = yfrda(j)%Get()
                    tempc%startsaddle = vf1(j) + nvinit ! per definition non-zero, so this is fine
                    tempc%endsaddle = vf2(j) + nvinit
                    call InsertTopologicalMeshContour(topomesh, magneticField, &
                        tempc, TMfaceradID, 0)
                end if 
            end do 

            ! Add to faces
            !do j = 1, size(vf1)
            !    call AddTopologicalMeshFace(topomesh, [vf1(j), vf2(j)], &
            !        xfrda(j), yfrda(j), TMfaceradID, 0)
            !end do 
        end do 

        ! Loop over fyp
        do i = 1, nfyc 
            ! Unpack
            teid = fypeid(i)%Get()
            tsid = fypsid(i)%Get()
            tsrid = fypsrid(i)%Get()
            allocate(sortind(size(tsrid)))
            call Sort(tsrid, ind=sortind)
            tsid = tsid(sortind)
            teid = teid(sortind)
            deallocate(sortind)

            ! Check
            if (size(teid) == 0) then 
                ! No intersections whatsoever, so no saddle points at 
                ! start and end - continue without adding
                cycle
            end if 
            
            ! Extract the faces from this polygon
            call ExtractTopologicalFacesFromPolygon(fyp(i), teid, &
                tsid, xe, ye, vf1, vf2, xfrda, yfrda)

            ! Insert into topological mesh as contour (need to rebuild contour)
            tempc = fyc(i)
            do j = 1, size(xfrda)
                ! Check
                if ((vf1(j) == 0) .or. (vf2(j) == 0)) then 
                    ! Skip
                    cycle
                end if 

                ! Check if we intersect a boundary face
                call temppol%Construct(xfrda(j)%Get(), yfrda(j)%Get())
                addsegment = .true. 
                do k = 1, topomesh%face%ntot
                    if (topomesh%face%type(k) == TMfacebndID) then 
                        ! Check for intersections
                        call PolygonIntersections(topomesh%face%pol(k), &
                            temppol, tempxint, tempyint, temps1, temps2)
                        if (size(tempxint) /= 0) then 
                            ! Skip and exit
                            addsegment = .false. 
                            exit 
                        end if 
                    end if 
                end do 

                ! Check if we need to add 
                if (addsegment) then 
                    tempc%x = xfrda(j)%Get()
                    tempc%y = yfrda(j)%Get()
                    tempc%startsaddle = vf1(j) + nvinit ! per definition non-zero, so this is fine
                    tempc%endsaddle = vf2(j) + nvinit
                    call InsertTopologicalMeshContour(topomesh, magneticField, &
                        tempc, TMfaceradID, 0)
                end if 
            end do 

        end do 


    end subroutine 

    ! Tangency points and vessel addition
    subroutine AddTopologicalMeshTangencyPoints(topomesh, tpx, tpy, tpf, &
            tptype, bndtracer)

        ! Description
        !============
        ! This routine adds the tangency points and the vessel segments they
        ! construct to the topological mesh. To this end, the tangency points
        ! should already be traced using the TraceTangencyPoints2DBox routine. The
        ! vessel polygon is retraced (using resolution resx, resy) and the tangency
        ! points are then inserted at the nearest edge. Note that this may result
        ! in slight vessel geometry changes. 

        ! Algorithm
        !==========
        ! Actually it is quite easy: we add the tangency point as 'saddle' points,
        ! and make sure the value of those saddle points is exactly zero. When
        ! tracing the field line starting from one of the tangency points, the
        ! TraceContourLineSTructured2DPoint routine will return a single contour
        ! with all tangency points properly included. To make sure we trace each
        ! separate contour, we loop until all tangency points have been found. Easy
        ! peasy lemon squeazy

        ! Note: it is assumed that the tangency points were already 
        ! added to the boundary tracer beforehand! 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT), intent(inout)   :: topomesh 
        class(ContourTracerUDT), intent(in) :: bndtracer
        real(R8), intent(in)                :: tpx(:), tpy(:), tpf(:)
        integer(I8), intent(in)             :: tptype(:)

        ! Auxiliary
        integer(I8)                         :: vcinit, ttype, ttp, &
            ftype, facevert(1:2), fsID
        logical, allocatable, dimension(:)  :: notfound 
        type(ContourUDT), allocatable       :: tc(:)
        type(RealDynamicArrayUDT)           :: tx, ty 

        ! Loop
        integer(I8)                 :: i 

        ! Initialize
        !===========
        ! Counters (relative to current state of topomesh)
        vcinit = topomesh%vert%ntot

        ! Add tangency points and construct mapping
        do i = 1, size(tpx)
            if (tptype(i) == 1) then 
                ttype = TMvertextp1ID 
            else
                ttype = TMvertextp2ID
            end if 
            call AddTopologicalMeshVertex(topomesh, tpx(i), tpy(i), &
                tpf(i), ttype, topomesh%nfs+1)
            topomesh%nfs = topomesh%nfs + 1
        end do

        ! Loop
        !=====
        allocate(notfound(size(tpx)))
        notfound = .true. 
        do while (any(notfound)) 
            ! Find first tangency point to trace from
            ttp = findloc(notfound, .true., 1)

            ! Sanity check
            if (ttp == 0) then 
                call gdErrorHandler('AddTopologicalMeshTangencyPoints: bug detected')
            end if 
            
            ! Trace the contour
            tc = bndtracer%TraceContours([tpx(ttp)], [tpy(ttp)])
            
            ! Add faces & set notfound
            do i = 1, size(tc)
                ! Check
                if ((tc(i)%startsaddle == 0) .or. (tc(i)%endsaddle == 0)) then 
                    call gdErrorHandler('AddTopologicalMeshTangencyPoints: ' // & 
                        'contour lines do not start and end in tangency ' // & 
                        'points for vessel, check input')
                end if 
                
                ! Make sure no segments are added that close upon themselves -
                ! shouldn't be possible at this stage, only in some 
                ! very exceptional cases
                if (tc(i)%startsaddle == tc(i)%endsaddle) then 
                    print *, 'AddTopologicalMeshTangencyPoints: segment detected that closes upon itself. Removing...'
                    cycle
                end if 

                ! Add the face
                ftype = TMfacebndID ! vessel boundary -> outer boundary
                facevert = [tc(i)%startsaddle, tc(i)%endsaddle];
                notfound(facevert) = .false.
                facevert = facevert + vcinit; ! adjust here to get the correct global vertex index!
                fsID = 0; ! no ID, because vessel boundary
                tx = ConstructRealDynamicArray(tc(i)%x)
                ty = ConstructRealDynamicArray(tc(i)%y)
                call AddTopologicalMeshFace(topomesh, facevert, tx, ty, ftype, fsID, 0.0_R8)
            end do 
        end do



    end subroutine

    ! 2D tangency points
    subroutine AddTopologicalMeshTangencyPoints2(topomesh, &
        boundarytracer, vessel, magneticField, options)

        ! Description
        !============
        ! This routine traces the tangency points on a prescribed boundary.
        ! We base ourselves solely on the
        ! discrete representation of the vessel geometry and compute the tangency
        ! points by looking at the field evaluated in the vertices of the vessel
        ! polygon(s). Here, minima and maxima are quite easily and rapidly found
        ! without having to compute any intersections between curves. To determine
        ! the type of tangency point (i.e. whether the curve bends inside or
        ! outside the domain), we check whether it's a minima or maxima along the
        ! vessel curve and whether the tangency point field value is found in the
        ! interior of the domain near the tangency point location. 

        ! The latter is a bit tricky here: since we don't explicitly know the
        ! magnetic field in the interior, we evaluate the derivative in the
        ! tangency point location and base ourselves on that. This, however, may be
        ! inaccurate since the tangency point location itself will not correspond
        ! with the location of the continuous optimum. Caution is therefore advised
        ! when the local magnetic field varies strongly locally (or when the mesh
        ! is too coarse). 


        ! Notes
        !======
        ! Note 1: the contour tracing algorithm used to trace dFdx = 0 and dFdy = 0
        ! contours is assumed to be able to deal with saddle points and should
        ! return only 'simple' polygons. Saddle points may exist even in the dFdx =
        ! 0 and dFdy = 0 fields. Therefore, we've implemented our own tracing
        ! routine in TraceContourLineStructured2D. 

        ! Note 2: it is assumed that the vessel polygon forms a closed (or a set of
        ! closed) surfaces. These polygons are oriented later on such that the face
        ! normal points inward the domain at all times. 


        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT), intent(inout)       :: topomesh
        class(ContourTracerUDT), intent(in)     :: boundarytracer
        type(MagneticFieldUDT), intent(in)      :: magneticField 
        type(TopomeshOptionsUDT), intent(in)    :: options
        type(VesselUDT), intent(in)             :: vessel

        ! Auxiliary
        integer(I8)                             :: flag, ntp
        integer(I8), allocatable                :: dvals(:), ddvals(:), &
            tv(:), extrlocind(:), tt(:), vf1(:), vf2(:), eID(:), sID(:), &
            sortind(:)
        real(R8)                                :: fdifftol 
        real(R8), allocatable, dimension(:)     :: val, dval, tx, ty, &
            tf, nxpe, nype, nxp, nyp, normprod, dFdx, dFdy
        logical                                 :: hasbeendeleted
        logical, allocatable                    :: extrloc(:), keepind(:), &
            doflip(:)
        type(RealDynamicArrayUDT), allocatable  :: xf(:), yf(:)
        type(ContourUDT), allocatable           :: bndcontours(:)
        type(PolygonUDT), allocatable           :: bndpol(:), realbndpol(:)
        type(PolygonSetUDT)                     :: bndps

        ! Loop
        integer(I8)                             :: i, j, k 

        ! Initialize
        !===========
        ! Initialize dynamic arrays

        ! Construct boundary polygonset
        !==============================
        ! Check if we just use the existing vessel polygon set or if
        ! we retrace 
        if (options%dotpvesselbased) then
            ! Unpack
            bndps = vessel%polygonset
            bndpol = bndps%polygons
        else
            ! Boundary contours
            bndcontours = boundarytracer%TraceContours([0.0_R8])
            
            ! Sanity check
            if (size(bndcontours) == 0) then 
                ! Throw error - no boundary polygon
                call gdErrorHandler('AddTopologicalMeshTangencyPoints2: could not trace ' // & 
                    'boundary contour, check input')
            end if 

            
            ! Construct polygons from boundary
            allocate(bndpol(size(bndcontours)), keepind(size(bndcontours)))
            keepind = .true. 
            do i = 1, size(bndcontours)
                ! Construct
                call bndpol(i)%Construct(bndcontours(i)%x, bndcontours(i)%y)

                ! Check for closedness, if not -> warning and remove
                if (.not. bndpol(i)%isclosed) then 
                    print *, 'AddTopologicalMeshTangencyPoints2: ' //&
                        'boundary polygon part is not closed and will not be ' // &
                        'added'
                    print *, 'part: ', i
                    keepind(i) = .false. 
                    !call gdErrorHandler('AddTopologicalMeshTangencyPoints2: boundary polygon ' // & 
                    !    'is not closed, not supported. Check input')
                end if 
            end do 

            ! Get only closed polygons
            allocate(realbndpol(count(keepind)))
            realbndpol = pack(bndpol, keepind)
            if (size(realbndpol) == 0) then 
                call gdErrorHandler('AddTopologicalMeshTangencyPoints2: no ' // & 
                    'closed boundary polygons found, check input')
            end if 
            bndpol = realbndpol

            ! Construct polygonset
            call bndps%Construct(bndpol)

            ! Orient
            call bndps%OrientNestedClosedPolygons(flag)

            ! Check if successful
            if (flag /= 0) then
                call gdErrorHandler('AddTopologicalMeshTangencyPoints2: could not orient ' // & 
                    'boundary polygons, check input')
            end if 
        end if

        ! Compute tangency points
        !========================
        do i = 1, size(bndps%polygons)
            ! Associate for ease
            associate(p         => bndps%polygons(i))

            ! Evaluate field at vertex locations
            allocate(val(size(p%vert)))
            call magneticField%interp%Evaluate(p%x(p%vert), p%y(p%vert), &
                0, 0, val) ! assumed start and end point the same
            val = [val, val(2)] ! extend to take next edge into account

            ! Take difference
            dval = val(2:size(val)) - val(1:size(val)-1)

            ! Check where this changes sign
            allocate(dvals(size(dval)))
            where (dval > 0) dvals = 1
            where (dval <= 0) dvals = -1
            ddvals = dvals(2:size(dvals)) - dvals(1:size(dvals)-1) 

            ! Find the location and value of extrema
            extrloc = [.false., ddvals /= 0]
            allocate(tv(count(extrloc)))
            tv = pack(p%vert, extrloc)
            tx = p%x(tv)
            ty = p%y(tv)
            tf = pack(val(1:size(val)-1), extrloc)
            
            ! Check if we should exclude extremum pairs based on field value
            ! difference
            fdifftol = (maxval(val) - minval(val))*tprelfieldtol
            k = 1
            hasbeendeleted = .false.
            do while (k < size(tf))
                ! Check difference
                if ((abs(tf(k+1)-tf(k)) < fdifftol) .or. (tv(k+1)-tv(k) == 1)) then ! 
                    ! Remove values, such that subsequent ones can be
                    ! checked too
                    hasbeendeleted = .true.
                    extrloc(tv(k:k+1)) = .false.
                    tv = [tv(:k-1), tv(k+1:)]
                    tf = [tf(:k-1), tf(k+1:)]
                    k = 1
                else
                    k = k + 1
                end if
            end do 
            
            ! Check last 'edge' - shouldn't do this, since boundary
            ! polygons are closed and therefore first and last 
            ! vertex are the same
            !if (abs(tf(size(tf))-tf(1)) < fdifftol) then 
            !    extrloc(tv([1, size(tf)])) = .false.
            !    tv = pack(p%vert, extrloc)
            !    tf = val(tv)
            !    hasbeendeleted = .true.
            !end if
            
            ! Issue message
            if (hasbeendeleted) then 
                print *, 'AddTopologicalMeshTangencyPoints2: some ' // & 
                    'tangency points were deleted based on their field ' // & 
                    'values as they are very close together'
            end if

            ! Recompute points
            tv = pack(p%vert, extrloc)
            extrlocind = tv 
            extrlocind = pack([(k, k = 1, size(extrloc))], extrloc)
            tx = p%x(tv)
            ty = p%y(tv)
            tf = val(extrlocind)

            ! Compute product between normal at vertex and magnetic field
            allocate(dFdx(size(tf)), dFdy(size(tf)))
            call magneticField%interp%Evaluate(tx, ty, 1, 0, dFdx)
            call magneticField%interp%Evaluate(tx, ty, 0, 1, dFdy)

            ! Get normals of edges
            nxpe = [p%nx, p%nx(1)] 
            nype = [p%ny, p%ny(1)]

            ! Check orientation - default polygon orientation is inward pointing,
            ! we need outward pointing
            doflip = [p%edges(:, 1) == p%vert(1:p%ne)] ! if equal, then the normal points inwards and we need to flip
            doflip = [doflip, doflip(1)]
            where (doflip)
                nxpe = -nxpe
                nype = -nype 
            end where

            ! Get normals at vertices
            nxp = 0.5*(nxpe(extrlocind-1) + nxpe(extrlocind))
            nyp = 0.5*(nype(extrlocind-1) + nype(extrlocind))

            ! Compute dot product
            normprod = dFdx*nxp + dFdy*nyp

            ! Determine type and add
            allocate(tt(size(tf)))
            tt = 0
            do j = 1, size(extrlocind)

                ! Determine type
                if ((val(extrlocind(j)-1) < tf(j)) .and. &
                    (val(extrlocind(j)+1) < tf(j))) then 

                    ! Local maximum
                    if (normprod(j) < 0) then 
                        ! Curve bends outwards of the domain
                        tt(j) = TMvertextp1ID
                    else
                        ! Curve bends inwards
                        tt(j) = TMvertextp2ID
                    end if
                else
                    ! Local minimum
                    if (normprod(j) > 0) then 
                        ! Curve bends outwards of the domain
                        tt(j) = TMvertextp1ID
                    else
                        ! Curve bends inwards
                        tt(j) = TMvertextp2ID
                    end if 
                end if
            end do 

            ! Add to topological mesh
            !========================
            ! Extract faces from polygon
            eID = [(k, k = 1, size(tv))]
            sID = tv 
            ntp = size(tv)
            allocate(sortind(size(sID)))
            call Sort(sID, ind=sortind)
            eID = eID(sortind)

            ! Check if the first vertex is a tangency point. If so, add 
            ! the first vertex as well
            if (sID(1) == 1) then 
                sID = [sID, bndpol(i)%ne]
                eID = [eID, eID(1)]
            end if 
            call ExtractTopologicalFacesFromPolygon(bndpol(i), eID, sID, &
                tx, ty, vf1, vf2, xf, yf)  
            deallocate (sortind) 

            ! Update face vertices to correct indices for after adding
            ! vertices
            vf1 = vf1 + topomesh%vert%ntot 
            vf2 = vf2 + topomesh%vert%ntot 

            ! Add vertices
            do j = 1, ntp 
                call AddTopologicalMeshVertex(topomesh, tx(j), ty(j), &
                    tf(j), tt(j), topomesh%nfs+1)
                topomesh%nfs = topomesh%nfs + 1
            end do 

            ! Add faces
            do j = 1, size(xf)
                call AddTopologicalMeshFace(topomesh, [vf1(j), vf2(j)], &
                    xf(j), yf(j), TMfacebndID, 0_I8, 0.0_R8)
            end do 

            ! Housekeeping
            deallocate(val, dvals, dFdx, dFdy, tt, tv)
            end associate

        end do 
    end subroutine

    ! Necessary contours
    subroutine AddTopologicalMeshContours(topomesh, magneticField, vessel, &
        fieldtracer, options)

        ! Description
        !============
        ! This routine adds all necessary remaining contours of the magnetic field
        ! and computes their intersection(s) with the existing boundaries. It is
        ! assumed that all radial lines (cuts, vessel boundaries, etc) are already
        ! properly added to the topological mesh and that sufficient information is
        ! available for which points to compute the contours. The contours of the
        ! following points should be traced, starting from that point:
        !
        ! - tangency points with contour in the vessel (type 5)
        ! - saddle points (type 2)
        ! - split points on boundaries 
        !
        ! At this point, no contours for maxima/minima (types 1, 3) or for other
        ! tangency points (type 4) are added. This can be done in a later step for
        ! gridding purposes (i.e. if one does not want a region to be gridded near
        ! a maximum or tangency point).

        ! Algorithm
        !==========
        ! 1) Compute all additional contour parts starting from the points of
        ! interest. Make sure to parse all possible saddle point locations to the
        ! tracing routine.
        ! 2) Compute all intersections between all current boundaries and the new
        ! contours (contours themselves shouldn't intersect!)
        ! 3) Add the new vertices and faces to the topological mesh
        ! 4) Remove all points outside of the boundary and delete any associated
        ! faces, vertices, ...
        ! 5) Check for any faces that are closed or duplicate. If there are, split
        ! them up into parts such that each face has its own unique vertex pair.

        ! Notes
        !======
        ! Note 1: it is assumed that saddle points which should be treated with
        ! equal field value have exactly the same field value.

        ! Note 2: for tangency points, we need to be careful when tracing contours.
        ! The initial or end part of the contour, which should start or end at the
        ! tangency point, should lie strictly inside the vessel, i.e. the first
        ! point should be inside the vessel domain. If this is not the case, these
        ! points are eliminated up to the first point that lies inside of the
        ! vessel again. Problems may arise if contour mesh resolution is (locally)
        ! not fine enough, and mesh refinement may be needed. However, this may
        ! require a very fine mesh in cases where the contour and vessel are almost
        ! aligned over a substantial part of the domain. Therefore, we propose a
        ! different approach: first, we compute the contours. Then, we check which
        ! boundary faces (type 3) the contours intersect. We then process the
        ! contours and boundary faces as follows: only segments are kept that start
        ! from the tangency point. However, if the first of such segments goes
        ! outside of the vessel, we have to treat it differently (this can happen
        ! due to different order of approximation of the contour (first order) and
        ! determination of tangency points (higher order)). Then, we remove this
        ! first segment from tangency point to the intersection and connect the
        ! tangency point to the first next point of the contour. We do the same for
        ! the vessel polygon, which ensures a proper and expected topology.

        ! Note 3: other boundaries and points may be inserted afterwards that still
        ! lie outside of the domain (this is typically the case with saddle point
        ! contours). These contours are not trimmed beforehand since we still need
        ! to find the proper intersections with the vessel wall. 

        ! Note 4: also for saddle points/cut-like boundaries, we need
        ! to be careful when computing intersections. It is possible that
        ! near the x-point, the separatrix intersects with a core boundary
        ! due to inaccurate tracing of the field line contour locally. 
        ! This needs to be hedged for by modifying the contours/lines 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        type(MagneticFieldUDT),intent(in)       :: magneticField 
        type(VesselUDT), intent(inout)          :: vessel
        class(ContourTracerUDT), intent(inout)  :: fieldtracer 
        type(TopomeshOptionsUDT), intent(in)    :: options 

        ! Auxiliary
        real(R8)                                :: dist, startsr, endsr
        real(R8), allocatable, dimension(:)     :: pspx, pspy, pspf, &
            xout, yout, iout, jout, tscr, tsfr, tx, ty, tsrfh, tsrsh
        integer(I8)                             :: npsp, &
            intersectind, nstc, ntpc, nint, startind, endind, indtpc
        integer(I8), allocatable, dimension(:)  :: psptype, &
            sortind, pspID, allcurvetypes, allfsIDs, &
            vindI, vindJ, notdelind, vindIfh, vindIsh, tsc, tfaceind, &
            tsf, tcontourind, tcstartind, tcendind
        logical                                 :: intfacestart, &
            intcstart, dointersect
        logical, allocatable, dimension(:)      :: tracepoints, keepind, &
            ispseudosaddlepoint, isnbface, rmind, isstartingcontour, &
            isendingcontour

        type(ContourUDT), allocatable           :: tc(:), allc(:), alltpc(:)
        type(PolygonUDT), allocatable           :: allpc(:)
        type(PolygonSetUDT)                     :: tempps
        type(IntegerDynamicArrayUDT)            :: curvetypes, fsIDs
        type(IntegerDynamicArrayUDT), allocatable       :: sc(:), &
            faceind(:), sf(:), contourind(:)
        type(RealDynamicArrayUDT), allocatable          :: scr(:), sfr(:)

        type(PolygonSetUDT)                             :: bndps
        type(PolygonUDT), allocatable, dimension(:)     :: bndpol
        class(PLF2DOptionsUDT), allocatable     :: bndplfoptions

        ! Loop 
        integer(I8)                             :: i, j, k, cc
        
        ! Initialize
        !===========
        ! Associate
        associate(nfs       => topomesh%nfs)
            
        ! Initialize the contour structure & dynamic arrays
        allocate(allc(0))
        curvetypes = ConstructIntegerDynamicArray()
        fsIDs = ConstructIntegerDynamicArray()

        ! Build 'saddle' points
        ispseudosaddlepoint = (topomesh%vert%type == TMvertexmaxID) .or. &
            (topomesh%vert%type == TMvertexminID) .or. &
            (topomesh%vert%type == TMvertexsaddleID) .or. & 
            (topomesh%vert%type == TMvertextp1ID) .or. & 
            (topomesh%vert%type == TMvertextp2ID) .or. &
            (topomesh%vert%type == TMvertexsplitID)
        npsp = count(ispseudosaddlepoint)
        allocate(pspx(npsp), pspy(npsp), pspf(npsp), psptype(npsp), &
            pspID(npsp))
        pspx    = pack(topomesh%vert%x, ispseudosaddlepoint)
        pspy    = pack(topomesh%vert%y, ispseudosaddlepoint) 
        pspf    = pack(topomesh%vert%fval, ispseudosaddlepoint) 
        psptype = pack(topomesh%vert%type, ispseudosaddlepoint)
        pspID   = pack(topomesh%vert%ID, ispseudosaddlepoint)

        ! Check for saddle points that are connected - don't trace twice
        allocate(tracepoints(npsp))
        tracepoints = .true.
        !do i = 1, npsp-1
        !    do j = i+1, npsp
        !        if ((topomesh%vert%fsID(i) == topomesh%vert%fsID(j)) .and. &
        !            (psptype(i) == TMvertexsaddleID) .and. (psptype(j) == TMvertexsaddleID)) then 
        !            tracepoints(j) = .false.
        !        end if 
        !    end do 
        !end do

        ! Add saddle points to field tracer - only tangency points and
        ! extrema, no other 'regular' vertices

        fieldtracer%xs = pspx 
        fieldtracer%ys = pspy
        fieldtracer%vs = pspf
        fieldtracer%IDs = pspID
        fieldtracer%order = 0*psptype   

        ! Trace contours
        !===============
        ! Add saddle point contours
        !-------------------------- 
        do i = 1, npsp
            if ((psptype(i) == TMvertexsaddleID) .and. tracepoints(i)) then 
                ! Trace
                tc = fieldtracer%TraceContours([pspx(i)], [pspy(i)])
                
                ! Process
                call CleanContours(tc)

                ! Check if any other saddle points were encountered 
                ! during tracing. If so, do not trace these anymore 
                ! (leads to duplicate faces)
                do j = 1, size(tc)
                    if (tc(j)%startsaddle /= 0) then
                        if (psptype(tc(j)%startsaddle) == TMvertexsaddleID) then  
                            tracepoints(tc(j)%startsaddle) = .false. 
                        end if 
                    end if 
                    if (tc(j)%endsaddle /= 0) then
                        if (psptype(tc(j)%endsaddle) == TMvertexsaddleID) then  
                            tracepoints(tc(j)%endsaddle) = .false. 
                        end if 
                    end if 
                end do 

                ! Check intersections with radial core boundaries
                do j = 1, topomesh%face%ntot
                    ! Check where exactly it (should) intersect
                    if ((topomesh%face%type(j) == TMfaceradID)) then 
                        
                        ! Loop over all contours
                        do k = 1, size(tc)
                            ! Check for intersections
                            dointersect = .true.
                            if ((topomesh%face%vert(j, 1) == tc(k)%startsaddle)) then
                                intfacestart = .true.
                                intcstart = .true.
                            elseif ((topomesh%face%vert(j, 2) == tc(k)%startsaddle)) then
                                intfacestart = .false.
                                intcstart = .true.
                            elseif ((topomesh%face%vert(j, 1) == tc(k)%endsaddle)) then
                                intfacestart = .true.
                                intcstart = .false.
                            elseif ((topomesh%face%vert(j, 2) == tc(k)%endsaddle)) then
                                intfacestart = .false.
                                intcstart = .false.
                            else 
                                ! No intersection found - cycle
                                dointersect = .false. 
                            end if 

                            if (dointersect) then 
                                call SimplePolygonIntersections(tc(k)%x, tc(k)%y, &
                                    topomesh%face%x(j)%Get(), topomesh%face%y(j)%Get(), &
                                    xout, yout, vindI, vindJ, iout, jout)

                                ! Check intersections - first one should be in x-point itself
                                ! Note: we don't need to hedge for closed contours here yet
                                nstc = size(tc(k)%x)-1 ! number of segments
                                if (allocated(iout)) then 
                                    if (size(iout) < 1) then 
                                        call gdErrorHandler('AddTopologicalMeshContour: ' // & 
                                            'could not find any intersection ' // & 
                                            'between radial line and contour ' // & 
                                            'though they should start or end ' // &
                                            'in the same point')
                                    end if 
                                else 
                                    call gdErrorHandler('AddTopologicalMeshContour: ' // & 
                                            'could not find any intersection ' // & 
                                            'between radial line and contour ' // & 
                                            'though they should start or end ' // &
                                            'in the same point')
                                end if 

                                ! If any intersections are found, check
                                if (intfacestart .and. intcstart) then 
                                    if (.not. (jout(1) == 0.0_R8) .and. .not. (iout(1) == 0.0_R8)) then 
                                        call gdErrorHandler('AddTopologicalMeshContour: ' // & 
                                            'first intersection between separatrix contour ' // & 
                                            'and radial line should be in first ' // & 
                                            'point of separatrix but this is not the case')
                                    end if
                                end if 
                                if (intfacestart .and. .not. intcstart) then 
                                    if (.not. (jout(1) == 0.0_R8) .and. .not. (iout(size(iout)) == nstc)) then 
                                        call gdErrorHandler('AddTopologicalMeshContour: ' // & 
                                            'first intersection between separatrix contour ' // & 
                                            'and radial line should be in first ' // & 
                                            'point of separatrix but this is not the case')
                                    end if
                                end if 
                                if (.not. intfacestart .and. intcstart) then 
                                    if (.not. (jout(size(jout)) == topomesh%face%x(j)%Size()-1) .and. &
                                        .not. (iout(1) == 0.0_R8)) then 
                                        call gdErrorHandler('AddTopologicalMeshContour: ' // & 
                                            'first intersection between separatrix contour ' // & 
                                            'and radial line should be in first ' // & 
                                            'point of separatrix but this is not the case')
                                    end if
                                end if 
                                if (.not. intfacestart .and. .not. intcstart) then 
                                    if (.not. (jout(size(jout)) == topomesh%face%x(j)%Size()-1) .and. &
                                        .not. (iout(size(iout)) == nstc)) then 
                                        call gdErrorHandler('AddTopologicalMeshContour: ' // & 
                                            'first intersection between separatrix contour ' // & 
                                            'and radial line should be in first ' // & 
                                            'point of separatrix but this is not the case')
                                    end if
                                end if

                                ! Sort the intersections according to contour
                                sortind = [(k, k = 1, size(xout))]
                                call Sort(iout, ind=sortind, ascend=.true.)
                                vindI = vindI(sortind)
                                iout = iout(sortind)
                                
                                ! Sort the intersections according to face
                                sortind = [(k, k = 1, size(xout))]
                                call Sort(jout, ind=sortind, ascend=.true.)
                                vindJ = vindJ(sortind)
                                jout = jout(sortind)

                                ! If there are other intersections, adjust the face and the contour
                                if (size(xout) > 1) then 
                                    ! Face
                                    !-----
                                    tx = topomesh%face%x(j)%Get()
                                    ty = topomesh%face%y(j)%Get()
                                    if (intfacestart) then 
                                        ! Face is oriented from start to end
                                        call DeleteCurveSegment(tx, ty, [maxval(jout)], 'start', &
                                            [distfrac, distfrac], .true., .true.)
                                        ! delind = [(cc, cc = 2, maxval(vindJ))]
                                    else
                                        call DeleteCurveSegment(tx, ty, [jout(1)], 'end', [distfrac, distfrac], .true., .true.)
                                        ! delind = [(cc, cc = vindJ(1)+1, topomesh%face%x(j)%Size()-1)]
                                    end if 

                                    ! Delete
                                    !call topomesh%face%x(j)%Remove(delind)
                                    !call topomesh%face%y(j)%Remove(delind)
                                    call topomesh%face%x(j)%Set(tx)
                                    call topomesh%face%y(j)%Set(ty)

                                    ! Reconvert to polygon
                                    call topomesh%face%pol(j)%Construct(&
                                        topomesh%face%x(j)%Get(), topomesh%face%y(j)%Get())

                                    ! Contour
                                    !--------
                                    ! Need to hedge for closed contour
                                    if (tc(k)%isclosed) then
                                        ! Sanity check
                                        if ((iout(1) /= 0.0_R8) .or. &
                                            (iout(size(iout)) /= size(tc(k)%x)-1)) then 
                                            ! Intersection should happen in both end and start point
                                            print *, 'vertex : ', i
                                            call gdErrorHandler('AddTopologicalMeshContours: ' // & 
                                                'closed contour does not intersect in both ' // &
                                                'start and en points, check input')
                                        end if

                                        ! Check if there are exactly two
                                        ! intersections - in that case, 
                                        ! this is simply a radial line
                                        ! intersecting with a closed 
                                        ! contour
                                        if (size(iout) /= 2) then 
                                            ! If not, issue warning
                                            print *, 'vertex: ', i 
                                            print *, 'AddTopologicalMeshContours: ' // & 
                                                'multiple intersections with ' // & 
                                                'closed separatrix found, attempting ' // & 
                                                'to adjust separatrix and radial line' 
                                        end if 

                                        

                                        ! Take only intersections that 
                                        ! are not in end points
                                        vindI = vindI(2:size(vindI)-1)

                                        ! Partition into intersections in 
                                        ! first half of curve and second half
                                        allocate(vindIfh(count(vindI < nstc/2_I8)), vindIsh(count(vindI >= nstc/2_I8)))
                                        vindIfh = pack(vindI, vindI < nstc/2_I8)
                                        vindIsh = pack(vindI, vindI >= nstc/2_I8)

                                        ! Check which part to keep
                                        if ((size(vindIfh) > 0) .and. (size(vindIsh) > 0)) then 
                                            ! Intersections at both sides, take 'middle' piece
                                            notdelind = [1, (cc, cc = maxval(vindIfh)+1, minval(vindIsh)-1), nstc]
                                        elseif (size(vindIfh) > 0) then 
                                            ! Intersection only at first half
                                            notdelind = [1, (cc, cc = maxval(vindIfh)+1, nstc)]
                                        elseif (size(vindIsh) > 0) then 
                                            ! Intersection only at second half
                                            notdelind = [(cc, cc = 1, minval(vindIsh)-1), nstc]
                                        else
                                            ! Only two intersections in 
                                            ! start and end, move along
                                            notdelind = [(cc, cc = 1, size(tc(k)%x))]

                                        end if 
                                        deallocate(vindIfh, vindIsh)

                                    elseif (intcstart) then 
                                        ! Contour is oriented from start to end
                                        notdelind = [1, (cc, cc = vindI(2)+1, size(tc(k)%x))]
                                    else
                                        notdelind = [(cc, cc = 1, vindI(1)), size(tc(k)%x)]
                                    end if  

                                    ! Remove
                                    tc(k)%x = tc(k)%x(notdelind)
                                    tc(k)%y = tc(k)%y(notdelind)
                                end if 
                            end if 
                        end do 
                    end if 
                end do

                ! Add
                allc = [allc, tc]
                call curvetypes%Append(spread(TMfacesepID, 1, size(tc)))
                
                ! Add flux surface ID
                nfs = nfs + 1
                call fsIDs%Append(spread(nfs, 1, size(tc)))
            end if 
        end do

        ! Add tangency point contours
        !----------------------------
        ! First, trace all contours
        allocate(alltpc(0))
        !$omp parallel do default(none) private(i, k, tc, keepind, dist) & 
        !$omp shared(fieldtracer, pspx, pspy, psptype, npsp, tracepoints, &
        !$omp alltpc, curvetypes, fsIDs, pspid) & 
        !$omp schedule(dynamic)
        do i = 1, npsp
            if ((psptype(i) == TMvertextp2ID) .and. tracepoints(i)) then 
                ! Trace contour
                tc = fieldtracer%TraceContours([pspx(i)], [pspy(i)])

                ! Process
                call CleanContours(tc)

                ! Checks
                allocate(keepind(size(tc)))
                keepind = .true.
                do k = 1, size(tc)
                    ! Is the starting point the actual given start point? (should
                    ! be exactly the same since added as starting point in the
                    ! contouring algorithm)
                    dist = sqrt((tc(k)%x(1) - pspx(i))**2 + (tc(k)%y(1) - pspy(i))**2)
                    if (dist > 0.0_R8) then 
                        print *, 'AddTopologicalMEshContours: tangency ' // & 
                            'contour segment found that does not start ' // & 
                            'in given tangency point. Removing...'
                        keepind(k) = .false.
                    end if
                    if (tc(k)%startsaddle /= pspID(i)) then
                        ! Normally this should come from the contour 
                        ! tracer, but we can add it afterwards as well 
                        print *, 'AddTopologicalMeshContours: tangency ' // &
                            'contour segment found that starts in given ' // & 
                            'tangency point, but that does not have the ' // &
                            'starting saddle point ID as tangency point. ' // & 
                            'adjusting starting ID...'
                        tc(k)%startsaddle = pspID(i)
                    end if 
                    if (tc(k)%isclosed .and. (tc(k)%endsaddle /= pspID(i))) then 
                        ! Ensure start and end saddle point are the same
                        tc(k)%endsaddle = pspID(i)
                    end if 
                end do

                ! Remove
                tc = pack(tc, keepind)
                deallocate(keepind)

                ! Add
                !$omp critical
                alltpc = [alltpc, tc]
                call curvetypes%Append(spread(TMfacepolID, 1, size(tc)))
                
                ! Add flux surface ID
                call fsIDs%Append(spread(i, 1, size(tc)))
                !$omp end critical
            end if 
        end do 
        !$omp end parallel do 

        ! Compute all intersections with other faces
        ntpc = size(alltpc)
        allocate(sc(ntpc), scr(ntpc), faceind(ntpc), sf(topomesh%face%ntot), &
            sfr(topomesh%face%ntot), contourind(topomesh%face%ntot))
        do i = 1, ntpc
            sc(i) = ConstructIntegerDynamicArray()
            scr(i) = ConstructRealDynamicArray()
            faceind(i) = ConstructIntegerDynamicArray()
        end do 
        do i = 1, topomesh%face%ntot
            sf(i) = ConstructIntegerDynamicArray()
            sfr(i) = ConstructRealDynamicArray()
            contourind(i) = ConstructIntegerDynamicArray()
        end do 

        ! Loop
        !$omp parallel do default(none) &
        !$omp schedule(dynamic) collapse(2) & 
        !$omp private(xout, yout, vindI, vindJ, iout, jout) & 
        !$omp shared(ntpc, topomesh, alltpc, sc, scr, faceind, sf, sfr, contourind)
        do i = 1, ntpc
            do j = 1, topomesh%face%ntot
                if (topomesh%face%type(j) == TMfacebndID) then 
                    ! Compute intersections
                    call SimplePolygonIntersections(alltpc(i)%x, alltpc(i)%y, &
                        topomesh%face%x(j)%Get(), topomesh%face%y(j)%Get(), &
                        xout, yout, vindI, vindJ, iout, jout)

                    ! Add, if any
                    if (allocated(xout)) then 
                        if (size(xout) > 0) then 
                            !$omp critical
                            ! Add intersection data to contour
                            call sc(i)%Append(vindI)
                            call scr(i)%Append(iout)
                            call faceind(i)%Append(spread(j, 1, size(vindI)))

                            ! Add intersection data to face
                            call sf(j)%Append(vindJ)
                            call sfr(j)%Append(jout)
                            call contourind(j)%Append(spread(i, 1, size(vindI))) 
                            !$omp end critical
                        end if 
                    end if
                end if
            end do 
        end do
        !$omp end parallel do 

        ! Process contours
        do i = 1, ntpc
            ! Unpack
            tsc = sc(i)%Get()
            tfaceind = faceind(i)%Get()
            tscr = scr(i)%Get()
            nint = size(tsc)
            nstc = size(alltpc(i)%x)-1

            ! Sort 
            sortind = [(k, k = 1, nint)]
            call Sort(tscr, ind=sortind, ascend=.true.)
            tsc = tsc(sortind)
            tfaceind = tfaceind(sortind)

            ! Check
            if (tscr(1) /= 0.0_R8) then 
                ! First intersection should always be in tangency point
                call gdErrorHandler('AddTopologicalMeshContours: first ' // & 
                    'intersection of tangency point contour is not in  ' // & 
                    'tangency point itself' )
            end if 
            if (alltpc(i)%isclosed) then 
                ! Last intersection should also be in tangency point
                if (tscr(nint) /= size(alltpc(i)%x)-1) then 
                    call gdErrorHandler('AddTopologicalMeshContours: last ' // & 
                    'intersection of closed tangency point contour is not in  ' // & 
                    'tangency point itself' )
                end if 
            end if 

            ! Determine which intersections were in a 'neighbouring' face
            ! (i.e. a boundary face that has the tangency point)
            allocate(isnbface(nint))
            isnbface = .false.
            do j = 1, nint
                if (any(topomesh%face%vert(tfaceind(j), :) == alltpc(i)%startsaddle)) then 
                    isnbface(j) = .true.
                end if 
            end do             

            ! Determine which part(s) of contour(s) to keep
            if (alltpc(i)%isclosed) then ! closed contour
                ! This is much more tricky, need to check additional cases
                if (all(isnbface)) then 
                    ! No intersections with other boundaries
                    if (nint == 2) then 
                        ! Simply in start and end - add full contour, so 
                        ! do nothing
                    else
                        ! Here, we need to look at which face IDs the 
                        ! intersections have, since it is 
                        ! Here, we don't have much more to go on than 
                        ! assuming that these intersections are somewhere
                        ! in the start/end nodes of the contour and that
                        ! we can partition based on that. This is of course
                        ! not a general way and may fail if the amount of
                        ! contour points is too low
                        ! Print a warning
                        print *, 'AddTopologicalMeshContours: multiple ' // & 
                            'intersections found for closed contour with ' // & 
                            'only neighbouring boundary faces, attempting to ' // &
                            'split contour by partitioning intersection into ' // & 
                            'intersections at start and end. This may not result ' // & 
                            'in desired behavior... (vertex: ', alltpc(i)%startsaddle, ')'

                        ! Determine start & end intersections
                        allocate(vindIfh(count(tsc < nstc/2)), &
                            vindIsh(count(tsc >= nstc/2)))
                        allocate(tsrfh(size(vindIfh)), tsrsh(size(vindIsh)))
                        vindIfh = pack(tsc, tsc < nstc/2)
                        vindIsh = pack(tsc, tsc >= nstc/2)
                        tsrfh = pack(tscr, tsc < nstc/2)
                        tsrsh = pack(tscr, tsc >= nstc/2)

                        ! Determine start index
                        if (size(vindIfh) > 0) then 
                            startind = maxval(vindIfh)
                            startsr = maxval(tsrfh)
                        else
                            startind = 2
                            startsr = 0.0_R8
                        end if 

                        ! Determine end index
                        if (size(vindIsh) > 0) then 
                            endind = minval(vindIsh)+1
                            endsr = minval(tsrsh)
                        else
                            endind = nstc
                            endsr = real(nstc, kind=R8)
                        end if 
                        
                        ! Adjust contour
                        if (startsr /= 0.0_R8 .and. endsr /= real(nstc, kind=R8)) then 
                            call DeleteCurveSegment(alltpc(i)%x, alltpc(i)%y, &
                                [startsr, endsr], 'both', [distfrac, distfrac], .true., .true.)
                        elseif (startsr /= 0.0_R8) then 
                            call DeleteCurveSegment(alltpc(i)%x, alltpc(i)%y, &
                                [startsr], 'start', [distfrac, distfrac], .true., .true.)
                        elseif (endsr /= real(nstc, kind=R8)) then 
                            call DeleteCurveSegment(alltpc(i)%x, alltpc(i)%y, &
                                [endsr], 'end', [distfrac, distfrac], .true., .true.)
                        end if 
                        !alltpc(i)%x = alltpc(i)%x([1, (k, k = startind, endind), nstc+1])
                        !alltpc(i)%y = alltpc(i)%y([1, (k, k = startind, endind), nstc+1])

                        ! Housekeeping
                        deallocate(vindIfh, vindIsh, tsrfh, tsrsh)
                    end if
                else
                    ! At least one intersection with another boundary. 
                    ! Note: here we do want to keep the intersection 
                    ! with these other boundaries in the contour!
                    ! Need to find first and second segment 

                    ! First segment
                    endind = findloc(isnbface, .false., 1, back=.false.)

                    ! Add this segment as additional contour
                    alltpc = [alltpc, alltpc(i)]
                    indtpc = size(alltpc)
                    alltpc(indtpc)%isclosed = .false.
                    alltpc(indtpc)%endsaddle = 0 ! doesn't end anymore in saddle point
                    if (tscr(endind-1) == 0.0_R8) then 
                        call DeleteCurveSegment(alltpc(indtpc)%x, alltpc(indtpc)%y, &
                            [tscr(endind)], 'end', [-distfrac], .true., .false.)
                    else
                        ! Also need to delete a first part
                        call DeleteCurveSegment(alltpc(indtpc)%x, alltpc(indtpc)%y, &
                            [tscr(endind-1:endind)], 'both', [distfrac, -distfrac], .true., .false.)
                    end if 
                    !call DeleteCurveSegment(alltpc(indtpc)%x, alltpc(indtpc)%y, &
                    !    [tscr(endind)], 'end', [0.0_R8], .true., .false.)
                    !alltpc(indtpc)%x = alltpc(i)%x([1, (k, k = tsc(startind-1)+1, tsc(startind)+1)])
                    !alltpc(indtpc)%y = alltpc(i)%y([1, (k, k = tsc(startind-1)+1, tsc(startind)+1)])

                    ! Append flux surface ID etc as well!
                    call curvetypes%Append(curvetypes%Get(size(allc) + i))
                    call fsIDs%Append(fsIDs%Get(size(allc) + i))

                    ! Second segment
                    startind = findloc(isnbface, .false., 1, back=.true.)

                    ! Add this segment by adjusting existing contour
                    alltpc(i)%isclosed = .false.
                    alltpc(i)%startsaddle = 0 ! doesn't start anymore in saddle point
                    if (tscr(startind+1) == real(nstc, kind=R8)) then 
                        call DeleteCurveSegment(alltpc(i)%x, alltpc(i)%y, &
                            [tscr(startind)], 'start', [-distfrac], .false., .true.)
                    else
                        ! Also need to delete last part
                        call DeleteCurveSegment(alltpc(i)%x, alltpc(i)%y, &
                            [tscr(startind:startind+1)], 'both', [-distfrac, distfrac], .false., .true.)
                    end if
                    !alltpc(i)%x = alltpc(i)%x([(k, k = tsc(endind), tsc(endind+1)), nstc+1])
                    !alltpc(i)%y = alltpc(i)%y([(k, k = tsc(endind), tsc(endind+1)), nstc+1])
                    
                end if 

            else ! open contour
                ! Check which intersection is the last intersection with
                ! the neighbouring boundary (should be first one)
                intersectind = findloc(isnbface, .true., 1, back=.true.)

                ! Sanity checks
                if (intersectind == 0) then 
                    call gdErrorHandler('AddTopologicalMeshContours: ' // & 
                        'tangency point contour intersects in tangency point ' // & 
                        'but not in neighbouring face - this is likely a bug')
                end if 
                if (any(.not. isnbface(1:intersectind))) then 
                    ! Issue warning - we got intersections of the contour
                    ! with non-neighbour faces inbetween - this is unexpected
                    print *, 'AddTopologicalMeshContours: intersections ' // & 
                        'found for tangency point ', alltpc(i)%startsaddle, &
                        'that occur inbetween intersections with contour and ' // & 
                        'neighbouring boundary face. Unexpected, intersections ' // & 
                        'are removed'
                end if 
                if (intersectind == size(isnbface)) then 
                    ! Unexpected, only intersections in neighbouring faces 
                    ! but not in other boundary face
                    call gdErrorHandler('AddTopologicalMeshContours: ' // & 
                        'tangency point contour intersects in tangency point ' // & 
                        'and neighbouring face but not in other boundary ' // &
                        'faces - unexpected')
                end if 

                ! Keep only this part of the contour coordinates
                call DeleteCurveSegment(alltpc(i)%x, alltpc(i)%y, &
                    [tscr(intersectind:intersectind+1)], 'both', [distfrac, -distfrac], .true., .false.)
                !notdelind = [(k, k = 2, tsc(intersectind))]
                !alltpc(i)%x = alltpc(i)%x(notdelind)
                !alltpc(i)%y = alltpc(i)%y(notdelind)
            end if 

            ! Housekeeping
            deallocate(isnbface)

        end do 

        ! Process faces
        do i = 1, topomesh%face%ntot 
            if (topomesh%face%type(i) == TMfacebndID) then 
                ! Unpack
                tsf = sf(i)%Get()
                tsfr = sfr(i)%Get()
                tcontourind = contourind(i)%Get()
                nint = size(tsf)

                ! Sort
                sortind = [(k, k = 1, nint)]
                call Sort(tsfr, ind=sortind, ascend=.true.)
                tsf = tsf(sortind)
                tcontourind = tcontourind(sortind)

                ! Check intersections
                if (nint > 0) then 
                    ! Check which intersections should be removed
                    allocate(rmind(nint))
                    rmind = .false.

                    do j = 1, nint
                        if (&
                            any(topomesh%face%vert(i, :) == alltpc(tcontourind(j))%startsaddle) .or. &
                            any(topomesh%face%vert(i, :) == alltpc(tcontourind(j))%endsaddle)) then 
                            rmind(j) = .true.
                        end if
                    end do 

                    ! Check how to remove intersections 
                    if (all(topomesh%vert%type(topomesh%face%vert(i, : )) == TMvertextp2ID)) then 
                        !print *, 'face ID: ', i
                        !print *, 'face vertices: ', topomesh%face%vert(i, :) 
                        !call gdErrorHandler('AddTopologicalMeshContours: ' // & 
                        !    'boundary detected with two type 2 tangency points, ' // & 
                        !    'this is likely a bug')

                        ! Sanity checks
                        if ((tsfr(1) /= 0.0_R8) .or. (tsfr(nint) /= topomesh%face%x(i)%Size()-1)) then 
                            print *, 'face ID: ', i
                            print *, 'face vertices: ', topomesh%face%vert(i, :) 
                            call gdErrorHandler('AddTopologicalMeshContours: ' //& 
                                'face should intersect with tangency point ' // & 
                                'at start and end, but not all intersections found')
                        end if 

                        ! Check which contours have the start and end
                        ! vertex as a starting/ending point
                        isstartingcontour = tsfr == 0.0_R8
                        isendingcontour = tsfr == topomesh%face%x(i)%Size()-1
                        allocate(tcstartind(count(isstartingcontour)), tcendind(count(isendingcontour)))
                        tcstartind = pack(tcontourind, isstartingcontour)
                        tcendind = pack(tcontourind, isendingcontour)

                        ! Check which intersections these are, these
                        ! give limits for face trimming
                        startind = 0
                        endind = nint+1 
                        k = 1
                        do while (k <= nint)
                            if (.not. any(tcontourind(k) == tcstartind)) then 
                                ! Exit the loop
                                exit 
                            else
                                ! Update counters
                                k = k + 1
                                startind = startind + 1
                            end if 
                        end do 
                        k = nint 
                        do while (k > 0)
                            if (.not. any(tcontourind(k) == tcendind)) then 
                                ! Exit the loop
                                exit
                            else
                                ! Update counters
                                k = k - 1
                                endind = endind - 1
                            end if 
                        end do 

                        ! Trim face
                        tx = topomesh%face%x(i)%Get()
                        ty = topomesh%face%y(i)%Get()
                        if (tsfr(startind) /= 0.0_R8 .and. tsfr(endind) /= real(topomesh%face%x(i)%Size()-1)) then 
                            call DeleteCurveSegment(tx, ty, tsfr([startind, endind]), &
                                'both', [distfrac, distfrac], .true., .true.)
                        elseif (tsfr(startind) /= 0.0_R8) then 
                            call DeleteCurveSegment(tx, ty, tsfr([startind]), &
                                'start', [distfrac, distfrac], .true., .true.)
                        elseif (tsfr(endind) /= real(topomesh%face%x(i)%Size()-1, kind=R8)) then 
                            call DeleteCurveSegment(tx, ty, tsfr([endind]), &
                                'end', [distfrac, distfrac], .true., .true.)
                        end if 

                            ! Reconstruct polygon
                        call topomesh%face%pol(i)%Construct(&
                            topomesh%face%x(i)%Get(), topomesh%face%y(i)%Get())
                        call topomesh%face%x(i)%Set(tx)
                        call topomesh%face%y(i)%Set(ty)
                        !call topomesh%face%x(i)%Remove([(k, k = 2, tsf(startind)+1), &
                        !    (k, k = tsf(endind+1), topomesh%face%x(i)%Size()-1)])
                        !call topomesh%face%y(i)%Remove([(k, k = 2, tsf(startind)+1), &
                        !    (k, k = tsf(endind+1), topomesh%face%y(i)%Size()-1)])

                        ! Reconstruct polygon
                        call topomesh%face%pol(i)%Construct(&
                            topomesh%face%x(i)%Get(), topomesh%face%y(i)%Get())

                        ! Housekeeping
                        deallocate(tcstartind, tcendind)

                    elseif (topomesh%vert%type(topomesh%face%vert(i, 1)) == TMvertextp2ID) then

                        ! Sanity checks
                        if (tsfr(1) /= 0.0_R8) then
                            print *, 'face ID: ', i
                            print *, 'face vertices: ', topomesh%face%vert(i, :) 
                            call gdErrorHandler('AddTopologicalMeshContours: ' //& 
                                'face should intersect with tangency point ' // & 
                                'initially, but no intersection found')
                        end if 

                        ! Find location
                        startind = findloc(rmind, .true., 1, back=.true.)

                        ! Sanity check
                        if (startind == 0) then 
                            ! This shouldn't happen
                            call gdErrorHandler('AddTopologicalMeshContours: ' // & 
                                'face has intersections but none to be removed - unexpected')
                        end if 

                        ! Trim face
                        if (tsfr(startind) /= 0.0_R8) then 
                            tx = topomesh%face%x(i)%Get()
                            ty = topomesh%face%y(i)%Get()
                            call DeleteCurveSegment(tx, ty, tsfr([startind]), &
                                'start', [distfrac, distfrac], .true., .true.)
                            call topomesh%face%x(i)%Set(tx)
                            call topomesh%face%y(i)%Set(ty)

                            ! Reconstruct polygon
                            call topomesh%face%pol(i)%Construct(&
                                topomesh%face%x(i)%Get(), topomesh%face%y(i)%Get())
                        end if 

                        !call topomesh%face%x(i)%Remove([(k, k = 2, tsf(startind)+1)])
                        !call topomesh%face%y(i)%Remove([(k, k = 2, tsf(startind)+1)])

                        

                    elseif (topomesh%vert%type(topomesh%face%vert(i, 2)) == TMvertextp2ID) then

                        ! Sanity checks
                        if (tsfr(nint) /= topomesh%face%x(i)%Size()-1) then
                            print *, 'face ID: ', i
                            print *, 'face vertices: ', topomesh%face%vert(i, :) 
                            call gdErrorHandler('AddTopologicalMeshContours: ' //& 
                                'face should intersect with tangency point ' // & 
                                'initially, but no intersection found')
                        end if 

                        ! Find location
                        startind = findloc(rmind, .true., 1, back=.false.)

                        ! Sanity check
                        if (startind == 0) then 
                            ! This shouldn't happen
                            call gdErrorHandler('AddTopologicalMeshContours: ' // & 
                                'face has intersections but none to be removed - unexpected')
                        end if 

                        ! Trim face
                        if (tsfr(startind) < real(topomesh%face%x(i)%Size()-1, kind=R8)) then 
                            tx = topomesh%face%x(i)%Get()
                            ty = topomesh%face%y(i)%Get()
                            call DeleteCurveSegment(tx, ty, tsfr([startind]), &
                                'end', [distfrac, distfrac], .true., .true.)
                            call topomesh%face%x(i)%Set(tx)
                            call topomesh%face%y(i)%Set(ty)

                            ! Reconstruct polygon
                            call topomesh%face%pol(i)%Construct(&
                                topomesh%face%x(i)%Get(), topomesh%face%y(i)%Get())
                        end if 
                        !call topomesh%face%x(i)%Remove([(k, k = tsf(startind), topomesh%face%x(i)%Size()-1)])
                        !call topomesh%face%y(i)%Remove([(k, k = tsf(startind), topomesh%face%y(i)%Size()-1)])

                        
                    else
                        ! This can happen, just skip since we shouldn't 
                        ! adjust anything here
                    end if

                    ! Housekeeping
                    deallocate(rmind)

                end if 

                
            end if 
        end do 

        ! Add contours (IDs etc should be done already...)
        allc = [allc, alltpc]

        


        ! Add to the topology mesh
        !=========================
        ! Clean the contours (just to be sure)
        call CleanContours(allc)

        allocate(allpc(size(allc)))
        do i = 1, size(allc)
            call allpc(i)%Construct(allc(i)%x, allc(i)%y)
        end do 
        call tempps%Construct(allpc)
        call tempps%WriteData('topocontours_all_beforeinsertion')

        ! Add the contours
        allcurvetypes = curvetypes%Get()
        allfsIDs = fsIDs%Get()
        do i = 1, size(allc)
            call InsertTopologicalMeshContour(topomesh, magneticField, &
                allc(i), allcurvetypes(i), allfsIDs(i))
        end do 

        ! Update the vessel description again
        allocate(bndpol(count(topomesh%face%type == TMfacebndID)))
        bndpol = pack(topomesh%face%pol, topomesh%face%type == TMfacebndID)
        call bndps%Construct(bndpol)
        call ConstructVesselPolygonSet(vessel, bndps)
        allocate(PLF2DClosedExactOptionsUDT::bndplfoptions)
        call InitializePolygonLevelsetFunction2D(vessel%plfvessel, vessel%polygonset, bndplfoptions)
        deallocate(bndpol)

        ! Trim the topological mesh
        call WriteTopologicalMesh(topomesh, 'topomesh_temp')
        call TrimTopologicalMesh(topomesh, magneticField, vessel)

        ! Split boundaries
        call WriteTopologicalMesh(topomesh, 'topomesh_temp')
        call SplitTopologicalMeshFaces(topomesh)   
        
        ! Boundary split vertex contours
        !===============================
        call AddBoundarySplitVertexContours(topomesh, magneticField, vessel, &
            fieldtracer)

        call SimplifyTopologicalMeshFaces(topomesh)
        
        ! Housekeeping
        !=============
        end associate

    end subroutine 

    ! Boundary split vertex contours
    subroutine AddBoundarySplitVertexContours(topomesh, magneticField, vessel, &
        fieldtracer)

        ! Description
        !============
        ! This routine adds contours that start from a split vertex 
        ! located on a boundary. These are necessary to deal with 
        ! boundaries that have been split and would result in 
        ! non-conforming flux tubes (perhaps expanding the flux tube
        ! concept will be necessary at some point...)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        type(MagneticFieldUDT),intent(in)       :: magneticField 
        type(VesselUDT), intent(inout)          :: vessel
        class(ContourTracerUDT), intent(inout)  :: fieldtracer 

        ! Auxiliary
        real(R8)                                :: dist
        integer(I8), allocatable, dimension(:)  :: allcurvetypes, allfsIDs
        logical, allocatable, dimension(:)      :: keepind

        type(ContourUDT), allocatable           :: tc(:), allc(:), &
            allbsvc(:)
        type(PolygonUDT), allocatable           :: allpc(:)
        type(PolygonSetUDT)                     :: tempps
        type(IntegerDynamicArrayUDT)            :: curvetypes, fsIDs


        ! Loop 
        integer(I8)                             :: i, k
        
        ! Initialize
        !===========
        ! Associate
        associate(nfs       => topomesh%nfs)

        ! Initialize the contour structure & dynamic arrays
        allocate(allc(0))
        curvetypes = ConstructIntegerDynamicArray()
        fsIDs = ConstructIntegerDynamicArray()

        ! First, trace all contours
        allocate(allbsvc(0))
        !$omp parallel do default(none) schedule(dynamic) &
        !$omp private(i, k, tc, keepind, dist) & 
        !$omp shared(topomesh, fieldtracer, allbsvc, curvetypes, &
        !$omp fsIDs) 
        do i = 1, topomesh%vert%ntot
            if ((topomesh%vert%type(i) == TMvertexsplitID) .and. &
                topomesh%vert%fsID(i) == 0) then 
                ! Trace contour
                tc = fieldtracer%TraceContours([topomesh%vert%x(i)], [topomesh%vert%y(i)])

                ! Process
                call CleanContours(tc)

                ! Checks
                allocate(keepind(size(tc)))
                keepind = .true.
                do k = 1, size(tc)
                    ! Is the starting point the actual given start point? (should
                    ! be exactly the same since added as starting point in the
                    ! contouring algorithm)
                    dist = sqrt((tc(k)%x(1) - topomesh%vert%x(i))**2 + &
                        (tc(k)%y(1) - topomesh%vert%y(i))**2)
                    if (dist > 0.0_R8) then 
                        print *, 'AddTopologicalMEshContours: split vertex ' // & 
                            'contour segment found that does not start ' // & 
                            'in given tangency point. Removing...'
                        keepind(k) = .false.
                    end if
                    if (tc(k)%startsaddle /= i) then
                        ! Normally this should come from the contour 
                        ! tracer, but we can add it afterwards as well 
                        print *, 'AddTopologicalMeshContours: split vertex ' // &
                            'contour segment found that starts in given ' // & 
                            'point, but that does not have the ' // &
                            'starting saddle point ID as tangency point. ' // & 
                            'adjusting starting ID...'
                        tc(k)%startsaddle = i
                    end if 
                    if (tc(k)%isclosed .and. (tc(k)%endsaddle /= i)) then 
                        ! Ensure start and end saddle point are the same
                        tc(k)%endsaddle = i
                    end if 
                end do

                ! Remove
                tc = pack(tc, keepind)
                deallocate(keepind)

                ! Add
                !$omp critical
                allbsvc = [allbsvc, tc]
                call curvetypes%Append(spread(TMfacepolID, 1, size(tc)))
                
                ! Add flux surface ID
                topomesh%nfs = topomesh%nfs + 1
                topomesh%vert%fsID(i) = topomesh%nfs
                topomesh%vert%type(i) = TMvertexbndID
                topomesh%vert%fval(i:i) = fieldtracer%Evaluate([topomesh%vert%x(i)], [topomesh%vert%y(i)])
                call fsIDs%Append(spread(topomesh%nfs, 1, size(tc)))
                !$omp end critical
            end if 
        end do 
        !$omp end parallel do 

        ! Add contours (IDs etc should be done already...)
        allc = [allc, allbsvc]

        ! Add to the topology mesh
        !=========================
        ! Clean the contours (just to be sure)
        if (size(allc) > 0) then 
            call CleanContours(allc)

            allocate(allpc(size(allc)))
            do i = 1, size(allc)
                call allpc(i)%Construct(allc(i)%x, allc(i)%y)
            end do 
            call tempps%Construct(allpc)
            call tempps%WriteData('topocontours_all_beforeinsertion2')

            ! Add the contours
            allcurvetypes = curvetypes%Get()
            allfsIDs = fsIDs%Get()
            do i = 1, size(allc)
                call InsertTopologicalMeshContour(topomesh, magneticField, &
                    allc(i), allcurvetypes(i), allfsIDs(i))
            end do 

            ! Trim the topological mesh
            call TrimTopologicalMesh(topomesh, magneticField, vessel)

            ! Split boundaries
            call SplitTopologicalMeshFaces(topomesh)   
        end if
        
        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Core boundary contours
    subroutine AddTopologicalMeshCoreBoundaries(topomesh, magneticField, &
        vessel, fieldtracer, options)

        ! Description
        !============
        ! This routine adds additional poloidal contours (type 2) for core parts.
        ! True limiter configurations are not yet supported (a warning is thrown). 
        ! The field value of this contour is determined as frac*(psiO - psiX) +
        ! psiX, where psiX is the value of the field at the X-point (or other 
        ! point), and psiO the value at the extremum (i.e. frac = 0 -> x-point,
        ! frac = 1 -> extremum). Note that values of frac close to 0 or 1 may lead
        ! to problems further downstream, as these contours may be very coarse
        ! depending on the tracing grid size. 

        ! It is assumed that all other necessary contours have been added (i.e.
        ! this routine should be invoked after 'AddTopologicalMeshContours).
        ! Otherwise, the introduced core boundary may not be properly determined.
        ! The remainder of this routine is much alike AddTopologicalMeshContours.
        
        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        type(magneticFieldUDT), intent(in)      :: magneticField
        type(VesselUDT), intent(inout)          :: vessel 
        class(ContourTracerUDT), intent(in)     :: fieldtracer 
        type(TopomeshOptionsUDT), intent(in)    :: options

        ! Auxiliary
        integer(I8)                             :: thisf
        integer(I8), allocatable, dimension(:)  :: tf, tfv, s1, s2, &
            contourtypes, fsIDs
        real(R8)                                :: thisfval, traceval 
        real(R8), allocatable, dimension(:)     :: xint, yint 
        logical, allocatable, dimension(:)      :: keepind 
        type(ContourUDT), allocatable           :: tc(:), allc(:)
        type(PolygonUDT)                        :: tcp 
        ! real(R8), allocatable, dimension(:)     :: 

        ! Loop
        integer(I8)                             :: i, j, k

        ! Initialize
        !============
        ! Associate
        associate(&
            nfs             => topomesh%nfs,    &
            vert            => topomesh%vert    &
        )

        ! Allocate
        allocate(allc(0), contourtypes(0), fsIDs(0))

        ! Trace contours
        !===============
        ! Loop
        do i = 1, vert%ntot 
            ! Check if we should trace for this extremum
            if (.not. ((vert%type(i) == TMvertexminID) .or. (vert%type(i) == TMvertexmaxID))) then 
                cycle 
            end if 

            ! Get the faces that have this extremum
            tf = pack([(k, k = 1, topomesh%face%ntot)], any(topomesh%face%vert == i, dim=2))
            !tf = findloc(any(topomesh%face%vert == i, dim=2), .true.)

            ! Sanity check
            if (size(tf) == 0) then 
                ! No faces found - limiter-like configuration, not yet 
                ! supported here
                print *, 'AddTopologicalMeshCoreBoundaries: limiter-like case ' // & 
                    'detected which is not yet supported. Not adding core ' // & 
                    'boundaries for extremum with coordinates (', & 
                    vert%x(i), ' , ', vert%y(i), ')'
                cycle 
            end if 

            ! Choose face with closest field value
            if (allocated(tfv)) then 
                deallocate(tfv)
            end if 
            allocate(tfv(size(tf)))
            do j = 1, size(tf)
                if (topomesh%face%vert(tf(j), 1) == i) then 
                    tfv(j) = topomesh%face%vert(tf(j), 2)
                else
                    tfv(j) = topomesh%face%vert(tf(j), 1) 
                end if 
            end do 
            thisf = tf(minloc(abs(vert%fval(i) - vert%fval(tfv)), dim=1))
            thisfval = vert%fval(tfv(minloc(abs(vert%fval(i) - vert%fval(tfv)), dim=1)))

            ! Compute the field value to trace
            traceval = options%coreboundariesfrac*(vert%fval(i) - thisfval) + thisfval

            ! Trace contours 
            tc = fieldtracer%TraceContours([traceval])
            call CleanContours(tc)

            ! Check which contour is closed and intersects with the 
            ! current face
            allocate(keepind(size(tc)))
            keepind = .true. 
            do j = 1, size(tc)
                ! Convert to polygon
                call tcp%Construct(tc(j)%x, tc(j)%y)

                ! Compute intersections
                call PolygonIntersections(topomesh%face%pol(thisf), tcp, &
                    xint, yint, s1, s2)

                ! Check
                if (size(xint) == 0) then 
                    keepind(j) = .false. 
                end if 
            end do

            ! Remove contours without intersections
            tc = pack(tc, keepind)
            deallocate(keepind)

            ! Check
            if (size(tc) == 0) then 
                ! No boundaries found
                print *, 'AddTopologicalMeshCoreBoundaries: no ' // & 
                    'contours found that intersect with the extremum face. ' // & 
                    'Not adding any core boundaries for extremum with ' // & 
                    'coordinates: (', vert%x(i), ', ', vert%y(i), ')'
                cycle 
            elseif (size(tc) > 1) then 
                ! Multiple boundaries found - also not expected, but
                ! may not be a problem 
                print *, 'AddTopologicalMeshCoreBoundaries: multiple' // & 
                    'contours found that intersect with the extremum face. ' // & 
                    'May lead to too many core boundaries for extremum with coordinates (', & 
                    vert%x(i), ' , ', vert%y(i), ')'
            end if 

            ! Concatenate
            allc = [allc, tc]
            contourtypes = [contourtypes, spread(TMfacecoreID, 1, size(tc))]

            ! Add flux surface ID
            nfs = nfs + 1
            fsIDs = [fsIDs, spread(nfs, 1, size(tc))]

        end do 

        ! End association
        end associate

        ! Add to the topology mesh
        !=========================
        ! Add contours
        do  i = 1, size(allc)
            call InsertTopologicalMeshContour(topomesh, magneticField, &
                allc(i), contourtypes(i), fsIDs(i))
        end do 

        ! Trim the topological mesh
        call TrimTopologicalMesh(topomesh, magneticField, vessel)

        ! Split boundaries
        call SplitTopologicalMeshFaces(topomesh)

    end subroutine 

    ! PF boundary contours
    subroutine AddTopologicalMeshPFBoundaries(topomesh, magneticField, &
        vessel, fieldtracer, options)

        ! Description
        !============
        ! This routine adds 'private flux' boundaries where applicable
        ! (here, the PF is defined as a region that has as a tangency 
        ! point with contour lying outside of the vessel, i.e. type 1 
        ! tangency points, and of which at least one of the two 
        ! neighbouring faces has a vertex that belongs to a separatrix
        ! segment). 
         
        ! The field value of this contour is determined as frac*(psiTP - psiX) +
        ! psiX, where psiX is the value of the field at the X-point (or other 
        ! point), and psiO the value at the tangency point (i.e. frac = 0 -> x-point,
        ! frac = 1 -> tp). Note that values of frac close to 0 or 1 may lead
        ! to problems further downstream, as these contours may be very coarse
        ! depending on the tracing grid size. 

        ! It is assumed that all other necessary contours have been added (i.e.
        ! this routine should be invoked after 'AddTopologicalMeshContours).
        ! Otherwise, the introduced boundary may not be properly determined.
        ! The remainder of this routine is much alike AddTopologicalMeshContours.
        
        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        type(magneticFieldUDT), intent(in)      :: magneticField
        type(VesselUDT), intent(inout)          :: vessel 
        class(ContourTracerUDT), intent(in)     :: fieldtracer 
        type(TopomeshOptionsUDT), intent(in)    :: options

        ! Auxiliary
        integer(I8)                             :: thisf
        integer(I8), allocatable, dimension(:)  :: tf, tfv, s1, s2, &
            contourtypes, fsIDs, sepfsIDs
        real(R8)                                :: thisfval, traceval
        real(R8), allocatable, dimension(:)     :: xint, yint
        logical                                 :: skipvertex
        logical, allocatable, dimension(:)      :: keepind, hasfv
        type(ContourUDT), allocatable           :: tc(:), allc(:)
        type(PolygonUDT)                        :: tcp 
        ! real(R8), allocatable, dimension(:)     :: 

        ! Loop
        integer(I8)                             :: i, j, k

        ! Initialize
        !============
        ! Associate
        associate(&
            nfs             => topomesh%nfs,    &
            vert            => topomesh%vert    &
        )

        ! Allocate
        allocate(allc(0), contourtypes(0), fsIDs(0))

        ! Get all separatrix flux surface IDs
        sepfsIDs = topomesh%GetSeparatrixFluxSurfaceIDs()

        ! Trace contours
        !===============
        ! Loop
        do i = 1, vert%ntot 
            ! Check if we should trace for this point
            if (.not. ((vert%type(i) == TMvertextp1ID))) then 
                cycle 
            end if 

            ! Get the faces that have this extremum
            hasfv = any(topomesh%face%vert == i, dim=2)
            if (allocated(tf)) then 
                deallocate(tf)
            end if 
            allocate(tf(count(hasfv)))
            tf = pack([(k, k = 1, topomesh%face%ntot)], hasfv)

            ! Sanity check
            if (size(tf) /= 2) then 
                ! Weird - normally, each tangency point of this type
                ! should only have two (radial) faces - ignore
                print *, 'AddTopologicalMeshPFBoundaries: tangeny point ' // & 
                    'vertex number: ', i, 'does not have exactly two ' // & 
                    'faces - unexpected. Skipping this tangency point...'
                cycle 
            end if 

            ! Check if both faces are target faces
            if (allocated(tfv)) then 
                deallocate(tfv)
            end if 
            allocate(tfv(size(tf)))
            do j = 1, size(tf)
                if (topomesh%face%vert(tf(j), 1) == i) then 
                    tfv(j) = topomesh%face%vert(tf(j), 2)
                else
                    tfv(j) = topomesh%face%vert(tf(j), 1) 
                end if 
            end do 

            skipvertex = .false.
            do j = 1, size(tfv)
                if (.not. any(topomesh%vert%fsID(tfv(j)) == sepfsIDs)) then 
                    skipvertex = .true. 
                    exit
                end if 
            end do 
            if (skipvertex) then 
                cycle
            end if 

            ! If we got here, check which value comes closest to tp value
            thisf = tf(minloc(abs(vert%fval(i) - vert%fval(tfv)), dim=1))
            thisfval = vert%fval(tfv(minloc(abs(vert%fval(i) - vert%fval(tfv)), dim=1)))

            ! Compute the field value to trace
            traceval = options%PFboundariesfrac*(vert%fval(i) - thisfval) + thisfval

            ! Trace contours 
            tc = fieldtracer%TraceContours([traceval])
            call CleanContours(tc)

            ! Check which contour intersects with both boundaries
            allocate(keepind(size(tc)))
            keepind = .true. 
            do j = 1, size(tc)
                ! Convert to polygon
                call tcp%Construct(tc(j)%x, tc(j)%y)

                ! Compute intersections with faces
                do k = 1, size(tf)
                    call PolygonIntersections(topomesh%face%pol(tf(1)), tcp, &
                        xint, yint, s1, s2)

                    ! Check
                    if (size(xint) == 0) then 
                        keepind(j) = .false. 
                    end if 
                end do 
            end do

            ! Remove contours without intersections
            tc = pack(tc, keepind)
            deallocate(keepind)

            ! Check
            if (size(tc) == 0) then 
                ! No boundaries found
                print *, 'AddTopologicalMeshPFBoundaries: no ' // & 
                    'contours found that intersect with the extremum face. ' // & 
                    'Not adding any core boundaries for extremum with ' // & 
                    'coordinates: (', vert%x(i), ', ', vert%y(i), ')'
                cycle 
            elseif (size(tc) > 1) then 
                ! Multiple boundaries found - also not expected, but
                ! may not be a problem 
                print *, 'AddTopologicalMeshPFBoundaries: multiple' // & 
                    'contours found that intersect with the extremum face. ' // & 
                    'May lead to too many core boundaries for extremum with coordinates (', & 
                    vert%x(i), ' , ', vert%y(i), ')'
            end if 

            ! Concatenate
            allc = [allc, tc]
            contourtypes = [contourtypes, spread(TMfacePFID, 1, size(tc))]

            ! Add flux surface ID
            nfs = nfs + 1
            fsIDs = [fsIDs, spread(nfs, 1, size(tc))]

            ! Housekeeping
            deallocate(tf, tfv)

        end do 

        ! End association
        end associate

        ! Add to the topology mesh
        !=========================
        ! Add contours
        do  i = 1, size(allc)
            call InsertTopologicalMeshContour(topomesh, magneticField, &
                allc(i), contourtypes(i), fsIDs(i))
        end do 

        ! Trim the topological mesh
        call TrimTopologicalMesh(topomesh, magneticField, vessel)

        ! Split boundaries
        call SplitTopologicalMeshFaces(topomesh)

    end subroutine 

#ifdef debug
    ! (X-point) streamlines
    subroutine AddTopologicalMeshStreamlines(topomesh, magneticField, &
        streamlinetracer, options)

        ! Description
        !============
        ! This routine adds any desired streamlines to the topological
        ! mesh. Currently, only streamlines originating from X-points 
        ! are implemented, as these are the only ones relevant for now.
        ! For X-points, we only trace the streamlines that do not go 
        ! to the core (these should already have been added in the
        ! extremum tracing routine). Furthermore, these streamlines are
        ! cut off at the first intersection with a poloidal face, or 
        ! they are not added if they first intersect with a non-poloidal
        ! face.
        
        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        type(magneticFieldUDT), intent(in)      :: magneticField
        class(StreamlineTracerUDT), intent(in)  :: streamlinetracer 
        type(TopomeshOptionsUDT), intent(in)    :: options

        ! Auxiliary
        integer(I8)                             :: thisf
        integer(I8), allocatable, dimension(:)  :: tf, contourtypes, &
            fsIDs, tfnb
        real(R8)                                :: thisfval, traceval
        real(R8), allocatable, dimension(:)     :: xint, yint, tx, ty, &
            tfval, xtrace, ytrace
        logical                                 :: skipvertex
        logical, allocatable, dimension(:)      :: keepind, hasfv
        type(ContourUDT), allocatable           :: tc(:), allc(:)
        type(PolygonUDT)                        :: tcp 
        ! real(R8), allocatable, dimension(:)     :: 

        ! Loop
        integer(I8)                             :: i, j, k

        ! Initialize
        !============
        ! Associate
        associate(&
            nfs             => topomesh%nfs,    &
            face            => topomesh%face,   &
            vert            => topomesh%vert    &
        )

        ! Initialize
        allocate(allc(0), contourtypes(0), fsIDs(0))

        ! X-points
        !=========
        if (options%addxpstreamlines) then 
            ! Loop over all vertices, check for saddle points
            do i = 1, nfs 
                ! Check if the current vertex is a saddle point
                if (.not. (vert%type(i) == TMvertexsaddleID)) then 
                    cycle
                end if 

                ! For this saddle point, check which faces it has
                tf = vert%GetFace(i)
                
                ! Loop over all faces and find face neighbours. If the
                ! left neighbour and current face are both separatrix
                ! parts, trace
                do j = 1, size(tf) 
                    ! Check
                    if (.not. face%type(tf(j)) == TMfacesepID) then 
                        cycle
                    end if

                    ! Get neighbours
                    tfnb = vert%GetFaceNeig(i, tf(j))

                    ! Sanity check
                    if (size(tfnb) /= 2) then 
                        print *, 'vertex number: ', i
                        call gdErrorHandler('AddTopologicalMeshStreamlines: ' // & 
                            'could not find two face neighbours for vertex')
                    end if 

                    ! Check next face
                    if (.not. face%type(tfnb(1)) == TMfacesepID) then 
                        cycle
                    end if 

                    ! Compute starting point
                end do 

                ! Check how many 
            end do 
        end if 

        ! End association
        end associate

        ! Add to the topology mesh
        !=========================
        ! Add contours
        do  i = 1, size(allc)
            call InsertTopologicalMeshContour(topomesh, magneticField, &
                allc(i), contourtypes(i), fsIDs(i))
        end do 

        ! Trim the topological mesh
        call TrimTopologicalMesh(topomesh, magneticField, vessel)

        ! Split boundaries
        call SplitTopologicalMeshFaces(topomesh)

    end subroutine
#endif 

    ! Limiter region removal
    subroutine RemoveTopologicalMeshLimiterRegions(topomesh, &
        magneticField, fieldtracer, streamlinetracer, options)

        ! Description
        !============
        ! This routine removes limiter-like regions (i.e. regions that 
        ! form disc cells, but that do not have any radial face, so 
        ! basically a minimum/maximum without radial line) by 
        ! introducing a connection between an extremum and a limiting
        ! tangency point. This connection is a streamline if it is
        ! found, otherwise it is a straight line between extremum and
        ! tangency point (this may perform poorly though, a warning will
        ! be issued). It is assumed that the necessary contours have 
        ! already been added, but that no core or other unnecessary 
        ! contours are present, nor that any regions have been removed. 
        ! We do assume that all remaining points are within the vessel
        ! boundary. 

        ! The algorithm is as follows: 
        ! 1) check which vertices are extrema without any faces 
        ! -> these are extrema in limiter-like configurations. 
        ! 2) retrace from these extrema a df/dx line (or df/dy line) 
        ! 3) check with which tangency point contour this line crosses 
        ! first -> that is the limiting tangency point
        ! 4) trace from the tangency point a streamline to the extremum
        ! 5) add that line as a face to the topological mesh

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh
        type(MagneticFieldUDT), intent(in)      :: magneticField
        class(ContourTracerUDT), intent(inout)  :: fieldtracer
        class(StreamlineTracerUDT), intent(in)  :: streamlinetracer
        class(TopomeshOptionsUDT), intent(in)   :: options 

        ! Auxiliary
        logical                                 :: hasfaces
        logical, allocatable, dimension(:)      :: islimiterextremum
        integer(I8)                             :: ngp, nle, faceID, &
            thistp
        integer(I8), allocatable, dimension(:)  :: sface, sfc, tfID, &
            tpfsID
        real(R8)                                :: xb(1:2), yb(1:2)
        real(R8), allocatable, dimension(:)     :: fx, fy, xg, yg, &
            xint, yint, srface, srfc, tsrfc
        class(ContourTracerUDT), allocatable    :: fxtracer, fytracer
        type(ContourUDT), allocatable           :: fxc(:), fyc(:), fc(:)
        type(StreamlineUDT), allocatable        :: streamlines(:)
        type(IntegerDynamicArrayUDT)            :: intfID
        type(RealDynamicArrayUDT)               :: intsrfc
        class(RealDynamicArrayUDT), allocatable :: tx, ty
        class(StreamlineTracerUDT), allocatable :: mystreamlinetracer

        ! Loop
        integer(I8)                             :: i, j

        ! Initialize
        !===========
        ! Associate for ease
        associate(&
            vert        => topomesh%vert, &
            face        => topomesh%face)        

        ! Initialize
        allocate(islimiterextremum(vert%ntot))
        islimiterextremum = .false. 

        ! Tracer bounds
        xb(1:2) = [minval(magneticField%R), maxval(magneticField%R)]
        yb(1:2) = [minval(magneticField%Z), maxval(magneticField%Z)]

        ! Copy and adjust tracer
        mystreamlinetracer = streamlinetracer 

        ! Set step size
        mystreamlinetracer%step = 0.1
        mystreamlinetracer%nsteps = 2000

        ! Determine limiter exterma
        !==========================
        ! Loop over all vertices
        do i = 1, vert%ntot 
            ! Check if extremum
            if (.not. any(topomesh%vert%type(i) == [TMvertexmaxID, TMvertexminID])) then 
                ! Skip
                cycle 
            end if 

            ! Check if this extremum has faces
            j = 1
            hasfaces = .false. 
            do while (j <= face%ntot)
                if (any(face%vert(j, :) == i)) then 
                    hasfaces = .true.
                    exit
                else
                    j = j + 1
                end if 
            end do 
            if (hasfaces) then 
                ! Skip
                cycle 
            end if 

            ! Extremum without faces found, set to true
            islimiterextremum(i) = .true.
        end do 

        ! Check if we should continue
        nle = count(islimiterextremum)
        if (nle == 0) then 
            return 
        end if 

        ! Initialize tracers
        !===================
        ! If we got here, limiter extrema were found. Initializing tracers
        ! Extract grid coordinates
        call fieldtracer%GetCoordinates(xg, yg)
        ngp = size(xg)

        ! Evaluate
        allocate(fx(ngp), fy(ngp))
        call magneticField%interp%Evaluate(xg, yg, 1, 0, fx)
        call magneticField%interp%Evaluate(xg, yg, 0, 1, fy)

        ! Construct new tracers - set only limiter extrema as saddle points
        fxtracer = fieldtracer
        fytracer = fieldtracer 
        if (allocated(fxtracer%xs)) then 
            deallocate(fxtracer%xs, fxtracer%ys, fxtracer%vs, &
                fxtracer%order, fxtracer%IDs)
        end if 
        allocate(fxtracer%xs(0), fxtracer%ys(0), fxtracer%vs(0), &
            fxtracer%order(0), fxtracer%IDs(0))
        if (allocated(fytracer%xs)) then 
            deallocate(fytracer%xs, fytracer%ys, fytracer%vs, &
            fytracer%order, fytracer%IDs)
        end if 
        allocate(fytracer%xs(0), fytracer%ys(0), fytracer%vs(0), &
            fytracer%order(0), fytracer%IDs(0))
        call fxtracer%SetValues(fx)
        call fytracer%SetValues(fy)

        ! Construct faces
        !================
        ! Get all flux surface IDs that belong to type 1 tangency points
        allocate(tpfsID(count(vert%type == TMvertextp2ID)))
        tpfsID = pack(vert%fsID, vert%type == TMvertextp2ID)

        ! Loop over all vertices
        do i = 1, vert%ntot
            ! Check if this is a limiter extremum
            if (.not. islimiterextremum(i)) then 
                ! Skip
                cycle
            end if 

            ! Trace the contour lines from this point
            fxc = fxtracer%TraceContours([vert%x(i)], [vert%y(i)]) ! not exactly df/dx = 0 but ok
            fyc = fytracer%TraceContours([vert%x(i)], [vert%y(i)])

            ! Hedge for possible non-allocation 
            if (.not. allocated(fxc)) then 
                allocate(fxc(0))
            end if 
            if (.not. allocated(fyc)) then 
                allocate(fyc(0))
            end if

            ! Keep only lines that start/end in this 

            ! Check
            if (size(fxc) == 0 .and. size(fyc) == 0) then 
                ! Strange, no lines found
                call gdErrorHandler('RemoveTopologicalMeshLimiterRegions: ' // &
                    'could not find any df/dx=0 or df/dy=0 lines in extremum, ' // &
                    'unexpected. Could be a bug')
            end if 
            fc = [fxc, fyc]

            ! For the first of these lines, find intersections with any 
            ! tangency point contour lines
            intfID = ConstructIntegerDynamicArray()
            intsrfc = ConstructRealDynamicArray()
            do j = 1, face%ntot
                if (any(face%fsID(j) == tpfsID)) then 
                    ! Look for intersections
                    call SimplePolygonIntersections(face%x(j)%Get(), &
                        face%y(j)%Get(), fc(1)%x, fc(1)%y, xint, yint, &
                        sface, sfc, srface, srfc)

                    ! Check if any found
                    if (size(xint) > 0) then 
                        ! Append 
                        call intfID%Append(j)
                        call intsrfc%Append(minval(srfc)) ! just take minimal value
                    end if 
                end if 
            end do 

            ! If no intersections found -> throw error
            if (intfID%Size() == 0) then 
                call gdErrorHandler('RemoveTopologicalMeshLimiterRegions: ' // & 
                    'no intersections with any tangency point contour found, this ' // & 
                    'may be a bug')
            end if 

            ! Find intersection closest to extremum (i.e. the one with 
            ! smallest segment index)
            tsrfc = intsrfc%Get()
            tfID = intfID%Get()
            faceID = tfID(minloc(tsrfc, 1))

            ! Find the corresponding tangency point
            thistp = findloc(vert%type == TMvertextp2ID .and. &
                vert%fsID == face%fsID(faceID), .true., 1)

            ! Trace a contour line from this point to the extremum
            if (vert%type(i) == TMvertexminID) then 
                ! Need to go downhill, so against the gradient
                streamlines = mystreamlinetracer%TraceStreamlines([vert%x(thistp)], &
                    [vert%y(thistp)], xb, yb, [-1_I8])
            else
                ! Need to go uphill, with the gradient
                streamlines = mystreamlinetracer%TraceStreamlines([vert%x(thistp)], &
                    [vert%y(thistp)], xb, yb, [1_I8])
            end if 

            ! Checks
            if (size(streamlines) /= 1) then 
                ! Print warning and simply connect tangency point and extremum
                print *, 'extremum vertex: ', i, 'tangeny vertex: ', thistp
                print *, 'RemoveTopologicalMeshLimiterRegions: ' // & 
                    'could not trace streamline from tangency point to ' // &
                    'extremum, adding straight line between these points'

                ! Add straight line
                tx = ConstructRealDynamicArray([vert%x(thistp), vert%x(i)])
                ty = ConstructRealDynamicArray([vert%y(thistp), vert%y(i)])
                call AddTopologicalMeshFace(topomesh, [thistp, i], &
                    tx, ty, TMfaceradID, 0_I8, 0.0_R8)
            else
                ! Append the extremum to the streamline
                tx = ConstructRealDynamicArray([streamlines(1)%x, vert%x(i)])
                ty = ConstructRealDynamicArray([streamlines(1)%y, vert%y(i)])

                ! Add
                call AddTopologicalMeshFace(topomesh, [thistp, i], &
                    tx, ty, TMfaceradID, 0_I8, 0.0_R8)
            end if 
        end do

        ! Housekeeping
        end associate

    end subroutine
    
    ! Core region removal
    subroutine RemoveTopologicalMeshCoreRegions(topomesh)

        ! Description
        !============
        ! This routine removes extrema in the topological mesh (vertices type 1 or
        ! 3) and removes all their faces and regions. This routine is useful if one wants to
        ! generate grids that do not extend to the extrema in the field, but only a
        ! certain fraction (if used in combination with
        ! AddTopologicalMeshCoreBoundaries). 

        ! This routine should be used after the topological mesh has been fully
        ! constructed. 

        ! Notes 
        !======

        ! Declare variables
        !==================
        ! Arguments 
        class(TopomeshUDT)                      :: topomesh 

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: tv
        logical, allocatable, dimension(:)      :: delv, delf, delc

        ! Loop
        integer(I8)                         :: i

        ! Initialize
        !===========
        ! Mark vertices for deletion
        delv = (topomesh%vert%type == TMvertexmaxID) .or. &
            (topomesh%vert%type == TMvertexminID)

        ! Mark faces for deletion
        allocate(delf(topomesh%face%ntot))
        delf = .false. 
        do i = 1, topomesh%face%ntot
            if (any(delv(topomesh%face%vert(i, :)))) then 
                delf(i) = .true.
            end if 
        end do

        ! Mark cells for deletion
        allocate(delc(topomesh%cell%ntot))
        delc = .false.
        do i = 1, topomesh%cell%ntot
            tv = GetTMCellVert(topomesh%cell, i)
            if (any(delv(tv))) then 
                delc(i) = .true.
            end if 
        end do 

        ! Delete
        !=======
        ! Vertices
        call RemoveTopologicalMeshVertexLogical(topomesh, delv)

        ! Faces
        call RemoveTopologicalMeshFaceLogical(topomesh, delf)

        ! Cells
        call RemoveTopologicalMeshCellLogical(topomesh, delc)

        ! Update again
        !=============
        ! Vertex faces
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Data
        call AddTopologicalMeshData(topomesh)

    end subroutine

    ! Wide grid region removal
    subroutine RemoveTopologicalMeshWideGridRegions(topomesh)

        ! Description
        !============
        ! This routine removes any regions that are classified as 
        ! 'wide grid' regions. In practice, this means any region that 
        ! does not have a separatrix boundary segment (note that if 
        ! core boundaries are inserted, also the core region will be 
        ! removed, if not removed already)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)      :: topomesh 

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: tcf, tcv
        logical, allocatable, dimension(:)      :: delv, delf, delc

        ! Loop
        integer(I8)                         :: i

        ! Initialize
        !===========
        ! Mark cells, vertices, faces for deletion
        allocate(delc(topomesh%cell%ntot), delv(topomesh%vert%ntot), &
            delf(topomesh%face%ntot))
        delc = .true. ! default true
        delv = .true.
        delf = .true.

        do i = 1, topomesh%cell%ntot 
            ! Get faces
            tcf = topomesh%cell%GetFace(i)

            ! Check if there is a separatrix face, then keep cell
            if (any(topomesh%face%type(tcf) == TMfacesepID)) then 
                delc(i)     = .false. 
                delf(tcf)   = .false.
                tcv         = topomesh%cell%GetVert(i)
                delv(tcv)   = .false.
            end if 
        end do

        ! Delete
        !=======
        ! Vertices
        call RemoveTopologicalMeshVertexLogical(topomesh, delv)

        ! Faces
        call RemoveTopologicalMeshFaceLogical(topomesh, delf)

        ! Cells
        call RemoveTopologicalMeshCellLogical(topomesh, delc)

        ! Update again
        !=============
        ! Vertex faces
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Data
        call AddTopologicalMeshData(topomesh)

    end subroutine

    ! Remove all-but core regions (for Anthony)
    subroutine RemoveTopologicalMeshNonCoreRegions(topomesh)

        ! Description
        !============
        ! This routine removes all regions that don't have either a 
        ! core boundary (inserted) or a maximum/minimum as boundary.
        ! This should yield the desired output for devices such as 
        ! e.g. TOMAS that are almost perfectly circular

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)              :: topomesh 

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: tcf, tcv
        logical, allocatable, dimension(:)      :: delv, delf, delc

        ! Loop
        integer(I8)                         :: i

        ! Initialize
        !===========
        ! Mark cells, vertices, faces for deletion
        allocate(delc(topomesh%cell%ntot), delv(topomesh%vert%ntot), &
            delf(topomesh%face%ntot))
        delc = .true. ! default true
        delv = .true.
        delf = .true.

        do i = 1, topomesh%cell%ntot 
            ! Get faces
            tcf = topomesh%cell%GetFace(i)

            ! Get vertices
            tcv = topomesh%cell%GetVert(i)

            ! Check if there is a core face or extremum, then keep cell
            if (any(topomesh%face%type(tcf) == TMfacecoreID) .or. &
                any(topomesh%vert%type(tcv) == TMvertexminID) .or. &
                any(topomesh%vert%type(tcv) == TMvertexmaxID)) then 
                delc(i)     = .false. 
                delf(tcf)   = .false.
                delv(tcv)   = .false.
            end if 
        end do

        ! Delete
        !=======
        ! Vertices
        call RemoveTopologicalMeshVertexLogical(topomesh, delv)

        ! Faces
        call RemoveTopologicalMeshFaceLogical(topomesh, delf)

        ! Cells
        call RemoveTopologicalMeshCellLogical(topomesh, delc)

        ! Update again
        !=============
        ! Vertex faces
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Data
        call AddTopologicalMeshData(topomesh)

    end subroutine

    ! Flux tube merging
    subroutine MergeTopologicalMeshFluxTubes(topomesh, magneticField, &
        vessel, fieldtracer, options)

        ! Description
        !============
        ! This routine merges topological mesh tubes if it is possible 
        ! and if it is desired by the user. Typically, only two adjacent
        ! tubes are merged at the same time. This operation results in 
        ! the deletion of cells, vertices, and flux surfaces, and the 
        ! construction of new cells. To apply this operation, the 
        ! topological mesh must be fully constructed with all additional
        ! interconnection information. This information is later 
        ! reconstructed, meaning that cells are retraced etc (this is 
        ! actually done after each tube merging, otherwise it's 
        ! impossible to proceed). 

        ! The currently supported types of merging are:
        ! - merging of tubes that share a tangency point contour as 
        !   boundary, and where one side has a single tube and the other
        !   side has two tubes, of which at least one is a single cell 
        !   tube with a tangency point (type 1) boundary (i.e. no real
        !   other flux surface boundary but a point). This sounds awfully
        !   specific, and it is, but it represents one of the most 
        !   common cases that requires coarsening (especially at nearly
        !   aligned walls etc that have many tangency points)

        ! See dedicated subroutines/documentation below for more information
        
        ! Algorithm
        !==========
        ! Merging of tangency point contours:
        !   1)  loop over all tubes until we reached the final tube and 
        !       did not apply a merging operation
        !   2)  Check if the current tube is eligible (a tube is eligible
        !       if the maximal difference in psi value of its vertices is
        !       below the specified threshold)
        !   3)  Check if we can merge the tube by deleting one of its 
        !       boundary face sets. This is possible if the tube only has
        !       two tube neighbours, and at least one neigbour is a 
        !       boundary tube (i.e. it has no boundary faces at one
        !       side)
        !   4)  To merge, the bounding faces are removed from the 
        !       topological mesh. The boundary vertices of these faces are
        !       kept, but are set to the regular boundary type. 
        !       Non-boundary vertices are deleted.
        !       If both tubes are boundary tubes, we can only keep 
        !       one of the tangency points that form the other boundary.
        !       The other one becomes a regular boundary vertex. This
        !       will allow to apply the SimplifyTopologicalMesh operation
        !       to simplify boundaries. The choice of which vertex to 
        !       take will 
        !   5)  All cell data etc is reconstructed after each tube merge

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT), intent(inout)       :: topomesh
        type(TopomeshOptionsUDT), intent(in)    :: options 
        type(MagneticFieldUDT),intent(in)       :: magneticField 
        type(VesselUDT), intent(inout)          :: vessel
        class(ContourTracerUDT), intent(inout)  :: fieldtracer

        ! Auxiliary
        logical                                 :: marked, &
            passedcheck
        logical, allocatable, dimension(:)      :: delf, delv
        integer(I8)                             :: maxind
        integer(I8), allocatable, dimension(:)  :: tf, tfv, tnb, &
            tfmerge, tnbmerge, tfvu, tfradmerge, tvfvID, tvfvIDu, &
            tvf, thisv
        real(R8)                                :: dpsi, dpsinb1, dpsinb2, &
            thisdeletedfval, lrad, lradnb1, lradnb2
        real(R8), allocatable, dimension(:)     :: fval, thisvfval

        ! Loop
        integer(I8)                             :: i, j

        ! Initialize
        !===========
        ! Associate
        associate(&
            vert        => topomesh%vert,   &
            face        => topomesh%face,   &
            cell        => topomesh%cell,   &
            tube        => topomesh%tube    &
            )

        

        ! Merge tubes
        !============
        ! Tangency point tubes
        !---------------------
        if (options%mergetangencypointtubes) then 
            ! Loop over all tubes to find eligible tubes for merging
            i = 1
            do while (i <= tube%ntot)
                
                ! Initialize
                marked = .false. 

                ! Set logicals
                if (allocated(delv)) then 
                    deallocate(delv, delf)
                end if
                allocate(delv(vert%ntot), delf(face%ntot))
                delv = .false.
                delf = .false.

                ! See if the tube can be merged on side one
                passedcheck = .true.
                if (.not. marked .and. (tube%ftneig1P(i, 2) > 0)) then 
                    ! Get the tube neighbours on that side
                    tnb = tube%GetNeig(i, 1_I8)

                    ! Check the amount of neighbours
                    if (size(tnb) == 2) then 

                        ! Check if at least one neighbour is a boundary 
                        ! tube (i.e. it doesn't have any other boundaries)
                        if ((tube%ftneig1P(tnb(1), 2) /= 0 .and. tube%ftneig2P(tnb(1), 2) /= 0) .and. &
                            (tube%ftneig1P(tnb(2), 2) /= 0 .and. tube%ftneig2P(tnb(2), 2) /= 0)) then
                            passedcheck = .false.  
                        end if

                        ! Check tube dpsi and radial length
                        tf = tube%GetFace(i)
                        dpsi = maxval(GetTMFaceDeltaPsi(topomesh, tf))
                        lrad = maxval(GetTMFaceRadialLength(topomesh, magneticField, tf))
                        tf = tube%GetFace(tnb(1))
                        dpsinb1 = maxval(GetTMFaceDeltaPsi(topomesh, tf))
                        lradnb1 = maxval(GetTMFaceRadialLength(topomesh, magneticField, tf))
                        tf = tube%GetFace(tnb(2))
                        dpsinb2 = maxval(GetTMFaceDeltaPsi(topomesh, tf))
                        lradnb2 = maxval(GetTMFaceRadialLength(topomesh, magneticField, tf))

                        ! Check if we should merge
                        if (all([dpsi, dpsinb1, dpsinb2] >= options%dpsimintangencypointtubes) .and. &
                            all([lrad, lradnb1, lradnb2] >= options%lradmintangencypointtubes)) then 
                            passedcheck = .false.
                        end if 

                        ! Check if we can merge
                        if (passedcheck) then 
                            ! Get merge data
                            tfmerge = tube%GetBndFace(i, 1_I8)
                            tnbmerge = tnb
                            tfradmerge = [tube%GetFace(tnb(1)), &
                                tube%GetFace(tnb(2)), tube%GetFace(i)]

                            ! Check for non-mergeable surfaces (separatrix basicall)
                            if (any(face%type(tfmerge) == TMfacesepID)) then 
                                marked = .false.
                            else
                                marked = .true.
                            end if

                        end if 

                    elseif (size(tnb) == 1) then 

                        ! Check if this neighbour has maximally one
                        ! other neighbour on each side
                        if ((tube%ftneig1P(tnb(1), 2) > 1 .or. tube%ftneig2P(tnb(1), 2) > 1)) then
                            passedcheck = .false.  
                        end if

                        ! Check tube dpsi and radial length
                        tf = tube%GetFace(i)
                        dpsi = maxval(GetTMFaceDeltaPsi(topomesh, tf))
                        lrad = maxval(GetTMFaceRadialLength(topomesh, magneticField, tf))
                        tf = tube%GetFace(tnb(1))
                        dpsinb1 = maxval(GetTMFaceDeltaPsi(topomesh, tf))
                        lradnb1 = maxval(GetTMFaceRadialLength(topomesh, magneticField, tf))

                        if (all([dpsi, dpsinb1] >= options%dpsimintangencypointtubes) .and. &
                            all([lrad, lradnb1] >= options%lradmintangencypointtubes)) then 
                            passedcheck = .false.
                        end if 

                        ! Check if we can merge
                        if (passedcheck) then 

                            ! Get merge data
                            tfmerge = tube%GetBndFace(i, 1_I8)
                            tnbmerge = tnb
                            tfradmerge = [tube%GetFace(tnb(1)), tube%GetFace(i)]

                            ! Check for non-mergeable surfaces (separatrix basicall)
                            if (any(face%type(tfmerge) == TMfacesepID)) then 
                                marked = .false.
                            else
                                marked = .true.
                            end if

                        end if 
                    end if 
                end if 

                ! See if the tube can be merged on side two
                if (.not. marked .and. (tube%ftneig2P(i, 2) > 0)) then 
                    ! Get the tube neighbours on that side
                    tnb = tube%GetNeig(i, 2_I8)

                    ! Check the amount of neighbours
                    if (size(tnb) == 2) then 

                        ! Check if at least one neighbour is a boundary 
                        ! tube (i.e. it doesn't have any other boundaries)
                        if ((tube%ftneig1P(tnb(1), 2) /= 0 .and. tube%ftneig2P(tnb(1), 2) /= 0) .and. &
                            (tube%ftneig1P(tnb(2), 2) /= 0 .and. tube%ftneig2P(tnb(2), 2) /= 0)) then
                            passedcheck = .false.  
                        end if

                        ! Check tube dpsi and radial length
                        tf = tube%GetFace(i)
                        dpsi = maxval(GetTMFaceDeltaPsi(topomesh, tf))
                        lrad = maxval(GetTMFaceRadialLength(topomesh, magneticField, tf))
                        tf = tube%GetFace(tnb(1))
                        dpsinb1 = maxval(GetTMFaceDeltaPsi(topomesh, tf))
                        lradnb1 = maxval(GetTMFaceRadialLength(topomesh, magneticField, tf))
                        tf = tube%GetFace(tnb(2))
                        dpsinb2 = maxval(GetTMFaceDeltaPsi(topomesh, tf))
                        lradnb2 = maxval(GetTMFaceRadialLength(topomesh, magneticField, tf))

                        ! Check if we should merge
                        if (all([dpsi, dpsinb1, dpsinb2] >= options%dpsimintangencypointtubes) .and. &
                            all([lrad, lradnb1, lradnb2] >= options%lradmintangencypointtubes)) then 
                            passedcheck = .false.
                        end if 

                        ! Check if we can merge
                        if (passedcheck) then 
                            ! Get merge data
                            tfmerge = tube%GetBndFace(i, 2_I8)
                            tnbmerge = tnb
                            tfradmerge = [tube%GetFace(tnb(1)), &
                                tube%GetFace(tnb(2)), tube%GetFace(i)]

                            ! Check for non-mergeable surfaces (separatrix basicall)
                            if (any(face%type(tfmerge) == TMfacesepID)) then 
                                marked = .false.
                            else
                                marked = .true.
                            end if

                        end if 

                    elseif (size(tnb) == 1) then 

                        ! Check if this neighbour has maximally one
                        ! other neighbour on each side
                        if ((tube%ftneig1P(tnb(1), 2) > 1 .or. tube%ftneig2P(tnb(1), 2) > 1)) then
                            passedcheck = .false.  
                        end if

                        ! Check tube dpsi and radial length
                        tf = tube%GetFace(i)
                        dpsi = maxval(GetTMFaceDeltaPsi(topomesh, tf))
                        lrad = maxval(GetTMFaceRadialLength(topomesh, magneticField, tf))
                        tf = tube%GetFace(tnb(1))
                        dpsinb1 = maxval(GetTMFaceDeltaPsi(topomesh, tf))
                        lradnb1 = maxval(GetTMFaceRadialLength(topomesh, magneticField, tf))

                        if (all([dpsi, dpsinb1] >= options%dpsimintangencypointtubes) .and. &
                            all([lrad, lradnb1] >= options%lradmintangencypointtubes)) then 
                            passedcheck = .false.
                        end if 

                        ! Check if we can merge
                        if (passedcheck) then 

                            ! Get merge data
                            tfmerge = tube%GetBndFace(i, 2_I8)
                            tnbmerge = tnb
                            tfradmerge = [tube%GetFace(tnb(1)), tube%GetFace(i)]

                            ! Check for non-mergeable surfaces (separatrix basicall)
                            if (any(face%type(tfmerge) == TMfacesepID)) then 
                                marked = .false.
                            else
                                marked = .true.
                            end if

                        end if 
                    end if 
                end if

                ! If marked for merging, apply merging operation
                if (marked) then 
                    ! Reset counter 
                    i = 1

                    ! Mark faces and vertices for removal
                    delf = .false. 
                    delv = .false.
                    delf(tfmerge) = .true. 
                    do j = 1, size(tfmerge)
                        ! Get face vertices
                        tfv = face%vert(tfmerge(j), :)

                        ! Check vertex type
                        if (any(vert%type(tfv(1)) == [TMvertexbndID, TMvertextp2ID])) then 
                            ! Don't delete, but reset type
                            vert%type(tfv(1)) = TMvertexbndID
                        else
                            ! Mark for deletion
                            delv(tfv(1)) = .true.
                        end if
                        if (any(vert%type(tfv(2)) == [TMvertexbndID, TMvertextp2ID])) then 
                            ! Don't delete, but reset type
                            vert%type(tfv(2)) = TMvertexbndID
                        else
                            ! Mark for deletion
                            delv(tfv(2)) = .true.
                        end if
                    end do 

                    ! Check if two (or more) type 1 tangency points were present.
                    ! If so, delete one of both (perhaps)
                    call Unique([face%vert(tfradmerge, 1), face%vert(tfradmerge, 2)], tfvu)
                    do while (count(vert%type(tfvu) == TMvertextp1ID) > 1)
                        ! Get all type 1 vertices
                        allocate(thisv(count(vert%type(tfvu) == TMvertextp1ID)))
                        thisv = pack(tfvu, vert%type(tfvu) == TMvertextp1ID)

                        ! Get field values
                        thisvfval = topomesh%fsfval%Get(vert%fsID(thisv))
                        thisdeletedfval = topomesh%fsfval%Get(face%fsID(tfmerge(1)))

                        ! Check which vertex is the furthest away in terms
                        ! of psi values - keep that one, delete the rest
                        maxind = maxloc(abs(thisvfval - thisdeletedfval), 1)
                        vert%type(thisv) = TMvertexbndID
                        vert%type(thisv(maxind)) = TMvertextp1ID

                        ! Housekeeping
                        deallocate(thisv)
                    end do 

                    ! Remove faces
                    call RemoveTopologicalMeshFaceLogical(topomesh, delf)

                    ! Remove vertices
                    call RemoveTopologicalMeshVertexLogical(topomesh, delv)

                    ! Simplify
                    call SimplifyTopologicalMeshFaces(topomesh)

                    !!! To be replaced by garbage point collector?

                    ! Check if the tangency points are currently located
                    ! in between two flux surfaces with different ID. 
                    ! In that case, 'delete' one of them
                    j = 1
                    do while (j <= vert%ntot)
                        if (vert%type(j) == TMvertextp1ID) then 
                            ! Check
                            tvf = vert%GetFace(j)

                            ! Sanity check
                            if (size(tvf) /= 2) then 
                                ! Write out some data for debugging
                                call WriteTopologicalMesh(topomesh, 'topomesh_error')
                                print *, 'vertex: ', j
                                call gdErrorHandler('MergeTopologicalMeshFluxTubes: ' // & 
                                    'type 1 tangency point does not have ' // & 
                                    'exactly two faces, unexpected. Check mesh in topomesh_error.dat')
                            end if 

                            ! Get all unique flux surface IDs of these vertices
                            tvfvID = vert%fsID([face%vert(tvf, 1), face%vert(tvf, 2)])
                            call Unique(tvfvID, tvfvIDu)
                            if (size(tvfvIDu) > 2) then 
                                ! Set vertex type to simple boundary vertex
                                ! for deletion later on in simplification 
                                ! step
                                vert%type(j) = TMvertexbndID

                                ! Simplify
                                call SimplifyTopologicalMeshFaces(topomesh)

                                ! Reset counter 
                                j = 0
                            end if 
                        end if 

                        ! Update counter
                        j = j + 1
                    end do 

                    ! Simplify again
                    call SimplifyTopologicalMeshFaces(topomesh)

                    ! Split faces if necessary
                    call SplitTopologicalMeshFaces(topomesh)   

                    ! Add contours for split boundary vertices in case they emerge
                    call AddBoundarySplitVertexContours(topomesh, &
                        magneticField, vessel, fieldtracer)

                    call WriteTopologicalMesh(topomesh, 'topomesh_temp')

                    ! Recompute all interconnections, cells, etc
                    ! Vertex faces
                    call AddTopologicalMeshVertexFaces(topomesh)

                    ! Data
                    call AddTopologicalMeshData(topomesh)

                    ! Add cells
                    call AddTopologicalMeshCells(topomesh)

                    ! Data (recompute)
                    call AddTopologicalMeshData(topomesh)

                    ! Compute interconnection data
                    call AddTopologicalMeshInterconnectionData(topomesh)

                    call WriteTopologicalMesh(topomesh, 'topomesh_temp')

                else
                    i = i + 1
                end if 

                ! Housekeeping
                deallocate(delv, delf)
                
            end do
        end if 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Contour insertion into topological mesh
    subroutine InsertTopologicalMeshContour(topomesh, magneticField, contour, &
        contourtype, contourfsID) 

        ! Description
        !============
        ! This routine provides a general way to insert a curve into an existing
        ! topological mesh. The topological mesh should contain vertices and faces
        ! already, but no cells yet (these are also not updated here). The curve
        ! should have at least x and y data and, if available, start and/or end
        ! vertex IDs of vertices that occur in the topological mesh. If new
        ! vertices should be introduced, this should be done in a separate routine
        ! beforehand! 

        ! We then compute all intersections of this curve with all already
        ! available faces in the domain. Faces that are intersected in this way
        ! will be deleted and the new faces and intersections will be added.

        ! IMPORTANT
        !==========
        ! It should be noted that not all exceptional intersection cases are
        ! checked (see also the notes below for more information) so use this
        ! routine wisely. 

        ! Notes
        !======
        ! Note 1: if a curve part does not start or end in a new intersection or
        ! existing vertex, it will still be added as a face here. Clean-up should
        ! be done later on, as it is possible that after adding more curves,
        ! additional intersections are found with these segments, which can't be
        ! found if we remove them prematurely. 

        ! Note 2: we assume that if the start or end vertex of a curve/face is
        ! present, the first/last x, y coordinate of the curve is exactly equal to
        ! that vertex's coordinates. 

        ! Note 3: it is assumed that the curves form simple polygons

        ! Note 4: there may still be exceptional intersection cases which are not
        ! properly captured by this routine. To avoid these as much as possible, it
        ! is good practice to compute these special intersection points beforehand.
        ! Typically these are saddle points, extrema, tangency points, ... When
        ! starting the curve from this point, which is assumed here, only checks
        ! need to be done whether the starting or ending point coincides with an
        ! existing one. 

        ! Algorithm
        !==========
        ! 
        ! 1) For each curve, we compute the intersections with all existing 
        !    topological faces. For each topological face, do:
        !       1.1) Compute intersections using standard polygon intersection routine
        !       1.2) All intersections that have been found are checked whether
        !       they are close (up to precision defined below) to an existing
        !       vertex. If they are not close, they are added as new vertices. If
        !       they are close, it is assumed that the intersection happens at one
        !       of the existing vertices and no new vertex is added. The
        !       intersection is then adjusted to be exactly this original vertex,
        !       and also the ID of the intersection is set to the ID of that
        !       vertex. 
        !       1.3) Now, we check the end points of the segments. If any end
        !       points coincide, we check the IDs of these points. If they are the
        !       same and non-zero, nothing must be done. If one is non-zero, then
        !       we update the zero ID to the non-zero one. If both are zero, we
        !       need to add a new vertex to the topology mesh. If they are
        !       non-zero, but not the same, we need to throw an error. 
        !       1.4) Add for each curve the start and end vertex also as
        !       intersection with updated ID. 
        ! 2) After computing all intersections and having partitioned the segments
        ! as stated above, we can remove all old faces and add the new ones. Note
        ! that, even if no intersections are found of an existing face, the
        ! algorithm above will add that face again as expected. This may result in
        ! some overhead, so it is best to call this function only once by
        ! precomputing all desired curves beforehand, if possible. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        type(MagneticFieldUDT), intent(in)      :: magneticField 
        type(ContourUDT), intent(inout)         :: contour
        integer(I8), intent(in)                 :: contourfsID, contourtype

        ! Auxiliary
        real(R8)                                :: dist, fsfval
        real(R8), allocatable, dimension(:)     :: xint, yint, s1r, s2r, &
            tfv, ts1r, tx, ty
        integer(I8)                             :: nint, fsID
        integer(I8), allocatable, dimension(:)  :: s1, s2, fID, vIDs, &
            vtypes, sortind, vf1, vf2, tvIDs, ts2, ts1
        logical                                 :: alreadyadded, &
            isinconsistent
        logical, allocatable, dimension(:)      :: iscoinciding, &
            delind, keepind
        type(PolygonUDT)                        :: cp 
        type(RealDynamicArrayUDT)               :: xda, yda, s1rda, &
            s2rda
        type(RealDynamicArrayUDT), allocatable  :: xfda(:), yfda(:)
        type(IntegerDynamicArrayUDT)            :: s1da, s2da, fda 

        ! Loop
        integer(I8)                             :: i, j, k, kold

        ! Initialize
        !===========
        ! Construct polygon from contour data
        call cp%Construct(contour%x, contour%y)

        ! Initialize dynamic arrays
        xda     = ConstructRealDynamicArray()
        yda     = ConstructRealDynamicArray()
        s1rda   = ConstructRealDynamicArray()
        s2rda   = ConstructRealDynamicArray()
        s1da    = ConstructIntegerDynamicArray()
        s2da    = ConstructIntegerDynamicArray()
        fda     = ConstructIntegerDynamicArray()

        ! Compute intersections 
        !======================
        ! Loop over all faces
        !$omp parallel do default(none) schedule(dynamic) &
        !$omp shared(topomesh, fda, xda, yda, s1da, s2da, s1rda, s2rda, cp) & 
        !$omp private(i, xint, yint, s1, s2, s1r, s2r)
        do i  = 1, topomesh%face%ntot 
            ! Associate current face polygon 
            associate(&
                fp          => topomesh%face%pol(i),    &
                fpx         => topomesh%face%pol(i)%x,  &
                fpy         => topomesh%face%pol(i)%y)

            ! Compute intersections
            call PolygonIntersections(fp, cp, xint, yint, s1, s2, &
                s1r=s1r, s2r=s2r)

            ! Store
            !$omp critical
            call fda%Append(spread(i, 1, size(xint)))
            call xda%Append(xint)
            call yda%Append(yint)
            call s1da%Append(s1)
            call s2da%Append(s2)
            call s1rda%Append(s1r)
            call s2rda%Append(s2r)
            !$omp end critical

            ! Housekeeping
            end associate
        end do
        !$omp end parallel do

        ! Extract
        fID     = fda%Get()
        xint    = xda%Get()
        yint    = yda%Get()
        s1      = s1da%Get()
        s2      = s2da%Get()
        s1r     = s1rda%Get()
        s2r     = s2rda%Get()
        nint    = size(xint)

        ! Add vertices
        !=============
        ! Initialize the vertex IDs and vertex types
        allocate(vIDs(nint), vtypes(nint))
        vIDs = 0
        vtypes = 0

        ! Loop over all intersections
        do i = 1, nint 
            ! Unpack
            associate(&
                xinti       => xint(i),     &
                yinti       => yint(i),     &
                fIDi        => fID(i),      &
                s1ri        => s1r(i),      &
                s2ri        => s2r(i)       &
                )

            ! Initialize logicals
            alreadyadded    = .false. 
            isinconsistent  = .false. 

            ! Set vertex type (default: 0)
            if ((contourtype == TMfacebndID) .or. &
                (topomesh%face%type(fIDi) == TMfacebndID)) then 
                ! Boundary vertex but no tangency point
                vtypes(i) = TMvertexbndID
            end if 

            ! Check if this vertex is already added
            iscoinciding = (abs(xinti - topomesh%vert%x) < disttol) .and. &
                (abs(yinti - topomesh%vert%y) < disttol)
            if (any (iscoinciding) ) then 
                ! Sanity check
                if (count(iscoinciding) > 1) then 
                    ! This shouldn't happen if the vertices in topomesh
                    ! are unique. Throw error
                    call gdErrorHandler('InsertTopologicalMeshContour: ' // & 
                        'duplicate vertices seem to appear in topomesh, ' // &
                        'check input')
                end if 

                ! Vertex coincides up to disttol precision -  do 
                ! sanity checks and adjustments of vertex ID
                alreadyadded = .true. 

                ! Was the intersection in a start point, if yes, do 
                ! vertex IDs correspond? 
                if ((s1ri == 0_R8) .and. (s2ri == 0_R8)) then 

                    ! Intersection in start of face and start of contour
                    if (contour%startsaddle /=0) then 
                        if (((topomesh%face%vert(fIDi, 1)) /= contour%startsaddle) .or. &
                            (topomesh%face%vert(fIDi, 1) == 0)) then 
                            ! Inconsistent - throw error later
                            isinconsistent = .true.
                        end if
                    else 
                        if (topomesh%face%vert(fIDi, 1) /= 0) then 
                            ! Not necessarily an error, but do print a warning
                            ! as this isn't expected
                            print *, 'vertex ID: ', topomesh%face%vert(fIDi, 1)
                            print *, 'InsertTopologicalMeshContour: ' // & 
                                'contour saddle point corresponds to ' // & 
                                'topomesh vertex but was not identified as ' // &
                                'such. Adding ID and continuing...'
                            contour%startsaddle = topomesh%face%vert(fIDi, 1)
                        else
                            ! Inconsistent
                            isinconsistent = .true.
                        end if 
                    end if 

                    ! Set ID 
                    vIDs(i) = contour%startsaddle

                elseif ((s1ri == 0_R8) .and. &
                    (s2ri == real(cp%ne, R8))) then 

                    ! Intersection in start of face and end of contour
                    if (contour%endsaddle /=0) then 
                        if (((topomesh%face%vert(fIDi, 1)) /= contour%endsaddle) .or. &
                            (topomesh%face%vert(fIDi, 1) == 0)) then 
                            ! Inconsistent - throw error later
                            isinconsistent = .true.
                        end if
                    else 
                        if (topomesh%face%vert(fIDi, 1) /= 0) then 
                            ! Not necessarily an error, but do print a warning
                            ! as this isn't expected
                            print *, 'vertex ID: ', topomesh%face%vert(fIDi, 1)
                            print *, 'InsertTopologicalMeshContour: ' // & 
                                'contour saddle point corresponds to ' // & 
                                'topomesh vertex but was not identified as ' // &
                                'such. Adding ID and continuing...'
                            contour%endsaddle = topomesh%face%vert(fIDi, 1)
                        else
                            ! Inconsistent
                            isinconsistent = .true.
                        end if 
                    end if 

                    ! Set ID 
                    vIDs(i) = contour%endsaddle

                elseif ((s1ri == real(topomesh%face%pol(fIDi)%ne, R8)) .and. &
                    (s2ri == 0_R8)) then 

                    ! Intersection in end of face and start of contour    
                    if (contour%startsaddle /=0) then 
                        if (((topomesh%face%vert(fIDi, 2)) /= contour%startsaddle) .or. &
                            (topomesh%face%vert(fIDi, 2) == 0)) then 
                            ! Inconsistent - throw error later
                            isinconsistent = .true.
                        end if
                    else 
                        if (topomesh%face%vert(fIDi, 2) /= 0) then 
                            ! Not necessarily an error, but do print a warning
                            ! as this isn't expected
                            print *, 'vertex ID: ', topomesh%face%vert(fIDi, 2)
                            print *, 'InsertTopologicalMeshContour: ' // & 
                                'contour saddle point corresponds to ' // & 
                                'topomesh vertex but was not identified as ' // &
                                'such. Adding ID and continuing...'
                            contour%startsaddle = topomesh%face%vert(fIDi, 2)
                        else
                            ! Inconsistent
                            isinconsistent = .true.
                        end if 
                    end if 

                    ! Set ID 
                    vIDs(i) = contour%startsaddle

                elseif ((s1ri == real(topomesh%face%pol(fIDi)%ne, R8)) .and. &
                    (s2ri == real(cp%ne, R8))) then 

                    ! Intersection in end of face and end of contour
                    if (contour%endsaddle /=0) then 
                        if (((topomesh%face%vert(fIDi, 2)) /= contour%endsaddle) .or. &
                            (topomesh%face%vert(fIDi, 2) == 0)) then 
                            ! Inconsistent - throw error later
                            isinconsistent = .true.
                        end if
                    else 
                        if (topomesh%face%vert(fIDi, 2) /= 0) then 
                            ! Not necessarily an error, but do print a warning
                            ! as this isn't expected
                            print *, 'vertex ID: ', topomesh%face%vert(fIDi, 2)
                            print *, 'InsertTopologicalMeshContour: ' // & 
                                'contour saddle point corresponds to ' // & 
                                'topomesh vertex but was not identified as ' // &
                                'such. Adding ID and continuing...'
                            contour%endsaddle = topomesh%face%vert(fIDi, 2)
                        else
                            ! Inconsistent
                            isinconsistent = .true.
                        end if 
                    end if 

                    ! Set ID 
                    vIDs(i) = contour%endsaddle

                else 

                    ! Apparently we're lucky and we get an intersection
                    ! in a vertex that already exists but that is 
                    ! not part of the start or end of a face or 
                    ! contour. This is very unlikely so we throw a 
                    ! warning yet continue and add the vertex ID 
                    print *, 'InsertTopologicalMeshContour: ' // & 
                        'found intersection that exactly coincides with ' // & 
                        'existing mesh vertex, yet was not identified ' // & 
                        'previously (i.e. the intersection happens to be )' // & 
                        'haphazardly together with an existing mesh vertex). ' // & 
                        'proceeding to add vertex, but results may be suprising'

                    ! Set the vertex ID 
                    vIDs(i) = findloc(iscoinciding, .true., 1)

                end if

                ! Set type
                vtypes(i) = topomesh%vert%type(vIDs(i))

            end if
            
            ! Check if we encountered an inconsistency, if so - call 
            ! error and exit. 
            if (isinconsistent) then 
                call gdErrorHandler('InsertTopologicalMeshContour: ' // &
                    'inconsistency encountered when adding intersection ' // & 
                    'with face number: ', fIDi)
            end if 

            ! Add the vertex (if not already present)
            if (.not. alreadyadded) then 
                ! Check flux surface ID
                fsID = 0
                if (contourfsID /= 0 ) then 
                    ! Sanity check
                    if (topomesh%face%fsID(fID(i)) /= 0) then 
                        ! Normally contours shouldn't intersect... 
                        print *, 'InsertTopologicalMeshContour: contour with fsID ', &
                            contourfsID, ' intersects with face ', fID(i), &
                            ' with fsID ', topomesh%face%fsID(fID(i)), ' - unexpected.' // & 
                            ' Taking face ID'
                        
                        ! Take face ID for contour
                        fsID = topomesh%face%fsID(fID(i))
                    else 
                        ! Take contour ID
                        fsID = contourfsID
                    end if 
                else 
                    ! Take face ID
                    fsID = topomesh%face%fsID(fID(i))
                end if 

                allocate(tfv(1))
                call magneticField%interp%Evaluate([xinti], [yinti], 0, 0, tfv)
                call AddTopologicalMeshVertex(topomesh, xinti, yinti, &
                    tfv(1), vtypes(i), fsID)
                ! Add ID as well 
                vIDs(i) = topomesh%vert%ntot
                deallocate(tfv)
            end if 
                
            ! Housekeeping
            end associate
        end do

        ! Add contour faces
        !==================
        ! Sort intersections according to contour coordinate
        allocate(sortind(size(s1r)))
        call Sort(s2r, ind=sortind)
        s1r = s1r(sortind)
        fID = fID(sortind)
        s1 = s1(sortind)
        s2 = s2(sortind)
        xint = xint(sortind)
        yint = yint(sortind)
        vIDs = vIDs(sortind)
        deallocate(sortind)

        ! Hedge for duplicate intersections (possible with closed polygons 
        ! or if multiple faces intersect in the same point)
        allocate(keepind(size(s1r)))
        keepind = .true. 
        do j = 1, size(keepind)-1
            if ((s2r(j+1) - s2r(j)) == 0_R8) then 
                keepind(j+1) = .false. 
            end if 
        end do

        ! Add start and end points as intersections if they have an 
        ! ID (and if that ID is not already present as an intersection)
        tvIDs = pack(vIDs, keepind) 
        ts2 = pack(s2, keepind)
        deallocate(keepind)
        if ((contour%startsaddle /= 0)) then 
            if (size(vIDs) > 0) then ! need to hedge for zero intersections
                if (contour%startsaddle /= vIDs(1)) then 
                    tvIDs = [contour%startsaddle, tvIDs]
                    ts2 = [1, ts2]
                end if 
            else
                tvIDs = [contour%startsaddle, tvIDs]
                ts2 = [1, ts2]
            end if 
        end if 
        if ((contour%endsaddle /= 0)) then 
            if (size(vIDs) > 0) then ! need to hedge for zero intersections
                if (contour%endsaddle /= vIDs(size(vIDs))) then 
                    tvIDs = [tvIDs, contour%endsaddle]
                    ts2 = [ts2, cp%ne]
                end if 
            else
                tvIDs = [tvIDs, contour%endsaddle]
                ts2 = [ts2, cp%ne]
            end if 
        end if 

        ! Extract faces
        call ExtractTopologicalFacesFromPolygon(cp, tvIDs, ts2, topomesh%vert%x, &
            topomesh%vert%y, vf1, vf2, xfda, yfda)

        ! Hedge for too small faces
        allocate(keepind(size(xfda)))
        keepind = .true. 
        do i = 1, size(xfda)
            tx = xfda(i)%Get()
            ty = yfda(i)%Get()
            dist = sum(sqrt((tx(2:size(tx)) - tx(1:size(tx)-1))**2) + &
                (ty(2:size(ty)) - ty(1:size(ty)-1))**2)
            if (dist <= disttol) then 
                ! Remove
                print *, 'local face ID: ', i 
                print *, 'vertices: ', vf1(i), vf2(i)
                print *, 'InsertTopologicalMeshContour: not adding face ' // & 
                    'with vertex indices as mentioned above as it is ' // & 
                    'smaller than distance tolerance'
                keepind(i) = .false. 
            end if 
        end do 
        
        ! Add to faces
        do i = 1, size(xfda)
            if (keepind(i)) then 
                call AddTopologicalMeshFace(topomesh, [vf1(i), vf2(i)], xfda(i), &
                    yfda(i), contourtype, contourfsID, contour%val)
            end if 
        end do 

        ! Housekeeping
        deallocate(keepind)
        
        ! Adjust existing faces
        !======================
        ! Sort intersections according to face index
        allocate(sortind(size(s1r)))
        call Sort(fID, ind=sortind)
        s2r = s2r(sortind)
        s1r = s1r(sortind)
        s1 = s1(sortind)
        s2 = s2(sortind)
        xint = xint(sortind)
        yint = yint(sortind)
        vIDs = vIDs(sortind)
        deallocate(sortind)

        ! Extract
        k = 0
        do while (k < size(fID)) 

            ! Update loop variables
            kold = k 
            k = findloc(fID, fID(kold+1), 1, back=.true.)

            ! Add start and end points as intersections if they have an 
            ! ID (and if that ID is not already present as an intersection)
            tvIDs = vIDs(kold+1:k) 
            ts1 = s1(kold+1:k)
            ts1r = s1r(kold+1:k)

            ! Sort along ts1r
            allocate(sortind(size(ts1r)))
            call Sort(ts1r, ind=sortind)
            tvIDs = tvIDs(sortind)
            ts1 = ts1(sortind)
            deallocate(sortind)

            ! Hedge for duplicate intersections
            allocate(keepind(size(tvIDs)))
            keepind = .true. 
            do j = 1, size(keepind)-1
                if ((ts1r(j+1) - ts1r(j)) == 0_R8) then 
                    keepind(j+1) = .false. 
                end if 
            end do 
            tvIDs = pack(tvIDs, keepind) ! can simply reduce here, not used afterwards
            ts1 = pack(ts1, keepind)
            deallocate(keepind)

            if ((topomesh%face%vert(fID(k), 1) /= 0) .and. (topomesh%face%vert(fID(k), 1) /= tvIDs(1))) then 
                tvIDs = [topomesh%face%vert(fID(k), 1), tvIDs]
                ts1 = [1, ts1]
            end if 
            if ((topomesh%face%vert(fID(k), 2) /= 0) .and. (topomesh%face%vert(fID(k), 2) /= tvIDs(size(tvIDs)))) then 
                tvIDs = [tvIDs, topomesh%face%vert(fID(k), 2)]
                ts1 = [ts1, topomesh%face%pol(fID(k))%ne]
            end if 

            ! Extract faces
            call ExtractTopologicalFacesFromPolygon(&
                topomesh%face%pol(fID(k)), tvIDs, ts1, topomesh%vert%x, &
                topomesh%vert%y, vf1, vf2, xfda, yfda)

            ! Add to faces
            do i = 1, size(xfda)
                if (topomesh%face%fsID(fID(k)) /= 0) then 
                    fsfval = topomesh%fsfval%Get(topomesh%face%fsID(fID(k)))
                else
                    fsfval = 0.0_R8
                end if 
                call AddTopologicalMeshFace(topomesh, [vf1(i), vf2(i)], xfda(i), &
                    yfda(i), topomesh%face%type(fID(k)), topomesh%face%fsID(fID(k)), fsfval)

                ! call WriteTopologicalMesh(topomesh, 'topomesh_temp')
            end do 

        end do 

        ! Remove adjusted faces
        allocate(delind(topomesh%face%ntot))
        delind = .false. 
        delind(fID) = .true. 
        call RemoveTopologicalMeshFaceLogical(topomesh, delind)

    end subroutine 

    ! Topological mesh face splitting
    subroutine SplitTopologicalMeshFaces(topomesh)

        ! Description
        !============
        ! This routine checks the topological mesh for the following faces:
        !
        ! - faces with same vertex indices (assumed still different face)
        ! - faces with the same start and end vertex indices (closed faces)
        !
        ! These faces are split up in resp. two and three parts in order to arrive
        ! at a conventional mesh format, where each face has a unique set of vertex
        ! indices (regardless the order) and no faces that close upon themselves
        ! exist. This is a prerequisite when mesh cells are determined. 

        ! To split up the faces, we add vertex nodes with ID -1 at the splitting
        ! points. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 

        ! Auxiliary 
        integer(I8)                             :: nncf, nndf, np, &
            ind(1:2), nfinit, tf
        integer(I8), allocatable, dimension(:)  :: faceID
        real(R8)                                :: newvcfx(1:2), newvcfy(1:2), &
            newvcff(1:2), newvdfx(1), newvdfy(1), newvdff(1), fsfval
        
        logical, allocatable, dimension(:)      :: isclosedface, &
            isduplicateface, delind

        type(RealDynamicArrayUDT)               :: xrda, yrda 

        ! Loop
        integer(I8)                             :: i, j, k

        ! Initialize
        !===========
        ! Set initial amount of faces
        nfinit = topomesh%face%ntot 

        ! Unpack for ease (only for determination of which faces to delete)
        associate(face      => topomesh%face, &
            ntot            => topomesh%face%ntot)

        ! Check faces
        !============
        ! Closed faces
        isclosedface = face%vert(:, 1) == face%vert(:, 2);

        ! Open faces
        allocate(isduplicateface(ntot), faceID(ntot))
        isduplicateface = .false. 
        faceID = 0
        do i = 1, ntot-1
            do j = i+1, ntot
                ! Skip closed faces - should be dealt with separately, even if
                ! multiple are present
                if ((.not. isclosedface(i)) .and. (.not. isclosedface(j))) then 
                    ! Check vertices
                    if (any(face%vert(i, 1) == face%vert(j, :)) .and. &
                        any(face%vert(i, 2) == face%vert(j, :))) then 
                        isduplicateface(j) = .true.
                        faceID(j) = i
                    end if 
                end if
            end do 
        end do 

        ! Count
        nncf = count(isclosedface)
        nndf = count(isduplicateface)

        ! Issue messages
        if (any(isclosedface)) then 
            print *, 'SplitTopologicalMeshFaces: ', nncf, &
                ' closed faces detected, splitting up ...'
        end if 
        if (any(isduplicateface)) then 
            print *, 'SplitTopologicalMeshFaces: ', nndf, &
            ' faces with the same vertices detected, splitting up ...'
        end if 

        ! Determine new vertices
        !=======================
        ! Closed faces
        do i = 1, ntot
            if (isclosedface(i)) then 
                
                ! Get number of points of this face
                np = face%x(i)%Size()
                if (np < 4) then  ! end points should be the same and duplicate
                    ! Issue message: we cannot split up this boundary
                    call gdErrorHandler('SplitTopologicalMeshFaces: ' // & 
                        'closed face found with only two coordinates, ' // & 
                        'cannot split up')
                end if 
                
                ! Split up the face into parts with approx. equal number of
                ! vertices
                if (np == 4) then 
                    ind = [2, 3]
                else
                    ind = [floor(real(np, R8)/3_R8), ceiling(real(2*np, R8)/3_R8)]
                end if
                
                ! Get vertex coordinates
                newvcfx = face%x(i)%Get(ind)
                newvcfy = face%y(i)%Get(ind)

                ! Get flux surface value
                if (face%fsID(i) /= 0) then 
                    fsfval = topomesh%fsfval%Get(face%fsID(i))
                else
                    fsfval = 0.0_R8
                end if 
                newvcff = fsfval
                
                ! Insert new vertices
                do j = 1, 2
                    call AddTopologicalMeshVertex(topomesh, &
                        newvcfx(j), newvcfy(j), newvcff(j), &
                        TMvertexsplitID, face%fsID(i))
                end do       
                
                ! Insert first face
                xrda = ConstructRealDynamicArray(face%x(i)%Get([(k, k = 1, ind(1))]))
                yrda = ConstructRealDynamicArray(face%y(i)%Get([(k, k = 1, ind(1))]))
                call AddTopologicalMeshFace(topomesh, &
                    [face%vert(i, 1), topomesh%vert%ntot-1], &
                    xrda, yrda, face%type(i), face%fsID(i), fsfval)

                ! Insert second face
                xrda = ConstructRealDynamicArray(face%x(i)%Get([(k, k = ind(1), ind(2))]))
                yrda = ConstructRealDynamicArray(face%y(i)%Get([(k, k = ind(1), ind(2))]))
                call AddTopologicalMeshFace(topomesh, &
                    [topomesh%vert%ntot-1, topomesh%vert%ntot], &
                    xrda, yrda, face%type(i), face%fsID(i), fsfval)

                ! Insert third face
                xrda = ConstructRealDynamicArray(face%x(i)%Get([(k, k = ind(2), face%x(i)%Size())]))
                yrda = ConstructRealDynamicArray(face%y(i)%Get([(k, k = ind(2), face%y(i)%Size())]))
                call AddTopologicalMeshFace(topomesh, &
                    [topomesh%vert%ntot, face%vert(i, 2)], &
                    xrda, yrda, face%type(i), face%fsID(i), fsfval)

            elseif (isduplicateface(i)) then 

                ! Check which face to split - preference to aligned faces
                ! instead of boundary faces (these may be contour parts
                ! etc)

                if (face%type(i) /= TMfacebndID) then 
                    tf = i
                elseif (face%type(faceID(i)) /= TMfacebndID) then 
                    tf = faceID(i)
                else
                    tf = i ! default
                end if 

                ! Get number of points of this face
                np = face%x(tf)%Size()
                if (np < 3) then  ! end points should be the same and duplicate
                    ! Issue message: we cannot split up this boundary
                    call gdErrorHandler('SplitTopologicalMeshFaces: ' // & 
                        'face with same vertices found with only ' // & 
                        'two coordinates, cannot split up')
                end if 
                
                ! Split up the face into parts with approx. equal number of
                ! vertices. 
                ind(1) = np/2+1
                
                ! Get vertex coordinates
                newvdfx = face%x(tf)%Get(ind(1))
                newvdfy = face%y(tf)%Get(ind(1))

                ! Get flux surface value
                if (face%fsID(tf) /= 0) then 
                    fsfval = topomesh%fsfval%Get(face%fsID(tf))
                else
                    fsfval = 0.0_R8
                end if 
                newvcff = fsfval
                
                ! Insert new vertex
                call AddTopologicalMeshVertex(topomesh, newvdfx(1), &
                    newvdfy(1), newvdff(1), TMvertexsplitID, face%fsID(tf))

                ! Insert first face
                xrda = ConstructRealDynamicArray(face%x(tf)%Get([(k, k = 1, ind(1))]))
                yrda = ConstructRealDynamicArray(face%y(tf)%Get([(k, k = 1, ind(1))]))
                call AddTopologicalMeshFace(topomesh, &
                    [face%vert(i, 1), topomesh%vert%ntot], &
                    xrda, yrda, face%type(tf), face%fsID(tf), fsfval)

                ! Insert second face
                xrda = ConstructRealDynamicArray(face%x(tf)%Get([(k, k = ind(1), face%x(tf)%Size())]))
                yrda = ConstructRealDynamicArray(face%y(tf)%Get([(k, k = ind(1), face%y(tf)%Size())]))
                
                call AddTopologicalMeshFace(topomesh, &
                    [topomesh%vert%ntot, face%vert(i, 2)], &
                    xrda, yrda, face%type(tf), face%fsID(tf), fsfval)

            end if 

        end do

        ! Remove old faces
        !=================
        ! Set IDs correctly
        allocate(delind(face%ntot))
        delind = .false.
        delind(1:nfinit) = isclosedface .or. isduplicateface
        call RemoveTopologicalMeshFaceLogical(topomesh, delind)

        ! Housekeeping
        !=============
        end associate

    end subroutine 

    ! Topomesh trimming
    subroutine TrimTopologicalMesh(topomesh, magneticField, vessel)

        ! Description
        !============
        ! This function removes boundaries and vertices from a topological mesh if
        ! those vertices and boundaries lie outside of the domain. To determine the
        ! latter, the vessel levelset function is checked, which is
        ! positive when points are outside of the boundary. Any points that lie on
        ! the boundary should be identified as tangency points (type 4 or 5) or
        ! other type of boundary points (type 6).

        ! All other points are removed. We assume that all intersections of faces with
        ! the boundary have been computed and that all segments are added as
        ! separate faces.

        ! Faces that do not have a start and end vertex are removed per definition.
        ! If a start or end vertex has been deleted, also the faces that have this
        ! start or end vertex will be deleted. For other faces, it is checked if
        ! the majority of points lies inside the vessel, excluding start and end
        ! vertices (they may lie exactly on the vessel or should be deleted already
        ! if they lie outside of the vessel). If all vertices (except end points)
        ! lie outside the vessel, there is no ambiguity and the face is deleted.
        ! If the first node near the start/end points are still inside
        ! the vessel, but the rest outside, it is assumed that this is due to mesh
        ! accuracy (a message will be displayed). If more than one node on each
        ! side is detected, a warning is issued and the face is not deleted. This
        ! may be due to not having computed and inserted all intersections, or due
        ! to misuse of this routine. 

        ! Notes
        !======
        ! Note 1: actually, we could also check based on the start and end vertex
        ! type whether we should consider a certain boundary for deletion (actually
        ! only boundaries with 'regular' intersections or tangency points). Now, we
        ! simply check all boundaries. 

        ! Note 2: we rely on the face identifiers to not consider boundary faces.
        ! This is necessary, since boundary faces may lie just on or off the vessel
        ! contour, depending on how accurate it was traced. Therefore, this routine
        ! is likely to fail for these boundaries (as expected...). Boundary faces
        ! are expected to be of type 3. 

        ! Note 3: separatrix segments (type 4) that do not have a 
        ! saddle point at either end are removed as well, since they
        ! are not strictly necessary in the topological mesh. Splitting
        ! points are also allowed of course. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        type(magneticFieldUDT), intent(in)      :: magneticField 
        type(VesselUDT), intent(inout)          :: vessel 

        ! Auxiliary
        real(R8), allocatable, dimension(:)     :: Vv
        logical, allocatable, dimension(:)      :: outbnd, rmvert, &
            rmface
        type(PolygonUDT), allocatable           :: bndpol(:)
        type(PolygonSetUDT)                     :: bndps 
        class(PLF2DOptionsUDT), allocatable     :: bndplfoptions 
        


        ! Loop
        integer(I8)                             :: i 

        ! Initialize
        !===========
        ! Unpack
        associate(&
            plf         => vessel%plfvessel,      &
            mfinterp    => magneticField%interp)

        ! Rebuild the vessel description to be sure
            allocate(bndpol(count(topomesh%face%type == TMfacebndID)))
            bndpol = pack(topomesh%face%pol, topomesh%face%type == TMfacebndID)
            call bndps%Construct(bndpol)
            call ConstructVesselPolygonSet(vessel, bndps)
            allocate(PLF2DClosedExactOptionsUDT::bndplfoptions)
            call InitializePolygonLevelsetFunction2D(vessel%plfvessel, &
                vessel%polygonset, bndplfoptions)

        ! Vertices
        !=========
        ! Check if in boundary
        allocate(Vv(topomesh%vert%ntot))
        call plf%Evaluate(topomesh%vert%x, topomesh%vert%y, 0, 0, Vv)
        outbnd = Vv >= 0

        ! Check if we should remove it
        rmvert = outbnd .and. (topomesh%vert%type /= TMvertextp1ID) .and. &
            (topomesh%vert%type /= TMvertextp2ID) .and. (topomesh%vert%type /= TMvertexbndID)

        ! Remove these vertices
        call RemoveTopologicalMeshVertexLogical(topomesh, rmvert)

        ! Faces
        !======
        ! Start by removing faces with zero start or end vertex
        rmface = (topomesh%face%vert(:, 1) == 0) .or. &
            (topomesh%face%vert(:, 2) == 0)

        ! Also remove separatrix faces without saddle points or splitting points
        do i = 1, topomesh%face%ntot
            if (topomesh%face%type(i) == TMfacesepID) then 
                ! Set to true, will be set to false if saddle point present
                rmface(i) = .true. 
                if (topomesh%face%vert(i, 1) /= 0) then 
                    if ((topomesh%vert%type(topomesh%face%vert(i, 1)) == TMvertexsaddleID) .or. &
                        (topomesh%vert%type(topomesh%face%vert(i, 1)) == TMvertexsplitID)) then 
                        rmface(i) = .false. 
                    end if 
                else 
                    ! Make sure is removed because of zero vertex
                    rmface(i) = .true.
                    
                end if 
                if (topomesh%face%vert(i, 2) /= 0) then 
                    if ((topomesh%vert%type(topomesh%face%vert(i, 2)) == TMvertexsaddleID) .or. &
                        (topomesh%vert%type(topomesh%face%vert(i, 2)) == TMvertexsplitID)) then 
                        rmface(i) = .false. 
                    end if 
                else 
                    ! Make sure is removed because of zero vertex
                    rmface(i) = .true.
                end if 
            end if 
        end do 

        ! Remove
        call RemoveTopologicalMeshFaceLogical(topomesh, rmface)
        call WriteTopologicalMesh(topomesh, 'topomesh_temp')
        call plf%ps%WriteData('vesselpolygon')
        ! call plf%Visualize('vesselplf', nxin=1000, nyin=1000)

        ! Check remaining faces 
        deallocate(rmface)
        allocate(rmface(topomesh%face%ntot))
        rmface = .false. 
        do i = 1, topomesh%face%ntot
            if (topomesh%face%type(i) /= TMfacebndID) then 
                ! Points in boundary (exclude end points)?
                deallocate(Vv)
                allocate(Vv(size(topomesh%face%x(i)%Get())))
                call plf%Evaluate(topomesh%face%x(i)%Get(), &
                    topomesh%face%y(i)%Get(), 0, 0, Vv)
                outbnd = Vv(2:size(Vv)-1) >= 0
                
                ! Check
                if (all(outbnd)) then 
                    ! Remove, no issue
                    rmface(i) = .true.
                elseif (all(outbnd(2:size(outbnd)-1)) .and. (size(outbnd) > 2)) then 
                    ! Remove, but display message
                    print *, 'face ID: ', i, 'face vertices: ', topomesh%face%vert(i, 1), topomesh%face%vert(i, 2)
                    print *, 'TrimTopologicalMesh: boundary removed ' // & 
                        'which still had first two points in domain'
                elseif ((.not. any(outbnd(2:size(outbnd)-1))) .and. (size(outbnd) > 2)) then 
                    ! Do not remove, nothing to do here
                else
                    ! Check if we simply had a very short boundary
                    if (size(outbnd) <= 2) then 
                        print *, 'TrimTopologicalMesh: face detected ' // & 
                            'with at most four points, may not be ' // & 
                            'properly diagnosed for removal. Try ' // & 
                            'increasing number of points by increasing ' //& 
                            'contour mesh resolution'
                    else
                        ! Points detected inside and outside - throw warning
                        print *, 'TrimTopologicalMesh: face detected ' // & 
                            'that has multiple points inside and outside ' // & 
                            'of domain. May result in unexpected outcome. Check input'

                        if (count(.not. outbnd) < count(outbnd)) then 
                            rmface(i) = .true.
                        end if 
                    end if 
                end if 
            end if 
        end do

        ! Remove
        call RemoveTopologicalMeshFaceLogical(topomesh, rmface)

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Topomesh simplification
    subroutine SimplifyTopologicalMeshFaces(topomesh)

        ! Description
        !============
        ! This routine attempts to simplify the topological mesh by 
        ! merging faces (and hence deleting vertices) if they would
        ! form a single, unique face again. To achieve this, we loop 
        ! over all vertices and check the following conditions:
        ! - does the vertex only appear in two separate boundaries?
        ! - are those two boundaries of the same type?
        ! - is the vertex a regular or boundary vertex?
        ! If all these conditions are met, the neighbouring faces may
        ! be safely merged into a single face. Normally, this merging 
        ! shouldn't be necessary, unless e.g. separatrix parts are 
        ! removed during topological mesh trimming. This routine should
        ! therefore be called after adding all boundaries to the 
        ! topological mesh, but before adding cells and other data. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 

        ! Auxiliary
        integer(I8)                             :: nforig, tf(1:2)
        integer(I8), allocatable, dimension(:)  :: rmvID, rmf1, rmf2, &
            fvert 
        logical, allocatable, dimension(:)      :: markv, markf, &
            appearstwice
        real(R8)                                :: fsfval
        real(R8), allocatable, dimension(:)     :: tempx, tempy 
        type(RealDynamicArrayUDT)               :: xda, yda 

        ! Loop 
        integer(I8)                             :: i, k 

        ! Initialize
        !===========
        ! Keep looping until exit condition reached
        do while (.true.) 

            ! Store original sizes
            nforig = topomesh%face%ntot

            ! Mark vertices and faces for removal
            allocate(markv(topomesh%vert%ntot), markf(topomesh%face%ntot), &
                rmvID(0), rmf1(0), rmf2(0))
            markv = .false. 
            markf = .false.
            do i = 1, topomesh%vert%ntot
                ! Is it a regular or boundary vertex?
                if ((.not. topomesh%vert%type(i) == TMvertexbndID) .and. &
                    (.not. topomesh%vert%type(i) == 0)) then 
                    cycle 
                end if 

                ! If it is, does it only appear twice in face%vert?
                appearstwice = ((topomesh%face%vert(:, 1) == i) .or. (topomesh%face%vert(:, 2) == i))
                if (.not. (count(appearstwice) == 2)) then 
                    cycle 
                end if 
                tf = pack([(k, k = 1, topomesh%face%ntot)], appearstwice)

                ! Are both faces not yet marked for deletion?
                if (any(markf(tf))) then 
                    cycle ! do in a next iteration
                end if 

                ! Do both boundaries have the same type and flux surface
                ! ID? 
                if (.not. (topomesh%face%type(tf(1)) == topomesh%face%type(tf(2)))) then 
                    cycle
                end if 
                if (.not. (topomesh%face%fsID(tf(1)) == topomesh%face%fsID(tf(2)))) then 
                    cycle
                end if 

                ! If we got here, we passed all checks. Mark for merging
                ! and deletion
                markv(i) = .true.
                markf(tf) = .true. 
                rmvID = [rmvID, i]
                rmf1 = [rmf1, tf(1)]
                rmf2 = [rmf2, tf(2)]

            end do 

            ! Check exit condition
            if (count(markv) == 0) then 
                exit 
            end if 

            ! Merge faces
            !============
            do i = 1, size(rmf1)
                ! Check vertices
                if (topomesh%face%vert(rmf1(i), 2) == rmvID(i)) then 
                    ! First face is in good position 
                    xda = topomesh%face%x(rmf1(i))
                    yda = topomesh%face%y(rmf1(i))
                    if (topomesh%face%vert(rmf2(i), 1) == rmvID(i)) then 
                        ! Simply append 
                        fvert = [topomesh%face%vert(rmf1(i), 1), &
                            topomesh%face%vert(rmf2(i), 2)]
                        call xda%Append(topomesh%face%x(rmf2(i))%Get())
                        call yda%Append(topomesh%face%y(rmf2(i))%Get())

                    elseif (topomesh%face%vert(rmf2(i), 2) == rmvID(i)) then 
                        ! Need to flip second face
                        fvert = [topomesh%face%vert(rmf1(i), 1), &
                            topomesh%face%vert(rmf2(i), 1)]
                        tempx = topomesh%face%x(rmf2(i))%Get()
                        tempy = topomesh%face%y(rmf2(i))%Get()
                        call xda%Append(tempx(size(tempx):1:-1))
                        call yda%Append(tempy(size(tempy):1:-1))

                    else 
                        ! This is a bug
                        call gdErrorHandler('SimplifyTopologicalMesh: ' // &
                            'bug detected')
                    end if 
                elseif (topomesh%face%vert(rmf1(i), 1) == rmvID(i)) then 
                    ! Need to flip first face
                    tempx = topomesh%face%x(rmf1(i))%Get()
                    tempy = topomesh%face%y(rmf1(i))%Get()
                    tempx = tempx(size(tempx):1:-1)
                    tempy = tempy(size(tempy):1:-1)
                    xda = ConstructRealDynamicArray(tempx)
                    yda = ConstructRealDynamicArray(tempy)
                    if (topomesh%face%vert(rmf2(i), 1) == rmvID(i)) then 
                        ! Simply append 
                        fvert = [topomesh%face%vert(rmf1(i), 2), &
                            topomesh%face%vert(rmf2(i), 2)]
                        call xda%Append(topomesh%face%x(rmf2(i))%Get())
                        call yda%Append(topomesh%face%y(rmf2(i))%Get())

                    elseif (topomesh%face%vert(rmf2(i), 2) == rmvID(i)) then 
                        ! Need to flip second face
                        fvert = [topomesh%face%vert(rmf1(i), 2), &
                            topomesh%face%vert(rmf2(i), 1)]
                        tempx = topomesh%face%x(rmf2(i))%Get()
                        tempy = topomesh%face%y(rmf2(i))%Get()
                        call xda%Append(tempx(size(tempx):1:-1))
                        call yda%Append(tempy(size(tempy):1:-1))
                        
                    else 
                        ! This is a bug
                        call gdErrorHandler('SimplifyTopologicalMesh: ' // &
                            'bug detected')
                    end if 
                else 
                    ! This is a bug
                    call gdErrorHandler('SimplifyTopologicalMesh: ' // &
                        'bug detected')
                end if 

                ! Add the new face
                if (topomesh%face%fsID(rmf1(i)) /= 0) then 
                    fsfval = topomesh%fsfval%Get(topomesh%face%fsID(rmf1(i)))
                else
                    fsfval = 0.0_R8
                end if 
                call AddTopologicalMeshFace(topomesh, fvert, &
                    xda, yda, topomesh%face%type(rmf1(i)), &
                    topomesh%face%fsID(rmf1(i)), fsfval)

            end do 

            ! Extend the face deletion vertices
            markf = [markf, spread(.false., 1, topomesh%face%ntot - nforig)]

            ! Delete the faces
            call RemoveTopologicalMeshFaceLogical(topomesh, markf)

            ! Delete the vertices
            call RemoveTopologicalMeshVertexLogical(topomesh, markv)

            ! Recompute some required data
            call AddTopologicalMeshVertexFaces(topomesh)

            ! Data
            call AddTopologicalMeshData(topomesh)

            ! Housekeeping
            deallocate(markv, markf, rmvID, rmf1, rmf2)

        end do 

    end subroutine 

    ! Garbage tangency point removal
    subroutine RemoveGarbageTangencyPoints(topomesh)

        ! Description
        !============
        ! This routine identifies 'garbage' tangency points and removes
        ! them by changing their type and applying the topomesh
        ! simplification algorithm (SimplifyTopologicalMeshFaces). 
        ! Garbage tangency points are defined as type 1 tangency points
        ! that are not limited on both sides by the same flux surface.
        ! This assumes that contours have been added to the topomesh!

        ! We now additionally remove boundary points that result in 
        ! garbage behavior related to tangency points. In particular, 
        ! we check type 2 tangency points and see if they have
        ! a face which is a boundary face that has a regular boundary
        ! point with the same flux surface ID. This should normally not
        ! happen and indicates that the tangency point contour originally
        ! intersected with the boundary, but that the intersection 
        ! removal was not successful. This happens very rarily though...

        ! Note: no additional interconnection data is updated

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh

        ! Auxiliary
        integer(I8)                             :: ntvfal, tv1, tv2
        integer(I8), allocatable, dimension(:)  :: tvf, tvfvID, tvfvIDu, &
            tvf1, tvf2
        logical, allocatable, dimension(:)      :: delf, delv

        ! Loop
        integer(I8)                             :: j, k

        ! Initialize
        !===========
        ! Associate
        associate(&
            vert        => topomesh%vert,   &
            face        => topomesh%face)

        ! Simplify to be sure
        call SimplifyTopologicalMeshFaces(topomesh)

        ! Reconstruct vertex faces
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Checks
        !=======
        ! Type 1 tangency point garbage
        !------------------------------
        ! Check which TPs should be removed
        do j = 1, vert%ntot
            if (vert%type(j) == TMvertextp1ID) then 
                ! Check
                tvf = vert%GetFace(j)

                ! Sanity check
                if (size(tvf) /= 2) then 
                    ! Write out some data for debugging
                    call WriteTopologicalMesh(topomesh, 'topomesh_error')
                    print *, 'vertex: ', j
                    call gdErrorHandler('RemoveGarbageTangencyPoints: ' // & 
                        'type 1 tangency point does not have ' // & 
                        'exactly two faces, unexpected. Check mesh in topomesh_error.dat')
                end if 

                ! Get all unique flux surface IDs of these vertices
                tvfvID = vert%fsID([face%vert(tvf, 1), face%vert(tvf, 2)])
                call Unique(tvfvID, tvfvIDu)
                if (size(tvfvIDu) > 2 .and. .not. any(tvfvIDu == 0_I8)) then 
                    ! Set vertex type to simple boundary vertex
                    ! for deletion later on in simplification 
                    ! step
                    vert%type(j) = TMvertexbndID
                    print *, 'RemoveGarbageTangencyPoints: removing vertex (type 1) ', j
                end if 
            end if 
        end do 

        ! Simplify 
        call SimplifyTopologicalMeshFaces(topomesh)

        ! Reconstruct vertex faces
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Type 2 tangency point garbage
        !------------------------------
        ! Initialize
        allocate(delv(vert%ntot), delf(face%ntot))
        delv = .false.
        delf = .false.

        ! Check which boundary tangency points should be removed
        do j = 1, vert%ntot
            if (vert%type(j) == TMvertextp2ID) then 
                ! Get the faces of this vertex
                tvf = vert%GetFace(j)

                ! Check if there are two faces with the same flux 
                ! surface ID. If not, we need to check boundary faces
                ntvfal = count(face%fsID(tvf) == vert%fsID(j))
                if (ntvfal < 2) then 
                    ! Get boundary faces, should be two
                    tvf = pack(tvf, face%type(tvf) == TMfacebndID)

                    ! Check
                    if (size(tvf) /= 2) then 
                        call WriteTopologicalMesh(topomesh, 'topomesh_error')
                        print *, 'vertex: ', j 
                        call gdErrorHandler('RemoveGarbageTangencyPoints: ' // & 
                            'type 2 tangency points does not have exactly ' // & 
                            'two boundaries, check topomesh in topomesh_error')
                    end if 

                    ! Check the flux surface ID of non-tangency point vertices
                    ! of both faces
                    if (face%vert(tvf(1), 1) == j) then 
                        tv1 = face%vert(tvf(1), 2)
                    else
                        tv1 = face%vert(tvf(1), 1)
                    end if 
                    if (face%vert(tvf(2), 1) == j) then 
                        tv2 = face%vert(tvf(2), 2)
                    else
                        tv2 = face%vert(tvf(2), 1)
                    end if 

                    ! Check first vertex
                    if (vert%fsID(tv1) == vert%fsID(j)) then 
                        ! Mark for deletion
                        delv(tv1) = .true.
                        delf(tvf(1)) = .true.

                        ! Print
                        print *, 'RemoveGarbageTangencyPoints: removing ', &
                            'vertex ', tv1

                        ! Get other faces of this vertex
                        tvf1 = vert%GetFace(tv1)
                        tvf1 = pack(tvf1, tvf1 /= tvf(1))

                        ! Adjust the vertex of these faces
                        do k = 1, size(tvf1)
                            ! Adjust the end point
                            if (face%vert(tvf1(k), 1) == tv1) then 
                                face%vert(tvf1(k), 1) = j 
                                call face%x(tvf1(k))%Set(1, vert%x(j))
                                call face%y(tvf1(k))%Set(1, vert%y(j))
                            else
                                face%vert(tvf1(k), 2) = j 
                                call face%x(tvf1(k))%Set(face%x(tvf1(k))%Size(), vert%x(j))
                                call face%y(tvf1(k))%Set(face%y(tvf1(k))%Size(), vert%y(j))
                            end if

                            ! Reconstruct the polygon
                            call face%pol(tvf1(k))%Construct(face%x(tvf1(k))%Get(), &
                                    face%y(tvf1(k))%Get())
                        end do 
                    end if

                    ! Check second vertex
                    if (vert%fsID(tv2) == vert%fsID(j)) then 
                        ! Mark for deletion
                        delv(tv2) = .true.
                        delf(tvf(2)) = .true.

                        ! Print
                        print *, 'RemoveGarbageTangencyPoints: removing ', &
                            'vertex ', tv2

                        ! Get other faces of this vertex
                        tvf2 = vert%GetFace(tv2)
                        tvf2 = pack(tvf2, tvf2 /= tvf(2))

                        ! Adjust the vertex of these faces
                        do k = 1, size(tvf2)
                            ! Adjust the end point
                            if (face%vert(tvf2(k), 1) == tv2) then 
                                face%vert(tvf2(k), 1) = j 
                                call face%x(tvf2(k))%Set(1, vert%x(j))
                                call face%y(tvf2(k))%Set(1, vert%y(j))
                            else
                                face%vert(tvf2(k), 2) = j 
                                call face%x(tvf2(k))%Set(face%x(tvf2(k))%Size(), vert%x(j))
                                call face%y(tvf2(k))%Set(face%y(tvf2(k))%Size(), vert%y(j))
                            end if

                            ! Reconstruct the polygon
                            call face%pol(tvf2(k))%Construct(face%x(tvf2(k))%Get(), &
                                    face%y(tvf2(k))%Get())
                        end do 
                    end if
                end if 
            end if 
        end do

        ! Remove
        call RemoveTopologicalMeshVertexLogical(topomesh, delv)
        call RemoveTopologicalMeshFaceLogical(topomesh, delf)


        ! Housekeeping
        !=============
        end associate

    end subroutine

    !------------------------------------------------------------------!
    !                   EQUILIBRIUM CHARACTERIZATION                   !
    !------------------------------------------------------------------!

    ! 2D extrema 
    subroutine TraceExtrema2D(xe, ye, fe, typee, fieldtracer, &
        magneticField, donewton)

        ! Description
        !============
        ! This routine traces all extrema for a 2D (magnetic) field.
        ! The routine then traces all points in the bounded domain that satisfy:
        !   
        !       dFdx = 0
        !       dFdy = 0
        !
        ! Depending on the higher order derivatives in those points, they can then
        ! be classified as 
        !       
        !       local minimum: [d2Fdx2 d2Fdxdy; d2Fdydx d2Fdy2] is positive
        !       definite
        !       saddle point: non-positive definite
        !       local maximum:  negative definite
        !
        ! This is checked by computing the eigenvalues of the local hessian.
        ! In this case, this is relatively trivial, since the matrix is 
        ! only 2D and hence the eigenvalues can be computed analytically.  

        ! Algorithm
        !==========
        ! Actually, the algorithm is rather simple since only 2 independent
        ! variables (x, y) are present. We trace contour lines of dFdx = 0 and
        ! dFdy = 0 and then look for intersections of these lines. At
        ! intersections, we call a Newton solver to compute the exact location of
        ! dFdx = 0, dFdy = 0, and we compute locally the hessian and its
        ! eigenvalues. 

        ! Notes
        !======
        ! Note 1: the contour tracing algorithm used to trace dFdx = 0 and dFdy = 0
        ! contours is assumed to be able to deal with saddle points and should
        ! return only 'simple' polygons. Saddle points may exist even in the dFdx =
        ! 0 and dFdy = 0 fields. Therefore, we've implemented our own tracing
        ! routine in TraceContourLineStructured2D. 

        ! Note 2: when computing the optimum with the Newton solver, it is assumed
        ! that no damping strategy is necessary since the initial point should be
        ! (very) close to the optimum. 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), allocatable, intent(out)      :: xe(:), ye(:), fe(:)
        integer(I8), allocatable, intent(out)   :: typee(:)
        class(ContourTracerUDT), intent(inout)  :: fieldtracer 
        type(MagneticFieldUDT), intent(in)      :: magneticField
        logical, intent(in)                     :: donewton

        ! Auxiliary
        integer(I8)                             :: ngp, nfxc, nfyc, &
            nx 
        integer(I8), allocatable, dimension(:)  :: ts1, ts2, tt
        real(R8)                                :: tempx, tempy
        real(R8), allocatable, dimension(:)     :: xg, yg, f, fx, fy, &
            tx, ty, thiseig, tf, tfxx, tfxy, tfyy
        logical                                 :: conv 
        type(RealDynamicArrayUDT)               :: xc, yc, fc 
        type(IntegerDynamicArrayUDT)            :: tc
        type(ContourUDT), allocatable           :: fxc(:), fyc(:)
        type(PolygonSetUDT)                     :: tempps
        type(PolygonUDT), allocatable           :: fxp(:), fyp(:)
        class(ContourTracerUDT), allocatable    :: fxtracer, fytracer

        ! Loop
        integer(I8)                             :: i, j, k

        ! Initialize
        !===========
        ! Check inputs
        if (.not. donewton) then 
            print *, 'TraceExtrema2DBox: not applying Newton solver, ' // & 
                'locations of extrema may be slightly inaccurate'
        end if

        ! Extract grid coordinates
        call fieldtracer%GetCoordinates(xg, yg)
        ngp = size(xg)

        ! Initialize dynamic arrays
        xc = ConstructRealDynamicArray()
        yc = ConstructRealDynamicArray()
        fc = ConstructRealDynamicArray()
        tc = ConstructIntegerDynamicArray()

        ! Trace contours
        !===============
        ! Evaluate derivatives
        allocate(f(ngp), fx(ngp), fy(ngp))
        call magneticField%interp%Evaluate(xg, yg, 0, 0, f)
        call magneticField%interp%Evaluate(xg, yg, 1, 0, fx)
        call magneticField%interp%Evaluate(xg, yg, 0, 1, fy)

        ! Construct new tracers
        fxtracer = fieldtracer
        fytracer = fieldtracer 
        if (allocated(fxtracer%xs)) then 
            deallocate(fxtracer%xs, fxtracer%ys, fxtracer%vs, &
                fxtracer%order, fxtracer%IDs)
        end if 
        allocate(fxtracer%xs(0), fxtracer%ys(0), fxtracer%vs(0), &
            fxtracer%order(0), fxtracer%IDs(0))
        if (allocated(fytracer%xs)) then 
            deallocate(fytracer%xs, fytracer%ys, fytracer%vs, &
            fytracer%order, fytracer%IDs)
        end if 
        allocate(fytracer%xs(0), fytracer%ys(0), fytracer%vs(0), &
            fytracer%order(0), fytracer%IDs(0))
        call fxtracer%SetValues(fx)
        call fytracer%SetValues(fy)
        

        ! Trace 
        fxc = fxtracer%TraceContours([0.0_R8])
        fyc = fytracer%TraceContours([0.0_R8])

        ! Clean
        call CleanContours(fxc)
        call CleanContours(fyc)

        ! Convert all to polygons
        nfxc = size(fxc)
        nfyc = size(fyc)
        allocate(fxp(nfxc), fyp(nfyc))

        do i = 1, nfxc 
            call fxp(i)%Construct(fxc(i)%x, fxc(i)%y)
        end do 
        do i = 1, nfyc 
            call fyp(i)%Construct(fyc(i)%x, fyc(i)%y)
        end do 

        ! Write out data
        call tempps%Construct(fxp)
        call tempps%WriteData('extrema_fx_lines')
        call tempps%Construct(fyp)
        call tempps%WriteData('extrema_fy_lines')
        

        ! Compute all intersections
        !==========================
        ! Compute intersections
        !$omp parallel do default(none) schedule(dynamic) collapse(2) &
        !$omp shared(nfxc, nfyc, fxc, fyc, magneticField, xc, yc, fc, tc) &
        !$omp private(i, j, tx, ty, ts1, ts2, nx, tt, k, tempx, tempy, &
        !$omp conv, tf, tfxx, tfxy, tfyy, thiseig)
        do i = 1, nfxc
            ! Only compute intersections with other polygons
            do j = 1, nfyc
                ! Compute intersections with next polygon
                call SimplePolygonIntersections(fxc(i)%x, fxc(i)%y, &
                    fyc(j)%x, fyc(j)%y, tx, ty, ts1, ts2)
                !call PolygonIntersections(fxp(i), fyp(j), tx, ty, ts1, ts2)
                
                ! Check if found
                nx = size(tx)
            
                ! Add intersections
                allocate(tt(nx))
                do k = 1, nx 
                    ! Refine
                    call TinyNewtonSolver(tempx, tempy, conv, &
                        tx(k), ty(k), magneticField)
                     
                    ! Check 
                    if (conv) then 
                        tx(k) = tempx
                        ty(k) = tempy
                    else
                        print *, 'TraceExtrema2DBox: newton solver did ' // & 
                            'not converge, taking original estimate'
                    end if
                end do 

                ! Compute value at location and second order derivatives
                allocate(tf(nx), tfxx(nx), tfxy(nx), tfyy(nx))
                call magneticField%interp%Evaluate(tx, ty, 0, 0, tf)
                call magneticField%interp%Evaluate(tx, ty, 2, 0, tfxx)
                call magneticField%interp%Evaluate(tx, ty, 1, 1, tfxy)
                call magneticField%interp%Evaluate(tx, ty, 0, 2, tfyy)

                ! Determine order
                do k = 1, nx

                    ! Compute value and check type of extremum
                    thiseig = ComputeEigenvaluesSymmetric2by2Matrix(&
                        tfxx(k), tfyy(k), tfxy(k))
                    if (any(thiseig == 0)) then 
                        print *, 'TraceExtrema2DBox: extrema with zero ' // & 
                            'eigenvalue detected, may be improperly identified'
                    end if 
                    if (all(thiseig > 0)) then 
                        ! Local minimum
                        tt(k) = -1
                    elseif (all(thiseig < 0)) then 
                        ! Local maximum
                        tt(k) = 1
                    else
                        ! Saddle point
                        tt(k) = 0
                    end if
                end do 

                ! Append
                !$omp critical
                call xc%Append(tx)
                call yc%Append(ty)
                call fc%Append(tf)
                call tc%Append(tt)
                !$omp end critical

                ! Housekeeping
                deallocate(tt, tf, tfxx, tfyy, tfxy)
            end do 
        end do
        !$omp end parallel do 
        
        ! Extract
        !========
        xe = xc%Get()
        ye = yc%Get()
        fe = fc%Get()
        typee = tc%Get()

    end subroutine 

    ! 2D tangency points
    subroutine TraceTangencyPoints2D(xtp, ytp, typetp, Ftp, &
        boundarytracer, magneticField)

        ! Description
        !============
        ! This routine traces the tangency points on a prescribed boundary.
        ! We base ourselves solely on the
        ! discrete representation of the vessel geometry and compute the tangency
        ! points by looking at the field evaluated in the vertices of the vessel
        ! polygon(s). Here, minima and maxima are quite easily and rapidly found
        ! without having to compute any intersections between curves. To determine
        ! the type of tangency point (i.e. whether the curve bends inside or
        ! outside the domain), we check whether it's a minima or maxima along the
        ! vessel curve and whether the tangency point field value is found in the
        ! interior of the domain near the tangency point location. 

        ! The latter is a bit tricky here: since we don't explicitly know the
        ! magnetic field in the interior, we evaluate the derivative in the
        ! tangency point location and base ourselves on that. This, however, may be
        ! inaccurate since the tangency point location itself will not correspond
        ! with the location of the continuous optimum. Caution is therefore advised
        ! when the local magnetic field varies strongly locally (or when the mesh
        ! is too coarse). 


        ! Notes
        !======
        ! Note 1: the contour tracing algorithm used to trace dFdx = 0 and dFdy = 0
        ! contours is assumed to be able to deal with saddle points and should
        ! return only 'simple' polygons. Saddle points may exist even in the dFdx =
        ! 0 and dFdy = 0 fields. Therefore, we've implemented our own tracing
        ! routine in TraceContourLineStructured2D. 

        ! Note 2: it is assumed that the vessel polygon forms a closed (or a set of
        ! closed) surfaces. These polygons are oriented later on such that the face
        ! normal points inward the domain at all times. 


        ! Declare variables
        !==================
        ! Arguments
        real(R8), allocatable, intent(out)      :: xtp(:), ytp(:), &
            Ftp(:)
        integer(I8), allocatable, intent(out)   :: typetp(:)
        class(ContourTracerUDT), intent(in)     :: boundarytracer
        type(MagneticFieldUDT), intent(in)      :: magneticField 

        ! Auxiliary
        integer(I8)                             :: flag
        integer(I8), allocatable                :: dvals(:), ddvals(:), &
            tv(:), extrlocind(:), tt(:)
        real(R8)                                :: fdifftol 
        real(R8), allocatable, dimension(:)     :: val, dval, tx, ty, &
            tf, nxpe, nype, nxp, nyp, normprod, dFdx, dFdy
        logical                                 :: hasbeendeleted
        logical, allocatable                    :: extrloc(:)
        type(RealDynamicArrayUDT)               :: xc, yc, fc
        type(IntegerDynamicArrayUDT)            :: tc
        type(ContourUDT), allocatable           :: bndcontours(:)
        type(PolygonUDT), allocatable           :: bndpol(:)
        type(PolygonSetUDT)                     :: bndps

        ! Loop
        integer(I8)                             :: i, j, k 

        ! Initialize
        !===========
        ! Initialize dynamic arrays
        xc = ConstructRealDynamicArray()
        yc = ConstructRealDynamicArray()
        fc = ConstructRealDynamicArray()
        tc = ConstructIntegerDynamicArray()

        ! Trace contours
        !---------------
        ! Boundary
        bndcontours = boundarytracer%TraceContours([0.0_R8])
        
        ! Sanity check
        if (size(bndcontours) == 0) then 
            ! Throw error - no boundary polygon
            call gdErrorHandler('TraceTangencyPoints2D: could not trace ' // & 
                'boundary contour, check input')
        end if 

        ! Construct boundary polygonset
        !==============================
        ! Construct polygons from boundary
        allocate(bndpol(size(bndcontours)))
        do i = 1, size(bndcontours)
            ! Construct
            call bndpol(i)%Construct(bndcontours(i)%x, bndcontours(i)%y)

            ! Check for closedness, if not -> error
            if (.not. bndpol(i)%isclosed) then 
                call gdErrorHandler('TraceTangencyPoints2D: boundary polygon ' // & 
                    'is not closed, not supported. Check input')
            end if 
        end do 

        ! Construct polygonset
        call bndps%Construct(bndpol)

        ! Orient
        call bndps%OrientNestedClosedPolygons(flag)

        ! Check if successful
        if (flag /= 0) then
            call gdErrorHandler('TraceTangencyPoints2D: could not orient ' // & 
                'boundary polygons, check input')
        end if 

        ! Compute tangency points
        !========================
        do i = 1, size(bndps%polygons)
            ! Associate for ease
            associate(p         => bndps%polygons(i))

            ! Evaluate field at vertex locations
            allocate(val(size(p%vert)))
            call magneticField%interp%Evaluate(p%x(p%vert), p%y(p%vert), &
                0, 0, val) ! assumed start and end point the same
            val = [val, val(2)] ! extend to take next edge into account

            ! Take difference
            dval = val(2:size(val)) - val(1:size(val)-1)

            ! Check where this changes sign
            allocate(dvals(size(dval)))
            where (dval > 0) dvals = 1
            where (dval <= 0) dvals = -1
            ddvals = dvals(2:size(dvals)) - dvals(1:size(dvals)-1) 

            ! Find the location and value of extrema
            extrloc = [.false., ddvals /= 0]
            allocate(tv(count(extrloc)))
            tv = pack(p%vert, extrloc)
            tx = p%x(tv)
            ty = p%y(tv)
            tf = val(tv) ! should be fine
            
            ! Check if we should exclude extremum pairs based on field value
            ! difference
            fdifftol = (maxval(val) - minval(val))*tprelfieldtol
            k = 1
            hasbeendeleted = .false.
            do while (k < size(tf))
                ! Check difference
                if (abs(tf(k+1)-tf(k)) < fdifftol) then 
                    ! Remove values, such that subsequent ones can be
                    ! checked too
                    hasbeendeleted = .true.
                    extrloc(tv(k:k+1)) = .false.
                    tv = [tv(:k-1), tv(k+1:)]
                    tf = [tf(:k-1), tf(k+1:)]
                    k = 1
                else
                    k = k + 1
                end if
            end do 
            
            ! Check last 'edge'
            if (abs(tf(size(tf))-tf(1)) < fdifftol) then 
                extrloc(tv([1, size(tf)])) = .false.
                tv = pack(p%vert, extrloc)
                tf = val(tv)
                hasbeendeleted = .true.
            end if
            
            ! Issue message
            if (hasbeendeleted) then 
                print *, 'TraceTangencyPoints2DBoxDiscrete: some ' // & 
                    'tangency points were deleted based on their field ' // & 
                    'values as they are very close together'
            end if

            ! Recompute points
            tv = pack(p%vert, extrloc)
            extrlocind = tv 
            extrlocind = pack([(k, k = 1, size(extrloc))], extrloc)
            tx = p%x(tv)
            ty = p%y(tv)
            tf = val(tv) ! should be fine

            ! Compute product between normal at vertex and magnetic field
            allocate(dFdx(size(tf)), dFdy(size(tf)))
            call magneticField%interp%Evaluate(tx, ty, 1, 0, dFdx)
            call magneticField%interp%Evaluate(tx, ty, 0, 1, dFdy)

            nxpe = [p%nx, p%nx(1)]
            nype = [p%ny, p%ny(1)]
            nxp = 0.5*(nxpe(extrlocind-1) + nxpe(extrlocind))
            nyp = 0.5*(nype(extrlocind-1) + nype(extrlocind))
            normprod = dFdx*nxp + dFdy*nyp

            ! Determine type and add
            call xc%Append(tx)
            call yc%Append(ty)
            call fc%Append(tf)
            allocate(tt(size(tf)))
            tt = 0
            do j = 1, size(extrlocind)

                ! Determine type
                if ((val(extrlocind(j)-1) < tf(j)) .and. &
                    (val(extrlocind(j)+1) < tf(j))) then 

                    ! Local maximum
                    if (normprod(j) < 0) then 
                        ! Curve bends outwards of the domain
                        tt(j) = 1
                    else
                        ! Curve bends inwards
                        tt(j) = -1
                    end if
                else
                    ! Local minimum
                    if (normprod(j) > 0) then 
                        ! Curve bends outwards of the domain
                        tt(j) = 1
                    else
                        ! Curve bends inwards
                        tt(j) = -1
                    end if 
                end if
            end do 

            ! Append
            call tc%Append(tt)

            ! Housekeeping
            deallocate(val, dvals, dFdx, dFdy, tt, tv)
            end associate

        end do 

        ! Extract values
        !===============
        xtp = xc%Get()
        ytp = yc%Get()
        Ftp = fc%Get()
        typetp = tc%Get()



    end subroutine

    !------------------------------------------------------------------!
    !                    TOPOLOGICAL MESH OPERATORS                    !
    !------------------------------------------------------------------!

    ! Initializers
    subroutine InitializeTopologicalMesh(topomesh)

        ! Description
        !============
        ! Initialize the topological mesh substructures

        ! Declare variables
        !==================
        class(TopomeshUDT)      :: topomesh 

        ! Initialize
        !===========
        ! number of distinct flux surfaces
        topomesh%nfs = 0
        topomesh%fsfval = ConstructRealDynamicArray()
        topomesh%fsID   = ConstructIntegerDynamicArray()

        ! Substructures (tubes are initialized empty)
        call topomesh%vert%Initialize()
        call topomesh%face%Initialize()
        call topomesh%cell%Initialize()
        call topomesh%tube%Initialize(0, 0, 0, 0, 0, 0, 0)

    end subroutine

    subroutine InitializeTopologicalMeshVertex(tpvert, nv)

        ! Description
        !============
        ! Initialize the topomesh vertex structure (simply empty arrays)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshVertUDT)      :: tpvert 
        integer(I8), intent(in), optional   :: nv 

        ! Initialize
        !===========
        if (present(nv)) then 
            tpvert%ntot = nv 
        else
            tpvert%ntot = 0
        end if 
        if (allocated(tpvert%ID)) then 
            ! Assume all allocated
            deallocate(tpvert%ID, tpvert%x, tpvert%y, tpvert%fval, &
                tpvert%type, tpvert%fsID, tpvert%BV)
        end if 
        associate(ntot  => tpvert%ntot)
        allocate(tpvert%ID(ntot), tpvert%x(ntot), tpvert%y(ntot), &
            tpvert%fval(ntot), tpvert%type(ntot), tpvert%fsID(ntot), &
            tpvert%BV(ntot))

        end associate

    end subroutine 

    subroutine InitializeTopologicalMeshFace(tpface, nf, nfc)

        ! Description
        !============
        ! Initialize the topomesh vertex structure (simply empty arrays)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshFaceUDT)      :: tpface
        integer(I8), intent(in), optional   :: nf, nfc

        ! Initialize
        !===========
        if (present(nf)) then 
            tpface%ntot = nf
        else
            tpface%ntot = 0
        end if
        if (present(nfc)) then 
            tpface%nfc = nfc
        else 
            tpface%nfc = 0
        end if 
        associate(ntot  => tpface%ntot, ntotc   => tpface%nfc)
        if (allocated(tpface%ID)) then 
            ! assume all allocated
            deallocate(tpface%ID, tpface%vert, tpface%cell, &
                tpface%fsID, tpface%x, tpface%y, &
                tpface%pol, tpface%type, tpface%cellP, tpface%label, &
                tpface%BF)
        end if 
        allocate(tpface%ID(ntot), tpface%vert(ntot, 2), tpface%cell(ntotc), &
            tpface%fsID(ntot), tpface%x(ntot), tpface%y(ntot), &
            tpface%pol(ntot), tpface%type(ntot), tpface%cellP(ntot, 2), &
            tpface%label(ntot), tpface%BF(ntot))

        ! Housekeeping
        end associate

    end subroutine

    subroutine InitializeTopologicalMeshCell(tpcell, nc, ncv, ncf)

        ! Description
        !============
        ! Initialize the topomesh cell structure (simply empty arrays)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshCellUDT)      :: tpcell
        integer(I8), intent(in), optional   :: nc, ncv, ncf

        ! Initialize
        !===========
        if (present(nc)) then 
            tpcell%ntot = nc
        else
            tpcell%ntot = 0
        end if 
        if (present(ncv)) then 
            tpcell%ntotv = ncv 
        else 
            tpcell%ntotv = 0
        end if 
        if (present(ncf)) then 
            tpcell%ntotf = ncf 
        else
            tpcell%ntotf = 0
        end if 

        if (allocated(tpcell%ID)) then 
            ! Assume all allocated
            deallocate(tpcell%ID, tpcell%vert, tpcell%vertP, tpcell%face,  &
                tpcell%faceP, tpcell%flags)
        end if 
        associate(ntot  => tpcell%ntot, ntotf   => tpcell%ntotf, ntotv  => tpcell%ntotv)
        allocate(tpcell%ID(ntot), tpcell%vert(ntotv), tpcell%vertP(ntot, 2), &
            tpcell%face(ntotf),  tpcell%faceP(ntot, 2), tpcell%flags(ntot))
        end associate

    end subroutine

    subroutine InitializeTopologicalMeshTube(tptube, ntot, nface, ncell, &
        nf1, nf2, nv1, nv2)

        ! Description
        !============
        ! Initialize the topomesh tube structure (simply empty arrays)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshTubeUDT)      :: tptube 
        integer(I8), intent(in)     :: ntot, nface, ncell, nv1, nv2, &
            nf1, nf2

        ! Initialize
        !===========
        tptube%ntot = ntot
        tptube%nface = nface
        tptube%ncell = ncell
        if (allocated(tptube%cell)) then 
            ! assume all allocated
            deallocate(tptube%cell, tptube%face, &
                tptube%cellP, tptube%faceP, tptube%isclosed, tptube%bndf1, tptube%bndf2, &
                tptube%bndv1, tptube%bndv2, tptube%bndf1P, &
                tptube%bndf2P, tptube%bndv1P, tptube%bndv2P)
        end if 

        allocate(tptube%cell(ncell), tptube%face(nface), &
            tptube%cellP(ntot, 2), tptube%faceP(ntot, 2), &
            tptube%isclosed(ntot), tptube%bndf1(nf1), tptube%bndf2(nf2), &
            tptube%bndv1(nv1), tptube%bndv2(nv2), tptube%bndf1P(ntot, 2), &
            tptube%bndf2P(ntot, 2), tptube%bndv1P(ntot, 2), tptube%bndv2P(ntot, 2))

    end subroutine

    ! Deallocators
    subroutine DeallocateTopologicalMeshTube(tptube)

        ! Description
        !============
        ! Soft deallocation

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshTubeUDT)          :: tptube 

        ! Deallocate
        !===========
        tptube%ntot = 0
        tptube%nface = 0
        tptube%ncell = 0
        if (allocated(tptube%face)) then 
            deallocate(tptube%face)
        end if 
        if (allocated(tptube%faceP)) then 
            deallocate(tptube%faceP)
        end if 
        if (allocated(tptube%cell)) then 
            deallocate(tptube%cell)
        end if 
        if (allocated(tptube%cellP)) then 
            deallocate(tptube%cellP)
        end if 
        if (allocated(tptube%isclosed)) then 
            deallocate(tptube%isclosed)
        end if 
        if (allocated(tptube%bndf1)) then 
            deallocate(tptube%bndf1)
        end if 
        if (allocated(tptube%bndf2)) then 
            deallocate(tptube%bndf2)
        end if 
        if (allocated(tptube%bndv1)) then 
            deallocate(tptube%bndv1)
        end if 
        if (allocated(tptube%bndv2)) then 
            deallocate(tptube%bndv2)
        end if 
        if (allocated(tptube%bndf1P)) then 
            deallocate(tptube%bndf1P)
        end if 
        if (allocated(tptube%bndf2P)) then 
            deallocate(tptube%bndf2P)
        end if 
        if (allocated(tptube%bndv1P)) then 
            deallocate(tptube%bndv1P)
        end if 
        if (allocated(tptube%bndv2P)) then 
            deallocate(tptube%bndv2P)
        end if 

    end subroutine

    ! Vertex addition
    subroutine AddTopologicalMeshVertex(topomesh, x, y, F, t, fsID)

        ! Description
        !============
        ! Add a vertex to the topological mesh. Not optimized for 
        ! memory usage (e.g. less allocate/deallocate), but 
        ! shouldn't be an issue since topomeshes are usually small. Note
        ! that no interconnection data is updated!

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                  :: topomesh 
        real(R8), intent(in)                :: x, y, F
        integer(I8), intent(in)             :: t, fsID

        ! Auxiliary
        real(R8)                            :: tF 

        ! Check 
        !======
        ! If the flux surface ID was already added, overwrite the flux
        ! value - this should be the same (but discretely this often
        ! doesn't hold) 
        tF = F 
        if (fsID /= 0) then 
            if ((topomesh%fsID%Size() > 0)) then 
                if (any(topomesh%fsID%Get() == fsID)) then 
                    ! Already exists, overwrite value
                    tF = topomesh%fsfval%Get(fsID)
                else
                    ! Does not yet exist, add 
                    call topomesh%fsID%Append(fsID)
                    call topomesh%fsfval%Set(fsID, F)
                end if
            else
                ! Does not yet exist, add 
                call topomesh%fsID%Append(fsID)
                call topomesh%fsfval%Set(fsID, F)
            end if 
        end if 


        ! Concatenate
        !============
        topomesh%vert%ntot = topomesh%vert%ntot + 1
        topomesh%vert%x = [topomesh%vert%x, x]
        topomesh%vert%y = [topomesh%vert%y, y]
        topomesh%vert%fval = [topomesh%vert%fval, tF]
        topomesh%vert%type = [topomesh%vert%type, t]
        topomesh%vert%ID = [topomesh%vert%ID, topomesh%vert%ntot]
        topomesh%vert%fsID = [topomesh%vert%fsID, fsID]

    end subroutine

    ! Vertex face addition
    subroutine AddTopologicalMeshVertexFaces(topomesh)

        ! Description
        !============
        ! This function adds the faces of each vertex and sorts them in a uniquely
        ! defined direction which is the same for all faces. It is assumed that all
        ! faces and vertices of the topology are adequately defined and that each
        ! face has a starting and end point etc (so basically after all
        ! intersections and topology faces have been added). 

        ! Algorithm
        !==========
        ! 1) For each vertex, find all faces that have this vertex and store them
        ! in vert.face and vert.faceP (list and pointer)
        ! 2) For each vertex, sort this list by doing the following steps:
        !   2.1) For each face, get the closest point of that face to the current
        !   vertex (but with distance > 1e-12 to hedge for numerical bullshit)
        !   2.2) Take one face as reference, compute the angle of all other faces
        !   w.r.t. that first face. 
        !   2.3) Sort the faces in sequence of increasing angle
        !   2.4) Check the cross product of each pair of consecutive faces. All
        !   cross-products should have the same sign. If it is positive, keep it
        !   like that. If it is negative, reverse the sorting direction. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: nfpv, tf, sortind, &
            tfsorted , fID
        real(R8)                                :: vx, vy
        real(R8), allocatable, dimension(:)     :: npx, npy, dx, dy,&
            theta
        logical, allocatable                    :: test(:)

        ! Loop
        integer(I8)                             :: i, j, k

        ! Initialize
        !===========
        ! Do some checks
        if (any(topomesh%face%vert(:, 1) == 0) .or. any(topomesh%face%vert(:, 2) == 0)) then  
            call gdErrorHandler('AddTopologicalMeshVertexFaces: faces without vertex ' // & 
                'encountered. First, clean up the topological mesh before calling this routine')
        end if 

        ! Get faces of each vertex
        !-========================
        ! Compute for each vertex how many faces they would have
        allocate(nfpv(topomesh%vert%ntot))
        nfpv = 0
        do i = 1, topomesh%face%ntot
            nfpv(topomesh%face%vert(i, :)) = nfpv(topomesh%face%vert(i, :)) + 1
        end do

        ! Initialize 
        if (allocated(topomesh%vert%face)) then 
            deallocate(topomesh%vert%face)
        end if 
        if (allocated(topomesh%vert%faceP)) then 
            deallocate(topomesh%vert%faceP)
        end if 
        allocate(topomesh%vert%face(sum(nfpv)))

        ! Construct the face pointer
        allocate(topomesh%vert%faceP(topomesh%vert%ntot, 2))
        topomesh%vert%faceP(:, 2) = nfpv 
        topomesh%vert%faceP(1, 1) = 1
        do i = 2, topomesh%vert%ntot 
            topomesh%vert%faceP(i, 1) = topomesh%vert%faceP(i-1, 1) + nfpv(i-1)
        end do 

        ! Get face list (unsorted)
        allocate(fID(topomesh%face%ntot))
        fID = [(k, k = 1, topomesh%face%ntot)]
        do i = 1, topomesh%vert%ntot
            ! Get face indices
            test = (topomesh%face%vert(:, 1) == i) .or. (topomesh%face%vert(:, 2) == i) 
            tf = pack(fID, test)
            
            ! Add
            topomesh%vert%face(topomesh%vert%faceP(i, 1):topomesh%vert%faceP(i, 1)+topomesh%vert%faceP(i, 2)-1) = tf
        end do

        ! Sort faces
        !===========
        ! Sort
        do i = 1, topomesh%vert%ntot
            ! Unpack this vertex
            vx = topomesh%vert%x(i)
            vy = topomesh%vert%y(i)

            ! Get faces
            tf = GetTMVertFace(topomesh%vert, i)
            
            ! Get the next point of each face 
            allocate(npx(size(tf)), npy(size(tf)))
            do j = 1, size(tf)
                if (topomesh%face%vert(tf(j), 1) == i) then  ! first vertex
                    npx(j) = topomesh%face%x(tf(j))%Get(2)
                    npy(j) = topomesh%face%y(tf(j))%Get(2)
                elseif (topomesh%face%vert(tf(j), 2) == i) then  ! second vertex
                    npx(j) = topomesh%face%x(tf(j))%Get(topomesh%face%x(tf(j))%size()-1)
                    npy(j) = topomesh%face%y(tf(j))%Get(topomesh%face%y(tf(j))%size()-1)
                else
                    call gdErrorHandler('AddTopologicalMeshVertexFaces: ' // & 
                        'This should not be happening and is a bug!')
                end if 
            end do 
            
            ! Compute angles (w.r.t. horizontal axis)
            dx = npx - vx
            dy = npy - vy
            theta = atan2(dy, dx)
            
            ! Sort
            allocate(sortind(size(theta)))
            call Sort(theta, ind=sortind)
            tfsorted = tf(sortind)
            deallocate(sortind)
            
            ! Add
            topomesh%vert%face(topomesh%vert%faceP(i, 1):&
                topomesh%vert%faceP(i, 1)+topomesh%vert%faceP(i, 2)-1) = tfsorted

            ! Housekeeping
            deallocate(npx, npy)
            
        end do


    end subroutine

    ! Face addition
    subroutine AddTopologicalMeshFace(topomesh, facevert, x, y, &
        t, fsID, F)

        ! Description
        !============
        ! Add a face to the topological mesh. Not optimized for 
        ! memory usage (e.g. less allocate/deallocate), but 
        ! shouldn't be an issue since topomeshes are usually small. Note
        ! that no interconnection data is updated!

        ! Note: to hedge for any issues in other routines, it is checked
        ! whether subsequent coordinates are coinciding to disttol 
        ! precision. These are removed from the x, y coordinates. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                  :: topomesh 
        type(RealDynamicArrayUDT)           :: x, y 
        integer(I8), intent(in)             :: t, facevert(1:2), fsID
        real(R8), intent(in)                :: F 

        ! Auxiliary
        integer(I8), allocatable            :: tmp(:, :)
        real(R8)                            :: tF 
        real(R8), allocatable               :: xt(:), yt(:), dxt(:), &
            dyt(:)
        logical, allocatable                :: delind(:)
        type(PolygonUDT)                    :: tp 

        ! Check if flux surface exists
        !=============================
        tF = F 
        if (fsID /= 0) then 
            if ((topomesh%fsID%Size() > 0)) then 
                if (any(topomesh%fsID%Get() == fsID)) then 
                    ! Already exists, overwrite value
                    tF = topomesh%fsfval%Get(fsID)
                else
                    ! Does not yet exist, add 
                    call topomesh%fsID%Append(fsID)
                    call topomesh%fsfval%Set(fsID, F)
                end if
            else
                ! Does not yet exist, add 
                call topomesh%fsID%Append(fsID)
                call topomesh%fsfval%Set(fsID, F)
            end if 
        end if 

        ! Hedge for duplicate points
        !===========================
        ! Detect duplicate points
        xt = x%Get()
        yt = y%Get()
        dxt = xt(2:size(xt)) - xt(1:size(xt)-1)
        dyt = yt(2:size(yt)) - yt(1:size(yt)-1)
        allocate(delind(size(xt)))
        delind = .false. 
        delind(2:size(delind)) = ((abs(dxt) <= disttol) .and. (abs(dyt)) <= disttol)
        if (any(delind)) then 
            x = ConstructRealDynamicArray(pack(xt, .not. delind))
            y = ConstructRealDynamicArray(pack(yt, .not. delind))
        end if 

        ! Make sure that final and start point are equal to vertex 
        ! coordinates if facevert is not zero
        if (x%size() < 2) then 
            print *, 'size is smaller than 2'
        end if 
        if (facevert(1) /= 0) then 
            call x%Set(1, topomesh%vert%x(facevert(1)))
            call y%Set(1, topomesh%vert%y(facevert(1)))
        end if 
        if (facevert(2) /= 0) then 
            call x%Set(x%Size(), topomesh%vert%x(facevert(2)))
            call y%Set(y%Size(), topomesh%vert%y(facevert(2)))
        end if 

        !  Construct new polygon
        call tp%Construct(x%Get(), y%Get())

        ! Concatenate
        !============
        topomesh%face%ntot = topomesh%face%ntot + 1
        tmp = topomesh%face%vert 
        deallocate(topomesh%face%vert)
        allocate(topomesh%face%vert(topomesh%face%ntot, 2))
        topomesh%face%vert(1:topomesh%face%ntot-1, :) = tmp 
        topomesh%face%vert(topomesh%face%ntot, :) = facevert 
        topomesh%face%x = [topomesh%face%x, x]
        topomesh%face%y = [topomesh%face%y, y]
        topomesh%face%fsID = [topomesh%face%fsID, fsID]
        topomesh%face%type = [topomesh%face%type, t]
        topomesh%face%ID = [topomesh%face%ID, topomesh%face%ntot]
        topomesh%face%pol = [topomesh%face%pol, tp]

    end subroutine

    ! Face cell addition
    subroutine AddTopologicalMeshFaceCells(topomesh)

        ! Description
        !============
        ! Add the cell neighbours of a face using the pointer way of working. It is
        ! assumed that all cells etc have been constructed correctly and that other
        ! basic topology information is available. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh   

        ! Auxiliary
        integer(I8)                             :: ind
        integer(I8), allocatable, dimension(:)  :: ncpf, tf, fc

        ! Loop
        integer(I8)                             :: i, j

        ! Add faces
        !==========
        ! Initialize
        if (allocated(topomesh%face%cellP)) then 
            deallocate(topomesh%face%cellP)
        end if 
        if (allocated(topomesh%face%cell)) then 
            deallocate(topomesh%face%cell)
        end if 
        allocate(topomesh%face%cellP(topomesh%face%ntot, 2))

        ! Compute number of cells per face
        allocate(ncpf(topomesh%face%ntot))
        ncpf = 0
        do i = 1, topomesh%cell%ntot
            ! Get cell faces
            tf = GetTMCellFace(topomesh%cell, i)
            ncpf(tf) = ncpf(tf) + 1
        end do

        ! Set cellP
        topomesh%face%cellP(:, 2) = ncpf
        topomesh%face%cellP(1, 1) = 1
        do i = 2, topomesh%face%ntot 
            topomesh%face%cellP(i, 1) = topomesh%face%cellP(i-1, 1) + & 
                topomesh%face%cellP(i-1, 2)
        end do 

        ! Set cells of each face
        allocate(topomesh%face%cell(sum(ncpf)), fc(topomesh%face%ntot))
        fc = 0
        topomesh%face%cell = 0
        do i = 1, topomesh%cell%ntot
            ! Get cell faces
            tf = GetTMCellFace(topomesh%cell, i);
            do j = 1, size(tf)
                ! Add
                ind = topomesh%face%cellP(tf(j), 1) + fc(tf(j))
                topomesh%face%cell(ind) = i
                
                ! Update counter
                fc(tf(j)) = fc(tf(j)) + 1
            end do 
        end do 

    end subroutine 

    ! Cell addition 
    subroutine AddTopologicalMeshCells(topomesh)

        ! Description
        !============
        ! This routine forms cells of the given topology mesh with consistent data
        ! of faces and vertices. It is assumed that all faces are 'simple', i.e.
        ! they do not start and end in the same point (in that case, points
        ! should've been introduced before). The cells can be of arbitrary topology
        ! (i.e. triangles, quads, ...) - if this is not desired, it should be
        ! checked afterwards.

        ! The following cell types are allowed: 
        ! - regular cells: cells with all unique boundaries that do not appear
        ! twice for this cell 
        ! - Disc-type cells: cells that have a single boundary appearing twice, but
        ! where this single boundary has one vertex that only has one face.

        ! All other cell types will not be properly identified by this routine and
        ! will very likely lead to undesired output or even crashing. If necessary,
        ! some other cell types may be supported in the future by making the stop
        ! condition for a cell more stringent (i.e. not just stopping if we find
        ! the same face, but if we also arrive at the face in the same way as we
        ! started - this would allow cells with duplicate boundaries as well). 

        ! Algorithm
        !==========
        ! We assume that the following data is readily available:
        ! - face indices of vertices, sorted either clockwise or counter-clockwise,
        ! but consistently the same for each vertex. This allows us to uniquely
        ! define the left and right neighbor face when determining the next one.
        ! This should be computed beforehand with the
        ! 'AddTopologicalMeshVertexFaces' routine.
        ! - Face identifiers that are accurate and up to date: at least boundary
        ! faces must be indicated to be type 3. No closed faces are assumed to
        ! exist. 

        ! Furthermore, we assume that the ensemble of vertices and faces results in
        ! a non-overlapping partition of the bounded 2D domain, and that the union
        ! of all faces with type 3 represents a closed boundary polygon (or
        ! multiple closed polygons) that bound the domain.

        ! Declare variables
        !==================
        ! Arguments 
        class(TopomeshUDT)                      :: topomesh 

        ! Auxiliary 
        integer(I8)                             :: tf, startface, & 
            turndirection, tfv(1:2), tvind, tv, startvert, &
            starttvind, nf, nfv(1:2)
        integer(I8), allocatable, dimension(:)  :: fc, disccellvert, &
            nfvfn, tcf, tcv, faceneig1, faceneig2, faceneig

        logical                                 :: istfv(1:2)
        logical, allocatable                    :: hasturned1(:, :), &
            hasturned2(:, :), donotstartfromface(:)

        type(IntegerDynamicArrayUDT), allocatable   :: cellvert(:), &
            cellface(:)
        type(IntegerDynamicArrayUDT)                :: thiscellvert, &
            thiscellface

        ! Loop 
        integer(I8)                             :: i, cc 

        ! Initialize
        !===========
        ! Unpack
        associate( & 
            vert        => topomesh%vert,   &
            face        => topomesh%face,   &
            cell        => topomesh%cell    &
            )

        ! Allocate
        allocate(fc(face%ntot), hasturned1(face%ntot, 2), &
            hasturned2(face%ntot, 2), donotstartfromface(face%ntot), &
            cellvert(0), cellface(0))

        ! Initialize face counters
        where (face%BF) 
            fc = 1 ! boundary faces: only once
        elsewhere 
            fc = 2
        end where 
            
        ! Initialize turn checkers
        hasturned1 = .false.
        hasturned2 = .false.

        ! Initialize cell counter
        cc = 0

        ! Check if there are any disc-type cells
        disccellvert = findloc(vert%faceP(:, 2), 1_I8)

        ! Check if there are any faces with only one adjacent face on each side -
        ! these faces shouldn't be started from, as one cannot determine the
        ! turning direction (and several other faces should remain that can be
        ! started from)
        donotstartfromface = .false. 
        do i = 1, face%ntot
            nfvfn = vert%faceP(face%vert(i, :), 2)
            if (all(nfvfn == 2)) then 
                donotstartfromface(i) = .true.
            end if 
        end do 

        ! Loop
        !=====
        do while (.true.)
            
            ! Initialize cell faces & vertices
            allocate(tcv(0), tcf(0))

            ! Find the next face (any next internal face)
            tf = findloc( (fc > 0) .and. (.not. face%BF) .and. &
                (.not. donotstartfromface), .true., 1)
            
            ! Check
            if (tf == 0) then 
                ! Check
                if (any(fc > 0)) then
                    ! Check if there are any internal faces to begin with
                    if (count(.not. face%BF) == 0) then 
                        print *, 'AddTopologicalMeshCells: case without any ' // & 
                            'internal faces detected'

                        ! Find the next face (any next face)
                        tf = findloc( (fc > 0)  , .true., 1)
                    else
                        call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                            'could not find next face')
                    end if 
                else
                    ! All faces added, exit
                    exit 
                end if 
            end if 
            
            ! Do not subtract a counter - we need to end up in this face again.
            ! Also, don't add, we do this later on
            
            ! Set starting face for this cell
            startface = tf
            
            ! Turning direction (0: not found, 1: first neighbour, 2: second
            ! neighbour)
            turndirection = 0
            
            ! Get neighbouring faces in correct order
            tfv = face%vert(tf, :)
            faceneig1 = GetTMVertFaceNeig(vert, tfv(1), tf);
            faceneig2 = GetTMVertFaceNeig(vert, tfv(2), tf);
            
            ! Sanity checks
            if ((size(faceneig1) == 0) .or. (size(faceneig2) == 0)) then 
                ! No neighbours found
                call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                    'could not find neighbouring faces, something wrong ' //& 
                    'with topological mesh construction. Check input')
            end if 
            if (all(faceneig1 == tf) .and. all(faceneig2 == tf)) then 
                ! Isolated face
                call gdErrorHandler('AddTopologicalMeshCells: isolated ' // & 
                    'face found, check input')
            end if
            
            ! Set the starting vertex index (if boundary vertices, needs checks)
            tvind = 0
            if ((.not. all(faceneig1 == tf)) .and. (.not. vert%BV(tfv(1)))) then 
                tvind = 1
            elseif ((.not. all(faceneig2 == tf)) .and. (.not. vert%BV(tfv(2)))) then 
                tvind = 2
            end if 

            ! If none found,  check if the first vertex has neighbours 
            ! with available faces
            if (any(fc(faceneig1) > 0) .and. (tvind == 0)) then 
                ! Check if both neighours are the same - in that case we can
                ! safely take this vertex as next vertex
                if (faceneig1(1) == tf) then 
                    ! Vertex has single face here, do not take it
                elseif (faceneig1(1) == faceneig1(2)) then 
                    tvind = 1
                else
                    ! Check if any face neighbour can be taken 
                    if (any((fc(faceneig1) > 0))) then 
                        tvind = 1
                    end if 
                end if 
            end if 

            ! If none found, check if the second vertex has neighbours 
            ! with available faces
            if (any(fc(faceneig2) > 0) .and. (tvind == 0)) then 
                ! Check if both neighours are the same - in that case we can
                ! safely take this vertex as next vertex
                if (faceneig2(1) == tf) then 
                    ! Vertex has single face here, do not take it
                elseif (faceneig2(1) == faceneig2(2)) then 
                    tvind = 2
                else
                    ! Check if any face neighbour can be taken
                    if (any((fc(faceneig2) > 0))) then  
                        tvind = 2
                    end if 
                end if 
            end if 
            
            ! Sanity check
            if (tvind == 0) then 
                ! No starting vertex found - possibly dangling face?
                call gdErrorHandler('AddTopologicalMeshCells: could not ' // & 
                    'find next cell although initial face was found - ' // & 
                    'possible dangling face detected, check input')
            end if
            
            ! Set current vertex
            tv = face%vert(tf, tvind)
            startvert = tv
            starttvind = tvind
        
            ! Loop 
            !=====
            do while (.true.)
                
                ! Get neighbouring faces in correct order
                faceneig = GetTMVertFaceNeig(vert, tv, tf)
                        
                ! Sanity checks
                if (size(faceneig) == 0) then 
                    ! No neighbours found
                    call gdErrorHandler('AddTopologicalMeshCells: could ' // & 
                        'not find neighbouring faces, something wrong ' // & 
                        'with topological mesh construction. Check input')
                end if 
                if (all(fc(faceneig) <= 0)) then 
                    ! No neighbours with counter left
                    call gdErrorHandler('AddTopologicalMeshCells: all ' // & 
                        'neighbouring faces cannot be taken anymore, ' // &
                        'faces do not seem to form cell')
                end if 
                
                ! Find the next face
                if (turndirection /= 0) then 
                    ! We have a turn direction, so we can only check if we should
                    ! throw errors
                    nf = faceneig(turndirection)
                    
                    ! Check counter
                    if (fc(nf) <= 0) then 
                        call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                            'next face is forced by turning direction ' // & 
                            'but is not available')
                    end if 
                    
                    ! If we passed this, we should check the vertices
                    nfv = face%vert(nf, :)
                    istfv = nfv == tv
                    if (.not.any(istfv)) then 
                        ! Current vertex is not found in the next face, this should
                        ! not be possible
                        call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                            'next face does not have current vertex, ' // & 
                            'check input')
                    end if 
                    if (all(istfv)) then 
                        ! Next face is a face that starts and ends in the same
                        ! vertex - not supported
                        call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                            'face detected with same start and end vertex, ' // & 
                            'not supported')
                    end if 
                    
                    ! Add and update
                    tcf = [tcf, nf]
                    tcv = [tcv, tv]
                    tf = nf
                    if (.not. istfv(1)) then 
                        tv = nfv(1)
                    else
                        tv = nfv(2)
                    end if 
                    
                    ! Update counter
                    fc(tf) = fc(tf) - 1
                    
                    ! Set that this direction can't be turned in anymore
                    ! from neither side for this face
                    if ((startface == nf) .and. (tv == startvert)) then 
                        ! Just break, turning direction etc already adjusted
                        ! before.
                        exit
                    else
                        if (turndirection == 1) then 
                            if (face%vert(tf, 1) == tv) then 
                                hasturned1(tf, 1) = .true. 
                                hasturned2(tf, 2) = .true. 
                            else 
                                hasturned1(tf, 2) = .true. 
                                hasturned2(tf, 1) = .true. 
                            end if 
                        elseif (turndirection == 2) then 
                            if (face%vert(tf, 1) == tv) then 
                                hasturned1(tf, 2) = .true. 
                                hasturned2(tf, 1) = .true. 
                            else 
                                hasturned1(tf, 1) = .true. 
                                hasturned2(tf, 2) = .true. 
                            end if 
                        else
                            call gdErrorHandler('AddTopologicalMEshCells: ' // & 
                                'bug detecetd when adjusting turning direction')
                        end if 
                    end if 
                    
                else
                    ! We don't have a turning direction yet. Check the current
                    ! neighbours
                    if (faceneig(1) == faceneig(2)) then 
                        ! We have to take this face if we can, or we should error.
                        ! Turning direction cannot be determined
                        if (fc(faceneig(1)) <= 0) then 
                            call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                                'next face is only possible face but ' // & 
                                'cannot be taken due to counter being zero')
                        end if             
                        
                        ! Take face
                        nf = faceneig(1)
                        
                        ! If we passed this, we should check the vertices
                        nfv = face%vert(nf, :)
                        istfv = nfv == tv
                        if (.not. any(istfv)) then 
                            ! Current vertex is not found in the next face, this should
                            ! not be possible
                            call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                                'next face does not have current vertex, ' // & 
                                'check input')
                        end if 
                        if (all(istfv)) then 
                            ! Next face is a face that starts and ends in the same
                            ! vertex - not supported
                            call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                                'face detected with same start and end ' // & 
                                'vertex, not supported')
                        end if

                        ! Add and update
                        tcf = [tcf, nf]
                        tcv = [tcv, tv]
                        tf = nf
                        if (.not. istfv(1)) then 
                            tv = nfv(1)
                        else
                            tv = nfv(2)
                        end if 
                        
                        ! Update counter
                        fc(tf) = fc(tf) - 1
                        
                        ! Is the next face the start face? If so, exit
                        if ((startface == nf) .and. (tv == startvert)) then 
                            exit 
                        end if 
                    else
                        ! Check which one has a positive counter (at least one,
                        ! cause we already passed a check for this)
                        if ((fc(faceneig(1)) > 0) .and. .not. hasturned1(startface, starttvind)) then 
                            ! Take face
                            nf = faceneig(1)
                            
                            ! If we passed this, we should check the vertices
                            nfv = face%vert(nf, :)
                            istfv = nfv == tv
                            if (.not. any(istfv)) then 
                                ! Current vertex is not found in the next face, this should
                                ! not be possible
                                call gdErrorHandler('AddTopologicalMeshCells:  ' // & 
                                    'next face does not have current vertex,  ' // & 
                                    'check input')
                            end if
                            if (all(istfv)) then 
                                ! Next face is a face that starts and ends in the same
                                ! vertex - not supported
                                call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                                    'face detected with same start and ' // & 
                                    'end vertex, not supported')
                            end if 
                            
                            ! Add and update
                            tcf = [tcf, nf]
                            tcv = [tcv, tv]
                            tf = nf
                            if (.not. istfv(1)) then 
                                tv = nfv(1)
                            else
                                tv = nfv(2)
                            end if 
                            
                            ! Update counter
                            fc(tf) = fc(tf) - 1
                            
                            ! Set turn direction
                            turndirection = 1
                            
                            ! Set that this direction can't be turned in anymore
                            ! from neither side
                            hasturned1(startface, starttvind) = .true.
                            if (starttvind == 1) then 
                                hasturned2(startface, 2) = .true.
                            else
                                hasturned2(startface, 1) = .true.
                            end if 
                            
                            ! Is the next face the start face? If so, exit
                            if ((startface == nf) .and. (tv == startvert)) then 
                                exit 
                            end if 
                            
                        elseif ((fc(faceneig(2)) > 0) .and. .not. hasturned2(startface, starttvind)) then 
                            
                            ! Take face
                            nf = faceneig(2)
                            
                            ! If we passed this, we should check the vertices
                            nfv = face%vert(nf, :)
                            istfv = nfv == tv;
                            if (.not. any(istfv)) then 
                                ! Current vertex is not found in the next face, this should
                                ! not be possible
                                call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                                    'next face does not have current vertex, ' // & 
                                    'check input')
                            end if
                            if (all(istfv)) then 
                                ! Next face is a face that starts and ends in the same
                                ! vertex - not supported
                                call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                                    'face detected with same start and ' // & 
                                    'end vertex, not supported')
                            end if
                            
                            ! Add and update
                            tcf = [tcf, nf]
                            tcv = [tcv, tv]
                            tf = nf
                            if (.not. istfv(1)) then 
                                tv = nfv(1)
                            else
                                tv = nfv(2)
                            end if 
                            
                            ! Update counter
                            fc(tf) = fc(tf) - 1
                            
                            ! Set turn direction
                            turndirection = 2
                            
                            ! Set that this direction can't be turned in anymore
                            ! from neither side
                            hasturned2(startface, starttvind) = .true.
                            if (starttvind == 1) then 
                                hasturned1(startface, 2) = .true.
                            else
                                hasturned1(startface, 1) = .true.
                            end if 
                            
                            ! Is the next face the start face? If so, exit
                            if ((startface == nf) .and. (tv == startvert)) then 
                                exit 
                            end if 
                            
                        else
                            call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                                'could not find next face since no ' // &
                                'combination of available start direction ' // &
                                'and free face found')
                        end if 
                    end if 
                end if 
            end do 
            
            ! Add found cell to the structure
            cc = cc + 1
            thiscellvert = ConstructIntegerDynamicArray(tcv)
            thiscellface = ConstructIntegerDynamicArray(tcf)
            cellvert = [cellvert, thiscellvert]
            cellface = [cellface, thiscellface]

            ! Housekeeping
            deallocate(tcv, tcf)
        end do 

        ! Add to the topology mesh
        if (allocated(topomesh%cell%vert)) then 
            deallocate(topomesh%cell%vert)
        end if 
        if (allocated(topomesh%cell%vertP)) then 
            deallocate(topomesh%cell%vertP)
        end if 
        if (allocated(topomesh%cell%face)) then 
            deallocate(topomesh%cell%face)
        end if 
        if (allocated(topomesh%cell%faceP)) then 
            deallocate(topomesh%cell%faceP)
        end if 

        topomesh%cell%ntot = cc
        topomesh%cell%ID = [(i, i = 1, cc)]
        if(allocated(topomesh%cell%flags)) then 
            deallocate(topomesh%cell%flags)
        end if 
        allocate(topomesh%cell%flags(cc))
        topomesh%cell%flags = 0_I8
        if (allocated(topomesh%cell%tube)) then 
            deallocate(topomesh%cell%tube)
        end if
        allocate(topomesh%cell%tube(cc))
        topomesh%cell%tube = 0_I8
        allocate(topomesh%cell%vertP(cc, 2))
        allocate(topomesh%cell%faceP(cc, 2))
        topomesh%cell%vertP(1, 1) = 1
        topomesh%cell%faceP(1, 1) = 1
        do i = 1, cc
            topomesh%cell%vertP(i, 2) = cellvert(i)%size()
            topomesh%cell%faceP(i, 2) = cellface(i)%size()
        end do 
        allocate(topomesh%cell%vert(sum(topomesh%cell%vertP(:, 2))))
        allocate(topomesh%cell%face(sum(topomesh%cell%faceP(:, 2))))
        do i = 2, cc
            topomesh%cell%vertP(i, 1) = topomesh%cell%vertP(i-1, 1) + &
                topomesh%cell%vertP(i-1, 2)
            topomesh%cell%faceP(i, 1) = topomesh%cell%faceP(i-1, 1) + &
                topomesh%cell%faceP(i-1, 2)
        end do 
        do i = 1, cc
            topomesh%cell%vert(&
                topomesh%cell%vertP(i, 1):topomesh%cell%vertP(i, 1)+topomesh%cell%vertP(i, 2)-1) = & 
                cellvert(i)%Get()
            topomesh%cell%face(&
                topomesh%cell%faceP(i, 1):topomesh%cell%faceP(i, 1)+topomesh%cell%faceP(i, 2)-1) = & 
                cellface(i)%Get()
        end do 
        topomesh%cell%ntotv = size(topomesh%cell%vert)
        topomesh%cell%ntotf = size(topomesh%cell%face)

        ! Housekeeping
        !=============
        end associate
        
    end subroutine 

    ! Data addition
    subroutine AddTopologicalMeshData(topomesh)

        ! Description
        !============
        ! Add data of the vertices, faces, and cells (if they are present) based on
        ! their type. The following fields are updated:
        ! - vert:
        !   * BV, IV: logical (nv-by-1) for boundary or internal vertex. Boundary
        !   vertices are those vertices with type 4, 5 or 6, all the rest are
        !   internal. (i.e. tangency points and regular boundary points)
        ! - face:
        !   * BF, IF: logical (nf-by-1) for boundary or internal face. Boundary
        !   faces are those faces with type 3, all the rest are internal (unless
        !   cells have already been computed and we can check the amount of 
        !   cells per face)

        ! Note: it is assumed that all necessary contours/faces have 
        ! been added already
        
        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 

        ! Loop
        integer(I8)                             :: i 

        ! Vertices
        !=========
        ! Logicals
        topomesh%vert%BV = (topomesh%vert%type == TMvertextp1ID) .or. &
            (topomesh%vert%type == TMvertextp2ID) .or. &
            (topomesh%vert%type == TMvertexbndID)

        ! Faces
        !======
        ! Logicals
        topomesh%face%BF = (topomesh%face%type == TMfacebndID)
        

        ! Check if cells are present
        if (allocated(topomesh%face%cellP)) then 
            if (size(topomesh%face%cellP, 1) == topomesh%face%ntot) then 
                topomesh%face%BF = topomesh%face%cellP(:, 2) == 1
            end if 
        end if 

        ! Labels
        if (allocated(topomesh%face%label)) then 
            deallocate(topomesh%face%label)
        end if 
        allocate(topomesh%face%label(topomesh%face%ntot))
        do i = 1, topomesh%face%ntot
            topomesh%face%label(i) = i
        end do 

    end subroutine

    ! Flux tube addition
    subroutine AddTopologicalMeshTubes(topomesh, types)

        ! Description
        !============
        ! This routine attempts to extract 'tubes' of cells (e.g. flux tubes for
        ! fusion applications). These tubes are sequences of cells that have
        ! exactly two faces of a certain topological type (or set of types). Tubes
        ! are either closed (starting and ending at the same cell) or open
        ! (starting and ending at boundary faces). If tubes can't be constructed,
        ! the output is the empty struct. If some cells remain (i.e. some tubes
        ! could not be constructed), a message is shown. 

        ! Note: 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        integer(I8), intent(in)                 :: types(:)

        ! Auxiliary
        integer(I8)                             :: thistf, thisc, &
            nextc, tc, ne, si1, si2, cpID
        integer(I8), allocatable, dimension(:)  :: tf, thesecells, &
            tfc, thiscf, tv, tfbnd, sortindex, sortedbndfaces, tf1, tf2, &
            tv1, tv2, polygonID
        integer(I8), allocatable, dimension(:, :)   :: bndfacevert, tfv
        logical                                 :: isclosed, addtube, &
            singlefacetube
        logical, allocatable, dimension(:)      :: hasfacetype, &
            hasfaces, iscellnotfound, ispolygonstart, isbranchingpolygon
        type(IntegerDynamicArrayUDT)            :: tubef, tubec, &
            ntubef, ntubec, isclosedtube, alltubef, alltubec, &
            tubebndf, alltubebndf1, alltubebndf2, ntubebndf1, ntubebndf2, &
            alltubebndv1, alltubebndv2, ntubebndv1, ntubebndv2, &
            temptf1, temptf2

        ! Loop
        integer(I8)                             :: i, k, ct

        ! Initialize
        !===========
        ! Unpack
        associate(&
            face            => topomesh%face,   &
            cell            => topomesh%cell    &
            )


        ! Get faces & cells
        !==================
        ! Check faces
        allocate(hasfacetype(face%ntot))
        hasfacetype = .false. 
        do i = 1, size(types)
            hasfacetype = hasfacetype .or. (face%type == types(i));
        end do

        ! Check cells
        allocate(hasfaces(cell%ntot))
        hasfaces = .false. 
        do i = 1, cell%ntot
            ! Get face cells
            tf = GetTMCellFace(cell, i)
            if (count(hasfacetype(tf)) == 2) then 
                hasfaces(i) = .true.
            end if 
        end do 

        ! Display warning if some cells were not found
        if (any(.not. hasfaces)) then 
            allocate(thesecells(count(.not. hasfaces)))
            thesecells = pack([(k, k = 1, cell%ntot)], .not. hasfaces)
            print *, 'The following cells do not have exactly two ' // & 
                'tube faces and will not be part of a tube: ',  thesecells
            print *, 'Vertices of these cells:'
            do i = 1, size(thesecells)
                tv = GetTMCellVert(topomesh%cell, thesecells(i))
                print *, 'vertices: ', tv 
            end do 
            deallocate(thesecells)
        end if 

        ! Construct tubes
        !================
        ! Initialize
        iscellnotfound = hasfaces
        ct = 0
        ntubef          = ConstructIntegerDynamicArray()
        ntubec          = ConstructIntegerDynamicArray()
        alltubef        = ConstructIntegerDynamicArray()
        alltubebndf1    = ConstructIntegerDynamicArray()
        alltubebndf2    = ConstructIntegerDynamicArray()
        alltubebndv1    = ConstructIntegerDynamicArray()
        alltubebndv2    = ConstructIntegerDynamicArray()
        ntubebndv1      = ConstructIntegerDynamicArray()
        ntubebndv2      = ConstructIntegerDynamicArray()
        ntubebndf1      = ConstructIntegerDynamicArray()
        ntubebndf2      = ConstructIntegerDynamicArray()
        alltubec        = ConstructIntegerDynamicArray()
        isclosedtube    = ConstructIntegerDynamicArray()

        ! Loop
        do while (.true.)
            
            ! Take the first cell
            tc = findloc(iscellnotfound, .true., 1)
            if (tc == 0) then 
                ! all found, exit
                exit 
            end if 
            iscellnotfound(tc) = .false.
            
            ! Initialize a new tube
            tubef           = ConstructIntegerDynamicArray()
            tubec           = ConstructIntegerDynamicArray()
            tubebndf        = ConstructIntegerDynamicArray()
            
            isclosed        = .false. 
            singlefacetube  = .false.
            addtube         = .true.

            ! Append the cell
            call tubec%Append(tc)
            
            ! Get the faces
            tf = GetTMCellFace(cell, tc)

            ! Get the bounding faces
            allocate(tfbnd(count(.not. hasfacetype(tf))))
            tfbnd = pack(tf, .not. hasfacetype(tf))
            call tubebndf%Append(tfbnd)
            deallocate(tfbnd)
            
            ! Get only faces that have the type
            tf = pack(tf, hasfacetype(tf))
            
            ! Sanity check
            if (size(tf) /= 2) then 
                call gdErrorHandler('AddTopologicalMeshTubes: did not ' // & 
                    'find two faces for cell which should have two ' // & 
                    'faces - this is a bug')
            end if 

            ! Check if both faces are the same - then don't trace, single cell tube
            if (tf(1) == tf(2)) then 
                isclosed = .true.
                singlefacetube = .true.
            end if 
            
            ! Take the first face and trace
            thistf = tf(1)
            thisc = tc
            
            ! Add if it is a single cell tube
            if (singlefacetube) then 
                call tubef%Append(tf)
            end if 
            
            ! Loop until all flux tube cells found from one side
            do while (.true. .and. (.not. singlefacetube)) 
                ! Add the current face before already found faces (to keep correct
                ! order)
                call tubef%Insert(thistf, 1)
                
                ! Is the face a boundary face? 
                if (face%BF(thistf)) then 
                    ! Exit
                    exit 
                end if 
                
                ! Get face cells (should be exactly two)
                tfc = GetTMFaceCell(face, thistf)
                if (size(tfc) /= 2 ) then 
                    call gdErrorHandler('Did not find two cells for an ' // & 
                        'inner face - this is a bug')
                end if 
                
                ! Check which cell is next
                if (tfc(1) /= thisc) then 
                    nextc = tfc(1)
                elseif (tfc(2) /= thisc) then 
                    nextc = tfc(2) 
                else 
                    call gdErrorHandler('Did not find next cell - this ' // &
                        'is a bug')
                end if
                
                ! Check if this cell is marked, otherwise we need to stop here and
                ! not add the flux tube
                if (((.not. hasfaces(nextc)) .or. (.not. iscellnotfound(nextc))) &
                    .and. (nextc /= tc))   then 
                    print *, 'AddTopologicalMeshTubes: could not extract a tube'
                    addtube = .false.
                    exit 
                end if 
                
                ! If we got here, we can safely add the cell and find the next face
                thisc = nextc
                if (nextc /= tc) then 
                    call tubec%Insert(thisc, 1)
                end if 
                iscellnotfound(thisc) = .false.
                
                ! Get cell faces
                thiscf = GetTMCellFace(cell, thisc)

                ! Add bounding faces
                allocate(tfbnd(count(.not. hasfacetype(thiscf))))
                tfbnd = pack(thiscf, .not. hasfacetype(thiscf))
                call tubebndf%Append(tfbnd)
                deallocate(tfbnd)

                ! Find next face
                thiscf = pack(thiscf, hasfacetype(thiscf) .and. thiscf /= thistf)
                
                ! Sanity check
                if (size(thiscf) > 1) then 
                    call gdErrorHandler('AddTopologicalMeshTubes: too ' // & 
                        'many faces found when trying to find next face')
                elseif (size(thiscf) < 1) then 
                    call gdErrorHandler('AddTopologicalMeshTubes: could ' // & 
                        'not find next face')
                end if 
                thistf = thiscf(1)
                
                ! Is the face the same face as we started from? (closed tube)
                if (thistf == tf(1)) then 
                    isclosed = .true.
                    call tubef%Append(tf(2)) ! add at end

                    ! Exit
                    exit 
                end if 
            end do 

            ! Check if we should cycle
            if (.not. addtube) then 
                cycle 
            end if 
            
            ! Take the second face and trace, unless we can't add the tube or the
            ! tube is closed
            thistf = tf(2)
            thisc = tc
            do while ((.true.) .and. (.not. isclosed))  
                ! Add the current face after already found faces (to keep correct
                ! order)
                call tubef%Append(thistf)

                ! Is the face a boundary face? 
                if (face%BF(thistf)) then 
                    ! Exit
                    exit 
                end if 

                ! Get face cells (should be exactly two)
                tfc = GetTMFaceCell(face, thistf)
                if (size(tfc) /= 2 ) then 
                    call gdErrorHandler('Did not find two cells for an ' // & 
                        'inner face - this is a bug')
                end if 
                
                ! Check which cell is next
                if (tfc(1) /= thisc) then 
                    nextc = tfc(1)
                elseif (tfc(2) /= thisc) then 
                    nextc = tfc(2) 
                else 
                    call gdErrorHandler('Did not find next cell - this ' // &
                        'is a bug')
                end if
                
                ! Check if this cell is marked, otherwise we need to stop here and
                ! not add the flux tube
                if (((.not. hasfaces(nextc)) .or. (.not. iscellnotfound(nextc))) &
                    .and. (nextc /= tc))   then 
                    print *, 'AddTopologicalMeshTubes: could not extract a tube'
                    addtube = .false.
                    exit 
                end if 
                
                ! If we got here, we can safely add the cell and find the next face
                thisc = nextc
                if (nextc /= tc) then 
                    call tubec%Append(thisc)
                end if 
                iscellnotfound(thisc) = .false.
                
                ! Get cell faces
                thiscf = GetTMCellFace(cell, thisc)

                ! Add bounding faces
                allocate(tfbnd(count(.not. hasfacetype(thiscf))))
                tfbnd = pack(thiscf, .not. hasfacetype(thiscf))
                call tubebndf%Append(tfbnd)
                deallocate(tfbnd)

                ! Find the next face
                thiscf = pack(thiscf, hasfacetype(thiscf) .and. thiscf /= thistf)
                
                ! Sanity check
                if (size(thiscf) > 1) then 
                    call gdErrorHandler('AddTopologicalMeshTubes: too ' // & 
                        'many faces found when trying to find next face')
                elseif (size(thiscf) < 1) then 
                    call gdErrorHandler('AddTopologicalMeshTubes: could ' // & 
                        'not find next face')
                end if 
                thistf = thiscf(1)
                
            end do 
            
            ! Add data to the tube if possible
            if (addtube) then 
                ! Update counter
                ct = ct + 1

                ! Add number of faces and cells
                call ntubef%Append(tubef%Size())
                call ntubec%Append(tubec%Size())
                if (isclosed) then 
                    call isclosedtube%Append(1_I8)
                else 
                    call isclosedtube%Append(0_I8)
                end if 

                ! Add faces to total faces etc
                call alltubef%Append(tubef%Get())
                call alltubec%Append(tubec%Get())

                ! Get 'polygons' from bounding faces
                ne = tubebndf%Size() 
                bndfacevert = face%vert(tubebndf%Get(), :)
                allocate(sortindex(ne), ispolygonstart(ne), isbranchingpolygon(ne))
                call SortPolygonEdges(bndfacevert, ne, sortindex, &
                    ispolygonstart, isbranchingpolygon, polygonID)
                bndfacevert(:, 1) = bndfacevert(sortindex, 1)
                bndfacevert(:, 2) = bndfacevert(sortindex, 2)
                sortedbndfaces = tubebndf%Get(sortindex)
                
                
                ! Checks
                if (any(isbranchingpolygon) .and. ne > 0) then 
                    ! Branching polygons, need to see which ones belong
                    ! together (max 2)
                    if (maxval(polygonID) > 2) then 
                        call gdErrorHandler('AddTopologicalMeshTubes: ' // & 
                            'tube has more than two branching polygons as ' // & 
                            'boundary, unexpected')
                    end if

                    ! Extract polygon faces
                    cpID = 0
                    temptf1 = ConstructIntegerDynamicArray()
                    temptf2 = ConstructIntegerDynamicArray()
                    do i = 1, ne
                        ! Get the polygon ID if it is a start
                        if (ispolygonstart(i)) then 
                            cpID = polygonID(i)
                        end if 

                        ! Add face
                        if (cpID == 1_I8) then 
                            call temptf1%Append(sortedbndfaces(i))
                        elseif (cpID == 2_I8) then 
                            call temptf2%Append(sortedbndfaces(i))
                        else
                            ! We shouldn't get here
                            call gdErrorHandler('Unexpected error')
                        end if 
                    end do 

                    ! At least the first polygon should exist
                    call alltubebndf1%Append(temptf1%Get())
                    call ntubebndf1%Append(temptf1%Size())

                    ! Add first set
                    tf1 = temptf1%Get()

                    ! Extract vertices of first set - unsorted here...
                    allocate(tv1(size(tf1)+1))
                    bndfacevert = face%vert(tf1, :)
                    call Unique([bndfacevert(:, 1), bndfacevert(:, 2)], tv1)
                    call alltubebndv1%Append(tv1)
                    call ntubebndv1%Append(size(tv1))
                    deallocate(tv1)

                    
                    ! Check if the second set exists
                    if (temptf2%Size() > 0) then 
                        ! Extract faces
                        call alltubebndf2%Append(temptf2%Get())
                        call ntubebndf2%Append(temptf2%Size())

                        ! Add second set
                        tf2 = temptf2%Get()

                        ! Extract vertices of second set 
                        allocate(tv2(size(tf2)+1))
                        bndfacevert = face%vert(tf2, :)
                        call Unique([bndfacevert(:, 1), bndfacevert(:, 2)], tv2)
                        call alltubebndv2%Append(tv2)
                        call ntubebndv2%Append(size(tv2))
                        deallocate(tv2)
                    else
                        ! Append a zero
                        call ntubebndf2%Append(0_I8)

                        ! Extract vertices
                        call ntubebndv2%Append(1_I8)
                        tf = tubef%Get()
                        tfv = face%vert(tf, :) 
                        if (all(tfv(1, 1) == tfv(:, 1) .or. tfv(1, 1) == tfv(:, 2))) then 
                            call alltubebndv2%Append(tfv(1, 1))
                        elseif (all(tfv(1, 2) == tfv(:, 1) .or. tfv(1, 2) == tfv(:, 2))) then
                            call alltubebndv2%Append(tfv(1, 2))
                        else
                            ! This shouldn't happen
                            call gdErrorHandler('AddTopologicalMeshTubes: ' // &
                                'found only one set of faces, but could not ' // & 
                                'find a common point of all tube faces as ' // & 
                                'other boundary side. Unexpected.')
                        end if 
                    end if 

                else 
                    ! Simple polygons, do sanity checks
                    if ((count(ispolygonstart) < 1) .and. (ne > 0)) then 
                        ! May be supported in the future
                        call gdErrorHandler('AddTopologicalMeshTubes: ' // & 
                            'tube does not have at least one aligned face, ' // & 
                            'not supported.')
                    elseif ((count(ispolygonstart) > 2) .and. (ne > 0)) then 
                        ! Weird 
                        call gdErrorHandler('AddTopologicalMeshTubes: ' // & 
                            'tube has more than two aligned face boundaries, ' // & 
                            'unexpected.')
                    end if 

                    ! Extract boundary faces and vertices
                    if (count(ispolygonstart) == 2) then 
                        ! Get starting indices (normally first one is always true)
                        si1 = 1
                        si2 = findloc(ispolygonstart, .true., 1, back=.true.)

                        ! Add first set
                        tf1 = sortedbndfaces(si1:si2-1)
                        call alltubebndf1%Append(tf1)
                        call ntubebndf1%Append(size(tf1))

                        ! Add second set
                        tf2 = sortedbndfaces(si2:ne)
                        call alltubebndf2%Append(tf2)
                        call ntubebndf2%Append(size(tf2))

                        ! Extract vertices of first set 
                        allocate(tv1(size(tf1)+1))
                        call ExtractPolygonVertices(bndfacevert(1:si2-1, 1:2), &
                            size(tf1), tv1)
                        call alltubebndv1%Append(tv1)
                        call ntubebndv1%Append(size(tv1))
                        deallocate(tv1)

                        ! Extract vertices of second set 
                        allocate(tv2(size(tf2)+1))
                        call ExtractPolygonVertices(bndfacevert(si2:ne, 1:2), &
                            size(tf2), tv2)
                        call alltubebndv2%Append(tv2)
                        call ntubebndv2%Append(size(tv2))
                        deallocate(tv2)

                    elseif (count(ispolygonstart) == 1) then 
                        ! Only one segment, need to search for tangency point
                        ! vertex (should be one vertex in common with all 
                        ! tube radial faces)

                        ! Add first set
                        tf1 = sortedbndfaces
                        call alltubebndf1%Append(tf1)
                        call ntubebndf1%Append(size(tf1))

                        ! 'Add' second set
                        call ntubebndf2%Append(0_I8)

                        ! Extract vertices of first set 
                        allocate(tv1(size(tf1)+1))
                        call ExtractPolygonVertices(bndfacevert, &
                            size(tf1), tv1)
                        call alltubebndv1%Append(tv1)
                        call ntubebndv1%Append(size(tv1))
                        deallocate(tv1)

                        ! Extract vertices of second set 
                        call ntubebndv2%Append(1_I8)
                        tf = tubef%Get()
                        tfv = face%vert(tf, :) 
                        if (all(tfv(1, 1) == tfv(:, 1) .or. tfv(1, 1) == tfv(:, 2))) then 
                            call alltubebndv2%Append(tfv(1, 1))
                        elseif (all(tfv(1, 2) == tfv(:, 1) .or. tfv(1, 2) == tfv(:, 2))) then
                            call alltubebndv2%Append(tfv(1, 2))
                        else
                            ! This shouldn't happen
                            call gdErrorHandler('AddTopologicalMeshTubes: ' // &
                                'found only one set of faces, but could not ' // & 
                                'find a common point of all tube faces as ' // & 
                                'other boundary side. Unexpected.')
                        end if 
                    end if 
                end if 

                ! Housekeeping
                deallocate(sortindex, ispolygonstart, isbranchingpolygon)
            end if 
            
        end do 

        ! Construct
        !==========
        ! Ensure deallocation
        call topomesh%tube%Deallocate()

        ! Reallocate
        call topomesh%tube%Initialize(ct, alltubef%Size(), alltubec%Size(), &
            alltubebndf1%Size(), alltubebndf2%Size(), alltubebndv1%Size(), &
            alltubebndv2%Size())

        ! Set fields
        topomesh%tube%face      = alltubef%Get()
        topomesh%tube%cell      = alltubec%Get()
        topomesh%tube%bndf1     = alltubebndf1%Get()
        topomesh%tube%bndf2     = alltubebndf2%Get()
        topomesh%tube%bndv1     = alltubebndv1%Get()
        topomesh%tube%bndv2     = alltubebndv2%Get()
        topomesh%tube%isclosed  = isclosedtube%Get() == 1_I8

        topomesh%tube%faceP(:, 2)   = ntubef%Get()
        topomesh%tube%cellp(:, 2)   = ntubec%Get()
        topomesh%tube%bndf1P(:, 2)  = ntubebndf1%Get()
        topomesh%tube%bndf2P(:, 2)  = ntubebndf2%Get()
        topomesh%tube%bndv1P(:, 2)  = ntubebndv1%Get()
        topomesh%tube%bndv2P(:, 2)  = ntubebndv2%Get()
        if (topomesh%tube%ntot > 0) then 
            topomesh%tube%faceP(1, 1)   = 1
            topomesh%tube%cellP(1, 1)   = 1
            topomesh%tube%bndf1P(1, 1)  = 1
            topomesh%tube%bndf2P(1, 1)  = 1
            topomesh%tube%bndv1P(1, 1)  = 1
            topomesh%tube%bndv2P(1, 1)  = 1
        end if 

        do i = 2, topomesh%tube%ntot 
            topomesh%tube%faceP(i, 1) = topomesh%tube%faceP(i-1, 1) + &     
                topomesh%tube%faceP(i-1, 2)
            topomesh%tube%cellP(i, 1) = topomesh%tube%cellP(i-1, 1) + &     
                topomesh%tube%cellP(i-1, 2)    
            topomesh%tube%bndf1P(i, 1) = topomesh%tube%bndf1P(i-1, 1) + &     
                topomesh%tube%bndf1P(i-1, 2)
            topomesh%tube%bndf2P(i, 1) = topomesh%tube%bndf2P(i-1, 1) + &     
                topomesh%tube%bndf2P(i-1, 2)    
            topomesh%tube%bndv1P(i, 1) = topomesh%tube%bndv1P(i-1, 1) + &     
                topomesh%tube%bndv1P(i-1, 2)
            topomesh%tube%bndv2P(i, 1) = topomesh%tube%bndv2P(i-1, 1) + &     
                topomesh%tube%bndv2P(i-1, 2)    
        end do 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Flux tube data addition
    subroutine AddTopologicalMeshTubeData(topomesh)

        ! Description
        !============
        ! This routine adds additional data to the topological mesh 
        ! tubes, such as which cells are in a tube (cell%tube) and which
        ! neighbours a tube has

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT), intent(inout)       :: topomesh 

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: ttc, tf1, tf2, tc, &
            tf1u, tf2u

        type(IntegerDynamicArrayUDT)            :: ftneig1, ftneig2

        ! Loop
        integer(I8)                             :: i, j, k, nnb1, nnb2

        ! Initialize
        !===========
        ! Associate for ease
        associate(&
            tube        => topomesh%tube, &
            cell        => topomesh%cell, &
            face        => topomesh%face, &
            vert        => topomesh%vert)

        ! Allocate 
        if (allocated(tube%ftneig1)) then 
            deallocate(tube%ftneig1)
        end if 
        if (allocated(tube%ftneig2)) then 
            deallocate(tube%ftneig2)
        end if 
        if (allocated(tube%ftneig1P)) then 
            deallocate(tube%ftneig1P)
        end if 
        if (allocated(tube%ftneig2P)) then 
            deallocate(tube%ftneig2P)
        end if 
        if (allocated(cell%tube)) then 
            deallocate(cell%tube)
        end if 
        allocate(cell%tube(cell%ntot))
        allocate(tube%ftneig1P(tube%ntot, 2), tube%ftneig2P(tube%ntot, 2))

        ! Initialize
        ftneig1     = ConstructIntegerDynamicArray()
        ftneig2     = ConstructIntegerDynamicArray()
        cell%tube       = 0
        tube%ftneig1P = 0
        tube%ftneig2P = 0
        if (tube%ntot > 1) then 
            tube%ftneig1P(1, 1) = 1_I8
            tube%ftneig2P(1, 1) = 1_I8
        end if  
        
        ! Cell tubes
        !===========
        do i = 1, tube%ntot
            ! Get tube cells
            ttc = tube%GetCell(i)

            ! Check
            if (any(cell%tube(ttc) /= 0)) then 
                print *, 'AddTopologicalMeshTubeData: some cells belong ' // &
                    'to multiple tubes, unexpected. Overwriting and continuing...'
            end if 

            ! Set cell tube
            cell%tube(ttc) = i
        end do 

        ! Tube neighbours
        !================
        ! Initialize
        nnb1 = 0
        nnb2 = 0

        ! Loop
        do i = 1, tube%ntot
            ! Get tube faces from first side
            tf1 = tube%GetBndFace(i, 1_I8)
            if (size(tf1) > 0) then 
                ! Loop over all face
                do j = 1, size(tf1)
                    ! Get the face cells
                    tc = face%GetCell(tf1(j))

                    ! Check if any cell tube IDs are not equal to the 
                    ! current ID. If so, add neighbour
                    do k = 1, size(tc)
                        if (cell%tube(tc(k)) /= i) then 
                            call ftneig1%Append(cell%tube(tc(k)))
                            tube%ftneig1P(i, 2) = tube%ftneig1P(i, 2) + 1
                        end if
                    end do
                end do
            end if 

            ! Get tube faces from second side
            tf2 = tube%GetBndFace(i, 2_I8)
            if (size(tf2) > 0) then 
                ! Loop over all face
                do j = 1, size(tf2)
                    ! Get the face cells
                    tc = face%GetCell(tf2(j))

                    ! Check if any cell tube IDs are not equal to the 
                    ! current ID. If so, add neighbour
                    do k = 1, size(tc)
                        if (cell%tube(tc(k)) /= i) then 
                            call ftneig2%Append(cell%tube(tc(k)))
                            tube%ftneig2P(i, 2) = tube%ftneig2P(i, 2) + 1
                        end if
                    end do
                end do
            end if 
        end do 

        ! Build (preliminary)
        tube%ftneig1 = ftneig1%Get()
        tube%ftneig2 = ftneig2%Get()
        do i = 2, tube%ntot
            tube%ftneig1P(i, 1) = tube%ftneig1P(i-1, 1) + tube%ftneig1P(i-1, 2)
            tube%ftneig2P(i, 1) = tube%ftneig2P(i-1, 1) + tube%ftneig2P(i-1, 2)
        end do 

        ! Reduce to unique set for each tube
        ftneig1 = ConstructIntegerDynamicArray()
        ftneig2 = ConstructIntegerDynamicArray()
        do i = 1, tube%ntot
            ! Get tube neighbours (this is possible now)
            tf1 = tube%GetNeig(i, 1_I8)
            tf2 = tube%GetNeig(i, 2_I8)

            ! Determine unique set
            call Unique(tf1, tf1u)
            call Unique(tf2, tf2u)

            ! Append
            call ftneig1%Append(tf1u)
            call ftneig2%Append(tf2u)

            ! Set pointer counter, but only of this tube
            tube%ftneig1P(i, 2) = size(tf1u)
            tube%ftneig2P(i, 2) = size(tf2u)
        end do 

        ! Build (final)
        tube%ftneig1 = ftneig1%Get()
        tube%ftneig2 = ftneig2%Get()
        do i = 2, tube%ntot
            tube%ftneig1P(i, 1) = tube%ftneig1P(i-1, 1) + tube%ftneig1P(i-1, 2)
            tube%ftneig2P(i, 1) = tube%ftneig2P(i-1, 1) + tube%ftneig2P(i-1, 2)
        end do 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Interconnection
    subroutine AddTopologicalMeshInterconnectionData(topomesh)

        ! Description
        !============
        ! Add additional interconnection data to the topological mesh. 
        ! Similar to AddGridInterconnections. Some fields may be 
        ! constructed again that already existed (e.g. cell faces etc), 
        ! but we keep it in here to have a general interconnections 
        ! routine.

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 

        ! Auxiliary
        integer(I8), allocatable                :: tempvcells(:, :)
        integer(I8), allocatable, dimension(:)  :: ncpv, vcount, tv

        ! Loop
        integer(I8)                             :: i, j

        ! Add face cells
        !===============
        ! Because separate routine
        call AddTopologicalMeshFaceCells(topomesh)
         
        ! Unpack & initialize
        !====================
        ! Data structures
        associate(&
            nv      => topomesh%vert%ntot,  &
            nf      => topomesh%face%ntot,  &
            nc      => topomesh%cell%ntot,  &
            v       => topomesh%vert,   &
            f       => topomesh%face,   &
            c       => topomesh%cell    &
            )
        
        ! Checks
        if (any((f%vert(:, 1) == 0) .or. (f%vert(:, 2) == 0))) then 
            call gdErrorHandler('AddTopologicalMeshInterconnectionData : ' // & 
                'Some vertex indices are zero in faces, not supported')
        end if 
        if (allocated(v%cell)) then 
            deallocate(v%cell)
        end if 
        if (allocated(v%cellP)) then 
            deallocate(v%cellP)
        end if 
        
        ! Basic interconnections
        !=======================
        ! Initialize
        allocate(ncpv(v%ntot))
        ncpv = 0
        do i = 1, size(c%vert)
            ncpv(c%vert(i)) = ncpv(c%vert(i))+1
        end do
        allocate(v%cell(sum(ncpv)), v%cellP(v%ntot, 2))
        v%cellP(1, 1) = 1
        v%cellP(:, 2) = ncpv
        do i = 2, v%ntot
            v%cellP(i, 1) = v%cellP(i-1, 1) + v%cellP(i-1, 2)
        end do 
        
        ! Sanity check
        if (any(v%cellP(:, 2) < 1)) then 
            call gdErrorHandler('AddTopologicalMeshInterconnectionData: ' // & 
                'vertex without any cells detected. Check grid interconnectivity')
        end if 
        
        ! Construct vertex cells (cells of faces and faces of cells 
        ! already found in constructor phase)
        ! Note: we construct first
        ! temporary arrays (nv-by-ncpv, nf-by-2) that are afterwards converted to
        ! cell and cellP arrays. 
        allocate(tempvcells(v%ntot, maxval(ncpv)), vcount(v%ntot))
        tempvcells = 0
        vcount = 0
        do i = 1, nc
            ! Get vertices of cell
            tv = GetTMCellVert(c, i) ! there may be doubles in here!
            
            ! Add vertex cells
            do j = 1, size(tv)
                ! Update counter
                vcount(tv(j)) = vcount(tv(j))+1

                ! Add cell
                tempvcells(tv(j), vcount(tv(j))) = i 
            end do 
        end do 
        
        ! Construct 
        do i = 1, v%ntot
            v%cell(v%cellP(i, 1):v%cellP(i, 1)+v%cellP(i, 2)-1) = & 
                tempvcells(i, 1:vcount(i))
        end do 
        
        ! Add logicals
        f%BF = f%cellP(:, 2) == 1

        ! Reconstruct boundary vertices based on boundary faces
        if (allocated(v%BV)) then 
            deallocate(v%BV)
        end if
        allocate(v%BV(v%ntot))
        v%BV = .false. 
        do i = 1, f%ntot
            if (f%BF(i)) then 
                v%BV(f%vert(i, :)) = .true. 
            end if 
        end do 

        ! Reconstruct fsID of vertices to be compliant to faces (note: 
        ! type 1 tangency points get their own flux surface ID)

                    

        ! Tubes
        !======
        ! Flux tubes 
        call AddTopologicalMeshTubes(topomesh, [TMfaceradID, TMfacebndID])

        ! Additional tube interconnection data
        call AddTopologicalMeshTubeData(topomesh)

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Vertex removal
    subroutine RemoveTopologicalMeshVertexLogical(topomesh, rmvert)

        ! Description
        !============
        ! Remove the vertices with IDs specified in 'IDs' from the topological
        ! mesh. Only vertex and face data is updated. Faces are not immediately
        ! removed, but the vertices of each face are updated (set to zero if
        ! vertex deleted, otherwise updated to new vertices). IDs is a 
        ! logical array of size 1:topomesh%vert%ntot

        ! Notes
        !======
        ! Note 1: Faces do not necessarily require start and end point yet (they
        ! can have zeros in face.vert)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        logical, intent(in)                     :: rmvert(:)

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: diffID, tv
        logical, allocatable, dimension(:)      :: keepvert

        ! Loop 
        integer(I8)                             :: i, j

        ! Initialize
        !===========
        ! Check
        if (size(rmvert) /= topomesh%vert%ntot) then 
            call gdErrorHandler('RemoveTopologicalMeshVertexLogical: ' // & 
                'rmvert has illegal size')
        end if 

        ! Initialize
        keepvert = .not. rmvert 

        ! Associate
        associate(&
            nv      => topomesh%vert%ntot,  &
            nf      => topomesh%face%ntot,  &
            nc      => topomesh%cell%ntot)

        ! Delete vertices
        !================
        ! Compute vertex ID adjuster
        allocate(diffID(nv))
        diffID = 0
        do i = 1, nv-1 
            if (rmvert(i)) then 
                diffID(i+1) = diffID(i) + 1
            else
                diffID(i+1) = diffID(i)
            end if
        end do

        ! Remove vertices
        topomesh%vert%ntot  = count(keepvert)
        topomesh%vert%ID    = topomesh%vert%ID - diffID 
        topomesh%vert%ID    = pack(topomesh%vert%ID, keepvert)
        topomesh%vert%x     = pack(topomesh%vert%x, keepvert)
        topomesh%vert%y     = pack(topomesh%vert%y, keepvert)
        topomesh%vert%fval  = pack(topomesh%vert%fval, keepvert)
        topomesh%vert%type  = pack(topomesh%vert%type, keepvert) 
        topomesh%vert%fsID  = pack(topomesh%vert%fsID, keepvert)

        ! Adjust faces
        !=============
        ! Loop over all faces
        do i = 1, nf
            ! First vertex
            if (topomesh%face%vert(i, 1) /= 0) then 
                if (rmvert(topomesh%face%vert(i, 1))) then 
                    topomesh%face%vert(i, 1) = 0
                else
                    topomesh%face%vert(i, 1) = topomesh%face%vert(i, 1) - diffID(topomesh%face%vert(i, 1))
                end if
            end if 
            if (topomesh%face%vert(i, 2) /= 0) then 
                if (rmvert(topomesh%face%vert(i, 2))) then 
                    topomesh%face%vert(i, 2) = 0
                else
                    topomesh%face%vert(i, 2) = topomesh%face%vert(i, 2) - diffID(topomesh%face%vert(i, 2))
                end if 
            end if
        end do

        ! Adjust cells
        !=============
        ! Loop over all cells
        do i = 1, nc 
            ! Get vertices
            tv = GetTMCellVert(topomesh%cell, i)
            
            ! Check
            do j = 1, size(tv)
                if (tv(j) /= 0) then
                    if (rmvert(tv(j))) then 
                        tv(j) = 0
                    else
                        tv(j) = tv(j) - diffID(tv(j))
                    end if 
                end if
            end do 
            !tv(delind(tv)) = 0;
            !tv(tv ~= 0) = tv(tv ~= 0) - diffID(tv(tv ~= 0));
            
            ! Reset
            topomesh%cell%vert(topomesh%cell%vertP(i, 1):topomesh%cell%vertP(i, 1)+topomesh%cell%vertP(i, 2)-1) = tv
        end do

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Face removal
    subroutine RemoveTopologicalMeshFaceLogical(topomesh, rmface)

        ! Description
        !============
        ! Remove the faces with IDs specified in rmface from the mesh. 
        ! IMPORTANT: no cell or vertex information is updated here! This should be
        ! recomputed afterwards. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                  :: topomesh 
        logical, intent(in)                 :: rmface(:)

        ! Auxiliary
        integer(I8), allocatable            :: tmp(:, :), diffIDf(:), &
            tcf(:)
        logical, allocatable, dimension(:)  :: keepface

        ! Loop
        integer(I8)                         :: i 

        ! Initialize
        !===========
        ! Unpack
        associate(&
            nf  => topomesh%face%ntot,    &
            nc  => topomesh%cell%ntot)
        keepface = .not. rmface 

        ! Check
        if (size(rmface) /= nf) then 
            call gdErrorHandler('RemoveTopologicalMeshFaceLogical: ' // & 
                'illegal size of rmface')
        end if 

        ! Determine ID shift
        allocate(diffIDf(nf))
        diffIDf = 0
        do i = 1, nf-1 
            if (rmface(i)) then 
                diffIDf(i+1) = diffIDf(i) + 1
            else
                diffIDf(i+1) = diffIDf(i)
            end if
        end do

        ! Delete
        !=======
        ! Update IDs
        topomesh%face%ID = topomesh%face%ID - diffIDf

        ! Remove data
        topomesh%face%ntot = count(keepface)
        tmp = topomesh%face%vert 
        deallocate(topomesh%face%vert)
        allocate(topomesh%face%vert(count(keepface), 2))
        topomesh%face%vert(:, 1) = pack(tmp(:, 1), keepface)
        topomesh%face%vert(:, 2) = pack(tmp(:, 2), keepface)
        topomesh%face%x         = pack(topomesh%face%x, keepface)
        topomesh%face%y         = pack(topomesh%face%y, keepface)
        topomesh%face%fsID      = pack(topomesh%face%fsID, keepface)
        topomesh%face%type      = pack(topomesh%face%type, keepface)
        topomesh%face%ID        = pack(topomesh%face%ID, keepface)
        topomesh%face%pol       = pack(topomesh%face%pol, keepface)  
        
        if (allocated(topomesh%face%cell)) then 
            deallocate(topomesh%face%cell)
        end if 
        if (allocated(topomesh%face%cellP)) then 
            deallocate(topomesh%face%cellP)
        end if 

        ! Update cell face IDs
        do i = 1, nc
            ! Get faces
            tcf = GetTMCellFace(topomesh%cell, i)
            
            ! Reset faces (cells are not removed here!)
            topomesh%cell%face(topomesh%cell%faceP(i, 1):topomesh%cell%faceP(i, 1)+topomesh%cell%faceP(i, 2)-1) = tcf - diffIDf(tcf)
        end do 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Cell removal 
    subroutine RemoveTopologicalMeshCellLogical(topomesh, rmcell)

        ! Description
        !============
        ! Remove the cells with IDs specified in IDs from the mesh. 
        ! IMPORTANT: no face or vertex information is updated here! This should be
        ! recomputed afterwards.

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        logical, intent(in)                     :: rmcell(:)

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: indv, indf, vp, &
            fp
        logical, allocatable, dimension(:)      :: delv, delf 

        ! Loop
        integer(I8)                             :: i, k 

        ! Delete
        !=======
        ! Initialize
        allocate(delv(size(topomesh%cell%vert)), delf(size(topomesh%cell%face)))
        delv = .false. 
        delf = .false. 

        ! Mark vertices and faces for removal
        do i = 1, topomesh%cell%ntot 
            if (rmcell(i)) then 
                indv = [(k, k = topomesh%cell%vertP(i, 1), &
                    topomesh%cell%vertP(i, 1)+topomesh%cell%vertP(i, 2)-1)]
                indf = [(k, k = topomesh%cell%faceP(i, 1), &
                    topomesh%cell%faceP(i, 1)+topomesh%cell%faceP(i, 2)-1)]
                delf(indf) = .true.
                delv(indv) = .true.
            end if 
        end do 

        ! Remove
        topomesh%cell%ntot = count(.not. rmcell)
        topomesh%cell%vert = pack(topomesh%cell%vert, .not. delv)
        topomesh%cell%face = pack(topomesh%cell%face, .not. delf)
        vp = pack(topomesh%cell%vertP(:, 2), .not. rmcell)
        fp = pack(topomesh%cell%faceP(:, 2), .not. rmcell)
        deallocate(topomesh%cell%vertP, topomesh%cell%faceP)
        allocate(topomesh%cell%vertP(topomesh%cell%ntot, 2), &
            topomesh%cell%faceP(topomesh%cell%ntot, 2))
        topomesh%cell%vertP(:, 2) = vp 
        topomesh%cell%faceP(:, 2) = fp 
        topomesh%cell%vertP(1, 1) = 1
        topomesh%cell%faceP(1, 1) = 1
        do i = 2, topomesh%cell%ntot 
            topomesh%cell%vertP(i, 1) = topomesh%cell%vertP(i-1, 1) + &
                topomesh%cell%vertP(i-1, 2)
            topomesh%cell%faceP(i, 1) = topomesh%cell%faceP(i-1, 1) + &
                topomesh%cell%faceP(i-1, 2)
        end do
        topomesh%cell%ntotv = size(topomesh%cell%vert)
        topomesh%cell%ntotf = size(topomesh%cell%face)

    end subroutine

    ! Topological face extraction
    subroutine ExtractTopologicalFacesFromPolygon(pol, eID, sID, &
            xint, yint, v1, v2, xf, yf)

            ! Description
            !============
            ! This function extracts faces from polygons given the segment intersection
            ! and vertex IDs. If the points don't exactly coincide with one of the
            ! segment's vertices, the vertex is added. xint and yint should contain the
            ! coordinates of all intersections, of which the right ones can be queried
            ! with xint(IDs). It is assumed that all intersections are unique. 
            ! Furthermore, it is assumed that the segment intersection 
            ! indices, sID, are already properly sorted (this is not
            ! explicitly checked!)

            ! Notes
            !======
            ! Note 1: it is assumed that the intersections are unique.

            ! Note 2: close polygons are accounted for. It is assumed that the last
            ! point in the polygon.vert structure is the same as the first one. This is
            ! not added twice in facedata.x, y.

            ! Declare variables
            !==================
            ! Arguments
            type(PolygonUDT)                        :: pol 
            integer(I8), allocatable, intent(out)   :: v1(:), v2(:)
            integer(I8), intent(inout)              :: eID(:), sID(:)
            type(RealDynamicArrayUDT), allocatable, intent(out) :: xf(:), yf(:)
            real(R8), intent(in)                    :: xint(:), yint(:)

            ! Auxiliary
            integer(I8)                             :: nf, fc, &
                seID, eeID, ssID, esID 
            integer(I8), allocatable, dimension(:)  :: sortedsID(:), &
                sortedeID(:)
            real(R8)                                :: sx, sy, ex, ey
            real(R8), allocatable, dimension(:)     :: tx, ty

            ! Loop
            integer(I8)                             :: i, si, ei  

            ! Initialize
            !===========
            ! Check maximal number of faces
            nf = size(eID) + 1

            ! Allocate & initialize
            allocate(v1(nf), v2(nf), xf(nf), yf(nf))
            v1 = 0
            v2 = 0
            
            ! Face counter
            fc = 0 

            ! Check
            if (nf == 1) then 
                ! No intersections, just return polygon as is
                xf(1) = ConstructRealDynamicArray(pol%x(pol%vert))
                yf(1) = ConstructRealDynamicArray(pol%y(pol%vert))
                return 
            end if 

            ! Initialize
            do i = 1, nf 
                xf(i) = ConstructRealDynamicArray()
                yf(i) = ConstructRealDynamicArray()
            end do 

            ! Associate
            associate( &
                nsID    =>      size(sID),          & ! size(sID) == size(eID) normally
                neID    =>      size(eID),          &
                px      =>      pol%x(pol%vert),    &
                py      =>      pol%y(pol%vert)     &
                )

            ! Construct segments
            !===================
            sortedsID = sID ! assumed already sorted!
            sortedeID = eID

            ! Treat start and end segments
            if (pol%isclosed) then 
                ! Closed polygon
                if (( (sortedsID(1) == 1 .and. px(1) == xint(sortedeID(1)))) .and. & 
                    ( (sortedsID(nsID) == pol%ne .and. px(pol%ne+1) == xint(sortedeID(nsID))))) then 
                    ! Start and end exactly in start point, do nothing
                elseif (.not. ( (sortedsID(1) == 1 .and. px(1) == xint(sortedeID(1)))) .and. & 
                    .not. ( (sortedsID(nsID) == pol%ne .and. px(pol%ne+1) == xint(sortedeID(nsID))))) then 
                    ! Start and end segment without start and end vertex,
                    ! so merge
                    fc = fc + 1

                    ! Add vertices
                    v1(fc) = sortedeID(neID)
                    v2(fc) = sortedeID(1)

                    ! Add face coordinates (hedge for exact duplicates)
                    tx = [xint(sortedeID(neID)), px(sortedsID(nsID)+1:pol%ne)]
                    ty = [yint(sortedeID(neID)), py(sortedsID(nsID)+1:pol%ne)]
                    tx = [tx, px(1:sortedsID(1)), xint(sortedeID(1))]
                    ty = [ty, py(1:sortedsID(1)), yint(sortedeID(1))]
                    if ((tx(1) == tx(2)) .and. (ty(1) == ty(2))) then 
                        tx = tx(2:size(tx))
                        ty = ty(2:size(ty))
                    end if 
                    if ((tx(size(tx)-1) == tx(size(tx))) .and. &
                        (ty(size(ty)-1) == ty(size(ty)))) then 
                        tx = tx(1:size(tx)-1)
                        ty = ty(1:size(ty)-1)
                    end if 

                    call xf(fc)%Append(tx)
                    call yf(fc)%Append(ty)
                else
                    ! Something weird here: intersection in exactly first node, but
                    ! only at first and not at last segment. Likely because
                    ! intersection with last segment was not added. Throw error. 
                    call gdErrorHandler(&
                        'ExtractTopologicalFacesFromPolygon: ' // & 
                        'closed polygon intersects in start point ' // & 
                        'but not in end point - check input '// & 
                        '(end point intersection should be added separately)')
                end if 
            else 
                ! Open polygon
                if ((sortedsID(1) == 1) .and. (px(1) == xint(sortedeID(1))) & 
                    .and. (py(1) == yint(sortedeID(1)))) then 
                    ! Start veretx at first intersection, no starting segment
                    ! without starting vertex 
                else
                    ! Start at second intersection as indicated, add starting
                    ! segment 
                    fc = fc + 1

                    ! Add vertices
                    v1(fc) = 0 ! no start vertex
                    v2(fc) = sortedeID(1) ! end vertex 

                    ! Add face coordinates
                    tx = [px(1:sortedsID(1)), xint(sortedeID(1))]
                    ty = [py(1:sortedsID(1)), yint(sortedeID(1))]
                    if ((tx(size(tx)-1) == tx(size(tx))) .and. &
                        (ty(size(ty)-1) == ty(size(ty)))) then 
                        tx = tx(1:size(tx)-1)
                        ty = ty(1:size(ty)-1)
                    end if 
                    call xf(fc)%Append(tx) 
                    call yf(fc)%Append(ty)
                
                end if 
                if ((sortedsID(nsID) == pol%ne) .and. (px(pol%ne+1) == xint(sortedeID(nsID))) & 
                    .and. (py(pol%ne+1) == yint(sortedeID(nsID)))) then 
                    ! Start veretx at first intersection, no starting segment
                    ! without starting vertex 
                else
                    ! Start at second intersection as indicated, add starting
                    ! segment 
                    fc = fc + 1

                    ! Add vertices
                    v1(fc) = sortedeID(neID) ! start vertex
                    v2(fc) = 0 ! no end vertex 

                    ! Add face coordinates
                    tx = [xint(sortedeID(neID)), px(sortedsID(nsID)+1:pol%ne+1)]
                    ty = [yint(sortedeID(neID)), py(sortedsID(nsID)+1:pol%ne+1)]
                    if ((tx(1) == tx(2)) .and. (ty(1) == ty(2))) then 
                        tx = tx(2:size(tx))
                        ty = ty(2:size(ty))
                    end if 
                    call xf(fc)%Append(tx) 
                    call yf(fc)%Append(ty)
                
                end if 
            end if 

            ! Set start and end indices for other segments
            si = 1
            ei = neID-1

            ! Loop over remaining segments
            do i = si, ei 
                ! All segments should be nicely closed so no special 
                ! checks needed

                ! Get start and end points & segments
                seID = sortedeID(i)
                eeID = sortedeID(i+1)
                ssID = sortedsID(i)
                esID = sortedsID(i+1)
                sx = xint(seID)
                sy = yint(seID)
                ex = xint(eeID)
                ey = yint(eeID)

                ! Update counter
                fc = fc + 1

                ! Add vertices
                v1(fc) = seID 
                v2(fc) = eeID 

                ! Add coordinates
                tx = px(ssID+1:esID)
                ty = py(ssID+1:esID)

                if (.not. allocated(tx)) then 
                    allocate(tx(0), ty(0))
                end if 

                ! Check
                if (size(tx) > 0) then 
                    ! Check if we need to add start and end points
                    if ((tx(1) /= sx) .or. (ty(1) /= sy)) then 
                        tx = [sx, tx]
                        ty = [sy, ty]
                    end if 
                    if ((tx(size(tx)) /= ex) .or. (ty(size(ty)) /= ey)) then 
                        tx = [tx, ex]
                        ty = [ty, ey]
                    end if 
                else 
                    ! Only add start and end points
                    tx = [sx, ex]
                    ty = [ex, ey]
                end if 

                ! Add
                call xf(fc)%Append(tx)
                call yf(fc)%Append(ty)

            end do

            ! Check if we need to shrink (may happen for closed polygons)
            if (fc < nf) then 
                xf = xf(1:fc)
                yf = yf(1:fc)
                V1 = V1(1:fc)
                V2 = V2(1:fc)
            end if 

            ! Housekeeping
            !=============
            end associate





    end subroutine 

    ! Getters
    function GetTMVertFace(vert, i) result(res)
        integer(I8)                 :: i 
        class(TopomeshVertUDT)      :: vert 
        integer(I8), allocatable    :: res(:)
        res = vert%face(vert%faceP(i, 1):(vert%faceP(i, 1) +  vert%faceP(i, 2) - 1))
    end function

    function GetTMCellVert(cell, i) result(res)
        integer(I8)                 :: i 
        class(TopomeshCellUDT)      :: cell 
        integer(I8), allocatable    :: res(:)
        res = cell%vert(cell%vertP(i, 1):(cell%vertP(i, 1) + cell%vertP(i, 2) - 1))
    end function

    function GetTMCellFace(cell, i) result(res)
        integer(I8)                 :: i 
        class(TopomeshCellUDT)      :: cell 
        integer(I8), allocatable    :: res(:)
        res = cell%face(cell%faceP(i, 1):(cell%faceP(i, 1) + cell%faceP(i, 2) - 1))
    end function

    function GetTMFaceCell(face, i) result(res)
        integer(I8)                 :: i 
        class(TopomeshFaceUDT)      :: face 
        integer(I8), allocatable    :: res(:)
        res = face%cell(face%cellP(i, 1):(face%cellP(i, 1) + face%cellP(i, 2) - 1))
    end function

    function GetTMTubeCell(tube, i) result(res)
        integer(I8)                 :: i 
        class(TopomeshTubeUDT)      :: tube 
        integer(I8), allocatable    :: res(:)
        res = tube%cell(tube%cellP(i, 1):(tube%cellP(i, 1) + tube%cellP(i, 2) - 1))
    end function

    function GetTMTubeFace(tube, i) result(res)
        integer(I8)                 :: i 
        class(TopomeshTubeUDT)      :: tube 
        integer(I8), allocatable    :: res(:)
        res = tube%face(tube%faceP(i, 1):(tube%faceP(i, 1) + tube%faceP(i, 2) - 1))
    end function

    function GetTMTubeBndFace(tube, i, j) result(res)
        class(TopomeshTubeUDT)      :: tube
        integer(I8), intent(in)     :: i, j
        integer(I8), allocatable    :: res(:)
        
        if (j == 1) then 
            res = tube%bndf1(tube%bndf1P(i, 1):(tube%bndf1P(i, 1) + tube%bndf1P(i, 2) - 1))
        else 
            res = tube%bndf2(tube%bndf2P(i, 1):(tube%bndf2P(i, 1) + tube%bndf2P(i, 2) - 1))
        end if 
    end function

    function GetTMTubeBndVert(tube, i, j) result(res)
        class(TopomeshTubeUDT)      :: tube
        integer(I8), intent(in)     :: i, j
        integer(I8), allocatable    :: res(:)

        if (j == 1) then 
            res = tube%bndv1(tube%bndv1P(i, 1):(tube%bndv1P(i, 1) + tube%bndv1P(i, 2) - 1))
        else 
            res = tube%bndv2(tube%bndv2P(i, 1):(tube%bndv2P(i, 1) + tube%bndv2P(i, 2) - 1))
        end if 
    end function

    function GetTMTubeNeig(tube, i, j) result(res)
        class(TopomeshTubeUDT)      :: tube
        integer(I8), intent(in)     :: i, j
        integer(I8), allocatable    :: res(:)

        if (j == 1) then 
            res = tube%ftneig1(tube%ftneig1P(i, 1):(tube%ftneig1P(i, 1) + tube%ftneig1P(i, 2) - 1))
        else 
            res = tube%ftneig2(tube%ftneig2P(i, 1):(tube%ftneig2P(i, 1) + tube%ftneig2P(i, 2) - 1))
        end if 
    end function

    function GetTMVertFaceNeig(vert, i, tf) result(res)

        ! Description
        !============
        ! Return the 'left' and 'right' neighbouring faces for vertex tv and
        ! current face tf. If the face is not found, fn is empty. Otherwise, fn is
        ! a 1-by-2 integer array containing the face indices.

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)     :: i, tf 
        class(TopomeshVertUDT)      :: vert 
        integer(I8), allocatable    :: res(:)

        ! Auxiliary
        integer(I8)                     :: tfind
        integer(I8), allocatable        :: tvf(:), tvfe(:)

        ! Get vertex faces
        tvf = GetTMVertFace(vert, i)

        ! Find where this face is located
        tfind = findloc(tvf, tf, 1)

        ! Find the neighbours 
        if (tfind == 0) then 
            allocate(res(0))
        else
            tvfe = [tvf(size(tvf)), tvf, tvf(1)]
            res = [tvfe(tfind), tvfe(tfind+2)]
        end if 
    end function

    ! ID getters
    function GetInternalFaceIDs(topomesh) result(ID)

        ! Description
        !============
        ! Get all internal face IDs (simply faces with two cell 
        ! neighbours)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)          :: topomesh 
        integer(I8), allocatable    :: ID(:)

        ! Loop
        integer(I8)                 :: k 

        ! Get
        !====
        allocate(ID(count(.not. topomesh%face%BF)))
        ID = pack([(k, k = 1, topomesh%face%ntot)], .not. topomesh%face%BF)

    end function

    function GetBoundaryFaceIDs(topomesh) result(ID)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)          :: topomesh 
        integer(I8), allocatable    :: ID(:)

        ! Loop
        integer(I8)                 :: k 

        ! Get
        !====
        allocate(ID(count(topomesh%face%BF)))
        ID = pack([(k, k = 1, topomesh%face%ntot)], topomesh%face%BF)

    end function

    function GetSeparatrixFaceIDs(topomesh) result(ID)

        ! Description
        !=============
        ! Get all face IDs that are separatrix parts

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)          :: topomesh 
        integer(I8), allocatable    :: ID(:)

        ! Auxiliary
        logical, allocatable        :: temp(:)

        ! Loop
        integer(I8)                 :: k 

        ! Get
        !====
        temp = topomesh%face%type == TMfacesepID
        allocate(ID(count(temp)))
        ID = pack([(k, k = 1, topomesh%face%ntot)], temp)

    end function

    function GetCoreFaceIDs(topomesh) result(ID)

        ! Description
        !=============
        ! Get all face IDs that are core parts

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)          :: topomesh 
        integer(I8), allocatable    :: ID(:)

        ! Auxiliary
        logical, allocatable        :: temp(:)

        ! Loop
        integer(I8)                 :: k 

        ! Get
        !====
        temp = topomesh%face%type == TMfacecoreID
        allocate(ID(count(temp)))
        ID = pack([(k, k = 1, topomesh%face%ntot)], temp)

    end function

    function GetVesselFaceIDs(topomesh) result(ID)

        ! Description
        !=============
        ! Get all face IDs that are vessel parts (i.e. type TMfacebndID)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)          :: topomesh 
        integer(I8), allocatable    :: ID(:)

        ! Auxiliary
        logical, allocatable        :: temp(:)

        ! Loop
        integer(I8)                 :: k 

        ! Get
        !====
        temp = topomesh%face%type == TMfacebndID
        allocate(ID(count(temp)))
        ID = pack([(k, k = 1, topomesh%face%ntot)], temp)

    end function

    function GetTargetFaceIDs(topomesh) result(ID)

        ! Description
        !=============
        ! Get all face IDs that are vessel parts (i.e. type TMfacebndID)
        ! and of which at least one vertex lies on a separatrix face.

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)          :: topomesh 
        integer(I8), allocatable    :: ID(:)

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: tempID, tv, tvf
        logical, allocatable, dimension(:)      :: temp

        ! Loop
        integer(I8)                 :: i, j

        ! Get
        !====
        ! Find vessel faces
        tempID = topomesh%GetVesselFaceIDs()

        ! Check which vessel faces have a separatrix vertex
        allocate(temp(size(tempID)))
        temp = .false. 
        do i = 1, size(tempID)
            ! Get face vertices
            tv = topomesh%face%vert(tempID(i), :)

            ! Check
            do j = 1, size(tv)
                tvf = topomesh%vert%GetFace(tv(j))
                if (any(topomesh%face%type(tvf) == TMfacesepID)) then 
                    temp(i) = .true.
                    exit 
                end if 
            end do 
        end do 

        ! Get indices
        allocate(ID(count(temp)))
        ID = pack(tempID, temp)

    end function

    function GetLastFluxSurfaceFaceIDs(topomesh) result(ID)

        ! Description
        !============
        ! Get faces that form last flux surfaces (e.g. PF or SOL
        ! boundary) - these are boundary faces that are not 
        ! vessel faces

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)          :: topomesh 
        integer(I8), allocatable    :: ID(:)

        ! Auxiliary
        logical, allocatable        :: temp(:)

        ! Loop
        integer(I8)                 :: k 

        ! Get
        !====
        temp = (topomesh%face%type /= TMfacebndID) .and. &
            topomesh%face%BF 
        allocate(ID(count(temp)))
        ID = pack([(k, k = 1, topomesh%face%ntot)], temp)
        
    end function

    function GetCoreCellIDs(topomesh) result(ID)

        ! Description
        !============
        ! Get the IDs of cells that are core regions. These either have
        ! a core boundary as a face, or they have a maximum/minimum as
        ! a vertex.

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                  :: topomesh 
        integer(I8), allocatable            :: ID(:)

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: tcf, tcv 
        logical, allocatable, dimension(:)      :: keepind

        ! Loop
        integer(I8)                         :: i, k

        ! Loop
        !=====
        ! Initialize
        allocate(keepind(topomesh%cell%ntot))
        keepind = .false. 

        ! Loop
        do i = 1, topomesh%cell%ntot 
            ! Check faces & vert
            tcf = topomesh%cell%GetFace(i)
            tcv = topomesh%cell%GetVert(i)

            ! Check
            if (any(topomesh%face%type(tcf) == TMfacecoreID)) then 
                keepind(i) = .true. 
            end if 
            if ((any(topomesh%vert%type(tcv) == TMvertexmaxID)) .or. &
                (any(topomesh%vert%type(tcv) == TMvertexminID))) then 
                    keepind(i) = .true.
            end if 
        end do 

        ! Extract IDs
        allocate(ID(count(keepind)))
        ID = pack([(k, k = 1, topomesh%cell%ntot)], keepind)

    end function

    function GetWideGridCellIDs(topomesh) result(ID)

        ! Description
        !============
        ! This function returns all the cells that are consider to be a 
        ! 'wide grid' cell. This includes basically all cells that 
        ! are not next to a separatrix, or cells that have an extremum.
        ! Note that if no core boundaries are added and the grid extends
        ! to an extremum, even the typical 'core' cell is a wide grid
        ! cell!

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        integer(I8), allocatable, dimension(:)  :: ID 

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: tcf, tcv
        logical, allocatable, dimension(:)      :: keepind

        ! Loop
        integer(I8)                             :: i 

        ! Initialize
        !===========
        ! Keeper
        allocate(keepind(topomesh%cell%ntot))
        keepind = .true. 

        ! Determine
        !==========
        do i = 1, topomesh%cell%ntot 
            ! Get faces and vertices of this cell
            tcf = topomesh%cell%GetFace(i)
            tcv = topomesh%cell%GetVert(i)

            ! Any extrema?
            if ((any(topomesh%vert%type(tcv) == TMvertexmaxID)) .or. & 
                (any(topomesh%vert%type(tcv) == TMvertexminID))) then 
                keepind(i) = .false. 
                cycle
            end if 

            ! Any separatrix parts? 
            if (.not. (any(topomesh%face%type(tcf) == TMfacesepID))) then 
                keepind(i) = .false.
            end if 
        end do 

        ! Get IDs
        allocate(ID(count(keepind)))
        ID = pack([(i, i = 1, topomesh%cell%ntot)], keepind)

    end function

    function GetSeparatrixFluxSurfaceIDs(topomesh) result(fsID)

        ! Description
        !============
        ! Get separatrix flux surface IDs by checking fsIDs of all 
        ! separatrix faces

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        integer(I8), allocatable, dimension(:)  :: fsID
        
        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: fsIDs
        logical, allocatable, dimension(:)      :: issepface

        ! Loop

        ! Loop over faces
        !================
        issepface = topomesh%face%type == TMfacesepID
        allocate(fsIDs(count(issepface)))
        fsIDs = pack(topomesh%face%fsID, issepface)
        call Unique(fsIDs, fsID)

    end function

    function GetStrikePointIDs(topomesh) result(ID)

        ! Description
        !============
        ! This function returns the vertex indices that are strike 
        ! points. Strike points are defined as boundary points of which
        ! at least one face is a separatrix segment.
        
        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)          :: topomesh 
        integer(I8), allocatable    :: ID(:)

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: tempID, tv
        logical, allocatable, dimension(:)      :: temp

        ! Loop
        integer(I8)                 :: i, j

        ! Get
        !====
        ! Find separatrix faces
        tempID = topomesh%GetSeparatrixFaceIDs()

        ! Check which separatrix faces have a boundary vertex
        allocate(temp(topomesh%vert%ntot))
        temp = .false. 
        do i = 1, size(tempID)
            ! Get face vertices
            tv = topomesh%face%vert(tempID(i), :)

            ! Check
            do j = 1, size(tv)
                if (topomesh%vert%type(tv(j)) == TMvertexbndID) then 
                    temp(tv(j)) = .true.
                end if 
            end do 
        end do 

        ! Get indices
        allocate(ID(count(temp)))
        ID = pack([(j, j = 1, topomesh%vert%ntot)], temp)

    end function

    function GetXPointIDs(topomesh) result(ID)

        ! Description
        !============
        ! This function returns all the X-point (saddle point) IDs

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh
        integer(I8), allocatable, dimension(:)  :: ID

        ! Auxiliary
        logical, allocatable, dimension(:)      :: isXP

        ! Loop
        integer(I8)                             :: k

        ! Find X-points
        !==============
        isXP = topomesh%vert%type == TMvertexsaddleID
        allocate(ID(count(isXP)))
        ID = pack([(k, k = 1, size(isXP))], isXP)

    end function

    function GetOPointIDs(topomesh) result(ID)

        ! Description
        !============
        ! This function returns all the O-points (extrema) point) IDs

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh
        integer(I8), allocatable, dimension(:)  :: ID

        ! Auxiliary
        logical, allocatable, dimension(:)      :: isOP

        ! Loop
        integer(I8)                             :: k

        ! Find X-points
        !==============
        isOP = (topomesh%vert%type == TMvertexminID) .or. &
            (topomesh%vert%type == TMvertexmaxID)
        allocate(ID(count(isOP)))
        ID = pack([(k, k = 1, size(isOP))], isOP)

    end function

    function GetPrimaryXPointIDs(topomesh) result(ID)

        ! Description
        !============
        ! This function returns the vertex IDs of all x-points that 
        ! lie close to a core (close in the sense that there is no other
        ! x-point between this x-point and the next o-point). This can 
        ! only occur if the x-point is present in a core region's 
        ! vertices. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh
        integer(I8), allocatable, dimension(:)  :: ID

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: corecells, tv

        ! Loop
        integer(I8)                             :: i

        ! Initialize
        !===========
        ! Allocate
        allocate(ID(0))

        ! Unpack
        associate(&
            cell        => topomesh%cell,   &
            vert        => topomesh%vert    &
            )

        ! Determine primary x-points
        !===========================
        ! Determine core cells
        corecells = topomesh%GetCoreCellIDs()
        do i = 1, size(corecells)
            ! Get vertices
            tv = cell%GetVert(corecells(i))

            ! Check which ones are x-points (if any) and add
            ID = [ID, pack(tv, vert%type(tv) == TMvertexsaddleID)]
        end do 

        ! Housekeeping
        !=============
        end associate

    end function

    function GetStrikePointXPointIDs(topomesh) result(ID)

        ! Description
        !============
        ! This function returns the x-point ID for each strike point 
        ! as returned by GetStrikePointIDs (a strike point can only 
        ! have one X-point normally speaking). To determine this, we 
        ! start from each strike point and keep tracing the aligned 
        ! field line it belongs to until an x-point is reached (normally, 
        ! this is already reached in the first face but ok)
        
        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)          :: topomesh 
        integer(I8), allocatable    :: ID(:)

        ! Auxiliary
        integer(I8)                             :: nsp, tspfloc, &
            nextv, tsp
        integer(I8), allocatable, dimension(:)  ::sp, tspf, tfvert

        ! Loop
        integer(I8)                 :: i

        ! Initialize
        !===========
        ! Get strike points
        sp = topomesh%GetStrikePointIDs()
        nsp = size(sp)
        
        ! Initialize array
        allocate(ID(nsp))
        ID = 0

        ! Associate
        associate(&
            vert    => topomesh%vert,   &
            face    => topomesh%face    &
            )

        ! Loop
        !=====
        do i = 1, nsp
            ! Unpack 
            tsp = sp(i)

            ! Get faces of this strike point
            tspf = vert%GetFace(tsp)

            ! Get separatrix face
            tspfloc = findloc(face%type(tspf), TMfacesepID, 1)

            ! Sanity check
            if (tspfloc == 0) then 
                call gdErrorHandler('GetStrikePointXPointIDs: could not ' // & 
                    'find separatrix face, unexpected')
            end if 

            ! Get the next vertex
            tfvert = face%vert(tspf(tspfloc), :)
            if (tfvert(1) == tsp) then 
                nextv = tfvert(2)
            elseif (tfvert(2) == tsp) then 
                nextv = tfvert(1)
            else
                call gdErrorHandler('GetStrikePointXPointIDs: something ' // & 
                    'wrong with topomesh interconnections, could not ' // & 
                    'find current vertex in face vertices')
            end if 

            ! Check the point type
            if (vert%type(nextv) == TMvertexsaddleID) then 
                ! Found, add and cycle
                ID(i) = nextv 
                cycle
            elseif (vert%type(nextv) == TMvertexbndID) then 
                ! Shouldn't happen, error
                call gdErrorHandler('GetStrikePointXPointIDs: could not ' // & 
                    'find X-point of strike point')
            end if 

            ! If we got here, we need to keep on walking along the 
            ! separatrix - this is not yet implemented
            call gdErrorHandler('GetStrikePointXPointIDs: strike points ' // & 
                'that do not directly connect to an x-point are not yet ' // & 
                'supported')

        end do

        ! Housekeeping
        !=============
        end associate

    end function

    function GetClosedContourTangencyPointIDs(topomesh) result(ID)

        ! Description
        !============
        ! This function returns the tangency point IDs (type 2) that 
        ! have a contour which closes upon itself. 
        
        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)          :: topomesh 
        integer(I8), allocatable    :: ID(:)

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: tv, &
            sortind, tf
        integer(I8), allocatable, dimension(:, :)   :: tfvert
        logical, allocatable, dimension(:)      :: isbranchingpolygon, &
            ispolygonstart

        ! Loop
        integer(I8)                 :: i, k

        ! Initialize
        !===========
        ! Allocate
        allocate(ID(0))

        ! Associate
        associate(&
            face        => topomesh%face,   &
            vert        => topomesh%vert    &
            )

        ! Loop
        !=====
        do i = 1, vert%ntot
            if (vert%type(i) == TMvertextp2ID) then 
                ! Get all topomesh faces with this flux surface ID
                allocate(tf(count(face%fsID == vert%fsID(i))))
                tf = pack([(k, k = 1, face%ntot)], face%fsID == vert%fsID(i))

                ! Check by sorting vertices 
                tfvert = face%vert(tf, :)
                allocate(sortind(size(tf)), ispolygonstart(size(tf)), &
                    isbranchingpolygon(size(tf)))
                call SortPolygonEdges(tfvert, size(tf), sortind, ispolygonstart, &
                    isbranchingpolygon)
                tfvert = tfvert(sortind, :)

                ! Sanity checks
                if (count(ispolygonstart) /= 1) then 
                    call gdErrorHandler('GetClosedContourTangencyPointIDs: ' // & 
                        'found no or multiple polygons, unexpected')
                end if 
                if (count(isbranchingpolygon) /= 0) then 
                    call gdErrorHandler('GetClosedContourTangencyPointIDs: ' // &
                        'found branching polygons, unexpected')
                end if 

                ! Check first and last edge
                call SetDiff(tfvert(1, :), tfvert(size(tf), :), tv)
                if (size(tv) == 0) then 
                    ! Found, add
                    ID = [ID, i]
                end if 

                ! Housekeeping
                deallocate(tf, sortind, ispolygonstart, isbranchingpolygon)

            end if 
        end do

        ! Housekeeping
        !=============
        end associate

    end function
 
    ! Topological mesh identification 
    function IdentifyTopologicalMeshType(topomesh) result(TMlabel)

        ! Description
        !============
        ! This function attempts to identify the topological mesh type
        ! and returns the topological mesh identification number 
        ! (see mod_definitions of the definition thereof). This is not 
        ! really required by goat itself, but is rather a convenience
        ! tool for later simulation etc and face label mapping. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)              :: topomesh
        integer(I8)                     :: TMlabel

        ! Auxiliary
        logical                                 :: issinglenull, &
            isdoublenull
        integer(I8)                             :: nxp, nop, nwgc, nsp, &
            ncc
        integer(I8), allocatable, dimension(:)  :: xp, op, wgc, sp, cc

        ! Initialize
        !===========
        ! Set initial label to default label
        TMlabel = TMTopGeneral

        ! Get all X-points
        xp = topomesh%GetXPointIDs()
        nxp = size(xp)

        ! Get all O-points
        op = topomesh%GetOPointIDs()
        nop = size(op)

        ! Strike points
        sp = topomesh%GetStrikePointIDs()
        nsp = size(sp)

        ! Get all wide grid cells
        wgc = topomesh%GetWideGridCellIDs()        
        nwgc = size(wgc)

        ! Core cells
        cc = topomesh%GetCoreCellIDs()
        ncc = size(cc)

        ! Initialize
        issinglenull = .true.
        isdoublenull = .true.

        ! Single null
        !============        
        ! X-point and O-point checks
        if (nxp /= 1 .or. nop /= 0) then 
            issinglenull = .false.
        end if 

        ! Strike point checks
        if (nsp /= 2) then 
            issinglenull = .false. 
        end if 

        ! Number of wide and narrow grid cells
        if (nwgc > 0) then 
            issinglenull = .false.
        elseif (topomesh%cell%ntot /= 3) then
            issinglenull = .false.
        end if 

        ! Core cells
        if (ncc /= 1) then 
            issinglenull = .false.
        end if 

        ! Double null
        !============
        ! X-point and O-point checks
        if (nxp /= 2 .or. nop /= 0) then 
            isdoublenull = .false.
        end if 

        ! Strike point checks
        if (nsp /= 4) then 
            isdoublenull = .false. 
        end if 

        ! Number of wide and narrow grid cells
        if (nwgc > 0) then 
            isdoublenull = .false.
        elseif (topomesh%cell%ntot /= 3) then
            isdoublenull = .false.
        end if 

        ! Core cells
        if (ncc /= 2) then 
            isdoublenull = .false.
        end if 

        ! Determine flag
        !===============
        ! Sanity checks
        if (count([issinglenull, isdoublenull]) > 1) then 
            ! Probably we missed something in the definition then
            print *, 'IdentifyTopologicalMeshType: multiple topologies ' // & 
                'appear valid, this is likely a bug. Setting flag to ' // & 
                'general flag...'

        elseif (issinglenull) then 
            TMlabel = TMTopSN
        elseif (isdoublenull) then 
            TMlabel = TMTopDN
        end if 

    end function

    ! Topological metric computations
    function GetTMFaceDeltaPsi(topomesh, faceID) result(dpsi)

        ! Description
        !============
        ! This function computes the (absolute) difference in psi value
        ! for the faces given in the array faceID. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT), intent(in)      :: topomesh 
        integer(I8), intent(in)             :: faceID(:) 
        real(R8), allocatable, dimension(:) :: dpsi 

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: tfv1, tfv2 
        real(R8), allocatable, dimension(:)     :: fval1, fval2

        ! Compute
        !========
        ! Get vertices of faces
        tfv1 = topomesh%face%vert(faceID, 1)
        tfv2 = topomesh%face%vert(faceID, 2)
        fval1 = topomesh%fsfval%Get(topomesh%vert%fsID(tfv1))
        fval2 = topomesh%fsfval%Get(topomesh%vert%fsID(tfv2))
        dpsi = abs(fval1 - fval2)

    end function

    function GetTMFaceRadialLength(topomesh, magneticField,faceID) &
        result(lrad)

        ! Description
        !============
        ! This function computes the (absolute) length of a face along
        ! the radial direction. This is determined as the sum of the 
        ! radial lengths of the face's edges. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT), intent(in)      :: topomesh 
        integer(I8), intent(in)             :: faceID(:) 
        type(MagneticFieldUDT), intent(in)  :: magneticField 
        real(R8), allocatable, dimension(:) :: lrad 

        ! Auxiliary
        integer(I8)                         :: nf    
        real(R8), allocatable, dimension(:) :: x, y, dx, dy, bx, by, &
            xf, yf, bn

        ! Loop
        integer(I8)                             :: i 

        ! Initialize
        !===========
        nf = size(faceID)
        allocate(lrad(nf))
        lrad = 0.0_R8

        ! Compute
        !========
        do i = 1, nf
            ! Get face coordinates
            x = topomesh%face%x(faceID(i))%Get()
            y = topomesh%face%y(faceID(i))%Get()

            ! Compute edge center coordinates and lengths, and psi gradient
            dx = x(2:) - x(1:size(x)-1)
            dy = y(2:) - y(1:size(y)-1)
            xf = 0.5*(x(2:) + x(1:size(x)-1))
            yf = 0.5*(y(2:) + y(1:size(y)-1))
            allocate(bx(size(xf)), by(size(xf)))
            call magneticField%interp%Evaluate(xf, yf, 1, 0, bx)
            call magneticField%interp%Evaluate(xf, yf, 0, 1, by)
            bn = sqrt(bx**2 + by**2)
            bx = bx/bn
            by = by/bn

            ! Compute length
            lrad(i) = sum(abs(bx*dx + by*dy))

            ! Housekeeping
            deallocate(bx, by)
        end do 

    end function

    !------------------------------------------------------------------!
    !                 TOPOLOGICAL MESH CELL OPERATORS                  !
    !------------------------------------------------------------------!
    ! Cell face sorter (CW or CCW)
    subroutine SortTopologicalMeshCellFaces(topomesh, cellind, &
        sortedfaces, sortedvertices, doflip)

        ! Description
        !============
        ! This routine sorts the faces of a cell in arbitrary direction
        ! (clockwise or counter-clockwise) and returns the sorted face
        ! IDs (sortedfaces) and an indication whether the face 
        ! coordinates should be flipped when querying the face from 
        ! the topological mesh. This may be useful to post-process 
        ! cells etc. Additionally, the vertices of the cell are also
        ! returned (closed form, so first and last are the same).

        ! Note: by construction, the cell faces should already be 
        ! sorted correctly! 

        ! Note: we should hedge for possible duplicate faces since 
        ! disc regions are allowed

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        integer(I8), intent(in)                 :: cellind 
        integer(I8), allocatable, intent(out)   :: sortedfaces(:), &
            sortedvertices(:)
        logical(I8), allocatable, intent(out)   :: doflip(:)

        ! Auxiliary
        integer(I8)                             :: startvert, sf
        integer(I8), allocatable, dimension(:)  :: tf, tfv, tdiff

        ! Loop
        integer(I8)                             :: i 

        ! Initialize
        !===========
        ! Associate
        associate(&
            cell        => topomesh%cell,   &
            face        => topomesh%face)

        ! Unpack 
        sortedfaces = GetTMCellFace(cell, cellind)
        sortedvertices = GetTMCellVert(cell, cellind)
        tf = sortedfaces

        ! Initialize
        if (allocated(doflip)) then 
            deallocate(doflip)
        end if 
        allocate(doflip(size(sortedfaces)))
        doflip = .false.

        ! Determine starting face as non-duplicate face (then all the rest
        ! should go fine in principle)
        tdiff = tf(2:size(tf)) - tf(1:size(tf)-1)
        sf = findloc(tdiff /= 0, .true., 1)
        if (sf == 1) then 
            ! Nothing to see here, move along
        elseif (sf == 2) then 
            ! Shift the first face backwards
            tf = [tf(2:size(tf)), tf(1)]
        else
            ! Unexpected, throw error
            call gdErrorHandler('SortTopologicalMeshCellFaces: more than two faces ' // & 
                'appear to be the same for this cell, unexpected')
        end if 

        ! Loop
        !=====
        if (any(face%vert(tf(1), 1) == face%vert(tf(2), :))) then 
            startvert = face%vert(tf(1), 1)
            doflip(1) = .true. 
        elseif (any(face%vert(tf(1), 2) == face%vert(tf(2), :))) then 
            startvert = face%vert(tf(1), 2)
            doflip(1) = .false. 
        else 
            call gdErrorHandler('SortTopologicalMeshCellFaces: unexpected error')
        end if

        do i = 2, size(tf)
            tfv = face%vert(tf(i), :)
            if (tfv(2) == startvert) then 
                doflip(i) = .true. 
                startvert = tfv(1)
            elseif (tfv(1) == startvert) then 
                doflip(i) = .false.
                startvert = tfv(2)
            else 
                call gdErrorHandler('SortTopologicalMeshCellFaces: unexpected error')
            end if 
        end do 
        


        ! Housekeeping
        !=============
        end associate 

    end subroutine

    ! Full cell polygon coordinate constructor
    subroutine ConstructTopologicalMeshCellPolygon(topomesh, cellind, &
        xcda, ycda)

        ! Description
        !============
        ! Construct the cell polygon by appending the face polygon
        ! coordinates. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        integer(I8), intent(in)                 :: cellind
        type(RealDynamicArrayUDT), intent(out)  :: xcda, ycda

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: sortedfaces, &
            sortedvertices
        real(R8), allocatable, dimension(:)     :: tempx, tempy 
        logical, allocatable, dimension(:)      :: doflip

        ! Loop
        integer(I8)                             :: i 

        ! Initialize
        !===========
        ! Initialize dynamic arrays
        xcda = ConstructRealDynamicArray()
        ycda = ConstructRealDynamicArray()

        ! Construct polygon
        !==================
        ! Sort the cell faces
        call SortTopologicalMeshCellFaces(topomesh, cellind, sortedfaces, &
            sortedvertices, doflip)

        ! Loop over the faces 
        do i = 1, size(sortedfaces)
            ! Get face coordinates
            tempx = topomesh%face%x(sortedfaces(i))%Get()
            tempy = topomesh%face%y(sortedfaces(i))%Get()

            ! Check
            if (doflip(i)) then 
                call xcda%Append(tempx(size(tempx):2:-1)) ! avoid double coordinates
                call ycda%Append(tempy(size(tempy):2:-1))
            else
                call xcda%Append(tempx(2:size(tempx))) ! avoid double coordinates
                call ycda%Append(tempy(2:size(tempy)))
            end if 
        end do 

    end subroutine 

    !------------------------------------------------------------------!
    !                            AUXILIARY                             !
    !------------------------------------------------------------------!
    ! Tiny newton solver
    subroutine TinyNewtonSolver(x, y, converged, x0, y0, magneticField)

        ! Description
        !============
        ! Tiny newton solver implementation to refine extrema. 
        ! No guarantee on convergence (also not strictly needed)

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(out)               :: x, y
        real(R8), intent(in)                :: x0, y0 
        logical, intent(out)                :: converged 
        type(magneticFieldUDT)              :: magneticField    
        
        ! Auxiliary
        integer(I8), parameter              :: maxit = 10
        real(R8), parameter                 :: tol = 1e-10 
        integer(I8)                         :: flag
        real(R8)                            :: xi(1), yi(1), dFdxi(1), &
            dFdyi(1), d2Fdx2(1), d2Fdy2(1), d2Fdxdy(1), A(1:2, 1:2), &
            res(1:2), d(1:2)

        ! Loop 
        integer(I8)                         :: it 

        ! Initialize
        !===========
        converged = .false.
        it = 1
        xi = x0
        yi = y0
        do while ((it <= maxit) .and. (.not. converged) )
            ! Compute residuals
            call magneticField%interp%Evaluate(xi, yi, 1, 0, dFdxi)
            call magneticField%interp%Evaluate(xi, yi, 1, 0, dFdyi)

            ! Check convergence 
            res = -[dFdxi, dFdyi]
            if (all(abs(res) <= tol)) then 
                converged = .true.
                exit 
            end if 
            if (any(isnan(res))) then 
                ! Probably out of bounds of the interpolator, exit
                converged = .false. 
                exit 
            end if 
            
            ! Compute update
            call magneticField%interp%Evaluate(xi, yi, 2, 0, d2Fdx2)
            call magneticField%interp%Evaluate(xi, yi, 1, 1, d2Fdxdy)
            call magneticField%interp%Evaluate(xi, yi, 0, 2, d2Fdy2)
            A(:, 1) = [d2Fdx2, d2Fdxdy] 
            A(:, 2) = [d2Fdxdy, d2Fdy2]
            call SolveDenseLinearSystemDI(A, res, d, flag)

            ! Check if solver converged, otherwise exit 
            if (flag /= 0) then 
                exit 
            end if 
            
            ! Apply updates
            xi = xi + d(1);
            yi = yi + d(2);
            
            ! Update counter
            it = it + 1
        end do 

        ! Unpack solution
        x = xi(1)
        y = yi(1)

    end subroutine 

    ! Eigenvalue computation for 2-by-2 matrices 
    function ComputeEigenvaluesSymmetric2by2Matrix(a, b, c) result(eig)
        
        ! Description
        !============
        ! Compute the eigenvalues of a 2-by-2 symmetric matrix, of which 
        ! the diagonal values should be given as the two first entries 
        ! and the off-diagonal value as the third entry. The 
        ! eigenvalue should be per definition real, since the matrix
        ! is symmetric. 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)        :: a, b, c
        real(R8)                    :: eig(1:2)

        ! Auxiliary
        real(R8)                    :: det 

        ! Compute
        !========
        ! Determinant
        det = sqrt((a+b)**2 - 4*(a*b-c**2))

        ! First eigenvalue
        eig(1) = ((a+b) + det)/2.0_R8
        eig(2) = ((a+b) - det)/2.0_R8

    end function 

    ! Coordinate trimmer
    subroutine DeleteCurveSegment(x, y, sr, loc, dloffset, keepstartpoint, &
        keependpoint)

        ! Description
        !============
        ! An auxiliary routine that deletes a segment of a curve, defined
        ! by its x, y coordinates. The segment that is deleted is defined 
        ! by 'loc', which can be 'start', 'end', or 'both'. 'sr' is then the
        ! length coordinate (curve length goes from 0 to number of segments) at which the 
        ! cutoff is defined. Typically, this length originates from an 
        ! intersection routine. an 'offset' can be defined (0 <= offset <= 1) 
        ! which gives the relative distance that should be kept between the
        ! intersection position and the next vertex. If 'loc' is 'both', 
        ! then offset and sr is expected to be a 1-by-2
        ! array. Finally, the boolean 'keependpoints' determines whether
        ! the start and end of the curve should be conserved or deleted.

        ! Declare variables
        !==================
        ! Arguments
        real(R8), allocatable, dimension(:), intent(inout)  :: x, y
        real(R8), dimension(:), intent(in)                  :: sr 
        character(*), intent(in)                            :: loc 
        real(R8), intent(in)                                :: dloffset(:) 
        logical, intent(in)                                 :: keepstartpoint 
        logical, intent(in)                                 :: keependpoint 

        ! Auxiliary
        real(R8), dimension(:), allocatable     :: offset
        real(R8), allocatable, dimension(:)     :: actualsr, &
            dlc, xint, yint
        integer(I8)                             :: ns
        integer(I8), allocatable, dimension(:)  :: keepind

        ! Loop
        integer(I8)                             :: i

        ! Delete segment
        !===============
        ! Initialize
        actualsr = sr 

        ! Project offset within bounds if necessary
        offset = dloffset
        !where (dloffset >= 1.0_R8) offset = 1.0_R8
        !where (dloffset <= 0.0_R8) offset = 0.0_R8
        !where (dloffset < 1.0_R8 .and. dloffset > 0.0_R8) offset = dloffset

        ! Compute length distribution
        dlc = [(i, i = 0, size(x)-1)]
        ns = size(dlc)-1

        ! Check case
        select case (loc)

        case ('start')

            ! Check sr
            if (sr(1) == real(ceiling(sr(1)), kind=R8)) then 
                ! Intersection exactly at end node, set segment index
                actualsr(1) = sr(1) + offset(1)
            else
                if (offset(1) == 1.0_R8) then
                    actualsr(1) = real(ceiling(sr(1)), kind=R8)
                else
                    actualsr(1) = sr(1) + (ceiling(sr(1)) - sr(1))*offset(1)
                end if 
            end if 

            ! Compute intersection coordinates
            call Interpolate1D(actualsr(1:1), xint, dlc, x)
            call Interpolate1D(actualsr(1:1), yint, dlc, y)

            ! Compute which points to keep 
            keepind = [(i, i = ceiling(actualsr(1))+1, ns+1)]

            ! Compute new coordinates
            if (keepstartpoint) then 
                x = [x(1), xint(1), x(keepind)]
                y = [y(1), yint(1), y(keepind)]
            else
                x = [xint(1), x(keepind)]
                y = [yint(1), y(keepind)]
            end if 

        case ('end')

            ! Check sr
            if (sr(1) == real(floor(sr(1)), kind=R8)) then 
                ! Intersection exactly at end node, set segment index
                actualsr(1) = sr(1) - offset(1)
            else
                actualsr(1) = sr(1) - (sr(1) - floor(sr(1)))*offset(1)
            end if 

            ! Compute intersection coordinates
            call Interpolate1D(actualsr(1:1), xint, dlc, x)
            call Interpolate1D(actualsr(1:1), yint, dlc, y)

            ! Compute which points to keep 
            keepind = [(i, i = 1, floor(actualsr(1)))]

            ! Compute new coordinates
            if (keependpoint) then 
                x = [x(keepind), xint, x(ns+1)]
                y = [y(keepind), yint, y(ns+1)]
            else
                x = [x(keepind), xint]
                y = [y(keepind), yint]
            end if 

        case ('both')

            ! Check sr
            if (sr(1) == real(ceiling(sr(1)), kind=R8)) then 
                ! Intersection exactly at end node, set segment index
                actualsr(1) = sr(1) + offset(1)
            else
                if (offset(1) == 1.0_R8) then
                    actualsr(1) = real(ceiling(sr(1)), kind=R8)
                else
                    actualsr(1) = sr(1) + (ceiling(sr(1)) - sr(1))*offset(1)
                end if 
            end if
            if (sr(2) == real(floor(sr(2)), kind=R8)) then 
                ! Intersection exactly at end node, set segment index
                actualsr(2) = sr(2) - offset(2)
            else
                actualsr(2) = sr(2) - (sr(2) - floor(sr(2)))*offset(2)
            end if 

            ! Compute intersection coordinates
            call Interpolate1D(actualsr, xint, dlc, x)
            call Interpolate1D(actualsr, yint, dlc, y)

            ! Compute which points to keep 
            keepind = [(i, i = ceiling(actualsr(1))+1, floor(actualsr(2)))]

            ! Compute new coordinates
            if (keepstartpoint .and. keependpoint) then 
                x = [x(1), xint(1), x(keepind), xint(2), x(ns+1)]
                y = [y(1), yint(1), y(keepind), yint(2), y(ns+1)]
            elseif (keepstartpoint) then 
                x = [x(1), xint(1), x(keepind), xint(2)]
                y = [y(1), yint(1), y(keepind), yint(2)]
            elseif (keependpoint) then 
                x = [xint(1), x(keepind), xint(2), x(ns+1)]
                y = [yint(1), y(keepind), yint(2), y(ns+1)]
            else
                x = [xint(1), x(keepind), xint(2)]
                y = [yint(1), y(keepind), yint(2)]
            end if 

        case default 

            call gdErrorHandler('DetermineCurveSegment: unknown location')

        end select

    end subroutine

    ! Tracer updater
    subroutine UpdateTracersFromTopomesh(topomesh, tracer, &
        magneticField, vessel, options)

        ! Description
        !============
        ! This routine updates the field/vessel tracers by adding the 
        ! topological mesh points as saddle points 

        ! Modules
        !========
        use mod_structured2Dgridding

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT), intent(in)          :: topomesh
        class(ContourTracerUDT), allocatable, intent(inout)  :: tracer
        type(TopomeshOptionsUDT), intent(in)    :: options
        type(MagneticFieldUDT), intent(in)      :: magneticField
        type(VesselUDT), intent(in)             :: vessel

        ! Auxiliary
        real(R8)                                :: dxfracmin, dyfracmin
        real(R8), allocatable, dimension(:)     :: xb, yb, &
            xg, yg, Vf, xgv, ygv, xps, yps
        real(R8), parameter                     :: emptyR8(0)= 0
        real(R8), allocatable, dimension(:)     :: xtp, ytp, Ftp
        integer(I8)                             :: ntp 
        integer(I8), allocatable, dimension(:)  :: IDs
        integer(I8), parameter                  :: emptyI8(0) = 0
        logical, allocatable, dimension(:)      :: includevert

        ! Grid
        !=====
        ! Initialize
        dxfracmin = 1e-4_R8 
        dyfracmin = 1e-4_R8

        ! Determine domain bounds based on vessel and magnetic field extent
        call vessel%plfvessel%ps%GetVertices(xps, yps)
        xb = [minval([xps, magneticField%interp%xgv]), maxval([xps, magneticField%interp%xgv])]
        yb = [minval([yps, magneticField%interp%ygv]), maxval([yps, magneticField%interp%ygv])]

        ! Construct refined grid based on tangency points and extrema
        includevert = (topomesh%vert%type == TMvertextp1ID) .or. &
            (topomesh%vert%type == TMvertextp2ID) .or. (topomesh%vert%type == TMvertexsaddleID)
        ! includevert = topomesh%vert%type == TMvertexsaddleID
        ntp = count(includevert)
        allocate(xtp(ntp), ytp(ntp), Ftp(ntp), IDs(ntp))
        xtp = pack(topomesh%vert%x, includevert)
        ytp = pack(topomesh%vert%y, includevert)
        Ftp = pack(topomesh%vert%fval, includevert)
        IDs = pack(topomesh%vert%ID, includevert)
        call ConstructRefined2DStructuredGrid(xg, yg, xgv, ygv, xb, yb, &
            options%vresx, options%vresy, xtp, ytp, 5, 5, dxfracmin, dyfracmin)  

        ! Evaluate magnetic field and vessel
        allocate(Vf(size(xg)))
        call magneticField%interp%Evaluate(xg, yg, 0, 0, Vf)
        call magneticField%interp%Evaluate(xtp, ytp, 0, 0, Ftp)

        ! Update the tracer 
        tracer = ConstructStructuredTracer(&
            reshape(Vf, [size(xgv), size(ygv)]), xgv, ygv, &
            xtp, ytp, Ftp, IDs, tracer%npmin, tracer%npmax, tracer%dl)

    end subroutine
    !------------------------------------------------------------------!
    !                               I/O                                !
    !------------------------------------------------------------------!
    ! Topological mesh data writer
    subroutine WriteTopologicalMesh(topomesh, filename)

        ! Description
        !============
        ! Routine to write topological mesh file in our own format. Also
        ! handy for reading in once constructed. The following data is
        ! written in a column-wise fashion:

        ! Vertex data
        !============
        ! 'vertices' 
        ! <vert%ntot> 
        ! 'ID, x, y, type, fsID, fval, BV'
        ! <ID, x, y, type, fsID, fval, BV (as zero or one)>
        ! 'vertex face pointer'
        ! <vertex%faceP(:, 1), vertex%faceP(:, 2)>
        ! 'vertex faces'
        ! <vertex%face>
        ! 'faces'
        ! <face%ntot> 
        ! 'ID, fsID, type, vert1, vert2, BF, nc'
        ! <ID, fsID, type, vert(:, 1), vert(:, 2), BF, face%x%size(), label>
        ! 'face <nf>' <repeated for each face including header>
        ! <x, y> 
        ! 'cells'
        ! <cell%ntot,  cells%nvert, cells%nface> 
        ! 'ID, cell polygon size>'
        ! <ID, nc>
        ! 'cell vertices'
        ! <cell%vert>
        ! 'cell vertex pointer'
        ! <cell%vertP>
        ! 'cell faces'
        ! <cell%face> 
        ! 'cell face pointer'
        ! <cell%faceP(:, 1), cell%faceP(:, 2)>
        ! 'cell polygon <nc>' (repeated nc times)
        ! <x, y> 
        ! 'tubes'
        ! <tubes%ntot,  tubes%nface, tubes%ncell> 
        ! 'tube faces'
        ! <tube%face>
        ! 'tube face pointer'
        ! <tube%faceP(:, 1), tube%faceP(:, 2)>
        ! 'tube cells'
        ! <tube%cell>
        ! 'tube cell pointer'
        ! <tube%cellP(:, 1), tube%cellP(:, 2)>
        
        ! Declare variables
        !==================
        ! Modules 
        use mod_plotter 
        use mod_specialchars, only : filesepchar

        ! Arguments
        class(TopomeshUDT)                      :: topomesh
        character(*), intent(in)                :: filename 

        ! Auxiliary
        integer                                 :: fu, BVval, cvsize, &
            cfsize
        integer(I8), allocatable, dimension(:)  :: tID
        real(R8), allocatable, dimension(:)     :: xf, yf, tfval
        character(:), allocatable               :: dir
        logical                                 :: isuptodate
        logical, allocatable, dimension(:)      :: BV, BF
        type(RealDynamicArrayUDT), allocatable  :: xcda(:), ycda(:)

        ! Loop
        integer(I8)                             :: i, j 

        ! Initialize
        !===========
        ! Unpack
        associate(&
            v       => topomesh%vert,   &
            f       => topomesh%face,   &
            c       => topomesh%cell,   &
            t       => topomesh%tube    &
        )

        ! Switch to check if full topomesh is up to date. If not, 
        ! certain (many) things will not be written out (cells, 
        ! some connectivity data, ...)
        isuptodate = .true.

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
        if (.not. allocated(v%BV)) then 
            allocate(BV(v%ntot))
            BV = .false.
        elseif (size(v%BV) /= v%ntot) then 
            allocate(BV(v%ntot))
            BV = .false.
        else
            BV = v%BV
        end if 
        if (.not. allocated(f%BF)) then 
            allocate(BF(f%ntot))
            BF = .false.
        elseif (size(f%BF) /= f%ntot) then 
            allocate(BF(f%ntot))
            BF = .false.
        else 
            BF = f%BF
        end if

        ! Write vertex data
        !==================
        ! Number of vertices
        write (fu, *) 'vertices'
        write (fu, *) v%ntot

        ! Basic vertex data
        write (fu, *) 'ID, x, y, type, fsID, fval, BV'
        do i = 1, v%ntot
            if (BV(i)) then 
                BVval = 1
            else 
                BVval = 0
            end if
            write (fu, *) v%ID(i), v%x(i), v%y(i), v%type(i), v%fsID(i), v%fval(i), BVval
        end do 

        ! Vertex face data
        if (allocated(v%faceP) .and. allocated(v%face)) then 
            if (size(v%faceP, 1) == v%ntot) then ! may not be up to date
                write (fu, *) 'vertex face pointer'
                do i = 1, v%ntot
                    write (fu, *) v%faceP(i, 1), v%faceP(i, 2)
                end do 

                write (fu, *) 'vertex facelist'
                do i = 1, size(v%face)
                    write (fu, *) v%face(i)
                end do 
            else
                isuptodate = .false.
            end if 
        else
            isuptodate = .false.
        end if 

        ! Write face data
        !================
        ! Number of faces
        write (fu, *) 'faces'
        write (fu, *) f%ntot 

        ! Basic face data
        write (fu, *) 'ID, fsID, type, vert1, vert2, BF, nc'
        do i = 1, f%ntot 
            if (BF(i)) then 
                BVval = 1
            else 
                BVval = 0
            end if 
            write (fu, *) f%ID(i), f%fsID(i), f%type(i), f%vert(i, 1), &
                f%vert(i, 2), BVval, f%x(i)%size()
        end do 
        
        ! Face coordinates
        write (fu, *) 'face coordinates'
        do i = 1, f%ntot 
            ! Write header
            write (fu, *) 'face ', i 

            ! Write coordinates
            xf = f%x(i)%Get()
            yf = f%y(i)%Get()
            do j = 1, size(xf)
                write(fu, *) xf(j), yf(j)
            end do 
        end do 

        ! Write cell data
        !================
        ! Number of cells, number of cell vertices, number of cell faces
        if (allocated(c%vert) .and. isuptodate) then 
            if ((size(c%vert) == c%ntotv) .and. .not. (any(c%vert == 0))) then 
                cvsize = size(c%vert)
            else
                cvsize = 0
                isuptodate = .false.
            end if 
        else
            cvsize = 0
            isuptodate = .false.
        end if 
        if (allocated(c%face) .and. isuptodate) then 
            if (size(c%face) == c%ntotf) then 
                cfsize = size(c%face)
            else
                cfsize = 0
                isuptodate = .false.
            end if 
        else
            cfsize = 0
            isuptodate = .false.
        end if 
        write (fu, *) 'cells'
        if (isuptodate) then 
            write (fu, *) c%ntot, cvsize, cfsize
        else
            write (fu, *) 0_I8, cvsize, cfsize
        end if 

        ! Print additional data
        if (isuptodate) then 
            ! Compute cell polygon data
            allocate(xcda(c%ntot), ycda(c%ntot))
            do i = 1, c%ntot
                call ConstructTopologicalMeshCellPolygon(topomesh, i, xcda(i), ycda(i))
            end do 

            ! Cell data
            write (fu, *) 'ID, nc'
            do i = 1, c%ntot 
                write (fu, *) i, xcda(i)%Size()
            end do 

            ! Cell vertices
            write (fu, *) 'cell vertexlist'
            do i = 1, cvsize
                write (fu, *) c%vert(i)
            end do 
            write (fu, *) 'cell vertex pointer'
            do i = 1, c%ntot 
                write (fu, *) c%vertP(i, 1), c%vertP(i, 2)
            end do 

            ! Cell faces
            write (fu, *) 'cell facelist'
            do i = 1, cfsize
                write (fu, *) c%face(i)
            end do 
            write (fu, *) 'cell face pointer'
            do i = 1, c%ntot 
                write (fu, *) c%faceP(i, 1), c%faceP(i, 2)
            end do 

            ! Cell polygons
            write (fu, *) 'cell polygons'
            do i = 1, c%ntot 
                ! Write header
                write (fu, *) 'cell ', i 

                ! Write coordinates
                xf = xcda(i)%Get()
                yf = ycda(i)%Get()
                do j = 1, size(xf)
                    write(fu, *) xf(j), yf(j)
                end do 
            end do 
        end if 

        ! Write flux surface data
        !========================
        if (isuptodate) then 
            ! Header and sizes
            write (fu, *) 'flux surfaces'
            write (fu, *) topomesh%nfs 
            
            ! Data
            write (fu, *) 'ID, fval'
            tID = topomesh%fsID%Get()
            tfval = topomesh%fsfval%Get()
            do i = 1, size(tID)
                write (fu, *) tID(i), tfval(i)
            end do 
        else
            ! Header and sizes
            write (fu, *) 'flux surfaces'
            write (fu, *) 0_I8 
        end if 

        ! Write tube data
        !================
        if (isuptodate) then 
            ! Header and sizes
            write (fu, *) 'tubes'
            write (fu, *) t%ntot, t%nface, t%ncell 

            ! Tube faces
            write (fu, *) 'tube facelist'
            do i = 1, t%nface
                write (fu, *) t%face(i)
            end do 
            write (fu, *) 'tube face pointer'
            do i = 1, t%ntot 
                write (fu, *) t%faceP(i, 1), t%faceP(i, 2)
            end do 

            ! Tube cells
            write (fu, *) 'tube cells'
            do i = 1, t%ncell
                write (fu, *) t%cell(i)
            end do 
            write (fu, *) 'tube cell pointer'
            do i = 1, t%ntot 
                write (fu, *) t%cellP(i, 1), t%cellP(i, 2)
            end do 
        else
            ! Header and sizes
            write (fu, *) 'tubes'
            write (fu, *) 0_I8, 0_I8, 0_I8
        end if 

        ! Housekeeping
        !=============
        ! Deallocate again

        ! Others
        end associate
        close(fu)

    end subroutine

    ! Topological mesh data reader 
    subroutine ReadTopologicalMesh(topomesh, filepath)

        ! Description
        !============
        ! This routine reads in the topological mesh data of an 
        ! existing topological mesh (which is assumed to be correctly
        ! set up) and adds any additional interconnection data that is 
        ! not directly available in the file itself.

        ! Declare variables
        !==================
        ! Modules 
        use mod_readwrite
        use mod_inputfileparser

        ! Arguments
        class(TopomeshUDT), intent(out)         :: topomesh
        character(*), intent(in)                :: filepath 

        ! Auxiliary
        integer                                 :: fu, BVval, &
            nfc, ncc
        real(R8), allocatable, dimension(:)     :: xf, yf, tfval 
        integer(I8), allocatable, dimension(:)  :: tID
        character(:), allocatable               :: thisline
        logical                                 :: reachedeof

        ! Loop
        integer(I8)                             :: i, j 

        ! Initialize
        !===========
        ! Unpack
        associate(&
            v       => topomesh%vert,   &
            f       => topomesh%face,   &
            c       => topomesh%cell,   &
            t       => topomesh%tube    &
        )

        ! Open file
        open (action='read', file=trim(filepath), newunit=fu, &
             status='old')

        ! Initialize
        call topomesh%Initialize()

        ! Vertices
        !=========
        ! Read until header found
        call ReadUntilFound(fu, 'vertices', reachedeof)
        if (reachedeof) then 
            call gdErrorHandler('ReadTopologicalMesh: could not find vertices')
        end if 

        ! Read number of vertices
        read(fu, *) v%ntot

        ! Initialize vertices 
        call v%Initialize(nv=v%ntot)

        ! Basic vertex data
        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, v%ntot
            read(fu, *) v%ID(i), v%x(i), v%y(i), v%type(i), v%fsID(i), v%fval(i), BVval
            if (BVval == 1) then 
                v%BV(i) = .true.
            else
                v%BV(i) = .false.
            end if 
        end do 

        ! Vertex face data
        call ReadSingleLine(fu, thisline, reachedeof)
        if (allocated(v%faceP)) then 
            deallocate(v%faceP)
        end if
        allocate(v%faceP(v%ntot, 2))
        v%faceP = 0
        do i = 1, v%ntot 
            read(fu, *) v%faceP(i, 1), v%faceP(i, 2)
        end do 
        if (allocated(v%face)) then 
            deallocate(v%face)
        end if
        allocate(v%face(sum(v%faceP(:, 2))))
        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, size(v%face)
            read(fu, *) v%face(i)
        end do 

        ! Faces
        !======
        ! Read until header found
        call ReadUntilFound(fu, 'faces', reachedeof)
        if (reachedeof) then 
            call gdErrorHandler('ReadTopologicalMesh: could not find faces')
        end if 

        ! Read number of faces
        read(fu, *) f%ntot

        ! Initialize vertices 
        call f%Initialize(nf=f%ntot)

        ! Read basic data
        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, f%ntot
            ! Read
            read(fu, *) f%ID(i), f%fsID(i), f%type(i), f%vert(i, 1), &
                f%vert(i, 2), BVval, nfc
            if (BVval == 1) then 
                f%BF(i) = .true.
            else
                f%BF(i) = .false.
            end if 

            ! Initialize
            f%x(i) = ConstructRealDynamicArray(spread(0.0_R8, 1, nfc))
            f%y(i) = ConstructRealDynamicArray(spread(0.0_R8, 1, nfc))
        end do 

        ! Read face coordinates
        call ReadSingleLine(fu, thisline, reachedeof) ! skip the header
        do i = 1, f%ntot 
            ! Read header
            call ReadSingleLine(fu, thisline, reachedeof) 

            ! Read coordinates
            xf = f%x(i)%Get()
            yf = f%y(i)%Get()
            do j = 1, f%x(i)%Size()
                read(fu, *) xf(j), yf(j)
            end do 
            call f%x(i)%Set(xf)
            call f%y(i)%Set(yf)
            call f%pol(i)%Construct(xf, yf)
        end do 

        ! Cells
        !======
        ! Read until header found
        call ReadUntilFound(fu, 'cells', reachedeof)
        if (reachedeof) then 
            call gdErrorHandler('ReadTopologicalMesh: could not find cells')
        end if 

        ! Read number of faces
        read(fu, *) c%ntot, c%ntotv, c%ntotf 

        ! Initialize vertices 
        call c%Initialize(nc=c%ntot, ncv=c%ntotv, ncf=c%ntotf)

        ! Basic data
        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, c%ntot 
            read(fu, *) c%ID(i), ncc
        end do 

        ! Cell vertices
        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, c%ntotv
            read (fu, *) c%vert(i)
        end do 

        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, c%ntot
            read(fu, *) c%vertP(i, 1), c%vertP(i, 2)
        end do 

        ! Cell faces
        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, c%ntotf
            read (fu, *) c%face(i)
        end do 

        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, c%ntot
            read(fu, *) c%faceP(i, 1), c%faceP(i, 2)
        end do 

        ! Flux surface data
        !==================
        ! Read until header found
        call ReadUntilFound(fu, 'flux surfaces', reachedeof)
        if (reachedeof) then 
            call gdErrorHandler('ReadTopologicalMesh: could not find flux surfaces')
        end if 

        ! Read number of surfaces
        read(fu, *) topomesh%nfs 

        ! Allocate
        allocate(tfval(topomesh%nfs), tID(topomesh%nfs))

        ! Read 
        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, topomesh%nfs
            read(fu, *) tID(i), tfval(i)
        end do 
        topomesh%fsID = ConstructIntegerDynamicArray(tID)
        topomesh%fsfval = ConstructRealDynamicArray(tfval)

        ! Add interconnection data
        !=========================
        call AddTopologicalMeshData(topomesh)
        call AddTopologicalMeshInterconnectionData(topomesh)

        ! Housekeeping
        !=============
        ! Deallocate again

        ! Others
        end associate
        close(fu)


    end subroutine


end module