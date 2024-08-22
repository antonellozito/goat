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
    use mod_definitions, only: TMvertexbndID, TMvertexmaxID, &
        TMvertexminID, TMvertexsaddleID, TMvertextp1ID, &
        TMvertextp2ID, TMfacepolID, TMfaceradID, TMfacebndID, &
        TMvertexsplitID
    use mod_linearsolverinterface, only: SolveDenseLinearSystemDI
    use goatmod_types, only : magneticFieldUDT, VesselUDT
    use goatmod_userinput, only : TopomeshOptionsUDT
    implicit none
    private 
    public :: TopomeshUDT, ConstructTopologicalMesh, TraceExtrema2D, &
        TraceTangencyPoints2D

    ! Module parameters
    real(R8), parameter, private        :: tprelfieldtol = 1e-2 ! relative field tolerance under which extrema are removed
    real(R8), parameter, private        :: disttol = 1e-12 ! distance tolerance

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
            flags(:), type(:), faceP(:, :), cellP(:, :)
        real(R8), allocatable           :: x(:), y(:), fval(:)
        logical, allocatable            :: BV(:)

    contains 

        ! Initializer
        procedure :: Initialize     => InitializeTopologicalMeshVertex

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

        integer(I8)                     :: ntot ! total number of vertices
        integer(I8), allocatable        :: ID(:), vert(:, :), cell(:), &
            label(:), fsID(:), type(:), cellP(:, :)
        type(RealDynamicArrayUDT), allocatable  :: x(:), y(:)
        real(R8), allocatable           :: fval(:)
        type(PolygonUDT), allocatable   :: pol(:)
        logical, allocatable            :: BF(:)

    contains 
    
        ! Initializer
        procedure :: Initialize     => InitializeTopologicalMeshFace

    end type

    ! Topological cell
    type :: TopomeshCellUDT

        ! Description
        !============
        ! Topological mesh cell

        integer(I8)                     :: ntot ! total number of vertices
        integer(I8), allocatable        :: ID(:), vert(:), vertP(:, :), &
            face(:), faceP(:, :), flags(:)
    contains 

        ! Initializer
        procedure :: Initialize     => InitializeTopologicalMeshCell

    end type

    ! General topological mesh type
    type :: TopomeshUDT 

        type(TopomeshVertUDT)   :: vert 
        type(TopomeshFaceUDT)   :: face 
        type(TopomeshCellUDT)   :: cell 
        integer(I8)             :: nfs 
    contains 

        ! Initialize
        procedure :: Initialize         => InitializeTopologicalMesh 
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
    function ConstructTopologicalMesh(vessel, magneticField, options) result(topomesh)

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

        ! Auxiliary
        class(ContourTracerUDT), allocatable    :: vesseltracer, fieldtracer
        real(R8), allocatable, dimension(:)     :: xb, yb, xps, &
            yps, xg, yg, Vf, Vv, xgv, ygv
        real(R8), parameter                     :: emptyR8(0)= 0
        real(R8), allocatable, dimension(:)     :: xtp, ytp, Ftp, &
            xe, ye, fe
        integer(I8)                             :: nv 
        integer(I8), allocatable, dimension(:)  :: typetp, typee

        ! Loop
        integer(I8)                             :: i, j

        ! Initialize
        !===========
        ! Initialize topomesh 
        call topomesh%Initialize()

        ! Determine domain bounds based on vessel and magnetic field extent
        call vessel%plfvessel%ps%GetVertices(xps, yps)
        xb = [minval([xps, magneticField%interp%xgv]), maxval([xps, magneticField%interp%xgv])]
        yb = [minval([yps, magneticField%interp%ygv]), maxval([yps, magneticField%interp%ygv])]

        ! Construct a 2D structured grid for tracing (may be extended
        ! in the future for different grid types)
        nv = options%vresx*options%vresy
        allocate(xg(nv), yg(nv), Vf(nv), Vv(nv), xgv(options%vresx), &
            ygv(options%vresy))
        call Construct2DStructuredUniformGrid(xg, yg, xgv, ygv, xb, yb, &
            options%vresx,  options%vresy, 0.0_R8, 0.0_R8)

        ! Evaluate the field and vessel values
        call vessel%plfvessel%Evaluate(xg, yg, 0, 0, Vv)
        call magneticField%interp%Evaluate(xg, yg, 0, 0, Vf)

        ! Construct the vessel and magnetic field tracers (currently in 
        ! structured format)
        vesseltracer = ConstructStructuredTracer(&
            reshape(Vv, [options%vresx, options%vresy]), xgv, ygv, &
            emptyR8, emptyR8, emptyR8)
        fieldtracer = ConstructStructuredTracer(&
            reshape(Vf, [options%vresx, options%vresy]), xgv, ygv, &
            emptyR8, emptyR8, emptyR8)

        ! Compute tangency points & extrema
        !==================================
        call TraceTangencyPoints2D(xtp, ytp, typetp, Ftp, &
            vesseltracer, magneticField)
        call TraceExtrema2D(xe, ye, fe, typee, fieldtracer, magneticField, &
            options%fdonewton)
        

        ! Update the tracers
        !===================
        ! Construct refined grid based on tangency points and extrema
        call ConstructRefined2DStructuredGrid(xg, yg, xgv, ygv, xb, yb, &
            options%vresx, options%vresy, [xtp, xe], [ytp, ye], 5, 5)  

        ! Evaluate magnetic field and vessel
        deallocate(Vv, Vf)
        allocate(Vv(size(xg)), Vf(size(xg)))
        call vessel%plfvessel%Evaluate(xg, yg, 0, 0, Vv)
        call magneticField%interp%Evaluate(xg, yg, 0, 0, Vf)

        ! Construct the vessel and magnetic field tracers (currently in 
        ! structured format)
        vesseltracer = ConstructStructuredTracer(&
            reshape(Vv, [size(xgv), size(ygv)]), xgv, ygv, &
            xtp, ytp, Ftp*0.0_R8) ! tangency points lie on vessel so zero value
        fieldtracer = ConstructStructuredTracer(&
            reshape(Vf, [size(xgv), size(ygv)]), xgv, ygv, &
            xtp, ytp, Ftp)

        ! Add extrema & related boundaries
        !=================================
        ! Extrema
        call AddTopologicalMeshExtrema(topomesh, fieldtracer, &
            magneticField, options)

        ! Tangency points
        call AddTopologicalMeshTangencyPoints(topomesh, xtp, ytp, ftp, &
            typetp, vesseltracer)

        ! Remove parts that do not lie inside the vessel
        call TrimTopologicalMesh(topomesh, magneticField, vessel)

        ! Process points
        !===============
        ! Set extrema with nearly identical values to be the same value
        do i = 1, topomesh%vert%ntot-1
            do j = i+1, topomesh%vert%ntot 
                if (abs(topomesh%vert%fval(i)-topomesh%vert%fval(j)) < options%ffieldtol) then 
                    topomesh%vert%fval(j) = topomesh%vert%fval(i)
                end if 
            end do 
        end do 

        ! Compute and add intersections with contours
        !============================================
        ! Necessary contours
        call AddTopologicalMeshContours(topomesh, magneticField, vessel, &
            fieldtracer, vesseltracer, options)

        ! Core boundaries? 
        if (options%addcoreboundaries) then 
            call AddTopologicalMeshCoreBoundaries(topomesh, magneticField, &
                vessel, fieldtracer, options)
        end if 

        ! Compute additional interconnnection data
        !=========================================
        ! Vertex faces
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Data 
        call AddTopologicalMeshData(topomesh)

        ! Add cells
        call AddTopologicalMeshCells(topomesh)

        ! Remove parts if desired
        if (options%removecoreregions) then 
            call RemoveTopologicalMeshCoreRegions(topomesh)
        end if 

        ! Compute interconnection data
        call AddTopologicalMeshInterconnectionData(topomesh)

        ! Write
        !======
        call WriteTopologicalMesh(topomesh, 'topomesh')

    end function

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

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT), intent(inout)       :: topomesh
        class(ContourTracerUDT), intent(inout)  :: fieldtracer 
        type(MagneticFieldUDT), intent(in)      :: magneticField 
        type(TopomeshOptionsUDT), intent(in)    :: options 

        ! Auxiliary
        integer(I8)                             :: ngp, nfxc, nfyc, &
            nx
        integer(I8), allocatable, dimension(:)  :: ts1, ts2, tt, typee, &
            vf1, vf2, teid, tsid, sortind
        real(R8)                                :: tempx, tempy
        real(R8), allocatable, dimension(:)     :: xg, yg, f, fx, fy, &
            tx, ty, thiseig, tf, tfxx, tfxy, tfyy, xe, ye, fe, tsr1, tsr2, &
            tsrid
        logical                                 :: conv, donewton 
        type(RealDynamicArrayUDT)               :: xc, yc, fc
        type(RealDynamicArrayUDT), allocatable  :: xfrda(:), yfrda(:), &
            fxpsrid(:), fypsrid(:) 
        type(IntegerDynamicArrayUDT), allocatable   :: fxpeid(:), &
            fxpsid(:), fypeid(:), fypsid(:)
        type(IntegerDynamicArrayUDT)            :: tc
        type(ContourUDT), allocatable           :: fxc(:), fyc(:)
        type(PolygonUDT), allocatable           :: fxp(:), fyp(:)
        class(ContourTracerUDT), allocatable    :: fxtracer, fytracer

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
            deallocate(fxtracer%xs, fxtracer%ys, fxtracer%vs, fxtracer%order)
        end if 
        allocate(fxtracer%xs(0), fxtracer%ys(0), fxtracer%vs(0), fxtracer%order(0))
        if (allocated(fytracer%xs)) then 
            deallocate(fytracer%xs, fytracer%ys, fytracer%vs, fytracer%order)
        end if 
        allocate(fytracer%xs(0), fytracer%ys(0), fytracer%vs(0), fytracer%order(0))
        call fxtracer%SetValues(fx)
        call fytracer%SetValues(fy)
        

        ! Trace 
        fxc = fxtracer%TraceContours([0.0_R8])
        fyc = fytracer%TraceContours([0.0_R8])

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
                    ec = ec + 1

                    ! Store intersection data
                    call fxpeid(i)%Append(ec)
                    call fxpsid(i)%Append(ts1(k))
                    call fxpsrid(i)%Append(tsr1(k))
                    call fypeid(j)%Append(ec)
                    call fypsid(j)%Append(ts2(k))
                    call fypsrid(j)%Append(tsr2(k))
                end do 

                ! Append
                call xc%Append(tx)
                call yc%Append(ty)
                call fc%Append(tf)
                call tc%Append(tt)

                ! Housekeeping
                deallocate(tt, tf, tfxx, tfyy, tfxy)
            end do 
        end do
        
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
            call AddTopologicalMeshVertex(topomesh, xe(i), ye(i), fe(i), &
                typee(i))
        end do

        ! Reconstruct faces
        !==================
        ! Loop over fxp
        do i = 1, nfxc 
            ! Extract the faces from this polygon
            teid = fxpeid(i)%Get()
            tsid = fxpsid(i)%Get()
            tsrid = fxpsrid(i)%Get()
            allocate(sortind(size(tsrid)))
            call Sort(tsrid, ind=sortind)
            tsid = tsid(sortind)
            teid = teid(sortind)
            deallocate(sortind)
            call ExtractTopologicalFacesFromPolygon(fxp(i), teid, &
                tsid, xe, ye, vf1, vf2, xfrda, yfrda)

            ! Add to faces
            do j = 1, size(vf1)
                call AddTopologicalMeshFace(topomesh, [vf1(j), vf2(j)], &
                    xfrda(j), yfrda(j), TMfaceradID, 0)
            end do 
        end do 

        ! Loop over fyp
        do i = 1, nfyc 
            ! Extract the faces from this polygon
            teid = fypeid(i)%Get()
            tsid = fypsid(i)%Get()
            tsrid = fypsrid(i)%Get()
            allocate(sortind(size(tsrid)))
            call Sort(tsrid, ind=sortind)
            tsid = tsid(sortind)
            teid = teid(sortind)
            deallocate(sortind)
            call ExtractTopologicalFacesFromPolygon(fyp(i), teid, &
                tsid, xe, ye, vf1, vf2, xfrda, yfrda)

            ! Add to faces
            do j = 1, size(vf1)
                call AddTopologicalMeshFace(topomesh, [vf1(j), vf2(j)], &
                    xfrda(j), yfrda(j), TMfaceradID, 0)
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
                tpf(i), ttype)
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
                call AddTopologicalMeshFace(topomesh, facevert, tx, ty, ftype, fsID)
            end do 
        end do



    end subroutine

    ! Necessary contours
    subroutine AddTopologicalMeshContours(topomesh, magneticField, vessel, &
        fieldtracer, bndtracer, options)

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

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        type(MagneticFieldUDT),intent(in)       :: magneticField 
        type(VesselUDT), intent(in)             :: vessel
        class(ContourTracerUDT), intent(inout)  :: fieldtracer, bndtracer 
        type(TopomeshOptionsUDT), intent(in)    :: options 

        ! Auxiliary
        real(R8)                                :: dist
        real(R8), allocatable, dimension(:)     :: pspx, pspy, pspf, &
            xout, yout, iout, jout, s1x, s1y, vesselval, s2x, s2y, &
            sortedI
        integer(I8)                             :: npsp, sizeI, sizetck, &
            ntc
        integer(I8), allocatable, dimension(:)  :: psptype, &
            sortind, pspID, delind, allcurvetypes, allfsIDs, &
            vindI, vindJ
        logical                                 :: deletefirstpart, &
            deletestart, deleteend
        logical, allocatable, dimension(:)      :: tracepoints, keepind, &
            isinpolyg

        type(ContourUDt)                        :: tempc
        type(ContourUDT), allocatable           :: tc(:), allc(:)
        type(PolygonUDT), allocatable           :: tcp(:)
        type(IntegerDynamicArrayUDT)            :: curvetypes, fsIDs, &
            intI 
        type(RealDynamicArrayUDT)               :: realI 

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
        pspx = topomesh%vert%x 
        pspy = topomesh%vert%y 
        pspf = topomesh%vert%fval 
        psptype = topomesh%vert%type 
        pspID = topomesh%vert%ID
        npsp = size(pspx)

        ! Check for saddle points with exactly the same field value - don't trace
        ! multiple times
        allocate(tracepoints(npsp))
        tracepoints = .true.
        do i = 1, npsp-1
            do j = i+1, npsp
                if ((pspf(i) == pspf(j)) .and. (psptype(i) == 2)) then 
                    tracepoints(j) = .false.
                end if 
            end do 
        end do

        ! Add saddle points to field tracer
        fieldtracer%xs = pspx 
        fieldtracer%ys = pspy
        fieldtracer%vs = pspf
        fieldtracer%order = 0*psptype   

        ! Trace contours
        !===============
        ! Add saddle point contours
        do i = 1, npsp
            if ((psptype(i) == 2) .and. tracepoints(i)) then 
                ! Trace
                tc = fieldtracer%TraceContours([pspx(i)], [pspy(i)])
                
                ! Process
                call CleanContours(tc)

                ! Add
                allc = [allc, tc]
                call curvetypes%Append([(k, k = 1, size(tc))]*4_I8)
                
                ! Add flux surface ID
                nfs = nfs + 1
                call fsIDs%Append([(k, k = 1, size(tc))]*nfs)
            end if 
        end do

        ! Add tangency point contours
        do i = 1, npsp
            if ((psptype(i) == 5) .and. tracepoints(i)) then 
                ! Trace
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
                    if (dist >= disttol) then 
                        print *, 'AddTopologicalMEshContours: tangency ' // & 
                            'contour segment found that does not start ' // & 
                            'in given tangency point. Removing...'
                        keepind(k) = .false.
                    end if
                end do 
                
                ! Remove
                tc = pack(tc, keepind)
                deallocate(keepind)

                ! Compute intersections with boundary faces
                ntc = size(tc)
                allocate(tcp(ntc))
                do k = 1, ntc
                    ! Grow dynamically because I'm lazy
                    intI = ConstructIntegerDynamicArray()
                    realI = ConstructRealDynamicArray()
                    deletefirstpart = .false.
                    
                    ! Construct polygon
                    call tcp(k)%Construct(tc(k)%x, tc(k)%y)

                    ! Loop over all faces
                    do j = 1, topomesh%face%ntot
                        if (topomesh%face%type(j) == 3) then 
                            ! If we have an intersection with a face that has this saddle point
                            ! as one of its vertices, we must remove parts of the intersection
                            ! (should always skip one boundary segment for intersection, since
                            ! adjacent points should be other tangency points right now
                            ! of type 4)
                            
                            ! We need to have a different approach for closed
                            ! contours!
                            
                            ! Compute intersections
                            call PolygonIntersections(tcp(k), topomesh%face%pol(j), &
                                xout, yout, vindI, vindJ, iout, jout)
                            
                            ! Check
                            if (size(xout) > 0) then 
                                ! Sort based on continuous polygon coordinate
                                allocate(sortind(size(jout)))
                                call Sort(iout, ind=sortind)
                                jout = jout(sortind)
                                vindI = vindI(sortind)
                                vindJ = vindJ(sortind)
                                deallocate(sortind)

                                ! Check
                                if (any(topomesh%face%vert(j, :) == pspID(i))) then 
                                    ! First intersection should be at starting point
                                    if (iout(1) > 0) then 
                                        call gdErrorHandler( & 
                                            'AddTopologicalMeshContours: ' // & 
                                            ' expected intersection in tangency ' // &
                                            'point but did not find it')
                                    end if
                                    
                                    ! Delete first intersection
                                    allocate(keepind(size(vindI)))
                                    keepind = .true. 
                                    keepind(1) = .false. 
                                    if (tcp(k)%isclosed) then 
                                        ! Also delete last intersection
                                        keepind(size(vindI)) = .false. 
                                    end if
                                    vindI = pack(vindI, keepind)
                                    vindJ = pack(vindJ, keepind)
                                    iout = pack(iout, keepind)
                                    jout = pack(jout, keepind)
                                    
                                    ! Check others
                                    if (size(vindI) > 0) then 
                                        
                                        ! Display a message
                                        print *, 'AddTopologicalMeshContours: ' // & 
                                            'tangency point contour of vertex nr ', &
                                            i, ' intersects prematurely with boundary, ' // & 
                                            ' adjusting contour and boundary segment'
                                        print *, 'AddTopologicalMeshContours: ' // & 
                                            'found ', size(vindI), &
                                            ' additional intersections, removing...'
                                        
                                        ! Set deletion to true
                                        deletefirstpart = .true.
                                        
                                        ! Issue warning: code not yet verified for
                                        ! more than one intersection (honestly
                                        ! don't think this should happen but you
                                        ! never know)
                                        if (size(vindI) > 1) then 
                                            print *, 'AddTopologicalMeshContours: ' // & 
                                                'code not yet verified ' // & 
                                                'for multiple intersections'
                                        end if 
                                        
                                        ! Check which part to delete from face data
                                        if (topomesh%face%vert(j, 1) == pspID(i)) then 
                                            delind = [1, (cc, cc = vindJ(size(vindJ))+1, topomesh%face%x(j)%Size() )]
                                        else
                                            delind = [(cc, cc = 1, vindJ(1)), topomesh%face%x(j)%Size()]
                                        end if 

                                        ! Delete
                                        call topomesh%face%x(j)%Remove(delind)
                                        call topomesh%face%y(j)%Remove(delind)
                                        
                                        ! Reconvert to polygon
                                        call topomesh%face%pol(j)%Construct(&
                                            topomesh%face%x(j)%Get(), topomesh%face%y(j)%Get())
                                    end if
                                end if
                            end if
                            
                            ! Add to intersections
                            if (size(vindI) > 0) then 
                                call realI%Append(iout(size(iout)))
                                call intI%Append(vindI(size(vindI)))
                            end if 
                        end if
                    end do
                    
                    ! Remove contour parts that lie outside of the vessel
                    sortedI = realI%Get()
                    vindI = intI%Get()
                    allocate(sortind(size(vindI)))
                    call Sort(sortedI, ind=sortind) ! Normally, first one will always be one to be deleted
                    vindI = vindI(sortind)
                    deallocate(sortind)

                    ! Trim contour data
                    sizeI = size(sortedI)
                    sizetck = sizetck
                    if (tcp(k)%isclosed) then 
                        ! If no other intersections than end points ->  continue
                        if (sizeI > 0) then 
                            
                            ! Need to determine how to trim closed contour.
                            ! going from both sides, we need to keep the first
                            ! segment inside of the vessel. We check using 
                            ! the vessel plf value, 
                            ! so this may be inaccurate when number of points is
                            ! too low per segment...
                            
                            ! Check first part
                            s1x = tc(k)%x(2:vindI(1))
                            s1y = tc(k)%y(2:vindI(1))
                            if (size(s1x) <= 3) then 
                                ! Issue warning
                                print *, 'AddTopologicalMeshContour: ' // & 
                                    'low amount of vertices for closed ' // & 
                                    'polygon tangency segment, may not ' // & 
                                    'be able to correctly identify segments'
                            end if
                            
                            ! Evaluate
                            allocate(vesselval(size(s1x)))
                            call vessel%plfvessel%Evaluate(s1x, s1y, 0, 0, vesselval)

                            ! Check if in vessel
                            isinpolyg = vesselval <= 0
                            if (all(isinpolyg)) then 
                                deletestart = .false.
                            else
                                deletestart = .true. 
                            end if 

                            ! Housekeeping
                            deallocate(vesselval)
                            
                            ! Check last part
                            s2x = tc(k)%x(vindI(sizeI)+1:sizetck-1)
                            s2y = tc(k)%y(vindI(sizeI)+1:sizetck-1)
                            if (size(s2x) <= 3) then 
                                ! Issue warning
                                print *, 'AddTopologicalMeshContour: ' // &
                                    'low amount of vertices for closed ' // &
                                    'polygon tangency segment, may not be ' // &
                                    'able to correctly identify segments'
                            end if 
                            
                            ! Evaluate
                            allocate(vesselval(size(s2x)))
                            call vessel%plfvessel%Evaluate(s2x, s2y, 0, 0, vesselval)

                            ! Check if in vessel
                            isinpolyg = vesselval <= 0
                            if (all(isinpolyg)) then 
                                deleteend = .false.
                            else
                                deleteend = .true.
                            end if
                            
                            ! Trim and add
                            if (sizeI == 1) then 
                                if (deletestart .and. .not. deleteend) then 
                                    s1x = [tc(k)%x(sizetck:vindI(1)+1:-1), &
                                        tc(k)%x(1)]
                                    s1y = [tc(k)%y(sizetck:vindI(1)+1:-1), &
                                        tc(k)%y(1)]
                                elseif (deleteend .and. .not. deletestart) then 
                                    s1x = [tc(k)%y(1:vindI(1)), &
                                        tc(k)%x(sizetck)]
                                    s1y = [tc(k)%y(1:vindI(1)), &
                                        tc(k)%y(sizetck)]
                                else
                                    ! This shouldn't be happening
                                    call gdErrorHandler('AddTopologicalMeshContours: ' // & 
                                        'unknown error when checking closed ' // & 
                                        'tangency contour (1 intersection), this is likely ' // & 
                                        'a bug') 
                                end if 
                                tc(k)%x = s1x;
                                tc(k)%y = s1y;
                            elseif (sizeI == 2) then 
                                ! Not yet properly verified, print warning
                                print *, 'AddTopologicalMeshContour: ' // & 
                                    'multiple intersections for closed ' // &
                                    'tangency contour detected. This part ' // &
                                    ' of the routine is not yet verified, ' // & 
                                    ' proceed with caution'
                                
                                ! Either both start and end to be removed, or both
                                ! to be added
                                if (deletestart .and. deleteend) then 
                                    s1x = tc(k)%x(vindI(1):vindI(2)+1)
                                    s1y = tc(k)%y(vindI(1):vindI(2)+1)
                                    tc(k)%x = s1x
                                    tc(k)%y = s1y
                                    tc(k)%startsaddle = 0
                                    tc(k)%endsaddle = 0
                                elseif (.not. deletestart .and. .not. deleteend) then 
                                    s1x = tc(k)%x(1:vindI(1)+1)
                                    s1y = tc(k)%y(1:vindI(1)+1)
                                    s2x = tc(k)%x(sizetck:vindI(sizeI):-1)
                                    s2y = tc(k)%y(sizetck:vindI(sizeI):-1)
                                    tc(k)%x = s1x
                                    tc(k)%y = s1y
                                    tc(k)%endsaddle = 0
                                    tempc = tc(k)
                                    tempc%x = s2x 
                                    tempc%y = s2y 
                                    tc = [tc, tempc]
                                    
                                else
                                    call gdErrorHandler('AddTopologicalMeshContours: ' // & 
                                        'unknown error when checking closed ' // & 
                                        'tangency contour (2 intersections), this is likely ' // & 
                                        'a bug') 
                                end if 
                            else
                                ! Multiple intersections, need at least 
                                ! an additional contour
                                tempc = tc(k) 

                                ! Check which parts to delete
                                if (deletestart) then 
                                    s1x = [tc(k)%x(1), tc(k)%x(vindI(1)+1:vindI(2)+1)]
                                    s1y = [tc(k)%y(1), tc(k)%y(vindI(1)+1:vindI(2)+1)]
                                    tc(k)%endsaddle = 0
                                else
                                    s1x = tc(k)%x(1:vindI(1)+1)
                                    s1y = tc(k)%y(1:vindI(1)+1)
                                    tc(k)%endsaddle = 0
                                end if
                                
                                if (deleteend) then 
                                    s2x = [tc(k)%x(sizetck), &
                                        tc(k)%x(vindI(sizeI):vindI(sizeI-1):-1)]
                                    s2y = [tc(k)%y(sizetck), &
                                        tc(k)%y(vindI(sizeI):vindI(sizeI-1):-1)]
                                    tempc%endsaddle = 0
                                    
                                else
                                    s2x = tc(k)%x(sizetck:vindI(sizeI):-1)
                                    s2y = tc(k)%y(sizetck:vindI(sizeI):-1)
                                    tempc%endsaddle = 0
                                end if
                                
                                ! Add
                                tc(k)%x = s1x
                                tc(k)%y = s1y
                                tempc%x = s2x
                                tempc%y = s2y
                                tc = [tc, tempc]
                                
                            end if 
                        end if 
                    else
                        ! Open contour
                        ! Check
                        if (size(vindI) < 1) then 
                            call gdErrorHandler('AddTopologicalMeshContour: '// & 
                                'at least one intersection of contour ' // & 
                                'with other boundary expected')
                        end if
                        
                        ! Delete
                        if (deletefirstpart) then 
                            tc(k)%x = [tc(k)%x(1), &
                                tc(k)%x(vindI(1)+1:vindI(2)+1)]
                            tc(k)%y = [tc(k)%y(1), &
                                tc(k)%y(vindI(1)+1:vindI(2)+1)]
                        else
                            tc(k)%x = tc(k)%x(1:vindI(1)+1)
                            tc(k)%y = tc(k)%y(1:vindI(1)+1)
                        end if 
                    end if 
                end do
           
                ! Add
                allc = [allc, tc]
                call curvetypes%Append(spread(2_I8, 1, size(tc)))
                
                ! Add flux surface ID
                nfs = nfs + 1;
                call fsIDs%Append(spread(nfs, 1, size(tc)))
            end if 
        end do 

        ! Add to the topology mesh
        !=========================
        ! Clean the contours (just to be sure)
        call CleanContours(allc)

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
        call SplitTopologicalMeshFaces(topomesh, magneticField, vessel)    
        
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
        type(VesselUDT), intent(in)             :: vessel 
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
        integer(I8)                             :: i, j

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
            tf = findloc(any(topomesh%face%vert == i, dim=2), .true.)

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
            allocate(tfv(size(tf)))
            do j = 1, size(tf)
                if (topomesh%face%vert(tf(j), 1) == i) then 
                    tfv(j) = topomesh%face%vert(tf(j), 2)
                else
                    tfv(j) = topomesh%face%vert(tf(j), 1) 
                end if 
            end do 
            thisf = tf(minloc(abs(vert%fval(i) - vert%fval(tfv)), dim=1))
            thisfval = minval(abs(vert%fval(i) - vert%fval(tfv)))

            ! Compute the field value to trace
            traceval = options%coreboundariesfrac*(vert%fval(i) - thisfval) + thisfval

            ! Trace contours 
            tc = fieldtracer%TraceContours([thisfval])
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
            contourtypes = [contourtypes, spread(2_I8, 1, size(tc))]

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
        call SplitTopologicalMeshFaces(topomesh, magneticField, vessel)

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
        type(ContourUDT), intent(in)            :: contour
        integer(I8), intent(in)                 :: contourfsID, contourtype

        ! Auxiliary
        real(R8), allocatable, dimension(:)     :: xint, yint, s1r, s2r, &
            tfv
        integer(I8)                             :: nint
        integer(I8), allocatable, dimension(:)  :: s1, s2, fID, vIDs, &
            vtypes, sortind, vf1, vf2, tvIDs, ts2
        logical                                 :: alreadyadded, &
            isinconsistent
        logical, allocatable, dimension(:)      :: iscoinciding, &
            delind
        type(PolygonUDT)                        :: cp 
        type(RealDynamicArrayUDT)               :: xda, yda, s1rda, &
            s2rda
        type(RealDynamicArrayUDT), allocatable  :: xfda(:), yfda(:)
        type(IntegerDynamicArrayUDT)            :: s1da, s2da, fda 

        ! Loop
        integer(I8)                             :: i, k, kold

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
            call fda%Append(spread(i, 1, size(xint)))
            call xda%Append(xint)
            call yda%Append(yint)
            call s1da%Append(s1)
            call s2da%Append(s2)
            call s1rda%Append(s1r)
            call s2rda%Append(s2r)

            ! Housekeeping
            end associate
        end do

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
                    if (((topomesh%face%vert(fIDi, 1)) /= contour%startsaddle) .or. &
                        (topomesh%face%vert(fIDi, 1) == 0) .or. &
                        (contour%startsaddle == 0)) then 
                        ! Inconsistent - throw error later
                        isinconsistent = .true.
                    end if 

                    ! Set ID 
                    vIDs(i) = contour%startsaddle

                elseif ((s1ri == 0_R8) .and. &
                    (s2ri == real(topomesh%face%pol(fIDi)%ne, R8))) then 

                    ! Intersection in start of face and end of contour
                    if (((topomesh%face%vert(fIDi, 1)) /= contour%endsaddle) .or. &
                        (topomesh%face%vert(fIDi, 1) == 0) .or. &
                        (contour%endsaddle == 0)) then 
                        ! Inconsistent - throw error later
                        isinconsistent = .true.
                    end if 

                    ! Set ID 
                    vIDs(i) = contour%endsaddle

                elseif ((s1ri == real(topomesh%face%pol(fIDi)%ne, R8)) .and. &
                    (s2ri == 0_R8)) then 

                    ! Intersection in end of face and start of contour    
                    if (((topomesh%face%vert(fIDi, 1)) /= contour%startsaddle) .or. &
                        (topomesh%face%vert(fIDi, 1) == 0) .or. &
                        (contour%startsaddle == 0)) then 
                        ! Inconsistent - throw error later
                        isinconsistent = .true.
                    end if 

                    ! Set ID 
                    vIDs(i) = contour%startsaddle

                elseif ((s1ri == real(topomesh%face%pol(fIDi)%ne, R8)) .and. &
                    (s2ri == real(topomesh%face%pol(fIDi)%ne, R8))) then 

                    ! Intersection in end of face and end of contour
                    if (((topomesh%face%vert(fIDi, 2)) /= contour%endsaddle) .or. &
                        (topomesh%face%vert(fIDi, 2) == 0) .or. &
                        (contour%endsaddle == 0)) then 
                        ! Inconsistent - throw error later
                        isinconsistent = .true.
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
                allocate(tfv(1))
                call magneticField%interp%Evaluate([xinti], [yinti], 0, 0, tfv)
                call AddTopologicalMeshVertex(topomesh, xinti, yinti, &
                    tfv(1), vtypes(i))
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

        ! Add start and end points as intersections if they have an 
        ! ID (and if that ID is not already present as an intersection)
        tvIDs = vIDs 
        ts2 = s2
        if ((contour%startsaddle /= 0) .and. (contour%startsaddle /= vIDs(1))) then 
            tvIDs = [contour%startsaddle, tvIDs]
            ts2 = [1, s2]
        end if 
        if ((contour%endsaddle /= 0) .and. (contour%endsaddle /= vIDs(size(vIDs)))) then 
            tvIDs = [tvIDs, contour%endsaddle]
            ts2 = [s2, cp%ne]
        end if 

        ! Extract faces
        call ExtractTopologicalFacesFromPolygon(cp, tvIDs, ts2, topomesh%vert%x, &
            topomesh%vert%y, vf1, vf2, xfda, yfda)
        
        ! Add to faces
        do i = 1, size(xfda)
            call AddTopologicalMeshFace(topomesh, [vf1(i), vf2(i)], xfda(i), &
                yfda(i), contourtype, contourfsID)
        end do 
        
        ! Adjust existing faces
        !======================
        ! Sort intersections according to face index
        allocate(sortind(size(s1r)))
        call Sort(s1r, ind=sortind)
        s1r = s1r(sortind)
        fID = fID(sortind)
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

            ! Mark face for deletion
            delind(fID(k)) = .true.

            ! Add start and end points as intersections if they have an 
            ! ID (and if that ID is not already present as an intersection)
            tvIDs = vIDs(kold+1:k) 
            ts2 = s2(kold+1:k)
            if ((topomesh%face%vert(k, 1) /= 0) .and. (topomesh%face%vert(k, 1) /= vIDs(1))) then 
                tvIDs = [topomesh%face%vert(k, 1), tvIDs]
                ts2 = [1, s2]
            end if 
            if ((topomesh%face%vert(k, 2) /= 0) .and. (topomesh%face%vert(k, 2) /= vIDs(size(vIDs)))) then 
                tvIDs = [tvIDs, topomesh%face%vert(k, 2)]
                ts2 = [s2, topomesh%face%pol(k)%ne]
            end if 

            ! Extract faces
            call ExtractTopologicalFacesFromPolygon(&
                topomesh%face%pol(fID(k)), tvIDs, ts2, topomesh%vert%x, &
                topomesh%vert%y, vf1, vf2, xfda, yfda)

            ! Add to faces
            do i = 1, size(xfda)
                call AddTopologicalMeshFace(topomesh, [vf1(i), vf2(i)], xfda(i), &
                    yfda(i), topomesh%face%type(k), topomesh%face%fsID(k))
            end do 

        end do 

        ! Remove adjusted faces
        allocate(delind(topomesh%face%ntot))
        delind = .false. 
        delind(fID) = .true. 
        call RemoveTopologicalMeshFaceLogical(topomesh, delind)

    end subroutine 

    ! Topological mesh face splitting
    subroutine SplitTopologicalMeshFaces(topomesh, magneticField, vessel)

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
        type(magneticFieldUDT)                  :: magneticField
        type(VesselUDT)                         :: vessel 

        ! Auxiliary 
        integer(I8)                             :: nncf, nndf, np, &
            ind(1:2), nfinit
        real(R8)                                :: newvcfx(1:2), newvcfy(1:2), &
            newvcff(1:2), newvdfx(1), newvdfy(1), newvdff(1)
        
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
        allocate(isduplicateface(ntot))
        isduplicateface = .false. 
        do i = 1, ntot-1
            do j = i+1, ntot
                ! Skip closed faces - should be dealt with separately, even if
                ! multiple are present
                if ((.not. isclosedface(i)) .and. (.not. isclosedface(j))) then 
                    ! Check vertices
                    if (any(face%vert(i, 1) == face%vert(j, :)) .and. &
                        any(face%vert(i, 2) == face%vert(j, :))) then 
                        isduplicateface(j) = .true.
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
            print *, 'SplitTopologicalMeshFaces: ', nncf, &
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
                call magneticField%interp%Evaluate(newvcfx, newvcfy, &
                    0, 0, newvcff)
                
                ! Insert new vertices
                do j = 1, 2
                    call AddTopologicalMeshVertex(topomesh, &
                        newvcfx(j), newvcfy(j), newvcff(j), &
                        TMvertexsplitID)
                end do 
                
                ! Insert first face
                xrda = ConstructRealDynamicArray(face%x(i)%Get([(k, k = 1, ind(1))]))
                yrda = ConstructRealDynamicArray(face%y(i)%Get([(k, k = 1, ind(1))]))
                call AddTopologicalMeshFace(topomesh, &
                    [face%vert(i, 1), topomesh%vert%ntot-1], &
                    xrda, yrda, face%type(i), face%fsID(i))

                ! Insert second face
                xrda = ConstructRealDynamicArray(face%x(i)%Get([(k, k = ind(1), ind(2))]))
                yrda = ConstructRealDynamicArray(face%y(i)%Get([(k, k = ind(1), ind(2))]))
                call AddTopologicalMeshFace(topomesh, &
                    [topomesh%vert%ntot-1, topomesh%vert%ntot], &
                    xrda, yrda, face%type(i), face%fsID(i))

                ! Insert third face
                xrda = ConstructRealDynamicArray(face%x(i)%Get([(k, k = ind(2), face%x(i)%Size())]))
                yrda = ConstructRealDynamicArray(face%y(i)%Get([(k, k = ind(2), face%y(i)%Size())]))
                call AddTopologicalMeshFace(topomesh, &
                    [topomesh%vert%ntot, face%vert(i, 2)], &
                    xrda, yrda, face%type(i), face%fsID(i))

            elseif (isduplicateface(i)) then 

                ! Get number of points of this face
                np = face%x(i)%Size()
                if (np < 3) then  ! end points should be the same and duplicate
                    ! Issue message: we cannot split up this boundary
                    call gdErrorHandler('SplitTopologicalMeshFaces: ' // & 
                        'face with same vertices found with only ' // & 
                        'two coordinates, cannot split up')
                end if 
                
                ! Split up the face into parts with approx. equal number of
                ! vertices. 
                ind(1) = np/2
                
                ! Get vertex coordinates
                newvdfx = face%x(i)%Get(ind)
                newvdfy = face%y(i)%Get(ind)
                call magneticField%interp%Evaluate(newvdfx, newvdfy, &
                    0, 0, newvdff)
                
                ! Insert new vertex
                call AddTopologicalMeshVertex(topomesh, newvdfx(1), &
                    newvdfy(1), newvdff(1), TMvertexsplitID);
                
                ! Insert first face
                xrda = ConstructRealDynamicArray(face%x(i)%Get([(k, k = 1, ind(1))]))
                yrda = ConstructRealDynamicArray(face%y(i)%Get([(k, k = 1, ind(1))]))
                call AddTopologicalMeshFace(topomesh, &
                    [face%vert(i, 1), topomesh%vert%ntot], &
                    xrda, yrda, face%type(i), face%fsID(i))

                ! Insert second face
                xrda = ConstructRealDynamicArray(face%x(i)%Get([(k, k = ind(1), face%x(i)%Size())]))
                yrda = ConstructRealDynamicArray(face%y(i)%Get([(k, k = ind(1), face%y(i)%Size())]))
                call AddTopologicalMeshFace(topomesh, &
                    [topomesh%vert%ntot, face%vert(i, 2)], &
                    xrda, yrda, face%type(i), face%fsID(i))

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
        type(PolygonUDT), allocatable           :: fxp(:), fyp(:)

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

        ! Trace 
        call fieldtracer%SetValues(fx)
        fxc = fieldtracer%TraceContours([0.0_R8])
        call fieldtracer%SetValues(fy)
        fyc = fieldtracer%TraceContours([0.0_R8])

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

        ! Compute intersections
        do i = 1, nfxc
            ! Only compute intersections with other polygons
            do j = 1, nfyc
                ! Compute intersections with next polygon
                call PolygonIntersections(fxp(i), fyp(j), tx, ty, ts1, ts2)
                
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
                end do 

                ! Append
                call xc%Append(tx)
                call yc%Append(ty)
                call fc%Append(tf)
                call tc%Append(tt)

                ! Housekeeping
                deallocate(tt, tf, tfxx, tfyy, tfxy)
            end do 
        end do
        
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
                    tv = pack(p%vert, extrloc)
                    tf = val(tv)
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

        ! Substructures
        call topomesh%vert%Initialize()
        call topomesh%face%Initialize()
        call topomesh%cell%Initialize()
        

    end subroutine

    subroutine InitializeTopologicalMeshVertex(tpvert)

        ! Description
        !============
        ! Initialize the topomesh vertex structure (simply empty arrays)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshVertUDT)      :: tpvert 

        ! Initialize
        !===========
        tpvert%ntot = 0
        allocate(tpvert%ID(0), tpvert%x(0), tpvert%y(0), tpvert%fval(0), &
            tpvert%type(0))

    end subroutine 

    subroutine InitializeTopologicalMeshFace(tpface)

        ! Description
        !============
        ! Initialize the topomesh vertex structure (simply empty arrays)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshFaceUDT)      :: tpface

        ! Initialize
        !===========
        tpface%ntot = 0
        allocate(tpface%ID(0), tpface%vert(0, 2), tpface%cell(0), &
            tpface%fsID(0), tpface%x(0), tpface%y(0), tpface%fval(0), &
            tpface%pol(0), tpface%type(0))

    end subroutine

    subroutine InitializeTopologicalMeshCell(tpcell)

        ! Description
        !============
        ! Initialize the topomesh cell structure (simply empty arrays)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshCellUDT)      :: tpcell

        ! Initialize
        !===========
        tpcell%ntot = 0
        allocate(tpcell%ID(0))

    end subroutine

    ! Vertex addition
    subroutine AddTopologicalMeshVertex(topomesh, x, y, F, t)

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
        integer(I8), intent(in)             :: t 

        ! Concatenate
        !============
        topomesh%vert%ntot = topomesh%vert%ntot + 1
        topomesh%vert%x = [topomesh%vert%x, x]
        topomesh%vert%y = [topomesh%vert%y, y]
        topomesh%vert%fval = [topomesh%vert%fval, F]
        topomesh%vert%type = [topomesh%vert%type, t]
        topomesh%vert%ID = [topomesh%vert%ID, topomesh%vert%ntot]

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
            tfsorted 
        real(R8)                                :: vx, vy
        real(R8), allocatable, dimension(:)     :: npx, npy, dx, dy,&
            theta

        ! Loop
        integer(I8)                             :: i, j

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
        do i = 1, topomesh%vert%ntot
            ! Get face indices
            tf = findloc(any(topomesh%face%vert == i, dim=2), .true.)
            
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
            
        end do


    end subroutine

    ! Face addition
    subroutine AddTopologicalMeshFace(topomesh, facevert, x, y, &
        t, fsID)

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

        ! Auxiliary
        integer(I8), allocatable            :: tmp(:, :)
        real(R8), allocatable               :: xt(:), yt(:), dxt(:), &
            dyt(:)
        logical, allocatable                :: delind(:)
        type(PolygonUDT)                    :: tp 

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
        integer(I8), allocatable, dimension(:)  :: ncpf, tf, fc, ind

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
                    call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                        'could not find next face')
                end if 
                ! All faces added, exit
                exit 
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
            tv = face%vert(tf, tvind);
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
            cc = cc + 1;
            thiscellvert = ConstructIntegerDynamicArray(tcv)
            thiscellface = ConstructIntegerDynamicArray(tcf)
            cellvert = [cellvert, thiscellvert]
            cellface = [cellface, thiscellface]
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

        topomesh%cell%ntot = cc;
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
            topomesh%cell%vert(topomesh%cell%vertP(i, 1):topomesh%cell%vertP(i, 2)) = & 
                cellvert(i)%Get()
            topomesh%cell%face(topomesh%cell%faceP(i, 1):topomesh%cell%faceP(i, 2)) = & 
                cellface(i)%Get()
        end do 

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
        !   faces are those faces with type 3, all the rest are internal. 

        ! Note: it is assumed that all necessary contours/faces have 
        ! been added already
        
        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                  :: topomesh 

        ! Loop
        integer(I8)                         :: i 

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

        ! Labels
        if (allocated(topomesh%face%label)) then 
            deallocate(topomesh%face%label)
        end if 
        allocate(topomesh%face%label(topomesh%face%ntot))
        do i = 1, topomesh%face%ntot
            topomesh%face%label(i) = i
        end do 

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
        logical, allocatable, dimension(:)      :: wasfound(:)

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
        allocate(tempvcells(v%ntot, maxval(ncpv)), vcount(v%ntot), &
            wasfound(v%ntot))
        tempvcells = 0
        vcount = 0
        do i = 1, nc
            ! Reset logical
            wasfound = .false. ! to hedge for doubly appearing vertices

            ! Get vertices of cell
            tv = GetTMCellVert(c, i) ! there may be doubles in here!
            
            ! Add vertex cells
            do j = 1, size(tv)
                if (.not. wasfound(tv(j))) then 
                    ! Update counter
                    vcount(tv(j)) = vcount(tv(j))+1

                    ! Update logical
                    wasfound(tv(j)) = .true.

                    ! Add cell
                    tempvcells(tv(j), vcount(tv(j))) = i 
                end if 
            end do 
        end do 
        
        ! Construct 
        do i = 1, v%ntot
            v%cell(v%cellP(i, 1):v%cellP(i, 1)+v%cellP(i, 2)) = & 
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
        integer(I8)                             :: i

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
            where (rmvert(tv)) 
                tv = 0
            elsewhere
                tv = tv - diffID(tv)
            end where
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
            ! Check number of faces
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

            ! Hedge for machine precision effects - shouldn't be necessary anymore?
            !if (sortedsID(size(sID)) > pol%ne) then 
            !    ! Check if it's only off by one -  otherwise it's an error
            !    if (sortedsID(size(sID)) == pol%ne+1) then 
            !        ! Shrink
            !        sortedsID(size(sID)) = pol%ne 
            !    else
            !        call gdErrorHandler('ExtractTopologicalFacesFromPolygon: ' // & 
            !            'segment index exceeds number of polygon edges, check input')
            !    end if 
            !end if 

            ! Treat start and end segments
            if (pol%isclosed) then 
                ! Closed polygon
                if ((sortedsID(1) /= 1) .and. (sortedsID(nsID) /= pol%ne)) then 
                    ! Start and end segment without start and end vertex,
                    ! so merge
                    fc = fc + 1

                    ! Add vertices
                    v1(fc) = sortedeID(neID)
                    v2(fc) = sortedeID(1)

                    ! Add face coordinates
                    call xf(fc)%Append([xint(sortedeID(neID)), &
                        px(sortedsID(nsID)+1:pol%ne), &
                        px(1:sortedsID(1)), xint(sortedeID(1))])
                    call yf(fc)%Append([yint(sortedeID(neID)), &
                        py(sortedsID(nsID)+1:pol%ne), &  
                        py(1:sortedsID(1)), yint(sortedeID(1))])
                elseif (((sortedsID(1) == 1) .and. (sortedsID(nsID) /= pol%ne)) .or. &
                    ((sortedsID(1) /= 1) .and. (sortedsID(nsID) == pol%ne))) then 
                        ! Something weird here: intersection in exactly first node, but
                        ! only at first and not at last segment. Likely because
                        ! intersection with last segment was not added. Throw error. 
                        call gdErrorHandler(&
                            'ExtractTopologicalFacesFromPolygon: ' // & 
                            'closed polygon intersects in start point ' // & 
                            'but not in end point - check input '// & 
                            '(end point intersection should be added separately)')

                else
                    ! Nothing to do: intersection exactly at start and end point. 
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
                    call xf(fc)%Append([px(1:sortedsID(1)), xint(sortedeID(1))]) 
                    call yf(fc)%Append([py(1:sortedsID(1)), yint(sortedeID(1))])
                
                end if 
                if ((sortedsID(nsID) == 1) .and. (px(1) == xint(sortedeID(nsID))) & 
                    .and. (py(1) == yint(sortedeID(nsID)))) then 
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
                    call xf(fc)%Append([xint(sortedeID(neID)), px(sortedsID(nsID)+1:pol%ne+1)]) 
                    call yf(fc)%Append([yint(sortedeID(neID)), py(sortedsID(nsID)+1:pol%ne+1)])
                
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

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        type(magneticFieldUDT), intent(in)      :: magneticField 
        type(VesselUDT), intent(in)             :: vessel 

        ! Auxiliary
        real(R8), allocatable, dimension(:)     :: Vv
        logical, allocatable, dimension(:)      :: outbnd, rmvert, &
            rmface

        ! Loop
        integer(I8)                             :: i 

        ! Initialize
        !===========
        ! Unpack
        associate(&
            plf         => vessel%plfvessel,      &
            mfinterp    => magneticField%interp)

        ! Vertices
        !=========
        ! Check if in boundary
        allocate(Vv(topomesh%vert%ntot))
        call plf%Evaluate(topomesh%vert%x, topomesh%vert%y, 0, 0, Vv)
        outbnd = Vv >= 0

        ! Check if we should remove it
        rmvert = outbnd .and. (topomesh%vert%type /= 4) .and. &
            (topomesh%vert%type /= 5) .and. (topomesh%vert%type /= 6)

        ! Remove these vertices
        call RemoveTopologicalMeshVertexLogical(topomesh, rmvert)

        ! Faces
        !======
        ! Start by removing faces with zero start or end vertex
        rmface = (topomesh%face%vert(:, 1) == 0) .or. &
            (topomesh%face%vert(:, 2) == 0)

        ! Remove
        call RemoveTopologicalMeshFaceLogical(topomesh, rmface)

        ! Check remaining faces 
        deallocate(rmface)
        allocate(rmface(topomesh%face%ntot))
        rmface = .false. 
        do i = 1, topomesh%face%ntot
            if (topomesh%face%type(i) /= 3) then 
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

    ! Getters
    function GetTMVertFace(vert, i) result(res)
        integer(I8)                 :: i 
        type(TopomeshVertUDT)       :: vert 
        integer(I8), allocatable    :: res(:)
        res = vert%face(vert%faceP(i, 1):(vert%faceP(i, 1) +  vert%faceP(i, 2) - 1))
    end function

    function GetTMCellVert(cell, i) result(res)
        integer(I8)                 :: i 
        type(TopomeshCellUDT)       :: cell 
        integer(I8), allocatable    :: res(:)
        res = cell%vert(cell%vertP(i, 1):(cell%vertP(i, 1) + cell%vertP(i, 2) - 1))
    end function

    function GetTMCellFace(cell, i) result(res)
        integer(I8)                 :: i 
        type(TopomeshCellUDT)       :: cell 
        integer(I8), allocatable    :: res(:)
        res = cell%face(cell%faceP(i, 1):(cell%faceP(i, 1) + cell%faceP(i, 2) - 1))
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
        type(TopomeshVertUDT)       :: vert 
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

    ! Contour clean-up
    subroutine CleanContours(contours)

        ! Description
        !============
        ! This routine cleans up the contours, i.e. it removes subsequent
        ! points that are up to disttol coinciding. May be necessary
        ! for later intersection computing etc. 

        ! Declare variables
        !==================
        ! Arguments
        class(ContourUDT), intent(inout)        :: contours(:)

        ! Auxiliary
        logical, allocatable                    :: delind(:)
        real(R8), allocatable                   :: dx(:), dy(:)

        ! Loop
        integer(I8)                             :: i 
        
        ! Clean
        !======
        do  i = 1, size(contours)
            dx = contours(i)%x(2:size(contours(i)%x)) - &
                contours(i)%x(1:size(contours(i)%x)-1)
            dy = contours(i)%y(2:size(contours(i)%y)) - &
                contours(i)%y(1:size(contours(i)%y)-1)
            allocate(delind(size(dx)))
            delind = (abs(dx) <= disttol) .and. (abs(dy) <= disttol)
            if (any(delind)) then 
                contours(i)%x = pack(contours(i)%x, [.true., .not. delind])
                contours(i)%y = pack(contours(i)%y, [.true., .not. delind])
            end if 
            deallocate(delind)
        end do 
    end subroutine 

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
        ! 'ID, x, y, fval, BV'
        ! <ID, x, y, fval, BV (as zero or one)>
        ! 'faces'
        ! <face%ntot> 
        ! 'ID, fsID, type, vert1, vert2, BF'
        ! <ID, fsID, type, vert(:, 1), vert(:, 2), BF>
        ! 'face <nf>' <repeated for each face including header>
        ! <x, y> 
        ! 'cells'
        ! <cell%ntot,  cells%nvert, cells%nface> 
        ! 'cell vertices'
        ! <cell%vert>
        ! 'cell vertex pointer'
        ! <cell%vertP>
        ! 'cell faces'
        ! <cell%face> 
        ! 'cell face pointer'
        ! <cell%faceP(:, 1), cell%faceP(:, 2)>
        
        ! Declare variables
        !==================
        ! Modules 
        use mod_plotter 
        use mod_specialchars, only : filesepchar

        ! Arguments
        class(TopomeshUDT)                      :: topomesh
        character(*), intent(in)                :: filename 

        ! Auxiliary
        integer                                 :: fu, BVval
        real(R8), allocatable, dimension(:)     :: xf, yf
        character(:), allocatable               :: dir

        ! Loop
        integer(I8)                             :: i, j 

        ! Initialize
        !===========
        ! Unpack
        associate(&
            v       => topomesh%vert,   &
            f       => topomesh%face,   &
            c       => topomesh%cell    &
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
        if (.not. allocated(v%BV)) then 
            allocate(v%BV(v%ntot))
            v%BV = .false.
        end if 
        if (.not. allocated(f%BF)) then 
            allocate(f%BF(f%ntot))
            f%BF = .false.
        end if

        ! Write vertex data
        !==================
        ! Number of vertices
        write (fu, *) 'vertices'
        write (fu, *) v%ntot

        ! Basic vertex data
        write (fu, *) 'ID, x, y, fval, BV'
        do i = 1, v%ntot
            if (v%BV(i)) then 
                BVval = 1
            else 
                BVval = 0
            end if
            write (fu, *) v%ID(i), v%x(i), v%y(i), v%fval(i), BVval
        end do 

        ! Write face data
        !================
        ! Number of faces
        write (fu, *) 'faces'
        write (fu, *) f%ntot 

        ! Basic face data
        write (fu, *) 'ID, fsID, type, vert1, vert2, BF'
        do i = 1, f%ntot 
            if (f%BF(i)) then 
                BVval = 1
            else 
                BVval = 0
            end if 
            write (fu, *) f%ID(i), f%fsID(i), f%type(i), f%vert(i, 1), &
                f%vert(i, 2), BVval 
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
        write (fu, *) 'cells'
        write (fu, *) c%ntot, size(c%vert), size(c%face)

        ! Cell vertices
        write (fu, *) 'cell vertices'
        do i = 1, size(c%vert)
            write (fu, *) c%vert(i)
        end do 
        write (fu, *) 'cell vertex pointer'
        do i = 1, c%ntot 
            write (fu, *) c%vertP(i, 1), c%vertP(i, 2)
        end do 

        ! Cell faces
        write (fu, *) 'cell faces'
        do i = 1, size(c%face)
            write (fu, *) c%face(i)
        end do 
        write (fu, *) 'cell face pointer'
        do i = 1, c%ntot 
            write (fu, *) c%faceP(i, 1), c%faceP(i, 2)
        end do 

        ! Housekeeping
        !=============
        end associate
        close(fu)

    end subroutine


end module