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
    use mod_linearsolverinterface, only: SolveDenseLinearSystemDI
    use goatmod_types, only : magneticFieldUDT, VesselUDT
    use goatmod_userinput, only : TopomeshOptionsUDT
    implicit none
    private 
    public :: TopomeshUDT, ConstructTopologicalMesh, TraceExtrema2D, &
        TraceTangencyPoints2D

    ! Module parameters
    real(R8), parameter, private        :: tprelfieldtol = 1e-2 ! relative field tolerance under which extrema are removed

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
            flags(:), type(:)
        real(R8), allocatable           :: x(:), y(:), fval(:)

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
            flags(:), fsID(:), type(:)
        type(RealDynamicArrayUDT), allocatable  :: x(:), y(:)
        real(R8), allocatable           :: fval(:)
        type(PolygonUDT), allocatable   :: pol(:)

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

    end type

    ! General topological mesh type
    type :: TopomeshUDT 

        type(TopomeshVertUDT)   :: vert 
        type(TopomeshFaceUDT)   :: face 
        type(TopomeshCellUDT)   :: cell 
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
        real(R8), allocatable, dimension(:)     :: xtp, ytp, Ftp
        integer(I8)                             :: nv 
        integer(I8), allocatable, dimension(:)  :: typetp

        ! Initialize
        !===========
        ! Initialize topomesh 
        call topomesh%Initialize()

        ! Determine domain bounds based on vessel extent
        call vessel%plfvessel%ps%GetVertices(xps, yps)
        xb = [minval(xps), maxval(xps)]
        yb = [minval(yps), maxval(yps)]

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

        ! Compute tangency points
        !========================
        call TraceTangencyPoints2D(xtp, ytp, typetp, Ftp, &
            vesseltracer, magneticField)

        ! Update the tracers
        !===================
        ! Construct refined grid based on tangency points
        call ConstructRefined2DStructuredGrid(xg, yg, xgv, ygv, xb, yb, &
            options%vresx, options%vresy, xtp, ytp, 5, 5)  

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
            vf1, vf2, teid, tsid
        real(R8)                                :: tempx, tempy
        real(R8), allocatable, dimension(:)     :: xg, yg, f, fx, fy, &
            tx, ty, thiseig, tf, tfxx, tfxy, tfyy, xe, ye, fe 
        logical                                 :: conv, donewton 
        type(RealDynamicArrayUDT)               :: xc, yc, fc
        type(RealDynamicArrayUDT), allocatable  :: xfrda(:), yfrda(:) 
        type(IntegerDynamicArrayUDT), allocatable   :: fxpeid(:), &
            fxpsid(:), fypeid(:), fypsid(:)
        type(IntegerDynamicArrayUDT)            :: tc
        type(ContourUDT), allocatable           :: fxc(:), fyc(:)
        type(PolygonUDT), allocatable           :: fxp(:), fyp(:)

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

        ! Initialize intersection trackers
        allocate(fxpeid(nfxc), fxpsid(nfxc), fypeid(nfyc), fypsid(nfyc))
        do i = 1, nfxc 
            fxpeid(i) = ConstructIntegerDynamicArray()
            fxpsid(i) = ConstructIntegerDynamicArray()
        end do 
        do i = 1, nfyc 
            fypeid(i) = ConstructIntegerDynamicArray()
            fypsid(i) = ConstructIntegerDynamicArray() 
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

                    ! Update counter
                    ec = ec + 1

                    ! Store intersection data
                    call fxpeid(i)%Append(ec)
                    call fxpsid(i)%Append(ts1(k))
                    call fypeid(j)%Append(ec)
                    call fypsid(j)%Append(ts2(k))
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
        typee = typee + 2 ! 1: minimum, 2: saddle, 3: maximum

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
            call ExtractTopologicalFacesFromPolygon(fxp(i), teid, &
                tsid, xe, ye, vf1, vf2, xfrda, yfrda)

            ! Add to faces
            do j = 1, size(vf1)
                call AddTopologicalMeshFace(topomesh, [vf1(j), vf2(j)], &
                    xfrda(j), yfrda(j), 1, 0)
            end do 
        end do 

        ! Loop over fyp
        do i = 1, nfyc 
            ! Extract the faces from this polygon
            teid = fypeid(i)%Get()
            tsid = fypsid(i)%Get()
            call ExtractTopologicalFacesFromPolygon(fyp(i), teid, &
                tsid, xe, ye, vf1, vf2, xfrda, yfrda)

            ! Add to faces
            do j = 1, size(vf1)
                call AddTopologicalMeshFace(topomesh, [vf1(j), vf2(j)], &
                    xfrda(j), yfrda(j), 1, 0)
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
                ttype = 4
            else
                ttype = 5
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
                ! shouldn't be possible at this stage
                if (tc(i)%startsaddle == tc(i)%endsaddle) then 
                    print *, 'AddTopologicalMeshTangencyPoints: segment detected that closes upon itself. Removing...'
                    cycle
                end if 

                ! Add the face
                ftype = 3 ! vessel boundary -> outer boundary
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
        ! Substructures
        call topomesh%vert%Initialize()
        call topomesh%face%Initialize()

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
        allocate(tpvert%ID(0), tpvert%x(0), tpvert%y(0), tpvert%fval(0))

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
            tpface%pol(0))

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

    ! Face addition
    subroutine AddTopologicalMeshFace(topomesh, facevert, x, y, &
        t, fsID)

        ! Description
        !============
        ! Add a face to the topological mesh. Not optimized for 
        ! memory usage (e.g. less allocate/deallocate), but 
        ! shouldn't be an issue since topomeshes are usually small. Note
        ! that no interconnection data is updated!

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                  :: topomesh 
        type(RealDynamicArrayUDT)           :: x, y 
        integer(I8), intent(in)             :: t, facevert(1:2), fsID

        ! Auxiliary
        integer(I8), allocatable            :: tmp(:, :)

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
        topomesh%face%ID = [topomesh%vert%ID, topomesh%vert%ntot]

    end subroutine

    ! Cell addition 

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

            ! The first output argument has the fields:
            ! - vert: nf-by-2 vector containing the IDs of the intersections. The end
            ! parts, if no intersections were exactly at the end of the polygon, will
            ! have at least one entry to be zero. The first is the start vertex, the
            ! second the end vertex. 

            ! The second output arguments has the fields:
            ! - x, y: coordinates of the face with added end points if applicable. The
            ! first entry corresponds to the start vertex, the second to the end vertex
            ! (should those be initialized). 

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
            integer(I8), allocatable, dimension(:)  :: sortind(:), &
                sortedsID(:), sortedeID(:)
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
                nsID    =>      size(sID),          &
                neID    =>      size(eID),          &
                px      =>      pol%x(pol%vert),    &
                py      =>      pol%y(pol%vert)     &
                )

            ! Construct segments
            !===================
            ! Sort 
            sortedsID = sID 
            call Sort(sortedsID, ind=sortind)
            sortedeID = eID(sortind)

            ! Hedge for machine precision effects 
            if (sortedsID(size(sID)) > pol%ne) then 
                ! Check if it's only off by one -  otherwise it's an error
                if (sortedsID(size(sID)) == pol%ne+1) then 
                    ! Shrink
                    sortedsID(size(sID)) = pol%ne 
                else
                    call gdErrorHandler('ExtractTopologicalFacesFromPolygon: ' // & 
                        'segment index exceeds number of polygon edges, check input')
                end if 
            end if 

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
                        px(sortedeID(neID)+1:pol%ne), &
                        px(1:sortedsID(1)), xint(sortedeID(1))])
                    call yf(fc)%Append([yint(sortedeID(neID)), &
                        py(sortedeID(neID)+1:pol%ne), &  
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
                if ((sortedsID(1) == 1) .and. (px(1) == xint(sortedsID(1))) & 
                    .and. (py(1) == yint(sortedsID(1)))) then 
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
                if ((sortedsID(nsID) == 1) .and. (px(1) == xint(sortedsID(nsID))) & 
                    .and. (py(1) == yint(sortedsID(nsID)))) then 
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
                    call xf(fc)%Append([xint(sortedeID(neID)), px(sortedsID(nsID)+1:pol%ne)]) 
                    call yf(fc)%Append([yint(sortedeID(neID)), py(sortedsID(nsID)+1:pol%ne)])
                
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
            end do

            ! Add
            call xf(fc)%Append(tx)
            call yf(fc)%Append(ty)


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

    ! eigenvalue computation for 2-by-2 matrices 
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


end module