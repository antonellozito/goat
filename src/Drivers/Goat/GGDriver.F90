subroutine GGDriver(goatoptions)

    ! Description
    !============
    ! This driver runs the grid deformation in standalone mode. The 
    ! goatoptions should be passed to this routine to identify which
    ! files to load etc. 

    ! Initialize
    !===========
    ! Modules
    use goatmod_types 
    use goatmod_userinput
    use ggmod_topology2D
    use ggmod_gridgeneration2D
    use mod_contour2D, only : ContourTracerUDT, ConstructStructuredTracer, &
        ContourUDT
    use mod_streamlinetracing2D, only: StreamlineTracerUDT, &
        ConstructStructuredStreamlineTracer
    use mod_structured2Dgridding

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    type(GoatoptionsUDT)        :: goatoptions
    type(GridUDT), allocatable  :: grid
    type(MagneticFieldUDT)      :: magneticField
    type(EnvironmentUDT)        :: environment

    ! Other options
    type(GGoptionsUDT)          :: ggoptions
    type(TopomeshOptionsUDT)    :: topomeshoptions

    ! Auxiliary
    type(TopomeshUDT), allocatable  :: topomesh
    class(ContourTracerUDT), allocatable    :: fieldtracer, vesseltracer
    class(StreamlineTracerUDT), allocatable :: streamlinetracer
    type(PolygonSetUDT)         :: voidps
    type(GGTMDataUDT), allocatable  :: ggtmdata

    real(R8), allocatable, dimension(:)     :: xb, yb, xps, &
        yps, xg, yg, Vf, Vv, xgv, ygv, Vfx, Vfy, defx, defy, newx, newy, &
        demx, demy, keepx, keepy
    logical, allocatable, dimension(:)      :: isnew, keepmask
    real(R8), parameter                     :: emptyR8(0)= 0
    integer(I8)                             :: nv, iter, ie
    integer(I8), parameter                  :: emptyI8(0) = 0
    integer(I8), parameter                  :: pinmaxiter = 20
    logical, parameter                      :: emptyL(0) = .false.

    ! Loop

    ! Initialize
    !===========
    ! Read and extract data
    call ExtractGGData(magneticField, environment, goatoptions)

    ! Set topological mesh options
    topomeshoptions%inputfilepath = goatoptions%inputfilepath 
    call topomeshoptions%Set()

    ! Set grid generation options
    ggoptions%inputfilepath = goatoptions%inputfilepath 
    call ggoptions%Set()

    ! Construct (initial) tracers
    !============================
    ! Determine domain bounds based on vessel and magnetic field extent
    call environment%vessel%plfvessel%ps%GetVertices(xps, yps)
    xb = [minval([xps, magneticField%interp%xgv(1:size(magneticField%interp%xgv))]), &
        maxval([xps, magneticField%interp%xgv(1:size(magneticField%interp%xgv))])]
    yb = [minval([yps, magneticField%interp%ygv(1:size(magneticField%interp%ygv))]), &
        maxval([yps, magneticField%interp%ygv(1:size(magneticField%interp%ygv))])]

    ! Construct a 2D structured grid for tracing (may be extended
    ! in the future for different grid types)
    nv = topomeshoptions%vresx*topomeshoptions%vresy
    allocate(xg(nv), yg(nv), Vf(nv), Vv(nv), xgv(topomeshoptions%vresx), &
        ygv(topomeshoptions%vresy))
    call Construct2DStructuredUniformGrid(xg, yg, xgv, ygv, xb, yb, &
        topomeshoptions%vresx,  topomeshoptions%vresy, 0.0_R8, 0.0_R8)

    ! Evaluate the field and vessel values
    allocate(Vfx(nv), Vfy(nv))
    call environment%vessel%plfvessel%Evaluate(xg, yg, 0, 0, Vv)
    call magneticField%interp%Evaluate(xg, yg, 0, 0, Vf)
    call magneticField%interp%Evaluate(xg, yg, 1, 0, Vfx)
    call magneticField%interp%Evaluate(xg, yg, 0, 1, Vfy)

    ! Build the grid, iterating to pin structure endpoints
    !======================================================
    ! Each pass (re)builds the tracers, the topological mesh and the grid. When
    ! endpoint pinning is on, structure endpoints still farther than the
    ! tolerance from a grid boundary node (measured along the wall contour) are
    ! appended to the forced-flux-surface list and the grid is rebuilt, until
    ! every endpoint is within tolerance or the iteration cap is reached. With
    ! pinning off the loop runs exactly once. The forced-surface list starts
    ! empty (set in the topomesh option defaults).
    do iter = 0, pinmaxiter

        ! Fresh derived types for this build
        if (allocated(topomesh)) deallocate(topomesh)
        if (allocated(grid))     deallocate(grid)
        if (allocated(ggtmdata)) deallocate(ggtmdata)
        allocate(topomesh, grid, ggtmdata)

        ! (Re)construct the tracers
        fieldtracer = ConstructStructuredTracer(&
            reshape(Vf, [topomeshoptions%vresx, topomeshoptions%vresy]), xgv, ygv, &
            emptyR8, emptyR8, emptyR8, emptyi8, emptyL, topomeshoptions%npmin, &
            topomeshoptions%npmax, topomeshoptions%dl)
        vesseltracer = ConstructStructuredTracer(&
            reshape(Vv, [topomeshoptions%vresx, topomeshoptions%vresy]), xgv, ygv, &
            emptyR8, emptyR8, emptyR8, emptyI8, emptyL, topomeshoptions%npmin, &
            topomeshoptions%npmax, topomeshoptions%dl)
        streamlinetracer = ConstructStructuredStreamlineTracer(&
            reshape(Vfx, [topomeshoptions%vresx, topomeshoptions%vresy]), &
            reshape(Vfy, [topomeshoptions%vresx, topomeshoptions%vresy]), &
            xgv, ygv, step=ggoptions%orthtracerstep, nsteps=ggoptions%orthtracernsteps)

        ! Topological mesh
        call ConstructTopologicalMesh(environment%vessel, magneticField, &
            topomeshoptions, topomesh, fieldtracer, vesseltracer, streamlinetracer)

        ! Grid
        call GenerateUnstructuredAlignedGrid(grid, topomesh, magneticField, &
            environment%vessel, fieldtracer, vesseltracer, streamlinetracer, &
            ggoptions, ggtmdataopt=ggtmdata)

        ! With pinning off, a single pass
        if (.not. ggoptions%pinstructureendpoints) exit

        ! Which structure endpoints are still beyond tolerance (measured along
        ! the wall contour)? This reports ALL deficient endpoints.
        call FindDeficientStructureEndpoints(grid, environment%vessel, &
            ggoptions%pinstructureendpointstol, defx, defy)

        ! RADIAL-ONLY pinning (the poloidal-fan path is disabled: it corrupted
        ! the face labels and its inserted nodes did not survive). New deficient
        ! endpoints are added to the radial list; endpoints already forced but
        ! still deficient are the fan/grazing corners the radial path cannot
        ! reach - reported as unpinnable, not retried, so the loop terminates.
        allocate(newx(0), newy(0))
        do ie = 1, size(defx)
            if (size(ggoptions%forcedx) > 0) then
                if (any((ggoptions%forcedx-defx(ie))**2 + &
                        (ggoptions%forcedy-defy(ie))**2 < 1e-12_R8)) cycle
            end if
            newx = [newx, defx(ie)]; newy = [newy, defy(ie)]
        end do

        if (size(newx) == 0) then
            print *, 'pinstructureendpoints: converged after', iter, &
                'iteration(s);', size(ggoptions%forcedx), 'pinned;', size(defx), &
                'endpoint(s) remain beyond tolerance (unpinnable fan/grazing)'
            exit
        end if

        if (iter == pinmaxiter) then
            print *, 'pinstructureendpoints: WARNING iteration cap reached;', &
                size(defx), 'endpoint(s) still beyond tolerance'
            exit
        end if

        ggoptions%forcedx = [ggoptions%forcedx, newx]
        ggoptions%forcedy = [ggoptions%forcedy, newy]
        print *, 'pinstructureendpoints: iteration', iter, '- forcing', &
            size(newx), 'more surface(s) (total', size(ggoptions%forcedx), ')'
        deallocate(newx, newy)

    end do

    ! Write data
    !===========
    ! Translate labels etc
    call TranslateGridLabels(grid, topomesh, environment%vessel, ggtmdata, &
        ggoptions, 'solps')

    ! Recompute topological data from grid for new face labels
    call ComputeTopologicalData(grid, topomesh)

    ! Grid data
    call WriteGOAT(goatoptions, grid, magneticField, environment)

    ! Fort.78 file with void regions
    call ComputeVoidRegionPolygonSet(grid, topomesh, environment%triangulationvessel, &
        ggtmdata, voidps)
    call WriteVoidRegionFile(voidps, grid, 'fort.78')
    call WriteVoidRegionFileGoat(voidps, grid, 'fort_goat.78')

    ! b2ag file
    !call Writeb2agdat(goatoptions, grid)



end subroutine