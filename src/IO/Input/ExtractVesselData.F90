subroutine ExtractVesselData(vessel, vesseloptions)

    ! Description
    !============
    ! Process the vessel structure to provide the necessary data for the
    ! grid deformation. Here, we construct a polygon set representation 
    ! of the vessel to easier manipulate the vessel structure. It is 
    ! assumed that the vessel structure contains segments that 
    ! properly intersect each other such that a single closed polygon
    ! can be extracted (nested polygons are supported as long as they 
    ! are simply nested, i.e. the interior of the domain can be 
    ! distinguished). 

    ! Notes
    !======
    ! Note 1: (20/03/2024) Added a separate polygon representation 
    ! for each target structure. Each structure that has to be a target
    ! should be identified by the field 'TP' in the vesseloptions 
    ! structure. Each integer i corresponds to the i-th vessel structure
    ! present in vessel%structures. Given that the aim of this 
    ! additional structure is to easily determine distributions etc
    ! around these structures, we artifically close the polygon (if not
    ! already closed), since this allows more types of distribution
    ! functions to be used. Note that self-intersecting polygons are not 
    ! allowed for targets and will cause an error. 

    ! Note 2: (27/05/2024) Added labels in polygon set representation with 
    ! the following format: label(:, 1), label(:, 2) are vessel structure
    ! ID(s) of the vertex (max. 2 allowed per vertex), label(:, 3)  is 
    ! the unique vertex ID (may, after intersections, not go from 1 to 
    ! number of vertices due to exclusion of vertices)

    ! Initialize
    !===========
    ! Declare modules
    use gdmod_types 
    use gdmod_userinput
    use gdmod_interfaces
    use gdmod_plots
    use PolygonLevelsetFunction2D
    use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    type(VesselUDT), intent(inout)      :: vessel
    type(VesselOptionsUDT), intent(in)  :: vesseloptions

    ! Loop variables
    integer(I8)                         :: i, k, status

    ! Auxiliary variables 
    type(PolygonSetUDT)                 :: tempps 
    type(PolygonUDT)                    :: temppol
    integer(I8)                         :: nvp, cc, nexcl, nTP, tne, &
        vID   
    integer(I8), allocatable            :: tv(:), templabels(:, :)
    real(R8), allocatable               :: tempx(:), &
        tempy(:), tx(:), ty(:), xp(:), yp(:), tvx(:), tvy(:), tnxp(:), &
        tnyp(:), tnnp(:), tnx(:), tny(:), tnn(:)

    class(PLF2DOptionsUDT), allocatable :: plfoptions

    ! Plotting

    ! Set NaN
    real(R8)                            :: NaN

    ! Initialize
    !===========
    ! Set NaN
    NaN = IEEE_VALUE(nan, IEEE_QUIET_NAN)
    
    ! Set vertex ID
    vID = 0

    ! Associate
    associate(&
        nvs         => vessel%nstructures,  &
        vs          => vessel%structures)

    ! Vessel polygon representation
    !==============================
    ! Count total number of points currently in the structures
    nvp = 0
    nexcl = 0
    do i = 1, nvs
        if (.not. any(vesseloptions%exclude == i)) then 
            nvp = vessel%structures(i)%np + nvp 
        else 
            nexcl = nexcl + 1
        end if 
    end do

    ! Allocate (account for NaNs)
    allocate(tempx(nvp+nvs-1-nexcl), tempy(nvp+nvs-1-nexcl))
    allocate(templabels(size(tempx), 3))
    templabels = 0
    
    ! Set vertices
    cc = 0
    do i = 1, nvs 
        if (.not. any(vesseloptions%exclude == i)) then 
            ! Coordinates
            tempx(cc+1:cc+vs(i)%np) = vs(i)%x
            tempy(cc+1:cc+vs(i)%np) = vs(i)%y
            
            ! Labels
            templabels(cc+1:cc+vs(i)%np, 1) = vs(i)%ID 
            templabels(cc+1:cc+vs(i)%np, 3) = [(k, k = vID+1, vID+vs(i)%np)]

            ! Update counters
            cc = cc + vs(i)%np
            vID = vID + vs(i)%np 

            ! Add NaN
            if (i .ne. nvs) then 
                tempx(cc+1) = NaN 
                tempy(cc+1) = NaN 
                templabels(cc+1, :) = 0
                cc = cc + 1
            end if 
        end if
    end do

    ! Construct temporary vessel polygon set
    call tempps%Construct(tempx, tempy, templabels)

    ! Construct vessel polygon set 
    call ConstructVesselPolygonSet(vessel, tempps)

    ! Target polygon representation
    !==============================
    ! Check how many targets there are
    nTP = size(vesseloptions%TP)

    ! If there are targets, check that TP does not exceed the number of
    ! structures
    if (nTP > 0) then 
        if (any(vesseloptions%TP > nvs)) then 
            ! Throw error
            call gdErrorHandler('ExtractVesselData: some target plate ' // &
                'indices are larger than number of vessel structures, ' // &
                'check input')
        end if
    end if

    ! Allocate
    allocate(vessel%targetpolygons(nTP))

    ! Construct polygons
    do i = 1, nTP

        ! Construct test polygon
        call temppol%Construct(vs(vesseloptions%TP(i))%x, vs(vesseloptions%TP(i))%y)

        ! Check
        if (temppol%selfintersecting) then 
            ! Throw error
            call gdErrorHandler('ExtractVesselData: target polygon is self intersecting, not supported')
        end if 

        ! Check if we need to close the polygon
        if (.not. temppol%isclosed) then 
            ! Display
            print *, 'ExtractVesselData: closing target polygon with vessel structure ID: ', vesseloptions%TP(i) 

            ! Get points in vertex order
            tne = temppol%ne 
            tv = temppol%vert 
            tvx = temppol%x(tv)
            tvy = temppol%y(tv)

            ! Compute normals in points
            tnx = -(tvy(2:tne+1) - tvy(1:tne))
            tny = (tvx(2:tne+1) - tvx(1:tne))
            tnn = sqrt(tnx**2 + tny**2)
            tnx = tnx/tnn 
            tny = tny/tnn
            tnxp = [tnx(1), 0.5*(tnx(1:tne-1)+tnx(2:tne)), tnx(tne)]
            tnyp = [tny(1), 0.5*(tny(1:tne-1)+tny(2:tne)), tny(tne)]
            tnnp = sqrt(tnxp**2 + tnyp**2)
            tnxp = tnxp/tnnp 
            tnyp = tnyp/tnnp 

            ! Shift the points slightly 
            tx = [tvx, tvx(size(tv)-1:1:-1)+1e-5*tnxp(size(tv)-1:1:-1), tvx(1)]
            ty = [tvy, tvy(size(tv)-1:1:-1)+1e-5*tnyp(size(tv)-1:1:-1), tvy(1)]

            call Write2DPolygonData(tx, ty, 'testpolyg')

            ! Close the polygon
            !tv2 = [tv, tv(size(tv)-1:1:-1)]

            ! Get new coordinates
            !tx = temppol%x(tv2)
            !ty = temppol%y(tv2)

            ! Construct the polygon
            call temppol%Deallocate()
            call temppol%Construct(tx, ty)

        end if 

        ! Assign
        vessel%targetpolygons(i) = temppol

        ! Deallocate test polygon
        call temppol%Deallocate()


    end do 

    ! Construct target polygonset
    call vessel%targetps%Construct(vessel%targetpolygons)

    ! Construct polygon representations
    !==================================
    ! Check how to construct
    select case (trim(vesseloptions%shapemeth))

    case ('polygon')

        ! Exact polygon representation
        allocate(PLF2DGeneralOptionsUDT::plfoptions)

        ! Set options (nothing to do here)

    case ('closedpolygon_exact')

        ! Exact representation of closed polygon
        allocate(PLF2DClosedExactOptionsUDT::plfoptions)

        ! Set options (nothing to do here)

    case ('closedpolygon_smoothapproximation')

        ! Approximate representation of closed polygon
        allocate(PLF2DClosedApproximationOptionsUDT::plfoptions, stat=status)


        ! Set options
        select type (plfoptions)

        type is (PLF2DClosedApproximationOptionsUDT) 

            ! Interpolation settings
            plfoptions%meth = 'uniformgrid'
            plfoptions%resx = vesseloptions%resx
            plfoptions%resy = vesseloptions%resy
            plfoptions%offsetx = vesseloptions%offsetfracx
            plfoptions%offsety = vesseloptions%offsetfracy
            plfoptions%C = vesseloptions%interpC
            plfoptions%M = vesseloptions%interpM
            plfoptions%xrange = vesseloptions%xrange
            plfoptions%yrange = vesseloptions%yrange

            if (size(plfoptions%xrange, 1) < 2) then 
                ! Reset
                call vessel%polygonset%GetVertices(xp, yp)
                plfoptions%xrange = xp 
                plfoptions%yrange = yp
            end if

        end select

    case default 

        ! Throw error
        call gdErrorHandler('ExtractVesselData: vessel shapemeth not implemented')

    end select

    ! Construct
    call InitializePolygonLevelsetFunction2D(vessel%plfvessel, vessel%polygonset, plfoptions)
    call InitializePolygonLevelsetFunction2D(vessel%plftarget, vessel%targetps, plfoptions)

    end associate 

end subroutine