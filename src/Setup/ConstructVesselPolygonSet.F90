subroutine ConstructVesselPolygonSet(vessel, ps)

    ! Description
    !============
    ! ConstructVesselPolygonSet constructs a closed vessel polygon set (possibly
    ! consisting of multiple closed polygons) starting from an initial set of
    ! polygons, given in the oldvessel input argument (see later for input
    ! specifications). The routine has support for multiple open and closed polygon
    ! boundaries, though the closed polygons (or the closed polygons resulting
    ! from combining multiple open polygons) should not intersect themselves
    ! with other open or closed polygons. Otherwise, an error will be thrown.
    ! See the notes for more explanation on how closed polygons are identified
    ! and how open polygons are merged to create a closed one. 

    ! Algorithm
    !==========
    ! 1) Check which structures are closed. For these closed structures, check
    ! if they have any intersections with other structures -> throw error if
    ! this is the case
    ! 2) Compute intersections of the remaining open polygons with all other
    ! polygons. If not exactly two intersections are found, an error is thrown,
    ! as this is indicates that the input assumptions are violated.
    ! 3) Construct the different closed polygon parts by looking at the
    ! intersections of each polygon.
    ! 4) Construct the new vessel polygon(s) by concatenating the polygon
    ! segments that were found. This is also the stage where target plate
    ! indices etc are constructed. 

    ! Notes
    !------
    ! Note 1: closed polygons are identified by comparing their end points. 
    ! These should be the same up to 100*macheps (to get rid of possible 
    ! numerical artefacts). This is measured in terms of Euclidean distance
    ! between the points (i.e. dist <= 100*macheps

    ! Note 2: open polygons can either truly intersect with each other at
    ! certain edges or coincide in one of their end vertices. The latter is
    ! checked by taking the distance between the end points and comparing
    ! to 100*macheps (same criterion as above for determining closed polygons).
    ! For polygon intersections, the end sections are cut off. Intersection
    ! computation is done using the 'intersections' routine

    ! Note 3: the resulting vessel polygon(s) is(are) arbitrarily sorted
    ! (either clockwise or counterclockwise). 

    ! Note 4: we propagate the structure number(s) of each vertex and 
    ! the vertex IDs as they were present in the original vessel 
    ! structure. Note that vertices are only allowed up to two structure
    ! IDs and only one vertex ID. Vertex IDs will be removed since 
    ! either the structures intersect exactly in the end points and only
    ! one node is retained, or the structures intersect somewhere along
    ! the polygon and the 'dangling' nodes are removed. It is assumed 
    ! that this information is already available through the polygon
    ! labels field, which should be setup in ExtractVesselData.F90. 

    ! Modules
    !========
    use gdmod_types 
    use gdmod_userinput
    use gdmod_plots

    implicit none

    ! Declare variables
    !==================
    ! Arguments
    type(VesselUDT)             :: vessel
    type(PolygonSetUDT)         :: ps 

    ! Auxiliary
    integer(I8)                 :: ni, nfinpol, thisp, firstpolygon, &
        c1, c2, nvest, nvv, tempnv, nextp, sv, ev, flag, vID
    logical                     :: polygonnotfound, doflip  
    real(R8)                    :: nan, xs, ys, xe, ye 

    real(R8), allocatable       :: xi(:), yi(:), tempx(:), tempy(:), &
        xv(:), yv(:)
    integer(I8)                 :: pis(1:2), pie(1:2)
    integer(I8), allocatable    :: p1(:), p2(:), s1(:), s2(:), &
        polcat(:, :), npol(:), remp(:), tempp(:), indi(:, :), &
        si(:, :), pi(:, :), ci(:), templabels(:, :), &
        labelsv(:, :)
    logical, allocatable        :: notfound(:)

    character(:), allocatable   :: vesselpath

    ! Loop
    integer(I8)                 :: i, j, k

    ! Initialize
    !===========
    ! Checks - seem unnecessary?
    !if (size(vesseloptions%TPind) == 0 ) then 
    !    call gdErrorHandler('ConstructVesselPolygon: no target plates are specified, check input')
    !elseif (size(vesseloptions%TPind) .ne. size(vesseloptions%TP)) then 
    !    call gdErrorHandler('ConstructVesselPolygon: number of ' &
    !        // 'elements in vesseloptions%TPind does not correspond' &
    !        // ' to number of structures in vesseloptions%TP')
    !end if 

    ! Set NaN
    nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)

    ! Get maximal vertex ID label - if available
    vID = 0
    if (size(ps%polygons(1)%labels, 2) >= 3) then 
        do i = 1, ps%np 
            vID = max(maxval(ps%polygons(i)%labels(:, 3)), vID)
        end do 
    else
        ! Throw warning
        print *, 'ConstructVesselPolygonSet: polygons do not have ' // & 
            'all labels, constructing vessel structure and vertex labels ' // & 
            'based on polygons'
        do i = 1, ps%np
            deallocate(ps%polygons(i)%labels)
            allocate(ps%polygons(i)%labels(ps%polygons(i)%nv, 3))
            ps%polygons(i)%labels(:, 1) = i 
            ps%polygons(i)%labels(:, 2) = 0
            ps%polygons(i)%labels(:, 3) = [(k, k = vID+1, vID+ps%polygons(i)%nv)]
            vID = vID + ps%polygons(i)%nv 
        end do 
    end if 

    ! Determine intersections
    !========================
    ! Compute all intersections
    call ps%SelfIntersections(xi, yi, p1, p2, s1, s2)
    
    ! Check if each polygon has either exactly zero or two intersections
    do i = 1, ps%np 
        ! Count
        ni = count(i == p1) + count(i == p2)

        ! Check
        if ( ni == 1 ) then 
            ! This shouldn't be happening with the right input format
            ! -> probably something wrong on the input side 
            print *, 'polygon number (excluded polygons not accounted):'
            print *, i
            call gdErrorHandler('ConstructVesselPolygon: only one ' &
                // 'intersection found for this polygon, check input')
        elseif ( ni > 2 ) then 
            ! This shouldn't be happening with the right input format
            ! -> probably something wrong on the input side 
            print *, 'polygon number (excluded polygons not accounted):'
            print *, i
            call gdErrorHandler('ConstructVesselPolygon: too many ' &
                // 'intersections found for this polygon, check input')
        elseif ( ni == 0 ) then 
            if (.not. ps%polygons(i)%isclosed) then 
                ! Open polygons should intersect twice, this shouldn't 
                ! be happening
                print *, 'polygon number (excluded polygons not accounted):'
                print *, i
                call gdErrorHandler('ConstructVesselPolygon: not enough ' &
                    // 'intersections found for this polygon, check input')
            end if 
        elseif ( (ni .ne. 0) .and. ps%polygons(i)%isclosed) then 
            ! Closed polygons should not intersect with itself (this is
            ! topologically possible but is not allowed here)
            print *, 'polygon number (excluded polygons not accounted):'
            print *, i
            print *, ps%polygons(i)%edges
            call gdErrorHandler('ConstructVesselPolygon: closed ' &
                // 'polygon found that intersects with itself, check input')

        end if

    end do

    ! Compute closed polygons
    !========================
    ! First, we check which polygons form one single closed polygon and
    ! how many closed polygons there are. 

    ! Allocate
    allocate(notfound(ps%np), polcat(ps%np, ps%np), &
        npol(ps%np))

    ! Initialize
    notfound(:) = .true. 
    polcat(:, :) = 0 ! concatenation structure: each column holds the polygon indices of a closed polygon
    npol(:) = 0 ! counter for how many polygons belong to each final polygon
    nfinpol = 0 ! number of final polygons

    ! Loop
    nfinpol = 0
    do while (any(notfound))
        ! Get the next polygon piece
        allocate(remp(count(notfound))) 
        remp = pack([(k, k=1, ps%np)], notfound)
        thisp = remp(1) 

        ! If polygon is closed -> add
        if (ps%polygons(thisp)%isclosed) then 

            ! Update total number of closed polygons
            nfinpol = nfinpol + 1
           
            ! Set logical to false 
            notfound(thisp) = .false.

            ! Update polygon structure counter
            npol(nfinpol) = npol(nfinpol) + 1 

            ! Add to polygon structure
            polcat(npol(nfinpol), nfinpol) = thisp 

        else 
            ! Polygon is not closed, need to loop untill we find all
            ! polygons for this closed polygon

            ! Update total number of closed polygons
            nfinpol = nfinpol + 1

            ! Set this polygon as the first
            firstpolygon = thisp 

            ! Loop
            polygonnotfound = .true. 
            do while (polygonnotfound) 
                
                ! Save the current polygon
                npol(nfinpol) = npol(nfinpol) + 1 
                polcat(npol(nfinpol), nfinpol) = thisp
                notfound(thisp) = .false. 

                ! Get next polygon candidates
                c2 = count( (thisp == p1) .and. (notfound(p2)) ) ! don't include polygons that were found already
                c1 = count( (thisp == p2) .and. (notfound(p1)) ) 

                ! Check the amount of candidates
                if ( (c2 + c1) == 1) then 
                    ! Just one candidate - ok 
                    if (c1 == 0) then 
                        ! Update thisp 
                        allocate(tempp(c2))
                        tempp = pack(p2, (thisp == p1) .and. (notfound(p2)))
                        thisp = tempp(1)
                        deallocate(tempp)
                    else 
                        ! Update thisp 
                        allocate(tempp(c1))
                        tempp = pack(p1, (thisp == p2) .and. (notfound(p1)))
                        thisp = tempp(1)
                        deallocate(tempp)
                    end if  
                elseif ( (c2 + c1) == 2) then 
                    ! Two candidates - only possible if this is the 
                    ! first polygon.
                    if (thisp .ne. firstpolygon) then  
                        ! Throw error
                        print *, 'polygon number: ', thisp
                        call gdErrorHandler('ConstructVesselPolygon: ' &
                            // 'multiple polygons found for polygon ' &
                            // 'that is not the starting polygon. ' &
                            // 'Check vessel input') 
                    end if 

                    ! Simply take one of the two
                    if (c1 > 0) then 
                        ! Update thisp 
                        allocate(tempp(c1))
                        tempp = pack(p1, (thisp == p2) .and. (notfound(p1)))
                        thisp = tempp(1)
                        deallocate(tempp)
                    else 
                        ! Update thisp 
                        allocate(tempp(c2))
                        tempp = pack(p2, (thisp == p1) .and. (notfound(p2)))
                        thisp = tempp(1)
                        deallocate(tempp)
                    end if 
                elseif ( (c2 + c1) == 0) then 
                    ! We should have reached the end of the polygon here
                    ! - need to check if the next polygon here would be 
                    ! the first one again
                    c2 = count( (thisp == p1) .and. (p2 == firstpolygon) )
                    c1 = count( (thisp == p2) .and. (p1 == firstpolygon) ) 
                    if ( (c1 + c2) == 1) then 
                        ! Ok, exit the loop
                        polygonnotfound = .false. 
                    else 
                        ! Not OK, throw error - this means an open 
                        ! vessel structure was found
                        call gdErrorHandler('ConstructVesselPolygon: ' &
                            // 'Open vessel structure detected, ' &
                            // 'check vessel input') 
                    end if 
                elseif ( (c2 + c1) > 2) then 
                    ! This should actually not happen anymore since this
                    ! is checked upfront, but still we handle this error
                    ! should it occur due to unknown causes. 
                    call gdErrorHandler('ConstructVesselPolygon: ' &
                        // 'more than two polygons found that ' &
                        // 'intersect with this polygon, ' &
                        // 'check vessel input')
                end if
            end do
        end if 

        ! Housekeeping
        deallocate(remp)
    end do 

    ! Construct new vessel
    !=====================
    ! Here, we 'merge' all the polygon pieces together for each closed
    ! polygon and reconstruct the vessel polygon set. To this end, we
    ! must make sure that the polygon pieces are properly trimmed at the
    ! intersection points and that the order of the points is correct as
    ! well. Since it is tedious (but not impossible) to precompute the 
    ! final total number of polygon vertices in order to allocate the 
    ! coordinate arrays, we make an overestimation and allocate too 
    ! much memory to trim later. 
    !
    ! The segments are sorted as follows: for each segment, we check
    ! which intersection corresponds to the intersection with the next
    ! polygon in the polcat structure. If that intersection has the 
    ! highest segment index of both intersections in that polygon, the 
    ! polygon is correctly oriented. Otherwise, we need to reverse the
    ! vertex order of this polygon before adding it to the vessel 
    ! structure. 

    ! First, cast the intersection data in easier to navigate 
    ! structures. We can safely assume here that there are either two 
    ! or zero intersections per polygon. We store:
    ! - pi: polygon indices of the intersections
    ! - si: section indices of the intersection of the polygon
    ! - indi: the index of the intersection (to retrieve the x- and y- 
    ! coordinates of this intersection as xi(indi) ) 
    ! - ci is a counter
    allocate(indi(2, ps%np), si(2, ps%np), pi(2, ps%np), ci(ps%np))
    indi(:, :)  = 0
    si(:, :)    = 0
    pi(:, :)    = 0
    ci(:)       = 0

    do i = 1, size(xi)
        ! First polygon contribution
        thisp = p1(i)
        ci(thisp) = ci(thisp) + 1
        indi(ci(thisp), thisp) = i
        si(ci(thisp), thisp) = s1(i)
        pi(ci(thisp), thisp) = p2(i)
        
        ! Second polygon contribution
        thisp = p2(i)
        ci(thisp) = ci(thisp) + 1
        indi(ci(thisp), thisp) = i
        si(ci(thisp), thisp) = s2(i)
        pi(ci(thisp), thisp) = p1(i)
    end do

    ! Make estimate of number of vessel vertices
    nvest = 0
    do i = 1, ps%np 
        nvest = nvest + ps%polygons(i)%nv 
    end do 
    nvest = nvest*2 ! factor 2 just to be sure 

    ! Allocate
    allocate(tempx(nvest), tempy(nvest), templabels(nvest, 3))

    ! Loop 
    nvv = 0 ! vessel vertex counter
    do i = 1, nfinpol 
        ! Check number of polygon pieces
        if (npol(i) == 1) then ! Closed polygon piece
            ! Get the polygon number
            thisp = polcat(1, i)

            ! Get number of vertices of this polygon
            ! Note: we need to work through the %vert structure, which 
            ! has size ne+1! Otherwise, results may be wrong 
            tempnv = ps%polygons(thisp)%ne+1

            ! Add coordinates
            tempx(nvv+1:nvv+tempnv) = ps%polygons(thisp)%x(ps%polygons(thisp)%vert) 
            tempy(nvv+1:nvv+tempnv) = ps%polygons(thisp)%y(ps%polygons(thisp)%vert)  

            ! Add labels
            templabels(nvv+1:nvv+tempnv, :) = ps%polygons(thisp)%labels(ps%polygons(thisp)%vert, :) 

            ! Update counter
            nvv = nvv + tempnv 

        else ! Loop over all segments
            do j = 1, npol(i)
                ! Get the current polygon
                thisp = polcat(j, i) 

                ! Check if the intersection with the largest si 
                ! corresponds to the intersection with the next polygon.
                ! If so, add the polygon as is. Otherwise, flip the 
                ! coordinates.

                ! Get the next polygon
                if (j < npol(i)) then 
                    nextp = polcat(j+1, i)
                else 
                    nextp = polcat(1, i)
                end if 

                ! Get start and end indices and check if we need to 
                ! flip. Set start and end coordinates 
                doflip = .false. 
                if (pi(1, thisp) == nextp) then 
                    ! First intersection, check si
                    if (si(1, thisp) >= si(2, thisp) ) then 
                        ! Don't flip
                        xs = xi(indi(2, thisp)) ! start 
                        ys = yi(indi(2, thisp)) 
                        pis = [p1(indi(2, thisp)), p2(indi(2, thisp))]
                        xe = xi(indi(1, thisp)) ! end 
                        ye = yi(indi(1, thisp)) 
                        pie = [p1(indi(1, thisp)), p2(indi(1, thisp))]
                        sv = si(2, thisp) + 1
                        ev = si(1, thisp) 
                    else 
                        ! Flip
                        doflip = .true. 
                        xs = xi(indi(2, thisp)) ! start 
                        ys = yi(indi(2, thisp)) 
                        pis = [p1(indi(2, thisp)), p2(indi(2, thisp))]
                        xe = xi(indi(1, thisp)) ! end 
                        ye = yi(indi(1, thisp)) 
                        pie = [p1(indi(1, thisp)), p2(indi(1, thisp))]
                        sv = si(1, thisp) + 1
                        ev = si(2, thisp) 
                    end if 
                elseif ( pi(2, thisp) == nextp) then 
                    ! Second intersection, check si
                    if (si(2, thisp) >= si(1, thisp) ) then 
                        ! Don't flip
                        xs = xi(indi(1, thisp)) ! start 
                        ys = yi(indi(1, thisp)) 
                        pis = [p1(indi(1, thisp)), p2(indi(1, thisp))]
                        xe = xi(indi(2, thisp)) ! end 
                        ye = yi(indi(2, thisp)) 
                        pie = [p1(indi(2, thisp)), p2(indi(2, thisp))]
                        sv = si(1, thisp) + 1 
                        ev = si(2, thisp) 
                    else 
                        ! Flip
                        doflip = .true. 
                        xs = xi(indi(1, thisp)) ! start 
                        ys = yi(indi(1, thisp)) 
                        pis = [p1(indi(1, thisp)), p2(indi(1, thisp))]
                        xe = xi(indi(2, thisp)) ! end 
                        ye = yi(indi(2, thisp))
                        pie = [p1(indi(2, thisp)), p2(indi(2, thisp))]
                        sv = si(2, thisp) + 1
                        ev = si(1, thisp) 
                    end if 
                else 
                    ! Unexpected error
                    call gdErrorHandler('ConstructVesselPolygon: ' &
                        // 'encountered unexpected error when ' &
                        // 'constructing full vessel polygon')
                end if 

                ! Add
                ! Starting point
                tempx(nvv+1) = xs 
                tempy(nvv+1) = ys 
                templabels(nvv+1, 1:2) = pis
                templabels(nvv+1, 3) = vID+1

                ! Update counters
                nvv     = nvv + 1
                vID     = vID + 1

                ! Polygon points
                tempnv = ev - sv + 1 
                
                if (doflip) then 
                    tempx(nvv+1:nvv+tempnv) = &
                        ps%polygons(thisp)%x(ps%polygons(thisp)%vert(ev:sv:-1))
                    tempy(nvv+1:nvv+tempnv) = &
                        ps%polygons(thisp)%y(ps%polygons(thisp)%vert(ev:sv:-1))
                    templabels(nvv+1:nvv+tempnv, :) = &
                        ps%polygons(thisp)%labels(ps%polygons(thisp)%vert(ev:sv:-1), :)
                else 
                    tempx(nvv+1:nvv+tempnv) = ps%polygons(thisp)%x(ps%polygons(thisp)%vert(sv:ev))
                    tempy(nvv+1:nvv+tempnv) = ps%polygons(thisp)%y(ps%polygons(thisp)%vert(sv:ev))
                    templabels(nvv+1:nvv+tempnv, :) = ps%polygons(thisp)%labels(ps%polygons(thisp)%vert(sv:ev), :)
                end if

                ! Hedge for zero/negative length
                if (tempnv < 0) then 
                    tempnv = 0
                end if
                nvv = nvv+tempnv

                ! End point
                tempx(nvv+1) = xe
                tempy(nvv+1) = ye 
                templabels(nvv+1, 1:2) = pie
                templabels(nvv+1, 3) = vID+1

                ! Update counters
                nvv = nvv + 1 
                vID = vID + 1

            end do 
        end if

        ! Add nan (unless if last polygon)
        if (i .ne. nfinpol) then 
            tempx(nvv+1) = nan 
            tempy(nvv+1) = nan
            templabels(nvv+1, :) = 0 
            nvv = nvv+1
        end if
    end do

    ! Construct vessel polygon set
    xv = tempx(1:nvv)
    yv = tempy(1:nvv)
    labelsv = templabels(1:nvv, :)
    call vessel%polygonset%Construct(xv, yv, labelsv)

    ! Test orientation
    call vessel%polygonset%OrientNestedClosedPolygons(flag)

    ! Check
    if (flag .ne. 0) then  
        ! Throw error
        print *, 'flag: ', flag
        call gdErrorHandler('ConstructVesselPolygon: could not orient polygons, OrientNestedClosedPolygons exited with flag above')
    end if 

    ! Write data
    vesselpath = 'vesselpolygon'
    call vessel%polygonset%WriteData(vesselpath)

    ! Housekeeping
    !=============
    deallocate(indi, si, pi, ci, tempx, tempy) 

end subroutine