subroutine SortPolygonEdges(pein,ne,sortindex,ispolygonstart)

    ! Description
    !===========
    ! This routine sorts the edges of a simply polygon. The definition
    ! that is used here for simple polygon demands that the polygon 
    ! either closes perfectly on itself, or is a non-branching polygon 
    ! (i.e. each vertex only has maximal two edges). This is checked 
    ! during the routine execution when determining the polygon edges. 
    ! Multiple open and closed polygons are supported. The logical 
    ! 'ispolygonstart' indicates which of the (sorted!) edges is the 
    ! start of a new polygon. 

    ! Arguments
    !==========
    !
    ! - pein:           (input) ne-by-2 array of vertex indices that has
    !                   to be sorted. 
    ! - ne:             (input) number of edges (integer)
    ! - sortindex:      (output) index that contains the sequence of 
    !                   sorted polygon edge indices, i.e. 
    !                   pein(sortindex,:) should sort the edges.
    ! - ispolygonstart: (output) ne-by-1 logical of which the i'th
    !                   element is true if sortindex(i) is the start 
    !                   index of a new polygon. The edges of this 
    !                   polygon are all edges between this true value 
    !                   and the next. 

    ! Algorithm
    !==========
    ! 1) Check if the polygon conforms to our demands

    ! Initialize
    !===========
    ! Modules
    use gdmod_types

    ! The usual
    implicit none

    ! Declare variables
    !==================
    ! Input
    integer(I8), dimension(ne,1:2)  :: pein ! polygon edges 
    integer                         :: ne
    
    ! Output
    integer(I8), dimension(ne)      :: sortindex 
    logical, dimension(ne)          :: ispolygonstart

    ! Mixed

    ! Loop
    logical                 :: allfound, startfound, polygonfound
    integer(I8)             :: i, k, spind

    ! Auxiliary
    integer(I8)                 :: nv, nremedges, nextremedge, &
                                tc1, tc2

    logical, allocatable        :: isedgesorted(:), isremedgesorted(:), &
                                mask(:)

    integer(I8), allocatable    :: remedges(:,:), edgeID(:), &
                                remedgeID(:), temparray(:)

    ! Main program
    !=============
    ! Check
    if (size(pein,2) /= 2) then
        ! Throw error
        call gdErrorHandler('SortPolygonEdges: input argument pein should be a ne-by-2 integer array')

    end if

    ! Initialize
    allocate(isedgesorted(ne)) ! logical to indicate if edge has been sorted
    allocate(edgeID(ne))

    ispolygonstart(:) = .false.
    isedgesorted(:) = .false.
    edgeID(:) = (/ (i, i=1,ne,1) /)
    allfound = .false. ! while loop variable
    spind = 1 ! sorted polygon index

    ! Loop
    do while (allfound .eqv. .false.)
        ! Set the polygon starting index
        ispolygonstart(spind) = .true.

        ! Allocate inner loop variables
        nremedges = count(isedgesorted .eqv. .false.)
        allocate(remedges(nremedges,2))
        allocate(remedgeID(nremedges))
        allocate(isremedgesorted(nremedges))
        allocate(mask(nremedges))

        ! Get the remaining edges
        remedges(:,1) = pack(pein(:,1), (isedgesorted .eqv. .false.))
        remedges(:,2) = pack(pein(:,2), (isedgesorted .eqv. .false.))
        remedgeID(:) = pack(edgeID, (isedgesorted .eqv. .false.))
        isremedgesorted(:) = .false.

        ! Find a starting vertex
        startfound = .false. 
        k = 1
        do while ((startfound .eqv. .false.) .and. (k <= nremedges))
            ! Count how many times the current edge vertices occur
            tc1 = count(remedges(:,1) == remedges(k,1)) & 
                + count(remedges(:,2) == remedges(k,1))
            tc2 = count(remedges(:,1) == remedges(k,2)) & 
                + count(remedges(:,2) == remedges(k,2))

            ! Check & add
            if (tc1 == 1) then
                ! Found a start vertex
                startfound = .true.
                nv = remedges(k,2) ! next vertex

                ! Add the current edge
                sortindex(spind) = remedgeID(k)

                ! Set edge as sorted
                isremedgesorted(k) = .true.

                ! Update indices
                k = k+1
                spind = spind+1
            else if (tc2 == 1) then
                ! Found a start vertex
                startfound = .true.
                nv = remedges(k,1) ! next vertex

                ! Add the current edge
                sortindex(spind) = remedgeID(k)

                ! Set edge as sorted
                isremedgesorted(k) = .true.

                ! Update counters
                k = k+1
                spind = spind+1
            else if ((tc1 > 2) .or. (tc2 > 2)) then
                ! Polygon branches, throw error
                call gdErrorHandler('SortPolygonEdges: branching polygon detected, not supported')
            else
                ! Next edge, update counter
                k = k+1
            end if
        end do

        ! Hedge for closed polygon(s)
        if (startfound .eqv. .false.) then
            ! No start index found, yet edges remain -> has to be closed
            ! polygon. Simply take the first vertex of the first edge
            startfound = .true.
            nv = remedges(1,2)

            ! Add the current edge
            sortindex(spind) = remedgeID(1)

            ! Set edge as sorted
            isremedgesorted(1) = .true.

            ! Update counters
            spind = spind+1
        end if

        ! Find the edges of this polygon
        polygonfound = .false.
        do while (polygonfound .eqv. .false.)
            ! The next edge has the current vertex and is not yet found.
            mask(:) = (isremedgesorted .eqv. .false.) ! first requirement
            mask = mask .and. & 
                ((remedges(:,1) == nv) .or. (remedges(:,2) == nv)) ! second requirement
            
            ! Check if this is the final edge
            if (count(mask) == 0) then
                ! All edges were found, exit
                polygonfound = .true. 
            else if (count(mask) > 1) then
                ! Unknown error, call error handler
                call gdErrorHandler('SortPolygonEdges: branching polygon detected, not supported')
            else
                ! Get the next edge
                allocate(temparray(1)) ! avoid rank conflicts
                temparray = pack((/ (i, i=1,nremedges,1) /),mask)
                nextremedge = temparray(1)

                ! Add edge
                sortindex(spind) = remedgeID(nextremedge)

                ! Set edge as sorted
                isremedgesorted(nextremedge) = .true.

                ! Get next vertex
                if (remedges(nextremedge,1) == nv) then
                    nv = remedges(nextremedge,2)
                else
                    nv = remedges(nextremedge,1)
                end if
                
                ! Update counters
                spind = spind+1

                ! Deallocate
                deallocate(temparray)
            end if 
        end do

        ! Update logicals
        isedgesorted(pack(remedgeID,isremedgesorted)) = .true. 

        ! Check if all edges have been found
        if (spind == ne+1) then ! one ahead due to updating rule here
            allfound = .true.
        end if

        ! Deallocate inner loop variables
        deallocate(remedges)
        deallocate(remedgeID)
        deallocate(isremedgesorted)
        deallocate(mask)

    end do

end subroutine