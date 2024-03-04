subroutine ExtractPolygonVertices(pe,ne,pv)

    ! Description
    !===========
    ! This routine extracts the vertices of a polygon from a set of 
    ! sorted polygon edges which are given by their vertex IDs in the 
    ! 'pe' (polygon edges) ne-by-2 array. Only single polygons (can be 
    ! closed or open) are supported, i.e. no branching or multiple 
    ! polygons. This is checked for during the routine. 

    ! Arguments
    !==========
    !
    ! - pe:             (input) ne-by-2 array of vertex indices that has
    !                   to be sorted. 
    ! - ne:             (input) number of edges (integer)
    ! - pv:             (output) ne+1-by-1 array of polygon vertices

    ! Algorithm
    !==========
    ! 0)    Initialize & check
    ! 1)    Take the first edge
    !
    ! 2)    Check which vertex is the start vertex by comparing with the 
    !       vertices of the next edge. Set this as the next vertex. 
    !       If (the vertex is not found) then
    !           throw error
    !       else 
    !           add the current vertex
    !           set the other vertex as the next vertex
    !       end if
    !
    ! 3)    for the remaining edges, repeat 2) but with the start vertex
    !       equal to the next vertex. If the last vertex is reached, add 
    !       the remaining vertex. 

    ! Initialize
    !===========
    ! Modules
    use mod_precision

    ! The usual
    implicit none

    ! Declare variables
    !==================
    ! Input
    integer(I8), dimension(ne,1:2)  :: pe ! polygon edges 
    integer                         :: ne
    
    ! Output
    integer(I8), dimension(ne+1)    :: pv

    ! Mixed

    ! Loop
    integer(I8)                     :: k, cvind, nv

    ! Auxiliary
    logical, dimension(1:4)         :: check

    ! Main program
    !=============
    ! Check
    if (size(pe,2) /= 2) then
        ! Throw error
        call gdErrorHandler('ExtractPolygonVertices: input argument pe should be a ne-by-2 integer array')

    end if

    ! Initialize
    check(:) = .false. 
    k = 1

    ! Hedge for the trivial case of one edge
    if (ne == 1) then 
        ! Check 
        if (pe(1,1) == pe(1,2)) then 
            call gdErrorHandler('ExtractPolygonVertices: starting polygon edge vertices are the same - check input')
        else
            pv(1:2) = pe(1,1:2)
        end if 

        ! Exit routine
        return 
    end if

    ! Get the starting vertex & do check on polygon
    check(1) = pe(1,1) == pe(2,1)
    check(2) = pe(1,1) == pe(2,2)
    check(3) = pe(1,2) == pe(2,1)
    check(4) = pe(1,2) == pe(2,2)
    if ( (check(1) .and. (.not. any(check(2:4))) ) .or. &
         (check(2) .and. (.not. any(check((/1, 3, 4/)))) ) ) then 
        ! Second vertex is starting vertex
        cvind = 2
    else if ( (check(4) .and. (.not. any(check(1:3))) ) .or. &
        (check(3) .and. (.not. any(check((/1, 2, 4/)))) ) ) then 
        ! First vertex is starting vertex
            cvind = 1
    else 
        ! Something wrong with input - throw error
        print *, pe(1,:)
        print *, pe(2,:)
        call gdErrorHandler('ExtractPolygonVertices: something wrong with input - please check input polygon')
    end if

    ! Loop
    do while (k < ne)
        ! Get the current vertex, add it to the polygon vertices
        pv(k) = pe(k,cvind)

        ! Get the next vertex
        if (cvind == 1) then
            nv = pe(k,2)
        else
            nv = pe(k,1)
        end if 

        ! Update cvind
        if (nv == pe(k+1,1)) then
            cvind = 1
        else if (nv == pe(k+1,2)) then
            cvind = 2
        else
            ! Throw error, next edge does not contain the vertex
            call gdErrorHandler('ExtractPolygonVertices: could not find next edge - please check input polygon')
        end if

        ! Update counter
        k = k+1

    end do

    ! Add last vertices
    pv(k) = pe(ne,cvind)
    if (cvind == 1) then
        pv(k+1) = pe(ne,2)
    else
        pv(k+1) = pe(ne,1)
    end if

end subroutine