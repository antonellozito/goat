subroutine MapVertexPairToFace(v1,v2,fvert,nfvert,faceID)

    ! Description
    !============
    ! This routine maps a vertex pair, given by the vertex indices v1 
    ! and v2, to the face ID that has these vertices (i.e. it retrieves
    ! the face ID with the same vertices). The sequence of the vertices 
    ! doesn't matter. Caution: the routine does not check for multiple 
    ! occurences and simply returns the first face id that matches!

    ! Notes
    !======
    ! Note 1: altough the routine works without sorting the face vertex
    ! array fvert, performance can be improved by doing this in a pre-
    ! processing step, e.g. setting fvert(:,1) = minval(fvert,2) and
    ! fvert(:,2) = maxval(fvert,2) - this is pseudocode!

    ! Note 2: it is assumed that v1 and v2 are scalar integers. 

    ! Note 3: if no face vertex pair matches, faceID is set to -1.

    ! Algorithm
    !==========
    ! First, we loop and check if the pair [v1 v2] is found to be equal 
    ! to any fvert(i,:). If none are found, then we loop to see if we 
    ! find [v2 v1]. If no pairs are found, faceID is set to -1. 

    ! Initialize
    !===========
    ! Declare modules
    use gdmod_types

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    integer(I8)                         :: v1, v2, nfvert
    integer(I8), dimension(nfvert,2)    :: fvert
    integer(I8)                         :: faceID

    ! Loop variables
    integer(I8)                         :: k
    logical                             :: notfound

    ! Auxiliary variables 

    ! Initialize
    !===========
    ! Default faceID output
    faceID = -1

    ! Iteration counter
    k = 1

    ! Loop variable
    notfound = .true.

    ! Loop
    !=====
    ! Is [v1 v2] a pair? 
    do while (notfound .and. k <= nfvert)

        ! Check the first vertex
        if (fvert(k,1) == v1) then 
            ! Check the second vertex
            if (fvert(k,2) == v2) then 
                ! Found edge, exit
                faceID = k
                notfound = .false.
            end if
        end if

        ! Update counter
        k = k+1;

    end do

    ! Reset k
    k = 1

    ! Is [v2 v1] a pair?
    do while (notfound .and. k <= nfvert)

        ! Check the second vertex
        if (fvert(k,2) == v1) then 
            ! Check the first vertex
            if (fvert(k,1) == v2) then 
                ! Found edge, exit
                faceID = k
                notfound = .false.
            end if
        end if

        ! Update counter
        k = k+1;

    end do





end subroutine