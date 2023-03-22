subroutine ConstructVesselPolygonSet(vessel, vesseloptions, ps)

    ! Description
    !============
    ! ConstructVesselPolygonSet constructs a closed vessel polygon set (possibly
    ! consisting of multiple closed polygons) starting from an initial set of
    ! polygons, given in the oldvessel input argument (see later for input
    ! specifications). Options can be parsed through the vesseloptions input
    ! structure. The routine has support for multiple open and closed polygon
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

    ! Modules
    !========
    use gdmod_types 
    use gdmod_userinput

    ! Declare variables
    !==================
    ! Arguments
    type(VesselUDT)             :: vessel
    type(VesselOptionsUDT)      :: vesseloptions 
    type(PolygonSetUDT)         :: ps 

    ! Auxiliary
    real(R8), allocatable       :: xi(:), yi(:)
    integer(I8), allocatable    :: p1(:), p2(:), s1(:), s2(:)

    ! Loop

    ! Initialize
    !===========
    ! Checks
    if (size(vesseloptions%TPind) == 0 ) then 
        call gdErrorHandler('ConstructVesselPolygon: no target plates are specified, check input')
    elseif (size(vesseloptions%TPind) .ne. size(vesseloptions%TP)) then 
        call gdErrorHandler('ConstructVesselPolygon: number of ' &
            // 'elements in vesseloptions%TPind does not correspond' &
            // ' to number of structures in vesseloptions%TP')
    end if 

    ! Determine intersections
    !========================
    call ps%SelfIntersections(xi, yi, p1, p2, s1, s2)
    print *, xi



end subroutine