subroutine ComputeGridInterconnections(grid)

    ! Description
    !============
    ! Compute additional grid topology data starting from basic grid
    ! information. The following fields are required (see gdmod_types
    ! for an explanation of the different fields)
    !
    ! Vertices: ntot, BV
    ! Faces: ntot, vert
    ! Cells: ntot, vertP, vert, nvert 
    !   
    ! All other fields are recomputed based on this basic 
    ! interconnection data.             

    ! Notes
    !======
    ! Note 1: guard cells are identified by checking which cells only 
    ! have two vertices. It

    ! Note 2: it is assumed that all fields of which the size is 
    ! predetermined by the number of vertices/faces/cells are already
    ! allocated. Other fields are allocated or overwritten. 

    ! Note 3: the original cell vertices are replaced by a sorted list.
    ! This list is generated based on the faces of the cell and their
    ! corresponding vertices. It is assumed that each cell has the same
    ! amount of faces as it has vertices.

    ! Note 4: the vertex neighbours and vertex cells are ordened in 
    ! clockwise or counter-clockwise manner. This sorting assumes that
    ! each vertex has more than two cells - otherwise, it doesn't really
    ! matter how it is sorted. If it has more than two cells, the cells
    ! have to be properly connected, in the sense that each cell should
    ! only have one neighbour. 

    ! Initialize
    !===========
    ! Declare modules
    use gdmod_types 
    use gdmod_interfaces

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    type(GridUDT)               :: grid

    ! Loop variables
    integer(I8)                 :: i, j, k, m, q

    ! Auxiliary variables 
    type(VertexUDT)             :: v 
    type(FaceUDT)               :: f
    type(CellUDT)               :: c

    integer(I8)                 :: nc, nv, nf, sp, ep, tcs, thiscell, &
        nvc, ntv, tf, ind, tfv(1:2), tc, tcn, ncf, nextcell, tcf, vc, &
        ntcf2, nb1, nb2, nfc, fcc
    integer(I8), allocatable    :: fcount(:), vcount(:), tv(:), &
        tempfcell(:, :), localID(:), fc(:), fn(:), allvertcells(:), &
        allfv(:), tcf2(:), nb(:), vcc(:), indv(:)

    logical, allocatable        :: cellfound(:), allfvind(:)
    logical                     :: startcellnotfound, allnotfound, &
        accountforGC

    ! Unpack & initialize
    !==================
    ! Data structures 
    v = grid%vert
    f = grid%face
    c = grid%cell 

    ! Checks
    if (size(f%vert,2) /= 2) then
        ! This should actually never happen, but just to be sure
        call gdErrorHandler('ComputeGridInterconnections: too many vertices per face')
    end if
    if (any(any(f%vert == 0,2),1)) then
        ! This could be the case when the face vertices are not correctly determined
        call gdErrorHandler('ComputeGridInterconnections: some vertex indices are zero in faces, not supported')
    end if

    ! Quantities
    nv = v%ntot
    nf = f%ntot
    nc = c%ntot

    ! Counters
    allocate(fcount(nf))
    allocate(vcount(nv))

    ! Face neighbours, vertex cells, cell faces
    !==========================================
    ! Check allocation
    if (allocated(c%face)) then
        ! Print out warning
        print *, 'ComputeGridInterconnections: recomputing cell%face'

        ! Deallocate
        deallocate(c%face)
    end if
    if (allocated(v%cell)) then 
        ! print out warning
        print *, 'ComputeGridInterconnections: recomputing vert%cell'

        ! Deallocate
        deallocate(v%cell)
    end if
    if (allocated(f%cell)) then 
        ! Print out warning
        print *, 'ComputeGridInterconnections: recomputing face%cell'

        ! Deallocate
        deallocate(f%cell)
    end if
    if (allocated(v%cell)) then 
        ! Print out warning
        print *, 'ComputeGridInterconnections: recomputing vert%cell'

        ! Deallocate
        deallocate(v%cell)
    end if

    ! Check which cells are guard cells
    c%GC(:) = .false. ! initialize to false
    do i = 1, nc 
        if (c%vertP(i,2) <= 2) then
            ! Guard cell found, set true
            c%GC(i) = .true.
        end if
    end do

    ! Check if we have to account for guard cells
    accountforGC = .false.
    if (any(c%GC)) then
        accountforGC = .true.
    end if 

    ! Make the pointer list for cell faces
    c%faceP(1,1) = 1
    if (c%GC(1)) then 
        ! Only one face
        c%faceP(1,2) = 1
    else
        ! Same amount of faces as vertices 
        c%faceP(1,2) = c%vertP(1,2) 
    end if
    

    do i = 2, nc
        ! Pointer
        c%faceP(i,1) = c%faceP(i-1,1) + c%faceP(i-1,2)

        ! Number of faces
        if (c%GC(i)) then 
            c%faceP(i,2) = 1
        else
            c%faceP(i,2) = c%vertP(i,2)
        end if

    end do

    ! Compute the total faces in cells%face
    c%nface = sum(c%faceP(:,2)) 

    ! Allocate
    allocate(c%face(c%nface))

    ! Make the pointer list for vertex cell neighbours
    v%cellP(:,:) = 0
    do i = 1, nc
        ! Get the number of vertices of the current cell
        ntv = c%vertP(i,2)

        ! Allocate
        allocate(tv(ntv))

        ! Get the vertices of the current cell
        tv = c%vert(c%vertP(i,1):(c%vertP(i,1)+ntv-1))

        ! Update the counter
        v%cellP(tv,2) = v%cellP(tv,2)+1

        ! Deallocate
        deallocate(tv)

    end do

    ! Compute the total number of cells in vert%cell
    v%ncell = sum(v%cellP(:,2))

    ! Allocate
    allocate(v%cell(v%ncell))

    ! Construct the pointer
    v%cellP(1,1) = 1
    do i = 2, nv 
        v%cellP(i,1) = v%cellP(i-1,1) + v%cellP(i-1,2)
    end do

    ! Set vertex cell and cell face and face neighbours
    fcount(:) = 1
    vcount(:) = 1
    allocate(tempfcell(f%ntot, 2))
    tempfcell(:, :) = 0
    do i = 1, nc ! loop over all cells
        ! Get vertices of this cell
        tv = GetCellVert(c, i)
        ntv = c%vertP(i,2)

        ! Loop over all vertices
        do j = 1, ntv

            ! Get the faces of the current cell
            ! Get the next face 
            if (j == ntv) then
                k = 1
            else
                k = j+1
            end if 

            call MapVertexPairToFace(tv(j), tv(k), f%vert, nf, tf)

            ! Check
            if (tf < 0) then
                ! This shouldn't happen
                print *, 'vertex pair: ', tv(j), tv(k)
                call gdErrorHandler(&
                    'ComputeGridInterconnections: could not find '&
                    // 'matching face for this vertex pair')
            else

                ! Add the current cell index as neighbour to this face
                if (.not. ((j > 1) .and. c%GC(i))) then ! hedge for guard cells
                    ! Sanity check
                    if (fcount(tf) > 2) then
                        print *, 'face ID: ', tf
                        print *, 'neighbours: ', tempfcell(tf,:), i
                        call gdErrorHandler(& 
                        'ComputeGridInterconnections: too many ' &
                            // 'neighbours for this face')
                    end if

                
                    tempfcell(tf,fcount(tf)) = i

                    ! Update fcount
                    fcount(tf) = fcount(tf)+1
                end if

                ! Add the current cell to the j'th vertex
                ind = v%cellP(tv(j),1) + vcount(tv(j)) - 1
                if (ind > v%ncell) then
                    call gdErrorHandler('unknown error')
                end if
                v%cell(ind) = i

                ! Update vcount
                vcount(tv(j)) = vcount(tv(j))+1

                ! Add cell face
                if (.not. ((j > 1) .and. c%GC(i))) then ! hedge for guard cells
                    ind = c%faceP(i,1) + j - 1
                    if (ind > c%nface) then
                        call gdErrorHandler('unknown error')
                    end if
                    c%face(ind) = tf
                end if
            end if

        end do

        ! Deallocate
        deallocate(tv)

    end do

    ! Construct cell arrays for faces and vertices
    fcount = fcount-1
    f%cellP(:, 2) = fcount 
    f%cellP(1, 1) = 1
    f%ncell = sum(f%cellP(:, 2))
    allocate(f%cell(f%ncell))
    fcc = 0
    nfc = fcount(1)
    f%cell(fcc+1:fcc+nfc) = tempfcell(1, 1:nfc)
    fcc = fcc + nfc
    do i = 2, f%ntot
        f%cellP(i, 1) = f%cellP(i-1, 2) + f%cellP(i-1, 1)
        nfc = fcount(i)
        f%cell(fcc+1:fcc+nfc) = tempfcell(i, 1:nfc)
        fcc = fcc + nfc
    end do

    ! Housekeeping
    deallocate(vcount, fcount)
    
    ! Logicals
    !=========
    ! Based on the interconnections computed above, the logical indices
    ! (boundary vertices => BV, boundary faces => BF) etc can be 
    ! computed. Boundary faces are either: 1) faces with a single 
    ! neighbour (if no guard cells are present), 2) faces with a guard
    ! cell neighbour (if guard cells are present)

    ! Compute BF & BV
    f%BF(:) = .false. ! initialize
    v%BV(:) = .false. ! initialize
    do i = 1, f%ntot

        ! Get face neighbours
        nb  = GetFaceCell(f, i)
        if (size(nb, 1) == 1) then 
            nb1 = nb(1)
            nb2 = 0
        elseif (size(nb, 1) == 2) then   
            nb1 = nb(1)
            nb2 = nb(2)
        else 
            ! This shouldn't happen
            call gdErrorHandler('More than two cells detected as face neighbours, check grid interconnection')
        end if 
        !deallocate(nb)

        ! Sanity checks
        if ( (nb1 == 0) .or. (accountforGC .and. (nb2 == 0)) ) then 
            ! Call error, this should not be the case
            print *, 'face ID: ', i
            call gdErrorHandler('ComputeGridInterconnections:' &
            // 'could not find neighbours for this face')
        end if

        ! Set BF & BV logical
        if (accountforGC) then
            if (c%GC(nb2) .or. c%GC(nb1) ) then
                ! This is a boundary face
                f%BF(i) = .true.
                v%BV(f%vert(i,:)) = .true.
            end if
        else
            if (nb2 == 0) then
                ! This is a boundary face
                f%BF(i) = .true.
                v%BV(f%vert(i,:)) = .true.
            end if
        end if
    end do


    ! Additional vertex interconnections
    !===================================
    ! Vertex neighbours
    ! Check allocation
    if (.not. allocated(v%neigP)) then
        allocate(v%neigP(v%ntot,2))
    end if

    ! Compute how much neighbours each vertex has by looping over faces
    v%neigP(:,:) = 0
    do i = 1, f%ntot
        ! Get face vertices
        tfv = f%vert(i,:)

        ! Update counters
        v%neigP(tfv(1),2) = v%neigP(tfv(1),2) + 1
        v%neigP(tfv(2),2) = v%neigP(tfv(2),2) + 1
    end do

    ! Compute the first column of neigP
    v%neigP(1,1) = 1
    do i = 2, v%ntot
        v%neigP(i,1) = v%neigP(i-1,1) + v%neigP(i-1,2)
    end do

    ! Check the neiglist allocation
    f%ncell = sum(v%neigP(:,2))
    if (allocated(v%neig)) then 
        ! This shouldn't be the case, but we hedge for it. Reallocate
        deallocate(v%neig)
        allocate(v%neig(f%ncell))
    else
        allocate(v%neig(f%ncell))
    end if

    ! Compute the neiglist
    allocate(vcount(v%ntot))
    vcount(:) = 0
    do i = 1, f%ntot
        ! Get the face vertices
        tfv = f%vert(i,:)

        ! Add the neighbours
        sp = v%neigP(tfv(1),1) + vcount(tfv(1))
        v%neig(sp) = tfv(2)
        sp = v%neigP(tfv(2),1) + vcount(tfv(2))
        v%neig(sp) = tfv(1)

        ! Update counter
        vcount(tfv) = vcount(tfv) + 1

    end do
    deallocate(vcount)

    ! Compute vertex faces
    do i = 1, f%ntot
        v%faceP(f%vert(i, :), 2) = v%faceP(f%vert(i, :), 2) + 1;
    end do
    v%nface = sum(v%faceP(:, 2))
    v%faceP(1, 1) = 1
    do i = 2, v%ntot 
        v%faceP(i, 1) = v%faceP(i-1, 1) + v%faceP(i-1, 2)
    end do 
    allocate(v%face(v%nface), vcc(v%ntot)) 
    vcc(:) = 0
    do i = 1, f%ntot
        tfv = f%vert(i, :)
        indv = vcc(tfv) + v%faceP(tfv, 1)
        v%face(indv) = i
        vcc(tfv) = vcc(tfv) + 1
    end do 


    ! Sort vertex neighbours and vertex cells
    !========================================
    
    ! Loop over all vertices
    do i = 1, v%ntot ! v%ntot

        ! Check how many distinct cell sequences there are by checking 
        ! the amount of faces of this vertex and the amount of cells. 
        ! If the difference is zero -> internal vertex (check if this is
        ! true, otherwise throw error), starting cell doesn't matter. 
        ! If one or higher -> only possible if boundary vertex 
        ! (check this). Difference indicates the amount of cell 
        ! sequences.

        ! Note: we have to account for guard cells here. Normally, if 
        ! guard cells are present, each boundary vertex should have a
        ! multiple of two amount of guard cells (more than 2 only holds
        ! for boundary vertices that connect two domains in a pointwise
        ! way (or more domains)).

        ! Number of cell sequences
        tcs = 0
        tcs = v%neigP(i, 2) - v%cellP(i, 2)

        ! Account for guard cells, if any
        if (accountforGC) then 
            tcs = tcs + count(c%GC(v%cell(v%cellP(i, 1):v%cellP(i, 2)+v%cellP(i, 1)-1)))
        end if

        ! Check
        if ( (tcs == 0) .and. (v%BV(i)) ) then 
            ! Throw error - this vertex should be an internal vertex,
            ! not a boundary one. Check the grid interconnection. 
            print *, 'vertex ID: ', i
            call gdErrorHandler(&
                'ComputeGridInterconnections: this vertex is ' &
                // 'falsely classified as boundary vertex')
        end if
        if ( (tcs > 0) .and. (.not. v%BV(i)) ) then 
            ! Throw error - this vertex should be a boundary vertex.
            ! Check the grid interconnection. 
            print *, 'vertex ID: ', i
            call gdErrorHandler(&
                'ComputeGridInterconnections: this vertex is ' &
                // 'falsely classified as internal vertex')
        end if

        ! Issue warning when tcs is higher than 1 - only occurs if two
        ! boundaries touch each other in one point, which should hardly
        ! ever be the case. 
        if (tcs > 1) then 
            print *, 'WARNING: boundary vertex ', i, &
                ' found with multiple cell sequences - check interconnection'
        end if 

        ! Set to one if equal to zero (for internal vertices)
        if (tcs == 0) then
            tcs = 1
        end if

        ! Loop over all cell sequences
        nvc = v%cellP(i,2)
        allocate(cellfound(nvc))
        allocate(localID(nvc))
        allocate(allvertcells(nvc))
        cellfound(:) = .false.
        sp = v%cellP(i,1)
        ep = v%cellP(i,1) + nvc-1
        allvertcells = v%cell(sp:ep)

        tcf = 0
        do j = 1, tcs

            ! Get the first cell and its faces
            startcellnotfound = .true.
            thiscell = 0
            if (.not. v%BV(i)) then 
                ! Internal cell, simply take the first cell, since only 
                ! one sequence
                sp = v%cellP(i,1)
                thiscell = v%cell(sp)


                ! Update 
                cellfound(1) = .true.
                startcellnotfound = .false.

            else
                ! Boundary cell, search for cell with only one (or none)
                ! faces in common with the other cells

                ! Find cell with one or zero common faces with the other
                ! cells
                
                m = 1
                do while (startcellnotfound .and. (m <= nvc))

                    ! Get the current cell
                    if (cellfound(m)) then
                        ! Update m and skip the rest
                        m = m + 1
                        cycle
                    end if

                    ! Get the faces of the current cell
                    ! tc = allremcells(m) ! current cell
                    tc = allvertcells(m)
                    fc = GetCellFace(c, tc)
                    !sp = c%faceP(tc,1)
                    !ep = c%faceP(tc,1) + c%faceP(tc,2) - 1 
                    !allocate(fc(c%faceP(tc,2)))
                    !fc = c%face(sp:ep)

                    ! Check how many common faces there are with faces
                    ! of the other remaining cells
                    ncf = 0 ! number of common faces
                    do k = 1, nvc
                        ! Skip if k == m or if cell is already found
                        if ((k == m) .or. cellfound(k)) then 
                            cycle 
                        end if

                        ! Get the faces of the next cell
                        tcn = allvertcells(k)
                        sp = c%faceP(tcn,1)
                        ep = c%faceP(tcn,1) + c%faceP(tcn,2) - 1 
                        allocate(fn(c%faceP(tcn,2)))
                        fn = c%face(sp:ep)

                        ! Count how many faces they have in common
                        do q = 1, c%faceP(tc,2)
                            ncf = ncf + count(fc(q) == fn)
                        end do

                        ! Deallocate
                        deallocate(fn)

                    end do

                    ! Deallocate
                    deallocate(fc)

                    ! Check how many common faces there are
                    !print *, ncf
                    if (ncf < 2) then 
                        ! Cell found
                        startcellnotfound = .false.

                        ! Set current cell
                        sp = v%cellP(i,1) + m - 1
                        thiscell = v%cell(sp)

                        ! Update logicals
                        ! cellfound(thislocalID(m)) = .true.
                        cellfound(m) = .true.

                        ! Set the remaining cells
                        !allocate(tcrem(nremcells-1))
                        !tcrem = allremcells((/(q, q = i, m-1),(q, q = m+1,nremcells)/))


                    else
                        ! Update m
                        m = m + 1

                    end if

                end do

                ! Check
                if (startcellnotfound) then
                    ! Should have found starting cell, throw error. 
                    ! Something is wrong then with the interconnection,
                    ! and it could be that this boundary vertex should 
                    ! in fact not be a boundary vertex. 
                    print *, 'vertex ID: ', i 
                    call gdErrorHandler('ComputeGridInterconnections: could not find starting cell for this boundary vertex')

                end if 


            end if

            ! Find the sequence
            !==================
            ! Now that the first cell of the sequence has been found, we
            ! need to extract the sequence itself. 

            ! Loop over the remaining cells
            k = 1
            vc = 0 ! vertex counter
            allnotfound = .true. 
            
            do while (k <= nvc)

                ! Set the faces of the current cell
                allocate(fc(c%faceP(thiscell,2)))
                sp = c%faceP(thiscell, 1)
                ep = sp + c%faceP(thiscell,2)-1
                fc = c%face(sp:ep)

                ! Take the next cell
                if (cellfound(k)) then
                    ! Update k, cycle
                    k = k + 1
                    deallocate(fc)
                    cycle
                else
                    nextcell = allvertcells(k)
                end if
                
                ! Get cell faces
                allocate(fn(c%faceP(nextcell,2)))
                sp = c%faceP(nextcell, 1)
                ep = sp + c%faceP(nextcell,2)-1
                fn = c%face(sp:ep)

                ! Check for common faces - should only be one
                ncf = 0
                do q = 1, c%faceP(thiscell, 2)
                    ncf = ncf + count(fc(q) == fn)
                    if ((count(fc(q) == fn)) == 1) then
                        ! Set the current face
                        tcf = fc(q)
                    end if
                end do

                ! Sanity check
                if (ncf > 1) then
                    ! Throw error: the next cell is connected with too
                    ! many faces to the previous one, which shouldn't be
                    ! possible in any decent grid. Though in some 
                    ! exceptional cases, this may be allowed, we don't 
                    ! support that here. 
                    call gdErrorHandler('ComputeGridInterconnections: '&
                        // 'too many common faces detected when ' &
                        // 'sorting vertex neighbours, check grid ' &
                        // 'interconnection')
                end if

                ! Check if a common face was found
                if (ncf == 1) then 
                    ! Add the first vertex if this is the first found 
                    ! cell
                    if (vc == 0) then
                        ! Get other face with this vertex

                        ! Get the faces that have the current vertex
                        if (c%GC(thiscell)) then
                            ! Hedge for guard cells: only one face
                            allocate(allfv(1))
                            allocate(allfvind(1))
                            allocate(tcf2(1))
                            tcf2(1) = fc(1)
                        else
                            allocate(allfvind(c%faceP(thiscell,2)))
                            allfvind = (f%vert(fc, 1) == i) .or. &
                                (f%vert(fc, 2) == i)
                            allocate(allfv(count(allfvind)))
                            allfv = pack(fc, allfvind)
                            allocate(tcf2(count(allfvind)))
                            tcf2 = pack(allfv, allfv .ne. tcf)
                        end if
                        
                        ! Check 
                        ntcf2 = size(tcf2)
                        if (ntcf2 .ne. 1) then
                            call gdErrorHandler( &
                            'ComputeGridInterconnections: could not ' &
                            // 'find a single common face, check grid' &
                            //' interconnections')
                        end if

                        ! Add cell and vertex neighbour - neighbour only
                        ! if the current cell is not a guard cell
                        
                        sp = v%cellP(i,1)
                        v%cell(sp) = thiscell 
                        if (.not. c%GC(thiscell)) then
                            sp = v%neigP(i,1)
                            tfv = f%vert(tcf2(1),:)
                            if (tfv(1) == i) then
                                v%neig(sp) = tfv(2)
                            else
                                v%neig(sp) = tfv(1)
                            end if

                            ! Update vc
                            vc = vc+1
                        end if

                        ! Housekeeping
                        deallocate(allfvind)
                        deallocate(allfv)
                        deallocate(tcf2)

                    end if 

                    ! Add the cell and vertex neighbour
                    sp = v%cellP(i,1) + vc
                    v%cell(sp) = nextcell
                    sp = v%neigP(i,1) + vc
                    tfv = f%vert(tcf,:)
                    if (tfv(1) == i) then 
                        v%neig(sp) = tfv(2)
                    else
                        v%neig(sp) = tfv(1)
                    end if

                    ! Update counter
                    vc = vc + 1

                    ! Update logicals
                    cellfound(k) = .true.

                    ! Reset iterator
                    k = 1

                    ! Set current cell
                    thiscell = nextcell

                else
                    ! Update the counter
                    k = k + 1
                end if

                ! Housekeeping
                deallocate(fc, fn)

            end do

            ! Add last vertex
            !================
            ! Only in case of boundary vertices
            if ( (v%BV(i)) .and. (v%cellP(i,2) > 1) ) then

                ! First, update vc - is now one-off
                ! vc = vc - 1

                ! Sanity check
                if (vc > v%neigP(i,2)) then
                    ! Throw error. This can happen if the supposed 
                    ! boundary vertex is no actual boundary vertex.
                    call gdErrorHandler('ComputeGridInterconnections: ' &
                        // 'supposed boundary vertex may be internal ' &
                        // 'vertex. Check grid interconnection')
                end if

                ! Set the faces of the current cell
                allocate(fc(c%faceP(thiscell,2)))
                sp = c%faceP(thiscell, 1)
                ep = sp + c%faceP(thiscell,2)-1
                fc = c%face(sp:ep)

                ! Add the vertex neighbour
                if (c%GC(thiscell)) then
                    ! Hedge for guard cells: only one face
                    allocate(tcf2(1))
                    allocate(allfv(1))
                    allocate(allfvind(1))
                    tcf2(1) = fc(1)
                else
                    allocate(allfvind(c%faceP(thiscell,2)))
                    allfvind = (f%vert(fc, 1) == i) .or. &
                        (f%vert(fc, 2) == i)
                    allocate(allfv(count(allfvind)))
                    allfv = pack(fc, allfvind)
                    tcf2 = pack(allfv, allfv .ne. tcf)

                end if
                
                ! Check 
                ntcf2 = size(tcf2)
                if (ntcf2 .ne. 1) then
                    call gdErrorHandler( &
                    'ComputeGridInterconnections: could not ' &
                    // 'find a single common face, check grid' &
                    //' interconnections')
                end if

                ! Add cell and vertex neighbour
                sp = v%neigP(i,1) + vc
                tfv = f%vert(tcf2(1),:)
                if (tfv(1) == i) then
                    v%neig(sp) = tfv(2)
                else
                    v%neig(sp) = tfv(1)
                end if

                ! Housekeeping
                deallocate(allfvind)
                deallocate(allfv)
                deallocate(fc, tcf2)
                
            end if
        end do

        ! Sanity check
        if (any(.not. cellfound)) then
            ! Not all cells were ordened, throw error
            print *, 'vertex ID: ', i
            call gdErrorHandler('ComputeGridInterconnections: could not sort all cells, check grid interconnection')
        end if 

        ! Housekeeping
        deallocate(cellfound)
        deallocate(localID)
        deallocate(allvertcells)

    end do

    !do i = 1, v%ntot 
    !    print *, i, v%BV(i), v%neig(v%neigP(i,1):v%neigP(i,1)+v%neigP(i,2)-1)
    !end do
    !    do i = 1, v%ntot 
    !    print *, i, v%BV(i), v%cell(v%cellP(i,1):v%cellP(i,1)+v%cellP(i,2)-1)
    !end do
    

    ! Add to grid
    !============
    grid%cell   = c 
    grid%face   = f 
    grid%vert    = v


end subroutine