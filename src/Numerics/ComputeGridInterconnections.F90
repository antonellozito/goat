subroutine ComputeGridInterconnections(grid)

    ! Description
    !============
    ! Compute additional grid topology data starting from basic grid
    ! information. The following fields are required (see gdmod_types
    ! for an explanation of the different fields)
    !
    ! Vertices: ntot, BV
    ! Faces: ntot, vert
    ! Cells: ntot, vertP, vertlist, nvertlist 
    !   
    ! All other fields are recomputed based on this basic 
    ! interconnection data.             

    ! Notes
    !======
    ! Note 1: it is assumed that the grid that is read does not contain
    ! any guard cells yet. Otherwise this routine has to be extended. 

    ! Note 2: it is assumed that all fields of which the size is 
    ! predetermined by the number of vertices/faces/cells are already
    ! allocated. Other fields are allocated or overwritten. 

    ! Note 3: the original cell vertices are replaced by a sorted list.
    ! This list is generated based on the faces of the cell and their
    ! corresponding vertices. It is assumed that each cell has the same
    ! amount of faces as it has vertices.

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
    type(GridUDT), intent(inout)    :: grid

    ! Loop variables
    integer(I8)                 :: i, j, k

    ! Auxiliary variables (
    type(VertexUDT)             :: v 
    type(FaceUDT)               :: f
    type(CellUDT)               :: c

    integer(I8)                 :: nc, nv, nf
    integer(I8), allocatable    :: fcount(:), vcount(:), tv(:)
    integer(I8)                 :: ntv, tf, ind

    ! Unpack & initialize
    !==================
    ! Data structures 
    v = grid%vert
    f = grid%faces
    c = grid%cells 

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
    if (allocated(c%facelist)) then
        ! Print out warning
        print *, 'ComputeGridInterconnections: recomputing cells%facelist'

        ! Deallocate
        deallocate(c%facelist)
    end if
    if (allocated(v%celllist)) then 
        ! print out warning
        print *, 'ComputeGridInterconnections: recomputing vert%celllist'

        ! Deallocate
        deallocate(v%celllist)
    end if

    ! Make the pointer list for cell faces
    c%faceP(1,1) = 1
    c%faceP(1,2) = c%vertP(1,2) ! same amount of faces as vertices
    do i = 2, nc
        ! Pointer
        c%faceP(i,1) = c%faceP(i-1,1) + c%faceP(i-1,2)

        ! Number of faces
        c%faceP(i,2) = c%vertP(i,2)

    end do

    ! Compute the total faces in cells%facelist
    c%nfacelist = sum(c%faceP(:,2)) 

    ! Allocate
    allocate(c%facelist(c%nfacelist))

    ! Make the pointer list for vertex cell neighbours
    v%cellP(:,:) = 0
    do i = 1, nc
        ! Get the number of vertices of the current cell
        ntv = c%vertP(i,2)

        ! Allocate
        allocate(tv(ntv))

        ! Get the vertices of the current cell
        tv = c%vertlist(c%vertP(i,1):c%vertP(i,1)+ntv-1)

        ! Update the counter
        v%cellP(tv,:) = v%cellP(tv,:)+1

        ! Deallocate
        deallocate(tv)

    end do

    ! Compute the total number of cells in vert%celllist
    v%ncelllist = sum(v%cellP(:,2))

    ! Allocate
    allocate(v%celllist(v%ncelllist))

    ! Construct the pointer
    v%cellP(1,1) = 1
    do i = 2, nv 
        v%cellP(i,1) = v%cellP(i-1,1) + v%cellP(i-1,2)
    end do

    ! Set vertex celllist and cell facelist and face neighbours
    fcount(:) = 1
    vcount(:) = 1
    do i = 1, nc ! loop over all cells
        ! Get the number of vertices of this cell
        ntv = c%vertP(i,2)

        ! Allocate
        allocate(tv(ntv))

        ! Get the vertices of the current cell
        tv(:) = c%vertlist(c%vertP(i,1):(c%vertP(i,1)+ntv))

        ! Loop over ntv-1 vertices
        do j = 1, ntv-1

            ! Get all faces of the j'th vertex
            do k = j+1, ntv
                ! Get the next face (if it exists)
                call MapVertexPairToFace(tv(j), tv(k), f%vert, nf, tf)

                ! Check
                if (tf < 0) then
                   ! do nothing, continue with next vertex
                else

                    ! Add the current cell index as neighbour to this face
                    f%neig(tf,fcount(tf)) = i

                    ! Update fcount
                    fcount(tf) = fcount(tf)+1

                    ! Add the current cell to the j'th vertex
                    ind = v%cellP(tv(j),1) + vcount(tv(j)) - 1
                    if (ind > v%ncelllist) then
                        call gdErrorHandler('unknown error')
                    end if
                    v%celllist(ind) = i

                    ! Update vcount
                    vcount(tv(j)) = vcount(tv(j))+1

                    ! Add cell face
                    ind = c%faceP(i,1) + j - 1
                    if (ind > c%nfacelist) then
                        call gdErrorHandler('unknown error')
                    end if
                    c%facelist(ind) = tf
                end if
            end do

        end do

        ! Deallocate
        deallocate(tv)

    end do

    ! Add to grid
    !============
    grid%cells   = c 
    grid%faces   = f 
    grid%vert    = v


end subroutine