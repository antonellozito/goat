subroutine ExtractGridData(grid, meth)

    ! Description
    !============
    ! Process the grid data given in 'grid' by adding the necessary 
    ! additional datastructures, checking for ghost vertices, ... 

    ! Notes:
    ! Note 1: it is assumed that the grid that is read does not contain
    ! any guard cells yet. Otherwise this routine has to be extended. 

    ! Note 2: it is checked whether 'ghost' vertices are present, which 
    ! are vertices that do not belong to any face. These tend to occur
    ! in files that were constructed based on a structured grid where
    ! typically the x-point was duplicated. These ghost vertices will be 
    ! removed and the grid structure will be updated accordingly. 

    ! Initialize
    !===========
    ! Declare modules
    use gdmod_types 

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    type(GridUDT)               :: grid
    character(len=*)            :: meth

    ! Loop variables
    integer(I8)                 :: i, j, iFT 

    ! Auxiliary variables (
    type(FluxDataUDT)           :: fluxdata
    type(VertexUDT)             :: newverts ! necessary if ghost vertices are present

    integer(I8)                 :: itf, ntf, ngv
    integer(I8), allocatable    :: tf(:), tfv(:,:), indgv(:), vdiff(:)

    logical, allocatable        :: isghostvert(:), mask(:)

    ! Check ghost vertices
    !=====================
    ! Treatment done later on, depending on 'meth'

    ! Allocate & initialize
    allocate(isghostvert(grid%vert%ntot))
    isghostvert(:) = .true.

    ! Loop
    do j = 1, 2
        do i = 1, grid%faces%ntot
            isghostvert(grid%faces%vert(i,j)) = .false.
        enddo
    enddo

    ! Count number of ghost vertices
    ngv = count(isghostvert .eqv. .true.) ! total number of ghost vertices

    ! Check method
    !=============
    select case(meth)

    case ('b2fgmtry')

        ! Process the grid
        !=================
        ! Extract the fieldline ID's of the vertices based on the flux 
        ! surface data

        ! Allocate
        fluxdata = grid%data%fluxdata
        allocate(tf(maxval(fluxdata%fluxsurfacefacesP(:,2),1)))
        allocate(tfv(maxval(fluxdata%fluxsurfacefacesP(:,2),1),2))

        ! Loop
        do iFT = 1, fluxdata%nFs
           ! Unpack
            itf = fluxdata%fluxsurfacefacesP(iFT,1); ! start index
            ntf = fluxdata%fluxsurfacefacesP(iFT,2); ! number of faces 
        
            ! Extract faces
            tf(1:ntf) = fluxdata%fluxsurfacefaces(itf:itf+ntf-1)

            ! Extract vertices of these faces
            tfv(1:ntf,:) = grid%faces%vert(tf(1:ntf),:)

            ! Set the flux tube index
            do j = 1, 2
                do i = 1, ntf 
                    grid%data%fluxdata%fluxsurfaceID(tfv(i,j)) = iFT
                    grid%vert%fieldlineID(tfv(i,j)) = iFT
                enddo
            enddo
        enddo

        ! Deallocate
        deallocate(tf)
        deallocate(tfv)

        ! Hedge for ghost vertices
        !=========================
        ! Remove the ghost vertex and update grid structure(s)
        if (ngv > 0) then
            ! Print out the amount of detected ghost vertices
            print *, 'detected ', ngv, ' ghost vertices, removing...'

            ! Allocate 
            allocate(indgv(ngv))
            allocate(vdiff(grid%vert%ntot))
            allocate(mask(grid%vert%ntot))
            newverts%ntot = grid%vert%ntot - ngv
            call AllocateVertices(newverts)

            ! Initialize
            vdiff(:) = 0

            ! Find the ghost vertex indices
            j = 1
            do i = 1, grid%vert%ntot
                if (isghostvert(i)) then
                    indgv(j) = i
                    j = j+1
                endif
            enddo

            ! Copy the coordinates
            mask(:)             = .true.
            mask(indgv)         = .false.
            newverts%x           = pack(grid%vert%x,mask)
            newverts%y           = pack(grid%vert%y,mask)
            newverts%fieldlineID = pack(grid%vert%fieldlineID,mask)

            ! Adjust the cell and face vertex IDs
            do  i = 1, ngv
                ! Check if any faces and cells actually contain this 
                ! vertex, if so -> throw error
                if (any(grid%faces%vert == indgv(i))) then
                    call gdErrorHandler('Ghost vertex detected yet present in face vertices')
                end if
                if (any(grid%faces%vert == indgv(i))) then
                    call gdErrorHandler('Ghost vertex detected yet present in cell vertices')
                end if

                ! Adjust mapping
                vdiff(indgv(i):grid%vert%ntot) = &
                    vdiff(indgv(i):grid%vert%ntot) + 1
                
            end do

            ! Adjust cell and faces
            do j = 1, 2
                grid%faces%vert(:,j) = grid%faces%vert(:,j) &
                    - vdiff(grid%faces%vert(:,j))
            end do

            grid%cells%vertlist = grid%cells%vertlist & 
                - vdiff(grid%cells%vertlist)

            ! Add the vertex structure
            call DeallocateVertices(grid%vert) ! first, deallocate
            grid%vert%ntot = newverts%ntot
            call AllocateVertices(grid%vert) ! allocate again
            grid%vert%x = newverts%x ! add the new structures
            grid%vert%y = newverts%y 
            grid%vert%fieldlineID = newverts%fieldlineID 

            ! Deallocate
            deallocate(indgv)
            deallocate(mask)

        endif

        ! Deallocate
        deallocate(isghostvert)

        

    case default

    end select 

end subroutine