subroutine ExtractGridData(grid, meth, gridoptions)

    ! Description
    !============
    ! Process the grid data given in 'grid' by adding the necessary 
    ! additional datastructures, checking for ghost vertices, ... 

    ! Notes
    !======
    ! Note 1: it is assumed that the grid that is read does not contain
    ! any guard cells yet. Otherwise this routine has to be extended. 

    ! Note 2: it is checked whether 'ghost' vertices are present, which 
    ! are vertices that do not belong to any face. These tend to occur
    ! in files that were constructed based on a structured grid where
    ! typically the x-point was duplicated. These ghost vertices will be 
    ! removed and the grid structure will be updated accordingly. 

    ! Note 3: it is assumed that for the b2fgmtry case the carre grid 
    ! generator is used (or that at least the face labels are set up 
    ! with the same convention as used in carre)

    ! Initialize
    !===========
    ! Declare modules
    use gdmod_types 
    use gdmod_interfaces
    use gdmod_plots
    use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    type(GridUDT), intent(inout)    :: grid
    type(GridOptionsUDT)            :: gridoptions
    character(len=*), intent(in)    :: meth

    ! Loop variables
    integer(I8)                 :: i, j, k, iFT, ib, il 

    ! Auxiliary variables (
    type(FluxDataUDT)           :: fluxdata
    type(VertexUDT)             :: newverts ! necessary if ghost vertices are present

    integer(I8)                 :: itf, ntf, ngv, nbnd, nfpb, nseg, &
        nlabels 
    integer(I8), allocatable    :: tf(:), tfv(:,:), indgv(:), vdiff(:)
    integer(I8), allocatable    :: gglabels(:), gdlabels(:), &
        sortindex(:), temparray(:,:), tempfaces(:), segstart(:)

    logical, allocatable        :: isghostvert(:), mask(:), &
                                ispolygonstart(:)

    integer(I8), allocatable    :: facevec(:) ! simply 1:grid%faces%ntot
    integer(I8)                 :: start, nvert, tv(1:4), newv(1:4), &
        segend

    ! Plotting
    real(R8)                    :: NaN
    integer(I8)                 :: ntemp
    real(R8), allocatable       :: tempx(:), tempy(:) ! coordinates
    logical                     :: makedebugplots = .false. 

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
            grid%vert = newverts

            ! Deallocate
            deallocate(indgv)
            deallocate(mask)
        endif 

        ! Reorden cell vertices
        !======================
        ! In a classical quadrilateral grid from CARRE(2), the vertices 
        ! are not ordened (counter)clockwise. Here, this should be the
        ! case. Therefore, for each cell that has 4 vertices (normally
        ! every cell here, otherwise we print a warning), the 3nd and 
        ! 4th vertex are switched. Guard cells without guard vertices 
        ! are skipped (these should, in fact, be removed) 
        if (any( (grid%cells%vertP(:, 2) .ne. 4) .and. (grid%cells%vertP(:,2) .ne. 2) )) then 
            print *, 'non-classical grid, vertex reordering may be wrong'
        end if 
        do i = 1, grid%cells%ntot
             if ((grid%cells%vertP(i,2) .ne. 4) .and. (grid%cells%vertP(i,2) .ne. 2))  then
                ! Do nothing
            else 
                start = grid%cells%vertP(i,1)
                nvert = grid%cells%vertP(i,2)
                tv = grid%cells%vertlist(start:start+nvert)
                newv = tv 
                newv(3) = tv(4)
                newv(4) = tv(3)
                grid%cells%vertlist(start:start+nvert) = newv
            end if
        end do

        ! Extract boundaries
        !===================
        ! Get the supported mapping between boundary labels 
        gglabels = gridoptions%facelabelmappingGG
        gdlabels = gridoptions%facelabelmappingGD 

        ! Substitute labels
        do i = 1, size(gridoptions%facelabelsubfrom)
            where (grid%data%regions%facelabel == gridoptions%facelabelsubfrom(i)) &
                grid%data%regions%facelabel = gridoptions%facelabelsubto(i)
        end do

        ! Loop over all face labels (not regions here!) to precompute
        ! number of grid boundaries (can be more/less)
        allocate(mask(grid%faces%ntot))
        allocate(facevec(grid%faces%ntot))
        facevec(:) = (/(i, i=1,grid%faces%ntot,1)/)
        nlabels = size(gglabels) 
        nbnd = 0 ! number of boundaries
        do il = 1, nlabels 
            ! Get the faces of this boundary
            mask(:) = grid%data%regions%facelabel == gglabels(il);
            nfpb = count(mask)

            ! Check
            if (nfpb == 0) then 
                ! this will not become a boundary, skip rest of loop
                cycle
            end if

            ! Extract faces
            allocate(tempfaces(nfpb))
            tempfaces = pack(facevec, mask) 

            ! Determine number of boundaries by sorting
            allocate(sortindex(nfpb))
            allocate(ispolygonstart(nfpb))
            allocate(temparray(nfpb, 2))
            temparray(:, :) = grid%faces%vert(tempfaces, :)
            call SortPolygonEdges(temparray, nfpb, sortindex, ispolygonstart)
            nbnd = nbnd + count(ispolygonstart)

            ! Housekeeping
            deallocate(tempfaces, sortindex, ispolygonstart, temparray)

        end do 

        ! Extract boundaries
        allocate(grid%bnd(nbnd))
        ib = 0 ! boundary counter
        do il = 1, nlabels
            

            ! Get the faces of this boundary
            mask(:) = grid%data%regions%facelabel == gglabels(il);
            nfpb = count(mask)

            ! Check
            if (nfpb == 0) then 
                ! Don't add as a boundary, skip rest of loop
                cycle 
            end if

            ! Sort the boundary faces
            allocate(sortindex(nfpb), ispolygonstart(nfpb), &
                tempfaces(nfpb), temparray(nfpb,2))

            tempfaces = pack(facevec, mask)
            temparray(:,:) = grid%faces%vert(tempfaces, :)
            call SortPolygonEdges(temparray, nfpb, sortindex, ispolygonstart)
            tempfaces = tempfaces(sortindex)
            
            ! Loop over all found boundary segments
            nseg = count(ispolygonstart)
            allocate(segstart(nseg))
            segstart = pack([(k, k = 1, nfpb)], ispolygonstart)
            do j = 1, nseg
                ! Update the boundary counter
                ib = ib + 1
                
                ! Compute end segment index
                if (j < nseg) then 
                    segend = segstart(j+1)-1
                else
                    segend = nfpb
                end if 

                ! Add the boundary ID
                grid%bnd(ib)%ID = gdlabels(il)
                grid%bnd(ib)%nfaces = segend-segstart(j)+1

                ! Allocate this boundary
                call AllocateBnd(grid%bnd(ib))

                ! Add faces (already sorted before)
                grid%bnd(ib)%faces(:) = tempfaces(segstart(j):segend)

                ! Extract vertices
                call ExtractPolygonVertices( & 
                    grid%faces%vert(grid%bnd(ib)%faces,:), &
                    grid%bnd(ib)%nfaces, grid%bnd(ib)%vert)
            end do 

            ! Deallocate
            deallocate(sortindex, ispolygonstart, temparray, tempfaces, &
                segstart)
        end do

        ! Make a plot to check
        if (makedebugplots) then
            ntemp = 0
            do ib = 1, nbnd
                ntemp = ntemp + grid%bnd(ib)%nvert + 1
            end do 
            allocate(tempx(ntemp))
            allocate(tempy(ntemp))
            NaN = IEEE_VALUE(nan, IEEE_QUIET_NAN)
            ntemp = 1
            do ib = 1, nbnd
                !print *, grid%bnd(ib)%vert
                tempx(ntemp:ntemp+grid%bnd(ib)%nvert-1) = &
                   grid%vert%x(grid%bnd(ib)%vert)
                tempy(ntemp:ntemp+grid%bnd(ib)%nvert-1) = &
                    grid%vert%y(grid%bnd(ib)%vert)
                tempx(ntemp+grid%bnd(ib)%nvert) = NaN
                tempy(ntemp+grid%bnd(ib)%nvert) = NaN
                ntemp = ntemp + grid%bnd(ib)%nvert + 1
            end do

            ntemp = size(tempx)
            call Plot2DPolygon(tempx, tempy, ntemp, '-p')
            deallocate(tempx)
            deallocate(tempy)
        end if


        ! Deallocate
        deallocate(isghostvert)
        deallocate(mask)

    case ('traduit')

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
            grid%vert = newverts

            ! Deallocate
            deallocate(indgv)
            deallocate(mask)
        endif 

        ! Extract boundaries
        !===================
        ! Get the supported mapping between boundary labels 
        gglabels = gridoptions%facelabelmappingGG
        gdlabels = gridoptions%facelabelmappingGD 

        ! Substitute labels
        do i = 1, size(gridoptions%facelabelsubfrom)
            where (grid%data%regions%facelabel == gridoptions%facelabelsubfrom(i)) &
                grid%data%regions%facelabel = gridoptions%facelabelsubto(i)
        end do

        ! Loop over all face labels (not regions here!) to precompute
        ! number of grid boundaries (can be more/less)
        allocate(mask(grid%faces%ntot))
        allocate(facevec(grid%faces%ntot))
        facevec(:) = (/(i, i=1,grid%faces%ntot,1)/)
        nlabels = size(gglabels) 
        nbnd = 0 ! number of boundaries
        do il = 1, nlabels 
            ! Get the faces of this boundary
            mask(:) = grid%data%regions%facelabel == gglabels(il);
            nfpb = count(mask)

            ! Check
            if (nfpb == 0) then 
                ! this will not become a boundary, skip rest of loop
                cycle
            end if

            ! Extract faces
            allocate(tempfaces(nfpb))
            tempfaces = pack(facevec, mask) 

            ! Determine number of boundaries by sorting
            allocate(sortindex(nfpb))
            allocate(ispolygonstart(nfpb))
            allocate(temparray(nfpb, 2))
            temparray(:, :) = grid%faces%vert(tempfaces, :)
            call SortPolygonEdges(temparray, nfpb, sortindex, ispolygonstart)
            nbnd = nbnd + count(ispolygonstart)

            ! Housekeeping
            deallocate(tempfaces, sortindex, ispolygonstart, temparray)

        end do 

        ! Extract boundaries
        allocate(grid%bnd(nbnd))
        ib = 0 ! boundary counter
        do il = 1, nlabels
            

            ! Get the faces of this boundary
            mask(:) = grid%data%regions%facelabel == gglabels(il);
            nfpb = count(mask)

            ! Check
            if (nfpb == 0) then 
                ! Don't add as a boundary, skip rest of loop
                cycle 
            end if

            ! Sort the boundary faces
            allocate(sortindex(nfpb), ispolygonstart(nfpb), &
                tempfaces(nfpb), temparray(nfpb,2))

            tempfaces = pack(facevec, mask)
            temparray(:,:) = grid%faces%vert(tempfaces, :)
            call SortPolygonEdges(temparray, nfpb, sortindex, ispolygonstart)
            tempfaces = tempfaces(sortindex)
            
            ! Loop over all found boundary segments
            nseg = count(ispolygonstart)
            allocate(segstart(nseg))
            segstart = pack([(k, k = 1, nfpb)], ispolygonstart)
            do j = 1, nseg
                ! Update the boundary counter
                ib = ib + 1
                
                ! Compute end segment index
                if (j < nseg) then 
                    segend = segstart(j+1)-1
                else
                    segend = nfpb
                end if 

                ! Add the boundary ID
                grid%bnd(ib)%ID = gdlabels(il)
                grid%bnd(ib)%nfaces = segend-segstart(j)+1

                ! Allocate this boundary
                call AllocateBnd(grid%bnd(ib))

                ! Add faces (already sorted before)
                grid%bnd(ib)%faces(:) = tempfaces(segstart(j):segend)

                ! Extract vertices
                call ExtractPolygonVertices( & 
                    grid%faces%vert(grid%bnd(ib)%faces,:), &
                    grid%bnd(ib)%nfaces, grid%bnd(ib)%vert)
            end do 

            ! Deallocate
            deallocate(sortindex, ispolygonstart, temparray, tempfaces, &
                segstart)
        end do

        ! Make a plot to check
        if (makedebugplots) then
            ntemp = 0
            do ib = 1, nbnd
                ntemp = ntemp + grid%bnd(ib)%nvert + 1
            end do 
            allocate(tempx(ntemp))
            allocate(tempy(ntemp))
            NaN = IEEE_VALUE(nan, IEEE_QUIET_NAN)
            ntemp = 1
            do ib = 1, nbnd
                !print *, grid%bnd(ib)%vert
                tempx(ntemp:ntemp+grid%bnd(ib)%nvert-1) = &
                   grid%vert%x(grid%bnd(ib)%vert)
                tempy(ntemp:ntemp+grid%bnd(ib)%nvert-1) = &
                    grid%vert%y(grid%bnd(ib)%vert)
                tempx(ntemp+grid%bnd(ib)%nvert) = NaN
                tempy(ntemp+grid%bnd(ib)%nvert) = NaN
                ntemp = ntemp + grid%bnd(ib)%nvert + 1
            end do

            ntemp = size(tempx)
            call Plot2DPolygon(tempx, tempy, ntemp, '-p')
            deallocate(tempx)
            deallocate(tempy)
        end if


        ! Deallocate
        deallocate(isghostvert)
        deallocate(mask)

    case default

        ! Unknown extraction case, throw error
        call gdErrorHandler('ExtractGridData: unknown method')

    end select 

end subroutine