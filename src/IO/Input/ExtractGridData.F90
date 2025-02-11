subroutine ExtractGridData(grid, meth, gridoptions)

    ! Description
    !============
    ! Process the grid data given in 'grid' by adding the necessary 
    ! additional datastructures, checking for ghost vertices, ... 

    ! Notes
    !======
    ! Note 1: it is assumed that the grid that is read does not contain
    ! any guard cells yet. Otherwise this routine has to be extended. 
    ! Additionally, it is assumed that no 'ghost' vertices (i.e. vertices 
    ! without any other connection to the grid) are present. 

    ! Initialize
    !===========
    ! Declare modules
    use goatmod_types 
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

    integer(I8)                 :: itf, ntf, nbnd, nfpb, nseg, &
        nlabels 
    integer(I8), allocatable    :: tf(:), tfv(:,:)
    integer(I8), allocatable    :: &
        sortindex(:), temparray(:,:), tempfaces(:), segstart(:)

    logical, allocatable        :: mask(:), ispolygonstart(:), &
        isbranchingpolygon(:)

    integer(I8), allocatable    :: facevec(:) ! simply 1:grid%face%ntot
    integer(I8)                 :: segend

    ! Plotting
    real(R8)                    :: NaN
    integer(I8)                 :: ntemp
    real(R8), allocatable       :: tempx(:), tempy(:) ! coordinates
    logical                     :: makedebugplots = .false. 

    ! Initialize
    !===========
    ! Associate
    associate(gglabels => gridoptions%facelabelmappingGG, &
        gdlabels => gridoptions%facelabelmappingGD)
    
    ! Initialize
    grid%vert%fieldlineID = 0
    grid%data%fluxdata%fluxsurfaceID = 0

    ! Check method
    !=============
    select case(meth)

    case ('traduitb2us', 'b2fgmtry_us')

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
            tfv(1:ntf,:) = grid%face%vert(tf(1:ntf),:)

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

        ! Extract boundaries
        !===================
        ! Get the supported mapping between boundary labels 
        gglabels = gridoptions%facelabelmappingGG
        gdlabels = gridoptions%facelabelmappingGD 

        ! Substitute labels
        do i = 1, size(gridoptions%facelabelsubfrom)
            where (grid%face%label == gridoptions%facelabelsubfrom(i)) &
                grid%face%label = gridoptions%facelabelsubto(i)
        end do

        ! Loop over all face labels (not regions here!) to precompute
        ! number of grid boundaries (can be more/less)
        allocate(mask(grid%face%ntot))
        allocate(facevec(grid%face%ntot))
        facevec(:) = (/(i, i=1,grid%face%ntot,1)/)
        nlabels = size(gglabels) 
        nbnd = 0 ! number of boundaries
        do il = 1, nlabels 
            ! Get the faces of this boundary
            mask(:) = grid%face%label == gglabels(il);
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
            allocate(ispolygonstart(nfpb), isbranchingpolygon(nfpb))
            allocate(temparray(nfpb, 2))
            temparray(:, :) = grid%face%vert(tempfaces, :)
            call SortPolygonEdges(temparray, nfpb, sortindex, ispolygonstart, &
                isbranchingpolygon)
            nbnd = nbnd + count(ispolygonstart)

            ! Housekeeping
            deallocate(tempfaces, sortindex, ispolygonstart, temparray, isbranchingpolygon)

        end do 

        ! Extract boundaries
        allocate(grid%bnd(nbnd))
        ib = 0 ! boundary counter
        do il = 1, nlabels
            

            ! Get the faces of this boundary
            mask(:) = grid%face%label == gglabels(il);
            nfpb = count(mask)

            ! Check
            if (nfpb == 0) then 
                ! Don't add as a boundary, skip rest of loop
                cycle 
            end if

            ! Sort the boundary faces
            allocate(sortindex(nfpb), ispolygonstart(nfpb), &
                tempfaces(nfpb), temparray(nfpb,2), isbranchingpolygon(nfpb))

            tempfaces = pack(facevec, mask)
            temparray(:,:) = grid%face%vert(tempfaces, :)
            call SortPolygonEdges(temparray, nfpb, sortindex, ispolygonstart, &
                isbranchingpolygon)
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
                grid%bnd(ib)%nface = segend-segstart(j)+1

                ! Allocate this boundary
                call AllocateBnd(grid%bnd(ib))

                ! Add faces (already sorted before)
                grid%bnd(ib)%face(:) = tempfaces(segstart(j):segend)

                ! Extract vertices
                call ExtractPolygonVertices( & 
                    grid%face%vert(grid%bnd(ib)%face,:), &
                    grid%bnd(ib)%nface, grid%bnd(ib)%vert)
            end do 

            ! Deallocate
            deallocate(sortindex, ispolygonstart, temparray, tempfaces, &
                segstart, isbranchingpolygon)
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


    case default

        ! Unknown extraction case, throw error
        call gdErrorHandler('ExtractGridData: unknown method')

    end select 

    end associate

end subroutine