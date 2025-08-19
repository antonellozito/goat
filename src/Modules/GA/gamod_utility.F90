!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains utility routines for the adapation modules


module gamod_utility

    ! Initialize
    !===========
    ! Load modules
    use mod_sort
    use mod_precision
    use goatmod_types
    use gdmod_utility_optimization


    ! The usual
    implicit none
    save
    public
    
    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================! 
    
    type :: QualityMetric

        ! Type that included several metrics
        real(R8), allocatable :: fcBias(:)
        real(R8), allocatable :: fcqalfc(:)
        real(R8), allocatable :: fcS(:)
        real(R8), allocatable :: cvS(:)
        real(R8), allocatable :: cvAR(:)
        real(R8), allocatable :: h_pol(:)
        real(R8), allocatable :: h_rad(:)
        real(R8), allocatable :: h_rad_psi(:)
        integer(I8) :: nCv

    contains
    
        procedure :: ComputeQM

    end type


    contains 

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!
    subroutine ComputeQM(qm,grid,options)

        ! Description
        !============
        ! Compute all quality metrics
        ! (Mirror of CalculateCvMetric.m)

        ! Declare variables
        !==================
        ! Arguments
        class(QualityMetric), intent(inout)     :: qm
        type(GridUDT), intent(in)               :: grid
        type(GAoptionsUDT), intent(in)          :: options

        ! Auxiliary
        real(R8) :: vec_n(grid%face%ntot,2), fcH(grid%face%ntot,2), &
         ncpf(grid%face%ntot), fccv(grid%face%ntot,2), &
         fcxx(grid%face%ntot)

        
    end subroutine

    subroutine CalculateQualityMetrics(grid,options,qm)
        ! Description
        !============
        ! Compute metric and criteria of cells necessary to execute grid adaptation.

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(in)           :: grid
        type(GAoptionsUDT), intent(in)      :: options
        type(QualityMetric), intent(inout)  :: qm

        ! Calculate cv metric
        call qm%ComputeQM(grid,options)


        ! Selecting splitting cell

        ! Selecting merging face
        


    end subroutine



    subroutine CheckVertOrder(grid,is_ordered,cells)

        ! Description
        !============
        ! Check whether the vertices and faces are ordered correctly,
        ! following the right turning scheme  
        ! is_ordered: indicates whether the vertices and faces of a cell are ordered
        ! cells: indicate which cells need to be checked

        ! Declare variables
        !==================
        type(GridUDT), intent(in)   :: grid
        logical, intent(inout)      :: is_ordered(grid%cell%ntot) &
        , cells(grid%cell%ntot)

        ! Auxiliary
        real(R8) :: fcX(grid%face%ntot), fcY(grid%face%ntot), sin1, sin2
        integer(I8) :: ic, v1, v2, f1, v1f, v2f, nv, i
        integer(I8), allocatable, dimension(:) :: tv, tf
        real(R8), allocatable, dimension(:) :: vec_vxx, vec_vxy, &
            vec_fcx, vec_fcy 

        ! Associate
        associate(&
            c  => grid%cell, &
            f  => grid%face, &
            v  => grid%vert &
            )

        ! Compute face centers
        fcX = 0.5_R8 * (v%x(f%vert(:,1)) + v%x(f%vert(:,2)))
        fcY = 0.5_R8 * (v%y(f%vert(:,1)) + v%y(f%vert(:,2)))

        ! Lopp over the cells
        is_ordered = .true.
        do ic = 1, c%ntot
            if (cells(ic)) then
                ! Get vertices and faces
                tv = GetCellVert(c, ic)
                tf = GetCellFace(c, ic)

                ! First check whether vertices and faces are connected according to correct ordering
                nv = size(tv)
                do i = 1, nv-1
                    v1 = tv(i)
                    v2 = tv(i+1)
                    f1 = tf(i)
                    v1f = f%vert(f1,1)
                    v2f = f%vert(f1,2)

                    if (((v1 /= v1f).and.(v1 /= v2f)).or. &
                          ((v2 /= v1f).and.(v2 /= v2f)) ) then
                            is_ordered(ic) = .false.
                            exit
                    end if 

                end do
                i = nv
                v1 = tv(i)
                f1 = tf(i)
                v1f = f%vert(f1,1)
                v2f = f%vert(f1,2)
                if (((v1 /= v1f).and.(v1 /= v2f)).or. &
                    ((tv(1) /= v1f).and.(tv(1) /= v2f)) ) then
                        is_ordered(ic) = .false.
                end if 


                ! If still ok, check the direction
                if (is_ordered(ic)) then
                    ! Vectors from cell center to vertices
                    vec_vxx = v%x(tv) - c%x(ic)
                    vec_vxy = v%y(tv) - c%y(ic)

                    ! Vectors from cell center to face centers
                    vec_fcx = fcX(tf) - c%x(ic)
                    vec_fcy = fcY(tf) - c%y(ic) 

                    do i = 1, nv-1
                        !check sinus between vector to vertex (a) and vector to face
                        !center (b)
                        !calculate angles (sin = |a x b| / norm(a)*norm(b))
                        !with |a x b | = ax*by - bx*ay
                        sin1 = vec_vxx(i)*vec_fcy(i) - vec_fcx(i) * vec_vxy(i)

                        ! Should be positive
                        if (sin1.lt.0.0_R8) then
                            is_ordered(ic) = .false.
                            exit
                        else 
                            !check sinus between vector to face center (b) and vector
                            !to second vertex (a)
                            sin2 = vec_fcx(i)*vec_vxy(i+1) - vec_vxx(i+1)*vec_fcy(i)
                            
                            if (sin2.lt.0.0_R8) then
                                is_ordered(ic) = .false.
                                exit
                            end if

                        end if

                    end do 
                end if

            end if
        end do

        end associate

    end subroutine

    subroutine ReorderCellConn(grid, is_ordered)

        ! Description
        !============
        ! Re-orders vertices and faces of cells in an order such that they form a
        ! chain as you would walk over the cell boundary 

        ! Declare variables
        !==================
        type(GridUDT), intent(inout)    :: grid
        logical, intent(in)             :: is_ordered(grid%cell%ntot)

        ! Declare variables
        !==================
        integer(I8) :: ic, nv, fcs(1:2), vs(1:2), s, n, i
        integer(I8), allocatable, dimension(:) :: tv, tf, ntv, ntf
        integer(I8), allocatable, dimension(:,:) :: vf
        real(R8) :: vec_start(1:2), vec_face1(1:2), vec_face2(1:2), &
            sin1, sin2
        real(R8), allocatable :: fcX(:), fcY(:)

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert  &
            )
        
        ! Compute face centers
        allocate(fcX(f%ntot), fcY(f%ntot))
        fcX = 0.5_R8 * (v%x(f%vert(:,1)) + v%x(f%vert(:,2)))
        fcY = 0.5_R8 * (v%y(f%vert(:,1)) + v%y(f%vert(:,2)))

        ! Loop over cells which are not ordered
        ! Could be implement to pack the cells and loop over them
        do ic = 1, c%ntot
            if (.not.is_ordered(ic)) then
                ! Get vertices and faces
                tv = GetCellVert(c, ic)
                tf = GetCellFace(c, ic)
                nv = size(tv)

                ! New vertices and faces
                ntv(1:nv) = 0
                ntf(1:nv) = 0

                ! Starting point
                ntv(1) = tv(1)

                vf = f%vert(tf,:)

                ! Find the two faces connected to that vertex and the other vertex connected to that face
            
                ! Find the first face
                fcs(1:2) = pack(vf, vf == ntv(1))

                ! Start right hand turning
                ! Define vector from cell centers to points
                vec_start(1) = v%x(ntv(1)) - c%x(ic)
                vec_start(2) = v%y(ntv(1)) - c%y(ic)

                ! Face1
                vec_face1(1) = fcX(fcs(1)) - c%x(ic)
                vec_face1(2) = fcY(fcs(1)) - c%y(ic)

                ! Face2
                vec_face2(1) = fcX(fcs(2)) - c%x(ic)
                vec_face2(2) = fcY(fcs(2)) - c%y(ic)
    
                ! Calculate angles (sin = |a x b| / norm(a)*norm(b))
                ! with |a x b | = ax*by - bx*ay
                ! dividing by norm is needed because norm is always positive and
                ! only the sign matters
                ! a = always vector from cell center to starting point
                sin1 = vec_start(1)*vec_face1(2) - vec_face1(1) * vec_start(2)
                sin2 = vec_start(1)*vec_face2(2) - vec_face2(1) * vec_start(2)
                
                ! Choose the one with positive sign, that right turning
                if (sin1 .gt. 0.0_R8) then
                    ntf(1) = fcs(1) ! First face
                elseif (sin2 .gt. 0.0_R8) then
                    ntf(1) = fcs(2)
                end if

                ! After the direction is fixed, implementation should be in a loop
                do i = 2, nv
                    ! Find the next vertex
                    vs = f%vert(ntf(i-1),:)

                    if (.not. any(ntv == vs(1))) then
                        ntv(i) = vs(1)
                    else
                        ntv(i) = vs(2)
                    end if

                    ! Find the next face
                    fcs(1:2) = pack(vf, vf == ntv(i))

                    if (.not. any(ntf == fcs(1))) then
                        ntf(i) = fcs(1)
                    else
                        ntf(i) = fcs(2)
                    endif

                end do

                ! Plug in the new verts and faces in cell%vert and cell%face
                s = c%vertP(ic,1)
                n = c%vertP(ic,2)
                c%vert(s:s+n-1) = ntv
                c%face(s:s+n-1) = ntf 

            end if

        end do

        ! Deallocate 
        deallocate(fcX, fcY)

        end associate


        
    end subroutine

    subroutine GetFsVxFromFsFc(grid)
        ! Description
        !============
        ! Get fsVx from fsFc

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(inout) :: grid

        ! Auxiliary
        integer(I8) :: fsVx(grid%vert%ntot), fsVxP(grid%data%fluxdata%nFs,2), &
            nv_counter, ifs, nv, nf , verts(1:grid%vert%ntot)
        integer(I8), allocatable, dimension(:) :: fcs, vxs 
         


        ! Initialize
        fsVx = 0
        fsVxP = 0
        nv_counter = 0
        verts = 0

        associate(&
            fd => grid%data%fluxdata, &
            f  => grid%face &
            )

        do ifs = 1, fd%nFs
            ! Get vertices from flux surface
            fcs = GetFSFace(fd, ifs)
            nf = size(fcs)
            verts(1:nf) = grid%face%vert(fcs,1) 
            verts(nf+1:nf*2) = grid%face%vert(fcs,2)
            call Unique(verts(1:nf*2), vxs)
            nv = size(vxs)

            ! Fill in into fsVx and fsVxP
            fsVxP(ifs,:) = [nv_counter+1, nv];
            fsVx(fsVxP(ifs,1):fsVxP(ifs,1)+fsVxP(ifs,2)-1) =  vxs;
            nv_counter = nv_counter + nv;        
        end do

        ! Trim - first implement GAGrid
        fd%fluxsurfaceverts = fsVx(1:nv_counter)
        fd%fluxsurfacevertsP = fsVxP(1:fd%nFs,:)

        end associate


    end subroutine

    subroutine GiveXpoint(grid,use_sep)

        ! Description
        !============
        ! Gives the Xpoint(s). Depending the grid information different methods are used. If all vertices have a fieldlineID the routine DetermineXPoints is used, otherwise the xpoint is determined by checking whether a vertex is surrounded by more then three regions and more than four cells. 

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(inout)   :: grid
        logical, intent(in)            :: use_sep              

        ! Auxiliary
        integer(I8)                         :: nxpind, i, j, iv, xpoints(1:100), counter, & 
            n , nc, s
        integer(I8), allocatable            :: xpind(:), order(:), vxs(:), cells(:), &
            regions(:), cvLookUp(:)
        logical :: use_sepID, start, use_fieldlineID, use_nsep
        

        ! Initialize
        counter = 0
        xpoints = 0
        

        ! Associate
        associate(&
            c   => grid%cell, &
            v   => grid%vert, &
            fd  => grid%data%fluxdata &
            )

            ! Find separatrix if not given
            if (.not.use_sep) then
                call GiveSeparatrices(grid,use_nsep,use_sepID,start,cvLookUp)
            end if

            ! Determine whether to use the vert.fieldlineID to determine the Xpoint
            use_fieldlineID = .false.
            if (allocated(v%fieldlineID)) then
                if (size(v%fieldlineID).eq.v%ntot) then
                    use_fieldlineID = .true.
                end if
            end if

            if (use_fieldlineID) then

                call DetermineXPoints(xpind, nxpind, order, grid)
                grid%data%xpointID  = xpind
                grid%data%nxp       = nxpind 

            elseif ((allocated(v%fieldlineID)).and.(allocated(v%cellP))) then

                ! Only check the vertices on the separatrices
                cvLookUp = GetCvLookUp(c)
                do i = 1, grid%data%nsep
                    s = fd%fluxsurfacevertsP(grid%data%sepID(i),1)
                    n = fd%fluxsurfacevertsP(grid%data%sepID(i),2)
                    vxs = fd%fluxsurfaceverts(s:s+n-1)

                    do j = 1, n
                        iv = vxs(j)
                        cells = GetVertCellGA(c, iv, cvLookUp)
                        call Unique(c%reg(cells), regions)
                        nc = size(cells)

                        if ((size(regions).ge.3) .and. (nc.gt.4)) then
                            counter = counter + 1
                            xpoints(counter) = iv
                        end if

                    end do 

                end do

                ! Saving 
                grid%data%xpointID  = xpoints(1:counter)
                grid%data%nxp       = counter

            else
                
                ! Check all vertices
                cvLookUp = GetCvLookUp(c)
                do iv = 1, v%ntot
                    cells = GetVertCellGA(c, iv, cvLookUp)
                    call Unique(c%reg(cells), regions)
                    nc = size(cells)

                    if ((size(regions).ge.3) .and. (nc.gt.4)) then
                        counter = counter + 1
                        xpoints(counter) = iv
                    end if

                end do 

                ! Saving 
                grid%data%xpointID  = xpoints(1:counter)
                grid%data%nxp       = counter

            end if

            if (grid%data%nxp == 0) then
                print *, 'GiveXpoints: Warning: no Xpoint found!'
            endif 

        end associate

    end subroutine 

    subroutine GiveSeparatrices(grid,use_nsep,use_sepID,start,cvLookUp)

        ! Description
        !============
        ! Determines the ID of the separatrices depending on the given information.
        ! The argument start is used when the separatrices need to be determined at the start of the program, so with little rest information.

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(inout)         :: grid
        logical, intent(in)                  :: use_nsep, use_sepID, start
        integer(I8), allocatable, optional :: cvLookUp(:)

        ! Initialize
        logical    :: use_xpointID     
        integer(I8), allocatable :: vxs(:), fcs(:), cvs(:)
        integer(I8) :: ifs, ix, nv, s, sepIDloc(1:100), &
            counter, counter_dummy, nf, step, i, ifc

        ! Initialize
        sepIDloc = 0
        counter = 0
        counter_dummy = 0
        

        ! Associate
        associate(&
            c   => grid%cell, &
            fd  => grid%data%fluxdata &
            )
            

            if (.not.present(cvLookUp)) then 
                cvLookUp = GetCvLookUp(c)
            end if
            
            ! Determine separatrices IDs
            if (start) then 
                if (.not.use_nsep) then

                    ! Checking if an Xpoint was already determined
                    if (allocated(grid%data%xpointID)) then
                        use_xpointID = .true.
                    end if

                    if (use_xpointID) then
                        ! Based on Xpoints
                        do ifs = 1, fd%nFs
                            nv = fd%fluxsurfacevertsP(ifs,1)
                            s = fd%fluxsurfacevertsP(ifs,2)
                            vxs = fd%fluxsurfaceverts(s:s+nv-1)

                            do ix = 1, grid%data%nxp
                                if (any(grid%data%xpointID(ix).eq.vxs)) then
                                    counter = counter + 1
                                    sepIDloc(counter) = ifs
                                end if
                            end do
                        end do

                    else

                        do ifs = 1, fd%nFs
                            nf = fd%fluxsurfacefacesP(ifs,1)
                            s = fd%fluxsurfacefacesP(ifs,2)
                            fcs = fd%fluxsurfacefaces(s:s+nf-1)

                            step = 1
                            sepIDloc = DetectSepID(ifs, c, fcs, sepIDloc, nf, step, cvLookUp, counter)
                            counter_dummy = counter + 1
                            if (sepIDloc(counter_dummy) /= 0) then
                                counter = counter + 1
                            end if 

                        end do 

                    end if

                    ! Save
                    grid%data%nsep = counter
                    grid%data%sepID = sepIDloc(1:counter)

                
                elseif (use_nsep .and. grid%data%nsep.eq.1) then
                    ! Determin separatrix for SN-case
                    ! Check whether previous sepID is still separatrix
                    if (use_sepID) then
                        nf = fd%fluxsurfacefacesP(grid%data%sepID(1),2)
                        s = fd%fluxsurfacefacesP(grid%data%sepID(1),1)
                        fcs = fd%fluxsurfacefaces(s:s+nf-1)
                        step = 1
                        do i = 1, step, nf
                            ifc = fcs(i)
                            cvs = GetFaceCellGA(c,ifc,cvLookUp)
                            if (size(cvs).eq.2) then
                                if (c%reg(cvs(1)) /= c%reg(cvs(2))) then
                                    ! Noting has changed and should not be possible to created new separatrices with grid adaptations
                                    return
                                end if
                            else
                                exit
                            end if
                        end do
                    end if 

                    ! Determine new separatrices (line 99)
                    sepIDloc = 0
                    do ifs = 1, fd%nFs
                        nf = fd%fluxsurfacefacesP(ifs,2)
                        s = fd%fluxsurfacefacesP(ifs,1)
                        fcs = fd%fluxsurfacefaces(s:s+nf-1)
                        step = 1

                        do i = 1, step, nf
                            ifc = fcs(i)
                            cvs = GetFaceCellGA(c,ifc,cvLookUp)
                            if (size(cvs).eq.2) then
                                if (c%reg(cvs(1)) /= c%reg(cvs(2))) then
                                    sepIDloc(1) = ifs
                                    exit
                                end if 
                            else
                                exit
                            end if
                        end do
                        if (sepIDloc(1) == ifs) then
                            grid%data%sepID = sepIDloc(1)
                            exit
                        end if

                    end do 

                    if (sepIDloc(1) .eq. 0) then
                        call gdErrorHandler('GiveSeparatrices: No separatrix found')
                    end if

                elseif (use_nsep.and.grid%data%nsep.gt.1) then
                    ! Determine separatrix for general case
                    sepIDloc = 0
                    counter = 0
                    ! Check whether previous is still correct
                    if (.not.any(grid%data%sepID /= 0)) then
                        do i = 1, grid%data%nsep
                            ifs = grid%data%sepID(i)
                            nf = fd%fluxsurfacefacesP(ifs, 2)
                            s = fd%fluxsurfacefacesP(ifs, 1)
                            fcs = fd%fluxsurfacefaces(s:s+nf-1)
                            step = 1
                            sepIDloc = DetectSepID(ifs, c, fcs, sepIDloc, nf, step, cvLookUp, counter)
                            counter_dummy = counter + 1
                            if (sepIDloc(counter_dummy) /= 0) then
                                counter = counter + 1
                            end if 
                        end do
                        if (counter.eq.grid%data%nsep) then
                            grid%data%sepID = sepIDloc(1:counter)
                            return
                        end if
                    end if

                    do ifs = 1, fd%nFs
                        if (.not.any(ifs == grid%data%sepID)) then
                            nf = fd%fluxsurfacefacesP(ifs, 2)
                            s = fd%fluxsurfacefacesP(ifs, 1)
                            fcs = fd%fluxsurfacefaces(s:s+nf-1)
                            step = 1
                            sepIDloc = DetectSepID(ifs, c, fcs, sepIDloc, nf, step, cvLookUp, counter)
                            counter_dummy = counter + 1
                            if (sepIDloc(counter_dummy) /= 0) then
                                counter = counter + 1
                            end if
                            if (sepIDloc(grid%data%nsep) /= 0 ) then ! all separatrices are found
                                exit
                            else 
                                call gdErrorHandler("GiveSeparatrices: not all separatrices were found.")
                                ! If above would not be sufficient, check give_iFs_sep.m line 194
                            end if
                        end if

                    end do

                    ! Save
                    grid%data%sepID = sepIDloc(1:counter)
                    !(nsep already there)

                end if 

            else ! no start, hopefully more efficient
                sepIDloc = 0
                counter = 0
                if  (.not.use_nsep) then
                    if (allocated(grid%data%xpointID)) then
                        do ifs = 1, fd%nFs
                            nv = fd%fluxsurfacevertsP(ifs,2)
                            s = fd%fluxsurfacevertsP(ifs,1)

                            vxs = fd%fluxsurfaceverts(s:s+nv-1)
                            do ix = 1, grid%data%nxp
                                if (any(grid%data%xpointID(ix)==vxs)) then
                                    counter = counter + 1
                                    sepIDloc(counter) = ifs
                                end if 
                            end do 
                        end do
                    else ! line 252
                        do ifs = 1, fd%nFs
                            nf = fd%fluxsurfacefacesP(ifs, 2)
                            s = fd%fluxsurfacefacesP(ifs, 1)
                            fcs = fd%fluxsurfacefaces(s:s+nf-1)
                            step = max(2,nint(nf/5.0_R8))
                            sepIDloc = DetectSepID(ifs,c, fcs, sepIDloc, nf, step, cvLookUp, counter)
                            counter_dummy = counter + 1
                            if (sepIDloc(counter_dummy) /= 0) then
                                    counter = counter + 1
                            end if

                        end do

                    end if

                    ! Save
                    grid%data%nsep = counter
                    grid%data%sepID = sepIDloc(1:counter)
                elseif (use_nsep.and.(grid%data%nsep.eq.1)) then
                    ! Check whether the previous is separatrix
                    if (grid%data%sepID(1) /= 0) then
                            nf = fd%fluxsurfacefacesP(ifs, 2)
                            s = fd%fluxsurfacefacesP(ifs, 1)
                            fcs = fd%fluxsurfacefaces(s:s+nf-1)
                            step = max(2,nint(nf/5.0_R8))   
                            do i = 3, step, nf-2
                                ifc = fcs(i)
                                cvs = GetFaceCellGA(c,ifc,cvLookUp)
                                if (size(cvs).eq.2) then
                                    if (c%reg(cvs(1))/=c%reg(cvs(2))) then
                                        return ! Separatrix is the same
                                    end if
                                else
                                    exit ! exit loop because this surface is boundary
                                end if

                            end do                     
                    end if

                    ! Determine new separatrix if necessary
                    sepIDloc = 0
                    counter = 0
                    do ifs = 1, fd%nFs
                        nf = fd%fluxsurfacefacesP(ifs, 2)
                        s = fd%fluxsurfacefacesP(ifs, 1)
                        fcs = fd%fluxsurfacefaces(s:s+nf-1)
                        step = 1   
                        do i = 3, step, nf-2 
                            ifc = fcs(i)
                            cvs = GetFaceCellGA(c, ifc, cvLookUp)
                            if (size(cvs).eq.2) then 
                                if (c%reg(cvs(1))/=c%reg(cvs(2))) then
                                    counter = counter + 1
                                    sepIDloc(counter) = ifs
                                    exit
                                end if
                            else
                                exit                               
                            end if  
                        end do
                        if (sepIDloc(counter).eq.ifs) then
                            exit
                        end if
                    end do 

                    ! Save
                    grid%data%sepID = sepIDloc(1)

                elseif (use_nsep.and.(grid%data%nsep.gt.1)) then
                    sepIDloc = 0
                    counter = 0
                    if (.not.any(grid%data%sepID.eq.0)) then
                        do i = 1, grid%data%nsep
                            ifs = grid%data%sepID(i)
                            nf = fd%fluxsurfacefacesP(ifs, 2)
                            s = fd%fluxsurfacefacesP(ifs, 1)
                            fcs = fd%fluxsurfacefaces(s:s+nf-1) 
                            step = 1
                            sepIDloc = DetectSepID(ifs,c, fcs, sepIDloc, nf, step, cvLookUp, counter)
                            counter_dummy = counter + 1
                            if (sepIDloc(counter_dummy) /= 0) then
                                    counter = counter + 1
                            end if                           
                        end do
                        if (counter .eq. grid%data%nsep) then
                            return ! All stays the same
                        end if
                    end if

                    ! Determine new separatrices if necessary
                    do ifs = 1, fd%nFs
                        if (.not.any(sepIDloc(1:grid%data%nsep).eq.ifs)) then
                            nf = fd%fluxsurfacefacesP(ifs, 2)
                            s = fd%fluxsurfacefacesP(ifs, 1)
                            fcs = fd%fluxsurfacefaces(s:s+nf-1) 
                            step = 1   
                            sepIDloc = DetectSepID(ifs,c, fcs, sepIDloc, nf, step, cvLookUp, counter)
                            counter_dummy = counter + 1
                            if (sepIDloc(counter_dummy) /= 0) then
                                    counter = counter + 1
                            end if  
                            if (sepIDloc(grid%data%nsep)/= 0) then
                                exit ! all separatrices are found
                            end if
                        end if
                    end do

                    ! If not all where found - extra check possible see line 390 in give_iFs_sep.m
                    if (sepIDloc(grid%data%nsep).eq.0) then
                       call gdErrorHandler("GiveSeparatrices: not all separatrices where found")
                    end if

                    ! Save
                    grid%data%sepID = sepIDloc(1:counter)

                end if

            end if 

        end associate

    end subroutine

    subroutine IdentifyAlignedFaces(grid,options,magneticField)
       
        ! Description
        !============
        ! Identify adhoc whether faces are aligned with the magnetic field or
        ! whether they were supposed to be aligned by the grid generator

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT) :: grid
        type(GAoptionsUDT) :: options
        type(MagneticFieldUDT) :: magneticField
        
        ! Auxiliary
        integer(I8) :: ifc, v1, v2
        integer(I8), allocatable :: facealigned(:)
        real(R8), allocatable :: fcX(:), fcY(:), &
            dpsidx(:), dpsidy(:), Bx(:), By(:), abs_cos(:), &
            t1x(:), t2x(:), t1y(:), t2y(:)

        ! Associate
        associate( &
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        ! Initialize
        allocate(facealigned(f%ntot), fcX(f%ntot), fcY(f%ntot), &
        Bx(f%ntot), By(f%ntot), tx(f%ntot,2), Bnorm(f%ntot), cos(f%ntot))
        fcX = 0.5_R8*(v%x(f%vert(:,1)) + v%x(f%vert(:,2)))
        fcY = 0.5_R8*(v%y(f%vert(:,1)) + v%y(f%vert(:,2)))

        ! Sort faces - Already done in GAInit

        ! Get Magneticfield
        call magneticField%interp%Evaluate(fcX,fcy,1,0,dpsidx)
        call magneticField%interp%Evaluate(fcX,fcy,0,1,dpsidy)  
        Bx = -dpsidy
        By = dpsidx
        Btot = sqrt(Bx**2 + By**2)

        abs_cos = 0
        ! Tangential vector
        t1x = v%x(f%vert(:,1))
        t2x = v%x(f%vert(:,2))
        t1y = v%y(f%vert(:,1))
        t2y = v%y(f%vert(:,2))

        t1x = t2x - t1x
        t1y = t2y - t1y

        Bnorm = sqrt(t1x**2 + t1y**2)*Btot
        cos = t1x*Bx + t1y*By
        
        ! Loop because of intrinsic abs an
        do ifc = 1, f%ntot
            abs_cos(ifc) = abs(max(-1.0_R8,min(1.0_R8,cos(ifc))) / Bnorm(ifc))
        end do

        ! Locally generalize face labels
        
              




        end associate

    end subroutine

    !==================================================================!
    !                                                                  !
    !                           FUNCTIONS                              !
    !                                                                  !
    !==================================================================!    

    ! Get cells of a vertex without using vert%cell
    function GetVertCellGA(cell, i, cvLookUp) result(res)
        integer(I8)                 :: i 
        type(CellUDT)               :: cell
        integer(I8), allocatable, optional    :: cvLookUp(:)
        integer(I8), allocatable    :: res(:)
        
        if (.not.present(cvLookUp)) then
            cvLookUp = GetCvLookUp(cell)
        end if
            
        res = pack(cvLookUp,cell%vert.eq.i)
    end function

    function GetFaceCellGA(cell, i, cvLookUp) result(res)
        integer(I8)                 :: i 
        type(CellUDT)               :: cell
        integer(I8), allocatable, optional    :: cvLookUp(:)
        integer(I8), allocatable    :: res(:)  
        
        if (.not.present(cvLookUp)) then
            cvLookUp = GetCvLookUp(cell)
        end if

        res = pack(cvLookUp,cell%face.eq.i)
    end function

    function GetFluxSurfaceFcs(fd, i) result(res)
        integer(I8) :: i, s, nf
        type(FluxDataUDT) :: fd
        integer(I8), allocatable :: res(:)

        nf = fd%fluxsurfacefacesP(i,2)
        s = fd%fluxsurfacefacesP(i,1)
        res = fd%fluxsurfacefaces(s:s+nf-1)

    end function

    function GetCvLookUp(cell) result(res)
        type(CellUDT)       :: cell
        integer(I8)         :: nc, ic, nv, s                
        integer(I8), allocatable :: res(:)

        nc = cell%ntot

        allocate(res(1:cell%vertP(nc,1)+cell%vertP(nc,2)-1))
        res = 0

        do ic = 1, nc
            s = cell%vertP(ic, 1)
            nv = cell%vertP(ic, 2)
            res(s:s+nv-1) = ic
        end do
    end function

    function DetectSepID(ifs, cell, faces, sepIDloc, nf, step, cvLookUp, counter) result(res)
        ! Description
        !============
        ! Check whethere is a core region on one side of the flux surface
        integer(I8) :: ifs, nf, counter, i, step, ifc, reg1, reg2, sepIDloc(1:100)
        type(CellUDT) ::  cell
        integer(I8), allocatable :: cvLookUp(:), faces(:), cvs(:), res(:)

        ! Initialize
        res = sepIDloc

        do i = 1, step, nf
            ifc = faces(i)
            cvs = GetFaceCellGA(cell,ifc, cvLookUp)
            if (size(cvs).eq.2) then
                reg1 = cell%reg(cvs(1))
                reg2 = cell%reg(cvs(2))
                if (reg1 /= reg2) then
                    if ((mod(reg1,4).eq.1) .or. (mod(reg2,4).eq.1)) then 
                        counter = counter + 1
                        res(counter) = ifs
                        exit
                    end if
                end if
            end if
        end do


        
    end function


end module 