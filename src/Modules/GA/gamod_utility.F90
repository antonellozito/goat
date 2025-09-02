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
    use mod_dynamicarrays
    use goatmod_types
    use goatmod_userinput
    !use gdmod_utility_optimization
    use gamod_types


    ! The usual
    implicit none
    save
    public
    
    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================! 
    
    contains 


    subroutine BuildFluxTubeData(grid, options, magneticField, ctype) 

        ! Description
        !============
        ! Build the flux tube data, necessary for b2fgmtry
        ! Probleem flux data also from guard cells => no here yet

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(inout)    :: grid
        type(GAoptionsUDT), intent(in)  :: options
        type(MagneticFieldUDT)          :: magneticField
        character(:), allocatable       :: ctype

        ! Auxiliary
        logical :: ft_open, ft_open_us, ft_closed
        integer(I8), allocatable :: cvs(:), cvs_rev(:), ind(:), fcs(:), fcs_rev(:), &
             ftCv(:), ftCvP(:,:), ftFc(:), ftFcP(:,:)
        integer(I8) :: nc, c1, c2, sv, nv, sf, nf, i, ift, ic, indc
        real(R8) :: c1_bx, c1_by, d, c1_bpolx, c1_bpoly, con_x, con_Y, d2, cosin, &
            dpsidx(grid%cell%ntot), dpsidy(grid%cell%ntot)

        ! Initialize
        ft_open = .false.
        ft_open_us = .false.
        ft_closed = .false.


        ! Determine method
        select case (ctype)

        case ('all')

            ft_open = .true.
            ft_closed = .true.

        case ('closed')

            ft_closed = .true.

        case ('open')

            ft_open = .true.

        case ('all_us')

            ft_open_us = .true.
            ft_closed = .true.

        case default

            call gdErrorHandler("BuildFluxTubeData: type not defined")

        end select

        ! Build open tubes
        !=================
        if (ft_open) then

            ! Build the open tubes connecting the target for structured grids   
            call BuildOpenTubes(grid, options, ftCv, ftCvP)

        else if (ft_open_us) then

            ! Build the open tubes connecting the target for unstructured grids
            call BuildOpenTubesUS(grid, options, ftCv, ftCvP, ftFc, ftFcP)

        end if

        ! Build closed tubes
        !===================
        if (ft_open .and. ft_closed) then

            ! Build the core flux tubes while open flux tubes were build
            call BuildClosedTubes(grid, options, ft_open_us, ftCv, ftCvP, ftFc, ftFcP)

        elseif (ft_open_us .and. ft_closed) then

            ! Build the core flux tubes while open flux tubes were build and add ftFc
            call BuildClosedTubes(grid,options,ft_open_us, ftCv, ftCvP, ftFc, ftFcP)

        elseif (.not.(ft_open .or. ft_open_us) .and. ft_closed) then 

            call gdErrorHandler('BuildFluxTubeData: only building the closed and the separatrix tube is not supported anymore')

            !call BuildClosedTubes(grid,options,ft_open_us, ftCv, ftCvP, ftFc, ftFcP, first_core_tube)
            !call AddSeparatrixTube(grid, ftCv, ftCvP, ftFc, ftFcP, first_core_tube)            

        end if

        ! Build ftFc - poloidal faces which lay into the flux tube
        if (.not.ft_open_us) then
            call BuildFtFc(grid)
        end if


        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert, &
            fd => grid%data%fluxdata &
        )

        ! Store
        fd%nFt = count(ftCvP(:,1) /= 0)
        fd%fluxtubecells = ftCv(1:ftCvP(fd%nFt,1)+ftCvP(fd%nFt,2)-1)
        fd%fluxtubecellsP = ftCvP(1:fd%nFt,:)
        fd%fluxtubefaces = ftFc(1:ftFcP(fd%nFt,1)+ftFcP(fd%nFt,2)-1)
        fd%fluxtubefacesP = ftFcP(1:fd%nFt,:)


        ! Orden the flux tubes along positive poloidal magnetic field
        call magneticField%interp%Evaluate(c%x,c%y,1,0,dpsidx)
        call magneticField%interp%Evaluate(c%x,c%y,0,1,dpsidy)
        do ift = 1,grid%data%fluxdata%nFt

            cvs = GetFTCell(fd, ift)
            nc = size(cvs)

            if (nc .gt. 1) then

                c1 = cvs(1)
                c2 = cvs(2)

                ! Compute direction
                c1_bx = dpsidx(c1)
                c1_by = dpsidy(c1)

                ! Get unit poloidal magneticField vector
                d = sqrt(c1_bx**2 + c1_by**2)
                c1_bpolx = -c1_by/d
                c1_bpoly = c1_bx/d

                ! Compute connector from first to second cell in the tube
                con_x = c%x(c2) - c%x(c1)
                con_y = c%y(c2) - c%y(c1)
                d2 = sqrt(con_x**2 + con_y**2)
                con_x = con_x/d2
                con_y = con_y/d2    

                ! Compute cosine between poloidal magnetic field and connector
                cosin = c1_bpolx*con_x + c1_bpoly*con_y 
                
                ! Reverse the ordening if cos < 0
                if (cosin .lt. 0.0_R8) then

                    ! Cells
                    ind = (/ (i, i = nc,1,-1) /)
                    cvs_rev = cvs(ind)

                    sv = fd%fluxtubecellsP(ift,1)
                    nv = fd%fluxtubecellsP(ift,2)  
                    
                    ind = (/ (i, i = sv, sv+nv-1)/)
                    fd%fluxtubecells(ind) = cvs_rev


                    ! Faces
                    sf = fd%fluxtubefacesP(ift,1)
                    nf = fd%fluxtubefacesP(ift,2)

                    fcs = fd%fluxtubefaces(sf:sf+nf-1)

                    ind = (/ (i, i = nf,1,-1) /)
                    fcs_rev = fcs(ind)

                    ind = (/ (i, i = sf, sf+nf-1)/)
                    fd%fluxtubefaces(ind) = fcs_rev
                    deallocate(ind)

                end if
            end if 
        end do

        ! Make cvFt -- inverse of ftCv
        c%ft = 0
        do ic = 1, c%ntot
            indc = findloc(fd%fluxtubecells, ic, 1)
            if (indc == 0) then
                print *, 'BuildFluxTubeData: cell without flux tube'
                c%ft(ic) = fd%nFt + 1 ! Convection to avoid zeros
            else 
                ift = GetFluxTubeFromCellIndex(fd, indc)
                c%ft(ic) = ift
            end if
        end do

    
        ! Make ftReg
        ! 1 = core
        ! 2 = SOL
        ! 3 = PF
        fd%fluxtuberegID = 0
        do ift = 1, fd%nFt
            cvs = GetFTCell(fd, ift)
            fd%fluxtuberegID(ift) = minval(c%reg(cvs))
        end do 

        end associate


    end subroutine

    subroutine BuildOpenTubes(grid, options, ftCv, ftCvP)

        ! Description
        !============
        ! Builds the flux tubes which are open.
        ! This means all the flux tubes expect the ones in the core region.
        
        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(in)           :: grid
        type(GAoptionsUDT), intent(in)      :: options
        integer(I8), allocatable, intent(out)  ::  ftCv(:), ftCvP(:,:)

        ! Auxiliary
        integer(I8), allocatable :: tube(:), fcLbl_loc(:), indf(:), &
             indFc(:), tf(:), fcs(:), cvs(:), ft_cells(:), fcs_next(:), &
             ind_face11(:)
        integer(I8) :: tube_count, start_cell, indmax, fc1, fc2, counter, &
            ic_next, ind_face1, ind_opp_face, s, k, i, cvs1(1:2), ic, ift, nf
        real(R8), allocatable :: dpsi_f(:)
        logical :: closed


        ! Associate
        associate(& 
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert, &
            fd => grid%data%fluxdata &
            )
        allocate(fcLbl_loc(f%ntot))
        allocate(indFc(f%ntot))
        allocate(ftCv(c%ntot))
        allocate(ftCvP(c%ntot,2)) 
        allocate(tube(c%ntot))
        ftCv = 0
        ftCvP = 0
        tube = 0
        tube_count = 0
        counter = 0

        ! Get GA face labels
        fcLbl_loc = f%label
        indFc = (/ (i, i=1,f%ntot) /) 
        do i = 1, size(options%facelabelmappingGG)
            fcLbl_loc(pack(indFc, f%label == options%facelabelmappingGG(i))) = options%facelabelmappingGA(i)
        end do

        ! Loop over all cells
        do ic = 1, c%ntot

            ! Flux tubes between the target
            ! Start at the boundary
            if (c%cflags(ic) == 3) then

                tf = GetCellFace(c, ic)
                nf = c%faceP(ic,2)

                start_cell = 0
                do i = 1, nf
                    if ((fcLbl_loc(tf(i)) .eq. 4) .or. &
                     (fcLbl_loc(tf(i)) .eq. 2)) then
                        start_cell = 1
                    end if
                end do
                if (.not.any(ic == ftCv) .and. start_cell.eq.1) then

                    ! Start tube
                    tube(1) = ic

                    ! Find other cells
                    fcs = c%face(c%faceP(ic,1):c%faceP(ic,1)+c%faceP(ic,2)-1)

                    dpsi_f = abs(v%psi(f%vert(fcs,1)) - v%psi(f%vert(fcs,2)))
                    indmax = maxloc(dpsi_f,1)
                    fc1 = fcs(indmax) ! poloidal face

                    if (c%faceP(ic,2) .eq. 4) then ! quad

                        dpsi_f(indmax) = 0 ! eliminate previous face
                        indmax = maxloc(dpsi_f,1)
                        fc2 = fcs(indmax) ! Get second poloidal face

                        if (fcLbl_loc(fc1) /= 0) then ! boundary face?
                            fc1 = fc2
                        end if 

                    else if (c%faceP(ic,2) .eq. 3) then

                        if (fcLbl_loc(fc1) /= 0) then
                            dpsi_f(indmax) = 0
                            indmax = maxloc(dpsi_f,1)
                            fc1 = fcs(indmax)

                        end if 

                    end if

                    ! Marching the direction of the max dpsi face
                    if (fcLbl_loc(fc1) == 0) then

                        cvs = f%cell(f%cellP(fc1,1):f%cellP(fc1,1)+f%cellP(fc1,2)-1)
                        do k = 1,2
                            if (.not.any(cvs(k) == tube(1:counter))) then
                                counter = counter + 1
                                tube(counter) = cvs(k)
                                ic_next = cvs(k)
                                exit                                 
                            end if 
                        end do 

                        closed = .false. ! Flag for core flux tubes
                        do while (fcLbl_loc(fc1) == 0) 
                            
                            ! Get opposite face
                            nf =  c%faceP(ic_next,2)
                            fcs_next = c%face(c%faceP(ic_next,1):c%faceP(ic_next,1)+nf-1);

                            indf = (/ (i, i=1,nf )/)
                            ind_face11 = pack(indf, fcs_next == fc1)
                            ind_face1 = ind_face11(1)
                            ind_opp_face = ind_face1 + 2

                            if (ind_opp_face .gt. nf) then
                                ind_opp_face = ind_opp_face - nf
                            end if

                            fc1 = fcs_next(ind_opp_face)
                            if (f%cellP(fc1,2) == 2) then

                                cvs = f%cell(f%cellP(fc1,1):f%cellP(fc1,1)+f%cellP(fc1,2)-1)
                                do k = 1, 2
                                    if (.not.any(cvs(k) == tube(1:counter))) then

                                        ! Add cell to tube
                                        counter = counter + 1
                                        tube(counter) = cvs(k)
                                        ic_next = cvs(k)
                                        exit  ! only add one of the cells
                                         
                                    end if
                                end do

                                if (fc1 == fc2) then

                                    closed = .true.

                                end if 


                            end if

                        end do 

                        ! Save the tube
                        tube_count = tube_count + 1
                        if (tube_count == 1) then

                            ftCvP(tube_count,1) = 1

                        else

                            ftCvP(tube_count,1) = ftCvP(tube_count-1,1) + ftCvP(tube_count-1,2)

                        end if

                        ftCvP(tube_count,2) = counter
                        s = ftCvP(tube_count,1)
                        ftCv(s:s+counter-1) = tube(1:counter)

                        ! Reset
                        tube = 0
                        counter = 0
                        closed = .false.

                    end if


                end if



            end if 

        end do

        ! Check for missing cells in the beginning of a flux tube
        do ift = 1, tube_count

            ! Get first cell of tube
            ft_cells = ftCv(ftCvP(ift,1):ftCvP(ift,1)+ftCvP(ift,2)-1)
            ic = ft_cells(1)

            ! Look on both sides in the  poloidal directions
            fcs = GetCellFace(c, ic)

            dpsi_f = abs(v%psi(f%vert(fcs,1)) - v%psi(f%vert(fcs,2)))
            indmax = maxloc(dpsi_f,1)
            fc1 = fcs(indmax) ! poloidal face1

            if ( c%faceP(ic,2).eq. 4) then ! Quad
                dpsi_f(indmax) = 0
                indmax = maxloc(dpsi_f,1)
                fc2 = fcs(indmax) ! poloidal face1
            end if

            ! Poloidal face 1 direction
            if (fcLbl_loc(fc1) == 0) then ! not a boundary face

                ! Get cells of face
                cvs1 = f%cell(f%cellP(fc1,1):f%cellP(fc1,1)+1)
                do k = 1, 2

                    if ((cvs1(k) /= ic) .and. (.not.any(cvs1(k) == ft_cells))) then

                        ! Add the cell in front of the tube
                        ! Push every value from ftCv(ift,1) one place futher
                        ftCv(ftCvP(ift,1)+1:ftCvP(tube_count,1)+ftCvP(tube_count,2)) = &
                          ftCv(ftCvP(ift,1):ftCvP(tube_count,1)+ftCvP(tube_count,2)-1)
                        ftCv(ftCvP(ift,1)) = cvs1(k)
                        ftCvP(ift,2) = ftCvP(ift,2)+1
                        ftCvP(ift+1:tube_count,1) = ftCvP(ift+1:tube_count,1) + 1

                        exit

                    end if

                end do

            end if

            if (c%faceP(ic,2) .eq. 4) then

                if (fcLbl_loc(fc2) == 0) then  ! not a boundary face
                    ! Get cells of face
                    cvs1 = f%cell(f%cellP(fc2,1):f%cellP(fc2,1)+1)

                    do k = 1, 2

                        if ((cvs1(k) /= ic) .and. (.not.any(cvs1(k) == ft_cells))) then

                            ! Add the cell in front of the tube
                            ! Push every value from ftCv(ift,1) one place futher
                            ftCv(ftCvP(ift,1)+1:ftCvP(tube_count,1)+ftCvP(tube_count,2)) = &
                            ftCv(ftCvP(ift,1):ftCvP(tube_count,1)+ftCvP(tube_count,2)-1)
                            ftCv(ftCvP(ift,1)) = cvs1(k)
                            ftCvP(ift,2) = ftCvP(ift,2)+1
                            ftCvP(ift+1:tube_count,1) = ftCvP(ift+1:tube_count,1) + 1

                            exit

                        end if
                
                    end do

                end if 
            
            end if 

        end do

        end associate
    end subroutine

    subroutine BuildOpenTubesUS(grid, options, ftCv, ftCvP, ftFc, ftFcP)

        ! Description
        !============
        ! Builds the flux tubes which are open.
        ! This means all the flux tubes expect the ones in the core region.
        
        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(in) :: grid
        type(GAoptionsUDT), intent(in) :: options
        integer(I8), allocatable, intent(out) :: ftCv(:), ftCvP(:,:), ftFc(:), ftFcP(:,:)

        ! Auxiliary
        integer(I8), allocatable, dimension(:) :: tube, tube_fc, fcs, & 
            pfaces, pfacesB, pfacesI, pfaces_dummy, pfaces_dummy2, cvs, faceB, &
            ic_dummy, ipface_dummy, ipface_dummy2
        integer(I8) :: tube_count, ic, ic1, npfacesB, i, tube_cv, tube_count_fc, &
            cv_tot_count, fc_tot_count, ifc, ind, next_ic, v1, v2, nipface, ipface, nt
        real(R8), allocatable :: fcX(:), fcY(:), vecX_i, vecY_i, vecX_B, vecY_B, &
            cosf(:)
        logical :: start_cell,  build_tube
        logical, allocatable :: in_ftCv(:), test(:)


        ! Associate
        associate(& 
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert, &
            fd => grid%data%fluxdata &
            )

        ! Initialize
        allocate(ftCv(c%ntot))
        allocate(ftCvP(c%ntot,2)) ! NFT not known probably
        allocate(ftFc(f%ntot))
        allocate(ftFcP(f%ntot,2))        
        allocate(tube(c%ntot))
        allocate(tube_fc(f%ntot))
        allocate(in_ftCv(c%ntot))
        allocate(fcX(f%ntot))
        allocate(fcY(f%ntot))
        ftCv = 0
        ftCvP = 0
        ftFc = 0
        ftFcP = 0
        tube = 0
        tube_count = 0
        tube_fc = 0
        tube_cv = 0
        in_ftCv = .false.
        cv_tot_count = 0
        fc_tot_count = 0
        fcX = 0.5_R8 * (v%x(f%vert(:,1)) + v%x(f%vert(:,2)))
        fcY = 0.5_R8 * (v%y(f%vert(:,1)) + v%y(f%vert(:,2)))

        ! Loop over all cells
        do ic = 1, c%ntot

            ! Find a starting cell
            if (c%cflags(ic) == 3) then

                if (.not.in_ftCv(ic)) then

                    ! Dummy variable to save the cell
                    ic1 = ic

                    ! Get faces of the cell
                    fcs = GetCellFace(c, ic)

                    ! Eliminate if cells is trapezoid
                    start_cell = .true.
                    if (allocated(pfaces)) deallocate(pfaces)
                    allocate(pfaces(count(f%aligned(fcs) == 0)))
                    pfaces = pack(fcs, f%aligned(fcs) == 0)

                    if (allocated(pfacesB)) deallocate(pfacesB)
                    allocate(pfacesB(count(f%label(pfaces) /= 0)))
                    pfacesB = pack(pfaces, f%label(pfaces) /= 0)
                    npfacesB = size(pfacesB)

                    ! Get the internal faces
                    if (npfacesB .gt. 1) then
                        if (npfacesB == 2) then
                            allocate(pfaces_dummy(count(pfaces /= pfacesB(1))))
                            pfaces_dummy = pack(pfaces, pfaces /= pfacesB(1)) 

                            allocate(pfaces_dummy2(count(pfaces_dummy /= pfacesB(2))))
                            pfaces_dummy2 = pack(pfaces_dummy, pfaces_dummy /= pfacesB(2))

                            allocate(pfacesI(size(pfaces_dummy2)))
                            pfacesI = pfaces_dummy2

                            deallocate(pfaces_dummy)
                            deallocate(pfaces_dummy2)
                        else
                            call gdErrorHandler('BuildOpenTubesUS: cell with three boundary faces not yet implemented')
                        end if

                    elseif (npfacesB == 1) then
                        if (allocated(pfacesI)) deallocate(pfacesI)
                        allocate(pfacesI(count(pfaces /= pfacesB(1))))
                        pfacesI = pack(pfaces, pfaces /= pfacesB(1))
                    end if

                    ! Do not start the tube if following conditions
                    if ((npfacesB == 0).or. (size(pfaces) == 3 .and. size(pfacesB) == 1)) then
                        start_cell = .false.
                    end if

                    ! Sticking out triangle tube (one cell)
                    if ((size(pfaces) == 2) .and. (size(pfacesB) == 2)) then

                        start_cell = .false.

                        ! Get the tube of one cell
                        tube(1) = ic1
                        tube_cv = 1

                        tube_fc(1) = pfaces(1)
                        tube_fc(2) = pfaces(2)
                        tube_count_fc = 2

                        ! Save
                        tube_count = tube_count + 1

                        ftCvP(tube_count,2) = tube_cv
                        cv_tot_count = cv_tot_count + 1
                        ftCv(cv_tot_count) = tube(1)


                        ftFcP(tube_count,2) = tube_count_fc
                        
                        ftFc(fc_tot_count+1:fc_tot_count + tube_count_fc) = tube_fc(1:tube_count_fc)
                        fc_tot_count = fc_tot_count + tube_count_fc

                        if (tube_count == 1) then
                            ftCvP(1,1) = 1
                            ftFcP(1,1) = 2
                        else
                            ftCvP(tube_count,2) = ftCvP(tube_count-1,1) + ftCvP(tube_count-1,2)
                            ftFcP(tube_count,2) = ftFcP(tube_count-1,1) + ftFcP(tube_count-1,2)
                        end if

                        in_ftCv(tube_fc(1:tube_count_fc)) = .true.
                        tube = 0
                        tube_fc = 0

                    end if

                    ! Regular tube 
                    if (start_cell) then

                        ! Add this cell to the tube
                        tube(1) = ic1
                        build_tube = .true.
                        tube_cv = 1

                        ! Add faces
                        npfacesB = size(pfacesB)
                        if (npfacesB .gt. 1) then ! Trapezoid at boundary

                            ! Pick the most aligned face with the internal face
                            allocate(cosf(npfacesB))
                            cosf = 0

                            ! Tangential vector to internal face
                            v1 = f%vert(pfacesI(1),1)
                            v2 = f%vert(pfacesI(1),2)
                            vecX_i = v%x(v2) - v%x(v1)
                            vecY_i = v%y(v2) - v%y(v1)
                            
                            ! Compute cosines of other faces
                            do i = 1, npfacesB

                                ifc = pfacesB(i)
                                v1 = f%vert(ifc,1)
                                v2 = f%vert(ifc,2)
                                vecX_B = v%x(v2) - v%x(v1)
                                vecY_B = v%y(v2) - v%y(v1)

                                cosf(i) = vecX_i*vecX_B + vecY_i*vecY_B

                            end do

                            ! Get most aligned face
                            ind = maxloc(abs(cosf),1)

                            ! Attritube faces in tube
                            tube_fc(1) = pfacesB(ind)
                            tube_fc(2) = pfacesI(1)

                            deallocate(cosf)
                        else

                            ! Attritube faces in tube
                            tube_fc(1) = pfacesB(1)
                            allocate(pfaces_dummy(count(pfaces /= pfacesB(1))))
                            pfaces_dummy = pack(pfaces,pfaces /= pfacesB(1))
                            tube_fc(2) = pfaces_dummy(1)
                            deallocate(pfaces_dummy)

                        end if 

                    
                        ! Housekeeping
                        deallocate(pfaces)
                        deallocate(pfacesB)

                        ! Update counter
                        tube_count_fc = 2

                        ! Marching algo
                        ! Get the internal poloidal faces
                        allocate(ipface_dummy(size(pfacesI)))
                        ipface_dummy = pfacesI
                        if (size(ipface_dummy) /= 1) then
                            build_tube = .false.
                        else 
                            ipface = ipface_dummy(1)
                        end if
                        deallocate(ipface_dummy)

                        if (build_tube) then

                            do while (f%label(ipface) == 0) 

                                ! Get cells of ipface
                                cvs = GetFaceCell(f,ipface)
                                allocate(ic_dummy(count(cvs /= ic1)))
                                ic_dummy = pack(cvs,cvs /= ic1)
                                next_ic = ic_dummy(1)
                                deallocate(ic_dummy)

                                if (any(next_ic == tube)) then

                                    exit

                                else

                                    tube_cv = tube_cv + 1
                                    tube(tube_cv) = next_ic

                                end if

                                ! Get ipface
                                ic1 = next_ic

                                fcs = GetCellFace(c, ic1)
                                if (allocated(pfaces)) deallocate(pfaces)
                                allocate(pfaces(count(f%aligned(fcs) == 0)))
                                pfaces = pack(fcs, f%aligned(fcs) == 0)

                                if (.not.any(pfaces == ipface)) then
                                    call gdErrorHandler('Something wrong')
                                end if

                                allocate(ipface_dummy(count(pfaces /= ipface)))
                                ipface_dummy = pack(pfaces, pfaces /= ipface)  

                                ! Eliminate possible poloidal boundary faces
                                nipface = size(ipface_dummy)
                                if (nipface .ge. 1) then
                                    if (nipface .gt. 1) then

                                        ! Eliminate a internal poloidal face
                                        if (count(f%label(ipface_dummy) /= 0) == nipface) then

                                                ! Poloidal faces are all boundary faces
                                                ipface = ipface_dummy(1)

                                        else
                                                ! Only keep the non-boundary ipface
                                                !if (count(f%label(ipface_dummy) == 0) /= 1) then
                                                test = isBoundaryFace1D(ipface_dummy,f)
                                                nt = count(.not.test)
                                                if (count(.not.test) /= 1) then 
                                                    call gdErrorHandler('BuildOpenTubeUS: not supported case')
                                                end if

                                                allocate(ipface_dummy2(count(f%label(ipface_dummy) == 0)))
                                                ipface_dummy2 = pack(ipface_dummy,f%label(ipface_dummy) == 0)
                                                ipface = ipface_dummy2(1)
                                                deallocate(ipface_dummy2)

                                        end if
                                    else if (nipface .eq. 1) then

                                        ipface = ipface_dummy(1)

                                    end if

                                    ! Add face to tube_fc
                                    tube_count_fc = tube_count_fc + 1
                                    tube_fc(tube_count_fc) = ipface

                                else if (nipface == 0) then

                                    print *, 'Warning: BuildOpenTubesUS: No ipface found, put boundary face in tube'
                                    faceB = pack(fcs,f%label(fcs) /= 0)
                                    ipface = faceB(1)

                                end if

                                ! Housekeeping
                                deallocate(ipface_dummy)
                                deallocate(pfaces)

                            end do ! Marching algo

                        end if

                        ! Save

                        tube_count = tube_count + 1
                        ftCvP(tube_count,1) = cv_tot_count + 1
                        ftCvP(tube_count,2) = tube_cv
                        ftCv(cv_tot_count+1:cv_tot_count+tube_cv) = tube(1:tube_cv)
                        cv_tot_count = cv_tot_count + tube_cv

                        ftFcP(tube_count,1) = fc_tot_count + 1
                        ftFcP(tube_count,2) = tube_count_fc
                        ftFc(fc_tot_count+1:fc_tot_count+tube_count_fc) = tube_fc(1:tube_count_fc)

                        in_ftCv(tube(1:tube_cv)) = .true.
                        tube = 0
                        tube_fc = 0

                    end if

                end if

            end if

        end do

    
        end associate
    end subroutine

    subroutine BuildClosedTubes(grid, options, ft_open_us, ftCv, ftCvP, ftFc, ftFcP)

        ! Description
        !============
        ! Builds the flux tubes in the core by tracing the tube

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT) :: grid
        type(GAoptionsUDT) :: options
        logical :: ft_open_us
        integer(I8), allocatable, intent(inout) :: ftCv(:), ftCvP(:,:), ftFc(:), ftFcP(:,:)

        ! Auxiliary
        integer(I8) :: i, j, nc, cv, tube_count, nbv, counterc, counter_tube, &
            fc1, ifs, iv, cv1
        integer(I8), allocatable, dimension(:) :: cvLookUp, fcLbl_loc, bf_core, bv_core, &
            bc_core, cvs, fs_closed, tube, indFc, tube1, bc_core1
        logical, allocatable :: in_tube(:)
        logical :: core_in_ftCv

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert, &
            fd => grid%data%fluxdata &
            )
        
        ! Initialize
        allocate(fcLbl_loc(f%ntot))
        allocate(bc_core(c%ntot))
        allocate(in_tube(c%ntot))
        allocate(fs_closed(fd%nFs))
        allocate(tube(c%ntot))
        tube_count = count(ftCvP(:,1) /= 0)
        tube = 0
        cvLookUp = GetCvLookUp(c)
        bc_core = 0
        in_tube = .false.

        ! Face labelling
        fcLbl_loc = GetfcLbl(f, options)

        ! Build the core tube
        if ((tube_count == 0) .or. (ft_open_us)) then

            ! Build the core tube
            indFc = (/ (i, i = 1, f%ntot)/)
            bf_core = pack(indFc, fcLbl_loc == 2) ! Boundary faces of the core
            bv_core = GetVxsFromFcs(f, bf_core)
            nbv = size(bv_core)
            counterc = 0

            do i = 1, nbv

                iv = bv_core(i)
                cvs = GetVertCell(v, iv)
                nc = size(cvs)

                do j = 1, nc

                    cv = cvs(j) 
                    if (.not.in_tube(cv)) then 

                        counterc = counterc + 1
                        bc_core(counterc) = cv
                        in_tube(cv) = .true.

                    end if

                end do

            end do


        end if

        ! Trim
        core_in_ftCv = any(bc_core(1) == ftCv)

        ! Save the tube
        if (tube_count == 0) then

            ftCvP(1,1) = 1
            ftCvP(1,2) = counterc
            tube_count = 1
            in_tube(bc_core(1:counterc)) = .true.

        else if ((ft_open_us).and.(.not.core_in_ftCv)) then

            tube_count = tube_count + 1
            ftCvP(tube_count,1) = ftCvP(tube_count-1,1) + ftCvP(tube_count-1,2)
            ftCvP(tube_count,2) = counterc
            ftCv(ftCvP(tube_count,1):ftCvP(tube_count,1)+counterc-1) = bc_core(1:counterc)
            in_tube(bc_core(1:counterc)) = .true.

            ! Also add faces in ftFc and ftFcP
            bc_core1 = bc_core(1:counterc)
            call GiveFtFc(f, c, v, ftFc, ftFcP, bc_core1, tube_count)

        end if
        
        ! Initiate ftFc and ftFcP for output of function
        if (.not.ft_open_us) then

            allocate(ftFc(f%ntot))
            allocate(ftFcP(f%ntot,2))
            ftFc = 0
            ftFcP = 0

        end if

        ! First determine the start and end cells of the tubes
        ! INFORMATION SAVED SHOULD BE CORRECT BEFORE!!!! CHECK THIS - TODO
        !call GiveXpoints(grid, .false., cvLookUp)
        !call GiveSeparatrices(grid, .false., .false., .true., cvLookUp)

        ! Recycle the in_tube array
        ! Try without requiring a core cut
        ! But base on the flux surface
        fs_closed = 0

        do ifs = 1, fd%nFs

            if ((fd%fluxsurfacevertsP(ifs,2) == fd%fluxsurfacefacesP(ifs,2)) &
                .and. (.not.any(ifs == grid%data%sepID))) then
                fs_closed(ifs) = 1

                ! Build the tube
                fc1 = fd%fluxsurfacefaces(fd%fluxsurfacefacesP(ifs, 1))
                cvs = GetFaceCell(f, fc1)

                do j = 1, size(cvs)
                    cv1 = cvs(j)

                    if (.not.in_tube(cv1)) then

                        ! Start a new tube
                        call TraceCloseFluxTube(grid, in_tube, cv1, tube, counter_tube)

                        ! Add tube to ftCv and ftCvP
                        tube_count = tube_count + 1
                        ftCvP(tube_count,1) = ftCvP(tube_count-1,1) + ftCvP(tube_count-1,2)
                        ftCvP(tube_count,2) = counter_tube
                        ftCv(ftCvP(tube_count,1):ftCvP(tube_count,1)+ftCvP(tube_count,2)-1) = tube(1:counter_tube)

                        ! Also add ftFc if ft_open_us through intersection of mean psi (should work for the core)
                        tube1 = tube(1:counter_tube)
                        if (ft_open_us) call GiveFtFc(f, c, v, ftFc, ftFcP, tube1, tube_count)                     

                    end if

                end do

            end if
        end do 

        end associate


    end subroutine

    ! Auxiliary subroutine for BuildClosedTubes
    subroutine GiveFtFc(f, c, v, ftFc, ftFcP, tube, tube_count)

        ! Description
        !============
        ! Add poloidal faces of a tube to ftFc and ftFcP

        ! Declare variables
        !==================
        ! Arguments
        type(FaceUDT), intent(in) :: f
        type(CellUDT), intent(in) :: c
        type(VertexUDT), intent(in) :: v
        integer(I8), allocatable, intent(inout) :: ftFc(:), ftFcP(:,:), tube(:)
        integer(I8), intent(in) :: tube_count

        ! Auxiliary
        integer(I8), allocatable :: fc_tube(:), cells(:), tf(:), tv(:), &
            pfaces(:), fc_tube_clean(:)
        integer(I8) :: fc_count, k, j, ic, np, ifc
        real(R8), allocatable :: psi_vert(:)
        real(R8) :: psic, v1p, v2p
        
        allocate(fc_tube(f%ntot)) 
        fc_tube = 0
        fc_count = 0
        cells = tube

        do k = 1, size(cells)

            ! Get faces and verts
            ic = cells(k)
            tf = GetCellFace(c, ic)
            tv = GetCellVert(c, ic)

            ! Mean psi
            psi_vert = v%psi(tv)
            psic = 0.5_R8 * (maxval(psi_vert) + minval(psi_vert))

            ! Get poloidal faces
            allocate(pfaces(count(f%aligned(tf) == 0)))
            pfaces = pack(tf, f%aligned(tf) == 0)
            np = size(pfaces)

            ! Loop over poloidal faces
            do j = 1, np

                ! Add the face if it gets intersection by the mean psi of the cell
                ifc = pfaces(j)
                v1p = v%psi(f%vert(ifc,1))
                v2p = v%psi(f%vert(ifc,2))
                if ((psic .gt. min(v1p, v2p)) .and. (psic .lt. max(v1p, v2p)))  then

                    fc_count = fc_count + 1
                    fc_tube(fc_count) = ifc

                end if

            end do

            deallocate(pfaces)

        end do

        ! Remove duplicates
        call Unique(fc_tube(1:fc_count),fc_tube_clean)
        
        ! Add tube to ftFc
        ftFcP(tube_count,1) = ftFcP(tube_count-1,1) + ftFcP(tube_count-1,2)
        ftFcP(tube_count,2) = size(fc_tube_clean)
        ftFc(ftFcP(tube_count,1):ftFcP(tube_count,1)+ftFcP(tube_count,2)-1) = fc_tube_clean

    end subroutine

    subroutine TraceCloseFluxTube(grid, in_tube, cv1, tube, counter)

        ! Description
        !============
        ! Trace a close flux tube
        ! cv1 = starting cell
        ! in_tube = array indicating whether a cell is already in a tube 
        
        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(in) :: grid
        logical, allocatable, intent(inout) :: in_tube(:)
        integer(I8), intent(in) :: cv1
        integer(I8), intent(out) :: counter
        integer(I8), allocatable, intent(out) :: tube(:)

        ! Auxiliary
        logical :: looped
        integer(I8), allocatable :: fcs(:), p_fcs(:), cvs(:)
        integer(I8) :: ic, opp_face
        character(:), allocatable :: meth

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        ! Initialize
        allocate(tube(c%ntot))
        tube = 0
        tube(1) = cv1
        counter = 1
        in_tube = .false.
        in_tube(cv1) = .true.
        looped = .false.

        ! Find common face
        fcs = GetCellFace(c, cv1)
        p_fcs = pack(fcs, f%aligned(fcs) == 0)
        if (size(p_fcs) /= 2) call gdErrorHandler('TraceClosedFluxTube: Something went wrong, cell has no two poloidal faces')

        ! Initialize while loop
        ic = cv1
        opp_face = p_fcs(1)

        ! Start the while loop
        do while (.not.looped)

            ! Find the oppostive face to the common face
            meth = 'pol'
            opp_face = GetOppositeFace(opp_face,ic, grid, meth)

            ! Add the second cell
            cvs = GetFaceCell(f, opp_face)

            if (size(cvs) /= 2) call gdErrorHandler('TraceClosedFluxTube: poloidal face in the core has no two cells')

            if ((cvs(1) /= cv1) .and. (.not.in_tube(cvs(1)))) then

                ! Add cvs(1)
                ic = cvs(1)
                counter = counter + 1
                tube(counter) = ic
                in_tube(ic) = .true.

            else if ((cvs(2) /= cv1) .and. (.not.in_tube(cvs(2)))) then

                ! Add cvs(2)
                ic = cvs(2)
                counter = counter + 1
                tube(counter) = ic
                in_tube(ic) = .true.
            
            else

                ! The tube is closed
                looped = .true.

            end if
            
        end do

        
        end associate


    end subroutine

    !subroutine AddSeparatrixTube(grid, ftCv, ftCvP, ftFc, ftFcP, first_core_tube) 

        ! Description
        !============
        ! Add the separatrix flux tube based on fsFc only where it is connected with
        ! the core region

        ! Declare variables
        !==================
        ! Arguments
        !type(GridUDT) :: grid
        !integer(I8), allocatable :: ftCv(:), ftCvP(:,:), ftFc(:), ftFcP(:,:)
        !integer(I8) :: first_core_tube(:)

        ! Auxiliary
        !logical, allocatable :: in_first_core_tube

    !end subroutine


    subroutine BuildFtFc(grid)

        ! Description
        !============
        ! Build the flux tube data ftFc based on ftCv. It contains for every flux tube
        ! the faces in the flux tube, being the poloidal faces. The poloidal faces don't
        ! need to have specific alignment to be part of ftFc. For triangles at the target boundary.

        ! Arguments
        !==========
        type(GridUDT), intent(inout)              :: grid

        ! Auxiliary
        integer(I8), allocatable :: fc_tube(:), cells(:), faces(:),  &
            fc_tube2(:), verts(:), pfaces(:), ind(:), ftFcP(:,:), ftFc(:)
        integer(I8) :: fc_count, ift, i, j, ic, s, n, np, nc, ifc, &
            v1, v2
        real(R8), allocatable :: psi_vert(:)
        real(R8) :: psic, v1p, v2p

        ! Initialize
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert, &
            fd => grid%data%fluxdata &
            )

        allocate(fc_tube(f%ntot))
        fc_tube = 0
        fc_count = 0
        allocate(ftFc(f%ntot))
        ftFc = 0
        allocate(ftFcP(fd%nFt, 2))
        ftFcP = 0

        ! Loop over flux tubes
        do ift = 1, fd%nFt

            ! Get cells for this tube
            s = fd%fluxtubecellsP(iFt,1)
            nc = fd%fluxtubecellsP(iFt,2)
            allocate(cells(nc))
            cells = fd%fluxtubecells(s:s+nc-1)

            ! Loop over cells in this tube
            do i = 1, nc
                ic = cells(i)

                n = c%faceP(ic,2)
                allocate(faces(n))
                allocate(verts(n))

                faces = GetCellFace(c, ic)
                verts = GetCellVert(c, ic)
                

                ! Mean psi
                allocate(psi_vert(n))
                psi_vert = v%psi(verts)
                psic = 0.5_R8 * (maxval(psi_vert) + minval(psi_vert))

                ! Get poloidal faces (those not aligned)
                np = count(f%aligned(faces) == 0)
                allocate(pfaces(np))
                pfaces = pack(faces, f%aligned(faces) == 0)

                ! Loop over poloidal faces
                do j = 1, np
                    ifc = pfaces(j)

                    v1  = f%vert(ifc,1)
                    v2  = f%vert(ifc,2)
                    v1p = v%psi(v1)
                    v2p = v%psi(v2)

                    if ( (psic .gt. min(v1p,v2p)) .and. (psic .gt. max(v1p,v2p)) ) then
                        fc_count = fc_count + 1
                        fc_tube(fc_count) = ifc
                    end if
                end do

                deallocate(faces, verts, psi_vert, pfaces)

            end do  ! cell loop

            ! Unique the fc_tube list
            if (fc_count .gt. 0) then
                call Unique(fc_tube(1:fc_count),fc_tube2)
            end if

            ! Store results in ftFc and ftFcP
            if (ift == 1) then
                ftFcP(ift,1) = 1
            else 
                ftFcP(ift,1) = ftFcP(ift-1,1) + ftFcP(ift-1,2) 
            end if
            ftFcP(ift,2) = size(fc_tube2)
            ind = (/(i, i= ftFcP(ift,1),ftFcP(ift,1)+ftFcP(ift,2)-1)/)
            ftFc(ind) = fc_tube2 

            !call AppendArray(ftFc, fc_tube(1:fc_count))

            ! Reset for next tube
            fc_tube = 0
            fc_count = 0
            deallocate(fc_tube2)
            deallocate(cells)

        end do

        ! Store
        fd%fluxtubefacesP = ftFcP
        fd%fluxtubefaces  = ftFc

        end associate

    end subroutine BuildFtFc


    subroutine DetermineMPs(grid, options)

        ! Description
        !============
        ! Determine the outer mid plane (OMP) and inner mid plane (IMP)

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(inout)    :: grid
        type(GAoptionsUDT), intent(inout)  :: options        

        ! Auxiliary
        integer(I8) :: mode
        integer(I8), allocatable, dimension(:) :: cvOMP, fcOMP, cvIMP, fcIMP

        ! OMP
        mode = 1
        call CvIntersectionsGOAT(options%OMP_r,options%OMP_z,grid,mode,cvOMP,fcOMP)
        
        ! IMP
        mode = 1
        call CvIntersectionsGOAT(options%IMP_r,options%IMP_z,grid,mode,cvIMP,fcIMP)

        ! Save
        grid%data%OMPcell   = cvOMP
        grid%data%nOMPcell  = size(cvOMP)
        grid%data%OMPface   = fcOMP
        grid%data%nOMPface  = size(fcOMP)
        grid%data%IMPcell   = cvIMP
        grid%data%nIMPcell  = size(cvIMP)
        grid%data%IMPface   = fcIMP
        grid%data%nIMPface  = size(fcIMP)

    end subroutine

    subroutine CvIntersectionsGOAT(segm_r, segm_z, grid, mode, listcv, listfc)

        ! Description
        !============
        ! Determine the cells intersected by a segment. 
        ! The r or x-coordinates of the segment are defined in segm_r
        ! The z or y-coordinates of the segment are defined in segm_z

        ! Declare variables
        ! =================
        ! Arguments
        real(R8), intent(inout) :: segm_r(1:2), segm_z(1:2)
        type(GridUDT) :: grid
        integer(I8) :: mode
        integer(I8), allocatable :: listcv(:), listfc(:)

        ! Auxiliary
        integer(I8) :: nn, ii, kk, cv, ifc, nc, nomp
        integer(I8), allocatable :: cvs(:), &
            listcv_dummy(:), listfc_dummy(:)
        logical :: no_intersections, found
        real(R8) :: m, a, b, cc
        real(R8), allocatable, dimension(:) :: dist_vert, p1, &
            q1, p2, q2


        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        ! Initialize
        allocate(listcv_dummy(200))
        allocate(listfc_dummy(200))
        allocate(dist_vert(v%ntot))
        nn = 1

        ! Check whether the segment intersects any vertices, if so shift the
        ! segment
        ! a line : a*x+ b*y + c = 0
        ! Combined with y-y_a = m (x - x_a)
        ! with m = (y_b - y_a)/(x_b - x_a)
        ! Gives : a = m, b = -1, c = y_a - m x_a   
        no_intersections = .false.

        do while (.not.no_intersections)

            ! Parameterize line
            m = (segm_z(2) - segm_z(1))/(segm_r(2) - segm_r(1))
            a = m 
            b = -1.0_R8
            cc = segm_z(1) - m*segm_r(1)

            ! Distance to vertices
            dist_vert = abs(a*v%x + b*v%y + cc)/sqrt(a**2 + b**2)

            if (minval(dist_vert) .lt. 1e-10_R8) then

                ! Adjust the segment and test again 
                segm_r = segm_r + 1e-6_R8
                segm_z = segm_z + 1e-6_R8
                no_intersections = .false.

            else

                ! Continue with current segment
                no_intersections = .true.

            end if

        end do

        ! Segment
        p1 = [segm_r(1), segm_z(1)]
        q1 = [segm_r(2), segm_z(2)]

        ! Loop over all faces
        do ifc = 1, f%ntot

            p2 = [ v%x(f%vert(ifc,1)), v%y(f%vert(ifc,1)) ]
            q2 = [ v%x(f%vert(ifc,2)), v%y(f%vert(ifc,2)) ]
            
            if (Intersects(p1,q1,p2,q2) == 1) then

                cvs = GetFaceCell(f, ifc)
                nc = size(cvs)
                do ii = 1, nc
                    cv = cvs(ii)
                    found = .false.
                    kk = 1
                    do while (.not.found .and. kk .le. nn)

                        if (listcv(kk) == cv) then 
                            found = .true.
                        end if
                        kk = kk + 1
                    end do
                    if (.not.found) then

                        listcv_dummy(nn) = cv
                        listfc_dummy(nn) = ifc
                        nn = nn + 1
                    end if

                end do

            end if

        end do

        ! Trim
        nn = nn - 1
        listcv = listcv_dummy(1:nn)
        listfc = listfc_dummy(1:nn)

        ! Sorting
        nomp = size(listcv)

        ! Sorting - direction does not matter for now
        call SortOmpListGOAT(grid,nomp,listcv,mode,segm_r,segm_z);



        end associate

    end subroutine 

    subroutine SortOmpListGOAT(grid,nomp,listcv,mode,segm_r,segm_z)

        ! Description
        !============
        ! Sort the cells in the intersection cell list

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(in) :: grid
        integer(I8), intent(in) :: nomp, mode
        integer(I8), allocatable, intent(inout) :: listcv(:)
        real(R8), intent(in) :: segm_r(2), segm_z(2)

        ! Auxiliary
        integer(I8) :: i, j, ic, nc, l, ic1, iv, ifc
        integer(I8), allocatable :: omp_cv_sorted(:), omp_cv(:), cvs(:), cvLookUp(:)
        logical :: found, in_list, in_list_sorted

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        ! Initialize
        cvLookUp = GetCvLookUp(c)
        omp_cv = listcv

        ! Find a boundary cell to start
        allocate(omp_cv_sorted(nomp))
        omp_cv_sorted = 0

        do i = nomp, 1, -1
            ic = omp_cv(i)

            if (omp_cv_sorted(1) == 0) then
                omp_cv_sorted(1) = ic
            elseif (mode == 1 .and. c%x(ic) .lt. c%x(omp_cv_sorted(1))) then
                omp_cv_sorted(1) = ic
            elseif (mode == 2 .and. c%y(ic) .lt. c%y(omp_cv_sorted(1))) then
                omp_cv_sorted(1) = ic
            end if

        end do

        ! Sorting itself
        ! First via poloidal faces
        found = .true.
        do i = 1, nomp-1

            ic = omp_cv_sorted(i)
            do j = c%faceP(ic,2), 1, -1

                ifc = c%face(c%faceP(ic,1)+j-1)
                if (f%aligned(ifc) == 0) then

                    cvs = GetFaceCell(f, ifc)
                    nc = size(cvs)
                    do l = 1, nc

                        in_list = .false.
                        in_list_sorted = .false.
                        ic1 = cvs(l)

                        if (ic1 /= ic .and. omp_cv_sorted(i+1) == 0) then

                            ! Check whether ic1 is in omp_cv and not in omp_cv_sorted
                            if (any(omp_cv == ic1) .and. .not. any(omp_cv_sorted(1:i) == ic1)) &
                                omp_cv_sorted(i+1) = ic1
                 
                        end if

                    end do

                end if
            end do
            
            ! Then through aligned faces
            if (omp_cv_sorted(i+1) == 0 ) then

                do j = c%faceP(ic,2), 1, -1

                    ifc = c%face(c%faceP(ic, 1)+j-1)
                    if (f%aligned(ifc) == 1) then

                        cvs = GetFaceCell(f, ifc)
                        nc = size(cvs)

                        do l = 1, nc

                            ic1 = cvs(l)
                            if (ic1 /= ic .and. omp_cv_sorted(i+1) == 0) then
                                if (any(omp_cv == ic1) .and. .not. any(omp_cv_sorted(1:i) == ic1)) &
                                omp_cv_sorted(i+1) = ic1

                            end if
                        end do
                    end if
                end do

                
            end if 

            ! if no next cell is found, allow to go over vertices
            ! if no next cell is found through face neigbours, look at vertex
            ! neighbours   
            if (omp_cv_sorted(i+1) == 0) then

                do j = 1, c%vertP(ic, 2)

                    iv = c%vert(c%vertP(ic,1)+j-1)
                    cvs = GetVertCell(v, i)
                    nc = size(cvs)

                    do l = 1, nc

                        ic = cvs(l)
                        if (ic1 /= ic .and. omp_cv_sorted(i+1) == 0) then
                            if (any(omp_cv == ic1) .and. .not. any(omp_cv_sorted(1:i) == ic1)) &
                            omp_cv_sorted(i+1) = ic1

                        end if                        
                    end do

                end do

            end if

            ! If there would be no next cell - error => change segment
            if (omp_cv_sorted(i+1) == 0) then
                found = .false.;
                exit
            end if          

        end do

        if (.not.found) then
            print *, "Warning: CvIntersectionsGOAT: possibly multiple grid intersections found, // &
             & check the position of the MP segments"
            omp_cv_sorted = 0
        end if

        ! Save
        listcv = omp_cv_sorted

        end associate

    end subroutine


    !==================================================================!
    !                                                                  !
    !                           FUNCTIONS                              !
    !                                                                  !
    !==================================================================!     

    function GetCvLookUp(cell) result(res)
        type(CellUDT)       :: cell
        integer(I8)         :: nc, ic, nv, s, i               
        integer(I8), allocatable :: res(:), range(:)

        nc = cell%ntot
        range = (/ (i, i = 1,(cell%vertP(nc,1)+cell%vertP(nc,2)-1))/)

        allocate(res(1:cell%vertP(nc,1)+cell%vertP(nc,2)-1))
        res = 0

        do ic = 1, nc
            s = cell%vertP(ic,1)
            nv = cell%vertP(ic,2)
            range = (/ (i, i = s, (s+nv-1)) /)
            res(range) = ic
        end do
    end function    

    function GetfcLbl(f,options) result(res)
        type(FaceUDT) :: f
        type(GAoptionsUDT) :: options
        integer(I8) :: res(1:f%ntot), indFc(1:f%ntot), i 

        res = f%label
        indFc = (/ (i, i=1,f%ntot) /)
        do i = 1, size(options%facelabelmappingGG)
            res(pack(indFc, f%label == options%facelabelmappingGG(i))) &
                = options%facelabelmappingGA(i)
        end do

    end function

    function GetVxsFromFcs(f,fcs) result(res)
        type(FaceUDT) :: f
        integer(I8), allocatable :: fcs(:), verts(:), res(:)
        integer(I8) :: nf

        nf = size(fcs)
        allocate(verts(1:nf*2))
        verts = 0
        verts(1:nf) = f%vert(fcs,1)
        verts(nf+1:nf*2) = f%vert(fcs,2)
        call Unique(verts, res)

    end function    

    function GetOppositeFace(ifc, ic, grid, meth) result(res)

        ! Description
        !============
        ! Give a face opposite to a given face in an aligned quad or pent.
        ! Assummes sequential ordening of the cell faces   
        
        ! Declare variables
        ! =================
        ! Arguments
        type(GridUDT) :: grid
        integer(I8) :: ifc, ic, res
        character(:), allocatable :: meth

        ! Auxiliary
        integer(I8), allocatable :: fcs(:), fcs1(:), fcs2(:), &
            indf(:), ind(:)
        integer(I8) :: vsF(1:2), i, ind_com_face, ind_opp_face, &
            vs(1:2), nf
        real(R8) :: psic

        res = 0
        fcs = GetCellFace(grid%cell, ic)

        select case (meth)

        case ('pol')
            
            vsF = grid%face%vert(ifc,:)
            psic = 0.5_R8*sum(grid%vert%psi(vsF))
            allocate(fcs1(count(fcs /= ifc)))
            fcs1 = pack(fcs, fcs /= ifc)
            allocate(fcs2(count(grid%face%aligned(fcs1) == 0)))
            fcs2 = pack(fcs1, grid%face%aligned(fcs1) == 0)
            do i = 1, size(fcs2)
                vs = grid%face%vert(fcs2(i),:)
                if ((psic .gt. minval(grid%vert%psi(vs))) &
                    .and. psic .lt. maxval(grid%vert%psi(vs))) then
                        res = fcs2(i)
                        exit
                end if

            end do

            if (res == 0) res = fcs2(1)

        case ('ordened_quad')

            nf = size(fcs)
            if (nf == 4) then
                indf = (/ (i, i = 1, nf)/)
                ind = pack(indf, fcs == ifc)
                ind_com_face = ind(1)
                ind_opp_face = ind_com_face + 2
                if (ind_opp_face .gt. nf) ind_opp_face = ind_opp_face - nf 
                res = fcs(ind_opp_face)

            end if
        case default

            call gdErrorHandler("GetOppositeFace: not yet implemented")

        end select

        if (res == 0) call gdErrorHandler('GetOppositeFace: no opposite face found')

    end function

    function Intersects(p1, q1, p2, q2) result(res)
        real(R8) :: p1(2), q1(2), p2(2), q2(2)
        integer(I8) :: res, o1, o2, o3, o4


        o1 = Orient(p1, q1, p2)
        o2 = Orient(p1, q1, q2)
        o3 = Orient(p2, q2, p1)
        o4 = Orient(p2, q2, q2)

        res = 1

        if (o1 /= o2 .and. o3 /= o4) then
            return
        end if

        !c1 = 
        if (o1 == 0 .and. OnSegment(p1, p2, q1) == 1) then 
            return
        end if
    
        if (o2 == 0 .and. OnSegment(p1, q2, q1) == 1) then
            return
        end if
    
        if (o3 == 0 .and. OnSegment(p2, p1, q2) == 1) then
            return 
        end if
    
        if (o4 == 0 .and. OnSegment(p2, q1, q2) == 1) then
            return 
        end if

        res = 0;

    end function

    function Orient(p, q, r)  result(res)
        real(R8) :: p(2), q(2), r(2), val
        integer(I8) :: res

        val = (q(2) - p(2)) * (r(1) - q(1)) - (q(1) - p(1)) * (r(2) - q(2))

        if (val > 0) then

            res = 1
            return

        else

            if (val < 0) then
                res = 2
                return
            else
                res = 0
                return
            end if

        end if

    end function

    function OnSegment(p, q, r) result(res)
        real(R8) :: p(2), q(2), r(2)
        integer(I8) :: res 

        res = 0
        if ( (q(1) <= max(p(1), r(1))) .and. (q(1) >= min(p(1), r(1))) .and. &
             (q(2) <= max(p(2), r(2))) .and.  (q(2) >= min(p(2), r(2))) ) then
            res = 1
        end if
        
        return               
    end function

    function isBoundaryFace0D(ifc, f) result(res)
        integer(I8) :: ifc
        type(FaceUDT) :: f
        logical :: res

        res = .true. 
        if (f%label(ifc) == 0) then
            res = .false.
        end if
    end function

    function isBoundaryFace1D(fc, f) result(res)
        integer(I8), allocatable :: fc(:)
        integer(I8) :: i
        type(FaceUDT) :: f
        logical, allocatable :: res(:)

        allocate(res(size(fc)))
        res = .true. 
        do i = 1, size(fc)
            if (f%label(fc(i)) == 0) then
                res(i) = .false.
            end if
        end do

    end function

    function GetFluxTubeFromCellIndex(fd, indc) result(res)
        type(FluxDataUDT) :: fd
        integer(I8) :: i, indc, n_el, res
        integer(I8), allocatable :: b(:), ind(:)

        n_el = count(fd%fluxtubecellsP(:,1).le.indc)
        allocate(b(n_el))
        ind = (/ (i, i = 1, n_el)/)
        b = pack(ind, fd%fluxtubecellsP(:,1).le.indc)
        res = b(n_el)
    end function


end module 