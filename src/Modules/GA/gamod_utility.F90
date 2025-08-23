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
           first_core_tube(:), ftCv(:), ftCvP(:,:), ftFc(:), ftFcP(:,:)
        integer(I8) :: nc, c1, c2, sv, nv, sf, nf, i, ift
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
        if (ft_open) then

            ! Build the open tubes connecting the target for structured grids   
            call BuildOpenTubes(grid, options, ftCv, ftCvP)

        else if (ft_open_us) then

            ! Build the open tubes connecting the target for unstructured grids
            call BuildOpenTubesUS(grid, options, ftCv, ftCvP, ftFc, ftFcP)

        end if

        ! Build closed tubes
        if (ft_open .and. ft_closed) then

            ! Build the core flux tubes while open flux tubes were build
            call BuildClosedTubes(grid, options, ft_open_us)

        elseif (ft_open_us .and. ft_closed) then

            ! Build the core flux tubes while open flux tubes were build and add ftFc
            call BuildClosedTubes(grid,options,ft_open_us)

        elseif (.not.(ft_open .or. ft_open_us) .and. ft_closed) then 

            call BuildClosedTubes(grid,options,ft_open_us, first_core_tube)
            call AddSeparatrixTube(grid, first_core_tube)            

        end if

        ! Build ftFc - poloidal faces which lay into the flux tube
        if (.not.ft_open_us) then
            call BuildFtFc(grid)
        end if

        ! Store
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert, &
            fd => grid%data%fluxdata &
        )

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
                    ind = (/ (i, i = nv,-1,1) /)
                    cvs_rev = cvs(ind)

                    sv = fd%fluxtubecellsP(ift,1)
                    nv = fd%fluxtubecellsP(ift,2)  
                    
                    ind = (/ (i, i = sv, sv+nv-1)/)
                    fd%fluxtubecells(ind) = cvs_rev
                    deallocate(ind)

                    ! Faces
                    sf = fd%fluxtubefacesP(ift,1)
                    nf = fd%fluxtubefacesP(ift,2)

                    fcs = fd%fluxtubefaces(sf:sf+nf-1)

                    fcs_rev = fcs(nf:-1:1)
                    fd%fluxtubefaces(sf:sf+nf-1) = fcs_rev

                end if
            end if 
        end do

        ! Make cvFt -- inverse of ftCv
    
        ! Make ftReg

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
        allocate(ftCvP(c%ntot,2)) ! fd%nFt not known !!!!! TODO
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
        integer(I8), allocatable :: tube(:)
        integer(I8) :: tube_count


        ! Associate
        associate(& 
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert, &
            fd => grid%data%fluxdata &
            )

        allocate(ftCv(c%ntot))
        allocate(ftCvP(fd%nFt,2))
        allocate(tube(c%ntot))
        ftCv = 0
        ftCvP = 0
        tube = 0
        tube_count = 0


    
        end associate
    end subroutine

    subroutine BuildClosedTubes(grid, options, ft_open_us, first_core_tube)

        ! Description
        !============
        ! Builds the flux tubes in the core by tracing the tube

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT) :: grid
        type(GAoptionsUDT) :: options
        logical :: ft_open_us
        integer(I8), allocatable, intent(out), optional :: first_core_tube(:)

        ! Auxiliary

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert, &
            fd => grid%data%fluxdata &
            )


        end associate


    end subroutine

    subroutine AddSeparatrixTube(grid, first_core_tube) 

        ! Description
        !============
        ! Add the separatrix flux tube based on fsFc only where it is connected with
        ! the core region

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT) :: grid
        integer(I8), allocatable :: first_core_tube(:)

    end subroutine


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
        type(GAoptionsUDT), intent(in)  :: options        

        ! Auxiliary

    end subroutine


    !==================================================================!
    !                                                                  !
    !                           FUNCTIONS                              !
    !                                                                  !
    !==================================================================!     



end module 