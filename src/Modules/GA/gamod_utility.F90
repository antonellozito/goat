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
        class(QualityMetric)    :: qm
        type(GridUDT)           :: grid
        type(GAoptionsUDT)      :: options

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
        logical, intent(in)          :: is_ordered(grid%cell%ntot)

        ! Declare variables
        !==================
        integer(I8) :: ic, nv, fcs(1:2), vs(1:2), s, n, i
        integer(I8), allocatable, dimension(:) :: tv, tf, ntv, ntf
        integer(I8), allocatable, dimension(:,:) :: vf
        real(R8) :: vec_start(1:2), vec_face1(1:2), vec_face2(1:2), &
            sin1, sin2, fcX(grid%face%ntot), fcY(grid%face%ntot)

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert  &
            )
        
        ! Compute face centers
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
        !fd%fluxsurfacevertices = fsVx(1:nv_counter)
        !fd%fluxsurfaceverticesP = fsVxP(1:fd%nFs,:)

        end associate





    end subroutine


end module 