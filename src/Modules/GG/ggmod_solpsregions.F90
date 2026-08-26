!======================================================================!
!                                                                      !
!                            DOCUMENTATION                             !
!                                                                      !
!======================================================================!

! Description
!============
! This module reworks the volume regions (cvReg) and face regions
! (fcReg) of a Goat grid into the region numbering that
! SOLPS-ITER expects, based on the catalog topology recognized by
! ClassifySOLPSCatalogTopology.
! One routine per catalog family. Currently implemented:
!  - SetSOLPSRegionsLIM: limiter topologies (NREG = 2, FREG = 6)
!  - SetSOLPSRegionsSN: lower/upper single null (NREG = 4, FREG = 13)

module ggmod_solpsregions

    ! Load modules
    use mod_precision
    use mod_errorhandler
    use mod_definitions
    use mod_constants, only: posinfval_R8
    use goatmod_types, only: GridUDT, VesselUDT, GetVertFace, &
        GetFaceCell, GetCellFace
    use ggmod_topology2D

    implicit none
    private
    public :: SetSOLPSRegionsLIM, SetSOLPSRegionsSN

contains

    ! SOLPS region conventions for the LIM (limiter) catalog topology
    subroutine SetSOLPSRegionsLIM(simgrid, topomesh, vessel)

        ! Description
        !============
        ! Assign the volume and face regions that SOLPS expects for a
        ! limiter case (NREG = 2, FREG = 6):
        !
        ! Volume regions (cvReg):
        !   1 = core (inside the separatrix)
        !   2 = SOL (everything else)
        !
        ! Face regions (fcReg):
        !   1 = left target: the part of the touched wall structure
        !       on the clockwise side (standard (R,z) view) of the
        !       point where the separatrix touches the wall.
        !   2 = right target: the counterclockwise part.
        !   3 = (empty) core connection.
        !   4 = core boundary (the inner grid boundary)
        !   5 = separatrix: internal aligned faces on the interface
        !       between core and SOL topomesh cells
        !   6 = the rest of the wall / outer boundary
        !
        ! The targets are bounded to the composite limiter surface
        ! formed by the SOLPS limiter target segments
        ! (goat.vessel.solps_limiter_target): the left target extends
        ! from the separatrix-wall touch point to one end of that
        ! composite surface, the right target to the other end. For a
        ! SOLPS limiter grid that list is REQUIRED to be non-empty
        ! and to form exactly one consecutive stretch of wall
        ! elements - otherwise this routine errors out (the same
        ! requirement is already enforced fail-fast, before any
        ! gridding, in InitializeGGTMLineRefiner).

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(inout)    :: simgrid
        type(TopomeshUDT), intent(in)   :: topomesh
        type(VesselUDT), intent(in)     :: vessel

        ! Auxiliary
        integer(I8), allocatable        :: cc(:), ccf(:), tp(:), &
            tcv(:), tfc(:), tf2(:)
        logical, allocatable            :: coremaskTM(:), &
            iscorefTM(:), issepTM(:), slmask(:)
        integer(I8)                     :: i, k, itp, istrtouch, lbl, &
            nleft, nright, nwall, ncoreb, nsep, nmix, iv, nrun, kp
        real(R8)                        :: xc, yc, xt, yt, xm, ym, &
            d, dmin, tht, thm, dth
        logical                         :: havetargets

        ! Volume regions
        !===============
        ! Core (marked with the SOLPS core region IDs during the
        ! label translation, or already 1/2 from the generic limiter
        ! remapping - the test is idempotent) -> 1, rest -> 2
        where (mod(simgrid%cell%reg - SOLPScoreregID, &
            SOLPScoreregIDincr) == 0)
            simgrid%cell%reg = 1
        elsewhere
            simgrid%cell%reg = 2
        end where

        ! Topomesh masks
        !===============
        allocate(coremaskTM(topomesh%cell%ntot), &
            iscorefTM(topomesh%face%ntot), &
            issepTM(topomesh%face%ntot))
        coremaskTM = .false.
        iscorefTM = .false.
        issepTM = .false.
        cc = topomesh%GetCoreCellIDs()
        coremaskTM(cc) = .true.
        ccf = topomesh%GetCoreFaceIDs()
        if (size(ccf) > 0) iscorefTM(ccf) = .true.

        ! Separatrix topomesh faces: interfaces between a core and a
        ! non-core topomesh cell
        do i = 1, topomesh%face%ntot
            tfc = topomesh%face%GetCell(i)
            if (size(tfc) == 2) then
                if (coremaskTM(tfc(1)) .neqv. coremaskTM(tfc(2))) then
                    issepTM(i) = .true.
                end if
            end if
        end do

        ! Core centroid (reference for the poloidal angle)
        !==================================================
        xc = 0.0_R8
        yc = 0.0_R8
        k = 0
        do i = 1, size(cc)
            tcv = topomesh%cell%GetVert(cc(i))
            xc = xc + sum(topomesh%vert%x(tcv))
            yc = yc + sum(topomesh%vert%y(tcv))
            k = k + size(tcv)
        end do
        xc = xc/real(k, kind=R8)
        yc = yc/real(k, kind=R8)

        ! Separatrix-wall touch point and touched wall structure
        !========================================================
        ! The touch point is the topomesh vertex shared by a
        ! separatrix (core/SOL interface) face and a wall (boundary)
        ! face. This is robust against the tangency-point typing and
        ! closed-contour checks of the topomesh construction.
        havetargets = .false.
        xt = 0.0_R8
        yt = 0.0_R8
        istrtouch = 0

        ! Requirement for SOLPS limiter grids: the SOLPS limiter
        ! target segments (goat.vessel.solps_limiter_target) must be
        ! given and must form exactly one consecutive stretch of wall
        ! elements (checked modulo the wall wrap-around)
        allocate(slmask(vessel%nstructures))
        slmask = .false.
        if (vessel%nsolpslim <= 0) then
            call gdErrorHandler('SetSOLPSRegionsLIM: ' // &
                'goat.vessel.solps_limiter_target is empty - a ' // &
                'SOLPS limiter grid requires the limiter ' // &
                'surface to be given there as one consecutive ' // &
                'stretch of wall elements')
        end if
        do i = 1, int(vessel%nsolpslim, kind=I8)
            k = int(vessel%solpslimind(i), kind=I8)
            if (slmask(k)) then
                call gdErrorHandler('SetSOLPSRegionsLIM: ' // &
                    'duplicate wall element in ' // &
                    'goat.vessel.solps_limiter_target')
            end if
            slmask(k) = .true.
        end do
        nrun = 0
        do k = 1, int(vessel%nstructures, kind=I8)
            kp = k - 1
            if (kp == 0) kp = int(vessel%nstructures, kind=I8)
            if (slmask(k) .and. .not. slmask(kp)) nrun = nrun + 1
        end do
        if ((nrun /= 1) .and. &
            (count(slmask) /= int(vessel%nstructures, kind=I8))) &
            then
            call gdErrorHandler('SetSOLPSRegionsLIM: the wall ' // &
                'elements in goat.vessel.solps_limiter_target ' // &
                'must form exactly one consecutive stretch ' // &
                '(the limiter surface) for a SOLPS limiter grid')
        end if

        allocate(tp(0))
        do i = 1, topomesh%face%ntot
            if (.not. issepTM(i)) cycle
            do k = 1, 2
                iv = topomesh%face%vert(i, k)
                if (size(tp) > 0) then
                    if (any(tp == iv)) cycle
                end if
                tf2 = topomesh%vert%GetFace(iv)
                if (any((topomesh%face%type(tf2) == TMfacebndID) &
                        .or. &
                        (topomesh%face%type(tf2) == TMfacealbndID))) &
                    then
                    tp = [tp, iv]
                end if
            end do
        end do
        if (size(tp) == 0) then
            print *, 'SetSOLPSRegionsLIM: no separatrix-wall ' // &
                'touch point found - leaving targets empty'
        else
            if (size(tp) > 1) then
                print *, 'SetSOLPSRegionsLIM: ', size(tp), &
                    ' separatrix-wall touch vertices found, ' // &
                    'using the first: ', tp
            end if
            itp = tp(1)
            xt = topomesh%vert%x(itp)
            yt = topomesh%vert%y(itp)

            ! Touched structure: nearest structure.dat sub-structure
            dmin = posinfval_R8()
            do k = 1, int(vessel%nstructures, kind=I8)
                d = DistToStructure(xt, yt, k)
                if (d < dmin) then
                    dmin = d
                    istrtouch = k
                end if
            end do
            havetargets = .true.
            print *, 'SetSOLPSRegionsLIM: separatrix touch point ', &
                xt, yt, ' on structure ', istrtouch
            if (.not. slmask(istrtouch)) then
                print *, 'SetSOLPSRegionsLIM: WARNING: the touch ' // &
                    'point does not lie on the limiter surface ' // &
                    '(goat.vessel.solps_limiter_target) - one of ' // &
                    'the two targets will be empty'
            end if
        end if

        ! Face regions
        !=============
        tht = atan2(yt - yc, xt - xc)
        simgrid%face%reg = 0
        nleft = 0
        nright = 0
        nwall = 0
        ncoreb = 0
        nsep = 0
        do i = 1, simgrid%face%ntot
            lbl = simgrid%face%TMfacelabel(i)
            if (simgrid%face%BF(i)) then

                ! Core boundary
                if ((lbl >= 1) .and. (lbl <= topomesh%face%ntot)) then
                    if (iscorefTM(lbl)) then
                        simgrid%face%reg(i) = 4
                        ncoreb = ncoreb + 1
                        cycle
                    end if
                end if

                ! Wall: split the touched structure into left/right
                ! targets, everything else is region 6
                simgrid%face%reg(i) = 6
                nwall = nwall + 1
                if (havetargets) then
                    xm = 0.5_R8*(simgrid%vert%x(simgrid%face%vert(i, 1)) &
                        + simgrid%vert%x(simgrid%face%vert(i, 2)))
                    ym = 0.5_R8*(simgrid%vert%y(simgrid%face%vert(i, 1)) &
                        + simgrid%vert%y(simgrid%face%vert(i, 2)))
                    dmin = posinfval_R8()
                    k = 0
                    do lbl = 1, int(vessel%nstructures, kind=I8)
                        d = DistToStructure(xm, ym, lbl)
                        if (d < dmin) then
                            dmin = d
                            k = lbl
                        end if
                    end do
                    if (slmask(k)) then
                        ! Face on the composite limiter surface:
                        ! clockwise side of the touch point = left
                        ! target (region 1), counterclockwise =
                        ! right target (region 2)
                        thm = atan2(ym - yc, xm - xc)
                        dth = atan2(sin(thm - tht), cos(thm - tht))
                        nwall = nwall - 1
                        if (dth < 0.0_R8) then
                            simgrid%face%reg(i) = 1
                            nleft = nleft + 1
                        else
                            simgrid%face%reg(i) = 2
                            nright = nright + 1
                        end if
                    end if
                end if

            else

                ! Internal faces: separatrix -> 5, rest stays 0
                if ((lbl >= 1) .and. (lbl <= topomesh%face%ntot)) then
                    if (issepTM(lbl)) then
                        simgrid%face%reg(i) = 5
                        nsep = nsep + 1
                    end if
                end if

            end if
        end do

        print *, 'SetSOLPSRegionsLIM: fcReg counts: left target ', &
            nleft, ', right target ', nright, ', core boundary ', &
            ncoreb, ', separatrix ', nsep, ', wall/outer ', nwall

        ! Flux tube regions
        !==================
        ! In a limiter topology every flux tube lies entirely in the
        ! core or entirely in the SOL: inherit the (rewritten) cell
        ! region -> 1 = core, 2 = SOL. No other values exist.
        nmix = 0
        associate(fd => simgrid%data%fluxdata)
        do i = 1, fd%nFt
            k = fd%fluxtubecells(fd%fluxtubecellsP(i, 1))
            fd%fluxtuberegID(i) = simgrid%cell%reg(k)
            do lbl = fd%fluxtubecellsP(i, 1), &
                fd%fluxtubecellsP(i, 1) + fd%fluxtubecellsP(i, 2) - 1
                if (simgrid%cell%reg(fd%fluxtubecells(lbl)) /= &
                    fd%fluxtuberegID(i)) nmix = nmix + 1
            end do
        end do
        if (nmix > 0) then
            print *, 'SetSOLPSRegionsLIM: WARNING: ', nmix, &
                ' flux tube cells with mixed core/SOL regions ' // &
                '(tube crossing the separatrix?)'
        end if
        print *, 'SetSOLPSRegionsLIM: ftReg: ', &
            count(fd%fluxtuberegID(1:fd%nFt) == 1), ' core tubes, ', &
            count(fd%fluxtuberegID(1:fd%nFt) == 2), ' SOL tubes'
        end associate

    contains

        ! Minimum distance from a point to vessel sub-structure ks
        function DistToStructure(px, py, ks) result(dres)
            real(R8), intent(in)        :: px, py
            integer(I8), intent(in)     :: ks
            real(R8)                    :: dres
            integer(I8)                 :: jj, npt
            real(R8)                    :: x1, y1, x2, y2
            dres = posinfval_R8()
            npt = vessel%structures(ks)%np
            do jj = 1, npt - 1
                x1 = vessel%structures(ks)%x(jj)
                y1 = vessel%structures(ks)%y(jj)
                x2 = vessel%structures(ks)%x(jj + 1)
                y2 = vessel%structures(ks)%y(jj + 1)
                dres = min(dres, DistToSegment(px, py, x1, y1, x2, y2))
            end do
            if (vessel%structures(ks)%isclosed .and. (npt > 2)) then
                x1 = vessel%structures(ks)%x(npt)
                y1 = vessel%structures(ks)%y(npt)
                x2 = vessel%structures(ks)%x(1)
                y2 = vessel%structures(ks)%y(1)
                dres = min(dres, DistToSegment(px, py, x1, y1, x2, y2))
            end if
        end function

        ! Minimum distance from a point to one segment
        function DistToSegment(px, py, x1, y1, x2, y2) result(dres)
            real(R8), intent(in)    :: px, py, x1, y1, x2, y2
            real(R8)                :: dres
            real(R8)                :: dx, dy, t
            dx = x2 - x1
            dy = y2 - y1
            if ((dx*dx + dy*dy) <= 0.0_R8) then
                dres = sqrt((px - x1)**2 + (py - y1)**2)
                return
            end if
            t = ((px - x1)*dx + (py - y1)*dy)/(dx*dx + dy*dy)
            t = max(0.0_R8, min(1.0_R8, t))
            dres = sqrt((px - (x1 + t*dx))**2 + (py - (y1 + t*dy))**2)
        end function

    end subroutine

    ! SOLPS region conventions for the SN (single null) catalog
    ! topologies
    subroutine SetSOLPSRegionsSN(simgrid, topomesh, vessel, isupper)

        ! Description
        !============
        ! Assign the volume and face regions that SOLPS expects for a
        ! single-null case (NREG = 4, FREG = 13). West/East follow the
        ! convention with its reversed target ordering: West = inboard
        ! for a lower single null, outboard for an upper single null.
        !
        ! The construction starts from the 8 "rays" leaving the
        ! X-point grid vertex: 4 are the separatrix branches (aligned
        ! faces), the other 4 are chains of non-aligned (radial) grid
        ! faces walked outwards from the X-point until a boundary
        ! vertex: into the core (core cut, face region 5), into the
        ! PFR (divertor cut, face region 6, splitting the inner and
        ! outer divertor), and into the two SOL sectors (divertor
        ! entrances, face regions 2/3, splitting main SOL from
        ! divertor SOL).
        !
        ! Volume regions (cvReg):
        !   1 = core
        !   2 = main SOL
        !   3 = Western divertor (divertor SOL + PFR, West side)
        !   4 = Eastern divertor (divertor SOL + PFR, East side)
        !
        ! Face regions (fcReg):
        !   1 = Western target: wall faces on the structures given in
        !       goat.vessel.solps_lower_divertor_inner_target (lower
        !       SN) / solps_upper_divertor_outer_target (upper SN)
        !   2 = entrance to the Western divertor (radial chain from
        !       the X-point to the wall through the West SOL)
        !   3 = entrance to the Eastern divertor
        !   4 = Eastern target (solps_lower_divertor_outer_target /
        !       solps_upper_divertor_inner_target)
        !   5 = core cut (radial chain X-point -> core boundary)
        !   6 = divertor cut (radial chain X-point -> wall through
        !       the PFR, splitting the two divertors)
        !   7 = Western divertor PFR wall: boundary faces adjacent to
        !       PFR cells of region 3 that are not target faces
        !   8 = core boundary (boundary faces adjacent to core cells)
        !   9 = Eastern divertor PFR wall (as 7, for region 4)
        !  10 = separatrix (internal faces on the separatrix,
        !       including the divertor legs)
        !  11 = Western divertor main-chamber wall: boundary faces
        !       adjacent to divertor-SOL cells of region 3 that are
        !       not target faces (empty when the entrance chain
        !       lands on the target itself)
        !  12 = main-chamber (main SOL) wall and any remaining wall
        !  13 = Eastern divertor main-chamber wall (as 11, region 4)
        !
        ! Flux tube regions (ftReg): 1 = core, 2 = SOL, 3 = PFR.
        !
        ! Empty divertor target inputs are tolerated (warning): the
        ! corresponding target face region stays empty and its wall
        ! faces fall to regions 7/9/11/13/12. Faces are attributed to
        ! a wall structure only when they lie within dtolwall of it,
        ! so flux-aligned domain boundaries far from the wall cannot
        ! be mistaken for targets.

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(inout)    :: simgrid
        type(TopomeshUDT), intent(in)   :: topomesh
        type(VesselUDT), intent(in)     :: vessel
        logical, intent(in)             :: isupper

        ! Parameters: cell kinds and wall attribution tolerance
        integer(I8), parameter          :: KCORE = 1, KSOL = 2, &
            KPFR = 3
        real(R8), parameter             :: dtolwall = 5.0e-3_R8

        ! Auxiliary
        integer(I8), allocatable        :: cc(:), xp(:), tcv(:), &
            tfc(:), vf(:), cf(:), cvkind(:), comp(:), queue(:), &
            chainreg(:), compreg(:), compn(:), zn(:), zkind(:)
        logical, allocatable            :: issepTM(:), ismain(:), &
            istarget1(:), istarget4(:), ispfrhalf(:)
        real(R8), allocatable           :: compx(:), compy(:), &
            zpsi(:)
        integer(I8)                     :: i, j, k, ich, iv, ivX, w, &
            f, fprev, lbl, c, c2, ncomp, nsolcomp, nqh, nqt, nmix, &
            nsepal, npfrhalf, fstart(4), chlen(4), termvx(4), &
            ncand(4), nfreg(0:13), ncvreg(4), nftreg(3)
        real(R8)                        :: xc, yc, ux, uy, xxp, yxp, &
            d, dmin, score, dirx, diry, dp, s, psiX, psicore, &
            xm, ym, vx, vy, cxr, cyr
        character(8)                    :: chname(4)

        chname = ['core    ', 'PFR     ', 'West SOL', 'East SOL']
        print *, 'SetSOLPSRegionsSN: applying SOLPS single-null ' // &
            'regions (orientation: ' // &
            trim(merge('upper', 'lower', isupper)) // ')'

        ! X-point and reference directions
        !=================================
        xp = topomesh%GetXPointIDs()
        if (size(xp) /= 1) then
            print *, 'SetSOLPSRegionsSN: expected exactly 1 ' // &
                'X-point in the topomesh, found ', size(xp), &
                ' - leaving the generic regions in place'
            return
        end if
        xxp = topomesh%vert%x(xp(1))
        yxp = topomesh%vert%y(xp(1))

        ! Grid vertex sitting on the X-point
        ivX = 1
        dmin = posinfval_R8()
        do i = 1, simgrid%vert%ntot
            d = (simgrid%vert%x(i) - xxp)**2 + &
                (simgrid%vert%y(i) - yxp)**2
            if (d < dmin) then
                dmin = d
                ivX = i
            end if
        end do
        if (sqrt(dmin) > 1.0e-6_R8) then
            print *, 'SetSOLPSRegionsSN: WARNING: nearest grid ' // &
                'vertex to the X-point is ', sqrt(dmin), ' m away'
        end if
        psiX = simgrid%vert%psi(ivX)

        ! Core centroid: defines the reference direction
        ! u = (X-point -> core), whose left side is SOLPS West
        cc = topomesh%GetCoreCellIDs()
        if (size(cc) == 0) then
            print *, 'SetSOLPSRegionsSN: no core topomesh cells - ' // &
                'leaving the generic regions in place'
            return
        end if
        xc = 0.0_R8
        yc = 0.0_R8
        k = 0
        do i = 1, size(cc)
            tcv = topomesh%cell%GetVert(cc(i))
            xc = xc + sum(topomesh%vert%x(tcv))
            yc = yc + sum(topomesh%vert%y(tcv))
            k = k + size(tcv)
        end do
        xc = xc/real(k, kind=R8)
        yc = yc/real(k, kind=R8)
        ux = xc - xxp
        uy = yc - yxp

        ! Separatrix topomesh faces (single X-point: all of them
        ! belong to its separatrix, legs included)
        issepTM = (topomesh%face%type == TMfacesepID)

        ! Cell kinds
        !===========
        ! Core from the region marks set during the label
        ! translation; PFR/SOL per ZONE (the pre-overwrite cvReg
        ! values, i.e. the topomesh cells, which never straddle the
        ! separatrix) from the mean flux side relative to the
        ! X-point: the PFR flux is always core-side of its null.
        ! The zone average makes the test immune to interpolation
        ! noise on the cells hugging the separatrix.
        allocate(zpsi(maxval(simgrid%cell%reg)), &
            zn(maxval(simgrid%cell%reg)), &
            zkind(maxval(simgrid%cell%reg)))
        zpsi = 0.0_R8
        zn = 0
        do i = 1, simgrid%cell%ntot
            k = simgrid%cell%reg(i)
            zpsi(k) = zpsi(k) + simgrid%cell%psi(i)
            zn(k) = zn(k) + 1
        end do
        psicore = 0.0_R8
        k = 0
        do j = 1, size(zkind)
            if (zn(j) == 0) cycle
            if (mod(j - SOLPScoreregID, SOLPScoreregIDincr) == 0) &
                then
                psicore = psicore + zpsi(j)
                k = k + zn(j)
            end if
        end do
        if (k == 0) then
            print *, 'SetSOLPSRegionsSN: no core grid cells - ' // &
                'leaving the generic regions in place'
            return
        end if
        s = sign(1.0_R8, psiX - psicore/real(k, kind=R8))
        do j = 1, size(zkind)
            zkind(j) = KSOL
            if (zn(j) == 0) cycle
            if (mod(j - SOLPScoreregID, SOLPScoreregIDincr) == 0) &
                then
                zkind(j) = KCORE
            elseif (s*(zpsi(j)/real(zn(j), kind=R8) - psiX) < &
                0.0_R8) then
                zkind(j) = KPFR
            end if
        end do
        allocate(cvkind(simgrid%cell%ntot))
        do i = 1, simgrid%cell%ntot
            cvkind(i) = zkind(simgrid%cell%reg(i))
        end do

        ! The four radial rays leaving the X-point vertex
        !================================================
        ! Non-aligned internal faces at the X-point vertex, classified
        ! by the sector they walk into: core/PFR by the flux side
        ! (told apart by the dot product with u), West/East SOL by
        ! the cross product with u
        fstart = 0
        ncand = 0
        nsepal = 0
        vf = GetVertFace(simgrid%vert, ivX)
        do j = 1, size(vf)
            f = vf(j)
            if (simgrid%face%BF(f)) cycle
            if (simgrid%face%aligned(f) == 1) then
                nsepal = nsepal + 1
                cycle
            end if
            if (simgrid%face%vert(f, 1) == ivX) then
                w = simgrid%face%vert(f, 2)
            else
                w = simgrid%face%vert(f, 1)
            end if
            dirx = simgrid%vert%x(w) - xxp
            diry = simgrid%vert%y(w) - yxp
            dp = s*(simgrid%vert%psi(w) - psiX)
            if (dp < 0.0_R8) then
                if ((dirx*ux + diry*uy) > 0.0_R8) then
                    ich = 1
                else
                    ich = 2
                end if
            else
                if ((ux*diry - uy*dirx) > 0.0_R8) then
                    ich = 3
                else
                    ich = 4
                end if
            end if
            ncand(ich) = ncand(ich) + 1
            if (fstart(ich) == 0) fstart(ich) = f
        end do
        print *, 'SetSOLPSRegionsSN: X-point grid vertex ', ivX, &
            ' at ', xxp, yxp, ': ', nsepal, ' separatrix rays, ' // &
            'radial ray candidates (core/PFR/West/East): ', ncand
        do ich = 1, 4
            if (fstart(ich) == 0) then
                print *, 'SetSOLPSRegionsSN: WARNING: no radial ' // &
                    'ray into the ' // trim(chname(ich)) // &
                    ' sector at the X-point'
            end if
        end do

        ! Walk the four radial chains from the X-point to a boundary
        ! vertex. Their faces are tagged with the face region of the
        ! cut they realize: 5 = core cut, 6 = divertor (PFR) cut,
        ! 2 = West divertor entrance, 3 = East divertor entrance
        allocate(chainreg(simgrid%face%ntot))
        chainreg = 0
        chlen = 0
        termvx = 0
        do ich = 1, 4
            if (fstart(ich) == 0) cycle
            iv = ivX
            f = fstart(ich)
            do
                select case (ich)
                case (1)
                    chainreg(f) = 5
                case (2)
                    chainreg(f) = 6
                case (3)
                    chainreg(f) = 2
                case (4)
                    chainreg(f) = 3
                end select
                chlen(ich) = chlen(ich) + 1
                if (chlen(ich) > simgrid%face%ntot) then
                    print *, 'SetSOLPSRegionsSN: WARNING: runaway ' // &
                        trim(chname(ich)) // ' chain walk - aborted'
                    exit
                end if
                if (simgrid%face%vert(f, 1) == iv) then
                    w = simgrid%face%vert(f, 2)
                else
                    w = simgrid%face%vert(f, 1)
                end if
                if (simgrid%vert%BV(w)) then
                    termvx(ich) = w
                    exit
                end if

                ! Continue on the unvisited radial face that best
                ! keeps the walking direction
                vx = simgrid%vert%x(w) - simgrid%vert%x(iv)
                vy = simgrid%vert%y(w) - simgrid%vert%y(iv)
                fprev = f
                f = 0
                score = -posinfval_R8()
                cf = GetVertFace(simgrid%vert, w)
                do j = 1, size(cf)
                    if (cf(j) == fprev) cycle
                    if (simgrid%face%BF(cf(j))) cycle
                    if (simgrid%face%aligned(cf(j)) == 1) cycle
                    if (chainreg(cf(j)) /= 0) cycle
                    if (simgrid%face%vert(cf(j), 1) == w) then
                        c2 = simgrid%face%vert(cf(j), 2)
                    else
                        c2 = simgrid%face%vert(cf(j), 1)
                    end if
                    dirx = simgrid%vert%x(c2) - simgrid%vert%x(w)
                    diry = simgrid%vert%y(c2) - simgrid%vert%y(w)
                    d = (dirx*vx + diry*vy)/ &
                        max(sqrt(dirx**2 + diry**2), tiny(1.0_R8))
                    if (d > score) then
                        score = d
                        f = cf(j)
                    end if
                end do
                if (f == 0) then
                    print *, 'SetSOLPSRegionsSN: WARNING: ' // &
                        trim(chname(ich)) // ' chain stopped at ' // &
                        'internal vertex ', w, ' after ', &
                        chlen(ich), ' faces'
                    termvx(ich) = w
                    exit
                end if
                iv = w
            end do
        end do
        do ich = 1, 4
            if (termvx(ich) == 0) cycle
            xm = simgrid%vert%x(termvx(ich))
            ym = simgrid%vert%y(termvx(ich))
            k = 0
            dmin = posinfval_R8()
            do j = 1, int(vessel%nstructures, kind=I8)
                d = StructureDistance(vessel, xm, ym, j)
                if (d < dmin) then
                    dmin = d
                    k = j
                end if
            end do
            if (ich == 1) then
                print *, 'SetSOLPSRegionsSN: ', trim(chname(ich)), &
                    ' cut: ', chlen(ich), ' faces, ends at ', xm, ym
            else
                print *, 'SetSOLPSRegionsSN: ', trim(chname(ich)), &
                    ' cut: ', chlen(ich), ' faces, intercepts ' // &
                    'the wall at ', xm, ym, ' (structure ', k, ')'
            end if
        end do

        ! Volume regions
        !===============
        ! Flood fill of the SOL cells with the two entrance chains as
        ! barriers (main SOL + the two divertor SOLs), then of the
        ! PFR cells with the divertor cut as barrier (the two PFR
        ! halves)
        allocate(comp(simgrid%cell%ntot), queue(simgrid%cell%ntot))
        comp = 0
        ncomp = 0
        call FloodFill(KSOL, 2_I8, 3_I8)
        nsolcomp = ncomp
        call FloodFill(KPFR, 6_I8, 6_I8)
        if (nsolcomp /= 3) then
            print *, 'SetSOLPSRegionsSN: WARNING: ', nsolcomp, &
                ' SOL components instead of 3 (incomplete ' // &
                'entrance chains?)'
        end if

        allocate(ismain(ncomp), compreg(ncomp), compx(ncomp), &
            compy(ncomp), compn(ncomp), ispfrhalf(ncomp))
        ismain = .false.
        compreg = 0

        ! The two PFR halves are the components adjacent to the
        ! divertor cut. Any other flux-side component is a far-field
        ! pocket whose flux folded back to the core side (e.g.
        ! beyond a clipped secondary separatrix): those belong to
        ! the main SOL, not to a divertor
        ispfrhalf = .false.
        do i = 1, simgrid%face%ntot
            if (chainreg(i) /= 6) cycle
            tfc = GetFaceCell(simgrid%face, i)
            do k = 1, size(tfc)
                if (cvkind(tfc(k)) == KPFR) then
                    ispfrhalf(comp(tfc(k))) = .true.
                end if
            end do
        end do
        npfrhalf = 0
        do c = nsolcomp + 1, ncomp
            if (ispfrhalf(c)) npfrhalf = npfrhalf + 1
        end do
        nmix = 0
        do i = 1, simgrid%cell%ntot
            if (cvkind(i) /= KPFR) cycle
            if (ispfrhalf(comp(i))) cycle
            cvkind(i) = KSOL
            ismain(comp(i)) = .true.
            nmix = nmix + 1
        end do
        if (nmix > 0) then
            print *, 'SetSOLPSRegionsSN: ', nmix, ' flux-side ' // &
                'cells in ', ncomp - nsolcomp - npfrhalf, &
                ' components not connected to the PFR (folded ' // &
                'far-field flux?) - assigned to the main SOL'
        end if
        if (npfrhalf /= 2) then
            print *, 'SetSOLPSRegionsSN: WARNING: ', npfrhalf, &
                ' PFR halves instead of 2 (incomplete divertor ' // &
                'cut?) - splitting the PFR by the side of the ' // &
                'X-point -> core axis instead'
        end if

        ! Main SOL = the SOL component(s) radially adjacent to the
        ! core across the separatrix; the divertor SOL components
        ! are identified by the entrance chain bounding them
        do i = 1, simgrid%face%ntot
            if (simgrid%face%BF(i)) cycle
            tfc = GetFaceCell(simgrid%face, i)
            if (size(tfc) /= 2) cycle
            if ((cvkind(tfc(1)) == KCORE) .and. &
                (cvkind(tfc(2)) == KSOL)) then
                ismain(comp(tfc(2))) = .true.
            elseif ((cvkind(tfc(2)) == KCORE) .and. &
                (cvkind(tfc(1)) == KSOL)) then
                ismain(comp(tfc(1))) = .true.
            end if
        end do
        do i = 1, simgrid%face%ntot
            if ((chainreg(i) /= 2) .and. (chainreg(i) /= 3)) cycle
            tfc = GetFaceCell(simgrid%face, i)
            do k = 1, size(tfc)
                if (cvkind(tfc(k)) /= KSOL) cycle
                c = comp(tfc(k))
                if (ismain(c)) cycle
                if (chainreg(i) == 2) then
                    compreg(c) = 3
                else
                    compreg(c) = 4
                end if
            end do
        end do

        ! PFR components: West (region 3) / East (region 4) by the
        ! side of the X-point -> core direction
        compx = 0.0_R8
        compy = 0.0_R8
        compn = 0
        do i = 1, simgrid%cell%ntot
            if (cvkind(i) /= KPFR) cycle
            compx(comp(i)) = compx(comp(i)) + simgrid%cell%x(i)
            compy(comp(i)) = compy(comp(i)) + simgrid%cell%y(i)
            compn(comp(i)) = compn(comp(i)) + 1
        end do
        do c = nsolcomp + 1, ncomp
            if (.not. ispfrhalf(c)) cycle
            if (compn(c) == 0) cycle
            cxr = compx(c)/real(compn(c), kind=R8) - xxp
            cyr = compy(c)/real(compn(c), kind=R8) - yxp
            if ((ux*cyr - uy*cxr) > 0.0_R8) then
                compreg(c) = 3
            else
                compreg(c) = 4
            end if
        end do

        ! Cell region assignment
        nmix = 0
        do i = 1, simgrid%cell%ntot
            select case (cvkind(i))
            case (KCORE)
                simgrid%cell%reg(i) = 1
            case (KSOL)
                if (ismain(comp(i))) then
                    simgrid%cell%reg(i) = 2
                elseif (compreg(comp(i)) /= 0) then
                    simgrid%cell%reg(i) = compreg(comp(i))
                else
                    simgrid%cell%reg(i) = 2
                    nmix = nmix + 1
                end if
            case default
                if ((npfrhalf == 2) .and. &
                    (compreg(comp(i)) /= 0)) then
                    simgrid%cell%reg(i) = compreg(comp(i))
                else
                    ! Fallback split by the X-point -> core axis
                    cxr = simgrid%cell%x(i) - xxp
                    cyr = simgrid%cell%y(i) - yxp
                    if ((ux*cyr - uy*cxr) > 0.0_R8) then
                        simgrid%cell%reg(i) = 3
                    else
                        simgrid%cell%reg(i) = 4
                    end if
                end if
            end select
        end do
        if (nmix > 0) then
            print *, 'SetSOLPSRegionsSN: WARNING: ', nmix, &
                ' SOL cells in components touching no entrance ' // &
                'chain and no core - assigned to the main SOL'
        end if
        ncvreg = 0
        do i = 1, simgrid%cell%ntot
            k = simgrid%cell%reg(i)
            if ((k >= 1) .and. (k <= 4)) ncvreg(k) = ncvreg(k) + 1
        end do
        print *, 'SetSOLPSRegionsSN: cvReg counts: core ', &
            ncvreg(1), ', main SOL ', ncvreg(2), ', West ' // &
            'divertor ', ncvreg(3), ', East divertor ', ncvreg(4)

        ! Face regions
        !=============
        ! Western/Eastern target structures from the Goat inputs,
        ! with the manual's reversed ordering: lower SN -> West =
        ! inner target, upper SN -> West = outer target
        allocate(istarget1(vessel%nstructures), &
            istarget4(vessel%nstructures))
        istarget1 = .false.
        istarget4 = .false.
        if (isupper) then
            do i = 1, int(vessel%nsolpsudo, kind=I8)
                istarget1(vessel%solpsudoind(i)) = .true.
            end do
            do i = 1, int(vessel%nsolpsudi, kind=I8)
                istarget4(vessel%solpsudiind(i)) = .true.
            end do
        else
            do i = 1, int(vessel%nsolpsldi, kind=I8)
                istarget1(vessel%solpsldiind(i)) = .true.
            end do
            do i = 1, int(vessel%nsolpsldo, kind=I8)
                istarget4(vessel%solpsldoind(i)) = .true.
            end do
        end if
        if (count(istarget1) == 0) then
            print *, 'SetSOLPSRegionsSN: WARNING: no Western ' // &
                'target structures given (goat.vessel.solps_' // &
                trim(merge('upper_divertor_outer', &
                'lower_divertor_inner', isupper)) // '_target) - ' // &
                'face region 1 will be empty'
        end if
        if (count(istarget4) == 0) then
            print *, 'SetSOLPSRegionsSN: WARNING: no Eastern ' // &
                'target structures given (goat.vessel.solps_' // &
                trim(merge('upper_divertor_inner', &
                'lower_divertor_outer', isupper)) // '_target) - ' // &
                'face region 4 will be empty'
        end if

        simgrid%face%reg = 0
        nfreg = 0
        do i = 1, simgrid%face%ntot
            if (simgrid%face%BF(i)) then
                tfc = GetFaceCell(simgrid%face, i)
                c = tfc(1)
                if (cvkind(c) == KCORE) then
                    ! Core boundary
                    simgrid%face%reg(i) = 8
                else
                    ! Nearest wall structure (trusted only when the
                    ! face actually sits on the wall)
                    xm = 0.5_R8*(simgrid%vert%x(simgrid%face%vert(i, 1)) &
                        + simgrid%vert%x(simgrid%face%vert(i, 2)))
                    ym = 0.5_R8*(simgrid%vert%y(simgrid%face%vert(i, 1)) &
                        + simgrid%vert%y(simgrid%face%vert(i, 2)))
                    k = 0
                    dmin = posinfval_R8()
                    do j = 1, int(vessel%nstructures, kind=I8)
                        d = StructureDistance(vessel, xm, ym, j)
                        if (d < dmin) then
                            dmin = d
                            k = j
                        end if
                    end do
                    if (dmin > dtolwall) k = 0
                    if (k > 0) then
                        if (istarget1(k)) then
                            simgrid%face%reg(i) = 1
                        elseif (istarget4(k)) then
                            simgrid%face%reg(i) = 4
                        end if
                    end if
                    if (simgrid%face%reg(i) == 0) then
                        if (cvkind(c) == KPFR) then
                            ! PFR wall, West/East half
                            if (simgrid%cell%reg(c) == 3) then
                                simgrid%face%reg(i) = 7
                            else
                                simgrid%face%reg(i) = 9
                            end if
                        else
                            ! Main-chamber wall: divertor SOL parts
                            ! keep their own regions
                            select case (simgrid%cell%reg(c))
                            case (3)
                                simgrid%face%reg(i) = 11
                            case (4)
                                simgrid%face%reg(i) = 13
                            case default
                                simgrid%face%reg(i) = 12
                            end select
                        end if
                    end if
                end if
            else
                ! Internal faces: separatrix, then the four cuts
                lbl = simgrid%face%TMfacelabel(i)
                if ((lbl >= 1) .and. (lbl <= topomesh%face%ntot)) then
                    if (issepTM(lbl)) simgrid%face%reg(i) = 10
                end if
                if ((simgrid%face%reg(i) == 0) .and. &
                    (chainreg(i) /= 0)) then
                    simgrid%face%reg(i) = chainreg(i)
                end if
            end if
            nfreg(simgrid%face%reg(i)) = &
                nfreg(simgrid%face%reg(i)) + 1
        end do
        print *, 'SetSOLPSRegionsSN: fcReg counts (regions ' // &
            '1..13): ', nfreg(1:13)

        ! Flux tube regions
        !==================
        ! 1 = core, 2 = SOL (main + divertor), 3 = PFR
        nmix = 0
        nftreg = 0
        associate(fd => simgrid%data%fluxdata)
        do i = 1, fd%nFt
            k = fd%fluxtubecells(fd%fluxtubecellsP(i, 1))
            fd%fluxtuberegID(i) = cvkind(k)
            nftreg(cvkind(k)) = nftreg(cvkind(k)) + 1
            do lbl = fd%fluxtubecellsP(i, 1), &
                fd%fluxtubecellsP(i, 1) + fd%fluxtubecellsP(i, 2) - 1
                if (cvkind(fd%fluxtubecells(lbl)) /= cvkind(k)) &
                    nmix = nmix + 1
            end do
        end do
        end associate
        if (nmix > 0) then
            print *, 'SetSOLPSRegionsSN: WARNING: ', nmix, &
                ' flux tube cells with mixed core/SOL/PFR kinds ' // &
                '(tube crossing the separatrix?)'
        end if
        print *, 'SetSOLPSRegionsSN: ftReg: ', nftreg(1), &
            ' core tubes, ', nftreg(2), ' SOL tubes, ', nftreg(3), &
            ' PFR tubes'

    contains

        ! Flood fill of the cells of kind kindval, without crossing
        ! the chain faces tagged breg1/breg2, appending the new
        ! components to comp/ncomp (host-associated)
        subroutine FloodFill(kindval, breg1, breg2)
            integer(I8), intent(in)     :: kindval, breg1, breg2
            integer(I8)                 :: i0, jj, kk, cq, cn
            integer(I8), allocatable    :: cfl(:), tfcl(:)
            do i0 = 1, simgrid%cell%ntot
                if (cvkind(i0) /= kindval) cycle
                if (comp(i0) /= 0) cycle
                ncomp = ncomp + 1
                comp(i0) = ncomp
                nqh = 1
                nqt = 1
                queue(1) = i0
                do while (nqh <= nqt)
                    cq = queue(nqh)
                    nqh = nqh + 1
                    cfl = GetCellFace(simgrid%cell, cq)
                    do jj = 1, size(cfl)
                        if (simgrid%face%BF(cfl(jj))) cycle
                        if ((chainreg(cfl(jj)) == breg1) .or. &
                            (chainreg(cfl(jj)) == breg2)) cycle
                        tfcl = GetFaceCell(simgrid%face, cfl(jj))
                        do kk = 1, size(tfcl)
                            cn = tfcl(kk)
                            if (cn == cq) cycle
                            if (cvkind(cn) /= kindval) cycle
                            if (comp(cn) /= 0) cycle
                            comp(cn) = ncomp
                            nqt = nqt + 1
                            queue(nqt) = cn
                        end do
                    end do
                end do
            end do
        end subroutine

    end subroutine

    ! Minimum distance from a point to vessel sub-structure ks
    function StructureDistance(vessel, px, py, ks) result(dres)
        type(VesselUDT), intent(in) :: vessel
        real(R8), intent(in)        :: px, py
        integer(I8), intent(in)     :: ks
        real(R8)                    :: dres
        integer(I8)                 :: jj, npt
        real(R8)                    :: x1, y1, x2, y2
        dres = posinfval_R8()
        npt = vessel%structures(ks)%np
        do jj = 1, npt - 1
            x1 = vessel%structures(ks)%x(jj)
            y1 = vessel%structures(ks)%y(jj)
            x2 = vessel%structures(ks)%x(jj + 1)
            y2 = vessel%structures(ks)%y(jj + 1)
            dres = min(dres, SegmentDistance(px, py, x1, y1, x2, y2))
        end do
        if (vessel%structures(ks)%isclosed .and. (npt > 2)) then
            x1 = vessel%structures(ks)%x(npt)
            y1 = vessel%structures(ks)%y(npt)
            x2 = vessel%structures(ks)%x(1)
            y2 = vessel%structures(ks)%y(1)
            dres = min(dres, SegmentDistance(px, py, x1, y1, x2, y2))
        end if
    end function

    ! Minimum distance from a point to one segment
    function SegmentDistance(px, py, x1, y1, x2, y2) result(dres)
        real(R8), intent(in)    :: px, py, x1, y1, x2, y2
        real(R8)                :: dres
        real(R8)                :: dx, dy, t
        dx = x2 - x1
        dy = y2 - y1
        if ((dx*dx + dy*dy) <= 0.0_R8) then
            dres = sqrt((px - x1)**2 + (py - y1)**2)
            return
        end if
        t = ((px - x1)*dx + (py - y1)*dy)/(dx*dx + dy*dy)
        t = max(0.0_R8, min(1.0_R8, t))
        dres = sqrt((px - (x1 + t*dx))**2 + (py - (y1 + t*dy))**2)
    end function


end module
