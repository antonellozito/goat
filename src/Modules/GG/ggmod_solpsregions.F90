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
!  - SetSOLPSRegionsDN: CDN and DDN bottom/top (NREG = 8,
!    FREG = 26/27)

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
    public :: SetSOLPSRegionsLIM, SetSOLPSRegionsSN, &
        SetSOLPSRegionsDN

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
        !       SN) / solps_upper_divertor_outer_target (upper SN).
        !       The West and East target lists may share a structure
        !       (both legs striking one flat wall run): its faces are
        !       then split between regions 1 and 4 at the divertor
        !       cut (via the volume regions), and the PFR wall
        !       regions 7/9 in between do not exist.
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
        !  10 = separatrix (internal faces on the core-enclosing
        !       part of the separatrix only - the divertor legs
        !       are excluded)
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

        ! Separatrix topomesh faces: face region 10
        ! is only the core-enclosing part of the
        ! separatrix, so the interfaces between a core and
        ! a non-core topomesh cell - the divertor legs
        ! (PFR/SOL interfaces) are excluded
        allocate(issepTM(topomesh%face%ntot))
        issepTM = .false.
        do i = 1, topomesh%face%ntot
            if (topomesh%face%type(i) /= TMfacesepID) cycle
            tfc = topomesh%face%GetCell(i)
            if (size(tfc) == 2) then
                if (any(tfc(1) == cc) .neqv. any(tfc(2) == cc)) then
                    issepTM(i) = .true.
                end if
            end if
        end do

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
        ! with reversed ordering: lower SN -> West =
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
                        if (istarget1(k) .and. istarget4(k)) then
                            ! Both targets share this wall structure
                            ! (e.g. both legs strike one flat wall
                            ! run): split it at the divertor (PFR)
                            ! cut, which the volume regions already
                            ! encode - West divertor cells (3) face
                            ! the West target, East cells (4) the
                            ! East one. Outside the divertor volumes
                            ! fall back to the side of the X-point ->
                            ! core axis. The PFR wall regions between
                            ! the targets then do not exist.
                            if (simgrid%cell%reg(c) == 3) then
                                simgrid%face%reg(i) = 1
                            elseif (simgrid%cell%reg(c) == 4) then
                                simgrid%face%reg(i) = 4
                            elseif ((ux*(ym - yxp) - &
                                uy*(xm - xxp)) > 0.0_R8) then
                                simgrid%face%reg(i) = 1
                            else
                                simgrid%face%reg(i) = 4
                            end if
                        elseif (istarget1(k)) then
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


    ! SOLPS region conventions for the DN (double null) catalog
    ! topologies: CDN, DDN bottom and DDN top
    subroutine SetSOLPSRegionsDN(simgrid, topomesh, vessel, isddn)

        ! Description
        !============
        ! Assign the volume and face regions that SOLPS expects for a
        ! double-null case (NREG = 8; FREG = 26 for CDN, 27 for DDN),
        ! 
        ! on the quadrilateral-equivalent picture of the unstructured
        ! Goat grid. The divertor numbering goes clockwise starting
        ! from the lower inner target. Only the innermost X-point per
        ! lower/upper side takes part in the construction (consistent
        ! with ClassifyBasicSOLPSCatalogTopology).
        !
        ! Each of the two X-points contributes 4 separatrix branches
        ! and 4 radial chains of non-aligned grid faces walked
        ! outwards from the X-point vertex: into the core-directed
        ! sector (core cuts 9/11; for a DDN the secondary X-point's
        ! chain first crosses the SOL band between the two
        ! separatrices - those faces are region 13, the connection
        ! between right and left SOL - and continues through the
        ! core), into the PFR (PFR cuts 12/10), and into the two SOL
        ! sectors (divertor entrances 2/7 below, 3/6 above).
        !
        ! Volume regions (cvReg):
        !   1 = left (inboard) core        5 = right (outboard) core
        !   2 = left SOL                   6 = right SOL
        !   3 = lower inner divertor       4 = upper inner divertor
        !   7 = upper outer divertor       8 = lower outer divertor
        ! (divertor = divertor SOL + PFR half). The cores are the
        ! flood-fill components of the core cells with the two core
        ! cuts as barriers; everything else comes from a single
        ! flood fill of the non-core cells with all chains as
        ! barriers: the components adjacent to the core across the
        ! separatrix are the left/right SOL, the remaining ones the
        ! four divertors (sorted by side and radius).
        !
        ! Face regions (fcReg), FREG = 26 (CDN) / 27 (DDN); the
        ! DDN-only region 13 shifts every region >= 13 up by one:
        !   1/4/5/8   = lower inner / upper inner / upper outer /
        !               lower outer target (wall faces on the
        !               goat.vessel.solps_*_divertor_*_target
        !               structures). The inner and outer list of a
        !               side may share a structure (both legs
        !               striking one flat wall run): its faces are
        !               then split between the two target regions at
        !               that side's PFR cut (via the volume regions),
        !               and the PFR wall regions in between do not
        !               exist.
        !   2/3/6/7   = entrances to the lower inner / upper inner /
        !               upper outer / lower outer divertor
        !   9/11      = lower/upper core cut
        !   10/12     = upper/lower PFR cut
        !   13        = (DDN only) connection between right and left
        !               SOL: the between-separatrix segment of the
        !               secondary X-point's core-directed chain
        !   13/15/20/22 (+1 for DDN) = lower inner / upper inner /
        !               upper outer / lower outer PFR wall boundary
        !   14/21 (+1) = left/right core boundary
        !   16/23 (+1) = left/right separatrix (the core-enclosing
        !               separatrix only, split by the X-points/cuts)
        !   17/19/24/26 (+1) = lower inner / upper inner / upper
        !               outer / lower outer divertor main wall
        !   18/25 (+1) = left/right SOL main wall boundary
        !
        ! Flux tube regions (ftReg): 1 = core, 2 = SOL, 3 = PFR
        ! (PFR = divertor cells on the core-side flux of their own
        ! X-point, decided per zone).
        !
        ! Empty target inputs are tolerated (warning): the target
        ! face region stays empty and its wall faces fall to the
        ! PFR/divertor/SOL wall regions. Faces are attributed to a
        ! wall structure only within dtolwall.

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(inout)    :: simgrid
        type(TopomeshUDT), intent(in)   :: topomesh
        type(VesselUDT), intent(in)     :: vessel
        logical, intent(in)             :: isddn

        ! Parameters: cell kinds and wall attribution tolerance
        integer(I8), parameter          :: KCORE = 1, KSOL = 2, &
            KPFR = 3
        real(R8), parameter             :: dtolwall = 5.0e-3_R8

        ! Auxiliary
        integer(I8), allocatable        :: xps(:), op(:), cc(:), &
            ccf(:), tmpi(:), tcv(:), tfc(:), vf(:), cf(:), zn(:), &
            zoneid(:), cvkind(:), comp(:), queue(:), chainreg(:), &
            compreg(:), compn(:)
        logical, allocatable            :: iscorezone(:), &
            issepTM(:), iscoreadj(:), istLI(:), istLO(:), &
            istUI(:), istUO(:)
        real(R8), allocatable           :: zpsi(:), compx(:), &
            compy(:)
        integer(I8)                     :: i, j, k, ich, iv, w, f, &
            fprev, lbl, c, c2, nxp, ncomp, ncorecomp, nqh, nqt, &
            ilow, iup, off, nmix, n13c, kside, ivX(2), fstart(2, 4), &
            chlen(2, 4), termvx(2, 4), ncand(2, 4), nfreg(0:27), &
            ncvreg(8), nftreg(3)
        real(R8)                        :: xc, yc, fcore, psicore, &
            s, d, dmin, score, dirx, diry, dp, xm, ym, vx, vy, zm, &
            xxp(2), yxp(2), psiX(2)
        character(14)                   :: chname(4)
        character(5)                    :: sidename(2)

        chname = ['core-directed ', 'PFR cut       ', &
            'inner entrance', 'outer entrance']
        sidename = ['lower', 'upper']
        off = merge(1_I8, 0_I8, isddn)
        print *, 'SetSOLPSRegionsDN: applying SOLPS double-null ' // &
            'regions (' // trim(merge('DDN', 'CDN', isddn)) // ')'

        ! Reference center and core-side flux value
        !==========================================
        ! From the O-point vertex when the axis is meshed, otherwise
        ! from the core-cell centroid and the core boundary faces
        cc = topomesh%GetCoreCellIDs()
        if (size(cc) == 0) then
            print *, 'SetSOLPSRegionsDN: no core topomesh cells - ' // &
                'leaving the generic regions in place'
            return
        end if
        op = topomesh%GetOPointIDs()
        if (size(op) == 1) then
            xc = topomesh%vert%x(op(1))
            yc = topomesh%vert%y(op(1))
            fcore = topomesh%vert%fval(op(1))
        else
            allocate(tmpi(0))
            do i = 1, size(cc)
                tmpi = [tmpi, topomesh%cell%GetVert(cc(i))]
            end do
            xc = sum(topomesh%vert%x(tmpi))/real(size(tmpi), kind=R8)
            yc = sum(topomesh%vert%y(tmpi))/real(size(tmpi), kind=R8)
            deallocate(tmpi)
            ccf = topomesh%GetCoreFaceIDs()
            if (size(ccf) == 0) then
                print *, 'SetSOLPSRegionsDN: no core-side flux ' // &
                    'reference - leaving the generic regions in place'
                return
            end if
            fcore = 0.0_R8
            do i = 1, size(ccf)
                fcore = fcore + &
                    0.5_R8*(topomesh%vert%fval( &
                            topomesh%face%vert(ccf(i), 1)) + &
                            topomesh%vert%fval( &
                            topomesh%face%vert(ccf(i), 2)))
            end do
            fcore = fcore/real(size(ccf), kind=R8)
        end if

        ! The two X-points: innermost per side
        !=====================================
        xps = topomesh%GetXPointIDs()
        nxp = size(xps)
        ilow = 0
        iup = 0
        vx = 0.0_R8
        vy = 0.0_R8
        do i = 1, nxp
            d = abs(topomesh%vert%fval(xps(i)) - fcore)
            if (topomesh%vert%y(xps(i)) >= yc) then
                if ((iup == 0) .or. (d < vy)) then
                    iup = i
                    vy = d
                end if
            else
                if ((ilow == 0) .or. (d < vx)) then
                    ilow = i
                    vx = d
                end if
            end if
        end do
        if ((ilow == 0) .or. (iup == 0)) then
            print *, 'SetSOLPSRegionsDN: could not find one ' // &
                'X-point per side - leaving the generic regions ' // &
                'in place'
            return
        end if
        xxp(1) = topomesh%vert%x(xps(ilow))
        yxp(1) = topomesh%vert%y(xps(ilow))
        xxp(2) = topomesh%vert%x(xps(iup))
        yxp(2) = topomesh%vert%y(xps(iup))

        ! Grid vertices sitting on the X-points
        do k = 1, 2
            ivX(k) = 1
            dmin = posinfval_R8()
            do i = 1, simgrid%vert%ntot
                d = (simgrid%vert%x(i) - xxp(k))**2 + &
                    (simgrid%vert%y(i) - yxp(k))**2
                if (d < dmin) then
                    dmin = d
                    ivX(k) = i
                end if
            end do
            if (sqrt(dmin) > 1.0e-6_R8) then
                print *, 'SetSOLPSRegionsDN: WARNING: nearest ' // &
                    'grid vertex to the ' // sidename(k) // &
                    ' X-point is ', sqrt(dmin), ' m away'
            end if
            psiX(k) = simgrid%vert%psi(ivX(k))
        end do

        ! Zones and core cells
        !=====================
        ! Snapshot of the pre-overwrite regions (the topomesh cells,
        ! which never straddle a separatrix): zone mean flux for the
        ! PFR tests, core zones from the region marks
        allocate(zoneid(simgrid%cell%ntot))
        zoneid = simgrid%cell%reg
        allocate(zpsi(maxval(zoneid)), zn(maxval(zoneid)), &
            iscorezone(maxval(zoneid)))
        zpsi = 0.0_R8
        zn = 0
        do i = 1, simgrid%cell%ntot
            zpsi(zoneid(i)) = zpsi(zoneid(i)) + simgrid%cell%psi(i)
            zn(zoneid(i)) = zn(zoneid(i)) + 1
        end do
        iscorezone = .false.
        psicore = 0.0_R8
        k = 0
        do j = 1, size(zn)
            if (zn(j) == 0) cycle
            if (mod(j - SOLPScoreregID, SOLPScoreregIDincr) == 0) &
                then
                iscorezone(j) = .true.
                psicore = psicore + zpsi(j)
                k = k + zn(j)
            end if
        end do
        if (k == 0) then
            print *, 'SetSOLPSRegionsDN: no core grid cells - ' // &
                'leaving the generic regions in place'
            return
        end if
        s = sign(1.0_R8, psiX(1) - psicore/real(k, kind=R8))
        allocate(cvkind(simgrid%cell%ntot))
        do i = 1, simgrid%cell%ntot
            cvkind(i) = merge(KCORE, KSOL, iscorezone(zoneid(i)))
        end do

        ! The four radial rays leaving each X-point vertex
        !=================================================
        ! Classified by the sector they walk into: core-side flux
        ! rays split into core-directed (towards the core centroid)
        ! and PFR (away from it); SOL-side rays split into inner and
        ! outer by the side of the X-point -> core direction
        fstart = 0
        ncand = 0
        do k = 1, 2
            vx = xc - xxp(k)
            vy = yc - yxp(k)
            vf = GetVertFace(simgrid%vert, ivX(k))
            do j = 1, size(vf)
                f = vf(j)
                if (simgrid%face%BF(f)) cycle
                if (simgrid%face%aligned(f) == 1) cycle
                if (simgrid%face%vert(f, 1) == ivX(k)) then
                    w = simgrid%face%vert(f, 2)
                else
                    w = simgrid%face%vert(f, 1)
                end if
                dirx = simgrid%vert%x(w) - xxp(k)
                diry = simgrid%vert%y(w) - yxp(k)
                dp = s*(simgrid%vert%psi(w) - psiX(k))
                if (dp < 0.0_R8) then
                    if ((dirx*vx + diry*vy) > 0.0_R8) then
                        ich = 1
                    else
                        ich = 2
                    end if
                else
                    ! Inner = left of u for the lower X-point,
                    ! right of u for the upper one
                    if (((vx*diry - vy*dirx) > 0.0_R8) .eqv. &
                        (k == 1)) then
                        ich = 3
                    else
                        ich = 4
                    end if
                end if
                ncand(k, ich) = ncand(k, ich) + 1
                if (fstart(k, ich) == 0) fstart(k, ich) = f
            end do
            print *, 'SetSOLPSRegionsDN: ' // sidename(k) // &
                ' X-point grid vertex ', ivX(k), ' at ', xxp(k), &
                yxp(k), ': radial ray candidates (core-directed/' // &
                'PFR/inner/outer): ', ncand(k, :)
            do ich = 1, 4
                if (fstart(k, ich) == 0) then
                    print *, 'SetSOLPSRegionsDN: WARNING: no ' // &
                        trim(chname(ich)) // ' ray at the ' // &
                        sidename(k) // ' X-point'
                end if
            end do
        end do

        ! Walk the eight radial chains
        !=============================
        ! Face tags: entrances 2/7 (lower inner/outer), 3/6 (upper
        ! inner/outer); PFR cuts 12 (lower), 10 (upper); the
        ! core-directed chains are tagged per face: 9/11 (lower/
        ! upper core cut) next to core cells, 13 (between the two
        ! separatrices, DDN only) next to SOL cells
        allocate(chainreg(simgrid%face%ntot))
        chainreg = 0
        chlen = 0
        termvx = 0
        do k = 1, 2
            do ich = 1, 4
                if (fstart(k, ich) == 0) cycle
                iv = ivX(k)
                f = fstart(k, ich)
                do
                    select case (ich)
                    case (1)
                        tfc = GetFaceCell(simgrid%face, f)
                        if (any(iscorezone(zoneid(tfc)))) then
                            chainreg(f) = merge(9_I8, 11_I8, k == 1)
                        else
                            chainreg(f) = 13
                        end if
                    case (2)
                        chainreg(f) = merge(12_I8, 10_I8, k == 1)
                    case (3)
                        chainreg(f) = merge(2_I8, 3_I8, k == 1)
                    case (4)
                        chainreg(f) = merge(7_I8, 6_I8, k == 1)
                    end select
                    chlen(k, ich) = chlen(k, ich) + 1
                    if (chlen(k, ich) > simgrid%face%ntot) then
                        print *, 'SetSOLPSRegionsDN: WARNING: ' // &
                            'runaway chain walk - aborted'
                        exit
                    end if
                    if (simgrid%face%vert(f, 1) == iv) then
                        w = simgrid%face%vert(f, 2)
                    else
                        w = simgrid%face%vert(f, 1)
                    end if
                    if (simgrid%vert%BV(w)) then
                        termvx(k, ich) = w
                        exit
                    end if
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
                            max(sqrt(dirx**2 + diry**2), &
                            tiny(1.0_R8))
                        if (d > score) then
                            score = d
                            f = cf(j)
                        end if
                    end do
                    if (f == 0) then
                        print *, 'SetSOLPSRegionsDN: WARNING: ' // &
                            trim(chname(ich)) // ' chain of the ' // &
                            sidename(k) // ' X-point stopped at ' // &
                            'internal vertex ', w, ' after ', &
                            chlen(k, ich), ' faces'
                        termvx(k, ich) = w
                        exit
                    end if
                    iv = w
                end do
                if (termvx(k, ich) /= 0) then
                    print *, 'SetSOLPSRegionsDN: ' // sidename(k) // &
                        ' ' // trim(chname(ich)) // ' chain: ', &
                        chlen(k, ich), ' faces, ends at ', &
                        simgrid%vert%x(termvx(k, ich)), &
                        simgrid%vert%y(termvx(k, ich))
                end if
            end do
        end do

        ! Volume regions
        !===============
        ! Flood fill of the core cells with the core cuts as
        ! barriers (left/right core), then of all the non-core
        ! cells with all the chains as barriers (left/right SOL +
        ! the four divertors)
        allocate(comp(simgrid%cell%ntot), queue(simgrid%cell%ntot))
        comp = 0
        ncomp = 0
        call FloodFill(.true.)
        ncorecomp = ncomp
        call FloodFill(.false.)
        if (ncorecomp /= 2) then
            print *, 'SetSOLPSRegionsDN: WARNING: ', ncorecomp, &
                ' core components instead of 2 (incomplete core ' // &
                'cuts?)'
        end if
        if ((ncomp - ncorecomp) /= 6) then
            print *, 'SetSOLPSRegionsDN: WARNING: ', &
                ncomp - ncorecomp, ' non-core components instead ' // &
                'of 6 (incomplete chains?)'
        end if

        ! Component classification: centroids, core adjacency
        allocate(compx(ncomp), compy(ncomp), compn(ncomp), &
            compreg(ncomp), iscoreadj(ncomp))
        compx = 0.0_R8
        compy = 0.0_R8
        compn = 0
        do i = 1, simgrid%cell%ntot
            compx(comp(i)) = compx(comp(i)) + simgrid%cell%x(i)
            compy(comp(i)) = compy(comp(i)) + simgrid%cell%y(i)
            compn(comp(i)) = compn(comp(i)) + 1
        end do
        iscoreadj = .false.
        do i = 1, simgrid%face%ntot
            if (simgrid%face%BF(i)) cycle
            if (chainreg(i) /= 0) cycle
            tfc = GetFaceCell(simgrid%face, i)
            if (size(tfc) /= 2) cycle
            if ((cvkind(tfc(1)) == KCORE) .neqv. &
                (cvkind(tfc(2)) == KCORE)) then
                if (cvkind(tfc(1)) == KCORE) then
                    iscoreadj(comp(tfc(2))) = .true.
                else
                    iscoreadj(comp(tfc(1))) = .true.
                end if
            end if
        end do

        ! Region of each component: cores 1/5 and SOLs 2/6 by the
        ! side of the reference center, divertors 3/4/7/8 by side
        ! (below/above the center) and radius (inboard/outboard of
        ! their X-point)
        do c = 1, ncomp
            if (compn(c) == 0) then
                compreg(c) = 0
                cycle
            end if
            xm = compx(c)/real(compn(c), kind=R8)
            ym = compy(c)/real(compn(c), kind=R8)
            if (c <= ncorecomp) then
                compreg(c) = merge(1_I8, 5_I8, xm < xc)
            elseif (iscoreadj(c)) then
                compreg(c) = merge(2_I8, 6_I8, xm < xc)
            else
                kside = merge(1_I8, 2_I8, ym < yc)
                if (kside == 1) then
                    compreg(c) = merge(3_I8, 8_I8, xm < xxp(1))
                else
                    compreg(c) = merge(4_I8, 7_I8, xm < xxp(2))
                end if
            end if
        end do
        do i = 1, simgrid%cell%ntot
            simgrid%cell%reg(i) = compreg(comp(i))
        end do
        ncvreg = 0
        do i = 1, simgrid%cell%ntot
            k = simgrid%cell%reg(i)
            if ((k >= 1) .and. (k <= 8)) ncvreg(k) = ncvreg(k) + 1
        end do
        print *, 'SetSOLPSRegionsDN: cvReg counts: left core ', &
            ncvreg(1), ', left SOL ', ncvreg(2), ', lower inner ' // &
            'div ', ncvreg(3), ', upper inner div ', ncvreg(4), &
            ', right core ', ncvreg(5), ', right SOL ', ncvreg(6), &
            ', upper outer div ', ncvreg(7), ', lower outer div ', &
            ncvreg(8)

        ! Cell kinds for the flux tube regions and the PFR walls:
        ! divertor cells on the core-side flux of their own X-point
        ! (per zone) are PFR
        do i = 1, simgrid%cell%ntot
            k = simgrid%cell%reg(i)
            if ((k == 3) .or. (k == 8)) then
                kside = 1
            elseif ((k == 4) .or. (k == 7)) then
                kside = 2
            else
                cycle
            end if
            zm = zpsi(zoneid(i))/real(zn(zoneid(i)), kind=R8)
            if (s*(zm - psiX(kside)) < 0.0_R8) cvkind(i) = KPFR
        end do

        ! Face regions
        !=============
        ! Separatrix topomesh faces: only the core-enclosing part
        ! (interfaces between a core and a non-core topomesh cell)
        allocate(issepTM(topomesh%face%ntot))
        issepTM = .false.
        do i = 1, topomesh%face%ntot
            if (topomesh%face%type(i) /= TMfacesepID) cycle
            tfc = topomesh%face%GetCell(i)
            if (size(tfc) == 2) then
                if (any(tfc(1) == cc) .neqv. any(tfc(2) == cc)) then
                    issepTM(i) = .true.
                end if
            end if
        end do

        ! Target structures from the Goat inputs
        allocate(istLI(vessel%nstructures), &
            istLO(vessel%nstructures), istUI(vessel%nstructures), &
            istUO(vessel%nstructures))
        istLI = .false.
        istLO = .false.
        istUI = .false.
        istUO = .false.
        do i = 1, int(vessel%nsolpsldi, kind=I8)
            istLI(vessel%solpsldiind(i)) = .true.
        end do
        do i = 1, int(vessel%nsolpsldo, kind=I8)
            istLO(vessel%solpsldoind(i)) = .true.
        end do
        do i = 1, int(vessel%nsolpsudi, kind=I8)
            istUI(vessel%solpsudiind(i)) = .true.
        end do
        do i = 1, int(vessel%nsolpsudo, kind=I8)
            istUO(vessel%solpsudoind(i)) = .true.
        end do
        if (count(istLI) == 0) print *, 'SetSOLPSRegionsDN: ' // &
            'WARNING: goat.vessel.solps_lower_divertor_inner_' // &
            'target empty - face region 1 will be empty'
        if (count(istUI) == 0) print *, 'SetSOLPSRegionsDN: ' // &
            'WARNING: goat.vessel.solps_upper_divertor_inner_' // &
            'target empty - face region 4 will be empty'
        if (count(istUO) == 0) print *, 'SetSOLPSRegionsDN: ' // &
            'WARNING: goat.vessel.solps_upper_divertor_outer_' // &
            'target empty - face region 5 will be empty'
        if (count(istLO) == 0) print *, 'SetSOLPSRegionsDN: ' // &
            'WARNING: goat.vessel.solps_lower_divertor_outer_' // &
            'target empty - face region 8 will be empty'

        simgrid%face%reg = 0
        nfreg = 0
        n13c = 0
        do i = 1, simgrid%face%ntot
            if (simgrid%face%BF(i)) then
                tfc = GetFaceCell(simgrid%face, i)
                c = tfc(1)
                if (cvkind(c) == KCORE) then
                    ! Left/right core boundary
                    simgrid%face%reg(i) = &
                        merge(14_I8, 21_I8, &
                        simgrid%cell%reg(c) == 1) + off
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
                        if (istLI(k) .and. istLO(k)) then
                            ! Both lower targets share this wall
                            ! structure: split it at the lower PFR
                            ! cut, which the volume regions already
                            ! encode (3 = lower inner divertor,
                            ! 8 = lower outer; left/right SOL fall
                            ! to the matching side). The PFR wall
                            ! regions between the targets then do
                            ! not exist.
                            select case (simgrid%cell%reg(c))
                            case (3)
                                simgrid%face%reg(i) = 1
                            case (8)
                                simgrid%face%reg(i) = 8
                            case (6)
                                simgrid%face%reg(i) = 8
                            case default
                                simgrid%face%reg(i) = 1
                            end select
                        elseif (istUI(k) .and. istUO(k)) then
                            ! Both upper targets share this wall
                            ! structure: split it at the upper PFR
                            ! cut (4 = upper inner divertor,
                            ! 7 = upper outer)
                            select case (simgrid%cell%reg(c))
                            case (4)
                                simgrid%face%reg(i) = 4
                            case (7)
                                simgrid%face%reg(i) = 5
                            case (6)
                                simgrid%face%reg(i) = 5
                            case default
                                simgrid%face%reg(i) = 4
                            end select
                        else
                            if (istLI(k)) simgrid%face%reg(i) = 1
                            if (istUI(k)) simgrid%face%reg(i) = 4
                            if (istUO(k)) simgrid%face%reg(i) = 5
                            if (istLO(k)) simgrid%face%reg(i) = 8
                        end if
                    end if
                    if (simgrid%face%reg(i) == 0) then
                        select case (simgrid%cell%reg(c))
                        case (2)
                            simgrid%face%reg(i) = 18 + off
                        case (6)
                            simgrid%face%reg(i) = 25 + off
                        case (3)
                            simgrid%face%reg(i) = merge(13_I8, &
                                17_I8, cvkind(c) == KPFR) + off
                        case (4)
                            simgrid%face%reg(i) = merge(15_I8, &
                                19_I8, cvkind(c) == KPFR) + off
                        case (7)
                            simgrid%face%reg(i) = merge(20_I8, &
                                24_I8, cvkind(c) == KPFR) + off
                        case (8)
                            simgrid%face%reg(i) = merge(22_I8, &
                                26_I8, cvkind(c) == KPFR) + off
                        case default
                            simgrid%face%reg(i) = 18 + off
                        end select
                    end if
                end if
            else
                ! Internal faces: left/right separatrix by the core
                ! component behind them, then the chains
                lbl = simgrid%face%TMfacelabel(i)
                if ((lbl >= 1) .and. (lbl <= topomesh%face%ntot)) &
                    then
                    if (issepTM(lbl)) then
                        tfc = GetFaceCell(simgrid%face, i)
                        c = tfc(1)
                        if (cvkind(c) /= KCORE) then
                            if (size(tfc) == 2) c = tfc(2)
                        end if
                        simgrid%face%reg(i) = merge(16_I8, 23_I8, &
                            simgrid%cell%reg(c) == 1) + off
                    end if
                end if
                if ((simgrid%face%reg(i) == 0) .and. &
                    (chainreg(i) /= 0)) then
                    if ((chainreg(i) == 13) .and. .not. isddn) then
                        ! No between-separatrix connection exists in
                        ! a CDN: leave unassigned and warn
                        n13c = n13c + 1
                    else
                        simgrid%face%reg(i) = chainreg(i)
                    end if
                end if
            end if
            nfreg(simgrid%face%reg(i)) = &
                nfreg(simgrid%face%reg(i)) + 1
        end do
        if (n13c > 0) then
            print *, 'SetSOLPSRegionsDN: WARNING: ', n13c, &
                ' between-separatrix chain faces in a CDN ' // &
                '(imperfectly connected null?) - left unassigned'
        end if
        print *, 'SetSOLPSRegionsDN: fcReg counts (regions 1..', &
            26 + off, '): ', nfreg(1:26 + off)

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
            print *, 'SetSOLPSRegionsDN: WARNING: ', nmix, &
                ' flux tube cells with mixed core/SOL/PFR kinds ' // &
                '(tube crossing a separatrix?)'
        end if
        print *, 'SetSOLPSRegionsDN: ftReg: ', nftreg(1), &
            ' core tubes, ', nftreg(2), ' SOL tubes, ', nftreg(3), &
            ' PFR tubes'

    contains

        ! Flood fill of the core (docore) or non-core (.not.docore)
        ! cells, without crossing any chain face, appending the new
        ! components to comp/ncomp (host-associated)
        subroutine FloodFill(docore)
            logical, intent(in)         :: docore
            integer(I8)                 :: i0, jj, kk, cq, cn
            integer(I8), allocatable    :: cfl(:), tfcl(:)
            do i0 = 1, simgrid%cell%ntot
                if ((cvkind(i0) == KCORE) .neqv. docore) cycle
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
                        if (chainreg(cfl(jj)) /= 0) cycle
                        tfcl = GetFaceCell(simgrid%face, cfl(jj))
                        do kk = 1, size(tfcl)
                            cn = tfcl(kk)
                            if (cn == cq) cycle
                            if ((cvkind(cn) == KCORE) .neqv. &
                                docore) cycle
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
