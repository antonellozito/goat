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

module ggmod_solpsregions

    ! Load modules
    use mod_precision
    use mod_errorhandler
    use mod_definitions
    use mod_constants, only: posinfval_R8
    use goatmod_types, only: GridUDT, VesselUDT
    use ggmod_topology2D

    implicit none
    private
    public :: SetSOLPSRegionsLIM

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


end module
