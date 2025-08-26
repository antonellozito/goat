!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains all routines related to the driver for the grid adaptation module, such as initialization and postprocessing routines.

module gamod_driver

    ! Initialize
    !============
    ! Load modules
    use goatmod_types
    use goatmod_userinput
    use gamod_utility
    use gamod_types

    ! The usual
    implicit none
    save
    public     
    
    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!  
    
    contains

    subroutine GridAdaptor(grid,environment,magneticField,options)

        ! Description
        !============
        ! Overarching driver for grid adaptation, one lever lower than GADriver. ! Adapts topology of unstructured grid with the goal to improve grid  quality based on grid metric and user inputs

        ! Initialize
        !===========
        ! Declare modules  
        
        ! Declare variables
        !==================
        ! Arguments
        type(GAGridUDT), intent(inout)              :: grid
        type(EnvironmentUDT), intent(in)            :: environment
        type(MagneticFieldUDT), intent(in)          :: magneticField     
        type(GAoptionsUDT), intent(in)              :: options

        ! Initialize grid adaptation
        !===========================
        call GAinit(grid,options,environment,magneticField)
        

        ! Driver Selection
        !=================
        select case (options%meth)

        case ('simple')

            ! Regular grid adaption
            call GAInternalDriver(grid,options,environment,magneticField)

        case ('aposteriori')

            ! Grid adaptation based on simulation information
            ! call GAapostDriver

        case default

            ! Call error handler
            call gdErrorHandler('GridAdaptor: unknown driver option')
        
        end select

        ! Postprocessing
        !===============
        call PostProcessGA(grid,magneticField,options)

    end subroutine

    subroutine GAinit(grid,options,environment,magneticField)

        ! Description
        !============

        ! Declare variables
        !==================
        ! Arguments
        type(GAGridUDT), intent(inout)          :: grid 
        type(GAoptionsUDT), intent(in)          :: options
        type(EnvironmentUDT), intent(in)        :: environment
        type(MagneticFieldUDT), intent(in)      :: magneticField

        ! Variables
        integer(I8) :: i
        integer(I8), allocatable, dimension(:) :: tv, cvLookUp
        logical :: cells(grid%cell%ntot), is_ordered(grid%cell%ntot), &
            use_nsep, use_sepID, start
        character(:), allocatable :: base_func


        ! Initialize
        !===========
 
        associate(&
            c  => grid%cell, &
            v  => grid%vert &
            )
        ! Recompute cell centers
        do i = 1, c%ntot
            ! Get cell vertices
            tv = GetCellVertGA(c, i)

            ! Compute coordinates
            call c%x%SetSingleElement(i,sum(v%x%GetMultipleElements(tv))/real(size(tv), kind=R8)) 
            call c%y%SetSingleElement(i,sum(v%y%GetMultipleElements(tv))/real(size(tv), kind=R8)) 

        end do


        ! Check order of vertices  (see GetGeo_usCouples.m)
        cells = .true.
        call grid%CheckVertOrder(is_ordered, cells) 

        if (.not. all(is_ordered)) then

            ! Do ordening
            call grid%ReorderCellConn(is_ordered)

        end if 

        ! Get fsVx from fsFc
        call grid%GetFsVxFromFsFc()

        ! Determine Xpoints and separatrices
        cvLookUp = GetCvLookUpGA(c)
        use_nsep = .false.
        call grid%GiveXpoints(use_nsep,cvLookUp)
        use_sepID = .false.
        start = .true.
        call grid%GiveSeparatrices(use_nsep, use_sepID, start, cvLookUp)

        ! Identify aligned faces
        call grid%IdentifyAlignedFaces(options,magneticField)

        ! Set up the distance functions
        if (options%dist_function) then

            ! Define base function to compute the distance function
            base_func = options%base_func !'exp(-dist/d)'

            ! Free distance function based on user-defined input
            grid%fun%d_char_type = options%d_char_type
            grid%fun%dist_type   = options%dist_type
            grid%fun%d_rescale   = options%d_rescale
            call grid%fun%ComputeDistanceFunction(grid,options,base_func)

            ! Distance function for high poloidal flux next to the separatrix
            grid%fun_r%d_char_type = options%d_char_type
            grid%fun_r%dist_type   = 'pol_flux_est'
            grid%fun_r%d_rescale   = options%d_rescale
            call grid%fun_r%ComputeDistanceFunction(grid,options,base_func)

            ! Distance function for wall proximity
            grid%fun_wall%d_char_type = options%d_char_type
            grid%fun_wall%dist_type = options%dist_type_wall
            grid%fun_wall%d_rescale = options%d_rescale_wall
            call grid%fun_wall%ComputeDistanceFunction(grid,options,base_func)

        end if 


        ! Check the connectivity
        if (options%debug) then
            call grid%CheckUnstructuredGrid(.false.)
        end if

        ! correct face labels - TODO

        ! Remove some connectivity field - TODO

        ! Identify farSOL cells - TODO

        end associate        


    end subroutine

    subroutine GAInternalDriver(grid,options,environment,magneticField)

        ! Description
        !============
        ! Internal driver for the grid adaptation where all real adaptation take place such as removal of small triangles, stacked triangles

        ! Declare variables
        !==================
        ! Argument
        type(GAGridUDT), intent(inout)       :: grid
        type(GAoptionsUDT), intent(in)       :: options
        type(EnvironmentUDT), intent(in)     :: environment
        type(MagneticFieldUDT), intent(in)   :: magneticField

        ! Calculate quality metric

        ! Remove Small triangles

        ! Remove flux tubes with only two triangles

        ! Stacked to cutcell

        ! Splitting non-alinged quads

        ! Splitting trapezoids in concave shaved-off flux tube

        ! Splitting  and merging

        ! Stacked triangles


        ! Remove sticking out triangles

        ! Remove boundary flux tubes with only two triangles

        ! Remove stickout quad

        ! Boundary layer grid




    end subroutine
 
    subroutine PostProcessGA(grid,magneticField,options)

        ! Description
        !============
        ! Postprocessing specifically after grid adaptations

        ! Declare variables
        !==================
        ! Arguments
        type(GAGridUDT), intent(inout)      :: grid
        type(MagneticFieldUDT), intent(in)  :: magneticField
        type(GAoptionsUDT), intent(in)         :: options

        ! Auxiliary
        logical :: is_ordered(grid%cell%ntot), cells(grid%cell%ntot), &
         use_nsep, use_sepID, start
        integer(I8), allocatable :: cvLookUp(:), fcs(:), f_ord(:,:), nf(:), &
            lbls(:), lbls2(:), fsVx(:)
        integer(I8) ::  i, iv, nl, nvi, lb, nind, fcReg(grid%face%ntot), &
            fcLbl_loc(grid%face%ntot), indFc(grid%face%ntot), &
            ind(grid%face%ntot), nflbl
        

        ! Check consistency
        if (options%debug) then
            call grid%CheckUnstructuredGrid(.false.)
        end if

        ! Order correctly
        cells = .true.
        call grid%CheckVertOrder(is_ordered,cells) 
        call grid%ReorderCellConn(is_ordered)

        ! Determine number of Xpoints and separatrices
        cvLookUp = GetCvLookUpGA(grid%cell)
        use_nsep = .true.
        use_sepID = .true. 
        call grid%GiveSeparatrices(use_nsep, use_sepID, start, cvLookUp)
        call grid%GiveXpoints(use_sepID,cvLookUp)

        ! Check flux surfaces for possible merging of open surfaces - only very specific cases for vesselmode
        ! TODO
        call grid%MergeFS()

        ! Recalculate bx, by
        call grid%RecalcMagn(magneticField)

        ! Check whether every vertices is an flux surface for the correct cases
        if (options%stacked_trias .and. .not.options%vesselmode) then

            ! Get all vertices belonging to a flux surface
            fsVx = grid%data%fluxdata%fluxsurfaceverts%GetAllElements()

            do iv = 1, grid%vert%ntot

                nvi = count(fsVx == iv)

                if (nvi /= 1) then

                    ! Visualize TODO
                    print *, grid%vert%x%Get(iv)
                    print *, grid%vert%y%Get(iv)
                    call gdErrorHandler('PostprocessGA: vertex does not occur once in fsVx')

                end if 

            end do

        end if

        ! fcReg
        !------
        ! inner main target => fcReg == 1 
        ! outer main target => fcReg == 4 
        ! inner secondary target => fcReg == 5
        ! inner secondary target => fcReg == 8  
        
        ! Locally 
        fcLbL_loc = GetfcLblGA(grid%face,options)

        call Unique(fcLbl_loc, lbls)
        lbls2 = pack(lbls,lbls /= 0)
        nl = size(lbls2)

        if (nl .gt. size(options%fcRegmappingGA)) then
            call gdErrorHandler('PostProcesGA: fcReg mapping not compitable for GA labels,(more than 2 divertors)')
        end if

        ! Reset fcReg to zero and apply at the right faces
        fcReg = 0
        indFc = (/ (i, i=1,grid%face%ntot) /)
        do i = 1, nl
            lb = lbls2(i)
            nind = count(fcLbl_loc == lb)
            ind(1:nind) = pack(indFc,fcLbl_loc == lb )
            fcReg(ind(1:nind)) = options%fcRegmappingGA(lb)
        end do

        ! Self-check if faces with fcReg label can be chained together
        do i = 1, size(options%fcRegmappingGA)

            if (options%fcRegmappingGA(i) /= 0) then
                nflbl = count(fcReg == options%fcRegmappingGA(i))
                allocate(fcs(nflbl))
                fcs = pack(indFc,fcReg == options%fcRegmappingGA(i))

                if (size(fcs) /= 0) then

                    ! Build chaines of faces
                    call grid%face%ChainFaces(fcs, f_ord, nf)
                    
                    ! One fcReg number should only have one chain
                    if (size(nf) > 1) then

                        ! Visualize TODO
                        call gdErrorHandler("PostprocessGA: more than one chain of faces detect with the fcReg value")

                    end if

                end if

                ! House keeping
                deallocate(fcs)


            end if
        end do

        call grid%face%reg%SetAllElementsArray(fcReg)

        ! Visualize end grid - TODO



    end subroutine

    subroutine PostProcessingGridInformation(grid,magneticField,options)

        ! Description
        !============
        ! Post-process grid information to receive correct fields for
        ! traduit.out.b2us-file. This includes:
        ! - psi, bx, by, ffbz at the vertex locations
        ! - x, y coordinates of cells
        ! - bp, bt according to b2ag

        ! Notes
        !======
        ! Note 1: as implemented now, the grid.psi values may be inconsistent with
        ! the new grid coordinates, as first grid.psi is computed based on the old
        ! cell coordinates.  A similar issue is with the t0 term when recalculating 
        ! bp and bt. Needs revision        

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(inout)        :: grid
        type(MagneticFieldUDT), intent(in)  :: magneticField
        type(GAoptionsUDT), intent(inout)      :: options 

        ! Auxiliary
        logical :: trias_present, pents_present
        integer(I8):: vxs(1:100), s, nv, ic, i, ifs
        real(R8) :: bp(grid%cell%ntot), bt(grid%cell%ntot), &
            bpvx(grid%vert%ntot), fsPsi(grid%data%fluxdata%nFs), &
            r
        character(:), allocatable :: ctype

        ! Initialize
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        ! Calculate additional connectivity
        call ComputeGridInterconnections(grid)

        ! Recalculate magnetic field
        call magneticField%interp%Evaluate(v%x,v%y,0,0,v%psi)
        call magneticField%interp%Evaluate(v%x,v%y,1,0,v%bx)
        call magneticField%interp%Evaluate(v%x,v%y,0,1,v%by)

        ! Recalculate ffbz (constant)
        v%ffbz = v%ffbz(1)
        if (magneticField%RBtor /= 0) then
            v%ffbz = 2*pi_R8*magneticField%RBtor
        end if 

        ! Cells
        !======
        ! Recalculate cell centers and bp and bt
        bp = 0
        bt = 0
        do ic = 1, c%ntot
            vxs = GetCellVert(c, ic)
            nv = c%vertP(ic,2)

            ! Compute cell centers as average of vertex coordinates
            c%x(ic) = sum(v%x(vxs(1:nv)))/real(nv, kind=R8)
            c%y(ic) = sum(v%y(vxs(1:nv)))/real(nv, kind=R8)

            ! Recompute bp and bt according to b2agbb (in b2ag (SOLPS-ITER))
            do i = 1, nv
                bpvx(i) = sqrt( v%bx(i)**2 + v%by(i)**2 )
            end do
            r = nv*2*pi_R8*c%x(ic)    ! For axissymmetry around Z = 0 ! TODO
            bp(ic) = -sum(bpvx(1:nv))/r
            bt(ic) = sum(v%ffbz(vxs(1:nv)))/r
            bpvx = 0
        end do

        ! Recalculate psi at cell centers
        call magneticField%interp%Evaluate(c%x,c%y,0,0,c%psi)

        ! Get number of guard cells
        c%ngc = count(f%label /= 0)

        ! Ordening should be fine - done in PostProcessGA

        ! Build flux tube data
        !=====================
        trias_present = any(c%vertP(:,2) == 3)
        pents_present = any(c%vertP(:,2) == 5)

        if (.not.options%rem_small_trias) then
        
            print *, "Warning: Postprocessing: BuildFluxTube goes to shit if small triangle present"

        else 

            if (.not.trias_present .and. .not.pents_present) then

                ! For grids that can be stored as structured grid
                ctype = 'all'
                call BuildFluxTubeData(grid,options,magneticField,ctype)

            elseif (trias_present .and. .not. pents_present) then

                ! For grids with triangles (so fully unstructured) but without pentagons
                ctype = 'all_us'
                call BuildFluxTubeData(grid,options,magneticField,ctype); 

            elseif (pents_present) then

                ! For unstructured grids it is sufficient to give the closed flux tubes
                ! in the core and the first tube in the SOL of cells neighbouring the core. 
                ctype = 'closed'
                call BuildFluxTubeData(grid,options,magneticField,ctype);             

            end if
        end if 

        ! Recalculate fsPsi
        fsPsi = 0
        vxs = 0
        do ifs = 1, grid%data%fluxdata%nFs
            s = grid%data%fluxdata%fluxsurfacevertsP(ifs,1)
            nv = grid%data%fluxdata%fluxsurfacevertsP(ifs,2)
            vxs(1:nv) = grid%data%fluxdata%fluxsurfaceverts(s:s+nv-1)
            fsPsi(ifs) = sum(v%psi(vxs(1:nv))) / real(nv, kind=R8)
        end do

        ! Determine OMP and IMP - TODO
        call DetermineMPs(grid, options)
        

        end associate

    end subroutine

    subroutine CarryOverOptions(goatoptions, gaoptions)

        ! Description
        !============
        ! Carry over some options from goatoptions to gaoptions

        ! Declare variables
        !==================
        ! Arguments
        type(GoatoptionsUDT), intent(in)    :: goatoptions
        type(GAoptionsUDT), intent(inout)   :: gaoptions
        
        gaoptions%vesselmode            = goatoptions%vesselmode 
        gaoptions%slab                  = goatoptions%slab
        gaoptions%debug                 = goatoptions%debug 
        gaoptions%facelabelmappingGG    = goatoptions%GGtoGAfacelabelmappingGG
        gaoptions%facelabelmappingGA    = goatoptions%GGtoGAfacelabelmappingGA
        gaoptions%OMP_r                 = goatoptions%OMP_r
        gaoptions%OMP_z                 = goatoptions%OMP_z
        gaoptions%IMP_r                 = goatoptions%IMP_r
        gaoptions%IMP_z                 = goatoptions%IMP_z

    end subroutine

end module