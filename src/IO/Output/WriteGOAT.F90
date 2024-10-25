subroutine WriteGOAT(goatoptions, grid, magneticField, environment)

    ! Description
    !============
    ! This routine writes out the grid data from GOAT to the output file
    ! specified in the goatoptions structure. The format is an 
    ! unstructured traduit.out.b2us file 

    ! Note: the format madness that you will find in this file is due to
    ! the overly pedantic read/write routines of SOLPS, which fail if 
    ! you dare to even add/remove a single whitespace. If, hopefully
    ! somewhere in the near future, this would get revision, one can
    ! get rid of all the 'fmt' business and just put 'write(fu, *)' 
    ! probably everywhere. 

    ! Modules
    !========
    use goatmod_types
    use goatmod_userinput
    use mod_std_formatspecs
    use mod_constants

    ! The usual
    !==========
    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    type(GoatoptionsUDT), intent(in)        :: goatoptions 
    type(GridUDT), intent(inout)            :: grid 
    type(MagneticFieldUDT), intent(in)      :: magneticField
    type(EnvironmentUDT), intent(in)        :: environment

    ! Auxiliary
    character(:), allocatable               :: version, tempstring, fmt
    integer                                 :: fu 

    ! Loop
    integer(I8)                             :: i 

    ! Initialize
    !===========
    ! Check if we need to write
    if (.not. goatoptions%write_final) then 
        ! Return
        return 
    end if 

    ! Open file, overwrite if existing
    open (action='write', file=trim(goatoptions%writefilepath), newunit=fu, &
        status='unknown')

    ! Recompute data
    call UpdateGridData(grid, magneticField, environment)

    ! Write grid (unstructured traduit format)
    !=========================================
    ! Associate
    associate(&
        mf              => magneticField%interp,    &
        xv              => grid%vert%x,         &
        yv              => grid%vert%y,         &
        nv              => grid%vert%ntot,      &
        psiv            => grid%vert%psi,       &
        Bxv             => grid%vert%bx,        &
        Byv             => grid%vert%by,        &
        ffbzv           => grid%vert%ffbz,      &
        nvertcell       => grid%vert%ncell,     &
        nvertface       => grid%vert%nface,     &
        Xpoint          => grid%data%xpointID,  &
        nX              => grid%data%nxp,       &
        nf              => grid%face%ntot,      &
        facevert        => grid%face%vert,      &
        labelf          => grid%face%label,     &
        regf            => grid%face%reg,       &
        alignedf        => grid%face%aligned,   &
        fcOMP           => grid%data%OMPface,   &
        nfcOMP          => grid%data%nOMPface,  &    
        fcIMP           => grid%data%IMPface,   &
        nfcIMP          => grid%data%nIMPface,  &
        nFmxCv          => grid%face%ncell,     &
        nc              => grid%cell%ntot,      &
        ngc             => grid%cell%ngc,       & 
        cellvert        => grid%cell%vert,      &
        cellvertP       => grid%cell%vertP,     &
        xc              => grid%cell%x,         &
        yc              => grid%cell%y,         &
        psic            => grid%cell%psi,       &
        bpc             => grid%cell%bp,        &
        btc             => grid%cell%bt,        &
        flagc           => grid%cell%cflags,    &
        regc            => grid%cell%reg,       &
        cellface        => grid%cell%face,      &
        ncellvert       => grid%cell%nvert,     &
        ncellface       => grid%cell%nface,     &
        cellft          => grid%cell%ft,        &    
        cvOMP           => grid%data%OMPcell,       &
        ncvOMP          => grid%data%nOMPcell,      &
        cvIMP           => grid%data%IMPcell,       &
        ncvIMP          => grid%data%nIMPcell,      &
        nft             => grid%data%fluxdata%nft,                  &
        ftcell          => grid%data%fluxdata%fluxtubecells,        &
        ftcellP         => grid%data%fluxdata%fluxtubecellsP,       &
        ftface          => grid%data%fluxdata%fluxtubefaces,        &
        ftfaceP         => grid%data%fluxdata%fluxtubefacesP,       &
        ftreg           => grid%data%fluxdata%fluxtuberegID,        &
        nfs             => grid%data%fluxdata%nfs,                  &
        fsface          => grid%data%fluxdata%fluxsurfacefaces,     &
        fsfaceP         => grid%data%fluxdata%fluxsurfacefacesP,    &
        fspsi           => grid%data%fluxdata%fluxsurfacepsi,       &
        sgnx            => grid%data%sglegacy%nx,       &
        sgny            => grid%data%sglegacy%ny,       &
        sgncut          => grid%data%sglegacy%nncut,    &
        isClassicalGrid => grid%data%sglegacy%isclassicalgrid       &
        )

    ! Version
    version = 'VERSION03.002.000 Matlab traduit.out.b2us'
    write(fu, '(a)') version 

    ! General grid information
    tempstring = '*cf:    int                6    nCi,nFc,nVx,nCg,nFs,nFt'
    write(fu, '(a)' ) tempstring
    fmt = '(6'//Ifm//')'
    write(fu, fmt) nc, nf, nv, ngc, nfs, nft

    tempstring = '*cf:    int                5    nCmxVx,nCmxFc,nFmxCv,nVmxCv,nVmxFc'
    write(fu, '(a)' ) tempstring

    fmt = '(5'//Ifm//')'
    write(fu, fmt) ncellvert, ncellface, nFmxCv, nvertcell, nvertface 
    tempstring = '*cf:    int                1    isClassicalGrid'
    write(fu, '(a)' ) tempstring

    fmt = '('//Ifm//')'
    write(fu, fmt) isClassicalGrid

    if (isClassicalGrid == 1) then 
        tempstring = '*cf:    int                3    nx,ny,nncut'
        write(fu, '(a)' ) tempstring
        fmt = '(3'//Ifm//')'
        write(fu, fmt) sgnx, sgny, sgncut
    end if 

    ! X-point(s) information
    if (goatoptions%write_Xpointdata) then 
        tempstring = '*cf:    int                1    nX'
        write(fu, '(a)' ) tempstring
        fmt = '('//Ifm//')'
        write(fu, fmt) nX
    end if 

    ! OMP and IMP lengths
    if (goatoptions%write_OMPdata) then 
        tempstring = '*cf:    int                4    ncvOMP,ncvIMP,nfcOMP,nfcIMP'
        write(fu, '(a)' ) tempstring
        fmt = '(4'//Ifm//')'
        write(fu, fmt) ncvOMP, ncvIMP, nfcOMP, nfcIMP 
    end if 

    ! Vertex information
    if (goatoptions%write_Xpointdata) then 
        tempstring = repeat(' ', 1000)
        write(tempstring, '(a)') '*cf:    int                ', nX, '    Xpoint(s)'
        fmt = '('//repeat(Ifm, nX)//')'
        write(fu, fmt) trim(tempstring), Xpoint
    end if

    tempstring = '*cf: Vx vxX vxY vxPsi vxBx vxBy vxFfbz'
    write(fu, '(a)' ) tempstring 
    fmt = '('//Ifm// ',' //repeat(spacefm // Rfm // ',', 5)// spacefm // Rfm //')'
    do i = 1, nv 
        write(fu, fmt) i, xv(i), yv(i), psiv(i), Bxv(i), Byv(i), ffbzv(i)
    end do

    ! Cell information
    tempstring = '*cf: cv cvVxP(:,1) cvVxP(:,2) cvX cvY psi bp bt cflags(:) cvReg cvFt'
    write(fu, '(a)' ) tempstring 
    fmt = '('//repeat(Ifm //',' //spacefm, 3)//&
        repeat(Rfm //',' //spacefm, 5)// &
        repeat(Ifm //',' //spacefm, 3)//')'
    do i = 1, nc 
        write(fu, fmt) i, cellvertP(i, 1), cellvertP(i, 2), xc(i), yc(i), &
            psic(i), bpc(i), btc(i), flagc(i), regc(i), cellft(i)
    end do
    fmt = repeat(' ', 1000)
    write(fmt, *) '(', 12, '(', Ifm, '))'
    write(fu, '(2a8,i12,4x,a4)') '*cf:    ','int     ', ncellvert, 'cvVx'
    write(fu, fmt) cellvert

    fmt = repeat(' ', 1000)
    write(fmt, *) '(', 12, '(', Ifm, '))'
    write(fu, '(2a8,i12,4x,a4)') '*cf:    ','int     ', ncellface, 'cvFc'
    write(fu, fmt) cellface

    if (goatoptions%write_OMPdata) then 
        fmt = repeat(' ', 1000)
        write(fu, '(2a8,i12,4x,a5)') '*cf:    ','int     ', ncvOMP, 'cvOMP'
        write(fmt, *) '(', 12, '(', Ifm, '))'
        write(fu, fmt) cvOMP
        write(fu, '(2a8,i12,4x,a5)') '*cf:    ','int     ', ncvIMP, 'cvIMP'
        write(fmt, *) '(', 12, '(', Ifm, '))'
        write(fu, fmt) cvIMP
    end if 

    ! Face information
    tempstring = '*cf: fc fcVx(:,1) fcVx(:,2) fcLbl fcReg fcAligned'
    write(fu, '(a)' ) tempstring
    fmt = '('//repeat(Ifm //',' //spacefm, 6) //')'
    do i = 1, nf 
          write(fu, fmt) i, facevert(i, 1), facevert(i, 2), labelf(i), &
             regf(i), alignedf(i)
    end do

    if (goatoptions%write_OMPdata) then 
        tempstring = repeat(' ', 1000)
        write(fu, '(2a8,i12,4x,a5)') '*cf:    ','int     ', nfcOMP, 'fcOMP'
        fmt = repeat(' ', 1000)
        write(fmt, *) '(', 12, '(', Ifm, '))'
        write(fu, fmt) fcOMP
        write(fu, '(2a8,i12,4x,a5)') '*cf:    ','int     ', nfcIMP, 'fcIMP'
        fmt = repeat(' ', 1000)
        write(fmt, *) '(', 12, '(', Ifm, '))'
        write(fu, fmt) fcIMP
    end if

    ! Flux tube information
    tempstring = '*cf: ft ftCvP(:,1) ftCvP(:,2) ftFcP(:,1) ftFcP(:,2) ftReg'
    write(fu, '(a)' ) tempstring
    fmt = '('//repeat(Ifm,  6) //')'
    do i = 1, nft
          write(fu, fmt) i, ftcellP(i, 1), ftcellP(i, 2), ftfaceP(i, 1), &
             ftfaceP(i, 2), ftreg(i)
    end do

    tempstring = repeat(' ', 1000)
    write(fu, '(2a8,i12,4x,a4)') '*cf:    ','int     ', ftcellP(nft, 1)+ftcellP(nft, 2)-1, 'ftCv'
    fmt = repeat(' ', 1000)
    write(fmt, *) '(', 12, '(', Ifm, '))'
    write(fu, fmt) ftcell(1:ftcellP(nft, 1)+ftcellP(nft, 2)-1) 
    write(fu, '(2a8,i12,4x,a4)') '*cf:    ','int     ', ftfaceP(nft, 1)+ftfaceP(nft, 2)-1, 'ftFc'
    fmt = repeat(' ', 1000)
    write(fmt, *) '(', 12, '(', Ifm, '))'  
    write(fu, fmt) ftface(1:ftfaceP(nft, 1)+ftfaceP(nft, 2)-1)

    ! Flux surface information
    tempstring = '*cf: fs fsFcP(:,1) fsFcP(:,2) fsPsi'
    write(fu, '(a)' ) tempstring
    fmt = '('//repeat(Ifm, 4) //')'
    do i = 1, nfs
          write(fu, *) i, fsfaceP(i, 1), fsfaceP(i, 2), fspsi(i)
    end do

    tempstring = repeat(' ', 1000)
    write(fu, '(2a8,i12,4x,a4)') '*cf:    ','int     ', fsfaceP(nfs, 1)+fsfaceP(nfs, 2)-1, 'fsFc'
    fmt = repeat(' ', 1000)
    write(fmt, *) '(', 12, '(', Ifm, '))'  
    write(fu, fmt) fsface(1:fsfaceP(nfs, 1)+fsfaceP(nfs, 2)-1)

    ! Housekeeping
    !=============
    ! Close file
    close (fu)

    ! Write vessel
    !=============
    call environment%vessel%polygonset%WriteData(goatoptions%writefilepath // '_vesselpolygonset')

    ! End associate
    end associate


end subroutine