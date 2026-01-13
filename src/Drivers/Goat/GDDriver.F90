subroutine GDDriver(goatoptions)

    ! Description
    !============
    ! This driver runs the grid deformation in standalone mode. The 
    ! goatoptions should be passed to this routine to identify which
    ! files to load etc. 

    ! Initialize
    !===========
    ! Modules
    use goatmod_types 
    use goatmod_userinput
    use ggmod_gridgeneration2D, only: WriteVoidRegionFile, WriteVoidRegionFileGoat, &
        ReadVoidRegionFileGoat, UpdateVoidRegionCoordinates

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    type(GoatoptionsUDT)        :: goatoptions 
    type(GridUDT)               :: grid 
    type(MagneticFieldUDT)      :: magneticField 
    type(EnvironmentUDT)        :: environment

    ! Other options
    type(GDoptionsUDT)          :: gdoptions
    
    ! Auxiliary
    logical                     :: fort78exists
    type(PolygonSEtUDT)         :: voidps

    ! Initialize
    !===========
    ! Read and extract data
    call ExtractGoatData(grid, magneticField, environment, goatoptions)

    ! Set grid deformation options
    gdoptions%inputfilepath = goatoptions%inputfilepath
    call gdoptions%Set()

    ! Run deformation
    !================
    call RunGridOptimization(grid, magneticField, environment, &
        gdoptions)

    ! Write data
    !===========
    ! Grid data
    call WriteGOAT(goatoptions, grid, magneticField, environment)

    ! Void polygon set (if available)
    inquire(file='fort_goat.78', exist=fort78exists) 
    if (fort78exists) then 
        ! Display message
        print *, 'GDDriver: fort_goat.78 file found, reading and adjusting ' // & 
            'vertex coordinates'

        ! Read void polygon set
        call ReadVoidRegionFileGoat(voidps, 'fort_goat.78')

        ! Update coordinates of polygons
        call UpdateVoidRegionCoordinates(voidps, grid)

        ! Write void polygon set
        print *, 'GDDriver: writing new void polygon to fort.78'
        call WriteVoidRegionFile(voidps, grid, 'fort.78')
        print *, 'GDDriver: writing new void polygon to fort_goat.78'
        call WriteVoidRegionFileGoat(voidps, grid, 'fort_goat.78')
    else
        ! Display message
        print *, 'GDDriver: fort_goat.78 file not found, continuing...'
    end if

    ! b2ag file
    !call Writeb2agdat(goatoptions, grid)

    






end subroutine