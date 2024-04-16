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
    call WriteGOAT(goatoptions, grid)

    ! b2ag file
    !call Writeb2agdat(goatoptions, grid)

    






end subroutine