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

    ! Run deformation
    !================
    call RunGridOptimization(grid, magneticField, environment, &
        gdoptions)


    






end subroutine