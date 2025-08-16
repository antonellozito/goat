subroutine GADriver(goatoptions)

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
    use gamod_driver

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
    type(GAoptionsUDT)          :: gaoptions
    
    ! Auxiliary

    ! Initialize
    !===========
    ! Read and extract data
    call ExtractGoatData(grid, magneticField, environment, goatoptions)

    ! Set grid adaptation options
    gaoptions%inputfilepath = goatoptions%inputfilepath
    call gaoptions%Set()

    ! Run adaptations
    !================
    call GridAdaptor(grid,environment,magneticField,gaoptions)


    ! Post-processing
    !================
    ! call PostProcessingGridInformation(grid,magneticField,goatoptions)

    ! Write data
    !===========
    ! Grid data
    call WriteGOAT(goatoptions, grid, magneticField, environment)

end subroutine

