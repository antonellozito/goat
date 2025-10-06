subroutine GAGDDriver(goatoptions)

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
    type(GDoptionsUDT)          :: gdoptions
    type(GridOptionsUDT)        :: gridoptions
    
    ! Auxiliary
    type(GAGridUDT)             :: GAgrid

    ! Initialize
    !===========
    ! Read and extract data
    call ExtractGoatData(grid, magneticField, environment, goatoptions)

    ! Set grid adaptation and deformation options
    gaoptions%inputfilepath         = goatoptions%inputfilepath
    gdoptions%inputfilepath         = goatoptions%inputfilepath
    gridoptions%inputfilepath       = goatoptions%inputfilepath 
    call gaoptions%Set()
    call gdoptions%Set()
    call gridoptions%Set()

    ! Grid adaptations
    !=================
    ! Carry-over options from goatoptions
    call CarryOverOptions(goatoptions, gaoptions)

    ! Translate Grid to GAGrid with dynamic arrays
    call TranslateGridTOGAGrid(grid,GAgrid)

    ! Run adaptations
    call GridAdaptor(GAgrid,environment,magneticField,gaoptions)

    ! Translate GAGrid to Grid
    call TranslateGAGridTOGrid(grid,GAgrid,gaoptions)

    ! Post-processing
    call PostProcessingGridInformation(grid,magneticField,gaoptions)

    ! Grid data
    call WriteGOAT(goatoptions, grid, magneticField, environment)

    ! Deallocate some grid properties
    if (allocated(grid%bnd)) deallocate(grid%bnd)

    ! Grid deformation
    !=================
    ! Extract the required grid data for grid deformation
    call ExtractGridData(grid, 'traduitb2us', gridoptions)

    ! Run deformation
    call RunGridOptimization(grid, magneticField, environment, &
        gdoptions)

    ! Grid data
    call WriteGOAT(goatoptions, grid, magneticField, environment)
    
end subroutine

