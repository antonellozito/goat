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
    type(GAGridUDT)             :: GAgrid

    ! Initialize
    !===========
    ! Read and extract data
    call ExtractGoatData(grid, magneticField, environment, goatoptions)

    ! Set grid adaptation options
    gaoptions%inputfilepath         = goatoptions%inputfilepath
    call gaoptions%Set()

    ! Override options from goatoptions
    gaoptions%vesselmode            = goatoptions%vesselmode 
    gaoptions%slab                  = goatoptions%slab
    gaoptions%debug                 = goatoptions%debug 
    gaoptions%facelabelmappingGG    = goatoptions%GGtoGAfacelabelmappingGG
    gaoptions%facelabelmappingGA    = goatoptions%GGtoGAfacelabelmappingGA
    gaoptions%OMPr                  = goatoptions%OMPr
    gaoptions%OMPz                  = goatoptions%OMPz
    gaoptions%IMPr                  = goatoptions%IMPr
    gaoptions%IMPz                  = goatoptions%IMPz

    ! Translate Grid to GAGrid with dynamic arrays
    call TranslateGridTOGAGrid(grid,GAgrid)

    ! Run adaptations
    !================
    call GridAdaptor(GAgrid,environment,magneticField,gaoptions)

    ! Translate GAGrid to Grid
    call TranslateGAGridTOGrid(grid,GAgrid,gaoptions)

    
    ! Post-processing
    !================
    call PostProcessingGridInformation(grid,magneticField,gaoptions)

    ! Write data
    !===========
    ! Grid data
    call WriteGOAT(goatoptions, grid, magneticField, environment)

end subroutine

