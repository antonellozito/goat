subroutine GGDriver(goatoptions)

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
    use ggmod_topology2D
    use ggmod_gridgeneration2D
    use mod_contour2D, only : ContourTracerUDT

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
    type(GGoptionsUDT)          :: ggoptions
    type(TopomeshOptionsUDT)    :: topomeshoptions
    
    ! Auxiliary
    type(TopomeshUDT)           :: topomesh
    class(ContourTracerUDT), allocatable    :: fieldtracer, vesseltracer

    ! Initialize
    !===========
    ! Read and extract data
    call ExtractGGData(magneticField, environment, goatoptions)

    ! Set topological mesh options
    topomeshoptions%inputfilepath = goatoptions%inputfilepath 
    call topomeshoptions%Set()

    ! Set grid generation options
    ggoptions%inputfilepath = goatoptions%inputfilepath 
    call ggoptions%Set()

    ! Generate the topological mesh
    !==============================
    call ConstructTopologicalMesh(environment%vessel, magneticField, &
        topomeshoptions, topomesh, fieldtracer, vesseltracer)

    ! Generate the grid
    !==================
    call GenerateUnstructuredAlignedGrid(topomesh, magneticField, &
        environment%vessel, fieldtracer, vesseltracer, ggoptions)

    ! Write data
    !===========
    ! Grid data
    ! call WriteGOAT(goatoptions, grid, magneticField, environment)

    ! b2ag file
    !call Writeb2agdat(goatoptions, grid)



end subroutine