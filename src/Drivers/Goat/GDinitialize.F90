subroutine GDinitialize(inputfilepath, optimizationdriver)

    ! Description
    !============
    ! This routine does an initial read-in of the grid, magnetic field, 
    ! and environment variables. 

    ! Declare modules
    !================
    use goatmod_types
    use goatmod_userinput
    use gdmod_optimizationengine

    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    character(*), intent(in)                    :: inputfilepath 
    type(OptimizationEngineGDUDT), intent(out)  :: optimizationdriver
    type(GridUDT)                               :: grid 
    type(MagneticFieldUDT)                      :: magneticField 
    type(EnvironmentUDT)                        :: environment
    type(GoatoptionsUDT)                        :: goatoptions
    type(GDoptionsUDT)                          :: gdoptions
    type(DesignOptionsUDT)                      :: designoptions

    ! Read the user input
    !====================
    ! fileID should always be GOAToptions.dat
    goatoptions%inputfilepath = inputfilepath 
    print *, 'Reading goat options from file: ', goatoptions%inputfilepath
    call goatoptions%Set()

    ! Read and extract data
    call ExtractGoatData(grid, magneticField, environment, goatoptions)

    ! Set grid deformation options
    gdoptions%inputfilepath = goatoptions%inputfilepath
    call gdoptions%Set()

    ! Set additional options
    !=======================
    ! Set paths from where options should be read
    designoptions%inputfilepath = gdoptions%designoptionsfile

    ! Set optimization options
    call designoptions%Set()

    ! Initialize
    !===========
    ! Initialize the grid design problem (as an optimization problem)
    call ConstructGridDesignProblem(optimizationdriver, designoptions, &
        grid, magneticField, environment)

    ! Further initialize the design problem
    call optimizationdriver%problem%Initialize()


end subroutine