subroutine GDinitialize(inputfilepath, grid, magneticField, &
    environment, goatoptions, gdoptions)

    ! Description
    !============
    ! This routine does an initial read-in of the grid, magnetic field, 
    ! and environment variables. 

    ! Declare modules
    !================
    use goatmod_types
    use goatmod_userinput

    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    character(*), intent(in)                :: inputfilepath 
    type(GridUDT), intent(out)              :: grid 
    type(MagneticFieldUDT), intent(out)     :: magneticField 
    type(EnvironmentUDT), intent(out)       :: environment
    type(GoatoptionsUDT), intent(out)       :: goatoptions
    type(GDoptionsUDT), intent(out)         :: gdoptions

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


end subroutine