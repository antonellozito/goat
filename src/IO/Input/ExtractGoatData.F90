subroutine ExtractGoatData(grid, magneticField, environment, options)

    ! Description
    !============
    ! Extract all necessary data to run the goat drivers. This includes
    ! the grid, magnetic field, and environment. The options structure
    ! should contain the necessary paths to read in these quantities. 

    ! Initialize
    !===========
    ! Modules
    use goatmod_types
    use goatmod_userinput 

    ! Declare variables
    !==================
    ! Arguments
    type(GridUDT)                       :: grid 
    type(MagneticFieldUDT)              :: magneticField 
    type(EnvironmentUDT)                :: environment 
    type(GoatoptionsUDT)                :: options 

    ! Auxiliary 
    type(GridOptionsUDT)                :: gridoptions
    type(MagneticFieldOptionsUDT)       :: mfoptions 
    type(EnvironmentOptionsUDT)         :: environmentoptions 

    ! Set options
    !============
    ! Set the path - assumed to be already set correctly in options
    gridoptions%inputfilepath           = options%gridoptionsfilepath 
    mfoptions%inputfilepath             = options%magneticfieldoptionsfilepath
    environmentoptions%inputfilepath    = options%environmentoptionsfilepath

    ! Set the options
    call gridoptions%Set()
    call mfoptions%Set() 
    call environmentoptions%Set()

    ! Read data
    !==========
    call ConstructGrid(grid, gridoptions)
    call ConstructMagneticField(magneticField, mfoptions) 
    call ConstructEnvironment(environment, environmentoptions) 


end subroutine