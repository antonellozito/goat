
subroutine ExtractGGData(magneticField, environment, options)

    ! Description
    !============
    ! Extract all necessary data to run the grid generator driver. This includes
    ! the magnetic field and environment. The options structure
    ! should contain the necessary paths to read in these quantities. 

    ! Initialize
    !===========
    ! Modules
    use goatmod_types
    use goatmod_userinput 

    ! Declare variables
    !==================
    ! Arguments
    type(MagneticFieldUDT)              :: magneticField 
    type(EnvironmentUDT)                :: environment 
    type(GoatoptionsUDT)                :: options 

    ! Auxiliary 
    type(MagneticFieldOptionsUDT)       :: mfoptions 
    type(EnvironmentOptionsUDT)         :: environmentoptions 

    integer                 :: filespecifier(0:2)
    
    ! Declare data
    data filespecifier /50, 51, 52/ ! numbers for file specification

    ! Initialize
    !===========
    ! Set the path - the same as for general goat options 
    environmentoptions%inputfilepath    = options%inputfilepath
    mfoptions%inputfilepath             = options%inputfilepath 
    
    ! Set the options
    call environmentoptions%Set()
    call mfoptions%Set()

    ! Reset vessel reading 
    environmentoptions%vesselfilepath = options%structurefilepath

    ! Read data
    !==========
    ! Read magnetic field
    call ReadMagneticField(magneticField, mfoptions, options%magneticfieldfilepath)

    ! Read additional data
    !=====================
    call ConstructMagneticField(magneticField, mfoptions) 
    call ConstructEnvironment(environment, environmentoptions) 

    ! Housekeeping
    !=============


end subroutine