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
    type(GDoptionsUDT)                  :: gdoptions

    ! Auxiliary 
    type(GridOptionsUDT)                :: gridoptions
    type(MagneticFieldOptionsUDT)       :: mfoptions 
    type(EnvironmentOptionsUDT)         :: environmentoptions 

    integer                 :: filespecifier(0:2)
    
    ! Declare data
    data filespecifier /50, 51, 52/ ! numbers for file specification

    ! Initialize
    !===========
    ! Set the path - the same as for general goat options 
    gdoptions%inputfilepath             = options%inputfilepath
    gridoptions%inputfilepath           = options%inputfilepath 
    environmentoptions%inputfilepath    = options%inputfilepath
    mfoptions%inputfilepath             = options%magneticfieldfilepath 
    
    ! Set the options
    call gdoptions%Set()
    call gridoptions%Set()
    call environmentoptions%Set()
    call mfoptions%Set()

    ! Set input filepaths of other options
    gdoptions%inputfilepath = options%gdinputfilepath

    ! Read data
    !==========
    ! Read grid
    select case (options%gridreadtype)

    case ('traduitb2us')

        ! Unstructured traduit file reading
        print *, options%gridfilepath
        call ReadTraduitUS(grid, options%gridfilepath)


    case default

        call gdErrorHandler('ExtractGoatData: unknown readtype option')

    end select 

    ! Read magnetic field
    call ReadMagneticField(magneticField, mfoptions, options%magneticfieldfilepath)

    ! Read additional data
    !=====================
    call ConstructGrid(grid, gridoptions)
    call ConstructMagneticField(magneticField, mfoptions) 
    call ConstructEnvironment(environment, environmentoptions) 

    ! Housekeeping
    !=============


end subroutine