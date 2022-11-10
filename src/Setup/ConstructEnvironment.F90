subroutine ConstructEnvironment(environment, environmentoptions)

    ! Description
    !============
    ! Construct the environment for the grid optimization. Here, all 
    ! other peripheral structures, vessel data, ... should be stored
    ! that may be needed in the grid optimization (apart from the 
    ! magnetic field, which, due to it's main role in the grid 
    ! generation, has it's own structure). It goes without saying that 
    ! this routine is heavily case dependent. 

    ! Notes
    !======
    ! Right now, only the vessel is added.

    ! Initialize
    !===========
    ! Declare modules
    use gdmod_types 
    use gdmod_userinput
    use gdmod_plots

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    type(EnvironmentOptionsUDT)       :: environmentoptions 
    type(EnvironmentUDT)              :: environment

    ! Loop variables

    ! Auxiliary variables 
    integer                             :: filespecifier
    logical                             :: debugplots

    ! Additional environment structures
    type(VesselOptionsUDT)              :: vesseloptions

    ! Data
    data filespecifier /60/
    data debugplots /.false./

    ! Initialize
    !===========
    ! Check which type of environment to construct
    select case (trim(environmentoptions%type))

    case ('vessel')

        ! Read in vessel options
        call SetVesselOptions(vesseloptions)

        ! Read in the vessel structure
        call ReadVessel(filespecifier, environment%vessel, vesseloptions)

        ! Extract the vessel data
        call ExtractVesselData(environment%vessel, vesseloptions)

        ! Visualize?
        if (debugplots) then 
            ! Visualize
            call PlotVesselPolygon(environment%vessel, '-p')
        end if

    case default

        ! Unknown case, throw error
        call gdErrorHandler('Unknown environment type')

    end select

end subroutine