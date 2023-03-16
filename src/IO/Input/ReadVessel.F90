subroutine ReadVessel(filespecifier, vessel, vesseloptions)

    ! Description
    !============
    ! Read in the vessel data. 

    ! Notes
    !======

    ! Initialize
    !===========
    ! Declare modules
    use gdmod_types 
    use gdmod_userinput
    use gdmod_interfaces

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    integer                         :: filespecifier
    type(VesselOptionsUDT)          :: vesseloptions 
    type(VesselUDT)                 :: vessel

    ! Loop variables

    ! Auxiliary variables 

    ! Main program
    !=============
    ! Check how to read the vessel structure
    select case(vesseloptions%readmeth)
        
    case ('read_structure')

        ! Read in the separate vessel structures
        call read_structure(filespecifier, vessel, vesseloptions)

        ! Reformat the structures of the vessel into a single vessel
        ! polygon
        ! call FormatVesselStructures(vessel)

    case default

        ! Throw error
        call gdErrorHandler('Unknown vessel reading method')

    end select

end subroutine