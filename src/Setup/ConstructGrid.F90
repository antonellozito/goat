subroutine ConstructGrid(grid, gridoptions)
    
    ! Description
    !============
    ! This function further processes the grid and constructs additional 
    ! data needed for most GOAT drivers. 

    ! Initialize
    !===========
    ! Declare modules
    use goatmod_types
    use goatmod_userinput

    ! The usual
    implicit none

    ! Arguments
    type(GridUDT)                           :: grid
    type(GridOptionsUDT)                    :: gridoptions

    ! Construct grid
    !===============
    select case (gridoptions%type)

    case ('plasma')

        ! Default plasma edge grid, check the provided format
        select case (gridoptions%readmeth)

        case ('b2fgmtry')

            ! Extract the grid data structures
            call ExtractGridData(grid, 'b2fgmtry', gridoptions)

        case ('traduitb2us')

            ! Extract the grid data structures 
            call ExtractGridData(grid, 'traduit', gridoptions)

        case default 

            call gdErrorHandler('unknown plasma grid input type')

        end select


    case default

        call gdErrorHandler('unknown grid type')

    end select

    ! Compute grid interconnections
    !==============================
    call ComputeGridInterconnections(grid)


end subroutine