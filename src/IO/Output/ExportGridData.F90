subroutine ExportGridData(grid,options)
    ! Description
    !============
    ! This routine exports the data of a given grid to a file. The 
    ! export options should be set in the SetExportOptions subroutine 
    ! in the gdmod_userinput module. 

    ! Initialize
    !===========
    ! Declare modules
    use gdmod_types
    use gdmod_userinput

    ! The usual
    implicit none

    ! Declare variables
    type(GridUDT)              :: grid
    type(ExportOptionsUDT)     :: options
    character(128)             :: errmsg
    

    ! Export
    !=======
    ! Check the format
    select case (options%gridformat)

    case ('structured')

        ! Export the grid in a structured fashion. 
        print *, 'gd: exporting grid...'

    case ('unstructured')

        ! Not yet implemented
        errmsg = 'export option unstructured not &
        & yet implemented'
        call gdErrorHandler(errmsg)

    case default 
        
        ! Unknown grid format, throw error
        errmsg = 'unknown export grid format: ' // trim(options%gridformat)
        call gdErrorHandler(errmsg)
        ! print *, 'gd: unknown export grid format: ', options%gridformat

    end select



end subroutine
    