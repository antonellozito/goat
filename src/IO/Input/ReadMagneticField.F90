subroutine ReadMagneticField(magneticField, mfoptions, filepath)

    ! Description
    !============
    ! Read in the magnetic field data. It is assumed that the data is 
    ! given as flux values on a structured (but not necessarily
    ! equidistant) grid. 

    ! Notes
    !======

    ! Initialize
    !===========
    ! Declare modules
    use goatmod_types 
    use goatmod_userinput
    
    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    integer                         :: filespecifier
    type(MagneticFieldOptionsUDT)   :: mfoptions 
    type(MagneticFieldUDT)          :: magneticField
    character(*)                    :: filepath

    ! Loop variables

    ! Auxiliary variables 

    ! Data
    data filespecifier /60/

    ! Main program
    !=============
    ! Open the file
    print *, 'reading magnetic field data from file: ' // filepath
    open(unit = filespecifier, file = filepath)

    ! Check how to read the file
    select case(mfoptions%readmeth)

    case ('rzpsi','default')

        ! Read in a classic rzpsi file 
        call Readrzpsi(filespecifier, magneticField)

    case ('equ')

        ! Read the equ file using CARRE subroutines - wrapper here
        call Readequ(filespecifier, magneticField)

    case default

        ! Unknown reading method, throw error
        call gdErrorHandler('ReadMagneticField: unknown reading method')

    end select

    ! Add additional data
    magneticField%RBtor = mfoptions%RBtor

    ! Housekeeping
    !=============
    close(filespecifier)


end subroutine