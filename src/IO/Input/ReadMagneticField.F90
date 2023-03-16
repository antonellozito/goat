subroutine ReadMagneticField(filespecifier, magneticField, mfoptions)

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
    use gdmod_types 
    use gdmod_userinput

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    integer                         :: filespecifier
    type(MagneticFieldOptionsUDT)   :: mfoptions 
    type(MagneticFieldUDT)          :: magneticField

    ! Loop variables

    ! Auxiliary variables 


    ! Main program
    !=============
    ! Check how to read the file
    select case(mfoptions%readmeth)

    case ('readrzpsi')

        ! Read the rzpsi file using CARRE subroutines - wrapper here
        call Readrzpsi(filespecifier, magneticField)

    case default

        ! Unknown reading method, throw error
        call gdErrorHandler('ReadMagneticField: unknown reading method')

    end select


end subroutine