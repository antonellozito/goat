subroutine ConstructMagneticField(mfoptions, magneticField)

    ! Description
    !============
    ! Construct a representation of the magnetic field based on given
    ! magnetic field flux data (or other data). 

    ! Notes
    !======

    ! Initialize
    !===========
    ! Declare modules
    use gdmod_types 

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    type(MagneticFieldOptionsUDT)   :: mfoptions 
    type(MagneticFieldUDT)          :: magneticField

    ! Loop variables

    ! Auxiliary variables 
    integer                 :: filespecifier

    ! Data
    data filespecifier /60/

    ! Main program
    !=============
    ! Open the file containing the grid data
    call cfopen(filespecifier,'inputfiles/rzpsi','old','un*formatted')

    ! Read 
    call ReadMagneticField(filespecifier, mfoptions, magneticField)

    ! Construct interpolant representation
    !call Make2DStructuredInterpolant(magneticField%R, magneticField%Z, &
    !    magneticField%nR, magneticField%nZ, magneticField%interp)
    

end subroutine
