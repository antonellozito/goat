subroutine ConstructMagneticField(magneticField, mfoptions)

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
    use gdmod_userinput
    use BicubicSplineInterpolant

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    type(MagneticFieldOptionsUDT)       :: mfoptions 
    type(MagneticFieldUDT)              :: magneticField

    ! Loop variables

    ! Auxiliary variables 
    integer                             :: filespecifier
    type(BicubicSplineInterpolantUDT)   :: interp

    ! Data
    data filespecifier /60/

    ! Main program
    !=============
    ! Open the file containing the grid data
    call cfopen(filespecifier, mfoptions%filepath,'old','un*formatted')

    ! Read 
    call ReadMagneticField(filespecifier, magneticField, mfoptions)

    ! Construct interpolant representation
    call ConstructBicubicSplineInterpolant(magneticField%Psi, &
        magneticField%R, magneticField%Z, magneticField%nR, &
        magneticField%nZ, interp)

    ! Add interpolant to the magnetic field
    magneticField%interp = interp

end subroutine
