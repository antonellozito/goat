subroutine ConstructMagneticField(magneticField, mfoptions)

    ! Description
    !============
    ! Construct a representation of the magnetic field based on given
    ! magnetic field flux data (or other data). It is assumed that the 
    ! basic magnetic field data is already available in the 
    ! magneticField structure

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

    ! Main program
    !=============
    ! Construct interpolant representation
    call ConstructBicubicSplineInterpolant(magneticField%Psi, &
        magneticField%R, magneticField%Z, magneticField%nR, &
        magneticField%nZ, interp)

    ! Add interpolant to the magnetic field
    magneticField%interp = interp

end subroutine
