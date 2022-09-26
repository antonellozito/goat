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
    call cfopen(filespecifier,'inputfiles/rzpsi_asdex.dat','old','un*formatted')

    ! Read 
    call ReadMagneticField(filespecifier, mfoptions, magneticField)

    ! Construct interpolant representation
    call ConstructBicubicSplineInterpolant(magneticField%Psi, &
        magneticField%R, magneticField%Z, magneticField%nR, &
        magneticField%nZ, interp)

    ! Add interpolant to the magnetic field
    magneticField%interp = interp

end subroutine
