subroutine Readequ(filespecifier, magneticField)

    ! Description
    !============
    ! Read in the equ data and add those to the magnetic field. This 
    ! routine serves as a wrapper for the CARRE reading files rdeqdg.F
    ! and rdeqlh.F. 

    ! Notes
    !======

    ! Initialize
    !===========
    ! Declare modules
    use gdmod_types 
    use gdmod_plots

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    integer                         :: filespecifier
    type(MagneticFieldUDT)          :: magneticField

    ! Loop variables

    ! Auxiliary variables 
    real(R8)  , allocatable     :: R(:), Z(:) ! coordinates of magnetic field grid
    real(R8)                    :: btf, rtf ! toroidal field strength and radius
    integer(I8)                 :: returncode ! return code from rdeqdg
    integer(I8)                 :: nr, nz ! real number of points in r, z directions
    integer(I8)                 :: maxnr, maxnz ! maximal number of points in r, z direction
    real(R8), allocatable       :: pfm(:, :) ! temporary psi value array

    ! Debug
    logical                     :: makedebugplots = .false.

    ! Data
    data maxnr /4000/
    data maxnz /4000/

    ! Main program
    !=============
    ! Allocate
    allocate(pfm(maxnr, maxnz))
    allocate(R(maxnr))
    allocate(Z(maxnz))

    ! Call file reader rdeqdg
    call rdeqdg(filespecifier, maxnr, maxnz, returncode, nr, nz, btf, rtf, &
            R, Z, pfm)

    ! Allocate the magnetic field 
    magneticField%nR = nr
    magneticField%nZ = nz
    call AllocateMagneticField(magneticField)

    ! Add data
    magneticField%R = R(1:nr)
    magneticField%Z = Z(1:nz)
    magneticField%Psi = pfm(1:nr, 1:nz)

    ! Make plots
    if (makedebugplots) then
        call PlotMagneticFlux(magneticField, '-p')
    end if

end subroutine