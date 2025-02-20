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
    use Interpolant

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    type(MagneticFieldOptionsUDT)       :: mfoptions 
    type(MagneticFieldUDT)              :: magneticField

    ! Auxiliary 
    real(R8)                                :: dR, dZ 
    real(R8), allocatable, dimension(:)     :: dx, dy, x, y, temp
    real(R8), allocatable, dimension(: ,:)  :: val, tempval  
    type(StructuredInterpolant2DUDT)    :: interp

    ! Loop
    integer(I8)                         :: i, k 


    ! Main program
    !=============
    ! Check if we want to reinterpolate on a different resolution
    if (mfoptions%reinterpolate) then 
        ! Compute field sizes
        dR = magneticField%R(magneticField%nR) - magneticField%R(1)
        dZ = magneticField%Z(magneticField%nZ) - magneticField%Z(1)

        ! Compute reinterpolation vectors 
        x = [(k, k = 0, mfoptions%resx)]*dR/real(mfoptions%resx, kind=R8) + magneticField%R(1)
        x(1) = magneticField%R(1) 
        x(size(x)) = magneticField%R(magneticField%nR)

        y = [(k, k = 0, mfoptions%resy)]*dZ/real(mfoptions%resy, kind=R8) + magneticField%Z(1)
        y(1) = magneticField%Z(1) 
        y(size(x)) = magneticField%Z(magneticField%nZ)

        ! Interpolate
        val = magneticField%Psi 
        allocate(tempval(size(x), magneticField%nZ))
        do i = 1, magneticField%nZ 
            call Interpolate1D(x, temp, magneticField%R, val(:, i))
            tempval(:, i) = temp
        end do 
        val = tempval 
        deallocate(tempval)
        allocate(tempval(size(x), size(y)))
        do i = 1, size(x) 
            call Interpolate1D(y, temp, magneticField%Z, val(i, :))
            tempval(i, :) = temp
        end do 
        val = tempval 
    else
        val = magneticField%Psi 
        x = magneticField%R 
        y = magneticField%Z 
    end if 

    ! Construct interpolant representation
    call interp%SetParameters(mfoptions%interpmeth, mfoptions%interpC, &
        mfoptions%interpM)
    call interp%Construct(x, y, val)

    ! Add interpolant to the magnetic field
    magneticField%interp = interp

    ! Visualize
    call magneticField%interp%Visualize('magneticfield_visualization')

end subroutine
