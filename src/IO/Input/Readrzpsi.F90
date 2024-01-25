subroutine Readrzpsi(filespecifier, magneticField)

    ! Description
    !============
    ! Read in the rzpsi data and add those to the magnetic field. It is 
    ! assumed that the file first gives the r coordinates, then z 
    ! coordinates, and finally the psi values in a 2D structured way. 
    ! The following format is assumed: 
    ! - first, we should encounter 'nr=xxx', where 'xxx' is an integer
    ! that specified the number of R-coordinates
    ! - then, the R-coordinates are specified as reals

    ! Notes
    !======

    ! Initialize
    !===========
    ! Declare modules
    use goatmod_types 
    use goatmod_userinput
    use mod_inputfileparser
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
    real(R8), allocatable       :: R(:), Z(:), psi(:, :), psivec(:), &
        valr(:)

    integer(I8)                 :: nr, nz ! real number of points in r, z directions
    integer(I8)                 :: n, k
    integer(I8), allocatable    :: val(:)
    
    character(:), allocatable   :: thisline ! temporary variable for line
    logical                     :: iseof = .false. 

    ! Debug
    logical                     :: makedebugplots = .true.

    ! Read R-coordinates
    !===================
    ! Read until we encounter 'nr'
    do while (.true.)
        ! Read in the next line
        call ReadSingleLine(filespecifier, thisline, iseof)

        ! Check if we reached the end of the file
        if (iseof) then 
            ! Throw error
            call gdErrorHandler('Readrzpsi: could not find r-coordinates in rzpsi file')

            ! Exit the loop
            exit 
        end if 

        ! Check if we encounter the specified string
        if (index(thisline, 'nr') .ne. 0) then 
            ! Found, extract nr and break the while loop
            call ReadIntegersFromString(thisline, val, n)

            ! Sanity check: only one value can be found
            if (n .ne. 1) then 
                call gdErrorHandler('Readrzpsi: found nr, but could not extract value')
            end if

            ! Set nr & exit
            nr = val(1)
            exit 
        end if
    end do

    ! Allocate
    allocate(R(nr))

    ! Read the R-coordinates line by line 
    k = 0
    do while (.true.)
        ! Read next line
        call ReadSingleLine(filespecifier, thisline, iseof)

        ! Check if we reached the end of the file prematurely
        if (iseof) then 
            ! Throw error
            call gdErrorHandler('Readrzpsi: could not find r-coordinates in rzpsi file')

            ! Exit the loop
            exit 
        end if 

        ! Extract reals
        call ReadRealsFromString(thisline, valr, n)

        ! Check 
        if ((k + n) > nr) then 
            ! This shouldn't be happening, throw error
            call gdErrorHandler('Readrzpsi: encountered more r-coordinate entries than given by nr')
        end if 

        ! Add 
        R(k+1:k+n) = valr 

        ! Update k
        k = k + n

        ! Check exit conditions
        if ( k == nr ) then 
            ! Normal exit
            exit
        elseif (index(thisline, 'z') .ne. 0) then 
            ! Abnormal exit
            call gdErrorHandler('Readrzpsi: encountered z-coordinates prematurely, probably less r-coordinate entries than nr')
        end if 
    end do 

    ! Read Z-coordinates
    !===================
    ! Read until we encounter 'nz'
    do while (.true.)
        ! Read in the next line
        call ReadSingleLine(filespecifier, thisline, iseof)

        ! Check if we reached the end of the file
        if (iseof) then 
            ! Throw error
            call gdErrorHandler('Readrzpsi: could not find r-coordinates in rzpsi file')

            ! Exit the loop
            exit 
        end if 

        ! Check if we encounter the specified string
        if (index(thisline, 'nz') .ne. 0) then 
            ! Found, extract nz and break the while loop
            call ReadIntegersFromString(thisline, val, n)

            ! Sanity check: only one value can be found
            if (n .ne. 1) then 
                call gdErrorHandler('Readrzpsi: found nz, but could not extract value')
            end if

            ! Set nz & exit
            nz = val(1)
            exit 
        end if
    end do

    ! Allocate
    allocate(Z(nz))

    ! Read the Z-coordinates line by line 
    k = 0
    do while (.true.)
        ! Read next line
        call ReadSingleLine(filespecifier, thisline, iseof)

        ! Check if we reached the end of the file prematurely
        if (iseof) then 
            ! Throw error
            call gdErrorHandler('Readrzpsi: could not find z-coordinates in rzpsi file')

            ! Exit the loop
            exit 
        end if 

        ! Extract reals
        call ReadRealsFromString(thisline, valr, n)

        ! Check 
        if ((k + n) > nz) then 
            ! This shouldn't be happening, throw error
            call gdErrorHandler('Readrzpsi: encountered more z-coordinate entries than given by nz')
        end if 

        ! Add 
        Z(k+1:k+n) = valr 

        ! Update k
        k = k + n

        ! Check exit conditions
        if ( k == nz ) then 
            exit
        elseif (index(thisline, 'psi') .ne. 0) then 
            ! Abnormal exit
            call gdErrorHandler('Readrzpsi: encountered psi values prematurely, probably less z-coordinate entries than nz')
        end if 
    end do 

    ! Read Psi values
    !================
    ! Read until we encounter 'psi'
    do while (.true.)
        ! Read in the next line
        call ReadSingleLine(filespecifier, thisline, iseof)

        ! Check if we reached the end of the file
        if (iseof) then 
            ! Throw error
            call gdErrorHandler('Readrzpsi: could not find psi values in rzpsi file')

            ! Exit the loop
            exit 
        end if 

        ! Check if we encounter the specified string
        if (index(thisline, 'psi') .ne. 0) then 
            ! Found, exit while loop
            exit 
        end if
    end do

    ! Allocate
    allocate(psi(nr, nz))
    allocate(psivec(nr*nz))

    ! Read the psi values line by line 
    k = 0
    do while (.true.)
        ! Read next line
        call ReadSingleLine(filespecifier, thisline, iseof)

        ! Check if we reached the end of the file prematurely
        if (iseof) then 
            ! Throw error
            call gdErrorHandler('Readrzpsi: could not find psi values in rzpsi file')

            ! Exit the loop
            exit 
        end if 

        ! Extract reals
        call ReadRealsFromString(thisline, valr, n)

        ! Check 
        if ((k + n) > nr*nz) then 
            ! This shouldn't be happening, throw error
            call gdErrorHandler('Readrzpsi: encountered more psi value entries than given by nr*nz')
        end if 

        ! Add 
        psivec(k+1:k+n) = valr 

        ! Update k
        k = k + n

        ! Check exit conditions
        if ( k == nr*nz ) then 
            exit
        end if 
    end do 

    ! Reshape psivec into psi
    psi = reshape(psivec, (/nr, nz/))

    ! Add to magnetic field
    !======================
    ! Allocate the magnetic field 
    magneticField%nR = nr
    magneticField%nZ = nz
    call AllocateMagneticField(magneticField)

    ! Add data
    magneticField%R = R
    magneticField%Z = Z
    magneticField%Psi = psi

    ! Make plots
    if (makedebugplots) then
        call PlotMagneticFlux(magneticField, '-p')
    end if

end subroutine