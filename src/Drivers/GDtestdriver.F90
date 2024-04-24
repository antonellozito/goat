subroutine GDtestdriver(goatoptions)

    ! Initialize
    !===========
    ! Modules
    use goatmod_types 
    use goatmod_userinput
    use gdmod_optimizationengine
    use mod_constants

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    type(GoatoptionsUDT), intent(inout)     :: goatoptions 
    type(OptimizationEngineGDUDT)           :: optimizationdriver
    
    ! Auxiliary
    real(R8)                        :: dist, angle
    real(R8), allocatable           :: xv(:), yv(:), dx(:), dy(:)

    integer(I8)                     :: nv, incr 

    character(:), allocatable       :: orig_writefilepath
    character(len=1024)             :: tempstring

    ! Loop
    integer(I8)                     :: i, cc

    ! Initialize
    !===========
    ! Set distance to move vessel vertices
    dist = 1e-1
    angle = 0.0*pi_R8/4
    incr = 100

    ! Set up optimization problem
    call GDinitialize(goatoptions%inputfilepath, optimizationdriver)

    ! Get the initial vessel coordinates
    associate(problem       => optimizationdriver%problem)
    select type(problem)

    type is (OptimizationProblemGDUDT)

        ! Associate for ease
        associate(vessel        => problem%environment%vessel)

        ! Compute total number of vertices
        nv = 0
        do i = 1, vessel%polygonset%np
            nv = nv + vessel%polygonset%polygons(i)%nv
        end do

        ! Extract
        cc = 0
        allocate(xv(nv), yv(nv), dx(nv), dy(nv))
        do i = 1, vessel%polygonset%np 
            xv(cc+1:cc+vessel%polygonset%polygons(i)%nv) = vessel%polygonset%polygons(i)%x
            yv(cc+1:cc+vessel%polygonset%polygons(i)%nv) = vessel%polygonset%polygons(i)%y
        end do

        ! Housekeeping
        end associate

    class default
    end select 
    end associate

    ! Solve for the initial geometry
    !===============================
    ! Solve
    call optimizationdriver%Driver()

    ! Loop over vessel geometry
    !==========================
    orig_writefilepath = goatoptions%writefilepath
    associate(problem       => optimizationdriver%problem)
    select type(problem)

    type is (OptimizationProblemGDUDT)

        do i = 1, incr 
            ! Print
            print *, '================================================='
            print *, '          vessel geometry iteration ', i 
            print *, '================================================='
            ! Update geometry
            dx = dist*cos(angle)*(real(1, kind=R8)/real(incr, kind=R8))
            dy = dist*sin(angle)*(real(1, kind=R8)/real(incr, kind=R8))
            xv = xv + dx 
            yv = yv + dy

            ! Update vessel description
            call problem%UpdateProblemParameters([xv, yv], 'vesselcoordinates')

            ! Solve for new geomeatry
            call optimizationdriver%solver%SolveOptimizationProblemKKT(problem)

            ! Write solution
            write(tempstring, '(I0.4, a)') i, trim(orig_writefilepath) 
            goatoptions%writefilepath = trim(tempstring)
            call WriteGOAT(goatoptions, problem%grid)

        end do

    end select  
    end associate

    ! b2ag file
    !call Writeb2agdat(goatoptions, grid)




end subroutine