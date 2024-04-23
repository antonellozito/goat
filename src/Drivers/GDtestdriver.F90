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
    real(R8), allocatable           :: xv(:), yv(:)

    integer(I8)                     :: nv, incr 

    character(:), allocatable       :: orig_writefilepath
    character(len=1024)             :: tempstring

    ! Loop
    integer(I8)                     :: i, cc

    ! Initialize
    !===========
    ! Set distance to move vessel vertices
    dist = 1e-2
    angle = pi_R8/4
    incr = 10

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
        allocate(xv(nv), yv(nv))
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
            ! Update geometry
            xv = xv + dist*cos(angle)*(1/incr)
            yv = yv + dist*sin(angle)*(1/incr)

            ! Update vessel description
            call problem%UpdateProblemParameters([xv, yv], 'vesselcoordinates')

            ! Solve for new geomeatry
            call optimizationdriver%solver%SolveOptimizationProblemKKT(optimizationdriver%problem)

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