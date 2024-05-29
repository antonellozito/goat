subroutine ShapeOptDriver(inputfilepath)

    ! Description
    !============
    ! This driver runs a shape optimization driver that solves a shape
    ! optimization problem (see somod_optimizationengine.F90 for more
    ! details).

    ! Initialize
    !===========
    ! Modules
    use somod_userinput
    use somod_optimizationengine
    use goatmod_types

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    character(:), allocatable       :: inputfilepath 

    ! Auxiliary
    type(OptimizationEngineSOUDT)   :: ShapeOptEngine
    type(GridUDT)                   :: grid
    type(GOAToptionsUDT)            :: goatoptions
    type(ShapeOptimizationOptionsUDT)   :: SOoptions

    ! Initialize
    !===========
    ! Initialize problem
    ShapeOptEngine%problem%inputfilepath = inputfilepath
    call ShapeOptEngine%problem%Initialize()

    ! Initialize options
    SOoptions%inputfilepath = inputfilepath 
    call SOoptions%Set()

    ! Initialize goat options (only for writing out data)
    goatoptions%inputfilepath = SOoptions%goatfilepath
    call goatoptions%Set()

    ! Solve
    !======
    ! Solve the problem
    call ShapeOptEngine%solver%SolveOptimizationProblem(ShapeOptEngine%problem)

    ! Write data
    !===========
    ! Extract
    associate(problem        => ShapeOptEngine%problem)

    select type (problem)

    type is (OptimizationProblemSOUDT)

        ! Extract grid
        grid = problem%goat%grid 

        ! Write
        call WriteGOAT(goatoptions, grid)

        ! b2ag file
        !call Writeb2agdat(goatoptions, grid)

    class default 
        
        call gdErrorHandler('ShapeOptDriver: unknown optimization problem type')

    end select

    ! Housekeeping
    !=============
    end associate

end subroutine