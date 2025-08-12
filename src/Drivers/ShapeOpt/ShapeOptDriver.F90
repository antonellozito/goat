subroutine ShapeOptDriver(inputfilepath)

    ! Description
    !============
    ! This driver runs a shape optimization driver that solves a shape
    ! optimization problem (see somod_optimizationengine.F90 for more
    ! details). Remeshing is supported using the goat grid generator.
    ! Note that a remeshing step means a full reinitialization of the 
    ! optimization problem, and is therefore a costly step. This is 
    ! necessary, since remeshing will change the number of state 
    ! and grid deformation equations. This means, however, that 
    ! remeshing in any subroutine of the optimization problem is 
    ! strictly speaking not allowed (or that it should be followed 
    ! by an optimization problem update, which can only be done at 
    ! the upper level right now)

    ! Initialize
    !===========
    ! Modules
    use somod_userinput
    use somod_optimizationengine
    use goatmod_types
    use ggmod_gridgenerator

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    character(*), intent(in)        :: inputfilepath 

    ! Auxiliary
    type(OptimizationEngineSOUDT)   :: ShapeOptEngine
    type(GridUDT)                   :: grid 
    type(GOAToptionsUDT)            :: goatoptions
    type(ShapeOptimizationOptionsUDT)   :: SOoptions
    type(GoatGridGenerator2DUDT)    :: gridgenerator 
    type(MagneticFieldUDT)          :: magneticField 
    type(EnvironmentUDT)            :: environment 
    type(GridOptionsUDT)            :: gridoptions 
    integer(I8), allocatable, dimension(:)  :: facelabelsGG, &
        facelabelsGD
    integer(I8), parameter                  :: emptyI8(0) = 0

    ! Loop
    integer(I8)                     :: itremesh 

    ! Initialize
    !===========
    ! Add the path 
    ShapeOptEngine%inputfilepath = inputfilepath
    ShapeOptEngine%inputfileprefix = 'so.'

    ! Initialize the driver
    call ShapeOptEngine%SetupOptimizationDriver()

    ! Initialize the problem
    ShapeOptEngine%problem%inputfilepath = inputfilepath
    !call ShapeOptEngine%problem%Initialize()

    ! Initialize options
    SOoptions%inputfilepath = inputfilepath 
    call SOoptions%Set()

    ! Initialize goat options (only for writing out data)
    goatoptions%inputfilepath = SOoptions%goatfilepath
    call goatoptions%Set()

    ! Read in grid generation data
    call ExtractGGData(magneticField, environment, goatoptions)

    ! Set grid options
    gridoptions%inputfilepath = goatoptions%inputfilepath 
    call gridoptions%Set()

    ! Initialize the grid generator
    call gridgenerator%Initialize(magneticField, environment, &
        SOoptions%goatfilepath)

    ! Solve
    !======
    ! Call driver for the first time
    call ShapeOptEngine%Driver()

    ! Check if we need to remesh and resolve
    associate(problem        => ShapeOptEngine%problem)

    select type (problem)

    type is (OptimizationProblemSOUDT)

        ! Initialize
        itremesh = 1

        ! Keep looping until exit conditions are met
        do while (SOoptions%doremesh .and. (itremesh <= SOoptions%itmaxremesh))

            ! Check if we need to remesh
            if (problem%doremesh) then 
                ! Reset remeshing
                problem%doremesh = .false. 

                ! Update remeshing iteration
                itremesh = itremesh + 1 
                
            else
                ! Converged or exited - stop
                exit 
            end if 

            ! Update the environment (mf not necessary)
            call gridgenerator%UpdateEnvironment(problem%goat%environment)

            ! Reconstruct the topomesh
            call gridgenerator%ConstructTopomesh()

            ! Reconstruct the grid
            call gridgenerator%ConstructGrid()

            ! Add to goat
            problem%goat%grid = gridgenerator%grid 

            ! Remap the face labels
            call GetGridFaceLabelMappingGD(gridgenerator%grid, &
                gridgenerator%topomesh, facelabelsGG, facelabelsGD)
            gridoptions%facelabelmappingGG = facelabelsGG 
            gridoptions%facelabelmappingGD = facelabelsGD 
            gridoptions%facelabelsubfrom   = emptyI8 
            gridoptions%facelabelsubto     = emptyI8 

            ! Recompute topological data from grid for new face labels
            call ComputeTopologicalData(problem%goat%grid, gridgenerator%topomesh)

            ! Extract the required grid data for grid deformation
            call ExtractGridData(problem%goat%grid, 'traduitb2us', gridoptions)

            ! Re-solve the optimization problem 
            call ShapeOptEngine%Driver()

        end do 

        ! Extract grid
        grid = problem%goat%grid 

        ! Write
        call WriteGOAT(goatoptions, grid, problem%goat%magneticField, &
            problem%goat%environment)

        ! b2ag file
        !call Writeb2agdat(goatoptions, grid)

    class default 
        
        call gdErrorHandler('ShapeOptDriver: unknown optimization problem type')

    end select

    ! Housekeeping
    !=============
    end associate

end subroutine