subroutine RunGridOptimization(grid, optimizationdriver, options)
    ! Description
    !============
    ! Main driver for optimization-based grid deformation. This routine
    ! is a wrapper where, depending on the options parsed in the input, 
    ! the initialization is done and the actual optimization driver is
    ! called. 

    ! Initialize
    !===========
    ! Declare modules
    use gdmod_types
    use gdmod_userinput 
    use gdmod_plots
    use mod_plotter
    use BicubicSplineInterpolant
    use gdmod_optimizationengine

    ! The usual
    implicit none

    ! Declare variables
    type(GridUDT)                   :: grid
    type(MagneticFieldUDT)          :: magneticField
    type(OptimizationEngineGDUDT)   :: optimizationdriver
    type(RunfileOptionsUDT)         :: options 
    type(GridOptionsUDT)            :: gridoptions
    type(DesignOptionsUDT)          :: designoptions
    type(MagneticFieldOptionsUDT)   :: mfoptions
    type(NumOptionsUDT)             :: num
    type(EnvironmentOptionsUDT)     :: environmentoptions
    type(EnvironmentUDT)            :: environment

    ! Debug
    logical                             :: makedebugplots = .false.
    real(R8), allocatable               :: xq(:), yq(:), vq(:), &
        xqmf(:,:), yqmf(:,:), vqmf(:,:)

    ! Set additional options
    !=======================
    ! Set grid options
    call SetGridOptions(gridoptions)

    ! Set optimization options
    call SetDesignOptions(designoptions)

    ! Set other numerical parameters
    call SetNumOptions(num);
    
    ! Magnetic field 
    call SetMagneticFieldOptions(mfoptions);

    ! Environment
    call SetEnvironmentOptions(environmentoptions)

    ! Initialize
    !===========
    ! Construct the initial grid
    call ConstructGrid(grid,gridoptions,options);

    ! Construct the initial magnetic field
    call ConstructMagneticField(mfoptions, magneticField)

    ! Construct the environment
    call ConstructEnvironment(environment, environmentoptions)

    ! Initialize the grid design problem (as an optimization problem)
    call ConstructGridDesignProblem(optimizationdriver, designoptions, &
        grid, magneticField, environment)

    ! Solve
    !======
    ! Simply call the solver ... 
    call optimizationdriver%Driver()

    ! ... and unpack the solution
    !grid = optimizationdriver%grid
    
    ! Do any post-processing if necessary

    ! Debug plots
    !============
    ! Make some plots
    if (makedebugplots) then
        !allocate(xq(grid%vert%ntot))
        !allocate(yq(grid%vert%ntot))
        !allocate(vq(grid%vert%ntot))
        !xq = grid%vert%x
        !yq = grid%vert%y
        
        allocate(xqmf(magneticField%nR-2, magneticField%nZ-2))
        allocate(yqmf(magneticField%nR-2, magneticField%nZ-2))
        allocate(vqmf(magneticField%nR-2, magneticField%nZ-2))
        allocate(xq((magneticField%nR-2)*(magneticField%nZ-2)))
        allocate(yq((magneticField%nR-2)*(magneticField%nZ-2)))
        allocate(vq((magneticField%nR-2)*(magneticField%nZ-2)))
        xqmf(:,:) = spread(magneticField%R(2:magneticField%nR-1), 2, magneticField%nZ-2)
        yqmf(:,:) = spread(magneticField%Z(2:magneticField%nZ-1), 1, magneticField%nR-2)
        xq = reshape(xqmf, (/((magneticField%nR-2)*(magneticField%nZ-2))/))
        yq = reshape(yqmf, (/((magneticField%nR-2)*(magneticField%nZ-2))/))
    
        call EvaluateBicubicSplineInterpolant(xq, yq, vq, magneticField%interp, '0', '0')
        vqmf = reshape(vq, (/magneticField%nR-2, magneticField%nZ-2/))
        call Plot2DStructuredField(magneticField%Psi, magneticField%R, magneticField%Z, magneticField%nR, magneticField%nZ, '-p')
        print *, size(magneticField%Psi,1), size(magneticField%Psi, 2)
        print *, size(vqmf,1), size(vqmf, 2)
        print *, size(yqmf, 1), size(yqmf, 2)
        print *, size(xq, 1), size(yq, 1), size(vq, 1)
        call Plot2DStructuredField(vqmf, magneticField%R(2:magneticField%nR-1), &
        magneticField%Z(2:magneticField%nZ-1), magneticField%nR-2, magneticField%nZ-2, '-p')
        
        ! call Plot2DUnstructuredField(vq, grid, 'v', '-p')
    end if



end subroutine