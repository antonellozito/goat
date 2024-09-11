!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains the implementation of the optimization engine for
! the grid deformation. It inherits from the base generic optimization
! structures in the optimization modules. 

! Notes
!======
! Note 1: Descriptions of what the deferred procedures should do are 
! provided in the abstract interface. 

module gdmod_optimizationengine
    
    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use optmod_optimizationengine
    use gdmod_designvariables
    use gdmod_costfunction
    use gdmod_constraints
    use gdmod_types

    ! The usual
    implicit none
    save
    public 

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    ! Optimization problem
    !=====================
    type, extends(OptimizationProblemUDT) :: OptimizationProblemGDUDT

        ! Description
        !============
        ! Overwrite the initial design routines with our own
        ! implemented and derived types. 

        ! The classical derived types
        class(DesignVariablesGDUDT), allocatable :: designvariables
        class(CostfunctionGDUDT), allocatable    :: costfunction
        type(ConstraintsGDUDT)                   :: constraints

        ! Additional fields that are needed
        type(GridUDT)                       :: grid 
        type(MagneticFieldUDT)              :: magneticField
        type(EnvironmentUDT)                :: environment
        type(DesignOptionsUDT)              :: designoptions 
        type(CGStructureUDT), allocatable   :: congroups(:) 
        type(DOFGStructureUDT), allocatable :: dofgroups(:)
        
    contains

        ! Overwrite 
        !==========
        ! Problem initialization
        procedure :: Initialize      => InitializeOptimizationProblemGD

        ! Dimension query
        procedure :: GetProblemDimensions   => GetProblemDimensionsGD 

        ! Design variables query
        procedure :: GetProblemDesignVariables => &
            GetProblemDesignVariablesGD

        ! Design update
        procedure :: UpdateDesign           => UpdateDesignGD

        ! Problem update
        procedure :: UpdateProblem          => UpdateProblemGD

        ! Problem parameter update
        procedure :: UpdateProblemParameters    => UpdateProblemParametersGD

        ! Cost function evaluation
        procedure :: EvaluateCostFunction   => EvaluateCostFunctionGD

        ! Equality constraints evaluation
        procedure :: EvaluateEqualityConstraints    &
                        => EvaluateEqualityConstraintsGD

        ! Inequality constraints evaluation
        procedure :: EvaluateInequalityConstraints &
                        => EvaluateInequalityConstraintsGD

        ! KKT relaxation
        procedure :: RelaxProblemKKTSystem => RelaxProblemKKTSystemGD

        ! Optimization iteration data writing
        procedure :: WriteIterationData => WriteIterationDataGD

        ! Additional routines
        !====================
        ! Initialization finalizer to account for cross-design/cfv/con
        ! initialization requirements
        procedure :: FinalizeInitialization  

        ! Inequality constraints group constructor
        procedure :: ConstructInequalityConstraintGroups

        ! Inequality constraint deactivator
        procedure :: DeactivateInequalityConstraints

        ! Linearizations
        !===============
        procedure :: EvaluateJacobian       => EvaluateJacobianGD
        procedure :: EvaluateJacobianGDGoatvariables
        procedure :: EvaluateJacobianGDVesselcoordinates

    end type 

    ! Optimization engine
    !====================
    type, extends(OptimizationEngineUDT) :: OptimizationEngineGDUDT 


    contains 

        ! Procedure for initialization of engine
        procedure :: SetupOptimizationDriver => SetupOptimizationDriverGD

    end type

    ! Abstract types
    !===============

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    ! Normal interfaces
    !==================
    interface SubtractDofs
        module procedure SubtractDofsReturnAttributed
        module procedure SubtractDofsDontReturnAttributed
    end interface

    contains

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                       OPTIMIZATION ENGINE                        !
    !------------------------------------------------------------------!
    ! Optimization engine initialization
    subroutine SetupOptimizationDriverGD(optimizationdriver) 

        ! Description
        !============
        ! This is only a dummy routine for setting up the optimization 
        ! engine and should be overwritten. Here, the design options etc
        ! are called to initialize the design. 

        ! Initialize
        !===========
        ! Declare modules
        !use gdmod_types
        !use gdmod_userinput

        ! The usual
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationEngineGDUDT)      :: optimizationdriver 

        ! Auxiliary
        type(OptimizationProblemGDUDT)      :: thisproblem
        type(OptimizationSolverKKTUDT)      :: thissolver

        ! Data

        ! Initialize
        !===========
        ! Allocate the problem type
        allocate(optimizationdriver%problem, source=thisproblem)

        ! Allocate the solver
        allocate(optimizationdriver%solver, source=thissolver)

        ! Propagate input filepath
        optimizationdriver%solver%inputfilepath = optimizationdriver%inputfilepath 
        optimizationdriver%solver%inputfileprefix = optimizationdriver%inputfileprefix
        optimizationdriver%problem%inputfilepath = optimizationdriver%inputfilepath

    end subroutine

    !------------------------------------------------------------------!
    !                       OPTIMIZATION PROBLEM                       !
    !------------------------------------------------------------------!
    ! Dimension query
    subroutine GetProblemDimensionsGD(problem, nphi, neq, nineq)

        ! Description
        !============
        ! Return the problem dimensions nphi, neq, nineq (number of 
        ! design variables, equality contraints, and inequality 
        ! constraints, resp. ). It is assumed that the problem is
        ! already properly initialized. 

        ! nphi is obtained from designvariables%nphi
        ! TO DO: change neq, nineq 

        ! Initialize
        !===========
        ! Declare modules

        ! The usual
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemGDUDT)     :: problem 
        integer(I8), intent(out)            :: nphi, neq, nineq

        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Output
        !=======
        ! Design variables
        nphi = problem%designvariables%nphi

        ! Constraints - call routine
        call problem%constraints%GetConstraintsDimensions(neq, nineq)

    end subroutine

    ! Design variables query
    subroutine GetProblemDesignVariablesGD(problem, phi)

        ! Description
        !============
        ! This routine returns the current design variable vector phi

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemGDUDT)         :: problem
        real(R8), allocatable, intent(out)      :: phi(:)

        ! Get!
        !=====
        phi = problem%designvariables%phi

    end subroutine

    ! Optimization problem initialization
    subroutine InitializeOptimizationProblemGD(problem) 

        ! Description
        !============
        ! This routine further initializes the design variables, cost
        ! function, and constraints. It is assumed that the grid, 
        ! magnetic field, and environment structures are properly 
        ! assigned. The allocatable components of the problem, namely 
        ! the cost function and design variables.

        ! Initialize
        !===========
        ! Declare modules
        use gdmod_types
        use gdmod_userinput

        ! The usual
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemGDUDT)      :: problem

        ! Loop variables

        ! Data

        ! Design options
        !===============
        ! Set the design options
        call problem%designoptions%Set()

        ! Design variables
        !=================
        ! Check if already allocated
        if (allocated(problem%designvariables)) then 
            ! Print and deallocate
            print *, 'InitializeOptimizationProblemGD: design variables ' // & 
                'already allocated, reinitializing...'

            deallocate(problem%designvariables)
        end if 

        ! Allocate the design variables, depending on the type.
        select case (trim(problem%designoptions%variables%type))

        case ('coordinates')

            ! Only coordinates
            allocate(DesignVariablesCoordinatesUDT::problem%designvariables)

        case ('desiredflux')

            ! Only flux values
            allocate(DesignVariablesFluxValuesUDT::problem%designvariables)

        case ('coordinates_desiredflux')

            ! Both coordinates and desired flux
            allocate(DesignVariablesCoordinatesFluxUDT::problem%designvariables)
            

        case default

            ! Throw error
            call gdErrorHandler('Unknown design variable type')

        end select

        ! Initialize the design variables
        call problem%designvariables%Initialize(problem%grid, &
            problem%magneticField, problem%environment)
        
        ! Cost function
        !==============
        if (allocated(problem%costfunction)) then 
            ! Print and deallocate
            print *, 'InitializeOptimizationProblemGD: cost function ' // & 
                'already allocated, reinitializing...'

            deallocate(problem%costfunction)
        end if
        
        ! Allocate the cost function, depending on the type
        select case (trim(problem%designoptions%costfunction%type))

        case ('lengthratio')

            ! Allocate
            allocate(CostFunctionLRUDT::problem%costfunction)

            ! Set type
            problem%costfunction%type = 'LR'

        case ('lengthratio2')

            ! Allocate
            allocate(CostFunctionLRUDT2::problem%costfunction)

            ! Set type
            problem%costfunction%type = 'LR2'

        case ('FAD')

            ! Allocate
            allocate(CostfunctionFADUDT::problem%costfunction)

            ! Set type
            problem%costfunction%type = 'FAD'

        case ('LR_FAD')

            ! Allocate
            allocate(CostfunctionLRFADUDT::problem%costfunction)

            ! Set type
            problem%costfunction%type = 'LR_FAD'

        case ('PRPB')

            ! Allocate
            allocate(CostfunctionPRPBUDT::problem%costfunction)

            ! Set type
            problem%costfunction%type = 'PRPB'
            
        case ('PRPB2')

            ! Allocate
            allocate(CostfunctionPRPB2UDT::problem%costfunction)

            ! Set type
            problem%costfunction%type = 'PRPB2'

        case ('general', 'LR_FAD_FA', 'LR_FAD_PRPB', 'LR_FAD_PRPB_FA', &
            'LR_FAD_PRPB_LRrad', 'LR_FAD_PRPB_LRrad_FA')

            ! Allocate
            allocate(CostfunctionGeneralUDT::problem%costfunction)

            ! Set type
            problem%costfunction%type = problem%designoptions%costfunction%type

        case default
            
            ! Throw error
            call gdErrorHandler('Unknown cost function type')

        end select

        ! Initialize the cost function
        call problem%costfunction%Initialize(problem%grid, &
            problem%magneticField, problem%environment, &
            problem%designoptions%costfunction)

        ! Constraints
        !============
        ! Given the (many) possible options for the constraints, the 
        ! constraints are set in its own initialization. 
        call problem%constraints%Initialize(problem%grid, &
            problem%magneticField, problem%environment, & 
            problem%designvariables, problem%designoptions%constraints)

        ! Set Lagrange multipliers
        if (allocated(problem%lambda)) then 
            deallocate(problem%lambda)
        end if 
        if (allocated(problem%mu)) then 
            deallocate(problem%mu)
        end if 

        allocate(problem%lambda(problem%constraints%eqcon%neqcon), &
            problem%mu(problem%constraints%ineqcon%nineqcon))
        problem%lambda = 0
        problem%mu = 0

        
        !=================
        ! Initialize design variables further for constraint/cfv 
        ! dependent fields
        call problem%FinalizeInitialization()


    end subroutine

    ! Finalize the problem initialization
    subroutine FinalizeInitialization(problem)

        ! Description
        !============
        ! This routine further initializes the optimization problem 
        ! after specific routines for initializing the design variables,
        ! constraints, and cost function have been called. As such, 
        ! interdependencies between these objects can be accounted for.
        ! For example, the amount of design variables for the 
        ! 'fluxvalue' (desired psi) type of design variables depends 
        ! on the number of flux function constraints. 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemGDUDT)     :: problem 

        ! Auxiliary

        ! Loop
        integer(I8)                         :: i, j

        ! Initialize
        !===========
        ! Associate
        associate(&
            nv                  => problem%grid%vert%ntot,      &
            x                   => problem%grid%vert%x,         &
            y                   => problem%grid%vert%y,         &
            npsi                => problem%constraints%eqcon%fluxfunction%nfluxsurfaces, &
            designvariables     => problem%designvariables,     &
            constraints         => problem%constraints,         &
            costfunction        => problem%costfunction         &
        )

        ! Update design
        !==============
        select type(designvariables)

        type is (DesignVariablesCoordinatesUDT)

            ! Nothing else to do here

        type is (DesignVariablesCoordinatesFluxUDT)

            ! Check
            if (.not. constraints%eqcon%dofluxfunction) then 
                ! Throw error - can't update if no constraints
                call gdErrorHandler('FinalizeInitialization: flux ' // &
                    'function constraints are not active and shoul be ' //&
                    ' for desiredpsi type design variables')
            end if 

            ! Need to initialize psi values further based on constraints
            if (allocated(designvariables%phi)) then 
                deallocate(designvariables%phi)
            end if

            ! Construct indices
            designvariables%xind = [(i, i = 1, nv)]
            designvariables%yind = [(i, i = nv+1, 2*nv)]
            designvariables%psiind = [(i, i = 2*nv+1, 2*nv+npsi)]

            ! Allocate & assign
            designvariables%nphi = npsi + 2*nv
            allocate(designvariables%phi(designvariables%nphi))
            designvariables%phi(designvariables%xind) = x
            designvariables%phi(designvariables%yind) = y
            designvariables%phi(designvariables%psiind) = constraints%eqcon%fluxfunction%fluxsurfaces%PsiD
            
            ! Set flux surface IDs
            designvariables%fsID = constraints%eqcon%fluxfunction%fluxsurfaces%fsID

        class default

            ! Throw error - unknown design variable type
            call gdErrorHandler('FinalizeInitialization: unknown design variable type')

        end select 

        ! Update cost function
        !=====================
        select type (costfunction)

        type is (CostfunctionPRPBUDT)

            ! Call initializer
            call costfunction%FinalizeInitialization(designvariables, &
                problem%grid, problem%magneticField, problem%environment, &
                problem%designoptions%costfunction)

        type is (CostfunctionPRPB2UDT)

            ! Call initializer
            call costfunction%FinalizeInitialization(designvariables, &
            problem%grid, problem%magneticField, problem%environment, &
            problem%designoptions%costfunction)

        type is (CostfunctionGeneralUDT)

            ! Check
            if (costfunction%doPRPB) then 
                call costfunction%cfv_prpb%FinalizeInitialization(designvariables, &
                    problem%grid, problem%magneticField, problem%environment, &
                    problem%designoptions%costfunction)
            end if

        class default 

            ! Do nothing

        end select

        ! Update constraints
        !===================
        ! Fixed flux values constraints
        if (constraints%eqcon%dofixedfluxvalues) then 

            ! Determine design variable indices for easier use later on
            select type (designvariables)

            type is (DesignVariablesCoordinatesFluxUDT) 

                do i = 1, size(constraints%eqcon%fixedfluxvalues%fsind, 1)
                    do j = 1, size(designvariables%psiind, 1)
                        if (constraints%eqcon%fixedfluxvalues%fsind(i) == designvariables%fsID(j)) then 
                            constraints%eqcon%fixedfluxvalues%psiind(i) = designvariables%psiind(j)
                        end if 
                    end do
                end do 

            class default
            
                call gdErrorHandler('FinalizeInitialization: fixed flux value constraint ' //&
                    'requires coordinates_desiredflux design variables' )
            end select
        end if 

        ! Construct constraint groups
        call problem%ConstructInequalityConstraintGroups()

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Problem update
    subroutine UpdateProblemGD(problem)

        ! Description
        !============
        ! Update the problem according to the design variables. Here,
        ! the grid, magnetic field, and environment should be updated,
        ! depending on the type of design variable. 

        ! Initialize
        !===========
        ! Declare modules

        ! The usual
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemGDUDT)     :: problem 

        ! Loop variables
        integer(I8)                         :: i

        ! Auxiliary variables 

        ! Data

        ! Initialize
        !===========
        ! Associate to use select type
        associate(&
            designvariables     => problem%designvariables, &
            constraints         => problem%constraints,     &
            grid                => problem%grid,            &
            magneticField       => problem%magneticField,   &
            environment         => problem%environment      &
            )

        ! Update design
        !==============
        ! Simply call the update routine from the design variables
        call designvariables%UpdateDesign(grid, magneticField, &
            environment)

        ! Update other fields
        !====================
        ! Actually this should be migrated to a different routine? E.G.
        ! UpdateConstraints or something? 
        select type (designvariables)

        type is (DesignVariablesCoordinatesFluxUDT)

            ! Update the flux function constraints - coordinates should
            ! be done already in UpdateDesign
            do i = 1, constraints%eqcon%fluxfunction%nfluxsurfaces
                constraints%eqcon%fluxfunction%fluxsurfaces(i)%PsiD = designvariables%phi(designvariables%psiind(i))
            end do
        
        class default 

            ! Do nothing

        end select

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Design update
    subroutine UpdateDesignGD(problem, dx)

        ! Description
        !============
        ! Update the design variables according to the update dx. 

        ! Initialize
        !===========
        ! Declare modules
        use gdmod_plots

        ! The usual
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemGDUDT)     :: problem 
        real(R8), intent(in)                :: dx(:)
    
        ! Loop variables

        ! Auxiliary variables 
        logical                             :: dodebugplots

        character(:), allocatable           :: vertpath, cellpath

        ! Data

        ! Do NOT set this to true unless you
        ! know what you're doing! prints out a figure, which may be in 
        ! each optimization iteration!
        
        data dodebugplots /.false./ 

        ! Update design
        !==============
        ! Simply update designvariables%phi
        problem%designvariables%phi = problem%designvariables%phi + dx 

        ! Debug plots?
        !=============
        if (dodebugplots) then
            ! Call grid visualization
            call PlotGridCells(problem%grid, '-p')
        end if
        if (problem%designoptions%writedata == 1) then 
            ! Call grid vertex writing routine
            allocate(character(len('vertices_iterate')) :: vertpath)
            allocate(character(len('cells_iterate')) :: cellpath)
            vertpath = 'vertices_iterate'
            cellpath = 'cells_iterate'
            call WriteGridVertices(problem%grid, vertpath) 
            call WriteGridCells(problem%grid, cellpath)
        end if

    end subroutine

    ! Problem parameters update
    subroutine UpdateProblemParametersGD(problem, values, updatemeth)

        ! Description
        !============
        ! Update the problem parameters (magnetic field &
        ! environment) and related cost function and constraint 
        ! parameters (design variables shouldn't change, as they are
        ! part of the problem definition and not a parameter of the 
        ! problem). Some examples of parameters:
        !
        ! - the shape of the vessel
        ! - the magnetic field (though topology shouldn't change!)
        ! 
        ! The following update methodes are currently implemented:
        !
        ! - vesselcoordinates: update vessel coordinates. Here, the 
        !   values are first all x-coordinates, then all y-coordinates,
        !   i.e. (x(i), y(i+np)) where np is the number of points forms
        !   the coordinates of the i-th point. 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemGDUDT)     :: problem 
        real(R8), intent(in)                :: values(:)
        character(*), intent(in)            :: updatemeth 

        ! Auxiliary
        logical                             :: upconbnd, upconlf, &
            upconiv

        integer(I8)                         :: np

        ! Initialize
        !===========
        ! Checks
        if (size(values, 1) == 0) then 
            ! Nothing to do here
            return 
        end if 

        ! Associate
        associate(&
            grid            => problem%grid,            &
            magneticField   => problem%magneticField,   &
            environment     => problem%environment)

        ! Set initial switches to update cost function and constraints
        upconbnd    = .false. 
        upconlf     = .false.
        upconiv     = .false.

        ! Update parameters
        !==================
        select case(trim(updatemeth))

        case ('vesselcoordinates')

            ! Update vessel coordinates
            np = size(values, 1)/2
            if (2*np /= size(values, 1)) then 
                ! This shouldn't happen
                call gdErrorHandler('UpdateProblemParameters: uneven number of coordinates detected, check input')
            end if 
            call problem%environment%vessel%UpdateVesselCoordinates(values(1:np), values(np+1:2*np))

            ! Select which constraints/cost function contributions to 
            ! update
            upconbnd    = .true. 
            upconlf     = .true.
            upconiv     = .true.

        case default

            ! Case not implemented, throw error
            call gdErrorHandler('UpdateProblemParameters: unknown updatemeth: "' // &
                trim(updatemeth) // '", exiting...')

        end select

        ! Update cost function
        !=====================

        ! Update constraints
        !===================
        ! Boundary function constraints
        if (upconbnd .and. problem%constraints%eqcon%doboundaryfunction) then 
            call problem%constraints%eqcon%boundaryfunction%Update(grid, &
                magneticField, environment)
        end if 

        ! Linefolding constraints
        if (upconlf .and. problem%constraints%ineqcon%dolinefolding) then 
            call problem%constraints%ineqcon%linefolding%Update(grid, &
                magneticField, environment)
        end if

        ! In vessel constraints
        if (upconiv .and. problem%constraints%ineqcon%doinvessel) then 
            call problem%constraints%ineqcon%invessel%Update(grid, &
                magneticField, environment)
        end if 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Cost function evaluation
    subroutine EvaluateCostFunctionGD(problem, J, gradJ, hessJ, &
        dogradient, dohessian)

        ! Description
        !============
        ! Evaluate the cost function, based on the current state of the
        ! optimization problem. Besides the scalar value itself, also
        ! the gradient and hessian are evaluated, if needed. This can 
        ! be set using the dogradient and dohessian logicals. 

        ! Initialize
        !===========
        ! Declare modules

        ! The usual
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemGDUDT)     :: problem 
        real(R8)                            :: J 
        real(R8), allocatable               :: gradJ(:)
        type(MySparseUDT)                   :: hessJ
        logical                             :: dogradient, dohessian                            

        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Compute the cost function
        !==========================
        ! Simply call the cost function computation routine
        call problem%costfunction%Evaluate(J, gradJ, hessJ, problem%grid, &
            problem%magneticField, problem%environment, dogradient, &
            dohessian, problem%designvariables)

    end subroutine

    ! Equality constraints evaluation
    subroutine EvaluateEqualityConstraintsGD(problem, G, gradG, hessG, &
        dogradient, dohessian, lambda)

        ! Description
        !============
        ! Evaluate the constraints, based on the current state of the
        ! problem. Note that the 'hessian' that is returned, is in fact
        ! the hessian-vector multiplication between the hessian and the
        ! vector lambda. Lambda should have been initialized with the 
        ! proper size. 

        ! Initialize
        !===========
        ! Declare modules

        ! The usual
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemGDUDT)     :: problem 
        real(R8), allocatable               :: G(:), lambda(:)
        type(MySparseUDT)                   :: gradG, hessG
        logical                             :: dogradient, dohessian                            

        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Compute the constraints
        !========================
        ! Simply call the constraint computation routine
        call problem%constraints%eqcon%Evaluate(G, gradG, hessG, problem%grid, &
            problem%magneticField, problem%environment, dogradient, &
            dohessian, problem%designvariables, lambda)

    end subroutine

    ! Inequality constraints evaluation
    subroutine EvaluateInequalityConstraintsGD(problem, H, gradH, &
        hessH, dogradient, dohessian, mu)

        ! Description
        !============
        ! Evaluate the constraints, based on the current state of the
        ! problem. Note that the 'hessian' that is returned, is in fact
        ! the hessian-vector multiplication between the hessian and the
        ! vector mu. mu should have been initialized with the 
        ! proper size. 

        ! Initialize
        !===========
        ! Declare modules

        ! The usual
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemGDUDT)     :: problem 
        real(R8), allocatable               :: H(:), mu(:)
        type(MySparseUDT)                   :: gradH, hessH
        logical                             :: dogradient, dohessian                            

        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Compute the constraints
        !========================
        ! Simply call the constraint computation routine
        call problem%constraints%ineqcon%Evaluate(H, gradH, hessH, &
            problem%grid, problem%magneticField, problem%environment, &
            dogradient, dohessian, problem%designvariables, mu)

        ! Deactivate constraints
        !=======================
        ! Only if necessary - determined in the routine itself
        ! Note: should we add switch here? And adjust gradient/Hessian
        ! as well? At least gradient should be tackled by activeness
        ! of constraint...
        call problem%DeactivateInequalityConstraints(H)

    end subroutine

    ! KKT relaxation 
    subroutine RelaxProblemKKTSystemGD(problem, KKT)

        ! Description
        !============
        ! Relax KKT system based on problem parameters. Note: this is not
        ! necessary/not supported yet for the general GD problem. If 
        ! used, an error will be thrown. 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemGDUDT)         :: problem 
        type(MySparseUDT)                       :: KKT 

        ! Call error handler
        call gdErrorHandler('RelaxProblemKKTSystemGD: method not yet implemented')

    end subroutine

    ! Optimization iteration data writing
    subroutine WriteIterationDataGD(problem, itopt)

        ! Description
        !============
        ! Write out any data of the optimization problem per iteration
        ! to a specific file. Uses the plotter module in the backend.
        ! Here, we don't use the itopt variable, which gives the 
        ! iteration number, so we overwrite each file.
        ! The following data is written out:
        ! - temp_gridcellsiterate.dat: grid cell coordinates with vertices for plotting (overwritten each time)
        ! - history.dat (created at iteration one, appended each iteration)

        ! Modules
        !========
        use mod_plotter, only: plotdir 
        use mod_specialchars, only: filesepchar

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemGDUDT)         :: problem 
        integer(I8)                             :: itopt 

        ! Auxiliary 
        integer, parameter                      :: fid = 70
        integer                                 :: tiostat
        logical                                 :: isfile 
        character(:), allocatable               :: filepath

        ! Initialize
        !===========
        ! Associate some fields for easier reading/writing
        associate(grid      => problem%grid)

        ! Set filepath
        filepath = plotdir // filesepchar // 'goat_optimization_history.dat'

        ! Write data
        !===========
        ! Cell vertex data
        call WriteGridCells(grid, 'temp_gridcellsiterate')

        ! Iteration data
        if (itopt == 1) then 
            ! Check if the file exists
            inquire(file=filepath, exist=isfile)
            if (isfile) then 
                ! Replace old file
                open(unit=fid, status='old', iostat=tiostat, file=filepath)
                rewind(fid)
            else
                ! Create new file
                open(unit=fid, status='new', iostat=tiostat, file=filepath)
            end if

            ! Write header
            call problem%monitor%WriteFileHeader(fid)
        else
            ! File should already by opened, append
            open(unit=fid, status='old', iostat=tiostat, file=filepath, access='sequential', position='append')
        end if

        ! Write
        call problem%monitor%WriteFileIterate(fid)

        ! Housekeeping
        !=============
        close(unit=fid)
        end associate

    end subroutine

    ! Constraint group construction
    subroutine ConstructInequalityConstraintGroups(problem)

        ! Description
        !============
        ! This routine sets up the groups and constraint groups to be used in the
        ! 'DeactiveInequalityConstraints' routine to (hopefully) render the
        ! optimization problem feasible in case of too many imposed constraints. It
        ! is imperative to realize that 'too many imposed constraints' is fully
        ! determined by the implementation of this routine! The following
        ! structures are added to the designParams.ineqcon structure:
        !
        ! - groups:     structure with 'cons' and 'dofs' fields. 'cons' holds the
        !               inequality constraint indices that affect the group, 'dofs'
        !               indicate how many inequality constraints can be active
        !               additional to the already imposed equality constraints. 
        ! - congroups:  structure with the 'groups' field, which holds the groups
        !               to which the inequality constraint can be assigned. 
        !
        ! Under 'groups' we understand either single vertices, flux surfaces or
        ! radial lines. The reasoning to determine the groups and their dofs is the
        ! following:
        !
        ! - flux surfaces are sets of vertices with the same (non-zero) fieldline 
        !   ID. A flux surface can have maximally one dof, namely their flux
        !   function value (if this is a design variable, otherwise it has no
        !   additional degree of freedom). 
        ! - radial lines are sets of vertices where each vertex has at least one,
        !   at most two faces in common with another vertex. All of these faces
        !   should have orthogonality constraints. Each radial line can therefore
        !   be seen as equivalent to a flux surface. They also have maximally one
        !   degree of freedom. 
        ! - Single vertices are simply vertices. They have two degrees of freedom,
        !   minus the amount of flux surfaces and radial lines they belong to.
        !   Normally, a vertex should not belong to more than one flux surface and
        !   one radial line (except for x-points). 
        !
        ! The attribution of the inequality constraints to each of these groups
        ! depends on the nature of the inequality constraints:
        !
        ! - poloidal linefolding: can only be attributed to single vertices or
        !   radial lines (changing psi value does not necessarily change the
        !   poloidal folding)
        ! - radial linefolding: same as poloidal linefolding, but now radial lines
        !   should't be used
        ! - vessel linefolding: all groups 
        ! - convexity: all groups
        ! - region: all groups 

        ! Algorithm
        !----------
        ! 1) Build groups:
        !   1.1) Flux surfaces should be present in the constraints already
        !   1.2) Radial lines must be reconstructed using GetOrthogonalLines
        !   1.3) All vertices are single vertices)
        ! 2) Determine dofs:
        !   2.1) Flux surfaces and radial lines start at one. Imposed psi value set
        !   it to zero. 
        !   2.2) Check which vertices are constrained to be on a poloidal/radial
        !   surface. Start at two, reduce by one for each constraint. 
        !   2.3) Check the remaining equality constraints. For each constraint,
        !   first try to attribute it to one of the single vertices that still have
        !   degrees of freedom. If this is not possible, attribute to a radial or
        !   poloidal surface:
        !       * Boundary constraints: both possible
        !       * Edge length constraints: poloidal if both same vertID (and
        !       nonzero), otherwise radial. If none possible -> probably overly
        !       constrained, issue warning
        !       * x-point: to be ignored, zero dofs (shouldn't move)
        !       * flux surface value constraints: obviously only flux surfaces
        !       * other constraints not supported yet
        ! 3) Determine congroups:
        !   3.1) Groups with already zero dofs can be eliminated a priori
        !   3.2) Any inequality constraint can be attributed to a single vertex.
        !   For others, see rules stated above. 


        ! Notes
        !------
        ! Note 1: currently, only linefolding, region, and convexity constraints 
        ! are considered. 

        ! Note 2: it is assumed (and this should anyhow be the case) that a vertex
        ! belongs maximally to one flux surface and one radial line (x-points have
        ! an exception for the latter)

        ! Note 3: flux surfaces are identified using the data of the flux surface
        ! constraints. Here, only vertices that belong to the flux surface
        ! substructure are considered - any vertices in the vpairs structure are
        ! not accounted for. If there are any present, a warning is thrown. 

        ! Note 4: it is assumed that for the equality constraints, the necessary
        ! constraint qualifications are satisfied. 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemGDUDT) :: problem 

        ! Auxiliary
        integer(I8)                         :: nfs, nrl, fsdof, nxpind, &
            gID, dummy_I8
        integer(I8), allocatable            :: vertrl(:, :), vertfs(:), &
            nvertrl(:), dofs(:), xpind(:), order(:), tv(:), &
            tvf(:, :), groupID(:), tID(:), groupindfs(:), groupindrl(:), &
            groupindvert(:)

        logical                             :: fslegal, rllegal, &
            throwerror, attributed, vertlegal, isattributed
        logical, allocatable                :: hasrl(:), hasfs(:), &
            hasmultrl(:), isineqconattributed(:)
        
        type(FFCStructureUDT), allocatable  :: fs(:) 
        type(PolygonSetUDT)                 :: rl
        

        ! Loop
        integer(I8)                         :: i, j, k, ng, cc

        ! Initialize
        !===========
        ! Associate
        associate(&
            nv              => problem%grid%vert%ntot,      &
            grid            => problem%grid,                &
            eqcon           => problem%constraints%eqcon,   &
            ineqcon         => problem%constraints%ineqcon, &
            designvariables => problem%designvariables      &
        )

        ! Check for flux surfaces
        if (eqcon%dofluxfunction) then 
            ! Unpack
            fs = eqcon%fluxfunction%fluxsurfaces
            nfs = size(fs, 1)
        else
            nfs = 0
        end if 

        ! Check for radial lines
        if (eqcon%doorthogonality) then 
            ! Unpack
            rl = eqcon%orthogonality%radiallines
            nrl = rl%np 
        else
            nrl = 0
        end if 

        ! Allocate
        if (allocated(problem%congroups)) then 
            deallocate(problem%congroups)
        end if 
        if (allocated(problem%dofgroups)) then 
            deallocate(problem%dofgroups)
        end if 
        allocate(hasrl(nv), hasfs(nv), nvertrl(nv), vertfs(nv))
        allocate(vertrl(nv, maxval(grid%vert%neigP(:, 2)))) ! vertex can belong to multiple radial lines if x-point... initialized too big here
        allocate(dofs(nv + nfs + nrl), problem%dofgroups(nv + nfs + nrl))
        allocate(problem%congroups(ineqcon%nineqcon), isineqconattributed(ineqcon%nineqcon))

        ! Associate further
        associate(&
            dofgroups       => problem%dofgroups,   &
            congroups       => problem%congroups    &
            )

        ! Initialize
        hasrl               = .false.
        hasfs               = .false.
        throwerror          = .true.
        isineqconattributed = .false.
        isattributed        = .false.

        vertrl  = 0
        nvertrl = 0
        vertfs  = 0
        dofs    = 0

        ! Check how many dofs a flux surface has
        fsdof = 0
        select case (designvariables%type)

        case ('coordinates_desiredflux')

            fsdof = 1

        case default

            fsdof = 0

        end select

        ! Build dof groups
        !=================
        ! Counters
        ng = 0 ! group counter

        ! Single vertices
        groupindvert = [(k, k = 1, nv)] + ng 
        do i = 1, nv
            ! Update counter
            ng = ng + 1
            
            ! Add dofs and vert 
            dofs(ng) = 2 ! default
            allocate(dofgroups(ng)%vert(1))
            dofgroups(ng)%vert = i 

        end do

        ! Flux surfaces
        groupindfs = [(k, k = 1, nfs)] + ng
        do i = 1, nfs 
            ! Update counter
            ng = ng + 1

            ! Add dofs and vert
            dofs(ng) = fsdof 
            dofgroups(ng)%vert = fs(i)%ID

            ! Check
            if (any(hasfs(fs(i)%ID))) then 
                ! Multiple flux surfaces with same vertex, not allowed
                call gdErrorHandler('ConstructInequalityConstraintGroups: ' // &
                    'vertices identified that belong to multiple flux surfaces, ' // &
                    'not supported') 
            end if
            hasfs(fs(i)%ID)     = .true. 
            vertfs(fs(i)%ID)    = i
        end do

        ! Radial lines
        groupindrl = [(k, k = 1, nrl)] + ng 
        do i = 1, nrl 
            ! Update counter
            ng = ng + 1

            ! Add dofs and vert
            dofs(ng) = 1
            tID = rl%polygons(i)%labels(:, 1)
            dofgroups(ng)%vert = tID

            ! Here, multiple radial lines are allowed
            hasrl(tID) = .true. 
            do j = 1, size(tID, 1)
                nvertrl(tID(j)) = nvertrl(tID(j)) + 1
                vertrl(tID(j), nvertrl(tID(j))) = i 
            end do
        end do

        ! Subtract dofs due to eqcons
        !============================
        ! Flux surfaces & radial lines
        do i = 1, nv 
            if (hasrl(i)) then 
                dofs(i) = dofs(i) - 1
            end if 
            if (hasfs(i)) then 
                dofs(i) = dofs(i) - 1
            end if 
        end do 

        ! Check if only X-points have multiple radial lines
        hasmultrl = nvertrl > 1 
        call DetermineXPoints(xpind, nxpind, order, grid)
        hasmultrl(xpind) = .false.
        if (any(hasmultrl)) then 
            ! Throw error, is not allowed
            call gdErrorHandler('ConstructInequalityConstraintGroups: ' // &
                'non X-points detected with multiple radial lines')
        end if 

        ! Boundary functions
        if (eqcon%doboundaryfunction) then 

            ! Set attribution rules
            fslegal     = .true.
            rllegal     = .true.
            vertlegal   = .true.

            ! Subtract dofs
            call SubtractDofs(eqcon%boundaryfunction%vert, dofs, &
                fslegal, rllegal, vertlegal, hasfs, hasrl, groupindrl, &
                groupindfs, vertrl, vertfs, nvertrl, throwerror)
            
        end if

        ! X-points
        if (eqcon%doxpoints) then 
            ! X-point indices already extracted upstream
            
            ! Set attribution rules
            rllegal = .true.
            fslegal = .true.
            vertlegal   = .true.
            call SubtractDofs(xpind, dofs, fslegal, rllegal, vertlegal, hasfs, hasrl, &
                groupindrl, groupindfs, vertrl, vertfs, nvertrl, throwerror)
        end if 

        ! Fluxfunction
        if (eqcon%dofluxfunction) then 
            ! Note: only special and fixed points are considered, not
            ! flux surfaces.
            
            ! Set attribution rules
            rllegal = .true.
            fslegal = .true.
            vertlegal   = .true.
            
            ! Fixed points
            do i = 1, eqcon%fluxfunction%nfixedpoints
                tv = eqcon%fluxfunction%fixedpoints(i)%ID
                call SubtractDofs(tv, dofs, fslegal, rllegal, vertlegal, hasfs, hasrl, &
                    groupindrl, groupindfs, vertrl, vertfs, nvertrl, throwerror)
            end do 
            
            ! Special points
            associate(sp => eqcon%fluxfunction%specialpoints)
            do i = 1, size(sp, 1)
                tv = sp(i)%ID
                call SubtractDofs(tv, dofs, fslegal, rllegal, vertlegal, hasfs, hasrl, &
                    groupindrl, groupindfs, vertrl, vertfs, nvertrl, throwerror)
            end do
            end associate

            ! Tangency points
            associate(tp => eqcon%fluxfunction%tangencypoints)
            do i = 1, size(tp, 1)
                tv = tp(i)%ID
                call SubtractDofs(tv, dofs, fslegal, rllegal, vertlegal, hasfs, hasrl, &
                    groupindrl, groupindfs, vertrl, vertfs, nvertrl, throwerror)
            end do
            end associate
        end if

        ! Orthogonality 
        ! (no additional steps)

        ! Edge lengths
        if (eqcon%doedgelengths) then 
            ! Extract vertices
            tvf = eqcon%edgelengths%edgevert

            ! Set attribution rules
            rllegal = .true.
            fslegal = .true.
            vertlegal   = .true.
            
            ! Attribute constraint
            do i = 1, size(tvf, 1)
                call AttributeConstraint(tvf(i, :), dofs, fslegal, rllegal, vertlegal, hasfs, hasrl, &
                    groupindrl, groupindfs, vertrl, vertfs, nvertrl, dummy_I8, isattributed, &
                    .true.)
            end do 
            
        end if

        ! Fixed flux values 
        if (eqcon%dofixedfluxvalues) then 

            ! Loop
            associate(tfsID => eqcon%fixedfluxvalues%fsind)
            do i = 1, size(tfsID, 1)
                ! Try to attribute first to flux surfaces only
                tv = fs(tfsID(i))%ID
                call AttributeConstraint(tv, dofs, .true., .false., .false., hasfs, hasrl, &
                    groupindrl, groupindfs, vertrl, vertfs, nvertrl, dummy_I8, &
                    isattributed, .false.)

                ! If that didn't work, try to attribute to vertices
                if (.not. isattributed) then 
                    call AttributeConstraint(tv, dofs, .true., .false., &
                        .true., hasfs, hasrl, &
                        groupindrl, groupindfs, vertrl, vertfs, nvertrl, dummy_I8, &
                        isattributed, .true.)
                end if 
                
            end do
            end associate 
        end if 

        ! Attribute inequality constraints
        !=================================
        ! Do this without updating the dofs! Note that all considered constraints
        ! already do some elimination, i.e. they are attributed already to a single
        ! node, so we only need to see if this node itself still has dofs, or if it
        ! has radial lines/flux surfaces with non-zero dofs. 

        ! IMPORTANT: the order of attributing constraints has to be exactly the
        ! same as the order in which the constraints are evaluated! 

        ! Constraint counter
        cc = 0

        ! Invessel
        if (ineqcon%doinvessel) then 
            
            ! Associate
            associate( &
                vert          => ineqcon%invessel%vert,         &
                nvert         => ineqcon%invessel%nvert         &
                )

            ! Set rules
            fslegal = .true. 
            rllegal = .true.
            vertlegal = .true. 

            ! Loop over all vertices
            do i = 1, nvert
                ! Update counter
                cc = cc + 1

                ! Attribute
                call DetermineConstraintDofgroups([vert(i)], dofs, fslegal, rllegal, &
                    vertlegal, hasfs, hasrl, groupindrl, groupindfs, vertrl, vertfs, nvertrl, &
                    attributed, groupID)
                
                ! Add
                congroups(cc)%dofgroups = groupID
                isineqconattributed(cc) = attributed

            end do 

            if (any(.not. isineqconattributed(cc-nvert+1:cc))) then 
                ! Throw warning, some constraints will not be set
                print *, 'ConstructInequalityConstraintGroups: some ' // & 
                    'invessel inequalities will never be active as ' // &
                    'they could not be attributed'
            end if
           
            ! Housekeeping
            end associate
        end if 


        ! Linefolding
        if (ineqcon%dolinefolding) then 
            
            ! Associate
            !----------
            associate( &
                vpairspol     => ineqcon%linefolding%vpairspol,     &
                vpairsrad     => ineqcon%linefolding%vpairsrad,     &
                vpairsves     => ineqcon%linefolding%vpairsves,     &
                nvpairspol    => ineqcon%linefolding%nvpairspol,    &
                nvpairsrad    => ineqcon%linefolding%nvpairsrad,    &
                nvpairsves    => ineqcon%linefolding%nvpairsves     &
                )

            ! Poloidal contributions
            !-----------------------
            ! Set rules
            fslegal = .false. 
            rllegal = .true.
            vertlegal = .true. 

            ! Loop over all vertex pairs
            do i = 1, nvpairspol
                ! Update counter
                cc = cc + 1

                ! Attribute
                call DetermineConstraintDofgroups(vpairspol(i, :), dofs, fslegal, rllegal, &
                    vertlegal, hasfs, hasrl, groupindrl, groupindfs, vertrl, vertfs, nvertrl, &
                    attributed, groupID)
                
                ! Add
                congroups(cc)%dofgroups = groupID
                isineqconattributed(cc) = attributed

            end do 

            if (any(.not. isineqconattributed(cc-nvpairspol+1:cc))) then 
                ! Throw warning, some constraints will not be set
                print *, 'ConstructInequalityConstraintGroups: some ' // & 
                    'poloidal linefolding inequalities will never be active as ' // &
                    'they could not be attributed'
            end if

            ! Radial contributions
            !---------------------
            ! Set rules
            fslegal = .true.
            rllegal = .false. 
            vertlegal = .true. 

            ! Loop over all vertex pairs
            do i = 1, nvpairsrad
                ! Update counter
                cc = cc + 1

                ! Attribute
                call DetermineConstraintDofgroups(vpairsrad(i, :), dofs, fslegal, rllegal, &
                    vertlegal, hasfs, hasrl, groupindrl, groupindfs, vertrl, vertfs, nvertrl, &
                    attributed, groupID)

                ! Add
                congroups(cc)%dofgroups = groupID
                isineqconattributed(cc) = attributed
            end do

            ! Check
            if (any(.not. isineqconattributed(cc-nvpairsrad+1:cc))) then 
                ! Throw warning, some constraints will not be set
                print *, 'ConstructInequalityConstraintGroups: some ' // & 
                    'radial linefolding inequalities will never be active as ' // &
                    'they could not be attributed'
            end if

            ! Vessel contributions
            !---------------------
            ! Set rules
            fslegal = .true.
            rllegal = .true.
            vertlegal = .true. 

            ! Loop over all vertex pairs
            do i = 1, nvpairsves
                ! Update counter
                cc = cc + 1

                ! Attribute
                call DetermineConstraintDofgroups(vpairsves(i, :), dofs, fslegal, rllegal, &
                    vertlegal, hasfs, hasrl, groupindrl, groupindfs, vertrl, vertfs, nvertrl, &
                    attributed, groupID)

                ! Add
                congroups(cc)%dofgroups = groupID
                isineqconattributed(cc) = attributed
            end do
        
                
            ! Check
            if (any(.not. isineqconattributed(cc-nvpairsves+1:cc))) then 
                ! Throw warning, some constraints will not be set
                print *, 'ConstructInequalityConstraintGroups: some ' // & 
                    'vessel linefolding inequalities will never be active as ' // &
                    'they could not be attributed'
            end if

            ! Housekeeping
            end associate
        end if 

        ! Convexity
        !if ineqcon.convexity
        !    ! Attribute constraints
        !    for j = 1:numel(ineqcon.parameters.cvert)
        !        cc = cc + 1;
        !        [attributed, groupID] = AttributeInequalityConstraints(ineqcon.parameters.cvert(j), dofs, fslegal, rllegal, hasfs, hasrl, ...
        !            groupindrl, groupindfs, vertrl, vertfs, nvertrl);
        !        congroups(cc).groups = groupID;
        !    end
        !    
        !    ! Check
        !    if any(~attributed)
        !        ! Throw warning, some constraints will not be set
        !        disp('ConstructInequalityConstraintGroups: some convexity inequalities will never be active as they could not be attributed')
        !    end
        !end

        ! Attribute constraints to groups
        !--------------------------------
        do i = 1, size(dofgroups)
            allocate(dofgroups(i)%cons(0))
        end do
        do i = 1, cc
            do j = 1, size(congroups(i)%dofgroups)
                gID = congroups(i)%dofgroups(j)
                dofgroups(gID)%cons = [dofgroups(gID)%cons, i]
            end do 
        end do 

        ! Add dofs
        do  i = 1, ng
            dofgroups(i)%dofs = dofs(i)
        end do 


        ! Housekeeping
        !=============
        end associate
        end associate
        

    end subroutine

    ! Inequality constraint deactivation
    subroutine DeactivateInequalityConstraints(problem, H)

        ! Description
        !============
        ! This function checks whether constraints should be deactivated in order
        ! to obtain a feasible problem. Note that this is not done based on the
        ! rank of the full Hessian of the problem (which could also be a
        ! possibility) but simply on insight in the problem. 
        ! This, of course, requires input from the user. Specifically, the
        ! following data should be given in ineqcon:
        !
        ! - groups:         ng-by-1 structure with fields 'dofs'. This structure
        !                   indicates for each (group of) vertices/... how many
        !                   constraints can additionally be imposed AFTER imposing
        !                   equality constraints (so the actual degrees of freedom
        !                   for that group of vertices). Additionally, the field
        !                   'cons' should be present which indicates all inequality
        !                   constraints that can be attributed tot that group. 
        ! - congroups:      nc-by-1 structure with fields 'groups'. Eeach ith
        !                   structure should correspond with the ith constraint,
        !                   and each 'groups' entry there should contain the
        !                   indices of which groups the constraint can be
        !                   contributed to. 
        !
        ! It is important to realize that everything depends on this user input: if
        ! this is not properly thought of, one may still encounter infeasible
        ! problems! In other words, this is a rather heuristic way of tackling
        ! possible infeasibility. 
        !
        ! Algorithm
        !==========
        ! The constraint attribution is done as follows:
        !   1) We check for each group how many constraints are active for that
        !   group. 
        !   2) Groups that only have one constraint are considered first. If they
        !   have at least one dof, the constraint is imposed. 
        !   3) Next, we consider the group(s) with the largest number of
        !   dofs and determine for that group the largest
        !   constraint value. This constraint is attributed. If each group either
        !   has no more dofs, or has no more active constraints, the algorithm is
        !   terminated. Any active constraints that are left are set to be
        !   inactive, as they are expected to lead to an infeasible problem. 
        
        ! Notes
        !======
        ! Note 1: observe that the 'groups' do not have to coincide with a group of
        ! vertices whatsoever and allows for a general implementation of which
        ! constraints to impose when and where. 
        
        ! Note 2: the algorithm mentioned above has no general convergence proof.
        ! However, since we always take the largest value of the constraint, one
        ! should eventually be able to identify the active set, provided that at
        ! the optimum the constraints satisfy the constraint qualifications (I
        ! think, not sure) 
        
        ! Note 3: if certain groups should be constrained first, they should appear
        ! first in the group structure. Because we start looping at the beginning
        ! of this structure, these will receive their constraints first. 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemGDUDT)     :: problem
        real(R8), intent(inout)             :: H(:)

        ! Auxiliary
        integer(I8)                         :: ng, nh, Hmaxind, &
            thisHind
        integer(I8), allocatable            :: dofs(:), nacg(:), &
            Hind(:), tg(:)

        real(R8)                            :: Hmaxval 
        real(R8), allocatable               :: Hval(:)
        logical, allocatable                :: A(:), nac(:)

        ! Loop
        integer(I8)                         :: i
        
        ! Initialize
        !===========
        ! Check
        if (size(H, 1) == 0) then 
            ! No inequality constraints present, return
            return
        end if 
        
        ! Determine active constraints
        A = H >= 0
        if (count(A) == 0) then 
            ! No active inequality constraints, return
            return
        end if 

        ! Associate
        associate(&
            congroups               => problem%congroups,   &
            dofgroups               => problem%dofgroups    &
            )
        
        ! Unpack
        ng          = size(dofgroups, 1)
        nh          = size(H, 1)
        dofs        = dofgroups%dofs ! degree of freedom counter
        nac         = A ! non-attributed constraints (nac)
        allocate(nacg(ng))
        nacg        = 0 ! number of active constraints per group
        
        ! Compute active constraints for each group
        !==========================================
        do i = 1, ng
            nacg(i) = count(A(dofgroups(i)%cons))
        end do
        
        ! Attribute constraints
        !======================
        ! dofgroups with only one constraint
        i = 1
        allocate(Hind(1))
        do while (i <= ng)
            if (nacg(i) == 1) then 
                
                ! Check
                if (count(nac(dofgroups(i)%cons)) /= 1) then 
                    ! This shouldn't happen and points to an implementation bug
                    call gdErrorHandler('DeactivateInequalityConstraints: implementation bug detected')
                end if 

                ! Get nacs
                Hind = pack(dofgroups(i)%cons, nac(dofgroups(i)%cons))
                
                ! Check dofs
                if (dofs(i) > 0) then 
                    ! Attribute
                    nac(Hind) = .false.
                    
                    ! Update dofs
                    dofs(i) = dofs(i) - 1
                    
                    ! Update constraint counter
                    tg = congroups(Hind(1))%dofgroups
                    nacg(tg) = nacg(tg) - 1 
                    
                    ! Reset counter to start at the beginning again to eliminate
                    ! all dofgroups with one constraint eventually
                    i = 1
                else
                    
                    ! Update counter
                    i = i + 1
                    
                end if
            else
                
                ! Update counter
                i = i + 1
                
            end if 
        end do
        deallocate(Hind)
        
        ! Sanity check
        if (any( (nacg == 1) .and. (dofs > 0))) then 
            ! This shouldn't be happening
            call gdErrorHandler('DeactivateInequalityConstraints: ' // & 
                'implementation bug: did not get all dofgroups with ' //& 
                'exactly one constraint')
        end if
        
        ! Check if there are any remaining constraints
        if (.not. any(nac)) then 
            ! None left, we can quit here
            return
        end if
        
        ! Attribute constraints
        !======================
        ! Remaining dofgroups with more than 1 dof and active constraints
        i = 1
        do while (i <= ng)
            ! Check 
            if ((nacg(i) > 0) .and. (dofs(i) > 1)) then 

                ! Check
                if (count(nac(dofgroups(i)%cons)) /= nacg(i)) then
                    ! Inconsistency
                    call gdErrorHandler('DeactivateInequalityConstraints: ' // &
                        'encountered bug')
                end if 

                ! Get all constraints
                allocate(Hind(nacg(i)), Hval(nacg(i)))
                Hind = pack(dofgroups(i)%cons, nac(dofgroups(i)%cons))
                Hval = H(Hind) 
                
                ! Get maximal value, attribute
                Hmaxval = maxval(Hval)
                Hmaxind = maxloc(Hval, 1)
                
                ! Sanity check
                if (Hmaxval < 0) then 
                    call gdErrorHandler('DeactiveInequalityConstraints: ' // &
                        'implementation bug: maximal active constraint ' // &
                        'value is strictly negative')
                end if
                
                ! Extract index 
                thisHind = Hind(Hmaxind)
                
                ! Attribute constraint
                nac(thisHind) = .false.
                
                ! Update dofs
                dofs(i) = dofs(i) - 1
                
                ! Update constraint counter
                tg = congroups(thisHind)%dofgroups
                nacg(tg) = nacg(tg) - 1
                
                ! Reset counter
                i = 1

                ! Housekeeping
                deallocate(Hind, Hval)
            else
                ! Update counter
                i = i + 1
            end if
                
        end do
        
        ! Sanity check
        if (any( (nacg > 0) .and. (dofs > 1))) then 
            call gdErrorHandler('DeactivateInequalityConstraints: ' // & 
                'implementation bug: did not get all dofgroups with ' // &
                'multiple constraints')
        end if
        
        ! Check if there are any remaining constraints
        if (.not. any(nac)) then 
            ! None left, we can quit here
            return
        end if
        
        ! Attribute constraints
        !======================
        ! Normally only dofgroups with active constraints with only one dof left
        i = 1
        do while (i <= ng)
            if ((nacg(i) > 0) .and. (dofs(i) == 1)) then 

                ! Check
                if (count(nac(dofgroups(i)%cons)) /= nacg(i)) then
                    ! Inconsistency
                    call gdErrorHandler('DeactivateInequalityConstraints: ' // &
                        'encountered bug')
                end if 

                ! Get all constraints
                allocate(Hind(nacg(i)), Hval(nacg(i)))
                Hind = pack(dofgroups(i)%cons, nac(dofgroups(i)%cons))
                Hval = H(Hind) 
                
                ! Get maximal value, attribute
                Hmaxval = maxval(Hval)
                Hmaxind = maxloc(Hval, 1)
                
                ! Sanity check
                if (Hmaxval < 0) then 
                    call gdErrorHandler('DeactiveInequalityConstraints: ' // &
                        'implementation bug: maximal active constraint ' // & 
                        'value is strictly negative')
                end if
                
                ! Extract index 
                thisHind = Hind(Hmaxind)
                
                ! Attribute constraint
                nac(thisHind) = .false.
                
                ! Update dofs
                dofs(i) = dofs(i) - 1
                
                ! Update constraint counter
                tg = congroups(thisHind)%dofgroups
                nacg(tg) = nacg(tg) - 1
                i = 1

                ! Housekeeping
                deallocate(Hind, Hval)
            else
                i = i + 1
            end if 
        end do 
        
        ! Final sanity checks
        if (any((nacg > 0) .and. (dofs > 0))) then 
            call gdErrorHandler('DeactiveInequalityConstraints: implementation bug: did not maximally attribute constraints')
        end if
        
        ! Deactivate constraints
        !=======================
        ! Set any remaining active constraint values to -inf
        if (count(nac) > 0) then 
            print *, 'DeactivateInequalityConstraints: ', count(nac), ' deactivated, ', count(A)-count(nac), ' remain'
            where (nac) H = -posinfval_R8()
            
        end if

        ! Housekeeping
        !=============
        end associate
               

    end subroutine

    ! Derivative calculation wrapper
    subroutine EvaluateJacobianGD(problem, var, values, jac)

        ! Description
        !============
        ! This routine computes the Jacobian of the optimization problem
        ! with respect to some predefined variable. Note that it is 
        ! assumed that all fields of the problem are up to date. The 
        ! variable(s) to differentiate to should be specified by the
        ! string 'var'. The Jacobian is returned in a sparse matrix 
        ! format and contains the linearization of the optimization 
        ! problem residuals w.r.t. the chosen variables. 

        ! Notes
        !======
        ! Note 1: this routine is a wrapper routine for specific 
        ! derivation computations which are defined in separate 
        ! routines. 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemGDUDT)     :: problem 
        character(*), intent(in)            :: var 
        real(R8), intent(in)                :: values(:)
        type(MySparseUDT), intent(out)      :: jac 

        ! Compute
        !========
        ! Check the variable name and call differentiation routine
        select case (var)

        case ('coordinates')

            ! Differentiate w.r.t. coordinates of the grid (x, y)
            call problem%EvaluateJacobianGDGoatvariables(jac, 'coordinates')

        case ('coordinates_flux')

            ! Differentiate w.r.t. coordinates of the grid (x, y) and 
            ! flux values that were optimized for (psi) 
            call problem%EvaluateJacobianGDGoatvariables(jac, 'coordinates_flux')

        case ('goatvariables')

            ! Differentiate w.r.t. all goat variables (design variables,
            ! lagrange multipliers)
            call problem%EvaluateJacobianGDGoatvariables(jac, 'all')

        case ('vesselcoordinates')

            ! Differentiate w.r.t. vessel coordinates
            call problem%EvaluateJacobianGDVesselcoordinates(values, jac)

        case default 

            ! Throw error, not implemented
            call gdErrorHandler('ComputeJacobianGD: variable not implemented')

        end select


    end subroutine

    ! Derivative calculation w.r.t. coordinates
    subroutine EvaluateJacobianGDGoatvariables(problem, jac, var)

        ! Description
        !============
        ! Evaluate the Jacobian of the optimization problem w.r.t. the 
        ! grid coordinates. This is done by calling the cost function, 
        ! constraint, and other evaluation routines. Then, only columns
        ! that correspond to coordinates are retained from the 
        ! computed linearization. 

        implicit none

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemGDUDT)         :: problem 
        type(MySparseUDT), intent(inout)        :: jac 
        character(*), intent(in)                :: var

        ! Loop
        integer(I8)                             :: k 

        ! Auxiliary
        ! Cost function 
        real(R8)                    :: J 
        real(R8), allocatable       :: gradJ(:)
        type(MySparseUDT)           :: hessJ 

        ! Equality constraints 
        real(R8), allocatable       :: G(:)
        type(MySparseUDT)           :: gradG, hessG  

        ! Inequality constraints 
        real(R8), allocatable       :: H(:), ncp(:)
        type(MySparseUDT)           :: gradH, hessH, gradncpmu, &
            gradncpphi   

        ! Lagrangian 
        real(R8)                    :: L 
        real(R8), allocatable       :: gradL(:), rhs(:)
        type(MySparseUDT)           :: hessL, lhs

        ! Logicals
        logical                     :: dogradient, dohessian
        logical, allocatable        :: tempind(:), I(:)

        ! Other
        integer(I8)                 :: nphi, neq, nineq
        integer(I8), allocatable    :: cind(:)

        type(OptimizationSolverKKTUDT)  :: kktsolver
        type(NumNCPUDT)                 :: num

        ! Initialize
        !===========
        ! Set numerics
        ! Set numerics
        num%ncpfun = 'max'
        num%alpha = 1

        ! Get problem dimensions
        call problem%GetProblemDimensions(nphi, neq, nineq)

        ! Cost function 
        allocate(gradJ(nphi))
        J = 0
        gradJ(:) = 0
        hessJ%nrow = nphi 
        hessJ%ncol = nphi 

        ! Equality constraint 
        allocate(G(neq))
        G(:) = 0
        gradG%nrow = nphi 
        gradG%ncol = neq 
        hessG%nrow = nphi 
        hessG%ncol = nphi

        ! Inequality constraint 
        allocate(H(nineq), ncp(nineq), I(nineq))
        H(:) = 0
        ncp = 0
        I = .true.
        gradH%nrow = nphi 
        gradH%ncol = nineq 
        hessH%nrow = nphi 
        hessH%ncol = nphi

        ! Lagrangian 
        allocate(gradL(nphi + neq + nineq))
        hessL%nrow = nphi + neq + nineq
        hessL%ncol = nphi + neq + nineq
        L = 0
        gradL(:) = 0
        rhs = gradL

        ! Set logicals
        dogradient  = .true. 
        dohessian   = .true.

        ! Associate
        associate(& 
            lambda          =>      problem%lambda, &
            mu              =>      problem%mu, &
            A               =>      problem%A)

        ! Evaluate
        !=========
        ! Update the optimization problem 
        call problem%UpdateProblem()

        ! Evaluate cost function
        call problem%EvaluateCostFunction(J, gradJ, & 
            hessJ, dogradient, dohessian)
        
        ! Evaluate the equality constraints
        call problem%EvaluateEqualityConstraints(G, gradG, &
            hessG, dogradient, dohessian, problem%lambda)
        
        ! Evaluate the inequality constraints
        call problem%EvaluateInequalityConstraints(H, gradH, &
            hessH, dogradient, dohessian, problem%mu)

        ! Nonlinear complementarity function
        call EvaluateNCPfunction(ncp, A, I, gradncpphi, gradncpmu, &
            H, gradH, problem%mu, num, dogradient)

        ! Evaluate lagrangian
        call problem%EvaluateLagrangian(L, gradL, hessL, J, gradJ, &
            hessJ, G, gradG, hessG, problem%lambda, H, gradH, hessH, problem%mu, A, &
            dogradient, dohessian)

        ! Build goat hessian & rhs
        call kktsolver%SetupCorrectionEquation(lhs, rhs, &
            gradJ, hessJ, G, gradG, hessG, problem%lambda, H, gradH, hessH, problem%mu, A, &
            ncp, gradncpphi, gradncpmu)

        ! Extract linearization
        !======================
        ! Get indices
        associate(designvariables  =>   problem%designvariables)
        select case (var)

        case ('coordinates')

            ! Only w.r.t. coordinates
            select type (designvariables)

            type is (DesignVariablesCoordinatesUDT)

                cind = [designvariables%xind, designvariables%yind]

            type is (DesignVariablesCoordinatesFluxUDT)

                cind = [designvariables%xind, designvariables%yind]

            class default

            call gdErrorHandler('EvaluateJacobianGDCoordinates: design variable type not implemented')

            end select
            

        case ('coordinates_flux')

            ! w.r.t. coordinates and flux values (only possible if they 
            ! are optimized)
            select type (designvariables)

            type is (DesignVariablesCoordinatesUDT)

                ! Throw error
                call gdErrorHandler('EvaluateJacobianGDGoatvariables: ' // & 
                    'cannot return derivatives w.r.t. flux as they are not ' // & 
                    'optimized for')

            type is (DesignVariablesCoordinatesFluxUDT)

                cind = [designvariables%xind, designvariables%yind, designvariables%psiind]

            class default

            call gdErrorHandler('EvaluateJacobianGDCoordinates: design variable type not implemented')

            end select

        case ('all')

            ! w.r.t. all variables, so including lagrange multipliers
            cind = [(k, k = 1, nphi+neq+nineq)]

        case default 

        end select
        end associate

        ! Extract linearization
        allocate(tempind(nphi+neq+nineq))
        tempind = .true.
        tempind(cind) = .false.
        jac = lhs%DeleteColumns(tempind)

        ! Housekeeping
        !=============
        end associate


    end subroutine

    ! Derivative calculation w.r.t. vessel coordinates
    subroutine EvaluateJacobianGDVesselcoordinates(problem, values, jac)

        ! Description
        !============
        ! Evaluate the jacobian of the optimization problem with 
        ! respect to the vessel polygon coordinates xv, yv. It is 
        ! assumed that the vessel polygon and its fields are up to date,
        ! and that the optimization problem is also fully solved and
        ! up to date. 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemGDUDT)     :: problem 
        real(R8), intent(in)                :: values(:)
        type(MySparseUDT)                   :: jac

        ! Auxiliary
        logical                             :: dogradient, dohessian
        logical, allocatable                :: A(:), I(:)
        real(R8)                            :: J
        real(R8), allocatable               :: gradJ(:), G(:), H(:), &
            dJdvar(:), ncp(:)
        type(MySparseUDT)                   :: hessJ, gradG, gradH, &
            hessG, hessH, dgradJdvar, dGdvar, dgradGdvar, dHdvar, &
            dgradHdvar, gradncpphi, gradncpmu
        type(NumNCPUDT)         :: num

        ! Initialize
        !===========
        ! Associate
        associate(&
            nphi        =>      problem%designvariables%nphi,   &
            neqcon      =>      problem%constraints%eqcon%neqcon, &
            nineqcon    =>      problem%constraints%ineqcon%nineqcon, &
            lambda      =>      problem%lambda, &
            mu          =>      problem%mu,     &
            grid        =>      problem%grid,   &
            magneticField   =>  problem%magneticField,  &
            environment     =>  problem%environment,    &
            designvariables =>  problem%designvariables)

        ! Allocate
        allocate(gradJ(nphi), G(neqcon), H(nineqcon), ncp(nineqcon), &
            A(nineqcon), I(nineqcon))
        
        ! Initialize - to be cleaned up...
        num%ncpfun = 'max'
        num%alpha = 1
        
        ! Evaluate
        !=========
        ! We don't need hessian evaluation w.r.t. coordinates
        dogradient  = .true. 
        dohessian   = .false.

        ! Cost function (gradient) sensitivities
        call problem%costfunction%Evaluate(J, gradJ, hessJ, grid, &
            magneticField, environment, dogradient, dohessian, &
            designvariables, 'vesselcoordinates', values, dJdvar, dgradJdvar)

        ! lambda*gradG and G sensitivities 
        call problem%constraints%eqcon%Evaluate(G, gradG, hessG, grid, &
            magneticField, environment, dogradient, dohessian, &
            designvariables, lambda, 'vesselcoordinates', values, dGdvar, dgradGdvar)

        ! mu*gradH and H sensitivities
        call problem%constraints%ineqcon%Evaluate(H, gradH, hessH, grid, &
            magneticField, environment, dogradient, dohessian, &
            designvariables, mu, 'vesselcoordinates', values, dHdvar, dgradHdvar)

        ! Nonlinear complementarity function
        call EvaluateNCPfunction(ncp, A, I, gradncpphi, gradncpmu, &
            H, gradH, mu, num, dogradient)

        ! Extract
        !========
        ! Full linearization of residuals should be:
        ! [hess + hessG + hessH, jacG, jacH] (nphi+neq+nineq-by-nvalues)
        jac = SpZeros(0, size(values, 1))

        ! 'Cost function' contribution
        jac = jac%Concatenate(dgradJdvar + dgradGdvar + dgradHdvar, 1)

        ! Equality constraint contribution
        jac = jac%Concatenate(dGdvar, 1)

        ! Inequality constraint contributions (need to set zero for inactive constraints)
        dHdvar = dHdvar%SetZeroRows(.not. A)
        jac = jac%Concatenate(dHdvar, 1)

        ! Housekeeping
        !=============
        end associate

    end subroutine

    !------------------------------------------------------------------!
    !                              AUXILIARY                           !
    !------------------------------------------------------------------!

    ! Routine for dof extraction when constructing inequality constraints
    subroutine SubtractDofsReturnAttributed(IDv, dofs, fslegal, rllegal, vertlegal, hasfs, hasrl, &
        groupindrl, groupindfs, vertrl, vertfs, nvertrl, throwerror, attributed, &
        dofgroupID)

        ! Description
        !============
        ! Remove a dof for each of the vertices specified in ID in the following
        ! way:
        ! - if vertex has dofs left -> do on vertex
        ! - if no dofs left on vertex:
        !   * elseif vertex belongs to rl, rllegal is true and at least one rl has 
        !   * dofs left -> put on rl
        !   * if vertex belongs to fs, fslegal is true, and fs of vert has dof left
        !   -> put on fs
        !   * else issue error

        ! The group ID of the attributed dof for each vertex is returned
        ! in 'dofgroupID' and has the same dimensions as IDv, i.e. the 
        ! vertices from which a dof needs to be subtracted. 
        
        ! Notes
        !======
        ! Note: if a point has multiple radial lines, this has to be an x-point. In
        ! that case, if it is attributed to a radial line, all other radial lines
        ! lose their degrees of freedom as well! 

        ! 

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)         :: IDv(:), groupindrl(:), &
            groupindfs(:), vertrl(:, :), vertfs(:), nvertrl(:)
        integer(I8), intent(inout)      :: dofs(:)
        logical, intent(in)             :: fslegal, rllegal, hasfs(:), &
            hasrl(:), throwerror, vertlegal
        logical, intent(out)            :: attributed(size(IDv, 1))
        integer(I8), intent(out)        :: dofgroupID(size(IDv, 1))

        ! Auxiliary
        integer(I8)                     :: ID, nrl, fsID 
        integer(I8), allocatable        :: rlID(:)

        ! Loop
        integer(I8)                     :: i, k 
        
        ! Initialize
        !===========
        attributed  = .false.
        dofgroupID     = 0
        
        ! Attribute
        !==========
        do i = 1, size(IDv, 1)
            ! Unpack ID
            ID = IDv(i)
            
            ! Subtract dof from vertex?
            if (vertlegal) then 
                if (dofs(ID) > 0) then  
                    ! Vertex has dof left
                    dofs(ID) = dofs(ID)-1
                    attributed(i) = .true.
                    dofgroupID(i) = ID
                end if 
            end if
            
            ! Subtract dof from rl? 
            if (rllegal .and. (.not. attributed(i))) then 
                if (hasrl(ID)) then 
                    ! Unpack radial line data
                    nrl = nvertrl(ID)
                    rlID = groupindrl(vertrl(ID, 1:nrl))

                    ! Set counter
                    k = 1
                    do while ((.not. attributed(i)) .and. (k <= nrl))
                        if (dofs(rlID(k)) > 0) then 
                            dofs(rlID(k)) = dofs(rlID(k)) - 1
                            attributed(i) = .true.
                            dofgroupID(i) = rlID(k)
                        end if
                        k = k + 1
                    end do 
                end if 
            end if
            
            ! Subtract dof from fs? 
            if (fslegal .and. (.not. attributed(i))) then 
                if (hasfs(ID)) then 
                    fsID = groupindfs(vertfs(ID))
                    if (dofs(fsID) > 0) then 
                        dofs(fsID) = dofs(fsID) - 1
                        attributed(i) = .true.
                        dofgroupID(i) = fsID
                    end if 
                end if 
            end if
            
            ! Check
            if ((.not. attributed(i)) .and. (throwerror)) then 
                call gdErrorHandler('ConstructInequalityConstraintGroups: ' // &
                    'could not subtract dof, problem is overly constrained')
            end if  
        end do 
    
    end subroutine

    subroutine SubtractDofsDontReturnAttributed(IDv, dofs, fslegal, rllegal, vertlegal, hasfs, hasrl, &
        groupindrl, groupindfs, vertrl, vertfs, nvertrl, throwerror)

        ! Description
        !============
        ! Same as previous, but now we don't return attributed and 
        ! group ID. Useful for some applications

        ! Arguments
        integer(I8), intent(in)         :: IDv(:), groupindrl(:), &
            groupindfs(:), vertrl(:, :), vertfs(:), nvertrl(:)
        integer(I8), intent(inout)      :: dofs(:)
        logical, intent(in)             :: fslegal, rllegal, hasfs(:), &
            hasrl(:), throwerror, vertlegal

        ! Auxiliary
        logical                         :: attributed(size(IDv, 1))
        integer(I8)                     :: dofgroupID(size(IDv, 1))

        ! Call
        !=====
        call SubtractDofs(IDv, dofs, fslegal, rllegal, vertlegal, hasfs, hasrl, &
            groupindrl, groupindfs, vertrl, vertfs, nvertrl, throwerror, &
            attributed, dofgroupID)

    end subroutine

    ! Routine for attributing a single constraint and subtracting dofs
    subroutine AttributeConstraint(IDv, dofs, fslegal, rllegal, vertlegal, hasfs, hasrl, &
        groupindrl, groupindfs, vertrl, vertfs, nvertrl, dofgroupID, isattributed, throwerror)

        ! Description
        !============
        ! Attribute a constraint that has multiple vertices involved to a single
        ! vertex and remove a dof from that vertex. To attribute the constraint, it
        ! is first checked if any vertex still has own dofs. If that's not the
        ! case, then the radial lines and flux surfaces are checked (if allowed). 

        ! Notes
        !======
        ! Note 1: this routine uses the SubtractDofs routine to determine
        ! to which dof a constraint can be attributed to
        
        ! Note 2: this routine should only be used if it is certain that
        ! the constraints will be active (i.e. primarily for equality 
        ! constraints), otherwise inequality constraint attribution may
        ! be overly conservative 

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)         :: IDv(:), groupindrl(:), &
            groupindfs(:), vertrl(:, :), vertfs(:), nvertrl(:)
        integer(I8), intent(inout)      :: dofs(:)
        logical, intent(in)             :: fslegal, rllegal, hasfs(:), &
            hasrl(:), vertlegal, throwerror
        logical, intent(out)            :: isattributed
        integer(I8)                     :: dofgroupID

        ! Auxiliary
        logical                         :: success, tfslegal, trllegal
        logical, allocatable            :: attributed(:)

        integer(I8), allocatable        :: tdofgroupID(:)
        
        ! Loop
        integer(I8)                     :: k 

        
        ! Initialize
        !===========
        allocate(attributed(1), tdofgroupID(1)) ! always a single constraint
        attributed = .false.
        success = .false.
        isattributed = .false. 
        
        ! Attribute
        !==========
        ! Try attributing to vertices only
        tfslegal = .false.
        trllegal = .false. 
        k = 1
        do while ((.not. success) .and. (k <= size(IDv, 1)))
            call SubtractDofs(spread(IDv(k), 1, 1), dofs, tfslegal, trllegal, vertlegal, hasfs, hasrl, &
                groupindrl, groupindfs, vertrl, vertfs, nvertrl, .false., attributed, &
                tdofgroupID)
            if (all(attributed)) then 
                success = .true.
                dofgroupID = tdofgroupID(1)
                isattributed = .true.
            else
                k = k + 1
                attributed = .false.
            end if 
            
        end do 
        
        
        ! Exit on success
        if (success) then    
            return
        end if
        
        ! Try attributing to radial lines or flux surfaces
        k = 1
        do while (( .not. success) .and. (k <= size(IDv, 1)))
            call SubtractDofs(spread(IDv(k), 1, 1), dofs, fslegal, rllegal, vertlegal, hasfs, hasrl, &
                groupindrl, groupindfs, vertrl, vertfs, nvertrl, .false., attributed, &
                tdofgroupID)
            if (all(attributed)) then 
                success = .true.
                dofgroupID = tdofgroupID(1)
                isattributed = .true.
            else
                k = k + 1
                attributed = .false.
            end if 
            
        end do 
        
        ! Check
        if ( (.not. success) .and. (throwerror)) then 
            call gdErrorHandler('ConstructInequalityConstraintGroups: could not attribute constraint to group')
        end if

        ! Housekeeping
        !=============
        deallocate(attributed, tdofgroupID)

    end subroutine

    ! Routine for attributing constraints to multiple dofs but not 
    ! subtracting them (for inequalities basically)
    subroutine DetermineConstraintDofgroups(IDv, dofs, fslegal, rllegal, vertlegal, hasfs, hasrl, &
        groupindrl, groupindfs, vertrl, vertfs, nvertrl, attributed, dofgroupID)

        ! Description
        !============
        ! Determine which  dofs can be used to attribute a certain 
        ! constraint to. Primarily useful for inequality constraints,
        ! for which it is a priori not known to which constraint they 
        ! should be attributed. IDv is supposed to be all vertices to 
        ! which this constraint can be attributed. 'attributed' is true
        ! if we could attribute the constraints, it is false if not. 
        ! The 'dofgroupID' field is then allocated to hold all group dofs 
        ! to which the constraint can be attributed to. 

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)         :: IDv(:), groupindrl(:), &
            groupindfs(:), vertrl(:, :), vertfs(:), nvertrl(:)
        integer(I8), intent(inout)      :: dofs(:)
        logical, intent(in)             :: fslegal, rllegal, hasfs(:), &
            hasrl(:), vertlegal
        logical, intent(out)            :: attributed
        integer(I8), intent(out), allocatable   :: dofgroupID(:)

        ! Auxiliary
        integer(I8)                     :: ID, nrl, fsID, nID, nIDv
        integer(I8), allocatable        :: rlID(:)

        ! Loop
        integer(I8)                     :: i, k 

        ! Initialize
        !===========
        attributed = .false. 
        nID = 0
        nIDv = size(IDv, 1)
        allocate(dofgroupID(3*nIDv)) ! too big, trim later
        dofgroupID = 0

        ! Determine all groups
        !=====================
        ! Loop over vertices
        do i = 1, nIDv
            ! Unpack ID
            ID = IDv(i)
            
            ! Subtract dof from vertex?
            if (vertlegal) then 
                if (dofs(ID) > 0 ) then 
                    ! Vertex has dof left
                    nID             = nID + 1
                    attributed      = .true.
                    dofgroupID(nID) = ID
                end if
            end if
            
            ! Subtract dof from rl? 
            if (rllegal) then 
                if (hasrl(ID)) then
                    nrl     = nvertrl(ID)
                    rlID    = groupindrl(vertrl(ID, 1:nrl))
                    if (all(dofs(rlID) > 0)) then ! .and. (nrl == 1)
                        attributed = .true.
                        do k = 1, nrl
                            nID = nID + 1
                            dofgroupID(nID) = rlID(k)
                        end do 
                    end if 
                end if 
            end if
            
            ! Subtract dof from fs? 
            if (fslegal) then
                if (hasfs(ID)) then 
                    fsID = groupindfs(vertfs(ID))
                    if (dofs(fsID) > 0) then 
                        nID = nID + 1
                        attributed = .true.
                        dofgroupID(nID) = fsID
                    end if
                end if
            end if
            
        end do

        ! Trim
        dofgroupID = dofgroupID(1:nID)

    end subroutine

    


end module
