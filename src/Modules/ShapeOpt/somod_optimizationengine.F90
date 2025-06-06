!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains the implementation of the optimization engine for
! shape optimization (without SOLPS for now). 

! Notes
!======
! Note 1: Descriptions of what the deferred procedures should do are 
! provided in the abstract interface. 

! Note 2: if shape optimization 

module somod_optimizationengine
    
    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use optmod_optimizationengine
    use somod_designvariables
    use somod_costfunction
    use somod_constraints
    use gdmod_optimizationengine
    use gdmod_designvariables
    use gdmod_costfunction
    use gdmod_constraints
    use gdmod_types

    ! Load SOLPS-specific modules if needed
#ifdef SOLPS 
    use sosmod_costfunction, only : CostfunctionGSRUDT
#endif

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
    type, extends(OptimizationProblemUDT) :: OptimizationProblemSOUDT

        ! Description
        !============
        ! Overwrite the initial design routines with our own
        ! implemented and derived types. 

        ! The classical derived types
        class(DesignVariablesSOUDT), allocatable :: designvariables
        class(CostfunctionSOUDT), allocatable    :: costfunction
        type(ConstraintsSOUDT)                   :: constraints

        ! Additional fields that are needed
        type(OptimizationProblemGDUDT)      :: goat 
        type(DesignOptionsSOUDT)            :: designoptions 
        type(CGStructureUDT), allocatable   :: congroups(:) 
        type(DOFGStructureUDT), allocatable :: dofgroups(:)
        logical                             :: doremesh 
        
    contains

        ! Overwrite 
        !==========
        ! Problem initialization
        procedure :: Initialize      => InitializeOptimizationProblemSO

        ! Dimension query
        procedure :: GetProblemDimensions   => GetProblemDimensionsSO 

        ! Design variables query
        procedure :: GetProblemDesignVariables => &
            GetProblemDesignVariablesSO

        ! Design update
        procedure :: UpdateDesign           => UpdateDesignSO

        ! Problem update
        procedure :: UpdateProblem          => UpdateProblemSO

        ! Cost function evaluation
        procedure :: EvaluateCostFunction   => EvaluateCostFunctionSO

        ! Equality constraints evaluation
        procedure :: EvaluateEqualityConstraints    &
                        => EvaluateEqualityConstraintsSO

        ! Inequality constraints evaluation
        procedure :: EvaluateInequalityConstraints &
                        => EvaluateInequalityConstraintsSO

        ! KKT relaxation
        procedure :: RelaxProblemKKTSystem => RelaxProblemKKTSystemSO

        ! Data output
        procedure :: WriteIterationData => WriteIterationDataSO 

        ! Additional routines
        !====================
        ! Initialization finalizer to account for cross-design/cfv/con
        ! initialization requirements
        procedure :: FinalizeInitialization => FinalizeInitializationSO  

    end type 

    ! Optimization engine
    !====================
    type, extends(OptimizationEngineUDT) :: OptimizationEngineSOUDT 


    contains 

        ! Procedure for initialization of engine
        procedure :: SetupOptimizationDriver => SetupOptimizationDriverSO

    end type

    ! Abstract types
    !===============

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

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
    subroutine SetupOptimizationDriverSO(optimizationdriver) 

        ! Description
        !============
        ! This routine should be used to initialize the solver and
        ! problem to the desired type. Further initialization is done
        ! through the problem and solver initialization routines. 

        ! The usual
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationEngineSOUDT)      :: optimizationdriver 

        ! Auxiliary
        type(OptimizationProblemSOUDT)      :: thisproblem
        type(OptimizationSolverKKTUDT)      :: thissolver

        ! Data

        ! Initialize
        !===========
        ! Allocate the problem type
        allocate(optimizationdriver%problem, source=thisproblem)

        ! Allocate the solver
        allocate(optimizationdriver%solver, source=thissolver)

        ! Propagate paths
        optimizationdriver%problem%inputfilepath = optimizationdriver%inputfilepath 
        optimizationdriver%solver%inputfilepath = optimizationdriver%inputfilepath
        optimizationdriver%solver%inputfileprefix = optimizationdriver%inputfileprefix

    end subroutine

    !------------------------------------------------------------------!
    !                       OPTIMIZATION PROBLEM                       !
    !------------------------------------------------------------------!
    ! Dimension query
    subroutine GetProblemDimensionsSO(problem, nphi, neq, nineq)

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
        class(OptimizationProblemSOUDT)     :: problem 
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
    subroutine GetProblemDesignVariablesSO(problem, phi)

        ! Description
        !============
        ! This routine returns the current design variable vector phi

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemSOUDT)         :: problem
        real(R8), allocatable, intent(out)      :: phi(:)

        ! Get!
        !=====
        phi = problem%designvariables%phi

    end subroutine

    ! Optimization problem initialization
    subroutine InitializeOptimizationProblemSO(problem) 

        ! Description
        !============
        ! This routine further initializes the design variables, cost
        ! function, and constraints. It is assumed that the auxiliary
        ! problem structures have been initialized beforehand. 

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemSOUDT)                 :: problem

        ! Auxiliary
        type(OptimizationEngineGDUDT)      :: goatengine  
        type(ShapeOptimizationOptionsUDT)  :: SOoptions

        ! Data

        ! Filepaths
        !==========
        ! Shape optimization problem filepaths
        SOoptions%inputfilepath = problem%inputfilepath

        ! Design options
        problem%designoptions%inputfilepath = problem%inputfilepath

        ! Remeshing
        !==========
        ! Initialize to false
        problem%doremesh = .false. 

        ! Design options
        !===============
        ! Set shape optimization options
        call SOoptions%Set()

        ! Set the design options
        call problem%designoptions%Set()

        ! GOAT
        !=====
        ! Initialize the GOAT optimization driver
        call GDinitialize(SOoptions%goatfilepath, goatengine)

        ! Assign problem
        associate(goat      => goatengine%problem)
        select type(goat)

        type is (OptimizationProblemGDUDT)

            problem%goat = goat

        class default

            call gdErrorHandler('InitializeOptimizationProblemSO: ' // & 
                'goat optimization problem type not supported')

        end select
        end associate

        ! Design variables
        !=================
        ! Allocate the design variables, depending on the type.
        select case (trim(problem%designoptions%variables%type))

        case ('vesselcoordinates')

            ! Vessel coordinates
            allocate(DesignVariablesVesselCoordinatesUDT::problem%designvariables)

        case ('vesselcoordinates_goat')

            ! Vessel coordinates & goat variables
            allocate(DesignVariablesVesselCoordinatesGoatUDT::problem%designvariables)

        case default

            ! Throw error
            call gdErrorHandler('Unknown design variable type')

        end select

        ! Initialize the design variables
        call problem%designvariables%Initialize(problem%goat)
        
        ! Cost function
        !==============
        ! Check some options
        if (.not. problem%designoptions%costfunction%dogoatreduction) then 
            ! Check if SOLPS is active - can only be active for reduced
            ! goat approach
            if (problem%designoptions%costfunction%includesolps) then 
                call gdErrorHandler('InitializeOptimizationProblemSO: '// & 
                    'solps contribution only allowed if goat is used in ' // & 
                    'reduced way')
            end if 
        end if 

        ! Allocate the cost function, depending on the type
        if (problem%designoptions%costfunction%dogoatreduction) then 

            ! Check if it's also SOLPS-based and if this is allowed
            if (problem%designoptions%costfunction%includesolps) then 
#ifdef SOLPS 
                ! Allowed
                allocate(CostFunctionGSRUDT::problem%costfunction)
#else 
                ! Not allowed, throw error
                call gdErrorHandler('InitializeOptimizationProblemSO: ' // & 
                    'cannot use solps-reduced cost function, since ' // & 
                    'solps is not available')
#endif
            else

                ! Goat-reduced cost function - initialization is done 
                ! further in the initialization routine of the cost function
                allocate(CostFunctionGRUDT::problem%costfunction)

            end if 
        else
            ! Classical cost function
            select case (trim(problem%designoptions%costfunction%type))

            case ('PLF')

                ! Allocate 
                allocate(CostFunctionPLFUDT::problem%costfunction)

                ! Set type
                problem%costfunction%type = 'PLF'

            case ('SOFA')

                ! Allocate
                allocate(CostFunctionSOFAUDT::problem%costfunction)

                ! Set type
                problem%costfunction%type = 'SOFA'

            case ('no')

                ! Zero cost function
                allocate(CostFunctionDummyUDT::problem%costfunction)

                ! Set type
                problem%costfunction%type = 'no'

            case ('general')
                
                ! General cost function
                allocate(CostFunctionGeneralSOUDT::problem%costfunction)

                ! Set type
                problem%costfunction%type = 'general'
                
            case default
                
                ! Throw error
                call gdErrorHandler('Unknown cost function type')

            end select
        end if

        ! Initialize the cost function
        call problem%costfunction%Initialize(problem%goat, &
            problem%designoptions%costfunction)

        ! Constraints
        !============
        ! Given the (many) possible options for the constraints, the 
        ! constraints are set in its own initialization. 
        call problem%constraints%Initialize(problem%goat, & 
            problem%designvariables, problem%designoptions%constraints)

        ! Set Lagrange multipliers
        allocate(problem%lambda(problem%constraints%eqcon%neqcon), &
            problem%mu(problem%constraints%ineqcon%nineqcon), &
            problem%A(problem%constraints%ineqcon%nineqcon))
        problem%lambda = 0
        problem%mu = 0
        problem%A = .false.

        ! Checks
        !=======
        ! Check some combinations of cfv/cons/designvariables
        select case (trim(problem%designoptions%variables%type))

        case ('vesselcoordinates')

            ! If goat constraints are active, the cost function should
            ! be reduced
            if (problem%constraints%eqcon%dogoat) then 
                call gdErrorHandler('Design variables "vesselcoordinates" '// & 
                    ' require reduced cost function formulation and ' // &
                    'inactive goat constraints (these are active now)')
            end if 
            

        case ('vesselcoordinates_goat')

            ! Check if goat constraints are on, otherwise throw error
            if (.not. problem%constraints%eqcon%dogoat) then 
                call gdErrorHandler('Design variables "vesselcoordinates_goat" ' // &
                    'require active goat constraints (are currently not active)')
            end if 

            ! Check if we don't have a goat-reduced cost function
            if (problem%designoptions%costfunction%dogoatreduction) then 
                call gdErrorHandler('Design variables "vesselcoordinates_goat" ' // &
                'cannot be used together with a goat-reduced cost function')
            end if 

        case default 

            ! Check cost function
            select case (problem%costfunction%type)

            case ('PLF', 'SOFA', 'general')

                ! Requires  vesselcoordinates(_goat)
                call gdErrorHandler('Cost function "' // problem%costfunction%type // &
                    '" requires design ' // &
                    'variables with vessel coordinates, check input')

            case default 

            end select

            ! Check constraints
            if (problem%constraints%eqcon%dofixedvesselpoints) then 
                ! Requires vesselcoordinates(_goat)
                call gdErrorHandler('"fixedvesselpoints" constraint requires design ' // &
                    'variables with vessel coordinates, check input')
            end if 

        end select

        ! Finalize initialization
        !========================
        ! Initialize design variables further for constraint/cfv 
        ! dependent fields
        call problem%FinalizeInitialization()


    end subroutine

    ! Finalize the problem initialization
    subroutine FinalizeInitializationSO(problem)

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
        class(OptimizationProblemSOUDT)                 :: problem 

        ! Auxiliary

        ! Loop

        ! Initialize
        !===========
        ! Associate

        ! Update cost function
        !=====================
        ! Hessian estimator
        ! Check which hessian estimator to construct
        associate(hopt => problem%designoptions%costfunction%hessapprox)
    
        select case (hopt%inithess)

        case ('diagonal')

            problem%costfunction%B = ConstructHessianApproximation(hopt%updatemethod, &
                problem%designvariables%nphi, hopt%diagind, hopt%diagval, hopt%storagetype)

        case default 

            call gdErrorHandler('FinalizeInitializationSO: ' // & 
                'unknown initial hessian option: "' // hopt%inithess // '"')

        end select
        end associate

        ! Update constraints
        !===================
        ! Construct constraint groups
        ! call problem%ConstructInequalityConstraintGroups()

        ! Write data
        !===========
        ! Write the original/initial vessel polygon structure
        call problem%goat%environment%vessel%polygonset%WriteData('vesselpolygon_orig_so')

        ! Write original grid cells
        call WriteGridCells(problem%goat%grid, 'cells_init')
        call WriteGridVertices(problem%goat%grid, 'vertices_init')

        ! Housekeeping
        !=============

    end subroutine

    ! Problem update
    subroutine UpdateProblemSO(problem)

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
        class(OptimizationProblemSOUDT)     :: problem 

        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Initialize
        !===========
        ! Associate to use select type
        associate(&
            designvariables     => problem%designvariables, &
            constraints         => problem%constraints,     &
            goat                => problem%goat             &
            )

        ! Update design
        !==============
        ! Simply call the update routine from the design variables
        call designvariables%UpdateDesign(goat)

        ! Update goat
        !===========
        ! Assumed coordinates already properly re-initialized during
        ! UpdateDesign
        call goat%UpdateProblem()

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Design update
    subroutine UpdateDesignSO(problem, dx)

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
        class(OptimizationProblemSOUDT)     :: problem 
        real(R8), intent(in)                :: dx(:)
    
        ! Loop variables

        ! Auxiliary variables 
        character(:), allocatable           :: vertpath, cellpath

        ! Data

        ! Update design
        !==============
        ! Simply update designvariables%phi
        problem%designvariables%phi = problem%designvariables%phi + dx 

        ! Debug plots?
        !=============
        ! Call grid vertex writing routine
        allocate(character(len('vertices_iterate')) :: vertpath)
        allocate(character(len('cells_iterate')) :: cellpath)
        vertpath = 'vertices_iterate'
        cellpath = 'cells_iterate'
        call WriteGridVertices(problem%goat%grid, vertpath) 
        call WriteGridCells(problem%goat%grid, cellpath)
        call problem%goat%environment%vessel%polygonset%WriteData('vesselpolygon_iterate_so')

    end subroutine

    ! Cost function evaluation
    subroutine EvaluateCostFunctionSO(problem, J, gradJ, hessJ, &
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
        class(OptimizationProblemSOUDT)     :: problem 
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
        call problem%costfunction%Evaluate(J, gradJ, hessJ, problem%goat, &
            dogradient, dohessian, problem%designvariables)

    end subroutine

    ! Equality constraints evaluation
    subroutine EvaluateEqualityConstraintsSO(problem, G, gradG, hessG, &
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
        class(OptimizationProblemSOUDT)     :: problem 
        real(R8), allocatable               :: G(:), lambda(:)
        type(MySparseUDT)                   :: gradG, hessG
        logical                             :: dogradient, dohessian                            

        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Compute the constraints
        !========================
        ! Simply call the constraint computation routine
        call problem%constraints%eqcon%Evaluate(G, gradG, hessG, &
            problem%goat, dogradient, &
            dohessian, problem%designvariables, lambda)

    end subroutine

    ! Inequality constraints evaluation
    subroutine EvaluateInequalityConstraintsSO(problem, H, gradH, &
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
        class(OptimizationProblemSOUDT)     :: problem 
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
            problem%goat, &
            dogradient, dohessian, problem%designvariables, mu)

        ! Deactivate constraints
        !=======================
        ! Only if necessary - determined in the routine itself
        ! Note: should we add switch here? And adjust gradient/Hessian
        ! as well? At least gradient should be tackled by activeness
        ! of constraint...
        ! call problem%DeactivateInequalityConstraints(H)

    end subroutine

    ! KKT system relaxation
    subroutine RelaxProblemKKTSystemSO(problem, KKT)

        ! Description
        !============
        ! Relax KKT system for shape optimization. It is recommended to
        ! use this routine if the goat is treated explicitly as a 
        ! constraint and therefore the optimization variables include
        ! lagrange multipliers of these constraints. It might make sense
        ! to do a problem-specific approach here, e.g. by applying 
        ! relaxation only on the hessian contributions of the design 
        ! variables in the upper and lower level problems. 

        ! In this implementation, we apply the same relaxation factor 
        ! as obtained from the problem monitor to the shape optimization
        ! and goat problem (this may be adjusted in the future). Only
        ! the cost function contribution to the hessian are retained, 
        ! and only those contributions that differentiate w.r.t. (only)
        ! the actual design variables (i.e. vessel coordinates and 
        ! grid coordinates/psi values)

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemSOUDT)         :: problem 
        type(MySparseUDT)                       :: KKT 

        ! Auxiliary
        integer(I8)                             :: nphi, neq, nineq
        integer(I8), allocatable, dimension(:)  :: phiind

        real(R8), allocatable, dimension(:)     :: diag
        
        type(MySparseUDT)                       :: KKTcopy, rxfmat

        ! Initialize
        !===========
        ! Associate
        associate(&
            rxf         => problem%monitor%rxf)

        ! Get the problem dimensions
        call problem%GetProblemDimensions(nphi, neq, nineq)

        ! Get the indices of the contributions to relax (i.e. shape opt
        ! design variables and, if present, goat design variables)
        associate(dv        => problem%designvariables)
        select type(dv)

        type is (DesignVariablesVesselCoordinatesUDT) 

            ! Only shape optimization design variables to account for
            phiind = [dv%xind, dv%yind]

        type is (DesignVariablesVesselCoordinatesGoatUDT)

            ! Shape and goat design variables
            phiind = [dv%xind, dv%yind, dv%phigoatind]

        class default

            ! Unknown
            call gdErrorHandler('RelaxProblemKKTSystemSO: unknown ' // & 
                'shape optimization design variable type')

        end select
        end associate

        ! Relax
        !======
        ! Compute diagonal entries
        KKTcopy = KKT 
        KKTcopy%val = abs(KKTcopy%val)
        call KKTcopy%SumColumnwiseFull(diag) ! diagonal matrix

        ! Set to one where diag is zero
        where(diag == 0.0) diag = 1.0

        ! Construct relaxation matrix
        rxfmat = ConstructMySparse(phiind, phiind, rxf*diag(phiind), &
            nphi+neq+nineq, nphi+neq+nineq)

        ! Add to KKT matrix 
        KKT = KKT + rxfmat

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Data writing
    subroutine WriteIterationDataSO(problem, itopt)

        ! Description
        !============
        ! Write out the optimization problem data in each optimization
        ! iteration. Here, we call the goat-specific data writing 
        ! routine and additionally write out the vessel polygon 
        ! coordinates.

        ! Modules
        !========
        use mod_plotter, only: plotdir 
        use mod_specialchars, only: filesepchar

        ! Declare variables
        !==================
        ! Arguments
        class(OptimizationProblemSOUDT)     :: problem 
        integer(I8)                         :: itopt 

        ! Auxiliary 
        integer, parameter                      :: fid = 70
        integer                                 :: iostat
        logical                                 :: isfile 
        character(:), allocatable               :: filepath 

        ! Initialize
        !===========
        associate(&
            vessel  => problem%goat%environment%vessel, &
            goat    => problem%goat                     &
            )

        ! Write data
        !===========
        ! Goat
        if (problem%designoptions%costfunction%dogoatreduction) then 
            call goat%WriteIterationData(itopt)
        end if 

        ! Vessel data
        call vessel%polygonset%WriteData('vesselpolygon_iterate_so')

        ! Optimization output
        !====================
        ! Associate some fields for easier reading/writing
        associate(grid      => goat%grid)

        ! Set filepath
        filepath = plotdir // filesepchar // 'so_optimization_history.dat'

        ! Cell vertex data
        call WriteGridCells(grid, 'cells_iterate')

        ! Iteration data
        if (itopt == 1) then 
            ! Check if the file exists
            inquire(file=filepath, exist=isfile)
            if (isfile) then 
                ! Replace old file
                open(unit=fid, status='old', iostat=iostat, file=filepath)
                rewind(fid)
            else
                ! Create new file
                open(unit=fid, status='new', iostat=iostat, file=filepath)
            end if

            ! Write header
            call problem%monitor%WriteFileHeader(fid)
        else
            ! File should already by opened, append
            open(unit=fid, status='old', iostat=iostat, file=filepath, access='sequential', position='append')
        end if

        ! Write
        call problem%monitor%WriteFileIterate(fid)


        ! Close file
        close(fid)
        

        ! Housekeeping
        !=============
        end associate
        end associate

    end subroutine


    


end module
