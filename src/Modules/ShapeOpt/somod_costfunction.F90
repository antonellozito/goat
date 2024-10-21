!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains the cost function implementation. Each cost
! function inherits from the basic cost function type for the grid 
! deformation, which itself inherits from the archetypical cost function 
! type defined in the optmod_costfunction module. The cost functions 
! defined here are specific for shape optimization. 

module somod_costfunction
    
    ! Initialize
    !============
    ! Load modules
    use optmod_costfunction
    use gdmod_optimizationengine
    use somod_userinput
    use somod_designvariables
    use PolygonLevelsetFunction2D
    use mod_linearsolverinterface

    ! The usual
    implicit none
    save
    public 

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    ! Abstract types
    !===============
    ! General cost function type
    type, abstract, extends(CostfunctionUDT) :: CostfunctionSOUDT

        ! Description
        !============
        ! Defines the basic cost function structure for the shape 
        ! optimization. The following general fields are added: 
        ! - J:          The cost function value (scalar)
        ! - type:       the cost function type (string)

        ! The following routines should be implemented for these cost
        ! functions (see also the interface below for a description of
        ! what the routines should do):
        ! - Initialize
        ! - Evaluate

        ! Cost function value
        real(R8)                        :: J 

        ! Cost function type
        character(:), allocatable       :: type

    contains

        ! Cost function initialization
        procedure(InitializeCostfunctionSOINT), deferred      :: Initialize

        ! Cost function evaluation
        procedure(EvaluateCostFunctionSOINT), deferred        :: Evaluate

    end type

    ! Derived types
    !==============
    ! Dummy (zero) cost function (e.g. for reduced formulation with solps)
    type, extends(CostFunctionSOUDT) :: CostFunctionDummyUDT 

        ! Description
        !============
        ! Dummy cost function that simply returns zero  for evaluation
        ! etc.

    contains 

        ! Initialization
        procedure :: Initialize         => InitializeCostFunctionDummy 

        ! Evaluation
        procedure :: Evaluate           => EvaluateCostFunctionDummy

    end type
    
    ! Levelset fitting cost function
    type, extends(CostfunctionSOUDT) :: CostfunctionPLFUDT

        ! Description
        !============
        ! This cost function aims to match the vessel coordinates to 
        ! a certain given levelset function. This levelset function 
        ! should be chosen wisely of course. Currently, the 
        ! initialization of this levelset function is based on an 
        ! additional vessel structure file that is read in and for which
        ! a polygon description is made. The following fields are 
        ! present:
        !
        !   targetplf       : target polygon levelset function
        !   lambda          : scaling parameter

        ! Fields
        class(PolygonLevelsetFunction2DUDT), allocatable    :: targetplf
        real(R8)                                            :: lambda

    contains

        ! Initialization
        procedure :: Initialize             => InitializeCostfunctionPLF

        ! Evaluation
        procedure :: Evaluate               => EvaluateCostFunctionPLF

    end type

    ! General goat-reduced cost function
    type, extends(CostfunctionSOUDT) :: CostFunctionGRUDT

        ! Description 
        !============
        ! General cost function for goat-reduced formulation. It has an
        ! allocatable cost function field which can take any of the 
        ! regular cost function options and the goat solver (the
        ! optimization problem is already saved as state variable 
        ! anyway). This allows for a nested formulation without 
        ! having to change solvers etc (i.e. KKT solver can be reused or
        ! any other if desired for the shape optimization problem). 
        ! Additionally, the full goat engine is added. Since the goat
        ! problem is already passed through to all the routines, the
        ! problem instance of the engine is simply equated to the 
        ! goat problem. 

        ! Note: one should *not* impose the goat constraints in 
        ! addition to the optimization problem 

        ! Fields
        class(CostFunctionSOUDT), allocatable       :: costfunction 
        type(OptimizationEngineGDUDT)               :: goatengine 

    contains 

        ! Initialization
        procedure :: Initialize             => InitializeCostfunctionGR

        ! Evaluation
        procedure :: Evaluate               => EvaluateCostFunctionGR

    end type

    
    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    ! Abstract interfaces
    !====================
    ! Cost function
    abstract interface

        ! Cost function initialization
        subroutine InitializeCostfunctionSOINT(costfunction, goat, options)

            ! Description
            !============
            ! This routine should initialize all additional parameters
            ! that are needed to evaluate the cost function (e.g. the 
            ! vertex indices where the cost function is defined).
            
            ! Import
            import :: CostfunctionSOUDT, OptimizationProblemGDUDT, &
                CostFunctionOptionsSOUDT

            ! Declare 
            class(CostfunctionSOUDT)        :: costfunction 
            type(OptimizationProblemGDUDT)  :: goat
            type(CostFunctionOptionsSOUDT)  :: options

        end subroutine

        ! Cost function evaluation
        subroutine EvaluateCostFunctionSOINT(costfunction, J, gradJ, &
            hessJ, goat, dogradient, dohessian, designvariables, &
            varin, valuesin, dJdvarin, dgradJdvarin)

            ! Description
            !============
            ! Main routine to evaluate the cost function and its 
            ! derivative and hessian w.r.t. design variables. 

            ! Import
            import :: CostfunctionSOUDT, MySparseUDT, R8, &
                OptimizationProblemGDUDT, DesignVariablesSOUDT
            
            ! Declare
            class(CostfunctionSOUDT)        :: costfunction 
            real(R8)                        :: J 
            real(R8), allocatable           :: gradJ(:)
            type(MySparseUDT)               :: hessJ 
            type(OptimizationProblemGDUDT)  :: goat
            logical                         :: dogradient, dohessian
            class(DesignVariablesSOUDT)     :: designvariables

            character(*), intent(in), optional  :: varin 
            real(R8), intent(in), optional      :: valuesin(:)
            real(R8), allocatable, optional     :: dJdvarin(:) 
            type(MySparseUDT), optional         :: dgradJdvarin

        end subroutine

    end interface

    contains

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                        DUMMY COST FUNCTION                       !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCostFunctionDummy(costfunction, goat, options)

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionDummyUDT)         :: costfunction
        type(OptimizationProblemGDUDT)      :: goat
        type(CostFunctionOptionsSOUDT)      :: options

        ! Nothing to do here
        
    end subroutine

    ! Cost function evaluation
    subroutine EvaluateCostFunctionDummy(costfunction, J, gradJ, hessJ, &
        goat, dogradient, dohessian, designvariables, &
        varin, valuesin, dJdvarin, dgradJdvarin)


        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionDummyUDT)     :: costfunction 
        real(R8)                        :: J
        real(R8), allocatable           :: gradJ(:) ! assumed initialized
        type(MySparseUDT)               :: hessJ ! assumed in initialized
        type(OptimizationProblemGDUDT)  :: goat
        logical                         :: dogradient, dohessian 
        class(DesignVariablesSOUDT)     :: designvariables

        ! Optional arguments
        character(*), intent(in), optional  :: varin 
        real(R8), intent(in), optional      :: valuesin(:)
        real(R8), allocatable, optional     :: dJdvarin(:) 
        type(MySparseUDT), optional         :: dgradJdvarin

        character(:), allocatable           :: var
        real(R8), allocatable               :: values(:)
        real(R8), allocatable               :: dJdvar(:) 
        type(MySparseUDT)                   :: dgradJdvar

        ! Initialize
        !===========
        ! Check inputs
        if (present(varin)) then 
            var = varin 
        else
            var = 'no'
        end if 
        if (present(valuesin)) then 
            values = valuesin 
        else
            allocate(values(0))
        end if 

        ! Set outputs to zero
        J = 0
        gradJ = 0
        hessJ = SpZeros(designvariables%nphi, designvariables%nphi)

        ! Other derivatives
        !==================
        ! Initialize
        allocate(dJdvar(size(values)))
        dJdvar = 0
        dgradJdvar = SpZeros(size(gradJ), size(values)) ! jacobian, not gradient

        ! Housekeeping
        !=============
        ! Optional arguments
        if (present(dJdvarin)) then 
            dJdvarin = dJdvar 
        end if 
        if (present(dgradJdvarin)) then 
            dgradJdvarin = dgradJdvar
        end if

    end subroutine

    !------------------------------------------------------------------!
    !                         LEVELSET FUNCTION                        !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCostFunctionPLF(costfunction, goat, options)

        ! Description
        !============
        ! Initialize the cost function and its parameters based on the 
        ! goat. The cost function is defined as:
        ! 
        !   cfv = lambda sum_i 0.5*V(xv(i), yv(i))^2,
        ! 
        ! where V is the polygon levelset function of the desired
        ! vessel location. This plf is set up by reading in another
        ! polygon structure from a file specified in the options and
        ! constructing a vessel from that file. 
        
        ! It is very important to realize that this cost function does 
        ! not map specific vessel vertices to a specific location! 
        ! So results may vary (e.g. vessel vertices may collide, or if 
        ! different closed structures exist, some vessel points may lie
        ! on one or the other)

        ! Modules

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionPLFUDT)           :: costfunction
        type(OptimizationProblemGDUDT)      :: goat
        type(CostFunctionOptionsSOUDT)      :: options

        ! Auxiliary 
        integer(I8)                         :: filespec 
        type(VesselOptionsUDT)              :: vesseloptions
        type(VesselUDT)                     :: newvessel

        ! Data
        
        ! Initialize
        !===========
        ! Set the scaling constant
        costfunction%lambda = options%plf%lambda

        ! Set the vessel options loading path (typically goat filepath)
        vesseloptions%inputfilepath = options%plf%vesselinputfilepath 

        ! Initialize the vessel options
        call vesseloptions%Set()

        ! Set the new vessel filepath
        vesseloptions%filepath = options%plf%newvesselfilepath

        ! Set file specifier
        filespec = 10
 
        ! Set up the vessel polygon
        !==========================
        ! Read in the vessel
        call ReadVessel(filespec, newvessel, vesseloptions)

        ! Extract data
        call ExtractVesselData(newvessel, vesseloptions)

        ! Set PLF
        !========
        costfunction%targetplf = newvessel%plfvessel
        
    end subroutine

    ! Cost function evaluation
    subroutine EvaluateCostFunctionPLF(costfunction, J, gradJ, hessJ, &
        goat, dogradient, dohessian, designvariables, &
        varin, valuesin, dJdvarin, dgradJdvarin)

        ! Description
        !============
        ! Evaluate the cost function, the gradient and its hessian. 
        ! Here, the target plf is evaluated at the vessel coordinates 
        ! (which may or may not be design variables...)

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionPLFUDT)       :: costfunction 
        real(R8)                        :: J
        real(R8), allocatable           :: gradJ(:) ! assumed initialized
        type(MySparseUDT)               :: hessJ ! assumed in initialized
        type(OptimizationProblemGDUDT)  :: goat
        logical                         :: dogradient, dohessian 
        class(DesignVariablesSOUDT)     :: designvariables

        ! Optional arguments
        character(*), intent(in), optional  :: varin 
        real(R8), intent(in), optional      :: valuesin(:)
        real(R8), allocatable, optional     :: dJdvarin(:) 
        type(MySparseUDT), optional         :: dgradJdvarin

        character(:), allocatable           :: var
        real(R8), allocatable               :: values(:)
        real(R8), allocatable               :: dJdvar(:) 
        type(MySparseUDT)                   :: dgradJdvar

        ! Auxiliary
        integer(I8)                             :: nv 
        integer(I8), allocatable, dimension(:)  :: row, col
        real(R8), allocatable, dimension(:)     :: xv, yv, val, valxx, &
            valxy, valyy, dvaldx, dvaldy, d2valdx2, d2valdxdy, d2valdy2

        ! Loop 
        integer(I8)                             :: k 

        ! Initialize
        !===========
        ! Associate
        associate(&
            ps          => goat%environment%vessel%polygonset,    &
            plf         => costfunction%targetplf)

        ! Check inputs
        if (present(varin)) then 
            var = varin 
        else
            var = 'no'
        end if 
        if (present(valuesin)) then 
            values = valuesin 
        else
            allocate(values(0))
        end if 

        ! Get vessel coordinates
        call ps%GetVertices(xv, yv)

        ! Initialize
        nv = size(xv)
        hessJ = SpZeros(designvariables%nphi, designvariables%nphi)

        ! Value
        !======
        ! Simply evaluate the levelset function 
        allocate(val(nv))
        call plf%Evaluate(xv, yv, 0, 0, val)

        ! Compute cost function contribution
        J = 0.5*costfunction%lambda*sum(val**2)

        ! Gradient & Hessian
        !===================
        ! Compute
        select case (designvariables%type)

        case ('vesselcoordinates', 'vesselcoordinates_goat')

            ! Same gradient contributions for both design variables, 
            ! since vessel coordinates come first in design variable 
            ! vector and since there are no contributions in terms of 
            ! goat coordinates

            ! Gradient
            !---------
            if (dogradient) then 
                
                ! Evaluate derivatives
                allocate(dvaldx(nv), dvaldy(nv))
                call plf%Evaluate(xv, yv, 1, 0, dvaldx)
                call plf%Evaluate(xv, yv, 0, 1, dvaldy)

                ! Evaluate contributions
                gradJ(1:nv)         = val*dvaldx 
                gradJ(nv+1:2*nv)    = val*dvaldy

            end if 

            ! Hessian
            !--------
            if (dohessian) then 

                ! Check 
                if (.not. dogradient) then 

                    ! Evaluate derivatives
                    allocate(dvaldx(nv), dvaldy(nv))
                    call plf%Evaluate(xv, yv, 1, 0, dvaldx)
                    call plf%Evaluate(xv, yv, 0, 1, dvaldy)

                end if 

                ! Evaluate second order derivatives
                allocate(d2valdx2(nv), d2valdxdy(nv), d2valdy2(nv))
                call plf%Evaluate(xv, yv, 2, 0, d2valdx2)
                call plf%Evaluate(xv, yv, 1, 1, d2valdxdy)
                call plf%Evaluate(xv, yv, 0, 2, d2valdy2)

                ! Evaluate contributions
                valxx = val*d2valdx2 + dvaldx**2
                valxy = val*d2valdxdy + dvaldx*dvaldy
                valyy = val*d2valdy2 + dvaldy**2

                ! Set hessian
                row = [(k, k = 1, nv), (k, k = 1, nv), &
                    (k, k = nv+1, 2*nv), (k, k = nv+1, 2*nv)]
                col = [(k, k = 1, nv), (k, k = nv+1, 2*nv), &
                    (k, k = 1, nv), (k, k = nv+1, 2*nv)]
                val = [valxx, valxy, valxy, valyy]
                hessJ = ConstructMySparse(row, col, val, designvariables%nphi, designvariables%nphi)

            end if 



        case default 

            ! Throw error
            call gdErrorHandler('EvaluateCostFunctionPLF: gradient ' // & 
                ' and hessian not implemented for design variables: ' // &
                designvariables%type)

        end select

        ! Scale
        !------
        gradJ = gradJ*costfunction%lambda 
        hessJ = hessJ*costfunction%lambda

        ! Other derivatives
        !==================
        ! Initialize
        allocate(dJdvar(size(values)))
        dJdvar = 0
        dgradJdvar = SpZeros(size(gradJ), size(values)) ! jacobian, not gradient

        ! Currently no derivatives w.r.t. any other variables
        select case (var)

        case ('goatvariables', 'no')

        case default 

            call gdErrorHandler('EvaluateCostFunctionPLF: unknown ' // & 
                'variable for derivative calculation: ' // var)

        end select

        ! Housekeeping
        !=============
        ! Optional arguments
        if (present(dJdvarin)) then 
            dJdvarin = dJdvar 
        end if 
        if (present(dgradJdvarin)) then 
            dgradJdvarin = dgradJdvar
        end if

        ! Assocation termination
        end associate

    end subroutine

    !------------------------------------------------------------------!
    !                     GOAT REDUCED COST FUNCTION                   !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCostFunctionGR(costfunction, goat, options)

        ! Description
        !============
        ! Initialize the goat reduced cost function. This basically means
        ! initializing the goat engine and the costfunction type. The 
        ! cost function is then allocated to the correct type and 
        ! initialized here as well. 

        ! Arguments
        class(CostfunctionGRUDT)            :: costfunction
        type(OptimizationProblemGDUDT)      :: goat
        type(CostFunctionOptionsSOUDT)      :: options

        ! Initialize goat engine
        !=======================
        ! Set the input filepaths for the driver
        costfunction%goatengine%inputfilepath = goat%inputfilepath
        costfunction%goatengine%inputfileprefix = 'gd.'
        
        ! Initialize the driver
        call costfunction%goatengine%SetupOptimizationDriver()

        ! Initialize the problem
        costfunction%goatengine%problem = goat

        ! Initialize the solver
        call costfunction%goatengine%solver%Initialize()

        ! Initialize costfunction
        !========================
        ! Check which cost function it is
        select case (options%type)

        case ('PLF')

            ! Allocate
            allocate(CostfunctionPLFUDT::costfunction%costfunction)

        case default 

            ! Throw error
            call gdErrorHandler('InitializeCostFunctionGR: unknown cost' // & 
                'function type: ' // options%type)

        end select

        ! Initialize the cost function further
        call costfunction%costfunction%Initialize(goat, options)

    end subroutine

    ! Evaluation
    subroutine EvaluateCostFunctionGR(costfunction, J, gradJ, hessJ, &
        goat, dogradient, dohessian, designvariables, &
        varin, valuesin, dJdvarin, dgradJdvarin)

        ! Description
        !============
        ! Evaluate the cost function, the gradient and its hessian. 
        ! First, the goat equations are solved by calling the driver (
        ! it is assumed that goat is up to date). Afterwards, the cost
        ! function is evaluated and its gradient is computed by 
        ! applying a discrete adjoint approach to account for the
        ! goat constraints. The contributions of the goat constraints to
        ! the hessian of the problem are not accounted for... 

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionGRUDT)        :: costfunction 
        real(R8)                        :: J
        real(R8), allocatable           :: gradJ(:) ! assumed initialized
        type(MySparseUDT)               :: hessJ ! assumed in initialized
        type(OptimizationProblemGDUDT)  :: goat
        logical                         :: dogradient, dohessian 
        class(DesignVariablesSOUDT)     :: designvariables

        ! Optional arguments
        character(*), intent(in), optional  :: varin 
        real(R8), intent(in), optional      :: valuesin(:)
        real(R8), allocatable, optional     :: dJdvarin(:) 
        type(MySparseUDT), optional         :: dgradJdvarin

        character(:), allocatable           :: var
        real(R8), allocatable               :: values(:)
        real(R8), allocatable               :: dJdvar(:) 
        type(MySparseUDT)                   :: dgradJdvar

        ! Auxiliary
        integer(I8)                                 :: flag

        real(R8), allocatable, dimension(:)         :: goatvariables, &
            gradJR, gradJgoat, lambdaG 

        type(MySparseUDT)                           :: hessJR, &
            jacGgoat, jacGdes, gradGdes, gradGgoat

        ! Initialize
        !===========
        ! Check inputs
        if (present(varin)) then 
            var = varin 
        else
            var = 'no'
        end if 
        if (present(valuesin)) then 
            values = valuesin 
        else
            allocate(values(0))
        end if 

        ! Initialize others
        allocate(gradJR(designvariables%nphi))
        gradJR = 0

        ! Construct goat variables
        goatvariables = [goat%designvariables%phi, goat%lambda, goat%mu]

        ! Initialize others
        allocate(gradJgoat(size(goatvariables)))

        ! Solve goat 
        !===========
        associate(goatproblem       => costfunction%goatengine%problem)
        select type (goatproblem)

        type is (OptimizationProblemGDUDT)

            ! Update the problem with the latest goat
            goatproblem = goat 

            ! Call the driver
            call costfunction%goatengine%solver%SolveOptimizationProblem(goatproblem)

            ! Update goat
            goat = goatproblem

        class default 

            ! Throw error
            call gdErrorHandler('EvaluateCostFunctionGR: unexpected goat type')

        end select
        end associate

        ! Compute reduced cost function
        !==============================
        ! Call evaluation routine
        call costfunction%costfunction%Evaluate(J, gradJR, hessJR, &
            goat, dogradient, dohessian, designvariables, &
            'goatvariables', goatvariables, gradJgoat)

        ! Compute gradient 
        !=================
        ! Evaluate goat jacobian w.r.t. goat variables
        call goat%EvaluateJacobian('goatvariables', goatvariables, jacGgoat)

        ! Evaluate goat jacobian w.r.t. design variables
        select case (designvariables%type)

        case ('vesselcoordinates')

            ! Derivatives w.r.t. vessel coordinate
            call goat%EvaluateJacobian('vesselcoordinates', &
                designvariables%phi, jacGdes)

        case ('vesselcoordinates_goat')

            ! Illegal, throw error
            call gdErrorHandler('EvaluteCostFunctionGR: goat variables' // & 
                ' cannot be present as explicit design variables, ' // & 
                ' not supported')

        case default 

            ! Unknown 
            call gdErrorHandler('EvaluateCostFunctionGR: design variable' // & 
                ' gradient w.r.t. design variables of type: ' // & 
                designvariables%type // ' are not implemented')

        end select

        ! Compute lagrange multipliers
        gradGgoat = jacGgoat%Transpose()
        call SolveSparseLinearSystemDI(gradGgoat, -gradJgoat, lambdaG, flag)

        ! Compute gradient
        gradGdes = jacGdes%Transpose()
        gradJ = gradJR + gradGdes%MatrixVectorProduct(lambdaG)

        ! Compute hessian
        !================
        ! Only take the raw cost function contribution as hessian - 
        ! perhaps replace by hessian estimator in the future? 
        hessJ = hessJR

        ! Housekeeping
        !=============
        ! Optional arguments
        if (present(dJdvarin)) then 
            dJdvarin = dJdvar 
        end if 
        if (present(dgradJdvarin)) then 
            dgradJdvarin = dgradJdvar
        end if




    end subroutine

end module