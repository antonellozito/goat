!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains the cost function implementation for shape
! optimization with SOLPS. It assumes that (at least) the b25 library
! with cost function (gradient) evaluation routines is available and
! linked to. It also makes use of the types defined in b2mod. 


module sosmod_costfunction
    
    ! Initialize
    !============
    ! Load modules
    use b2mod_switches_diff
    use b2us_map_diff
    use b2mod_ad_diff, only : nncf 
    use b2mod_main_diff, only : b2mn_init_b
    use b2us_geo_diff
    use b2us_plasma_diff
    use b2mod_agdr_diff, only : b2agdr_opt, b2agdr_opt_b, b2agdr_init, &
        b2agdr_init_b, b2agdr_write, b2agdr_fin, b2agdr_fin_b
    use b2mod_par_opt_diff
    use somod_costfunction 
    use gdmod_optimizationengine

    ! Rename precision ... 
    use b2mod_types, only: R8_B25 => R8
    use mod_precision, only: R8_G => R8

    

    

    ! The usual
    implicit none
    save
    public 


    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    ! Derived types
    !==============
    ! Basic solps cost function (not inherited, so has to be included
    ! in another cost function)
    type :: CostFunctionSOLPSUDT

        ! Description
        !============
        ! This cost function is a reduced cost function that is 
        ! fully based on a cost function definition from SOLPS. It is 
        ! assumed that there is a routine called 
        ! 'EvaluateCostFunctionGradient' that returns the cost 
        ! function value and its gradient. Our routines wrap around it
        ! and this UDT defines the required data structures from SOLPS
        ! to initialize, evaluate and update the cost function. 

        ! Note: we follow the same naming convention as defned in the
        ! SOLPS routines for (hopefully) clarity

        ! Fields
        type(geometry)          :: g ! geometry and grid
        type(geometry_diff)     :: gb ! geometry derivatives 
        type(b2state)           :: st ! state
        type(b2state_diff)      :: stb ! state derivatives
        type(switches)          :: switches 
        type(switches_diff)     :: switches_diff
        type(b2stateext)        :: state_ext 
        type(b2stateext_diff)   :: state_extb
        type(b2average)         :: state_avg 
        type(b2average)         :: state_avgb
        type(switches)          :: switch 
        type(switches_diff)     :: switchb 
        type(mapping)           :: mpg 
        type(mapping_diff)      :: mpgb

    contains 

        ! Initialization
        procedure :: Initialize             => InitializeCostfunctionSOLPS

        ! Evaluation
        procedure :: Evaluate               => EvaluateCostFunctionSOLPS
        
    end type

    ! Goat and SOLPS reduced cost function
    type, extends(CostFunctionGRUDT)    :: CostFunctionGSRUDT

        ! Description
        !============
        ! This cost function type inherits from the goat reduced
        ! optimization cost function and should be used to incorporate
        ! the SOLPS cost function together with the goat. Additional 
        ! cost function contributions from the purely goat side may be
        ! added in the future. Here, we need to overwrite the GR 
        ! initialize and evaluation routines.

        ! Additional fields
        type(CostFunctionSOLPSUDT)      :: cfvsolps 
        
    contains 

        ! Initialization
        procedure :: Initialize         => InitializeCostFunctionGSR

        ! Evaluation
        procedure :: Evaluate           => EvaluateCostFunctionGSR


    end type


    contains

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                        SOLPS COST FUNCTION                       !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCostFunctionSOLPS(costfunction)

        ! Description
        !============
        ! Initialize the SOLPS cost function. This is simply a wrapper
        ! for the InitializeCostFunction.F90 routine declared on the
        ! SOLPS side. 

        ! Declare variables
        !==================
        ! Arguments
        class(CostFunctionSOLPSUDT)         :: costfunction 

        ! Call initializer
        !=================
        call InitializeCostFunction(costfunction%switch, &
            costfunction%switchb, costfunction%g, costfunction%gb, &
            costfunction%mpg, costfunction%mpgb, costfunction%st, &
            costfunction%stb, costfunction%state_ext, costfunction%state_extb, &
            costfunction%state_avg, costfunction%state_avgb)
    
    end subroutine 

     ! Evaluation
    subroutine EvaluateCostFunctionSOLPS(costfunction, J, gradJ)

        ! Description
        !============
        ! Evaluate the SOLPS cost function and its gradient (hessian is
        ! not returned). The gradient is assumed to be with respect to
        ! the grid coordinates (vxx, vxy), the psi values at the 
        ! vertices (vxfpsi), the derivatives of the psi values at the
        ! vertices (vxbx, vxby) and the ffbz value. This routine is a 
        ! wrapper for the EvaluateCostFunctionGradient routine

        ! Use modules
        !============
        

        ! Declare variables
        !==================
        ! Arguments
        class(CostFunctionSOLPSUDT)                 :: costfunction 
        real(R8_G), intent(out)                  :: J 
        real(R8_G), allocatable, intent(out)     :: gradJ(:)

        ! Auxiliary
        real(R8_B25)                        :: J1(nncf), J1b(nncf)

        ! Initialize
        !===========
        ! Associate
        associate(&
            switch          => costfunction%switch,         &
            switchb         => costfunction%switchb,        &
            g               => costfunction%g,              &
            gb              => costfunction%gb,             &
            mpg             => costfunction%mpg,            &
            mpgb            => costfunction%mpgb,           &
            st              => costfunction%st,             &
            state_ext       => costfunction%state_ext,      &
            state_extb      => costfunction%state_extb,     &
            state_avg       => costfunction%state_avg,      &
            state_avgb      => costfunction%state_avgb,     &
            stb             => costfunction%stb             &
            )

        ! Evaluate
        !=========
        ! Call routine
        call EvaluateCostFunctionGradient(switch, switchb, g, gb, &
            mpg, mpgb, st, stb, state_ext, state_extb, state_avg, &
            state_avgb, J1, J1b)

        ! Extract and cast into our precision format
        J = real(J1(0), kind=R8_G) ! assumed first entry is total cost function
        gradJ = real([gb%vxx, gb%vxy, gb%vxfpsi, gb%vxbx, gb%vxby, &
            gb%vxffbz], kind=R8_G)
        
        ! Housekeeping
        !=============
        end associate

    end subroutine

    !------------------------------------------------------------------!
    !                GOAT AND SOLPS REDUCED COST FUNCTION              !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCostFunctionGSR(costfunction, goat, options)

        ! Description
        !============
        ! Initialize the goat reduced cost function. This basically means
        ! initializing the goat engine and the costfunction type. The 
        ! cost function is then allocated to the correct type and 
        ! initialized here as well. 

        ! Arguments
        class(CostfunctionGSRUDT)           :: costfunction
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

        ! Initialize GOAT-based costfunction
        !===================================
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

        ! Initialize SOLPS cost function
        !===============================
        ! Simply call initializer
        call costfunction%cfvsolps%Initialize()

    end subroutine

    ! Evaluation
    subroutine EvaluateCostFunctionGSR(costfunction, J, gradJ, hessJ, &
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
        class(CostfunctionGSRUDT)       :: costfunction 
        real(R8_G)                      :: J
        real(R8_G), allocatable         :: gradJ(:) ! assumed initialized
        type(MySparseUDT)               :: hessJ ! assumed in initialized
        type(OptimizationProblemGDUDT)  :: goat
        logical                         :: dogradient, dohessian 
        class(DesignVariablesSOUDT)     :: designvariables

        ! Optional arguments
        character(*), intent(in), optional  :: varin 
        real(R8_G), intent(in), optional    :: valuesin(:)
        real(R8_G), allocatable, optional   :: dJdvarin(:) 
        type(MySparseUDT), optional         :: dgradJdvarin

        character(:), allocatable           :: var
        real(R8_G), allocatable             :: values(:)
        real(R8_G), allocatable             :: dJdvar(:) 
        type(MySparseUDT)                   :: dgradJdvar

        ! Auxiliary
        integer(I8)                         :: flag

        real(R8_G)                                  :: Js, Jg
        real(R8_G), allocatable, dimension(:)       :: goatvariables, &
            gradJR, gradJgoat, lambdaG, gradJs, gradJg, dpsidx, dpsidy, &
            d2psidx2, d2psidxdy, d2psidy2, gradJsolps

        type(MySparseUDT)                           :: hessJg, &
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
        ! GOAT side
        call costfunction%costfunction%Evaluate(Jg, gradJg, hessJg, &
            goat, dogradient, dohessian, designvariables, &
            'goatvariables', goatvariables, gradJgoat)

        ! SOLPS side (gradient is w.r.t. coordinates, psi, dpsidx, dpsidy, ffbz)
        call costfunction%cfvsolps%Evaluate(Js, gradJs)

        ! Compute goat linearization 
        !===========================
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

        ! Compute gradients
        !==================
        ! Associate for ease
        associate(nv            => goat%grid%vert%ntot, &
            x                   => goat%grid%vert%x,    &
            y                   => goat%grid%vert%y,    &
            MFinterp            => goat%magneticField%interp)

        ! Compute partial derivatives of psi and dpsidx, dpsidy w.r.t. 
        ! grid coordinates (simply 'diagonal' matrix linearization)
        allocate(dpsidx(nv), dpsidy(nv), d2psidx2(nv), d2psidy2(nv), &
            d2psidxdy(nv))
        call MFinterp%Evaluate(x, y, 1, 0, dpsidx)
        call MFinterp%Evaluate(x, y, 0, 1, dpsidy)
        call MFinterp%Evaluate(x, y, 2, 0, d2psidx2)
        call MFinterp%Evaluate(x, y, 1, 1, d2psidxdy)
        call MFinterp%Evaluate(x, y, 0, 2, d2psidy2)

        ! Compute solps gradient w.r.t. goat variables (assumed no contribution of ffbz)
        gradJsolps = gradJgoat ! initialize
        gradJsolps = 0
        gradJsolps(1:2*nv) = gradJs(1:2*nv) ! dJdx, dJdy 
        gradJsolps(1:nv) = gradJsolps(1:nv) + gradJs(2*nv+1:3*nv)*dpsidx &
            + gradJs(3*nv+1:4*nv)*d2psidx2 + gradJs(4*nv+1:5*nv)*d2psidxdy ! dJdpsi*dpsidx + dJd(dpsidx)*d(dpsidx)dx + dJd(dpsidy)*d(dpsidy)dx 
        gradJsolps(nv+1:2*nv) = gradJsolps(nv+1:2*nv) + gradJs(2*nv+1:3*nv)*dpsidy &
            + gradJs(3*nv+1:4*nv)*d2psidxdy + gradJs(4*nv+1:5*nv)*d2psidy2 ! dJdpsi*dpsidy + dJd(dpsidx)*d(dpsidx)dy + dJd(dpsidy)*d(dpsidy)dy 

        ! Compute lagrange multipliers
        gradGgoat = jacGgoat%Transpose()
        call SolveSparseLinearSystemDI(gradGgoat, -(gradJgoat + gradJsolps), lambdaG, flag)

        ! Compute gradient
        gradGdes = jacGdes%Transpose()
        gradJ = gradJR + gradGdes%MatrixVectorProduct(lambdaG)

        ! End association
        end associate

        ! Compute hessian
        !================
        ! Only take the raw cost function contribution as hessian - 
        ! perhaps replace by hessian estimator in the future? 
        hessJ = hessJg

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