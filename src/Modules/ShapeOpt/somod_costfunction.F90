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
            hessJ, goat, dogradient, dohessian, designvariables)

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

        end subroutine

    end interface

    contains

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

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
        goat, dogradient, dohessian, designvariables)

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

        ! Housekeeping
        !=============
        end associate

    end subroutine

end module