!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains the constraint classes specific for the grid 
! deformation. It relies on user input defined in the gdmod_userinput
! module, and on the gdmod_user types. 

! The constraints are structured as follows:
! - All derived constraint types inherit from the 'mother' type 
!   'GenericConstraintsGDUDT', which contains the field 'ncon', the 
!   number of constraints, and the initialization and evaluation 
!   routines. 
! - The overarching constraints structure contains an object 'eqcon' and 
!   'ineqcon', which are objects that contain the specific equality and
!   inequality constraints, respectively. Both objects have general
!   initialization, evaluation, and destruction routines that should be
!   used in the optimizer. 
! - The equality constraints contain different fields (e.g. fluxfunction
!   ) with a logical (e.g. dofluxfunction) that indicates whether the
!   constraint should be considered. 
! - Each specific constraints (e.g. fluxfunction) has its own evaluation
!   , initialization, and destruction routines (type-bound). 

module gdmod_constraints
    
    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_sparseinterface
    use optmod_constraints
    use gdmod_types
    use gdmod_designvariables
    

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
    ! Generic constraint type
    type, abstract, extends(ConstraintsUDT) :: GenericConstraintsGDUDT

        ! Description
        !============
        ! Generic type for grid deformation constraints. Inherits from 
        ! the mother constraint type defined in gdmod_constraints.

        ! Fields:
        integer(I8)                 :: ncon = 0 ! number of constraints

    contains 

        ! Initialization
        procedure(InitializeConstraintsINT), deferred :: Initialize 

        ! Evaluation
        procedure(EvaluateConstraintsINT), deferred :: Evaluate

    end type

    ! Specific constraint types
    !==========================
    ! Flux function constraints
    type, extends(GenericConstraintsGDUDT) :: FluxfunctionConstraintsUDT

        ! Description
        !============
        ! Flux function constraints. Fixes the flux function values 
        ! for a set of desired nodes. The following fields are added:
        ! - vert:       the vertices to consider (1 constraint per 
        !               entry of this array, so ncon-by-1 dimension)
        ! - PsiD:       Desired psi value of each vertex (ncon-by-1)
        ! - ncon:       (inherited) number of constraints 

        ! No other routines than the standard initialization, evaluation
        ! and destruction routines are implemented nor needed. 

        ! Fields: 
        integer(I8), allocatable        :: vert(:)
        real(R8), allocatable           :: PsiD(:)

    contains

        ! Initialization
        procedure :: Initialize     => InitializeFluxfunctionConstraints

        ! Evaluation
        procedure :: Evaluate       => EvaluateFluxfunctionConstraints

        ! Destructor
        final :: DestroyFluxfunctionConstraints

    end type

    ! Boundary function constraints
    type, extends(ConstraintsUDT) :: BoundaryFunctionConstraintsUDT

    end type

    ! Overarching types
    !==================
    ! Equality constraints
    type, extends(ConstraintsUDT) :: EqConGDUDT

        ! Description
        !============
        ! This type contains all the different constraints as different
        ! derived types. For each type of constraint, a different
        ! type is defined. 

        ! Total number of constraints 
        integer(I8)                             :: neqcon = 0

        ! Flux function constraint
        logical                                 :: dofluxfunction = .false.
        type(FluxfunctionConstraintsUDT)        :: fluxfunction 
        

    contains

        ! Procedure to initialize constraints
        procedure :: Initialize         => InitializeEqCon

        ! Procedure to evaluate constraints
        procedure :: Evaluate           => EvaluateEqCon

    end type 

    ! Inequality constraints
    type, extends(ConstraintsUDT) :: IneqConGDUDT

        ! Total number of inequality constraints
        integer(I8)             :: nineqcon = 0

    contains

        ! Constraints initialization
        procedure :: Initialize         => InitializeIneqCon

        ! Constraints evaluation
        procedure :: Evaluate           => EvaluateIneqCon

    end type

    ! All constraints for the grid deformation
    type, extends(ConstraintsUDT) :: ConstraintsGDUDT

        ! Description
        !============
        ! Defines the basic optimization problem: it has a set of 
        ! design variables, constraints, and a cost function. 

        ! Fields: 

        ! Equality constraints 
        type(EqConGDUDT)        :: eqcon

        ! Inequality constraints 
        type(IneqConGDUDT)      :: ineqcon

    contains

        ! Initialization
        procedure :: Initialize         => InitializeConstraints

        ! Number of constraints getter
        procedure :: GetConstraintsDimensions  

        ! Evaluation
        ! procedure :: Evaluate           => EvaluateEqualityConstraintsGD
        
        ! Housekeeping

    end type

    ! Monitor
    !========
    type ConstraintsMonitorUDT

        ! Description
        !============
        ! This object can be used to keep track of the amount of 
        ! constraints that are imposed per vertex. As such, it can be
        ! checked whether the problem will be overly constrained and to 
        ! prevent constraints to be imposed on certain vertices if 
        ! a maximum number is reached. For now, the equality constraints
        ! per vertex are counted in eqvcc, the inequality in ineqvcc. 
        ! The maximum number of equality and inequality constraints per 
        ! vertex are set in maxeqvcc and maxineqvcc. The methods are 
        ! the following:

        ! - Initialize: initialization routine that takes the grid, 
        !   magnetic field, and environment as input (in case of 
        !   extension of this object in the future)
        ! - CheckConstraints: will check if the problem is already 
        !   overly constrained, or may become overly constrained if 
        !   inequality constraints become active. 

        ! Fields:
        integer(I8), allocatable           :: eqvcc(:), ineqvcc(:)
        integer(I8)                        :: maxeqvcc, maxineqvcc 
        
    contains 

        ! Routines
        procedure :: Initialize         => InitializeMonitor
        ! procedure :: CheckConstraints   => CheckConstraintsMonitor

        ! Finalization
        final :: DestroyMonitor

    end type

    

    ! Derived types
    !==============
    

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    ! Abstract interfaces
    !====================
    abstract interface

        ! Constraint initialization
        subroutine InitializeConstraintsINT(constraints, grid, &
            magneticField, environment, monitor)

            ! Description
            !============
            ! This routine serves as general initialization routine for
            ! a generic grid deformation constraint (that inherits from
            ! the GenericConstraintUDT type)

            ! Import
            import :: GenericConstraintsGDUDT, GridUDT, &
                MagneticFieldUDT, EnvironmentUDT, ConstraintsMonitorUDT

            ! Declare
            class(GenericConstraintsGDUDT)      :: constraints 
            type(GridUDT)                       :: grid 
            type(MagneticFieldUDT)              :: magneticField 
            type(EnvironmentUDT)                :: environment 
            type(ConstraintsMonitorUDT)         :: monitor

        end subroutine

        ! Constraint evaluation
        subroutine EvaluateConstraintsINT(constraints, G, gradG, & 
            hessG, grid, magneticField, environment, &
            dogradient, dohessian, designvariables, lambda)

            ! Description
            !============
            ! This reoutine serves as a general evaluation routine for 
            ! a generic grid deformation constraint. 

            ! Import
            import :: GenericConstraintsGDUDT, MySparseUDT, GridUDT, &
                R8, MagneticFieldUDT, EnvironmentUDT, & 
                DesignVariablesGDUDT
            
            ! Declare
            class(GenericConstraintsGDUDT)  :: constraints 
            real(R8), allocatable           :: G(:), lambda(:)
            type(MySparseUDT)               :: hessG, gradG 
            type(GridUDT)                   :: grid 
            type(MagneticFieldUDT)          :: magneticField 
            type(EnvironmentUDT)            :: environment 
            logical                         :: dogradient, dohessian
            class(DesignVariablesGDUDT)     :: designvariables

        end subroutine

    end interface

    contains

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                           CONSTRAINTS MONITOR                    !
    !------------------------------------------------------------------!
    ! Initialization
    subroutine InitializeMonitor(monitor, grid, magneticField, &
        environment)

        ! Description
        !============
        ! Initializes the constraints monitor structure. It is assumed
        ! that the grid, magnetic field and environment are properly 
        ! allocated and initialized. 

        ! Declare variables
        !==================
        ! Arguments
        class(ConstraintsMonitorUDT)        :: monitor 
        type(GridUDT)                       :: grid 
        type(MagneticFieldUDT)              :: magneticField 
        type(EnvironmentUDT)                :: environment

        ! Loop variables

        ! Auxiliary 

        ! Initialize
        !===========
        ! Set the maximal number of constraints
        monitor%maxeqvcc = 2 ! for most cases this is fine
        monitor%maxineqvcc = 1000 ! a stupid large number - can impose any number

        ! Allocate
        allocate(monitor%eqvcc(grid%vert%ntot))
        allocate(monitor%ineqvcc(grid%vert%ntot))

        ! Initialize
        monitor%eqvcc(:)    = 0
        monitor%ineqvcc(:)  = 0

    end subroutine

    ! Destructor
    subroutine DestroyMonitor(monitor)

        ! Description
        !============
        ! Finalization for the constraints monitor object

        ! Declare variables
        !==================
        ! Arguments
        type(ConstraintsMonitorUDT)         :: monitor 

        ! Destroy
        !========
        deallocate(monitor%eqvcc, monitor%ineqvcc)

    end subroutine

    !------------------------------------------------------------------!
    !                           GENERAL CONSTRAINTS                    !
    !------------------------------------------------------------------!
    ! Initialization
    subroutine InitializeConstraints(constraints, grid, magneticField, &
        environment, constraintoptions)

        ! Description
        !============
        ! Routine that initializes the equality and inequality 
        ! constraints, using the initialization routines of those 
        ! objects. 

        ! Declare variables
        !==================
        ! Arguments
        class(ConstraintsGDUDT)     :: constraints 
        type(GridUDT)               :: grid 
        type(MagneticFieldUDT)      :: magneticField 
        type(EnvironmentUDT)        :: environment 
        type(ConstraintOptionsUDT)  :: constraintoptions
        type(ConstraintsMonitorUDT) :: monitor

        ! Loop variables

        ! Auxiliary variables

        ! Initialize monitor
        !===================
        call monitor%Initialize(grid, magneticField, environment)

        ! Initialize constraints
        !=======================
        ! Equality constraints
        call constraints%eqcon%Initialize(grid, magneticField, &
            environment, constraintoptions, monitor)

        ! Inequality constraints
        call constraints%ineqcon%Initialize(grid, magneticField, &
            environment, constraintoptions, monitor)

    end subroutine

    ! Dimension getter
    subroutine GetConstraintsDimensions(constraints, neqcon, nineqcon)

        ! Description
        !============
        ! Return the current dimensions of the equality and inequality
        ! constraints. Can be used for initialization of other 
        ! quantities at higher levels. 

        ! Declare variables
        !==================
        ! Arguments
        class(ConstraintsGDUDT)         :: constraints 
        integer(I8), intent(out)        :: neqcon, nineqcon

        ! Extract dimensions
        !===================
        ! Stored in eqcon, ineqcon
        neqcon = constraints%eqcon%neqcon 
        nineqcon = constraints%ineqcon%nineqcon

    end subroutine

    

    !------------------------------------------------------------------!
    !                           EQUALITY CONSTRAINTS                   !
    !------------------------------------------------------------------!
    ! Initialization
    subroutine InitializeEqCon(constraints, grid, magneticField, &
        environment, constraintoptions, monitor)

        ! Description
        !============
        ! Routine that initializes the desired constraints, based on the
        ! user options d

        ! Declare variables
        !==================
        ! Arguments
        class(EqConGDUDT)           :: constraints 
        type(GridUDT)               :: grid 
        type(MagneticFieldUDT)      :: magneticField 
        type(EnvironmentUDT)        :: environment 
        type(ConstraintOptionsUDT)  :: constraintoptions
        type(ConstraintsMonitorUDT) :: monitor

        ! Loop variables

        ! Auxiliary variables
        
        ! Initialize constraints
        !=======================
        constraints%neqcon = 0

        ! Flux function
        if (constraintoptions%fluxfunction) then 
            ! Set the logical
            constraints%dofluxfunction = .true.

            ! Initialize
            call constraints%fluxfunction%Initialize(grid, &
                magneticField, environment, monitor)

            ! Add constraints number
            constraints%neqcon = constraints%neqcon + &
                constraints%fluxfunction%ncon 

        else
            ! Set to false, don't initialize
            constraints%dofluxfunction = .false.

        end if

    end subroutine

    ! Constraint evaluation
    subroutine EvaluateEqCon(constraints, G, gradG, hessG, &
        grid, magneticField, environment, dogradient, dohessian, & 
        designvariables, lambda)

        ! Description
        !============
        ! This routine evaluates the constraints G and the corresponding
        ! gradient and hessian. To do so, every type of constraint is 
        ! checked whether it is imposed, and the contributions are 
        ! added by calling the evaluation routine of each constraint. 

        ! Declare variables
        !==================
        ! Arguments
        class(EqConGDUDT)               :: constraints
        real(R8), intent(inout)         :: G(:)
        real(R8), intent(in)            :: lambda(:)
        type(MySparseUDT)               :: gradG, hessG 
        type(GridUDT)                   :: grid
        type(MagneticFieldUDT)          :: magneticField 
        type(EnvironmentUDT)            :: environment
        logical                         :: dogradient, dohessian 
        class(DesignVariablesGDUDT)     :: designvariables 

        ! Loop
        integer(I8)                     :: ic, ivg, ivh, k
        integer(I8), allocatable        :: conindex(:)

        ! Auxiliary
        real(R8), allocatable           :: G_flux(:), lambda_flux(:)
        type(MySparseUDT)               :: gradG_flux, hessG_flux

        ! Initialize
        !===========
        ! Set the constraint counter
        ic = 0

        ! Flux function constraints
        !==========================
        if (constraints%dofluxfunction) then 
            ! Construct the constraint index
            allocate(conindex(constraints%fluxfunction%ncon))
            conindex = [(k, k = ic+1, ic+constraints%fluxfunction%ncon)]

            ! Allocate & initialize
            allocate(lambda_flux(constraints%fluxfunction%ncon))
            lambda_flux = lambda(conindex)

            ! Call the evaluation routine
            call constraints%fluxfunction%Evaluate(G_flux, &
                gradG_flux, hessG_flux, &
                grid, magneticField, environment, dogradient, &
                dohessian, designvariables, &
                lambda_flux)

            ! Assign
            G(conindex) = G_flux

            ! Update the gradient column indices
            if (dogradient) then 
                ! For easier concatenation later on
                gradG_flux%col = gradG_flux%col + ic

            end if

            ! Update the constraint counter
            ic = ic + constraints%fluxfunction%ncon

            ! Housekeeping
            deallocate(conindex, lambda_flux) 
            if (allocated(G_flux)) then
                deallocate(G_flux)
            end if


        end if

        ! Concatenate gradient
        !=====================
        if (dogradient) then 

            ! Determine sizes
            !----------------
            ! Size of the gradient
            gradG%ncol = constraints%neqcon 
            gradG%nrow = designvariables%nphi

            ! Allocate
            if (.not. allocated(gradG%val)) then 
                ! Number of values (to be determined)
                gradG%nval = 0

                ! Add values of each constraint, if used
                if (constraints%dofluxfunction) then 
                    gradG%nval = gradG%nval + gradG_flux%nval  
                end if 

                ! Allocate
                call gradG%Allocate()
            end if

            ! Add contributions
            !------------------
            ! Initialize counter
            ivg = 0
            
            ! Flux function
            if (constraints%dofluxfunction) then 
                ! Associate 
                associate(&
                    nc      => constraints%fluxfunction%ncon, &
                    nval    => gradG_flux%nval)

                ! Add values
                gradG%row(ivg+1:ivg+nval) = gradG_flux%row 
                gradG%col(ivg+1:ivg+nval) = gradG_flux%col
                gradG%val(ivg+1:ivg+nval) = gradG_flux%val

                ! Update counter
                ivg = ivg + nval 

                ! End associate
                end associate

            end if

        end if

        ! Concatenate the hessian
        !========================
        if (dohessian) then 

            ! Determine sizes
            !----------------
            ! Size of the hessian
            hessG%ncol = designvariables%nphi
            hessG%nrow = designvariables%nphi

            ! Allocate
            if (.not. allocated(hessG%val)) then 
                ! Number of values (to be determined)
                hessG%nval = 0

                ! Add values of each constraint, if used
                if (constraints%dofluxfunction) then 
                    hessG%nval = hessG%nval + hessG_flux%nval  
                end if 

                ! Allocate
                call hessG%Allocate()
            end if

            ! Add contributions
            !------------------
            ! Initialize counter
            ivh = 0
            
            ! Flux function
            if (constraints%dofluxfunction) then 
                ! Associate 
                associate(&
                    nc      => constraints%fluxfunction%ncon, &
                    nval    => hessG_flux%nval)

                ! Add values
                hessG%row(ivh+1:ivh+nval) = hessG_flux%row 
                hessG%col(ivh+1:ivh+nval) = hessG_flux%col
                hessG%val(ivh+1:ivh+nval) = hessG_flux%val

                ! Update counter
                ivh = ivh + nval 

                ! End associate
                end associate

            end if

        end if


    end subroutine


    !------------------------------------------------------------------!
    !                          INEQUALITY CONSTRAINTS                  !
    !------------------------------------------------------------------!
    ! Initialization
    subroutine InitializeIneqCon(constraints, grid, magneticField, &
        environment, constraintoptions, monitor)

        ! Description
        !============
        ! Routine that initializes the desired constraints, based on the
        ! user options d

        ! Declare variables
        !==================
        ! Arguments
        class(IneqConGDUDT)         :: constraints 
        type(GridUDT)               :: grid 
        type(MagneticFieldUDT)      :: magneticField 
        type(EnvironmentUDT)        :: environment 
        type(ConstraintOptionsUDT)  :: constraintoptions
        type(ConstraintsMonitorUDT) :: monitor

        ! Loop variables

        ! Auxiliary variables

        ! Initialize
        !===========

        ! Initialize constraints
        !=======================
        constraints%nineqcon = 0

    end subroutine

    ! Constraint evaluation
    subroutine EvaluateIneqCon(constraints, G, gradG, hessG, &
        grid, magneticField, environment, dogradient, dohessian, & 
        designvariables, lambda)

        ! Description
        !============
        ! This routine evaluates the constraints G and the corresponding
        ! gradient and hessian. To do so, every type of constraint is 
        ! checked whether it is imposed, and the contributions are 
        ! added by calling the evaluation routine of each constraint. 

        ! Declare variables
        !==================
        ! Arguments
        class(IneqConGDUDT)             :: constraints
        real(R8), intent(inout)         :: G(:)
        real(R8), intent(in)            :: lambda(:)
        type(MySparseUDT)               :: gradG, hessG 
        type(GridUDT)                   :: grid
        type(MagneticFieldUDT)          :: magneticField 
        type(EnvironmentUDT)            :: environment
        logical                         :: dogradient, dohessian 
        class(DesignVariablesGDUDT)     :: designvariables 

        ! Loop
        integer(I8)                     :: ic, ivg, ivh, k
        integer(I8), allocatable        :: conindex(:)

        ! Auxiliary

        ! Initialize
        !===========
        ! Set the constraint counter
        ic = 0

        ! Initialize
        G(:) = 0
    
        ! Concatenate gradient
        !=====================
        if (dogradient) then 

            ! Determine sizes
            !----------------
            ! Size of the gradient
            gradG%ncol = constraints%nineqcon 
            gradG%nrow = designvariables%nphi

            ! Allocate
            if (.not. allocated(gradG%val)) then 
                ! Number of values (to be determined)
                gradG%nval = 0

                ! Add values of each constraint, if used

                ! Allocate
                call gradG%Allocate()
            end if

            ! Add contributions
            !------------------
            ! Initialize counter
            ivg = 0
            
        end if

        ! Concatenate the hessian
        !========================
        if (dohessian) then 

            ! Determine sizes
            !----------------
            ! Size of the hessian
            hessG%ncol = designvariables%nphi
            hessG%nrow = designvariables%nphi

            ! Allocate
            if (.not. allocated(hessG%val)) then 
                ! Number of values (to be determined)
                hessG%nval = 0

                ! Add values of each constraint, if used

                ! Allocate
                call hessG%Allocate()
            end if

            ! Add contributions
            !------------------
            ! Initialize counter
            ivh = 0

        end if


    end subroutine

    !------------------------------------------------------------------!
    !                           FLUX FUNCTION                          !
    !------------------------------------------------------------------!
    ! Initialize
    subroutine InitializeFluxfunctionConstraints(constraints, grid, &
        magneticField, environment, monitor)

        ! Description
        !============
        ! Initialize the required fields of the flux function 
        ! constraints. These constraints impose for each vertex that 
        ! lies on a flux surface (i.e. has a flux surface ID which is
        ! non-zero) by fixing its flux value. The flux values can be
        ! determined in different ways. Here, we compute the initial
        ! psi values by averaging the current psi values of the vertices
        ! that lie on a flux surface. This hedges a bit for 
        ! discretization errors originating from the re-interpretation 
        ! of the magnetic field data as a bicubic spline interpolant 
        ! (here) instead of a linear interpolant in the grid generator. 

        ! Notes
        !======
        ! Note 1: the last argument of this function is a derived type 
        ! used to monitor the constraints. It can be used to make sure
        ! that the problem is - at least not by the equality constraints 
        ! - is overly constrained. It is UP TO THE DEVELOPER to use this
        ! monitor properly!

        ! Note 2: the constraint options are passed to this function in 
        ! order to determine whether e.g. boundary nodes should be 
        ! considered for the constraints. 

        ! Note 3: for vessel mode grids (grids reaching up to the 
        ! vessel wall), some vertices may not belong to a flux surface.
        ! To fix them anyway on the vessel wall, the logical 
        ! 'fixfarvesselflux' can be set to 'true'. This is the 
        ! recommended default value. 

        ! Note 4: for target mode grids, the flux value at the 'corners'
        ! is typically also determined by averaging over the field line.
        ! However, this gives issues in some cases where the targets are
        ! nearly flux-aligned. To avoid this, set the logical 
        ! 'fixfluxalignedtargets' to true. This is the recommended 
        ! default value. 

        ! Note 5: it is assumed that the grid contains the flux surface
        ! data and that the flux surfaces are numbered from 1 to nFs

        ! Initialize
        !===========
        ! Modules
        use BicubicSplineInterpolant
        use gdmod_plots
        
        ! Declare variables
        !==================
        ! Arguments 
        class(FluxfunctionConstraintsUDT)       :: constraints 
        type(GridUDT)                           :: grid 
        type(MagneticFieldUDT)                  :: magneticField 
        type(EnvironmentUDT)                    :: environment 
        type(ConstraintsMonitorUDT)             :: monitor

        ! Loop variables
        integer(I8)                 :: i

        ! Auxiliary variables
        integer(I8), allocatable    :: vert_tmp(:), vertID(:) 
        real(R8), allocatable       :: PsiD_tmp(:) 
        logical, allocatable        :: delind(:), mask(:)
        real(R8)                    :: tpsi

        ! Data
        logical                     :: fixfluxalignedtargets = .true. 
        logical                     :: fixfarvesselflux = .true. 

        ! Initialize
        !===========
        ! Number of constraints
        constraints%ncon = 0

        ! Allocate temporary arrays
        allocate(vert_tmp(grid%vert%ntot))
        allocate(PsiD_tmp(grid%vert%ntot))

        ! Allocate auxiliary arrays
        allocate(mask(grid%vert%ntot))
        allocate(delind(grid%vert%ntot))
        allocate(vertID(grid%vert%ntot))

        ! Initialize
        PsiD_tmp(:) = 0
        mask(:) = .false.
        delind(:) = .false.
        vert_tmp = [(i, i = 1, grid%vert%ntot)]

        ! Associate
        associate(&
            vert    => grid%vert,       x       => grid%vert%x,     &
            y       => grid%vert%y,     cc      => monitor%eqvcc,   &
            maxcc   => monitor%maxeqvcc)
        
        ! Determine flux values to impose
        !================================
        ! Evaluate flux at all nodes
        call EvaluateBicubicSplineInterpolant(x, y, PsiD_tmp, &
            magneticField%interp, '0', '0')

        ! call Plot2DUnstructuredField(PsiD_tmp, grid, 'v', '-p')

        ! Loop over all the flux surfaces to compute desired flux
        do i = 1, grid%data%fluxdata%nFs
            ! Get all vertices with this ID
            mask = (vert%fieldlineID == i)

            ! Get the flux values
            if (any(pack(vert%BV,mask))) then 
                ! Average only over boundary vertices
                tpsi = sum(pack(PsiD_tmp, (mask .and. vert%BV))) & 
                    /count((mask .and. vert%BV))
            else
                ! Average over all vertices
                tpsi = sum(pack(PsiD_tmp, mask))/count(mask)
            end if

            ! Adjust PsiD
            where(mask) PsiD_tmp = tpsi 

        end do

        ! Compensate for flux aligned targets?
        if (fixfluxalignedtargets) then
            ! This is still to do
            print *, 'fix for flux aligned targets is not yet available'
        end if

        ! Fix flux values of vessel boundaries? 
        if (fixfarvesselflux) then 
            ! This is still to do
            print *, 'fix for vessel boundaries is not yet available'
        end if

        ! Set the deletion vector
        where(cc >= maxcc) delind = .true. ! don't constrain

        ! Update monitor
        !===============
        where (.not. delind) cc = cc + 1

        ! Allocate and assign
        !====================
        ! Allocate
        constraints%ncon = count( .not. delind)
        allocate(constraints%vert(constraints%ncon))
        allocate(constraints%PsiD(constraints%ncon))

        ! Assign
        !=======
        constraints%vert = pack(vert_tmp, (.not. delind))
        constraints%PsiD = pack(PsiD_tmp, (.not. delind))

        ! Housekeeping
        !=============
        ! End associate
        end associate

        ! Deallocate temporary arrays
        deallocate(vert_tmp)
        deallocate(PsiD_tmp)

        ! Deallocate auxiliary arrays
        deallocate(mask)
        deallocate(delind)
        deallocate(vertID)
        
    end subroutine

    ! Evaluation
    subroutine EvaluateFluxfunctionConstraints(constraints, G, gradG, & 
        hessG, grid, magneticField, environment, dogradient, &
        dohessian, designvariables, lambda)

        ! Description
        !============
        ! Evaluate the flux function constraints imposed on the 
        ! vertices. For each  vertex considered (see InitDesign), the 
        ! flux function is imposed mathematically as:
        ! 
        !       G_i = Psi(x_i,y_i) - Psi_D,
        !
        ! where x_i and y_i are the i-th vertex's coordinates, Psi is 
        ! the underlying flux function, characterized by a bicubic 
        ! spline interpolant, and Psi_D is a vector containing the
        ! desired flux function values. 

        ! The Hessian of each ith constraint is:
        !
        !       Hjk,i = d^2F/dx^2 if j == k == i
        !       Hjk,i = d^2F/dy^2 if j == k == i+numel(x)
        !       Hjk,i = d^2F/dxdy if j == k+numel(x)
        !       Otherwise zero
        !
        ! Therefore, the multiplication Hjk,i lambda_i is equal to:
        !
        !       Hjk,i lambda_i = d^2F/dx^2 lambda_i (idem for y)

        ! Notes
        !======
        ! Note 1: not the true hessian of the constraint vector is 
        ! returned, but the hessian-vector multiplication with the 
        ! vector lambda, which should be of suitable size. 

        ! Note 2: the row and column indices for the constraints are 
        ! local, meaning that in no way other constraints are accounted
        ! for in positioning the elements in the matrix. This should be
        ! done in an overarching routine. 

        ! Note 3: at first, we compute the linearization of the 
        ! constraints, meaning that we actually compute the Jacobian.
        ! Afterwards, we switch the row and column indices (i.e. 
        ! transpose) to obtain the gradient. 

        ! Initialize
        !===========
        ! Modules
        use BicubicSplineInterpolant
        
        ! Declare variables
        !==================
        ! Arguments 
        class(FluxfunctionConstraintsUDT)   :: constraints 
        real(R8), allocatable               :: G(:) 
        real(R8), allocatable               :: lambda(:)
        type(MySparseUDT)                   :: hessG, gradG, jacG 
        type(GridUDT)                       :: grid 
        type(MagneticFieldUDT)              :: magneticField 
        type(EnvironmentUDT)                :: environment 
        logical                             :: dogradient, dohessian
        class(DesignVariablesGDUDT)         :: designvariables                

        ! Loop variables
        integer(I8)                         :: ic, ivg, ivh, k
        integer(I8), allocatable            :: valindex(:), conindex(:)

        ! Auxiliary variables
        real(R8), allocatable               :: psival(:), dpsidx(:), &
            dpsidy(:), valxx(:), valxy(:), valyy(:)
        integer(I8)                         :: ntv

        ! Data

        ! Initialize
        !===========
        ! Checks
        if ( (.not. allocated(lambda)) .and. dohessian) then
            ! Throw error
            call gdErrorHandler('When evaluating the hessian vector' &
                // ' multiplication, lambda must be given')
        end if

        if (size(lambda) .ne. constraints%ncon) then
            ! Lambda should have the same size as the constraints
            call gdErrorHandler('Lambda should have the same size ' &
                // 'as the constraint vector')
        end if

        ! Counters
        ic = 0 ! constraint counter (local)
        ivg = 0 ! value index for gradient
        ivh = 0 ! value index for hessian

        ! Associate
        associate(&
            nc      => constraints%ncon,        &
            psiD    => constraints%PsiD,        &
            tv      => constraints%vert,        &
            Psifun  => magneticField%interp,    & 
            x       => grid%vert%x,             & 
            y       => grid%vert%y              & 
            )

        ! Constraint value
        !=================
        ! Initialize
        ntv = size(tv)

        ! Allocate
        allocate(psival(nc))
        allocate(G(nc))
        psival(:) = 0

        ! Evaluate
        call EvaluateBicubicSplineInterpolant(x(tv), y(tv), &
            psival, Psifun, '0', '0')
        G(:) = psival - psiD

        ! Constraint gradient
        !====================
        if (dogradient) then 
            ! Initialize
            jacG%nrow = nc 
            jacG%ncol = designvariables%nphi

            ! Check design variables
            select case(designvariables%type)

            case ('coordinates')

                ! Order in jacobian: first x, then y. Has as many 
                ! non-zero elements as there are design variables. 

                ! Allocate
                jacG%nval = 2*ntv
                call jacG%Allocate() 
                allocate(dpsidx(designvariables%nphi))
                allocate(dpsidy(designvariables%nphi))
                allocate(conindex(ntv))
                allocate(valindex(ntv))

                ! Compute the derivative values
                call EvaluateBicubicSplineInterpolant(&
                    x(tv), y(tv), dpsidx, Psifun, '1', '0')
                call EvaluateBicubicSplineInterpolant(&
                    x(tv), y(tv), dpsidy, Psifun, '0', '1')

                ! x-contribution
                !---------------
                ! Build indices
                conindex = [(k, k = ic+1, ic+ntv)]
                valindex = [(k, k = ivg+1, ivg+ntv)]

                ! Add values
                jacG%row(valindex) = conindex  
                jacG%col(valindex) = tv
                jacG%val(valindex) = dpsidx 

                ! y-contribution
                !---------------
                ! Build indices
                ivg = ivg + ntv
                valindex = valindex + ntv

                ! Add values
                jacG%row(valindex) = conindex 
                jacG%col(valindex) = tv + grid%vert%ntot 
                jacG%val(valindex) = dpsidy 

                ! Build gradient
                gradG%nrow = jacG%ncol 
                gradG%ncol = jacG%nrow 
                gradG%nval = jacG%nval 
                
                call gradG%Allocate()
                gradG%row = jacG%col 
                gradG%col = jacG%row
                gradG%val = jacG%val

                ! Housekeeping
                call jacG%Deallocate()

            case default

                ! Unknown, throw error
                call gdErrorHandler('Gradient not implemented for ' &
                    // 'this type of design variable')

            end select

        end if

        ! Constraint hessian
        !===================
        if (dohessian) then 

            ! Initialize
            hessG%nrow = designvariables%nphi 
            hessG%ncol = designvariables%nphi 

            ! Check design variables
            select case(designvariables%type)

            case ('coordinates')
            
                ! Allocate
                hessG%nval = 4*ntv
                if (.not. allocated(valindex)) then
                    allocate(valindex(ntv))
                end if
                if (.not. allocated(conindex)) then 
                    allocate(conindex(ntv))
                end if 
                if (.not. allocated(hessG%val)) then
                    call hessG%Allocate()
                end if
                allocate(valxx(ntv))
                allocate(valxy(ntv))
                allocate(valyy(ntv))

                ! Initialize
                valxx(:) = 0
                valyy(:) = 0
                valxy(:) = 0

                ! Compute contributions
                call EvaluateBicubicSplineInterpolant(&
                    x(tv), y(tv), valxx, Psifun, '2', '0')
                call EvaluateBicubicSplineInterpolant(&
                    x(tv), y(tv), valyy, Psifun, '0', '2')
                call EvaluateBicubicSplineInterpolant(&
                    x(tv), y(tv), valxy, Psifun, '1', '1') ! 

                ! xx-contribution
                !----------------
                k = 1
                ! Build indices
                valindex = [(k, k = ivh+1, ivh+ntv)] 

                ! Add values
                hessG%row(valindex) = tv 
                hessG%col(valindex) = tv 
                hessG%val(valindex) = valxx*lambda ! element-wise mult. 
                
                ! xy-contribution
                !----------------
                ! Build indices
                ivh = ivh + ntv
                valindex = valindex + ntv

                ! Add values
                hessG%row(valindex) = tv
                hessG%col(valindex) = tv + grid%vert%ntot 
                hessG%val(valindex) = valxy*lambda ! element-wise mult. 

                ! yx-contribution
                !----------------
                ! symmetric with xy
                ! Build indices
                ivh = ivh + ntv
                valindex = valindex + ntv 

                ! Add values
                hessG%row(valindex) = tv + grid%vert%ntot 
                hessG%col(valindex) = tv 
                hessG%val(valindex) = valxy*lambda ! element-wise mult. 

                ! yy-contribution
                !----------------
                ! Build indices
                ivh = ivh + ntv
                valindex = valindex + ntv

                ! Add values
                hessG%row(valindex) = tv + grid%vert%ntot 
                hessG%col(valindex) = tv + grid%vert%ntot 
                hessG%val(valindex) = valyy*lambda ! element-wise mult. 

            case default

                ! Unknown, throw error
                call gdErrorHandler('Gradient not implemented for ' &
                    // 'this type of design variable')

            end select

        end if

        ! Housekeeping
        !=============
        ! End associate
        end associate

        ! Deallocate
        deallocate(psival)

        if (dogradient) then 
            deallocate(dpsidx, dpsidy, valindex, conindex)
        end if 

        if (dohessian) then 
            if (allocated(valindex)) then 
                deallocate(valindex, conindex)
            end if 
            deallocate(valxx, valxy, valyy)
        end if

    end subroutine

    ! Destructor
    subroutine DestroyFluxfunctionConstraints(constraints)

        ! Description
        !============
        ! Destructor of the flux function constraints

        ! Declare variables
        !==================
        ! Arguments
        type(FluxfunctionConstraintsUDT)        :: constraints 

        ! Destroy
        !========
        if (allocated(constraints%vert)) then
            deallocate(constraints%vert)
        end if
        if (allocated(constraints%PsiD)) then
            deallocate(constraints%PsiD)
        end if

    end subroutine

end module