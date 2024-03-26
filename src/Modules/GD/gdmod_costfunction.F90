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
! type defined in the optmod_costfunction module. 

module gdmod_costfunction
    
    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_sparseinterface
    use gdmod_types 
    use gdmod_designvariables
    use gdmod_utility_optimization
    use optmod_costfunction
    use gdmod_plots
    use DistributionFunction

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
    type, abstract, extends(CostfunctionUDT) :: CostfunctionGDUDT

        ! Description
        !============
        ! Defines the basic cost function structure for the grid 
        ! deformation. The following general fields are added: 
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
        procedure(InitializeCostfunctionINT), deferred :: Initialize

        ! Cost function evaluation
        procedure(EvaluateCostFunctionINT), deferred :: Evaluate

    end type

    ! Derived types
    !==============
    ! Length ratio cost function
    type, extends(CostfunctionGDUDT) :: CostfunctionLRUDT 

        ! Description
        !============
        ! Cost function based on the length ratio distribution along
        ! coordinate lines. Here, the coordinate lines are lines of
        ! constant magnetic flux. It is assumed that the vertex neighbours
        ! are sorted either clockwise or counter clockwise. The 
        ! allocatable arrays are defined as nv-by-1 (or by 2 for the
        ! vertex pairs) arrays, where nv is the number of grid vertices.
        ! If certain vertices should not play a role in the cost 
        ! function, set the weight wt to zero for that vertex. 

        ! IMPORTANT: this cost function shouldn't be used, as it gives
        ! quite shitty results. It is, however, used in a better way
        ! by the CostfunctionLRUDT2 cost function, where the bias is 
        ! also considered in the other direction (see the description
        ! there)
    
        ! Cost function formula
        !======================
        !
        !   J_i = wt(i) * (d1/d2 - b)**2 
        !   J = lambda * sum(J_i), i = 1, vert%ntot

        ! Notes
        !======
        ! Note 1: it is assumed that each vertex only has one vertex 
        ! pair 

        ! Fields
        real(R8)                    :: lambda ! scaling constant
        real(R8)                    :: eta ! length scale to prevent NaN/Inf
        logical                     :: dovessel ! include vessel edges in cost function?
        real(R8), allocatable       :: b0(:) ! desired length ratio per vertex
        real(R8), allocatable       :: wt(:) ! weight factor per vertex
        integer(I8), allocatable    :: vpairs(:, :), nvpairs(:) ! vertex pairs

    contains

        ! Initialization
        procedure :: Initialize             => InitializeCostfunctionLR

        ! Evaluation
        procedure :: Evaluate               => EvaluateCostFunctionLR

        ! Data output
        procedure :: WriteData              => WriteCostFunctionDataLR

        ! Housekeeping
        procedure :: Allocate               => AllocateCostfunctionLR
        procedure :: Deallocate             => DeallocateCostFunctionLR
        final :: DestroyCostFunctionLR

    end type

    ! Length ratio cost function, radial
    type, extends(CostfunctionLRUDT) :: CostfunctionLRradUDT 

        ! Description
        !============
        ! Length ratio cost function, fully analogous to LR cost
        ! function, but different initialization routine. Therefore, 
        ! extends directly the CostfunctionLRUDT type.

        ! No additional fields needed
    contains

        ! Initialization to be overwritten
        procedure :: Initialize         => InitializeCostfunctionLRrad

        ! Data writing to be overwritten
        procedure :: WriteData          => WriteCostFunctionDataLRrad

    end type

    ! Length ratio 2 cost function
    type, extends(CostfunctionGDUDT) :: CostfunctionLRUDT2 

        ! Description
        !============
        ! Cost function based on the length ratio distribution along
        ! coordinate lines. Similar to the lengthratio cost function, 
        ! but now we also consider the 'inverse' of this cost function 
        ! (this is not exactly the inverse though, we simply flip the 
        ! edge order and of course the bias). Reason for this  is to not
        ! put excessive weight on smaller edges d2 (which happens in the 
        ! other cost function). This cost function is recommended over 
        ! the original length ratio cost function. 

        ! To evaluate it, we simply add the original cost function as an
        ! object and, in evaluation, we adjust the cost function 
        ! slightly. 
    
        ! Cost function formula
        !======================
        !
        !   J_i = wt(i) * ( (d1/d2 - b)**2  + (d2/d1 - 1/b)**2)
        !   J = lambda * sum(J_i), i = 1, vert%ntot

        ! Notes
        !======
        ! Note 1: it is assumed that each vertex only has one vertex 
        ! pair 

        ! Fields
        type(CostfunctionLRUDT)     :: cfv_lr

    contains

        ! Initialization
        procedure :: Initialize             => InitializeCostfunctionLR2

        ! Evaluation
        procedure :: Evaluate               => EvaluateCostFunctionLR2

        ! Housekeeping
        procedure :: Allocate               => AllocateCostFunctionLR2
        procedure :: Deallocate             => DeallocateCostFunctionLR2
        final :: DestroyCostFunctionLR2

    end type

    ! Length ratio, radial 2 cost function
    type, extends(CostfunctionGDUDT) :: CostfunctionLRrad2UDT 

        ! Description
        !============
        ! Cost function based on the length ratio in radial direction.
        ! Similar to LR2, but now for radial length ratio. 

        ! Fields
        type(CostfunctionLRradUDT)     :: cfv_lrrad

    contains

        ! Initialization
        procedure :: Initialize             => InitializeCostfunctionLRrad2

        ! Evaluation
        procedure :: Evaluate               => EvaluateCostFunctionLRrad2

        ! Housekeeping
        procedure :: Allocate               => AllocateCostFunctionLRrad2
        procedure :: Deallocate             => DeallocateCostFunctionLRrad2
        final :: DestroyCostFunctionLRrad2

    end type

    ! Face angle difference cost function
    type, extends(CostfunctionGDUDT) :: CostfunctionFADUDT 

        ! Description
        !============
        ! Cost function based on the angle between magnetic field and face
        ! (and the difference of this angle between two vertex pairs). This
        ! cost function will, ideally, lead to more smooth grids in the
        ! radial direction.
        
        ! Note: this cost function was based on quadrilateral grids, where
        ! the ideal angle difference between two non-aligned faces is zero,
        ! i.e. the faces are perfectly perpendicular to the local magnetic
        ! field. However, for other cell geometries (triangles etc) this
        ! kind of angle is impossible to achieve. However, this cost
        ! function will lead in the ideal case to the lowest variation of
        ! the angle. 

        ! Notes
        !======
        ! Note 1: it is technically not yet possible to impose a desired 
        ! angle difference dtheta0, so the default desired value is 
        ! zero. This may be further extended in the future if needed. 

        ! Note 2: currently (18/11/22), we do not include vertices with
        ! more than one vertex pair to consider. This means in practice
        ! that x-points and vertices that have other cells than 
        ! quadrilateral cells as neighbours are not included. 
    
        ! Cost function formula
        !======================
        !
        !   J_i = wt(i) * ( (theta1 - theta2 )**2 )
        !   J = lambda * sum(J_i), i = 1, vert%ntot
        !   theta1 = atan( cross(d1, B1 )/dot(d1, B1) ) (theta2 similar)

        ! Notes
        !======
        ! Note 1: it is assumed that each vertex only has one vertex 
        ! pair 

        ! Fields
        real(R8)                    :: lambda ! scaling constant
        integer(I8)                 :: nvpairs ! number of vertex pairs
        real(R8), allocatable       :: wt(:) ! weight factor per vertex
        integer(I8), allocatable    :: vpairs(:, :)! vertex pairs

    contains

        ! Initialization
        procedure :: Initialize             => InitializeCostfunctionFAD

        ! Evaluation
        procedure :: Evaluate               => EvaluateCostFunctionFAD

        ! Data output
        procedure :: WriteData              => WriteCostFunctionDataFAD

        ! Housekeeping
        procedure :: Allocate               => AllocateCostFunctionFAD
        procedure :: Deallocate             => DeallocateCostFunctionFAD
        final :: DestroyCostFunctionFAD

    end type

    ! Face angle cost function
    type, extends(CostfunctionGDUDT)    :: CostFunctionFAUDT 

        ! Description
        !============
        ! This cost function penalizes the angle of non-aligned faces
        ! w.r.t the magnetic field  (orthogonality is promoted). This 
        ! can be useful when orthogonality constraints are not possible
        ! to apply, but one still wants faces as orthogonal as possible,
        ! relative to other cost function contributions of course. Only
        ! faces that have:
        ! - a different vertex ID for each vertex (can be zero, but not twice zero)
        ! - at most one vessel vertex (i.e. the face is not a vessel face)
        ! are included in the cost function. Each contribution is stored
        ! as a vertex pair that is used for evaluation later on. 

        ! Notes:
        !=======
        ! Note 1: this cost function can work perfectly in combination
        ! with the FAD cost function: then strong face angle differences
        ! but also strong non-orthogonality is prevented. Using only the
        ! FA cost function will not necessarily lead to better FAD, 
        ! though it may help a bit. 

        ! Fields
        real(R8)                    :: lambda ! scaling constant
        real(R8), allocatable       :: wt(:) ! weight 
        integer(I8), allocatable    :: vpairs(:, :) ! vertex pairs
        integer(I8)                 :: nvpairs ! total number of vertex pairs


    contains

        ! Initialization
        procedure :: Initialize             => InitializeCostfunctionFA

        ! Evaluation
        procedure :: Evaluate               => EvaluateCostFunctionFA

        ! Data output
        procedure :: WriteData              => WriteCostFunctionDataFA

        ! Housekeeping
        procedure :: Allocate               => AllocateCostFunctionFA
        procedure :: Deallocate             => DeallocateCostFunctionFA
        final :: DestroyCostFunctionFA

    end type

    ! Psi ratio, psi based
    type, extends(CostfunctionGDUDT)    :: CostfunctionPRPBUDT 

        ! Description
        !============
        ! This cost function penalizes differences in flux value 
        ! between adjacent flux surfaces. It acts immediately upon the
        ! desired psi values of the flux function constraints, and can
        ! therefore obviously only be included if these constraints are
        ! active. 

        ! The format is very much alike the length ratio cost function. 
        ! Likewise, one best uses the PRPB2 cost function to avoid any
        ! length effects.

        ! Fields
        real(R8)                        :: lambda ! scaling constant
        real(R8), allocatable           :: wt(:) ! weigths 
        real(R8), allocatable           :: b0(:) ! desired bias
        integer(I8), allocatable        :: psipairs(:, :) ! pairs of psi values
        integer(I8)                     :: npsipairs ! number of psi value pairs


    contains

        ! Initialization
        procedure :: Initialize             => InitializeCostfunctionPRPB
        procedure :: FinalizeInitialization => &
            FinalizeInitializationCostFunctionPRPB

        ! Evaluation
        procedure :: Evaluate               => EvaluateCostFunctionPRPB

        ! Data output
        procedure :: WriteData              => WriteCostFunctionDataPRPB

        ! Housekeeping
        procedure :: Allocate               => AllocateCostFunctionPRPB
        procedure :: Deallocate             => DeallocateCostFunctionPRPB
        final :: DestroyCostFunctionPRPB
    end type

    ! Psi ratio, psi based 2 (analogous to LR2)
    type, extends(CostfunctionGDUDT)    :: CostfunctionPRPB2UDT 

        ! Description
        !============
        ! This cost function simply wraps around the PRPB cost function
        ! and symmetrizes it, similar to the lengthratio2 cost function,
        ! to prevent size effects. All real implementation is in the 
        ! PRPB cost function definition.

        ! Fields
        type(CostfunctionPRPBUDT)       :: cfv_prpb


    contains

        ! Initialization
        procedure :: Initialize             => InitializeCostfunctionPRPB2
        procedure :: FinalizeInitialization => &
            FinalizeInitializationCostFunctionPRPB2

        ! Evaluation
        procedure :: Evaluate               => EvaluateCostFunctionPRPB2

        ! Data output
        procedure :: WriteData              => WriteCostFunctionDataPRPB2

        ! Housekeeping
        procedure :: Allocate               => AllocateCostFunctionPRPB2
        procedure :: Deallocate             => DeallocateCostFunctionPRPB2
        final :: DestroyCostFunctionPRPB2
    end type

    ! LRFAD cost function
    type, extends(CostfunctionGDUDT) :: CostfunctionLRFADUDT

        ! Description
        !============
        ! Cost function that combines the length ratio 2 and face angle
        ! difference cost function. The total cost function is simply
        ! compute as the sum of both cost functions (no additional 
        ! weighing is applied here, so all weighing is to be done 
        ! through the lambda values of each cost function separately - 
        ! see also the initialization routine of this cost function). 

        ! Currently, this is the best performing cost function for 
        ! practical purposes (02/12/2022)
    
        ! Cost function formula
        !======================
        ! See formulas of separate cost functions. 

        ! Notes
        !======

        ! Fields
        type(CostfunctionLRUDT2)    :: cfv_lr
        type(CostFunctionFADUDT)    :: cfv_fad

    contains

        ! Initialization
        procedure :: Initialize         => InitializeCostfunctionLRFAD

        ! Evaluation
        procedure :: Evaluate           => EvaluateCostFunctionLRFAD

        ! Housekeeping
        procedure :: Allocate           => AllocateCostFunctionLRFAD
        procedure :: Deallocate         => DeallocateCostFunctionLRFAD
        final :: DestroyCostFunctionLRFAD

    end type

        ! Cost function with all possible contributions
    type, extends(CostfunctionGDUDT) :: CostfunctionGeneralUDT

        ! Description
        !============
        ! Cost function that accounts for all possible combinations of 
        ! length ratio(s), angles, differences, ... The inclusion of a 
        ! cost function value is determined based on the value of the 
        ! scaling coefficient lambda. If this is zero or negative, the 
        ! contribution is not included. One should beware that if the 
        ! lambda values are not properly set in the input file, 
        ! contributions may be unexpectedly included since the default
        ! value for these contributions is non-zero. If no contributions
        ! would be included, the system is likely underdetermined, 
        ! leading to NaNs/divergence of the solver. 

        ! Fields
        type(CostfunctionLRUDT2)        :: cfv_lr
        type(CostFunctionFADUDT)        :: cfv_fad
        type(CostFunctionFAUDT)         :: cfv_fa
        type(CostfunctionPRPB2UDT)      :: cfv_prpb
        type(CostfunctionLRrad2UDT)     :: cfv_lrrad

        ! Switches
        logical                         :: doLR, doFA, doFAD, doPRPB, &
            doLRrad

    contains 

        ! Initialization
        procedure :: Initialize         => InitializeCostfunctionGeneral

        ! Evaluation
        procedure :: Evaluate           => EvaluateCostFunctionGeneral

        ! Housekeeping
        procedure :: Allocate           => AllocateCostFunctionGeneral
        procedure :: Deallocate         => DeallocateCostFunctionGeneral
        final :: DestroyCostFunctionGeneral

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
        subroutine InitializeCostfunctionINT(costfunction, grid, &
            magneticField, environment, options)

            ! Description
            !============
            ! This routine should initialize all additional parameters
            ! that are needed to evaluate the cost function (e.g. the 
            ! vertex indices where the cost function is defined).
            
            ! Import
            import :: CostfunctionGDUDT, GridUDT, MagneticFieldUDT, &
                EnvironmentUDT, CostFunctionOptionsUDT

            ! Declare 
            class(CostfunctionGDUDT)        :: costfunction 
            type(GridUDT)                   :: grid
            type(MagneticFieldUDT)          :: magneticField 
            type(EnvironmentUDT)            :: environment
            type(CostFunctionOptionsUDT)    :: options

        end subroutine

        ! Cost function evaluation
        subroutine EvaluateCostFunctionINT(costfunction, J, gradJ, &
            hessJ, grid, magneticField, environment, dogradient, &
            dohessian, designvariables)

            ! Description
            !============
            ! Main routine to evaluate the cost function and its 
            ! derivative and hessian w.r.t. design variables. 

            ! Import
            import :: CostfunctionGDUDT, MySparseUDT, GridUDT, R8, &
                MagneticFieldUDT, EnvironmentUDT, DesignVariablesGDUDT
            
            ! Declare
            class(CostfunctionGDUDT)        :: costfunction 
            real(R8)                        :: J 
            real(R8), allocatable           :: gradJ(:)
            type(MySparseUDT)               :: hessJ 
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
    !                           LENGTH RATIO                           !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCostFunctionLR(costfunction, grid, &
        magneticField, environment, options)

        ! Description
        !============
        ! Initialize the cost function and its parameters based on the 
        ! grid, magnetic field, and environment structures. Here, the 
        ! length ratio cost function is initialized, which requires

        ! Modules

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionLRUDT)            :: costfunction
        type(GridUDT)                       :: grid
        type(MagneticFieldUDT)              :: magneticField 
        type(EnvironmentUDT)                :: environment 
        type(CostFunctionOptionsUDT)        :: options

        ! Loop variables
        integer(I8)                         :: i, j

        ! Auxiliary variables
        type(StructuredPLF2DDistanceDFUDT)  :: dfbias

        integer                             :: sgn2, sgn3
        integer(I8)                         :: sp, ep, tID, v2, v3, tv
        real(R8)                            :: Btx2, Btx3, Bty2, Bty3, &
            dx2, dy2, dx3, dy3, a0

        integer(I8), allocatable            :: tvn(:), temptvn(:), &
            vesselvert(:), tvf(:), tvfv(:, :)
        real(R8), allocatable               :: Btx(:), Bty(:)
        logical, allocatable                :: cID(:), isvesselvertex(:), &
            isvesselface(:)

        ! Data
        
        ! Initialize
        !===========
        ! Set the scaling constant
        costfunction%lambda = options%LR%lambda
        
        ! Set the small length parameter eta
        costfunction%eta    = options%LR%eta

        ! Check if we need to include vessel edges
        costfunction%dovessel   = options%LR%dovessel

        ! Allocate
        call costfunction%Allocate(grid%vert%ntot, 4)
        allocate(Btx(grid%vert%ntot))
        allocate(Bty(grid%vert%ntot))

        ! Associate some fields
        associate(vert => grid%vert, x => grid%vert%x, y => grid%vert%y, &
            b0 => costfunction%b0, wt => costfunction%wt, &
            nvpairs => costfunction%nvpairs, &
            vpairs => costfunction%vpairs, &
            dovessel => costfunction%dovessel)

        ! Set the initial weighting factors
        wt(:) = 1

        ! Initialize
        vpairs(:, :) = 0
        nvpairs(:) = 0

        ! Compute the magnetic field vectors at the vertex locations
        call magneticField%interp%Evaluate(x, y, 0, 1, Btx)
        call magneticField%interp%Evaluate(x, y, 1, 0, Bty)
        Btx = -Btx ! adjust sign: Btx = - dPsidy

        ! Compute desired length ratio
        !=============================
        ! Desired bias far away from vessel (set to one)
        a0 = 1

        ! Construct polygon set based on target plates only
        !call targetps%Construct(environment%vessel%targetpolygons)

        ! Construct interpolant
        !call targetinterp%SetParameters('uniformgrid', 3, 6)
        !call targetinterp%Construct()

        ! Construct
        call dfbias%Initialize(magneticField%interp, &
            environment%vessel%plftarget, environment%vessel%plfvessel, options%LR%biasatvessel, a0, &
            options%LR%lengthparam, 'signed')
        
        ! Evaluate
        call dfbias%Evaluate(x, y, b0)

        ! Visualize
        call environment%vessel%plftarget%Visualize('costfunctionLR_vesselcontours')
        call dfbias%Visualize([minval(x), maxval(x)], [minval(y), maxval(y)], 100, 100, 'costfunctionLR_desiredbias')

        ! Write
        !call Write3DCoordinateData(x, y, b0, 'costfunctionLR_desiredbias')

        ! Compute the vertex pairs
        !=========================
        do i = 1, grid%vert%ntot 
            ! Housekeeping
            allocate(temptvn(vert%neigP(i, 2)))
            allocate(cID(vert%neigP(i, 2)))

            ! Get the vertex neighbours of this vertex
            temptvn = GetvertNeig(vert, i)
            sp = vert%neigP(i, 1)
            ep = vert%neigP(i, 1) + vert%neigP(i, 2)-1
            temptvn = vert%neig(sp:ep)
            
            ! Get the ID of the coordinate line
            tID = vert%fieldlineID(i)
            
            ! Check which vertices have the same ID
            if (tID /= 0) then 
                cID = (tID == vert%fieldlineID(temptvn))
            else
                ! Determine later
                cID(:) = .true. 
            end if
            
            ! Extract
            allocate(tvn(count(cID)))
            tvn = pack(temptvn, cID)
            
            ! Assemble
            if (size(tvn) > 0) then
                ! Update counter
                nvpairs(i) = (size(tvn)/2)
                
                ! If multiple pairs, loop
                if(tID /= 0) then 
                    do j = 1, nvpairs(i)
                        ! Normally, multiple pairs only occur at x-points, and,
                        ! since the coordinates are sorted, the corresponding
                        ! pairs should be tvn(j) and tvn(j+nvpairs(i))
                        
                        ! Get vertices
                        v2 = tvn(j)
                        v3 = tvn(j+nvpairs(i))
                    
                    ! Get vectors
                        dx2 = x(v2) - x(i); dy2 = y(v2) - y(i)
                        dx3 = x(v3) - x(i); dy3 = y(v3) - y(i)
                        
                        ! Check if we're dealing with an x-point
                        if (nvpairs(i) > 1) then
                            ! Here, the gradient *should* vanish. For now, we
                            ! cope with this by setting the desired ratio to 1,
                            ! such that it does not matter which length is
                            ! considered first.
                            b0(i) = 1
                            sgn2 = -1
                            wt(i) = 0
                        else
                            ! Evaluate sign of dot product of magnetic field
                            ! coordinates with vector
                            Btx2 = 0.5*(Btx(i) + Btx(v2))
                            Bty2 = 0.5*(Bty(i) + Bty(v2))
                            Btx3 = 0.5*(Btx(i) + Btx(v3))
                            Bty3 = 0.5*(Bty(i) + Bty(v3))
                            if ( (dx2*Btx2 + dy2*Bty2) < 0 ) then 
                                sgn2 = -1
                            else 
                                sgn2 = 1
                            end if
                            if ( (dx3*Btx3 + dy3*Bty3) < 0 ) then 
                                sgn3 = -1
                            else 
                                sgn3 = 1
                            end if
                            
                            ! Consistency check: normally, one positive and one
                            ! negative sign
                            if ( ((sgn2 < 0) .and. (sgn3 < 0)) &
                                .or. ((sgn2 > 0) .and. (sgn3 > 0)) ) then
                                ! Most likely we're near an x-point here, so
                                ! the magnetic field is off. Ignore these
                                ! vertices
                                
                                b0(i) = 1
                                wt(i) = 0
                                sgn2 = -1
                            
                            end if
                        end if
                        
                        ! Add vertices in the direction along the coordinate
                        ! line
                        if (sgn2 < 0) then
                            vpairs(i,2*j-1:2*j) = (/v2, v3/)
                        else
                            vpairs(i,2*j-1:2*j) = (/v3, v2/)
                        end if
                        
                    end do
                else
                    ! Vertex without fieldline ID - don't include. 
                    ! It is assumed that these vertices only appear on 
                    ! the boundary, and these vertices are considered
                    ! later on if vessel edges are considered.
                    

                end if
            end if

            ! Housekeeping
            deallocate(tvn, temptvn, cID)

        end do

        ! Include vessel vertices?
        if (dovessel) then 
            ! Get all vessel vertices
            call DetermineVesselVertices(isvesselvertex, isvesselface, grid)
            allocate(vesselvert(count(isvesselvertex)))
            vesselvert = pack([(i, i=1, grid%vert%ntot)], isvesselvertex)

            ! Overwrite potential other vertex pairs (ordering doesn't 
            ! matter because we set bias to one anyway)
            do i = 1, size(vesselvert, 1)
                ! Unpack
                tv = vesselvert(i)

                ! Get the faces of this vertex
                tvf = GetVertFace(vert, tv)

                ! Check
                if (count(isvesselface(tvf)) == 2) then
                    ! Get the other two vertices
                    allocate(tvfv(2, 2))
                    tvfv = grid%face%vert(pack(tvf, isvesselface(tvf)), :)
                    allocate(tvn(count(tvfv /= tv)))
                    tvn = pack(tvfv, tvfv /= tv)

                    ! Check
                    if (size(tvn, 1) /= 2) then 
                        ! Shouldn't happen, throw error
                        call gdErrorHandler('InitializeCostFunctionLR: ' // &
                            'unknown error, something seems wrong in grid interconnection')
                    end if 

                    ! Add
                    vpairs(tv, 1:2) = tvn
                    b0(tv) = 1
                    nvpairs(tv) = 1 

                    ! Deallocate 
                    deallocate(tvfv, tvf, tvn)
                end if 

            end do

            ! Housekeeping
            deallocate(vesselvert)


        end if 

        ! Housekeeping
        deallocate(Btx, Bty)

        ! End associate
        end associate

        ! Write data
        !===========
        if (options%writedata == 1) then 
            call costfunction%WriteData(grid)
        end if 

    end subroutine

    ! Cost function evaluation
    subroutine EvaluateCostFunctionLR(costfunction, J, gradJ, hessJ, &
        grid, magneticField, environment, dogradient, dohessian, &
        designvariables)

        ! Description
        !============
        ! Evaluate the cost function, the gradient and its hessian. 

        ! Notes:
        !=======
        ! Possible future performance improvements:
        ! - Allocating hessian stuff only once and storing indices, 
        ! since they don't change
        ! - Instead of recomputing auxiliary variables, store them. May
        ! not actually be better in terms of computational time, but 
        ! may lead to shorter and hence better maintainable code. 

        ! Note: we've added a relaxation term in the denominator (so for
        ! d2) to prevent NaN/Inf behavior when points coincide. This does
        ! lead to higher gradients, but ensure differentiability. Additionally,
        ! if the parameter in the denominator, eta, is relatively large 
        ! compared to d2, the obtained ratio will (perhaps strongly) differ
        ! from the desired one locally. However, if chosen too small, the 
        ! gradient may become extremely large leading to difficult to solve
        ! subproblems

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionLRUDT)        :: costfunction 
        real(R8)                        :: J
        real(R8), allocatable           :: gradJ(:) 
        type(MySparseUDT)               :: hessJ 
        type(GridUDT)                   :: grid 
        type(MagneticFieldUDT)          :: magneticField 
        type(EnvironmentUDT)            :: environment
        logical                         :: dogradient, dohessian 
        class(DesignVariablesGDUDT)     :: designvariables

        ! Loop variables
        integer(I8)                     :: i, k, cc 

        ! Auxiliary
        integer(I8)                     :: v1, v2, v3
        real(R8)                        :: x1, x2, x3, y1, y2, y3, wti, b0i
        real(R8)                        :: dx1, dx2, dy1, dy2, d1, d2, rat
        integer(I8), allocatable        :: row(:), col(:) 
        real(R8), allocatable           :: valxx(:),  valxy(:), &
                                        valyx(:), valyy(:)
                                        
        ! Associate
        !==========
        associate(&
            vert    => grid%vert, &
            vpairs  => costfunction%vpairs, &
            nvpairs => costfunction%nvpairs, &
            b0      => costfunction%b0, &
            x       => grid%vert%x, &
            y       => grid%vert%y, &
            lambda  => costfunction%lambda, &
            wt      => costfunction%wt, &
            eta     => costfunction%eta )

        ! Initialize
        !===========
        ! Cost function
        J = 0

        ! Gradient
        gradJ(:) = 0

        ! Compute cost function
        !======================
        ! Loop over all vertices
        do i = 1, vert%ntot
            ! Set the current vertex
            v1 = i 

            ! Loop over all vertex pairs of this vertex
            do k = 1, nvpairs(i)
                ! Get the current pair
                v2 = vpairs(i, 2*k-1)
                v3 = vpairs(i, 2*k)

                ! Compute intermediate quantities
                x1 = x(v1)
                x2 = x(v2)
                x3 = x(v3)

                y1 = y(v1)
                y2 = y(v2)
                y3 = y(v3)

                dx1 = x2 - x1
                dx2 = x3 - x1 

                dy1 = y2 - y1
                dy2 = y3 - y1 

                ! Compute lengths
                d1 = sqrt(dx1**2 + dy1**2)
                d2 = sqrt(dx2**2 + dy2**2)

                ! Compute cost function contribution
                J = J + 0.5*wt(i)*(d1/(d2 + eta) - b0(i))**2

            end do 
        end do

        ! Scale
        J = lambda*J

        ! Compute gradient
        !=================
        if (dogradient) then 

            ! Check the design variables
            select case (trim(designvariables%type))

            case ('coordinates', 'coordinates_desiredflux')

                ! Loop over all vertices
                do i = 1, vert%ntot
                    ! Set the current vertex
                    v1 = i 

                    ! Unpack
                    wti = wt(i)
                    b0i = b0(i)

                    ! Loop over all vertex pairs of this vertex
                    do k = 1, nvpairs(i)
                        ! Get the current pair
                        v2 = vpairs(i, 2*k-1)
                        v3 = vpairs(i, 2*k)

                        ! Compute intermediate quantities
                        x1 = x(v1)
                        x2 = x(v2)
                        x3 = x(v3)

                        y1 = y(v1)
                        y2 = y(v2)
                        y3 = y(v3)

                        dx1 = x2 - x1
                        dx2 = x3 - x1 

                        dy1 = y2 - y1
                        dy2 = y3 - y1 

                        ! Compute lengths
                        d1 = sqrt(dx1**2 + dy1**2)
                        d2 = sqrt(dx2**2 + dy2**2)

                        ! Compute ratio
                        rat = d1/(d2 + eta)

                        ! Compute cost function contribution
                        gradJ(v1) = gradJ(v1) + & 
                            wti*(b0i - rat)*(dx1/(d1*(d2 + eta)) - &
                            (dx2*rat)/(d2*(d2 + eta))) ! x1
                        gradJ(v2) = gradJ(v2) + &
                            -(dx1*wti*(b0i - rat))/(d1*(d2 + eta)) ! x2
                        gradJ(v3) = gradJ(v3) + &
                            (dx2*rat*wti*(b0i - rat))/(d2*(d2 + eta)) ! x3
                        
                        gradJ(v1+vert%ntot) = gradJ(v1+vert%ntot) + &
                            wti*(b0i - rat)*(dy1/(d1*(d2 + eta)) &
                            - (dy2*rat)/(d2*(d2 + eta))) ! y1
                        gradJ(v2+vert%ntot) = gradJ(v2+vert%ntot) + &
                            -(dy1*wti*(b0i - rat))/(d1*(d2 + eta)) ! y2
                        gradJ(v3+vert%ntot) = gradJ(v3+vert%ntot) + &
                            (dy2*rat*wti*(b0i - rat))/(d2*(d2 + eta)) ! y3

                    end do 
                end do

            case default

                ! Not implemented, throw error
                call gdErrorHandler('EvaluateCostFunctionLR: gradient' &
                    // ' not yet implemented for this design variable' &
                    // ' type')

            end select

            ! Scale
            gradJ = lambda*gradJ

        end if

        ! Compute hessian
        !================
        ! Allocate the hessian (if not already done so)
        if (.not. allocated(hessJ%row)) then
            ! Allocate the sparse matrix
            select case (trim(designvariables%type))

            case ('coordinates', 'coordinates_desiredflux')

                hessJ%nval = 36*sum(nvpairs) ! this should be exact and constant

            end select 
            call hessJ%Allocate()

        end if

        if (dohessian) then

            ! Check the design variables
            select case (trim(designvariables%type))

            case ('coordinates', 'coordinates_desiredflux')

                ! Allocate auxiliary variables
                allocate(row(9*sum(nvpairs)))
                allocate(col(9*sum(nvpairs)))
                allocate(valxx(9*sum(nvpairs)))
                allocate(valxy(9*sum(nvpairs)))
                allocate(valyx(9*sum(nvpairs)))
                allocate(valyy(9*sum(nvpairs)))

                ! Initialize counter
                cc = 1

                ! Loop over all vertices
                do i = 1, vert%ntot
                    ! Set the current vertex
                    v1 = i 

                    ! Unpack
                    b0i = b0(i)
                    wti = wt(i)

                    ! Loop over all vertex pairs of this vertex
                    do k = 1, nvpairs(i)
                        ! Get the current pair
                        v2 = vpairs(i, 2*k-1)
                        v3 = vpairs(i, 2*k)

                        ! Compute intermediate quantities
                        x1 = x(v1)
                        x2 = x(v2)
                        x3 = x(v3)

                        y1 = y(v1)
                        y2 = y(v2)
                        y3 = y(v3)

                        dx1 = x2 - x1
                        dx2 = x3 - x1 

                        dy1 = y2 - y1
                        dy2 = y3 - y1 

                        ! Compute lengths
                        d1 = sqrt(dx1**2 + dy1**2)
                        d2 = sqrt(dx2**2 + dy2**2)

                        ! Compute ratio
                        rat = d1/(d2 + eta) 

                        ! Compute Hessian contributions - ordened per
                        ! vertex pair (e.g. v1, v1), split up in xx, xy,
                        ! yx, yy. 

                        ! d2J/dx1**2, d2J/dy1**2, d2J/dx1dy1, d2J/dy1dx1
                        row(cc) = v1; col(cc) = v1
                        valxx(cc) = wti*(dx1/(d1*(d2 + eta)) - &
                            (dx2*rat)/(d2*(d2 + eta)))**2 - &
                            wti*(b0i - rat)*(1/(d1*(d2 + eta)) - &
                            rat/(d2*(d2 + eta)) - &
                            dx1**2/(d1**3*(d2 + eta)) + &
                            (2*dx2**2*rat)/(d2**2*(d2 + eta)**2) + &
                            (dx2**2*rat)/(d2**3*(d2 + eta)) - &
                            (2*dx1*dx2)/(d1*d2*(d2 + eta)**2)) ! x1x1

                        valyy(cc) = wti*(dy1/(d1*(d2 + eta)) - &
                            (dy2*rat)/(d2*(d2 + eta)))**2 - &
                            wti*(b0i - rat)*(1/(d1*(d2 + eta)) &
                            - rat/(d2*(d2 + eta)) - &
                            dy1**2/(d1**3*(d2 + eta)) + &
                            (2*dy2**2*rat)/(d2**2*(d2 + eta)**2) + &
                            (dy2**2*rat)/(d2**3*(d2 + eta)) - &
                            (2*dy1*dy2)/(d1*d2*(d2 + eta)**2)) !y1y1

                        valxy(cc) = wti*(b0i - rat)*&
                            ((dx1*dy1)/(d1**3*(d2 + eta)) + &
                            (dx1*dy2)/(d1*d2*(d2 + eta)**2) + &
                            (dx2*dy1)/(d1*d2*(d2 + eta)**2) - &
                            (2*dx2*dy2*rat)/(d2**2*(d2 + eta)**2) - &
                            (dx2*dy2*rat)/(d2**3*(d2 + eta))) + &
                            wti*(dx1/(d1*(d2 + eta)) - &
                            (dx2*rat)/(d2*(d2 + eta)))*&
                            (dy1/(d1*(d2 + eta)) - &
                            (dy2*rat)/(d2*(d2 + eta))) !x1y1

                        valyx(cc) = valxy(cc)
                        cc = cc+1
                        
                        ! d2J/dx1dx2, d2J/dy1dy2, d2J/dx1dy2, d2J/dy1dx2
                        row(cc) = v1; col(cc) = v2
                        valxx(cc) = - wti*(b0i - rat)*&
                            (dx1**2/(d1**3*(d2 + eta)) - &
                            1/(d1*(d2 + eta)) + &
                            (dx1*dx2)/(d1*d2*(d2 + eta)**2)) - &
                            (dx1*wti*(dx1/(d1*(d2 + eta)) - &
                            (dx2*rat)/(d2*(d2 + eta))))/(d1*(d2 + eta)) !x1x2
                        valyy(cc) = - wti*(b0i - rat)*&
                            (dy1**2/(d1**3*(d2 + eta)) - &
                            1/(d1*(d2 + eta)) + (dy1*dy2)&
                            /(d1*d2*(d2 + eta)**2)) - &
                            (dy1*wti*(dy1/(d1*(d2 + eta)) - &
                            (dy2*rat)/(d2*(d2 + eta))))/(d1*(d2 + eta)) !y1y2
                        valxy(cc) = - wti*((dx1*dy1)/(d1**3*(d2 + eta)) &
                            + (dx2*dy1)/(d1*d2*(d2 + eta)**2))*(b0i - rat) &
                            - (dy1*wti*(dx1/(d1*(d2 + eta)) - &
                            (dx2*rat)/(d2*(d2 + eta))))/(d1*(d2 + eta)) !x1y2
                        valyx(cc) = - wti*((dx1*dy1)/(d1**3*(d2 + eta)) &
                            + (dx1*dy2)/(d1*d2*(d2 + eta)**2))*(b0i - rat) &
                            - (dx1*wti*(dy1/(d1*(d2 + eta)) &
                            - (dy2*rat)/(d2*(d2 + eta))))/(d1*(d2 + eta)) !y1x2
                        cc = cc+1
                        
                        row(cc) = v2; col(cc) = v1
                        valxx(cc) = valxx(cc-1) !x2x1
                        valyy(cc) = valyy(cc-1) !y2y1
                        valxy(cc) = valyx(cc-1) !x1y2
                        valyx(cc) = valxy(cc-1) !y2x1
                        cc = cc+1
                        
                        row(cc) = v2; col(cc) = v2
                        valxx(cc) = (dx1**2*wti)/(d1**2*(d2 + eta)**2) &
                            - (wti*(b0i - rat))/(d1*(d2 + eta)) &
                            + (dx1**2*wti*(b0i - rat))/(d1**3*(d2 + eta)) ! x2x2
                        valyy(cc) = (dy1**2*wti)/(d1**2*(d2 + eta)**2) &
                            - (wti*(b0i - rat))/(d1*(d2 + eta)) &
                            + (dy1**2*wti*(b0i - rat))/(d1**3*(d2 + eta)) ! y2y2
                        valxy(cc) = (dx1*dy1*wti)/(d1**2*(d2 + eta)**2) &
                            + (dx1*dy1*wti*(b0i - rat))/(d1**3*(d2 + eta)) !x2y2
                        valyx(cc) = valxy(cc)
                        cc = cc+1
                        
                        row(cc) = v1; col(cc) = v3
                        valxx(cc) = (dx2*rat*wti*(dx1/(d1*(d2 + eta)) - &
                            (dx2*rat)/(d2*(d2 + eta))))/(d2*(d2 + eta)) &
                            - wti*(b0i - rat)*(rat/(d2*(d2 + eta)) &
                            - (2*dx2**2*rat)/(d2**2*(d2 + eta)**2) &
                            - (dx2**2*rat)/(d2**3*(d2 + eta)) &
                            + (dx1*dx2)/(d1*d2*(d2 + eta)**2)) !x1x3
                        valyy(cc) = (dy2*rat*wti*(dy1/(d1*(d2 + eta)) &
                            - (dy2*rat)/(d2*(d2 + eta))))/(d2*(d2 + eta)) &
                            - wti*(b0i - rat)*(rat/(d2*(d2 + eta)) &
                            - (2*dy2**2*rat)/(d2**2*(d2 + eta)**2) &
                            - (dy2**2*rat)/(d2**3*(d2 + eta)) &
                            + (dy1*dy2)/(d1*d2*(d2 + eta)**2)) !y1y3
                        valxy(cc) = wti*(b0i - rat)*&
                            ((2*dx2*dy2*rat)/(d2**2*(d2 + eta)**2) &
                            - (dx1*dy2)/(d1*d2*(d2 + eta)**2) &
                            + (dx2*dy2*rat)/(d2**3*(d2 + eta))) &
                            + (dy2*rat*wti*(dx1/(d1*(d2 + eta)) &
                            - (dx2*rat)/(d2*(d2 + eta))))/(d2*(d2 + eta)) !x1y3
                        valyx(cc) = wti*(b0i - rat)*&
                            ((2*dx2*dy2*rat)/(d2**2*(d2 + eta)**2) &
                            - (dx2*dy1)/(d1*d2*(d2 + eta)**2) &
                            + (dx2*dy2*rat)/(d2**3*(d2 + eta))) &
                            + (dx2*rat*wti*(dy1/(d1*(d2 + eta)) &
                            - (dy2*rat)/(d2*(d2 + eta))))/(d2*(d2 + eta)) !y1x3
                        cc = cc+1
                        
                        row(cc) = v3; col(cc) = v1
                        valxx(cc) = valxx(cc-1)
                        valyy(cc) = valyy(cc-1)
                        valxy(cc) = valyx(cc-1)
                        valyx(cc) = valxy(cc-1)
                        cc = cc+1
                        
                        row(cc) = v2; col(cc) = v3
                        valxx(cc) = (dx1*dx2*wti*(b0i - rat))/&
                            (d1*d2*(d2 + eta)**2) - (dx1*dx2*wti)/(d2*(d2 + eta)**3) !x2x3
                        valyy(cc) = (dy1*dy2*wti*(b0i - rat))/(d1*d2*(d2 + eta)**2) &
                            - (dy1*dy2*wti)/(d2*(d2 + eta)**3) !y2y3
                        valxy(cc) = (dx1*dy2*wti*(b0i - rat))/(d1*d2*(d2 + eta)**2) &
                            - (dx1*dy2*wti)/(d2*(d2 + eta)**3) !x2y3
                        valyx(cc) = (dx2*dy1*wti*(b0i - rat))/(d1*d2*(d2 + eta)**2) &
                            - (dx2*dy1*wti)/(d2*(d2 + eta)**3) !y2x3
                        cc = cc+1
                        
                        row(cc) = v3; col(cc) = v2
                        valxx(cc) = valxx(cc-1)
                        valyy(cc) = valyy(cc-1)
                        valxy(cc) = valyx(cc-1)
                        valyx(cc) = valxy(cc-1)
                        cc = cc+1
                        
                        row(cc) = v3; col(cc) = v3
                        valxx(cc) = (rat*wti*(b0i - rat))/(d2*(d2 + eta)) &
                            + (dx2**2*wti*(dx1**2 + dy1**2))/(d2**2*(d2 + eta)**4) &
                            - (2*dx2**2*rat*wti*(b0i - rat))/(d2**2*(d2 + eta)**2) &
                            - (dx2**2*rat*wti*(b0i - rat))/(d2**3*(d2 + eta)) !x3x3
                        valyy(cc) = (rat*wti*(b0i - rat))/(d2*(d2 + eta)) &
                            + (dy2**2*wti*(dx1**2 + dy1**2))/(d2**2*(d2 + eta)**4) &
                            - (2*dy2**2*rat*wti*(b0i - rat))/(d2**2*(d2 + eta)**2) &
                            - (dy2**2*rat*wti*(b0i - rat))/(d2**3*(d2 + eta)) !y3y3
                        valxy(cc) = (dx2*dy2*wti*(dx1**2 + dy1**2))/(d2**2*(d2 + eta)**4) &
                            - (2*dx2*dy2*rat*wti*(b0i - rat))/(d2**2*(d2 + eta)**2) &
                            - (dx2*dy2*rat*wti*(b0i - rat))/(d2**3*(d2 + eta)) !x3y3
                        valyx(cc) = valxy(cc)
                        cc = cc+1

                    end do 
                end do

                ! Build full hessian
                hessJ%row = [row, row, row+vert%ntot, row+vert%ntot]
                hessJ%col = [col, col+vert%ntot, col, col+vert%ntot]
                hessJ%val = [valxx, valxy, valyx, valyy]

                ! Scale
                hessJ%val = lambda*hessJ%val

                ! Housekeeping
                deallocate(row)
                deallocate(col)
                deallocate(valxx)
                deallocate(valxy)
                deallocate(valyx)
                deallocate(valyy)


            case default

                ! Not implemented, throw error
                call gdErrorHandler('EvaluateCostFunctionLR: hessian' &
                    // ' not yet implemented for this design variable' &
                    // ' type')

            end select

        end if

        

        ! Deassociate
        !============
        end associate

    end subroutine

    ! Cost function data writing 
    subroutine WriteCostFunctionDataLR(costfunction, grid)

        ! Description
        !============
        ! Write out the cost function data for the LR cost function.
        ! Here, this consists of the vertex pair data in IDn, xn, yn 
        ! format

        ! Declare variables
        !==================
        ! Arguments
        class(CostFunctionLRUDT)        :: costfunction 
        type(GridUDT)                   :: grid

        ! Auxiliary
        integer(I8)                     :: ncol, nrow 

        integer(I8), allocatable        :: IDn(:, :) 
        real(R8), allocatable           :: xn(:, :), yn(:, :)
        character(:), allocatable       :: filename 

        ! Loop
        integer(I8)                     :: j 

        ! Initialize
        !===========
        ! Set filename
        allocate(character(len('costfunction_vertexpairs_LR')) :: filename)
        filename = 'costfunction_vertexpairs_LR'

        ! Allocate
        nrow = size(costfunction%vpairs, 1)
        ncol = size(costfunction%vpairs, 2)
        allocate(IDn(nrow, ncol), xn(nrow, ncol), yn(nrow, ncol))

        ! Unpack
        associate(&
            vpairs      => costfunction%vpairs,         &
            x           => grid%vert%x,                 &
            y           => grid%vert%y)

        ! Loop
        do j = 1, ncol 
            IDn(:, j) = vpairs(:, j) 
            xn(:, j) = x(vpairs(:, j)) 
            yn(:, j) = y(vpairs(:, j)) 
        end do

        ! Call writer
        !============
        call WriteVertexPairData(IDn, xn, yn, filename)

        ! Housekeeping
        !=============
        end associate
        deallocate(IDn, xn, yn)
        


    end subroutine

    ! Housekeeping
    subroutine AllocateCostFunctionLR(costfunction, nv, nvn)

        ! Description
        !============
        ! Allocate, assumed that costfunction%nvpairs is given

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionLRUDT)        :: costfunction
        integer(I8)                     :: nv, nvn

        ! Allocate
        !=========
        allocate(costfunction%nvpairs(nv))
        allocate(costfunction%vpairs(nv,2*nvn))
        allocate(costfunction%b0(nv))
        allocate(costfunction%wt(nv))

    end subroutine

    subroutine DeallocateCostFunctionLR(costfunction)

        ! Description
        !============
        ! Deallocate

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionLRUDT)        :: costfunction

        ! Deallocate
        !===========
        if (allocated(costfunction%vpairs)) then 
            deallocate(costfunction%vpairs)
            deallocate(costfunction%b0)
            deallocate(costfunction%wt)
        end if

    end subroutine

    subroutine DestroyCostFunctionLR(costfunction)

        ! Description
        !============
        ! Deallocate

        ! Declare variables
        !==================
        ! Arguments
        type(CostfunctionLRUDT)        :: costfunction

        ! Destroy
        !========
        call costfunction%Deallocate()

    end subroutine

    !------------------------------------------------------------------!
    !                        LENGTH RATIO, RADIAL                      !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCostFunctionLRrad(costfunction, grid, &
        magneticField, environment, options)

        ! Description
        !============
        ! Initialize the cost function and its parameters based on the 
        ! grid, magnetic field, and environment structures. Here, the 
        ! length ratio cost function is initialized, which requires

        ! Modules

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionLRradUDT)         :: costfunction
        type(GridUDT)                       :: grid
        type(MagneticFieldUDT)              :: magneticField 
        type(EnvironmentUDT)                :: environment 
        type(CostFunctionOptionsUDT)        :: options

        ! Loop variables
        integer(I8)                         :: i, j

        ! Auxiliary variables
        integer                             :: sgn2, sgn3
        integer(I8)                         :: sp, ep, tID, v2, v3, tv, &
            tempID(1:2)
        real(R8)                            :: Btx2, Btx3, Bty2, Bty3, &
            dx2, dy2, dx3, dy3

        integer(I8), allocatable            :: tvn(:), temptvn(:), &
            vesselvert(:), tvf(:), tvfv(:, :), tvnID(:)
        real(R8), allocatable               :: Btx(:), Bty(:)

        logical                             :: includecutcellfaces, &
            excludedomainfaces
        logical, allocatable                :: cID(:), isvesselvertex(:), &
            isvesselface(:)

        ! Data
        
        ! Initialize
        !===========
        ! Set the scaling constant
        costfunction%lambda = options%LRrad%lambda
        
        ! Set the small length parameter eta
        costfunction%eta    = options%LRrad%eta

        ! Include cut cell faces?
        includecutcellfaces = options%LRrad%includecutcellfaces 

        ! Exclude domain faces?
        excludedomainfaces  = options%LRrad%excludedomainfaces

        ! Allocate
        call costfunction%Allocate(grid%vert%ntot, 4)
        allocate(Btx(grid%vert%ntot))
        allocate(Bty(grid%vert%ntot))

        ! Associate some fields
        associate(vert => grid%vert, x => grid%vert%x, y => grid%vert%y, &
            b0 => costfunction%b0, wt => costfunction%wt, &
            nvpairs => costfunction%nvpairs, &
            vpairs => costfunction%vpairs, &
            dovessel => costfunction%dovessel)

        ! Set the initial weighting factors
        wt(:) = 1

        ! Initialize
        vpairs(:, :) = 0
        nvpairs(:) = 0

        ! Compute the magnetic field vectors at the vertex locations
        call magneticField%interp%Evaluate(x, y, 0, 1, Btx)
        call magneticField%interp%Evaluate(x, y, 1, 0, Bty)
        Btx = -Btx ! adjust sign: Btx = - dPsidy

        ! Compute desired length ratio
        !=============================
        ! This has to be replaced with an actual computation based on 
        ! the magnetic field... Right now, simply ones
        b0(:) = 1

        ! Compute the vertex pairs
        !=========================
        do i = 1, grid%vert%ntot 
            ! Housekeeping
            allocate(tvn(vert%neigP(i, 2)))
            allocate(cID(vert%neigP(i, 2)))

            ! Get the vertex neighbours of this vertex
            tvn = GetvertNeig(vert, i)
            sp = vert%neigP(i, 1)
            ep = vert%neigP(i, 1) + vert%neigP(i, 2)-1
            tvn = vert%neig(sp:ep)
            
            ! Get the ID of the coordinate line
            tID = vert%fieldlineID(i)
            tvnID = vert%fieldlineID(tvn)
            
            ! Assemble
            if ( (size(tvn) > 0) .and. (tID /= 0)) then
                ! Vertices with different ID?
                if (includecutcellfaces) then 
                    ! Zeros allowed if boundary vertex
                    cID = ( (tID /= tvnID) .and. ( (tvnID /= 0) .or. (vert%BV(tvn))) )
                else 
                    ! No zeros allowed
                    cID = ( (tID /= tvnID) .and. (tvnID /= 0) )
                end if

                ! Are two vertices remaining?
                if (count(cID) == 2) then 
                    ! Do both vertices have a different ID?
                    tempID = pack(tvnID, cID)
                    if (tempID(1) == tempID(2)) then 
                        ! Don't include
                        cID = .false. 
                    elseif ( (.not. any(tempID == 0)) .and. excludedomainfaces) then 
                        ! Don't include
                        cID = .false.
                    end if
                else
                    ! Skip
                    cID = .false.
                end if
            else
                ! Determine later
                if (allocated(cID)) then 
                    deallocate(cID)
                end if
                allocate(cID(size(tvnID)))
                cID = .true.
            end if
            
            ! Extract
            allocate(temptvn(count(cID)))
            temptvn = pack(tvn, cID)
            tvn = temptvn

            ! Assemble
            if (size(tvn) > 0) then 

                ! Update counter
                nvpairs(i) = (size(tvn)/2)

                if ((tID /= 0) .and. (.not. vert%BV(i))) then ! regular vertex with fieldline ID
                    do j = 1, nvpairs(i)
                        
                        ! Normally, multiple pairs only occur at x-points, and,
                        ! since the coordinates are sorted, the corresponding
                        ! pairs should be tvn(j) and tvn(j+nvpairs(i))
                        
                        ! Get vertices
                        v2 = tvn(j)
                        v3 = tvn(j+nvpairs(i))
                    
                        ! Get vectors
                        dx2 = x(v2) - x(i); dy2 = y(v2) - y(i)
                        dx3 = x(v3) - x(i); dy3 = y(v3) - y(i)
                        
                        ! Check if we're dealing with an x-point
                        if (nvpairs(i) > 1) then
                            ! Here, the gradient *should* vanish. For now, we
                            ! cope with this by setting the desired ratio to 1,
                            ! such that it does not matter which length is
                            ! considered first.
                            b0(i) = 1
                            sgn2 = -1
                            wt(i) = 0
                        else
                            ! Evaluate sign of dot product of magnetic field
                            ! coordinates with vector
                            Btx2 = 0.5*(Btx(i) + Btx(v2))
                            Bty2 = 0.5*(Bty(i) + Bty(v2))
                            Btx3 = 0.5*(Btx(i) + Btx(v3))
                            Bty3 = 0.5*(Bty(i) + Bty(v3))
                            if ( (dx2*Btx2 + dy2*Bty2) < 0 ) then 
                                sgn2 = -1
                            else 
                                sgn2 = 1
                            end if
                            if ( (dx3*Btx3 + dy3*Bty3) < 0 ) then 
                                sgn3 = -1
                            else 
                                sgn3 = 1
                            end if
                            
                            ! Consistency check: normally, one positive and one
                            ! negative sign
                            if ( ((sgn2 < 0) .and. (sgn3 < 0)) &
                                .or. ((sgn2 > 0) .and. (sgn3 > 0)) ) then
                                ! Most likely we're near an x-point here, so
                                ! the magnetic field is off. Ignore these
                                ! vertices
                                
                                b0(i) = 1
                                wt(i) = 0
                                sgn2 = -1
                            
                            end if
                        end if
                        
                        ! Add vertices in the direction along the coordinate
                        ! line
                        if (sgn2 < 0) then
                            vpairs(i,2*j-1:2*j) = (/v2, v3/)
                        else
                            vpairs(i,2*j-1:2*j) = (/v3, v2/)
                        end if
                        
                    end do
                else
                    ! Vertex without fieldline ID - don't include. 
                    ! It is assumed that these vertices only appear on 
                    ! the boundary, and these vertices are considered
                    ! later on if vessel edges are considered.
                    nvpairs(i)  = 0
                    wt(i)       = 0

                end if
            end if

            ! Housekeeping
            deallocate(tvn, temptvn, cID)

        end do

        ! Include vessel vertices?
        if (dovessel) then 
            ! Get all vessel vertices
            call DetermineVesselVertices(isvesselvertex, isvesselface, grid)
            allocate(vesselvert(count(isvesselvertex)))
            vesselvert = pack([(i, i=1, grid%vert%ntot)], isvesselvertex)

            ! Overwrite potential other vertex pairs (ordering doesn't 
            ! matter because we set bias to one anyway)
            do i = 1, size(vesselvert, 1)
                ! Unpack
                tv = vesselvert(i)

                ! Get the faces of this vertex
                tvf = GetVertFace(vert, tv)

                ! Check
                if (count(isvesselface(tvf)) == 2) then
                    ! Get the other two vertices
                    allocate(tvfv(2, 2))
                    tvfv = grid%face%vert(pack(tvf, isvesselface(tvf)), :)
                    allocate(tvn(count(tvfv /= tv)))
                    tvn = pack(tvfv, tvfv /= tv)

                    ! Check
                    if (size(tvn, 1) /= 2) then 
                        ! Shouldn't happen, throw error
                        call gdErrorHandler('InitializeCostFunctionLR: ' // &
                            'unknown error, something seems wrong in grid interconnection')
                    end if 

                    ! Add
                    vpairs(tv, 1:2) = tvn
                    b0(tv) = 1
                    nvpairs(tv) = 1 

                    ! Deallocate 
                    deallocate(tvfv, tvf, tvn)
                end if 

            end do

            ! Housekeeping
            deallocate(vesselvert)


        end if 

        ! Housekeeping
        deallocate(Btx, Bty)

        ! End associate
        end associate

        ! Write data
        !===========
        if (options%writedata == 1) then 
            call costfunction%WriteData(grid)
        end if 

    end subroutine

    ! Cost function data writing 
    subroutine WriteCostFunctionDataLRrad(costfunction, grid)

        ! Description
        !============
        ! Write out the cost function data for the LR cost function.
        ! Here, this consists of the vertex pair data in IDn, xn, yn 
        ! format

        ! Declare variables
        !==================
        ! Arguments
        class(CostFunctionLRradUDT)     :: costfunction 
        type(GridUDT)                   :: grid

        ! Auxiliary
        integer(I8)                     :: ncol, nrow 

        integer(I8), allocatable        :: IDn(:, :) 
        real(R8), allocatable           :: xn(:, :), yn(:, :)
        character(:), allocatable       :: filename 

        ! Loop
        integer(I8)                     :: j 

        ! Initialize
        !===========
        ! Set filename
        allocate(character(len('costfunction_vertexpairs_LRrad')) :: filename)
        filename = 'costfunction_vertexpairs_LRrad'

        ! Allocate
        nrow = size(costfunction%vpairs, 1)
        ncol = size(costfunction%vpairs, 2)
        allocate(IDn(nrow, ncol), xn(nrow, ncol), yn(nrow, ncol))

        ! Unpack
        associate(&
            vpairs      => costfunction%vpairs,         &
            x           => grid%vert%x,                 &
            y           => grid%vert%y)

        ! Loop
        do j = 1, ncol 
            IDn(:, j) = vpairs(:, j) 
            xn(:, j) = x(vpairs(:, j)) 
            yn(:, j) = y(vpairs(:, j)) 
        end do

        ! Call writer
        !============
        call WriteVertexPairData(IDn, xn, yn, filename)

        ! Housekeeping
        !=============
        end associate
        deallocate(IDn, xn, yn)
        


    end subroutine

    !------------------------------------------------------------------!
    !                        FACE ANGLE DIFFERENCE                     !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCostFunctionFAD(costfunction, grid, &
        magneticField, environment, options)

        ! Description
        !============
        ! Initialize the cost function and its parameters based on the 
        ! grid, magnetic field, and environment structures. To determine 
        ! which faces should be considered in the FAD cost function, the
        ! vertex ID is compared with the ID of its neighbours. If it is 
        ! the same (and non-zero), the face is NOT considered as it is 
        ! aligned. If it is not the same and both are non-zero, it is 
        ! a potential face to consider. 

        ! Modules
        !========

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionFADUDT)           :: costfunction
        type(GridUDT)                       :: grid
        type(MagneticFieldUDT)              :: magneticField 
        type(EnvironmentUDT)                :: environment 
        type(CostFunctionOptionsUDT)        :: options
        
        ! Loop variables
        integer(I8)                         :: i, j, k, vpc

        ! Auxiliary variables
        type(StructuredPLF2DDistanceDFUDT)  :: dfwt

        integer(I8)                         :: tID, vp1(1:2), &
            vp2(1:2), ntvn, ntemptvn, nvp

        real(R8)                            :: dxf, dyf, txf, tyf, tnf,&
            bxf, byf, bnf, dpf

        logical                             :: thischeck

        integer(I8), allocatable            :: tvn(:), temptvn(:), &
            vpairs(:,:), nvpairs(:), reverse(:), tempvpairs(:,:)

        real(R8), allocatable               :: bx(:), by(:), xv(:,:), &
            yv(:,:), xf(:,:), yf(:,:), gxf(:,:), gyf(:,:), dx(:,:), &
            dy(:,:), dotprod(:,:), wt(:)

        logical, allocatable                :: cID(:), mask(:), &
            isaligned(:)

        ! Data
        real(R8)                            :: epsalignment 
        data epsalignment /0.1/
        
        ! Initialize
        !===========
        ! Set the scaling constant
        costfunction%lambda = options%FAD%lambda

        ! Allocate (initialize too big)
        nvp = grid%vert%ntot*maxval(grid%vert%neigP(:,2)) ! maximal number of pairs
        allocate(vpairs(nvp,4), nvpairs(grid%vert%ntot))
        ! call costfunction%Allocate(grid%vert%ntot, 4)

        ! Compute magnetic field in grid points
        allocate(bx(grid%vert%ntot), by(grid%vert%ntot))
        call magneticField%interp%Evaluate(grid%vert%x, grid%vert%y, 0, 1, bx)
        call magneticField%interp%Evaluate(grid%vert%x, grid%vert%y, 1, 0, by)
        !call EvaluateBicubicSplineInterpolant( &
        !    grid%vert%x, grid%vert%y, bx, magneticField%interp, '0', '1')
        !call EvaluateBicubicSplineInterpolant( &
        !    grid%vert%x, grid%vert%y, by, magneticField%interp, '1', '0')
        bx = -bx

        ! Associate some fields
        associate(&
            vert    => grid%vert, &
            x       => grid%vert%x, &
            y       => grid%vert%y &
            )

        ! Determine vertex pairs
        !=======================
        ! Initialize
        vpairs(:, :) = 0
        nvpairs(:) = 0
        vpc = 1 ! vertex pair index & counter

        ! Loop over all vertices
        do i = 1, vert%ntot
            ! Skip if the current vertex is a boundary vertex 
            if (vert%BV(i)) then 
                cycle ! skip the remainder of the do-loop for this iteration
            end if

            ! Dimensions
            ntemptvn = vert%neigP(i, 2)

            ! Allocate
            allocate(temptvn(ntemptvn))
            allocate(cID(ntemptvn))

            ! Get the vertex neighbours of this vertex
            temptvn = vert%neig(& 
                vert%neigP(i, 1):vert%neigP(i, 1)+vert%neigP(i, 2)-1)
            
            ! Get the ID of the coordinate line
            tID = vert%fieldlineID(i)

            ! Check if there is an ID
            ntvn = 0
            if (tID == 0) then ! no vertex ID - check all vertex neighbours
                ! Check which faces are aligned to know which vertex 
                ! pairs to (not) include. Each vertex pair consists of 
                ! the current node and one of its neighbours in temptvn. 
                ! Alignment checking is based on the dot product between 
                ! the normalized face and magnetic field vectors. The 
                ! magnetic field vector is interpolated from the 
                ! precomputed magnetic field on the vertices. 

                ! Allocate & initialize aligned indicator
                allocate(isaligned(ntemptvn))
                isaligned(:) = .false. 

                ! Loop
                do j = 1, ntemptvn
                    ! Get normalized face vector
                    dxf = x(temptvn(j)) - x(i)
                    dyf = y(temptvn(j)) - y(i)
                    tnf = sqrt(dxf**2 + dyf**2)
                    txf = dxf/tnf 
                    tyf = dyf/tnf 

                    ! Get normalized magnetic field vector, approximated
                    ! by simple arithmetic average
                    bxf = (bx(temptvn(j)) + bx(i))*0.5
                    byf = (by(temptvn(j)) + by(i))*0.5
                    bnf = sqrt(bxf**2 + byf**2)
                    bxf = bxf/bnf 
                    byf = byf/bnf

                    ! Compute dot product
                    dpf = bxf*txf + byf*tyf 

                    ! Check
                    if (abs( abs(dpf) - 1 ) < epsalignment) then 
                        ! Consider as aligned 
                        isaligned(j) = .true.
                    end if
                end do

                ! Extract the non-aligned nodes if two non-aligned faces
                ! remain. Otherwise, set to zero. 
                ntvn = count(.not. isaligned)
                if (ntvn .ne. 2) then 
                    ntvn = 0
                    ! Allocate - to prevent deallocation errors downstream
                    allocate(tvn(ntvn))
                else
                    ! Allocate
                    allocate(tvn(ntvn))

                    ! Extract the vertices
                    tvn = pack(temptvn, .not. isaligned)
                end if
    
                ! Deallocate
                deallocate(isaligned)

            else
                ! Check which vertices have the same ID
                cID = tID == vert%fieldlineID(temptvn)

                ! Check which faces are aligned
                allocate(isaligned(ntemptvn))

                ! Loop
                do j = 1, ntemptvn
                    ! Get normalized face vector
                    dxf = x(temptvn(j)) - x(i)
                    dyf = y(temptvn(j)) - y(i)
                    tnf = sqrt(dxf**2 + dyf**2)
                    txf = dxf/tnf 
                    tyf = dyf/tnf 

                    ! Get normalized magnetic field vector, approximated
                    ! by simple arithmetic average
                    bxf = (bx(temptvn(j)) + bx(i))*0.5
                    byf = (by(temptvn(j)) + by(i))*0.5
                    bnf = sqrt(bxf**2 + byf**2)
                    bxf = bxf/bnf 
                    byf = byf/bnf

                    ! Compute dot product
                    dpf = bxf*txf + byf*tyf 

                    ! Check
                    if (abs( abs(dpf) - 1 ) < epsalignment) then 
                        ! Consider as aligned 
                        isaligned(j) = .true.
                    end if
                end do
            
                ! Extract vertices that have NOT the same ID and are not
                ! aligned. 
                allocate(mask(ntemptvn))
                mask = (.not. cID) ! .and. (.not. isaligned)
                allocate(tvn(count(mask)))
                tvn = pack(temptvn, mask)

                ! Only keep if there are two vertices 
                if  (size(tvn) .ne. 2) then
                    ntvn = 0
                else
                    ntvn = 2
                end if


                ! Deallocate
                deallocate(mask, isaligned)
            end if
            
            ! Constrain each pair (more generally written in case 
            ! multiple pairs are allowed in the future)
            nvpairs(i) = (ntvn/2) ! this automatically floors
            do j = 1, nvpairs(i)

                ! Get the vertex pairs
                vp1 = [i, tvn(2*j-1)]
                vp2 = [i, tvn(2*j)]
                
                ! Check if these pairs should be added. For now, we 
                ! don't consider vertices where there are more than two
                ! pairs, as here the cost function does not make much
                ! sense. 
                thischeck = .true.
                thischeck = thischeck .and. (.not. (nvpairs(i) > 1)) 
                
                ! Add the pair if allowed
                if (thischeck) then 
                    ! Add the pairs
                    vpairs(vpc,:) = [vp1, vp2]
                    
                    ! Update counter
                    vpc = vpc+1

                end if
            end do

            ! Deallocate
            deallocate(tvn, temptvn, cID)
        end do

        ! Check order
        !============
        ! Check the vertex pair order by computing the dot product with 
        ! the local magnetic field (now evaluated with the actual face
        ! coordinates)

        ! Update counter (-1 to get actual number of edges)
        vpc = vpc-1

        ! Allocate
        allocate(xv(vpc, 4), yv(vpc, 4), xf(vpc, 2), &
            yf(vpc, 2), gxf(vpc, 2), gyf(vpc, 2), &
            dx(vpc, 2), dy(vpc, 2), dotprod(vpc, 2), &
            tempvpairs(size(vpairs,1), size(vpairs, 2)))
        call costfunction%Allocate(vpc, 4)

        ! Compute vectors
        do i = 1, 4
            xv(:, i) = x(vpairs(1:vpc, i)) !reshape(x(vpairs),[vert%ntot,4])
            yv(:, i) = y(vpairs(1:vpc, i)) !reshape(y(vpairs),[vert%ntot,4])
        end do

        ! Compute edge faces
        xf(:, 1) = xv(:, 1) + xv(:, 2)
        xf(:, 2) = xv(:, 3) + xv(:, 4)
        xf = 0.5*xf
        yf(:, 1) = yv(:, 1) + yv(:, 2)
        yf(:, 2) = yv(:, 3) + yv(:, 4)
        yf = 0.5*yf

        ! Compute the vector perpendicular on the magnetic field vector
        do i = 1, 2
            call magneticField%interp%Evaluate(xf(:, i), yf(:, i ), 1, 0, gxf(:, i))
            call magneticField%interp%Evaluate(xf(:, i), yf(:, i ), 0, 1, gyf(:, i))
            !call EvaluateBicubicSplineInterpolant( &
            !    xf(:,i), yf(:,i), gxf(:,i), magneticField%interp, '1', '0')
            !call EvaluateBicubicSplineInterpolant( &
            !    xf(:,i), yf(:,i), gyf(:,i), magneticField%interp, '0', '1')
        end do

        ! Compute distances
        dx(:, 1) = xv(:, 2) - xv(:, 1)
        dx(:, 2) = xv(:, 4) - xv(:, 3)
        dy(:, 1) = yv(:, 2) - yv(:, 1)
        dy(:, 2) = yv(:, 4) - yv(:, 3)

        ! Compute dot product
        dotprod = dx*gxf + dy*gyf

        ! Compute faces to reverse 
        tempvpairs = vpairs
        allocate(mask(vpc))
        do i = 1, 2
            mask = (dotprod(:,i) < 0)
            allocate(reverse(count(mask)))
            reverse = pack([(k, k = 1, vpc)], mask)
            vpairs(reverse, 2*i-1)  = tempvpairs(reverse, 2*i)
            vpairs(reverse, 2*i)    = tempvpairs(reverse, 2*i-1)
            deallocate(reverse)
        end do

        ! Determine weigths
        !==================
        ! Initialize (magnetic field as dummy since unsigned anyway)
        call dfwt%Initialize(magneticField%interp, &
            environment%vessel%plfvessel, environment%vessel%plfvessel, &
            options%FAD%weightatvessel, options%FAD%weightatinf, &
            options%FAD%decaylength, 'unsigned')

        ! Evaluate
        allocate(wt(size(xv, 1)))
        call dfwt%Evaluate(xv(:, 1), yv(:, 1), wt)

        ! Visualize
        call dfwt%Visualize([minval(x), maxval(x)], &
            [minval(y), maxval(y)], 100, 100, 'costfunctionFAD_weights')

        ! Assign to cost function
        !========================
        costfunction%vpairs = vpairs(1:vpc,:)
        costfunction%nvpairs = vpc 
        costfunction%wt = wt

        ! Visualize? 
        !allocate(xplot(size(xf)), yplot(size(yf)))
        !xplot = reshape(xf, [size(xf)])
        !yplot = reshape(yf, [size(yf)])
        !call PlotGridWithPoints(grid, xplot, yplot, '-p')
        !deallocate(xplot, yplot)
        

        ! End associate
        end associate

        ! Write data
        !===========
        if (options%writedata == 1) then 
            call costfunction%WriteData(grid) 
        end if 

        ! Deallocate
        deallocate(tempvpairs, dx, dy, gxf, gyf, xv, yv, xf, &
            yf, dotprod, bx, by, vpairs, nvpairs, mask)
        

        

    end subroutine

    ! Cost function evaluation
    subroutine EvaluateCostFunctionFAD(costfunction, J, gradJ, hessJ, &
        grid, magneticField, environment, dogradient, dohessian, &
        designvariables)

        ! Description
        !============
        ! Evaluate the cost function, the gradient and its hessian. 
        ! Here, we simply call the same cost function twice, but switch
        ! the order of the indices and recompute the bias. 

        ! Notes:
        !=======
        ! Possible future performance improvements:
        ! - Allocating hessian stuff only once and storing indices, 
        ! since they don't change
        ! - Instead of recomputing auxiliary variables, store them. May
        ! not actually be better in terms of computational time, but 
        ! may lead to shorter and hence better maintainable code. 

        ! Modules
        !========

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionFADUDT)       :: costfunction 
        real(R8)                        :: J
        real(R8), allocatable           :: gradJ(:)
        type(MySparseUDT)               :: hessJ
        type(GridUDT)                   :: grid 
        type(MagneticFieldUDT)          :: magneticField 
        type(EnvironmentUDT)            :: environment
        logical                         :: dogradient, dohessian 
        class(DesignVariablesGDUDT)     :: designvariables

        ! Loop variables
        integer(I8)                     :: i, cc

        ! Auxiliary
        real(R8)                        :: gxf1, gxf2, gyf1, gyf2, dx1, &
            dx2, dy1, dy2, sp1, sp2, cp1, cp2, rat1, rat2, gxxf1, &
            gxxf2, gxyf1, gxyf2, gyxf1, gyxf2, gyyf1, gyyf2, gxxxf1, &
            gxyxf1, gxyyf1, gyxxf1, gyyxf1, gyyyf1, gxxxf2, &
            gxyxf2, gxyyf2, gyxxf2, gyyxf2, gyyyf2, wt

        integer(I8), allocatable        :: row(:), col(:)

        real(R8), allocatable           :: xv(:,:), yv(:,:), xfv1(:), &
            yfv1(:), xfv2(:), yfv2(:), gxfv1(:), gyfv1(:), gxfv2(:), &
            gyfv2(:), dxv1(:), dyv1(:), dxv2(:), dyv2(:), spv1(:), &
            spv2(:), cpv1(:), cpv2(:), ratv1(:), ratv2(:), gxxfv1(:), &
            gxxfv2(:), gxyfv1(:), gxyfv2(:), gyxfv1(:), gyxfv2(:), &
            gyyfv1(:), gyyfv2(:), gxxxfv1(:), gxxxfv2(:), gxxyfv1(:), &
            gxxyfv2(:), gxyxfv1(:), gxyxfv2(:), gxyyfv1(:), gxyyfv2(:), &
            gyxxfv1(:), gyxxfv2(:), gyxyfv1(:), gyxyfv2(:), gyyyfv1(:), &
            gyyyfv2(:), valxx(:), valxy(:), valyx(:), valyy(:) 

        ! Initialize
        !===========
        ! Cost function
        J = 0

        ! Gradient
        gradJ(:) = 0

        ! Associate
        associate(&
            np          => costfunction%nvpairs,    &
            vpairs      => costfunction%vpairs,     &
            lambda      => costfunction%lambda,     &
            wtv         => costfunction%wt,         &
            x           => grid%vert%x,             &
            y           => grid%vert%y,             &
            vert        => grid%vert)

        ! Allocate auxiliary arrays
        allocate(xv(np, 4), yv(np, 4), xfv1(np), xfv2(np), yfv1(np), &
            yfv2(np), gxfv1(np), gxfv2(np), gyfv1(np), gyfv2(np), &
            dxv1(np), dxv2(np), dyv1(np), dyv2(np), &
            spv1(np), spv2(np), cpv1(np), cpv2(np), ratv1(np), &
            ratv2(np))

        ! Coordinates
        do i = 1, 4
            xv(:, i) = x(vpairs(:, i))
            yv(:, i) = y(vpairs(:, i))
        end do

        ! Face coordinates
        xfv1 = 0.5*( xv(:, 1) + xv(:, 2) )  
        xfv2 = 0.5*( xv(:, 3) + xv(:, 4) ) 
        yfv1 = 0.5*( yv(:, 1) + yv(:, 2) ) 
        yfv2 = 0.5*( yv(:, 3) + yv(:, 4) ) 

        ! Coordinate line vectors at face coordinates
        call magneticField%interp%Evaluate(xfv1, yfv1, 1, 0, gxfv1)
        call magneticField%interp%Evaluate(xfv2, yfv2, 1, 0, gxfv2)
        call magneticField%interp%Evaluate(xfv1, yfv1, 0, 1, gyfv1)
        call magneticField%interp%Evaluate(xfv2, yfv2, 0, 1, gyfv2)
        
        ! Face vectors
        dxv1 = xv(:, 2) - xv(:, 1)
        dxv2 = xv(:, 4) - xv(:, 3)
        dyv1 = yv(:, 2) - yv(:, 1)
        dyv2 = yv(:, 4) - yv(:, 3)

        spv1 = gxfv1*dxv1 + gyfv1*dyv1
        cpv1 = gxfv1*dyv1 - gyfv1*dxv1
        spv2 = gxfv2*dxv2 + gyfv2*dyv2
        cpv2 = gxfv2*dyv2 - gyfv2*dxv2

        ratv1 = cpv1/spv1
        ratv2 = cpv2/spv2

        ! Compute cost function
        !======================
        ! Sum contributions
        J = sum( 0.5*( wtv*(atan(ratv1) - atan(ratv2) )**2 ) )

        ! Scale
        J = lambda*J

        ! Precompute values
        !==================
        ! Precompute values needed for gradient and hessian calculation,
        ! based on the design variable type. 
        if (dohessian) then

            select case (trim(designvariables%type))

            case ('coordinates', 'coordinates_desiredflux')

                ! Allocate
                allocate(gxxfv1(np), gxxfv2(np), gxyfv1(np), gxyfv2(np), &
                    gyxfv1(np), gyxfv2(np), gyyfv1(np), gyyfv2(np))
                allocate(gxxxfv1(np), gxxxfv2(np), gxxyfv1(np), &
                    gxxyfv2(np), gxyxfv1(np), gxyxfv2(np), &
                    gxyyfv1(np), gxyyfv2(np), gyxxfv1(np), gyxxfv2(np), &
                    gyxyfv1(np), gyxyfv2(np), gyyyfv1(np), gyyyfv2(np))

                ! Precompute gradient quantities
                call magneticField%interp%Evaluate(xfv1, yfv1, 2, 0, gxxfv1)
                call magneticField%interp%Evaluate(xfv1, yfv1, 1, 1, gxyfv1)
                call magneticField%interp%Evaluate(xfv1, yfv1, 0, 2, gyyfv1)
                gyxfv1 = gxyfv1 ! symmetric, done for ease here

                call magneticField%interp%Evaluate(xfv2, yfv2, 2, 0, gxxfv2)
                call magneticField%interp%Evaluate(xfv2, yfv2, 1, 1, gxyfv2)
                call magneticField%interp%Evaluate(xfv2, yfv2, 0, 2, gyyfv2)
                gyxfv2 = gxyfv2 ! symmetric, done for ease here

                ! Precompute hessian quantities
                call magneticField%interp%Evaluate(xfv1, yfv1, 3, 0, gxxxfv1)
                call magneticField%interp%Evaluate(xfv1, yfv1, 2, 1, gxxyfv1)
                call magneticField%interp%Evaluate(xfv1, yfv1, 1, 2, gxyyfv1)
                call magneticField%interp%Evaluate(xfv1, yfv1, 0, 3, gyyyfv1)
                gyxyfv1 = gxyyfv1 ! symmetric, repeated for ease
                gyxxfv1 = gxxyfv1 

                call magneticField%interp%Evaluate(xfv2, yfv2, 3, 0, gxxxfv2)
                call magneticField%interp%Evaluate(xfv2, yfv2, 2, 1, gxxyfv2)
                call magneticField%interp%Evaluate(xfv2, yfv2, 1, 2, gxyyfv2)
                call magneticField%interp%Evaluate(xfv2, yfv2, 0, 3, gyyyfv2)
                gyxyfv2 = gxyyfv2 ! symmetric, repeated for ease
                gyxxfv2 = gxxyfv2 

            case default

                ! Not implemented, throw error
                call gdErrorHandler('EvaluateCostFunctionFAD: gradient' &
                    // ' not yet implemented for this design variable' &
                    // ' type')

            end select 

        elseif (dogradient) then

            select case (trim(designvariables%type))

            case ('coordinates', 'coordinates_desiredflux')

                ! Allocate
                allocate(gxxfv1(np), gxxfv2(np), gxyfv1(np), gxyfv2(np), &
                    gyxfv1(np), gyxfv2(np), gyyfv1(np), gyyfv2(np))

                ! Precompute gradient quantities
                call magneticField%interp%Evaluate(xfv1, yfv1, 2, 0, gxxfv1)
                call magneticField%interp%Evaluate(xfv1, yfv1, 1, 1, gxyfv1)
                call magneticField%interp%Evaluate(xfv1, yfv1, 0, 2, gyyfv1)
                gyxfv1 = gxyfv1 ! symmetric, done for ease here

                call magneticField%interp%Evaluate(xfv2, yfv2, 2, 0, gxxfv2)
                call magneticField%interp%Evaluate(xfv2, yfv2, 1, 1, gxyfv2)
                call magneticField%interp%Evaluate(xfv2, yfv2, 0, 2, gyyfv2)
                gyxfv2 = gxyfv2 ! symmetric, done for ease here

            end select

        end if

        ! Compute gradient
        !=================
        if (dogradient) then 

            ! Check the design variables
            select case (trim(designvariables%type))

            case ('coordinates', 'coordinates_desiredflux')

                ! Loop over all pairs
                do i = 1, np
                    ! Unpack
                    wt = wtv(i) 
                    
                    ! Coordinate line vectors at face coordinates
                    gxf1 = gxfv1(i) 
                    gyf1 = gyfv1(i) 
                    gxf2 = gxfv2(i) 
                    gyf2 = gyfv2(i) 
                    
                    ! Face vectors
                    dx1 = dxv1(i) 
                    dx2 = dxv2(i) 
                    dy1 = dyv1(i) 
                    dy2 = dyv2(i) 
                    
                    sp1 = spv1(i) 
                    cp1 = cpv1(i) 
                    sp2 = spv2(i) 
                    cp2 = cpv2(i) 
                    
                    rat1 = ratv1(i) 
                    rat2 = ratv2(i) 
                    
                    ! Derivatives
                    gxxf1 = gxxfv1(i) 
                    gxyf1 = gxyfv1(i) 
                    gyxf1 = gyxfv1(i) 
                    gyyf1 = gyyfv1(i) 
                    
                    gxxf2 = gxxfv2(i) 
                    gxyf2 = gxyfv2(i) 
                    gyxf2 = gyxfv2(i) 
                    gyyf2 = gyyfv2(i) 
                    
                    ! Evaluate gradient
                    gradJ(vpairs(i,1)) = gradJ(vpairs(i,1)) + &
                        (wt*(atan(rat1) - atan(rat2))* &
                        ((gyf1 - 0.5*dx1*gyxf1 + 0.5*dy1*gxxf1)/sp1 - &
                        (rat1*(0.5*dx1*gxxf1 - gxf1 + 0.5*dy1*gyxf1))/sp1))/(rat1**2 + 1) 
                    gradJ(vpairs(i,2)) = gradJ(vpairs(i,2)) + &
                        -(wt*(atan(rat1) - atan(rat2))* &
                        ((gyf1 + 0.5*dx1*gyxf1 - 0.5*dy1*gxxf1)/sp1 + &
                        (rat1*(gxf1 + 0.5*dx1*gxxf1 + 0.5*dy1*gyxf1))/sp1))/(rat1**2 + 1) 
                    gradJ(vpairs(i,3)) = gradJ(vpairs(i,3)) + &
                        -(wt*(atan(rat1) - atan(rat2))*&
                        ((gyf2 - 0.5*dx2*gyxf2 + 0.5*dy2*gxxf2)/sp2 - &
                        (rat2*(0.5*dx2*gxxf2 - gxf2 + 0.5*dy2*gyxf2))/sp2))/(rat2**2 + 1) 
                    gradJ(vpairs(i,4)) = gradJ(vpairs(i,4)) + &
                        (wt*(atan(rat1) - atan(rat2))* &
                        ((gyf2 + 0.5*dx2*gyxf2 - 0.5*dy2*gxxf2)/sp2 + &
                        (rat2*(gxf2 + 0.5*dx2*gxxf2 + 0.5*dy2*gyxf2))/sp2))/(rat2**2 + 1) 
                        
                    
                    gradJ(vpairs(i,1)+size(x)) = gradJ(vpairs(i,1)+size(x)) + &
                        -(wt*(atan(rat1) - atan(rat2))* &
                        ((gxf1 + 0.5*dx1*gyyf1 - 0.5*dy1*gxyf1)/sp1 + &
                        (rat1*(0.5*dx1*gxyf1 - gyf1 + 0.5*dy1*gyyf1))/sp1))/(rat1**2 + 1) 
                    gradJ(vpairs(i,2)+size(x)) = gradJ(vpairs(i,2)+size(x)) + &
                        (wt*(atan(rat1) - atan(rat2))* &
                        ((gxf1 - 0.5*dx1*gyyf1 + 0.5*dy1*gxyf1)/sp1 - & 
                        (rat1*(gyf1 + 0.5*dx1*gxyf1 + 0.5*dy1*gyyf1))/sp1))/(rat1**2 + 1) 
                    gradJ(vpairs(i,3)+size(x)) = gradJ(vpairs(i,3)+size(x)) + &
                        (wt*(atan(rat1) - atan(rat2))* &
                        ((gxf2 + 0.5*dx2*gyyf2 - 0.5*dy2*gxyf2)/sp2 + &
                        (rat2*(0.5*dx2*gxyf2 - gyf2 + 0.5*dy2*gyyf2))/sp2))/(rat2**2 + 1) 
                    gradJ(vpairs(i,4)+size(x)) = gradJ(vpairs(i,4)+size(x)) + &
                        -(wt*(atan(rat1) - atan(rat2))* &
                        ((gxf2 - 0.5*dx2*gyyf2 + 0.5*dy2*gxyf2)/sp2 - &
                        (rat2*(gyf2 + 0.5*dx2*gxyf2 + 0.5*dy2*gyyf2))/sp2))/(rat2**2 + 1) 
                    
                end do

            case default

                ! Not implemented, throw error
                call gdErrorHandler('EvaluateCostFunctionFAD: gradient' &
                    // ' not yet implemented for this design variable' &
                    // ' type')

            end select

            ! Scale
            gradJ = lambda*gradJ

        end if

        ! Compute hessian
        !================
        ! Allocate the hessian (if not already done so)
        if (.not. allocated(hessJ%row)) then
            ! Allocate the sparse matrix
            select case (trim(designvariables%type))

            case ('coordinates', 'coordinates_desiredflux')

                hessJ%nval = 64*np ! this should be exact and constant

            end select
            call hessJ%Allocate()

        end if

        if (dohessian) then

            ! Check the design variables
            select case (trim(designvariables%type))

            case ('coordinates', 'coordinates_desiredflux')

                ! Allocate auxiliary variables
                allocate(row(16*np))
                allocate(col(16*np))
                allocate(valxx(16*np))
                allocate(valxy(16*np))
                allocate(valyx(16*np))
                allocate(valyy(16*np))

                ! Initialize counter
                cc = 1

                ! Loop over all pairs
                do i = 1, np
                    ! Unpack
                    wt = wtv(i) 
                    
                    ! Coordinate line vectors at face coordinates
                    gxf1 = gxfv1(i) 
                    gyf1 = gyfv1(i) 
                    gxf2 = gxfv2(i) 
                    gyf2 = gyfv2(i) 
                    
                    ! Face vectors
                    dx1 = dxv1(i) 
                    dx2 = dxv2(i) 
                    dy1 = dyv1(i) 
                    dy2 = dyv2(i) 
                    
                    sp1 = spv1(i) 
                    cp1 = cpv1(i) 
                    sp2 = spv2(i) 
                    cp2 = cpv2(i) 
                    
                    rat1 = ratv1(i) 
                    rat2 = ratv2(i) 
                    
                    ! Derivatives
                    gxxf1 = gxxfv1(i) 
                    gxyf1 = gxyfv1(i) 
                    gyxf1 = gyxfv1(i) 
                    gyyf1 = gyyfv1(i) 
                    
                    gxxf2 = gxxfv2(i) 
                    gxyf2 = gxyfv2(i) 
                    gyxf2 = gyxfv2(i) 
                    gyyf2 = gyyfv2(i) 

                    gxxxf1 = gxxxfv1(i)  
                    gxyxf1 = gxxyfv1(i)  
                    gxyyf1 = gxyyfv1(i)  
                    gyxxf1 = gyxxfv1(i)  
                    gyyxf1 = gyxyfv1(i)  
                    gyyyf1 = gyyyfv1(i) 
                    
                    gxxxf2 = gxxxfv2(i)  
                    gxyxf2 = gxxyfv2(i)  
                    gxyyf2 = gxyyfv2(i)  
                    gyxxf2 = gyxxfv2(i)  
                    gyyxf2 = gyxyfv2(i)  
                    gyyyf2 = gyyyfv2(i) 

                    ! Compute contributions
                    row(cc) = vpairs(i,1)  
                    col(cc) = vpairs(i,1) 
                    valxx(cc) = wt*(((1.0*gyxf1 - 0.25*dx1*gyxxf1 + &
                        0.25*dy1*gxxxf1)/sp1 + &
                        (2*rat1*(0.5*dx1*gxxf1 - &
                        gxf1 + 0.5*dy1*gyxf1)**2)/sp1**2 - &
                        (2*(gyf1 - 0.5*dx1*gyxf1 + 0.5*dy1*gxxf1)*&
                        (0.5*dx1*gxxf1 - gxf1 + 0.5*dy1*gyxf1))&
                        /sp1**2 - (rat1*(0.25*dx1*gxxxf1 - 1.0*gxxf1 + &
                        0.25*dy1*gyxxf1))/sp1)/(rat1**2 + 1) + &
                        (((2*rat1**2*(0.5*dx1*gxxf1 - gxf1 + &
                        0.5*dy1*gyxf1))/sp1 - (2*rat1*(gyf1 - &
                        0.5*dx1*gyxf1 + 0.5*dy1*gxxf1))/sp1)*&
                        ((gyf1 - 0.5*dx1*gyxf1 + 0.5*dy1*gxxf1)&
                        /sp1 - (rat1*(0.5*dx1*gxxf1 - gxf1 + &
                        0.5*dy1*gyxf1))/sp1))/(rat1**2 + 1)**2)*&
                        (atan(rat1) - atan(rat2)) + &
                        (wt*((gyf1 - 0.5*dx1*gyxf1 + 0.5*dy1*gxxf1)/sp1 - &
                        (rat1*(0.5*dx1*gxxf1 - gxf1 + 0.5*dy1*gyxf1))&
                        /sp1)**2)/(rat1**2 + 1)**2 ! x1x1
                    valxy(cc) = wt*(atan(rat1) - atan(rat2))*&
                        ((((gxf1 + 0.5*dx1*gyyf1 - 0.5*dy1*gxyf1)*&
                        (0.5*dx1*gxxf1 - gxf1 + 0.5*dy1*gyxf1))/sp1**2 - &
                        (0.5*gxxf1 - 0.5*gyyf1 + 0.25*dx1*gyyxf1 - &
                        0.25*dy1*gxyxf1)/sp1 - ((gyf1 - 0.5*dx1*gyxf1 + &
                        0.5*dy1*gxxf1)*(0.5*dx1*gxyf1 - gyf1 + 0.5*dy1*gyyf1))&
                        /sp1**2 + (rat1*(0.5*gxyf1 + 0.5*gyxf1 - &
                        0.25*dx1*gxyxf1 - 0.25*dy1*gyyxf1))/sp1 + &
                        (2*rat1*(0.5*dx1*gxxf1 - gxf1 + 0.5*dy1*gyxf1)*&
                        (0.5*dx1*gxyf1 - gyf1 + 0.5*dy1*gyyf1))/sp1**2)&
                        /(rat1**2 + 1) + (((2*rat1**2*(0.5*dx1*gxyf1 - &
                        gyf1 + 0.5*dy1*gyyf1))/sp1 + (2*rat1*(gxf1 + &
                        0.5*dx1*gyyf1 - 0.5*dy1*gxyf1))/sp1)*&
                        ((gyf1 - 0.5*dx1*gyxf1 + 0.5*dy1*gxxf1)/sp1 - &
                        (rat1*(0.5*dx1*gxxf1 - gxf1 + 0.5*dy1*gyxf1))&
                        /sp1))/(rat1**2 + 1)**2) - (wt*((gyf1 - &
                        0.5*dx1*gyxf1 + 0.5*dy1*gxxf1)/sp1 - &
                        (rat1*(0.5*dx1*gxxf1 - gxf1 + 0.5*dy1*gyxf1))&
                        /sp1)*((gxf1 + 0.5*dx1*gyyf1 - 0.5*dy1*gxyf1)&
                        /sp1 + (rat1*(0.5*dx1*gxyf1 - gyf1 + 0.5*dy1*&
                        gyyf1))/sp1))/(rat1**2 + 1)**2  !x1y1
                    valyx(cc) = valxy(cc)  !y1x1
                    valyy(cc) = (wt*((gxf1 + 0.5*dx1*gyyf1 - &
                    0.5*dy1*gxyf1)/sp1 + (rat1*(0.5*dx1*gxyf1 - &
                    gyf1 + 0.5*dy1*gyyf1))/sp1)**2)/(rat1**2 + 1)**2 -&
                        wt*(((1.0*gxyf1 + 0.25*dx1*gyyyf1 - &
                        0.25*dy1*gxyyf1)/sp1 - (2*rat1*(0.5*dx1*gxyf1 - &
                        gyf1 + 0.5*dy1*gyyf1)**2)/sp1**2 - &
                        (2*(gxf1 + 0.5*dx1*gyyf1 - 0.5*dy1*gxyf1)*&
                        (0.5*dx1*gxyf1 - gyf1 + 0.5*dy1*gyyf1))/&
                        sp1**2 + (rat1*(0.25*dx1*gxyyf1 - 1.0*gyyf1 + &
                        0.25*dy1*gyyyf1))/sp1)/(rat1**2 + 1) + &
                        (((2*rat1**2*(0.5*dx1*gxyf1 - gyf1 + &
                        0.5*dy1*gyyf1))/sp1 + (2*rat1*(gxf1 + &
                        0.5*dx1*gyyf1 - 0.5*dy1*gxyf1))/sp1)*&
                        ((gxf1 + 0.5*dx1*gyyf1 - 0.5*dy1*gxyf1)/sp1 +&
                        (rat1*(0.5*dx1*gxyf1 - gyf1 + 0.5*dy1*gyyf1))&
                        /sp1))/(rat1**2 + 1)**2)*(atan(rat1) - atan(rat2))  !y1y1
                    cc = cc+1 
                    
                    row(cc) = vpairs(i,1); col(cc) = vpairs(i,2) 
                    valxx(cc) = - wt*(atan(rat1) - atan(rat2))*&
                        (((0.25*dx1*gyxxf1 - 0.25*dy1*gxxxf1)/sp1 + &
                        ((gxf1 + 0.5*dx1*gxxf1 + 0.5*dy1*gyxf1)*&
                        (gyf1 - 0.5*dx1*gyxf1 + 0.5*dy1*gxxf1))&
                        /sp1**2 - ((gyf1 + 0.5*dx1*gyxf1 - 0.5*&
                        dy1*gxxf1)*(0.5*dx1*gxxf1 - gxf1 + 0.5*dy1*&
                        gyxf1))/sp1**2 + (rat1*(0.25*dx1*gxxxf1 + 0.25*&
                        dy1*gyxxf1))/sp1 - (2*rat1*(gxf1 + 0.5*dx1*gxxf1 &
                        + 0.5*dy1*gyxf1)*(0.5*dx1*gxxf1 - gxf1 + 0.5*dy1 &
                        *gyxf1))/sp1**2)/(rat1**2 + 1) - (((2*rat1*(gyf1 &
                        + 0.5*dx1*gyxf1 - 0.5*dy1*gxxf1))/sp1 + (2*rat1**2*&
                        (gxf1 + 0.5*dx1*gxxf1 + 0.5*dy1*gyxf1))/sp1)*&
                        ((gyf1 - 0.5*dx1*gyxf1 + 0.5*dy1*gxxf1)/sp1 - &
                        (rat1*(0.5*dx1*gxxf1 - gxf1 + 0.5*dy1*gyxf1))&
                        /sp1))/(rat1**2 + 1)**2) - (wt*((gyf1 + 0.5*dx1&
                        *gyxf1 - 0.5*dy1*gxxf1)/sp1 + (rat1*(gxf1 + &
                        0.5*dx1*gxxf1 + 0.5*dy1*gyxf1))/sp1)*&
                        ((gyf1 - 0.5*dx1*gyxf1 + 0.5*dy1*gxxf1)/sp1 -&
                        (rat1*(0.5*dx1*gxxf1 - gxf1 + 0.5*dy1*gyxf1))&
                        /sp1))/(rat1**2 + 1)**2  ! x1x2
                    valxy(cc) = (wt*((gxf1 - 0.5*dx1*gyyf1 + &
                        0.5*dy1*gxyf1)/sp1 - (rat1*(gyf1 + 0.5*dx1*gxyf1 &
                        + 0.5*dy1*gyyf1))/sp1)*((gyf1 - 0.5*dx1*gyxf1 + &
                        0.5*dy1*gxxf1)/sp1 - (rat1*(0.5*dx1*gxxf1 - gxf1 &
                        + 0.5*dy1*gyxf1))/sp1))/(rat1**2 + 1)**2 - &
                        wt*(atan(rat1) - atan(rat2))*((((gyf1 - &
                        0.5*dx1*gyxf1 + 0.5*dy1*gxxf1)*(gyf1 + &
                        0.5*dx1*gxyf1 + 0.5*dy1*gyyf1))/sp1**2 - &
                        (0.5*gxxf1 + 0.5*gyyf1 - 0.25*dx1*gyyxf1 + &
                        0.25*dy1*gxyxf1)/sp1 + ((gxf1 - 0.5*dx1*gyyf1 +&
                        0.5*dy1*gxyf1)*(0.5*dx1*gxxf1 - gxf1 + &
                        0.5*dy1*gyxf1))/sp1**2 + (rat1*(0.5*gyxf1 - &
                        0.5*gxyf1 + 0.25*dx1*gxyxf1 + 0.25*dy1*gyyxf1)) &
                        /sp1 - (2*rat1*(gyf1 + 0.5*dx1*gxyf1 + &
                        0.5*dy1*gyyf1)*(0.5*dx1*gxxf1 - gxf1 + &
                        0.5*dy1*gyxf1))/sp1**2)/(rat1**2 + 1) + &
                        (((2*rat1*(gxf1 - 0.5*dx1*gyyf1 + &
                        0.5*dy1*gxyf1))/sp1 - (2*rat1**2*(gyf1 + &
                        0.5*dx1*gxyf1 + 0.5*dy1*gyyf1))/sp1)* &
                        ((gyf1 - 0.5*dx1*gyxf1 + 0.5*dy1*gxxf1)/sp1 - &
                        (rat1*(0.5*dx1*gxxf1 - gxf1 + 0.5*dy1*gyxf1)) &
                        /sp1))/(rat1**2 + 1)**2)  !x1y2
                    valyx(cc) = wt*(atan(rat1) - atan(rat2))*&
                        ((((gxf1 + 0.5*dx1*gxxf1 + 0.5*dy1*gyxf1)&
                        *(gxf1 + 0.5*dx1*gyyf1 - 0.5*dy1*gxyf1))/&
                        sp1**2 - (0.5*gxxf1 + 0.5*gyyf1 + 0.25*dx1&
                        *gyyxf1 - 0.25*dy1*gxyxf1)/sp1 + &
                        ((gyf1 + 0.5*dx1*gyxf1 - 0.5*dy1*gxxf1)*&
                        (0.5*dx1*gxyf1 - gyf1 + 0.5*dy1*gyyf1))/sp1**2 &
                         - (rat1*(0.5*gxyf1 - 0.5*gyxf1 + 0.25*dx1*&
                         gxyxf1 + 0.25*dy1*gyyxf1))/sp1 + (2*rat1*&
                         (gxf1 + 0.5*dx1*gxxf1 + 0.5*dy1*gyxf1)*&
                         (0.5*dx1*gxyf1 - gyf1 + 0.5*dy1*gyyf1))&
                         /sp1**2)/(rat1**2 + 1) - (((gyf1 + &
                         0.5*dx1*gyxf1 - 0.5*dy1*gxxf1)/sp1 + &
                         (rat1*(gxf1 + 0.5*dx1*gxxf1 + 0.5*dy1*&
                         gyxf1))/sp1)*((2*rat1**2*(0.5*dx1*gxyf1 - &
                         gyf1 + 0.5*dy1*gyyf1))/sp1 + (2*rat1*&
                         (gxf1 + 0.5*dx1*gyyf1 - 0.5*dy1*gxyf1))/sp1))&
                         /(rat1**2 + 1)**2) + (wt*((gyf1 + &
                         0.5*dx1*gyxf1 - 0.5*dy1*gxxf1)/sp1 + &
                         (rat1*(gxf1 + 0.5*dx1*gxxf1 + &
                         0.5*dy1*gyxf1))/sp1)*((gxf1 + 0.5*dx1*gyyf1 - &
                         0.5*dy1*gxyf1)/sp1 + (rat1*(0.5*dx1*gxyf1 - &
                         gyf1 + 0.5*dy1*gyyf1))/sp1))/(rat1**2 + 1)**2  !y1x2
                    valyy(cc) = - wt*(atan(rat1) - atan(rat2))*&
                        (((0.25*dx1*gyyyf1 - 0.25*dy1*gxyyf1)/sp1 - &
                        ((gxf1 + 0.5*dx1*gyyf1 - 0.5*dy1*gxyf1)*(gyf1 + &
                        0.5*dx1*gxyf1 + 0.5*dy1*gyyf1))/sp1**2 + &
                        ((gxf1 - 0.5*dx1*gyyf1 + 0.5*dy1*gxyf1)*&
                        (0.5*dx1*gxyf1 - gyf1 + 0.5*dy1*gyyf1))/sp1**2 + &
                        (rat1*(0.25*dx1*gxyyf1 + 0.25*dy1*gyyyf1))/sp1 - &
                        (2*rat1*(gyf1 + 0.5*dx1*gxyf1 + 0.5*dy1*gyyf1)*&
                        (0.5*dx1*gxyf1 - gyf1 + 0.5*dy1*gyyf1))/sp1**2)/&
                        (rat1**2 + 1) - (((gxf1 - 0.5*dx1*gyyf1 +&
                        0.5*dy1*gxyf1)/sp1 - (rat1*(gyf1 + 0.5*dx1*gxyf1 &
                        + 0.5*dy1*gyyf1))/sp1)*((2*rat1**2*(0.5*dx1*gxyf1 &
                        - gyf1 + 0.5*dy1*gyyf1))/sp1 + (2*rat1*(gxf1 + &
                        0.5*dx1*gyyf1 - 0.5*dy1*gxyf1))/sp1))/(rat1**2 + 1)**2) &
                        - (wt*((gxf1 - 0.5*dx1*gyyf1 + 0.5*dy1*gxyf1)/sp1 &
                        - (rat1*(gyf1 + 0.5*dx1*gxyf1 + 0.5*dy1*gyyf1))&
                        /sp1)*((gxf1 + 0.5*dx1*gyyf1 - 0.5*dy1*gxyf1)/sp1 &
                        + (rat1*(0.5*dx1*gxyf1 - gyf1 + 0.5*dy1*gyyf1))/sp1))&
                        /(rat1**2 + 1)**2  !y1y2
                    cc = cc+1 
                    
                    row(cc) = vpairs(i,2); col(cc) = vpairs(i,1) 
                    valxx(cc) = valxx(cc-1)  ! x2x1
                    valxy(cc) = valyx(cc-1)  !x2y1
                    valyx(cc) = valxy(cc-1)  !y2x1
                    valyy(cc) = valyy(cc-1)  !y2y1
                    cc = cc+1 
                    
                    row(cc) = vpairs(i,1); col(cc) = vpairs(i,3) 
                    valxx(cc) = -(wt*((gyf1 - 0.5*dx1*gyxf1 + &
                        0.5*dy1*gxxf1)/sp1 - (rat1*(0.5*dx1*gxxf1 - &
                        gxf1 + 0.5*dy1*gyxf1))/sp1)*((gyf2 - 0.5*dx2*gyxf2 &
                        + 0.5*dy2*gxxf2)/sp2 - (rat2*(0.5*dx2*gxxf2 - &
                        gxf2 + 0.5*dy2*gyxf2))/sp2))/((rat1**2 + 1)*&
                        (rat2**2 + 1))  ! x1x3
                    valxy(cc) = (wt*((gyf1 - 0.5*dx1*gyxf1 + &
                        0.5*dy1*gxxf1)/sp1 - (rat1*(0.5*dx1*gxxf1 - &
                        gxf1 + 0.5*dy1*gyxf1))/sp1)*((gxf2 + 0.5*dx2*gyyf2 &
                        - 0.5*dy2*gxyf2)/sp2 + (rat2*(0.5*dx2*gxyf2 - &
                        gyf2 + 0.5*dy2*gyyf2))/sp2))/((rat1**2 + 1)*&
                        (rat2**2 + 1))  !x1y3
                    valyx(cc) = (wt*((gyf2 - 0.5*dx2*gyxf2 + &
                        0.5*dy2*gxxf2)/sp2 - (rat2*(0.5*dx2*gxxf2 - &
                        gxf2 + 0.5*dy2*gyxf2))/sp2)*((gxf1 + &
                        0.5*dx1*gyyf1 - 0.5*dy1*gxyf1)/sp1 + &
                        (rat1*(0.5*dx1*gxyf1 - gyf1 + 0.5*dy1*gyyf1))&
                        /sp1))/((rat1**2 + 1)*(rat2**2 + 1))  !y1x3
                    valyy(cc) = -(wt*((gxf1 + 0.5*dx1*gyyf1 - &
                        0.5*dy1*gxyf1)/sp1 + (rat1*(0.5*dx1*gxyf1 - &
                        gyf1 + 0.5*dy1*gyyf1))/sp1)*((gxf2 + 0.5*dx2*gyyf2 &
                        - 0.5*dy2*gxyf2)/sp2 + (rat2*(0.5*dx2*gxyf2 - &
                        gyf2 + 0.5*dy2*gyyf2))/sp2))/((rat1**2 + 1)*(rat2**2 + 1))  !y1y3
                    cc = cc+1 
                    
                    row(cc) = vpairs(i,3); col(cc) = vpairs(i,1) 
                    valxx(cc) = valxx(cc-1)  ! x3x1
                    valxy(cc) = valyx(cc-1)  !x3y1
                    valyx(cc) = valxy(cc-1)  !y3x1
                    valyy(cc) = valyy(cc-1)  !y3y1
                    cc = cc+1 
                    
                    row(cc) = vpairs(i,1); col(cc) = vpairs(i,4) 
                    valxx(cc) = (wt*((gyf2 + 0.5*dx2*gyxf2 - &
                        0.5*dy2*gxxf2)/sp2 + (rat2*(gxf2 + 0.5*dx2*gxxf2 + &
                        0.5*dy2*gyxf2))/sp2)*((gyf1 - 0.5*dx1*gyxf1 + &
                        0.5*dy1*gxxf1)/sp1 - (rat1*(0.5*dx1*gxxf1 - &
                        gxf1 + 0.5*dy1*gyxf1))/sp1))/((rat1**2 + 1)*(rat2**2 + 1))  ! x1x4
                    valxy(cc) = -(wt*((gxf2 - 0.5*dx2*gyyf2 + &
                        0.5*dy2*gxyf2)/sp2 - (rat2*(gyf2 + 0.5*dx2*gxyf2 + &
                        0.5*dy2*gyyf2))/sp2)*((gyf1 - 0.5*dx1*gyxf1 +&
                        0.5*dy1*gxxf1)/sp1 - (rat1*(0.5*dx1*gxxf1 - &
                        gxf1 + 0.5*dy1*gyxf1))/sp1))/((rat1**2 + 1)&
                        *(rat2**2 + 1))  !x1y4
                    valyx(cc) = -(wt*((gyf2 + 0.5*dx2*gyxf2 - &
                        0.5*dy2*gxxf2)/sp2 + (rat2*(gxf2 + 0.5*dx2*gxxf2 &
                        + 0.5*dy2*gyxf2))/sp2)*((gxf1 + 0.5*dx1*gyyf1 &
                        - 0.5*dy1*gxyf1)/sp1 + (rat1*(0.5*dx1*gxyf1 - &
                        gyf1 + 0.5*dy1*gyyf1))/sp1))/((rat1**2 + 1)&
                        *(rat2**2 + 1))  !y1x4
                    valyy(cc) = (wt*((gxf2 - 0.5*dx2*gyyf2 + &
                        0.5*dy2*gxyf2)/sp2 - (rat2*(gyf2 + 0.5*dx2*gxyf2 + &
                        0.5*dy2*gyyf2))/sp2)*((gxf1 + 0.5*dx1*gyyf1 - &
                        0.5*dy1*gxyf1)/sp1 + (rat1*(0.5*dx1*gxyf1 - &
                        gyf1 + 0.5*dy1*gyyf1))/sp1))/((rat1**2 + 1)&
                        *(rat2**2 + 1))  !y1y4
                    cc = cc+1 
                    
                    row(cc) = vpairs(i,4); col(cc) = vpairs(i,1) 
                    valxx(cc) = valxx(cc-1)  ! x4x1
                    valxy(cc) = valyx(cc-1)  !x4y1
                    valyx(cc) = valxy(cc-1)  !y4x1
                    valyy(cc) = valyy(cc-1)  !y4y1
                    cc = cc+1 
                    
                    row(cc) = vpairs(i,2); col(cc) = vpairs(i,2) 
                    valxx(cc) = (wt*((gyf1 + 0.5*dx1*gyxf1 - &
                        0.5*dy1*gxxf1)/sp1 + (rat1*(gxf1 + 0.5*dx1*gxxf1 &
                        + 0.5*dy1*gyxf1))/sp1)**2)/(rat1**2 + 1)**2 - &
                        wt*(atan(rat1) - atan(rat2))*(((1.0*gyxf1 + &
                        0.25*dx1*gyxxf1 - 0.25*dy1*gxxxf1)/sp1 - &
                        (2*(gxf1 + 0.5*dx1*gxxf1 + 0.5*dy1*gyxf1)*&
                        (gyf1 + 0.5*dx1*gyxf1 - 0.5*dy1*gxxf1))/sp1**2 - &
                        (2*rat1*(gxf1 + 0.5*dx1*gxxf1 + 0.5*dy1*gyxf1)**2)&
                        /sp1**2 + (rat1*(1.0*gxxf1 + 0.25*dx1*gxxxf1 +&
                        0.25*dy1*gyxxf1))/sp1)/(rat1**2 + 1) + &
                        (((2*rat1*(gyf1 + 0.5*dx1*gyxf1 - &
                        0.5*dy1*gxxf1))/sp1 + (2*rat1**2*(gxf1 + &
                        0.5*dx1*gxxf1 + 0.5*dy1*gyxf1))/sp1)*&
                        ((gyf1 + 0.5*dx1*gyxf1 - 0.5*dy1*gxxf1)&
                        /sp1 + (rat1*(gxf1 + 0.5*dx1*gxxf1 + &
                        0.5*dy1*gyxf1))/sp1))/(rat1**2 + 1)**2)  ! x2x2
                    valxy(cc) = wt*(atan(rat1) - atan(rat2))*&
                        (((0.5*gxxf1 - 0.5*gyyf1 - 0.25*dx1*gyyxf1 + &
                        0.25*dy1*gxyxf1)/sp1 - ((gxf1 + 0.5*dx1*gxxf1 + &
                        0.5*dy1*gyxf1)*(gxf1 - 0.5*dx1*gyyf1 + &
                        0.5*dy1*gxyf1))/sp1**2 + ((gyf1 + 0.5*dx1*gyxf1 - &
                        0.5*dy1*gxxf1)*(gyf1 + 0.5*dx1*gxyf1 + &
                        0.5*dy1*gyyf1))/sp1**2 - (rat1*(0.5*gxyf1 + &
                        0.5*gyxf1 + 0.25*dx1*gxyxf1 + 0.25*dy1*gyyxf1))&
                        /sp1 + (2*rat1*(gxf1 + 0.5*dx1*gxxf1 + &
                        0.5*dy1*gyxf1)*(gyf1 + 0.5*dx1*gxyf1 + &
                        0.5*dy1*gyyf1))/sp1**2)/(rat1**2 + 1) + &
                        (((2*rat1*(gxf1 - 0.5*dx1*gyyf1 + &
                        0.5*dy1*gxyf1))/sp1 - (2*rat1**2*(gyf1 + &
                        0.5*dx1*gxyf1 + 0.5*dy1*gyyf1))/sp1)*&
                        ((gyf1 + 0.5*dx1*gyxf1 - 0.5*dy1*gxxf1)/sp1 + &
                        (rat1*(gxf1 + 0.5*dx1*gxxf1 + 0.5*dy1*gyxf1))/sp1))&
                        /(rat1**2 + 1)**2) - (wt*((gyf1 + &
                        0.5*dx1*gyxf1 - 0.5*dy1*gxxf1)/sp1 + &
                        (rat1*(gxf1 + 0.5*dx1*gxxf1 + &
                        0.5*dy1*gyxf1))/sp1)*((gxf1 - 0.5*dx1*gyyf1 + &
                        0.5*dy1*gxyf1)/sp1 - (rat1*(gyf1 + 0.5*dx1*gxyf1 + &
                        0.5*dy1*gyyf1))/sp1))/(rat1**2 + 1)**2  !x2y2
                    valyx(cc) = valxy(cc)  !y2x2
                    valyy(cc) = (wt*((gxf1 - 0.5*dx1*gyyf1 + &
                        0.5*dy1*gxyf1)/sp1 - (rat1*(gyf1 + 0.5*dx1*gxyf1 + &
                        0.5*dy1*gyyf1))/sp1)**2)/(rat1**2 + 1)**2 + &
                        wt*(atan(rat1) - atan(rat2))*(((1.0*gxyf1 - &
                        0.25*dx1*gyyyf1 + 0.25*dy1*gxyyf1)/sp1 - &
                        (2*(gxf1 - 0.5*dx1*gyyf1 + 0.5*dy1*gxyf1)*&
                        (gyf1 + 0.5*dx1*gxyf1 + 0.5*dy1*gyyf1))/sp1**2 + &
                        (2*rat1*(gyf1 + 0.5*dx1*gxyf1 + &
                        0.5*dy1*gyyf1)**2)/sp1**2 - &
                        (rat1*(1.0*gyyf1 + 0.25*dx1*gxyyf1 + &
                        0.25*dy1*gyyyf1))/sp1)/(rat1**2 + 1) - &
                        (((2*rat1*(gxf1 - 0.5*dx1*gyyf1 + &
                        0.5*dy1*gxyf1))/sp1 - (2*rat1**2*(gyf1 + &
                        0.5*dx1*gxyf1 + 0.5*dy1*gyyf1))/sp1)*((gxf1 - &
                        0.5*dx1*gyyf1 + 0.5*dy1*gxyf1)/sp1 - &
                        (rat1*(gyf1 + 0.5*dx1*gxyf1 + 0.5*dy1*gyyf1))&
                        /sp1))/(rat1**2 + 1)**2)  !y2y2
                    cc = cc+1 
                    
                    row(cc) = vpairs(i,2); col(cc) = vpairs(i,3) 
                    valxx(cc) = (wt*((gyf1 + 0.5*dx1*gyxf1 - &
                        0.5*dy1*gxxf1)/sp1 + (rat1*(gxf1 + &
                        0.5*dx1*gxxf1 + 0.5*dy1*gyxf1))/sp1)*&
                        ((gyf2 - 0.5*dx2*gyxf2 + 0.5*dy2*gxxf2)/sp2 - &
                        (rat2*(0.5*dx2*gxxf2 - gxf2 + 0.5*dy2*gyxf2))&
                        /sp2))/((rat1**2 + 1)*(rat2**2 + 1))  ! x2x3
                    valxy(cc) = -(wt*((gyf1 + 0.5*dx1*gyxf1 - &
                        0.5*dy1*gxxf1)/sp1 + (rat1*(gxf1 + &
                        0.5*dx1*gxxf1 + 0.5*dy1*gyxf1))/sp1)*&
                        ((gxf2 + 0.5*dx2*gyyf2 - 0.5*dy2*gxyf2)/sp2 + &
                        (rat2*(0.5*dx2*gxyf2 - gyf2 + 0.5*dy2*gyyf2))&
                        /sp2))/((rat1**2 + 1)*(rat2**2 + 1))  !x2y3
                    valyx(cc) = -(wt*((gxf1 - 0.5*dx1*gyyf1 + &
                        0.5*dy1*gxyf1)/sp1 - (rat1*(gyf1 + &
                        0.5*dx1*gxyf1 + 0.5*dy1*gyyf1))/sp1)*&
                        ((gyf2 - 0.5*dx2*gyxf2 + 0.5*dy2*gxxf2)/sp2 - &
                        (rat2*(0.5*dx2*gxxf2 - gxf2 + 0.5*dy2*gyxf2))&
                        /sp2))/((rat1**2 + 1)*(rat2**2 + 1))  !y2x3
                    valyy(cc) = (wt*((gxf1 - 0.5*dx1*gyyf1 + &
                        0.5*dy1*gxyf1)/sp1 - (rat1*(gyf1 + &
                        0.5*dx1*gxyf1 + 0.5*dy1*gyyf1))/sp1)*&
                        ((gxf2 + 0.5*dx2*gyyf2 - 0.5*dy2*gxyf2)/sp2 + &
                        (rat2*(0.5*dx2*gxyf2 - gyf2 + &
                        0.5*dy2*gyyf2))/sp2))/((rat1**2 + 1)*&
                        (rat2**2 + 1))  !y2y3
                    cc = cc+1 
                    
                    row(cc) = vpairs(i,3); col(cc) = vpairs(i,2) 
                    valxx(cc) = valxx(cc-1)  ! x3x3
                    valxy(cc) = valyx(cc-1)  !x3y2
                    valyx(cc) = valxy(cc-1)  !y3x2
                    valyy(cc) = valyy(cc-1)  !y3y2
                    cc = cc+1 
                    
                    row(cc) = vpairs(i,2); col(cc) = vpairs(i,4) 
                    valxx(cc) = -(wt*((gyf1 + 0.5*dx1*gyxf1 - &
                        0.5*dy1*gxxf1)/sp1 + (rat1*(gxf1 + &
                        0.5*dx1*gxxf1 + 0.5*dy1*gyxf1))/sp1)*&
                        ((gyf2 + 0.5*dx2*gyxf2 - 0.5*dy2*gxxf2)/sp2 + &
                        (rat2*(gxf2 + 0.5*dx2*gxxf2 + &
                        0.5*dy2*gyxf2))/sp2))/((rat1**2 + 1)&
                        *(rat2**2 + 1))  ! x2x4
                    valxy(cc) = (wt*((gyf1 + 0.5*dx1*gyxf1 - &
                        0.5*dy1*gxxf1)/sp1 + (rat1*(gxf1 + &
                        0.5*dx1*gxxf1 + 0.5*dy1*gyxf1))/sp1)*&
                        ((gxf2 - 0.5*dx2*gyyf2 + 0.5*dy2*gxyf2)/sp2 - &
                        (rat2*(gyf2 + 0.5*dx2*gxyf2 + 0.5*dy2*gyyf2))&
                        /sp2))/((rat1**2 + 1)*(rat2**2 + 1))  !x2y4
                    valyx(cc) = (wt*((gyf2 + 0.5*dx2*gyxf2 - &
                        0.5*dy2*gxxf2)/sp2 + (rat2*(gxf2 + &
                        0.5*dx2*gxxf2 + 0.5*dy2*gyxf2))/sp2)*&
                        ((gxf1 - 0.5*dx1*gyyf1 + 0.5*dy1*gxyf1)/sp1 - &
                        (rat1*(gyf1 + 0.5*dx1*gxyf1 + &
                        0.5*dy1*gyyf1))/sp1))/((rat1**2 + 1)&
                        *(rat2**2 + 1))  !y2x4
                    valyy(cc) = -(wt*((gxf1 - 0.5*dx1*gyyf1 + &
                        0.5*dy1*gxyf1)/sp1 - (rat1*(gyf1 + &
                        0.5*dx1*gxyf1 + 0.5*dy1*gyyf1))/sp1)*&
                        ((gxf2 - 0.5*dx2*gyyf2 + 0.5*dy2*gxyf2)/sp2 - &
                        (rat2*(gyf2 + 0.5*dx2*gxyf2 + 0.5*dy2*gyyf2))&
                        /sp2))/((rat1**2 + 1)*(rat2**2 + 1))  !y2y4
                    cc = cc+1 
                    
                    row(cc) = vpairs(i,4); col(cc) = vpairs(i,2) 
                    valxx(cc) = valxx(cc-1)  ! x4x2
                    valxy(cc) = valyx(cc-1)  !x4y2
                    valyx(cc) = valxy(cc-1)  !y4x2
                    valyy(cc) = valyy(cc-1)  !y4y2
                    cc = cc+1 
                    
                    row(cc) = vpairs(i,3); col(cc) = vpairs(i,3) 
                    valxx(cc) = (wt*((gyf2 - 0.5*dx2*gyxf2 + &
                        0.5*dy2*gxxf2)/sp2 - (rat2*(0.5*dx2*gxxf2 -  &
                        gxf2 + 0.5*dy2*gyxf2))/sp2)**2)/(rat2**2 &
                        + 1)**2 - wt*(((1.0*gyxf2 - 0.25*dx2*gyxxf2 &
                        + 0.25*dy2*gxxxf2)/sp2 + &
                        (2*rat2*(0.5*dx2*gxxf2 - gxf2 + &
                        0.5*dy2*gyxf2)**2)/sp2**2 - (2*(gyf2 - &
                        0.5*dx2*gyxf2 + 0.5*dy2*gxxf2)*(0.5*dx2*gxxf2 - &
                        gxf2 + 0.5*dy2*gyxf2))/sp2**2 - &
                        (rat2*(0.25*dx2*gxxxf2 - 1.0*gxxf2 + &
                        0.25*dy2*gyxxf2))/sp2)/(rat2**2 + 1) + &
                        (((2*rat2**2*(0.5*dx2*gxxf2 - gxf2 + &
                        0.5*dy2*gyxf2))/sp2 - (2*rat2*(gyf2 - &
                        0.5*dx2*gyxf2 + 0.5*dy2*gxxf2))/sp2)*&
                        ((gyf2 - 0.5*dx2*gyxf2 + 0.5*dy2*gxxf2)/sp2 - &
                        (rat2*(0.5*dx2*gxxf2 - gxf2 + 0.5*dy2*gyxf2))&
                        /sp2))/(rat2**2 + 1)**2)*(atan(rat1) - atan(rat2))  ! x3x3
                    valxy(cc) = - wt*(atan(rat1) - atan(rat2))*&
                        ((((gxf2 + 0.5*dx2*gyyf2 - 0.5*dy2*gxyf2)*&
                        (0.5*dx2*gxxf2 - gxf2 + 0.5*dy2*gyxf2))/sp2**2&
                         - (0.5*gxxf2 - 0.5*gyyf2 + 0.25*dx2*gyyxf2 - &
                        0.25*dy2*gxyxf2)/sp2 - ((gyf2 - 0.5*dx2*gyxf2 + &
                        0.5*dy2*gxxf2)*(0.5*dx2*gxyf2 - gyf2 + &
                        0.5*dy2*gyyf2))/sp2**2 + (rat2*(0.5*gxyf2 + &
                        0.5*gyxf2 - 0.25*dx2*gxyxf2 - &
                        0.25*dy2*gyyxf2))/sp2 + (2*rat2*(0.5*dx2*gxxf2 &
                        - gxf2 + 0.5*dy2*gyxf2)*(0.5*dx2*gxyf2 - gyf2 + &
                        0.5*dy2*gyyf2))/sp2**2)/(rat2**2 + 1) + &
                        (((2*rat2**2*(0.5*dx2*gxyf2 - gyf2 + &
                        0.5*dy2*gyyf2))/sp2 + (2*rat2*(gxf2 + &
                        0.5*dx2*gyyf2 - 0.5*dy2*gxyf2))/sp2)*&
                        ((gyf2 - 0.5*dx2*gyxf2 + 0.5*dy2*gxxf2)/sp2 - &
                        (rat2*(0.5*dx2*gxxf2 - gxf2 + &
                        0.5*dy2*gyxf2))/sp2))/(rat2**2 + 1)**2) - &
                        (wt*((gyf2 - 0.5*dx2*gyxf2 + &
                        0.5*dy2*gxxf2)/sp2 - (rat2*(0.5*dx2*gxxf2 - &
                        gxf2 + 0.5*dy2*gyxf2))/sp2)*((gxf2 + &
                        0.5*dx2*gyyf2 - 0.5*dy2*gxyf2)/sp2 + &
                        (rat2*(0.5*dx2*gxyf2 - gyf2 + &
                        0.5*dy2*gyyf2))/sp2))/(rat2**2 + 1)**2  !x3y3
                    valyx(cc) = valxy(cc)  !y3x3
                    valyy(cc) = wt*(((1.0*gxyf2 + 0.25*dx2*gyyyf2 - &
                        0.25*dy2*gxyyf2)/sp2 - (2*rat2*(0.5*dx2*gxyf2 - &
                        gyf2 + 0.5*dy2*gyyf2)**2)/sp2**2 - &
                        (2*(gxf2 + 0.5*dx2*gyyf2 - 0.5*dy2*gxyf2)&
                        *(0.5*dx2*gxyf2 - gyf2 + 0.5*dy2*gyyf2))&
                        /sp2**2 + (rat2*(0.25*dx2*gxyyf2 - 1.0*gyyf2 + &
                        0.25*dy2*gyyyf2))/sp2)/(rat2**2 + 1) + &
                        (((2*rat2**2*(0.5*dx2*gxyf2 - gyf2 + &
                        0.5*dy2*gyyf2))/sp2 + (2*rat2*(gxf2 + &
                        0.5*dx2*gyyf2 - 0.5*dy2*gxyf2))/sp2)*&
                        ((gxf2 + 0.5*dx2*gyyf2 - 0.5*dy2*gxyf2)/sp2 + &
                        (rat2*(0.5*dx2*gxyf2 - gyf2 + &
                        0.5*dy2*gyyf2))/sp2))/(rat2**2 + 1)**2)*&
                        (atan(rat1) - atan(rat2)) + (wt*((gxf2 + &
                        0.5*dx2*gyyf2 - 0.5*dy2*gxyf2)/sp2 + &
                        (rat2*(0.5*dx2*gxyf2 - gyf2 + &
                        0.5*dy2*gyyf2))/sp2)**2)/(rat2**2 + 1)**2  !y3y3
                    cc = cc+1 
                    
                    row(cc) = vpairs(i,3); col(cc) = vpairs(i,4) 
                    valxx(cc) = wt*(atan(rat1) - atan(rat2))*&
                        (((0.25*dx2*gyxxf2 - 0.25*dy2*gxxxf2)/sp2 + &
                        ((gxf2 + 0.5*dx2*gxxf2 + 0.5*dy2*gyxf2)*&
                        (gyf2 - 0.5*dx2*gyxf2 + 0.5*dy2*gxxf2))/sp2**2&
                         - ((gyf2 + 0.5*dx2*gyxf2 - 0.5*dy2*gxxf2)*&
                        (0.5*dx2*gxxf2 - gxf2 + 0.5*dy2*gyxf2))/sp2**2 &
                        + (rat2*(0.25*dx2*gxxxf2 + 0.25*dy2*gyxxf2))&
                        /sp2 - (2*rat2*(gxf2 + 0.5*dx2*gxxf2 + &
                        0.5*dy2*gyxf2)*(0.5*dx2*gxxf2 - gxf2 + &
                        0.5*dy2*gyxf2))/sp2**2)/(rat2**2 + 1) - &
                        (((2*rat2*(gyf2 + 0.5*dx2*gyxf2 - &
                        0.5*dy2*gxxf2))/sp2 + (2*rat2**2*(gxf2 + &
                        0.5*dx2*gxxf2 + 0.5*dy2*gyxf2))/sp2)*((gyf2 - &
                        0.5*dx2*gyxf2 + 0.5*dy2*gxxf2)/sp2 - &
                        (rat2*(0.5*dx2*gxxf2 - gxf2 + 0.5*dy2*gyxf2))&
                        /sp2))/(rat2**2 + 1)**2) - (wt*((gyf2 + &
                        0.5*dx2*gyxf2 - 0.5*dy2*gxxf2)/sp2 + &
                        (rat2*(gxf2 + 0.5*dx2*gxxf2 + 0.5*dy2*gyxf2))&
                        /sp2)*((gyf2 - 0.5*dx2*gyxf2 + 0.5*dy2*gxxf2)&
                        /sp2 - (rat2*(0.5*dx2*gxxf2 - gxf2 + &
                        0.5*dy2*gyxf2))/sp2))/(rat2**2 + 1)**2  ! x3x4
                    valxy(cc) = wt*(atan(rat1) - atan(rat2))*&
                        ((((gyf2 - 0.5*dx2*gyxf2 + 0.5*dy2*gxxf2)*&
                        (gyf2 + 0.5*dx2*gxyf2 + 0.5*dy2*gyyf2))/sp2**2 &
                        - (0.5*gxxf2 + 0.5*gyyf2 - 0.25*dx2*gyyxf2 + &
                        0.25*dy2*gxyxf2)/sp2 + ((gxf2 - 0.5*dx2*gyyf2 &
                        + 0.5*dy2*gxyf2)*(0.5*dx2*gxxf2 - gxf2 + &
                        0.5*dy2*gyxf2))/sp2**2 + (rat2*(0.5*gyxf2 - &
                        0.5*gxyf2 + 0.25*dx2*gxyxf2 + 0.25*dy2*gyyxf2))&
                        /sp2 - (2*rat2*(gyf2 + 0.5*dx2*gxyf2 + &
                        0.5*dy2*gyyf2)*(0.5*dx2*gxxf2 - gxf2 + &
                        0.5*dy2*gyxf2))/sp2**2)/(rat2**2 + 1) + &
                        (((2*rat2*(gxf2 - 0.5*dx2*gyyf2 + &
                        0.5*dy2*gxyf2))/sp2 - (2*rat2**2*(gyf2 + &
                        0.5*dx2*gxyf2 + 0.5*dy2*gyyf2))/sp2)*&
                        ((gyf2 - 0.5*dx2*gyxf2 + 0.5*dy2*gxxf2)/sp2 - &
                        (rat2*(0.5*dx2*gxxf2 - gxf2 + 0.5*dy2*gyxf2))&
                        /sp2))/(rat2**2 + 1)**2) + (wt*((gxf2 - &
                        0.5*dx2*gyyf2 + 0.5*dy2*gxyf2)/sp2 - &
                        (rat2*(gyf2 + 0.5*dx2*gxyf2 + 0.5*dy2*gyyf2))&
                        /sp2)*((gyf2 - 0.5*dx2*gyxf2 + &
                        0.5*dy2*gxxf2)/sp2 - (rat2*(0.5*dx2*gxxf2 - &
                        gxf2 + 0.5*dy2*gyxf2))/sp2))/(rat2**2 + 1)**2  !x3y4
                    valyx(cc) = (wt*((gyf2 + 0.5*dx2*gyxf2 - &
                        0.5*dy2*gxxf2)/sp2 + (rat2*(gxf2 + &
                        0.5*dx2*gxxf2 + 0.5*dy2*gyxf2))/sp2)*&
                        ((gxf2 + 0.5*dx2*gyyf2 - 0.5*dy2*gxyf2)/sp2 + &
                        (rat2*(0.5*dx2*gxyf2 - gyf2 + 0.5*dy2*gyyf2))&
                        /sp2))/(rat2**2 + 1)**2 - wt*(atan(rat1) - &
                        atan(rat2))*((((gxf2 + 0.5*dx2*gxxf2 + &
                        0.5*dy2*gyxf2)*(gxf2 + 0.5*dx2*gyyf2 - &
                        0.5*dy2*gxyf2))/sp2**2 - (0.5*gxxf2 + &
                        0.5*gyyf2 + 0.25*dx2*gyyxf2 - 0.25*dy2*gxyxf2)&
                        /sp2 + ((gyf2 + 0.5*dx2*gyxf2 - 0.5*dy2*gxxf2)&
                        *(0.5*dx2*gxyf2 - gyf2 + 0.5*dy2*gyyf2))&
                        /sp2**2 - (rat2*(0.5*gxyf2 - 0.5*gyxf2 + &
                        0.25*dx2*gxyxf2 + 0.25*dy2*gyyxf2))/sp2 + &
                        (2*rat2*(gxf2 + 0.5*dx2*gxxf2 + 0.5*dy2*gyxf2)&
                        *(0.5*dx2*gxyf2 - gyf2 + 0.5*dy2*gyyf2))&
                        /sp2**2)/(rat2**2 + 1) - (((2*rat2*(gyf2 + &
                        0.5*dx2*gyxf2 - 0.5*dy2*gxxf2))/sp2 + &
                        (2*rat2**2*(gxf2 + 0.5*dx2*gxxf2 + &
                        0.5*dy2*gyxf2))/sp2)*((gxf2 + 0.5*dx2*gyyf2 - &
                        0.5*dy2*gxyf2)/sp2 + (rat2*(0.5*dx2*gxyf2 - &
                        gyf2 + 0.5*dy2*gyyf2))/sp2))/(rat2**2 + 1)**2)  !y3x4
                    valyy(cc) = wt*(atan(rat1) - atan(rat2))*&
                        (((0.25*dx2*gyyyf2 - 0.25*dy2*gxyyf2)/sp2 - &
                        ((gxf2 + 0.5*dx2*gyyf2 - 0.5*dy2*gxyf2)*&
                        (gyf2 + 0.5*dx2*gxyf2 + 0.5*dy2*gyyf2))&
                        /sp2**2 + ((gxf2 - 0.5*dx2*gyyf2 + 0.5*dy2*&
                        gxyf2)*(0.5*dx2*gxyf2 - gyf2 + 0.5*dy2*gyyf2))&
                        /sp2**2 + (rat2*(0.25*dx2*gxyyf2 + &
                        0.25*dy2*gyyyf2))/sp2 - (2*rat2*(gyf2 + &
                        0.5*dx2*gxyf2 + 0.5*dy2*gyyf2)*(0.5*dx2*gxyf2 &
                        - gyf2 + 0.5*dy2*gyyf2))/sp2**2)/(rat2**2 + 1) &
                        - (((2*rat2*(gxf2 - 0.5*dx2*gyyf2 + &
                        0.5*dy2*gxyf2))/sp2 - (2*rat2**2*(gyf2 + &
                        0.5*dx2*gxyf2 + 0.5*dy2*gyyf2))/sp2)*((gxf2 + &
                        0.5*dx2*gyyf2 - 0.5*dy2*gxyf2)/sp2 + &
                        (rat2*(0.5*dx2*gxyf2 - gyf2 + 0.5*dy2*gyyf2))&
                        /sp2))/(rat2**2 + 1)**2) - (wt*((gxf2 - &
                        0.5*dx2*gyyf2 + 0.5*dy2*gxyf2)/sp2 - &
                        (rat2*(gyf2 + 0.5*dx2*gxyf2 + 0.5*dy2*gyyf2))&
                        /sp2)*((gxf2 + 0.5*dx2*gyyf2 - 0.5*dy2*gxyf2)&
                        /sp2 + (rat2*(0.5*dx2*gxyf2 - gyf2 + 0.5*dy2&
                        *gyyf2))/sp2))/(rat2**2 + 1)**2  !y3y4
                    cc = cc+1 
                    
                    row(cc) = vpairs(i,4); col(cc) = vpairs(i,3) 
                    valxx(cc) = valxx(cc-1)  ! x4x3
                    valxy(cc) = valyx(cc-1)  !x4y3
                    valyx(cc) = valxy(cc-1)  !y4x3
                    valyy(cc) = valyy(cc-1)  !y4y3
                    cc = cc+1 
                    
                    row(cc) = vpairs(i,4); col(cc) = vpairs(i,4) 
                    valxx(cc) = (wt*((gyf2 + 0.5*dx2*gyxf2 - &
                        0.5*dy2*gxxf2)/sp2 + (rat2*(gxf2 + &
                        0.5*dx2*gxxf2 + 0.5*dy2*gyxf2))/sp2)**2)&
                        /(rat2**2 + 1)**2 + wt*(atan(rat1) - &
                        atan(rat2))*(((1.0*gyxf2 + 0.25*dx2*gyxxf2 - &
                        0.25*dy2*gxxxf2)/sp2 - (2*(gxf2 + 0.5*dx2*gxxf2 &
                        + 0.5*dy2*gyxf2)*(gyf2 + 0.5*dx2*gyxf2 - &
                        0.5*dy2*gxxf2))/sp2**2 - (2*rat2*(gxf2 + &
                        0.5*dx2*gxxf2 + 0.5*dy2*gyxf2)**2)/sp2**2 + &
                        (rat2*(1.0*gxxf2 + 0.25*dx2*gxxxf2 + &
                        0.25*dy2*gyxxf2))/sp2)/(rat2**2 + 1) + &
                        (((2*rat2*(gyf2 + 0.5*dx2*gyxf2 - &
                        0.5*dy2*gxxf2))/sp2 + (2*rat2**2*(gxf2 + &
                        0.5*dx2*gxxf2 + 0.5*dy2*gyxf2))/sp2)*((gyf2 + &
                        0.5*dx2*gyxf2 - 0.5*dy2*gxxf2)/sp2 + &
                        (rat2*(gxf2 + 0.5*dx2*gxxf2 + 0.5*dy2*gyxf2))&
                        /sp2))/(rat2**2 + 1)**2)  ! x4x4
                    valxy(cc) = - wt*(atan(rat1) - atan(rat2))*&
                        (((0.5*gxxf2 - 0.5*gyyf2 - 0.25*dx2*gyyxf2 + &
                        0.25*dy2*gxyxf2)/sp2 - ((gxf2 + 0.5*dx2*gxxf2 &
                        + 0.5*dy2*gyxf2)*(gxf2 - 0.5*dx2*gyyf2 + &
                        0.5*dy2*gxyf2))/sp2**2 + ((gyf2 + &
                        0.5*dx2*gyxf2 - 0.5*dy2*gxxf2)*(gyf2 + &
                        0.5*dx2*gxyf2 + 0.5*dy2*gyyf2))/sp2**2 - &
                        (rat2*(0.5*gxyf2 + 0.5*gyxf2 + &
                        0.25*dx2*gxyxf2 + 0.25*dy2*gyyxf2))/sp2 + &
                        (2*rat2*(gxf2 + 0.5*dx2*gxxf2 + &
                        0.5*dy2*gyxf2)*(gyf2 + 0.5*dx2*gxyf2 + &
                        0.5*dy2*gyyf2))/sp2**2)/(rat2**2 + 1) + &
                        (((2*rat2*(gxf2 - 0.5*dx2*gyyf2 + &
                        0.5*dy2*gxyf2))/sp2 - (2*rat2**2*(gyf2 + &
                        0.5*dx2*gxyf2 + 0.5*dy2*gyyf2))/sp2)*((gyf2 + &
                        0.5*dx2*gyxf2 - 0.5*dy2*gxxf2)/sp2 + &
                        (rat2*(gxf2 + 0.5*dx2*gxxf2 + 0.5*dy2*gyxf2))&
                        /sp2))/(rat2**2 + 1)**2) - (wt*((gyf2 + &
                        0.5*dx2*gyxf2 - 0.5*dy2*gxxf2)/sp2 + &
                        (rat2*(gxf2 + 0.5*dx2*gxxf2 + &
                        0.5*dy2*gyxf2))/sp2)*((gxf2 - 0.5*dx2*gyyf2 + &
                        0.5*dy2*gxyf2)/sp2 - (rat2*(gyf2 + &
                        0.5*dx2*gxyf2 + 0.5*dy2*gyyf2))/sp2))&
                        /(rat2**2 + 1)**2  !x4y4
                    valyx(cc) = valxy(cc)  !y4x4
                    valyy(cc) = (wt*((gxf2 - 0.5*dx2*gyyf2 + &
                        0.5*dy2*gxyf2)/sp2 - (rat2*(gyf2 + &
                        0.5*dx2*gxyf2 + 0.5*dy2*gyyf2))/sp2)**2)&
                        /(rat2**2 + 1)**2 - wt*(atan(rat1) - &
                        atan(rat2))*(((1.0*gxyf2 - &
                        0.25*dx2*gyyyf2 + 0.25*dy2*gxyyf2)/sp2 - &
                        (2*(gxf2 - 0.5*dx2*gyyf2 + 0.5*dy2*gxyf2)*&
                        (gyf2 + 0.5*dx2*gxyf2 + 0.5*dy2*gyyf2))/sp2**2 &
                        + (2*rat2*(gyf2 + 0.5*dx2*gxyf2 + 0.5*dy2*&
                        gyyf2)**2)/sp2**2 - (rat2*(1.0*gyyf2 + 0.25*&
                        dx2*gxyyf2 + 0.25*dy2*gyyyf2))/sp2)/(rat2**2 + &
                        1) - (((2*rat2*(gxf2 - 0.5*dx2*gyyf2 + &
                        0.5*dy2*gxyf2))/sp2 - (2*rat2**2*(gyf2 + &
                        0.5*dx2*gxyf2 + 0.5*dy2*gyyf2))/sp2)*&
                        ((gxf2 - 0.5*dx2*gyyf2 + 0.5*dy2*gxyf2)/sp2 - &
                        (rat2*(gyf2 + 0.5*dx2*gxyf2 + 0.5*dy2*gyyf2))&
                        /sp2))/(rat2**2 + 1)**2)  !y4y4
                    cc = cc+1 

                end do

                ! Build full hessian
                hessJ%row = [row, row, row+vert%ntot, row+vert%ntot]
                hessJ%col = [col, col+vert%ntot, col, col+vert%ntot]
                hessJ%val = [valxx, valxy, valyx, valyy]

                ! Scale
                hessJ%val = lambda*hessJ%val

                ! Housekeeping
                deallocate(row)
                deallocate(col)
                deallocate(valxx)
                deallocate(valxy)
                deallocate(valyx)
                deallocate(valyy)

            case default

                ! Not implemented, throw error
                call gdErrorHandler('EvaluateCostFunctionLR: hessian' &
                    // ' not yet implemented for this design variable' &
                    // ' type')

            end select

        end if

        ! Housekeeping
        !=============
        ! Deallocate auxiliary arrays
        deallocate(xv, yv, xfv1 , xfv2 , yfv1 , &
            yfv2 , gxfv1 , gxfv2 , gyfv1 , gyfv2 , &
            dxv1 , dxv2 , dyv1 , dyv2 , &
            spv1 , spv2 , cpv1 , cpv2 , ratv1 , &
            ratv2)

        if (dohessian) then

            select case (trim(designvariables%type))

            case ('coordinates', 'coordinates_desiredflux')

                ! Deallocate
                deallocate(gxxfv1 , gxxfv2 , gxyfv1 , gxyfv2 , &
                    gyxfv1 , gyxfv2 , gyyfv1 , gyyfv2 )
                deallocate(gxxxfv1 , gxxxfv2 , gxxyfv1 , &
                    gxxyfv2 , gxyxfv1 , gxyxfv2 , &
                    gxyyfv1 , gxyyfv2 , gyxxfv1 , gyxxfv2 , &
                    gyxyfv1 , gyxyfv2 , gyyyfv1 , gyyyfv2 )

            case default

                ! Not implemented, throw error
                call gdErrorHandler('EvaluateCostFunctionFAD: gradient' &
                    // ' not yet implemented for this design variable' &
                    // ' type')

            end select 

        elseif (dogradient) then

            select case (trim(designvariables%type))

            case ('coordinates', 'coordinates_desiredflux')

                ! Deallocate
                deallocate(gxxfv1 , gxxfv2 , gxyfv1 , gxyfv2 , &
                    gyxfv1 , gyxfv2 , gyyfv1 , gyyfv2 )

            end select

        end if

        ! Deassociate
        end associate
        

    end subroutine

    ! Cost function data writing 
    subroutine WriteCostFunctionDataFAD(costfunction, grid)

        ! Description
        !============
        ! Write out the cost function data for the FAD cost function.
        ! Here, this consists of the vertex pair data in IDn, xn, yn 
        ! format

        ! Declare variables
        !==================
        ! Arguments
        class(CostFunctionFADUDT)       :: costfunction 
        type(GridUDT)                   :: grid

        ! Auxiliary
        integer(I8)                     :: ncol 

        integer(I8), allocatable        :: IDn(:, :) 
        real(R8), allocatable           :: xn(:, :), yn(:, :)
        character(:), allocatable       :: filename 

        ! Loop
        integer(I8)                     :: j

        ! Initialize
        !===========
        ! Set filename
        allocate(character(len('costfunction_vertexpairs_FAD')) :: filename)
        filename = 'costfunction_vertexpairs_FAD'

        ! Unpack
        associate(&
            vpairs      => costfunction%vpairs,         &
            nvpairs     => costfunction%nvpairs,        &
            x           => grid%vert%x,                 &
            y           => grid%vert%y)
        
        ! Allocate
        ncol = size(vpairs, 2)
        allocate(IDn(nvpairs, ncol), xn(nvpairs, ncol), &
            yn(nvpairs, ncol))

        ! Loop
        do j = 1, ncol 
            IDn(:, j) = vpairs(:, j) 
            xn(:, j) = x(vpairs(:, j)) 
            yn(:, j) = y(vpairs(:, j)) 
        end do

        ! Call writer
        !============
        call WriteVertexPairData(IDn, xn, yn, filename)

        ! Housekeeping
        !=============
        end associate
        deallocate(IDn, xn, yn)
        


    end subroutine

    ! Housekeeping
    subroutine AllocateCostFunctionFAD(costfunction, nv, nvn)

        ! Description
        !============
        ! Allocate, assume number of vpairs given

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionFADUDT)       :: costfunction
        integer(I8)                     :: nv, nvn

        ! Allocate
        !=========
        allocate(costfunction%vpairs(nv, nvn))
        allocate(costfunction%wt(nv))

    end subroutine

    subroutine DeallocateCostFunctionFAD(costfunction)

        ! Description
        !============
        ! Deallocate

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionFADUDT)       :: costfunction

        ! Deallocate
        !===========
        if (allocated(costfunction%vpairs)) then 
            deallocate(costfunction%vpairs)
            deallocate(costfunction%wt)
        end if

    end subroutine

    subroutine DestroyCostFunctionFAD(costfunction)

        ! Description
        !============
        ! Deallocate

        ! Declare variables
        !==================
        ! Arguments
        type(CostfunctionFADUDT)       :: costfunction

        ! Destroy
        !========
        call costfunction%Deallocate()

    end subroutine

    !------------------------------------------------------------------!
    !                             FACE ANGLE                           !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCostFunctionFA(costfunction, grid, &
        magneticField, environment, options)

        ! Description
        !============
        ! Initialize the cost function and its parameters based on the 
        ! grid, magnetic field, and environment structures. Here, the 
        ! length ratio cost function is initialized, which requires

        ! Modules

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionFAUDT)            :: costfunction
        type(GridUDT)                       :: grid
        type(MagneticFieldUDT)              :: magneticField 
        type(EnvironmentUDT)                :: environment 
        type(CostFunctionOptionsUDT)        :: options

        ! Loop variables
        integer(I8)                         :: i

        ! Auxiliary variables
        type(StructuredPLF2DDistanceDFUDT)  :: dfwt
        integer(I8), allocatable            :: vpairs(:, :), tv(:), &
            tID(:) 

        real(R8), allocatable               :: wt(:), xv(:, :), &
            yv(:, :), dx(:), dy(:), gxf(:), gyf(:), dp(:), xf(:), yf(:)

        logical, allocatable                :: isvesselvertex(:), &
            isvesselface(:)

        ! Data
        
        ! Initialize
        !===========
        ! Set the scaling constant
        costfunction%lambda = options%FA%lambda

        ! Initialize temporary arrays (too big for now, trim later)
        allocate(vpairs(grid%face%ntot, 2))

        ! Associate
        associate(&
            face        => grid%face,   &
            x           => grid%vert%x, &
            y           => grid%vert%y, &
            fieldlineID => grid%vert%fieldlineID,   &
            nvpairs     => costfunction%nvpairs)

        ! Initialize counter
        nvpairs = 0

        ! Determine vertex pairs
        !=======================
        ! Check which faces are on vessel
        call DetermineVesselVertices(isvesselvertex, isvesselface, &
            grid)

        ! Loop
        do i = 1, face%ntot
            ! Skip if vessel face
            if (isvesselface(i)) then 
                cycle 
            end if

            ! Get vertices
            tv = face%vert(i, :)

            ! Get fieldline ID
            tID = fieldlineID(tv)

            ! Check if the face can be added
            if (tID(1) /= tID(2) ) then 
                ! Update counter
                nvpairs = nvpairs + 1

                ! Add pair
                vpairs(nvpairs, :) = tv 
            end if
        end do

        ! Allocate and add
        call costfunction%Allocate(nvpairs)
        costfunction%vpairs = vpairs(1:nvpairs, :)

        ! Check orientation
        !==================
        ! Flip vertex pairs if dotproduct is smaller than zero
        allocate(xv(nvpairs, 2), yv(nvpairs, 2), gxf(nvpairs), gyf(nvpairs))
        allocate(dy(nvpairs), dx(nvpairs))
        do i = 1, 2
            xv(:, i) = x(costfunction%vpairs(:, i))
            yv(:, i) = y(costfunction%vpairs(:, i))
        end do
        dx = xv(:, 2) - xv(:, 1)
        dy = yv(:, 2) - yv(:, 1)
        xf = 0.5*(xv(:, 1) + xv(:, 2))
        yf = 0.5*(yv(:, 1) + yv(:, 2))
        call magneticField%interp%Evaluate(xf, yf, 1, 0, gxf)
        call magneticField%interp%Evaluate(xf, yf, 0, 1, gyf)
        dp = dx*gxf + dy*gyf 
        do i = 1, nvpairs
            if (dp(i) < 0) then 
                ! flip
                costfunction%vpairs(i, :) = costfunction%vpairs(i, 2:1:-1)
            end if
        end do 

        ! Determine weigths
        !==================
        ! Initialize (magnetic field as dummy since unsigned anyway)
        call dfwt%Initialize(magneticField%interp, &
            environment%vessel%plfvessel, environment%vessel%plfvessel, &
            options%FA%weightatvessel, options%FA%weightatinf, &
            options%FA%decaylength, 'unsigned')

        ! Evaluate
        allocate(wt(nvpairs))
        call dfwt%Evaluate(xf, yf, wt)

        ! Visualize
        call dfwt%Visualize([minval(x), maxval(x)], &
            [minval(y), maxval(y)], 100, 100, 'costfunctionFA_weights')

        ! Add
        costfunction%wt = wt

        ! Housekeeping
        !=============
        ! End associate
        end associate

        ! Write data
        !===========
        if (options%writedata == 1) then 
            call costfunction%WriteData(grid)
        end if 

    end subroutine

    ! Cost function evaluation
    subroutine EvaluateCostFunctionFA(costfunction, J, gradJ, hessJ, &
        grid, magneticField, environment, dogradient, dohessian, &
        designvariables)

        ! Description
        !============
        ! Evaluate the cost function, the gradient and its hessian. The 
        ! cost function penalizes the angle between the magnetic field 
        ! and the face normal. It is assumed that the optimal angle is
        ! zero. 

        ! Notes:
        !=======
        ! Note 1: for ease, we compute the dot product between the 
        ! tangent and the magnetic field normal, which somewhat boils 
        ! down to the same thing. 

        ! Note 2: we could simplify the cost function expression by 
        ! penalizing the dot product instead of the actual angle. 
        ! However, in the future one may desire to have a non-zero 
        ! desired angle, which is now easily adjusted by adding a 
        ! desired theta value. Anyway, this implementation may be 
        ! improved. 

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionFAUDT)        :: costfunction 
        real(R8)                        :: J
        real(R8), allocatable           :: gradJ(:) 
        type(MySparseUDT)               :: hessJ 
        type(GridUDT)                   :: grid 
        type(MagneticFieldUDT)          :: magneticField 
        type(EnvironmentUDT)            :: environment
        logical                         :: dogradient, dohessian 
        class(DesignVariablesGDUDT)     :: designvariables

        ! Loop variables
        integer(I8)                     :: i, cc 

        ! Auxiliary
        integer(I8), allocatable        :: row(:), col(:) 

        real(R8)                        :: wti, gxf, gyf, dx, dy, dp, cp, &
            rat, theta, gxxf, gxyf, gyxf, gyyf, gxxxf, gxyyf, &
            gyxxf, gyyyf, gyyxf, gxyxf
        real(R8), allocatable           :: valxx(:),  valxy(:), &
            valyx(:), valyy(:), xv(:, :), yv(:, :), xfv(:), yfv(:), &
            gxfv(:), gyfv(:), dxv(:), dyv(:), dpv(:), cpv(:), ratv(:), &
            thetav(:), gxxfv(:), gxyfv(:), gyxfv(:), gyyfv(:), &
            gxxxfv(:), gxxyfv(:), gxyyfv(:), gyxxfv(:), gyxyfv(:), &
            gyyyfv(:)
                                        
        ! Associate
        !==========
        associate(&
            vert    => grid%vert, &
            vpairs  => costfunction%vpairs, &
            nvpairs => costfunction%nvpairs, &
            x       => grid%vert%x, &
            y       => grid%vert%y, &
            nv      => grid%vert%ntot, &
            lambda  => costfunction%lambda, &
            wt      => costfunction%wt)

        ! Initialize
        !===========
        ! Cost function
        J = 0

        ! Gradient
        gradJ(:) = 0

        ! Precompute
        !===========
        ! Coordinates
        allocate(xv(nvpairs, 2), yv(nvpairs, 2), gxfv(nvpairs), gyfv(nvpairs))
        allocate(dxv(nvpairs), dyv(nvpairs), dpv(nvpairs), cpv(nvpairs), &
            ratv(nvpairs), thetav(nvpairs))
        do i = 1, 2
            xv(:, i) = x(vpairs(:, i))
            yv(:, i) = y(vpairs(:, i))
        end do
        xfv = 0.5*(xv(:, 1) + xv(:, 2))
        yfv = 0.5*(yv(:, 1) + yv(:, 2))

        ! Magnetic field
        call magneticField%interp%Evaluate(xfv, yfv, 1, 0, gxfv)
        call magneticField%interp%Evaluate(xfv, yfv, 0, 1, gyfv)

        ! Face vectors
        dxv = xv(:, 2) - xv(:, 1)
        dyv = yv(:, 2) - yv(:, 1)

        dpv = dxv*gxfv + dyv*gyfv
        cpv = dxv*gyfv - dyv*gxfv 

        ratv = cpv/dpv 
        thetav = atan(ratv)

        ! Compute cost function
        !======================
        ! Compute
        J = sum(0.5*wt*thetav**2)

        ! Scale
        J = lambda*J

        ! Compute derivatives
        !====================
        select case (designvariables%type) 

        case ('coordinates', 'coordinates_desiredflux') ! no flux contributions

            ! Initialize
            if (.not. allocated(hessJ%row)) then 
                ! Allocate
                hessJ%nval = 16*nvpairs 
                call hessJ%Allocate()

            end if 
            
            ! Precompute
            if (dogradient .or. dohessian) then 
                ! Magnetic field vector derivatives
                allocate(gxxfv(nvpairs), gxyfv(nvpairs), gyyfv(nvpairs))
                call magneticField%interp%Evaluate(xfv, yfv, 2, 0, gxxfv)
                call magneticField%interp%Evaluate(xfv, yfv, 1, 1, gxyfv)
                call magneticField%interp%Evaluate(xfv, yfv, 0, 2, gyyfv)
                gyxfv = gxyfv ! symmetric, for ease

            end if
            if (dohessian) then 
                ! Additional derivatives
                allocate(gxxxfv(nvpairs), gxxyfv(nvpairs), &
                    gxyyfv(nvpairs), gyyyfv(nvpairs))
                call magneticField%interp%Evaluate(xfv, yfv, 3, 0, gxxxfv)
                call magneticField%interp%Evaluate(xfv, yfv, 2, 1, gxxyfv)
                call magneticField%interp%Evaluate(xfv, yfv, 1, 2, gxyyfv)
                call magneticField%interp%Evaluate(xfv, yfv, 0, 3, gyyyfv)
                gyxxfv = gxxyfv 
                gyxyfv = gxyyfv 

                ! Allocate local variables
                allocate(valxx(4*nvpairs), valxy(4*nvpairs), &
                    valyx(4*nvpairs), valyy(4*nvpairs), row(4*nvpairs), &
                    col(4*nvpairs))
                
            end if 

            ! Gradient
            if (dogradient) then 
                do i = 1, nvpairs
                    ! Unpack
                    wti     = wt(i)
                    gxf     = gxfv(i)
                    gyf     = gyfv(i)
                    dx      = dxv(i)
                    dy      = dyv(i)
                    dp      = dpv(i)
                    cp      = cpv(i)
                    rat     = ratv(i)
                    theta   = thetav(i)
                    gxxf    = gxxfv(i)
                    gxyf    = gxyfv(i)
                    gyxf    = gyxfv(i)
                    gyyf    = gyyfv(i)

                    ! Evaluate gradient
                    gradJ(vpairs(i, 1)) = gradJ(vpairs(i, 1)) + &
                        -(theta*wti*((gyf - 0.5*dx*gyxf + 0.5*dy*gxxf)/dp &
                        + (rat*(0.5*dx*gxxf - gxf + 0.5*dy*gyxf))/dp))/(rat**2 + 1) !x1
                    gradJ(vpairs(i, 2)) = gradJ(vpairs(i, 2)) + &
                        (theta*wti*((gyf + 0.5*dx*gyxf - 0.5*dy*gxxf)/dp &
                        - (rat*(gxf + 0.5*dx*gxxf + 0.5*dy*gyxf))/dp))/(rat**2 + 1) ! x2

                    gradJ(vpairs(i,1)+nv) = gradJ(vpairs(i,1)+nv) + &
                        (theta*wti*((gxf + 0.5*dx*gyyf - 0.5*dy*gxyf)/dp &
                        - (rat*(0.5*dx*gxyf - gyf + 0.5*dy*gyyf))/dp))/(rat**2 + 1) !y1
                    gradJ(vpairs(i,2)+nv) = gradJ(vpairs(i,2)+nv) + &
                        -(theta*wti*((gxf - 0.5*dx*gyyf + 0.5*dy*gxyf)/dp &
                        + (rat*(gyf + 0.5*dx*gxyf + 0.5*dy*gyyf))/dp))/(rat**2 + 1) !y2

                end do
            end if 

            ! Scale
            gradJ = lambda*gradJ

            ! Hessian
            if (dohessian) then 
                ! Initialize
                cc = 1
                do i = 1, nvpairs
                    ! Unpack
                    wti     = wt(i)
                    gxf     = gxfv(i)
                    gyf     = gyfv(i)
                    dx      = dxv(i)
                    dy      = dyv(i)
                    dp      = dpv(i)
                    cp      = cpv(i)
                    rat     = ratv(i)
                    theta   = thetav(i)
                    gxxf    = gxxfv(i)
                    gxyf    = gxyfv(i)
                    gyxf    = gyxfv(i)
                    gyyf    = gyyfv(i)

                    gxxxf = gxxxfv(i)
                    gxyxf = gxxyfv(i)
                    gxyyf = gxyyfv(i)
                    gyxxf = gyxxfv(i)
                    gyyxf = gyxyfv(i)
                    gyyyf = gyyyfv(i)

                    ! v1 v1
                    row(cc) = vpairs(i, 1)  
                    col(cc) = vpairs(i, 1) 
                    valxx(cc) = (wti*((gyf - 0.5*dx*gyxf + 0.5*dy*gxxf)/dp &
                        + (rat*(0.5*dx*gxxf - gxf + 0.5*dy*gyxf))/dp)**2) &
                        /(rat **2 + 1) **2 - (theta*wti*((1.0*gyxf - &
                        0.25*dx*gyxxf + 0.25*dy*gxxxf)/dp - (2*(gyf &
                        - 0.5*dx*gyxf + 0.5*dy*gxxf)*(0.5*dx*gxxf &
                        - gxf + 0.5*dy*gyxf))/dp **2 + (rat*(0.25*dx*gxxxf &
                        - 1.0*gxxf + 0.25*dy*gyxxf))/dp - (2*rat*(0.5*dx*gxxf &
                        - gxf + 0.5*dy*gyxf) **2)/dp **2))/(rat **2 + 1) &
                        - (theta*wti*((2*rat*(gyf - 0.5*dx*gyxf + 0.5*dy*gxxf))/dp &
                        + (2*rat **2*(0.5*dx*gxxf - gxf + 0.5*dy*gyxf))/dp)*((gyf &
                        - 0.5*dx*gyxf + 0.5*dy*gxxf)/dp + (rat*(0.5*dx*gxxf &
                        - gxf + 0.5*dy*gyxf))/dp))/(rat **2 + 1) **2  !x1x1
                    valxy(cc) = (theta*wti*((0.5*gxxf - 0.5*gyyf &
                        + 0.25*dx*gyyxf - 0.25*dy*gxyxf)/dp - ((gxf &
                        + 0.5*dx*gyyf - 0.5*dy*gxyf)*(0.5*dx*gxxf &
                        - gxf + 0.5*dy*gyxf))/dp **2 + ((gyf - 0.5*dx*gyxf &
                        + 0.5*dy*gxxf)*(0.5*dx*gxyf - gyf + 0.5*dy*gyyf))/dp**2 &
                        + (rat*(0.5*gxyf + 0.5*gyxf - 0.25*dx*gxyxf &
                        - 0.25*dy*gyyxf))/dp + (2*rat*(0.5*dx*gxxf &
                        - gxf + 0.5*dy*gyxf)*(0.5*dx*gxyf - gyf &
                        + 0.5*dy*gyyf))/dp **2))/(rat **2 + 1) &
                        - (wti*((gyf - 0.5*dx*gyxf + 0.5*dy*gxxf)/dp &
                        + (rat*(0.5*dx*gxxf - gxf + 0.5*dy*gyxf))/dp)&
                        *((gxf + 0.5*dx*gyyf - 0.5*dy*gxyf)/dp &
                        - (rat*(0.5*dx*gxyf - gyf + 0.5*dy*gyyf))/dp))&
                        /(rat **2 + 1) **2 + (theta*wti*((2*rat*(gxf &
                        + 0.5*dx*gyyf - 0.5*dy*gxyf))/dp - (2*rat**2*&
                        (0.5*dx*gxyf - gyf + 0.5*dy*gyyf))/dp)*((gyf &
                        - 0.5*dx*gyxf + 0.5*dy*gxxf)/dp + &
                        (rat*(0.5*dx*gxxf - gxf + 0.5*dy*gyxf))/dp))&
                        /(rat **2 + 1) **2  !x1y1
                    valyx(cc) = valxy(cc)  !y1x1
                    valyy(cc) = (wti*((gxf + 0.5*dx*gyyf - 0.5*dy*gxyf)/dp &
                        - (rat*(0.5*dx*gxyf - gyf + 0.5*dy*gyyf))/dp) **2)&
                        /(rat **2 + 1) **2 + (theta*wti*((1.0*gxyf &
                        + 0.25*dx*gyyyf - 0.25*dy*gxyyf)/dp &
                        - (2*(gxf + 0.5*dx*gyyf - 0.5*dy*gxyf)&
                        *(0.5*dx*gxyf - gyf + 0.5*dy*gyyf))/dp **2 &
                        - (rat*(0.25*dx*gxyyf - 1.0*gyyf + 0.25*dy*gyyyf))/dp &
                        + (2*rat*(0.5*dx*gxyf - gyf + 0.5*dy*gyyf)**2)/dp**2))&
                        /(rat **2 + 1) - (theta*wti*((2*rat*(gxf + 0.5*dx*gyyf &
                        - 0.5*dy*gxyf))/dp - (2*rat **2*(0.5*dx*gxyf &
                        - gyf + 0.5*dy*gyyf))/dp)*((gxf + 0.5*dx*gyyf &
                        - 0.5*dy*gxyf)/dp - (rat*(0.5*dx*gxyf - gyf &
                        + 0.5*dy*gyyf))/dp))/(rat **2 + 1) **2  !y1y1
                    cc = cc+1 
                    
                    ! v1 v2
                    row(cc) = vpairs(i, 1)  
                    col(cc) = vpairs(i, 2) 
                    valxx(cc) = (theta*wti*((0.25*dx*gyxxf - 0.25*dy*gxxxf)/dp &
                        + ((gxf + 0.5*dx*gxxf + 0.5*dy*gyxf)*(gyf - 0.5*dx*gyxf &
                        + 0.5*dy*gxxf))/dp **2 - ((gyf + 0.5*dx*gyxf &
                        - 0.5*dy*gxxf)*(0.5*dx*gxxf - gxf + 0.5*dy*gyxf))/dp **2 &
                        - (rat*(0.25*dx*gxxxf + 0.25*dy*gyxxf))/dp + (2*rat*(gxf &
                        + 0.5*dx*gxxf + 0.5*dy*gyxf)*(0.5*dx*gxxf - gxf &
                        + 0.5*dy*gyxf))/dp **2))/(rat **2 + 1) - (wti*((gyf &
                        + 0.5*dx*gyxf - 0.5*dy*gxxf)/dp - (rat*(gxf &
                        + 0.5*dx*gxxf + 0.5*dy*gyxf))/dp)*((gyf - 0.5*dx*gyxf &
                        + 0.5*dy*gxxf)/dp + (rat*(0.5*dx*gxxf - gxf &
                        + 0.5*dy*gyxf))/dp))/(rat **2 + 1) **2 &
                        + (theta*wti*((2*rat*(gyf + 0.5*dx*gyxf - 0.5*dy*gxxf))/dp &
                        - (2*rat **2*(gxf + 0.5*dx*gxxf + 0.5*dy*gyxf))/dp)*((gyf &
                        - 0.5*dx*gyxf + 0.5*dy*gxxf)/dp + (rat*(0.5*dx*gxxf - gxf &
                        + 0.5*dy*gyxf))/dp))/(rat **2 + 1) **2  !x1x2
                    valxy(cc) = (wti*((gxf - 0.5*dx*gyyf + 0.5*dy*gxyf)/dp &
                        + (rat*(gyf + 0.5*dx*gxyf + 0.5*dy*gyyf))/dp)*((gyf &
                        - 0.5*dx*gyxf + 0.5*dy*gxxf)/dp + (rat*(0.5*dx*gxxf &
                        - gxf + 0.5*dy*gyxf))/dp))/(rat **2 + 1) **2 &
                        + (theta*wti*(((gyf - 0.5*dx*gyxf + 0.5*dy*gxxf)*(gyf &
                        + 0.5*dx*gxyf + 0.5*dy*gyyf))/dp **2 - (0.5*gxxf &
                        + 0.5*gyyf - 0.25*dx*gyyxf + 0.25*dy*gxyxf)/dp &
                        + ((gxf - 0.5*dx*gyyf + 0.5*dy*gxyf)*(0.5*dx*gxxf &
                        - gxf + 0.5*dy*gyxf))/dp **2 - (rat*(0.5*gyxf - 0.5*gxyf &
                        + 0.25*dx*gxyxf + 0.25*dy*gyyxf))/dp + (2*rat*(gyf &
                        + 0.5*dx*gxyf + 0.5*dy*gyyf)*(0.5*dx*gxxf - gxf &
                        + 0.5*dy*gyxf))/dp **2))/(rat **2 + 1) &
                        - (theta*wti*((2*rat*(gxf - 0.5*dx*gyyf + 0.5*dy*gxyf))/dp &
                        + (2*rat **2*(gyf + 0.5*dx*gxyf + 0.5*dy*gyyf))/dp)&
                        *((gyf - 0.5*dx*gyxf + 0.5*dy*gxxf)/dp + (rat*(0.5*dx*gxxf &
                        - gxf + 0.5*dy*gyxf))/dp))/(rat **2 + 1) **2  !x1y2
                    valyx(cc) = (wti*((gyf + 0.5*dx*gyxf - 0.5*dy*gxxf)/dp &
                        - (rat*(gxf + 0.5*dx*gxxf + 0.5*dy*gyxf))/dp)*((gxf &
                        + 0.5*dx*gyyf - 0.5*dy*gxyf)/dp - (rat*(0.5*dx*gxyf &
                        - gyf + 0.5*dy*gyyf))/dp))/(rat **2 + 1) **2 &
                        - (theta*wti*(((gxf + 0.5*dx*gxxf + 0.5*dy*gyxf)*(gxf &
                        + 0.5*dx*gyyf - 0.5*dy*gxyf))/dp **2 - (0.5*gxxf &
                        + 0.5*gyyf + 0.25*dx*gyyxf - 0.25*dy*gxyxf)/dp &
                        + ((gyf + 0.5*dx*gyxf - 0.5*dy*gxxf)*(0.5*dx*gxyf &
                        - gyf + 0.5*dy*gyyf))/dp **2 + (rat*(0.5*gxyf - 0.5*gyxf &
                        + 0.25*dx*gxyxf + 0.25*dy*gyyxf))/dp - (2*rat*(gxf &
                        + 0.5*dx*gxxf + 0.5*dy*gyxf)*(0.5*dx*gxyf - gyf &
                        + 0.5*dy*gyyf))/dp **2))/(rat **2 + 1) &
                        - (theta*wti*((2*rat*(gyf + 0.5*dx*gyxf - 0.5*dy*gxxf))/dp &
                        - (2*rat **2*(gxf + 0.5*dx*gxxf + 0.5*dy*gyxf))/dp)*((gxf &
                        + 0.5*dx*gyyf - 0.5*dy*gxyf)/dp - (rat*(0.5*dx*gxyf - gyf &
                        + 0.5*dy*gyyf))/dp))/(rat **2 + 1) **2  !y1x2
                    valyy(cc) = (theta*wti*((0.25*dx*gyyyf - 0.25*dy*gxyyf)/dp &
                        - ((gxf + 0.5*dx*gyyf - 0.5*dy*gxyf)*(gyf + 0.5*dx*gxyf &
                        + 0.5*dy*gyyf))/dp **2 + ((gxf - 0.5*dx*gyyf &
                        + 0.5*dy*gxyf)*(0.5*dx*gxyf - gyf + 0.5*dy*gyyf))/dp **2 &
                        - (rat*(0.25*dx*gxyyf + 0.25*dy*gyyyf))/dp &
                        + (2*rat*(gyf + 0.5*dx*gxyf + 0.5*dy*gyyf)*(0.5*dx*gxyf &
                        - gyf + 0.5*dy*gyyf))/dp **2))/(rat **2 + 1) &
                        - (wti*((gxf - 0.5*dx*gyyf + 0.5*dy*gxyf)/dp + (rat*(gyf &
                        + 0.5*dx*gxyf + 0.5*dy*gyyf))/dp)*((gxf + 0.5*dx*gyyf &
                        - 0.5*dy*gxyf)/dp - (rat*(0.5*dx*gxyf - gyf &
                        + 0.5*dy*gyyf))/dp))/(rat **2 + 1) **2 &
                        + (theta*wti*((2*rat*(gxf - 0.5*dx*gyyf &
                        + 0.5*dy*gxyf))/dp + (2*rat **2*(gyf + 0.5*dx*gxyf &
                        + 0.5*dy*gyyf))/dp)*((gxf + 0.5*dx*gyyf - 0.5*dy*gxyf)/dp &
                        - (rat*(0.5*dx*gxyf - gyf + 0.5*dy*gyyf))/dp))&
                        /(rat **2 + 1) **2  !y1y2
                    cc = cc+1 
                    
                    ! v2 v1
                    row(cc) = vpairs(i,2)  
                    col(cc) = vpairs(i,1) 
                    valxx(cc) = valxx(cc-1)  ! x2x1
                    valxy(cc) = valyx(cc-1)  !x2y1
                    valyx(cc) = valxy(cc-1)  !y2x1
                    valyy(cc) = valyy(cc-1)  !y2y1
                    cc = cc+1 
                    
                    ! v2 v2
                    row(cc) = vpairs(i,2)  
                    col(cc) = vpairs(i,2) 
                    valxx(cc) = (wti*((gyf + 0.5*dx*gyxf - 0.5*dy*gxxf)/dp &
                        - (rat*(gxf + 0.5*dx*gxxf + 0.5*dy*gyxf))/dp) **2)/(rat **2 &
                        + 1) **2 + (theta*wti*((1.0*gyxf + 0.25*dx*gyxxf &
                        - 0.25*dy*gxxxf)/dp - (2*(gxf + 0.5*dx*gxxf &
                        + 0.5*dy*gyxf)*(gyf + 0.5*dx*gyxf - 0.5*dy*gxxf))/dp **2 &
                        + (2*rat*(gxf + 0.5*dx*gxxf + 0.5*dy*gyxf) **2)/dp **2 &
                        - (rat*(1.0*gxxf + 0.25*dx*gxxxf + 0.25*dy*gyxxf))/dp))&
                        /(rat **2 + 1) - (theta*wti*((2*rat*(gyf + 0.5*dx*gyxf &
                        - 0.5*dy*gxxf))/dp - (2*rat **2*(gxf + 0.5*dx*gxxf &
                        + 0.5*dy*gyxf))/dp)*((gyf + 0.5*dx*gyxf - 0.5*dy*gxxf)/dp &
                        - (rat*(gxf + 0.5*dx*gxxf + 0.5*dy*gyxf))/dp))/(rat **2 + 1) **2  !x2x2
                    valxy(cc) = (theta*wti*((2*rat*(gxf - 0.5*dx*gyyf + 0.5*dy*gxyf))/dp &
                        + (2*rat **2*(gyf + 0.5*dx*gxyf + 0.5*dy*gyyf))/dp)&
                        *((gyf + 0.5*dx*gyxf - 0.5*dy*gxxf)/dp - (rat*(gxf &
                        + 0.5*dx*gxxf + 0.5*dy*gyxf))/dp))/(rat **2 + 1) **2 &
                        - (wti*((gyf + 0.5*dx*gyxf - 0.5*dy*gxxf)/dp &
                        - (rat*(gxf + 0.5*dx*gxxf + 0.5*dy*gyxf))/dp)*((gxf &
                        - 0.5*dx*gyyf + 0.5*dy*gxyf)/dp + (rat*(gyf + 0.5*dx*gxyf &
                        + 0.5*dy*gyyf))/dp))/(rat **2 + 1) **2 &
                        - (theta*wti*((0.5*gxxf - 0.5*gyyf - 0.25*dx*gyyxf &
                        + 0.25*dy*gxyxf)/dp - ((gxf + 0.5*dx*gxxf + 0.5*dy*gyxf)&
                        *(gxf - 0.5*dx*gyyf + 0.5*dy*gxyf))/dp **2 &
                        + ((gyf + 0.5*dx*gyxf - 0.5*dy*gxxf)*(gyf + 0.5*dx*gxyf &
                        + 0.5*dy*gyyf))/dp **2 + (rat*(0.5*gxyf + 0.5*gyxf &
                        + 0.25*dx*gxyxf + 0.25*dy*gyyxf))/dp - (2*rat*(gxf &
                        + 0.5*dx*gxxf + 0.5*dy*gyxf)*(gyf + 0.5*dx*gxyf &
                        + 0.5*dy*gyyf))/dp **2))/(rat **2 + 1)  !x2y2
                    valyx(cc) = valxy(cc)  !y2x2
                    valyy(cc) = (wti*((gxf - 0.5*dx*gyyf + 0.5*dy*gxyf)/dp + &
                        (rat*(gyf + 0.5*dx*gxyf + 0.5*dy*gyyf))/dp) **2)/(rat **2 + 1) **2 &
                        - (theta*wti*((1.0*gxyf - 0.25*dx*gyyyf + 0.25*dy*gxyyf)/dp &
                        - (2*(gxf - 0.5*dx*gyyf + 0.5*dy*gxyf)*(gyf + 0.5*dx*gxyf &
                        + 0.5*dy*gyyf))/dp **2 - (2*rat*(gyf + 0.5*dx*gxyf &
                        + 0.5*dy*gyyf) **2)/dp **2 + (rat*(1.0*gyyf + 0.25*dx*gxyyf &
                        + 0.25*dy*gyyyf))/dp))/(rat **2 + 1) &
                        - (theta*wti*((2*rat*(gxf - 0.5*dx*gyyf + 0.5*dy*gxyf))/dp &
                        + (2*rat **2*(gyf + 0.5*dx*gxyf + 0.5*dy*gyyf))/dp)&
                        *((gxf - 0.5*dx*gyyf + 0.5*dy*gxyf)/dp + (rat*(gyf &
                        + 0.5*dx*gxyf + 0.5*dy*gyyf))/dp))/(rat **2 + 1) **2  !y2y2
                    cc = cc+1 

                end do

                ! Build full hessian
                hessJ%row = [row, row, row+vert%ntot, row+vert%ntot]
                hessJ%col = [col, col+vert%ntot, col, col+vert%ntot]
                hessJ%val = [valxx, valxy, valyx, valyy]

                ! Scale
                hessJ%val = lambda*hessJ%val

            end if 

        case default 

            ! Throw error
            call gdErrorHandler('design variable type "' // designvariables%type &
                // '" not yet implemented for face angle cost function')

        end select

       
        ! Deassociate
        !============
        end associate

    end subroutine

    ! Cost function data writing 
    subroutine WriteCostFunctionDataFA(costfunction, grid)

        ! Description
        !============
        ! Write out the cost function data for the LR cost function.
        ! Here, this consists of the vertex pair data in IDn, xn, yn 
        ! format

        ! Declare variables
        !==================
        ! Arguments
        class(CostFunctionFAUDT)        :: costfunction 
        type(GridUDT)                   :: grid

        ! Auxiliary
        integer(I8)                     :: ncol, nrow 

        integer(I8), allocatable        :: IDn(:, :) 
        real(R8), allocatable           :: xn(:, :), yn(:, :)
        character(:), allocatable       :: filename 

        ! Loop
        integer(I8)                     :: j 

        ! Initialize
        !===========
        ! Set filename
        allocate(character(len('costfunction_vertexpairs_FA')) :: filename)
        filename = 'costfunction_vertexpairs_FA'

        ! Allocate
        nrow = size(costfunction%vpairs, 1)
        ncol = size(costfunction%vpairs, 2)
        allocate(IDn(nrow, ncol), xn(nrow, ncol), yn(nrow, ncol))

        ! Unpack
        associate(&
            vpairs      => costfunction%vpairs,         &
            x           => grid%vert%x,                 &
            y           => grid%vert%y)

        ! Loop
        do j = 1, ncol 
            IDn(:, j) = vpairs(:, j) 
            xn(:, j) = x(vpairs(:, j)) 
            yn(:, j) = y(vpairs(:, j)) 
        end do

        ! Call writer
        !============
        call WriteVertexPairData(IDn, xn, yn, filename)

        ! Housekeeping
        !=============
        end associate
        deallocate(IDn, xn, yn)
        


    end subroutine

    ! Housekeeping
    subroutine AllocateCostFunctionFA(costfunction, nvp)

        ! Description
        !============
        ! Allocate, assumed that costfunction%nvpairs is given

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionFAUDT)        :: costfunction
        integer(I8)                     :: nvp

        ! Allocate
        !=========
        allocate(costfunction%vpairs(nvp, 2))
        allocate(costfunction%wt(nvp))

    end subroutine

    subroutine DeallocateCostFunctionFA(costfunction)

        ! Description
        !============
        ! Deallocate

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionFAUDT)        :: costfunction

        ! Deallocate
        !===========
        if (allocated(costfunction%vpairs)) then 
            deallocate(costfunction%vpairs)
            deallocate(costfunction%wt)
        end if

    end subroutine

    subroutine DestroyCostFunctionFA(costfunction)

        ! Description
        !============
        ! Deallocate

        ! Declare variables
        !==================
        ! Arguments
        type(CostfunctionFAUDT)        :: costfunction

        ! Destroy
        !========
        call costfunction%Deallocate()

    end subroutine

    !------------------------------------------------------------------!
    !                           LENGTH RATIO 2                         !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCostFunctionLR2(costfunction, grid, &
        magneticField, environment, options)

        ! Description
        !============
        ! Initialize the cost function and its parameters based on the 
        ! grid, magnetic field, and environment structures. 

        ! Simply call the initialization of the original lenght ratio
        ! cost function. 

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionLRUDT2)           :: costfunction
        type(GridUDT)                       :: grid
        type(MagneticFieldUDT)              :: magneticField 
        type(EnvironmentUDT)                :: environment 
        type(CostFunctionOptionsUDT)        :: options
        
        ! Initialize
        !===========
        call costfunction%cfv_lr%Initialize(grid, magneticField, &
            environment, options)

        ! (Re)set the scaling constant
        costfunction%cfv_lr%lambda = options%LR%lambda ! seems to agree well with most grids

    end subroutine

    ! Cost function evaluation
    subroutine EvaluateCostFunctionLR2(costfunction, J, gradJ, hessJ, &
        grid, magneticField, environment, dogradient, dohessian, &
        designvariables)

        ! Description
        !============
        ! Evaluate the cost function, the gradient and its hessian. 
        ! Here, we simply call the same cost function twice, but switch
        ! the order of the indices and recompute the bias. 

        ! Notes:
        !=======
        ! Possible future performance improvements:
        ! - Allocating hessian stuff only once and storing indices, 
        ! since they don't change
        ! - Instead of recomputing auxiliary variables, store them. May
        ! not actually be better in terms of computational time, but 
        ! may lead to shorter and hence better maintainable code. 

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionLRUDT2)       :: costfunction 
        real(R8)                        :: J, J1, J2
        real(R8), allocatable           :: gradJ(:), gradJ1(:), &
            gradJ2(:) 
        type(MySparseUDT)               :: hessJ, hessJ1, hessJ2 
        type(GridUDT)                   :: grid 
        type(MagneticFieldUDT)          :: magneticField 
        type(EnvironmentUDT)            :: environment
        logical                         :: dogradient, dohessian 
        class(DesignVariablesGDUDT)     :: designvariables

        ! Loop variables
        integer(I8)                     :: i

        ! Auxiliary
        integer(I8), allocatable        :: tempvpairs(:,:)
        real(R8), allocatable           :: tempb0(:) 
                                        
        ! Initialize
        !===========
        ! Store original vertex pairs and bias
        allocate(tempvpairs, source=costfunction%cfv_lr%vpairs)
        allocate(tempb0, source=costfunction%cfv_lr%b0)
        tempvpairs = costfunction%cfv_lr%vpairs 
        tempb0  = costfunction%cfv_lr%b0
         
        ! Cost function
        J = 0
        J1 = 0
        J2 = 0

        ! Gradient
        gradJ(:) = 0
        allocate(gradJ1(size(gradJ)), gradJ2(size(gradJ)))

        ! Hessian
        hessJ1%nrow = hessJ%nrow 
        hessJ2%nrow = hessJ%nrow 
        hessJ1%ncol = hessJ%ncol 
        hessJ2%ncol = hessJ%ncol 

        ! Compute cost function
        !======================
        ! First contribution
        call costfunction%cfv_lr%Evaluate(J1, gradJ1, &
            hessJ1, grid, magneticField, environment, dogradient, &
            dohessian, designvariables)
        
        ! Adjust vertex pairs and bias
        do i = 1, maxval(costfunction%cfv_lr%nvpairs)
            costfunction%cfv_lr%vpairs(:, 2*i-1) = tempvpairs(:, 2*i)
            costfunction%cfv_lr%vpairs(:, 2*i) = tempvpairs(:, 2*i-1)
        end do
        costfunction%cfv_lr%b0(:) = 1/tempb0

        ! Second contribution
        call costfunction%cfv_lr%Evaluate(J2, gradJ2, &
            hessJ2, grid, magneticField, environment, dogradient, &
            dohessian, designvariables)

        ! Reset vertex pairs and bias
        costfunction%cfv_lr%b0(:) = tempb0
        costfunction%cfv_lr%vpairs(:,:) = tempvpairs

        ! Add
        J = J1 + J2 
        gradJ = gradJ1 + gradJ2
        hessJ = hessJ1 + hessJ2

        ! Housekeeping
        !=============
        deallocate(gradJ1, gradJ2)

    end subroutine

    ! Housekeeping
    subroutine AllocateCostFunctionLR2(costfunction, nv, nvn)

        ! Description
        !============
        ! Allocate, assumed that costfunction%nvpairs is given

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionLRUDT2)       :: costfunction
        integer(I8)                     :: nv, nvn

        ! Allocate
        !=========
        call costfunction%cfv_lr%Allocate(nv, nvn)

    end subroutine

    subroutine DeallocateCostFunctionLR2(costfunction)

        ! Description
        !============
        ! Deallocate

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionLRUDT2)       :: costfunction

        ! Deallocate
        !===========
        call costfunction%cfv_lr%Deallocate()

    end subroutine

    subroutine DestroyCostFunctionLR2(costfunction)

        ! Description
        !============
        ! Deallocate

        ! Declare variables
        !==================
        ! Arguments
        type(CostfunctionLRUDT2)       :: costfunction

        ! Destroy
        !========
        call costfunction%cfv_lr%Deallocate()

    end subroutine

    !------------------------------------------------------------------!
    !                       LENGTH RATIO, RADIAL 2                     !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCostFunctionLRrad2(costfunction, grid, &
        magneticField, environment, options)

        ! Description
        !============
        ! Initialize the cost function and its parameters based on the 
        ! grid, magnetic field, and environment structures. 

        ! Simply call the initialization of the original lenght ratio
        ! cost function. 

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionLRrad2UDT)           :: costfunction
        type(GridUDT)                       :: grid
        type(MagneticFieldUDT)              :: magneticField 
        type(EnvironmentUDT)                :: environment 
        type(CostFunctionOptionsUDT)        :: options
        
        ! Initialize
        !===========
        call costfunction%cfv_lrrad%Initialize(grid, magneticField, &
            environment, options)

        ! (Re)set the scaling constant
        costfunction%cfv_lrrad%lambda = options%LRrad%lambda ! seems to agree well with most grids

    end subroutine

    ! Cost function evaluation
    subroutine EvaluateCostFunctionLRrad2(costfunction, J, gradJ, hessJ, &
        grid, magneticField, environment, dogradient, dohessian, &
        designvariables)

        ! Description
        !============
        ! Evaluate the cost function, the gradient and its hessian. 
        ! Here, we simply call the same cost function twice, but switch
        ! the order of the indices and recompute the bias. 

        ! Notes:
        !=======
        ! Possible future performance improvements:
        ! - Allocating hessian stuff only once and storing indices, 
        ! since they don't change
        ! - Instead of recomputing auxiliary variables, store them. May
        ! not actually be better in terms of computational time, but 
        ! may lead to shorter and hence better maintainable code. 

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionLRrad2UDT)    :: costfunction 
        real(R8)                        :: J, J1, J2
        real(R8), allocatable           :: gradJ(:), gradJ1(:), &
            gradJ2(:) 
        type(MySparseUDT)               :: hessJ, hessJ1, hessJ2 
        type(GridUDT)                   :: grid 
        type(MagneticFieldUDT)          :: magneticField 
        type(EnvironmentUDT)            :: environment
        logical                         :: dogradient, dohessian 
        class(DesignVariablesGDUDT)     :: designvariables

        ! Loop variables
        integer(I8)                     :: i

        ! Auxiliary
        integer(I8), allocatable        :: tempvpairs(:,:)
        real(R8), allocatable           :: tempb0(:) 
                                        
        ! Initialize
        !===========
        ! Store original vertex pairs and bias
        allocate(tempvpairs, source=costfunction%cfv_lrrad%vpairs)
        allocate(tempb0, source=costfunction%cfv_lrrad%b0)
        tempvpairs = costfunction%cfv_lrrad%vpairs 
        tempb0  = costfunction%cfv_lrrad%b0
         
        ! Cost function
        J = 0
        J1 = 0
        J2 = 0

        ! Gradient
        gradJ(:) = 0
        allocate(gradJ1(size(gradJ)), gradJ2(size(gradJ)))

        ! Hessian
        hessJ1%nrow = hessJ%nrow 
        hessJ2%nrow = hessJ%nrow 
        hessJ1%ncol = hessJ%ncol 
        hessJ2%ncol = hessJ%ncol 

        ! Compute cost function
        !======================
        ! First contribution
        call costfunction%cfv_lrrad%Evaluate(J1, gradJ1, &
            hessJ1, grid, magneticField, environment, dogradient, &
            dohessian, designvariables)
        
        ! Adjust vertex pairs and bias
        do i = 1, maxval(costfunction%cfv_lrrad%nvpairs)
            costfunction%cfv_lrrad%vpairs(:, 2*i-1) = tempvpairs(:, 2*i)
            costfunction%cfv_lrrad%vpairs(:, 2*i) = tempvpairs(:, 2*i-1)
        end do
        costfunction%cfv_lrrad%b0(:) = 1/tempb0

        ! Second contribution
        call costfunction%cfv_lrrad%Evaluate(J2, gradJ2, &
            hessJ2, grid, magneticField, environment, dogradient, &
            dohessian, designvariables)

        ! Reset vertex pairs and bias
        costfunction%cfv_lrrad%b0(:) = tempb0
        costfunction%cfv_lrrad%vpairs(:,:) = tempvpairs

        ! Add
        J = J1 + J2 
        gradJ = gradJ1 + gradJ2
        hessJ = hessJ1 + hessJ2

        ! Housekeeping
        !=============
        deallocate(gradJ1, gradJ2)

    end subroutine

    ! Housekeeping
    subroutine AllocateCostFunctionLRrad2(costfunction, nv, nvn)

        ! Description
        !============
        ! Allocate, assumed that costfunction%nvpairs is given

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionLRrad2UDT)    :: costfunction
        integer(I8)                     :: nv, nvn

        ! Allocate
        !=========
        call costfunction%cfv_lrrad%Allocate(nv, nvn)

    end subroutine

    subroutine DeallocateCostFunctionLRrad2(costfunction)

        ! Description
        !============
        ! Deallocate

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionLRrad2UDT)       :: costfunction

        ! Deallocate
        !===========
        call costfunction%cfv_lrrad%Deallocate()

    end subroutine

    subroutine DestroyCostFunctionLRrad2(costfunction)

        ! Description
        !============
        ! Deallocate

        ! Declare variables
        !==================
        ! Arguments
        type(CostfunctionLRrad2UDT)       :: costfunction

        ! Destroy
        !========
        call costfunction%cfv_lrrad%Deallocate()

    end subroutine


    !------------------------------------------------------------------!
    !                       PSI RATIO, PSI BASED                       !
    !------------------------------------------------------------------!

    ! Initialization (somewhat dummy)
    subroutine InitializeCostFunctionPRPB(costfunction, grid, &
        magneticField, environment, options)

        ! Description
        !============
        ! Initialize the cost function and its parameters based on the 
        ! grid, magnetic field, and environment structures. This 
        ! function is actually a dummy function, since the 
        ! initialization can only be done after the constraints have 
        ! been set up. Therefore, no real implementation is provided 
        ! here. The true initialization is done in 
        ! FinalizeInitializationCostFunctionPRPB, which is specific 
        ! for this cost function and should be called in the 
        ! GD optimization engine after initializing the constraints 
        ! etc. 

        ! Note: the only thing set here, is the scaling constant...

        ! Modules

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionPRPBUDT)          :: costfunction
        type(GridUDT)                       :: grid
        type(MagneticFieldUDT)              :: magneticField 
        type(EnvironmentUDT)                :: environment 
        type(CostFunctionOptionsUDT)        :: options

        ! Loop variables

        ! Auxiliary variables

        ! Data
        
        ! Initialize
        !===========
        ! Set the scaling constant
        costfunction%lambda = options%PRPB%lambda
        
    end subroutine

    ! True cost function initialization
    subroutine FinalizeInitializationCostFunctionPRPB(costfunction, &
        designvariables, grid, magneticField, environment, options)

        ! Description
        !============
        ! This is the actual initialization for the PRPB cost function.
        ! Here, we check whether the design variables are compatible 
        ! with the cost function. 

        ! Declare
        !========
        ! Arguments
        class(CostfunctionPRPBUDT)      :: costfunction
        class(DesignVariablesGDUDT)     :: designvariables
        type(GridUDT)                   :: grid 
        type(MagneticFieldUDT)          :: magneticField
        type(EnvironmentUDT)            :: environment
        type(CostFunctionOptionsUDT)    :: options

        ! Auxiliary
        type(Coordinates1DFieldDistanceDFUDT)    :: dfbias, dfwt
        integer(I8)                     :: nfsID, nxpind, tID
        integer(I8), allocatable        :: map2fsind(:), fsID(:), &
            xpind(:), order(:), psipairs(:, :), tvn(:), fsIDcounter(:), &
            tvnID(:), allfsIDs(:), tvID(:)

        real(R8), allocatable           :: fsPsi(:), b0v(:), wtv(:), &
            wtp(:), b0p(:), thispsival(:, :), sgn1(:), sgn2(:)

        logical, allocatable            :: mask(:), doflip(:), &
            issepvert(:)

        ! Loop
        integer(I8)                     :: npp, i, j, k

        ! Initialize
        !===========
        ! Unpack
        associate(&
            x               => grid%vert%x,             &
            y               => grid%vert%y,             &
            nv              => grid%vert%ntot,            &
            fieldlineID     => grid%vert%fieldlineID    &
            )

        ! Get psi values
        select type (designvariables)

        type is (DesignVariablesCoordinatesUDT) 

            ! Throw error
            call gdErrorHandler('PRPB cost function cannot be used for ' // &
                ' design variable type "coordinates"')

        type is (DesignVariablesCoordinatesFluxUDT) 

            ! Get psi values and indices of flux surfaces
            allocate(fsPsi(size(designvariables%psiind)))
            allocate(fsID(size(fsPsi)))
            fsPsi   = designvariables%phi(designvariables%psiind)
            fsID    = designvariables%fsID
            

        class default 

            ! Throw error
            call gdErrorHandler('PRPB cost function: unknown design variable type')

        end select

        ! Construct mapping
        nfsID = size(fsID, 1)
        allfsIDs = [(k, k = 1, maxval(fsID))]
        allocate(map2fsind(maxval(fsID)))
        do i = 1, nfsID
            map2fsind(fsID(i)) = i 
        end do 

        ! Initialize fsID counter
        allocate(fsIDcounter(maxval(fsID)))

        ! Compute desired distributions
        !==============================
        ! Get all separatrix vertices
        call DetermineXPoints(xpind, nxpind, order, grid)
        
        allocate(issepvert(grid%vert%ntot))
        issepvert = .false.
        do i = 1, nxpind 
            where (fieldlineID == fieldlineID(xpind(i))) issepvert = .true.
        end do 
        allocate(tvID(count(issepvert)))
        tvID = pack([(k, k = 1, grid%vert%ntot)], issepvert)

        ! Initialize distributions
        call dfbias%Initialize(magneticField%interp, x(tvID), y(tvID), &
            options%PRPB%biasatsep, options%PRPB%biasatinf, &
            options%PRPB%biasdecaylength, 'signed')
        call dfwt%Initialize(magneticField%interp, x(tvID), y(tvID), &
            options%PRPB%weightatsep, options%PRPB%weightatinf, &
            options%PRPB%weightdecaylength, 'unsigned')

        ! Evaluate
        allocate(b0v(size(x, 1)), wtv(size(x, 1)))
        call dfbias%Evaluate(x, y, b0v)
        call dfwt%Evaluate(x, y, wtv)

        ! Visualize
        call dfbias%Visualize([minval(x), maxval(x)], [minval(y), maxval(y)], &
            100, 100, 'costfunctionPRPB_desiredbias')
        call dfwt%Visualize([minval(x), maxval(x)], [minval(y), maxval(y)], &
            100, 100, 'costfunctionPRPB_weight')

        ! Set desired bias to one at separatrix nodes
        do i = 1, nxpind 
            where (fieldlineID == fieldlineID(xpind(i))) b0v = 1
        end do

        ! Initialize pairs too big
        allocate(psipairs(nv, 3), wtp(size(wtv)), b0p(size(b0v)))
        wtp = wtv 
        b0p = b0v 

        ! Set counter
        npp = 0

        ! Determine pairs
        !================
        ! Loop over all vertices
        do i = 1, nv
            ! Get vertex neighbours
            tvn = GetVertNeig(grid%vert, i)

            ! Get the ID of the coordinate line
            tID = fieldlineID(i)

            ! Check which vertices have a different ID (no zeros 
            ! allowed) 
            if ((tID /= 0) .and. (any(tID == fsID)) ) then 
                ! Get vertices with different ID
                mask = (fieldlineID(tvn) /= tID) .and. (fieldlineID(tvn) /= 0)
                allocate(tvnID(count(mask)))
                tvnID = pack(fieldlineID(tvn), mask)

                ! Check if only two unique flux surface IDs remain
                fsIDcounter = 0
                do j = 1, size(tvnID, 1)
                    fsIDcounter(tvnID(j)) = fsIDcounter(tvnID(j)) + 1
                end do
                if ( (maxval(fsIDcounter) < 2) .and. (count(fsIDcounter > 0) == 2)) then 
                    ! Update counter
                    npp = npp + 1

                    ! Add
                    psipairs(npp, :) = [tID, pack(allfsIDs, fsIDcounter > 0)]
                    wtp(npp) = wtv(i)
                    b0p(npp) = b0v(i)

                end if 

                ! Housekeeping
                deallocate(tvnID)
            end if 
        end do

        ! Map psi pairs to local indices and get current psi values
        allocate(thispsival(npp, 3))
        do i = 1, 3
            psipairs(1:npp, i)      = map2fsind(psipairs(1:npp, i))
            thispsival(1:npp, i)    = fsPsi(psipairs(1:npp, i))
        end do

        ! Check and switch
        sgn1 = thispsival(:, 1) - thispsival(:, 2)
        sgn2 = thispsival(:, 3) - thispsival(:, 1)
        if (any(sgn1*sgn2 < 0)) then 
            ! Shouldn't happen, throw error
            call gdErrorHandler('Something wrong with initial psi values')
        end if 
        doflip = (sgn1 < 0) .and. (sgn2 < 0)
        do i = 1, npp
            if (doflip(i)) then 
                psipairs(i, :) = psipairs(i, [1, 3, 2])
            end if 
        end do 

        ! Add to cost function
        costfunction%b0         = b0p(1:npp)
        costfunction%wt         = wtp(1:npp)
        costfunction%psipairs   = psipairs(1:npp, :)
        costfunction%npsipairs  = npp
        
        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Cost function evaluation
    subroutine EvaluateCostFunctionPRPB(costfunction, J, gradJ, hessJ, &
        grid, magneticField, environment, dogradient, dohessian, &
        designvariables)

        ! Description
        !============
        ! Evaluate the cost function, the gradient and its hessian. 

        ! Notes:
        !=======
        ! Possible future performance improvements:
        ! - Allocating hessian stuff only once and storing indices, 
        ! since they don't change
        ! - Instead of recomputing auxiliary variables, store them. May
        ! not actually be better in terms of computational time, but 
        ! may lead to shorter and hence better maintainable code. 

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionPRPBUDT)      :: costfunction
        real(R8)                        :: J
        real(R8), allocatable           :: gradJ(:) 
        type(MySparseUDT)               :: hessJ 
        type(GridUDT)                   :: grid 
        type(MagneticFieldUDT)          :: magneticField 
        type(EnvironmentUDT)            :: environment
        logical                         :: dogradient, dohessian 
        class(DesignVariablesGDUDT)     :: designvariables

        ! Loop variables
        integer(I8)                     :: i, cc 

        ! Auxiliary
        integer(I8)                     :: v1, v2, v3
        integer(I8), allocatable        :: row(:), col(:), psiind(:)

        real(R8)                        :: wti, &
            b0i, d1, d2, psi1, psi2, psi3
        real(R8), allocatable           ::  psi(:), psipairsval(:, :), &
            d1v(:), d2v(:), valpp(:) 
                            
        ! Associate
        !==========
        associate(&
            vert        => grid%vert, &
            psipairs    => costfunction%psipairs, &
            npsipairs   => costfunction%npsipairs, &
            b0          => costfunction%b0, &
            x           => grid%vert%x, &
            y           => grid%vert%y, &
            lambda      => costfunction%lambda, &
            wt          => costfunction%wt)

        ! Initialize
        !===========
        ! Cost function
        J = 0

        ! Gradient
        gradJ(:) = 0

        ! Compute cost function
        !======================
        ! Precompute
        select type (designvariables)

        type is (DesignVariablesCoordinatesFluxUDT)

            ! Get psi values
            allocate(psiind(size(designvariables%psiind)))
            allocate(psi(size(psiind)))
            psiind = designvariables%psiind
            Psi = designvariables%phi(psiind) 
            
        class default

            ! Throw error
            call gdErrorHandler('EvaluateCostFunctionPRPB: incompatible design variable type')

        end select

        ! Psi values
        allocate(psipairsval(npsipairs, 3), d1v(npsipairs), d2v(npsipairs))
        do i = 1, 3
            psipairsval(:, i) = psi(psipairs(:, i))
        end do

        ! Cost function contribution
        d1v = psipairsval(:, 1) - psipairsval(:, 2)
        d2v = psipairsval(:, 3) - psipairsval(:, 1)
        J = J + 0.5*sum(wt*(d1v/d2v - b0)**2)

        ! Scale
        J = lambda*J

        ! Compute gradient
        !=================
        if (dogradient) then 

            ! Check the design variables
            select case (trim(designvariables%type))

            case ('coordinates_desiredflux')

                ! Only flux value contributions
                do i = 1, npsipairs 
                    ! Unpack
                    psi1    = psipairsval(i, 1)
                    psi2    = psipairsval(i, 2)
                    psi3    = psipairsval(i, 3)
                    d1      = d1v(i)
                    d2      = d2v(i)
                    wti     = wt(i)
                    b0i     = b0(i)

                    ! Get indices
                    v1      = psiind(psipairs(i, 1))
                    v2      = psiind(psipairs(i, 2))
                    v3      = psiind(psipairs(i, 3))

                    ! Compute gradient contribution
                    gradJ(v1) = gradJ(v1) + &
                        -wti*(b0i - d1/d2)*(d1/d2**2 + 1/d2) !psi1
                    gradJ(v2) = gradJ(v2) + &
                        (wti*(b0i - d1/d2))/d2 !psi2
                    gradJ(v3) = gradJ(v3) + &
                        (d1*wti*(b0i - d1/d2))/d2**2 !psi3

                end do

            case default

                ! Not implemented, throw error
                call gdErrorHandler('EvaluateCostFunctionLR: gradient' &
                    // ' not yet implemented for this design variable' &
                    // ' type')

            end select

            ! Scale
            gradJ = lambda*gradJ

        end if

        ! Compute hessian
        !================
        ! Allocate the hessian (if not already done so)
        if (.not. allocated(hessJ%row)) then
            ! Allocate the sparse matrix
            select case (trim(designvariables%type))

            case ('coordinates_desiredflux')

                hessJ%nval = 9*npsipairs ! this should be exact and constant
                allocate(valpp(hessJ%nval))
                allocate(row(hessJ%nval), col(hessJ%nval))

            end select 
            call hessJ%Allocate()

        end if

        if (dohessian) then
            
            ! Initialize counter
            cc = 1

            ! Check the design variables
            select case (trim(designvariables%type))

            case ('coordinates_desiredflux')

                ! Only flux value contributions
                do i = 1, npsipairs 
                    ! Unpack
                    psi1    = psipairsval(i, 1)
                    psi2    = psipairsval(i, 2)
                    psi3    = psipairsval(i, 3)
                    d1      = d1v(i)
                    d2      = d2v(i)
                    wti     = wt(i)
                    b0i     = b0(i)

                    ! Get indices
                    v1      = psiind(psipairs(i, 1))
                    v2      = psiind(psipairs(i, 2))
                    v3      = psiind(psipairs(i, 3))

                    ! Compute hessian contributions
                    row(cc) = v1
                    col(cc) = v1
                    valpp(cc) = wti*(d1/d2**2 + 1/d2)**2 &
                        - wti*((2*d1)/d2**3 + 2/d2**2)*(b0i - d1/d2) !psi1psi1
                    cc = cc+1
                    
                    row(cc) = v1
                    col(cc) = v2
                    valpp(cc) = (wti*(b0i - d1/d2))/d2**2 &
                        - (wti*(d1/d2**2 + 1/d2))/d2 !psi1psi2
                    cc = cc+1
                    
                    row(cc) = v2
                    col(cc) = v1
                    valpp(cc) = valpp(cc-1) !psi2psi1
                    cc = cc+1
                    
                    row(cc) = v2
                    col(cc) = v2
                    valpp(cc) = wti/d2**2 !psi2psi2
                    cc = cc+1
                    
                    row(cc) = v1
                    col(cc) = v3
                    valpp(cc) = wti*(b0i - d1/d2)*((2*d1)/d2**3 + 1/d2**2) &
                        - (d1*wti*(d1/d2**2 + 1/d2))/d2**2 !psi1psi3
                    cc = cc+1
                    
                    row(cc) = v3
                    col(cc) = v1
                    valpp(cc) = valpp(cc-1)
                    cc = cc+1
                    
                    row(cc) = v2
                    col(cc) = v3
                    valpp(cc) = (d1*wti)/d2**3 - (wti*(b0i - d1/d2))/d2**2 !psi2psi3
                    cc = cc+1
                    
                    row(cc) = v3
                    col(cc) = v2
                    valpp(cc) = valpp(cc-1)
                    cc = cc+1
                    
                    row(cc) = v3
                    col(cc) = v3
                    valpp(cc) = (d1**2*wti)/d2**4 &
                        - (2*d1*wti*(b0i - d1/d2))/d2**3 !psi3psi3
                    cc = cc+1

                end do

                ! Scale
                valpp = lambda*valpp

                ! Construct hessian contributions
                hessJ%val = valpp 
                hessJ%col = col 
                hesSJ%row = row
                
            case default

                ! Not implemented, throw error
                call gdErrorHandler('EvaluateCostFunctionLR: hessian' &
                    // ' not yet implemented for this design variable' &
                    // ' type')

            end select

        end if

        ! Deassociate
        !============
        end associate

    end subroutine

    ! Cost function data writing 
    subroutine WriteCostFunctionDataPRPB(costfunction, grid)

        ! Description
        !============
        ! Write out the cost function data for the LR cost function.
        ! Here, this consists of the vertex pair data in IDn, xn, yn 
        ! format

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionPRPBUDT)      :: costfunction
        type(GridUDT)                   :: grid

        ! Auxiliary
        integer(I8)                     :: ncol, nrow 

        integer(I8), allocatable        :: IDn(:, :) 
        real(R8), allocatable           :: xn(:, :), yn(:, :)
        character(:), allocatable       :: filename 

        ! Loop
        integer(I8)                     :: j 

        ! Initialize
        !===========
        ! Set filename
        allocate(character(len('costfunction_psipairs_PRPB')) :: filename)
        filename = 'costfunction_psipairs_PRPB'

        ! Allocate
        nrow = size(costfunction%psipairs, 1)
        ncol = size(costfunction%psipairs, 2)
        allocate(IDn(nrow, ncol), xn(nrow, ncol), yn(nrow, ncol))

        ! Unpack
        associate(&
            vpairs      => costfunction%psipairs,         &
            x           => grid%vert%x,                 &
            y           => grid%vert%y)

        ! Loop
        do j = 1, ncol 
            IDn(:, j) = vpairs(:, j) 
            xn(:, j) = x(vpairs(:, j)) 
            yn(:, j) = y(vpairs(:, j)) 
        end do

        ! Call writer
        !============
        call WriteVertexPairData(IDn, xn, yn, filename)

        ! Housekeeping
        !=============
        end associate
        deallocate(IDn, xn, yn)
        


    end subroutine

    ! Housekeeping
    subroutine AllocateCostFunctionPRPB(costfunction, npp)

        ! Description
        !============
        ! Allocate, assumed that costfunction%nvpairs is given

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionPRPBUDT)      :: costfunction
        integer(I8)                     :: npp

        ! Allocate
        !=========
        allocate(costfunction%psipairs(npp, 3))
        allocate(costfunction%b0(npp))
        allocate(costfunction%wt(npp))

    end subroutine

    subroutine DeallocateCostFunctionPRPB(costfunction)

        ! Description
        !============
        ! Deallocate

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionPRPBUDT)      :: costfunction

        ! Deallocate
        !===========
        if (allocated(costfunction%psipairs)) then 
            deallocate(costfunction%psipairs)
            deallocate(costfunction%b0)
            deallocate(costfunction%wt)
        end if

    end subroutine

    subroutine DestroyCostFunctionPRPB(costfunction)

        ! Description
        !============
        ! Deallocate

        ! Declare variables
        !==================
        ! Arguments
        type(CostfunctionPRPBUDT)      :: costfunction

        ! Destroy
        !========
        call costfunction%Deallocate()

    end subroutine

    !------------------------------------------------------------------!
    !                       PSI RATIO, PSI BASED 2                     !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCostFunctionPRPB2(costfunction, grid, &
        magneticField, environment, options)

        ! Description
        !============
        ! Initialize the cost function and its parameters based on the 
        ! grid, magnetic field, and environment structures. 

        ! Simply call the initialization of the original lenght ratio
        ! cost function. 

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionPRPB2UDT)         :: costfunction
        type(GridUDT)                       :: grid
        type(MagneticFieldUDT)              :: magneticField 
        type(EnvironmentUDT)                :: environment 
        type(CostFunctionOptionsUDT)        :: options
        
        ! Initialize
        !===========
        call costfunction%cfv_prpb%Initialize(grid, magneticField, &
            environment, options)

        ! (Re)set the scaling constant
        costfunction%cfv_prpb%lambda = options%PRPB%lambda

    end subroutine

    ! True initialization
    subroutine FinalizeInitializationCostFunctionPRPB2(costfunction, &
        designvariables, grid, magneticField, environment, options)

        ! Description
        !============
        ! This is the actual initialization for the PRPB cost function.
        ! Here, we check whether the design variables are compatible 
        ! with the cost function. 

        ! Declare
        !========
        ! Arguments
        class(CostfunctionPRPB2UDT)     :: costfunction
        class(DesignVariablesGDUDT)     :: designvariables
        type(GridUDT)                   :: grid 
        type(MagneticFieldUDT)          :: magneticField
        type(EnvironmentUDT)            :: environment
        type(CostfunctionOptionsUDT)    :: options

        ! Call subroutine
        !================
        call costfunction%cfv_prpb%FinalizeInitialization(designvariables, &
            grid, magneticField, environment, options)

    end subroutine

    ! Cost function evaluation
    subroutine EvaluateCostFunctionPRPB2(costfunction, J, gradJ, hessJ, &
        grid, magneticField, environment, dogradient, dohessian, &
        designvariables)

        ! Description
        !============
        ! Evaluate the cost function, the gradient and its hessian. 
        ! Here, we simply call the same cost function twice, but switch
        ! the order of the indices and recompute the bias. 

        ! Notes:
        !=======
        ! Possible future performance improvements:
        ! - Allocating hessian stuff only once and storing indices, 
        ! since they don't change
        ! - Instead of recomputing auxiliary variables, store them. May
        ! not actually be better in terms of computational time, but 
        ! may lead to shorter and hence better maintainable code. 

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionPRPB2UDT)     :: costfunction
        real(R8)                        :: J, J1, J2
        real(R8), allocatable           :: gradJ(:), gradJ1(:), &
            gradJ2(:) 
        type(MySparseUDT)               :: hessJ, hessJ1, hessJ2 
        type(GridUDT)                   :: grid 
        type(MagneticFieldUDT)          :: magneticField 
        type(EnvironmentUDT)            :: environment
        logical                         :: dogradient, dohessian 
        class(DesignVariablesGDUDT)     :: designvariables

        ! Loop variables

        ! Auxiliary
        integer(I8), allocatable        :: tempvpairs(:,:)
        real(R8), allocatable           :: tempb0(:) 
                                        
        ! Initialize
        !===========
        ! Store original vertex pairs and bias
        allocate(tempvpairs, source=costfunction%cfv_prpb%psipairs)
        allocate(tempb0, source=costfunction%cfv_prpb%b0)
        tempvpairs = costfunction%cfv_prpb%psipairs 
        tempb0  = costfunction%cfv_prpb%b0
         
        ! Cost function
        J = 0
        J1 = 0
        J2 = 0

        ! Gradient
        gradJ(:) = 0
        allocate(gradJ1(size(gradJ)), gradJ2(size(gradJ)))

        ! Hessian
        hessJ1%nrow = hessJ%nrow 
        hessJ2%nrow = hessJ%nrow 
        hessJ1%ncol = hessJ%ncol 
        hessJ2%ncol = hessJ%ncol 

        ! Compute cost function
        !======================
        ! First contribution
        call costfunction%cfv_prpb%Evaluate(J1, gradJ1, &
            hessJ1, grid, magneticField, environment, dogradient, &
            dohessian, designvariables)
        
        ! Adjust vertex pairs and bias
        costfunction%cfv_prpb%psipairs(:, :) = tempvpairs(:, [1, 3, 2])
        costfunction%cfv_prpb%b0(:) = 1/tempb0

        ! Second contribution
        call costfunction%cfv_prpb%Evaluate(J2, gradJ2, &
            hessJ2, grid, magneticField, environment, dogradient, &
            dohessian, designvariables)

        ! Reset vertex pairs and bias
        costfunction%cfv_prpb%b0(:) = tempb0
        costfunction%cfv_prpb%psipairs(:, :) = tempvpairs

        ! Add
        J = J1 + J2 
        gradJ = gradJ1 + gradJ2
        hessJ = hessJ1 + hessJ2

        ! Housekeeping
        !=============
        deallocate(gradJ1, gradJ2)

    end subroutine

    ! Data writing
    subroutine WriteCostFunctionDataPRPB2(costfunction, grid)

        ! Description
        !============
        ! This routine calls the data writer of the PRPB cost function

        ! Declare
        !========
        ! Arguments
        class(CostfunctionPRPB2UDT)     :: costfunction
        type(GridUDT)                   :: grid 

        ! Call
        !=====
        call costfunction%cfv_prpb%WriteData(grid)

    end subroutine

    ! Housekeeping
    subroutine AllocateCostFunctionPRPB2(costfunction, npp)

        ! Description
        !============
        ! Allocate, assumed that costfunction%nvpairs is given

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionPRPB2UDT)     :: costfunction
        integer(I8)                     :: npp

        ! Allocate
        !=========
        call costfunction%cfv_prpb%Allocate(npp)

    end subroutine

    subroutine DeallocateCostFunctionPRPB2(costfunction)

        ! Description
        !============
        ! Deallocate

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionPRPB2UDT)     :: costfunction

        ! Deallocate
        !===========
        call costfunction%cfv_prpb%Deallocate()

    end subroutine

    subroutine DestroyCostFunctionPRPB2(costfunction)

        ! Description
        !============
        ! Deallocate

        ! Declare variables
        !==================
        ! Arguments
        type(CostfunctionPRPB2UDT)     :: costfunction

        ! Destroy
        !========
        call costfunction%cfv_prpb%Deallocate()

    end subroutine

    !------------------------------------------------------------------!
    !                             LR-FAD                               !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCostFunctionLRFAD(costfunction, grid, &
        magneticField, environment, options)

        ! Description
        !============
        ! Initialize the cost function and its parameters based on the 
        ! grid, magnetic field, and environment structures. 

        ! Simply call the initialization of the original lenght ratio
        ! cost function. 

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionLRFADUDT)         :: costfunction
        type(GridUDT)                       :: grid
        type(MagneticFieldUDT)              :: magneticField 
        type(EnvironmentUDT)                :: environment 
        type(CostFunctionOptionsUDT)        :: options
        
        ! Initialize
        !===========
        call costfunction%cfv_lr%Initialize(grid, magneticField, &
            environment, options)
        call costfunction%cfv_fad%Initialize(grid, magneticField, &
            environment, options)

        ! (Re)set the scaling constant
        !costfunction%cfv_lr%lambda = 1e4 ! seems to agree well with most grids
        !costfunction%cfv_fad%lambda = 1e3

    end subroutine

    ! Cost function evaluation
    subroutine EvaluateCostFunctionLRFAD(costfunction, J, gradJ, hessJ, &
        grid, magneticField, environment, dogradient, dohessian, &
        designvariables)

        ! Description
        !============
        ! Evaluate the cost function, the gradient and its hessian. 
        ! Here, we simply call the same cost function twice, but switch
        ! the order of the indices and recompute the bias. 

        ! Notes:
        !=======
        ! Possible future performance improvements:
        ! - Allocating hessian stuff only once and storing indices, 
        ! since they don't change
        ! - Instead of recomputing auxiliary variables, store them. May
        ! not actually be better in terms of computational time, but 
        ! may lead to shorter and hence better maintainable code. 

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionLRFADUDT)     :: costfunction
        real(R8)                        :: J, J1, J2
        real(R8), allocatable           :: gradJ(:), gradJ1(:), &
            gradJ2(:) 
        type(MySparseUDT)               :: hessJ, hessJ1, hessJ2 
        type(GridUDT)                   :: grid 
        type(MagneticFieldUDT)          :: magneticField 
        type(EnvironmentUDT)            :: environment
        logical                         :: dogradient, dohessian 
        class(DesignVariablesGDUDT)     :: designvariables

        ! Loop variables

        ! Auxiliary
                                        
        ! Initialize
        !===========
        ! Cost function
        J = 0
        J1 = 0
        J2 = 0

        ! Gradient
        gradJ(:) = 0
        allocate(gradJ1(size(gradJ)), gradJ2(size(gradJ)))

        ! Hessian
        hessJ1%nrow = hessJ%nrow 
        hessJ2%nrow = hessJ%nrow 
        hessJ1%ncol = hessJ%ncol 
        hessJ2%ncol = hessJ%ncol 

        ! Compute cost function
        !======================
        ! First contribution
        call costfunction%cfv_lr%Evaluate(J1, gradJ1, &
            hessJ1, grid, magneticField, environment, dogradient, &
            dohessian, designvariables)
        
        ! Second contribution
        call costfunction%cfv_fad%Evaluate(J2, gradJ2, &
            hessJ2, grid, magneticField, environment, dogradient, &
            dohessian, designvariables)

        ! Add
        J = J1 + J2 
        gradJ = gradJ1 + gradJ2
        hessJ = hessJ1 + hessJ2

        ! Housekeeping
        !=============
        deallocate(gradJ1, gradJ2)

    end subroutine

    ! Housekeeping
    subroutine AllocateCostFunctionLRFAD(costfunction)

        ! Description
        !============
        ! Dummy function, nothing to be done here

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionLRFADUDT)     :: costfunction

    end subroutine

    subroutine DeallocateCostFunctionLRFAD(costfunction)

        ! Description
        !============
        ! Deallocate

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionLRFADUDT)       :: costfunction

        ! Deallocate
        !===========
        call costfunction%cfv_lr%Deallocate()
        call costfunction%cfv_fad%Deallocate()

    end subroutine

    subroutine DestroyCostFunctionLRFAD(costfunction)

        ! Description
        !============
        ! Deallocate

        ! Declare variables
        !==================
        ! Arguments
        type(CostfunctionLRFADUDT)       :: costfunction

        ! Destroy
        !========
        call costfunction%cfv_lr%Deallocate()
        call costfunction%cfv_fad%Deallocate()

    end subroutine

    !------------------------------------------------------------------!
    !                             GENERAL                              !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCostFunctionGeneral(costfunction, grid, &
        magneticField, environment, options)

        ! Description
        !============
        ! Initialize the cost function and its parameters based on the 
        ! grid, magnetic field, and environment structures. 

        ! Simply call the initialization of the original lenght ratio
        ! cost function. 

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionGeneralUDT)       :: costfunction
        type(GridUDT)                       :: grid
        type(MagneticFieldUDT)              :: magneticField 
        type(EnvironmentUDT)                :: environment 
        type(CostFunctionOptionsUDT)        :: options
        
        ! Initialize
        !===========
        ! Set evaluation switches
        costfunction%doLR       = .false.
        costfunction%doFAD      = .false.
        costfunction%doFA       = .false.
        costfunction%doPRPB     = .false. 
        costfunction%doLRrad    = .false.

        ! Check based on cost function type
        select case (costfunction%type)

        case ('LR_FAD')

            ! Set lambda of other contributions 
            options%FA%lambda       = -1
            options%PRPB%lambda     = -1
            options%LRrad%lambda    = -1
            
        case ('LR_FAD_FA')

            ! Set lambda of other contributions
            options%PRPB%lambda     = -1
            options%LRrad%lambda    = -1

        case ('LR_FAD_PRPB')

            ! Set lambda of other contributions
            options%FA%lambda       = -1
            options%LRrad%lambda    = -1

        case ('LR_FAD_PRPB_FA')

            ! Set lambda of other contributions
            options%LRrad%lambda    = -1

        case ('LR_FAD_PRPB_LRrad')

            ! Set lambda of other contributions
            options%FA%lambda    = -1

        case ('LR_FAD_PRPB_LRrad_FA')

            ! Set lambda of other contributions

        case ('general')

            ! Include all contributions

        case default 

            ! Throw error
            call gdErrorHandler('Unknown cost function type')

        end select

        ! Initialize if necessary
        if (options%LR%lambda > 0) then
            costfunction%doLR = .true.
            call costfunction%cfv_lr%Initialize(grid, magneticField, &
                environment, options)
        end if 
        if (options%FAD%lambda > 0) then 
            costfunction%doFAD = .true.
            call costfunction%cfv_fad%Initialize(grid, magneticField, &
                environment, options)
        end if 
        if (options%FA%lambda > 0) then 
            costfunction%doFA = .true.
            call costfunction%cfv_fa%Initialize(grid, magneticField, &
                environment, options)
        end if
        if (options%PRPB%lambda > 0) then 
            costfunction%doPRPB = .true.
            call costfunction%cfv_prpb%Initialize(grid, magneticField, &
                environment, options)
        end if
        if (options%LRrad%lambda > 0) then 
            costfunction%doLRrad = .true.
            call costfunction%cfv_lrrad%Initialize(grid, magneticField, &
                environment, options)
        end if

    end subroutine

    ! Cost function evaluation
    subroutine EvaluateCostFunctionGeneral(costfunction, J, gradJ, hessJ, &
        grid, magneticField, environment, dogradient, dohessian, &
        designvariables)

        ! Description
        !============
        ! Evaluate the cost function, the gradient and its hessian. 
        ! Here, we simply call the same cost function twice, but switch
        ! the order of the indices and recompute the bias. 

        ! Notes:
        !=======
        ! Possible future performance improvements:
        ! - Allocating hessian stuff only once and storing indices, 
        ! since they don't change
        ! - Instead of recomputing auxiliary variables, store them. May
        ! not actually be better in terms of computational time, but 
        ! may lead to shorter and hence better maintainable code. 

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionGeneralUDT)       :: costfunction
        real(R8)                        :: J, Jtemp
        real(R8), allocatable           :: gradJ(:), gradJtemp(:)
        type(MySparseUDT)               :: hessJ, hessJtemp
        type(GridUDT)                   :: grid 
        type(MagneticFieldUDT)          :: magneticField 
        type(EnvironmentUDT)            :: environment
        logical                         :: dogradient, dohessian 
        class(DesignVariablesGDUDT)     :: designvariables

        ! Loop variables

        ! Auxiliary
                                        
        ! Initialize
        !===========
        ! Cost function
        J = 0
        Jtemp = 0

        ! Gradient
        gradJ(:) = 0
        allocate(gradJtemp(size(gradJ)))

        ! Hessian
        hessJtemp%nrow = hessJ%nrow 
        hessJtemp%ncol = hessJ%ncol 
        
        ! Allocate initially to avoid errors 
        if (.not. allocated(hessJ%row)) then 
            hessJ%nval = 0
            call hessJ%Allocate()
        else 
            ! Reset hessian
            call hessJ%Deallocate()
            hessJ%nval = 0
            call hessJ%Allocate()
        end if 
            

        ! Compute cost function
        !======================
        ! Length ratio
        if (costfunction%doLR) then 
            ! Compute
            call costfunction%cfv_lr%Evaluate(Jtemp, gradJtemp, &
                hessJtemp, grid, magneticField, environment, dogradient, &
                dohessian, designvariables)
            
            ! Add
            J       = J + Jtemp 
            gradJ   = gradJ + gradJtemp
            hessJ   = hessJ + hessJtemp

            ! Deallocate
            call hessJtemp%Deallocate()
        end if 
        
        ! Face angle difference
        if (costfunction%doFAD) then 
            ! Compute
            call costfunction%cfv_fad%Evaluate(Jtemp, gradJtemp, &
                hessJtemp, grid, magneticField, environment, dogradient, &
                dohessian, designvariables)

            ! Add
            J       = J + Jtemp 
            gradJ   = gradJ + gradJtemp
            hessJ   = hessJ + hessJtemp

            ! Deallocate
            call hessJtemp%Deallocate()
        end if 

        ! Face angle
        if (costfunction%doFA) then 
            ! Compute
            call costfunction%cfv_fa%Evaluate(Jtemp, gradJtemp, &
                hessJtemp, grid, magneticField, environment, dogradient, &
                dohessian, designvariables)

            ! Add
            J       = J + Jtemp 
            gradJ   = gradJ + gradJtemp
            hessJ   = hessJ + hessJtemp

            ! Deallocate
            call hessJtemp%Deallocate()
        end if 

        ! Psi ratio, psi based
        if (costfunction%doPRPB) then 
            ! Compute
            call costfunction%cfv_prpb%Evaluate(Jtemp, gradJtemp, &
                hessJtemp, grid, magneticField, environment, dogradient, &
                dohessian, designvariables)

            ! Add
            J       = J + Jtemp 
            gradJ   = gradJ + gradJtemp
            hessJ   = hessJ + hessJtemp

            ! Deallocate
            call hessJtemp%Deallocate()
        end if 

        ! Length ratio, radial
        if (costfunction%doLRrad) then 
            ! Compute
            call costfunction%cfv_lrrad%Evaluate(Jtemp, gradJtemp, &
                hessJtemp, grid, magneticField, environment, dogradient, &
                dohessian, designvariables)

            ! Add
            J       = J + Jtemp 
            gradJ   = gradJ + gradJtemp
            hessJ   = hessJ + hessJtemp

            ! Deallocate
            call hessJtemp%Deallocate()
        end if 

        ! Housekeeping
        !=============
        deallocate(gradJtemp)

    end subroutine

    ! Housekeeping
    subroutine AllocateCostFunctionGeneral(costfunction)

        ! Description
        !============
        ! Dummy function, nothing to be done here

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionGeneralUDT)       :: costfunction

    end subroutine

    subroutine DeallocateCostFunctionGeneral(costfunction)

        ! Description
        !============
        ! Deallocate

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionGeneralUDT)       :: costfunction

        ! Deallocate
        !===========
        call costfunction%cfv_lr%Deallocate()
        call costfunction%cfv_fad%Deallocate()
        call costfunction%cfv_fa%Deallocate()
        call costfunction%cfv_prpb%Deallocate()
        call costfunction%cfv_lrrad%Deallocate()

    end subroutine

    subroutine DestroyCostFunctionGeneral(costfunction)

        ! Description
        !============
        ! Deallocate

        ! Declare variables
        !==================
        ! Arguments
        type(CostfunctionGeneralUDT)       :: costfunction

        ! Destroy
        !========
        call costfunction%cfv_lr%Deallocate()
        call costfunction%cfv_fad%Deallocate()
        call costfunction%cfv_fa%Deallocate()
        call costfunction%cfv_prpb%Deallocate()
        call costfunction%cfv_lrrad%Deallocate()

    end subroutine
    

end module