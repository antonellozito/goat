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
    use gdmod_types 
    use optmod_costfunction

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

        ! The following routines should be implemented for these cost
        ! functions (see also the interface below for a description of
        ! what the routines should do):
        ! - Initialize
        ! - Evaluate

        ! Cost function value
        real(R8)        :: J 

    contains

        ! Cost function initialization
        procedure(InitializeCostfunctionINT), deferred :: Initialize

        ! Cost function evaluation
        !procedure(EvaluateCostFunctionINT), deferred :: Evaluate

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
    
        ! Notes
        !======
        ! Note 1: it is assumed that each vertex only has one vertex 
        ! pair 

        ! Fields
        real(R8)                    :: lambda ! scaling constant
        real(R8), allocatable       :: b0(:) ! desired length ratio per vertex
        real(R8), allocatable       :: wt(:) ! weight factor per vertex
        integer(I8), allocatable    :: vpairs(:, :), nvpairs(:) ! vertex pairs

    contains

        ! Initialization
        procedure :: Initialize             => InitializeCostfunctionLR

        ! Housekeeping
        procedure :: Allocate               => AllocateCostfunctionLR
        procedure :: Deallocate             => DeallocateCostFunctionLR
        final :: DestroyCostFunctionLR

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
            magneticField, environment)

            ! Description
            !============
            ! This routine should initialize all additional parameters
            ! that are needed to evaluate the cost function (e.g. the 
            ! vertex indices where the cost function is defined).
            
            ! Import
            import :: CostfunctionGDUDT, GridUDT, MagneticFieldUDT, &
                EnvironmentUDT

            ! Declare 
            class(CostfunctionGDUDT)        :: costfunction 
            type(GridUDT)                   :: grid
            type(MagneticFieldUDT)          :: magneticField 
            type(EnvironmentUDT)            :: environment

        end subroutine

        ! Cost function evaluation
        !subroutine EvaluateCostFunctionINT(costfunction, grid, &
         !   magneticField, environment, dogradient, dohessian)

    end interface

    contains

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    ! Length ratio cost function
    !===========================
    ! Initialization
    subroutine InitializeCostFunctionLR(costfunction, grid, &
        magneticField, environment)

        ! Description
        !============
        ! Initialize the cost function and its parameters based on the 
        ! grid, magnetic field, and environment structures. Here, the 
        ! length ratio cost function is initialized, which requires

        ! Modules
        use BicubicSplineInterpolant

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionLRUDT)            :: costfunction
        type(GridUDT)                       :: grid
        type(MagneticFieldUDT)              :: magneticField 
        type(EnvironmentUDT)                :: environment 

        ! Loop variables
        integer(I8)                         :: i, j

        ! Auxiliary variables
        type(VertexUDT)                       :: vert

        integer                             :: sgn2, sgn3
        integer(I8)                         :: sp, ep, tID, v2, v3
        real(R8)                            :: Btx2, Btx3, Bty2, Bty3, &
            dx2, dy2, dx3, dy3

        integer(I8), allocatable            :: tvn(:), temptvn(:)
        real(R8), allocatable               :: Btx(:), Bty(:)
        logical, allocatable                :: cID(:)

        ! Data
        
        ! Initialize
        !===========

        ! Set the scaling constant
        costfunction%lambda = 1e4 ! seems to agree well with most grids

        ! Allocate
        call costfunction%Allocate(grid%vert%ntot)
        allocate(Btx(grid%vert%ntot))
        allocate(Bty(grid%vert%ntot))

        ! Associate some fields
        associate(vert => grid%vert, x => grid%vert%x, y => grid%vert%y, &
            b0 => costfunction%b0, wt => costfunction%wt, &
            nvpairs => costfunction%nvpairs, &
            vpairs => costfunction%vpairs)

        ! Set the initial weighting factors
        wt(:) = 1

        ! Compute the magnetic field vectors at the vertex locations
        call EvaluateBicubicSplineInterpolant(x, y, Btx, &
            magneticField%interp, '0', '1') 
        call EvaluateBicubicSplineInterpolant(x, y, Bty, &
            magneticField%interp, '1', '0') 
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
            allocate(temptvn(vert%neigP(i, 2)))
            allocate(cID(vert%neigP(i, 2)))

            ! Get the vertex neighbours of this vertex
            sp = vert%neigP(i, 1)
            ep = vert%neigP(i, 1) + vert%neigP(i, 2)-1
            temptvn = vert%neiglist(sp:ep)
            
            ! Get the ID of the coordinate line
            tID = vert%fieldlineID(i)
            
            ! Check which vertices have the same ID
            cID = (tID == vert%fieldlineID(temptvn))
            
            ! Extract
            allocate(tvn(count(cID)))
            tvn = pack(temptvn, cID)
            
            ! Assemble
            if (size(tvn) > 0) then
                ! Update counter
                nvpairs(i) = (size(tvn)/2)
                
                ! If multiple pairs, loop
                do j = 1, nvpairs(i)
                    ! Normally, multiple pairs only occur at x-points, and,
                    ! since the coordinates are sorted, the corresponding
                    ! pairs should be tvn(j) and tvn(j+nvpairs(i))
                    
                    ! Get vertices
                    v2 = tvn(j);
                    v3 = tvn(j+nvpairs(i));
                   
                   ! Get vectors
                    dx2 = x(v2) - x(i); dy2 = y(v2) - y(i);
                    dx3 = x(v3) - x(i); dy3 = y(v3) - y(i);
                    
                    ! Compute weight
                    wt(i) = wt(i) + &
                        1/(dx2**2 + dy2**2) + 1/(dx3**2 + dy3**2);
                    
                    ! Check if we're dealing with an x-point
                    if (nvpairs(i) > 1) then
                        ! Here, the gradient *should* vanish. For now, we
                        ! cope with this by setting the desired ratio to 1,
                        ! such that it does not matter which length is
                        ! considered first.
                        b0(i) = 1;
                        sgn2 = -1;
                        wt(i) = 0;
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
                            costfunction%wt(i) = 0
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
            end if

            ! Housekeeping
            deallocate(tvn, temptvn, cID)

        end do

        ! Housekeeping
        deallocate(Btx, Bty)

        ! End associate
        end associate


    end subroutine

    ! Housekeeping
    subroutine AllocateCostFunctionLR(costfunction, nv)

        ! Description
        !============
        ! Allocate, assumed that costfunction%nvpairs is given

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionLRUDT)        :: costfunction
        integer(I8)                     :: nv

        ! Allocate
        !=========
        print *, nv
        allocate(costfunction%nvpairs(nv))
        allocate(costfunction%vpairs(nv,2))
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
        deallocate(costfunction%vpairs)
        deallocate(costfunction%b0)
        deallocate(costfunction%wt)

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

end module