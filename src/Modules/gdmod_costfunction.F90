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
    use optmod_costfunction
    use gdmod_plots

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
        real(R8), allocatable       :: b0(:) ! desired length ratio per vertex
        real(R8), allocatable       :: wt(:) ! weight factor per vertex
        integer(I8), allocatable    :: vpairs(:, :), nvpairs(:) ! vertex pairs

    contains

        ! Initialization
        procedure :: Initialize             => InitializeCostfunctionLR

        ! Evaluation
        procedure :: Evaluate               => EvaluateCostFunctionLR

        ! Housekeeping
        procedure :: Allocate               => AllocateCostfunctionLR
        procedure :: Deallocate             => DeallocateCostFunctionLR
        final :: DestroyCostFunctionLR

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

        ! Housekeeping
        procedure :: Allocate               => AllocateCostFunctionFAD
        procedure :: Deallocate             => DeallocateCostFunctionFAD
        final :: DestroyCostFunctionFAD

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
        call costfunction%Allocate(grid%vert%ntot, 4)
        allocate(Btx(grid%vert%ntot))
        allocate(Bty(grid%vert%ntot))

        ! Associate some fields
        associate(vert => grid%vert, x => grid%vert%x, y => grid%vert%y, &
            b0 => costfunction%b0, wt => costfunction%wt, &
            nvpairs => costfunction%nvpairs, &
            vpairs => costfunction%vpairs)

        ! Set the initial weighting factors
        wt(:) = 1

        ! Initialize
        vpairs(:, :) = 0
        nvpairs(:) = 0

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
                    v2 = tvn(j)
                    v3 = tvn(j+nvpairs(i))
                   
                   ! Get vectors
                    dx2 = x(v2) - x(i); dy2 = y(v2) - y(i)
                    dx3 = x(v3) - x(i); dy3 = y(v3) - y(i)
                    
                    ! Compute weight
                    wt(i) = wt(i) + &
                        1/(dx2**2 + dy2**2) + 1/(dx3**2 + dy3**2)
                    
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
            end if

            ! Housekeeping
            deallocate(tvn, temptvn, cID)

        end do

        ! Housekeeping
        deallocate(Btx, Bty)

        ! End associate
        end associate


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
        real(R8)                        :: x1, x2, x3, y1, y2, y3
        real(R8)                        :: dx1, dx2, dy1, dy2, d1, d2
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
            lambda  => costfunction%lambda )

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
                J = J + 0.5*(d1/d2 - b0(i))**2

            end do 
        end do

        ! Scale
        J = lambda*J

        ! Compute gradient
        !=================
        if (dogradient) then 

            ! Check the design variables
            select case (trim(designvariables%type))

            case ('coordinates')

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
                        gradJ(v1) = gradJ(v1) + & 
                            (b0(i) - d1/d2) * &
                            (dx1/(d1*d2) - (d1*dx2)/d2**3)
                        gradJ(v2) = gradJ(v2) + &
                            -(dx1*(b0(i) - d1/d2))/(d1*d2)
                        gradJ(v3) = gradJ(v3) + &
                            (d1*dx2*(b0(i) - d1/d2))/d2**3
                        
                        gradJ(v1+vert%ntot) = gradJ(v1+vert%ntot) + &
                            (b0(i) - d1/d2) * &
                            (dy1/(d1*d2) - (d1*dy2)/d2**3)
                        gradJ(v2+vert%ntot) = gradJ(v2+vert%ntot) + &
                            -(dy1*(b0(i) - d1/d2))/(d1*d2)
                        gradJ(v3+vert%ntot) = gradJ(v3+vert%ntot) + &
                            (d1*dy2*(b0(i) - d1/d2))/d2**3

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
        if (dohessian) then

            ! Check the design variables
            select case (trim(designvariables%type))

            case ('coordinates')

                ! Allocate the hessian (if not already done so)
                if (.not. allocated(hessJ%row)) then
                    ! Allocate the sparse matrix
                    hessJ%nval = 36*sum(nvpairs) ! this should be exact and constant
                    call hessJ%Allocate()
                end if

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

                        ! Compute Hessian contributions - ordened per
                        ! vertex pair (e.g. v1, v1), split up in xx, xy,
                        ! yx, yy. 

                        ! d2J/dx1**2, d2J/dy1**2, d2J/dx1dy1, d2J/dy1dx1
                        row(cc) = v1; col(cc) = v1
                        valxx(cc) = (b0(i) - d1/d2) & 
                            * (d1/d2**3 - 1/(d1*d2) - & 
                            (3*d1*dx2**2)/d2**5 + dx1**2/(d1**3*d2) & 
                            + (2*dx1*dx2)/(d1*d2**3)) + (dx1/(d1*d2) & 
                            - (d1*dx2)/d2**3)**2
                        valyy(cc) = (b0(i) - d1/d2) &
                            * (d1/d2**3 - 1/(d1*d2) - & 
                            (3*d1*dy2**2)/d2**5 + dy1**2/(d1**3*d2) & 
                            + (2*dy1*dy2)/(d1*d2**3)) + & 
                            (dy1/(d1*d2) - (d1*dy2)/d2**3)**2
                        valxy(cc) = (b0(i) - d1/d2) & 
                            *((dx1*dy1)/(d1**3*d2) - & 
                            (3*d1*dx2*dy2)/d2**5 + & 
                            (dx1*dy2)/(d1*d2**3) + & 
                            (dx2*dy1)/(d1*d2**3)) + (dx1/(d1*d2) & 
                            - (d1*dx2)/d2**3)*(dy1/(d1*d2) &
                            - (d1*dy2)/d2**3)
                        valyx(cc) = valxy(cc)
                        cc = cc+1
                        
                        ! d2J/dx1dx2, d2J/dy1dy2, d2J/dx1dy2, d2J/dy1dx2
                        row(cc) = v1; col(cc) = v2
                        valxx(cc) = - (b0(i) - d1/d2) &
                            * (dx1**2/(d1**3*d2) - 1/(d1*d2) + & 
                            (dx1*dx2)/(d1*d2**3)) - (dx1*(dx1/(d1*d2) - & 
                            (d1*dx2)/d2**3))/(d1*d2) !x1x2
                        valyy(cc) = - (b0(i) - d1/d2)&
                            *(dy1**2/(d1**3*d2) - 1/(d1*d2) &
                            + (dy1*dy2)/(d1*d2**3)) - &
                            (dy1*(dy1/(d1*d2) - &
                            (d1*dy2)/d2**3))/(d1*d2) !y1y2
                        valxy(cc) = - (b0(i) - d1/d2)&
                            *((dx1*dy1)/(d1**3*d2) + &
                            (dx2*dy1)/(d1*d2**3)) - &
                            (dy1*(dx1/(d1*d2) - &
                            (d1*dx2)/d2**3))/(d1*d2) !x1y2
                        valyx(cc) = - (b0(i) - d1/d2) &
                            *((dx1*dy1)/(d1**3*d2) + &
                            (dx1*dy2)/(d1*d2**3)) - &
                            (dx1*(dy1/(d1*d2) - &
                            (d1*dy2)/d2**3))/(d1*d2) !y1x2
                        cc = cc+1
                        
                        row(cc) = v2; col(cc) = v1
                        valxx(cc) = valxx(cc-1) !x2x1
                        valyy(cc) = valyy(cc-1) !y2y1
                        valxy(cc) = valyx(cc-1) !x1y2
                        valyx(cc) = valxy(cc-1) !y2x1
                        cc = cc+1
                        
                        row(cc) = v2; col(cc) = v2
                        valxx(cc) = dx1**2/(d1**2*d2**2) - &
                            (b0(i) - d1/d2)/(d1*d2) + &
                            (dx1**2*(b0(i) - d1/d2))/(d1**3*d2)
                        valyy(cc) = dy1**2/(d1**2*d2**2) - &
                            (b0(i) - d1/d2)/(d1*d2) + &
                            (dy1**2*(b0(i) - d1/d2))/(d1**3*d2)
                        valxy(cc) = (dx1*dy1)/(d1**2*d2**2) &
                            + (dx1*dy1*(b0(i) - d1/d2))/(d1**3*d2)
                        valyx(cc) = valxy(cc)
                        cc = cc+1
                        
                        row(cc) = v1; col(cc) = v3
                        valxx(cc) = (d1*dx2*(dx1/(d1*d2) - &
                            (d1*dx2)/d2**3))/d2**3 - &
                            (b0(i) - d1/d2)*(d1/d2**3 - &
                            (3*d1*dx2**2)/d2**5 + (dx1*dx2)/(d1*d2**3))
                        valyy(cc) = (d1*dy2*(dy1/(d1*d2) - &
                            (d1*dy2)/d2**3))/d2**3 - &
                            (b0(i) - d1/d2)*(d1/d2**3 - &
                            (3*d1*dy2**2)/d2**5 + (dy1*dy2)/(d1*d2**3))
                        valxy(cc) = (b0(i) - d1/d2) &
                            *((3*d1*dx2*dy2)/d2**5 - &
                            (dx1*dy2)/(d1*d2**3)) + &
                            (d1*dy2*(dx1/(d1*d2) - &
                            (d1*dx2)/d2**3))/d2**3
                        valyx(cc) = (b0(i) - d1/d2)&
                            *((3*d1*dx2*dy2)/d2**5 - &
                            (dx2*dy1)/(d1*d2**3)) + &
                            (d1*dx2*(dy1/(d1*d2) - &
                            (d1*dy2)/d2**3))/d2**3
                        cc = cc+1
                        
                        row(cc) = v3; col(cc) = v1
                        valxx(cc) = valxx(cc-1)
                        valyy(cc) = valyy(cc-1)
                        valxy(cc) = valyx(cc-1)
                        valyx(cc) = valxy(cc-1)
                        cc = cc+1
                        
                        row(cc) = v2; col(cc) = v3
                        valxx(cc) = (dx1*dx2*(b0(i) - d1/d2)) &
                            /(d1*d2**3) - (dx1*dx2)/d2**4
                        valyy(cc) = (dy1*dy2*(b0(i) - d1/d2)) &
                            /(d1*d2**3) - (dy1*dy2)/d2**4
                        valxy(cc) = (dx1*dy2*(b0(i) - d1/d2)) &
                            /(d1*d2**3) - (dx1*dy2)/d2**4
                        valyx(cc) = (dx2*dy1*(b0(i) - d1/d2)) &
                            /(d1*d2**3) - (dx2*dy1)/d2**4
                        cc = cc+1
                        
                        row(cc) = v3; col(cc) = v2
                        valxx(cc) = valxx(cc-1)
                        valyy(cc) = valyy(cc-1)
                        valxy(cc) = valyx(cc-1)
                        valyx(cc) = valxy(cc-1)
                        cc = cc+1
                        
                        row(cc) = v3; col(cc) = v3
                        valxx(cc) = (dx2**2*(dx1**2 + dy1**2))/d2**6 &
                            + (d1*(b0(i) - d1/d2))/d2**3 - &
                            (3*d1*dx2**2*(b0(i) - d1/d2))/d2**5
                        valyy(cc) = (dy2**2*(dx1**2 + dy1**2))/d2**6 &
                            + (d1*(b0(i) - d1/d2))/d2**3 - &
                            (3*d1*dy2**2*(b0(i) - d1/d2))/d2**5
                        valxy(cc) = (dx2*dy2*(dx1**2 + dy1**2)) &
                            /d2**6 - (3*d1*dx2*dy2*(b0(i) - d1/d2)) &
                            /d2**5
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

    !------------------------------------------------------------------!
    !                        FACE ANGLE DIFFERENCE                     !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCostFunctionFAD(costfunction, grid, &
        magneticField, environment)

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
        use BicubicSplineInterpolant

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionFADUDT)           :: costfunction
        type(GridUDT)                       :: grid
        type(MagneticFieldUDT)              :: magneticField 
        type(EnvironmentUDT)                :: environment 
        
        ! Loop variables
        integer(I8)                         :: i, j, k, vpc

        ! Auxiliary variables
        integer(I8)                         :: tID, vp1(1:2), &
            vp2(1:2), ntvn, ntemptvn, nvp

        real(R8)                            :: dxf, dyf, txf, tyf, tnf,&
            bxf, byf, bnf, dpf

        logical                             :: thischeck

        integer(I8), allocatable            :: tvn(:), temptvn(:), &
            vpairs(:,:), nvpairs(:), reverse(:), tempvpairs(:,:)

        real(R8), allocatable               :: bx(:), by(:), xv(:,:), &
            yv(:,:), xf(:,:), yf(:,:), gxf(:,:), gyf(:,:), dx(:,:), &
            dy(:,:), dotprod(:,:), xplot(:), yplot(:)

        logical, allocatable                :: cID(:), mask(:), &
            isaligned(:)

        ! Data
        real(R8)                            :: epsalignment 
        data epsalignment /0.1/
        
        ! Initialize
        !===========
        ! Set the scaling constant
        costfunction%lambda = 1e-3 ! seems to agree well with most grids

        ! Allocate (initialize too big)
        nvp = grid%vert%ntot*maxval(grid%vert%neigP(:,2)) ! maximal number of pairs
        allocate(vpairs(nvp,4), nvpairs(grid%vert%ntot))
        ! call costfunction%Allocate(grid%vert%ntot, 4)

        ! Compute magnetic field in grid points
        allocate(bx(grid%vert%ntot), by(grid%vert%ntot))
        call EvaluateBicubicSplineInterpolant( &
            grid%vert%x, grid%vert%y, bx, magneticField%interp, '0', '1')
        call EvaluateBicubicSplineInterpolant( &
            grid%vert%x, grid%vert%y, by, magneticField%interp, '1', '0')
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
            temptvn = vert%neiglist(& 
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
            call EvaluateBicubicSplineInterpolant( &
                xf(:,i), yf(:,i), gxf(:,i), magneticField%interp, '1', '0')
            call EvaluateBicubicSplineInterpolant( &
                xf(:,i), yf(:,i), gyf(:,i), magneticField%interp, '0', '1')
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

        ! Assign to cost function
        !========================
        costfunction%vpairs = vpairs(1:vpc,:)
        costfunction%nvpairs = vpc 
        costfunction%wt(:) = 1

        ! Visualize? 
        !allocate(xplot(size(xf)), yplot(size(yf)))
        !xplot = reshape(xf, [size(xf)])
        !yplot = reshape(yf, [size(yf)])
        !call PlotGridWithPoints(grid, xplot, yplot, '-p')
        !deallocate(xplot, yplot)
        

        ! End associate
        end associate

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
        use BicubicSplineInterpolant

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
        call EvaluateBicubicSplineInterpolant(xfv1, yfv1, gxfv1, &
            magneticField%interp, '1', '0')
        call EvaluateBicubicSplineInterpolant(xfv2, yfv2, gxfv2, &
            magneticField%interp, '1', '0')
        call EvaluateBicubicSplineInterpolant(xfv1, yfv1, gyfv1, &
            magneticField%interp, '0', '1')
        call EvaluateBicubicSplineInterpolant(xfv2, yfv2, gyfv2, &
            magneticField%interp, '0', '1')
        
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

            case ('coordinates')

                ! Allocate
                allocate(gxxfv1(np), gxxfv2(np), gxyfv1(np), gxyfv2(np), &
                    gyxfv1(np), gyxfv2(np), gyyfv1(np), gyyfv2(np))
                allocate(gxxxfv1(np), gxxxfv2(np), gxxyfv1(np), &
                    gxxyfv2(np), gxyxfv1(np), gxyxfv2(np), &
                    gxyyfv1(np), gxyyfv2(np), gyxxfv1(np), gyxxfv2(np), &
                    gyxyfv1(np), gyxyfv2(np), gyyyfv1(np), gyyyfv2(np))

                ! Precompute gradient quantities
                call EvaluateBicubicSplineInterpolant(xfv1, yfv1, &
                    gxxfv1, magneticField%interp, '2', '0')
                call EvaluateBicubicSplineInterpolant(xfv1, yfv1, &
                    gxyfv1, magneticField%interp, '1', '1')
                call EvaluateBicubicSplineInterpolant(xfv1, yfv1, &
                    gyyfv1, magneticField%interp, '0', '2')
                gyxfv1 = gxyfv1 ! symmetric, done for ease here

                call EvaluateBicubicSplineInterpolant(xfv2, yfv2, &
                    gxxfv2, magneticField%interp, '2', '0')
                call EvaluateBicubicSplineInterpolant(xfv2, yfv2, &
                    gxyfv2, magneticField%interp, '1', '1')
                call EvaluateBicubicSplineInterpolant(xfv2, yfv2, &
                    gyyfv2, magneticField%interp, '0', '2')
                gyxfv2 = gxyfv2 ! symmetric, done for ease here

                ! Precompute hessian quantities
                call EvaluateBicubicSplineInterpolant(xfv1, yfv1, &
                    gxxxfv1, magneticField%interp, '3', '0')
                call EvaluateBicubicSplineInterpolant(xfv1, yfv1, &
                    gxxyfv1, magneticField%interp, '2', '1')
                call EvaluateBicubicSplineInterpolant(xfv1, yfv1, &
                    gxyyfv1, magneticField%interp, '1', '2')
                call EvaluateBicubicSplineInterpolant(xfv1, yfv1, &
                    gyyyfv1, magneticField%interp, '0', '3')
                gyxyfv1 = gxyyfv1 ! symmetric, repeated for ease
                gyxxfv1 = gxxyfv1 

                call EvaluateBicubicSplineInterpolant(xfv2, yfv2, &
                    gxxxfv2, magneticField%interp, '3', '0')
                call EvaluateBicubicSplineInterpolant(xfv2, yfv2, &
                    gxxyfv2, magneticField%interp, '2', '1')
                call EvaluateBicubicSplineInterpolant(xfv2, yfv2, &
                    gxyyfv2, magneticField%interp, '1', '2')
                call EvaluateBicubicSplineInterpolant(xfv2, yfv2, &
                    gyyyfv2, magneticField%interp, '0', '3')
                gyxyfv2 = gxyyfv2 ! symmetric, repeated for ease
                gyxxfv2 = gxxyfv2 

            case default

                ! Not implemented, throw error
                call gdErrorHandler('EvaluateCostFunctionFAD: gradient' &
                    // ' not yet implemented for this design variable' &
                    // ' type')

            end select 

        elseif (dogradient) then

            ! Allocate
            allocate(gxxfv1(np), gxxfv2(np), gxyfv1(np), gxyfv2(np), &
                gyxfv1(np), gyxfv2(np), gyyfv1(np), gyyfv2(np))

            ! Precompute gradient quantities
            call EvaluateBicubicSplineInterpolant(xfv1, yfv1, &
                gxxfv1, magneticField%interp, '2', '0')
            call EvaluateBicubicSplineInterpolant(xfv1, yfv1, &
                gxyfv1, magneticField%interp, '1', '1')
            call EvaluateBicubicSplineInterpolant(xfv1, yfv1, &
                gyyfv1, magneticField%interp, '0', '2')
            gyxfv1 = gxyfv1 ! symmetric, done for ease here

            call EvaluateBicubicSplineInterpolant(xfv2, yfv2, &
                gxxfv2, magneticField%interp, '2', '0')
            call EvaluateBicubicSplineInterpolant(xfv2, yfv2, &
                gxyfv2, magneticField%interp, '1', '1')
            call EvaluateBicubicSplineInterpolant(xfv2, yfv2, &
                gyyfv2, magneticField%interp, '0', '2')
            gyxfv2 = gxyfv2 ! symmetric, done for ease here

        end if

        ! Compute gradient
        !=================
        if (dogradient) then 

            ! Check the design variables
            select case (trim(designvariables%type))

            case ('coordinates')

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
        if (dohessian) then

            ! Check the design variables
            select case (trim(designvariables%type))

            case ('coordinates')

                ! Allocate the hessian (if not already done so)
                if (.not. allocated(hessJ%row)) then
                    ! Allocate the sparse matrix
                    hessJ%nval = 64*np ! this should be exact and constant
                    call hessJ%Allocate()
                end if

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

            case ('coordinates')

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

            ! Deallocate
            deallocate(gxxfv1 , gxxfv2 , gxyfv1 , gxyfv2 , &
                gyxfv1 , gyxfv2 , gyyfv1 , gyyfv2 )

        end if

        


        ! Deassociate
        end associate
        

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
        deallocate(costfunction%vpairs)
        deallocate(costfunction%wt)

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
    !                           LENGTH RATIO 2                         !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCostFunctionLR2(costfunction, grid, &
        magneticField, environment)

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
        
        ! Initialize
        !===========
        call costfunction%cfv_lr%Initialize(grid, magneticField, &
            environment)

        ! (Re)set the scaling constant
        costfunction%cfv_lr%lambda = 1e4 ! seems to agree well with most grids

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

    

end module