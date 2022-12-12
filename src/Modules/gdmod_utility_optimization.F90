!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains utility routines for the optimization modules of
! the grid deformation framework. These routines include for example 
! routines to compute faces on which orthogonality constraints should be
! imposed, routines that allow to derive the vertices that are x-points,
! and so on. In general, any routine that is necessary to set up the 
! constraints/costfunction, but that is not needed anymore for the 
! remainder of the code, can be added here. 

module gdmod_utility_optimization

    ! Initialize
    !============
    ! Load modules
    !use mod_plotter
    use gdmod_types
    use gdmod_userinput 
    !use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN

    ! The usual
    implicit none
    save
    public 

    contains 

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                            Grid quantities                       !
    !------------------------------------------------------------------!

    ! Determination of X-point(s)
    subroutine DetermineXPoints(xpind, nxpind, order, grid)

        ! Description
        !============
        ! This routine determines the indices (and their number) based 
        ! on the information given in the grid structure. Note that 
        ! the location of the x-point should already be available 
        ! through output of the grid generator, but we provide a routine 
        ! here to recompute the location if it wasn't saved explicitly
        ! in the grid structure. Additionally, we determine the order of
        ! the x-point (see algorithm description below)

        ! Algorithm
        !==========
        ! Basically, we simply check if a vertex has 4 or more 
        ! neighbours with the same flux surface ID (which is non-zero)
        ! as the current vertex. The order is then simply determined by 
        ! the number of separatrix legs, i.e. 
        !
        !       o = n/2 - 1,
        !
        ! where n is the number of vertices with the same ID and o is 
        ! the order. 

        ! Initialize
        !===========
        ! The usual
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        integer(I8)                         :: nxpind
        integer(I8), allocatable            :: xpind(:), order(:)
        type(GridUDT), intent(in)           :: grid

        ! Loop variables
        integer(I8)                         :: i

        ! Auxiliary variables 
        integer(I8)                         :: tfID, ncIDs
        integer(I8), allocatable            :: temporder(:), &
            tempxpind(:), tvn(:), tvnfID(:)

        ! Data

        ! Initialize
        !===========
        ! Associate
        associate(&
            fID         => grid%vert%fieldlineID,   &
            vert        => grid%vert)

        ! Check allocation - shouldn't be the case as size unknown a 
        ! priori
        if (allocated(xpind)) then
            ! Deallocate
            deallocate(xpind)
        end if
        if (allocated(order)) then 
            ! Deallocate
            deallocate(order)
        end if

        ! Allocate temporary arrays (too big, trim later)
        nxpind = 0
        allocate(tempxpind(vert%ntot))
        allocate(temporder(vert%ntot))

        ! Determine x-points
        !===================
        ! Loop over all vertices
        do i = 1, vert%ntot 
            ! Get the current field line ID, skip if zero
            tfID = fID(i)
            if (tfID == 0) then
                cycle ! skip rest of the loop for this index
            end if

            ! Get the vertex neighbours
            allocate(tvn(vert%neigP(i, 2)))
            tvn = vert%neiglist(vert%neigP(i, 1):&
                (vert%neigP(i, 1) + vert%neigP(i, 2)-1))

            ! Get their IDs
            allocate(tvnfID(vert%neigP(i, 2)))
            tvnfID = fID(tvn)

            ! Get the number of common IDs
            ncIDs = count(tvnfID == tfID)

            ! If equal or larger than 4, add as x-point
            if (ncIDs >= 4) then 
                ! Check if it is a multiple of 2
                if (modulo(ncIDs, 2) .ne. 0) then 
                    ! Uneven number, throw error
                    call gdErrorHandler('DetermineXPoints: '&
                    // 'supposed x-point has uneven number' &
                    // 'of vertices with same ID, check grid' &
                    // ' consistency')
                end if

                ! Update the counter
                nxpind = nxpind + 1

                ! Add the x-point
                tempxpind(nxpind) = i 

                ! Compute the order
                temporder(nxpind) = ncIDs/2 - 1

            end if
            
            ! Deallocate
            deallocate(tvn, tvnfID)

        end do

        ! Output
        !=======
        ! Allocate
        allocate(xpind(nxpind), order(nxpind))

        ! Set output
        xpind = tempxpind(1:nxpind)
        order = temporder(1:nxpind)

        ! Housekeeping
        !=============
        ! Deallocate 
        deallocate(tempxpind, temporder)

        ! Deassociate
        end associate


    end subroutine

end module