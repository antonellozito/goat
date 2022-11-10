subroutine ExtractVesselData(vessel, vesseloptions)

    ! Description
    !============
    ! Process the vessel structure to provide the necessary data for the
    ! grid deformation. Here, we mainly construct the vessel edges and
    ! indicate the vessel segments. We also store the (unique) vessel 
    ! coordinates. The actual vessel data should be read in already. 

    ! Notes
    !======

    ! Initialize
    !===========
    ! Declare modules
    use gdmod_types 
    use gdmod_interfaces
    use gdmod_plots
    use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    type(VesselUDT), intent(inout)      :: vessel
    type(VesselOptionsUDT), intent(in)  :: vesseloptions

    ! Loop variables
    integer(I8)                         :: i, k 

    ! Auxiliary variables 
    integer(I8)                         :: pec, nc, netot, nv, thisne 
    integer(I8), allocatable            :: polygonstarts(:), pe(:, :)
    real(R8), allocatable               :: xv(:), yv(:)

    ! Plotting

    ! Vessel edges
    !=============
    ! Allocate
    allocate(polygonstarts(vessel%nstructures))

    ! Counters
    pec = 1 ! polygon edge counter
    nc = 1 ! 

    ! Compute the final and total number of edges
    netot = 0
    nv = 0
    do i = 1, vessel%nstructures
        ! Compute the number of edges 
        thisne = vessel%structures(i)%np-1 
        if (vessel%structures(i)%isclosed) then 
            thisne = thisne + 1 
        end if

        ! Compute the number of points
        nv = nv + vessel%structures(i)%np

        ! Set the polygon start
        polygonstarts(i) = pec

        ! Update the counters
        pec = pec + thisne
        netot = netot + thisne
    end do

    ! Allocate
    allocate(xv(nv), yv(nv), pe(netot, 2))

    ! Reset counters
    pec = 1 
    nc = 1  

    ! Loop over all structures
    do i = 1, vessel%nstructures
        ! Compute the current number of edges 
        thisne = vessel%structures(i)%np-1 
        if (vessel%structures(i)%isclosed) then 
            thisne = thisne + 1 
        end if
        
        ! Set the edges
        pe(pec:pec+thisne-1, 1) = [(k, k = nc, nc+thisne-1)]
        pe(pec:pec+thisne-1, 2) = [(k, k = nc+1, nc+thisne)]
        if (vessel%structures(i)%isclosed) then 
            ! Adjust the last point of the last edge to be the first 
            ! point
            pe(pec+thisne-1, 2) = pe(pec, 1)
        end if
        
        ! Add the coordinates
        xv(nc:nc+vessel%structures(i)%np-1) = vessel%structures(i)%x 
        yv(nc:nc+vessel%structures(i)%np-1) = vessel%structures(i)%y

        ! Update counters
        pec = pec + thisne 
        nc = nc + vessel%structures(i)%np

    end do 

    ! Allocate vessel
    if (allocated(vessel%x)) then 
        ! Deallocate
        print *, 'ExtractVesselData: reinitializing vessel coordinates'
        deallocate(vessel%x, vessel%y)
    end if
    allocate(vessel%x(nv), vessel%y(nv))
    allocate(vessel%edges(netot, 2), &
        vessel%polygonstart(vessel%nstructures))

    ! Assign to vessel structure
    vessel%x            = xv 
    vessel%y            = yv 
    vessel%nv           = nv
    vessel%edges        = pe 
    vessel%nedges       = size(pe, 1)
    vessel%polygonstart = polygonstarts 
    vessel%np           = vessel%nstructures
    
    ! Deallocate
    deallocate(polygonstarts, pe, xv, yv)

end subroutine