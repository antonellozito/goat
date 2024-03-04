subroutine ExtractVesselData(vessel, vesseloptions)

    ! Description
    !============
    ! Process the vessel structure to provide the necessary data for the
    ! grid deformation. Here, we construct a polygon set representation 
    ! of the vessel to easier manipulate the vessel structure. It is 
    ! assumed that the vessel structure contains segments that 
    ! properly intersect each other such that a single closed polygon
    ! can be extracted (nested polygons are supported as long as they 
    ! are simply nested, i.e. the interior of the domain can be 
    ! distinguished). 

    ! Notes
    !======

    ! Initialize
    !===========
    ! Declare modules
    use gdmod_types 
    use gdmod_userinput
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
    integer(I8)                         :: i

    ! Auxiliary variables 
    type(PolygonSetUDT)                 :: tempps 
    integer(I8)                         :: nvp, cc, nexcl   
!    integer(I8), allocatable            :: polygonstarts(:)
    real(R8), allocatable               :: tempx(:), &
        tempy(:)

    ! Plotting

    ! Set NaN
    real(R8)                            :: NaN

    ! Initialize
    !===========
    ! Set NaN
    NaN = IEEE_VALUE(nan, IEEE_QUIET_NAN)
    
    ! Associate
    associate(&
        nvs         => vessel%nstructures,  &
        vs          => vessel%structures)

    ! Vessel polygon representation
    !==============================
    ! Count total number of points currently in the structures
    nvp = 0
    nexcl = 0
    do i = 1, nvs
        if (.not. any(vesseloptions%exclude == i)) then 
            nvp = vessel%structures(i)%np + nvp 
        else 
            nexcl = nexcl + 1
        end if 
    end do

    ! Allocate (account for NaNs)
    allocate(tempx(nvp+nvs-1-nexcl), tempy(nvp+nvs-1-nexcl))
    
    ! Set vertices
    cc = 0
    do i = 1, nvs 
        if (.not. any(vesseloptions%exclude == i)) then 
            ! Coordinates
            tempx(cc+1:cc+vs(i)%np) = vs(i)%x
            tempy(cc+1:cc+vs(i)%np) = vs(i)%y
            cc = cc + vs(i)%np 

            ! Add NaN
            if (i .ne. nvs) then 
                tempx(cc+1) = NaN 
                tempy(cc+1) = NaN 
                cc = cc + 1
            end if 
        end if
    end do

    ! Construct temporary vessel polygon set
    call tempps%Construct(tempx, tempy)

    ! Construct vessel polygon set 
    call ConstructVesselPolygonSet(vessel, vesseloptions, tempps)

    end associate 

end subroutine