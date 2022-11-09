!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! The purpose of this module is to interface with other programs (e.g 
! the grid generator) and to translate their data types and variables to
! those used in the grid deformation modules. 

module gdmod_interfaces

    ! Initialize
    !============
    ! Load modules
    use gdmod_types
    use gdmod_userinput 

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
    !                              Grid data                           !
    !------------------------------------------------------------------!

    ! Boundaries
    !===========
    subroutine InterfaceBoundaryMapping(from, labelsfrom, labelsto, &
        bndmapping)

        ! Description
        !============
        ! This routine simply returns the boundary ID bndmapping to 
        ! communicate between e.g. a grid generator and its output flags 
        ! and the grid deformation routine. 
    
        ! Notes
        !======
        ! Note 1: right now, the 'from' cases are assumed to be based on the
        ! grid generator output, e.g. carre. 
    
        ! Boundary definitions: see gdmod_types.F90, BndUDT implementation
    
        ! Initialize
        !===========
        ! Declare modules
        use gdmod_types 
    
        ! The usual
        implicit none 
    
        ! Declare variables
        !==================
        ! Arguments
        character(*)                :: from ! input
        integer(I8), allocatable    :: labelsfrom(:), labelsto(:), bndmapping(:,:)
        integer(I8)                 :: nlf, nlt
    
        ! Set the mapping
        !================
        select case(from)
    
        case ('carre')
    
            ! Mapping description
            !====================
            ! CARRE         GRID DEFORMATION            PHYSICAL MEANING
            !
            ! -13           1                           target plate (inner)
            ! -34           2                           target plate (outer)
            ! -23           3                           private flux
            ! -24           3                           private flux
            ! -21           4                           core boundary
            ! -42           5                           outermost flux surf.
            ! -43           5                           outermost flux surf.
            ! -44           5                           outermost flux surf.
            !
            ! Note that information may be lost in this bndmapping, but this 
            ! is information that is not used by the grid deformation 
            ! module. 
            
            ! Define mapping
            !===============
            ! Initialize and allocate
            nlf = 8
            nlt = 5
            allocate(labelsfrom(nlf))
            allocate(labelsto(nlt))
            allocate(bndmapping(nlf,2))

            ! Set labels
            labelsfrom      = (/-13, -34, -23, -24, -21, -42, -43, -44/)
            labelsto        = (/1,  2, 3, 4, 5/)
    
            ! Set bndmapping
            bndmapping(:,1) = labelsfrom
            bndmapping(:,2) = (/1, 2, 3,   3,   4,   5,   5,   5/)
            
    
        end select
    
    end subroutine

    !------------------------------------------------------------------!
    !                              Vessel data                         !
    !------------------------------------------------------------------!
    ! Given the different vessel input formats, several reading routines
    ! are provided here. Also, the vessel data extraction routines are
    ! given here, since these can be very specific for each use case. 

    ! Reading vessel data from structure.dat file
    subroutine read_structure(filespecifier, vessel, vesseloptions)

        ! Description
        !============
        ! Read the vessel data from a structure.dat file. This file 
        ! should be strictly formatted (though not checked) as follows:
        ! - The first line should contain the number of structures 
        ! - The next two lines are skipped
        ! - Starting from the fourth line, the structures follow. Each 
        ! structure should first start with a header (which is ignored)
        ! and subsequently with the number of points in the 
        ! structure, where the sign of this number indicates whether the 
        ! polygon should be closed (negative) or open (positive).

        ! Notes:
        ! The vessel should have an allocatable array of substructures
        ! of the type 'VesselStructure'. This array is allocated while
        ! reading in the number of vessel structures.

        ! Declare variables
        !==================
        ! Arguments
        type(VesselUDT)                         :: vessel
        type(VesselOptionsUDT), intent(in)      :: vesseloptions
        integer(I8), intent(in)                 :: filespecifier

        ! Loop variables
        integer(I8)                             :: i, j

        ! Auxiliary variables
        integer(I4)                             :: nstruct, npoints
        character(C32)                          :: dummy

        ! Initialize
        !===========
        ! Print from where we're reading
        print *, 'reading vessel from file: ', vesseloptions%dir

        ! Open the file
        open(unit = filespecifier, file = vesseloptions%dir)

        ! Read
        !=====
        ! Read the amount of structures
        read(filespecifier, *) nstruct
        print *, 'there are ', nstruct, ' structures present'
        print *, 'vessel structure ID | number of points | is closed'

        ! Allocate
        vessel%nstructures = nstruct
        allocate(vessel%structures(nstruct))

        ! Skip the next line
        read(filespecifier, *)

        ! Read in structures
        do i = 1, nstruct
            ! Read the header (structure <structureID>)
            read(filespecifier, *) dummy, vessel%structures(i)%ID

            ! Read the number of points
            read(filespecifier, *) npoints

            ! Allocate
            vessel%structures(i)%np = abs(npoints)
            vessel%structures(i)%isclosed = (npoints .gt. 0)
            call AllocateVesselStructure(vessel%structures(i))

            ! Print
            print *, vessel%structures(i)%ID,  vessel%structures(i)%np, vessel%structures(i)%isclosed

            ! Read coordinates
            do j = 1, vessel%structures(i)%np
                ! First x, then y coordinate
                read(filespecifier, *) vessel%structures(i)%x(j), vessel%structures(i)%y(j)
            end do
        end do

        ! Close the file
        !===============
        close(filespecifier)
        

    end subroutine
    


end module