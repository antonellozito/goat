!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! The purpose of this model is to interface with other programs (e.g 
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
    


end module