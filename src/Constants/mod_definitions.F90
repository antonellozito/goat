!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! General definitions to be used in other modules. These definitions can
! be for example the integer that is linked to a specific type of 
! vessel boundary (e.g. target plate). These definitions should replace
! hard-coded constants etc. 

module mod_definitions

    use mod_precision

    implicit none
    public 

    ! Vessel parts definitions
    integer(I8), parameter  :: targetID = 1, coreID = 2, outerboundaryID = 3, &
        vesselID = 4, interiorID = 5

   
end module