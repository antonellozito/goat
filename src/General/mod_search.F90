!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains some basic search routines, mainly for integer
! arrays (character/string search often has fortran intrinsics such 
! as 'index', 'verify', and 'scan' - if anything else is required, it
! can be defined here). Wherever possible, we make use of (standardized)
! fortran intrinsics, such as findloc. 

module mod_search

    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_errorhandler
    
    ! The usual
    implicit none 
    private 

    public :: Findloc1D

    !==================================================================!
    !                                                                  !
    !                            TYPES                                 !
    !                                                                  !
    !==================================================================!

    ! Currently no additional types to define

    !==================================================================!
    !                                                                  !
    !                          INTERFACES                              !
    !                                                                  !
    !==================================================================!

    ! Findloc1D
    interface Findloc1D 
        module procedure Findloc1D_I8
    end interface

    contains 

    !==================================================================!
    !                                                                  !
    !                           ROUTINES                               !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                          INTEGER SEARCHING                       !
    !------------------------------------------------------------------!

    ! Find first occurrence of 1D subarray in 1D array
    function Findloc1D_I8(mainarray, subarray, back) result(ind)

        ! Description
        !============
        ! Find the first start location of the 1D subarray in the main 
        ! array. If 'back' is specified and true, then we find the last
        ! occurrence. If the subarray is not present in mainarray, then 
        ! 0 is returned. If the subarray or mainarray are empty, the 
        ! result is also 0. 

        ! Algorithm
        !==========
        ! 0) Do preliminary checks to escape trivial cases
        ! 1) Set the initial index to 1
        ! 2) Find the first element that equals the first element of 
        ! subarray in mainarray starting from the current index + 1 (if 
        ! that is not out of bounds). If this is not found -> exit with 
        ! ind = 0. Otherwise, check if the remaining entries correspond
        ! to those in subarray. If this is the case -> exit. Otherwise,
        ! update the index to the found location and repeat.

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), dimension(:), intent(in)   :: mainarray, subarray 
        logical, optional, intent(in)           :: back 
        integer(I8)                             :: ind 

        ! Auxiliary
        integer(I8)                             :: nmain, nsub, oldind 
        logical                                 :: doback, isfound 

        ! Loop
        integer(I8)                             :: i

        ! Initialize
        !===========
        ! Default value
        ind = 0

        ! Optional arguments
        if (present(back)) then 
            doback = back 
        else 
            doback = .false. 
        end if 

        ! Sanity checks
        nmain = size(mainarray)
        nsub = size(subarray)
        if (nsub > nmain) then 
            return 
        end if 
        if ((nsub == 0) .or. (nmain == 0)) then 
            return 
        end if 

        ! Find array
        !===========
        isfound = .false. 
        if (doback) then 
            ! Backward mode
            do while (.true.)
                ! Check
                if (ind == 1) then 
                    exit 
                end if 

                ! Find location of last element
                oldind = ind 
                ind = findloc(mainarray(1:oldind-1), subarray(1), 1, back=.true.) + (nmain - oldind)

                ! Check
                if (ind == (nmain - oldind)) then 
                    ind = 0
                    exit 
                end if 
                if ((ind + nsub - 1) > nmain) then 
                    ! Subarray size will exceed remaining mainarray elements, skip
                    cycle 
                end if 

                ! We found an element, check the remainder
                isfound = .true. 
                do i = 2, nsub
                    if (subarray(i) /= mainarray(ind + i - 1)) then 
                        isfound = .false. 
                        exit 
                    end if 
                end do 

                ! Check exit conditions
                if (isfound) then 
                    exit 
                end if 
            end do 
        else
            ! Forward mode
            do while (.true.)
                ! Check
                if (ind == size(mainarray)) then 
                    exit 
                end if 

                ! Find location of first element
                oldind = ind 
                ind = findloc(mainarray(oldind+1:), subarray(1), 1, back=.false.) + oldind

                ! Check
                if (ind == oldind) then 
                    ind = 0
                    exit 
                end if 
                if ((ind + nsub - 1) > nmain) then 
                    ! Subarray size will exceed remaining mainarray elements, exit
                    ind = 0
                    exit 
                end if 

                ! We found an element, check the remainder
                isfound = .true. 
                do i = 2, nsub
                    if (subarray(i) /= mainarray(ind + i - 1)) then 
                        isfound = .false. 
                        exit 
                    end if 
                end do 

                ! Check exit conditions
                if (isfound) then 
                    exit 
                end if 
            end do 
        end if 


    end function
    

end module