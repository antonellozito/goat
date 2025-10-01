!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains math routines for the adapation modules

module gamod_math

    ! Initialize
    !===========
    ! Load modules
    use mod_precision


    ! The usual
    implicit none
    save
    public

    !==================================================================!
    !                                                                  !
    !                          INTERFACES                              !
    !                                                                  !
    !==================================================================!

    ! General norm
    interface Norm
        module procedure Norm0D, Norm1D
    end interface

    contains

    function Intersects(p1, q1, p2, q2) result(res)
        real(R8) :: p1(2), q1(2), p2(2), q2(2)
        integer(I8) :: res, o1, o2, o3, o4


        o1 = Orient(p1, q1, p2)
        o2 = Orient(p1, q1, q2)
        o3 = Orient(p2, q2, p1)
        o4 = Orient(p2, q2, q1)

        res = 1

        if (o1 /= o2 .and. o3 /= o4) then
            return
        end if

        !c1 = 
        if (o1 == 0 .and. OnSegment(p1, p2, q1) == 1) then 
            return
        end if
    
        if (o2 == 0 .and. OnSegment(p1, q2, q1) == 1) then
            return
        end if
    
        if (o3 == 0 .and. OnSegment(p2, p1, q2) == 1) then
            return 
        end if
    
        if (o4 == 0 .and. OnSegment(p2, q1, q2) == 1) then
            return 
        end if

        res = 0

    end function
    
    function Orient(p, q, r)  result(res)
        real(R8) :: p(2), q(2), r(2), val
        integer(I8) :: res

        val = (q(2) - p(2)) * (r(1) - q(1)) - (q(1) - p(1)) * (r(2) - q(2))

        if (val > 0) then

            res = 1
            return

        else

            if (val < 0) then
                res = 2
                return
            else
                res = 0
                return
            end if

        end if

    end function

    function OnSegment(p, q, r) result(res)
        real(R8) :: p(2), q(2), r(2)
        integer(I8) :: res 

        res = 0
        if ( (q(1) <= max(p(1), r(1))) .and. (q(1) >= min(p(1), r(1))) .and. &
             (q(2) <= max(p(2), r(2))) .and.  (q(2) >= min(p(2), r(2))) ) then
            res = 1
        end if
        
        return               
    end function  
    
    function TriangleArea(x0, y0, x1, y1, x2, y2) result(res)
        real(R8) :: x0, y0, x1, y1, x2, y2, res

        res = 0.5_R8 * abs( (x1-x0)*(y2-y0) - (y1-y0)*(x2-x0) )

    end function 
    
    function Norm0D(x0, y0) result(res)
        real(R8) :: x0, y0, res
        res = sqrt(x0**2 + y0**2)
    end function

    function Norm1D(x0, y0) result(res)
        real(R8) :: x0(:), y0(:)
        real(R8), allocatable :: res(:)
        res = sqrt(x0**2 + y0**2)
    end function    
    
    
end module 