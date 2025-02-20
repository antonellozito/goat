!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! Simple module that provides 1D interpolation of data

module Interpolant1D

    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_constants, only: nanval_R8
    use mod_errorhandler
    use Interpolant2D_auxiliaries, only: GetBinIndex

    implicit none 
    private 
    
    public :: Interpolate1D

    !==================================================================!
    !                                                                  !
    !                            TYPES                                 !
    !                                                                  !
    !==================================================================!

    contains 

    !==================================================================!
    !                                                                  !
    !                          ROUTINES                                !
    !                                                                  !
    !==================================================================!

    ! Simple 1D interpolation routine
    subroutine Interpolate1D(xq, yq, xv, yv)

        ! Description
        !============
        ! Simple 1D interpolation of real data yv that is given on 
        ! datapoints xv. It is assumed that xv is monotonically 
        ! increasing (this is not explicitly checked). If the query
        ! point coordinates are out of bounds, NaNs are returned. 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)                :: xv(:), yv(:), xq(:)
        real(R8), allocatable, intent(out)  :: yq(:)

        ! Auxiliary
        real(R8)                            :: dx, dy 
        integer(I8)                         :: ind 
        
        ! Loop
        integer(I8)                         :: i 

        ! Initialize
        !===========
        ! Checks
        if (size(xv) /= size(yv)) then 
            call gdErrorHandler('Interpolate1D: incompatible dimensions ' // & 
                'of xv and yv')
        end if 
        if (allocated(yq)) then 
            if (size(yq) /= size(xq)) then 
                deallocate(yq)
                allocate(yq(size(xq)))
            end if 
        else
            allocate(yq(size(xq)))
        end if 


        ! Interpolate
        !============
        do i = 1, size(xq)
            ! Get segment index
            call GetBinIndex(xq(i), xv, size(xv), ind)

            ! Check 
            if (ind == 0) then 
                yq(i) = nanval_R8()
            else
                dx = xv(ind+1) - xv(ind)
                dy = yv(ind+1) - yv(ind)
                if (dy == 0_R8) then 
                    yq(i) = yv(ind)
                elseif (dx == 0_r8) then
                    ! This actually means there's a discontinuity, so we
                    ! just take the half of dy... (discontinuity only if 
                    ! dy /= 0)
                    yq(i) =  0.5*dy + yv(ind)
                else
                    ! Standard interpolation
                    yq(i) = dy/dx*(xq(i) - xv(ind)) + yv(ind)
                end if
            end if 
        end do 

    end subroutine 

end module 
