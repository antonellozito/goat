!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module provides vertex distributors that distribute vertices 
! along a certain field. 

module ggmod_vertexdistribution2D

    ! Load modules
    use mod_precision
    use mod_errorhandler
    use mod_dynamicarrays
    use mod_contour2D
    use mod_polygon
    use mod_sort
    use Interpolant1D
    use mod_lagrangefunctions
    use mod_constants, only: nanval_R8
    use mod_linearsolverinterface, only: SolveDenseLinearSystemDI
    use goatmod_types, only : magneticFieldUDT, VesselUDT, GridUDT
    use goatmod_userinput, only : TopomeshOptionsUDT
    implicit none
    private 

    ! Module parameters
    real(R8), parameter, private        :: tprelfieldtol = 1e-10 ! relative field tolerance under which extrema are removed
    real(R8), parameter, private        :: disttol = 1e-12 ! distance tolerance

    !==================================================================!
    !                                                                  !
    !                            TYPES                                 !
    !                                                                  !
    !==================================================================!

    ! General vertex distributor
    type, abstract :: VertexDistributor2D 

        ! Description
        !============
        ! General abstract type that mainly stipulates which routines
        ! a vertex distributor should have. This is currently only the
        ! 'DistributeOverCurve' subroutine, which yields vertex 
        ! coordinates that lie on a given curve.

    contains 

        ! Distribution over given coordinates
        procedure(DistributeOverCurveINT), deferred :: DistributeOverCurve

    end type 


    ! Uniform vertex distributor
    type, extends(VertexDistributor2D) :: UniformVertexDistributor2D

        ! Description
        !============
        ! Distributor that distributes vertices uniformly with a spacing
        ! equal to 'd' (this is an input parameter of the distribution
        ! function)

        real(R8)                        :: d 

    contains 

        ! Distribution over given coordinates
        procedure :: DistributeOverCurve    => DistributeUniformOverCurve

    end type

    !==================================================================!
    !                                                                  !
    !                         INTERFACES                               !
    !                                                                  !
    !==================================================================!

    abstract interface 

        ! Distribution over curve
        subroutine DistributeOverCurveINT(vd, xc, yc, xv, yv, nv)

            import :: I8, R8, VertexDistributor2D
            class(VertexDistributor2D)          :: vd 
            real(R8), intent(in)                :: xc(:), yc(:)
            real(R8), allocatable, intent(out)  :: xv(:), yv(:)
            integer(I8), intent(out)            :: nv 

        end subroutine 

    end interface 


    contains 

    !==================================================================!
    !                                                                  !
    !                           ROUTINES                               !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                           CONSTRUCTORS                           !
    !------------------------------------------------------------------!
    
    ! Uniform distributor
    function ConstructUniformVertexDistributor(d) result(vd)

        ! Description
        !============
        ! Construct a uniform vertex distributor. Options are directly
        ! passed through the constructor

        ! Declare variables
        !==================
        ! Arguments
        class(VertexDistributor2D), allocatable :: vd 
        real(R8), intent(in)                    :: d 

        ! Initialize
        !===========
        ! Allocate
        allocate(UniformVertexDistributor2D::vd)

        ! Set the fields
        select type (vd)

        type is (UniformVertexDistributor2D)

            ! Set parameters
            vd%d = d ! distance 

        end select 

    end function 

    !------------------------------------------------------------------!
    !                       UNIFORM DISTRIBUTOR                        !
    !------------------------------------------------------------------!

    ! Distributor
    subroutine DistributeUniformOverCurve(vd, xc, yc, xv, yv, nv)

        ! Description
        !============
        ! Distribute the vertex coordinates xv, yv uniformly over the
        ! curve given by the coordinates xc, yc. 

        ! Declare variables
        !==================
        ! Arguments
        class(UniformVertexDistributor2D)   :: vd 
        real(R8), intent(in)                :: xc(:), yc(:)
        real(R8), allocatable, intent(out)  :: xv(:), yv(:)
        integer(I8), intent(out)            :: nv 

        ! Auxiliary
        real(R8)                            :: l 
        real(R8), allocatable, dimension(:) :: dx, dy, dl, distr

        ! Loop
        integer(I8)                         :: k 


        ! Initialize
        !===========
        ! Checks
        if (size(xc) /= size(yc)) then 
            call gdErrorHandler('DistributeUniformOverCurve: curve ' // & 
                'coordinates have incompatible dimensions')
        end if 
        if (allocated(xv)) then 
            deallocate(xv)
        end if 
        if (allocated(yv)) then 
            deallocate(yv)
        end if 

        ! Precompute
        !===========
        ! Compute curve quantities
        dx = xc(2:) - xc(1:size(xc)-1)
        dy = yc(2:) - yc(1:size(yc)-1)  
        dl = sqrt(dx**2 + dy**2)
        l = sum(dl)

        ! Compute number of vertices
        nv = ceiling(l/vd%d)+1

        ! Compute distribution
        distr = [(k, k = 0, nv)]*(l/nv)
        distr(size(distr)) = l ! just to make sure

        ! Distribute
        call DistributeVerticesLine(xc, yc, dl, distr, xv, yv)

    end subroutine 


    !------------------------------------------------------------------!
    !                    THE ONES THAT DO THE WORK                     !
    !------------------------------------------------------------------!

    ! Line distributor
    subroutine DistributeVerticesLine(xl, yl, dl, distr, xv, yv)

        ! Description
        !============
        ! Return a set of vertices xv, yv, distributed so that the length between
        ! these vertices is approximately equal to the desired length distribution
        ! specified in distr. The reason this is approximate, is that we need to
        ! interpolate on the line edges, and that if there are (sharp) angles
        ! between the line vertices, the actual obtained length can differ
        ! (substantially) from the desired one. If you want to capture these sharp
        ! features, you need to add points/refine the distribution...

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)    :: xl(:), yl(:), dl(:), distr(:)
        real(R8), allocatable, intent(out)  :: xv(:), yv(:)

        ! Auxiliary
        real(R8)                :: dlc(1:size(xl))
        real(R8), allocatable   :: dx(:), dy(:)

        ! Loop
        integer(I8)             :: i 

        ! Initialize
        !===========
        ! Construct the length axis
        dx = xl(2:) - xl(1:size(xl)-1)
        dy = yl(2:) - yl(1:size(yl)-1)
        dlc(1) = 0
        do i = 1, size(dl) ! cumulative sum
            dlc(i+1) = dl(i) + dlc(i)
        end do 

        ! Interpolate
        !============
        call Interpolate1D(distr, xv, dlc, xl)
        call Interpolate1D(distr, yv, dlc, yl)

    end subroutine

end module 