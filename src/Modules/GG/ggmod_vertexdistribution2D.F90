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
    use Interpolant1D
    use mod_lagrangefunctions
    use mod_constants, only: nanval_R8
    use DistributionFunction, only: DistributionFunctionUDT
    implicit none
    private 

    public :: VertexDistributor2DUDT, ConstructUniformVertexDistributor, &
        FieldDistributor1DUDT, ConstructUniformFieldDistributor

    ! Module parameters
    real(R8), parameter, private        :: tprelfieldtol = 1e-10 ! relative field tolerance under which extrema are removed
    real(R8), parameter, private        :: disttol = 1e-12 ! distance tolerance

    !==================================================================!
    !                                                                  !
    !                            TYPES                                 !
    !                                                                  !
    !==================================================================!

    ! Vertex distribution
    !====================
    ! General vertex distributor
    type, abstract :: VertexDistributor2DUDT

        ! Description
        !============
        ! General abstract type that mainly stipulates which routines
        ! a vertex distributor should have. This is currently only the
        ! 'DistributeOverCurve' subroutine, which yields vertex 
        ! coordinates that lie on a given curve.

    contains 

        ! Distribution over given coordinates
        procedure(DistributeVerticesOverCurveINT), deferred :: DistributeOverCurve

        ! Distribution over field
        procedure(DistributeVerticesOverFieldINT), deferred :: DistributeOverField

    end type 

    


    ! Uniform vertex distributor
    type, extends(VertexDistributor2DUDT) :: UniformVertexDistributor2DUDT

        ! Description
        !============
        ! Distributor that distributes vertices uniformly with a spacing
        ! equal to 'd' (this is an input parameter of the distribution
        ! function) in the case of a curve based distribution or with
        ! field values separated by a distance 'fd' in the case of a
        ! field based distribution. These distances can only be 
        ! approximated, since the start and end points are always taken
        ! to be on the given 

        real(R8)                        :: d, fd

    contains 

        ! Distribution over given coordinates
        procedure :: DistributeOverCurve    => DistributeVerticesUniformOverCurve

        ! Distribution over given field
        procedure :: DistributeOverField    => DistributeVerticesUniformOverField

    end type

    ! Field distribution
    !===================
    ! General field distributor
    type, abstract :: FieldDistributor1DUDT

        ! Description
        !============
        ! General abstract type that mainly stipulates which routines
        ! a field distributor should have. This is currently only the
        ! 'DistributeOverCurve' subroutine, which yields vertex 
        ! coordinates that lie on a given curve. The evaluation routine
        ! expects a distribution from the class 'DistributionFunctionUDT'
        ! that can be evaluated at x, y coordinates

    contains 

        ! Distribution of field values over given coordinates
        procedure(DistributeFieldOverCurveINT), deferred :: DistributeOverCurve

    end type 

    ! Uniform field distributor
    type, extends(FieldDistributor1DUDT) :: UniformFieldDistributor1DUDT

        ! Description
        !============
        ! Distributor that distributes vertices uniformly with a spacing
        ! equal to 'd' (this is an input parameter of the distribution
        ! function) in the case of a curve based distribution or with
        ! field values separated by a distance 'fd' in the case of a
        ! field based distribution. These distances can only be 
        ! approximated, since the start and end points are always taken
        ! to be on the given 

        real(R8)                        :: d, fd

    contains 

        ! Distribution over given coordinates
        procedure :: DistributeOverCurve    => DistributeFieldUniformOverCurve

    end type


    !==================================================================!
    !                                                                  !
    !                         INTERFACES                               !
    !                                                                  !
    !==================================================================!

    abstract interface 

        ! Vertex distribution over curve
        subroutine DistributeVerticesOverCurveINT(vd, xc, yc, xv, yv, nv)

            import :: I8, R8, VertexDistributor2DUDT
            class(VertexDistributor2DUDT)          :: vd 
            real(R8), intent(in)                :: xc(:), yc(:)
            real(R8), allocatable, intent(out)  :: xv(:), yv(:)
            integer(I8), intent(out)            :: nv 

        end subroutine 

        ! Vertex distribution over field
        subroutine DistributeVerticesOverFieldINT(vd, xc, yc, field, &
            xv, yv, nv)

            import :: I8, R8, VertexDistributor2DUDT, DistributionFunctionUDT
            class(VertexDistributor2DUDT)       :: vd 
            class(DistributionFunctionUDT), intent(in)  :: field
            real(R8), intent(in)                :: xc(:), yc(:)
            real(R8), allocatable, intent(out)  :: xv(:), yv(:)
            integer(I8), intent(out)            :: nv 

        end subroutine 

        ! Field distribution over curve
        subroutine DistributeFieldOverCurveINT(vd, xc, yc, field, nv, fv)

            import :: I8, R8, FieldDistributor1DUDT, &
                DistributionFunctionUDT
            class(FieldDistributor1DUDT)            :: vd
            real(R8), intent(in)                    :: xc(:), yc(:)
            real(R8), allocatable, intent(out)      :: fv(:)
            integer(I8), intent(out)                :: nv 
            class(DistributionFunctionUDT), intent(in)  :: field    

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
    
    ! Uniform vertex distributor
    function ConstructUniformVertexDistributor(d, fd) result(vd)

        ! Description
        !============
        ! Construct a uniform vertex distributor. Options are directly
        ! passed through the constructor. these are:
        ! - d: desired distance between vertices using the curve based
        ! distribution function
        ! - fd: desired distance measured in the field width that is 
        ! given using the field based distribution function

        ! Declare variables
        !==================
        ! Arguments
        class(VertexDistributor2DUDT), allocatable :: vd 
        real(R8), intent(in)                    :: d, fd

        ! Initialize
        !===========
        ! Allocate
        allocate(UniformVertexDistributor2DUDT::vd)

        ! Set the fields
        select type (vd)

        type is (UniformVertexDistributor2DUDT)

            ! Set parameters
            vd%d = d ! distance 
            vd%fd = fd ! field distance

        end select 

    end function 

    ! Uniform field distributor
    function ConstructUniformFieldDistributor(d, fd) result(vd)

        ! Description
        !============
        ! Construct a uniform vertex distributor. Options are directly
        ! passed through the constructor. these are:
        ! - d: desired distance between vertices using the curve based
        ! distribution function
        ! - fd: desired distance measured in the field width that is 
        ! given using the field based distribution function

        ! Declare variables
        !==================
        ! Arguments
        class(FieldDistributor1DUDT), allocatable   :: vd 
        real(R8), intent(in)                        :: d, fd

        ! Initialize
        !===========
        ! Allocate
        allocate(UniformFieldDistributor1DUDT::vd)

        ! Set the fields
        select type (vd)

        type is (UniformFieldDistributor1DUDT)

            ! Set parameters
            vd%d = d ! distance 
            vd%fd = fd ! field distance

        end select 

    end function 

    !------------------------------------------------------------------!
    !                    UNIFORM VERTEX DISTRIBUTOR                    !
    !------------------------------------------------------------------!

    ! Curve based distributor
    subroutine DistributeVerticesUniformOverCurve(vd, xc, yc, xv, yv, nv)

        ! Description
        !============
        ! Distribute the vertex coordinates xv, yv uniformly over the
        ! curve given by the coordinates xc, yc. 

        ! Declare variables
        !==================
        ! Arguments
        class(UniformVertexDistributor2DUDT):: vd 
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

    ! Field based distributor
    subroutine DistributeVerticesUniformOverField(vd, xc, yc, field, &
        xv, yv, nv)

        ! Description
        !============
        ! Distribute the vertices over the given curve xc, yc based on 
        ! uniform distance in terms of the field values fc that are defined
        ! on the curve coordinates. fc is made monotonic by building the
        ! interpolation coordinate as the sum of the absolute value of
        ! the difference in fc (which is therefore always monotonic). 

        ! Declare variables
        !==================
        ! Arguments
        class(UniformVertexDistributor2DUDT)        :: vd 
        class(DistributionFunctionUDT), intent(in)  :: field 
        real(R8), intent(in)                :: xc(:), yc(:)
        real(R8), allocatable, intent(out)  :: xv(:), yv(:)
        integer(I8), intent(out)            :: nv 

        ! Auxiliary
        real(R8)                            :: l 
        real(R8), allocatable, dimension(:) :: dx, dy, dl, distr, fc

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
        allocate(fc(size(xc)))
        call field%Evaluate(xc, yc, fc) 
        dl = abs(fc(2:) - fc(1:size(fc)-1))
        l = sum(dl)

        ! Compute number of vertices
        nv = ceiling(l/vd%fd)+1

        ! Compute distribution
        distr = [(k, k = 0, nv)]*(l/nv)
        distr(size(distr)) = l ! just to make sure

        ! Distribute
        call DistributeVerticesLine(xc, yc, dl, distr, xv, yv)

    end subroutine
    
    
    !------------------------------------------------------------------!
    !                     UNIFORM FIELD DISTRIBUTOR                    !
    !------------------------------------------------------------------!

    ! Curve based distributor
    subroutine DistributeFieldUniformOverCurve(vd, xc, yc, field, nv, fv)

        ! Declare variables
        !==================
        ! Arguments
        class(UniformFieldDistributor1DUDT)         :: vd 
        real(R8), intent(in)                        :: xc(:), yc(:)
        real(R8), allocatable, intent(out)          :: fv(:)
        integer(I8), intent(out)                    :: nv 
        class(DistributionFunctionUDT), intent(in)  :: field
        
        ! Auxiliary
        real(R8)                        :: fval(1:2), df 

        ! Loop
        integer(I8)                     :: k 

        ! Distribute
        !===========
        ! Very simple
        call field%Evaluate([xc(1), xc(ubound(xc))], [yc(1), yc(ubound(yc))], &
            fval)
        df = fval(2) - fval(1)
        nv = ceiling(abs(df)/vd%fd) + 1 
        fv = [(k, k = 0, nv-1)]*df/(nv-1)

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

        ! Loop
        integer(I8)             :: i 

        ! Initialize
        !===========
        ! Construct the length axis
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