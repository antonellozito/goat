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
    use StructuredInterpolant2D
    use mod_lagrangefunctions
    use mod_constants, only: nanval_R8
    use DistributionFunction, only: DistributionFunctionUDT
    implicit none
    private 

    public :: VertexDistributor2DUDT, ConstructUniformVertexDistributor, &
        FieldDistributor1DUDT, ConstructUniformFieldDistributor, &
        ConstructDensityBasedVertexDistributor

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

    ! Uniform vertex distributor
    type, extends(VertexDistributor2DUDT) :: DensityBasedVertexDistributor2DUDT

        ! Description
        !============
        ! Construct a distributor that works density based. Here, it is assumed
        ! that a 'density' distribution (e.g. in #/m) is given in 'options'. The
        ! distribution is then integrated along the polygon, which determines the
        ! number of vertices to be distributed thereon. The location is then
        ! obtained by finding the location at which the accumulated mass equals
        ! 1/N times the total mass of the system. Hence, where the density is
        ! higher, the vertices will be concentrated.
        
        ! For integration of the density, we work per segment and use a polynomial
        ! representation of order 'vd.order'. In higher than linear order cases, we
        ! introduce additional nodes in an equidistant fashion. By parameterization
        ! per segment, we obtain a single coefficient matrix to be inverted for
        ! each segment. This can all be done in pre-processing and should therefore
        ! be quite efficient (only one matrix inversion for all to be interpolated
        ! quantities). As basis function, we us Lagrange polynomials. These have
        ! the nice feature that at any given point, all basis functions vanish
        ! except for one, which is equal to one, such that we can immediately see
        ! that the coefficients should be equal to the density at those points.
        ! Hence, the inverse of this system is trivial. In other words:
        !
        !   rho(s) = sum(rho_i l_i)
        !
        ! where l_i is the i-th Lagrange polynomial. This polynomial is given by:
        !
        !   l_i(s) = prod( (x - x_m)(x_i - x_m) ) for m = 0..k, m != i
        !
        ! If desired, one can also write the Lagrange polynomial in the form of
        ! l_i = prod( 1/(x_i - x_m) ) sum_i a_i x^i, which is convenient for
        ! integration/differentiation. Here, the a_i are computed recursively
        ! (perhaps there's a closed and elegant formulation but haven't found that
        ! one yet). If we start with the 'constant' part (m = 0) and move further
        ! to m = M, then, we can find each set of coefficients of m+1 as function
        ! of the coefficients of m (note: we assume here that i is skipped, but do
        ! not reflect this in notation):
        !
        !   a_j^m+1 = a_j-1^m (-x_m+1) + (a_j^m)
        !
        ! For example:
        !
        !   m   |   a_0^m   a_1^m-1     ...
        !   0   |   1       0     0     ...
        !   1   |   1       -x_1    0
        !   2   |   1       -(x_1 + x_2)    x_1 x_2
        !   3   |   1       -(x_1 + x_2 + x_3)  x_1 x_2 + x_1 x_3 + x_2 x_3     -x_1 x_2 x_3
        !
        ! For the definition of l_i above, this would result for example for l_1
        ! and M = 3 in:
        !
        !   m   |   a_0^m   a_1^m-1     ...
        !   0   |   1       0     0     ...
        !   ----------------------------------
        !   2   |   1       -(x_2)
        !   3   |   1       -(x_2 + x_3)  x_2 x_3
        !
        ! Finally, we note that here, the coefficients will immediately include the
        ! product prefactor (i.e. we cosider b_i = a_i prod( 1/(x_i - x_m) ) )

        integer(I8)                                     :: order 
        class(DistributionFunctionUDT), allocatable     :: densityfunction
        real(R8), allocatable, dimension(:)             :: xi
        real(R8), allocatable, dimension(:, :)          :: lagcoef, intlagcoef

        real(R8)                        :: d, fd

    contains 

        ! Distribution over given coordinates
        procedure :: DistributeOverCurve    => DistributeVerticesDensityBasedOverCurve

        ! Distribution over given field
        procedure :: DistributeOverField    => DistributeVerticesDensityBasedOverField2

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
        subroutine DistributeVerticesOverCurveINT(vd, xc, yc, nv, xv, yv, &
            ldistr)

            import :: I8, R8, VertexDistributor2DUDT
            class(VertexDistributor2DUDT)          :: vd 
            real(R8), intent(in)                :: xc(:), yc(:)
            real(R8), allocatable, intent(out), optional  :: xv(:), yv(:)
            integer(I8), intent(out)            :: nv 
            real(R8), allocatable, intent(out), optional :: ldistr(:)

        end subroutine 

        ! Vertex distribution over field
        subroutine DistributeVerticesOverFieldINT(vd, xc, yc, field, nv, &
            xv, yv, ldistr)

            import :: I8, R8, VertexDistributor2DUDT, DistributionFunctionUDT
            class(VertexDistributor2DUDT)       :: vd 
            class(DistributionFunctionUDT), intent(in)  :: field
            real(R8), intent(in)                :: xc(:), yc(:)
            real(R8), allocatable, intent(out), optional  :: xv(:), yv(:), ldistr(:)
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

    ! Density-based vertex distributor
    function ConstructDensityBasedVertexDistributor(distribution, order) result(vd)

        ! Description
        !============
        ! Construct a density-based vertex distribution function based
        ! on the density distribution which is given in 'distribution'. 
        ! For proper performance, the distribution function should be
        ! (strictly) positive everywhere. This cannot be checked 
        ! explicitly beforehand, but during evaluation this may be 
        ! checked. 
 
        ! Declare variables
        !==================
        ! Arguments
        class(VertexDistributor2DUDT), allocatable  :: vd 
        class(DistributionFunctionUDT), intent(in)  :: distribution
        integer(I8), intent(in)                     :: order

        
        ! Loop
        integer(I8)                                 :: k

        ! Initialize
        !===========
        ! Allocate
        allocate(DensityBasedVertexDistributor2DUDT::vd)

        ! Initialize
        select type(vd)

        type is (DensityBasedVertexDistributor2DUDT)

            ! Set some fields
            vd%order = order
            vd%densityfunction = distribution

            ! Construct the lagrangian basis functions
            vd%xi = real([(k, k = 0, vd%order)], kind=R8)/(real(vd%order, kind=R8))
            call ConstructLagrangianBasisFunctions(vd%order, vd%xi, &
                vd%lagcoef, vd%intlagcoef)

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
    subroutine DistributeVerticesUniformOverCurve(vd, xc, yc, nv, xv, yv, &
        ldistr)

        ! Description
        !============
        ! Distribute the vertex coordinates xv, yv uniformly over the
        ! curve given by the coordinates xc, yc. 

        ! Declare variables
        !==================
        ! Arguments
        class(UniformVertexDistributor2DUDT):: vd 
        real(R8), intent(in)                :: xc(:), yc(:)
        real(R8), allocatable, intent(out), optional  :: xv(:), yv(:)
        integer(I8), intent(out)            :: nv 
        real(R8), allocatable, intent(out), optional :: ldistr(:)

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
        if (present(xv)) then 
            if (allocated(xv)) then 
                deallocate(xv)
            end if 
        end if 
        if (present(yv)) then 
            if (allocated(yv)) then 
                deallocate(yv)
            end if 
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
        distr = [(k, k = 0, nv-1)]*(l/(nv-1))
        distr(size(distr)) = l ! just to make sure

        ! Distribute
        if (present(xv) .and. present(yv)) then 
            call DistributeVerticesLine(xc, yc, dl, distr, xv, yv)
        end if 

        ! Return length distribution
        if (present(ldistr)) then 
            ldistr = distr 
        end if 

    end subroutine

    ! Field based distributor
    subroutine DistributeVerticesUniformOverField(vd, xc, yc, field, nv, &
        xv, yv, ldistr)

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
        integer(I8), intent(out)            :: nv 
        real(R8), allocatable, intent(out), optional    :: ldistr(:), xv(:), yv(:)

        ! Auxiliary
        real(R8)                            :: l 
        real(R8), allocatable, dimension(:) :: dx, dy, dl, distr, fc, &
            dllc, dlc, dll

        ! Loop
        integer(I8)                         :: k 

        ! Initialize
        !===========
        ! Checks
        if (size(xc) /= size(yc)) then 
            call gdErrorHandler('DistributeUniformOverCurve: curve ' // & 
                'coordinates have incompatible dimensions')
        end if 
        if (present(xv)) then 
            if (allocated(xv)) then 
                deallocate(xv)
            end if 
        end if 
        if (present(yv)) then 
            if (allocated(yv)) then 
                deallocate(yv)
            end if 
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
        distr = [(k, k = 0, nv-1)]*(l/(nv-1))
        distr(size(distr)) = l ! just to make sure
        if (present(ldistr)) then 
            ! This has to be the length distribution!
            dll = sqrt(dx**2 + dy**2)
            dllc = spread(0, 1, size(xc))
            dlc = spread(0, 1, size(xc))
            do k = 1, size(xc)-1
                dllc(k+1) = dllc(k) + dll(k)
            end do 
            do k = 1, size(xc)-1
                dlc(k+1) = dlc(k) + dl(k)
            end do 
            ldistr = distr 
            call Interpolate1D(distr, ldistr, dlc, dllc)
        end if 

        ! Distribute
        if (present(xv) .and. (present(yv))) then 
            call DistributeVerticesLine(xc, yc, dl, distr, xv, yv)
        end if 

    end subroutine

    !------------------------------------------------------------------!
    !                 DENSITY BASED VERTEX DISTRIBUTOR                 !
    !------------------------------------------------------------------!

    ! Curve based distributor
    subroutine DistributeVerticesDensityBasedOverCurve(vd, xc, yc, nv, xv, yv, &
        ldistr)

        ! Description
        !============
        ! Distribute the vertex coordinates xv, yv non-uniformly over the
        ! curve given by the coordinates xc, yc. 

        ! Declare variables
        !==================
        ! Arguments
        class(DensityBasedVertexDistributor2DUDT)   :: vd 
        real(R8), intent(in)                :: xc(:), yc(:)
        real(R8), allocatable, intent(out), optional  :: xv(:), yv(:)
        integer(I8), intent(out)            :: nv 
        real(R8), allocatable, intent(out), optional :: ldistr(:)

        ! Auxiliary
        real(R8)                            :: l, Mtot
        real(R8), allocatable, dimension(:) :: dx, dy, dl, distr, temp, &
            Mi, Mdistr, dll, dllc, dlc
        real(R8), allocatable, dimension(:, :)  :: xi, yi, rhoi 
        integer(I8), allocatable, dimension(:)  :: pxi

        ! Loop
        integer(I8)                         :: i, k 


        ! Initialize
        !===========
        ! Checks
        if (size(xc) /= size(yc)) then 
            call gdErrorHandler('DistributeDensityBasedOverCurve: curve ' // & 
                'coordinates have incompatible dimensions')
        end if 
        if (present(xv)) then 
            if (allocated(xv)) then 
                deallocate(xv)
            end if 
        end if 
        if (present(yv)) then 
            if (allocated(yv)) then 
                deallocate(yv)
            end if 
        end if 

        ! Associate
        associate(&
            xb          => vd%xi, &
            rho         => vd%densityfunction, &
            coef        => vd%lagcoef,  &
            intcoef     => vd%intlagcoef) 

        ! Precompute
        !===========
        ! Compute curve quantities
        dx = xc(2:) - xc(1:size(xc)-1)
        dy = yc(2:) - yc(1:size(yc)-1)  
        dl = sqrt(dx**2 + dy**2)
        l = sum(dl)

        ! Compute lagrange interpolation 
        xi = spread(xb, 1, size(dx))
        yi = spread(xb, 1, size(dy))
        do i = 1, size(dx)
            xi(i, :) = xi(i, :)*dx(i) + xc(i)
            yi(i, :) = yi(i, :)*dy(i) + yc(i)
        end do 
        allocate(temp(size(xi)))
        call rho%Evaluate(reshape(xi, [size(xi)]), reshape(yi, [size(yi)]), temp)
        rhoi = reshape(temp, [size(xi, 1), size(xi, 2)])
        pxi = [(k, k = 0, vd%order)]

        ! Compute number of vertices
        !===========================
        ! Compute mass of each segment
        allocate(Mi(size(dx)))
        do i = 1, size(dx)
            Mi(i) = sum(rhoi(i, :)*sum(coef, 2)*dl(i))
        end do

        ! Compute total mass
        Mtot = sum(Mi)

        ! Determine total number of vertices (at least two)
        nv = max(ceiling(Mtot), 2)

        ! Determine vertex distribution
        !==============================
        ! Determine mass distribution 
        allocate(Mdistr(size(xc)))
        Mdistr = 0
        do i = 1, size(xc)-1
            Mdistr(i+1) = Mdistr(i) + Mi(i)
        end do 
        Mdistr(size(Mdistr)) = Mtot 

        ! Compute distribution
        distr = real([(k, k = 0, nv-1)], kind=R8)*(Mtot/(nv-1))
        distr(size(distr)) = Mtot ! just to make sure
        distr(1) = 0

        ! Distribute - simply interpolate linearly for now
        if (present(xv) .and. present(yv)) then 
            call DistributeVerticesLine(xc, yc, Mdistr, distr, xv, yv)
        end if 

        ! Return length distribution
        if (present(ldistr)) then 
            ! This has to be the length distribution!
            dll = sqrt(dx**2 + dy**2)
            dllc = spread(0, 1, size(xc))
            dlc = spread(0, 1, size(xc))
            do k = 1, size(xc)-1
                dllc(k+1) = dllc(k) + dll(k)
            end do 
            do k = 1, size(xc)-1
                dlc(k+1) = dlc(k) + dl(k)
            end do 
            ldistr = distr 
            call Interpolate1D(distr, ldistr, Mdistr, dllc)
        end if 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Field based distributor
    subroutine DistributeVerticesDensityBasedOverField(vd, xc, yc, field, nv, &
        xv, yv, ldistr)

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
        class(DensityBasedVertexDistributor2DUDT)   :: vd 
        class(DistributionFunctionUDT), intent(in)  :: field 
        real(R8), intent(in)                :: xc(:), yc(:)
        integer(I8), intent(out)            :: nv 
        real(R8), allocatable, intent(out), optional    :: ldistr(:), xv(:), yv(:)

        ! Auxiliary
        real(R8)                            :: l, signcor, Mtot 
        real(R8), allocatable, dimension(:) :: dx, dy, dl, distr, fc, &
            dllc, dlc, dll, temp, df, Mi, Mdistr
        real(R8), allocatable, dimension(:, :)  :: xi, yi, rhoi 
        integer(I8), allocatable, dimension(:)  :: pxi 

        ! Loop
        integer(I8)                         :: i, k 

        ! Initialize
        !===========
        ! Checks
        if (size(xc) /= size(yc)) then 
            call gdErrorHandler('DistributeUniformOverCurve: curve ' // & 
                'coordinates have incompatible dimensions')
        end if 
        if (present(xv)) then 
            if (allocated(xv)) then 
                deallocate(xv)
            end if 
        end if 
        if (present(yv)) then 
            if (allocated(yv)) then 
                deallocate(yv)
            end if 
        end if

        ! Associate
        associate(&
            xb          => vd%xi, &
            rho         => vd%densityfunction, &
            coef        => vd%lagcoef,  &
            intcoef     => vd%intlagcoef) 

        ! Precompute
        !===========
        ! Compute curve quantities
        dx = xc(2:) - xc(1:size(xc)-1)
        dy = yc(2:) - yc(1:size(yc)-1) 
        allocate(fc(size(xc)))
        call field%Evaluate(xc, yc, fc) 
        df = abs(fc(2:) - fc(1:size(fc)-1))
        dl = sqrt(dx**2 + dy**2)
        l = sum(dl)

        ! Check
        if (all(df <= 0)) then 
            signcor = -1.0_R8
        elseif (all(df >= 0)) then 
            signcor = 1.0_R8
        else
            ! Non-monotonic, print warning
            print *, 'DistributeVerticesDensityBasedOverField: ' // & 
                'field distribution is non-monotonic'
            signcor = 1.0_R8
        end if 

        ! Compute lagrange interpolation 
        xi = spread(xb, 1, size(dx))
        yi = spread(xb, 1, size(dy))
        do i = 1, size(dx)
            xi(i, :) = xi(i, :)*dx(i) + xc(i)
            yi(i, :) = yi(i, :)*dy(i) + yc(i)
        end do 
        allocate(temp(size(xi)))
        call rho%Evaluate(reshape(xi, [size(xi)]), reshape(yi, [size(yi)]), temp)
        rhoi = reshape(temp, [size(xi, 1), size(xi, 2)])
        pxi = [(k, k = 0, vd%order)]

        ! Compute number of vertices
        !===========================
        ! Compute mass of each segment
        allocate(Mi(size(dx)))
        do i = 1, size(dx)
            Mi(i) = sum(rhoi(i, :)*reshape(coef, [size(coef)])*df(i))/sum(df)
        end do

        ! Correct for sign
        Mi = Mi*signcor

        ! Compute total mass
        Mtot = sum(Mi)

        ! Determine total number of vertices (at least two)
        nv = max(ceiling(abs(Mtot)), 2)

        ! Determine vertex distribution
        !==============================
        ! Determine mass distribution 
        allocate(Mdistr(size(xc)))
        Mdistr = 0
        do i = 1, size(xc)-1
            Mdistr(i+1) = Mdistr(i) + Mi(i)
        end do 

        ! Compute distribution
        distr = real([(k, k = 0, nv-1)], kind=R8)*(Mtot/(nv-1))
        distr(size(distr)) = Mtot ! just to make sure

        ! Distribute - simply interpolate linearly for now
        if (present(xv) .and. present(yv)) then 
            call DistributeVerticesLine(xc, yc, Mdistr, distr, xv, yv)
        end if 

        ! Return length distribution
        if (present(ldistr)) then 
            ! This has to be the length distribution!
            dll = sqrt(dx**2 + dy**2)
            dllc = spread(0, 1, size(xc))
            dlc = spread(0, 1, size(xc))
            do k = 1, size(xc)-1
                dllc(k+1) = dllc(k) + dll(k)
            end do 
            do k = 1, size(xc)-1
                dlc(k+1) = dlc(k) + dl(k)
            end do 
            ldistr = distr 
            call Interpolate1D(distr, ldistr, Mdistr, dllc)
        end if 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    subroutine DistributeVerticesDensityBasedOverField2(vd, xc, yc, field, nv, &
        xv, yv, ldistr)

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
        class(DensityBasedVertexDistributor2DUDT)   :: vd 
        class(DistributionFunctionUDT), intent(in)  :: field 
        real(R8), intent(in)                :: xc(:), yc(:)
        integer(I8), intent(out)            :: nv 
        real(R8), allocatable, intent(out), optional    :: ldistr(:), xv(:), yv(:)

        ! Auxiliary
        integer(I8)                         :: nc 
        real(R8)                            :: l, Mtot 
        real(R8), allocatable, dimension(:) :: dx, dy, dl, distr, fc, &
            dllc, dlc, dll, temp, df, Mi, Mdistr, dfdx, dfdy, fx, fy, fn, &
            xf, yf
        real(R8), allocatable, dimension(:, :)  :: xi, yi, rhoi 
        integer(I8), allocatable, dimension(:)  :: pxi 
        character(:), allocatable               :: lengthtype

        ! Loop
        integer(I8)                         :: i, k 

        ! Initialize
        !===========
        nc = size(xc)

        ! Checks
        if (size(xc) /= size(yc)) then 
            call gdErrorHandler('DistributeUniformOverCurve: curve ' // & 
                'coordinates have incompatible dimensions')
        end if 
        if (present(xv)) then 
            if (allocated(xv)) then 
                deallocate(xv)
            end if 
        end if 
        if (present(yv)) then 
            if (allocated(yv)) then 
                deallocate(yv)
            end if 
        end if

        ! Associate
        associate(&
            xb          => vd%xi, &
            rho         => vd%densityfunction, &
            coef        => vd%lagcoef,  &
            intcoef     => vd%intlagcoef) 

        ! Precompute
        !===========
        ! Compute curve quantities
        dx = xc(2:) - xc(1:nc-1)
        dy = yc(2:) - yc(1:nc-1) 
        dl = sqrt(dx**2 + dy**2)
        l = sum(dl)

        ! Compute field quantities
        allocate(fc(nc), dfdx(nc), dfdy(nc), fx(nc), fy(nc))
        call field%Evaluate(xc, yc, fc)
        call field%EvaluateDerivative(xc, yc, 1, 0, dfdx)
        call field%EvaluateDerivative(xc, yc, 0, 1, dfdy)
        df = abs(fc(2:) - fc(1:size(fc)-1))
        fn = sqrt(dfdx**2 + dfdy**2) 
        fx = 0
        fy = 0
        where( fn > disttol) 
            fx = dfdx/fn
            fy = dfdy/fn 
        end where 

        ! Recompute the length coordinate
        fc = 0
        lengthtype = 'radial'
        select case (lengthtype)

        case ('psi')

            do i = 2, nc
                fc(i) = fc(i-1) + df(i-1)
            end do
            
        case ('euler')

            do i = 2, nc
                fc(i) = fc(i-1) + dl(i-1)
            end do
            df = dl 

        case ('radial')

            do i = 2, nc
                df(i-1) = abs(dx(i-1)*fx(i-1) + dy(i-1)*fy(i-1))
                fc(i) = fc(i-1) + df(i-1)
            end do
            

        case default 

            call gdErrorHandler('DistributeVerticesDensityBasedOverField:' // & 
                'unknown length type')

        end select

        ! Compute density at face centers
        xf = 0.5*(xc(1:nc-1) + xc(2:nc))
        yf = 0.5*(yc(1:nc-1) + yc(2:nc))
        xi = spread(xb, 1, size(dx))
        yi = spread(xb, 1, size(dy))
        do i = 1, size(dx)
            xi(i, :) = xi(i, :)*dx(i) + xc(i)
            yi(i, :) = yi(i, :)*dy(i) + yc(i)
        end do 
        allocate(temp(nc-1))
        call rho%Evaluate(xf, yf, temp)

        ! Compute number of vertices
        !===========================
        ! Compute mass of each segment
        allocate(Mi(nc-1))
        do i = 1, nc-1
            Mi(i) = temp(i)*df(i)
        end do

        ! Compute total mass
        Mtot = sum(Mi)

        ! Determine total number of vertices (at least two)
        nv = max(ceiling(abs(Mtot)), 2)

        ! Determine vertex distribution
        !==============================
        ! Determine mass distribution 
        allocate(Mdistr(nc))
        Mdistr = 0
        do i = 1, nc-1
            Mdistr(i+1) = Mdistr(i) + Mi(i)
        end do 

        ! Compute distribution
        distr = real([(k, k = 0, nv-1)], kind=R8)*(Mtot/(nv-1))
        distr(size(distr)) = Mtot ! just to make sure

        ! Distribute - simply interpolate linearly for now
        if (present(xv) .and. present(yv)) then 
            call DistributeVerticesLine(xc, yc, Mdistr, distr, xv, yv)
        end if 

        ! Return length distribution
        if (present(ldistr)) then 
            ! This has to be the length distribution!
            dll = sqrt(dx**2 + dy**2)
            dllc = spread(0, 1, nc)
            dlc = spread(0, 1, nc)
            do k = 1, nc-1
                dllc(k+1) = dllc(k) + dll(k)
            end do 
            do k = 1, nc-1
                dlc(k+1) = dlc(k) + dl(k)
            end do 
            ldistr = distr 
            call Interpolate1D(distr, ldistr, Mdistr, dllc)
        end if 

        ! Housekeeping
        !=============
        end associate

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