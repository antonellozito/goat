!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains some small and easy subroutines to create the
! grid coordinates of structured 2D grids. 

module mod_structured2Dgridding

    ! Declare modules
    !================
    use mod_precision 
    use mod_errorhandler
    use mod_sort

    implicit none 
    private 
    save
    public Construct2DStructuredGrid, Construct2DStructuredUniformGrid, &
        ConstructRefined2DStructuredGrid

    ! General precision parameters
    real(R8), parameter, private    :: disttol = 1e-12

    ! No derived types

    contains 

    ! Structured non-uniform 2D grid, no refinement
    subroutine Construct2DStructuredGrid(xgv, ygv, nx, ny, xg, yg)

        ! Description
        !============ 
        ! This routine returns the grid coordinates (in vector form!) of a
        ! 2D structured grid defined by the gridding vectors xgv and ygv.
        ! It is assumed that these gridding vectors are monotonic, though
        ! this is not explicitly checked. 
    
        ! Initialize
        !===========
        ! Modules
        use mod_precision
    
        ! The usual
        implicit none
    
        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)             :: nx, ny 
        real(R8), intent(in)                :: xgv(1:nx), ygv(1:ny)
        real(R8), intent(out)               :: xg(1:nx*ny), yg(1:nx*ny)
    
        ! Auxiliary
        real(R8), allocatable               :: xaux(:, :), yaux(:, :)
    
        ! Compute grid
        !=============
        ! Allocate
        allocate(xaux(nx, ny))
        allocate(yaux(nx, ny))
    
        ! Replicate
        xaux(:,:) = spread(xgv, 2, ny)
        yaux(:,:) = spread(ygv, 1, nx)
    
        ! Reshape
        xg = reshape(xaux, (/(nx*ny)/) ) 
        yg = reshape(yaux, (/(nx*ny)/) )
    
        ! Deallocate
        deallocate(xaux, yaux)
    
    end subroutine

    ! Structured non-uniform 2D grid with possible refinement near given
    ! points
    subroutine ConstructRefined2DStructuredGrid(xg, yg, xgv, ygv, &
        xb, yb, resx, resy, xp, yp, padx, pady, dxfracmin, dyfracmin)

        ! Description
        !============
        ! This function constructs a 2D grid using meshgrid starting from the
        ! bounds given in xb, yb, the desired resolution in x and y direction resx,
        ! resy, and any specific points to be introduced into the gridding vectors.
        ! It returns the coordinate arrays x, y as outputted by meshgrid and the
        ! constructing coordinate vectors xgv, ygv. padx and pady should be
        ! positive integer numbers that indicate how many additional coordinates
        ! should be added before and after an inserted point.

        ! This function is particularly useful to deal with insufficient resolution
        ! near saddle points when tracing contour lines without having to (overly)
        ! refine the domain (depending of course on the number of saddle points and
        ! locations). 

        ! It is, however, important to note that not exactly resx and resy amount
        ! of vertices will be added in both directions, but that this will be
        ! larger depending on the number of saddle points (coordinates xp and yp)
        ! and the desired padding padx, pady. 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), allocatable, intent(out)      :: xg(:), yg(:), &
            xgv(:), ygv(:)
        real(R8), intent(in)                    :: xb(1:2), yb(1:2), &
            xp(:), yp(:)
        real(R8), intent(inout)                 :: dxfracmin, dyfracmin
        integer(I8), intent(in)                 :: resx, resy, padx, pady 

        ! Auxiliary
        real(R8)                                :: minx, miny, maxx, &
            maxy, sx, sy, dxmin, dymin
        real(R8), allocatable, dimension(:)     :: xgvp, &
            ygvp, dxgvp, dygvp, xgvpnew, ygvpnew
        integer(I8)                             :: nxp, nyp
        integer(I8), allocatable, dimension(:)  :: nsx, nsy
        logical, allocatable, dimension(:)      :: keep 

        ! Loop
        integer(I8)                             :: i, k, xc, yc 

        ! Initialize
        !===========
        ! Check input
        nxp = size(xp)
        nyp = size(yp)

        if (dxfracmin > 0.1_R8) then 
            ! This seems much to high, cap
            print *, 'given dxfracmin: ', dxfracmin
            print *, 'ConstructRefined2DStructuredGrid: dxfracmin ' // & 
                'unrealistically high, capping to 0.1'
            dxfracmin = 0.1_R8
        end if 
        if (dyfracmin > 0.1_R8) then 
            ! This seems much to high, cap
            print *, 'given dyfracmin: ', dyfracmin
            print *, 'ConstructRefined2DStructuredGrid: dyfracmin ' // & 
                'unrealistically high, capping to 0.1'
            dyfracmin = 0.1_R8
        end if 
            
    

        ! Construct gridding vectors
        !===========================
        ! Compute domain metrics
        minx = minval(xb) 
        maxx = maxval(xb)
        miny = minval(yb) 
        maxy = maxval(yb)

        ! Range
        sx = maxx - minx
        sy = maxy - miny

        ! Minimal distance
        dxmin = sx*dxfracmin 
        dymin = sy*dyfracmin

        ! Add additional points & padding in x-direction
        if (size(xp) > 0) then 
            ! Initialize
            xgvp = [minx, xp, maxx]
            call Sort(xgvp)
            
            ! Compute distance between seeding points
            dxgvp = xgvp(2:nxp+2) - xgvp(1:nxp+1)

            ! Check for almost zero distance, eliminate
            keep = dxgvp > dxmin
            dxgvp = pack(dxgvp, keep)  
            xgvp = pack(xgvp, [keep, .true.])

            ! Add padding points around additional points (except minx, maxx)
            allocate(xgvpnew(3*(size(xgvp)-2)+2))
            xgvpnew(1) = xgvp(1)
            xgvpnew(size(xgvpnew)) = xgvp(size(xgvp))
            k = 1
            do i = 2, size(xgvp)-1
                xgvpnew(k+1:k+3) = [-2*padx*dxmin, 0.0_R8, 2*padx*dxmin] + xgvp(i)
                k = k + 3
            end do  
            xgvp = xgvpnew
            call Sort(xgvp, ascend=.true.)
            dxgvp = xgvp(2:size(xgvp)) - xgvp(1:size(xgvp)-1)
            
            ! Compute number of segments
            nsx = floor(dxgvp/(sx/resx))+1
            
            ! Adjust for desired padding
            where(nsx < 2*padx) nsx = 2*padx 
            
            ! Recompute gridding vector
            allocate(xgv(sum(nsx)+1))
            xgv = 0
            xc = 0 !counter
            do i = 1, size(nsx)
                xgv(xc+1:xc+nsx(i)+1) = [(k, k = 0, nsx(i))]*&
                    (xgvp(i+1) - xgvp(i))/(nsx(i)) + xgvp(i)
                xc = xc + nsx(i)
            end do 
        else 
            ! Uniform gridding vector
            xgv = real([(k, k = 0, resx)], kind=R8)*sx + minx 
            
        end if

        ! Add additional points & padding in y-direction
        if (size(yp) > 0) then 
            ! Initialize
            ygvp = [miny, yp, maxy]
            call Sort(ygvp)
            
            ! Compute distance between seeding points
            dygvp = ygvp(2:nyp+2) - ygvp(1:nyp+1)
            
            ! Check for almost zero distance, eliminate
            keep = dygvp > dymin
            dygvp = pack(dygvp, keep)  
            ygvp = pack(ygvp, [keep, .true.])

            ! Add padding points around additional points (except minx, maxx)
            allocate(ygvpnew(3*(size(ygvp)-2)+2))
            ygvpnew(1) = ygvp(1)
            ygvpnew(size(ygvpnew)) = ygvp(size(ygvp))
            k = 1
            do i = 2, size(ygvp)-1
                ygvpnew(k+1:k+3) = [-2*pady*dymin, 0.0_R8, 2*pady*dymin] + ygvp(i)
                k = k + 3
            end do  
            ygvp = ygvpnew
            call Sort(ygvp, ascend=.true.)
            dygvp = ygvp(2:size(ygvp)) - ygvp(1:size(ygvp)-1)

            ! Compute number of segments
            nsy = floor(dygvp/(sy/resy))+1
            
            ! Adjust for desired padding
            where(nsy < 2*pady) nsy = 2*pady 
            
            ! Recompute gridding vector
            allocate(ygv(sum(nsy)+1))
            ygv = 0
            yc = 0 !counter
            do i = 1, size(nsy)
                ygv(yc+1:yc+nsy(i)+1) = [(k, k = 0, nsy(i))]*&
                    (ygvp(i+1) - ygvp(i))/(nsy(i)) + ygvp(i)
                yc = yc + nsy(i)
            end do 
        else 
            ! Uniform gridding vector
            ygv = real([(k, k = 0, resy)], kind=R8)*sy + miny
            
        end if

        ! Construct grid
        allocate(xg(size(xgv)*size(ygv)), yg(size(xgv)*size(ygv)))
        call Construct2DStructuredGrid(xgv, ygv, size(xgv), size(ygv), &
            xg, yg)

    end subroutine 

    ! Structured uniform 2D grid, no refinement
    subroutine Construct2DStructuredUniformGrid(xg, yg, xgv, ygv, &
        xrange, yrange, nx, ny, offsetx, offsety)
    
        ! Description
        !============ 
        ! This routine returns the grid coordinates (in vector form!) of a
        ! 2D structured uniform grid. The extent of the grid is determined
        ! by the minimal and maximal values in xrange and yrange and by the
        ! offsets offsetx and offsety. The latter are relative offsets 
        ! compared to the initial length and width of the domain. 
    
        ! Initialize
        !===========
        ! Modules
        use mod_precision
    
        ! The usual
        implicit none
    
        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)                :: xrange(:), yrange(:), &
            offsetx, offsety
        integer(I8), intent(in)             :: nx, ny 
        real(R8), intent(out)               :: xg(1:nx*ny), yg(1:nx*ny), &
            xgv(1:nx), ygv(1:ny)
    
        ! Auxiliary
        real(R8)                            :: maxx, minx, maxy, miny, &
            dx, dy 
    
        ! Loop
        integer(I8)                         :: k 
    
        ! Compute gridding vectors
        !=========================
        ! Compute domain extent
        maxx = maxval(xrange)
        minx = minval(xrange)
        maxy = maxval(yrange)
        miny = minval(yrange)
        dx   = maxx - minx 
        dy   = maxy - miny 
    
        ! Extend using offset
        maxx = maxx + offsetx*dx 
        minx = minx - offsetx*dx 
        maxy = maxy + offsety*dy 
        miny = miny - offsety*dy 
    
        ! Recompute extent
        dx   = maxx - minx 
        dy   = maxy - miny
    
        ! Construct gridding vectors
        xgv = [(k, k = 0, nx-1)]*dx/(nx-1) + minx 
        ygv = [(k, k = 0, ny-1)]*dy/(ny-1) + miny 
    
        ! Compute grid
        !=============
        call Construct2DStructuredGrid(xgv, ygv, nx, ny, xg, yg)
    
    end subroutine



end module 