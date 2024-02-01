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