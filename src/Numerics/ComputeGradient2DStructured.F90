subroutine ComputeGradient2DStructured(v, xg, yg, nx, ny, vx, vy)

    ! Description
    !============
    ! This routine computes the gradient on a structured (non-uniform)
    ! quadrilateral grid. xg and yg are assumed to be the grid 
    ! coordinates, where nx and ny are the grid dimensions (i.e. xg is a
    ! nx-by-ny array). v are the values (same size as xg) of which the
    ! gradient has to be computed. vx and vy are the gradients in x-
    ! and y-direction, respectively, and are of the same size as xg.

    ! Note: in the bulk of the domain, a central difference scheme is 
    ! used, while at the edges, a forward or backward difference scheme
    ! is applied. 

    ! Initialize
    !===========
    ! Declare modules
    use mod_precision 

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    integer(I8), intent(in)         :: nx, ny
    real(R8), intent(in)            :: xg(1:nx, 1:ny), yg(1:nx, 1:ny), &
        v(1:nx, 1:ny)
    real(R8), intent(out)           :: vx(1:nx, 1:ny), vy(1:nx, 1:ny)

    ! Loop variables

    ! Auxiliary variables 
    real(R8)                        :: vxFW(1:nx, 1:ny), &
        vxBW(1:nx, 1:ny), vyFW(1:nx, 1:ny), vyBW(1:nx, 1:ny)

    ! Initialize
    !===========
    vxFW(:,:) = 0
    vxBW(:,:) = 0
    vyFW(:,:) = 0
    vyBW(:,:) = 0

    ! Gradient, x-direction
    !======================
    ! Forward difference
    vxFW(1:nx-1,:) = (v(2:nx,:) - v(1:nx-1,:))/(xg(2:nx,:) - xg(1:nx-1,:))

    ! Backward difference
    vxBW(2:nx,:) = (v(1:nx-1,:) - v(2:nx,:))/(xg(1:nx-1,:) - xg(2:nx,:))

    ! Central difference (bulk, simple average)
    vx(2:nx-1,:) = 0.5*(vxFW(2:nx-1,:) + vxBW(2:nx-1,:))
    vx(1,:) = vxFW(1,:)
    vx(nx,:) = vxBW(nx,:)

    ! Gradient, y-direction
    !======================
    ! Forward differenc-e
    vyFW(:,1:ny-1) = (v(:,2:ny) - v(:,1:ny-1))/(yg(:,2:ny) - yg(:,1:ny-1))

    ! Backward difference
    vyBW(:,2:ny) = (v(:,1:ny-1) - v(:,2:ny))/(yg(:,1:ny-1) - yg(:,2:ny))

    ! Central difference (bulk, simple average)
    vy(:,2:ny-1) = 0.5*(vyFW(:,2:ny-1) + vyBW(:,2:ny-1))
    vy(:,1) = vyFW(:,1)
    vy(:,ny) = vyBW(:,ny)

end subroutine