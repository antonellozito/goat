!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains some functionality to make plots during runtime 
! using the unix GNUplot functionality. Most (if not all) plots are made
! by writing the data to a text file (depending on which plot to make)
! which may become computtionally costly if done e.g. every iteration
! in a loop. The  plotting routines from the mod_plotter mod are 
! used to do the basic plotting.

! The idea is that every plotting routine in this module has a .plt 
! counterpart in /src/Visualization, which can be called using gnuplot. 
! In order to not close the window after gnuplot has finished 
! (in case the x11 terminal output is chosen), the 'persistent' option 
! should be parsed. The data for the plot has to be written to a data 
! file, also in Visualization, which should have the same name as the
! function, but with a .dat extension. 

! IMPORTANT: every time this function is called, a new window or 
! figure is created - make sure not to call these in a loop in the main 
! code! 

module gdmod_plots

    ! Initialize
    !============
    ! Load modules
    use mod_plotter
    use gdmod_types
    use gdmod_userinput 
    use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN

    ! The usual
    implicit none
    save
    public 

    contains 

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                            Main plotters                         !
    !------------------------------------------------------------------!

    ! 2D unstructured state field
    subroutine Plot2DUnstructuredField(field, grid, loc, gnuplotoptions)

        ! Description
        !============
        ! Plot a cell-centered 2D field. The field should be a nc-by-1 
        ! array (with nc the number of cells in grid%cells). This 
        ! routine calls the Patchplot2D routine to do the plotting. 

        ! The usual
        implicit none

        ! Arguments
        real(R8)                            :: field(:)
        type(GridUDT), intent(in)           :: grid 
        character(*)                        :: gnuplotoptions
        character(1)                        :: loc ! location of where to plot 

        ! Loop variables
        integer(I8)                         :: j, k, ind

        ! Auxiliary variables
        integer(I8)                         :: np, vertind
        real(R8), allocatable               :: xp(:), yp(:), val(:)

        ! Set NaN
        real(R8)                            :: NaN

        ! Initialize
        !===========
        ! Set NaN
        NaN = IEEE_VALUE(nan, IEEE_QUIET_NAN)

        ! Construct plot data
        !====================
        ! Check the location 
        select case (loc)

        case ('c')

            ! Cell based quantities

            ! Length of xp, yp
            np = sum(grid%cells%vertP(:,2)) + grid%cells%ntot

            ! Allocate
            allocate(xp(np))
            allocate(yp(np))
            allocate(val(np))

            ! Loop over all cells
            ind = 1
            do k = 1, grid%cells%ntot
                ! Loop over all vertices of the current cell
                do j = 1, grid%cells%vertP(k,2)
                    ! Add the coordinates and value
                    vertind = grid%cells%vertP(k,1)+j-1
                    xp(ind) = grid%vert%x(vertind)
                    yp(ind) = grid%vert%y(vertind)
                    val(ind) = field(k)

                    ! Update ind
                    ind = ind+1
                end do
                
                ! Add NaN
                xp(ind) = NaN
                yp(ind) = NaN
                val(ind) = NaN

                ! Update ind
                ind = ind+1
            end do

        case ('v')
            
            ! Vertex based quantities (still plotted per cell)

            ! Length of xp, yp
            np = sum(grid%cells%vertP(:,2)) + grid%cells%ntot

            ! Allocate
            allocate(xp(np))
            allocate(yp(np))
            allocate(val(np))

            ! Loop over all cells
            ind = 1
            do k = 1, grid%cells%ntot
                ! Loop over all vertices of the current cell
                do j = 1, grid%cells%vertP(k,2)
                    ! Add the coordinates and value
                    vertind = grid%cells%vertlist(grid%cells%vertP(k,1)+j-1)
                    xp(ind) = grid%vert%x(vertind)
                    yp(ind) = grid%vert%y(vertind)
                    val(ind) = field(vertind)

                    ! Update ind
                    ind = ind+1
                end do
                
                ! Add NaN
                xp(ind) = NaN
                yp(ind) = NaN
                val(ind) = NaN

                ! Update ind
                ind = ind+1
            end do

        case default

            call gdErrorHandler('Plotter: plot option not implemented')

        end select

        ! Call the plotter
        !=================
        print *, size(xp,1), size(yp,1), size(val,1), val(1), np, ind
        call Patchplot2D(xp, yp, val, np, gnuplotoptions)


    end subroutine

    !------------------------------------------------------------------!
    !                              Grid plots                          !
    !------------------------------------------------------------------!

    ! Grid 
    !=====
    ! Grid (nodes)
    subroutine PlotGrid(grid,gnuplotoptions)

        ! Description
        !============
        ! Make a plot of the grid for inspection. 

        ! The usual
        implicit none

        ! Declare variables
        type(GridUDT)                       :: grid
        integer                             :: i, fu
        character(*)                        :: gnuplotoptions

        ! Initialize
        !===========
        ! Set the correct directories
        call SetGnuplotNames(plotfile,datafile,'plotgrid')

        ! Write the data file
        !====================
        ! Write vertex coordinates to file
        open (action='write', file=trim(datafile), newunit=fu, &
             status='replace')
    
        do i = 1, grid%vert%ntot
            write (fu, *) grid%vert%x(i), grid%vert%y(i)
        end do
    
        close (fu)

        ! Call plotter
        !=============
        call gnuplotexe(gnuplotoptions,trim(plotfile))

    end subroutine

    ! Cells
    subroutine PlotGridCells(grid,gnuplotoptions)

        ! Description
        !============
        ! Make a plot of the grid cells for inspection. 

        ! The usual
        implicit none

        ! Declare variables
        type(GridUDT)                       :: grid
        integer                             :: fu
        character(*)                        :: gnuplotoptions

        integer(I8)                         :: i, j, nvpc, si
        integer(I8), allocatable            :: tcv(:)

        ! Initialize
        !===========
        ! Set the correct directories
        call SetGnuplotNames(plotfile,datafile,'plotgridcells')

        ! Write the data file
        !====================
        ! Write vertex coordinates to file
        open (action='write', file=trim(datafile), newunit=fu, &
             status='replace')
    
        ! Loop over all cells
        do i = 1, grid%cells%ntot
            ! Write vertex coordinates as ordened in the cell
            si = grid%cells%vertP(i,1) ! start index
            nvpc = grid%cells%vertP(i,2) ! number of vertices per cell

            ! Allocate
            allocate(tcv(nvpc))

            ! Get the vertex indices of the current cell
            tcv = grid%cells%vertlist(si:si+nvpc-1)
            do j = 1, nvpc
                ! Print
                write (fu, *) grid%vert%x(tcv(j)), grid%vert%y(tcv(j))
            end do 
            write(fu, *) ! leave blank line between each cell

            ! Deallocate
            deallocate(tcv)

        end do
    
        close (fu)

        ! Call plotter
        !=============
        call gnuplotexe(gnuplotoptions,trim(plotfile))

    end subroutine

    ! Faces
    subroutine PlotGridFaces(grid,gnuplotoptions)

        ! Description
        !============
        ! Make a plot of the grid faces for inspection. 

        ! The usual
        implicit none

        ! Declare variables
        type(GridUDT)                       :: grid
        integer                             :: fu
        character(*)                        :: gnuplotoptions

        integer(I8)                         :: i

        ! Initialize
        !===========
        ! Set the correct directories
        call SetGnuplotNames(plotfile,datafile,'plotgridfaces')

        ! Write the data file
        !====================
        ! Write vertex coordinates to file
        open (action='write', file=trim(datafile), newunit=fu, &
             status='replace')
    
        ! Loop over all faces
        do i = 1, grid%faces%ntot
            ! Write vertex coordinates of the face
            write(fu, *) grid%vert%x(grid%faces%vert(i,1)), &
                grid%vert%y(grid%faces%vert(i,1))
            write(fu, *) grid%vert%x(grid%faces%vert(i,2)), &
                grid%vert%y(grid%faces%vert(i,2))
            write(fu, *) ! blank line

        end do
    
        close (fu)

        ! Call plotter
        !=============
        call gnuplotexe(gnuplotoptions,trim(plotfile))

    end subroutine

    ! Flux surfaces as vertices and numbers
    subroutine PlotFluxSurfaces(grid,gnuplotoptions)

        ! Description
        !============
        ! Plot all the flux surfaces, including the ID of each flux 
        ! surface vertex

        ! The usual
        implicit none

        ! Declare variables
        type(GridUDT)                       :: grid
        integer                             :: i, fu
        character(*)                        :: gnuplotoptions

        ! Initialize
        !===========
        ! Set the correct directories
        call SetGnuplotNames(plotfile,datafile,'plotfluxsurfaces')

        ! Write the data file
        !====================
        ! Write vertex coordinates and their surface indices to file
        open (action='write', file=trim(datafile), newunit=fu, &
             status='replace')
    
        do i = 1, grid%vert%ntot
            write (fu, *) grid%vert%x(i), grid%vert%y(i), &
                grid%data%fluxdata%fluxsurfaceID(i)
        end do
    
        close (fu)

        ! Call plotter
        !=============
        call gnuplotexe(gnuplotoptions,trim(plotfile))

    end subroutine

    ! Boundaries
    subroutine PlotGridBoundaries(grid,gnuplotoptions)

        ! Description
        !============
        ! Plot all the different boundaries according to the provided
        ! labels. Will loop over boundaries and 
        ! plot the faces. Face coordinates are constructed
        ! locally for plotting starting from the vertex coordinates and
        ! taking a simple arithmetic average. 

        ! The usual
        implicit none

        ! Declare variables
        type(GridUDT)                       :: grid
        integer                             :: i, j, fu, tnf, nb
        real(R8), allocatable               :: fx(:), fy(:), tfx(:), tfy(:)
        logical, allocatable                :: mask(:)
        character(*)                        :: gnuplotoptions

        ! Initialize
        !===========
        ! Set the correct directories
        call SetGnuplotNames(plotfile,datafile,'plotgridboundaries')

        ! Write the data file
        !====================
        ! Get number of boundaries
        nb = size(grid%bnd)

        ! Allocate
        allocate(fx(grid%faces%ntot))
        allocate(fy(grid%faces%ntot))
        allocate(mask(grid%faces%ntot))

        ! Compute coordinates
        fx(:) = (grid%vert%x(grid%faces%vert(:,1)) & 
            + grid%vert%x(grid%faces%vert(:,2)))*0.5
        fy(:) = (grid%vert%y(grid%faces%vert(:,1)) & 
            + grid%vert%y(grid%faces%vert(:,2)))*0.5

        ! Write vertex coordinates and their surface indices to file
        open (action='write', file=trim(datafile), newunit=fu, &
             status='replace')
    
        do i = 1, nb
            ! Get all faces 
            tnf = size(grid%bnd(i)%faces)

            ! Allocate coordinate vectors, add coordinates
            allocate(tfx(tnf))
            allocate(tfy(tnf))
            print *, tnf 

            tfx = fx(grid%bnd(i)%faces)
            tfy = fy(grid%bnd(i)%faces)

            ! Write
            do j = 1, tnf
                write (fu, *) tfx(j), tfy(j), grid%bnd(i)%ID
            end do

            ! Write blank line
            write (fu, *) 

            ! Deallocate coordinates
            deallocate(tfx)
            deallocate(tfy)
            
        end do
    
        close (fu)

        ! Call plotter
        !=============
        call gnuplotexe(gnuplotoptions,trim(plotfile))

    end subroutine

    ! Face labels
    subroutine PlotFaceLabels(grid,labels,gnuplotoptions)

        ! Description
        !============
        ! Plot all the different boundaries according to the provided
        ! face labels. Will loop over face labels and check which faces
        ! have the specified label. Face coordinates are constructed
        ! locally for plotting starting from the vertex coordinates and
        ! taking a simple arithmetic average. 

        ! The usual
        implicit none

        ! Declare variables
        type(GridUDT)                       :: grid
        integer                             :: i, j, fu, nl, tnf
        integer                             :: labels(:)
        real(R8), allocatable               :: fx(:), fy(:), tfx(:), tfy(:)
        logical, allocatable                :: mask(:)
        character(*)                        :: gnuplotoptions

        ! Initialize
        !===========
        ! Set the correct directories
        call SetGnuplotNames(plotfile,datafile,'plotfacelabels')

        ! Write the data file
        !====================
        ! Get number of labels
        nl = size(labels)

        ! Allocate
        allocate(fx(grid%faces%ntot))
        allocate(fy(grid%faces%ntot))
        allocate(mask(grid%faces%ntot))

        ! Compute coordinates
        fx(:) = (grid%vert%x(grid%faces%vert(:,1)) & 
            + grid%vert%x(grid%faces%vert(:,2)))*0.5
        fy(:) = (grid%vert%y(grid%faces%vert(:,1)) & 
            + grid%vert%y(grid%faces%vert(:,2)))*0.5

        ! Write vertex coordinates and their surface indices to file
        open (action='write', file=trim(datafile), newunit=fu, &
             status='replace')
    
        do i = 1, nl
            ! Get all faces with this label
            mask = grid%data%regions%facelabel == labels(i)
            tnf = count(mask)

            ! Allocate coordinate vectors, add coordinates
            allocate(tfx(tnf))
            allocate(tfy(tnf))

            tfx = pack(fx, mask)
            tfy = pack(fy, mask)

            ! Write
            do j = 1, tnf
                write (fu, *) tfx(j), tfy(j), labels(i)
            end do

            ! Write blank line
            write (fu, *) 

            ! Deallocate coordinates
            deallocate(tfx)
            deallocate(tfy)
            
        end do
    
        close (fu)

        ! Call plotter
        !=============
        call gnuplotexe(gnuplotoptions,trim(plotfile))

    end subroutine

    ! Data with grid on background
    subroutine PlotGridWithPoints(grid, xp, yp, gnuplotoptions)

        ! Description
        !============
        ! Plot data points on a grid background. First, the 
        ! PlotGridCells routine is called to construct the grid file 
        ! (and a plot that is not persisted - may be cleaned up in the
        ! future), afterwards, the points are plotted. This routine uses
        ! therefore two separate files to store the data. 

        ! The usual
        implicit none

        ! Declare variables
        type(GridUDT)                       :: grid
        integer                             :: i, fu
        character(*)                        :: gnuplotoptions

        real(R8), allocatable, intent(in)   :: xp(:), yp(:)

        ! Plot grid
        !==========
        call PlotGridCells(grid, '')

        ! Initialize
        !===========
        ! Set the correct directories
        call SetGnuplotNames(plotfile,datafile,'plotgridwithpoints')

        ! Write the data file
        !====================
        ! Write point coordinates to file
        open (action='write', file=trim(datafile), newunit=fu, &
             status='replace')
    
        do i = 1, size(xp)
            write (fu, *) xp(i), yp(i)
        end do
    
        close (fu)

        ! Call plotter
        !=============
        call gnuplotexe(gnuplotoptions,trim(plotfile))

    end subroutine

    !------------------------------------------------------------------!
    !                            Optimization                          !
    !------------------------------------------------------------------!

    !------------------------------------------------------------------!
    !                               State                              !
    !------------------------------------------------------------------!

    ! Magnetic flux
    subroutine PlotMagneticFlux(magneticField, gnuplotoptions)

        ! Description
        !============
        ! Plot the magnetic flux values (psi), which are in this case 
        ! given on a 2D (possibly nonuniform) mesh. The general 
        ! patchplot routine is used to visualize the data. 

        ! The usual
        implicit none

        ! Declare variables
        type(MagneticFieldUDT)              :: magneticField
        integer                             :: i, j, k, np
        real(R8), allocatable               :: R(:), Z(:), Psi(:, :), &
                                            x(:), y(:), val(:)
        character(*)                        :: gnuplotoptions

        ! Set NaN
        real(R8)                            :: NaN

        ! Initialize
        !===========
        ! Set NaN
        NaN = IEEE_VALUE(nan, IEEE_QUIET_NAN)

        ! Write the data 
        !===============
        ! Number of polygons to be plotted (one cell each)
        np = (magneticField%nR - 1) * (magneticField%nZ  - 1)

        ! Allocate
        allocate(R(magneticField%nR))
        allocate(Z(magneticField%nZ))
        allocate(Psi(magneticField%nR, magneticField%nZ))

        allocate(x(5*np))
        allocate(y(5*np))
        allocate(val(5*np))

        ! Unpack
        R = magneticField%R
        Z = magneticField%Z
        Psi = magneticField%Psi

        ! Loop over both coordinate directions
        print *, 
        k = 1
        do i = 1, magneticField%nZ-1
            do j = 1, magneticField%nR-1
                ! Write
                x(k) = R(j)
                y(k) = Z(i)
                val(k) = Psi(j, i)
                k = k+1

                x(k) = R(j+1)
                y(k) = Z(i)
                val(k) = Psi(j+1, i)
                k = k+1

                x(k) = R(j+1)
                y(k) = Z(i+1)
                val(k) = Psi(j+1, i+1)
                k = k+1

                x(k) = R(j)
                y(k) = Z(i+1)
                val(k) = Psi(j, i+1)
                k = k+1

                x(k) = NaN
                y(k) = NaN
                val(k) = NaN
                k = k+1

            end do
        end do

        ! Call the plotter
        !=================
        call Patchplot2D(x, y, val, 5*np, gnuplotoptions)

    end subroutine

    !------------------------------------------------------------------!
    !                             Environment                          !
    !------------------------------------------------------------------!
    
    ! Vessel polygon
    subroutine PlotVesselPolygon(vessel, gnuplotoptions)

        ! Description
        !============
        ! Plot the vessel polygon 

        ! The usual
        implicit none

        ! Declare variables
        type(VesselUDT)                     :: vessel
        character(*)                        :: gnuplotoptions
        integer(I8)                         :: i
        integer                             :: fu

        ! Initialize
        !===========
        ! Set the correct directories
        call SetGnuplotNames(plotfile,datafile,'plot2dpolygon')

        ! Write vertex coordinates and their surface indices to file
        open (action='write', file=trim(datafile), newunit=fu, &
             status='replace')

        ! Write the data file
        !====================
        ! Loop over all polygon parts
        do i = 1, vessel%nedges

            if (any(i == vessel%polygonstart)) then 
                ! Write white line
                write(fu, *)
            end if

            ! Write coordinates
            write (fu, *) vessel%x(vessel%edges(i, 1)), vessel%y(vessel%edges(i, 1))
            write (fu, *) vessel%x(vessel%edges(i, 2)), vessel%y(vessel%edges(i, 2))

            ! Insert white line
            !write(fu, *)

        end do
    
        ! Close file
        close (fu)

        ! Call plotter
        !=============
        call gnuplotexe(gnuplotoptions,trim(plotfile))

    end subroutine

    !------------------------------------------------------------------!
    !                         Auxiliary routines                       !
    !------------------------------------------------------------------!
    
end module 