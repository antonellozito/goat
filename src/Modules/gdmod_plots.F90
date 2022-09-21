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
! in a loop. 

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
    use gdmod_types
    use gdmod_userinput 
    use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN

    ! The usual
    implicit none
    save
    public 

    ! Define the (relative) file directory where the plot are located
    character(len=*), parameter     :: plotdir = './src/Visualization/'

    ! Define the (general) plot filename and data filename
    character(C128)                 :: plotfile
    character(C128)                 :: datafile

    contains 

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                            Main plotters                         !
    !------------------------------------------------------------------!

    ! Generic
    !========
    ! 2D polygon plotter (only points)
    subroutine Plot2DPolygon(x,y,n,gnuplotoptions)

        ! Description
        !============
        ! Plot a generic 2D polygon based on x and y coordinates. 
        ! Different polygons can be separated using NaNs, for which a 
        ! blank line will be written to disconnect the polygons in the
        ! plot. 

        ! The usual
        implicit none

        ! Declare variables
        integer                             :: n, i, fu
        real(R8), dimension(n)              :: x, y
        character(*)                        :: gnuplotoptions

        ! Initialize
        !===========
        ! Set the correct directories
        call SetGnuplotNames(plotfile,datafile,'plot2dpolygon')

        ! Write the data file
        !====================
        ! Write vertex coordinates to file
        open (action='write', file=trim(datafile), newunit=fu, &
             status='replace')
    
        do i = 1, n
            if (isnan(x(i))) then 
                ! Write blank line
                write(fu, *)
            else
                ! Write coordinates
                write (fu, *) x(i), y(i)
            end if
        end do
    
        close (fu)

        ! Call plotter
        !=============
        call gnuplotexe(gnuplotoptions,trim(plotfile))

    end subroutine

    ! Patchplot (filled polygons/curves)
    subroutine Patchplot2D(x, y, z, n, gnuplotoptions)

        ! Description
        !============
        ! Make a patchplot of the polygon with coordinates x, y where 
        ! the color value of the polygon points is given in z. N should 
        ! give the total number of points. Polygons should be separated
        ! with NaNs.

        ! The usual
        implicit none

        ! Declare variables
        integer                             :: n, i, fu
        real(R8), dimension(n)              :: x, y, z
        character(*)                        :: gnuplotoptions

        ! Initialize
        !===========
        ! Set the correct directories
        call SetGnuplotNames(plotfile,datafile,'patchplot')

        ! Write the data file
        !====================
        ! Write vertex coordinates to file
        open (action='write', file=trim(datafile), newunit=fu, &
             status='replace')
    
        do i = 1, n
            if (isnan(x(i))) then 
                ! Write blank line
                write(fu, *)
            else
                ! Write coordinates
                write (fu, *) x(i), y(i), z(i)
            end if
        end do
    
        close (fu)

        ! Call plotter
        !=============
        call gnuplotexe(gnuplotoptions,trim(plotfile))

    end subroutine

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

    ! Optimization
    !=============

    ! State
    !======
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
        integer                             :: i, j, k, fu, np
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
    !                         Auxiliary routines                       !
    !------------------------------------------------------------------!
    
    ! gnuplot execution
    subroutine gnuplotexe(gnuplotoptions,plotfile)

        ! Description
        !============
        ! Wrapper to execute gnuplot from command line. 'plotfile' 
        ! should refer to the .plt plot file, gnuplotoptions should be 
        ! a string containing more plotting options. 

        ! The usual
        implicit none

        ! Declare variables
        character(*)             :: gnuplotoptions
        character(*)             :: plotfile

        print *, trim(plotfile)
        ! Execute gnuplot
        call execute_command_line('gnuplot ' // trim(gnuplotoptions) &
            // ' ' // trim(plotfile))

    end subroutine

    ! Plot name and data file determination
    subroutine SetGnuplotNames(plt,dat,funname)

        ! Description
        !============
        ! Wrapper to set the plot name (.plt), data file (.dat) for the 
        ! specified function 'funname'.

        ! The usual
        implicit none

        ! Declare variables
        character(C128)         :: plt
        character(C128)         :: dat
        character(*)         :: funname

        ! Set output
        plt = trim(plotdir) // trim(funname) // '.plt'
        dat = trim(plotdir) // trim(funname) // '.dat'

    end subroutine

end module 