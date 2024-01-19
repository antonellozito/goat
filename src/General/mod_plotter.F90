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

module mod_plotter

    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN

    ! The usual
    implicit none
    save
    public 

    ! Define the (relative) file directory where the plot are located
    character(len=*), parameter     :: plotdir = './src/Visualization/'

    ! Define the (general) plot filename and data filename
    character(:), allocatable                 :: plotfile
    character(:), allocatable                 :: datafile

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

    ! 2D quiverplot
    subroutine Quiverplot2D(x, y, dx, dy, n, gnuplotoptions)

        ! Description
        !============
        ! Plot a generic 2D quiver plot with arrows at the locations
        ! (x, y) with x- and y-component sizes dx and dy. 

        ! The usual
        implicit none

        ! Declare variables
        integer                             :: n, i, fu
        real(R8), dimension(n)              :: x, y, dx, dy
        character(*)                        :: gnuplotoptions

        ! Initialize
        !===========
        ! Set the correct directories
        call SetGnuplotNames(plotfile,datafile,'quiverplot2d')

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
                write (fu, *) x(i), y(i), dx(i), dy(i)
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

    ! Spy plot
    subroutine SpyPlot(x, y, n, gnuplotoptions)

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
        integer(I8), dimension(n)              :: x, y
        character(*)                        :: gnuplotoptions

        ! Initialize
        !===========
        ! Set the correct directories
        call SetGnuplotNames(plotfile,datafile,'spyplot')

        ! Write the data file
        !====================
        ! Write vertex coordinates to file
        open (action='write', file=trim(datafile), newunit=fu, &
             status='replace')
    
        do i = 1, n
            ! Write coordinates
            write (fu, *) x(i), y(i)
        end do
    
        close (fu)

        ! Call plotter
        !=============
        call gnuplotexe(gnuplotoptions,trim(plotfile))

    end subroutine

    ! 2D structured state field
    subroutine Plot2DStructuredField(v, x, y, nx, ny, gnuplotoptions)

        ! Description
        !============
        ! Routine to plot 2D unstructured data, e.g. the magnetic field.
        ! x, y should be nx-by-1 and ny-by-1 arrays. v should be an 
        ! nx-by-ny array. Plotting is taken care of by the Patchplot2D 
        ! routine. 

        ! The usual
        implicit none

        ! Declare variables
        integer(I8), intent(in)             :: nx, ny
        integer(I8)                         :: i, j, k, np
        real(R8)                            :: v(1:nx,1:ny), & 
            x(1:nx), y(1:nx)
        character(*)                        :: gnuplotoptions
        real(R8), allocatable               :: xp(:), yp(:), val(:)

        ! Set NaN
        real(R8)                            :: NaN

        ! Initialize
        !===========
        ! Set NaN
        NaN = IEEE_VALUE(nan, IEEE_QUIET_NAN)

        ! Write the data 
        !===============
        ! Number of polygons to be plotted (one cell each)
        np = (nx - 1) * (ny  - 1)

        ! Allocate
        allocate(xp(5*np))
        allocate(yp(5*np))
        allocate(val(5*np))

        ! Loop over both coordinate directions
        k = 1
        do i = 1, ny-1
            do j = 1, nx-1
                ! Write
                xp(k) = x(j)
                yp(k) = y(i)
                val(k) = v(j, i)
                k = k+1

                xp(k) = x(j+1)
                yp(k) = y(i)
                val(k) = v(j+1, i)
                k = k+1

                xp(k) = x(j+1)
                yp(k) = y(i+1)
                val(k) = v(j+1, i+1)
                k = k+1

                xp(k) = x(j)
                yp(k) = y(i+1)
                val(k) = v(j, i+1)
                k = k+1

                xp(k) = NaN
                yp(k) = NaN
                val(k) = NaN
                k = k+1

            end do
        end do

        ! Call the plotter
        !=================
        call Patchplot2D(xp, yp, val, 5*np, gnuplotoptions)

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
        character(:), allocatable       :: plt
        character(:), allocatable       :: dat
        character(*)         :: funname

        ! Set output
        plt = trim(plotdir) // trim(funname) // '.plt'
        dat = trim(plotdir) // trim(funname) // '.dat'

    end subroutine

    !------------------------------------------------------------------!
    !                          Writing routines                        !
    !------------------------------------------------------------------!

    ! Vertex writing routine
    subroutine WriteVertexData(ID, x, y, filepath)

        ! Description
        !============
        ! General vertex data writing routine. writes IDs, x, y

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), allocatable, intent(in)        :: ID(:)
        real(R8), allocatable, intent(in)           :: x(:), y(:) 
        integer                                     :: i, fu
        character(:), allocatable, intent(in)       :: filepath 

        ! Initialize
        !===========
        ! Set the correct directories
        call SetGnuplotNames(plotfile, datafile, filepath)

        ! Write
        !======
        open (action='write', file=trim(datafile), newunit=fu, &
             status='replace')
            
        ! Header
        write (fu, *) 'ID, x, y ', 'nrow', size(ID, 1), 'ncol', 3
    
        ! Data
        do i = 1, size(ID)
            write (fu, *) ID(i), x(i), y(i)
        end do
    
        close (fu)

    end subroutine

    ! Vertex pair writing routine
    subroutine WriteVertexPairData(ID, x, y, filepath)

        ! Description
        !============
        ! General vertex data writing routine. writes IDs, x, y in the 
        ! format ID1, ID2, ... IDn, x1, x2, ... xn, y1, y2, ... yn

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), allocatable, intent(in)        :: ID(:, :)
        real(R8), allocatable, intent(in)           :: x(:, :), y(:, :) 
        integer                                     :: i, fu
        character(:), allocatable, intent(in)       :: filepath 

         ! Initialize
        !===========
        ! Set the correct directories
        call SetGnuplotNames(plotfile, datafile, filepath)

        ! Write
        !======
        open (action='write', file=trim(datafile), newunit=fu, &
             status='replace')
            
        ! Header
        write (fu, *) 'IDn, xn, yn ', 'nrow', size(ID, 1), 'ncol', 3*size(ID, 2)
    
        ! Data
        do i = 1, size(ID, 1)
            write (fu, *) ID(i, :), x(i, :), y(i, :)
        end do
    
        close (fu)

    end subroutine

    ! Polygon writing routine
    subroutine Write2DPolygonData(x, y, filepath)

        ! Description
        !============
        ! General polygon data writing routine. writes x, y coordinates.
        ! Data can contain NaNs - these are replaced by white lines, and
        ! should be used to distinguish different polygons. 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), allocatable, intent(in)           :: x(:), y(:) 
        integer                                     :: i, fu
        character(:), allocatable, intent(in)       :: filepath 


        
        ! Initialize
        !===========
        ! Set the correct directories
        call SetGnuplotNames(plotfile, datafile, filepath)

        ! Write
        !======
        open (action='write', file=trim(datafile), newunit=fu, &
             status='replace')
        print *, datafile
            
        ! Header
        write (fu, *) 'x, y'
    
        ! Data
        do i = 1, size(x)
            if (isnan(x(i))) then 
                write (fu, *)
            else 
                write (fu, *) x(i), y(i)
            end if
        end do
    
        close (fu)

    end subroutine

end module