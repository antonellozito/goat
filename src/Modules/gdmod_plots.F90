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

    ! Grid 
    !=====
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

    ! Optimization
    !=============


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