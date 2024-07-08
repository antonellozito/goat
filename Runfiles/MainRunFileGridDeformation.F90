
!    ______     _     __
!   / ____/____(_)___/ /
!  / / __/ ___/ / __  /
! / /_/ / /  / / /_/ /
! \____/_/  /_/\__,_/
!        __     ____                           __  _
!   ____/ /__  / __/___  _________ ___  ____ _/ /_(_)___  ____
!  / __  / _ \/ /_/ __ \/ ___/ __ `__ \/ __ `/ __/ / __ \/ __ \
! / /_/ /  __/ __/ /_/ / /  / / / / / / /_/ / /_/ / /_/ / / / /
! \__,_/\___/_/  \____/_/  /_/ /_/ /_/\__,_/\__/_/\____/_/ /_/
!    __            __  __                    __
!   / /____  _____/ /_/ /_  ___  ____  _____/ /_
!  / __/ _ \/ ___/ __/ __ \/ _ \/ __ \/ ___/ __ \
! / /_/  __(__  ) /_/ /_/ /  __/ / / / /__/ / / /
! \__/\___/____/\__/_.___/\___/_/ /_/\___/_/ /_/

!======================================================================!
!                                                                      !
!                            DOCUMENTATION                             !
!                                                                      !
!======================================================================!

! Short description
!==================
! This is the main driver to run the grid deformation routines, developed
! for plasma edge grids (and possibly others). It allows to redistribute
! the grid nodes in different ways (optimization based or others, if
! available), starting from a known initial grid. 

! Running the code
!=================
! Normally, everything in the run file is set such that a single case can
! be run. The regular work flow can be found in the driver routines, e.g.
! RunGridOptimization.m, and typically consists of the following steps:
!
! 1) Set the initial options for the grid type, grid deformation method,
! and possibly others (e.g. constraints on the grid nodes, such as
! alignment with the magnetic field).
! 2) Construct the initial grid, either from scratch or based on a saved
! grid (see also ConstructGrid.m)
! 3) Redistribute the nodes
! 4) Do post process if necessary
! 
! Below, a list of the most important functions is given, including a short
! description. 

program GridDeformation
    ! Initialization
    !===============
    ! Declare modules to be used
    use gdmod_types
    use gdmod_userinput
    use gdmod_optimizationengine
    
    ! The usual
    implicit none

    ! Declare variables
    type(RunfileOptionsUDT)    :: runfileOptions ! options for mainfile
    type(ExportOptionsUDT)     :: exportOptions ! options for exporting
    type(GridUDT)              :: grid ! grid structure
    type(OptimizationEngineGDUDT)      :: optimizationdriver 

    ! Set the main options
    !=====================
    ! Set the reading path
    runfileOptions%inputfilepath = './Examples/TCV/GOAToptions.dat'

    ! Set the default values
    call runfileOptions%Set()

    ! Main driver
    !============
    select case (runfileOptions%runtype)

    case ('optimize')

        ! Run the main optimization driver
        call RunGridOptimization(grid, optimizationdriver, runfileOptions)

    case default 

        ! Unknown case - print error
        print *, 'gd: unknown mainOptions%runtype value: ', runfileOptions%meth

    end select

    ! Post-processing
    !================

    ! Export
    !=======
    if (runfileOptions%export == 1) then
        
        ! Set the export options
        call SetExportOptions(exportOptions)

        ! Export
        call ExportGridData(grid,exportOptions)

    endif


end program GridDeformation