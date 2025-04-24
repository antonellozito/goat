! 
!    _____ __                                                  
!   / ___// /_  ____ _____  ___                                
!   \__ \/ __ \/ __ `/ __ \/ _ \
!  ___/ / / / / /_/ / /_/ /  __/                               
! /____/_/ /_/\__,_/ .___/\___/                                
!    ____        _/_/_           _             __  _           
!   / __ \____  / /_(_)___ ___  (_)___  ____ _/ /_(_)___  ____ 
!  / / / / __ \/ __/ / __ `__ \/ /_  / / __ `/ __/ / __ \/ __ \
! / /_/ / /_/ / /_/ / / / / / / / / /_/ /_/ / /_/ / /_/ / / / /
! \____/ .___/\__/_/_/ /_/ /_/_/ /___/\__,_/\__/_/\____/_/ /_/ 
!     /_/                                                      

!======================================================================!
!                                                                      !
!                            DOCUMENTATION                             !
!                                                                      !
!======================================================================!

! Short description
!==================
! This is the main driver to run shape optimization with the goat 
! grid deformation as a constraint.

program ShapeOptimization

    ! Modules
    !========
    use somod_userinput
    use mod_global_environment, only: solps, SolpsPreamble
    use mod_plotter, only: plotdir


    implicit none

    ! Declare variables
    !==================
    ! Arguments

    ! Loop variables

    ! Auxiliary
    character(:), allocatable           :: filepath

    ! Initialize
    !===========
    allocate(character(len('./SOoptions.dat')) :: filepath)
    filepath = './SOoptions.dat'

    ! Call solps preamble
    if (solps) then 
        call SolpsPreamble('goat')
    else
        call execute_command_line('mkdir ' // plotdir)
    end if

    ! Read the user input
    !====================
    ! fileID should always be SOoptions.dat
    print *, 'Reading shape optimization options from file: ', filepath

    ! Run driver
    !===========
    call ShapeOptDriver(filepath)

    ! Print out the error stack
    !==========================
    call ErrorStack%Print()


end program ShapeOptimization
                                    