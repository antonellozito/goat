!   ______     _     __                                                  
!  / ____/____(_)___/ /                                                  
! / / __/ ___/ / __  /                                                   
!/ /_/ / /  / / /_/ /                                                    
!\____/_/  /_/\__,_/                                                     
!   ____        __  _           _             __  _                ___   
!  / __ \____  / /_(_)___ ___  (_)___  ____ _/ /_(_)___  ____     ( _ )  
! / / / / __ \/ __/ / __ `__ \/ /_  / / __ `/ __/ / __ \/ __ \   / __ \/|
!/ /_/ / /_/ / /_/ / / / / / / / / /_/ /_/ / /_/ / /_/ / / / /  / /_/  < 
!\____/ .___/\__/_/_/ /_/ /_/_/ /___/\__,_/\__/_/\____/_/ /_/   \____/\/ 
!    /_/       __            __        __  _                             
!   /   | ____/ /___ _____  / /_____ _/ /_(_)___  ____                   
!  / /| |/ __  / __ `/ __ \/ __/ __ `/ __/ / __ \/ __ \
! / ___ / /_/ / /_/ / /_/ / /_/ /_/ / /_/ / /_/ / / / /                  
!/_/  |_\__,_/\__,_/ .___/\__/\__,_/\__/_/\____/_/ /_/                   
!  ______         /_/____                                                
! /_  __/___  ____  / / /_  ____  _  __                                  
!  / / / __ \/ __ \/ / __ \/ __ \| |/_/                                  
! / / / /_/ / /_/ / / /_/ / /_/ />  <                                    
!/_/  \____/\____/_/_.___/\____/_/|_|  

!======================================================================!
!                                                                      !
!                            DOCUMENTATION                             !
!                                                                      !
!======================================================================!

! Short description
!==================
! This is the main driver to run the grid optimization and adaptation 
! toolbox (goat). To run a case, a .dat file should be provided with the 
! necessary user input (see also gdmod_userinput and the manual). 
! This program is used to call the desired grid driver to adapt and/or
! deform the grid. See the separate drivers for more info. 

program Goat

    ! Modules
    !========
    use goatmod_userinput
    use mod_global_environment, only: solps, SolpsPreamble
    use mod_plotter, only: plotdir

#if (defined(MUMPS) || defined(USE_MPI))
    use mpi
#endif


    ! Declare variables
    !==================
    ! Arguments

    ! Loop variables

    ! Auxiliary
    type(GoatoptionsUDT)        :: goatoptions
    character(:), allocatable   :: filepath

#if (defined(MUMPS) || defined(USE_MPI))
    ! MPI
    integer                     :: rank, nrank, ierror
#endif 

    ! Initialize
    !===========
#if (defined(MUMPS) || defined(USE_MPI))
    ! MPI call
    call mpi_init(ierror)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierror)
    call mpi_comm_size(MPI_COMM_WORLD, nrank, ierror)
    if (rank == 0) then 
        print *, 'mpi initialized with: ', nrank, 'processes'
        print *, MPI_COMM_WORLD
        if (nrank > 1) then 
            print *, 'warning: goat not yet capable of using multiple processes, mpi only intended for MUMPS solver'
        end if 
    end if 
#endif 

    ! Set the filepath
    allocate(character(len('./GOAToptions.dat')) :: filepath)
    filepath = './GOAToptions.dat'

    ! Call solps preamble
    if (solps) then 
        call SolpsPreamble('goat')
    else
        call execute_command_line('mkdir ' // plotdir)
    end if

    ! Read the user input
    !====================
    ! fileID should always be GOAToptions.dat
    goatoptions%inputfilepath = filepath 
    print *, 'Reading goat options from file: ', goatoptions%inputfilepath
    call goatoptions%Set()
    

    ! Run driver
    !===========
    select case (goatoptions%meth)

    case ('GD')

        ! Grid deformation only
        call GDDriver(goatoptions)

    case ('GG')

        ! Grid generation only 
        call GGDriver(goatoptions)

    case ('GGGD')

        ! Grid generation with subsequent grid deformation
        call GGGDDriver(goatoptions)

    case ('GDtest')

        call GDtestdriver(goatoptions)

    case default 

        ! Call error handler
        call gdErrorHandler('Goat: unknown driver option')

    end select 

    ! Finalize
    !=========
    ! Error handler
    call ErrorStack%Print()

#if (defined(MUMPS) || defined(USE_MPI))
    ! MPI call
    call mpi_finalize(ierror)
#endif 



end program Goat
                                    