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

    ! Declare variables
    !==================
    ! Arguments

    ! Loop variables

    ! Auxiliary
    type(GoatoptionsUDT)    :: goatoptions
    character(:), allocatable    :: filepath

    ! Initialize
    !===========
    allocate(character(len('./GOAToptions.dat')) :: filepath)
    filepath = './GOAToptions.dat'

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

    case ('GDtest')

        call GDtestdriver(goatoptions)

    case ('ShapeOpt')

        ! Move to different program?
        call ShapeOptDriver(goatoptions)

    case default 

        ! Call error handler
        call gdErrorHandler('Goat: unknown driver option')

    end select 

end program Goat
                                    