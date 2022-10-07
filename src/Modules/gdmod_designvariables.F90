!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains all the necessary routines to manipulate the 
! design variables (except for the actual optimization routines). The
! design variable UDT should be set in gdmod_types, user input in 
! gdmod_userinput. 

! Notes
!======
! Note 1: this module is specifically built for grid optimization. Most
! routines therefore expect at least the input arguments 'grid', 
! 'magneticField', and 'environment'. 

! Note 2: although the code structure (should) allow(s) for any type of 
! design variable, currently only the 'coordinates' type is implemented.

module gdmod_designvariables
    
    ! Initialize
    !============
    ! Load modules
    use gdmod_types
    use gdmod_userinput 
    use optmod_designvariables
    use gdmod_state

    ! The usual
    implicit none
    save
    public 

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    ! Derived types
    !==============
    ! Should all extend the base type DesignVariablesUDT, defined in 
    ! optmod_designvariables.

    ! Design variables for coordinates
    type, extends(DesignVariablesUDT) :: DesignVariablesCoordinatesUDT 

    contains 

        procedure ::  InitializeDesign => InitializeDesignCoordinates

    end type

    type, abstract :: DesignVariablesParUDT

    contains

        procedure(InitializeDesignVariablesParametersINT), deferred :: &
            InitializeDesignVariablesParameters 

    end type  


    ! Parameter type for coordinates
    type, extends(DesignVariablesParUDT) :: DesignVariablesCoordinatesParUDT

        ! Indices for x, y coordinates in phi
        integer(I8), allocatable        :: xind(:), yind(:)

    contains 

        ! Comment for select type behavior
        procedure :: InitializeDesignVariablesParameters => InitializeDesignParametersCoordinates

    end type

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    ! Parameters
    abstract interface
    
        ! Design initialization
        subroutine InitializeDesignVariablesParametersINT(parameters, &
            grid, magneticField, environment)
            import :: DesignVariablesParUDT, &
                gridUDT, MagneticFieldUDT, EnvironmentUDT 
            class(DesignVariablesParUDT)        :: parameters 
            type(GridUDT)                       :: grid
            type(MagneticFieldUDT)              :: magneticField
            type(EnvironmentUDT)                :: environment
        end subroutine

    end interface


    contains

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                              Coordinates                         !
    !------------------------------------------------------------------!

    ! Design initialization
    subroutine InitializeDesignCoordinates(designvariables, state)

        ! Description
        !============
        ! Initialize the design variables, which are in this case only 
        ! the x and y coordinates. 

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariablesCoordinatesUDT)           :: designvariables
        class(StateUDT)                                :: state 
    
        ! Loop variables

        ! Auxiliary variables 

        ! Data

    end subroutine

    ! Design initialization (main routine)
   ! subroutine SetupDesign(designvariables, grid, magneticField, &
    !    environment)

        ! Description
        !============
        ! Initialize the design variables. It is assumed the design 
        ! variable type is already set. Based on that, the variables phi
        ! are initialized differently. 

        ! Notes
        !======

        ! Initialize
        !===========
        ! Declare modules

        ! The usual

        ! Declare variables
        !==================
        ! Arguments
      !  type(DesignVariablesUDT)            :: designvariables 
      !  type(GridUDT), intent(in)           :: grid 
      !  type(MagneticFieldUDT), intent(in)  :: magneticField
      !  type(EnvironmentUDT), intent(in)    :: environment

        ! Loop variables

        ! Auxiliary variables 
      !  class(DesignVariablesCoordinatesParUDT), pointer  :: coordinatespar 
      !  class(DesignVariablesCoordinatesParUDT), allocatable :: mycoordinatespar

        ! Data

        ! Initialize design
        !==================
        ! Check the design type
      !  select case(trim(designvariables%type))

      !  case ('coordinates')

            ! Design variables are x, y coordinates (first x, then y)

            ! Compute total number of design variables
      !      designvariables%nphi = 2*grid%vert%ntot 

            ! Allocate
      !      allocate(designvariables%phi(designvariables%nphi))

            ! Set parameters
      !      allocate(DesignVariablesCoordinatesParUDT::designvariables%parameters)

            ! Comment for select type beahvior
      !      call designvariables%parameters%InitializeDesignParameters(grid, magneticField, environment)

            ! This should not work here
            ! call designvariables%parameters%mydummyroutineforclasses(grid)
       !     select type(associate=>designvariables%parameters)

      !      type is (DesignVariablesCoordinatesParUDT)

                ! Though we get here... 
       !         print *, 'we are here in the select type while running classes'

                ! ... None of this works - perfect 'encapsulation'... 
                ! print *, designvariables%parameters%xind(1)
                !parameters%mydummyroutineforclasses(grid)
                !coordinatespar => designvariables%parameters

                ! Though If we really want, we can use other objects NOT related to the original designvariables%parameters...
         !       allocate(DesignVariablesCoordinatesParUDT::mycoordinatespar)

         !       call mycoordinatespar%InitializeDesignParameters(grid, magneticField, environment)
         !       call mycoordinatespar%mydummyroutineforclasses(grid)
!
                ! Note that this should not work
                ! call mycoordinatespar%InitializeDesignParametersCoordinates(grid, magneticField, environment)


          !  end select

            ! Uncomment for select type behavior
            ! call InitializeDesignParametersCoordinates(designvariables%parameters, grid, magneticField, environment)

            ! Initialize design variables
            ! call SetDesignVariables(designvariables%parameters, grid, magneticField, environment)


       ! case default

            ! Throw error
       !     call gdErrorHandler('Design variable type not implemented')

       ! end select


   !end subroutine

    ! Design variable value initialization
    !=====================================
    subroutine SetDesignVariablesCoordinates(parameters, &
        designvariables, grid, magneticField, environment)

        ! Description
        !============
        ! Set the design variable parameters for the coordinates.

        ! Notes
        !======

        ! Initialize
        !===========
        ! Declare modules

        ! The usual

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariablesUDT)                :: designvariables
        class(DesignVariablesCoordinatesParUDT)  :: parameters 
        type(GridUDT)              :: grid 
        type(MagneticFieldUDT)     :: magneticField
        type(EnvironmentUDT)       :: environment

        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Initialize values
        !==================
        ! Simply x and y coordinates here
        !designvariables%phi(parameters%xind) = grid%vert%x
        !designvariables%phi(parameters%yind) = grid%vert%y


    end subroutine

    ! Parameter initialization
    !=========================
    subroutine InitializeDesignParametersCoordinates(parameters, &
        grid, magneticField, environment)

        ! Description
        !============
        ! Set the design variable parameters for the coordinates.

        ! Notes
        !======

        ! Initialize
        !===========
        ! Declare modules

        ! The usual

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariablesCoordinatesParUDT)  :: parameters 
        type(GridUDT)              :: grid 
        type(MagneticFieldUDT)     :: magneticField
        type(EnvironmentUDT)       :: environment

        ! Loop variables

        ! Auxiliary variables 

        ! Data        

    end subroutine

end module