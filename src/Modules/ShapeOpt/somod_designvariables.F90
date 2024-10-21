!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains all the necessary routines to manipulate the 
! design variables (except for the actual optimization routines). There 
! are currently two types of design variables supported: vesselcoordinates
! and vesselcoordinates_goat. The first only considers the vessel 
! coordinates as true design variables. This can only be useful in a 
! reduced setting where the goat optimization problem is solved in each
! optimization iteration (this effect should be included in the cost 
! function actually). The latter considers also the goat optimization 
! variables (design variables, but also lagrange multipliers) and should
! only be used in a fully coupled monolithic solver. 

module somod_designvariables
    
    ! Initialize
    !============
    ! Load modules
    use optmod_designvariables
    use gdmod_optimizationengine

    ! The usual
    implicit none
    save
    public 

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    ! Abstract types
    !===============
    ! Main design variable type for shape optimization. Others should 
    ! be derived from this type. 
    type, abstract, extends(DesignVariablesUDT) :: DesignVariablesSOUDT

        character(:), allocatable                :: type 

    contains 

        ! Design initialization
        procedure(InitializeSOINT), deferred      :: Initialize

        ! Design update
        procedure(UpdateDesignSOINT), deferred      :: UpdateDesign

    end type

    ! Derived types
    !==============
    ! Should all extend the base type DesignVariablesSOUDT

    ! Design variables for (all) vessel coordinates
    type, extends(DesignVariablesSOUDT) :: DesignVariablesVesselCoordinatesUDT 

        ! Coordinate indices
        integer(I8), allocatable        :: xind(:), yind(:) 

    contains 
    
        ! Design initialization
        procedure :: Initialize     => InitializeDesignVesselCoordinates

        ! Design update
        procedure :: UpdateDesign   => UpdateDesignVesselCoordinates

    end type

    ! Design variables for vessel coordinates and goat design variables
    type, extends(DesignVariablesSOUDT) :: DesignVariablesVesselCoordinatesGoatUDT

        ! Coordinate & goat indices
        integer(I8), allocatable        :: xind(:), yind(:), &
            phigoatind(:), lambdagoatind(:), mugoatind(:)

    contains 

        ! Design initialization
        procedure :: Initialize     => InitializeDesignVesselCoordinatesGoat

        ! Design update
        procedure :: UpdateDesign   => UpdateDesignVesselCoordinatesGoat

    end type

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    ! Abstract interfaces
    !====================
    abstract interface

        ! Design initialization, specific for our purposes
        subroutine InitializeSOINT(designvariables, goat)

            ! Description
            !============
            ! Interface for the design initialization routine that each
            ! derived type should have. Since this is specific for the
            ! shape optimization routines, the goat has to be passed. 
            ! En avant avec la chèvre, copain

            ! Import
            import :: DesignVariablesSOUDT, OptimizationProblemGDUDT

            ! Declare
            class(DesignVariablesSOUDT)                 :: designvariables 
            type(OptimizationProblemGDUDT), intent(in)  :: goat

        end subroutine

        ! Design update
        subroutine UpdateDesignSOINT(designvariables, goat)

            ! Description
            !============
            ! Interface for the design update routine that each
            ! derived type should have. Since this is specific for the
            ! shape optimization routines, we have to update the goat
            ! here. 

            ! Import
            import :: DesignVariablesSOUDT, OptimizationProblemGDUDT

            ! Declare
            class(DesignVariablesSOUDT)                     :: designvariables 
            type(OptimizationProblemGDUDT), intent(inout)   :: goat

        end subroutine

    end interface

    contains

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                          Vessel coordinates                      !
    !------------------------------------------------------------------!

    ! Design initialization
    subroutine InitializeDesignVesselCoordinates(designvariables, goat)

        ! Description
        !============
        ! Initialize the design variables, which are in this case only 
        ! the x and y coordinates. 

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariablesVesselCoordinatesUDT)  :: designvariables
        type(OptimizationProblemGDUDT), intent(in)  :: goat

    
        ! Loop variables
        integer(I8)                                 :: i

        ! Auxiliary variables 
        integer(I8)                 :: nv 
        real(R8), allocatable       :: xv(:), yv(:)

        ! Data

        ! Initialize
        !===========
        ! Set the type
        designvariables%type = 'vesselcoordinates' 

        ! Get the vessel coordinates
        call goat%environment%vessel%GetVesselCoordinates(xv, yv)
        nv = size(xv)
        
        ! Set the number of design variables
        designvariables%nphi = 2*nv

        ! Set the design variables
        designvariables%phi = [xv, yv]

        ! Set other fields
        designvariables%xind = (/(i, i = 1, nv)/)
        designvariables%yind = designvariables%xind + nv

    end subroutine

    ! Design update
    subroutine UpdateDesignVesselCoordinates(designvariables, goat)

        ! Description
        !============
        ! Update the design coordinates according to the phi values. 

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariablesVesselCoordinatesUDT)      :: designvariables
        type(OptimizationProblemGDUDT), intent(inout)   :: goat
    
        ! Auxiliary variables 

        ! Data

        ! Update
        !=======
        ! Call updater (x and y already in correct order)
        call goat%UpdateProblemParameters(designvariables%phi, 'vesselcoordinates')

    end subroutine

    !------------------------------------------------------------------!
    !                       Vessel coordinates & Goat                  !
    !------------------------------------------------------------------!

    ! Design initialization
    subroutine InitializeDesignVesselCoordinatesGoat(designvariables, goat)

        ! Description
        !============
        ! Initialize the design variables, which are in this case  
        ! the x and y coordinates of the vessel and all variables of the
        ! goat optimization problem. 

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariablesVesselCoordinatesGoatUDT)  :: designvariables
        type(OptimizationProblemGDUDT), intent(in)  :: goat

    
        ! Loop variables
        integer(I8)                                 :: i

        ! Auxiliary variables 
        integer(I8)                         :: nv 
        real(R8), allocatable, dimension(:) :: xv, yv, phigoat, &
            lambdagoat, mugoat

        ! Data

        ! Initialize
        !===========
        ! Set the type
        designvariables%type = 'vesselcoordinates_goat' 

        ! Get the vessel coordinates
        call goat%environment%vessel%GetVesselCoordinates(xv, yv)
        nv = size(xv)

        ! Get the goat variables
        call goat%GetProblemDesignVariables(phigoat)
        lambdagoat  = goat%lambda 
        mugoat      = goat%mu 
        
        ! Set the number of design variables
        designvariables%nphi = 2*nv + size(phigoat) + size(lambdagoat) & 
            + size(mugoat)

        ! Set the design variables
        designvariables%phi = [xv, yv, phigoat, lambdagoat, mugoat]

        ! Set other fields
        designvariables%xind = (/(i, i = 1, nv)/)
        designvariables%yind = designvariables%xind + nv
        designvariables%phigoatind = [(i, i = 2*nv+1, size(phigoat)+2*nv)]
        designvariables%lambdagoatind = [(i, i = 1, size(lambdagoat))] + 2*nv+size(phigoat)
        designvariables%mugoatind = [(i, i = 1, size(mugoat))] + size(lambdagoat) + size(phigoat) + 2*nv

    end subroutine

    ! Design update
    subroutine UpdateDesignVesselCoordinatesGoat(designvariables, goat)

        ! Description
        !============
        ! Update the design coordinates according to the phi values. 

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariablesVesselCoordinatesGoatUDT)  :: designvariables
        type(OptimizationProblemGDUDT), intent(inout)   :: goat
    
        ! Auxiliary variables 

        ! Data

        ! Update
        !=======
        ! Call updater for vessel coordinates
        call goat%UpdateProblemParameters(&
            designvariables%phi([designvariables%xind, designvariables%yind]), &
            'vesselcoordinates')

        ! Update goat design variables
        goat%designvariables%phi = designvariables%phi(designvariables%phigoatind)

        ! Update goat lagrange multipliers
        goat%lambda = designvariables%phi(designvariables%lambdagoatind)
        goat%mu     = designvariables%phi(designvariables%mugoatind)
        
    end subroutine

end module