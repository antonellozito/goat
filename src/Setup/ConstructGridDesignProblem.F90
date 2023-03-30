subroutine ConstructGridDesignProblem(optimizationdriver, &
    designoptions, grid, magneticField, environment)

    ! Description
    !============
    ! Construct the initial grid design problem by initializing the 
    ! following quantities:
    ! 
    ! - Cost function: type, parameters
    ! - design variables: type, values
    ! - constraints: 

    ! Notes
    !======

    ! Initialize
    !===========
    ! Declare modules
    use gdmod_types 
    use gdmod_userinput
    use gdmod_designvariables
    use gdmod_optimizationengine
    ! use gdmod_costfunction

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    type(OptimizationEngineGDUDT)       :: optimizationdriver
    type(DesignOptionsUDT), intent(in)  :: designoptions
    type(GridUDT), intent(in)           :: grid
    type(MagneticFieldUDT), intent(in)  :: magneticField
    type(EnvironmentUDT), intent(in)    :: environment

    ! Loop variables

    ! Auxiliary variables 

    ! Data

    ! State
    !======
    optimizationdriver%inputfilepath = designoptions%inputfilepath
    call optimizationdriver%SetupOptimizationDriver()
    
    ! Associate in order to execute select type...
    associate(thisproblem => optimizationdriver%problem) 

        select type(thisproblem)

        type is (OptimizationProblemGDUDT)
    
            ! This should be the only possible type 
            thisproblem%grid            = grid 
            thisproblem%magneticField   = magneticField 
            thisproblem%environment     = environment
            thisproblem%designoptions   = designoptions

        class default

            ! Unknown, throw error
            stop 'Unknown optimization problem type'

        end select

    end associate

end subroutine
