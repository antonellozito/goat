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

    ! Design variables
    !=================
    ! Setup the design variables
    print *, 'we are here'
    call optimizationdriver%SetupOptimizationDriver()
    
    ! Set the type
    !optimizationproblem%designvariables%type = &
    !    designoptions%variables%type 

    ! Initialize the design
    !call InitializeDesign(optimizationproblem%designvariables, grid, magneticField, &
    !    environment)
    !call InitializeDesignParameters(optimizationproblem%designvariables%parameters, &
    !    optimizationproblem%designvariables, grid, magneticField, environment)

    ! Initialize design variables
    !call SetDesignVariables(optimizationproblem%designvariables%parameters, grid, &
    !       magneticField, environment)



end subroutine
