subroutine RunGridOptimization(grid,designParams,options)
    ! Description
    !============
    ! Main driver for optimization-based grid deformation. This routine
    ! is a wrapper where, depending on the options parsed in the input, 
    ! the initialization is done and the actual optimization driver is
    ! called. 

    ! Initialize
    !===========
    ! Declare modules
    use gdmod_types
    use gdmod_userinput 

    ! The usual
    implicit none

    ! Declare variables
    type(GridUDT)              :: grid
    type(DesignParamsUDT)      :: designParams
    type(RunfileOptionsUDT)    :: options 

    ! Set additional options
    !=======================
    ! 

end subroutine