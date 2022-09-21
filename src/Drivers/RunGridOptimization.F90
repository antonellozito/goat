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
    type(GridUDT)                   :: grid
    type(MagneticFieldUDT)          :: magneticField
    type(DesignParamsUDT)           :: designParams
    type(RunfileOptionsUDT)         :: options 
    type(GridOptionsUDT)            :: gridoptions
    type(DesignOptionsUDT)          :: designoptions
    type(MagneticFieldOptionsUDT)   :: mfoptions
    type(NumOptionsUDT)             :: num

    ! Set additional options
    !=======================
    ! Set grid options
    call SetGridOptions(gridoptions)

    ! Set optimization options
    call SetDesignOptions(designoptions)

    ! Set other numerical parameters
    call SetNumOptions(num);
    
    ! Magnetic field 
    call SetMagneticFieldOptions(mfoptions);

    ! Initialize
    !===========
    ! Construct the initial grid
    call ConstructGrid(grid,gridoptions,options);

    ! Construct the initial magnetic field
    call ConstructMagneticField(mfoptions, magneticField)


end subroutine