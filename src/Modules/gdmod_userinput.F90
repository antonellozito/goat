!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains all the user-defined input parameters and related
! routines for the grid deformation module. 

module gdmod_userinput

    ! Initialize
    !============
    ! Load modules
    use gdmod_types

    ! The usual
    implicit none
    save
    private 

    ! Read namelist from file
    ! namelist /gridsmoothing/ 

    ! All functions
    public SetRunfileOptions
    public SetExportOptions
    public SetDesignOptions
    public SetGridOptions
    public SetNumOptions
    public SetMagneticFieldOptions
    public SetEnvironmentOptions 
    public SetVesselOptions

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!
    
    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    contains 

    subroutine SetRunfileOptions(options)
        
        ! Description
        !============
        ! Set the run options. The following fields have to be set:
        !
        ! - runtype: define how to run the grid deformation. 'optimize'
        ! will use the optimization drivers (which can be further 
        ! specified using 'meth') - other options are not yet 
        ! implemented. 
        !
        ! - gridtype: type of grid to be considered. Can only be 
        ! 'plasma' for now, though in the future e.g. EIRENE grids could
        ! also be considered.
        !
        ! - meth: method to be used for the specific runtype. For 
        ! 'optimize', only 'KKT' is currently available. Here, the full
        ! KKT conditions of the grid deformation optimization problem 
        ! are set up and solved.
        !
        ! - export: export the results in a certain format to a 
        ! specified file (can be true or false). See also the
        ! SetExportOptions subroutine in this file. 

        ! Declaration
        type (RunfileOptionsUDT), intent(inout)    :: options

        ! Default options
        options%runtype     = 'optimize'
        options%gridtype    = 'plasma'
        options%meth        = 'KKT'
        options%export      = .true.
        
        ! Overwrite user-defined options
        ! options%runtype     = runtype

    end subroutine

    subroutine SetExportOptions(options)
        
        ! Description
        !============
        ! Set the export options. The following fields have to be set:
        !
        ! - gridformat: the format in which the grid needs to be written
        ! either 'structured' or 'unstructured'.

        ! Declaration
        type (ExportOptionsUDT), intent(inout)    :: options

        ! Default options
        options%gridformat     = 'structured'

    end subroutine

    subroutine SetDesignOptions(options)

        ! Description
        !============
        ! Set the options for the design. The following fields have to
        ! be set:
        !
        ! - costfunction%type: the type of cost function to be 
        ! considered. See EvaluateCostFunction.F for the different 
        ! options. A good default option is the 'LR_FAD' cost function,
        ! which combines the length ratio along the magnetic field lines
        ! with a contribution of the difference in angle between 
        ! magnetic field vector and (radial) face vector.
        ! 
        ! - variables%type: the type of design variables to be 
        ! considered. Currently, this is (in general) all the 2D
        ! plasma edge grid coordinates (x,y). Note that all constraints 
        ! and cost function assume that at least the coordinates are 
        ! free to be optimized. 
        !
        ! - constraints%fluxfunction: -logical which decides whether
        ! the grid nodes should be aligned to the magnetic field flux 
        ! surfaces.
        !
        ! - constraints%xpoints: logical which decides whether the 
        ! x-point location(s) should be fixed. This assumes that the
        ! x-point location(s) are known beforehand.
        !
        ! - constraints%orthogonality: logical which decides whether the
        ! radial surfaces should remain orthogonal to the field. Only 
        ! possible for a limited set of faces, usually in the core
        ! TO BE SET IN DIFFERENT INPUT FILE
        !
        ! - constraints%linefolding: prevent 'line folding', i.e. that 
        ! grid nodes move along the magnetic field lines over their 
        ! neighbours. Only possible if the flux function constraint is 
        ! set to true. 
        !

        ! Declaration
        type (DesignOptionsUDT), intent(inout)    :: options

        ! Default options, design variables
        options%variables%type                  = 'coordinates'

        ! Default options, cost function
        options%costfunction%type               = 'lengthratio'

        ! Default options, equality constraints
        options%constraints%boundaryfunctions   = .true.
        options%constraints%fluxfunction        = .true.
        options%constraints%xpoints             = .true.
        options%constraints%edgelengths         = .true.
        options%constraints%orthogonality       = .true.
        options%constraints%neq                 = 5 

        ! Default options, inequality constraints
        options%constraints%linefolding         = .true.
        options%constraints%nineq               = 1

    end subroutine

    subroutine SetGridOptions(options)
        ! Description
        !============
        ! Set the options for the grid. The following fields have to be
        ! set: 
        !
        ! - inputtype: string containing the type of input 
        !
        ! - vesselrefine: logical to decide whether the original vessel
        ! polygon should be further refined (increased resolution)
        !
        ! - vesselmaxdist: real number that gives the maximal distance 
        ! between vessel nodes to be considered when increasing the 
        ! resolution of the vessel polygon. Only active when vessel
        ! refinement is active. 

        ! Declaration
        type(GridOptionsUDT)        ::  options

        ! Default options
        options%inputtype           = 'b2fgmtry'
    
    end subroutine

    subroutine SetNumOptions(options)
        ! Description
        !============
        ! Set the numerical parameters and values. The following fields 
        ! have to be set:
        !

        ! Declaration
        type(NumOptionsUDT)         :: options

        ! Default options
        options%itmax = 1000

    end subroutine

    subroutine SetMagneticFieldOptions(options)
        ! Description
        !============
        ! Set the options for the magnetic field. The following fields
        ! have to be set:

        ! Declaration
        type(MagneticFieldOptionsUDT)       :: options

        ! Default options   
        options%readmeth                = 'readrzpsi'

    end subroutine

    subroutine SetVesselOptions(options)

        ! Description
        !============
        ! Set the options for the vessel.

        ! Declare variables
        !==================
        type(VesselOptionsUDT)          :: options

        ! Set options
        !============
        options%readmeth    = 'read_structure'
        options%geom        = 'ASDEX_Nathan'
        options%dir         = 'inputfiles/structure_tcv.dat'
        options%maxdist     = 0.01
        options%refine      = .true.

    end subroutine

    subroutine SetEnvironmentOptions(options)

        ! Description
        !============
        ! Set the options to read in the environment.

        ! Declare variables
        !==================
        type(EnvironmentOptionsUDT)         :: options

        ! Set options
        !============
        options%type = 'vessel'

    end subroutine


end module