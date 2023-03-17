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
    use goatmod_userinput
    use mod_inputfileparser

    ! The usual
    implicit none
    save
    public 

    ! Read namelist from file
    ! namelist /gridsmoothing/ 

    ! All functions

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    ! Options for the main runfile type - extend from general options 
    type, extends(OptionsUDT) :: RunfileOptionsUDT

        ! Fields
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

        ! - ...filepath: filepaths for the different other options. In
        ! general, different files for the different options are 
        ! allowed, but default assumes the same filepath as for the 
        ! runoptions. 

        character(:), allocatable   :: runtype ! type of run: 'optimize' or 'test'
        character(:), allocatable   :: gridtype ! type of grid: 'plasma' 
        character(:), allocatable   :: meth ! method for grid deformation: 'KKT'
        integer                     :: export ! do export? 

        character(:), allocatable   :: gridfilepath, &
            environmentfilepath, designfilepath, magneticfieldfilepath

    contains 

        procedure :: SetDefaults    => SetDefaultRunfileOptions
        procedure :: Read           => ReadRunfileOptions

    end type

    ! Options for exporting
    type ExportOptionsUDT
        character(C32)       :: gridformat ! format of grid to export
    end type

    ! Options for design variables
    type, extends(OptionsUDT) :: DesignVariableOptionsUDT

        ! Fields
        ! - type:   type of design variables. Currently, only 
        !           'coordinates' is available
        character(:), allocatable    :: type 

    contains 

        procedure :: SetDefaults    => SetDefaultDesignVariableOptions
        procedure :: Read           => ReadDesignVariableOptions

    end type

    ! Options for the cost function
    type, extends(OptionsUDT) :: CostFunctionOptionsUDT

        ! Fields
        ! - type:   cost function type (e.g. LR_FAD, LR, LR2, ...) See 
        !           gdmod_costfunction to see which cost functions are 
        !           available

        ! Notes: parameters for the cost function are type specific and 
        ! are therefore not defined here. They should be read in by 
        ! dedicated routines defined in gdmod_costfunction
        character(:), allocatable   :: type 

    contains 

        procedure :: SetDefaults    => SetDefaultCostFunctionOptions
        procedure :: Read           => ReadCostFunctionOptions

    end type

    ! Options for the constraints
    type, extends(OptionsUDT) :: ConstraintOptionsUDT

        ! Fields
        ! - Per type of constraint, a field to indicate to see if the
        !   constraint should be imposed (then should be equal to 1)
        ! - total number of equality and inequality constraints

        ! Fields for equality constraints
        integer(I8)         :: boundaryfunctions ! impose boundary functions
        integer(I8)         :: fluxfunction ! impose constraints on flux
        integer(I8)         :: xpoints ! impose x-point location
        integer(I8)         :: edgelengths ! impose edge length cons
        integer(I8)         :: orthogonality 

        ! Fields for inequality constraints
        integer             :: linefolding ! prevent flux line folding

        ! Number of (continuous) constraints
        integer(I8)         :: neq ! number of equality constraints
        integer(I8)         :: nineq ! number of inequality constraints 

    contains 

        procedure :: SetDefaults    => SetDefaultConstraintOptions
        procedure :: Read           => ReadConstraintOptions

    end type

    ! Options for design optimization
    type, extends(OptionsUDT) :: DesignOptionsUDT

        ! Contains the different design option structures
        type(CostFunctionOptionsUDT)       :: costfunction
        type(DesignVariableOptionsUDT)     :: variables
        type(ConstraintOptionsUDT)         :: constraints

    contains 

        procedure :: SetDefaults    => SetDefaultDesignOptions
        procedure :: Read           => ReadDesignOptions 

    end type

    !==================================================================!
    !                                                                  !
    !                              INTERFACES                          !
    !                                                                  !
    !==================================================================!
    
    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    contains 

    !------------------------------------------------------------------!
    !                            Option setters                        !
    !------------------------------------------------------------------!

    subroutine SetDefaultRunfileOptions(options)
        
        ! Description
        !============
        ! Set the run options. 

        ! Declare variables
        !==================
        ! Arguments
        class(RunfileOptionsUDT)    :: options

        ! Default options
        !================
        options%runtype     = 'optimize'
        options%gridtype    = 'plasma'
        options%meth        = 'KKT'
        options%export      = 1

        options%gridfilepath            = options%inputfilepath 
        options%magneticfieldfilepath   = options%inputfilepath 
        options%designfilepath          = options%inputfilepath 
        options%environmentfilepath     = options%inputfilepath 

    end subroutine

    subroutine SetDefaultDesignVariableOptions(options)

        ! Description
        !============
        ! Set default design variable options

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariableOptionsUDT) :: options 

        ! Default options
        !================
        options%type        = 'coordinates' 

    end subroutine

    subroutine SetDefaultCostFunctionOptions(options)

        ! Description
        !============
        ! Set default cost function options

        ! Declare variables
        !==================
        ! Arguments
        class(CostFunctionOptionsUDT) :: options 

        ! Default options
        !================
        options%type        = 'LR_FAD' 

    end subroutine

    subroutine SetDefaultConstraintOptions(options)

        ! Description
        !============
        ! Set default constraint options

        ! Declare variables
        !==================
        ! Arguments
        class(ConstraintOptionsUDT) :: options 

        ! Default options
        !================
        options%boundaryfunctions   = 1
        options%fluxfunction        = 1
        options%xpoints             = 1
        options%edgelengths         = 0
        options%orthogonality       = 1

        options%linefolding         = 0

        options%neq                 = 4
        options%nineq               = 0

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

    subroutine SetDefaultDesignOptions(options)

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

        ! Declare variables
        !==================
        class(DesignOptionsUDT)                :: options

        ! Set defaults
        !=============
        ! Set inputfile paths 
        options%variables%inputfilepath = options%inputfilepath 
        options%costfunction%inputfilepath = options%inputfilepath 
        options%constraints%inputfilepath = options%inputfilepath

        ! Default options
        call options%variables%Set()
        call options%costfunction%Set() 
        call options%constraints%Set()

    end subroutine

    !------------------------------------------------------------------!
    !                            Option readers                        !
    !------------------------------------------------------------------!

    subroutine ReadRunfileOptions(options)

        ! Description
        !============
        ! Read in user-specified runfile options

        ! Declare variables
        !==================
        ! Arguments
        class(RunfileOptionsUDT)  :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: thisline, field
        integer, parameter              :: fid = 10 
        logical                         :: reachedeof

        ! Initialize
        !===========
        ! Variables
        reachedeof = .false. 

        ! Open the file, check if it exists
        open(unit=fid, file=options%inputfilepath, status='old', &
            iostat=openstatus)

        if (openstatus > 0) then 
            ! Something wrong when reading file - continue with default
            ! values
            print *, 'ReadRunfileOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadRunfileOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Run parameters
        field = 'gd.main.runtype'
        call ExtractOptionValueCharacter(fid, field, options%runtype)
        field = 'gd.main.gridtype'
        call ExtractOptionValueCharacter(fid, field, options%gridtype)
        field = 'gd.main.meth'
        call ExtractOptionValueCharacter(fid, field, options%meth)

        ! File paths
        field = 'gd.main.gridoptionsfile'
        call ExtractOptionValueCharacter(fid, field, options%gridfilepath)
        field = 'gd.main.magneticfieldoptionsfile'
        call ExtractOptionValueCharacter(fid, field, options%magneticfieldfilepath)
        field = 'gd.main.designoptionsfile'
        call ExtractOptionValueCharacter(fid, field, options%designfilepath)
        field = 'gd.main.environmentoptionsfile'
        call ExtractOptionValueCharacter(fid, field, options%environmentfilepath)

        ! Export
        field = 'gd.main.export'
        call ExtractOptionValueInteger0D(fid, field, options%export)

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)

    end subroutine

    subroutine ReadDesignVariableOptions(options)

        ! Description
        !============
        ! Read in user-specified design variable options

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariableOptionsUDT) :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: thisline, field
        integer, parameter              :: fid = 10 
        logical                         :: reachedeof

        ! Initialize
        !===========
        ! Variables
        reachedeof = .false. 

        ! Open the file, check if it exists
        open(unit=fid, file=options%inputfilepath, status='old', &
            iostat=openstatus)

        if (openstatus > 0) then 
            ! Something wrong when reading file - continue with default
            ! values
            print *, 'ReadDesignVariableOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadDesignVariableOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Type
        field = 'gd.design.dv.type'
        call ExtractOptionValueCharacter(fid, field, options%type)
        
        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)

    end subroutine

    subroutine ReadCostFunctionOptions(options)

        ! Description
        !============
        ! Read in user-specified cost function options

        ! Declare variables
        !==================
        ! Arguments
        class(CostFunctionOptionsUDT)   :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: thisline, field
        integer, parameter              :: fid = 10 
        logical                         :: reachedeof

        ! Initialize
        !===========
        ! Variables
        reachedeof = .false. 

        ! Open the file, check if it exists
        open(unit=fid, file=options%inputfilepath, status='old', &
            iostat=openstatus)

        if (openstatus > 0) then 
            ! Something wrong when reading file - continue with default
            ! values
            print *, 'ReadCostFunctionOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadCostFunctionOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Type
        field = 'gd.design.cfv.type'
        call ExtractOptionValueCharacter(fid, field, options%type)
        
        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)

    end subroutine

    subroutine ReadConstraintOptions(options)

        ! Description
        !============
        ! Read in user-specified constraints options

        ! Declare variables
        !==================
        ! Arguments
        class(ConstraintOptionsUDT)     :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: thisline, field
        integer, parameter              :: fid = 10 
        logical                         :: reachedeof

        ! Initialize
        !===========
        ! Variables
        reachedeof = .false. 

        ! Open the file, check if it exists
        open(unit=fid, file=options%inputfilepath, status='old', &
            iostat=openstatus)

        if (openstatus > 0) then 
            ! Something wrong when reading file - continue with default
            ! values
            print *, 'ReadCostFunctionOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadCostFunctionOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Equality constraints
        field = 'gd.design.ec.boundaryfunctions'
        call ExtractOptionValueInteger0D(fid, field, options%boundaryfunctions)
        field = 'gd.design.ec.fluxfunction'
        call ExtractOptionValueInteger0D(fid, field, options%fluxfunction)
        field = 'gd.design.ec.xpoints'
        call ExtractOptionValueInteger0D(fid, field, options%xpoints)
        field = 'gd.design.ec.edgelengths'
        call ExtractOptionValueInteger0D(fid, field, options%edgelengths)
        field = 'gd.design.ec.orthogonality'
        call ExtractOptionValueInteger0D(fid, field, options%orthogonality)

        ! Inequality constraints
        field = 'gd.design.inec.linefolding'
        call ExtractOptionValueInteger0D(fid, field, options%linefolding)

        ! Set number to the correct value
        options%neq     = 0
        options%nineq   = 0
        if (options%boundaryfunctions == 1) then 
            options%neq = options%neq + 1
        end if
        if (options%fluxfunction == 1) then 
            options%neq = options%neq + 1
        end if
        if (options%xpoints == 1) then 
            options%neq = options%neq + 1
        end if
        if (options%edgelengths == 1) then 
            options%neq = options%neq + 1
        end if
        if (options%orthogonality == 1) then 
            options%neq = options%neq + 1
        end if

        if (options%linefolding == 1) then 
            options%nineq = options%nineq + 1
        end if
        
        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)

    end subroutine

    subroutine ReadDesignOptions(options)

        ! Description
        !============
        ! Read in user-specified design options - currently nothing to 
        ! read in, since this is a superstructure without any specific 
        ! fields. 

        ! Declare variables
        !==================
        ! Arguments
        class(DesignOptionsUDT)         :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: thisline, field
        integer, parameter              :: fid = 10 
        logical                         :: reachedeof

        ! Initialize
        !===========
        ! Variables
        reachedeof = .false. 

        ! Open the file, check if it exists
        open(unit=fid, file=options%inputfilepath, status='old', &
            iostat=openstatus)

        if (openstatus > 0) then 
            ! Something wrong when reading file - continue with default
            ! values
            print *, 'ReadCostFunctionOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadCostFunctionOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Do nothing for now...
        
        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)

    end subroutine

end module