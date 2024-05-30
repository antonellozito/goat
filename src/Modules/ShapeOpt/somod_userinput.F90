!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains all options and structures to process user input
! for shape optimization.

module somod_userinput

    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_readwrite
    use mod_inputfileparser

    ! The usual
    implicit none
    save
    public 

    
    !------------------------------------------------------------------!
    !                         DESIGN VARIABLES                         !
    !------------------------------------------------------------------!

    ! Vessel coordinates design variable options
    type, extends(optionsUDT) :: DesignVariableOptionsVCUDT

        ! Desciption
        !===========
        ! Design variable options for the vessel coordinates. Here, not
        ! too many options are available, since we assume all vessel 
        ! coordinates are design variables (they may be constrained to 
        ! reduce the design space)
        
        ! No actual options needed ...
         
    contains

        procedure :: SetDefaults    => SetDefaultDesignVariableOptionsVC
        procedure :: Read           => ReadDesignVariableOptionsVC

    end type  

    ! Design variable options
    type, extends(optionsUDT) :: DesignVariableOptionsSOUDT

        ! Description
        !============
        ! General option structure for all design variable options. 
        ! Has different substructures for different other design 
        ! variable options

        character(:), allocatable   :: type ! design variable type
        type(DesignVariableOptionsVCUDT)    :: vc

    contains 

        procedure :: SetDefaults    => SetDefaultDesignVariableOptionsSO
        procedure :: Read           => ReadDesignVariableOptionsSO

    end type

    !------------------------------------------------------------------!
    !                           CONSTRAINTS                            !
    !------------------------------------------------------------------!

    ! Fixed vessel points constraints
    type, extends(optionsUDT)   :: FixedVesselPointsConOptionsUDT

        ! Description
        !============
        ! Constraint options for fixed vessel points. The following
        ! fields are present:
        !   structureIDs:   structure IDs of structures of which the 
        !                   vertices should be constrained 
        !   vertIDs:        IDs of specific vessel vertices (according
        !                   to the original numbering in the 
        !                   structure.dat file) to be constrained
        ! 
        ! Note that vertices are constrained to their original location

        ! Fields
        integer(I8), allocatable    :: structureIDs(:), vertIDs(:)

    contains 

        procedure :: SetDefaults    => SetDefaultFixedVesselPointsConOptions
        procedure :: Read           => ReadFixedVesselPointsConOptions

    end type

    ! General constraints
    type, extends(optionsUDT)   :: ConstraintOptionsSOUDT 

        ! Description
        !============
        ! General constraint options structure with different subfields
        ! for constraint specific options. Contains switches
        ! on/off (in)equality constraints. Note: for the goat constraints,
        ! no separate options structure is given since these should be
        ! read in separately anyway. 
        
        ! Switches
        integer(I8)                 :: fixedvesselpoints, goat

        ! Constraint options
        type(FixedVesselPointsConOptionsUDT)    :: fvpoptions

    contains 

        procedure :: SetDefaults    => SetDefaultConstraintOptionsSO
        procedure :: Read           => ReadConstraintOptionsSO

    end type

    !------------------------------------------------------------------!
    !                          COST FUNCTION                           !
    !------------------------------------------------------------------!

    ! Polygon levelset function cost function
    type, extends(optionsUDT)   :: CostfunctionOptionsPLFUDT

        ! Description
        !============
        ! Cost function options for the levelset function cost function.
        ! This includes:
        !
        !   vesselinputfilepath:    path to the inputfile for vessel 
        !                           options. This is likely the same as
        !                           the options file that defines the 
        !                           original vessel options
        !   newvesselfilepath:      filepath to the new vessel structure
        !                           (i.e. the new structure.dat file)
        !   lambda:                 scaling constant for cost function
        
        ! Fields
        real(R8)                    :: lambda
        character(:), allocatable   :: vesselinputfilepath, &
            newvesselfilepath

    contains 

        procedure :: SetDefaults    => SetDefaultCostfunctionOptionsPLF
        procedure :: Read           => ReadCostfunctionOptionsPLF

    end type

    ! General cost function options
    type, extends(optionsUDT)   :: CostfunctionOptionsSOUDT 

        ! Description
        !============
        ! General cost function options structure that contains fields
        ! for each possible cost function contribution. 

        ! Fields
        character(:), allocatable           :: type
        type(CostFunctionOptionsPLFUDT)     :: plf 

    contains 

        procedure :: SetDefaults    => SetDefaultCostfunctionOptionsSO
        procedure :: Read           => ReadCostfunctionOptionsSO

    end type

    !------------------------------------------------------------------!
    !                       SHAPE OPTIMIZATION                         !
    !------------------------------------------------------------------!

    ! General shape optimization options
    type, extends(OptionsUDT) :: ShapeOptimizationOptionsUDT

        ! Structure containing the options for shape optimization. The following 
        ! fields are present:

        ! General
        ! - filepath:   filepath to option file from which all options 
        !               are read in 
        ! - optionprefix:   prefix preceding any option key in the 
        !                   options key-value pairs
        
        ! Paths and options for I/O
        ! - goatfilepath:   input filepath to read goat options
        ! - writefilepath:  file path to which output should be written 
        !                   (without any extension - added separately)

        ! General
        character(:), allocatable   :: optionprefix, &
            goatfilepath, writefilepath  
    
    contains

        ! Routines to manipulate the options
        procedure   :: Read             => ReadShapeOptimizationOptions
        procedure   :: SetDefaults      => SetDefaultShapeOptimizationOptions
        

    end type

    ! Overarching design options
    type, extends(OptionsUDT) :: DesignOptionsSOUDT 

        ! Contains the different design option structures
        type(CostFunctionOptionsSOUDT)       :: costfunction
        type(DesignVariableOptionsSOUDT)     :: variables
        type(ConstraintOptionsSOUDT)         :: constraints

    contains 

        procedure :: SetDefaults    => SetDefaultDesignOptionsSO
        procedure :: Read           => ReadDesignOptionsSO

    end type

    contains 
        
    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                            Option setters                        !
    !------------------------------------------------------------------!

    ! Shape optimization 
    subroutine SetDefaultShapeOptimizationOptions(options)

        ! Description
        !============
        ! Default setter

        ! Declare variables
        !==================
        ! Arguments
        class(ShapeOptimizationOptionsUDT)  :: options
        
        ! Set
        !====
        ! Filepaths
        options%inputfilepath   = 'SOoptions.dat'
        options%writefilepath   = 'so'
        options%optionprefix    = 'so'
        options%goatfilepath    = 'GOAToptions.dat'

    end subroutine

    ! General design variables
    subroutine SetDefaultDesignVariableOptionsSO(options)

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariableOptionsSOUDT)   :: options 
        
        ! Set defaults
        !=============
        ! Own fields
        options%type = 'vesselcoordinates'

        ! Propagate inputfilepath
        options%vc%inputfilepath = options%inputfilepath

        ! Initialize
        !===========
        call options%vc%SetDefaults()

    end subroutine

    ! Vessel coordinates design variables
    subroutine SetDefaultDesignVariableOptionsVC(options)

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariableOptionsVCUDT)   :: options 

        ! Set defaults
        !=============

    end subroutine

    ! Fixed vessel points constraints
    subroutine SetDefaultFixedVesselPointsConOptions(options)

        ! Declare variables
        !==================
        ! Arguments
        class(FixedVesselPointsConOptionsUDT)   :: options 

        ! Set defaults
        !=============
        if (allocated(options%structureIDs)) then 
            deallocate(options%structureIDs)
        end if
        if (allocated(options%vertIDs)) then 
            deallocate(options%vertIDs)
        end if 
        allocate(options%structureIDs(0), options%vertIDs(0))

    end subroutine

    ! General constraints
    subroutine SetDefaultConstraintOptionsSO(options)

        ! Declare variables
        !==================
        ! Arguments
        class(ConstraintOptionsSOUDT)   :: options 

        ! Set defaults
        !=============
        ! Own fields
        options%fixedvesselpoints       = 0
        options%goat                    = 0

        ! Propagate filepaths
        options%fvpoptions%inputfilepath = options%inputfilepath

        ! Other constraints
        call options%fvpoptions%SetDefaults()

    end subroutine

    ! PLF cost function
    subroutine SetDefaultCostfunctionOptionsPLF(options)

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionOptionsPLFUDT)    :: options 

        ! Set defaults
        !=============
        options%lambda                  = 1
        options%vesselinputfilepath     = 'GOAToptions.dat'
        options%newvesselfilepath       = 'structure_new.dat'

    end subroutine

    ! General cost function
    subroutine SetDefaultCostfunctionOptionsSO(options)

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionOptionsSOUDT)   :: options 

        ! Set defaults
        !=============
        ! Set type
        options%type                = 'PLF'

        ! Propagate filepaths
        options%plf%inputfilepath = options%inputfilepath

        ! Contributions
        call options%plf%SetDefaults()

    end subroutine

    ! General design options
    subroutine SetDefaultDesignOptionsSO(options)

        ! Declare variables
        !==================
        ! Arguments
        class(DesignOptionsSOUDT)   :: options 

        ! Set defaults
        !=============
        ! Propagate filepaths
        options%costfunction%inputfilepath = options%inputfilepath
        options%constraints%inputfilepath = options%inputfilepath
        options%variables%inputfilepath = options%inputfilepath

        ! Contributions
        call options%costfunction%SetDefaults()
        call options%constraints%SetDefaults()
        call options%variables%SetDefaults()

    end subroutine


    !------------------------------------------------------------------!
    !                            Option readers                        !
    !------------------------------------------------------------------!

    ! Shape optimization 
    subroutine ReadShapeOptimizationOptions(options)

        ! Description
        !============
        ! Read in options from the specified inputfile

        ! Declare variables
        !==================
        ! Arguments
        class(ShapeOptimizationOptionsUDT)      :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: field
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
            print *, 'ReadShapeOptimizationOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadShapeOptimizationOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Filepaths etc
        field = 'so.main.goatfilepath'
        call ExtractOptionValueCharacter(fid, field, options%goatfilepath)
        field = 'so.main.writefilepath'
        call ExtractOptionValueCharacter(fid, field, options%writefilepath)
        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)


    end subroutine

    ! General design variables
    subroutine ReadDesignVariableOptionsSO(options)

        ! Arguments
        class(DesignVariableOptionsSOUDT)  :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: field
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
            print *, 'ReadDesignVariableOptionsSO: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadDesignVariableOptionsSO: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Design variable type
        field = 'so.dv.type'
        call ExtractOptionValueCharacter(fid, field, options%type)

        ! Read other structures
        call options%vc%Read()

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)


    end subroutine

    ! Vessel coordinates design variables
    subroutine ReadDesignVariableOptionsVC(options)

        ! Declare variables
        !==================
        ! Arguments
        class(DesignVariableOptionsVCUDT)   :: options 

        ! Read
        !=====

    end subroutine

    ! Fixed vessel points constraints
    subroutine ReadFixedVesselPointsConOptions(options)

        ! Declare variables
        !==================
        ! Arguments
        class(FixedVesselPointsConOptionsUDT)  :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: field
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
            print *, 'ReadFixedVesselPointsConOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadFixedVesselPointsConOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Structure & vertex IDs
        field = 'so.ec.par.fvp.structureIDs'
        call ExtractOptionValueInteger1D(fid, field, options%structureIDs)
        field = 'so.ec.par.fvp.vertIDs'
        call ExtractOptionValueInteger1D(fid, field, options%vertIDs)

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)


    end subroutine

    ! General constraints
    subroutine ReadConstraintOptionsSO(options)

        ! Declare variables
        !==================
        ! Arguments
        class(ConstraintOptionsSOUDT)   :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: field
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
            print *, 'ReadConstraintOptionsSO: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadConstraintOptionsSO: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Constraint switches
        field = 'so.ec.fixedvesselpoints'
        call ExtractOptionValueInteger0D(fid, field, options%fixedvesselpoints)
        field = 'so.ec.goat'
        call ExtractOptionValueInteger0D(fid, field, options%goat)

        ! Other constraint options
        call options%fvpoptions%Read()

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)


    end subroutine

    ! PLF cost function
    subroutine ReadCostFunctionOptionsPLF(options)

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionOptionsPLFUDT)    :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: field
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
            print *, 'ReadCostFunctionOptionsPLF: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadCostFunctionOptionsPLF: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Scaling constant
        field = 'so.cfv.par.PLF.lambda'
        call ExtractOptionValueReal0D(fid, field, options%lambda)

        ! Vessel path
        field  = 'so.cfv.par.PLF.vesselinputfilepath'
        call ExtractOptionValueCharacter(fid, field, options%vesselinputfilepath)
        field  = 'so.cfv.par.PLF.newvesselfilepath'
        call ExtractOptionValueCharacter(fid, field, options%newvesselfilepath)

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)

    end subroutine

    ! General cost function
    subroutine ReadCostfunctionOptionsSO(options)

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionOptionsSOUDT)   :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: field
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
            print *, 'ReadCostFunctionOptionsSO: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadCostFunctionOptionsSO: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Cost function type
        field = 'so.cfv.type'
        call ExtractOptionValueCharacter(fid, field, options%type)

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)
        
        ! Read other options
        !===================
        ! Contributions
        call options%plf%Read()

    end subroutine

    ! General design
    subroutine ReadDesignOptionsSO(options)

        ! Declare variables
        !==================
        ! Arguments
        class(DesignOptionsSOUDT)   :: options 
        
        ! Read options
        !=============
        ! Contributions
        call options%costfunction%Read()
        call options%variables%Read()
        call options%constraints%Read()

    end subroutine




end module somod_userinput