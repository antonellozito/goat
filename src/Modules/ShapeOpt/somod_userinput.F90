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
    use mod_constants, only: Pi_R8

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

    ! Fixed vessel flux constraints
    type, extends(optionsUDT)   :: FixedVesselFluxConOptionsUDT

        ! Description
        !============
        ! Constraint options for fixed vessel flux. The following
        ! fields are present:
        !   structureIDs:   structure IDs of structures of which the 
        !                   vertices should be constrained 
        !   vertIDs:        IDs of specific vessel vertices (according
        !                   to the original numbering in the 
        !                   structure.dat file) to be constrained
        ! 
        ! Note that vertices are constrained to their original psi value

        ! Fields
        integer(I8), allocatable    :: structureIDs(:), vertIDs(:)

    contains 

        procedure :: SetDefaults    => SetDefaultFixedVesselFluxConOptions
        procedure :: Read           => ReadFixedVesselFluxConOptions

    end type

    ! Vessel distance options
    type, extends(optionsUDT)   :: VesselDistanceConOptionsUDT

        ! Description
        !============
        ! Constraint options for vessel distance. The following
        ! fields are present:
        !   structureIDs:   structure IDs of structures of which the 
        !                   vertices should be constrained 
        !   vertIDs:        IDs of specific vessel vertices (according
        !                   to the original numbering in the 
        !                   structure.dat file) to be constrained
        !   d:              desired distance from polygon
        !   xp, yp:         polygon coordinates
        !   plftype:        type of polygon levelset function 
        !   resx, resy      options for closed polygon approximation
        !   C, M, meth
        !   offsetx, 
        !   offsety
        ! 

        ! Fields
        integer(I8)                 :: resx, resy, C, M
        integer(I8), allocatable    :: structureIDs(:), vertIDs(:)
        real(R8)                    :: offsetx, offsety, d 
        real(R8), allocatable       :: xp(:), yp(:)
        character(:), allocatable   :: meth, plftype
        
    contains 

        procedure :: SetDefaults    => SetDefaultVesselDistanceConOptions
        procedure :: Read           => ReadVesselDistanceConOptions

    end type

    ! Vessel upper bound options
    type, extends(VesselDistanceConOptionsUDT) :: VesselUBOptionsUDT

    contains 

        procedure :: SetDefaults    => SetDefaultVesselUpperBoundConOptions
        procedure :: Read           => ReadVesselUpperBoundConOptions

    end type 

    ! Vessel lower bound options
    type, extends(VesselDistanceConOptionsUDT) :: VesselLBOptionsUDT

    contains 

        procedure :: SetDefaults    => SetDefaultVesselLowerBoundConOptions
        procedure :: Read           => ReadVesselLowerBoundConOptions

    end type 

    ! Vessel angle difference options
    type, extends(optionsUDT)   :: VesselAngleDifferenceConOptionsUDT

        ! Description
        !============
        ! Constraint options for vessel distance. The following
        ! fields are present:
        !   structureIDs:   structure IDs of structures of which the 
        !                   vertices should be constrained 
        !   vertIDs:        IDs of specific vessel vertices (according
        !                   to the original numbering in the 
        !                   structure.dat file) to be constrained
        !   dthetas         bound for angle deviation (structure based)
        !   dthetav,        same, but vertID based (should have same size)
        !   alpha           minimal angle (scalar, one value)
        ! IMPORTANT: angles should be set in degrees in the input! They
        ! are converted to radians when reading in data.
        
        ! Fields
        integer(I8), allocatable    :: structureIDs(:), vertIDs(:)
        real(R8)                    :: alpha
        real(R8), allocatable       :: dthetas(:), dthetav(:)
        
    contains 

        procedure :: SetDefaults    => SetDefaultVesselAngleDifferenceConOptions
        procedure :: Read           => ReadVesselAngleDifferenceConOptions

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
        integer(I8)                 :: fixedvesselpoints, &
            fixedvesselflux, goat, vesselupperbound, vessellowerbound, &
            vesselangledifference

        ! Constraint options
        type(FixedVesselPointsConOptionsUDT)    :: fvpoptions
        type(FixedVesselFluxConOptionsUDT)      :: fvfoptions
        class(VesselDistanceConOptionsUDT), allocatable  :: vdoptions
        type(VesselUBOptionsUDT)                :: vduboptions 
        type(VesselLBOptionsUDT)                :: vdlboptions
        type(VesselAngleDifferenceConOptionsUDT)    :: vadoptions

    contains 

        procedure :: SetDefaults    => SetDefaultConstraintOptionsSO
        procedure :: Read           => ReadConstraintOptionsSO

    end type

    !------------------------------------------------------------------!
    !                       HESSIAN APPROXIMATION                      !
    !------------------------------------------------------------------!
    ! General hessian approximation options
    type, extends(optionsUDT)   :: HessianApproximationOptionsUDT

        ! Description
        !============
        ! Options for the hessian approximation of the cost function.
        ! See also the optmod_hessapprox module for more information on
        ! which hessian approximations etc are currently available. 

        ! The following fields are present:
        ! - storagetype:    type of hessian approximation storage ('dense', 
        !                   'sparse', or others)
        ! - updatemethod:   method to update the hessian (may depend on 
        !                   storage type, typically 'BFGS')
        ! - inithess:       initial hessian type ('diagonal')

        ! Options for 'diagonal' inithess:
        ! - diagind:        diagonal indices (relative to main diagonal) 
        !                   with non-zeros
        ! - diagval:        values on the diagonal (for uniform values on
        !                   separate diagonals, equal amount of values 
        !                   as diagonal indices, for non-uniform: 
        !                   specify equal to total number of elements)

        character(:), allocatable       :: storagetype, updatemethod, &
            inithess 
        integer(I8), allocatable        :: diagind(:)
        real(R8), allocatable           :: diagval(:)

    contains 

        procedure :: SetDefaults    => SetDefaultHessianApproximationOptionsSO
        procedure :: Read           => ReadHessianApproximationOptionsSO

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

    ! Face angle cost function
    type, extends(optionsUDT)   :: CostfunctionOptionsSOFAUDT

        ! Description
        !============
        ! Cost function options for the face angle cost function.
        ! This includes:
        !
        !   structureIDs:           IDs of the structure where to
        !                           impose the costfunction
        !   vertexIDs:              IDs of the vertices where to impose
        !                           the cost function (note: at least 
        !                           three subsequent IDs have to be 
        !                           present to have a cost function
        !                           contribution)
        !   lambda:                 scaling constant for cost function
        
        ! Fields
        real(R8)                    :: lambda
        integer(I8), allocatable    :: structureIDs(:), vertIDs(:)

    contains 

        procedure :: SetDefaults    => SetDefaultCostfunctionOptionsSOFA
        procedure :: Read           => ReadCostfunctionOptionsSOFA

    end type

    ! General cost function options
    type, extends(optionsUDT)   :: CostfunctionOptionsSOUDT 

        ! Description
        !============
        ! General cost function options structure that contains fields
        ! for each possible cost function contribution. Additionally, 
        ! there is the option to make a reduced cost function by 
        ! solving the goat equations for each cost function evaluation. 
        ! This is set by the logical 'dogoatreduction'. If a solps 
        ! contribution should be included, this should be set by
        ! the 'includesolps' switch. Note that this switch requires 
        ! also the goat reduction to be active

        ! Fields
        character(:), allocatable           :: type
        logical                             :: dogoatreduction, includesolps
        type(CostFunctionOptionsPLFUDT)     :: plf 
        type(CostfunctionOptionsSOFAUDT)    :: fa 
        type(HessianApproximationOptionsUDT)    :: hessapprox

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

        ! Remeshing options
        ! - doremesh    switch to turn remeshing off or on. If off, then 
        !               the program will simply exit on a remeshing 
        !               request. Hence, by default, this is turned on. 
        ! - itmaxremesh maximal number of remeshing steps that are 
        !               allowed. 
        
        ! Paths and options for I/O
        ! - goatfilepath:   input filepath to read goat options
        ! - writefilepath:  file path to which output should be written 
        !                   (without any extension - added separately)

        ! General
        character(:), allocatable   :: optionprefix, &
            goatfilepath, writefilepath  

        ! Remeshing
        logical                     :: doremesh 
        integer(I8)                 :: itmaxremesh
    
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

        ! Remeshing options
        options%doremesh        = .true. 
        options%itmaxremesh     = 10 

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

    ! Fixed vessel flux constraints
    subroutine SetDefaultFixedVesselFluxConOptions(options)

        ! Declare variables
        !==================
        ! Arguments
        class(FixedVesselFluxConOptionsUDT)   :: options 

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

    ! Fixed vessel distance constraints
    subroutine SetDefaultVesselDistanceConOptions(options)

        ! Declare variables
        !==================
        ! Arguments
        class(VesselDistanceConOptionsUDT)   :: options 

        ! Set defaults
        !=============
        if (allocated(options%structureIDs)) then 
            deallocate(options%structureIDs)
        end if
        if (allocated(options%vertIDs)) then 
            deallocate(options%vertIDs)
        end if 
        if (allocated(options%xp)) then 
            deallocate(options%xp)
        end if 
        if (allocated(options%yp)) then 
            deallocate(options%yp)
        end if
        allocate(options%structureIDs(0), options%vertIDs(0), &
            options%xp(0), options%yp(0))

        options%plftype = 'closedpolygon_exact'
        options%resx = 100
        options%resy = 100
        options%M = 6
        options%C = 3
        options%meth = 'uniform'
        options%offsetx = 0.05
        options%offsety = 0.05

    end subroutine

    ! Fixed vessel distance constraints, upper bound
    subroutine SetDefaultVesselUpperBoundConOptions(options)

        ! Declare variables
        !==================
        ! Arguments
        class(VesselUBOptionsUDT)   :: options 

        ! Set defaults
        !=============
        if (allocated(options%structureIDs)) then 
            deallocate(options%structureIDs)
        end if
        if (allocated(options%vertIDs)) then 
            deallocate(options%vertIDs)
        end if 
        if (allocated(options%xp)) then 
            deallocate(options%xp)
        end if 
        if (allocated(options%yp)) then 
            deallocate(options%yp)
        end if
        allocate(options%structureIDs(0), options%vertIDs(0), &
            options%xp(0), options%yp(0))

        options%plftype = 'closedpolygon_exact'
        options%resx = 100
        options%resy = 100
        options%M = 6
        options%C = 3
        options%meth = 'uniform'
        options%offsetx = 0.05
        options%offsety = 0.05

    end subroutine

    ! Fixed vessel distance constraints, lower bound
    subroutine SetDefaultVesselLowerBoundConOptions(options)

        ! Declare variables
        !==================
        ! Arguments
        class(VesselLBOptionsUDT)   :: options 

        ! Set defaults
        !=============
        if (allocated(options%structureIDs)) then 
            deallocate(options%structureIDs)
        end if
        if (allocated(options%vertIDs)) then 
            deallocate(options%vertIDs)
        end if 
        if (allocated(options%xp)) then 
            deallocate(options%xp)
        end if 
        if (allocated(options%yp)) then 
            deallocate(options%yp)
        end if
        allocate(options%structureIDs(0), options%vertIDs(0), &
            options%xp(0), options%yp(0))

        options%plftype = 'closedpolygon_exact'
        options%resx = 100
        options%resy = 100
        options%M = 6
        options%C = 3
        options%meth = 'uniform'
        options%offsetx = 0.05
        options%offsety = 0.05

    end subroutine

    ! Fixed vessel angle difference constraints
    subroutine SetDefaultVesselAngleDifferenceConOptions(options)

        ! Declare variables
        !==================
        ! Arguments
        class(VesselAngleDifferenceConOptionsUDT)   :: options 

        ! Set defaults
        !=============
        if (allocated(options%structureIDs)) then 
            deallocate(options%structureIDs)
        end if
        if (allocated(options%vertIDs)) then 
            deallocate(options%vertIDs)
        end if 
        if (allocated(options%dthetas)) then 
            deallocate(options%dthetas)
        end if 
        if (allocated(options%dthetav)) then 
            deallocate(options%dthetav)
        end if
        allocate(options%structureIDs(0), options%vertIDs(0), &
            options%dthetas(0), options%dthetav(0))

        options%alpha = 5.0_R8*Pi_R8/180.0_R8 ! in radians!

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
        options%fixedvesselflux         = 0 
        options%goat                    = 0
        options%vesselupperbound        = 0
        options%vessellowerbound        = 0 
        options%vesselangledifference   = 0

        ! Propagate filepaths
        options%fvpoptions%inputfilepath    = options%inputfilepath
        options%fvfoptions%inputfilepath    = options%inputfilepath
        options%vduboptions%inputfilepath   = options%inputfilepath 
        options%vdlboptions%inputfilepath   = options%inputfilepath
        options%vadoptions%inputfilepath    = options%inputfilepath

        ! Other constraints
        call options%fvpoptions%SetDefaults()
        call options%fvfoptions%SetDefaults()
        call options%vduboptions%SetDefaults()
        call options%vdlboptions%SetDefaults()
        call options%vadoptions%SetDefaults()

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

    ! PLF cost function
    subroutine SetDefaultCostfunctionOptionsSOFA(options)

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionOptionsSOFAUDT)    :: options 

        ! Set defaults
        !=============
        options%lambda                  = 1
        allocate(options%structureIDs(0), options%vertIDs(0))

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

        ! Set goat reduction switch
        options%dogoatreduction = .false. 
        options%includesolps    = .false.

        ! Propagate filepaths
        options%plf%inputfilepath = options%inputfilepath
        options%fa%inputfilepath = options%inputfilepath
        options%hessapprox%inputfilepath = options%inputfilepath

        ! Contributions
        call options%plf%SetDefaults()
        call options%fa%SetDefaults()
        call options%hessapprox%SetDefaults()

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
    
    ! General cost function
    subroutine SetDefaultHessianApproximationOptionsSO(options)

        ! Declare variables
        !==================
        ! Arguments
        class(HessianApproximationOptionsUDT)   :: options 

        ! Set defaults
        !=============
        ! Initial hessian (identity matrix)
        options%inithess            = 'diagonal'
        options%diagind             = [0]
        options%diagval             = [1.0_R8]
        options%updatemethod        = 'no' ! no updates
        options%storagetype         = 'sparse'

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

        ! Remeshing options
        field = 'so.remesh.doremesh'
        call ExtractOptionValueLogical0D(fid, field, options%doremesh)
        field = 'so.remesh.itmax'
        call ExtractOptionValueInteger0D(fid, field, options%itmaxremesh)
        
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

    ! Fixed vessel flux constraints
    subroutine ReadFixedVesselFluxConOptions(options)

        ! Declare variables
        !==================
        ! Arguments
        class(FixedVesselFluxConOptionsUDT)  :: options 

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
            print *, 'ReadFixedVesselFluxConOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadFixedVesselFluxConOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Structure & vertex IDs
        field = 'so.ec.par.fvf.structureIDs'
        call ExtractOptionValueInteger1D(fid, field, options%structureIDs)
        field = 'so.ec.par.fvf.vertIDs'
        call ExtractOptionValueInteger1D(fid, field, options%vertIDs)

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)


    end subroutine

    ! Vessel distance constraints
    subroutine ReadVesselDistanceConOptions(options)

        ! Declare variables
        !==================
        ! Arguments
        class(VesselDistanceConOptionsUDT)  :: options 

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
            print *, 'ReadVesselDistanceConOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadVesselDistanceConOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Structure & vertex IDs
        field = 'so.ec.par.vd.structureIDs'
        call ExtractOptionValueInteger1D(fid, field, options%structureIDs)
        field = 'so.ec.par.vd.vertIDs'
        call ExtractOptionValueInteger1D(fid, field, options%vertIDs)

        ! Plf options
        field = 'so.ec.par.vd.plftype'
        call ExtractOptionValueCharacter(fid, field, options%plftype)
        field = 'so.ec.par.vd.meth'
        call ExtractOptionValueCharacter(fid, field, options%meth)
        field = 'so.ec.par.vdlb.d'
        call ExtractOptionValueReal0D(fid, field, options%d)
        field = 'so.ec.par.vd.resx'
        call ExtractOptionValueInteger0D(fid, field, options%resx)
        field = 'so.ec.par.vd.resy'
        call ExtractOptionValueInteger0D(fid, field, options%resy)
        field = 'so.ec.par.vd.offsetx'
        call ExtractOptionValueReal0D(fid, field, options%offsetx)
        field = 'so.ec.par.vd.offsety'
        call ExtractOptionValueReal0D(fid, field, options%offsety)
        field = 'so.ec.par.vd.M'
        call ExtractOptionValueInteger0D(fid, field, options%M)
        field = 'so.ec.par.vd.C'
        call ExtractOptionValueInteger0D(fid, field, options%C)
        field = 'so.ec.par.vd.xp'
        call ExtractOptionValueReal1D(fid, field, options%xp)
        field = 'so.ec.par.vd.yp'
        call ExtractOptionValueReal1D(fid, field, options%yp)


        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)


    end subroutine

    ! Vessel upper boundary constraints
    subroutine ReadVesselUpperBoundConOptions(options)

        ! Declare variables
        !==================
        ! Arguments
        class(VesselUBOptionsUDT)  :: options 

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
            print *, 'ReadVesselUpperBoundaryConOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadVesselUpperBoundaryConOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Structure & vertex IDs
        field = 'so.inec.par.vdub.structureIDs'
        call ExtractOptionValueInteger1D(fid, field, options%structureIDs)
        field = 'so.inec.par.vdub.vertIDs'
        call ExtractOptionValueInteger1D(fid, field, options%vertIDs)

        ! Plf options
        field = 'so.inec.par.vdub.plftype'
        call ExtractOptionValueCharacter(fid, field, options%plftype)
        field = 'so.inec.par.vdub.meth'
        call ExtractOptionValueCharacter(fid, field, options%meth)
        field = 'so.inec.par.vdlb.d'
        call ExtractOptionValueReal0D(fid, field, options%d)
        field = 'so.inec.par.vdub.resx'
        call ExtractOptionValueInteger0D(fid, field, options%resx)
        field = 'so.inec.par.vdub.resy'
        call ExtractOptionValueInteger0D(fid, field, options%resy)
        field = 'so.inec.par.vdub.offsetx'
        call ExtractOptionValueReal0D(fid, field, options%offsetx)
        field = 'so.inec.par.vdub.offsety'
        call ExtractOptionValueReal0D(fid, field, options%offsety)
        field = 'so.inec.par.vdub.M'
        call ExtractOptionValueInteger0D(fid, field, options%M)
        field = 'so.inec.par.vdub.C'
        call ExtractOptionValueInteger0D(fid, field, options%C)
        field = 'so.inec.par.vdub.xp'
        call ExtractOptionValueReal1D(fid, field, options%xp)
        field = 'so.inec.par.vdub.yp'
        call ExtractOptionValueReal1D(fid, field, options%yp)

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)


    end subroutine

     ! Vessel upper boundary constraints
    subroutine ReadVesselLowerBoundConOptions(options)

        ! Declare variables
        !==================
        ! Arguments
        class(VesselLBOptionsUDT)  :: options 

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
            print *, 'ReadVesselLowerBoundaryConOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadVesselLowerBoundaryConOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Structure & vertex IDs
        field = 'so.inec.par.vdlb.structureIDs'
        call ExtractOptionValueInteger1D(fid, field, options%structureIDs)
        field = 'so.inec.par.vdlb.vertIDs'
        call ExtractOptionValueInteger1D(fid, field, options%vertIDs)

        ! Plf options
        field = 'so.inec.par.vdlb.plftype'
        call ExtractOptionValueCharacter(fid, field, options%plftype)
        field = 'so.inec.par.vdlb.meth'
        call ExtractOptionValueCharacter(fid, field, options%meth)
        field = 'so.inec.par.vdlb.resx'
        call ExtractOptionValueInteger0D(fid, field, options%resx)
        field = 'so.inec.par.vdlb.resy'
        call ExtractOptionValueInteger0D(fid, field, options%resy)
        field = 'so.inec.par.vdlb.d'
        call ExtractOptionValueReal0D(fid, field, options%d)
        field = 'so.inec.par.vdlb.offsetx'
        call ExtractOptionValueReal0D(fid, field, options%offsetx)
        field = 'so.inec.par.vdlb.offsety'
        call ExtractOptionValueReal0D(fid, field, options%offsety)
        field = 'so.inec.par.vdlb.M'
        call ExtractOptionValueInteger0D(fid, field, options%M)
        field = 'so.inec.par.vdlb.C'
        call ExtractOptionValueInteger0D(fid, field, options%C)


        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)


    end subroutine

    ! Vessel distance constraints
    subroutine ReadVesselAngleDifferenceConOptions(options)

        ! Declare variables
        !==================
        ! Arguments
        class(VesselAngleDifferenceConOptionsUDT)  :: options 

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
            print *, 'ReadVesselAngleDifferenceConOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadVesselAngleDifferenceConOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Structure & vertex IDs
        field = 'so.inec.par.vad.structureIDs'
        call ExtractOptionValueInteger1D(fid, field, options%structureIDs)
        field = 'so.inec.par.vad.vertIDs'
        call ExtractOptionValueInteger1D(fid, field, options%vertIDs)

        ! Angles
        field  = 'so.inec.par.vad.dthetas'
        call ExtractOptionValueReal1D(fid, field, options%dthetas)
        field  = 'so.inec.par.vad.dthetav'
        call ExtractOptionValueReal1D(fid, field, options%dthetav)
        field  = 'so.inec.par.vad.alpha'
        call ExtractOptionValueReal0D(fid, field, options%alpha)

        ! Convert to radians
        options%dthetas = options%dthetas*Pi_R8/180_R8
        options%dthetav = options%dthetav*Pi_R8/180_R8
        options%alpha = options%alpha*Pi_R8/180_R8

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
        field = 'so.ec.fixedvesselflux'
        call ExtractOptionValueInteger0D(fid, field, options%fixedvesselflux)
        field = 'so.ec.goat'
        call ExtractOptionValueInteger0D(fid, field, options%goat)
        field = 'so.inec.vesselupperbound'
        call ExtractOptionValueInteger0D(fid, field, options%vesselupperbound)
        field = 'so.inec.vessellowerbound'
        call ExtractOptionValueInteger0D(fid, field, options%vessellowerbound)
        field = 'so.inec.vesselangledifference'
        call ExtractOptionValueInteger0D(fid, field, options%vesselangledifference)

        ! Other constraint options
        call options%fvpoptions%Read()
        call options%fvfoptions%Read()
        call options%vduboptions%Read()
        call options%vdlboptions%Read()
        call options%vadoptions%Read()

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

    ! PLF cost function
    subroutine ReadCostFunctionOptionsSOFA(options)

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionOptionsSOFAUDT)    :: options 

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
            print *, 'ReadCostFunctionOptionsSOFA: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadCostFunctionOptionsSOFA: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Scaling constant
        field = 'so.cfv.par.FA.lambda'
        call ExtractOptionValueReal0D(fid, field, options%lambda)

        ! Structure & vertex IDs
        field = 'so.cfv.par.FA.structureIDs'
        call ExtractOptionValueInteger1D(fid, field, options%structureIDs)
        field = 'so.cfv.par.FA.vertIDs'
        call ExtractOptionValueInteger1D(fid, field, options%vertIDs)

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

        ! Goat reduction
        field = 'so.cfv.dogoatreduction'
        call ExtractOptionValueLogical0D(fid, field, options%dogoatreduction)

        ! Solps inclusion
        field = 'so.cfv.includesolps'
        call ExtractOptionValueLogical0D(fid, field, options%includesolps)

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)
        
        ! Read other options
        !===================
        ! Contributions
        call options%plf%Read()
        call options%fa%Read()
        call options%hessapprox%Read()

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

    ! General design
    subroutine ReadHessianApproximationOptionsSO(options)

         ! Declare variables
        !==================
        ! Arguments
        class(HessianApproximationOptionsUDT)    :: options 

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
            print *, 'ReadHessianApproximationOptionsSO: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadHessianApproximationOptionsSO: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Storage type etc
        field = 'so.cfv.par.hess.storagetype'
        call ExtractOptionValueCharacter(fid, field, options%storagetype)
        field = 'so.cfv.par.hess.inithess'
        call ExtractOptionValueCharacter(fid, field, options%inithess)
        field = 'so.cfv.par.hess.updatemethod'
        call ExtractOptionValueCharacter(fid, field, options%updatemethod)

        ! Initial hessian options for diagonal 
        field = 'so.cfv.par.hess.diagind'
        call ExtractOptionValueInteger1D(fid, field, options%diagind)
        field = 'so.cfv.par.hess.diagval'
        call ExtractOptionValueReal1D(fid, field, options%diagval)

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)

    end subroutine




end module somod_userinput