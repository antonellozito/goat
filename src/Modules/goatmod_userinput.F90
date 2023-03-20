!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains all options and structures to process user input.
! It uses the modules mod_readwrite and mod_inputfileparser to process
! the input provided in the input file. 

module goatmod_userinput

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

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                               Abstract                           !
    !------------------------------------------------------------------!


    ! Abstract option type
    type, abstract :: OptionsUDT  

        ! General abstract type for options. Should at least contain the
        ! path from where the options should be read. 
        character(:), allocatable       :: inputfilepath

    contains 

        procedure(ReadOptionsINT), deferred     :: Read 
        procedure(SetDefaultsINT), deferred     :: SetDefaults
        procedure                               :: Set => SetOptions
        procedure                               :: SetInputFile 

    end type 

    !------------------------------------------------------------------!
    !                               Goat                               !
    !------------------------------------------------------------------!

    ! General goat options
    type, extends(OptionsUDT) :: GoatoptionsUDT

        ! Structure containing the options for goat. The following 
        ! fields are present:
        ! - driver:     specify the driver to be used in the goat 
        !               program (see Goat.F90 for the options)
        ! - facelabelsubfrom:   labels in the original grid file to 
        !                       replace by other files
        ! - facelabelsubto:     replacement labels

        character(:), allocatable   :: driver ! driver to be taken for goat
        integer                     :: itmax 
        character(:), allocatable   :: filepath ! file path to options file

        ! Mappings
        integer(I8), allocatable    :: facelabelsubfrom(:) ! labels

    contains

        ! Routines to manipulate the options
        procedure   :: Read             => ReadGoatOptions
        procedure   :: SetDefaults      => SetDefaultGoatOptions
        

    end type

    ! Options for the grid
    type, extends(OptionsUDT) :: GridOptionsUDT

        ! Structure containing options for grid manipulation (mainly
        ! reading and writing). The following fields are present:
        ! - inputtype:      the type of inputfile. Can be 'b2fgmtry' 
        !                   or 'traduit'. Both are unstructured!
        ! - filepath:       file path indicating where the grid file
        !                   data is stored
        ! - facelabelsubfrom/facelabelsubto: substitute the labels 
        !                   in '..from' with the labels defined in 
        !                   '..to'. This can be used to substitute some
        !                   boundary labels with other labels if 
        !                   desired. 
        ! - facelabelmappingGG/facelabelmappingGD: 
        !                   defines the mapping between the 
        !                   face labels of the grid generator (GG) and 
        !                   the face labels of the grid deformation (GD)
        !                   see below for a mapping description 
        !                   
        ! Mapping description
        !====================
        ! GRID DEFORMATION            PHYSICAL MEANING
        !
        ! 1                           target plate (inner)
        ! 2                           target plate (outer)
        ! 3                           private flux
        ! 4                           core boundary
        ! 5                           outermost flux surf.
        !
        ! Note that information may be lost in this bndmapping, but this 
        ! is information that is not used by the grid deformation 
        ! module. 

        character(:), allocatable   :: inputtype
        character(:), allocatable   :: filepath 

        ! Mappings
        integer(I8), allocatable    :: facelabelsubfrom(:), &
            facelabelsubto(:), facelabelmappingGG(:), &
            facelabelmappingGD(:)
        

    contains 

        procedure :: Read               => ReadGridOptions
        procedure :: SetDefaults        => SetDefaultGridOptions
        
    end type  

    ! Options for magnetic field
    type, extends(OptionsUDT) :: MagneticFieldOptionsUDT

        ! Structure containing options for reading in the magnetic
        ! field. The following fields are present:
        ! - readmeth:   method to read in the magnetic field data. Only
        !               'readrzpsi' and 'default' supported for now. 
        !               These use the readrzpsi routine
        ! - filepath:   path where the magnetic field data is stored
        ! - interpC:    desired continuity of the interpolant 
        ! - interpM     order of the interpolant describing the values 
        !               at the grid nodes 
        ! - interpmeth  interpolant method ('uniformgrid' for data given 
        !               on a uniform grid, 'centered' for default 
        !               non-uniform grids. 

        
        character(:), allocatable   :: readmeth
        character(:), allocatable   :: filepath
        character(:), allocatable   :: interpmeth 

        integer(I8)                 :: interpC, interpM

    contains
    
        procedure :: Read           => ReadMagneticFieldOptions
        procedure :: SetDefaults    => SetDefaultMagneticFieldOptions
        
    end type  

    ! Options for the vessel
    type, extends(OptionsUDT) :: VesselOptionsUDT

        ! Fields:
        ! - readmeth:   'read_structure' for structure.dat files. Only
        !               method that is currently supported.
        ! - filepath:   path to the file to be read
        ! - refine:     set to 1 to do refinement of vessel (insert 
        !               more nodes)
        ! - maxdist:    maximum distance between two nodes. If larger,
        !               nodes will be added in between when refine == 1
        ! - TP:         integer giving the indices of which part of the
        !               structure should be considered as target plates.
        !               This should not account for any structures being
        !               excluded using 'exclude' (see below)
        ! - TPind:      target plate enumerators 
        ! - exclude:    list of structures in the vessel to be excluded
        !               from the vessel 

        character(:), allocatable       :: readmeth
        character(:), allocatable       :: filepath

        real(R8)                        :: maxdist
        integer(I8)                     :: refine
        integer(I8), allocatable        :: TP(:), TPind(:), exclude(:)

    contains 

        procedure :: Read           => ReadVesselOptions
        procedure :: SetDefaults    => SetDefaultVesselOptions
        
    end type 

    ! Options for the environment
    type, extends(OptionsUDT) :: EnvironmentOptionsUDT

        ! Structure used to keep all other possible necessary 
        ! structures that do not immediately fall under the grid, vessel
        ! or magnetic field. Currently empty. 

        character(:), allocatable   :: type
        character(:), allocatable   :: vesselfilepath

    contains 

        procedure :: Read           => ReadEnvironmentOptions
        procedure :: SetDefaults    => SetDefaultEnvironmentOptions
        
    end type

    !------------------------------------------------------------------!
    !                         Grid deformation                         !
    !------------------------------------------------------------------!

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    ! Abstract interfaces
    abstract interface

        subroutine ReadOptionsINT(options)
            import :: OptionsUDT
            class(OptionsUDT) :: options 

        end subroutine

        subroutine SetDefaultsINT(options)
            import :: OptionsUDT
            class(OptionsUDT) :: options 

        end subroutine

    end interface

    contains 

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                            Option setters                        !
    !------------------------------------------------------------------!

    ! Main option setter
    subroutine SetOptions(options)

        ! Description
        !============
        ! Set the options by first setting the defaults and overriding
        ! them later on with the user-defined values from the input 
        ! file. It is assumed that the inputfile is set beforehand. This
        ! inputfilepath should not be overridden in the defaults! Use 
        ! the 'SetInputFile' routine in the setup to determine the
        ! input file to be read. 

        ! Declare variables
        !==================
        class(OptionsUDT)       :: options 

        ! Set options
        !============
        ! Defaults
        call options%SetDefaults()

        ! Read
        call options%Read()

    end subroutine

    ! Main option inputfile path setter
    subroutine SetInputFile(options, path)

        ! Description
        !============ 
        ! Set the option inputfilepath

        ! Declare variables
        !==================
        ! Arguments
        class(OptionsUDT)                       :: options
        character(:), allocatable, intent(in)   :: path 

        ! Set path
        !=========
        options%inputfilepath = path

    end subroutine

    ! Goat options routines
    subroutine SetDefaultGoatOptions(options)

        ! Declare variables
        !==================
        ! Arguments
        class(GoatoptionsUDT)       :: options
        character(:), allocatable   :: filepath
        
        ! Set default values
        !===================
        ! Input file
        options%filepath    = './Examples/TCV/GOAToptions.dat'

        ! Numerics
        options%itmax       = 5
        options%driver      = 'GD'

        ! Mappings
        allocate(options%facelabelsubfrom(0))

    end subroutine

    ! Grid option routines
    subroutine SetDefaultGridOptions(options)

        ! Description
        !============
        ! Set the options for the grid

        ! Declaration
        class(GridOptionsUDT)        ::  options

        ! Default options
        options%gridtype            = 'plasma'
        options%inputtype           = 'b2fgmtry'

        ! Default grid location /data/le/Examples/TCV/b2fgmtry_us_tcv
        options%filepath            = './Examples/TCV/b2fgmtry_us_tcv'
    
        ! Default mappings
        allocate(options%facelabelsubfrom(0), options%facelabelsubto(0))
        allocate(options%facelabelmappingGG(8), options%facelabelmappingGD(8))
        options%facelabelmappingGG = [-13, -34, -23, -24, -21, -42, -43, -44]
        options%facelabelmappingGD = [1, 2, 3,   3,   4,   5,   5,   5]

    
    end subroutine

    ! Magnetic field
    subroutine SetDefaultMagneticFieldOptions(options)
        ! Description
        !============
        ! Set the options for the magnetic field. The following fields
        ! have to be set:

        ! Declaration
        class(MagneticFieldOptionsUDT)       :: options

        ! Default options   
        options%readmeth                = 'default'
        options%filepath                = './Examples/TCV/rzpsi_tcv.dat'

        options%interpmeth              = 'uniformgrid' 
        options%interpC                 = 3
        options%interpM                 = 6

    end subroutine

    ! Vessel
    subroutine SetDefaultVesselOptions(options)

        ! Description
        !============
        ! Set the options for the vessel.

        ! Declare variables
        !==================
        class(VesselOptionsUDT)          :: options

        ! Set options
        !============
        ! Input
        options%readmeth    = 'read_structure'
        options%filepath    = './Examples/TCV/structure_tcv.dat'

        ! Refinement options
        options%refine      = 1
        options%maxdist     = 0.01

        ! Target plates
        allocate(options%TP(2))
        allocate(options%TPind(2))
        allocate(options%exclude(0))
        options%TP      = [1, 2]
        options%TPind   = [1, 2]
        

    end subroutine

    ! Environment
    subroutine SetDefaultEnvironmentOptions(options)

        ! Description
        !============
        ! Set the options to read in the environment.

        ! Declare variables
        !==================
        class(EnvironmentOptionsUDT)         :: options

        ! Set options
        !============
        options%type = 'vessel'
        options%vesselfilepath = options%inputfilepath

    end subroutine

    !------------------------------------------------------------------!
    !                            Option readers                        !
    !------------------------------------------------------------------!

    ! Goat options reader
    subroutine ReadGoatOptions(options)

        ! Description
        !============
        ! This routine reads in the goat options from a file of which 
        ! the full path should be given in options%inputfilepath. The default
        ! options should have already been set at this point, as this 
        ! routine will only overwrite options that are present in the 
        ! user-specified options file. If no options file is present, 
        ! nothing is read in and a message will be shown. 

        ! Declare variables
        !==================
        ! Arguments
        class(GoatoptionsUDT)            :: options 

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
            print *, 'ReadGoatOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadGoatOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Driver
        field = 'GOAToptions.driver'
        call ExtractOptionValueCharacter(fid, field, options%driver)

        ! Max number of iterations
        field = 'GOAToptions.itmax'
        call ExtractOptionValueInteger0D(fid, field, options%itmax)

        ! Mappings
        field = 'GOAToptions.GDtoGA.facelabelsubfrom'
        call ExtractOptionValueInteger1D(fid, field, &
            options%facelabelsubfrom)
        

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)


    end subroutine

    ! Grid options reader
    subroutine ReadGridOptions(options)

        ! Description
        !============
        ! Read in grid options from file. It is assumed that the 
        ! filepath has been set correctly. 

        ! Declare variables
        !==================
        ! Arguments
        class(GridOptionsUDT)            :: options 

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
            print *, 'ReadGridOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadGridOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Inputtype
        field = 'gd.grid.inputtype'
        call ExtractOptionValueCharacter(fid, field, options%inputtype)

        ! Filepath
        field = 'gd.grid.gridfilepath'
        call ExtractOptionValueCharacter(fid, field, options%filepath)

        ! Mappings
        field = 'gd.grid.facelabelsubfrom'
        call ExtractOptionValueInteger1D(fid, field, options%facelabelsubfrom)
        field = 'gd.grid.facelabelsubto'
        call ExtractOptionValueInteger1D(fid, field, options%facelabelsubto)
        field = 'gd.grid.facelabelmappingGG'
        call ExtractOptionValueInteger1D(fid, field, options%facelabelmappingGG)
        field = 'gd.grid.facelabelmappingGD'
        call ExtractOptionValueInteger1D(fid, field, options%facelabelmappingGD)

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)

    end subroutine

    ! Magnetic field options reader
    subroutine ReadMagneticFieldOptions(options)

        ! Description
        !============
        ! Read in grid options from file. It is assumed that the 
        ! filepath has been set correctly. 

        ! Declare variables
        !==================
        ! Arguments
        class(MagneticFieldOptionsUDT)            :: options 

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
            print *, 'ReadMagneticFieldOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadMagneticFieldOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! I/O
        field = 'gd.mf.readmeth'
        call ExtractOptionValueCharacter(fid, field, options%readmeth)
        field = 'gd.mf.filepath'
        call ExtractOptionValueCharacter(fid, field, options%filepath)

        ! Interpolant
        field = 'gd.mf.interpC'
        call ExtractOptionValueInteger0D(fid, field, options%interpC)
        field = 'gd.mf.interpM'
        call ExtractOptionValueInteger0D(fid, field, options%interpM)
        field = 'gd.mf.interpmeth'
        call ExtractOptionValueCharacter(fid, field, options%interpmeth)

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)

    end subroutine

    ! Vessel options reader
    subroutine ReadVesselOptions(options)

        ! Description
        !============
        ! Read in vessel options from file. It is assumed that the 
        ! filepath has been set correctly. 

        ! Declare variables
        !==================
        ! Arguments
        class(VesselOptionsUDT)         :: options 

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
            print *, 'ReadVesselOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadVesselOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! I/O
        field = 'gd.grid.vesselloadmeth'
        call ExtractOptionValueCharacter(fid, field, options%readmeth)
        field = 'gd.grid.vesselfilepath'
        call ExtractOptionValueCharacter(fid, field, options%filepath)
        field = 'gd.grid.exclude'
        call ExtractOptionValueInteger1D(fid, field, options%exclude)

        ! Refinement
        field = 'gd.grid.refinevessel'
        call ExtractOptionValueInteger0D(fid, field, options%refine)
        field = 'gd.grid.maxvesseldist'
        call ExtractOptionValueReal0D(fid, field, options%maxdist)
        
        ! Target plates
        field = 'gd.grid.TP'
        call ExtractOptionValueInteger1D(fid, field, options%TP)
        field = 'gd.grid.TPind'
        call ExtractOptionValueInteger1D(fid, field, options%TPind)

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)

    end subroutine

    ! Environment options reader
    subroutine ReadEnvironmentOptions(options)

        ! Description
        !============
        ! Read in environment options from file. It is assumed that the 
        ! filepath has been set correctly. 

        ! Declare variables
        !==================
        ! Arguments
        class(EnvironmentOptionsUDT)    :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: thisline, field, &
            loaddir,loadfile
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
            print *, 'ReadEnvironmentOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadEnvironmentOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Nothing to be read in currently

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)

    end subroutine


end module goatmod_userinput