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
    !                               GOAT                               !
    !------------------------------------------------------------------!

    ! General goat options
    type, extends(OptionsUDT) :: GoatoptionsUDT

        ! Structure containing the options for goat. The following 
        ! fields are present:

        ! General fields:
        ! - debug: general debug plots on the goat level
        ! - meth: goat running method. Currently, only 'GD' is supported
        ! - gridreadtype: type of grid input file to read. Can be 
        ! 'traduitb2us' or b2fgmtry (also in unstructured variant)
        ! - magneticfieldreadtype: type of magnetic field input to read.
        ! Can be 'rzpsi' or 'equ' for rzpsi.dat or equ.dat files.
        ! - filepath: path towards the file where options are defined
        ! - gdfilepath: path towards the file where options for grid 
        ! deformation are defined

        ! Input filenames:
        ! - gridfilepath: file path to file with grid data (e.g. traduit.out.b2us file)
        ! - structurefilepath: structure.dat file to read
        ! - magneticfieldfilepath: magnetic field file to read
        ! - writefilepath: path where to write output traduit file

        ! Output options
        ! - write_final: write final output
        ! - write_traduitb2us: write unstructured traduit file
        ! - write_b2agdat:  write final b2ag.dat file for use in b2ag
        ! - write_Xpointdata: write out X-point data in traduit file
        ! - write_OMPdata: write OMP data in traduit file

        ! Case identification options
        ! - vesselmode: set to true if the case is a vessel mode grid
        ! - slab: true if slab grid
        ! - artificial_slab: true if artificial slab

        ! Face label mappings
        ! - GAtoGDfacelabelmappingGG: labels as defined in grid
        ! generator for interfacing between GA and GD
        ! - GAtoGDfacelabelmappingGD: corresponding labels for GD (so 
        ! first GG label is mapped to first GD label here)
        ! - GAtoGDfacelabelsubfrom: substitution of this face label ... 
        ! - GAtoGDfacelabelsubto: ... to this face label in GA to GD 
        ! interface

        ! - GGtoGDfacelabelmappingGG: idem above but for GG to GD 
        ! - GGtoGDfacelabelmappingGA
        ! - GGtoGDfacelabelsubfrom
        ! - GGtoGDfacelabelsubto

        ! Structure options
        ! - TP: structure numbers that are target plates
        ! - TPind: indices of target plates
        ! - exclude: indices of structures in structure.dat that should 
        ! be excluded when generating the bounding polygon

        ! OMP and IMP definition
        ! - OMP_r, OMP_z: R, Z coordinates that define the outer mid 
        ! plane line segment
        ! - IMP_r, IMP_z: R, Z coordinates that define the inner mid 
        ! plane line segment

        ! General
        logical                     :: debug     
        character(:), allocatable   :: meth 
        character(:), allocatable   :: magneticfieldreadtype
        character(:), allocatable   :: filepath
        character(:), allocatable   :: gdinputfilepath

        ! Specify input filenames
        character(:), allocatable   :: gridfilepath
        character(:), allocatable   :: structurefilepath
        character(:), allocatable   :: magneticfieldfilepath
        character(:), allocatable   :: writefilepath

        ! Output options
        logical                     :: write_final 
        logical                     :: write_traduitb2us
        logical                     :: write_b2agdat
        logical                     :: write_Xpointdata 
        logical                     :: write_OMPdata

        ! Case identification options
        logical                     :: vesselmode 
        logical                     :: slab 
        logical                     :: artificial_slab

        ! Face label mappings
        integer(I8), allocatable    :: GDtoGAfacelabelmappingGG(:)
        integer(I8), allocatable    :: GDtoGAfacelabelmappingGD(:) 
        integer(I8), allocatable    :: GDtoGAfacelabelsubfrom(:) 
        integer(I8), allocatable    :: GDtoGAfacelabelsubto(:) 

        integer(I8), allocatable    :: GGtoGDfacelabelmappingGG(:)
        integer(I8), allocatable    :: GGtoGDfacelabelmappingGD(:) 
        integer(I8), allocatable    :: GGtoGDfacelabelsubfrom(:) 
        integer(I8), allocatable    :: GGtoGDfacelabelsubto(:) 

        ! Structure options
        integer(I8), allocatable    :: TP(:)
        integer(I8), allocatable    :: TPind(:) 
        integer(I8), allocatable    :: exclude(:)

        ! OMP and IMP
        real(R8), allocatable       :: OMP_r(:), OMP_z(:), IMP_r(:), &
            IMP_z(:)

    contains

        ! Routines to manipulate the options
        procedure   :: Read             => ReadGoatOptions
        procedure   :: SetDefaults      => SetDefaultGoatOptions
        

    end type

    !------------------------------------------------------------------!
    !                          GRID DEFORMATION                        !
    !------------------------------------------------------------------!

    ! General grid deformation options
    type, extends(OptionsUDT) :: GDoptionsUDT

        ! Structure containing the options for the grid deformation part 
        ! of GOAT. The following fields are present:

        ! - runtype: type of run for grid deformation. currently, only 
        ! 'optimize' is available, though may be extended in the future
        ! - gridtype: only available options is plasma edge grid 
        ! ('plasma')
        ! - meth: only supported option is 'KKT', though augmented 
        ! lagrangian is also available. Method to solve optimization 
        ! problem when using 'optimize' as runtype

        ! - designoptionsfile: file to read design options
        ! - gridoptionsfile: similar to above, but for grid
        ! - magneticfieldoptionsfile: for MF
        ! - numparamsoptionsfile: for numerical parameters
        ! - environmentoptionsfile: for other environment stuff

        ! 
        ! - driver:     specify the driver to be used in the goat 
        !               program (see Goat.F90 for the options)
        ! - facelabelsubfrom:   labels in the original grid file to 
        !                       replace by other files
        ! - facelabelsubto:     replacement labels

        ! General deformation options
        character(:), allocatable   :: runtype 
        character(:), allocatable   :: gridtype 
        character(:), allocatable   :: meth 
        character(:), allocatable   :: filepath 

        ! Files to read other options
        character(:), allocatable   :: designoptionsfile 
        character(:), allocatable   :: gridoptionsfile 
        character(:), allocatable   :: magneticfieldoptionsfile 
        character(:), allocatable   :: numparamsoptionsfile
        character(:), allocatable   :: environmentoptionsfile

    contains

        ! Routines to manipulate the options
        procedure   :: Read             => ReadGDOptions
        procedure   :: SetDefaults      => SetDefaultGDOptions
        

    end type

    ! Options for the grid
    type, extends(OptionsUDT) :: GridOptionsUDT

        ! Structure containing options for grid manipulation (mainly
        ! reading and writing). The following fields are present:
        ! - type:       the type of grid that is considered. Only
        !                   'plasma' is currently available. 
        ! - readmeth:       the type of inputfile. Can be 'b2fgmtry' 
        !                   or 'traduit'. Both are unstructured!
        ! - filepath:       file path indicating where the grid file
        !                   data is stored

        ! Input filepaths and types
        character(:), allocatable   :: type 
        character(:), allocatable   :: readmeth 
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
        ! - shapemeth:  method used to represent (and possibly 
        ! approximate) the shape of the vessel. 
        ! - resx:       resolution in x-direction for shape rep
        ! - resy:       resolution in y-direction for shape rep
        ! - offsetfracx:    fractional offset to be taken from minimal 
        ! and maximal x-value of original vessel polygon
        ! - offsetfracy:    same as offsetfracx, but for y-direction
        ! - interpC:    only used for interpolant based methods. Order
        ! of interpolant
        ! - interpM:    only used for interpolant based methods. Order 
        ! used to approximate derivatives to construct interpolant 
        ! function. 

        character(:), allocatable       :: readmeth
        character(:), allocatable       :: filepath

        real(R8)                        :: maxdist
        integer(I8)                     :: refine
        integer(I8), allocatable        :: TP(:), TPind(:), exclude(:)

        ! Vessel representation options
        character(:), allocatable       :: shapemeth 
        integer(I8)                     :: resx, resy, interpC, interpM 
        real(R8)                        :: offsetfracx, offsetfracy

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
        character(:), allocatable   :: filepath
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
        
        ! Input file
        options%inputfilepath    = './GOAToptions.dat'

        ! General
        options%debug           = .false. 
        options%meth            = 'GD'
        options%gdinputfilepath = './GOAToptions.dat'

        ! Specify input filenames
        options%gridfilepath            = './traduit.out.b2us'
        options%structurefilepath       = './structure.dat'
        options%magneticfieldfilepath   = './rzpsi.dat'

        ! Output options
        options%writefilepath       = 'traduit.out.b2us_smoothed'
        options%write_final         = .true. 
        options%write_traduitb2us   = .true.
        options%write_b2agdat       = .true. 
        options%write_Xpointdata    = .false. 
        options%write_OMPdata       = .false. 

        ! Case identification options
        options%vesselmode          = .false. 
        options%slab                = .false.
        options%artificial_slab     = .false.
        
        ! Face label mappings
        allocate(options%GDtoGAfacelabelmappingGG(0), &
            options%GDtoGAfacelabelmappingGD(0), &
            options%GDtoGAfacelabelsubfrom(0), &
            options%GDtoGAfacelabelsubto(0), &
            options%GGtoGDfacelabelmappingGG(0), &
            options%GGtoGDfacelabelmappingGD(0), &
            options%GGtoGDfacelabelsubfrom(0), &
            options%GGtoGDfacelabelsubto(0))
        
        ! Structure options
        allocate(options%TP(0), options%TPind(0), options%exclude(0))

        ! OMP and IMP
        allocate(options%OMP_r(1:2), options%OMP_z(1:2), &
            options%IMP_r(1:2), options%IMP_z(1:2))
        options%OMP_r(:)   = 0
        options%OMP_z(:)   = 0
        options%IMP_r(:)   = 0
        options%IMP_z(:)   = 0

    end subroutine

    ! Grid deformation options routines
    subroutine SetDefaultGDOptions(options)

        ! Declare variables
        !==================
        ! Arguments
        class(GDoptionsUDT)       :: options        
        
        ! Input file
        options%filepath    = './GOAToptions.dat'

        ! General
        options%runtype     = 'optimize'
        options%gridtype    = 'plasma'
        options%meth        = 'KKT'

        ! Files to read other options
        options%designoptionsfile           = './GOAToptions.dat'
        options%gridoptionsfile             = './GOAToptions.dat'
        options%magneticfieldoptionsfile    = './GOAToptions.dat'
        options%numparamsoptionsfile        = './GOAToptions.dat'
        options%environmentoptionsfile      = './GOAToptions.dat'
        
    end subroutine

    ! Grid option routines
    subroutine SetDefaultGridOptions(options)

        ! Description
        !============
        ! Set the options for the grid

        ! Declaration
        class(GridOptionsUDT)        ::  options

        ! Default options
        options%type                = 'plasma'
        options%readmeth            = 'b2fgmtry'

        ! Default grid location
        options%filepath            = 'traduit.out.b2us'

        ! Default mappings
        allocate(options%facelabelsubfrom(0), options%facelabelsubto(0))
        allocate(options%facelabelmappingGG(1:8), options%facelabelmappingGD(1:8))
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
        options%filepath    = './structure.dat'

        ! Refinement options
        options%refine      = 1
        options%maxdist     = 0.01

        ! Target plates
        allocate(options%TP(2))
        allocate(options%TPind(2))
        allocate(options%exclude(0))
        options%TP      = [1, 2]
        options%TPind   = [1, 2]

        ! Vessel shape representation
        allocate(character(len('closedpolygon_smoothapproximation')) :: &
            options%shapemeth)
        options%shapemeth = 'closedpolygon_smoothapproximation'
        options%resx            = 100
        options%resy            = 100 
        options%offsetfracx     = 0.1
        options%offsetfracy     = 0.1
        options%interpC         = 3
        options%interpM         = 6
        

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

        ! Set data

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
            print *, 'ReadGoatOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadGoatOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! General
        field = 'goat.debug'
        call ExtractOptionValueLogical0D(fid, field, options%debug)
        field = 'goat.meth'
        call ExtractOptionValueCharacter(fid, field, options%meth)

        ! Input filenames
        field = 'goat.gridfilepath'
        call ExtractOptionValueCharacter(fid, field, options%gridfilepath)
        field = 'goat.structurefilepath'
        call ExtractOptionValueCharacter(fid, field, options%structurefilepath)
        field = 'goat.magneticfieldfilepath'
        call ExtractOptionValueCharacter(fid, field, options%magneticfieldfilepath)
        field = 'goat.GDinputfilepath'
        call ExtractOptionValueCharacter(fid, field, options%gdinputfilepath)

        ! Output options
        field = 'goat.writefilepath'
        call ExtractOptionValueCharacter(fid, field, options%writefilepath)
        field = 'goat.write_final'
        call ExtractOptionValueLogical0D(fid, field, options%write_final)
        field = 'goat.write_traduitb2us'
        call ExtractOptionValueLogical0D(fid, field, options%write_traduitb2us)
        field = 'goat.write_b2agdat'
        call ExtractOptionValueLogical0D(fid, field, options%write_b2agdat)
        field = 'goat.write_Xpointdata'
        call ExtractOptionValueLogical0D(fid, field, options%write_Xpointdata)
        field = 'goat.write_OMPdata'
        call ExtractOptionValueLogical0D(fid, field, options%write_OMPdata)

        ! Case identification options
        field = 'goat.vesselmode'
        call ExtractOptionValueLogical0D(fid, field, options%vesselmode)
        field = 'goat.slab'
        call ExtractOptionValueLogical0D(fid, field, options%slab)
        field = 'goat.artificial_slab'
        call ExtractOptionValueLogical0D(fid, field, options%artificial_slab)

        ! Face label mappings
        field = 'goat.GDtoGA.facelabelmappingGG'
        call ExtractOptionValueInteger1D(fid, field, &
            options%GDtoGAfacelabelmappingGG)
        field = 'goat.GDtoGA.facelabelmappingGD'
        call ExtractOptionValueInteger1D(fid, field, &
            options%GDtoGAfacelabelmappingGD) 
        field = 'goat.GDtoGA.facelabelsubfrom'
        call ExtractOptionValueInteger1D(fid, field, &
            options%GDtoGAfacelabelsubfrom)
        field = 'goat.GDtoGA.facelabelsubto'
        call ExtractOptionValueInteger1D(fid, field, &
            options%GDtoGAfacelabelsubto)

        field = 'goat.GGtoGD.facelabelmappingGG'
        call ExtractOptionValueInteger1D(fid, field, &
            options%GGtoGDfacelabelmappingGG)
        field = 'goat.GGtoGD.facelabelmappingGD'
        call ExtractOptionValueInteger1D(fid, field, &
            options%GGtoGDfacelabelmappingGD)
        field = 'goat.GGtoGD.facelabelsubfrom'
        call ExtractOptionValueInteger1D(fid, field, &
            options%GGtoGDfacelabelsubfrom)
        field = 'goat.GGtoGD.facelabelsubto'
        call ExtractOptionValueInteger1D(fid, field, &
            options%GGtoGDfacelabelsubto)

        ! OMP and IMP
        field = 'goat.OMP_r'
        call ExtractOptionValueReal1D(fid, field, &
            options%OMP_r)
        field = 'goat.OMP_z'
        call ExtractOptionValueReal1D(fid, field, &
            options%OMP_z)
        field = 'goat.IMP_r'
        call ExtractOptionValueReal1D(fid, field, &
            options%IMP_r)
        field = 'goat.IMP_z'
        call ExtractOptionValueReal1D(fid, field, &
            options%IMP_z)

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)


    end subroutine

    ! Grid deformation options reader
    subroutine ReadGDOptions(options)

        ! Description
        !============
        ! This routine reads in the grid deformation options from a file of which 
        ! the full path should be given in options%inputfilepath. The default
        ! options should have already been set at this point, as this 
        ! routine will only overwrite options that are present in the 
        ! user-specified options file. If no options file is present, 
        ! nothing is read in and a message will be shown. 

        ! Declare variables
        !==================
        ! Arguments
        class(GDoptionsUDT)             :: options 

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
            print *, 'ReadGDOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadGDOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! General
        field = 'gd.main.runtype'
        call ExtractOptionValueCharacter(fid, field, options%runtype)
        field = 'gd.main.gridtype'
        call ExtractOptionValueCharacter(fid, field, options%gridtype)
        field = 'gd.main.meth'
        call ExtractOptionValueCharacter(fid, field, options%meth)

        ! Files to read other options
        field = 'gd.main.designoptionsfile'
        call ExtractOptionValueCharacter(fid, field, options%designoptionsfile)
        field = 'gd.main.gridoptionsfile'
        call ExtractOptionValueCharacter(fid, field, options%gridoptionsfile)
        field = 'gd.main.magneticfieldoptionsfile'
        call ExtractOptionValueCharacter(fid, field, options%magneticfieldoptionsfile)
        field = 'gd.main.numparamsoptionsfile'
        call ExtractOptionValueCharacter(fid, field, options%numparamsoptionsfile)
        field = 'gd.main.environmentoptionsfile'
        call ExtractOptionValueCharacter(fid, field, options%environmentoptionsfile)
        
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
            print *, 'ReadGridOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadGridOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Gridtype
        field = 'goat.grid.type'
        call ExtractOptionValueCharacter(fid, field, options%type)

        ! Reading method
        field = 'goat.grid.readmeth'
        call ExtractOptionValueCharacter(fid, field, options%readmeth)

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
        character(:), allocatable       :: field
        integer, parameter              :: fid = 10 
        logical                         :: reachedeof

        ! Initialize
        !===========
        ! Variables
        reachedeof = .false. 

        ! Open the file, check if it exists
        print *, options%inputfilepath
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
        field = 'goat.mf.readmeth'
        call ExtractOptionValueCharacter(fid, field, options%readmeth)

        ! Interpolant
        field = 'goat.mf.interpC'
        call ExtractOptionValueInteger0D(fid, field, options%interpC)
        field = 'goat.mf.interpM'
        call ExtractOptionValueInteger0D(fid, field, options%interpM)
        field = 'goat.mf.interpmeth'
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
        field = 'goat.vessel.readmeth'
        call ExtractOptionValueCharacter(fid, field, options%readmeth)
        field = 'goat.vessel.exclude'
        call ExtractOptionValueInteger1D(fid, field, options%exclude)

        ! Refinement
        field = 'goat.vessel.refinevessel'
        call ExtractOptionValueInteger0D(fid, field, options%refine)
        field = 'goat.vessel.maxvesseldist'
        call ExtractOptionValueReal0D(fid, field, options%maxdist)
        
        ! Target plates
        field = 'goat.vessel.TP'
        call ExtractOptionValueInteger1D(fid, field, options%TP)
        field = 'goat.vessel.TPind'
        call ExtractOptionValueInteger1D(fid, field, options%TPind)

        ! Shape representation options
        field = 'goat.vessel.shapemeth'
        call ExtractOptionValueCharacter(fid, field, options%shapemeth)
        field = 'goat.vessel.resx'
        call ExtractOptionValueInteger0D(fid, field, options%resx)
        field = 'goat.vessel.resy'
        call ExtractOptionValueInteger0D(fid, field, options%resy)
        field = 'goat.vessel.interpC'
        call ExtractOptionValueInteger0D(fid, field, options%interpC)
        field = 'goat.vessel.interpM'
        call ExtractOptionValueInteger0D(fid, field, options%interpM)
        field = 'goat.vessel.offsetfracx'
        call ExtractOptionValueReal0D(fid, field, options%offsetfracx)
        field = 'goat.vessel.offsetfracy'
        call ExtractOptionValueReal0D(fid, field, options%offsetfracy)

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