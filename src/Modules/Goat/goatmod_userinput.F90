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
    use mod_global_environment, only: solps, solps_inputfilepath, &
        solps_writefilepath, solps_gridfilepath, solps_magneticfieldfilepath, &
        solps_structurefilepath

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
        ! - RBtor       product of torodial magnetic field and major 
        !               radius (assumed constants here)

        
        character(:), allocatable   :: readmeth
        character(:), allocatable   :: filepath
        character(:), allocatable   :: interpmeth 

        integer(I8)                 :: interpC, interpM
        real(R8)                    :: RBtor

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
        ! - xrange, yrange: only used for interpolant based methods. 
        ! Sets range over which interpolant can be evaluated (range is 
        ! determined by max and min value in array, may be overwritten))

        character(:), allocatable       :: readmeth
        character(:), allocatable       :: filepath

        real(R8)                        :: maxdist
        integer(I8)                     :: refine
        integer(I8), allocatable        :: TP(:), TPind(:), exclude(:)

        ! Vessel representation options
        character(:), allocatable       :: shapemeth 
        integer(I8)                     :: resx, resy, interpC, interpM 
        real(R8)                        :: offsetfracx, offsetfracy
        real(R8), allocatable           :: xrange(:), yrange(:)

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
    !                          GRID GENERATION                         !
    !------------------------------------------------------------------!

    ! Options for topological mesh generation
    type, extends(OptionsUDT) :: TopomeshOptionsUDT 

        ! Description
        !============
        ! Options for the topological mesh construction. The 
        ! following fields are set:
        ! - fresx, fresy    : resolution for 2D tracing of contours and
        !                   extrema of the field
        ! - fdonewton       : option to refine extrema with newton 
        !                   solver (may not always converge!)
        ! - v(...)          : same options but for vessel 
        ! - ffieldtol       : tolerance on field value of extrema (if 
        !                   difference is below tolerance, two extrema 
        !                   are considered to have exactly the same field
        !                   value)
        ! - doadaptations       : general switch to apply or not apply
        !                       adaptations to the basic topological 
        !                       mesh (all things like core/PF boundaries, 
        !                       optional region removal, tube merging, ...)
        ! - addcoreboundaries   : adds additional core boundaries
        ! - addPFboundaries     : adds additional private flux like 
        !                       boundaries (at tangency points of which 
        !                       contour doesn't enter vessel)
        ! - removecoreregions   : remove the innermost core regions
        ! - coreboundariesfrac  : fraction in field value between core 
        !                       and connecting X-point (1: tangency point, 0: 
        !                       other PF boundary, e.g. separatrix)
        ! - PFboundariesfrac    : fraction in field value between outer 
        !                       PF boundary and tangency point
        ! - removewidegridregions: remove all regions that are not next    
        !                           to a separatrix 
        ! - npmin, npmax, dl    : minimal and maximal number of points
        !                       of contours (when doing coarsening) and
        !                       desired uniform edge length
        ! - readexistingTM:     read in an existing topomesh file, 
        !                       for which the full path is defined in  
        !                       TMfilepath
        ! - mergetangencypointtubes     merge tubes that are too small 
        !                       and that have tangency point tubes 
        !                       as neighbours. 'too small' is based on 
        !                       the (absolute) difference in flux values
        !                       of the tube's radial face vertices
        ! - dpsimintangencypointtubes   minimal delta psi for tangency 
        !                       point tubes (if below, we attempt to 
        !                       merge)

        integer(I8)             :: fresx, fresy, vresx, vresy, npmin, &
            npmax
        logical                 :: addcoreboundaries, removecoreregions, &
            fdonewton, vdonewton, removewidegridregions, addPFboundaries, &
            readexistingTM, removenoncoreregions, mergetangencypointtubes, &
            doadaptations
        real(R8)                :: coreboundariesfrac, ffieldtol, dl, &
            PFboundariesfrac, dpsimintangencypointtubes
        character(:), allocatable   :: TMfilepath
    contains 

        procedure :: Read           => ReadTopomeshOptions
        procedure :: SetDefaults    => SetDefaultTopomeshOptions

    end type

    ! Options for grid generation
    type, extends(OptionsUDT) :: GGOptionsUDT 

        ! Description
        !============
        ! Options for the grid generator, including vertex distribution 
        ! options (poloidal and radial), flux surface removal options, 
        ! boundary triangle removal options, and small face removal 
        ! options. The following fields are present: 

        ! General grid generation
        ! - verbosity       the higher, the more information is 
        !                   printed out (default: 1, 0 suppresses all)
        ! - ggmethod        method how to construct the grid. Can be 
        !                   'independent' (treating each flux surface
        !                    independently) or 'orthogonal' (
        !                   yields largely orthogonal grid, but needs 
        !                   dependency between mesh cells)
        ! - cellconstructionmethod: method how to determine grid cells
        !                   starting from given vertex distribution. 
        !                   'quads_triangles' is recommended one
        ! - TMcellgriddingorder:    order in which to grid the 
        !                   topological mesh cells. Can be 'independent'
        !                   (then it is  pretty much random) or 
        !                   'sequential' -> that one is the recommended
        !                   one. Here, initial distributions are 
        !                   propagated through
        ! 
        ! Vertex distribution, poloidal direction:
        ! - vdptype:        'uniform' for uniform distribution, 
        !                   'densitybased' for non-uniform distribution 
        !                   based on predefined node density distribution
        !                   function 
        ! - vdpdfacelength: length of face for 'uniform' distribution 
        ! - vdpddecaylength:    decay length from vessel for 
        !                       'densitybased' option
        ! - vdpdx, y, d, val:   parameters to add attraction points for
        !                       densitybased option: x, y are 
        !                       coordinates, d is a decay length, and 
        !                       val is the desired density at the point
        ! - vdpddensityatvessel:    desired density at the vessel boundary
        ! - vdpddensityatinf:       desired density far from the vessel 
        !                           boundary 

        ! Vertex distribution, radial direction:
        ! largely the same as poloidal one (but then vdr instead of vdp)
        ! but some differences
        ! - vdrtype:        (same as vdptype)
        ! - vdrdfieldwidth: desired field with for uniform distribution
        ! - vdrddecaylength:    decay length parameter 
        ! - vdrddensityatseparatrix:    desired density at separatrix
        ! - vdrddensityatinf:           density far from separatrix
 
        ! Options for extending flux tubes with vessel parts (so-called 
        ! 'cut cells')
        ! - extendtptubes:          extend tubes at the tangency point
        !                           side to include (part of) the 
        !                           vessel from both side of type 1
        !                           tangency points

        ! Options for flux surface removal
        ! - removefluxsurfaces:     switch to remove or not
        ! - remfspsitol:            absolute tolerance for difference in 
        !                           flux value between flux surfaces 
        !                           (if below, one is removed)
        ! - remfspsirattol:         maximal ratio between flux surface 
        !                           difference (and inverse), if above/
        !                           below: removed

        ! Options for boundary triangle removal
        ! - removenarrowboundarytriangles:  switch
        ! - rembndtriacriterion:    criterion for removal, typically 
        !                           'angle' (if too small, remove triangle)
        ! - rmbndtriaminangle:      minimal angle [rad] (too small - removed)
        !                           input is in degrees!
        
        ! Options for small face removal
        ! - removefaces             switch 
        ! - remfacescriterion       criterion to remove faces. Typically
        !                           'facelength_radial_bnd' to remove
        !                           too narrow radial faces at the 
        !                           boundary 
        ! - remfacesminlength       minimal length of these faces  [m] 

        ! Options for refinement
        ! - refmeth:                method for refinement. Can be 'no' 
        !                           (no additional refinement), 
        !                           'lengthbased' (ref based on min and
        !                           max length distributions) 
        
        ! (poloidal) Refinement options for lengthbased option (this is currently 
        ! based on exponential decay functions defined in points):
        ! - refLBlmininf    minimal length at infinity
        ! - refLBLmaxinf    maximal length at infinity
        ! - refLBdoxp:      refine near x-points (x-points are added)
        ! - refLBLminxp     minimal length on x-point
        ! - refLBLmaxxp     maximal length on x-point
        ! - refLBdecaylengthxp  decaylength on x-point (larger - wider influence)
        ! - refLBdovessel   refine near (specified) vessel vertices
        ! - refLBstructureIDs   structure IDs of vessel structures to include
        ! - refLBvertIDs:   vessel vertex IDs of vessel vertices to include, 
        !                   as specified in structure.dat file (some 
        !                   may not be included if deleted during vessel 
        !                   polygon construction)
        ! - refLBLminstructure  minimal length on structure i (etc, 
        !                   similar for vertices)

        ! - refBLdotarget   do BL refinement at targets
        ! - refBLdovessel   do BL refinement at far vessel boundaries
        ! - refBLnctarget   number of desired boundary layer cells at 
        !                   the target (similar for vessel)
        ! - refBLdltarget   desired lengths for these cells 

        ! (Radial) refinement options for lengthbased refiner:
        !   mostly the same refinement options as the poloidal direction,
        !   but names now have 'rad' in front. 
        ! - radrefLBlmininf    minimal length at infinity (in [m]!)
        ! - radrefLBLmaxinf    maximal length at infinity
        ! - radrefLBdosp:      refine near strike points (x-points are added)
        ! - radrefLBLminsp     minimal length on strike point
        ! - radrefLBLmaxsp     maximal length on strike point
        ! - radrefLBdecaylengthsp  decaylength on strike point (larger - wider influence)
        
        ! - radrefBLdosp    do BL refinement at strike points
        ! - radrefBLncsp   number of desired boundary layer cells at 
        !                   the strike point 
        ! - radrefBLdlsp   desired lengths for these cells

        
        logical                     :: removefluxsurfaces, &
            removenarrowboundarytriangles, removefaces, refLBdoxp, &
            refLBdovessel, vdpdincludexp, coarsencontours, refBLdotarget, &
            refBLdovessel, readexistingrefdata, radrefBLdosp, radrefLBdosp, &
            extendtptubes
        integer(I8)                 :: gcresx, gcresy, &
            verbosity, orthtracernsteps, refBLnctarget, refBLncvessel, &
            radrefBLncsp
        integer(I8), allocatable, dimension(:)  :: refLBstructureIDs, &
            refLBvertIDs
        real(R8)                    :: vdpdfacelength, vdpddecaylengthplf, &
            vdpddecaylengthxp, vdpddensityatvessel, vdpddensityatxp, &
            vdpddensityatinf, vdrdfieldwidth, &
            vdrddecaylength, vdrddensityatseparatrix, vdrddensityatinf, &
            remfspsitol, remfspsirattol, rembndtriaminangle, &
            remfacesminlength, refLBLmininf, refLBLmaxinf, refLBLminxp, &
            refLBLmaxxp, refLBdecaylengthxp, orthtracerstep, &
            radrefLBLmininf, radrefLBLmaxinf, radrefLBLminsp, &
            radrefLBLmaxsp, radrefLBdecaylengthsp
        real(R8), allocatable, dimension(:)     :: vdpdx, vdpdy, vdpdd, &
            vdpdval, refLBLminstructure, refLBLminvert, refLBLmaxstructure, &
            refLBLmaxvert, refLBdecaylengthstructure, refLBdecaylengthvert, &
            refBLdltarget, refBLdlvessel, radrefBLdlsp
        character(:), allocatable   :: vdptype, vdpdtype, vdrtype, &
            vdrdtype, rembndtriacriterion, remfacescriterion, ggmethod, &
            cellconstructionmethod, TMcellgriddingorder, refmeth, vdpplftype, &
            refdatafile, radrefmeth
    contains 

        procedure :: Read           => ReadGGOptions
        procedure :: SetDefaults    => SetDefaultGGOptions

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

    ! Goat options routines
    subroutine SetDefaultGoatOptions(options)

        ! Declare variables
        !==================
        ! Arguments
        class(GoatoptionsUDT)       :: options        
        
        ! Input file
        if (solps) then 
            ! SOLPS defaults
            options%inputfilepath    = solps_inputfilepath
        else
            options%inputfilepath    = './GOAToptions.dat'
        end if 

        ! General
        options%debug           = .false. 
        options%meth            = 'GD'
        options%gdinputfilepath = './GOAToptions.dat'

        ! Specify input filenames
        if (solps) then 
            ! SOLPS defaults
            options%gridfilepath            = solps_gridfilepath
            options%structurefilepath       = solps_structurefilepath
            options%magneticfieldfilepath   = solps_magneticfieldfilepath
        else
            options%gridfilepath            = './traduit.out.b2us'
            options%structurefilepath       = './structure.dat'
            options%magneticfieldfilepath   = './rzpsi.dat'
        end if 

        ! Output options
        if (solps) then 
            ! SOLPS defaults
            options%writefilepath   = solps_writefilepath
        else
            options%writefilepath       = 'traduit.out.b2us_smoothed'
        end if 
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

        options%RBtor                   = 0

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
        allocate(options%xrange(0), options%yrange(0))
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

    ! Topomesh
    subroutine SetDefaultTopomeshOptions(options)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshOptionsUDT)   :: options

        ! Set defaults
        !=============
        ! I/O
        options%readexistingTM = .false. 
        options%TMfilepath = 'topomesh.dat'
        
        ! Contouring (field)
        options%fresx = 100
        options%fresy = 100

        ! Contouring (vessel)
        options%vresx = 100
        options%vresy = 100

        ! Contouring (general)
        options%npmin = 10
        options%npmax = 10000
        options%dl = 1e-3

        ! Refining options for extrema (field)
        options%fdonewton = .true. 
        options%ffieldtol = 1e-8 
        
        ! Refining options for extrema (vessel)
        options%vdonewton = .true.

        ! Additional options
        options%doadaptations               = .true.
        options%addcoreboundaries           = .true. 
        options%coreboundariesfrac          = 0.2
        options%addPFboundaries             = .true. 
        options%PFboundariesfrac            = 0.2
        options%removecoreregions           = .true. 
        options%removewidegridregions       = .true. 
        options%removenoncoreregions        = .false.
        options%mergetangencypointtubes     = .false.
        options%dpsimintangencypointtubes   = 100_R8 ! some absurd large value

    end subroutine 

    ! Grid generation
    subroutine SetDefaultGGOptions(options)

        ! Modules
        !========
        use mod_constants, only : pi_R8

        ! Declare variables
        !==================
        ! Arguments
        class(GGOptionsUDT)   :: options

        ! Set defaults
        !=============
        ! Contouring options in grid generator
        options%gcresx = 100 
        options%gcresy = 100
        options%coarsencontours = .false.

        ! Options for poloidal vertex distribution
        options%vdptype             = 'densitybased'
        options%vdpplftype          = 'target'
        options%vdpdfacelength      = 4e-2 
        options%vdpddecaylengthplf  = 0.05
        options%vdpddecaylengthxp   = 0.05
        options%vdpdincludexp       = .false.
        allocate(options%vdpdx(0), options%vdpdy(0), options%vdpdd(0), &
            options%vdpdval(0))
        options%vdpddensityatvessel = 250.0_R8
        options%vdpddensityatxp     = 250.0_R8
        options%vdpddensityatinf    = 10.0_R8

        ! Grid generation approach
        options%verbosity           = 1
        options%ggmethod            = 'independent'
        options%cellconstructionmethod  = 'quads_triangles'
        options%TMcellgriddingorder = 'sequential'
        options%readexistingrefdata = .false. 
        options%refdatafile         = './output/refdataTM.dat'

        ! Poloidal refinement options ('lengthbased' refinement options only)
        options%refmeth         = 'no'      
        options%refLBdoxp       = .true. 
        options%refLBdovessel   = .false. 
        options%refLBLmininf    = 0.0_R8
        options%refLBLmaxinf    = 100 ! some absurd big number
        options%refLBLminxp     = 0.0_R8
        options%refLBLmaxxp     = 100 ! some absurd big number  
        options%refLBdecaylengthxp = 0.1 
        allocate(options%refLBdecaylengthstructure(0), &
            options%refLBdecaylengthvert(0), options%refLBstructureIDs(0), &
            options%refLBvertIDs(0), options%refLBLminstructure(0), &
            options%refLBLmaxstructure(0))
        
        ! Radial refinement options
        options%radrefmeth         = 'no'      
        options%radrefLBdosp       = .true. 
        options%radrefLBLmininf    = 0.0_R8
        options%radrefLBLmaxinf    = 100_R8 ! some absurd big number
        options%radrefLBLminsp     = 0.0_R8
        options%radrefLBLmaxsp     = 100_R8 ! some absurd big number  
        options%radrefLBdecaylengthsp = 0.1_R8 

        ! Poloidal boundary layer options
        options%refBLdotarget   = .false. 
        options%refBLdovessel   = .false. 
        options%refBLnctarget   = 0
        options%refBLncvessel   = 0
        allocate(options%refBLdltarget(options%refBLnctarget), &
            options%refBLdlvessel(options%refBLncvessel))
        options%refBLdltarget   = 0.001_R8 ! in m 
        options%refBLdlvessel   = 0.01_R8 ! in m 

        ! Radial boundary layer options
        options%radrefBLdosp    = .false. 
        options%radrefBLncsp    = 0
        allocate(options%radrefBLdlsp(options%radrefBLncsp))
        options%radrefBLdlsp       = 0.001_R8 ! in m 

        ! Streamline tracer options
        options%orthtracerstep = 0.5
        options%orthtracernsteps = 2000

        ! Options for radial vertex distribution
        options%vdrtype             = 'densitybased'
        options%vdrdfieldwidth      = 4e-3
        options%vdrddecaylength     = 0.005
        options%vdrddensityatseparatrix     = 2500.0_R8
        options%vdrddensityatinf            = 250.0_R8

        ! Options for flux tube extensions
        options%extendtptubes       = .true. 

        ! Options for flux surface removal 
        options%removefluxsurfaces = .true.
        options%remfspsitol = 1e-4 
        options%remfspsirattol = 1e-1 
        
        ! Options for boundary triangle removal
        options%removenarrowboundarytriangles = .true. 
        options%rembndtriacriterion = 'angle' 
        options%rembndtriaminangle = 15.0_R8*pi_R8/180.0_R8

        ! Options for small face removal
        options%removefaces = .false. 
        options%remfacescriterion = 'facelength_radial_bnd'
        options%remfacesminlength = 2e-3        

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
        logical                         :: reachedeof, throwerror

        ! Initialize
        !===========
        ! Variables
        reachedeof = .false. 
        throwerror = .false. 

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
        field = 'goat.OMPr'
        call ExtractOptionValueReal1D(fid, field, &
            options%OMP_r)
        field = 'goat.OMPz'
        call ExtractOptionValueReal1D(fid, field, &
            options%OMP_z)
        field = 'goat.IMPr'
        call ExtractOptionValueReal1D(fid, field, &
            options%IMP_r)
        field = 'goat.IMPz'
        call ExtractOptionValueReal1D(fid, field, &
            options%IMP_z)

        ! Checks
        !=======
        ! Check for name clashes when using solps
        if (solps) then 
            ! Check input file names and assert they are equal to 
            ! assumed solps filenames
            if (options%gridfilepath /= solps_gridfilepath) then 
                print *, 'ReadGoatOptions: set goat.gridfilepath value ' // &
                    'to "' // solps_gridfilepath // '" in GOAToptions.dat'
                throwerror = .true. 
            end if 
            if (options%magneticfieldfilepath /= solps_magneticfieldfilepath) then 
                print *, 'ReadGoatOptions: set goat.magneticfieldfilepath value ' // &
                    'to "' // solps_magneticfieldfilepath // '" in GOAToptions.dat'
                throwerror = .true. 
            end if 
            if (options%structurefilepath /= solps_structurefilepath) then 
                print *, 'ReadGoatOptions: set goat.structurefilepath value ' // &
                    'to "' // solps_structurefilepath // '" in GOAToptions.dat'
                throwerror = .true. 
            end if 
            if (options%writefilepath /= solps_writefilepath) then 
                print *, 'ReadGoatOptions: set goat.writefilepath value ' // &
                    'to "' // solps_writefilepath // '" in GOAToptions.dat'
                throwerror = .true. 
            end if 

            ! Throw error
            if (throwerror) then 
                call gdErrorHandler('ReadGoatOptions: exiting due to ' // &
                    'wrong I/O filenames for SOLPS (see messages above)')
            end if

        end if 

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

        ! RBtor
        field  = 'goat.mf.RBtor'
        call ExtractOptionValueReal0D(fid, field, options%RBtor)

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
        field = 'goat.vessel.xrange'
        call ExtractOptionValueReal1D(fid, field, options%xrange)
        field = 'goat.vessel.yrange'
        call ExtractOptionValueReal1D(fid, field, options%yrange)

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

    ! Topomesh options reader
    subroutine ReadTopomeshOptions(options)

        ! Description
        !============
        ! Read in Topomesh options from file. It is assumed that the 
        ! filepath has been set correctly. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshOptionsUDT)       :: options 

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
            print *, 'ReadTopomeshOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadTopomeshOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! I/O
        field = 'gg.tm.readexistingTM'
        call ExtractOptionValueLogical0D(fid, field, options%readexistingTM)
        field = 'gg.tm.TMfilepath'
        call ExtractOptionValueCharacter(fid, field, options%TMfilepath)

        ! Resolution 
        field = 'gg.tm.field.resx'
        call ExtractOptionValueInteger0D(fid, field, options%fresx)
        field = 'gg.tm.field.resy'
        call ExtractOptionValueInteger0D(fid, field, options%fresy)
        field = 'gg.tm.vessel.resx'
        call ExtractOptionValueInteger0D(fid, field, options%vresx)
        field = 'gg.tm.vessel.resy'
        call ExtractOptionValueInteger0D(fid, field, options%vresy)
        field = 'gg.tm.contour.npmin'
        call ExtractOptionValueInteger0D(fid, field, options%npmin)
        field = 'gg.tm.contour.npmax'
        call ExtractOptionValueInteger0D(fid, field, options%npmax)
        field = 'gg.tm.contour.dl'
        call ExtractOptionValueReal0D(fid, field, options%dl)

        ! Refinement
        field = 'gg.tm.field.donewton'
        call ExtractOptionValueLogical0D(fid, field, options%fdonewton)
        field = 'gg.tm.vessel.donewton'
        call ExtractOptionValueLogical0D(fid, field, options%vdonewton)
        field = 'gg.tm.fieldtol'
        call ExtractOptionValueReal0D(fid, field, options%ffieldtol)

        ! Additional options
        field = 'gg.tm.doadaptations'
        call ExtractOptionValueLogical0D(fid, field, options%doadaptations)
        field = 'gg.tm.addcoreboundaries'
        call ExtractOptionValueLogical0D(fid, field, options%addcoreboundaries)
        field = 'gg.tm.addPFboundaries'
        call ExtractOptionValueLogical0D(fid, field, options%addPFboundaries)
        field = 'gg.tm.removecoreregions'
        call ExtractOptionValueLogical0D(fid, field, options%removecoreregions)
        field = 'gg.tm.removewidegridregions'
        call ExtractOptionValueLogical0D(fid, field, options%removewidegridregions)
        field = 'gg.tm.cbnd.frac'
        call ExtractOptionValueReal0D(fid, field, options%coreboundariesfrac)
        field = 'gg.tm.PFbnd.frac'
        call ExtractOptionValueReal0D(fid, field, options%PFboundariesfrac)
        field = 'gg.tm.removenoncoreregions'
        call ExtractOptionValueLogical0D(fid, field, options%removenoncoreregions)
        field = 'gg.tm.mergetangencypointtubes'
        call ExtractOptionValueLogical0D(fid, field, options%mergetangencypointtubes)
        field = 'gg.tm.dpsimintangencypointtubes'
        call ExtractOptionValueReal0D(fid, field, options%dpsimintangencypointtubes)

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)

    end subroutine

    ! Grid generator options reader
    subroutine ReadGGOptions(options)

        ! Description
        !============
        ! Read in GG options from file. It is assumed that the 
        ! filepath has been set correctly. 

        ! Declare variables
        !==================
        ! Arguments
        class(GGOptionsUDT)             :: options 

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
            print *, 'ReadGGOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadGGOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! General options
        field  = 'gg.verbosity'
        call ExtractOptionValueInteger0D(fid, field, options%verbosity)
        field  = 'gg.vertexconstructionmethod'
        call ExtractOptionValueCharacter(fid, field, options%ggmethod)
        field  = 'gg.cellconstructionmethod'
        call ExtractOptionValueCharacter(fid, field, options%cellconstructionmethod)
        field = 'gg.TMcellgriddingorder'
        call ExtractOptionValueCharacter(fid, field, options%TMcellgriddingorder)

        ! Refinement options (general)
        field = 'gg.ref.meth'
        call ExtractOptionValueCharacter(fid, field, options%refmeth) 
        field = 'gg.radref.meth'
        call ExtractOptionValueCharacter(fid, field, options%radrefmeth) 
        field = 'gg.ref.readexistingrefdata'
        call ExtractOptionValueLogical0D(fid, field, options%readexistingrefdata) 
        field = 'gg.ref.refdatafile'
        call ExtractOptionValueCharacter(fid, field, options%refdatafile) 

        ! Length-based refinement options (poloidal)
        field  = 'gg.ref.LB.doxp'
        call ExtractOptionValueLogical0D(fid, field, options%refLBdoxp)
        field  = 'gg.ref.LB.dovessel'
        call ExtractOptionValueLogical0D(fid, field, options%refLBdovessel)
        field  = 'gg.ref.LB.Lmininf'
        call ExtractOptionValueReal0D(fid, field, options%refLBLmininf)
        field  = 'gg.ref.LB.Lmaxinf'
        call ExtractOptionValueReal0D(fid, field, options%refLBLmaxinf)
        
        field  = 'gg.ref.LB.Lminxp'
        call ExtractOptionValueReal0D(fid, field, options%refLBLminxp)
        field  = 'gg.ref.LB.Lmaxxp'
        call ExtractOptionValueReal0D(fid, field, options%refLBLmaxxp)
        field  = 'gg.ref.LB.decaylengthxp'
        call ExtractOptionValueReal0D(fid, field, options%refLBdecaylengthxp)

        field  = 'gg.ref.LB.Lminstructure'
        call ExtractOptionValueReal1D(fid, field, options%refLBLminstructure)
        field  = 'gg.ref.LB.Lmaxstructure'
        call ExtractOptionValueReal1D(fid, field, options%refLBLmaxstructure)
        field  = 'gg.ref.LB.decaylengthstructure'
        call ExtractOptionValueReal1D(fid, field, options%refLBdecaylengthstructure)
        field  = 'gg.ref.LB.structureIDs'   
        call ExtractOptionValueInteger1D(fid, field, options%refLBstructureIDs)
        field  = 'gg.ref.LB.vertIDs'   
        call ExtractOptionValueInteger1D(fid, field, options%refLBvertIDs)

        ! Refinement options (radial)
        field  = 'gg.radref.LB.dosp'
        call ExtractOptionValueLogical0D(fid, field, options%radrefLBdosp)
        field  = 'gg.radref.LB.Lmininf'
        call ExtractOptionValueReal0D(fid, field, options%radrefLBLmininf)
        field  = 'gg.radref.LB.Lmaxinf'
        call ExtractOptionValueReal0D(fid, field, options%radrefLBLmaxinf)
        
        field  = 'gg.radref.LB.Lminsp'
        call ExtractOptionValueReal0D(fid, field, options%radrefLBLminsp)
        field  = 'gg.radref.LB.Lmaxsp'
        call ExtractOptionValueReal0D(fid, field, options%radrefLBLmaxsp)
        field  = 'gg.radref.LB.decaylengthsp'
        call ExtractOptionValueReal0D(fid, field, options%radrefLBdecaylengthsp)

        ! Boundary layer options (only for length-based ref, poloidal)
        field = 'gg.ref.BL.dotarget'
        call ExtractOptionValueLogical0D(fid, field, options%refBLdotarget)
        field = 'gg.ref.BL.dovessel'
        call ExtractOptionValueLogical0D(fid, field, options%refBLdovessel)
        field = 'gg.ref.BL.nctarget'
        call ExtractOptionValueInteger0D(fid, field, options%refBLnctarget)
        field = 'gg.ref.BL.ncvessel'
        call ExtractOptionValueInteger0D(fid, field, options%refBLncvessel)
        field = 'gg.ref.BL.dltarget'
        call ExtractOptionValueReal1D(fid, field, options%refBLdltarget)
        field = 'gg.ref.BL.dlvessel'
        call ExtractOptionValueReal1D(fid, field, options%refBLdlvessel)

        ! Boundary layer options (only for length-based ref, radial)
        field = 'gg.radref.BL.dosp'
        call ExtractOptionValueLogical0D(fid, field, options%radrefBLdosp)
        field = 'gg.radref.BL.ncsp'
        call ExtractOptionValueInteger0D(fid, field, options%radrefBLncsp)
        field = 'gg.radref.BL.dlsp'
        call ExtractOptionValueReal1D(fid, field, options%radrefBLdlsp)

        ! Contouring options in grid generator
        field = 'gg.vd.contouring.resx'
        call ExtractOptionValueInteger0D(fid, field, options%gcresx)
        field = 'gg.vd.contouring.resy'
        call ExtractOptionValueInteger0D(fid, field, options%gcresy)
        field = 'gg.vd.contouring.coarsen'
        call ExtractOptionValueLogical0D(fid, field, options%coarsencontours)

        ! Orthogonal line tracer options
        field = 'gg.vd.orthlinetracing.step'
        call ExtractOptionValueReal0D(fid, field, options%orthtracerstep)
        field = 'gg.vd.orthlinetracing.nsteps'
        call ExtractOptionValueInteger0D(fid, field, options%orthtracernsteps)

        ! Options for poloidal vertex distribution
        field = 'gg.vd.pd.type'
        call ExtractOptionValueCharacter(fid, field, options%vdptype)
        field = 'gg.vd.pd.distribution.plftype'
        call ExtractOptionValueCharacter(fid, field, options%vdpplftype)
        field = 'gg.vd.pd.distribution.facelength'
        call ExtractOptionValueReal0D(fid, field, options%vdpdfacelength)
        field = 'gg.vd.pd.distribution.decaylengthplf'
        call ExtractOptionValueReal0D(fid, field, options%vdpddecaylengthplf)
        field = 'gg.vd.pd.distribution.decaylengthxp'
        call ExtractOptionValueReal0D(fid, field, options%vdpddecaylengthxp)
        field = 'gg.vd.pd.distribution.vdpdincludexp'
        call ExtractOptionValueLogical0D(fid, field, options%vdpdincludexp)
        field = 'gg.vd.pd.distribution.points.x'
        call ExtractOptionValueReal1D(fid, field, options%vdpdx)
        field = 'gg.vd.pd.distribution.points.y'
        call ExtractOptionValueReal1D(fid, field, options%vdpdy)
        field = 'gg.vd.pd.distribution.points.d'
        call ExtractOptionValueReal1D(fid, field, options%vdpdd)
        field = 'gg.vd.pd.distribution.points.val'
        call ExtractOptionValueReal1D(fid, field, options%vdpdval)
        field = 'gg.vd.pd.distribution.densityatvessel'
        call ExtractOptionValueReal0D(fid, field, options%vdpddensityatvessel)
        field = 'gg.vd.pd.distribution.densityatxp'
        call ExtractOptionValueReal0D(fid, field, options%vdpddensityatxp)
        field = 'gg.vd.pd.distribution.densityatinf'
        call ExtractOptionValueReal0D(fid, field, options%vdpddensityatinf)

        ! Options for radial vertex distribution
        field = 'gg.vd.rd.type'
        call ExtractOptionValueCharacter(fid, field, options%vdrtype)
        field = 'gg.vd.rd.distribution.fieldwidth'
        call ExtractOptionValueReal0D(fid, field, options%vdrdfieldwidth)
        field = 'gg.vd.rd.distribution.decaylength'
        call ExtractOptionValueReal0D(fid, field, options%vdrddecaylength)
        field = 'gg.vd.rd.distribution.densityatseparatrix'
        call ExtractOptionValueReal0D(fid, field, options%vdrddensityatseparatrix)
        field = 'gg.vd.rd.distribution.densityatinf'
        call ExtractOptionValueReal0D(fid, field, options%vdrddensityatinf)

        ! Options for extending flux tubes
        field = 'gg.adap.extendtptubes'
        call ExtractOptionValueLogical0D(fid, field, options%extendtptubes)

        ! Options for flux surface removal 
        field = 'gg.vd.removefluxsurfaces'
        call ExtractOptionValueLogical0D(fid, field, options%removefluxsurfaces)
        field = 'gg.vd.rfs.mark.psitol'
        call ExtractOptionValueReal0D(fid, field, options%remfspsitol)
        field = 'gg.vd.rfs.mark.psirattol'
        call ExtractOptionValueReal0D(fid, field, options%remfspsirattol)
        
        ! Options for boundary triangle removal
        field = 'gg.vd.removenarrowboundarytriangles'
        call ExtractOptionValueLogical0D(fid, field, options%removenarrowboundarytriangles)
        field = 'gg.vd.rnbt.mark.criterion'
        call ExtractOptionValueCharacter(fid, field, options%rembndtriacriterion)
        field = 'gg.vd.rnbt.mark.minangle'
        call ExtractOptionValueReal0D(fid, field, options%rembndtriaminangle)

        ! Options for small face removal
        field = 'gg.vd.removefaces'
        call ExtractOptionValueLogical0D(fid, field, options%removefaces)
        field = 'gg.vd.rf.criterion'
        call ExtractOptionValueCharacter(fid, field, options%remfacescriterion)
        field = 'gg.vd.rf.minlength'
        call ExtractOptionValueReal0D(fid, field, options%remfacesminlength)       

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)

    end subroutine


end module goatmod_userinput