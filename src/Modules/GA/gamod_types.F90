!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains all types that are only used in the grid adaptation module


module gamod_types

    ! Initialize
    !===========
    ! Load modules
    use mod_sort   
    use mod_precision
    use mod_dynamicarrays    
    use goatmod_types 
    use DistributionFunction
    !use goatmod_userinput
    !use gdmod_utility_optimization


    ! The usual
    implicit none
    save
    public

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!   

    type :: QualityMetric

        ! Type that included several metrics
        real(R8), allocatable :: fcBias(:)
        real(R8), allocatable :: fcqalfc(:)
        real(R8), allocatable :: fcS(:)
        real(R8), allocatable :: cvS(:)
        real(R8), allocatable :: cvAR(:)
        real(R8), allocatable :: h_pol(:)
        real(R8), allocatable :: h_rad(:)
        real(R8), allocatable :: h_rad_psi(:)
        integer(I8) :: nCv

    contains
    
        procedure :: ComputeQM

    end type 
    
    ! Vertex structure
    type GAVertexUDT

        ! Description
        !============
        ! Fields:
        ! - x, y            : coordinates [m]
        ! - BV              : logical index to indicate if vertex lies
        !                   on boundary
        ! - fieldlineID     : ID of the field line the vertex lies on
        ! - ntot            : total number of vertices  
        ! - faceP           : ntot-by-2 array containing in the first 
        !                   column the starting index and in the second
        !                   column the number of faces of the cell 
        !                   (for querying the faces of the cell
        !                   stored in vert%facelist)            
        ! - facelist        : list of vertex faces, to be queried as 
        !                   vert%facelist(vert%faceP(i,1):
        !                   vert%faceP(i,1)+vert%faceP(i,2)-1) 
        ! - nfacelist       : length of facelist   
        ! - cellP           : similar to faceP, but for vertex cells
        ! - celllist        : similar to facelist, but for vertex cells
        ! - ncellist        : length of celllist    
        ! - neiglist        : list of neighbouring vertices (nneiglist-
        !                   by-1)
        ! - nneiglist       : dimension of neiglist (scalar)
        ! - neigP           : ntot-by-2 array, analogous to faceP and 
        !                   cellP, but for neiglist.
        ! - psi             : psi values at vertex locations    
        ! - bx, by          : magnetic field vector at vertex locations  
        ! - ffbz            : ???

        ! Coordinates
        class(RealDynamicArrayUDT), allocatable        :: x, y

        ! Logicals and indices
        !logical, allocatable                          :: BV(:)
        class(IntegerDynamicArrayUDT), allocatable        :: fieldlineID 
        integer(I8)                                   :: ntot = 0

        !type(IntegerDynamicArrayUDT), allocatable     :: faceP(:,:)
        !type(IntegerDynamicArrayUDT), allocatable     :: face(:)
        !integer(I8)                                   :: nface = 0

        !type(IntegerDynamicArrayUDT), allocatable     :: cellP(:,:)
        !type(IntegerDynamicArrayUDT), allocatable     :: cell(:)
        !integer(I8)                                   :: ncell = 0

        !type(IntegerDynamicArrayUDT), allocatable     :: neigP(:,:)
        !type(IntegerDynamicArrayUDT), allocatable     :: neig(:)
        !integer(I8)                                   :: nneig = 0

        ! Other data
        class(RealDynamicArrayUDT), allocatable :: bx, &
            by, psi, ffbz 
    contains

        ! Initializer
        procedure :: Initialize     => InitializeGAVertex

    end type 

    ! Face structure
    type GAFaceUDT
        ! Description
        !============
        ! Fields:
        ! - ntot            : total number of faces
        ! - vert            : set of (two) vertices belonging to that 
        !                   face
        ! - cell            : list of cell neighbours (max 2)
        ! - cellP           : pointer for list of cell neighbours
        ! - neig            : cell neighbours of face
        ! - BF              : logical index that is true if the face
        !                   is a boundary face
        ! - label           : face labels
        ! - TMfacelabel     : face label corresponding to topomesh ID (only when GG is used)
        ! - reg             : face regions
        ! - aligned         : integer that is 1 if the face is aligned


        ! Logicals and indices
        class(IntegerDynamicArrayUDT), allocatable         :: vert1
        class(IntegerDynamicArrayUDT), allocatable         :: vert2

        class(IntegerDynamicArrayUDT), allocatable         :: &
            label, reg, aligned

        !logical :: aligned(:)
        !type(IntegerDynamicArrayUDT), allocatable         :: cell(:)    
        !type(IntegerDynamicArrayUDT), allocatable         :: cellP(:,:)
        !integer(I8)                                       :: ncell = 0

        integer(I8)                                       :: ntot = 0
        !type(IntegerDynamicArrayUDT), allocatable         :: BF(:)
        !logical, allocatable                                 :: BF(:)

    contains

        ! Initialize
        procedure :: Initialize     => InitializeGAFace
        procedure :: ChainFaces

    end type

    ! Cell structure
    type GACellUDT

        ! Description
        !============
        ! Fields:   
        ! - ntot            : total number of cells
        ! - vertP           : ntot-by-2 array containing in the first 
        !                   column the starting index and in the second
        !                   column the number of vertices of the cell 
        !                   (for querying the vertices of the cell
        !                   stored in cells%vert)            
        ! - vert            : list of cell vertices, to be queried as 
        !                   cells%vert(cells%vertP(i,1):
        !                   cells%vertP(i,1)+cells%vertP(i,2)-1) 
        ! - nvert          : length of vert
        ! - faceP           : similar to vertP, but for faces
        ! - face            : similar to vert, but for faces
        ! - nface          : similar to nvert, but for faces
        ! - GC              : ntot-by-1 logical vector indicating if the
        !                   cell is a guard cell


        ! Logicals and indices
        class(IntegerDynamicArrayUDT), allocatable        :: vertP1
        class(IntegerDynamicArrayUDT), allocatable        :: vertP2
        class(IntegerDynamicArrayUDT), allocatable        :: vert
        integer(I8)                                       :: nvert = 0

        class(IntegerDynamicArrayUDT), allocatable        :: faceP1
        class(IntegerDynamicArrayUDT), allocatable        :: faceP2
        class(IntegerDynamicArrayUDT), allocatable        :: face
        integer(I8)                                       :: nface = 0

        !logical, allocatable                             :: GC(:)

        integer(I8)                                       :: ntot = 0, ngc
        
        class(RealDynamicArrayUDT), allocatable           :: psi, &
            bp, bt, x, y
        class(IntegerDynamicArrayUDT), allocatable        :: cflags, reg, ft
    contains

        ! Initialize
        procedure :: Initialize     => InitializeGACell

    end type

    ! Flux data
    type GAFluxDataUDT

        ! Description
        !============
        ! This data type contains all information regarding flux tubes 
        ! and flux surfaces, and the cells associated with these. 
        ! Fields:
        !
        ! - nFt:            : total number of flux tubes (scalar, to be 
        !                   read in)
        ! - nFs:            : total number of flux surfaces (scalar, to be
        !                   read in)
        ! - fluxtubecellsP  : nFt-by-2 array where the first index 
        !                   is the start index in the fluxtubecells 
        !                   array, and the second the amount of cells of
        !                   the flux tube.
        ! - fluxtubecells   : nCv(number of cells)-by-1 array containing 
        !                   the cell numbers that correspond to flux 
        !                   tubes. 
        ! - fluxtubefacesP  : nFt-by-2 array where the first index 
        !                   is the start index in the fluxtubefaces 
        !                   array, and the second the amount of faces of
        !                   the flux tube.
        ! - fluxtubefsIDs   : flux surface IDs that bound the flux tube
        !                   (nFt-by-2)
        ! - fluxtubefaces   : nFv(number of faces)-by-1 array containing 
        !                   the face numbers that correspond to flux 
        !                   tubes. 
        ! - fluxsurfacefacesP  : nFs-by-2 array where the first index 
        !                   is the start index in the fluxsurfacefaces
        !                   array, and the second the amount of faces of
        !                   the flux surface.
        ! - fluxsurfacefaces   : nFv(number of faces)-by-1 array containing 
        !                   the face numbers that correspond to flux 
        !                   surfaces. 
        ! - fluxsurfaceID   : nv(number of vertices)-by-1 array 
        !                   containing the flux tube ID of each vertex

        ! Logicals and indices
        integer(I8)                         :: nFt = 0
        integer(I8)                         :: nFs = 0

        ! Arrays, flux tube data
        !type(IntegerDynamicArrayUDT), allocatable :: &
        !    fluxtubecellsP(:,:)
        !type(IntegerDynamicArrayUDT), allocatable :: fluxtubecells(:)
        !type(IntegerDynamicArrayUDT), allocatable :: &
        !    fluxtubefacesP(:,:)
        !type(IntegerDynamicArrayUDT), allocatable :: fluxtubefaces(:)
        !type(IntegerDynamicArrayUDT), allocatable :: &
        !    fluxtubefsIDs(:, :), fluxtuberegID(:)
        !logical, allocatable                :: isclosedft(:)

        ! Arrays, flux surface data
        class(IntegerDynamicArrayUDT), allocatable :: fluxsurfacefacesP1
        class(IntegerDynamicArrayUDT), allocatable :: fluxsurfacefacesP2
        class(IntegerDynamicArrayUDT), allocatable :: fluxsurfacefaces
        class(IntegerDynamicArrayUDT), allocatable :: fluxsurfaceID
        !type(IntegerDynamicArrayUDT), allocatable :: &
        !    fluxsurfaceneig(:), fluxsurfaceneigP(:, :)
        class(RealDynamicArrayUDT), allocatable    :: fluxsurfacepsi
        class(IntegerDynamicArrayUDT), allocatable :: fluxsurfacevertsP1
        class(IntegerDynamicArrayUDT), allocatable :: fluxsurfacevertsP2
        class(IntegerDynamicArrayUDT), allocatable :: fluxsurfaceverts

    contains

        ! Initialize
        procedure :: Initialize     => InitializeGAFluxData

    end type  
    
    ! Grid data
    type GAGridDataUDT

        ! Description
        !===========
        ! Data type to collect all other grid related data which is not 
        ! linked to the core data, i.e. the cells, vertices, and faces.
        ! This can be additional info from the grid generator on flux 
        ! tubes, for example, and interconnection data thereof (e.g. 
        ! the cells that belong to a certain flux tube). 
        ! Fields:
        !
        ! - fluxdata            : UDT with all flux data such as flux
        !                       flux tube data, flux surfaces, ... 
        ! - sglegacy            : data from legacy structured grids
        ! - OMPcell, OMPface    : cells and faces belonging to outer mid
        !                       plane
        ! - IMPcell, IMPface    : same, but inner mid plane
        ! - OMPr, OMPz          : points defining line segment of OMP
        ! - IMPr, IMPz          : same but for inner mid plane
        ! - topoflag            : flag indicating the topological mesh 
        !                       type (see also mod_definitions)
        ! - xpointID            : array containing all X-point vertex IDs
        ! - isprimaryxp         : integer array which is 1 if the 
        !                       x-point is a primary x-point (i.e. it 
        !                       lies on a separatrix that is connected
        !                       directly to a core region), otherwise 0
        ! - nxp                 : number of x-points
        ! - spointID            : array containing all strike point vertex IDs
        ! - nsp                 : number of strike points
        ! - opointID            : O point IDs (in the grid)
        ! - nop                 : number of o points
        ! - tpointID            : Tangency point IDs (but here tangency 
        !                       points are only those tangency points 
        !                       that have a closed contour in the grid!)
        ! - ntp                 : number of t points
        ! - spointxpID          : x-point vertex number of the x-point 
        !                       the strike point belongs to
        ! - ndiv                : number of divertor targets/plates
        ! - spointdivID         : list of targets that belong to strike 
        !                       points (one target each)
        ! - tpointdivID         : same as above but for tangency points
        ! - ndivFc              : total number of divertor faces
        ! - divFc               : list of divertor faces
        ! - divFcP              : pointer for the list of divertor faces

        ! Flux data
        type(GAFluxDataUDT)           :: fluxdata

        ! Legacy data of structured grid
        !type(StructuredGridDataUDT) :: sglegacy

        ! Topological mesh type
        !integer(I8)                             :: topoflag
        
        ! OMP & IMP
        !integer(I8), allocatable, dimension(:)  :: OMPcell, OMPface, &
        !    IMPcell, IMPface
        !integer(I8)                             :: nOMPcell, nOMPface, &
        !    nIMPcell, nIMPface
        !real(R8), dimension(1:2)                :: OMPr, OMPz, IMPr, &
        !    IMPz

        ! X-point(s), strike points, o points, ...
        integer(I8), allocatable, dimension(:)  :: xpointID, sepID
        !integer(I8), allocatable, dimension(:)  :: spointID, opointID, tpointID, isprimaryxp, spointxpID, &
        !    divFc, spointdivID, tpointdivID
        !integer(I8), allocatable, dimension(:, :)   :: divFcP
        integer(I8)                             :: nxp, nsep
        !integer(I8)                             :: nsp, nop, ntp, &
        !    ndiv, ndivFc  
            
    contains

        ! Initialize
        procedure :: Initialize     => InitializeGAGridData
        
    end type   
    
    ! Distance function type
    type DistFunctUDT

        ! Description
        !============
        ! Data type to save an interpolant with certain property used to indicate grid regions spatially.

        ! - d_char:  characteristic length
        ! - dist_type: type of distance function
        ! - d_rescale: 
        ! - d_char_type: option to choose characteristic length   

        ! Properties
        character(:), allocatable           :: dist_type
        character(:), allocatable           :: d_char_type
        real(R8)                            :: d_rescale

        ! Interpolant
        type(Coordinates2DDistanceDFUDT)       :: distr
    
    contains

        procedure :: ComputeDistanceFunction


    end type
    
    ! Main grid structure
    !====================
    type GAGridUDT
        ! Description
        !============
        ! Data structure containing all the grid data substructures for grid adaptation. 
        ! Fields:
        !
        ! - vert            : see type definition for description
        ! - face           
        ! - cell
        ! - data
        ! - bnd
        !
        ! Note: the bnd substructure has to be allocated separately

        ! Vertices
        type(GAVertexUDT)                     :: vert

        ! Faces
        type(GAFaceUDT)                       :: face

        ! Cells
        type(GACellUDT)                       :: cell

        ! Additional data
        type(GAGridDataUDT)                   :: data

        ! Boundaries
        !type(BndUDT), allocatable           :: bnd(:)

        ! Distance functions
        type(DistFunctUDT)                    :: fun 
        type(DistFunctUDT)                    :: fun_r
        type(DistFunctUDT)                    :: fun_wall      

    contains

        ! Initialize
        procedure :: Initialize         => InitializeGAGrid


        ! Various operations
        procedure :: CheckVertOrder
        procedure :: ReorderCellConn
        procedure :: GetFsVxFromFsFc
        procedure :: GiveXpoints
        procedure :: GiveSeparatrices
        procedure :: IdentifyAlignedFaces
        procedure :: CheckUnstructuredGrid
        procedure :: RecalcMagn
        procedure :: MergeFS

    end type  

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    contains 

    !------------------------------------------------------------------!
    !                  QUALITY METRIC COMPUTATIONS                     !
    !------------------------------------------------------------------!
    subroutine ComputeQM(qm,grid,options)

        ! Description
        !============
        ! Compute all quality metrics
        ! (Mirror of CalculateCvMetric.m)

        ! Declare variables
        !==================
        ! Arguments
        class(QualityMetric), intent(inout)     :: qm
        type(GridUDT), intent(in)               :: grid
        type(GAoptionsUDT), intent(in)          :: options

        ! Auxiliary
        real(R8) :: vec_n(grid%face%ntot,2), fcH(grid%face%ntot,2), &
         ncpf(grid%face%ntot), fccv(grid%face%ntot,2), &
         fcxx(grid%face%ntot)

        
    end subroutine

    subroutine CalculateQualityMetrics(grid,options,qm)
        ! Description
        !============
        ! Compute metric and criteria of cells necessary to execute grid adaptation.

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(in)           :: grid
        type(GAoptionsUDT), intent(in)      :: options
        type(QualityMetric), intent(inout)  :: qm

        ! Calculate cv metric
        call qm%ComputeQM(grid,options)


        ! Selecting splitting cell

        ! Selecting merging face
        


    end subroutine

    !------------------------------------------------------------------!
    !                        GAGRID OPERATIONS                         !
    !------------------------------------------------------------------!  

    subroutine InitializeGAGrid(grid)

        ! Description
        !============
        ! Initialize the GAgrid substructures
        
        ! Declare variables
        !==================
        class(GAGridUDT) :: grid

        ! Initialize
        !===========
        ! Substructures
        call grid%vert%Initialize()
        call grid%face%Initialize()
        call grid%cell%Initialize()
        call grid%data%Initialize()

    end subroutine

    subroutine InitializeGAVertex(GAvert)

        ! Description
        !============
        ! Initialize the grid adaptation vertex structure (simply empty arrays)        

        ! Declare variables
        !==================
        ! Arguments
        class(GAVertexUDT) :: GAvert

        ! Initialize
        !===========
        if (allocated(GAvert%x)) then 
            ! Assume all allocated
            deallocate(GAvert%x, GAvert%y, GAvert%fieldlineID, &
                GAvert%psi, GAvert%bx, GAvert%by, GAvert%ffbz)
        end if 
        allocate(RealDynamicArrayUDT::GAvert%x, GAvert%y, &
            GAvert%psi, GAvert%bx, GAvert%by, GAvert%ffbz)
        allocate(IntegerDynamicArrayUDT:: GAvert%fieldlineID)



    end subroutine

    subroutine InitializeGAFace(GAface)

        ! Description
        !============
        ! Initialize the grid adapation face structure (simply empty arrays)

        ! Declare variables
        !==================
        ! Arguments
        class(GAFaceUDT)      :: GAface

        ! Initialize
        !===========
        if (allocated(GAface%vert1)) then 
            ! assume all allocated
            deallocate(GAface%vert1, GAface%vert2, GAface%label, GAface%reg, &
                GAface%aligned)
        end if 
        allocate(IntegerDynamicArrayUDT:: GAface%vert1, GAface%vert2, GAface%label, &
             GAface%reg, GAface%aligned)

    end subroutine

    subroutine InitializeGACell(GAcell)

        ! Description
        !============
        ! Initialize the grid adapation cell structure (simply empty arrays)

        ! Declare variables
        !==================
        ! Arguments
        class(GACellUDT)      :: GAcell


        ! Initialize
        !===========
        if (allocated(GAcell%vertP1)) then 
            ! assume all allocated
            deallocate(GAcell%vertP1, GAcell%vertP2, GAcell%vert, &
                GAcell%faceP1, GAcell%faceP2, GAcell%face, GAcell%psi, &
                GAcell%bp, GAcell%bt, &
                GAcell%x, GAcell%y, GAcell%cflags, GAcell%reg, &
                GAcell%ft )
        end if 
        allocate(RealDynamicArrayUDT:: GAcell%psi, GAcell%bp, &
                 GAcell%bt, GAcell%x, GAcell%y )
        allocate(IntegerDynamicArrayUDT:: GAcell%vertP1, GAcell%vertP2, &
                 GAcell%vert, GAcell%faceP1, GAcell%faceP2, GAcell%face, &
                 GAcell%cflags, GAcell%reg, GAcell%ft )

    end subroutine

    subroutine InitializeGAGridData(GAgriddata)

        ! Description
        !============
        ! Initialize the grid adapation data structure (simply empty arrays)

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridDataUDT)      :: GAgriddata

        ! Initialize
        call GAgriddata%fluxdata%Initialize()     

    end subroutine

    subroutine InitializeGAFluxData(GAfluxdata)

        ! Description
        !============
        ! Initialize the grid adapation data structure (simply empty arrays)

        ! Declare variables
        !==================
        ! Arguments
        class(GAFluxdataUDT)      :: GAfluxdata

        ! Initialize
        !===========
        if (allocated(GAfluxdata%fluxsurfacefaces)) then 
            ! assume all allocated
            deallocate(GAfluxdata%fluxsurfacefacesP1, GAfluxdata%fluxsurfacefacesP2, &
                GAfluxdata%fluxsurfacefaces, GAfluxdata%fluxsurfaceID, &
                GAfluxdata%fluxsurfacevertsP1, GAfluxdata%fluxsurfacevertsP2, &
                GAfluxdata%fluxsurfaceverts, GAfluxdata%fluxsurfacepsi)
        end if 
        allocate(RealDynamicArrayUDT:: GAfluxdata%fluxsurfacepsi)
        allocate(IntegerDynamicArrayUDT:: GAfluxdata%fluxsurfacefacesP1, &
                GAfluxdata%fluxsurfacefacesP2, GAfluxdata%fluxsurfacefaces, &
                GAfluxdata%fluxsurfaceID, GAfluxdata%fluxsurfacevertsP1, &
                GAfluxdata%fluxsurfacevertsP2, GAfluxdata%fluxsurfaceverts)


    end subroutine

    subroutine TranslateGridTOGAGrid(grid,GAgrid)
        ! Description
        ! ===========
        ! Translating the information in the GridUDT type to a GAGridUDT type with dynamic arrays

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(in)       :: grid
        type(GAGridUDT), intent(out)    :: GAgrid

        ! Initialize GAGrid
        call GAgrid%Initialize()

        ! Associate
        associate(&
            gc => grid%cell, &
            gf => grid%face, &
            gv => grid%vert, &
            gfd => grid%data%fluxdata, &
            GAc => GAgrid%cell, &
            GAf => GAgrid%face, &
            GAv => GAgrid%vert, &
            GAfd => GAgrid%data%fluxdata &
            )

        ! Give information in grid to GAgrid
        ! Vertex infromation
        GAv%x           = ConstructRealDynamicArray(gv%x)
        GAv%y           = ConstructRealDynamicArray(gv%y)
        GAv%bx          = ConstructRealDynamicArray(gv%bx)
        GAv%by          = ConstructRealDynamicArray(gv%by)
        GAv%psi         = ConstructRealDynamicArray(gv%psi)
        GAv%ffbz        = ConstructRealDynamicArray(gv%ffbz)
        GAv%fieldlineID = ConstructIntegerDynamicArray(gv%fieldlineID)
        GAv%ntot        = gv%ntot

        ! Face information
        GAf%vert1   = ConstructIntegerDynamicArray(gf%vert(:,1))
        GAf%vert2   = ConstructIntegerDynamicArray(gf%vert(:,2))
        GAf%label   = ConstructIntegerDynamicArray(gf%label)
        GAf%reg     = ConstructIntegerDynamicArray(gf%reg)
        GAf%aligned = ConstructIntegerDynamicArray(gf%aligned)
        GAf%ntot    = gf%ntot

        ! Cell information
        GAc%vertP1  = ConstructIntegerDynamicArray(gc%vertP(:,1))
        GAc%vertP2  = ConstructIntegerDynamicArray(gc%vertP(:,2))
        GAc%vert    = ConstructIntegerDynamicArray(gc%vert)
        GAc%faceP1  = ConstructIntegerDynamicArray(gc%faceP(:,1))
        GAc%faceP2  = ConstructIntegerDynamicArray(gc%faceP(:,2))
        GAc%face    = ConstructIntegerDynamicArray(gc%face)
        GAc%cflags  = ConstructIntegerDynamicArray(gc%cflags)
        GAc%reg     = ConstructIntegerDynamicArray(gc%reg)
        GAc%ft      = ConstructIntegerDynamicArray(gc%ft)
        GAc%psi     = ConstructRealDynamicArray(gc%psi)
        GAc%bp      = ConstructRealDynamicArray(gc%bp)
        GAc%bt      = ConstructRealDynamicArray(gc%bt)
        GAc%x       = ConstructRealDynamicArray(gc%x)
        GAc%y       = ConstructRealDynamicArray(gc%y)
        GAc%ntot    = gc%ntot
        GAc%ngc     = gc%ngc
        GAc%nvert   = gc%nvert
        GAc%nface   = gc%nface
        

        ! Grid data - flux surface data
        GAgrid%data%xpointID    = grid%data%xpointID
        GAgrid%data%nxp         = grid%data%nxp
        GAgrid%data%sepID       = grid%data%sepID
        GAgrid%data%nsep        = grid%data%nsep
        GAfd%fluxsurfacefacesP1 = ConstructIntegerDynamicArray(gfd%fluxsurfacefacesP(:,1))
        GAfd%fluxsurfacefacesP2 = ConstructIntegerDynamicArray(gfd%fluxsurfacefacesP(:,2))
        GAfd%fluxsurfacefaces   = ConstructIntegerDynamicArray(gfd%fluxsurfacefaces)
        GAfd%fluxsurfacevertsP1 = ConstructIntegerDynamicArray(gfd%fluxsurfacevertsP(:,1))
        GAfd%fluxsurfacevertsP2 = ConstructIntegerDynamicArray(gfd%fluxsurfacevertsP(:,2))
        GAfd%fluxsurfaceverts   = ConstructIntegerDynamicArray(gfd%fluxsurfaceverts)
        GAfd%fluxsurfaceID      = ConstructIntegerDynamicArray(gfd%fluxsurfaceID)
        GAfd%fluxsurfacepsi     = ConstructRealDynamicArray(gfd%fluxsurfacepsi)
        GAfd%nFs                = gfd%nFs
        GAfd%nFt                = gfd%nFt

    
        end associate

        ! Deallocate grid 
        call DeallocateGrid(grid)


    end subroutine

    subroutine TranslateGAGridTOGrid(grid,GAgrid,options)
        ! Description
        ! ===========
        ! Translating the information in the GridUDT type to a GAGridUDT type with dynamic arrays

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(out)      :: grid
        type(GAGridUDT), intent(inout)  :: GAgrid
        type(GAoptionsUDT), intent(in)  :: options

        ! Auxiliary
        integer(I8) :: i, indFc(GAgrid%face%ntot), fcLbl_loc(GAgrid%face%ntot), nfb
        integer(I8), allocatable :: fcsbnd(:), vxsbnd(:)


        ! Give information in GAgrid to grid
        associate(&
            gc => grid%cell, &
            gf => grid%face, &
            gv => grid%vert, &
            gfd => grid%data%fluxdata, &
            GAc => GAgrid%cell, &
            GAf => GAgrid%face, &
            GAv => GAgrid%vert, &
            GAfd => GAgrid%data%fluxdata &
            )
        ! Initialize grid
        ! Correct sizes for allocation
        gv%ntot     = GAv%ntot    
        gf%ntot     = GAf%ntot
        gc%ntot     = GAc%ntot
        gc%nvert    = GAc%vertP1%Get(GAc%ntot) + Gac%vertP2%Get(GAc%ntot) - 1 
        gc%nface    = GAc%faceP1%Get(GAc%ntot) + Gac%faceP2%Get(GAc%ntot) - 1 
        gfd%nFs     = GAfd%nFs
        gfd%nFt     = GAfd%nFt

        ! Boundary information - TODO
        !fcLbl_loc = GetfcLblGA(GAf, options)
        !nfb = count(fcLbL_loc /= 0)
        !grid%bnd%nface = nfb
        !indFc = (/ (i, i = 1, GAf%ntot) /)
        !allocate(fcsbnd(nfb))
        !fcsbnd = pack(indFc, fcLbL_loc /= 0 )
        !vxsbnd = GetVxsFromFcs(GAf, fcsbnd)
        !grid%bnd%nvert = size(vxsbnd)

        call AllocateGrid(grid)


        ! Vertex information
        gv%x            = GAv%x%GetAllElements()
        gv%y            = GAv%y%GetAllElements()
        gv%bx           = GAv%bx%GetAllElements()
        gv%by           = GAv%by%GetAllElements()
        gv%psi          = GAv%psi%GetAllElements()
        gv%ffbz         = GAv%ffbz%GetAllElements()
        gv%fieldlineID  = GAv%fieldlineID%GetAllElements()

        ! Face information
        gf%vert(:,1)    = GAf%vert1%GetAllElements()
        gf%vert(:,2)    = GAf%vert2%GetAllElements()
        gf%label        = GAf%label%GetAllElements()
        gf%reg          = GAf%reg%GetAllElements()
        gf%aligned      = GAf%aligned%GetAllElements()

        ! Cell information
        gc%vertP(:,1)   = GAc%vertP1%GetAllElements() 
        gc%vertP(:,2)   = GAc%vertP2%GetAllElements() 
        gc%vert         = GAc%vert%GetAllElements() 
        gc%faceP(:,1)   = GAc%faceP1%GetAllElements()
        gc%faceP(:,2)   = GAc%faceP2%GetAllElements()
        gc%face         = GAc%face%GetAllElements()
        gc%cflags       = GAc%cflags%GetAllElements()
        gc%ft           = GAc%ft%GetAllElements()
        gc%psi          = GAc%psi%GetAllElements()
        gc%bp           = GAc%bp%GetAllElements()
        gc%bt           = GAc%bt%GetAllElements()
        gc%x            = GAc%x%GetAllElements()
        gc%y            = GAc%y%GetAllElements()

        ! Grid data - flux surface data
        grid%data%xpointID          = GAgrid%data%xpointID
        grid%data%nxp               = GAgrid%data%nxp
        grid%data%sepID             = GAgrid%data%sepID
        grid%data%nsep              = GAgrid%data%nsep
        gfd%fluxsurfacefacesP(:,1)  = GAfd%fluxsurfacefacesP1%GetAllElements()
        gfd%fluxsurfacefacesP(:,2)  = GAfd%fluxsurfacefacesP2%GetAllElements()
        gfd%fluxsurfacefaces        = GAfd%fluxsurfacefaces%GetAllElements()
        gfd%fluxsurfacevertsP(:,1)  = GAfd%fluxsurfacevertsP1%GetAllElements()
        gfd%fluxsurfacevertsP(:,2)  = GAfd%fluxsurfacevertsP2%GetAllElements()
        gfd%fluxsurfaceverts        = GAfd%fluxsurfaceverts%GetAllElements()
        gfd%fluxsurfaceID           = GAfd%fluxsurfaceID%GetAllElements()
        gfd%fluxsurfacepsi          = GAfd%fluxsurfacepsi%GetAllElements()

        ! Problem need to compute some extra fields for grid 
        ! See what is needed for WriteGOAT: TODO



        end associate
        ! Deallocate GAgrid - maybe not problematic




    end subroutine 

    subroutine CheckVertOrder(grid,is_ordered,cells)

        ! Description
        !============
        ! Check whether the vertices and faces are ordered correctly,
        ! following the right turning scheme  
        ! is_ordered: indicates whether the vertices and faces of a cell are ordered
        ! cells: indicate which cells need to be checked

        ! Declare variables
        !==================
        class(GAGridUDT)            :: grid
        logical, intent(inout)      :: is_ordered(grid%cell%ntot) &
        , cells(grid%cell%ntot)

        ! Auxiliary
        real(R8) ::  sin1, sin2
        integer(I8) :: ic, v1, v2, f1, v1f, v2f, nv, i
        integer(I8), allocatable, dimension(:) :: tv, tf, v1n, v2n
        real(R8), allocatable, dimension(:) :: vec_vxx, vec_vxy, &
            vec_fcx, vec_fcy, fcX, fcY, vx, vy, cx, cy

        ! Associate
        associate(&
            c  => grid%cell, &
            f  => grid%face, &
            v  => grid%vert &
            )

        ! Initialize
        allocate(fcX(f%ntot), fcY(f%ntot), v1n(f%ntot), &
            v2n(f%ntot), cx(c%ntot), cy(c%ntot), vx(v%ntot), &
            vy(v%ntot)) 
        v1n = f%vert1%GetAllElements()
        v2n = f%vert2%GetAllElements()
        vx = v%x%GetAllElements()
        vy = v%y%GetAllElements()
        cx = c%x%GetAllElements()
        cy = c%y%GetAllElements()
        fcX = 0.5_R8 * (vx(v1n) + vx(v1n))
        fcY = 0.5_R8 * (vy(v2n) + vy(v2n))

        ! Lopp over the cells
        is_ordered = .true.
        do ic = 1, c%ntot
            if (cells(ic)) then
                ! Get vertices and faces
                tv = GetCellVertGA(c, ic)
                tf = GetCellFaceGA(c, ic)

                ! First check whether vertices and faces are connected according to correct ordering
                nv = size(tv)
                do i = 1, nv-1
                    v1 = tv(i)
                    v2 = tv(i+1)
                    f1 = tf(i)
                    v1f = v1n(f1)
                    v2f = v2n(f1)

                    if (((v1 /= v1f).and.(v1 /= v2f)).or. &
                          ((v2 /= v1f).and.(v2 /= v2f)) ) then
                            is_ordered(ic) = .false.
                            exit
                    end if 

                end do
                i = nv
                v1 = tv(i)
                f1 = tf(i)
                v1f = v1n(f1)
                v2f = v2n(f1)
                if (((v1 /= v1f).and.(v1 /= v2f)).or. &
                    ((tv(1) /= v1f).and.(tv(1) /= v2f)) ) then
                        is_ordered(ic) = .false.
                end if 


                ! If still ok, check the direction
                if (is_ordered(ic)) then
                    ! Vectors from cell center to vertices
                    vec_vxx = vx(tv) - cx(ic)
                    vec_vxy = vy(tv) - cy(ic)

                    ! Vectors from cell center to face centers
                    vec_fcx = fcX(tf) - cx(ic)
                    vec_fcy = fcY(tf) - cy(ic) 

                    do i = 1, nv-1
                        !check sinus between vector to vertex (a) and vector to face
                        !center (b)
                        !calculate angles (sin = |a x b| / norm(a)*norm(b))
                        !with |a x b | = ax*by - bx*ay
                        sin1 = vec_vxx(i)*vec_fcy(i) - vec_fcx(i) * vec_vxy(i)

                        ! Should be positive
                        if (sin1.lt.0.0_R8) then
                            is_ordered(ic) = .false.
                            exit
                        else 
                            !check sinus between vector to face center (b) and vector
                            !to second vertex (a)
                            sin2 = vec_fcx(i)*vec_vxy(i+1) - vec_vxx(i+1)*vec_fcy(i)
                            
                            if (sin2.lt.0.0_R8) then
                                is_ordered(ic) = .false.
                                exit
                            end if

                        end if

                    end do 
                end if

            end if
        end do

        end associate

    end subroutine 

    subroutine ReorderCellConn(grid, is_ordered)

        ! Description
        !============
        ! Re-orders vertices and faces of cells in an order such that they form a
        ! chain as you would walk over the cell boundary 

        ! Declare variables
        !==================
        class(GAGridUDT), intent(inout)    :: grid
        logical, intent(in)                :: is_ordered(grid%cell%ntot)

        ! Declare variables
        !==================
        integer(I8) :: ic, nv, fcs(1:2), vs(1:2), s, n, i, j, ntv(1:20), &
            ntf(1:20), indf(1:2)
        integer(I8), allocatable, dimension(:) :: tv, tf,  &
            v1n, v2n, range
        integer(I8), allocatable, dimension(:,:) :: vf, indv
        real(R8) :: vec_start(1:2), vec_face1(1:2), vec_face2(1:2), &
            sin1, sin2
        real(R8), allocatable :: fcX(:), fcY(:), vx(:), vy(:), cx(:), cy(:)
        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert  &
            )
        
        ! Compute face centers
        allocate(fcX(f%ntot), fcY(f%ntot), v1n(f%ntot), &
            v2n(f%ntot), vx(v%ntot), vy(v%ntot))
        v1n = f%vert1%GetAllElements()
        v2n = f%vert2%GetAllElements()        
        vx = v%x%GetAllElements()
        vy = v%y%GetAllElements()
        cx = c%x%GetAllElements()
        cy = c%y%GetAllElements()        
        fcX = 0.5_R8 * (vx(v1n) + vx(v2n))
        fcY = 0.5_R8 * (vy(v1n) + vy(v2n))

        ! Loop over cells which are not ordered
        ! Could be implement to pack the cells and loop over them
        do ic = 1, c%ntot
            if (.not.is_ordered(ic)) then
                ! Get vertices and faces
                tv = GetCellVertGA(c, ic)
                tf = GetCellFaceGA(c, ic)
                nv = size(tv)

                ! New vertices and faces
                ntv(1:nv) = 0
                ntf(1:nv) = 0

                ! Starting point
                ntv(1) = tv(1)
                
                allocate(vf(nv,2))
                vf(:,1) = v1n(tf)
                vf(:,2) = v2n(tf)

                ! Find the two faces connected to that vertex and the other vertex connected to that face
            
                ! Find the first face
                allocate(indv(nv,2))
                indv(:,1) = (/ (i, i = 1, nv)/)
                indv(:,2) = indv(:,1)
                if (count(vf == ntv(1)).eq.1) then
                    call gdErrorHandler('ReorderCellConn: vertex with one occurence')
                end if 
                indf(1:2) = pack(indv, vf == ntv(1))
                fcs(1:2) = tf(indf)
                !fcs2 = pack(indv, vf == ntv(1))
                !fcs(1) = fcs1(1)
                !fcs(2) = fcs2(1)


                ! Start right hand turning
                ! Define vector from cell centers to points
                vec_start(1) = vx(ntv(1)) - cx(ic)
                vec_start(2) = vy(ntv(1)) - cy(ic)

                ! Face1
                vec_face1(1) = fcX(fcs(1)) - cx(ic)
                vec_face1(2) = fcY(fcs(1)) - cy(ic)

                ! Face2
                vec_face2(1) = fcX(fcs(2)) - cx(ic)
                vec_face2(2) = fcY(fcs(2)) - cy(ic)
    
                ! Calculate angles (sin = |a x b| / norm(a)*norm(b))
                ! with |a x b | = ax*by - bx*ay
                ! dividing by norm is needed because norm is always positive and
                ! only the sign matters
                ! a = always vector from cell center to starting point
                sin1 = vec_start(1)*vec_face1(2) - vec_face1(1) * vec_start(2)
                sin2 = vec_start(1)*vec_face2(2) - vec_face2(1) * vec_start(2)
                
                ! Choose the one with positive sign, that right turning
                if (sin1 .gt. 0.0_R8) then
                    ntf(1) = fcs(1) ! First face
                elseif (sin2 .gt. 0.0_R8) then
                    ntf(1) = fcs(2)
                end if

                ! After the direction is fixed, implementation should be in a loop
                do i = 2, nv
                    ! Find the next vertex
                    vs(1) = v1n(ntf(i-1)) 
                    vs(2) = v2n(ntf(i-1))

                    if (.not. any(ntv(1:nv) == vs(1))) then
                        ntv(i) = vs(1)
                    else
                        ntv(i) = vs(2)
                    end if

                    ! Find the next face
                    indf(1:2) = pack(indv, vf == ntv(i))
                    fcs(1:2) = tf(indf)

                    if (.not. any(ntf(1:nv) == fcs(1))) then
                        ntf(i) = fcs(1)
                    else
                        ntf(i) = fcs(2)
                    endif

                end do

                ! Plug in the new verts and faces in cell%vert and cell%face
                s = c%vertP1%GetSingleElement(ic)
                n = c%vertP2%GetSingleElement(ic)
                allocate(range(n))
                range = (/ (j, j = s,s+n-1)/)
                call c%vert%SetMultipleElements(range,ntv(1:nv))
                call c%face%SetMultipleElements(range,ntf(1:nv))
                deallocate(range) 
                deallocate(vf) 
                deallocate(indv)

            end if

        end do

        ! Deallocate 
        deallocate(fcX, fcY)

        end associate

    end subroutine

    subroutine GetFsVxFromFsFc(grid)
        ! Description
        !============
        ! Get fsVx from fsFc

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout) :: grid

        ! Auxiliary
        integer(I8) :: nv_counter, ifs, nv, nf
        integer(I8), allocatable, dimension(:) :: fcs, vxs, verts, fsVx, v1n, v2n
        integer(I8), allocatable :: fsVxP(:,:)
         



        ! Associate
        associate(&
            fd => grid%data%fluxdata, &
            f  => grid%face, &
            v  => grid%vert &
            )
        
        ! Initialize
        allocate(fsVx(v%ntot), fsVxP(fd%nFs,2),verts(v%ntot), &
            v1n(f%ntot), v2n(f%ntot))
        v1n = f%vert1%GetAllElements()
        v2n = f%vert2%GetAllElements()
        fsVx = 0
        fsVxP = 0
        nv_counter = 0
        verts = 0
        

        do ifs = 1, fd%nFs
            ! Get vertices from flux surface
            fcs = GetFluxSurfaceFcsGA(fd, ifs)
            nf = size(fcs)
            verts(1:nf) = v1n(fcs) 
            verts(nf+1:nf*2) = v2n(fcs)
            call Unique(verts(1:nf*2), vxs)
            nv = size(vxs)

            ! Fill in into fsVx and fsVxP
            fsVxP(ifs,:) = [nv_counter+1, nv];
            fsVx(fsVxP(ifs,1):fsVxP(ifs,1)+fsVxP(ifs,2)-1) =  vxs;
            nv_counter = nv_counter + nv;        
        end do

        ! Trim - first implement GAGrid
        call fd%fluxsurfaceverts%SetAllElementsArray(fsVx(1:nv_counter))
        call fd%fluxsurfacevertsP1%SetAllElementsArray(fsVxP(1:fd%nFs,1))
        call fd%fluxsurfacevertsP2%SetAllElementsArray(fsVxP(1:fd%nFs,2))

        end associate


    end subroutine

    subroutine GiveXpoints(grid,use_sep,cvLookUp)

        ! Description
        !============
        ! Gives the Xpoint(s). Depending the grid information different methods are used. If all vertices have a fieldlineID the routine DetermineXPointsGA is used, otherwise the xpoint is determined by checking whether a vertex is surrounded by more then three regions and more than four cells. 

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)  :: grid
        logical, intent(in)             :: use_sep  
        integer(I8), allocatable, optional :: cvLookUp(:)            

        ! Auxiliary
        integer(I8)                         :: i, j, iv, xpoints(1:100), counter, & 
            n , nc
        integer(I8), allocatable            :: vxs(:), cells(:), &
            regions(:), creg(:)
        logical :: use_sepID, start, use_fieldlineID, use_nsep
        

        ! Initialize
        counter = 0
        xpoints = 0
        

        ! Associate
        associate(&
            c   => grid%cell, &
            v   => grid%vert, &
            fd  => grid%data%fluxdata &
            )

            ! Initialize
            creg = c%reg%GetAllElements()

            ! Find separatrix if not given
            if (.not.use_sep) then
                use_nsep = .false.
                use_sepID = .false.
                start = .true.
                call GiveSeparatrices(grid,use_nsep,use_sepID,start,cvLookUp)
            end if

            ! Determine whether to use the vert.fieldlineID to determine the Xpoint
            use_fieldlineID = .false.
            if (allocated(v%fieldlineID)) then
                if (v%fieldlineID%Size().eq.v%ntot) then
                    use_fieldlineID = .true.
                end if
            end if

            !if (use_fieldlineID) then

                !!! NOT USED BECAUSE ALGO USES vert%neig WHICH IS NOT PRESENT IN GAGRID !!
            !    call DetermineXPointsGA(xpind, nxpind, order, grid)
            !    grid%data%xpointID  = xpind
            !    grid%data%nxp       = nxpind 

            if ((allocated(v%fieldlineID))) then!.and.(allocated(v%cellP1))) then

                ! Only check the vertices on the separatrices
                cvLookUp = GetCvLookUpGA(c)
                do i = 1, grid%data%nsep
                    n = fd%fluxsurfacevertsP2%Get(grid%data%sepID(i))
                    vxs = GetFluxSurfaceVxsGA(fd, grid%data%sepID(i))

                    do j = 1, n
                        iv = vxs(j)
                        cells = GetVertCellGA(c, iv, cvLookUp)
                        call Unique(creg(cells), regions)
                        nc = size(cells)

                        if ((size(regions).ge.3) .and. (nc.gt.4)) then
                            counter = counter + 1
                            xpoints(counter) = iv
                        end if

                    end do 

                end do

                ! Saving 
                grid%data%xpointID  = xpoints(1:counter)
                grid%data%nxp       = counter

            else
                
                ! Check all vertices
                cvLookUp = GetCvLookUpGA(c)
                do iv = 1, v%ntot
                    cells = GetVertCellGA(c, iv, cvLookUp)
                    call Unique(creg(cells), regions)
                    nc = size(cells)

                    if ((size(regions).ge.3) .and. (nc.gt.4)) then
                        counter = counter + 1
                        xpoints(counter) = iv
                    end if

                end do 

                ! Saving 
                grid%data%xpointID  = xpoints(1:counter)
                grid%data%nxp       = counter

            end if

            if (grid%data%nxp == 0) then
                print *, 'GiveXpoints: Warning: no Xpoint found!'
            endif 

        end associate

    end subroutine 

    !subroutine DetermineXPointsGA(xpind, nxpind, order, grid)

        ! Description
        !============
        ! This routine determines the indices (and their number) based 
        ! on the information given in the GAgrid structure. Note that 
        ! the location of the x-point should already be available 
        ! through output of the grid generator, but we provide a routine 
        ! here to recompute the location if it wasn't saved explicitly
        ! in the grid structure. Additionally, we determine the order of
        ! the x-point (see algorithm description below)

        ! Algorithm
        !==========
        ! Basically, we simply check if a vertex has 4 or more 
        ! neighbours with the same flux surface ID (which is non-zero)
        ! as the current vertex. The order is then simply determined by 
        ! the number of separatrix legs, i.e. 
        !
        !       o = n/2 - 1,
        !
        ! where n is the number of vertices with the same ID and o is 
        ! the order. 

        ! Initialize
        !===========
        ! The usual
        !implicit none 

        ! Declare variables
        !==================
        ! Arguments
        !integer(I8)                         :: nxpind
        !integer(I8), allocatable            :: xpind(:), order(:)
        !class(GAGridUDT), intent(in)           :: grid

        ! Loop variables
        !integer(I8)                         :: i

        ! Auxiliary variables 
        !integer(I8)                         :: tfID, ncIDs
        !integer(I8), allocatable            :: temporder(:), &
        !    tempxpind(:), tvn(:), tvnfID(:)

        ! Data

        ! Initialize
        !===========
        ! Associate
        !associate(&
        !    vert        => grid%vert)
        !fID = vert%fieldlineID%GetAllElements()

        ! Check allocation - shouldn't be the case as size unknown a 
        ! priori
        !if (allocated(xpind)) then
        !    ! Deallocate
        !    deallocate(xpind)
        !end if
        !if (allocated(order)) then 
        !    ! Deallocate
        !    deallocate(order)
        !end if

        ! Allocate temporary arrays (too big, trim later)
        !nxpind = 0
        !allocate(tempxpind(vert%ntot))
        !allocate(temporder(vert%ntot))

        ! Determine x-points
        !===================
        ! Loop over all vertices
        !do i = 1, vert%ntot 
            ! Get the current field line ID, skip if zero
        !    tfID = fID(i)
        !    if (tfID == 0) then
        !        cycle ! skip rest of the loop for this index
        !    end if

            ! Get the vertex neighbours
        !    allocate(tvn(vert%neigP(i, 2)))
        !    tvn = vert%neig(vert%neigP(i, 1):&
        !        (vert%neigP(i, 1) + vert%neigP(i, 2)-1))

        !    ! Get their IDs
        !    allocate(tvnfID(vert%neigP(i, 2)))
        !    tvnfID = fID(tvn)

        !    ! Get the number of common IDs
        !    ncIDs = count(tvnfID == tfID)

            ! If equal or larger than 4, add as x-point
        !    if (ncIDs >= 4) then 
                ! Check if it is a multiple of 2
        !        if (modulo(ncIDs, 2) .ne. 0) then
                    ! Uneven number, throw error
        !            call gdErrorHandler('DetermineXPoints: '&
        !            // 'supposed x-point has uneven number' &
        !            // 'of vertices with same ID, check grid' &
        !            // ' consistency')
        !        end if

                ! Update the counter
        !        nxpind = nxpind + 1

                ! Add the x-point
        !        tempxpind(nxpind) = i 

                ! Compute the order
        !        temporder(nxpind) = ncIDs/2 - 1

        !    end if
            
            ! Deallocate
        !    deallocate(tvn, tvnfID)

        !end do

        ! Output
        !=======
        ! Allocate
        !allocate(xpind(nxpind), order(nxpind))

        ! Set output
        !xpind = tempxpind(1:nxpind)
        !order = temporder(1:nxpind)

        ! Housekeeping
        !=============
        ! Deallocate 
        !deallocate(tempxpind, temporder)

        ! Deassociate
        !end associate

    !end subroutine

    subroutine GiveSeparatrices(grid,use_nsep,use_sepID,start,cvLookUp)

        ! Description
        !============
        ! Determines the ID of the separatrices depending on the given information.
        ! The argument start is used when the separatrices need to be determined at the start of the program, so with little rest information.

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT)        :: grid
        logical, intent(in)                  :: use_nsep, use_sepID, start
        integer(I8), allocatable, optional :: cvLookUp(:)

        ! Initialize
        logical    :: use_xpointID     
        integer(I8), allocatable :: vxs(:), fcs(:), cvs(:), creg(:)
        integer(I8) :: ifs, ix, nv, sepIDloc(1:100), &
            counter, counter_dummy, nf, step, i, ifc

        ! Initialize
        sepIDloc = 0
        counter = 0
        counter_dummy = 0
        
        

        ! Associate
        associate(&
            c   => grid%cell, &
            fd  => grid%data%fluxdata &
            )

            ! Initialize
            creg = c%reg%GetAllElements()
            use_xpointID = .false.
            

            if (.not.present(cvLookUp)) then 
                cvLookUp = GetCvLookUpGA(c)
            end if
            
            ! Determine separatrices IDs
            if (start) then 
                if (.not.use_nsep) then

                    ! Checking if an Xpoint was already determined
                    if (allocated(grid%data%xpointID)) then
                        if (size(grid%data%xpointID) /= 0) then
                            use_xpointID = .true.
                        end if
                    end if

                    if (use_xpointID) then
                        ! Based on Xpoints
                        do ifs = 1, fd%nFs
                            nv = fd%fluxsurfacevertsP2%Get(ifs)
                            vxs = GetFluxSurfaceVxsGA(fd, ifs)

                            do ix = 1, grid%data%nxp
                                if (any(grid%data%xpointID(ix).eq.vxs)) then
                                    counter = counter + 1
                                    sepIDloc(counter) = ifs
                                end if
                            end do
                        end do

                    else

                        do ifs = 1, fd%nFs
                            fcs = GetFluxSurfaceFcsGA(fd, ifs)
                            nf = fd%fluxsurfacefacesP2%Get(ifs)

                            step = 1
                            sepIDloc = DetectSepID(ifs, c, fcs, sepIDloc, nf, step, cvLookUp, counter)
                            counter_dummy = counter + 1
                            if (sepIDloc(counter_dummy) /= 0) then
                                counter = counter + 1
                            end if 

                        end do 

                    end if

                    ! Save
                    grid%data%nsep = counter
                    grid%data%sepID = sepIDloc(1:counter)

                
                elseif (use_nsep .and. grid%data%nsep.eq.1) then
                    ! Determin separatrix for SN-case
                    ! Check whether previous sepID is still separatrix
                    if (use_sepID) then
                        ifs = grid%data%sepID(1)
                        nf = fd%fluxsurfacefacesP2%Get(ifs)
                        fcs = GetFluxSurfaceFcsGA(fd, ifs)
                        step = 1
                        do i = 1, nf, step
                            ifc = fcs(i)
                            cvs = GetFaceCellGA(c,ifc,cvLookUp)
                            if (size(cvs).eq.2) then
                                if (creg(cvs(1)) /= creg(cvs(2))) then
                                    ! Noting has changed and should not be possible to created new separatrices with grid adaptations
                                    return
                                end if
                            else
                                exit
                            end if
                        end do
                    end if 

                    ! Determine new separatrices (line 99)
                    sepIDloc = 0
                    do ifs = 1, fd%nFs
                        nf = fd%fluxsurfacefacesP2%Get(ifs)
                        fcs = GetFluxSurfaceFcsGA(fd, ifs)
                        step = 1

                        do i = 1, nf, step
                            ifc = fcs(i)
                            cvs = GetFaceCellGA(c,ifc,cvLookUp)
                            if (size(cvs).eq.2) then
                                if (creg(cvs(1)) /= creg(cvs(2))) then
                                    sepIDloc(1) = ifs
                                    exit
                                end if 
                            else
                                exit
                            end if
                        end do
                        if (sepIDloc(1) == ifs) then
                            grid%data%sepID = sepIDloc(1)
                            exit
                        end if

                    end do 

                    if (sepIDloc(1) .eq. 0) then
                        call gdErrorHandler('GiveSeparatrices: No separatrix found')
                    end if

                elseif (use_nsep.and.grid%data%nsep.gt.1) then
                    ! Determine separatrix for general case
                    sepIDloc = 0
                    counter = 0
                    ! Check whether previous is still correct
                    if (.not.any(grid%data%sepID /= 0)) then
                        do i = 1, grid%data%nsep
                            ifs = grid%data%sepID(i)
                            nf = fd%fluxsurfacefacesP2%Get(ifs)
                            fcs = GetFluxSurfaceFcsGA(fd, ifs)
                            step = 1
                            sepIDloc = DetectSepID(ifs, c, fcs, sepIDloc, nf, step, cvLookUp, counter)
                            counter_dummy = counter + 1
                            if (sepIDloc(counter_dummy) /= 0) then
                                counter = counter + 1
                            end if 
                        end do
                        if (counter.eq.grid%data%nsep) then
                            grid%data%sepID = sepIDloc(1:counter)
                            return
                        end if
                    end if

                    do ifs = 1, fd%nFs
                        if (.not.any(ifs == grid%data%sepID)) then
                            nf = fd%fluxsurfacefacesP2%Get(ifs)
                            fcs = GetFluxSurfaceFcsGA(fd, ifs)
                            step = 1
                            sepIDloc = DetectSepID(ifs, c, fcs, sepIDloc, nf, step, cvLookUp, counter)
                            counter_dummy = counter + 1
                            if (sepIDloc(counter_dummy) /= 0) then
                                counter = counter + 1
                            end if
                            if (sepIDloc(grid%data%nsep) /= 0 ) then ! all separatrices are found
                                exit
                            else 
                                call gdErrorHandler("GiveSeparatrices: not all separatrices were found.")
                                ! If above would not be sufficient, check give_iFs_sep.m line 194
                            end if
                        end if

                    end do

                    ! Save
                    grid%data%sepID = sepIDloc(1:counter)
                    !(nsep already there)

                end if 

            else ! no start, hopefully more efficient
                sepIDloc = 0
                counter = 0
                if  (.not.use_nsep) then
                    if (allocated(grid%data%xpointID).and. size(grid%data%xpointID) /= 0) then
                        do ifs = 1, fd%nFs
                            nv = fd%fluxsurfacevertsP2%Get(ifs)
                            vxs = GetFluxSurfaceVxsGA(fd, ifs)
                            do ix = 1, grid%data%nxp
                                if (any(grid%data%xpointID(ix)==vxs)) then
                                    counter = counter + 1
                                    sepIDloc(counter) = ifs
                                end if 
                            end do 
                        end do
                    else ! line 252
                        do ifs = 1, fd%nFs
                            nf = fd%fluxsurfacefacesP2%Get(ifs)
                            fcs = GetFluxSurfaceFcsGA(fd, ifs)
                            step = max(2,nint(nf/5.0_R8))
                            sepIDloc = DetectSepID(ifs,c, fcs, sepIDloc, nf, step, cvLookUp, counter)
                            counter_dummy = counter + 1
                            if (sepIDloc(counter_dummy) /= 0) then
                                    counter = counter + 1
                            end if

                        end do

                    end if

                    ! Save
                    grid%data%nsep = counter
                    grid%data%sepID = sepIDloc(1:counter)
                elseif (use_nsep.and.(grid%data%nsep.eq.1)) then
                    ! Check whether the previous is separatrix
                    if (grid%data%sepID(1) /= 0) then
                            nf = fd%fluxsurfacefacesP2%Get(ifs)
                            fcs = GetFluxSurfaceFcsGA(fd, ifs)
                            step = max(2,nint(nf/5.0_R8))   
                            do i = 3, nf-2, step
                                ifc = fcs(i)
                                cvs = GetFaceCellGA(c,ifc,cvLookUp)
                                if (size(cvs).eq.2) then
                                    if (creg(cvs(1))/=creg(cvs(2))) then
                                        return ! Separatrix is the same
                                    end if
                                else
                                    exit ! exit loop because this surface is boundary
                                end if

                            end do                     
                    end if

                    ! Determine new separatrix if necessary
                    sepIDloc = 0
                    counter = 0
                    do ifs = 1, fd%nFs
                        nf = fd%fluxsurfacefacesP2%Get(ifs)
                        fcs = GetFluxSurfaceFcsGA(fd, ifs)
                        step = 1   
                        do i = 3, nf-2, step 
                            ifc = fcs(i)
                            cvs = GetFaceCellGA(c, ifc, cvLookUp)
                            if (size(cvs).eq.2) then 
                                if (creg(cvs(1))/=creg(cvs(2))) then
                                    counter = counter + 1
                                    sepIDloc(counter) = ifs
                                    exit
                                end if
                            else
                                exit                               
                            end if  
                        end do
                        if (sepIDloc(counter).eq.ifs) then
                            exit
                        end if
                    end do 

                    ! Save
                    grid%data%sepID = sepIDloc(1)

                elseif (use_nsep.and.(grid%data%nsep.gt.1)) then
                    sepIDloc = 0
                    counter = 0
                    if (.not.any(grid%data%sepID.eq.0)) then
                        do i = 1, grid%data%nsep
                            ifs = grid%data%sepID(i)
                            nf = fd%fluxsurfacefacesP2%Get(ifs)
                            fcs = GetFluxSurfaceFcsGA(fd, ifs) 
                            step = 1
                            sepIDloc = DetectSepID(ifs,c, fcs, sepIDloc, nf, step, cvLookUp, counter)
                            counter_dummy = counter + 1
                            if (sepIDloc(counter_dummy) /= 0) then
                                    counter = counter + 1
                            end if                           
                        end do
                        if (counter .eq. grid%data%nsep) then
                            return ! All stays the same
                        end if
                    end if

                    ! Determine new separatrices if necessary
                    do ifs = 1, fd%nFs
                        if (.not.any(sepIDloc(1:grid%data%nsep).eq.ifs)) then
                            nf = fd%fluxsurfacefacesP2%Get(ifs)
                            fcs = GetFluxSurfaceFcsGA(fd, ifs) 
                            step = 1   
                            sepIDloc = DetectSepID(ifs,c, fcs, sepIDloc, nf, step, cvLookUp, counter)
                            counter_dummy = counter + 1
                            if (sepIDloc(counter_dummy) /= 0) then
                                    counter = counter + 1
                            end if  
                            if (sepIDloc(grid%data%nsep)/= 0) then
                                exit ! all separatrices are found
                            end if
                        end if
                    end do

                    ! If not all where found - extra check possible see line 390 in give_iFs_sep.m
                    if (sepIDloc(grid%data%nsep).eq.0) then
                       call gdErrorHandler("GiveSeparatrices: not all separatrices where found")
                    end if

                    ! Save
                    grid%data%sepID = sepIDloc(1:counter)

                end if

            end if 

        end associate

    end subroutine

    subroutine IdentifyAlignedFaces(grid,options,magneticField)
       
        ! Description
        !============
        ! Identify adhoc whether faces are aligned with the magnetic field or
        ! whether they were supposed to be aligned by the grid generator

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)            :: grid
        type(GAoptionsUDT), intent(in)          :: options
        type(MagneticFieldUDT), intent(in)   :: magneticField
        
        ! Auxiliary
        integer(I8) :: ifc, lim, i, ic, nf, indmax, indmin, &
            rface3, pface3(1:2), ind3(1:3), ind4(1:4), f1, f2, n_al, &
            fcs_al3(1:3), fcs_al2(1:2) , indmax2
        integer(I8), allocatable :: facealigned(:), fcLbl_loc(:), indFc(:), &
         tf(:), fcs(:), fb(:), v1n(:), v2n(:), ccflags(:)
        logical, allocatable :: b_flag(:)
        real(R8) :: abs_cos_loc3(1:3), abs_cos_loc4(1:4), dpsi_f(1:4), &
            dpsi_r, cos2, dpsi3(1:3), dpsi2(1:2)
        real(R8), allocatable, dimension(:) :: fcX, fcY, &
            dpsidx, dpsidy, Bx, By, Btot, abs_cos, &
            t1x, t2x, t1y, t2y, Bnorm, cosB, vpsi

        ! Associate
        associate( &
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        ! Initialize
        allocate(facealigned(f%ntot), fcX(f%ntot), fcY(f%ntot), &
        dpsidx(f%ntot), dpsidy(f%ntot), Bx(f%ntot), By(f%ntot), &
        Btot(f%ntot), Bnorm(f%ntot), cosB(f%ntot), t1x(f%ntot),&
        t2x(f%ntot), t1y(f%ntot), t2y(f%ntot), &
        fcLbl_loc(f%ntot), b_flag(f%ntot), indFc(f%ntot), abs_cos(f%ntot), &
        v1n(f%ntot), v2n(f%ntot), vpsi(v%ntot), ccflags(c%ntot))

        v1n = f%vert1%GetAllElements()
        v2n = f%vert2%GetAllElements()
        fcX = 0.5_R8*(v%x%GetMultipleElements(v1n) + v%x%GetMultipleElements(v2n))
        fcY = 0.5_R8*(v%y%GetMultipleElements(v1n) + v%y%GetMultipleElements(v2n))
        facealigned = 0

        ! Sort faces - Already done in GAInit

        ! Get Magneticfield
        call magneticField%interp%Evaluate(fcX,fcy,1,0,dpsidx)
        call magneticField%interp%Evaluate(fcX,fcy,0,1,dpsidy)  
        Bx = -dpsidy
        By = dpsidx
        Btot = sqrt(Bx**2 + By**2)

        abs_cos = 0
        ! Tangential vector
        t1x = v%x%GetMultipleElements(v1n)
        t2x = v%x%GetMultipleElements(v2n)
        t1y = v%y%GetMultipleElements(v1n)
        t2y = v%y%GetMultipleElements(v2n)

        t1x = t2x - t1x
        t1y = t2y - t1y

        Bnorm = sqrt(t1x**2 + t1y**2)*Btot
        cosB = t1x*Bx + t1y*By
        
        ! Loop because of intrinsic abs an
        do ifc = 1, f%ntot
            abs_cos(ifc) = abs(max(-1.0_R8,min(1.0_R8,cosB(ifc))) / Bnorm(ifc))
        end do

        ! Locally generalize face labels
        fcLbl_loc = GetfcLblGA(f,options)

        ! Boundary faces
        b_flag = fcLbl_loc /= 0

        ! Criteria for trapezoids
        if (options%vesselmode) then
            lim = 10
        else
            lim = 100
        end if

        ! Loop over internal cells - these are 
        !   Quads with two aligned faces
        !   Triangles with one aligned face
        ! Initialize
        ind3 = (/ (i, i=1,3) /)
        ind4 = (/ (i, i=1,4) /)
        vpsi = v%psi%GetAllElements()

        do ic = 1, c%ntot
            tf = GetCellFaceGA(c, ic)
            nf = size(tf)

            ! Get the type of cell
            if (nf.eq.3) then ! Triangle
                ! Minimal angles with magneticfield
                abs_cos_loc3 = abs_cos(tf)

                ! Maximal absolute cosine
                indmax = maxloc(abs_cos_loc3,1)
     
                rface3 = tf(indmax)
                pface3 = tf(pack(ind3, tf/=rface3))

                facealigned(rface3) = 1
                facealigned(pface3) = 0

            elseif (nf.eq.4) then
                ! Minimal angles with magneticfield
                abs_cos_loc4 = abs_cos(tf)

                ! Maximal absolute cosine
                indmax = maxloc(abs_cos_loc4,1)
                f1 = tf(indmax)

                ! Get dpsi_f for opposite faces in increasing psi value
                dpsi_f = abs(vpsi(v1n(tf)) - vpsi(v2n(tf)));
                indmin = minloc(dpsi_f,1);

                ! Put f1 on aligned
                facealigned(f1) = 1;

                ! Get opposite face
                abs_cos_loc4(indmax) = 0.0_R8;
                indmax2 = maxloc(abs_cos_loc4,1);
                f2 = tf(indmax2);

                ! Criteria for f2 being aligned with psi => possibility of trapezoidal cell
                dpsi_r = dpsi_f(indmax2)/dpsi_f(indmin);

                cos2 = abs_cos_loc4(indmax2);

                ! Decide alignment based on case and criteria
                if (options%vesselmode) then

                    if (b_flag(f2)) then
                        if ((abs(cos2).gt.0.999_R8) .and. (dpsi_r.lt.lim)) then
                            facealigned(f2) = 1
                        end if 
                    else ! so no boundary face, so should be aligned with fieldlines by GG
                        facealigned(f2) = 1
                    end if

                else ! not vesselmode

                    if (b_flag(f2)) then
                        if (fcLbl_loc(f2).lt.4) then
                            if (abs(cos2).gt. 0.9 .and. dpsi_r.lt.lim*10) then
                                facealigned(f2) = 1
                            end if
                        else  ! target face (fclabel >= 4)
                            if (abs(cos2).gt.0.999_R8 .and. dpsi_r.lt.lim) then
                                facealigned(f2) = 1
                            end if
                        end if
                    end if

                end if

            else 
                call gdErrorHandler("IdentifyAlignedFace: cell is quad nor triangle, not supported")
            end if

        end do


        ! Put all core boundary faces as aligned
        indFc = (/ (i, i = 1, f%ntot)/)
        facealigned(pack(indFc,fcLbl_loc==2)) = 1

        ! Check for quads with three aligned faces and triangles with two aligned faces
        ccflags = c%cflags%GetAllElements()
        do ic = 1, c%ntot
            tf = GetCellFaceGA(c, ic)
            n_al = sum(facealigned(tf))
            nf = size(tf)

            if ((nf.eq.4).and.(n_al.eq.3)) then
                ! Find quad with more than two aligned faces
                ! Find triangle with more than one aligned face 
                ! Kick out the least aligned face via psi values
                fcs_al3 = tf(pack(ind4,facealigned(tf).eq.1));
                dpsi3 = abs(vpsi(v1n(fcs_al3))  - vpsi(v2n(fcs_al3)));
                indmax = maxloc(dpsi3,1);
                facealigned(fcs_al3(indmax)) = 0;

            elseif ((nf.eq.3) .and. (n_al.eq.2)) then
                if (ccflags(ic).eq.3) then

                    ! Make the boundary face non-algined
                    fb = tf(pack(ind3,b_flag(tf)));
                    facealigned(fb) = 0;     

                else
                    ! Kick out the least aligned face via psi values
                    fcs_al2 = tf(pack(ind3,facealigned(tf).eq.1));;
                    dpsi2 = abs(vpsi(v1n(fcs_al2))  - vpsi(v2n(fcs_al2)));
                    indmax = maxloc(dpsi2, 1);
                    facealigned(fcs_al2(indmax)) = 0;

                end if
            
            elseif (((nf.eq.4) .and. (n_al.eq.4)).or.((nf.eq.3).and.(n_al.eq.3))) then
                call gdErrorHandler("IdentifyAlignedFaces: not yet implemented")
            elseif (nf.gt.4) then
                call gdErrorHandler("IdentifyAlignedFaces: higher than quad is not supported")
            end if

        end do

        ! Override values for faces that are part of a flux surface (maybe put this upfront, or better not rely on correctness of fluxsurfaces????)

        if (allocated(grid%data%fluxdata%fluxsurfacefaces)) then
            if (grid%data%fluxdata%fluxsurfacefaces%Size() /= 0) then
                facealigned = 0
                call Unique(grid%data%fluxdata%fluxsurfacefaces%GetAllElements(), fcs)
                facealigned(fcs(2:size(fcs))) = 1  
            end if
        end if

        ! Save
        call f%aligned%SetAllElementsArray(facealigned)

        end associate

    end subroutine
    
    subroutine CheckUnstructuredGrid(grid,check_extra_conn)

        ! Description
        !============
        ! Check whether the connectivity of the grid is correct through 
        ! Multiple tests

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(in) :: grid
        logical :: check_extra_conn

        ! Auxiliary
        integer(I8), allocatable, dimension(:) :: cvLookUp, v1n, v2n, cf, &
            cellnumbers, verts_of_cell, cvertP2, cfaceP2, fcs, ccflags, nvxs, &
            b_cells, indCv
        integer(I8) :: verts_of_face(1:2), nc, nf, l, ic, nb, i, j, counter, &
            nface, nvert

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert, &
            fd => grid%data%fluxdata &
            )

            ! Initialize
            nface = c%faceP1%Get(c%ntot)+c%faceP2%Get(c%ntot)-1
            allocate(v1n(f%ntot), v2n(f%ntot), cf(nface), &
                cvertP2(c%ntot), cfaceP2(c%ntot), ccflags(c%ntot), &
                indCv(c%ntot))
            v1n = f%vert1%GetAllElements()
            v2n = f%vert2%GetAllElements()
            cf  = c%face%GetAllElements()
            cvertP2 = c%vertP2%GetAllElements()
            cfaceP2 = c%faceP2%GetAllElements()
            ccflags = c%cflags%GetAllElements()

            ! Check 1 - Equal lengths of data
            if ((c%vertP1%Size() /= c%ntot) .or. (c%faceP1%Size() /=  c%ntot) &
                    .or. (f%vert1%Size() /= f%ntot) .or. (v%x%Size() /= v%ntot)) then
                call gdErrorHandler('CheckUnstructuredGrid: check 1, length of data is not matching')
            end if

            ! Check2 -every cells has an equal amount of faces and vertices
            if (any(cfaceP2 /= cvertP2)) then
               call gdErrorHandler('CheckUnstructuredGrid: check 2, not every cell has equal amount of face and vertices')
            end if

            ! Check 3 - faces and vertices connection

            !voor alle faces, minstens 1 cell, 
            ! de twee vertices van het face moeten voorkomen in de verts van de (beide) cell

            !first make the cv_look_up
            nvert = c%vertP1%Get(c%ntot)+c%vertP2%Get(c%ntot)-1
            allocate(cvLookUp(nvert))
            cvLookUp = GetCvLookUpGA(c)

            do j = 1,f%ntot
                verts_of_face(1) = v1n(j)
                verts_of_face(2) = v2n(j)

                ! find a cell where the face is attached
                nc = count(cf == j)
                allocate(cellnumbers(nc))
                cellnumbers = pack(cvLookUp,cf == j)
                if (nc .gt. 2) then
                    call gdErrorHandler('CheckUnstructuredGrid: check 3, more than two // &
                     & cells connected to a face')
                end if

                !get vertices of cells
                do l = 1, nc
                    ic = cellnumbers(l);
                    verts_of_cell = GetCellVertGA(c,ic)
                    !verts_of_face moeten voorkomen in verts_of_cells
                    do i = 1,2
                        if (.not.any(verts_of_cell == verts_of_face(i))) then
                                call gdErrorHandler('CheckUnstructuredGrid: check 3, vertices // &
                                & of face not in vertices of the cells of that face')
                        end if
                    end do
                end do
                
                deallocate(cellnumbers)

            end do

            ! Check 4 - whether all internal cells have faces which are connected to two cells
            do ic = 1, c%ntot
                nf = cfaceP2(ic);
                fcs = GetCellFaceGA(c, ic)
                if (ccflags(ic) == 1) then
                    do i = 1, nf
                        if (.not.( size( pack(cvLookUp,cf == fcs(i)) ) .eq. 2 ) ) then
                            call gdErrorHandler('CheckUnstructuredGrid: check 4, not all internal faces have two cells')
                        end if

                    end do
                end if
            
            end do

            ! Check 5 duplicate face
            ![dup_faces,grid] = DetectOverlappingFaces(grid);

            !if ~isempty(dup_faces)
            !    fcX = 0.5*sum(vert.x(face.vert),2);
            !    fcY = 0.5*sum(vert.y(face.vert),2);
            !    figure,plotgeo_us(grid,'fast'),hold on
            !    plot(fcX(dup_faces),fcY(dup_faces),'g*')
            !    error('CheckUnstructuredGrid: Check 5, overlapping faces')
            !end
            ! Check 6 every vertex is only once in a flux surface
            if (allocated(fd%fluxsurfaceverts)) then
                if (fd%fluxsurfaceverts%Size() /= 0) then
                    call Unique(fd%fluxsurfaceverts%GetAllElements(), nvxs)
                    if (fd%fluxsurfaceverts%Size() /= size(nvxs)) then
                        !figure,plotgeo_us(grid,'fast'),hold on
                        !fs = grid.fs;
                        !for iFs = 1:fs.ntot
                        !    s = fs.vertP(iFs,1);
                        !    verts = fs.vert(s:s+fs.vertP(iFs,2)-1);
                        !    plot(grid.vert.x(verts),grid.vert.y(verts),'-*'),hold on;
                        !end
                        !for i = 1:length(fs.vert)
                        !    vx = fs.vert(i);
                        !    if sum(ismember(fs.vert,vx)) > 1
                        !        plot(grid.vert.x(vx),grid.vert.y(vx),'g*')
                        !         break
                        !    end
                        !end
                    
                        call gdErrorHandler('CheckUnstructuredGrid: check 6, vertex multiple times in fs.vert')

                    end if 
                end if
            end if

            ! Check 7 - wether all boundary cells have one faces which is only connected to that cell
            nb = count(ccflags .eq. 3)
            allocate(b_cells(nb))
            indCv = (/ (i, i = 1, c%ntot) /)
            b_cells = pack(indCv,ccflags .eq. 3) ! boundary cells
            do ic = 1, nb
                nf = cfaceP2(ic);
                fcs = GetCellFaceGA(c, ic)

                counter = 0;
                do j = 1,nf
                    if (size(pack(cvLookUp,cf == fcs(j)) ) .eq. 1 )  then
                        counter = counter + 1;
                    end if
                end do

                if ((counter .eq. 0) .or. (counter .eq. nf)) then!only one disconnected face == boundary face
                    !this statement is not closing as for al regular BC cells,
                    !count should be one and only for the corner cells it needs to
                    !be two
                    !figure, plotgeo_us(grid), hold on
                    !plot(cell.x(iCv),cell.y(iCv),'r*')
                    call gdErrorHandler('CheckUnstructuredGrid: check 7, boundary cell has no boundary face')
                end if
            end do

            ! Extra check on the extra connectivity fields
            if (check_extra_conn) then

                ! Not necessary for GAgrid I believe, and these field are not present anyway

                !%% Check 7
                !if any(face.cellP(:,2)==0)
                !    flags(7) = 1;
                !    error('CheckUnstructuredGrid: check 7, face(s) with zero cells')
                !end

                !%% Check 8
                !if any(vert.cellP(:,2)==0)
                !    flags(8) = 1;
                !    figure,plotgeo_us(grid,'fast'),hold on
                !    plot(vert.x(vert.cellP(:,2)==0),vert.y(vert.cellP(:,2)==0),'g*')
                !    error('CheckUnstructuredGrid: check 8, vertex with zero cells')
                    
                !end

                !%% Check 9
                !if any(vert.faceP(:,2)==0)
                !    flags(9) = 1;
                !        figure,plotgeo_us(grid,'fast'),hold on
                !    plot(vert.x(vert.cellP(:,2)==0),vert.y(vert.cellP(:,2)==0),'g*')
                !    error('CheckUnstructuredGrid: check 9, vertex with zero faces')
                !end


            end if 

        end associate


    end subroutine

    subroutine RecalcMagn(grid, magneticField)

        ! Description
        !============
        ! Evaluates the magneticField interpolants

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout) :: grid
        type(MagneticFieldUDT)          :: magneticField

        ! Auxiliary
        real(R8) :: dpsidx(grid%vert%ntot), dpsidy(grid%vert%ntot), &
            psi_cv(grid%cell%ntot), psi_vx(grid%vert%ntot),  &
            vx(grid%vert%ntot), vy(grid%vert%ntot), cx(grid%cell%ntot), &
            cy(grid%cell%ntot)

        ! Initialize
        vx = grid%vert%x%GetAllElements()
        vy = grid%vert%y%GetAllElements()
        cx = grid%cell%x%GetAllElements()
        cy = grid%cell%y%GetAllElements()

        ! Evaluate magneticField
        call magneticField%interp%Evaluate(cx,cy,0,0,psi_cv)
        call magneticField%interp%Evaluate(vx,vy,0,0,psi_vx)
        call magneticField%interp%Evaluate(vx,vy,1,0,dpsidx)
        call magneticField%interp%Evaluate(vx,vy,0,1,dpsidy) 

        call grid%cell%psi%SetAllElementsArray(psi_cv) 
        call grid%vert%psi%SetAllElementsArray(psi_vx) 
        call grid%vert%bx%SetAllElementsArray(dpsidx) 
        call grid%vert%by%SetAllElementsArray(dpsidy) 

    end subroutine

    subroutine MergeFS(grid)

        ! Description
        !============
        ! Checks whether flux surfaces should be merged to one and performs the
!        mergeing.   

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT) :: grid

        ! Auxiliary
        integer(I8), allocatable ::  mFS(:,:), mFS_update(:,:)
        integer(I8), allocatable, dimension(:) :: fieldlineID, tfv, ar, tv1, &
            tf1, tv2, tf2, new_faces, new_verts, rem_ind_v, rem_ind_f, fID, &
            hasID, range
        integer(I8) :: i, j, ifs, counter, ifs1, ifs2, nv1, nf1, nv2, &
            nf2, counterv, counterf

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert, &
            fd => grid%data%fluxdata &
            )

        !Initialize
        allocate(mFS(fd%nFs,2))
        allocate(hasID(v%ntot))
        fieldlineID = 0
        hasID = 0
        mFS = 0
        counter = 0

        ! Augment to mergeFS and mergeFSP
        do ifs = 1, fd%nFs

            ! Get vertices of flux surface
            tfv = GetFluxSurfaceVxsGA(fd, ifs)

            ! Check if already attributed
            if (any(hasID(tfv).eq.1)) then

                ! vertices detected that belong to multiple flux surfaces
                ! Indicate the surfaces to merge
                fID = fieldlineID(tfv)
                ar = fID(hasID(tfv))
                ifs1 = ar(1)
                if (size(ar) .gt. 1) &
                    print *, 'Warning: MergeFS: merging of multiple flux surfaces not yet supported'

                counter = counter + 1
                mFS(counter,1) = ifs1
                mFS(counter,2) = ifs
            end if

            ! Set the fieldline ID
            fieldlineID(tfv) = ifs
            hasID(tfv) = 1

        end do

        ! Merge indicated flux surfaces - augment to n-fluxsurface - TODO
        do i = 1, counter

            ! Get flux surface indices
            ifs1 = mFS(i,1)
            ifs2 = mFS(i,2)

            ! Get vertices and faces
            tv1 = GetFluxSurfaceVxsGA(fd, ifs1)
            nv1 = fd%fluxsurfacevertsP2%Get(ifs1)
            tf1 = GetFluxSurfaceFcsGA(fd, ifs1)
            nf1 = fd%fluxsurfacefacesP2%Get(ifs1)
            tv2 = GetFluxSurfaceVxsGA(fd, ifs2)
            nv2 = fd%fluxsurfacevertsP2%Get(ifs2)
            tf2 = GetFluxSurfaceFcsGA(fd, ifs2)
            nf2 = fd%fluxsurfacefacesP2%Get(ifs2)

            ! Remove both of them and add one at the back
            ! Get the indices to remove out of fsVx and fsFc
            allocate(rem_ind_v(v%ntot))
            allocate(rem_ind_f(f%ntot))
            rem_ind_v = 0
            rem_ind_f = 0

            counterv = 0
            counterf = 0

            ! First ifs
            rem_ind_v(counterv + 1: counterv+nv1) = (/ (j, j = fd%fluxsurfacevertsP1%Get(ifs1),&
                fd%fluxsurfacevertsP1%Get(ifs1)+nv1-1) /)
            counterv = counterv + nv1
            rem_ind_f(counterf + 1: counterf+nf1) = (/ (j, j = fd%fluxsurfacefacesP1%Get(ifs1),&
                fd%fluxsurfacefacesP1%Get(ifs1)+nf1-1) /)
            counterf = counterf + nf1

            ! Second ifs
            rem_ind_v(counterv + 1: counterv+nv2) = (/ (j, j = fd%fluxsurfacevertsP1%Get(ifs2),&
                fd%fluxsurfacevertsP1%Get(ifs2)+nv2-1) /)
            counterv = counterv + nv2
            rem_ind_f(counterf + 1: counterf+nf2) = (/ (j, j = fd%fluxsurfacefacesP1%Get(ifs2),&
                fd%fluxsurfacefacesP1%Get(ifs2)+nf2-1) /)
            counterf = counterf + nf2

            ! Remove in fsVx and fsFc
            call fd%fluxsurfaceverts%RemoveMultipleElements(rem_ind_v(1:counterv))
            call fd%fluxsurfacefaces%RemoveMultipleElements(rem_ind_f(1:counterf))

            ! Remove and adjust pointer arrays
            call fd%fluxsurfacevertsP1%RemoveSingleElement(ifs1)
            call fd%fluxsurfacevertsP2%RemoveSingleElement(ifs1)
            call fd%fluxsurfacefacesP1%RemoveSingleElement(ifs1)
            call fd%fluxsurfacefacesP2%RemoveSingleElement(ifs1)
            fd%nFs = fd%nFs - 1

            range = (/ (j, j = ifs1, fd%nFs)/)
            call fd%fluxsurfacevertsP1%SetMultipleElements(range, fd%fluxsurfacevertsP1%GetMultipleElements(range) - nv1)
            call fd%fluxsurfacefacesP1%SetMultipleElements(range, fd%fluxsurfacefacesP1%GetMultipleElements(range) - nf1)

            ifs2 = ifs2 - 1
            call fd%fluxsurfacevertsP1%RemoveSingleElement(ifs2)
            call fd%fluxsurfacevertsP2%RemoveSingleElement(ifs2)
            call fd%fluxsurfacefacesP1%RemoveSingleElement(ifs2)
            call fd%fluxsurfacefacesP2%RemoveSingleElement(ifs2)
            fd%nFs = fd%nFs - 1

            range = (/ (j, j = ifs2, fd%nFs)/)
            call fd%fluxsurfacevertsP1%SetMultipleElements(range,fd%fluxsurfacevertsP1%GetMultipleElements(range) - nv2)
            call fd%fluxsurfacefacesP1%SetMultipleElements(range,fd%fluxsurfacefacesP1%GetMultipleElements(range) - nf2)

            ! Add new flux surface
            fd%nFs = fd%nFs + 1
            call Unique([tv1, tv2], new_verts)
            call Unique([tf1, tf2], new_faces)

            call fd%fluxsurfacevertsP1%AppendSingleElement(fd%fluxsurfacevertsP1%Get(fd%nFs-1)+ &
                fd%fluxsurfacevertsP2%Get(fd%nFs-1))
            call fd%fluxsurfacevertsP2%AppendSingleElement(size(new_verts))
            call fd%fluxsurfaceverts%AppendMultipleElements(new_verts)
            call fd%fluxsurfacefacesP1%AppendSingleElement(fd%fluxsurfacefacesP1%Get(fd%nFs-1)+ &
                fd%fluxsurfacefacesP2%Get(fd%nFs-1))
            call fd%fluxsurfacefacesP2%AppendSingleElement(size(new_faces))
            call fd%fluxsurfacefaces%AppendMultipleElements(new_faces)

            ! Update other flux surface indices in mFS
            allocate(mFS_update(size(mFS(:,1)),2))
            mFS_update = 0
            do j = i+1, counter

                if (mFS(j,1) .gt. ifs1) mFS_update(j,1) = mFS_update(j,1) + 1 
                if (mFS(j,2) .gt. ifs1) mFS_update(j,2) = mFS_update(j,2) + 1 
                if (mFS(j,1) .gt. ifs2+1) mFS_update(j,1) = mFS_update(j,1) + 1 
                if (mFS(j,2) .gt. ifs2+1) mFS_update(j,2) = mFS_update(j,2) + 1 

            end do

            ! Update
            mFS = mFs - mFS_update
            deallocate(mFS_update)
            deallocate(rem_ind_f)
            deallocate(rem_ind_v)

        end do

        end associate

        
        
        
    end subroutine

    !------------------------------------------------------------------!
    !                        GAFACE ROUTINES                           !
    !------------------------------------------------------------------!    

    subroutine ChainFaces(face, f_list, f_ord, nf)

        ! Description
        !============
        ! Order a list of faces by chaining then together. The faces need to be an open chain
        ! - f_list  : input list of faces
        ! - f_ord   : ordered list of faces
        ! - nf      : number of faces per detected chains

        ! Declare variables
        !==================
        ! Arguments
        class(GAFaceUDT), intent(in)            :: face
        integer(I8), allocatable, intent(in)    :: f_list(:)
        integer(I8), allocatable, intent(out)   :: f_ord(:,:), nf(:)

        ! Auxiliary
        integer(I8), allocatable :: verts_list(:,:), indv(:,:), &
            a(:)
        integer(I8) :: i, j, k, l, ends, nfcs, v1n(face%ntot), &
            v2n(face%ntot), start_fcs, com_vert(1:2), na, f

        ! Initialize
        nfcs = size(f_list)
        allocate(verts_list(nfcs,2), indv(nfcs,2))
        v1n = face%vert1%GetAllElements()
        v2n = face%vert2%GetAllElements()
        verts_list(:,1) = v1n(f_list)
        verts_list(:,2) = v2n(f_list)
        ends = 0
        indv(:,1) = (/(i, i = 1, nfcs)/)
        indv(:,2) = (/(i, i = 1, nfcs)/)

        ! Check the occurence of vertices to determine the starting and endings of a chain, and save one end face for algo efficienvy
        do i = 1, nfcs
            do j = 1, 2
                if (count(verts_list == verts_list(i,j)) .eq. 1) then
                    start_fcs = f_list(i)
                    ends = ends + 1 
                end if
            end do
        end do

        ! Get correct sizes for allocatable output arrays
        if (ends .gt. 2) then
            allocate(f_ord(nfcs, ends/2))
            f_ord = 0
            allocate(nf(ends/2))
            nf = 0
        else
            allocate(f_ord(nfcs,1))
            f_ord = 0
            allocate(nf(1))
            nf = size(f_list)
        end if

        ! Construct the chaines
        do k = 1, ends/2

            ! Determine the first face
            if (k /= 1) then 
                do i = 1, nfcs
                    do j = 1, 2 
                        if (count(verts_list == verts_list(i,j)) .eq. 1) then
                            f = f_list(i)
                            if (.not.any(f_ord == f)) then
                                f_ord(1, k) = f_list(i)
                            end if
                        end if
                    end do
                end do

            else 
                f_ord(1, 1) = start_fcs
            end if

            ! Get next faces
            do i = 2, nfcs

                ! Check whether to continue
                if (f_ord(i-1,k) == 0) then
                    nf(k) = i-2
                    exit
                end if

                com_vert(1) = v1n(f_ord(i-1,k))
                com_vert(2) = v2n(f_ord(i-1,k))
                do l = 1, 2 ! Check two vertices

                    na = count(verts_list == com_vert(l))
                    allocate(a(na))
                    a = pack(indv,verts_list == com_vert(l))
                    
                    ! Check if vertex is end of chain or not
                    if (na .eq. 2) then

                        ! If one of the faces is not in the list, add it
                        if (.not.any(f_list(a(1)) == f_ord(:,k))) then
                            f_ord(i,k) = f_list(a(1))
                        elseif (.not.any(f_list(a(2)) == f_ord(:,k))) then
                            f_ord(i,k) = f_list(a(2))
                        end if

                    elseif (na .eq. 1) then

                        if (.not.any(f_list(a) == f_ord(:,k))) then
                            f_ord(i,k) = f_list(a(1))
                        end if

                    else 
                        call gdErrorHandler('ChainFaces: Something went wont, probably faces can not form a chain')
                    end if

                    deallocate(a)
                end do

            end do

        end do

    end subroutine

    !------------------------------------------------------------------!
    !                        DISTANCE FUNCTIONS                        !
    !------------------------------------------------------------------! 
    subroutine ComputeDistanceFunction(DistFunc, grid, options, base_func)

        ! Description
        !============
        ! Construct a distance function to use as criterium for splitting and merging. 
        !The base function is normal distribution b = exp(-1/d (dist)) with d the chararcteristic length.

        ! Inputs
        !======
        ! - dist_type: type of distance function
        ! - d_rescale: 
        ! - d_char_type: option to chose characteristic length 
        ! - base_func: base function for distance calculation

        ! Declare variables
        !==================
        ! Arguments 
        class(DistFunctUDT), intent(inout)      :: DistFunc
        type(GAGridUDT), intent(inout)          :: grid 
        type(GAoptionsUDT), intent(in)          :: options
        character(:), allocatable, intent(in)   :: base_func

        ! Auxiliary
        integer(I8), allocatable, dimension(:) :: fcs_t1,  &
            fcs1, fcs, fclbl_fcs1, indFc, v1n, v2n, cvLookUp, fcLbL_loc, & 
            vxs_t1, faces_sep, vxs_sep, vx1, vx2
        real(R8), allocatable, dimension(:) :: b, dFun, fcX, fcY, &
            vx, vy, vxFun, vyFun, distX, fcLength, vertsX, vertsY
        integer(I8) :: i, nv, nf, resx, resy
        real(R8) :: d, xx, xy, SOLwidth, xrange(1:2), yrange(1:2), valinf
        logical :: use_nsep, use_sepID, start
        character(:), allocatable :: savefilepath

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert, &
            data => grid%data &
            )

            ! Initialize
            allocate(fcX(f%ntot),fcY(f%ntot), vx(v%ntot), vy(v%ntot), &
                indFc(f%ntot), fcLbL_loc(f%ntot))
            v1n = f%vert1%GetAllElements()
            v2n = f%vert2%GetAllElements()
            vx  = v%x%GetAllElements()
            vy  = v%y%GetAllElements()

            ! Get GA facelabels
            fcLbl_loc = GetfcLblGA(f,options)

            ! Determine charateristic decay length
            !=====================================
            ! Take the target vertices
            indFc = (/ (i, i = 1, f%ntot)/)
            fcs_t1 = pack(indFc, fcLbl_loc.ge.4 )

            vxs_t1 = GetVxsFromFcsGA(f,fcs_t1)
            nv = size(vxs_t1)

            ! Get vertex coordinates
            vertsX = vx(vxs_t1)
            vertsY = vy(vxs_t1)

            cvLookUp = GetCvLookUpGA(c)
            use_nsep = .true.
            use_sepID = .true. 
            start = .true.
            call grid%GiveSeparatrices(use_nsep,use_sepID,start,cvLookUp)
            call grid%GiveXpoints(use_nsep,cvLookUp)
            allocate(distX(1:data%nxp))
            distX = 0 

            select case (DistFunc%d_char_type)

            case ('min_Xpoint_dist')

                do i = 1, data%nxp
                    xx = vx(data%xpointID(i))
                    xy = vy(data%xpointID(i))
                    distX(i) = minval(sqrt((vertsX - xx)**2 + (vertsY - xy)**2))
                end do

            case ('max_Xpoint_dist')

                do i = 1, data%nxp
                    xx = vx(data%xpointID(i))
                    xy = vy(data%xpointID(i))
                    distX(i) = maxval(sqrt((vertsX - xx)**2 + (vertsY - xy)**2))
                end do 
            case default

                call gdErrorHandler("ComputeDistanceFunction: d_char_type not implemented")
                
            end select

            ! Rescale characteristic distance
            d = minval(distX) / DistFunc%d_rescale

            ! Chose type of distance function
            select case (DistFunc%dist_type)

            case ('target_single_null')
                ! Reuse target vertices
                vXfun = vertsX
                vYfun = vertsY

                ! Characteristic length
                allocate(dFun(1:nv))
                dFun = d 

            case ('target_to_vessel')

                ! Take boundary faces based on face labels
                allocate(fcs1(count(fcLbl_loc /= 0)))
                fcs1 = pack(indFc,fcLbl_loc /= 0)

                fclbl_fcs1 = fcLbl_loc(fcs1)

                nf = count(fclbl_fcs1 /= 2)
                allocate(fcs(nf))
                fcs = pack(fcs1,fclbl_fcs1 /= 2)

                fcX = 0.5_R8*(vx(v1n)+vy(v2n))
                fcY = 0.5_R8*(vy(v1n)+vy(v2n))

                allocate(vxFun(nf), vyFun(nf))
                vxFun = fcX(fcs)
                vyFun = fcY(fcs)

                ! Compute length of faces
                allocate(fcLength(nf))
                fcLength = sqrt( (vx(v2n(fcs)) - vx(v1n(fcs)))**2 &
                          + (vy(v2n(fcs)) - vy(v1n(fcs)))**2 )
                
                allocate(dFun(1:nf))
                dFun = sum(fcLength)*2/(real(nf, kind=R8))
                nv = nf

                print *, 'Warning: ComputeDistanceFunction:  problem with target to vessel option ' // &
                    'that the distance function is lower where boundary faces are larger. Should be' // &
                    'changed to distance function based on line segments.' 


            
            case ('pol_flux_est')

                print *, 'Warning: ComputeDistanceFunction: dist_type "pol_flux_est" only implemented for SN-case' // &
                         'Other wise the assumption is made that sepID(1) is the main separatrix.'
                 
                ! Get faces of separatrix
                faces_sep = GetFluxSurfaceFcsGA(data%fluxdata,data%sepID(1))

                vxs_sep = GetVxsFromFcsGA(f, faces_sep)
                nv = size(vxs_sep)

                ! Refinement if needed
                vx1 = v1n(faces_sep)
                vx2 = v2n(faces_sep)
                fcLength = sqrt( (vx(v2n(faces_sep)) - vx(v1n(faces_sep)))**2 &
                               + (vy(v2n(faces_sep)) - vy(v1n(faces_sep)))**2)
                               
                ! Vertex coordinates
                allocate(vxFun(1:nv), vyFun(1:nv), dFun(1:nv))
                vxFun = vx(vxs_sep)
                vyFun = vy(vxs_sep)

                ! Very rough calculation solwidth - could be improved by computing distance between min-max positions
                ! But good enough if the core is large and divertor legs are short
                SOLwidth = 0.25_R8 * ( (maxval(vx)- maxval(vxFun)) + &
                            (minval(vx - minval(vxFun))) + (maxval(vy) - maxval(vyFun)) & 
                            + (minval(vy) - minval(vyFun))) 

                dFun = SOLwidth

            case ('separatrix')

                call gdErrorHandler("ComputeDistanceFunction: dist_type 'separatrix' not implemented")

            case default

                call gdErrorHandler("ComputeDistanceFunction: dist_type not implemented")
            
            end select

            ! Desired values at vertices
            allocate(b(1:nv))
            b = 1.0_R8

            ! Desired values at infinity
            valinf = 0.0_R8


            ! Initialize distance function
            call DistFunc%distr%Initialize(vxFun,vyFun,b,valinf,dFun)

            ! Visualize
            if (options%plt_dist_func) then
                xrange(1) = minval(vx)
                xrange(2) = maxval(vx) 
                yrange(1) = minval(vy)
                yrange(2) = maxval(vy) 
                resx = nint(xrange(2) - xrange(1)*100)
                resy = nint(yrange(2) - yrange(1)*100)
                savefilepath = DistFunc%dist_type
                call DistFunc%distr%Visualize(xrange, yrange, resx, resy, savefilepath)
            end if
                

        end associate

    end subroutine



    !==================================================================!
    !                                                                  !
    !                           FUNCTIONS                              !
    !                                                                  !
    !==================================================================! 
    
    ! Get vertices of a cell and dynamic arrays
    function GetCellVertGA(cell, i) result(res)
        integer(I8)         :: i, j, s
        type(GACellUDT)     :: cell
        integer(I8), allocatable :: res(:), range(:)
        s = cell%vertP1%Get(i)
        range = (/ (j, j = s, (s + cell%vertP2%Get(i) - 1)) /)
        res = cell%vert%GetMultipleElements(range)   
    end function

    ! Get vertices of a cell and dynamic arrays
    function GetCellFaceGA(cell, i) result(res)
        integer(I8)         :: i, j, s 
        type(GACellUDT)     :: cell
        integer(I8), allocatable :: res(:), range(:)
        s = cell%faceP1%Get(i)        
        range = (/ (j, j = s, (s + cell%faceP2%Get(i) - 1)) /)
        res = cell%face%GetMultipleElements(range)   
    end function
    ! Get cells of a face with dynamic arrays
    function GetFaceCellGA(cell, i, cvLookUp) result(res)
        integer(I8)                 :: i 
        type(GACellUDT)               :: cell
        integer(I8), allocatable, optional    :: cvLookUp(:)
        integer(I8), allocatable    :: res(:)  
        
        if (.not.present(cvLookUp)) then
            cvLookUp = GetCvLookUpGA(cell)
        end if

        res = pack(cvLookUp,cell%face%GetAllElements().eq.i)
    end function

    function GetFluxSurfaceFcsGA(fd, i) result(res)
        integer(I8) :: i, s, nf, j
        type(GAFluxDataUDT) :: fd
        integer(I8), allocatable :: res(:), range(:)

        nf = fd%fluxsurfacefacesP2%Get(i)
        s = fd%fluxsurfacefacesP1%Get(i)
        range = (/(j, j = s, (s+nf-1) )/)
        res = fd%fluxsurfacefaces%GetMultipleElements(range)

    end function   
    
    function GetFluxSurfaceVxsGA(fd, i) result(res)
        integer(I8) :: i, s, nf, j
        type(GAFluxDataUDT) :: fd
        integer(I8), allocatable :: res(:), range(:)

        nf = fd%fluxsurfacevertsP2%Get(i)
        s = fd%fluxsurfacevertsP1%Get(i)
        range = (/(j, j = s, (s+nf-1) )/)
        res = fd%fluxsurfaceverts%GetMultipleElements(range)
    end function

    ! Get cells of a vertex without using vert%cell and dynamic arrays
    function GetVertCellGA(cell, i, cvLookUp) result(res)
        integer(I8)                 :: i 
        type(GACellUDT)               :: cell
        integer(I8), allocatable, optional    :: cvLookUp(:)
        integer(I8), allocatable    :: res(:)
        
        if (.not.present(cvLookUp)) then
            cvLookUp = GetCvLookUpGA(cell)
        end if
            
        res = pack(cvLookUp,cell%vert%GetAllElements().eq.i)
    end function

    function GetCvLookUpGA(cell) result(res)
        type(GACellUDT)       :: cell
        integer(I8)         :: nc, ic, nv, s, i               
        integer(I8), allocatable :: res(:), range(:)

        nc = cell%ntot
        range = (/ (i, i = 1,(cell%vertP1%Get(nc)+cell%vertP2%Get(nc)-1))/)

        allocate(res(1:cell%vertP1%Get(nc)+cell%vertP2%Get(nc)-1))
        res = 0

        do ic = 1, nc
            s = cell%vertP1%Get(ic)
            nv = cell%vertP2%Get(ic)
            range = (/ (i, i = s, (s+nv-1)) /)
            res(range) = ic
        end do
    end function

    function GetfcLblGA(f,options) result(res)
        type(GAFaceUDT) :: f
        type(GAoptionsUDT) :: options
        integer(I8) :: res(1:f%ntot), indFc(1:f%ntot), i 

        res = f%label%GetAllElements()
        indFc = (/ (i, i=1,f%ntot) /)
        do i = 1, size(options%facelabelmappingGG)
            res(pack(indFc, f%label%GetAllElements() == options%facelabelmappingGG(i))) &
                = options%facelabelmappingGA(i)
        end do

    end function

    function GetVxsFromFcsGA(f,fcs) result(res)
        type(GAFaceUDT) :: f
        integer(I8), allocatable :: fcs(:), verts(:), res(:)
        integer(I8) :: nf

        nf = size(fcs)
        allocate(verts(1:nf*2))
        verts = 0
        verts(1:nf) = f%vert1%GetMultipleElements(fcs)
        verts(nf+1:nf*2) = f%vert2%GetMultipleElements(fcs)
        call Unique(verts, res)

    end function

    function DetectSepID(ifs, cell, faces, sepIDloc, nf, step, cvLookUp, counter) result(res)
        ! Description
        !============
        ! Check whethere is a core region on one side of the flux surface
        integer(I8) :: ifs, nf, counter, i, step, ifc, reg1, reg2, sepIDloc(1:100)
        type(GACellUDT) ::  cell
        integer(I8), allocatable :: cvLookUp(:), faces(:), cvs(:), res(:), creg(:)

        ! Initialize
        res = sepIDloc
        creg = cell%reg%GetAllElements()

        do i = 1, nf, step
            ifc = faces(i)
            cvs = GetFaceCellGA(cell,ifc, cvLookUp)
            if (size(cvs).eq.2) then
                reg1 = creg(cvs(1))
                reg2 = creg(cvs(2))
                if (reg1 /= reg2) then
                    if ((mod(reg1,4).eq.1) .or. (mod(reg2,4).eq.1)) then 
                        counter = counter + 1
                        res(counter) = ifs
                        exit
                    end if
                end if
            end if
        end do


        
    end function


end module 