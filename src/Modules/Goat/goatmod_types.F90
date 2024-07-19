!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains all the type definitions used in the goat. 
! See the option setting routines in goatmod_userinput
! for an explanation of the different fields and options. Additionally,
! precision types are set here that can (should) be used in all 
! subroutines. 

! Note: default accessibility is set to public, since all these types 
! and corresponding functions are free to be used. 

module goatmod_types

    ! Initialize
    !============
    use mod_precision
    use mod_polygon
    use Interpolant
    use PolygonShapeFunction
    use PolygonLevelsetFunction2D

    ! The usual
    implicit none
    save
    public 

    ! All private types
    
    ! All private functions

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                               Grid                               !
    !------------------------------------------------------------------!

    ! Main grid substructures
    !========================
    ! Vertex structure
    type VertexUDT

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
        real(R8), allocatable               :: x(:), y(:)

        ! Logicals and indices
        logical, allocatable                :: BV(:)
        integer(I8), allocatable            :: fieldlineID(:) 
        integer(I8)                         :: ntot = 0

        integer(I8), allocatable            :: faceP(:,:)
        integer(I8), allocatable            :: face(:)
        integer(I8)                         :: nface = 0

        integer(I8), allocatable            :: cellP(:,:)
        integer(I8), allocatable            :: cell(:)
        integer(I8)                         :: ncell = 0

        integer(I8), allocatable            :: neigP(:,:)
        integer(I8), allocatable            :: neig(:)
        integer(I8)                         :: nneig = 0

        ! Other data
        real(R8), allocatable, dimension(:) :: bx, by, psi, ffbz 

    end type

    ! Face structure
    type FaceUDT

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
        ! - reg             : face regions
        ! - aligned         : integer that is 1 if the face is aligned


        ! Logicals and indices
        integer(I8), allocatable            :: vert(:,:)

        integer(I8), allocatable            :: cell(:), label(:), &
            reg(:), aligned(:)
        integer(I8), allocatable            :: cellP(:, :)
        integer(I8)                         :: ncell = 0

        integer(I8)                         :: ntot = 0
        logical, allocatable                :: BF(:)

    end type

    ! Cell structure
    type CellUDT

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
        integer(I8), allocatable            :: vertP(:,:)
        integer(I8), allocatable            :: vert(:)
        integer(I8)                         :: nvert = 0

        integer(I8), allocatable            :: faceP(:,:)
        integer(I8), allocatable            :: face(:)
        integer(I8)                         :: nface = 0

        logical, allocatable                :: GC(:)

        integer(I8)                         :: ntot = 0, ngc
        
        real(R8), allocatable, dimension(:) :: psi, bp, bt, x, y
        integer(I8), allocatable, dimension(:)  :: cflags(:), reg(:), &
            ft(:)
    end type

    ! Boundary structure
    type BndUDT

        ! Description
        !============
        ! Fields:   
        ! - nvert:                      total number of vertices
        ! - vert:                       all vertices belonging to that 
        !                               boundary
        ! - nface:                      total number of faces 
        ! - face:                       all face indices 
        ! - ID:                         boundary ID (see below) that 
        !                               identifies the boundary type.
        ! 
        ! See also the interface routine InterfaceBoundaryMapping in 
        ! gdmod_interfaces on how the boundary IDs are mapped given 
        ! input from the grid generator. 

        ! Note: the amount of vertices is per definition the amount of
        ! faces minus one, also for closed boundaries. In the latter 
        ! case, the begin and end vertex index are the same. This 
        ! assumes that the boundaries should be polygons that either 
        ! close perfectly on themselves, or are simple polygons (i.e. 
        ! no polygon vertex should have more than three faces it belongs
        ! to).

        ! Boundary ID meaning
        !====================
        ! GRID DEFORMATION              PHYSICAL MEANING
        !
        ! 1                             target plate (inner)
        ! 2                             target plate (outer)
        ! 3                             private flux
        ! 4                             core boundary
        ! 5                             outermost flux surf.


        ! Boundary vertices
        integer(I8)                         :: nvert = 0 ! simply nfaces-1 actually
        integer(I8), allocatable            :: vert(:)   

        ! Boundary faces
        integer(I8)                         :: nface = 0
        integer(I8), allocatable            :: face(:)             
        
        ! Boundary ID
        integer(I8)                         :: ID

    end type

    ! Additional grid data structures
    !================================
    ! Used to collect derived grid data, e.g. which link the grid cells 
    ! flux tubes. Other output data from the grid generator that does 
    ! not fit within the classical cell/faces/vert structures can be 
    ! added as well through user defined types. 

    ! Flux data
    type FluxDataUDT

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
        integer(I8), allocatable            :: fluxtubecellsP(:,:)
        integer(I8), allocatable            :: fluxtubecells(:)
        integer(I8), allocatable            :: fluxtubefacesP(:,:)
        integer(I8), allocatable            :: fluxtubefaces(:)

        ! Arrays, flux surface data
        integer(I8), allocatable            :: fluxsurfacefacesP(:,:)
        integer(I8), allocatable            :: fluxsurfacefaces(:)
        integer(I8), allocatable            :: fluxsurfaceID(:)
        real(R8), allocatable               :: fluxsurfacepsi(:)


    end type

    ! Region data
    type RegionDataUDT

        ! Description
        !============
        ! Data type to collect all information on which cells/verts/face
        ! belongs to which grid region. 
        ! Fields:
        !
        ! - fluxtuberegID       : fluxdata%nFt-by-1 array containing 
        !                       the region IDs for each flux tube

        ! Arrays
        integer(I8), allocatable            :: fluxtuberegID(:)

    end type

    ! Structured grid data (to be removed in the future)
    type StructuredGridDataUDT

        ! Description
        !============
        ! Data type to collect all data related to a possibly initial 
        ! structured grid that served as basis for the current 
        ! unstructured grid data. Only saved for backward compatibility 
        ! reasons. Hopefully deleted in the future. 
        ! Fields:
        !
        ! - isClassicalGrid         : integer indicating whether the
        !                           whether the grid was a classic 
        !                           structured grid
        ! - nx, ny                  : dimensions of the original 
        !                           structured grid
        ! - nncut                   : number of cuts

        ! Logicals & scalars
        integer(I8)             :: nx, ny, nncut
        integer(I4)             :: isClassicalGrid

    end type

    ! Grid data 
    type GridDataUDT

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
        ! - regions             : UDT with all data to which region
        !                       cells, faces, ... belong
        ! - sglegacy            : data from legacy structured grids
        ! - OMPcell, OMPface    : cells and faces belonging to outer mid
        !                       plane
        ! - IMPcell, IMPface    : same, but inner mid plane
        ! - OMPr, OMPz          : points defining line segment of OMP
        ! - IMPr, IMPz          : same but for inner mid plane
        ! - xpointID            : array containing all X-point vertex IDs
        ! - nxp                 : number of x-points

        ! Flux data
        type(FluxDataUDT)           :: fluxdata

        ! Region data
        type(RegionDataUDT)         :: regions

        ! Legacy data of structured grid
        type(StructuredGridDataUDT) :: sglegacy

        ! OMP & IMP
        integer(I8), allocatable, dimension(:)  :: OMPcell, OMPface, &
            IMPcell, IMPface
        integer(I8)                             :: nOMPcell, nOMPface, &
            nIMPcell, nIMPface
        real(R8), dimension(1:2)                :: OMPr, OMPz, IMPr, &
            IMPz

        ! X-point(s)
        integer(I8), allocatable                :: xpointID(:)
        integer(I8)                             :: nxp

    end type

    ! Main grid structure
    !====================
    type GridUDT

        ! Description
        !============
        ! Data structure containing all the grid data substructures. 
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
        type(VertexUDT)                     :: vert

        ! Faces
        type(FaceUDT)                       :: face

        ! Cells
        type(CellUDT)                       :: cell

        ! Additional data
        type(GridDataUDT)                   :: data

        ! Boundaries
        type(BndUDT), allocatable           :: bnd(:)

    end type

    !------------------------------------------------------------------!
    !                               Numerics                           !
    !------------------------------------------------------------------!

    ! Interpolation
    !==============

    !------------------------------------------------------------------!
    !                           Magnetic field                         !
    !------------------------------------------------------------------!

    ! Magnetic field
    !===============
    type MagneticFieldUDT

        ! Description
        !============
        ! Data structure containing all the magnetic field data. The 
        ! following fields should be present:
        !
        ! - interpolant:    see type definition for explanation
        ! - nR:             number of sample points in the R-direction
        ! - nZ:             number of sample points in the Z-direction
        ! - R:              nR-by-1 array with R-coordinates
        ! - Z:              nZ-by-1 array with Z-coordinates
        ! - Psi:            nR-by-nZ array with magnetic flux values

        ! Coordinates
        integer(I8)                 :: nR = 0 
        integer(I8)                 :: nZ = 0
        real(R8), allocatable       :: R(:)
        real(R8), allocatable       :: Z(:)
        real(R8), allocatable       :: Psi(:,:)
        
        ! Interpolant
        type(StructuredInterpolant2DUDT)    :: interp

    end type

    !------------------------------------------------------------------!
    !                            Environment                           !
    !------------------------------------------------------------------!
    ! Vessel structures
    !==================
    type VesselStructureUDT

        ! Description
        !============
        ! User defined type for the different substructures of the 
        ! vessel. These are actually polygons (so with x, y coordinates)
        ! that are used to build up the vessel structure afterwards. 
        ! The following fields are present:
        !
        ! - np:         number of points in the polygon
        ! - x, y:       polygon coordinates
        ! - isclosed:   logical to indicate if the polygon should close
        !               upon itself
        ! - id:         identifier (integer number)

        ! Coordinates
        integer(I8)                         :: np = 0
        real(R8), allocatable               :: x(:), y(:)

        ! Logicals
        logical                             :: isclosed

        ! ID
        integer(I4)                         :: ID

    end type

    ! Vessel
    !=======
    type VesselUDT
        
        ! Description
        !============
        ! User defined structure to keep track of all vessel related 
        ! data that is required for the grid optimization. Currently, 
        ! this includes the following structures and fields:
        !
        ! - np:         number of points in the boundary polygon
        ! - x, y:       x, y-coordinates (np-by-1 array)
        ! - TPind:      integer array (np-by-1) indicating on which TP
        !               a node is located
        ! - ntp:        number of all target plate indices
        ! - allTPind:   a list (ntp-by-1) of all target plate indices
        ! - geom:       the vessel geometry being considered (character)
        ! - nstructures:number of sub-structures of the vessel
        ! - structures : nstructures-by-1 array of structures
        ! - targetpolygons:     closed polygon representation of each 
        !                       vessel target plate

        ! Target plates
        integer(I4)                         :: ntp = 0
        integer(I4), allocatable            :: allTPind(:)
        integer(I4), allocatable            :: TPind(:)

        ! Structures
        integer(I4)                         :: nstructures = 0
        type(VesselStructureUDT), allocatable       :: structures(:)
        type(PolygonSetUDT)                 :: polygonset, targetps 
        type(PolygonUDT), allocatable       :: targetpolygons(:)
        
        ! Shape function
        class(PolygonShapeFunctionUDT), allocatable     :: psf
        class(PolygonLevelsetFunction2DUDT), allocatable   :: plfvessel, &
            plftarget

    contains 

        ! Update vessel description using coordinates
        procedure :: UpdateVesselCoordinates

        ! Vessel coordinate getter
        procedure :: GetVesselCoordinates

    end type

    ! Environment
    !============
    type EnvironmentUDT

        ! Description
        !============
        ! Overarching type that stores all other structures etc which 
        ! may be needed for grid optimization, and which are not 
        ! related to the grid or the magnetic field. Currently, only
        ! the vessel structure is stored here. 

        ! Note: the routine to set up the vessel is currently a 
        ! standalone routine. Should we include it here as a 
        ! method of the vessel structure?

        type(VesselUDT)                 :: vessel

    contains

    end type
    
    !==================================================================!
    !                                                                  !
    !                             INTERFACES                           !
    !                                                                  !
    !==================================================================!

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    contains 


    !------------------------------------------------------------------!
    !                               Grid                               !
    !------------------------------------------------------------------!

    ! Allocation 
    !===========
    ! Main grid structure
    subroutine AllocateGrid(grid)

        ! Description
        !============
        ! Allocate the fields in the vertex, faces, and cell structures
        ! of the grid. 
        !
        ! Note: the bnd structure has to be allocated separately, since 
        ! the size of this structure may vary, depending on the rest of
        ! the grid. 

        ! The usual
        implicit none

        ! Declare variables
        type(GridUDT)       :: grid

        ! Return if already allocated
        if (allocated(grid%vert%x)) return

        ! Allocate vertex data
        call AllocateVertices(grid%vert)
        
        ! Allocate face data
        call AllocateFaces(grid%face)

        ! Allocate cell data
        call AllocateCells(grid%cell)

        ! Allocate other data structures
        call AllocateGridData(grid%data,grid)


    end subroutine

    ! Vertex substructure
    subroutine AllocateVertices(vert)

        ! Description
        !============
        ! Allocate the fields in the vertex structure. At least the 
        ! following scalar fields have to be present:
        !
        ! - ntot        : total number of vertices

        ! Note: derived quantities, such as the vertex neighbours, 
        ! cannot be allocated beforehand. These should be done 
        ! separately, e.g. in the ComputeGridInterconnections routine


        ! The usual
        implicit none

        ! Declare variables
        type(VertexUDT)       :: vert

        ! Allocate
        !=========
        ! Vertex data
        allocate(vert%x(vert%ntot))
        allocate(vert%y(vert%ntot))
        allocate(vert%BV(vert%ntot))
        allocate(vert%fieldlineID(vert%ntot))
        allocate(vert%bx(vert%ntot), vert%by(vert%ntot), &
            vert%ffbz(vert%ntot), vert%psi(vert%ntot))

        ! Face data
        allocate(vert%faceP(vert%ntot,2))
        !allocate(vert%face(vert%nface)) 

        ! Cell data
        allocate(vert%cellP(vert%ntot,2))
        !allocate(vert%cell(vert%ncell))

        ! Neighbour data
        allocate(vert%neigP(vert%ntot, 2))
        !allocate(vert%neig(vert%nneig))

    end subroutine

    ! Face substructure
    subroutine AllocateFaces(face)

        ! Description
        !============
        ! Allocate the fields in the faces structure. At least the 
        ! following scalar fields have to be present:
        !
        ! - ntot        : total number of faces

        ! Note: the cell structure cannot be allocated a priori - see 
        ! ComputeGridInterconnections


        ! The usual
        implicit none

        ! Declare variables
        type(FaceUDT)       :: face

        ! Allocate
        !=========
        ! Vertices
        allocate(face%vert(face%ntot,2))

        ! Cells
        allocate(face%cellP(face%ntot,2)) 
        ! allocate(face%cell(face%ncell))
        allocate(face%BF(face%ntot))

        ! Other
        allocate(face%aligned(face%ntot), face%label(face%ntot), &
            face%reg(face%ntot))

    end subroutine

    ! Cell substruture
    subroutine AllocateCells(cell)

        ! Description
        !============
        ! Allocate the fields in the cells structure. At least the 
        ! following scalar fields have to be present:
        !
        ! - ntot        : total number of cells
        ! - nvert       : length of cells%vert
        ! - nfacelist   : length of cells%face
        ! 
        ! The following fields are allocated: 

        ! Note: roughly the same conventions on data structures. fields
        ! 


        ! The usual
        implicit none

        ! Declare variables
        type(CellUDT)       :: cell

        ! Allocate
        !=========
        ! Guard cells
        allocate(cell%GC(cell%ntot))

        ! Vertex values
        allocate(cell%vertP(cell%ntot,2))
        allocate(cell%vert(cell%nvert)) 

        ! Face values
        allocate(cell%faceP(cell%ntot,2))
        allocate(cell%face(cell%nface))

        ! Others
        allocate(cell%psi(cell%ntot), cell%bp(cell%ntot), &
            cell%bt(cell%ntot), cell%x(cell%ntot), cell%y(cell%ntot), &
            cell%cflags(cell%ntot), cell%reg(cell%ntot), cell%ft(cell%ntot))

    end subroutine

    ! Boundary substructure
    subroutine AllocateBnd(bnd)

        ! Description
        !============
        ! Allocate the boundary fields. It is assumed that the following
        ! fields are present:
        !
        ! - nface:          number of boundary faces
        !
        ! Note: the amount of vertices is per definition the amount of
        ! faces plus one, also for closed boundaries. In the latter 
        ! case, the begin and end vertex index are the same. This 
        ! assumes that the boundaries should be polygons that either 
        ! close perfectly on themselves, or are simple polygons (i.e. 
        ! no polygon vertex should have more than three faces it belongs
        ! to)
        
        ! The usual
        implicit none

        ! Declare variables
        type(BndUDT)       :: bnd

        ! Allocate
        !=========
        bnd%nvert = bnd%nface+1
        allocate(bnd%face(bnd%nface))
        allocate(bnd%vert(bnd%nvert))
        

    end subroutine
 

    ! Additional grid data structures
    !================================
    ! Data substructure
    subroutine AllocateGridData(data,grid)

        ! Description
        !============
        ! Allocate the structures inside the grid data structures. It 
        ! is assumed that the main grid structures (cells/faces/vert) 
        ! are already properly initialized. 

        ! The usual
        implicit none

        ! Declare variables
        type(GridDataUDT)       :: data
        type(GridUDT)           :: grid

        ! Flux data
        call AllocateFluxData(data%fluxdata,grid)

        ! Region data
        call AllocateRegionData(data%regions,grid)

    end subroutine

    ! Flux data substructure
    subroutine AllocateFluxData(fluxdata,grid)

        ! Description
        !============
        ! Allocate the flux data structure. At least the following 
        ! fields have to be initialized (see definition of FluxDataUDT
        ! for explanation): 
        !
        ! - nFs
        ! - nFt
        ! - grid%cell%ntot
        ! - grid%face%ntot

        ! The usual
        implicit none

        ! Declare variables
        type(FluxDataUDT)       :: fluxdata
        type(GridUDT)           :: grid

        ! Allocate
        !=========
        ! Flux tube data
        allocate(fluxdata%fluxtubecellsP(fluxdata%nFt,2))
        allocate(fluxdata%fluxtubecells(grid%cell%ntot)) 
        allocate(fluxdata%fluxtubefacesP(fluxdata%nFt,2))
        allocate(fluxdata%fluxtubefaces(grid%face%ntot)) 

        ! Flux surface data
        allocate(fluxdata%fluxsurfacefacesP(fluxdata%nFs,2))
        allocate(fluxdata%fluxsurfacefaces(grid%face%ntot))
        allocate(fluxdata%fluxsurfaceID(grid%vert%ntot))
        allocate(fluxdata%fluxsurfacepsi(fluxdata%nFs))

    end subroutine

    ! Region data substrucure
    subroutine AllocateRegionData(regions,grid)

        ! Description
        !============
        ! Allocat the region data grid substructure. The following 
        ! fields should be present:
        !
        ! - grid%cell%ntot
        ! - grid%face%ntot
        ! - grid%data%fluxdata%nFt

        ! The usual
        implicit none

        ! Declare variables
        type(RegionDataUDT)         :: regions
        type(GridUDT)               :: grid

        ! Allocate
        !=========
        allocate(regions%fluxtuberegID(grid%data%fluxdata%nFt))

    end subroutine

    ! Deallocation
    !=============
    ! Vertex substructure
    subroutine DeallocateVertices(vert)

        ! Description
        !============
        ! Deallocate the fields in the vert structure

        ! The usual
        implicit none

        ! Declare variables
        type(VertexUDT)       :: vert

        ! Deallocate
        !=========
        ! Vertex data
        deallocate(vert%x)
        deallocate(vert%y)
        deallocate(vert%BV)
        deallocate(vert%fieldlineID)

        ! Face data
        deallocate(vert%faceP)
        deallocate(vert%face) 

        ! Cell data
        deallocate(vert%cellP)
        deallocate(vert%cell)

        ! Neighbour data
        deallocate(vert%neigP)
        deallocate(vert%neig)

    end subroutine

    ! Basic operations
    !=================
    ! Get cells of a vertex
    function GetVertCell(vert, i) result(res)
        integer(I8)                 :: i 
        type(VertexUDT)             :: vert 
        integer(I8), allocatable    :: res(:)
        res = vert%cell(vert%cellP(i, 1):(vert%cellP(i, 1) + vert%cellP(i, 2) - 1))
    end function

    ! Get faces of a vertex
    function GetVertFace(vert, i) result(res)
        integer(I8)                 :: i 
        type(VertexUDT)             :: vert 
        integer(I8), allocatable    :: res(:)
        res = vert%face(vert%faceP(i, 1):(vert%faceP(i, 1) + vert%faceP(i, 2) - 1))
    end function

    ! Get neig of a vertex
    function GetVertNeig(vert, i) result(res)
        integer(I8)                 :: i 
        type(VertexUDT)             :: vert 
        integer(I8), allocatable    :: res(:)
        res = vert%neig(vert%neigP(i, 1):(vert%neigP(i, 1) + vert%neigP(i, 2) - 1))
    end function

    ! Get cells of a face
    function GetFaceCell(face, i) result(res)
        integer(I8)                 :: i 
        type(FaceUDT)               :: face 
        integer(I8), allocatable    :: res(:)
        res = face%cell(face%cellP(i, 1):(face%cellP(i, 1) + face%cellP(i, 2) - 1))
    end function

    ! Get faces of a cell
    function GetCellFace(cell, i) result(res)
        integer(I8)                 :: i 
        type(CellUDT)               :: cell 
        integer(I8), allocatable    :: res(:)
        res = cell%face(cell%faceP(i, 1):(cell%faceP(i, 1) + cell%faceP(i, 2) - 1))
    end function

    ! Get vertices of a cell
    function GetCellVert(cell, i) result(res)
        integer(I8)                 :: i 
        type(CellUDT)               :: cell 
        integer(I8), allocatable    :: res(:)
        res = cell%vert(cell%vertP(i, 1):(cell%vertP(i, 1) + cell%vertP(i, 2) - 1))
    end function

    ! Get faces of a flux surface
    function GetFSFace(fd, i) result(res)
        integer(I8)                 :: i 
        type(FluxDataUDT)           :: fd 
        integer(I8), allocatable    :: res(:)
        res = fd%fluxsurfacefaces(fd%fluxsurfacefacesP(i, 1):&
            (fd%fluxsurfacefacesP(i, 1) + fd%fluxsurfacefacesP(i, 2) - 1))
    end function

    ! Get cells of a flux tube
    function GetFTCell(fd, i) result(res)
        integer(I8)                 :: i 
        type(FluxDataUDT)           :: fd 
        integer(I8), allocatable    :: res(:)
        res = fd%fluxtubecells(fd%fluxtubecellsP(i, 1):&
            (fd%fluxtubecellsP(i, 1) + fd%fluxtubecellsP(i, 2) - 1))
    end function

    ! Get faces of a flux tube
    function GetFTFace(fd, i) result(res)
        integer(I8)                 :: i 
        type(FluxDataUDT)           :: fd 
        integer(I8), allocatable    :: res(:)
        res = fd%fluxtubefaces(fd%fluxtubefacesP(i, 1):&
            (fd%fluxtubefacesP(i, 1) + fd%fluxtubefacesP(i, 2) - 1))
    end function

    ! Interconnections
    subroutine ComputeGridInterconnections(grid)

        ! Description
        !============
        ! Compute additional grid topology data starting from basic grid
        ! information. The following fields are required (see gdmod_types
        ! for an explanation of the different fields)
        !
        ! Vertices: ntot, BV
        ! Faces: ntot, vert
        ! Cells: ntot, vertP, vert, nvert 
        !   
        ! All other fields are recomputed based on this basic 
        ! interconnection data.             
    
        ! Notes
        !======
        ! Note 1: guard cells are identified by checking which cells only 
        ! have two vertices. It
    
        ! Note 2: it is assumed that all fields of which the size is 
        ! predetermined by the number of vertices/faces/cells are already
        ! allocated. Other fields are allocated or overwritten. 
    
        ! Note 3: the original cell vertices are replaced by a sorted list.
        ! This list is generated based on the faces of the cell and their
        ! corresponding vertices. It is assumed that each cell has the same
        ! amount of faces as it has vertices.
    
        ! Note 4: the vertex neighbours and vertex cells are ordened in 
        ! clockwise or counter-clockwise manner. This sorting assumes that
        ! each vertex has more than two cells - otherwise, it doesn't really
        ! matter how it is sorted. If it has more than two cells, the cells
        ! have to be properly connected, in the sense that each cell should
        ! only have one neighbour. 
    
        ! Initialize
        !===========
        ! The usual
        implicit none 
    
        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT)               :: grid
    
        ! Loop variables
        integer(I8)                 :: i, j, k, m, q
    
        ! Auxiliary variables 
        type(VertexUDT)             :: v 
        type(FaceUDT)               :: f
        type(CellUDT)               :: c
    
        integer(I8)                 :: nc, nv, nf, sp, ep, tcs, thiscell, &
            nvc, ntv, tf, ind, tfv(1:2), tc, tcn, ncf, nextcell, tcf, vc, &
            ntcf2, nb1, nb2, nfc, fcc
        integer(I8), allocatable    :: fcount(:), vcount(:), tv(:), &
            tempfcell(:, :), localID(:), fc(:), fn(:), allvertcells(:), &
            allfv(:), tcf2(:), nb(:), vcc(:), indv(:)
    
        logical, allocatable        :: cellfound(:), allfvind(:)
        logical                     :: startcellnotfound, allnotfound, &
            accountforGC
    
        ! Unpack & initialize
        !==================
        ! Data structures 
        v = grid%vert
        f = grid%face
        c = grid%cell 
    
        ! Checks
        if (size(f%vert,2) /= 2) then
            ! This should actually never happen, but just to be sure
            call gdErrorHandler('ComputeGridInterconnections: too many vertices per face')
        end if
        if (any(any(f%vert == 0,2),1)) then
            ! This could be the case when the face vertices are not correctly determined
            call gdErrorHandler('ComputeGridInterconnections: some vertex indices are zero in faces, not supported')
        end if
    
        ! Quantities
        nv = v%ntot
        nf = f%ntot
        nc = c%ntot
    
        ! Counters
        allocate(fcount(nf))
        allocate(vcount(nv))
    
        ! Face neighbours, vertex cells, cell faces
        !==========================================
        ! Check allocation
        if (allocated(c%face)) then
            ! Print out warning
            print *, 'ComputeGridInterconnections: recomputing cell%face'
    
            ! Deallocate
            deallocate(c%face)
        end if
        if (allocated(v%cell)) then 
            ! print out warning
            print *, 'ComputeGridInterconnections: recomputing vert%cell'
    
            ! Deallocate
            deallocate(v%cell)
        end if
        if (allocated(f%cell)) then 
            ! Print out warning
            print *, 'ComputeGridInterconnections: recomputing face%cell'
    
            ! Deallocate
            deallocate(f%cell)
        end if
        if (allocated(v%cell)) then 
            ! Print out warning
            print *, 'ComputeGridInterconnections: recomputing vert%cell'
    
            ! Deallocate
            deallocate(v%cell)
        end if
    
        ! Check which cells are guard cells
        c%GC(:) = .false. ! initialize to false
        do i = 1, nc 
            if (c%vertP(i,2) <= 2) then
                ! Guard cell found, set true
                c%GC(i) = .true.
            end if
        end do
    
        ! Check if we have to account for guard cells
        accountforGC = .false.
        if (any(c%GC)) then
            accountforGC = .true.
        end if 
    
        ! Make the pointer list for cell faces
        c%faceP = 0
        c%faceP(1,1) = 1
        if (c%GC(1)) then 
            ! Only one face
            c%faceP(1,2) = 1
        else
            ! Same amount of faces as vertices 
            c%faceP(1,2) = c%vertP(1,2) 
        end if
    
        do i = 2, nc
            ! Pointer
            c%faceP(i,1) = c%faceP(i-1,1) + c%faceP(i-1,2)
    
            ! Number of faces
            if (c%GC(i)) then 
                c%faceP(i,2) = 1
            else
                c%faceP(i,2) = c%vertP(i,2)
            end if
    
        end do
    
        ! Compute the total faces in cells%face
        c%nface = sum(c%faceP(:,2)) 
    
        ! Allocate
        allocate(c%face(c%nface))
        c%face = 0
    
        ! Make the pointer list for vertex cell neighbours
        v%cellP = 0
        do i = 1, nc
            ! Get the number of vertices of the current cell
            ntv = c%vertP(i,2)
    
            ! Allocate
            allocate(tv(ntv))
    
            ! Get the vertices of the current cell
            tv = c%vert(c%vertP(i,1):(c%vertP(i,1)+ntv-1))
    
            ! Update the counter
            v%cellP(tv,2) = v%cellP(tv,2)+1
    
            ! Deallocate
            deallocate(tv)
    
        end do
    
        ! Compute the total number of cells in vert%cell
        v%ncell = sum(v%cellP(:,2))
    
        ! Allocate
        allocate(v%cell(v%ncell))
        v%cell = 0
    
        ! Construct the pointer
        v%cellP(1,1) = 1
        do i = 2, nv 
            v%cellP(i,1) = v%cellP(i-1,1) + v%cellP(i-1,2)
        end do
    
        ! Set vertex cell and cell face and face neighbours
        fcount = 1
        vcount = 1
        allocate(tempfcell(f%ntot, 2))
        tempfcell = 0
        do i = 1, nc ! loop over all cells
            ! Get vertices of this cell
            tv = GetCellVert(c, i)
            ntv = c%vertP(i,2)
    
            ! Loop over all vertices
            do j = 1, ntv
    
                ! Get the faces of the current cell
                ! Get the next face 
                if (j == ntv) then
                    k = 1
                else
                    k = j+1
                end if 
    
                call MapVertexPairToFace(tv(j), tv(k), f%vert, nf, tf)
    
                ! Check
                if (tf < 0) then
                    ! This shouldn't happen
                    print *, 'vertex pair: ', tv(j), tv(k)
                    call gdErrorHandler(&
                        'ComputeGridInterconnections: could not find '&
                        // 'matching face for this vertex pair')
                else
    
                    ! Add the current cell index as neighbour to this face
                    if (.not. ((j > 1) .and. c%GC(i))) then ! hedge for guard cells
                        ! Sanity check
                        if (fcount(tf) > 2) then
                            print *, 'face ID: ', tf
                            print *, 'neighbours: ', tempfcell(tf,:), i
                            call gdErrorHandler(& 
                            'ComputeGridInterconnections: too many ' &
                                // 'neighbours for this face')
                        end if
    
                    
                        tempfcell(tf,fcount(tf)) = i
    
                        ! Update fcount
                        fcount(tf) = fcount(tf)+1
                    end if
    
                    ! Add the current cell to the j'th vertex
                    ind = v%cellP(tv(j),1) + vcount(tv(j)) - 1
                    if (ind > v%ncell) then
                        call gdErrorHandler('unknown error')
                    end if
                    v%cell(ind) = i
    
                    ! Update vcount
                    vcount(tv(j)) = vcount(tv(j))+1
    
                    ! Add cell face
                    if (.not. ((j > 1) .and. c%GC(i))) then ! hedge for guard cells
                        ind = c%faceP(i,1) + j - 1
                        if (ind > c%nface) then
                            call gdErrorHandler('unknown error')
                        end if
                        c%face(ind) = tf
                    end if
                end if
    
            end do
    
            ! Deallocate
            deallocate(tv)
    
        end do
    
        ! Construct cell arrays for faces and vertices
        fcount = fcount-1
        f%cellP(:, 2) = fcount 
        f%cellP(1, 1) = 1
        f%ncell = sum(f%cellP(:, 2))
        allocate(f%cell(f%ncell))
        fcc = 0
        nfc = fcount(1)
        f%cell(fcc+1:fcc+nfc) = tempfcell(1, 1:nfc)
        fcc = fcc + nfc
        do i = 2, f%ntot
            f%cellP(i, 1) = f%cellP(i-1, 2) + f%cellP(i-1, 1)
            nfc = fcount(i)
            f%cell(fcc+1:fcc+nfc) = tempfcell(i, 1:nfc)
            fcc = fcc + nfc
        end do
    
        ! Housekeeping
        deallocate(vcount, fcount)
        
        ! Logicals
        !=========
        ! Based on the interconnections computed above, the logical indices
        ! (boundary vertices => BV, boundary faces => BF) etc can be 
        ! computed. Boundary faces are either: 1) faces with a single 
        ! neighbour (if no guard cells are present), 2) faces with a guard
        ! cell neighbour (if guard cells are present)
    
        ! Compute BF & BV
        f%BF(:) = .false. ! initialize
        v%BV(:) = .false. ! initialize
        do i = 1, f%ntot
    
            ! Get face neighbours
            nb  = GetFaceCell(f, i)
            nb1 = 0 ! Initialize to be sure
            nb2 = 0
            if (size(nb, 1) == 1) then 
                nb1 = nb(1)
                nb2 = 0
            elseif (size(nb, 1) == 2) then   
                nb1 = nb(1)
                nb2 = nb(2)
            else 
                ! This shouldn't happen
                call gdErrorHandler('More than two cells detected as face neighbours, check grid interconnection')
            end if 
            !deallocate(nb)
    
            ! Sanity checks
            if ( (nb1 == 0) .or. (accountforGC .and. (nb2 == 0)) ) then 
                ! Call error, this should not be the case
                print *, 'face ID: ', i
                call gdErrorHandler('ComputeGridInterconnections:' &
                // 'could not find neighbours for this face')
            end if
    
            ! Set BF & BV logical
            if (accountforGC) then
                if (c%GC(nb2) .or. c%GC(nb1) ) then
                    ! This is a boundary face
                    f%BF(i) = .true.
                    v%BV(f%vert(i,:)) = .true.
                end if
            else
                if (nb2 == 0) then
                    ! This is a boundary face
                    f%BF(i) = .true.
                    v%BV(f%vert(i,:)) = .true.
                end if
            end if
        end do
    
    
        ! Additional vertex interconnections
        !===================================
        ! Vertex neighbours
        ! Check allocation
        if (.not. allocated(v%neigP)) then
            allocate(v%neigP(v%ntot,2))
        end if
    
        ! Compute how much neighbours each vertex has by looping over faces
        v%neigP = 0
        do i = 1, f%ntot
            ! Get face vertices
            tfv = f%vert(i,:)
    
            ! Update counters
            v%neigP(tfv(1),2) = v%neigP(tfv(1),2) + 1
            v%neigP(tfv(2),2) = v%neigP(tfv(2),2) + 1
        end do
    
        ! Compute the first column of neigP
        v%neigP(1,1) = 1
        do i = 2, v%ntot
            v%neigP(i,1) = v%neigP(i-1,1) + v%neigP(i-1,2)
        end do
    
        ! Check the neiglist allocation
        f%ncell = sum(v%neigP(:,2))
        if (allocated(v%neig)) then 
            ! This shouldn't be the case, but we hedge for it. Reallocate
            deallocate(v%neig)
            allocate(v%neig(f%ncell))
        else
            allocate(v%neig(f%ncell))
        end if
        v%neig = 0
    
        ! Compute the neiglist
        allocate(vcount(v%ntot))
        vcount = 0
        do i = 1, f%ntot
            ! Get the face vertices
            tfv = f%vert(i,:)
    
            ! Add the neighbours
            sp = v%neigP(tfv(1),1) + vcount(tfv(1))
            v%neig(sp) = tfv(2)
            sp = v%neigP(tfv(2),1) + vcount(tfv(2))
            v%neig(sp) = tfv(1)
    
            ! Update counter
            vcount(tfv) = vcount(tfv) + 1
    
        end do
        deallocate(vcount)
    
        ! Compute vertex faces
        v%faceP = 0
        do i = 1, f%ntot
            v%faceP(f%vert(i, :), 2) = v%faceP(f%vert(i, :), 2) + 1;
        end do
        v%nface = sum(v%faceP(:, 2))
        v%faceP(1, 1) = 1
        do i = 2, v%ntot 
            v%faceP(i, 1) = v%faceP(i-1, 1) + v%faceP(i-1, 2)
        end do 
        if (allocated(v%face)) then 
            deallocate(v%face)
        end if 
        allocate(v%face(v%nface), vcc(v%ntot)) 
        vcc(:) = 0
        do i = 1, f%ntot
            tfv = f%vert(i, :)
            indv = vcc(tfv) + v%faceP(tfv, 1)
            v%face(indv) = i
            vcc(tfv) = vcc(tfv) + 1
        end do 
    
    
        ! Sort vertex neighbours and vertex cells
        !========================================
        
        ! Loop over all vertices
        do i = 1, v%ntot ! v%ntot
    
            ! Check how many distinct cell sequences there are by checking 
            ! the amount of faces of this vertex and the amount of cells. 
            ! If the difference is zero -> internal vertex (check if this is
            ! true, otherwise throw error), starting cell doesn't matter. 
            ! If one or higher -> only possible if boundary vertex 
            ! (check this). Difference indicates the amount of cell 
            ! sequences.
    
            ! Note: we have to account for guard cells here. Normally, if 
            ! guard cells are present, each boundary vertex should have a
            ! multiple of two amount of guard cells (more than 2 only holds
            ! for boundary vertices that connect two domains in a pointwise
            ! way (or more domains)).
    
            ! Number of cell sequences
            tcs = 0
            tcs = v%neigP(i, 2) - v%cellP(i, 2)
    
            ! Account for guard cells, if any
            if (accountforGC) then 
                tcs = tcs + count(c%GC(v%cell(v%cellP(i, 1):v%cellP(i, 2)+v%cellP(i, 1)-1)))
            end if
    
            ! Check
            if ( (tcs == 0) .and. (v%BV(i)) ) then 
                ! Throw error - this vertex should be an internal vertex,
                ! not a boundary one. Check the grid interconnection. 
                print *, 'vertex ID: ', i
                call gdErrorHandler(&
                    'ComputeGridInterconnections: this vertex is ' &
                    // 'falsely classified as boundary vertex')
            end if
            if ( (tcs > 0) .and. (.not. v%BV(i)) ) then 
                ! Throw error - this vertex should be a boundary vertex.
                ! Check the grid interconnection. 
                print *, 'vertex ID: ', i
                call gdErrorHandler(&
                    'ComputeGridInterconnections: this vertex is ' &
                    // 'falsely classified as internal vertex')
            end if
    
            ! Issue warning when tcs is higher than 1 - only occurs if two
            ! boundaries touch each other in one point, which should hardly
            ! ever be the case. 
            if (tcs > 1) then 
                print *, 'WARNING: boundary vertex ', i, &
                    ' found with multiple cell sequences - check interconnection'
            end if 
    
            ! Set to one if equal to zero (for internal vertices)
            if (tcs == 0) then
                tcs = 1
            end if
    
            ! Loop over all cell sequences
            nvc = v%cellP(i,2)
            allocate(cellfound(nvc))
            allocate(localID(nvc))
            allocate(allvertcells(nvc))
            cellfound(:) = .false.
            sp = v%cellP(i,1)
            ep = v%cellP(i,1) + nvc-1
            allvertcells = v%cell(sp:ep)
    
            tcf = 0
            do j = 1, tcs
    
                ! Get the first cell and its faces
                startcellnotfound = .true.
                thiscell = 0
                if (.not. v%BV(i)) then 
                    ! Internal cell, simply take the first cell, since only 
                    ! one sequence
                    sp = v%cellP(i,1)
                    thiscell = v%cell(sp)
    
    
                    ! Update 
                    cellfound(1) = .true.
                    startcellnotfound = .false.
    
                else
                    ! Boundary cell, search for cell with only one (or none)
                    ! faces in common with the other cells
    
                    ! Find cell with one or zero common faces with the other
                    ! cells
                    
                    m = 1
                    do while (startcellnotfound .and. (m <= nvc))
    
                        ! Get the current cell
                        if (cellfound(m)) then
                            ! Update m and skip the rest
                            m = m + 1
                            cycle
                        end if
    
                        ! Get the faces of the current cell
                        ! tc = allremcells(m) ! current cell
                        tc = allvertcells(m)
                        fc = GetCellFace(c, tc)
                        !sp = c%faceP(tc,1)
                        !ep = c%faceP(tc,1) + c%faceP(tc,2) - 1 
                        !allocate(fc(c%faceP(tc,2)))
                        !fc = c%face(sp:ep)
    
                        ! Check how many common faces there are with faces
                        ! of the other remaining cells
                        ncf = 0 ! number of common faces
                        do k = 1, nvc
                            ! Skip if k == m or if cell is already found
                            if ((k == m) .or. cellfound(k)) then 
                                cycle 
                            end if
    
                            ! Get the faces of the next cell
                            tcn = allvertcells(k)
                            sp = c%faceP(tcn,1)
                            ep = c%faceP(tcn,1) + c%faceP(tcn,2) - 1 
                            allocate(fn(c%faceP(tcn,2)))
                            fn = c%face(sp:ep)
    
                            ! Count how many faces they have in common
                            do q = 1, c%faceP(tc,2)
                                ncf = ncf + count(fc(q) == fn)
                            end do
    
                            ! Deallocate
                            deallocate(fn)
    
                        end do
    
                        ! Deallocate
                        deallocate(fc)
    
                        ! Check how many common faces there are
                        !print *, ncf
                        if (ncf < 2) then 
                            ! Cell found
                            startcellnotfound = .false.
    
                            ! Set current cell
                            sp = v%cellP(i,1) + m - 1
                            thiscell = v%cell(sp)
    
                            ! Update logicals
                            ! cellfound(thislocalID(m)) = .true.
                            cellfound(m) = .true.
    
                            ! Set the remaining cells
                            !allocate(tcrem(nremcells-1))
                            !tcrem = allremcells((/(q, q = i, m-1),(q, q = m+1,nremcells)/))
    
    
                        else
                            ! Update m
                            m = m + 1
    
                        end if
    
                    end do
    
                    ! Check
                    if (startcellnotfound) then
                        ! Should have found starting cell, throw error. 
                        ! Something is wrong then with the interconnection,
                        ! and it could be that this boundary vertex should 
                        ! in fact not be a boundary vertex. 
                        print *, 'vertex ID: ', i 
                        call gdErrorHandler('ComputeGridInterconnections: could not find starting cell for this boundary vertex')
    
                    end if 
    
    
                end if
    
                ! Find the sequence
                !==================
                ! Now that the first cell of the sequence has been found, we
                ! need to extract the sequence itself. 
    
                ! Loop over the remaining cells
                k = 1
                vc = 0 ! vertex counter
                allnotfound = .true. 
                
                do while (k <= nvc)
    
                    ! Set the faces of the current cell
                    allocate(fc(c%faceP(thiscell,2)))
                    sp = c%faceP(thiscell, 1)
                    ep = sp + c%faceP(thiscell,2)-1
                    fc = c%face(sp:ep)
    
                    ! Take the next cell
                    if (cellfound(k)) then
                        ! Update k, cycle
                        k = k + 1
                        deallocate(fc)
                        cycle
                    else
                        nextcell = allvertcells(k)
                    end if
                    
                    ! Get cell faces
                    allocate(fn(c%faceP(nextcell,2)))
                    sp = c%faceP(nextcell, 1)
                    ep = sp + c%faceP(nextcell,2)-1
                    fn = c%face(sp:ep)
    
                    ! Check for common faces - should only be one
                    ncf = 0
                    do q = 1, c%faceP(thiscell, 2)
                        ncf = ncf + count(fc(q) == fn)
                        if ((count(fc(q) == fn)) == 1) then
                            ! Set the current face
                            tcf = fc(q)
                        end if
                    end do
    
                    ! Sanity check
                    if (ncf > 1) then
                        ! Throw error: the next cell is connected with too
                        ! many faces to the previous one, which shouldn't be
                        ! possible in any decent grid. Though in some 
                        ! exceptional cases, this may be allowed, we don't 
                        ! support that here. 
                        call gdErrorHandler('ComputeGridInterconnections: '&
                            // 'too many common faces detected when ' &
                            // 'sorting vertex neighbours, check grid ' &
                            // 'interconnection')
                    end if
    
                    ! Check if a common face was found
                    if (ncf == 1) then 
                        ! Add the first vertex if this is the first found 
                        ! cell
                        if (vc == 0) then
                            ! Get other face with this vertex
    
                            ! Get the faces that have the current vertex
                            if (c%GC(thiscell)) then
                                ! Hedge for guard cells: only one face
                                allocate(allfv(1))
                                allocate(allfvind(1))
                                allocate(tcf2(1))
                                tcf2(1) = fc(1)
                            else
                                allocate(allfvind(c%faceP(thiscell,2)))
                                allfvind = (f%vert(fc, 1) == i) .or. &
                                    (f%vert(fc, 2) == i)
                                allocate(allfv(count(allfvind)))
                                allfv = pack(fc, allfvind)
                                allocate(tcf2(count(allfvind)))
                                tcf2 = pack(allfv, allfv .ne. tcf)
                            end if
                            
                            ! Check 
                            ntcf2 = size(tcf2)
                            if (ntcf2 .ne. 1) then
                                call gdErrorHandler( &
                                'ComputeGridInterconnections: could not ' &
                                // 'find a single common face, check grid' &
                                //' interconnections')
                            end if
    
                            ! Add cell and vertex neighbour - neighbour only
                            ! if the current cell is not a guard cell
                            
                            sp = v%cellP(i,1)
                            v%cell(sp) = thiscell 
                            if (.not. c%GC(thiscell)) then
                                sp = v%neigP(i,1)
                                tfv = f%vert(tcf2(1),:)
                                if (tfv(1) == i) then
                                    v%neig(sp) = tfv(2)
                                else
                                    v%neig(sp) = tfv(1)
                                end if
    
                                ! Update vc
                                vc = vc+1
                            end if
    
                            ! Housekeeping
                            deallocate(allfvind)
                            deallocate(allfv)
                            deallocate(tcf2)
    
                        end if 
    
                        ! Add the cell and vertex neighbour
                        sp = v%cellP(i,1) + vc
                        v%cell(sp) = nextcell
                        sp = v%neigP(i,1) + vc
                        tfv = f%vert(tcf,:)
                        if (tfv(1) == i) then 
                            v%neig(sp) = tfv(2)
                        else
                            v%neig(sp) = tfv(1)
                        end if
    
                        ! Update counter
                        vc = vc + 1
    
                        ! Update logicals
                        cellfound(k) = .true.
    
                        ! Reset iterator
                        k = 1
    
                        ! Set current cell
                        thiscell = nextcell
    
                    else
                        ! Update the counter
                        k = k + 1
                    end if
    
                    ! Housekeeping
                    deallocate(fc, fn)
    
                end do
    
                ! Add last vertex
                !================
                ! Only in case of boundary vertices
                if ( (v%BV(i)) .and. (v%cellP(i,2) > 1) ) then
    
                    ! First, update vc - is now one-off
                    ! vc = vc - 1
    
                    ! Sanity check
                    if (vc > v%neigP(i,2)) then
                        ! Throw error. This can happen if the supposed 
                        ! boundary vertex is no actual boundary vertex.
                        call gdErrorHandler('ComputeGridInterconnections: ' &
                            // 'supposed boundary vertex may be internal ' &
                            // 'vertex. Check grid interconnection')
                    end if
    
                    ! Set the faces of the current cell
                    allocate(fc(c%faceP(thiscell,2)))
                    sp = c%faceP(thiscell, 1)
                    ep = sp + c%faceP(thiscell,2)-1
                    fc = c%face(sp:ep)
    
                    ! Add the vertex neighbour
                    if (c%GC(thiscell)) then
                        ! Hedge for guard cells: only one face
                        allocate(tcf2(1))
                        allocate(allfv(1))
                        allocate(allfvind(1))
                        tcf2(1) = fc(1)
                    else
                        allocate(allfvind(c%faceP(thiscell,2)))
                        allfvind = (f%vert(fc, 1) == i) .or. &
                            (f%vert(fc, 2) == i)
                        allocate(allfv(count(allfvind)))
                        allfv = pack(fc, allfvind)
                        tcf2 = pack(allfv, allfv .ne. tcf)
    
                    end if
                    
                    ! Check 
                    ntcf2 = size(tcf2)
                    if (ntcf2 .ne. 1) then
                        call gdErrorHandler( &
                        'ComputeGridInterconnections: could not ' &
                        // 'find a single common face, check grid' &
                        //' interconnections')
                    end if
    
                    ! Add cell and vertex neighbour
                    sp = v%neigP(i,1) + vc
                    tfv = f%vert(tcf2(1),:)
                    if (tfv(1) == i) then
                        v%neig(sp) = tfv(2)
                    else
                        v%neig(sp) = tfv(1)
                    end if
    
                    ! Housekeeping
                    deallocate(allfvind)
                    deallocate(allfv)
                    deallocate(fc, tcf2)
                    
                end if
            end do
    
            ! Sanity check
            if (any(.not. cellfound)) then
                ! Not all cells were ordened, throw error
                print *, 'vertex ID: ', i
                call gdErrorHandler('ComputeGridInterconnections: could not sort all cells, check grid interconnection')
            end if 
    
            ! Housekeeping
            deallocate(cellfound)
            deallocate(localID)
            deallocate(allvertcells)
    
        end do
    
        !do i = 1, v%ntot 
        !    print *, i, v%BV(i), v%neig(v%neigP(i,1):v%neigP(i,1)+v%neigP(i,2)-1)
        !end do
        !    do i = 1, v%ntot 
        !    print *, i, v%BV(i), v%cell(v%cellP(i,1):v%cellP(i,1)+v%cellP(i,2)-1)
        !end do
        
    
        ! Add to grid
        !============
        grid%cell   = c 
        grid%face   = f 
        grid%vert    = v
    
    
    end subroutine


    !------------------------------------------------------------------!
    !                             Numerics                             !
    !------------------------------------------------------------------!

    ! Grid data updating
    subroutine UpdateGridData(grid, magneticField, environment)

        ! Description
        !============
        ! Update grid data such as metrics, magnetic field at vertices,
        ! ... that depends on the magnetic field and possibly on 
        ! the environment. Assumes all other fields are properly 
        ! initialized and allocated. 

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(inout)        :: grid 
        type(MagneticFieldUDT), intent(in)  :: magneticField 
        type(EnvironmentUDT), intent(in)    :: environment

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: tv, tf, tfv 
        real(R8), allocatable, dimension(:)     :: tpsi
        
        ! Loop
        integer(I8)                             :: i

        ! Initialize
        !===========
        ! Associate
        associate(mf        => magneticField%interp,    &
                nc          => grid%cell%ntot,          &
                nfs         => grid%data%fluxdata%nfs)

        ! Compute
        !========
        ! Magnetic field at vertices
        call mf%Evaluate(grid%vert%x, grid%vert%y, 0, 0, grid%vert%psi)
        call mf%Evaluate(grid%vert%x, grid%vert%y, 1, 0, grid%vert%bx)
        call mf%Evaluate(grid%vert%x, grid%vert%y, 0, 1, grid%vert%by)

        ! Cell centers
        do i = 1, nc 
            tv = GetCellVert(grid%cell, i)
            grid%cell%x(i) = sum(grid%vert%x(tv))/size(tv)
            grid%cell%y(i) = sum(grid%vert%y(tv))/size(tv)
        end do 

        ! Magnetic field at cell centers
        call mf%Evaluate(grid%cell%x, grid%cell%y, 0, 0, grid%cell%psi)
        do i = 1, nc
            tv = GetCellVert(grid%cell, i)
            grid%cell%bp(i) = -sum(sqrt(grid%vert%bx(tv)**2 + grid%vert%by(tv)**2))/(size(tv)*2*pi_R8*grid%cell%x(i))
            grid%cell%bt(i) = sum(grid%vert%ffbz(tv))/(size(tv)*2*pi_R8*grid%cell%x(i))
        end do

        ! Flux surface psi values
        do i = 1, nfs 

            ! Get flux surface faces
            tf = GetFSFace(grid%data%fluxdata, i)

            ! Get vertices
            tfv = reshape(grid%face%vert(tf, :), [size(tf)*2])

            ! Get psi values
            tpsi = grid%vert%psi(tfv)

            ! Average and add
            grid%data%fluxdata%fluxsurfacepsi = sum(tpsi)/size(tpsi)

        end do

        ! Housekeeping
        !=============
        end associate


    end subroutine

    !------------------------------------------------------------------!
    !                          Magnetic field                          !
    !------------------------------------------------------------------!

    ! Magnetic field
    !===============
    subroutine AllocateMagneticField(magneticField)

        ! Description
        !============
        ! Allocate the magnetic field properties. At least the fields nR
        ! and nZ should be present. 

        ! The usual
        implicit none

        ! Declare variables
        type(MagneticFieldUDT)          :: magneticField

        ! Allocate
        !=========
        allocate(magneticField%R(magneticField%nR))
        allocate(magneticField%Z(magneticField%nZ))
        allocate(magneticField%Psi(magneticField%nR, magneticField%nZ))

    end subroutine

    !------------------------------------------------------------------!
    !                            Environment                           !
    !------------------------------------------------------------------!
    ! Vessel structures
    !==================
    ! Allocation
    subroutine AllocateVesselStructure(vesselstructure)

        ! Description
        !============
        ! Routine to allocate the vessel structure quantities assuming
        ! that np is given.

        ! The usual
        implicit none

        ! Declare variables
        type(VesselStructureUDT)          :: vesselstructure

        ! Allocate
        !=========
        allocate(vesselstructure%x(vesselstructure%np))
        allocate(vesselstructure%y(vesselstructure%np))

    end subroutine

    ! Vessel
    !=======
    ! Allocation
    subroutine AllocateVessel(vessel)

        ! Description
        !============
        ! Allocate the vessel quantities assuming that np and ntp are 
        ! given. 

        ! The usual
        implicit none

        ! Declare variables
        type(VesselUDT)          :: vessel

        ! Allocate
        !=========
        !allocate(vessel%x(vessel%nv))
        !allocate(vessel%y(vessel%nv))
        !allocate(vessel%TPind(vessel%nv))
        !allocate(vessel%allTPind(vessel%ntv))

    end subroutine

    ! Deallocation
    subroutine DeallocateVessel(vessel)

        ! Description
        !============
        ! Deallocate the vessel quantities.

        ! The usual
        implicit none

        ! Declare variables
        type(VesselUDT)          :: vessel

        ! Allocate
        !=========
        deallocate(vessel%TPind)
        deallocate(vessel%allTPind)

    end subroutine

    ! Getting coordinates (assumed already initialized)
    subroutine GetVesselCoordinates(vessel, xv, yv)

        ! Description
        !============
        ! This routine returns the vessel coordinates in the same order
        ! as they are updated. Useful for e.g. shape optimization 
        ! purposes. 

        ! Declare variables
        !==================
        ! Arguments
        class(VesselUDT)                    :: vessel 
        real(R8), allocatable, intent(out)  :: xv(:), yv(:)

        ! Auxiliary
        integer(I8)                     :: nv

        ! Loop 
        integer(I8)                     :: i, cc

        ! Compute total number of vertices
        nv = 0
        do i = 1, vessel%polygonset%np
            nv = nv + vessel%polygonset%polygons(i)%nv
        end do

        ! Extract
        cc = 0
        if (allocated(xv)) then 
            deallocate(xv)
        end if 
        if (allocated(yv)) then 
            deallocate(yv)
        end if 
        allocate(xv(nv), yv(nv))
        do i = 1, vessel%polygonset%np 
            xv(cc+1:cc+vessel%polygonset%polygons(i)%nv) = vessel%polygonset%polygons(i)%x
            yv(cc+1:cc+vessel%polygonset%polygons(i)%nv) = vessel%polygonset%polygons(i)%y
        end do


    end subroutine

    ! Updating (assumed already initialized)
    subroutine UpdateVesselCoordinates(vessel, xv, yv)

        ! Description
        !============
        ! This routine updates the vessel coordinates and all other 
        ! dependent structures according to the new vessel coordinates
        ! given by the xv, yv pairs. It is assumed that xv, yv has an 
        ! equal number of entries as there are vertices in the vessel
        ! polygon structure. 

        ! Notes
        !======
        ! Note 1: the targets are NOT updated yet!
        ! Note 2: we must construct the polygon sets using the polygon
        ! vertices, not simply passing the new set of coordinates! This
        ! will lead to an open polygon by default...

        ! Declare variables
        !==================
        ! Arguments
        class(VesselUDT)                :: vessel 
        real(R8), intent(in)            :: xv(:), yv(:)

        ! Auxiliary
        integer(I8)                     :: flag, npv, npvtot
        real(R8), allocatable           :: xvp(:), yvp(:)
        character(:), allocatable       :: vesselpath 

        ! Loop 
        integer(I8)                     :: i, k 

        ! Check
        !======
        npvtot = 0
        do i = 1, vessel%polygonset%np 
            npvtot = npvtot + vessel%polygonset%polygons(i)%nv 
        end do
        if ( (npvtot /= size(xv, 1)) .or. (size(xv, 1) /= size(yv, 1)) ) then 
            ! Incompatible dimensions
            call gdErrorHandler('UpdateVesselCoordinates: incompatible ' // &
                'dimensions of new coordinates and original vessel polygon')
        end if  

        ! Adjust coordinates
        !===================
        ! Vessel polygon vertices
        k = 0
        do i = 1, vessel%polygonset%np
            ! Get number of vertices of this polygon
            npv = vessel%polygonset%polygons(i)%nv 

            ! Get the coordinates
            xvp = xv(k+1:k+npv)
            yvp = yv(k+1:k+npv)

            ! Assign
            call vessel%polygonset%polygons(i)%Construct(&
                xv(vessel%polygonset%polygons(i)%vert), &
                yv(vessel%polygonset%polygons(i)%vert), &
                vessel%polygonset%polygons(i)%labels(vessel%polygonset%polygons(i)%vert,:))
            
            ! Update counter
            k = k + npv 
        end do

        ! Reconstruct polygonset (just in case, shouldn't be necessary)
        call vessel%polygonset%Construct(vessel%polygonset%polygons)

        ! Test orientation
        call vessel%polygonset%OrientNestedClosedPolygons(flag)

        ! Check
        if (flag .ne. 0) then  
            ! Throw error
            print *, 'flag: ', flag
            call gdErrorHandler('UpdateVesselCoordinates: could not orient ' // &
                'polygons, OrientNestedClosedPolygons exited with flag above')
        end if 

        ! Write data
        vesselpath = 'vesselpolygon'
        call vessel%polygonset%WriteData(vesselpath)

        ! Adjust vessel description
        !==========================
        call vessel%plfvessel%Initialize(vessel%polygonset)

    end subroutine
    

end module