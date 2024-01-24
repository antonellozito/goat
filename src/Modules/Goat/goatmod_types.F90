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

        ! Coordinates
        real(R8), allocatable               :: x(:),y(:) 

        ! Logicals and indices
        logical, allocatable                :: BV(:)
        integer(I8), allocatable            :: fieldlineID(:) 
        integer(I8)                         :: ntot = 0

        integer(I8), allocatable            :: faceP(:,:)
        integer(I8), allocatable            :: facelist(:)
        integer(I8)                         :: nfacelist = 0

        integer(I8), allocatable            :: cellP(:,:)
        integer(I8), allocatable            :: celllist(:)
        integer(I8)                         :: ncelllist = 0

        integer(I8), allocatable            :: neigP(:,:)
        integer(I8), allocatable            :: neiglist(:)
        integer(I8)                         :: nneiglist = 0

    end type

    ! Face structure
    type FaceUDT

        ! Description
        !============
        ! Fields:
        ! - ntot            : total number of faces
        ! - vert            : set of (two) vertices belonging to that 
        !                   face
        ! - neig            : cell neighbours of face
        ! - BF              : logical index that is true if the face
        !                   is a boundary face

        ! Logicals and indices
        integer(I8), allocatable            :: vert(:,:)
        integer(I8), allocatable            :: neig(:,:)
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
        !                   stored in cells%vertlist)            
        ! - vertlist        : list of cell vertices, to be queried as 
        !                   cells%vertlist(cells%vertP(i,1):
        !                   cells%vertP(i,1)+cells%vertP(i,2)-1) 
        ! - nvertlist       : length of vertlist
        ! - faceP           : similar to vertP, but for faces
        ! - facelist        : similar to vertlist, but for faces
        ! - nfacelist       : similar to nvertlist, but for faces
        ! - GC              : ntot-by-1 logical vector indicating if the
        !                   cell is a guard cell


        ! Logicals and indices
        integer(I8), allocatable            :: vertP(:,:)
        integer(I8), allocatable            :: vertlist(:)
        integer(I8)                         :: nvertlist = 0

        integer(I8), allocatable            :: faceP(:,:)
        integer(I8), allocatable            :: facelist(:)
        integer(I8)                         :: nfacelist = 0

        logical, allocatable                :: GC(:)

        integer(I8)                         :: ntot = 0                      
    end type

    ! Boundary structure
    type BndUDT

        ! Description
        !============
        ! Fields:   
        ! - nvert:                      total number of vertices
        ! - vert:                       all vertices belonging to that 
        !                               boundary
        ! - nfaces:                     total number of faces 
        ! - faces:                      all face indices 
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
        integer(I8)                         :: nfaces = 0
        integer(I8), allocatable            :: faces(:)             
        
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
        ! - cellfluxtubeID  : nCv(number of cells)-by-1 array containing
        !                   the flux tube number (ID) for each cell. 
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
        integer(I8), allocatable            :: cellfluxtubeID(:)

        ! Arrays, flux surface data
        integer(I8), allocatable            :: fluxsurfacefacesP(:,:)
        integer(I8), allocatable            :: fluxsurfacefaces(:)
        integer(I8), allocatable            :: fluxsurfaceID(:)

    end type

    ! Region data
    type RegionDataUDT

        ! Description
        !============
        ! Data type to collect all information on which cells/verts/face
        ! belongs to which grid region. 
        ! Fields:
        !
        ! - cellregID           : grid%cells%ntot-by-1 array containing
        !                       the region IDs for each cell
        ! - faceregID           : grid%faces%ntot-by-1 array containing
        !                       the region IDs for each face
        ! - fluxtuberegID       : fluxdata%nFt-by-1 array containing 
        !                       the region IDs for each flux tube

        ! Arrays
        integer(I8), allocatable            :: cellregID(:)
        integer(I8), allocatable            :: faceregID(:)
        integer(I8), allocatable            :: fluxtuberegID(:)
        integer(I8), allocatable            :: facelabel(:)
        integer(I8), allocatable            :: celllabel(:)

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

        ! Logicals & scalars
        integer(I8)             :: nx, ny
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

        ! Flux data
        type(FluxDataUDT)           :: fluxdata

        ! Region data
        type(RegionDataUDT)         :: regions

        ! Legacy data of structured grid
        type(StructuredGridDataUDT) :: sglegacy

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
        ! - faces           
        ! - cells
        ! - data
        ! - bnd
        !
        ! Note: the bnd substructure has to be allocated separately

        ! Vertices
        type(VertexUDT)                     :: vert

        ! Faces
        type(FaceUDT)                       :: faces

        ! Cells
        type(CellUDT)                       :: cells

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
    ! Bicubic spline interpolant
    type BicubicSplineInterpolantUDT

        ! Description
        !============
        ! Bicubic spline interpolation structure. Stores all necessary 
        ! data to compute interpolated quantities (data has to be given
        ! on a structured non-uniform grid) and their spatial 
        ! derivatives (only in 2D)

        ! Dimensions 
        integer(I8)                 :: nx, ny ! number of points in x, y direction
        integer(I8)                 :: nc ! number of 'cells' 

        ! Coordinates
        real(R8), allocatable       :: x(:), y(:) ! coordinate vectors

        ! Interpolation coefficients
        integer(I8), allocatable    :: cellindex(:) ! index of interpolant grid cells
        real(R8), allocatable       :: a(:,:) ! actual coefficient matrix

        ! Scaling constants of the interpolant
        real(R8), allocatable       :: refx(:), refy(:), refdx(:), refdy(:) 

    end type

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
        type(BicubicSplineInterpolantUDT)        :: interp

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
        ! - strcutures : nstructures-by-1 array of structures

        ! Target plates
        integer(I4)                         :: ntp = 0
        integer(I4), allocatable            :: allTPind(:)
        integer(I4), allocatable            :: TPind(:)

        ! Structures
        integer(I4)                         :: nstructures = 0
        type(VesselStructureUDT), allocatable       :: structures(:)
        type(PolygonSetUDT)                 :: polygonset 


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

        type(VesselUDT)                 :: vessel

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
        call AllocateFaces(grid%faces)

        ! Allocate cell data
        call AllocateCells(grid%cells)

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

        ! Face data
        allocate(vert%faceP(vert%ntot,2))
        allocate(vert%facelist(vert%nfacelist)) 

        ! Cell data
        allocate(vert%cellP(vert%ntot,2))
        allocate(vert%celllist(vert%ncelllist))

    end subroutine

    ! Face substructure
    subroutine AllocateFaces(faces)

        ! Description
        !============
        ! Allocate the fields in the faces structure. At least the 
        ! following scalar fields have to be present:
        !
        ! - ntot        : total number of faces


        ! The usual
        implicit none

        ! Declare variables
        type(FaceUDT)       :: faces

        ! Allocate
        !=========
        allocate(faces%vert(faces%ntot,2))
        allocate(faces%neig(faces%ntot,2)) 
        allocate(faces%BF(faces%ntot))

    end subroutine

    ! Cell substruture
    subroutine AllocateCells(cells)

        ! Description
        !============
        ! Allocate the fields in the cells structure. At least the 
        ! following scalar fields have to be present:
        !
        ! - ntot        : total number of cells
        ! - nvertlist   : length of cells%vertlist
        ! - nfacelist   : length of cells%facelist
        ! 
        ! The following fields are allocated: 

        ! Note: roughly the same conventions on data structures


        ! The usual
        implicit none

        ! Declare variables
        type(CellUDT)       :: cells

        ! Allocate
        !=========
        ! Guard cells
        allocate(cells%GC(cells%ntot))
        ! Vertex values
        allocate(cells%vertP(cells%ntot,2))
        allocate(cells%vertlist(cells%nvertlist)) 

        ! Face values
        allocate(cells%faceP(cells%ntot,2))
        allocate(cells%facelist(cells%nfacelist))

        

    end subroutine

    ! Boundary substructure
    subroutine AllocateBnd(bnd)

        ! Description
        !============
        ! Allocate the boundary fields. It is assumed that the following
        ! fields are present:
        !
        ! - nfaces:         number of boundary faces
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
        bnd%nvert = bnd%nfaces+1
        allocate(bnd%faces(bnd%nfaces))
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
        ! - grid%cells%ntot
        ! - grid%faces%ntot

        ! The usual
        implicit none

        ! Declare variables
        type(FluxDataUDT)       :: fluxdata
        type(GridUDT)           :: grid

        ! Allocate
        !=========
        ! Flux tube data
        allocate(fluxdata%fluxtubecellsP(fluxdata%nFt,2))
        allocate(fluxdata%fluxtubecells(grid%cells%ntot)) 
        allocate(fluxdata%fluxtubefacesP(fluxdata%nFt,2))
        allocate(fluxdata%fluxtubefaces(grid%faces%ntot)) 
        allocate(fluxdata%cellfluxtubeID(grid%cells%ntot))

        ! Flux surface data
        allocate(fluxdata%fluxsurfacefacesP(fluxdata%nFs,2))
        allocate(fluxdata%fluxsurfacefaces(grid%faces%ntot))
        allocate(fluxdata%fluxsurfaceID(grid%vert%ntot))

    end subroutine

    ! Region data substrucure
    subroutine AllocateRegionData(regions,grid)

        ! Description
        !============
        ! Allocat the region data grid substructure. The following 
        ! fields should be present:
        !
        ! - grid%cells%ntot
        ! - grid%faces%ntot
        ! - grid%data%fluxdata%nFt

        ! The usual
        implicit none

        ! Declare variables
        type(RegionDataUDT)         :: regions
        type(GridUDT)               :: grid

        ! Allocate
        !=========
        allocate(regions%cellregID(grid%cells%ntot))
        allocate(regions%faceregID(grid%faces%ntot))
        allocate(regions%fluxtuberegID(grid%data%fluxdata%nFt))
        allocate(regions%facelabel(grid%faces%ntot))
        allocate(regions%celllabel(grid%cells%ntot))

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
        deallocate(vert%facelist) 

        ! Cell data
        deallocate(vert%cellP)
        deallocate(vert%celllist)

    end subroutine

    !------------------------------------------------------------------!
    !                             Numerics                             !
    !------------------------------------------------------------------!
    subroutine AllocateBicubicSplineInterpolant(interp)

        ! Description
        !============
        ! Allocate the fields of the interpolant. It is assumed that nx
        ! and ny are set correctly. 

        ! The usual
        implicit none

        ! Declare variables
        type(BicubicSplineInterpolantUDT)   :: interp

        ! Allocate
        !=========
        interp%nc = (interp%nx-1)*(interp%ny-1)
        allocate(interp%x(interp%nx))
        allocate(interp%y(interp%ny))
        allocate(interp%a(interp%nc, 16))
        allocate(interp%refx(interp%nc))
        allocate(interp%refy(interp%nc))
        allocate(interp%refdx(interp%nc))
        allocate(interp%refdy(interp%nc))



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
    

end module