!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains all the type definitions used in the grid 
! deformation module. See the option setting routines in gdmod_userinput
! for an explanation of the different fields and options. Additionally,
! precision types are set here that can (should) be used in all 
! subroutines. 

! Note: default accessibility is set to public, since all these types 
! and corresponding functions are free to be used. 

module gdmod_types

    ! Initialize
    !============
    ! The usual
    implicit none
    save
    public 

    ! Precision types
    integer, parameter       :: R4      = selected_real_kind(6)
    integer, parameter       :: R8      = selected_real_kind(14)
    integer, parameter       :: R16     = selected_real_kind(33)

    integer, parameter       :: C32     = 32
    integer, parameter       :: C64     = 64
    integer, parameter       :: C132    = 132

    integer, parameter       :: I4      = selected_int_kind(4)
    integer, parameter       :: I8      = selected_int_kind(8)
    integer, parameter       :: I16     = selected_int_kind(16)

    ! All private types
    
    ! All private functions

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                               Options                            !
    !------------------------------------------------------------------!
    ! Each option type has a setting routine called 
    ! SetDefault<optionname>, see subroutines after the 'contain' 
    ! statement. 

    ! Options for the main runfile type
    type RunfileOptionsUDT
        character(C32)          :: runtype ! type of run: 'optimize' or 'test'
        character(C32)          :: gridtype ! type of grid: 'plasma' 
        character(C32)          :: meth ! method for grid deformation: 'KKT'
        logical                 :: export ! do export? 
    end type

    ! Options for exporting
    type ExportOptionsUDT
        character(C32)       :: gridformat ! format of grid to export
    end type

    ! Options for design variables
    type DesignVariableOptionsUDT
        character(C32)       :: type ! type of design variables
    end type

    ! Options for the cost function
    type CostFunctionOptionsUDT
        character(C32)       :: type ! the cost function format
    end type

    ! Options for the constraints
    type ConstraintOptionsUDT
        ! Fields for equality constraints
        logical             :: fluxfunction ! impose constraints on flux
        logical             :: xpoints ! impose x-point location
        logical             :: edgelengths ! impose edge length cons
        logical             :: orthogonality 

        ! Fields for inequality constraints
        logical             :: linefolding ! prevent flux line folding
    end type

    ! Options for design optimization
    type DesignOptionsUDT
        type(CostFunctionOptionsUDT)       :: costfunction
        type(DesignVariableOptionsUDT)     :: variables
        type(ConstraintOptionsUDT)         :: constraints
    end type

    ! Options for the grid
    type GridOptionsUDT
        character(C32)          :: inputtype
        logical                 :: vesselrefine
        real(R8)                :: vesselmaxdist
    end type  

    ! Options for numerics
    type NumOptionsUDT
        integer(R8)     :: itmax 
    end type

    ! Options for magnetic field
    type MagneticFieldOptionsUDT
        character(C32)       :: readmeth
    end type


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

        ! Coordinates
        real(R8), allocatable               :: x(:),y(:) 

        ! Logicals and indices
        logical, allocatable                :: BV(:)
        integer(I8), allocatable            :: fieldlineID(:) 
        integer(I8)                         :: ntot

        integer(I8), allocatable            :: faceP(:,:)
        integer(I8), allocatable            :: facelist(:)
        integer(I8)                         :: nfacelist

        integer(I8), allocatable            :: cellP(:,:)
        integer(I8), allocatable            :: celllist(:)
        integer(I8)                         :: ncelllist

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

        ! Logicals and indices
        integer(I8), allocatable            :: vert(:,:)
        integer(I8), allocatable            :: neig(:,:)
        integer(I8)                         :: ntot
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


        ! Logicals and indices
        integer(I8), allocatable            :: vertP(:,:)
        integer(I8), allocatable            :: vertlist(:)
        integer(I8)                         :: nvertlist

        integer(I8), allocatable            :: faceP(:,:)
        integer(I8), allocatable            :: facelist(:)
        integer(I8)                         :: nfacelist

        integer(I8)                         :: ntot                      
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

        ! Logicals and indices
        integer(I8)                         :: nFt
        integer(I8)                         :: nFs

        ! Arrays, flux tube data
        integer(I8), allocatable            :: fluxtubecellsP(:,:)
        integer(I8), allocatable            :: fluxtubecells(:)
        integer(I8), allocatable            :: fluxtubefacesP(:,:)
        integer(I8), allocatable            :: fluxtubefaces(:)
        integer(I8), allocatable            :: cellfluxtubeID(:)

        ! Arrays, flux surface data
        integer(I8), allocatable            :: fluxsurfacefacesP(:,:)
        integer(I8), allocatable            :: fluxsurfacefaces(:)

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

        ! Vertices
        type(VertexUDT)                     :: vert

        ! Faces
        type(FaceUDT)                       :: faces

        ! Cells
        type(CellUDT)                       :: cells

        ! Additional data
        type(GridDataUDT)                   :: data

    end type

    !------------------------------------------------------------------!
    !                            Optimization                          !
    !------------------------------------------------------------------!

    ! Design parameter structure
    type DesignParamsUDT
        ! Design variables
        real, allocatable   :: phi(:) 

    end type

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
        ! Vertex values
        allocate(cells%vertP(cells%ntot,2))
        allocate(cells%vertlist(cells%nvertlist)) 

        ! Face values
        allocate(cells%faceP(cells%ntot,2))
        allocate(cells%facelist(cells%nfacelist))

        

    end subroutine

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

    end subroutine

    ! Deallocation


end module