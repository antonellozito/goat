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
    use mod_precision
    use goatmod_types 
    !use goatmod_userinput
    use mod_dynamicarrays

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
        integer(I8), allocatable, dimension(:)  :: xpointID
        !integer(I8), allocatable, dimension(:)  :: spointID, opointID, tpointID, isprimaryxp, spointxpID, &
        !    divFc, spointdivID, tpointdivID
        !integer(I8), allocatable, dimension(:, :)   :: divFcP
        integer(I8)                             :: nxp
        !integer(I8)                             :: nsp, nop, ntp, &
        !    ndiv, ndivFc  
            
    contains

        ! Initialize
        procedure :: Initialize     => InitializeGAGridData
        
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

    contains

        ! Initialize
        procedure :: Initialize         => InitializeGAGrid

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
                GAvert%bx, GAvert%by, GAvert%ffbz)
        end if 
        allocate(RealDynamicArrayUDT::GAvert%x, GAvert%y, &
            GAvert%bx, GAvert%by, GAvert%ffbz)
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

        ! Deallocate grid - maybe no problem
        !call DeallocateVertices(grid)
        !call DeallocateFaces(grid)



    end subroutine

    subroutine TranslateGAGridTOGrid(grid,GAgrid)
        ! Description
        ! ===========
        ! Translating the information in the GridUDT type to a GAGridUDT type with dynamic arrays

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT) :: grid
        type(GAGridUDT) :: GAgrid

        ! Initialize grid

        ! Give information in GAgrid to grid

        ! Deallocate GAgrid




    end subroutine 

    !==================================================================!
    !                                                                  !
    !                           FUNCTIONS                              !
    !                                                                  !
    !==================================================================! 
    
    ! Get vertices of a cell and dynamic arrays
    function GetCellVertGA(cell, i) result(res)
        integer(I8)         :: i, j 
        type(GACellUDT)     :: cell
        integer(I8), allocatable :: res(:), range(:)
        range = (/ (j, j = cell%vertP1%Get(i), (cell%vertP1%Get(i) + cell%vertP2%Get(i) - 1)) /)
        res = cell%vert%GetMultipleElements(range)   
    end function

    ! Get cells of a face with dynamic arrays
    function GetFaceCellGA(cell, i, cvLookUp) result(res)
        integer(I8)                 :: i 
        type(GACellUDT)               :: cell
        integer(I8), allocatable, optional    :: cvLookUp(:)
        integer(I8), allocatable    :: res(:)  
        
        if (.not.present(cvLookUp)) then
            cvLookUp = GetCvLookUp(cell)
        end if

        res = pack(cvLookUp,cell%face%GetAllElements().eq.i)
    end function

    function GetFluxSurfaceFcs(fd, i) result(res)
        integer(I8) :: i, s, nf, j
        type(GAFluxDataUDT) :: fd
        integer(I8), allocatable :: res(:), range(:)

        nf = fd%fluxsurfacefacesP2%Get(i)
        s = fd%fluxsurfacefacesP1%Get(i)
        range = (/(j, j = s, (s+nf-1) )/)
        res = fd%fluxsurfacefaces%GetMultipleElements(range)

    end function    

    ! Get cells of a vertex without using vert%cell and dynamic arrays
    function GetVertCellGA(cell, i, cvLookUp) result(res)
        integer(I8)                 :: i 
        type(GACellUDT)               :: cell
        integer(I8), allocatable, optional    :: cvLookUp(:)
        integer(I8), allocatable    :: res(:)
        
        if (.not.present(cvLookUp)) then
            cvLookUp = GetCvLookUp(cell)
        end if
            
        res = pack(cvLookUp,cell%vert%GetAllElements.eq.i)
    end function

    function GetCvLookUp(cell) result(res)
        type(GACellUDT)       :: cell
        integer(I8)         :: nc, ic, nv, s, i               
        integer(I8), allocatable :: res(:), range(:), range2(:)

        nc = cell%ntot
        range = (/ (i, i = 1,(cell%vertP1%Get(nc)+cell%vertP2%Get(nc)-1))/)

        allocate(res(range))
        res = 0

        do ic = 1, nc
            s = cell%vertP1%Get(ic)
            nv = cell%vertP2%Get(ic)
            range = (/ (i, i = s, (s+nv-1)) /)
            res(range) = ic
        end do
    end function

end module 