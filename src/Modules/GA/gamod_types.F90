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
        type(RealDynamicArrayUDT), allocatable        :: x(:), y(:)

        ! Logicals and indices
        logical, allocatable                          :: BV(:)
        type(IntegerDynamicArrayUDT), allocatable        :: fieldlineID(:) 
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
        type(RealDynamicArrayUDT), allocatable, dimension(:) :: bx, &
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
        type(IntegerDynamicArrayUDT), allocatable         :: vert(:,:)

        type(IntegerDynamicArrayUDT), allocatable         :: &
            label(:), reg(:), aligned(:)

        !logical :: aligned(:)
        !type(IntegerDynamicArrayUDT), allocatable         :: cell(:)    
        !type(IntegerDynamicArrayUDT), allocatable         :: cellP(:,:)
        integer(I8)                                       :: ncell = 0

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
        type(IntegerDynamicArrayUDT), allocatable        :: vertP(:,:)
        type(IntegerDynamicArrayUDT), allocatable        :: vert(:)
        integer(I8)                                      :: nvert = 0

        type(IntegerDynamicArrayUDT), allocatable        :: faceP(:,:)
        type(IntegerDynamicArrayUDT), allocatable        :: face(:)
        integer(I8)                                      :: nface = 0

        !logical, allocatable                             :: GC(:)

        integer(I8)                                      :: ntot = 0, ngc
        
        type(RealDynamicArrayUDT), allocatable, dimension(:) :: psi, &
            bp, bt, x, y
        type(IntegerDynamicArrayUDT), allocatable, dimension(:)  :: &
            cflags, reg, ft
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
        type(IntegerDynamicArrayUDT), allocatable :: &
            fluxsurfacefacesP(:,:)
        type(IntegerDynamicArrayUDT), allocatable :: &
            fluxsurfacefaces(:)
        type(IntegerDynamicArrayUDT), allocatable :: fluxsurfaceID(:)
        !type(IntegerDynamicArrayUDT), allocatable :: &
        !    fluxsurfaceneig(:), fluxsurfaceneigP(:, :)
        type(RealDynamicArrayUDT), allocatable    :: fluxsurfacepsi(:)
        type(IntegerDynamicArrayUDT), allocatable :: fluxsurfacevertsP(:,:)
        type(IntegerDynamicArrayUDT), allocatable :: fluxsurfaceverts(:)
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
        type(StructuredGridDataUDT) :: sglegacy

        ! Topological mesh type
        integer(I8)                             :: topoflag
        
        ! OMP & IMP
        integer(I8), allocatable, dimension(:)  :: OMPcell, OMPface, &
            IMPcell, IMPface
        integer(I8)                             :: nOMPcell, nOMPface, &
            nIMPcell, nIMPface
        real(R8), dimension(1:2)                :: OMPr, OMPz, IMPr, &
            IMPz

        ! X-point(s), strike points, o points, ...
        integer(I8), allocatable, dimension(:)  :: xpointID, &
            spointID, opointID, tpointID, isprimaryxp, spointxpID, &
            divFc, spointdivID, tpointdivID
        integer(I8), allocatable, dimension(:, :)   :: divFcP
        integer(I8)                             :: nxp, nsp, nop, ntp, &
            ndiv, ndivFc  
            
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

    subroutine InitializeGAGrid(grid, nv, nf, nc, ncv, ncf, &
        nfs,nft,nsv,nsf)

        ! Description
        !============
        ! Initialize the GAgrid substructures
        
        ! Declare variables
        !==================
        class(GAGridUDT) :: grid
        integer(I8) :: nv, nf, nc, ncv, ncf, nfs, nft, nsv, nsf

        ! Aux

        ! Initialize
        !===========

        ! Substructures
        call grid%vert%Initialize(nv)
        call grid%face%Initialize(nf)
        call grid%cell%Initialize(nc, ncv, ncf)
        call grid%data%Initialize(nfs,nft,nv,nsv,nsf)

    end subroutine

    subroutine InitializeGAVertex(GAvert, nv)

        ! Description
        !============
        ! Initialize the grid adaptation vertex structure (simply empty arrays)        

        ! Declare variables
        !==================
        ! Arguments
        class(GAVertexUDT) :: GAvert
        integer(I8), intent(in), optional :: nv
        ! Initialize
        !===========
        if (present(nv)) then 
            GAvert%ntot = nv 
        else
            GAvert%ntot = 0
        end if 
        if (allocated(GAvert%x)) then 
            ! Assume all allocated
            deallocate(GAvert%x, GAvert%y, GAvert%fieldlineID, &
                GAvert%bx, GAvert%by, GAvert%ffbz)
        end if 
        associate(ntot  => GAvert%ntot)
        allocate(GAvert%x(ntot), GAvert%y(ntot), GAvert%fieldlineID(ntot), &
                GAvert%bx(ntot), GAvert%by(ntot), GAvert%ffbz(ntot))

        end associate

    end subroutine

    subroutine InitializeGAFace(GAface,nf)

        ! Description
        !============
        ! Initialize the grid adapation face structure (simply empty arrays)

        ! Declare variables
        !==================
        ! Arguments
        class(GAFaceUDT)      :: GAface
        integer(I8), intent(in), optional   :: nf

        ! Initialize
        !===========
        if (present(nf)) then 
            GAface%ntot = nf
        else
            GAface%ntot = 0
        end if
        associate(ntot  => GAface%ntot)
        if (allocated(GAface%vert)) then 
            ! assume all allocated
            deallocate(GAface%vert, GAface%label, GAface%reg, &
                GAface%aligned)
        end if 
        allocate(GAface%vert(ntot,2), GAface%label(ntot), GAface%reg(ntot), &
                GAface%aligned(ntot))

        ! Housekeeping
        end associate

    end subroutine

    subroutine InitializeGACell(GAcell,nc,ncv,ncf)

        ! Description
        !============
        ! Initialize the grid adapation cell structure (simply empty arrays)

        ! Declare variables
        !==================
        ! Arguments
        class(GACellUDT)      :: GAcell
        integer(I8), intent(in), optional   :: nc, ncv, ncf

        ! Initialize
        !===========
        if (present(nc)) then 
            GAcell%ntot = nc
        else
            GAcell%ntot = 0
        end if
        if (present(ncv)) then 
            GAcell%nvert = ncv 
        else 
            GAcell%nvert = 0
        end if 
        if (present(ncf)) then 
            GAcell%nface = ncf 
        else
            GAcell%nface = 0
        end if 
        associate(ntot  => GAcell%ntot, ntotf => GAcell%nface, &
            ntotv => GAcell%nvert &
            )
        if (allocated(GAcell%vertP)) then 
            ! assume all allocated
            deallocate(GAcell%vertP, GAcell%vert, GAcell%faceP, &
                GAcell%face, GAcell%psi, GAcell%bp, GAcell%bt, &
                GAcell%x, GAcell%y, GAcell%cflags, GAcell%reg, &
                GAcell%ft )
        end if 
        allocate(GAcell%vertP(ntot,2), GAcell%vert(ntotv), GAcell%faceP(ntot,2), &
                GAcell%face(ntotf), GAcell%psi(ntot), GAcell%bp(ntot), &
                 GAcell%bt(ntot), GAcell%x(ntot), GAcell%y(ntot), &
                 GAcell%cflags(ntot), GAcell%reg(ntot), GAcell%ft(ntot) )

        ! Housekeeping
        end associate

    end subroutine

    subroutine InitializeGAGridData(GAgriddata,nfs,nft,nv,nsv,nsf)

        ! Description
        !============
        ! Initialize the grid adapation data structure (simply empty arrays)

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridDataUDT)      :: GAgriddata
        integer(I8)               :: nfs,nft,nv,nsv,nsf

        ! Initialize
        call GAgriddata%fluxdata%Initialize(nfs,nft,nv,nsv,nsf)     

    end subroutine

    subroutine InitializeGAFluxData(GAfluxdata,nfs,nft,nv,nsv,nsf)

        ! Description
        !============
        ! Initialize the grid adapation data structure (simply empty arrays)

        ! Declare variables
        !==================
        ! Arguments
        class(GAFluxdataUDT)      :: GAfluxdata
        integer(I8), intent(in), optional   :: nfs, nft, nv, nsv, nsf 

        ! Auxiliary
        integer(I8) :: vntot, nsv1, nsf1

        ! Initialize
        !===========
        if (present(nfs)) then 
            GAfluxdata%nFs = nfs
        else
            GAfluxdata%nFs = 0
        end if
        if (present(nft)) then
            GAfluxdata%nFt = nft
        else
            GAfluxdata%nft = 0
        end if
        if (present(nv)) then
            vntot = nv
        else
            vntot = 0
        end if
        if (present(nsv)) then
            nsv1 = nsv
        else
            nsv1 = 0
        end if 
        if (present(nsf)) then
            nsf1 = nsf
        else
            nsf1 = 0
        end if
        associate(nfs  => GAfluxdata%nfs, nft => GAfluxdata%nft)
        if (allocated(GAfluxdata%fluxsurfacefaces)) then 
            ! assume all allocated
            deallocate(GAfluxdata%fluxsurfacefacesP, GAfluxdata%fluxsurfacefaces, &
                GAfluxdata%fluxsurfaceID, GAfluxdata%fluxsurfacevertsP, &
                GAfluxdata%fluxsurfaceverts, GAfluxdata%fluxsurfacepsi)
        end if 
        allocate(GAfluxdata%fluxsurfacefacesP(nfs,2), GAfluxdata%fluxsurfacefaces(nsf1), &
                GAfluxdata%fluxsurfaceID(vntot), GAfluxdata%fluxsurfacevertsP(nfs,2), &
                GAfluxdata%fluxsurfaceverts(nsv1), GAfluxdata%fluxsurfacepsi(nfs))

        ! Housekeeping
        end associate

    end subroutine

    subroutine TranslateGridTOGAGrid(grid,GAgrid)
        ! Description
        ! ===========
        ! Translating the information in the GridUDT type to a GAGridUDT type with dynamic arrays

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT) :: grid
        type(GAGridUDT) :: GAgrid

        ! Auxiliary
        integer(I8) :: nv, nf, nc, ncv, ncf, nfs, nft, nsv, nsf

        ! Initialize GAGrid
        nv = grid%vert%ntot
        nf = grid%face%ntot
        nc = grid%cell%ntot
        ncv = grid%cell%vertP(nc,1) + grid%cell%vertP(nc,2) - 1
        ncf = grid%cell%faceP(nc,1) + grid%cell%faceP(nc,2) - 1
        nfs = grid%data%fluxdata%nFs
        nft = grid%data%fluxdata%nFt
        nsv = grid%data%fluxdata%fluxsurfacevertsP(nfs,1) + grid%data%fluxdata%fluxsurfacevertsP(nfs,2) - 1
        nsf = grid%data%fluxdata%fluxsurfacefacesP(nfs,1) + grid%data%fluxdata%fluxsurfacefacesP(nfs,2) - 1
        call GAgrid%Initialize(nv, nf, nc, ncv, ncf, &
        nfs, nft, nsv, nsf)

        ! Give information in grid to GAgrid
        GAgrid%vert%x = ConstructRealDynamicArray(grid%vert%x)


        ! Deallocate grid


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

end module 