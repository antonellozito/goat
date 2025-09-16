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
    use mod_polygon
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

    type :: QualityMetricUDT

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

        ! Indicator for merging and splitting
        integer(I8)           :: merge_fc
        integer(I8)           :: split_cv

    contains

        procedure :: Initialize
        procedure :: CalculateQualityMetrics
        procedure :: ComputeQM
        procedure :: SelectMergingFace
        procedure :: SelectSplitCell

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
        class(RealDynamicArrayBufferedUDT), allocatable        :: x, y

        ! Logicals and indices
        !logical, allocatable                          :: BV(:)
        !class(IntegerDynamicArrayUDT), allocatable    :: fieldlineID 
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
        class(RealDynamicArrayBufferedUDT), allocatable :: bx, &
            by, psi
        real(R8) :: ffbz 
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
        class(IntegerDynamicArrayBufferedUDT), allocatable         :: vert1
        class(IntegerDynamicArrayBufferedUDT), allocatable         :: vert2
        class(IntegerDynamicArrayBufferedUDT), allocatable         :: &
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
        procedure :: CheckFace

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
        class(IntegerDynamicArrayBufferedUDT), allocatable        :: vertP1
        class(IntegerDynamicArrayBufferedUDT), allocatable        :: vertP2
        class(IntegerDynamicArrayBufferedUDT), allocatable        :: vert
        integer(I8)                                       :: nvert = 0

        class(IntegerDynamicArrayBufferedUDT), allocatable        :: faceP1
        class(IntegerDynamicArrayBufferedUDT), allocatable        :: faceP2
        class(IntegerDynamicArrayBufferedUDT), allocatable        :: face
        integer(I8)                                       :: nface = 0

        !logical, allocatable                             :: GC(:)

        integer(I8)                                       :: ntot = 0, ngc
        
        class(RealDynamicArrayBufferedUDT), allocatable           :: psi, x, y
        class(IntegerDynamicArrayBufferedUDT), allocatable        :: cflags, reg
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
        class(IntegerDynamicArrayBufferedUDT), allocatable :: fluxsurfacefacesP1
        class(IntegerDynamicArrayBufferedUDT), allocatable :: fluxsurfacefacesP2
        class(IntegerDynamicArrayBufferedUDT), allocatable :: fluxsurfacefaces
        !class(IntegerDynamicArrayUDT), allocatable :: fluxsurfaceID
        !type(IntegerDynamicArrayUDT), allocatable :: &
        !    fluxsurfaceneig(:), fluxsurfaceneigP(:, :)
        !class(RealDynamicArrayUDT), allocatable    :: fluxsurfacepsi
        class(IntegerDynamicArrayBufferedUDT), allocatable :: fluxsurfacevertsP1
        class(IntegerDynamicArrayBufferedUDT), allocatable :: fluxsurfacevertsP2
        class(IntegerDynamicArrayBufferedUDT), allocatable :: fluxsurfaceverts

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


        ! Pre- and PostProcessing
        procedure :: CheckVertOrder
        procedure :: ReorderCellConn
        procedure :: GetFsVxFromFsFc
        procedure :: GiveXpoints
        procedure :: GiveSeparatrices
        procedure :: IdentifyAlignedFaces
        procedure :: CheckUnstructuredGrid
        procedure :: RecalcMagn
        procedure :: MergeFS

        ! Grid adaptation
        !================

        ! Removal Small Triangles
        procedure :: RemoveSmallTriangle
        procedure :: DetectSmallTrias
        procedure :: DetermineRemovalMethod
        procedure :: LocalSmallTriangleRemoval

        ! Stacked triangles
        procedure :: StackedTrias
        procedure :: DetectCutCell
        procedure :: GetFirstTrap
        procedure :: GetNextTraps
        procedure :: GetConnectionVertex
        procedure :: CheckStackedTriaOverlap
        procedure :: StackAdaptation

        ! Stacked to cutcell
        procedure :: StackedToCutcell
        procedure :: DetectStackedTrias

        ! Merging
        procedure :: DoMerging
        procedure :: MergeCells
        procedure :: OneMerge
        procedure :: MergeRec
        procedure :: DetermineMergeCaseID
        procedure :: Merge4444
        procedure :: Merge3443
        procedure :: Merge4443
        procedure :: Merge4433
        procedure :: Merge334
        procedure :: MakeMergePent
        procedure :: MakeNewThreeFace
        procedure :: AdaptNeigThreeVert

        ! Splitting 
        procedure :: DoSplitting
        procedure :: Splitting
        procedure :: OneSplit
        procedure :: SplitPentsRec
        procedure :: DetermineTcaseID
        procedure :: DetermineQcaseID
        procedure :: DeterminePcaseID
        procedure :: DetermineT4caseID
        procedure :: SelectSplitPent
        procedure :: GetRadSplitFacePent
        procedure :: GetNeigTypeFromFace
        procedure :: GetIntersectedPsiFaces
        procedure :: GetPents
        procedure :: SplitTQ
        procedure :: SplitTP
        procedure :: SplitQTB
        procedure :: SplitBTB
        procedure :: SplitTTT
        procedure :: SplitTTQ
        procedure :: SplitQQQ
        !procedure :: SplitQQT
        !procedure :: SplitBQT
        !procedure :: SplitTQT
        !procedure :: SplitQshaved
        procedure :: SplitPQ
        procedure :: SplitPT
        procedure :: SplitT4Q
        procedure :: SplitT4B
        procedure :: SplitFace
        procedure :: SplitTriaStacked
        procedure :: SplitCenterTria
        procedure :: SplitCenterQuad
        procedure :: SplitT4
        procedure :: SplitT4Stacked
        procedure :: TriaToQuad
        procedure :: QuadToPent
        procedure :: SplitCenterPent
        procedure :: DeterminePerpFaceTria
        procedure :: DeterminePerpFaceQuad
        procedure :: DeterminePerpFacePent
        procedure :: DetermineFreeVertTria
        procedure :: DetermineHangingNodePent


        ! Boundary layer grid
        procedure :: BoundaryLayerGrid


        ! Grid information utility
        procedure :: GetForbiddenMergeFaces
        procedure :: GetCutsXpoints
        procedure :: GetTangencyPoints
        procedure :: GetRadLineFaces
        procedure :: RecursiveGridMarching
        procedure :: DetermineCflags
        procedure :: AreaConstraintPents
        procedure :: GetXcellsGeometric

        ! Grid operation
        procedure :: RemoveCells
        procedure :: RemoveFaces
        procedure :: RemoveVertices
        procedure :: GetFaceNumber
        procedure :: AddFaceToFsFc
        procedure :: AddVertToFsVx
        procedure :: AddCell
        procedure :: AddVert
        procedure :: AddFluxSurface
        procedure :: AddIntoExistingSurface

        ! Computing
        !===========
        procedure :: CalcHpol0D
        generic   :: CalcHpol => CalcHpol0D  
        procedure :: CalcHrad0D
        generic   :: CalcHrad => CalcHrad0D         
        procedure :: CalcCentroid0DGA
        procedure :: CalcCentroid1DGA
        generic   :: CalcCentroidGA => CalcCentroid0DGA, CalcCentroid1DGA

        ! Visualization
        !==============
        procedure :: WriteData => WriteGAGridData



    end type  

    !==================================================================!
    !                                                                  !
    !                          INTERFACES                              !
    !                                                                  !
    !==================================================================!

    ! General isBoundaryCell
    interface isBoundaryCellGA 
        module procedure isBoundaryCell0DGA, isBoundaryCell1DGA
    end interface

    ! General isBoundaryFace
    interface isBoundaryFaceGA
        module procedure isBoundaryFace0DGA, isBoundaryFace1DGA
    end interface

    ! General isBoundaryCell
    interface isBoundaryVertGA
        module procedure isBoundaryVert0DGA, isBoundaryVert1DGA
    end interface

    ! General norm
    interface Norm
        module procedure Norm0D, Norm1D
    end interface



    contains 
    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                  QUALITY METRIC COMPUTATIONS                     !
    !------------------------------------------------------------------!
    subroutine Initialize(qm, grid)

        ! Description
        !============
        ! Initialize quality metric type

        ! Declare variables
        !==================
        class(QualityMetricUDT)    :: qm
        type(GAGridUDT)         :: grid

        ! Initialize the qm arrays?? - TODO
        if (allocated(qm%fcBias)) then
            deallocate(qm%fcBias)
            deallocate(qm%fcqalfc)
            deallocate(qm%fcS)
            deallocate(qm%cvS)
            deallocate(qm%cvAR)
            deallocate(qm%h_pol)
            deallocate(qm%h_rad)
            deallocate(qm%h_rad_psi)
        end if
        allocate(qm%fcBias(grid%face%ntot))
        allocate(qm%fcqalfc(grid%face%ntot))
        allocate(qm%fcS(grid%face%ntot))
        allocate(qm%cvS(grid%cell%ntot))
        allocate(qm%cvAR(grid%cell%ntot))
        allocate(qm%h_pol(grid%cell%ntot))
        allocate(qm%h_rad(grid%cell%ntot))
        allocate(qm%h_rad_psi(grid%cell%ntot))

        ! Set to zero
        qm%fcBias       = 0
        qm%fcqalfc      = 0
        qm%fcS          = 0
        qm%cvS          = 0
        qm%cvAR         = 0
        qm%h_pol        = 0
        qm%h_rad        = 0
        qm%h_rad_psi    = 0

    end subroutine

    subroutine ComputeQM(qm,grid,options, magneticField)

        ! Description
        !============
        ! Compute all quality metrics
        ! (Mirror of CalculateCvMetric.m)

        ! Declare variables
        !==================
        ! Arguments
        class(QualityMetricUDT), intent(inout)  :: qm
        type(GAGridUDT), intent(in)             :: grid
        type(GAoptionsUDT), intent(in)          :: options
        type(MagneticFieldUDT), intent(in)      :: magneticField

        ! Auxiliary
        real(R8) :: fcH(grid%face%ntot,2), psic, &
            isx(4), isy(4), min_vpsi, max_vpsi, &
            h_rad_int, v1p, v2p, cvx1, cvy1, cvx2, cvy2, t0
        real(R8), allocatable, dimension(:) :: v1x, v1y, &
            v2x, v2y, fcX, fcY, fcS, fcBx, fcBy, fcs_fcs, &
            vx, vy, vertsX, vertsY, vpsi, cx, cy, h_pol, &
            h_rad, h_rad_psi, cvS,  cvAR, fcBias, psi_verts,  &
            fcqalfc, fcxx
        real(R8), allocatable :: fcBb(:,:), vec_n(:,:)
        integer(I8), allocatable, dimension(:) :: vx1, vx2, &
            tv, tf, tc, indCv, ncpf
        integer(I8) :: i, ic, nv, ii, ifc, ir, v1, v2
        integer(I8), allocatable :: fccv(:,:) 
        logical :: not_aligned_f(grid%face%ntot), &
            comp_range(grid%cell%ntot)



        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        ! Initialize
        allocate(vx1(f%ntot), vx2(f%ntot), v1x(f%ntot), v1y(f%ntot), &
            v2x(f%ntot), v2y(f%ntot), fcX(f%ntot), fcY(f%ntot), &
            vec_n(f%ntot,2), vx(v%ntot), vy(v%ntot), vpsi(v%ntot), &
            h_pol(c%ntot), h_rad(c%ntot), h_rad_psi(c%ntot), &
            cvS(c%ntot), cvAR(c%ntot), fcBias(f%ntot), fccv(f%ntot,2), &
            indCv(c%ntot), fcqalfc(f%ntot), ncpf(f%ntot), fcxx(f%ntot), &
            fcBx(f%ntot), fcBy(f%ntot),fcBb(f%ntot,2))
        ncpf = 0
        fccv = 0
        fcxx = 0

        vx1 = f%vert1%Get()
        vx2 = f%vert2%Get()
        vx = v%x%Get()
        vy = v%y%Get()
        vpsi = v%psi%Get()
        cx = c%x%Get()
        cy = c%y%Get()

        ! Calculate face properties
        v1x = v%x%Get(vx1)
        v1y = v%y%Get(vx1)
        v2x = v%x%Get(vx2)
        v2y = v%y%Get(vx2)

        fcX = 0.5_R8 * (v1x + v2x)
        fcY = 0.5_R8 * (v1y + v2y)

        fcS = sqrt( (v2x - v1x)**2 + (v2y - v1y)**2 )

        ! Normal vectors
        vec_n(:,1) = (v2y - v1y) / fcS
        vec_n(:,2) = -(v2x - v1x) / fcS

        call magneticField%interp%Evaluate(fcX,fcY,1,0,fcBx)
        call magneticField%interp%Evaluate(fcX,fcY,0,1,fcBy)

        fcBb(:,1) = -fcBy / sqrt(fcBx**2 + fcBy**2)
        fcBb(:,2) =  fcBx / sqrt(fcBx**2 + fcBy**2)

        fcqalfc = (fcBb(:,1)*vec_n(:,1) + fcBb(:,2)*vec_n(:,2));

        ! Pre-process non aligned faces
        !==============================
        not_aligned_f = (f%aligned%Get() == 0)

        indCv = (/ (i, i = 1, c%ntot )/)
        comp_range = indCv .le. qm%nCv

        ! Cell loop 
        !==========
        do ic = 1, c%ntot

            tv = GetCellVertGA(c, ic)
            tf = GetCellFaceGA(c, ic)
            nv = size(tv)

            if (any(tv .gt. size(vpsi)) .or. any(tv .lt. 1)) then
                print * , tv
            end if
            psi_verts = vpsi(tv)
            max_vpsi = maxval(psi_verts)
            min_vpsi = minval(psi_verts)

            ! Add faces to an array
            ncpf(tf) = ncpf(tf) + 1
            do i = 1, nv
                fccv(tf(i), ncpf(tf(i))) = ic
            end do

            ! Aspect ratio
            fcs_fcs = fcs(tf)
            cvAR(ic) = maxval(fcs_fcs) / minval(fcs_fcs)

            if (comp_range(ic) .and. .not.options%slab) then

                if (.not.(cvAR(ic) == qm%cvAR(ic))) then

                    !TODO - group in :: call CalcCvMetric(ic,tv,vx,vy,max_vpsi,min_vpsi,tf,vx1,vx2)

                    ! Area
                    !=====
                    vertsX = vx(tv)
                    vertsY = vy(tv)
                    if (nv == 4) then ! Quad

                        cvS(ic) = TriangleArea(vertsX(1),vertsY(1),vertsX(2),vertsY(2),vertsX(3),vertsY(3))  &
                        + TriangleArea(vertsX(1),vertsY(1),vertsX(3),vertsY(3),vertsX(4),vertsY(4))

                    elseif (nv == 3) then

                        cvS(ic) = triangleArea(vertsX(1),vertsY(1),vertsX(2),vertsY(2),vertsX(3),vertsY(3))

                    elseif (nv == 2) then

                        cvS(ic) = 0

                    elseif (nv == 5) then

                        cvS(ic) = TriangleArea(vertsX(1),vertsY(1),vertsX(2),vertsY(2),vertsX(3),vertsY(3)) &
                        + TriangleArea(vertsX(3),vertsY(3),vertsX(4),vertsY(4),vertsX(5),vertsY(5)) & 
                        + TriangleArea(vertsX(3),vertsY(3),vertsX(1),vertsY(1),vertsX(5),vertsY(5))

                    else

                        call gdErrorHandler('ComputeQM: cell type not implemented')

                    end if

                    ! H pol
                    !======
                    psic = 0.5_R8* (max_vpsi + min_vpsi) 

                    isx = 0
                    isy = 0
                    ii = 0
                    ir = 0
                    h_rad_int = 0
                    do i = 1, nv
                        ifc = tf(i)
                        v1 = vx1(ifc)
                        v2 = vx2(ifc)
                        v1p = vpsi(v1)
                        v2p = vpsi(v2)
                        if ((psic .gt. min(v1p, v2p)) &
                            .and. (psic .lt. max(v1p, v2p)) ) then
                            
                            t0 = (psic - v1p) / (v2p - v1p)
                            ii = ii + 1
                            isx(ii) = vx(v1) + t0 *(vx(v2) - vx(v1))
                            isy(ii) = vy(v1) + t0 *(vy(v2) - vy(v1))

                            if (not_aligned_f(ifc)) then
                                ir = ir + 1
                                h_rad_int = h_rad_int + fcS(ifc)*abs(fcqalfc(ifc))
                            end if
                            
                        end if 
                    end do

                    if (ii /= 2) call gdErrorHandler('ComputeQM: insufficient intersection with mean psi value found')

                    h_pol(ic) = sqrt( ( isx(2) - isx(1) )**2 + ( isy(2) - isy(1) )**2 )
                    h_rad(ic) = h_rad_int / ir
                    h_rad_psi(ic) = abs(max_vpsi - min_vpsi)

                else

                    ! Is aspect ratios are the same, the metric are not recalculated
                    h_pol(ic) = qm%h_pol(ic);
                    h_rad(ic) = qm%h_rad(ic);
                    h_rad_psi(ic) = qm%h_rad_psi(ic);
                    cvS(ic) = qm%cvS(ic);


                end if

            else 

                ! For new cells
                ! Area
                !=====
                vertsX = vx(tv)
                vertsY = vy(tv)
                if (nv == 4) then ! Quad

                    cvS(ic) = TriangleArea(vertsX(1),vertsY(1),vertsX(2),vertsY(2),vertsX(3),vertsY(3))  &
                    + TriangleArea(vertsX(1),vertsY(1),vertsX(3),vertsY(3),vertsX(4),vertsY(4))

                elseif (nv == 3) then

                    cvS(ic) = triangleArea(vertsX(1),vertsY(1),vertsX(2),vertsY(2),vertsX(3),vertsY(3))

                elseif (nv == 2) then

                    cvS(ic) = 0

                elseif (nv == 5) then

                    cvS(ic) = TriangleArea(vertsX(1),vertsY(1),vertsX(2),vertsY(2),vertsX(3),vertsY(3)) &
                    + TriangleArea(vertsX(3),vertsY(3),vertsX(4),vertsY(4),vertsX(5),vertsY(5)) & 
                    + TriangleArea(vertsX(3),vertsY(3),vertsX(1),vertsY(1),vertsX(5),vertsY(5))

                else

                    call gdErrorHandler('ComputeQM: cell type not implemented')

                end if

                ! H pol
                !======
                psic = 0.5_R8* (max_vpsi + min_vpsi) 

                isx = 0
                isy = 0
                ii = 0
                ir = 0
                h_rad_int = 0
                do i = 1, nv
                    ifc = tf(i)
                    v1 = vx1(ifc)
                    v2 = vx2(ifc)
                    v1p = vpsi(v1)
                    v2p = vpsi(v2)
                    if ((psic .gt. min(v1p, v2p)) &
                        .and. (psic .lt. max(v1p, v2p)) ) then
                        
                        t0 = (psic - v1p) / (v2p - v1p)
                        ii = ii + 1
                        isx(ii) = vx(v1) + t0 *(vx(v2) - vx(v1))
                        isy(ii) = vy(v1) + t0 *(vy(v2) - vy(v1))

                        if (not_aligned_f(ifc)) then
                            ir = ir + 1
                            h_rad_int = h_rad_int + fcS(ifc)*abs(fcqalfc(ifc))
                        end if
                        
                    end if 
                end do

                if (ii /= 2) call gdErrorHandler('ComputeQM: insufficient intersection with mean psi value found')

                h_pol(ic) = sqrt( ( isx(2) - isx(1) )**2 + ( isy(2) - isy(1) )**2 )
                h_rad(ic) = h_rad_int / ir
                h_rad_psi(ic) = abs(max_vpsi - min_vpsi)
            end if


        end do

        ! Face loop
        do ifc = 1, f%ntot

            tc = fccv(ifc,:)
            cvx1 = cx(tc(1))
            cvy1 = cy(tc(1))
            fcH(ifc,1) = sqrt( (fcX(ifc) - cvx1)**2 + (fcY(ifc) - cvy1)**2  )
            if (ncpf(ifc) == 2) then
                cvx2 = cx(tc(2));
                cvy2 = cy(tc(2))
                fcH(ifc,2) = sqrt( (fcX(iFc) - cvx2)**2 + (fcY(iFc) - cvy2)**2  )
                fcxx(ifc) = sqrt( (cvx1 - cvx2)**2 + (cvy1 - cvy2)**2  )
            else
                fcxx(ifc) = sqrt( (cvx1 - fcX(ifc))**2 + (cvy1 - fcY(ifc))**2  )
            end if

            fcBias(ifc) = maxval(fcH(ifc,:)) / minval(fcH(ifc,:))

        end do

        ! Initialize the qm arrays?? - TODO
        if (allocated(qm%fcBias)) then
            deallocate(qm%fcBias)
            deallocate(qm%fcqalfc)
            deallocate(qm%fcS)
            deallocate(qm%cvS)
            deallocate(qm%cvAR)
            deallocate(qm%h_pol)
            deallocate(qm%h_rad)
            deallocate(qm%h_rad_psi)
        end if
        allocate(qm%fcBias(grid%face%ntot))
        allocate(qm%fcqalfc(grid%face%ntot))
        allocate(qm%fcS(grid%face%ntot))
        allocate(qm%cvS(grid%cell%ntot))
        allocate(qm%cvAR(grid%cell%ntot))
        allocate(qm%h_pol(grid%cell%ntot))
        allocate(qm%h_rad(grid%cell%ntot))
        allocate(qm%h_rad_psi(grid%cell%ntot))

        ! Save computed data
        qm%fcBias   = fcBias
        qm%fcqalfc  = fcqalfc
        qm%fcS      = fcS
        qm%cvS      = cvS
        qm%cvAR     = cvAR
        qm%h_pol    = h_pol
        qm%h_rad    = h_rad
        qm%h_rad_psi= h_rad_psi
        qm%nCv      = c%ntot




        
        end associate



        
    end subroutine

    subroutine CalculateQualityMetrics(qm,grid,options,magneticField, select_split, select_merge)
        ! Description
        !============
        ! Compute metric and criteria of cells necessary to execute grid adaptation.

        ! Declare variables
        !==================
        ! Arguments
        class(QualityMetricUDT), intent(inout)  :: qm
        type(GAGridUDT), intent(inout)          :: grid
        type(GAoptionsUDT), intent(in)          :: options
        type(MagneticFieldUDT), intent(in)      :: magneticField
        logical :: select_split, select_merge


        ! Calculate cv metric
        call qm%ComputeQM(grid,options,magneticField)


        ! Selecting splitting cell
        if (select_split) then
            call qm%SelectSplitCell(grid, options)
        end if

        ! Selecting merging face
        if (select_merge) then
            call qm%SelectMergingFace(grid, options)
        end if
        


    end subroutine

    subroutine SelectMergingFace(qm, grid, options)

        ! Description
        !============
        ! Select a face to merge the two cell of

        ! Declare variables
        !==================
        ! Arguments
        class(QualityMetricUDT)     :: qm
        type(GAGridUDT)             :: grid
        type(GAoptionsUDT)          :: options

        ! Auxiliary 
        integer(I8) :: i, j, k, neig, nf, fcs_al1, fcs_al2, ind, &
            common_face, counter
        integer(I8), allocatable, dimension(:) :: forbidden_fcs, &
            indsort, cells, small_cells, fcs, cvs, cvLookUp, indcv, &
            no_cells, trias, cells2, indfc, pol_faces, fcs1, fcs2, cellsD, &
            ind_sort, neigs, vxs, cvs_sep, cvs_sepU, core_faces, Xcells, &
            core_facesD, rad_facesD, rad_faces
        real(R8) :: crit, dfunv(grid%cell%ntot), h_pol_no_cells_crit, &
            mean_pol_flux, bench
        real(R8), allocatable, dimension(:) :: area_small_cells, &
            h_pol_no_cells_sorted, h_pol_cells, h_pol_cvs, bias, &
            pol_fluxdens_est, pol_fluxdens_estD, h_rad_cells
        logical, allocatable :: log(:), log2(:), trias_log(:)

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            fd => grid%data%fluxdata &
            )

        ! Determine forbidden merge faces iFs of separatrix
        if (.not.options%slab) call grid%GetForbiddenMergeFaces(forbidden_fcs)

        ! Initialize
        qm%merge_fc = 0
        cvLookUp = GetCvLookUpGA(c)

        ! Criteria
        select case (options%merge_crit)

        case ('tria_to_quad')

            ! NOT IN MATLAB
            ! Determine merge face (case for when a quad consists of two triangles)
            ! Aligned face in other flux surface, or aligned faces have no common vertex
            indcv = (/ (i, i = 1, c%ntot)/)
            allocate(trias(count(c%faceP2%Get() == 3)))
            trias = pack(indcv, c%faceP2%Get() == 3)
            do i = 1, size(trias)

                ! Get faces and check if it is not aligned and has another triangle as neighbor
                fcs = GetCellFaceGA(c, trias(i))
                do j = 1, 3

                    if (f%aligned%Get(fcs(j)) == 0) then

                        cvs = GetFaceCellGA(c, fcs(j), cvLookUp)
                        if (size(cvs) == 2) then
                            if ((c%faceP2%Get(cvs(1)) == 3) .and. (c%faceP2%Get(cvs(2)) == 3)) then
                                ! Every triangle should have one aligned face
                                ! In this case the merge can happend if those aligned faces 
                                ! do not have a common vertex
                                fcs1 = GetCellFaceGA(c, cvs(1))
                                fcs2 = GetCellFaceGA(c, cvs(2))
                                fcs_al1 = 0
                                fcs_al2 = 0
                                do k = 1, 3
                                    if (f%aligned%Get(fcs1(k)) == 1) fcs_al1 = fcs1(k)
                                    if (f%aligned%Get(fcs2(k)) == 1) fcs_al2 = fcs2(k)
                                end do

                                ! If faces are found and have no common vertex, the trias can be merged
                                if (fcs_al1 /= 0 .and. fcs_al2 /= 0) then

                                    if (.not.HaveCommonVert(f, fcs_al1, fcs_al2)) then
                                        qm%merge_fc = fcs(j)
                                        exit
                                    end if

                                end if

                            end if

                        end if

                    end if

                end do

                if (qm%merge_fc /= 0) exit

            end do


        case ('min_area')

            ! Initialize criteria
            crit = sum(qm%cvS) / c%ntot

            ! Distance function 
            call grid%fun%distr%Evaluate(c%x%Get(), c%y%Get(), dfunv)
            log = (dfunv .lt. options%dist_function_threshold_merge)
            indcv = (/ (i, i = 1, c%ntot)/)
            allocate(cells(count(log)))
            cells = pack(indcv, log)

            ! Get small cells
            log2 = ((qm%cvS(cells) < crit) .and. (c%reg%Get(cells) /= 3) .and. (c%reg%Get(cells) /= 4))
            allocate(small_cells(count(log2)))
            small_cells = pack(cells, log2)

            ! Sort small cells for area
            area_small_cells = qm%cvS(small_cells)
            call Sort(area_small_cells, indsort, .true.)
            small_cells = small_cells(indsort)

            ! Find merge face
            if (size(small_cells) /= 0) then

                do i = 1, size(small_cells)

                    ! Get faces
                    fcs = GetCellFaceGA(c, small_cells(i))

                    do j = 1, size(fcs)

                        cvs = GetFaceCellGA(c, fcs(j), cvLookUp)
                        if (size(cvs) == 2) then

                            if (cvs(1) /= small_cells(i)) then
                                neig = cvs(1)
                            else if (cvs(2) /= small_cells(i)) then
                                neig = cvs(2)
                            end if
                            
                            if (isPoloidal(grid, fcs(j), cvLookUp) &
                                .and. .not.any(fcs(j) == forbidden_fcs) &
                                .and. any(neig == small_cells)) then

                                    if (c%reg%Get(cvs(1)) == c%reg%Get(cvs(2))) then

                                        qm%merge_fc = fcs(j)

                                    end if

                            end if

                        end if
                        
                    end do

                    ! Exit when face is found
                    if (qm%merge_fc /= 0) exit

                end do

            end if
            
        case ('minimal_grid')

            call gdErrorHandler('SelectMergingFace: merge criterium minimal grid not implemented') ! TODO

        case ('h_pol')

            ! Remove merge cells with too small h_pol

            ! Distance function 
            call grid%fun%distr%Evaluate(c%x%Get(), c%y%Get(), dfunv)
            log = (dfunv .lt. options%dist_function_threshold_merge)
            indcv = (/ (i, i = 1, c%ntot)/)
            allocate(cells(count(log)))
            cells = pack(indcv, log)

            ! Second distance function use
            log = (dfunv .gt. options%dist_function_threshold_merge + 0.01_R8)
            indcv = (/ (i, i = 1, c%ntot)/)
            allocate(no_cells(count(log)))
            no_cells = pack(indcv, log)

            !  Take the 10 percentile biggest cell to avoid large Xcells
            ! influencing
            h_pol_no_cells_sorted = qm%h_pol(no_cells)
            allocate(indsort(size(h_pol_no_cells_sorted)))
            call Sort(h_pol_no_cells_sorted,indsort,.false.)
            deallocate(indsort)
            h_pol_no_cells_crit = h_pol_no_cells_sorted(nint(size(no_cells)/real(10, kind=R8)))

            ! Multiplication of criterium = criterium for good size
            ! All cells smaller should be merged
            crit = h_pol_no_cells_crit * options%merge_h_pol_factor

            ! Give triangles artifical larger h_pol
            trias_log = (c%vertP2%Get() ==3)
            allocate(trias(count(trias_log)))
            trias = pack(indcv, trias_log)
            qm%h_pol(trias) = qm%h_pol(trias)*2


            allocate(cells2(count(qm%h_pol(cells) .lt. crit)))
            cells2 = pack(cells, qm%h_pol(cells) .lt. crit)  
            h_pol_cells = qm%h_pol(cells2)
            allocate(indsort(size(h_pol_cells)))
            call Sort(h_pol_cells, indsort)    
            cells2 = cells2(indsort)

            ! Find the most interesting merge face
            if (size(cells2) /= 0) then

                do i = 1, size(cells2)

                    ! Get faces
                    fcs = GetCellFaceGA(c, cells2(i))

                    do j = 1, size(fcs)

                        cvs = GetFaceCellGA(c, fcs(j), cvLookUp)
                        if (size(cvs) == 2) then

                            if (cvs(1) /= cells2(i)) then
                                neig = cvs(1)
                            else if (cvs(2) /= cells2(i)) then
                                neig = cvs(2)
                            end if
                            
                            if (isPoloidal(grid, fcs(j), cvLookUp) &
                                .and. .not.any(fcs(j) == forbidden_fcs) &
                                .and. any(neig == cells2)) then

                                    if (c%reg%Get(cvs(1)) == c%reg%Get(cvs(2))) then

                                        qm%merge_fc = fcs(j)

                                    end if

                            end if

                        end if
                        
                    end do

                    ! Exit when face is found
                    if (qm%merge_fc /= 0) exit

                end do

            end if


        case ('bias')

            ! Merge the cells with strong bias
            ! Get non-aligned faces
            indfc = (/ (i, i = 1, f%ntot) /)
            log = (f%aligned%Get() == 0 .and. .not.isMember(indfc,forbidden_fcs))
            allocate(pol_faces(count(log)))
            pol_faces = pack(indfc, log)

            nf = size(pol_faces)
            allocate(bias(nf))
            bias = 0
            do i = 1, nf

                cvs = GetFaceCellGA(c, pol_faces(i), cvLookUp)
                if (size(cvs) == 2) then
                    if (c%reg%Get(cvs(1)) == c%reg%Get(cvs(2))) then

                        ! Get faces
                        fcs1 = GetCellFaceGA(c, cvs(1))
                        fcs2 = GetCellFaceGA(c, cvs(2))

                        if (sum(f%aligned%Get(fcs1)) == 2 .and. sum(f%aligned%Get(fcs2)) == 2) then
                            h_pol_cvs = qm%h_pol(cvs)
                            bias(i) = maxval(h_pol_cvs) / minval(h_pol_cvs)
                        end if

                    end if

                end if

            end do

            ! Sort for descending bias
            ! Get pol_faces with largest bias
            ind = maxloc(bias,1)
            if (bias(ind) .gt. options%merge_bias_limit) then
                qm%merge_fc = pol_faces(ind)
            end if

            !call Sort(bias, indsort, .false.)
            !pol_faces = pol_faces(indsort)
            !if (bias(1) .gt. options%merge_bias_limit) then
            !    qm%merge_fc = pol_faces(1)
            !end if

        case ('pol_flux')

            ! Distance function fun_r
            ! Only use cells in SOL
            log = (c%reg%Get() == 2)
            allocate(cellsD(count(log)))
            cellsD = pack(indcv, log)

            allocate(pol_fluxdens_estD(size(cellsD)))
            call grid%fun_r%distr%Evaluate(c%x%Get(cellsD), c%y%Get(cellsD), pol_fluxdens_estD)

            log = (pol_fluxdens_estD .lt. 0.5_R8)
            allocate(cells(count(log)))
            cells = pack(cellsD, log)

            allocate(pol_fluxdens_est(count(log)))
            pol_fluxdens_est = pack(pol_fluxdens_est, log)

            pol_fluxdens_est = pol_fluxdens_est*qm%h_rad(cells)
            mean_pol_flux = sum(pol_fluxdens_est)/size(pol_fluxdens_est)

            call Sort(pol_fluxdens_est, ind_sort)
            cells = cells(ind_sort)

            do i = 1, size(pol_fluxdens_est)

                if (pol_fluxdens_est(i) .lt. mean_pol_flux * 0.2_R8) then

                    neigs = GetCellNeigsGA(grid, cells(i), cvLookUp)
                    bench = 1000.0_R8
                    do j = 1, size(neigs)
                        if (any(neigs(j) == cells)) then
                            ind = findloc(cells, neigs(j), 1)
                            if (pol_fluxdens_est(ind) .lt. bench) then
                                common_face = GetCommonFace(c, cells(i), neigs(j))
                                if ((f%aligned%Get(common_face) == 1) .and.(f%label%Get(common_face) == 0)) then
                                    qm%merge_fc = common_face
                                    bench = pol_fluxdens_est(ind)
                                end if
                            end if
                        end if
                    end do

                    if (qm%merge_fc /= 0) exit

                end if

            end do

        case ('h_rad')

            ! Initialize
            ! Get all cells around the separatrices
            call grid%GiveSeparatrices(.true.,.true.,.false.,cvLookUp)
            counter = 0
            allocate(cvs_sep(c%ntot))
            do k = 1, grid%data%nsep
                vxs = GetFluxSurfaceVxsGA(fd, grid%data%sepID(k))
                cvs_sep = 0 
                do i = 1, size(vxs)
                    cvs = GetVertCellGA(c, vxs(i), cvLookUp)
                    do j = 1, size(cvs)
                        cvs_sep(counter + j) = cvs(j)
                    end do
                    counter = counter + size(cvs)
                end do
            end do

            call Unique(cvs_sep(1:counter), cvs_sepU) 

            ! Find cells with smallest psi width
            ! First only pick the cells smaller than h_rad_threshold and from creg 1 or 2
            cellsD = (/(i, i = 1, c%ntot)/)
            log = (qm%h_rad .lt. options%h_rad_threshold .and. c%reg%Get() .lt. 3)
            allocate(cells(count(log)))
            cells = pack(cellsD, log)

            ! Sort cells for smallest h_rad_psi
            h_rad_cells = qm%h_rad_psi(cells)
            call Sort(h_rad_cells, ind_sort)
            cells = cells(ind_sort)

            do i = 1, size(cells)

                if ((.not.any(cells(i) == cvs_sepU))) then
                    neigs = GetCellNeigsGA(grid, cells(i), cvLookUp)
                    bench = 1e10_R8

                    ! Loop over neigs
                    do j = 1, size(neigs)

                        ! If neig is also in cells, not next to separatrix and smaller than threshold
                        if (any(neigs(j) == cells) &
                            .and. (.not.any(neigs(j) == cvs_sepU)) &
                            .and. (qm%h_rad(neigs(j)) .lt. options%h_rad_threshold) ) then

                                if (qm%h_rad_psi(neigs(j)) .lt. bench) then

                                    common_face = GetCommonFace(c, cells(i), neigs(j))
                                    if ((f%aligned%Get(common_face) == 1) .and. (f%label%Get(common_face) == 0)) then

                                        qm%merge_fc = common_face
                                        bench = qm%h_rad_psi(neigs(j))

                                    end if
                                end if

                        end if
                    end do


                end if

                if (qm%merge_fc /= 0) exit

            end do

        case ('h_rad_core')

            ! Initialize
            call grid%GetCutsXpoints(cvLookUp, core_facesD)
            Xcells = GetVertCellGA(c, grid%data%xpointID(1), cvLookUp)
            allocate(core_faces(count(f%aligned%Get(core_facesD)==0)))
            core_faces = pack(core_facesD,f%aligned%Get(core_facesD)==0)

            allocate(cells(size(core_faces)*2))
            cells = 0
            do i = 1, size(core_faces)
                cvs = GetFaceCellGA(c, core_faces(i), cvLookUp)
                cells(2*i-1) = cvs(1)
                cells(2*i) = cvs(2)
            end do

            ! Sort to smallest h_rad
            h_rad_cells = qm%h_rad(cells)
            call Sort(h_rad_cells, ind_sort)
            cells = cells(ind_sort)
            
            do i = 1, size(cells)
                if (.not.any(cells(i)==Xcells) .and. (qm%h_rad(cells(i)) .lt. options%h_rad_core_threshold)) then

                    neigs = GetCellNeigsGA(grid, cells(i), cvLookUp)
                    bench = 1e10_R8
                    do j = 1, size(neigs)
                        if (any(neigs(j) == cells) .and. (.not.any(neigs(j) == Xcells)) &
                            .and. (qm%h_rad(neigs(j)) .lt. options%h_rad_core_threshold)) then
                            if (qm%h_rad(neigs(j)) .lt. bench) then
                                common_face = GetCommonFace(c, cells(i), neigs(j))
                                if (f%aligned%Get(common_face) == 1 .and. f%label%Get(common_face)== 0) then
                                    qm%merge_fc = common_face
                                    bench = qm%h_rad(neigs(j))
                                end if
                            end if
                        end if
                    end do

                end if

                if (qm%merge_fc /= 0) exit

            end do




        case ('bias_rad_farSOL')

            ! Remove cells with strong bias in the radial direction in the farSOL region - TODO after detect farSOL


        case ('bias_rad')

            ! Remove cells with strong bias in the radial direction

            ! Get aligned faces
            indfc = (/(i, i = 1, f%ntot)/)
            allocate(rad_facesD(count(f%aligned%Get() == 1)))
            rad_facesD = pack(indfc, f%aligned%Get() == 1)

            log = (.not.isMember(rad_facesD, forbidden_fcs))
            allocate(rad_faces(count(log)))
            rad_faces = pack(rad_facesD, log)

            ! Compute bias for these faces
            allocate(bias(size(rad_faces)))
            bias = 0
            do i = 1, size(rad_faces)

                cvs = GetFaceCellGA(c, rad_faces(i), cvLookUp)
                if (size(cvs) == 2) then
                    if (c%reg%Get(cvs(1)) == c%reg%Get(cvs(2))) then

                        ! Both only quads
                        fcs1 = GetCellFaceGA(c, cvs(1))
                        fcs2 = GetCellFaceGA(c, cvs(2))

                        if ((count(f%aligned%Get(fcs1) == 1) .gt. 1) .and. (count(f%aligned%Get(fcs2) == 1) .gt. 1)) then

                            bias(i) = maxval(qm%h_rad(cvs)) / minval(qm%h_rad(cvs))
    
                        end if

                    end if

                end if

            end do

            ! Pick the largest bias
            call Sort(bias, indsort, .false.)
            rad_faces = rad_faces(indsort)
            if (bias(1) .gt. options%merge_bias_limit) then
                qm%merge_fc = rad_faces(1)
            end if


        case ('skew_tria')

            ! TODO 

        case ('manual')

            call gdErrorHandler('SelectMergingFace: manual merging via input not possible in precompile code')
        
        case default

            call gdErrorHandler('SelectMergingFace: merge criterium not implemented')

        end select

        end associate

    end subroutine

    subroutine SelectSplitCell(qm, grid, options)

        ! Description
        !============
        ! Select a cell to split based on the criteria chosen in options

        ! Declare variables
        !==================
        ! Arguments
        class(QualityMetricUDT), intent(inout)  :: qm
        type(GAGridUDT), intent(inout)          :: grid
        type(GAoptionsUDT), intent(in)          :: options

        ! Auxiliary
        integer(I8) :: i, k, ind, nal1, nal2, fcs1_al, fcs2_al, indmin, &
            splitfaces(2), n_sub
        integer(I8), allocatable, dimension(:) :: indcv, cells,  &
            cellsD, indsort, cells2, fcs, fcs_sep, indf, b_faces, b_verts, &
            cvLookUp, cvs, cvsD, fcs1, fcs2, tria_cells, subset_tria
        real(R8) :: mean_tot_flux, threshold, dfunv(grid%cell%ntot)
        real(R8), allocatable, dimension(:) :: h_rad_cells, pol_fluxdens_est, &
            dpsi, cvAR_tria, cvS_tria, h_pol_sol_cells
        logical, allocatable :: log(:), log2(:)

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert, &
            fd => grid%data%fluxdata &
            )

        ! Initialize
        qm%split_cv = 0
        indcv = (/ (i, i = 1, c%ntot )/)

        select case (options%splittype)

        case ('rad')

            ! Apply distance function or area constraint
            if (options%dist_function) then

                call grid%fun%distr%Evaluate(c%x%Get(), c%y%Get(), dfunv)
                log = (dfunv .lt. options%dist_function_threshold_split)
                allocate(cellsD(count(log)))
                cellsD = pack(indcv, log)

            else if (options%no_pents_area_split) then

                call grid%AreaConstraintPents(options, options%dist_function_threshold_split, cellsD)

            end if

            ! Only select quads to start radial splitting
            allocate(cells(count(c%faceP2%Get(cellsD) == 4)))
            cells = pack(cellsD, c%faceP2%Get(cellsD) == 4)

            ! Select split criterium
            select case (options%rad_type)

            case ('h_rad_psi')

                ! Sort for decreasing h_rad_psi
                h_rad_cells = qm%h_rad_psi(cells)
                call Sort(h_rad_cells, indsort, .false.)
                cells = cells(indsort)

                ! Find split cell
                i = 1
                do while (qm%split_cv == 0)

                    if (c%reg%Get(cells(i)) == 2 .and. c%cflags%Get(cells(i)) /= 3) then
                        qm%split_cv = cells(i)
                    else 
                        i = i + 1
                    end if

                end do
 
            case ('h_rad')

                ! Sort for decreasing h_rad_psi
                h_rad_cells = qm%h_rad(cells)
                call Sort(h_rad_cells, indsort, .false.)
                cells = cells(indsort)

                ! Find split cell
                i = 1
                do while (qm%split_cv == 0)

                    if (c%reg%Get(cells(i)) /= 1 .and. c%cflags%Get(cells(i)) /= 3) then
                        qm%split_cv = cells(i)
                    else 
                        i = i + 1
                    end if

                end do

            case ('pol_flux')

                ! Distance function fun_r
                allocate(cells2(count(c%reg%Get(cells) == 2)))
                cells2 = pack(cells, c%reg%Get(cells) == 2)
                call grid%fun_r%distr%Evaluate(c%x%Get(cells2),c%y%Get(cells2), pol_fluxdens_est)
                pol_fluxdens_est = pol_fluxdens_est * qm%h_rad(cells2)
                mean_tot_flux = sum(pol_fluxdens_est) / size(pol_fluxdens_est)

                ! Sort for descending absolute poloidal flux
                call Sort(pol_fluxdens_est, indsort, .false.)
                cells2 = cells2(indsort)

                ! Get separatrix
                cvLookUp = GetCvLookUpGA(c)
                call grid%GiveSeparatrices(.true., .true., .false., cvLookUp)
                if (grid%data%nsep == 1) then
                    fcs_sep = GetFluxSurfaceFcsGA(fd, grid%data%sepID(1))
                else
                    call gdErrorHandler('SelectSplitCell: pol_flux criterion, multiple separatrix not implemented yet')
                end if

                do i = 1, size(fcs_sep)
                    if (pol_fluxdens_est(i) .gt. mean_tot_flux*2) then

                        fcs = GetCellFaceGA(c, cells2(i))
                        if (.not.any(fcs == fcs_sep)) then
                            qm%split_cv = cells2(i)
                            exit
                        end if

                    end if

                end do



            case ('farSOL')

                ! TODO

            case ('farSOLrefinement')

                ! TODO

            case ('farSOLrefinement_targets')

                ! TODO

            case ('no_aligned_faces')

                ! Loop over cells
                do i = 1, c%ntot

                    fcs = GetCellFaceGA(c, i)

                    if (.not.any(f%aligned%Get(fcs) == 1) .and. size(fcs) == 4) then

                        qm%split_cv = i
                        exit

                    end if
                
                end do

            case ('shaved_off_tubes')

                ! Automatic option, used before splitting and merging
                ! Get boundary vertices
                indf = (/ (i, i = 1, f%ntot)/)
                log = isBoundaryFaceGA(f, indf)
                allocate(b_faces(count(log)))
                b_faces = pack(indf, log)
                b_verts = GetVxsFromFcsGA(f, b_faces)

                cvLookUp = GetCvLookUpGA(c)
                threshold = 0.1_R8
                do i = 1, size(b_verts)

                    ! Get boundary cells of vertex
                    cvsD = GetVertCellGA(c, b_verts(i), cvLookUp)
                    log2 = (isBoundaryCellGA(grid, cvsD))
                    allocate(cvs(count(log2)))
                    cvs = pack(cvsD, log2)

                    if (size(cvs) == 2) then ! Boundary vertices has two boundary cells

                        ! Check if both cells are trapezoids
                        fcs1 = GetCellFaceGA(c, cvs(1))
                        fcs2 = GetCellFaceGA(c, cvs(2))

                        nal1 = count(f%aligned%Get(fcs1) == 1)
                        nal2 = count(f%aligned%Get(fcs2) == 1)

                        ! Check triangles and quad combinations - futher
                        ! checks not required if both triangles
                        if (size(fcs1) == 4 .and. size(fcs2) == 4) then

                            ! Determine specific case
                            if (nal1 == 1 .and. nal2 == 1) then

                                ! Two trapezoids - check the trapezoids are neighbors with a common face
                                if (HaveCommonFace(c, cvs(1), cvs(2))) then

                                    ! Look at the psi values of cells and the vert
                                    ! Get vertices of boundary faces to check
                                    ! the inclination via the psi values
                                    if (CheckCellPsiIntersection(grid,b_verts(i),cvs(1)) &
                                        .and. CheckCellPsiIntersection(grid,b_verts(i),cvs(2)) ) then

                                        ! Check whether position of intersection is not to close to another vertex
                                        if (CompareCellPsiIntersection(grid, cvs(1), fcs1, b_verts(i), threshold)) then
                                            qm%split_cv = cvs(1)
                                            exit
                                        else if (CompareCellPsiIntersection(grid, cvs(2), fcs2, b_verts(i), threshold)) then
                                            qm%split_cv = cvs(2)
                                            exit
                                        end if

                                    end if

                                end if

                            else if (nal1 == 1 .and. nal2 /= 1) then
                                ! One trapezoids and one aligned quad, so
                                ! split the trapezoids if concave 
                                ! Check for concave boundary
                                ! configuration
                                if (CheckCellPsiIntersection(grid,b_verts(i), cvs(1))) then
                                    if (CompareCellPsiIntersection(grid, cvs(1), fcs1, b_verts(i), threshold)) then
                                        qm%split_cv = cvs(1)
                                        exit
                                    end if
                                end if

                            else if (nal1 /= 1 .and. nal2 == 1) then
                                ! One trapezoids and one aligned quad, so
                                ! split the trapezoids
                                if (CheckCellPsiIntersection(grid,b_verts(i), cvs(2))) then
                                    if (CompareCellPsiIntersection(grid, cvs(2), fcs2, b_verts(i), threshold)) then
                                        qm%split_cv = cvs(2)
                                        exit
                                    end if
                                end if

                            end if

                        else if (size(fcs1) == 3 .and. size(fcs2) == 4) then
                            ! Triangle and quad
                            if (nal2 == 1) then
                                ! Check if the aligned faces of both cells
                                ! are on the same flux surface
                                ind = findloc(f%aligned%Get(fcs1), 1, 1)
                                fcs1_al = fcs1(ind)
                                ind = findloc(f%aligned%Get(fcs2), 1, 1)
                                fcs2_al = fcs2(ind)

                                if (.not.HaveCommonVert(f, fcs1_al, fcs2_al) &
                                    .and. CheckCellPsiIntersection(grid, b_verts(i), cvs(2))) then

                                        if (CompareCellPsiIntersection(grid, cvs(2), fcs2, b_verts(i), threshold)) then
                                            qm%split_cv = cvs(2)
                                            exit
                                        end if

                                end if

                            end if

                        else if (size(fcs1) == 4 .and. size(fcs2) == 3) then

                            if (nal1 == 1) then ! Trapezoid
                                ! Check if the aligned faces of both cells
                                ! are on the same flux surface
                                ind = findloc(f%aligned%Get(fcs1), 1, 1)
                                fcs1_al = fcs1(ind)
                                ind = findloc(f%aligned%Get(fcs2), 1, 1)
                                fcs2_al = fcs2(ind)

                                if (.not.HaveCommonVert(f, fcs1_al, fcs2_al) &
                                    .and. CheckCellPsiIntersection(grid, b_verts(i), cvs(1))) then

                                        if (CompareCellPsiIntersection(grid, cvs(1), fcs1, b_verts(i), threshold)) then
                                            qm%split_cv = cvs(1)
                                            exit
                                        end if

                                end if
                            end if

                        end if

                    end if

                end do


            case ('manual')

                call gdErrorHandler('SelectSplitCell: manual selection method is not supported in pre-compiled framework')

            case default

                call gdErrorHandler('SelectSplitCell: split criterion not implemented')

            end select


        case ('pol')

            select case (options%pol_type)

            case ('trias_cvS')

                ! Get triangles with largest area
                log = (c%faceP2%Get(indcv) == 3)
                allocate(tria_cells(count(log)))
                tria_cells = pack(indcv, log)

                ! Sort for descending surface area
                ! Get triangle with largest area
                ind = maxloc(qm%cvS(tria_cells),1)
                qm%split_cv = tria_cells(ind)

            case ('h_pol')
                
                ! Restrict to area where pents are allowed
                if (options%no_pents) then
                    call grid%AreaConstraintPents(options, options%dist_function_threshold_split, cells)
                end if

                ! Criterion
                h_pol_sol_cells = qm%h_pol(cells)
                allocate(indsort(size(cells)))
                call Sort(h_pol_sol_cells, indsort, .false.)
                cells = cells(indsort)
                i = 1
                cvLookUp = GetCvLookUpGA(c)
                do while (qm%split_cv == 0) 

                    qm%split_cv = cells(i)

                    if (c%faceP2%Get(qm%split_cv) == 4) then

                        ! Get faces and check whether the cell has two split faces
                        fcs = GetCellFaceGA(c, qm%split_cv)

                        dpsi = abs(v%psi%Get(f%vert1%Get(fcs)) - v%psi%Get(f%vert2%Get(fcs)))
                        indmin = minloc(dpsi,1)
                        splitfaces(1) = fcs(indmin)
                        
                        dpsi(indmin) = 1e99_R8
                        indmin = minloc(dpsi,1)
                        splitfaces(2) = fcs(indmin)
                        
                        do k = 1, 2
                            cvs = GetFaceCellGA(c, splitfaces(k), cvLookUp)
                            if (any(c%faceP2%Get(cvs) == 5)) then
                                qm%split_cv = 0
                                i = i + 1
                                exit
                            end if
                        end do

                    else if (c%faceP2%Get(qm%split_cv) == 3 .or. c%cflags%Get(qm%split_cv) == 3) then

                        qm%split_cv = 0
                        i = i + 1

                    end if

                end do

            case ('trias_farSOL')

                ! Get triangular cells
                log = (c%faceP2%Get() ==3)
                allocate(tria_cells(count(log)))
                tria_cells = pack(indcv, log)
                
                ! Get largest triangles by area
                cvS_tria = qm%cvS(tria_cells)
                n_sub = max(nint(size(tria_cells)/1.5_R8), 1)

                allocate(indsort(size(tria_cells)))
                call Sort(cvS_tria, indsort, .false.)
                tria_cells = tria_cells(indsort)
                deallocate(indsort)
                subset_tria = tria_cells(1:n_sub)

                ! Get triangles with the highest AR
                cvAR_tria = qm%cvAR(subset_tria)
                allocate(indsort(size(subset_tria)))
                call Sort(cvAR_tria, indsort, .false.)
                subset_tria = subset_tria(indsort)

                do i = 1, n_sub
                    if (cvAR_tria(i) .gt. 3) then
                        qm%split_cv = subset_tria(i)
                        exit
                    end if
                end do



            case ('farSOLrefinement_hpol')

                ! TODO

            case ('farSOLrefinement_targets')

                ! TODO

            case ('manual')

                call gdErrorHandler('SelectSplitCell: manual selection method is not supported in pre-compiled framework')

            case default

                call gdErrorHandler('SelectSplitCell: split criterion not implemented')

            end select

        end select

        end associate

    end subroutine

    function CheckCellPsiIntersection(grid, iv, ic) result(res)
        type(GAGridUDT) :: grid
        integer(I8) :: iv, ic
        logical :: res

        integer(I8), allocatable :: vxs(:)
        real(R8) :: psi_iv
        real(R8), allocatable :: psi_cv(:)

        res = .false.
        vxs = GetCellVertGA(grid%cell, iv)
        psi_cv = grid%vert%psi%Get(vxs)
        psi_iv = grid%vert%psi%Get(iv)
        if (psi_iv .lt. maxval(psi_cv) .and. psi_iv .gt. minval(psi_cv)) then
            res = .true.
        end if

    end function

    function CompareCellPsiIntersection(grid, ic, fcs, iv, threshold) result(res)
        type(GAGridUDT) :: grid
        integer(I8) :: iv, ic
        integer(I8), allocatable :: fcs(:), int_faces(:)
        real(R8) :: threshold, psic, min_dpsi, dpsi_v(2)
        real(R8), allocatable :: psi_f1(:), psi_f2(:), dpsi_f(:)
        logical :: res 
        
        res = .false.
        psic = grid%vert%psi%Get(iv)
        call grid%GetIntersectedPsiFaces(ic, fcs, int_faces, psic)
        if (size(int_faces) .gt. 1) then
            call gdErrorHandler('CompareCellPsiIntersection: more sections than expected (one expected)')
        else
            psi_f1 = grid%vert%psi%Get(grid%face%vert1%Get(int_faces))
            psi_f2 = grid%vert%psi%Get(grid%face%vert2%Get(int_faces))
            dpsi_f = abs(psi_f1(1) - psi_f2(1))
            dpsi_v(1) = abs(psi_f1(1) - psic)
            dpsi_v(2) = abs(psi_f2(1) - psic)
            min_dpsi = minval(dpsi_v)

            if ((min_dpsi/dpsi_f(1)) .gt. threshold) then
                res = .true.
            end if
        end if
        
    end function


    !------------------------------------------------------------------!
    !                        GAGRID ROUTINES                           !
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
            deallocate(GAvert%x, GAvert%y, &
                GAvert%psi, GAvert%bx, GAvert%by)
        end if 
        allocate(RealDynamicArrayBufferedUDT::GAvert%x, GAvert%y, &
            GAvert%psi, GAvert%bx, GAvert%by)


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
            deallocate(GAface%vert1, GAface%vert2, GAface%label,  &
                GAface%reg, GAface%aligned)
        end if 
        allocate(IntegerDynamicArrayBufferedUDT:: GAface%vert1, GAface%vert2, GAface%label, &
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
                GAcell%x, GAcell%y, GAcell%cflags, GAcell%reg)
        end if 
        allocate(RealDynamicArrayBufferedUDT:: GAcell%psi, GAcell%x, GAcell%y )
        allocate(IntegerDynamicArrayBufferedUDT:: GAcell%vertP1, GAcell%vertP2, &
                 GAcell%vert, GAcell%faceP1, GAcell%faceP2, GAcell%face, &
                 GAcell%cflags, GAcell%reg)

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
                GAfluxdata%fluxsurfacefaces,  GAfluxdata%fluxsurfacevertsP1, &
                GAfluxdata%fluxsurfacevertsP2, GAfluxdata%fluxsurfaceverts)
        end if 
        allocate(IntegerDynamicArrayBufferedUDT:: GAfluxdata%fluxsurfacefacesP1, &
                GAfluxdata%fluxsurfacefacesP2, GAfluxdata%fluxsurfacefaces, &
                GAfluxdata%fluxsurfacevertsP1, GAfluxdata%fluxsurfacevertsP2, &
                GAfluxdata%fluxsurfaceverts)


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
        GAv%x           = ConstructRealDynamicArrayBuffered(gv%x)
        GAv%y           = ConstructRealDynamicArrayBuffered(gv%y)
        GAv%bx          = ConstructRealDynamicArrayBuffered(gv%bx)
        GAv%by          = ConstructRealDynamicArrayBuffered(gv%by)
        GAv%psi         = ConstructRealDynamicArrayBuffered(gv%psi)
        GAv%ffbz        = gv%ffbz(1)
        GAv%ntot        = gv%ntot

        ! Face information
        GAf%vert1   = ConstructIntegerDynamicArrayBuffered(gf%vert(:,1))
        GAf%vert2   = ConstructIntegerDynamicArrayBuffered(gf%vert(:,2))
        GAf%label   = ConstructIntegerDynamicArrayBuffered(gf%label)
        GAf%reg     = ConstructIntegerDynamicArrayBuffered(gf%reg)
        GAf%aligned = ConstructIntegerDynamicArrayBuffered(gf%aligned)
        GAf%ntot    = gf%ntot

        ! Cell information
        GAc%vertP1  = ConstructIntegerDynamicArrayBuffered(gc%vertP(:,1))
        GAc%vertP2  = ConstructIntegerDynamicArrayBuffered(gc%vertP(:,2))
        GAc%vert    = ConstructIntegerDynamicArrayBuffered(gc%vert)
        GAc%faceP1  = ConstructIntegerDynamicArrayBuffered(gc%faceP(:,1))
        GAc%faceP2  = ConstructIntegerDynamicArrayBuffered(gc%faceP(:,2))
        GAc%face    = ConstructIntegerDynamicArrayBuffered(gc%face)
        GAc%cflags  = ConstructIntegerDynamicArrayBuffered(gc%cflags)
        GAc%reg     = ConstructIntegerDynamicArrayBuffered(gc%reg)
        GAc%psi     = ConstructRealDynamicArrayBuffered(gc%psi)
        GAc%x       = ConstructRealDynamicArrayBuffered(gc%x)
        GAc%y       = ConstructRealDynamicArrayBuffered(gc%y)
        GAc%ntot    = gc%ntot
        GAc%ngc     = gc%ngc
        GAc%nvert   = gc%nvert
        GAc%nface   = gc%nface
        

        ! Grid data - flux surface data
        GAgrid%data%xpointID    = grid%data%xpointID
        GAgrid%data%nxp         = grid%data%nxp
        GAgrid%data%sepID       = grid%data%sepID
        GAgrid%data%nsep        = grid%data%nsep
        GAfd%fluxsurfacefacesP1 = ConstructIntegerDynamicArrayBuffered(gfd%fluxsurfacefacesP(:,1))
        GAfd%fluxsurfacefacesP2 = ConstructIntegerDynamicArrayBuffered(gfd%fluxsurfacefacesP(:,2))
        GAfd%fluxsurfacefaces   = ConstructIntegerDynamicArrayBuffered(gfd%fluxsurfacefaces)
        GAfd%fluxsurfacevertsP1 = ConstructIntegerDynamicArrayBuffered(gfd%fluxsurfacevertsP(:,1))
        GAfd%fluxsurfacevertsP2 = ConstructIntegerDynamicArrayBuffered(gfd%fluxsurfacevertsP(:,2))
        GAfd%fluxsurfaceverts   = ConstructIntegerDynamicArrayBuffered(gfd%fluxsurfaceverts)
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
        !integer(I8) :: i, indFc(GAgrid%face%ntot), fcLbl_loc(GAgrid%face%ntot), nfb
        !integer(I8), allocatable :: fcsbnd(:), vxsbnd(:)


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
        gv%ffbz         = GAv%ffbz

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
        gc%reg          = Gac%reg%GetAllElements()
        gc%psi          = GAc%psi%GetAllElements()
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

        ! Problem need to compute some extra fields for grid 
        ! See what is needed for WriteGOAT: TODO
        ! To give to grid deformation
        ! => fd%fluxsurfaceID, vert%fieldlineID



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
        logical, intent(out)        :: is_ordered(grid%cell%ntot)
        logical, intent(in)         :: cells(grid%cell%ntot)

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

        ! Loop over the cells
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
                else
                    print *, 'Error on cell: ', ic, ' ', c%cflags%Get(ic), ' ', c%faceP2%Get(ic)
                    print *, 'cx: ', c%x%Get(ic), ' cy: ', c%y%Get(ic)
                    print *, 'vx :', v%x%Get(tv)
                    print *, 'vy :', v%y%Get(tv)
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

    subroutine GetFsVxFromFsFc(grid, options)
        ! Description
        !============
        ! Get fsVx from fsFc

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout) :: grid
        type(GAoptionsUDT) :: options

        ! Auxiliary
        integer(I8) :: nv_counter, ifs, nv
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

        if (options%debug) call CheckUniqueness(fd%fluxsurfacefaces)
        
        do ifs = 1, fd%nFs
            ! Get vertices from flux surface
            fcs = GetFluxSurfaceFcsGA(fd, ifs)
            vxs = GetVxsFromFcsGA(f, fcs)
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

        ! Check uniqueness
        if (options%debug) call CheckUniqueness(fd%fluxsurfaceverts)

        end associate


    end subroutine

    subroutine GiveXpoints(grid,use_sep,cvLookUp)

        ! Description
        !============
        ! Gives the Xpoint(s). Depending the grid information different methods are used. If all vertices have a fieldlineID the routine DetermineXPointsGA is used, otherwise the xpoint is determined by checking whether a vertex is surrounded by more then three regions and more than four cells. 

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)     :: grid
        logical, intent(in)                 :: use_sep  
        integer(I8), allocatable, optional  :: cvLookUp(:)            

        ! Auxiliary
        integer(I8)                         :: i, j, iv, xpoints(1:100), counter, & 
            n , nc
        integer(I8), allocatable            :: vxs(:), cells(:), &
            regions(:), creg(:)
        logical :: use_sepID, start, use_nsep
        

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
                call grid%GiveSeparatrices(use_nsep,use_sepID,start,cvLookUp)
            end if

            ! Determine whether to use the vert.fieldlineID to determine the Xpoint
            !use_fieldlineID = .false.
            !if (allocated(v%fieldlineID)) then
            !    if (v%fieldlineID%Size().eq.v%ntot) then
            !        use_fieldlineID = .true.
            !    end if
            !end if

            !if (use_fieldlineID) then

                !!! NOT USED BECAUSE ALGO USES vert%neig WHICH IS NOT PRESENT IN GAGRID !!
            !    call DetermineXPointsGA(xpind, nxpind, order, grid)
            !    grid%data%xpointID  = xpind
            !    grid%data%nxp       = nxpind 

            if ((allocated(fd%fluxsurfaceverts))) then!.and.(allocated(v%cellP1))) then

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
        class(GAGridUDT), intent(inout)     :: grid
        logical, intent(in)                 :: use_nsep, use_sepID, start
        integer(I8), allocatable, optional  :: cvLookUp(:)

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
                            nf = fd%fluxsurfacefacesP2%Get(grid%data%sepID(1))
                            fcs = GetFluxSurfaceFcsGA(fd, grid%data%sepID(1))
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
        class(GAGridUDT), intent(inout)      :: grid
        type(GAoptionsUDT), intent(in)       :: options
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
                call Unique(grid%data%fluxdata%fluxsurfacefaces%Get(), fcs)
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
        class(GAGridUDT), intent(in)    :: grid
        logical                         :: check_extra_conn

        ! Auxiliary
        integer(I8), allocatable, dimension(:) :: cvLookUp, v1n, v2n, cf, &
            cellnumbers, verts_of_cell, cvertP2, cfaceP2, fcs, ccflags, nvxs, &
            b_cells, indCv, cvertP1, cfaceP1
        integer(I8) :: verts_of_face(1:2), nf, l, ic, nb, i, j, counter, &
            nface, nvert, ncpf(grid%face%ntot)
        logical, allocatable :: logcf(:)   

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
            cvertP1 = c%vertP1%GetAllElements()
            cvertP2 = c%vertP2%GetAllElements()
            cfaceP1 = c%faceP1%GetAllElements()
            cfaceP2 = c%faceP2%GetAllElements()
            ccflags = c%cflags%GetAllElements()

            ! Check 0 - Check pointer consistency
            if (any(cvertP1(2:c%ntot) /= (cvertP1(1:c%ntot-1)+cvertP2(1:c%ntot-1)))) then

                ! Get where the inconsistency is 
                do ic = 2, c%ntot
                    if (cvertP1(ic) /= cvertP1(ic-1)+cvertP2(ic-1)) then
                        print *, 'Cell:', ic 
                        call gdErrorHandler('CheckUnstructuredGrid: check 0, point cvertP not consistent')
                    endif
                end do 

            end if 
            if (any(cfaceP1(2:c%ntot) /= cfaceP1(1:c%ntot-1)+cfaceP2(1:c%ntot-1))) then

                ! Get where the inconsistency is 
                do ic = 2, c%ntot
                    if (cfaceP1(ic) /= cfaceP1(ic-1)+cfaceP2(ic-1)) then
                        print *, 'Cell:', ic 
                        call gdErrorHandler('CheckUnstructuredGrid: check 0, point cfaceP not consistent')
                    endif                    
                end do

            end if

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

            ! Initialize
            nvert = c%vertP1%Get(c%ntot)+c%vertP2%Get(c%ntot)-1
            allocate(cvLookUp(nvert))
            cvLookUp = GetCvLookUpGA(c)
            allocate(cellnumbers(10))
            cellnumbers = 0
            allocate(logcf(size(cf))) 

            do j = 1,f%ntot
                verts_of_face(1) = v1n(j)
                verts_of_face(2) = v2n(j)

                ! Find a cell where the face is attached
                logcf = (cf == j)
                ncpf(j) = count(logcf)
                !allocate(cellnumbers(nc))
                cellnumbers(1:ncpf(j)) = pack(cvLookUp,logcf)
                if (ncpf(j) .gt. 2) then
                    print *, 'Face: ', j
                    call gdErrorHandler('CheckUnstructuredGrid: check 3, more than two' // &
                     & 'cells connected to face')
                end if

                !get vertices of cells
                do l = 1, ncpf(j)
                    verts_of_cell = GetCellVertGA(c,cellnumbers(l))
                    !verts_of_face moeten voorkomen in verts_of_cells
                    do i = 1,2
                        if (.not.any(verts_of_cell == verts_of_face(i))) then
                                print *, 'Face: ', j
                                print *, 'V1 :', verts_of_face(1)
                                print *, 'x coor: ', v%x%Get(verts_of_face(1))
                                print *, 'y coor: ', v%y%Get(verts_of_face(1))
                                print *, 'V2 :', verts_of_face(2)
                                print *, 'x coor: ', v%x%Get(verts_of_face(2))
                                print *, 'y coor: ', v%y%Get(verts_of_face(2))    
                                call gdErrorHandler('CheckUnstructuredGrid: check 3, vertices' // &
                                & ' of face not in vertices of the cells of that face')
                        end if
                    end do
                end do

            end do

            ! Check 4 - whether all internal cells have faces which are connected to two cells
            do ic = 1, c%ntot

                ! Only internal faces
                if (ccflags(ic) == 1) then
                    fcs = GetCellFaceGA(c, ic)
                    do i = 1, size(fcs)
                        if (.not.( ncpf(fcs(i)) .eq. 2 ) ) then
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
                    call Unique(fd%fluxsurfaceverts%Get(), nvxs)
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
            do i = 1, nb

                nf = cfaceP2(b_cells(i))
                fcs = GetCellFaceGA(c, b_cells(i))

                counter = 0
                do j = 1,nf
                    if ( ncpf(fcs(j)) .eq. 1 )  then
                        counter = counter + 1
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

            ! Extra check on the extra connectivity fields - TODO
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
        ! mergeing.   

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
        allocate(fieldlineID(v%ntot))
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
                allocate(ar(count(hasID(tfv) == 1)))
                ar = pack(fID,hasID(tfv) == 1)
                ifs1 = ar(1)
                if (size(ar) .gt. 1) &
                    print *, 'Warning: MergeFS: merging of multiple flux surfaces not yet supported'

                counter = counter + 1
                mFS(counter,1) = ifs1
                mFS(counter,2) = ifs

                deallocate(ar)

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
            call fd%fluxsurfacevertsP1%SumMask(range, -nv1)
            call fd%fluxsurfacefacesP1%SumMask(range, -nf1)
            !call fd%fluxsurfacevertsP1%SetMultipleElements(range, fd%fluxsurfacevertsP1%GetMultipleElements(range) - nv1)
            !call fd%fluxsurfacefacesP1%SetMultipleElements(range, fd%fluxsurfacefacesP1%GetMultipleElements(range) - nf1)

            ifs2 = ifs2 - 1
            call fd%fluxsurfacevertsP1%RemoveSingleElement(ifs2)
            call fd%fluxsurfacevertsP2%RemoveSingleElement(ifs2)
            call fd%fluxsurfacefacesP1%RemoveSingleElement(ifs2)
            call fd%fluxsurfacefacesP2%RemoveSingleElement(ifs2)
            fd%nFs = fd%nFs - 1

            range = (/ (j, j = ifs2, fd%nFs)/)
            call fd%fluxsurfacevertsP1%SumMask(range, -nv2)
            call fd%fluxsurfacefacesP1%SumMask(range, -nf2)
            !call fd%fluxsurfacevertsP1%SetMultipleElements(range,fd%fluxsurfacevertsP1%GetMultipleElements(range) - nv2)
            !call fd%fluxsurfacefacesP1%SetMultipleElements(range,fd%fluxsurfacefacesP1%GetMultipleElements(range) - nf2)

            ! Add new flux surface
            fd%nFs = fd%nFs + 1
            call Unique([tv1, tv2], new_verts)
            call Unique([tf1, tf2], new_faces)

            call fd%fluxsurfacevertsP1%Append(fd%fluxsurfacevertsP1%Get(fd%nFs-1)+ &
                fd%fluxsurfacevertsP2%Get(fd%nFs-1))
            call fd%fluxsurfacevertsP2%Append(size(new_verts))
            call fd%fluxsurfaceverts%Append(new_verts)
            call fd%fluxsurfacefacesP1%Append(fd%fluxsurfacefacesP1%Get(fd%nFs-1)+ &
                fd%fluxsurfacefacesP2%Get(fd%nFs-1))
            call fd%fluxsurfacefacesP2%Append(size(new_faces))
            call fd%fluxsurfacefaces%Append(new_faces)

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
    !                       GAGRID ADAPTATIONS                         !
    !------------------------------------------------------------------!
    
    ! Removing small triangles
    !=========================
    subroutine RemoveSmallTriangle(grid, magneticField, qm, options)

        ! Description
        !============
        ! Removing of very small triangles at the boundary of the grid. 
        ! There are two methods: 
        ! 1) The merging method
        !       This method is used for small triangles and merges a whole
        !       row of cells to eliminate the small triangle.
        ! 2) Local method
        !       Eliminates mini triangles by only changing the cells in 
        !       the vicinty of the mini triangle.



        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)         :: grid
        type(QualityMetricUDT), intent(inout)   :: qm
        type(GAoptionsUDT), intent(in)          :: options
        type(MagneticFieldUDT)                  :: magneticField

        ! Auxiliary
        integer(I8) :: small_tria, faceA, faceB, faceC
        logical :: plain_merging
        type(GAoptionsUDT) :: options2

        ! Show progress
        print *, 'Remove small triangles'
    
        ! Recalculate areas and poloidal distances ! should be done!!!

        ! Detect small triangles
        !=======================
        small_tria = 0
        call grid%DetectSmallTrias(qm, options, small_tria)


        ! Remove the triangles
        !=====================
        do while (small_tria /= 0)

            ! Determine removal method
            call grid%DetermineRemovalMethod(qm, small_tria, faceA, faceB, faceC, plain_merging)

            if (plain_merging) then

                ! Method 1:
                print *, 'Apply cell merging'
                qm%merge_fc = faceC
                options2 = options
                options2%n_merge = 1
                options2%merging = .true.
                
                call grid%MergeCells(qm, options2)
                if (.not.options%slab) call grid%RecalcMagn(magneticField)

            else 

                ! Method2 2: Local method
                print *, 'Apply local removal method'

                ! Visualize - TODO

                call grid%LocalSmallTriangleRemoval(small_tria, faceA, faceB, faceC, options)

                if (.not.options%slab) call grid%RecalcMagn(magneticField)

            end if

            if (options%debug) call grid%CheckUnstructuredGrid(.false.)

            ! Redection
            small_tria = 0            
            call qm%ComputeQM(grid, options, magneticField)
            call grid%DetectSmallTrias(qm, options, small_tria)

        end do

        ! Show progress
        print *, 'Ended removing small triangles'



    end subroutine

    subroutine DetectSmallTrias(grid, qm, options, small_tria)

        ! Description
        !============
        ! Detects small triangles according to criteria, cut_off_pol
        ! and cut_off_surf. One triangle is select at a time because
        ! cell numbers change when doing adaptation, so detecting
        ! all the small triangles at the same time would require
        ! elaborate operations on the small_trias array

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(in)        :: grid
        type(QualityMetricUDT), intent(in)  :: qm
        type(GAoptionsUDT), intent(in)      :: options
        integer(I8), intent(out)            :: small_tria

        ! Auxiliary
        integer(I8) :: i, j, nb, ic, indmax, int_face, ic_neig
        integer(I8), allocatable :: cvLookUp(:), indCv(:), b_trias(:), &
            tf_int(:), cvs(:), neigs(:), tf(:)
        real(R8) :: ratio_s, ratio_psi, ratio_surf
        real(R8), allocatable :: dpsi(:)


        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        ! Initialize
        small_tria = 0
        cvLookUp = GetCvLookUpGA(c)
        indCv = (/ (i, i= 1,c%ntot)/)

        ! Find all boundary triangles

        nb = count(c%vertP2%Get() == 3 .and. c%cflags%Get() == 3)
        allocate(b_trias(nb))
        b_trias = pack(indCv,c%vertP2%Get() == 3 .and. c%cflags%Get() == 3)

        do i = 1, nb

            ! Get cell number
            ic = b_trias(i)

            ! Find the internal face of the triangle which is no boundary face
            tf = GetCellFaceGA(c, ic)

            allocate(tf_int(count(.not.isBoundaryFaceGA(f, tf)))) ! Note: in matlab using c%face to determine this, not the face labels
            tf_int = pack(tf,.not.isBoundaryFaceGA(f, tf))

            ! Find the poloidal face
            dpsi = abs( v%psi%Get(f%vert1%Get(tf_int)) - v%psi%Get(f%vert2%Get(tf_int)))
            indmax = maxloc(dpsi,1)
            int_face = tf_int(indmax)

            ! Get the cells and its neighbouring cells to compute criteria
            cvs = GetFaceCellGA(c, int_face, cvLookUp)

            do j = 1, size(cvs)
                if (cvs(j) /= ic) ic_neig = cvs(j)
            end do

            neigs = GetCellNeigsGA(grid, ic, cvLookUp)

            ! Criteria
            if ((c%reg%Get(ic) == c%reg%Get(ic_neig)) &
                .and. (minval(c%vertP2%Get(neigs)) .gt. 3)) then

                ! Size citerion
                ratio_s = qm%h_pol(ic) / qm%h_pol(ic_neig)

                ! Width criterion - if tria spans the whole flux tube => other removal method needed
                ratio_psi = qm%h_rad_psi(ic) / qm%h_rad_psi(ic_neig)

                ! Surface criterion
                ratio_surf = qm%cvS(ic) / sum(qm%cvS(neigs))

                ! Selection
                if (ratio_psi .lt. 0.9) then

                    if ((ratio_s .lt. options%cut_off_pol) &
                        .or. (ratio_surf .lt. options%cut_off_surf) &
                        .or. (ratio_psi .lt. 0.4)) then

                            small_tria = ic
                            return
                    end if

                end if

            end if

            ! Housekeeping
            deallocate(tf_int)

        end do

        end associate

    end subroutine

    subroutine DetermineRemovalMethod(grid, qm, small_tria, faceA, faceB, faceC, plain_merging)

        ! Description
        !============
        ! Determines which method to use for small triangle removal.

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)     :: grid
        type(QualityMetricUDT), intent(in)  :: qm
        integer(I8), intent(in)             :: small_tria
        integer(I8), intent(out)            :: faceA, faceB, faceC
        logical, intent(out)                :: plain_merging

        ! Auxiliary
        integer(I8) :: i, ic,  indmin, neigA
        integer(I8), allocatable, dimension(:) :: faceA_dummy, faceB_dummy, &
            tf, cvLookUp, cvs, forbidden_fcs
        real(R8), allocatable :: dpsi_f(:)
        logical, allocatable :: is_aligned(:), is_boundary(:)


        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        ! Initialize
        ic  = small_tria
        cvLookUp = GetCvLookUpGA(c)

        ! Get faces
        tf = GetCellFaceGA(c, ic)

        ! Find align face of the tria
        is_aligned = f%aligned%Get(tf) == 1;

        if ((count(is_aligned)) /= 1) then

            ! Most aligned faces
            dpsi_f = abs(v%psi%Get(f%vert1%Get(tf)) - v%psi%Get(f%vert2%Get(tf)))
            indmin = minloc(dpsi_f,1)
            faceA = tf(indmin)

        else

            allocate(faceA_dummy(count(is_aligned)))
            faceA_dummy = pack(tf, is_aligned)
            faceA = faceA_dummy(1)
            deallocate(faceA_dummy)

        end if

        ! Find boundary face of the triangle
        is_boundary = isBoundaryFaceGA(f, tf)

        ! Get boundary face
        allocate(faceB_dummy(count(is_boundary)))
        faceB_dummy = pack(tf, is_boundary)
        faceB = faceB_dummy(1)
        deallocate(faceB_dummy)

        ! Get other face
        do i = 1, size(tf)
            if ((tf(i) /= faceA) .and.(tf(i) /= faceB)) faceC = tf(i)
        end do

        ! Get neighboring cells at the aligned face
        cvs = GetFaceCellGA(c, faceA, cvLookUp)
        do i = 1, size(cvs)
            if (cvs(i) /= ic) neigA = cvs(i)
        end do

        ! Face length criterion
        plain_merging = .false.
        if ((qm%fcS(faceA) * 1.5_R8) .gt. qm%h_pol(neigA)) then

            ! Check whether faceC is not a region boundary face
            cvs = GetFaceCellGA(c, faceC, cvLookUp)

            if (c%reg%Get(cvs(1)) == c%reg%Get(cvs(2))) then

                ! Check whether faceC is not a forbidden face
                call grid%GetForbiddenMergeFaces(forbidden_fcs)
                if (.not.any(faceC == forbidden_fcs)) then
                    plain_merging = .true.
                end if

            end if

        end if

        end associate

    end subroutine

    subroutine LocalSmallTriangleRemoval(grid, cell, faceA, faceB, faceC, options)

        ! Description 
        !============
        ! Removes the small triangle on the classical method
        ! - faceA is the aligned face
        ! - faceB is the boundary face
        ! - faceC is the other (usually poloidal) face of the triangle
        ! checks

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT)        :: grid
        integer(I8), intent(in) ::  cell, faceA, faceB, faceC
        type(GAoptionsUDT)      :: options

        ! Auxiliary
        integer(I8) :: j, k, n, n_oc, s, vx_common, indexf, cellA, &
            vertsA(2), vertsB(2), cellC, indf, indv
        integer(I8), allocatable :: cell_rem(:), &
            tv_cellA(:), tv_cellC(:), new_verts_cellC(:),&
            range(:), new_verts_cellA(:), faces(:), vx_rem(:)

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert, &
            fd => grid%data%fluxdata &
            )

        ! Check if flux surface surfaces are unique
        if (options%debug) call CheckUniqueness(fd%fluxsurfacefaces)
        

        ! Find vertex which is both in faceA and faceB
        vertsA = [f%vert1%Get(faceA), f%vert2%Get(faceA)]
        vertsB = [f%vert1%Get(faceB), f%vert2%Get(faceB)]

        allocate(vx_rem(2))
        vx_rem = 0

        if (any(vertsB == vertsA(1))) then
            vx_common = vertsA(1)
            vx_rem(1) = vertsA(2)
        else
            vx_common = vertsA(2)
            vx_rem(1) = vertsA(1)
        end if

        do j = 1, 2
            if (vertsB(j) /= vx_common) vx_rem(2) = vertsB(j)
        end do

        ! META state starts here
        allocate(cell_rem(1))
        cell_rem = cell
        call grid%RemoveCells(cell_rem)

        ! Aligned cell
        !-------------
        ! Find the cell at the aligned face
        indexf = findloc(c%face%Get(), faceA, 1)
        cellA = GetCellFromFaceIndex(c, indexf)

        !n_f = indexf - c%faceP1%Get(cellA) + 1

        !tf_cellA = GetCellFaceGA(c, cellA)

        ! Remove the aligned face in that cell and set the remaining faces => become triangle
        call c%face%Remove(indexf)
        call c%faceP2%SumMask(cellA, -1)
        range = (/ (j, j = cellA+1, c%faceP1%Size())/)
        call c%faceP1%SumMask(range, -1)
        !allocate(new_face_cellA(3))
        !new_face_cellA = [tf_cellA(1:nf_1), tf_cellA(n_f+1:size(tf_cellA))]
        !s = grid%cell%faceP1%Get(cellA)
        !n = grid%cell%faceP2%Get(cellA)
        !call grid%cell%face%Replace(s,s+n-1,new_face_cellA)

        ! Remove the vertex of the aligned face in that cell
        tv_cellA = GetCellVertGA(c, cellA)
        indv = findloc(tv_cellA, vx_rem(1), 1) ! ordering of vertices matters!!

        new_verts_cellA = [tv_cellA(1:indv-1), tv_cellA(indv+1:size(tv_cellA)) ]
        s = c%vertP1%Get(cellA)
        n = c%vertP2%Get(cellA)
        call c%vert%Replace(s,s+n-1,new_verts_cellA) 

        call c%vertP2%SumMask(cellA, -1)
        range = (/ (j, j = cellA+1, c%vertP1%Size())/)
        call c%vertP1%SumMask(range, -1)

        ! Cell at non-aligned and non-boundary face
        !------------------------------------------
        ! Find the cell at the non-boundary faces and non-aligned face
        indexf = findloc(c%face%Get(),faceC, 1)
        cellC = GetCellFromFaceIndex(c, indexf)

        call c%face%Remove(indexf)
        call c%faceP2%SumMask(cellC, -1)
        range = (/ (j, j = cellC+1, c%faceP1%Size())/)
        call c%faceP1%SumMask(range, -1)


        ! Remove the vertex of the aligned face in that cell
        tv_cellC = GetCellVertGA(c, cellC)
        indv = findloc(tv_cellC, vx_rem(1), 1) 

        new_verts_cellC = [tv_cellC(1:indv-1), tv_cellC(indv+1:size(tv_cellC))]
        s = c%vertP1%Get(cellC)
        n = c%vertP2%Get(cellC)        
        call c%vert%Replace(s,s+n-1,new_verts_cellC)
        call c%vertP2%SumMask(cellC, -1)
        range = (/ (j, j = cellC+1, c%vertP1%Size())/)
        call c%vertP1%SumMask(range, -1)

        ! Find cells where vx_rem(1:2) is and replace is by the common vertex (vx_com)
        do j = 1, 2

            ! cell%vert
            indv = findloc(c%vert%Get(), vx_rem(j), 1)
            call c%vert%Set(indv, vx_common)
            
            ! face%vert
            n_oc = count(f%vert1%Get() == vx_rem(j))
            do k = 1,n_oc
                indf = findloc(f%vert1%Get(), vx_rem(j), 1)
                call f%vert1%Set(indf, vx_common)
            enddo

            n_oc = count(f%vert2%Get() == vx_rem(j))
            do k = 1,n_oc
                indf = findloc(f%vert2%Get(), vx_rem(j), 1)
                call f%vert2%Set(indf, vx_common)
            enddo

        end do
        
        ! Check flux surface faces
        if (options%debug) call CheckUniqueness(fd%fluxsurfacefaces)
        
        ! Now faces of small triangle can get removed
        faces = [faceA, faceB, faceC]
        call grid%RemoveFaces(faces)

        ! Now vertices can get removed
        call grid%RemoveVertices(vx_rem)


        ! Check flux surface faces
        if (options%debug) call CheckUniqueness(fd%fluxsurfacefaces)

        call grid%GetFsVxFromFsFc(options)


        end associate

    end subroutine

    ! Stacked Triangles
    !==================
    subroutine StackedTrias(grid, magneticField, qm, options)

        ! Description
        !============
        ! Transforms a cut-cell triangles to stacked triangles
        ! Con_vert is the vertex to which all triangles are attached

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)         :: grid
        type(MagneticFieldUDT), intent(in)      :: magneticField
        type(QualityMetricUDT), intent(inout)   :: qm
        type(GAoptionsUDT), intent(in)          :: options

        ! Auxiliary
        integer(I8) :: i, tria, con_vert
        integer(I8), allocatable, dimension(:) :: cctria, cctraps, traps, ctest1, ctest2
        integer(I8), allocatable :: cctrapsP(:,:)
        logical :: approved

        ! Printing
        print *, 'Apply Stacked Triangle adaptation'

        ! Calculating quality metric
        call qm%CalculateQualityMetrics(grid, options, magneticField, .false., .false.)

        ! Detection
        !==========
        ! Number of cells should not change, so try to detect whole group of cells
        call grid%DetectCutCell(qm, options, cctria, cctraps, cctrapsP)

        ! For visualization
        if (options%debug) then
            call WriteArray(cctria, 'cctria')
            call WriteArray(cctraps, 'cctraps')
            ctest1 = cctrapsP(:,1)
            call WriteArray(ctest1, 'cctrapsP1')
            ctest2 = cctrapsP(:,2)
            call WriteArray(ctest2, 'cctrapsP2')
        end if


        ! Adaptation
        !===========
        do i = 1, size(cctria)

            ! Unpack
            tria = cctria(i)
            traps = cctraps(cctrapsP(i,1):cctrapsP(i,1)+cctrapsP(i,2)-1)

            ! Get the connection vertex
            call grid%GetConnectionVertex(tria, traps, con_vert)

            ! Check cell overlap when transformation would be done
            if (con_vert /= 0) &
                call grid%CheckStackedTriaOverlap(tria, traps, con_vert, approved)
            
            !if (options%debug) then
            !    print *, con_vert
            !    print *, 'x-coor: ' , grid%vert%x%Get(con_vert)
            !    print *, 'y-coor: ' , grid%vert%y%Get(con_vert)
            !    print *, approved
            !end if

            ! Do the adaptation if oké
            if (con_vert /= 0 .and. approved) &
                call grid%StackAdaptation(tria, traps, con_vert)

        end do

        ! Merge some triangles if to shallow angle
        if (options%merge_stacked_trias) then
            ! call grid%MergeStackedTrias() - TODO

        end if

        ! Check from boundary trapezoids which should be included in the stacked triangle
        if (options%merge_trap_into_stacked) then
            ! call grid%MergeTrapIntoStacked() - TODO
        end if


        ! Recalculate the magneticfield
        if (.not. options%slab) call grid%RecalcMagn(magneticField)

        ! Printing
        print *, 'Ended Stacked Triangle adaptation'
        
    end subroutine

    subroutine DetectCutCell(grid, qm, options, cctria, cctraps, cctrapsP)

        ! Description
        !============
        ! Detecting cutcells in the grid to be transformed to stacked triangles

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)         :: grid
        type(QualityMetricUDT), intent(in)      :: qm
        type(GAoptionsUDT)                      :: options 
        integer(I8), allocatable, intent(out)   :: cctria(:), cctraps(:), cctrapsP(:,:)
        
        ! Auxiliary
        integer(I8) :: i, ic, counter, trap1, counter_tria, counter_trapgroups, s
        integer(I8), allocatable :: cctrapsPD(:,:)
        integer(I8), allocatable, dimension(:) :: cctriaD, cctrapsD, trias, &
            triasD, cvLookUp, vxs, indCv, traps
        real(R8) :: dpsi_max
        real(R8), allocatable :: rpolrad(:)
        logical, allocatable :: in_cutcell(:)

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        ! Initialize
        allocate(cctriaD(c%ntot))
        allocate(cctrapsD(c%ntot))
        allocate(cctrapsPD(c%ntot,2))
        allocate(in_cutcell(c%ntot))
        allocate(traps(c%ntot))
        counter = 0
        counter_tria = 0
        counter_trapgroups = 0
        cctriaD = 0
        cctrapsD = 0 
        cctrapsPD = 0 
        traps = 0
        in_cutcell = .false.
        cvLookUp = GetCvLookUpGA(c)

        rpolrad = qm%h_pol / qm%h_rad

        ! Get triangles 
        indCv = (/ (i, i = 1, c%ntot)/)
        if (options%stacked_trias_checkAR) then

            allocate(triasD(count(c%faceP2%Get() == 3)))
            triasD = pack(indCv, c%faceP2%Get() == 3)
            allocate(trias(count(rpolrad(triasD) .lt. options%stacked_trias_maxAR)))
            trias = pack(triasD,rpolrad(triasD) .lt. options%stacked_trias_maxAR )
            deallocate(triasD)

        else

            allocate(trias(count(c%faceP2%Get() == 3)))
            trias = pack(indCv, c%faceP2%Get() == 3)

        end if

        ! Loop over triangles
        do i = 1, size(trias)

            ! Get cell
            ic = trias(i)
            
            ! Only boundary cells
            if (c%cflags%Get(ic) == 3) then

                ! Get vertices and dpsi_max
                vxs = GetCellVertGA(c, ic)
                dpsi_max = abs(maxval(v%psi%Get(vxs)) - minval(v%psi%Get(vxs)))

                ! Get the first trapezoidal cells on top
                call grid%GetFirstTrap(ic, trap1, cvLookUp)

                ! Continue if found
                if (trap1 /= 0) then

                    ! Add to traps array
                    counter = counter + 1
                    traps(1) = trap1

                    ! Get the other traps above
                    call grid%GetNextTraps(ic, cvLookUp, traps, counter)

                    ! Check if there is no Xpoint to near
                    if (any(isMember(vxs,grid%data%xpointID))) then
                        traps = 0 ! Reset
                        counter = 0
                    end if

                    ! Check if some trapezoids where already in another 
                    ! group of cutcells
                    if (any(in_cutcell(traps(1:counter)))) then
                        traps = 0
                        counter = 0
                    end if


                    if (counter /= 0) then

                        ! Save 
                        counter_tria = counter_tria + 1
                        cctriaD(counter_tria) = ic

                        counter_trapgroups = counter_trapgroups + 1
                        if (counter_trapgroups == 1) then
                            cctrapsPD(1,1) =   1
                        else
                            cctrapsPD(counter_trapgroups,1) =   &
                                cctrapsPD(counter_trapgroups-1,1) + cctrapsPD(counter_trapgroups-1,2)
                        end if
                        cctrapsPD(counter_trapgroups,2) = counter
                        s = cctrapsPD(counter_trapgroups,1)
                        cctrapsD(s:s+counter-1) = traps(1:counter)


                        ! Save in_cutcell
                        in_cutcell(ic) = .true.
                        in_cutcell(traps(1:counter)) = .true.
                        
                    end if

                end if 

            end if

        end do

        ! Trim
        cctria = cctriaD(1:counter_trapgroups)
        cctrapsP = cctrapsPD(1:counter_trapgroups,:)
        cctraps = cctrapsD(1:cctrapsPD(counter_trapgroups,1)+cctrapsPD(counter_trapgroups,2)-1)

        end associate

    end subroutine

    subroutine GetFirstTrap(grid, ic, trap1, cvLookUp)

        ! Description
        !============
        ! Find the first trapezoidal cells of a triangle

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(in)            :: grid
        integer(I8), intent(in)                 :: ic
        integer(I8), intent(in)                 :: cvLookUp(:)
        integer(I8), intent(out)                :: trap1

        ! Auxiliary
        integer(I8) :: j, counter_b, common_face, indmin
        integer(I8), allocatable, dimension(:) :: vxs, neigs, b_neig, b_neig2, &
            fcs_neig, b_face_neig, bvxs, bface_neig_verts
        real(R8) :: dpsi_max, dpsi_b
        real(R8), allocatable :: dpsi(:) 

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        ! Get vertices and dpsi_max
        vxs = GetCellVertGA(c, ic)
        dpsi_max = abs(maxval(v%psi%Get(vxs)) - minval(v%psi%Get(vxs)))   

        ! Get the first trapezoidal cells on top
        neigs = GetCellNeigsGA(grid, ic, cvLookUp)
        allocate(b_neig(count(c%cflags%Get(neigs) == 3)))
        b_neig = pack(neigs, c%cflags%Get(neigs) == 3) ! Boundary cells only

        counter_b = 0 
        allocate(b_neig2(size(b_neig)))
        b_neig2 = 0
        do j = 1, size(b_neig)

            ! No triangles
            if (c%vertP2%Get(b_neig(j)) .gt. 3) then

                common_face = GetCommonFace(c, b_neig(j), ic)
                fcs_neig = GetCellFaceGA(c, b_neig(j))
                allocate(b_face_neig(count(isBoundaryFaceGA(f, fcs_neig))))
                b_face_neig = pack(fcs_neig, isBoundaryFaceGA(f, fcs_neig))

                ! Check for common vertex with triangle
                allocate(bface_neig_verts(2))
                bface_neig_verts(1) = f%vert1%Get(b_face_neig(1))
                bface_neig_verts(2) = f%vert2%Get(b_face_neig(1))


                ! Dpsi of the neighbor
                bvxs = GetCellVertGA(c, b_neig(j))
                dpsi_b = abs(maxval(v%psi%Get(bvxs)) - minval(v%psi%Get(bvxs)))

                if ((f%aligned%Get(common_face) == 0) &
                    .and. (f%aligned%Get(b_face_neig(1)) == 0) &
                    .and. (any(isMember(bface_neig_verts, vxs))) &
                    .and. (dpsi_b .gt. (dpsi_max* 1.001_R8)) ) then
                    
                    ! Keep the cell
                    counter_b = counter_b + 1
                    b_neig2(counter_b) = b_neig(j)

                end if

            end if

        end do

        ! If multiple neighbors are left, choose the one with the most
        ! similar psi value
        if (counter_b .ge. 2) then

            dpsi = abs(c%psi%Get(ic) - c%psi%Get(b_neig2(1:counter_b)))
            indmin = minloc(dpsi,1)
            trap1 = b_neig2(indmin) 

        elseif (counter_b == 1) then

            trap1 = b_neig2(1)

        else if (counter_b == 0) then

            trap1 = 0

        end if


        end associate

    end subroutine
    
    subroutine GetNextTraps(grid, ic, cvLookUp, traps, counter)

        ! Description
        !============
        ! Find the next trapezoidal cells, given the triangle and the first trapezoid

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(in)            :: grid
        integer(I8), intent(in)                 :: ic
        integer(I8), intent(in)                 :: cvLookUp(:)
        integer(I8), intent(inout)              :: traps(:)
        integer(I8), intent(inout)              :: counter
 
        ! Auxiliary
        integer(I8), allocatable, dimension(:) :: neigs, b_neig, b_verts, b_faces, &
            trap_verts, trap_faces, bface_trap, bvert_trap, bvert_neig, bfaces_neig, &
            verts_last, trap_cells
        integer(I8) :: trap_next, nc_vert, trap_last, indmin
        real(R8), allocatable, dimension(:) :: dpsi
        real(R8) :: dpsi_max, dpsi_max_b, b_neig_psi, min_psi, max_psi
        logical, allocatable :: log(:)

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        ! Initialize
        trap_last = traps(1)
        verts_last = GetCellVertGA(c, trap_last)
        dpsi_max = abs(maxval(v%psi%Get(verts_last)-minval(v%psi%Get(verts_last))))
        counter = 1
        do while (trap_last /= 0)

            ! Get neigs
            neigs = GetCellNeigsGA(grid, trap_last, cvLookUp)
            trap_cells = traps(1:counter)

            ! Get b_neig
            log = (c%cflags%Get(neigs) == 3 &
                .and. c%faceP2%Get(neigs) == 4 &
                .and. .not.isMember(neigs, trap_cells))

            allocate(b_neig(count(log)))
            b_neig = pack(neigs,log)

            if (size(b_neig) .ge. 2) then

                dpsi = abs(c%psi%Get(trap_last) - c%psi%Get(b_neig))
                indmin = minloc(dpsi,1)
                trap_next = b_neig(indmin)

            else if (size(b_neig) == 1) then

                trap_next = b_neig(1)
            
            else 

                trap_next = 0

            end if


            ! Check extra criteria
            if (trap_next /= 0) then

                ! Criteria based on psi values
                b_verts = GetCellVertGA(c, trap_next)
                b_faces = GetCellFaceGA(c, trap_next)
                dpsi_max_b = abs(maxval(v%psi%Get(b_verts)) - minval(v%psi%Get(b_verts)))
    
                b_neig_psi = c%psi%Get(trap_next)
                trap_verts = GetCellVertGA(c, trap_last)
                min_psi = minval(v%psi%Get(trap_verts))
                max_psi = maxval(v%psi%Get(trap_verts))

                ! Criteria based on connectivity - boundary faces need 
                ! to have at least a common vertex
                trap_faces = GetCellFaceGA(c, trap_last)
                allocate(bface_trap(count(f%label%Get(trap_faces) /= 0)))
                bface_trap = pack(trap_faces, f%label%Get(trap_faces) /= 0)
                bvert_trap = GetVxsFromFcsGA(f, bface_trap) 

                allocate(bfaces_neig(count(f%label%Get(b_faces) /= 0)))
                bfaces_neig = pack(b_faces, f%label%Get(b_faces) /= 0)
                bvert_neig = GetVxsFromFcsGA(f, bfaces_neig)
                nc_vert = count(isMember(bvert_neig, bvert_trap))

                ! Apply criteria
                if ( (dpsi_max_b .gt. (dpsi_max*1.01_R8)) &
                     .and. (b_neig_psi .gt. min_psi) &
                     .and. (b_neig_psi .lt. max_psi) &
                     .and. (nc_vert .gt. 0)) then

                        ! Add cell
                        counter = counter + 1
                        traps(counter) = trap_next
                        trap_last = trap_next
                        dpsi_max = dpsi_max_b

                else
                    
                    ! Exit the loop
                    trap_last = 0
                    
                end if

                ! Housekeeping
                deallocate(bface_trap)
                deallocate(bfaces_neig)
                
            else

                ! No trap_next found so exit the loop
                trap_last = 0

            end if

            deallocate(b_neig)

        end do

        end associate

    end subroutine

    subroutine GetConnectionVertex(grid, tria, traps, con_vert)

        ! Description
        !============
        ! Finds the connection vertex to which all stacked triangles
        ! need to be connected, based on information on cutcells

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(in)  :: grid
        integer(I8), intent(in)       :: tria
        integer(I8), intent(in)       :: traps(:)
        integer(I8), intent(out)      :: con_vert

        ! Auxiliary
        integer(I8) :: i, k, big_trap, indmin, indmax, lt, vxs2(2)
        integer(I8), allocatable, dimension(:) :: vxs, cvs, cv2, &
             vx, bfcs, ts, vxs_lt, fcs

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        ! Initialize
        con_vert = 0
        big_trap = traps(size(traps))
        vxs = GetCellVertGA(c, big_trap)
        indmin = minloc(v%psi%Get(vxs),1)
        indmax = maxloc(v%psi%Get(vxs),1)

        do i = 1, size(vxs)

            ! Must be a boundary vertex
            if (isBoundaryVertGA(grid, vxs(i))) then

                ! Get cells
                cvs = GetVertCellGA(c, vxs(i))
                if ((size(cvs) .gt. 2) .and. (i == indmax .or. i == indmin)) then
                    con_vert = vxs(i)
                    exit
                end if

            end if

        end do

        ! To capture con_vert of cut cells in the outermost fluxsurface flux tube
        if (con_vert == 0) then

            ts = [ tria, traps ]
            lt = ts(size(ts)-1)
            vxs_lt = GetCellVertGA(c, lt)

            do i = 1, size(vxs)

                if (isBoundaryVertGA(grid, vxs(i))) then

                    ! Get cells
                    cvs = GetVertCellGA(c, vxs(i))
                    allocate(cv2(count(cvs /= big_trap)))
                    cv2 = pack(cvs, cvs /= big_trap)

                    fcs = GetVertFaceGA(f, vxs(i))
                    allocate(bfcs(count(isBoundaryFaceGA(f, fcs))))
                    bfcs = pack(fcs, isBoundaryFaceGA(f, fcs))

                    ! If the other cells are not part of the cutcells
                    if (.not.any(isMember(cv2, ts))) then

                        do k = 1, size(bfcs)
                            vxs2(1) = f%vert1%Get(bfcs(k)) 
                            vxs2(2) = f%vert2%Get(bfcs(k)) 
                            allocate(vx(count(vxs2 /= vxs(i))))
                            vx = pack(vxs2, vxs2 /= vxs(i))

                            if (any(isMember(vx, vxs_lt))) then
                                con_vert = vxs(i)
                                exit
                            end if

                            deallocate(vx)

                        end do
                    end if

                    if (con_vert /= 0) exit

                    deallocate(cv2)
                    deallocate(bfcs)

                end if

            end do

        end if

        end associate

    end subroutine

    subroutine CheckStackedTriaOverlap(grid, tria, traps, con_vert, approved)

        ! Description
        !============
        ! % Check whether it is allow to transform to stacked triangle
        ! Method
        ! 1) make the vector between the tria_tip and the con_vert
        ! 2) get the vertices of the trapezoids which are not boundary vertices
        ! 3) make vectors from the tria_tip and these vertices
        ! 4) the sine between the vector from tria_tip and con_vert, and the
        ! vectors to the trapezoid vertices, should be bigger than 0

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(in)    :: grid
        integer(I8), intent(in)         :: tria, con_vert
        integer(I8), intent(in)         :: traps(:)
        logical, intent(out)            :: approved

        ! Auxiliary
        integer(I8) :: nt, counter_v, j, prev_icv, ind, nvix, nvt
        integer(I8), allocatable, dimension(:) :: verts_trap, vxs, fcs_al, &
            vix, tria_tip, vxs_trap1, vxs_tria, verts_trapUD, fcs, vtrap
        real(R8) :: vec_outx, vec_outy, norm_vec_out
        real(R8), allocatable :: vec_inx(:), vec_iny(:), norm_vec_in(:), sint(:)
        logical, allocatable :: log(:)
        type(IntegerDynamicArrayUDT), allocatable :: verts_trapU

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        ! Initialize
        nt = size(traps)
        allocate(verts_trap(nt*2))
        counter_v = 0
        approved = .true.

        ! Loop over the trapezoids
        do j = 1, nt

            ! Get internal vertices
            vxs = GetCellVertGA(c, traps(j))
            log = (.not.isBoundaryVertGA(grid, vxs))
            nvix = count(log)

            if (nvix /= 2) then

                ! Take the vertices of the aligned faces
                fcs = GetCellFaceGA(c, traps(j))
                allocate(fcs_al(count(f%aligned%Get(fcs) == 1)))
                fcs_al = pack(fcs, f%aligned%Get(fcs) == 1)

                allocate(vix(2))
                vix(1) = f%vert1%Get(fcs_al(1))
                vix(2) = f%vert2%Get(fcs_al(1))


                deallocate(fcs_al)

            else 

                ! Take the non-boundary vertices
                allocate(vix(count(log)))
                vix = pack(vxs, log)

            end if

            ! Add vertices to verts_trap
            verts_trap(counter_v + 1: counter_v + size(vix)) = vix
            counter_v = counter_v + nvix

            deallocate(vix)

        end do

        ! Unique on vertstrap
        call Unique(verts_trap, verts_trapUD)
        allocate(IntegerDynamicArrayUDT :: verts_trapU)
        verts_trapU = ConstructIntegerDynamicArray(verts_trapUD)

        ! Loop over the trapezoids
        do j = 1, nt

            if (j == 1) then

                ! Get the tria vert = tip_vert
                vxs_tria = GetCellVertGA(c, tria)
                vxs_trap1 = GetCellVertGA(c, traps(1))
                
                log = (.not.isMember(vxs_tria, vxs_trap1))
                allocate(tria_tip(count(log)))
                tria_tip = pack(vxs_tria, log)

            else if (j == 2) then

                ! Tria tip has to become the next vertex
                prev_icv = tria
                vxs = GetCellVertGA(c, prev_icv)
                vtrap = verts_trapU%Get()
                allocate(tria_tip(count(isMember(vxs,vtrap))))
                tria_tip = pack(vxs, isMember(vxs,vtrap))

                ! And remove it from the verts_trap
                ind = findloc(verts_trapU%Get(), tria_tip(1), 1)
                call verts_trapU%Remove(ind) 
                
            else

                ! Tria tip has to become the next vertex
                prev_icv = traps(j-2)
                vxs = GetCellVertGA(c, prev_icv)
                vtrap = verts_trapU%Get()
                allocate(tria_tip(count(isMember(vxs,vtrap))))
                tria_tip = pack(vxs, isMember(vxs,vtrap))

                ! And remove it from the verts_trap
                ind = findloc(verts_trapU%Get(), tria_tip(1), 1)
                call verts_trapU%Remove(ind) 

            end if

            ! Make outer vector between the tria_tip and 
            vec_outx = v%x%Get(con_vert) -  v%x%Get(tria_tip(1))
            vec_outy = v%y%Get(con_vert) -  v%y%Get(tria_tip(1))
            norm_vec_out = sqrt(vec_outx**2 + vec_outy**2)

            ! Make vectors
            vec_inx = v%x%Get(verts_trapU%Get()) -  v%x%Get(tria_tip(1))
            vec_iny = v%y%Get(verts_trapU%Get()) -  v%y%Get(tria_tip(1))
            norm_vec_in = sqrt(vec_inx**2 + vec_iny**2)

            ! Compute sine for each case
            ! Calculate angles (sin = |a x b| / norm(a)*norm(b))
            sint = (vec_outx * vec_iny - vec_inx * vec_outy) / (norm_vec_out * norm_vec_in )

            ! If all sines are positive or all are negative, and first triangle
            ! not to skewed
            nvt = verts_trapU%Size()
            if (.not.((count(sint .lt. 0) == nvt) .or. (count(sint .gt. 0) == nvt))) then
                approved = .false.
                exit
            end if  
            
            ! Housekeeping
            deallocate(tria_tip)

        end do

        end associate

    end subroutine 

    subroutine StackAdaptation(grid, tria, traps, con_vert)

        ! Description
        !============
        ! Perform the adaptation from cutcell to stacked triangles

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout) :: grid
        integer(I8), intent(in)         :: tria, con_vert
        integer(I8), intent(in)         :: traps(:)

        ! Auxiliary
        integer(I8) :: i, j, ic, common_face, b_vert, counterv, &
            counterf, nv, nvT
        integer(I8), allocatable, dimension(:) :: trapsC, facesQ, &
            b_facesQD, b_facesQ, vertsQ, indbvert, indbvert1, indbvert2, &
            indcv, indfc, ar, ar1, ar2, vx_remD, vx_rem, fc_remD, fc_rem, &
            vertsT, loc1, facesQD, vertsQD
        real(R8) :: cx, cy
        logical, allocatable :: log(:)

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        ! Initialize
        trapsC = [tria, traps]
        allocate(vx_remD(v%ntot))
        allocate(fc_remD(f%ntot))
        vx_remD = 0
        fc_remD = 0
        counterv = 0
        counterf = 0

        ! Begin with the last trap

        do i = 1, (size(trapsC)-1)

            ic = trapsC(size(trapsC)+1-i)
            facesQ = GetCellFaceGA(c, ic)
            allocate(b_facesQD(count(isBoundaryFaceGA(f, facesQ))))
            b_facesQD = pack(facesQ, isBoundaryFaceGA(f, facesQ))

            if (size(b_facesQD) .gt. 1) then
                allocate(b_facesQ(count(f%aligned%Get(b_facesQD) == 0)))
                b_facesQ = pack(b_facesQD, f%aligned%Get(b_facesQD) == 0)
            else
                allocate(b_facesQ(size(b_facesQD)))
                b_facesQ = b_facesQD
            end if

            ! Get common face with next trap
            common_face = GetCommonFace(c, ic, trapsC(size(trapsC)-i))

            b_vert = GetCommonVert(f, common_face, b_facesQ(1))

            ! Replace the occurrence of bvert with con_vert in cell%vert
            log = (c%vert%Get() == b_vert)
            allocate(indbvert(count(log)))
            indcv = (/ (j, j = 1, c%vertP1%Get(c%ntot)+c%vertP2%Get(c%ntot)-1) /)
            indbvert = pack(indcv,log)
            ar = (/ (con_vert, j = 1, size(indbvert) ) /)
            call c%vert%Set(indbvert, ar)

            ! Replace the occurrence of bvert with con_vert in face%vert
            allocate(indbvert1(count(f%vert1%Get() == b_vert)))
            allocate(indbvert2(count(f%vert2%Get() == b_vert)))
            indfc = (/ (j, j = 1, f%ntot) /)
            indbvert1 = pack(indfc, f%vert1%Get() == b_vert)
            indbvert2 = pack(indfc, f%vert2%Get() == b_vert)
            ar1 = (/ (con_vert, j = 1, size(indbvert1) ) /)
            ar2 = (/ (con_vert, j = 1, size(indbvert2) ) /)
            call f%vert1%Set(indbvert1, ar1)
            call f%vert2%Set(indbvert2, ar2)

            ! Remove bvert out of cell%vert
            vertsQD = GetCellVertGA(c, ic)
            call Unique(vertsQD, vertsQ)
            call c%vert%Replace(c%vertP1%Get(ic),c%vertP1%Get(ic)+c%vertP2%Get(ic)-1, vertsQ)
            nv = c%vertP2%Get(ic) - 1
            call c%vertP2%Set(ic, nv)
            loc1 = (/ (j, j = ic+1 , c%vertP1%Size()) /)
            call c%vertP1%SumMask(loc1, -1)

            ! Adjust centroid
            cx = sum(v%x%Get(vertsQ))/(real(nv, kind=R8))
            cy = sum(v%y%Get(vertsQ))/(real(nv, kind=R8))
            call c%x%Set(ic, cx)
            call c%y%Set(ic, cy)

            ! Remove boundaryFace out of cell%face
            allocate(facesQD(count(facesQ /= b_facesQ(1))))
            facesQD = pack(facesQ, facesQ /= b_facesQ(1))
            call c%face%Replace(c%faceP1%Get(ic),c%faceP1%Get(ic)+c%faceP2%Get(ic)-1, facesQD)
            call c%faceP2%SumMask(ic, -1)
            loc1 = (/ (j, j = ic+1 , c%faceP1%Size()) /)            
            call c%faceP1%SumMask(loc1, -1)

            ! Save vertices and face to remove
            counterv = counterv + 1
            counterf = counterf + 1
            vx_remD(counterv) = b_vert
            fc_remD(counterf) = b_facesQ(1)

            deallocate(b_facesQ)
            deallocate(b_facesQD)
            deallocate(indbvert)
            deallocate(indbvert1)
            deallocate(indbvert2)
            deallocate(facesQD)

        end do

        ! Trim
        fc_rem = fc_remD(1:counterf)
        vx_rem = vx_remD(1:counterv)

        ! Adjust centroid van tria
        nvT = c%vertP2%Get(tria)
        vertsT = GetCellVertGA(c, tria)
        cx = sum(v%x%Get(vertsT))/real(nvT, kind=R8)
        cy = sum(v%y%Get(vertsT))/real(nvT, kind=R8)
        call c%x%Set(tria, cx)
        call c%y%Set(tria, cy)

        ! Remove faces & vertices
        call grid%RemoveFaces(fc_rem)
        call grid%RemoveVertices(vx_rem)

        ! Determine cflags
        call grid%DetermineCflags(trapsC)

        end associate

    end subroutine

    ! Stacked to Cutcell
    !===================
    subroutine StackedToCutcell(grid, magneticField, options)

        ! Description
        !============
        ! Convert stacked triangle back to cutcells. This will give better
        ! possiblities for merging and splitting.
        ! The boundary face will be uniformly splitted in the right number of cells

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)     :: grid
        type(MagneticFieldUDT), intent(in)  :: magneticField
        type(GAoptionsUDT), intent(in)      :: options

        ! Auxiliary
        integer(I8) :: i, ic, nc, nv, v1, v2, start_vertex, end_vertex, ind, &
            bfcs, lbl, vxs(2), face_num, new_verts3(3), new_faces3(3), &
            new_verts4(4), new_faces4(4), fcsA, fcst_up_correct, fcst_al, &
            c_fcs
        integer(I8), allocatable, dimension(:) :: cvs, fcs, &
             cvs_v1, cvs_v2, new_vertex, vecV, ar,  &
            cvs1, cvs2, fcst_up,  new_bfcs, i_verts, new_ifcs, regs, &
            rem_faces, rem_cells, cells, cvLookUp, fcst,  &
            new_faces, new_verts
        real(R8) :: vec_x, vec_y
        real(R8), allocatable, dimension(:) :: v1_nx, v1_ny, v1_psi, &
            v1_bx, v1_by, fcA_length, fcA_length_int, &
            Vdistribution, inVdistribution, fcA_dist

        logical :: found

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        print *, 'Apply StackedToCutcell'

        ! Detection of stacked triangles
        call grid%DetectStackedTrias(found, cvs, nc)

        ! Adaptation
        do while (found)

            ! Get boundary face
            fcs = GetCellFaceGA(c, cvs(1))
            ind = findloc(isBoundaryFaceGA(f, fcs), .true. ,1)
            bfcs = fcs(ind)
            lbl = f%label%Get(bfcs)

            ! Save the regions
            regs = c%reg%Get()

            ! New vertices
            nv = nc - 1

            ! Initialize some arrays
            allocate(new_bfcs(nc))            
            allocate(i_verts(nv))
            allocate(new_ifcs(nv))

            ! Vector tangential to boundary face
            cvLookUp = GetCvLookUpGA(c)

            ! Start vertex
            v1 = f%vert1%Get(bfcs)
            v2 = f%vert2%Get(bfcs)

            cvs_v1 = GetVertCellGA(c, v1, cvLookUp)
            cvs_v2 = GetVertCellGA(c, v2, cvLookUp)

            if (count(isMember(cvs_v1,cvs)) == 1) then
                start_vertex = v1
                end_vertex = v2               
            else if (count(isMember(cvs_v2,cvs)) == 1) then
                start_vertex = v2
                end_vertex = v1                  
            end if

            vec_x = v%x%Get(end_vertex) - v%x%Get(start_vertex);
            vec_y = v%y%Get(end_vertex) - v%y%Get(start_vertex);

            new_vertex = (/ (i, i = v%ntot+1, v%ntot+nv) /)

            ! Get the vertex distribution along the line
            if (options%stacked_to_cutcell_nonuniform) then

                ! Non-uniform distribution based on sizes of aligned faces
                allocate(fcA_length(nc))
                allocate(fcA_length_int(nc))
                allocate(fcA_dist(nc))
                fcA_length = 0
                fcA_length_int = 0

                do i = 1, nc
                    fcs = GetCellFaceGA(c, cvs(i))
                    ind = findloc(f%aligned%Get(fcs), 1, 1)
                    fcsA = fcs(ind)

                    v1 = f%vert1%Get(fcsA)
                    v2 = f%vert2%Get(fcsA)

                    fcA_length(i) = sqrt( (v%x%Get(v1) - v%x%Get(v2))**2 &
                        + (v%y%Get(v1) - v%y%Get(v2))**2)

                    if (i == 1) then
                        fcA_length_int(i) = fcA_length(i)
                    else
                        fcA_length_int(i) = fcA_length(i) + fcA_length_int(i-1)
                    end if

                end do

                ! Save non-uniform distribution
                fcA_dist = fcA_length_int / sum(fcA_length)
                Vdistribution = fcA_dist(1:nv)
                inVdistribution = 1 - Vdistribution

                deallocate(fcA_length) 
                deallocate(fcA_length_int) 
                deallocate(fcA_dist)

            else

                ! Uniform distribution
                vecV = (/(i, i = 1, nv)/)
                Vdistribution = vecV / nc
                inVdistribution = 1 - Vdistribution

            end if

            ! Append the vertices
            v1_nx = v%x%Get(start_vertex) + vec_x * Vdistribution
            v1_ny = v%y%Get(start_vertex) + vec_y * Vdistribution
            v1_psi = v%psi%Get(start_vertex)*inVdistribution + v%psi%Get(end_vertex)*Vdistribution
            allocate(v1_bx(size(v1_nx)))
            allocate(v1_by(size(v1_nx)))
            call magneticField%interp%Evaluate(v1_nx, v1_ny, 1, 0, v1_bx)
            call magneticField%interp%Evaluate(v1_nx, v1_ny, 0, 1, v1_by)

            call v%x%Append(v1_nx)
            call v%y%Append(v1_ny)
            call v%psi%Append(v1_psi)
            call v%bx%Append(v1_bx)
            call v%by%Append(v1_by)
            v%ntot = v%ntot + nv

            deallocate(v1_bx)
            deallocate(v1_by)

            ! Make the new boundary faces
            new_bfcs = 0
            do i = 1, nc
                if (i == 1) then
                    call grid%GetFaceNumber(start_vertex, new_vertex(i), 3, face_num)
                else if (i == nc) then
                    call grid%GetFaceNumber(new_vertex(i-1), end_vertex, 3, face_num)
                else
                    call grid%GetFaceNumber(new_vertex(i-1), new_vertex(i), 3, face_num)
                end if
                new_bfcs(i) = face_num
            end do

            ! Give face label
            ar = (/ (lbl, i = 1, nc)/)
            call f%label%Set(new_bfcs,ar)

            ! Make new internal faces
            ! Get the internal vertixes 
            i_verts = 0
            do i = 1, nv

                ! Get common face
                c_fcs = GetCommonFace(c, cvs(i), cvs(i+1))

                vxs = [f%vert1%Get(c_fcs), f%vert2%Get(c_fcs)]

                if (vxs(1) /= end_vertex) then
                    i_verts(i) = vxs(1)
                else if (vxs(2) /= end_vertex) then
                    i_verts(i) = vxs(i)
                end if

            end do

            ! Connect

            new_ifcs = 0
            do i = 1, nv
                call grid%GetFaceNumber(new_vertex(i),i_verts(i), 3, face_num)
                new_ifcs(i) = face_num
            end do

            ! Make the cells
            do i = 1, nc

                ! Get aligned face
                fcst = GetCellFaceGA(c, cvs(i))
                ind = findloc(f%aligned%Get(fcst), 1, 1)
                fcst_al = fcst(ind)

                if (i == 1) then
                    
                    ! Bottom triangle
                    new_verts3 = 0
                    new_faces3 = 0

                    new_faces3(1) = fcst_al
                    new_faces3(2) = new_ifcs(i)
                    new_faces3(3) = new_bfcs(i)

                    new_verts3(1) = f%vert1%Get(fcst_al)
                    new_verts3(2) = f%vert2%Get(fcst_al)
                    new_verts3(3) = new_vertex(i)

                    new_faces = new_faces3
                    new_verts = new_verts3

                else if (i == nc) then

                    ! Last trapezoid
                    new_verts4 = 0
                    new_faces4 = 0

                    ! Get upper face
                    allocate(fcst_up(count(fcst /= fcst_al)))
                    fcst_up = pack(fcst, fcst /= fcst_al)

                    ! Test with cells of fcst_up
                    cvs1 = GetFaceCellGA(c, fcst_up(1), cvLookUp)
                    cvs2 = GetFaceCellGA(c, fcst_up(2), cvLookUp)

                    if (count(isMember(cvs1, cvs)) == 1) then
                        fcst_up_correct = fcst_up(1)
                    elseif (count(isMember(cvs2, cvs)) == 1) then
                        fcst_up_correct = fcst_up(2)
                    end if

                    new_faces4(1) = fcst_al
                    new_faces4(2) = fcst_up_correct
                    new_faces4(3) = new_bfcs(i)
                    new_faces4(4) = new_ifcs(i-1)

                    new_verts4(1) = f%vert1%Get(fcst_al)
                    new_verts4(2) = f%vert2%Get(fcst_al)
                    new_verts4(3) = f%vert1%Get(new_bfcs(i))
                    new_verts4(4) = f%vert2%Get(new_bfcs(i))    
                    
                    new_faces = new_faces4
                    new_verts = new_verts4

                    deallocate(fcst_up)
                    
                else

                    ! Inner trapezoid
                    new_verts4 = 0
                    new_faces4 = 0      

                    new_faces4(1) = fcst_al;
                    new_faces4(2) = new_ifcs(i)
                    new_faces4(3) = new_bfcs(i)
                    new_faces4(4) = new_ifcs(i-1)

                    new_verts4(1) = f%vert1%Get(fcst_al)
                    new_verts4(2) = f%vert2%Get(fcst_al)
                    new_verts4(3) = f%vert1%Get(new_bfcs(i))
                    new_verts4(4) = f%vert2%Get(new_bfcs(i))

                    new_faces = new_faces4
                    new_verts = new_verts4

                end if

                ! Add the cell
                call grid%AddCell(new_faces, new_verts, regs(i), ic)

                ! Adjust centroid
                call grid%CalcCentroidGA(c%ntot)

                ! Give correct cflag
                allocate(cells(1))
                cells = c%ntot
                call grid%DetermineCflags(cells)

                ! Housekeeping
                deallocate(cells)


            end do
            
            ! House keeping
            deallocate(i_verts)
            deallocate(new_ifcs)
            deallocate(new_bfcs)

            ! Remove faces - boundary faces and all non-aligned faces
            allocate(rem_faces(nc))
            rem_faces(1) = bfcs
            do i = 2, nc
                rem_faces(i) = GetCommonFace(c, cvs(i-1), cvs(i))
            end do

            rem_cells = cvs

            call grid%RemoveFaces(rem_faces)
            call grid%RemoveCells(rem_cells)

            ! Detection
            call grid%DetectStackedTrias(found, cvs, nc)

            ! Housekeeping
            deallocate(rem_faces)

        end do

        ! Recalculate the magneticfield
        if (.not. options%slab) call grid%RecalcMagn(magneticField)

        print *, 'Ended StackedToCutcell'

        ! Order - maybe not necessary - check TODO
        !cells = .true.
        !call grid%CheckVertOrder(is_ordered,cells)
        !call grid%ReorderCellConn(is_ordered)

        end associate

    end subroutine

    subroutine DetectStackedTrias(grid, found, cvs, counter)

        ! Description
        !============
        ! Detect stacked triangle. If a group of triangle goes from one boundary
        ! face to another only a part of the triangle group will be considered
        ! depending on the length ratio of the boundary face.

        ! Declare variables
        !==================
        class(GAGridUDT), intent(in)    :: grid
        logical, intent(out)            :: found
        integer(I8), intent(out)        :: counter
        integer(I8), allocatable, intent(out) :: cvs(:)

        ! Auxiliary
        integer(I8) :: i, v11, v12, v21, v22, nf, fcs, ifc, prev_cv, &
            prev_face
        integer(I8), allocatable, dimension(:) :: cvsD, bfaces, cv, fcsD, &
            fcsD2, fcs3, fcs3_nal, cvs2, cvLookUp, indfc
        real(R8) :: fcs_b1, fcs_b2, r 
        logical, allocatable :: log(:)


        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        ! Initialize
        found = .false.
        allocate(cvsD(c%ntot))
        counter = 0
        cvLookUp = GetCvLookUpGA(c)

        ! Get boundary faces
        indfc = (/ (i, i = 1, f%ntot)/)
        log = isBoundaryFaceGA(f, indfc)
        allocate(bfaces(count(log)))
        bfaces = pack(indfc,log)

        ! Loop over boundary faces
        do i = 1, size(bfaces)

            ! Get boundary face
            ifc = bfaces(i)
            fcs = bfaces(i)

            ! Get boundary cells
            cv = GetFaceCellGA(c, ifc, cvLookUp)

            ! Continue if triangle
            nf = c%faceP2%Get(cv(1))

            do while (nf == 3) 

                ! Store
                counter = counter + 1
                cvsD(counter) = cv(1)

                ! Boundary face is prev face
                prev_face = fcs
                prev_cv = cv(1)

                ! Take new poloidal faces
                fcsD = GetCellFaceGA(c, cv(1))
                log = ((fcsD /= prev_face) .and. (f%aligned%Get(fcsD) == 0))
                allocate(fcsD2(count(log)))
                fcsD2 = pack(fcsD, log)
                fcs = fcsD2(1)

                ! Get next cell
                if (.not.isBoundaryFaceGA(f, fcs)) then

                    ! Get the next face
                    cvs2 = GetFaceCellGA(c, fcs, cvLookUp)
                    cv = pack(cvs2, cvs2 /= prev_cv)

                    nf = c%faceP2%Get(cv(1))

                    ! Check to continue
                    if (nf == 3) then

                        ! Get faces
                        fcs3 = GetCellFaceGA(c, cv(1))
                        
                        ! Do not continue if any non-aligned face is a boundary face
                        fcs3_nal = pack(fcs3, f%aligned%Get(fcs3) == 0)
                        if (any(isBoundaryFaceGA(f, fcs3_nal)))  nf = 0

                    end if

                else

                    ! Exit
                    nf = 0

                end if

                ! Housekeeping
                deallocate(fcsD2)


            end do

            if (counter .gt. 1) then

                exit

            else

                ! Reinitialization
                cvsD = 0
                counter = 0 

            end if

        end do

        ! Post-process
        if (counter .gt. 1) then

            if (nf == 0) then

                ! Other boundary face reached
                ! Compute ratio of boundary faces
                v11 = f%vert1%Get(ifc)
                v12 = f%vert1%Get(ifc)
                fcs_b1 = sqrt( (v%x%Get(v11) - v%x%Get(v12))**2 + (v%y%Get(v11) - v%y%Get(v12))**2 )

                v21 = f%vert1%Get(fcs)
                v22 = f%vert1%Get(fcs)
                fcs_b2 = sqrt( (v%x%Get(v21) - v%x%Get(v22))**2 + (v%y%Get(v21) - v%y%Get(v22))**2 )

                r = fcs_b1 / (fcs_b1 + fcs_b2)

                ! Rescale counter
                counter = nint(counter*r)

            end if

            ! Flag
            found = .true.

        end if

        ! Trim
        cvs = cvsD(1:counter)

        end associate

    end subroutine

    ! Merging
    !========
    subroutine DoMerging(grid, magneticField, qm, options)

        ! Description
        !============
        ! Wrapper function for merging operations

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT)        :: grid
        type(MagneticFieldUDT)  :: magneticField
        type(QualityMetricUDT)  :: qm
        type(GAoptionsUDT)      :: options

        ! Auxiliary
        integer(I8) :: i

        ! Printing
        print *, 'Merging: ', options%merge_crit


        ! Calculate metrics
        call qm%CalculateQualityMetrics(grid, options, magneticField,.false.,.true.)

        i = 0
        do while ((qm%merge_fc /= 0) .and. (i .lt. options%n_merge))

            ! Merge
            call grid%MergeCells(qm, options)

            ! Update counter
            i = i + 1

            ! Printing
            print *, 'Applied cell merging ', i

            ! Recalculate magneticField
            if (.not.options%slab) call grid%RecalcMagn(magneticField)

            ! Recompute metrics
            call qm%CalculateQualityMetrics(grid, options, magneticField,.false.,.true.)  
            
            ! Check grid
            if (options%debug) call grid%CheckUnstructuredGrid(.false.)

        end do

        ! Transform remaining pentagons into triangles - TODO

        ! Printing
        print *, 'Ended merging: ', options%merge_crit



    end subroutine

    subroutine MergeCells(grid, qm, options)

        ! Description
        !============
        ! Wrapper for merging routine which can merge two cells together.

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)         :: grid
        type(QualityMetricUDT), intent(inout)   :: qm
        type(GAoptionsUDT)                      :: options

        ! Auxiliary
        integer(I8) :: i
        integer(I8), allocatable :: empty_surf(:), indfs(:), cvLookUp(:), &
            cells(:)
        logical, allocatable :: b_flag(:), log(:)
       

        ! Get separatrix
        cvLookUp = GetCvLookUpGA(grid%cell)
        call grid%GiveSeparatrices(.true., .true., .false., cvLookUp)

        ! Get merging case and do the merge for face fc
        call grid%OneMerge(qm%merge_fc,.true.,.false.,.false.)

        ! If hexagonal and pentagonal cells are not allowed
        if (options%no_pents  .or. options%no_hex) then
            call grid%MergeRec(options)
        end if

        ! Determine cflags
        b_flag = (grid%face%label%Get() /= 0)
        cells = (/ (i, i = 1, grid%cell%ntot )/)
        call grid%DetermineCflags(cells, b_flag)

        ! Remove empty flux surfaces
        log = (grid%data%fluxdata%fluxsurfacefacesP2%Get() == 0)
        allocate(empty_surf(count(log)))
        indfs = (/ (i, i= 1, grid%data%fluxdata%fluxsurfacefacesP2%Size() )/)
        empty_surf = pack(indfs, log)

        call grid%data%fluxdata%fluxsurfacefacesP1%Remove(empty_surf)
        call grid%data%fluxdata%fluxsurfacefacesP2%Remove(empty_surf)

        ! Update number of flux surfaces
        grid%data%fluxdata%nFs = grid%data%fluxdata%fluxsurfacefacesP1%Size()

        ! Check uniqueness
        if (options%debug) call CheckUniqueness(grid%data%fluxdata%fluxsurfacefaces)

        ! Recover flux surface vertices
        call grid%GetFsVxFromFsFc(options)

    end subroutine

    subroutine OneMerge(grid, fc, starter,pent_to_tria,special_case)

        ! Description
        !============
        ! Performs one merging operation. The method is base on the grid 
        ! surrounding the selected face

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout) :: grid
        integer(I8), intent(in)         :: fc
        logical                         :: starter, pent_to_tria, special_case

        ! Auxiliary
        integer(I8), allocatable :: cvs(:)
        character(:), allocatable :: caseID


        ! Determine the merge case
        call grid%DetermineMergeCaseID(fc, starter, pent_to_tria, special_case, caseID, cvs)

        ! Select the case
        select case (caseID)

        case ('4444')

            ! This is the case in a internal regular grid. Probably this is the starting point
            call grid%Merge4444(fc, cvs)

        case ('4443')

            ! This case is most often called in a grid with quads, so optimize this
            call grid%Merge4443(fc, cvs)

        case ('4433')

            ! Merge two quads to a quad at a boundary
            call grid%Merge4433(fc, cvs)

        case ('3443')

            ! Merging a quad and triangle at a boundary
            call grid%Merge3443(fc, cvs, starter)

        case ('3433')
        case ('334')

            ! Merging to triangles 
            call grid%Merge334(fc, cvs)

        case ('333')
        case ('53')
        case ('53B')
        case ('4433B1')
        case ('5T3')
        case ('5spec')
        case ('99')

            print *, 'Merge can not be done because not allow, try to smooth first'

        case default

            call gdErrorHandler('OneMerge: caseID not implemented')

        end select
    end subroutine
    
    recursive subroutine MergeRec(grid, options)

        ! Description
        !============
        ! Merges cells recursively to avoid pentagonal and hexagonal cells.
        ! Free verts (= hanging nodes) are vertices from pentagonal or
        ! hexagonal cells which should be removed to obtain only
        ! quadrilateral or triangular cells.

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout) :: grid
        type(GAoptionsUDT), intent(in)  :: options

        ! Auxiliary
        logical :: pent_to_tria, special_case 
        integer(I8) :: i, j, k, lim, fc
        integer(I8), allocatable, dimension(:) :: cellsD2, cellsD, cells, indCv, &
            verts, cvs, fcs, bfcs, cvP, cvs_b, vs1, vs2, cvLookUp

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        ! Initialize
        allocate(indCv(c%ntot))
        indCv = (/ (i, i = 1, c%ntot)/)
        pent_to_tria = .false.
        special_case = .false.
        cvLookUp = GetCvLookUpGA(c)        

        ! Get cells to merge
        if (options%no_pents) then
            if (options%no_pents_area_merge) then
                call grid%AreaConstraintPents(options,options%dist_function_threshold_merge,cellsD)
            else
                cellsD = (/ (i, i = 1, c%ntot )/)
            end if
            allocate(cellsD2(count(c%faceP2%Get(cellsD) .ge. 5)))
            cellsD2 = pack(cellsD, c%faceP2%Get(cellsD) .ge. 5)
            lim = 4

        else if (options%no_hex) then

            allocate(cellsD2(count(c%cflags%Get() .gt. 5)))
            cellsD2 = pack(indCv, c%cflags%Get() .gt. 5)
            lim = 5

        end if

        ! Also ad tria4 cases
        allocate(cells(size(cellsD2) + count(c%cflags%Get() == 4)))
        cells = [cellsD2, pack(indCv, c%cflags%Get() == 4)]
        deallocate(cellsD2)

        ! Get the correct face to merge cells by finding handing nodes
        if (size(cells) /= 0) then

            ! Loop over chosen cells
            do k = 1, size(cells)

                ! Get vertices
                verts = GetCellVertGA(c, cells(k))

                ! Find a merge face
                fc = 0
                do i = 1, size(verts)

                    ! Check faces of vertices
                    fcs = GetVertFaceGA(f, verts(i))

                    ! Hanging nodes inside the mesh have three faces or less
                    if ((size(fcs) .le. 3) .and. .not.isBoundaryVertGA(grid, verts(i))) then

                        ! Get the correct merge face
                        do j = 1, size(fcs)

                            cvs = GetFaceCellGA(c, fcs(j),cvLookUp)
                            if (.not.any(cells(k) == cvs)) then

                                fc = fcs(j)
                                exit

                            end if

                        end do

                    end if
                end do

                ! If no merge face was found yet
                ! Allow boundary vertices for special case where hanging node of
                ! pentagon lays on the boundary  
                if (fc == 0) then

                    fcs = GetCellFaceGA(c, cells(k))
                    if (size(fcs) == 5 .and. count(isBoundaryFaceGA(f, fcs)) == 2) then

                        pent_to_tria = .true.
                        allocate(bfcs(2))
                        bfcs = pack(fcs, isBoundaryFaceGA(f, fcs))
                        fc = bfcs(1)

                        deallocate(bfcs)
                    else

                        ! Loop over vertices of cell
                        do i = 1, size(verts)

                            fcs = GetVertFaceGA(f, verts(i))

                            if (size(fcs) .lt. 3) then

                                do j = 1, size(fcs)

                                    cvs = GetFaceCellGA(c, fcs(j), cvLookUp)
                                    if (size(cvs) == 2) then

                                        cvP = c%faceP2%Get(cvs)

                                        if (maxval(cvP) == 5) then
                                            if ((cvP(1)==3 .and. cvP(2)==5) &
                                                .or.(cvP(1)==5 .and. cvP(1)==3)) then
                                                
                                                fc = fcs(j)
                                                exit
                                            end if
                                        else

                                            fc = fcs(j)
                                            exit

                                        end if

                                    end if

                                end do

                            end if

                        end do

                    end if



                end if

                ! Special case where two following cells to merge only have one
                ! vertex in common  
                if (fc == 0) then

                    cvs = GetCellNeigsGA(grid, cells(k), cvLookUp)

                    allocate(cvs_b(count(isBoundaryCellGA(grid, cvs))))
                    cvs_b = pack(cvs, isBoundaryCellGA(grid, cvs))

                    if (size(cvs_b) == 2) then

                        ! Check whether they only have one common vert
                        vs1 = GetCellVertGA(c, cvs_b(1))
                        vs2 = GetCellVertGA(c, cvs_b(2))

                        if (count(isMember(vs1, vs2)) == 1) then

                            special_case = .true.

                            ! Get common face
                            fc = GetCommonFace(c,cells(k), cvs_b(1))

                        end if
                    else

                        call gdErrorHandler('MergeRec: case not yet implemented')

                    end if

                end if

                if ((fc /= 0) .and. .not.special_case) then

                    ! Avoid merging a pentagon with a quad
                    cvs = GetFaceCellGA(c, fc, cvLookUp)
                    if (size(cvs) == 2) then
                        
                        cvP = c%vertP2%Get(cvs)
                        if (((cvP(1) == 5) .and. (cvP(2) == 4)) &
                            .or. ((cvP(1)==4) .and. (cvP(2)==5))) then

                            ! Set an indicator to not merge this cell
                            fc = 0                                   
                            call c%cflags%Set(cells(k), 7)

                        end if

                    end if

                end if

                if (fc /= 0) exit

            end do

            ! Do a one merge with fc
            if (fc /= 0) then

                call grid%OneMerge(fc, .false., pent_to_tria, special_case)

            end if

        end if

        ! Housekeeping
        deallocate(cells)
        
        ! Recursive call
        ! Determine cells
        if (options%no_pents_area_merge) then

            call grid%AreaConstraintPents(options,options%dist_function_threshold_merge,cellsD)

            ! Detect pents in the no_pents_area which are not mergeable
            allocate(cells(count(c%cflags%Get(cellsD) /= 7 )))
            cells = pack(cellsD, c%cflags%Get(cellsD) /= 7 )
        else

            cells = (/ (i, i = 1, c%ntot)/)

        end if 

        if ((maxval(c%faceP2%Get(cells)) .gt. lim) &
            .or. (any(c%cflags%Get() == 4)) ) then

                call grid%MergeRec(options)

        end if

        end associate
    end subroutine

    subroutine DetermineMergeCaseID(grid, fc, starter, pent_to_tria, special_case, caseID, cvs)

        ! Description
        !============
        ! Determines merge case

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)         :: grid
        integer(I8), intent(in)                 :: fc
        logical, intent(in)                     :: starter, pent_to_tria, special_case
        character(:), allocatable, intent(out)  :: caseID
        integer(I8), allocatable, intent(out)   :: cvs(:)

        ! Auxiliary
        integer(I8) :: i, j, nfc1, nfc2, faces_sepD(grid%face%ntot), &
            s, nf, counter, vB, ncvs, cvB, vs_fcB(2), tria
        integer(I8), allocatable :: cvLookUp(:), fc1(:), fc2(:), loc(:), &
            faces_sep(:), vc(:), tf(:), tfB(:), vertsT(:), cvtypes(:)
        logical, allocatable :: v1(:), b1(:)

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            fd => grid%data%fluxdata &
            )
        
        ! Initialize
        cvLookUp = GetCvLookUpGA(c)

        ! Special cases
        if (pent_to_tria) then

            caseID = '5T3'
            cvs = GetFaceCellGA(c, fc, cvLookUp)
            return

        end if

        if (special_case) then
            caseID = '5spec'
            cvs = GetFaceCellGA(c, fc, cvLookUp)
        end if

        ! More regular cases
        !===================
        ! Get information on surroundings
        allocate(vc(2))
        vc(1) = f%vert1%Get(fc)
        vc(2) = f%vert2%Get(fc)

        fc1 = GetVertFaceGA(f, vc(1))
        nfc1 = size(fc1)
        if (nfc1 .gt. 4) nfc1 = 4   ! Clipping to 4 because more than 4 doesn't change the situation

        fc2 = GetVertFaceGA(f, vc(2))
        nfc2 = size(fc2)
        if (nfc2 .gt. 4) nfc2 = 4

        cvs = GetFaceCellGA(c, fc, cvLookUp)
        ncvs = size(cvs)

        ! Check for boundary face and other condition where merging is not possible
        if (ncvs == 1) then

            caseID = '99'
            print *, 'Merging not possible because face is boundary face'
            return

        end if

        if (c%reg%Get(cvs(1)) /= c%reg%Get(cvs(2))) then

            caseID = '99'
            print *, 'Cannot merge because neighboring cells are from' // &
             & 'different regions. This would result in issues with the X-point.' // &
             & 'This should be detect earlier in the code!'

            return 

        end if

        ! Check for face on separatrix
        faces_sepD = 0
        counter = 0
        do i = 1, grid%data%nsep
            nf = fd%fluxsurfacefacesP2%Get(grid%data%sepID(i))
            s = fd%fluxsurfacefacesP1%Get(grid%data%sepID(i))
            loc = (/ (j, j = s, s+nf-1)/)
            faces_sepD(counter+1:counter+nf) = fd%fluxsurfacefaces%Get(loc)
            counter = counter + nf
        end do

        ! Trim
        faces_sep = faces_sepD(1:counter)

        if (any(faces_sep == fc)) then

            caseID = '99'
            print *, 'Merging not possible as face is on flux surface of separatrix'
        end if

        ! Determine actual case
        allocate(cvtypes(size(cvs)))
        cvtypes = c%faceP2%Get(cvs)

        if (cvtypes(1) == 4 .and. cvtypes(2) == 4) then
        
            if ( nfc1 == 4 .and. nfc2 == 4) then

                caseID = '4444'

            elseif ( (nfc1 == 3 .and. nfc2 == 4) &
                   .or. (nfc1 == 4 .and. nfc2 == 3)) then
                ! For special case where on quad is a boundary cell and one of the
                ! vertices of the merge face is also a boundary cell
                caseID = '4443'
                v1 = isBoundaryVertGA(grid, vc)
                if (count(v1) == 1) then

                    ! Get boundary vertex
                    if (v1(1)) then
                        vB = vc(1)
                    else if (v1(2)) then
                        vB = vc(2)
                    end if

                    b1 = isBoundaryCellGA(grid, cvs)
                    if (count(b1) == 1) then

                        ! Get boundary cell
                        if (b1(1)) then
                            cvB = cvs(1)
                        else if (v1(2)) then
                            cvB = cvs(2)
                        end if       
                        tf = GetCellFaceGA(c, cvB)   
                        allocate(tfB(count( isBoundaryFaceGA(f, tf)) ))
                        tfB = pack(tf, isBoundaryFaceGA(f, tf))

                        vs_fcB(1) = f%vert1%Get(tfB(1))
                        vs_fcB(2) = f%vert2%Get(tfB(1))

                        if (any(vB == vs_fcB)) caseID = '443B1'

                        if (size(tfB) .gt. 1) call gdErrorHandler('DetermineMergeCaseID:extra implementation needed')

                    end if

                end if
                


            elseif ( nfc1 == 3 .and. nfc2 == 3) then 

                caseID = '4433'

            else

                call gdErrorHandler('DetermineMergeCaseID: undefined case')

            end if

        else if ((cvtypes(1) == 3 .and. cvtypes(2) == 4) &
                .or. (cvtypes(1) == 4 .and. cvtypes(2) == 3)) then

            if (nfc1 == 4 .and. nfc2 == 4) then
                caseID = '4444'
            else if ((nfc1 == 4 .and. nfc2 == 3) .or. (nfc1 == 3 .and. nfc2 == 4)) then
                caseID = '3443'
            else if (nfc1 == 3 .and. nfc2 == 3) then
                caseID = '3433'
            end if

        else if (cvtypes(1) == 3 .and. cvtypes(2) == 3) then

            if (min(nfc1, nfc2) == 3 .and. starter) then

                caseID = '334'

            elseif (min(nfc1, nfc2) == 3 .and. starter) then

                caseID = '333'

            elseif (min(nfc1, nfc2) == 4) then

                caseID = '334'

            end if

        else if ((cvtypes(1) == 5 .and. cvtypes(2) == 3) &
                .or. (cvtypes(1) == 3 .and. cvtypes(2) == 5)) then

            caseID = '53'
            if (cvtypes(1) == 3) then
                tria = cvs(1)
            else if (cvtypes(2) == 3) then
                tria = cvs(2)
            end if

            vertsT = GetCellVertGA(c, tria)
            if (count(isBoundaryVertGA(grid, vertsT)) == 3) caseID = '53B'

        else if (maxval(cvtypes) .gt. 4) then

            caseID = '99'
            call gdErrorHandler('DetermineMergeCaseID: merging pents or hex not implemented')

        end if

        end associate
 
    end subroutine

    subroutine Merge4444(grid, fc, cvs)

        ! Description
        !============
        ! Merges two quads to a hex in starting position

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT)    :: grid
        integer(I8)         :: fc, cvs(2)

        ! Auxiliary
        integer(I8) :: ic
        integer(I8), allocatable, dimension(:) :: vx1, vx2, fc1, fc2, new_verts, new_facesD, &
            new_faces, fc_rem, cv_rem, vertsD, facesD

        ! Make hexagon
        vx1 = GetCellVertGA(grid%cell, cvs(1))
        vx2 = GetCellVertGA(grid%cell, cvs(2))
        fc1 = GetCellFaceGA(grid%cell, cvs(1))
        fc2 = GetCellFaceGA(grid%cell, cvs(2))

        ! Verts
        vertsD = [vx1, vx2]
        call Unique(vertsD, new_verts)

        ! Faces
        facesD = [fc1, fc2]
        call Unique(facesD, new_facesD)
        allocate(new_faces(count(new_facesD /= fc)))
        new_faces = pack(new_facesD, new_facesD /= fc)

        ! Add new cell
        call grid%AddCell(new_faces, new_verts, grid%cell%reg%Get(cvs(1)), ic)

        ! Remove face
        fc_rem = fc
        call grid%RemoveFaces(fc_rem)

        ! Remove cells
        cv_rem = cvs
        call grid%RemoveCells(cv_rem)

    end subroutine

    subroutine Merge3443(grid, fc, cvs, starter)

        ! Description
        !============
        ! Merging a quad and a triangle to a pent as one of the
        ! common face vertices has one three faces attached

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout) :: grid
        integer(I8)                     :: fc, cvs(2)
        logical                         :: starter

        ! Auxiliary
        integer(I8) :: i, ic, vxF(2), nv1cvs, nv2cvs, three_vert, f1n
        integer(I8), allocatable, dimension(:) :: v1cvs, v2cvs, fc23, &
            cells_rem, faces_rem, verts_ic, fcs_v, verts_rem        

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face)

        ! Initialize
        vxF = [f%vert1%Get(fc),  f%vert2%Get(fc)]

        ! Three vert
        v1cvs = GetVertCellGA(c, vxF(1))
        v2cvs = GetVertCellGA(c, vxF(2))
        nv1cvs = size(v1cvs)
        nv2cvs = size(v2cvs)

        ! Determine three vert
        if (starter) then
            if (nv1cvs .lt. 4) then
                three_vert = vxF(1)
            else if (nv2cvs .lt. 4) then
                three_vert = vxF(2)
            end if
        else
            if ((nv1cvs .lt. 4) .and. .not.isBoundaryVertGA(grid, vxF(1))) then
                three_vert = vxF(1)
            else if ((nv2cvs .lt. 4) .and. .not.isBoundaryVertGA(grid, vxF(2))) then
                three_vert = vxF(2)
            end if
        end if

        ! Make merge cell
        !================
        call grid%MakeMergePent(three_vert, fc, cvs, f1n, fc23, ic)

        ! Adapt the neighbor, if there is any
        !====================================
        call grid%AdaptNeigThreeVert(three_vert, cvs, fc23, f1n, cells_rem)

        ! Remove fc
        faces_rem = [fc, fc23]
        call grid%RemoveFaces(faces_rem)

        ! Remove cells
        call grid%RemoveCells(cells_rem)
        
        ! Mark the new cells as tria4 if necessary
        ic = ic - size(cells_rem)
        verts_ic = GetCellVertGA(c, ic)
        do i = 1, size(verts_ic)
            fcs_v = GetVertFaceGA(f, verts_ic(i))
            if (size(fcs_v) == 3 .and. .not.isBoundaryVertGA(grid, verts_ic(i))) then
                call c%cflags%Set(ic, 4)
            end if
        end do

        ! Remove vertex
        allocate(verts_rem(1))
        verts_rem = three_vert
        call grid%RemoveVertices(verts_rem)

        end associate

    end subroutine

    subroutine Merge4443(grid, fc, cvs)

        ! Description
        !============
        ! Merges two quads to a pent as one of  
        ! the commonface vertices has only three faces attached

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT) :: grid
        integer(I8)      :: fc, cvs(2)

        ! Auxiliary
        integer(I8) :: three_vert, vxF(2), nv1cvs, nv2cvs, f1n, ic 
        integer(I8), allocatable, dimension(:) :: v1cvs, v2cvs, fc23, faces_rem, verts_rem, &
            cells_rem


        associate(&
            c => grid%cell, &
            f => grid%face &
            )
            
        ! Verts
        vxF = [f%vert1%Get(fc), f%vert2%Get(fc)]
        ! Three vert
        v1cvs = GetVertCellGA(c, vxF(1))
        v2cvs = GetVertCellGA(c, vxF(2))
        nv1cvs = size(v1cvs)
        nv2cvs = size(v2cvs)

        ! Determine three vert
        if (nv1cvs .lt. 4) then
            three_vert = vxF(1)
        else if (nv2cvs .lt. 4) then
            three_vert = vxF(2)
        end if

        ! Make new three face and build new pentagon cell
        call grid%MakeMergePent(three_vert, fc, cvs, f1n, fc23, ic)

        ! Adapt neighbor of three_vert
        call grid%AdaptNeigThreeVert(three_vert, cvs, fc23, f1n, cells_rem)

        ! Remove fc
        faces_rem = [fc, fc23]
        call grid%RemoveFaces(faces_rem)

        ! Remove cells
        call grid%RemoveCells(cells_rem)
        
        ! Remove vertex
        allocate(verts_rem(1))
        verts_rem = three_vert
        call grid%RemoveVertices(verts_rem)

        end associate

    end subroutine

    subroutine Merge4433(grid, fc, cvs)

        ! Description
        !============
        ! Merges two quads to a quad at a boundary 

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT) :: grid
        integer(I8)      :: fc, cvs(2)

        ! Auxiliary
        integer(I8) :: i, ic, vxF(2), nv1cvs, nv2cvs, three_vert, f1n, f2n, counter, b_vert
        integer(I8), allocatable, dimension(:) :: v1cvs, v2cvs, fc3, fc23, fc3b, &
            fc23b , fc1, fc2, new_verts, new_vertsUD, new_vertsD, new_vertsD2, fc_rem, &
            new_facesUD, new_facesL, new_facesD, fc_rem2, cells_rem1, cells_rem2, cells_remD, &
            cells_rem, vx_rem, new_faces, vx1, vx2

        associate(&
            c => grid%cell, &
            f => grid%face &
            )
            
        ! Verts
        vxF = [f%vert1%Get(fc), f%vert2%Get(fc)]
        ! Three vert
        v1cvs = GetVertCellGA(c, vxF(1))
        v2cvs = GetVertCellGA(c, vxF(2))
        nv1cvs = size(v1cvs)
        nv2cvs = size(v2cvs)

        ! Determine three vert
        if ((nv1cvs == 3).and.(nv2cvs == 2)) then
            three_vert = vxF(1)
            b_vert = vxF(2)
        else if ((nv2cvs == 3) .and. (nv1cvs == 2)) then
            three_vert = vxF(2)
            b_vert = vxF(1)
        else ! does not matter, do check later
            three_vert = vxF(1)
            b_vert = vxF(2)
        end if

        ! Make new face a three_vert
        call grid%MakeNewThreeFace(three_vert, fc, fc3, fc23, f1n)
        call grid%MakeNewThreeFace(b_vert, fc, fc3b, fc23b, f2n)

        ! Make merge cell at the boundary
        vx1 = GetCellVertGA(grid%cell, cvs(1))
        vx2 = GetCellVertGA(grid%cell, cvs(2))

        fc1 = GetCellFaceGA(grid%cell, cvs(1))
        fc2 = GetCellFaceGA(grid%cell, cvs(2))

        ! Verts
        new_vertsD = [vx1, vx2]
        call Unique(new_vertsD, new_vertsUD)
        allocate(new_vertsD2(count(new_vertsUD /= three_vert)))
        new_vertsD2 = pack(new_vertsUD,new_vertsUD /= three_vert)

        allocate(new_verts(count(new_vertsD2 /= b_vert)))
        new_verts = pack(new_vertsD2, new_vertsD2 /= b_vert)

        ! Faces
        new_facesD = [fc1, fc2, f1n, f2n]
        call Unique(new_facesD, new_facesUD)
        
        allocate(new_facesL(size(new_facesUD)))
        new_facesL = 0
        counter = 0
        fc_rem  = [fc3, fc23b]
        do i = 1, size(new_facesUD)
            if (.not.any(new_facesUD(i) == fc_rem)) then
                counter = counter + 1
                new_facesL(counter) = new_facesUD(i)
            end if
        end do

        ! Trim 
        new_faces = new_facesL(1:counter)

        ! Add new cell
        call grid%AddCell(new_faces, new_verts, grid%cell%reg%Get(cvs(1)), ic) 
        call grid%cell%cflags%Set(ic, maxval(grid%cell%cflags%Get(cvs)))

        ! Adapt neighboring cell
        call grid%AdaptNeigThreeVert(three_vert, cvs, fc23, f1n, cells_rem1)
        call grid%AdaptNeigThreeVert(b_vert, cvs, fc23b, f2n, cells_rem2)

        ! Remove fc 
        fc_rem2 = [fc_rem, fc23]
        call grid%RemoveFaces(fc_rem2)

        ! Remove cells
        cells_remD = [cells_rem1, cells_rem2]
        call Unique(cells_remD, cells_rem)
        call grid%RemoveCells(cells_rem)

        ! Remove vertices
        vx_rem = vxF
        call grid%RemoveVertices(vx_rem)

        end associate

    end subroutine

    subroutine Merge334(grid, fc, cvs)

        ! Description
        !============
        ! Merges two triangles without merged cells around

        ! Declare variables
        !==================
        class(GAGridUDT) :: grid
        integer(I8) :: fc, cvs(2)

        ! Auxiliary
        integer(I8) :: ic, n_al
        integer(I8), allocatable, dimension(:) :: vx1, vx2, fc1, fc2, vertsD, new_verts, facesD, &
            new_faces, fc_rem, cell_rem, fcs_al

        ! Make merge cell
        vx1 = GetCellVertGA(grid%cell, cvs(1))
        vx2 = GetCellVertGA(grid%cell, cvs(2))

        fc1 = GetCellFaceGA(grid%cell, cvs(1))
        fc2 = GetCellFaceGA(grid%cell, cvs(2))

        ! Verts
        vertsD = [vx1, vx2]
        call Unique(vertsD, new_verts)

        ! Faces
        facesD = [fc1, fc2]
        allocate(new_faces(count(facesD /= fc)))
        new_faces = pack(facesD, facesD /= fc)

        ! Add new cell
        call grid%AddCell(new_faces, new_verts, grid%cell%reg%Get(cvs(1)), ic)

        ! Determine cflag
        n_al = count(grid%face%aligned%Get(new_faces) == 1)
        if (n_al == 2) then
            allocate(fcs_al(n_al))
            fcs_al = pack(new_faces, grid%face%aligned%Get(new_faces) == 1)
            if (HaveCommonVert(grid%face, fcs_al(1),fcs_al(2)))  call grid%cell%cflags%Set(ic, 4) 
            deallocate(fcs_al)
        else 
            call grid%cell%cflags%Set(ic, 4)
        end if

        ! Remove fc
        fc_rem = fc
        call grid%RemoveFaces(fc_rem)

        ! Remove cells
        cell_rem = cvs
        call grid%RemoveCells(cell_rem)

    end subroutine 

    subroutine MakeMergePent(grid,three_vert, fc, cvs, f1n, fc23, ic)

        ! Description
        !============
        ! Commonly used piece of code to make a new face at across faces connected by a 
        ! hanging node or vertex with only three faces attached.

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT)                        :: grid
        integer(I8), intent(in)                 :: three_vert, fc, cvs(2)
        integer(I8), intent(out)                :: f1n, ic 
        integer(I8), allocatable, intent(out)   :: fc23(:)

        ! Auxiliary
        integer(I8) :: counter, i
        integer(I8), allocatable, dimension(:) :: fc3, vx1, vx2, &
            fc1, fc2, new_vertsD, new_vertsUD, new_verts, new_facesD, new_facesUD, &
            new_facesL, new_faces

        ! Make new three_face
        call grid%MakeNewThreeFace(three_vert, fc, fc3, fc23, f1n)

        ! Make new cell with hanging node
        vx1 = GetCellVertGA(grid%cell, cvs(1))
        vx2 = GetCellVertGA(grid%cell, cvs(2))

        fc1 = GetCellFaceGA(grid%cell, cvs(1))
        fc2 = GetCellFaceGA(grid%cell, cvs(2))

        ! Verts
        new_vertsD = [vx1, vx2]
        call Unique(new_vertsD, new_vertsUD)
        allocate(new_verts(count(new_vertsUD /= three_vert)))
        new_verts = pack(new_vertsUD,new_vertsUD /= three_vert)

        ! Faces
        new_facesD = [fc1, fc2, f1n]
        call Unique(new_facesD, new_facesUD)
        
        allocate(new_facesL(size(new_facesUD)))
        new_facesL = 0
        counter = 0
        do i = 1, size(new_facesUD)
            if (.not.any(new_facesUD(i) == fc3)) then
                counter = counter + 1
                new_facesL(counter) = new_facesUD(i)
            end if
        end do

        ! Trim 
        new_faces = new_facesL(1:counter)

        ! Add new cell
        call grid%AddCell(new_faces, new_verts, grid%cell%reg%Get(cvs(1)), ic)        

        
    end subroutine

    subroutine MakeNewThreeFace(grid, three_vert, fc, fc3, fc23, f1n)

        ! Description
        !============
        ! Make a new face from two faces connect with hanging node

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT)                        :: grid
        integer(I8), intent(in)                 :: three_vert, fc
        integer(I8), intent(out)                :: f1n
        integer(I8), allocatable, intent(out)   :: fc3(:), fc23(:)

        ! Auxiliary
        integer(I8), allocatable, dimension(:) :: vf1n, vc3

        fc3 = GetVertFaceGA(grid%face, three_vert)
        allocate(fc23(count(fc3 /= fc)))
        fc23 = pack(fc3, fc3 /= fc)

        allocate(vc3(size(fc23)*2))
        vc3 = [grid%face%vert1%Get(fc23(1)), grid%face%vert1%Get(fc23(2)), &
               grid%face%vert2%Get(fc23(1)), grid%face%vert2%Get(fc23(2)) ]
        allocate(vf1n(count(vc3 /= three_vert)))
        vf1n = pack(vc3, vc3 /= three_vert)

        ! Create new face
        call grid%GetFaceNumber(vf1n(1), vf1n(2), 3, f1n) ! Fast method
        call grid%face%label%Set(f1n,grid%face%label%Get(fc23(1))) 
        if (grid%face%label%Get(fc23(1)) /= grid%face%label%Get(fc23(2))) &
            print *, 'MakeNewThreeFace: face%label of merged faces were not equal!'
        if ( grid%face%aligned%Get(fc23(1)) == 1 .and. grid%face%aligned%Get(fc23(2)) == 1 ) &
            call grid%face%aligned%Set(f1n,1)
        call grid%AddFaceToFsFc(f1n, fc23)

    end subroutine

    subroutine AdaptNeigThreeVert(grid, three_vert, cvs, fc23, f1n, cells_rem)

        ! Description
        !============
        ! Commonly used piece of code to adapt the neighbor to the merged cells

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT)                        :: grid
        integer(I8), intent(in)                 :: three_vert, f1n, cvs(2)
        integer(I8), intent(in)                 :: fc23(:)
        integer(I8), allocatable, intent(out)   :: cells_rem(:)

        ! Auxiliary
        integer(I8) :: counter, ic, i
        integer(I8), allocatable, dimension(:) :: cells_n, neigD, verts_n, new_vertsN, &
            faces_n, faces_nD, new_faces_n

        ! Get neighbors in radial direction
        cells_n = GetVertCellGA(grid%cell, three_vert)

        ! Eliminate cvs
        allocate(neigD(size(cells_n)))
        neigD = 0
        counter = 0
        do i = 1, size(cells_n)
            if ((cells_n(i) /= cvs(1)) .and. (cells_n(i) /= cvs(2))) then
                counter = counter + 1
                neigD(counter) = cells_n(i)
            end if
        end do

        if (counter == 1) then

            ! Verts
            verts_n = GetCellVertGA(grid%cell, neigD(1))
            allocate(new_vertsN(count(verts_n /= three_vert)))
            new_vertsN = pack(verts_n, verts_n /= three_vert)

            ! Face 
            faces_n = GetCellFaceGA(grid%cell, neigD(1))

            ! Eliminate fc23
            allocate(faces_nD(size(faces_n)))
            faces_nD = 0
            counter = 0
            do i = 1, size(faces_n)
                if (faces_n(i) /= fc23(1) .and. faces_n(i) /= fc23(2)) then
                    counter = counter + 1
                    faces_nD(counter) = faces_n(i)
                end if
            end do

            allocate(new_faces_n(counter+1))
            new_faces_n = [faces_nD(1:counter), f1n]

            ! Add new cell
            call grid%AddCell(new_faces_n, new_vertsN, grid%cell%reg%Get(neigD(1)), ic)

            ! To remove
            allocate(cells_rem(3))
            cells_rem = [cvs, neigD(1)]


        elseif (counter == 0) then

            if (grid%face%label%Get(fc23(1)) /= grid%face%label%Get(fc23(2))) then
                print *, 'AdaptNeigThreeVert: face%label of merged faces were not equal!'
            end if 
            cells_rem = cvs
        else if (counter .gt. 1) then

            call gdErrorHandler('AdaptNeigThreeVert: case not implemented')

        end if


    end subroutine

    ! Splitting
    !==========
    subroutine DoSplitting(grid, magneticField, qm, options)

        ! Description
        !============
        ! Wrapper function for splitting operation

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)         :: grid
        type(MagneticFieldUDT), intent(in)      :: magneticField
        type(QualityMetricUDT), intent(inout)   :: qm
        type(GAoptionsUDT), intent(inout)       :: options

        ! Auxiliary
        integer(I8) :: j
        character(:), allocatable :: t1, t2, t3
        logical :: rad_splitting


        ! Printing
        select case (options%splittype)
        case ('rad')
            t1 = 'radial'
            t2 = options%rad_type
            t3 = 'Applied radial splitting'
            rad_splitting = .true.
        case ('pol')
            t1 = 'poloidal'
            t2 = options%pol_type
            t3 = 'Applied poloidal splitting'
            rad_splitting  = .false.
        end select
        print *,  'Splitting: ', t1, ' ', t2
        
        if (options%slab) call grid%vert%psi%Set(grid%vert%x%Get())

        ! Calculate metric and find splitting cell
        call qm%CalculateQualityMetrics(grid, options, magneticField, .true., .false.)
        j = 0

        ! While a splitting cell is found
        do while (qm%split_cv /= 0 .and. j .lt. options%n_split)

            call grid%Splitting(qm, options, magneticField)

            if (.not.options%slab) call grid%RecalcMagn(magneticField)

            ! Update counter
            j = j + 1

            ! Printing 
            print *, t3, ' ', j

            ! Check grid
            if (options%debug) call grid%CheckUnstructuredGrid(.false.)

            if (rad_splitting) then

                ! Remove small triangle which were possibly created
                ! Calculate quality metric
                call qm%CalculateQualityMetrics(grid, options, magneticField, .false., .false.)

                ! Remove Small triangles
                if (options%rem_small_trias) then

                    print *, 'Removing small triangles'
                    call grid%RemoveSmallTriangle(magneticField, qm, options)
                    print *, 'Ended removing small triangles'

                    if (options%debug) call grid%CheckUnstructuredGrid(.false.)

                end if

            end if

            ! Recalculate quality metric and find a split cell
            call qm%CalculateQualityMetrics(grid, options, magneticField, .true., .false.)

            
        end do

        ! Transform remaining pentagons into triangles
        !if (options%pents_to_tria) call grid%TransPentsToTrias() - TODO

        ! Ordening not necessary probably


    end subroutine

    subroutine Splitting(grid, qm, options, magneticField)

        ! Description
        !============
        ! Splits cells based on indication of certain cell

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)         :: grid
        type(QualityMetricUDT), intent(inout)   :: qm
        type(GAoptionsUDT), intent(inout)       :: options
        type(MagneticFieldUDT), intent(in)      :: magneticField

        ! Auxiliary
        integer(I8) :: cv

        ! Get splitting case and do the split on cell cv
        call grid%OneSplit(qm, magneticField,  options, qm%split_cv)

        ! If pentagons are not allowed

        if (options%no_pents) then
            cv = 0
            call grid%SplitPentsRec(qm, magneticField, options, cv)
        end if

        ! Do ordering for nice plots -  not necessary

    end subroutine

    subroutine OneSplit(grid, qm, magneticField, options, cv)

        ! Description
        !============
        ! Performs splitting on one cell and adjusts neighbouring cells
        ! topology

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)         :: grid
        type(QualityMetricUDT), intent(inout)   :: qm
        type(MagneticFieldUDT), intent(in)      :: magneticField
        type(GAoptionsUDT), intent(inout)       :: options
        integer(I8), intent(inout)              :: cv

        ! Auxiliary
        integer(I8) :: rface1, rface2, neig1, neig2, Qface, neig, common_vert
        integer(I8), allocatable :: fcs(:), Xcells(:)
        character(:), allocatable :: caseID, type
        logical, allocatable :: is_ordered(:), cells(:)



        ! Get faces
        fcs = GetCellFaceGA(grid%cell, cv)

        ! Determine if cv in Xradcells for geometric splitting strategy
        call grid%GetXcellsGeometric(Xcells)
        options%XpointSplitting = .false.
        if (any(cv == Xcells)) then
            options%XpointSplitting = .true.
        end if

        ! Determine the correct split case
        type = options%splittype
        select case (grid%cell%cflags%Get(cv))

        case (1:3)

            select case (size(fcs))

            case (3)

                call grid%DetermineTcaseID(cv, fcs, type, options%typeT, options, caseID, rface1, rface2, neig1, neig2)
                if (options%debug) print *, caseID

                select case (caseID)

                case ('34', '30')
                    call grid%SplitTQ(magneticField, cv, rface1, neig1, type, options)
                case ('35')
                    call grid%SplitTP(magneticField, cv, rface1, neig1, type, options)
                case ('430','034')
                    call grid%SplitQTB(magneticField, cv, rface1, rface2, neig1, neig2, type, options%typeT)
                case ('030')
                    call grid%SplitBTB(magneticField, cv, rface1, rface2, type, options%typeT)
                case ('333')
                    call grid%SplitTTT(magneticField, cv, rface1, rface2, neig1, neig2, type) 
                case ('334','433','033','330')
                    call grid%SplitTTQ(magneticField, cv, rface1, rface2, neig1, neig2, type)
                case default

                    print *, caseID
                    call gdErrorHandler('OneSplit: Tcase, case not implemented')

                end select

            case (4)

                call grid%DetermineQcaseID(cv, fcs, type, options, caseID, rface1, rface2, neig1, neig2)
                if (options%debug) print *, caseID

                select case (caseID)

                case ('444', '440', '044')
                    call grid%SplitQQQ(magneticField, cv, rface1, rface2, neig1, neig2, type, options)
                case ('344', '443')
                    !call grid%SplitQQT(magneticField, cv, rface1, rface2, neig1, neig2, type, options%typeT, options%QTtype, options)
                case ('043', '340')
                    !call grid%SplitBQT(magneticField, cv, rface1, rface2, neig1, neig2, type, options%QTtype)
                case ('343')
                    !call grid%SplitTQT(magneticField, cv, rface1, rface2, neig1, neig2, type, options%typeT, options%QTtype)
                case('4shaved')
                    !call grid%SpltQshaved(magneticField, cv, rface1,rface2,neig1,neig2)
                case default

                    print * , caseID
                    call gdErrorHandler('OneSplit: QcaseID not implemented')

                end select

            case (5)

                call grid%DeterminePcaseID(cv, fcs, type, caseID, Qface, neig, common_vert)
                if (options%debug) print *, caseID

                select case (caseID)

                !case ('55')
                case ('54','50')
                    call grid%SplitPQ(magneticField, cv, Qface, neig, common_vert, type, options)
                case ('53')
                    call grid%SplitPT(magneticField, cv, Qface, neig, common_vert, type, options%typeT)

                !case ('534')
                !case ('56')
                !case ('573')
                !case ('574')
                !case ('570')
                !case ('580')
                case default

                    print *, caseID
                    call gdErrorHandler('OneSplit: PcaseID not implemented')

                end select

            end select

        case (4) ! triangles that where turned into quad

            call grid%DetermineT4caseID(cv, fcs, options%splittype, caseID, rface1, neig1, common_vert)
            if (options%debug) print *, caseID

            select case (caseID)
            case ('30')
                call grid%SplitT4B(qm, magneticField, cv, rface1,common_vert, type, options)
            case ('31')
                !call grid%SplitT4BV()
            case ('34')
                call grid%SplitT4Q(magneticField, cv, rface1, common_vert, neig1, type)
            case ('35')
            case ('33')
            case default

                print *, caseID
                call gdErrorHandler('OneSplit: T4caseID not implemented')                

            end select

        case default

            call gdErrorHandler('OneSplit: not determined')

        end select

        if (options%debug) then 
            call grid%CheckUnstructuredGrid(.false.)
            allocate(cells(grid%cell%ntot), is_ordered(grid%cell%ntot))
            cells = .true.
            call grid%CheckVertOrder(is_ordered, cells)
            call grid%ReorderCellConn(is_ordered)
        end if

    end subroutine

    recursive subroutine SplitPentsRec(grid, qm, magneticField, options, cv)

        ! Description
        !============
        ! Splits pentagons recursively, so by splitting a cell.
        ! Also included triangles with one split face (triangle4)

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)         :: grid
        type(QualityMetricUDT), intent(inout)   :: qm
        type(MagneticFieldUDT), intent(in)      :: magneticField
        type(GAoptionsUDT), intent(inout)       :: options
        integer(I8), intent(inout)              :: cv
        
        ! Auxiliary
        integer(I8) :: i
        integer(I8), allocatable :: cells(:), pents(:)

        ! Select the cell to split
        ! Area constraints
        if (options%no_pents_area_split) then
            call AreaConstraintPents(grid, options, options%dist_function_threshold_split, cells)
        else
            cells = (/ (i, i = 1, grid%cell%ntot) /)
        end if

        ! Get the cell
        if (cv == 0) then

            ! Search pents to split
            call grid%GetPents(cells, pents)

            if (size(pents) /= 0) then

                call grid%SelectSplitPent(pents, options, cv)

                if (cv /= 0) then

                    call grid%OneSplit(qm, magneticField, options, cv)

                end if

            end if

        else 

            call grid%OneSplit(qm, magneticField, options, cv)

        end if

        ! Recursive part
        if (options%no_pents_area_split) then

            call grid%AreaConstraintPents(options, options%dist_function_threshold_split, cells)

        else 

            cells = (/ (i, i = 1, grid%cell%ntot)/)

        end if

        ! Get pentagonal and triangle4 cells
        call grid%GetPents(cells, pents)

        ! Eliminate cells which should not be split
        if (size(pents) /= 0) then

            call grid%SelectSplitPent(pents, options, cv)

            if (cv /= 0) then

                call grid%SplitPentsRec(qm, magneticField, options, cv)

            end if

        end if



    end subroutine

    subroutine SelectSplitPent(grid, pents, options, cv)

        ! Description
        !============
        ! Selects a pentagon or a triangle4 to be split based on surrounding cells

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(in)    :: grid
        integer(I8), intent(in)         :: pents(:)
        type(GAoptionsUDT), intent(in)  :: options
        integer(I8), intent(inout)      :: cv
 
        ! Auxiliary
        integer(I8) :: i, j, ic, indmax, splitface, counter, splitface_n
        integer(I8), allocatable :: fcs(:), neigs5(:), fcs_n(:), neigs(:)
        real(R8), allocatable :: dpsi(:)
        logical, allocatable :: log(:)


        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        cv = 0
        if (options%split_out) return

        do i = 1, size(pents)

            ic = pents(i)
            if (c%faceP2%Get(ic) == 5) then

                select case (options%splittype)

                case ('rad')

                    ! If the poloidal faces that should be split is also in an
                    ! existing pentagon, the splitting is not possible
                    neigs = GetCellNeigsGA(grid, ic)

                    if (any(c%faceP2%Get(neigs) == 5) .or. any(c%cflags%Get(neigs) == 4)) then

                        fcs = GetCellFaceGA(c, ic)
                        dpsi = abs(v%psi%Get(f%vert1%Get(fcs)) - v%psi%Get(f%vert2%Get(fcs)))
                        indmax = maxloc(dpsi,1)
                        splitface = fcs(indmax)

                        log = (c%faceP2%Get(neigs) == 5)
                        allocate(neigs5(count(log)))
                        neigs5 = pack(neigs, log)

                        counter = 0
                        do j = 1, size(neigs5)

                            fcs_n = GetCellFaceGA(c, neigs5(j))
                            if (any(splitface == fcs_n)) then

                                dpsi = abs(v%psi%Get(f%vert1%Get(fcs_n)) - v%psi%Get(f%vert2%Get(fcs_n)))
                                indmax = maxloc(dpsi, 1)
                                if (fcs_n(indmax) == splitface) then
                                    cv = ic
                                    exit
                                end if

                            else

                                counter = counter + 1

                            end if

                        end do

                        if (counter == size(neigs5)) then ! So for a neighboring pent the split face ws not a common face

                            cv = ic
                            return

                        end if
                            

                        ! Housekeeping
                        deallocate(neigs5)
                        
                    else

                        cv = ic
                        return

                    end if


                case ('pol')

                    ! If the radial faces that should be split is also in an
                    ! existing pentagon, the splitting is not possible, expect
                    ! when it is a correct splitting face for the pentagon
                    neigs = GetCellNeigsGA(grid, ic)
                    if (any(c%faceP2%Get(neigs) ==5) .or. any(c%cflags%Get(neigs) == 4)) then

                        fcs = GetCellFaceGA(c, ic)
                        call grid%GetRadSplitFacePent(fcs, splitface)

                        allocate(neigs5(count(c%faceP2%Get(neigs) == 5)))
                        neigs5 = pack(neigs, c%faceP2%Get(neigs) == 5)

                        counter = 0
                        do j = 1, size(neigs5)
                            fcs_n = GetCellFaceGA(c, neigs5(j))
                            if (any(splitface == fcs_n)) then
                                call grid%GetRadSplitFacePent(fcs, splitface_n)
                                if (splitface_n == splitface) then
                                    cv = ic
                                    exit
                                else
                                    cv = neigs5(j)
                                end if

                            else 

                                counter = counter + 1

                            end if
                        end do

                        if (counter == size(neigs5)) then ! So for the neigboring pents the split face was not a common face
                            cv = ic
                            return
                        end if
                    
                    else

                        cv = ic
                        return

                    end if

                end select

            else if (c%cflags%Get(ic) == 4) then

                cv = ic
                return

            else

                call gdErrorHandler('SelectSplitPent: not yet implemented')

            end if

        end do

        end associate

    end subroutine

    subroutine GetRadSplitFacePent(grid, faces, splitface)

        ! Description
        !============
        ! Give the radial split face of a pentagon (for poloidal splitting), the
        ! pentagon has good alignment. 
        ! Assumption splitted face is a radial face

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(in)    :: grid
        integer(I8), intent(in)         :: faces(:)
        integer(I8), intent(out)        :: splitface

        ! Auxiliary
        integer(I8) :: rface1, rface2, rface3, indmin
        real(R8), allocatable :: dpsi(:)


        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
        )

        ! Get radial faces
        dpsi = abs(v%psi%Get(f%vert1%Get(faces)) - v%psi%Get(f%vert2%Get(faces)))
        indmin = minloc(dpsi,1)
        rface1 = faces(indmin) ! Radial face1
        dpsi(indmin) = 1e99_R8

        indmin = minloc(dpsi,1)
        rface2 = faces(indmin) ! Radial face2
        dpsi(indmin) = 1e99_R8

        indmin = minloc(dpsi,1)
        rface3 = faces(indmin) ! Radial face3

        ! Radial faces with common vertex are faces on the side of
        ! the two quads
        if (HaveCommonVert(f, rface1, rface2)) then
            splitface = rface3
        else if (HaveCommonVert(f, rface1, rface3)) then
            splitface = rface2
        else
            splitface = rface1
        end if

        end associate

    end subroutine

    subroutine DetermineTcaseID(grid, cv, fcs, type, typeT, options, caseID, rface1, rface2, neig1, neig2)

        ! Description
        !============
        ! Determine the splitting case for a (aligned) triangle. Because in plasma
        ! grid always field aligned, a triangle at least has one aligned face.

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(in)            :: grid
        integer(I8), intent(in)                 :: cv, fcs(:)
        character(:), allocatable, intent(in)   :: type, typeT
        type(GAoptionsUDT), intent(in)          :: options
        character(:), allocatable, intent(out)  :: caseID
        integer(I8), intent(out)                :: rface1, rface2, neig1, neig2

        ! Auxiliary
        integer(I8), allocatable, dimension(:) :: rface1D, faceB, cells1, rfaces, int_faces
        integer(I8) :: neig1type, neig2type
        character(len=10) :: x1, x2

        ! Determine the case
        select case (type)

        case ('pol') ! Poloidal splitting

            select case (typeT)

            case ('stacked')

                ! Get radial face
                allocate(rface1D(count(grid%face%aligned%Get(fcs) == 1)))
                rface1D = pack(fcs, grid%face%aligned%Get(fcs) == 1)
                
                if (size(rface1D) == 0) call gdErrorHandler('DetermineTcaseID: Triangle should have an aligned face')

                rface1 = rface1D(1)

                call grid%GetNeigTypeFromFace(rface1, cv, neig1, neig1type)


                write (x1, '(I0)') neig1type
                caseID = '3' // trim(x1) 

            case ('cutcell')

                ! Get radial face
                allocate(rface1D(count(grid%face%aligned%Get(fcs) == 1)))
                rface1D = pack(fcs, grid%face%aligned%Get(fcs) == 1)
                
                if (size(rface1D) == 0) call gdErrorHandler('DetermineTcaseID: Triangle should have an aligned face')

                rface1 = rface1D(1)

                if (isBoundaryFaceGA(grid%face, rface1)) then

                    neig1 = 0
                    neig1type = 0

                    if (count(isBoundaryFaceGA(grid%face, fcs)) == 2) then

                        allocate(faceB(count(isBoundaryFaceGA(grid%face, fcs))))
                        faceB = pack(fcs, isBoundaryFaceGA(grid%face, fcs))
                        rface2 = faceB(1) 
                        neig2 = 0
                        neig2type = 0

                    else

                        call gdErrorHandler('DetermineTcaseID: geometry not considered')

                    end if

                else

                    ! Get radial neighbors
                    cells1 = GetFaceCellGA(grid%cell, rface1)
                    if (cells1(1) /= cv) then
                        neig1 = cells1(1)
                    else if (cells1(2) /= cv) then
                        neig1 = cells1(2)
                    end if
                    neig1type = grid%cell%faceP2%Get(neig1)

                    ! Get boundary face
                    allocate(faceB(count(isBoundaryFaceGA(grid%face, fcs))))
                    faceB = pack(fcs, isBoundaryFaceGA(grid%face, fcs))
                    rface2 = faceB(1)
                    neig2 = 0
                    neig2type = 0

                end if

                ! Make caseID
                write (x1, '(I0)') neig1type
                write (x2, '(I0)') neig2type
                caseID =  trim(x1) // '3' // trim(x2)
                

            end select

        case ('rad')

            ! Get the rfaces
            if (options%rad_type == 'no_aligned_faces') then

                ! Check which faces are intersected
                call grid%GetIntersectedPsiFaces(cv, fcs, int_faces)

                ! Extract faces
                rface1 = int_faces(1)
                rface2 = int_faces(2)

                call gdErrorHandler('DetermineTcaseID: triangle without aligned faces are not supported yet')

            else 

                ! Get the poloidal faces
                allocate(rfaces(count(grid%face%aligned%Get(fcs) == 0)))
                rfaces = pack(fcs, grid%face%aligned%Get(fcs) == 0)

                rface1 = rfaces(1)
                rface2 = rfaces(2)


            end if

            ! Get the neighbors and neighbor types
            call grid%GetNeigTypeFromFace(rface1, cv, neig1, neig1type)
            call grid%GetNeigTypeFromFace(rface2, cv, neig2, neig2type)

            ! Make caseID
            write (x1, '(I0)') neig1type
            write (x2, '(I0)') neig2type
            if (options%rad_type == 'no_aligned_faces') then

                caseID = trim(x1) // '3' // trim(x2) // 'A'

            else

                caseID  = trim(x1) // '3' // trim(x2)

            end if

        end select

    end subroutine

    subroutine DetermineQcaseID(grid, cv, fcs, type, options, caseID, rface1, rface2, neig1, neig2)

        ! Description
        !============
        ! Determine splitting case for quadrilateral cell
        
        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(in)            :: grid
        integer(I8), intent(in)                 :: cv, fcs(:)
        character(:), allocatable, intent(in)   :: type
        type(GAoptionsUDT), intent(in)          :: options
        character(:), allocatable, intent(out)  :: caseID
        integer(I8), intent(out)                :: rface1, rface2, neig1, neig2

        ! Auxiliary
        character(len=10) :: x1, x2
        integer(I8) :: ind, face_al, neig1type, neig2type, rface3
        integer(I8), allocatable, dimension(:) :: rfaces, rfaces2D, int_faces, cells1, cells2

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face &
            )

        select case (type)

        case ('pol')

            ! Get poloidal faces
            allocate(rfaces(count(f%aligned%Get(fcs) == 1)))
            rfaces = pack(fcs, f%aligned%Get(fcs) == 1)
            rface1 = rfaces(1)
            if (size(rfaces) == 2) then
                rface2 = rfaces(2)
            elseif (size(rfaces) == 1) then
                allocate(rfaces2D(count(isBoundaryFaceGA(f, fcs)))) 
                rfaces2D = pack(fcs, isBoundaryFaceGA(f, fcs))
                rface2 = rfaces2D(1)
            else 
                call gdErrorHandler('DetermineQcaseID: polooidal splitting of quad without aligned faces not implemented')
            end if


        case ('rad')

            ! Get poloidal faces
            allocate(rfaces(count(f%aligned%Get(fcs) == 0)))
            rfaces = pack(fcs, f%aligned%Get(fcs) == 0)
            rface1 = rfaces(1)
            rface2 = rfaces(2)
            if (size(rfaces) == 3) then !trapezoid
                ! Third non-aligned face
                rface3 = rfaces(3)

                ! Get the aligned face
                ind = findloc(f%aligned%Get(fcs), 1, 1)
                face_al = fcs(ind)

                if (.not.HaveCommonVert(f, rface1, face_al)) then
                    rface1 = rfaces(2)
                    rface2 = rfaces(3)
                elseif (.not.HaveCommonVert(f, rface2, face_al)) then
                    rface1 = rfaces(1)
                    rface2 = rfaces(3)
                elseif (.not.HaveCommonVert(f, rface3, face_al)) then
                    rface1 = rfaces(1)
                    rface2 = rfaces(2) 
                end if

            elseif (size(rfaces) == 4) then
                ! Find faces that are intersected by cv_psi
                call grid%GetIntersectedPsiFaces(cv, fcs, int_faces)

                rface1 = int_faces(1)
                rface2 = int_faces(2)

            end if

        end select

        ! Determine case
        if (isBoundaryFaceGA(f, rface1)) then
            neig1 = 0 ! indication that rface1 is boundary
            neig1type = 0
        else
            ! Get radial neighbours
            cells1 = GetFaceCellGA(c, rface1)
            if (cells1(1) /= cv) then
                neig1 = cells1(1)
            else if (cells1(2) /= cv) then
                neig1 = cells1(2)
            end if
            neig1type = c%faceP2%Get(neig1)             
        end if

        if (rface2 /= 0) then
            if (isBoundaryFaceGA(f, rface2)) then
                neig2 = 0
                neig2type = 0
            else
                ! Get radial neighbours
                cells2 = GetFaceCellGA(c, rface2)
                if (cells2(1) /= cv) then
                    neig2 = cells2(1)
                else if (cells2(2) /= cv) then
                    neig2 = cells2(2)
                end if
                neig2type = c%faceP2%Get(neig2)
            end if
        else
            ! Trapezoid so a boundary normally
            neig2 = 0
            neig2type = 0
        end if

        ! Make caseID
        write (x1, '(I0)') neig1type
        write (x2, '(I0)') neig2type
        caseID = trim(x1) // '4' // trim(x2) 

        ! Special case correction
        if (type == 'rad') then
            if (options%rad_type == 'shaved-off_tubes') then
                caseID = '4shaved'
            end if 
        end if

        end associate

    end subroutine

    subroutine DeterminePcaseID(grid, cv, fcs, type, caseID, Qface, neig, common_vert)

        ! Description
        !============
        ! Determine splitting case for pentagonal cells

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(in)            :: grid
        integer(I8), intent(in)                 :: cv, fcs(:)
        character(:), allocatable, intent(in)   :: type
        character(:), allocatable, intent(out)  :: caseID
        integer(I8), intent(out)                :: Qface, neig, common_vert

        ! Auxiliary
        integer(I8) :: neigtype, face_al, ind, i, Qface_old, v1, v2
        integer(I8), allocatable, dimension(:) :: rfaces, cells, cvs, vxs, &
            cvsv, vxsD
        real(R8) :: dpsi_al, dpsi_max, dpsi_min, h_pol, h_rad, AR, psic
        real(R8), allocatable :: dpsi(:), dpsi_rel(:)
        character(len=10) :: x1

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        ! Determine hanging node and Qface (opposite face)
        call grid%DetermineHangingNodePent(cv, fcs, type, common_vert, Qface, rfaces)

        ! Determine case
        if (isBoundaryFaceGA(f, Qface)) then
            neig = 0
            neigtype = 0
        else
            cells = GetFaceCellGA(c, Qface)
            if (cells(1) /= cv) then
                neig = cells(1) 
            else if (cells(2) /= cv) then
                neig = cells(2)
            end if
            neigtype = c%faceP2%Get(neig)
            if (c%cflags%Get(neig) == 4) then
                neigtype = 34
            end if
        end if

        ! Check whether pentagonal is trapezoidal because than
        ! above data can be incorrect

        if (type == 'rad') then

            ! For radial splitting it can be better to split pent different if it is a trapezoidal pent

            ! Other method to get common vert
            vxs = GetCellVertGA(c, cv)
            do i = 1, size(vxs)
                cvsv = GetVertCellGA(c, vxs(i))
                if (size(cvsv) == 3 .and. .not.isBoundaryVertGA(grid, vxs(i))) then
                    common_vert = vxs(i)
                    exit
                end if
            end do

            ! Determine whether it is a trapezoid 
            allocate(vxsD(count(vxs /= common_vert)))
            vxsD = pack(vxs, vxs /= common_vert)

            if (count(f%aligned%Get(fcs)==1) == 1) then
                ! So it is a trapezoid
                ind = findloc(f%aligned%Get(fcs), 1, 1)
                face_al = fcs(ind)
                dpsi_al = abs(v%psi%Get(f%vert1%Get(face_al)) - v%psi%Get(f%vert2%Get(face_al)))

                dpsi = abs(v%psi%Get(common_vert) - v%psi%Get(vxs))
                dpsi_max = maxval(dpsi)
                dpsi_min = minval(dpsi)
                dpsi_rel = dpsi/dpsi_max

                ! Recalculate because qm outdated at this point
                call grid%CalcHpol(cv, h_pol)
                call grid%CalcHrad(cv, h_rad)

                AR = h_pol/ h_rad

                if ((dpsi_min/dpsi_al .gt. 2) .or. ((AR .gt. 3) .and. (dpsi_min/dpsi_al .lt. 10))) then

                    ! using Ptrap = splitting in triangle and quad
                    neigtype = 6 ! just convention, will connect common face with a vertex
                    Qface = 0
                    common_vert = 0
                    neig = 0   
                else

                    ! else just go for psi splitting
                    Qface_old = Qface
                    psic = v%psi%Get(common_vert)

                    ! Determine psic crossing with a face
                    Qface = 0
                    do i = 1, size(fcs) 
                        v1 = f%vert1%Get(fcs(i))
                        v2 = f%vert2%Get(fcs(i))
                        if (psic .gt. min(v%psi%Get(v1),v%psi%Get(v2)) &
                            .and. psic .lt. max(v%psi%Get(v1),v%psi%Get(v2))) then
                                Qface = fcs(i)
                                exit
                        end if
                    end do

                    if (Qface /= 0) then
                        call gdErrorHandler('DeterminePcaseID: No psic crossing found, which should be impossible!')
                    end if

                    if (isBoundaryFaceGA(f, Qface)) then
                        if (Qface_old /= Qface) then

                            ! Here the psi splitting
                            neigtype = 80
                            neig = 0

                        else 

                            ! Need psi splitting here
                            neigtype = 70
                            neig = 0

                        end if
                    else 

                        ! Determine neig 
                        cvs = GetFaceCellGA(c, Qface)
                        if (cvs(1) /= cv) then
                            neig = cvs(1)
                        else if (cvs(2) /= cv) then
                            neig = cvs(2)
                        end if
                        if (c%vertP2%Get(neig) == 3) then
                            neigtype = 73
                        elseif (c%vertP2%Get(neig) == 4) then
                            neigtype = 74
                        end if

                    end if

                end if

            end if

        end if



        ! Construct caseID string
        write (x1, '(I0)') neigtype
        caseID = '5' // trim(x1)

        end associate

    end subroutine

    subroutine DetermineT4caseID(grid, cv, fcs, type, caseID, rface1, neig, common_vert)

        ! Description
        !============
        ! Determine the splitting case for a triangle with one face split in two
        ! which is indicated with a cflags == 4
        !
        ! cv: cellnumber
        ! fcs: faces of the cell
        ! type: pol or rad

        ! caseID: point to correct splitting routine
        ! rface1: that will be splitted by the splitting routine
        !         for radial splitting this is the longest poloidal face opposite to the common vert
        !         for poloidal splitting this is the longest radial face opposite
        !         to the common vert
        ! neig: the cell number of a neigboring cell if any, 0 if none
        ! common vert: is the common vertex of the earlier splitting face

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(in)            :: grid
        integer(I8), intent(in)                 :: cv, fcs(:)
        character(:), allocatable, intent(in)   :: type
        character(:), allocatable, intent(out)  :: caseID
        integer(I8), intent(out)                :: rface1, neig, common_vert
        
        ! Auxiliary     
        integer(I8) :: ind, pface1, pface2, neigtype, vp1(2), vp2(2)
        integer(I8), allocatable, dimension(:) :: split_faces, facesD, pfaces, cvs, vxs
        real(R8), allocatable :: dpsi_f(:)
        real(R8) :: fcS1, fcS2
        character(len=10) :: x1

        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )
 
        ! rface1
        dpsi_f = abs(v%psi%Get(f%vert1%Get(fcs)) - v%psi%Get(f%vert2%Get(fcs)))
        ind = maxloc(dpsi_f,1)

        select case (type)
        case ('pol')

            ! rface1 => shoulde be the longest face of the poloidal faces
            pface1 = fcs(ind) !poloidal face1

            dpsi_f(ind) = 0
            ind = maxloc(dpsi_f,1)
            pface2 = fcs(ind) !poloidal face1

            ! Take the longest of the two faces
            vp1 = [f%vert1%Get(pface1), f%vert2%Get(pface1)]
            vp2 = [f%vert1%Get(pface2), f%vert2%Get(pface2)]
            fcS1 = sqrt( (v%x%Get(vp1(2)) - v%x%Get(vp1(1)))**2 + (v%y%Get(vp1(2)) - v%y%Get(vp1(1)))**2)
            fcS2 = sqrt( (v%x%Get(vp2(2)) - v%x%Get(vp2(1)))**2 + (v%y%Get(vp2(2)) - v%y%Get(vp2(1)))**2)

            if (fcS1 .gt. fcS2) then
                rface1 = pface1
            else if (fcS2 .gt. fcS1) then
                rface1 = pface2
            else
                call gdErrorHandler('DetermineT4caseID: Both poloidal faces in the triangle4' // &
                        & 'are even long which is geometrically not possible!')
            end if

            ! Common vert
            allocate(split_faces(count(fcs /= pface2)))
            split_faces = pack(fcs, fcs /= pface2)
            common_vert = GetCommonVert(f, split_faces(1), split_faces(2))

            ! Determine CaseID
            vxs = GetCellVertGA(c, cv)

            if (isBoundaryFaceGA(f, rface1)) then
                neig = 0
                neigtype = 0
                ! SplitT4B
            else if (count(isBoundaryVertGA(grid, vxs)) == 1) then
                neig = 0
                neigtype = 1
            else
                cvs = GetFaceCellGA(c, rface1)
                if (cvs(1) /= cv) then
                    neig = cvs(1)
                else if (cvs(2) /= cv) then
                    neig = cvs(2)
                end if
                neigtype = c%faceP2%Get(neig)
            end if

        case ('rad')

            ! rface1
            ! Get face that will be split
            pface1 = fcs(ind) !poloidal face1
            
            ! common_vert
            ! Determine the common vertex (vertex on split face)
            ! lay between two faces with maximaal dpsi expect the rface1
            allocate(facesD(count(fcs /= rface1)))
            facesD = pack(fcs, fcs /= rface1)   
            
            ! Get poloidal faces
            allocate(pfaces(count(f%aligned%Get(facesD) == 0)))
            pfaces = pack(facesD, f%aligned%Get(facesD) == 0)

            if (size(pfaces) /= 2) call gdErrorHandler('DetermineT4caseID: Something wrong with face%aligned')

            pface1 = pfaces(1)
            pface2 = pfaces(2)

            ! Common vert
            if (HaveCommonVert(f, pface1, pface2)) then
                common_vert = GetCommonVert(f, pface1, pface2)
            else
                call gdErrorHandler('DetermineT4caseID: No common vertex: probably the triangle a no aligned face.')
            end if

            ! Neig
            if (isBoundaryFaceGA(f, rface1)) then
                neig = 0
                neigtype = 0
            else 
                cvs = GetFaceCellGA(c, rface1)
                if (cvs(1) /= cv) then
                    neig = cvs(1)
                else if (cvs(2) /= cv) then
                    neig = cvs(2)
                end if
                neigtype = c%faceP2%Get(neig)
            end if

        end select

        write (x1, '(I0)') neigtype
        caseID = '3' // trim(x1)

        end associate
    end subroutine

    subroutine GetNeigTypeFromFace(grid, rface, cv, neig, neigtype)

        ! Description
        !============
        ! Getting neighbor cell and neigtype of a radial face from a split cell

        ! Declare variables
        !==================
        class(GAGridUDT), intent(in)    :: grid
        integer(I8), intent(in)         :: cv, rface
        integer(I8), intent(out)        :: neig, neigtype

        ! Auxiliary
        integer(I8), allocatable :: cvs(:)
        
        if (isBoundaryFaceGA(grid%face, rface)) then

            neig = 0
            neigtype = 0

        else

            ! Get radial neighbors
            cvs = GetFaceCellGA(grid%cell, rface)
            if (cvs(1) /= cv) then
                neig = cvs(1)
            else if (cvs(2) /= cv) then
                neig = cvs(2)
            end if
            neigtype = grid%cell%faceP2%Get(neig)

        end if

    end subroutine

    subroutine GetIntersectedPsiFaces(grid, cv, fcs, int_faces, psic)

        ! Description
        !============
        ! Returns the two faces which were intersecting by the average psi of
        ! cell cv of a specified psic. The intersection points are stored in isx and isy.

        ! Declare variables
        !==================
        class(GAGridUDT), intent(in)            :: grid
        integer(I8), intent(in)                 :: cv, fcs(:)
        integer(I8), allocatable, intent(out)   :: int_faces(:)
        real(R8), optional                      :: psic

        ! Auxiliary
        integer(I8) :: ii, int_facesD(4), v1, v2, i
        integer(I8), allocatable :: vxs(:)
        real(R8) :: isxD(4), isyD(4), v1p, v2p, t0
        real(R8), allocatable :: isx(:), isy(:)

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        ! Initialize
        int_facesD = 0
        isx = 0
        isy = 0
        ii = 0

        ! Get the intersected faces
        if (.not.present(psic)) then
            vxs = GetCellVertGA(c, cv)
            psic = 0.5_R8*(maxval(v%psi%Get(vxs)) + minval(v%psi%Get(vxs)))
        end if 

        do i = 1, c%vertP2%Get(cv)
            v1 = f%vert1%Get(fcs(i))
            v2 = f%vert2%Get(fcs(i))
            v1p = v%psi%Get(v1)
            v2p = v%psi%Get(v2)

            if ((psic .gt. min(v1p, v2p)) .and. (psic .lt. max(v1p, v2p))) then
                t0 = (psic - v1p) / (v2p - v1p)
                ii = ii + 1
                isxD(ii) = v%x%Get(v1) + t0 * (v%x%Get(v2) - v%x%Get(v1))
                isyD(ii) = v%y%Get(v1) + t0 * (v%y%Get(v2) - v%y%Get(v1))
                int_facesD(ii) = fcs(i)
            end if

        end do

        if (.not.present(psic)) then
            if (ii /= 2) then
                call gdErrorHandler('GetIntersectedPsiFaces: no 2 intersections found')
            end if
        end if

        ! Trim
        isx = isxD(1:ii)
        isy = isyD(1:ii)
        int_faces = int_facesD(1:ii)


        end associate

    end subroutine  

    subroutine GetPents(grid, cells, pents)

        ! Description
        !============
        ! Get the pentagonal and triangle4 cells

        ! Declare variables
        !==================
        class(GAGridUDT), intent(in)            :: grid
        integer(I8), allocatable, intent(in)    :: cells(:)
        integer(I8), allocatable, intent(out)   :: pents(:)

        ! Auxiliary
        logical, allocatable :: log(:), log2(:)
        integer(I8) :: i
        integer(I8), allocatable :: pentsD(:), triangle4(:), indcv(:)

        log = (grid%cell%faceP2%Get(cells) == 5)
        allocate(pentsD(count(log)))
        pentsD = pack(cells, log)
        log2 = (grid%cell%cflags%Get() == 4)
        allocate(triangle4(count(log2)))
        indcv = (/ (i, i = 1, grid%cell%ntot)/)
        triangle4 = pack(indcv, log2)
        pents = [pentsD, triangle4]

    end subroutine

    subroutine SplitTQ(grid, magneticField, cv, tface, neig1, type, options)

        ! Description
        !============
        ! Splits a triangle 
        ! Tface is the aligned faces of the triangle

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)     :: grid 
        type(MagneticFieldUDT), intent(in)  :: magneticField
        integer(I8), intent(in)             :: cv, tface, neig1
        character(:), allocatable           :: type
        type(GAoptionsUDT), intent(in)      :: options

        ! Auxiliary
        integer(I8) :: v1n, free_vert, f1n, f2n, &
            f3n, cv1, cv2, face_old
        integer(I8), allocatable, dimension(:) :: fcs, vxs, face_rem, cells

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face &
            )

        ! Split tface
        call grid%SplitFace(magneticField, tface, v1n, f1n, f2n)

        ! Determine free vertex of the triangle by eliminate vertices of tface
        call grid%DetermineFreeVertTria(cv, tface, free_vert)

        ! Make new face     
        call grid%GetFaceNumber(v1n, free_vert, 3, f3n)

        ! Add new vert to correct flux surface
        if (.not.options%slab) then
            select case (type)
            case ('pol')
                allocate(vxs(1))
                vxs = v1n
                fcs = [f1n, f2n]
                face_old = tface
            case ('rad')
                allocate(fcs(1))
                vxs = [v1n, free_vert]
                fcs = f3n
                face_old = 0
            end select
        end if
        call grid%AddVertToFsVx(vxs, fcs, face_old, type)

        ! Build new cells

        ! Split cv in two triangles
        call grid%SplitTriaStacked(cv, Tface, v1n, free_vert, cv1, cv2)

        ! Neig1 => becomes pentagonal cell
        if (neig1 /= 0) call grid%QuadToPent(neig1, Tface, f1n, f2n, v1n)

        ! Remove tface
        face_rem = Tface
        call grid%RemoveFaces(face_rem)

        ! Determine cflags
        call grid%DetermineCflags(cells)
        
        end associate

    end subroutine

    subroutine SplitTP(grid, magneticField, cv, tface, neig1, type, options)

        ! Description
        !============
        ! Splits a triangle neighbored by a pentagon in half

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)     :: grid
        type(MagneticFieldUDT), intent(in)  :: magneticField         
        integer(I8), intent(in)             :: cv, tface, neig1
        character(:), allocatable           :: type
        type(GAoptionsUDT), intent(in)      :: options

        ! Auxiliary
        integer(I8) :: i, counter, free_vert, v1n, common_vert, qface, &
            pfaces(5), face_old, cvT1, cvT2, cvP1, cvP2, f1n, f2n, f3n, f4n 
        integer(I8), allocatable, dimension(:) :: fcs, rfaces, face_rem, cells, &
            vxs

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face &
            )

        ! Make new vertices
        call grid%SplitFace(magneticField, tface, v1n, f1n, f2n)

        ! Determine free vertex of triangle
        call grid%DetermineFreeVertTria(cv, tface, free_vert)

        ! Determine common vert of pentagon
        fcs = GetCellFaceGA(c, neig1)
        call grid%DetermineHangingNodePent(neig1, fcs, type, common_vert, qface, rfaces)

        ! Make new faces
        call grid%GetFaceNumber(v1n, free_vert, 3, f3n)
        call grid%GetFaceNumber(v1n, common_vert, 3, f4n)

        ! For alignement of f4_n, determine the two faces which are not rfaces
        pfaces = 0
        counter = 0
        do i = 1, size(fcs)
            if (.not.any(fcs(i) == rfaces)) then
                counter = counter + 1
                pfaces(counter) = fcs(i)
            end if
        end do
        if (counter .lt. 2) call gdErrorHandler('SplitTP: not enough poloidal faces in the pentagon')
        if (f%aligned%Get(pfaces(1)) == 1 .and.  f%aligned%Get(pfaces(2)) == 1) call f%aligned%Set(f4n, 1)

        ! Add new vert to correct flux sruface
        select case (type)
        case ('pol')
            allocate(vxs(1))
            vxs = v1n
            fcs = [f1n, f2n]
            face_old = tface
        case ('rad')
            vxs = [v1n, free_vert, common_vert]
            fcs = [f3n, f4n]
            face_old = 0
        end select
        call grid%AddVertToFsVx(vxs, fcs, face_old, type)

        ! Build new cell
        !===============
        ! Split cv in two triangles
        call grid%SplitTriaStacked(cv, tface, v1n, free_vert, cvT1, cvT2)

        ! Split pent in two quads
        call grid%SplitCenterPent(neig1,tface, v1n, common_vert, f4n, f1n, f2n, cvP1, cvP2)

        ! Remove tface
        face_rem = tface
        call grid%RemoveFaces(face_rem)

        ! Determine cflags
        cells = [cvT1, cvT2, cvP1, cvP2]
        call grid%DetermineCflags(cells)

        end associate

    end subroutine   

    subroutine SplitQTB(grid, magneticField, cv, rface1, rface2, neig1, neig2, type, typeT)

        ! Description
        !============
        ! Does the radial or poloidal splitting of a starting triangle which is a boundary cell
        ! and has a quad neigbour.

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)     :: grid 
        type(MagneticFieldUDT), intent(in)  :: magneticField        
        integer(I8), intent(in)             :: cv, rface1, rface2, neig1, neig2
        character(:), allocatable           :: type, typeT

        ! Auxiliary
        integer(I8) :: Bface, Iface, ic, v1n, v2n, f11n, f12n, &
            f21n, f22n, perp_face, face_old1, cv1, cv2, f3n
        integer(I8), allocatable :: vxs1(:), fcs1(:), &
            vxs2(:), fcs2(:), face_rem(:), cells(:)

        ! Check case
        if (type == 'pol' .and. typeT == 'stacked') &
            call gdErrorHandler('PolSplitQTB: wrong caseID for poloidal stacked triangle splitting')

        ! Determine which is the boundary face and with the internal face
        if (neig1 == 0) then
            Bface = rface1
            Iface = rface2
            ic = neig2 ! internal cell
        else 
            Bface = rface2
            Iface = rface1
            ic = neig1
        end if

        ! Determine perpendicular face
        call grid%DeterminePerpFaceTria(cv, Bface, Iface, perp_face)

        ! Make new vertices and faces - just geometric middle
        call grid%SplitFace(magneticField, Bface, v1n, f11n, f12n)
        call grid%SplitFace(magneticField, Iface, v2n, f21n, f22n)

        call grid%GetFaceNumber(v1n, v2n, 3, f3n)
        call grid%face%aligned%Set(f3n, grid%face%aligned%Get(perp_face))

        ! Add new vertices to correct flux sruface
        select case (type)
        case ('pol')
            allocate(vxs1(1), vxs2(1))
            vxs1 = v1n
            fcs1 = [f11n, f12n]
            vxs2 = v2n
            fcs2 = [f21n, f22n]        
            call grid%AddVertToFsVx(vxs1, fcs1, Bface, type)
            call grid%AddVertToFsVx(vxs2, fcs2, Iface, type)
        case ('rad')
            allocate(fcs1(1))
            vxs1 = [v1n, v2n]
            fcs1 = f3n
            face_old1 = 0
            call grid%AddVertToFsVx(vxs1, fcs1, face_old1, type)
        end select

        ! Build new cells
        call grid%QuadToPent(ic, Iface, f21n, f22n, v2n)

        ! Split the center tria
        call grid%SplitCenterTria(cv, Bface, v1n, v2n, perp_face, cv1, cv2)

        ! Remove Iface and Bface
        face_rem = [Bface, Iface]
        call grid%RemoveFaces(face_rem)

        ! Determine cflags
        cells = [cv1, cv2]
        call grid%DetermineCflags(cells)
         
    end subroutine   

    subroutine SplitBTB(grid, magneticField, cv, rface1, rface2, type, typeT)

        ! Description
        !============
        ! Does the radial or poloidal  splitting of a starting triangle which is a boundary cell
        ! and has two boundary faces.

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)     :: grid 
        type(MagneticFieldUDT), intent(in)  :: magneticField        
        integer(I8), intent(in)             :: cv, rface1, rface2
        character(:), allocatable           :: type, typeT

        ! Auxiliary
        integer(I8) :: perp_face, cv1, cv2, Bface, Iface, &
            f11n, f12n, f21n, f22n, f3n, v1n, v2n 
        integer(I8), allocatable, dimension(:) :: vxs1, &
            vxs2, fcs1, fcs2, cells, face_rem


        if (type == 'pol' .and. typeT == 'cutcell') then

            ! Determine which is the boundary face and with the internal face
            Bface = rface1
            Iface = rface2

            ! Determine perpendicular face
            call grid%DeterminePerpFaceTria(cv, Bface, Iface, perp_face)

            ! Split faces
            call grid%SplitFace(magneticField, Bface, v1n, f11n, f12n)
            call grid%SplitFace(magneticField, Iface, v2n, f21n, f22n)

            ! Make new face
            call grid%GetFaceNumber(v1n, v2n, 3, f3n)
            call grid%face%aligned%Set(f3n, grid%face%aligned%Get(perp_face))

            ! Add new vert to flux surface.
            allocate(vxs1(1), vxs2(1))
            vxs1 = v1n
            vxs2 = v2n 
            fcs1 = [f11n, f12n]
            fcs2 = [f21n, f22n]
            call grid%AddVertToFsVx(vxs1, fcs1, Bface, type)
            call grid%AddVertToFsVx(vxs2, fcs2, Iface, type)

            ! Build new cells
            ! cv => becomes pent

            ! Split the cneter tria
            call grid%SplitCenterTria(cv, Bface, v1n, v2n, perp_face, cv1, cv2)

            ! Remove faces
            face_rem = [Bface, Iface]
            call grid%RemoveFaces(face_rem)

            ! Determine cflags
            cells = [cv1, cv2]
            call grid%DetermineCflags(cells)

        else if (type ==  'rad') then

            call gdErrorHandler('PolSplitBTB: radial splitting not implemented')

        else if (typeT == 'stacked') then

            call gdErrorHandler('PolSplitBTB: wrong caseID for poloidal stacked triangle splitting')

        end if



    end subroutine   
 
    subroutine SplitTTT(grid, magneticField, cv, rface1, rface2, neig1, neig2, type)

        ! Description
        !============
        !  Split a triangle with a triangle as neighbour along both sides. Only
        ! available for radial splitting

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)     :: grid
        type(MagneticFieldUDT), intent(in)  :: magneticField         
        integer(I8), intent(in)             :: cv, rface1, rface2, neig1, neig2
        character(:), allocatable           :: type

        ! Auxiliary
        integer(I8) :: v1n, v2n, f11n, f12n, f21n, f22n, perp_face, cv1, cv2, f3n
        integer(I8), allocatable, dimension(:) :: vxs1, fcs1, &
            face_rem, cells
        real(R8) :: psic
        character(:), allocatable :: typeV

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        ! Check type
        if (type == 'rad') then

            ! Determine perpendicular face
            call grid%DeterminePerpFaceTria(cv, rface1, rface2, perp_face)
            
            ! Split faces
            psic = 0.5_R8 * (v%psi%Get(f%vert1%Get(rface1)) + v%psi%Get(f%vert2%Get(rface1)))
            typeV = 'psi'
            call grid%SplitFace(magneticField, rface1, v1n, f11n, f12n, typeV, psic)
            call grid%SplitFace(magneticField, rface2, v2n, f21n, f22n, typeV, psic)

            ! Make new face
            call grid%GetFaceNumber(v1n, v2n, 3, f3n)
            call grid%face%aligned%Set(f3n, grid%face%aligned%Get(perp_face))

            ! Add new vertices to flux surface
            allocate(fcs1(1))
            vxs1 = [v1n, v2n]
            fcs1 = f3n
            call grid%AddVertToFsVx(vxs1, fcs1, 0, type)

            ! Turn neighboring triangles to triangle4
            call grid%TriaToQuad(neig1, rface1, f11n, f12n, v1n)
            call grid%TriaToQuad(neig2, rface2, f21n, f22n, v2n)

            ! Split the center tria
            call grid%SplitCenterTria(cv, rface1, v1n, v2n, perp_face, cv1, cv2)

            ! Remove faces
            face_rem = [rface1, rface2]
            call grid%RemoveFaces(face_rem)

            ! Determine cflags
            cells = [cv1, cv2]
            call grid%DetermineCflags(cells)


        else if (type == 'pol') then

            call gdErrorHandler('SplitTTT for poloidal splitting: Not supported')

        end if

        end associate

    end subroutine 

    subroutine SplitTTQ(grid, magneticField, cv, rface1, rface2, neig1, neig2, type)

        ! Description
        !============
        ! Splits a triangle neigbored by a triangle and a quad in the splitting
        ! direction

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)     :: grid
        type(MagneticFieldUDT), intent(in)  :: magneticField         
        integer(I8), intent(in)             :: cv, rface1, rface2, neig1, neig2
        character(:), allocatable           :: type

        ! Auxiliary
        integer(I8) :: Tface, Qface, cvT, cvQ, nf1, nf2, cv1, cv2, &
            f11n, f12n, f21n, f22n, v1n, v2n, f3n, perp_face, v1B, v2B
        integer(I8), allocatable :: vxs1(:), fcs1(:), face_rem(:), cells(:)
        real(R8) :: psic
        character(:), allocatable :: typeV

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        if (type == 'rad') then

            ! Determine which is the boundary face and with the internal face
            if (neig1 == 0) then
                Tface = rface2
                Qface = rface1
                cvT = neig2
                cvQ = neig1 
            else if (neig2 == 0) then
                Tface = rface1
                Qface = rface2
                cvT = neig1
                cvQ = neig2                
            else
                nf1 = c%faceP2%Get(neig1)
                nf2 = c%faceP2%Get(neig2)
                if (nf1 == 3 .and. nf2 == 4) then
                    Tface = rface1
                    Qface = rface2
                    cvT = neig1
                    cvQ = neig2
                else if (nf1 == 4 .and. nf2 == 3) then
                    Tface = rface2
                    Qface = rface1
                    cvT = neig2
                    cvQ = neig1         
                else
                    call gdErrorHandler('Something wrong')
                end if
            end if 

            ! Determine perpendicular face
            call grid%DeterminePerpFaceTria(cv, Qface, Tface, perp_face)

            ! Split faces
            v1B = f%vert1%Get(Qface)
            v2B = f%vert2%Get(Qface)
            psic = 0.5_R8 * (v%psi%Get(v1B) + v%psi%Get(v2B))
            typeV = 'psi'
            call grid%SplitFace(magneticField, Tface, v1n, f11n, f12n, typeV, psic)
            call grid%SplitFace(magneticField, Qface, v2n, f21n, f22n, typeV, psic)

            ! Make new face
            call grid%GetFaceNumber(v1n, v2n, 3, f3n)
            call f%aligned%Set(f3n, f%aligned%Get(perp_face))

            ! Add new vertices to flux surface
            allocate(fcs1(1))
            vxs1 = [v1n, v2n]
            fcs1 = f3n
            call grid%AddVertToFsVx(vxs1, fcs1, 0, type)

            ! Build new cells
            ! cvQ => becomes pent
            if (cvQ /= 0) call grid%QuadToPent(cvQ, Qface, f21n, f22n, v2n)

            ! Split the center tria
            call grid%SplitCenterTria(cv, Tface, v1n, v2n, perp_face, cv1, cv2)

            ! Turn neighbouring triangle to triangle4
            call grid%TriaToQuad(cvT, Tface, f11n, f12n, v1n)

            ! Remove faces
            face_rem = [Tface, Qface]
            call grid%RemoveFaces(face_rem)

            ! Determine cflags
            cells = [cv1, cv2]
            call grid%DetermineCflags(cells)

            
        else if (type == 'pol') then

            call gdErrorHandler('SplitTTQ: poloidal splitting not supported here')

        end if

        end associate

    end subroutine 

    subroutine SplitQQQ(grid, magneticField, cv, rface1, rface2, neig1, neig2, type, options)

        ! Description
        !============
        !  Splitting of a quadrilateral cell cv where faces rface1 and
        ! rface2 are split in half. This results in a configuration of pentagonal - 2x quad - pentagonal
        ! cell.

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)         :: grid
        type(MagneticFieldUDT), intent(in)      :: magneticField
        integer(I8), intent(in)                 :: cv, rface1, rface2
        integer(I8), intent(inout)              ::  neig1, neig2
        character(:), allocatable, intent(in)   :: type
        type(GAoptionsUDT), intent(in)          :: options

        ! Auxiliary
        integer(I8) :: i, vxs1r(2), vxs2r(2), v1n, v2n, f11n, f12n, f21n, f22n, f3n, &
            f1n, f2n, f4n, f5n, Qface1, Qface2, vxsT, &
            vxsQ1, vxsQ2, cv1, s, cv2, neig, split_face
        integer(I8), allocatable, dimension(:) :: verts, face_rem, cells, &
            vxs1, vxs2, fcs1, fcs2, cvsQ1, range, vxs_newp, fcs_newp, vxsD, fcsD, &
            vertsT, facesT, fcs_al, fcs, vxs, perp_faces, int_faces, cvs, fcs_cv1, &
            vxs_fcs, ind, fcs_vxsQ1, query_faces
        real(R8) :: psic, vp1(2), vp2(2)
        character(:), allocatable :: typeV, type2
        logical :: normal_flag
        logical, allocatable :: is_ordered(:), cells_log(:)


        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )
          
        ! Get perpendicular faces
        call grid%DeterminePerpFaceQuad(cv, rface1, rface2, perp_faces)

        ! Flag to switch method
        normal_flag = .true.

        ! Check if center quad is a trapezoid
        psic = 0_R8
        if (.not.options%XpointSplitting) then

            if (options%splittype == 'rad') then

                if (count(f%aligned%Get(perp_faces) == 1) .lt. 2) then

                    verts = GetCellVertGA(c, cv)
                    psic = 0.5_R8 * (maxval(v%psi%Get(verts)) + minval(v%psi%Get(verts)))

                    ! Check intersect of rface1 and rface2 by psic
                    vxs1r = [f%vert1%Get(rface1), f%vert2%Get(rface1)]
                    vp1 = v%psi%Get(vxs1r)
                    vxs2r = [f%vert1%Get(rface2), f%vert2%Get(rface2)]
                    vp2 = v%psi%Get(vxs2r)
                    if (.not.( (psic > min(vp1(1),vp1(2))) .and. (psic < max(vp1(1),vp1(2)))) &
                        .or. .not.((psic > min(vp2(1),vp2(2))) .and. (psic < max(vp2(1),vp2(2))))) then
                        ! One of the faces is not intersected, so split the
                        ! trapezoid in a triangle and a pentagon, which is
                        ! later splitted poloidally
                        normal_flag = .false.

                    else

                        typeV = 'psi'

                    end if
                else

                    typeV = 'geometric'

                end if 
            else if (options%splittype == 'pol') then

                typeV = 'pol'

            end if
        else

            typeV = 'geometric'

        end if

        if (normal_flag) then

            ! Split faces
            call grid%SplitFace(magneticField, rface1, v1n, f11n, f12n, typeV, psic)
            call grid%SplitFace(magneticField, rface2, v2n, f21n, f22n, typeV, psic)

            ! Make new face
            call grid%GetFaceNumber(v1n, v2n, 3, f3n)
            if (grid%face%aligned%Get(perp_faces(1)) == 1 .or. grid%face%aligned%Get(perp_faces(1)) == 1) &
                call grid%face%aligned%Set(f3n, 1)
            if (type == 'pol') then
                call grid%face%aligned%Set(f3n, 0)
            end if

            ! Add vert to flux surface
            select case (options%splittype)
            case ('pol')
                allocate(vxs1(1), vxs2(1))
                vxs1 = v1n
                vxs2 = v2n
                fcs1 = [f11n, f12n]
                fcs2 = [f21n, f22n]
                call grid%AddVertToFsVx(vxs1, fcs1, rface1, options%splittype)
                call grid%AddVertToFsVx(vxs2, fcs2, rface2, options%splittype)
            case ('rad')
                allocate(fcs1(1))
                vxs1 = [v1n, v2n]
                fcs1 = f3n  
                call grid%AddVertToFsVx(vxs1, fcs1, rface1, options%splittype)                         
            end select

            ! Build new cells
            ! neig1 => becomes pent, neig2 => becomes pent
            if (neig1 /= 0) call grid%QuadToPent(neig1, rface1, f11n, f12n, v1n)
            if (neig2 /= 0) call grid%QuadToPent(neig2, rface2, f21n, f22n, v2n)

            ! Split center quad
            call grid%SplitCenterQuad(cv, rface1, v1n, v2n, cv1, cv2)

            ! Remove faces
            face_rem = [rface1, rface2]
            call grid%RemoveFaces(face_rem)

            ! Determine cflags
            cells = [cv1, cv2]
            call grid%DetermineCflags(cells)

        else

            ! Not regular case
            ! Other splitting method:  pure one psi, split trapezoid in a triangle
            ! and a pentagon (ref. SplitPtrapPSIB2.m)

            ! Find intersection faces
            fcs = GetCellFaceGA(c, cv)
            allocate(fcs_al(count(f%aligned%Get(fcs) == 1)))
            fcs_al = pack(fcs, f%aligned%Get(fcs) == 1)
            vxs = GetCellVertGA(c, cv)
            call grid%GetIntersectedPsiFaces(cv, fcs, int_faces)

            ! Identify boundary face, if any
            if (isBoundaryFaceGA(f, int_faces(1))) then
                Qface2 = int_faces(1)
                Qface1 = int_faces(2)
            else if (isBoundaryFaceGA(f, int_faces(2))) then
                Qface2 = int_faces(2)
                Qface1 = int_faces(1)
            else 
                call gdErrorHandler('SplitQQQ: internal trapezoid (= 1 aligned face); splitting not implemented')
            end if

            ! Identify quad neighbour
            cvsQ1 = GetFaceCellGA(c, Qface1)
            if (cvsQ1(1) /= cv) then
                neig1 = cvsQ1(1)
            else if (cvsQ1(2) /= cv) then
                neig1 = cvsQ1(2)
            endif

            ! Start the splitting
            typeV = 'psi'
            call grid%SplitFace(magneticField, Qface1, v1n, f1n, f2n, typeV, psic)
            call grid%SplitFace(magneticField, Qface2, v2n, f3n, f4n, typeV, psic)

            vxsT = GetCommonVert(f, Qface1, Qface2)
            if (f%vert1%Get(Qface1) /= vxsT) then
                vxsQ1 = f%vert1%Get(Qface1)
            else if (f%vert2%Get(Qface1) /= vxsT) then
                vxsQ1 = f%vert2%Get(Qface1)
            end if 
            if (f%vert1%Get(Qface2) /= vxsT) then
                vxsQ2 = f%vert1%Get(Qface2)
            else if (f%vert2%Get(Qface2) /= vxsT) then
                vxsQ2 = f%vert2%Get(Qface2)
            end if 

            ! Make new face
            call grid%GetFaceNumber(v1n, v2n, 3, f5n)

            ! Add vertex to flux surface
            allocate(fcs1(1))
            vxs1 = [v1n, v2n]
            fcs1 = f5n
            call grid%AddVertToFsVx(vxs1, fcs1, 0, options%splittype)

            ! Build new cells
            !================
            ! Quad becomes Pentagon
            allocate(vxsD(count(vxs /= vxsT)))
            vxsD = pack(vxs, vxs /= vxsT)
            vxs_newp = [vxsD, v1n, v2n]
            if (size(vxs_newp) /= 5) call gdErrorHandler('SplitQQQ: Something wrong')

            allocate(fcsD(count(fcs /= Qface1)))
            fcsD = pack(fcs, fcs /= Qface1)
            fcs_newp = [fcsD, f2n, f4n, f5n]
            if (size(fcs_newp) /= 5) call gdErrorHandler('SplitQQQ: Something wrong')

            ! Face
            cv1 = cv
            s = c%faceP1%Get(cv1)
            call c%face%Replace(s, s + c%faceP2%Get(cv1) - 1, fcs_newp)
            call c%faceP2%Set(cv1, 5)
            range = (/ (i, i = cv1+1, c%faceP1%Size())/)
            call c%faceP1%SumMask(range, -1)

            ! Verts
            call c%vert%Replace(s, s + c%faceP2%Get(cv1) - 1, vxs_newp)
            call c%vertP2%Set(cv1, 5)
            call c%vertP1%SumMask(range, -1)

            ! Turn quad on the Qface1 side into a pent
            ! neig1 => become pent
            call grid%QuadToPent(neig1, Qface1, f1n, f2n, v1n)

            ! Triangle
            vertsT = [v1n, v2n, vxsT]
            facesT = [f1n, f3n, f5n]
            call grid%AddCell(facesT, vertsT, c%reg%Get(cv1), cv2)
            
            ! Remove Qface
            face_rem = [Qface1, Qface2]
            call grid%RemoveFaces(face_rem)

            ! Determine cflags
            cells = [cv1, cv2]
            call grid%DetermineCflags(cells)

            ! Split further in the radial direction (to reduce poloidal length)
            ! If there is still an aligned faces 
            if (size(fcs_al) /= 0) then
                ! Update counter
                if (fcs_al(1) .gt. Qface1) fcs_al(1) = fcs_al(1) - 1
                if (Qface2 .gt. Qface1) Qface2 = Qface2 - 1
                if (fcs_al(1) .gt. Qface2) fcs_al(1) = fcs_al(1) - 1

                ! Get cells across the aligned face
                cvs = GetFaceCellGA(c, fcs_al(1))
                if (cvs(1) /= cv1) then
                    neig = cvs(1)
                else if (cvs(2) /= cv1) then
                    neig = cvs(2)
                end if
                type2 = 'pol'
                ! call grid%SplitPQpol_rad(cv1, fcs_al, neig, v2n, type2)

            else 

 
                ! If there is no aligned face, so no neighbor there
                neig = 0
                type2 = 'pol'
                f1n = f1n - 2  !Update f1
                
                ! Get the face of the pentagon to split
                ! Update face number of face f1_n => new faces to high number, so will not be affect by removed faces
                fcs_cv1 = GetCellFaceGA(c, cv1)
                
                ! Get faces which have vxsQ1 as vertex
                vxs_fcs = [f%vert1%Get(fcs_cv1), f%vert2%Get(fcs_cv1)]
                query_faces = [(/(i, i = 1, size(fcs_cv1))/), (/(i, i = 1, size(fcs_cv1))/)]
                allocate(ind(count(vxs_fcs /= vxsQ1)))
                ind = pack(query_faces, vxs_fcs /= vxsQ1)
                fcs_vxsQ1 = fcs_cv1(ind)
                if (fcs_vxsQ1(1) /= f1n) then
                    split_face = fcs_vxsQ1(1) 
                else if (fcs_vxsQ1(2) /= f1n) then
                    split_face = fcs_vxsQ1(2)                     
                end if
                
                !call grid%SplitPQpol_rad(grid,cv1,split_face,neig,v2n,type2)            
            end if

            ! Extract fsVx from fsFc
            if (.not.options%slab) then
                call grid%GetFsVxFromFsFc(options)
            end if

            ! If pentagon are not allowed
            ! Make sure only the splitting in radial direction is continiued in this
            ! routine
            if (options%no_pents) then
                !call grid%PropagatePolSplitting(options)    
            end if            


            ! Do ordening of vert and face for nice plot
            if (options%debug) then
                if (.not.options%slab) then
                    allocate(cells_log(c%ntot))
                    cells_log = .true.
                    call grid%CheckVertOrder(is_ordered, cells_log)
                    call grid%ReOrderCellConn(is_ordered)
                else
                    !%Only order cells which are not in the cuts
                    !cells = 1:grid.cell.ntot
                    !cells = cells(logical(~ismember(cells,grid.cellsCut)));
                    !is_ordered = CheckVertOrder(grid,cells);
                    !grid = ReOrderCellConn(grid,is_ordered,cells); 
                end if
            end if




        end if

        end associate

    end subroutine 

    subroutine SplitPQ(grid, magneticField, cv, Qface, neig, common_vert, type, options)

        ! Description
        !============
        !  Splits a pentagonal cell with a quads as neighbor.
        ! This routine is used very often to try to optimize local implementation.

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)     :: grid
        type(MagneticFieldUDT), intent(in)  :: magneticField
        integer(I8), intent(in)             :: cv, Qface, neig, common_vert
        character(:), allocatable           :: type
        type(GAoptionsUDT), intent(in)      :: options

        ! Auxiliary
        integer(I8) :: face_old, v1n, f1n, f2n, f3n, cv1, cv2
        integer(I8), allocatable, dimension(:) :: perp_faces, vxs1, fcs1, cells, face_rem
        real(R8) :: psic
        character(:), allocatable :: typeV

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        psic = 0_R8
        if (.not.options%XpointSplitting) then
            if (type == 'pol') then
                typeV = 'pol'
            else if (type == 'rad') then
                typeV = 'psi'
                psic = v%psi%Get(common_vert)
            end if
        else 
            typeV = 'geometric'
        end if 

        ! Split faces
        call grid%SplitFace(magneticField, Qface, v1n, f1n, f2n, typeV, psic)

        ! Determine perpendicular faces
        call grid%DeterminePerpFacePent(cv, Qface, perp_faces)

        ! Make new face
        call grid%GetFaceNumber(v1n, common_vert, 3, f3n)

        if (all(f%aligned%Get(perp_faces) == 1)) call f%aligned%Set(f3n, 1)

        ! Add new vert to flux surface
        select case (type)
        case ('pol')
            allocate(vxs1(1))
            vxs1 = v1n
            fcs1 = [f1n, f2n]
            face_old = Qface
        case ('rad')
            allocate(fcs1(1))
            vxs1 = [v1n, common_vert]
            fcs1 = f3n
            face_old = 0
        end select
        call grid%AddVertToFsVx(vxs1, fcs1, face_old, type)

        ! Build new cells
        ! neig1 => become pent
        if (neig /= 0) call grid%QuadToPent(neig, Qface, f1n, f2n, v1n)

        ! Split pent in two quad
        call grid%SplitCenterPent(cv, Qface, v1n, common_vert, f3n, f1n, f2n, cv1, cv2)

        ! Remove Qface
        allocate(face_rem(1))
        face_rem = Qface  
        call grid%RemoveFaces(face_rem)

        ! Determine cflags
        cells = [cv1, cv2]
        call grid%DetermineCflags(cells)

        end associate

    end subroutine 

    subroutine SplitPT(grid, magneticField, cv, Tface, Tneig, common_vert, type, typeT)

        ! Description
        !============
        ! Splits a pentagonal cell with a triangle as neigbor in the poloidal or direction
        
        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)     :: grid
        type(MagneticFieldUDT), intent(in)  :: magneticField
        integer(I8), intent(in)             :: cv, Tface, Tneig, common_vert
        character(:), allocatable           :: type, typeT

        ! Auxiliary
        integer(I8) :: v1n, f1n, f2n, f3n, f4n, free_vert, face_old, cv1, cv2, Tneig1, Tneig2
        integer(I8), allocatable, dimension(:) ::  vxs1, fcs1, face_rem, perp_faces, cells
        character(:), allocatable :: typeV
        real(R8) :: psic

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        psic = 0
        if (type == 'pol') then
            typeV = 'pol'
        else if (type == 'rad') then
            typeV = 'psic'
            psic = v%psi%Get(common_vert)
        end if

        ! Split Tface
        call grid%SplitFace(magneticField, Tface, v1n, f1n, f2n, typeV, psic)

        ! Determine free_vert and make new face
        if (type == 'pol' .and. typeT == 'stacked') then
                call grid%DetermineFreeVertTria(Tneig, Tface, free_vert)
                call grid%GetFaceNumber(v1n, free_vert, 3, f4n)
        end if

        ! Make new face in pentagon
        call grid%GetFaceNumber(v1n, common_vert, 3, f3n)
        if (type == 'rad') then
            call grid%DeterminePerpFacePent(cv, Tface, perp_faces)
            if (all(f%aligned%Get(perp_faces) == 1)) call f%aligned%Set(f3n, 1)
        end if

        ! Add new vert to fsVx
        if (type == 'pol') then
            allocate(vxs1(1))
            vxs1 = v1n
            fcs1 = [f1n, f2n]
            face_old = Tface
        else if (type == 'rad') then
            allocate(fcs1(1))
            vxs1 = [v1n, common_vert]
            fcs1 = f3n
            face_old = 0
        end if
        call grid%AddVertToFsVx(vxs1, fcs1, face_old, type)

        ! Build new cells
        ! Neighboring triangle
        if (type == 'pol' .and. typeT == 'stacked') then
            call grid%SplitTriaStacked(Tneig, Tface, v1n, free_vert, Tneig1, Tneig2)
        else
            call grid%TriaToQuad(Tneig, Tface, f1n, f2n, v1n)
        end if

        ! Split pent in two quads
        call grid%SplitCenterPent(cv, Tface, v1n, common_vert, f3n, f1n, f2n, cv1, cv2)

        ! Remove Tface
        allocate(face_rem(1))
        face_rem = Tface
        call grid%RemoveFaces(face_rem)

        ! Determine cflags
        if (type == 'pol' .and. typeT == 'stacked') then
            cells = [Tneig1, Tneig2, cv1, cv2]
        else
            cells = [cv1, cv2]
        end if
        call grid%DetermineCflags(cells)

        end associate

    end subroutine

    !subroutine SplitPB(grid, magneticField, cv, Qface, common_vert, type)

        ! Description
        !============
        ! Splitting of pentagonal boundary cells with the internal
        ! cell as the splitted face. Qface is the face the should be splitted.
        ! Common_vert is the vertex of the pentagon that splits the internal radial
        ! face.

        ! Declare variables
        !==================
        ! Arguments
        !class(GAGridUDT), intent(inout)     :: grid
        !type(MagneticFieldUDT), intent(in)  :: magneticField
        !integer(I8), intent(in)             :: cv, Qface, common_vert
        !character(:), allocatable           :: type

        ! Auxiliary
        !integer(I8), allocatable :: perp_faces(:)
        !real(R8) :: psic
        !character(:), allocatable :: typeV

        ! Determine method
        !psic = 0
        !if (type == 'pol') then
        !    typeV = 'pol'
        !else if (type == 'rad') then
        !    typeV = 'psi'
        !    psic = v%psi%Get(common_vert)
        !end if

        ! Determine perp faces
        !call grid%DeterminePerpFacePent(cv, Qface, perp_faces)

        ! Split face



    !end subroutine
    
    subroutine SplitT4B(grid, qm, magneticField, cv, rface1, common_vert, type, options)

        ! Description
        !============
        !  Splitting of triangle4 boundary cell

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)         :: grid
        type(QualityMetricUDT), intent(inout)   :: qm
        type(MagneticFieldUDT), intent(in)      :: magneticField
        integer(I8), intent(in)                 :: cv, rface1, common_vert
        character(:), allocatable               :: type
        type(GAoptionsUDT), intent(in)          :: options

        ! Auxiliary
        integer(I8) :: v1n, f1n, f2n, cv1, cv2, f3n, face_old, b_fcs, free_vert
        integer(I8), allocatable, dimension(:) :: perp_faces, face_rem, cells, &
            fcs1, vxs1, b_fcsD, fcs_perp, fcs_al, fcs
        real(R8) :: psic
        character(:), allocatable :: typeV
        type(GAoptionsUDT) :: options2

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        select case (type)
        case ('rad')

            ! Split face
            psic = v%psi%Get(common_vert)
            typeV = 'psi'
            call grid%SplitFace(magneticField, rface1, v1n, f1n, f2n, typeV, psic)

            ! Determine perpendicular faces
            call grid%DeterminePerpFacePent(cv, rface1, perp_faces)

            ! Make new face
            call grid%GetFaceNumber(v1n, common_vert, 3, f3n)

            ! Add new vertex in flux surface
            allocate(fcs1(1))
            vxs1 = [v1n, common_vert]
            fcs1 = f3n
            face_old = 0
            call grid%AddVertToFsVx(vxs1, fcs1, face_old, type)

            ! Build new cells
            call grid%SplitT4(cv, rface1, v1n, common_vert, cv1, cv2)

            ! Remove face rface1
            allocate(face_rem(1))
            face_rem = rface1
            call grid%RemoveFaces(face_rem)

            ! Determine cflags
            cells = [cv1, cv2]
            call grid%DetermineCflags(cells)

            ! Do a poloidal split
            options2 = options
            options2%splittype = 'pol'
            options2%QTtype = 'pol-rad'
            qm%split_cv = cv2
            call grid%Splitting(qm, options2, magneticField)

        case ('pol')

            ! Get non-aligned faces
            fcs = GetCellFaceGA(c, cv)
            allocate(fcs_perp(count(f%aligned%Get(fcs) == 0)))
            fcs_perp = pack(fcs, f%aligned%Get(fcs) == 0)
            allocate(fcs_al(count(f%aligned%Get(fcs) == 1)))
            fcs_al = pack(fcs, f%aligned%Get(fcs) == 1)

            if (size(fcs_perp) /= 2) call gdErrorHandler('SplitT4B: something wrong')

            select case (options%typeT)
            case ('stacked')

                ! Get free vert
                free_vert = GetCommonVert(f, fcs_perp(1), fcs_perp(2))

                ! Make new face
                call grid%GetFaceNumber(free_vert, common_vert, 3, f1n)

                ! Build new cells
                call grid%SplitT4Stacked(cv, f1n, free_vert, fcs_perp, fcs_al, cv1, cv2)

                ! cflags
                cells = [cv1, cv2]
                call grid%DetermineCflags(cells)


            case ('cutcell')

                ! Get boundary faces
                allocate(b_fcsD(count(isBoundaryFaceGA(f, fcs))))
                b_fcsD = pack(fcs, isBoundaryFaceGA(f, fcs))
                b_fcs = b_fcsD(1)

                ! Split face
                typeV = 'pol'
                psic = 0_R8
                call grid%SplitFace(magneticField, b_fcs, v1n, f1n, f2n, typeV, psic)

                ! Make new face
                call grid%GetFaceNumber(v1n, common_vert, 3, f3n)

                ! Add new vert to flux surface
                allocate(vxs1(1))
                vxs1 = v1n
                fcs1 = [f1n, f2n]
                call grid%AddVertToFsVx(vxs1, fcs1, b_fcs, type)

                ! Build new cells
                call grid%SplitT4(cv, b_fcs, v1n, common_vert, cv1, cv2)

                ! Remove face
                allocate(face_rem(1))
                face_rem = b_fcs
                call grid%RemoveFaces(face_rem)

                ! Determine cflags
                cells = [cv1, cv2]
                call grid%DetermineCflags(cells)

            end select

        end select

        end associate

    end subroutine

    subroutine SplitT4Q(grid, magneticField, cv, rface1, common_vert, neig1, type)

        ! Description
        !============
        ! Radial splitting of triangle4 internal cell with a quad neighboring

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)     :: grid
        type(MagneticFieldUDT), intent(in)  :: magneticField
        integer(I8), intent(in)             :: cv, rface1, neig1, common_vert
        character(:), allocatable           :: type

        ! Auxiliary 
        integer(I8) :: face_old, v1n, f1n, f2n, f3n, cv1, cv2
        integer(I8), allocatable, dimension(:) :: vxs1, fcs1, perp_faces, &
            face_rem, cells  
        character(:), allocatable :: typeV
        real(R8) :: psic
         

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )
            
        psic = 0_R8
        if (type == 'pol') then
            typeV = 'pol'
        else if (type == 'rad') then
            typeV = 'psi'
            psic = v%psi%Get(common_vert)
        end if

        ! Determine perpendicular faces
        call grid%DeterminePerpFacePent(cv, rface1, perp_faces) ! cv is not a pent but works as well

        ! Split faces
        call grid%SplitFace(magneticField, rface1, v1n, f1n, f2n)

        ! Make new face
        call grid%GetFaceNumber(v1n, common_vert, 3, f3n)
        call f%aligned%Set(f3n, f%aligned%Get(perp_faces(1)))

        ! Add new vertices to flux surface
        select case (type)
        case ('pol')
            allocate(vxs1(1))
            vxs1 = v1n
            fcs1 = [f1n, f2n]
            face_old = rface1
        case ('rad')
            allocate(fcs1(1))
            vxs1 = [v1n, common_vert]
            fcs1 = f3n
            face_old = 0
        end select
        call grid%AddVertToFsVx(vxs1, fcs1, face_old, type)

        ! Build new cells
        call grid%SplitT4(cv, rface1, v1n, common_vert, cv1, cv2)

        if (neig1 /= 0) call grid%QuadToPent(neig1, rface1, f1n, f2n, v1n)

        ! Remove rface1
        allocate(face_rem(1))
        face_rem = rface1
        call grid%RemoveFaces(face_rem)

        ! Determine cflags
        cells = [cv1, cv2]
        call grid%DetermineCflags(cells)

        end associate

    end subroutine

    subroutine SplitFace(grid, magneticField, face_in, v1n, f1n, f2n, type, psic)

        ! Description
        !============
        ! Splitting a face into two providing two new faces and the new vertex

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)                     :: grid
        type(MagneticFieldUDT), intent(in)                  :: magneticField
        integer(I8), intent(in)                             :: face_in
        integer(I8), intent(out)                            :: v1n, f1n, f2n
        character(:), allocatable, optional, intent(in)     :: type
        real(R8), intent(in), optional                      :: psic

        ! Make new vertices
        if (present(type)) then
            call grid%AddVert(magneticField, face_in, v1n, type, psic)
        else 
            call grid%AddVert(magneticField, face_in, v1n)
        end if
        
        ! Make new faces
        call grid%GetFaceNumber(v1n, grid%face%vert1%Get(face_in),3, f1n)
        call grid%GetFaceNumber(v1n, grid%face%vert2%Get(face_in),3, f2n)
        call grid%face%aligned%Set(f1n, grid%face%aligned%Get(face_in))     
        call grid%face%aligned%Set(f2n, grid%face%aligned%Get(face_in))     
        call grid%face%label%Set(f1n, grid%face%label%Get(face_in))     
        call grid%face%label%Set(f2n, grid%face%label%Get(face_in)) 


    end subroutine

    subroutine SplitTriaStacked(grid, cv, Tface, new_v, free_vert, cv1, cv2)

        ! Description
        !============
        ! Splits the triangular cv where the aligned face is splitted in two
        ! new_v is the new vertex
        ! Tface is the original face which is split in half and is
        ! split by new_v.
        ! Free_vert is the vertex which is not part of the Tface.
        ! cv1 and cv2 are the new triangles

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout) :: grid
        integer(I8), intent(in)         :: cv, Tface, new_v, free_vert
        integer(I8), intent(out)        :: cv1, cv2 

        ! Auxiliary
        integer(I8) :: s, f1_T1, f2_T1, f3_T1, f1_T2, f2_T2, f3_T2
            
        integer(I8), allocatable :: verts_T1(:), verts_T2(:), &
            faces_T1(:), faces_T2(:)

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        ! The original triangle is split into two new triangles Tneig1 and Tneig2.

        ! cv => split in two trianlges
        cv1 = cv
        verts_T1 = [new_v , f%vert1%Get(Tface) , free_vert]
        s = c%vertP1%Get(cv1)
        call c%vert%Replace(s,s+c%vertP2%Get(cv1)-1, verts_T1 )

        call grid%GetFaceNumber(new_v,f%vert1%Get(Tface),1, f1_T1)
        call grid%GetFaceNumber(f%vert1%Get(Tface),free_vert,1, f2_T1)
        call grid%GetFaceNumber(free_vert,new_v,1, f3_T1)

        faces_T1 = [f1_T1 , f2_T1 , f3_T1]
        s = c%faceP1%Get(cv1)
        call c%face%Replace(s, s + c%faceP2%Get(cv1) - 1, faces_T1)

        ! Tneig2
        verts_T2 = [new_v , f%vert2%Get(Tface) , free_vert]

        call grid%GetFaceNumber(new_v,f%vert2%Get(Tface),1, f1_T2)
        call grid%GetFaceNumber(f%vert2%Get(Tface),free_vert, 1, f2_T2)
        call grid%GetFaceNumber(free_vert,new_v, 1, f3_T2)

        faces_T2 = [f1_T2 , f2_T2 , f3_T2]

        ! Add new cell Tneig2
        call grid%AddCell(faces_T2, verts_T2, c%reg%Get(cv1), cv2)

        end associate

    end subroutine

    subroutine SplitCenterTria(grid, cv, face, v1n, v2n, perp_face, cv1, cv2)

        ! Description
        !============
        ! Splits the center tria cv
        ! new_v1 and new_v2 are the new vertices
        ! face is one of the original two faces which are split in half and is split by new_v1.
        ! Perp_face is the face where the free_vert is not part of. 

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout) :: grid
        integer(I8), intent(in)         :: cv, face, v1n, v2n, perp_face
        integer(I8), intent(out)        :: cv1, cv2

        ! Auxiliary
        integer(I8) :: free_vert, s, bface_vert, f1_cv1, f2_cv1, f3_cv1, &
            f2_cv2, f3_cv2, f4_cv2, vertp
        integer(I8), allocatable :: verts_cv(:), verts_cv1(:), faces_cv1(:), &
            verts_cv2(:), faces_cv2(:)

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face &
            )

        ! Build the triangle
        !===================
        cv1 = cv
        verts_cv = GetCellVertGA(c, cv)

        ! Verts
        ! Get free vert
        call grid%DetermineFreeVertTria(cv,perp_face,free_vert)
        verts_cv1 = [v1n, v2n, free_vert]

        ! Faces
        call grid%GetFaceNumber(v1n, v2n, 1, f1_cv1)
        call grid%GetFaceNumber(v1n, free_vert, 1, f2_cv1)
        call grid%GetFaceNumber(v1n, free_vert, 1, f3_cv1)
        faces_cv1 = [f1_cv1, f2_cv1, f3_cv1]

        s = c%faceP1%Get(cv)         
        call grid%cell%vert%Replace(s, s + c%faceP1%Get(cv) - 1, verts_cv1)
        call grid%cell%face%Replace(s, s + c%faceP1%Get(cv) - 1, faces_cv1)
        
        ! Recalculate centroid
        call grid%CalcCentroidGA(cv1)

        ! cv2 => the quad
        ! Verts
        ! Get bface_vert which is not the free_vert
        if (f%vert1%Get(face) /= free_vert) then
            bface_vert = f%vert1%Get(face)
        else if (f%vert2%Get(face) /= free_vert) then
            bface_vert = f%vert2%Get(face)          
        end if 
        if (f%vert1%Get(perp_face) /= bface_vert) then
            vertp = f%vert1%Get(perp_face)
        else if (f%vert2%Get(perp_face) /= bface_vert) then
            vertp = f%vert2%Get(perp_face)
        end if

        verts_cv2 = [v1n, v2n, bface_vert, vertp]

        ! Face
        call grid%GetFaceNumber(v1n, bface_vert, 1, f2_cv2)
        call grid%GetFaceNumber(vertp, bface_vert, 1, f3_cv2)  
        call grid%GetFaceNumber(vertp, v2n, 1, f4_cv2)  
        
        faces_cv2 = [f1_cv1, f2_cv2, f3_cv2, f4_cv2]

        ! Add new cell cv2
        call grid%AddCell(faces_cv2, verts_cv2, c%reg%Get(cv1), cv2)

        end associate

    end subroutine

    subroutine SplitCenterQuad(grid, cv, face, v1n, v2n, cv1, cv2)

        ! Description
        !============
        ! Splits the center quad cv
        ! new_v1 and new_v2 are the new vertices
        ! face is one of the original two faces which are split in half and is
        ! split by new_v1.

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout) :: grid
        integer(I8), intent(in)         :: cv, face, v1n, v2n
        integer(I8), intent(out)        :: cv1, cv2

        ! Auxiliary
        integer(I8) :: s, pface, v1, f1_cv1, f2_cv1, f3_cv1, f4_cv1, vertp, &
            f2_cv2, f3_cv2, f4_cv2, v2
        integer(I8), allocatable, dimension(:) :: faces_cv, verts_cv, query_faces, fcs_a, &
            faces_cv1, verts_cv1, fcs_b, faces_cv2, verts_cv2

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face &
            )

        ! cv1 
        cv1 = cv
        v1 = f%vert1%Get(face)
        faces_cv = GetCellFaceGA(c, cv1)
        verts_cv = [ f%vert1%Get(faces_cv), f%vert2%Get(faces_cv)]
        query_faces = [faces_cv, faces_cv]
        allocate(fcs_a(count(verts_cv == v1)))
        fcs_a = pack(query_faces, verts_cv == v1)

        if (fcs_a(1) /= face) then
            pface = fcs_a(1)
        else if (fcs_a(2) /= face) then
            pface = fcs_a(2)
        end if

        if (f%vert1%Get(pface) /= v1) then
            vertp = f%vert1%Get(pface)
        else if (f%vert2%Get(pface) /= v1) then
            vertp = f%vert2%Get(pface)
        end if
        
        verts_cv1 = [v2n, v1n, v1, vertp]

        call grid%GetFaceNumber(v1n, v2n, 1, f1_cv1)
        call grid%GetFaceNumber(v1n, v1, 1, f2_cv1)
        call grid%GetFaceNumber(vertp, v1, 1, f3_cv1)
        call grid%GetFaceNumber(vertp, v2n, 1, f4_cv1)

        faces_cv1 = [f1_cv1, f2_cv1, f3_cv1, f4_cv1]

        ! Add in connectivity
        s = c%faceP1%Get(cv1)
        call c%face%Replace(s, s + c%faceP2%Get(cv1) - 1, faces_cv1)
        call c%vert%Replace(s, s + c%vertP2%Get(cv1) - 1, verts_cv1)

        ! Recalculate centroid
        call grid%CalcCentroidGA(cv1)

        ! Second cell cv2
        v2 = f%vert2%Get(face)
        allocate(fcs_b(count(verts_cv == v2)))
        fcs_b = pack(query_faces, verts_cv == v2)

        if (fcs_b(1) /= face) then
            pface = fcs_b(1)
        else if (fcs_b(2) /= face) then
            pface = fcs_b(2)
        end if

        if (f%vert1%Get(pface) /= v2) then
            vertp = f%vert1%Get(pface)
        else if (f%vert2%Get(pface) /= v2) then
            vertp = f%vert2%Get(pface)
        end if

        verts_cv2 = [v2n, v1n, v2, vertp]

        call grid%GetFaceNumber(v1n, v2, 1, f2_cv2)
        call grid%GetFaceNumber(vertp, v2, 1, f3_cv2)
        call grid%GetFaceNumber(vertp, v2n, 1, f4_cv2)

        faces_cv2 = [f1_cv1, f2_cv2, f3_cv2, f4_cv2]

        ! Add cell cv2
        call grid%AddCell(faces_cv2, verts_cv2, c%reg%Get(cv1), cv2)

        end associate

    end subroutine

    subroutine SplitT4(grid, cv, face, v1n, common_vert, cv1, cv2)

        ! Description
        !============
        ! Splits a triangle in half resulting in a quad and triangle
        ! The original cv number is used for the triangle, a new cell is made for
        ! the triangle

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)     :: grid
        integer(I8), intent(in)             :: cv, face, v1n, common_vert
        integer(I8), intent(out)            :: cv1, cv2

        ! Auxiliary
        integer(I8) :: i, s, f1, tria_vert, vert1, vert2, f1_cv1, f2_cv1, f3_cv1, &
            f2_cv2, f3_cv2, f4_cv2
        integer(I8), allocatable, dimension(:) :: verts_cv, test_verts, range, &
            verts_cv2_rest, faces_cv2, verts_cv2, faces_cv1, verts_cv1

        associate(&
            c => grid%cell, &
            f => grid%face &
            )

        ! Save verts_cv for later
        verts_cv = GetCellVertGA(c, cv)

        ! cv1 => triangle becomes other triangle
        cv1 = cv

        ! verts of triangle are v1_n, common_vert, ... and the vertex which is in
        !face and also has already a face together with the common vert 
        test_verts = [f%vert1%Get(face), f%vert2%Get(face)]
        call f%CheckFace(common_vert, test_verts(1), f1)

        if (f1 /= 0) then
            tria_vert = test_verts(1)
        else 
            tria_vert = test_verts(2)
        end if

        verts_cv1 = [common_vert, v1n, tria_vert]

        ! Verts
        s = c%vertP1%Get(cv1)
        call c%vert%Replace(s, s + c%vertP2%Get(cv1) - 1, verts_cv1)
        call c%vertP2%Set(cv1, size(verts_cv1))
        range = (/ (i, i = cv+1, c%vertP1%Size())/)
        call c%vertP1%SumMask(range, -1)

        ! Faces
        call grid%GetFaceNumber(common_vert, v1n, 1, f1_cv1)
        call grid%GetFaceNumber(v1n, tria_vert, 1, f2_cv1)
        call grid%GetFaceNumber(tria_vert, common_vert, 1, f3_cv1)

        faces_cv1 = [f1_cv1, f2_cv1, f3_cv1]
        call c%face%Replace(s, s + c%faceP2%Get(cv1) - 1, faces_cv1)
        call c%faceP2%Set(cv1, size(faces_cv1))
        call c%faceP1%SumMask(range, -1)


        ! cv2
        allocate(verts_cv2_rest(count(verts_cv /= tria_vert .and. verts_cv /= common_vert)))
        verts_cv2_rest = pack(verts_cv, verts_cv /= tria_vert .and. verts_cv /= common_vert)

        ! Get vert connect to v1n
        if (test_verts(1) /= tria_vert) then
            vert1 = test_verts(1)
        else if (test_verts(2) /= tria_vert) then
            vert1 = test_verts(2)
        end if
        if (verts_cv2_rest(1) /= common_vert) then
            vert2 = verts_cv2_rest(1)
        else if (verts_cv2_rest(2) /= common_vert) then
            vert2 = verts_cv2_rest(2)
        end if

        verts_cv2 = [v1n, common_vert, vert2, vert1]

        call grid%GetFaceNumber(common_vert, vert2, 1, f2_cv2)
        call grid%GetFaceNumber(vert2, vert1, 1, f3_cv2)
        call grid%GetFaceNumber(vert1, common_vert, 1, f4_cv2)

        faces_cv2 = [f1_cv1, f2_cv2, f3_cv2, f4_cv2]

        ! Add cell cv2
        call grid%AddCell(faces_cv2, verts_cv2, c%reg%Get(cv1), cv2)

        end associate
    end subroutine
    
    subroutine SplitT4Stacked(grid, cv, f1n, free_vert, fcs_perp, fcs_al, cv1, cv2)

        ! Description
        !============
        ! Splitting a 4 triangle poloidally with stacked method

        ! Declare variables
        !==================
        ! Arguments
        class(GAgridUDT), intent(inout) :: grid
        integer(I8), intent(in)         :: cv, f1n, free_vert, fcs_perp(:), fcs_al(:)
        integer(I8), intent(out)        :: cv1, cv2

        ! Auxiliary
        integer(I8) :: fcP1, fcA1, fcP2, fcA2, s, i
        integer(I8), allocatable :: fcs(:), range(:), faces_cv1(:), verts_cv1(:), &
            faces_cv2(:), verts_cv2(:)

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face &
            )

        fcs = GetCellFaceGA(c, cv)

        ! Group correctly
        if (HaveCommonVert(f, fcs_perp(1), fcs_al(1))) then
            fcP1 = fcs_perp(1)
            fcA1 = fcs_al(1)
            fcP2 = fcs_perp(2)
            fcA2 = fcs_al(2)
        else if (HaveCommonVert(f, fcs_perp(1), fcs_al(2))) then
            fcP1 = fcs_perp(1)
            fcA1 = fcs_al(2)
            fcP2 = fcs_perp(2)
            fcA2 = fcs_al(1)           
        end if


        ! cv1 => 4 triangle become real triangle
        cv1 = cv
        faces_cv1 = [ fcA1, fcP1, f1n]
        s = c%faceP1%Get(cv1)
        call c%face%Replace(s, s + c%faceP2%Get(cv1) - 1, faces_cv1)
        call c%faceP2%Set(cv1, size(faces_cv1))
        range = (/ (i, i = cv1+1, c%faceP1%Size() )/)
        call c%faceP1%SumMask(range, -1)

        verts_cv1 = [f%vert1%Get(fcA1), f%vert2%Get(fcA1), free_vert]
        call c%vert%Replace(s, s + c%vertP2%Get(cv1) - 1, verts_cv1)
        call c%vertP2%Set(cv1, size(verts_cv1))
        range = (/ (i, i = cv1+1, c%vertP1%Size() )/)
        call c%vertP1%SumMask(range, -1)

        ! cv2, new cell
        faces_cv2 = [fcA2, fcP2, f1n]
        verts_cv2 = [f%vert1%Get(fcA2), f%vert2%Get(fcA2), free_vert]

        ! Add new cell
        call grid%AddCell(faces_cv2, verts_cv2, c%reg%Get(cv1), cv2)

        end associate

    end subroutine

    subroutine TriaToQuad(grid, cv, face, new_f1, new_f2, new_v)

        ! Description
        !============
        ! Turn a triangle in a quadrilateral cell by splitting one face. cv is cell
        ! number. Face is the original face which is split in two
        ! new_f1 and new_f2 are the new faces
        ! new_v is the new vertex

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout) :: grid
        integer(I8), intent(in)         :: cv, face, new_f1, new_f2, new_v

        ! Auxiliary
        integer(I8) :: i, s
        integer(I8), allocatable :: fcs1D(:), fcs1(:), faces_cv(:), range(:), &
            vxs1D(:), verts_cv(:)

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face &
            )

        ! Get faces for the quadrilateral
        fcs1D = GetCellFaceGA(c, cv)
        allocate(fcs1(count(fcs1D /= face)))
        fcs1 = pack(fcs1D, fcs1D /= face)
        if (size(fcs1) /= 2) then
            call gdErrorHandler('TriaToQuad: face not is faces1')
        end if
        faces_cv = [fcs1, new_f1, new_f2]

        ! Replace faces
        s = c%faceP1%Get(cv)
        call c%face%Replace(s, s + c%faceP2%Get(cv) - 1,  faces_cv)
        call c%faceP2%Set(cv, size(faces_cv))
        range = (/(i, i = cv+1, c%faceP1%Size() )/)
        call c%faceP1%SumMask(range, 1)

        ! Vertices
        vxs1D = GetCellVertGA(c, cv)
        verts_cv = [vxs1D, new_v]

        call c%vert%Replace(s, s + c%vertP2%Get(cv) - 1,  verts_cv)
        call c%vertP2%Set(cv, size(verts_cv))
        call c%vertP1%SumMask(range, 1)
        
        ! To signal this cell (because should be removed)
        call c%cflags%Set(cv, 4) ! triangle with 4 vertices

        end associate

    end subroutine

    subroutine QuadToPent(grid, cv, ifc, new_f1, new_f2, new_v)

        ! Description
        !============
        ! Turns a quad in a pentagonal cell. cv is the cell number. 
        ! ifc is the original face which is split in two
        ! new_f1 and new_f2 are the new faces
        ! new_v is the new vertex

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout) :: grid
        integer(I8), intent(in)         :: cv, ifc, new_f1, new_f2, new_v

        ! Auxiliary
        integer(I8) :: s, i, ind
        integer(I8), allocatable :: faces1D(:), faces1(:), verts1D(:), verts1(:), &
            range(:)

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face &
            )

        ! Construct thepentagon
        ! Faces
        faces1D = GetCellFaceGA(c, cv)
        ind = findloc(faces1D, ifc, 1)
        if (ind == 0) call gdErrorHandler('QuadToPent: old face not in quads faces')
        faces1D(ind) = new_f1
        faces1 = [faces1D, new_f2]

        s = c%faceP1%Get(cv)
        call c%face%Replace(s, s + c%faceP2%Get(cv) - 1, faces1)
        call c%faceP2%Set(cv, size(faces1)) ! should be five
        range = (/(i, i = cv+1, c%faceP1%Size() )/)
        call c%faceP1%SumMask(range, 1)

        ! Vertices
        verts1D = GetCellVertGA(c, cv)
        verts1 = [verts1D, new_v]

        s = c%vertP1%Get(cv)
        call c%vert%Replace(s, s + c%vertP2%Get(cv) - 1, verts1)
        call c%vertP2%Set(cv, size(verts1))
        call c%vertP1%SumMask(range, 1)
        
        end associate

    end subroutine

    subroutine SplitCenterPent(grid, cv, face, new_v1, common_vert, fc_new, f1_new, f2_new, cv1, cv2)

        ! Description
        !============
        ! Splits the center pentagonal cv
        ! new_v1 is the new vertex
        ! face is  the original face which is split in half and is split by new_v1
        ! The common vert is the vertex of the pentagon that splits a face is two.
        ! fc_new is the new face created between the hanging node and new_v1.
        ! f1_new and f2_new are the new faces created by splitting 'face'.

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout) :: grid
        integer(I8), intent(in) :: cv, face, new_v1, common_vert, &
            fc_new, f1_new, f2_new
        integer(I8), intent(out) :: cv1, cv2

        ! Auxiliary
        integer(I8) :: i, s, v1, v2, f1_cv1, f2_cv1, f3_cv1, f4_cv1, &
            f1_cv2, f2_cv2, f3_cv2, f4_cv2, vertp
        integer(I8), allocatable, dimension(:) :: verts_faces_cv, query_faces, &
            ind, pface, verts_pface, faces_cv, verts_cv1, range, query5, &
            faces_cv1, ind1, ind2, ind3, faces_cv2, verts_cv2
        logical, allocatable, dimension(:) :: log, log1, log2, log3, log4, log5, log6

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face &
            )

        ! cv1 => pent becomes quad, on the side of grid%face%vert1%Get(face)
        cv1 = cv
        faces_cv = GetCellFaceGA(c, cv)
        allocate(verts_cv1(4))
        verts_cv1 = 0
        v1 = f%vert1%Get(face)
        verts_cv1(1:3) = [common_vert, new_v1, v1]

        ! Get faces where grid.face.vert(face,1) is vertex of to find all correct vertices
        verts_faces_cv = [ f%vert1%Get(faces_cv), f%vert2%Get(faces_cv) ]
        query_faces = [ faces_cv, faces_cv ]
        
        allocate(ind(count(verts_faces_cv == v1)))
        ind = pack(query_faces, verts_faces_cv == v1)

        allocate(pface(count(ind /= face)))
        pface = pack(ind, ind /= face)
        verts_pface = [f%vert1%Get(pface(1)), f%vert2%Get(pface(1))]
        deallocate(pface)
        if (verts_pface(1) /= v1) then
            verts_cv1(4) = verts_pface(1)
        else if (verts_pface(2) /= v1) then
            verts_cv1(4) = verts_pface(2)
        end if
        vertp = verts_cv1(4)

        ! Put verts in c%vert
        s = c%vertP1%Get(cv1)
        call c%vert%Replace(s, s+c%vertP2%Get(cv1)-1, verts_cv1)
        call c%vertP2%Set(cv1, size(verts_cv1))
        range = (/ (i, i = cv1+1, c%vertP1%Size())/)
        call c%vertP1%SumMask(range, -1)

        ! Get face for cv1
        f1_cv1 = fc_new
        f2_cv1 = f1_new

        ! Get third and fourth face
        query5 = (/(i, i = 1, 5)/)
        log1 = (verts_faces_cv == vertp)
        log2 = (log1(1:5) .or. log1(6:10))
        log1 = (verts_faces_cv == v1)
        log3 = (log1(1:5) .or. log1(6:10))
        log = (log2 .and. log3)
        allocate(ind1(count(log)))
        ind1 = pack(query5, log)
        f3_cv1 = faces_cv(ind1(1))
        deallocate(ind1)

        log1 = (verts_faces_cv ==  common_vert)
        log4 = (log1(1:5) .or. log1(6:10))
        log = (log2 .and. log4)
        allocate(ind2(count(log)))
        ind2 = pack(query5, log)
        f4_cv1 = faces_cv(ind2(1))
        deallocate(ind2)

        faces_cv1 = [f1_cv1, f2_cv1, f3_cv1, f4_cv1]

        ! Put faces in c%face
        s = c%faceP1%Get(cv1)
        call c%face%Replace(s, s+c%faceP2%Get(cv1)-1, faces_cv1)
        call c%faceP2%Set(cv1, size(faces_cv1))
        call c%faceP1%SumMask(range, -1)

        ! Recalculate centroid
        call grid%CalcCentroidGA(cv1)

        ! cv2 => new quad, on the side of f%vert2%Get(face)
        !==================================================
        v2 = f%vert2%Get(face)
        allocate(verts_cv2(4))
        verts_cv2 = 0
        verts_cv2(1:3) = [common_vert, new_v1, v2]

        allocate(ind3(count(verts_faces_cv == v2)))
        ind3 = pack(query_faces, verts_faces_cv == v2)

        allocate(pface(count(ind3 /= face)))
        pface = pack(ind3, ind3 /= face)
        verts_pface = [f%vert1%Get(pface(1)), f%vert2%Get(pface(1))]
        if (verts_pface(1) /= v2) then
            verts_cv2(4) = verts_pface(1)
        else if (verts_pface(2) /= v2) then
            verts_cv2(4) = verts_pface(2)
        end if 
        vertp = verts_cv2(4)
        
        ! Faces
        f1_cv2 = fc_new
        f2_cv2 = f2_new
        log1 = (verts_faces_cv == vertp)
        log5 = (log1(1:5) .or. log1(6:10))        
        log1 = (verts_faces_cv == v2)
        log6 = (log1(1:5) .or. log1(6:10))
        log = (log5 .and. log6)
        allocate(ind1(count(log)))
        ind1 = pack(query5, log)
        f3_cv2 = faces_cv(ind1(1))

        log = (log5 .and. log4)
        allocate(ind2(count(log)))
        ind2 = pack(query5, log)
        f4_cv2 = faces_cv(ind2(1))

        faces_cv2 = [f1_cv2, f2_cv2, f3_cv2, f4_cv2]

        ! Add new cv2
        call grid%AddCell(faces_cv2, verts_cv2, c%reg%Get(cv1), cv2)

        end associate

    end subroutine

    subroutine DeterminePerpFaceTria(grid, cv, face1, face2, perp_face)

        ! Description
        !============
        ! Get the thrid face of a triangle not equal to face1 and face2.
        
        ! Declare variables
        !==================
        ! Argument
        class(GAGridUDT), intent(in)    :: grid
        integer(I8), intent(in)         :: cv, face1, face2
        integer(I8), intent(out)        :: perp_face

        ! Auxiliary
        integer(I8), allocatable :: fcs(:), perp_faces(:)


        ! Get perpendicular face
        fcs = GetCellFaceGA(grid%cell, cv)
        allocate(perp_faces(count(fcs /= face1 .and. fcs /= face2)))
        perp_faces = pack(fcs, fcs /= face1 .and. fcs /= face2)
        perp_face = perp_faces(1) 

    end subroutine

    subroutine DeterminePerpFaceQuad(grid, cv, face1, face2, perp_faces)

        ! Description
        !============
        ! Get the perpendicular faces from a quad by eliminating radial faces face1 and face2

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(in)            :: grid
        integer(I8), intent(in)                 :: cv, face1, face2
        integer(I8), allocatable, intent(out)   :: perp_faces(:)

        ! Auxiliary
        integer(I8), allocatable :: fcs(:)

        ! Get perpendicular faces
        fcs = GetCellFaceGA(grid%cell, cv)
        allocate(perp_faces(count(fcs /= face1 .and. fcs /= face2)))
        perp_faces = pack(fcs, fcs /= face1 .and. fcs /= face2)
        
    end subroutine

    subroutine DeterminePerpFacePent(grid, cv, Qface, perp_faces)

        ! Description
        !============
        ! Get the perpendicular faces from a pentagon by eliminating radial faces face1 and face2

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(in)            :: grid
        integer(I8), intent(in)                 :: cv, Qface
        integer(I8), allocatable, intent(out)   :: perp_faces(:)

        ! Auxiliary
        integer(I8) :: n
        integer(I8), allocatable :: fcs(:)

        ! Get perpendicular faces
        fcs = GetCellFaceGA(grid%cell, cv)
        n = 1
        if (grid%face%aligned%Get(Qface) == 1)  n = 0
        allocate(perp_faces(count(grid%face%aligned%Get(fcs) == n)))
        perp_faces = pack(fcs, grid%face%aligned%Get(fcs) == n)

    end subroutine

    subroutine DetermineFreeVertTria(grid, cv, tface, free_vert)

        ! Description
        !============
        ! Determine the vertex of a triangle that is nog in the aligned face tface

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(in) :: grid
        integer(I8), intent(in)      :: cv, tface
        integer(I8), intent(out)     :: free_vert

        ! Auxiliary
        integer(I8) :: i
        integer(I8), allocatable :: vxs(:)

        vxs = GetCellVertGA(grid%cell, cv)
        do i = 1, size(vxs)
            if (vxs(i) /= grid%face%vert1%Get(tface) .and. vxs(i) /= grid%face%vert2%Get(tface)) then
                free_vert = vxs(i)
                exit
            end if
        end do
        
    end subroutine

    subroutine DetermineHangingNodePent(grid, cv, faces, type, common_vert, qface, rfaces)

        ! Description
        !============
        ! Determine the hanging node of a pentagon

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(in)            :: grid
        integer(I8), intent(in)                 :: cv, faces(:)
        character(:), allocatable, intent(in)   :: type
        integer(I8), intent(out)                :: common_vert, qface
        integer(I8), allocatable, intent(out)   :: rfaces(:)    

        ! Auxiliary
        integer(I8) :: i, iv, rface1, rface2, rface3, nr, &
                vx, indmin, pfaces(5), counter
        integer(I8), allocatable, dimension(:) :: vxs, rfaces_na, &
            vxs_faces1, vxs_faces2, vxs_faces, indf, query_facesD, query_faces
        real(R8), allocatable, dimension(:) :: fcX, fcY, sin
        real(R8) :: vec_vf1_x, vec_vf1_y, vec_vf2_x, vec_vf2_y

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        select case (type)

        case ('pol')

            ! ASSUMPTION: SPLITTED FACES IS A RADIAL FACE
            ! Get radial faces
            allocate(rfaces(count(f%aligned%Get(faces) == 1)))
            rfaces = pack(faces, f%aligned%Get(faces) == 1)
            rface1 = rfaces(1)
            rface2 = rfaces(2)
            nr = size(rfaces)
            if (nr == 3) then
                rface3 = rfaces(3)
            else if (nr .gt. 3) then
                call gdErrorHandler('DetermineHangingNodePent: too many radial faces in a pentagon')
            end if

        case ('rad')

            ! Get poloidal faces
            allocate(rfaces_na(count(f%aligned%Get(faces) == 0)))
            rfaces_na = pack(faces, f%aligned%Get(faces) == 0)

            ! Normal procedure
            allocate(rfaces_na(count(.not.isBoundaryFaceGA(f, faces))))
            rfaces = pack(rfaces_na, .not.isBoundaryFaceGA(f, faces))
            rface1 = rfaces(1)
            rface2 = rfaces(2)          
            nr = size(rfaces)
            if (c%vertP2%Get(cv) == 5) then
                if (size(rfaces_na) == 4) nr = 4   ! Trapezoid with splitted poloidal face
            end if

            if (nr == 3) then
                rface3 = rfaces(3)
            else if (nr .gt. 3) then

                ! Compute angle of check if it is a trapezoidal pentagon
                vxs = GetCellVertGA(c, cv)
                vxs_faces1 = f%vert1%Get(faces)
                vxs_faces2 = f%vert2%Get(faces)
                fcX = 0.5_R8 * (v%x%Get(vxs_faces1) + v%x%Get(vxs_faces2))
                fcY = 0.5_R8 * (v%y%Get(vxs_faces1) + v%y%Get(vxs_faces2))

                ! Get faces of vert
                vxs_faces = [vxs_faces1, vxs_faces2]
                query_facesD = (/ (i, i = 1, size(faces))/)
                query_faces = [ query_faces, query_faces]
                ! MAKE A QUERY
                allocate(sin(size(vxs)))
                
                sin = 0
                do i = 1, size(vxs)
                    iv = vxs(i)
                    allocate(indf(count(vxs_faces == iv)))
                    indf = pack(query_faces, vxs_faces == iv)
                    if (size(indf) .lt. 2) call gdErrorHandler('DetermineHangingNodePent: vertex only once in faces of cells')

                    vec_vf1_x = fcX(indf(1)) - v%x%Get(iv)
                    vec_vf1_y = fcY(indf(1)) - v%y%Get(iv)
                    vec_vf2_x = fcX(indf(2)) - v%x%Get(iv)
                    vec_vf2_y = fcY(indf(2)) - v%y%Get(iv)

                    ! Calculate angles (sin = |a x b| / norm(a)*norm(b))
                    ! with |a x b | = ax*by - bx*ay
                    sin(i) = (vec_vf1_x*vec_vf2_y - vec_vf2_x*vec_vf1_y) &
                             / (Norm(vec_vf1_x,vec_vf1_y)*Norm(vec_vf2_x,vec_vf2_y))
                    
                    ! Housekeeping
                    deallocate(indf)

                end do

                ! Common vert
                indmin = minloc(abs(sin),1)
                vx = vxs(indmin)

                ! Faces of common vert
                allocate(indf(count(vxs_faces == vx)))
                indf = pack(query_faces, vxs_faces == vx)
                rface1 = faces(indf(1))
                rface2 = faces(indf(2))
                nr = 2

            end if

        end select

        ! Radial faces with common vertex are faces on the side of the two quads
        if (nr == 2) then

            common_vert = GetCommonVert(f, rface1, rface2)

            ! Determine Qface = face without common vert
            pfaces = 0
            counter = 0
            do i = 1, size(faces)
                if (faces(i) /= rface1 .and. faces(i) /= rface2) then
                    counter = counter + 1
                    pfaces(counter) = faces(i)
                end if
            end do

            do i = 1, counter
                if (.not.HaveCommonVert(f, pfaces(i), rface1) .and. .not.HaveCommonVert(f, pfaces(i), rface2)) then
                    qface = pfaces(i)
                end if
            end do

        else if (nr == 3) then

            if (HaveCommonVert(f, rface1, rface2)) then

                qface = rface3
                common_vert = GetCommonVert(f, rface1, rface2)

            else if (HaveCommonVert(f, rface1, rface3)) then

                qface = rface2
                common_vert = GetCommonVert(f, rface1, rface3)

            else 

                qface = rface1
                common_vert = GetCommonVert(f, rface2, rface3)

            end if
        end if

        end associate

    end subroutine

    ! Boundary layer grid
    !====================
    subroutine BoundaryLayerGrid(grid, magneticField, options)

        ! Description
        !============
        ! 

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)     :: grid
        type(MagneticFieldUDT), intent(in)  :: magneticField
        type(GAoptionsUDT), intent(in)      :: options
  
    end subroutine

    !------------------------------------------------------------------!
    !                   GAGRID INFORMATION UTILITY                     !
    !------------------------------------------------------------------! 

    subroutine GetForbiddenMergeFaces(grid, forbidden_fcs)

        ! Description
        !============
        ! Determine forbidden merging faces
        ! The forbidden merging faces are
        ! - the separatrix faces
        ! - the faces of the core cut
        ! - boundary faces
        ! Important note: faces which are a region boundary can also not be used
        ! for merging. So always check whether the regions of the cells of a face
        ! are from the same region.

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)         :: grid
        integer(I8), allocatable, intent(out)   :: forbidden_fcs(:)

        ! Auxiliary
        integer(I8) :: i, indFc(grid%face%ntot), nb, nc, nr
        integer(I8), allocatable, dimension(:) :: cvLookUp, b_faces, tang_points, &
            rad_faces, fcsAll, core_faces
        
        ! Precompute
        cvLookUp = GetCvLookUpGA(grid%cell)

        ! 1. faces iFs of separatrix
        ! [~,~,grid] = give_iFs_sep(grid,grid.n_sep,cv_look_up,grid.iFs_sep);
        ! s = grid.fs.faceP(grid.iFs_sep,1);
        ! faces_sep = grid.fs.face(s:s+grid.fs.faceP(grid.iFs_sep,2)-1);

        ! 2. cuts around xpoint
        call grid%GetCutsXpoints(cvLookUp, core_faces)

        ! 3. boundary faces
        indFc = (/ (i, i = 1,grid%face%ntot) /)
        allocate(b_faces(count(grid%face%label%Get() /= 0)))
        b_faces = pack(indFc, grid%face%label%Get() /= 0);

        ! 4. region boundary - too expensive to do with every merge proposition - do in the determinecase
        !     cv_look_up = give_cv_look_up(grid);
        !     int_faces = faces(grid.fcLbl == 0);
        !    face_regions = zeros(round(grid.nFc/2),1); %over estimation
        !     counter = 0;
        !     for i = 1:length(int_faces)
        !         iFc = int_faces(i);
        !         cvs = give_cells_of_face(iFc,grid.cvFc,grid.cvFcP,cv_look_up);
        !        if grid.cvReg(cvs(1)) ~= grid.cvReg(cvs(2))
        !             counter = counter + 1;
        !             face_regions(counter) = iFc;
        !         end
        !     end
        !     face_regions = face_regions(1:counter);

        ! 5. tangency points - radial line faces
        call grid%GetTangencyPoints(b_faces, tang_points)
        call grid%GetRadLineFaces(tang_points,rad_faces)
    
        ! Compose forbidden faces array
        nc = size(core_faces)
        nb = size(b_faces)
        nr = size(rad_faces)
        allocate(fcsAll(nc+nb+nr))
        fcsAll(1:nc) = core_faces
        fcsAll(nc+1:nc+nb) = b_faces
        fcsAll(nc+nb+1:nc+nb+nr) = rad_faces
        call Unique(fcsAll,forbidden_fcs)
    

    end subroutine

    subroutine GetCutsXpoints(grid, cvLookUp, fcs)

        ! Description
        !============
        ! Give the poloidal and radial faces out of an Xpoint which should not be used to merge

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)         :: grid
        integer(I8), allocatable, intent(in)    :: cvLookUp(:)
        integer(I8), allocatable, intent(out)   :: fcs(:)

        ! Auxiliary
        integer(I8) :: fcsD(grid%face%ntot), counterf, i, &
            counter1, counter2, fcs1(grid%face%ntot), fcs2(grid%face%ntot)
        logical :: b_flag(grid%face%ntot), in_flag1(grid%face%ntot), &
            in_flag2(grid%face%ntot)


        ! Initialize
        call grid%GiveXpoints(.true.,cvLookUp)
        fcsD = 0
        counterf = 0
        counter1 = 0
        counter2 = 0
        b_flag = (grid%face%label%Get() == 0)
        in_flag1 = .false.
        in_flag2 = .false.
        fcs1 = 0
        fcs2 = 0

        ! Loop over the Xpoints
        do i = 1, grid%data%nxp

            ! Get the poloidal and radial faces
            call grid%RecursiveGridMarching(grid%data%xpointID(i),fcs1,0,counter1,b_flag,in_flag1)
            call grid%RecursiveGridMarching(grid%data%xpointID(i),fcs2,1,counter2,b_flag,in_flag2)

            ! Add fcs1
            fcsD(counterf+1:counterf+counter1) = fcs1(1:counter1)
            counterf = counterf + counter1

            ! Add fcs2
            fcsD(counterf+1:counterf+counter2) = fcs2(1:counter2)
            counterf = counterf + counter2

        end do

        ! Trim
        fcs = fcsD(1:counterf)

    end subroutine

    subroutine GetTangencyPoints(grid, b_Faces, tang_points)

        ! Description
        !============
        ! Give tangency points of a grid

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(in)            :: grid
        integer(I8),  intent(in)                :: b_faces(:)
        integer(I8), allocatable, intent(out)   :: tang_points(:)

        ! Auxiliary
        integer(I8) :: i, j, counter, nbv, iv, ifs, nfa
        integer(I8), allocatable, dimension(:) :: vxsB, cvLookUp, fsvLookUp, &
            tf, tang_pointsD, tf_aligned, cvs, int_face


        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face &
            )
        
        ! Initialize 
        allocate(tang_pointsD(nint(grid%vert%ntot/10.0_R8)))
        tang_pointsD = 0
        counter = 0

        ! Get boundary faces locally - check if necessary TODO

        
        ! Get boundary vertices
        vxsB = GetVxsFromFcsGA(f, b_faces)
        nbv = size(vxsB)

        ! Precompute
        cvLookUp = GetCvLookUpGA(c)
        fsvLookUp = GetFsvLookUpGA(grid%data%fluxdata)

        ! Loop over boundary vertices
        do i = 1, nbv
            
            ! Get vertex number
            iv = vxsB(i)

            ! Get flux surface
            ifs = GetVertFsvGA(grid%data%fluxdata, iv, fsvLookUp)

            if (ifs /= 0) then

                ! Get faces of vertex
                tf = GetVertFaceGA(f, iv)

                ! Get aligned faces
                nfa = count(f%aligned%Get(tf) == 1)
                allocate(tf_aligned(nfa))
                tf_aligned = pack(tf, f%aligned%Get(tf) == 1)

                if (nfa == 2) then ! Otherwise just end point of surface

                    ! Check whether faces are external of internal
                    allocate(int_face(nfa))
                    int_face = 0

                    do j = 1, nfa

                        cvs = GetFaceCellGA(c, tf_aligned(j), cvLookUp)
                        if (size(cvs) == 2) then

                            int_face(j) = 1

                        end if

                    end do

                    ! Add tangency point
                    if (sum(int_face) .gt. 0) then

                        counter = counter + 1
                        tang_pointsD(counter) = iv

                    end if

                    ! Housekeeping
                    deallocate(int_face)

                end if

                ! Housekeeping
                deallocate(tf_aligned)


            end if

        end do 

        ! Trim
        tang_points = tang_pointsD(1:counter)

        end associate

    end subroutine

    subroutine GetRadLineFaces(grid, verts, faces)

        ! Description
        !============
        ! Gives the faces of radial lines of the give vertices

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(in)            :: grid
        integer(I8), intent(in)                 :: verts(:)
        integer(I8), allocatable, intent(out)   :: faces(:)

        ! Auxiliary
        integer(I8) :: i, nv, faces1(grid%face%ntot), counter
        logical :: b_flag(grid%face%ntot), in_flag(grid%face%ntot)

        ! Initialize
        faces1 = 0
        counter = 0
        nv = size(verts)
        b_flag = (grid%face%label%Get() /= 0)
        in_flag = .false.

        ! Loop over vertices
        do i = 1, nv

            call grid%RecursiveGridMarching(verts(i), faces1, 0, counter, b_flag, in_flag)

        end do 

        ! Trim
        faces = faces1(1:counter)


    end subroutine

    recursive subroutine RecursiveGridMarching(grid, iv, faces, aligned, counter, b_flag, in_flag)

        ! Description
        !============
        ! Recursive Grid Marching starting from a vertex.
        ! Input:
        ! - starting vertex
        ! - grid
        ! - faces: should be empty to start
        ! - aligned: poloidal marching => aligned = 0
        !            radial marching   => aligned = 1
        ! - counter: should be 0 to start
        ! - b_flag = logical for boundary faces (1 == boundary face)

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(in)    :: grid
        integer(I8), intent(in)         :: iv, aligned
        integer(I8), intent(inout)      :: counter 
        integer(I8), intent(inout)      :: faces(grid%face%ntot)
        logical, intent(in)             :: b_flag(grid%face%ntot)
        logical, intent(inout)          :: in_flag(grid%face%ntot)

        ! Auxiliary
        integer(I8) :: i, nf, vs(1:2), iv_next
        integer(I8), allocatable :: fcs1(:), fcs2(:), fcs3(:), fcs4(:)

        ! Give faces
        fcs1 = GetVertFaceGA(grid%face, iv)

        ! Eliminate faces based on marching direction - TODO replace all pack by 1 pack by combining the criteria
        allocate(fcs2(count(grid%face%aligned%Get(fcs1) == aligned)))
        fcs2 = pack(fcs1,grid%face%aligned%Get(fcs1) == aligned)

        ! Eliminate boundary faces
        allocate(fcs3(count(.not.b_flag(fcs2))))
        fcs3 = pack(fcs2,.not.b_flag(fcs2))

        ! Eliminate faces which are in array faces
        allocate(fcs4(count(.not.in_flag(fcs3))))
        fcs4 = pack(fcs3,.not.in_flag(fcs3))

        ! Add remaining faces - if fcs would be empty, RecursiveGridMarching just
        ! returns faces and counter
        nf = size(fcs4)
        do i = 1, nf
            counter = counter + 1
            faces(counter) = fcs4(i)
            in_flag(fcs4(i)) = .true.

            ! Select new vertex
            vs(1) = grid%face%vert1%Get(fcs4(i))
            vs(2) = grid%face%vert2%Get(fcs4(i))

            if (vs(1) /= iv) then
                iv_next = vs(1)
            else if (vs(2) /= iv) then
                iv_next = vs(2)
            else

                call gdErrorHandler('RecursiveMarching: something went wrong')

            end if

            ! Recursive Marching on new vertex
            call grid%RecursiveGridMarching(iv_next, faces, aligned, counter, b_flag, in_flag)


        end do


    end subroutine

    subroutine AreaConstraintPents(grid, options, threshold, cells)

        ! Description
        !============
        ! Give cells which are not allow to be a pentagon and should be selected for
        ! splitting/merging

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT)                        :: grid
        type(GAoptionsUDT)                      :: options
        real(R8)                                :: threshold
        integer(I8), allocatable, intent(out)   :: cells(:)

        ! Auxiliary
        integer(I8) :: i
        integer(I8), allocatable :: indCv(:)
        real(R8) ::  val(grid%cell%ntot)

        ! Initiliaze
        indCv = (/ (i, i = 1, grid%cell%ntot)/)

        select case (options%no_pents_area_type)
            case ('coordinates')

                allocate(cells(count(grid%cell%x%Get() .lt. options%no_pents_area_maxR &
                        .and. grid%cell%x%Get() .lt. options%no_pents_area_minR &
                        .and. grid%cell%y%Get() .lt. options%no_pents_area_maxZ &
                        .and. grid%cell%y%Get() .lt. options%no_pents_area_minZ)))

                cells = pack(indCv, grid%cell%x%Get() .lt. options%no_pents_area_maxR &
                        .and. grid%cell%x%Get() .lt. options%no_pents_area_minR &
                        .and. grid%cell%y%Get() .lt. options%no_pents_area_maxZ &
                        .and. grid%cell%y%Get() .lt. options%no_pents_area_minZ)

            case ('dist_function')
                
                call grid%fun%distr%Evaluate(grid%cell%x%Get(),grid%cell%y%Get(), val)
                allocate(cells(count(val > threshold)))
                cells = pack(indCv, val > threshold)

            case default

                call gdErrorHandler('AreaConstraintPents: method not implemented')

        end select

        if (size(cells) == 0) call gdErrorHandler('AreaConstraintPents: no cells in area')    

    end subroutine

    subroutine DetermineCflags(grid, cells, b_flag)

        !  Description
        !=============
        ! Computes the cflags for cells. 
        ! A speed-up method is to pass precomputed data b_flag which is an logical array of
        ! all faces indicating whether a face is a boundary face.

        ! Declare variables
        !==================
        ! Argument
        class(GAGridUDT), intent(inout) :: grid
        integer(I8), intent(in)         :: cells(:)
        logical , optional              :: b_flag(grid%face%ntot)

        ! Auxiliary
        integer(I8) :: i, j, nb, indFc(grid%face%ntot)
        integer(I8), allocatable :: b_faces(:), cvs(:), cvLookUp(:), ar(:)

        ! Initialize - TODO - method without using face label - check if necessary
        if (.not.present(b_flag)) then

            ! Loop over cells
            do i = 1, size(cells)

                if (isBoundaryCellGA(grid,cells(i))) then
                    call grid%cell%cflags%Set(cells(i), 3)
                else
                    call grid%cell%cflags%Set(cells(i), 1)
                end if

            end do

        else 

            ! Initiliaze
            cvLookUp = GetCvLookUpGA(grid%cell)
            ar = (/ (1, j = 1, size(cells))/)
            call grid%cell%cflags%SetMultipleElements(cells,ar)
            indFc = (/ (i, i= 1, grid%face%ntot)/)
            nb = count(b_flag)
            allocate(b_faces(nb))
            b_faces = pack(indFc, b_flag)

            ! Loop over boundary cells
            do i = 1, nb
                cvs = GetFaceCellGA(grid%cell, b_faces(i), cvLookUp)
                ar = (/ (3, j = 1, size(cvs))/)
                call grid%cell%cflags%SetMultipleElements(cvs, ar)
            end do

        end if



    end subroutine

    subroutine GetXcellsGeometric(grid, Xcells)
        ! Description
        !============
        ! Determine the cells which should be splitted via an geometric method to
        ! avoid problems near the X-point due to noisy magnetic field.
        ! For the moment just cells around the Xpoints

        class(GAGridUDT) :: grid
        integer(I8) :: counter, i, XcellsD(grid%cell%ntot)
        integer(I8), allocatable :: Xcells(:), cvLookUp(:), cvs(:) 

        cvLookUp = GetCvLookUpGA(grid%cell) 

        ! Cells directly around the Xpoint
        do i = 1, grid%data%nxp
            cvs = GetVertCellGA(grid%cell, grid%data%xpointID(i), cvLookUp)
            XcellsD(counter+1:counter+size(cvs)) = cvs
        end do 

        ! Trim
        Xcells = XcellsD(1:counter)

    end subroutine

    !------------------------------------------------------------------!
    !                       GAGRID OPERATIONS                          !
    !------------------------------------------------------------------!  

    subroutine RemoveCells(grid, cells)
        
        ! Description
        !============
        ! Remove cells from the grid data
        ! But does not change the geometry of neighboring cells!
        ! Therefore, only cells that are obsolete should be removed with this routine        
        ! Input
        !-----
        ! grid%cell: struct with unstructured grid data as in Make_traduitoutb2us.m
        ! cells_to_remove: list of cells to remove ordered from small to larg cell number 

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)         :: grid
        integer(I8), intent(in)                 :: cells(:)

        ! Auxiliary
        integer(I8) :: i, j, s, n, nv, ic, nc, counter, c_num, &
            rem_ind_dummy(grid%cell%ntot)
        integer(I8), allocatable :: cellsU(:), rem_ind(:), range(:)


        ! Make unique
        call Unique(cells,cellsU) 

        associate(&
            c => grid%cell &
            )

        nc = size(cellsU)
        if (nc /= 0) then

            call c%x%Remove(cellsU)
            call c%y%Remove(cellsU)
            call c%psi%Remove(cellsU)
            call c%cflags%Remove(cellsU)
            call c%reg%Remove(cellsU)

            ! Adjust the c%vert and c%face array
            rem_ind_dummy = 0
            counter = 0
            do i = 1, nc
                ic = cellsU(i)
                s = c%faceP1%Get(ic)
                n = c%faceP2%Get(ic)
                rem_ind_dummy(counter+1:counter+n) = (/ (j, j=s, s+n-1)/)
                counter = counter + n
            end do

            rem_ind = rem_ind_dummy(1:counter)

            call c%vert%Remove(rem_ind)
            call c%face%Remove(rem_ind)

            ! Adjust the pointers
            do i = 1, nc

                ! Determine correct cell number
                c_num = cellsU(i) - (i-1)

                ! Get nv
                nv = c%vertP2%Get(c_num)

                ! Remove and adjust
                call c%vertP1%Remove(c_num)
                call c%vertP2%Remove(c_num)
                range = (/ (j, j = c_num, c%vertP1%Size())/)
                call c%vertP1%SumMask(range, -nv)

                ! Copy for c%face
                call c%faceP1%SetAllElementsArray(c%vertP1%Get())
                call c%faceP2%SetAllElementsArray(c%vertP2%Get())

            end do

            ! Change number of cells
            c%ntot = c%ntot - nc

        end if

        end associate

    end subroutine

    subroutine RemoveFaces(grid, faces)

        ! Description
        !============
        ! Removes faces out of face%vert and adapts the face number in cell%face

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)         :: grid
        integer(I8), intent(in)                 :: faces(:)

        ! Auxiliary
        integer(I8) :: i, nf, ind, ifs, face_num
        integer(I8), allocatable :: facesU(:)

        ! Make unique
        associate(&
            c => grid%cell, &
            f => grid%face, &
            fd => grid%data%fluxdata &
            )
        
        call Unique(faces, facesU)

        nf = size(facesU)
        if (nf /= 0) then

            call f%vert1%Remove(facesU)
            call f%vert2%Remove(facesU)
            call f%label%Remove(facesU)
            call f%aligned%Remove(facesU)

            do i = 1, nf

                ! Adjusting all face number in cell%face
                face_num = facesU(i) - (i-1)

                call c%face%UpdateArray(face_num)

                ! Remove face from fluxsurface if in a flux tube
                ind = findloc(fd%fluxsurfacefaces%Get(), face_num, 1)

                if (ind /= 0) then

                    ifs = GetFsFcFromFaceIndex(fd, ind)
                    call fd%fluxsurfacefaces%Remove(ind)
                    call fd%fluxsurfacefacesP2%SumMask(ifs, -1)
                    call fd%fluxsurfacefacesP1%UpdateArray(fd%fluxsurfacefacesP1%Get(ifs))

                end if

                ! Change other facenumber
                call fd%fluxsurfacefaces%UpdateArray(face_num)


            end do

            ! Adjust 
            f%ntot = f%ntot - nf

        end if

        end associate

    end subroutine

    subroutine RemoveVertices(grid, verts)

        ! Description
        !============
        !  Removes vertices out of coordinate and magnetic field data
        ! Adapts also the vertexnumbering in cell.vert and face.vert
        ! change coordinate data

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout) :: grid
        integer(I8), intent(in)         :: verts(:)

        ! Auxiliary
        integer(I8) :: i, j, nv, vx_num, ind, ifs
        integer(I8), allocatable :: vertsU(:)

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert, &
            fd => grid%data%fluxdata &
        )

        call Unique(verts, vertsU)

        nv = size(vertsU)
        if (nv /= 0) then

            call v%x%Remove(vertsU)
            call v%y%Remove(vertsU)
            call v%psi%Remove(vertsU)
            call v%bx%Remove(vertsU)
            call v%by%Remove(vertsU)

        end if

        do i = 1, nv

            ! Determine correct vertex number
            vx_num = vertsU(i) - (i-1)

            ! Cell%vert
            call c%vert%UpdateArray(vx_num)       

            ! Change face%vert
            call f%vert1%UpdateArray(vx_num)
            call f%vert2%UpdateArray(vx_num)

            ! Remove vertex from fsVx
            ind = findloc(fd%fluxsurfaceverts%Get(), vx_num, 1)

            if (ind /= 0) then

                ifs = GetFsVxFromVertInd(fd, ind)
                call fd%fluxsurfaceverts%Remove(ind)
                call fd%fluxsurfacevertsP2%SumMask(ifs, -1)
                call fd%fluxsurfacevertsP1%UpdateArray(fd%fluxsurfacevertsP1%Get(ifs))

            end if

            ! Update fsVx
            call fd%fluxsurfaceverts%UpdateArray(vx_num)

            ! Update xpointID
            do j = 1, grid%data%nxp
                if (grid%data%xpointID(j) .gt. vx_num) &
                    grid%data%xpointID(j) = grid%data%xpointID(j) - 1
            end do
           
        end do

        ! Change v%ntot
        v%ntot = v%ntot - nv

        ! Remove duplicated Face - TODO - probably not necessary

        end associate

    end subroutine

    subroutine GetFaceNumber(grid, v1, v2, meth, face_num)

        ! Description
        !============
        ! Give face number between two vertices. There are three methods:
        ! 1) Just check whether there exists a face and provide the number (flag1 = 0)
        ! 2) Check whether the face exist, if not, create a new face (flag1=1, flag2=0)
        ! 3) Immediatelly create a new face (flag1=1, flag2=1)     
        
        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)     :: grid
        integer(I8), intent(in)             :: v1, v2, meth
        integer(I8), intent(out)            :: face_num

        ! Auxiliary
        integer(I8) :: i, j
        integer(I8), allocatable, dimension(:) :: fc1, fc2

        ! Checks
        if (v1 == v2) call gdErrorHandler('GetFaceNumber: Same vertices')

        ! Initialize
        face_num = 0
        
        if (meth == 1 .or. meth == 2) then

            ! Check wheter face exist
            fc1 = GetVertFaceGA(grid%face, v1)
            fc2 = GetVertFaceGA(grid%face, v2)

            ! Search common face
            do i = 1, size(fc1)
                do j = 1, size(fc2)
                    if (fc1(i) == fc2(j)) then
                        face_num = fc1(i)
                        exit
                    end if
                end do
                if (face_num /= 0) return
            end do

            ! Create new face if necessary
            if (meth == 2 .and. face_num == 0) then

                ! Set default face properties
                face_num = grid%face%ntot + 1

                call grid%face%vert1%Append(v1)
                call grid%face%vert2%Append(v2)
                call grid%face%label%Append(0)
                call grid%face%aligned%Append(0)
                grid%face%ntot = face_num

            else

                call gdErrorHandler('GetFaceNumber: no face found or created')

            end if 

        else if (meth == 3) then

            ! Set default face properties
            face_num = grid%face%ntot + 1

            call grid%face%vert1%Append(v1)
            call grid%face%vert2%Append(v2)
            call grid%face%label%Append(0)
            call grid%face%aligned%Append(0)
            grid%face%ntot = face_num

        else

            call gdErrorHandler('GetFaceNumber: meth not implemented')

        end if

    end subroutine

    subroutine AddFaceToFsFc(grid, new_f, old_fs)

        ! Description
        !============
        ! Adds a new faces the fsFc to the correct flux surface based on an old face
        ! of the same flux surface. Old_fs are face on the same flux surface
        ! as the new face

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout) :: grid
        integer(I8), intent(in)         :: new_f, old_fs(:)    

        ! Auxiliary
        integer(I8) :: i, old_nf, s, ifs, nf
        integer(I8), allocatable :: ind(:), ifss(:), ifssU(:), fcs(:), &
            new_faces(:), new_facesU(:), loc(:)

        associate(&
            fd => grid%data%fluxdata &
            )

        ! Find the flux surface of old faces
        nf = size(old_fs)
        allocate(ind(nf))
        ind = 0
        do i = 1, nf
            ind(i) = findloc(fd%fluxsurfacefaces%Get(), old_fs(i),1)
        end do

        if (sum(ind) /= 0) then

            ! Get the correct flux surface
            allocate(ifss(nf))
            ifss = 0
            do i = 1, nf
                if (ind(i) /= 0) ifss(i) = GetFsFcFromFaceIndex(fd,ind(i))
            end do

            ! Unique 
            call Unique(ifss, ifssU)

        end if

        if (allocated(ifssU)) then

            ! Add the face to that flux surface
            ifs = ifssU(1)
            old_nf = fd%fluxsurfacefacesP2%Get(ifs)
            s = fd%fluxsurfacefacesP1%Get(ifs)
            fcs = GetFluxSurfaceFcsGA(fd, ifs)

            ! Build new faces
            allocate(new_faces(old_nf + 1))
            new_faces =  [fcs , new_f]
            call Unique(new_faces, new_facesU)

            ! Put in the data 
            call fd%fluxsurfacefaces%Replace(s,s+old_nf-1,new_faces)
            call fd%fluxsurfacefacesP2%Set(ifs, size(new_faces))
            loc = (/ (i, i = ifs+1, fd%nFs)/)
            call fd%fluxsurfacefacesP1%SumMask(loc,1)

        end if
            
        end associate

    end subroutine

    subroutine AddCell(grid, faces, verts, reg, ic)

        ! Description
        !============
        ! Add cell in connectivity cvFc and cvVx at the back.

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout) :: grid
        integer(I8), intent(in)         :: reg, faces(:), verts(:)
        integer(I8), intent(out)        :: ic

        ! Auxiliary
        integer(I8) :: ind, nv

        associate(&
            c => grid%cell &
            )

        ! Check
        if (size(faces) /= size(verts)) then
            call gdErrorHandler('AddCell: size faces and size verts is not the same')
        end if

        ! Adjust and add the cell
        c%ntot = c%ntot + 1
        ic = c%ntot

        ind = c%vertP1%Get(ic-1) + c%vertP2%Get(ic-1)
        nv = size(verts)

        ! Verts
        call c%vertP1%Append(ind)
        call c%vertP2%Append(nv)
        call c%vert%Append(verts)

        ! Face
        call c%faceP1%Append(ind)
        call c%faceP2%Append(nv)
        call c%face%Append(faces) 
        
        ! Cell
        call c%cflags%Append(0)
        call c%reg%Append(reg)        
        call c%psi%Append(0.0_R8)


        ! Calculate centroid
        call c%x%Append(sum(grid%vert%x%Get(verts))/real(nv, kind=R8))
        call c%y%Append(sum(grid%vert%y%Get(verts))/real(nv, kind=R8))




        end associate

    end subroutine

    subroutine AddVert(grid, magneticField, face, vert, type, psic, common_vert)

        ! Description
        !============
        ! Add a new vertex on the segment defined by the face

        ! Declare variables
        !==================
        class(GAGridUDT), intent(inout)     :: grid
        type(MagneticFieldUDT), intent(in)  :: magneticField
        integer(I8), intent(in)             :: face
        integer(I8), intent(out)            :: vert
        character(*), optional              :: type  
        real(R8), optional                  :: psic
        integer(I8), optional               :: common_vert

        ! Auxiliary
        integer(I8) :: v1, v2
        real(R8) :: t0, v1x, v1y, v2x, v2y, fcX, fcY, d, &
            t1, vec_v1(2), vec_v2(2), X1(2), Y1(2), X2(2), Y2(2), &
            x_int, y_int, vcx, vcy, d_v1(2), d_v2(2), d2, dn, &
            fcX_int, fcY_int, d1, vec_c(2)
        real(R8), allocatable, dimension(:) :: v1_nx, v1_ny, v1_psi, v1_bx, v1_by


        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        allocate(v1_nx(1), v1_ny(1), v1_psi(1), v1_bx(1), v1_by(1))

        if (.not.present(type)) then

            ! Make new vertices at the crossing with psic
            v1 = f%vert1%Get(face)
            v2 = f%vert2%Get(face)
            v1_nx = 0.5_R8 * (v%x%Get(v1) + v%x%Get(v2))
            v1_ny = 0.5_R8 * (v%y%Get(v1) + v%y%Get(v2))
            call magneticField%interp%Evaluate(v1_nx, v1_ny, 0, 0, v1_psi)
            call magneticField%interp%Evaluate(v1_nx, v1_ny, 1, 0, v1_bx)
            call magneticField%interp%Evaluate(v1_nx, v1_ny, 0, 1, v1_by)

        else 

            select case (type)

            case ('psi')

                ! Make new vertices at the crossing with psic
                v1 = f%vert1%Get(face)
                v2 = f%vert2%Get(face)
                t0 = (psic - v%psi%Get(v1)) / (v%psi%Get(v2) - v%psi%Get(v1))

                ! Get vertex properties
                v1_nx = v%x%Get(v1) + t0 * (v%x%Get(v2) - v%x%Get(v1))
                v1_ny = v%y%Get(v1) + t0 * (v%y%Get(v2) - v%y%Get(v1))
                v1_psi = psic
                call magneticField%interp%Evaluate(v1_nx, v1_ny, 1, 0, v1_bx)
                call magneticField%interp%Evaluate(v1_nx, v1_ny, 0, 1, v1_by)
            
            case ('pol')

                ! Make new vertices on the average of intersection of constant psi vectors 
                ! and the face center if face is aligned - TODO upgrade to construct Cubic Hermite polynomal for psi-curve
                v1 = f%vert1%Get(face)
                v2 = f%vert2%Get(face)
                if (f%aligned%Get(face) == 1) then

                    v1x = v%x%Get(v1)
                    v1y = v%y%Get(v1)
                    v2x = v%x%Get(v2)
                    v2y = v%y%Get(v2)
                    fcX = 0.5_R8 * (v1x + v2x)
                    fcY = 0.5_R8 * (v1y + v2y)
                    d = sqrt((v1x-v2x)**2 + (v1y-v2y)**2)

                    ! Normalized vector out of vert1
                    d1 = sqrt(v%bx%Get(v1)**2 + v%by%Get(v1)**2)
                    vec_v1(1) = -v%by%Get(v1)/d1
                    vec_v1(2) = v%bx%Get(v1)/d1

                    ! Face normal
                    t0 = v2y - v1y
                    t1 = -(v2x - v1x)
                    dn = sqrt(t0**2 + t1**2)
                    vec_v2(1) = t0/dn
                    vec_v2(2) = t1/dn

                    ! Normalized vector out of vert2
                    X1 = [v1x + vec_v1(1)*2*d, v1x - vec_v1(1)*2*d]
                    Y1 = [v1y + vec_v1(2)*2*d, v1y - vec_v1(2)*2*d]
                    X2 = [fcX + vec_v2(1)*2*d, fcX - vec_v2(1)*2*d]
                    Y2 = [fcY + vec_v2(2)*2*d, fcY - vec_v2(2)*2*d]

                    call SegmentIntersections(x_int, y_int, X1(1), Y1(1), X1(2), Y1(2), X2(1), Y2(1), X2(2), Y2(2))

                    if (.not.isnan(x_int)) then

                        ! Approximate constant psi curve - average of intersection
                        ! point and old face center
                        v1_nx = 0.5_R8 *(x_int + 0.5_R8 * (v1x + v2x))
                        v1_ny = 0.5_R8 *(y_int + 0.5_R8 * (v1y + v2y))
                        call magneticField%interp%Evaluate(v1_nx, v1_ny, 0, 0, v1_psi)
                        call magneticField%interp%Evaluate(v1_nx, v1_ny, 0, 0, v1_psi)
                        call magneticField%interp%Evaluate(v1_nx, v1_ny, 0, 0, v1_psi)

                    else 

                        call gdErrorHandler('AddVert: no intersection found')

                    end if

                else 

                    ! Geometric splitting
                    v1_nx = 0.5_R8 * (v%x%Get(v1) + v%x%Get(v2))
                    v1_ny = 0.5_R8 * (v%y%Get(v1) + v%y%Get(v2))
                    call magneticField%interp%Evaluate(v1_nx, v1_ny, 0, 0, v1_psi)
                    call magneticField%interp%Evaluate(v1_nx, v1_ny, 1, 0, v1_bx)
                    call magneticField%interp%Evaluate(v1_nx, v1_ny, 0, 1, v1_by)                   

                end if

            case ('pol_ortho')

                ! Make new vertices on intersection of constant psi vectors
                ! if face is aligned
                v1 = f%vert1%Get(face)
                v2 = f%vert2%Get(face)
                if (f%aligned%Get(face) == 1) then

                    v1x = v%x%Get(v1)
                    v1y = v%y%Get(v1)
                    v2x = v%x%Get(v2)
                    v2y = v%y%Get(v2)
                    vcx = v%x%Get(common_vert)
                    vcy = v%y%Get(common_vert)
                    fcX = 0.5_R8 * (v1x + v2x)
                    fcY = 0.5_R8 * (v1y + v2y)
                    d = sqrt((v1x-v2x)**2 + (v1y-v2y)**2)

                    d_v1(1) = v1x - fcX
                    d_v1(2) = v1y - fcY

                    d_v2(1) = v2x - fcX
                    d_v2(2) = v2y - fcY

                    ! Normalized vector out of vert1
                    d1 = sqrt(v%bx%Get(v1)**2 + v%by%Get(v1)**2) 
                    vec_v1(1) = - v%by%Get(v1) / d1
                    vec_v1(2) = v%bx%Get(v1) / d1

                    ! Face normal
                    t0 = v2y - v1y
                    t1 = -(v2x - v1x)
                    dn = sqrt(t0**2 + t1**2)
                    vec_v2(1) = t0/dn
                    vec_v2(2) = t1/dn

                    ! ORTHOGONAL INTERSECTION
                    ! FcX_int and fcY_int is the intersection point of orthogonal line
                    ! out of the hanging node of the pentagon with opposite face in that pentagon
                    ! The orthogonal line starts in the common vert
                    ! Normalized vector out of common vert
                    d2 = sqrt(v%bx%Get(common_vert)**2 + v%by%Get(common_vert)**2)
                    vec_c(1) = v%bx%Get(common_vert) / d2
                    vec_c(2) = v%by%Get(common_vert) / d2  
                    
                    ! Define segment 75% of original face centers around the face center
                    X1 = [fcX + 0.75_R8 * d_v2(1), fcX + 0.75_R8 * d_v1(1)] 
                    Y1 = [fcX + 0.75_R8 * d_v2(2), fcX + 0.75_R8 * d_v1(2)]

                    ! Orthogonal vector out of the common vertex (hanging node of pentagon)
                    X2 = [vcx + vec_c(1)*10*d, vcx - vec_c(1)*10*d]
                    Y2 = [vcy + vec_c(2)*10*d, vcy - vec_c(2)*10*d]

                    call SegmentIntersections(fcX_int, fcY_int, X1(1), Y1(1), X1(2), Y1(2), X2(1), Y2(1), X2(2), Y2(2))

                    if (isnan(fcX_int)) call gdErrorHandler('AddVert: No intersection found')

                    ! Average face center with new face intersection
                    fcX = 0.3*fcX + 0.7*fcX_int
                    fcY = 0.3*fcY + 0.7*fcY_int

                    ! INTERSECTION FOR NEW POINT
                    ! intersection between face normal out if new (fcX,fcY)
                    ! and the magnetic vector out of v1

                    X1 = [v1x + vec_v1(1)*d, v1x - vec_v1(1)*d]
                    Y1 = [v1y + vec_v1(2)*d, v1y - vec_v1(2)*d]
                    X2 = [fcX + vec_v2(1)*d, fcX - vec_v2(1)*d]
                    Y2 = [fcY + vec_v2(2)*d, fcY - vec_v2(2)*d]

                    call SegmentIntersections(x_int, y_int, X1(1), Y1(1), X1(2), Y1(2), X2(1), Y2(1), X2(2), Y2(2))

                    if (.not.isnan(x_int)) then

                        ! Approximate constant psi curve - average of intersection and old face center
                        v1_nx = 0.5_R8 * (x_int + fcX)
                        v1_ny = 0.5_R8 * (y_int + fcY)
                        call magneticField%interp%Evaluate(v1_nx, v1_ny, 0, 0, v1_psi)
                        call magneticField%interp%Evaluate(v1_nx, v1_ny, 1, 0, v1_bx)
                        call magneticField%interp%Evaluate(v1_nx, v1_ny, 0, 1, v1_by)

                    else

                        call gdErrorHandler('AddVert: no intersection found')

                    end if

                else 

                    ! Geometric splitting
                    v1_nx = 0.5_R8 * (v%x%Get(v1) + v%x%Get(v2))
                    v1_ny = 0.5_R8 * (v%y%Get(v1) + v%y%Get(v2))
                    call magneticField%interp%Evaluate(v1_nx, v1_ny, 0, 0, v1_psi)
                    call magneticField%interp%Evaluate(v1_nx, v1_ny, 1, 0, v1_bx)
                    call magneticField%interp%Evaluate(v1_nx, v1_ny, 1, 0, v1_by) 

                end if
            case ('geometric')

                ! Geometric splitting
                v1 = f%vert1%Get(face)
                v2 = f%vert2%Get(face)
                v1_nx = 0.5_R8 * (v%x%Get(v1) + v%x%Get(v2))
                v1_ny = 0.5_R8 * (v%y%Get(v1) + v%y%Get(v2))
                call magneticField%interp%Evaluate(v1_nx, v1_ny, 0, 0, v1_psi)
                call magneticField%interp%Evaluate(v1_nx, v1_ny, 1, 0, v1_bx)
                call magneticField%interp%Evaluate(v1_nx, v1_ny, 1, 0, v1_by)  

            case default

                call gdErrorHandler('AddVert: other type of face splitting not implemented.')

            end select

        end if

        ! Append the vertex
        call v%x%Append(v1_nx)
        call v%y%Append(v1_ny)
        call v%psi%Append(v1_psi)
        call v%bx%Append(v1_bx)
        call v%by%Append(v1_by)

        ! Increase total number of vertices and output the vertex number
        v%ntot = v%ntot + 1
        vert = v%ntot

        end associate

    end subroutine 

    subroutine AddVertToFsVx(grid, new_v, new_f, old_f, type)

        ! Description
        !============
        ! Adds the new vertices to the correct flux surface if possible based
        ! on information from the existing faces.

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)         :: grid
        integer(I8), intent(in)                 :: old_f, new_v(:), new_f(:)
        character(:), allocatable, intent(in)   :: type

        ! Auxiliary
        integer(I8) :: i, ind, ifs, ind_good

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert, &
            fd => grid%data%fluxdata &
            )

        select case (type)

        case ('pol')

            ! Find the ifs of the face, only for aligned faces
            if (f%aligned%Get(old_f) == 1) then

                if (count(fd%fluxsurfacefaces%Get() == old_f) /= 1) then
                    call gdErrorHandler('AddVertToFsVx: multiple possible flux surface to add vertex to')
                else
                    ind = findloc(fd%fluxsurfacefaces%Get(), old_f, 1)
                end if

                ! If old face not in a surface, add new surface
                if (ind == 0) then

                    call grid%AddFluxSurface(new_v, new_f)

                !  Or at the same surface of the old face
                else 

                    ! Get the flux surface id
                    ifs = GetFsFcFromFaceIndex(fd, ind)
                    call grid%AddIntoExistingSurface(ifs, new_v, new_f)
                             
                end if
                
            end if

        case ('rad')

            ! Find ifs of the given vertex
            ind_good = 0
            do i = 1, size(new_v)
                ind = findloc(fd%fluxsurfaceverts%Get(), new_v(i), 1)
                if (ind /= 0) then
                    ind_good = ind
                    exit
                end if
            end do

            if (ind_good == 0) then

                call grid%AddFluxSurface(new_v, new_f)

            else

                ifs = GetFsVxFromVertInd(fd, ind_good)
                call grid%AddIntoExistingSurface(ifs, new_v, new_f)

            end if


        end select

        end associate

    end subroutine

    subroutine AddFluxSurface(grid, vxs, fcs)

        ! Description
        !============
        ! Adding a new flux surface in fsFc and fsVx

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)         :: grid
        integer(I8), intent(in)                 :: vxs(:), fcs(:)

        ! Auxiliary

        associate(&
            fd => grid%data%fluxdata &
            )

        ! Increase number of flux surfaces
        fd%nFs = fd%nFs + 1

        ! Verts
        call fd%fluxsurfacevertsP1%Append(fd%fluxsurfacevertsP1%Get(fd%nFs-1) + fd%fluxsurfacevertsP2%Get(fd%nFs-1))
        call fd%fluxsurfacevertsP2%Append(size(vxs))
        call fd%fluxsurfaceverts%Append(vxs)

        ! Faces
        call fd%fluxsurfacefacesP1%Append(fd%fluxsurfacefacesP1%Get(fd%nFs-1) + fd%fluxsurfacefacesP2%Get(fd%nFs-1))
        call fd%fluxsurfacefacesP2%Append(size(fcs))
        call fd%fluxsurfaceverts%Append(fcs)

        end associate

    end subroutine

    subroutine AddIntoExistingSurface(grid, ifs, vxs, fcs)

        ! Description
        !============
        ! Add vxs and fcs into an existing flux surface

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout)         :: grid
        integer(I8), intent(in)                 :: ifs, vxs(:), fcs(:)

        ! Auxiliary
        integer(I8) :: s, j
        integer(I8), allocatable :: vxs_current(:), fcs_current(:), &
            new_verts(:), new_vertsD(:), new_faces(:), new_facesD(:), range(:)

        associate(&
            fd => grid%data%fluxdata &
            )

        ! Add the vertex to that flux surface
        vxs_current = GetFluxSurfaceVxsGA(fd, ifs)
        new_vertsD = [vxs_current, vxs]
        call Unique(new_vertsD, new_verts)
        s = fd%fluxsurfacevertsP1%Get(ifs)
        call fd%fluxsurfaceverts%Replace(s,s+fd%fluxsurfacevertsP2%Get(ifs)-1, new_verts)
        call fd%fluxsurfacevertsP2%Set(ifs, size(new_verts))
        range = (/ (j, j = ifs+1, fd%nFs)/)
        call fd%fluxsurfacefacesP1%SumMask(range, size(vxs))

        ! Add the faces to that flux surfaces
        fcs_current = GetFluxSurfaceFcsGA(fd, ifs)
        new_facesD = [fcs_current, fcs]
        call Unique(new_facesD, new_faces)
        s = fd%fluxsurfacefacesP1%Get(ifs)
        call fd%fluxsurfacefaces%Replace(s,s+fd%fluxsurfacefacesP2%Get(ifs)-1, new_faces)
        call fd%fluxsurfacefacesP2%Set(ifs, size(new_faces))
        call fd%fluxsurfacefacesP1%SumMask(range, size(fcs)) 
        
        end associate

    end subroutine

    !------------------------------------------------------------------!
    !                         VISUALIZATION                            !
    !------------------------------------------------------------------!    

    subroutine WriteGAGridData(grid, filename)

        ! Description
        !============
        ! Writing out the GAgrid to plot it.

        ! 'vertices'
        ! <vert%ntot>
        ! 'ID, x, y'
        ! <ID, x, y>
        ! 'faces'
        ! <face%ntot> 
        ! 'ID, v1, v2, label'
        ! <ID, v1, v2, label>
        ! 'cells'
        ! <cell%ntot, cell%nvert> 
        ! 'ID, vp1, vp2, region>'
        ! <ID, vp1, vp2, region>
        ! 'cell vertices'
        ! <cell%vert> 

        ! Declare variables
        !==================
        ! Modules 
        use mod_plotter 
        use mod_specialchars, only : filesepchar


        ! Arguments
        class(GAGridUDT)                        :: grid
        character(*), intent(in)                :: filename 

        ! Auxiliary
        integer                                 :: fu
        real(R8), allocatable, dimension(:)     :: x, y, cx, cy
        integer(I8), allocatable, dimension(:)  :: v1, v2, region, &   
            label, vc, aligned
        character(:), allocatable               :: dir

        ! Loop
        integer(I8)                             :: i


        ! Initialize
        !===========
        ! Unpack
        associate(&
            f       => grid%face,   &
            c       => grid%cell,   &
            v       => grid%vert    &
        )

        ! Construct writing directory
        dir = plotdir // filesepchar // filename // '.dat'

        ! Open file
        open (action='write', file=trim(dir), newunit=fu, &
             status='unknown')

        ! Write header
        write(fu, *) 'VERSION3.00.00'

        ! Write vertex data
        !==================
        ! Unpack
        x = v%x%Get()
        y = v%y%Get()

        ! Number of vertices
        write (fu, *) 'vertices'
        write (fu, *) v%ntot 

        ! Vertex data
        write (fu, *) 'ID, x, y'
        do i = 1, v%ntot 
            write (fu, *) i, x(i), y(i)
        end do 

        ! Write face data
        !================
        ! Unpack
        v1 = f%vert1%Get()
        v2 = f%vert2%Get()
        label = f%label%Get()
        aligned = f%aligned%Get()

        ! Number of faces
        write (fu, *) 'faces'
        write (fu, *) f%ntot

        ! Face data
        write (fu, *) 'ID, v1, v2, label, aligned'
        do i = 1, f%ntot
            write (fu, *) i, v1(i), v2(i), label(i), aligned(i)
        end do 

        ! Write cell data
        !================
        ! Unpack
        vc = c%vert%Get()
        v1 = c%vertP1%Get()
        v2 = c%vertP2%Get()
        region = c%reg%Get()
        cx = c%x%Get()
        cy = c%y%Get()

        ! Number of cells
        write (fu, *) 'cells'
        write (fu, *) c%ntot, size(vc)

        ! Cell data
        write (fu, *) 'ID, vp1, vp2, region, x, y'
        do i = 1, c%ntot
            write (fu, *) i, v1(i), v2(i), region(i), cx(i), cy(i)
        end do 

        ! Cell vertices
        write (fu, *) 'cell vertices'
        do i = 1, size(vc)
            write (fu, *) vc(i)
        end do 

        ! Housekeeping
        !=============
        ! Deallocate again

        ! Others
        end associate
        close(fu)

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
        integer(I8), intent(in)                 :: f_list(:)
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

                    ! Housekeeping
                    deallocate(a)

                end do

            end do

        end do

    end subroutine

    subroutine CheckFace(f, v1, v2, face_num)

        ! Description
        !============
        ! Gives face number of face between two verts if any

        ! Declare variables
        !==================
        ! Arguments
        class(GAFaceUDT), intent(in)    :: f
        integer(I8), intent(in)         :: v1, v2
        integer(I8), intent(out)        :: face_num

        ! Auxiliary
        integer(I8) :: i, j
        integer(I8), allocatable :: fc1(:), fc2(:)


        ! Checks
        if (v1 == v2) call gdErrorHandler('CheckFace: Same vertices')

        ! Initialize
        face_num = 0

        ! Check wheter face exist
        fc1 = GetVertFaceGA(f, v1)
        fc2 = GetVertFaceGA(f, v2)

        ! Search common face
        do i = 1, size(fc1)
            do j = 1, size(fc2)
                if (fc1(i) == fc2(j)) then
                    face_num = fc1(i)
                    exit
                end if
            end do
            if (face_num /= 0) exit
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
            allocate(fcs_t1(count(fcLbl_loc.ge.4)))
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

    ! Get neighboring cells of a cell
    function GetCellNeigsGA(g, ic, cvLookUp) result(res)
        type(GAGridUDT) :: g
        integer(I8) :: ic, i, j, s, ifc, counter, nf
        integer(I8), optional :: cvLookUp(:)
        integer(I8), allocatable :: cvs(:), res(:), res_dummy(:)

        if (.not.present(cvLookUp)) then

            s = g%cell%faceP1%Get(ic)
            nf = g%cell%faceP2%Get(ic)
            counter = 0
            allocate(res_dummy(nf))
            do i = 1, nf
                ifc = g%cell%face%Get(s+i-1)
                cvs = GetFaceCellGA(g%cell,ifc)

                ! Add cell if not equal to ic
                do j = 1,size(cvs)
                    if (cvs(j) /= ic) then
                        counter = counter + 1
                        res_dummy(counter) = cvs(j) 
                    end if
                end do

            end do

            res = res_dummy(1:counter)

        else 

            s = g%cell%faceP1%Get(ic)
            nf = g%cell%faceP2%Get(ic)
            counter = 0
            allocate(res_dummy(nf))
            do i = 1, nf
                ifc = g%cell%face%Get(s+i-1)
                cvs = GetFaceCellGA(g%cell,ifc,cvLookUp)

                ! Add cell if not equal to ic
                do j = 1,size(cvs)
                    if (cvs(j) /= ic) then
                        counter = counter + 1
                        res_dummy(counter) = cvs(j) 
                    end if
                end do

            end do

            res = res_dummy(1:counter)

        end if

        ! House keeping
        deallocate(res_dummy)

        s = g%cell%faceP1%Get(ic)
        nf = g%cell%faceP2%Get(ic)
        counter = 0
        allocate(res_dummy(nf))
        do i = 1, nf
            ifc = g%cell%face%Get(s+i-1)
            cvs = GetFaceCellGA(g%cell,ifc,cvLookUp)

            ! Add cell if not equal to ic
            do j = 1,size(cvs)
                if (cvs(j) /= ic) then
                    counter = counter + 1
                    res_dummy(counter) = cvs(j) 
                end if
            end do

        end do

        res = res_dummy(1:counter)

    end function

    function GetCellFromFaceIndex(c, indf) result(res)
        type(GACellUDT) :: c
        integer(I8) :: i, indf, n_el, res, faceP1(c%ntot)
        integer(I8), allocatable :: b(:), ind(:)

        faceP1 = c%faceP1%Get()
        n_el = count(faceP1.le.indf)
        allocate(b(n_el))
        ind = (/ (i, i = 1, n_el)/)
        b = pack(ind, faceP1.le.indf)
        res = b(n_el)
    end function

    function GetCellFromVertIndex(c, indv) result(res)
        type(GACellUDT) :: c
        integer(I8) :: i, indv, n_el, res, vertP1(c%ntot)
        integer(I8), allocatable :: b(:), ind(:)

        vertP1 = c%vertP1%Get()
        n_el = count(vertP1.le.indv)
        allocate(b(n_el))
        ind = (/ (i, i = 1, n_el)/)
        b = pack(ind, vertP1.le.indv)
        res = b(n_el)        
    end function

    ! Get cells of a face with dynamic arrays
    function GetFaceCellGA(cell, i, cvLookUp) result(res)
        integer(I8)                 :: i, j
        type(GACellUDT)             :: cell
        integer(I8), optional       :: cvLookUp(:)
        integer(I8), allocatable    :: res(:), indcf(:), ind(:)
        logical, allocatable        :: log(:) 
        
        if (present(cvLookUp)) then
            allocate(res(count(cell%face%GetAllElements().eq.i)))
            res = pack(cvLookUp,cell%face%GetAllElements().eq.i)
        else 
            log = (cell%face%Get() == i)
            allocate(ind(count(log)))
            indcf = (/ (j, j = 1, cell%face%Size()) /)
            ind = pack(indcf, log)

            allocate(res(size(ind)))
            do j = 1, size(ind)
                res(j) = GetCellFromFaceIndex(cell, ind(j))
            end do

        end if

    end function

    ! Get faces of a vertex with dynamic arrays
    function GetVertFaceGA(face, i) result(res)
        type(GAFaceUDT) :: face
        integer(I8) :: i, j, n1, n2, &
            fv1(face%ntot), fv2(face%ntot), indf(face%ntot)
        integer(I8), allocatable :: res(:), res1(:), &
            res2(:), res3(:)

        fv1 = face%vert1%Get()
        fv2 = face%vert2%Get()

        allocate(res1(count(fv1 == i)))
        allocate(res2(count(fv2 == i)))
        indf = (/ (j, j = 1, face%ntot)/)
        res1 = pack(indf,fv1 == i)
        res2 = pack(indf,fv2 == i)
        n1 = size(res1)
        n2 = size(res2)
        allocate(res3(n1+n2))
        res3(1:n1) = res1 
        res3(n1+1:n1+n2) = res2 

        call Unique(res3, res)

    end function

    ! Get the fluxsurface of a vertex
    function GetVertFsvGA(fd, i, fsvLookUp) result(res)
        type(GAFluxDataUDT) :: fd
        integer(I8) :: i, res
        integer(I8), allocatable :: res1(:)
        integer(I8), allocatable, optional :: fsvLookUp(:)

        res = 0

        if (.not.present(fsvLookUp)) fsvLookUp = GetFsvLookUpGA(fd)
        allocate(res1(count(fd%fluxsurfaceverts%Get() == i)))
        res1 = pack(fsvLookUp,fd%fluxsurfaceverts%Get() == i)
        if (count(fd%fluxsurfaceverts%Get() == i) .lt. 1) return
        res = res1(1)

    end function
 
    function GetFsFcFromFaceIndex(fd, indf) result(res)
        type(GAFluxDataUDT) :: fd
        integer(I8) :: i, n_el, indf, res, fdP1(fd%nFs)
        integer(I8), allocatable :: b(:), ind(:)

        fdP1 = fd%fluxsurfacefacesP1%Get()
        n_el = count(fdP1.le.indf)
        allocate(b(n_el))
        ind = (/ (i, i = 1, n_el)/)
        b = pack(ind, fdP1.le.indf)
        res = b(n_el)

    end function

    function GetFsVxFromVertInd(fd, indv) result(res)
        type(GAFluxDataUDT) :: fd 
        integer(I8) :: i, n_el, indv, res, fdP1(fd%nFs)
        integer(I8), allocatable :: b(:), ind(:)

        fdP1 = fd%fluxsurfacevertsP1%Get()
        n_el = count(fdP1.le.indv)
        allocate(b(n_el))
        ind = (/ (i, i = 1, n_el) /)
        b = pack(ind, fdP1.le.indv)
        res = b(n_el)       
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
        integer(I8)                 :: i, j, k, counter, cell_num
        type(GACellUDT)               :: cell
        integer(I8), optional    :: cvLookUp(:)
        integer(I8), allocatable    :: res(:), indc(:), ind(:), covD(:)
        logical, allocatable :: log(:) 
        
        if (.not.present(cvLookUp)) then
            ! Do without cvLookUp
            !cvLookUp = GetCvLookUpGA(cell)
            allocate(covD(cell%ntot))
            covD = 0
            counter = 0

            indc = (/( j, j = 1, cell%vertP1%Get(cell%ntot) + cell%vertP2%Get(cell%ntot)-1)/)
            log = (cell%vert%Get() == i)
            allocate(ind(count(log)))
            ind = pack(indc, log)

            do k = 1, size(ind)
                cell_num = GetCellFromVertIndex(cell, ind(k))
                counter = counter + 1
                covD(counter) = cell_num
            end do

            res = covD(1:counter)
        else
            log = (cell%vert%Get().eq.i)
            allocate(res(count(log)))
            res = pack(cvLookUp,log)
        end if

    end function

    function GetCvLookUpGA(cell) result(res)
        type(GACellUDT)       :: cell
        integer(I8)         :: nc, ic, nv, s               
        integer(I8), allocatable :: res(:)! , range(:)

        nc = cell%ntot
        !range = (/ (i, i = 1,(cell%vertP1%Get(nc)+cell%vertP2%Get(nc)-1))/)

        allocate(res(cell%vertP1%Get(nc)+cell%vertP2%Get(nc)-1))
        res = 0

        do ic = 1, nc
            s = cell%vertP1%Get(ic)
            nv = cell%vertP2%Get(ic)
            res(s:s+nv-1) = ic
            !range = (/ (i, i = s, (s+nv-1)) /)
            !res(range) = ic
        end do
    end function

    function GetFsvLookUpGA(fd) result(res)
        type(GAFluxDataUDT) :: fd
        integer(I8) :: ifs, nv, sv
        integer(I8), allocatable :: res(:)
        allocate(res(fd%fluxsurfacevertsP1%Get(fd%nFs) +  &
            fd%fluxsurfacevertsP2%Get(fd%nFs)-1)) 
        res = 0
        do ifs = 1, fd%nFs
            nv = fd%fluxsurfacevertsP2%Get(ifs)
            sv = fd%fluxsurfacevertsP1%Get(ifs)
            res(sv:sv+nv-1) = ifs
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
        integer(I8) :: fcs(:)
        integer(I8), allocatable :: verts(:), res(:)
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
        integer(I8) :: ifs, nf, counter, i, step, ifc, reg1, reg2, sepIDloc(1:100), cvLookUp(:)
        type(GACellUDT) ::  cell
        integer(I8), allocatable :: faces(:), cvs(:), res(:), creg(:)

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

    function isBoundaryFace0DGA(f, ifc) result(res)
        integer(I8) :: ifc
        type(GAFaceUDT) :: f
        logical :: res

        res = .true. 
        if (f%label%Get(ifc) == 0) then
            res = .false.
        end if
    end function

    function isBoundaryFace1DGA(f, tf) result(res)
        integer(I8) :: tf(:)
        type(GAFaceUDT) :: f
        logical, allocatable :: res(:)

        allocate(res(size(tf)))
        res = (f%label%Get(tf) /= 0)

        !do i = 1, nf
        !    if (f%label%Get(i) == 0) then
        !        res = .false.
        !    end if
        !end do

    end function

    function isBoundaryCell0DGA(grid, cell) result(res)
        integer(I8), allocatable :: fcs(:)
        type(GAGridUDT) :: grid
        logical, allocatable :: b_flag(:)
        integer(I8) :: cell
        logical :: res

        res = .false.
        b_flag = (grid%face%label%Get() /= 0)
        fcs = GetCellFaceGA(grid%cell, cell)
        if (any(b_flag(fcs))) res = .true.
    
    end function 

    function isBoundaryCell1DGA(grid, cell) result(res)
        integer(I8), allocatable :: fcs(:)
        type(GAGridUDT) :: grid
        logical, allocatable :: b_flag(:)
        integer(I8) :: cell(:)
        logical, allocatable :: res(:)
        integer(I8) :: i

        allocate(res(size(cell)))
        res = .false.
        b_flag = (grid%face%label%Get() /= 0)
        do i = 1, size(cell)
            fcs = GetCellFaceGA(grid%cell, cell(i))
            if (any(b_flag(fcs))) res(i) = .true.
        end do
    
    end function
    function isBoundaryVert0DGA(grid, iv) result(res)
        type(GAGridUDT) :: grid
        integer(I8) :: iv
        logical :: res
        integer(I8), allocatable :: fcs(:)

        res = .false.
        fcs = GetVertFaceGA(grid%face, iv)
        if (any(isBoundaryFaceGA(grid%face,fcs))) then
            res = .true.
        end if

    end function

    function isBoundaryVert1DGA(grid, verts) result(res)
        type(GAGridUDT) :: grid
        integer(I8) :: i, verts(:)
        logical, allocatable :: res(:)
        integer(I8), allocatable :: fcs(:) 

        allocate(res(size(verts)))
        res = .false.
        do i = 1, size(verts)
            fcs = GetVertFaceGA(grid%face, verts(i))
            if (any(isBoundaryFaceGA(grid%face, fcs))) then
                res(i) = .true.
            end if   
        end do
     
    end function

    function isPoloidal(grid, fcs, cvLookUp) result(res)
        type(GAGridUDT) :: grid
        integer(I8) :: fcs, cvLookUp(:)
        logical :: res
        integer(I8), allocatable :: fcs1(:), fcs2(:), cvs(:), &
            faces(:)

        res = .false.
        if (grid%face%aligned%Get(fcs) == 0) then
            cvs = GetFaceCellGA(grid%cell, fcs, cvLookUp)
            if (size(cvs) == 2) then
                fcs1 = GetCellFaceGA(grid%cell, cvs(1))
                fcs2 = GetCellFaceGA(grid%cell, cvs(2))

                if (size(fcs1) + size(fcs1) == 8) then

                    faces = [fcs1, fcs2]
                    if ((sum(grid%face%aligned%Get(faces)) == 4) &
                        .or. (sum(grid%face%aligned%Get(faces)) == 3)) &
                        res = .true.

                else if (size(fcs1) == 3 .and. size(fcs2) == 4) then

                    if ((sum(grid%face%aligned%Get(fcs1)) == 1) &
                        .and. (sum(grid%face%aligned%Get(fcs2)) == 2)) &
                        res = .true.

                else if (size(fcs1) == 4 .and. size(fcs2) == 3) then
                    if ((sum(grid%face%aligned%Get(fcs1)) == 2) &
                        .and. (sum(grid%face%aligned%Get(fcs2)) == 1)) &
                        res = .true.                    
                end if

            else ! boundary face
                fcs1 = GetCellFaceGA(grid%cell, cvs(1))
                if (size(fcs1) == 3) then
                   if (sum(grid%face%aligned%Get(fcs1))== 1) res = .true.
                else if (size(fcs1) == 4) then
                   if (sum(grid%face%aligned%Get(fcs1))== 2) res = .true. 
                end if
            end if
        end if

    end function
    
    function GetCommonFace(cell, ic1, ic2) result(res)
        type(GACellUDT) :: cell
        integer(I8) :: ic1, ic2, res, i
        integer(I8), allocatable :: fcs1(:), fcs2(:)

        res = 0

        fcs1 = GetCellFaceGA(cell, ic1)
        fcs2 = GetCellFaceGA(cell, ic2)

        do i = 1, size(fcs1)
            if (any(fcs1(i) == fcs2)) then
                res = fcs1(i)
                return
            end if
        end do

        if (res == 0) call gdErrorHandler('GetCommonFace: cells have no common face')

    end function

    function HaveCommonFace(cell, ic1, ic2) result(res)
        type(GACellUDT) :: cell
        integer(I8) :: ic1, ic2, i
        integer(I8), allocatable :: fcs1(:), fcs2(:)
        logical :: res

        res = .false.

        fcs1 = GetCellFaceGA(cell, ic1)
        fcs2 = GetCellFaceGA(cell, ic2)

        do i = 1, size(fcs1)
            if (any(fcs1(i) == fcs2)) then
                res = .true.
                return
            end if
        end do

    end function

    function GetCommonVert(f, f1, f2) result(res)
        type(GAFaceUDT) :: f
        integer(I8) :: i, f1, f2, res, vx1(2), vx2(2)

        res = 0

        if (HaveCommonVert(f,f1,f2)) then
            vx1(1) = f%vert1%Get(f1) 
            vx1(2) = f%vert2%Get(f1) 
            vx2(1) = f%vert1%Get(f2) 
            vx2(2) = f%vert2%Get(f2) 

            do i = 1, 2
                if (vx1(i) == vx2(1)) then
                    res = vx2(1)
                elseif (vx1(i) == vx2(2)) then
                    res = vx2(2)
                end if
            end do
        else

            call gdErrorHandler('GetCommonVert: face1 and face2 have no common vertex')
        end if

    end function

    function HaveCommonVert(f, f1, f2) result(res)
        type(GAFaceUDT) :: f
        integer(I8) :: i, f1, f2, vx1(2), vx2(2)
        logical :: res
        res = .false.
        if ((f1 == 0) .or. (f2 == 0)) return

        vx1(1) = f%vert1%Get(f1) 
        vx1(2) = f%vert2%Get(f1) 
        vx2(1) = f%vert1%Get(f2) 
        vx2(2) = f%vert2%Get(f2) 

        do i = 1, 2
            if (any(vx1(i) == vx2)) res = .true.
        end do

    end function

    function TriangleArea(x0, y0, x1, y1, x2, y2) result(res)
        real(R8) :: x0, y0, x1, y1, x2, y2, res

        res = 0.5_R8 * abs( (x1-x0)*(y2-y0) - (y1-y0)*(x2-x0) )

    end function 

    function Norm0D(x0, y0) result(res)
        real(R8) :: x0, y0, res
        res = sqrt(x0**2 + y0**2)
    end function

    function Norm1D(x0, y0) result(res)
        real(R8) :: x0(:), y0(:)
        real(R8), allocatable :: res(:)
        res = sqrt(x0**2 + y0**2)
    end function

    subroutine CalcCentroid0DGA(grid, ic)
        class(GAGridUDT) :: grid
        integer(I8) :: ic
        integer(I8), allocatable :: vxs(:)
        vxs = GetCellVertGA(grid%cell, ic)
        call grid%cell%x%Set(ic, sum(grid%vert%x%Get(vxs)/real(size(vxs), kind=R8)))
        call grid%cell%y%Set(ic, sum(grid%vert%y%Get(vxs)/real(size(vxs), kind=R8)))

    end subroutine

    subroutine CalcCentroid1DGA(grid, cells)
        class(GAGridUDT) :: grid
        integer(I8) :: i, cells(:)
        integer(I8), allocatable :: vxs(:)
        do i = 1, size(cells)
            vxs = GetCellVertGA(grid%cell, cells(i))
            call grid%cell%x%Set(cells(i), sum(grid%vert%x%Get(vxs)/real(size(vxs), kind=R8)))
            call grid%cell%y%Set(cells(i), sum(grid%vert%y%Get(vxs)/real(size(vxs), kind=R8)))
        end do
    end subroutine

    subroutine CalcHpol0D(grid, ic, h_pol)

        ! Description
        !============
        ! Calculate poloidal length of a cell. 
        ! The poloidal length is defined as the distance between the two intersection
        ! point of the mean psi line of the cell

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(in)    :: grid
        integer(I8), intent(in)         :: ic
        real(R8), intent(out)           :: h_pol

        ! Auxiliary
        integer(I8) :: ii, v1, v2, i
        integer(I8), allocatable, dimension(:) :: vxs, fcs
        real(R8) :: psic, v1p, v2p, isx(4), isy(4), t0


        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        h_pol = 0
        
        ! Get vertices
        vxs = GetCellVertGA(c, ic)

        ! mean psi
        psic = 0.5_R8 * ( maxval(v%psi%Get(vxs)) + minval(v%psi%Get(vxs)))

        ! Search faces intersections with psic
        fcs = GetCellFaceGA(c, ic)
        isx = 0
        isy = 0
        ii = 1
        do i = 1, size(vxs)
            v1 = f%vert1%Get(fcs(i))
            v2 = f%vert2%Get(fcs(i))
            v1p = v%psi%Get(v1)
            v2p = v%psi%Get(v2)
            if (psic .gt. min(v1p, v2p) &
                .and. psic .lt. max(v1p, v2p)) then
                    t0 = (psic - v1p) / (v2p - v1p)
                    isx(ii) = v%x%Get(v1) + t0 * (v%x%Get(v2) - v%x%Get(v1))
                    isy(ii) = v%y%Get(v1) + t0 * (v%y%Get(v2) - v%y%Get(v1))
                    ii = ii + 1 
            end if
        end do

        h_pol = sqrt( (isx(2) - isx(1))**2 + (isy(2) - isy(1))**2)

        end associate

    end subroutine

    subroutine CalcHrad0D(grid, ic, h_rad)

        ! Description
        !============
        ! Calculate radial length of a cell
        
        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(in)    :: grid
        integer(I8), intent(in)         :: ic
        real(R8), intent(out)           :: h_rad

        ! Auxiliary
        integer(I8) :: i, vs(2)
        integer(I8), allocatable, dimension(:) :: fcs, verts1, verts2
        real(R8) :: dx, dy, ds, vec_v1x, vec_v1y, d1, c1, dist

        ! Associate
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )
        
        fcs = GetCellFaceGA(c, ic)
        verts1 = f%vert1%Get(fcs)
        verts2 = f%vert2%Get(fcs)

        do i = 1, size(fcs)
            vs = [verts1(i), verts2(i)]

            ! Face vector
            dx = v%x%Get(vs(2)) - v%x%Get(vs(1))
            dy = v%y%Get(vs(2)) - v%y%Get(vs(1))
            ds = Norm(dx, dy)

            ! Magnetic field vector
            vec_v1x = v%bx%Get(vs(1))
            vec_v1y = v%by%Get(vs(1))
            d1 = Norm(vec_v1x, vec_v1y)

            ! Project distance
            c1 = abs(dx*vec_v1x + dy*vec_v1y) / (ds*d1)
            dist = c1*ds

            h_rad = h_rad + dist
        end do

        ! Gone all around the cell 
        h_rad = h_rad / 2


        end associate

    end subroutine

    subroutine CheckUniqueness(ida)
        type(IntegerDynamicArrayBufferedUDT) :: ida
        integer(I8), allocatable :: ida0(:), &
            ida_loc(:), ida0U(:)

        ida_loc = ida%Get()
        allocate(ida0(count(ida_loc /= 0)))
        ida0 = pack(ida_loc,ida_loc /= 0)
        call Unique(ida0,ida0U)
        if (size(ida0) /= size(ida0U)) &
            call gdErrorHandler('CheckUniqueness: integer array not unique')

    end subroutine

    subroutine WriteArray(a, filename)

        ! Description
        !============
        ! Write an array in a file

        ! Declare variables
        !==================
        ! Modules 
        use mod_plotter 
        use mod_specialchars, only : filesepchar

        ! Arguments
        integer(I8), allocatable :: a(:)
        character(*), intent(in) :: filename 

        ! Auxiliary
        integer :: fu     
        integer(I8) :: i
        character(:), allocatable :: dir

        ! Construct writing directory
        dir = plotdir // filesepchar // filename // '.dat'

        ! Open file
        open (action='write', file=trim(dir), newunit=fu, &
             status='unknown')

        ! Size data
        write (fu, *) 'Elements'
        write (fu, *) size(a)

        ! Array
        write (fu, *) 'ID val(ID)'
        do i = 1, size(a)
            write(fu, *) i, a(i)
        end do

        close(fu)

    end subroutine


end module 