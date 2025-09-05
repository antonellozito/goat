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
        class(IntegerDynamicArrayUDT), allocatable    :: fieldlineID 
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

        ! Grid information utility
        procedure :: GetForbiddenMergeFaces
        procedure :: GetCutsXpoints
        procedure :: GetTangencyPoints
        procedure :: GetRadLineFaces
        procedure :: RecursiveGridMarching
        procedure :: DetermineCflags
        procedure :: AreaConstraintPents

        ! Grid operation
        procedure :: RemoveCells
        procedure :: RemoveFaces
        procedure :: RemoveVertices
        procedure :: GetFaceNumber
        procedure :: AddFaceToFsFc
        procedure :: AddCell

        ! Computing
        !===========
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

        vx1 = f%vert1%GetAllElements()
        vx2 = f%vert2%GetAllElements()
        vx = v%x%GetAllElements()
        vy = v%y%GetAllElements()
        vpsi = v%psi%GetAllElements()
        cx = c%x%GetAllElements()
        cy = c%y%GetAllElements()

        ! Calculate face properties
        v1x = v%x%GetMultipleElements(vx1)
        v1y = v%y%GetMultipleElements(vx1)
        v2x = v%x%GetMultipleElements(vx2)
        v2y = v%y%GetMultipleElements(vx2)

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
        type(GAGridUDT), intent(in)             :: grid
        type(GAoptionsUDT), intent(in)          :: options
        type(MagneticFieldUDT), intent(in)      :: magneticField
        logical :: select_split, select_merge


        ! Calculate cv metric
        call qm%ComputeQM(grid,options,magneticField)


        ! Selecting splitting cell
        if (select_split) then
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
        integer(I8) :: i, j, neig, nf
        integer(I8), allocatable, dimension(:) :: forbidden_fcs, &
            indsort, cells, small_cells, fcs, cvs, cvLookUp, indcv, &
            no_cells, trias, cells2, indfc, pol_faces, fcs1, fcs2
        real(R8) :: crit, dfunv(grid%cell%ntot), h_pol_no_cells_crit
        real(R8), allocatable, dimension(:) :: area_small_cells, &
            h_pol_no_cells_sorted, h_pol_cells, h_pol_cvs, bias
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
            ! Determine merge face
            ! Aligned face in other flux surface, or aligned faces have no common vertex



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

            call gdErrorHandler('SelectMergingFace: merge criterium minimal grid not implemented')

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
            call Sort(bias, indsort, .false.)
            pol_faces = pol_faces(indsort)
            if (bias(1) .gt. options%merge_bias_limit) then
                qm%merge_fc = pol_faces(1)
            end if

        case ('pol_flux')

            ! Distance function fun_r
            ! Only use cells in SOL
            log = (c%reg%Get() == 2)
            allocate(cells(count(log)))
            cells = pack(indcv, log)

            allocate(pol_fluxdens_est(size(cells)))
            call grid%fun_r%distr%Evaluate(c%x%Get(cells), c%y%Get(cells), pol_fluxdens_est)

            log = (pol_fluxdens_est .lt. 0.5_R8)
            



        case ('h_rad')
        case ('h_rad_core')
        case ('bias_rad_farSOL')
        case ('bias_rad')
        case ('skew_tria')
        case ('manual')

            call gdErrorHandler('SelectMergingFace: manual merging via input not possible in precompile code')
        
        case default

            call gdErrorHandler('SelectMergingFace: merge criterium not implemented')

        end select












        end associate

    end subroutine


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
        class(GAGridUDT) :: grid
        integer(I8), intent(in) ::  cell, faceA, faceB, faceC
        type(GAoptionsUDT) :: options

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
        class(GAGridUDT), intent(inout)     :: grid
        type(QualityMetricUDT), intent(in)  :: qm
        type(GAoptionsUDT)                  :: options 
        integer(I8), allocatable, intent(out) :: cctria(:), cctraps(:), cctrapsP(:,:)
        
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
        class(GAGridUDT), intent(in) :: grid
        integer(I8), intent(in) :: ic
        integer(I8), allocatable, intent(in) :: cvLookUp(:)
        integer(I8), intent(out) :: trap1

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
        class(GAGridUDT), intent(in) :: grid
        integer(I8), intent(in) :: ic
        integer(I8), allocatable, intent(in) :: cvLookUp(:)
        integer(I8), intent(inout) :: traps(:)
        integer(I8), intent(inout) :: counter
 
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
        class(GAGridUDT), intent(in) :: grid
        integer(I8), intent(in) :: tria
        integer(I8), allocatable, intent(in) :: traps(:)
        integer(I8), intent(out) :: con_vert

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
        class(GAGridUDT), intent(in) :: grid
        integer(I8), intent(in) :: tria, con_vert
        integer(I8), allocatable, intent(in) :: traps(:)
        logical, intent(out) :: approved

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
        integer(I8), intent(in) :: tria, con_vert
        integer(I8), allocatable, intent(in) :: traps(:)

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
            v1_ffbz, v1_bx, v1_by, fcA_length, fcA_length_int, &
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
            v1_ffbz = v%ffbz%Get(start_vertex)*inVdistribution + v%ffbz%Get(end_vertex)*Vdistribution
            allocate(v1_bx(size(v1_nx)))
            allocate(v1_by(size(v1_nx)))
            call magneticField%interp%Evaluate(v1_nx, v1_ny, 1, 0, v1_bx)
            call magneticField%interp%Evaluate(v1_nx, v1_ny, 0, 1, v1_by)

            call v%x%Append(v1_nx)
            call v%y%Append(v1_ny)
            call v%psi%Append(v1_psi)
            call v%ffbz%Append(v1_ffbz)
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
                call grid%AddCell(new_faces, new_verts, ic)

                ! Adjust centroid
                call grid%CalcCentroidGA(c%ntot)

                ! Give region
                call c%reg%Set(c%ntot, regs(i))

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
        class(GAGridUDT), intent(inout) :: grid
        type(QualityMetricUDT), intent(inout)  :: qm
        type(GAoptionsUDT) :: options

        ! Auxiliary
        integer(I8) :: i, ne
        integer(I8), allocatable :: empty_surf(:), indfs(:), cvLookUp(:), &
            cells(:)
        logical, allocatable :: b_flag(:)
       

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
        ne = count(grid%data%fluxdata%fluxsurfacefacesP2%Get() == 0)
        allocate(empty_surf(ne))
        indfs = (/ (i, i= 1, grid%data%fluxdata%fluxsurfacefacesP2%Size() )/)
        empty_surf = pack(indfs, grid%data%fluxdata%fluxsurfacefacesP2%Get() == 0)

        call grid%data%fluxdata%fluxsurfacefacesP1%Remove(empty_surf)
        call grid%data%fluxdata%fluxsurfacefacesP2%Remove(empty_surf)

        ! Update number of flux surfaces
        grid%data%fluxdata%nFs = grid%data%fluxdata%fluxsurfacefacesP1%Size()

        ! Check uniqueness
        call CheckUniqueness(grid%data%fluxdata%fluxsurfacefaces)


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
        integer(I8), intent(in) :: fc
        logical :: starter, pent_to_tria, special_case

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
        class(GAGridUDT), intent(inout) :: grid
        integer(I8), intent(in) :: fc
        logical, intent(in) :: starter, pent_to_tria, special_case
        character(:), allocatable, intent(out) :: caseID
        integer(I8), allocatable, intent(out) :: cvs(:)

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
        class(GAGridUDT) :: grid
        integer(I8) :: fc, cvs(2)

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
        call grid%AddCell(new_faces, new_verts, ic)
        call grid%cell%reg%Set(ic, grid%cell%reg%Get(cvs(1)))

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
        integer(I8) :: fc, cvs(2)
        logical :: starter

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
        integer(I8) :: fc, cvs(2)

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
        integer(I8) :: fc, cvs(2)

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
        call grid%AddCell(new_faces, new_verts, ic)
        call grid%cell%reg%Set(ic,grid%cell%reg%Get(cvs(1))) 
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
        call grid%AddCell(new_faces, new_verts, ic)

        call grid%cell%reg%Set(ic, grid%cell%reg%Get(cvs(1)))

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
        class(GAGridUDT) :: grid
        integer(I8), intent(in) :: three_vert, fc, cvs(2)
        integer(I8), intent(out) :: f1n, ic 
        integer(I8), allocatable, intent(out) :: fc23(:)

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
        call grid%AddCell(new_faces, new_verts, ic)
        call grid%cell%reg%Set(ic,grid%cell%reg%Get(cvs(1)))        

        
    end subroutine

    subroutine MakeNewThreeFace(grid, three_vert, fc, fc3, fc23, f1n)

        ! Description
        !============
        ! Make a new face from two faces connect with hanging node

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT) :: grid
        integer(I8), intent(in) :: three_vert, fc
        integer(I8), intent(out) :: f1n
        integer(I8), allocatable, intent(out) :: fc3(:), fc23(:)

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
        class(GAGridUDT) :: grid
        integer(I8), intent(in) :: three_vert, f1n, cvs(2)
        integer(I8), allocatable, intent(in) :: fc23(:)
        integer(I8), allocatable, intent(out) :: cells_rem(:)

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
            call grid%AddCell(new_faces_n, new_vertsN, ic)

            call grid%cell%reg%Set(ic, grid%cell%reg%Get(neigD(1)))

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
            fcs(counterf+1:counterf+counter1) = fcs1(1:counter1)
            counterf = counterf + counter1

            ! Add fcs2
            fcs(counterf+1:counterf+counter2) = fcs2(1:counter2)
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
        class(GAGridUDT), intent(in) :: grid
        integer(I8), allocatable, intent(in) :: b_faces(:)
        integer(I8), allocatable, intent(out) :: tang_points(:)

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
        class(GAGridUDT), intent(in) :: grid
        integer(I8), allocatable, intent(in) :: verts(:)
        integer(I8), allocatable, intent(out) :: faces(:)

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
        class(GAGridUDT)    :: grid
        type(GAoptionsUDT)  :: options
        real(R8) :: threshold
        integer(I8), allocatable, intent(out) :: cells(:)

        ! Auxiliary
        integer(I8) :: i
        integer(I8), allocatable :: indCv(:)
        real(R8), allocatable ::  val(:)

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
        integer(I8), allocatable :: cells(:)
        logical , optional :: b_flag(grid%face%ntot)

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
        integer(I8), allocatable, intent(in)    :: cells(:)

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
            call c%bp%Remove(cellsU)
            call c%bt%Remove(cellsU)
            call c%cflags%Remove(cellsU)
            call c%reg%Remove(cellsU)
            call c%ft%Remove(cellsU)

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
        integer(I8), allocatable, intent(in)    :: faces(:)

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
            call f%reg%Remove(facesU)
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
        integer(I8), allocatable, intent(in) :: verts(:)

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
            call v%ffbz%Remove(vertsU)

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
                if (face_num /= 0) exit
            end do

            ! Create new face if necessary
            if (meth == 2 .and. face_num == 0) then

                ! Set default face properties
                face_num = grid%face%ntot + 1

                call grid%face%vert1%Append(v1)
                call grid%face%vert2%Append(v2)
                call grid%face%label%Append(0)
                call grid%face%reg%Append(0)
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
            call grid%face%reg%Append(0)
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
        integer(I8), intent(in) :: new_f
        integer(I8), allocatable :: old_fs(:)

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

    subroutine AddCell(grid, faces, verts, ic)

        ! Description
        !============
        ! Add cell in connectivity cvFc and cvVx at the back.

        ! Declare variables
        !==================
        ! Arguments
        class(GAGridUDT), intent(inout) :: grid
        integer(I8), allocatable, intent(in) :: faces(:), verts(:)
        integer(I8), intent(out) :: ic

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
        call c%psi%Append(0.0_R8)
        call c%bp%Append(0.0_R8)
        call c%bt%Append(0.0_R8)
        call c%ft%Append(0) ! Give the correct number in postprocessing

        ! Calculate centroid
        call c%x%Append(sum(grid%vert%x%Get(verts))/real(nv, kind=R8))
        call c%y%Append(sum(grid%vert%y%Get(verts))/real(nv, kind=R8))




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
        ! 'ID, x, y, fieldlineID'
        ! <ID, x, y, fieldlineID>
        ! 'faces'
        ! <face%ntot> 
        ! 'ID, v1, v2, label, region'
        ! <ID, v1, v2, label, region>
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
        integer(I8), allocatable, dimension(:)  :: fID, v1, v2, region, &   
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
        fID = v%fieldlineID%Get()

        ! Number of vertices
        write (fu, *) 'vertices'
        write (fu, *) v%ntot 

        ! Vertex data
        write (fu, *) 'ID, x, y, fieldlineID'
        do i = 1, v%ntot 
            write (fu, *) i, x(i), y(i), fID(i)
        end do 

        ! Write face data
        !================
        ! Unpack
        v1 = f%vert1%Get()
        v2 = f%vert2%Get()
        region = f%reg%Get()
        label = f%label%Get()
        aligned = f%aligned%Get()

        ! Number of faces
        write (fu, *) 'faces'
        write (fu, *) f%ntot

        ! Face data
        write (fu, *) 'ID, v1, v2, label, region, aligned'
        do i = 1, f%ntot
            write (fu, *) i, v1(i), v2(i), label(i), region(i), aligned(i)
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

                    ! Housekeeping
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
        integer(I8), allocatable, optional :: cvLookUp(:)
        integer(I8), allocatable :: cvs(:), res(:), res_dummy(:)

        if (.not.present(cvLookUp)) then
            cvLookUp = GetCvLookUpGA(g%cell)
        end if

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
        integer(I8)                 :: i 
        type(GACellUDT)               :: cell
        integer(I8), allocatable, optional  :: cvLookUp(:)
        integer(I8), allocatable    :: res(:)  
        
        if (.not.present(cvLookUp)) then
            cvLookUp = GetCvLookUpGA(cell)
        end if
        allocate(res(count(cell%face%GetAllElements().eq.i)))
        res = pack(cvLookUp,cell%face%GetAllElements().eq.i)
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
        integer(I8), allocatable, optional    :: cvLookUp(:)
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
        integer(I8), allocatable :: tf(:)
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
        integer(I8), allocatable :: cell(:)
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
        integer(I8) :: i
        logical, allocatable :: res(:)
        integer(I8), allocatable :: fcs(:), verts(:)

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
        integer(I8) :: fcs
        logical :: res
        integer(I8), allocatable :: fcs1(:), fcs2(:), cvs(:), &
            faces(:), cvLookUp(:)

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
        integer(I8) :: i
        integer(I8), allocatable :: vxs(:), cells(:)
        do i = 1, size(cells)
            vxs = GetCellVertGA(grid%cell, cells(i))
            call grid%cell%x%Set(cells(i), sum(grid%vert%x%Get(vxs)/real(size(vxs), kind=R8)))
            call grid%cell%y%Set(cells(i), sum(grid%vert%y%Get(vxs)/real(size(vxs), kind=R8)))
        end do
    end subroutine

    subroutine CheckUniqueness(ida)
        type(IntegerDynamicArrayUDT) :: ida
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
        character(:), allocatable               :: dir

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