!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module provides functionality to generate 2D grids (potentially in 
! different ways) starting from a topological mesh. The main drivers 
! are made public and can be found under the 'drivers' subroutines tab. 

module ggmod_gridgeneration2D

    ! Load modules
    use mod_precision
    use mod_errorhandler
    use mod_dynamicarrays
    use mod_contour2D
    use mod_polygon
    use mod_sort
    use mod_definitions, only: TMvertexbndID, TMvertexmaxID, &
        TMvertexminID, TMvertexsaddleID, TMvertextp1ID, &
        TMvertextp2ID, TMfacepolID, TMfaceradID, TMfacebndID, &
        TMvertexsplitID, TMfacesepID
    use mod_linearsolverinterface, only: SolveDenseLinearSystemDI
    use goatmod_types, only : magneticFieldUDT, VesselUDT, GridUDT
    use goatmod_userinput, only : GGoptionsUDT
    use ggmod_topology2D
    use ggmod_vertexdistribution2D
    use DistributionFunction
    implicit none
    private 
    public :: GenerateUnstructuredAlignedGrid

    ! Module parameters
    real(R8), parameter, private        :: tprelfieldtol = 1e-10 ! relative field tolerance under which extrema are removed
    real(R8), parameter, private        :: disttol = 1e-12 ! distance tolerance

    !==================================================================!
    !                                                                  !
    !                            TYPES                                 !
    !                                                                  !
    !==================================================================!

    ! Grid data related to topological mesh
    !======================================
    ! Vertex data
    type :: GGTMVertexDataUDT

    contains 

        ! Initializer
        procedure :: Initialize         => InitializeGGTMVertexData

    end type 
    
    ! Face data
    type :: GGTMFaceDataUDT

        ! Description
        !============
        ! Contains additional face data, including the grid vertices
        ! xv, yv and the field values at these vertices, fv.
        type(RealDynamicArrayUDT), allocatable  :: xv(:), yv(:), fv(:) 

    contains 

        ! Initializer
        procedure :: Initialize         => InitializeGGTMFaceData

    end type 

    ! Cell data
    type :: GGTMCellDataUDT

    contains 

        ! Initializer
        procedure :: Initialize         => InitializeGGTMCellData

    end type

    ! Tube data
    type :: GGTMTubeDataUDT

        ! Description
        !============
        ! Contains additional flux tube data, such as field value 
        ! data etc 

        type(RealDynamicArrayUDT), allocatable  :: fval(:)
        integer(I8), allocatable, dimension(:)  :: distributionface
    
    contains 

        ! Initializer
        procedure :: Initialize         => InitializeGGTMTubeData

    end type 


    ! Topological mesh grid generator data
    type :: GGTMDataUDT

        ! Description
        !============
        ! This type contains additional information related to the 
        ! topological mesh to construct a grid. 

        type(GGTMVertexDataUDT)        :: vert
        type(GGTMFaceDataUDT)          :: face 
        type(GGTMCellDataUDT)          :: cell 
        type(GGTMTubeDataUDT)          :: tube 

    contains 

        ! Initializer
        procedure :: Initialize         => InitializeGGTMData

    end type 


    contains 

    !==================================================================!
    !                                                                  !
    !                           ROUTINES                               !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                     GRID GENERATION DRIVERS                      !
    !------------------------------------------------------------------!

    ! Unstructured aligned grid generator
    subroutine GenerateUnstructuredAlignedGrid(topomesh, magneticField, &
        vessel, fieldtracer, boundarytracer, options)

        ! Description
        !============
        ! Generate a grid from scratch for a topological mesh given in topomesh. It
        ! is assumed that each cell in the topological mesh has exactly two
        ! 'radial' boundaries (type 1) or one radial and one vessel boundary (type
        ! 3), and an arbitrary number of poloidal boundaries. First, we construct
        ! based on this topology mesh the tubes (list of cell and faces) through
        ! which one bundle of field coordinates should pass. Then, we distribute
        ! vertices on all faces of the topological mesh, where we assume that for
        ! each 'poloidal' face, we have the freedom to distribute however we like
        ! regardless of the other faces, but that the distribution of the radial
        ! faces is determined based on the tubes. Afterwards, we trace the field
        ! lines of each tube. Then, we are ready to distribute for each grid cell
        ! the grid nodes, which should be based on a given distribution (determined
        ! through the options). Finally, we connect all cells and vertices together
        ! and call all other necessary grid initialization routines. 

        ! Notes
        !======
        ! Note 1: all vertices of the topological mesh will become vertices of the
        ! grid

        ! Note 2: vertex distributions of each face are added on the topological
        ! mesh level (facedata)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)          :: topomesh 
        type(MagneticFieldUDT)      :: magneticField
        type(VesselUDT)             :: vessel 
        class(ContourTracerUDT)     :: fieldtracer, boundarytracer 
        type(GGoptionsUDT)          :: options 

        ! Auxiliary
        type(GGTMDataUDT)           :: ggtmdata 
        class(VertexDistributor2DUDT), allocatable      :: &
            poloidalvertexdistributor, radialvertexdistributor
        class(DistributionFunctionUDT), allocatable     :: & 
            magneticFieldDF 

        ! Initialize
        !===========
        ! Required data of topomesh for grid generator
        call ggtmdata%Initialize(topomesh)

        ! Magnetic field distribution function
        magneticFieldDF = ConstructStructured2DDF(magneticField%interp)

        ! Initial grid structure
        !vert = struct('x', zeros(0, 1), 'y', zeros(0, 1), 'BV', zeros(0, 1), 'ntot', 0);
        !face = struct('vert', zeros(0, 2) , 'ntot', 0, 'labels', zeros(0, 1), 'region', zeros(0, 1));
        !cell = struct('vert', zeros(0, 1), 'vertP', zeros(0, 2), 'region', zeros(0, 1), 'ntot', 0, 'nvert', 0);
        !fs = struct('ntot', 0);
        !grid = struct('vert', vert, 'face', face, 'cell', cell, 'fs', fs);

        ! Set up vertex distribution
        !===========================
        ! Poloidal vertex distributor
        select case (options%vdptype)

        case ('uniform')

            ! Construct uniform distributor with facelength 'options%vdpdfacelength'
            poloidalvertexdistributor = ConstructUniformVertexDistributor(&
                options%vdpdfacelength, options%vdrdfieldwidth)

        case default 

            ! Unknown option
            call gdErrorHandler('GenerateUnstructuredAlignedGrid: ' // & 
                'poloidal vertex distribution option: ' // options%vdptype // &
                ' not implemented')

        end select

        ! Radial vertex distributor
        select case (options%vdrtype)

        case ('uniform')

            ! Construct uniform distributor 
            radialvertexdistributor = ConstructUniformVertexDistributor(&
                options%vdpdfacelength, options%vdrdfieldwidth)

        case Default

            ! Unknown option
            call gdErrorHandler('GenerateUnstructuredAlignedGrid: ' // & 
                'radial vertex distribution option: ' // options%vdrtype // &
                ' not implemented')

        end select

        !pdoptions = options.poloidaldistributor;
        !pdoptions = SetPoloidalDistributorOptions(pdoptions);
        !pdoptions.distribution.distrfield = field;
        !pdoptions.distribution.bnd = bnd;
        !poloidaldistributor = ConstructVertexDistributor(pdoptions);

        ! Extract x-point data
        !isxp = cat(1, topomesh.vertdata.type) == 2;
        !xpdata = cat(1, topomesh.vertdata(isxp).F);

        ! Radial
        !rdoptions = options.radialdistributor;
        !rdoptions = SetRadialDistributorOptions(rdoptions);
        !rdoptions.distribution.distrfield = field;
        !rdoptions.distribution.bnd = bnd;
        !rdoptions.distribution.xpointdata = xpdata;
        !radialdistributor = ConstructFieldDistributor(rdoptions, field);

        ! Distribute vertices on topological faces
        !=========================================
        ! Poloidal faces
        call DistributeVerticesTopologicalMeshFaces(ggtmdata, topomesh, &
            poloidalvertexdistributor, [TMfacepolID, TMfacesepID])

        ! Radial faces and tubes
        call DistributeFieldValuesTopologicalMeshTubes(ggtmdata, &
            topomesh, radialvertexdistributor, magneticFieldDF, [TMfaceradID, TMfacebndID])

        ! Generate grid
        !--------------
        ![grid] = DistributeVerticesTopologicalMeshCells(grid, topomesh, tubes, poloidaldistributor, field, bnd, options.gridding);

        ! Add additional interconnections
        !--------------------------------
        ![grid] = ComputeGridInterconnections(grid);

        ! Visualize
        !----------
        !VisualizeGrid(grid, 1);


    end subroutine

    !------------------------------------------------------------------!
    !                  TOPOMESH VERTEX DISTRIBUTION                    !
    !------------------------------------------------------------------!

    ! Distribution of vertices along topological mesh faces

    ! Distribution of field values on topological flux tubes

    !------------------------------------------------------------------!
    !                  TOPOMESH GRID DATA ROUTINES                     !
    !------------------------------------------------------------------!

    ! Initializers
    subroutine InitializeGGTMVertexData(vert, topomesh)

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMVertexDataUDT)            :: vert 
        class(TopomeshUDT),intent(in)       :: topomesh
        
        ! Initialize
        !===========

    end subroutine

    subroutine InitializeGGTMFaceData(face, topomesh)

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMFaceDataUDT)              :: face 
        class(TopomeshUDT),intent(in)       :: topomesh
        
        ! Initialize
        !===========
        ! Associate
        associate(nf        => topomesh%face%ntot)

        ! Allocate coordinates
        allocate(face%xv(nf), face%yv(nf))

        ! Housekeeping
        !=============
        end associate
    end subroutine 

    subroutine InitializeGGTMCellData(cell, topomesh)

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMCellDataUDT)              :: cell 
        class(TopomeshUDT),intent(in)       :: topomesh
        
        ! Initialize
        !===========

    end subroutine 

    subroutine InitializeGGTMTubeData(tube, topomesh)

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMTubeDataUDT)              :: tube 
        class(TopomeshUDT),intent(in)       :: topomesh
        
        ! Initialize
        !===========
        ! Associate
        associate(nt        => topomesh%tube%ntot)

        ! Allocate
        allocate(tube%fval(nt), tube%distributionface(nt))

        ! Initial value where applicable
        tube%distributionface = 0

        ! Housekeeping
        end associate

    end subroutine 

    subroutine InitializeGGTMData(ggtmdata, topomesh)

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMDataUDT)                  :: ggtmdata 
        class(TopomeshUDT),intent(in)       :: topomesh
        
        ! Initialize
        !===========
        ! Substructures
        call ggtmdata%vert%Initialize(topomesh)
        call ggtmdata%face%Initialize(topomesh)
        call ggtmdata%cell%Initialize(topomesh)
        call ggtmdata%tube%Initialize(topomesh)

    end subroutine 

    ! Vertex distribution over faces
    subroutine DistributeVerticesTopologicalMeshFaces(ggtmdata, topomesh, &
        vd, facetypes)

        ! Description
        !============
        ! Distribute vertices over topological mesh faces of the types
        ! defined in 'facetypes'. The vertices etc are stored in the
        ! ggtmdata structure. The vertex distribution is done based on 
        ! the vertexdistributor that is passed 

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMDataUDT)                      :: ggtmdata 
        class(TopomeshUDT), intent(in)          :: topomesh 
        class(VertexDistributor2DUDT), intent(in)   :: vd 
        integer(I8), intent(in)                 :: facetypes(:)

        ! Auxiliary
        integer(I8)                             :: nv 
        real(R8), allocatable, dimension(:)     :: tx, ty 

        ! Loop
        integer(I8)                             :: i 

        ! Initialize
        !===========
        ! Associate
        associate(&
            vert        => topomesh%vert,   &
            face        => topomesh%face,   &
            facedata    => ggtmdata%face    &
            )

        ! Loop over all faces and distribute
        !===================================
        do i = 1, face%ntot
            if (any(face%type(i) == facetypes)) then 
                ! Distribute
                call vd%DistributeOverCurve(face%x(i)%Get(), face%y(i)%Get(), tx, ty, nv)
                
                ! Adjust start and end to be sure
                tx(1) = vert%x(face%vert(i, 1))
                ty(1) = vert%y(face%vert(i, 1))
                tx(nv) = vert%x(face%vert(i, 2))
                ty(nv) = vert%y(face%vert(i, 2))
                
                ! Add
                facedata%xv(i) = ConstructRealDynamicArray(tx)
                facedata%yv(i) = ConstructRealDynamicArray(ty)
                !facedata%faceflags(i) = ConstructIntegerDynamicArray(spread())
                !facedata(i).faceflags = ones(numel(facedata(i).vx)-1, size(face.flags, 2)).*face.flags(i, :);
            end if 
        end do 

        ! Housekeeping
        !=============
        end associate




    end subroutine 

    ! Field value distribution over tubes
    subroutine DistributeFieldValuesTopologicalMeshTubes(ggtmdata, topomesh, &
        vd, field, facetypes)

        ! Description
        !============
        ! Distribute field values over topological mesh faces of the types
        ! defined in 'facetypes'. The field values are computed per tube,
        ! first based on an initial distribution on the faces in 
        ! 'facetypes'. This chosen distribution is simply 
        ! taken as the distribution that gives the maximal amount 
        ! of faces, as this is likely the desired one. We store the 
        ! chosen topological face ID for each flux tube (we don't 
        ! propagate the distribution to other faces yet, since we have 
        ! to compute intersections anyway)

        ! Declare variables
        !==================
        ! Arguments
        class(GGTMDataUDT)                      :: ggtmdata 
        class(TopomeshUDT), intent(in)          :: topomesh 
        class(VertexDistributor2DUDT), intent(in)    :: vd 
        class(DistributionFunctionUDT), intent(in)  :: field 
        integer(I8), intent(in)                 :: facetypes(:)

        ! Auxiliary
        integer(I8)                             :: nv, nfl, nflmax, &
            tfmax
        integer(I8), allocatable, dimension(:)  :: tf
        real(R8), allocatable, dimension(:)     :: fc, xc, yc, tx, &
            ty

        ! Loop
        integer(I8)                             :: i, j  

        ! Initialize
        !===========
        ! Associate
        associate(&
            vert        => topomesh%vert,   &
            face        => topomesh%face,   &
            facedata    => ggtmdata%face,   &
            tube        => topomesh%tube,   &
            tubedata    => ggtmdata%tube    &
            )

        ! Distribute over faces
        !======================
        do i = 1, face%ntot 
            if (any(face%type(i) == facetypes)) then 
                ! Unpack
                xc = face%x(i)%Get()
                yc = face%y(i)%Get()
                allocate(fc(face%x(i)%Size()))

                ! Distribute
                call vd%DistributeOverField(xc, yc, field, tx, ty, nv)
                tx(1) = vert%x(face%vert(i, 1))
                ty(1) = vert%y(face%vert(i, 1))
                tx(nv) = vert%x(face%vert(i, 2))
                ty(nv) = vert%y(face%vert(i, 2))

                ! Add data
                facedata%xv(i) = ConstructRealDynamicArray(tx)
                facedata%yv(i) = ConstructRealDynamicArray(ty)

            end if 
        end do 

        ! Determine tube distribution
        !============================
        do i = 1, tube%ntot 
            ! Get tube faces
            tf = tube%GetFace(i)

            ! Initialize
            nfl = 0
            nflmax = 0

            ! Loop
            do j = 1, size(tf)
                ! Determine number of field lines
                nfl = facedata%fv(tf(j))%Size()
                if (nfl > nflmax) then 
                    nflmax = nfl 
                    tfmax = tf(j)
                end if 
            end do 

            ! Add maximal distribution to tube
            tubedata%distributionface(i) = tfmax

        end do 

        ! Overwrite
        !==========

        ! Housekeeping
        !=============
        end associate

    end subroutine 


end module 