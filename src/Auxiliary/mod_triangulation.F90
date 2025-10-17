!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains functionality for triangulation.

module mod_triangulation

    ! Load modules
    use mod_precision
    use mod_errorhandler
    use mod_gradient

    implicit none

    !==================================================================!
    !                                                                  !
    !                            TYPES                                 !
    !                                                                  !
    !==================================================================!
    
    ! Triangulation type
    type :: TriangulationUDT

        ! Description
        !============
        ! The triangulation type contains the following data describing a triangle mesh.
        ! - x, y    coordinates of the vertices (simple array)
        ! - cvert   connectivity of vertices to cells. The array has size
        !           (number of cells, 3). 

        real(R8), allocatable       :: x(:), y(:)
        integer(I8), allocatable    :: cvert(:,:)

    contains

        ! Constructors (unstructured)
        procedure :: ConstructTriaFromUnstructuredData
        generic :: Construct => ConstructTriaFromUnstructuredData

        ! Visualization
        procedure :: Visualize          =>  VisualizeTriangulation

    end type

    type, extends(GradientReconstructionUDT) :: GradientReconstructionTriaUDT

        ! Description
        !============
        ! Gradient reconstrunction implementation for triangulated grid

    contains

        ! Set parameters
        procedure :: SetParameters      => SetParametersGRTria

        ! Constructor
        procedure :: Construct          => ConstructGRTria

        ! Evaluate
        procedure :: Evaluate           => EvaluateGRTria

    end type

    contains

    !==================================================================!
    !                                                                  !
    !                             ROUTINES                             !
    !                                                                  !
    !==================================================================!
    
    !------------------------------------------------------------------!
    !                           CONSTRUCTORS                           !
    !------------------------------------------------------------------!
    
    subroutine ConstructTriaFromUnstructuredData(triangulation, xv, yv, vertlist, vertP1, vertP2)

        ! Description
        !============

        ! Declare variables
        !==================
        ! Arguments
        class(TriangulationUDT)                 :: triangulation
        real(R8), intent(in)                    :: xv(:), yv(:)
        integer(I8), intent(in)                 :: vertlist(:), vertP1(:), vertP2(:)


        ! Auxiliary
        integer(I8) :: i

        ! Checks
        if (size(vertP1) /= size(vertP2)) &
            call gdErrorHandler("ConstructFromUnstructuredData: vertP1 and vertP2, should have same size")
        if (size(xv) /= size(yv)) &
            call gdErrorHandler("ConstructFromUnstructuredData: xv, yv should have same size")
        if (.not.all(vertP2 == 3)) &
            call gdErrorHandler("ConstructFromUnstructuredData: grid does not only consists of triangular cells")

        ! Assign vertex coordinate
        triangulation%x = xv
        triangulation%y = yv

        ! Assign connecitivity
        allocate(triangulation%cvert(size(vertP1),3))
        do i = 1, size(vertP1)
            triangulation%cvert(i,:) = vertlist(vertP1(i) : vertP1(i) + vertP2(i) - 1)
        end do

    end subroutine

    subroutine VisualizeTriangulation(triangulation, filename)
        
        ! Description
        !============
        ! Writes out triangulation.dat file 

        ! 'vertices'
        ! <vert%ntot>
        ! 'ID, x, y'
        ! <ID, x, y>
        ! 'cell vertices'
        ! <ID, v1, v2, v3> 

        ! Declare variables
        !==================
        ! Modules 
        use mod_plotter 
        use mod_specialchars, only : filesepchar

        ! Arguments
        class(TriangulationUDT)                 :: triangulation
        character(*), intent(in)                :: filename

        ! Auxiliary
        integer                                 :: fu, nv, nc, i
        integer(I8), allocatable                :: cvert(:,:)
        real(R8), allocatable, dimension(:)     :: x, y
        character(:), allocatable               :: dir

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
        x = triangulation%x
        y = triangulation%y
        nv = size(x)

        ! Number of vertices
        write (fu, *) 'vertices'
        write (fu, *) nv 

        ! Vertex data
        write (fu, *) 'ID, x, y'
        do i = 1, nv 
            write (fu, *) i, x(i), y(i)
        end do 

        ! Write cell data
        !================
        ! Unpack
        cvert = triangulation%cvert
        nc = size(cvert, 1)

        ! Number of cells
        write (fu, *) 'cells'
        write (fu, *) nc

        ! Cell vertices
        write (fu, *) 'ID, v1, v2, v3'
        do i = 1, nc
            write (fu, *) i, cvert(i,1), cvert(i,2), cvert(i,3)
        end do 

        ! Housekeeping
        close(fu)

    end subroutine

    !------------------------------------------------------------------!
    !                     GRADIENT RECONSTRUCTION                      !
    !------------------------------------------------------------------!    

        subroutine SetParametersGRTria(GR, type1, type2, meth)

        ! Description
        !============
        ! Set parameters for gradient reconstruction

        ! Declare variables
        !==================
        ! Arguments
        class(GradientReconstructionTriaUDT)  :: GR
        character(:), allocatable           :: type1, type2, meth

        GR%type1 = type1
        GR%type2 = type2
        GR%meth = meth

    end subroutine

    subroutine ConstructGRTria(GR)

        ! Description
        !============
        ! Constructor 

        ! Declare variables
        !==================
        ! Arguments
        class(GradientReconstructionTriaUDT) :: GR

        call gdErrorHandler('ConstructGRTria: not implemented yet')

    end subroutine

    subroutine EvaluateGRTria(GR, v)

        ! Description
        !============
        ! Evaluator

        ! Declare variables
        !==================
        ! Arguments
        class(GradientReconstructionTriaUDT) :: GR
        real(R8), intent(in)               :: v(:)

        call gdErrorHandler('EvaluateGRTria: not implemented yet')

    end subroutine

end module