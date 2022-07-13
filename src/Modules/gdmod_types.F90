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

    ! Option types
    !=============
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


    ! Grid types
    !===========
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

    ! Main grid structure
    type GridUDT
        ! Vertices
        type(VertexUDT)                     :: vert

        ! Faces
        type(FaceUDT)                       :: faces

        ! Cells
        type(CellUDT)                       :: cells
    end type

    ! Optimization
    !=============
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

    subroutine SetDefaultRunfileOptions(options)
        ! Description
        !============
        ! Set the default runfile options

        ! Declaration
        type (RunfileOptionsUDT), intent(inout)    :: options

        ! Default options
        options%runtype     = 'optimize'
        options%gridtype    = 'plasma'
        options%meth        = 'KKT'
        options%export      = .true.  

    end subroutine

    ! Grid routines
    !==============
    ! Allocation
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


    end subroutine

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

    ! Deallocation


end module