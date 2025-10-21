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
    use mod_readwrite
    use mod_gradient
    use mod_sort

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

        ! Stencil construction
        procedure :: AddVxsToStencil

        ! Visualization
        procedure :: Visualize          =>  VisualizeTriangulation
        procedure :: WriteErrorData     =>  WriteErrorDataTria

    end type

    type, extends(GradientReconstructionUDT) :: GradientReconstructionTriaUDT

        ! Description
        !============
        ! Gradient reconstrunction implementation for triangulated grid

        ! - deriv: indicate order of derivative required

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

    subroutine WriteErrorDataTria(tria, verts, flag)

        ! Description
        !============
        ! Write data for an error plot        
        ! Declare variables
        !==================
        ! Arguments
        class(TriangulationUDT), intent(in) :: tria
        integer(I8), intent(in)         :: verts(:), flag

        ! Write grid and array for vertices to indicate problematic area
        call tria%Visualize('tria_error')
        call WriteArray(verts, 'vertices_error')

        if (flag == 1) print *, 'Error occurs: use pgaerror to visualize'

    end subroutine

    !------------------------------------------------------------------!
    !                     GRADIENT RECONSTRUCTION                      !
    !------------------------------------------------------------------!    

    subroutine SetParametersGRTria(GR, type1, type2, meth, deriv)

        ! Description
        !============
        ! Set parameters for gradient reconstruction

        ! Declare variables
        !==================
        ! Arguments
        class(GradientReconstructionTriaUDT)    :: GR
        character(*)                            :: type1, type2, meth
        integer(I8)                             :: deriv

        GR%type1 = type1
        GR%type2 = type2
        GR%meth = meth
        GR%deriv = deriv

    end subroutine

    subroutine ConstructGRTria(GR, tria)

        ! Description
        !============
        ! Constructor 

        ! Declare variables
        !==================
        ! Arguments
        class(GradientReconstructionTriaUDT)    :: GR
        type(TriangulationUDT), intent(in)      :: tria

        ! Auxiliary
        integer(I8) :: iv, j, k, l, m, n, cNvP(size(tria%x),2), &
            vxs(200), counterv, counter, tv(3), tv2(3), stencil_est, min_stencil
        integer(I8), allocatable :: cvs(:), cvs2(:), cNv(:), ar(:)
        real(R8), allocatable :: ATA_loc(:,:), distx(:), disty(:), ATA(:,:), w(:), weight(:)


        ! Set parameters depending on deriv
        ! Compute number of derivatives
        if (GR%meth /= 'global') call gdErrorHandler('ConstructGRTria: no other meth than "global" implemented')
        if (GR%deriv .gt. 6) call gdErrorHandler('ConstructGRTria: deriv > 6 not yet implement')
        ar = (/(j, j = 2, 7)/)
        n = sum(ar(1:GR%deriv))
        min_stencil = n + 1
        stencil_est = min_stencil + 10
        !if (GR%deriv == 1) then
        !    stencil_est = 15
        !else if (GR%deriv == 2) then
        !    stencil_est = 20
        !else if (GR%deriv == 3) then
        !    stencil_est = 25
        !else if (GR%deriv == 4) then
        !    stencil_est = 30
        !else if (GR%deriv == 5) then
        !    stencil_est = 40
        !else if (GR%deriv == 6) then
        !    stencil_est = 50
        !end if

           
            
        ! Initialize
        allocate(cNv(size(tria%x)*stencil_est), ATA(size(tria%x)*stencil_est,n), w(size(tria%x)*stencil_est))
        counter = 0
        cNv = 0
        cNvP = 0
        ATA = 0
        w = 0  
        select case (GR%type1)
        case ('vert')

            select case (GR%type2)
            case ('vert')

                ! Loop over vertices
                do iv = 1, size(tria%x)

                    ! Reset
                    vxs = 0

                    ! Get cells
                    vxs(1) = iv
                    counterv = 1
                    cvs = GetVertCellTria(tria, iv)
                    do j = 1, size(cvs)
                            tv = tria%cvert(cvs(j),:)
                            do k = 1, 3
                                if (.not. any(tv(k) == vxs(1:counterv))) then
                                    counterv = counterv + 1
                                    vxs(counterv) = tv(k)
                                end if

                            end do

                    end do 
                    
                    ! Need at least min_stencil data points for second order
                    if (counterv .lt. min_stencil) then
                        do j = 1, size(cvs)
                            tv = tria%cvert(cvs(j),:)
                            do k = 1, 3
                                cvs2 = GetVertCellTria(tria, tv(k))
                                do l = 1, size(cvs2)
                                    tv2 = tria%cvert(cvs2(l),:)
                                    do m = 1, 3
                                        if (.not. any(tv2(m) == vxs(1:counterv))) then
                                            counterv = counterv + 1
                                            vxs(counterv) = tv2(m)
                                            if (counterv .gt. min_stencil) exit
                                        end if
                                    end do 
                                    if (counterv .gt. min_stencil + 5) exit                                       
                                end do
                            end do
                        end do
                    end if

                    call tria%AddVxsToStencil(min_stencil, vxs, counterv)

                    if (counterv .lt. min_stencil) then
                        call tria%WriteErrorData(vxs(1:counterv), 1)
                        print *, 'Number of vertices: ', counterv
                        print *, 'Minimum stencil: ', min_stencil
                        print *, 'Vertex (iv): ', iv
                        call gdErrorHandler('ConstructGRTria: stencil insufficient, implement this!')
                    end if

                    ! Add to cNv
                    cNvP(iv, 1) = counter + 1
                    cNvP(iv, 2) = counterv
                    cNv(counter+1:counter+counterv) = vxs(1:counterv)    

                    ! Compute coefficients
                    distx = tria%x(vxs(1:counterv)) - tria%x(iv)
                    disty = tria%y(vxs(1:counterv)) - tria%y(iv)

                    ! Compute weight is necessary
                    if (GR%deriv .gt. 1) then
                        weight = 1/(distx**2 + disty**2)
                        weight(1) = maxval(weight(2:size(distx)))*2
                    else 
                        allocate(weight(size(distx)))
                        weight = 1
                    end if
                    w(counter+1:counter+counterv) = weight 
                    if (allocated(weight)) deallocate(weight) 

                    call ComputeATA(distx, disty, GR%deriv, .false., ATA_loc)
                    ATA(counter+1:counter+counterv,:) = ATA_loc
                    counter = counter + counterv

                end do

            case default

                    call gdErrorHandler('ConstructGRTria: type1 == vert, type2 not implemented')

            end select

        case default

            call gdErrorHandler('ConstructGRTria: type1 not implemented')

        end select        


        ! Save in GR type
        GR%cNv = cNv(1:counter)
        GR%cNvP = cNvP   
        GR%invA = ATA(1:counter,:)
        GR%w = w(1:counter)


    end subroutine

    subroutine EvaluateGRTria(GR, v, deriv_vals)

        ! Description
        !============
        ! Evaluator of gradient reconstruction on triangulated grid
        ! The resulting matrix deriv_vals contains the derivatives in the different 
        ! columns:
        ! 1) phi
        ! 2) dphidx
        ! 3) dphidy
        ! 4) d2phidx2
        ! 5) d2phidy2
        ! 6) d2phidxdy
        ! 7) d3phidx3
        ! 8) d3phidy3
        ! 9) d3phidx2dy
        ! 10) d3phidxdy2
        ! 11) d4phidx4
        ! 12) d4phidy4
        ! 13) d4phidx3dy
        ! 14) d4phidx2dy2
        ! 15) d4phidxdy3
        ! 16) d5phidx5
        ! 17) d5phidy5
        ! 18) d5phidx4dy
        ! 19) d5phidx3dy2
        ! 20) d5phidx2dy3
        ! 21) d5phidxdy4
        ! 22) d6phidx6
        ! 23) d6phidy6
        ! 24) d6phidx5dy
        ! 25) d6phidx4dy2
        ! 26) d6phidx3dy3
        ! 27) d6phidx2dy4
        ! 28) d6phidxdy5

        ! Declare variables
        !==================
        ! Arguments
        class(GradientReconstructionTriaUDT) :: GR
        real(R8), intent(in)               :: v(:)
        real(R8), intent(out), allocatable :: deriv_vals(:,:)
        ! Auxiliary
        integer(I8) :: j, nv, iv, s, n, ne
        integer(I8), allocatable :: vxs(:), ar(:)
        real(R8), allocatable :: b(:), c(:), coef(:,:)
        
        if (GR%deriv .gt. 6) call gdErrorHandler('ConstructGRTria: deriv > 6 not yet implement')
        ar = (/(j, j = 2, 7)/)
        ne = sum(ar(1:GR%deriv)) + 1

        select case (GR%type1)
        case ('vert')

            select case (GR%type2)
            case ('vert')

                nv = size(GR%cNvP,1)
                if (size(v) /= nv) call gdErrorHandler('EvaluateGRGA: incompatible v')
                allocate(deriv_vals(nv,ne), c(size(GR%invA,2)))
                deriv_vals = 0
                deriv_vals(:,1) = v
                do iv = 1, nv
                    s = GR%cNvP(iv,1)
                    n = GR%cNvP(iv,2)
                    vxs = GR%cNv(s:s+n-1)
                    coef = GR%invA(s:s+n-1,:)
                    b = GR%w(s:s+n-1) *(v(vxs) - v(iv))

                    c = matmul(transpose(coef), b)
                    deriv_vals(iv,2:ne) = c
                    !deriv_vals(iv,2) = c(1)
                    !deriv_vals(iv,3) = c(2)

                end do
            case default

                call gdErrorHandler('ConstructGRTria: type1 == vert, type2 not implemented')

            end select

        case default

            call gdErrorHandler('EvaluateGRTria: type1 not implemented')

        end select

    end subroutine

    subroutine AddVxsToStencil(tria, min_stencil, vxs, counterv)

        ! Description
        !============
        ! Keep adding vertices to the stencil till minimal stencil size is reached

        ! Declare variables
        !==================
        ! Arguments
        class(TriangulationUDT), intent(in) :: tria
        integer(I8), intent(in)             :: min_stencil
        integer(I8), intent(inout)          :: counterv, vxs(:)

        ! Auxiliary
        integer(I8) :: j, l, m, counter_loc
        integer(I8), allocatable :: cvs(:), tv(:)

        counter_loc = counterv
        do while (counter_loc .lt. min_stencil)

            ! Loop over current stencil
            do j = 2, counter_loc
                cvs = GetVertCellTria(tria, vxs(j))
                do l = 1, size(cvs)
                        tv = tria%cvert(cvs(l),:)
                        do m = 1, 3
                            if (.not. any(tv(m) == vxs(1:counterv))) then
                                counterv = counterv + 1
                                vxs(counterv) = tv(m)
                                if (counterv .gt. min_stencil) exit
                            end if
                        end do 
                        if (counterv .gt. min_stencil + 5) exit                                       
                end do
            end do

            ! Update
            counter_loc = counterv

        end do

    end subroutine

    !==================================================================!
    !                                                                  !
    !                           FUNCTIONS                              !
    !                                                                  !
    !==================================================================!  
    
    ! Get cells of vertex
    function GetVertCellTria(tria, iv) result(res)
        type(TriangulationUDT) :: tria
        integer(I8) :: iv, j, n1, n2, n3
        integer(I8), allocatable :: indc(:), res(:), &
            res1(:), res2(:), res3(:), res4(:), cv1(:), &
            cv2(:), cv3(:)

        cv1 = tria%cvert(:,1)
        cv2 = tria%cvert(:,2)
        cv3 = tria%cvert(:,3)

        allocate(res1(count(cv1 == iv)))
        allocate(res2(count(cv2 == iv)))
        allocate(res3(count(cv3 == iv)))
        indc = (/ (j, j = 1, size(cv1))/)
        res1 = pack(indc,cv1 == iv)
        res2 = pack(indc,cv2 == iv)
        res3 = pack(indc,cv3 == iv)
        n1 = size(res1)
        n2 = size(res2)
        n3 = size(res3)
        allocate(res4(n1+n2+n3))
        res4(1:n1) = res1 
        res4(n1+1:n1+n2) = res2 
        res4(n1+n2+1:n1+n2+n3) = res3 

        call Unique(res4, res)
    end function

end module