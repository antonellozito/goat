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
    use mod_utility

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
        ! - nv      number of vertices
        ! - nc      number of cells
        ! - vcell   cells of vertices
        ! - vcellP  pointer array for vcell

        real(R8), allocatable       :: x(:), y(:), dx(:,:), dy(:,:)
        integer(I8), allocatable    :: cvert(:,:)
        integer(I8)                 :: nv, nc
        integer(I8), allocatable    :: vcellP(:,:), vcell(:)

    contains

        ! Constructors (unstructured)
        procedure :: ConstructTriaFromUnstructuredData
        generic :: Construct => ConstructTriaFromUnstructuredData

        ! Stencil construction
        procedure :: ConstructStencil

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
        integer(I8) :: i, iv, counter, n
        integer(I8), allocatable :: query(:,:), vcellP(:,:), vcell(:), cvs(:)
        integer(I8), allocatable, dimension(:) :: v1, v2, v3

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

        ! Save number of vertices and cells
        triangulation%nv = size(triangulation%x)
        triangulation%nc = size(vertP1)

        ! Construct inverse of cvert
        allocate(query(triangulation%nc,3), vcellP(triangulation%nv,2), vcell(triangulation%nv*20))
        query = spread((/ (i, i = 1, triangulation%nc)/), 2, 3)
        counter = 0
        do iv = 1, triangulation%nv
            cvs = GetVertCellTriaQ(triangulation, iv, query)
            n = size(cvs)
            vcellP(iv,1) = counter + 1
            vcellP(iv,2) = n
            vcell(counter+1:counter+n) = cvs
            counter = counter + n

        end do
        triangulation%vcellP = vcellP
        triangulation%vcell = vcell(1:counter)

        ! Compute vectors
        allocate(triangulation%dx(triangulation%nc,3))
        allocate(triangulation%dy(triangulation%nc,3))
        v1 = triangulation%cvert(:,1)
        v2 = triangulation%cvert(:,2)
        v3 = triangulation%cvert(:,3)
        triangulation%dx(:,1) = triangulation%x(v2) - triangulation%x(v1)
        triangulation%dx(:,2) = triangulation%x(v3) - triangulation%x(v2)
        triangulation%dx(:,3) = triangulation%x(v1) - triangulation%x(v3)
        triangulation%dy(:,1) = triangulation%y(v2) - triangulation%y(v1)
        triangulation%dy(:,2) = triangulation%y(v3) - triangulation%y(v2)
        triangulation%dy(:,3) = triangulation%y(v1) - triangulation%y(v3)



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
        integer(I8) :: iv, j, n, n1, cNvP(size(tria%x),2), &
            vxs(200), counterv, counter, stencil_est, min_stencil, &
            info, min_stencil_loc
        integer(I8), allocatable :: cNv(:)
        real(R8), allocatable :: ATA_loc(:,:), distx(:), disty(:), ATA(:,:), w(:), weight(:)
        real(R8) :: t_start, t_end
        !real(R8) :: delta_t1, t_start1, t_end1, t_start2, t_end2, delta_t2
        logical :: int


        ! Set parameters depending on deriv
        ! Compute number of derivatives
        if (GR%meth /= 'global') call gdErrorHandler('ConstructGRTria: no other meth than "global" implemented')
        if (GR%deriv .gt. 6) call gdErrorHandler('ConstructGRTria: deriv > 6 not yet implement')
        print *, 'Constructing gradient reconstruction of ', GR%deriv,'th order on triangulated grid'
        n = sum((/(j, j = 2, GR%deriv + 1)/))
        min_stencil = n + 1
        stencil_est = min_stencil + 10
        int = .false.

        ! Timing
        call wall_time(t_start) 
        !delta_t1 = 0   
        !delta_t2 = 0       
            
        ! Initialize
        if (int) then
            n1 = n + 1
        else 
            n1 = n 
        end if
        allocate(cNv(tria%nv*stencil_est), ATA(tria%nv*stencil_est,n1), &
            w(tria%nv*stencil_est))
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
                do iv = 1, tria%nv

                    ! Timing
                    !call wall_time(t_start1)
                    info = 1

                    ! While loop in case solver could not converge
                    min_stencil_loc = min_stencil
                    do while (info /= 0)

                        ! ConstructStencil
                        call tria%ConstructStencil(iv, min_stencil_loc, vxs, counterv)
                        
                        ! Timing
                        !call wall_time(t_end1)
                        !delta_t1 = max(t_end1-t_start1, delta_t1)

                        ! Timing
                        !call wall_time(t_start2)

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

                        ! Compute coefficients
                        call ComputeATA(distx, disty, GR%deriv, int, ATA_loc, info)

                        ! Check convergence
                        if (info /= 0) then
                            print *, 'Could not converge at vertex: ', iv
                            print *, 'Trying with larger stencil'
                            call tria%WriteErrorData(vxs(1:counterv), 0)
                            min_stencil_loc = min_stencil_loc + 5
                            vxs = 0
                            counterv = 0
                        else 
                            ! Add to cNv
                            cNvP(iv, 1) = counter + 1
                            cNvP(iv, 2) = counterv
                            cNv(counter+1:counter+counterv) = vxs(1:counterv)

                            w(counter+1:counter+counterv) = weight 
                            if (allocated(weight)) deallocate(weight) 

                            ! Save coefficients
                            ATA(counter+1:counter+counterv,:) = ATA_loc
                            counter = counter + counterv                            
                        end if

                    end do 

                    ! Timing
                    !call wall_time(t_end2)
                    !delta_t2 = max(t_end2-t_start2, delta_t2)

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
        GR%coef = ATA(1:counter,:)
        GR%w = w(1:counter)

        ! Check
        if (size(GR%cNv) /= GR%cNvP(tria%nv,1)+GR%cNvP(tria%nv,2)-1) then
            call gdErrorHandler('ConstructGRTria: incompatible size of cNv and cNvP')
        end if

        ! Timing
        call wall_time(t_end)

        ! Display
        !print *, 'Time to construct GR coefficients: ', t_end - t_start, ' seconds'

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
        ! 5) d2phidxdy
        ! 6) d2phidy2
        ! 7) d3phidx3
        ! 8) d3phidx2dy
        ! 9) d3phidxdy2
        ! 10) d3phidy3
        ! 11) d4phidx4
        ! 12) d4phidx3dy
        ! 13) d4phidx2dy2
        ! 14) d4phidxdy3
        ! 15) d4phidy4
        ! 16) d5phidx5
        ! 17) d5phidx4dy
        ! 18) d5phidx3dy2
        ! 19) d5phidx2dy3
        ! 20) d5phidxdy4
        ! 21) d5phidy5
        ! 22) d6phidx6
        ! 23) d6phidx5dy
        ! 24) d6phidx4dy2
        ! 25) d6phidx3dy3
        ! 26) d6phidx2dy4
        ! 27) d6phidxdy5
        ! 28) d6phidy6

        ! Declare variables
        !==================
        ! Arguments
        class(GradientReconstructionTriaUDT) :: GR
        real(R8), intent(in)               :: v(:)
        real(R8), intent(out), allocatable :: deriv_vals(:,:)
        ! Auxiliary
        integer(I8) :: j, nv, iv, s, n, ne
        integer(I8), allocatable :: vxs(:)
        real(R8), allocatable :: b(:), c(:), coef(:,:)
        
        if (GR%deriv .gt. 6) call gdErrorHandler('ConstructGRTria: deriv > 6 not yet implement')
        ne = sum((/(j, j = 2, GR%deriv+1)/)) + 1

        select case (GR%type1)
        case ('vert')

            select case (GR%type2)
            case ('vert')

                nv = size(GR%cNvP,1)
                if (size(v) /= nv) call gdErrorHandler('EvaluateGRGA: incompatible v')
                allocate(deriv_vals(nv,ne), c(size(GR%coef,2)))
                deriv_vals = 0
                deriv_vals(:,1) = v
                do iv = 1, nv
                    s = GR%cNvP(iv,1)
                    n = GR%cNvP(iv,2)
                    vxs = GR%cNv(s:s+n-1)
                    coef = GR%coef(s:s+n-1,:)
                    b = GR%w(s:s+n-1) *(v(vxs) - v(iv))

                    c = matmul(transpose(coef), b)
                    deriv_vals(iv,2:ne) = c

                end do
            case default

                call gdErrorHandler('ConstructGRTria: type1 == vert, type2 not implemented')

            end select

        case default

            call gdErrorHandler('EvaluateGRTria: type1 not implemented')

        end select

    end subroutine

    subroutine ConstructStencil(tria, iv, min_stencil, vxs, counterv)

        ! Description
        !============
        ! Keep adding vertices to the stencil till minimal stencil size is reached

        ! Declare variables
        !==================
        ! Arguments
        class(TriangulationUDT), intent(in) :: tria
        integer(I8), intent(in)             :: iv, min_stencil
        integer(I8), intent(inout)          :: counterv, vxs(:)

        ! Auxiliary
        integer(I8) :: j, k, l, m, counter_loc
        integer(I8), allocatable, dimension(:) :: cvs, cvs2, tv, tv2
        logical, allocatable :: cvs_done(:)

        ! Construct query
        allocate(cvs_done(tria%nc))
        cvs_done = .false.

        ! Get loop over cells of vertex
        vxs = 0
        vxs(1) = iv
        counterv = 1
        cvs = GetVertCellTria(tria, vxs(1))
        do j = 1, size(cvs)

                ! Get vertices of cell
                tv = tria%cvert(cvs(j),:)
                do k = 1, 3
                    if (.not. any(tv(k) == vxs(1:counterv))) then
                        counterv = counterv + 1
                        vxs(counterv) = tv(k)
                    end if
                end do
        end do 
        cvs_done(cvs) = .true.
                    
        ! Need at least min_stencil data points for second order
        if (counterv .lt. min_stencil) then
            do j = 1, size(cvs)
                tv = tria%cvert(cvs(j),:)
                do k = 1, 3
                    cvs2 = GetVertCellTria(tria, tv(k))
                    do l = 1, size(cvs2)
                        if (.not.cvs_done(cvs2(l))) then
                            tv2 = tria%cvert(cvs2(l),:)
                            do m = 1, 3
                                if (.not. any(tv2(m) == vxs(1:counterv))) then
                                    counterv = counterv + 1
                                    vxs(counterv) = tv2(m)
                                    if (counterv .gt. min_stencil) exit
                                end if
                            end do 
                            cvs_done(cvs2(l)) = .true.
                        end if
                        if (counterv .gt. min_stencil + 5) exit
                    end do
                end do
            end do
        end if


        counter_loc = counterv
        do while (counter_loc .lt. min_stencil)

            ! Loop over current stencil
            do j = 2, counter_loc
                cvs = GetVertCellTria(tria, vxs(j))
                do l = 1, size(cvs)
                    if (.not.cvs_done(cvs(l))) then
                        tv = tria%cvert(cvs(l),:)
                        do m = 1, 3
                            if (.not. any(tv(m) == vxs(1:counterv))) then
                                counterv = counterv + 1
                                vxs(counterv) = tv(m)
                                if (counterv .gt. min_stencil) exit
                            end if
                        end do 
                        cvs_done(cvs(l)) = .true.     
                    end if
                    if (counterv .gt. min_stencil + 5) exit
                end do
            end do

            ! Update
            counter_loc = counterv

        end do

        if (counterv .lt. min_stencil) then
            call tria%WriteErrorData(vxs(1:counterv), 1)
            print *, 'Number of vertices: ', counterv
            print *, 'Minimum stencil: ', min_stencil
            print *, 'Vertex (iv): ', iv
            call gdErrorHandler('ConstructGRTria: stencil insufficient, implement this!')
        end if        

    end subroutine

    !==================================================================!
    !                                                                  !
    !                           FUNCTIONS                              !
    !                                                                  !
    !==================================================================!  
    
    ! Get cells of vertex not using inverse connectivity
    function GetVertCellTriaQ(tria, iv, query) result(res)
        type(TriangulationUDT) :: tria
        integer(I8) :: iv, j
        integer(I8), allocatable :: res(:)
        integer(I8), allocatable, optional :: query(:,:)
        logical, allocatable :: log(:,:)

        if (.not.present(query)) then

            allocate(query(tria%nc,3))
            query = spread((/ (j, j = 1, tria%nc)/), 2, 3)

        end if

        ! Get 
        log = tria%cvert == iv
        allocate(res(count(log)))
        res = pack(query, log)        

    end function

    function GetVertCellTria(tria, iv) result(res)
        type(TriangulationUDT) :: tria
        integer(I8) :: iv, s
        integer(I8), allocatable :: res(:)
        
        s = tria%vcellP(iv,1)
        res = tria%vcell(s:s+tria%vcellP(iv,2)-1)
    end function

end module