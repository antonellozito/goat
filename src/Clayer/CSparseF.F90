!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module provides the interface to routines in the C layer that
! interface to the CXSparse code of Tim Alden Davis. 


module CSparseF 

    ! Modules
    !========
    use, intrinsic :: iso_c_binding
    use CUtilities

    implicit none 
    public 

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    ! CSparse type for interfacing own CSparse type in interlayer
    type, bind(c) :: CSparseUDT 

        type(c_ptr) :: row
        type(c_ptr) :: col
        type(c_ptr) :: val
        integer(c_int) :: nrow
        integer(c_int)  :: ncol 
        integer(c_int)  :: nval
        

    end type 

    !==================================================================!
    !                                                                  !
    !                             INTERFACES                           !
    !                                                                  !
    !==================================================================!
    
    interface 
        ! Sparse Matrix Matrix muliplication 
        function SpMMF(A, B) result(C) &
            bind(c, name='SpMMV')
            use, intrinsic :: iso_c_binding
            import :: CSparseUDT
            type(CSparseUDT), value               :: A, B
            type(CSparseUDT)                     :: C
        end function

        ! Sparse Matrix Matrix multiplication, only int/doubles as input
        function SpMMFBasic(nr1, nc1, nv1, nr2, nc2, nv2, r1, c1, v1, r2, c2, v2) result(C) & 
            bind(c, name='SpMMBasic')
            use, intrinsic :: iso_c_binding
            import :: CSparseUDT 
            integer(c_int), value           :: nr1
            integer(c_int), value           :: nc1
            integer(c_int), value           :: nv1
            integer(c_int), value           :: nr2
            integer(c_int), value           :: nc2
            integer(c_int), value           :: nv2
            integer(c_int)                  :: r1(*) 
            integer(c_int)                  :: c1(*)
            real(c_double)                  :: v1(*)
            integer(c_int)                  :: r2(*)
            integer(c_int)                  :: c2(*)
            real(c_double)                  :: v2(*)
            type(CSparseUDT)                :: C
        end function 

        ! Sparse constructor
        function ConstructMyCSparseF(nrow, ncol, nval, row, col, val) result(A) &
            bind(c, name='ConstructMyCSparse')
            use, intrinsic :: iso_c_binding
            import :: CSparseUDT
            integer(c_int), intent(in)      :: row(*), col(*)
            real(c_double), intent(in)      :: val(*)
            integer(c_int), value           :: nrow, ncol, nval
            type(c_ptr)                     :: A
        end function

        function ParseTestF(A) result(B) &
            bind(c, name='ParseTest')
            use, intrinsic :: iso_c_binding
            import CSparseUDT 
            type(CSparseUDT), value     :: A
            type(CSparseUDT)            :: B 
        end function 

        function ParseTestFP2(A) result(B) &
            bind(c, name='ParseTestP2')
            use, intrinsic :: iso_c_binding
            import CSparseUDT 
            type(c_ptr), value          :: A
            type(c_ptr)                 :: B 
        end function 
    end interface 

    contains 

    !==================================================================!
    !                                                                  !
    !                              ROUTINES                            !
    !                                                                  !
    !==================================================================!

    ! Free memory usage of sparse matrix
    subroutine FreeCSparse(cs)

        ! Description
        !============
        ! Free the c-pointer data of a CSparseUDT 
        ! Important: it is assumed that these components are allocated!
        ! Otherwise, memory issues may occur... 

        ! Declare variables
        !==================
        ! Arguments
        type(CSparseUDT), intent(inout)     :: cs 

        ! Free memory
        !============
        call FreeCptr(cs%row)
        call FreeCptr(cs%col)
        call FreeCptr(cs%val)

    end subroutine

    

end module
