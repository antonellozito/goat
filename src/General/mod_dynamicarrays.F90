!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module provides several implementations for arrays that 
! have to change size (i.e. dynamic arrays). These arrays have some
! basic operators, just as classic fortran arrays, but may not have 
! all functionality that basic arrays have. The idea of this module
! is to provide several derived dynamic array types with different 
! optimized routines depending on the use case at hand. The arrays 
! should be constructed using the dedicated array constructors. 

! Notes
!======
! Note 1: focus lies on 1D arrays (not explicitly mentioned in the name
! ) since 2D arrays are likely dealt with better by sparse matrices or
! coarrays. 

! Note 2: most (if not all) functionality of each derived type should 
! also be available for other derived types such that switching 
! between them in the code becomes as easy as changing the type and
! constructor name in the declaration/initialization part of the code.
! See also Note 3 on why this couldn't be done (yet) in a more generic
! way.

! Note 3: It is a great disappointment, but at the time of writing this 
! code, generic features such as parameterized types with type-bound
! procedures of which the type can be determined at runtime are still 
! unavailable in the fortran standard (let alone available in popular
! compilers such as gfortran/ifort). Therefore, we need to fix the 
! precision of the arrays in this modules (in this case R8). In the hope
! that in the future this constraint is alleviated, this dependency has
! been programmed as an integer module parameter (private to this 
! module) and having all the rest implemented as derived
! types with type-bound procedures where applicable. In the future, one
! might be able to remove the precision parameter and have it as a truly
! parameterized type. 

module mod_dynamicarrays
    
    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_errorhandler
    use mod_sparseinterface

    ! The usual
    implicit none
    save
    public

    ! Precision definition
    integer, parameter, private :: rk = R8
    integer, parameter, private :: ik = I8

    ! Memory MGMT settings
    integer, parameter, private :: grow = 2 ! factor by which to grow memory of array 
    integer, parameter, private :: init_array_size = grow**8 ! initial array size 

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    ! Real dynamic array class
    type :: RealDynamicArrayUDT

        ! Description
        !============
        ! This is the base dynamic array type for general purpose 
        ! dynamic array operations. Only for real values (integer has 
        ! its own UDT) and not optimized for anything (simple 
        ! straightforward implementation of methods). All implementation
        ! is basically done as one would do with a standard array. 
        ! This can be sufficient in some cases where the array changes
        ! only few times, in which case the overhead for other 
        ! other operations such as addition etc in optimized types may
        ! lead to poor performance. 

        real(kind=rk), allocatable   :: val(:) ! values 

    contains 

        ! Element insertion
        procedure   :: InsertSingleElement      => InsertSingleElementRDA
        procedure   :: InsertMultipleElements   => InsertMultipleElementsRDA
        generic     :: Insert                   => &
            InsertSingleElement, InsertMultipleElements

        ! Element appending
        procedure   :: AppendSingleElement      => AppendSingleElementRDA
        procedure   :: AppendMultipleElements   => AppendMultipleElementsRDA
        generic     :: Append                   => &
            AppendSingleElement, AppendMultipleElements

        ! Element removal
        procedure   :: RemoveSingleElement      => RemoveSingleElementRDA
        procedure   :: RemoveMultipleElements   => RemoveMultipleElementsRDA
        generic     :: Remove                   => &
            RemoveSingleElement, RemoveMultipleElements

        ! Element getter
        procedure   :: GetSingleElement      => GetSingleElementRDA
        procedure   :: GetMultipleElements   => GetMultipleElementsRDA
        procedure   :: GetAllElements        => GetAllElementsRDA
        generic     :: Get                   => &
            GetSingleElement, GetMultipleElements, GetAllElements

        ! Element setter
        procedure   :: SetSingleElement      => SetSingleElementRDA
        procedure   :: SetMultipleElements   => SetMultipleElementsRDA
        procedure   :: SetAllElements        => SetAllElementsRDA
        generic     :: Set                   => &
            SetSingleElement, SetMultipleElements, SetAllElements

    end type 

    ! Type optimized for appending values
    !type, extends(RealDynamicArrayUDT) :: RealAppendingOptimizedArrayUDT
    !contains
    !end type

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    ! Summation overriding
    interface operator(+)
        module procedure SumRDA, SumRDAscalar, SumRDAarray
    end interface

    ! Subtraction overriding
    interface operator(-)
        module procedure SubtractRDA, SubtractRDAscalar, SubtractRDAarray
    end interface

    ! Multiplication overriding
    interface operator(*)
        module procedure MultiplyRDA, MultiplyRDAscalar, MultiplyRDAarray
    end interface

    ! Division overriding
    interface operator(/)
        module procedure DivideRDA, DivideRDAscalar, DivideRDAarray
    end interface


contains

    
    !==================================================================!
    !                                                                  !
    !                             ROUTINES                             !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                        REAL DYNAMIC ARRAY                        !
    !------------------------------------------------------------------!

    ! Constructors
    !=============
    ! General constructor
    function ConstructRealDynamicArray(val) result(rda)

        ! Description
        !============
        ! Arrays are constructed empty by default. Adding values should
        ! be done later on (if the optional argument 'val' is given, 
        ! then these values are constructed directly)

        ! Declare variables
        !==================
        ! Arguments
        type(RealDynamicArrayUDT)           :: rda 
        real(R8), intent(in), optional      :: val(:)

        ! Construct
        !==========
        if (present(val)) then 
            rda%val = val 
        else 
            ! empty array
            allocate(rda%val(0))
        end if 
        
    end function

    ! Array manipulators
    !===================
    ! Single element appending 
    subroutine AppendSingleElementRDA(rda, val)

        ! Description
        !============
        ! Append a single element to the end of the array.

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayUDT)          :: rda 
        real(rk), intent(in)                :: val 

        ! Append
        !=======
        rda%val = [rda%val, val]

    end subroutine 

    ! Multiple element appending
    subroutine AppendMultipleElementsRDA(rda, val)

        ! Description
        !============
        ! Append a multiple elements to the end of the array.

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayUDT)          :: rda 
        real(rk), intent(in)                :: val(:) 

        ! Append
        !=======
        rda%val = [rda%val, val]

    end subroutine

    ! Single element insertion
    subroutine InsertSingleElementRDA(rda, val, loc)

        ! Description
        !============
        ! Add a single element to an arbitrary location to the array. 
        ! Note that we don't explicitly check for a feasible location
        ! value! In case it is negative, then the element will simply be
        ! added to the beginning of the array. If it's higher than the 
        ! array size, it will be appended. 

        ! Declare variables
        !==================
        ! Arguments 
        class(RealDynamicArrayUDT)      :: rda 
        real(rk), intent(in)            :: val 
        integer(ik), intent(in)         :: loc 

        ! Auxiliary
        real(rk), allocatable           :: temp1(:), temp2(:)

        ! Insert
        !=======
        ! Split up array
        temp1 = rda%val(1:loc-1)
        temp2 = rda%val(loc:size(rda%val))

        ! Reconstruct
        rda%val = [temp1, val, temp2]

    end subroutine 

    ! Multiple element insertion
    subroutine InsertMultipleElementsRDA(rda, val, loc)

        ! Description
        !============
        ! Add multiple elements to an arbitrary location to the array. 
        ! Here, we simply call the single insertion routine multiple
        ! times and update the location as needed. Probably not the 
        ! most efficient way for many elements or when called many times

        ! Declare variables
        !==================
        ! Arguments 
        class(RealDynamicArrayUDT)      :: rda 
        real(rk), intent(in)            :: val(:) 
        integer(ik), intent(in)         :: loc(:) 

        ! Auxiliary
        integer(ik)                     :: nval
        integer(ik)                     :: tloc(size(val))

        ! Loop
        integer(ik)                     :: i

        ! Check
        !======
        if (size(loc) /= size(val)) then 
            call gdErrorHandler('InsertMultipleElementsRDA: incompatible ' // &
                'sizes of val and loc')
        end if 

        ! Add
        !====
        tloc = loc
        nval = size(val)
        do i = 1, nval
            ! Add the next element
            call rda%Insert(val(i), tloc(i))

            ! Update the location
            where (tloc >= tloc(i)) tloc = tloc + 1
        end do 

    end subroutine 

    ! Single element deletion
    subroutine RemoveSingleElementRDA(rda, loc)

        ! Description
        !============
        ! Remove a single element at a certain array location. We do 
        ! not explicitly check if the location is valid. 

        ! Declare variables
        !==================
        ! Arguments 
        class(RealDynamicArrayUDT)      :: rda 
        integer(ik), intent(in)         :: loc 

        ! Auxiliary
        real(rk), allocatable           :: temp1(:), temp2(:)

        ! Remove
        !=======
        ! Split up array
        temp1 = rda%val(1:loc-1)
        temp2 = rda%val(loc+1:size(rda%val))

        ! Reconstruct
        rda%val = [temp1, temp2]

    end subroutine 

    ! Multiple element deletion
    subroutine RemoveMultipleElementsRDA(rda, loc)

        ! Description
        !============
        ! Remove multiple elements at a certain array location. We do 
        ! not explicitly check if the location is valid. 

        ! Declare variables
        !==================
        ! Arguments 
        class(RealDynamicArrayUDT)      :: rda 
        integer(ik), intent(in)         :: loc(:)

        ! Auxiliary
        logical, allocatable            :: mask(:)

        ! Remove
        !=======
        ! Set mask
        allocate(mask(size(rda%val)))
        mask = .true.
        mask(loc) = .false. 

        ! Remove
        rda%val = pack(rda%val, mask)

    end subroutine 

    ! Single element getter
    function GetSingleElementRDA(rda, loc) result(val)

        ! Description
        !============
        ! Get the value of a single element on a specified location. It
        ! is not checked whether this location is valid. 

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayUDT)  :: rda 
        integer(ik), intent(in)     :: loc 
        real(rk)                    :: val 

        ! Get 
        !====
        val = rda%val(loc)

    end function

    ! Multiple element getter
    function GetMultipleElementsRDA(rda, loc) result(val)

        ! Description
        !============
        ! Get the value of multiple elements on a specified location. It
        ! is not checked whether this location is valid. 

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayUDT)  :: rda 
        integer(ik), intent(in)     :: loc(:) 
        real(rk)                    :: val(size(loc))

        ! Get 
        !====
        val = rda%val(loc)
        
    end function

    ! All element getter
    function GetAllElementsRDA(rda) result(val)

        ! Description
        !============
        ! Get the value of all elements as an array. 

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayUDT)  :: rda 
        real(rk)                    :: val(size(rda%val))

        ! Get 
        !====
        val = rda%val
        
    end function

    ! Single element setter
    subroutine SetSingleElementRDA(rda, loc, val)

        ! Description
        !============
        ! Set the value of a single element on a specified location. It
        ! is not checked whether this location is valid. 

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayUDT)  :: rda 
        integer(ik), intent(in)     :: loc 
        real(rk)                    :: val 

        ! Set 
        !====
        rda%val(loc) = val

    end subroutine

    ! Multiple element setter
    subroutine SetMultipleElementsRDA(rda, loc, val)

        ! Description
        !============
        ! Set the value of multiple elements on a specified location. It
        ! is not checked whether this location is valid. 

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayUDT)  :: rda 
        integer(ik), intent(in)     :: loc(:) 
        real(rk)                    :: val(size(loc))

        ! Set 
        !====
        rda%val(loc) = val
        
    end subroutine

    ! All element setter
    subroutine SetAllElementsRDA(rda, val) 

        ! Description
        !============
        ! Set the value of all elements to a single scalar value

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayUDT)  :: rda 
        real(rk)                    :: val

        ! Set 
        !====
        rda%val = val
        
    end subroutine

    ! Elementray array operations
    !============================
    ! Summation
    function SumRDA(rda1, rda2) result(rda3)

        ! Description
        !============
        ! Sum two real dynamic arrays elementwise

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayUDT), intent(in)  :: rda1, rda2 
        class(RealDynamicArrayUDT), allocatable :: rda3 

        ! Sum values
        !===========
        rda3 = rda1 
        rda3%val = rda1%val + rda2%val

    end function

    function SumRDAscalar(rda1, scalar) result(rda3)

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayUDT), intent(in)  :: rda1 
        real(R8), intent(in)                    :: scalar 
        class(RealDynamicArrayUDT), allocatable :: rda3

        ! Sum
        !====
        rda3 = rda1 
        rda3%val = rda3%val + scalar

    end function 

    function SumRDAarray(rda1, array) result(rda3)

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayUDT), intent(in)  :: rda1 
        real(R8), intent(in)                    :: array(:)
        class(RealDynamicArrayUDT), allocatable :: rda3

        ! Sum
        !====
        rda3 = rda1 
        rda3%val = rda3%val + array

    end function

    ! Subtraction
    function SubtractRDA(rda1, rda2) result(rda3)

        ! Description
        !============
        ! Subtract two real dynamic arrays elementwise

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayUDT), intent(in)  :: rda1, rda2 
        class(RealDynamicArrayUDT), allocatable :: rda3 

        ! Subtract values
        !===========
        rda3 = rda1 
        rda3%val = rda1%val - rda2%val

    end function

    function SubtractRDAscalar(rda1, scalar) result(rda3)

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayUDT), intent(in)  :: rda1 
        real(R8), intent(in)                    :: scalar 
        class(RealDynamicArrayUDT), allocatable :: rda3

        ! Subtract
        !====
        rda3 = rda1 
        rda3%val = rda3%val - scalar

    end function 

    function SubtractRDAarray(rda1, array) result(rda3)

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayUDT), intent(in)  :: rda1 
        real(R8), intent(in)                    :: array(:)
        class(RealDynamicArrayUDT), allocatable :: rda3

        ! Subtract
        !====
        rda3 = rda1 
        rda3%val = rda3%val - array

    end function

    ! Multiplication
    function MultiplyRDA(rda1, rda2) result(rda3)

        ! Description
        !============
        ! Multiply two real dynamic arrays elementwise

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayUDT), intent(in)  :: rda1, rda2 
        class(RealDynamicArrayUDT), allocatable :: rda3 

        ! Multiply values
        !===========
        rda3 = rda1 
        rda3%val = rda1%val*rda2%val

    end function

    function MultiplyRDAscalar(rda1, scalar) result(rda3)

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayUDT), intent(in)  :: rda1 
        real(R8), intent(in)                    :: scalar 
        class(RealDynamicArrayUDT), allocatable :: rda3

        ! Multiply
        !====
        rda3 = rda1 
        rda3%val = rda3%val*scalar

    end function 

    function MultiplyRDAarray(rda1, array) result(rda3)

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayUDT), intent(in)  :: rda1 
        real(R8), intent(in)                    :: array(:)
        class(RealDynamicArrayUDT), allocatable :: rda3

        ! Multiply
        !====
        rda3 = rda1 
        rda3%val = rda3%val*array

    end function

    ! Division
    function DivideRDA(rda1, rda2) result(rda3)

        ! Description
        !============
        ! Divide two real dynamic arrays elementwise

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayUDT), intent(in)  :: rda1, rda2 
        class(RealDynamicArrayUDT), allocatable :: rda3 

        ! Divide values
        !===========
        rda3 = rda1 
        rda3%val = rda1%val/rda2%val

    end function

    function DivideRDAscalar(rda1, scalar) result(rda3)

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayUDT), intent(in)  :: rda1 
        real(R8), intent(in)                    :: scalar 
        class(RealDynamicArrayUDT), allocatable :: rda3

        ! Divide
        !====
        rda3 = rda1 
        rda3%val = rda3%val/scalar

    end function 

    function DivideRDAarray(rda1, array) result(rda3)

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayUDT), intent(in)  :: rda1 
        real(R8), intent(in)                    :: array(:)
        class(RealDynamicArrayUDT), allocatable :: rda3

        ! Divide
        !====
        rda3 = rda1 
        rda3%val = rda3%val/array

    end function
    

end module 