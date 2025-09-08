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
        procedure   :: SetAllElementsScalar  => SetAllElementsScalarRDA
        procedure   :: SetAllElementsArray   => SetAllElementsArrayRDA
        generic     :: Set                   => &
            SetSingleElement, SetMultipleElements, SetAllElementsScalar, &
            SetAllElementsArray
        
        ! Number of element getter
        procedure   :: Size                  => GetSizeRDA

        ! Sum with Mask
        procedure   :: SumMask               => SumRDAScalarMask

    end type 

    ! Integer dynamic array class
    type :: IntegerDynamicArrayUDT

        ! Description
        !============
        ! This is the base dynamic array type for general purpose 
        ! dynamic array operations. Only for Integer values (integer has 
        ! its own UDT) and not optimized for anything (simple 
        ! straightforward implementation of methods). All implementation
        ! is basically done as one would do with a standard array. 
        ! This can be sufficient in some cases where the array changes
        ! only few times, in which case the overhead for other 
        ! other operations such as addition etc in optimized types may
        ! lead to poor performance. 

        integer(kind=ik), allocatable   :: val(:) ! values 

    contains 

        ! Element insertion
        procedure   :: InsertSingleElement      => InsertSingleElementIDA
        procedure   :: InsertMultipleElements   => InsertMultipleElementsIDA
        generic     :: Insert                   => &
            InsertSingleElement, InsertMultipleElements

        ! Element appending
        procedure   :: AppendSingleElement      => AppendSingleElementIDA
        procedure   :: AppendMultipleElements   => AppendMultipleElementsIDA
        generic     :: Append                   => &
            AppendSingleElement, AppendMultipleElements

        ! Element removal
        procedure   :: RemoveSingleElement      => RemoveSingleElementIDA
        procedure   :: RemoveMultipleElements   => RemoveMultipleElementsIDA
        generic     :: Remove                   => &
            RemoveSingleElement, RemoveMultipleElements

        ! Element getter
        procedure   :: GetSingleElement      => GetSingleElementIDA
        procedure   :: GetMultipleElements   => GetMultipleElementsIDA
        procedure   :: GetAllElements        => GetAllElementsIDA
        generic     :: Get                   => &
            GetSingleElement, GetMultipleElements, GetAllElements

        ! Element setter
        procedure   :: SetSingleElement      => SetSingleElementIDA
        procedure   :: SetMultipleElements   => SetMultipleElementsIDA
        procedure   :: SetAllElementsScalar  => SetAllElementsScalarIDA
        procedure   :: SetAllElementsArray   => SetAllElementsArrayIDA
        generic     :: Set                   => &
            SetSingleElement, SetMultipleElements, SetAllElementsScalar, &
            SetAllElementsArray

        ! Number of element getter
        procedure   :: Size                  => GetSizeIDA

        ! Sum with Mask
        procedure   :: SumMask1D               => SumIDAScalarMask1D
        procedure   :: SumMask0D               => SumIDAScalarMask0D
        generic     :: SumMask                 => &
            SumMask0D, SumMask1D

        ! Replace
        procedure   :: Replace               => ReplaceIDA

        ! Update after entry removal
        procedure   :: UpdateArray

    end type
    
    ! Real dynamic array memory efficient class
    type, extends(RealDynamicArrayUDT) :: RealDynamicArrayBufferedUDT

        ! Description
        !============
        ! This is the base dynamic array type for general purpose 
        ! dynamic array operations. Only for Integer values (integer has 
        ! its own UDT). This type of array is typically larger than
        ! the amount of values to have less events where the 
        ! array size need to change when new elements are appended.

        integer(kind=ik)                :: nel    ! number of element in use

    contains

        ! Expanding array
        procedure   :: Expand                   => ExpandRDABuf 

        ! Element insertion
        procedure   :: InsertSingleElement      => InsertSingleElementRDABuf
        procedure   :: InsertMultipleElements   => InsertMultipleElementsRDABuf

        ! Element appending
        procedure   :: AppendSingleElement      => AppendSingleElementRDABuf
        procedure   :: AppendMultipleElements   => AppendMultipleElementsRDABuf

        ! Element removal
        procedure   :: RemoveSingleElement      => RemoveSingleElementRDABuf
        procedure   :: RemoveMultipleElements   => RemoveMultipleElementsRDABuf

        ! Element getter
        procedure   :: GetSingleElement      => GetSingleElementRDABuf
        procedure   :: GetMultipleElements   => GetMultipleElementsRDABuf
        procedure   :: GetAllElements        => GetAllElementsRDABuf

        ! Element setter
        procedure   :: SetSingleElement      => SetSingleElementRDABuf
        procedure   :: SetMultipleElements   => SetMultipleElementsRDABuf
        procedure   :: SetAllElementsScalar  => SetAllElementsScalarRDABuf
        procedure   :: SetAllElementsArray   => SetAllElementsArrayRDABuf
        
        ! Number of element getter
        procedure   :: Size                  => GetSizeRDABuf

        ! Sum with Mask
        procedure   :: SumMask1D               => SumRDAScalarMask1DBuf
        procedure   :: SumMask0D               => SumRDAScalarMask0DBuf

    end type

    ! Integer dynamic array memory efficient class
    type, extends(IntegerDynamicArrayUDT) :: IntegerDynamicArrayBufferedUDT

        ! Description
        !============
        ! This is the base dynamic array type for general purpose 
        ! dynamic array operations. Only for Integer values (integer has 
        ! its own UDT). This type of array is typically larger than
        ! the amount of values to have less events where the 
        ! array size need to change when new elements are appended.

        integer(kind=ik)                :: nel    ! number of element in use
        
    contains

        ! Expand array
        procedure   :: Expand                   => ExpandIDABuf

        ! Element insertion
        procedure   :: InsertSingleElement      => InsertSingleElementIDABuf
        procedure   :: InsertMultipleElements   => InsertMultipleElementsIDABuf

        ! Element appending
        procedure   :: AppendSingleElement      => AppendSingleElementIDABuf
        procedure   :: AppendMultipleElements   => AppendMultipleElementsIDABuf

        ! Element removal
        procedure   :: RemoveSingleElement      => RemoveSingleElementIDABuf
        procedure   :: RemoveMultipleElements   => RemoveMultipleElementsIDABuf

        ! Element getter
        procedure   :: GetSingleElement      => GetSingleElementIDABuf
        procedure   :: GetMultipleElements   => GetMultipleElementsIDABuf
        procedure   :: GetAllElements        => GetAllElementsIDABuf

        ! Element setter
        procedure   :: SetSingleElement      => SetSingleElementIDABuf
        procedure   :: SetMultipleElements   => SetMultipleElementsIDABuf
        procedure   :: SetAllElementsScalar  => SetAllElementsScalarIDABuf
        procedure   :: SetAllElementsArray   => SetAllElementsArrayIDABuf

        ! Number of element getter
        procedure   :: Size                  => GetSizeIDABuf

        ! Sum with Mask
        procedure   :: SumMask1D               => SumIDAScalarMask1DBuf
        procedure   :: SumMask0D               => SumIDAScalarMask0DBuf

        ! Replace
        procedure   :: Replace               => ReplaceIDABuf

        ! Update after entry removal
        procedure   :: UpdateArray           => UpdateArrayIDABuf

    end type
    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    ! Summation overriding
    interface operator(+)
        module procedure SumRDA, SumRDAscalar, SumRDAarray
        module procedure SumIDA, SumIDAscalar, SumIDAarray
    end interface

    ! Subtraction overriding
    interface operator(-)
        module procedure SubtractRDA, SubtractRDAscalar, SubtractRDAarray
        module procedure SubtractIDA, SubtractIDAscalar, SubtractIDAarray
    end interface

    ! Multiplication overriding
    interface operator(*)
        module procedure MultiplyRDA, MultiplyRDAscalar, MultiplyRDAarray
        module procedure MultiplyIDA, MultiplyIDAscalar, MultiplyIDAarray
    end interface

    ! Division overriding
    interface operator(/)
        module procedure DivideRDA, DivideRDAscalar, DivideRDAarray
        module procedure DivideIDA, DivideIDAscalar, DivideIDAarray
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
        ! Split 
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
        ! Set the value of a single element on a specified location. If
        ! the location is larger than the current array size, the array
        ! is extended with zeros up to the required array size. 

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayUDT)  :: rda 
        integer(ik), intent(in)     :: loc 
        real(rk)                    :: val 

        ! Set 
        !====
        if (rda%Size() < loc) then 
            call rda%Append(spread(0.0_R8, 1, loc - rda%Size()))
        end if 
        rda%val(loc) = val

    end subroutine

    ! Multiple element setter
    subroutine SetMultipleElementsRDA(rda, loc, val)

        ! Description
        !============
        ! Set the value of multiple elements on a specified location. If
        ! the location is larger than the current array size, the array
        ! is extended with zeros up to the required array size. 

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayUDT)  :: rda 
        integer(ik), intent(in)     :: loc(:) 
        real(rk)                    :: val(size(loc))

        ! Set 
        !====
        if (rda%Size() < maxval(loc)) then 
            call rda%Append(spread(0.0_R8, 1, maxval(loc) - rda%Size()))
        end if 
        rda%val(loc) = val
        
    end subroutine

    ! All element setter
    subroutine SetAllElementsScalarRDA(rda, val) 

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

    subroutine SetAllElementsArrayRDA(rda, val) 

        ! Description
        !============
        ! Set the value of all elements to the current array 

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayUDT)  :: rda 
        real(rk), dimension(:)      :: val

        ! Set 
        !====
        rda%val = val
        
    end subroutine

    ! Size getter
    function GetSizeRDA(rda) result(s)

        ! Declare variables
        !==================
        class(RealDynamicArrayUDT)  :: rda 
        integer(ik)                 :: s 

        ! Determine size
        !===============
        ! easy here
        s = size(rda%val)

    end function 

    ! Sum with Mask
    subroutine SumRDAScalarMask(rda,loc,val) 

        ! Description
        !============
        ! Increment the value of multiple elements on a specified location with val. If
        ! the location is larger than the current array size, the array
        ! is extended with zeros up to the required array size.
        
        ! Declare variables
        !==================
        class(RealDynamicArrayUDT)      :: rda 
        integer(ik), intent(in)         :: loc(:) 
        integer(rk)                     :: val

        if (rda%Size() < maxval(loc)) then 
            call rda%Append(spread(0.0_R8, 1, maxval(loc) - rda%Size()))
        end if 

        rda%val(loc) = rda%val(loc) + val
        
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

    !------------------------------------------------------------------!
    !                      INTEGER DYNAMIC ARRAY                       !
    !------------------------------------------------------------------!

    ! Constructors
    !=============
    ! General constructor
    function ConstructIntegerDynamicArray(val) result(ida)

        ! Description
        !============
        ! Arrays are constructed empty by default. Adding values should
        ! be done later on (if the optional argument 'val' is given, 
        ! then these values are constructed directly)

        ! Declare variables
        !==================
        ! Arguments
        type(IntegerDynamicArrayUDT)        :: ida 
        integer(I8), intent(in), optional   :: val(:)

        ! Construct
        !==========
        if (present(val)) then 
            ida%val = val 
        else 
            ! empty array
            allocate(ida%val(0))
        end if 
        
    end function

    ! Array manipulators
    !===================
    ! Single element appending 
    subroutine AppendSingleElementIDA(ida, val)

        ! Description
        !============
        ! Append a single element to the end of the array.

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayUDT)          :: ida 
        integer(ik), intent(in)                :: val 

        ! Append
        !=======
        ida%val = [ida%val, val]

    end subroutine 

    ! Multiple element appending
    subroutine AppendMultipleElementsIDA(ida, val)

        ! Description
        !============
        ! Append a multiple elements to the end of the array.

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayUDT)          :: ida 
        integer(ik), intent(in)                :: val(:) 

        ! Append
        !=======
        ida%val = [ida%val, val]

    end subroutine

    ! Single element insertion
    subroutine InsertSingleElementIDA(ida, val, loc)

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
        class(IntegerDynamicArrayUDT)       :: ida 
        integer(ik), intent(in)             :: val 
        integer(ik), intent(in)             :: loc 

        ! Auxiliary
        integer(ik), allocatable           :: temp1(:), temp2(:)

        ! Insert
        !=======
        ! Split up array
        temp1 = ida%val(1:loc-1)
        temp2 = ida%val(loc:size(ida%val))

        ! Reconstruct
        ida%val = [temp1, val, temp2]

    end subroutine 

    ! Multiple element insertion
    subroutine InsertMultipleElementsIDA(ida, val, loc)

        ! Description
        !============
        ! Add multiple elements to an arbitrary location to the array. 
        ! Here, we simply call the single insertion routine multiple
        ! times and update the location as needed. Probably not the 
        ! most efficient way for many elements or when called many times

        ! Declare variables
        !==================
        ! Arguments 
        class(IntegerDynamicArrayUDT)       :: ida 
        integer(ik), intent(in)             :: val(:) 
        integer(ik), intent(in)             :: loc(:) 

        ! Auxiliary
        integer(ik)                     :: nval
        integer(ik)                     :: tloc(size(val))

        ! Loop
        integer(ik)                     :: i

        ! Check
        !======
        if (size(loc) /= size(val)) then 
            call gdErrorHandler('InsertMultipleElementsIDA: incompatible ' // &
                'sizes of val and loc')
        end if 

        ! Add
        !====
        tloc = loc
        nval = size(val)
        do i = 1, nval
            ! Add the next element
            call ida%Insert(val(i), tloc(i))

            ! Update the location
            where (tloc >= tloc(i)) tloc = tloc + 1
        end do 

    end subroutine 

    ! Single element deletion
    subroutine RemoveSingleElementIDA(ida, loc)

        ! Description
        !============
        ! Remove a single element at a certain array location. We do 
        ! not explicitly check if the location is valid. 

        ! Declare variables
        !==================
        ! Arguments 
        class(IntegerDynamicArrayUDT)       :: ida 
        integer(ik), intent(in)             :: loc 

        ! Auxiliary
        integer(ik), allocatable           :: temp1(:), temp2(:)

        ! Remove
        !=======
        if (loc > size(ida%val)) then 
            call gdErrorHandler('RemoveSingleElementIDA: index is out of bounds')
        end if 
        ! Split up array
        temp1 = ida%val(1:loc-1)
        temp2 = ida%val(loc+1:size(ida%val))

        ! Reconstruct
        ida%val = [temp1, temp2]

    end subroutine 

    ! Multiple element deletion
    subroutine RemoveMultipleElementsIDA(ida, loc)

        ! Description
        !============
        ! Remove multiple elements at a certain array location. We do 
        ! not explicitly check if the location is valid. 

        ! Declare variables
        !==================
        ! Arguments 
        class(IntegerDynamicArrayUDT)       :: ida 
        integer(ik), intent(in)             :: loc(:)

        ! Auxiliary
        logical, allocatable            :: mask(:)

        ! Remove
        !=======
        ! Set mask
        if (any(loc > size(ida%val))) then 
            call gdErrorHandler('RemoveSingleElementIDA: index is out of bounds')
        end if 
        allocate(mask(size(ida%val)))
        mask = .true.
        mask(loc) = .false. 

        ! Remove
        ida%val = pack(ida%val, mask)

    end subroutine 

    ! Single element getter
    function GetSingleElementIDA(ida, loc) result(val)

        ! Description
        !============
        ! Get the value of a single element on a specified location. It
        ! is not checked whether this location is valid. 

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayUDT)   :: ida 
        integer(ik), intent(in)         :: loc 
        integer(ik)                     :: val 

        ! Get 
        !====
        val = ida%val(loc)

    end function

    ! Multiple element getter
    function GetMultipleElementsIDA(ida, loc) result(val)

        ! Description
        !============
        ! Get the value of multiple elements on a specified location. It
        ! is not checked whether this location is valid. 

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayUDT)   :: ida 
        integer(ik), intent(in)         :: loc(:) 
        integer(ik)                     :: val(size(loc))

        ! Get 
        !====
        val = ida%val(loc)
        
    end function

    ! All element getter
    function GetAllElementsIDA(ida) result(val)

        ! Description
        !============
        ! Get the value of all elements as an array. 

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayUDT)  :: ida 
        integer(ik)                    :: val(size(ida%val))

        ! Get 
        !====
        val = ida%val
        
    end function

    ! Single element setter
    subroutine SetSingleElementIDA(ida, loc, val)

        ! Description
        !============
        ! Set the value of a single element on a specified location. If
        ! the location is larger than the current array size, the array
        ! is extended with zeros up to the required array size. 

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayUDT)   :: ida 
        integer(ik), intent(in)         :: loc 
        integer(ik)                     :: val 

        ! Set 
        !====
        if (ida%Size() < loc) then 
            call ida%Append(spread(0_I8, 1, loc - ida%Size()))
        end if 
        ida%val(loc) = val

    end subroutine

    ! Multiple element setter
    subroutine SetMultipleElementsIDA(ida, loc, val)

        ! Description
        !============
        ! Set the value of multiple elements on a specified location. If
        ! the location is larger than the current array size, the array
        ! is extended with zeros up to the required array size. 

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayUDT)   :: ida 
        integer(ik), intent(in)         :: loc(:) 
        integer(ik)                     :: val(size(loc))

        ! Set 
        !====
        if (ida%Size() < maxval(loc)) then 
            call ida%Append(spread(0_I8, 1, maxval(loc) - ida%Size()))
        end if 
        ida%val(loc) = val
        
    end subroutine

    ! All element setter
    subroutine SetAllElementsScalarIDA(ida, val) 

        ! Description
        !============
        ! Set the value of all elements to a single scalar value

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayUDT)  :: ida 
        integer(ik)                    :: val

        ! Set 
        !====
        ida%val = val
        
    end subroutine

    subroutine SetAllElementsArrayIDA(ida, val) 

        ! Description
        !============
        ! Set the value of all elements to a single scalar value

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayUDT)  :: ida 
        integer(ik), dimension(:)      :: val

        ! Set 
        !====
        ida%val = val
        
    end subroutine
    
    ! Size getter
    function GetSizeIDA(ida) result(s)

        ! Declare variables
        !==================
        class(IntegerDynamicArrayUDT)   :: ida 
        integer(ik)                     :: s 

        ! Determine size
        !===============
        ! easy here
        s = size(ida%val)

    end function 

    ! Replace
    subroutine ReplaceIDA(ida,loc_begin,loc_end,val)
        
        ! Description
        !============
        ! Replaces the array elements located from index loc_begin
        ! till index loc_end, with the array val. This operation
        ! is a short cut for removing and inserting elements
        ! with different sizes at the same location.
        ! Note: it is not check whether loc_begin and loc_end are
        ! in the array size!

        ! Declare variables
        !==================
        class(IntegerDynamicArrayUDT)   :: ida
        integer(ik), intent(in)         :: loc_begin
        integer(ik), intent(in)         :: loc_end
        integer(ik), allocatable        :: val(:)

        ! Auxiliary
        integer(ik), allocatable           :: temp1(:), temp2(:)

        ! Insert val
        !===========
        ! Split up array
        temp1 = ida%val(1:loc_begin-1)
        temp2 = ida%val(loc_end+1:size(ida%val))

        ! Reconstruct
        ida%val = [temp1, val, temp2]        

    end subroutine

    ! Sum with Mask as single location
    subroutine SumIDAScalarMask0D(ida,loc,val) 

        ! Description
        !============
        ! Increment the value of one elements on a specified location with val.

        ! Declare variables
        !==================
        class(IntegerDynamicArrayUDT)   :: ida 
        integer(ik), intent(in)         :: loc 
        integer(ik)                     :: val

        !if (ida%Size() < loc) then 
        !    call ida%Append(spread(0_I8, 1, loc - ida%Size()))
        !end if 

        ida%val(loc) = ida%val(loc) + val
        
    end subroutine

    ! Sum with Mask as 1D array
    subroutine SumIDAScalarMask1D(ida,loc,val) 

        ! Description
        !============
        ! Increment the value of multiple elements on a specified location with val. If
        ! the location is larger than the current array size, the array
        ! is extended with zeros up to the required array size.

        ! Declare variables
        !==================
        class(IntegerDynamicArrayUDT)           :: ida 
        integer(ik), intent(in), allocatable    :: loc(:) 
        integer(ik)                             :: val

        if (ida%Size() < maxval(loc)) then 
            call ida%Append(spread(0_I8, 1, maxval(loc) - ida%Size()))
        end if 

        ida%val(loc) = ida%val(loc) + val
        
    end subroutine

    subroutine UpdateArray(ida,val)

        ! Description
        !============
        ! In the case where a entry of an array got removed
        ! and all the entry of the array which are larger than
        ! the removed array need to be decreased by one
        ! F.e. a vertex got removed out of vert%x etc. So 
        ! all vertices with larger number need to decrease by one
        ! in cell%vert and face%vert. This routine bundels the 
        ! needed operations

        ! Declare variables
        !==================
        class(IntegerDynamicArrayUDT) :: ida
        integer(I8) :: val

        ! Auxiliary
        integer(I8) :: j 
        integer(I8), allocatable :: ind(:), loc(:)

        ind = (/ (j, j = 1, size(ida%val)) /)
        loc = pack(ind, ida%val .gt. val)
        ida%val(loc) = ida%val(loc) - 1


    end subroutine

    ! Elementray array operations
    !============================
    ! Summation
    function SumIDA(ida1, ida2) result(ida3)

        ! Description
        !============
        ! Sum two real dynamic arrays elementwise

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayUDT), intent(in)  :: ida1, ida2 
        class(IntegerDynamicArrayUDT), allocatable :: ida3 

        ! Sum values
        !===========
        ida3 = ida1 
        ida3%val = ida1%val + ida2%val

    end function

    function SumIDAscalar(ida1, scalar) result(ida3)

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayUDT), intent(in)  :: ida1 
        integer(I8), intent(in)                    :: scalar 
        class(IntegerDynamicArrayUDT), allocatable :: ida3

        ! Sum
        !====
        ida3 = ida1 
        ida3%val = ida3%val + scalar

    end function 

    function SumIDAarray(ida1, array) result(ida3)

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayUDT), intent(in)  :: ida1 
        integer(I8), intent(in)                    :: array(:)
        class(IntegerDynamicArrayUDT), allocatable :: ida3

        ! Sum
        !====
        ida3 = ida1 
        ida3%val = ida3%val + array

    end function

    ! Subtraction
    function SubtractIDA(ida1, ida2) result(ida3)

        ! Description
        !============
        ! Subtract two real dynamic arrays elementwise

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayUDT), intent(in)  :: ida1, ida2 
        class(IntegerDynamicArrayUDT), allocatable :: ida3 

        ! Subtract values
        !===========
        ida3 = ida1 
        ida3%val = ida1%val - ida2%val

    end function

    function SubtractIDAscalar(ida1, scalar) result(ida3)

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayUDT), intent(in)  :: ida1 
        integer(I8), intent(in)                    :: scalar 
        class(IntegerDynamicArrayUDT), allocatable :: ida3

        ! Subtract
        !====
        ida3 = ida1 
        ida3%val = ida3%val - scalar

    end function 

    function SubtractIDAarray(ida1, array) result(ida3)

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayUDT), intent(in)  :: ida1 
        integer(I8), intent(in)                    :: array(:)
        class(IntegerDynamicArrayUDT), allocatable :: ida3

        ! Subtract
        !====
        ida3 = ida1 
        ida3%val = ida3%val - array

    end function

    ! Multiplication
    function MultiplyIDA(ida1, ida2) result(ida3)

        ! Description
        !============
        ! Multiply two real dynamic arrays elementwise

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayUDT), intent(in)  :: ida1, ida2 
        class(IntegerDynamicArrayUDT), allocatable :: ida3 

        ! Multiply values
        !===========
        ida3 = ida1 
        ida3%val = ida1%val*ida2%val

    end function

    function MultiplyIDAscalar(ida1, scalar) result(ida3)

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayUDT), intent(in)  :: ida1 
        integer(I8), intent(in)                    :: scalar 
        class(IntegerDynamicArrayUDT), allocatable :: ida3

        ! Multiply
        !====
        ida3 = ida1 
        ida3%val = ida3%val*scalar

    end function 

    function MultiplyIDAarray(ida1, array) result(ida3)

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayUDT), intent(in)  :: ida1 
        integer(I8), intent(in)                    :: array(:)
        class(IntegerDynamicArrayUDT), allocatable :: ida3

        ! Multiply
        !====
        ida3 = ida1 
        ida3%val = ida3%val*array

    end function

    ! Division
    function DivideIDA(ida1, ida2) result(ida3)

        ! Description
        !============
        ! Divide two real dynamic arrays elementwise - this is 
        ! integer division though!

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayUDT), intent(in)  :: ida1, ida2 
        class(IntegerDynamicArrayUDT), allocatable :: ida3 

        ! Divide values
        !===========
        ida3 = ida1 
        ida3%val = ida1%val/ida2%val

    end function

    function DivideIDAscalar(ida1, scalar) result(ida3)

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayUDT), intent(in)  :: ida1 
        integer(I8), intent(in)                    :: scalar 
        class(IntegerDynamicArrayUDT), allocatable :: ida3

        ! Divide
        !====
        ida3 = ida1 
        ida3%val = ida3%val/scalar

    end function 

    function DivideIDAarray(ida1, array) result(ida3)

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayUDT), intent(in)  :: ida1 
        integer(I8), intent(in)                    :: array(:)
        class(IntegerDynamicArrayUDT), allocatable :: ida3

        ! Divide
        !====
        ida3 = ida1 
        ida3%val = ida3%val/array

    end function

    !------------------------------------------------------------------!
    !                      REAL DYNAMIC ARRAY BUFFERED                 !
    !------------------------------------------------------------------!

    ! Constructors
    !=============
    ! General constructor
    function ConstructRealDynamicArrayBuffered(val) result(rda)

        ! Description
        !============
        ! Arrays are constructed empty by default. Adding values should
        ! be done later on (if the optional argument 'val' is given, 
        ! then these values are constructed directly)

        ! Declare variables
        !==================
        ! Arguments
        type(RealDynamicArrayBufferedUDT)      :: rda 
        real(R8), intent(in), optional   :: val(:)

        ! Auxiliary


        ! Construct
        !==========
        if (present(val)) then 
            rda%nel = size(val)
            allocate(rda%val(nint(rda%nel*1.1_R8)))
            rda%val(1:rda%nel) = val 
        else 
            ! empty array
            allocate(rda%val(0))
        end if 
        
    end function   

    ! Array manipulators
    !===================
    ! Expanding array
    subroutine ExpandRDABuf(rda, n)

        ! Description
        !============
        ! Expands an array with 10% extra size
        ! above the minimal needed expansion

        ! Declare variables
        !==================
        class(RealDynamicArrayBufferedUDT) :: rda
        integer(ik), intent(in)         :: n

        ! Auxiliary
        real(rk), allocatable        :: val1(:)

        ! Determine new size
        val1 = rda%val
        deallocate(rda%val)
        allocate(rda%val(nint(size(val1)*1.1_R8) + n))
        rda%val = 0
        rda%val(1:size(val1)) = val1

    end subroutine

    ! Single element appending 
    subroutine AppendSingleElementRDABuf(rda, val)

        ! Description
        !============
        ! Append a single element to the end of the array.

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayBufferedUDT)   :: rda 
        real(rk), intent(in)                    :: val 

        ! Append
        !=======
        ! Update counter
        rda%nel = rda%nel + 1
            
        if (rda%nel .gt. size(rda%val)) then
            ! Make the array larger
            call rda%Expand(1)
        end if

        rda%val(rda%nel) = val

    end subroutine 

    ! Multiple element appending
    subroutine AppendMultipleElementsRDABuf(rda, val)

        ! Description
        !============
        ! Append a multiple elements to the end of the array.

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayBufferedUDT)      :: rda 
        real(rk), intent(in)                    :: val(:) 

        ! Auxiliary
        integer(ik)                             :: old_size

        ! Append
        !=======
        ! Update counter
        old_size = rda%nel
        rda%nel = rda%nel + size(val)
            
        if (rda%nel .gt. size(rda%val)) then

            ! Make the array larger
            call rda%Expand(size(val))

        end if

        rda%val(old_size+1:rda%nel) = val

    end subroutine   
    
    ! Single element insertion
    subroutine InsertSingleElementRDABuf(rda, val, loc)

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
        class(RealDynamicArrayBufferedUDT)      :: rda 
        real(rk), intent(in)                    :: val 
        integer(ik), intent(in)                 :: loc 

        ! Auxiliary
        real(rk), allocatable           :: temp1(:)

        ! Insert
        !=======
        ! Split up array
        temp1 = rda%val(loc:rda%nel)

        ! Update counter
        rda%nel = rda%nel + 1

        if (rda%nel .gt. size(rda%val)) then
            call rda%Expand(1)
        end if

        ! Reconstruct
        rda%val(loc:rda%nel) = [val, temp1]

    end subroutine

    ! Multiple element insertion
    subroutine InsertMultipleElementsRDABuf(rda, val, loc)

        ! Description
        !============
        ! Add multiple elements to an arbitrary location to the array. 
        ! Here, we simply call the single insertion routine multiple
        ! times and update the location as needed. Probably not the 
        ! most efficient way for many elements or when called many times

        ! Declare variables
        !==================
        ! Arguments 
        class(RealDynamicArrayBufferedUDT)      :: rda 
        real(rk), intent(in)                    :: val(:) 
        integer(ik), intent(in)                 :: loc(:) 

        ! Auxiliary
        integer(ik)                     :: nval
        integer(ik)                     :: tloc(size(val))

        ! Loop
        integer(ik)                     :: i

        ! Check
        !======
        if (size(loc) /= size(val)) then 
            call gdErrorHandler('InsertMultipleElementsRDABuf: incompatible ' // &
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
    subroutine RemoveSingleElementRDABuf(rda, loc)

        ! Description
        !============
        ! Remove a single element at a certain array location. We do 
        ! not explicitly check if the location is valid. 

        ! Declare variables
        !==================
        ! Arguments 
        class(RealDynamicArrayBufferedUDT)      :: rda 
        integer(ik), intent(in)                 :: loc 

        ! Auxiliary
        real(rk), allocatable                   :: temp1(:)

        ! Remove
        !=======
        if (loc > size(rda%val)) then 
            call gdErrorHandler('RemoveSingleElementGArda: index is out of bounds')
        end if 

        ! Split up array - take only slice that moves
        temp1 = rda%val(loc+1:rda%nel)

        ! Reconstruct
        rda%nel = rda%nel - 1
        rda%val(loc:rda%nel) = temp1
        rda%val(rda%nel + 1:size(rda%val) ) = 0

    end subroutine

    ! Multiple element deletion
    subroutine RemoveMultipleElementsRDABuf(rda, loc)

        ! Description
        !============
        ! Remove multiple elements at a certain array location. We do 
        ! not explicitly check if the location is valid. 

        ! Declare variables
        !==================
        ! Arguments 
        class(RealDynamicArrayBufferedUDT)      :: rda 
        integer(ik), intent(in)                 :: loc(:)

        ! Auxiliary
        logical, allocatable                :: mask(:)
        integer(ik)                         :: old_size

        ! Remove
        !=======
        ! Set mask
        old_size = rda%nel
        if (any(loc > size(rda%val))) then 
            call gdErrorHandler('RemoveSingleElementRDABuf: index is out of bounds')
        end if 
        allocate(mask(old_size))
        mask = .true.
        mask(loc) = .false. 

        ! Remove
        rda%nel = rda%nel - size(loc)
        rda%val(1:rda%nel) = pack(rda%val(1:old_size), mask)
        rda%val(rda%nel+1:size(rda%val)) = 0

        ! Maybe shrink the array if size is big difference - TODO


    end subroutine 

    ! Single element getter
    function GetSingleElementRDABuf(rda, loc) result(val)

        ! Description
        !============
        ! Get the value of a single element on a specified location. It
        ! is checked whether this location is valid. 

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayBufferedUDT)  :: rda 
        integer(ik), intent(in)             :: loc 
        real(rk)                            :: val 

        ! Get 
        !====
        val = rda%val(loc)

        if (loc .gt. rda%nel) then
            call gdErrorHandler('GetSingleElementRDABuf: invald location')
        end if

    end function

    ! Multiple element getter
    function GetMultipleElementsRDABuf(rda, loc) result(val)

        ! Description
        !============
        ! Get the value of multiple elements on a specified location. It
        ! is not checked whether this location is valid. 

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayBufferedUDT)  :: rda 
        integer(ik), intent(in)             :: loc(:) 
        real(rk)                            :: val(size(loc))

        ! Get 
        !====
        val = rda%val(loc)

        if (maxval(loc) .gt. rda%nel) then
            call gdErrorHandler('GetMultipleElementsRDABuf: invalid location')
        end if
        
    end function

    ! All element getter
    function GetAllElementsRDABuf(rda) result(val)

        ! Description
        !============
        ! Get the value of all elements as an array. 

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayBufferedUDT)  :: rda 
        real(rk)                            :: val(rda%nel)

        ! Get 
        !====
        val = rda%val(1:rda%nel)
        
    end function

    ! Single element setter
    subroutine SetSingleElementRDABuf(rda, loc, val)

        ! Description
        !============
        ! Set the value of a single element on a specified location. If
        ! the location is larger than the current array size, the array
        ! is extended with zeros up to the required array size and
        ! 10 percent longer.

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayBufferedUDT)      :: rda 
        integer(ik), intent(in)                 :: loc 
        real(rk)                                :: val 

        ! Set 
        !====
        if (size(rda%val) < loc) then 
            call rda%Expand(loc - size(rda%val))
            rda%nel = loc
        else if (loc > rda%nel) then
            rda%nel = loc
        end if
        rda%val(loc) = val

    end subroutine
    ! Multiple element setter
    subroutine SetMultipleElementsRDABuf(rda, loc, val)

        ! Description
        !============
        ! Set the value of multiple elements on a specified location. If
        ! the location is larger than the current array size, the array
        ! is extended with zeros up to the required array size and
        ! 10 procent more. 

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayBufferedUDT)  :: rda 
        integer(ik), intent(in)             :: loc(:) 
        real(rk)                            :: val(size(loc))

        ! Set 
        !====
        if (size(rda%val) < maxval(loc)) then 
            call rda%Expand(maxval(loc) - size(rda%val))
            rda%nel = maxval(loc)
        else if (maxval(loc) .gt. rda%nel) then
            rda%nel = maxval(loc)
        end if
        rda%val(loc) = val
        
    end subroutine

    ! All element setter
    subroutine SetAllElementsScalarRDABuf(rda, val) 

        ! Description
        !============
        ! Set the value of all elements to a single scalar value

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayBufferedUDT)  :: rda 
        real(rk)                            :: val

        ! Set 
        !====
        rda%val(1:rda%nel) = val
        
    end subroutine

    ! All element setter
    subroutine SetAllElementsArrayRDABuf(rda, val) 

        ! Description
        !============
        ! Set the value of all elements to a single scalar value

        ! Declare variables
        !==================
        ! Arguments
        class(RealDynamicArrayBufferedUDT)  :: rda 
        real(rk), dimension(:)              :: val

        ! Set 
        !====
        rda%val(1:rda%nel) = val
        
    end subroutine
    
    ! Size getter
    function GetSizeRDABuf(rda) result(s)

        ! Declare variables
        !==================
        class(RealDynamicArrayBufferedUDT)   :: rda 
        integer(ik)                  :: s 

        ! Determine size
        !===============
        ! easy here
        s = rda%nel

    end function 

    ! Replace
    subroutine ReplaceRDABuf(rda,loc_begin,loc_end,val)
        
        ! Description
        !============
        ! Replaces the array elements located from index loc_begin
        ! till index loc_end, with the array val. This operation
        ! is a short cut for removing and inserting elements
        ! with different sizes at the same location.
        ! Note: if the size loc_begin to loc_end is longer than
        ! the size of val, the rest of the element will be zeros

        ! Declare variables
        !==================
        class(RealDynamicArrayBufferedUDT)      :: rda
        integer(ik), intent(in)                 :: loc_begin
        integer(ik), intent(in)                 :: loc_end
        real(rk), allocatable                   :: val(:)

        ! Auxiliary
        integer(ik)                         :: size_place

        ! Insert val
        !===========
        ! Check size of locations
        size_place = loc_end - loc_begin + 1
        if (size_place .lt. size(val)) then

            call gdErrorHandler('ReplaceRDABuf: not able to replace')

        end if 

        ! Check whether appending is needed
        if (loc_end .gt. size(rda%val)) then
            call rda%Expand(loc_end - size(rda%val))
            rda%nel = loc_end
        else if (loc_end .gt. rda%nel) then
            rda%nel = loc_end
        end if 

        ! Assign
        rda%val(loc_begin:loc_end) = 0
        rda%val(loc_begin:loc_begin+size(val)-1) = val        

    end subroutine

    ! Sum with Mask as single location
    subroutine SumRDAScalarMask0DBuf(rda,loc,val) 

        ! Description
        !============
        ! Increment the value of one elements on a specified location with val.

        ! Declare variables
        !==================
        class(RealDynamicArrayBufferedUDT)      :: rda 
        integer(ik), intent(in)                 :: loc 
        real(rk)                                :: val


        if (size(rda%val) < loc) then 
            call rda%Append(spread(0.0_R8, 1, loc - size(rda%val)))
        end if 
        
        ! Sum
        rda%val(loc) = rda%val(loc) + val
        
    end subroutine  
    
    ! Sum with Mask as 1D array
    subroutine SumRDAScalarMask1DBuf(rda,loc,val) 

        ! Description
        !============
        ! Increment the value of multiple elements on a specified location with val. If
        ! the location is larger than the current array size, the array
        ! is extended with zeros up to the required array size.

        ! Declare variables
        !==================
        class(RealDynamicArrayBufferedUDT)      :: rda 
        integer(ik), intent(in), allocatable    :: loc(:) 
        real(rk)                                :: val

        if (size(rda%val) < maxval(loc)) then 
            call rda%Append(spread(0.0_R8, 1, maxval(loc) - size(rda%val)))
        end if 

        rda%val(loc) = rda%val(loc) + val
        
    end subroutine 
    
    subroutine UpdateArrayRDABuf(rda,val)

        ! Description
        !============
        ! In the case where a entry of an array got removed
        ! and all the entry of the array which are larger than
        ! the removed array need to be decreased by one
        ! F.e. a vertex got removed out of vert%x etc. So 
        ! all vertices with larger number need to decrease by one
        ! in cell%vert and face%vert. This routine bundels the 
        ! needed operations

        ! Declare variables
        !==================
        class(RealDynamicArrayBufferedUDT) :: rda
        real(rk) :: val

        ! Auxiliary
        integer(I8) :: j 
        integer(I8), allocatable :: ind(:), loc(:)

        ind = (/ (j, j = 1, rda%nel) /)
        loc = pack(ind, rda%val(1:rda%nel) .gt. val)
        rda%val(loc) = rda%val(loc) - 1

    end subroutine  

    !------------------------------------------------------------------!
    !                  INTEGER DYNAMIC ARRAY BUFFERED                  !
    !------------------------------------------------------------------!

    ! Constructors
    !=============
    ! General constructor
    function ConstructIntegerDynamicArrayBuffered(val) result(ida)

        ! Description
        !============
        ! Arrays are constructed empty by default. Adding values should
        ! be done later on (if the optional argument 'val' is given, 
        ! then these values are constructed directly)

        ! Declare variables
        !==================
        ! Arguments
        type(IntegerDynamicArrayBufferedUDT)      :: ida 
        integer(I8), intent(in), optional   :: val(:)

        ! Auxiliary


        ! Construct
        !==========
        if (present(val)) then 
            ida%nel = size(val)
            allocate(ida%val(nint(ida%nel*1.1_R8)))
            ida%val(1:ida%nel) = val 
        else 
            ! empty array
            allocate(ida%val(0))
        end if 
        
    end function   
    
    ! Array manipulators
    !===================
    ! Expanding array
    subroutine ExpandIDABuf(ida, n)

        ! Description
        !============
        ! Expands an array with 10% extra size
        ! above the minimal needed expansion

        ! Declare variables
        !==================
        class(IntegerDynamicArrayBufferedUDT) :: ida
        integer(ik), intent(in)         :: n

        ! Auxiliary
        integer(ik), allocatable        :: val1(:)

        ! Determine new size
        val1 = ida%val
        deallocate(ida%val)
        allocate(ida%val(nint(size(val1)*1.1_R8) + n))
        ida%val = 0
        ida%val(1:size(val1)) = val1

    end subroutine

    ! Single element appending 
    subroutine AppendSingleElementIDABuf(ida, val)

        ! Description
        !============
        ! Append a single element to the end of the array.

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayBufferedUDT)          :: ida 
        integer(ik), intent(in)                :: val 

        ! Append
        !=======
        ! Update counter
        ida%nel = ida%nel + 1
            
        if (ida%nel .gt. size(ida%val)) then
            ! Make the array larger
            call ida%Expand(1)
        end if

        ida%val(ida%nel) = val

    end subroutine 

    ! Multiple element appending
    subroutine AppendMultipleElementsIDABuf(ida, val)

        ! Description
        !============
        ! Append a multiple elements to the end of the array.

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayBufferedUDT)          :: ida 
        integer(ik), intent(in)                :: val(:) 

        ! Auxiliary
        integer(ik)                             :: old_size

        ! Append
        !=======
        ! Update counter
        old_size = ida%nel
        ida%nel = ida%nel + size(val)
            
        if (ida%nel .gt. size(ida%val)) then

            ! Make the array larger
            call ida%Expand(size(val))

        end if

        ida%val(old_size+1:ida%nel) = val

    end subroutine   
    
    ! Single element insertion
    subroutine InsertSingleElementIDABuf(ida, val, loc)

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
        class(IntegerDynamicArrayBufferedUDT)       :: ida 
        integer(ik), intent(in)             :: val 
        integer(ik), intent(in)             :: loc 

        ! Auxiliary
        integer(ik), allocatable           :: temp1(:)

        ! Insert
        !=======
        ! Split up array
        temp1 = ida%val(loc:ida%nel)

        ! Update counter
        ida%nel = ida%nel + 1

        if (ida%nel .gt. size(ida%val)) then
            call ida%Expand(1)
        end if

        ! Reconstruct
        ida%val(loc:ida%nel) = [val, temp1]

    end subroutine

    ! Multiple element insertion
    subroutine InsertMultipleElementsIDABuf(ida, val, loc)

        ! Description
        !============
        ! Add multiple elements to an arbitrary location to the array. 
        ! Here, we simply call the single insertion routine multiple
        ! times and update the location as needed. Probably not the 
        ! most efficient way for many elements or when called many times

        ! Declare variables
        !==================
        ! Arguments 
        class(IntegerDynamicArrayBufferedUDT)       :: ida 
        integer(ik), intent(in)             :: val(:) 
        integer(ik), intent(in)             :: loc(:) 

        ! Auxiliary
        integer(ik)                     :: nval
        integer(ik)                     :: tloc(size(val))

        ! Loop
        integer(ik)                     :: i

        ! Check
        !======
        if (size(loc) /= size(val)) then 
            call gdErrorHandler('InsertMultipleElementsIDABuf: incompatible ' // &
                'sizes of val and loc')
        end if 

        ! Add
        !====
        tloc = loc
        nval = size(val)
        do i = 1, nval
            ! Add the next element
            call ida%Insert(val(i), tloc(i))

            ! Update the location
            where (tloc >= tloc(i)) tloc = tloc + 1
        end do 

    end subroutine 

    ! Single element deletion
    subroutine RemoveSingleElementIDABuf(ida, loc)

        ! Description
        !============
        ! Remove a single element at a certain array location. We do 
        ! not explicitly check if the location is valid. 

        ! Declare variables
        !==================
        ! Arguments 
        class(IntegerDynamicArrayBufferedUDT)     :: ida 
        integer(ik), intent(in)             :: loc 

        ! Auxiliary
        integer(ik), allocatable            :: temp1(:)

        ! Remove
        !=======
        if (loc > size(ida%val)) then 
            call gdErrorHandler('RemoveSingleElementGAIDA: index is out of bounds')
        end if 

        ! Split up array - take only slice that moves
        temp1 = ida%val(loc+1:ida%nel)

        ! Reconstruct
        ida%nel = ida%nel - 1
        ida%val(loc:ida%nel) = temp1
        ida%val(ida%nel + 1:size(ida%val) ) = 0

    end subroutine

    ! Multiple element deletion
    subroutine RemoveMultipleElementsIDABuf(ida, loc)

        ! Description
        !============
        ! Remove multiple elements at a certain array location. We do 
        ! not explicitly check if the location is valid. 

        ! Declare variables
        !==================
        ! Arguments 
        class(IntegerDynamicArrayBufferedUDT)     :: ida 
        integer(ik), intent(in)             :: loc(:)

        ! Auxiliary
        logical, allocatable                :: mask(:)
        integer(ik)                         :: old_size

        ! Remove
        !=======
        ! Set mask
        old_size = ida%nel
        if (any(loc > size(ida%val))) then 
            call gdErrorHandler('RemoveSingleElementIDABuf: index is out of bounds')
        end if 
        allocate(mask(old_size))
        mask = .true.
        mask(loc) = .false. 

        ! Remove
        ida%nel = ida%nel - size(loc)
        ida%val(1:ida%nel) = pack(ida%val(1:old_size), mask)
        ida%val(ida%nel+1:size(ida%val)) = 0

        ! Maybe shrink the array if size is big difference - TODO


    end subroutine 

    ! Single element getter
    function GetSingleElementIDABuf(ida, loc) result(val)

        ! Description
        !============
        ! Get the value of a single element on a specified location. It
        ! is checked whether this location is valid. 

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayBufferedUDT) :: ida 
        integer(ik), intent(in)         :: loc 
        integer(ik)                     :: val 

        ! Get 
        !====
        val = ida%val(loc)

        if (loc .gt. ida%nel) then
            call gdErrorHandler('GetSingleElementIDABuf: invald location')
        end if

    end function

    ! Multiple element getter
    function GetMultipleElementsIDABuf(ida, loc) result(val)

        ! Description
        !============
        ! Get the value of multiple elements on a specified location. It
        ! is not checked whether this location is valid. 

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayBufferedUDT)   :: ida 
        integer(ik), intent(in)         :: loc(:) 
        integer(ik)                     :: val(size(loc))

        ! Get 
        !====
        val = ida%val(loc)

        if (maxval(loc) .gt. ida%nel) then
            call gdErrorHandler('GetMultipleElementsIDABuf: invalid location')
        end if
        
    end function

    ! All element getter
    function GetAllElementsIDABuf(ida) result(val)

        ! Description
        !============
        ! Get the value of all elements as an array. 

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayBufferedUDT)     :: ida 
        integer(ik)                         :: val(ida%nel)

        ! Get 
        !====
        val = ida%val(1:ida%nel)
        
    end function

    ! Single element setter
    subroutine SetSingleElementIDABuf(ida, loc, val)

        ! Description
        !============
        ! Set the value of a single element on a specified location. If
        ! the location is larger than the current array size, the array
        ! is extended with zeros up to the required array size and
        ! 10 percent longer.

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayBufferedUDT)   :: ida 
        integer(ik), intent(in)         :: loc 
        integer(ik)                     :: val 

        ! Set 
        !====
        if (size(ida%val) < loc) then 
            call ida%Expand(loc - size(ida%val))
            ida%nel = loc
        else if (loc > ida%nel) then
            ida%nel = loc
        end if
        ida%val(loc) = val

    end subroutine
    ! Multiple element setter
    subroutine SetMultipleElementsIDABuf(ida, loc, val)

        ! Description
        !============
        ! Set the value of multiple elements on a specified location. If
        ! the location is larger than the current array size, the array
        ! is extended with zeros up to the required array size and
        ! 10 procent more. 

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayBufferedUDT) :: ida 
        integer(ik), intent(in)         :: loc(:) 
        integer(ik)                     :: val(size(loc))

        ! Set 
        !====
        if (size(ida%val) < maxval(loc)) then 
            call ida%Expand(maxval(loc) - size(ida%val))
            ida%nel = maxval(loc)
        else if (maxval(loc) .gt. ida%nel) then
            ida%nel = maxval(loc)
        end if
        ida%val(loc) = val
        
    end subroutine

    ! All element setter
    subroutine SetAllElementsScalarIDABuf(ida, val) 

        ! Description
        !============
        ! Set the value of all elements to a single scalar value

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayBufferedUDT)  :: ida 
        integer(ik)                    :: val

        ! Set 
        !====
        ida%val(1:ida%nel) = val
        
    end subroutine

    ! All element setter
    subroutine SetAllElementsArrayIDABuf(ida, val) 

        ! Description
        !============
        ! Set the value of all elements to a single scalar value

        ! Declare variables
        !==================
        ! Arguments
        class(IntegerDynamicArrayBufferedUDT)  :: ida 
        integer(ik), dimension(:)      :: val

        ! Set 
        !====
        ida%val(1:ida%nel) = val
        
    end subroutine
    
    ! Size getter
    function GetSizeIDABuf(ida) result(s)

        ! Declare variables
        !==================
        class(IntegerDynamicArrayBufferedUDT)   :: ida 
        integer(ik)                     :: s 

        ! Determine size
        !===============
        ! easy here
        s = ida%nel

    end function 

    ! Replace
    subroutine ReplaceIDABuf(ida,loc_begin,loc_end,val)
        
        ! Description
        !============
        ! Replaces the array elements located from index loc_begin
        ! till index loc_end, with the array val. This operation
        ! is a short cut for removing and inserting elements
        ! with different sizes at the same location.
        ! Note: if the size loc_begin to loc_end is longer than
        ! the size of val, the rest of the element will be zeros

        ! Declare variables
        !==================
        class(IntegerDynamicArrayBufferedUDT)     :: ida
        integer(ik), intent(in)             :: loc_begin
        integer(ik), intent(in)             :: loc_end
        integer(ik), allocatable            :: val(:)

        ! Auxiliary
        integer(ik)                         :: size_place

        ! Insert val
        !===========
        ! Check size of locations
        size_place = loc_end - loc_begin + 1
        if (size_place .lt. size(val)) then

            call gdErrorHandler('ReplaceIDABuf: not able to replace')

        end if 

        ! Check whether appending is needed
        if (loc_end .gt. size(ida%val)) then
            call ida%Expand(loc_end - size(ida%val))
            ida%nel = loc_end
        else if (loc_end .gt. ida%nel) then
            ida%nel = loc_end
        end if 

        ! Assign
        ida%val(loc_begin:loc_end) = 0
        ida%val(loc_begin:loc_begin+size(val)-1) = val        

    end subroutine

    ! Sum with Mask as single location
    subroutine SumIDAScalarMask0DBuf(ida,loc,val) 

        ! Description
        !============
        ! Increment the value of one elements on a specified location with val.

        ! Declare variables
        !==================
        class(IntegerDynamicArrayBufferedUDT)   :: ida 
        integer(ik), intent(in)         :: loc 
        integer(ik)                     :: val

        ! Sum
        ida%val(loc) = ida%val(loc) + val
        
    end subroutine  
    
    ! Sum with Mask as 1D array
    subroutine SumIDAScalarMask1DBuf(ida,loc,val) 

        ! Description
        !============
        ! Increment the value of multiple elements on a specified location with val. If
        ! the location is larger than the current array size, the array
        ! is extended with zeros up to the required array size.

        ! Declare variables
        !==================
        class(IntegerDynamicArrayBufferedUDT)           :: ida 
        integer(ik), intent(in), allocatable    :: loc(:) 
        integer(ik)                             :: val

        if (size(ida%val) < maxval(loc)) then 
            call ida%Append(spread(0_I8, 1, maxval(loc) - size(ida%val)))
        end if 

        ida%val(loc) = ida%val(loc) + val
        
    end subroutine 
    
    subroutine UpdateArrayIDABuf(ida,val)

        ! Description
        !============
        ! In the case where a entry of an array got removed
        ! and all the entry of the array which are larger than
        ! the removed array need to be decreased by one
        ! F.e. a vertex got removed out of vert%x etc. So 
        ! all vertices with larger number need to decrease by one
        ! in cell%vert and face%vert. This routine bundels the 
        ! needed operations

        ! Declare variables
        !==================
        class(IntegerDynamicArrayBufferedUDT) :: ida
        integer(I8) :: val

        ! Auxiliary
        integer(I8) :: j 
        integer(I8), allocatable :: ind(:), loc(:)

        ind = (/ (j, j = 1, ida%nel) /)
        loc = pack(ind, ida%val(1:ida%nel) .gt. val)
        ida%val(loc) = ida%val(loc) - 1

    end subroutine  

end module 