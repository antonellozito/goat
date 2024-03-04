!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains all 2D interpolant types and implementations. 
! Each interpolant inherits from the generic 2D interpolant type, which 
! defines which routines etc should be present to evaluate. Depends on
! the 'precision' module. 

module Interpolant2D

    ! Initialize
    !============
    ! Load modules
    use mod_precision

    !==================================================================!
    !                                                                  !
    !                            TYPES                                 !
    !                                                                  !
    !==================================================================!
    
    ! Generic interpolant type
    !=========================
    type, abstract :: GenericInterpolant2DUDT

        ! Description
        !============
        ! Generic 2D interpolant. Contains the following fields
        ! - xv, yv: coordinates at which the values are given
        ! - v:      the values at the coordinates
        ! Any other required fields are to be added by the specific
        ! interpolants themselves. 

        ! The following routines should be present in any 2D interpolant
        ! - Construct(xv, yv, v): interpolant constructor
        ! - Deconstruct()               : deconstructor
        ! - Evaluate(xq, yq, derivx, derivy) : evaluator, derivx and 
        !   derivy indicate the derivative order in x and y direction
        
        real(R8), allocatable   :: xv(:), yv(:), v(:)

    contains 

        ! Constructors (structured and unstructured)
        procedure :: ConstructStructuredGEN => ConstructFromStructuredData
        procedure :: ConstructUnstructuredGEN => ConstructFromUnstructuredData
        generic :: Construct => ConstructStructuredGEN, ConstructUnstructuredGEN

        procedure(ConstructSINT), deferred      :: ConstructStructured 
        procedure(ConstructUSINT), deferred     :: ConstructUnstructured  

        ! Evaluation
        procedure(EvaluateINT), deferred        :: Evaluate

    end type

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!
    ! Abstract interface
    abstract interface 
        
        ! Structured constructor
        subroutine ConstructSINT(interp, xg, yg, v)
            import :: GenericInterpolant2DUDT, R8
            class(GenericInterpolant2DUDT)      :: interp 
            real(R8), allocatable               :: xg(:), yg(:), v(:, :)   
        end subroutine  
        
        ! Unstructured constructor
        subroutine ConstructUSINT(interp, xg, yg, v)
            import :: GenericInterpolant2DUDT, R8
            class(GenericInterpolant2DUDT)      :: interp 
            real(R8), allocatable               :: xg(:), yg(:), v(:)   
        end subroutine  

        ! Evaluator
        subroutine EvaluateINT(interp, xq, yq, derivx, derivy, vq)
            import :: GenericInterpolant2DUDT, R8, I8 
            class(GenericInterpolant2DUDT)      :: interp 
            real(R8), intent(in)                :: xq(:), yq(:)
            real(R8), intent(out)               :: vq(:)
            integer(I8), intent(in)             :: derivx, derivy
        end subroutine  
        
    end interface        
    
    contains 

    ! Structured constructor routine (wrapper)
    subroutine ConstructFromStructuredData(interp, xg, yg, v) 

        class(GenericInterpolant2DUDT)      :: interp 
        real(R8), allocatable               :: xg(:), yg(:), v(:, :)

        call interp%ConstructStructured(xg, yg, v) 

    end subroutine 

    ! Unstructured constructor routine (wrapper)
    subroutine ConstructFromUnstructuredData(interp, xg, yg, v) 
        
        class(GenericInterpolant2DUDT)      :: interp 
        real(R8), allocatable               :: xg(:), yg(:), v(:)

        call interp%ConstructUnstructured(xg, yg, v)

    end subroutine 
    

end module