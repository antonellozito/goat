!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains general diagnostics routines and derived types, 
! such as for example routines to check gradients and hessians with 
! finite differences. This module should be used as a sort of 'template'
! for one's own implementation, in the sense that the objects and 
! routines here provide the main machinery for the diagnostics, but that
! specific implementation (e.g. how to evaluate a certain function)
! is left to the user. An example:
!
! One wants to use a finite difference checker to see if the gradient
! is correctly computed. The user-specific field here is the 'fun' field
! where the function evaluation etc should be defined. The user then 
! has to create a module where he defines a derived type that inherits
! from the basic 'fun' (abstract) type, and where he implements the 
! required methods that adhere to the interface described here. Then, 
! the user simply writes a program, where the finite difference checker
! is used and assigns the 'fun' field to be of his derived type. 

! This allows modular application of this functionality at the penalty
! of implementing one's own evaluation routines. 

! Notes
!======
! Note 1: currently, only a finite differencing checker is implemented. 
! It can cope with gradient and hessian estimates. See also the 
! documentation of these routines.

! Note 2: this module relies on the mod_precision and 
! mod_sparseinterface modules for some derived types and precision 
! specifications. 

module mod_diagnostics
    
    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_sparseinterface

    ! The usual
    implicit none
    save
    private
    
    public DiagnosticsFunctionUDT, FDcheckerUDT

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    ! Abstract types
    !===============
    ! Function handle for finite difference checker
    type, abstract :: DiagnosticsFunctionUDT

        ! Description
        !============
        ! This abstract type defines the requirements for a diagnostics
        ! function that is used by the FDchecker. No fields are 
        ! explicitly required - it is assumed that the necessary data
        ! for function evaluation is fully stored in the derived type
        ! which is to be constructed by the user. 

    contains

        ! Function evaluation
        procedure(EvaluateDiagnosticsFunctionINT), deferred :: &
            Evaluate 

        ! Dimension getter
        procedure(GetDimensionsDiagnosticsFunctionINT), deferred :: &
            GetDimensions

        ! Argument getter
        procedure(GetArgumuentsDiagnosticsFunctionINT), deferred :: & 
            GetArguments 

    end type

    ! Derived types
    !==============
    ! Finite difference checker
    type :: FDcheckerUDT

        ! Description
        !============
        ! Object that is used to see if the obtained gradient and/or 
        ! hessian of a function is accurate, compared to the finite 
        ! difference one. The evaluation etc is encapsulated in the 
        ! 'fun' field, which contains the function to be evaluated. 
        ! Note: the 'fun' field must be allocated by the user! 

        ! Fields:
        ! - fun:        function to be evaluated and checked
        ! - nvars:      number of variables to check
        ! - vars:       variable indices that have to be checked 
        ! - d:          step size vector
        ! - nd:         number of steps

        class(DiagnosticsFunctionUDT), allocatable      :: fun 
        integer(I8)                                     :: nvars, nd
        integer(I8), allocatable                        :: vars(:)
        real(R8), allocatable                           :: d(:)

    contains

        ! Constructor
        procedure :: Initialize     => InitializeFDchecker

        ! Gradient checker
        procedure :: CheckGradient

        ! Hessian checker
        procedure :: CheckHessian

        ! Printing
        procedure :: PrintGradientHeader
        procedure :: PrintGradientIterate
        procedure :: PrintHessianHeader
        procedure :: PrintHessianIterate

        ! Housekeeping procedures
        procedure :: Allocate       => AllocateFDchecker
        procedure :: Deallocate     => DeallocateFDchecker

        ! Destructor
        final :: DestroyFDchecker

    end type

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                       DIAGNOSTIC FUNCTION                        !
    !------------------------------------------------------------------!
    abstract interface

        ! Evaluation
        subroutine EvaluateDiagnosticsFunctionINT(fun, x, f, df, d2f, &
            dogradient, dohessian)

            ! Description
            !============
            ! This subroutine evaluates the function, given the input 
            ! argument x, and returns the value, f, the gradient, df, 
            ! and the hessian, d2f. dogradient and dohessian
            ! are logicals that can be set to true or false in order to
            ! specify whether gradients and/or hessians should actually
            ! be computed. 

            ! Import
            import DiagnosticsFunctionUDT, R8, MySparseUDT 

            ! Declare
            class(DiagnosticsFunctionUDT)       :: fun
            real(R8)                            :: f
            real(R8), allocatable               :: x(:), df(:)
            type(MySparseUDT)                   :: d2f 
            logical                             :: dogradient, dohessian

        end subroutine

        ! Dimension getter
        subroutine GetDimensionsDiagnosticsFunctionINT(fun, dimx)

            ! Description
            !============
            ! This subroutine returns the dimension of x, which also 
            ! fixes the gradient and hessian dimensions. 

            ! Import
            import :: DiagnosticsFunctionUDT, I8

            ! Declare
            class(DiagnosticsFunctionUDT)       :: fun 
            integer(I8)                         :: dimx 

        end subroutine

        ! Argument getter
        subroutine GetArgumuentsDiagnosticsFunctionINT(fun, x)

            ! Description
            !============
            ! This function returns the arguments of the function, x. 

            ! Import
            import :: DiagnosticsFunctionUDT, R8

            ! Declare
            class(DiagnosticsFunctionUDT)       :: fun 
            real(R8), allocatable               :: x(:)


        end subroutine

    end interface


    contains

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!    

    !------------------------------------------------------------------!
    !                              FD CHECKER                          !
    !------------------------------------------------------------------!

    ! Constructor
    subroutine InitializeFDchecker(FDchecker, nvars, vars) 

        ! Description
        !============
        ! Initialization routine for the finite difference checker, 
        ! which sets the nvars and vars values. These have to be 
        ! specified by the user. Externally, the user should also set 
        ! the FDchecker%fun field by allocation with his own type. 

        ! Initialize
        !===========
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        class(FDcheckerUDT)             :: FDchecker 
        integer(I8), intent(in)         :: nvars
        integer(I8), intent(in)         :: vars(1:nvars)

        ! Initialize
        !===========
        ! Initialize
        FDchecker%nvars = nvars 
        FDchecker%nd = 5

        ! Allocate
        call FDchecker%Allocate()

        ! Set values
        FDchecker%vars = vars  

        ! Set the step sizes (hard coded - these values are fine)
        FDchecker%d = [1e-2, 1e-4, 1e-6, 1e-8, 1e-10]

    end subroutine

    ! Gradient checker
    subroutine CheckGradient(FDchecker)

        ! Description
        !============
        ! This routine computes the gradient and checks by doing a 
        ! finite difference step study whether the error reduces as
        ! expected. The results are printed out to the terminal. A 
        ! central difference scheme is used to compute the derivatives. 

        ! Declare variables
        !==================
        ! Arguments
        class(FDcheckerUDT)         :: FDchecker 

        ! Loop
        integer(I8)                 :: i, j

        ! Auxiliary
        real(R8)                    :: fref, fFW, fBW, gradFW, gradBW, &
            gradC, reFW, reBW, reC
        real(R8), allocatable       :: gradref(:), xref(:), xFW(:), &
            xBW(:)

        integer(I8)                 :: dimx, tv

        logical                     :: dogradient, dohessian

        ! Data

        ! Dummy
        type(MySparseUDT)           :: dummyhess
        real(R8), allocatable       :: dummygrad(:)

        ! Initialize
        !===========
        ! Set logicals
        dogradient = .false.
        dohessian = .false.

        ! Get the dimensions of x
        call FDchecker%fun%GetDimensions(dimx)

        ! Allocate
        allocate(xref(dimx), xFW(dimx), xBW(dimx), &
            gradref(dimx), dummygrad(dimx))

        ! Compute reference and its gradient
        call FDchecker%fun%GetArguments(xref) ! current point
        call FDchecker%fun%Evaluate(xref, fref, gradref, dummyhess, &
            .true., .false.)

        ! Print the header
        call FDchecker%PrintGradientHeader()

        ! Compute the finite difference gradients
        !========================================
        ! Associate
        associate(&
            d       => FDchecker%d,     &
            nd      => FDchecker%nd,    &
            vars    => FDchecker%vars,  &
            nvars   => FDchecker%nvars  &
            )

        ! Loop over all variables
        do j = 1, nvars
            ! Get current variable
            tv = vars(j)
            ! Loop over all step sizes
            do i = 1, FDchecker%nd
                ! Initialize
                xFW     = xref 
                xBW     = xref 

                ! Perturb
                xFW(tv) = xref(tv) + d(i)
                xBW(tv) = xref(tv) - d(i)

                ! Compute function values
                call FDchecker%fun%Evaluate(xFW, fFW, dummygrad, &
                    dummyhess, .false., .false.)
                call FDchecker%fun%Evaluate(xBW, fBW, dummygrad, &
                    dummyhess, .false., .false.)

                ! Compute FD gradients
                gradFW  = (fFW - fref)/d(i)
                gradBW  = (fref - fBW)/d(i)
                gradC   = 0.5*(gradFW + gradBW)

                ! Compute relative error
                call ComputeScalarRelativeError(gradref(tv), gradFW, &
                    reFW)
                call ComputeScalarRelativeError(gradref(tv), gradBW, &
                    reBW)
                call ComputeScalarRelativeError(gradref(tv), gradC, &
                    reC)

                ! Print
                call FDchecker%PrintGradientIterate(tv, d(i), &
                    gradref(tv), gradFW, gradBW, gradC, reFW, reBW, &
                    reC)
            end do
        end do

        ! End associate
        end associate

        ! Deallocate
        deallocate(xref, xFW, xBW, gradref)


    end subroutine

    ! Hessian checker
    subroutine CheckHessian(FDchecker)

        ! Description
        !============
        ! This routine computes the hessian and checks by doing a 
        ! finite difference step study whether the error reduces as
        ! expected. The results are printed out to the terminal. A 
        ! central difference scheme is used to compute the derivatives. 
        ! 
        ! Notes
        !======
        ! Note 1: since the 'hessian' here will be found by computing a
        ! finite difference of the gradient, computed using the user 
        ! provided funciton, the outcome will depend on whether these
        ! gradients are correctly computed! 
        
        ! Note 2: what is actually checked, are the columns of the 
        ! hessian that correspond to the current variable. Note that it 
        ! is NOT explicitly checked whether the hessian is symmetric 
        ! (i.e. whether the row is the same as the column transposed)! 
        ! The true error criterium is the maximal difference between the
        ! non-zero values of the hessian column. 

        ! Declare variables
        !==================
        ! Arguments
        class(FDcheckerUDT)         :: FDchecker 

        ! Loop
        integer(I8)                 :: i, j

        ! Auxiliary
        real(R8)                    :: fref, fFW, fBW, &
            gradC, reFW, reBW, reC
        integer(I8)                 ::  indexFW(1), indexBW(1), indexC(1), index
        real(R8), allocatable       :: gradref(:), xref(:), xFW(:), &
            xBW(:), hesscolFW(:), hesscolBW(:), hesscolC(:), &
            hesscolref(:), gradFW(:), gradBW(:)
        type(MySparseUDT)           ::  hessref

        integer(I8)                 :: dimx, tv

        ! Data
        real(R8)                    :: macheps = 1e-15

        ! Dummy
        type(MySparseUDT)           :: dummyhess

        ! Initialize
        !===========
        ! Get the dimensions of x
        call FDchecker%fun%GetDimensions(dimx)

        ! Allocate
        allocate(xref(dimx), xFW(dimx), xBW(dimx), &
            gradref(dimx), hesscolFW(dimx), hesscolBW(dimx), &
            hesscolC(dimx), hesscolref(dimx), gradFW(dimx), &
            gradBW(dimx))

        ! Compute reference and its gradient
        call FDchecker%fun%GetArguments(xref) ! current point
        call FDchecker%fun%Evaluate(xref, fref, gradref, hessref, &
            .true., .true.)

        ! Print the header
        call FDchecker%PrintHessianHeader()

        ! Compute the finite difference gradients
        !========================================
        ! Associate
        associate(&
            d       => FDchecker%d,     &
            nd      => FDchecker%nd,    &
            vars    => FDchecker%vars,  &
            nvars   => FDchecker%nvars  &
            )

        ! Loop over all variables
        do j = 1, nvars
            ! Get current variable
            tv = vars(j)

            ! Extract the reference
            call hessref%ExtractColumnFull(hesscolref, tv)

            ! Loop over all step sizes
            do i = 1, FDchecker%nd
                ! Initialize
                xFW     = xref 
                xBW     = xref 

                ! Perturb
                xFW(tv) = xref(tv) + d(i)
                xBW(tv) = xref(tv) - d(i)

                ! Compute function gradients
                gradFW(:) = 0 ! very important to re-initialize!
                gradBW(:) = 0
                call FDchecker%fun%Evaluate(xFW, fFW, gradFW, &
                    dummyhess, .true., .false.)
                call FDchecker%fun%Evaluate(xBW, fBW, gradBW, &
                    dummyhess, .true., .false.)

                ! Compute FD hessian columns
                hesscolFW(:) = 0
                hesscolBW(:) = 0
                hesscolC(:) = 0
                hesscolFW  = (gradFW - gradref)/d(i)
                hesscolBW  = (gradref - gradBW)/d(i)
                hesscolC   = 0.5*(hesscolBW + hesscolFW)

                ! Compute relative error
                call ComputeVectorRelativeError(hesscolref, hesscolFW, &
                    reFW, indexFW)
                call ComputeVectorRelativeError(hesscolref, hesscolBW, &
                    reBW, indexBW)
                call ComputeVectorRelativeError(hesscolref, hesscolC, &
                    reC, indexC)

                ! Print - we only print the (tv, tv) component!
                index = indexC(1)
                call FDchecker%PrintHessianIterate(tv, d(i), &
                    hesscolref(index), hesscolFW(index), hesscolBW(index), &
                    hesscolC(index), reFW, reBW, reC)
                
            end do
        end do

        ! End associate
        end associate

        ! Deallocate
        deallocate(xref, xFW, xBW, gradref, hesscolFW, hesscolBW, &
            hesscolC, hesscolref, gradFW, gradBW)

    end subroutine

    ! Allocation
    subroutine AllocateFDchecker(FDchecker)

        ! Description
        !============
        ! Allocate the fields (except fun) of the FD checker. Assumed 
        ! that nvars is known and set. 

        ! Declare variables
        !==================
        ! Arguments
        class(FDcheckerUDT)             :: FDchecker 

        ! Allocate
        !=========
        allocate(FDchecker%vars(FDchecker%nvars))
        allocate(FDchecker%d(FDchecker%nd))


    end subroutine

    ! Deallocation
    subroutine DeallocateFDchecker(FDchecker)

        ! Description
        !============
        ! Deallocator

        ! Declare variables
        !==================
        ! Arguments 
        class(FDcheckerUDT)             :: FDchecker 

        ! Deallocate
        !===========
        if (allocated(FDchecker%vars)) then
            deallocate(FDchecker%vars)
        end if
        if (allocated(FDchecker%d)) then
            deallocate(FDchecker%d)
        end if

    end subroutine

    ! Destructor
    subroutine DestroyFDchecker(FDchecker)

        ! Description
        !============
        ! Deallocator

        ! Declare variables
        !==================
        ! Arguments 
        type(FDcheckerUDT)              :: FDchecker 

        ! Deallocate
        !===========
        call FDchecker%Deallocate()

    end subroutine

    ! Print gradient header
    subroutine PrintGradientHeader(FDchecker)

        ! Description
        !============
        ! Print out the header of the gradient checker

        ! Declare variables
        !==================
        ! Arguments
        class(FDcheckerUDT)             :: FDchecker 

        ! Print
        !======
        print *, '!===================================================!'
        print *, '!         Gradient Finite Difference checking       !'
        print *, '!===================================================!'
        print "(4x, a4, 4x, 4x, a4, 4x, 4x, a4, 4x, 4x, a4, 4x, &
            & 4x, a4, 4x, 4x, a4, 4x, 4x, a4, 4x, 4x, a4, 4x, 4x, a4, 4x)" &
            ,   'var ', 'step',    'gref ',    'gFW ', 'gBW ', 'gC ', &
            'reFW', 'reBW', 'reC '

    end subroutine

    ! Print Hessian header
    subroutine PrintHessianHeader(FDchecker)

        ! Description
        !============
        ! Print out the header of the Hessian checker

        ! Declare variables
        !==================
        ! Arguments
        class(FDcheckerUDT)             :: FDchecker 

        ! Print
        !======
        print *, '!===================================================!'
        print *, '!         Hessian Finite Difference checking        !'
        print *, '!===================================================!'
        print "(4x, a4, 4x, 4x, a4, 4x, 4x, a4, 4x, 4x, a4, 4x, &
            & 4x, a4, 4x, 4x, a4, 4x, 4x, a4, 4x, 4x, a4, 4x, 4x, a4, 4x)" &
            ,   'var ', 'step',    'href ',    'hFW ', 'hBW ', 'hC ', &
            'reFW', 'reBW', 'reC '

    end subroutine

    ! Print gradient iterate
    subroutine PrintGradientIterate(FDchecker, var, step, gref, gFW, &
        gBW, gC, reFW, reBW, reC)

        ! Description
        !============
        ! Print an iterate of the gradient checker

        ! Declare variables
        !==================
        ! Arguments
        class(FDcheckerUDT)                 :: FDchecker 
        integer(I8), intent(in)             :: var 
        real(R8), intent(in)                :: gref, gFW, gBW, gC, &
            reFW, reBW, reC, step 
        
        ! Print
        !======
        print "(i8, 4x, e8.2, 4x, e8.2, 4x, e8.2, 4x, e8.2,&
         & 4x, e8.2, 4x, e8.2, 4x, e8.2, 4x, e8.2, 4x)", &
            var, step, gref, gFW, gBW, gC, reFW, reBW, reC

    end subroutine

    ! Print Hessian iterate
    subroutine PrintHessianIterate(FDchecker, var, step, gref, gFW, &
        gBW, gC, reFW, reBW, reC)

        ! Description
        !============
        ! Print an iterate of the gradient checker

        ! Declare variables
        !==================
        ! Arguments
        class(FDcheckerUDT)                 :: FDchecker 
        integer(I8), intent(in)             :: var 
        real(R8), intent(in)                :: gref, gFW, gBW, gC, &
            reFW, reBW, reC, step 
        
        ! Print
        !======
        print "(i8, 4x, e8.2, 4x, e8.2, 4x, e8.2, 4x, e8.2,&
         & 4x, e8.2, 4x, e8.2, 4x, e8.2, 4x, e8.2, 4x)", &
            var, step, gref, gFW, gBW, gC, reFW, reBW, reC

    end subroutine

    !------------------------------------------------------------------!
    !                        ERROR COMPUTATION                         !
    !------------------------------------------------------------------!

    ! Absolute error computation 
    subroutine ComputeScalarAbsoluteError(a, b, err)

        ! Description
        !============
        ! Compute the absolute error between the scalars a and b

        ! Declare variables
        !===================
        ! Arguments
        real(R8), intent(in)        :: a, b 
        real(R8), intent(out)       :: err
        
        ! Compute
        !========
        err = abs(a - b)

    end subroutine

    subroutine ComputeVectorAbsoluteError(a, b, err)

        ! Description
        !============
        ! Compute the absolute error between the vectors a and b based 
        ! on the infinity norm. 

        ! Declare variables
        !===================
        ! Arguments
        real(R8), intent(in), allocatable       :: a(:), b(:) 
        real(R8), intent(out)                   :: err
        
        ! Compute
        !========
        err = maxval(abs(a - b))
        

    end subroutine

    ! Relative error computation 
    subroutine ComputeScalarRelativeError(a, b, err)

        ! Description
        !============
        ! Compute the relative error between the scalars a and b w.r.t.
        ! the reference a

        ! Declare variables
        !===================
        ! Arguments
        real(R8), intent(in)        :: a, b 
        real(R8), intent(out)       :: err
        
        ! Compute
        !========
        ! Compute absolute error
        call ComputeScalarAbsoluteError(a, b, err)

        ! Checks
        if (a == 0) then 
            if (b == 0) then 
                ! perfectly fine
                err = 0
            else
                ! Should be infinity, return huge
                err = huge(R8)
            end if
        else
            ! Compute
            err = err/a 
        end if

    end subroutine

    subroutine ComputeVectorRelativeError(a, b, err, location)

        ! Description
        !============
        ! Compute the relative error between the scalars a and b w.r.t.
        ! the reference a

        ! Declare variables
        !===================
        ! Arguments
        real(R8), intent(in), allocatable       :: a(:), b (:)
        real(R8), intent(out)                   :: err
        integer(I8)                             :: location(1)

        ! Auxiliary
        real(R8), allocatable                   :: errvec(:)
        
        ! Compute
        !========
        ! Allocate
        allocate(errvec(size(a)))

        ! Compute absolute error vector
        where ((a == 0) .and. (b == 0)) errvec = 0 ! perfectly fine
        where ((a == 0) .and. (b .ne. 0)) errvec = huge(R8) ! infinite
        where ((a .ne. 0) .and. (b .ne. 0)) errvec = (a - b)/a

        ! Compute maximal value
        err = maxval(errvec)
        location = maxloc(errvec)

        ! Deallocate
        deallocate(errvec)

    end subroutine


end module