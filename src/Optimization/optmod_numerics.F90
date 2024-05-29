!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains the numerical parameters for the optimization. 

module optmod_numerics
    
    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_readwrite
    use mod_inputfileparser

    ! The usual
    implicit none
    save
    public 

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    ! Abstract types
    !===============
    ! General numerical parameter object
    type :: NumUDT

        ! Description
        !============
        ! Abstract type that contains the general numerical parameters
        ! and routines (if any) for optimization solvers. 
        ! Fields:
        ! - tol:            tolerance to which the solver has to solve
        ! - maxit:          maximum number of iterations
        ! - verbosity:      verbosity level of printing
        ! - inputfilepath:  filepath from which user-specified data 
        !                   can be read in.
        ! - fieldprefix: all numerical options are read in assuming a 
        ! format '<fieldprefix>opt.num.<field>'. Fieldprefix can be 
        ! empty or given (default is empty)

        character(:), allocatable   :: inputfilepath
        character(:), allocatable   :: fieldprefix

        real(R8)            :: tol 
        integer(I8)         :: maxit 
        integer(I8)         :: verbosity

    contains 

        procedure :: SetDefaultNumParams => SetDefaultNumParamsINT
        procedure :: InitializeNumParams => InitializeGeneralNumParamsINT

    end type

    ! Derived types
    !==============
    ! KKT numerics
    type, extends(NumUDT) :: NumKKTUDT 

        ! Description
        !============
        ! Numerical parameters for the KKT solver. The following fields
        ! are added:
        ! - rxf:    relaxation factor on the cost function hessian 
        !           (sum(abs(hessJ,2))*unity matrix*rxf is added). 
        !           Higher rxf -> more relaxation -> (normally) better 
        !           convergence. Set to zero for no relaxation.
        ! - rxfdesign:  relaxation factor on the design update. Rescales
        !               the design update as 
        !               x_new = x_old + rxfdesign*delta x. 
        !               smaller values -> slower convergence, yet more
        !               stable. Set to one for no relaxation. 
        ! - rxfdec: relaxation reduction factor per iteration, i.e. the
        !           next iteration, rxf_new = rxf_old*rxfdec. Set to one
        !           for constant relaxation factor. Set smaller to one 
        !           to reduce relaxation factor (0.98 works good in most
        !           cases)
        ! - rxfmin: minimal value of relaxation factor
        ! - fieldprefix: all numerical options are read in assuming a 
        ! format '<fieldprefix>opt.num.<field>'. Fieldprefix can be 
        ! empty or given (default is empty)

        ! Relaxation factors
        real(R8)            :: rxf
        real(R8)            :: rxfdec 
        real(R8)            :: rxfmin
        real(R8)            :: rxfdesign

    contains 

        procedure :: SetDefaultNumParamsKKT => SetDefaultNumParamsKKTINT
        procedure :: InitializeNumParams => InitializeNumParamsKKTINT
        procedure :: Read               => ReadNumKKTOptions

    end type

    ! Line search numerics
    type, extends(NumUDT) :: NumLSUDT

        ! Description
        !============
        ! This type contains all numerical options for linesearch 
        ! methods. The following fields are present:
        ! - type:   type of linesearch, e.g. 'wolfe' or 'backtracking'
        ! - meritfunction:  type of merit function to use, e.g. 'l1'
        ! - dolinesearch:   logical to decide whether to do a linesearch
        ! - maxit:  maximal number of linesearch iterations
        ! - dec:    decrease factor if cost function reduction is insufficient 
        ! - inc:    increase factor if gradient condition is insufficient
        ! - c1:     constant for first condition
        ! - c2:     constant for second condition
        ! - mfdelta:    constant for merit function 
        ! - mfpenfactol:    tolerance under which the default penalization 
        !                   factor becomes active for the merit function
        ! - mfpenfacdef:    default penalization factor, used in
        !                   initial iterate


        character(:), allocatable       :: type, meritfunction 
        logical                         :: dolinesearch
        real(R8)                        :: dec, inc, c1, c2, mfdelta, &
            mfpenfactol, mfpenfacdef



    contains

        procedure :: SetDefaultNumParamsLS  
        procedure :: InitializeNumParams    => InitializeNumParamsLS 
        procedure :: Read                   => ReadNumLSOptions

    end type

    ! Quasi-Newton numerics
    type, extends(NumUDT) :: NumQNUDT 

        
    contains 

        procedure :: SetDefaultNumParamsQN 
        procedure :: InitializeNumParams    => InitializeNumParamsQN 
        procedure :: Read                   => ReadNumQNOptions

    end type

    ! Nonlinear complementarity problem numerics
    type :: NumNCPUDT 

        ! Description
        !============
        ! Numerics structure for ncp function. Not extended from 
        ! general numerics type, because does not need those fields. 
        ! Two fields are present, alpha and ncpfun. The latter determines
        ! which type of nonlinear complementarity problem function is
        ! taken (can be 'max', 'FB' (Fisher-Burmeister), or 'FBsmooth').
        ! The first is a numerical parameter that is used for the 'max'
        ! and 'FBsmooth' ncpfun options. 

        ! Fields
        real(R8)                    :: alpha 
        character(:), allocatable   :: ncpfun 
        character(:), allocatable   :: inputfilepath
        character(:), allocatable   :: fieldprefix

    contains

        procedure :: SetDefaultNumParams => SetDefaultNumParamsNCP
        procedure :: InitializeNumParams => InitializeNumParamsNCP
        procedure :: Read                => ReadNumNCPOptions

    end type

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    contains

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                           GENERAL OPTIONS                        !
    !------------------------------------------------------------------!
    ! Set the general numerical parameters
    subroutine SetDefaultNumParamsINT(num)

        ! Description
        !============
        ! Set default numerical parameters tol, itmax, verbosity

        ! Declare variables
        !==================
        ! Arguments
        class(NumUDT)            :: num
    
        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Set defaults
        !=============
        num%tol         = 1e-6
        num%maxit       = 5
        num%verbosity   = 1

    end subroutine

    ! Initialize the numerics
    subroutine InitializeGeneralNumParamsINT(num)

        ! Description
        !============
        ! Initialize the general numerical parameters

        ! Declare variables
        !==================
        ! Arguments
        class(NumUDT)            :: num
    
        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Set defaults
        call num%SetDefaultNumParams()

        ! Override with user settings (to be implemented)

    end subroutine

    !------------------------------------------------------------------!
    !                             KKT OPTIONS                          !
    !------------------------------------------------------------------!
    ! Set the KKT numerical parameters
    subroutine SetDefaultNumParamsKKTINT(num)

        ! Description
        !============
        ! Set default numerical parameters tol, itmax, verbosity

        ! Declare variables
        !==================
        ! Arguments
        class(NumKKTUDT)            :: num
    
        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Set defaults
        !=============
        ! General numerics
        call num%SetDefaultNumParams()

        ! Specifics for KKT numerics
        num%rxf             = 1e0
        num%rxfdesign       = 1
        num%rxfdec          = 0.98
        num%rxfmin          = 2e-2

    end subroutine

    ! Initialize the numerics
    subroutine InitializeNumParamsKKTINT(num)

        ! Description
        !============
        ! Initialize the general numerical parameters

        ! Declare variables
        !==================
        ! Arguments
        class(NumKKTUDT)            :: num
    
        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Set defaults 
        call num%SetDefaultNumParamsKKT()

        ! Override with user settings (to be implemented)
        call num%Read()

    end subroutine

    ! Read user data
    subroutine ReadNumKKTOptions(num)

        ! Description
        !============
        ! Read in grid options from file. It is assumed that the 
        ! filepath has been set correctly. 

        ! Declare variables
        !==================
        ! Arguments
        class(NumKKTUDT)                 :: num 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: field
        integer, parameter              :: fid = 10 
        logical                         :: reachedeof

        ! Initialize
        !===========
        ! Variables
        reachedeof = .false. 

        ! Open the file, check if it exists
        open(unit=fid, file=num%inputfilepath, status='old', &
            iostat=openstatus)

        if (openstatus > 0) then 
            ! Something wrong when reading file - continue with default
            ! values
            print *, 'ReadNumKKTOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadNumKKTOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! General
        print *, num%fieldprefix 
        field = num%fieldprefix // 'opt.num.itmax'
        call ExtractOptionValueInteger0D(fid, field, num%maxit)
        field = num%fieldprefix // 'opt.num.verbosity'
        call ExtractOptionValueInteger0D(fid, field, num%verbosity)
        field = num%fieldprefix // 'opt.num.tol'
        call ExtractOptionValueReal0D(fid, field, num%tol)
        
        ! Relaxation factors
        field = num%fieldprefix // 'opt.num.rxf'
        call ExtractOptionValueReal0D(fid, field, num%rxf)
        field = num%fieldprefix // 'opt.num.rxfdec'
        call ExtractOptionValueReal0D(fid, field, num%rxfdec)
        field = num%fieldprefix // 'opt.num.rxfmin'
        call ExtractOptionValueReal0D(fid, field, num%rxfmin)
        field = num%fieldprefix // 'opt.num.rxfdesign'
        call ExtractOptionValueReal0D(fid, field, num%rxfdesign)

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)

    end subroutine

    !------------------------------------------------------------------!
    !                         LINESEARCH OPTIONS                       !
    !------------------------------------------------------------------!
    ! Set the KKT numerical parameters
    subroutine SetDefaultNumParamsLS(num)

        ! Description
        !============
        ! Set default numerical parameters tol, itmax, verbosity

        ! Declare variables
        !==================
        ! Arguments
        class(NumLSUDT)            :: num
    
        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Set defaults
        !=============
        ! General numerics
        call num%SetDefaultNumParams()

        ! Specifics for linesearch numerics (recommended)
        num%type            = 'backtracking_soc' 
        num%meritfunction   = 'l1'
        num%dolinesearch    = .true.
        num%c1              = 1e-4
        num%c2              = 0.9
        num%dec             = 0.5
        num%inc             = 2
        num%maxit           = 8 
        num%mfdelta         = 1e-4 
        num%mfpenfactol     = 1e-12
        num%mfpenfacdef     = 1

    end subroutine

    ! Initialize the numerics
    subroutine InitializeNumParamsLS(num)

        ! Description
        !============
        ! Initialize the general numerical parameters

        ! Declare variables
        !==================
        ! Arguments
        class(NumLSUDT)            :: num
    
        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Set defaults 
        call num%SetDefaultNumParamsLS()

        ! Override with user settings (to be implemented)
        call num%Read()

    end subroutine

    ! Read user data
    subroutine ReadNumLSOptions(num)

        ! Description
        !============
        ! Read in grid options from file. It is assumed that the 
        ! filepath has been set correctly. 

        ! Declare variables
        !==================
        ! Arguments
        class(NumLSUDT)                 :: num 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: field
        integer, parameter              :: fid = 10 
        logical                         :: reachedeof

        ! Initialize
        !===========
        ! Variables
        reachedeof = .false. 

        ! Open the file, check if it exists
        open(unit=fid, file=num%inputfilepath, status='old', &
            iostat=openstatus)

        if (openstatus > 0) then 
            ! Something wrong when reading file - continue with default
            ! values
            print *, 'ReadNumLSOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadNumLSOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! General
        print *, num%fieldprefix 
        field = num%fieldprefix // 'opt.num.ls.type'
        call ExtractOptionValueCharacter(fid, field, num%type)
        field = num%fieldprefix // 'opt.num.ls.meritfunction'
        call ExtractOptionValueCharacter(fid, field, num%meritfunction)
        field = num%fieldprefix // 'opt.num.ls.dolinesearch'
        call ExtractOptionValueLogical0D(fid, field, num%dolinesearch)
        field = num%fieldprefix // 'opt.num.ls.itmax'
        call ExtractOptionValueInteger0D(fid, field, num%maxit)
        field = num%fieldprefix // 'opt.num.ls.verbosity'
        call ExtractOptionValueInteger0D(fid, field, num%verbosity)
        field = num%fieldprefix // 'opt.num.ls.dec'
        call ExtractOptionValueReal0D(fid, field, num%dec)
        field = num%fieldprefix // 'opt.num.ls.inc'
        call ExtractOptionValueReal0D(fid, field, num%inc)
        field = num%fieldprefix // 'opt.num.ls.c1'
        call ExtractOptionValueReal0D(fid, field, num%c1)
        field = num%fieldprefix // 'opt.num.ls.c2'
        call ExtractOptionValueReal0D(fid, field, num%c2)
        field = num%fieldprefix // 'opt.num.ls.mf.delta'
        call ExtractOptionValueReal0D(fid, field, num%mfdelta)
        field = num%fieldprefix // 'opt.num.ls.mf.penfactol'
        call ExtractOptionValueReal0D(fid, field, num%mfpenfactol)
        field = num%fieldprefix // 'opt.num.ls.mf.penfacdef'
        call ExtractOptionValueReal0D(fid, field, num%mfpenfacdef)

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)

    end subroutine

    !------------------------------------------------------------------!
    !                        NCP FUNCTION OPTIONS                      !
    !------------------------------------------------------------------!
    ! Set the KKT numerical parameters
    subroutine SetDefaultNumParamsNCP(num)

        ! Description
        !============
        ! Set default numerical parameters tol, itmax, verbosity

        ! Declare variables
        !==================
        ! Arguments
        class(NumNCPUDT)            :: num
    
        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Set defaults
        !=============
        num%alpha = 1
        num%ncpfun = 'max'

    end subroutine

    ! Initialize the numerics
    subroutine InitializeNumParamsNCP(num)

        ! Description
        !============
        ! Initialize the general numerical parameters

        ! Declare variables
        !==================
        ! Arguments
        class(NumNCPUDT)            :: num
    
        ! Loop variables

        ! Auxiliary variables 

        ! Data

        ! Set defaults 
        call num%SetDefaultNumParams()

        ! Override with user settings (to be implemented)
        call num%Read()

    end subroutine

    ! Read user data
    subroutine ReadNumNCPOptions(num)

        ! Description
        !============
        ! Read in grid options from file. It is assumed that the 
        ! filepath has been set correctly. 

        ! Declare variables
        !==================
        ! Arguments
        class(NumNCPUDT)                 :: num 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: field
        integer, parameter              :: fid = 10 
        logical                         :: reachedeof

        ! Initialize
        !===========
        ! Variables
        reachedeof = .false. 

        ! Open the file, check if it exists
        open(unit=fid, file=num%inputfilepath, status='old', &
            iostat=openstatus)

        if (openstatus > 0) then 
            ! Something wrong when reading file - continue with default
            ! values
            print *, 'ReadNumNCPOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadNumNCPOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! General
        print *, num%fieldprefix 
        field = num%fieldprefix // 'opt.num.ncp.ncpfun'
        call ExtractOptionValueCharacter(fid, field, num%ncpfun)
        field = num%fieldprefix // 'opt.num.ncp.alpha'
        call ExtractOptionValueReal0D(fid, field, num%alpha)

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)

    end subroutine

    !------------------------------------------------------------------!
    !                         QN SOLVER OPTIONS                        !
    !------------------------------------------------------------------!
    ! Set numerical parameters
    subroutine SetDefaultNumParamsQN(num)

        ! Description
        !============
        ! Set default numerical parameters of quasi newton solver

        ! Declare variables
        !==================
        class(NumQNUDT)         :: num 

    end subroutine

    ! Initialize
    subroutine InitializeNumParamsQN(num)

        ! Description
        !============
        ! Initialization routine that calls the default parameter setter
        ! and the reader

        ! Declare variables
        !==================
        class(NumQNUDT)         :: num
        
        ! Set parameters
        !===============
        ! Set defaults 
        call num%SetDefaultNumParams()

        ! Override with user settings (to be implemented)
        call num%Read()

    end subroutine 

    ! Read
    subroutine ReadNumQNOptions(num)

        ! Declare variables
        !==================
        class(NumQNUDT)         :: num 

    end subroutine

end module