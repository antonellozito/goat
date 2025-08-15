!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains all the user-defined input parameters and related
! routines for the grid adaptation module. 

module gamod_userinput

    ! Initialize
    !============
    ! Load modules
    use goatmod_userinput
    use mod_inputfileparser
    
    
    ! The usual
    implicit none
    save
    public 

    ! Read namelist from file

    ! All functions

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!    
    
    ! Options to active the different possible adaptation operations
    type, extends(OptionsUDT) :: OperationsOptionsUDT


        ! TODO: explanation on all options!!!!
        logical         :: rem_small_trias
        real(R8)        :: cut_off_pol
        real(R8)         :: cut_off_surf

        logical         :: stacked_trias
        logical         :: stacked_trias_checkAR
        real(R8)        :: stacked_trias_maxAR
        logical         :: merge_stacked_trias
        real(R8)        :: merge_stacked_trias_angle_threshold
        logical         :: merge_trap_into_stacked

        logical         :: stacked_to_cutcell
        logical         :: stacked_to_cutcell_uniform

        logical         :: split_shaved_off_tube


        logical         :: splitting
        logical         :: merging
        logical         :: pents_to_tria
        real(R8)        :: h_rad_threshold
        real(R8)        :: h_rad_core_threshold

        logical         :: rem_stickout_trias
        logical         :: rem_trias_flux
        real(R8)        :: rem_tube_outershell_threshold
        character(:), allocatable   :: outershell_handling
        logical         :: remove_stickoutquad
        logical         :: split_noalignedquads


    contains

        procedure :: SetDefaults     => SetDefaultOperationOptions
        procedure :: Read           => ReadOperationOptions

    end type

    ! Options to specify the splitting operations
    type, extends(OptionsUDT) :: SplittingOptionsUDT
    contains
  
        procedure :: SetDefaults     => SetDefaultSplittingOptions
        procedure :: Read           => ReadSplittingOptions

    end type    

    ! Options to specify the merging operations
    type, extends(OptionsUDT) :: MergingOptionsUDT
    contains
  
        procedure :: SetDefaults     => SetDefaultMergingOptions
        procedure :: Read           => ReadMergingOptions

    end type
    
    ! Options to specify the pentagon setting
    type, extends(OptionsUDT) :: PentagonOptionsUDT
    contains
  
        procedure :: SetDefaults     => SetDefaultPentagonOptions
        procedure :: Read           => ReadPentagonOptions

    end type
    
    ! Options to specify the distance function settings
    type, extends(OptionsUDT) :: DistancefunctionOptionsUDT
    contains
  
        procedure :: SetDefaults     => SetDefaultDistancefunctionOptions
        procedure :: Read           => ReadDistancefunctionOptions

    end type     
    !==================================================================!
    !                                                                  !
    !                              INTERFACES                          !
    !                                                                  !
    !==================================================================!
    
    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    contains 

    !------------------------------------------------------------------!
    !                            Option setters                        !
    !------------------------------------------------------------------!

    subroutine SetDefaultOperationOptions(options)

        ! Description
        !============
        ! Set the operation options. 
        
        ! Declare variables
        !==================
        ! Arguments
        class(OperationsOptionsUDT)    :: options

        ! Default options
        !================
        options%rem_small_trias                     = .false.
        options%cut_off_pol                         = 0.3
        options%cut_off_surf                        = 0.05
    
        options%stacked_trias                       = .false.
        options%stacked_trias_checkAR               = .false.
        options%stacked_trias_maxAR                 = 10
        options%merge_stacked_trias                 = .false.
        options%merge_stacked_trias_angle_threshold = 5
        options%merge_trap_into_stacked             = .true.

        options%stacked_to_cutcell                  = .false.
        options%stacked_to_cutcell_uniform          = .true.

        options%split_shaved_off_tube               = .false.

        options%splitting                           = .false.
        options%pents_to_tria                       = .false.
        options%merging                             = .false.
        options%h_rad_threshold                     = 0.01
        options%h_rad_core_threshold                = 0.04

        options%rem_stickout_trias                  = .false.
        options%rem_trias_flux                      = .false.
        options%rem_tube_outershell_threshold       = 2
        options%outershell_handling                 = 'merge' 
        options%remove_stickoutquad                 = .false.
        options%split_noalignedquads                = .true.   
        
    end subroutine 
    
    subroutine SetDefaultSplittingOptions(options)
        ! Description
        !============
        ! Set the splitting options. 
        
        ! Declare variables
        !==================
        ! Arguments
        class(SplittingOptionsUDT)    :: options

        ! Default options
        !================        
    end subroutine

    subroutine SetDefaultMergingOptions(options)
        ! Description
        !============
        ! Set the merging options. 
        
        ! Declare variables
        !==================
        ! Arguments
        class(MergingOptionsUDT)    :: options

        ! Default options
        !================        
    end subroutine

    subroutine SetDefaultPentagonOptions(options)
        ! Description
        !============
        ! Set the pentagon options. 
        
        ! Declare variables
        !==================
        ! Arguments
        class(PentagonOptionsUDT)    :: options

        ! Default options
        !================        
    end subroutine

    subroutine SetDefaultDistancefunctionOptions(options)
        ! Description
        !============
        ! Set the operation options. 
        
        ! Declare variables
        !==================
        ! Arguments
        class(DistancefunctionOptionsUDT)    :: options

        ! Default options
        !================        
    end subroutine

    !------------------------------------------------------------------!
    !                            Option readers                        !
    !------------------------------------------------------------------!   
    
    subroutine ReadOperationOptions(options)

        ! Description
        !============
        ! Read in user-specified adaption operation options 
        
        ! Declare variables
        !==================
        ! Arguments
        class(OperationsOptionsUDT)  :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: field
        integer, parameter              :: fid = 10        
        
        ! Initialize
        !===========
                
        ! Open the file, check if it exists
        open(unit=fid, file=options%inputfilepath, status='old', &
            iostat=openstatus)

        if (openstatus > 0) then 
            ! Something wrong when reading file - continue with default
            ! values
            print *, 'ReadOperationOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadOperationOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if     
        
        ! Read options
        !=============
        field = 'ga.rem_small_trias'
        call ExtractOptionValueLogical0D(fid, field, options%rem_small_trias)        
        field = 'ga.cut_off_pol'
        call ExtractOptionValueReal0D(fid, field, options%cut_off_pol) 
        field = 'ga.cut_off_surf'
        call ExtractOptionValueReal0D(fid, field, options%cut_off_surf)
        field = 'ga.stacked_trias'
        call ExtractOptionValueLogical0D(fid, field, options%stacked_trias)          
        

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)        
        
    end subroutine

    subroutine ReadSplittingOptions(options)
        ! Description
        !============
        ! Read in user-specified adaption operation options 
        
        ! Declare variables
        !==================
        ! Arguments
        class(SplittingOptionsUDT)  :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: field
        integer, parameter              :: fid = 10        
        
        ! Initialize
        !===========        

    end subroutine

    subroutine ReadMergingOptions(options)
        ! Description
        !============
        ! Read in user-specified adaption operation options 
        
        ! Declare variables
        !==================
        ! Arguments
        class(MergingOptionsUDT)  :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: field
        integer, parameter              :: fid = 10        
        
        ! Initialize
        !===========        
    end subroutine

    subroutine ReadPentagonOptions(options)
        ! Description
        !============
        ! Read in user-specified adaption operation options 
        
        ! Declare variables
        !==================
        ! Arguments
        class(PentagonOptionsUDT)  :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: field
        integer, parameter              :: fid = 10        
        
        ! Initialize
        !===========        
    end subroutine

    subroutine ReadDistancefunctionOptions(options)
        ! Description
        !============
        ! Read in user-specified adaption operation options 
        
        ! Declare variables
        !==================
        ! Arguments
        class(DistancefunctionOptionsUDT)  :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: field
        integer, parameter              :: fid = 10        
        
        ! Initialize
        !===========        
    end subroutine

end module