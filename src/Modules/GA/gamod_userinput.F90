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
        logical                     :: rem_small_trias
        real(R8)                    :: cut_off_pol
        real(R8)                    :: cut_off_surf

        logical                     :: stacked_trias
        logical                     :: stacked_trias_checkAR
        real(R8)                    :: stacked_trias_maxAR
        logical                     :: merge_stacked_trias
        real(R8)                    :: merge_stacked_trias_angle_threshold
        logical                     :: merge_trap_into_stacked

        logical                     :: stacked_to_cutcell
        logical                     :: stacked_to_cutcell_uniform

        logical                     :: split_shaved_off_tube


        logical                     :: splitting
        logical                     :: merging
        logical                     :: pents_to_tria
        real(R8)                    :: h_rad_threshold
        real(R8)                    :: h_rad_core_threshold

        logical                     :: rem_stickout_trias
        logical                     :: rem_trias_flux
        real(R8)                    :: rem_tube_outershell_threshold
        character(:), allocatable   :: outershell_handling
        logical                     :: remove_stickoutquad
        logical                     :: split_noalignedquads


    contains

        procedure :: SetDefaults     => SetDefaultOperationOptions
        procedure :: Read           => ReadOperationOptions

    end type

    ! Options to specify the splitting operations
    type, extends(OptionsUDT) :: SplittingOptionsUDT

        ! TODO explanation on all options
        logical                     :: no_pents
        character(:), allocatable   :: QTtype
        logical                     :: split_out
        character(:), allocatable   :: splittype
        integer(I8)                 :: n_split
        character(:), allocatable   :: typeT
        character(:), allocatable   :: rad_type
        character(:), allocatable   :: pol_type
        real(R8)                    :: dist_function_threshold_split
        real(R8)                    :: dist_function_threshold_split_Wall
        
        
    contains
  
        procedure :: SetDefaults     => SetDefaultSplittingOptions
        procedure :: Read           => ReadSplittingOptions

    end type    

    ! Options to specify the merging operations
    type, extends(OptionsUDT) :: MergingOptionsUDT

        ! TODO explain all options
        logical                     :: no_hex
        character(:), allocatable   :: merge_crit
        real(R8)                    :: merge_h_pol_factor 
        integer(I8)                 :: n_merge
        real(R8)                    :: merge_bias_limit
        real(R8)                    :: dist_function_threshold_merge

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
        options%no_pents                            = .true.
        options%QTtype                              = 'regular'
        options%split_out                           = .false.
        options%splittype                           = 'rad'
        options%n_split                             = 20
        options%typeT                               = 'cutcell'
        options%rad_type                            = 'h_rad'
        options%pol_type                            = 'trias'
        options%dist_function_threshold_split       = 0.9
        options%dist_function_threshold_split_wall  = 0.6    
        
        
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
        options%no_hex                          = .true.
        options%merge_crit                      = 'h_pol'
        options%merge_h_pol_factor              = 1
        options%n_merge                          = 20
        options%merge_bias_limit                = 5
        options%dist_function_threshold_merge   = 0.6
        
        
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

        ! Small triangles
        !================
        field = 'ga.rem_small_trias'
        call ExtractOptionValueLogical0D(fid, field, options%rem_small_trias)        
        field = 'ga.cut_off_pol'
        call ExtractOptionValueReal0D(fid, field, options%cut_off_pol) 
        field = 'ga.cut_off_surf'
        call ExtractOptionValueReal0D(fid, field, options%cut_off_surf)

        ! Stacked triangles
        !==================
        field = 'ga.stacked_trias'
        call ExtractOptionValueLogical0D(fid, field, options%stacked_trias)          
        field = 'ga.stacked_trias_checkAR'
        call ExtractOptionValueLogical0D(fid, field, options%stacked_trias_checkAR)     
        field = 'ga.stacked_trias_maxAR'
        call ExtractOptionValueReal0D(fid, field, options%stacked_trias_maxAR)
        field = 'ga.merge_stacked_trias'
        call ExtractOptionValueLogical0D(fid, field, options%merge_stacked_trias)
        field = 'ga.merge_stacked_trias_angle_threshold'
        call ExtractOptionValueReal0D(fid, field, options%merge_stacked_trias_angle_threshold)
        field = 'ga.merge_trap_into_stacked'
        call ExtractOptionValueLogical0D(fid, field, options%merge_trap_into_stacked)    
        
        ! Stacked to cutcell
        !===================
        field = 'ga.stacked_to_cutcell'
        call ExtractOptionValueLogical0D(fid, field, options%stacked_to_cutcell)
        field = 'ga.stacked_to_cutcell_uniform'
        call ExtractOptionValueLogical0D(fid, field, options%stacked_to_cutcell_uniform)

        ! Splitting and merging
        !======================
        field = 'ga.splitting'
        call ExtractOptionValueLogical0D(fid, field, options%splitting)    
        field = 'ga.pents_to_tria' 
        call ExtractOptionValueLogical0D(fid, field, options%pents_to_tria)
        field = 'ga.merging'
        call ExtractOptionValueLogical0D(fid, field, options%merging)
        field = 'ga.h_rad_threshold'
        call ExtractOptionValueReal0D(fid, field, options%h_rad_threshold)                
        field = 'ga.h_rad_core_threshold'
        call ExtractOptionValueReal0D(fid, field, options%h_rad_core_threshold)  
        
        ! Special operations
        !===================
        field = 'ga.rem_stickout_trias'
        call ExtractOptionValueLogical0D(fid, field, options%rem_stickout_trias)        
        field = 'ga.rem_trias_flux'
        call ExtractOptionValueLogical0D(fid, field, options%rem_trias_flux)        
        field = 'ga.rem_tube_outershell_threshold'
        call ExtractOptionValueReal0D(fid, field, options%rem_tube_outershell_threshold)  
        field = 'ga.outershell_handling'
        call ExtractOptionValueCharacter(fid, field, options%outershell_handling)
        field = 'ga.remove_stickoutquad'        
        call ExtractOptionValueLogical0D(fid, field, options%remove_stickoutquad)
        field = 'ga.split_noalignedquads'
        call ExtractOptionValueLogical0D(fid, field, options%split_noalignedquads)        
        
        
        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)        
        
    end subroutine

    subroutine ReadSplittingOptions(options)
        ! Description
        !============
        ! Read in user-specified splitting options 
        
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

        ! Open the file, check if it exists
        open(unit=fid, file=options%inputfilepath, status='old', &
            iostat=openstatus)

        if (openstatus > 0) then 
            ! Something wrong when reading file - continue with default
            ! values
            print *, 'ReadSplittingOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadSplittingOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if     
        
        ! Read options
        !=============  
        field = 'ga.no_pents'
        call ExtractOptionValueLogical0D(fid, field, options%no_pents)         
        field = 'ga.QTtype'
        call ExtractOptionValueCharacter(fid, field, options%QTtype) 
        field = 'ga.split_out'
        call ExtractOptionValueLogical0D(fid, field, options%split_out)
        field = 'ga.splittype'
        call ExtractOptionValueCharacter(fid, field, options%splittype)
        field = 'ga.n_split'
        call ExtractOptionValueInteger0D(fid, field, options%n_split)
        field = 'ga.typeT'
        call ExtractOptionValueCharacter(fid, field, options%typeT)                                   
        field = 'ga.rad_type'
        call ExtractOptionValueCharacter(fid, field, options%rad_type)                                   
        field = 'ga.pol_type'
        call ExtractOptionValueCharacter(fid, field, options%pol_type)
        field = 'ga.dist_function_threshold_split'
        call ExtractOptionValueReal0D(fid, field, options%dist_function_threshold_split)                                    
        field = 'ga.dist_function_threshold_split_wall'
        call ExtractOptionValueReal0D(fid, field, options%dist_function_threshold_split_wall)                                           

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid) 


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

        ! Open the file, check if it exists
        open(unit=fid, file=options%inputfilepath, status='old', &
            iostat=openstatus)

        if (openstatus > 0) then 
            ! Something wrong when reading file - continue with default
            ! values
            print *, 'ReadMergingOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadMergingOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if     
        
        ! Read options
        !=============      
        field = 'ga.no_hex'
        call ExtractOptionValueLogical0D(fid, field, options%no_hex)
        field = 'ga.merge_crit'
        call ExtractOptionValueCharacter(fid, field, options%merge_crit) 
        field = 'ga.merge_h_pol_factor'
        call ExtractOptionValueReal0D(fid, field, options%merge_h_pol_factor)
        field = 'ga.n_merge'
        call ExtractOptionValueInteger0D(fid, field, options%n_merge)
        field = 'ga.merge_bias_limit'
        call ExtractOptionValueReal0D(fid, field, options%merge_bias_limit)        
        field = 'ga.dist_function_threshold_merge'
        call ExtractOptionValueReal0D(fid, field, options%dist_function_threshold_merge)        


        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)


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