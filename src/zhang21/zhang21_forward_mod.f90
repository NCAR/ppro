module zhang21_forward_mod
!--------------------------
! Zhang21 Polarimetric Radar Operator
!
! Parameterized polarimetric radar operator following
! Zhang, G., J. Gao, and M. Du, 2021: Parameterized forward operators 
!   for simulation and assimilation of polarimetric radar data with 
!   numerical weather predictions. Adv. Atmos. Sci., 38(5), 737−754.
!
! Authors:
!   Zhiquan (Jake) Liu (NCAR/MMM) - Initial implementation, July 2022
!   Hejun Xie (NCAR/MMM) - Integration into JEDI-UFO framework, TL/AD development (not testing yet), August 2024
!   Tao Sun (NCAR/MMM) - Adding hail categories and extended microphysics support, April 2025
!   Rong Kong (NCAR/MMM) - Bug fixes, operator tuning and testing; Modularization and optimization, October 2025
!--------------------------------------

  implicit none

  ! Physical constants
  real (kind=8), parameter :: PI = 3.141592653589793238462643383279502884197_8
  real (kind=8), parameter :: dpi = PI  ! Keep for backward compatibility

  ! Hydrometeor densities [g/cm^3]
  real (kind=8), target :: density_rain    = 1.0    ! Pure water
  real (kind=8), target :: density_snow    = 0.1    ! Dry snow
  real (kind=8), target :: density_graupel = 0.5    ! Dry graupel
  real (kind=8), target :: density_hail    = 0.917  ! Dry hail (ice)
                        
  ! Water fraction (ratio of melting, 0=dry, 1=wet)
  real (kind=8), target :: ratio_rain     = 1.0
  real (kind=8), target :: ratio_snow     = 0.0
  real (kind=8), target :: ratio_graupel  = 0.0
  real (kind=8), target :: ratio_hail     = 0.0

  ! Lookup table coefficients
  real (kind=8), dimension(0:12,0:3), target :: sband_snow_coefs, sband_graupel_coefs, sband_hail_coefs
  real (kind=8), dimension(0:12,0:3), target :: cband_snow_coefs, cband_graupel_coefs, cband_hail_coefs
  real (kind=8), dimension(1:4,0:4), target  :: sband_rain_coefs, cband_rain_coefs

  real (kind=8), dimension(0:12), target :: sband_snow_a, sband_graupel_a, sband_hail_a
  real (kind=8), dimension(0:12), target :: cband_snow_a, cband_graupel_a, cband_hail_a

  ! Initialization flag
  logical, save :: coefs_initialized = .false.
  character(len=512), save :: coef_path_saved = ''

  ! ============================================================================
  ! Melting skip threshold (HARD-CODED)
  ! ============================================================================
  ! Skip melting calculation if qr + qx < qsum_min_threshold
  ! This avoids computing melting for negligible mixing ratios
  ! ============================================================================
  real(kind=8), parameter :: qsum_min_threshold = 1.0d-3  ! g/kg (0.001 g/kg)
  
  ! ============================================================================
  ! Skip small qx (configurable via YAML)
  ! ============================================================================
  ! When enabled, skips all radar variable calculations for hydrometeor species
  ! with very small mixing ratios. Sets ZH/KDP/dm/w/z=0, ZDR=1, PhV=1.
  ! When disabled, calculates radar variables for all non-zero mixing ratios.
  ! ============================================================================
  logical, save :: skip_small_qx = .true.  ! Default: enabled
  public :: skip_small_qx
  
  ! ============================================================================
  ! Melting transition control (configurable via YAML)
  ! ============================================================================
  ! When enabled, uses smooth transition function instead of hard cutoff
  ! Transition is based on balance = min(qx/qr, qr/qx):
  !   - balance < ratio_low  -> factor = 0 (no melting)
  !   - balance > ratio_high -> factor = 1 (full melting)
  !   - in between           -> linear transition
  ! ============================================================================
  logical, save :: enable_melting_transition = .false.  ! Default disabled, can be overridden via YAML
  
  ! Transition parameters (configurable via YAML)
  real(kind=8), save :: snow_ratio_low    = 0.01d0    ! factor=0 below this
  real(kind=8), save :: snow_ratio_high   = 0.05d0    ! factor=1 above this
  real(kind=8), save :: graupel_ratio_low  = 0.01d0
  real(kind=8), save :: graupel_ratio_high = 0.05d0
  real(kind=8), save :: hail_ratio_low    = 0.01d0
  real(kind=8), save :: hail_ratio_high   = 0.05d0
  
  public :: enable_melting_transition
  public :: snow_ratio_low, snow_ratio_high
  public :: graupel_ratio_low, graupel_ratio_high
  public :: hail_ratio_low, hail_ratio_high

  ! ============================================================================
  ! Melting water fraction scaling factor (configurable via YAML)
  ! ============================================================================
  ! Scales the melting fraction (rats/ratg/rath) before mass partitioning:
  !   rats = min(melting_water_fraction * rats, 1.0)
  ! This affects both mass partitioning (qmsr = qms * rats) and downstream
  ! scattering coefficients (coef_a) and melting density.
  !   1.0 (default) = no scaling (original melting fraction)
  !   0.3           = reduce melting water fraction to 30%
  ! ============================================================================
  real(kind=8), save :: melting_water_fraction = 1.0d0
  
  public :: melting_water_fraction

  ! ============================================================================
  ! Dm-based melting limit (configurable via YAML)
  ! ============================================================================
  ! When enabled, prevents melting of large pure ice particles based on
  ! temperature-dependent Dm threshold. This is physically motivated:
  ! - At low temperatures (< -5°C), large ice particles (> 1 mm) are less likely to melt
  ! - At higher temperatures (> 0°C), larger particles (> 2 mm) can still melt
  ! - Uses smooth transition function to avoid discontinuities
  ! ============================================================================
  logical, save :: enable_dm_melting_limit = .false.  ! Default disabled
  real(kind=8), save :: dm_melting_transition_width = 0.5d0  ! mm, controls smoothness of transition
  
  public :: enable_dm_melting_limit
  public :: dm_melting_transition_width

  ! ============================================================================
  ! Dm regularization for small ntx/qx (configurable via YAML)
  ! ============================================================================
  ! When enabled, applies regularization in dm_z_2moment to handle small ntx/qx:
  !   - ntx < 1000 AND qx < 1e-5: dm = 0 (effectively no contribution)
  !   - qx > 1e-5 AND ntx < 2000: smooth transition using effective qx/ntx
  !   - otherwise: standard calculation
  ! This prevents extreme Dm values when number concentration is very small.
  ! When disabled (default): uses standard dm calculation for all ntx/qx values,
  ! which matches the original codeall behavior.
  ! ============================================================================
  logical, save :: enable_dm_regularization = .true.    ! Default enabled
  
  public :: enable_dm_regularization

  ! ============================================================================
  ! Dm tuning coefficients for each hydrometeor type (configurable via YAML)
  ! ============================================================================
  ! Multiplicative scaling factors applied to mean diameter (Dm) after computation:
  !   Dm_new = tuning_dm_xxx * Dm_old
  !   Z_new  = Z_old * tuning_dm_xxx^3   (since Z ∝ Dm^3)
  ! Default: 1.0 (no scaling). Values < 1 reduce Dm, > 1 increase it.
  ! Each coefficient affects all downstream radar variables (ZH, ZDR, KDP)
  ! for the corresponding hydrometeor species.
  ! ============================================================================
  real(kind=8), save :: tuning_dm_rain           = 1.0d0
  real(kind=8), save :: tuning_dm_melting_snow    = 1.0d0
  real(kind=8), save :: tuning_dm_melting_graupel = 1.0d0
  real(kind=8), save :: tuning_dm_melting_hail    = 1.0d0
  real(kind=8), save :: tuning_dm_pure_snow       = 1.0d0
  real(kind=8), save :: tuning_dm_pure_graupel    = 1.0d0
  real(kind=8), save :: tuning_dm_pure_hail       = 1.0d0

  public :: tuning_dm_rain
  public :: tuning_dm_melting_snow, tuning_dm_melting_graupel, tuning_dm_melting_hail
  public :: tuning_dm_pure_snow, tuning_dm_pure_graupel, tuning_dm_pure_hail

  ! ============================================================================
  ! Melting rain exponent (legacy, kept for backward compatibility)
  ! ============================================================================
  real(kind=8), save :: melting_rain_exponent = 0.5d0
  public :: melting_rain_exponent

  ! ============================================================================
  ! Melting power-law exponents (configurable via YAML, used in liu24)
  ! ============================================================================
  ! Controls the power-law formula for melting species:
  !   qmx  = qr^a  * qx^b  * coeff
  !   ntmx = ntr^a * ntx^b * coeff
  ! Default: a = 0.5, b = 0.5 => geometric mean sqrt(qr*qx)
  ! a and b are independent (do not need to sum to 1.0).
  ! ============================================================================
  real(kind=8), save :: melting_exp_rain = 0.5d0  ! exponent for rain (qr, ntr)
  real(kind=8), save :: melting_exp_ice  = 0.5d0  ! exponent for ice  (qx, ntx)
  public :: melting_exp_rain, melting_exp_ice

  ! ============================================================================
  ! Melting fraction exponents (configurable via YAML)
  ! ============================================================================
  ! Controls the power-law formula for melting fraction (water content ratio):
  !   rats = qr^a / (qr^a + qs^b)
  !   ratg = qr^a / (qr^a + qg^b)
  !   rath = qr^a / (qr^a + qh^b)
  ! Default: a = b = 1.0 => standard ratio qr/(qr+qx)
  ! ============================================================================
  real(kind=8), save :: fraction_exp_rain = 1.0d0  ! exponent for rain in fraction
  real(kind=8), save :: fraction_exp_ice  = 1.0d0  ! exponent for ice in fraction
  public :: fraction_exp_rain, fraction_exp_ice

  ! ============================================================================
  ! Melting fraction tuning coefficients (configurable via YAML)
  ! ============================================================================
  ! Multiplicative scaling factors applied to melting fraction
  ! after it is computed by the melting scheme:
  !   ratio_new = min(tuning_melt_frac_xxx * ratio_old, 1.0)
  ! Default: 1.0 (no scaling). Values < 1 reduce the liquid fraction (more ice-like),
  ! > 1 increase it (more rain-like). Clamped to [0, 1].
  ! Affects: melting density, scattering coefficients (coef_a), and mass partitioning.
  ! ============================================================================
  real(kind=8), save :: tuning_melt_frac_snow    = 1.0d0
  real(kind=8), save :: tuning_melt_frac_graupel = 1.0d0
  real(kind=8), save :: tuning_melt_frac_hail    = 1.0d0

  public :: tuning_melt_frac_snow, tuning_melt_frac_graupel, tuning_melt_frac_hail

  ! ============================================================================
  ! Melting model selection (configurable via YAML)
  ! ============================================================================
  ! Melting model based on Liu et al. 2024 with modifications by Kong:
  !   - Configurable power-law exponents (melting_exp_rain, melting_exp_ice)
  !   - Dm-based melting inhibition for large ice particles
  !   - Smooth melting transition function
  !   - Melting water fraction scaling
  !   - Melting coefficient tuning
  ! Default exponents: a = b = 0.5 (geometric mean)
  ! Exponents are configurable via melting_exp_rain and melting_exp_ice
  character(len=16), save :: melting_scheme_option = 'liu24'
  public :: melting_scheme_option

  ! ============================================================================
  ! dmmax configuration (configurable via YAML)
  ! ============================================================================
  ! Maximum mean diameter limits for each hydrometeor species [mm]
  ! Each species can be individually configured via YAML.
  ! If not specified, uses default values.
  !
  ! YAML keys:
  !   dmmax rain:            pure rain (default: 5.0)
  !   dmmax pure snow:       pure snow (default: 10.0)
  !   dmmax melting snow:    melting snow (default: 5.0)
  !   dmmax pure graupel:    pure graupel (default: 10.0)
  !   dmmax melting graupel: melting graupel (default: 5.0)
  !   dmmax pure hail:       pure hail (default: 10.0)
  !   dmmax melting hail:    melting hail (default: 5.0)
  ! ============================================================================
  real(kind=8), save :: dmmax_rain = 5.0d0            ! pr
  real(kind=8), save :: dmmax_pure_snow = 10.0d0      ! ps
  real(kind=8), save :: dmmax_melting_snow = 5.0d0    ! ms
  real(kind=8), save :: dmmax_pure_graupel = 10.0d0   ! pg
  real(kind=8), save :: dmmax_melting_graupel = 5.0d0 ! mg
  real(kind=8), save :: dmmax_pure_hail = 10.0d0      ! ph
  real(kind=8), save :: dmmax_melting_hail = 4.0d0    ! mh
  
  public :: dmmax_rain
  public :: dmmax_pure_snow, dmmax_melting_snow
  public :: dmmax_pure_graupel, dmmax_melting_graupel
  public :: dmmax_pure_hail, dmmax_melting_hail

  ! ============================================================================
  ! Treat hail as graupel option (configurable via YAML)
  ! ============================================================================
  ! When enabled, uses graupel coefficients and formula for hail calculations
  ! This can be useful for microphysics schemes that don't distinguish between
  ! hail and graupel well, or for testing purposes
  ! Default: false (uses normal hail coefficients and C-band special formula)
  ! ============================================================================
  logical, save :: treat_hail_as_graupel = .false.   ! Default disabled
  public :: treat_hail_as_graupel

  ! Public interfaces
  public :: zhang21_compute_point
  public :: zhang21_init_coefs
  public :: zhang21_finalize
  public :: melting_model_liu24
  
  ! Make module variables accessible to zhang21_tlad_mod
  ! These are needed by the TL/AD module
  public :: density_rain, density_snow, density_graupel, density_hail
  public :: ratio_rain, ratio_snow, ratio_graupel, ratio_hail
  public :: sband_snow_a, sband_graupel_a, sband_hail_a
  public :: cband_snow_a, cband_graupel_a, cband_hail_a
  public :: sband_rain_coefs, sband_snow_coefs, sband_graupel_coefs, sband_hail_coefs
  public :: cband_rain_coefs, cband_snow_coefs, cband_graupel_coefs, cband_hail_coefs
  
  ! Public utility functions (needed by zhang21_tlad_mod)
  public :: n0_lambda_wsm6
  public :: n0_lambda_gceop

contains

  ! ------------------------------------------------------------------------------
  subroutine read_coefs_icephase(filename, coefs)
  !------------------------------------------------------------------------
  ! Read coef C as defined in Table 1 (snow), 2 (hail) and 3 (graupel) 
  ! of Zhang et al., 2021, 13 x 4 array
  !------------------------------------------------------------------------
    implicit none
    character (len=*), intent(in)      :: filename
    real (kind=8),     intent(inout)   :: coefs(0:12,0:3)

    integer :: i, ios, unit_num
    logical :: file_exists

    ! Check if file exists
    inquire(file=trim(filename), exist=file_exists)
    if (.not. file_exists) then
      write(*,'(A)') '=========================================='
      write(*,'(A)') 'ERROR: PPRO coefficient file not found:'
      write(*,'(A)') '  ' // trim(filename)
      write(*,'(A)') 'Please ensure coefficient files are in:'
      write(*,'(A)') '  1. Current directory'
      write(*,'(A)') '  2. ./fix/ subdirectory'
      write(*,'(A)') '  3. Path specified by PPRO_COEF_PATH'
      write(*,'(A)') '=========================================='
      stop 1
    endif

    ! Open file with error checking
    open(newunit=unit_num, file=trim(filename), status='old', &
         form='formatted', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,I0)') 'ERROR: Cannot open file: ' // trim(filename) // ', iostat=', ios
      stop 1
    endif

    ! Read coefficients
    do i=0,12
       read(unit_num, *, iostat=ios) coefs(i,:)
       if (ios /= 0) then
         write(*,'(A,I0,A,I0)') 'ERROR: Reading line ', i, ' from ' // trim(filename), ios
         close(unit_num)
         stop 1
       endif
    end do
    close(unit_num)

  end subroutine read_coefs_icephase

  ! ------------------------------------------------------------------------------
  subroutine read_coefs_rain(filename, coefs)
  !------------------------------------------------------------------------
  ! Read coef in Eq. (13) - (16) of Zhang et al., 2021
  ! 4 x 5 array
  !------------------------------------------------------------------------
    implicit none
    character (len=*), intent(in)      :: filename
    real (kind=8),     intent(inout)   :: coefs(1:4,0:4)

    integer :: i, ios, unit_num
    logical :: file_exists

    ! Check if file exists
    inquire(file=trim(filename), exist=file_exists)
    if (.not. file_exists) then
      write(*,'(A)') '=========================================='
      write(*,'(A)') 'ERROR: PPRO coefficient file not found:'
      write(*,'(A)') '  ' // trim(filename)
      write(*,'(A)') '=========================================='
      stop 1
    endif

    ! Open file with error checking
    open(newunit=unit_num, file=trim(filename), status='old', &
         form='formatted', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,I0)') 'ERROR: Cannot open file: ' // trim(filename) // ', iostat=', ios
      stop 1
    endif

    ! Read coefficients
    do i=1,4
       read(unit_num, *, iostat=ios) coefs(i,:)
       if (ios /= 0) then
         write(*,'(A,I0,A)') 'ERROR: Reading line ', i, ' from ' // trim(filename)
         close(unit_num)
         stop 1
       endif
    end do
    close(unit_num)

  end subroutine read_coefs_rain

  ! ------------------------------------------------------------------------------
  !> @brief Public interface to initialize PPRO coefficients
  !>
  !> This should be called once before using ppro_compute_point.
  !> Subsequent calls will be skipped if already initialized.
  !>
  !> @param[in] coef_path Optional path to coefficient files directory
  subroutine zhang21_init_coefs(coef_path)
    implicit none
    character(len=*), intent(in), optional :: coef_path
    
    if (coefs_initialized) then
      write(*,'(A)') 'PPRO: Coefficients already initialized, skipping'
      return
    endif
    
    if (present(coef_path)) then
      coef_path_saved = trim(coef_path)
    else
      ! Try to get from environment variable
      call get_environment_variable('PPRO_COEF_PATH', coef_path_saved)
      if (len_trim(coef_path_saved) == 0) then
        coef_path_saved = '.'  ! Default to current directory
      endif
    endif
    
    write(*,'(A)') '=========================================='
    write(*,'(A)') 'PPRO: Initializing coefficients'
    write(*,'(A)') 'Coefficient path: ' // trim(coef_path_saved)
    write(*,'(A)') '=========================================='
    
    call init_coefs_internal(coef_path_saved)
    coefs_initialized = .true.
    
    write(*,'(A)') 'PPRO: Coefficients initialized successfully'
    
  end subroutine zhang21_init_coefs

  ! ------------------------------------------------------------------------------
  subroutine init_coefs_internal(base_path)
  !------------------------------------------------------------------------
  ! Internal subroutine to read coefficient files
  ! Read coef of rain/snow/graupel/hail of Zhang et al., 2021
  !------------------------------------------------------------------------
    implicit none
    character(len=*), intent(in) :: base_path
    character(len=512) :: filepath
    integer :: i

    ! Helper function to construct full path
    call build_coef_filepath(base_path, 'sband_rain_coefs.txt', filepath)
    call read_coefs_rain(filepath, sband_rain_coefs)
    
    call build_coef_filepath(base_path, 'sband_snow_coefs.txt', filepath)
    call read_coefs_icephase(filepath, sband_snow_coefs)
    
    call build_coef_filepath(base_path, 'sband_graupel_coefs.txt', filepath)
    call read_coefs_icephase(filepath, sband_graupel_coefs)
    
    call build_coef_filepath(base_path, 'sband_hail_coefs.txt', filepath)
    call read_coefs_icephase(filepath, sband_hail_coefs)

    ! Compute a coefficients for dry hydrometeors
    do i = 0, 12
      call coef_a(sband_snow_coefs(i,0:3), ratio_snow, sband_snow_a(i))
      call coef_a(sband_graupel_coefs(i,0:3), ratio_graupel, sband_graupel_a(i))
      call coef_a(sband_hail_coefs(i,0:3), ratio_hail, sband_hail_a(i))
    end do

    call build_coef_filepath(base_path, 'cband_rain_coefs.txt', filepath)
    call read_coefs_rain(filepath, cband_rain_coefs)
    
    call build_coef_filepath(base_path, 'cband_snow_coefs.txt', filepath)
    call read_coefs_icephase(filepath, cband_snow_coefs)
    
    call build_coef_filepath(base_path, 'cband_graupel_coefs.txt', filepath)
    call read_coefs_icephase(filepath, cband_graupel_coefs)
    
    call build_coef_filepath(base_path, 'cband_hail_coefs.txt', filepath)
    call read_coefs_icephase(filepath, cband_hail_coefs)

    do i = 0, 12
      call coef_a(cband_snow_coefs(i,0:3), ratio_snow, cband_snow_a(i))
      call coef_a(cband_graupel_coefs(i,0:3), ratio_graupel, cband_graupel_a(i))
      call coef_a(cband_hail_coefs(i,0:3), ratio_hail, cband_hail_a(i))
    end do

  end subroutine init_coefs_internal

  ! ------------------------------------------------------------------------------
  subroutine build_coef_filepath(base_path, filename, fullpath)
    !------------------------------------------------------------------------
    ! Build full path to coefficient file, trying multiple locations
    !------------------------------------------------------------------------
    implicit none
    character(len=*), intent(in) :: base_path, filename
    character(len=*), intent(out) :: fullpath
    logical :: file_exists
    
    ! Try 1: base_path/filename
    fullpath = trim(base_path) // '/' // trim(filename)
    inquire(file=trim(fullpath), exist=file_exists)
    if (file_exists) return
    
    ! Try 2: base_path/fix/zhang21_coefs/filename  
    fullpath = trim(base_path) // '/fix/zhang21_coefs/' // trim(filename)
    inquire(file=trim(fullpath), exist=file_exists)
    if (file_exists) return
    
    ! Try 3: base_path/fix/filename (backward compatibility)
    fullpath = trim(base_path) // '/fix/' // trim(filename)
    inquire(file=trim(fullpath), exist=file_exists)
    if (file_exists) return
    
    ! Try 4: Just filename (current directory)
    fullpath = trim(filename)
    inquire(file=trim(fullpath), exist=file_exists)
    if (file_exists) return
    
    ! Try 5: ./fix/zhang21_coefs/filename
    fullpath = './fix/zhang21_coefs/' // trim(filename)
    
  end subroutine build_coef_filepath

  ! ------------------------------------------------------------------------------
  !> @brief Finalize PPRO and reset initialization flag
  subroutine zhang21_finalize()
    implicit none
    coefs_initialized = .false.
    write(*,'(A)') 'PPRO: Finalized'
  end subroutine zhang21_finalize

  ! ------------------------------------------------------------------------------
  ! Backward compatibility wrapper
  subroutine init_coefs()
    call zhang21_init_coefs()
  end subroutine init_coefs

  ! ------------------------------------------------------------------------------
  subroutine dualpol_op_rain(w, dm, a, zh, zdr, kdp, phv)
  !--------------------------------------------------------
  ! Eqs. (13) - (16) of Zhang et al., 2021
  !------------------------------------------
    implicit none
    real(kind=8), intent(in)   :: w   ! density_air * mixing ratio , g/m^3
    real(kind=8), intent(in)   :: dm  ! mass/volume-weighted mean diameter, mm
    real(kind=8), intent(in)   :: a(1:4,0:4) ! coefs for rain
    real(kind=8), intent(out)  :: zh  ! horizontal reflectivity, mm^6/m^3
    real(kind=8), intent(out)  :: zdr ! differential reflectivity, -
    real(kind=8), intent(out)  :: kdp ! specific differential phase, degree/km^-1
    real(kind=8), intent(out)  :: phv ! co-polar correlation coefficient, 0-1

    zh  = w*(a(1,0) + a(1,1)*dm + a(1,2)*dm**2 + a(1,3)*dm**3 + a(1,4)*dm**4)**2

    zdr =  a(2,0) + a(2,1)*dm + a(2,2)*dm**2 + a(2,3)*dm**3 + a(2,4)*dm**4
    
    kdp = w*(a(3,0) + a(3,1)*dm + a(3,2)*dm**2 + a(3,3)*dm**3 + a(3,4)*dm**4)

    phv =    a(4,0) + a(4,1)*dm + a(4,2)*dm**2 + a(4,3)*dm**3 + a(4,4)*dm**4

  end subroutine dualpol_op_rain

  ! ------------------------------------------------------------------------------
  subroutine dualpol_op_icephase(z, w, rho, dm, a, zh, zdr, kdp, phv, iband, ptype)
  !-----------------------------------------------------------------
  ! Eqs. (17) - (20) of Zhang et al., 2021 for snow/graupel/hail
  !---------------------------------------------------------------
    implicit none
    integer, intent(in), optional   :: iband  
    character(len=*), intent(in), optional :: ptype
    real(kind=8), intent(in)   :: z   ! 6th moment of PSD, Eq. (6)
    real(kind=8), intent(in)   :: w   ! density_air * mixing ratio , g/m^3
    real(kind=8), intent(in)   :: rho ! density of melting snow/graupel/hail, g/cm^3
    real(kind=8), intent(in)   :: dm  ! mass/volume-weighted mean diameter, mm
    real(kind=8), intent(in)   :: a(0:12) ! coef

    real(kind=8), intent(out)  :: zh  ! horizontal reflectivity, mm^6/m^3
    real(kind=8), intent(out)  :: zdr ! differential reflectivity, -
    real(kind=8), intent(out)  :: kdp ! specific differential phase, degree/km^-1
    real(kind=8), intent(out)  :: phv ! co-polar correlation coefficient, 0-1

  ! --- Calculate reflectivity zh ---
    !C-band, the coefficient for C band can not fit the S-band formula well for hail
    if(present(iband) .and. iband ==2 .and. present(ptype) .and. &
       (ptype == "mh" .or. ptype == "ph")) then 
       !due to strong non-Rayleigh scattering effects for hail
       if(dm > 1.0D-3)then !in mm
          zh  = z*((a(0) + a(1)*dm + a(2)*dm**2 + a(3)*dm**3)/dm)**2
       else
          zh = 0.0d0
       endif
    else 
       zh  = z*(a(0) + a(1)*dm + a(2)*dm**2 + a(3)*dm**3)**2
    endif
    zdr =    a(4) + a(5)*dm + a(6)*dm**2

    if(rho > 1.0D-3)then
       kdp = w*(a(7) + a(8)*dm + a(9)*dm**2)/rho
    else
       kdp = 0.0d0
    endif
    phv =   a(10) + a(11)*dm + a(12)*dm**2

  end subroutine dualpol_op_icephase


  ! --------------------------------------------------------------
  subroutine dualpol_op_total ( zh_prain,  zh_psnow,  zh_pgraupel, &
                               zdr_prain, zdr_psnow, zdr_pgraupel, &
                               kdp_prain, kdp_psnow, kdp_pgraupel, &
                               phv_prain, phv_psnow, phv_pgraupel, &
                               zh, zdr, kdp, phv, &
                               zh_msnow,  zh_mgraupel, &                  !optional
                               zdr_msnow,  zdr_mgraupel, &                !optional
                               kdp_msnow,  kdp_mgraupel, &                !optional
                               phv_msnow,  phv_mgraupel, &                !optional
                               zh_mhail, zdr_mhail, kdp_mhail, phv_mhail, &   !optional
                               zh_phail, zdr_phail, kdp_phail, phv_phail, &
                               dmms, dmmg, dmmh, dm, temperature)    !optional
  !---------------------------------------------------------------
  ! compute total zh, zdr, kdp, phv
  ! Eqs. 22 - 25 of Zhang et al., 2021
  !---------------------------------------------------------------

    implicit none

    real(kind=8), intent(in)  :: zh_prain, zh_psnow, zh_pgraupel    ! horizontal reflectivity, mm^6/m^3
    real(kind=8), intent(in)  :: zdr_prain, zdr_psnow, zdr_pgraupel ! differential reflectivity, -
    real(kind=8), intent(in)  :: kdp_prain, kdp_psnow, kdp_pgraupel ! specific differential phase, degree/km^-1
    real(kind=8), intent(in)  :: phv_prain, phv_psnow, phv_pgraupel ! co-polar correlation coefficient, 0-1

    real(kind=8), intent(out)  :: zh  ! horizontal reflectivity, mm^6/m^3
    real(kind=8), intent(out)  :: zdr ! differential reflectivity, -
    real(kind=8), intent(out)  :: kdp ! specific differential phase, degree/km^-1
    real(kind=8), intent(out)  :: phv ! co-polar correlation coefficient, 0-1

    ! Melting snow and graupel (optional)
    real(kind=8), intent(in), optional   :: zh_msnow, zh_mgraupel    ! horizontal reflectivity, mm6 m-3
    real(kind=8), intent(in), optional   :: zdr_msnow, zdr_mgraupel   ! differential reflectivity, -
    real(kind=8), intent(in), optional   :: kdp_msnow, kdp_mgraupel ! specific differential phase, degree/km^-1
    real(kind=8), intent(in), optional   :: phv_msnow, phv_mgraupel ! co-polar correlation coefficient, 0-1

    ! Hail and pure hail (optional)
    real(kind=8), intent(in), optional :: zh_mhail, zh_phail    ! horizontal reflectivity, mm6 m-3
    real(kind=8), intent(in), optional :: zdr_phail              ! differential reflectivity, -
    real(kind=8), intent(in), optional :: zdr_mhail              ! differential reflectivity, -
    real(kind=8), intent(in), optional :: kdp_mhail, kdp_phail  ! specific differential phase, degree/km^-1
    real(kind=8), intent(in), optional :: phv_mhail, phv_phail  ! co-polar correlation coefficient, 0-1
    real(kind=8), intent(in), optional :: dmms, dmmg, dmmh
    real(kind=8), intent(in), optional :: dm(7)  ! Array of all mean diameters [pr,ms,mg,mh,ps,pg,ph]
    real(kind=8), intent(in), optional :: temperature  ! Temperature in Celsius

    real(kind=8) :: zbar_prain, zbar_msnow, zbar_mgraupel, zbar_mhail
    real(kind=8) :: zv_prain, zv_msnow, zv_mgraupel, zv_mhail
    real(kind=8) :: zbar_psnow, zbar_pgraupel, zbar_phail
    real(kind=8) :: zv_psnow, zv_pgraupel, zv_phail

    real(kind=8) :: zv, alpha
    real(kind=8) :: zh_thresh = 5.0
 
    logical :: include_hail = .false.
    logical :: pure_hydro = .false.

    if (present(zh_mhail) .and. present(zh_phail)) then
      include_hail = .true.
    endif

    ! Here the zdr is dimensionless (linear domain, before log10 conversion)
    ! For rain, zdr >= 1.0 is physically expected (oblate drops)
    ! When zdr < 1.0 (polynomial extrapolation artifact at small Dm),
    ! assume isotropic scattering (zdr = 1.0) as a safe fallback
    if(zdr_prain >= 1.0d0)then
       zv_prain = zh_prain / zdr_prain
       zbar_prain = zh_prain / sqrt(zdr_prain)
    else
       ! Isotropic assumption: zv = zh (equivalent to zdr = 1.0, i.e., 0 dB)
       zv_prain = zh_prain
       zbar_prain = zh_prain
    endif

    if (present(zh_msnow) .and. present(zdr_msnow)) then
      if(zdr_msnow >= 1.0d0)then
         zv_msnow = zh_msnow / zdr_msnow
         zbar_msnow = zh_msnow / sqrt(zdr_msnow)
      else
         ! Isotropic assumption: zdr < 1.0 is non-physical for melting snow
         zv_msnow = zh_msnow
         zbar_msnow = zh_msnow
      endif
    else
      zv_msnow = 0.0d0
      zbar_msnow = 0.0d0
    endif

    if (present(zh_mgraupel) .and. present(zdr_mgraupel)) then
      if(zdr_mgraupel >= 1.0d0)then
        zv_mgraupel = zh_mgraupel / zdr_mgraupel
        zbar_mgraupel = zh_mgraupel / sqrt(zdr_mgraupel)
      else
        ! Isotropic assumption: zdr < 1.0 is non-physical for melting graupel
        zv_mgraupel = zh_mgraupel
        zbar_mgraupel = zh_mgraupel
      endif
    else
      zv_mgraupel = 0.0d0
      zbar_mgraupel = 0.0d0
    endif

    if(zdr_psnow >= 1.0d0)then
       zv_psnow = zh_psnow / zdr_psnow
       zbar_psnow = zh_psnow / sqrt(zdr_psnow)
    else
       ! Isotropic assumption: zdr < 1.0 is non-physical for pure snow
       zv_psnow = zh_psnow
       zbar_psnow = zh_psnow
    endif

    if(zdr_pgraupel >= 1.0d0)then
       zv_pgraupel = zh_pgraupel / zdr_pgraupel
       zbar_pgraupel = zh_pgraupel / sqrt(zdr_pgraupel)
    else
       ! Isotropic assumption: zdr < 1.0 is non-physical for pure graupel
       zv_pgraupel = zh_pgraupel
       zbar_pgraupel = zh_pgraupel
    endif

    if (include_hail) then
       if (present(zh_mhail) .and. present(zdr_mhail)) then
         if(zdr_mhail >= 1.0d0)then
            zv_mhail = zh_mhail / zdr_mhail
            zbar_mhail = zh_mhail / sqrt(zdr_mhail)
         else
            ! Isotropic assumption: zdr < 1.0 is non-physical for melting hail
            zv_mhail = zh_mhail
            zbar_mhail = zh_mhail
         endif
       else
         zv_mhail = 0.0d0
         zbar_mhail = 0.0d0
       endif

       if (present(zh_phail) .and. present(zdr_phail)) then
         if(zdr_phail >= 1.0d0)then
            zv_phail = zh_phail / zdr_phail
            zbar_phail = zh_phail / sqrt(zdr_phail)
         else
            ! Isotropic assumption: zdr < 1.0 is non-physical for pure hail
            zv_phail = zh_phail
            zbar_phail = zh_phail
         endif
       else
         zv_phail = 0.0d0
         zbar_phail = 0.0d0
       endif
    endif

    if (include_hail) then
       if(pure_hydro)then
          zh  = zh_prain + zh_psnow + zh_pgraupel + zh_phail ! Eq. 22
          zv  = zv_prain + zv_psnow + zv_pgraupel + zv_phail
          if(zv > 1.0D-3)then
             zdr = zh/zv ! Eq. 23
          else
             zdr = 1.0d0  ! No signal: ZDR = 1 (0 dB)
          endif
          kdp = kdp_prain + kdp_psnow + kdp_pgraupel + kdp_phail ! Eq. 24
          if(zbar_prain + zbar_psnow + zbar_pgraupel + zbar_phail > 1.0D-3)then
             phv = (zbar_prain*phv_prain + &
                    zbar_psnow*phv_psnow + zbar_pgraupel*phv_pgraupel + zbar_phail*phv_phail) / &
                    (zbar_prain + zbar_psnow + zbar_pgraupel + zbar_phail) ! Eq. 25
          else
             phv = 1.0d0  ! No signal: PhV = 1 (perfect correlation)
          endif
       else
          zh  = zh_prain + zh_msnow + zh_mgraupel + zh_psnow + zh_pgraupel + zh_mhail + zh_phail ! Eq. 22
          zv  = zv_prain + zv_msnow + zv_mgraupel + zv_psnow + zv_pgraupel + zv_mhail + zv_phail
          if(zv > 1.0D-3)then
             zdr = zh/zv ! Eq. 23
          else
             zdr = 1.0d0  ! No signal: ZDR = 1 (0 dB)
          endif
          kdp = kdp_prain + kdp_msnow + kdp_mgraupel + kdp_psnow + kdp_pgraupel + kdp_mhail + kdp_phail ! Eq. 24
          if(zbar_prain + zbar_msnow + zbar_mgraupel + zbar_psnow + zbar_pgraupel + zbar_mhail + zbar_phail > 1.0D-3)then
             phv = (zbar_prain*phv_prain + zbar_msnow*phv_msnow + zbar_mgraupel*phv_mgraupel + zbar_mhail*phv_mhail + &
                    zbar_psnow*phv_psnow + zbar_pgraupel*phv_pgraupel + zbar_phail*phv_phail) / &
                    (zbar_prain + zbar_msnow + zbar_mgraupel + zbar_psnow + zbar_pgraupel + zbar_mhail + zbar_phail) ! Eq. 25
          else
             phv = 1.0d0  ! No signal: PhV = 1 (perfect correlation)
          endif    
       endif
    else
       ! No hail case
       zh  = zh_prain + zh_msnow + zh_mgraupel + zh_psnow + zh_pgraupel ! Eq. 22
       zv  = zv_prain + zv_msnow + zv_mgraupel + zv_psnow + zv_pgraupel

       if(zv > 1.0D-3)then
          zdr = zh/zv ! Eq. 23
       else
          zdr = 1.0d0  ! No signal: ZDR = 1 (0 dB)
       endif
      
       kdp = kdp_prain + kdp_msnow + kdp_mgraupel + kdp_psnow + kdp_pgraupel ! Eq. 24

       if(zbar_prain + zbar_msnow + zbar_mgraupel + zbar_psnow + zbar_pgraupel > 1.0D-3)then
          phv = (zbar_prain*phv_prain + zbar_msnow*phv_msnow + zbar_mgraupel*phv_mgraupel + &
                 zbar_psnow*phv_psnow + zbar_pgraupel*phv_pgraupel) / &
                 (zbar_prain + zbar_msnow + zbar_mgraupel + zbar_psnow + zbar_pgraupel) ! Eq. 25
       else
          phv = 1.0d0  ! No signal: PhV = 1 (perfect correlation)
       endif
    endif

    ! KDP non-negative protection: negative KDP is non-physical in normal precipitation
    kdp = max(kdp, 0.0d0)

  end subroutine dualpol_op_total

  ! ------------------------------------------------------------------------------
  subroutine coef_a(c, gamma, a)
  !---------------------------------------
  ! Eq. (21) of Zhang et al., 2021
  !---------------------------------------
    implicit none
    real(kind=8), intent(in)   :: c(0:3)  ! coef c in Table
    real(kind=8), intent(in)   :: gamma ! coef for Zdr
    real(kind=8), intent(out)  :: a  ! horizontal reflectivity, dBZ

    integer  :: i

    a = c(0)
    do i = 1, 3
      a = a + c(i) * gamma**i
    end do

  end subroutine coef_a

  !-------------------------------------------------------------------
  subroutine rainfrac(qr, qx, ratio)
  !------------------------------------------------------
    implicit none
    real(kind=8), intent(in)   :: qr    ! mixing ratio of rain, kg/kg
    real(kind=8), intent(in)   :: qx    ! mixing ratio of hydrometeor, kg/kg
    real(kind=8), intent(out)  :: ratio ! qr/(qr+qx)
    ! Calculate ratio directly, consistent with melting_model_liu24
    ! This function is kept for backward compatibility but is not used in melting model
    ! Melting model calculates ratio directly: rats = qr^a / (qr^a + qs^b)
    ! If qr + qx is very small, still calculate ratio correctly (unless both are zero)
    if(qr + qx > 1.0d-10)then  ! Changed from 1.0E-5 to 1.0d-10 to be more permissive
       ratio = qr/(qr+qx)
    else
      ! If sum is essentially zero, set ratio to 0 (no rain component)
      ! Note: This is consistent with melting_model_liu24 behavior when qr + qx < threshold
      ratio = 0.0d0
    endif

  end subroutine rainfrac

  !-------------------------------------------------------------------
  function get_dm_threshold_melting(temperature_celsius, species) result(dm_thresh)
  !-------------------------------------------------------------------
  ! Returns a temperature-dependent Dm threshold [mm] for melting 
  ! inhibition. Used by get_dm_melting_factor() to decide:
  !   Dm > dm_thresh  -->  melting is suppressed  (dm_factor -> 0)
  !   Dm < dm_thresh  -->  melting proceeds normally (dm_factor -> 1)
  !
  ! Temperature interpolation:
  !   T < thresh1 (-5C):  dm_thresh = dm_cold  (stricter, less melting)
  !   T > thresh2 ( 0C):  dm_thresh = dm_warm  (relaxed, more melting)
  !   thresh1 <= T <= thresh2:  linear interpolation between the two
  !
  ! Species-dependent values [mm]:
  !   Species   dm_cold  dm_warm   (denser ice melts slower -> larger threshold)
  !   Snow       1.0      2.0
  !   Graupel    2.0      4.0
  !   Hail       3.0      6.0
  !
  ! Fig 1: dm_thresh vs Temperature (snow example)
  !
  !   dm_thresh
  !    2.0 |     *-- dm_warm
  !        |    /
  !    1.5 |   /
  !        |  /  <- linear interpolation
  !    1.0 |-*  dm_cold
  !        |
  !        +--+--+-- T(C)
  !          -5  0
  !
  ! Fig 2: dm_factor vs Dm (sigmoid, for a given dm_thresh)
  !
  !   factor
  !    1.0 |----,
  !        |     \
  !    0.5 |      *  <- dm_thresh (midpoint)
  !        |       \
  !    0.0 |        '----
  !        +------+---- Dm
  !               ^
  !           dm_thresh
  !
  !   Dm < dm_thresh: factor ~ 1.0 (melting proceeds)
  !   Dm = dm_thresh: factor = 0.5 (transition midpoint)
  !   Dm > dm_thresh: factor ~ 0.0 (melting suppressed)
  !   Transition width controlled by dm_melting_transition_width
  !-------------------------------------------------------------------
    implicit none
    real(kind=8), intent(in) :: temperature_celsius
    character(len=*), intent(in), optional :: species
    real(kind=8) :: dm_thresh
    real(kind=8) :: dm_cold, dm_warm
    real(kind=8), parameter :: thresh1 = -5.0d0  ! Cold-end temperature [C]
    real(kind=8), parameter :: thresh2 =  0.0d0  ! Warm-end temperature [C]
    
    ! Species-dependent Dm thresholds at cold/warm ends
    ! hrad to melt, when dm > dm_thresh
    if (present(species)) then
      select case (trim(species))
        case ('snow')
          dm_cold = 1.0d0;  dm_warm = 2.0d0   ! Fluffy, melts fast
        case ('graupel')
          !dm_cold = 2.0d0;  dm_warm = 4.0d0   ! Medium density
          dm_cold = 1.5d0;  dm_warm = 3.0d0   ! Medium density
        case ('hail')
          dm_cold = 2.0d0;  dm_warm = 4.0d0   ! Dense, melts slow
        case default
          dm_cold = 1.0d0;  dm_warm = 2.0d0   ! Default: same as snow
      end select
    else
      dm_cold = 1.0d0;  dm_warm = 2.0d0       ! Default: same as snow
    endif
    
    if (temperature_celsius < thresh1) then
      dm_thresh = dm_cold
    else if (temperature_celsius > thresh2) then
      dm_thresh = dm_warm
    else
      ! Linear interpolation: thresh1 -> dm_cold, thresh2 -> dm_warm
      dm_thresh = dm_cold + (temperature_celsius - thresh1) / (thresh2 - thresh1) * (dm_warm - dm_cold)
    endif
    
  end function get_dm_threshold_melting

  !-------------------------------------------------------------------
  function get_dm_melting_factor(dm, dm_threshold, transition_width) result(factor)
  !-------------------------------------------------------------------
  ! Calculate smooth melting inhibition factor based on Dm
  !
  ! Uses sigmoid function to provide smooth transition:
  !   factor = 1 / (1 + exp((dm - dm_threshold) / transition_width))
  !
  ! Behavior:
  !   - When dm << dm_threshold: factor -> 1.0 (normal melting)
  !   - When dm >> dm_threshold: factor -> 0.0 (no melting)
  !   - Smooth transition around dm_threshold
  !
  ! Parameters:
  !   dm: Mean diameter of pure ice particle [mm]
  !   dm_threshold: Temperature-dependent Dm threshold [mm]
  !   transition_width: Width of transition region [mm]
  !
  ! Returns:
  !   factor: Melting inhibition factor [0-1]
  !-------------------------------------------------------------------
    implicit none
    real(kind=8), intent(in) :: dm, dm_threshold, transition_width
    real(kind=8) :: factor
    
    ! Sigmoid function for smooth transition
    ! factor = 1 / (1 + exp((dm - dm_threshold) / transition_width))
    if (transition_width > 0.0d0) then
      factor = 1.0d0 / (1.0d0 + exp((dm - dm_threshold) / transition_width))
    else
      ! If transition_width is zero or negative, use hard cutoff
      if (dm < dm_threshold) then
        factor = 1.0d0
      else
        factor = 0.0d0
      endif
    endif
    
  end function get_dm_melting_factor

  !-------------------------------------------------------------------
  subroutine weticephase_density(precip_type, ratio, rho)
  !------------------------------------------------------
  ! Eq. (7) of Zhang et al., 2021
  !--------------------------------
    implicit none
    character (len=*), intent(in) :: precip_type
    real(kind=8), intent(in)   :: ratio  ! qr/(qr+qx)
    real(kind=8), intent(out)  :: rho    ! density of melting precip, g/cm^3

    select case (precip_type)
       case ('snow')
          rho =  density_snow * (1.0 - ratio**2) + ratio**2
       case ('graupel')
          rho =  density_graupel * (1.0 - ratio**2) + ratio**2
       case ('hail')
          rho =  density_hail * (1.0 - ratio**2) + ratio**2
       case default
          write(*,*) 'Unknown precipitation type: ', precip_type
          stop
    end select

  end subroutine weticephase_density

  !-------------------------------------------------------------------
  subroutine watercontent(air_density, qx, w)
  !------------------------------------------------------
    implicit none
    real(kind=8), intent(in)   :: air_density   ! Air Density, 1.225 kg/m^3 
    real(kind=8), intent(in)   :: qx            ! mixing ratio of precip, g/kg
    real(kind=8), intent(out)  :: w             ! water content, g/m^3
    !print*,'air_density=',air_density
    w = air_density * qx

  end subroutine watercontent

  !------------------------------------------------------------------------------
  subroutine dm_z_2moment(air_density, qx, rhox, ntx, dm, z, dmmax, dm_min)
  !-----------------------------------------------------------
  ! Eqs. (5) and (6) of Zhang et al., 2021
  ! Assume 2-moment MP scheme with ntx and qx as prognostic variables
  ! Exponential PSD: N(D) = N0 * exp(-lambda D)
  !
  ! Dm = 4 * ( (rho_air * qx) / (pi * rhox * ntx) )^(1/3)
  ! Z  = 11250 * (rho_air * qx / (pi * rhox)) * Dm^3
  !
  ! Optional bounds:
  !   dmmax: upper bound on Dm (caps unrealistically large diameters)
  !   dm_min: lower bound on Dm; if Dm < dm_min, set Dm=0 and Z=0.
  !           Use 0.3 mm for ice species (snow/graupel/hail) since
  !           sub-0.3mm particles are cloud-ice scale. Do NOT pass
  !           for rain (drizzle Dm ~ 0.1-0.3 mm is physical).
  !
  ! Regularization (enable_dm_regularization = .true.):
  !   When both qx and ntx are very small (qx < 1e-3 g/kg AND
  !   ntx < 1000 /m^3), the particle population is physically
  !   meaningless and dm is set to 0. Otherwise, standard formula
  !   is used (caller ensures ntx >= ntx_min=10 for zero protection).
  !
  ! Without regularization:
  !   Standard formula with basic ntx > 0 protection.
  !-------------------------------------------------------------------
    implicit none
    real(kind=8), intent(in)   :: air_density   ! Air Density, 1.225 kg/m^3 
    real(kind=8), intent(in)   :: rhox          ! density of wet hydrometeor,g/cm^3 
                                                ! rain density = 1 g/cm^3
    real(kind=8), intent(in)   :: qx            ! mixing ratio of hydrometeor, g/kg
    real(kind=8), intent(in)   :: ntx           ! number concentration, /m^3
    real(kind=8), intent(in), optional   :: dmmax         ! maximum mean diameter, mm
    real(kind=8), intent(in), optional   :: dm_min        ! minimum meaningful Dm, mm
    real(kind=8), intent(out)  :: dm            ! mean diameter, mm
    real(kind=8), intent(out)  :: z             ! 6th moment of PSD

    if (enable_dm_regularization) then
      ! When both qx and ntx are small, particle population is negligible
      ! qx < 1e-3 g/kg: very small mixing ratio
      ! ntx < 1000 /m^3 (= 1/L): statistically unreliable sample
      if (qx < 1.0d-3 .and. ntx < 1000.0d0) then
        dm = 0.0d0
      else
        ! Standard calculation (caller guarantees ntx >= ntx_min=10)
        dm = 4.0d0 * ((air_density * qx) / (dpi * rhox * ntx)) ** (1.0d0 / 3.0d0)
        dm = dm * 10.0d0 ! cm -> mm
      endif
    else
      ! No regularization: standard calculation with basic protection
      if (ntx < 1.0d0) then
        ! Prevent division by zero / near-zero
        dm = 0.0d0
      else
        dm = 4.0d0 * ((air_density * qx) / (dpi * rhox * ntx)) ** (1.0d0 / 3.0d0)
        dm = dm * 10.0d0 ! cm -> mm
      endif
    endif

    if (present(dmmax)) then
      dm = min(dm, dmmax) ! physical upper bound on particle diameter
    endif

    ! Optional lower bound: Dm below dm_min has no physical meaning
    ! (e.g., 0.3 mm for ice species = cloud-ice scale, not precipitation)
    if (present(dm_min)) then
      if (dm < dm_min) dm = 0.0d0
    endif

    z = 11250.0d0 * ((air_density * qx) / (dpi * rhox)) * dm**3

  end subroutine dm_z_2moment

   !---------------------------------------------------------------------------
  subroutine melting_model_liu24( &
    qr, qs, qg, &
    qms, qmg, qpr, qps, qpg, rats, ratg, &
    ntr, nts, ntg, &
    ntms, ntmg, ntpr, ntps, ntpg, &
    qh, nth, rath, qmh, qph, ntmh, ntph, qx_min, ntx_min, coeff_melt, &
    density_air, temperature )
  !---------------------------------------------------------------------------
  ! Melting model based on Liu et al., 2024 with modifications by Kong (NCAR/MMM)
  !
  ! Original Liu24 formula: qmx = sqrt(qr * qx), i.e. geometric mean
  ! Modifications:
  !   - Configurable power-law exponents via melting_exp_rain / melting_exp_ice
  !     (default 0.5/0.5 = geometric mean)
  !   - Temperature-dependent Dm-based melting inhibition for large ice particles
  !   - Smooth melting transition function (balance ratio based)
  !   - Melting water fraction scaling (rats/ratg/rath)
  !   - Tunable melting coefficient (coeff_melt)
  !
  ! Configurable power-law formula:
  ! qmx  = qr^melting_exp_rain  * qx^melting_exp_ice  * coeff  (default a=b=0.5 => sqrt)
  ! ntmx = ntr^melting_exp_rain * ntx^melting_exp_ice * coeff
  ! Eq. (15-18) of Zhang et al., 2024
  ! Assume 2-moment MP scheme with ntx and qx as prognostic variables
  ! Exponential PSD: N(D) = N0 * exp(-lambda D)
  !
  ! Optional parameters:
  !   density_air: Air density [kg/m^3], needed for Dm calculation
  !   temperature: Air temperature [°C], needed for Dm-based melting limit
  !
  ! Special treatment and QC added by Rong Kong (NCAR/MMM):
  !   - Temperature-dependent Dm-based melting inhibition: Large pure ice
  !     particles (Dm > threshold) are less likely to melt, with threshold
  !     depending on temperature (1.0 mm at T < -5°C, 2.0 mm at T > 0°C).
  !   - Uses smooth sigmoid transition function to avoid discontinuities.
  !   - Quality control: Only calculates Dm when mixing ratios are non-negligible
  !     (qs/qg/qh > 1.0d-6) to avoid numerical issues.
  !   - Important: Without the qsum_min_threshold and Dm-based threshold
  !     checks, melting species (qms, qmg, qmh) would be overestimated,
  !     leading to unrealistic melting of small or large particles.
  !---------------------------------------------------------------------------
    implicit none
    real(kind=8), intent(in)   :: qr            ! mixing ratio of raw rain, g/kg
    real(kind=8), intent(in)   :: qs            ! mixing ratio of raw snow, g/kg
    real(kind=8), intent(in)   :: qg            ! mixing ratio of raw graupel, g/kg
    real(kind=8), intent(out)  :: qms           ! mixing ratio of melting snow, g/kg
    real(kind=8), intent(out)  :: qmg           ! mixing ratio of melting graupel, g/kg
    real(kind=8), intent(out)  :: qpr           ! mixing ratio of pure rain, g/kg
    real(kind=8), intent(out)  :: qps           ! mixing ratio of pure snow, g/kg
    real(kind=8), intent(out)  :: qpg           ! mixing ratio of pure graupel, g/kg
    real(kind=8), intent(out)  :: rats          ! mass water fraction of melting snow, [0-1]
    real(kind=8), intent(out)  :: ratg          ! mass water fraction of melting graupel, [0-1]
    real(kind=8), intent(in), optional   :: ntr           ! number concentration of raw rain, -/m3
    real(kind=8), intent(in), optional   :: nts           ! number concentration of raw snow, -/m3
    real(kind=8), intent(in), optional   :: ntg           ! number concentration of raw graupel, -/m3
    real(kind=8), intent(out), optional  :: ntms          ! number concentration of melting snow, -/m3
    real(kind=8), intent(out), optional  :: ntmg          ! number concentration of melting graupel, -/m3
    real(kind=8), intent(out), optional  :: ntpr          ! number concentration of pure rain, -/m3
    real(kind=8), intent(out), optional  :: ntps          ! number concentration of pure snow, -/m3
    real(kind=8), intent(out), optional  :: ntpg          ! number concentration of pure graupel, -/m3
    ! For hail
    real(kind=8), intent(in), optional   :: qh            ! mixing ratio of raw hail, g/kg
    real(kind=8), intent(in), optional   :: nth           ! number concentration of raw hail, -/m3
    real(kind=8), intent(out), optional  :: rath          ! mass water fraction of melting hail, [0-1]
    real(kind=8), intent(out), optional  :: qmh           ! mixing ratio of melting hail, g/kg
    real(kind=8), intent(out), optional  :: qph           ! mixing ratio of pure hail, g/kg
    real(kind=8), intent(out), optional  :: ntmh          ! number concentration of melting hail, -/m3
    real(kind=8), intent(out), optional  :: ntph          ! number concentration of pure hail, -/m3
    real(kind=8), intent(in), optional  ::  qx_min(7), ntx_min(7)
    real(kind=8), intent(in), optional  ::  coeff_melt    ! tuning coefficient for melting
    real(kind=8), intent(in), optional  ::  density_air   ! Air density [kg/m^3], for Dm calculation
    real(kind=8), intent(in), optional  ::  temperature   ! Air temperature [°C], for Dm-based melting limit


    real(kind=8) :: qmsr, qmgr, qmsd, qmgd, qmhr, qmhd
    real(kind=8) :: melt_coeff
    real(kind=8) :: balance, factor_s, factor_g, factor_h  ! For melting transition
    real(kind=8) :: eff_coeff_s, eff_coeff_g, eff_coeff_h  ! Effective melt coefficient = melt_coeff * factor
    logical      :: double_moment ! Flag for 2-moment microphysics scheme
    
    ! Variables for Dm-based melting limit
    real(kind=8) :: dm_pure_snow, dm_pure_graupel, dm_pure_hail
    real(kind=8) :: dm_thresh_snow, dm_thresh_graupel, dm_thresh_hail
    real(kind=8) :: dm_factor_s, dm_factor_g, dm_factor_h
    real(kind=8) :: zs_temp, zg_temp, zh_temp  ! Temporary variables for Dm calculation

   ! Check if double moment (2-moment scheme has number concentrations)
    double_moment = .false.
    if (present(ntr) .or. present(nts) .or. present(ntg) .or. present(nth)) then
        double_moment = .true.
    endif
    
    ! Initialize factors to 1.0 (no transition adjustment by default)
    factor_s = 1.0d0
    factor_g = 1.0d0
    factor_h = 1.0d0
    
    ! Initialize effective coefficients (will be set in melting calculation)
    eff_coeff_s = 0.0d0
    eff_coeff_g = 0.0d0
    eff_coeff_h = 0.0d0
    ! Note: WSM6 (1-moment) doesn't pass any number concentrations, double_moment=.false. is valid
    
    ! Initialize melting coefficient (default 0.3, or from parameter)
    if (present(coeff_melt)) then
      melt_coeff = coeff_melt
    else
      melt_coeff = 0.3d0
    endif
    
    ! Initialize Dm-based melting factors (default to 1.0, no inhibition)
    ! Special treatment added by Rong Kong (NCAR/MMM): Temperature-dependent
    ! Dm-based melting inhibition to prevent large pure ice particles from melting
    dm_factor_s = 1.0d0
    dm_factor_g = 1.0d0
    dm_factor_h = 1.0d0
    
    ! Calculate Dm-based melting inhibition factors if enabled
    ! QC: Only calculate when enable_dm_melting_limit is true and both
    ! temperature and density_air are provided
    IF (enable_dm_melting_limit .and. present(temperature) .and. present(density_air)) THEN
      ! Calculate Dm threshold based on temperature (species-dependent)
      dm_thresh_snow = get_dm_threshold_melting(temperature, 'snow')
      dm_thresh_graupel = get_dm_threshold_melting(temperature, 'graupel')
      dm_thresh_hail = get_dm_threshold_melting(temperature, 'hail')
      
      ! Calculate pure ice Dm before melting (for snow)
      ! QC: Only calculate when qs is non-negligible to avoid numerical issues
      if (qs > 1.0d-6) then
        if (present(nts) .and. nts > 1.0d-6) then
          ! 2-moment: calculate Dm from qs and nts
          call dm_z_2moment(density_air, qs, density_snow, nts, dm_pure_snow, zs_temp, dmmax=dmmax_pure_snow)
        else
          ! 1-moment: use WSM6 method (temperature is guaranteed to be present due to outer IF condition)
          call dm_z_wsm6('snow', density_air, temperature + 273.15d0, qs, density_snow, dm_pure_snow, zs_temp, dmmax_pure_snow)
        endif
        
        if (dm_pure_snow >= 0.3d0) then
          dm_factor_s = get_dm_melting_factor(dm_pure_snow, dm_thresh_snow, dm_melting_transition_width)
        else
          ! Dm < 0.3 mm: sub-cloud-particle scale, not physically meaningful
          ! as a melting species -> suppress melting entirely
          dm_factor_s = 0.0d0
        endif
      endif
      
      ! Calculate pure ice Dm before melting (for graupel)
      ! QC: Only calculate when qg is non-negligible to avoid numerical issues
      if (qg > 1.0d-6) then
        if (present(ntg) .and. ntg > 1.0d-6) then
          ! 2-moment: calculate Dm from qg and ntg
          call dm_z_2moment(density_air, qg, density_graupel, ntg, dm_pure_graupel, zg_temp, dmmax=dmmax_pure_graupel)
        else
          ! 1-moment: use WSM6 method (temperature is guaranteed to be present due to outer IF condition)
          call dm_z_wsm6('graupel', density_air, temperature + 273.15d0, qg, density_graupel, dm_pure_graupel, zg_temp, dmmax_pure_graupel)
        endif
        
        if (dm_pure_graupel >= 0.3d0) then
          dm_factor_g = get_dm_melting_factor(dm_pure_graupel, dm_thresh_graupel, dm_melting_transition_width)
        else
          dm_factor_g = 0.0d0
        endif
      endif
      
      ! Calculate pure ice Dm before melting (for hail)
      ! QC: Only calculate when qh is present and non-negligible to avoid numerical issues
      ! Same treatment as graupel: use standard threshold (no special qh-based adjustment)
      if (present(qh) .and. qh > 1.0d-6) then
        if (present(nth) .and. nth > 1.0d-6) then
          ! 2-moment: calculate Dm from qh and nth
          call dm_z_2moment(density_air, qh, density_hail, nth, dm_pure_hail, zh_temp, dmmax=dmmax_pure_hail)
        else
          ! 1-moment: use WSM6 method (temperature is guaranteed to be present due to outer IF condition)
          call dm_z_wsm6('hail', density_air, temperature + 273.15d0, qh, density_hail, dm_pure_hail, zh_temp, dmmax_pure_hail)
        endif
        
        if (dm_pure_hail >= 0.3d0) then
          dm_factor_h = get_dm_melting_factor(dm_pure_hail, dm_thresh_hail, dm_melting_transition_width)
        else
          dm_factor_h = 0.0d0
        endif
      endif
    ENDIF
    
    ! Compute melting mixing ratios with tuning coefficient
    ! When enable_melting_transition is true, use smooth transition function
    ! based on balance = min(qx/qr, qr/qx):
    !   balance < ratio_low  -> factor = 0 (no melting)
    !   balance > ratio_high -> factor = 1 (full melting)
    !   in between           -> linear transition
    
    ! =========================================================================
    ! Snow melting
    ! =========================================================================
    ! QC: Skip melting if qr + qs is below qsum_min_threshold to avoid
    ! overestimation of melting species from negligible mixing ratios
    IF (qr + qs < qsum_min_threshold) THEN
      ! Skip melting if qr + qs is essentially zero (avoid division by zero)
      qms = 0.0d0
      rats = 0.0d0
    ELSE
      ! Calculate transition factor first (if enabled)
      IF (enable_melting_transition) THEN
        balance = min(qs/qr, qr/qs)
        IF (balance < snow_ratio_low) THEN
          factor_s = 0.0d0
        ELSEIF (balance > snow_ratio_high) THEN
          factor_s = 1.0d0
        ELSEIF (snow_ratio_high <= snow_ratio_low) THEN
          factor_s = 1.0d0
        ELSE
          factor_s = (balance - snow_ratio_low) / (snow_ratio_high - snow_ratio_low)
        ENDIF
      ENDIF
      
      ! Apply factor to melt_coeff for consistent qms and ntms calculation
      ! Also apply Dm-based melting inhibition factor to prevent overestimation
      ! of melting species for large particles (Dm > threshold)
      eff_coeff_s = melt_coeff * factor_s * dm_factor_s
      qms = (qr**melting_exp_rain) * (qs**melting_exp_ice) * eff_coeff_s
      rats = qr**fraction_exp_rain / (qr**fraction_exp_rain + qs**fraction_exp_ice)
      qms = max(qms, qx_min(2))
    ENDIF

    ! =========================================================================
    ! Graupel melting
    ! =========================================================================
    ! QC: Skip melting if qr + qg is below qsum_min_threshold to avoid
    ! overestimation of melting species from negligible mixing ratios
    IF (qr + qg < qsum_min_threshold) THEN
      ! Skip melting if qr + qg is essentially zero (avoid division by zero)
      qmg = 0.0d0
      ratg = 0.0d0
    ELSE
      ! Calculate transition factor first (if enabled)
      IF (enable_melting_transition) THEN
        balance = min(qg/qr, qr/qg)
        IF (balance < graupel_ratio_low) THEN
          factor_g = 0.0d0
        ELSEIF (balance > graupel_ratio_high) THEN
          factor_g = 1.0d0
        ELSEIF (graupel_ratio_high <= graupel_ratio_low) THEN
          factor_g = 1.0d0
        ELSE
          factor_g = (balance - graupel_ratio_low) / (graupel_ratio_high - graupel_ratio_low)
        ENDIF
      ENDIF
      
      ! Apply factor to melt_coeff for consistent qmg and ntmg calculation
      ! Also apply Dm-based melting inhibition factor to prevent overestimation
      ! of melting species for large particles (Dm > threshold)
      eff_coeff_g = melt_coeff * factor_g * dm_factor_g
      qmg = (qr**melting_exp_rain) * (qg**melting_exp_ice) * eff_coeff_g
      ratg = qr**fraction_exp_rain / (qr**fraction_exp_rain + qg**fraction_exp_ice)
      qmg = max(qmg, qx_min(3))
    ENDIF

    ! =========================================================================
    ! Hail melting
    ! =========================================================================
    if (present(qh)) then
      ! QC: Skip melting if qr + qh is below qsum_min_threshold to avoid
      ! overestimation of melting species from negligible mixing ratios
      IF (qr + qh < qsum_min_threshold) THEN
        ! Skip melting if qr + qh is essentially zero (avoid division by zero)
        qmh = 0.0d0
        rath = 0.0d0
      ELSE
        ! Calculate transition factor first (if enabled)
        IF (enable_melting_transition) THEN
          balance = min(qh/qr, qr/qh)
          IF (balance < hail_ratio_low) THEN
            factor_h = 0.0d0
          ELSEIF (balance > hail_ratio_high) THEN
            factor_h = 1.0d0
          ELSEIF (hail_ratio_high <= hail_ratio_low) THEN
            factor_h = 1.0d0
          ELSE
            factor_h = (balance - hail_ratio_low) / (hail_ratio_high - hail_ratio_low)
          ENDIF
        ENDIF
        
        ! Apply factor to melt_coeff for consistent qmh and ntmh calculation
        ! Also apply Dm-based melting inhibition factor to prevent overestimation
        ! of melting species for large particles (Dm > threshold)
        eff_coeff_h = melt_coeff * factor_h * dm_factor_h
        qmh = (qr**melting_exp_rain) * (qh**melting_exp_ice) * eff_coeff_h
        rath = qr**fraction_exp_rain / (qr**fraction_exp_rain + qh**fraction_exp_ice)
        qmh = max(qmh, qx_min(4))
      ENDIF
    endif

    ! Scale melting fractions by melting_water_fraction (default 1.0 = no scaling)
    ! This controls both mass partitioning (qmsr) and downstream scattering coefficients (coef_a)
    rats = min(melting_water_fraction * rats, 1.0d0)
    ratg = min(melting_water_fraction * ratg, 1.0d0)
    if (present(qh)) rath = min(melting_water_fraction * rath, 1.0d0)

    qmsr = qms * rats
    qmgr = qmg * ratg
   
    if (present(qh)) then
       qmhr = qmh * rath
    endif

    ! Dry component: qmsd = qms - qmsr (mass conservation)
    qmsd = qms - qmsr
    qmgd = qmg - qmgr
    if (present(qh)) qmhd = qmh - qmhr

    qpr = qr - qmsr - qmgr
    qps = qs - qmsd
    qpg = qg - qmgd
    if (present(qh)) then
       qpr = qpr - qmhr
       qph = qh - qmhd
    endif

    qpr = max(qpr, qx_min(1))
    qps = max(qps, qx_min(5))
    qpg = max(qpg, qx_min(6))
    if (present(qh)) qph = max(qph, qx_min(7))



    if (double_moment) then
      ! Number concentration for melting species (same exponents a, b as qmx)
      if(present(nts) .and. present(ntr) .and. present(ntms)) then
         if (qms > 1.0d-3) then
           ntms = (ntr**melting_exp_rain) * (nts**melting_exp_ice) * eff_coeff_s
           ntms = max(ntms, ntx_min(2))
         else
           ntms = ntx_min(2)
         endif
      endif

      if(present(ntg) .and. present(ntr) .and. present(ntmg)) then
         if (qmg > 1.0d-3) then
           ntmg = (ntr**melting_exp_rain) * (ntg**melting_exp_ice) * eff_coeff_g
           ntmg = max(ntmg, ntx_min(3))
         else
           ntmg = ntx_min(3)
         endif
      endif

      if (present(nth) .and. present(ntr) .and. present(ntmh)) then
         if (present(qh) .and. qmh > 1.0d-3) then
           ntmh = (ntr**melting_exp_rain) * (nth**melting_exp_ice) * eff_coeff_h
           ntmh = max(ntmh, ntx_min(4))
         else
           ntmh = ntx_min(4)
         endif
      endif

      if(present(ntr)) ntpr = max(ntr,ntx_min(1))
      if(present(nts)) ntps = max(nts,ntx_min(5))
      if(present(ntg)) ntpg = max(ntg,ntx_min(6))
      if(present(nth)) ntph = max(nth,ntx_min(7))

    endif

  end subroutine melting_model_liu24

  !---------------------------------------------------------------------------
  subroutine melting_model_jung08( &
    qr, qs, qg, &
    qms, qmg, qpr, qps, qpg, rats, ratg, &
    ntr, nts, ntg, &
    ntms, ntmg, ntpr, ntps, ntpg, &
    qh, nth, rath, qmh, qph, ntmh, ntph)
  !---------------------------------------------------------------------------
  ! Classic melting scheme based on Jung et al., 2008
  ! Uses min-ratio criterion and power-law formula:
  !   qmx = (qr + qx) * 0.5 * minratio^0.3
  ! where minratio = min(qr/qx, qx/qr), threshold > 0.01
  !---------------------------------------------------------------------------
    implicit none
    real(kind=8), intent(in)   :: qr            ! mixing ratio of raw rain, g/kg
    real(kind=8), intent(in)   :: qs            ! mixing ratio of raw snow, g/kg
    real(kind=8), intent(in)   :: qg            ! mixing ratio of raw graupel, g/kg
    real(kind=8), intent(out)  :: qms           ! mixing ratio of melting snow, g/kg
    real(kind=8), intent(out)  :: qmg           ! mixing ratio of melting graupel, g/kg
    real(kind=8), intent(out)  :: qpr           ! mixing ratio of pure rain, g/kg
    real(kind=8), intent(out)  :: qps           ! mixing ratio of pure snow, g/kg
    real(kind=8), intent(out)  :: qpg           ! mixing ratio of pure graupel, g/kg
    real(kind=8), intent(out)  :: rats          ! mass water fraction of melting snow, [0-1]
    real(kind=8), intent(out)  :: ratg          ! mass water fraction of melting graupel, [0-1]
    real(kind=8), intent(in), optional   :: ntr           ! number concentration of raw rain, -/m3
    real(kind=8), intent(in), optional   :: nts           ! number concentration of raw snow, -/m3
    real(kind=8), intent(in), optional   :: ntg           ! number concentration of raw graupel, -/m3
    real(kind=8), intent(out), optional  :: ntms          ! number concentration of melting snow, -/m3
    real(kind=8), intent(out), optional  :: ntmg          ! number concentration of melting graupel, -/m3
    real(kind=8), intent(out), optional  :: ntpr          ! number concentration of pure rain, -/m3
    real(kind=8), intent(out), optional  :: ntps          ! number concentration of pure snow, -/m3
    real(kind=8), intent(out), optional  :: ntpg          ! number concentration of pure graupel, -/m3
    ! For hail
    real(kind=8), intent(in), optional   :: qh            ! mixing ratio of raw hail, g/kg
    real(kind=8), intent(in), optional   :: nth           ! number concentration of raw hail, -/m3
    real(kind=8), intent(out), optional  :: rath          ! mass water fraction of melting hail, [0-1]
    real(kind=8), intent(out), optional  :: qmh           ! mixing ratio of melting hail, g/kg
    real(kind=8), intent(out), optional  :: qph           ! mixing ratio of pure hail, g/kg
    real(kind=8), intent(out), optional  :: ntmh          ! number concentration of melting hail, -/m3
    real(kind=8), intent(out), optional  :: ntph          ! number concentration of pure hail, -/m3

    real(kind=8) :: qmsr, qmgr, qmsd, qmgd, qmhr, qmhd
    logical      :: double_moment ! Flag for 2-moment microphysics scheme
    real(kind=8) :: minratio

    double_moment = .false.
    if (present(ntr) .and. present(nts) .and. present(ntg)) then
      if (present(ntms) .and. present(ntmg) .and. present(ntpr) .and. present(ntps) .and. present(ntpg)) then
        double_moment = .true.
      else
        print *, "Please specify output for number concentrations of the melting scheme"
        stop
      endif
    endif

    IF (qs .GT. 0.0 .and. qr .GT. 0.0) then
      minratio = MIN(qr/qs, qs/qr)
      IF (minratio .GT. 0.01) then 
      qms = (qr + qs) * 0.5 * minratio ** 0.3
      rats = qr / (qr + qs)
      ELSE
        qms = 0.0
        rats = 0.0
      ENDIF
    ELSE
      qms = 0.0
      rats = 0.0
    ENDIF

    IF (qg .GT. 0.0 .and. qr .GT. 0.0) then
      minratio = MIN(qr/qg, qg/qr)
      IF (minratio .GT. 0.01) then 
        qmg = (qr + qg) * 0.5 * minratio ** 0.3
        ratg = qr / (qr + qg)
      ELSE
        qmg = 0.0
        ratg = 0.0
      ENDIF
    ELSE
      qmg = 0.0
      ratg = 0.0
    ENDIF

    if (present(qh)) then
       if (qh .gt. 0. .and. qr .gt. 0.) then
          minratio = MIN(qr/qh, qh/qr)
          if (minratio .GT. 0.01) then
             qmh = (qr + qh) * 0.5 * minratio ** 0.3
             rath = qr / (qr + qh)
          else
             qmh = 0.0
             rath = 0.0
          endif
       endif
    endif

    qmsr = qms * rats
    qmgr = qmg * ratg
    qmsd = qms * (1.0 - rats)
    qmgd = qmg * (1.0 - ratg)

    if (present(qh)) then
       qmhr = qmh * rath
       qmhd = qmh * (1.0 - rath)
    endif

    qpr = qr - qmsr - qmgr
    qps = qs - qmsd
    qpg = qg - qmgd

    if (present(qh)) then
       qpr = qpr - qmhr
       qph = qh - qmhd
    endif

    if (double_moment) then
      ntms = sqrt(ntr * nts)
      ntmg = sqrt(ntr * ntg)
      ntpr = ntr
      ntps = nts
      ntpg = ntg
      if (present(nth)) then
         ntmh = sqrt(ntr * nth)
         ntph = nth
      endif 
    endif
  
  end subroutine melting_model_jung08

  subroutine n0_lambda_wsm6(precip_type, air_density, temp, qx, rhox, n0, lambda)
   !------------------------------------------------------------------------
   ! Single moment scheme assuming Exponential PSD: N(D) = N0 * exp(-lambda D)
   ! N0 is fixed except for snow, then calcuate lambda from N0
   !------------------------------------------------------------------------
     implicit none
     character (len=*), intent(in) :: precip_type  ! rain/snow/graupel
     real(kind=8), intent(in)   :: air_density   ! Air Density, 1.225 kg/m^3 
     real(kind=8), intent(in)   :: rhox          ! density of wet hydrometeor,g/cm^3 
                                                 ! rain density = 1 g/cm^3
     real(kind=8), intent(in)   :: temp          ! temperature, Kelvin
     real(kind=8), intent(in)   :: qx            ! mixing ratio of hydrometeor, g/kg
     real(kind=8), intent(out)  :: n0            ! intercept parameter of exponential PSD, m-4
     real(kind=8), intent(out)  :: lambda        ! slope parameter of exponential PSD, m-1
 
     select case (precip_type)
        case ('rain')
            n0 = 8.0e6 ! m^(-4)
        case ('snow', 'pure snow')
            n0 = 2.0e6*exp(0.12*(273.15-temp))   !Eq (6) in Hong et al. 2004
        case ('graupel', 'pure graupel','hail','pure hail')
            n0 = 4.0e6 ! m^(-4)
        case default
           write(*,*) 'Unknow precipitation type: ', precip_type
           stop
     end select
 
       !rhox = rhox * 1000.0 ! g/cm^3 = 10^-3 kg/ (10-2 m)^3 = 1000 kg/m^3
       !  qx = qx * 0.001   ! g/kg -> kg/kg
     lambda = ((dpi*1000.0*rhox*n0)/(air_density*qx*0.001))**0.25 ! m^-1, Year2-report
   
   end subroutine n0_lambda_wsm6

  !------------------------------------------------------------------------------
  subroutine dm_z_wsm6(precip_type, air_density, temp, qx, rhox, dm, z, dmmax, dm_min)
  !------------------------------------------------------------------------
  ! Single moment scheme assuming Exponential PSD: N(D) = N0 * exp(-lambda D)
  ! N0 is fixed except for snow, then calcuate lambda from N0
  ! Optional dm_min: lower bound on Dm (use 0.3 for ice, omit for rain)
  !-------------------------------------------------------------------
    implicit none
    character (len=*), intent(in) :: precip_type  ! rain/snow/graupel
    real(kind=8), intent(in)   :: air_density   ! Air Density, 1.225 kg/m^3 
    real(kind=8), intent(in)   :: rhox          ! density of wet hydrometeor,g/cm^3 
                                                ! rain density = 1 g/cm^3
    real(kind=8), intent(in)   :: temp          ! temperature, Kelvin
    real(kind=8), intent(in)   :: qx            ! mixing ratio of hydrometeor, g/kg
    real(kind=8), intent(in),optional   :: dmmax         ! maximum mean diameter
    real(kind=8), intent(in),optional   :: dm_min        ! minimum meaningful Dm, mm
    real(kind=8), intent(out)  :: dm            ! mean diameter, mm
    real(kind=8), intent(out)  :: z             ! 6th moment of PSD

    real(kind=8)  :: n0  ! intercept parameter of exponential PSD
    real(kind=8)  :: lambda ! slope parameter of exponential PSD

    call n0_lambda_wsm6(precip_type, air_density, temp, qx, rhox, n0, lambda)
   
    dm = 4000.0 / lambda  ! Eq. (5), should be in mm, *1000?

    if(present(dmmax))then
       dm = min(dm, dmmax) !to avoid extremely large diameter due to very small ntx
    endif

    ! Optional lower bound for ice species
    if (present(dm_min)) then
      if (dm < dm_min) dm = 0.0d0
    endif

    z = 11250.0 * ( (air_density*qx)/(dpi*rhox) )* dm**3 ! Eq. (6)

    if(qx < 1.0D-3)then
      dm = 0.0
      z = 0.0
    endif

  end subroutine dm_z_wsm6

  subroutine n0_lambda_gceop(precip_type, air_density, qx, rhox, n0, lambda)
    !------------------------------------------------------------------------
    ! Single moment scheme assuming Exponential PSD: N(D) = N0 * exp(-lambda D)
    ! N0 is fixed, calcuate lambda from N0
    ! Goddard Cumulus Ensemble scheme used in CWB WRF operation
    ! only difference from wsm6: n0 for snow is a fixed value
    !------------------------------------------------------------------------
      implicit none
      character (len=*), intent(in) :: precip_type  ! rain/snow/graupel
      real(kind=8), intent(in)   :: air_density   ! Air Density, 1.225 kg/m^3 
      real(kind=8), intent(in)   :: rhox          ! density of wet hydrometeor,g/cm^3 
                                                  ! rain density = 1 g/cm^3
      real(kind=8), intent(in)   :: qx            ! mixing ratio of hydrometeor, g/kg
      real(kind=8), intent(out)  :: n0            ! intercept parameter of exponential PSD, m-4
      real(kind=8), intent(out)  :: lambda        ! slope parameter of exponential PSD, m-1
  
      select case (precip_type)
          case ('rain')
              n0 = 8.0e6 ! m^(-4)
          case ('snow', 'pure snow')
              n0 = 1.6e7 ! m^(-4)
          case ('graupel', 'pure graupel')
              n0 = 4.0e6 ! m^(-4)
          case default
            write(*,*) 'Unknow precipitation type: ', precip_type
            stop
      end select
  
        !rhox = rhox * 1000.0 ! g/cm^3 = 10^-3 kg/ (10-2 m)^3 = 1000 kg/m^3
        !  qx = qx * 0.001   ! g/kg -> kg/kg
      lambda = ((dpi*1000.0*rhox*n0)/(air_density*qx*0.001))**0.25 ! m^-1, Year2-report
      
  end subroutine n0_lambda_gceop

  !------------------------------------------------------------------------------
  subroutine dm_z_gceop(precip_type, air_density, qx, rhox, dm, z)
  !------------------------------------------------------------------------
  ! Single moment scheme assuming Exponential PSD: N(D) = N0 * exp(-lambda D)
  ! N0 is fixed, calcuate lambda from N0
  ! Goddard Cumulus Ensemble scheme used in CWB WRF operation
  ! only difference from wsm6: n0 for snow is a fixed value
  !-------------------------------------------------------------------
    implicit none
    character (len=*), intent(in) :: precip_type  ! rain/snow/graupel
    real(kind=8), intent(in)   :: air_density   ! Air Density, 1.225 kg/m^3 
    real(kind=8), intent(in)   :: rhox          ! density of wet hydrometeor,g/cm^3 
                                                ! rain density = 1 g/cm^3
    real(kind=8), intent(in)   :: qx            ! mixing ratio of hydrometeor, g/kg
    real(kind=8), intent(out)  :: dm            ! mean diameter, mm
    real(kind=8), intent(out)  :: z             ! 6th moment of PSD

    real(kind=8)  :: n0  ! intercept parameter of exponential PSD
    real(kind=8)  :: lambda ! slope parameter of exponential PSD

    call n0_lambda_gceop(precip_type, air_density, qx, rhox, n0, lambda)

    dm = 4000.0 / lambda  ! Eq. (5), should be in mm, *1000?
      z = 11250.0 * ( (air_density*qx)/(dpi*rhox) )* dm**3 ! Eq. (6)

  end subroutine dm_z_gceop

  subroutine dm_den_gcecwb(precip_type, logw, logt, frac, rhox, dm)
    implicit none
    character (len=*), intent(in) :: precip_type  ! rain/snow/graupel
    real(kind=8), intent(in)      :: logw         ! ln(water content), kg/m^3
    real(kind=8), intent(in)      :: logt         ! ln(temperature)
    real(kind=8), intent(in)      :: frac         ! melting fraction=qr/(qr+qx)
    real(kind=8), intent(out)     :: rhox         ! mean diameter, mm
    real(kind=8), intent(out)     :: dm           ! mean diameter, mm

    ! Diagnose density
    select case (precip_type)
       case ('rain')
          rhox = 1000.0 ! kg/m^3
            dm = exp(-20.189916 + 6.888453*logt + 0.25*logw  &
                         - 0.3132*logt**2 + 1.43e-18*logw**2  &
                         - 5.27932e-16*logt*logw) ! in micron meter
       case ('snow', 'pure snow')
          rhox = exp(65.84696 - 17.124*logt + 0.84975*logw  &
                         + 1.069132*logt**2 - 0.00614*logw**2 &
                         - 0.208331*logt*logw) ! kg/m^3
            dm = exp(-87.705062 + 26.74956*logt - 1.327406*logw  &
                         - 1.6701*logt**2 + 0.009598*logw**2  &
                         + 0.32543632*logt*logw) ! in micron meter
       case ('graupel', 'pure graupel')
          rhox = 1000.0*frac + &  ! kg/m^3
                 (1-frac)*max(250.0,exp(0.056433*logw+6.608048))
            dm = exp(-74.804847 + 16.08519*logt - 2.530305*logw  &
                         - 0.19094*logt**2 + 0.005846*logw**2  &
                         + 0.51917156*logt*logw) ! in micron meter
       case default
          write(*,*) 'Unknown precipitation type: ', precip_type
          stop
    end select

    dm = dm * 0.001   ! um -> mm
    rhox = rhox * 0.001 ! kg/m^3 -> g/cm^3

  end subroutine dm_den_gcecwb

  !------------------------------------------------------------------------------
  subroutine dm_z_gcecwb(precip_type, air_density, temp, qx, frac, dm, z, rhox)
  !-------------------------------------------------------------------------------
  ! Single moment scheme assuming Gamma PSD: N(D) = N0 * D^alpha * exp(-lambda D)
  ! Goddard Cumulus Ensemble scheme modified by CWB WRF
  ! hydrometeor particle density and mean diameter are diagnosed.
  !
  ! T.-C. Tsai, S.-Y. Jiang, J.-P. Chen, H.-L. Huang, Yi.-C. Lo, J.-S. Hong, 2022:
  !  A Gamma-type Single Moment Bulk Microphysics Scheme for Operational Forecasting.
  !----------------------------------------------------------------------------------
    implicit none
    character (len=*), intent(in) :: precip_type  ! rain/snow/graupel
    real(kind=8), intent(in)   :: air_density   ! Air Density, 1.225 kg/m^3 
    real(kind=8), intent(in)   :: temp          ! temperature, Kelvin
    real(kind=8), intent(in)   :: qx            ! mixing ratio of hydrometeor, g/kg
    real(kind=8), intent(in)   :: frac          ! melting fraction=qr/(qr+qx)
    real(kind=8), intent(out)  :: dm            ! mean diameter, mm
    real(kind=8), intent(out)  :: z             ! 6th moment of PSD
    real(kind=8), intent(out)  :: rhox          ! density of melting snow/graupel

    real(kind=8)  :: logw  ! ln(water content), kg/m^3
    real(kind=8)  :: logt  ! ln(temperature)

    
    logw = log(air_density * (0.001*qx)) ! kg/m^3
    logt = log(temp)  ! temp in Kelvin

    call dm_den_gcecwb(precip_type, logw, logt, frac, rhox, dm)

    ! air_density(kg/m^3), qx(g/kg), rhox(g/cm^3), dm(mm)  
    z = 11250.0 * ( (air_density*qx)/(dpi*rhox) ) * dm**3 ! Eq. (6)

  end subroutine dm_z_gcecwb

function lin_interp(h, h1, h2, v1, v2) result(vout)
  real(kind=8), intent(in) :: h, h1, h2, v1, v2
  real(kind=8) :: vout
  if (h <= h1) then
    vout = v1
  elseif (h >= h2) then
    vout = v2
  else
    vout = v1 + (v2 - v1) * (h - h1) / (h2 - h1)
  endif
end function lin_interp

  ! ==============================================================================
  !> @brief High-level interface: Compute dual-pol variables for one observation
  !>
  !> This is the main API for ppro library. It encapsulates all physical computations
  !> including melting layer treatment and aggregation of all hydrometeor types.
  !>
  !> @param[in]  iband         Radar band: 1=S-band, 2=C-band
  !> @param[in]  scheme_type   Microphysics scheme: 'WSM6', 'THOMPSON', 'NSSL', 'TCWA2'
  !> @param[in]  density_air   Dry air density [kg/m^3]
  !> @param[in]  temp_air      Air temperature [K]
  !> @param[in]  qr,qs,qg      Rain/snow/graupel mixing ratio [g/kg]
  !> @param[out] zh,zdr,kdp,phv  Dual-pol variables (linear units)
  !> @param[in]  qh,nr,ns,ng,nh,vg,vh  Optional for different schemes
  !>
  !> @note set_trajectory functionality moved to separate TL/AD interface
  ! ==============================================================================
  subroutine zhang21_compute_point(iband, scheme_type, density_air, temp_air, &
                                qr, qs, qg, zh, zdr, kdp, phv, &
                                qh, nr, ns, ng, nh, vg, vh, &
                                height, zhobs, zdrobs, kdpobs, coeff_melt, temperature, iobs)
    implicit none
    
    ! Inputs
    integer, intent(in) :: iband
    character(len=*), intent(in) :: scheme_type
    real(kind=8), intent(in) :: density_air, temp_air
    real(kind=8), intent(in) :: qr, qs, qg
    
    ! Outputs
    real(kind=8), intent(out) :: zh, zdr, kdp, phv
    
    ! Optional inputs
    real(kind=8), intent(in), optional :: qh, nr, ns, ng, nh, vg, vh
    real(kind=8), intent(in), optional :: height, zhobs, zdrobs, kdpobs
    real(kind=8), intent(in), optional :: coeff_melt, temperature
    integer, intent(in), optional :: iobs
    
    ! Local variables
    real(kind=8) :: qrreg, qsreg, qgreg, qhreg
    real(kind=8) :: nrreg, nsreg, ngreg, nhreg
    real(kind=8) :: qpr, qms, qmg, qps, qpg, qmh, qph
    real(kind=8) :: ntpr, ntms, ntmg, ntps, ntpg, ntmh, ntph
    real(kind=8) :: zhpr, zdrpr, kdppr, phvpr
    real(kind=8) :: zhms, zdrms, kdpms, phvms
    real(kind=8) :: zhmg, zdrmg, kdpmg, phvmg
    real(kind=8) :: zhps, zdrps, kdpps, phvps
    real(kind=8) :: zhpg, zdrpg, kdppg, phvpg
    real(kind=8) :: zhmh, zdrmh, kdpmh, phvmh
    real(kind=8) :: zhph, zdrph, kdpph, phvph
    real(kind=8) :: wpr, wms, wmg, wps, wpg, wmh, wph
    real(kind=8) :: zpr, zms, zmg, zps, zpg, zmh, zph
    real(kind=8) :: dmpr, dmms, dmmg, dmps, dmpg, dmmh, dmph
    real(kind=8) :: denms, denmg, denmh, denpg, denph
    real(kind=8) :: rats, ratg, rath
    real(kind=8) :: as(0:12), ag(0:12), ah(0:12)
    real(kind=8), dimension(:,:), pointer :: rain_coefs => NULL()
    real(kind=8), dimension(:,:), pointer :: snow_coefs => NULL()
    real(kind=8), dimension(:,:), pointer :: graupel_coefs => NULL()
    real(kind=8), dimension(:,:), pointer :: hail_coefs => NULL()
    real(kind=8), dimension(:), pointer :: snow_a => NULL()
    real(kind=8), dimension(:), pointer :: graupel_a => NULL()
    real(kind=8), dimension(:), pointer :: hail_a => NULL()
    integer :: i
    real(kind=8) :: qx_min(7), ntx_min(7)
    real(kind=8) :: dmmax(7)
    real(kind=8) :: dm(7)  ! Array of mean diameters: [pr, ms, mg, mh, ps, pg, ph]
    real(kind=8) :: temperature_celsius  ! Temperature in Celsius for melting scheme
    
    ! Convert temperature from Kelvin to Celsius if provided
    if (present(temperature)) then
      temperature_celsius = temperature  ! Assume temperature is already in Celsius
    else
      ! If temperature not provided, convert from temp_air (Kelvin) to Celsius
      temperature_celsius = temp_air - 273.15d0
    endif
    
    ! Initialize parameters
    qx_min = 0.0d0
    ntx_min = 10.0d0  ! Same as code_develop: nx_min = 10.0
    ! Use individual dmmax values (configurable via YAML, with defaults)
    ! Index: 1=pr, 2=ms, 3=mg, 4=mh, 5=ps, 6=pg, 7=ph
    dmmax(1) = dmmax_rain
    dmmax(2) = dmmax_melting_snow
    dmmax(3) = dmmax_melting_graupel
    dmmax(4) = dmmax_melting_hail
    dmmax(5) = dmmax_pure_snow
    dmmax(6) = dmmax_pure_graupel
    dmmax(7) = dmmax_pure_hail
    !----------------------------------------------------------------
    !1_pr; 2_ms; 3_mg; 4_mh; 5_ps; 6_pg; 7_ph 
    !----------------------------------------------------------------

    ! Select coefficients based on radar band
    if (iband == 1) then
      rain_coefs => sband_rain_coefs
      snow_coefs => sband_snow_coefs
      graupel_coefs => sband_graupel_coefs
      hail_coefs => sband_hail_coefs
      snow_a => sband_snow_a
      graupel_a => sband_graupel_a
      hail_a => sband_hail_a
    else
      rain_coefs => cband_rain_coefs
      snow_coefs => cband_snow_coefs
      graupel_coefs => cband_graupel_coefs
      hail_coefs => cband_hail_coefs
      snow_a => cband_snow_a
      graupel_a => cband_graupel_a
      hail_a => cband_hail_a
    endif
    
    ! Regularize inputs
    qrreg = max(qx_min(1), qr)
    qsreg = max(qx_min(5), qs)
    qgreg = max(qx_min(6), qg)
    if (present(qh)) qhreg = max(qx_min(7), qh)
    if (present(nr)) nrreg = max(ntx_min(1), nr)
    if (present(ns)) nsreg = max(ntx_min(5), ns)
    if (present(ng)) ngreg = max(ntx_min(6), ng)
    if (present(nh)) nhreg = max(ntx_min(7), nh)
    
    ! Call melting scheme based on microphysics type
    ! ============================================================================
    ! Currently available in MPAS:
    !   - WSM6      : 1-moment (no number concentrations)
    !   - Thompson  : Partial 2-moment rain only (nr)
    !   - TCWA2     : Full 2-moment without hail (nr,ns,ng)
    !   - NSSL      : Full 2-moment with hail (nr,ns,ng,nh,qh)
    !
    ! Placeholder for future MPAS availability:
    !   - Lin, WSM5, etc.    : 1-moment (same as WSM6 mode)
    !   - Morrison           : Full 2-moment without hail (same as TCWA2 mode)
    !   - Milbrandt-Yau, MY2 : Full 2-moment with hail (same as NSSL mode)
    !   - (nr,ns,ng,qh)      : 2-moment + 1-moment hail (placeholder)
    !   - (nr,ns)            : Partial 2-moment rain+snow (placeholder)
    !   - (nr,qh)            : Partial 2-moment rain + hail (placeholder)
    !   - (qh only)          : 1-moment with hail (placeholder)
    ! ============================================================================
    
    ! Call melting model (dispatch based on melting_scheme_option)
    select case (trim(melting_scheme_option))
    case ('jung08')
      ! Jung et al. 2008 classic melting scheme
      if (present(nr) .and. present(ns) .and. present(ng) .and. present(nh) .and. present(qh)) then
        call melting_model_jung08(qrreg, qsreg, qgreg, qms, qmg, qpr, qps, qpg, rats, ratg, &
             ntr=nrreg, nts=nsreg, ntg=ngreg, ntms=ntms, ntmg=ntmg, ntpr=ntpr, ntps=ntps, ntpg=ntpg, &
             qh=qhreg, nth=nhreg, rath=rath, qmh=qmh, qph=qph, ntmh=ntmh, ntph=ntph)
      elseif (present(nr) .and. present(ns) .and. present(ng) .and. present(qh)) then
        call melting_model_jung08(qrreg, qsreg, qgreg, qms, qmg, qpr, qps, qpg, rats, ratg, &
             ntr=nrreg, nts=nsreg, ntg=ngreg, ntms=ntms, ntmg=ntmg, ntpr=ntpr, ntps=ntps, ntpg=ntpg, &
             qh=qhreg, rath=rath, qmh=qmh, qph=qph)
      elseif (present(nr) .and. present(ns) .and. present(ng)) then
        call melting_model_jung08(qrreg, qsreg, qgreg, qms, qmg, qpr, qps, qpg, rats, ratg, &
             ntr=nrreg, nts=nsreg, ntg=ngreg, ntms=ntms, ntmg=ntmg, ntpr=ntpr, ntps=ntps, ntpg=ntpg)
      elseif (present(nr) .and. present(ns) .and. present(qh)) then
        call melting_model_jung08(qrreg, qsreg, qgreg, qms, qmg, qpr, qps, qpg, rats, ratg, &
             ntr=nrreg, nts=nsreg, ntms=ntms, ntpr=ntpr, ntps=ntps, &
             qh=qhreg, rath=rath, qmh=qmh, qph=qph)
      elseif (present(nr) .and. present(ns)) then
        call melting_model_jung08(qrreg, qsreg, qgreg, qms, qmg, qpr, qps, qpg, rats, ratg, &
             ntr=nrreg, nts=nsreg, ntms=ntms, ntpr=ntpr, ntps=ntps)
      elseif (present(nr) .and. present(qh)) then
        call melting_model_jung08(qrreg, qsreg, qgreg, qms, qmg, qpr, qps, qpg, rats, ratg, &
             ntr=nrreg, ntpr=ntpr, &
             qh=qhreg, rath=rath, qmh=qmh, qph=qph)
      elseif (present(nr)) then
        call melting_model_jung08(qrreg, qsreg, qgreg, qms, qmg, qpr, qps, qpg, rats, ratg, &
             ntr=nrreg, ntpr=ntpr)
      elseif (present(qh)) then
        call melting_model_jung08(qrreg, qsreg, qgreg, qms, qmg, qpr, qps, qpg, rats, ratg, &
             qh=qhreg, rath=rath, qmh=qmh, qph=qph)
      else
        call melting_model_jung08(qrreg, qsreg, qgreg, qms, qmg, qpr, qps, qpg, rats, ratg)
      endif

    case default
      ! 'liu24' (default): Liu et al. 2024 with modifications by Kong
      if (present(nr) .and. present(ns) .and. present(ng) .and. present(nh) .and. present(qh)) then
        call melting_model_liu24(qrreg, qsreg, qgreg, qms, qmg, qpr, qps, qpg, rats, ratg, &
             ntr=nrreg, nts=nsreg, ntg=ngreg, ntms=ntms, ntmg=ntmg, ntpr=ntpr, ntps=ntps, ntpg=ntpg, &
             qh=qhreg, nth=nhreg, rath=rath, qmh=qmh, qph=qph, ntmh=ntmh, ntph=ntph, &
             qx_min=qx_min, ntx_min=ntx_min, coeff_melt=coeff_melt, &
             density_air=density_air, temperature=temperature_celsius)
      elseif (present(nr) .and. present(ns) .and. present(ng) .and. present(qh)) then
        call melting_model_liu24(qrreg, qsreg, qgreg, qms, qmg, qpr, qps, qpg, rats, ratg, &
             ntr=nrreg, nts=nsreg, ntg=ngreg, ntms=ntms, ntmg=ntmg, ntpr=ntpr, ntps=ntps, ntpg=ntpg, &
             qh=qhreg, rath=rath, qmh=qmh, qph=qph, &
             qx_min=qx_min, ntx_min=ntx_min, coeff_melt=coeff_melt, &
             density_air=density_air, temperature=temperature_celsius)
      elseif (present(nr) .and. present(ns) .and. present(ng)) then
        call melting_model_liu24(qrreg, qsreg, qgreg, qms, qmg, qpr, qps, qpg, rats, ratg, &
             ntr=nrreg, nts=nsreg, ntg=ngreg, ntms=ntms, ntmg=ntmg, ntpr=ntpr, ntps=ntps, ntpg=ntpg, &
             qx_min=qx_min, ntx_min=ntx_min, coeff_melt=coeff_melt, &
             density_air=density_air, temperature=temperature_celsius)
      elseif (present(nr) .and. present(ns) .and. present(qh)) then
        call melting_model_liu24(qrreg, qsreg, qgreg, qms, qmg, qpr, qps, qpg, rats, ratg, &
             ntr=nrreg, nts=nsreg, ntms=ntms, ntpr=ntpr, ntps=ntps, &
             qh=qhreg, rath=rath, qmh=qmh, qph=qph, &
             qx_min=qx_min, ntx_min=ntx_min, coeff_melt=coeff_melt, &
             density_air=density_air, temperature=temperature_celsius)
      elseif (present(nr) .and. present(ns)) then
        call melting_model_liu24(qrreg, qsreg, qgreg, qms, qmg, qpr, qps, qpg, rats, ratg, &
             ntr=nrreg, nts=nsreg, ntms=ntms, ntpr=ntpr, ntps=ntps, &
             qx_min=qx_min, ntx_min=ntx_min, coeff_melt=coeff_melt, &
             density_air=density_air, temperature=temperature_celsius)
      elseif (present(nr) .and. present(qh)) then
        call melting_model_liu24(qrreg, qsreg, qgreg, qms, qmg, qpr, qps, qpg, rats, ratg, &
             ntr=nrreg, ntpr=ntpr, &
             qh=qhreg, rath=rath, qmh=qmh, qph=qph, &
             qx_min=qx_min, ntx_min=ntx_min, coeff_melt=coeff_melt, &
             density_air=density_air, temperature=temperature_celsius)
      elseif (present(nr)) then
        call melting_model_liu24(qrreg, qsreg, qgreg, qms, qmg, qpr, qps, qpg, rats, ratg, &
             ntr=nrreg, ntpr=ntpr, qx_min=qx_min, ntx_min=ntx_min, coeff_melt=coeff_melt, &
             density_air=density_air, temperature=temperature_celsius)
      elseif (present(qh)) then
        call melting_model_liu24(qrreg, qsreg, qgreg, qms, qmg, qpr, qps, qpg, rats, ratg, &
             qh=qhreg, rath=rath, qmh=qmh, qph=qph, &
             qx_min=qx_min, ntx_min=ntx_min, coeff_melt=coeff_melt, &
             density_air=density_air, temperature=temperature_celsius)
      else
        call melting_model_liu24(qrreg, qsreg, qgreg, qms, qmg, qpr, qps, qpg, rats, ratg, &
             qx_min=qx_min, ntx_min=ntx_min, coeff_melt=coeff_melt, &
             density_air=density_air, temperature=temperature_celsius)
      endif

    end select

    ! -------------------------------------------------------------------------
    ! Apply melting fraction tuning coefficients (default 1.0 = no change)
    ! -------------------------------------------------------------------------
    ! NOTE: At this point, rats/ratg/rath have already been scaled by
    ! melting_water_fraction inside the melting scheme (affecting mass
    ! partitioning: qmsr, qmgr, qpr, etc.).
    !
    ! This additional tuning further scales rats/ratg/rath for downstream
    ! scattering calculations only (does NOT re-partition mass):
    !   1. Scattering coefficients a(0:12) via coef_a(coefs, rats, a):
    !        a = c0 + c1*rats + c2*rats^2 + c3*rats^3   (Eq. 21, Zhang 2021)
    !   2. Melting species density via weticephase_density(rats, denms):
    !        Higher rats -> higher density (more liquid-like).
    !
    ! Physical meaning: tuning > 1 makes melting particles appear more
    ! liquid-like (higher density, rain-like scattering); tuning < 1 makes
    ! them appear more ice-like.
    ! -------------------------------------------------------------------------
    if (tuning_melt_frac_snow /= 1.0d0) then
      rats = min(tuning_melt_frac_snow * rats, 1.0d0)
    endif
    if (tuning_melt_frac_graupel /= 1.0d0) then
      ratg = min(tuning_melt_frac_graupel * ratg, 1.0d0)
    endif
    if (present(qh) .and. tuning_melt_frac_hail /= 1.0d0) then
      rath = min(tuning_melt_frac_hail * rath, 1.0d0)
    endif
    
    ! Regularize melting scheme outputs
    qpr = max(qx_min(1), qpr)
    qps = max(qx_min(5), qps)
    qpg = max(qx_min(6), qpg)
    if (present(qh) .and. present(nh)) then
      ! NSSL scheme - full regularization
      qph = max(qx_min(7), qph)
      ntpr = max(ntx_min(1), ntpr)
      ntms = max(ntx_min(2), ntms)
      ntmg = max(ntx_min(3), ntmg)
      ntmh = max(ntx_min(4), ntmh)
      ntps = max(ntx_min(5), ntps)
      ntpg = max(ntx_min(6), ntpg)
      ntph = max(ntx_min(7), ntph)
    elseif (present(nr) .and. present(ns) .and. present(ng)) then
      ! TCWA2 scheme - regularize all except hail
      ntpr = max(ntx_min(1), ntpr)
      ntms = max(ntx_min(2), ntms)
      ntmg = max(ntx_min(3), ntmg)
      ntps = max(ntx_min(5), ntps)
      ntpg = max(ntx_min(6), ntpg)
    elseif (present(nr)) then
      ! Thompson scheme - only rain
      ntpr = max(ntx_min(1), ntpr)
    endif
    
    ! Process all hydrometeor types
    ! 1. Pure rain
    if (skip_small_qx .and. qpr < 1.0d-3) then
      wpr = 0.0d0;  zpr = 0.0d0;  dmpr = 0.0d0
      zhpr = 0.0d0;  zdrpr = 1.0d0;  kdppr = 0.0d0;  phvpr = 1.0d0
    else
      call watercontent(density_air, qpr, wpr)
      if (present(nr)) then
        call dm_z_2moment(density_air, qpr, density_rain, ntpr, dmpr, zpr, dmmax=dmmax(1))
      else
        call dm_z_wsm6('rain', density_air, temp_air, qpr, density_rain, dmpr, zpr, dmmax=dmmax(1))
      endif
      ! Apply rain Dm tuning coefficient (default 1.0 = no change)
      if (tuning_dm_rain /= 1.0d0) then
        dmpr = tuning_dm_rain * dmpr
        zpr  = zpr * tuning_dm_rain**3  ! z ∝ dm^3
      endif
      call dualpol_op_rain(wpr, dmpr, rain_coefs, zhpr, zdrpr, kdppr, phvpr)
    endif
    
    ! 2. Melting snow
    if (skip_small_qx .and. qms < 0.001d0) then
      wms = 0.0d0;  zms = 0.0d0;  dmms = 0.0d0
      denms = density_snow
      zhms = 0.0d0;  zdrms = 1.0d0;  kdpms = 0.0d0;  phvms = 1.0d0
    else
      call watercontent(density_air, qms, wms)
      call weticephase_density('snow', rats, denms)
      ! Clamp melting snow density to [0.1, 0.5] g/cm³
      denms = min(0.5d0, max(0.1d0, denms))
      if (present(ns)) then
        call dm_z_2moment(density_air, qms, denms, ntms, dmms, zms, dmmax=dmmax(2), dm_min=0.3d0)
      else
        call dm_z_wsm6('snow', density_air, temp_air, qms, denms, dmms, zms, dmmax=dmmax(2), dm_min=0.3d0)
      endif
      ! Apply melting snow Dm tuning coefficient (default 1.0 = no change)
      if (tuning_dm_melting_snow /= 1.0d0) then
        dmms = tuning_dm_melting_snow * dmms
        zms  = zms * tuning_dm_melting_snow**3
      endif
      do i = 0,12
        call coef_a(snow_coefs(i,0:3), rats, as(i))
      end do
      call dualpol_op_icephase(zms, wms, denms, dmms, as(0:12), zhms, zdrms, kdpms, phvms)
    endif
    
    ! 3. Pure snow
    if (skip_small_qx .and. qps < 1.0d-4) then
      wps = 0.0d0;  zps = 0.0d0;  dmps = 0.0d0
      zhps = 0.0d0;  zdrps = 1.0d0;  kdpps = 0.0d0;  phvps = 1.0d0
    else
      call watercontent(density_air, qps, wps)
      if (present(ns)) then
        call dm_z_2moment(density_air, qps, density_snow, ntps, dmps, zps, dmmax=dmmax(5), dm_min=0.3d0)
      else
        call dm_z_wsm6('snow', density_air, temp_air, qps, density_snow, dmps, zps, dmmax=dmmax(5), dm_min=0.3d0)
      endif
      ! Apply pure snow Dm tuning coefficient (default 1.0 = no change)
      if (tuning_dm_pure_snow /= 1.0d0) then
        dmps = tuning_dm_pure_snow * dmps
        zps  = zps * tuning_dm_pure_snow**3
      endif
      call dualpol_op_icephase(zps, wps, density_snow, dmps, snow_a, zhps, zdrps, kdpps, phvps)
    endif
    
    ! 4. Melting graupel
    if (skip_small_qx .and. qmg < 0.01d0) then
      wmg = 0.0d0;  zmg = 0.0d0;  dmmg = 0.0d0
      denmg = density_graupel
      zhmg = 0.0d0;  zdrmg = 1.0d0;  kdpmg = 0.0d0;  phvmg = 1.0d0
    else
      call watercontent(density_air, qmg, wmg)
      call weticephase_density('graupel', ratg, denmg)
      if (present(ng)) then
        call dm_z_2moment(density_air, qmg, denmg, ntmg, dmmg, zmg, dmmax=dmmax(3), dm_min=0.3d0)
      else
        call dm_z_wsm6('graupel', density_air, temp_air, qmg, denmg, dmmg, zmg, dmmax=dmmax(3), dm_min=0.3d0)
      endif
      ! Apply melting graupel Dm tuning coefficient (default 1.0 = no change)
      if (tuning_dm_melting_graupel /= 1.0d0) then
        dmmg = tuning_dm_melting_graupel * dmmg
        zmg  = zmg * tuning_dm_melting_graupel**3
      endif
      do i = 0,12
        call coef_a(graupel_coefs(i,0:3), ratg, ag(i))
      end do
      call dualpol_op_icephase(zmg, wmg, denmg, dmmg, ag(0:12), zhmg, zdrmg, kdpmg, phvmg)
    endif
    
    ! 5. Pure graupel
    if (skip_small_qx .and. qpg < 0.001d0) then
      wpg = 0.0d0;  zpg = 0.0d0;  dmpg = 0.0d0
      denpg = density_graupel
      zhpg = 0.0d0;  zdrpg = 1.0d0;  kdppg = 0.0d0;  phvpg = 1.0d0
    else
      call watercontent(density_air, qpg, wpg)
      if (present(vg) .and. present(ng)) then
        if (qpg > 0.0 .and. vg > 0.0 .and. ng > 0) then
          denpg = density_air * qpg / vg * 1.0E-6
          denpg = min(max(denpg, 0.2_8), 0.6_8)
        else
          denpg = density_graupel
        endif
      else
        denpg = density_graupel
      endif
      if (present(ng)) then
        call dm_z_2moment(density_air, qpg, denpg, ntpg, dmpg, zpg, dmmax=dmmax(6), dm_min=0.3d0)
      else
        call dm_z_wsm6('graupel', density_air, temp_air, qpg, denpg, dmpg, zpg, dmmax=dmmax(6), dm_min=0.3d0)
      endif
      ! Apply pure graupel Dm tuning coefficient (default 1.0 = no change)
      if (tuning_dm_pure_graupel /= 1.0d0) then
        dmpg = tuning_dm_pure_graupel * dmpg
        zpg  = zpg * tuning_dm_pure_graupel**3
      endif
      call dualpol_op_icephase(zpg, wpg, denpg, dmpg, graupel_a, zhpg, zdrpg, kdppg, phvpg)
    endif
    
    ! 6-7. Hail (if present)
    if (present(nh) .and. present(qh)) then
      ! 6. Melting hail
      if (skip_small_qx .and. qmh < 0.01d0) then
        wmh = 0.0d0;  zmh = 0.0d0;  dmmh = 0.0d0
        if (treat_hail_as_graupel) then; denmh = density_graupel; else; denmh = density_hail; endif
        zhmh = 0.0d0;  zdrmh = 1.0d0;  kdpmh = 0.0d0;  phvmh = 1.0d0
      else
        call watercontent(density_air, qmh, wmh)
        if (treat_hail_as_graupel) then
          denmh = density_graupel
          call weticephase_density('graupel', rath, denmh)
        else
          denmh = density_hail
          call weticephase_density('hail', rath, denmh)
        endif
        if (present(nh)) then
          call dm_z_2moment(density_air, qmh, denmh, ntmh, dmmh, zmh, dmmax=dmmax(4), dm_min=0.3d0)
        else
          if (treat_hail_as_graupel) then
            call dm_z_wsm6('graupel', density_air, temp_air, qmh, denmh, dmmh, zmh, dmmax=dmmax(4), dm_min=0.3d0)
          else
            call dm_z_wsm6('hail', density_air, temp_air, qmh, denmh, dmmh, zmh, dmmax=dmmax(4), dm_min=0.3d0)
          endif
        endif
        ! Apply melting hail Dm tuning coefficient (default 1.0 = no change)
        if (tuning_dm_melting_hail /= 1.0d0) then
          dmmh = tuning_dm_melting_hail * dmmh
          zmh  = zmh * tuning_dm_melting_hail**3
        endif
        do i = 0,12
          if (treat_hail_as_graupel) then
            call coef_a(graupel_coefs(i,0:3), rath, ah(i))
          else
            call coef_a(hail_coefs(i,0:3), rath, ah(i))
          endif
        end do
        if (treat_hail_as_graupel) then
          call dualpol_op_icephase(zmh, wmh, denmh, dmmh, ah(0:12), zhmh, zdrmh, kdpmh, phvmh, iband=iband)
        else
          call dualpol_op_icephase(zmh, wmh, denmh, dmmh, ah(0:12), zhmh, zdrmh, kdpmh, phvmh, iband=iband, ptype='mh')
        endif
      endif
      
      ! 7. Pure hail
      if (skip_small_qx .and. qph < 0.001d0) then
        wph = 0.0d0;  zph = 0.0d0;  dmph = 0.0d0
        if (treat_hail_as_graupel) then; denph = density_graupel; else; denph = density_hail; endif
        zhph = 0.0d0;  zdrph = 1.0d0;  kdpph = 0.0d0;  phvph = 1.0d0
      else
        call watercontent(density_air, qph, wph)
        if (present(vh) .and. present(nh)) then
          if (qph > 0.0 .and. vh > 0.0 .and. nh > 0) then
            denph = density_air * qph / vh * 1.0E-6
            denph = min(max(denph, 0.7_8), 1.0_8)
          else
            if (treat_hail_as_graupel) then; denph = density_graupel; else; denph = density_hail; endif
          endif
        else
          if (treat_hail_as_graupel) then; denph = density_graupel; else; denph = density_hail; endif
        endif
        if (present(nh)) then
          call dm_z_2moment(density_air, qph, denph, ntph, dmph, zph, dmmax=dmmax(7), dm_min=0.3d0)
        else
          if (treat_hail_as_graupel) then
            call dm_z_wsm6('graupel', density_air, temp_air, qph, denph, dmph, zph, dmmax=dmmax(7), dm_min=0.3d0)
          else
            call dm_z_wsm6('hail', density_air, temp_air, qph, denph, dmph, zph, dmmax=dmmax(7), dm_min=0.3d0)
          endif
        endif
        ! Apply pure hail Dm tuning coefficient (default 1.0 = no change)
        if (tuning_dm_pure_hail /= 1.0d0) then
          dmph = tuning_dm_pure_hail * dmph
          zph  = zph * tuning_dm_pure_hail**3  ! z ∝ dm^3
        endif
        if (treat_hail_as_graupel) then
          do i = 0,12
            call coef_a(graupel_coefs(i,0:3), 0.0d0, ah(i))
          end do
          call dualpol_op_icephase(zph, wph, denph, dmph, ah(0:12), zhph, zdrph, kdpph, phvph, iband=iband)
        else
          call dualpol_op_icephase(zph, wph, denph, dmph, hail_a, zhph, zdrph, kdpph, phvph, &
                                   iband=iband, ptype='ph')
        endif
      endif
      
      ! Prepare dm array: [pr, ms, mg, mh, ps, pg, ph]
      dm(1) = dmpr
      dm(2) = dmms
      dm(3) = dmmg
      dm(4) = dmmh
      dm(5) = dmps
      dm(6) = dmpg
      dm(7) = dmph
      
      ! Total with hail
      call dualpol_op_total(zhpr, zhps, zhpg, zdrpr, zdrps, zdrpg, &
                            kdppr, kdpps, kdppg, phvpr, phvps, phvpg, &
                            zh, zdr, kdp, phv, &
                            zh_msnow=zhms, zh_mgraupel=zhmg, &
                            zdr_msnow=zdrms, zdr_mgraupel=zdrmg, &
                            kdp_msnow=kdpms, kdp_mgraupel=kdpmg, &
                            phv_msnow=phvms, phv_mgraupel=phvmg, &
                            zh_mhail=zhmh, zdr_mhail=zdrmh, kdp_mhail=kdpmh, phv_mhail=phvmh, &
                            zh_phail=zhph, zdr_phail=zdrph, kdp_phail=kdpph, phv_phail=phvph, &
                            dmms=dmms, dmmg=dmmg, dmmh=dmmh, dm=dm, temperature=temperature)
    else
      ! Prepare dm array without hail
      dm(1) = dmpr
      dm(2) = dmms
      dm(3) = dmmg
      dm(4) = 0.0d0  ! No hail
      dm(5) = dmps
      dm(6) = dmpg
      dm(7) = 0.0d0  ! No hail
      
      ! Total without hail
      call dualpol_op_total(zhpr, zhps, zhpg, zdrpr, zdrps, zdrpg, &
                            kdppr, kdpps, kdppg, phvpr, phvps, phvpg, &
                            zh, zdr, kdp, phv, &
                            zh_msnow=zhms, zh_mgraupel=zhmg, &
                            zdr_msnow=zdrms, zdr_mgraupel=zdrmg, &
                            kdp_msnow=kdpms, kdp_mgraupel=kdpmg, &
                            phv_msnow=phvms, phv_mgraupel=phvmg, &
                            dm=dm, temperature=temperature)
    endif
    
  end subroutine zhang21_compute_point

  ! ------------------------------------------------------------------------------

end module zhang21_forward_mod
