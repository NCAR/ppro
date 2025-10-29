module zhang21_core_mod
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

  ! Public interfaces
  public :: zhang21_compute_point
  public :: zhang21_init_coefs
  public :: zhang21_finalize
  
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
    
    ! Try 2: base_path/fix/filename  
    fullpath = trim(base_path) // '/fix/' // trim(filename)
    inquire(file=trim(fullpath), exist=file_exists)
    if (file_exists) return
    
    ! Try 3: Just filename (current directory)
    fullpath = trim(filename)
    inquire(file=trim(fullpath), exist=file_exists)
    if (file_exists) return
    
    ! Try 4: ./fix/filename
    fullpath = './fix/' // trim(filename)
    
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
    if(present(iband) .and. iband ==2 .and. present(ptype) .and. ptype == "ph") then 
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
                               dmms, dmmg, dmmh) !optional
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

    ! Pure snow and graupel (optional)
    real(kind=8), intent(in), optional   :: zh_msnow, zh_mgraupel    ! horizontal reflectivity, mm6 m-3
    real(kind=8), intent(inout), optional   :: zdr_msnow, zdr_mgraupel ! differential reflectivity, -
    real(kind=8), intent(in), optional   :: kdp_msnow, kdp_mgraupel ! specific differential phase, degree/km^-1
    real(kind=8), intent(in), optional   :: phv_msnow, phv_mgraupel ! co-polar correlation coefficient, 0-1

    ! Hail and pure hail (optional)
    real(kind=8), intent(in), optional :: zh_mhail, zh_phail    ! horizontal reflectivity, mm6 m-3
    real(kind=8), intent(in), optional :: zdr_phail  ! differential reflectivity, -
    real(kind=8), intent(inout), optional :: zdr_mhail ! differential reflectivity, -
    real(kind=8), intent(in), optional :: kdp_mhail, kdp_phail  ! specific differential phase, degree/km^-1
    real(kind=8), intent(in), optional :: phv_mhail, phv_phail  ! co-polar correlation coefficient, 0-1
    real(kind=8), intent(in), optional :: dmms, dmmg, dmmh 

    real(kind=8) :: zbar_prain, zbar_msnow, zbar_mgraupel, zbar_mhail
    real(kind=8) :: zv_prain, zv_msnow, zv_mgraupel, zv_mhail
    real(kind=8) :: zbar_psnow, zbar_pgraupel, zbar_phail
    real(kind=8) :: zv_psnow, zv_pgraupel, zv_phail

    real(kind=8) :: zh_qc_prain, zh_qc_psnow, zh_qc_pgraupel, zh_qc_phail 
    real(kind=8) :: zv_qc_prain, zv_qc_psnow, zv_qc_pgraupel, zv_qc_phail 

    real(kind=8) :: zh_qc_msnow, zh_qc_mgraupel, zh_qc_mhail
    real(kind=8) :: zv_qc_msnow, zv_qc_mgraupel, zv_qc_mhail

    real(kind=8) :: zhqc, zv, zvqc, alpha
    real(kind=8) :: zh_thresh = 5.0

    real(kind=8) :: zdr_dB_limit_msnow, zdr_dB_limit_mgraupel, zdr_dB_limit_mhail
    real(kind=8) :: zdr_limit_msnow, zdr_limit_mgraupel, zdr_limit_mhail
  
 
    logical :: include_hail = .false.
    logical :: pure_hydro = .false.

    if (present(zh_mhail) .and. present(zh_phail)) then
      include_hail = .true.
    endif

    ! Here the zdr is dimensionless, not in unit of dBZ, therefore always positive
    if(zdr_prain > 1.0D-3)then
       zv_prain = zh_prain/zdr_prain
    else
       zv_prain = 0.0d0
    endif

    if(zdr_prain > 1.0-3)then
       zbar_prain = zh_prain / sqrt(zdr_prain)
    else
       zbar_prain = 0.0d0
    endif

    if(zdr_msnow > 1.0D-3)then
       zv_msnow = zh_msnow/zdr_msnow
    else
       zv_msnow = 0.0d0
    endif

    if(zdr_msnow > 1.0D-3)then
       zbar_msnow = zh_msnow / sqrt(zdr_msnow)
    else
       zbar_msnow = 0.0d0
    endif

    if(zdr_mgraupel>1.0D-3)then
      zv_mgraupel = zh_mgraupel/zdr_mgraupel
      zbar_mgraupel = zh_mgraupel / sqrt(zdr_mgraupel)
    else
      zv_mgraupel = 0.0d0
      zbar_mgraupel = 0.0d0
    endif

    if(zdr_psnow > 1.0D-3)then
       zv_psnow = zh_psnow/zdr_psnow
       zbar_psnow = zh_psnow / sqrt(zdr_psnow)
    else
       zv_psnow = 0.0d0
       zbar_psnow = 0.0d0
    endif

    if(zdr_pgraupel > 1.0D-3)then
       zv_pgraupel = zh_pgraupel/zdr_pgraupel
       zbar_pgraupel = zh_pgraupel / sqrt(zdr_pgraupel)
    else
       zv_pgraupel = 0.0d0
       zbar_pgraupel = 0.0d0
    endif


    zh_qc_prain    = zh_prain
    zh_qc_psnow    = zh_psnow
    zh_qc_pgraupel = zh_pgraupel
    zh_qc_phail    = zh_phail

    zv_qc_prain    = zv_prain
    zv_qc_psnow    = zv_psnow
    zv_qc_pgraupel = zv_pgraupel
    zv_qc_phail    = zv_phail

    zh_qc_msnow    = zh_msnow
    zv_qc_msnow    = zv_msnow
    zh_qc_mgraupel = zh_mgraupel
    zv_qc_mgraupel = zv_mgraupel

    if (include_hail) then
       if(zdr_mhail > 1.0D-3)then
          zv_mhail = zh_mhail/zdr_mhail
          zbar_mhail = zh_msnow / sqrt(zdr_mhail)
       else
          zv_mhail = 0.0d0
          zbar_mhail = 0.0d0
       endif

       if(zdr_phail > 1.0D-3)then
          zv_phail = zh_phail/zdr_phail
          zbar_phail = zh_phail / sqrt(zdr_phail)
       else
          zv_phail = 0.0d0
          zbar_phail = 0.0d0
       endif
       zh_qc_mhail    = zh_mhail
       zv_qc_mhail    = zv_mhail
    endif


    if (include_hail) then
       if(pure_hydro)then
          zh  = zh_prain + zh_psnow + zh_pgraupel + zh_phail ! Eq. 22
          zv  = zv_prain + zv_psnow + zv_pgraupel + zv_phail
          if(zv > 1.0D-3)then
             zdr = zh/zv ! Eq. 23
          else
             zdr = 0.0d0
          endif
          kdp = kdp_prain + kdp_psnow + kdp_pgraupel + kdp_phail ! Eq. 24
          if(zbar_prain + zbar_psnow + zbar_pgraupel + zbar_phail > 1.0D-3)then
             phv = (zbar_prain*phv_prain + &
                    zbar_psnow*phv_psnow + zbar_pgraupel*phv_pgraupel + zbar_phail*phv_phail) / &
                    (zbar_prain + zbar_psnow + zbar_pgraupel + zbar_phail) ! Eq. 25
          else
             phv = 0.0d0
          endif
       else
          zh  = zh_prain + zh_msnow + zh_mgraupel + zh_psnow + zh_pgraupel + zh_mhail + zh_phail ! Eq. 22
          zv = zv_prain + zv_msnow + zv_mgraupel + zv_psnow + zv_pgraupel + zv_mhail + zv_phail

          if(zv > 1.0D-3)then
             zdr = zh/zv ! Eq. 23
          else
             zdr = 0.0d0
          endif
          kdp = kdp_prain + kdp_msnow + kdp_mgraupel + kdp_psnow + kdp_pgraupel + kdp_mhail + kdp_phail ! Eq. 24
          if(zbar_prain + zbar_msnow + zbar_mgraupel + zbar_psnow + zbar_pgraupel + zbar_mhail + zbar_phail > 1.0D-3)then
             phv = (zbar_prain*phv_prain + zbar_msnow*phv_msnow + zbar_mgraupel*phv_mgraupel + zbar_mhail*phv_mhail + &
                    zbar_psnow*phv_psnow + zbar_pgraupel*phv_pgraupel + zbar_phail*phv_phail) / &
                    (zbar_prain + zbar_msnow + zbar_mgraupel + zbar_psnow + zbar_pgraupel + zbar_mhail + zbar_phail) ! Eq. 25
          else
             phv = 0.0d0
          endif    
       endif
    else
       zh  = zh_prain + zh_msnow + zh_mgraupel + zh_psnow + zh_pgraupel ! Eq. 22

       if(zv_prain + zv_msnow + zv_mgraupel + zv_psnow + zv_pgraupel > 1.0D-3)then
          zdr = zh/(zv_prain + zv_msnow + zv_mgraupel + zv_psnow + zv_pgraupel) ! Eq. 23
       else
          zdr = 0.0d0
       endif
      
       kdp = kdp_prain + kdp_msnow + kdp_mgraupel + kdp_psnow + kdp_pgraupel ! Eq. 24

       if(zbar_prain + zbar_msnow + zbar_mgraupel + zbar_psnow + zbar_pgraupel > 1.0D-3)then
          phv = (zbar_prain*phv_prain + zbar_msnow*phv_msnow + zbar_mgraupel*phv_mgraupel + &
                 zbar_psnow*phv_psnow + zbar_pgraupel*phv_pgraupel) / &
                 (zbar_prain + zbar_msnow + zbar_mgraupel + zbar_psnow + zbar_pgraupel) ! Eq. 25
       else
          phv = 0.0d0
       endif
    endif

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
    if(qr + qx > 1.0E-5)then
       ratio = qr/(qr+qx)
    else
      ratio = 0.0d0
    endif

  end subroutine rainfrac

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
  subroutine dm_z_2moment(air_density, qx, rhox, ntx, dm, z, dmmax)
  !-----------------------------------------------------------
  ! Eqs. (5) and (6) of Zhang et al., 2021
  ! Assume 2-moment MP scheme with ntx and qx as prognostic variables
  ! Exponential PSD: N(D) = N0 * exp(-lambda D)
  !-------------------------------------------------------------------
    implicit none
    real(kind=8), intent(in)   :: air_density   ! Air Density, 1.225 kg/m^3 
    real(kind=8), intent(in)   :: rhox          ! density of wet hydrometeor,g/cm^3 
                                                ! rain density = 1 g/cm^3
    real(kind=8), intent(in)   :: qx            ! mixing ratio of hydrometeor, g/kg
    real(kind=8), intent(in)   :: ntx           ! number concentration, ?/m^3
    real(kind=8), intent(in), optional   :: dmmax         ! maximum mean diameter
    real(kind=8), intent(out)  :: dm            ! mean diameter, mm
    real(kind=8), intent(out)  :: z             ! 6th moment of PSD

    dm = 4.0 * ((air_density * qx) / (dpi * rhox * ntx)) ** (1.0 / 3.0)
    dm = dm * 10.0 ! cm -> mm   
    if(present(dmmax))then
       dm = min(dm, dmmax) !to avoid extremely large diameter due to very small ntx
    endif

     z = 11250.0 * ( (air_density*qx)/(dpi*rhox) )* dm**3

  end subroutine dm_z_2moment

   !---------------------------------------------------------------------------
  subroutine melting_scheme_zhang24( &
    qr, qs, qg, &
    qms, qmg, qpr, qps, qpg, rats, ratg, &
    ntr, nts, ntg, &
    ntms, ntmg, ntpr, ntps, ntpg, &
    qh, nth, rath, qmh, qph, ntmh, ntph, qx_min, ntx_min )
  !---------------------------------------------------------------------------
  ! Eq. (15-18) of Zhang et al., 2024
  ! Assume 2-moment MP scheme with ntx and qx as prognostic variables
  ! Exponential PSD: N(D) = N0 * exp(-lambda D)
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
    logical      :: dm ! double moment
    real(kind=8), intent(in), optional  ::  qx_min(7), ntx_min(7)

   ! Check if double moment

    dm = .false.
    if (present(ntr) .or. present(nts) .or. present(ntg) .or. present(nth)) then
        dm = .true.
    else
        print *, "Please specify output for number concentrations of the melting scheme"
        stop
    endif
    !print*,'dm=',dm
    qms  = sqrt(qr * qs)
    qmg  = sqrt(qr * qg)
    if (present(qh)) qmh = sqrt(qr * qh)
    qms = max(qms, qx_min(2))
    qmg = max(qmg, qx_min(3))
    if (present(qh)) qmh = max(qmh, qx_min(4))

    IF (qr + qs .le. 1.0D-3   .or. ntr + nts < 100.0 ) THEN
      rats = 0.0
    ELSE
      rats = qr / (qr + qs)
    ENDIF

    IF (qr + qg .le. 1.0D-3 .or. ntr + ntg < 100.0 ) THEN
      ratg = 0.0
    ELSE
      ratg = qr / (qr + qg)
    ENDIF

    if (present(qh)) then
      if (qr + qh .le. 1.0D-3 .or. ntr + nth < 100.0 ) THEN
        rath = 0.0
      else
        rath = qr / (qr + qh)
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

    qpr = max(qpr, qx_min(1))
    qps = max(qps, qx_min(5))
    qpg = max(qpg, qx_min(6))
    if (present(qh)) qph = max(qph, qx_min(7))



    if (dm) then
      if(present(nts) .and. present(ntr) .and. present(ntms)) then
         ntms = sqrt(ntr * nts)
         ntms = max(ntms, ntx_min(2))
      endif

      if(present(ntg) .and. present(ntr) .and. present(ntmg)) then
         ntmg = sqrt(ntr * ntg)
         ntmg = max(ntmg, ntx_min(3))
      endif

      if (present(nth) .and. present(ntr) .and. present(ntmh)) then
         ntmh = sqrt(ntr * nth)
         ntmh = max(ntmh, ntx_min(4))
      endif

      if(present(ntr)) ntpr = max(ntr,ntx_min(1))
      if(present(nts)) ntps = max(nts,ntx_min(5))
      if(present(ntg)) ntpg = max(ntg,ntx_min(6))
      if(present(nth)) ntph = max(nth,ntx_min(7))

    endif

  end subroutine melting_scheme_zhang24

  !---------------------------------------------------------------------------
  subroutine melting_scheme_zhang08( &
    qr, qs, qg, &
    qms, qmg, qpr, qps, qpg, rats, ratg, &
    ntr, nts, ntg, &
    ntms, ntmg, ntpr, ntps, ntpg, &
    qh, nth, rath, qmh, qph, ntmh, ntph)
  !---------------------------------------------------------------------------
  ! Eq. (11-14) of Zhang et al., 2024
  ! Assume 2-moment MP scheme with ntx and qx as prognostic variables
  ! Exponential PSD: N(D) = N0 * exp(-lambda D)
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
    logical      :: dm ! double moment
    real(kind=8) :: minratio

    dm = .false.
    if (present(ntr) .and. present(nts) .and. present(ntg)) then
      if (present(ntms) .and. present(ntmg) .and. present(ntpr) .and. present(ntps) .and. present(ntpg)) then
        dm = .true.
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

    if (dm) then
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
  
  end subroutine melting_scheme_zhang08

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
  subroutine dm_z_wsm6(precip_type, air_density, temp, qx, rhox, dm, z, dmmax)
  !------------------------------------------------------------------------
  ! Single moment scheme assuming Exponential PSD: N(D) = N0 * exp(-lambda D)
  ! N0 is fixed except for snow, then calcuate lambda from N0
  !-------------------------------------------------------------------
    implicit none
    character (len=*), intent(in) :: precip_type  ! rain/snow/graupel
    real(kind=8), intent(in)   :: air_density   ! Air Density, 1.225 kg/m^3 
    real(kind=8), intent(in)   :: rhox          ! density of wet hydrometeor,g/cm^3 
                                                ! rain density = 1 g/cm^3
    real(kind=8), intent(in)   :: temp          ! temperature, Kelvin
    real(kind=8), intent(in)   :: qx            ! mixing ratio of hydrometeor, g/kg
    real(kind=8), intent(in),optional   :: dmmax         ! maximum mean diameter
    real(kind=8), intent(out)  :: dm            ! mean diameter, mm
    real(kind=8), intent(out)  :: z             ! 6th moment of PSD

    real(kind=8)  :: n0  ! intercept parameter of exponential PSD
    real(kind=8)  :: lambda ! slope parameter of exponential PSD

    call n0_lambda_wsm6(precip_type, air_density, temp, qx, rhox, n0, lambda)
   
    dm = 4000.0 / lambda  ! Eq. (5), should be in mm, *1000?

    if(present(dmmax))then
       dm = min(dm, dmmax) !to avoid extremely large diameter due to very small ntx
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
                                qh, nr, ns, ng, nh, vg, vh)
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
    
    ! Initialize parameters
    qx_min = 0.0d0
    ntx_min = 1.0D-3
    dmmax = 25.0d0
    
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
    if (present(nr) .and. present(ns) .and. present(ng) .and. present(nh) .and. present(qh)) then
      ! NSSL scheme
      call melting_scheme_zhang24(qrreg, qsreg, qgreg, qms, qmg, qpr, qps, qpg, rats, ratg, &
           ntr=nrreg, nts=nsreg, ntg=ngreg, ntms=ntms, ntmg=ntmg, ntpr=ntpr, ntps=ntps, ntpg=ntpg, &
           qh=qhreg, nth=nhreg, rath=rath, qmh=qmh, qph=qph, ntmh=ntmh, ntph=ntph, &
           qx_min=qx_min, ntx_min=ntx_min)
    elseif (present(nr) .and. present(ns) .and. present(ng)) then
      ! TCWA2 scheme
      call melting_scheme_zhang24(qrreg, qsreg, qgreg, qms, qmg, qpr, qps, qpg, rats, ratg, &
           ntr=nrreg, nts=nsreg, ntg=ngreg, ntms=ntms, ntmg=ntmg, ntpr=ntpr, ntps=ntps, ntpg=ntpg, &
           qx_min=qx_min, ntx_min=ntx_min)
    elseif (present(nr) .and. .not. present(ns) .and. .not. present(ng)) then
      ! Thompson scheme
      call melting_scheme_zhang24(qrreg, qsreg, qgreg, qms, qmg, qpr, qps, qpg, rats, ratg, &
           ntr=nrreg, ntpr=ntpr, qx_min=qx_min, ntx_min=ntx_min)
    else
      ! WSM6 scheme
      call melting_scheme_zhang24(qrreg, qsreg, qgreg, qms, qmg, qpr, qps, qpg, rats, ratg, &
           qx_min=qx_min, ntx_min=ntx_min)
    endif
    
    ! Regularize melting scheme outputs
    qpr = max(qx_min(1), qpr)
    qps = max(qx_min(5), qps)
    qpg = max(qx_min(6), qpg)
    if (present(qh) .and. present(nh)) then
      qph = max(qx_min(7), qph)
      ntpr = max(ntx_min(1), ntpr)
      ntms = max(ntx_min(2), ntms)
      ntmg = max(ntx_min(3), ntmg)
      ntmh = max(ntx_min(4), ntmh)
      ntps = max(ntx_min(5), ntps)
      ntpg = max(ntx_min(6), ntpg)
      ntph = max(ntx_min(7), ntph)
    elseif (present(nr)) then
      ntpr = max(ntx_min(1), ntpr)
    endif
    
    ! Process all hydrometeor types
    ! 1. Pure rain
    call watercontent(density_air, qpr, wpr)
    if (present(nr)) then
      call dm_z_2moment(density_air, qpr, density_rain, ntpr, dmpr, zpr, dmmax=dmmax(1))
    else
      call dm_z_wsm6('rain', density_air, temp_air, qpr, density_rain, dmpr, zpr, dmmax=dmmax(1))
    endif
    call dualpol_op_rain(wpr, dmpr, rain_coefs, zhpr, zdrpr, kdppr, phvpr)
    
    ! 2. Melting snow
    call watercontent(density_air, qms, wms)
    call weticephase_density('snow', rats, denms)
    if (present(ns)) then
      call dm_z_2moment(density_air, qms, denms, ntms, dmms, zms, dmmax=dmmax(2))
    else
      call dm_z_wsm6('snow', density_air, temp_air, qms, denms, dmms, zms, dmmax=dmmax(2))
    endif
    do i = 0,12
      call coef_a(snow_coefs(i,0:3), rats, as(i))
    end do
    call dualpol_op_icephase(zms, wms, denms, dmms, as(0:12), zhms, zdrms, kdpms, phvms)
    
    ! 3. Pure snow
    call watercontent(density_air, qps, wps)
    if (present(ns)) then
      call dm_z_2moment(density_air, qps, density_snow, ntps, dmps, zps, dmmax=dmmax(5))
    else
      call dm_z_wsm6('snow', density_air, temp_air, qps, density_snow, dmps, zps, dmmax=dmmax(5))
    endif
    call dualpol_op_icephase(zps, wps, density_snow, dmps, snow_a, zhps, zdrps, kdpps, phvps)
    
    ! 4. Melting graupel
    call watercontent(density_air, qmg, wmg)
    call weticephase_density('graupel', ratg, denmg)
    if (present(ng)) then
      call dm_z_2moment(density_air, qmg, denmg, ntmg, dmmg, zmg, dmmax=dmmax(3))
    else
      call dm_z_wsm6('graupel', density_air, temp_air, qmg, denmg, dmmg, zmg, dmmax=dmmax(3))
    endif
    do i = 0,12
      call coef_a(graupel_coefs(i,0:3), ratg, ag(i))
    end do
    call dualpol_op_icephase(zmg, wmg, denmg, dmmg, ag(0:12), zhmg, zdrmg, kdpmg, phvmg)
    
    ! 5. Pure graupel
    call watercontent(density_air, qpg, wpg)
    if (present(vg)) then
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
      call dm_z_2moment(density_air, qpg, denpg, ntpg, dmpg, zpg, dmmax=dmmax(6))
    else
      call dm_z_wsm6('graupel', density_air, temp_air, qpg, denpg, dmpg, zpg, dmmax=dmmax(6))
    endif
    call dualpol_op_icephase(zpg, wpg, density_graupel, dmpg, graupel_a, zhpg, zdrpg, kdppg, phvpg)
    
    ! 6-7. Hail (if present)
    if (present(nh) .and. present(qh)) then
      ! 6. Melting hail
      call watercontent(density_air, qmh, wmh)
      call weticephase_density('hail', rath, denmh)
      if (present(nh)) then
        call dm_z_2moment(density_air, qmh, denmh, ntmh, dmmh, zmh, dmmax=dmmax(4))
      else
        call dm_z_wsm6('hail', density_air, temp_air, qmh, denmh, dmmh, zmh, dmmax=dmmax(4))
      endif
      do i = 0,12
        call coef_a(hail_coefs(i,0:3), rath, ah(i))
      end do
      call dualpol_op_icephase(zmh, wmh, denmh, dmmh, ah(0:12), zhmh, zdrmh, kdpmh, phvmh, iband=iband)
      
      ! 7. Pure hail
      call watercontent(density_air, qph, wph)
      if (present(vh)) then
        if (qph > 0.0 .and. vh > 0.0 .and. nh > 0) then
          denph = density_air * qph / vh * 1.0E-6
          denph = min(max(denph, 0.7_8), 1.0_8)
        else
          denph = density_hail
        endif
      else
        denph = density_hail
      endif
      if (present(nh)) then
        call dm_z_2moment(density_air, qph, denph, ntph, dmph, zph, dmmax=dmmax(7))
      else
        call dm_z_wsm6('hail', density_air, temp_air, qph, denph, dmph, zph, dmmax=dmmax(7))
      endif
      call dualpol_op_icephase(zph, wph, denph, dmph, hail_a, zhph, zdrph, kdpph, phvph, &
                               iband=iband, ptype='ph')
      
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
                            dmms=dmms, dmmg=dmmg, dmmh=dmmh)
    else
      ! Total without hail
      call dualpol_op_total(zhpr, zhps, zhpg, zdrpr, zdrps, zdrpg, &
                            kdppr, kdpps, kdppg, phvpr, phvps, phvpg, &
                            zh, zdr, kdp, phv, &
                            zh_msnow=zhms, zh_mgraupel=zhmg, &
                            zdr_msnow=zdrms, zdr_mgraupel=zdrmg, &
                            kdp_msnow=kdpms, kdp_mgraupel=kdpmg, &
                            phv_msnow=phvms, phv_mgraupel=phvmg)
    endif
    
  end subroutine zhang21_compute_point

  ! ------------------------------------------------------------------------------

end module zhang21_core_mod
