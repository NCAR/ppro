module dualpol_op_mod
!--------------------------
! PPRO Library Main Interface Module
!
! This module provides the unified interface for PPRO library, supporting
! multiple polarimetric radar operators:
!   - Zhang21: Zhang et al. 2021 (lookup table based)
!   - TCWA2: Tsai et al. (empirically-fitted polynomials based on gamma PSD)
!
! Authors:
!   Rong Kong (NCAR/MMM) - Interface design and multi-operator support, 2025
!--------------------------------------
  use zhang21_core_mod
  use tcwa2_core_mod
  implicit none

  ! Operator type enumeration
  integer, parameter, public :: OPERATOR_ZHANG21 = 1
  integer, parameter, public :: OPERATOR_TCWA2   = 2
  
  ! Current operator selection
  integer, save :: current_operator = OPERATOR_ZHANG21
  
  ! Initialization flag
  logical, save :: initialized = .false.

  ! Public interfaces
  public :: ppro_init_coefs
  public :: ppro_compute_point
  public :: ppro_finalize
  public :: ppro_set_operator
  public :: ppro_get_operator
  
  ! Re-export needed variables from zhang21 for backward compatibility
  public :: density_rain, density_snow, density_graupel, density_hail
  public :: ratio_rain, ratio_snow, ratio_graupel, ratio_hail
  public :: sband_snow_a, sband_graupel_a, sband_hail_a
  public :: cband_snow_a, cband_graupel_a, cband_hail_a
  public :: sband_rain_coefs, sband_snow_coefs, sband_graupel_coefs, sband_hail_coefs
  public :: cband_rain_coefs, cband_snow_coefs, cband_graupel_coefs, cband_hail_coefs
  
  ! Re-export melting scheme parameters (configurable via YAML in UFO layer)
  public :: melting_scheme_option
  public :: enable_melting_transition
  public :: snow_ratio_low, snow_ratio_high
  public :: graupel_ratio_low, graupel_ratio_high
  public :: hail_ratio_low, hail_ratio_high
  public :: melting_water_fraction
  public :: dmmax_rain
  public :: dmmax_pure_snow, dmmax_melting_snow
  public :: dmmax_pure_graupel, dmmax_melting_graupel
  public :: dmmax_pure_hail, dmmax_melting_hail
  public :: treat_hail_as_graupel
  public :: enable_dm_melting_limit, dm_melting_transition_width
  public :: skip_small_qx
  public :: enable_dm_regularization  ! Dm regularization for small ntx/qx
  public :: tuning_dm_rain
  public :: tuning_dm_melting_snow, tuning_dm_melting_graupel, tuning_dm_melting_hail
  public :: tuning_dm_pure_snow, tuning_dm_pure_graupel, tuning_dm_pure_hail
  public :: melting_rain_exponent
  public :: tuning_melt_frac_snow, tuning_melt_frac_graupel, tuning_melt_frac_hail

contains

  ! ------------------------------------------------------------------------------
  !> @brief Initialize PPRO library with specified operator
  !>
  !> @param[in] coef_path       Optional path to coefficient files (Zhang21 only)
  !> @param[in] operator_name   Optional operator name: 'Zhang21' or 'TCWA2'
  ! ------------------------------------------------------------------------------
  subroutine ppro_init_coefs(coef_path, operator_name)
    implicit none
    character(len=*), intent(in), optional :: coef_path
    character(len=*), intent(in), optional :: operator_name
    
    ! Set operator type if provided
    if (present(operator_name)) then
      call ppro_set_operator(operator_name)
    endif
    
    if (initialized) then
      write(*,'(A)') 'PPRO: Already initialized, skipping'
      return
    endif
    
    ! Initialize operator-specific components
    if (current_operator == OPERATOR_ZHANG21) then
      call zhang21_init_coefs(coef_path)
    else if (current_operator == OPERATOR_TCWA2) then
      write(*,'(A)') '=========================================='
      write(*,'(A)') 'PPRO: Using TCWA2 operator'
      write(*,'(A)') 'No coefficient files needed (built-in fitted polynomials)'
      write(*,'(A)') '=========================================='
    endif
    
    initialized = .true.
    
  end subroutine ppro_init_coefs

  ! ------------------------------------------------------------------------------
  !> @brief Main computation interface - dispatches to selected operator
  ! ------------------------------------------------------------------------------
  subroutine ppro_compute_point(iband, scheme_type, density_air, temp_air, &
                                qr, qs, qg, zh, zdr, kdp, phv, &
                                qh, nr, ns, ng, nh, vg, vh, qi, ni, qc, smlf, gmlf, &
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
    real(kind=8), intent(in), optional :: qi, ni, qc, smlf, gmlf
    real(kind=8), intent(in), optional :: height, zhobs, zdrobs, kdpobs
    real(kind=8), intent(in), optional :: coeff_melt, temperature
    integer, intent(in), optional :: iobs
    
    ! Dispatch to selected operator
    if (current_operator == OPERATOR_ZHANG21) then
      ! Zhang21: Use lookup tables with melting layer treatment
      ! Note: vg, vh not passed to disable dynamic density calculation (use fixed density)
      call zhang21_compute_point(iband, scheme_type, density_air, temp_air, &
                                 qr, qs, qg, zh, zdr, kdp, phv, &
                                 qh, nr, ns, ng, nh, &
                                 height=height, zhobs=zhobs, zdrobs=zdrobs, kdpobs=kdpobs, &
                                 coeff_melt=coeff_melt, temperature=temperature, iobs=iobs)
      
    else if (current_operator == OPERATOR_TCWA2) then
      ! TCWA2: Use empirically-fitted polynomial formulation
      if (.not. (present(qi) .and. present(ni) .and. present(qc) .and. &
                 present(smlf) .and. present(gmlf) .and. present(nr) .and. &
                 present(ns) .and. present(ng))) then
        write(*,*) 'ERROR: TCWA2 operator requires qi, ni, qc, smlf, gmlf, nr, ns, ng'
        zh = 0.0_8; zdr = 0.0_8; kdp = 0.0_8; phv = 1.0_8
        return
      endif
      
      ! Call TCWA2 operator (need to implement aggregation wrapper)
      call tcwa2_compute_point(iband, density_air, temp_air, qr, qs, qg, qc, qi, &
                               nr, ns, ng, ni, smlf, gmlf, zh, zdr, kdp, phv)
    endif
    
  end subroutine ppro_compute_point

  ! ------------------------------------------------------------------------------
  !> @brief Set the polarimetric operator type
  ! ------------------------------------------------------------------------------
  subroutine ppro_set_operator(operator_name)
    implicit none
    character(len=*), intent(in) :: operator_name
    
    select case (trim(operator_name))
      case ('Zhang21', 'ZHANG21', 'zhang21')
        current_operator = OPERATOR_ZHANG21
        write(*,*) 'PPRO: Selected Zhang21 operator (lookup table based)'
      case ('TCWA2', 'tcwa2', 'Tcwa2')
        current_operator = OPERATOR_TCWA2
        write(*,*) 'PPRO: Selected TCWA2 operator (fitted polynomial formulation)'
      case default
        write(*,*) 'WARNING: Unknown operator "', trim(operator_name), '". Using Zhang21.'
        current_operator = OPERATOR_ZHANG21
    end select
  end subroutine ppro_set_operator

  ! ------------------------------------------------------------------------------
  !> @brief Get current operator type
  ! ------------------------------------------------------------------------------
  function ppro_get_operator() result(operator_type)
    implicit none
    integer :: operator_type
    operator_type = current_operator
  end function ppro_get_operator

  ! ------------------------------------------------------------------------------
  !> @brief Finalize and cleanup
  ! ------------------------------------------------------------------------------
  subroutine ppro_finalize()
    implicit none
    
    if (current_operator == OPERATOR_ZHANG21) then
      call zhang21_finalize()
    endif
    
    initialized = .false.
    write(*,'(A)') 'PPRO: Finalized'
  end subroutine ppro_finalize

  ! ------------------------------------------------------------------------------
  !> @brief TCWA2 aggregation wrapper (calls individual TCWA2 subroutines)
  ! ------------------------------------------------------------------------------
  subroutine tcwa2_compute_point(iband, density_air, temp_air, qr, qs, qg, qc, qi, &
                                 nr, ns, ng, ni, smlf, gmlf, zh, zdr, kdp, phv)
    implicit none
    integer, intent(in) :: iband
    real(kind=8), intent(in) :: density_air, temp_air
    real(kind=8), intent(in) :: qr, qs, qg, qc, qi
    real(kind=8), intent(in) :: nr, ns, ng, ni
    real(kind=8), intent(in) :: smlf, gmlf
    real(kind=8), intent(out) :: zh, zdr, kdp, phv
    
    ! Local variables for individual hydrometeor contributions
    real(kind=8) :: zh_rain, zv_rain, kdp_rain
    real(kind=8) :: zh_ice, zv_ice, kdp_ice
    real(kind=8) :: zh_snow, zv_snow, kdp_snow
    real(kind=8) :: zh_graup, zv_graup, kdp_graup
    real(kind=8) :: zdr_temp
    
    ! Call TCWA2 subroutines for each hydrometeor type
    call dualpol_op_rain_tcwa2(iband, density_air, qr, nr, zh_rain, zv_rain, kdp_rain)
    call dualpol_op_ice_tcwa2(iband, temp_air, density_air, qi, ni, zh_ice, zv_ice, kdp_ice)
    call dualpol_op_snow_tcwa2(iband, temp_air, density_air, qc, qr, qs, ns, smlf, &
                               zh_snow, zv_snow, kdp_snow)
    call dualpol_op_graup_tcwa2(iband, temp_air, density_air, qc, qr, qg, ng, gmlf, &
                                zh_graup, zv_graup, kdp_graup)
    
    ! Aggregate contributions (linear sum for zh, zv, kdp)
    zh = zh_rain + zh_ice + zh_snow + zh_graup
    kdp = kdp_rain + kdp_ice + kdp_snow + kdp_graup
    
    ! Compute ZDR from zh and total zv
    zdr_temp = (zv_rain + zv_ice + zv_snow + zv_graup) / max(zh, 1.0e-20_8)
    zdr = zdr_temp  ! Already in linear units, not dB
    
    ! Compute correlation coefficient (simplified, TCWA2 doesn't provide phv directly)
    phv = 0.95_8  ! Placeholder - TCWA2 operator doesn't compute phv explicitly
    
  end subroutine tcwa2_compute_point

end module dualpol_op_mod
