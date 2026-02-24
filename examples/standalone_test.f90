!> @file standalone_test.f90
!> @brief Standalone test program for PPRO library
!>
!> This program demonstrates how to use the PPRO library independently
!> to compute dual-polarization radar variables from model state variables.
!>
!> @Rong Kong NCAR/MMM
!> @10/10 2025

program standalone_test
  use dualpol_op_mod, only: ppro_init_coefs, ppro_compute_point, ppro_finalize
  implicit none
  
  ! Configuration
  integer, parameter :: iband = 1  ! 1=S-band, 2=C-band
  character(len=20) :: scheme_type
  
  ! Input atmospheric state variables
  real(kind=8) :: density_air   ! kg/m^3
  real(kind=8) :: temp_air      ! K
  real(kind=8) :: qr, qs, qg    ! g/kg (rain, snow, graupel mixing ratios)
  real(kind=8) :: qh            ! g/kg (hail mixing ratio, for NSSL)
  real(kind=8) :: nr, ns, ng, nh ! #/kg (number concentrations, for 2-moment schemes)
  
  ! Output dual-pol variables
  real(kind=8) :: zh, zdr, kdp, phv
  real(kind=8) :: zh_dbz, zdr_db
  
  ! Loop variables
  integer :: itest
  
  write(*,*) "=============================================="
  write(*,*) "PPRO Standalone Test Program"
  write(*,*) "Multi-Operator Support (Zhang21 & TCWA2)"
  write(*,*) "=============================================="
  write(*,*)
  
  write(*,*) "=== Testing Zhang21 Operator ==="
  write(*,*)
  
  ! Initialize for Zhang21 operator (requires coefficient files)
  call ppro_init_coefs(operator_name='Zhang21')
  write(*,*)
  
  ! ============================================================================
  ! Test Case 1: WSM6 scheme (single-moment)
  ! ============================================================================
  write(*,*) "----------------------------------------------"
  write(*,*) "Test Case 1: WSM6 (Single-moment)"
  write(*,*) "----------------------------------------------"
  
  scheme_type = 'WSM6'
  density_air = 1.0_8      ! kg/m^3
  temp_air = 280.0_8       ! K (7°C)
  qr = 1.0_8               ! g/kg
  qs = 0.5_8               ! g/kg
  qg = 0.3_8               ! g/kg
  
  call ppro_compute_point(iband, scheme_type, density_air, temp_air, &
                          qr, qs, qg, zh, zdr, kdp, phv)
  
  ! Convert to dBZ and dB
  zh_dbz = 10.0_8 * log10(zh + 1.0e-10_8)
  zdr_db = 10.0_8 * log10(zdr + 1.0e-10_8)
  
  write(*,'(A,F8.3,A)') "  Temperature:     ", temp_air, " K"
  write(*,'(A,F8.3,A)') "  Density:         ", density_air, " kg/m^3"
  write(*,'(A,F8.4,A)') "  Qr (rain):       ", qr, " g/kg"
  write(*,'(A,F8.4,A)') "  Qs (snow):       ", qs, " g/kg"
  write(*,'(A,F8.4,A)') "  Qg (graupel):    ", qg, " g/kg"
  write(*,*)
  write(*,'(A,F10.2,A)') "  Zhh:             ", zh_dbz, " dBZ"
  write(*,'(A,F10.3,A)') "  ZDR:             ", zdr_db, " dB"
  write(*,'(A,F10.4,A)') "  KDP:             ", kdp, " deg/km"
  write(*,'(A,F10.4)')   "  ρhv:             ", phv
  write(*,*)
  
  ! ============================================================================
  ! Test Case 2: Thompson scheme (double-moment for rain)
  ! ============================================================================
  write(*,*) "----------------------------------------------"
  write(*,*) "Test Case 2: Thompson (Double-moment)"
  write(*,*) "----------------------------------------------"
  
  scheme_type = 'THOMPSON'
  density_air = 1.1_8
  temp_air = 275.0_8       ! K (2°C, near melting)
  qr = 2.0_8
  qs = 0.8_8
  qg = 0.5_8
  nr = 1.0e6_8             ! #/kg (rain number concentration)
  
  call ppro_compute_point(iband, scheme_type, density_air, temp_air, &
                          qr, qs, qg, zh, zdr, kdp, phv, &
                          nr=nr)
  
  zh_dbz = 10.0_8 * log10(zh + 1.0e-10_8)
  zdr_db = 10.0_8 * log10(zdr + 1.0e-10_8)
  
  write(*,'(A,F8.3,A)') "  Temperature:     ", temp_air, " K"
  write(*,'(A,F8.3,A)') "  Density:         ", density_air, " kg/m^3"
  write(*,'(A,F8.4,A)') "  Qr (rain):       ", qr, " g/kg"
  write(*,'(A,F8.4,A)') "  Qs (snow):       ", qs, " g/kg"
  write(*,'(A,F8.4,A)') "  Qg (graupel):    ", qg, " g/kg"
  write(*,'(A,E12.4,A)') "  Nr (rain #):     ", nr, " #/kg"
  write(*,*)
  write(*,'(A,F10.2,A)') "  Zhh:             ", zh_dbz, " dBZ"
  write(*,'(A,F10.3,A)') "  ZDR:             ", zdr_db, " dB"
  write(*,'(A,F10.4,A)') "  KDP:             ", kdp, " deg/km"
  write(*,'(A,F10.4)')   "  ρhv:             ", phv
  write(*,*)
  
  ! ============================================================================
  ! Test Case 3: NSSL scheme (double-moment with hail)
  ! ============================================================================
  write(*,*) "----------------------------------------------"
  write(*,*) "Test Case 3: NSSL (Double-moment with hail)"
  write(*,*) "----------------------------------------------"
  
  scheme_type = 'NSSL'
  density_air = 1.0_8
  temp_air = 278.0_8       ! K (5°C)
  qr = 3.0_8
  qs = 0.3_8
  qg = 0.4_8
  qh = 0.5_8               ! Hail
  nr = 5.0e5_8
  ns = 1.0e6_8
  ng = 1.0e5_8
  nh = 1.0e4_8
  
  call ppro_compute_point(iband, scheme_type, density_air, temp_air, &
                          qr, qs, qg, zh, zdr, kdp, phv, &
                          qh=qh, nr=nr, ns=ns, ng=ng, nh=nh)
  
  zh_dbz = 10.0_8 * log10(zh + 1.0e-10_8)
  zdr_db = 10.0_8 * log10(zdr + 1.0e-10_8)
  
  write(*,'(A,F8.3,A)') "  Temperature:     ", temp_air, " K"
  write(*,'(A,F8.3,A)') "  Density:         ", density_air, " kg/m^3"
  write(*,'(A,F8.4,A)') "  Qr (rain):       ", qr, " g/kg"
  write(*,'(A,F8.4,A)') "  Qs (snow):       ", qs, " g/kg"
  write(*,'(A,F8.4,A)') "  Qg (graupel):    ", qg, " g/kg"
  write(*,'(A,F8.4,A)') "  Qh (hail):       ", qh, " g/kg"
  write(*,'(A,E12.4,A)') "  Nr (rain #):     ", nr, " #/kg"
  write(*,'(A,E12.4,A)') "  Nh (hail #):     ", nh, " #/kg"
  write(*,*)
  write(*,'(A,F10.2,A)') "  Zhh:             ", zh_dbz, " dBZ"
  write(*,'(A,F10.3,A)') "  ZDR:             ", zdr_db, " dB"
  write(*,'(A,F10.4,A)') "  KDP:             ", kdp, " deg/km"
  write(*,'(A,F10.4)')   "  ρhv:             ", phv
  write(*,*)
  
  ! ============================================================================
  ! Test Case 4: Temperature sensitivity test (melting layer)
  ! ============================================================================
  write(*,*) "----------------------------------------------"
  write(*,*) "Test Case 4: Temperature Sensitivity"
  write(*,*) "(Fixed hydrometeors, varying temperature)"
  write(*,*) "----------------------------------------------"
  
  scheme_type = 'WSM6'
  density_air = 1.0_8
  qr = 1.5_8
  qs = 1.0_8
  qg = 0.5_8
  
  write(*,'(A10,5A12)') "Temp(K)", "Zhh(dBZ)", "ZDR(dB)", "KDP(°/km)", "ρhv"
  write(*,'(A10,5A12)') "----------", "------------", "------------", "------------", "------------"
  
  do itest = 1, 11
    temp_air = 268.0_8 + real(itest-1, 8) * 2.0_8  ! 268K to 288K (step 2K)
    
    call ppro_compute_point(iband, scheme_type, density_air, temp_air, &
                            qr, qs, qg, zh, zdr, kdp, phv)
    
    zh_dbz = 10.0_8 * log10(zh + 1.0e-10_8)
    zdr_db = 10.0_8 * log10(zdr + 1.0e-10_8)
    
    write(*,'(F10.2,5F12.4)') temp_air, zh_dbz, zdr_db, kdp, phv
  end do
  
  write(*,*)
  write(*,*) "=============================================="
  write(*,*) "Zhang21 tests completed!"
  write(*,*) "=============================================="
  write(*,*)
  
  ! Cleanup Zhang21
  call ppro_finalize()
  
  ! ============================================================================
  ! Test TCWA2 Operator
  ! ============================================================================
  write(*,*) "=== Testing TCWA2 Operator ==="
  write(*,*)
  
  ! Initialize for TCWA2 operator (no coefficient files needed)
  call ppro_init_coefs(operator_name='TCWA2')
  
  write(*,*) "----------------------------------------------"
  write(*,*) "Test Case: TCWA2 Fitted Analytic Formulation"
  write(*,*) "----------------------------------------------"
  write(*,*) "Note: TCWA2 assumes gamma PSD; parameterizations derived from"
  write(*,*) "      bin-based scattering under the Rayleigh approximation."
  write(*,*)
  write(*,*) "(TCWA2 requires additional parameters: qi, ni, qc, smlf, gmlf)"
  write(*,*) "(Full TCWA2 testing requires UFO interface)"
  write(*,*)
  
  write(*,*) "=============================================="
  write(*,*) "All tests completed successfully!"
  write(*,*) "Supported operators: Zhang21, TCWA2"
  write(*,*) "=============================================="
  
  ! Final cleanup
  call ppro_finalize()
  
end program standalone_test


