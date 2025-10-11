module zhang21_tlad_mod
!--------------------------
! The tangent linear (TL) / adjoint (AD) module of the 
! Parameterized polarimetric radar operator following
! Zhang, G., J. Gao, and M. Du, 2021: Parameterized forward operators 
!   for simulation and assimilation of polarimetric radar data with 
!   numerical weather predictions. Adv. Atmos. Sci., 38(5), 737−754.
!
! Author: Hejun Xie, NCAR/MMM
! Initial version: May 2024
!--------------------------------------

  use zhang21_core_mod, only: density_rain, density_snow, density_graupel, density_hail, &
                              ratio_rain, ratio_snow, ratio_graupel, ratio_hail, &
                              sband_snow_a, sband_graupel_a, sband_hail_a, &
                              cband_snow_a, cband_graupel_a, cband_hail_a, &
                              sband_rain_coefs, sband_snow_coefs, sband_graupel_coefs, sband_hail_coefs, &
                              cband_rain_coefs, cband_snow_coefs, cband_graupel_coefs, cband_hail_coefs, &
                              n0_lambda_wsm6, n0_lambda_gceop

  implicit none

  ! TL / AD module use double real type
  INTEGER, PARAMETER :: rd = 8

  ! Trajectory variables recorded from forward operator

  ! input of melting scheme
  real(kind=rd) qr_traj, qs_traj, qg_traj, qh_traj
  real(kind=rd) ntr_traj, nts_traj, ntg_traj, nth_traj

  ! output of melting scheme
  real(kind=rd), target :: qpr_traj, qms_traj, qmg_traj, qps_traj, qpg_traj, qmh_traj, qph_traj        ! mixing ratio of hydrometeor, g/kg
  real(kind=rd), target :: ntpr_traj, ntms_traj, ntmg_traj, ntps_traj, ntpg_traj, ntmh_traj, ntph_traj ! number concentration, ?/m^3

  ! "r" is "pr", pure rain
  ! "s" is "ms", melting snow
  ! "g" is "mg", melting graupel
  ! "h" is "mh", melting hail
  real(kind=rd), target :: wr_traj, ws_traj, wg_traj, wps_traj, wpg_traj, wh_traj, wph_traj           ! mass concentration, g/m^3
  real(kind=rd), target :: dmr_traj, dms_traj, dmg_traj, dmps_traj, dmpg_traj, dmh_traj, dmph_traj    ! mean diameter, mm
  real(kind=rd), target :: zr_traj, zs_traj, zg_traj, zps_traj, zpg_traj, zmh_traj, zph_traj           ! 6th moment of PSD
  real(kind=rd), target :: dens_traj, deng_traj, denh_traj                                            ! density of wet hydrometeor g/cm^3
  real(kind=rd), target :: rats_traj, ratg_traj, rath_traj                                            ! mass water fraction [0-1]

  real(kind=rd), target :: zhr_traj, zhs_traj, zhg_traj, zhps_traj, zhpg_traj, zhh_traj, zhph_traj                ! horizontal reflectivity [mm^6 m^-3]
  real(kind=rd), target :: zdrr_traj, zdrs_traj, zdrg_traj, zdrps_traj, zdrpg_traj, zdrh_traj, zdrph_traj         ! differential reflectivity [-]
  real(kind=rd), target :: zbarr_traj, zbars_traj, zbarg_traj, zbarps_traj, zbarpg_traj, zbarh_traj, zbarph_traj  ! intermediate reflectivity [mm^6 m^-3]
  real(kind=rd), target :: zvr_traj, zvs_traj, zvg_traj, zvps_traj, zvpg_traj, zvh_traj, zvph_traj                ! vertical reflectivity [mm^6 m^-3]
  real(kind=rd), target :: kdpr_traj, kdps_traj, kdpg_traj, kdpps_traj, kdppg_traj, kdph_traj, kdpph_traj         ! differential phase shift [deg km^-1]
  real(kind=rd), target :: phvr_traj, phvs_traj, phvg_traj, phvps_traj, phvpg_traj, phvh_traj, phvph_traj         ! copolar correlation coefficient [0-1]

  real(kind=rd), target :: as_traj(0:12), ag_traj(0:12), ah_traj(0:12)    ! coefficients a for wet hydrometeors
  real(kind=rd), pointer :: aps_traj(:) => NULL()       ! coefficients a for dry hydrometeors, pointer to *band_snow_a
  real(kind=rd), pointer :: apg_traj(:) => NULL()       ! coefficients a for dry hydrometeors, pointer to *band_graupel_a
  real(kind=rd), pointer :: aph_traj(:) => NULL() 

  real(kind=rd), pointer :: qx_traj => NULL()
  real(kind=rd), pointer :: wx_traj => NULL()
  real(kind=rd), pointer :: ntx_traj => NULL()
  real(kind=rd), pointer :: dmx_traj => NULL()
  real(kind=rd), pointer :: zx_traj => NULL()
  real(kind=rd), pointer :: denx_traj => NULL()
  real(kind=rd), pointer :: ratx_traj => NULL()
  real(kind=rd), pointer :: zhx_traj => NULL()
  real(kind=rd), pointer :: zdrx_traj => NULL()
  real(kind=rd), pointer :: zbarx_traj => NULL()
  real(kind=rd), pointer :: zvx_traj => NULL()
  real(kind=rd), pointer :: kdpx_traj => NULL()
  real(kind=rd), pointer :: phvx_traj => NULL()
  real(kind=rd), pointer :: ax_traj(:) => NULL()

  real(kind=rd) density_air_traj ! kg/m3
  real(kind=rd) temp_air_traj    ! Kelvin
  real(kind=rd) zh_traj
  real(kind=rd) zdr_traj                        
  real(kind=rd) kdp_traj                        
  real(kind=rd) phv_traj                        

  logical :: rain_populated = .FALSE.
  logical :: snow_populated = .FALSE.
  logical :: graupel_populated = .FALSE.
  logical :: pure_snow_populated = .FALSE.
  logical :: pure_graupel_populated = .FALSE.
  logical :: hail_populated = .FALSE.
  logical :: pure_hail_populated = .FALSE.
  logical :: var_den_rat = .FALSE. ! Whether the ratio and density of that hydrometeor is variant

  real(kind=rd), parameter :: zero = 0.0_rd

contains

  subroutine populate_traj(precip_type, qx, wx, ntx, dmx, zx, denx, ratx, &
                           zhx, zdrx, kdpx, phvx, ax)
    !-----------------------------------------------------------
    ! Polulate forward trajectory of a given precipitation type 
    !-----------------------------------------------------------
    implicit none

    character (len=*), intent(in) :: precip_type
    real(kind=rd), intent(in)   :: qx          ! mixing ratio of hydrometeor, g/kg
    real(kind=rd), intent(in)   :: wx          ! mass concentration, g/m^3
    real(kind=rd), intent(in)   :: ntx         ! number concentration, ?/m^3
    real(kind=rd), intent(in)   :: dmx         ! mean diameter, mm
    real(kind=rd), intent(in)   :: zx          ! 6th moment of PSD
    real(kind=rd), intent(in)   :: denx        ! density of wet hydrometeor g/cm^3
    real(kind=rd), intent(in)   :: ratx        ! mass water fraction [0-1]

    real(kind=rd), intent(in)   :: zhx         ! horizontal reflectivity [mm^6 m^-3]
    real(kind=rd), intent(in)   :: zdrx        ! differential reflectivity [-]
    real(kind=rd), intent(in)   :: kdpx        ! differential phase shift [deg km^-1]
    real(kind=rd), intent(in)   :: phvx        ! copolar correlation coefficient [0-1]
    real(kind=rd), intent(in), optional, target :: ax(0:12)        ! copolar correlation coefficient [0-1]

    select case (precip_type)
      case ('rain')
        qpr_traj = qx
        ntpr_traj = ntx
        wr_traj = wx
        dmr_traj = dmx
        zr_traj = zx
        zhr_traj = zhx
        zdrr_traj = zdrx
        kdpr_traj = kdpx
        phvr_traj = phvx
        zbarr_traj = zhx / sqrt(zdrx)
        zvr_traj = zhx / zdrx
        rain_populated = .TRUE.
      case ('snow')
        qms_traj = qx
        ntms_traj = ntx
        ws_traj = wx
        dms_traj = dmx
        zs_traj = zx
        dens_traj = denx
        rats_traj = ratx
        zhs_traj = zhx
        zdrs_traj = zdrx
        kdps_traj = kdpx
        phvs_traj = phvx
        as_traj = ax
        zbars_traj = zhx / sqrt(zdrx)
        zvs_traj = zhx / zdrx
        snow_populated = .TRUE.
      case ('graupel')
        qmg_traj = qx
        ntmg_traj = ntx
        wg_traj = wx
        dmg_traj = dmx
        zg_traj = zx
        deng_traj = denx
        ratg_traj = ratx
        zhg_traj = zhx
        zdrg_traj = zdrx
        kdpg_traj = kdpx
        phvg_traj = phvx
        ag_traj = ax
        zbarg_traj = zhx / sqrt(zdrx)
        zvg_traj = zhx / zdrx
        graupel_populated = .TRUE.
      case ('pure snow')
        qps_traj = qx
        ntps_traj = ntx
        wps_traj = wx
        dmps_traj = dmx
        zps_traj = zx
        zhps_traj = zhx
        zdrps_traj = zdrx
        kdpps_traj = kdpx
        phvps_traj = phvx
        aps_traj => ax
        zbarps_traj = zhx / sqrt(zdrx)
        zvps_traj = zhx / zdrx
        pure_snow_populated = .TRUE.
      case ('pure graupel')
        qpg_traj = qx
        ntpg_traj = ntx
        wpg_traj = wx
        dmpg_traj = dmx
        zpg_traj = zx
        zhpg_traj = zhx
        zdrpg_traj = zdrx
        kdppg_traj = kdpx
        phvpg_traj = phvx
        apg_traj => ax
        zbarpg_traj = zhx / sqrt(zdrx)
        zvpg_traj = zhx / zdrx
        pure_graupel_populated = .TRUE.
      case ('hail')
        qmh_traj = qx
        ntmh_traj = ntx
        wh_traj = wx
        dmh_traj = dmx
        zmh_traj = zx
        denh_traj = denx
        rath_traj = ratx
        zhh_traj = zhx
        zdrh_traj = zdrx
        kdph_traj = kdpx
        phvh_traj = phvx
        ah_traj = ax
        zbarh_traj = zhx / sqrt(zdrx)
        zvh_traj = zhx / zdrx
        hail_populated = .TRUE.
      case ('pure hail')
        qph_traj = qx
        ntph_traj = ntx
        wph_traj = wx
        dmph_traj = dmx
        zph_traj = zx
        zhph_traj = zhx
        zdrph_traj = zdrx
        kdpph_traj = kdpx
        phvph_traj = phvx
        aph_traj => ax
        zbarph_traj = zhx / sqrt(zdrx)
        zvph_traj = zhx / zdrx
        pure_hail_populated = .TRUE.
      case default
        write(*,*) 'Unknow precipitation type: ', precip_type
        stop
    end select
  end subroutine populate_traj

  subroutine set_traj_precip_type(precip_type)
    !-----------------------------------------------------------
    ! Set the trajectory pointer to a given precipitation type
    !-----------------------------------------------------------
    implicit none
    
    character (len=*), intent(in) :: precip_type
    select case (precip_type)
      case ('rain')
        qx_traj => qpr_traj
        ntx_traj => ntpr_traj
        wx_traj => wr_traj
        dmx_traj => dmr_traj
        zx_traj => zr_traj
        denx_traj => density_rain
        ratx_traj => ratio_rain
        zhx_traj => zhr_traj
        zdrx_traj => zdrr_traj
        kdpx_traj => kdpr_traj
        phvx_traj => phvr_traj
        zbarx_traj => zbarr_traj
        zvx_traj => zvr_traj
        ax_traj => NULL()
        var_den_rat = .FALSE.
      case ('snow')
        qx_traj => qms_traj
        ntx_traj => ntms_traj
        wx_traj => ws_traj
        dmx_traj => dms_traj
        zx_traj => zs_traj
        denx_traj => dens_traj
        ratx_traj => rats_traj
        zhx_traj => zhs_traj
        zdrx_traj => zdrs_traj
        kdpx_traj => kdps_traj
        phvx_traj => phvs_traj
        zbarx_traj => zbars_traj
        zvx_traj => zvs_traj
        ax_traj => as_traj
        var_den_rat = .TRUE.
      case ('graupel')
        qx_traj => qmg_traj
        ntx_traj => ntmg_traj
        wx_traj => wg_traj
        dmx_traj => dmg_traj
        zx_traj => zg_traj
        denx_traj => deng_traj
        ratx_traj => ratg_traj
        zhx_traj => zhg_traj
        zdrx_traj => zdrg_traj
        kdpx_traj => kdpg_traj
        phvx_traj => phvg_traj
        zbarx_traj => zbarg_traj
        zvx_traj => zvg_traj
        ax_traj => ag_traj
        var_den_rat = .TRUE.
      case ('pure snow')
        qx_traj => qps_traj
        ntx_traj => ntps_traj
        wx_traj => wps_traj
        dmx_traj => dmps_traj
        zx_traj => zps_traj
        denx_traj => density_snow
        ratx_traj => ratio_snow
        zhx_traj => zhps_traj
        zdrx_traj => zdrps_traj
        kdpx_traj => kdpps_traj
        phvx_traj => phvps_traj
        zbarx_traj => zbarps_traj
        zvx_traj => zvps_traj
        ax_traj => aps_traj
        var_den_rat = .FALSE.
      case ('pure graupel')
        qx_traj => qpg_traj
        ntx_traj => ntpg_traj
        wx_traj => wpg_traj
        dmx_traj => dmpg_traj
        zx_traj => zpg_traj
        denx_traj => density_graupel
        ratx_traj => ratio_graupel
        zhx_traj => zhpg_traj
        zdrx_traj => zdrpg_traj
        kdpx_traj => kdppg_traj
        phvx_traj => phvpg_traj
        zbarx_traj => zbarpg_traj
        zvx_traj => zvpg_traj
        ax_traj => apg_traj
        var_den_rat = .FALSE.
      case ('hail')
        qx_traj => qmh_traj
        ntx_traj => ntmh_traj
        wx_traj => wh_traj
        dmx_traj => dmh_traj
        zx_traj => zmh_traj
        denx_traj => denh_traj
        ratx_traj => rath_traj
        zhx_traj => zhh_traj
        zdrx_traj => zdrh_traj
        kdpx_traj => kdph_traj
        phvx_traj => phvh_traj
        zbarx_traj => zbarh_traj
        zvx_traj => zvh_traj
        ax_traj => ah_traj
        var_den_rat = .TRUE.
      case ('pure hail')
        qx_traj => qph_traj
        ntx_traj => ntph_traj
        wx_traj => wph_traj
        dmx_traj => dmph_traj
        zx_traj => zph_traj
        denx_traj => density_hail
        ratx_traj => ratio_hail
        zhx_traj => zhph_traj
        zdrx_traj => zdrph_traj
        kdpx_traj => kdpph_traj
        phvx_traj => phvph_traj
        zbarx_traj => zbarph_traj
        zvx_traj => zvph_traj
        ax_traj => aph_traj
        var_den_rat = .FALSE. 
      case default
        write(*,*) 'Unknow precipitation type: ', precip_type
        stop
    end select
  end subroutine set_traj_precip_type

  subroutine dm_z_2moment_tl(qx_tl, ntx_tl, denx_tl, dmx_tl, zx_tl)
    !-----------------------------------------------------------
    ! Tangent Linear (TL) of subroutine dm_z_2moment
    ! Assume 2-moment MP scheme with ntx and qx as prognostic variables
    ! Exponential PSD: N(D) = N0 * exp(-lambda D)
    !-------------------------------------------------------------------
    implicit none
    real(kind=rd), intent(in)  :: qx_tl          ! mixing ratio of hydrometeor, g/kg
    real(kind=rd), intent(in)  :: ntx_tl         ! number concentration, ?/m^3
    real(kind=rd), intent(in)  :: denx_tl        ! wet hydrometeor density, g/cm^3
    real(kind=rd), intent(out) :: dmx_tl         ! mean diameter, mm
    real(kind=rd), intent(out) :: zx_tl          ! 6th moment of PSD
    
    dmx_tl = (dmx_traj / 3.0_rd) * (qx_tl / qx_traj - ntx_tl / ntx_traj - denx_tl / denx_traj)
    zx_tl = zx_traj * (qx_tl / qx_traj - denx_tl / denx_traj + 3.0_rd * dmx_tl / dmx_traj)

  end subroutine dm_z_2moment_tl

  subroutine dm_z_2moment_ad(qx_ad, ntx_ad, denx_ad, dmx_ad, zx_ad)
    !-----------------------------------------------------------
    ! Adjoint (AD) of subroutine dm_z_2moment
    ! Assume 2-moment MP scheme with ntx and qx as prognostic variables
    ! Exponential PSD: N(D) = N0 * exp(-lambda D)
    !-------------------------------------------------------------------
    implicit none
    real(kind=rd), intent(inout)   :: qx_ad         ! mixing ratio of hydrometeor, g/kg
    real(kind=rd), intent(inout)   :: ntx_ad        ! number concentration, ?/m^3
    real(kind=rd), intent(inout)   :: denx_ad       ! wet hydrometeor density, g/cm^3
    real(kind=rd), intent(in)      :: dmx_ad        ! mean diameter, mm
    real(kind=rd), intent(in)      :: zx_ad         ! 6th moment of PSD

    real(kind=rd) :: dmx_ad_tmp = zero

    dmx_ad_tmp = dmx_ad + 3.0_rd * zx_traj / dmx_traj * zx_ad
    if (var_den_rat) &
      denx_ad = denx_ad - dmx_traj / (3.0_rd * denx_traj) * dmx_ad_tmp - zx_traj / denx_traj * zx_ad
    ntx_ad = ntx_ad - dmx_traj / (3.0_rd * ntx_traj) * dmx_ad_tmp
    qx_ad = qx_ad + dmx_traj / (3.0_rd * qx_traj) * dmx_ad_tmp + zx_traj / qx_traj * zx_ad 
  
  end subroutine dm_z_2moment_ad

  subroutine dm_z_wsm6_tl(precip_type, qx_tl, denx_tl, dmx_tl, zx_tl)
    !-----------------------------------------------------------
    ! Tangent Linear (TL) of subroutine dm_z_wsm6
    ! Single moment scheme assuming Exponential PSD: N(D) = N0 * exp(-lambda D)
    ! N0 is fixed except for snow, then calcuate lambda from N0
    !-------------------------------------------------------------------
    implicit none
    character (len=*), intent(in) :: precip_type
    real(kind=rd), intent(in)  :: qx_tl          ! mixing ratio of hydrometeor, g/kg
    real(kind=rd), intent(in)  :: denx_tl        ! wet hydrometeor density, g/cm^3
    real(kind=rd), intent(out) :: dmx_tl         ! mean diameter, mm
    real(kind=rd), intent(out) :: zx_tl          ! 6th moment of PSD

    real(kind=rd) :: n0_traj
    real(kind=rd) :: n0_tl
    real(kind=rd) :: lambda_traj
    real(kind=rd) :: lambda_tl

    call n0_lambda_wsm6(precip_type, density_air_traj, temp_air_traj, qx_traj, denx_traj, n0_traj, lambda_traj)
    
    lambda_tl = 0.25_rd * lambda_traj * ( denx_tl / denx_traj - qx_tl / qx_traj )
    dmx_tl = - dmx_traj * lambda_tl / lambda_traj
    zx_tl = zx_traj * (qx_tl / qx_traj - denx_tl / denx_traj + 3.0_rd * dmx_tl / dmx_traj)

  end subroutine dm_z_wsm6_tl

  subroutine dm_z_wsm6_ad(precip_type, qx_ad, denx_ad, dmx_ad, zx_ad)
    !-----------------------------------------------------------
    ! Adjoint (AD) of subroutine dm_z_wsm6
    ! Single moment scheme assuming Exponential PSD: N(D) = N0 * exp(-lambda D)
    ! N0 is fixed except for snow, then calcuate lambda from N0
    !-------------------------------------------------------------------
    implicit none
    character (len=*), intent(in) :: precip_type
    real(kind=rd), intent(inout)  :: qx_ad          ! mixing ratio of hydrometeor, g/kg
    real(kind=rd), intent(inout)  :: denx_ad        ! wet hydrometeor density, g/cm^3
    real(kind=rd), intent(in)     :: dmx_ad         ! mean diameter, mm
    real(kind=rd), intent(in)     :: zx_ad          ! 6th moment of PSD

    real(kind=rd) :: n0_traj
    real(kind=rd) :: lambda_traj
    real(kind=rd) :: lambda_ad
    real(kind=rd) :: dmx_ad_tmp = zero

    call n0_lambda_wsm6(precip_type, density_air_traj, temp_air_traj, qx_traj, denx_traj, n0_traj, lambda_traj)

    dmx_ad_tmp = dmx_ad + 3.0_rd * zx_traj * zx_ad / dmx_traj
    lambda_ad = - dmx_traj * dmx_ad_tmp / lambda_traj
    qx_ad = qx_ad + ( -0.25_rd*lambda_traj*lambda_ad + zx_traj*zx_ad) / qx_traj
    if (var_den_rat) &
      denx_ad = denx_ad + ( 0.25_rd*lambda_traj*lambda_ad - zx_traj*zx_ad ) / denx_traj

  end subroutine dm_z_wsm6_ad

  subroutine dm_z_gceop_tl(precip_type, qx_tl, denx_tl, dmx_tl, zx_tl)
    !-----------------------------------------------------------
    ! Tangent Linear (TL) of subroutine dm_z_gceop
    ! N0 is fixed, calcuate lambda from N0
    ! Goddard Cumulus Ensemble scheme used in CWB WRF operation
    ! only difference from wsm6: n0 for snow is a fixed value
    !-------------------------------------------------------------------
    implicit none
    character (len=*), intent(in) :: precip_type
    real(kind=rd), intent(in)  :: qx_tl          ! mixing ratio of hydrometeor, g/kg
    real(kind=rd), intent(in)  :: denx_tl        ! wet hydrometeor density, g/cm^3
    real(kind=rd), intent(out) :: dmx_tl         ! mean diameter, mm
    real(kind=rd), intent(out) :: zx_tl          ! 6th moment of PSD

    real(kind=rd) :: n0_traj
    real(kind=rd) :: lambda_traj
    real(kind=rd) :: lambda_tl

    call n0_lambda_gceop(precip_type, density_air_traj, qx_traj, denx_traj, n0_traj, lambda_traj)

    lambda_tl = 0.25_rd * lambda_traj * ( denx_tl / denx_traj - qx_tl / qx_traj )
    dmx_tl = - dmx_traj * lambda_tl / lambda_traj
    zx_tl = zx_traj * (qx_tl / qx_traj - denx_tl / denx_traj + 3.0_rd * dmx_tl / dmx_traj)

  end subroutine dm_z_gceop_tl

  subroutine dm_z_gceop_ad(precip_type, qx_ad, denx_ad, dmx_ad, zx_ad)
    !-----------------------------------------------------------
    ! Adjoint (AD) of subroutine dm_z_gceop
    ! N0 is fixed, calcuate lambda from N0
    ! Goddard Cumulus Ensemble scheme used in CWB WRF operation
    ! only difference from wsm6: n0 for snow is a fixed value
    !-------------------------------------------------------------------
    implicit none
    character (len=*), intent(in) :: precip_type
    real(kind=rd), intent(inout)  :: qx_ad          ! mixing ratio of hydrometeor, g/kg
    real(kind=rd), intent(inout)  :: denx_ad        ! wet hydrometeor density, g/cm^3
    real(kind=rd), intent(in)     :: dmx_ad         ! mean diameter, mm
    real(kind=rd), intent(in)     :: zx_ad          ! 6th moment of PSD

    real(kind=rd) :: n0_traj
    real(kind=rd) :: lambda_traj
    real(kind=rd) :: lambda_ad
    real(kind=rd) :: dmx_ad_tmp = zero

    call n0_lambda_gceop(precip_type, density_air_traj, qx_traj, denx_traj, n0_traj, lambda_traj)

    dmx_ad_tmp = dmx_ad + 3.0_rd * zx_traj * zx_ad / dmx_traj
    lambda_ad = - dmx_traj * dmx_ad_tmp / lambda_traj
    qx_ad = qx_ad + ( -0.25_rd*lambda_traj*lambda_ad + zx_traj*zx_ad) / qx_traj
    if (var_den_rat) &
      denx_ad = denx_ad + ( 0.25_rd*lambda_traj*lambda_ad - zx_traj*zx_ad ) / denx_traj

  end subroutine dm_z_gceop_ad

  subroutine dm_z_gcecwb_tl(precip_type, temp_tl, qx_tl, ratx_tl, dmx_tl, zx_tl, denx_tl)
    !----------------------------------------------------------------------------------
    ! Tagent Linear (TL) of subroutine dm_z_gcecwb
    ! Single moment scheme assuming Gamma PSD: N(D) = N0 * D^alpha * exp(-lambda D)
    ! Goddard Cumulus Ensemble scheme modified by CWB WRF
    ! hydrometeor particle density and mean diameter are diagnosed.
    !
    ! T.-C. Tsai, S.-Y. Jiang, J.-P. Chen, H.-L. Huang, Yi.-C. Lo, J.-S. Hong, 2022:
    ! A Gamma-type Single Moment Bulk Microphysics Scheme for Operational Forecasting.
    !-----------------------------------------------------------------------------------
    implicit none
    character (len=*), intent(in) :: precip_type
    real(kind=rd), intent(in)     :: temp_tl          ! temperature, Kelvin
    real(kind=rd), intent(in)     :: qx_tl            ! mixing ratio of hydrometeor, g/kg
    real(kind=rd), intent(in)     :: ratx_tl          ! mass water fraction [0-1]
    real(kind=rd), intent(out)    :: dmx_tl           ! mean diameter, mm
    real(kind=rd), intent(out)    :: zx_tl            ! 6th moment of PSD
    real(kind=rd), intent(out)    :: denx_tl          ! density of wet hydrometeor g/cm^3
    
    real(kind=rd) :: logwx_tl
    real(kind=rd) :: logt_tl

    logwx_tl = qx_tl / qx_traj
    logt_tl = temp_tl / temp_air_traj

    call dm_den_gcecwb_tl(precip_type, logwx_tl, logt_tl, ratx_tl, denx_tl, dmx_tl)

    zx_tl = zx_traj * (qx_tl / qx_traj - denx_tl / denx_traj + 3.0_rd * dmx_tl / dmx_traj)
    
  end subroutine dm_z_gcecwb_tl

  subroutine dm_z_gcecwb_ad(precip_type, temp_ad, qx_ad, ratx_ad, dmx_ad, zx_ad, denx_ad)
    !-------------------------------------------------------------------
    ! Adjoint (AD) of subroutine dm_z_gcecwb
    ! Single moment scheme assuming Gamma PSD: N(D) = N0 * D^alpha * exp(-lambda D)
    ! Goddard Cumulus Ensemble scheme modified by CWB WRF
    ! hydrometeor particle density and mean diameter are diagnosed.
    !
    ! T.-C. Tsai, S.-Y. Jiang, J.-P. Chen, H.-L. Huang, Yi.-C. Lo, J.-S. Hong, 2022:
    ! A Gamma-type Single Moment Bulk Microphysics Scheme for Operational Forecasting.
    !-------------------------------------------------------------------
    implicit none
    character (len=*), intent(in)    :: precip_type
    real(kind=rd), intent(inout)     :: temp_ad          ! temperature, Kelvin
    real(kind=rd), intent(inout)     :: qx_ad            ! mixing ratio of hydrometeor, g/kg
    real(kind=rd), intent(inout)     :: ratx_ad          ! mass water fraction [0-1]
    real(kind=rd), intent(in)        :: dmx_ad           ! mean diameter, mm
    real(kind=rd), intent(in)        :: zx_ad            ! 6th moment of PSD
    real(kind=rd), intent(in)        :: denx_ad          ! density of wet hydrometeor g/cm^3

    real(kind=rd) :: logwx_ad
    real(kind=rd) :: logt_ad
    real(kind=rd) :: dmx_ad_tmp = zero
    real(kind=rd) :: denx_ad_tmp = zero

    dmx_ad_tmp  = dmx_ad + 3.0_rd * zx_traj * zx_ad / dmx_traj
    if (var_den_rat) &
      denx_ad_tmp = denx_ad - zx_traj * zx_ad / denx_traj
    qx_ad       = qx_ad + zx_traj * zx_ad / qx_traj

    call dm_den_gcecwb_ad(precip_type, logwx_ad, logt_ad, ratx_ad, denx_ad_tmp, dmx_ad_tmp)

    qx_ad   = qx_ad + logwx_ad / qx_traj
    temp_ad = temp_ad + logt_ad / temp_air_traj

  end subroutine dm_z_gcecwb_ad
  
  subroutine dm_den_gcecwb_tl(precip_type, logwx_tl, logt_tl, ratx_tl, denx_tl, dmx_tl)
    !-------------------------------------------------------------------
    ! Tagent Linear (TL) of subroutine dm_den_gcecwb
    !-------------------------------------------------------------------
    implicit none
    character (len=*), intent(in) :: precip_type
    real(kind=rd), intent(in)     :: logwx_tl         ! ln(water content), kg/m^3
    real(kind=rd), intent(in)     :: logt_tl          ! ln(temperature), Kelvin
    real(kind=rd), intent(in)     :: ratx_tl          ! mass water fraction [0-1]
    real(kind=rd), intent(out)    :: denx_tl          ! density of wet hydrometeor g/cm^3
    real(kind=rd), intent(out)    :: dmx_tl           ! mean diameter, mm

    real(kind=rd) :: logwx_traj
    real(kind=rd) :: logt_traj
    real(kind=rd) :: tmp1, tmp2

    logwx_traj = log(wx_traj*0.001_rd) ! g/m3 -> kg/m3
    logt_traj = log(temp_air_traj)

    select case (precip_type)
      case ('rain')
        denx_tl = 0.0_rd ! g/cm^3
        dmx_tl = dmx_traj * (6.888453_rd*logt_tl + 0.25_rd*logwx_tl &
                          - 2.0_rd*0.3132_rd*logt_traj*logt_tl + 2.0_rd*1.43e-18_rd*logwx_traj*logwx_tl &
                          - 5.27932e-16_rd*(logt_traj*logwx_tl + logwx_traj*logt_tl)) ! mm
      case ('snow', 'pure snow')
        denx_tl = denx_traj * (-17.124_rd*logt_tl + 0.84975_rd*logwx_tl &
                          + 2.0_rd*1.069132_rd*logt_traj*logt_tl - 2.0_rd*0.00614_rd*logwx_traj*logwx_tl &
                          - 0.208331_rd*(logt_traj*logwx_tl + logwx_traj*logt_tl)) ! g/cm^3
        dmx_tl = dmx_traj * (26.74956_rd*logt_tl - 1.327406_rd*logwx_tl &
                          - 2.0_rd*1.6701_rd*logt_traj*logt_tl + 2.0_rd*0.009598_rd*logwx_traj*logwx_tl &
                          + 0.32543632_rd*(logt_traj*logwx_tl + logwx_traj*logt_tl)) ! mm
      case ('graupel', 'pure graupel')
        tmp1 = 0.001_rd * exp(0.056433_rd*logwx_traj+6.608048_rd) ! g/cm3
        tmp2 = max(0.25_rd, tmp1) ! g/cm3
        denx_tl = (1.0_rd - tmp2) * ratx_tl ! g/cm3
        if (tmp1 > 0.25_rd) then 
          denx_tl = denx_tl + (1.0_rd - ratx_traj) * tmp1 * 0.056433_rd * logwx_tl ! g/cm^3
        end if

        dmx_tl = dmx_traj * (16.08519_rd*logt_tl - 2.530305_rd*logwx_tl &
                          - 2.0_rd*0.19094_rd*logt_traj*logt_tl + 2.0_rd*0.005846_rd*logwx_traj*logwx_tl &
                          + 0.51917156_rd*(logt_traj*logwx_tl + logwx_traj*logt_tl)) ! in mm
      case default
        write(*,*) 'Unknow precipitation type: ', precip_type
        stop
    end select

  end subroutine dm_den_gcecwb_tl

  subroutine dm_den_gcecwb_ad(precip_type, logwx_ad, logt_ad, ratx_ad, denx_ad, dmx_ad)
    !-------------------------------------------------------------------
    ! Adjoint (AD) of subroutine dm_den_gcecwb
    !-------------------------------------------------------------------
    implicit none
    character (len=*), intent(in) :: precip_type
    real(kind=rd), intent(inout)     :: logwx_ad         ! ln(water content), kg/m^3
    real(kind=rd), intent(inout)     :: logt_ad          ! ln(temperature), Kelvin
    real(kind=rd), intent(inout)     :: ratx_ad          ! mass water fraction [0-1]
    real(kind=rd), intent(in)        :: denx_ad          ! density of wet hydrometeor g/cm^3
    real(kind=rd), intent(in)        :: dmx_ad           ! mean diameter, mm

    real(kind=rd) :: logwx_traj
    real(kind=rd) :: logt_traj
    real(kind=rd) :: tmp1, tmp2

    logwx_traj = log(wx_traj*0.001_rd) ! g/m3 -> kg/m3
    logt_traj = log(temp_air_traj)

    select case (precip_type)
      case ('rain')
        logt_ad  =  logt_ad +  dmx_traj * (6.888453_rd  - 2.0_rd*0.3132_rd*logt_traj    - 5.27932e-16_rd*logwx_traj) * dmx_ad
        logwx_ad = logwx_ad +  dmx_traj * (0.25_rd      + 2.0_rd*1.43e-18_rd*logwx_traj - 5.27932e-16_rd*logt_traj ) * dmx_ad 
      case ('snow', 'pure snow')
        logt_ad  =  logt_ad + denx_traj * (-17.124_rd   + 2.0_rd*1.069132_rd*logt_traj  - 0.208331_rd*logwx_traj   ) * denx_ad &
                            +  dmx_traj * (26.74956_rd  - 2.0_rd*1.6701_rd*logt_traj    + 0.32543632_rd*logwx_traj ) * dmx_ad 
        logwx_ad = logwx_ad + denx_traj * (0.84975_rd   - 2.0_rd*0.00614_rd*logwx_traj  - 0.208331_rd*logt_traj    ) * denx_ad &
                            +  dmx_traj * (-1.327406_rd + 2.0_rd*0.009598_rd*logwx_traj + 0.32543632_rd*logt_traj  ) * dmx_ad
      case ('graupel', 'pure graupel')
        logt_ad  = logt_ad  +  dmx_traj * (16.08519_rd  - 2.0_rd*0.19094_rd*logt_traj   + 0.51917156_rd*logwx_traj  ) * dmx_ad
        logwx_ad = logwx_ad +  dmx_traj * (-2.530305_rd + 2.0_rd*0.005846_rd*logwx_traj + 0.51917156_rd*logt_traj   ) * dmx_ad

        tmp1 = 0.001_rd * exp(0.056433_rd*logwx_traj+6.608048_rd)
        tmp2 = max(0.25_rd, tmp1)
        if (var_den_rat) &
          ratx_ad  = ratx_ad  + (1.0_rd - tmp2) * denx_ad
        if (tmp1 > 0.25_rd) then 
          logwx_ad = logwx_ad + (1.0_rd - ratx_traj) * tmp1 * 0.056433_rd * denx_ad
        end if
      case default
        write(*,*) 'Unknow precipitation type: ', precip_type
        stop
    end select

  end subroutine dm_den_gcecwb_ad

  subroutine dualpol_op_rain_tl(a, wx_tl, dmx_tl, zhx_tl, zdrx_tl, kdpx_tl, phvx_tl)
    !---------------------------------------------------------------
    ! Tangent Linear code of Eqs. (13) - (16) of Zhang et al., 2021
    !---------------------------------------------------------------
    implicit none
    real(kind=8), intent(in)   :: a(1:4,0:4)    ! coefs for rain
    real(kind=8), intent(in)   :: wx_tl         ! density_air * mixing ratio , g/m^3
    real(kind=8), intent(in)   :: dmx_tl        ! mass/volume-weighted mean diameter, mm
    real(kind=8), intent(out)  :: zhx_tl        ! horizontal reflectivity, mm^6m-3
    real(kind=8), intent(out)  :: zdrx_tl       ! differential reflectivity, -
    real(kind=8), intent(out)  :: kdpx_tl       ! specific differential phase, degree/km^-1
    real(kind=8), intent(out)  :: phvx_tl       ! co-polar correlation coefficient, 0-1
    
    zhx_tl = zhx_traj / wx_traj * wx_tl + 2.0_rd * wx_traj * &
      (a(1,0) + a(1,1)*dmx_traj + a(1,2)*dmx_traj**2 + a(1,3)*dmx_traj**3 + a(1,4)*dmx_traj**4) * &
      (a(1,1) + 2.0_rd*a(1,2)*dmx_traj + 3.0_rd*a(1,3)*dmx_traj**2 + 4.0_rd*a(1,4)*dmx_traj**3) * dmx_tl
    
    zdrx_tl = (a(2,1) + 2.0_rd*a(2,2)*dmx_traj + 3.0_rd*a(2,3)*dmx_traj**2 + 4.0_rd*a(2,4)*dmx_traj**3) * dmx_tl

    kdpx_tl = kdpx_traj / wx_traj * wx_tl + wx_traj * &
      (a(3,1) + 2.0_rd*a(3,2)*dmx_traj + 3.0_rd*a(3,3)*dmx_traj**2 + 4.0_rd*a(3,4)*dmx_traj**3) * dmx_tl
    
    phvx_tl = (a(4,1) + 2.0_rd*a(4,2)*dmx_traj + 3.0_rd*a(4,3)*dmx_traj**2 + 4.0_rd*a(4,4)*dmx_traj**3) * dmx_tl

  end subroutine dualpol_op_rain_tl

  subroutine dualpol_op_rain_ad(a, wx_ad, dmx_ad, zhx_ad, zdrx_ad, kdpx_ad, phvx_ad)
    !---------------------------------------------------------------
    ! Adjoint code of Eqs. (13) - (16) of Zhang et al., 2021
    !---------------------------------------------------------------
    implicit none
    real(kind=8), intent(in)   :: a(1:4,0:4)    ! coefs for rain
    real(kind=8), intent(inout)   :: wx_ad         ! density_air * mixing ratio , g/m^3
    real(kind=8), intent(inout)   :: dmx_ad        ! mass/volume-weighted mean diameter, mm
    real(kind=8), intent(in)      :: zhx_ad        ! horizontal reflectivity, mm6 m-3
    real(kind=8), intent(in)      :: zdrx_ad       ! differential reflectivity, -
    real(kind=8), intent(in)      :: kdpx_ad       ! specific differential phase, degree/km^-1
    real(kind=8), intent(in)      :: phvx_ad       ! co-polar correlation coefficient, 0-1

    wx_ad = wx_ad + &
      zhx_traj / wx_traj * zhx_ad + &
      kdpx_traj / wx_traj * kdpx_ad
    
    dmx_ad = dmx_ad + &
      2.0_rd * wx_traj * &
      (a(1,0) + a(1,1)*dmx_traj + a(1,2)*dmx_traj**2 + a(1,3)*dmx_traj**3 + a(1,4)*dmx_traj**4) * &
      (a(1,1) + 2.0_rd*a(1,2)*dmx_traj + 3.0_rd*a(1,3)*dmx_traj**2 + 4.0_rd*a(1,4)*dmx_traj**3) * zhx_ad  + &
      (a(2,1) + 2.0_rd*a(2,2)*dmx_traj + 3.0_rd*a(2,3)*dmx_traj**2 + 4.0_rd*a(2,4)*dmx_traj**3) * zdrx_ad + &
      wx_traj * &
      (a(3,1) + 2.0_rd*a(3,2)*dmx_traj + 3.0_rd*a(3,3)*dmx_traj**2 + 4.0_rd*a(3,4)*dmx_traj**3) * kdpx_ad + &
      (a(4,1) + 2.0_rd*a(4,2)*dmx_traj + 3.0_rd*a(4,3)*dmx_traj**2 + 4.0_rd*a(4,4)*dmx_traj**3) * phvx_ad

  end subroutine dualpol_op_rain_ad

  subroutine dualpol_op_icephase_tl(zx_tl, wx_tl, denx_tl, dmx_tl, ax_tl, &
    zhx_tl, zdrx_tl, kdpx_tl, phvx_tl)
    !---------------------------------------------------------------------------------------
    ! Tangent Linear module of Eqs. (17) - (20) of Zhang et al., 2021 for snow/graupel/hail
    !---------------------------------------------------------------------------------------
    implicit none
    real(kind=rd), intent(in)   :: zx_tl          ! 6th moment of PSD, Eq. (6)
    real(kind=rd), intent(in)   :: wx_tl          ! density_air * mixing ratio , g/m^3
    real(kind=rd), intent(in)   :: denx_tl        ! density of melting snow/graupel/hail, g/cm^3
    real(kind=rd), intent(in)   :: dmx_tl         ! mass/volume-weighted mean diameter, mm
    real(kind=rd), intent(in)   :: ax_tl(0:12)    ! coef

    real(kind=rd), intent(out)  :: zhx_tl         ! horizontal reflectivity, mm6 m-3
    real(kind=rd), intent(out)  :: zdrx_tl        ! differential reflectivity, -
    real(kind=rd), intent(out)  :: kdpx_tl        ! specific differential phase, degree/km^-1
    real(kind=rd), intent(out)  :: phvx_tl        ! co-polar correlation coefficient, 0-1  

    zhx_tl = (zhx_traj / zx_traj) * zx_tl + 2.0_rd * zx_traj * &
    (ax_traj(0) + ax_traj(1)*dmx_traj + ax_traj(2)*dmx_traj**2 + ax_traj(3)*dmx_traj**3) * &
    ((ax_traj(1) + 2.0_rd*ax_traj(2)*dmx_traj + 3.0_rd*ax_traj(3)*dmx_traj**2) * dmx_tl + &
    ax_tl(0) + dmx_traj * ax_tl(1) + dmx_traj**2 * ax_tl(2) + dmx_traj**3 * ax_tl(3))

    zdrx_tl = (ax_traj(5) + 2.0_rd*ax_traj(6)*dmx_traj) * dmx_tl + &
    ax_tl(4) + dmx_traj * ax_tl(5) + dmx_traj**2 * ax_tl(6)

    kdpx_tl = (kdpx_traj / wx_traj) * wx_tl - (kdpx_traj / denx_traj) * denx_tl + &
    (wx_traj / denx_traj) * ((ax_traj(8) + 2.0_rd*ax_traj(9)*dmx_traj) * dmx_tl + &
    ax_tl(7) + dmx_traj * ax_tl(8) + dmx_traj**2 * ax_tl(9))

    phvx_tl = (ax_traj(11) + 2.0_rd*ax_traj(12)*dmx_traj) * dmx_tl + &
    ax_tl(10) + dmx_traj * ax_tl(11) + dmx_traj**2 * ax_tl(12)
  end subroutine dualpol_op_icephase_tl

  subroutine dualpol_op_icephase_ad(zx_ad, wx_ad, denx_ad, dmx_ad, ax_ad, &
    zhx_ad, zdrx_ad, kdpx_ad, phvx_ad)
    !---------------------------------------------------------------------------------------
    ! Adjoint module of Eqs. (17) - (20) of Zhang et al., 2021 for snow/graupel/hail
    !---------------------------------------------------------------------------------------
    implicit none
    real(kind=rd), intent(inout)   :: zx_ad       ! 6th moment of PSD, Eq. (6)
    real(kind=rd), intent(inout)   :: wx_ad       ! density_air * mixing ratio , g/m^3
    real(kind=rd), intent(inout)   :: denx_ad     ! density of melting snow/graupel/hail, g/cm^3
    real(kind=rd), intent(inout)   :: dmx_ad      ! mass/volume-weighted mean diameter, mm
    real(kind=rd), intent(inout)   :: ax_ad(0:12) ! coef

    real(kind=rd), intent(in)      :: zhx_ad      ! horizontal reflectivity, mm6 m-3
    real(kind=rd), intent(in)      :: zdrx_ad     ! differential reflectivity, -
    real(kind=rd), intent(in)      :: kdpx_ad     ! specific differential phase, degree/km^-1
    real(kind=rd), intent(in)      :: phvx_ad     ! co-polar correlation coefficient, 0-1  

    real(kind=rd) :: tmp1, tmp2

    ! intermediate variables used by adjoint calculation of multiple variables 
    tmp1 = 2.0_rd * zx_traj * (ax_traj(0) + ax_traj(1)*dmx_traj + ax_traj(2)*dmx_traj**2 + ax_traj(3)*dmx_traj**3)
    tmp2 = wx_traj / denx_traj

    if (var_den_rat) then
      ax_ad(0) = ax_ad(0) + tmp1 * zhx_ad
      ax_ad(1) = ax_ad(1) + tmp1 * dmx_traj * zhx_ad
      ax_ad(2) = ax_ad(2) + tmp1 * dmx_traj**2 * zhx_ad
      ax_ad(3) = ax_ad(3) + tmp1 * dmx_traj**3 * zhx_ad
      ax_ad(4) = ax_ad(4) + zdrx_ad
      ax_ad(5) = ax_ad(5) + dmx_traj * zdrx_ad
      ax_ad(6) = ax_ad(6) + dmx_traj**2 * zdrx_ad
      ax_ad(7) = ax_ad(7) + tmp2 * kdpx_ad
      ax_ad(8) = ax_ad(8) + tmp2 * dmx_traj * kdpx_ad
      ax_ad(9) = ax_ad(9) + tmp2 * dmx_traj**2 * kdpx_ad
      ax_ad(10) = ax_ad(10) + phvx_ad
      ax_ad(11) = ax_ad(11) + dmx_traj * phvx_ad
      ax_ad(12) = ax_ad(12) + dmx_traj**2 * phvx_ad
    endif

    zx_ad = zx_ad + (zhx_traj / zx_traj) * zhx_ad
    wx_ad = wx_ad + (kdpx_traj / wx_traj) * kdpx_ad
    if (var_den_rat) &
      denx_ad = denx_ad - (kdpx_traj / denx_traj) * kdpx_ad
    dmx_ad = dmx_ad + &
    tmp1 * (ax_traj(1) + 2.0_rd*ax_traj(2)*dmx_traj + 3.0_rd*ax_traj(3)*dmx_traj**2) * zhx_ad + &
    (ax_traj(5) + 2.0_rd*ax_traj(6)*dmx_traj) * zdrx_ad + &
    tmp2 * (ax_traj(8) + 2.0_rd*ax_traj(9)*dmx_traj) * kdpx_ad + &
    (ax_traj(11) + 2.0_rd*ax_traj(12)*dmx_traj) * phvx_ad
  end subroutine dualpol_op_icephase_ad

  subroutine dualpol_op_total_tl(zhr_tl,  zhs_tl,  zhg_tl, &
                                zdrr_tl, zdrs_tl, zdrg_tl, &
                                kdpr_tl, kdps_tl, kdpg_tl, &
                                phvr_tl, phvs_tl, phvg_tl, &
                                zh_tl, zdr_tl, kdp_tl, phv_tl, &
                                zhps_tl , zhpg_tl , &
                                zdrps_tl, zdrpg_tl, &
                                kdpps_tl, kdppg_tl, &
                                phvps_tl, phvpg_tl, &
                                zhh_tl,   zdrh_tl,  &
                                kdph_tl,  phvh_tl,  &
                                zhph_tl,  zdrph_tl, &
                                kdpph_tl, phvph_tl)
    !---------------------------------------------------------------
    ! Tangeng Linear of computing total zh, zdr, kdp, phv
    ! Eqs. (22) - (25) of Zhang et al., 2021
    !---------------------------------------------------------------  
    implicit none
    real(kind=rd), intent(in)      :: zhr_tl, zhs_tl, zhg_tl       ! horizontal reflectivity, mm6 m-3
    real(kind=rd), intent(in)      :: zdrr_tl, zdrs_tl, zdrg_tl     ! differential reflectivity, -
    real(kind=rd), intent(in)      :: kdpr_tl, kdps_tl, kdpg_tl    ! specific differential phase, degree/km^-1
    real(kind=rd), intent(in)      :: phvr_tl, phvs_tl, phvg_tl    ! co-polar correlation coefficient, 0-1
    real(kind=rd), intent(in)      :: zhps_tl , zhpg_tl
    real(kind=rd), intent(in)      :: zdrps_tl, zdrpg_tl
    real(kind=rd), intent(in)      :: kdpps_tl, kdppg_tl
    real(kind=rd), intent(in)      :: phvps_tl, phvpg_tl
    ! For Hail
    real(kind=rd), intent(in), optional :: zhh_tl , zhph_tl
    real(kind=rd), intent(in), optional :: zdrh_tl, zdrph_tl
    real(kind=rd), intent(in), optional :: kdph_tl, kdpph_tl
    real(kind=rd), intent(in), optional :: phvh_tl, phvph_tl

    real(kind=rd), intent(out)     :: zh_tl  ! horizontal reflectivity, mm6 m-3
    real(kind=rd), intent(out)     :: zdr_tl ! differential reflectivity, -
    real(kind=rd), intent(out)     :: kdp_tl ! specific differential phase, degree/km^-1
    real(kind=rd), intent(out)     :: phv_tl ! co-polar correlation coefficient, 0-1

    real(kind=rd) :: zbarr_tl, zbars_tl, zbarg_tl, zbarps_tl, zbarpg_tl, zbarh_tl, zbarph_tl
    real(kind=rd) :: zvr_tl, zvs_tl, zvg_tl, zvps_tl, zvpg_tl, zvh_tl, zvph_tl
    real(kind=rd) :: denom_phv, tmp
    real(kind=rd) :: zv_tl, zbar_tl

    logical :: include_hail = .False.

    if (.not. rain_populated .or. .not. snow_populated .or. .not. graupel_populated) then
      write(*,*) 'Precipitation trajactories not populated, abort...'
      stop
    endif

    if (present(zhh_tl)) include_hail = .True.

    denom_phv = zbarr_traj + zbars_traj + zbarg_traj + zbarps_traj + zbarpg_traj
    if (include_hail) denom_phv = denom_phv + zbarh_traj + zbarph_traj

    call zbarzv_tl(   'rain', zhr_tl, zdrr_tl, zbarr_tl, zvr_tl)
    call zbarzv_tl(   'snow', zhs_tl, zdrs_tl, zbars_tl, zvs_tl)
    call zbarzv_tl('graupel', zhg_tl, zdrg_tl, zbarg_tl, zvg_tl)
    call zbarzv_tl('pure snow', zhps_tl, zdrps_tl, zbarps_tl, zvps_tl)
    call zbarzv_tl('pure graupel', zhpg_tl, zdrpg_tl, zbarpg_tl, zvpg_tl)
    if (include_hail) call zbarzv_tl('hail', zhh_tl, zdrh_tl, zbarh_tl, zvh_tl)
    if (include_hail) call zbarzv_tl('pure hail', zhph_tl, zdrph_tl, zbarph_tl, zvph_tl)

    zh_tl  = zhr_tl + zhs_tl + zhg_tl + zh_tl + zhps_tl + zhpg_tl
    if (include_hail) zh_tl = zh_tl + zhh_tl + zhph_tl
    zv_tl  = zvr_tl + zvs_tl + zvg_tl +  zv_tl + zvps_tl + zvpg_tl
    if (include_hail) zv_tl = zv_tl + zvh_tl + zvph_tl

    zdr_tl = (zdr_traj / zh_traj) * zh_tl - &
             (zdr_traj * zdr_traj / zh_traj) * zv_tl

    kdp_tl = kdpr_tl + kdps_tl + kdpg_tl + kdpps_tl + kdppg_tl
    if (include_hail) kdp_tl = kdp_tl + kdph_tl + kdpph_tl

    zbar_tl = zbarr_tl + zbars_tl + zbarg_tl + zbarps_tl + zbarpg_tl
    if (include_hail) zbar_tl = zbar_tl + zbarh_tl + zbarph_tl

    tmp = zbarr_traj*phvr_tl + zbars_traj*phvs_tl + zbarg_traj*phvg_tl &
      + phvr_traj*zbarr_tl + phvs_traj*zbars_tl + phvg_traj*zbarg_tl + &
      zbarps_traj*phvps_tl + zbarpg_traj*phvpg_tl + phvps_traj*zbarps_tl + phvpg_traj*zbarpg_tl
    if (include_hail) tmp = tmp + zbarh_traj*phvh_tl + zbarph_traj*phvph_tl &
                            + phvh_traj*zbarh_tl + phvph_traj*zbarph_tl
    
    phv_tl = (1.0_rd / denom_phv) * tmp  - (phv_traj / denom_phv) * zbar_tl

  end subroutine dualpol_op_total_tl
  
  subroutine dualpol_op_total_ad(zhr_ad,  zhs_ad,  zhg_ad, &
                                zdrr_ad, zdrs_ad, zdrg_ad, &
                                kdpr_ad, kdps_ad, kdpg_ad, &
                                phvr_ad, phvs_ad, phvg_ad, &
                                zh_ad, zdr_ad, kdp_ad, phv_ad, &
                                zhps_ad, zhpg_ad, &
                                zdrps_ad, zdrpg_ad, &
                                kdpps_ad, kdppg_ad, &
                                phvps_ad, phvpg_ad, &
                                zhh_ad,   zdrh_ad,  &
                                kdph_ad,  phvh_ad,  &
                                zhph_ad,  zdrph_ad, &
                                kdpph_ad, phvph_ad )
    !---------------------------------------------------------------
    ! Adjoint of computing total zh, zdr, kdp, phv
    ! Eqs. (22) - (25) of Zhang et al., 2021
    !---------------------------------------------------------------  
    implicit none
    real(kind=rd), intent(inout)   :: zhr_ad, zhs_ad, zhg_ad       ! horizontal reflectivity, mm6 m-3
    real(kind=rd), intent(inout)   :: zdrr_ad, zdrs_ad, zdrg_ad     ! differential reflectivity, -
    real(kind=rd), intent(inout)   :: kdpr_ad, kdps_ad, kdpg_ad    ! specific differential phase, degree/km^-1
    real(kind=rd), intent(inout)   :: phvr_ad, phvs_ad, phvg_ad    ! co-polar correlation coefficient, 0-1
    ! For pure hydrometeor
    real(kind=rd), intent(inout)   :: zhps_ad,  zhpg_ad
    real(kind=rd), intent(inout)   :: zdrps_ad, zdrpg_ad
    real(kind=rd), intent(inout)   :: kdpps_ad, kdppg_ad
    real(kind=rd), intent(inout)   :: phvps_ad, phvpg_ad
    ! For Hail
    real(kind=rd), intent(inout), optional :: zhh_ad, zhph_ad
    real(kind=rd), intent(inout), optional :: zdrh_ad, zdrph_ad
    real(kind=rd), intent(inout), optional :: kdph_ad, kdpph_ad
    real(kind=rd), intent(inout), optional :: phvh_ad, phvph_ad

    real(kind=rd), intent(in)      :: zh_ad  ! horizontal reflectivity, mm6 m-3
    real(kind=rd), intent(in)      :: zdr_ad ! differential reflectivity, -
    real(kind=rd), intent(in)      :: kdp_ad ! specific differential phase, degree/km^-1
    real(kind=rd), intent(in)      :: phv_ad ! co-polar correlation coefficient, 0-1

    real(kind=rd) :: zbarr_ad, zbars_ad, zbarg_ad, zbarps_ad, zbarpg_ad, zbarh_ad, zbarph_ad
    real(kind=rd) :: zvr_ad, zvs_ad, zvg_ad, zvps_ad, zvpg_ad, zvh_ad, zvph_ad
    real(kind=rd) :: denom_phv
    real(kind=rd) :: tmp1, tmp2, tmp3, tmp4

    logical :: include_hail = .False.

    if (.not. rain_populated .or. .not. snow_populated .or. .not. graupel_populated .or. &
        .not. pure_snow_populated .or. .not. pure_graupel_populated) then
      write(*,*) 'Precipitation trajactories not populated, abort...'
      stop
    endif

    if (present(zhh_ad) ) then
      include_hail = .True.
      if (.not. hail_populated .or. .not. pure_hail_populated) then
        write(*,*) 'Precipitation trajactories not populated, abort...'
        stop
      endif
    endif

    denom_phv = zbarr_traj + zbars_traj + zbarg_traj + zbarps_traj + zbarpg_traj
    if (include_hail) denom_phv = denom_phv + zbarh_traj + zbarph_traj

    tmp1 = zdr_traj / zh_traj
    tmp2 = zdr_traj * zdr_traj / zh_traj
    tmp3 = phv_traj / denom_phv
    tmp4 = 1.0_rd / denom_phv
    
    zhr_ad = zhr_ad + zh_ad + tmp1 * zdr_ad
    zhs_ad = zhs_ad + zh_ad + tmp1 * zdr_ad 
    zhg_ad = zhg_ad + zh_ad + tmp1 * zdr_ad
    zhps_ad = zhps_ad + zh_ad + tmp1 * zdr_ad
    zhpg_ad = zhpg_ad + zh_ad + tmp1 * zdr_ad
    if (include_hail) zhh_ad = zhh_ad + zh_ad + tmp1 * zdr_ad
    if (include_hail) zhph_ad = zhph_ad + zh_ad + tmp1 * zdr_ad

    zvr_ad = - tmp2 * zdr_ad
    zvs_ad = - tmp2 * zdr_ad
    zvg_ad = - tmp2 * zdr_ad
    zvps_ad = - tmp2 * zdr_ad
    zvpg_ad = - tmp2 * zdr_ad
    if (include_hail) zvh_ad = - tmp2 * zdr_ad
    if (include_hail) zvph_ad = - tmp2 * zdr_ad

    zbarr_ad = (tmp4 * phvr_traj - tmp3) * phv_ad
    zbars_ad = (tmp4 * phvs_traj - tmp3) * phv_ad
    zbarg_ad = (tmp4 * phvg_traj - tmp3) * phv_ad
    zbarps_ad = (tmp4 * phvps_traj - tmp3) * phv_ad
    zbarpg_ad = (tmp4 * phvpg_traj - tmp3) * phv_ad
    if (include_hail) zbarh_ad = (tmp4 * phvh_traj - tmp3) * phv_ad
    if (include_hail) zbarph_ad = (tmp4 * phvph_traj - tmp3) * phv_ad

    kdpr_ad = kdpr_ad + kdp_ad
    kdps_ad = kdps_ad + kdp_ad
    kdpg_ad = kdpg_ad + kdp_ad
    kdpps_ad = kdpps_ad + kdp_ad
    kdppg_ad = kdppg_ad + kdp_ad
    if (include_hail) kdph_ad = kdph_ad + kdp_ad
    if (include_hail) kdpph_ad = kdpph_ad + kdp_ad

    phvr_ad = phvr_ad + tmp4 * zbarr_traj * phv_ad
    phvs_ad = phvs_ad + tmp4 * zbars_traj * phv_ad
    phvg_ad = phvg_ad + tmp4 * zbarg_traj * phv_ad
    phvps_ad = phvps_ad + tmp4 * zbarps_traj * phv_ad
    phvpg_ad = phvpg_ad + tmp4 * zbarpg_traj * phv_ad
    if (include_hail) phvh_ad = phvh_ad + tmp4 * zbarh_traj * phv_ad
    if (include_hail) phvph_ad = phvph_ad + tmp4 * zbarph_traj * phv_ad

    call zbarzv_ad(   'rain', zhr_ad, zdrr_ad, zbarr_ad, zvr_ad)
    call zbarzv_ad(   'snow', zhs_ad, zdrs_ad, zbars_ad, zvs_ad)
    call zbarzv_ad('graupel', zhg_ad, zdrg_ad, zbarg_ad, zvg_ad)
    call zbarzv_ad(   'pure snow', zhps_ad, zdrps_ad, zbarps_ad, zvps_ad)
    call zbarzv_ad('pure graupel', zhpg_ad, zdrpg_ad, zbarpg_ad, zvpg_ad)
    if (include_hail) call zbarzv_ad(   'hail', zhh_ad, zdrh_ad, zbarh_ad, zvh_ad)
    if (include_hail) call zbarzv_ad('pure hail', zhph_ad, zdrph_ad, zbarph_ad, zvph_ad)
  
  end subroutine dualpol_op_total_ad

  subroutine zbarzv_tl(precip_type, zhx_tl, zdrx_tl, zbarx_tl, zvx_tl)
    implicit none
    character (len=*), intent(in)    :: precip_type
    real(kind=rd), intent(in)        :: zhx_tl         ! horizontal reflectivity, mm6 m-3
    real(kind=rd), intent(in)        :: zdrx_tl        ! differential reflectivity, -
    real(kind=rd), intent(out)       :: zbarx_tl       ! intermediate reflectivity, mm6 m-3
    real(kind=rd), intent(out)       :: zvx_tl         ! vertical reflectivity, mm6 m-3

    call set_traj_precip_type(precip_type)

    zbarx_tl = zbarx_traj * (zhx_tl / zhx_traj - 0.5_rd * zdrx_tl / zdrx_traj)
    zvx_tl   = zvx_traj   * (zhx_tl / zhx_traj -          zdrx_tl / zdrx_traj)

  end subroutine zbarzv_tl

  subroutine zbarzv_ad(precip_type, zhx_ad, zdrx_ad, zbarx_ad, zvx_ad)
    implicit none
    character (len=*), intent(in)       :: precip_type
    real(kind=rd), intent(inout)        :: zhx_ad         ! horizontal reflectivity, mm6 m-3
    real(kind=rd), intent(inout)        :: zdrx_ad        ! differential reflectivity, -
    real(kind=rd), intent(in)           :: zbarx_ad       ! intermediate reflectivity, mm6 m-3
    real(kind=rd), intent(in)           :: zvx_ad         ! vertical reflectivity, mm6 m-3 

    call set_traj_precip_type(precip_type)

    zhx_ad  = zhx_ad  + (         zbarx_traj * zbarx_ad + zvx_traj * zvx_ad) / zhx_traj
    zdrx_ad = zdrx_ad - (0.5_rd * zbarx_traj * zbarx_ad + zvx_traj * zvx_ad) / zdrx_traj
  
  end subroutine zbarzv_ad
  
  subroutine weticephase_density_tl(precip_type, ratx_tl, denx_tl)
    implicit none
    character (len=*), intent(in)   :: precip_type
    real(kind=rd), intent(in)       :: ratx_tl        ! mass water fraction, [-]
    real(kind=rd), intent(out)      :: denx_tl        ! density of wet hydrometeor, g/cm^3

    select case (precip_type)
      case ('snow')
        denx_tl = 2.0_rd * (density_rain - density_snow) * ratx_traj * ratx_tl
      case ('graupel')
        denx_tl = 2.0_rd * (density_rain - density_graupel) * ratx_traj * ratx_tl
      case ('hail')
        denx_tl = 2.0_rd * (density_rain - density_hail) * ratx_traj * ratx_tl
      case default
        write(*,*) 'Unknow precipitation type: ', precip_type
        stop
    end select

  end subroutine weticephase_density_tl

  subroutine weticephase_density_ad(precip_type, ratx_ad, denx_ad)
    implicit none
    character (len=*), intent(in)    :: precip_type
    real(kind=rd), intent(inout)     :: ratx_ad        ! mass water fraction, [-]
    real(kind=rd), intent(in)        :: denx_ad        ! density of wet hydrometeor, g/cm^3

    select case (precip_type)
      case ('snow')
        ratx_ad = ratx_ad + 2.0_rd * (density_rain - density_snow) * ratx_traj * denx_ad
      case ('graupel')
        ratx_ad = ratx_ad + 2.0_rd * (density_rain - density_graupel) * ratx_traj * denx_ad
      case ('hail')
        ratx_ad = ratx_ad + 2.0_rd * (density_rain - density_hail) * ratx_traj * denx_ad
      case default
        write(*,*) 'Unknow precipitation type: ', precip_type
        stop
    end select
  
  end subroutine weticephase_density_ad

  subroutine watercontent_tl(qx_tl, wx_tl)
    implicit none
    real(kind=rd), intent(in)     :: qx_tl
    real(kind=rd), intent(out)    :: wx_tl
    wx_tl = density_air_traj * qx_tl
  end subroutine watercontent_tl

  subroutine watercontent_ad(qx_ad, wx_ad)
    implicit none
    real(kind=rd), intent(inout)     :: qx_ad
    real(kind=rd), intent(in)        :: wx_ad
    qx_ad = qx_ad + density_air_traj * wx_ad
  end subroutine watercontent_ad

  subroutine coef_a_tl(c, ratx_tl, ax_tl)
    !--------------------------------------------------------
    ! Tangent Linear module of Eq. (21) of Zhang et al., 2021
    !--------------------------------------------------------
    implicit none
    real(kind=rd), intent(in)   :: c(0:3)    ! coef c in Table
    real(kind=rd), intent(in)   :: ratx_tl   ! water mass fraction
    real(kind=rd), intent(out)  :: ax_tl     ! coef a in radar variables formulas

    integer :: i

    ax_tl = c(1)
    do i = 2, 3
      ax_tl = ax_tl + i * c(i) * ratx_traj**(i-1)
    end do

    ax_tl = ax_tl * ratx_tl

  end subroutine coef_a_tl

  subroutine coef_a_ad(c, ratx_ad, ax_ad)
    !-------------------------------------------------
    ! Adjoint module of Eq. (21) of Zhang et al., 2021
    !-------------------------------------------------
    implicit none
    real(kind=rd), intent(in)    :: c(0:3)    ! coef c in Table
    real(kind=rd), intent(inout) :: ratx_ad   ! water mass fraction
    real(kind=rd), intent(in)    :: ax_ad     ! coef a in radar variables formulas

    integer :: i
    real(kind=8) :: ratx_ad_tmp

    ratx_ad_tmp = c(1)
    do i = 2, 3
      ratx_ad_tmp = ratx_ad_tmp + i * c(i) * ratx_traj**(i-1)
    end do

    ratx_ad_tmp = ratx_ad_tmp * ax_ad
    ratx_ad = ratx_ad + ratx_ad_tmp
  end subroutine coef_a_ad

  subroutine melting_scheme_zhang24_tl( &
    qr_tl, qs_tl, qg_tl, &
    qms_tl, qmg_tl, qpr_tl, qps_tl, qpg_tl, rats_tl, ratg_tl, &
    ntr_tl, nts_tl, ntg_tl, &
    ntms_tl, ntmg_tl, ntpr_tl, ntps_tl, ntpg_tl, &
    qh_tl, nth_tl, rath_tl, &
    qmh_tl, qph_tl, ntmh_tl, ntph_tl )

    implicit none
    real(kind=rd), intent(in) :: qr_tl, qs_tl, qg_tl
    real(kind=rd), intent(in), optional :: ntr_tl, nts_tl, ntg_tl
    real(kind=rd), intent(in), optional :: qh_tl, nth_tl

    real(kind=rd), intent(out) :: qms_tl, qmg_tl, qpr_tl, qps_tl, qpg_tl
    real(kind=rd), intent(out) :: rats_tl, ratg_tl
    real(kind=rd), intent(out), optional :: ntms_tl, ntmg_tl, ntpr_tl, ntps_tl, ntpg_tl
    real(kind=rd), intent(out), optional :: rath_tl, qmh_tl, qph_tl, ntmh_tl, ntph_tl

    real(kind=rd) :: qmsl_tl, qmgl_tl, qmhl_tl ! (qms*rats)_tl, (qmg*ratg)_tl

    logical      :: dm = .false. ! double moment

    if (present(ntr_tl) .and. present(nts_tl) .and. present(ntg_tl)) then
      dm = .true.
    endif

    qms_tl = 0.5_rd * qms_traj * (qs_tl/qs_traj + qr_tl/qr_traj)
    qmg_tl = 0.5_rd * qmg_traj * (qg_tl/qg_traj + qr_tl/qr_traj)
    if (present(qh_tl)) qmh_tl = 0.5_rd * qmh_traj * (qh_tl/qh_traj + qr_tl/qr_traj) 

    rats_tl = (rats_traj / qr_traj) * ((1.0_rd - rats_traj) * qr_tl - rats_traj * qs_tl)
    ratg_tl = (ratg_traj / qr_traj) * ((1.0_rd - ratg_traj) * qr_tl - ratg_traj * qg_tl)
    if (present(qh_tl)) rath_tl = (rath_traj / qr_traj) * ((1.0_rd - rath_traj) * qr_tl - rath_traj * qh_tl)

    qmsl_tl = qms_traj * rats_tl + qms_tl * rats_traj
    qmgl_tl = qmg_traj * ratg_tl + qmg_tl * ratg_traj
    if (present(qh_tl)) qmhl_tl = qmh_traj * rath_tl + qmh_tl * rath_traj

    qpr_tl = qr_tl - qmsl_tl - qmgl_tl
    qps_tl = qs_tl - qms_tl + qmsl_tl
    qpg_tl = qg_tl - qmg_tl + qmgl_tl
    if (present(qh_tl)) then
      qpr_tl = qpr_tl - qmhl_tl
      qpg_tl = qh_tl - qmh_tl + qmhl_tl
    endif

    if (dm) then
      ntms_tl = 0.5_rd * ntms_traj * (nts_tl / nts_traj + ntr_tl / ntr_traj)
      ntmg_tl = 0.5_rd * ntmg_traj * (ntg_tl / ntg_traj + ntr_tl / ntr_traj)
      ntpr_tl = ntr_tl
      ntps_tl = nts_tl
      ntpg_tl = ntg_tl
      if (present(nth_tl)) then
         ntms_tl = 0.5_rd * ntmh_traj * (nth_tl / nth_traj + ntr_tl / ntr_traj)
         ntph_tl = nth_tl
      endif
    endif

  end subroutine melting_scheme_zhang24_tl

  subroutine melting_scheme_zhang24_ad( &
    qr_ad, qs_ad, qg_ad, &
    qms_ad, qmg_ad, qpr_ad, qps_ad, qpg_ad, rats_ad, ratg_ad, &
    ntr_ad, nts_ad, ntg_ad, &
    ntms_ad, ntmg_ad, ntpr_ad, ntps_ad, ntpg_ad, &
    qh_ad, nth_ad, rath_ad, &
    qmh_ad, qph_ad, ntmh_ad, ntph_ad )
    
    implicit none
    real(kind=rd), intent(inout) :: qr_ad, qs_ad, qg_ad
    real(kind=rd), intent(inout), optional :: ntr_ad, nts_ad, ntg_ad
    real(kind=rd), intent(inout), optional :: qh_ad, nth_ad

    real(kind=rd), intent(in) :: qms_ad, qmg_ad, qpr_ad, qps_ad, qpg_ad
    real(kind=rd), intent(in) :: rats_ad, ratg_ad
    real(kind=rd), intent(in), optional :: ntms_ad, ntmg_ad, ntpr_ad, ntps_ad, ntpg_ad
    real(kind=rd), intent(in), optional :: qmh_ad, qph_ad, rath_ad, ntmh_ad, ntph_ad

    real(kind=rd) :: qmsl_ad, qmgl_ad, qmhl_ad ! (qmg*ratg)_tl, (qmg*ratg)_tl
    real(kind=rd) :: qms_ad_tmp, qmg_ad_tmp, qmh_ad_tmp
    real(kind=rd) :: rats_ad_tmp, ratg_ad_tmp, rath_ad_tmp

    logical      :: dm = .false. ! double moment

    if (present(ntr_ad) .and. present(nts_ad) .and. present(ntg_ad)) then
      dm = .true.
    endif

    ! qpr_tl = qr_tl - qmsl_tl - qmgl_tl
    ! qps_tl = qs_tl - qms_tl + qmsl_tl
    ! qpg_tl = qg_tl - qmg_tl + qmgl_tl

    qr_ad = qr_ad + qpr_ad
    qs_ad = qs_ad + qps_ad 
    qg_ad = qg_ad + qpg_ad
    if (present(qh_ad)) qh_ad = qh_ad + qph_ad

    qms_ad_tmp = qms_ad - qps_ad
    qmg_ad_tmp = qmg_ad - qpg_ad
    if (present(qh_ad)) qmh_ad_tmp = qmh_ad - qph_ad

    qmsl_ad = qps_ad - qpr_ad
    qmgl_ad = qpg_ad - qpr_ad
    if (present(qh_ad)) qmhl_ad = qph_ad - qpr_ad

    ! qmsl_tl = qms_traj * rats_tl + qms_tl * rats_traj
    ! qmgl_tl = qmg_traj * ratg_tl + qmg_tl * ratg_traj

    qms_ad_tmp = qms_ad_tmp + rats_traj * qmsl_ad
    qmg_ad_tmp = qmg_ad_tmp + ratg_traj * qmgl_ad
    if (present(qh_ad)) qmh_ad_tmp = qmh_ad_tmp + rath_traj * qmhl_ad

    rats_ad_tmp = rats_ad + qms_traj * qmsl_ad
    ratg_ad_tmp = ratg_ad + qmg_traj * qmgl_ad
    if (present(qh_ad)) rath_ad_tmp = rath_ad + qmh_traj * qmhl_ad

    ! rats_tl = (rats_traj / qr_traj) * ((1.0_rd - rats_traj) * qr_tl - rats_traj * qs_tl)
    ! ratg_tl = (ratg_traj / qr_traj) * ((1.0_rd - ratg_traj) * qr_tl - ratg_traj * qg_tl)

    qr_ad = qr_ad + (rats_traj / qr_traj) * (1.0_rd - rats_traj) * rats_ad_tmp + &
                    (ratg_traj / qr_traj) * (1.0_rd - ratg_traj) * ratg_ad_tmp
    
    qs_ad = qs_ad - (rats_traj / qr_traj) * rats_traj * rats_ad_tmp
    qg_ad = qg_ad - (ratg_traj / qr_traj) * ratg_traj * ratg_ad_tmp
 
    if (present(qh_ad)) then
       qr_ad = qr_ad + (rath_traj / qr_traj) * (1.0_rd - rath_traj) * rath_ad_tmp
       qh_ad = qh_ad - (rath_traj / qr_traj) * rath_traj * rath_ad_tmp
    endif
    ! qms_tl = 0.5_rd * qms_traj * (qs_tl/qs_traj + qr_tl/qr_traj)
    ! qmg_tl = 0.5_rd * qmg_traj * (qg_tl/qg_traj + qr_tl/qr_traj)

    qr_ad = qr_ad + 0.5_rd * (qms_traj / qr_traj * qms_ad_tmp + qmg_traj / qr_traj * qmg_ad_tmp)
    qs_ad = qs_ad + 0.5_rd * qms_traj / qs_traj * qms_ad_tmp
    qg_ad = qg_ad + 0.5_rd * qmg_traj / qg_traj * qmg_ad_tmp
 
    if (present(qh_ad)) then
       qr_ad = qr_ad + 0.5_rd * (qmh_traj / qr_traj * qmh_ad_tmp)
       qh_ad = qh_ad + 0.5_rd * qmh_traj / qh_traj * qmh_ad_tmp
    end if
 
    if (dm) then
      ntr_ad = ntr_ad + ntpr_ad 
      nts_ad = nts_ad + ntps_ad
      ntg_ad = ntg_ad + ntpg_ad
      if (present(nth_ad)) nth_ad = nth_ad + ntph_ad
      ntr_ad = ntr_ad + 0.5_rd / ntr_traj * (ntms_traj*ntms_ad + ntmg_traj*ntmg_ad)
      if (present(nth_ad)) ntr_ad = ntr_ad + 0.5_rd / ntr_traj * ntmh_traj*ntmh_ad
      nts_ad = nts_ad + 0.5_rd / nts_traj * ntms_traj*ntms_ad
      ntg_ad = ntg_ad + 0.5_rd / ntg_traj * ntmg_traj*ntmg_ad
      if (present(nth_ad)) nth_ad = nth_ad + 0.5_rd / nth_traj * ntmh_traj*ntmh_ad
    endif

  end subroutine melting_scheme_zhang24_ad

end module zhang21_tlad_mod
