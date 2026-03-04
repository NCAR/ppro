module tcwa2_forward_mod
!--------------------------
! TCWA2 polarimetric radar operator
! 
! This module contains the TCWA2-specific dual-polarization radar forward
! operator subroutines. These assume a gamma PSD and provide empirically fitted
! analytic parameterizations as functions of gamma-PSD parameters (λ, α)
! (or equivalently bulk moments qx, Ntx with derived size metrics such as Dx).
! The parameterizations are derived by fitting results from offline
! particle-size-bin scattering calculations integrated over the PSD under the
! Rayleigh approximation. Coefficients are embedded in the source (no external
! lookup/coeff files required at runtime).
!
! Reference: Tsai, T.-C., J.-P. Chen, Z. Liu, S.-Y. Jiang, R. Kong, Y.-J. Wu, 
!            J. Ban, L.-F. Hsiao, Y.-S. Tang, P.-L. Chang, and J.-S. Hong, 2026:
!            Development of the TCWA2 Bulk Cloud Microphysics Scheme and Its 
!            Integration with a Dual-Polarization Radar Operator for Forecasting 
!            Applications. to be submitted.
! 
! Author: Tzu-Chin Tsai (CWA, original implementation)
! Integrated into PPRO library: Rong Kong (NCAR/MMM), 2025
!--------------------------------------
  implicit none

  private
  
  ! Public interfaces
  public :: dualpol_op_rain_tcwa2
  public :: dualpol_op_ice_tcwa2
  public :: dualpol_op_snow_tcwa2
  public :: dualpol_op_graup_tcwa2
  public :: GAMLN

contains

!---------------------------------------TCWA2------------------------------------------
!--------------------------------------------------------------------------------------
      subroutine dualpol_op_rain_tcwa2(iband,rho,qc0,qr0,nr0,zh,zv,kdp)
!--------------------------------------------------------------------------------------
      implicit none
      integer, intent(in) :: iband       ! radar scan band (wavelength)
      real(kind=8), intent(in)  :: rho   ! air density, kg/m^3
      real(kind=8), intent(in)  :: qc0   ! mass mixing ratio, g/kg
      real(kind=8), intent(in)  :: qr0   ! mass mixing ratio, g/kg
      real(kind=8), intent(in)  :: nr0   ! number concentration, 1/m^3
      real(kind=8), intent(out) :: zh    ! horizontal reflectivity, mm^6/m^3
      real(kind=8), intent(out) :: zv    ! vertical reflectivity, mm^6/m^3
      real(kind=8), intent(out) :: kdp   ! specific differential phase, degree/km^-1
      real, parameter :: drmin = 1.E-4, drmax = 6.E-3
      real :: qc,qr,nr,mvrr,afar,lamr,lamrmin,lamrmax,gr1,gr4,mvdr,    &
              llmr,llmr2,llmr3,lar2,lar22,lar23,fdbzr,dbzr,fzdrr,zdr,  &
              fkdpr,dbzwl
      real(kind=8), parameter :: thrd = 1./3.,         zxmin = 1.e-21
      real(kind=8), parameter :: afamin = 0.,          afamax = 30.
      real(kind=8), parameter :: afa1 = -0.35897435,   afa2 = 0.010041278
      real(kind=8), parameter :: afa3 = 1.9349749E-3,  afa4 = -0.067453607
      real(kind=8), parameter :: afa5 = -5.5693923E-3, afa6 = 0.06376154
      real(kind=8), parameter :: afb1 = -6.5100404,    afb2 = 0.039597245
      real(kind=8), parameter :: afb3 = 2.3169628E-3,  afb4 = -0.45856646
      real(kind=8), parameter :: afb5 = -0.04537631,   afb6 = 0.25713479

      zh = zxmin
      zv = zxmin
      kdp = 0.
      qc = qc0/1.e3      ! g/kg -> kg/kg
      qr = qr0/1.e3      ! g/kg -> kg/kg
      nr = nr0/rho       ! 1/m^3 -> 1/kg
      if (iband.eq.1) then
         dbzwl = 0.11    ! for S band
      else
         dbzwl = 0.053   ! for C band
      endif
      if (qr.ge.1.e-14.and.nr.ge.1.e-2) then
         mvrr = MIN(drmax/2.,MAX(drmin/2.,(qr/nr/4.18879E3)**thrd))
         if (qc.ge.1.e-14.or.nr.lt.0.5) then
            afar = MIN(afamax,MAX(afamin,(afa1+afa2*LOG(nr)+afa3*      &
                   LOG(nr)**2.+afa4*LOG(mvrr))/(1.+afa5*LOG(nr)+afa6*  &
                   LOG(mvrr))))
         else
            afar = MIN(afamax,MAX(afamin,(afb1+afb2*LOG(nr)+afb3*      &
                   LOG(nr)**2.+afb4*LOG(mvrr))/(1.+afb5*LOG(nr)+afb6*  &
                   LOG(mvrr))))
         endif
         if ((nr*mvrr**6.).ge.1.E-16) then
            if (nr.ge.3.8e3) then
               afar = afamax
            else
               afar = 0.5*afar
            endif
         endif
         gr1     = GAMLN(afar+1.)
         gr4     = GAMLN(afar+4.)
         lamr    = (EXP(gr4-gr1)*nr*523.5987756/qr)**thrd
         lamrmin = (EXP(gr4-gr1))**thrd/drmax
         lamrmax = (EXP(gr4-gr1))**thrd/drmin
         IF (lamr.LT.lamrmin) THEN
            lamr = lamrmin
            nr = qr*1.909859E-3*EXP(gr1-gr4+3.*LOG(lamr))
         ELSEIF (lamr.GT.lamrmax) THEN
            lamr = lamrmax
            nr = qr*1.909859E-3*EXP(gr1-gr4+3.*LOG(lamr))
         ENDIF
         mvdr = (EXP(gr4-gr1))**thrd/lamr
      if (qr.ge.1.e-6) then
         llmr = LOG(lamr)
         lar2 = LOG(afar+2.)
         llmr2 = llmr*llmr
         llmr3 = llmr2*llmr
         lar22 = lar2*lar2
         lar23 = lar22*lar2
         fdbzr = 1950.2082+29.685692*lar2-17.761215*SQRT(lar2)+        &
                 328.61414*LOG(llmr)**2.-1122.8963*SQRT(llmr)
         zh = nr*10.**((fdbzr-200.)/10.-18.)
         dbzr = 10.*LOG10(1.E18*zh)
         IF (afar.GE.3.) THEN
            fzdrr = EXP(27.639873+3.7680268*lar2-4.3382556*llmr+       &
                    0.36498072*lar22+0.29286255*llmr2-0.58159197*lar2* &
                    llmr+9.968596E-3*lar23-7.1806791E-3*llmr3+         &
                    0.023642573*lar2*llmr2-0.028767891*lar22*llmr)/1.E4
         ELSE
            fzdrr = EXP(10.243103+0.46561642*lar2+1.7270962*llmr+      &
                    5.8974032E-3*lar22-0.34268382*llmr2+0.10193409*    &
                    lar2*llmr-0.11005788*lar23+0.015622039*llmr3-      &
                    0.038623263*lar2*llmr2+0.10413007*lar22*llmr)/1.E4
         ENDIF
         zv = 10.**(dbzr/10.-18.-fzdrr/10.)
         IF (mvdr.LT.7.E-4) THEN
            fkdpr = MIN(1.5/MAX(1.,nr),EXP(14.837602-18.786496/lar2+   &
                    72.211021/lar22-124.16879/lar23+98.39131/lar2**4.- &
                    28.355229/lar2**5.-0.092813834*llmr)/1.E10)
         ELSE
            IF (afar.LE.3.) THEN
               fkdpr = EXP(35.450715+2.4505311*lar2-1.4651261*llmr-    &
                       0.22974058*lar22-0.22139098*llmr2+0.29382343*   &
                       lar2*llmr)/1.E10
            ELSE
               fkdpr = EXP((38.35738+8.8170679*lar2-8.971269*llmr+     &
                       0.49634922*lar22+0.51897478*llmr2-1.0162854*    &
                       lar2*llmr)/(1.+0.16321213*lar2-0.1676006*llmr+  &
                       4.0072333E-3*lar22+4.9701831E-3*llmr2-          &
                       9.1512549E-3*lar2*llmr))/1.E10
            ENDIF
            fkdpr = MIN(10./MAX(1.,nr),fkdpr)
         ENDIF
         kdp = nr*fkdpr*0.11/dbzwl
         IF (zv.GT.zh) zv = 0.9999*zh
         zdr = 10.*LOG10(zh/zv)
      endif
      endif

      end subroutine dualpol_op_rain_tcwa2
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
      subroutine dualpol_op_ice_tcwa2(iband,tk,rho,qc0,qi0,ni0,zh,zv,kdp)
!-------------------------------------------------------------------------------------
      implicit none
      integer, intent(in) :: iband      ! radar scan band (wavelength)
      real(kind=8), intent(in)  :: tk   ! air temperature, K
      real(kind=8), intent(in)  :: rho  ! air density, kg/m^3
      real(kind=8), intent(in)  :: qc0  ! mass mixing ratio, g/kg
      real(kind=8), intent(in)  :: qi0  ! mass mixing ratio, g/kg
      real(kind=8), intent(in)  :: ni0  ! number concentration, 1/m^3
      real(kind=8), intent(out) :: zh   ! horizontal reflectivity, mm^6/m^3
      real(kind=8), intent(out) :: zv   ! vertical reflectivity, mm^6/m^3
      real(kind=8), intent(out) :: kdp  ! specific differential phase, degree/km^-1
      integer :: hid
      real, parameter :: dimin = 6.E-6, dimax = 5.E-3
      real :: qc,qi,ni,dbzwl,rhoi,mvri,afai,gi1,gi4,lami,zdr,lamimin,  &
              lamimax,mvdi,ltk,lbdi,adagr,zeta,dbzi,fdbzi,llmi,llmi2,  &
              lai2,lai22,lai23,lroi,lroi2,lroi3,lmdi,lmdi2,lmdi3,zeta2,&
              zeta3,zhrox,zhzet,zdrox,zdzet,fzdri,fkdpi,kdrox,kdzet,   &
              i3m,tc,ttr,inhgr
      real(kind=8), parameter :: thrd = 1./3.,         zxmin = 1.e-21
      real(kind=8), parameter :: afamin = 0.,          afamax = 30.
      real(kind=8), parameter :: afa1 = -0.35897435,   afa2 = 0.010041278
      real(kind=8), parameter :: afa3 = 1.9349749E-3,  afa4 = -0.067453607
      real(kind=8), parameter :: afa5 = -5.5693923E-3, afa6 = 0.06376154
      real(kind=8), parameter :: afb1 = -6.5100404,    afb2 = 0.039597245
      real(kind=8), parameter :: afb3 = 2.3169628E-3,  afb4 = -0.45856646
      real(kind=8), parameter :: afb5 = -0.04537631,   afb6 = 0.25713479
      real(kind=8), dimension(0:120) :: itble
      data itble /1.000000,0.979490,0.959401,0.939723,0.920450,0.899498,&
                  0.879023,0.857038,0.833681,0.810961,0.783430,0.755092,&
                  0.703072,0.537032,0.467735,0.524807,0.630957,0.812831,&
                  1.096478,1.479108,1.905461,2.089296,2.290868,2.398833,&
                  2.454709,2.426610,2.371374,2.290868,2.137962,1.995262,&
                  1.862087,1.737801,1.621810,1.513561,1.396368,1.288250,&
                  1.188502,1.096478,1.000000,0.922571,0.851138,0.785236,&
                  0.724436,0.668344,0.616595,0.575440,0.537032,0.501187,&
                  0.467735,0.436516,0.407380,0.380189,0.354813,0.331131,&
                  0.316228,0.301995,0.291743,0.285102,0.281838,0.278612,&
                  0.275423,0.278612,0.281838,0.285102,0.291743,0.298538,&
                  0.309030,0.319890,0.331131,0.346737,0.367282,0.393550,&
                  0.426580,0.457088,0.489779,0.524807,0.562341,0.609537,&
                  0.660693,0.716143,0.785236,0.860994,0.954993,1.047129,&
                  1.148154,1.258925,1.380384,1.496236,1.603245,1.698244,&
                  1.778279,1.840772,1.883649,1.905461,1.905461,1.883649,&
                  1.862087,1.840772,1.798871,1.737801,1.698244,1.640590,&
                  1.584893,1.548817,1.513561,1.475707,1.452112,1.428894,&
                  1.412538,1.393157,1.377209,1.361445,1.348963,1.336596,&
                  1.327394,1.318257,1.309182,1.303167,1.294196,1.288250,&
                  1.279381/

      zh = zxmin
      zv = zxmin
      kdp = 0.
      qc = qc0/1.e3      ! g/kg -> kg/kg
      qi = qi0/1.e3      ! g/kg -> kg/kg
      ni = ni0/rho       ! 1/m^3 -> 1/kg
      if (iband.eq.1) then
         dbzwl = 0.11    ! for S band
      else
         dbzwl = 0.053   ! for C band
      endif
      if (qi.ge.1.e-14.and.ni.ge.1.e-2) then
         rhoi = 900.
         mvri = MIN(dimax/2.,MAX(dimin/2.,(1.90986*qi/ni/rhoi)**thrd))
         if (qc.ge.1.e-14.or.ni.lt.0.5) then
            afai = MIN(afamax,MAX(afamin,(afa1+afa2*LOG(ni)+afa3*      &
                   LOG(ni)**2.+afa4*LOG(mvri))/(1.+afa5*LOG(ni)+afa6*  &
                   LOG(mvri))))
         else
            afai = MIN(afamax,MAX(afamin,(afb1+afb2*LOG(ni)+afb3*      &
                   LOG(ni)**2.+afb4*LOG(mvri))/(1.+afb5*LOG(ni)+afb6*  &
                   LOG(mvri))))
         endif
         gi1  = GAMLN(afai+1.)
         gi4  = GAMLN(afai+4.)
         i3m  = 1.90986*qi/rhoi
         lami = (EXP(gi4-gi1)*ni/i3m)**thrd
         lamimin = (EXP(gi4-gi1))**thrd/dimax
         lamimax = (EXP(gi4-gi1))**thrd/dimin
         IF (lami.LT.lamimin) THEN
            lami = lamimin
            ni = i3m*EXP(gi1-gi4+3.*LOG(lami))
         ELSEIF (lami.GT.lamimax) THEN
            lami = lamimax
            ni = i3m*EXP(gi1-gi4+3.*LOG(lami))
         ENDIF
         mvdi  = (EXP(gi4-gi1))**thrd/lami
         IF (mvdi.GE.6.e-6) THEN
            lbdi  = LOG(mvdi*1.E6)
            ttr   = MIN(1.,MAX(0.01,9.207E-3*lbdi**3.-0.150212*lbdi**2.+&
                    0.565632*lbdi+0.401116))
            tc    = tk-273.15
            hid   = MAX(MIN(NINT(ABS(tc)/0.25),120),0)
            inhgr = itble(hid)
            adagr = inhgr**ttr
         ELSE
            adagr = 1.
         ENDIF
         zeta = (adagr-1.)/(adagr+2.)
      if (qi.ge.1.e-6) then
         llmi  = LOG(lami)
         lai2  = LOG(afai+2.)
         lroi  = LOG(rhoi)
         lmdi  = LOG(1.E6*mvdi)
         llmi2 = llmi*llmi
         lai22 = lai2*lai2
         lai23 = lai22*lai2
         lmdi2 = lmdi*lmdi
         lmdi3 = lmdi2*lmdi
         lroi2 = lroi*lroi
         lroi3 = lroi2*lroi
         zeta2 = zeta*zeta
         zeta3 = zeta2*zeta
         IF ((adagr-1.).GE.1.e-2) THEN
            fdbzi = 389.04082+18.161704*lai2+1.3401764*lai22-          &
                    0.07360894*lai23-26.137103*llmi
            zhrox = 1.9382E-6*EXP(1.9260732373*lroi)
            zhzet = 1.
            zh    = ni*10.**((fdbzi-200.)/10.-18.)*zhrox*zhzet
            dbzi  = 10.*LOG10(1.E18*zh)
            zdrox = 0.0585053*lroi3-0.80223023*lroi2+3.75756297*lroi-  &
                    5.86526003
            zdzet = MAX(0.,-0.01644484+0.010454643*lmdi+9.4367973*zeta-&
                    0.0026137203*lmdi2+5.7453232*zeta2+0.3234453*lmdi* &
                    zeta+2.0622912E-4*lmdi3-1.3945604*zeta3-2.2879493* &
                    lmdi*zeta2-0.013856582*lmdi2*zeta)
            fzdri = EXP((15.181145+0.34794935*lai2+0.11677701*lai22-   &
                    6.3942628E-3*lai23-1.0984758*llmi)/(1.+0.022930408*&
                    lai2+7.5639466E-3*lai22-4.146688E-4*lai23-         &
                    0.071648235*llmi))/1.E6*zdzet*zdrox
            zv    = 10.**(dbzi/10.-18.+fzdri/10.)
            fkdpi = EXP((54.728473+7.4097498*lai2-7.5413134*llmi+      &
                    0.25580463*lai22+0.25257362*llmi2-0.50222586*lai2* &
                    llmi)/(1.+0.081702748*lai2-0.085895343*llmi+       &
                    7.6712618E-4*lai22+4.5412963E-4*llmi2-9.3521894E-4*&
                    lai2*llmi))/2.E18
            kdrox = 1.75E-6*EXP(1.93990711*lroi)
            kdzet = MAX(0.,-0.0067191022+0.030489666*lmdi+7.6470878*   &
                    zeta-0.01097138*lmdi2+9.8497711*zeta2+0.83287327*  &
                    lmdi*zeta+9.1590015E-4*lmdi3-7.2917335*zeta3-      &
                    2.2235781*lmdi*zeta2-0.053116593*lmdi2*zeta)
            kdp   = ni*fkdpi*kdrox*kdzet*0.11/dbzwl
            IF (zv.LT.zh) zv = 1.0001*zh
            zdr   = 10.*LOG10(zv/zh)
         ELSEIF ((1.-adagr).GE.1.e-2) THEN
            fdbzi = 390.47265+18.157805*lai2+1.3558212*lai22-          &
                    0.074452587*lai23-26.228751*llmi
            zhrox = 1.4172E-6*EXP(1.9687416987*lroi)
            zhzet = 1.1531946-0.10474946*lmdi+1.3268439*zeta+          &
                    0.013483871*lmdi2+3.0239757*zeta2-0.53774565*lmdi* &
                    zeta-7.6839476E-4*lmdi3+4.3808336*zeta3+           &
                    0.080046674*lmdi*zeta2+0.015810946*lmdi2*zeta
            zh    = ni*10.**((fdbzi-200.)/10.-18.)*zhrox*zhzet
            dbzi  = 10.*LOG10(1.E18*zh)
            zdrox = 0.05683556*lroi3-0.77585835*lroi2+3.62207807*lroi- &
                    5.63826259
            zdzet = MAX(0.,0.12567123-0.077409338*lmdi-8.3245922*zeta+ &
                    0.013440196*lmdi2+8.7068519*zeta2-0.53878869*lmdi* &
                    zeta-7.0609359E-4*lmdi3+6.371616*zeta3-1.6563927*  &
                    lmdi*zeta2+0.032682807*lmdi2*zeta)
            fzdri = EXP((15.587817+0.23578544*lai2+0.045497555*lai22-  &
                    1.0862071*llmi+3.4688369E-3*llmi2)/(1.+0.015556043*&
                    lai2+2.1779895E-3*lai22+7.2733931E-5*lai23-        &
                    0.064401818*llmi))/1.E6*zdzet*zdrox
            zv    = 10.**(dbzi/10.-18.-fzdri/10.)
            fkdpi = EXP((55.096057+7.4047265*lai2-7.6277374*llmi+      &
                    0.25324833*lai22+0.25735457*llmi2-0.50528828*lai2* &
                    llmi)/(1.+0.078240442*lai2-0.083742581*llmi+       &
                    5.4080817E-4*lai22+2.9610344E-4*llmi2-5.9621769E-4*&
                    lai2*llmi))/2.E18
            kdrox = 1.89E-6*EXP(1.92963337*lroi)
            kdzet = MAX(0.,0.044182106-0.020578987*lmdi-8.8659688*zeta+&
                    0.0022174675*lmdi2+5.2415098*zeta2-0.40918383*lmdi*&
                    zeta-3.696556E-5*lmdi3+3.5504491*zeta3-1.4353113*  &
                    lmdi*zeta2+0.022796706*lmdi2*zeta)
            kdp   = ni*fkdpi*kdrox*kdzet*0.11/dbzwl
            IF (zv.GT.zh) zv = 0.9999*zh
            zdr   = 10.*LOG10(zh/zv)
         ELSEIF (ABS(adagr-1.).LT.1.e-2) THEN
            fdbzi = 390.47265+18.157805*lai2+1.3558212*lai22-          &
                    0.074452587*lai23-26.228751*llmi
            zhrox = 1.9382E-6*EXP(1.9260732373*lroi)
            zhzet = 1.
            zh    = ni*10.**((fdbzi-200.)/10.-18.)*zhrox*zhzet
            dbzi  = 10.*LOG10(1.E18*zh)
            zv    = 10.**(dbzi/10.-18.)
            kdp   = 0.
            zdr   = 0.
         ENDIF
      endif
      endif

  end subroutine dualpol_op_ice_tcwa2
!---------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
      subroutine dualpol_op_snow_tcwa2(iband,tk,rho,qc0,qr0,qs0,ns0,   &
                                       smlf,zh,zv,kdp)
!-------------------------------------------------------------------------------------
      implicit none
      integer, intent(in) :: iband      ! radar scan band (wavelength)
      real(kind=8), intent(in)  :: tk   ! air temperature, K
      real(kind=8), intent(in)  :: rho  ! air density, kg/m^3
      real(kind=8), intent(in)  :: qc0  ! mass mixing ratio, g/kg
      real(kind=8), intent(in)  :: qr0  ! mass mixing ratio, g/kg
      real(kind=8), intent(in)  :: qs0  ! mass mixing ratio, g/kg
      real(kind=8), intent(in)  :: ns0  ! number concentration, 1/m^3
      real(kind=8), intent(in)  :: smlf ! melted fraction
      real(kind=8), intent(out) :: zh   ! horizontal reflectivity, mm^6/m^3
      real(kind=8), intent(out) :: zv   ! vertical reflectivity, mm^6/m^3
      real(kind=8), intent(out) :: kdp  ! specific differential phase, degree/km^-1
      real, parameter :: dsmin = 2.E-5, dsmax = 1.E-2
      real :: qc,qs,ns,dbzwl,ltk,ltk2,lqs,lqs2,rhos,s3m,afas,gs1,gs4,  &
              lams,lamsmin,lamsmax,mvds,saspr,smlv,dsmm,sca0,scala,    &
              scalb,scafa,scafb,tc,ora0,orar,ora1,ora2,ora3,scaw1,rads,&
              zdr,llms,epsin,epss,alphe,lamdz,nerel,scawa,scawb,scaw0, &
              scakw,las2,lros,lasp,fdbzs,zhrox,dbzs,zdrox,zdasp,fzdrs, &
              fkdps,kdrox,kdasp,qr,mvrs,lbds
      real(kind=8), parameter :: thrd = 1./3.,         zxmin = 1.e-21
      real(kind=8), parameter :: afamin = 0.,          afamax = 30.
      real(kind=8), parameter :: afa1 = -0.35897435,   afa2 = 0.010041278
      real(kind=8), parameter :: afa3 = 1.9349749E-3,  afa4 = -0.067453607
      real(kind=8), parameter :: afa5 = -5.5693923E-3, afa6 = 0.06376154
      real(kind=8), parameter :: afb1 = -6.5100404,    afb2 = 0.039597245
      real(kind=8), parameter :: afb3 = 2.3169628E-3,  afb4 = -0.45856646
      real(kind=8), parameter :: afb5 = -0.04537631,   afb6 = 0.25713479

      zh = zxmin
      zv = zxmin
      kdp = 0.
      qs = qs0/1.e3      ! g/kg -> kg/kg
      qc = qc0/1.e3      ! g/kg -> kg/kg
      qr = qr0/1.e3
      ns = ns0/rho       ! 1/m^3 -> 1/kg
      if (iband.eq.1) then
         dbzwl = 0.11    ! for S band
      else
         dbzwl = 0.053   ! for C band
      endif
      if (qs.ge.1.e-14.and.ns.ge.1.e-2) then
         ltk  = LOG(tk)
         ltk2 = ltk*ltk
         lqs  = -1.*LOG(qs)
         lqs2 = lqs*lqs
         IF (tk.LT.273.15) THEN
            rhos = MIN(900.,EXP(75.401034-25.730216*ltk+0.81467519*lqs+&
                   2.339118*ltk2+1.9976667E-3*lqs2-0.14086434*ltk*lqs))
         ELSE
            rhos = MAX(100.,EXP(-64808.666+23113.508*ltk-36.46632*lqs- &
                   2060.6024*ltk2-0.005729458*lqs2+6.5057411*ltk*lqs))
         ENDIF
         rhos = MIN(MAX(rhos,100.),910.)
!         smlf = MIN(1.-rhos/1000.,MAX(0.,(MIN(qr,qs)/MAX(qr,qs))**0.3))
         mvrs = MIN(dsmax/2.,MAX(dsmin/2.,(1.90986*qs/ns/rhos)**thrd))
         if (qc.ge.1.e-14.or.ns.lt.0.5) then
            afas = MIN(afamax,MAX(afamin,(afa1+afa2*LOG(ns)+afa3*      &
                   LOG(ns)**2.+afa4*LOG(mvrs))/(1.+afa5*LOG(ns)+afa6*  &
                   LOG(mvrs))))
         else
            afas = MIN(afamax,MAX(afamin,(afb1+afb2*LOG(ns)+afb3*      &
                   LOG(ns)**2.+afb4*LOG(mvrs))/(1.+afb5*LOG(ns)+afb6*  &
                   LOG(mvrs))))
         endif
         if (tk.gt.273.15) afas = MIN(afamax,0.01*EXP(-1.*LOG(mvrs)))
         gs1  = GAMLN(afas+1.)
         gs4  = GAMLN(afas+4.)
         s3m  = 1.90986*qs/rhos
         lams = (EXP(gs4-gs1)*ns/s3m)**thrd
         lamsmin = (EXP(gs4-gs1))**thrd/dsmax
         lamsmax = (EXP(gs4-gs1))**thrd/dsmin
         IF (lams.LT.lamsmin) THEN
            lams = lamsmin
            ns = s3m*EXP(gs1-gs4+3.*LOG(lams))
         ELSEIF (lams.GT.lamsmax) THEN
            lams = lamsmax
            ns = s3m*EXP(gs1-gs4+3.*LOG(lams))
         ENDIF
         mvds  = (EXP(gs4-gs1))**thrd/lams
      if (qs.ge.1.e-6) then
         llms  = LOG(lams)
         lbds  = LOG(mvds*1.E6)
         lros  = LOG(rhos)
         saspr = MAX(0.2,MIN(1.,1.0207483-0.40684481*lros+0.043489251* &
                 lros**2.+0.020462134*lbds+0.014926001*lbds**2.-       &
                 1.1252149E-3*lbds**3.))
         IF (tk.LT.238.16) saspr = 0.7
         IF (smlf.GE.0.01.and.tk.ge.273.16) THEN
            tc    = tk-273.15
            epsin = 5.27137+2.16474E-2*tc-1.31198E-3*tc**2.
            epss  = 78.54*(1.-4.579E-3*(tc-25.)+1.19E-5*(tc-25.)**2.-  &
                    2.8E-8*(tc-25.)**3.)
            alphe = -16.8129/tk+6.09265E-2
            lamdz = 3.3836E-4*EXP(2513.98/tk)*1.E-2
            nerel = 1.+2.*(lamdz/dbzwl)**(1.-alphe)*SIN(1.570796327*   &
                    alphe)+(lamdz/dbzwl)**(2.-2.*alphe)
            scawa = epsin+((epss-epsin)*((lamdz/dbzwl)**(1.-alphe)*    &
                    SIN(1.570796327*alphe)+1.))/nerel
            scawb = ((epss-epsin)*((lamdz/dbzwl)**(1.-alphe)*          &
                    COS(1.570796327*alphe)))/nerel+1.25664*dbzwl/1.88496
            scaw0 = scawa**2.+4.*scawa+4.+scawb**2.
            scakw = ((scawa**2.+scawa-2.+scawb**2.)/scaw0)**2.+(3.*    &
                    scawb/scaw0)**2.
            s3m   = 1.90986*qs/rhos
            mvrs  = MIN(dsmax/2.,MAX(dsmin/2.,(1.90986*qs/ns/          &
                    rhos)**thrd))
            afas  = MIN(afamax,MAX(afamin,EXP(-18.420681-0.5*LOG(ns)-  &
                    2.94*LOG(mvrs))))
            gs1   = GAMLN(afas+1.)
            gs4   = GAMLN(afas+4.)
            lams  = (EXP(gs4-gs1)*ns/s3m)**thrd
            lamsmin = (EXP(gs4-gs1))**thrd/dsmax
            lamsmax = (EXP(gs4-gs1))**thrd/dsmin
            IF (lams.LT.lamsmin) THEN
               lams = lamsmin
               ns = s3m*EXP(gs1-gs4+3.*LOG(lams))
            ELSEIF (lams.GT.lamsmax) THEN
               lams = lamsmax
               ns = s3m*EXP(gs1-gs4+3.*LOG(lams))
            ENDIF
            mvds  = (EXP(gs4-gs1))**thrd/lams
            smlv  = (rhos/(1.e3/smlf-1.e3+rhos))**0.3
            dsmm  = mvds*1.E3
            saspr = smlf*MAX(0.4,0.9951+2.51E-2*dsmm-3.644E-2*dsmm**2.+&
                    5.303E-3*dsmm**3.-2.492E-4*dsmm**4.)+(1.-smlf)*saspr
            rhos  = MIN(910.,MAX(50.,smlf*1.e3+(1.-smlf)*rhos))
            ora0  = 10.*smlf+(1.-smlf)*40.
            orar  = EXP(-2.*(ora0*0.017453293)**2.)
            ora1  = (3.+4.*orar+orar**4.)/8.
            ora2  = (3.-4.*orar+orar**4.)/8.
            ora3  = (1.-orar**4.)/8.
            sca0  = MAX(1.E-6,SQRT(saspr**(-2.)-1.))
            scala = (1.+sca0**2.)/sca0**2.*(1.-ATAN(sca0)/sca0)
            scalb = (1.-scala)/2.
            scaw1 = -1.0892332E-2+2.0820143E-2/smlv+4.8553456E-6*rhos- &
                    7.6529057E-5/smlv**2.-1.1783334E-9*rhos**2.-       &
                    2.9599665E-6*rhos/smlv
            scafa = 1./(scala+scaw1)
            scafb = 1./(scalb+scaw1)
            rads  = ns*EXP(GAMLN(afas+7.)-gs1-6.*llms)
            zh    = MAX(1.E-21,rads*(ora1*scafb**2.+2.*ora3*scafa*     &
                    scafb+ora2*scafa**2.)/(9.*scakw))
            zh    = zh*(MAX(0.,-0.1*afas+0.5)*smlf+1.)
            zv    = MAX(1.E-21,rads*(ora1*scafa**2.+2.*ora3*scafa*     &
                    scafb+ora2*scafb**2.)/(9.*scakw))
            kdp   = 94247.77961/dbzwl*ns*EXP(gs4-gs1-3.*llms)*         &
                    ABS(scafa-scafb)*orar
            IF (zv.GT.zh) zv = 0.9999*zh
            zdr = 10.*LOG10(zh/zv)
         ELSE
            las2  = LOG(afas+2.)
            lros  = LOG(rhos)
            lasp  = LOG(200.*MAX(0.01,saspr))
            fdbzs = 428.62235+18.226687*las2+1.299724*las2**2.-        &
                    0.070558679*las2**3.-26.070123*llms
            zhrox = 1.544E-5*EXP(1.95714566*lros)
            zh    = ns*10.**((fdbzs-250.)/10.-18.)*zhrox
            dbzs  = 10.*LOG10(1.E18*zh)
            zdrox = 0.14542522*lros**3.-1.97576772*lros**2.+           &
                    9.33469066*lros-14.95051206
            zdasp = MAX(0.,-8.85535E-3*lasp**3.-0.10877652*lasp**2.+   &
                    0.24567081*lasp+3.06007047)
            fzdrs = 0.384*zdrox*zdasp/0.37728
            zv    = 10.**(dbzs/10.-18.-fzdrs/10.)
            fkdps = EXP(44.30884+2.9903578*las2-0.25943237/las2-       &
                    3.0036608*llms)/1.E15
            kdrox = 1.974E-5*EXP(1.90820312*lros)
            kdasp = MAX(0.,-0.00958782*lasp**3.-0.10043861*lasp**2.+   &
                    0.22095032*lasp+3.06677159)
            kdp   = ns*fkdps*kdrox*kdasp*0.11/dbzwl
            IF (zv.GT.zh) zv = 0.9999*zh
            zdr   = 10.*LOG10(zh/zv)
         ENDIF
      endif
      endif

  end subroutine dualpol_op_snow_tcwa2
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
      subroutine dualpol_op_graup_tcwa2(iband,tk,rho,qc0,qr0,qg0,ng0,  &
                                        gmlf,zh,zv,kdp)
!---------------------------------------------------------------------------------
      implicit none
      integer, intent(in) :: iband      ! radar scan band (wavelength)
      real(kind=8), intent(in)  :: tk   ! air temperature, K
      real(kind=8), intent(in)  :: rho  ! air density, kg/m^3
      real(kind=8), intent(in)  :: qc0  ! mass mixing ratio, g/kg
      real(kind=8), intent(in)  :: qr0  ! mass mixing ratio, g/kg
      real(kind=8), intent(in)  :: qg0  ! mass mixing ratio, g/kg
      real(kind=8), intent(in)  :: ng0  ! number concentration, 1/m^3
      real(kind=8), intent(in)  :: gmlf ! melted fraction
      real(kind=8), intent(out) :: zh   ! horizontal reflectivity, mm^6/m^3
      real(kind=8), intent(out) :: zv   ! vertical reflectivity, mm^6/m^3
      real(kind=8), intent(out) :: kdp  ! specific differential phase, degree/km^-1
      real, parameter :: dgmin = 5.E-5, dgmax = 2.E-2
      integer :: ihail
      real :: qg,qc,qr,ng,dbzwl,ltk,ltk2,lqg,rhog,tc,mvg0,dsll,g3m,bdg,&
              lbdg,afag,gg1,gg4,lamg,lamgmin,lamgmax,mvdg,gaspr,llmg,  &
              llmg2,gmlv,dgmm,lag5,lag52,zdr,epsin,epss,alphe,lamdz,   &
              nerel,scawa,scawb,scaw0,scakw,ora0,orar,ora1,ora2,ora3,  &
              sca0,scala,scalb,scaw1,scafa,scafb,radg,lag2,lrog,lasp,  &
              zhrox,fdbzg,dbzg,zdrox,zdasp,fzdrg,fkdpg,kdrox,kdasp,mvrg
      real(kind=8), parameter :: thrd = 1./3.,         zxmin = 1.e-21
      real(kind=8), parameter :: afamin = 0.,          afamax = 30.
      real(kind=8), parameter :: afa1 = -0.35897435,   afa2 = 0.010041278
      real(kind=8), parameter :: afa3 = 1.9349749E-3,  afa4 = -0.067453607
      real(kind=8), parameter :: afa5 = -5.5693923E-3, afa6 = 0.06376154
      real(kind=8), parameter :: afb1 = -6.5100404,    afb2 = 0.039597245
      real(kind=8), parameter :: afb3 = 2.3169628E-3,  afb4 = -0.45856646
      real(kind=8), parameter :: afb5 = -0.04537631,   afb6 = 0.25713479

      zh = zxmin
      zv = zxmin
      kdp = 0.
      qg = qg0/1.e3      ! g/kg -> kg/kg
      qc = qc0/1.e3      ! g/kg -> kg/kg
      qr = qr0/1.e3      ! g/kg -> kg/kg
      ng = ng0/rho       ! 1/m^3 -> 1/kg
      if (iband.eq.1) then
         dbzwl = 0.11    ! for S band
      else
         dbzwl = 0.053   ! for C band
      endif

      if (qg.ge.1.e-14.and.ng.ge.1.e-2) then
         ltk   = LOG(tk)
         ltk2  = ltk*ltk
         ihail = 0
         rhog  = 400.
         tc    = tk-273.15
         mvg0  = EXP(8.8283+0.25*LOG(qg))*1.E-6
         dsll  = MAX(1.E-3,MIN(1.,2.*1.E-2*(EXP(MIN(20.,-tc/MAX(0.1,(  &
                 1.1E4*(qc+qr)-1.3E3*qg+1.))))-1.)))
         IF (mvg0.GT.dsll) THEN
            ihail = 1
            rhog = 900.
         ELSE
            ihail = 0
            lqg  = -1.*LOG(qg)
            rhog = MIN(900.,MAX(400.,EXP(-9.57E-5*lqg**3.+3.077E-3*    &
                   lqg**2.-6.800923E-2*lqg+6.8175231)))
         ENDIF
!         gmlf = MIN(1.-rhog/1000.,MAX(0.,(MIN(qr,qg)/MAX(qr,qg))**0.3))
         mvrg = MIN(dgmax/2.,MAX(dgmin/2.,(1.90986*qg/ng/rhog)**thrd))
         if (qc.ge.1.e-14.or.ng.lt.1.) then
            afag = MIN(afamax,MAX(afamin,(afa1+afa2*LOG(ng)+afa3*      &
                   LOG(ng)**2.+afa4*LOG(mvrg))/(1.+afa5*LOG(ng)+afa6*  &
                   LOG(mvrg))))
         else
            afag = MIN(afamax,MAX(afamin,(afb1+afb2*LOG(ng)+afb3*      &
                   LOG(ng)**2.+afb4*LOG(mvrg))/(1.+afb5*LOG(ng)+afb6*  &
                   LOG(mvrg))))
         endif
         if (tk.gt.273.15) afag = MIN(afamax,0.01*EXP(-1.*LOG(mvrg)))
         gg1  = GAMLN(afag+1.)
         gg4  = GAMLN(afag+4.)
         g3m  = 1.90986*qg/rhog
         lamg = (EXP(gg4-gg1)*ng/g3m)**thrd
         lamgmin = (EXP(gg4-gg1))**thrd/dgmax
         lamgmax = (EXP(gg4-gg1))**thrd/dgmin
         IF (lamg.LT.lamgmin) THEN
            lamg = lamgmin
            ng = g3m*EXP(gg1-gg4+3.*LOG(lamg))
         ELSEIF (lamg.GT.lamgmax) THEN
            lamg = lamgmax
            ng = g3m*EXP(gg1-gg4+3.*LOG(lamg))
         ENDIF
         mvdg  = (EXP(gg4-gg1))**thrd/lamg
      if (qg.ge.1.e-6) then
         llmg  = LOG(lamg)
         llmg2 = llmg*llmg
         IF (ihail.EQ.0) THEN
            lbdg  = LOG(mvdg*1.E6)
            lrog  = LOG(rhog)
            gaspr = MAX(0.2,MIN(1.,1.0207483-0.40684481*lrog+          &
                    0.043489251*lrog**2.+0.020462134*lbdg+0.014926001* &
                    lbdg**2.-1.1252149E-3*lbdg**3.))
         ELSEIF (ihail.EQ.1) THEN
            lag5  = LOG(afag+5.)
            lag52 = lag5*lag5
            gaspr = MIN(1.,((5.6387098+1.5626893*lag5-1.5481826*llmg+  &
                    0.12481435*lag52+0.11244742*llmg2-0.2336531*lag5*  &
                    llmg)/(1.+0.28272803*lag5-0.2780308*llmg+          &
                    2.1750225E-2*lag52+1.9960971E-2*llmg2-4.1317594E-2*&
                    lag5*llmg))/10.)
         ENDIF
         IF (gmlf.GE.0.01.and.tk.ge.273.16) THEN
            epsin = 5.27137+2.16474E-2*tc-1.31198E-3*tc**2.
            epss  = 78.54*(1.-4.579E-3*(tc-25.)+1.19E-5*(tc-25.)**2.-  &
                    2.8E-8*(tc-25.)**3.)
            alphe = -16.8129/tk+6.09265E-2
            lamdz = 3.3836E-4*EXP(2513.98/tk)*1.E-2
            nerel = 1.+2.*(lamdz/dbzwl)**(1.-alphe)*SIN(1.570796327*   &
                    alphe)+(lamdz/dbzwl)**(2.-2.*alphe)
            scawa = epsin+((epss-epsin)*((lamdz/dbzwl)**(1.-alphe)*    &
                    SIN(1.570796327*alphe)+1.))/nerel
            scawb = ((epss-epsin)*((lamdz/dbzwl)**(1.-alphe)*          &
                    COS(1.570796327*alphe)))/nerel+1.25664*dbzwl/1.88496
            scaw0 = scawa**2.+4.*scawa+4.+scawb**2.
            scakw = ((scawa**2.+scawa-2.+scawb**2.)/scaw0)**2.+(3.*    &
                    scawb/scaw0)**2.
            g3m   = 1.90986*qg/rhog
            mvrg  = MIN(dgmax/2.,MAX(dgmin/2.,(1.90986*qg/ng/          &
                    rhog)**thrd))
            afag  = MIN(afamax,MAX(afamin,EXP(-18.420681-0.5*LOG(ng)-  &
                    2.94*LOG(mvrg))))
            gg1   = GAMLN(afag+1.)
            gg4   = GAMLN(afag+4.)
            lamg  = (EXP(gg4-gg1)*ng/g3m)**thrd
            lamgmin = (EXP(gg4-gg1))**thrd/dgmax
            lamgmax = (EXP(gg4-gg1))**thrd/dgmin
            IF (lamg.LT.lamgmin) THEN
               lamg = lamgmin
               ng = g3m*EXP(gg1-gg4+3.*LOG(lamg))
            ELSEIF (lamg.GT.lamgmax) THEN
               lamg = lamgmax
               ng = g3m*EXP(gg1-gg4+3.*LOG(lamg))
            ENDIF
            mvdg  = (EXP(gg4-gg1))**thrd/lamg
            gmlv  = (rhog/(1.e3/gmlf-1.e3+rhog))**0.3
            dgmm  = mvdg*1.E3
            gaspr = gmlf*MAX(0.4,0.9951+2.51E-2*dgmm-3.644E-2*dgmm**2.+&
                    5.303E-3*dgmm**3.-2.492E-4*dgmm**4.)+(1.-gmlf)*gaspr
            rhog  = MIN(900.,MAX(100.,gmlf*1.e3+(1.-gmlf)*rhog))
            ora0  = 10.*gmlf+(1.-gmlf)*40.
            orar  = EXP(-2.*(ora0*0.017453293)**2.)
            ora1  = (3.+4.*orar+orar**4.)/8.
            ora2  = (3.-4.*orar+orar**4.)/8.
            ora3  = (1.-orar**4.)/8.
            sca0  = MAX(1.E-6,SQRT(gaspr**(-2.)-1.))
            scala = (1.+sca0**2.)/sca0**2.*(1.-ATAN(sca0)/sca0)
            scalb = (1.-scala)/2.
            scaw1 = -1.0892332E-2+2.0820143E-2/gmlv+4.8553456E-6*rhog- &
                    7.6529057E-5/gmlv**2.-1.1783334E-9*rhog**2.-       &
                    2.9599665E-6*rhog/gmlv
            scafa = 1./(scala+scaw1)
            scafb = 1./(scalb+scaw1)
            radg  = ng*EXP(GAMLN(afag+7.)-gg1-6.*llmg)
            zh    = MAX(1.E-21,radg*(ora1*scafb**2.+2.*ora3*scafa*     &
                    scafb+ora2*scafa**2.)/(9.*scakw))
            zh    = zh*(MAX(0.,-0.1*afag+0.5)*gmlf+1.)
            zv    = MAX(1.E-21,radg*(ora1*scafa**2.+2.*ora3*scafa*     &
                    scafb+ora2*scafb**2.)/(9.*scakw))
            kdp   = 94247.77961/dbzwl*ng*EXP(gg4-gg1-3.*llmg)*         &
                    ABS(scafa-scafb)*orar
            IF (zv.GT.zh) zv = 0.9999*zh
            zdr   = 10.*LOG10(zh/zv)
         ELSE
            lag2  = LOG(afag+2.)
            lrog  = LOG(rhog)
            lasp  = LOG(200.*MAX(0.01,gaspr))
            zhrox = 1.544E-5*EXP(1.95714566*lrog)
            fdbzg = 428.62235+18.226687*lag2+1.299724*lag2**2.-        &
                    0.070558679*lag2**3.-26.070123*llmg
            zh    = ng*10.**((fdbzg-250.)/10.-18.)*zhrox
            dbzg  = 10.*LOG10(1.E18*zh)
            zdrox = 0.14542522*lrog**3.-1.97576772*lrog**2.+9.33469066*&
                    lrog-14.95051206
            zdasp = MAX(0.,-8.85535E-3*lasp**3.-0.10877652*lasp**2.+   &
                    0.24567081*lasp+3.06007047)
            fzdrg = 0.384*zdrox*zdasp
            zv    = 10.**(dbzg/10.-18.-fzdrg/10.)
            fkdpg = EXP(44.30884+2.9903578*lag2-0.25943237/lag2-       &
                    3.0036608*llmg)/1.E15
            kdrox = 1.974E-5*EXP(1.90820312*lrog)
            kdasp = MAX(0.,-0.00958782*lasp**3.-0.10043861*lasp**2.+   &
                    0.22095032*lasp+3.06677159)
            kdp   = ng*fkdpg*kdrox*kdasp*0.11/dbzwl
            IF (zv.GT.zh) zv = 0.9999*zh
            zdr = 10.*LOG10(zh/zv)
         ENDIF
      endif
      endif

  end subroutine dualpol_op_graup_tcwa2
!---------------------------------------------------------------------------------------
!======================================================================
  REAL FUNCTION GAMLN(XX)
!======================================================================
!     --- RETURNS THE VALUE LN(GAMMA(XX)) FOR XX > 0.
    IMPLICIT NONE
    REAL, INTENT(IN):: XX
    DOUBLE PRECISION, PARAMETER:: STP = 2.5066282746310005D0
    DOUBLE PRECISION, DIMENSION(6), PARAMETER:: &
             COF = (/76.18009172947146D0, -86.50532032941677D0, &
                     24.01409824083091D0, -1.231739572450155D0, &
                    .1208650973866179D-2, -.5395239384953D-5/)
    DOUBLE PRECISION:: SER,TMP,X,Y
    INTEGER:: J

    X = XX
    Y = X
    TMP = X+5.5D0
    TMP = (X+0.5D0)*LOG(TMP)-TMP
    SER = 1.000000000190015D0
    DO 11 J = 1,6
      Y = Y+1.D0
      SER = SER+COF(J)/Y
11  CONTINUE
    GAMLN = TMP+LOG(STP*SER/X)

  END FUNCTION GAMLN
!  (C) Copr. 1986-92 Numerical Recipes Software 2.02
!======================================================================

end module tcwa2_forward_mod

