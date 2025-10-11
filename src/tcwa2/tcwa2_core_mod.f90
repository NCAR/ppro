module tcwa2_core_mod
!--------------------------
! TCWA2 polarimetric radar operator
! 
! This module contains the TCWA2-specific dual-polarization radar forward
! operator subroutines. These use analytical formulations based on gamma
! distribution parameters rather than lookup tables.
!
! Reference: Tsai et al. TCWA2 microphysics scheme
! 
! Author: Tzu-Chin Tsai (original implementation)
! Integrated into PPRO library: 2025
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
  subroutine dualpol_op_rain_tcwa2(iband,rho,qr0,nr0,zh,zv,kdp)
!--------------------------------------------------------------------------------------
    implicit none
      integer, intent(in) :: iband       ! radar scan band (wavelength)
      real(kind=8), intent(in)  :: rho   ! air density, kg/m^3
      real(kind=8), intent(in)  :: qr0   ! mass mixing ratio, g/kg
      real(kind=8), intent(in)  :: nr0   ! number concentration, 1/m^3
      real(kind=8), intent(out) :: zh    ! horizontal reflectivity, mm^6/m^3
      real(kind=8), intent(out) :: zv    ! vertical reflectivity, mm^6/m^3
      real(kind=8), intent(out) :: kdp   ! specific differential phase, degree/km^-1
      real, parameter :: thrd = 1./3., drmin = 1.E-4, drmax = 6.E-3
      real :: qr,nr,mvrr,afar,lamr,lamrmin,lamrmax,gr1,gr4,mvdr,llmr,  &
              llmr2,llmr3,lar2,lar22,lar23,fdbzr,dbzr,fzdrr,zdr,fkdpr, &
              dbzwl

      zh = 1.e-21
      zv = 1.e-21
      kdp = 0.
      qr = qr0/1.e3      ! g/kg -> kg/kg
      nr = nr0/rho       ! 1/m^3 -> 1/kg
      if (iband.eq.1) then
         dbzwl = 0.11    ! for S band
      else
         dbzwl = 0.053   ! for C band
      endif
      if (qr.ge.1.e-14.and.nr.ge.1.e-2) then
         mvrr = MIN(drmax/2.,MAX(drmin/2.,(qr/nr/4.18879E3)**thrd))
         afar = MIN(20.,MAX(1.E-4,(-0.69914266+2.3818105E-3*LOG(nr)+   &
                1.6103668E-3*LOG(nr)**2.-0.12277336*LOG(mvrr))/(1.-    &
                2.3054262E-3*LOG(nr)+0.092294568*LOG(mvrr))))
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
         IF (nr.ge.1.) THEN
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
         ELSE
         fzdrr = 0.
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
  subroutine dualpol_op_ice_tcwa2(iband,tk,rho,qi0,ni0,zh,zv,kdp)
!-------------------------------------------------------------------------------------
    implicit none
      integer, intent(in) :: iband      ! radar scan band (wavelength)
      real(kind=8), intent(in)  :: tk   ! air temperature, K
      real(kind=8), intent(in)  :: rho  ! air density, kg/m^3
      real(kind=8), intent(in)  :: qi0  ! mass mixing ratio, g/kg
      real(kind=8), intent(in)  :: ni0  ! number concentration, 1/m^3
      real(kind=8), intent(out) :: zh   ! horizontal reflectivity, mm^6/m^3
      real(kind=8), intent(out) :: zv   ! vertical reflectivity, mm^6/m^3
      real(kind=8), intent(out) :: kdp  ! specific differential phase, degree/km^-1
      integer :: hid
      real, parameter :: thrd = 1./3., dimin = 6.E-6, dimax = 5.E-3
      real :: qi,ni,dbzwl,rhoi,i3m,bdi,fdi,kdx,afai,gi1,gi4,lami,zdr,  &
              lamimin,lamimax,mvdi,ltk,lbdi,adagr,zeta,dbzi,fdbzi,llmi,&
              llmi2,lai2,lai22,lai23,lroi,lroi2,lroi3,lmdi,lmdi2,lmdi3,&
              zeta2,zeta3,zhrox,zhzet,zdrox,zdzet,fzdri,fkdpi,kdrox,   &
              kdzet

      zh = 1.e-21
      zv = 1.e-21
      kdp = 0.
      qi = qi0/1.e3      ! g/kg -> kg/kg
      ni = ni0/rho       ! 1/m^3 -> 1/kg
      if (iband.eq.1) then
         dbzwl = 0.11    ! for S band
      else
         dbzwl = 0.053   ! for C band
      endif
      if (qi.ge.1.e-14.and.ni.ge.1.e-2) then
         rhoi = 900.
         i3m  = 1.909859317*qi/rhoi
         bdi  = MIN(MAX((i3m/ni)**thrd*1.E3,dimin*1.E3),dimax*1.E3)
         fdi  = 7.4015986E-2+7.9866676E-1*bdi-9.4468892E-3*LOG(ni)+    &
                3.8235092E-1*bdi**2.+2.9811542E-4*LOG(ni)**2.+         &
                1.9052614E-2*bdi*LOG(ni)
         kdx  = MAX(0.556,MIN(0.999,(bdi/fdi)**3.))
         afai = (6.*kdx-3.+SQRT(8.*kdx+1.))/(2.-2.*kdx)
         afai = MIN(MAX(afai,0.),3.E4)
         gi1  = GAMLN(afai+1.)
         gi4  = GAMLN(afai+4.)
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
         IF (mvdi.GT.6.e-6) THEN
            ltk   = LOG(tk)
            lbdi  = LOG(mvdi*1.E6)
            adagr = MIN(1.6,MAX(0.2,831.26148-299.86635*ltk+23.257362/ &
                    lbdi+27.074238*ltk**2.+0.004494653/lbdi**2.-       &
                    4.1838773*ltk/lbdi))
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
            kdp   = ni*fkdpi*kdrox*kdzet*dbzwl/0.11
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
            kdp   = ni*fkdpi*kdrox*kdzet*dbzwl/0.11
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
      real, parameter :: thrd = 1./3., dsmin = 2.E-5, dsmax = 1.E-2
      real :: qs,ns,dbzwl,ltk,ltk2,lqs,lqs2,rhos,s3m,bds,lbds,fds,kdx, &
              afas,gs1,gs4,lams,lamsmin,lamsmax,mvds,saspr,smlv,dsmm,  &
              sca0,scala,scalb,scafa,scafb,tc,ora0,orar,ora1,ora2,ora3,&
              scaw1,rads,zdr,llms,qc,epsin,epss,alphe,lamdz,nerel,     &
              scawa,scawb,scaw0,scakw,las2,lros,lasp,fdbzs,zhrox,dbzs, &
              zdrox,zdasp,fzdrs,fkdps,kdrox,kdasp,qr

      zh = 1.e-21
      zv = 1.e-21
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
            rhos = MIN(650.,EXP(75.401034-25.730216*ltk+0.81467519*lqs+&
                   2.339118*ltk2+1.9976667E-3*lqs2-0.14086434*ltk*lqs))
         ELSE
            rhos = EXP(-64808.666+23113.508*ltk-36.46632*lqs-2060.6024*&
                   ltk2-0.005729458*lqs2+6.5057411*ltk*lqs)
         ENDIF
         rhos = MIN(MAX(rhos,50.),910.)
!         smlf = MIN(1.-rhos/1000.,MAX(0.,(MIN(qr,qs)/MAX(qr,qs))**0.3))
         s3m  = 1.909859317*qs/rhos
         bds  = MIN(MAX((s3m/ns)**thrd*1.E3,dsmin*1.E3),dsmax*1.E3)
         IF (tk.GE.273.15) THEN
            fds = -0.21911541+1.2739845*bds+0.10141003*LOG(ns)+        &
                  0.30063818*bds**2.-4.3857765E-3*LOG(ns)**2.-         &
                  7.8801732E-2*bds*LOG(ns)
         ELSE
            IF (qc.GE.1.E-8) THEN
               fds = -1.1527014+2.9067645*bds+0.25316062*LOG(ns)-      &
                     0.17768557*bds**2.-0.013117292*LOG(ns)**2.-       &
                     0.17020429*bds*LOG(ns)
            ELSE
               fds = -0.2813929+1.7275463*bds+0.045550156*LOG(ns)-     &
                     0.16526226*bds**2.-1.7699916E-3*LOG(ns)**2.-      &
                     4.6441257E-2*bds*LOG(ns)
            ENDIF
         ENDIF
         kdx  = MAX(0.223,MIN(0.999,(bds/fds)**3.))
         afas = (6.*kdx-3.+SQRT(8.*kdx+1.))/(2.-2.*kdx)
         afas = MIN(MAX(afas,0.),3.E4)
         gs1  = GAMLN(afas+1.)
         gs4  = GAMLN(afas+4.)
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
         saspr = MAX(0.2,MIN(1.,-5.598876-4.5087011*ltk+5.1416616*lbds+&
                 1.0361366*ltk2+0.017687308*lbds**2.-0.9649671*ltk*lbds))
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
            zh    = MAX(1.E-28,rads*(ora1*scafb**2.+2.*ora3*scafa*     &
                    scafb+ora2*scafa**2.)/(9.*scakw))
            zv    = MAX(1.E-28,rads*(ora1*scafa**2.+2.*ora3*scafa*     &
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
            kdp   = ns*fkdps*kdrox*kdasp*dbzwl/0.11
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
      real, parameter :: thrd = 1./3., dgmin = 5.E-5, dgmax = 2.E-2
      integer :: ihail
      real :: qg,qc,qr,ng,dbzwl,ltk,ltk2,lqg,rhog,tc,mvg0,dsll,g3m,bdg,&
              lbdg,fdg,kdx,afag,gg1,gg4,lamg,lamgmin,lamgmax,mvdg,     &
              gaspr,llmg,llmg2,gmlv,dgmm,lag5,lag52,zdr,epsin,epss,    &
              alphe,lamdz,nerel,scawa,scawb,scaw0,scakw,ora0,orar,ora1,&
              ora2,ora3,sca0,scala,scalb,scaw1,scafa,scafb,radg,lag2,  &
              lrog,lasp,zhrox,fdbzg,dbzg,zdrox,zdasp,fzdrg,fkdpg,kdrox,&
              kdasp

      zh = 1.e-21
      zv = 1.e-21
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
            rhog = MIN(900.,MAX(100.,EXP(-9.57E-5*lqg**3.+3.077E-3*    &
                   lqg**2.-6.800923E-2*lqg+6.8175231)))
         ENDIF
!         gmlf = MIN(1.-rhog/1000.,MAX(0.,(MIN(qr,qg)/MAX(qr,qg))**0.3))
         g3m = 1.909859317*qg/rhog
         bdg = MIN(MAX((g3m/ng)**thrd*1.E3,dgmin*1.E3),dgmax*1.E3)
         IF (tk.GE.273.15) THEN
            fdg = 0.58006354+0.79661229*bdg-0.18394382*LOG(ng)+        &
                  0.067371044*bdg**2.+9.832945E-3*LOG(ng)**2.+         &
                  0.12433055*bdg*LOG(ng)
         ELSE
            IF (qc.GE.1.E-8) THEN
               fdg = 0.17363469+1.5044291*bdg-0.050639722*LOG(ng)+     &
                     0.015101052*bdg**2.+2.5974719E-3*LOG(ng)**2.+     &
                     0.01961464*bdg*LOG(ng)
            ELSE
               fdg = -4.8667704E-2+1.0504692*bdg+3.1159905E-3*LOG(ng)+ &
                     0.2509613*bdg**2.+1.8369028E-3*LOG(ng)**2.-       &
                     2.0083465E-2*bdg*LOG(ng)
            ENDIF
         ENDIF
         kdx  = MAX(0.223,MIN(0.999,(bdg/fdg)**3.))
         afag = (6.*kdx-3.+SQRT(8.*kdx+1.))/(2.-2.*kdx)
         afag = MIN(MAX(afag,0.),3.E4)
         gg1  = GAMLN(afag+1.)
         gg4  = GAMLN(afag+4.)
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
            gaspr = MAX(0.4,MIN(1.,-5.598876-4.5087011*ltk+5.1416616*  &
                    lbdg+1.0361366*ltk2+0.017687308*lbdg**2.-0.9649671*&
                    ltk*lbdg))
            IF (tk.LT.238.16) gaspr = 0.7
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
            zh    = MAX(1.E-28,radg*(ora1*scafb**2.+2.*ora3*scafa*     &
                    scafb+ora2*scafa**2.)/(9.*scakw))
            zv    = MAX(1.E-28,radg*(ora1*scafa**2.+2.*ora3*scafa*     &
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
            kdp   = ng*fkdpg*kdrox*kdasp*dbzwl/0.11
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

end module tcwa2_core_mod
