C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/monod/monod_forcing.F,v
C 1.6 2012/05/31 21:08:25 stephd Exp $
C $Name:  $

#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "DARWIN_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_MONOD

c=============================================================
c subroutine MONOD_forcing
c step forward bio-chemical tracers in time
C==============================================================
         SUBROUTINE MONOD_FORCING(
     U                  Ptr,
     I                  bi,bj,imin,imax,jmin,jmax,
     I                  myTime,myIter,myThid)
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#ifdef USE_QSW
#include "FFIELDS.h"
#endif
#ifdef ALLOW_LONGSTEP
#include "LONGSTEP.h"
#endif
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "GCHEM.h"
#include "MONOD_SIZE.h"
#include "MONOD.h"
#include "DARWIN_IO.h"
#include "DARWIN_FLUX.h"
#include "MONOD_FIELDS.h"

c ANNA include wavebands_params.h
#ifdef WAVEBANDS
#include "SPECTRAL_SIZE.h"
#include "SPECTRAL.h"
#include "WAVEBANDS_PARAMS.h"
#endif

C     === Global variables ===
c tracers
       _RL Ptr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,nDarwin)
       INTEGER bi,bj,imin,imax,jmin,jmax
       INTEGER myIter
       _RL myTime
       INTEGER myThid

C     !FUNCTIONS:
C     ==  Functions ==
#ifdef ALLOW_PAR_DAY
      LOGICAL  DIFF_PHASE_MULTIPLE
      EXTERNAL DIFF_PHASE_MULTIPLE
#endif

C============== Local variables
C============================================
c plankton arrays
      _RL  ZooP(nzmax)
      _RL  ZooN(nzmax)
      _RL  ZooFe(nzmax)
      _RL  ZooSi(nzmax)
      _RL  Phy(npmax)
      _RL  Phy_k(npmax,Nr)
      _RL  Phyup(npmax)
      _RL  part_k(Nr)
#ifdef ALLOW_CDOM
      _RL cdom_k(Nr)
#endif
c iron partitioning
      _RL  freefe(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
c some working variables
      _RL  sumpy
      _RL  sumpyup
c light variables
      _RL  PAR(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  sfac(1-OLy:sNy+OLy)
      _RL  atten,lite
      _RL  newtime     ! for sub-timestepping
      _RL  runtim      ! time from tracer initialization


c ANNA define variables for wavebands 
#ifdef WAVEBANDS
       integer ilam
       _RL PARw_k(tlam,Nr)
       _RL PARwup(tlam)
       _RL acdom_k(Nr,tlam)
#ifdef DAR_RADTRANS
       integer iday,iyr,imon,isec,lp,wd,mydate(4)
       _RL Edwsf(tlam),Eswsf(tlam)
       _RL Edz(tlam,Nr),Esz(tlam,Nr),Euz(tlam,Nr),Eutop(tlam,Nr)
       _RL tirrq(nr)
       _RL tirrwq(tlam,nr)
       _RL solz
       _RL rmud
       _RL actot,bctot,bbctot
       _RL apart_k(Nr,tlam),bpart_k(Nr,tlam),bbpart_k(Nr,tlam)
       _RL bt_k(Nr,tlam), bb_k(Nr,tlam)
#else
       _RL PARwdn(tlam)
#endif
C      always need for diagnostics
       _RL a_k(Nr,tlam)
#endif /* WAVEBANDS */

#ifdef DAR_DIAG_IRON
      _RL FeTbgcSource(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif

#ifdef DAR_DIAG_DIVER
      _RL  Diver1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  Diver2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  Diver3(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  Diver4(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL Shannon(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL Simpson(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)

      _RL  tmpphy(npmax)
      _RL  totphy, biotot, maxphy, phymax
#endif

#ifdef GEIDER
      _RL phychl(npmax)
      _RL phychl_k(npmax,Nr)
#ifdef DYNAMIC_CHL
      _RL dphychl(npmax)
      _RL chlup(npmax)
#endif
#endif
#ifdef ALLOW_CDOM
      _RL cdoml  
      _RL dcdoml
#endif

#ifdef ALLOW_DIAGNOSTICS
COJ for diagnostics
      _RL  PParr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  Nfixarr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
c ANNA_TAVE
#ifdef WAVES_DIAG_PCHL
      _RL  Pchlarr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,npmax)
#endif
c ANNA end TAVE
#ifdef DAR_DIAG_RSTAR
      _RL  Rstararr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,npmax)
#endif
#ifdef ALLOW_DIAZ
#ifdef DAR_DIAG_NFIXP
      _RL  NfixParr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,npmax)
#endif
#endif
#endif


      _RL  totphyC
#ifdef ALLOW_PAR_DAY
      LOGICAL itistime
      INTEGER PARiprev, PARiaccum, iperiod, nav
      _RL phase
      _RL dtsubtime
#endif
#ifdef DAR_DIAG_CHL
      _RL ChlGeiderlocal, ChlDoneylocal, ChlCloernlocal
#ifdef ALLOW_DIAGNOSTICS
      _RL GeiderChlarr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL GeiderChl2Carr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL DoneyChlarr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL DoneyChl2Carr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL CloernChlarr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL CloernChl2Carr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif
#endif
c  ballast 
      _RL freefu
      _RL inputFel
      _RL inputDustl
      _RL dusthinl
      _RL dustsinl
      _RL pocqahinl
      _RL pocqasinl
      _RL pocqainl
      INTEGER lukevar
c some local variables
       _RL PO4l
       _RL NO3l
       _RL FeTl
       _RL Sil
       _RL DOPl
       _RL DONl
       _RL DOFel
       _RL POPl
       _RL PONl
       _RL POFel
       _RL PSil
       _RL POPupl
       _RL PONupl
       _RL POFeupl
       _RL PSiupl
       _RL Tlocal
       _RL Slocal
       _RL Qswlocal
       _RL NH4l
       _RL NO2l
       _RL PARl
       _RL dzlocal
       _RL ztoplocal
       _RL zbotlocal
       _RL dz_k(Nr)
       _RL dtplankton
       _RL bottom
       _RL PP
       _RL Nfix
       _RL denit
       _RL Chl
       _RL Rstarl(npmax)
       _RL RNstarl(npmax)
#ifdef DAR_DIAG_GROW
       _RL Growl(npmax)
       _RL Growsql(npmax)
#endif
#ifdef ALLOW_DIAZ
#ifdef DAR_DIAG_NFIXP
       _RL NfixPl(npmax)
#endif
#endif

c local tendencies
       _RL dphy(npmax)
       _RL dzoop(nzmax)
       _RL dzoon(nzmax)
       _RL dzoofe(nzmax)
       _RL dzoosi(nzmax)
       _RL dPO4l
       _RL dNO3l
       _RL dFeTl
       _RL dSil
       _RL dDOPl
       _RL dDONl
       _RL dDOFel
       _RL dPOPl
       _RL dPONl
       _RL dPOFel
       _RL dPSil
       _RL dNH4l
       _RL dNO2l

#ifdef ALLOW_CARBON
       _RL dicl
       _RL docl
       _RL pocl
       _RL picl
       _RL alkl
       _RL o2l
       _RL ZooCl(nzmax)
       _RL pocupl
       _RL picupl
       _RL psihl
       _RL psihupl
       _RL pichl
       _RL pichupl
       _RL pochpl
       _RL pochpupl
       _RL pochsl
       _RL pochsupl
       _RL pocspl
       _RL pocspupl
       _RL pocssl
       _RL pocssupl
       _RL pochdl
       _RL pochdupl
       _RL pocsdl
       _RL pocsdupl
       _RL dustsl 
       _RL dustsupl
       _RL dusthl
       _RL dusthupl
c ballast production terms
       _RL dustprodinl
       _RL pocnfl
       _RL pocfl
       _RL poclprod
       _RL psilprod
       _RL piclprod
       _RL pocuplprod
       _RL psiuplprod
       _RL picuplprod
c tendencies
       _RL ddicl
       _RL ddocl
       _RL dpocl
       _RL dpicl
       _RL dalkl
       _RL do2l
       _RL dZooCl(nzmax)
       _RL dpsihl
       _RL dpichl
       _RL dpochpl
       _RL dpochsl
       _RL dpocspl 
       _RL dpocssl
       _RL dpocqastmpl
       _RL dpocqahtmpl
       _RL dpocqatmpl
       _RL dpocqasl
       _RL dpocqahl
       _RL dpochdl
       _RL dpocsdl
       _RL ddusthl
       _RL ddustsl
c ballast local tendencies
       _RL dpocprodl
       _RL dpicprodl
       _RL dpsiprodl
       _RL dpocfreel
       _RL dpocnonfreel
       _RL dpcaco3l
       _RL ddustsoftl
       _RL ddusthardl
       _RL dpocnfprodl
       _RL dpocfprodl
c ballast local up tendencies
       _RL dpocavaill
c ballast upper fluxes
       _RL pocnfupl
       _RL pocfupl 
       _RL ddustsfluxoupl
       _RL ddusthfluxoupl
c air-sea fluxes
       _RL flxCO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
       _RL flxALK(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
       _RL flxO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#endif

       _RL tot_Nfix

       _RL tmp
       _RL tmpflxpocfree
       _RL tmpflxpicsoft
       _RL tmpflxpichard
       _RL tmpflxsisoft
       _RL tmpflxsihard
       _RL tmpflxdsoft
       _RL tmpflxdhard

       _RL flxpocfree
       _RL flxpicsoft
       _RL flxpichard
       _RL flxsisoft
       _RL flxsihard
       _RL flxdustsoft
       _RL flxdusthard

       _RL phytmp, chltmp

       _RL flx, tmpflx
      INTEGER kflx

      INTEGER i,j,k,it, ktmp
      INTEGER np, nz, np2, npsave
      INTEGER debug
      CHARACTER*8 diagname

c
c
c
      DO j=1-OLy,sNy+OLy
      DO i=1-OLx,sNx+OLx
       do k=1,Nr
           freefe(i,j,k)=0. _d 0
           PAR(i,j,k) = 0. _d 0
#ifdef DAR_DIAG_DIVER
           Diver1(i,j,k)=0. _d 0
           Diver2(i,j,k)=0. _d 0
           Diver3(i,j,k)=0. _d 0
           Diver4(i,j,k)=0. _d 0
           Shannon(i,j,k)=0. _d 0
           Simpson(i,j,k)=1. _d 0
#endif

#ifdef ALLOW_DIAGNOSTICS
COJ for diagnostics
           PParr(i,j,k) = 0. _d 0
           Nfixarr(i,j,k) = 0. _d 0
#ifdef DAR_DIAG_CHL
           GeiderChlarr(i,j,k) = 0. _d 0
           GeiderChl2Carr(i,j,k) = 0. _d 0
           DoneyChlarr(i,j,k) = 0. _d 0
           DoneyChl2Carr(i,j,k) = 0. _d 0
           CloernChlarr(i,j,k) = 0. _d 0
           CloernChl2Carr(i,j,k) = 0. _d 0
#endif
c ANNA_TAVE
#ifdef WAVES_DIAG_PCHL
           DO np=1,npmax
             Pchlarr(i,j,k,np) = 0. _d 0
           ENDDO
#endif
c ANNA end TAVE
#ifdef DAR_DIAG_RSTAR
           DO np=1,npmax
             Rstararr(i,j,k,np) = 0. _d 0
           ENDDO
#endif
COJ
#ifdef ALLOW_DIAZ
#ifdef DAR_DIAG_NFIXP
           DO np=1,npmax
             NfixParr(i,j,k,np) = 0. _d 0
           ENDDO
#endif
#endif
#endif
        enddo
       ENDDO
       ENDDO
c
c bio-chemical time loop
c--------------------------------------------------
      DO it=1,nsubtime
c -------------------------------------------------
         tot_Nfix=0. _d 0
COJ cannot use dfloat because of adjoint
COJ division will be double precision anyway because of dTtracerLev
         newtime=myTime-dTtracerLev(1)+
     &           float(it)*dTtracerLev(1)/float(nsubtime)
c        print*,'it  ',it,newtime,nsubtime,myTime
         runtim=myTime-float(PTRACERS_Iter0)*dTtracerLev(1)

c determine iron partitioning  - solve for free iron
c ---------------------------
         call darwin_fe_chem(bi,bj,iMin,iMax,jMin,jMax, 
     &                       Ptr(1-OLx,1-OLy,1,bi,bj,iFeT), freefe,
     &                       myIter, mythid)
c --------------------------
#ifdef ALLOW_CARBON
c air-sea flux and dilution of CO2
         call dic_surfforcing(Ptr(1-OLx,1-OLy,1,bi,bj,iDIC),
     &                        Ptr(1-OLx,1-OLy,1,bi,bj,iALK),
     &                        Ptr(1-OLx,1-OLy,1,bi,bj,iPO4),
     &                        Ptr(1-OLx,1-OLy,1,bi,bj,iSi),
     &                        flxCO2,
     &                    bi,bj,imin,imax,jmin,jmax,
     &                    myIter,myTime,myThid)
c air-sea flux of O2
         call dic_o2_surfforcing(Ptr(1-OLx,1-OLy,1,bi,bj,iO2),
     &                       flxO2,
     &                    bi,bj,imin,imax,jmin,jmax,
     &                    myIter,myTime,myThid)
c dilusion of alkalinity
         call dic_alk_surfforcing(Ptr(1-OLx,1-OLy,1,bi,bj,iALK),
     &                       flxALK,
     &                    bi,bj,imin,imax,jmin,jmax,
     &                    myIter,myTime,myThid)
#endif

 
c find light in each grid cell
c ---------------------------
c determine incident light
#ifndef READ_PAR
#ifndef USE_QSW
         DO j=1-OLy,sNy+OLy
          sfac(j)=0. _d 0
         ENDDO
         call darwin_insol(newTime,sfac,bj)
#endif /* not USE_QSW */
#endif /* not READ_PAR */

#ifdef ALLOW_PAR_DAY
C find out which slot of PARday has previous day's average
         dtsubtime = dTtracerLev(1)/float(nsubtime)
C        running index of averaging period
C        myTime has already been incremented in this iteration,
C        go back half a substep to avoid roundoff problems
         iperiod = FLOOR((newtime-0.5 _d 0*dtsubtime)
     &                   /darwin_PARavPeriod)
C        0 -> 1, 1->2, 2->0, ...
         PARiprev = MOD(iperiod, 2) + 1
         
#ifdef   ALLOW_DIAGNOSTICS
C        always fill; this will be the same during PARavPeriod, but this
C        way it won't blow up for weird diagnostics periods.
C        we fill before updating, so the diag is the one used in this
C        time
C        step
         CALL DIAGNOSTICS_FILL(
     &         PARday(1-Olx,1-Oly,1,bi,bj,PARiprev),'PARday  ',
     &         0,Nr,2,bi,bj,myThid )
#endif
#endif /* ALLOW_PAR_DAY */

#ifdef DAR_RADTRANS
#ifndef DAR_RADTRANS_USE_MODEL_CALENDAR
#ifdef ALLOW_CAL
C get current date and time of day: iyr/imon/iday+isec
         CALL CAL_GETDATE( myIter, newtime, mydate, mythid )
         CALL CAL_CONVDATE( mydate,iyr,imon,iday,isec,lp,wd,mythid )
#else
         STOP 'need cal package or DAR_RADTRANS_USE_MODEL_CALENDAR'
#endif
#endif
#endif

C.................................................................
C.................................................................


C ========================== i,j loops =================================
       DO j=1,sNy 
        DO i=1,sNx

c ------------ these are convenient ------------------------------------
         DO k=1,Nr
          part_k(k)  = max(Ptr(i,j,k,bi,bj,iPOP),0. _d 0)
#ifdef ALLOW_CDOM
          cdom_k(k)  = max(Ptr(i,j,k,bi,bj,iCDOM),0. _d 0)
#endif
          DO np = 1,npmax
            Phy_k(np,k) = max(Ptr(i,j,k,bi,bj,iPhy+np-1),0. _d 0)
#ifdef GEIDER
#ifdef DYNAMIC_CHL
            phychl_k(np,k) = max(Ptr(i,j,k,bi,bj,iChl+np-1),0. _d 0)
#else
            phychl_k(np,k) = max(Chl_phy(i,j,k,bi,bj,np), 0. _d 0)
#endif
#endif
          ENDDO
         ENDDO

c ------------ GET CDOM_k FOR WAVEBANDS_3D and RADTRANS ----------------
#ifdef WAVEBANDS
#if defined(DAR_CALC_ACDOM) || defined(DAR_RADTRANS)
#ifdef ALLOW_CDOM
        call MONOD_ACDOM(cdom_k,
     O                     acdom_k,
     I                     myThid)
#else
         call MONOD_ACDOM(phychl_k,aphy_chl,aw,
     O                     acdom_k,
     I                     myThid)
#endif
#else
         DO k=1,Nr
          DO ilam = 1,tlam
            acdom_k(k,ilam) = acdom(ilam)
          ENDDO
         ENDDO
#endif /* DAR_CALC_ACDOM or DAR_RADTRANS */
#endif /* WAVEBANDS */

c ------------ GET INCIDENT NON-SPECTRAL LIGHT -------------------------
#if !(defined(WAVEBANDS) && defined(OASIM))
#ifdef  READ_PAR

         lite = sur_par(i,j,bi,bj)

#else /* not READ_PAR */
#ifdef USE_QSW

#ifdef ALLOW_LONGSTEP
         Qswlocal=LS_Qsw(i,j,bi,bj)
#else
         Qswlocal=Qsw(i,j,bi,bj)
#endif
         lite = -parfrac*Qswlocal*parconv*maskC(i,j,1,bi,bj)

#else /* not USE_QSW */

         lite = sfac(j)*maskC(i,j,1,bi,bj)/86400*1 _d 6

#endif /* not USE_QSW */
#endif /* not READ_PAR */

c take ice coverage into account
c unless already done in seaice package
#if !(defined (ALLOW_SEAICE) && defined (USE_QSW))
         lite = lite*(1. _d 0-fice(i,j,bi,bj))
#endif
#endif /* not(WAVEBANDS and OASIM) */

c ------------ LIGHT ATTENUATION: --------------------------------------
#ifndef WAVEBANDS
c ------------ SINGLE-BAND ATTENUATION ---------------------------------
         atten=0. _d 0
         do k=1,Nr
          if (HFacC(i,j,k,bi,bj).gt.0. _d 0) then
             sumpyup = sumpy
             sumpy = 0. _d 0
             do np=1,npmax
#ifdef GEIDER
                sumpy = sumpy + phychl_k(np,k)
#else
                sumpy = sumpy + Phy_k(np,k)
#endif
             enddo
             atten= atten + (k0 + kc*sumpy)*5. _d -1*drF(k)
             if (k.gt.1)then
               atten = atten + (k0+kc*sumpyup)*5. _d -1*drF(k-1)
             endif
             PAR(i,j,k) = lite*exp(-atten)
          endif
         enddo

#else /* WAVEBANDS */
#ifndef DAR_RADTRANS
c ------------ WAVEBANDS W/O RADTRANS ----------------------------------
         do ilam = 1,tlam
#ifdef OASIM
c add direct and diffuse, convert to uEin/m2/s/nm
          PARwup(ilam) = WtouEins(ilam)*(oasim_ed(i,j,ilam,bi,bj)+
     &                      oasim_es(i,j,ilam,bi,bj))
c and take ice fraction into account
c         PARwup(ilam) = PARwup(ilam)*(1 _d 0 - fice(i,j,bi,bj))
#else
c sf is per nm; convert to per waveband
          PARwup(ilam) = wb_width(ilam)*sf(ilam)*lite
#endif
         enddo

         do k=1,Nr
           if (HFacC(i,j,k,bi,bj).gt.0. _d 0) then
             do ilam = 1,tlam
              sumpy = 0.
              do np = 1,npmax
c get total attenuation (absorption) by phyto at each wavelength
               sumpy = sumpy + (phychl_k(np,k)*aphy_chl(np,ilam))
              enddo
c             for diagnostic
              a_k(k,ilam) = aw(ilam) + sumpy + acdom_k(k,ilam)
              atten = a_k(k,ilam)*drF(k)
              PARwdn(ilam) = PARwup(ilam)*exp(-atten)
             enddo

c find for the midpoint of the gridcell (gridcell mean)
             do ilam = 1,tlam
C             PARw_k(ilam,k)=exp((log(PARwup(ilam))+log(PARwdn(ilam)))*0.5)
              PARw_k(ilam,k)=sqrt(PARwup(ilam)*PARwdn(ilam))
             enddo

c cycle
             do ilam=1,tlam
              PARwup(ilam) = PARwdn(ilam)
             enddo
           else
             do ilam=1,tlam
              PARw_k(ilam,k) = 0. _d 0
             enddo
           endif

c sum wavebands for total PAR at the mid point of the gridcell (PARl)
           PAR(i,j,k) = 0.
           do ilam = 1,tlam
             PAR(i,j,k) = PAR(i,j,k) + PARw_k(ilam,k)
           enddo
         enddo

#else /* DAR_RADTRANS */
c ------------ FULL RADIATIVE TRANSFER CODE ----------------------------
         do ilam = 1,tlam
          Edwsf(ilam) = oasim_ed(i,j,ilam,bi,bj)
          Eswsf(ilam) = oasim_es(i,j,ilam,bi,bj)
         enddo

#ifdef DAR_RADTRANS_USE_MODEL_CALENDAR
C simplified solar zenith angle for 360-day year and daily averaged
C light
C cos(solz) is average over daylight period
         call darwin_solz360(newtime, YC(i,j,bi,bj),
     O                       solz)

#else /* not DAR_RADTRANS_USE_MODEL_CALENDAR */
C use calendar date for full solar zenith angle computation
C oj: average light effective at noon?
         solz = 0.0 _d 0
         isec = 12*3600       
         call radtrans_sfcsolz(rad,iyr,imon,iday,isec,
     I                         XC(i,j,bi,bj),YC(i,j,bi,bj),
     O                         solz)
#endif /* not DAR_RADTRANS_USE_MODEL_CALENDAR */

c have Ed,Es below surface - no need for this adjustment on Ed Es for surface affects
c         do ilam=1,tlam
c           rod(ilam) = 0.0 _d 0
c           ros(ilam) = 0.0 _d 0
c         enddo

c compute 1/cos(zenith) for direct light below surface
         call radtrans_sfcrmud(rad,solz,
     O                         rmud)

C compute absorption/scattering coefficients for radtrans
         DO k=1,Nr
          dz_k(k) = drF(k)*HFacC(i,j,k,bi,bj)
          DO ilam = 1,tlam
c           absorption by phyto
            actot = 0.0
            bctot = 0.0
            bbctot = 0.0
            DO np = 1,npmax
             actot  = actot  + phychl_k(np,k)*aphy_chl(np,ilam)
             bctot  = bctot  + phychl_k(np,k)*bphy_chl(np,ilam)
             bbctot = bbctot + phychl_k(np,k)*bbphy_chl(np,ilam)
            ENDDO
c           particulate
            apart_k(k,ilam) = part_k(k)*apart_P(ilam)
            bpart_k(k,ilam) = part_k(k)*bpart_P(ilam)
            bbpart_k(k,ilam) = part_k(k)*bbpart_P(ilam)
c           add water and CDOM
            a_k(k,ilam) = aw(ilam)+acdom_k(k,ilam)+actot+apart_k(k,ilam)
            bt_k(k,ilam) = bw(ilam) + bctot + bpart_k(k,ilam)
            bb_k(k,ilam) = darwin_bbw*bw(ilam)+bbctot+bbpart_k(k,ilam)
            bb_k(k,ilam) = MAX(darwin_bbmin, bb_k(k,ilam))
          ENDDO
         ENDDO

#ifdef DAR_RADTRANS_ITERATIVE
         call MONOD_RADTRANS_ITER(
     I                dz_k,rmud,Edwsf,Eswsf,a_k,bt_k,bb_k,
     I                darwin_radtrans_kmax,darwin_radtrans_niter,
     O                Edz,Esz,Euz,Eutop,
     O                tirrq,tirrwq,
     I                myThid)
#else
c dzlocal ?????
         call MONOD_RADTRANS(
     I                drF,rmud,Edwsf,Eswsf,a_k,bt_k,bb_k,
     O                Edz,Esz,Euz,Eutop,
     O                tirrq,tirrwq,
     I                myThid)
#endif
c
c uses chl from prev timestep (as wavebands does)
c keep like this in case need to consider upwelling irradiance as affecting the grid box above
c will pass to plankton: PARw only, but will be for this timestep for RT and prev timestep for WAVBANDS
c
c now copy
         DO k=1,Nr
           PAR(i,j,k) = tirrq(k)
           DO ilam = 1,tlam
             PARw_k(ilam,k) = tirrwq(ilam,k)
           ENDDO
         ENDDO
#endif /* DAR_RADTRANS */

c oj: ???
c so PARw and PARwup from WAVEBANDS_1D are from previous timestep (attenuation done in plankton)
c but PARw and PARwup from WAVEBANDS_3D and RADTRANS are for the current timestep

#endif /* WAVEBANDS */

C ============================ k loop ==================================
c for each layer ...   
        do k= 1, NR   
         if (HFacC(i,j,k,bi,bj).gt.0. _d 0) then

c make sure we only deal with positive definite numbers
c brute force...
             po4l  = max(Ptr(i,j,k,bi,bj,iPO4  ),0. _d 0)
             no3l  = max(Ptr(i,j,k,bi,bj,iNO3  ),0. _d 0)
             fetl  = max(Ptr(i,j,k,bi,bj,iFeT  ),0. _d 0)
             sil   = max(Ptr(i,j,k,bi,bj,iSi   ),0. _d 0)
             dopl  = max(Ptr(i,j,k,bi,bj,iDOP  ),0. _d 0)
             donl  = max(Ptr(i,j,k,bi,bj,iDON  ),0. _d 0)
             dofel = max(Ptr(i,j,k,bi,bj,iDOFe ),0. _d 0)
             DO nz = 1,nzmax
               ZooP(nz)  = max(Ptr(i,j,k,bi,bj,iZooP (nz)),0. _d 0)
               ZooN(nz)  = max(Ptr(i,j,k,bi,bj,iZooN (nz)),0. _d 0)
               ZooFe(nz) = max(Ptr(i,j,k,bi,bj,iZooFe(nz)),0. _d 0)
               ZooSi(nz) = max(Ptr(i,j,k,bi,bj,iZooSi(nz)),0. _d 0)
             ENDDO
             popl  = max(Ptr(i,j,k,bi,bj,iPOP  ),0. _d 0)
             ponl  = max(Ptr(i,j,k,bi,bj,iPON  ),0. _d 0)
             pofel = max(Ptr(i,j,k,bi,bj,iPOFe ),0. _d 0)
             psil  = max(Ptr(i,j,k,bi,bj,iPOSi ),0. _d 0)
             NH4l  = max(Ptr(i,j,k,bi,bj,iNH4  ),0. _d 0)
             NO2l  = max(Ptr(i,j,k,bi,bj,iNO2  ),0. _d 0)
#ifdef ALLOW_CDOM
             cdoml = max(Ptr(i,j,k,bi,bj,iCDOM  ),0. _d 0)
#endif
#ifdef ALLOW_CARBON
             dicl  = max(Ptr(i,j,k,bi,bj,iDIC  ),0. _d 0)
             docl  = max(Ptr(i,j,k,bi,bj,iDOC  ),0. _d 0)
             pocl  = max(Ptr(i,j,k,bi,bj,iPOC  ),0. _d 0) 
             picl  = max(Ptr(i,j,k,bi,bj,iPIC  ),0. _d 0)
             alkl  = max(Ptr(i,j,k,bi,bj,iALK  ),0. _d 0)
             o2l   = max(Ptr(i,j,k,bi,bj,iO2   ),0. _d 0)
             DO nz = 1,nzmax
               ZooCl(nz)  = max(Ptr(i,j,k,bi,bj,iZooC (nz)),0. _d 0)
             ENDDO
             psihl  = max(Ptr(i,j,k,bi,bj,iPOSih  ),0. _d 0)
             pichl  = max(Ptr(i,j,k,bi,bj,iPICh   ),0. _d 0)
             pochpl = max(Ptr(i,j,k,bi,bj,iPOChp  ),0. _d 0)
             pochsl = max(Ptr(i,j,k,bi,bj,iPOChs  ),0. _d 0)
             pocspl = max(Ptr(i,j,k,bi,bj,iPOCsp  ),0. _d 0)
             pocssl = max(Ptr(i,j,k,bi,bj,iPOCss  ),0. _d 0)
             pochdl = max(Ptr(i,j,k,bi,bj,iPOChd  ),0. _d 0)
             pocsdl = max(Ptr(i,j,k,bi,bj,iPOCsd  ),0. _d 0)
             dustsl = max(Ptr(i,j,k,bi,bj,iDUSTs  ),0. _d 0)
             dusthl = max(Ptr(i,j,k,bi,bj,iDUSTh  ),0. _d 0)
c ballast local production
             pocnfl = max(Ptr(i,j,k,bi,bj,iPOC ), 0. _d 0)
             pocfl  = max(Ptr(i,j,k,bi,bj,iPOC ), 0. _d 0)
             poclprod = max(Ptr(i,j,k,bi,bj,iPOC  ),0. _d 0)
             psilprod = max(Ptr(i,j,k,bi,bj,iPOC ),0. _d 0)
             piclprod = max(Ptr(i,j,k,bi,bj,iPOC ),0. _d 0)
#endif

             totphyC = 0. _d 0
             DO np=1,npmax
              totphyC = totphyC + R_PC(np)*Ptr(i,j,k,bi,bj,iPhy+np-1)
             ENDDO

             DO np = 1,npmax
               Phy(np) = Phy_k(np,k)
#ifdef GEIDER
               phychl(np) = phychl_k(np,k)
#endif
             ENDDO

#ifdef DAR_DIAG_DIVER
            Diver1(i,j,k)=0. _d 0
            Diver2(i,j,k)=0. _d 0
            Diver3(i,j,k)=0. _d 0
            Diver4(i,j,k)=0. _d 0
            totphy=0. _d 0
            do np=1,npmax
              totphy=totphy + Phy(np)
              tmpphy(np)=Phy(np)
            enddo
            if (totphy.gt.diver_thresh0) then
              do np=1,npmax
c simple threshhold
               if (Phy(np).gt.diver_thresh1) then
                 Diver1(i,j,k)=Diver1(i,j,k)+1. _d 0
               endif
c proportion of total biomass
               if (Phy(np)/totphy.gt.diver_thresh2) then
                 Diver2(i,j,k)=Diver2(i,j,k)+1. _d 0
               endif
              enddo
c majority of biomass by finding rank order
               biotot=0. _d 0
               do np2=1,npmax
                 phymax=0. _d 0
                 do np=1,npmax
                   if (tmpphy(np).gt.phymax) then
                     phymax=tmpphy(np)
                     npsave=np
                   endif
                 enddo
                 if (biotot.lt.totphy*diver_thresh3) then
                    Diver3(i,j,k)=Diver3(i,j,k)+1. _d 0
                 endif
                 biotot=biotot+tmpphy(npsave)
                 tmpphy(npsave)=0. _d 0
                 if (np2.eq.1) then
                    maxphy=phymax
                 endif
               enddo
c ratio of maximum species
               do np=1,npmax
                 if (Phy(np).gt.diver_thresh4*maxphy) then
                    Diver4(i,j,k)=Diver4(i,j,k)+1. _d 0
                 endif
               enddo
c           totphy > thresh0
            endif
c Shannon and Simpson indices
            Shannon(i,j,k) = 0. _d 0
c note: minimal valid value is 1, but we set to zero below threshold
            Simpson(i,j,k) = 0. _d 0
            if (totphy.gt.shannon_thresh) then 
              do np=1,npmax
               if (Phy(np) .gt. 0. _d 0) then
                tmpphy(np) = Phy(np)/totphy
                Shannon(i,j,k)=Shannon(i,j,k)+tmpphy(np)*LOG(tmpphy(np))
                Simpson(i,j,k)=Simpson(i,j,k)+tmpphy(np)*tmpphy(np)
               endif
              enddo
              Shannon(i,j,k) = -Shannon(i,j,k)
              Simpson(i,j,k) = 1./Simpson(i,j,k)
            endif
#endif

c..........................................................
c find local light
c..........................................................

             PARl = PAR(i,j,k)
c..........................................................

c for explicit sinking of particulate matter and phytoplankton
             if (k.eq.1) then
                popupl =0. _d 0
                ponupl =0. _d 0
                pofeupl = 0. _d 0
                psiupl = 0. _d 0
                do np=1,npmax
                  Phyup(np)=0. _d 0
#ifdef DYNAMIC_CHL
                  chlup(np)=0. _d 0
#endif
                enddo
#ifdef ALLOW_CARBON
                pocupl   = 0. _d 0
                picupl   = 0. _d 0
                psihupl  = 0. _d 0 
                pichupl  = 0. _d 0
                pochpupl = 0. _d 0
                pochsupl = 0. _d 0
                pocspupl = 0. _d 0
                pocssupl = 0. _d 0
                pochdupl = 0. _d 0 
                pocsdupl = 0. _d 0
                dustsupl = 0. _d 0
                dusthupl = 0. _d 0
c ballast intitalize production
                pocuplprod = 0. _d 0
                psiuplprod = 0. _d 0
                picuplprod = 0. _d 0
                pocnfupl = 0. _d 0
                pocfupl  = 0. _d 0  
#endif
             endif

#ifdef ALLOW_LONGSTEP
             Tlocal = LS_theta(i,j,k,bi,bj)
             Slocal = LS_salt(i,j,k,bi,bj)
#else
             Tlocal = theta(i,j,k,bi,bj)
             Slocal = salt(i,j,k,bi,bj)
#endif

c             dzlocal = drF(k)*HFacC(i,j,k,bi,bj)
             dzlocal = drF(k)
             ztoplocal = rF(k)
             zbotlocal = rF(k+1)
             freefu = max(freefe(i,j,k),0. _d 0)
             if (k.eq.1) then
               inputFel = inputFe(i,j,bi,bj)
c----------------------------------------------------------
c Initialize dust deposiation as production in top cell
c used in ballast model - L. gloege 2016
c----------------------------------------------------------
               lukevar     = 1  ! 1 if dustPOC > POC avail 0 if dustPOC < POCavail
               inputDustl  = inputDust(i,j,bi,bj)*1000    ! Dust deposition (mg/m^2/s)
               dustprodinl = (inputDustl)/dzlocal   ! Dust "production"
               dusthinl    = (1-f_dust_bio)*(f_dust_hard) *dustprodinl 
               dustsinl    = (1-f_dust_bio)*(1-f_dust_hard)*dustprodinl 
               pocqahinl   = omega_dust * dusthinl  ! POC qual. assoc. with h dust
               pocqasinl   = omega_dust * dustsinl  ! POC qual. assoc. with s dust
               pocqainl    = pocqasinl  + pocqahinl
c               WRITE( *,'(A,D)') 'lukeDUSThin ', dusthinl
c               WRITE( *,'(A,D)') 'lukeDUSTsin ', dustsinl
c               WRITE( *,'(A,D)') 'lukePOCqahin ', pocqahinl
c               WRITE( *,'(A,D)') 'lukePOCqasin ', pocqasinl
c               WRITE( *,'(A,D)') 'lukePOCqain ', pocqainl
c----------------------------------------------------------
             else
               inputFel = 0. _d 0
               inputDustl = 0. _d 0
             endif

c            WRITE( *,'(A,D)') 'lukeDZ ', dzlocal
            

c             dzlocal = drF(k)*HFacC(i,j,k,bi,bj)
c             ztoplocal = rF(k)
c             zbotlocal = rF(k+1)

c set bottom=1.0 if the layer below is not ocean
             ktmp=min(nR,k+1)
             if(hFacC(i,j,ktmp,bi,bj).eq.0. _d 0.or.k.eq.Nr) then
               bottom = 1.0 _d 0
             else
               bottom = 0.0 _d 0
             endif

c set tendencies to 0
             do np=1,npmax
               dphy(np)=0. _d 0
             enddo
             do nz=1,nzmax
               dzoop(nz)=0. _d 0
               dzoon(nz)=0. _d 0
               dzoofe(nz)=0. _d 0
               dzoosi(nz)=0. _d 0
             enddo
             dPO4l=0. _d 0
             dNO3l=0. _d 0
             dFeTl=0. _d 0
             dSil=0. _d 0
             dDOPl=0. _d 0
             dDONl=0. _d 0
             dDOFel=0. _d 0
             dPOPl=0. _d 0
             dPONl=0. _d 0
             dPOFel=0. _d 0
             dPSil=0. _d 0
             dNH4l=0. _d 0
             dNO2l=0. _d 0
#ifdef DYNAMIC_CHL
             do np=1,npmax
               dphychl(np)=0. _d 0
             enddo
#endif
#ifdef ALLOW_CDOM
            dcdoml=0. _d 0
#endif
#ifdef ALLOW_CARBON
             ddicl=0. _d 0
             ddocl=0. _d 0
             dpocl=0. _d 0
             dpicl=0. _d 0
             dalkl=0. _d 0
             do2l=0. _d 0
             do nz=1,nzmax
               dzoocl(nz)=0. _d 0
             enddo
             dpsihl=0. _d 0
             dpichl=0. _d 0
             dpochpl=0. _d 0
             dpochsl=0. _d 0
             dpocspl=0. _d 0
             dpocssl=0. _d 0
             dpocqastmpl=0. _d 0
             dpocqahtmpl=0. _d 0
             dpocqatmpl=0. _d 0
             dpocqasl=0. _d 0
             dpocqahl=0. _d 0
             dpochdl=0. _d 0
             dpocsdl=0. _d 0
             ddusthl=0. _d 0
             ddustsl=0. _d 0
c ballast stuff
             dpocprodl=0. _d 0
             dpicprodl=0. _d 0
             dpsiprodl=0. _d 0
             dpcaco3l=0. _d 0
             ddustsoftl=0. _d 0
             ddusthardl=0. _d 0
             dpocfreel=0. _d 0
             dpocnonfreel=0. _d 0
c 
             dpocnfprodl=0. _d 0
             dpocfprodl=0. _d 0
             dpocavaill=0. _d 0
#endif
c set other arguments to zero
             PP=0. _d 0
             Nfix=0. _d 0
             denit=0. _d 0
             do np=1,npmax
                Rstarl(np)=0. _d 0
                RNstarl(np)=0. _d 0
#ifdef DAR_DIAG_GROW
                Growl(np)=0.  _d 0
                Growsql(np)=0. _d 0
#endif
#ifdef ALLOW_DIAZ
#ifdef DAR_DIAG_NFIXP
                NfixPl(np)=0. _d 0
#endif
#endif
             enddo


            debug=0
c           if (i.eq.20.and.j.eq.20.and.k.eq.1) debug=8
c           if (i.eq.10.and.j.eq.10.and.k.eq.1) debug=100
c           if (i.eq.1.and.j.eq.10.and.k.eq.1) debug=10
c           if (i.eq.1.and.j.eq.1.and.k.eq.10) debug=14
 
            if (debug.eq.7) print*,'PO4, DOP, POP, ZooP',
     &                       PO4l,  DOPl, POPl, zooP
            if (debug.eq.7) print*,'NO3, NO2, NH4, DON, PON, ZooN', 
     &                       NO3l,NO2l,NH4l, DONl, PONl, ZooN
            if (debug.eq.7) print*,'FeT, DOFe, POFe, Zoofe', 
     &                       FeTl,  DOFel, POFel, zooFe
            if (debug.eq.7) print*,'Si, Psi, zooSi', 
     &                       Sil, PSil, zooSi
            if (debug.eq.7) print*,'Total Phy', sumpy, PARl, lite
            if (debug.eq.7) print*,'Phy', Phy

            if (debug.eq.8) print*,'k, PARl, inputFel, dzlocal',
     &                        PARl, inputFel, dzlocal

c           if (NO3l.eq.0. _d 0.or.NO2l.eq.0. _d 0
c    &          .or.NH4l.eq.0. _d 0) then
c            print*,'QQ N zeros',i,j,k,NO3l,NO2l,NH4l
c           endif


c ANNA pass extra variables if WAVEBANDS
             CALL MONOD_PLANKTON(
     U                       Phy,
     I                       zooP, zooN, zooFe, zooSi,
     O                       PP, Chl, Nfix, denit, 
     I                       PO4l, NO3l, FeTl, Sil,
     I                       NO2l, NH4l, 
     I                       DOPl, DONl, DOFel,
     I                       POPl, PONl, POFel, PSil,
     I                       phyup, popupl, ponupl, 
     I                       pofeupl, psiupl,
     I                       PARl,
     I                       Tlocal, Slocal,
     I                       freefu, inputFel, inputDustl,
     I                       bottom, dzlocal, 
     I                       ztoplocal, zbotlocal,
     I                       dusthinl, dustsinl,
     I                       pocqahinl, pocqasinl,
     I                       pocqainl,
     I                       lukevar,
     O                       Rstarl, RNstarl,
#ifdef DAR_DIAG_GROW
     O                       Growl, Growsql,
#endif
#ifdef ALLOW_DIAZ
#ifdef DAR_DIAG_NFIXP
     O                       NfixPl,
#endif
#endif
     O                       dphy, dzooP, dzooN, dzooFe,
     O                       dzooSi,
     O                       dPO4l, dNO3l, dFeTl, dSil,
     O                       dNH4l, dNO2l,
     O                       dDOPl, dDONl, dDOFel, 
     O                       dPOPl, dPONl, dPOFel, dPSil,
#ifdef ALLOW_CARBON
     I                       dicl, docl, pocl, picl,
     I                       pocnfl, pocfl, poclprod,
     I                       alkl, o2l, zoocl,
     I                       pocupl, picupl, pocuplprod,
     I                       psilprod, piclprod,
     I                       psiuplprod, picuplprod,
     I                       pocnfupl,
     I                       pocfupl,
     I                       ddustsfluxoupl,
     I                       ddusthfluxoupl,
     I                       psihl, psihupl,
     I                       pichl, pichupl,
     I                       pochpl, pochpupl,
     I                       pochsl, pochsupl,
     I                       pocspl, pocspupl,
     I                       pocssl, pocssupl,
     I                       pochdl, pochdupl,
     I                       pocsdl, pocsdupl,
     I                       dustsl, dustsupl,
     I                       dusthl, dusthupl,
     O                       ddicl, ddocl, dpocl, dpicl,
     O                       dalkl, do2l, dzoocl, 
     O                       dpocnonfreel, phy1pocl, phy2pocl,
     O                       dpocfreel, dpsio3l, dpcaco3l,
     O                       ddustsoftl, ddusthardl, ddustl,
     O                       dpocprodl, dpicprodl, dpsiprodl,
     O                       dpocnfprodl, dpocfprodl,
     O                       dpocavaill,
     O                       dpsihl, dpichl,
     O                       dpochpl, dpochsl,
     O                       dpocspl, dpocssl,
     O                       dpocqastmpl, dpocqahtmpl, dpocqatmpl,
     O                       dpocqasl, dpocqahl,
     O                       dpochdl, dpocsdl,
     O                       ddusthl, ddustsl,
#endif
#ifdef GEIDER
     O                       phychl,
#ifdef DYNAMIC_CHL
     I                       dphychl,
     I                       chlup,
#endif/
#ifdef ALLOW_CDOM
     O                       dcdoml,
     I                       cdoml,
#endif
#ifdef WAVEBANDS
     I                       PARw_k(1,k),
#endif
#endif
#ifdef ALLOW_PAR_DAY
     I                       PARday(i,j,k,bi,bj,PARiprev),
#endif
#ifdef DAR_DIAG_CHL
     O                       ChlGeiderlocal, ChlDoneylocal,
     O                       ChlCloernlocal,
#endif
     I                       debug,
     I                       runtim,
     I                       MyThid) 

c
c            if (i.eq.1.and.k.eq.1.and.j.eq.5) then
c              print*,i,j,k
c              print*,'NO3,No2,NH4',  NO3l, NO2l, NH4l
c              print*,'dNO3 etc',dNO3l,dNH4l, dNO2l
c              print*,'PO4',PO4l,dPO4l
c            endif
c
#ifdef IRON_SED_SOURCE
c only above minimum depth (continental shelf)
             if (rF(k).gt.-depthfesed) then
c only if bottom layer
               if (bottom.eq.1.0 _d 0) then
#ifdef IRON_SED_SOURCE_VARIABLE
c calculate sink of POP into bottom layer
                tmp=(wp_sink*POPupl)/(dzlocal)
c convert to dPOCl
                dFetl=dFetl+fesedflux_pcm*(tmp*106. _d 0)
#else
                 dFetl=dFetl+fesedflux/
     &                   (drF(k)*hFacC(i,j,k,bi,bj))
#endif
               endif
             endif
#endif


             popupl = POPl
             ponupl = PONl
             pofeupl = POFel
             psiupl = PSil
             do np=1,npmax
               Phyup(np) = Phy(np) 
#ifdef DYNAMIC_CHL
               chlup(np) = phychl(np)
#endif
             enddo

c
c these are the up local terms
c note: not case senstive (POCl == pocl)
c you added the production terms
c luke gloege (ballast)
c
#ifdef ALLOW_CARBON
             pocupl   = POCl
             picupl   = PICl
             psihupl  = PSIhl
             pichupl  = PIChl
             pochpupl = POChpl
             pochsupl = POChsl
             pocspupl = POCspl
             pocssupl = POCssl
             pocuplprod = POClprod
             psiuplprod = PSIlprod
             picuplprod = PIClprod
             pochdupl = POChdl
             pocsdupl = POCsdl

             if (k.eq.1) then
                dustsupl = dustsinl
                dusthupl = dusthinl
             else
                dustsupl = DUSTsl
                dusthupl = DUSThl
             endif

c include surface forcing
             if (k.eq.1) then
              ddicl  =  ddicl  + flxCO2(i,j)
              dalkl  =  dalkl  + flxALK(i,j)
              do2l   =  do2l   + flxO2(i,j)
             endif
             pocnfupl = POCnfl
             pocfupl  =  POCfl
c ballast fluxes tendencies
             ddustsfluxoupl = ddustsoftl
             ddusthfluxoupl = ddusthardl
#endif
c
#ifdef CONS_SUPP
c only works for two layer model
            if (k.eq.2) then
              dpo4l=0. _d 0
              dno3l=0. _d 0
              dfetl=0. _d 0
              dsil=0. _d 0
            endif
#endif
#ifdef RELAX_NUTS
#ifdef DENIT_RELAX
           if (rF(k).lt.-depthdenit) then
            if (darwin_relaxscale.gt.0. _d 0) then
             IF ( darwin_NO3_RelaxFile .NE. ' '  ) THEN
c Fanny's formulation
              tmp=(Ptr(i,j,k,bi,bj,iNO3 )-no3_obs(i,j,k,bi,bj))
              if (tmp.gt.0. _d 0) then
                dno3l=dno3l-(tmp/
     &                      darwin_relaxscale)
                denit=tmp/
     &                      darwin_relaxscale
              else
                denit=0. _d 0
              endif
c --- end fanny's formulation
             ENDIF
c steph's alternative
c            tmp=(Ptr(i,j,k,bi,bj,iNO3 )-
c    &                16. _d 0 * Ptr(i,j,k,bi,bj,iPO4 ))
c            if (tmp.gt.0. _d 0) then
c               dno3l=dno3l-(tmp/
c    &                      darwin_relaxscale)
c               denit=tmp/
c    &                      darwin_relaxscale
c             else
c               denit=0. _d 0
c             endif
c ---- end steph's alternative
            endif
           endif
#else
c -------------------------
c hard relax back to observed
c if deeper than euphotic depth
c 
c L. Gloege May 2015
c --------------------------

            if (rF(k).lt.zeu) then 
            if (darwin_relaxscale.gt.0. _d 0) then
             IF ( darwin_PO4_RelaxFile .NE. ' '  ) THEN
              tmp=(Ptr(i,j,k,bi,bj,iPO4 )-po4_obs(i,j,k,bi,bj))
              if (tmp.ne.0. _d 0) then
                dpo4l=dpo4l-(tmp/
     &                      darwin_relaxscale)
              endif
             ENDIF
             IF ( darwin_NO3_RelaxFile .NE. ' '  ) THEN
              tmp=(Ptr(i,j,k,bi,bj,iNO3 )-no3_obs(i,j,k,bi,bj))
              if (tmp.ne.0. _d 0) then
                dno3l=dno3l-(tmp/
     &                      darwin_relaxscale)
              endif
             ENDIF
             IF ( darwin_Fet_RelaxFile .NE. ' '  ) THEN
              tmp=(Ptr(i,j,k,bi,bj,iFeT )-fet_obs(i,j,k,bi,bj))
              if (tmp.ne.0. _d 0) then
                dfetl=dfetl-(tmp/
     &                      darwin_relaxscale)
              endif
             ENDIF
             IF ( darwin_Si_RelaxFile .NE. ' '  ) THEN
              tmp=( Ptr(i,j,k,bi,bj,iSi  )-si_obs(i,j,k,bi,bj))
              if (tmp.ne.0. _d 0) then
                dsil=dsil-(tmp/
     &                    darwin_relaxscale)
              endif
             ENDIF
            endif
            endif

c --------------------
c do not go below observed 
c between surface and 
c euphotic depth
c
c tmp.lt.0 (semi relax)
c tmp.ne.0 (hard relax)
c
c L. Gloege May 2015
c ---------------------

            if (rF(k).gt.zeu) then 
            if (darwin_relaxscale.gt.0. _d 0) then
             IF ( darwin_PO4_RelaxFile .NE. ' '  ) THEN
              tmp=(Ptr(i,j,k,bi,bj,iPO4 )-po4_obs(i,j,k,bi,bj))
              if (tmp.lt.0. _d 0) then
                dpo4l=dpo4l-(tmp/
     &                      darwin_relaxscale)
              endif
             ENDIF
             IF ( darwin_NO3_RelaxFile .NE. ' '  ) THEN
              tmp=(Ptr(i,j,k,bi,bj,iNO3 )-no3_obs(i,j,k,bi,bj))
              if (tmp.lt.0. _d 0) then
                dno3l=dno3l-(tmp/
     &                      darwin_relaxscale)
              endif
             ENDIF
             IF ( darwin_Fet_RelaxFile .NE. ' '  ) THEN
              tmp=(Ptr(i,j,k,bi,bj,iFeT )-fet_obs(i,j,k,bi,bj))
              if (tmp.lt.0. _d 0) then
                dfetl=dfetl-(tmp/
     &                      darwin_relaxscale)
              endif
             ENDIF
             IF ( darwin_Si_RelaxFile .NE. ' '  ) THEN
              tmp=( Ptr(i,j,k,bi,bj,iSi  )-si_obs(i,j,k,bi,bj))
              if (tmp.lt.0. _d 0) then
                dsil=dsil-(tmp/
     &                    darwin_relaxscale)
              endif
             ENDIF
            endif
            endif
#endif
#endif
#ifdef FLUX_NUTS
                dpo4l=dpo4l+po4_flx(i,j,k,bi,bj)
                dno3l=dno3l+no3_flx(i,j,k,bi,bj)
                dfetl=dfetl+fet_flx(i,j,k,bi,bj)
                dsil=dsil+si_flx(i,j,k,bi,bj)
#endif
c
c now update main tracer arrays
          dtplankton = PTRACERS_dTLev(k)/float(nsubtime)
          Ptr(i,j,k,bi,bj,iPO4 ) = Ptr(i,j,k,bi,bj,iPO4) +
     &                                  dtplankton*dpo4l
          Ptr(i,j,k,bi,bj,iNO3 ) = Ptr(i,j,k,bi,bj,iNO3) +
     &                                  dtplankton*dno3l
          Ptr(i,j,k,bi,bj,iFeT ) = Ptr(i,j,k,bi,bj,iFeT) +
     &                                  dtplankton*dfetl
          Ptr(i,j,k,bi,bj,iSi  ) = Ptr(i,j,k,bi,bj,iSi ) +
     &                                  dtplankton*dsil
          Ptr(i,j,k,bi,bj,iDOP ) = Ptr(i,j,k,bi,bj,iDOP) +
     &                                  dtplankton*ddopl
          Ptr(i,j,k,bi,bj,iDON ) = Ptr(i,j,k,bi,bj,iDON) +
     &                                  dtplankton*ddonl
          Ptr(i,j,k,bi,bj,iDOFe) = Ptr(i,j,k,bi,bj,iDOFe) +
     &                                  dtplankton*ddofel
          Ptr(i,j,k,bi,bj,iPOP ) = Ptr(i,j,k,bi,bj,iPOP ) +
     &                                  dtplankton*dpopl
          Ptr(i,j,k,bi,bj,iPON ) = Ptr(i,j,k,bi,bj,iPON ) +
     &                                  dtplankton*dponl
          Ptr(i,j,k,bi,bj,iPOFe) = Ptr(i,j,k,bi,bj,iPOFe) +
     &                                  dtplankton*dpofel
          Ptr(i,j,k,bi,bj,iPOSi) = Ptr(i,j,k,bi,bj,iPOSi) +
     &                                  dtplankton*dpsil
          Ptr(i,j,k,bi,bj,iNH4 ) = Ptr(i,j,k,bi,bj,iNH4 ) +
     &                                  dtplankton*dnh4l
          Ptr(i,j,k,bi,bj,iNO2 ) = Ptr(i,j,k,bi,bj,iNO2 ) +
     &                                  dtplankton*dno2l
          DO nz = 1,nzmax
           Ptr(i,j,k,bi,bj,iZooP (nz)) = Ptr(i,j,k,bi,bj,iZooP (nz)) +
     &                                        dtplankton*dzoop (nz)
           Ptr(i,j,k,bi,bj,iZooN (nz)) = Ptr(i,j,k,bi,bj,iZooN (nz)) +
     &                                        dtplankton*dzoon (nz)
           Ptr(i,j,k,bi,bj,iZooFe(nz)) = Ptr(i,j,k,bi,bj,iZooFe(nz)) +
     &                                        dtplankton*dzoofe(nz)
           Ptr(i,j,k,bi,bj,iZooSi(nz)) = Ptr(i,j,k,bi,bj,iZooSi(nz)) +
     &                                        dtplankton*dzoosi(nz)
          ENDDO
          DO np = 1,npmax
           Ptr(i,j,k,bi,bj,iPhy+np-1) = Ptr(i,j,k,bi,bj,iPhy+np-1) +
     &                                       dtplankton*dPhy(np)
#ifdef GEIDER
#ifdef DYNAMIC_CHL
           if (np.eq.1) Chl=0. _d 0
           Ptr(i,j,k,bi,bj,iChl+np-1) = Ptr(i,j,k,bi,bj,iChl+np-1) +
     &                                       dtplankton*dphychl(np)
c          chltmp=Ptr(i,j,k,bi,bj,iChl+np-1)
c          phytmp=Ptr(i,j,k,bi,bj,iPhy+np-1)
c          Ptr(i,j,k,bi,bj,iChl+np-1)=
c    &         max(chltmp,phytmp*R_PC(np)*chl2cmin(np))
c          if (np.eq.1.and.i.eq.1.and.j.eq.1.and.k.eq.1)
c    &     print*,chltmp,phytmp,phytmp*R_PC(np)*chl2cmin(np),
c    &                   phytmp*R_PC(np)*chl2cmax(np)
c in darwin_plankton this is stored for previous timestep. Reset here.
           Chl=Chl+Ptr(i,j,k,bi,bj,iChl+np-1)
#else
           Chl_phy(i,j,k,bi,bj,np)=phychl(np)
#endif
#endif
         ENDDO
#ifdef ALLOW_CDOM
          Ptr(i,j,k,bi,bj,iCDOM ) = Ptr(i,j,k,bi,bj,iCDOM ) +
     &                                  dtplankton*dcdoml
#endif
#ifdef ALLOW_CARBON
          Ptr(i,j,k,bi,bj,iDIC ) = Ptr(i,j,k,bi,bj,iDIC ) +
     &                                  dtplankton*ddicl
          Ptr(i,j,k,bi,bj,iDOC ) = Ptr(i,j,k,bi,bj,iDOC ) +
     &                                  dtplankton*ddocl
          Ptr(i,j,k,bi,bj,iPOC ) = Ptr(i,j,k,bi,bj,iPOC ) +
     &                                  dtplankton*dpocl

          Ptr(i,j,k,bi,bj,iPIC ) = Ptr(i,j,k,bi,bj,iPIC ) +
     &                                  dtplankton*dpicl
          Ptr(i,j,k,bi,bj,iALK ) = Ptr(i,j,k,bi,bj,iALK ) +
     &                                  dtplankton*dalkl
          Ptr(i,j,k,bi,bj,iO2  ) = Ptr(i,j,k,bi,bj,iO2  ) +
     &                                  dtplankton*do2l
          DO nz = 1,nzmax
           Ptr(i,j,k,bi,bj,iZooC (nz)) = Ptr(i,j,k,bi,bj,iZooC (nz)) +
     &                                        dtplankton*dzoocl (nz)
          ENDDO

c----------------------------------------------------------------------
c ballast portion of code
c the term after dtplankton needs to be changed!!
c----------------------------------------------------------------------

          Ptr(i,j,k,bi,bj,iPOSih) = Ptr(i,j,k,bi,bj,iPOSih) +
     &                                  dtplankton*dpsihl
     
          Ptr(i,j,k,bi,bj,iPICh) = Ptr(i,j,k,bi,bj,iPICh) +
     &                                  dtplankton*dpichl

          Ptr(i,j,k,bi,bj,iPOChp ) = Ptr(i,j,k,bi,bj,iPOChp ) +
     &                                  dtplankton*dpochpl
 
          Ptr(i,j,k,bi,bj,iPOChs ) = Ptr(i,j,k,bi,bj,iPOChs ) +
     &                                  dtplankton*dpochsl

          Ptr(i,j,k,bi,bj,iPOCsp ) = Ptr(i,j,k,bi,bj,iPOCsp ) +
     &                                  dtplankton*dpocspl

          Ptr(i,j,k,bi,bj,iPOCss ) = Ptr(i,j,k,bi,bj,iPOCss ) +
     &                                  dtplankton*dpocssl

c dust stuff
          Ptr(i,j,k,bi,bj,iPOChd ) = Ptr(i,j,k,bi,bj,iPOChd ) +
     &                                  dtplankton*dpochdl

          Ptr(i,j,k,bi,bj,iPOCsd ) = Ptr(i,j,k,bi,bj,iPOCsd ) +
     &                                  dtplankton*dpocsdl

          Ptr(i,j,k,bi,bj,iDUSTh ) = Ptr(i,j,k,bi,bj,iDUSTh ) +
     &                                  dtplankton*ddusthl

          Ptr(i,j,k,bi,bj,iDUSTs ) = Ptr(i,j,k,bi,bj,iDUSTs ) +
     &                                  dtplankton*ddustsl

          Ptr(i,j,k,bi,bj,iPOCflux ) = 
     &                       (Ptr(i,j,k,bi,bj,iPOChp) * wpic_sink) +
     &                       (Ptr(i,j,k,bi,bj,iPOChs) * wsi_sink ) +
     &                       (Ptr(i,j,k,bi,bj,iPOCsp) * wpic_sink) +
     &                       (Ptr(i,j,k,bi,bj,iPOCss) * wsi_sink ) +
     &                       (Ptr(i,j,k,bi,bj,iPOCsd) * wc_sink  ) + 
     &                       (Ptr(i,j,k,bi,bj,iPOChd) * wc_sink  ) +
     &                       (Ptr(i,j,k,bi,bj,iPOC  ) * wc_sink  ) 
#endif
c
#ifdef ALLOW_MUTANTS
cQQQQTEST
          if (debug.eq.11) then
          if (k.lt.8) then
          do np=1,60
          if(mod(np,4).eq. 1. _d 0)then
          np2=np+1
          np4=np+3

Coj: couldn't test this part after change Phynp -> Ptr(...,iPhy+np-1)
Coj: used to be many copies of this:
C           if (dPhy(2).gt.dPhy(4).and.dPhy(4).gt.0. _d 0) then
C             print*,'QQQ dphy2 > dphy4',i,j,k,Phy2(i,j,k),
C    &             Phy4(i,j,k), dPhy(2), dPhy(4)
C           endif
C           if (Phy2(i,j,k).gt.Phy4(i,j,k).and.
C    &                           Phy4(i,j,k).gt.0. _d 0) then
C             print*,'QQ phy02 > phy04',i,j,k,Phy2(i,j,k),
C    &             Phy4(i,j,k), dPhy(2), dPhy(4)
C           endif

          if (dPhy(np2).gt.dPhy(np4).and.dPhy(np4).gt.0. _d 0) then
             print*,'QQQ dphy',np2,' > dphy',np4,i,j,k,Phy2(i,j,k),
     &             Ptr(i,j,k,bi,bj,iPhy+np4-1), dPhy(2), dPhy(4)
          endif
          if (Ptr(i,j,k,bi,bj,iphy+np2-1).gt.Ptr(i,j,k,bi,bj,iPhy+np4-1)
     &        .and. Ptr(i,j,k,bi,bj,iPhy+np4-1).gt.0. _d 0) then
            print*,'QQ phy',np2,' > ',np4,i,j,k,
     &             Ptr(i,j,k,bi,bj,iPhy+np2-1),
     &             Ptr(i,j,k,bi,bj,iPhy+np4-1), dPhy(2), dPhy(4)
          endif

          endif
          enddo  ! np
          endif  ! k
          endif
#endif

#ifdef ALLOW_DIAGNOSTICS
COJ for diagnostics
             PParr(i,j,k) = PP
             Nfixarr(i,j,k) = Nfix
c ANNA_TAVE
#ifdef WAVES_DIAG_PCHL
             DO np = 1,npmax
               Pchlarr(i,j,k,np) = phychl(np)
             ENDDO
#endif
c ANNA end TAVE
#ifdef DAR_DIAG_RSTAR
             DO np = 1,npmax
               Rstararr(i,j,k,np) = Rstarl(np)
             ENDDO
#endif
#ifdef ALLOW_DIAZ
#ifdef DAR_DIAG_NFIXP
             DO np = 1,npmax
               NfixParr(i,j,k,np) = NfixPl(np)
             ENDDO
#endif
#endif
#ifdef DAR_DIAG_CHL
             GeiderChlarr(i,j,k) = ChlGeiderlocal
             DoneyChlarr(i,j,k) = ChlDoneylocal
             CloernChlarr(i,j,k) = ChlCloernlocal
             IF (totphyC .NE. 0. _d 0) THEN
               GeiderChl2Carr(i,j,k) = ChlGeiderlocal/totphyC
               DoneyChl2Carr(i,j,k) = ChlDoneylocal/totphyC
               CloernChl2Carr(i,j,k) = ChlCloernlocal/totphyC
             ELSE
               GeiderChl2Carr(i,j,k) = 0. _d 0
               DoneyChl2Carr(i,j,k) = 0. _d 0
               CloernChl2Carr(i,j,k) = 0. _d 0
             ENDIF
#endif
COJ
#ifdef DAR_DIAG_IRON
             FeTbgcSource(i,j,k) = dfetl
#endif
#endif /* ALLOW_DIAGNOSTICS */

c total fixation (NOTE - STILL NEEDS GLOB SUM)
             tot_Nfix=tot_Nfix+ 
     &             Nfix*rA(i,j,bi,bj)*rF(k)*hFacC(i,j,k,bi,bj)

#ifdef ALLOW_TIMEAVE
c save averages
c            Phygrow1ave(i,j,k,bi,bj)=Phygrow1ave(i,j,k,bi,bj)+
c    &                             mu1*py1*deltaTclock
c    &                       /float(nsubtime)
c            Phygrow2ave(i,j,k,bi,bj)=Phygrow2ave(i,j,k,bi,bj)+
c    &                             mu2*py2*deltaTclock
c    &                       /float(nsubtime)
c            Zoograzave(i,j,k,bi,bj)=Zoograzave(i,j,k,bi,bj)+
c    &                      (gampn1*graz1*zo +gampn2*graz2*zo)*
c    &                       deltaTclock/float(nsubtime)
#ifdef GEIDER
             Chlave(i,j,k,bi,bj)=Chlave(i,j,k,bi,bj)+
     &                           Chl*dtplankton
#endif
             PARave(i,j,k,bi,bj)=PARave(i,j,k,bi,bj)+
     &                           PARl*dtplankton
             PPave(i,j,k,bi,bj)=PPave(i,j,k,bi,bj)+
     &                           PP*dtplankton
             Nfixave(i,j,k,bi,bj)=Nfixave(i,j,k,bi,bj)+
     &                           Nfix*dtplankton
             Denitave(i,j,k,bi,bj)=Denitave(i,j,k,bi,bj)+
     &                           denit*dtplankton
#ifdef WAVES_DIAG_PCHL
             do np=1,npmax
               Pchlave(i,j,k,bi,bj,np)=Pchlave(i,j,k,bi,bj,np)+
     &                           phychl(np)*dtplankton
             enddo
#endif
#ifdef DAR_DIAG_ACDOM
c            print*,'acdom',k,acdom_k(k,darwin_diag_acdom_ilam)
             aCDOMave(i,j,k,bi,bj)=aCDOMave(i,j,k,bi,bj)+
     &                    acdom_k(k,darwin_diag_acdom_ilam)*dtplankton
#endif
#ifdef DAR_DIAG_IRR
             do ilam = 1,tlam
              if (k.EQ.1) then
               Edave(i,j,k,bi,bj,ilam)=Edave(i,j,k,bi,bj,ilam)+
     &                                 Edwsf(ilam)*dtplankton
               Esave(i,j,k,bi,bj,ilam)=Esave(i,j,k,bi,bj,ilam)+
     &                                 Eswsf(ilam)*dtplankton
Coj            no Eu at surface (yet)
              else
               Edave(i,j,k,bi,bj,ilam)=Edave(i,j,k,bi,bj,ilam)+
     &                                 Edz(ilam,k-1)*dtplankton
               Esave(i,j,k,bi,bj,ilam)=Esave(i,j,k,bi,bj,ilam)+
     &                                 Esz(ilam,k-1)*dtplankton
               Euave(i,j,k,bi,bj,ilam)=Euave(i,j,k,bi,bj,ilam)+
     &                                 Euz(ilam,k-1)*dtplankton
              endif
              Eutave(i,j,k,bi,bj,ilam)=Eutave(i,j,k,bi,bj,ilam)+
     &                                 Eutop(ilam,k)*dtplankton
             enddo
#endif
#ifdef DAR_DIAG_ABSORP
             do ilam = 1,tlam
               aave(i,j,k,bi,bj,ilam)=aave(i,j,k,bi,bj,ilam)+
     &                                a_k(k,ilam)*dtplankton
             enddo
#endif
#ifdef DAR_DIAG_SCATTER
             do ilam = 1,tlam
               btave(i,j,k,bi,bj,ilam)=btave(i,j,k,bi,bj,ilam)+
     &                                 bt_k(k,ilam)*dtplankton
               bbave(i,j,k,bi,bj,ilam)=bbave(i,j,k,bi,bj,ilam)+
     &                                 bb_k(k,ilam)*dtplankton
             enddo
#endif
#ifdef DAR_DIAG_PART_SCATTER
             do ilam = 1,tlam
               apartave(i,j,k,bi,bj,ilam)=apartave(i,j,k,bi,bj,ilam)+
     &                                apart_k(k,ilam)*dtplankton
               btpartave(i,j,k,bi,bj,ilam)=btpartave(i,j,k,bi,bj,ilam)+
     &                                 bpart_k(k,ilam)*dtplankton
               bbpartave(i,j,k,bi,bj,ilam)=bbpartave(i,j,k,bi,bj,ilam)+
     &                                 bbpart_k(k,ilam)*dtplankton
             enddo
#endif
#ifdef DAR_DIAG_RSTAR
             do np=1,npmax
               Rstarave(i,j,k,bi,bj,np)=Rstarave(i,j,k,bi,bj,np)+
     &                           Rstarl(np)*dtplankton
               RNstarave(i,j,k,bi,bj,np)=RNstarave(i,j,k,bi,bj,np)+
     &                           RNstarl(np)*dtplankton
             enddo
#endif
#ifdef DAR_DIAG_DIVER
             Diver1ave(i,j,k,bi,bj)=Diver1ave(i,j,k,bi,bj)+
     &                           Diver1(i,j,k)*dtplankton
             Diver2ave(i,j,k,bi,bj)=Diver2ave(i,j,k,bi,bj)+
     &                           Diver2(i,j,k)*dtplankton
             Diver3ave(i,j,k,bi,bj)=Diver3ave(i,j,k,bi,bj)+
     &                           Diver3(i,j,k)*dtplankton
             Diver4ave(i,j,k,bi,bj)=Diver4ave(i,j,k,bi,bj)+
     &                           Diver4(i,j,k)*dtplankton
#endif
#ifdef DAR_DIAG_GROW
             do np=1,npmax
               Growave(i,j,k,bi,bj,np)=Growave(i,j,k,bi,bj,np)+
     &                           Growl(np)*dtplankton
               Growsqave(i,j,k,bi,bj,np)=Growsqave(i,j,k,bi,bj,np)+
     &                           Growsql(np)*dtplankton
             enddo
#endif

#ifdef ALLOW_DIAZ
#ifdef DAR_DIAG_NFIXP
             do np=1,npmax
               NfixPave(i,j,k,bi,bj,np)=NfixPave(i,j,k,bi,bj,np)+
     &                           NfixPl(np)*dtplankton
             enddo
#endif
#endif
#endif

#ifdef ALLOW_CARBON
             if (k.eq.1) then
              SURave(i,j,bi,bj)    =SURave(i,j,bi,bj)+
     &                              flxCO2(i,j)*dtplankton
              SURCave(i,j,bi,bj)    =SURCave(i,j,bi,bj)+
     &                              FluxCO2(i,j,bi,bj)*dtplankton
              SUROave(i,j,bi,bj)   =SUROave(i,j,bi,bj)+
     &                              flxO2(i,j)*dtplankton
              pCO2ave(i,j,bi,bj)   =pCO2ave(i,j,bi,bj)+
     &                              pCO2(i,j,bi,bj)*dtplankton
              pHave(i,j,bi,bj)     =pHave(i,j,bi,bj)+
     &                              pH(i,j,bi,bj)*dtplankton
             endif
#endif
          endif  
c end if hFac>0

         enddo ! k
c end layer loop
c

       ENDDO ! i
       ENDDO ! j

#ifdef ALLOW_PAR_DAY
C      1 <-> 2
       PARiaccum  = 3 - PARiprev

       DO k=1,nR
        DO j=1,sNy
         DO i=1,sNx
           PARday(i,j,k,bi,bj,PARiaccum) =
     &     PARday(i,j,k,bi,bj,PARiaccum) + PAR(i,j,k)
         ENDDO
        ENDDO
       ENDDO

       phase = 0. _d 0
       itistime = DIFF_PHASE_MULTIPLE( phase, darwin_PARavPeriod,
     &              newtime, dtsubtime)

       IF ( itistime ) THEN
C       compute average
        nav = darwin_PARnav
        IF (newtime - baseTime .LT. darwin_PARavPeriod) THEN
C         incomplete period at beginning of run
          nav = NINT((newtime-baseTime)/dtsubtime)
        ENDIF
        DO k=1,nR
         DO j=1,sNy
          DO i=1,sNx
            PARday(i,j,k,bi,bj,PARiaccum) =
     &      PARday(i,j,k,bi,bj,PARiaccum) / nav
          ENDDO
         ENDDO
        ENDDO
C       reset the other slot for averaging
        DO k=1,nR
         DO j=1,sNy
          DO i=1,sNx
            PARday(i,j,k,bi,bj,PARiprev) = 0. _d 0
          ENDDO
         ENDDO
        ENDDO
       ENDIF
C itistime
#endif

COJ fill diagnostics
#ifdef ALLOW_DIAGNOSTICS
       IF ( useDiagnostics ) THEN
        diagname = '        '
        WRITE(diagname,'(A8)') 'PAR     '
        CALL DIAGNOSTICS_FILL( PAR(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
        WRITE(diagname,'(A8)') 'PP      '
        CALL DIAGNOSTICS_FILL( PParr(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
        WRITE(diagname,'(A8)') 'Nfix    '
        CALL DIAGNOSTICS_FILL( Nfixarr(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
#ifdef DAR_DIAG_IRON
        WRITE(diagname,'(A8)') 'FeTbgc  '
        CALL DIAGNOSTICS_FILL( FeTbgcSource(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
C        WRITE(diagname,'(A8)') 'UDIAG4  '
C        CALL DIAGNOSTICS_FILL( consumpFeT(1-Olx,1-Oly,1), diagname,
C     &                         0,Nr,2,bi,bj,myThid )
C        WRITE(diagname,'(A8)') 'UDIAG5  '
C        CALL DIAGNOSTICS_FILL( reminFeT(1-Olx,1-Oly,1), diagname,
C     &                         0,Nr,2,bi,bj,myThid )
C        WRITE(diagname,'(A8)') 'UDIAG6  '
C        CALL DIAGNOSTICS_FILL( scavFeT(1-Olx,1-Oly,1), diagname,
C     &                         0,Nr,2,bi,bj,myThid )
#endif
c ANNA_TAVE
#ifdef WAVES_DIAG_PCHL
        DO np=1,MIN(99,npmax)
         WRITE(diagname,'(A5,I2.2,A1)') 'Pchl',np,' '
         CALL DIAGNOSTICS_FILL( Pchlarr(1-Olx,1-Oly,1,np), diagname,
     &                          0,Nr,2,bi,bj,myThid )
        ENDDO
#endif
c ANNA end TAVE
#ifdef DAR_DIAG_RSTAR
        DO np=1,MIN(99,npmax)
         WRITE(diagname,'(A5,I2.2,A1)') 'Rstar',np,' '
         CALL DIAGNOSTICS_FILL( Rstararr(1-Olx,1-Oly,1,np), diagname,
     &                          0,Nr,2,bi,bj,myThid )
        ENDDO
#endif
#ifdef DAR_DIAG_DIVER
        WRITE(diagname,'(A8)') 'Diver1  '
        CALL DIAGNOSTICS_FILL( Diver1(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
        WRITE(diagname,'(A8)') 'Diver2  '
        CALL DIAGNOSTICS_FILL( Diver2(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
        WRITE(diagname,'(A8)') 'Diver3  '
        CALL DIAGNOSTICS_FILL( Diver3(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
        WRITE(diagname,'(A8)') 'Diver4  '
        CALL DIAGNOSTICS_FILL( Diver4(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
        WRITE(diagname,'(A8)') 'Shannon '
        CALL DIAGNOSTICS_FILL( Shannon(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
        WRITE(diagname,'(A8)') 'Simpson '
        CALL DIAGNOSTICS_FILL( Simpson(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
#endif
#ifdef ALLOW_DIAZ
#ifdef DAR_DIAG_NFIXP
        DO np=1,MIN(99,npmax)
         WRITE(diagname,'(A5,I2.2,A1)') 'NfixP',np,' '
         CALL DIAGNOSTICS_FILL( NfixParr(1-Olx,1-Oly,1,np), diagname,
     &                          0,Nr,2,bi,bj,myThid )
        ENDDO
#endif
#endif
#ifdef DAR_DIAG_CHL
        CALL DIAGNOSTICS_FILL( GeiderChlarr(1-Olx,1-Oly,1), 'ChlGeide',
     &                         0,Nr,2,bi,bj,myThid )
        CALL DIAGNOSTICS_FILL( GeiderChl2Carr(1-Olx,1-Oly,1),'Chl2CGei',
     &                         0,Nr,2,bi,bj,myThid )
        CALL DIAGNOSTICS_FILL( DoneyChlarr(1-Olx,1-Oly,1), 'ChlDoney',
     &                         0,Nr,2,bi,bj,myThid )
        CALL DIAGNOSTICS_FILL( DoneyChl2Carr(1-Olx,1-Oly,1), 'Chl2CDon',
     &                         0,Nr,2,bi,bj,myThid )
        CALL DIAGNOSTICS_FILL( CloernChlarr(1-Olx,1-Oly,1), 'ChlCloer',
     &                         0,Nr,2,bi,bj,myThid )
        CALL DIAGNOSTICS_FILL( CloernChl2Carr(1-Olx,1-Oly,1),'Chl2CClo',
     &                         0,Nr,2,bi,bj,myThid )
#endif
#ifdef ALLOW_CARBON
        CALL DIAGNOSTICS_FILL( flxCO2(1-Olx,1-Oly), 'DICTFLX ',
     &                         0,1,2,bi,bj,myThid )
        CALL DIAGNOSTICS_FILL( FluxCO2(1-Olx,1-Oly,bi,bj), 'DICCFLX ',
     &                         0,1,2,bi,bj,myThid )
        CALL DIAGNOSTICS_FILL( flxO2(1-Olx,1-Oly), 'DICOFLX ',
     &                         0,1,2,bi,bj,myThid )
        CALL DIAGNOSTICS_FILL( pCO2(1-Olx,1-Oly,bi,bj), 'DICPCO2 ',
     &                         0,1,2,bi,bj,myThid )
        CALL DIAGNOSTICS_FILL( pH(1-Olx,1-Oly,bi,bj), 'DICPHAV ',
     &                         0,1,2,bi,bj,myThid )
#endif /* ALLOW_CARBON */
       ENDIF
#endif /* ALLOW_DIAGNOSTICS */
COJ

c determine iron partitioning  - solve for free iron
         call darwin_fe_chem(bi,bj,iMin,iMax,jMin,jMax,
     &                       Ptr(1-OLx,1-OLy,1,bi,bj,iFeT), freefe,
     &                       myIter, mythid)
c
#ifdef ALLOW_TIMEAVE
c save averages
       do k=1,nR
         dar_timeave(bi,bj,k)=dar_timeave(bi,bj,k)
     &                         +dtplankton
#ifdef ALLOW_CARBON
         dic_timeave(bi,bj,k)=dic_timeave(bi,bj,k)
     &                         +dtplankton
#endif
       enddo
#endif
c
c -----------------------------------------------------
       ENDDO ! it
c -----------------------------------------------------
c end of bio-chemical time loop
c 
        RETURN
        END
#endif  /*MONOD*/
#endif  /*ALLOW_PTRACERS*/

C============================================================================

