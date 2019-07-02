C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/monod/monod_plankton.F,v 1.6 2012/05/31 21:08:25 stephd Exp $
C $Name:  $

#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "DARWIN_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_MONOD

c ====================================================================
c SUBROUTINE MONOD_PLANKTON
c 1. Local ecological interactions for models with many phytoplankton
c    "functional groups"
c 2. Timestep plankton and nutrients locally
c 3. Includes explicit DOM and POM
c 4. Remineralization of detritus also determined in routine
c 5. Sinking particles and phytoplankton
c 6. NOT in this routine: iron chemistry
c
c       Mick Follows, Scott Grant, Fall/Winter 2005
c       Stephanie Dutkiewicz Spring/Summer 2006
c
c - add extra diagnostics, including R* (#define DAR_DIAG_RSTAR) - Stephanie, Spring 2007
c - add check for conservation (#define CHECK_CONS) - Stephanie, Spring 2007
c - improve grazing (#undef OLD_GRAZING) - Stephanie, Spring 2007
c - add diazotrophy (#define ALLOW_DIAZ) - Stephanie, Spring 2007
c - add mutation code (#define ALLOW_MUTANTS) - Jason Bragg, Spring/Summer 2007
c - new nitrogen limiting scheme (#undef OLD_NSCHEME) - Jason Bragg, Summer 2007
c - fix bug in diazotroph code  - Stephanie, Fall 2007
c - add additional r* diagnostic for (no3+no2) - Stephanie, Winter 2007
c - add diversity diagnostics - Stephanie, Winter 2007
c - add geider chl:c ratio and growth rate dependence,
c   though has no photo-inhibtion at this point - Stephanie, Spring 2008
c - add waveband dependence of light attenuation and absorption,
c   NOTE: need to have geider turned on too  - Anna Hickman, Summer 2008
c ====================================================================

c ANNA pass extra variables if WAVEBANDS
         SUBROUTINE MONOD_PLANKTON(
     U                       phyto,
     I                       zooP, zooN, zooFe, zooSi,
     O                       PP, Chl, Nfix, denit,
     I                       PO4local, NO3local, FeTlocal, Silocal,
     I                       NO2local, NH4local, 
     I                       DOPlocal, DONlocal, DOFelocal,
     I                       POPlocal, PONlocal, POFelocal, PSilocal,
     I                       phytoup, popuplocal, ponuplocal, 
     I                       pofeuplocal, psiuplocal,
     I                       PARlocal,Tlocal, Slocal,
     I                       freefelocal, inputFelocal, inputDustlocal,
     I                       bottom, dzlocal,
     I                       ztoplocal, zbotlocal,
     I                       ddusthindt, ddustsindt,
     I                       dPOCqahindt, dPOCqasindt,
     I                       dPOCqadt,
     I                       lukevar,
     O                       Rstarlocal, RNstarlocal,
#ifdef DAR_DIAG_GROW
     O                       Growlocal, Growsqlocal,
#endif
#ifdef ALLOW_DIAZ
#ifdef DAR_DIAG_NFIXP
     O                       NfixPlocal,
#endif
#endif
     O                       dphytodt, dzooPdt, dzooNdt, dzooFedt,
     O                       dzooSidt,
     O                       dPO4dt, dNO3dt, dFeTdt, dSidt,
     O                       dNH4dt, dNO2dt,
     O                       dDOPdt, dDONdt, dDOFedt, 
     O                       dPOPdt, dPONdt, dPOFedt, dPSidt,
#ifdef ALLOW_CARBON
     I                       DIClocal, DOClocal, POClocal,PIClocal,
     I                       POCnflocal, POCflocal, POClocalprod,
     I                       ALKlocal, O2local, ZooClocal,
     I                       POCuplocal, PICuplocal, POCuplocalprod,
     I                       pSilocalprod, PIClocalprod,
     I                       pSiuplocalprod, PICuplocalprod,
     I                       POCnfuplocal,
     I                       POCfuplocal,
     I                       dDUSTsfluxouplocal,
     I                       dDUSThfluxouplocal,
     I                       PSihlocal, PSihuplocal,
     I                       PIChlocal, PIChuplocal,
     I                       POChplocal, POChpuplocal,
     I                       POChslocal, POChsuplocal,
     I                       POCsplocal, POCspuplocal,
     I                       POCsslocal, POCssuplocal,
     I                       POChdlocal, POChduplocal,
     I                       POCsdlocal, POCsduplocal,
     I                       dustslocal, dustsuplocal,
     I                       dusthlocal, dusthuplocal,
     O                       dDICdt, dDOCdt, dPOCdt, dPICdt,
     O                       dALKdt, dO2dt, dZOOCdt,
     O                       dPOCnonfreedt, phy1_poc, phy2_poc,
     O                       dPOCfreedt, dPSiO3dt, dPCaCO3dt,
     O                       dDUSTsfluxodt, dDUSThfluxodt, dDUSTfluxodt,
     O                       dPOCproddt, dPICproddt, dpSiproddt,
     O                       dPOCnfproddt, dPOCfproddt,
     O                       dPOCavaildt,
     O                       dpSihdt, dPIChdt,
     O                       dPOChpdt, dPOChsdt,
     O                       dPOCspdt, dPOCssdt,
     O                       dPOCqastmpdt, dPOCqahtmpdt, dPOCqatmpdt,
     O                       dPOCqasdt, dPOCqahdt,
     O                       dPOChddt, dPOCsddt,
     O                       dDUSThdt, dDUSTsdt, 
#endif
#ifdef GEIDER
     I                       phychl, 
#ifdef DYNAMIC_CHL
     O                       dphychl, Chlup,
#endif
#ifdef ALLOW_CDOM
     O                       dcdomdt   , 
     I                       cdomlocal, 
#endif
#ifdef WAVEBANDS
     I                       PARwlocal,
#endif
#endif
#ifdef ALLOW_PAR_DAY
     I                       PARdaylocal,
#endif
#ifdef DAR_DIAG_CHL
     O                       ChlGeiderlocal, ChlDoneylocal,
     O                       ChlCloernlocal,
#endif
     I                       debug,
     I                       runtim,
     I                       MyThid) 
                  

         implicit none
#include "MONOD_SIZE.h"
#include "MONOD.h"
#include "DARWIN_PARAMS.h"

c ANNA set wavebands params
#ifdef WAVEBANDS
#include "SPECTRAL_SIZE.h"
#include "WAVEBANDS_PARAMS.h"
#endif


C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid
CEOP
c === GLOBAL VARIABLES =====================
c npmax = no of phyto functional groups
c nzmax = no of grazer species
c phyto = phytoplankton
c zoo   = zooplankton
         _RL phyto(npmax)
         _RL zooP(nzmax)
         _RL zooN(nzmax)
         _RL zooFe(nzmax)
         _RL zooSi(nzmax)
         _RL PP
         _RL Nfix
         _RL denit
         _RL Chl
         _RL PO4local
         _RL NO3local
         _RL FeTlocal
         _RL Silocal
         _RL NO2local
         _RL NH4local
         _RL DOPlocal
         _RL DONlocal
         _RL DOFelocal
         _RL POPlocal 
         _RL PONlocal
         _RL POFelocal
         _RL PSilocal
         _RL phytoup(npmax)
         _RL POPuplocal
         _RL PONuplocal
         _RL POFeuplocal
         _RL PSiuplocal 
         _RL PARlocal
         _RL Tlocal
         _RL Slocal
         _RL freefelocal
         _RL inputFelocal
         _RL inputDustlocal
         _RL bottom
         _RL dzlocal
         _RL ztoplocal
         _RL zbotlocal
         _RL Rstarlocal(npmax)
         _RL RNstarlocal(npmax)
#ifdef DAR_DIAG_GROW
         _RL Growlocal(npmax)
         _RL Growsqlocal(npmax)
#endif
#ifdef ALLOW_DIAZ
#ifdef DAR_DIAG_NFIXP
         _RL NfixPlocal(npmax)
#endif
#endif
         INTEGER debug
         _RL dphytodt(npmax)
         _RL dzooPdt(nzmax)
         _RL dzooNdt(nzmax)
         _RL dzooFedt(nzmax)
         _RL dzooSidt(nzmax)
         _RL dPO4dt
         _RL dNO3dt
         _RL dNO2dt
         _RL dNH4dt
         _RL dFeTdt
         _RL dSidt
         _RL dDOPdt
         _RL dDONdt
         _RL dDOFedt
         _RL dPOPdt
         _RL dPONdt
         _RL dPOFedt
         _RL dPSidt
#ifdef ALLOW_CARBON
         _RL DIClocal
         _RL DOClocal
         _RL POClocal
         _RL POClocalprod
         _RL PIClocal
         _RL POCnflocal
         _RL POCflocal
         _RL ALKlocal
         _RL O2local
         _RL ZooClocal(nzmax)
         _RL POCuplocal
         _RL POCuplocalprod
         _RL PICuplocal
         _RL dDICdt
         _RL dDOCdt
         _RL dPOCdt
         _RL dPOCproddt 
         _RL dPICproddt 
         _RL dpSiproddt
         _RL dPOCfreedt
         _RL dPOCnonfreedt
         _RL phy1_poc
         _RL phy2_poc
         _RL dPSiO3dt
         _RL dPCaCO3dt
         _RL dDUSTsfluxodt
         _RL dDUSThfluxodt
         _RL dDUSTfluxodt
         _RL dPICdt
         _RL dALKdt
         _RL dO2dt
         _RL dZOOCdt(nzmax)
         _RL pSilocalprod
         _RL PIClocalprod
         _RL pSiuplocalprod
         _RL PICuplocalprod
         _RL dPOCnfproddt
         _RL dPOCfproddt
         _RL dPOCavaildt
         _RL POCnfuplocal
         _RL POCfuplocal
         _RL dDUSTsfluxouplocal
         _RL dDUSThfluxouplocal
         _RL dpSihdt 
         _RL dPIChdt
         _RL dPOChpdt 
         _RL dPOChsdt
         _RL dPOCspdt
         _RL dPOCssdt
         _RL dPOCqastmpdt
         _RL dPOCqahtmpdt
         _RL dPOCqatmpdt
         _RL dPOCqasdt
         _RL dPOCqahdt
         _RL dPOChddt
         _RL dPOCsddt
         _RL dDUSThdt
         _RL dDUSTsdt
         _RL PSihlocal 
         _RL PSihuplocal
         _RL PIChlocal 
         _RL PIChuplocal
         _RL POChplocal
         _RL POChpuplocal
         _RL POChslocal
         _RL POChsuplocal
         _RL POCsplocal
         _RL POCspuplocal
         _RL POCsslocal 
         _RL POCssuplocal
         _RL POChdlocal
         _RL POChduplocal
         _RL POCsdlocal
         _RL POCsduplocal
         _RL dustslocal 
         _RL dustsuplocal
         _RL dusthlocal
         _RL dusthuplocal
#endif
#ifdef GEIDER
         _RL phychl(npmax)
#ifdef DYNAMIC_CHL
         _RL dphychl(npmax)
         _RL Chlup(npmax)
#endif
#endif
#ifdef ALLOW_CDOM
         _RL cdomlocal
         _RL dcdomdt  
#ifdef ALLOW_CARBON
         _RL cdomclocal, dcdomcdt    
#endif
#endif
#ifdef ALLOW_PAR_DAY
         _RL PARdaylocal
#endif
#ifdef DAR_DIAG_CHL
         _RL ChlGeiderlocal, ChlDoneylocal, ChlCloernlocal
#endif
         _RL runtim

c ANNA Global variables for WAVEBANDS
c ANNA these variables are passed in/out of darwin_forcing.F
#ifdef WAVEBANDS
         _RL PARwlocal(tlam)      !PAR at midpoint of previous(in) and local(out) gridcell
#endif
c ANNA endif





c LOCAL VARIABLES
c -------------------------------------------------------------

c WORKING VARIABLES
c np = phytoplankton index
         integer np
c nz = zooplankton index
         integer nz

c variables for ballast
c luke gloege
         _RL ddusthindt
         _RL ddustsindt
         _RL dPOCqahindt
         _RL dPOCqasindt
         _RL dPOCqadt
         INTEGER lukevar
         _RL new_DUST_def

c ballastlocal
         _RL dPIChproddt
         _RL dPICsproddt
         _RL dpSihproddt
         _RL dpSisproddt   
c         _RL dPOCnfproddt
c         _RL dPOCavaildt
         _RL dPOChpproddt
         _RL dPOChsproddt
         _RL dPOCspproddt
         _RL dPOCssproddt

c variables for phytoplankton growth rate/nutrient limitation
c phytoplankton specific nutrient limitation term
         _RL limit(npmax)
c phytoplankton light limitation term
         _RL ilimit(npmax)
         _RL ngrow(npmax)
         _RL grow(npmax)
         _RL PspecificPO4(npmax)
         _RL phytoTempFunction(npmax)
         _RL dummy
         _RL Ndummy
         _RL Nsourcelimit(npmax)
         _RL Nlimit(npmax)
         _RL NO3limit(npmax)
         _RL NO2limit(npmax)
         _RL NH4limit(npmax)

c for check N limit scheme
         _RL Ndiff
         _RL NO3limcheck
         _RL NO2limcheck
         _RL Ndummy1
         LOGICAL check_nlim

#ifndef OLD_NSCHEME
c [jbmodif] some new N terms
         integer N2only
         integer noNOdadv
         integer NOreducost
         _RL NO2zoNH4
         _RL NOXzoNH4
#endif

c varible for mimumum phyto
         _RL phytomin(npmax)

#ifdef OLD_GRAZE
c  variables for zooplankton grazing rates
         _RL zooTempFunction(nzmax)
         _RL grazing_phyto(npmax)
         _RL grazingP(nzmax)
         _RL grazingN(nzmax)
         _RL grazingFe(nzmax)
         _RL grazingSi(nzmax)
#else
c  variables for zooplankton grazing rates
         _RL zooTempFunction(nzmax)
         _RL allphyto(nzmax)
         _RL denphyto(nzmax)
         _RL grazphy(npmax,nzmax)
         _RL sumgrazphy(npmax)
         _RL sumgrazzoo(nzmax)
         _RL sumgrazzooN(nzmax)
         _RL sumgrazzooFe(nzmax)
         _RL sumgrazzooSi(nzmax)
         _RL sumgrazloss(nzmax)
         _RL sumgrazlossN(nzmax)
         _RL sumgrazlossFe(nzmax)
         _RL sumgrazlossSi(nzmax)
#endif

#ifdef GEIDER
        _RL alpha_I(npmax)
        _RL pcarbon(npmax)
        _RL pcm(npmax)
        _RL chl2c(npmax)
#ifdef DYNAMIC_CHL
        _RL acclim(npmax)
        _RL psinkchl(npmax)
        _RL rhochl(npmax)
#endif
#endif

#ifdef ALLOW_CDOM
        _RL cdomp_degrd, cdomn_degrd, cdomfe_degrd
        _RL preminP_cdom, preminN_cdom, preminFe_cdom
#ifdef ALLOW_CARBON
        _RL cdomc_degrd
        _RL preminC_cdom
#endif
#endif

#ifdef DAR_DIAG_CHL
         _RL tmppcm
         _RL tmpchl2c
#endif
c variables for nutrient uptake
         _RL consumpPO4
         _RL consumpNO3
         _RL consumpNO2
         _RL consumpNH4
         _RL consumpFeT
         _RL consumpSi

c variables for reminerlaization of DOM and POM
         _RL reminTempFunction
         _RL DOPremin
         _RL DONremin
         _RL DOFeremin
         _RL preminP
         _RL preminN
         _RL preminFe
         _RL preminSi
         _RL preminSih

c for sinking matter
         _RL psinkP
         _RL psinkN
         _RL psinkFe
         _RL psinkSi
         _RL psinkSih
         _RL psinkphy(npmax)

#ifdef ALLOW_CARBON
         _RL consumpDIC
         _RL consumpDIC_PIC
         _RL preminC
         _RL preminChp
         _RL preminChs
         _RL preminCsp
         _RL preminCss
         _RL preminChd
         _RL preminCsd
         _RL preminDUSTh
         _RL preminDUSTs
         _RL DOCremin
         _RL totphy_doc
         _RL totzoo_doc
         _RL totphy_poc
         _RL totzoo_poc
         _RL totphy_pic
         _RL totzoo_pic
         _RL psinkC
         _RL psinkChp
         _RL psinkChs
         _RL psinkCsp
         _RL psinkCss
         _RL psinkChd
         _RL psinkCsd
         _RL psinkPIC
         _RL psinkPICh
         _RL disscPIC
         _RL disscPICh
         _RL psinkDUSTh
         _RL psinkDUSTs
#ifdef OLD_GRAZE
         _RL grazingC(nzmax)
#else
c  variables for zooplankton grazing rates
         _RL sumgrazzooC(nzmax)
         _RL sumgrazlossC(nzmax)
         _RL sumgrazlossPIC(nzmax)
#endif

#endif

c variables for conversions from phyto and zoo to DOM and POM
         _RL totphy_dop
         _RL totphy_pop
         _RL totphy_don
         _RL totphy_pon
         _RL totphy_dofe
         _RL totphy_pofe
         _RL totphy_dosi
         _RL totphy_posi

         _RL totzoo_dop
         _RL totzoo_pop
         _RL totzoo_don
         _RL totzoo_pon
         _RL totzoo_dofe
         _RL totzoo_pofe
         _RL totzoo_posi

#ifdef ALLOW_CDOM
         _RL totphy_cdomp
         _RL totphy_cdomn
         _RL totphy_cdomfe
         _RL totzoo_cdomp
         _RL totzoo_cdomn
         _RL totzoo_cdomfe
#ifdef ALLOW_CARBON
         _RL totphy_cdomc
         _RL totzoo_cdomc
#endif
#endif

         _RL NO2prod
         _RL NO3prod
  
         _RL facpz

         _RL kpar, kinh

         _RL tmpr,tmpz, tmpgrow, tmp1, tmp2

         integer ITEST

#ifdef PART_SCAV
        _RL scav_part
        _RL scav_poc
#endif


c ANNA local variables for WAVEBANDS
#ifdef WAVEBANDS 
         integer i,ilam
         integer nl

c ANNA for interpolation
         _RL cu_area
C         _RL waves_diff
C         _RL light_diff
C         _RL alphaI_diff
C         _RL squ_part
C         _RL tri_part
C         _RL seg_area

c ANNA inportant but local variables that can be fogotten
         _RL PARwdn(tlam)         !light at bottom of local gridcell
         _RL attenwl(tlam)        !attenuation (m-1)
         _RL sumaphy_nl(tlam)     !total phyto absorption at each wavelength
#endif
c ANNA endif

c ANNA - for inhib
         _RL Ek
         _RL EkoverE

c.................................................................

#ifdef ALLOW_MUTANTS
c -m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-
c mutation variables [jbmodif]
         INTEGER nsisone
         INTEGER nsistwo  
         INTEGER nsisthree
         INTEGER nsisfour
         INTEGER npro
         INTEGER taxind         
         _RL mutfor, mutback
         _RL grow1
         _RL grow2
         _RL grow3
         _RL grow4
#endif

         INTEGER numtax
         _RL oneyr,threeyr

#ifdef  ALLOW_MUTANTS
c compile time options -- could maybe be moved to
c run time and set in data.gchem???
c QQQQQQQ
c Initialize sister taxon mutation scheme 
c if numtax = 1, mutation is off
         numtax = 4
c number of plankton types to assign for 
c wild and mutants types
         npro = 60
#else
         numtax=1
#endif

         oneyr = 86400.0 _d 0*360.0 _d 0 
         threeyr = oneyr*3. _d 0

c end mutation variables [jbmodif]
c -m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-

#ifndef OLD_NSCHEME
c [jbmodif] init new N terms
c  if those not using NO3 has 
c  N limit with denominator with NO3 or not: 0=NO3 in denom; 1=NO2 only
         N2only = 1
c ??
         noNOdadv = 1
c energetic disadvantage of using NO2/No3: off=0, on=1
         NOreducost =0
#endif

#ifdef GEIDER
         do np=1,npmax
           pcarbon(np) = 0. _d 0
           pcm(np)=0. _d 0
           chl2c(np)=0. _d 0
#ifdef DYNAMIC_CHL
           acclim(np)=0. _d 0
           psinkChl(np)=0. _d 0
#endif
         enddo
#endif
         

c set sum totals to zero
         totphy_pop = 0. _d 0
         totphy_dop = 0. _d 0
         totphy_don = 0. _d 0
         totphy_pon = 0. _d 0
         totphy_dofe = 0. _d 0
         totphy_pofe = 0. _d 0
         totphy_posi = 0. _d 0

         totzoo_dop = 0. _d 0
         totzoo_pop = 0. _d 0
         totzoo_don = 0. _d 0
         totzoo_pon = 0. _d 0
         totzoo_dofe = 0. _d 0
         totzoo_pofe = 0. _d 0
         totzoo_posi = 0. _d 0

         consumpPO4 = 0.0 _d 0
         consumpNO3 = 0.0 _d 0
         consumpNO2 = 0.0 _d 0
         consumpNH4 = 0.0 _d 0
         consumpFeT = 0.0 _d 0
         consumpSi  = 0.0 _d 0

#ifdef ALLOW_CARBON
         totphy_doc = 0. _d 0
         totphy_poc = 0. _d 0
         totphy_pic = 0. _d 0
         totzoo_doc = 0. _d 0
         totzoo_poc = 0. _d 0
         totzoo_pic = 0. _d 0
         consumpDIC = 0.0 _d 0
         consumpDIC_PIC = 0.0 _d 0
#endif

c zeros for diagnostics
        PP=0. _d 0
        Nfix=0. _d 0
        denit=0. _d 0
        Chl=0. _d 0

c set up phtyoplankton array to be used for grazing and mortality
c set up other variable used more than once to zero
         do np = 1, npmax
           dummy = phyto(np)-phymin
           phytomin(np)=max(dummy,0. _d 0)
           NH4limit(np)=0. _d 0
           NO2limit(np)=0. _d 0
           NO3limit(np)=0. _d 0
#ifdef ALLOW_DIAZ
#ifdef DAR_DIAG_NFIXP
           NfixPlocal(np)=0. _d 0
#endif
#endif
         enddo


#ifdef ALLOW_MUTANTS
c SWD if parent population is zero (ie. negative) treat all mutants
c as zeros too
       if(runtim .gt. threeyr) then
       if(numtax .gt. 1)then
        do np=1,npro
          if(mod(np,numtax).eq. 1. _d 0)then
             nsisone   = np
             nsistwo   = np+1
             nsisthree = np+2
             nsisfour  = np+3

             if (phyto(nsisone).le.0. _d 0) then
                if (numtax.gt.1) phyto(nsistwo)=0. _d 0
                if (numtax.gt.2) phyto(nsisthree)=0. _d 0
                if (numtax.gt.3) phyto(nsisfour)=0. _d 0
             endif
           endif
        enddo
       endif
       endif
ccccccccccccccccccccccccccccccc
#endif


c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          call MONOD_TEMPFUNC(Tlocal,phytoTempFunction,
     &            zooTempFunction, reminTempFunction, myThid)
          if (debug.eq.1) print*,'phytoTempFunction',
     &                          phytoTempFunction, Tlocal
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c ******************** GROWTH OF PHYTO ****************************
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifndef GEIDER
c ANNA also if not wavebands
#ifndef WAVEBANDS
c Determine phytoplantkon light limitation: will affect growth rate
c using Platt-like equations with inhibition
          do np = 1, npmax
           if (PARlocal.gt.1. _d 0) then
             kpar=ksatPAR(np)/10. _d 0
             kinh=kinhib(np)/1000. _d 0
             ilimit(np)=(1.0 _d 0 - EXP(-PARlocal*kpar))
     &                  *(EXP(-PARlocal*kinh)) / 
     &       ( kpar/(kpar+kinh)*EXP(kinh/kpar*LOG(kinh/(kpar+kinh))) )
             ilimit(np)=min(ilimit(np),1. _d 0)
           else
             ilimit(np)=0. _d 0
           endif
          enddo
          if (debug.eq.1) print*,'ilimit',ilimit, PARlocal
#endif
#endif
c ANNA endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Determine phytoplankton nutrient limitation as mimimum of 
c P,N,Si,Fe. However N can be utilized in several forms, so
c also determine which is used
          do np=1, npmax
            limit(np) = 1.0 _d 0
c P limitation
            if (ksatPO4(np).gt.0. _d 0) then
              dummy = PO4local/(PO4local+ksatPO4(np))
              if (dummy .lt. limit(np)) limit(np) = dummy
            endif
c Fe limitation
            if (ksatFeT(np).gt.0. _d 0) then
              dummy = FeTlocal/(FeTlocal+ksatFeT(np))
              if (dummy .lt. limit(np))limit(np) = dummy
            endif
c Si limiation
            if (R_SiP(np) .ne. 0. _d 0.and.ksatSi(np).gt.0. _d 0) then
              dummy = Silocal/(Silocal+ksatSi(np))
              if (dummy .lt. limit(np))limit(np) = dummy
            endif

c N limitation [jbmodif]
c nsource: genetic preference for {1:NH4&NO2 2:NH4 3:ALL Sources}
c Nsourcelimit marker for which nsource will be consumed {1:NO3 2:NO2 3:NH4}
c (Note: very different to way 1-D model does this)
            if(diazotroph(np) .ne. 1.0 _d 0)then

c NH4, all nsource
               if (ksatNH4(np).gt.0. _d 0) then
                   NH4limit(np) = NH4local/(NH4local+ksatNH4(np))
               endif

#ifdef OLD_NSCHEME
               if (ksatNO2(np).gt.0. _d 0) then
c NO2, if nsource is 1 or 3
                   NO2limit(np) = NO2local/(NO2local+ksatNO2(np))*
     &                    EXP(-sig1*NH4local)
                   NO2limcheck = NO2local/(NO2local+ksatNO2(np))
               endif
c NO3, if nsource is 3
               if (ksatNO3(np).gt.0. _d 0) then
                   NO3limit(np) = NO3local/(NO3local+ksatNO3(np))*
     &                    EXP(-sig2*NH4local - sig3*NO2local)
                   NO3limcheck = NO3local/(NO3local+ksatNO3(np))
               endif
#else
c [jbmodif]
c NO2, if nsource is 1 or 3
               if (ksatNO2(np).gt.0. _d 0 .and. nsource(np).ne.2) then
                 if (N2only.eq.1 .and. nsource(np).eq.1) then
c                if (nsource(np).eq.1) then
                  NO2limit(np) = NO2local/(NO2local+ksatNO2(np))
     &                    *EXP(-sig1*NH4local)
                  NO2limcheck = NO2local/(NO2local+ksatNO2(np))
                 else
                  if (ksatNO3(np).gt.0. _d 0) then
                   NO2limit(np)=NO2local/(NO3local+NO2local+ksatNO3(np))
     &                    *EXP(-sig1*NH4local)
                   NO2limcheck=NO2local/(NO3local+NO2local+ksatNO3(np))
                  endif
                 endif
               endif
c NO3, if nsource is 3
               if (ksatNO3(np).gt.0. _d 0 .and. nsource(np).eq.3) then
                  NO3limit(np)=NO3local/(NO3local+NO2local+ksatNO3(np))
     &                    *EXP(-sig1*NH4local)
                  NO3limcheck=NO3local/(NO3local+NO2local+ksatNO3(np))
               endif

#endif

               if (nsource(np).eq.2) then
                  NO2limit(np) = 0. _d 0
                  NO3limit(np) = 0. _d 0
                  NO2limcheck  = 0. _d 0
                  NO3limcheck  = 0. _d 0
               endif
               if (nsource(np).eq.1) then
                 NO3limit(np) = 0. _d 0
                 NO3limcheck  = 0. _d 0
               endif
               if (nsource(np).eq.3) then
c don't do anything
               endif

               Ndummy = NO3limit(np)+NO2limit(np)+NH4limit(np)
c
c make sure no Nlim disadvantage;
c check that limit doesn't decrease at high NH4 levels
               check_nlim=.FALSE.
               if (check_nlim) then
                  Ndummy1=NO3limcheck+NO2limcheck
                  if (Ndummy.gt.0. _d 0.and.Ndummy.lt.Ndummy1) then
c                   print*,'QQ N limit WARNING',Ndummy, Ndummy1,
c    &                             NO3local,NO2local,NH4local
                    Ndiff=Ndummy1-NH4limit(np)
                    NO2limit(np)=Ndiff * 
     &                     NO2limit(np)/(NO2limit(np)+NO3limit(np))
                    NO3limit(np)=Ndiff *
     &                     NO3limit(np)/(NO2limit(np)+NO3limit(np))
                    Ndummy = NO3limit(np)+NO2limit(np)+NH4limit(np)
                  endif
               endif

               if (Ndummy.gt.1. _d 0) then
                 NO3limit(np) = NO3limit(np)/Ndummy
                 NO2limit(np) = NO2limit(np)/Ndummy
                 NH4limit(np) = NH4limit(np)/Ndummy 
               endif
               Nlimit(np)=NO3limit(np)+NO2limit(np)+NH4limit(np)
               if (Nlimit(np).gt.1.01 _d 0) then
                print*,'QQ Nlimit', Nlimit(np), NO3limit(np), 
     &                              NO2limit(np), NH4limit(np)
               endif
               if (Nlimit(np).le.0. _d 0) then
c                if (np.eq.1) then
c                 print*,'QQ Nlimit', Nlimit(np), NO3limit(np),
c    &                              NO2limit(np), NH4limit(np)
c                 print*,'QQ limit',limit(np), np
c                endif
                 Nlimit(np)=0. _d 0   !1 _d -10
               endif

#ifdef OLD_NSCHEME
c lower growth for higher NO3 consumption at higher light
               if (Nlimit(np).le.0. _d 0) then
                 ngrow(np)=1. _d 0
               else
                if (parlocal.gt.ilight) then
                  ngrow(np)=ngrowfac+(1. _d 0-ngrowfac)*
     &                      (NH4limit(np)+NO2limit(np))/Nlimit(np)
                else
                  ngrow(np)=1. _d 0
                endif
                ngrow(np)=min(ngrow(np),1. _d 0)
               endif
#else
c disadvantage of oxidized inorganic N
c for now, ignore - a first attempt is included below
              ngrow(np) = 1.0 _d 0

cc lower growth for higher NO3 consumption at higher light
c one possible way of counting cost of reducing NOX
               if (NOreducost .eq. 1)then 
                if (Nlimit(np).le.0. _d 0) then
                  ngrow(np)=1. _d 0
                else
                 ngrow(np)= (10. _d 0*4. _d 0 +2. _d 0)   /
     &              (10. _d 0*4. _d 0 +2. _d 0*NH4limit(np)/Nlimit(np)
     &                           +8. _d 0*NO2limit(np)/Nlimit(np)
     &                           +10. _d 0*NO3limit(np)/Nlimit(np))
                 ngrow(np)=min(ngrow(np),1. _d 0)
                endif
               endif
c
c might consider other costs, too
c             if (NOironcost .eq. 1)then 
c
c             endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#endif

c Now Check Against General Nutrient Limiting Tendency
               if (ksatNH4(np).gt.0. _d 0.or.ksatNO2(np).gt.0. _d 0
     &              .or.ksatNO3(np).gt.0. _d 0) then
                 if(Nlimit(np) .lt. limit(np)) limit(np) = Nlimit(np)
               endif
             else
                ngrow(np)=1. _d 0
                Nlimit(np)=1. _d 0
                NO3limit(np)=0. _d 0
                NO2limit(np)=0. _d 0
                NH4limit(np)=0. _d 0
             endif   ! diaz
             limit(np)=min(limit(np),1. _d 0)
           enddo    !np
          if (debug.eq.1) print*,'nut limit',
     &                          limit, PO4local, FeTlocal, Silocal
          if (debug.eq.1) print*,'Nlimit',
     &                 Nlimit
          if (debug.eq.1) print*,'NH4limit',
     &                  NH4limit, NH4local
          if (debug.eq.1) print*,'NO2limit',
     &                   NO2limit, NO2local
          if (debug.eq.1) print*,'NO3limit',
     &                  NO3limit, NO3local
          if (debug.eq.1) print*,'ngrow',
     &                  ngrow


#ifdef GEIDER

#ifdef WAVEBANDS
c ANNA if wavebands then uses spectral alphachl derived from spectral alpha * I
c so first get value for alphachl_nl * PARwlocal
c value will depend on matchup between spectra of alphachl_nl (ie. aphy_chl) and PARwlocal
c integrate alpha*PAR over wavebands
          do np = 1,npmax
           alpha_I(np) = 0 _d 0
           do nl = 1,tlam
            alpha_I(np) = alpha_I(np) + alphachl_nl(np,nl)*PARwlocal(nl)
           end do
          end do
c Geider growth (and chl2c) now depends on this (sinlge) value of alpha_chl * I

c alpha_mean now precomputed in darwin_init_vari
#else
c ANNA if not wavebands uses alphachl derived from mQyield * aphy_chl_ave
c for use with generic geider equation need to use alpha_I (ie. alphachl*PARlocal)
          do np = 1, npmax
            alpha_I(np)=alphachl(np)*PARlocal
          enddo
c ANNA endif
#endif

          do np = 1, npmax
           pcm(np)=pcmax(np)*limit(np)*phytoTempFunction(np)
#ifdef DYNAMIC_CHL
           if (phyto(np).gt. 0. _d 0) then
             chl2c(np)=phychl(np)/(phyto(np)*R_PC(np))
           else
             chl2c(np)= 0. _d 0
           endif
#endif
           if (pcm(np).gt.0.d0) then
#ifndef DYNAMIC_CHL
c assumes balanced growth, eq A14 in Geider et al 1997
             chl2c(np)=chl2cmax(np)/
     &            (1+(chl2cmax(np)*alpha_I(np))/
     &               (2*pcm(np)))
             chl2c(np)=min(chl2c(np),chl2cmax(np))
             chl2c(np)=max(chl2c(np),chl2cmin(np))
#endif
             if (PARlocal.gt.1. _d -1) then
c Eq A1 in Geider et al 1997
               pcarbon(np)=pcm(np)*( 1 -
     &          exp((-alpha_I(np)*chl2c(np))/(pcm(np))) )
c for inhibition
               if (inhibcoef_geid(np).gt.0. _d 0) then
#ifdef WAVEBANDS 
                 Ek = pcm(np)/(chl2c(np)*alpha_mean(np))
#else
                 Ek = pcm(np)/(chl2c(np)*alphachl(np))
#endif
                 EkoverE = Ek / PARlocal
                 if (PARlocal .ge. Ek) then !photoinhibition begins
                  pcarbon(np) = pcarbon(np)*(EkoverE*inhibcoef_geid(np))
                 endif
               endif
c end inhib
               if (pcarbon(np).lt. 0. _d 0) 
     &                 print*,'QQ ERROR pc=',np,pcarbon(np)
               if (pcm(np).gt.0. _d 0) then
                 ilimit(np)=pcarbon(np)/pcm(np)
               else
                 ilimit(np)= 0. _d 0
               endif
             else
               ilimit(np)=0. _d 0
               pcarbon(np)=0. _d 0
             endif
           else ! if pcm 0
             pcm(np)=0.d0
#ifndef DYNAMIC_CHL
             chl2c(np)=chl2cmin(np)
#endif
             pcarbon(np)=0.d0
             ilimit(np)=0.d0
           endif
#ifdef DYNAMIC_CHL
c Chl:C acclimated to current conditions 
c (eq A14 in Geider et al 1997)
          acclim(np)=chl2cmax(np)/
     &            (1+(chl2cmax(np)*alpha_I(np))/
     &               (2*pcm(np)))
          acclim(np)=min(acclim(np),chl2cmax(np))
c         acclim(np)=max(acclim(np),chl2cmin(np))
#else
           phychl(np)=phyto(np)*R_PC(np)*chl2c(np)
#endif
          enddo
          if (debug.eq.14) print*,'ilimit',ilimit, PARlocal
          if (debug.eq.14) print*,'chl:c',chl2c
          if (debug.eq.14) print*,'chl',phychl
#ifdef DYNAMIC_CHL
          if (debug.eq.14) print*,'acclim',acclim
#endif
#endif /* GEIDER */

#ifdef DAR_DIAG_CHL
c diagnostic version of the above that does not feed back to growth
          ChlGeiderlocal = 0. _d 0
          do np = 1, npmax
           tmppcm = mu(np)*limit(np)*phytoTempFunction(np)
           if (tmppcm.gt.0.d0) then
             tmpchl2c = Geider_chl2cmax(np)/
     &         (1+(Geider_chl2cmax(np)*Geider_alphachl(np)*PARdaylocal)/
     &            (2*tmppcm))
             tmpchl2c = min(tmpchl2c, Geider_chl2cmax(np))
             tmpchl2c = max(tmpchl2c, Geider_chl2cmin(np))
           else
             tmpchl2c = Geider_chl2cmin(np)
           endif
           ChlGeiderlocal = ChlGeiderlocal + phyto(np)*R_PC(np)*tmpchl2c
          enddo
C Chl a la Doney
          ChlDoneylocal = 0. _d 0
          do np = 1, npmax
            tmpchl2c = (Doney_Bmax - (Doney_Bmax-Doney_Bmin)*
     &                          MIN(1. _d 0,PARdaylocal/Doney_PARstar))
     &                    *limit(np)
            ChlDoneylocal = ChlDoneylocal +
     &                      tmpchl2c*R_PC(np)*phyto(np)
          enddo
C Chl a la Cloern
          ChlCloernlocal = 0. _d 0
          do np = 1, npmax
            tmpchl2c = Cloern_chl2cmin +
     &                     Cloern_A*exp(Cloern_B*Tlocal)
     &                             *exp(-Cloern_C*PARdaylocal)
     &                    *limit(np)
            ChlCloernlocal = ChlCloernlocal +
     &                       tmpchl2c*R_PC(np)*phyto(np)
          enddo
#endif /* DAR_DIAG_CHL */


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ******************* END GROWTH PHYTO ******************************* 


#ifdef OLD_GRAZE
c------------------------------------------------------------------------
c GRAZING sum contributions of all zooplankton
           do np=1,npmax
              grazing_phyto(np) = 0.0 _d 0
              do nz = 1, nzmax
                   grazing_phyto(np) = grazing_phyto(np)
     &               + graze(np,nz)*zooP(nz)*zooTempFunction(nz)
              enddo
           enddo
           if (debug.eq.2) print*,'grazing_phyto',grazing_phyto
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#else
c------------------------------------------------------------------------
c sum all palatability*phyto and find phyto specific grazing rate
          do nz=1,nzmax
           allphyto(nz)=0. _d 0
           do np=1,npmax
            allphyto(nz)=allphyto(nz)+palat(np,nz)*phyto(np)
           enddo
           if (allphyto(nz).le.0. _d 0) allphyto(nz)=phygrazmin
#ifdef SER_GRAZ
           denphyto(nz)=0. _d 0
           do np=1,npmax
            denphyto(nz)=denphyto(nz)+
     &            (palat(np,nz)*phyto(np)/allphyto(nz))*phyto(np)
           enddo
           if (denphyto(nz).le.0. _d 0) denphyto(nz)=phygrazmin
#endif
           do np=1,npmax
            tmpz=max(0. _d 0,(allphyto(nz)-phygrazmin) )
            grazphy(np,nz)=grazemax(nz)*
#ifdef SER_GRAZ
c as in Vallina et al, 2011
     &          (((palat(np,nz)*phyto(np)/allphyto(nz))*phyto(np))/
     &           denphyto(nz)) *
#else
c as in Dutkiewicz et al, GBC, 2009
     &           (palat(np,nz)*phyto(np)/allphyto(nz))*
#endif
     &           ( tmpz**hollexp/
     &             (tmpz**hollexp+kgrazesat**hollexp) )
           enddo
          enddo
          if (debug.eq.2) print*,'allphyto',allphyto
c         if (debug.eq.2) print*,'grazephy',grazphy
c sum over zoo for impact on phyto
          do np=1,npmax
           sumgrazphy(np)=0. _d 0
           do nz=1,nzmax
             sumgrazphy(np)=sumgrazphy(np)+
     &           grazphy(np,nz)*zooP(nz)
           enddo
          enddo
          if (debug.eq.2) print*,'sumgrazephy',sumgrazphy
c sum over phy for impact on zoo, and all remainder to go to POM
          do nz=1,nzmax
           sumgrazzoo(nz)=0. _d 0
           sumgrazzooN(nz)=0. _d 0
           sumgrazzooFe(nz)=0. _d 0
           sumgrazzooSi(nz)=0. _d 0
           sumgrazloss(nz)=0. _d 0
           sumgrazlossN(nz)=0. _d 0
           sumgrazlossFe(nz)=0. _d 0
           sumgrazlossSi(nz)=0. _d 0
#ifdef ALLOW_CARBON
           sumgrazzooC(nz)=0. _d 0
           sumgrazlossC(nz)=0. _d 0
           sumgrazlossPIC(nz)=0. _d 0
#endif
           do np=1,npmax
              sumgrazzoo(nz)=sumgrazzoo(nz)+
     &             asseff(np,nz)*grazphy(np,nz)*zooP(nz)
              sumgrazloss(nz)=sumgrazloss(nz)+
     &             (1. _d 0-asseff(np,nz))*grazphy(np,nz)*zooP(nz)
              sumgrazzooN(nz)=sumgrazzooN(nz)+
     &             asseff(np,nz)*grazphy(np,nz)*zooP(nz)*R_NP(np)
              sumgrazlossN(nz)=sumgrazlossN(nz)+
     &             (1. _d 0-asseff(np,nz))*grazphy(np,nz)*
     &                                      zooP(nz)*R_NP(np)
              sumgrazzooFe(nz)=sumgrazzooFe(nz)+
     &             asseff(np,nz)*grazphy(np,nz)*
     &                                      zooP(nz)*R_FeP(np)
              sumgrazlossFe(nz)=sumgrazlossFe(nz)+
     &             (1. _d 0-asseff(np,nz))*grazphy(np,nz)*
     &                                      zooP(nz)*R_FeP(np)
              sumgrazzooSi(nz)=sumgrazzooSi(nz)+
     &             asseff(np,nz)*grazphy(np,nz)*
     &                                     zooP(nz)*R_SiP(np)
              sumgrazlossSi(nz)=sumgrazlossSi(nz)+
     &             (1. _d 0-asseff(np,nz))*grazphy(np,nz)*
     &                                     zooP(nz)*R_SiP(np)
#ifdef ALLOW_CARBON
              sumgrazzooC(nz)=sumgrazzooC(nz)+
     &             asseff(np,nz)*grazphy(np,nz)*zooP(nz)*R_PC(np)
              sumgrazlossC(nz)=sumgrazlossC(nz)+
     &             (1. _d 0-asseff(np,nz))*grazphy(np,nz)*
     &                                      zooP(nz)*R_PC(np)
              sumgrazlossPIC(nz)=sumgrazlossPIC(nz)+
     &             (1. _d 0)*grazphy(np,nz)*
     &                           zooP(nz)*R_PC(np)*R_PICPOC(np)
#endif
           enddo
          enddo
          if (debug.eq.2) print*,'sumgrazzoo',sumgrazzoo
          if (debug.eq.2) print*,'sumgrazloss',sumgrazloss
          if (debug.eq.2) print*,'sumgrazzooN',sumgrazzooN
          if (debug.eq.2) print*,'sumgrazlossN',sumgrazlossN
          if (debug.eq.2) print*,'sumgrazzooFe',sumgrazzooFe
          if (debug.eq.2) print*,'sumgrazlossFe',sumgrazlossFe
          if (debug.eq.2) print*,'sumgrazzooSi',sumgrazzooSi
          if (debug.eq.2) print*,'sumgrazlossSi',sumgrazlossSi
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c accumulate particulate and dissolved detritus
           do np=1, npmax
               totphy_pop=totphy_pop+
     &                ExportFracP(np)*mortphy(np)*phytomin(np)
               totphy_dop=totphy_dop+
     &                (1. _d 0-ExportFracP(np))*mortphy(np)*phytomin(np)
               totphy_pon=totphy_pon+ R_NP(np)*
     &                ExportFracP(np)*mortphy(np)*phytomin(np)
               totphy_don=totphy_don+ R_NP(np)*
     &               (1. _d 0-ExportFracP(np))*mortphy(np)*phytomin(np)
               totphy_pofe=totphy_pofe+ R_FeP(np)*
     &                ExportFracP(np)*mortphy(np)*phytomin(np)
               totphy_dofe=totphy_dofe+ R_FeP(np)*
     &               (1. _d 0-ExportFracP(np))*mortphy(np)*phytomin(np)
               totphy_posi=totphy_posi+ R_SiP(np)*
     &                mortphy(np)*phytomin(np)
#ifdef ALLOW_CARBON
C commented out by Luke want to separate POC into small and large 
C Becuase we are assuming a portion of large phyto is CaCO3
C               totphy_poc=totphy_poc+ R_PC(np)*
C     &                ExportFracP(np)*mortphy(np)*phytomin(np)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c POC and PIC for large (np=1) and small (np=2) species
c
c small : 7% POC routed to PIC remainder is PIC
c
c Large: all is routed to POC
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c large phytoplankton (diatoms)
              if (np.eq.1) then
c               large POC
                 phy1_poc=R_PC(np)*ExportFracP(np)*
     &                    mortphy(np)*phytomin(np)

c               large POC added to total pool
                  totphy_poc = totphy_poc + phy1_poc
              endif

c small phytoplankton
               if (np.eq.2) then
C               small POC
                  phy2_poc=(1-0.07)*R_PC(np)*ExportFracP(np)*
     &                     mortphy(np)*phytomin(np)

C               small POC added to total pool
                  totphy_poc = totphy_poc + phy2_poc

C               small PIC (route 7% to PIC pool)
                  totphy_pic=totphy_pic + 0.07*phy2_poc
               endif
c ====================================================================

               totphy_doc=totphy_doc+ R_PC(np)*
     &               (1. _d 0-ExportFracP(np))*mortphy(np)*phytomin(np)

C Commented out by Luke (this is the original code)
C               totphy_pic=totphy_pic+ R_PC(np)*R_PICPOC(np)*
C     &                mortphy(np)*phytomin(np)
#endif
            enddo 
#ifdef ALLOW_CDOM
               totphy_cdomp=(fraccdom)*totphy_dop
               totphy_dop=totphy_dop-totphy_cdomp
               totphy_cdomn=rnp_cdom*totphy_cdomp
               totphy_don=totphy_don-totphy_cdomn
               totphy_cdomfe=rfep_cdom*totphy_cdomp
               totphy_dofe=totphy_dofe-totphy_cdomfe
#ifdef ALLOW_CARBON
               totphy_cdomc=rcp_cdom*totphy_cdomp
               totphy_doc=totphy_doc-totphy_cdomc
#endif
#endif
            if (debug.eq.3) print*,'tot_phy_pop',totphy_pop
            if (debug.eq.3) print*,'tot_phy_dop',totphy_dop
            if (debug.eq.3) print*,'tot_phy_pon',totphy_pon
            if (debug.eq.3) print*,'tot_phy_don',totphy_don
            if (debug.eq.3) print*,'tot_phy_pofe',totphy_pofe
            if (debug.eq.3) print*,'tot_phy_dofe',totphy_dofe
            if (debug.eq.3) print*,'tot_phy_posi',totphy_posi


c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


#ifdef OLD_GRAZE
c ****************** ZOO GRAZING RATE ****************************
c determine zooplankton grazing rates 
            do nz = 1, nzmax
c grazing: sum contribution from all phytoplankton
              grazingP(nz) = 0.0 _d 0
              grazingN(nz) = 0.0 _d 0
              grazingFe(nz) = 0.0 _d 0
              grazingSi(nz) = 0.0 _d 0
#ifdef ALLOW_CARBON
              grazingC(nz) = 0.0 _d 0
#endif
              do np = 1, npmax
               facpz = (phytomin(np)/(phytomin(np) + kgrazesat))
     &                  *zooTempFunction(nz)
               grazingP(nz) = grazingP(nz) + 
     &                         graze(np,nz)*facpz
               grazingN(nz) = grazingN(nz) + 
     &                         graze(np,nz)*R_NP(np)*facpz
               grazingFe(nz) = grazingFe(nz) + 
     &                          graze(np,nz)*R_FeP(np)*facpz
               grazingSi(nz) = grazingSi(nz) + 
     &                          graze(np,nz)*R_SiP(np)*facpz
#ifdef ALLOW_CARBON
               grazingC(nz) = grazingC(nz) +
     &                         graze(np,nz)*R_PC(np)*facpz
#endif
              enddo
            enddo
            if (debug.eq.4) print*,'grazingP', grazingP
            if (debug.eq.4) print*,'grazingN', grazingN
            if (debug.eq.4) print*,'grazingFe', grazingFe
            if (debug.eq.4) print*,'grazingSi', grazingSi
c *************  END ZOO GRAZING *********************************
#endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c accumulate particulate and dissolved detritus
            do nz=1, nzmax
               totzoo_pop=totzoo_pop+
     &                     ExportFracZ(nz)*( mortzoo(nz)*zooP(nz)
     &                     +  mortzoo2(nz)*zooP(nz)**2 )
               totzoo_dop=totzoo_dop+
     &                    (1. _d 0-ExportFracZ(nz))*(
     &                     mortzoo(nz)*zooP(nz)+
     &                     mortzoo2(nz)*zooP(nz)**2 )
               totzoo_pon=totzoo_pon+
     &                     ExportFracZ(nz)*( mortzoo(nz)*zooN(nz)
     &                     +  mortzoo2(nz)*zooN(nz)**2 )
               totzoo_don=totzoo_don+
     &                    (1. _d 0-ExportFracZ(nz))*(
     &                     mortzoo(nz)*zooN(nz)+
     &                     mortzoo2(nz)*zooN(nz)**2 )
               totzoo_pofe=totzoo_pofe+
     &                     ExportFracZ(nz)*( mortzoo(nz)*zooFe(nz)
     &                     +  mortzoo2(nz)*zooFe(nz)**2 )
               totzoo_dofe=totzoo_dofe+
     &                   (1. _d 0-ExportFracZ(nz))*(
     &                    mortzoo(nz)*zooFe(nz) +
     &                    mortzoo2(nz)*zooFe(nz)**2 )
               totzoo_posi=totzoo_posi+
     &                     ( mortzoo(nz)*zooSi(nz)+
     &                       mortzoo2(nz)*zooSi(nz)**2 )
#ifdef ALLOW_CARBON
               totzoo_poc=totzoo_poc+
     &                    ExportFracZ(nz)*( mortzoo(nz)*zooClocal(nz)
     &                      +  mortzoo2(nz)*zooClocal(nz)**2 )
               totzoo_doc=totzoo_doc+
     &           (1. _d 0-ExportFracZ(nz))*( mortzoo(nz)*zooClocal(nz)
     &                      +  mortzoo2(nz)*zooClocal(nz)**2 )
#endif
            enddo

#ifndef OLD_GRAZE
            do nz=1, nzmax
               totzoo_pop=totzoo_pop+
     &                    ExportFracGraz(nz)*sumgrazloss(nz)
               totzoo_dop=totzoo_dop+
     &                 (1. _d 0-ExportFracGraz(nz))*sumgrazloss(nz)
               totzoo_pon=totzoo_pon+
     &                    ExportFracGraz(nz)*sumgrazlossN(nz)
               totzoo_don=totzoo_don+
     &                 (1. _d 0-ExportFracGraz(nz))*sumgrazlossN(nz)
               totzoo_pofe=totzoo_pofe+
     &                    ExportFracGraz(nz)*sumgrazlossFe(nz)
               totzoo_dofe=totzoo_dofe+
     &                 (1. _d 0-ExportFracGraz(nz))*sumgrazlossFe(nz)
               totzoo_posi=totzoo_posi+
     &                     sumgrazlossSi(nz)
#ifdef ALLOW_CARBON
               totzoo_poc=totzoo_poc+
     &                    ExportFracGraz(nz)*sumgrazlossC(nz)
               totzoo_doc=totzoo_doc+
     &                 (1. _d 0-ExportFracGraz(nz))*sumgrazlossC(nz)
               totzoo_pic=totzoo_pic+
     &                    sumgrazlossPIC(nz)
#endif
            enddo
#endif

#ifdef ALLOW_CDOM
               totzoo_cdomp=(fraccdom)*totzoo_dop
               totzoo_dop=totzoo_dop-totzoo_cdomp
               totzoo_cdomn=rnp_cdom*totzoo_cdomp
               totzoo_don=totzoo_don-totzoo_cdomn
               totzoo_cdomfe=rfep_cdom*totzoo_cdomp
               totzoo_dofe=totzoo_dofe-totzoo_cdomfe
#ifdef ALLOW_CARBON
               totzoo_cdomc=rcp_cdom*totzoo_cdomp
               totzoo_doc=totzoo_doc-totzoo_cdomc
#endif
#endif
            if (debug.eq.5) print*,'totzoo_pop',totzoo_pop
            if (debug.eq.5) print*,'totzoo_dop',totzoo_dop
            if (debug.eq.5) print*,'totzoo_pon',totzoo_pon
            if (debug.eq.5) print*,'totzoo_don',totzoo_don
            if (debug.eq.5) print*,'totzoo_pofe',totzoo_pofe
            if (debug.eq.5) print*,'totzoo_dofe',totzoo_dofe
            if (debug.eq.5) print*,'totzoo_posi',totzoo_posi
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c ********************* NUTRIENT UPTAKE *******************************
c determine nutrient uptake
c consumption - sum of phytoplankton contributions
            do np = 1, npmax
c phospate uptake by each phytoplankton
#ifndef GEIDER
               grow(np)=ngrow(np)*mu(np)*limit(np)*ilimit(np)*
     &                            phytoTempFunction(np)
#endif
#ifdef GEIDER
               grow(np)=ngrow(np)*pcarbon(np)
               if (debug.eq.1) print*,'grow', grow(np), pcarbon(np)
               if (debug.eq.14) print*,'grow', grow(np), pcarbon(np)
#ifdef DYNAMIC_CHL
c geider 97 for dChl/dt (source part) Eq. 3
               if (acclim(np).gt. 0. _d 0) then
                rhochl(np)=chl2cmax(np) *
     &                      (grow(np)/(alpha_I(np)*acclim(np)) )
               else
                rhochl(np)= 0. _d 0
               endif
               if (debug.eq.14) print*,'rhochl',rhochl(np)
#endif
#endif
               PspecificPO4(np) = grow(np)*phyto(np)
c        write(6,*)'np =',np, ' PspecificPO4 ='
c     &              ,PspecificPO4(np)
               consumpPO4 = consumpPO4 + PspecificPO4(np)
               consumpFeT = consumpFeT + PspecificPO4(np)*R_FeP(np)
               consumpSi  = consumpSi  + PspecificPO4(np)*R_SiP(np)
cswd should have O2prod as function of np?
c New Way of doing Nitrogen Consumption .......................
               if(diazotroph(np) .ne. 1.0 _d 0)then
                if (Nlimit(np).le.0. _d 0) then
                 consumpNO3 = consumpNO3
                 consumpNO2 = consumpNO2
                 consumpNH4 = consumpNH4
                else 
                 consumpNO3 = consumpNO3 + 
     &              NO3limit(np)/Nlimit(np)*PspecificPO4(np)*R_NP(np)
                 consumpNO2 = consumpNO2 +
     &              NO2limit(np)/Nlimit(np)* PspecificPO4(np)*R_NP(np)
                 consumpNH4 = consumpNH4 + 
     &              NH4limit(np)/Nlimit(np)*PspecificPO4(np)*R_NP(np)
                endif
               else
                 consumpNO3 = consumpNO3
                 consumpNO2 = consumpNO2
                 consumpNH4 = consumpNH4
                 Nfix=Nfix+PspecificPO4(np)*R_NP(np)
#ifdef ALLOW_DIAZ
#ifdef DAR_DIAG_NFIXP
                 NfixPlocal(np)=PspecificPO4(np)*R_NP(np)
#endif
#endif
               endif
#ifdef ALLOW_CARBON
               consumpDIC  = consumpDIC  + PspecificPO4(np)*R_PC(np)
               consumpDIC_PIC = consumpDIC_PIC + 
     &                           PspecificPO4(np)*R_PC(np)*R_PICPOC(np)
#endif
            enddo
           if (debug.eq.7) print*,'local', parlocal,tlocal,po4local,
     &                  no3local, no2local,nh4local,fetlocal,silocal
           if (debug.eq.7) print*,'grow',grow
           if (debug.eq.6) print*,'pspecificpo4', PspecificPO4
           if (debug.eq.6) print*,'consumpPO4', consumpPO4
           if (debug.eq.6) print*,'consumpFeT', consumpFeT
           if (debug.eq.6) print*,'consumpSi ', consumpsi 
           if (debug.eq.6) print*,'consumpNO3', consumpNO3
           if (debug.eq.6) print*,'consumpNO2', consumpNO2
           if (debug.eq.6) print*,'consumpNH4', consumpNH4
c ****************** END NUTRIENT UPTAKE ****************************

c sinking phytoplankton and POM
             if(bottom .eq. 1.0 _d 0)then
              psinkP   = (wp_sink*POPuplocal)/(dzlocal)
              psinkN   = (wn_sink*PONuplocal)/(dzlocal)
              psinkFe  = (wfe_sink*POFeuplocal)/(dzlocal)
              psinkSi  = (wsi_sink*PSiuplocal)/(dzlocal)
              psinkSih = (wsi_sink*PSihuplocal)/(dzlocal)
              do np=1,npmax
               psinkPhy(np) =
     &                     (wsink(np)*Phytoup(np))/(dzlocal)
              enddo
#ifdef DYNAMIC_CHL
              do np=1,npmax
               psinkChl(np) =
     &                     (wsink(np)*Chlup(np))/(dzlocal)
              enddo
#endif
#ifdef ALLOW_CARBON
              psinkC    = (wc_sink*POCuplocal)/(dzlocal)    ! labile POC
              psinkChp  = (wc_sink*POChpuplocal)/(dzlocal)  ! POC assoc with hard PIC 
              psinkChs  = (wc_sink*POChsuplocal)/(dzlocal)  ! POC assoc with hard opal
              psinkCsp  = (wc_sink*POCspuplocal)/(dzlocal)  ! POC assoc with soft PIC
              psinkCss  = (wc_sink*POCssuplocal)/(dzlocal)  ! POC assoc with soft opal
              psinkPIC  = (wpic_sink*PICuplocal)/(dzlocal)  ! soft PIC
              psinkPICh = (wpic_sink*PIChuplocal)/(dzlocal) ! hard PIC
              psinkChd  = (wc_sink*POChduplocal)/(dzlocal)  ! POC assoc with hard dust
              psinkCsd  = (wc_sink*POCsduplocal)/(dzlocal)  ! POC assoc with soft dust
              psinkDUSTh = (wc_sink*dusthuplocal)/(dzlocal) ! hard dust
              psinkDUSTs = (wc_sink*dustsuplocal)/(dzlocal) ! soft dust
#endif
            else
              psinkP   = (wp_sink*(POPuplocal-POPlocal))/(dzlocal)
              psinkN   = (wn_sink*(PONuplocal-PONlocal))/(dzlocal)
              psinkFe  = (wfe_sink*(POFeuplocal-POFelocal))/(dzlocal)
              psinkSi  = (wsi_sink*(PSiuplocal-PSilocal))/(dzlocal)
              psinkSih = (wsi_sink*(PSihuplocal-PSihlocal))/(dzlocal)
              do np=1,npmax
               psinkPhy(np) =
     &                  (wsink(np))*(Phytoup(np)-Phyto(np))/(dzlocal)
              enddo
#ifdef DYNAMIC_CHL
              do np=1,npmax
               psinkChl(np) =
     &                  (wsink(np))*(Chlup(np)-phychl(np))/(dzlocal)
              enddo
#endif
#ifdef ALLOW_CARBON
c             labile POC
              psinkC     = (wc_sink*(POCuplocal-POClocal))/(dzlocal)   
c             POC assoc with hard PIC and hard opal 
              psinkChp   = (wc_sink*(POChpuplocal-POChplocal))/(dzlocal)
              psinkChs   = (wc_sink*(POChsuplocal-POChslocal))/(dzlocal)
c             POC assoc with soft PIC and soft opal
              psinkCsp   = (wc_sink*(POCspuplocal-POCsplocal))/(dzlocal)
              psinkCss   = (wc_sink*(POCssuplocal-POCsslocal))/(dzlocal)
c             soft and hard PIC
              psinkPIC   = (wpic_sink*(PICuplocal-PIClocal))/(dzlocal)
              psinkPICh  = (wpic_sink*(PIChuplocal-PIChlocal))/(dzlocal)
c             POC assoc with hard and soft dust
              psinkChd   = (wc_sink*(POChduplocal-POChdlocal))/(dzlocal)
              psinkCsd   = (wc_sink*(POCsduplocal-POCsdlocal))/(dzlocal)
c             hard and soft dust
              psinkDUSTh = (wc_sink*(dusthuplocal-dusthlocal))/(dzlocal)
              psinkDUSTs = (wc_sink*(dustsuplocal-dustslocal))/(dzlocal)
#endif
            endif

c DOM remineralization rates
            DOPremin =  reminTempFunction * Kdop * DOPlocal
            DONremin =  reminTempFunction * Kdon * DONlocal
            DOFeremin =  reminTempFunction * KdoFe * DOFelocal

c remineralization of sinking particulate
            preminP   = reminTempFunction * Kpremin_P   * POPlocal
            preminN   = reminTempFunction * Kpremin_N   * PONlocal
            preminFe  = reminTempFunction * Kpremin_Fe  * POFelocal

c remineralization of soft and hard opal
            preminSi  = reminTempFunction * Kpremin_Si  * PSilocal
            preminSih = Kpremin_hard* PSihlocal
#ifdef ALLOW_CDOM
            preminP_cdom = fraccdom*preminP
            preminP=preminP-preminP_cdom
            preminN_cdom = rnp_cdom*preminP_cdom
            preminN=preminN-preminN_cdom
            preminFe_cdom = rfep_cdom*preminP_cdom
            preminFe=preminFe-preminFe_cdom
#endif

#ifdef ALLOW_CARBON
           DOCremin   = reminTempFunction * Kdoc        * DOClocal
           preminC    = reminTempFunction * Kpremin_C   * POClocal
c          POC assoc with hard PIC and hard opal
           preminChp  = Kpremin_hard* POChplocal 
           preminChs  = Kpremin_hard* POChslocal
c          POC assoc with soft PIC and soft opal 
           preminCsp  = Kdissc * POCsplocal
           preminCss  = reminTempFunction * Kpremin_Si  * POCsslocal
c          POC assoc with hard and soft dust
           preminChd = Kpremin_DUSTh * POChdlocal
           preminCsd = Kpremin_DUSTs * POCsdlocal
c          remineralization of hard and soft dust
           preminDUSTh = Kpremin_DUSTh * dusthlocal
           preminDUSTs = Kpremin_DUSTs * dustslocal
#ifdef ALLOW_CDOM
            preminC_cdom = rcp_cdom*preminP_cdom
            preminC=preminC-preminC_cdom
#endif

c dissolution of soft and hard PIC
           disscPIC  = Kdissc*PIClocal
           disscPICh = Kdissc_hard*PIChlocal
#endif

#ifdef ALLOW_CDOM
c degradation of  CDOM - high when bleached by light
           cdomp_degrd = reminTempFunction * cdomlocal
     &            *(cdomdegrd+cdombleach*min(PARlocal/PARcdom,1. _d 0) )
           cdomn_degrd = rnp_cdom * cdomp_degrd
           cdomfe_degrd= rfep_cdom * cdomp_degrd
#ifdef ALLOW_CARBON
           cdomc_degrd = rcp_cdom  * cdomp_degrd
#endif
#endif

c chemistry
c NH4 -> NO2 -> NO3 by bacterial action
            NO2prod = knita*( 1. _d 0-min(PARlocal/PAR0,1. _d 0) )
     &                *NH4local
            NO3prod = knitb*( 1. _d 0-min(PARlocal/PAR0,1. _d 0) )
     &                *NO2local 
c           NO2prod = knita*NH4local
c           NO3prod = knitb*NO2local
c
#ifdef PART_SCAV
            scav_poc=POPlocal/1.1321 _d -4
c scav rate 
            scav_part=scav_rat*scav_inter*(scav_poc**scav_exp)
#endif
c -------------------------------------------------------------------
c calculate tendency terms (and some diagnostics)
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c phytoplankton
            do np=1,npmax
              dphytodt(np) =   PspecificPO4(np)
#ifdef OLD_GRAZE
     &                        - grazing_phyto(np)*
     &              (phytomin(np)/(phytomin(np) + kgrazesat))
#else
     &                        - sumgrazphy(np)
#endif
     &                        - mortphy(np)*phytomin(np)
     &                        + psinkphy(np)
#ifdef GEIDER
#ifdef DYNAMIC_CHL
              dphychl(np) = acclim(np)*PspecificPO4(np)*R_PC(np)
c             dphychl(np) = rhochl(np)*PspecificPO4(np)*R_PC(np) 
     &                     + acclimtimescl *
     &                      (acclim(np)-chl2c(np))*phyto(np)*R_PC(np)
     &                     +(
#ifdef OLD_GRAZE
     &                        - grazing_phyto(np)*
     &              (phytomin(np)/(phytomin(np) + kgrazesat))
#else
     &                   - sumgrazphy(np)
#endif
     &                   - mortphy(np)*phytomin(np)) *chl2c(np)*R_PC(np)
     &                        + psinkChl(np) 
#endif
              Chl=Chl + phychl(np)
#endif
c %% diagnostics
              PP = PP + PspecificPO4(np)
c%%% 
#ifdef OLD_GRAZE
              tmpr=grazing_phyto(np)*
     &              (phytomin(np)/(phytomin(np) + kgrazesat)) 
     &                        + mortphy(np)*phytomin(np)
     &                        - psinkphy(np)
#else
              tmpr=sumgrazphy(np)
     &                        + mortphy(np)*phytomin(np)
     &                        - psinkphy(np)
#endif
#ifdef DAR_DIAG_RSTAR
#ifndef GEIDER
              tmpgrow=ngrow(np)*mu(np)*ilimit(np)*
     &              phytoTempFunction(np)
#endif
#ifdef GEIDER
              tmpgrow=grow(np)/limit(np)
#endif
              tmp1=tmpgrow*phyto(np)-tmpr
              tmp2=tmpgrow*phyto(np)*(exp(-sig1*nh4local)+NH4limit(np))
     &                 -tmpr 
              if (tmp1.ne.0. _d 0) then 
                Rstarlocal(np)=ksatPO4(np)*tmpr/tmp1
              else
                Rstarlocal(np)=-9999. _d 0
              endif
              if (tmp2.ne.0. _d 0) then 
                RNstarlocal(np)=ksatNO3(np)*
     &                 (tmpr-tmpgrow*NH4limit(np)*phyto(np))/tmp2
              else
                RNstarlocal(np)=-9999. _d 0
              endif
#endif
#ifdef DAR_DIAG_GROW
c include temp, light, nutrients
c             Growlocal(np)=grow(np)
c include temp and light, but not nutrients
              Growlocal(np)=ngrow(np)*mu(np)*ilimit(np)*
     &              phytoTempFunction(np)
c include temp, but not nutrients or light
c             Growlocal(np)=ngrow(np)*mu(np)*
c    &              phytoTempFunction(np)
              Growsqlocal(np)=Growlocal(np)**2
#endif
            enddo
c end np loop
            if (debug.eq.10) print*,'dphytodt',dphytodt
c
#ifdef OLD_GRAZE
c zooplankton growth by grazing
            do nz=1,nzmax
c zoo in P currency
              dzooPdt(nz)  = grazingP(nz)*zooP(nz) 
C zooplankton stoichiometry varies according to food source
              dzooNdt(nz)  = grazingN(nz)*zooP(nz) 
              dzooFedt(nz) = grazingFe(nz)*zooP(nz) 
              dzooSidt(nz) = grazingSi(nz)*zooP(nz) 
            enddo
#else
            do nz=1,nzmax
c zoo in P currency
              dzooPdt(nz)  = sumgrazzoo(nz)
C zooplankton stoichiometry varies according to food source
              dzooNdt(nz)  = sumgrazzooN(nz)
              dzooFedt(nz) = sumgrazzooFe(nz)
              dzooSidt(nz) = sumgrazzooSi(nz)
            enddo
#endif
            if (debug.eq.10) print*,'dZooPdt',dZooPdt

c zooplankton mortality
            do nz=1,nzmax
c zoo in P currency
              dzooPdt(nz)  = dzooPdt(nz)
     &                           - mortzoo(nz)*zooP(nz)
     &                           - mortzoo2(nz)*zooP(nz)**2
c zooplankton in other currencies
C zooplankton stoichiometry varies according to food source
              dzooNdt(nz)  = dzooNdt(nz)
     &                           - mortzoo(nz)*zooN(nz)
     &                           - mortzoo2(nz)*zooN(nz)**2
              dzooFedt(nz) = dzooFedt(nz)
     &                           - mortzoo(nz)*zooFe(nz)
     &                           - mortzoo2(nz)*zooFe(nz)**2
              dzooSidt(nz) = dzooSidt(nz)
     &                           - mortzoo(nz)*zooSi(nz)
     &                           - mortzoo2(nz)*zooSi(nz)**2
            enddo


c sum contributions to inorganic nutrient tendencies
            dPO4dt =  - consumpPO4 +  
#ifdef ALLOW_CDOM
     &                 DOPremin
#else
     &                 preminP  + DOPremin
#endif
            dNH4dt =  - consumpNH4 +  
#ifdef ALLOW_CDOM
     &                 DONremin
#else
     &                 preminN  +  DONremin 
#endif
     &                                   - NO2prod
            dNO2dt =  - consumpNO2                        
     &                                   + NO2prod - NO3prod
            dNO3dt =  - consumpNO3             
     &                                   + NO3prod
c-ONLYNO3   dNO3dt =  - consumpNO3 +  preminN  + DONremin
#ifdef ALLOW_DENIT
            if (O2local.le.O2crit) then
              denit =  denit_np
#ifdef ALLOW_CDOM
     &                *(DOPremin)
#else 
     &                *(preminP  + DOPremin)
#endif
              dNO3dt = dNO3dt - denit
              dNH4dt = dNH4dt -  
#ifdef ALLOW_CDOM
     &                 (DONremin)
#else
     &                 (preminN  + DONremin)
#endif
            endif
#endif
            dFeTdt =  - consumpFeT +  
#ifdef ALLOW_CDOM
     &                 DOFeremin
#else
     &                 preminFe + DOFeremin
#endif
#ifdef PART_SCAV
     &                 - scav_part*freefelocal +
#else
     &                 - scav*freefelocal + 
#endif          
     &                 alpfe*inputFelocal/dzlocal
            dSidt  =  - consumpSi  +  preminSi

c tendency of dissolved organic pool
            dDOPdt  = totphy_dop  + totzoo_dop  - DOPremin 
#ifdef ALLOW_CDOM
     &                +preminP + cdomp_degrd
#endif
            dDONdt  = totphy_don  + totzoo_don  - DONremin 
#ifdef ALLOW_CDOM
     &                +preminN + cdomn_degrd
#endif
            dDOFedt = totphy_dofe + totzoo_dofe - DOFeremin 
#ifdef ALLOW_CDOM
     &                +preminFe + cdomfe_degrd
#endif

c tendency of particulate detritus pools
            dpopdt  = totphy_pop  + totzoo_pop   - preminP + psinkP
#ifdef ALLOW_CDOM
     &                -preminP_cdom
#endif
            dpondt  = totphy_pon  + totzoo_pon   - preminN + psinkN
#ifdef ALLOW_CDOM
     &                -preminN_cdom
#endif
            dpofedt = totphy_pofe + totzoo_pofe  - preminFe + psinkFe
#ifdef ALLOW_CDOM
     &                -preminFe_cdom
#endif

c ========================================================================
c ballast hypothesis: 
c based on Armstrong 2002. Followed Moore et al. 2004 and Lima et al. 2014
c L. Gloege 2015
c ========================================================================

c-----------------------------------------------------------------------
c production tendencies
c (7% of small phtyo. produce PIC, set in generatate_phyto)
c-----------------------------------------------------------------------
            dPOCproddt  = totphy_poc  + totzoo_poc    ! production POC
            dPICproddt  = totphy_pic  + totzoo_pic    ! production PIC
            dpSiproddt  = totphy_posi + totzoo_posi   ! production opal

c PIC hard,soft production
            dPIChproddt = f_pcaco3_hard     * dPICproddt
            dPICsproddt = (1-f_pcaco3_hard) * dPICproddt

c OPAL hard,soft production
            dpSihproddt = f_psio3_hard     * dpSiproddt
            dpSisproddt = (1-f_psio3_hard) * dpSiproddt

c-----------------------------------------------------------------------
c ballast associated POC (non-free component
c----------------------------------------------------------------------
            dPOCnfproddt = (omega_pcaco3 * dPICproddt) +
     &                     (omega_psio3  * dpSiproddt)

c-----------------------------------------------------------------------
c POC free production
c make sure there enough produced POC to be associated with ballast minerals
c-----------------------------------------------------------------------
            if ( dPOCproddt.ge.dPOCnfproddt ) then
               dPOCavaildt  = dPOCproddt   - dPOCnfproddt  ! free POC
               dPOChpproddt = omega_pcaco3 * dPIChproddt   ! h POC-PIC
               dPOChsproddt = omega_psio3  * dpSihproddt   ! h POC-opal
               dPOCspproddt = omega_pcaco3 * dPICsproddt   ! s POC-PIC
               dPOCssproddt = omega_psio3  * dpSisproddt   ! s POC-opal
            else
               dPOCavaildt  = dPOCproddt
               dPOCnfproddt = 0. _d 0
               dPOChpproddt = 0. _d 0
               dPOChsproddt = 0. _d 0
               dPOCspproddt = 0. _d 0
               dPOCssproddt = 0. _d 0
            endif


c           WRITE( *,'(A,D)') 'lukedPOCqas1 ', dPOCqasdt

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c------------------------------------------------------------------------
c           set outgoing POC associated with dust
c------------------------------------------------------------------------
c           -------------------------------------------------
c           DUST ASSOCIATED POC .less than. AVAILABLE
c           -------------------------------------------------
            if (dPOCqadt.le.dPOCavaildt) then 
               dPOCavaildt  = dPOCavaildt - dPOCqadt  ! leftover is labile
               dPOCqahdt    = dPOCqahindt             ! POC associated with hard dust
               dPOCqasdt    = dPOCqasindt             ! POC associated with soft dust
            endif

c           -------------------------------------------------
c           DUST ASSOCIATED POC is ZERO
c           -------------------------------------------------
            if (dPOCqadt.eq.0. _d 0) then
               dPOCqahdt = 0   ! POC associated with hard dust
               dPOCqasdt = 0   ! POC associated with soft dust
            endif

c           -------------------------------------------------
c           DUST ASSOCIATED POC .greater than. AVAILABLE
c           then split available POC between hard and soft 
c           class associated with dust and leave none as labile
c           -------------------------------------------------
            if (dPOCqadt.gt.dPOCavaildt) then
c              split available across hard and soft
               dPOCqahindt  = (f_dust_hard) * dPOCavaildt 
               dPOCqasindt  = (1 - f_dust_hard) * dPOCavaildt
               dPOCqahdt    = dPOCqahindt   ! POC associated with hard dust
               dPOCqasdt    = dPOCqasindt   ! POC associated with soft dust
c              now set available POC to zero
               dPOCavaildt = 0
            endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c------------------------------------------------------------------------
c tendency of hard and soft dust
c------------------------------------------------------------------------
            dDUSThdt = -preminDUSTh + psinkDUSTh
            dDUSTsdt = -preminDUSTs + psinkDUSTs

c------------------------------------------------------------------------
c tendency of hard and soft opal
c------------------------------------------------------------------------
            dpSidt   = dpSisproddt  - preminSi  + psinkSi
            dpSihdt  = dpSihproddt  - preminSih + psinkSih

#ifdef ALLOW_CARBON
            dDICdt =  - consumpDIC - consumpDIC_PIC 
#ifdef ALLOW_CDOM
     &                                         +  DOCremin
#else
     &                                         +  preminC  + DOCremin
     &                                         + preminChp + preminChs 
     &                                         + preminCsp + preminCss
     &                                         + preminChd + preminCsd
#endif
     &                                         + disscPIC
     &                                         + disscPICh ! gloege added this
            dDOCdt  = totphy_doc  + totzoo_doc  - DOCremin
#ifdef ALLOW_CDOM
     &                    +preminC + cdomc_degrd
#endif
c------------------------------------------------------------------------
c tendency of POC, labile and ballast associated
c------------------------------------------------------------------------
            dPOChpdt = dPOChpproddt - preminChp + psinkChp ! hard PIC POC
            dPOChsdt = dPOChsproddt - preminChs + psinkChs ! soft PIC POC
            dPOCspdt = dPOCspproddt - preminCsp + psinkCsp ! hard opal POC
            dPOCssdt = dPOCssproddt - preminCss + psinkCss ! soft opal POC
            dPOChddt = dPOCqahdt    - preminChd + psinkChd ! hard dust POC
            dPOCsddt = dPOCqasdt    - preminCsd + psinkCsd ! soft dust POC
            dPOCdt   = dPOCavaildt  - preminC   + psinkC   ! labile POC

c original implementation, exponetial decay
c            dPOCdt  = totphy_poc  + totzoo_poc  - preminC + psinkC
#ifdef ALLOW_CDOM
     &                -preminC_cdom
#endif
c------------------------------------------------------------------------
c tendency of hard and soft PIC
c------------------------------------------------------------------------
            dPICdt  = dPICsproddt  - disscPIC  + psinkPIC
            dPIChdt = dPIChproddt  - disscPICh + psinkPICh

c NOTHING BELOW IS RELATED TO BALLAST HYPTOHESIS (L. Gloege)
c======================================================================== 

c L. Gloege, added on the dissc PICh hard term
C
            dALKdt  = - dNO3dt - 2.d0 * (consumpDIC_PIC - disscPIC 
     &                                                  - disscPICh)
c should be = O2prod - preminP - DOPremin?
c OLD WAY
c            dO2dt   = - R_OP*dPO4dt
c production of O2 by photosynthesis
             dO2dt   = R_OP*consumpPO4
c loss of O2 by remineralization
             if (O2local.gt.O2crit) then
               dO2dt   = dO2dt - R_OP
#ifdef ALLOW_CDOM
     &                      *(DOPremin)
#else
     &                      *(preminP  + DOPremin)
#endif
             endif
#ifdef OLD_GRAZE
            do nz=1,nzmax
              dzooCdt(nz)  = grazingC(nz)*zooClocal(nz)
     &                           - mortzoo(nz)*zooClocal(nz)
     &                           - mortzoo2(nz)*zooClocal(nz)**2
            enddo
#else
            do nz=1,nzmax
              dzooCdt(nz)  = sumgrazzooc(nz)
     &                           - mortzoo(nz)*zooClocal(nz)
     &                           - mortzoo2(nz)*zooClocal(nz)**2
            enddo
#endif

#ifdef ALLOW_CDOM
            dcdomdt = totphy_cdomp + totzoo_cdomp +  preminP_cdom
     &                  -cdomp_degrd
#endif

#endif

            if (debug.eq.10) print*,'dDOPdt', dDOPdt
            if (debug.eq.10) print*,'dpopdt',dpopdt
            if (debug.eq.10) print*,'dDONdt',dDONdt
            if (debug.eq.10) print*,'dpondt',dpondt
c
c -------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c --------------------------------------------------------------------------

c -m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-
c  Mutation - apply mutation to tendencies [jbmodif]

#ifdef ALLOW_MUTANTS
c apply to all sisters when first sister is encountered
       if(runtim .gt. threeyr) then
       mutfor=1 _d -8
       mutback=1 _d -12
       if(numtax .gt. 1)then
        do np=1,npro
          if(mod(np,numtax).eq. 1. _d 0)then
             nsisone   = np
             nsistwo   = np+1 
             nsisthree = np+2             
             nsisfour  = np+3 
              
             grow1 = PspecificPO4(nsisone)
             grow2 = PspecificPO4(nsistwo)

             if(numtax.eq.2)grow3 = 0.0 _d 0
             if(numtax.eq.2)grow4 = 0.0 _d 0

             if(numtax.eq.3)grow4 = 0.0 _d 0
             if(numtax.ge.3)grow3 = PspecificPO4(nsisthree)
             
             if(numtax.eq.4)grow4 = PspecificPO4(nsisfour)
                        


                dphytodt(nsisone)  = dphytodt(nsisone)
     &                - grow1 *1.4427 _d 0*mutfor
     &                - grow1 *1.4427 _d 0*mutfor
     &                - grow1 *1.4427 _d 0*mutfor
     &                + grow2 *1.4427 _d 0*mutback
     &                + grow3 *1.4427 _d 0*mutback
     &                + grow4 *1.4427 _d 0*mutback

                dphytodt(nsistwo)  = dphytodt(nsistwo)
     &                - grow2 *1.4427 _d 0*mutback
     &                + grow1 *1.4427 _d 0*mutfor

                if(numtax .ge. 3)then
                  dphytodt(nsisthree)  = dphytodt(nsisthree)
     &                - grow3 *1.4427 _d 0*mutback
     &                + grow1 *1.4427 _d 0*mutfor
                endif
        
                if(numtax .eq. 4)then
                  dphytodt(nsisfour)  = dphytodt(nsisfour)
     &                - grow4 *1.4427 _d 0*mutback
     &                + grow1 *1.4427 _d 0*mutfor
c QQQQQQQQQQ FIX FOR NIT RUNS ONLY!!!
                  if (phyto(nsisfour).eq.0. _d 0) then
                    if (phyto(nsistwo).eq.0. _d 0) then
                     if (dphytodt(nsistwo).gt.dphytodt(nsisfour)) then
                         dphytodt(nsisfour)=dphytodt(nsistwo)
                     endif
                    endif
                    if (phyto(nsisthree).eq.0. _d 0) then
                     if (dphytodt(nsisthree).gt.dphytodt(nsisfour)) then
                         dphytodt(nsisfour)=dphytodt(nsisthree)
                     endif
                    endif
                  endif
c QQQQQQQQQQQQQ
                endif

c QQQQQQQQQQQQTEST
         if (debug.eq.11) then
         if (PARlocal.gt.1. _d 0) then
          if (dphytodt(nsistwo).gt.dphytodt(nsisfour).and.
     &                             dphytodt(nsisfour).gt.0. _d 0) then
                print*,'QQQQ nsistwo>nsisfour',nsistwo,nsisfour, 
     &               dphytodt(nsistwo),  dphytodt(nsisfour),
     &               phyto(nsistwo), phyto(nsisfour),  
     &               phyto(nsisone)
          endif
          if (dphytodt(nsisthree).gt.dphytodt(nsisfour).and.
     &                             dphytodt(nsisfour).gt.0. _d 0) then
                 print*,'QQQQ nsisthree>nsisfour',nsisthree,nsisfour,
     &                dphytodt(nsisthree), dphytodt(nsisfour),
     &               phyto(nsisthree), phyto(nsisfour),  
     &               phyto(nsisone)
          endif
          if (dphytodt(nsisfour).gt.dphytodt(nsisone).and.
     &                            dphytodt(nsisone).gt.0. _d 0) then
                 print*,' BIG QQQQ nsisfour>nsisone',nsisone,nsisfour,
     &                dphytodt(nsisfour), dphytodt(nsisone),
     &               phyto(nsisfour), phyto(nsisone)
          endif
         endif
         endif
c QQQQQQQQQTEST
             endif
          enddo
        endif
        endif 

c mutation is finished
c -m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-m-
#endif



         RETURN
         END
#endif  /*MONOD*/
#endif  /*ALLOW_PTRACERS*/
c ==================================================================

