C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/darwin/darwin_mnc_init.F,v 1.1 2011/04/13 18:56:24 jahn Exp $
C $Name:  $

#include "DARWIN_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: DARWIN_MNC_INIT

C     !INTERFACE:
      SUBROUTINE DARWIN_MNC_INIT( myThid )

C     !DESCRIPTION:
C     Initialize MNC for DARWIN

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "DARWIN_SIZE.h"
#include "SPECTRAL_SIZE.h"

C     !INPUT PARAMETERS:
C     myThid               :: thread number
      INTEGER myThid
CEOP
      INTEGER np
      INTEGER ii,jj,ilnb
      CHARACTER*(MAX_LEN_MBUF) suff
      CHARACTER*(1) pf
      CHARACTER*(MAX_LEN_FNAM) name
#ifdef ALLOW_QUOTA
      INTEGER iPhy
#endif
C     Functions
      integer ILNBLNK

#ifdef ALLOW_DARWIN
#ifdef ALLOW_MNC

      IF ( useMNC ) THEN

#ifdef ALLOW_QUOTA
      iPhy=ibiomass
#endif

      CALL MNC_CW_ADD_VNAME(
     &     'dar_PAR_ave', 'Cen_xy_Hn__C__t', 4,5, myThid)
      CALL MNC_CW_ADD_VATTR_TEXT(
     &     'dar_PAR_ave','units','--', myThid)
      CALL MNC_CW_ADD_VATTR_TEXT('dar_PAR_ave','description',
     &     '', myThid)

      CALL MNC_CW_ADD_VNAME(
     &     'dar_PP_ave', 'Cen_xy_Hn__C__t', 4,5, myThid)
      CALL MNC_CW_ADD_VATTR_TEXT(
     &     'dar_PP_ave','units','--', myThid)
      CALL MNC_CW_ADD_VATTR_TEXT('dar_PP_ave','description',
     &     '', myThid)
#ifndef ALLOW_QUOTA
      CALL MNC_CW_ADD_VNAME(
     &     'dar_Chl_ave', 'Cen_xy_Hn__C__t', 4,5, myThid)
      CALL MNC_CW_ADD_VATTR_TEXT(
     &     'dar_Chl_ave','units','--', myThid)
      CALL MNC_CW_ADD_VATTR_TEXT('dar_Chl_ave','description',
     &     '', myThid)
#endif
      CALL MNC_CW_ADD_VNAME(
     &     'dar_Nfix_ave', 'Cen_xy_Hn__C__t', 4,5, myThid)
      CALL MNC_CW_ADD_VATTR_TEXT(
     &     'dar_Nfix_ave','units','--', myThid)
      CALL MNC_CW_ADD_VATTR_TEXT('dar_Nfix_ave','description',
     &     '', myThid)

      CALL MNC_CW_ADD_VNAME(
     &     'dar_Denit_ave', 'Cen_xy_Hn__C__t', 4,5, myThid)
      CALL MNC_CW_ADD_VATTR_TEXT(
     &     'dar_Denit_ave','units','--', myThid)
      CALL MNC_CW_ADD_VATTR_TEXT('dar_Denit_ave','description',
     &     '', myThid)

#ifdef DAR_DIAG_ACDOM
      CALL MNC_CW_ADD_VNAME(
     &     'dar_aCDOM_ave', 'Cen_xy_Hn__C__t', 4,5, myThid)
      CALL MNC_CW_ADD_VATTR_TEXT(
     &     'dar_aCDOM_ave','units','--', myThid)
      CALL MNC_CW_ADD_VATTR_TEXT('dar_aCDOM_ave','description',
     &     '', myThid)
#endif

#ifdef DAR_DIAG_IRR
      do ii=1,tlam
          DO jj = 1,MAX_LEN_FNAM
             name(jj:jj) = ' '
          ENDDO
          WRITE(name,'(a,i2.2)') 'dar_Ed',ii
          CALL MNC_CW_ADD_VNAME(name,'Cen_xy_Hn__C__t',4,5,myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'units','--', myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'description','', myThid)
          WRITE(name,'(a,i2.2)') 'dar_Es',ii
          CALL MNC_CW_ADD_VNAME(name,'Cen_xy_Hn__C__t',4,5,myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'units','--', myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'description','', myThid)
          WRITE(name,'(a,i2.2)') 'dar_Eu',ii
          CALL MNC_CW_ADD_VNAME(name,'Cen_xy_Hn__C__t',4,5,myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'units','--', myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'description','', myThid)
          WRITE(name,'(a,i2.2)') 'dar_Eutop',ii
          CALL MNC_CW_ADD_VNAME(name,'Cen_xy_Hn__C__t',4,5,myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'units','--', myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'description','', myThid)
      enddo
#endif

#ifdef DAR_DIAG_ABSORP
      do ii=1,tlam
          DO jj = 1,MAX_LEN_FNAM
             name(jj:jj) = ' '
          ENDDO
          WRITE(name,'(a,i2.2)') 'dar_a',ii
          CALL MNC_CW_ADD_VNAME(name,'Cen_xy_Hn__C__t',4,5,myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'units','--', myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'description','', myThid)
      enddo
#endif
#ifdef DAR_DIAG_SCATTER
      do ii=1,tlam
          DO jj = 1,MAX_LEN_FNAM
             name(jj:jj) = ' '
          ENDDO
          WRITE(name,'(a,i2.2)') 'dar_bt',ii
          CALL MNC_CW_ADD_VNAME(name,'Cen_xy_Hn__C__t',4,5,myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'units','--', myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'description','', myThid)
          WRITE(name,'(a,i2.2)') 'dar_bb',ii
          CALL MNC_CW_ADD_VNAME(name,'Cen_xy_Hn__C__t',4,5,myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'units','--', myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'description','', myThid)
      enddo
#endif
#ifdef DAR_DIAG_PART_SCATTER
      do ii=1,tlam
          DO jj = 1,MAX_LEN_FNAM
             name(jj:jj) = ' '
          ENDDO
          WRITE(name,'(a,i2.2)') 'dar_apart',ii
          CALL MNC_CW_ADD_VNAME(name,'Cen_xy_Hn__C__t',4,5,myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'units','--', myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'description','', myThid)
          WRITE(name,'(a,i2.2)') 'dar_btpart',ii
          CALL MNC_CW_ADD_VNAME(name,'Cen_xy_Hn__C__t',4,5,myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'units','--', myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'description','', myThid)
          WRITE(name,'(a,i2.2)') 'dar_bbpart',ii
          CALL MNC_CW_ADD_VNAME(name,'Cen_xy_Hn__C__t',4,5,myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'units','--', myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'description','', myThid)
      enddo
#endif

c ANNA_TAVE
#ifdef WAVES_DIAG_PCHL
      do np=1,npmax
          ii=iPhy+np-1
          DO jj = 1,MAX_LEN_FNAM
             name(jj:jj) = ' '
          ENDDO
          ilnb = ILNBLNK(PTRACERS_names(ii))
          WRITE(name,'(a,a)') 'dar_Pchl_',PTRACERS_names(ii)(1:ilnb)
          CALL MNC_CW_ADD_VNAME(name,'Cen_xy_Hn__C__t',4,5,myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'units','--', myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'description','', myThid)
      enddo
#endif
c ANNA end TAVE
#ifdef DAR_DIAG_RSTAR
      do np=1,npmax
          ii=iPhy+np-1
          DO jj = 1,MAX_LEN_FNAM
             name(jj:jj) = ' '
          ENDDO
          ilnb = ILNBLNK(PTRACERS_names(ii))
          WRITE(name,'(a,a)') 'dar_Rstar_',PTRACERS_names(ii)(1:ilnb)
          CALL MNC_CW_ADD_VNAME(name,'Cen_xy_Hn__C__t',4,5,myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'units','--', myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'description','', myThid)
      enddo
      do np=1,npmax
          ii=iPhy+np-1
          DO jj = 1,MAX_LEN_FNAM
             name(jj:jj) = ' '
          ENDDO
          ilnb = ILNBLNK(PTRACERS_names(ii))
          WRITE(name,'(a,a)') 'dar_RNstar_',PTRACERS_names(ii)(1:ilnb)
          CALL MNC_CW_ADD_VNAME(name,'Cen_xy_Hn__C__t',4,5,myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'units','--', myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'description','', myThid)
      enddo
#endif

#ifdef DAR_DIAG_DIVER
      CALL MNC_CW_ADD_VNAME(
     &     'dar_Diver1_ave', 'Cen_xy_Hn__C__t', 4,5, myThid)
      CALL MNC_CW_ADD_VATTR_TEXT(
     &     'dar_Diver1_ave','units','--', myThid)
      CALL MNC_CW_ADD_VATTR_TEXT('dar_Diver1_ave','description',
     &     '', myThid)

      CALL MNC_CW_ADD_VNAME(
     &     'dar_Diver2_ave', 'Cen_xy_Hn__C__t', 4,5, myThid)
      CALL MNC_CW_ADD_VATTR_TEXT(
     &     'dar_Diver2_ave','units','--', myThid)
      CALL MNC_CW_ADD_VATTR_TEXT('dar_Diver2_ave','description',
     &     '', myThid)

      CALL MNC_CW_ADD_VNAME(
     &     'dar_Diver3_ave', 'Cen_xy_Hn__C__t', 4,5, myThid)
      CALL MNC_CW_ADD_VATTR_TEXT(
     &     'dar_Diver3_ave','units','--', myThid)
      CALL MNC_CW_ADD_VATTR_TEXT('dar_Diver3_ave','description',
     &     '', myThid)

      CALL MNC_CW_ADD_VNAME(
     &     'dar_Diver4_ave', 'Cen_xy_Hn__C__t', 4,5, myThid)
      CALL MNC_CW_ADD_VATTR_TEXT(
     &     'dar_Diver4_ave','units','--', myThid)
      CALL MNC_CW_ADD_VATTR_TEXT('dar_Diver4_ave','description',
     &     '', myThid)

#endif
#ifdef DAR_DIAG_GROW
      do np=1,npmax
          ii=iPhy+np-1
          DO jj = 1,MAX_LEN_FNAM
             name(jj:jj) = ' '
          ENDDO
          ilnb = ILNBLNK(PTRACERS_names(ii))
          WRITE(name,'(a,a)') 'dar_Grow_',PTRACERS_names(ii)(1:ilnb)
          CALL MNC_CW_ADD_VNAME(name,'Cen_xy_Hn__C__t',4,5,myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'units','--', myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'description','', myThid)
      enddo
      do np=1,npmax
          ii=iPhy+np-1
          DO jj = 1,MAX_LEN_FNAM
             name(jj:jj) = ' '
          ENDDO
          ilnb = ILNBLNK(PTRACERS_names(ii))
          WRITE(name,'(a,a)') 'dar_Growsq_',PTRACERS_names(ii)(1:ilnb)
          CALL MNC_CW_ADD_VNAME(name,'Cen_xy_Hn__C__t',4,5,myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'units','--', myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'description','', myThid)
      enddo
#endif

#ifdef ALLOW_DIAZ
#ifdef DAR_DIAG_NFIXP
      do np=1,npmax
          ii=iPhy+np-1
          DO jj = 1,MAX_LEN_FNAM
             name(jj:jj) = ' '
          ENDDO
          ilnb = ILNBLNK(PTRACERS_names(ii))
          WRITE(name,'(a,a)') 'dar_NfixP_',PTRACERS_names(ii)(1:ilnb)
          CALL MNC_CW_ADD_VNAME(name,'Cen_xy_Hn__C__t',4,5,myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'units','--', myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(name,'description','', myThid)
      enddo
#endif
#endif

      

c     CALL MNC_CW_ADD_VNAME(
c    &     'dic_pH_ave', 'Cen_xy_Hn__-__t', 3,4, myThid)
c     CALL MNC_CW_ADD_VATTR_TEXT(
c    &     'dic_pH_ave','units','--', myThid)
c     CALL MNC_CW_ADD_VATTR_TEXT('dic_pH_ave','description',
c    &     '', myThid)

      ENDIF

#endif /* ALLOW_MNC */
#endif /* ALLOW_DARWIN */

      RETURN
      END
