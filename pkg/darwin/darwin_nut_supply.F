C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/darwin/darwin_nut_supply.F,v 1.1 2011/04/13 18:56:24 jahn Exp $
C $Name:  $

#include "DARWIN_OPTIONS.h"

CBOP
C !ROUTINE: DARWIN_NUT_SUPPLY
C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_NUT_SUPPLY(Nut,
     &                             bi,bj,imin,imax,jmin,jmax,
     &                             myIter,myTime, myThid )

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_FIELDS.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      _RL  Nut(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      INTEGER bi,bj,imin,imax,jmin,jmax
      INTEGER myThid, myIter
      _RL myTime

c!LOCAL VARIABLES: ====================================================
CEOP

#ifdef ALLOW_DARWIN

#ifdef ALLOW_QUOTA
#endif /* ALLOW_QUOTA */

#ifdef ALLOW_MONOD
c articficial supply of nutrients
      call MONOD_Nut_Supply( Nut,
     &                       bi,bj,imin,imax,jmin,jmax,
     &                       myTime,myIter,myThid)
#endif /* ALLOW_MONOD */

#endif /* ALLOW_DARWIN */

      RETURN
      END
