C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/darwin/darwin_forcing.F,v 1.1 2011/04/13 18:56:24 jahn Exp $
C $Name:  $

#include "DARWIN_OPTIONS.h"

CBOP
C !ROUTINE: DARWIN_FORCING
C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_FORCING(Ptr,
     &                          bi,bj,imin,imax,jmin,jmax,
     &                          myIter,myTime, myThid )

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_FIELDS.h"
#include "DARWIN_SIZE.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      _RL Ptr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,nDarwin)
      INTEGER bi,bj,imin,imax,jmin,jmax
      INTEGER myThid, myIter
      _RL myTime

c!LOCAL VARIABLES: ====================================================
CEOP

#ifdef ALLOW_DARWIN

#ifdef ALLOW_QUOTA
      call QUOTA_Forcing( Ptr,
     &                    bi,bj,imin,imax,jmin,jmax,
     &                    myTime,myIter,myThid)
#endif /* ALLOW_QUOTA */

#ifdef ALLOW_MONOD
C darwin_forcing operates on bi,bj part only, but needs to get full
C array because of last (iPtr) index
      call MONOD_Forcing(  Ptr,
     &                     bi,bj,imin,imax,jmin,jmax,
     &                     myTime,myIter,myThid)
#endif /* ALLOW_MONOD */

#endif /* ALLOW_DARWIN */

      RETURN
      END
