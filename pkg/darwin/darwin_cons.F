C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/darwin/darwin_cons.F,v 1.1 2011/04/13 18:56:24 jahn Exp $
C $Name:  $

#include "DARWIN_OPTIONS.h"

CBOP
C !ROUTINE: DARWIN_CONS
C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_CONS( myIter,myTime, myThid )

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_FIELDS.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER bi,bj,imin,imax,jmin,jmax
      INTEGER myThid, myIter
      _RL myTime

c!LOCAL VARIABLES: ====================================================
CEOP

#ifdef ALLOW_DARWIN

#ifdef ALLOW_QUOTA
      call QUOTA_cons(myTime,myIter,myThid)
#endif /* ALLOW_QUOTA */

#ifdef ALLOW_MONOD
      call MONOD_cons(myTime,myIter,myThid)
#endif /* ALLOW_MONOD */

#endif /* ALLOW_DARWIN */

      RETURN
      END
