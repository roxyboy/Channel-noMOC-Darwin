C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/darwin/darwin_check.F,v 1.1 2011/04/13 18:56:24 jahn Exp $
C $Name:  $

#include "DARWIN_OPTIONS.h"

      SUBROUTINE DARWIN_CHECK( myThid )
C     /==========================================================\
C     | SUBROUTINE DARWIN_CHECK                                  |
C     | o Check dependences with other packages:                 |
C     |   Print darwin indices (iPO4,...) and ptracer names      |
C     |   to allow for manual check                              |
C     \==========================================================/
      IMPLICIT NONE

C     === Global variables ===
c#include "EEPARAMS.h"

C     === Routine arguments ===
C     myThid -  Number of this instances
      INTEGER myThid

C     === Local variables ===

#ifdef ALLOW_DARWIN
C     run checks specific to darwin
#ifdef ALLOW_QUOTA
      CALL QUOTA_CHECK( myThid )
#endif
#ifdef ALLOW_MONOD
      CALL MONOD_CHECK( myThid )
#endif
#endif /* ALLOW_DARWIN */
       
      RETURN
      END

