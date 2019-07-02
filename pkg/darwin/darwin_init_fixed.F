C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/darwin/darwin_init_fixed.F,v 1.1 2011/04/13 18:56:24 jahn Exp $
C $Name:  $

#include "DARWIN_OPTIONS.h"

C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_INIT_FIXED(myThid )

C !DESCRIPTION:
C calls subroutines that initialized fixed variables for any tracer
C experiment

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#ifdef ALLOW_SPECTRAL
#include "SPECTRAL_SIZE.h"
#include "SPECTRAL_PARAMS.h"
#endif

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid
CEOP

#ifdef ALLOW_DARWIN

C     !LOCAL VARIABLES:
#ifdef ALLOW_SPECTRAL
      INTEGER ilam
      _RL planck, c, hc, oavo, hcoavo, rlamm
#endif

#ifdef ALLOW_SPECTRAL
c   Quanta conversion
      planck = 6.6256 _d -34   !Plancks constant J sec
      c = 2.998 _d 8                 !speed of light m/sec
      hc = 1.0/(planck*c)
      oavo = 1.0/6.023 _d 23   ! 1/Avogadros number
      hcoavo = hc*oavo
      do ilam = 1,tlam
       rlamm = darwin_waves(ilam)*1 _d -9  !lambda in m
       WtouEins(ilam) = 1 _d 6*rlamm*hcoavo      !Watts to uEin/s conversion
      enddo
#endif

#ifdef ALLOW_QUOTA
      call QUOTA_INIT_FIXED(myThid)
#endif
#ifdef ALLOW_MONOD
      call MONOD_INIT_FIXED(myThid)
#endif

#endif /* ALLOW_DARWIN */

      RETURN
      END
