C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/darwin/DARWIN_SIZE.h,v 1.1 2011/04/13 18:56:23 jahn Exp $
C $Name:  $

CBOP
C !ROUTINE: DARWIN_SIZE.h

C !INTERFACE:
C #include "DARWIN_SIZE.h"

C !DESCRIPTION:
C Includes _SIZE header of appropriate subpackage
C to make location of certain tracers within ptracers
C as well as number of phytoplankton, ... available to
C top-level darwin package
C
C !USES:

#ifdef ALLOW_MONOD
#include "MONOD_SIZE.h"
#endif

#ifdef ALLOW_QUOTA
#include "QUOTA_SIZE.h"
#endif

C !DEFINED PARAMETERS:

CBOP

