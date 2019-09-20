C $Header: /u/gcmpack/MITgcm_contrib/darwin2/verification/monod_1d_carbon/code/DARWIN_OPTIONS.h,v 1.1 2011/04/13 19:25:26 jahn Exp $
C $Name:  $

#ifndef DARWIN_OPTIONS_H
#define DARWIN_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#ifdef ALLOW_DARWIN

#include "CPP_OPTIONS.h"

CBOP
C    !ROUTINE: DARWIN_OPTIONS.h
C    !INTERFACE:

C    !DESCRIPTION:
c options for darwin package
CEOP

#define READ_PAR
#define MINFE
#undef  NUT_SUPPLY
#undef  FLUX_NUTS
#define RELAX_NUTS
#undef  CONS_SUPP
#undef  OLD_GRAZE
#undef  ALLOW_DIAZ
#undef  OLD_NSCHEME
#undef  ALLOW_MUTANTS
#define PORT_RAND
#undef OLDSEED

#undef NOTEMP
#define TEMP_VERSION 2
#define TEMP_RANGE

#define TWO_SPECIES_SETUP
#undef NINE_SPECIES_SETUP

#define CALC_RATE_TOTALS

#undef IRON_SED_SOURCE
#undef IRON_SED_SOURCE_VARIABLE
#define PART_SCAV

#define ALLOW_DIAZ
#define ALLOW_DENIT
#undef  DENIT_RELAX

#define ALLOW_CARBON
#define ALLOW_OLD_VIRTUALFLUX


#undef  GEIDER

c ANNA turn wavebands on/off
#undef WAVEBANDS 
c#define WAVEBANDS

#undef  CHECK_CONS
#undef  DAR_DIAG_RSTAR
#undef  DAR_DIAG_DIVER
#define DAR_DIAG_IRON

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C dependencies
c if two or nine species setup we don't want specific temperature ranges
#ifdef  TWO_SPECIES_SETUP
#undef TEMP_RANGE
#endif
#ifdef  NINE_SPECIES_SETUP
#undef TEMP_RANGE
#endif

#ifdef DAR_DIAG_CHL
#define ALLOW_PAR_DAY
#endif
 
#endif /* ALLOW_DARWIN */
#endif /* DARWIN_OPTIONS_H */

