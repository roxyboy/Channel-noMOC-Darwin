#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "DARWIN_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_DARWIN

#ifdef ALLOW_CARBON


      SUBROUTINE DIC_READ_PICKUP(
     O                            pH_isLoaded,
     I                            myIter, myThid )

      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DARWIN_FLUX.h"

C     == Routine arguments ==
C     myThid  :: my Thread Id number
      LOGICAL pH_isLoaded
      INTEGER myIter
      INTEGER myThid


C     !FUNCTIONS

C     !LOCAL VARIABLES:
C     == Local variables ==
      CHARACTER*(MAX_LEN_FNAM) fn, filNam
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      LOGICAL useCurrentDir, fileExist
      INTEGER fp, ioUnit
CEOP

      pH_isLoaded =.FALSE.
      ioUnit = errorMessageUnit

C--   pickup file name :
      IF (pickupSuff.EQ.' ') THEN
        WRITE(fn,'(A,I10.10)') 'pickup_dic.', myIter
      ELSE
        WRITE(fn,'(A,A10)') 'pickup_dic.', pickupSuff
      ENDIF
      fp = precFloat64

C--   First check if pickup file exist
#ifdef ALLOW_MDSIO
      useCurrentDir = .FALSE.
      CALL MDS_CHECK4FILE(
     I                     fn, '.data', 'DIC_READ_PICKUP',
     O                     filNam, fileExist,
     I                     useCurrentDir, myThid )
#else
      STOP 'ABNORMAL END: S/R DIC_READ_PICKUP: Needs MDSIO pkg'
#endif

      IF ( fileExist ) THEN
C--   Read pickup file
        CALL READ_REC_3D_RL( fn, fp, 1, pH, 1, myIter, myThid )
        pH_isLoaded = .TRUE.

        _EXCH_XY_RL( pH, myThid )
      ELSE
        pH_isLoaded = .FALSE.
        IF ( pickupStrictlyMatch ) THEN
          WRITE(msgBuf,'(4A)') 'DIC_READ_PICKUP: ',
     &      'try with " pickupStrictlyMatch=.FALSE.,"',
     &      ' in file: "data", NameList: "PARM03"'
          CALL PRINT_MESSAGE( msgBuf, ioUnit, SQUEEZE_RIGHT, myThid )
          STOP 'ABNORMAL END: S/R DIC_READ_PICKUP'
        ELSE
          WRITE(msgBuf,'(2A)') 'WARNING >> DIC_READ_PICKUP: ',
     &      'will restart from approximated pH'
          CALL PRINT_MESSAGE( msgBuf, ioUnit, SQUEEZE_RIGHT, myThid )
        ENDIF
      ENDIF


      RETURN
      END
#endif  /*ALLOW_CARBON*/

#endif  /*DARWIN*/
#endif  /*ALLOW_PTRACERS*/
c ==================================================================
