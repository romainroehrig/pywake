SUBROUTINE SUPHY2(KULOUT,PTSPHY)

!**** *SUPHY2 * - Routine to initialize physics parameters.

!     Purpose.
!     --------
!           Initialize and print the common YOMPHY2

!**   Interface.
!     ----------
!        *CALL* *SUPHY2 (..)

!        Explicit arguments :
!        --------------------


!        Implicit arguments :
!        --------------------
!        COMMON YOMPHY2

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!      R. Roehrig 
!      Original : 2019-05-02

!     Modifications.
!     --------------
!     ------------------------------------------------------------------

USE PARKIND1, ONLY : JPIM, JPRB

USE YOMPHY2,  ONLY : TSPHY

IMPLICIT NONE    

INTEGER(KIND=JPIM), INTENT(IN) :: KULOUT

REAL(KIND=JPRB), INTENT(IN)    :: PTSPHY

!     ------------------------------------------------------------------

TSPHY=PTSPHY

WRITE(KULOUT,*) '---------------------------------------------------'
WRITE(KULOUT,*) '          PHYSICS PARAMETERS'

WRITE(KULOUT,'('' TSPHY = '',E13.7,'' -'')')TSPHY

END SUBROUTINE SUPHY2 
