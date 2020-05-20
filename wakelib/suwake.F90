SUBROUTINE SUWAKE(KULOUT,                                                        &
              &   PWAPECUT,    PSIGMAD,     PHWMIN,    PSIGMAW_MAX,  PDENS_RATE, &
              &   PWDENSMIN,   PSTARK,      PALPK,     PWDENS_REF,   PCOEFGW,    &
              &   PTAU_CV,     PCREP_UPPER, PCREP_SOL, PDELTA_T_MIN, PRZERO)

!**** *SUWAKE * - Routine to initialize wake parameters.

!     Purpose.
!     --------
!           Initialize and print the common YOMWAKE

!**   Interface.
!     ----------
!        *CALL* *SUWAKE (..)

!        Explicit arguments :
!        --------------------


!        Implicit arguments :
!        --------------------
!        COMMON YOMWAKE

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

USE YOMWAKE,  ONLY : WAPECUT,    SIGMAD,     HWMIN,    SIGMAW_MAX,  DENS_RATE, &
                 &   WDENSMIN,   STARK,      ALPK,     WDENS_REF,   COEFGW,    &
                 &   TAU_CV,     CREP_UPPER, CREP_SOL, DELTA_T_MIN, RZERO


IMPLICIT NONE    

INTEGER(KIND=JPIM), INTENT(IN) :: KULOUT

REAL(KIND=JPRB), INTENT(IN)    :: PWAPECUT
REAL(KIND=JPRB), INTENT(IN)    :: PSIGMAD
REAL(KIND=JPRB), INTENT(IN)    :: PHWMIN
REAL(KIND=JPRB), INTENT(IN)    :: PSIGMAW_MAX
REAL(KIND=JPRB), INTENT(IN)    :: PDENS_RATE
REAL(KIND=JPRB), INTENT(IN)    :: PWDENSMIN
REAL(KIND=JPRB), INTENT(IN)    :: PSTARK
REAL(KIND=JPRB), INTENT(IN)    :: PALPK
REAL(KIND=JPRB), INTENT(IN)    :: PWDENS_REF(2)
REAL(KIND=JPRB), INTENT(IN)    :: PCOEFGW
REAL(KIND=JPRB), INTENT(IN)    :: PTAU_CV
REAL(KIND=JPRB), INTENT(IN)    :: PCREP_UPPER
REAL(KIND=JPRB), INTENT(IN)    :: PCREP_SOL
REAL(KIND=JPRB), INTENT(IN)    :: PDELTA_T_MIN
REAL(KIND=JPRB), INTENT(IN)    :: PRZERO

!     ------------------------------------------------------------------

WAPECUT      = PWAPECUT
SIGMAD       = PSIGMAD
HWMIN        = PHWMIN
SIGMAW_MAX   = PSIGMAW_MAX
DENS_RATE    = PDENS_RATE
WDENSMIN     = PWDENSMIN
STARK        = PSTARK
ALPK         = PALPK
WDENS_REF(1) = PWDENS_REF(1)
WDENS_REF(2) = PWDENS_REF(2)
COEFGW       = PCOEFGW
TAU_CV       = PTAU_CV
CREP_UPPER   = PCREP_UPPER
CREP_SOL     = PCREP_SOL
DELTA_T_MIN  = PDELTA_T_MIN
RZERO        = PRZERO

WRITE(KULOUT,*) '---------------------------------------------------'
WRITE(KULOUT,*) '          WAKE PARAMETERS'

WRITE(KULOUT,'('' WAPECUT     = '',E13.7,'' -'')')WAPECUT
WRITE(KULOUT,'('' SIGMAD      = '',E13.7,'' -'')')SIGMAD
WRITE(KULOUT,'('' HWMIN       = '',E13.7,'' -'')')HWMIN
WRITE(KULOUT,'('' SIGMAW_MAX  = '',E13.7,'' -'')')SIGMAW_MAX
WRITE(KULOUT,'('' DENS_RATE   = '',E13.7,'' -'')')DENS_RATE
WRITE(KULOUT,'('' WDENSMIN    = '',E13.7,'' -'')')WDENSMIN
WRITE(KULOUT,'('' STARK       = '',E13.7,'' -'')')STARK
WRITE(KULOUT,'('' ALPK        = '',E13.7,'' -'')')ALPK
WRITE(KULOUT,'('' WDENS_REF(1)= '',E13.7,'' -'')')WDENS_REF(1)
WRITE(KULOUT,'('' WDENS_REF(2)= '',E13.7,'' -'')')WDENS_REF(2)
WRITE(KULOUT,'('' COEFGW      = '',E13.7,'' -'')')COEFGW
WRITE(KULOUT,'('' TAU_CV      = '',E13.7,'' -'')')TAU_CV
WRITE(KULOUT,'('' CREP_UPPER  = '',E13.7,'' -'')')CREP_UPPER
WRITE(KULOUT,'('' CREP_SOL    = '',E13.7,'' -'')')CREP_SOL
WRITE(KULOUT,'('' DELTA_T_MIN = '',E13.7,'' -'')')DELTA_T_MIN
WRITE(KULOUT,'('' RZERO       = '',E13.7,'' -'')')RZERO

END SUBROUTINE SUWAKE      
