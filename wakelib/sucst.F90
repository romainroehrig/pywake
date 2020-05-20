SUBROUTINE SUCST(KULOUT,KPRINTLEV)

!**** *SUCST * - Routine to initialize the constants of the model.

!     Purpose.
!     --------
!           Initialize and print the common YOMCST + initialize
!         date and time of YOMRIP.

!**   Interface.
!     ----------
!        *CALL* *SUCST (..)

!        Explicit arguments :
!        --------------------


!        Implicit arguments :
!        --------------------
!        COMMON YOMCST

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!      R. Roehrig (from sucst.F90)
!      Original : 2019-05-02

!     Modifications.
!     --------------
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE YOMCST   , ONLY : RPI      ,RCLUM    ,RHPLA    ,RKBOL    ,&
 & RNAVO    ,RDAY     ,REA      ,REPSM    ,RSIYEA   ,&
 & RSIDAY   ,ROMEGA   ,RA       ,RG       ,R1SA     ,&
 & RSIGMA   ,RI0      ,R        ,RMD      ,RMV      ,&
 & RMO3     ,RD       ,RV       ,RCPD     ,RCPV     ,&
 & RMCO2    ,RMCH4    ,RMN2O    ,RMCO     ,RMHCHO   ,&
 & RMSO2    ,RMSO4    ,RMNO2    ,RMSF6    ,&
 & RCVD     ,RCVV     ,RKAPPA   ,RETV     ,RCW      ,&
 & RCS      ,RLVTT    ,RLSTT    ,RLVZER   ,RLSZER   ,&
 & RLMLT    ,RTT      ,RATM     ,RDT      ,RESTT    ,&
 & RALPW    ,RBETW    ,RGAMW    ,RALPS    ,RBETS    ,&
 & RGAMS    ,RALPD    ,RBETD    ,RGAMD    ,RMCFC11  ,&
 & RMCFC12  ,RMHCFC22 ,RMCCl4   ,REXC     ,RPERW    ,&
 & REOM     ,RDAYSOX  ,RUNDEF

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT
INTEGER(KIND=JPIM),INTENT(IN)    :: KPRINTLEV
INTEGER(KIND=JPIM) :: IA, ID, IDAT, IM, ISSS, J, IA2

REAL(KIND=JPRB) :: ZDE, ZET, ZJU, ZRS, ZRSREL, ZTETA, ZTI, ZTI2

LOGICAL :: LAQUA

!#include "fcttrm.h"
!#include "fcttim.h"

!      -----------------------------------------------------------------

LAQUA=.FALSE.

!*       1.    DEFINE FUNDAMENTAL CONSTANTS.
!              -----------------------------

RPI=2.0_JPRB*ASIN(1.0_JPRB)
RCLUM=299792458._JPRB
RHPLA=6.6260755E-34_JPRB
RKBOL=1.380658E-23_JPRB
RNAVO=6.0221367E+23_JPRB

!     ------------------------------------------------------------------

!*       2.    DEFINE ASTRONOMICAL CONSTANTS.
!              ------------------------------

RDAY=86400._JPRB
REA=149597870000._JPRB
REXC=0.016704_JPRB
! perihelie longitude (rad ~ 102.938°)
RPERW=1.79660702_JPRB
! Mean Earth orbital angular speed around the Sun(rad.yr-1)
REOM=6.283076_JPRB
! Day of the spring equinox (~ 21 march) ???? To be checked?
RDAYSOX=79.
!     ------------------------------------------------------------------

!*       3.    DEFINE RADIATION CONSTANTS.
!              ---------------------------

RSIGMA=2.0_JPRB * RPI**5 * RKBOL**4 /(15._JPRB* RCLUM**2 * RHPLA**3)
IF( LAQUA ) THEN
  ! aqua-planet special, obliquity == 0.0
  RI0=1365._JPRB
  REPSM=0.0_JPRB
ELSE
  RI0=1366._JPRB
  REPSM=0.409093_JPRB
ENDIF


RSIYEA=365.25_JPRB*RDAY*2.0_JPRB*RPI/REOM
RSIDAY=RDAY/(1.0_JPRB+RDAY/RSIYEA)
ROMEGA=2.0_JPRB*RPI/RSIDAY

!     ------------------------------------------------------------------

!*       4.    DEFINE GEOIDE.
!              --------------

RA=6371229._JPRB
IF( LAQUA ) THEN
  RG=9.79764_JPRB
! intercomparison RA=6371000.
ELSE
  RG=9.80665_JPRB
ENDIF
R1SA=REAL(1.0_JPRB/REAL(RA,KIND(1.0_JPRB)),KIND(R1SA))


!     ------------------------------------------------------------------

!*       5.    DEFINE THERMODYNAMIC CONSTANTS, GAS PHASE.
!              ------------------------------------------

R=RNAVO*RKBOL
RMD=28.9644_JPRB
RMV=18.0153_JPRB
RMO3=47.9942_JPRB
RD=1000._JPRB*R/RMD
RV=1000._JPRB*R/RMV
RCPD=3.5_JPRB*RD
RCVD=RCPD-RD
RCPV=4._JPRB *RV
RCVV=RCPV-RV
RKAPPA=RD/RCPD
RETV=RV/RD-1.0_JPRB
RMCO2=44.0095_JPRB
RMCH4=16.04_JPRB
RMN2O=44.013_JPRB
RMSF6=146.05_JPRB
RMCO=28.01_JPRB
RMHCHO=30.03_JPRB
RMNO2=46.01_JPRB
RMSO2=64.056_JPRB
RMSO4=96.052_JPRB
RMCFC11=137.37_JPRB
RMCFC12=120.91_JPRB
RMHCFC22=86.48_JPRB
RMCCl4=153.8336_JPRB
!     ------------------------------------------------------------------

!*       6.    DEFINE THERMODYNAMIC CONSTANTS, LIQUID PHASE.
!              ---------------------------------------------

RCW=4218._JPRB

!     ------------------------------------------------------------------

!*       7.    DEFINE THERMODYNAMIC CONSTANTS, SOLID PHASE.
!              --------------------------------------------

RCS=2106._JPRB

!     ------------------------------------------------------------------

!*       8.    DEFINE THERMODYNAMIC CONSTANTS, TRANSITION OF PHASE.
!              ----------------------------------------------------

RTT=273.16_JPRB
RDT=11.82_JPRB
RLVTT=2.5008E+6_JPRB
RLSTT=2.8345E+6_JPRB
RLVZER=RLVTT+RTT*(RCW-RCPV)
RLSZER=RLSTT+RTT*(RCS-RCPV)
RLMLT=RLSTT-RLVTT
RATM=100000._JPRB

!     ------------------------------------------------------------------

!*       9.    SATURATED VAPOUR PRESSURE.
!              --------------------------

RESTT=611.14_JPRB
RGAMW=(RCW-RCPV)/RV
RBETW=RLVTT/RV+RGAMW*RTT
RALPW=LOG(RESTT)+RBETW/RTT+RGAMW*LOG(RTT)
RGAMS=(RCS-RCPV)/RV
RBETS=RLSTT/RV+RGAMS*RTT
RALPS=LOG(RESTT)+RBETS/RTT+RGAMS*LOG(RTT)
RGAMD=RGAMS-RGAMW
RBETD=RBETS-RBETW
RALPD=RALPS-RALPW

! Missing value

RUNDEF = 1.E20_JPRB

!     ------------------------------------------------------------------

!*      10.    PRINTS

IF (KPRINTLEV >= 1) THEN
  WRITE(KULOUT,'(''0*** Constants of the ICM   ***'')')
  WRITE(KULOUT,'('' *** Fundamental constants ***'')')
  WRITE(KULOUT,'(''           PI = '',E13.7,'' -'')')RPI
  WRITE(KULOUT,'(''            c = '',E13.7,''m s-1'')')RCLUM
  WRITE(KULOUT,'(''            h = '',E13.7,''J s'')')RHPLA
  WRITE(KULOUT,'(''            K = '',E13.7,''J K-1'')')RKBOL
  WRITE(KULOUT,'(''            N = '',E13.7,''mol-1'')')RNAVO
  WRITE(KULOUT,'('' *** Astronomical constants ***'')')
  WRITE(KULOUT,'(''          day = '',E13.7,'' s'')')RDAY
  WRITE(KULOUT,'('' half g. axis = '',E13.7,'' m'')')REA
  WRITE(KULOUT,'('' mean anomaly = '',E13.7,'' rad (obliquity)'')')REPSM
  WRITE(KULOUT,'('' sideral year = '',E13.7,'' s'')')RSIYEA
  WRITE(KULOUT,'(''  sideral day = '',E13.7,'' s'')')RSIDAY
  WRITE(KULOUT,'(''        omega = '',E13.7,'' s-1'')')ROMEGA
  WRITE(KULOUT,'('' excentricity = '',E13.7,'' - '')')REXC
  WRITE(KULOUT,'(''    perihelie = '',E13.7,''  '')')RPERW

  WRITE(KULOUT,'('' The initial date of the run is :'')')
  WRITE(KULOUT,'(1X,I8,1X,I5,5X,I4,1X,I2,1X,I2)')IDAT,ISSS,IA,IM,ID
  WRITE(KULOUT,'('' The Julian date is : '',F11.2)') ZJU
  WRITE(KULOUT,'('' Time of the model  : '',F15.2,'' s'')')ZTI
  WRITE(KULOUT,'('' Distance Earth-Sun : '',E13.7,'' m'')')ZRS
  WRITE(KULOUT,'('' Relative Dist. E-S : '',E13.7,'' m'')')ZRSREL
  WRITE(KULOUT,'('' Declination        : '',F12.5)') ZDE
  WRITE(KULOUT,'('' Eq. of time        : '',F12.5,'' s'')')ZET
  WRITE(KULOUT,'('' ***         Geoide         ***'')')
  WRITE(KULOUT,'(''      Gravity = '',E13.7,'' m s-2'')')RG
  WRITE(KULOUT,'('' Earth radius = '',E13.7,'' m'')')RA
  WRITE(KULOUT,'('' Inverse E.R. = '',E13.7,'' m'')')R1SA
  WRITE(KULOUT,'('' ***        Radiation       ***'')')
  WRITE(KULOUT,'('' Stefan-Bol.  = '',E13.7,'' W m-2 K-4'')')  RSIGMA
  WRITE(KULOUT,'('' Solar const. = '',E13.7,'' W m-2'')')RI0
!  WRITE(KULOUT,'('' Correction of solar drift = '',L)') LASTRF
  WRITE(KULOUT,'('' *** Thermodynamic, gas     ***'')')
  WRITE(KULOUT,'('' Perfect gas  = '',e13.7)') R
  WRITE(KULOUT,'('' Dry air mass = '',e13.7)') RMD
  WRITE(KULOUT,'('' Vapour  mass = '',e13.7)') RMV
  WRITE(KULOUT,'('' Ozone   mass = '',e13.7)') RMO3
  WRITE(KULOUT,'('' Dry air cst. = '',e13.7)') RD
  WRITE(KULOUT,'('' Vapour  cst. = '',e13.7)') RV
  WRITE(KULOUT,'(''         Cpd  = '',e13.7)') RCPD
  WRITE(KULOUT,'(''         Cvd  = '',e13.7)') RCVD
  WRITE(KULOUT,'(''         Cpv  = '',e13.7)') RCPV
  WRITE(KULOUT,'(''         Cvv  = '',e13.7)') RCVV
  WRITE(KULOUT,'(''      Rd/Cpd  = '',e13.7)') RKAPPA
  WRITE(KULOUT,'(''     Rv/Rd-1  = '',e13.7)') RETV
  WRITE(KULOUT,'('' *** Thermodynamic, liquid  ***'')')
  WRITE(KULOUT,'(''         Cw   = '',E13.7)') RCW
  WRITE(KULOUT,'('' *** thermodynamic, solid   ***'')')
  WRITE(KULOUT,'(''         Cs   = '',E13.7)') RCS
  WRITE(KULOUT,'('' *** Thermodynamic, trans.  ***'')')
  WRITE(KULOUT,'('' Fusion point  = '',E13.7)') RTT
  WRITE(KULOUT,'('' RTT-Tx(ew-ei) = '',E13.7)') RDT
  WRITE(KULOUT,'(''        RLvTt  = '',E13.7)') RLVTT
  WRITE(KULOUT,'(''        RLsTt  = '',E13.7)') RLSTT
  WRITE(KULOUT,'(''        RLv0   = '',E13.7)') RLVZER
  WRITE(KULOUT,'(''        RLs0   = '',E13.7)') RLSZER
  WRITE(KULOUT,'(''        RLMlt  = '',E13.7)') RLMLT
  WRITE(KULOUT,'('' Normal press. = '',E13.7)') RATM
  WRITE(KULOUT,'('' Latent heat :  '')')
  WRITE(KULOUT,'(10(1X,E10.4))') (10._JPRB*J,J=-4,4)
!  WRITE(KULOUT,'(10(1X,E10.4))') (RLV(RTT+10._JPRB*J),J=-4,4)
!  WRITE(KULOUT,'(10(1X,E10.4))') (RLS(RTT+10._JPRB*J),J=-4,4)
  WRITE(KULOUT,'('' *** Thermodynamic, satur.  ***'')')
  WRITE(KULOUT,'('' Fusion point = '',E13.7)') RTT
  WRITE(KULOUT,'(''      es(Tt)  = '',e13.7)') RESTT
  WRITE(KULOUT,'('' es(T) :  '')')
  WRITE(KULOUT,'(10(1X,E10.4))') (10._JPRB*J,J=-4,4)
!  WRITE(KULOUT,'(10(1X,E10.4))') (ESW(RTT+10._JPRB*J),J=-4,4)
!  WRITE(KULOUT,'(10(1X,E10.4))') (ESS(RTT+10._JPRB*J),J=-4,4)
!  WRITE(KULOUT,'(10(1X,E10.4))') (ES (RTT+10._JPRB*J),J=-4,4)
  WRITE(KULOUT,'(''RUNDEF = '',E13.7)') RUNDEF
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE SUCST

