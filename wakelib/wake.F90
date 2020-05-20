!OPTIONS XOPT(NOEVAL)
SUBROUTINE wake(KIDIA, KFDIA, KLON,KLEV,znatsurf, p, ph, pi, dtime, &
                te0, qe0, omgb, &
                dtdwn, dqdwn, amdwn, amup, dta, dqa, wgen, &
                sigd_con, Cin, &
                deltatw, deltaqw, sigmaw, awdens, wdens, &                  ! state variables
                dth, hw, wape, fip, gfl, &
                dtls, dqls, ktopw, omgbdth, dp_omgb, tu, qu, &
                dtke, dqke, omg, dp_deltomg, spread, cstar, &
                d_deltat_gw, &
                d_deltatw2, d_deltaqw2, d_sigmaw2, d_awdens2, d_wdens2)     ! tendencies


  ! **************************************************************
  ! *
  ! WAKE                                                        *
  ! retour a un Pupper fixe                                *
  ! *
  ! written by   :  GRANDPEIX Jean-Yves   09/03/2000            *
  ! modified by :   ROEHRIG Romain        01/29/2007            *
  ! **************************************************************

  
USE PARKIND1, ONLY : JPIM, JPRB
!USE YOMHOOK,  ONLY : LHOOK,   DR_HOOK

USE YOMCST, ONLY : RG, RD, RETV

USE YOMWAKE, ONLY : WAPECUT, SIGMAD, HWMIN, SIGMAW_MAX, DENS_RATE, WDENSMIN, &
                 &  STARK, ALPK, WDENS_REF, COEFGW, TAU_CV, CREP_UPPER,CREP_SOL, DELTA_T_MIN, RZERO

USE YOMPHY2  , ONLY : TSPHY
  
IMPLICIT NONE
  ! ============================================================================


  ! But : Decrire le comportement des poches froides apparaissant dans les
  ! grands systemes convectifs, et fournir l'energie disponible pour
  ! le declenchement de nouvelles colonnes convectives.

  ! State variables : 
  ! deltatw    : temperature difference between wake and off-wake regions
  ! deltaqw    : specific humidity difference between wake and off-wake regions
  ! sigmaw     : fractional area covered by wakes.
  ! wdens      : number of wakes per unit area

  ! Variable de sortie :

  ! wape : WAke Potential Energy
  ! fip  : Front Incident Power (W/m2) - ALP
  ! gfl  : Gust Front Length per unit area (m-1)
  ! dtls : large scale temperature tendency due to wake
  ! dqls : large scale humidity tendency due to wake
  ! hw   : wake top hight (given by hw*deltatw(1)/2=wape)
  ! dp_omgb : vertical gradient of large scale omega
  ! awdens  : densite de poches actives
  ! wdens   : densite de poches
  ! omgbdth: flux of Delta_Theta transported by LS omega
  ! dtKE   : differential heating (wake - unpertubed)
  ! dqKE   : differential moistening (wake - unpertubed)
  ! omg    : Delta_omg =vertical velocity diff. wake-undist. (Pa/s)
  ! dp_deltomg  : vertical gradient of omg (s-1)
  ! spread  : spreading term in d_t_wake and d_q_wake
  ! deltatw     : updated temperature difference (T_w-T_u).
  ! deltaqw     : updated humidity difference (q_w-q_u).
  ! sigmaw      : updated wake fractional area.
  ! d_deltat_gw : delta T tendency due to GW

  ! Variables d'entree :

  ! aire : aire de la maille
  ! te0  : temperature dans l'environnement  (K)
  ! qe0  : humidite dans l'environnement     (kg/kg)
  ! omgb : vitesse verticale moyenne sur la maille (Pa/s)
  ! dtdwn: source de chaleur due aux descentes (K/s)
  ! dqdwn: source d'humidite due aux descentes (kg/kg/s)
  ! dta  : source de chaleur due courants satures et detrain  (K/s)
  ! dqa  : source d'humidite due aux courants satures et detra (kg/kg/s)
  ! wgen : number of wakes generated per unit area and per sec (/m^2/s)
  ! amdwn: flux de masse total des descentes, par unite de
  !        surface de la maille (kg/m2/s)
  ! amup : flux de masse total des ascendances, par unite de
  !        surface de la maille (kg/m2/s)
  ! sigd_con: 
  ! Cin  : convective inhibition
  ! p    : pressions aux milieux des couches (Pa)
  ! ph   : pressions aux interfaces (Pa)
  ! pi  : (p/p_0)**kapa (adim)
  ! dtime: increment temporel (s)

  ! Variables internes :

  ! rhow : masse volumique de la poche froide
  ! rho  : environment density at P levels
  ! rhoh : environment density at Ph levels
  ! te   : environment temperature | may change within
  ! qe   : environment humidity    | sub-time-stepping
  ! the  : environment potential temperature
  ! thu  : potential temperature in undisturbed area
  ! tu   :  temperature  in undisturbed area
  ! qu   : humidity in undisturbed area
  ! dp_omgb: vertical gradient og LS omega
  ! omgbw  : wake average vertical omega
  ! dp_omgbw: vertical gradient of omgbw
  ! omgbdq : flux of Delta_q transported by LS omega
  ! dth  : potential temperature diff. wake-undist.
  ! th1  : first pot. temp. for vertical advection (=thu)
  ! th2  : second pot. temp. for vertical advection (=thw)
  ! q1   : first humidity for vertical advection
  ! q2   : second humidity for vertical advection
  ! d_deltatw   : terme de redistribution pour deltatw
  ! d_deltaqw   : terme de redistribution pour deltaqw
  ! deltatw0   : deltatw initial
  ! deltaqw0   : deltaqw initial
  ! hw0    : wake top hight (defined as the altitude at which deltatw=0)
  ! amflux : horizontal mass flux through wake boundary
  ! wdens_ref: initial number of wakes per unit area (3D) or per
  ! unit length (2D), at the beginning of each time step
  ! Tgw    : 1 sur la période de onde de gravité
  ! Cgw    : vitesse de propagation de onde de gravité
  ! LL     : distance entre 2 poches

  ! -------------------------------------------------------------------------
  ! Déclaration de variables
  ! -------------------------------------------------------------------------

!  include "YOMCST.h"
!  include "cvthermo.h"

  ! Arguments en entree
  ! --------------------

  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  INTEGER(KIND=JPIM), INTENT(IN) :: KLON

  INTEGER(KIND=JPIM), INTENT(IN) :: znatsurf(klon)
  REAL(KIND=JPRB), INTENT(IN)    :: p(klon, klev)
  REAL(KIND=JPRB), INTENT(IN)    :: pi(klon, klev)
  REAL(KIND=JPRB), INTENT(IN)    :: ph(klon, klev+1)
  REAL(KIND=JPRB), INTENT(IN)    :: omgb(klon, klev)
  REAL(KIND=JPRB), INTENT(IN)    :: dtime
  REAL(KIND=JPRB), INTENT(IN)    :: te0(klon, klev)
  REAL(KIND=JPRB), INTENT(IN)    :: qe0(klon, klev)
  REAL(KIND=JPRB), INTENT(IN)    :: dtdwn(klon, klev)
  REAL(KIND=JPRB), INTENT(IN)    :: dqdwn(klon, klev)
  REAL(KIND=JPRB), INTENT(IN)    :: amdwn(klon, klev)
  REAL(KIND=JPRB), INTENT(IN)    :: amup(klon, klev)
  REAL(KIND=JPRB), INTENT(IN)    :: dta(klon, klev)
  REAL(KIND=JPRB), INTENT(IN)    :: dqa(klon, klev)
  REAL(KIND=JPRB), INTENT(IN)    :: wgen(klon)
  REAL(KIND=JPRB), INTENT(IN)    :: sigd_con(klon)
  REAL(KIND=JPRB), INTENT(IN)    :: Cin(klon)

  ! Entrees/Sorties
  ! --------

  REAL(KIND=JPRB), INTENT(INOUT)    :: deltatw(klon, klev)
  REAL(KIND=JPRB), INTENT(INOUT)    :: deltaqw(klon, klev)
  REAL(KIND=JPRB), INTENT(INOUT)    :: sigmaw(klon)
  REAL(KIND=JPRB), INTENT(INOUT)    :: awdens(klon)
  REAL(KIND=JPRB), INTENT(INOUT)    :: wdens(klon)

  ! Sorties
  ! --------

  REAL(KIND=JPRB), INTENT(OUT)    :: dth(klon, klev)
  REAL(KIND=JPRB), INTENT(OUT)    :: tu(klon, klev)
  REAL(KIND=JPRB), INTENT(OUT)    :: qu(klon, klev)
  REAL(KIND=JPRB), INTENT(OUT)    :: dtls(klon, klev)
  REAL(KIND=JPRB), INTENT(OUT)    :: dqls(klon, klev)
  REAL(KIND=JPRB), INTENT(OUT)    :: dtke(klon, klev)
  REAL(KIND=JPRB), INTENT(OUT)    :: dqke(klon, klev)
  REAL(KIND=JPRB), INTENT(OUT)    :: spread(klon, klev)
  REAL(KIND=JPRB), INTENT(OUT)    :: omgbdth(klon, klev)
  REAL(KIND=JPRB), INTENT(OUT)    :: omg(klon, klev)
  REAL(KIND=JPRB), INTENT(OUT)    :: dp_omgb(klon, klev)
  REAL(KIND=JPRB), INTENT(OUT)    :: dp_deltomg(klon, klev)
  REAL(KIND=JPRB), INTENT(OUT)    :: d_deltat_gw(klon, klev)  
  REAL(KIND=JPRB), INTENT(OUT)    :: hw(klon)
  REAL(KIND=JPRB), INTENT(OUT)    :: wape(klon)
  REAL(KIND=JPRB), INTENT(OUT)    :: fip(klon)
  REAL(KIND=JPRB), INTENT(OUT)    :: gfl(klon)
  REAL(KIND=JPRB), INTENT(OUT)    :: cstar(klon)
  INTEGER(KIND=JPIM), INTENT(OUT) :: ktopw(klon)
  REAL(KIND=JPRB), INTENT(OUT)    :: d_deltatw2(klon, klev)
  REAL(KIND=JPRB), INTENT(OUT)    :: d_deltaqw2(klon, klev)
  REAL(KIND=JPRB), INTENT(OUT)    :: d_sigmaw2(klon)
  REAL(KIND=JPRB), INTENT(OUT)    :: d_awdens2(klon)
  REAL(KIND=JPRB), INTENT(OUT)    :: d_wdens2(klon)

  ! Variables internes
  ! -------------------

  ! Variables à fixer
  !REAL(KIND=JPRB) :: alon
  INTEGER(KIND=JPIM) :: igout
!  LOGICAL, SAVE :: first = .TRUE.
!  REAL, SAVE ::  stark, wdens_ref, coefgw, alpk, crep_upper, crep_sol  
!  REAL(KIND=JPRB)    :: crep_upper, crep_sol, wdens_ref(2)
!  REAL(KIND=JPRB)    :: tau_cv
!  REAL(KIND=JPRB)    :: rzero ! minimal wake radius
  REAL(KIND=JPRB)    :: aa0 ! minimal wake area
  REAL(KIND=JPRB)    :: cstart
  LOGICAL            :: flag_wk_check_trgl
  INTEGER(KIND=JPIM) :: iflag_wk_check_trgl
  INTEGER(KIND=JPIM) :: iflag_wk_pop_dyn
!  REAL(KIND=JPRB)    :: delta_t_min
  INTEGER(KIND=JPIM) :: nsub
  REAL(KIND=JPRB)    :: dtimesub
!  REAL(KIND=JPRB)    :: wdensmin
!  REAL(KIND=JPRB) :: sigmad, hwmin, wapecut, cstart
!  REAL(KIND=JPRB)    :: sigmaw_max
!  REAL(KIND=JPRB)    :: dens_rate
  REAL(KIND=JPRB)    :: wdens0
  ! IM 080208
  LOGICAL            :: gwake(KLON)

  ! Variables de sauvegarde
  REAL(KIND=JPRB)    :: deltatw0(klon, klev)
  REAL(KIND=JPRB)    :: deltaqw0(klon, klev)
  REAL(KIND=JPRB)    :: te(klon, klev), qe(klon, klev)
!!  REAL(KIND=JPRB) :: sigmaw0(klon), sigmaw1(klon)

  ! Variables liees a la dynamique de population
  REAL(KIND=JPRB)    :: act(klon)
  REAL(KIND=JPRB)    :: rad_wk(klon), tau_wk_inv(klon)
  REAL(KIND=JPRB)    :: f_shear(klon)
  REAL(KIND=JPRB)    :: drdt(klon)
  REAL(KIND=JPRB)    :: d_sig_gen(klon), d_sig_death(klon), d_sig_col(klon)
  REAL(KIND=JPRB)    :: wape1_act(klon), wape2_act(klon)
  LOGICAL            :: kill_wake(klon)
  INTEGER(KIND=JPIM) :: iflag_wk_act

  REAL(KIND=JPRB)    :: drdt_pos
  REAL(KIND=JPRB)    :: tau_wk_inv_min

  ! Variables pour les GW
  REAL(KIND=JPRB)    :: ll(klon)
  REAL(KIND=JPRB)    :: n2(klon, klev)
  REAL(KIND=JPRB)    :: cgw(klon, klev)
  REAL(KIND=JPRB)    :: tgw(klon, klev)

  ! Variables liées au calcul de hw
  REAL(KIND=JPRB)    :: ptop_provis(klon), ptop(klon), ptop_new(klon)
  REAL(KIND=JPRB)    :: sum_dth(klon)
  REAL(KIND=JPRB)    :: dthmin(klon)
  REAL(KIND=JPRB)    :: z(klon), dz(klon), hw0(klon)
  INTEGER(KIND=JPIM) :: ktop(klon), kupper(klon)

  ! Variables liees au test de la forme triangulaire du profil de Delta_theta
  REAL(KIND=JPRB)    :: sum_half_dth(klon)
  REAL(KIND=JPRB)    :: dz_half(klon)

  ! Sub-timestep tendencies and related variables
  REAL(KIND=JPRB)    :: d_deltatw(klon, klev), d_deltaqw(klon, klev)
  REAL(KIND=JPRB)    :: d_te(klon, klev), d_qe(klon, klev)
  REAL(KIND=JPRB)    :: d_awdens(klon), d_wdens(klon)
  REAL(KIND=JPRB)    :: d_sigmaw(klon), alpha(klon)
  REAL(KIND=JPRB)    :: q0_min(klon), q1_min(klon)
  LOGICAL            :: wk_adv(klon), ok_qx_qw(klon)
  REAL(KIND=JPRB)    :: epsilon
!  DATA epsilon/1.E-15/
  REAL(KIND=JPRB) :: epsim1

  ! Autres variables internes
  INTEGER(KIND=JPIM) :: isubstep, k, i
  INTEGER(KIND=JPIM) :: prt_level

  REAL(KIND=JPRB)    :: wdens_targ
  REAL(KIND=JPRB)    :: sigmaw_targ

  REAL(KIND=JPRB)    :: sum_thu(klon), sum_tu(klon), sum_qu(klon), sum_thvu(klon)
  REAL(KIND=JPRB)    :: sum_dq(klon), sum_rho(klon)
  REAL(KIND=JPRB)    :: sum_dtdwn(klon), sum_dqdwn(klon)
  REAL(KIND=JPRB)    :: av_thu(klon), av_tu(klon), av_qu(klon), av_thvu(klon)
  REAL(KIND=JPRB)    :: av_dth(klon), av_dq(klon), av_rho(klon)
  REAL(KIND=JPRB)    :: av_dtdwn(klon), av_dqdwn(klon)

  REAL(KIND=JPRB)    :: rho(klon, klev), rhow(klon, klev)
  REAL(KIND=JPRB)    :: rhoh(klon, klev+1) 
  REAL(KIND=JPRB)    :: rhow_moyen(klon, klev)
  REAL(KIND=JPRB)    :: zh(klon, klev)
  REAL(KIND=JPRB)    :: zhh(klon, klev+1) 
  REAL(KIND=JPRB)    :: epaisseur1(klon, klev), epaisseur2(klon, klev)

  REAL(KIND=JPRB)    :: the(klon, klev), thu(klon, klev)

  REAL(KIND=JPRB)    :: omgbw(klon, klev)
  REAL(KIND=JPRB)    :: pupper(klon)
  REAL(KIND=JPRB)    :: omgtop(klon)
  REAL(KIND=JPRB)    :: dp_omgbw(klon, klev)
  REAL(KIND=JPRB)    :: ztop(klon), dztop(klon)
  REAL(KIND=JPRB)    :: alpha_up(klon, klev)

  REAL(KIND=JPRB)    :: rre1(klon), rre2(klon)
  REAL(KIND=JPRB)    :: rrd1, rrd2
  REAL(KIND=JPRB)    :: th1(klon, klev), th2(klon, klev), q1(klon, klev), q2(klon, klev)
  REAL(KIND=JPRB)    :: d_th1(klon, klev), d_th2(klon, klev), d_dth(klon, klev)
  REAL(KIND=JPRB)    :: d_q1(klon, klev), d_q2(klon, klev), d_dq(klon, klev)
  REAL(KIND=JPRB)    :: omgbdq(klon, klev)

  REAL(KIND=JPRB)    :: ff(klon), gg(klon)
  REAL(KIND=JPRB)    :: wape2(klon), cstar2(klon), heff(klon)

  REAL(KIND=JPRB)    :: crep(klon, klev)

  REAL(KIND=JPRB)    :: ppi(klon, klev)

  ! cc nrlmd
  REAL(KIND=JPRB)    :: death_rate(klon)
!  REAL(KIND=JPRB) :: nat_rate(klon)
  REAL(KIND=JPRB)    :: entr(klon, klev)
  REAL(KIND=JPRB)    :: detr(klon, klev)

  REAL(KIND=JPRB)    :: sigmaw_in(klon)             ! pour les prints
  REAL(KIND=JPRB)    :: awdens_in(klon), wdens_in(klon)   ! pour les prints

  REAL(KIND=JPRB) :: ZHOOK_HANDLE

  ! -------------------------------------------------------------------------
  ! Initialisations
  ! -------------------------------------------------------------------------

!  IF (LHOOK) CALL DR_HOOK('WAKE',0,ZHOOK_HANDLE)

  prt_level = 0

  IF (prt_level >= 10) THEN
    print*, 'wake initialisations'
  ENDIF
  epsilon = 1.E-15_JPRB

  epsim1 = RETV

  ! Essais d'initialisation avec sigmaw = 0.02 et hw = 10.
  ! -------------------------------------------------------------------------

!!  DATA wapecut, sigmad, hwmin/5., .02, 10./
!  DATA wapecut, sigmad, hwmin/1., .02, 10./
!!  DATA wdensmin/1.e-12/
!  DATA wdensmin/1.e-14/
!  ! cc nrlmd
!  DATA sigmaw_max/0.4/
!  DATA dens_rate/0.1/
!  DATA rzero /5000./
  ! cc
  ! Longueur de maille (en m)
  ! -------------------------------------------------------------------------

  ! ALON = 3.e5
  ! alon = 1.E6

  !  Provisionnal; to be suppressed when f_shear is parameterized
  f_shear(:) = 1.       ! 0. for strong shear, 1. for weak shear


  ! Configuration de coefgw,stark,wdens (22/02/06 by YU Jingmei)

  ! coefgw : Coefficient pour les ondes de gravité
  ! stark : Coefficient k dans Cstar=k*sqrt(2*WAPE)
  ! wdens : Densité surfacique de poche froide
  ! -------------------------------------------------------------------------

  ! cc nrlmd      coefgw=10
  ! coefgw=1
  ! wdens0 = 1.0/(alon**2)
  ! cc nrlmd      wdens = 1.0/(alon**2)
  ! cc nrlmd      stark = 0.50
  ! CRtest
  ! cc nrlmd      alpk=0.1
  ! alpk = 1.0
  ! alpk = 0.5
  ! alpk = 0.05

! if (first) then

  igout = klon/2+1/klon

!  crep_upper = 0.9
!  crep_sol = 1.0

  aa0 = 3.14*rzero*rzero

!  tau_cv = 4000.

  ! cc nrlmd Lecture du fichier wake_param.data
!  stark=0.33
!  CALL getin_p('stark',stark)
  cstart = stark*sqrt(2.*wapecut)

!  alpk=0.25
!  CALL getin_p('alpk',alpk)
!jyg<
!!  wdens_ref=8.E-12
!!  CALL getin_p('wdens_ref',wdens_ref)
!  wdens_ref(1)=8.E-12
!  wdens_ref(2)=8.E-12
!  CALL getin_p('wdens_ref_o',wdens_ref(1))    !wake number per unit area ; ocean
!  CALL getin_p('wdens_ref_l',wdens_ref(2))    !wake number per unit area ; land
!>jyg
  iflag_wk_pop_dyn = 0
!  CALL getin_p('iflag_wk_pop_dyn',iflag_wk_pop_dyn) ! switch between wdens prescribed 
                                                    ! and wdens prognostic
  iflag_wk_act = 0
!  CALL getin_p('iflag_wk_act',iflag_wk_act) ! 0: act(:)=0.
                                            ! 1: act(:)=1.
                                            ! 2: act(:)=f(Wape)
!  coefgw=4.
!  CALL getin_p('coefgw',coefgw)

!  WRITE(*,*) 'stark=', stark
!  WRITE(*,*) 'alpk=', alpk
!jyg<
!!  WRITE(*,*) 'wdens_ref=', wdens_ref
!  WRITE(*,*) 'wdens_ref_o=', wdens_ref(1)
!  WRITE(*,*) 'wdens_ref_l=', wdens_ref(2)
!>jyg
!  WRITE(*,*) 'iflag_wk_pop_dyn=',iflag_wk_pop_dyn
!  WRITE(*,*) 'iflag_wk_act',iflag_wk_act
!  WRITE(*,*) 'coefgw=', coefgw 

  flag_wk_check_trgl=.false.
!  CALL getin_p('flag_wk_check_trgl ', flag_wk_check_trgl)
!  WRITE(*,*) 'flag_wk_check_trgl=', flag_wk_check_trgl
!  WRITE(*,*) 'flag_wk_check_trgl OBSOLETE. Utilisr iflag_wk_check_trgl plutot'
  iflag_wk_check_trgl=0 ; IF (flag_wk_check_trgl) iflag_wk_check_trgl=1
!  CALL getin_p('iflag_wk_check_trgl ', iflag_wk_check_trgl)
!  WRITE(*,*) 'iflag_wk_check_trgl=', iflag_wk_check_trgl
!
!  first=.false.
! endif

 IF (iflag_wk_pop_dyn == 0) THEN
  ! Initialisation de toutes des densites a wdens_ref.
  ! Les densites peuvent evoluer si les poches debordent
  ! (voir au tout debut de la boucle sur les substeps)
  !jyg<
  !!  wdens(:) = wdens_ref
   DO i = KIDIA,KFDIA
     wdens(i) = wdens_ref(znatsurf(i)+1)
   ENDDO 
  !>jyg
 ENDIF  ! (iflag_wk_pop_dyn == 0)

  ! print*,'stark',stark
  ! print*,'alpk',alpk
  ! print*,'wdens',wdens
  ! print*,'coefgw',coefgw
  ! cc
  ! Minimum value for |T_wake - T_undist|. Used for wake top definition
  ! -------------------------------------------------------------------------

!  delta_t_min = 0.2

  ! 1. - Save initial values, initialize tendencies, initialize output fields
  ! ------------------------------------------------------------------------

!jyg<
!!  DO k = 1, klev
!!    DO i = 1, klon
!!      ppi(i, k) = pi(i, k)
!!      deltatw0(i, k) = deltatw(i, k)
!!      deltaqw0(i, k) = deltaqw(i, k)
!!      te(i, k) = te0(i, k)
!!      qe(i, k) = qe0(i, k)
!!      dtls(i, k) = 0.
!!      dqls(i, k) = 0.
!!      d_deltat_gw(i, k) = 0.
!!      d_te(i, k) = 0.
!!      d_qe(i, k) = 0.
!!      d_deltatw(i, k) = 0.
!!      d_deltaqw(i, k) = 0.
!!      ! IM 060508 beg
!!      d_deltatw2(i, k) = 0.
!!      d_deltaqw2(i, k) = 0.
!!      ! IM 060508 end
!!    END DO
!!  END DO
      ppi(:,:) = pi(:,:)
      deltatw0(:,:) = deltatw(:,:)
      deltaqw0(:,:) = deltaqw(:,:)
      te(:,:) = te0(:,:)
      qe(:,:) = qe0(:,:)
      dtls(:,:) = 0.
      dqls(:,:) = 0.
      d_deltat_gw(:,:) = 0.
      d_te(:,:) = 0.
      d_qe(:,:) = 0.
      d_deltatw(:,:) = 0.
      d_deltaqw(:,:) = 0.
      d_deltatw2(:,:) = 0.
      d_deltaqw2(:,:) = 0.

      IF (iflag_wk_act == 0) THEN
        act(:) = 0.
      ELSEIF (iflag_wk_act == 1) THEN
        act(:) = 1.
      ENDIF

!!  DO i = 1, klon
!!   sigmaw_in(i) = sigmaw(i)
!!  END DO
   sigmaw_in(:) = sigmaw(:)
!>jyg

  ! sigmaw1=sigmaw
  ! IF (sigd_con.GT.sigmaw1) THEN
  ! print*, 'sigmaw,sigd_con', sigmaw, sigd_con
  ! ENDIF
  IF (iflag_wk_pop_dyn >=1) THEN
    DO i = KIDIA, KFDIA
      wdens_targ = max(wdens(i),wdensmin)
      d_wdens2(i) = wdens_targ - wdens(i)
      wdens(i) = wdens_targ
    END DO
  ELSE
    DO i = KIDIA, KFDIA
      d_awdens2(i) = 0.
      d_wdens2(i) = 0.
    END DO
  ENDIF  ! (iflag_wk_pop_dyn >=1)
!
  DO i = KIDIA, KFDIA
    ! c      sigmaw(i) = amax1(sigmaw(i),sigd_con(i))
!jyg<
!!    sigmaw(i) = amax1(sigmaw(i), sigmad)
!!    sigmaw(i) = amin1(sigmaw(i), 0.99)
    sigmaw_targ = min(max(sigmaw(i), sigmad),0.99)
    d_sigmaw2(i) = sigmaw_targ - sigmaw(i)
    sigmaw(i) = sigmaw_targ
!>jyg
  END DO

!
  IF (iflag_wk_pop_dyn >= 1) THEN
    awdens_in(:) = awdens(:)
    wdens_in(:) = wdens(:)
!!    wdens(:) = wdens(:) + wgen(:)*dtime
!!    d_wdens2(:) = wgen(:)*dtime
!!  ELSE
  ENDIF  ! (iflag_wk_pop_dyn >= 1)

  wape(:) = 0.
  wape2(:) = 0.
  d_sigmaw(:) = 0.
  ktopw(:) = 0
!
!<jyg
dth(:,:) = 0.
tu(:,:) = 0.
qu(:,:) = 0.
dtke(:,:) = 0.
dqke(:,:) = 0.
spread(:,:) = 0.
omgbdth(:,:) = 0.
omg(:,:) = 0.
dp_omgb(:,:) = 0.
dp_deltomg(:,:) = 0.
hw(:) = 0.
wape(:) = 0.
fip(:) = 0.
gfl(:) = 0.
cstar(:) = 0.
ktopw(:) = 0
!
!  Vertical advection local variables
omgbw(:,:) = 0.
omgtop(:) = 0
dp_omgbw(:,:) = 0.
omgbdq(:,:) = 0.
!>jyg
!
  IF (prt_level>=10) THEN
    PRINT *, 'wake-1, sigmaw(igout) ', sigmaw(igout)
    PRINT *, 'wake-1, deltatw(igout,k) ', (k,deltatw(igout,k), k=1,klev)
    PRINT *, 'wake-1, deltaqw(igout,k) ', (k,deltaqw(igout,k), k=1,klev)
    PRINT *, 'wake-1, deltatw0(igout,k) ', (k,deltatw0(igout,k), k=1,klev)
    PRINT *, 'wake-1, deltaqw0(igout,k) ', (k,deltaqw0(igout,k), k=1,klev)    
    PRINT *, 'wake-1, dowwdraughts, amdwn(igout,k) ', (k,amdwn(igout,k), k=1,klev)
    PRINT *, 'wake-1, dowwdraughts, dtdwn(igout,k) ', (k,dtdwn(igout,k), k=1,klev)
    PRINT *, 'wake-1, dowwdraughts, dqdwn(igout,k) ', (k,dqdwn(igout,k), k=1,klev)
    PRINT *, 'wake-1, updraughts, amup(igout,k) ', (k,amup(igout,k), k=1,klev)
    PRINT *, 'wake-1, updraughts, dta(igout,k) ', (k,dta(igout,k), k=1,klev)
    PRINT *, 'wake-1, updraughts, dqa(igout,k) ', (k,dqa(igout,k), k=1,klev)
  ENDIF

  ! 2. - Prognostic part
  ! --------------------


  ! 2.1 - Undisturbed area and Wake integrals
  ! ---------------------------------------------------------

  DO i = KIDIA, KFDIA
    z(i) = 0.
    ktop(i) = 0
    kupper(i) = 0
    sum_thu(i) = 0.
    sum_tu(i) = 0.
    sum_qu(i) = 0.
    sum_thvu(i) = 0.
    sum_dth(i) = 0.
    sum_dq(i) = 0.
    sum_rho(i) = 0.
    sum_dtdwn(i) = 0.
    sum_dqdwn(i) = 0.

    av_thu(i) = 0.
    av_tu(i) = 0.
    av_qu(i) = 0.
    av_thvu(i) = 0.
    av_dth(i) = 0.
    av_dq(i) = 0.
    av_rho(i) = 0.
    av_dtdwn(i) = 0.
    av_dqdwn(i) = 0.
  END DO

  ! Distance between wakes
  DO i = KIDIA, KFDIA
    ll(i) = (1-sqrt(sigmaw(i)))/sqrt(wdens(i))
  END DO
  ! Potential temperatures and humidity
  ! ----------------------------------------------------------
  DO k = 1, klev
    DO i = KIDIA, KFDIA
      ! write(*,*)'wake 1',i,k,rd,te(i,k)
      rho(i, k) = p(i, k)/(rd*te(i,k))
      ! write(*,*)'wake 2',rho(i,k)
      IF (k==1) THEN
        ! write(*,*)'wake 3',i,k,rd,te(i,k)
        rhoh(i, k) = ph(i, k)/(rd*te(i,k))
        ! write(*,*)'wake 4',i,k,rd,te(i,k)
        zhh(i, k) = 0
      ELSE
        ! write(*,*)'wake 5',rd,(te(i,k)+te(i,k-1))
        rhoh(i, k) = ph(i, k)*2./(rd*(te(i,k)+te(i,k-1)))
        ! write(*,*)'wake 6',(-rhoh(i,k)*RG)+zhh(i,k-1)
        zhh(i, k) = (ph(i,k)-ph(i,k-1))/(-rhoh(i,k)*rg) + zhh(i, k-1)
      END IF
      ! write(*,*)'wake 7',ppi(i,k)
      the(i, k) = te(i, k)/ppi(i, k)
      thu(i, k) = (te(i,k)-deltatw(i,k)*sigmaw(i))/ppi(i, k)
      tu(i, k) = te(i, k) - deltatw(i, k)*sigmaw(i)
      qu(i, k) = qe(i, k) - deltaqw(i, k)*sigmaw(i)
      ! write(*,*)'wake 8',(rd*(te(i,k)+deltatw(i,k)))
      rhow(i, k) = p(i, k)/(rd*(te(i,k)+(1.-sigmaw(i))*deltatw(i,k)))  !RR bizarre, should be p(i, k)/(rd*(te(i,k)+(1-sigmaw(i)*deltatw(i,k)))
      dth(i, k) = deltatw(i, k)/ppi(i, k)
    END DO
  END DO

  DO k = 1, klev - 1
    DO i = KIDIA, KFDIA
      IF (k==1) THEN
        n2(i, k) = 0
      ELSE
        n2(i, k) = amax1(0., -rg**2/the(i,k)*rho(i,k)*(the(i,k+1)-the(i,k-1))/ &
                             (p(i,k+1)-p(i,k-1)))
      END IF
      zh(i, k) = (zhh(i,k)+zhh(i,k+1))/2

      cgw(i, k) = sqrt(n2(i,k))*zh(i, k)
      tgw(i, k) = coefgw*cgw(i, k)/ll(i)
    END DO
  END DO

  DO i = KIDIA, KFDIA
    n2(i, klev) = 0
    zh(i, klev) = 0
    cgw(i, klev) = 0
    tgw(i, klev) = 0
  END DO

  ! Calcul de la masse volumique moyenne de la colonne   (bdlmd)
  ! -----------------------------------------------------------------

  DO k = 1, klev
    DO i = KIDIA, KFDIA
      epaisseur1(i, k) = 0.
      epaisseur2(i, k) = 0.
    END DO
  END DO

  DO i = KIDIA, KFDIA
    epaisseur1(i, 1) = -(ph(i,2)-ph(i,1))/(rho(i,1)*rg) + 1.
    epaisseur2(i, 1) = -(ph(i,2)-ph(i,1))/(rho(i,1)*rg) + 1.
    rhow_moyen(i, 1) = rhow(i, 1)
  END DO

  DO k = 2, klev
    DO i = KIDIA, KFDIA
      epaisseur1(i, k) = -(ph(i,k+1)-ph(i,k))/(rho(i,k)*rg) + 1.
      epaisseur2(i, k) = epaisseur2(i, k-1) + epaisseur1(i, k)
      rhow_moyen(i, k) = (rhow_moyen(i,k-1)*epaisseur2(i,k-1)+rhow(i,k)* &
        epaisseur1(i,k))/epaisseur2(i, k)
    END DO
  END DO


  ! Choose an integration bound well above wake top
  ! -----------------------------------------------------------------

  ! Pupper = 50000.  ! melting level
  ! Pupper = 60000.
  ! Pupper = 80000.  ! essais pour case_e
  DO i = KIDIA, KFDIA
    pupper(i) = 0.6*ph(i, 1)
    pupper(i) = max(pupper(i), 45000.)
    ! cc        Pupper(i) = 60000.
  END DO

  IF (prt_level>=10) THEN
    PRINT *, 'wake-1.1, pupper(igout) ', pupper(igout)
  ENDIF


  ! Determine Wake top pressure (Ptop) from buoyancy integral
  ! --------------------------------------------------------

  ! -1/ Pressure of the level where dth becomes less than delta_t_min.

!<rr
  !DO i = KIDIA, KFDIA
  !  ptop_provis(i) = ph(i, 1)
  !END DO
  !DO k = 2, klev
  !  DO i = KIDIA, KFDIA
  !    IF (dth(i,k)>-delta_t_min .AND. dth(i,k-1)<-delta_t_min .AND. &
  !        ptop_provis(i)==ph(i,1)) THEN
  !      ptop_provis(i) = ((dth(i,k)+delta_t_min)*p(i,k-1)- &
  !                        (dth(i,k-1)+delta_t_min)*p(i,k))/(dth(i,k)-dth(i,k-1))
  !    END IF
  !  END DO
  !END DO

  DO i = KIDIA, KFDIA
    ptop_provis(i) = ph(i, klev)
  END DO
  DO k = klev-1,1,-1
    DO i = KIDIA, KFDIA
      IF (dth(i,k+1)>-delta_t_min .AND. dth(i,k)<-delta_t_min .AND. &
          ptop_provis(i)==ph(i,1)) THEN
        ptop_provis(i) = ((dth(i,k+1)+delta_t_min)*p(i,k)-(dth(i, &
          k)+delta_t_min)*p(i,k+1))/(dth(i,k+1)-dth(i,k))
      END IF
    END DO
  END DO
!>rr

  IF (prt_level>=10) THEN
    PRINT *, 'wake-1.2, ptop_provis(igout) ', ptop_provis(igout)
  ENDIF

  ! -2/ dth integral

  DO i = KIDIA, KFDIA
    sum_dth(i) = 0.
    dthmin(i) = -delta_t_min
    z(i) = 0.
  END DO

  DO k = 1, klev
    DO i = KIDIA, KFDIA
      dz(i) = -(amax1(ph(i,k+1),ptop_provis(i))-ph(i,k))/(rho(i,k)*rg)
      IF (dz(i)>0) THEN
        z(i) = z(i) + dz(i)
        sum_dth(i) = sum_dth(i) + dth(i, k)*dz(i)
        dthmin(i) = amin1(dthmin(i), dth(i,k))
      END IF
    END DO
  END DO

  ! -3/ height of triangle with area= sum_dth and base = dthmin

  DO i = KIDIA, KFDIA
    hw0(i) = 2.*sum_dth(i)/amin1(dthmin(i), -0.5)
    hw0(i) = amax1(hwmin, hw0(i))
  END DO

  ! -4/ now, get Ptop

  DO i = KIDIA, KFDIA
    z(i) = 0.
    ptop(i) = ph(i, 1)
  END DO

  DO k = 1, klev
    DO i = KIDIA, KFDIA
      dz(i) = amin1(-(ph(i,k+1)-ph(i,k))/(rho(i,k)*rg), hw0(i)-z(i))
      IF (dz(i)>0) THEN
        z(i) = z(i) + dz(i)
        ptop(i) = ph(i, k) - rho(i, k)*rg*dz(i)
      END IF
    END DO
  END DO

  IF (prt_level>=10) THEN
    PRINT *, 'wake-2, ptop_provis(igout), ptop(igout) ', ptop_provis(igout), ptop(igout)
  ENDIF


  ! -5/ Determination de ktop et kupper

  DO k = klev, 1, -1
    DO i = KIDIA, KFDIA
      IF (ph(i,k+1)<ptop(i)) ktop(i) = k
      IF (ph(i,k+1)<pupper(i)) kupper(i) = k
    END DO
  END DO

  ! On evite kupper = 1 et kupper = klev
  DO i = KIDIA, KFDIA
    kupper(i) = max(kupper(i), 2)
    kupper(i) = min(kupper(i), klev-1)
  END DO


  ! -6/ Correct ktop and ptop

  DO i = KIDIA, KFDIA
    ptop_new(i) = ptop(i)
  END DO
  DO k = klev, 2, -1
    DO i = KIDIA, KFDIA
      IF (k<=ktop(i) .AND. ptop_new(i)==ptop(i) .AND. &
          dth(i,k)>-delta_t_min .AND. dth(i,k-1)<-delta_t_min) THEN
        ptop_new(i) = ((dth(i,k)+delta_t_min)*p(i,k-1)-(dth(i, &
          k-1)+delta_t_min)*p(i,k))/(dth(i,k)-dth(i,k-1))
      END IF
    END DO
  END DO

  DO i = KIDIA, KFDIA
    ptop(i) = ptop_new(i)
  END DO

  DO k = klev, 1, -1
    DO i = KIDIA, KFDIA
      IF (ph(i,k+1)<ptop(i)) ktop(i) = k
    END DO
  END DO

  IF (prt_level>=10) THEN
    PRINT *, 'wake-3, ktop(igout), kupper(igout), ptop(igout) ', ktop(igout), kupper(igout), ptop(igout)
  ENDIF

  ! -5/ Set deltatw & deltaqw to 0 above kupper

  DO k = 1, klev
    DO i = KIDIA, KFDIA
      IF (k>=kupper(i)) THEN
        deltatw(i, k) = 0.
        deltaqw(i, k) = 0.
        d_deltatw2(i,k) = -deltatw0(i,k)
        d_deltaqw2(i,k) = -deltaqw0(i,k)        
      END IF
    END DO
  END DO
  IF (prt_level>=10) WRITE(*,*) 'WAKE deltatw', deltatw(1,:)


  ! Vertical gradient of LS omega

  DO k = 1, klev
    DO i = KIDIA, KFDIA
      IF (k<=kupper(i)) THEN
        dp_omgb(i, k) = (omgb(i,k+1)-omgb(i,k))/(ph(i,k+1)-ph(i,k))
      END IF
    END DO
  END DO

  ! Integrals (and wake top level number)
  ! --------------------------------------

  ! Initialize sum_thvu to 1st level virt. pot. temp.

  DO i = KIDIA, KFDIA
    z(i) = 1.
    dz(i) = 1.
    sum_thvu(i) = thu(i, 1)*(1.+epsim1*qu(i,1))*dz(i)
    sum_dth(i) = 0.
  END DO

  DO k = 1, klev
    DO i = KIDIA, KFDIA
      dz(i) = -(amax1(ph(i,k+1),ptop(i))-ph(i,k))/(rho(i,k)*rg)
      IF (dz(i)>0) THEN
        z(i) = z(i) + dz(i)
        sum_thu(i) = sum_thu(i) + thu(i, k)*dz(i)
        sum_tu(i) = sum_tu(i) + tu(i, k)*dz(i)
        sum_qu(i) = sum_qu(i) + qu(i, k)*dz(i)
        sum_thvu(i) = sum_thvu(i) + thu(i, k)*(1.+epsim1*qu(i,k))*dz(i)
        sum_dth(i) = sum_dth(i) + dth(i, k)*dz(i)
        sum_dq(i) = sum_dq(i) + deltaqw(i, k)*dz(i)
        sum_rho(i) = sum_rho(i) + rhow(i, k)*dz(i)
        sum_dtdwn(i) = sum_dtdwn(i) + dtdwn(i, k)*dz(i)
        sum_dqdwn(i) = sum_dqdwn(i) + dqdwn(i, k)*dz(i)
      END IF
    END DO
  END DO

  DO i = KIDIA, KFDIA
    hw0(i) = z(i)
  END DO
  IF (prt_level>=10) WRITE(*,*) 'WAKE hw0', hw0(1)


  ! 2.1 - WAPE and mean forcing computation
  ! ---------------------------------------

  ! ---------------------------------------

  ! Means

  DO i = KIDIA, KFDIA
    av_thu(i) = sum_thu(i)/hw0(i)
    av_tu(i) = sum_tu(i)/hw0(i)
    av_qu(i) = sum_qu(i)/hw0(i)
    av_thvu(i) = sum_thvu(i)/hw0(i)
    ! av_thve = sum_thve/hw0
    av_dth(i) = sum_dth(i)/hw0(i)
    av_dq(i) = sum_dq(i)/hw0(i)
    av_rho(i) = sum_rho(i)/hw0(i)
    av_dtdwn(i) = sum_dtdwn(i)/hw0(i)
    av_dqdwn(i) = sum_dqdwn(i)/hw0(i)

    wape(i) = -rg*hw0(i)*(av_dth(i)+ &
        epsim1*(av_thu(i)*av_dq(i)+av_dth(i)*av_qu(i)+av_dth(i)*av_dq(i)))/av_thvu(i)

  END DO
  IF (prt_level>=10) THEN
    WRITE(*,*) 'WAKE WAPE', wape(1)
    WRITE(*,*) 'WAKE av_dth', av_dth(1)
    WRITE(*,*) 'WAKE hw0', hw0(1)
  ENDIF

  ! 2.2 Prognostic variable update
  ! ------------------------------

  ! Filter out bad wakes

  DO k = 1, klev
    DO i = KIDIA, KFDIA
      IF (wape(i)<0.) THEN
        deltatw(i, k) = 0.
        deltaqw(i, k) = 0.
        dth(i, k) = 0.
        d_deltatw2(i,k) = -deltatw0(i,k)
        d_deltaqw2(i,k) = -deltaqw0(i,k)        
      END IF
    END DO
  END DO

  DO i = KIDIA, KFDIA
    IF (wape(i)<0.) THEN
      wape(i) = 0.
      cstar(i) = 0.
      hw(i) = hwmin
!jyg<
!!      sigmaw(i) = amax1(sigmad, sigd_con(i))
      sigmaw_targ = max(sigmad, sigd_con(i))
      d_sigmaw2(i) = d_sigmaw2(i) + sigmaw_targ - sigmaw(i)
      sigmaw(i) = sigmaw_targ
!>jyg
      fip(i) = 0.
      gwake(i) = .FALSE.
    ELSE
      hw(i) = hw0(i)
      cstar(i) = stark*sqrt(2.*wape(i))
      gwake(i) = .TRUE.
    END IF
  END DO


  ! Check qx and qw positivity
  ! --------------------------
  DO i = KIDIA, KFDIA
    q0_min(i) = min((qe(i,1)-sigmaw(i)*deltaqw(i,1)),  &
                    (qe(i,1)+(1.-sigmaw(i))*deltaqw(i,1)))
  END DO
  DO k = 2, klev
    DO i = KIDIA, KFDIA
      q1_min(i) = min((qe(i,k)-sigmaw(i)*deltaqw(i,k)), &
                      (qe(i,k)+(1.-sigmaw(i))*deltaqw(i,k)))
      IF (q1_min(i)<=q0_min(i)) THEN
        q0_min(i) = q1_min(i)
      END IF
    END DO
  END DO
  IF (prt_level>=10) WRITE(*,*) 'WAKE q0_min', q0_min

  DO i = KIDIA, KFDIA
    ok_qx_qw(i) = q0_min(i) >= 0.
    alpha(i) = 1.
  END DO

  IF (prt_level>=10) THEN
    PRINT *, 'wake-4, sigmaw(igout), cstar(igout), wape(igout), ktop(igout) ', &
                      sigmaw(igout), cstar(igout), wape(igout), ktop(igout)
  ENDIF


  ! C -----------------------------------------------------------------
  ! Sub-time-stepping
  ! -----------------

  nsub = 10
  dtimesub = dtime/nsub

  ! ------------------------------------------------------------
  DO isubstep = 1, nsub
    ! ------------------------------------------------------------

    ! wk_adv is the logical flag enabling wake evolution in the time advance
    ! loop
    DO i = KIDIA, KFDIA
      wk_adv(i) = ok_qx_qw(i) .AND. alpha(i) >= 1.
    END DO
    IF (prt_level>=10) THEN
      PRINT *, 'wake-4.1, isubstep,wk_adv(igout),cstar(igout),wape(igout), ptop(igout) ', &
                          isubstep,wk_adv(igout),cstar(igout),wape(igout), ptop(igout)
    ENDIF

    ! cc nrlmd   Ajout d'un recalcul de wdens dans le cas d'un entrainement
    ! négatif de ktop à kupper --------
    ! cc           On calcule pour cela une densité wdens0 pour laquelle on
    ! aurait un entrainement nul ---
    !jyg<
    ! Dans la configuration avec wdens prognostique, il s'agit d'un cas ou 
    ! les poches sont insuffisantes pour accueillir tout le flux de masse 
    ! des descentes unsaturees. Nous faisons alors l'hypothese que la 
    ! convection profonde cree directement de nouvelles poches, sans passer 
    ! par les thermiques. La nouvelle valeur de wdens est alors imposée.

    DO i = KIDIA, KFDIA
      ! c       print *,' isubstep,wk_adv(i),cstar(i),wape(i) ',
      ! c     $           isubstep,wk_adv(i),cstar(i),wape(i)
      IF (wk_adv(i) .AND. cstar(i)>0.01) THEN
        omg(i, kupper(i)+1) = -rg*amdwn(i, kupper(i)+1)/sigmaw(i) + &
                               rg*amup(i, kupper(i)+1)/(1.-sigmaw(i))
        wdens0 = (sigmaw(i)/(4.*3.14))* &
          ((1.-sigmaw(i))*omg(i,kupper(i)+1)/((ph(i,1)-pupper(i))*cstar(i)))**(2)
        IF (prt_level >= 10) THEN
             print*,'omg(i,kupper(i)+1),wdens0,wdens(i),cstar(i), ph(i,1)-pupper(i)', &
                     omg(i,kupper(i)+1),wdens0,wdens(i),cstar(i), ph(i,1)-pupper(i)
        ENDIF
        IF (wdens(i)<=wdens0*1.1) THEN
          IF (iflag_wk_pop_dyn >= 1) THEN
             d_wdens2(i) = d_wdens2(i) + wdens0 - wdens(i)
          ENDIF
          wdens(i) = wdens0
        END IF
      END IF
    END DO

    DO i = KIDIA, KFDIA
      IF (wk_adv(i)) THEN
        gfl(i) = 2.*sqrt(3.14*wdens(i)*sigmaw(i))
        rad_wk(i) = sqrt(sigmaw(i)/(3.14*wdens(i)))
!jyg<
!!        sigmaw(i) = amin1(sigmaw(i), sigmaw_max)
        sigmaw_targ = min(sigmaw(i), sigmaw_max)
        d_sigmaw2(i) = d_sigmaw2(i) + sigmaw_targ - sigmaw(i)
        sigmaw(i) = sigmaw_targ
!>jyg
      END IF
    END DO

    IF (iflag_wk_pop_dyn >= 1) THEN
!    The variable "death_rate" is significant only when iflag_wk_pop_dyn = 0.
!    Here, it has to be set to zero.
      death_rate(:) = 0.

      IF (iflag_wk_act ==2) THEN
      DO i = KIDIA, KFDIA
        IF (wk_adv(i)) THEN
          wape1_act(i) = abs(cin(i))
          wape2_act(i) = 2.*wape1_act(i) + 1.
          act(i) = min(1., max(0., (wape(i)-wape1_act(i)) / (wape2_act(i)-wape1_act(i)) ))
        ENDIF  ! (wk_adv(i))
      ENDDO
      ENDIF  ! (iflag_wk_act ==2)


      DO i = KIDIA, KFDIA
        IF (wk_adv(i)) THEN
!!          tau_wk(i) = max(rad_wk(i)/(3.*cstar(i))*((cstar(i)/cstart)**1.5 - 1), 100.)
          tau_wk_inv(i) = max( (3.*cstar(i))/(rad_wk(i)*((cstar(i)/cstart)**1.5 - 1)), 0.)
          tau_wk_inv_min = min(tau_wk_inv(i), 1./dtimesub)
          drdt(i) = (cstar(i) - wgen(i)*(sigmaw(i)/wdens(i)-aa0)/gfl(i)) / &
                    (1 + 2*f_shear(i)*(2.*sigmaw(i)-aa0*wdens(i)) - 2.*sigmaw(i))
!!                    (1 - 2*sigmaw(i)*(1.-f_shear(i)))
          drdt_pos=max(drdt(i),0.)

!!          d_wdens(i) = ( wgen(i)*(1.+2.*(sigmaw(i)-sigmad)) &
!!                     - wdens(i)*tau_wk_inv_min &
!!                     - 2.*gfl(i)*wdens(i)*Cstar(i) )*dtimesub
          d_awdens(i) = ( wgen(i) - (1./tau_cv)*(awdens(i) - act(i)*wdens(i)) )*dtimesub
          d_wdens(i) = ( wgen(i) - (wdens(i)-awdens(i))*tau_wk_inv_min -  &
                         2.*wdens(i)*gfl(i)*drdt_pos )*dtimesub
          d_wdens(i) = max(d_wdens(i), wdensmin-wdens(i))

!!          d_sigmaw(i) = ( (1.-2*f_shear(i)*sigmaw(i))*(gfl(i)*Cstar(i)+wgen(i)*sigmad/wdens(i)) &
!!                      + 2.*f_shear(i)*wgen(i)*sigmaw(i)**2/wdens(i) &
!!                      - sigmaw(i)*tau_wk_inv_min )*dtimesub
          d_sig_gen(i) = wgen(i)*aa0
          d_sig_death(i) = - sigmaw(i)*(1.-awdens(i)/wdens(i))*tau_wk_inv_min
!!          d_sig_col(i) = - 2*f_shear(i)*sigmaw(i)*gfl(i)*drdt_pos
          d_sig_col(i) = - 2*f_shear(i)*(2.*sigmaw(i)-wdens(i)*aa0)*gfl(i)*drdt_pos
          d_sigmaw(i) = ( d_sig_gen(i) + d_sig_death(i) + d_sig_col(i) + gfl(i)*cstar(i) )*dtimesub
          d_sigmaw(i) = max(d_sigmaw(i), sigmad-sigmaw(i))
        ENDIF
      ENDDO

      IF (prt_level >= 10) THEN
        print *,'wake, cstar(1), cstar(1)/cstart, rad_wk(1), tau_wk_inv(1), drdt(1) ', &
                       cstar(1), cstar(1)/cstart, rad_wk(1), tau_wk_inv(1), drdt(1)
        print *,'wake, wdens(1), awdens(1), act(1), d_awdens(1) ', &
                       wdens(1), awdens(1), act(1), d_awdens(1)
        print *,'wake, wgen, -(wdens-awdens)*tau_wk_inv, -2.*wdens*gfl*drdt_pos, d_wdens ', &
                       wgen(1), -(wdens(1)-awdens(1))*tau_wk_inv(1), -2.*wdens(1)*gfl(1)*drdt_pos, d_wdens(1)
        print *,'wake, d_sig_gen(1), d_sig_death(1), d_sig_col(1), d_sigmaw(1) ', &
                       d_sig_gen(1), d_sig_death(1), d_sig_col(1), d_sigmaw(1)
      ENDIF
    
    ELSE  !  (iflag_wk_pop_dyn >= 1)

    ! cc nrlmd

      DO i = KIDIA, KFDIA
        IF (wk_adv(i)) THEN
          ! cc nrlmd          Introduction du taux de mortalité des poches et
          ! test sur sigmaw_max=0.4
          ! cc         d_sigmaw(i) = gfl(i)*Cstar(i)*dtimesub
          IF (sigmaw(i)>=sigmaw_max) THEN
            death_rate(i) = gfl(i)*cstar(i)/sigmaw(i)
          ELSE
            death_rate(i) = 0.
          END IF
    
          d_sigmaw(i) = gfl(i)*cstar(i)*dtimesub - death_rate(i)*sigmaw(i)* &
            dtimesub
          ! $              - nat_rate(i)*sigmaw(i)*dtimesub
          ! c        print*, 'd_sigmaw(i),sigmaw(i),gfl(i),Cstar(i),wape(i),
          ! c     $  death_rate(i),ktop(i),kupper(i)',
          ! c     $	         d_sigmaw(i),sigmaw(i),gfl(i),Cstar(i),wape(i),
          ! c     $  death_rate(i),ktop(i),kupper(i)
    
          ! sigmaw(i) =sigmaw(i) + gfl(i)*Cstar(i)*dtimesub
          ! sigmaw(i) =min(sigmaw(i),0.99)     !!!!!!!!
          ! wdens = wdens0/(10.*sigmaw)
          ! sigmaw =max(sigmaw,sigd_con)
          ! sigmaw =max(sigmaw,sigmad)
        END IF
      END DO

    ENDIF   !  (iflag_wk_pop_dyn >= 1)
    IF (prt_level>=10) THEN
      WRITE(*,*) 'WAKE sigmaw', sigmaw(1)
      WRITE(*,*) 'WAKE d_sigmaw', d_sigmaw(1)
      WRITE(*,*) 'WAKE cstar', cstar(1)
      WRITE(*,*) 'WAKE death_rate', death_rate(1)
    ENDIF


    ! calcul de la difference de vitesse verticale poche - zone non perturbee
    ! IM 060208 differences par rapport au code initial; init. a 0 dp_deltomg
    ! IM 060208 et omg sur les niveaux de 1 a klev+1, alors que avant l'on definit
    ! IM 060208 au niveau k=1..?
    !JYG 161013 Correction : maintenant omg est dimensionne a klev.
    DO k = 1, klev
      DO i = KIDIA, KFDIA
        IF (wk_adv(i)) THEN !!! nrlmd
          dp_deltomg(i, k) = 0.
        END IF
      END DO
    END DO
    DO k = 1, klev
      DO i = KIDIA, KFDIA
        IF (wk_adv(i)) THEN !!! nrlmd
          omg(i, k) = 0.
        END IF
      END DO
    END DO

    DO i = KIDIA, KFDIA
      IF (wk_adv(i)) THEN
        z(i) = 0.
        omg(i, 1) = 0.
        dp_deltomg(i, 1) = -(gfl(i)*cstar(i))/(sigmaw(i)*(1-sigmaw(i)))
      END IF
    END DO

    DO k = 2, klev
      DO i = KIDIA, KFDIA
        IF (wk_adv(i) .AND. k<=ktop(i)) THEN
          dz(i) = -(ph(i,k)-ph(i,k-1))/(rho(i,k-1)*rg)
          z(i) = z(i) + dz(i)
          dp_deltomg(i, k) = dp_deltomg(i, 1)
          omg(i, k) = dp_deltomg(i, 1)*z(i)
        END IF
      END DO
    END DO
    IF (prt_level>=10) THEN
      WRITE(*,*) 'WAKE dp_deltomg', dp_deltomg(1,:)
      WRITE(*,*) 'WAKE omg', omg(1,:)
    ENDIF

    DO i = KIDIA, KFDIA
      IF (wk_adv(i)) THEN
        dztop(i) = -(ptop(i)-ph(i,ktop(i)))/(rho(i,ktop(i))*rg)
        ztop(i) = z(i) + dztop(i)
        omgtop(i) = dp_deltomg(i, 1)*ztop(i)
      END IF
    END DO

    IF (prt_level>=10) THEN
      PRINT *, 'wake-4.2, omg(igout,k) ', (k,omg(igout,k), k=1,klev)
      PRINT *, 'wake-4.2, omgtop(igout), ptop(igout), ktop(igout) ', &
                          omgtop(igout), ptop(igout), ktop(igout)
    ENDIF

    ! -----------------
    ! From m/s to Pa/s
    ! -----------------

    DO i = KIDIA, KFDIA
      IF (wk_adv(i)) THEN
        omgtop(i) = -rho(i, ktop(i))*rg*omgtop(i)
        dp_deltomg(i, 1) = omgtop(i)/(ptop(i)-ph(i,1))
      END IF
    END DO
    IF (prt_level>=10) WRITE(*,*) 'WAKE omgtop', omgtop(1)

    DO k = 1, klev
      DO i = KIDIA, KFDIA
        IF (wk_adv(i) .AND. k<=ktop(i)) THEN
          omg(i, k) = -rho(i, k)*rg*omg(i, k)
          dp_deltomg(i, k) = dp_deltomg(i, 1)
        END IF
      END DO
    END DO
    IF (prt_level>=10) THEN
      WRITE(*,*) 'WAKE omg', omg(1,:)
      WRITE(*,*) 'WAKE dp_deltomg', dp_deltomg(1,:)
    ENDIF

    ! raccordement lineaire de omg de ptop a pupper

    DO i = KIDIA, KFDIA
      IF (wk_adv(i) .AND. kupper(i)>ktop(i)) THEN
        omg(i, kupper(i)+1) = -rg*amdwn(i, kupper(i)+1)/sigmaw(i) + &
          rg*amup(i, kupper(i)+1)/(1.-sigmaw(i))
        dp_deltomg(i, kupper(i)) = (omgtop(i)-omg(i,kupper(i)+1))/ &
          (ptop(i)-pupper(i))
      END IF
    END DO
    IF (prt_level>=10) THEN
      WRITE(*,*) 'WAKE omg', omg(1,:)
      WRITE(*,*) 'WAKE dp_deltomg', dp_deltomg(1,:)
    ENDIF

    ! c      DO i=1,klon
    ! c        print*,'Pente entre 0 et kupper (référence)'
    ! c     $   	,omg(i,kupper(i)+1)/(pupper(i)-ph(i,1))
    ! c        print*,'Pente entre ktop et kupper'
    ! c     $  	,(omg(i,kupper(i)+1)-omgtop(i))/(pupper(i)-ptop(i))
    ! c      ENDDO
    ! c
    DO k = 1, klev
      DO i = KIDIA, KFDIA
        IF (wk_adv(i) .AND. k>ktop(i) .AND. k<=kupper(i)) THEN
          dp_deltomg(i, k) = dp_deltomg(i, kupper(i))
          omg(i, k) = omgtop(i) + (ph(i,k)-ptop(i))*dp_deltomg(i, kupper(i))
        END IF
      END DO
    END DO
!!    print *,'omg(igout,k) ', (k,omg(igout,k),k=1,klev)
    IF (prt_level>=10) THEN
      WRITE(*,*) 'WAKE omg', omg(1,:)
      WRITE(*,*) 'WAKE dp_deltomg', dp_deltomg(1,:)
    ENDIF
    ! cc nrlmd
    ! c      DO i=1,klon
    ! c      print*,'deltaw_ktop,deltaw_conv',omgtop(i),omg(i,kupper(i)+1)
    ! c      END DO
    ! cc


    ! --    Compute wake average vertical velocity omgbw


    DO k = 1, klev
      DO i = KIDIA, KFDIA
        IF (wk_adv(i)) THEN
          omgbw(i, k) = omgb(i, k) + (1.-sigmaw(i))*omg(i, k)
        END IF
      END DO
    END DO
    IF (prt_level>=10) WRITE(*,*) 'WAKE omgbw', omgbw(1,:)
    ! --    and its vertical gradient dp_omgbw

    DO k = 1, klev-1
      DO i = KIDIA, KFDIA
        IF (wk_adv(i)) THEN
          dp_omgbw(i, k) = (omgbw(i,k+1)-omgbw(i,k))/(ph(i,k+1)-ph(i,k))
        END IF
      END DO
    END DO
    DO i = 1, klon
      IF (wk_adv(i)) THEN
          dp_omgbw(i, klev) = 0.
      END IF
    END DO
    IF (prt_level>=10) WRITE(*,*) 'WAKE dp_omgbw', dp_omgbw(1,:)

    ! --    Upstream coefficients for omgb velocity
    ! --    (alpha_up(k) is the coefficient of the value at level k)
    ! --    (1-alpha_up(k) is the coefficient of the value at level k-1)
    DO k = 1, klev
      DO i = KIDIA, KFDIA
        IF (wk_adv(i)) THEN
          alpha_up(i, k) = 0.
          IF (omgb(i,k)>0.) alpha_up(i, k) = 1.
        END IF
      END DO
    END DO

    ! Matrix expressing [The,deltatw] from  [Th1,Th2]

    DO i = KIDIA, KFDIA
      IF (wk_adv(i)) THEN
        rre1(i) = 1. - sigmaw(i)
        rre2(i) = sigmaw(i)
      END IF
    END DO
    rrd1 = -1.
    rrd2 = 1.

    ! --    Get [Th1,Th2], dth and [q1,q2]

    DO k = 1, klev
      DO i = KIDIA, KFDIA
        IF (wk_adv(i) .AND. k<=kupper(i)+1) THEN
          dth(i, k) = deltatw(i, k)/ppi(i, k)
          th1(i, k) = the(i, k) - sigmaw(i)*dth(i, k) ! undisturbed area
          th2(i, k) = the(i, k) + (1.-sigmaw(i))*dth(i, k) ! wake
          q1(i, k) = qe(i, k) - sigmaw(i)*deltaqw(i, k) ! undisturbed area
          q2(i, k) = qe(i, k) + (1.-sigmaw(i))*deltaqw(i, k) ! wake
        END IF
      END DO
    END DO
    IF (prt_level>=10) WRITE(*,*) 'WAKE dth', dth(1,:)

    DO i = KIDIA, KFDIA
      IF (wk_adv(i)) THEN !!! nrlmd
        d_th1(i, 1) = 0.
        d_th2(i, 1) = 0.
        d_dth(i, 1) = 0.
        d_q1(i, 1) = 0.
        d_q2(i, 1) = 0.
        d_dq(i, 1) = 0.
      END IF
    END DO

    DO k = 2, klev
      DO i = KIDIA, KFDIA
        IF (wk_adv(i) .AND. k<=kupper(i)+1) THEN
          d_th1(i, k) = th1(i, k-1) - th1(i, k)
          d_th2(i, k) = th2(i, k-1) - th2(i, k)
          d_dth(i, k) = dth(i, k-1) - dth(i, k)
          d_q1(i, k) = q1(i, k-1) - q1(i, k)
          d_q2(i, k) = q2(i, k-1) - q2(i, k)
          d_dq(i, k) = deltaqw(i, k-1) - deltaqw(i, k)
        END IF
      END DO
    END DO

    DO i = KIDIA, KFDIA
      IF (wk_adv(i)) THEN
        omgbdth(i, 1) = 0.
        omgbdq(i, 1) = 0.
      END IF
    END DO

    DO k = 2, klev
      DO i = KIDIA, KFDIA
        IF (wk_adv(i) .AND. k<=kupper(i)+1) THEN !   loop on interfaces
          omgbdth(i, k) = omgb(i, k)*(dth(i,k-1)-dth(i,k))
          omgbdq(i, k) = omgb(i, k)*(deltaqw(i,k-1)-deltaqw(i,k))
        END IF
      END DO
    END DO

    IF (prt_level>=10) THEN
      PRINT *, 'wake-4.3, th1(igout,k) ', (k,th1(igout,k), k=1,klev)
      PRINT *, 'wake-4.3, th2(igout,k) ', (k,th2(igout,k), k=1,klev)
      PRINT *, 'wake-4.3, dth(igout,k) ', (k,dth(igout,k), k=1,klev)
      PRINT *, 'wake-4.3, omgbdth(igout,k) ', (k,omgbdth(igout,k), k=1,klev)
    ENDIF

    ! -----------------------------------------------------------------
    DO k = 1, klev-1
      DO i = KIDIA, KFDIA
        IF (wk_adv(i) .AND. k<=kupper(i)-1) THEN
          ! -----------------------------------------------------------------

          ! Compute redistribution (advective) term

          d_deltatw(i, k) = dtimesub/(ph(i,k)-ph(i,k+1))* &
            (rrd1*omg(i,k)*sigmaw(i)*d_th1(i,k) - &
             rrd2*omg(i,k+1)*(1.-sigmaw(i))*d_th2(i,k+1)- &
             (1.-alpha_up(i,k))*omgbdth(i,k)- &
             alpha_up(i,k+1)*omgbdth(i,k+1))*ppi(i, k)
!           print*,'d_deltatw=', k, d_deltatw(i,k)

          d_deltaqw(i, k) = dtimesub/(ph(i,k)-ph(i,k+1))* &
            (rrd1*omg(i,k)*sigmaw(i)*d_q1(i,k)- &
             rrd2*omg(i,k+1)*(1.-sigmaw(i))*d_q2(i,k+1)- &
             (1.-alpha_up(i,k))*omgbdq(i,k)- &
             alpha_up(i,k+1)*omgbdq(i,k+1))
!           print*,'d_deltaqw=', k, d_deltaqw(i,k)

          ! and increment large scale tendencies




          ! C
          ! -----------------------------------------------------------------
          d_te(i, k) = dtimesub*((rre1(i)*omg(i,k)*sigmaw(i)*d_th1(i,k)- &
                                  rre2(i)*omg(i,k+1)*(1.-sigmaw(i))*d_th2(i,k+1))/ &
                                 (ph(i,k)-ph(i,k+1)) &
                                 -sigmaw(i)*(1.-sigmaw(i))*dth(i,k)*(omg(i,k)-omg(i,k+1))/ &
                                 (ph(i,k)-ph(i,k+1)) )*ppi(i, k)

          d_qe(i, k) = dtimesub*((rre1(i)*omg(i,k)*sigmaw(i)*d_q1(i,k)- &
                                  rre2(i)*omg(i,k+1)*(1.-sigmaw(i))*d_q2(i,k+1))/ &
                                 (ph(i,k)-ph(i,k+1)) &
                                 -sigmaw(i)*(1.-sigmaw(i))*deltaqw(i,k)*(omg(i,k)-omg(i,k+1))/ &
                                 (ph(i,k)-ph(i,k+1)) )
        ELSE IF (wk_adv(i) .AND. k==kupper(i)) THEN
          d_te(i, k) = dtimesub*(rre1(i)*omg(i,k)*sigmaw(i)*d_th1(i,k)/(ph(i,k)-ph(i,k+1)))*ppi(i, k)

          d_qe(i, k) = dtimesub*(rre1(i)*omg(i,k)*sigmaw(i)*d_q1(i,k)/(ph(i,k)-ph(i,k+1)))

        END IF
        ! cc
      END DO
    END DO
    IF (prt_level>=10) THEN
      WRITE(*,*) 'WAKE d_deltatw', d_deltatw(1,:)
      WRITE(*,*) 'WAKE d_te', d_te(1,:)
    ENDIF
    ! ------------------------------------------------------------------

    IF (prt_level>=10) THEN
      PRINT *, 'wake-4.3, d_deltatw(igout,k) ', (k,d_deltatw(igout,k), k=1,klev)
      PRINT *, 'wake-4.3, d_deltaqw(igout,k) ', (k,d_deltaqw(igout,k), k=1,klev)
    ENDIF

    ! Increment state variables
!jyg<
    IF (iflag_wk_pop_dyn >= 1) THEN
      DO k = 1, klev
        DO i = KIDIA, KFDIA
          IF (wk_adv(i) .AND. k<=kupper(i)) THEN
            detr(i,k) = - d_sig_death(i) - d_sig_col(i)      
            entr(i,k) = d_sig_gen(i)
          ENDIF
        ENDDO
      ENDDO
      ELSE  ! (iflag_wk_pop_dyn >= 1)
      DO k = 1, klev
        DO i = KIDIA, KFDIA
          IF (wk_adv(i) .AND. k<=kupper(i)) THEN
            detr(i, k) = 0.
   
            entr(i, k) = 0.
          ENDIF
        ENDDO
      ENDDO
    ENDIF  ! (iflag_wk_pop_dyn >= 1)

    IF (prt_level>=10) PRINT *, 'wake-4.4'
    

    DO k = 1, klev
      DO i = KIDIA, KFDIA
        ! cc nrlmd       IF( wk_adv(i) .AND. k .LE. kupper(i)-1) THEN
        IF (wk_adv(i) .AND. k<=kupper(i)) THEN
          ! cc



          ! Coefficient de répartition

          crep(i, k) = crep_sol*(ph(i,kupper(i))-ph(i,k))/ &
            (ph(i,kupper(i))-ph(i,1))
          crep(i, k) = crep(i, k) + crep_upper*(ph(i,1)-ph(i,k))/ &
            (p(i,1)-ph(i,kupper(i)))


          ! Reintroduce compensating subsidence term.

          ! dtKE(k)=(dtdwn(k)*Crep(k))/sigmaw
          ! dtKE(k)=dtKE(k)-(dtdwn(k)*(1-Crep(k))+dta(k))
          ! .                   /(1-sigmaw)
          ! dqKE(k)=(dqdwn(k)*Crep(k))/sigmaw
          ! dqKE(k)=dqKE(k)-(dqdwn(k)*(1-Crep(k))+dqa(k))
          ! .                   /(1-sigmaw)

          ! dtKE(k)=(dtdwn(k)*Crep(k)+(1-Crep(k))*dta(k))/sigmaw
          ! dtKE(k)=dtKE(k)-(dtdwn(k)*(1-Crep(k))+dta(k)*Crep(k))
          ! .                   /(1-sigmaw)
          ! dqKE(k)=(dqdwn(k)*Crep(k)+(1-Crep(k))*dqa(k))/sigmaw
          ! dqKE(k)=dqKE(k)-(dqdwn(k)*(1-Crep(k))+dqa(k)*Crep(k))
          ! .                   /(1-sigmaw)

          dtke(i, k) = (dtdwn(i,k)/sigmaw(i)-dta(i,k)/(1.-sigmaw(i)))
          dqke(i, k) = (dqdwn(i,k)/sigmaw(i)-dqa(i,k)/(1.-sigmaw(i)))
          ! print*,'dtKE= ',dtKE(i,k),' dqKE= ',dqKE(i,k)

!

          ! cc nrlmd          Prise en compte du taux de mortalité
          ! cc               Définitions de entr, detr
!jyg<
!!            detr(i, k) = 0.
!!   
!!            entr(i, k) = detr(i, k) + gfl(i)*cstar(i) + &
!!              sigmaw(i)*(1.-sigmaw(i))*dp_deltomg(i, k)
!!
            entr(i, k) = entr(i,k) + gfl(i)*cstar(i) + &
                         sigmaw(i)*(1.-sigmaw(i))*dp_deltomg(i, k)   
!>jyg
            spread(i, k) = (entr(i,k)-detr(i,k))/sigmaw(i)

          ! cc        spread(i,k) =
          ! (1.-sigmaw(i))*dp_deltomg(i,k)+gfl(i)*Cstar(i)/
          ! cc     $  sigmaw(i)


          ! ajout d'un effet onde de gravité -Tgw(k)*deltatw(k) 03/02/06 YU
          ! Jingmei

          ! write(lunout,*)'wake.F ',i,k, dtimesub,d_deltat_gw(i,k),
          ! &  Tgw(i,k),deltatw(i,k)
          d_deltat_gw(i, k) = d_deltat_gw(i, k) - tgw(i, k)*deltatw(i, k)* &
            dtimesub
          ! write(lunout,*)'wake.F ',i,k, dtimesub,d_deltatw(i,k)
          ff(i) = d_deltatw(i, k)/dtimesub

          ! Sans GW

          ! deltatw(k)=deltatw(k)+dtimesub*(ff+dtKE(k)-spread(k)*deltatw(k))

          ! GW formule 1

          ! deltatw(k) = deltatw(k)+dtimesub*
          ! $         (ff+dtKE(k) - spread(k)*deltatw(k)-Tgw(k)*deltatw(k))

          ! GW formule 2

          IF (dtimesub*tgw(i,k)<1.E-10) THEN
            d_deltatw(i, k) = dtimesub*(ff(i)+dtke(i,k) - & 
               entr(i,k)*deltatw(i,k)/sigmaw(i) - &
               (death_rate(i)*sigmaw(i)+detr(i,k))*deltatw(i,k)/(1.-sigmaw(i)) - & ! cc
               tgw(i,k)*deltatw(i,k) )
          ELSE
            d_deltatw(i, k) = 1/tgw(i, k)*(1-exp(-dtimesub*tgw(i,k)))* &
               (ff(i)+dtke(i,k) - &
                entr(i,k)*deltatw(i,k)/sigmaw(i) - &
                (death_rate(i)*sigmaw(i)+detr(i,k))*deltatw(i,k)/(1.-sigmaw(i)) - &
                tgw(i,k)*deltatw(i,k) )
          END IF

          dth(i, k) = deltatw(i, k)/ppi(i, k)

          gg(i) = d_deltaqw(i, k)/dtimesub

          d_deltaqw(i, k) = dtimesub*(gg(i)+dqke(i,k) - & 
            entr(i,k)*deltaqw(i,k)/sigmaw(i) - &
            (death_rate(i)*sigmaw(i)+detr(i,k))*deltaqw(i,k)/(1.-sigmaw(i)))
          ! cc

          ! cc nrlmd
          ! cc       d_deltatw2(i,k)=d_deltatw2(i,k)+d_deltatw(i,k)
          ! cc       d_deltaqw2(i,k)=d_deltaqw2(i,k)+d_deltaqw(i,k)
          ! cc
        END IF
      END DO
    END DO
    IF (prt_level>=10) THEN
      WRITE(*,*) 'WAKE d_deltatw', d_deltatw(1,:)
      WRITE(*,*) 'WAKE dtke', dtke(1,:)
      WRITE(*,*) 'WAKE entr', entr(1,:)
      WRITE(*,*) 'WAKE detr', detr(1,:)
      WRITE(*,*) 'WAKE deltatw', deltatw(1,:)
      WRITE(*,*) 'WAKE dth', dth(1,:)
      WRITE(*,*) 'WAKE sigmaw', sigmaw(1)
      WRITE(*,*) 'WAKE tgw', tgw(1,:)
    ENDIF

    IF (prt_level>=10) PRINT *, 'wake-4.5'
    ! Scale tendencies so that water vapour remains positive in w and x.

    CALL wake_vec_modulation(klon, klev, wk_adv(KIDIA:KFDIA), epsilon, &
        qe(KIDIA:KFDIA,1:KLEV), d_qe(KIDIA:KFDIA,1:KLEV), &
        deltaqw(KIDIA:KFDIA,1:KLEV), d_deltaqw(KIDIA:KFDIA,1:KLEV), &
        sigmaw(KIDIA:KFDIA), d_sigmaw(KIDIA:KFDIA), alpha(KIDIA:KFDIA))
    
    IF (prt_level>=10) THEN
      PRINT *, 'wake-4.6'
      WRITE(*,*) 'WAKE alpha', alpha(1)
    ENDIF

    ! cc nrlmd
    ! c      print*,'alpha'
    ! c      do i=1,klon
    ! c         print*,alpha(i)
    ! c      end do
    ! cc
    DO k = 1, klev
      DO i = KIDIA, KFDIA
        IF (wk_adv(i) .AND. k<=kupper(i)) THEN
          d_te(i, k) = alpha(i)*d_te(i, k)
          d_qe(i, k) = alpha(i)*d_qe(i, k)
          d_deltatw(i, k) = alpha(i)*d_deltatw(i, k)
          d_deltaqw(i, k) = alpha(i)*d_deltaqw(i, k)
          d_deltat_gw(i, k) = alpha(i)*d_deltat_gw(i, k)
        END IF
      END DO
    END DO
    DO i = KIDIA, KFDIA
      IF (wk_adv(i)) THEN
        d_sigmaw(i) = alpha(i)*d_sigmaw(i)
      END IF
    END DO

    IF (prt_level>=10) PRINT *, 'wake-4.7'    
    ! Update large scale variables and wake variables
    ! IM 060208 manque DO i + remplace DO k=1,kupper(i)
    ! IM 060208     DO k = 1,kupper(i)
    DO k = 1, klev
      DO i = KIDIA, KFDIA
        IF (wk_adv(i) .AND. k<=kupper(i)) THEN
          dtls(i, k) = dtls(i, k) + d_te(i, k)
          dqls(i, k) = dqls(i, k) + d_qe(i, k)
          ! cc nrlmd
          d_deltatw2(i, k) = d_deltatw2(i, k) + d_deltatw(i, k)
          d_deltaqw2(i, k) = d_deltaqw2(i, k) + d_deltaqw(i, k)
          ! cc
        END IF
      END DO
    END DO
    IF (prt_level>=10) PRINT *, 'wake-4.8'    
    DO k = 1, klev
      DO i = KIDIA, KFDIA
        IF (wk_adv(i) .AND. k<=kupper(i)) THEN
          te(i, k) = te0(i, k) + dtls(i, k)
          qe(i, k) = qe0(i, k) + dqls(i, k)
          the(i, k) = te(i, k)/ppi(i, k)
          deltatw(i, k) = deltatw(i, k) + d_deltatw(i, k)
          deltaqw(i, k) = deltaqw(i, k) + d_deltaqw(i, k)
          dth(i, k) = deltatw(i, k)/ppi(i, k)
          ! c      print*,'k,qx,qw',k,qe(i,k)-sigmaw(i)*deltaqw(i,k)
          ! c     $        ,qe(i,k)+(1-sigmaw(i))*deltaqw(i,k)
        END IF
      END DO
    END DO
    IF (prt_level>=10) THEN
      PRINT *, 'wake-4.9'    
      WRITE(*,*) 'WAKE sigmaw', sigmaw(1)
      WRITE(*,*) 'WAKE d_sigmaw', d_sigmaw(1)
    ENDIF
    DO i = KIDIA, KFDIA
      IF (wk_adv(i)) THEN
        sigmaw(i) = sigmaw(i) + d_sigmaw(i)
        d_sigmaw2(i) = d_sigmaw2(i) + d_sigmaw(i)
      END IF
    END DO
!jyg<
    IF (prt_level>=10) PRINT *, 'wake-4.10'
    IF (iflag_wk_pop_dyn >= 1) THEN
      DO i = KIDIA, KFDIA
        IF (wk_adv(i)) THEN
          awdens(i) = awdens(i) + d_awdens(i)
          wdens(i) = wdens(i) + d_wdens(i)
          d_awdens2(i) = d_awdens2(i) + d_awdens(i)
          d_wdens2(i) = d_wdens2(i) + d_wdens(i)
        END IF
      END DO
      DO i = KIDIA, KFDIA
        IF (wk_adv(i)) THEN
          wdens_targ = max(wdens(i),wdensmin)
          d_wdens2(i) = d_wdens2(i) + wdens_targ - wdens(i)
          wdens(i) = wdens_targ
!
          wdens_targ = min( max(awdens(i),0.), wdens(i) )
          d_awdens2(i) = d_awdens2(i) + wdens_targ - awdens(i)
          awdens(i) = wdens_targ
        END IF
      END DO
      DO i = KIDIA, KFDIA
        IF (wk_adv(i)) THEN
          sigmaw_targ = max(sigmaw(i),sigmad)
          d_sigmaw2(i) = d_sigmaw2(i) + sigmaw_targ - sigmaw(i)
          sigmaw(i) = sigmaw_targ
        END IF
      END DO
    ENDIF  ! (iflag_wk_pop_dyn >= 1)
!>jyg

    IF (prt_level>=10) PRINT *, 'wake-4.11'
    ! Determine Ptop from buoyancy integral
    ! ---------------------------------------

    ! -     1/ Pressure of the level where dth changes sign.

!<rr
    !DO i = KIDIA, KFDIA
    !  IF (wk_adv(i)) THEN
    !    ptop_provis(i) = ph(i, 1)
    !  END IF
    !END DO

    !DO k = 2, klev
    !  DO i = KIDIA, KFDIA
    !    IF (wk_adv(i) .AND. ptop_provis(i)==ph(i,1) .AND. &
    !        dth(i,k)>-delta_t_min .AND. dth(i,k-1)<-delta_t_min) THEN
    !      ptop_provis(i) = ((dth(i,k)+delta_t_min)*p(i,k-1) - &
    !                        (dth(i,k-1)+delta_t_min)*p(i,k))/(dth(i,k)-dth(i,k-1))
    !    END IF
    !  END DO
    !END DO

    DO i = KIDIA, KFDIA
      IF (wk_adv(i)) THEN
        ptop_provis(i) = ph(i, klev)
      END IF
    END DO

    DO k = klev-1,1,-1
      DO i = KIDIA, KFDIA
        IF (wk_adv(i) .AND. ptop_provis(i)==ph(i,klev) .AND. &
            dth(i,k+1)>-delta_t_min .AND. dth(i,k)<-delta_t_min) THEN
          ptop_provis(i) = ((dth(i,k+1)+delta_t_min)*p(i,k)-(dth(i, &
            k)+delta_t_min)*p(i,k+1))/(dth(i,k+1)-dth(i,k))
        END IF
      END DO
    END DO
!>rr
    IF (prt_level>=10) THEN
      WRITE(*,*) 'WAKE ptop_provis', ptop_provis(1)
      PRINT *, 'wake-4.12'
    ENDIF
    ! -     2/ dth integral

    DO i = KIDIA, KFDIA
      IF (wk_adv(i)) THEN !!! nrlmd
        sum_dth(i) = 0.
        dthmin(i) = -delta_t_min
        z(i) = 0.
      END IF
    END DO

    DO k = 1, klev
      DO i = KIDIA, KFDIA
        IF (wk_adv(i)) THEN
          dz(i) = -(amax1(ph(i,k+1),ptop_provis(i))-ph(i,k))/(rho(i,k)*rg)
          IF (dz(i)>0) THEN
            z(i) = z(i) + dz(i)
            sum_dth(i) = sum_dth(i) + dth(i, k)*dz(i)
            dthmin(i) = amin1(dthmin(i), dth(i,k))
          END IF
        END IF
      END DO
    END DO
    IF (prt_level>=10) PRINT *, 'wake-4.13'
    ! -     3/ height of triangle with area= sum_dth and base = dthmin

    DO i = KIDIA, KFDIA
      IF (wk_adv(i)) THEN
        hw(i) = 2.*sum_dth(i)/amin1(dthmin(i), -0.5)
        hw(i) = amax1(hwmin, hw(i))
      END IF
    END DO
    IF (prt_level>=10) WRITE(*,*) 'WAKE hw', hw(1)

    ! -     4/ now, get Ptop
    IF (prt_level>=10) PRINT *, 'wake-4.14'
    DO i = KIDIA, KFDIA
      IF (wk_adv(i)) THEN !!! nrlmd
        ktop(i) = 0
        z(i) = 0.
      END IF
    END DO

    DO k = 1, klev
      DO i = KIDIA, KFDIA
        IF (wk_adv(i)) THEN
          dz(i) = amin1(-(ph(i,k+1)-ph(i,k))/(rho(i,k)*rg), hw(i)-z(i))
          IF (dz(i)>0) THEN
            z(i) = z(i) + dz(i)
            ptop(i) = ph(i, k) - rho(i, k)*rg*dz(i)
            ktop(i) = k
          END IF
        END IF
      END DO
    END DO
    IF (prt_level>=10) THEN
      PRINT *, 'wake-4.15'   
      WRITE(*,*) 'WAKE ptop ktop', ptop(1), ktop(1)
    ENDIF

    ! 4.5/Correct ktop and ptop

    DO i = KIDIA, KFDIA
      IF (wk_adv(i)) THEN
        ptop_new(i) = ptop(i)
      END IF
    END DO

!rr<
    !DO k = klev, 2, -1
    !  DO i = KIDIA, KFDIA
    !    ! IM v3JYG; IF (k .GE. ktop(i)
    !    IF (wk_adv(i) .AND. k<=ktop(i) .AND. ptop_new(i)==ptop(i) .AND. &
    !        dth(i,k)>-delta_t_min .AND. dth(i,k-1)<-delta_t_min) THEN
    !      ptop_new(i) = ((dth(i,k)+delta_t_min)*p(i,k-1) - &
    !                     (dth(i,k-1)+delta_t_min)*p(i,k))/(dth(i,k)-dth(i,k-1))
    !    END IF
    !  END DO
    !END DO

    DO k = klev-1, 1, -1
      DO i = KIDIA, KFDIA
        ! IM v3JYG; IF (k .GE. ktop(i)
        IF (wk_adv(i) .AND. k<=ktop(i)+1 .AND. ptop_new(i)==ptop(i) .AND. &
            dth(i,k+1)>-delta_t_min .AND. dth(i,k)<-delta_t_min) THEN
          ptop_new(i) = ((dth(i,k+1)+delta_t_min)*p(i,k)-(dth(i, &
            k)+delta_t_min)*p(i,k+1))/(dth(i,k+1)-dth(i,k))
        END IF
      END DO
    END DO
!>rr

    IF (prt_level>=10) PRINT *, 'wake-4.16'    

    DO i = KIDIA, KFDIA
      IF (wk_adv(i)) THEN
        ptop(i) = ptop_new(i)
      END IF
    END DO

    DO k = klev, 1, -1
      DO i = KIDIA, KFDIA
        IF (wk_adv(i)) THEN !!! nrlmd
          IF (ph(i,k+1)<ptop(i)) ktop(i) = k
        END IF
      END DO
    END DO

    IF (prt_level>=10) WRITE(*,*) 'WAKE ptop ktop', ptop(1), ktop(1)

    ! 5/ Set deltatw & deltaqw to 0 above kupper

    DO k = 1, klev
      DO i = KIDIA, KFDIA
        IF (wk_adv(i) .AND. k>=kupper(i)) THEN
          deltatw(i, k) = 0.
          deltaqw(i, k) = 0.
          d_deltatw2(i,k) = -deltatw0(i,k)
          d_deltaqw2(i,k) = -deltaqw0(i,k)
        END IF
      END DO
    END DO
    IF (prt_level>=10) PRINT *, 'wake-4.17'

    ! -------------Cstar computation---------------------------------
    DO i = KIDIA, KFDIA
      IF (wk_adv(i)) THEN !!! nrlmd
        sum_thu(i) = 0.
        sum_tu(i) = 0.
        sum_qu(i) = 0.
        sum_thvu(i) = 0.
        sum_dth(i) = 0.
        sum_dq(i) = 0.
        sum_rho(i) = 0.
        sum_dtdwn(i) = 0.
        sum_dqdwn(i) = 0.

        av_thu(i) = 0.
        av_tu(i) = 0.
        av_qu(i) = 0.
        av_thvu(i) = 0.
        av_dth(i) = 0.
        av_dq(i) = 0.
        av_rho(i) = 0.
        av_dtdwn(i) = 0.
        av_dqdwn(i) = 0.
      END IF
    END DO

    ! Integrals (and wake top level number)
    ! --------------------------------------

    ! Initialize sum_thvu to 1st level virt. pot. temp.

    DO i = KIDIA, KFDIA
      IF (wk_adv(i)) THEN !!! nrlmd
        z(i) = 1.
        dz(i) = 1.
        sum_thvu(i) = thu(i, 1)*(1.+epsim1*qu(i,1))*dz(i)
        sum_dth(i) = 0.
      END IF
    END DO
    IF (prt_level>=10) PRINT *, 'wake-4.18'
    DO k = 1, klev
      DO i = KIDIA, KFDIA
        IF (wk_adv(i)) THEN !!! nrlmd
          dz(i) = -(max(ph(i,k+1),ptop(i))-ph(i,k))/(rho(i,k)*rg)
          IF (dz(i)>0) THEN
            z(i) = z(i) + dz(i)
            sum_thu(i) = sum_thu(i) + thu(i, k)*dz(i)
            sum_tu(i) = sum_tu(i) + tu(i, k)*dz(i)
            sum_qu(i) = sum_qu(i) + qu(i, k)*dz(i)
            sum_thvu(i) = sum_thvu(i) + thu(i, k)*(1.+epsim1*qu(i,k))*dz(i)
            sum_dth(i) = sum_dth(i) + dth(i, k)*dz(i)
            sum_dq(i) = sum_dq(i) + deltaqw(i, k)*dz(i)
            sum_rho(i) = sum_rho(i) + rhow(i, k)*dz(i)
            sum_dtdwn(i) = sum_dtdwn(i) + dtdwn(i, k)*dz(i)
            sum_dqdwn(i) = sum_dqdwn(i) + dqdwn(i, k)*dz(i)
          END IF
        END IF
      END DO
    END DO

    DO i = KIDIA, KFDIA
      IF (wk_adv(i)) THEN !!! nrlmd
        hw0(i) = z(i)
      END IF
    END DO
    IF (prt_level>=10) PRINT *, 'wake-4.19'

    ! - WAPE and mean forcing computation
    ! ---------------------------------------

    ! ---------------------------------------

    ! Means

    DO i = KIDIA, KFDIA
      IF (wk_adv(i)) THEN !!! nrlmd
        av_thu(i) = sum_thu(i)/hw0(i)
        av_tu(i) = sum_tu(i)/hw0(i)
        av_qu(i) = sum_qu(i)/hw0(i)
        av_thvu(i) = sum_thvu(i)/hw0(i)
        av_dth(i) = sum_dth(i)/hw0(i)
        av_dq(i) = sum_dq(i)/hw0(i)
        av_rho(i) = sum_rho(i)/hw0(i)
        av_dtdwn(i) = sum_dtdwn(i)/hw0(i)
        av_dqdwn(i) = sum_dqdwn(i)/hw0(i)

        wape(i) = -rg*hw0(i)*(av_dth(i)+epsim1*(av_thu(i)*av_dq(i) + &
                              av_dth(i)*av_qu(i)+av_dth(i)*av_dq(i)))/av_thvu(i)
      END IF
    END DO
    IF (prt_level>=10) WRITE(*,*) 'WAKE wape', wape(1)

    ! Filter out bad wakes

    DO k = 1, klev
      DO i = KIDIA, KFDIA
        IF (wk_adv(i)) THEN !!! nrlmd
          IF (wape(i)<0.) THEN
            deltatw(i, k) = 0.
            deltaqw(i, k) = 0.
            dth(i, k) = 0.
            d_deltatw2(i,k) = -deltatw0(i,k)
            d_deltaqw2(i,k) = -deltaqw0(i,k)
          END IF
        END IF
      END DO
    END DO
    IF (prt_level>=10) PRINT *, 'wake-4.20'
    DO i = KIDIA, KFDIA
      IF (wk_adv(i)) THEN !!! nrlmd
        IF (wape(i)<0.) THEN
          wape(i) = 0.
          cstar(i) = 0.
          hw(i) = hwmin
!jyg<
!!          sigmaw(i) = max(sigmad, sigd_con(i))
          sigmaw_targ = max(sigmad, sigd_con(i))
          d_sigmaw2(i) = d_sigmaw2(i) + sigmaw_targ - sigmaw(i)
          sigmaw(i) = sigmaw_targ
!>jyg
          fip(i) = 0.
          gwake(i) = .FALSE.
        ELSE
          cstar(i) = stark*sqrt(2.*wape(i))
          gwake(i) = .TRUE.
        END IF
      END IF
    END DO
  IF (prt_level>=10) WRITE(*,*) 'WAKE deltatw', deltatw(1,:)
  END DO ! end sub-timestep loop

  IF (prt_level>=10) THEN
    PRINT *, 'wake-5, sigmaw(igout), cstar(igout), wape(igout), ptop(igout) ', &
                      sigmaw(igout), cstar(igout), wape(igout), ptop(igout)
  ENDIF


  ! ----------------------------------------------------------
  ! Determine wake final state; recompute wape, cstar, ktop;
  ! filter out bad wakes.
  ! ----------------------------------------------------------

  ! 2.1 - Undisturbed area and Wake integrals
  ! ---------------------------------------------------------

  DO i = KIDIA, KFDIA
    ! cc nrlmd       if (wk_adv(i)) then !!! nrlmd
    IF (ok_qx_qw(i)) THEN
      ! cc
      z(i) = 0.
      sum_thu(i) = 0.
      sum_tu(i) = 0.
      sum_qu(i) = 0.
      sum_thvu(i) = 0.
      sum_dth(i) = 0.
      sum_half_dth(i) = 0.
      sum_dq(i) = 0.
      sum_rho(i) = 0.
      sum_dtdwn(i) = 0.
      sum_dqdwn(i) = 0.

      av_thu(i) = 0.
      av_tu(i) = 0.
      av_qu(i) = 0.
      av_thvu(i) = 0.
      av_dth(i) = 0.
      av_dq(i) = 0.
      av_rho(i) = 0.
      av_dtdwn(i) = 0.
      av_dqdwn(i) = 0.

      dthmin(i) = -delta_t_min
    END IF
  END DO
  ! Potential temperatures and humidity
  ! ----------------------------------------------------------

  DO k = 1, klev
    DO i = KIDIA, KFDIA
      ! cc nrlmd       IF ( wk_adv(i)) THEN
      IF (ok_qx_qw(i)) THEN
        ! cc
        rho(i, k) = p(i, k)/(rd*te(i,k))
        IF (k==1) THEN
          rhoh(i, k) = ph(i, k)/(rd*te(i,k))
          zhh(i, k) = 0
        ELSE
          rhoh(i, k) = ph(i, k)*2./(rd*(te(i,k)+te(i,k-1)))
          zhh(i, k) = (ph(i,k)-ph(i,k-1))/(-rhoh(i,k)*rg) + zhh(i, k-1)
        END IF
        the(i, k) = te(i, k)/ppi(i, k)
        thu(i, k) = (te(i,k)-deltatw(i,k)*sigmaw(i))/ppi(i, k)
        tu(i, k) = te(i, k) - deltatw(i, k)*sigmaw(i)
        qu(i, k) = qe(i, k) - deltaqw(i, k)*sigmaw(i)
        rhow(i, k) = p(i, k)/(rd*(te(i,k)+(1.-sigmaw(i))*deltatw(i,k)))
        dth(i, k) = deltatw(i, k)/ppi(i, k)
      END IF
    END DO
  END DO

  ! Integrals (and wake top level number)
  ! -----------------------------------------------------------

  ! Initialize sum_thvu to 1st level virt. pot. temp.

  DO i = KIDIA, KFDIA
    ! cc nrlmd       IF ( wk_adv(i)) THEN
    IF (ok_qx_qw(i)) THEN
      ! cc
      z(i) = 1.
      dz(i) = 1.
      dz_half(i) = 1.
      sum_thvu(i) = thu(i, 1)*(1.+epsim1*qu(i,1))*dz(i)
      sum_dth(i) = 0.
    END IF
  END DO

  DO k = 1, klev
    DO i = KIDIA, KFDIA
      ! cc nrlmd       IF ( wk_adv(i)) THEN
      IF (ok_qx_qw(i)) THEN
        ! cc
        dz(i) = -(amax1(ph(i,k+1),ptop(i))-ph(i,k))/(rho(i,k)*rg)
        dz_half(i) = -(amax1(ph(i,k+1),0.5*(ptop(i)+ph(i,1)))-ph(i,k))/(rho(i,k)*rg)
        IF (dz(i)>0) THEN
          z(i) = z(i) + dz(i)
          sum_thu(i) = sum_thu(i) + thu(i, k)*dz(i)
          sum_tu(i) = sum_tu(i) + tu(i, k)*dz(i)
          sum_qu(i) = sum_qu(i) + qu(i, k)*dz(i)
          sum_thvu(i) = sum_thvu(i) + thu(i, k)*(1.+epsim1*qu(i,k))*dz(i)
          sum_dth(i) = sum_dth(i) + dth(i, k)*dz(i)
          sum_dq(i) = sum_dq(i) + deltaqw(i, k)*dz(i)
          sum_rho(i) = sum_rho(i) + rhow(i, k)*dz(i)
          sum_dtdwn(i) = sum_dtdwn(i) + dtdwn(i, k)*dz(i)
          sum_dqdwn(i) = sum_dqdwn(i) + dqdwn(i, k)*dz(i)
!
          dthmin(i) = min(dthmin(i), dth(i,k))
        END IF
        IF (dz_half(i)>0) THEN
          sum_half_dth(i) = sum_half_dth(i) + dth(i, k)*dz_half(i)
        END IF
      END IF
    END DO
  END DO

  DO i = KIDIA, KFDIA
    ! cc nrlmd       IF ( wk_adv(i)) THEN
    IF (ok_qx_qw(i)) THEN
      ! cc
      hw0(i) = z(i)
    END IF
  END DO

  ! - WAPE and mean forcing computation
  ! -------------------------------------------------------------

  ! Means

  DO i = KIDIA, KFDIA
    ! cc nrlmd       IF ( wk_adv(i)) THEN
    IF (ok_qx_qw(i)) THEN
      ! cc
      av_thu(i) = sum_thu(i)/hw0(i)
      av_tu(i) = sum_tu(i)/hw0(i)
      av_qu(i) = sum_qu(i)/hw0(i)
      av_thvu(i) = sum_thvu(i)/hw0(i)
      av_dth(i) = sum_dth(i)/hw0(i)
      av_dq(i) = sum_dq(i)/hw0(i)
      av_rho(i) = sum_rho(i)/hw0(i)
      av_dtdwn(i) = sum_dtdwn(i)/hw0(i)
      av_dqdwn(i) = sum_dqdwn(i)/hw0(i)

      wape2(i) = -rg*hw0(i)*(av_dth(i)+epsim1*(av_thu(i)*av_dq(i) + &
                             av_dth(i)*av_qu(i)+av_dth(i)*av_dq(i)))/av_thvu(i)
    END IF
  END DO
  IF (prt_level>=10) WRITE(*,*) 'WAKE wape2', wape2(1)


  ! Prognostic variable update
  ! ------------------------------------------------------------

  ! Filter out bad wakes

  IF (iflag_wk_check_trgl>=1) THEN
    ! Check triangular shape of dth profile
    DO i = 1, klon
      IF (ok_qx_qw(i)) THEN
        IF (prt_level>=10) THEN
          print *,'wake, hw0(i), dthmin(i) ', hw0(i), dthmin(i)
          print *,'wake, 2.*sum_dth(i)/(hw0(i)*dthmin(i)) ', &
                         2.*sum_dth(i)/(hw0(i)*dthmin(i))
          print *,'wake, sum_half_dth(i), sum_dth(i) ', &
                         sum_half_dth(i), sum_dth(i)
        ENDIF
        IF ((hw0(i) < 1.) .or. (dthmin(i) >= -delta_t_min) ) THEN
          wape2(i) = -1.
          IF (prt_level>=10) print *,'wake, rej 1'
        ELSE IF (iflag_wk_check_trgl==1.AND.abs(2.*sum_dth(i)/(hw0(i)*dthmin(i)) - 1.) > 0.5) THEN
          wape2(i) = -1.
          IF (prt_level>=10) print *,'wake, rej 2'
        ELSE IF (abs(sum_half_dth(i)) < 0.5*abs(sum_dth(i)) ) THEN
          wape2(i) = -1.
          IF (prt_level>=10) print *,'wake, rej 3'
        END IF
      END IF
    END DO
  END IF


  IF (prt_level>=10) WRITE(*,*) 'ok_qx_qw, wape2',ok_qx_qw(1), wape2(1)
  DO k = 1, klev
    DO i = KIDIA, KFDIA
      ! cc nrlmd        IF ( wk_adv(i) .AND. wape2(i) .LT. 0.) THEN
      IF (ok_qx_qw(i) .AND. wape2(i)<0.) THEN
        ! cc
        deltatw(i, k) = 0.
        deltaqw(i, k) = 0.
        dth(i, k) = 0.
        d_deltatw2(i,k) = -deltatw0(i,k)
        d_deltaqw2(i,k) = -deltaqw0(i,k)
      END IF
    END DO
  END DO


  DO i = KIDIA, KFDIA
    ! cc nrlmd       IF ( wk_adv(i)) THEN
    IF (ok_qx_qw(i)) THEN
      ! cc
      IF (wape2(i)<0.) THEN
        wape2(i) = 0.
        cstar2(i) = 0.
        hw(i) = hwmin
!jyg<
!!      sigmaw(i) = amax1(sigmad, sigd_con(i))
      sigmaw_targ = max(sigmad, sigd_con(i))
      d_sigmaw2(i) = d_sigmaw2(i) + sigmaw_targ - sigmaw(i)
      sigmaw(i) = sigmaw_targ
!>jyg
        fip(i) = 0.
        gwake(i) = .FALSE.
      ELSE
        IF (prt_level>=10) PRINT *, 'wape2>0'
        cstar2(i) = stark*sqrt(2.*wape2(i))
        gwake(i) = .TRUE.
      END IF
    END IF
  END DO

  DO i = KIDIA, KFDIA
    ! cc nrlmd       IF ( wk_adv(i)) THEN
    IF (ok_qx_qw(i)) THEN
      ! cc
      ktopw(i) = ktop(i)
    END IF
  END DO

  DO i = KIDIA, KFDIA
    ! cc nrlmd       IF ( wk_adv(i)) THEN
    IF (ok_qx_qw(i)) THEN
      ! cc
      IF (ktopw(i)>0 .AND. gwake(i)) THEN

        ! jyg1     Utilisation d'un h_efficace constant ( ~ feeding layer)
        ! cc       heff = 600.
        ! Utilisation de la hauteur hw
        ! c       heff = 0.7*hw
        heff(i) = hw(i)

        fip(i) = 0.5*rho(i, ktopw(i))*cstar2(i)**3*heff(i)*2* &
          sqrt(sigmaw(i)*wdens(i)*3.14)
        fip(i) = alpk*fip(i)
        ! jyg2
      ELSE
        fip(i) = 0.
      END IF
    END IF
  END DO

  ! Limitation de sigmaw

  ! cc nrlmd
  ! DO i=1,klon
  ! IF (OK_qx_qw(i)) THEN
  ! IF (sigmaw(i).GE.sigmaw_max) sigmaw(i)=sigmaw_max
  ! ENDIF
  ! ENDDO
  ! cc

  !jyg<
  IF (iflag_wk_pop_dyn >= 1) THEN
    DO i = KIDIA, KFDIA
      kill_wake(i) = ((wape(i)>=wape2(i)) .AND. (wape2(i)<=wapecut)) .OR. (ktopw(i)<=2) .OR. &
          .NOT. ok_qx_qw(i) .OR. (wdens(i) < 2.*wdensmin)
    ENDDO
  ELSE  ! (iflag_wk_pop_dyn >= 1)
    DO i = KIDIA, KFDIA
      kill_wake(i) = ((wape(i)>=wape2(i)) .AND. (wape2(i)<=wapecut)) .OR. (ktopw(i)<=2) .OR. &
          .NOT. ok_qx_qw(i)
    ENDDO
  ENDIF  ! (iflag_wk_pop_dyn >= 1)
  !>jyg

  IF (prt_level>=10) THEN
    WRITE(*,*) 'kill_wake:', kill_wake(1)
    WRITE(*,*) 'wape,wape2,wapecut, ktopw,ok_qx_qw:', wape(1),wape2(1),wapecut,ktopw(1),ok_qx_qw(1)
  ENDIF

  DO k = 1, klev
    DO i = KIDIA, KFDIA
!!jyg      IF (((wape(i)>=wape2(i)) .AND. (wape2(i)<=wapecut)) .OR. (ktopw(i)<=2) .OR. &
!!jyg          .NOT. ok_qx_qw(i)) THEN
      IF (kill_wake(i)) THEN
        ! cc
        dtls(i, k) = 0.
        dqls(i, k) = 0.
        deltatw(i, k) = 0.
        deltaqw(i, k) = 0.
        d_deltatw2(i,k) = -deltatw0(i,k)
        d_deltaqw2(i,k) = -deltaqw0(i,k)
      END IF  ! (kill_wake(i))
    END DO
  END DO

  DO i = KIDIA, KFDIA
!!jyg    IF (((wape(i)>=wape2(i)) .AND. (wape2(i)<=wapecut)) .OR. (ktopw(i)<=2) .OR. &
!!jyg        .NOT. ok_qx_qw(i)) THEN
      IF (kill_wake(i)) THEN
      ktopw(i) = 0
      wape(i) = 0.
      cstar(i) = 0.
!!jyg   Outside subroutine "Wake" hw, wdens and sigmaw are zero when there are no wakes
!!      hw(i) = hwmin                       !jyg
!!      sigmaw(i) = sigmad                  !jyg
      hw(i) = 0.                            !jyg
      fip(i) = 0.
!!      sigmaw(i) = 0.                        !jyg
      sigmaw_targ = 0.
      d_sigmaw2(i) = d_sigmaw2(i) + sigmaw_targ - sigmaw(i)
      sigmaw(i) = sigmaw_targ
      IF (iflag_wk_pop_dyn >= 1) THEN
!!        awdens(i) = 0.
!!        wdens(i) = 0.
        wdens_targ = 0.
        d_wdens2(i) = wdens_targ - wdens(i)
        wdens(i) = wdens_targ
        wdens_targ = 0.
        d_awdens2(i) = wdens_targ - awdens(i)
        awdens(i) = wdens_targ
      ENDIF  ! (iflag_wk_pop_dyn >= 1)
    ELSE  ! (kill_wake(i))
      wape(i) = wape2(i)
      cstar(i) = cstar2(i)
    END IF  ! (kill_wake(i))
    ! c        print*,'wape wape2 ktopw OK_qx_qw =',
    ! c     $          wape(i),wape2(i),ktopw(i),OK_qx_qw(i)
  END DO

  IF (prt_level>=10) THEN
    PRINT *, 'wake-6, wape wape2 ktopw OK_qx_qw =', &
                      wape(igout),wape2(igout),ktopw(igout),OK_qx_qw(igout)
  ENDIF


  ! -----------------------------------------------------------------
  ! Get back to tendencies per second

  DO k = 1, klev
    DO i = KIDIA, KFDIA

      ! cc nrlmd        IF ( wk_adv(i) .AND. k .LE. kupper(i)) THEN
!jyg<
!!      IF (ok_qx_qw(i) .AND. k<=kupper(i)) THEN
      IF (ok_qx_qw(i)) THEN
!>jyg
        ! cc
        dtls(i, k) = dtls(i, k)/dtime
        dqls(i, k) = dqls(i, k)/dtime
        d_deltatw2(i, k) = d_deltatw2(i, k)/dtime
        d_deltaqw2(i, k) = d_deltaqw2(i, k)/dtime
        d_deltat_gw(i, k) = d_deltat_gw(i, k)/dtime
        ! c      print*,'k,dqls,omg,entr,detr',k,dqls(i,k),omg(i,k),entr(i,k)
        ! c     $         ,death_rate(i)*sigmaw(i)
      END IF
    END DO
  END DO
!jyg<
  DO i = KIDIA, KFDIA
    d_sigmaw2(i) = d_sigmaw2(i)/dtime
    d_awdens2(i) = d_awdens2(i)/dtime
    d_wdens2(i) = d_wdens2(i)/dtime
  ENDDO
!>jyg

!  IF (LHOOK) CALL DR_HOOK('WAKE',1,ZHOOK_HANDLE)
END SUBROUTINE wake

SUBROUTINE wake_vec_modulation(nlon, nl, wk_adv, epsilon, qe, d_qe, deltaqw, &
    d_deltaqw, sigmaw, d_sigmaw, alpha)
  ! ------------------------------------------------------
  ! Dtermination du coefficient alpha tel que les tendances
  ! corriges alpha*d_G, pour toutes les grandeurs G, correspondent
  ! a une humidite positive dans la zone (x) et dans la zone (w).
  ! ------------------------------------------------------

  USE PARKIND1, ONLY : JPIM, JPRB
!  USE YOMHOOK,  ONLY : LHOOK,   DR_HOOK

  IMPLICIT NONE

  ! Input
  INTEGER(KIND=JPIM), INTENT(IN) :: nl
  INTEGER(KIND=JPIM), INTENT(IN) :: nlon

  REAL(KIND=JPRB), INTENT(IN) :: qe(nlon, nl)
  REAL(KIND=JPRB), INTENT(IN) :: d_qe(nlon, nl)
  REAL(KIND=JPRB), INTENT(IN) :: deltaqw(nlon, nl)
  REAL(KIND=JPRB), INTENT(IN) :: d_deltaqw(nlon, nl)
  REAL(KIND=JPRB), INTENT(IN) :: sigmaw(nlon)
  REAL(KIND=JPRB), INTENT(IN) :: d_sigmaw(nlon)
  REAL(KIND=JPRB), INTENT(IN) :: epsilon
  LOGICAL, INTENT(IN) :: wk_adv(nlon)
  ! Output
  REAL(KIND=JPRB), INTENT(OUT) :: alpha(nlon)
  ! Internal variables
  REAL(KIND=JPRB) zeta(nlon, nl)
  REAL(KIND=JPRB) alpha1(nlon)
  REAL(KIND=JPRB) x, a, b, c, discrim
  INTEGER(KIND=JPIM) i,k

  REAL(KIND=JPRB) :: ZHOOK_HANDLE

!  IF (LHOOK) CALL DR_HOOK('wake_vec_modulation',0,ZHOOK_HANDLE)

  DO k = 1, nl
    DO i = 1, nlon
      IF (wk_adv(i)) THEN
        IF ((deltaqw(i,k)+d_deltaqw(i,k))>=0.) THEN
          zeta(i, k) = 0.
        ELSE
          zeta(i, k) = 1.
        END IF
      END IF
    END DO
    DO i = 1, nlon
      IF (wk_adv(i)) THEN
        x = qe(i, k) + (zeta(i,k)-sigmaw(i))*deltaqw(i, k) + d_qe(i, k) + &
          (zeta(i,k)-sigmaw(i))*d_deltaqw(i, k) - d_sigmaw(i) * &
          (deltaqw(i,k)+d_deltaqw(i,k))
        a = -d_sigmaw(i)*d_deltaqw(i, k)
        b = d_qe(i, k) + (zeta(i,k)-sigmaw(i))*d_deltaqw(i, k) - &
          deltaqw(i, k)*d_sigmaw(i)
        c = qe(i, k) + (zeta(i,k)-sigmaw(i))*deltaqw(i, k) + epsilon
        discrim = b*b - 4.*a*c
        ! print*, 'x, a, b, c, discrim', x, a, b, c, discrim
        IF (a+b>=0.) THEN !! Condition suffisante pour la positivité de ovap
          alpha1(i) = 1.
        ELSE
          IF (x>=0.) THEN
            alpha1(i) = 1.
          ELSE
            IF (a>0.) THEN
              alpha1(i) = 0.9*min( (2.*c)/(-b+sqrt(discrim)),  &
                                   (-b+sqrt(discrim))/(2.*a) )
            ELSE IF (a==0.) THEN
              alpha1(i) = 0.9*(-c/b)
            ELSE
              ! print*,'a,b,c discrim',a,b,c discrim
              alpha1(i) = 0.9*max( (2.*c)/(-b+sqrt(discrim)),  &
                                   (-b+sqrt(discrim))/(2.*a))
            END IF
          END IF
        END IF
        alpha(i) = min(alpha(i), alpha1(i))
      END IF
    END DO
  END DO

  RETURN

!  IF (LHOOK) CALL DR_HOOK('wake_vec_modulation',1,ZHOOK_HANDLE)
END SUBROUTINE wake_vec_modulation
