import sys, os

import numpy as np

import ctypes as ct
_solib = ct.CDLL("wakelib/wake.so")

outfile=1

lverbose = True

def init_cst(out=outfile,printlev=1):
    """
    Initialize a few constants
    """

    _solib.sucst_.argtypes = [\
                              ct.POINTER(ct.c_int32),\
                              ct.POINTER(ct.c_int32),\
                              ]

    _solib.sucst_.restype  = ct.c_void_p    

    _solib.sucst_(\
                  ct.byref(ct.c_int32(out)),\
                  ct.byref(ct.c_int32(printlev)),\
                  )

def init_phy2(out=outfile,timestep=300.):
    """
    Initialize some environment parameter including timestep (default = 300 s)
    """

    _solib.suphy2_.argtypes = [\
                               ct.POINTER(ct.c_int32),\
                               ct.POINTER(ct.c_double),\
                               ]

    _solib.suphy2_.restype  = ct.c_void_p    

    _solib.suphy2_(\
                   ct.byref(ct.c_int32(out)),\
                   ct.byref(ct.c_double(timestep)),\
                  )


def init_wake(out=outfile,wapecut=1.,sigmad=0.02,hwmin=10.,sigmaw_max=0.4,\
              dens_rate=0.1,wdensmin=1.e-14,stark=0.33,alpk=0.25,\
              wdens_ref=(8.e-12,8.e-12),coefgw=4.,tau_cv=4000.,\
              crep_upper=0.9,crep_sol=1.0,delta_t_min=0.2,rzero=5000.):
    """
    Initialize wake parameters
    """

    wdens_ref_loc = np.copy(np.asfortranarray(wdens_ref,dtype=np.float64))

    _solib.suwake_.argtypes = [\
                               ct.POINTER(ct.c_int32),\
                               ct.POINTER(ct.c_double),\
                               ct.POINTER(ct.c_double),\
                               ct.POINTER(ct.c_double),\
                               ct.POINTER(ct.c_double),\
                               ct.POINTER(ct.c_double),\
                               ct.POINTER(ct.c_double),\
                               ct.POINTER(ct.c_double),\
                               ct.POINTER(ct.c_double),\
                               np.ctypeslib.ndpointer(dtype=np.float64,ndim=len(wdens_ref_loc.shape),flags='F_CONTIGUOUS'),\
                               ct.POINTER(ct.c_double),\
                               ct.POINTER(ct.c_double),\
                               ct.POINTER(ct.c_double),\
                               ct.POINTER(ct.c_double),\
                               ct.POINTER(ct.c_double),\
                               ct.POINTER(ct.c_double),\
                               ]
    _solib.suwake_.restype  = None    

    _solib.suwake_(\
                   ct.byref(ct.c_int32(out)),\
                   ct.byref(ct.c_double(wapecut)),\
                   ct.byref(ct.c_double(sigmad)),\
                   ct.byref(ct.c_double(hwmin)),\
                   ct.byref(ct.c_double(sigmaw_max)),\
                   ct.byref(ct.c_double(dens_rate)),\
                   ct.byref(ct.c_double(wdensmin)),\
                   ct.byref(ct.c_double(stark)),\
                   ct.byref(ct.c_double(alpk)),\
                   wdens_ref_loc,\
                   ct.byref(ct.c_double(coefgw)),\
                   ct.byref(ct.c_double(tau_cv)),\
                   ct.byref(ct.c_double(crep_upper)),\
                   ct.byref(ct.c_double(crep_sol)),\
                   ct.byref(ct.c_double(delta_t_min)),\
                   ct.byref(ct.c_double(rzero)),\
                   )

def init(kwargs):
    """
    Initialize what is needed for running the wake parameterization
    """

    init_cst()
    if kwargs.has_key('timestep'):
        init_phy2(timestep=kwargs['timestep'])
        kwargs_loc = kwargs.copy()
        del(kwargs_loc['timestep'])
    else:
        init_phy2()
        kwargs_loc = kwargs.copy()

    init_wake(**kwargs_loc)

def wake(znatsurf, p, ph, pi, dtime,te0, qe0, omgb,\
         dtdwn, dqdwn, amdwn, amup, dta, dqa, wgen,\
         sigd_con, cin,\
         deltatw, deltaqw, sigmaw, awdens, wdens):

    """
    Main call to wake Fortran routine
    """

    KIDIA = 1
    KLON  = p.shape[0]
    KFDIA = KIDIA + KLON - 1
    KLEV  = p.shape[1]

    znatsurf_loc = np.copy(np.asfortranarray(znatsurf,dtype=np.int32))
    p_loc        = np.copy(np.asfortranarray(p,       dtype=np.float64))
    ph_loc       = np.copy(np.asfortranarray(ph,      dtype=np.float64))
    pi_loc       = np.copy(np.asfortranarray(pi,      dtype=np.float64))
    te0_loc      = np.copy(np.asfortranarray(te0,     dtype=np.float64))
    qe0_loc      = np.copy(np.asfortranarray(qe0,     dtype=np.float64))
    omgb_loc     = np.copy(np.asfortranarray(omgb,    dtype=np.float64))
    dtdwn_loc    = np.copy(np.asfortranarray(dtdwn,   dtype=np.float64))
    dqdwn_loc    = np.copy(np.asfortranarray(dqdwn,   dtype=np.float64))
    amdwn_loc    = np.copy(np.asfortranarray(amdwn,   dtype=np.float64))
    amup_loc     = np.copy(np.asfortranarray(amup,    dtype=np.float64))
    dta_loc      = np.copy(np.asfortranarray(dta,     dtype=np.float64))
    dqa_loc      = np.copy(np.asfortranarray(dqa,     dtype=np.float64))
    wgen_loc     = np.copy(np.asfortranarray(wgen,    dtype=np.float64))
    sigd_con_loc = np.copy(np.asfortranarray(sigd_con,dtype=np.float64))
    cin_loc      = np.copy(np.asfortranarray(cin,     dtype=np.float64))
    deltatw_loc  = np.copy(np.asfortranarray(deltatw, dtype=np.float64))
    deltaqw_loc  = np.copy(np.asfortranarray(deltaqw, dtype=np.float64))
    sigmaw_loc   = np.copy(np.asfortranarray(sigmaw,  dtype=np.float64))
    awdens_loc   = np.copy(np.asfortranarray(awdens,  dtype=np.float64))
    wdens_loc    = np.copy(np.asfortranarray(wdens,   dtype=np.float64))

    # IN/INOUT
    _solib.wake_.argtypes = [\
                              ct.POINTER(ct.c_int32),# INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
                              ct.POINTER(ct.c_int32),# INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
                              ct.POINTER(ct.c_int32),# INTEGER(KIND=JPIM), INTENT(IN) :: KLON
                              ct.POINTER(ct.c_int32),# INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
                              np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1,flags='F_CONTIGUOUS'),# INTEGER(KIND=JPIM), INTENT(IN)  :: znatsurf(klon)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(IN)     :: p(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(IN)     :: ph(klon, klev+1)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(IN)     :: pi(klon, klev)
                              ct.POINTER(ct.c_double),# REAL(KIND=JPRB), INTENT(IN)    :: dtime
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(IN)     :: te0(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(IN)     :: qe0(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(IN)     :: omgb(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(IN)     :: dtdwn(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(IN)     :: dqdwn(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(IN)     :: amdwn(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(IN)     :: amup(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(IN)     :: dta(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(IN)     :: dqa(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(IN)     :: wgen(klon)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(IN)     :: sigd_con(klon)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(IN)     :: Cin(klon)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(INOUT)  :: deltatw(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(INOUT)  :: deltaqw(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(INOUT)  :: sigmaw(klon)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(INOUT)  :: awdens(klon)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(INOUT)  :: wdens(klon)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(OUT)    :: dth(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(OUT)    :: hw(klon)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(OUT)    :: wape(klon)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(OUT)    :: fip(klon)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(OUT)    :: gfl(klon)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(OUT)    :: dtls(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(OUT)    :: dqls(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1,flags='F_CONTIGUOUS'),# INTEGER(KIND=JPIM), INTENT(OUT) :: ktopw(klon)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(OUT)    :: omgbdth(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(OUT)    :: dp_omgb(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(OUT)    :: tu(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(OUT)    :: qu(klon, klev)          
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(OUT)    :: dtke(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(OUT)    :: dqke(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(OUT)    :: omg(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(OUT)    :: dp_deltomg(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(OUT)    :: spread(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(OUT)    :: cstar(klon)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(OUT)    :: d_deltat_gw(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(OUT)    :: d_deltatw2(klon, klev) 
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(OUT)    :: d_deltaqw2(klon, klev)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(OUT)    :: d_sigmaw2(klon)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(OUT)    :: d_awdens2(klon)
                              np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='F_CONTIGUOUS'),# REAL(KIND=JPRB), INTENT(OUT)    :: d_wdens2(klon)
                              ]

    _solib.wake_.restype  = None

    # OUT
    dth         = np.ndarray((KLON,KLEV),dtype=np.float64, order='F')
    tu          = np.ndarray((KLON,KLEV),dtype=np.float64, order='F')
    qu          = np.ndarray((KLON,KLEV),dtype=np.float64, order='F')
    dtls        = np.ndarray((KLON,KLEV),dtype=np.float64, order='F')
    dqls        = np.ndarray((KLON,KLEV),dtype=np.float64, order='F')
    dtke        = np.ndarray((KLON,KLEV),dtype=np.float64, order='F')
    dqke        = np.ndarray((KLON,KLEV),dtype=np.float64, order='F')
    omgbdth     = np.ndarray((KLON,KLEV),dtype=np.float64, order='F')
    omg         = np.ndarray((KLON,KLEV),dtype=np.float64, order='F')
    dp_omgb     = np.ndarray((KLON,KLEV),dtype=np.float64, order='F')
    dp_deltomg  = np.ndarray((KLON,KLEV),dtype=np.float64, order='F')
    d_deltat_gw = np.ndarray((KLON,KLEV),dtype=np.float64, order='F')
    hw          = np.ndarray((KLON,),    dtype=np.float64, order='F')
    wape        = np.ndarray((KLON,),    dtype=np.float64, order='F')
    fip         = np.ndarray((KLON,),    dtype=np.float64, order='F')
    gfl         = np.ndarray((KLON,),    dtype=np.float64, order='F')
    spread      = np.ndarray((KLON,KLEV),dtype=np.float64, order='F')
    cstar       = np.ndarray((KLON,),    dtype=np.float64, order='F')
    ktopw       = np.ndarray((KLON,),    dtype=np.int32,   order='F')
    d_deltatw2  = np.ndarray((KLON,KLEV),dtype=np.float64, order='F')
    d_deltaqw2  = np.ndarray((KLON,KLEV),dtype=np.float64, order='F')
    d_sigmaw2   = np.ndarray((KLON,),    dtype=np.float64, order='F')
    d_awdens2   = np.ndarray((KLON,),    dtype=np.float64, order='F')
    d_wdens2    = np.ndarray((KLON,),    dtype=np.float64, order='F')


    _solib.wake_(\
                  ct.byref(ct.c_int32(KIDIA)),\
                  ct.byref(ct.c_int32(KFDIA)),\
                  ct.byref(ct.c_int32(KLON)),\
                  ct.byref(ct.c_int32(KLEV)),\
                  znatsurf_loc,\
                  p_loc,\
                  ph_loc,\
                  pi_loc,\
                  ct.byref(ct.c_double(dtime)),\
                  te0_loc,\
                  qe0_loc,\
                  omgb_loc,\
                  dtdwn_loc,\
                  dqdwn_loc,\
                  amdwn_loc,\
                  amup_loc,\
                  dta_loc,\
                  dqa_loc,\
                  wgen_loc,\
                  sigd_con_loc,\
                  cin_loc,\
                  deltatw_loc,\
                  deltaqw_loc,\
                  sigmaw_loc,\
                  awdens_loc,\
                  wdens_loc,\
                  dth,\
                  hw,\
                  wape,\
                  fip,\
                  gfl,\
                  dtls,\
                  dqls,\
                  ktopw,\
                  omgbdth,\
                  dp_omgb,\
                  tu,\
                  qu,\
                  dtke,\
                  dqke,\
                  omg,\
                  dp_deltomg,\
                  spread,\
                  cstar,\
                  d_deltat_gw,\
                  d_deltatw2,\
                  d_deltaqw2,\
                  d_sigmaw2,\
                  d_awdens2,\
                  d_wdens2\
                  )

    return deltatw_loc, deltaqw_loc, sigmaw_loc,awdens_loc, wdens_loc,\
           dth, hw, wape, fip, gfl,\
           dtls, dqls, ktopw, omgbdth, dp_omgb, tu, qu,\
           dtke, dqke, omg, dp_deltomg, spread, cstar,\
           d_deltat_gw,\
           d_deltatw2, d_deltaqw2, d_sigmaw2, d_awdens2, d_wdens2

# Order is a priori important
wake_in = ['znatsurf','p','ph','pi','te','qe','omgbe','dtdwn','dqdwn','amdwn','amup','dta','dqa','wgen','sigd0','cin','dtime']
wake_inout = ['deltatw','deltaqw','sigmaw','awdens','wdens']
wake_out = ['dth','hw','wape','fip','gfl','dtls','dqls','ktopw','omgbdth','dp_omgb','tu','qu','dtke','dqke','omg','dp_deltomg','spread','cstar','d_deltat_gw','d_deltatw','d_deltaqw','d_sigmaw2','d_awdens2','d_wdens2']

def execute(data_in,data_inout):
    """
    Execute the wake main routine
        data_in: dictionnary containing input data
            znatsurf: type of surface (-)
            p:        full-level pressure (Pa)
            ph:       half-level pressure (Pa)
            pi:       Exner function (p/p0)**kappa (-)
            te:       environnement (grid-scale) temperature (K)
            qe:       environnement (grid-scale) specific humidity (kg kg-1)
            omgbe:    large-scale vertical velocity (Pa s-1)
            dtdwn:    temperature tendency due to downdraft (K s-1)
            dqdwn:    specific humidity tendency due to downdraft (kg kg-1 s-1)
            amdwn:    downdraft mass flux (per unit area of grid cell) (kg m-2 s-1) 
            amup:     updraft mass flux (per unit area of grid cell) (kg m-2 s-1)
            dta:      temperature tendency due to updraft (K s-1)
            dqa:      specific humidity tendency due to updraft (kg kg-1 s-1)
            wgen:     number of wakes generated per unit area and per second (m-2 s-1)
            sigd_con: fractional surface of downdraft (-)
            cin:      convective inhibition (J kg-1)
            dtime:    time step (s)

        data_inout: dictionnary containing in/out data
            deltatw: temperature difference between wake and off-wake regions (K)
            deltaqw: specific humidity difference between wake and off-wake regions (kg kg-1)
            sigmaw:  fractional area covered by wakes (-)
            awdens:  density of active wakes (m-2)
            wdens:   wake surface density (m-2)

    Return
        data_inout: dictionnary containing new in/out data
            deltatw: temperature difference between wake and off-wake regions (K)
            deltaqw: specific humidity difference between wake and off-wake regions (kg kg-1)
            sigmaw:  fractional area covered by wakes (-)
            awdens:  density of active wakes (m-2)
            wdens:   wake surface density (m-2)
        
        data_out : dictionnary containing output data
            dth:         potential temperature difference between wake and off-wake regions (K)
            hw:          wake depth (m)
            wape:        wake available potential energy (J kg-1)
            fip:         wake available lifting power (W m-2)
            gfl:         gust front length per unit area (m-1)
            dtls:        large-scale temperature tendency due to wakes (K s-1)
            dqls:        large-scale specific humidity tendency due to wakes (kg kg-1 s-1)
            ktopw:       index level of wake top (-)
            omgbdth:     flux of delta theta transported by large-scale omega (K Pa s-1)
            dp_omgb:     vertical gradient of large-scale omega (s-1)
            tu:          wake temperature (K)
            qu:          wake specific humidity (kg kg-1)
            dtke:        differential heating between wake and off-wake regions (K s-1)
            dqke:        differential moistening between wake and off-wake regions (kg kg s-1)
            omg:         vertical velocity difference between wake and off-wake regions (Pa s-1)
            dp_deltomg:  vertical gradient of omg (s-1)
            spread:      spreading term in wake tendencies
            cstar:       wake spreading velocity (m s-1)
            d_deltat_gw: delta T tendency due to gravity waves (K s-1)
            d_deltatw:   delta T tendency (K s-1)
            d_deltaqw:   delta qv tendency (kg kg-1 s-1)
            d_sigmaw2:   fractional wake area tendency (s-1)
            d_awdens2:   active wake density tendency (m-2 s-1)
            d_wdens2':   wake surface density tendency (m-2 s-1)
    """

    nlev, = data_in['p'].shape
    nlevp1 = nlev + 1

    znatsurf = np.zeros((1,),dtype=np.int32) + data_in['znatsurf']
    p        = np.reshape(data_in['p']    ,(1,nlev))
    ph       = np.reshape(data_in['ph']   ,(1,nlevp1))
    pi       = np.reshape(data_in['pi']   ,(1,nlev))
    te0      = np.reshape(data_in['te']   ,(1,nlev))
    qe0      = np.reshape(data_in['qe']   ,(1,nlev))
    omgb     = np.reshape(data_in['omgbe'],(1,nlev))
    dtdwn    = np.reshape(data_in['dtdwn'],(1,nlev))
    dqdwn    = np.reshape(data_in['dqdwn'],(1,nlev))
    amdwn    = np.reshape(data_in['amdwn'],(1,nlev))
    amup     = np.reshape(data_in['amup'] ,(1,nlev))
    dta      = np.reshape(data_in['dta']  ,(1,nlev))
    dqa      = np.reshape(data_in['dqa']  ,(1,nlev))
    wgen     = np.zeros((1,),dtype=np.float64) + data_in['wgen']
    sigd_con = np.zeros((1,),dtype=np.float64) + data_in['sigd0']
    cin      = np.zeros((1,),dtype=np.float64) + data_in['cin']

    deltatw  = np.reshape(data_inout['deltatw'],(1,nlev))
    deltaqw  = np.reshape(data_inout['deltaqw'],(1,nlev))
    sigmaw   = np.zeros((1,),dtype=np.float64) + data_inout['sigmaw']
    awdens   = np.zeros((1,),dtype=np.float64) + data_inout['awdens']
    wdens    = np.zeros((1,),dtype=np.float64) + data_inout['wdens']

    dtime = data_in['dtime']

    result = wake(znatsurf, p, ph, pi, dtime,te0, qe0, omgb, dtdwn, dqdwn, amdwn, amup, dta, dqa, wgen, sigd_con, cin, deltatw, deltaqw, sigmaw, awdens, wdens)

    i = 0
    data_inout_new = {}
    for v in wake_inout:
        data_inout_new[v] = result[i]
        i += 1

    data_out = {}
    for v in wake_out:
        data_out[v] = result[i]
        i += 1

    return data_inout_new, data_out

def itere(data_in,nstep=1,update=False):
    """
    Execute the wake main routine
        data_in: dictionnary containing input data
            znatsurf: type of surface (-)
            p:        full-level pressure (Pa)
            ph:       half-level pressure (Pa)
            pi:       Exner function (p/p0)**kappa (-)
            te:       environnement (grid-scale) temperature (K)
            qe:       environnement (grid-scale) specific humidity (kg kg-1)
            omgbe:    large-scale vertical velocity (Pa s-1)
            dtdwn:    temperature tendency due to downdraft (K s-1)
            dqdwn:    specific humidity tendency due to downdraft (kg kg-1 s-1)
            amdwn:    downdraft mass flux (per unit area of grid cell) (kg m-2 s-1) 
            amup:     updraft mass flux (per unit area of grid cell) (kg m-2 s-1)
            dta:      temperature tendency due to updraft (K s-1)
            dqa:      specific humidity tendency due to updraft (kg kg-1 s-1)
            wgen:     number of wakes generated per unit area and per second (m-2 s-1)
            sigd_con: fractional surface of downdraft (-)
            cin:      convective inhibition (J kg-1)
            timestep: time step (s)

        nstep: number of timestep to be iterated (default = 1)

        update: if True, prognostic variables of the wake parameterization are updated after each timestep
                if False, the parameterization is purely diagnostic at each timestep
                (default = False)

    Return
        data_out: dictionnary containing output data
            deltatw:     temperature difference between wake and off-wake regions (K)
            deltaqw:     specific humidity difference between wake and off-wake regions (kg kg-1)
            sigmaw:      fractional area covered by wakes (-)
            awdens:      density of active wakes (m-2)
            wdens:       wake surface density (m-2)
            dth:         potential temperature difference between wake and off-wake regions (K)
            hw:          wake depth (m)
            wape:        wake available potential energy (J kg-1)
            fip:         wake available lifting power (W m-2)
            gfl:         gust front length per unit area (m-1)
            dtls:        large-scale temperature tendency due to wakes (K s-1)
            dqls:        large-scale specific humidity tendency due to wakes (kg kg-1 s-1)
            ktopw:       index level of wake top (-)
            omgbdth:     flux of delta theta transported by large-scale omega (K Pa s-1)
            dp_omgb:     vertical gradient of large-scale omega (s-1)
            tu:          wake temperature (K)
            qu:          wake specific humidity (kg kg-1)
            dtke:        differential heating between wake and off-wake regions (K s-1)
            dqke:        differential moistening between wake and off-wake regions (kg kg s-1)
            omg:         vertical velocity difference between wake and off-wake regions (Pa s-1)
            dp_deltomg:  vertical gradient of omg (s-1)
            spread:      spreading term in wake tendencies
            cstar:       wake spreading velocity (m s-1)
            d_deltat_gw: delta T tendency due to gravity waves (K s-1)
            d_deltatw:   delta T tendency (K s-1)
            d_deltaqw:   delta qv tendency (kg kg-1 s-1)
            d_sigmaw2:   fractional wake area tendency (s-1)
            d_awdens2:   active wake density tendency (m-2 s-1)
            d_wdens2':   wake surface density tendency (m-2 s-1)
        
    """
    
    loc_in = {}
    loc_inout = {}

    res = {}
    for v in wake_inout + wake_out:
        res[v] = {}

    for it in range(0,nstep):
        if lverbose: print 'step {0}/{1}'.format(it+1,nstep)
        for v in wake_in:
            if v == 'dtime':
                loc_in[v] = data_in['timestep']
            else:
                loc_in[v] = np.copy(data_in[v][it])
        if update:
            if it == 0:
                for v in wake_inout:
                    loc_inout[v] = np.copy(data_in[v][it])
        else:
            for v in wake_inout:
                loc_inout[v] = np.copy(data_in[v][it])

        tmp_inout, tmp_out = execute(loc_in,loc_inout)

        if update:
            loc_inout = tmp_inout.copy()

        for v in wake_inout:
            res[v][it] = np.copy(tmp_inout[v])
        for v in wake_out:
            res[v][it] = np.copy(tmp_out[v])

    data_out = {}
    for v in res.keys():
        if len(res[v][0].shape) == 1:
            data_out[v] = np.zeros(nstep,dtype=np.float)
            l3d = False
        elif len(res[v][0].shape) == 2:
            nl,nlev = res[v][0].shape
            data_out[v] = np.zeros((nstep,nlev),dtype=np.float)
            l3d = True
        else:
            print 'shape length unexpected:',len(dataout[v][0].shape), 'for var:', var
            sys.exit()

        for it in range(0,nstep):
            data_out[v][it] = res[v][it]

    return data_out

def clean(outfile=outfile):
    """
    Some cleaning at the end of a run
    """
    
    os.remove('fort.{0}'.format(outfile))

