# -*- coding: utf-8 -*-
"""
ROMS default parameters

This module stores the default parameters for ROMS Bering Sea runs.

Created on Wed May  4 13:59:46 2016

@author: kakearney
"""

from collections import OrderedDict

# The default Bering 10K ocean.in (ROMS circa 2009) dictionary
# Most parameters are set to the values we use for Bering  10K runs.  However,
# all file names are set to defaults and need to be changed.
def ocean():
    """
    Populate dictionary with default ocean.in ROMS parameters

    Stores parameters related to the physical ROMS simulation that are written
    in the default standard input file (usually named ocean.in, or a variant
    thereof)

    Note that although the ROMS documentation states that parameters may be
    passed in any order, that's not entirely true.  Certain restrictions (such
    as NtileI coming before NtileJ) exist. So this function returns an ordered
    dictionary to preserve the order of keys for writing later.

    Returns:
        d: ROMS parameter dictionary.  Keys correspond to ROMS variables.
    """
    nat = 2
    npt = 167

    d = OrderedDict((
        # Application title
        ('TITLE', 'Bering Sea 10 km Grid'),
        # C-preprocessing Flag
        ('MyAppCPP', 'NEP5'),
        # Input variable information file name.
        ('VARNAME', 'varinfo.dat'),
        # Grid dimension parameters.
        ('Lm', 180),
        ('Mm', 256),
        ('N',   10),
        ('Nbed', 0),
        ('NAT', nat),
        ('NPT', npt),
        ('NCS', 0),
        ('NNS', 0),
        # Domain decomposition parameters
        ('NtileI', 7),
        ('NtileJ', 20),
        # Time-Stepping parameters
        ('NTIMES', 0),
        ('DT', 600),
        ('NDTFAST', 40),
        # Model iteration loops parameters.
        ('ERstr', 1),
        ('ERend', 1),
        ('Nouter', 1),
        ('Ninner', 1),
        ('Nintervals', 1),
        # Number of eigenvalues and eigenvectors for GST
        ('NEV', 2),
        ('NCV', 10),
        # Input/Output parameters.
        ('NRREC', 0),
        ('LcycleRST', True),
        ('NRST', 1008),
        ('NSTA', 6),
        ('NFLT', 144),
        ('NINFO', 1),
        # Output history, average, diagnostic files parameters.
        ('LDEFOUT', True),
        ('NHIS', 1008),
        ('NDEFHIS', 10080),
        ('NTSAVG', 1),
        ('NAVG', 1008),
        ('NDEFAVG', 10080),
        ('NTSDIA', 1),
        ('NDIA', 10),
        ('NDEFDIA', 0),
        # Output tangent linear and adjoint models parameters.
        ('LcycleTLM' , False),
             ('NTLM' , 72),
          ('NDEFTLM' , 0),
        ('LcycleADJ' , False),
             ('NADJ' , 72),
          ('NDEFADJ' , 0),
        # Output check pointing GST restart parameters.
           ('LrstGST' ,  False),
        ('MaxIterGST' ,  500),
              ('NGST' ,  10),
        # Relative accuracy of the Ritz values computed in the GST analysis.
        ('Ritz_tol' ,  1.0e-15),
        # Harmonic/biharmonic horizontal diffusion of tracers
        ('TNU2' , (nat+npt)*[25.0]),
        ('TNU4' , (nat+npt)*[0.0]),
        # Harmononic/biharmonic, horizontal viscosity coefficient
        ('VISC2' , 25.0),
        ('VISC4' , 0.0),
        # Vertical mixing coefficients for active tracers
        ('AKT_BAK' , (nat+npt)*[1.0e-6]),
        # Vertical mixing coefficient for momentum
        ('AKV_BAK' , 1.0e-5),
        # Turbulent closure parameters
        ('AKK_BAK' , 5.0e-6),
        ('AKP_BAK' , 5.0e-6),
         ('TKENU2' , 0.0),
         ('TKENU4' , 0.0),
        # Generic length-scale turbulence closure parameters
            ('GLS_P' , -1.0),
            ('GLS_M' , 0.5),
            ('GLS_N' , -1.0),
         ('GLS_Kmin' , 7.6e-6),
         ('GLS_Pmin' , 1.0e-12),
         ('GLS_CMU0' , 0.5477),
           ('GLS_C1' , 0.555),
           ('GLS_C2' , 0.833),
          ('GLS_C3M' , -0.6),
          ('GLS_C3P' , 1.0),
         ('GLS_SIGK' , 2.0),
         ('GLS_SIGP' , 2.0),
        # Constants used in momentum stress computation
          ('RDRG' , 3.0e-04),
         ('RDRG2' , 3.0e-03),
           ('Zob' , 0.02),
           ('Zos' , 0.02),
        # Height (m) of atmospheric measurements for Bulk fluxes parameterization
        ('BLK_ZQ' , 10.0),
        ('BLK_ZT' , 10.0),
        ('BLK_ZW' , 10.0),
        # Minimum depth for wetting and drying
        ('DCRIT' , 0.50),
        # Various parameters
          ('WTYPE' , 5),
        ('LEVSFRC' , 15),
        ('LEVBFRC' , 1),
        # Vertical S-coordinates parameters
        ('THETA_S' , 5.0),
        ('THETA_B' , 0.4),
         ('TCLINE' , 10.0),
        # Mean Density and Brunt-Vaisala frequency
            ('RHO0' ,  1025.0),
         ('BVF_BAK' ,  1.0e-4),
        # Time-stamps
        ('DSTART', ''),
        ('TIDE_START', -693962.0),
        ('TIME_REF', 19000101.0),
        # Nudging/relaxation time scales
         ('TNUDG' , 360.0),
         ('ZNUDG' , 360.0),
        ('M2NUDG' , 360.0),
        ('M3NUDG' , 360.0),
        # Factor between passive and active open boundary condition
        ('OBCFAC' , 120.0),
        # Linear equation of State parameters
             ('R0' , 1027.0),
             ('T0' , 10.0),
             ('S0' , 35.0),
          ('TCOEF' , 1.7e-4),
          ('SCOEF' , 7.6e-4),
        # Slipperiness parameter
        ('GAMMA2' , 1.0),
        # Starting  and ending day for adjoint sensitivity forcing
        ('DstrS' , 0.0),
        ('DendS' , 0.0),
        # Starting and ending vertical levels of the 3D adjoint state variables
        ('KstrS' , 1),
        ('KendS' , 1),
        # Specify the adjoint variables whose sensitivity is required
        ('Lstate' , {
            'isFsur' : False,
            'isUbar' : False,
            'isVbar' : False,
            'isUvel' : False,
            'isVvel' : False,
            'isTvar' : [False, False]
        }),
        # Stochastic optimals time decorrelation scale
        ('SO_decay' , 2.0),
        # Specify the surface forcing variables whose stochastic optimals are required.
        ('SOstate' , {
            'isUstr' : True,
            'isVstr' : True,
            'isTsur' : [False, False]
        }),
        # Stochastic optimals surface forcing standard deviation for dimensionalization.
        ('SO_sdev', {
            'isUstr' : 1.0,
            'isVstr' : 1.0,
            'isTsur' : [1.0, 1.0]
        }),
        # Activate writing of fields into HISTORY output file.
        ('Hout', {
            'idUvel'   : True,
            'idVvel'   : True,
            'idWvel'   : True,
            'idOvel'   : True,
            'idUbar'   : True,
            'idVbar'   : True,
            'idFsur'   : True,
            'idTvar'   : [True, True],
            'idUair'   : False,
            'idVair'   : False,
            'idUsms'   : True,
            'idVsms'   : True,
            'idUbms'   : True,
            'idVbms'   : True,
            'idUbrs'   : False,
            'idVbrs'   : False,
            'idUbws'   : False,
            'idVbws'   : False,
            'idUbcs'   : False,
            'idVbcs'   : False,
            'idUbot'   : False,
            'idVbot'   : False,
            'idUbur'   : False,
            'idVbvr'   : False,
            'idTsur'   : [True, True],
            'idLhea'   : True,
            'idShea'   : True,
            'idLrad'   : True,
            'idSrad'   : True,
            'idevap'   : False,
            'idrain'   : False,
            'idDano'   : False,
            'idVvis'   : True,
            'idTdif'   : True,
            'idSdif'   : False,
            'idHsbl'   : True,
            'idHbbl'   : True,
            'idMtke'   : False,
            'idMtls'   : False,
            'idUice'   : True,
            'idVice'   : True,
            'idAice'   : True,
            'idHice'   : True,
            'idTice'   : True,
            'idHsno'   : True,
            'idTimid'  : True,
            'idSfwat'  : True,
            'idTauiw'  : True,
            'idChuiw'  : True,
            'idAgeice' : True,
            'idSig11'  : True,
            'idSig12'  : True,
            'idSig22'  : True,
            'idS0mk'   : True,
            'idT0mk'   : True,
            'inert'    : npt*[True],
            'idBott'   : [True, True, True, True, True, True, True, True, True, False, False, False, False, False, False, False]
        }),
        # Generic User parameters
        ('NUSER' , 0),
         ('USER' , 0.0),
        # Input NetCDF file names
        ('GRDNAME', 'grd.nc'),
        ('ININAME', 'ini.nc'),
        ('IRPNAME', 'irp.nc'),
        ('IADNAME', 'iad.nc'),
        ('BRYNAME', 'bry.nc'),
        ('ADSNAME', 'ads.nc'),
        # Input forcing NetCDF file name(s).
        ('NFFILES', 1),
        ('FRCNAME', 'frc.nc'),
        # Output NetCDF file names
        ('RSTNAME', 'rst.nc'),
        ('HISNAME', 'his.nc'),
        ('AVGNAME', 'avg.nc'),
        ('STANAME', 'sta.nc'),
        ('FLTNAME', 'flt.nc'),
        ('GSTNAME', 'gst.nc'),
        ('TLMNAME', 'tlm.nc'),
        ('TLFNAME', 'tlf.nc'),
        ('ADJNAME', 'adj.nc'),
        # Input ASCII parameter filenames
        ('APARNAM', 'assimilation.in'),
        ('SPOSNAM', 'stations.in'),
        ('FPOSNAM', 'floats.in'),
        ('IPARNAM', 'ice.in'),
        ('BPARNAM', 'bio.in'),
        ('SPARNAM', 'sediment.in'),
        ('USRNAME', 'myfile.dat')
    ))
    return d

def bestnpz():
    """
    Populate dictionary with BESTNPZ biological parameters

    Stores parameters related to the BESTNPZ (Bering Sea NPZ) biological
    module.  These parameters are input through the file indicated by the
    BPARNAM variable.

    Returns:
        d: ROMS parameter dictionary.  Keys correspond to ROMS variables.
    """
    d = OrderedDict((
    # General
    ('Lbiology' , True),      # Logical flag for biology
    ('BioIter'  , 1),         # Number of iterations for nonlinear convergence.
    # Light
    ('PARfrac'  , 0.5),       # Fraction of irradiance that is photosynthetically available (PAR)
    ('k_ext'    , 0.046),     # Extinction coefficient due to seawater (1/m)
    ('k_chl'    , 0.121),     # Extinction coefficient due to Phy. (1/m)
    ('k_extZ'   , 0),         # Cokelet light parameter, TODO: need to get value from Georgina
    ('k_chlA'   , 0),         # Cokelet light parameter, TODO: need to get value from Georgina
    ('k_chlB'   , 0),         # Cokelet light parameter, TODO: need to get value from Georgina
    ('k_chlC'   , 0),         # Cokelet light parameter, TODO: need to get value from Georgina
    ('a_frac'   , 0.58),      # TODO: what is this
    ('a_mu1'    , 0.35),      # TODO: what is this
    # Biological conversions
    ('xi'       , 0.0126),    # Nitrogen,Carbon ratio (mmol N / mg C)
    ('ccr'      , 65.0),      # Carbon,Chlorophyll ratio (mg C / mg Chl-a), small phyto
    ('ccrPhL'   , 25.0),      # Carbon,Chlorophyll ratio (mg C / mg Chl-a), large phyto
    ('FeC'      , 1.667e-4),  # Fe(umol),Carbon(mg) is 2 umol Fe , mol C
    # Phytoplankton growth parameters
    ('DiS'      , 0.5),       # Doubling rate parameter
    ('DiL'      , 1.0),       # Doubling rate parameter
    ('DpS'      , 0.0275),    # Doubling rate exponent
    ('DpL'      , 0.0275),    # Doubling rate exponent
    ('k1PhS'    , 1.0),       # Half-saturation constant for NO3 limitation
    ('k1PhL'    , 2.0),       # Half-saturation constant for NO3 limitation
    ('k2PhS'    , 0.5),       # Half-saturation constant for NH4 limitation
    ('k2PhL'    , 2.0),       # Half-saturation constant for NH4 limitation
    ('FeCritPS' , 2.0),       # TODO: get description
    ('FeCritPL' , 2.0),       # TODO: get description
    ('kfePhS'   , 0.3),       # half-saturation const. PhS (umol per m-3)
    ('kfePhL'   , 1.0),       # half-saturation const. PhL (umol per m-3)
    # Feeding preference
    ('fpPhSMZL' , 1.0),       # PhS->MZL  Feeding preference
    ('fpPhLMZL' , 0.2),       # PhL->MZL  Feeding preference
    ('fpMZSMZL' , 0.0),       # PhS->Cop  Feeding preference
    ('fpPhSCop' , 0.8),       # PhL->Cop  Feeding preference
    ('fpPhLCop' , 0.7),       # MZL->Cop  Feeding preference
    ('fpMZLCop' , 0.5),       # PhS->NCa  Feeding preference
    ('fpPhSNCa' , 0.1),       # PhL->NCa  Feeding preference
    ('fpPhLNCa' , 1.0),       # MZL->NCa  Feeding preference
    ('fpMZLNCa' , 1.0),       # PhS->Eup  Feeding preference
    ('fpPhLEup' , 1.0),       # PhL->Eup  Feeding preference
    ('fpMZLEup' , 1.0),       # MZL->Eup  Feeding preference
    ('fpCopEup' , 0.2),       # Cop->Eup  Feeding preference
    ('fpDetEup' , 0.0),       # Det->Eup  Feeding preference TODO: get value
    ('fpDetEupO', 0.0),       # Det->EupO Feeding preference TODO: get value
    ('fpCopJel' , 1.0),       # Cop->Jel  Feeding preference
    ('fpNCaJel' , 1.0),       # NCa->Jel  Feeding preference
    ('fpEupJel' , 1.0),       # Eup->Jel  Feeding preference
    # Zooplankton growth and feeding
    ('eMZL'     ,  0.4),      # maximum specific ingestion rate (mg C/mg C/d)
    ('eCop'     ,  0.4),
    ('eNCa'     ,  0.3),
    ('eEup'     ,  0.3),
    ('eJel'     ,  0.069),
    ('Q10MZL'   ,  2.0),      # Q10 for growth rate
    ('Q10Cop'   ,  1.7),
    ('Q10NCa'   ,  1.6),
    ('Q10Eup'   ,  1.50),
    ('Q10Jele'  ,  2.4),
    ('Q10MZLT'  ,  5.0),      # Temperature coefficient for Q10 (deg. C)
    ('Q10CopT'  ,  5.0),
    ('Q10NCaT'  ,  5.0),
    ('Q10EupT'  ,  5.0),
    ('Q10JelTe' , 10.0),
    ('fMZL'     , 20.0),      # Half-saturation constant for grazing (mg C/m3)
    ('fCop'     , 30.0),
    ('fNCa'     , 30.0),
    ('fEup'     , 40.0),
    ('fJel'     ,  0.01),
    ('gammaMZL' ,  0.7),      # Growth efficiency
    ('gammaCop' ,  0.7),
    ('gammaNCa' ,  0.7),
    ('gammaEup' ,  0.7),
    ('gammaJel' ,  1.0),
    # Phytoplankton senescence
    ('mPhS'     , 0.01),      # daily linear mortality rate (1/d)
    ('mPhL'     , 0.01),      # daily linear mortality rate (1/d)
    # Predation closure
    ('mpredMZL' , 0.010),     # Daily mortality for Large Microzoo. (1/d)
    ('mpredCop' , 0.05),      # Daily mortality for Copepods (1/d)
    ('mpredNCa' , 0.05),      # Daily mortality for Neocalanus (1/d)
    ('mpredEup' , 0.05),      # Daily mortality for Euphausiids (1/d)
    ('mpredJel' , 0.006),
    # Sinking
    ('wPhS'     , 0.05),      # Sinking rate for Small Phytoplankton (m/d)
    ('wPhL'     , 1.0),       # Sinking rate for Large Phytoplankton (m/d)
    ('wDet'     , 1.0),       # Sinking rate for Detritus (m/d)
    ('wDetF'    , 10.0),      # Sinking rate for Detritus (m/d)
    # Respiration
    ('respPhS'  , 0.02),      # Specific respiration rate for PhS
    ('respPhL'  , 0.02),      # Specific respiration rate for PhL
    ('respMZL'  , 0.08),      # Specific respiration rate for MZL  0.1d0
    ('respCop'  , 0.04),      # Basal metabolic rate day**-1
    ('respNCa'  , 0.03),      # Basal metabolic rate day**-1
    ('respEup'  , 0.02),      # Basal metabolic rate day**-1 !0.044d0
    ('respJel'  , 0.02),      # Basal metabolic rate day**-1 0- try 2 - see above
    ('Q10Jelr'  , 2.8),       # Jel Q10 for respiration rate
    ('Q10JelTr' ,10.0),       # Jel Temperature coefficient for resp Q10 (deg. C)
    ('KtBm_PhS' , 0.03),      # temperature coefficient for respiration (1/deg C)
    ('KtBm_PhL' , 0.03),
    ('KtBm_MZL' , 0.069),
    ('ktbmC'    , 0.05),      # Temperature response degrees C**-1
    ('ktbmN'    , 0.05),      # Temperature response degrees C**-1
    ('ktbmE'    , 0.069),     # Temperature response degrees C**-1
    ('TmaxPhS'  ,10.0),       # reference temperature for respiration (deg C)
    ('TmaxPhL'  ,10.0),
    ('TmaxMZL'  , 8.0),
    ('TrefC'    ,15.0),       # Reference temperature degrees C
    ('TrefN'    , 5.0),       # Reference temperature degrees C
    ('TrefE'    , 5.0),       # Reference temperature degrees C
    ('TrefJ'    , 10.),       # Reference temperature degrees C
    # Iron climatology
    ('Feinlo'   ,   2.0),     # inshore/surface (micromol Fe m-3 or nM)
    ('Feinhi'   ,   4.0),     # inshore/deep    (micromol Fe m-3 or nM)
    ('Feinh'    ,  20.0),     # inshore isobath of transition (m)
    ('Feofflo'  ,   0.01),    # offshore/surface (micromol Fe m-3 or nM)
    ('Feoffhi'  ,   2.0),     # offshore/deep    (micromol Fe m-3 or nM)
    ('Feoffh'   , 100.0),     # offshore isobath of transition (m)
    # Diapause
    ('wNCrise'  ,  12.0),     # upward velocity (m/day), tuned not data
    ('wNCsink'  ,  11.0),     # downward velocity (m/day), tuned not data
    ('RiseStart',   0.0),     # Date NC begin to move upward (Day of Year)
    ('RiseEnd'  ,  60.0),     # Date NC stop moving upward (Day of Year)
    ('SinkStart', 155.0),     # Date NC begin to move downward (Day of Year)
    ('SinkEnd'  , 366.0),     # Date NC stop moving downward (Day of Year)
    # Reminaralization and nitrification
    ('Pv0'      , 0.1),       # PON dicompositon at 0 deg C (d-1)
    ('PvT'      , 0.069),     # Temperature coefficient (deg C-1)
    ('Nitr0'    , 0.0107),
    ('ktntr'    , 0.002),
    ('KNH4Nit'  , 0.057),     # Half Sat Con mg N/m3/day 0.08d0
    ('tI0'      , 0),         # Threshold for light limitation of nitrification (W m^-2) TODO: get value
    ('KI'       , 0),         # Half saturation light intensity for nitrification (w m^-2) TODO: get value
    # Benthos
    ('iremin'   , 0.800),       # related to nitrification bflx TODO: define
    ('q10r'     , 1.5),         # Ben Q10 for feeding
    ('Rup'      , 0.05),        # Ben Feeding rate (1/d)
    ('KupD'     , 2000.),       # Ben Half-saturation constant for feeding on pelagic food (mg C/m^2)
    ('KupP'     , 10.),         # Ben Half-saturation constant for feeding on benthic detritus (mg C/m^2)
    ('LupD'     , 292.),        # Ben Lower threshold for feeding on benthic detritus (mg C/m^2)
    ('LupP'     , 1.),          # Ben Lower threshold for feeding on pelagic food (mg C/m^2)
    ('Qres'     , 0.25),        # Ben Active metabolism reference? (1/d)
    ('Rres'     , 0.0027),      # Ben Basal metabolism reference? (1/d)
    ('rmort'    , 0.0021),      # Ben mortality rate (1/d)
    ('eex'      , [0.3, 0.3]),  # Ben fraction of pelagic food excreted
    ('eexD'     , [0.5, 0.7]),  # Ben fraction of detrital food excreted
    ('prefD'    , [1.0, 0.1]),  # BenDet->Ben feeding preference TODO: rename to match others?
    ('prefPL'   , [1.0, 0.1]),  # PhL->Ben feeding preference TODO: rename to match others?
    ('prefPS'   , [0.1, 0.1]),  # PhS->Ben feeding preference TODO: rename to match others?
    ('T0benr'   , 5.0),         # Ben Q10 Reference temperature for feeding
    ('BenPred'  , 0.000001),    # Ben mortality rate due to undefined predation (1/d)
    # Ice biology
    ('alphaIb'  , 0.80),        # IcePhL Chl-a specific attenuation coefficient (W^-1 m^-2)
    ('betaI'    , 0.018),       # IcePhL photosynthetic efficiency (W^-1 m^-2)
    ('inhib'    , 1.46),        # IcePhL NH4 inhibition on NO3 uptake (m^3/mmol N)
    ('ksnut1'   , 1.0),         # IcePhL half-saturation constant for NO3 (mmolN/m^3)
    ('ksnut2'   , 4.0),         # IcePhL half-saturation constant for NH4 (mmolN/m^3)
    ('mu0'      , 2.4),         # IcePhL maximum growth rate at 0 deg C (1/d)
    ('R0i'      , 0.05),        # IcePhL respiration rate (1/d)
    ('rg0'      , 0.01),        # IcePhL mortality rate at 0 deg C (1/d)
    ('rg'       , 0.03),        # IcePhL temperature coefficient for mortality (1/deg C)
    ('annit'    , 0.0149),      # IcePhL nitrification factor (1/d)
    ('aidz'     , 0.02),        # Ice thickness (m)
    # Horizontal diffusion of biological tracers
    ('TNU2'     , 15*[25.0] + [5.0]), # m2/s
    ('TNU4'     , 15*[0.0] + [2.0]),  # m2/s
    # Vertical mixing coefficients for biological tracers
    ('AKT_BAK'  , 15*[1.0e-6]),
    # Nudging/relaxation time scales
    ('TNUDG'    , 2*[360.0] + 12*[36000.0] + [360.0]),
    # Activate writing of biological tracers into history file
    ('Hout', OrderedDict((
        ('idTvar'    , 15*[True] + [True]), # State variables
        ('idTSvar'   , 16*[True]),          # STATIONARY
        ('idTS2var'  , 8*[True]),           # STATIONARY2
        ('idPT3var'  , 10*[True]),          # PROD3
        ('idPT2var'  , 3*[True]),           # PROD2
        ('idBvar'    , 3*[True]),           # BENTHIC
        ('idIceBvar' , 4*[True]),           # ICE_BIO
        ('idIcePhL'  , True),
        ('idIceNO3'  , True),
        ('idIceNH4'  , True),
        ('idIceLog'  , True),
        ('idTBFvar'  , True)                # BIOFLUX
        )))
    ))
    return d

def feast(flag=False):
    """
    Populate dictionary with FEAST biological parameters

    Stores parameters related to the FEAST (Bering Sea fish) biological
    module.  These parameters are input through the file indicated by the
    BPARNAM variable.  Includes BESTNPZ params for lower trophic levels.

    Returns:
        d: ROMS parameter dictionary.  Keys correspond to ROMS variables.
    """
    d = OrderedDict((
    ('feast_test'     , 1.0),    # dummy variable to test i/o
    ('force_ration'   , -0.7),   # set ration to <0 to do rations by functional re
    ('force_vonB'     , -1.0),   # >0 forces to vonB (pollock params used).  Overrides
    ('feast_mixed'    , 0),      # 0 : no mixed layer (all one water column) 1 : mixed
    ('feast_coupled'     , 1),   # 0 : no fish to zoop feedback 1 : feedback
    ('feast_mort'        , 1),   # 0 : no mzero or mtwo
    ('feast_fishing'     , 1),   # 0 : no fishing
    ('feast_growth'      , 1),   # 0 : no growth
    ('feast_recruitment' , 1),   # 0 : no recruitment
    ('feast_movement'    , 1),   # 0 : no movement
    ('feast_promotion'   , 0),   # 0 : no promotion
    ('feast_sp_view'     , 1),   # species to output detailed statistics
    ('feast_age_view'    , 2),   # age to output detailed statistics
    ('mpredCop'          , 0.035),  # Daily mortality for Copepods (1/d)
    ('mpredNCa'          , 0.025),  # Daily mortality for Neocalanus (1/d)
    ('mpredEup'          , 0.025),
    ('fpredCop'  ,   0.015),
    ('fpredNcaS'  ,  0.025),
    ('fpredNcaO'  ,  0.025),
    ('fpredEupS'  ,  0.025),
    ('fpredEupO'  ,  0.025),
    # Age offset (rows = species, cols = Species, Age 0, 1, 2, ... 10)
    ('fsh_age_offset' , [
        [1,  0,  0,  2,  3,  5,  6,  6,  7,  7,  8,  9],  # POL
        [2,  0,  0,  2,  4,  8,  9, 11, 12, 13, 14, 15],  # COD
        [3,  0,  0,  0,  3,  4,  5,  6,  7,  7,  8,  9],  # ATF
    ]),
    # For all below, cols = predators: POL, COD, ATF, HER, CAP, EUL, SAN, MYC, SAL1, SAL2
    # Size bin info
    ('fsh_Lsize'     ,  [4 ,  4,  4,  2,  2,  2,  2,  2,  2,  4]),
    ('fsh_base_prey' ,  [16, 16, 16, 16, 16, 16, 16, 16, 16, 16]),
    # Vertical water column preference
    ('fsh_a_T' , [0.45,  0.5, 0.45, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]),
    ('fsh_b_T' , [ 0.1, 0.15,  0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
    ('fsh_c_T' , [0.45,  0.8, 0.45, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
    # L-W relationship from VonB spreadsheets
    ('fsh_A_L' , [0.00553096, 0.00411781, 0.00443866, 0.008593946, 0.00033996, 0.00443866, 0.00581592, 0.00581592, 0.00581592, 0.00581592]),
    ('fsh_B_L' , [3.044172  , 3.25325765, 3.19894001, 3.107793351,     4.2304, 3.19894001,     3.0294,     3.0294,     3.0294,     3.0294]),
    # Functional response TOFIT
    ('fsh_Bv_min'  , [  0.01,   0.01,   0.01,   0.01,   0.01,   0.01,    0.01,   0.01,    0.01,    0.01]),
    ('fsh_B_Lzero' , [1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0,  1000.0, 1000.0,  1000.0,  1000.0]),
    ('fsh_B_Lone'  , [  30.0,   30.0,   30.0,   30.0,   30.0,   30.0,    30.0,   30.0,    30.0,    30.0]),
    ('fsh_B_Lpow'  , [   2.0,    2.0,    2.0,    2.0,    2.0,    2.0,     2.0,    2.0,     2.0,     2.0]),
    ('fsh_A_enc'   , [ 0.001,  0.002,  0.001,  0.001,  0.001,  0.001,  0.001 ,  0.001,   0.001,   0.001]),
    ('fsh_B_enc'   , [   2.0,    2.0,    2.0,    2.0,    2.0,   2.0 ,    2.0 ,   2.0 ,     2.0,     2.0]),
    # Using max stomach size in samples and digestion rates from REEM Samples
    ('fsh_A_S' , [0.072024, 0.01, 0.005,   0.00056460, 0.00056460, 0.00056460, 0.00056460, 0.00056460, 0.00056460, 0.00056460]),
    ('fsh_B_S' , [     1.2,  2.0,   1.7,   3.15843577, 3.15843577, 3.15843577, 3.15843577, 3.15843577, 3.15843577, 3.15843577]),
    # Standard bionenergetics consumption subtract 1 for standard g/g/day units on B_C
    ('fsh_C_TM' , [ 15,   15, 15 , 15 , 15 , 15 , 15 , 15 , 15 ,  15]),
    ('fsh_C_T0' , [ 10,   10, 10 , 10 , 10 , 10 , 10 , 10 , 10 ,  10]),
    ('fsh_C_Q'  , [2.6, 1.88, 5.5, 4.6, 4.6, 4.6, 4.6, 4.6, 4.6, 4.6]),
    # Caloric Density (in joules/g ww)
    ('fsh_ED_m' , [38.75, 38.75, 38.75,     0,     0,     0,     0,     0,     0,     0]),
    ('fsh_ED_b' , [ 2500,  2500, 2500 , 4499 , 4499 , 4499 , 4499 , 4499 , 4499 , 4499 ]),
    # Respiration subtract 1 for standard g/g/day units on B_R
    ('fsh_A_R'  , [0.0075, 0.008, 0.0057, 0.0195,  0.0195, 0.0195, 0.0195, 0.0195, 0.0195, 0.0195]),
    ('fsh_B_R'  , [ 0.749, 0.828,  0.656,   0.74,    0.74,   0.74,   0.74,   0.74,   0.74,   0.74]),
    ('fsh_F_A'  , [  0.15,  0.17,    0.2,    0.2,     0.2,    0.2,    0.2,    0.2,    0.2,    0.2]),
    ('fsh_U_A'  , [  0.11,  0.09,  0.111,  0.111,   0.111,  0.111,  0.111,  0.111,  0.111,  0.111]),
    ('fsh_SDA'  , [ 0.125,  0.17,  0.161,  0.125,   0.125,  0.125,  0.125,  0.125,  0.125,  0.125]),
    ('fsh_R_TM' , [    18,    24,   24.9,     18,      18,     18,     18,     18,     18,     18]),
    ('fsh_R_T0' , [    13,    21,   20.9,     15,      15,     15,     15,     15,     15,     15]),
    ('fsh_R_Q'  , [   2.6,  1.88,    5.5,    4.6,     4.6,    4.6,    4.6,    4.6,    4.6,    4.6]),
    # Mortality TOFIT
       ('fsh_omega' , [0.1, 0.2, 0.03, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]),
          ('fsh_mu' , [  1,   1,   1 ,  1 ,  1 ,  1 ,  1 ,  1 ,  1 ,  1 ]),
        ('fsh_zeta' , [  0,   0,   0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ]),
      ('fsh_s_mega' , [  0,   0,   0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ]),
    # Growth partitioning
    ('fsh_g_W' , 10*[1]),
    # Ivonne's recruitment functions
    ('fsh_mat_a'    , [    -7.8 ,    -12.4,       -9.4 ,     -12.5 ,        -6 ,       -10 ,      -7.5  ,      -10  ,    0.  ,    0.]),
    ('fsh_mat_b'    , [     0.2 ,      0.2,        0.2 ,       0.5 ,       0.5 ,       0.5 ,       0.5  ,      0.5  ,    0.  ,    0.]),
    ('fsh_fec_a'    , [     -15 ,    -13.4,     -19.65 ,     -12.5 ,        -6 ,       -10 ,      -7.5  ,      -10  ,    0.  ,    0.]),
    ('fsh_fec_b'    , [     0.3 ,      0.2,        0.3 ,       0.5 ,       0.5 ,       0.5 ,       0.5  ,      0.5  ,    0.  ,    0.]),
    ('fsh_fec_max'  , [    0.20 ,    0.045,  0.0730576 ,       0.1 ,       0.1 ,       0.1 ,       0.1  ,      0.1  ,  0.1  ,  0.1]),
    ('fsh_fem_prop' , [     0.5 ,      0.5,        0.5 ,       0.5 ,       0.5 ,       0.5 ,       0.5  ,      0.5  ,    0]), # short 1?
    ('fsh_rec_prop' , [  0.28384,    0.29 ,        1.36,     14.192,     14.192,     14.192,      14.192,     14.192,      0.,      0.]),
     ('fsh_z_muL'   , [     3.0 ,      2.4,        3.0 ,    1.3355 ,    1.3355 ,       3.0 ,    1.3355  ,   1.3355  ,     1. ,     1.]),
     ('fsh_z_sdL'   , [0.222717 , 0.455848,   0.222717 ,  0.222717 ,  0.222717 ,      0.25 ,  0.222717  , 0.222717  ,     1. ,     1.]),
    ('fsh_sp_sday'  , [     28  ,      30 ,        30  ,       930 ,       930 ,       930 ,       930  ,       930 ,   930 ,   930]),
    ('fsh_sp_eday'  , [     60  ,      60 ,        60  ,       960 ,       960 ,       960 ,       960  ,       960 ,   960 ,   960]),
    ('fsh_z_sday'   , [    105  ,      15 ,       120  ,      9100 ,      9150 ,      9120 ,      9210  ,      9150 ,    0  ,    0]),
    ('fsh_z_eday'   , [    273  ,     210 ,       240  ,      9180 ,      9240 ,      9240 ,      9300  ,      9210 ,    0  ,    0]),
    ('fsh_max_speed' , [ 1.0,  1.0,  1.0,  1.0,  1.0 , 1.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0]),
    ('fsh_happy_01'  , [-0.1, -0.1, -0.1, -0.1, -0.1 ,-0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1]),
    ('fsh_happy_99'  , [ 0.1,  0.1,  0.1,  0.1,  0.1 , 0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1]),
    # pred/prey matrices: SHR     SQU     EPI     CRA     OTH
    ('fsh_simple_len' , [   2,   10,    1,    2,    2]),  # mean length in cm
    ('fsh_simple_JperG', [3000, 3000, 3000, 3000, 3000]),  # cal per g wet weight
    ('fsh_simple_wt'   , [   5,   20,    5,   20,    5]),  # NOT USED g per N wet weight
    # j/gWW zoop dense             COP         NCAS       NCAO        EUPS        EUPO    BEN
    ('fsh_zoop_len'   , [       0.10,        0.35,        0.60,        1.70,        1.70,    1.0]),  # mean length in cm
    ('fsh_zoop_JperG' , [     2500.0,      3000.0,      3500.0,      4000.0,      4000.0, 2929.0]),  # cal per g wet weight
    ('fsh_zoop_wt'    , [0.000974293, 0.005333067, 0.014366590, 0.055869750, 0.055869750,    1.0]),  # g per N wet weight NOT USED
    # Predator/prey preferences
    ('fsh_q_G' , [
        [ 1,  1.0  ,  0.1  ,   0.01,  0.4  ,  0.4  ,  0.4  ,  0.1  ,  3.0  ,  2.0  ,  2.0  ,  1.9  , 4.0  ,   0.06 ,  0.001 ,  0.2],
        [ 2,  0.2  ,  0.3  ,   0.2 ,  1.0  ,  1.0  ,  1.0  ,  0.25 ,  1.0  ,  1.0  ,  1.0  ,  0.1  , 1.0  ,   0.5  ,  2.0   ,  0.5],
        [ 3,  1.2  ,  0.4  ,   0.2 ,  1.0  ,  1.0  ,  1.0  ,  0.25 ,  1.0  ,  1.0  ,  1.0  ,  0.1  , 1.0  ,   0.01 ,  0.0001,  0.3],
        [ 4,  0.1  ,  0.1  ,   0.1 ,  0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.05 , 2.0  ,   0.001,  0.0001,  0.2],
        [ 5,  0.1  ,  0.1  ,   0.1 ,  0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.05 , 0.5  ,   0.001,  0.0001,  0.2],
        [ 6,  0.1  ,  0.1  ,   0.1 ,  0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.05 , 0.5  ,   0.001,  0.0001,  0.2],
        [ 7,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001, 0.001,   1.0  ,  0.0001,  0.2],
        [ 8,  0.3  ,  0.3  ,   0.3 ,  2.0  ,  2.0  ,  2.0  ,  2.0  ,  2.0  ,  2.0  ,  2.0  ,  0.05 , 2.0  ,   0.001,  0.0001,  0.2],
        [ 9,  1.0  ,  1.0  ,   1.0 ,  2.0  ,  2.0  ,  2.0  ,  2.0  ,  2.0  ,  2.0  ,  2.0  ,  0.05 , 2.0  ,   0.001,  0.0001,  0.2],
        [10,  1.0  ,  1.0  ,   1.0 ,  2.0  ,  2.0  ,  2.0  ,  2.0  ,  2.0  ,  2.0  ,  2.0  ,  0.05 , 2.0  ,   0.001,  0.0001,  0.2],
        [11,  0.1  ,  0.1  ,   0.1 ,  1.0  , 1.0   , 1.0   ,  1.0  ,  1.0  ,  1.0  ,  1.0  ,  1.0  , 1.0  ,   1.0  ,  1.0   ,  1.0],
        [12,  0.4  ,  0.3  ,   0.1 ,  1.0  , 1.0   , 1.0   ,  1.0  ,  1.0  ,  1.0  ,  1.0  ,  1.0  , 1.0  ,   1.0  ,  1.0   ,  1.0],
        [13,  1.0  ,  1.0  ,   1.0 ,  1.0  , 1.0   , 1.0   ,  1.0  ,  1.0  ,  1.0  ,  1.0  ,  1.0  , 1.0  ,   1.0  ,  1.0   ,  1.0],
        [14,  1.0  ,  1.0  ,   1.0 ,  1.0  , 1.0   , 1.0   ,  1.0  ,  1.0  ,  1.0  ,  1.0  ,  1.0  , 1.0  ,   1.0  ,  1.0   ,  1.0],
        [15,  1.0  ,  1.0  ,   1.0 ,  1.0  , 1.0   , 1.0   ,  1.0  ,  1.0  ,  1.0  ,  1.0  ,  1.0  , 1.0  ,   1.0  ,  1.0   ,  1.0]
    ]),
    ('fsh_q_Gz' , [
        [ 1,  5.0,  3.0,  3.0,  6.0,  6.0,  0.0001],
        [ 2,  2.0,  2.0,  2.0,  0.5,  0.5,  0.001 ],
        [ 3,  2.0,  2.0,  2.0,  2.0,  2.0,  0.002 ],
        [ 4,  5.0,  2.0,  2.0,  4.0,  4.0,  0.0001],
        [ 5,  5.0,  2.0,  2.0,  4.0,  4.0,  0.0001],
        [ 6,  5.0,  2.0,  2.0,  4.0,  4.0,  0.0001],
        [ 7,  5.0,  2.0,  2.0,  4.0,  4.0,  0.0001],
        [ 8,  5.0,  2.0,  2.0,  4.0,  4.0,  0.0001],
        [ 9,  5.0,  2.0,  2.0,  4.0,  4.0,  0.0001],
        [10,  5.0,  2.0,  2.0,  4.0,  4.0,  0.0001],
        [11,  1.0,  1.0,  1.0,  1.0,  1.0,  0.0],
        [12,  1.0,  1.0,  1.0,  1.0,  1.0,  0.0],
        [13,  1.0,  1.0,  1.0,  1.0,  1.0,  0.0],
        [14,  1.0,  1.0,  1.0,  1.0,  1.0,  0.0],
        [15,  1.0,  1.0,  1.0,  1.0,  1.0,  0.0]
    ]),
    ('fsh_alpha_G', [
        [ 1,  16.975, 12.975, 12.975, 16.975, 16.975, 12.975, 16.975, 160.73, 12.975, 12.975, 234.73, 140.73, 50.975, 12.975, 228.62],
        [ 2,  5.4260, 5.4260, 5.4260, 5.4260, 5.4260, 5.4260, 5.4260, 5.4260, 5.4260, 5.4260, 41.372, 23.502, 80.135, 133.02, 12.505],
        [ 3,  3.0   , 6.6279, 6.6279, 6.6279, 6.6279, 6.6279, 6.6279, 6.6279, 6.6279, 6.6279, 27.298, 3.4421, 427.95, 21.797, 5.6126],
        [ 4,  18.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 194.73, 3.3176, 12.975, 12.975, 228.62],
        [ 5,  18.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 194.73, 3.3176, 12.975, 12.975, 228.62],
        [ 6,  18.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 194.73, 3.3176, 12.975, 12.975, 228.62],
        [ 7,  18.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 194.73, 3.3176, 12.975, 12.975, 228.62],
        [ 8,  18.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 194.73, 3.3176, 12.975, 12.975, 228.62],
        [ 9,  18.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 194.73, 3.3176, 12.975, 12.975, 228.62],
        [10,  18.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 12.975, 194.73, 3.3176, 12.975, 12.975, 228.62]
    ]),
    ('fsh_alpha_Gz', [
         [ 1,  30.805, 41.922, 35.662, 22.503, 22.503, 12.975],
         [ 2,  30.805, 41.922, 35.662, 19.322, 19.322, 98.658],
         [ 3,  289.89, 144.95, 97.844, 6.3573, 6.3573, 61.732],
         [ 4,  27    , 35    , 34.662, 26    , 26    , 12.975],
         [ 5,  27    , 35    , 34.662, 26    , 26    , 12.975],
         [ 6,  27    , 35    , 34.662, 26    , 26    , 12.975],
         [ 7,  27    , 35    , 34.662, 26    , 26    , 12.975],
         [ 8,  27    , 35    , 34.662, 26    , 26    , 12.975],
         [ 9,  27    , 35    , 34.662, 26    , 26    , 12.975],
         [10,  27    , 35    , 34.662, 26    , 26    , 12.975]
    ]),
    ('fsh_beta_G', [
        [ 1,  0.1647, 0.1647, 0.1647,  0.1647,  0.1647,  0.1647, 0.1647, 0.0154, 0.1647, 0.1647, 0.0174, 0.0174, 0.0847, 0.1647,  0.0162],
        [ 2,  0.2560, 0.2560, 0.2560,  0.2560,  0.2560,  0.2560, 0.2560, 0.2560, 0.2560, 0.2560, 0.0716, 0.1104, 0.0484, 0.0257,  0.1669],
        [ 3,  0.5   , 0.1863, 0.1863,  0.1863,  0.1863,  0.1863, 0.1863, 0.1863, 0.1863, 0.1863, 0.0813, 0.6590, 0.0090, 0.1635,  0.3045],
        [ 4,  0.1   , 0.1647, 0.1647,  0.1647,  0.1647,  0.1647, 0.1647, 0.1647, 0.1647, 0.1647, 0.0174, 0.7731, 0.1647, 0.1647,  0.0162],
        [ 5,  0.1   , 0.1647, 0.1647,  0.1647,  0.1647,  0.1647, 0.1647, 0.1647, 0.1647, 0.1647, 0.0174, 0.7731, 0.1647, 0.1647,  0.0162],
        [ 6,  0.1   , 0.1647, 0.1647,  0.1647,  0.1647,  0.1647, 0.1647, 0.1647, 0.1647, 0.1647, 0.0174, 0.7731, 0.1647, 0.1647,  0.0162],
        [ 7,  0.1   , 0.1647, 0.1647,  0.1647,  0.1647,  0.1647, 0.1647, 0.1647, 0.1647, 0.1647, 0.0174, 0.7731, 0.1647, 0.1647,  0.0162],
        [ 8,  0.1   , 0.1647, 0.1647,  0.1647,  0.1647,  0.1647, 0.1647, 0.1647, 0.1647, 0.1647, 0.0174, 0.7731, 0.1647, 0.1647,  0.0162],
        [ 9,  0.1   , 0.1647, 0.1647,  0.1647,  0.1647,  0.1647, 0.1647, 0.1647, 0.1647, 0.1647, 0.0174, 0.7731, 0.1647, 0.1647,  0.0162],
        [10,  0.1   , 0.1647, 0.1647,  0.1647,  0.1647,  0.1647, 0.1647, 0.1647, 0.1647, 0.1647, 0.0174, 0.7731, 0.1647, 0.1647,  0.0162]
    ]),
    ('fsh_beta_Gz' , [
         [ 1, 0.1425, 0.1088, 0.1246, 0.2116, 0.2116, 0.1647],
         [ 2, 0.1425, 0.1088, 0.1246, 0.1253, 0.1253, 0.0366],
         [ 3, 0.0157, 0.0229, 0.0284, 0.2931, 0.2931, 0.0481],
         [ 4, 0.1375, 0.11  , 0.11  , 0.13  , 0.13  , 0.1647],
         [ 5, 0.1375, 0.11  , 0.11  , 0.13  , 0.13  , 0.1647],
         [ 6, 0.1375, 0.11  , 0.11  , 0.13  , 0.13  , 0.1647],
         [ 7, 0.1375, 0.11  , 0.11  , 0.13  , 0.13  , 0.1647],
         [ 8, 0.1375, 0.11  , 0.11  , 0.13  , 0.13  , 0.1647],
         [ 9, 0.1375, 0.11  , 0.11  , 0.13  , 0.13  , 0.1647],
         [10, 0.1375, 0.11  , 0.11  , 0.13  , 0.13  , 0.1647]
    ]),
    # Catch
    #  Rows:  1) CP_PLCK_TWL  2) CP_PCOD_TWL  3) CP_PCOD_HAL 4)  CP_PCOD_POT
    #         5) CP_OTHR_TWL  6) CP_OTHR_HAL  7) CP_OTHR_POT 8)  CV_PLCK_TWL
    #         9) CV_PCOD_TWL 10) CV_PCOD_HAL 11) CV_PCOD_POT 12) CV_OTHR_TWL
    #        13) CV_OTHR_HAL 14) CV_OTHR_POT 15) HERR_GILL   16) HERR_SEINE
    #  Cols: GEAR POL  COD  ATF  HER  CAP  EUL  SAN  MYC  SAL SAL2  SHR  SQU  EPI  CRA  OTH
    ('fsh_catch_sel' , [
        [ 1,   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1],
        [ 2,   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1],
        [ 3,   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1],
        [ 4,   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1],
        [ 5,   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1],
        [ 6,   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1],
        [ 7,   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1],
        [ 8,   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1],
        [ 9,   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1],
        [10,   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1],
        [11,   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1],
        [12,   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1],
        [13,   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1],
        [14,   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1],
        [15,   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1],
        [16,   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1]
      ]),
      ('fsh_catch_01' , [
        [ 1,  26,   30,   23,   19,   19,   19,   19,   19,   19,   19,   19,   19,   19,   19,   19],
        [ 2,  26,   30,   23,   19,   19,   19,   19,   19,   19,   19,   19,   19,   19,   19,   19],
        [ 3,  44,   40,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20],
        [ 4,  47,   46,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20],
        [ 5,  26,   30,   23,   19,   19,   19,   19,   19,   19,   19,   19,   19,   19,   19,   19],
        [ 6,  44,   40,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20],
        [ 7,  47,   46,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20],
        [ 8,  30,   33,   23,   23,   23,   23,   23,   23,   23,   23,   23,   23,   23,   23,   23],
        [ 9,  30,   33,   23,   23,   23,   23,   23,   23,   23,   23,   23,   23,   23,   23,   23],
        [10,  44,   43,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20],
        [11,  50,   45,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20],
        [12,  30,   33,   23,   23,   23,   23,   23,   23,   23,   23,   23,   23,   23,   23,   23],
        [13,  44,   43,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20],
        [14,  50,   45,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20],
        [15,  20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20],
        [16,  20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20]
      ]),
      ('fsh_catch_99' , [
        [ 1,  61,  108,   75,   36,   36,   36,   36,   36,   36,   36,   36,   36,   36,   36,   36],
        [ 2,  61,  108,   75,   36,   36,   36,   36,   36,   36,   36,   36,   36,   36,   36,   36],
        [ 3,  83,  101,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60],
        [ 4,  86,  106,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60],
        [ 5,  61,  108,   75,   36,   36,   36,   36,   36,   36,   36,   36,   36,   36,   36,   36],
        [ 6,  83,  101,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60],
        [ 7,  86,  106,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60],
        [ 8,  63,   70,   84,   35,   35,   35,   35,   35,   35,   35,   35,   35,   35,   35,   35],
        [ 9,  63,   70,   84,   35,   35,   35,   35,   35,   35,   35,   35,   35,   35,   35,   35],
        [10,  83,  101,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60],
        [11,  84,  101,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60],
        [12,  63,   70,   84,   35,   35,   35,   35,   35,   35,   35,   35,   35,   35,   35,   35],
        [13,  83,  101,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60],
        [14,  84,  101,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60],
        [15,  60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60],
        [16,  60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60]
      ])
    ))
    # In some of these matrices, float values are typed as integers... quick
    # correction of that

    checkfloats(d)

    dnpz = bestnpz()

    if flag:
        return dict(merge_dicts(dnpz, d))
    else:
        return merge_dicts(dnpz, d)

def checkfloats(x):
    """
    Within dict, change lists of mixed floats and integers to all floats
    """
    for ky in x.keys():
        if isinstance(x[ky], list):
            if isinstance(x[ky][0], list):
                x[ky] = [list2floats(i) for i in x[ky]]
            else:
                x[ky] = list2floats(x[ky])
        elif isinstance(x[ky], dict):
            x[ky] = checkfloats(x[ky])
    return x


def list2floats(x):
    """
    Change list of mixed floats and integers to all floats
    """
    if isinstance(x[0], (float,int)):
        if not (all(isinstance(y, float) for y in x) or
                all(isinstance(y, int) for y in x)):
            x = [float(i) for i in x]
    return x


def merge_dicts(*dict_args):
    '''
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    '''
    result = OrderedDict()
    for dictionary in dict_args:
        result.update(dictionary)
    return result

def ice():
    """
    Populate dictionary with ice model parameters

    Stores parameters related to the ice model.  These parameters are input
    through the file indicated by the IPARNAM variable.

    Returns:
        d: ROMS parameter dictionary.  Keys correspond to ROMS variables.
    """
    d = OrderedDict((
    # Logical flag to turn on ice
    ('Lice', True),
    # Time step for sea ice model [s]
    ('DTICE', 600.0),
    # Time step for updating (equilibrium) ice drift in sea ice model [s]
    ('DTICE_EQ', 10800.0),
    # Maximum number of iterations to achieve convergence of the nonlinear
    # ice dynamics solution
    ('nstrs', 1),
    # Number of elastic time steps per viscous step in EVP dynamics
    ('nevp', 60),
    # Density of sea ice [kg/m3]
    ('rhoice', 900.0),
    # Ice-water bulk drag coefficient (nondimensional).
    ('cdiw', 7.5e-3),
    # Air-ice bulk drag coefficient (nondimensional).
    ('cdai', 3.0e-3),
    # Slipperiness parameter for ice: 1.0 (free slip) or -1.0 (no slip)
    ('GAMMA2', -1.0),
    # Density of air [kg/m3]
    ('rho_air', 1.267),
    # Density of dry snow [kg/m3]
    ('rhosnow_dry', 330.0),
    # Density of wet snow [kg/m3]
    ('rhosnow_wet', 450.0),
    # Ice strength coefficient [N/m2]
    ('pstar', 5.0e03),
    # Ice strength exponential weighting coefficient on concentration
    # (nondimensional)
    ('astren', 20.0),
    # Maximum shear strength [N/m2]
    ('zetamax', 1.0e14),
    # Minimum shear strength [N/m2]
    ('zetamin', 0.0),
    # Ellipticity squared of yield curve
    # (nondimensional)
    ('ellip_sq', 4.0),
    # Mohr-Coulomb stress angle [degrees]
    ('alphai', 45.0),
    # Tolerance against zero-divide
    ('tol', 1.e-16),
    # Minimum average ice thickness [m]
    ('min_h', 0.00),
    # Minimum ice concentration
    # (nondimensional)
    ('min_a', 0.00),
    # Maximum ice concentration
    # (nondimensional)
    ('max_a', 1.0),
    # Ice-water turning angle [degrees]
    ('stressang', 0.0),
    # Ice emissivity
    # (nondimensional)
    ('ice_emiss', 0.97),
    # Specific heat of air
    ('spec_heat_air', 1004.0),
    # Transfer coefficient
    # (nondimensional)
    ('trans_coeff', 1.75e-3),
    # Latent heat of sublimation
    ('sublim_latent_heat', 2.834e6),
    # Zero degree Celsius in Kelvin
    ('t0deg', 273.15)
    ))
    return d

def stations():
    """
    Populate dictionary with station output parameters

    Stores parameters related to the stations output.  These parameters are
    input through the file indicated by the SPOSNAM variable.

    Returns:
        d: ROMS parameter dictionary.  Keys correspond to ROMS variables.
    """

    d = OrderedDict((
    # Logical flag to turn on writing station data
    ('Lstations', True),
    # Activate writing of fields
    ('Sout', OrderedDict((
        ('idFsur', True),             # free-surface
        ('idUbar', True),             # 2D U-velocity
        ('idVbar', True),             # 2D V-velocity
        ('idUvel', True),             # 3D U-velocity
        ('idVvel', True),             # 3D V-velocity
        ('idWvel', True),             # 3D W-velocity
        ('idOvel', True),             # 3D omega vertical velocity
        ('idTvar', 17*[True]),        # all (NT) tracers
        ('idUsms', True),             # surface U-stress
        ('idVsms', True),             # surface V-stress
        ('idUbms', True),             # bottom U-stress
        ('idVbms', True),             # bottom V-stress
        ('idUbws', False),            # bottom U-wave stress
        ('idVbws', False),            # bottom V-wave stress
        ('idUbed', False),            # bed wave orbital U-velocity
        ('idVbed', False),            # bed wave orbital V-velocity
        ('idUbot', False),            # bottom U-momentum above bed
        ('idVbot', False),            # bottom V-momentum above bed
        ('idTsur', [True, True]),     # surface net heat and salt flux
        ('idLhea', True),             # latent heat flux
        ('idShea', True),             # sensible heat flux
        ('idLrad', True),             # longwave radiation flux
        ('idSrad', True),             # shortwave radiation flux
        ('idevap', False),            # evaporation rate
        ('idrain', False),            # precipitation rate
        ('idDano', True),             # density anomaly
        ('idVvis', True),             # vertical viscosity
        ('idTdif', True),             # vertical T-diffusion
        ('idSdif', True),             # vertical Salinity diffusion
        ('idHsbl', True),             # depth of surface boundary layer
        ('idHbbl', True),             # depth of bottom boundary layer
        ('idMtke', False),            # turbulent kinetic energy
        ('idMtls', False),            # turbulent length scale
        ('idUice'  , True),           # Ice prognostic variables
        ('idVice'  , True),
        ('idAice'  , True),
        ('idHice'  , True),
        ('idHice'  , True),
        ('idHsnbc' , True),
        ('idTice'  , True),
        ('idHsno'  , True),
        ('idTimid' , True),
        ('idSfwat' , True),
        ('idTauiw' , False),
        ('idChuiw' , False),
        ('idAgeice', True),
        ('idSig11' , True),
        ('idSig12' , True),
        ('idSig22' , True),
        ('idWfr'   , True),
        ('idS0mk'  , True),
        ('idT0mk'  , True),
        ('idTSvar', 14*[True]),      # Stationary, Biology
        ('idTS2var', 8*[True]),
        ('idPT3var', 10*[True]),     # PROD3
        ('idPT2var', 3*[True]),      # PROD2
        ('idBvar', [True, True]),    # Benthic variables
        ('idIcePhL', True),          # Ice Biology variables
        ('idIceNO3', True),
        ('idIceNH4', True),
        ('idBott', 16*[False])       # exposed sediment layer properties
        ))),
    ('NSTATION', 3),
    ('POS', {
        'GRID': [1, 1, 1],
        'FLAG': [1, 1, 1],
        'X-POS': [195.9500, 195.0, 181.9],
        'Y-POS': [56.8636, 63.855, 56.597]
    })
    ))
    return d



