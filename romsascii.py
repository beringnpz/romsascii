# -*- coding: utf-8 -*-
"""
ROMS utility function module

This module provides a few utility functions to assist in running ROMS
simulations.  It is primarily focused on setting up the primary ROMS 
ocean.in input file, and on the parsing of standard output from a ROMS
simulation.

At the moment, this is tailored to the Bering 10K domain simulations.

"""

import datetime
import os.path
import subprocess
from string import Template
import copy
from decimal import Decimal
import numpy as np
import re

# The default Bering 10K ocean.in dictionary
# Most parameters are set to the values we use for Bering  10K runs.  However, 
# all file names are set to defaults and need to be changed. 
def defaultparams():
    """
    Populate dictionary with default ROMS parameters
    
    Returns:
        d: ROMS parameter dictionary.  Keys correspond to ROMS variables.
    """
    nat = 2
    npt = 167
    
    d = {
        # Application title
        'TITLE': 'Bering Sea 10 km Grid',
        # C-preprocessing Flag
        'MyAppCPP': 'NEP5',
        # Input variable information file name.
        'VARNAME': 'varinfo.dat',
        # Grid dimension parameters.
        'Lm': 180,
        'Mm': 256,
        'N':   10,
        'Nbed': 0,
        'NAT': nat,
        'NPT': npt,
        'NCS': 0,
        'NNS': 0,
        # Domain decomposition parameters
        'NtileI': 7,
        'NtileJ': 20,
        # Time-Stepping parameters
        'NTIMES': 0,
        'DT': 600,
        'NDTFAST': 40,
        # Model iteration loops parameters.
        'ERstr': 1,
        'ERend': 1,
        'Nouter': 1,
        'Ninner': 1,
        'Nintervals': 1,
        # Number of eigenvalues and eigenvectors for GST
        'NEV': 2,
        'NCV': 10,
        # Input/Output parameters.
        'NRREC': 0,
        'LcycleRST': True,
        'NRST': 1008,
        'NSTA': 6,
        'NFLT': 144,
        'NINFO': 1,
        # Output history, average, diagnostic files parameters.
        'LDEFOUT': True,
        'NHIS': 1008,
        'NDEFHIS': 10080,
        'NTSAVG': 1,
        'NAVG': 1008,
        'NDEFAVG': 10080,
        'NTSDIA': 1, 
        'NDIA': 10,
        'NDEFDIA': 0,
        # Output tangent linear and adjoint models parameters.
        'LcycleTLM' : False,
             'NTLM' : 72,
          'NDEFTLM' : 0,
        'LcycleADJ' : False,
             'NADJ' : 72,
          'NDEFADJ' : 0,
        # Output check pointing GST restart parameters.
           'LrstGST' :  False,   
        'MaxIterGST' :  500, 
              'NGST' :  10,  
        # Relative accuracy of the Ritz values computed in the GST analysis.
        'Ritz_tol' :  1.0e-15,
        # Harmonic/biharmonic horizontal diffusion of tracers
        'TNU2' : (nat+npt)*[25.0],
        'TNU4' : (nat+npt)*[0.0], 
        # Harmononic/biharmonic, horizontal viscosity coefficient
        'VISC2' : 25.0,  
        'VISC4' : 0.0,   
        # Vertical mixing coefficients for active tracers
        'AKT_BAK' : (nat+npt)*[1.0e-6], 
        # Vertical mixing coefficient for momentum
        'AKV_BAK' : 1.0e-5,  
        # Turbulent closure parameters
        'AKK_BAK' : 5.0e-6, 
        'AKP_BAK' : 5.0e-6, 
         'TKENU2' : 0.0,  
         'TKENU4' : 0.0,
        # Generic length-scale turbulence closure parameters
            'GLS_P' : -1.0,  
            'GLS_M' : 0.5,
            'GLS_N' : -1.0,
         'GLS_Kmin' : 7.6e-6,
         'GLS_Pmin' : 1.0e-12,
         'GLS_CMU0' : 0.5477,
           'GLS_C1' : 0.555,
           'GLS_C2' : 0.833,
          'GLS_C3M' : -0.6,
          'GLS_C3P' : 1.0,
         'GLS_SIGK' : 2.0,
         'GLS_SIGP' : 2.0,
        # Constants used in momentum stress computation
          'RDRG' : 3.0e-04,
         'RDRG2' : 3.0e-03,
           'Zob' : 0.02, 
           'Zos' : 0.02, 
        # Height (m) of atmospheric measurements for Bulk fluxes parameterization
        'BLK_ZQ' : 10.0,
        'BLK_ZT' : 10.0,
        'BLK_ZW' : 10.0,
        # Minimum depth for wetting and drying
        'DCRIT' :0.50,
        # Various parameters
          'WTYPE' : 5,
        'LEVSFRC' : 15,
        'LEVBFRC' : 1,
        # Vertical S-coordinates parameters
        'THETA_S' : 5.0,  
        'THETA_B' : 0.4,  
         'TCLINE' : 10.0, 
        # Mean Density and Brunt-Vaisala frequency
            'RHO0' :  1025.0, 
         'BVF_BAK' :  1.0e-4,   
        # Time-stamps
        'DSTART': '',
        'TIDE_START': -693962.0,
        'TIME_REF': 19000101.0,
        # Nudging/relaxation time scales
         'TNUDG' : 360.0, 
         'ZNUDG' : 360.0,
        'M2NUDG' : 360.0,
        'M3NUDG' : 360.0,
        # Factor between passive and active open boundary condition
        'OBCFAC' : 120.0,
        # Linear equation of State parameters
             'R0' : 1027.0,
             'T0' : 10.0,  
             'S0' : 35.0,  
          'TCOEF' : 1.7e-4,   
          'SCOEF' : 7.6e-4,   
        # Slipperiness parameter
        'GAMMA2' : 1.0,
        # Starting  and ending day for adjoint sensitivity forcing
        'DstrS' : 0.0,
        'DendS' : 0.0,
        # Starting and ending vertical levels of the 3D adjoint state variables
        'KstrS' : 1, 
        'KendS' : 1,
        # Specify the adjoint variables whose sensitivity is required
        'Lstate' : {
            'isFsur' : False,
            'isUbar' : False,
            'isVbar' : False,
            'isUvel' : False,
            'isVvel' : False,
            'isTvar' : [False, False]
        },
        # Stochastic optimals time decorrelation scale
        'SO_decay' : 2.0, 
        # Specify the surface forcing variables whose stochastic optimals are required.
        'SOstate' : {
            'isUstr' : True,
            'isVstr' : True,
            'isTsur' : [False, False]
        },
        # Stochastic optimals surface forcing standard deviation for dimensionalization.
        'SO_sdev': {
            'isUstr' : 1.0,      
            'isVstr' : 1.0,      
            'isTsur' : [1.0, 1.0]
        },
        # Activate writing of fields into HISTORY output file.
        'Hout': {
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
        },
        # Generic User parameters
        'NUSER' : 0,
         'USER' : 0.0,
        # Input NetCDF file names
        'GRDNAME': 'grd.nc',
        'ININAME': 'ini.nc',
        'IRPNAME': 'irp.nc',
        'IADNAME': 'iad.nc',
        'BRYNAME': 'bry.nc',
        'ADSNAME': 'ads.nc',
        # Input forcing NetCDF file name(s).
        'NFFILES': 1,
        'FRCNAME': 'frc.nc',
        # Output NetCDF file names
        'RSTNAME': 'rst.nc',
        'HISNAME': 'his.nc',
        'AVGNAME': 'avg.nc',
        'STANAME': 'sta.nc',
        'FLTNAME': 'flt.nc',
        'GSTNAME': 'gst.nc',
        'TLMNAME': 'tlm.nc',
        'TLFNAME': 'tlf.nc',
        'ADJNAME': 'adj.nc',
        # Input ASCII parameter filenames
        'APARNAM': 'assimilation.in',
        'SPOSNAM': 'stations.in',
        'FPOSNAM': 'floats.in',
        'IPARNAM': 'ice.in',
        'BPARNAM': 'bio.in',
        'SPARNAM': 'sediment.in',
        'USRNAME': 'myfile.dat'
    }    
    return d    
    
# Fill in the appropriate initial, forcing, and boundary files based on either
# CORE or CFSR data
    
def ini_frc_bry(d, indir, group):
    """
    Set up input files to use either CORE or CFSR sets
    
    This function alters the ININAME, BRYNAME, FRCNAME, and NFFILE values 
    in a ROMS parameter dictionary to match predefined sets for CORE or 
    CFSR data.
    
    Args:
        d:      ROMS parameter dictionary
        indir:  folder where input files are stored
        group:  string, 'core' or 'cfsr'
    
    Returns:
        d:      ROMS parameter dictionary
    """
    
    same = ('tides_OTBS.nc',
            'riverrunoff-usgs-fsu.1950-2015.efol20km.nc',
            'sss.clim.nc'
            )
    
    core = ('Pair.1948-2006.Bering.nc',
            'Qair.1948-2006.Bering.nc',
            'tair_all.Bering.nc',
            'lwrad.1948-2006.Bering.nc',
            'swrad.1948-2006.Bering.nc',
            'Uwind.1948-2006.Bering.nc',
            'Vwind.1948-2006.Bering.nc',
            'rain.1948-2006.Bering.nc'
            )
            
    cfsr = ('roms-cfs-atmos-Pair-regridded-2002-2014.nc',
            'roms-cfs-atmos-Qair-regridded-2002-2014.nc',
            'roms-cfs-atmos-Tair-regridded-2002-2014.nc',
            'roms-cfs-atmos-lwrad-reduced-by-3-percent-regridded-2002-2014.nc',
            'roms-cfs-atmos-swradave-reduced-by-10-percent-regridded-2002-2014.nc',
            'roms-cfs-atmos-Uwind-regridded-2002-2014.nc',
            'roms-cfs-atmos-Vwind-regridded-2002-2014.nc',
            'roms-cfs-atmos-rain-regridded-2002-2014.nc'
            )
            
    if group == 'core':
        frc = core + same
        bry = 'core_bry.1969-2004.nc'
        ini = 'ini.1969.nc'
    elif group == 'cfsr':
        frc = cfsr + same
        bry = 'cfsr_bc_2002-2014_5day_bio_10layer.nc'
        ini = 'ini.2004.nc'
        
    
    frc = list(map(lambda x: os.path.join(*(indir,'frc',x)), frc))      
    d['FRCNAME'] = '\n               '.join(frc)
    d['NFFILES'] = len(frc)
    d['BRYNAME'] = os.path.join(*(indir,'bry',bry))   
    d['ININAME'] = os.path.join(*(indir,'ini',ini)) 
    
    return d    
    
# Fill in all input files except forcing, boundary, and initial    
    
def other_in(d, indir):
    """
    Set input files for typical Bering 10K run
    
    This function alters the VARNAME, GRDNAME, SPOSNAM, BPARNAM, and 
    IPARNAM values in a ROMS parameter dictionary to match predefined 
    sets for our typical Bering 10K runs.
    
    Args:
        d:      ROMS parameter dictionary
        indir:  folder where input files are located
    
    Returns:
        d:      ROMS parameter dictionary
    """
    
    d['VARNAME'] = os.path.join(*(indir,'var','GK_varinfo.dat'))
    d['GRDNAME'] = os.path.join(*(indir,'grd','Bering_grid_withFeast.nc'))
    d['SPOSNAM'] = os.path.join(*(indir,'spos','stations_bering_10k.in'))
    d['BPARNAM'] = os.path.join(*(indir,'bpar','sebs_bio_ajh_08_26_11.in'))
    d['IPARNAM'] = os.path.join(*(indir,'ipar','ice.in'))
    
    return d

# Create file based on template but filling in values from d
    
def filltemplate(d, templatefile, outfile):
    """
    Fills in ocean.in template file with dictionary values
    
    This function fills in an ocean.in template file.  The template file
    should have variable placeholders (e.g $var) wherever d['var'] will 
    be filled in.
    
    Args:
        d:              ROMS parameter dictionary
        templatefile:   template file name
        outfile:        name of output file to be created
    
    """    
  
    with open(templatefile, 'r') as fin:
        with open(outfile, 'w') as fout:
            src = Template( fin.read() )
            result = src.substitute(d)
            fout.write(result)
    
def bool2str(x):
    if not isinstance(x, bool):
        return x
    y = '{}'.format(x)[0]
    return y

        
    
def float2str(x):  
    if not isinstance(x, float):
        return 
    y = '{}'.format(x).replace('e','d')
    if not 'd' in y  :
        y = '{}d0'.format(y)
    return y
        

def consecutive(data, stepsize=0):
    data = np.array(data)
    tmp = np.split(data, np.where(np.diff(data) != stepsize)[0]+1)
    tmp = [x.tolist() for x in tmp]
    return tmp
    

def formatforin(d):
    newdict = copy.copy(d)
    for x in newdict:
        if isinstance(newdict[x], float):
            newdict[x] = float2str(newdict[x])
        elif isinstance(newdict[x], bool):
            newdict[x] = bool2str(newdict[x])
        elif isinstance(newdict[x], list):
            tmp = newdict[x]
            consec = consecutive(tmp)
            if isinstance(tmp[0], float):
                y = map(lambda x: '{num}*{val}'.format(num=len(x),val=float2str(x[0])), consec)
            elif isinstance(newdict[x][0], bool):
                y = map(lambda x: '{num}*{val}'.format(num=len(x),val=bool2str(x[0])), consec)
            y = ' '.join(y)
            y = re.sub('\s1\*', ' ', y)
            newdict[x] = y
            
    return newdict

# Set all the time-related variables

def filltimevars(d, tstep, datestart, dateend, tref, dthis, dtavg, dtsta, dtdefhis, dtdefavg):   
    """
    Calculate time-related parameters for ROMS dictionary
    
    This function calculates DT, DSTART, NTIMES, NHIS, NAVG, NSTA, 
    NDEFHIS, and NDEFAVG variables based on desired calendar dates and 
    time intervals.
    
    Args:
        d:          ROMS parameter dictionary
        tstep:      timedelta object, simulation time step
        datestart:  datetime object, date of simulation start (date of 
                    step = 0, actual start of computation set by 
                    initialization file time)
        dateend:    datetime object, date of simulation end
        tref:       datetime object, reference date
        dthis:      timedelta object, time between history file writes
        dtavg:      timedelta object, time between average file writes
        dtsta:      timedelta object, time between station file writes
        dtdefhis:   timedelta object, time between new history file 
                    creation
        dtdefavg:   timedelta object, time between new average file 
                    creation
    """ 
    
    d['DT'] = int(tstep.total_seconds())
    d['DSTART'] = '{:d}.0d0'.format((datestart - tref).days)
    d['NTIMES'] = int((dateend - datestart).total_seconds()/tstep.total_seconds())
    d['NHIS'] = int(dthis.total_seconds()/tstep.total_seconds())
    d['NAVG'] = int(dtavg.total_seconds()/tstep.total_seconds())
    d['NSTA'] = int(dtsta.total_seconds()/tstep.total_seconds())
    d['NDEFHIS'] = int(dtdefhis.total_seconds()/tstep.total_seconds())
    d['NDEFAVG'] = int(dtdefhis.total_seconds()/tstep.total_seconds())
    
    return d
    
def parseromslog(fname):
    """
    Parse ROMS standard output log for some details
    
    This function extracts details about the success (or not) of a ROMS 
    simulation
    
    Args
        fname:  log file name
    
    Returns:
        dictionary object with the following keys:
            cleanrun:   True if simulation ran ran with errors
            blowup:     True if simulation blew up
            laststep:   Index of last step recorded
            lasthis:    Name of last history file defined
    """
    
    with open(fname, 'r') as f:
        lines = f.read()
        lnnum = lines.find('ROMS/TOMS: DONE')
        cleanrun = lnnum != -1
    
        lnnum = lines.find('Blowing-up: Saving latest model state into  RESTART file')    
        blowup = lnnum != -1
    
    step = []
    lasthis = []
    if cleanrun:
        with open(fname, 'r') as f:
        
            datablock = False
    
            for line in f:
                if line.find('STEP   Day HH:MM:SS  KINETIC_ENRG   POTEN_ENRG    TOTAL_ENRG    NET_VOLUME') != -1:
                    datablock = True
                elif line.find('Elapsed CPU time (seconds):') != -1:
                    datablock = False
                elif datablock:
                    tmp = line.split() #string.split(line.strip())
                    if len(tmp) == 7 and tmp[0].isdigit():
                        step = int(tmp[0])  
                    if len(tmp) == 6 and tmp[0] == 'DEF_HIS':
                        lasthis = tmp[-1]     
    
    return {'cleanrun': cleanrun, 'blowup': blowup, 'laststep': step, 'lasthis':lasthis}
    
def runroms(d, templatefile, outdir, logdir, outbase, logbase, mpivars, dryrun=False):
    """
    Run a ROMS simulation
    
    This function runs the ROMS executable with the specified options and 
    inputs
    
    Args:
        d:              ROMS parameter dictionary
        templatefile:   ocean.in template file
        outdir:         path to folder where output files will be placed
        logdir:         path to folder where log files (standard output 
                        and standard error from ROMS) will be placed
        outbase:        string, base for all output file names
        logbase:        string, base for all log file names
        mpivars:        dictionary with running options
                        mpiexe:     path to mpirun executable
                        np:         number of processors to use
                        hostfile:   host file name specifying cores to 
                                    use
                        romsexe:    path to roms executable
                    
    Returns:
        dictionary object with the following keys:
                        log:        path to log file with ROMS standard output
                        err:        path to log file with ROMS standard error
    """

    # Set output file names  

    if not os.path.exists(outdir):
        os.makedirs(outdir) 
        

    d['RSTNAME'] = os.path.join(*(outdir, '{}_rst.nc'.format(outbase)))
    d['HISNAME'] = os.path.join(*(outdir, '{}_his.nc'.format(outbase)))
    d['AVGNAME'] = os.path.join(*(outdir, '{}_avg.nc'.format(outbase)))
    d['STANAME'] = os.path.join(*(outdir, '{}_sta.nc'.format(outbase)))
    d['FLTNAME'] = os.path.join(*(outdir, '{}_flt.nc'.format(outbase)))
        
    filltemplate(d, templatefile, 'ocean.tmp.in')
    
    if not os.path.exists(logdir):
        os.makedirs(logdir)     
    
    logfile = os.path.join(*(logdir, 'out_{}.txt'.format(logbase)))
    errfile = os.path.join(*(logdir, 'err_{}.txt'.format(logbase))) 
    
    cmd = [mpivars['mpiexe'], 
               '-np', str(mpivars['np']), 
               '--hostfile', mpivars['hostfile'], 
               mpivars['romsexe'], 
               'ocean.tmp.in']
    if dryrun:
        print('{} >{} 2>{}'.format(' '.join(cmd), logfile, errfile))
    else:
        with open(logfile, 'w') as fout, open(errfile, 'w') as ferr:
            subprocess.run(cmd, stdout=fout, stderr=ferr)
        
    return {'log': logfile, 'err': errfile}

    
def runromsthroughblowup(d, templatefile, logdir, outdir, outbase, faststep, slowstep, timevars, mpivars):
    """
    Run a ROMS simulation, attempting to get past blowups
    
    This function runs the ROMS executable with the specified options and 
    inputs.  If the simulation flows up, it decreases the simulation step
    size and runs for a month before returning to the original step size.
    
    Args:
        d:              ROMS parameter dictionary
        templatefile:   ocean.in template file
        logdir:         path to folder where log files (standard output 
                        and standard error from ROMS) will be placed
        outdir:         path to folder where output files will be placed
        outbase:        string, base used for both output and log files 
                        (will be modified with counters indicating 
                        restarts)
        faststep:       timedelta object, simulation time step for 
                        initial run
        slowstep:       timedelta object, smaller time step used to get 
                        past blowups
        timevars:       dictionary with desired time-related variables
                        datestart:  datetime object, date of simulation 
                                    start (date of step = 0, actual start 
                                    of computation set by initialization 
                                    file time)
                        dateend:    datetime object, date of simulation 
                                    end
                        tref:       datetime object, reference date
                        dthis:      timedelta object, time between 
                                    history file writes
                        dtavg:      timedelta object, time between 
                                    average file writes
                        dtsta:      timedelta object, time between 
                                    station file writes
                        dtdefhis:   timedelta object, time between new 
                                    history file creation
                        dtdefavg:   timedelta object, time between new 
                                    average file creation
        mpivars:        dictionary with running options
                        mpiexe:     path to mpirun executable
                        np:         number of processors to use
                        hostfile:   host file name specifying cores to 
                                    use
                        romsexe:    path to roms executable
    
    """
    
    #---------------------
    # Run initial attempt    
    #---------------------  
    
    count = 1  
    
    # Start with input time bounds and the fast time step
    
    filltimevars(d, faststep, timevars['datestart'], timevars['dateend'], 
                 timevars['tref'], timevars['dthis'], timevars['dtavg'], 
                 timevars['dtsta'], timevars['dtdefhis'], timevars['dtdefavg'])
                
    # Run ROMS
                
    outbasetmp = '{}_{:02d}'.format(outbase, count)
    logbasetmp = '{}_{:02d}_fast'.format(outbase, count)
    
    print('Running {} (Init)'.format(outbasetmp))
    s = runroms(d, templatefile, outdir, logdir, outbasetmp, logbasetmp, mpivars)
    
    # Parse the log file to make sure ROMS ran cleanly, and to check for a 
    # blowup
    
    r = parseromslog(s['log'])  
    
    if not r['cleanrun']:
        print('ROMS crashed')
        return

    # If it blew up...
    
    while r['blowup']:
        
        #---------------------
        # Run a slow-step run    
        #---------------------

        # Increment the count, for record-keeping
        
        count+=1     
        
        # Restart run, using last good history file for initialization
        
        d['NRREC'] = -1  
        d['ININAME'] = r['lasthis']

        # Calculate new end date, 30 days from blow up point, and change time 
        # variables to use this and the slow timestep        
        
        dateblewup = timevars['datestart'] + datetime.timedelta(seconds=r['laststep']*faststep.total_seconds())
        newend = dateblewup + datetime.timedelta(days=30)
        
        filltimevars(d, slowstep, timevars['datestart'], newend, 
                     timevars['tref'], timevars['dthis'], timevars['dtavg'], 
                     timevars['dtsta'], timevars['dtdefhis'], timevars['dtdefavg'])
                     
        # Run ROMS

        outbasetmp = '{}_{:02d}'.format(outbase, count)
        logbasetmp = '{}_{:02d}_slow'.format(outbase, count)
        
        print('Running {} (slow)'.format(outbasetmp))
        s = runroms(d, templatefile, outdir, logdir, outbasetmp, logbasetmp, mpivars)
    
        # Parse this run's log to make sure it finished cleanly
    
        r = parseromslog(s['log'])  
        if not r['cleanrun']:
            print('ROMS crashed')
            return

        # Switch back to fast time step and try to run until the end

        count+=1     

        d['NRREC'] = -1  
        d['ININAME'] = r['lasthis']      
        
        filltimevars(d, faststep, timevars['datestart'], timevars['dateend'], 
                     timevars['tref'], timevars['dthis'], timevars['dtavg'], 
                     timevars['dtsta'], timevars['dtdefhis'], timevars['dtdefavg'])
                     
        outbasetmp = '{}_{:02d}'.format(outbase, count)
        logbasetmp = '{}_{:02d}_fast'.format(outbase, count)
        
        print('Running {} (fast)'.format(outbasetmp))
        s = runroms(d, templatefile, outdir, logdir, outbasetmp, logbasetmp, mpivars)
            
        # Parse this run's log to check for another blowup
    
        r = parseromslog(s['log'])  
        if not r['cleanrun']:
            print('ROMS crashed')
            return
            
    print('Done')
    



    













