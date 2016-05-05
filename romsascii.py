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
#from decimal import Decimal
import numpy as np
import re




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
    """
    Formats input boolean as string 'T' or 'F'
    """
    if not isinstance(x, bool):
        return x
    y = '{}'.format(x)[0]
    return y
    

def float2str(x):
    """
    Formats input float as Fortran-style double-precision string
    """
    if not isinstance(x, float):
        return
    y = '{}'.format(x).replace('e','d')
    if not 'd' in y:
        y = '{}d0'.format(y)
    return y


def consecutive(data, stepsize=0):
    """
    Groups values in list based on difference between consecutive elements
    
    Args:
        data:       a list
        stepsize:   difference between consecutive elements to use for
                    grouping. Default = 0, i.e. identical values grouped
    
    Returns:
        list of lists, values grouped.
    
    Example:
        consecutive([1, 1, 1, 2, 2, 4, 5]) -> [[1, 1, 1], [2, 2], [4], [5]]
        consecutive([1, 1, 1, 2, 2, 4, 5], 1) -> [[1], [1], [1, 2], [2], [4, 5]]
    """
    data = np.array(data)
    tmp = np.split(data, np.where(np.diff(data) != stepsize)[0]+1)
    tmp = [x.tolist() for x in tmp]
    return tmp

def list2str(tmp):
    """
    Convert list of bools, floats, or integers to string
    """    
    if not (all(isinstance(x, float) for x in tmp) or 
            all(isinstance(x, bool)  for x in tmp) or
            all(isinstance(x, int)   for x in tmp)):  
        return tmp
                
    consec = consecutive(tmp)

    if isinstance(tmp[0], float):
        y = map(lambda x: '{num}*{val}'.format(num=len(x),val=float2str(x[0])), consec)
    elif isinstance(tmp[0], bool):
        y = map(lambda x: '{num}*{val}'.format(num=len(x),val=bool2str(x[0])), consec)
    else:
        y = map(lambda x: '{num}*{val}'.format(num=len(x),val=x[0]), consec)
    y = ' '.join(y)
    y = re.sub('\s1\*', ' ', y)
    y = re.sub('^1\*', '', y)
    return y    

def checkforstring(x, prefix=''):
    """
    Check that all dictionary entries have been stringified
    """
    for ky in x.keys():
        if isinstance(x[ky], dict):
            checkforstring(x[ky], ky)
        else:
            if not isinstance(x[ky], (str)):
                print('{}{}'.format(prefix, ky))
    

def formatforascii(d):
    """
    Formats values in ROMS parameter dictionary for ascii input files
    
    This function converts values to the Fortran-ish syntax used in ROMS ascii
    input files.  Floats are converted to the double-precision format of
    Fortran read/write statements, booleans are converted to T/F, and lists are
    converted to space-delimited strings of the above (compressed using * for
    repeated values where applicable)
    
    Args:
        d:  ROMS parameter dictionary
    
    Returns:
        dictionary with the same keys as the input dictionary, but with values
        replaced by the new strings (or integers) that will be used in printing
        to file.
    
    """
    newdict = copy.copy(d)
    for x in newdict:
        if isinstance(newdict[x], float):
            newdict[x] = float2str(newdict[x])
        elif isinstance(newdict[x], bool):
            newdict[x] = bool2str(newdict[x])
        elif isinstance(newdict[x], int):
            newdict[x] = '{}'.format(newdict[x])
        elif isinstance(newdict[x], list):                  
            tmp = newdict[x]
            if isinstance(tmp[0], list):
                newdict[x] = map(float2str, tmp)
            else:
                newdict[x] = list2str(tmp)                
        elif isinstance(newdict[x], dict):
            newdict[x] = formatforascii(newdict[x])
    
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



    














