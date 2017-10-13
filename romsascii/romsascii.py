# -*- coding: utf-8 -*-
"""
ROMS utility function module

This module provides a few utility functions to assist in running ROMS
simulations.  It is primarily focused on setting up the primary ROMS
ocean.in input file, and on the parsing of standard output from a ROMS
simulation.

At the moment, this is tailored to the Bering 10K domain simulations.

"""

from datetime import datetime, timedelta
import os
import subprocess
from string import Template
import copy
import numpy as np
import re
import sys
import glob

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
    # if not any(x in y for x in ['d','.']):
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
    
    consec = consecutive(tmp, stepsize=-99999)

    consecstr = [None]*len(consec)

    for ii in range(0, len(consec)):
        n = len(consec[ii]) 
        if isinstance(tmp[0], float):
            sampleval = float2str(consec[ii][0])
        elif isinstance(tmp[0], bool): 
            sampleval = bool2str(consec[ii][0])
        else:
            sampleval = consec[ii][0]
    
        if n > 1:
            consecstr[ii] = '{num}*{val}'.format(num=n,val=sampleval)
        else:
            consecstr[ii] = '{val}'.format(val=sampleval)
        
    y = ' '.join(consecstr)
    return y

def checkforstring(x, prefix=''):
    """
    Check that all dictionary entries have been stringified
    """
    for ky in x.keys():
        if isinstance(x[ky], dict):
            checkforstring(x[ky], ky)
        else:
            
            if not (isinstance(x[ky], (str)) or
                    (isinstance(x[ky],list) and
                    (all(isinstance(i,str) for i in x[ky])))):
                print('{}{}'.format(prefix, ky))


def formatforascii(d):
    """
    Formats values in ROMS parameter dictionary for ascii input files
    
    This function converts values to the Fortran-ish syntax used in ROMS ascii
    input files.  Floats are converted to the double-precision format of
    Fortran read/write statements, booleans are converted to T/F, integers are
    converted straight to strings, and lists are converted to space-delimited
    strings of the above (compressed using * for repeated values where
    applicable).  A few specific lists are converted to the appropriate tables.
    
    Args:
        d:  ROMS parameter dictionary
    
    Returns:
        dictionary with the same keys as the input dictionary, but with values
        replaced by the new strings (or lists or dicts of strings) that will be
        used in printing to file.
    
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
            if x == 'FRCNAME':
                # Forcing files: 1 per line
                newdict[x] = '\n{:18s}'.format('').join(tmp)
            elif isinstance(tmp[0], list):
                newdict[x] = [list2str(i) for i in tmp]
            else:
                newdict[x] = list2str(tmp)
        elif isinstance(newdict[x], dict):
            if x == 'POS':
                # Stations table
                tmp = newdict[x]
                tablestr = ['{:4s} {:4s} {:12s} {:12s} {:12s}'.format('GRID','FLAG', 'X-POS', 'Y-POS', 'COMMENT')]
                for ii in range(len(tmp['GRID'])):
                    if tmp['FLAG'][ii] == 1: # lat, lon grid pairs
                        tablestr.append('{space:17s}{grid:4d} {flag:4d} {xpos:12f} {ypos:12f}'.format(space='',grid=tmp['GRID'][ii], flag=tmp['FLAG'][ii], xpos=tmp['X-POS'][ii], ypos=tmp['Y-POS'][ii]))
                    elif tmp['FLAG'][ii] == 0: # I,J grid pairs
                        tablestr.append('{space:17s}{grid:4d} {flag:4d} {xpos:12d} {ypos:12d}'.format(space='',grid=tmp['GRID'][ii], flag=tmp['FLAG'][ii], xpos=tmp['X-POS'][ii], ypos=tmp['Y-POS'][ii]))
                newdict[x] = '\n'.join(tablestr)
            else:
                newdict[x] = formatforascii(newdict[x])
    
    return newdict

def writeromsascii(d, fname, filetype='phys'):
    """
    Create ROMS ascii file based on values in a dictionary
    
    Assuming the keys in d correspond to ROMS keyword input arguments, this
    function builds a ROMS ascii input file using those keywords and their
    corresponding values.
    
    Args:
        d:          ROMS parameter dictionary
        fname:      name of file to create
    
    Optional keyword arguments:
        filetype:   type of file (used to determine which parameters need
                    = vs ==).  Can be 'phys', 'bio', 'ice', or 'stations'.  If 
                    not included, default is 'phys'.
    
    """
    dstr = formatforascii(d)
    
    with open(fname, 'w') as f:
        for ky in dstr:
            if isinstance(dstr[ky], list):
                for i in dstr[ky]:
                    f.write(formatline(ky, i, filetype))
            elif isinstance(dstr[ky], dict):
                for i in dstr[ky]:
                    newkey = '{}({})'.format(ky,i)
                    f.write(formatline(newkey, dstr[ky][i], filetype))
            else:
                f.write(formatline(ky, dstr[ky], filetype))
                

def formatline(kw, val, filetype='phys'):
    """
    Format an parameter definition line for a ROMS file
    
    ROMS input parameter statements follow the format of
      KEY1 =  VALUE1
      KEY2 == VALUE2
    The single equals sign indicates singular assingment (i.e. no nested 
    grid dependency), while the double equals indicates that the 
    parameter expects Ngrids replicates of the parameter values, one per 
    grid.  
    
    This function formats a parameter-assignment line, comparing the 
    input key to a saved collection of variables that I know expect 
    singular assignment.  Any keys not in that collection will be 
    formatted with a plural assignment.
    
    Args:
        kw:         ROMS input variable name
        val:        stringified version of ROMS variable value(s)
    
    Optional keyword arguments:
        filetype:   type of file (used to determine which parameters need
                    = vs ==).  Can be 'phys', 'bio', 'ice', or 'stations'.  
                    If not included, default is 'phys'.
    
    """
    if filetype == 'phys':
        singular = ['TITLE', 'MyAppCPP', 'VARNAME', 'Nbed', 'NAT', 'NPT', 'NCS',
                    'NNS', 'ERstr', 'ERend', 'Nouter', 'Ninner', 'Nintervals',
                    'NEV', 'NCV', 'NRREC', 'LrstGST', 'MaxIterGST', 'NGST',
                    'Ritz_tol', 'RHO0', 'BVF_BAK', 'DSTART', 'TIDE_START',
                    'TIME_REF', 'NUSER', 'USER', 'APARNAM', 'SPOSNAM', 'IPARNAM',
                    'BPARNAM', 'SPARNAM', 'USRNAME']
    elif filetype == 'bio':
        singular = []
    elif filetype == 'ice':
        singular = ['GAMMA2', 'rho_air', 'tol', 'stressang', 'ice_emiss',
                    'spec_heat_air', 'trans_coeff', 'sublim_latent_heat',
                    't0deg']
    elif filetype == 'stations':
        singular = ['POS']
    else:
        raise Exception("filetype input must be either 'phys', 'bio', 'ice', or 'stations'")
        
    
    if kw in singular:
        return '{:>14s} = {}\n'.format(kw,val)
    else:
        return '{:>14s} == {}\n'.format(kw,val)


# Set all the time-related variables

def filltimevars(d, tstep, datestart, dateend, tref, 
                 dthis=timedelta(weeks=1), dtavg=timedelta(weeks=1), 
                 dtsta=timedelta(weeks=1), dtrst=timedelta(weeks=1), 
                 dtflt=timedelta(weeks=1), dtdefhis=timedelta(weeks=10), 
                 dtdefavg=timedelta(weeks=10)):
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
                 
    Optional keyword arguments:
        dthis:      timedelta object, time between history file writes
                    default = 1 week
        dtavg:      timedelta object, time between average file writes
                    default = 1 week
        dtsta:      timedelta object, time between station file writes
                    default = 1 week
        dtrst:      timedelta object, time between restart file writes
                    default = 1 week
        dtdefhis:   timedelta object, time between new history file
                    creation
                    default = 10 weeks
        dtdefavg:   timedelta object, time between new average file
                    creation
                    default = 10 weeks
    """
    
    d['DT'] = int(tstep.total_seconds())
    d['DSTART'] = '{:d}.0d0'.format((datestart - tref).days)
    d['NTIMES'] = int((dateend - datestart).total_seconds()/tstep.total_seconds())
    d['NHIS'] = int(dthis.total_seconds()/tstep.total_seconds())
    d['NAVG'] = int(dtavg.total_seconds()/tstep.total_seconds())
    d['NSTA'] = int(dtsta.total_seconds()/tstep.total_seconds())
    d['NRST'] = int(dtrst.total_seconds()/tstep.total_seconds())
    d['NDEFHIS'] = int(dtdefhis.total_seconds()/tstep.total_seconds())
    d['NDEFAVG'] = int(dtdefhis.total_seconds()/tstep.total_seconds())
    
    dfrac = (tref - datetime(tref.year, tref.month, tref.day)).total_seconds()/86400.0
    datefloat = float('{year}{month:02d}{day:02d}'.format(year=tref.year, month=tref.month, day=tref.day))
    
    d['TIME_REF'] = datefloat + dfrac
    
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
            cleanrun:   True if simulation ran without errors
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
    
def reportstatus(logoutfile):
    """
    Parse ROMS simulation log and report whether simulation completed.
    
    This function parses a ROMS simulation log file.  If that simulation
    finished cleanly, it prints a simple statement to standard output.  
    If not, it terminates execution of the calling program along with a
    print to standard output.
    """
    s = parseromslog(logoutfile)
    if s['cleanrun']:
        print('Done: completed cleanly')
    else:
        print('Done: crashed, see {}'.format(logoutfile))
        sys.exit()

def runroms(rundata, mpiexe="mpirun", romsexe="oceanM", hostfile="", 
            dryrun=False):
    """
    Run a ROMS simulation (without a system scheduler)
    
    This function runs the ROMS executable with the specified options and
    inputs, assuming mpirun can be called directly on the computer system 
    (i.e. no system scheduler for batch job submission)
    
    Args:
        rundata:        dictionary of ROMS simulation info (see output 
                        from createinputfiles)
    
    Optional keyword arguments:
    
        mpiexe:         location of mpirun executable
                        default = 'mpirun'
        romsexe:        location of ROMS compiled executable
                        default = 'oceanM'
        hostfile:       path to host file (telling mpirun which 
                        processors to use
        dryrun:         true to run in dry run mode (prints mpirun 
                        subprocess command to screen), false to actually 
                        run it
                        default = False
    
    Returns:
        dictionary object with the following keys:
        log:            path to log file with ROMS standard output
        err:            path to log file with ROMS standard error
    """
    
    # Set up main command
    
    if hostfile:
        cmd = [mpiexe, 
               '-np', str(rundata['np']),
               '--hostfile', hostfile,
               romsexe, 
               rundata['in']
           ]
    else:
        cmd = [mpiexe, 
               '-np', str(rundata['np']),
               romsexe, 
               rundata['in']
           ]
                  
    # Run ROMS via subprocess call (or print equivalent command to screen):
    
    if dryrun:
        print('{} >{} 2>{}'.format(' '.join(cmd), rundata['out'], rundata['err']))
    else:
        with open(rundata['out'], 'w') as fout, open(rundata['err'], 'w') as ferr:
            subprocess.run(cmd, stdout=fout, stderr=ferr)
            
    return {'log': rundata['out'], 'err': rundata['err']}
    
def createinputfiles(d, outbase, logbase, outdir='.', logdir='.', 
        indir = '.', bio={}, ice={}, stations={}, 
        oceanfile='ocean.tmp.in', bparfile='bio.tmp.in', 
        iparfile='ice.tmp.in', sposfile='stations.tmp.in'):
    """
    Create ascii input files for a ROMS simulation
        
    This function creates all the files necessary to run a ROMS 
    simulation.  At minimum this includes the standard input file (i.e.
    ocean.in).  It will also create the appropriate biology, ice, and 
    stations files as needed, and set up names for files to catch 
    standard output and standard error.
        
    Args:
        d:              ROMS parameter dictionary
        outbase:        string, base for all output file names
        logbase:        string, base for all log file names
    
    Optional keyword arguments:
        outdir:         path to folder where output files will be placed
                        default = current directory
        logdir:         path to folder where log files (standard output
                        and standard error from ROMS) will be placed
                        default = current directory
        indir:          path to folder where any dynamically-generated 
                        input files are saved
                        default = current directory
        bio:            dictionary of biological parameters.  If not 
                        empty, the BPARNAM file will be dynamically 
                        generated based on the values in this dictionary.
                        default = {}
        ice:            dictionary of ice parameters.  If not empty, the 
                        IPARNAM file will be dynamically generated based 
                        on the values in this dictionary.
                        default = {}
        stations:       dictionary of station parameters.  If not empty, 
                        the SPOSNAM file will be dynamically generated 
                        based on the values in this dictionary.
                        default = {}
        oceanfile:      filename for ROMS standard input file
                        default = 'ocean.tmp.in'
        bparfile:       filename for dynamically-generated biological 
                        parameter file (if bio passed as input)
                        default = 'bio.tmp.in'
        iparfile:       filename for dynamically-generated ice parameter 
                        file (if ice passed as input)
                        default = 'ice.tmp.in'
        sposfile:       filename for dynamically-generated stations 
                        parameter file (if stations passed as input)
                        default = 'stations.tmp.in'
        
    Returns:
        dictionary object with the following keys:
        in:             path to ROMS standard input file
        out:            path to ROMS standard output file
        err:            path to ROMS standard error file
        np:             number of processors required to run this 
                        simulation in parallel
        
    """
    # Set output file names
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    d['RSTNAME'] = os.path.join(*(outdir, '{}_rst.nc'.format(outbase)))
    d['HISNAME'] = os.path.join(*(outdir, '{}_his.nc'.format(outbase)))
    d['AVGNAME'] = os.path.join(*(outdir, '{}_avg.nc'.format(outbase)))
    d['STANAME'] = os.path.join(*(outdir, '{}_sta.nc'.format(outbase)))
    d['FLTNAME'] = os.path.join(*(outdir, '{}_flt.nc'.format(outbase)))
    
    # Create ascii input files
    
    if not os.path.exists(indir):
        os.makedirs(indir)    
    
    if bio:
        bparfullfile =  os.path.join(*(indir, bparfile))
        writeromsascii(bio, bparfullfile, filetype='bio')
        d['BPARNAM'] = bparfullfile
    
    if ice:
        iparfullfile = os.path.join(*(indir, iparfile))
        writeromsascii(ice, iparfullfile, filetype='ice')
        d['IPARNAM'] = iparfullfile
        
    if stations:
        sposfullfile = os.path.join(*(indir, sposfile))
        writeromsascii(stations, sposfullfile, filetype='stations') 
        d['SPOSNAM'] = sposfullfile   
    
    oceanfullfile = os.path.join(*(indir, oceanfile))
    writeromsascii(d, oceanfullfile, filetype='phys')
    
    # Set up log file names
    
    if not os.path.exists(logdir):
        os.makedirs(logdir)
    
    logfile = os.path.join(*(logdir, 'out_{}.txt'.format(logbase)))
    errfile = os.path.join(*(logdir, 'err_{}.txt'.format(logbase)))
        
    # Return info necessary to run this ROMS simulation
    
    return {'in': oceanfullfile, 'out': logfile, 'err': errfile, 
            'np': d['NtileI'] * d['NtileJ']}
    
def parserst(filebase):
    """
    Parse restart counters from ROMS simulation restart files
    
    This function finds the name of, and parses the simulation counter, 
    from a series of ROMS restart files.  It assumes that those files 
    were using the naming scheme from runromssmart, i.e.
        filebase_XX_rst.nc
    where XX is the counter for number of restarts.
    
    Args:
        filebase:   base name for restart files (can include full path)
    
    Returns:
        d:          dictionary object with the following keys:
                    lastfile:   full path to last restart file
                    cnt:        restart counter of last file incremented 
                                by 1 (i.e. count you would want to 
                                restart with in runromssmart)
    """
    allrst = sorted(glob.glob(os.path.join(filebase + "*rst*.nc")))
    if len(allrst) == 0:
        rst = []
        cnt = 1
    else:
        rst = allrst[-1]

        pattern = filebase + "_(\d+)_rst.nc"
        m = re.search(pattern, rst)
        cnt = int(m.group(1)) + 1
    
    return {'lastfile': rst, 'count': cnt}


def runromssmart(d, outbase, timevars, mpivars, logdir='.', outdir='.',
                indir='.', faststep=[], slowstep=[], count=1, bio={}, 
                ice={}, stations={}, dryrun=False):
    """
    Run a ROMS simulation, attempting to deal with crashes
    
    This function runs the ROMS executable with the specified options and
    inputs.  Before running, it checks to see if a restart file for the 
    simulation already exists, and if so, resumes from there. If the 
    simulation blows up, it decreases the simulation step size and runs 
    for a month before returning to the original step size.
    
    Args:
        d:              ROMS parameter dictionary
        outbase:        string, base name used for input, output, and log 
                        files (will be modified with counters indicating
                        restarts where appropriate)
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
                        hostfile:   host file name specifying cores to
                                    use
                        romsexe:    path to roms executable
    
    Optional keyword arguments:
        logdir:         path to folder where log files (standard output
                        and standard error from ROMS) will be placed
                        default = current directory
        outdir:         path to folder where output files will be placed
                        default = current directory
        indir:          path to folder where any dynamically-generated 
                        input files are saved
        faststep:       timedelta object, simulation time step for
                        initial run.  If empty, the value from d['DT'] 
                        will be used.
                        default = []
        slowstep:       timedelta object, smaller time step used to get
                        past blowups.  If empty, will be half the faststep
                        value.
                        default = []
        count:          index assigned to initial run, assuming no 
                        restart files are found
                        default = 1
        bio:            dictionary of biological parameters.  If not empty, the
                        BPARNAM file will be dynamically generated based on the
                        values in this dictionary.
                        default = {}
        ice:            dictionary of ice parameters.  If not empty, the 
                        IPARNAM file will be dynamically generated based on the 
                        values in this dictionary.
                        default = {}
        stations:       dictionary of station parameters.  If not empty, the 
                        SPOSNAM file will be dynamically generated based on the 
                        values in this dictionary.
                
    Returns:
        cleanexit:      logical scalar, True if ROMS exited the 
                        simulation cleanly (this can include a blowup, as
                        long as the code made it to "ROMS/TOMS Done")
        
    """
    
    #---------------------
    # Run initial attempt
    #---------------------
    
    # Start with input time bounds and the fast time step
    
    filltimevars(d, **timevars)
    if faststep:
        d['DT'] = int(faststep.total_seconds())
    else:
        faststep = timedelta(seconds=d['DT'])
    
    if not slowstep:
        slowstep = faststep/2
    
    # Create ascii input files that will be reused across blowups
    # (No need to keep making these every time)
    
    if not os.path.exists(indir):
        os.makedirs(indir)
    
    if bio:
        bparfullfile =  os.path.join(*(indir, '{}.bio.in'.format(outbase)))
        writeromsascii(bio, bparfullfile, filetype='bio')
        d['BPARNAM'] = bparfullfile
        
    if ice:  
        iparfullfile =  os.path.join(*(indir, '{}.ice.in'.format(outbase)))
        writeromsascii(ice, iparfullfile, filetype='ice')
        d['IPARNAM'] = iparfullfile
        
    if stations:
        sposfullfile =  os.path.join(*(indir, '{}.stations.in'.format(outbase)))
        writeromsascii(stations, sposfullfile, filetype='stations')
        d['SPOSNAM'] = sposfullfile
        
    # Check for restart files
    
    rstinfo = parserst(os.path.join(*(outdir, outbase)))
    if rstinfo['lastfile']:
        count = rstinfo['count']
        d['ININAME'] = rstinfo['lastfile']
        d['NRREC'] = -1
    
    # Run ROMS
    
    outbasetmp = '{}_{:02d}'.format(outbase, count)
    logbasetmp = '{}_{:02d}'.format(outbase, count)
    oceantmp = '{}_{:02d}.ocean.in'.format(outbase, count)
    
    rundata = createinputfiles(d, outbasetmp, logbasetmp, 
                               outdir=outdir, logdir=logdir, 
                               indir=indir, oceanfile=oceantmp)
    
    print('Running {}'.format(outbasetmp))
    if dryrun:
        runroms(rundata, **mpivars, dryrun=True)
        cleanexit = True
        return cleanexit
        
    runroms(rundata, **mpivars)

    # Parse the log file to make sure ROMS ran cleanly, and to check for a
    # blowup
    
    r = parseromslog(rundata['out'])
    
    cleanexit = True
    
    if not r['cleanrun']:
        print('ROMS crashed')
        cleanexit = False
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
        
        dateblewup = timevars['datestart'] + timedelta(seconds=r['laststep']*faststep.total_seconds())
        newend = max([dateblewup + timedelta(days=30), timevars['dateend']])
        
        filltimevars(d, slowstep, timevars['datestart'], newend,
                     timevars['tref'], dthis=timevars['dthis'], dtavg=timevars['dtavg'],
                     dtsta=timevars['dtsta'], dtrst=timevars['dtrst'],
                     dtdefhis=timevars['dtdefhis'], dtdefavg=timevars['dtdefavg'])
        
        # Run ROMS
        
        outbasetmp = '{}_{:02d}'.format(outbase, count)
        logbasetmp = '{}_{:02d}_slow'.format(outbase, count)
        oceantmp = 'ocean.tmp{:02d}.in'.format(count)
        
        rundata = createinputfiles(d, outbasetmp, logbasetmp, 
                                   outdir=outdir, logdir=logdir, 
                                   indir=indir, oceanfile=oceantmp)
    
        print('Running {} (slow)'.format(outbasetmp))        
        runroms(rundata, **mpivars)
        
        # Parse this run's log to make sure it finished cleanly
        
        r = parseromslog(rundata['out'])
        if not r['cleanrun']:
            print('ROMS crashed')
            cleanexit = False
            return cleanexit
            
        if r['blowup']:
            print('ROMS blew up with slow step; stopping simulation')
            return cleanexit
        
        # Switch back to fast time step and try to run until the end
        
        count+=1
        
        d['NRREC'] = -1
        d['ININAME'] = r['lasthis']
        
        filltimevars(d, faststep, timevars['datestart'], timevars['dateend'],
                     timevars['tref'], dthis=timevars['dthis'], dtavg=timevars['dtavg'],
                     dtsta=timevars['dtsta'], dtdefhis=timevars['dtdefhis'], 
                     dtdefavg=timevars['dtdefavg'], dtrst=timevars['dtrst'])
        
        outbasetmp = '{}_{:02d}'.format(outbase, count)
        logbasetmp = '{}_{:02d}'.format(outbase, count)
        oceantmp = 'ocean.tmp{:02d}.in'.format(count)
        
        rundata = createinputfiles(d, outbasetmp, logbasetmp, 
                                   outdir=outdir, logdir=logdir, 
                                   indir=indir, oceanfile=oceantmp)
    
        print('Running {}'.format(outbasetmp))        
        runroms(rundata, **mpivars)
    
        # Parse this run's log to check for another blowup
    
        r = parseromslog(rundata['out'])
        if not r['cleanrun']:
            print('ROMS crashed')
            cleanexit = False
            return cleanexit
    
    print('Done')
    return cleanexit


    















