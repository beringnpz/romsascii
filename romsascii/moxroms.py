from romsascii import romsascii as r
import os
import sys
from datetime import datetime, timedelta
import netCDF4 as nc
import csv
import glob

"""
Bering10K-on-mox module

This module provides wrappers around some common tasks used when running the Bering10K 
ROMS model on the mox-hyak (UW) computer.  It's very application-specific and not intended
for general use.
"""

def setinfiles(d, ininame, ncinputfolder,
               nlayer=10,
               breakyr=1995,
               setvinfo=True,
               tide='tides_OTBS.nc',
               river='runoff.kearney.efol20.updated201809.nc',
               sss='sss.clim.nc',
               coreadjust=True):
    """
    Modify input file name variables in parameter dict based on year
    This function is applicable to all hindcast runs
    """
    # Grid file
    d['GRDNAME'] = os.path.join(ncinputfolder,'grd','Bering_grid_withFeast.nc')

    # BESTNPZ/FEAST varinfo file
    
    if setvinfo:
        d['VARNAME'] = os.path.join(ncinputfolder,'var','varinfo_bestnpzfeast_new.dat')
        
    # Forcing files used for the entire timeperiod: SSS, tides,
    # runoff
    frcclim = (tide,
               sss
                )
    frcclim = list(map(lambda x: os.path.join(ncinputfolder,'generic',x), frcclim))

    # Initialization file
    d['ININAME'] = ininame
    
    # In what year does this initialization file start?
    
    f = nc.Dataset(ininame, 'r')
    tunit = f.variables['ocean_time'].units
    tcal = f.variables['ocean_time'].calendar
    tini = max(nc.num2date(f.variables['ocean_time'][:], units=tunit, calendar=tcal))
    
    yr = tini.year
    if (datetime(yr+1,1,1) - tini) < timedelta(hours=6):
        yr = yr+1 # in overlap period... start next year

    # Year-specific (and layer-specific) boundary file
    
    if (nlayer == 10):
        bryfol = 'bry'
    else:
        bryfol = 'bry{:d}'.format(nlayer)
        
    if yr < breakyr:
        d['BRYNAME'] = os.path.join(ncinputfolder,'hindcast_core', '{}'.format(yr), 'roms-core-bry-N{}-{}.nc'.format(nlayer,yr))
        # d['BRYNAME'] = os.path.join(ncinputfolder,bryfol,'roms-core-bry-{}.nc'.format(yr))
    else:
        d['BRYNAME'] = os.path.join(ncinputfolder,'hindcast_cfs', '{}'.format(yr), 'roms-cfs-bry-N{}-{}.nc'.format(nlayer,yr))
        # d['BRYNAME'] = os.path.join(ncinputfolder,bryfol,'roms-cfs-bry-{}.nc'.format(yr))

    # Add year-specific forcing files to full-sim ones
    
    if yr < breakyr:
        if coreadjust:
            frc = [os.path.join(ncinputfolder,'hindcast_core', '{}'.format(yr), 'roms-core-atmos-lwrad-increased-{}.nc'.format(yr)),
                   os.path.join(ncinputfolder,'hindcast_core', '{}'.format(yr), 'roms-core-atmos-Pair-{}.nc'.format(yr)),
                   os.path.join(ncinputfolder,'hindcast_core', '{}'.format(yr), 'roms-core-atmos-Qair-{}.nc'.format(yr)),
                   os.path.join(ncinputfolder,'hindcast_core', '{}'.format(yr), 'roms-core-atmos-rain-{}.nc'.format(yr)),
                   os.path.join(ncinputfolder,'hindcast_core', '{}'.format(yr), 'roms-core-atmos-swrad-increased-{}.nc'.format(yr)),
                   os.path.join(ncinputfolder,'hindcast_core', '{}'.format(yr), 'roms-core-atmos-Tair-{}.nc'.format(yr)),
                   os.path.join(ncinputfolder,'hindcast_core', '{}'.format(yr), 'roms-core-atmos-Uwind-{}.nc'.format(yr)),
                   os.path.join(ncinputfolder,'hindcast_core', '{}'.format(yr), 'roms-core-atmos-Vwind-{}.nc'.format(yr))]
        else:
           frc = [os.path.join(ncinputfolder,'hindcast_core', '{}'.format(yr), 'roms-core-atmos-lwrad-{}.nc'.format(yr)),
                  os.path.join(ncinputfolder,'hindcast_core', '{}'.format(yr), 'roms-core-atmos-Pair-{}.nc'.format(yr)),
                  os.path.join(ncinputfolder,'hindcast_core', '{}'.format(yr), 'roms-core-atmos-Qair-{}.nc'.format(yr)),
                  os.path.join(ncinputfolder,'hindcast_core', '{}'.format(yr), 'roms-core-atmos-rain-{}.nc'.format(yr)),
                  os.path.join(ncinputfolder,'hindcast_core', '{}'.format(yr), 'roms-core-atmos-swrad-{}.nc'.format(yr)),
                  os.path.join(ncinputfolder,'hindcast_core', '{}'.format(yr), 'roms-core-atmos-Tair-{}.nc'.format(yr)),
                  os.path.join(ncinputfolder,'hindcast_core', '{}'.format(yr), 'roms-core-atmos-Uwind-{}.nc'.format(yr)),
                  os.path.join(ncinputfolder,'hindcast_core', '{}'.format(yr), 'roms-core-atmos-Vwind-{}.nc'.format(yr))]
    else: 
        frc = [os.path.join(ncinputfolder,'hindcast_cfs', '{}'.format(yr), 'roms-cfs-atmos-lwrad-{}.nc'.format(yr)),
               os.path.join(ncinputfolder,'hindcast_cfs', '{}'.format(yr), 'roms-cfs-atmos-Pair-{}.nc'.format(yr)),
               os.path.join(ncinputfolder,'hindcast_cfs', '{}'.format(yr), 'roms-cfs-atmos-Qair-{}.nc'.format(yr)),
               os.path.join(ncinputfolder,'hindcast_cfs', '{}'.format(yr), 'roms-cfs-atmos-rain-{}.nc'.format(yr)),
               os.path.join(ncinputfolder,'hindcast_cfs', '{}'.format(yr), 'roms-cfs-atmos-swrad-{}.nc'.format(yr)),
               os.path.join(ncinputfolder,'hindcast_cfs', '{}'.format(yr), 'roms-cfs-atmos-Tair-{}.nc'.format(yr)),
               os.path.join(ncinputfolder,'hindcast_cfs', '{}'.format(yr), 'roms-cfs-atmos-Uwind-{}.nc'.format(yr)),
               os.path.join(ncinputfolder,'hindcast_cfs', '{}'.format(yr), 'roms-cfs-atmos-Vwind-{}.nc'.format(yr))]

    d['FRCNAME'] = frcclim + frc + [os.path.join(ncinputfolder, 'river', 'runoff_kearney_{}.nc'.format(yr))]
    d['NFFILES'] = len(d['FRCNAME'])
    
    return tini
    
def runhindcast(ocean, simdir, simname, inifile, enddate, mpivars, timevars, fast, slow,
                nrrec=0, 
                ncinputfolder='/gscratch/bumblereem/bering10k/input/',
                dryrunflag=False,
                bio={}, ice={}, stations={}, 
                setvinfo=True,
                river='runoff.kearney.efol20.updated201809.nc',
                breakyr=1995,
                coreadjust=True):
    """
    Run the hindcast.  This takes care of all the messiness of switching up input files, 
    checking for blowups and reducing timesteps, and resuming partially-completed runs.
                
    Returns:
                status: 'success', 'error', 'blowup', or 'dryrun'
    """
                
    # Set up input, output, and log folders

    outdir = os.path.join(simdir, "Out")
    indir  = os.path.join(simdir, "In")
    logdir = os.path.join(simdir, "Log")
                
    # Initialization file: check for any existing restart files, and if not 
    # found, start from the initialization file

    rstinfo = r.parserst(os.path.join(outdir, simname))
    if rstinfo['lastfile']:
        cnt = rstinfo['count']
        ocean['ININAME'] = rstinfo['lastfile']
        ocean['NRREC'] = -1
    else:
        cnt = 1
        ocean['ININAME'] = inifile
        ocean['NRREC'] = nrrec
        
    # Set input files based on the initialization time

    tini = setinfiles(ocean, ocean['ININAME'], ncinputfolder, nlayer=ocean['N'], setvinfo=setvinfo, river=river, breakyr=breakyr, coreadjust=coreadjust)
    
    # Create ascii input files that will be reused across all restarts

    if not os.path.exists(indir):
        os.makedirs(indir, 0o775)
        os.chmod(indir, 0o775)
    if not os.path.exists(outdir):
        os.makedirs(outdir, 0o775)
        os.chmod(outdir, 0o775)
    if not os.path.exists(logdir):
        os.makedirs(logdir, 0o775)
        os.chmod(logdir, 0o775)

    if bio:
        bparfullfile =  os.path.join(*(indir, '{}.bio.in'.format(simname)))
        r.writeromsascii(bio, bparfullfile, filetype='bio', consecstep=0) # need to compress the Hout flags to be readable past 80-ish characters
        ocean['BPARNAM'] = bparfullfile
        os.chmod(bparfullfile, 0o644)

    if ice:
        iparfullfile =  os.path.join(*(indir, '{}.ice.in'.format(simname)))
        r.writeromsascii(ice, iparfullfile, filetype='ice')
        ocean['IPARNAM'] = iparfullfile
        os.chmod(iparfullfile, 0o644)
        
    if stations:
        sposfullfile =  os.path.join(*(indir, '{}.stations.in'.format(simname)))
        r.writeromsascii(stations, sposfullfile, filetype='stations')
        ocean['SPOSNAM'] = sposfullfile
        os.chmod(sposfullfile, 0o644)

    # Create log file to document slow-stepping time periods

    steplog = os.path.join(logdir, 'step_{}.txt'.format(simname))

    if not os.path.isfile(steplog):
        fstep = open(steplog, "w+")
        fstep.close()

    # Run sim

    while tini < (enddate - timedelta(hours=6)):
    
        # Set end date as furthest point we can run.  This will be either 
        # the simulation end date, the end of the current set of 
        # forcing/boundary files (i.e. end of the year), or the end of the 
        # slow-stepping period (if we are in one), whichever comes first
    
        yr = tini.year
        if (datetime(yr+1,1,1) - tini) < timedelta(hours=6):
            yr = yr+1 # in overlap period...
        
        endyear = datetime(yr+1,1,1)
    
        timevars['tstep'] = fast
        endslow = endyear
        with open(steplog) as fstep:
            readCSV = csv.reader(fstep, delimiter=',')
            for row in readCSV:
                t1 = datetime.strptime(row[0], '%Y-%m-%d-%H-%M:%S')
                t2 = datetime.strptime(row[1], '%Y-%m-%d-%H-%M:%S')
                if (tini >= t1) & (tini <= (t2-timedelta(hours=4))): # in a slow-step period
                    endslow = t2
                    timevars['tstep'] = slow
        
        timevars['dateend'] = min(enddate, endyear, endslow)
        r.filltimevars(ocean, **timevars)
    
        # Create ascii input files
    
        oceanfile = '{}_{:02d}.ocean.in'.format(simname, cnt)
        logbase = simname
    
        rundata = r.createinputfiles(ocean, simname, logbase, 
                                   outdir=outdir, logdir=logdir, 
                                   indir=indir, oceanfile=oceanfile, 
                                   addrstcount=True, addstacount=True, 
                                   addlogcount=True, addhiscount=True, 
                                   addavgcount=True, count=cnt)
                               
        # Run simulation bit
    
        print('Running ROMS hindcast/nowcast')
        print('  Counter block: {}'.format(cnt))
        print('  Start date: {}'.format(tini.strftime('%Y-%m-%d %H:%M:%S')))
        print('  End date:   {}'.format(timevars['dateend'].strftime('%Y-%m-%d %H:%M:%S')))
        print('  Time step:  {}'.format(str(timevars['tstep'])))
        print('  Executable: {}'.format(mpivars['romsexe']))
        print('  Input file: {}'.format(rundata['in']))
        print('  Standard output: {}'.format(rundata['out']))
        print('  Standard error:  {}'.format(rundata['err']))
        r.runroms(rundata, **mpivars, dryrun=dryrunflag, addnp=False)
    
        if dryrunflag:
            return 'dryrun'
        
        rsim = r.parseromslog(rundata['out'])
    
        # Did the run crash (i.e. anything but successful end or blowup)? If 
        # so, we'll exit now
    
        if (not rsim['cleanrun']) & (not rsim['blowup']): 
            print('  Similation block terminated with error')
            return 'error'
        
        # Did it blow up?  If it did so during a slow-step period, we'll exit
        # now.  If it blew up during a fast-step period, set up a new 
        # slow-step period and reset input to start with last history file.  
        # If it ran to completion, reset input to start with last restart 
        # file
    
        rstinfo = r.parserst(os.path.join(outdir, simname))
        cnt = rstinfo['count']
    
        if rsim['blowup']:
            if timevars['tstep'] == slow: 
                print('  Simulation block blew up in a slow-step period')
                return 'blowup'
        
            # Find the most recent history file written to
            hisfile = rsim['lasthis']
        
            if not hisfile: # non-clean blowup, no his file defined
                allhis = sorted(glob.glob(os.path.join(outdir, simname + "*his*.nc")))
                hisfile = allhis[-1]
                    
            fhis = nc.Dataset(hisfile)
            if len(fhis.variables['ocean_time']) == 0:
                allhis = glob.glob(os.path.join(outdir, simname + "*his*.nc"))
                allhis = sorted(list(set(allhis) - set([hisfile])))
                hisfile = allhis[-1]
        
            tini = setinfiles(ocean, hisfile, ncinputfolder, nlayer=ocean['N'], setvinfo=setvinfo, river=river, breakyr=breakyr, coreadjust=coreadjust)
            ocean['ININAME'] = hisfile
            ocean['NRREC'] = -1
        
            # dateblewup = timevars['datestart'] + timedelta(seconds=rsim['laststep']*timevars['tstep'].total_seconds())
        
            t1 = tini.strftime('%Y-%m-%d-%H-%M:%S')
            t2 = (tini + timedelta(days=30)).strftime('%Y-%m-%d-%H-%M:%S')
            fstep = open(steplog, "a+")
            fstep.write('{},{}\n'.format(t1,t2))
            fstep.close()
        
        else:
            tini = setinfiles(ocean, rstinfo['lastfile'], ncinputfolder, nlayer=ocean['N'], setvinfo=setvinfo, river=river, breakyr=breakyr, coreadjust=coreadjust)
            ocean['ININAME'] = rstinfo['lastfile']
            ocean['NRREC'] = -1
        
    # Print completion status message
        
    print('Simulation completed through specified end date')
    return 'success'
    
def runforecast(ocean, simdir, simname, inifile, enddate, mpivars, timevars, fast, slow,
                nrrec=0, 
                dryrunflag=False,
                bio={}, ice={}, stations={}):
    """
    Run a forecast.  Same as runhindcast but without the need to switch input files as it 
    goes (assumes all input files have already been set up in ocean dictionary)
                
    Returns:
                status: 'success', 'error', 'blowup', or 'dryrun'
    """
                
    # Set up input, output, and log folders

    outdir = os.path.join(simdir, "Out")
    indir  = os.path.join(simdir, "In")
    logdir = os.path.join(simdir, "Log")
                
    # Initialization file: check for any existing restart files, and if not 
    # found, start from the initialization file

    rstinfo = r.parserst(os.path.join(outdir, simname))
    if rstinfo['lastfile']:
        cnt = rstinfo['count']
        ocean['ININAME'] = rstinfo['lastfile']
        ocean['NRREC'] = -1
    else:
        cnt = 1
        ocean['ININAME'] = inifile
        ocean['NRREC'] = nrrec
        
    # Check that all input files exist (better to do this here than let ROMS try and fail)
    
    for fl in ocean['FRCNAME']+[ocean['GRDNAME']]+[ocean['VARNAME']]+[ocean['BRYNAME']]+[ocean['ININAME']]:
        if not os.path.isfile(fl):
            print('WARNING!: Cannot find file {}'.format(fl))
            sys.exit()    
        
    # Get starting time from initialization file

    f = nc.Dataset(ocean['ININAME'], 'r')
    tunit = f.variables['ocean_time'].units
    tcal = f.variables['ocean_time'].calendar
    tini = max(nc.num2date(f.variables['ocean_time'][:], units=tunit, calendar=tcal))
    
    # Create ascii input files that will be reused across all restarts

    if not os.path.exists(indir):
        os.makedirs(indir, 0o775)
        os.chmod(indir, 0o775)
    if not os.path.exists(outdir):
        os.makedirs(outdir, 0o775)
        os.chmod(outdir, 0o775)
    if not os.path.exists(logdir):
        os.makedirs(logdir, 0o775)
        os.chmod(logdir, 0o775)

    if bio:
        bparfullfile =  os.path.join(*(indir, '{}.bio.in'.format(simname)))
        r.writeromsascii(bio, bparfullfile, filetype='bio', consecstep=0) # need to compress the Hout flags to be readable past 80-ish characters
        ocean['BPARNAM'] = bparfullfile
        os.chmod(bparfullfile, 0o644)

    if ice:
        iparfullfile =  os.path.join(*(indir, '{}.ice.in'.format(simname)))
        r.writeromsascii(ice, iparfullfile, filetype='ice')
        ocean['IPARNAM'] = iparfullfile
        os.chmod(iparfullfile, 0o644)
        
    if stations:
        sposfullfile =  os.path.join(*(indir, '{}.stations.in'.format(simname)))
        r.writeromsascii(stations, sposfullfile, filetype='stations')
        ocean['SPOSNAM'] = sposfullfile
        os.chmod(sposfullfile, 0o644)

    # Create log file to document slow-stepping time periods

    steplog = os.path.join(logdir, 'step_{}.txt'.format(simname))

    if not os.path.isfile(steplog):
        fstep = open(steplog, "w+")
        fstep.close()

    # Run sim

    while tini < (enddate - timedelta(hours=6)):
    
        # Set end date as furthest point we can run.  This will be either 
        # the simulation end date, the end of the current set of 
        # forcing/boundary files (i.e. end of the year), or the end of the 
        # slow-stepping period (if we are in one), whichever comes first
    
        timevars['tstep'] = fast
        endslow = enddate
        with open(steplog) as fstep:
            readCSV = csv.reader(fstep, delimiter=',')
            for row in readCSV:
                t1 = datetime.strptime(row[0], '%Y-%m-%d-%H-%M:%S')
                t2 = datetime.strptime(row[1], '%Y-%m-%d-%H-%M:%S')
                if (tini >= t1) & (tini <= (t2-timedelta(hours=4))): # in a slow-step period
                    endslow = t2
                    timevars['tstep'] = slow
        
        timevars['dateend'] = min(enddate, endslow)
        r.filltimevars(ocean, **timevars)
    
        # Create ascii input files
    
        oceanfile = '{}_{:02d}.ocean.in'.format(simname, cnt)
        logbase = simname
    
        rundata = r.createinputfiles(ocean, simname, logbase, 
                                   outdir=outdir, logdir=logdir, 
                                   indir=indir, oceanfile=oceanfile, 
                                   addrstcount=True, addstacount=True, 
                                   addlogcount=True, addhiscount=True, 
                                   addavgcount=True, count=cnt)
                               
        # Run simulation bit
    
        print('Running ROMS forecast')
        print('  Counter block: {}'.format(cnt))
        print('  Start date: {}'.format(tini.strftime('%Y-%m-%d %H:%M:%S')))
        print('  End date:   {}'.format(timevars['dateend'].strftime('%Y-%m-%d %H:%M:%S')))
        print('  Time step:  {}'.format(str(timevars['tstep'])))
        print('  Executable: {}'.format(mpivars['romsexe']))
        print('  Input file: {}'.format(rundata['in']))
        print('  Standard output: {}'.format(rundata['out']))
        print('  Standard error:  {}'.format(rundata['err']))
        r.runroms(rundata, **mpivars, dryrun=dryrunflag, addnp=False)
    
        if dryrunflag:
            return 'dryrun'
        
        rsim = r.parseromslog(rundata['out'])
    
        # Did the run crash (i.e. anything but successful end or blowup)? If 
        # so, we'll exit now
    
        if (not rsim['cleanrun']) & (not rsim['blowup']): 
            print('  Similation block terminated with error')
            return 'error'
        
        # Did it blow up?  If it did so during a slow-step period, we'll exit
        # now.  If it blew up during a fast-step period, set up a new 
        # slow-step period and reset input to start with last history file.  
        # If it ran to completion, reset input to start with last restart 
        # file
    
        rstinfo = r.parserst(os.path.join(outdir, simname))
        cnt = rstinfo['count']
    
        if rsim['blowup']:
            if timevars['tstep'] == slow: 
                print('  Simulation block blew up in a slow-step period')
                return 'blowup'
        
            # Find the most recent history file written to
            hisfile = rsim['lasthis']
        
            if not hisfile: # non-clean blowup, no his file defined
                allhis = sorted(glob.glob(os.path.join(outdir, simname + "*his*.nc")))
                hisfile = allhis[-1]
                    
            fhis = nc.Dataset(hisfile)
            if len(fhis.variables['ocean_time']) == 0:
                allhis = glob.glob(os.path.join(outdir, simname + "*his*.nc"))
                allhis = sorted(list(set(allhis) - set([hisfile])))
                hisfile = allhis[-1]
        
            # tini = setinfiles(ocean, hisfile, ncinputfolder, nlayer=ocean['N'])
            ocean['ININAME'] = hisfile
            ocean['NRREC'] = -1
        
            # dateblewup = timevars['datestart'] + timedelta(seconds=rsim['laststep']*timevars['tstep'].total_seconds())
            f = nc.Dataset(ocean['ININAME'], 'r')
            tunit = f.variables['ocean_time'].units
            tcal = f.variables['ocean_time'].calendar
            tini = max(nc.num2date(f.variables['ocean_time'][:], units=tunit, calendar=tcal))
        
            t1 = tini.strftime('%Y-%m-%d-%H-%M:%S')
            t2 = (tini + timedelta(days=30)).strftime('%Y-%m-%d-%H-%M:%S')
            fstep = open(steplog, "a+")
            fstep.write('{},{}\n'.format(t1,t2))
            fstep.close()
        
        else:
            ocean['ININAME'] = rstinfo['lastfile']
            ocean['NRREC'] = -1
            
            f = nc.Dataset(ocean['ININAME'], 'r')
            tunit = f.variables['ocean_time'].units
            tcal = f.variables['ocean_time'].calendar
            tini = max(nc.num2date(f.variables['ocean_time'][:], units=tunit, calendar=tcal))
        
    # Print completion status message
        
    print('Simulation completed through specified end date')
    return 'success'
