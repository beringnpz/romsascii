# -*- coding: utf-8 -*-
"""
Created on Mon May  2 10:21:05 2016

@author: kelly.kearney
"""

from romsascii import romsascii as r
from romsascii import defaultparams
from romsascii import bering10kinput as b10k
from datetime import datetime, timedelta

#--------------------
# Setup
#--------------------

# Default parameters: physical, NPZ, ice, stations 
    
ocean = defaultparams.ocean()  
npz = defaultparams.bestnpz()
ice = defaultparams.ice()
stat = defaultparams.stations()

# Fill in Bering 10K input files, using CORE forcing
  
b10k.ini_frc_bry(ocean, '/share/storage2/Bering10k_input', 'core')
b10k.other_in(ocean, '/share/storage2/Bering10k_input')

# Fill in time variables
# Run a short period that crosses a blowup spot, to test the 
#   runromsacrossblowup routine
# Weekly archiving, 10 weeks per file
# Time step = 600 s, halved to get across blowups

timevars = {
'tstep'     : timedelta(seconds=600),
'datestart' : datetime(1969,1,16),
'dateend'   : datetime(1971,3,1),
'tref'      : datetime(1900,1,1),
'dthis'     : timedelta(weeks=1),
'dtavg'     : timedelta(weeks=1),
'dtsta'     : timedelta(hours=1),
'dtdefhis'  : timedelta(weeks=10),
'dtdefavg'  : timedelta(weeks=10)}

fast = timedelta(seconds=600)
slow = timedelta(seconds=300)

# MPI variables

mpivars = {
'mpiexe': '/usr/mpi/intel/openmpi-1.4.1/bin/mpirun',
'np':   140,
'hostfile': 'hostfile_0-7',
'romsexe': 'oceanM_npz_srb'}

#--------------------
# Test sims
#--------------------

# Test 1: A simple call to the ROMS executable.  This run should stop when it 
# blows up.  Start from a previous history file

logdir = 'testlog/'
outdir = 'testout/'
indir  = 'testin/'

ocean['NRREC'] = -1
ocean['ININAME'] = '../Out/HindcastYKChinook/core_01_his_00011.nc'

ocean = r.filltimevars(ocean, **timevars)

print('Single run test')
r.runroms(ocean, 'pythontest', 'pythontest', mpivars, 
          outdir=outdir, logdir=logdir, indir=indir, 
          bio=npz, ice=ice, stations=stat, dryrun=False)

# Test 2: Same as above run, but run through the blowup by slow-stepping for a 
# month when we hit it.

indir = 'testin2/'

print('Blowup run test')
r.runromsthroughblowup(ocean, 'pythontest_bu', timevars, mpivars, 
                       logdir=logdir, outdir=outdir, indir=indir, faststep=fast, 
                       slowstep=slow, bio=npz, ice=ice, stations=stat)
