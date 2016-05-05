# -*- coding: utf-8 -*-
"""
Created on Mon May  2 10:21:05 2016

@author: kelly.kearney
"""

import romsascii as r
import defaultparams
import bering10kinput as b10k
from datetime import datetime, timedelta

#--------------------
# Setup
#--------------------

# Default parameters: physical and NPZ 
    
ocean = defaultparams.ocean()  
npz = defaultparams.bestnpz()

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

logdir = '../Log'
outdir = '../Out/PythonTest/'

ocean['NRREC'] = -1
ocean['ININAME'] = '../Out/HindcastYKChinook/core_01_his_00011.nc'

ocean = r.filltimevars(ocean, **timevars)

print('Single run test')
r.runroms(ocean, 'pythontest', 'pythontest', mpivars, 
          outdir=outdir, logdir=logdir, bio=npz, dryrun=True)

# Test2: Same as above run, bun run through the blowup by slow-stepping for a 
# month when we hit it.

logdir = '../Log'
outdir = '../Out/PythonTestBlowup/'

print('Blowup run test')
r.runromsthroughblowup(ocean, 'pythontest_bu', timevars, mpivars, 
                       logdir=logdir, outdir=outdir, faststep=fast, 
                       slowstep=slow, bio=npz)
