# -*- coding: utf-8 -*-
"""
Created on Mon May  2 10:21:05 2016

@author: kelly.kearney
"""

import romsascii as r
import defaultparams
import bering10kinput as b10k
import datetime

#--------------------
# Setup
#--------------------

# Default physical parameters   
    
ocean = defaultparams.ocean()  

# Fill in Bering 10K input files
  
b10k.ini_frc_bry(ocean, '/share/storage2/Bering10k_input', 'core')
b10k.other_in(ocean, '/share/storage2/Bering10k_input')

# Fill in time variables
# Run a short period that crosses a blowup spot, to test the 
#   runromsacrossblowup routine
# Weekly archiving, 10 weeks per file
# Time step = 600 s, halved to get across blowups

timevars = {
'tstep'     : datetime.timedelta(seconds=600),
'datestart' : datetime.datetime(1969,1,16),
'dateend'   : datetime.datetime(1971,3,1),
'tref'      : datetime.datetime(1900,1,1),
'dthis'     : datetime.timedelta(weeks=1),
'dtavg'     : datetime.timedelta(weeks=1),
'dtsta'     : datetime.timedelta(hours=1),
'dtdefhis'  : datetime.timedelta(weeks=10),
'dtdefavg'  : datetime.timedelta(weeks=10)}

fast = datetime.timedelta(seconds=600)
slow = datetime.timedelta(seconds=300)

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

test['NRREC'] = -1
test['ININAME'] = '../Out/HindcastYKChinook/core_01_his_00011.nc'

r.filltimevars(test, timevars['tstep'], timevars['datestart'], timevars['dateend'], 
                 timevars['tref'], timevars['dthis'], timevars['dtavg'], 
                 timevars['dtsta'], timevars['dtdefhis'], timevars['dtdefavg'])

 print('Single run test')
# r.runroms(ocean, 'ocean_NBSnpz.template.in', outdir, logdir, 'pythontest', 'pythontest', mpivars, dryrun=True)

# Test running through blowup

logdir = '../Log'
outdir = '../Out/PythonTestBlowup/'

print('Blowup run test')
#r.runromsthroughblowup(test, 'ocean_NBSnpz.template.in', logdir, outdir, 'pythontest_bu', fast, slow, timevars, mpivars)
