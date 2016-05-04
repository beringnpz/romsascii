# -*- coding: utf-8 -*-
"""
Created on Mon May  2 10:21:05 2016

@author: kelly.kearney
"""

import romsutilities as r
import datetime

# Setup for a run    
    
test = r.defaultparams()    
test = r.ini_frc_bry(test, '/share/storage2/Bering10k_input', 'core')
test = r.other_in(test, '/share/storage2/Bering10k_input')

# A short run, across a blowup spot, enough to stop, slow step, and resume

timevars = {
'tstep'     : datetime.timedelta(seconds=600),
'datestart' : datetime.datetime(1969,1,16),
#'dateend'   : datetime.datetime(2004,1,16),'
'dateend'   : datetime.datetime(1971,3,1),
'tref'      : datetime.datetime(1900,1,1),
'dthis'     : datetime.timedelta(weeks=1),
'dtavg'     : datetime.timedelta(weeks=1),
'dtsta'     : datetime.timedelta(hours=1),
'dtdefhis'  : datetime.timedelta(weeks=10),
'dtdefavg'  : datetime.timedelta(weeks=10)}

fast = datetime.timedelta(seconds=600)
slow = datetime.timedelta(seconds=300)

logdir = '../Log'
outdir = '../Out/PythonTest/'

# MPI variables

mpivars = {
'mpiexe': '/usr/mpi/intel/openmpi-1.4.1/bin/mpirun',
'np':   140,
'hostfile': 'hostfile_0-7',
'romsexe': 'oceanM_npz_srb'}


# Test a single run, using an initialization time right before a blowup

test['NRREC'] = -1
test['ININAME'] = '../Out/HindcastYKChinook/core_01_his_00011.nc'

r.filltimevars(test, timevars['tstep'], timevars['datestart'], timevars['dateend'], 
                 timevars['tref'], timevars['dthis'], timevars['dtavg'], 
                 timevars['dtsta'], timevars['dtdefhis'], timevars['dtdefavg'])

# print('Single run test')
# r.runroms(test, 'ocean_NBSnpz.template.in', outdir, logdir, 'pythontest', 'pythontest', mpivars, dryrun=False)

# Test running through blowup

logdir = '../Log'
outdir = '../Out/PythonTestBlowup/'

print('Blowup run test')
r.runromsthroughblowup(test, 'ocean_NBSnpz.template.in', logdir, outdir, 'pythontest_bu', fast, slow, timevars, mpivars)
