# -*- coding: utf-8 -*-
"""
Module to store file locations for Bering 10K ROMS App

This module holds the locations of input files used in the Bering10K ROMS 
simulations, specifically, the files used for the Yukon-Kuskokwim Chinook 
project simulations.

It also provides some functions to create initialization files for the 
BEST_NPZ and FEAST biological modules, and to parse history and restart 
files under my typical simulation folder tree.

This is extremely application- and computer-specific... not intended for 
distribution.

Created on Wed May  4 15:21:46 2016

@author: kakearney
"""

import os

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
    
    cfsr1 = ('cfs-Pair-2002-2014.nc',
             'cfs-Qair-2002-2014.nc',
             'cfs-Tair-2002-2014.nc',
             'cfs-lwrad-reduced-3-percent-2002-2014.nc',
             'cfs-swradave-reduced-10-percent-2002-2014.nc',
             'cfs-Uwind-2002-2014.nc',
             'cfs-Vwind-2002-2014.nc',
             'cfs-rain-2002-2014.nc'
            )
            
    cfsr2 = ('cfs-Pair-2015-2016.nc',
             'cfs-Qair-2015-2016.nc',
             'cfs-Tair-2015-2016.nc',
             'cfs-lwrad-reduced-by-3-percent-2015-2016.nc',
             'cfs-swradave-reduced-by-10-percent-2015-2016.nc',
             'cfs-Uwind-2015-2016.nc',
             'cfs-Vwind-2015-2016.nc',
             'cfs-rain-2015-2016.nc'
            )
    
    if group == 'core':    # 1969/01/07 - 2002/01/02 
        frc = core + same
        bry = 'core_bry.1969-2004.nc'
        ini = 'ini.1969.nc'
    elif group == 'cfsr1': # 2002/01/02 - 2013-01-04
        frc = cfsr1 + same
        bry = 'cfsr_bc_2002-2014_5day_bio_10layer.nc'
        ini = 'ini.2004.nc'
    elif group == 'cfsr2': # 2013/01/04 - 2014/12/31
        frc = cfsr1 + same
        bry = 'bering10k_bc_2013-2016_with_bio.nc'
        ini = 'ini.2004.nc'
    elif group == 'cfsr3': # 2014/12/31 - 2016/08/25
        frc = cfsr2 + same
        bry = 'bering10k_bc_2013-2016_with_bio.nc'
        ini = 'ini.2004.nc'
    
    
    frc = list(map(lambda x: os.path.join(*(indir,'frc',x)), frc))
    d['FRCNAME'] = frc
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
    
def parsehis(simdir, shortname):
    """
    Parse restart counter from ROMS simulation history files
    
    This function finds the name of, and parses the simulation counter, 
    from a series of ROMS history files.  It assumes that those files 
    were using my standard naming scheme, where each file is named
        simdir/Out/shortname_XX_his_YYYY.nc
    where XX is the counter for number of restarts (assigned by 
    runromsthroughblowup or similar), and YYYY is the counter for number 
    of history files (assigned by ROMS.)
    
    Args:
        simdir:     Name of simulation folder
        shortname:  base name for history files
    
    Returns:
        d:          dictionary object with the following keys:
                    lastfile:   full path to last history file
                    cnt:        restart counter of last file incremented 
                                by 1 (i.e. count you would want to 
                                restart with in runromsthroughblowup)
      
    """
    allhis = sorted(glob.glob(os.path.join(simdir, "Out", shortname + "*his*.nc")))
    his = allhis[-1]

    pattern = shortname + "_(\d+)_his_(\d+)"
    m = re.search(pattern, his)
    cnt = int(m.group(1)) + 1
    
    return {'lastfile': his, 'count': cnt}
    
def parserst(simdir, shortname):
    """
    Parse restart counters from ROMS simulation restart files
    
    This function finds the name of, and parses the simulation counter, 
    from a series of ROMS restart files.  It assumes that those files 
    were using my standard naming scheme, where each file is named
        simdir/Out/shortname_XX_rst.nc
    where XX is the counter for number of restarts (assigned by 
    runromsthroughblowup or similar.)
    
    Args:
        simdir:     Name of simulation folder
        shortname:  base name for history files
    
    Returns:
        d:          dictionary object with the following keys:
                    lastfile:   full path to last restart file
                    cnt:        restart counter of last file incremented 
                                by 1 (i.e. count you would want to 
                                restart with in runromsthroughblowup)
    """
    allrst = sorted(glob.glob(os.path.join(simdir, "Out", shortname + "*rst*.nc")))
    if len(allrst) == 0:
        rst = []
        cnt = 1
    else:
        rst = allrst[-1]

        pattern = shortname + "_(\d+)_rst.nc"
        m = re.search(pattern, rst)
        cnt = int(m.group(1)) + 1
    
    return {'lastfile': rst, 'count': cnt}
    
def buildinifile(iniphys, ininpz):
    """
    Create ROMS initialization file for a BEST_NPZ simulation
    
    This function adds biological variables to a physics-only history 
    file.  It is intended to make it easy to start a new BEST_NPZ 
    simulation from any time point in a physics run that used the same 
    ROMS domain.  The benthic variables (Ben and DetBen) are initialized
    to specific values; all others are initialized to 0 (with the 
    assumption that analytical initial conditions will be used.)
    
    Args:
        iniphys:    full file name of physics-only history file
        ininpz:     full file name of new BEST_NPZ initialization file
    """
    shutil.copyfile(iniphys, ininpz)

    f = nc.Dataset(ininpz, 'r+')
    sz = (f.dimensions['ocean_time'].size, f.dimensions['s_rho'].size,
          f.dimensions['eta_rho'].size,    f.dimensions['xi_rho'].size)

    vars2D = [['IcePhL', "Ice algae concentration",                "ice algae conc",         "mgC/m2"],
              ['IceNO3', "Ice nitrate concentration",              "ice nitrate conc",       "mgC/m2"],
              ['IceNH4', "Ice Ammonium concentration",             "ice Ammonium conc",      "mgC/m2"],
              ['IceLog', "Logical Ice Counter",                    "Logical Ice Counter",    "unitless"],
              ['DetBen', "benthic detritus concentration",         "benthic detritus",       "milligram carbon meter-2"],
              ['Ben',    "benthos concentration",                  "benthos",                "milligram carbon meter-2"]
             ]

    vars3D = [['NO3',    "nitrate concentration",                  "nitrate",                "millimole nitrogen meter-3"],
              ['NH4',    "ammonia concentration",                  "ammonia",                "millimole nitrogen meter-3"],
              ['PhS',    "small phytoplankton concentration",      "small phytoplankton",    "milligram carbon meter-3"],
              ['PhL',    "large phytoplankton concentration",      "large phytoplankton",    "milligram carbon meter-3"],
              ['MZS',    "small microzooplankton concentration",   "small microzooplankton", "milligram carbon meter-3"],
              ['MZL',    "large microzooplankton concentration",   "large microzooplankton", "milligram carbon meter-3"],
              ['Cop',    "small coastal copepod concentration",    "copepod",                "milligram carbon meter-3"],
              ['NCaS',   "neocalanus spp. concentration",          "neocalanus",             "milligram carbon meter-3"],
              ['NCaO',   "Offshore neocalanus spp. concentration", "neocalanus",             "milligram carbon meter-3"],
              ['EupS',   "euphausiid concentration",               "euphausiid",             "milligram carbon meter-3"],
              ['EupO',   "Offshore euphausiid concentration",      "euphausiid",             "milligram carbon meter-3"],
              ['Det',    "detritus concentration",                 "detritus",               "milligram carbon meter-3"],
              ['DetF',   "Fast sinking detritus concentration",    "detritus",               "milligram carbon meter-3"],
              ['Jel',    "Jellyfish concentration",                "jellyfish",              "milligram carbon meter-3"],
              ['Iron',   "iron concentration",                     "iron",                   "micromol Fe m-3"]
             ]

    for var in vars2D:
        dvar = f.createVariable(var[0], 'f4', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=1e37)
        dvar.long_name = var[1]
        dvar.units = var[3]
        dvar.time = "ocean_time"
        dvar.coordinates = "lon_rho lat_rho ocean_time"
        dvar.field = '{}, scalar, series'.format(var[2])
        if var[0] == 'Ben':
            dvar[:] = np.ones([sz[i] for i in [0,2,3]])*8000.0;
        elif var[0] == 'DetBen':
            dvar[:] = np.ones([sz[i] for i in [0,2,3]])*500.0;
        else:
            dvar[:] = np.zeros([sz[i] for i in [0,2,3]])

    for var in vars3D:
        dvar = f.createVariable(var[0], 'f4', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=1e37)
        dvar.long_name = var[1]
        dvar.units = var[3]
        dvar.time = "ocean_time"
        dvar.coordinates = "lon_rho lat_rho s_rho ocean_time"
        dvar.field = '{}, scalar, series'.format(var[2])
        dvar[:] = np.zeros(sz)

    f.close()
    


