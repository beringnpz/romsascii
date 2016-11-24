# -*- coding: utf-8 -*-
"""
Module to store file locations for Bering 10K ROMS App

This module holds the locations of input files used in the Bering10K ROMS 
simulations, specifically, the files used for the Yukon-Kuskokwim Chinook 
project simulations.

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
             'cfs-lwrad-reduced-3-percent-2015-2016.nc',
             'cfs-swradave-reduced-10-percent-2015-2016.nc',
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

