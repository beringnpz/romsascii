"""
FEAST helper functions module

This module provides several function to streamline the setup and running
of ROMS simulations with the FEAST biological module.  

The FEAST module itself is still in development, and therefore some of 
these functions (most notably, the buildfeastini function) are still very 
specific to our development environment and file structure.  We hope to
make these utilities more generic as a more stable version of FEAST 
emerges.
"""

import netCDF4 as nc
import numpy as np
import os
import glob
import shutil
import os
import subprocess


def findclosesttime(folder, targetdate):
    """
    Search folder for history file with time closest to the target date
    
    Args:
        folder:     folder holding output of a BESTNPZ ROMS simulation,
                    with history files matching the pattern *his*.nc.  
                    Alternatively, can be a list of history filenames 
                    (useful if you want to include a smaller subset from 
                    within a folder)
        targetdate: datenumber, target date
    
    Returns: 
        dictionary with the following keys/values:
            filename:   full path to history file including nearest date
            idx:        time index within that file (0-based) of nearest 
                        date
            dt:         timedelta, time between nearest date and target 
                        date
            unit:       time units used in history file
            cal:        calendar used by history file
    """
    
    if (type(folder) is str) and os.path.isdir(folder):
        hisfiles = glob.glob(os.path.join(*(folder, '*his*.nc')))
    else:
        hisfiles = folder

    f = nc.Dataset(hisfiles[0], 'r')
    tunit = f.variables['ocean_time'].units
    tcal = f.variables['ocean_time'].calendar

    dtmin = []
    
    d = {}
    for fn in hisfiles:
        try:
            f = nc.Dataset(fn, 'r')
            time = nc.num2date(f.variables['ocean_time'][:], units=tunit, calendar=tcal)

            dt = abs(time - targetdate)
    
            if not d:
                d['filename'] = fn
                d['idx'] = np.argmin(dt)
                d['dt'] = dt[d['idx']]
                d['time'] = time[d['idx']]
            else:
                if min(dt) < d['dt']:
                    d['filename'] = fn
                    d['idx'] = np.argmin(dt)
                    d['dt'] = dt[d['idx']]
                    d['time'] = time[d['idx']]
        except:
            pass
                
    d['unit'] = tunit
    d['cal'] = tcal
                
    return d
    
def fishtable(lsize, age_offset):
    """
    Creates a dictionary with dye parameters
    
    This function simply builds a dictionary that can be used as a lookup 
    table for various fish-related properies, such as match dye tracers 
    to fish species.
    
    Args:
        lsize:      list, size of length bins (cm) for each lengthed fish 
                    species (feast input variable: fsh_Lsize)
        age_offset: nested list of length bin offsets associated with 
                    each age of the age/length species (feast input 
                    variable: fsh_age_offset)
    
    Returns:
        dictionary with the following keys:
        gtype:      group type, either 'Age/length', 'Length', or 
                    'Simple'
        groupidx:   species index.  Indices follow python conventions, 
                    running from 0-14 for the 15 fish species
        name:       3-letter code for the species name
        lengthbin:  index of length bin (0-based).  A value of None 
                    implies that this property does not apply to this 
                    group.
        agebin:     index of age bin (0-based), currently equal to actual 
                    age (also in increments of 1 year starting at 0)  A 
                    value of None implies that this property does not 
                    apply to this group.
        variable:   variable being tracked by a given state variable.  
                    Can be 'Number', 'Condition factor', 'Caloric 
                    content', or 'Biomass'
        dye:        index of dye tracer (0-based).  Due to the 0-based 
                    nature of python indices, these will be off my one 
                    from the Fortran-like name (e.g. dye = 0 corresponds 
                    to dye_001)
        layer:      index of depth layer within a dye variable (0-based).
    """
                            
    # Number of each species type: aged-and-lengthed, lengthed, simple
                            
    naged = 3  # age-and-lengthed
    nleng = 7  # length-only
    nsimp = 5  # simple

    name_aged = ['POL', 'COD', 'ATF'] # TODO: would be nice to get full names
    name_leng = ['HER', 'CAP', 'EUL', 'SAN', 'MYC', 'SAL1', 'SAL2' ]
    name_simp = ['SHR', 'SQU', 'EPI', 'CRA', 'OTH']

    name = name_aged + name_leng + name_simp
    
    # Number of age and length bins

    nlength_leng = 20
    nlength_aged = 14

    nage = 11

    # Length bins for each fish-age combo
    
    ages = []
    for ii in range(naged):
        tmp = []
        for ia in range(nage):
            if ia == 0:
                lsz = lsize[ii]/2
            else:
                lsz = lsize[ii]
            tmp.append(np.arange(0,lsz*(nlength_aged+1),lsz) + age_offset[ii][ia+1]*lsz)
        ages.append(tmp)
        
    for ii in range(nleng):
        ages.append(np.arange(0,lsize[ii]*(nlength_leng+1),lsize[ii]))    
                        
    # Plankton

    nplank = 6 # TODO: COP, NCAS, NCA0, EUPS, EUPO, BEN...

    # Number of variables tracked per species type

    nvar_aged = 3 # Numbers, condition factor, caloric content
    nvar_leng = 3 # Numbers, condition factor, caloric content
    nvar_simp = 2 # Biomass, caloric content

    # Max indices in group types

    idxmax_aged = nvar_aged*naged*nlength_aged*nage
    idxmax_leng = nvar_leng*nleng*nlength_leng
    idxmax_simp = nvar_simp*nsimp

    # Dye variables and layers in ROMS

    ndye = 190
    nlayer = 10

    # Build table of details
    # Column 0: group type (0 = aged-lengthed, 1 = lengthed, 2 = simple)
    # Column 1: group index
    # Column 2: group name
    # Column 3: length bin index
    # Column 4: age bin index
    # Column 5: variable type 
    # Column 6: age bin edges

    vartype = ['Number', 'Condition factor', 'Calories', 'Biomass']

    vartable = []
    for ii in range(idxmax_aged):
        tmp = np.unravel_index(ii, (nvar_aged, naged, nlength_aged, nage), order='F')
        vidx = tmp[0]
        gidx = tmp[1]
        lidx = tmp[2]
        aidx = tmp[3]
        vartable.append(['Age/length', gidx, name_aged[gidx], lidx, aidx, vartype[vidx], list(ages[gidx][aidx][lidx:lidx+2])])
    
    for ii in range(idxmax_leng):
        tmp = np.unravel_index(ii, (nvar_leng, nleng, nlength_leng), order='F')
        vidx = tmp[0]
        gidx = tmp[1]
        lidx = tmp[2]
        vartable.append(['Length', gidx+naged, name_leng[gidx], lidx, None, vartype[vidx], list(ages[gidx+naged][lidx:lidx+2])])
    
    for ii in range(idxmax_simp):
        tmp = np.unravel_index(ii, (nvar_simp, nsimp), order='F')
        vidx = tmp[0]
        gidx = tmp[1]
        if vidx == 0:
            vtype = 'Biomass'
        else:
            vtype = 'Calories'    
        vartable.append(['Simple', gidx+naged+nleng, name_simp[gidx], None, None, vtype, None])
    
    for ii in range(len(vartable)):
        vartable[ii][7:8] = np.unravel_index(ii, (nlayer, ndye), order='F')
    
    # Convert to easier-to-query dictionary
    
    keys = ['gtype', 'groupidx', 'name', 'lengthbin', 'agebin', 'variable', 'length', 'layer', 'dye']
    d = {}
    for ii in range(len(keys)):
        d[keys[ii]] = [x[ii] for x in vartable]
    
    return d
    
def fishlookup(tbl, gtype=[], groupidx=[], name=[], lengthbin=[], agebin=[], variable=[], dye=[], layer=[], length=[]):
    """
    Return entries from the table that match all properties
    
    Args:
        tbl:    the fish lookup table (see fishtable())
    
    Optional key/value inputs:
        keys correspond to keys of tbl.  The value includes a list of 
        entries to match against.  For example, dye=[0,1] would return 
        the subset of tbl where dye values match 0 or 1.
    """
    
    ismatch = [all(tup) for tup in zip(keymatches(tbl, 'gtype', gtype), 
                                       keymatches(tbl, 'groupidx', groupidx),
                                       keymatches(tbl, 'name', name),
                                       keymatches(tbl, 'lengthbin', lengthbin), 
                                       keymatches(tbl, 'agebin', agebin),
                                       keymatches(tbl, 'variable', variable),
                                       keymatches(tbl, 'dye', dye), 
                                       keymatches(tbl, 'layer', layer),
                                       keymatches(tbl, 'dye', dye), 
                                       keymatches(tbl, 'layer', layer),
                                       keymatches(tbl, 'length', length)
                                        )]
    
    d = {}
    for ky in tbl:
        d[ky] = [tbl[ky][ii] for ii in range(len(ismatch)) if ismatch[ii]]
    return d
    
    
def keymatches(tbl, key, vals):
    if not isinstance(vals, list):
        vals = [vals]
    if vals:
        x = [ii in vals for ii in tbl[key]]
    else:
        x = len(tbl[key]) * [True]
    return x

def addfeastdyes(fname, ndye=190):
    """
    Add zeroed-out dye variables to a BESTNPZ initialization file
    
    Args:
        fname:  name of file to add dyes to
        ndye:   number of dyes.  Default is 190
    """
    
    f = nc.Dataset(fname, 'r+')
    sz = (f.dimensions['ocean_time'].size, f.dimensions['s_rho'].size, 
          f.dimensions['eta_rho'].size,    f.dimensions['xi_rho'].size)
    for id in range(1,ndye+1):
        dname = 'dye_{:03d}'.format(id)
        dvar = f.createVariable(dname, 'f4', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=1e37)
        dvar.units = 'kilogram meter-3'
        dvar.long_name = 'time-averaged dye concentration, type {}'.format(id)
        dvar.time = 'ocean_time'
        dvar.coordinates = 'lon_rho lat_rho s_rho ocean_time'
        dvar.field = '{}, scalar, series'.format(dname)
        dvar.missing_value = 1e37
        dvar[:] = np.zeros(sz)
    f.close()
        
def buildfeastbry(oldbry, newbry):
    """
    Creates boundary condition file for a FEAST run
    
    This function adds FEAST dye variables to a ROMS BEST_NPZ boundary 
    condition file.  The single dye variable will be used for all fish 
    boundary conditions, and is set to 0.
    
    Args:
        oldbry: full file name to BEST_NPZ boundary condition file
        newbry: full file name to new FEAST boundary condition file
    """

    shutil.copyfile(oldbry, newbry)
    f = nc.Dataset(newbry, 'r+')
    
    dyew = f.createVariable('dye_west_001', 'f4', ('ocean_time', 's_rho', 'eta_rho'))
    dyew.units = 'kilogram meter-3'
    dyew.missing_value = 1e36
    dyew.long_name = 'dye boundary condition'
    
    dyes = f.createVariable('dye_south_001', 'f4', ('ocean_time', 's_rho', 'xi_rho'))
    dyes.units = 'kilogram meter-3'
    dyes.missing_value = 1e36
    dyes.long_name = 'dye boundary condition'
     
    dyew[:] = np.zeros((f.dimensions['ocean_time'].size, f.dimensions['s_rho'].size, f.dimensions['eta_rho'].size))
    dyes[:] = np.zeros((f.dimensions['ocean_time'].size, f.dimensions['s_rho'].size, f.dimensions['xi_rho'].size))

    f.close()
    
def buildfeastini(years, hisfiles, inibase):
    """
    Creates initialization file for a FEAST simulation
    
    This function creates initialization files for single year runs 
    of ROMS with FEAST by adding dye variables to a BEST_NPZ history 
    file.  A FEAST simulation always simulates a single year, with a 
    5-month lead spinup period; this function looks through a set of 
    history files to find the appropriate starting time slice (closest to 
    July 1 of the year before the requested simulation years) and then 
    runs Kerim Aydin's addFishToIni.r R script.
    
    Note that this process currently relies on a very computer-specific
    path to the R script, which in turn relies on a computer-specific
    version of the netCDF library.)
    
    Args:
        years:      array of years for which to set up initialization 
                    files (one file per year)
        hisfiles:   array of full filenames of history file names from a 
                    BEST_NPZ simulation
        inibase:    base name for new FEAST initialization name (created 
                    files will be named inibase_YYYY.nc, where YYYY 
                    corresponds to the requested years.)
    """
    for yr in years:
        
        # Find history file from the designated BESTNPZ sim that is 
        # closest to July 1 of the previous year
        
        print(yr)
        print('  Extracting bestnpz time slice for July-before-Jan 1, {}'.format(yr))
        day1 = datetime(yr,1,1)
    
        prevjuly = datetime(day1.year-1, 7, 1)
        closest = findclosesttime(hisfiles, prevjuly)
        
        # Use ncks to slice out the appropriate time...

        inifile = '{}_{}.nc'.format(inibase, day1.year)

        cmd1 = ['ncks', '-O', '-d', 'ocean_time,{:d}'.format(closest['idx']), 
               closest['filename'], inifile]
        subprocess.run(cmd1)
        print('  Extracted to {}'.format(inifile))

        # ... and append the zeroed-out dyes

        print('  Appending dye variables')
        addfeastdyes(inifile)

    # Use Kerim's R script to add in all the appropriate fish data

    print('All years\n  Adding fish data to dye variables')

    env = os.environ.copy();
    env['LD_LIBRARY_PATH'] = ':'.join([env['LD_LIBRARY_PATH'], '/home/aydink/lib'])
    env['NETCDFHOME'] = '/home/aydink'

    cmd3 = ['/home/aydink/bin/Rscript', 'addFishToIni.r', ocean['GRDNAME'], inibase] + [str(x) for x in years]
     
    subprocess.call(cmd3, env=env)
    
    
        
    