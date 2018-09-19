#romsascii python package

This package is designed to allow for programmatic organization of the many, many input variables associated with a ROMS simulation.  In typical ROMS use, these parameters are stored in an input file passed to the ROMS executable (typically named ocean.in, or some variant on that).  Additional input files are defined by variables within the ocean.in file itself (for example, the SPOSNAM variable holds the path to a file with all the station-related input variables.)

All of the modules in this package are based on the concept of storing these input variables in a python dictionary rather than in static input files, and dynamically generating the properly-formatted input files on a per-simulation basis.  The goal is to allow all setup and running of ROMS simulation to be done in a single python script, allowing clear documentation of any modifications to input parameters, and easy replication of simulations without having to keep track of too many separate files.

Note that this package focuses only on ascii input variables, *NOT* netCDF inputs.  There are other tools available (e.g. Pyroms: https://github.com/ESMG/pyroms) that are designed to help with that part of the setup process.

This is still very much a work in progress.

## Modules

###romsascii

Generic functions to create input files, run roms, and parse log files. This module includes functions not tailored to any particular ROMS domain.  Primary functions that will be accessed by the user are:

- **writeromsascii**: Dynamically generate ascii input files based on the parameters in a ROMS parameter dictionary
- **runroms**: This function dynamically generates the ascii input files needed for a ROMS simulation (i.e. ocean.in) based on a ROMS parameter dictionary and then issues the mpirun command needed to start the run.
- **filltimevars**: Fill in time-related variables based on  dates and time interval values, eliminating the need for manual calculation of number of time steps in variables such as NTIMES, NDEFHIS, etc. 
- **parseromslog**: Parse the standard output file from a ROMS simulation to see if it ran cleanly, if it hit a blowup, what the last timestep recorded was, and what the last history file written to was.
- **parserst**: Finds the name of and parses the simulation counter from a set of restart files, assuming files were created using the naming convention in runromssmart and its successors (i.e. filebase\_XX\_rst.nc).
- **runromssmart**: DEPRECATED This was intended to be a wrapper around runroms in order to take care of some of the setup and cleanup work for a ROMS simulation. It turned out that this sort of thing was pretty computer- and simulation-specific (i.e. different for hindcasts, forecasts, spinup loops, etc.), so instead I've replaced this with the more application-specific functions in the moxroms module.

###defaultparams

This module includes functions to create default parameter dictionaries.  These parameter dictionaries form the basis for all the functions in this package.  Each dictionary key corresponds to a ROMS input variable.  These functions are tailored to our specific version of ROMS, but they could be expanded for newer versions and/or more generic parameters in the future.

- **ocean**: standard input file (ocean.in) input parameters
- **bestnpz**: biological parameters for the BEST_NPZ biology
- **feast**: biological parameters for the FEAST biology
- **ice**: ice parameters
- **stations**: station parameters

###bering10kinput

DEPRECATED  File naming conventions vary between simulation types and computers, so I've replaced this with functions in the moxroms module for now,

### moxroms

This module handles some of the very computer- and simulation-specific stuff that I tend to do over and over.  Right now, that includes:

- **setinfiles**: sets the GRDNAME, VARNAME, ININAME, FRCNAME, BRYNAME, and NFFILES parameters for a Bering 10K simulation on hyak-mox, automatically choosing the appropriate ones based on the time in the specified initialization file and a few other parameters (# layers, break year for CORE-to-CFS transition).
- **runhindcast**: Runs a simulation with CORE-CFS hindcast input, taking care of all the messiness of switching up input files, checking for blowups, reducing and increasing timesteps, and resuming partially-completed runs.


## To build

To rebuild in development mode:

`python setup.py develop`

(or, if you want/need to write to a user-specific location rather than the system default location)

`python setup.py develop --user`

To rebuild in final mode:

`python setup.py install [--user]`