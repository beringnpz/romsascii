#romsascii python package

This package includes a few simple tools to create the ascii input files needed for a ROMS simulation, as well as functions to run the simulation and parse the standard output files produced.  The goal is to allow all setup and running of a ROMS simulation to be done in a single python script, allowing clear documentation of any modifications to input parameters, and easy replication of simulations.

This is still very much a work in progress.

## Modules

###romsascii

Generic functions to create input files, run roms, and parse log files.  This module includes functions not tailored to any particular ROMS domain.  Primary functions that will be accessed by the user are:

- **filltimevars**: Fill in time-related variables based on  dates and time interval values
- **runroms**: Run ROMS, using the parameters in the input dictionaries
- **runromsthroughblowup**: Run ROMS, using the parameters in the input dictionaries, and try to get past blowups when they occur
- **writeromsascii**: Create ascii input file based on parameters in a dictionary
- **parseromslog**: Parse the standard output file from a ROMS simulation to see if it ran cleanly, if it hit a blowup, what the last timestep recorded was, and what the last history file written to was.

###defaultparams

This module includes functions to create default parameter dictionaries.  These parameter dictionaries form the basis for all the functions in this package.  Each dictionary key corresponds to a ROMS input variable.  These functions are tailored to our specific version of ROMS, though they could be expanded for newer versions and/or more generic parameters in the future.

- **ocean**:		standard input file (ocean.in) input parameters
- **bestnpz**:	biological parameters for the BEST_NPZ biology
- **feast**:		biological parameters for the FEAST biology

###bering10kinput

This module includes some very Bering 10K-specific functions.  They mostly serve as shortcuts so I can update .nc file locations and names all in one place.

- **ini\_frc\_bry**: locations of CORE and CFSR datasets
- **other\_in**: locations of varinfo.dat, grid file, and default stations, biology, and ice (the latter three will often be overwritten via parameter dictionaries in user scripts, but I provide this default just in case).
                         