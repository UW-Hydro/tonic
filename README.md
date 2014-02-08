VICpy
=====

A pre/post processing toolbox for the [VIC model](https://github.com/UW-Hydro/VIC)

Truely a work in progress...

This is a library of python scripts to go along with the pre/post
processing of VIC model forcings/parameters/outputs.

Right now this is just a collection of a few scripts but as the toolbox grows,
a more formal structure may be put together.

Most of the tools contained here are run from the command line and use a combination of command line arguments and configuration files.  For tool specific usage try raising the help flag, -h (e.g. vic2netcdf.py -h).

-----

### vic2netcdf.py
This tool converts binary or ASCII VIC output to a netCDF format.  Options and settings are input via configuration file (e.g. VIC2nc_example.cfg).

Usage: `vic2netcdf.py some_config_file.cfg`

-----

### netcdf2vic.py
This tool converts netCDF meteorological forcings to VIC binary or ASCII format.  Options and settings are input via configuration file (e.g. example_netcdf2vic.cfg).

Usage: `netcdf2vic.py some_config_file.cfg`

-----

### grid_params.py
This tool grids VIC soil, snow, and vegetation parameter files and stores them in a single netCDF file.  This allows for easy viewing of parameters using tools like [ncview](http://meteora.ucsd.edu/~pierce/ncview_home_page.html).  This tool take a list of command line arguments:

Usage:
For soil only: `grid_params.py --soil_file=global_soil_param_new -o outfile.nc`

For soil and veg: `grid_params.py -s global_soil_param_new -o outfile.nc -v global_lai_0.25deg.txt`

-----

Questions? Sure, Joe Hamman - jhamman at hydro.washington.edu

