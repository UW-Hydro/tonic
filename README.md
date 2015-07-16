TONIC
=====

`tonic` is a toolkit for working with distributed hydrologic models and their output.

This is truely a work in progress...

## Models:
1. [the Vairable Infiltration Capacity (VIC) model](https://github.com/UW-Hydro/VIC)
1.  [Community Land Model (CLM)](http://www.cgd.ucar.edu/tss/clm/)
1.  [Unified Land Model (ULM)](https://github.com/UW-Hydro/ULM)
1.  [Precipitation Runoff Modeling System (PRMS)](http://wwwbrr.cr.usgs.gov/projects/SW_MoWS/PRMS.html)
1.  [Noah Land Surface Model](http://www.ral.ucar.edu/research/land/technology/lsm.php)
1.  [Structure for Unifying Multiple Modeling Alternatives (SUMMA)](http://www.ral.ucar.edu/projects/summa/)
1.  [FLO-2D hydrodynamic floodplain model](http://www.flo-2d.com/)
1.  [SNOW17 accumulation and ablation model](http://www.nws.noaa.gov/oh/hrl/nwsrfs/users_manual/part2/_pdf/22snow17.pdf)

## Scripts:
`tonic` currently has 1 script available for use on the command line.

**vic_utils**: is a utility script that runs a number of VIC related processing tools.  Once `tonic` is installed, `vic_utils` will be available in your path.  Run `vic_utils -h` for a description of the utilities available.

## Install:
Dependencies:
- python 3
- netCDF4
- xray
- matplotlib
- basemap
- pandas

To install `tonic`, run `python setup.py install` from this directory.

Questions? Find us on [![Join the chat at https://gitter.im/UW-Hydro/tonic](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/UW-Hydro/tonic?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge).
