TONIC
=====

[![Build Status](https://travis-ci.org/UW-Hydro/tonic.svg?branch=develop)](https://travis-ci.org/UW-Hydro/tonic)

`tonic` is a toolkit for working with distributed hydrologic models and their output.

This is truely a work in progress...

## Models:
1. [the Vairable Infiltration Capacity (VIC) model](https://github.com/UW-Hydro/VIC)
1.  [the Community Land Model (CLM)](http://www.cgd.ucar.edu/tss/clm/)
1.  [the Unified Land Model (ULM)](https://github.com/UW-Hydro/ULM)
1.  [the Precipitation Runoff Modeling System (PRMS)](http://wwwbrr.cr.usgs.gov/projects/SW_MoWS/PRMS.html)
1.  [the Noah Land Surface Model](http://www.ral.ucar.edu/research/land/technology/lsm.php)
1.  [the Structure for Unifying Multiple Modeling Alternatives (SUMMA)](http://www.ral.ucar.edu/projects/summa/)
1.  [FLO-2D](http://www.flo-2d.com/)

## Scripts:
`tonic` currently has 1 script available for use on the command line.

**vic_utils**: is a utility script that runs a number of VIC related processing tools.  Once `tonic` is installed, `vic_utils` will be available in your path.  Run `vic_utils -h` for a description of the utilities available.

## Install:
Dependencies:
- python 2.7 or greater
- netCDF4
- xray
- matplotlib
- basemap
- pandas

To install `tonic`, run `python setup.py install` from this directory.

Questions? Sure, Joe Hamman - jhamman at hydro.washington.edu

