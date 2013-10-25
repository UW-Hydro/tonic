#!/usr/bin/python
"""
This script calculates the minimum and maximum daily temperature from an input netcdf file
"""
import os
import glob
import tarfile
import numpy as np
from cdo import *
cdo = Cdo()

def main():

def process_command_line():
    """
    Parse arguments and assign flags for further loading of variables, for
    information on input arguments, type rout.py -h
    """
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", type=str, help="Input file")
    parser.add_argument("--outdir", type=str, help="Output file to append tmin/tmax to, if not specified a new file will be written")
    args = parser.parse_args()

    infile = args.infile
    try:
        outfile = args.outfile
    else:
        outfile = False
    return infile, outfile


def read_netcdf(ncFile,vars = [],coords = False, verbose = False):
    """
    Read data from input netCDF. Will read all variables if none provided.
    Will also return all variable attributes.
    Both variables (data and attributes) are returned as dictionaries named by variable
    """
    if verbose:
        print 'Reading input data vars:', vars, ', from file:',ncFile
    f = Dataset(ncFile,'r')
    if vars==[]: vars = f.variables.keys()
    d={}
    a={}
    g={}
    if coords:
        if isinstance(vars,str):
            d[vars] = f.variables[vars][coords]
            a[vars] = f.variables[vars].__dict__
        else:
            for var in vars:
                d[var] = f.variables[var][coords]
                a[var] = f.variables[var].__dict__
    else:
        if isinstance(vars,str):
            d[vars] = f.variables[vars][:]
            a[vars] = f.variables[vars].__dict__
        else:
            for var in vars:
                d[var] = f.variables[var][:]
                a[var] = f.variables[var].__dict__
    
    for attr in f.ncattrs():
        g[attr] = getattr(f,attr)
    
    f.close()
    
    return d,a,g

##################################################################################
# Run Program
##################################################################################
if __name__ == "__main__":
    main()
