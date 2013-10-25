#!/usr/local/bin/python

import numpy as np
from netCDF4 import Dataset
import argparse

def main():

    targetFile,outFile = process_command_line()
    d,a = read_netcdf(targetFile)

    cols,rows = np.nonzero(d['mask']>0)

    cell_ids = np.arange(d['xc'][cols,rows].size)
    lons = np.ravel(d['xc'][cols,rows])
    lons[lons>180] -= 360
    lats = np.ravel(d['yc'][cols,rows])

    out = np.empty((len(cell_ids),5))
    out[:,0] = cell_ids
    out[:,1] = rows
    out[:,2] = cols
    out[:,3] = lats
    out[:,4] = lons

    print 'writing to outfile:',outFile
    np.savetxt(outFile,out,delimiter=' ', fmt=['%1.1i','%1.1i','%1.1i','%1.3f','%1.3f'])

    return

def process_command_line():
    """
    Process command line arguments. Must have target grid (in netcdf format) and inputPath where VIC outputs are located
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("targetFile", type=str, help="Input netCDF target grid")
    parser.add_argument("outFile", type=str, help="Input netCDF target grid")

    args = parser.parse_args()
    targetFile = args.targetFile
    outFile = args.outFile

    return targetFile,outFile

def read_netcdf(ncFile,vars = [],coords = False, verbose = True):
    """
    Read data from input netCDF. Will read all variables if none provided.
    Will also return all variable attributes.
    Both variables (data and attributes) are returned as dictionaries named by variable
    """
    if verbose: print 'Reading input data vars:', vars, ', from file:',ncFile
    f = Dataset(ncFile,'r')
    if vars==[]: vars = f.variables.keys()
    d={}
    a={}
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
    f.close()

    return d,a

if __name__ == "__main__":
    main()
