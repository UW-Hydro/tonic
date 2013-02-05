#!/usr/local/bin/python

import numpy as np
import struct
from netCDF4 import Dataset
import time as tm
import sys
import datetime

def main():
    outFile = 'data_'+??
    precipFile
    tempFile
    windFile
    outPrefix = 'data_'

def process_command_line():
    """
    Parse arguments and assign flags for further loading of variables, for
    information on input arguments, type rout.py -h
    """
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-A","--allFile", type=str, help="Input netCDF containing precip data")
    parser.add_argument("-P","--precipFile", type=str, help="Input netCDF containing precip data")
    parser.add_argument("-TH","--tmaxFile", type=str, help="Input netCDF containing tmax data")
    parser.add_argument("-TL","--tminFile", type=str, help="Input netCDF containing tmin data")
    parser.add_argument("-W","--windFile", type=str, help="Input netCDF containing wind data")
    parser.add_argument("-LON","--longitude", type=float, help="Input longitude of point to access data")
    parser.add_argument("-LAT","--latitude", type=float, help="Input latitude of point to access data")
    parser.add_argument("-v", "--verbose", help="Increase output verbosity", action="store_true")
    parser.add_argument("-C", "--config", help="Use Configuration File", action="store_true")
    args = parser.parse_args()
    if args.config:
        print 'Using configuration file'
    else:
        if args.allfile: precipFile = tmaxFile = tminFile = windFile = args.allFile
        if args.precipFile: precipFile = args.precipFile
        if args.tmaxFile: tmaxFile = args.tmaxFile
        if args.tminFile: tminFile = args.tminFile
        if args.windFile: windFile = args.windFile
        if (args.longitude and args.latitude):
            lon = args.longitude
            lat = args.latitude
            Full_Grid = False
            print 'Processing data from grid cell at lon/lat:',str(lon),',',str(lat)
        elif args.longitude:
            raise NameError("Need both -LON & -LAT args if one is specified, if neither are specified, the full grid will be processed into individual files")
        elif args.latitude:
            raise NameError("Need both -LON & -LAT args if one is specified,if neither are specified, the full grid will be processed into individual files")
        else:
            print "Processing all grid cells in input data set(s)"
            Full_Grid = True
    print 'Reading Precip data from:',precipFile
    print 'Reading Tmax data from:',tmaxFile
    print 'Reading Tmin data from:',tminFile
    print 'Reading Wind data from:',windFile
    
    
def process_config():
    """
    Parse arguments and assign flags for further loading of variables, for
    information on input arguments, type rout.py -h
    """

##################################################################################
## Read netCDF Inputs
## Read data from input netCDF.
##################################################################################
def read_netcdf(ncFile,vars,verbose):
    """
    Read data from input netCDF. Will read all variables if none provided.
    """
    f = Dataset(ncFile,'r')
    if vars==[]: vars = f.variables.keys()
    if verbose: print 'Reading input data vars:', vars, ',from file:',ncFile
    d={}
    for var in vars:
        d[var] = f.variables[var][:]
        f.close()
    return d

def write_ascii(outFile,array):
    """
    Write an array to standard VIC ASCII output.
    """
    np.savetxt(outFile,array)

def write_binary(outFile,array):
    """
    Write a given array to standard binary short int format.  
    """
    ---> add offsets
    np.save(outFile, array)
##################################################################################
##  Find Indicies
##  Given an input lat or lon, the function returns the nearest index location
##################################################################################
def find_nearest(array,value):
    """ Find the index location in (array) with value nearest to (value)"""
    idx = (np.abs(array-value)).argmin()
    return idx
    

if __name__ == "__main__":
    main()    
