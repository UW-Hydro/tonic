#!/usr/local/bin/python

import numpy as np
import struct
from netCDF4 import Dataset
import time as tm
import sys, os
import datetime
import ConfigParser

def main(Files,Points,Full_Grid,output,Binary_Mult,Binary_Type,Paths,outPrefix,verbose):
    if Full_Grid:

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
    parser.add_argument("-v","--verbose", help="Increase output verbosity", default = False,action="store_true")
    parser.add_argument("-C","--config", help="Use Configuration File", action="store_true")
    parser.add_argument("-B","--Binary", help="output in binary format (default)",action="store_true")
    parser.add_argument("-T","--ASCII", help="output in ASCII format",action="store_true")
    args = parser.parse_args()
    if args.config:
        print 'Using configuration file'
        (Files,Points,Full_Grid,output,Binary_Mult,Binary_Type,Paths,verbose) = process_config(configFile)
    else:
        verbose = args.verbose
        if args.allfile: 
            Files = {'precip':[args.allFile],'tmax':[args.allFile],'tmin':[args.allFile],'wind':[args.allFile]}
        if args.precipFile: files['precip'] = [args.precipFile]
        if args.tmaxFile: files['tmax'] = [args.tmaxFile]
        if args.tminFile: files['tmin'] = [args.tminFile]
        if args.windFile: files['wind'] = [args.windFile]
        if (args.longitude and args.latitude):
            Points = [(args.longitude,args.latitude)]
            Full_Grid = False
            if verbose: print 'Processing data from grid cell at lon/lat:',str(lon),',',str(lat)
        elif args.longitude:
            raise NameError("Need both -LON & -LAT args if one is specified, if neither are specified, the full grid will be processed into individual files")
        elif args.latitude:
            raise NameError("Need both -LON & -LAT args if one is specified,if neither are specified, the full grid will be processed into individual files")
        else:
            print "Processing all grid cells in input data set(s)"
            Full_Grid = True
            Points = [()]
        if (args.Binary==False and args.ASCII==False):
            Output = {'Binary':True,'ASCII':False}
            if verbose: print 'No output data type was specified (-B or -T) so default will be used (Binary).'
        else: Output = {'Binary':args.Binary,'ASCII':args.ASCII}
        outPrefix = 'data_'
        Paths = {'inPath': os.getcwd, 'outPath':os.getcwd()}
        Binary_Mult = {'precip':40,'tmax':100,'tmin':100,'wind':100}
        Binary_Type = {'precip':np.uint16,'tmax':np.int16,'tmin':np.int16,'wind':np.int16}
    if verbose: print 'Reading Precip data from:',files
    return (Files,Points,Full_Grid,output,Binary_Mult,Binary_Type,Paths,outPrefix,verbose)
        
def process_config(configFile):
    """
    Parse arguments and assign flags for further loading of files.  Configuration
    flag must be raised on command line (-C) to configuration file must be provided.
    Usage:  netcdf2vic.py -C netcdf2vic.cfg
    """
    config = ConfigParser.ConfigParser()
    config.read(configFile)
    config.sections()
    keys = config.get('Basics','Var_Keys').split(',')
    verbose = config.getboolean('Basics','verbose')
    Full_Grid = config.getboolean('Basics','Full_Grid')
    output = {'Binary':config.getboolean('Basics','Binary'),
              'ASCII' = config.getboolean('Basics','ASCII')}
    Paths = {'inPath':config.get('Paths','inPath'),
             'outPath': config.get('Paths','outPath')}
    Binary_Mult = {}
    Binary_Type = {}
    Files = {}
    for key in keys:
        Binary_Mult[key] = config.getint(key,'Binary_Mult')
        Binary_Type[key] = config.get(key,'Binary_Type')
        Files[key] = config.get(key,'Files').split(',') 
    if Full_Grid == True:
        Points = [()]
    else:
        Points = zip(config.get('Points','Lons'),config.get('Points','Lats'))
    return (Files,Points,Full_Grid,output,Binary_Mult,Binary_Type,Paths,outPrefix,verbose)

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
