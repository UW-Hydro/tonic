#!/usr/local/bin/python

import numpy as np
import struct
from netCDF4 import Dataset
import time as tm
import sys, os
import datetime
import argparse
import ConfigParser

"""
left to do...
1. Figure out how to write to binary
2.  

"""

##################################################################################
def main():
    (Files,Coord_Keys,Var_Keys,Points,Full_Grid,output,Binary_Mult,Binary_Type,Paths,outPrefix,verbose) = process_command_line()
    file = Paths['inPath']+Files[Files.keys()[0]][0]
    xs,ys = get_points_list(file,Coord_Keys)
    for y in xrange(len(ys)):
        for x in xrange(len(xs)):
            point = (xs[y,x],ys[y,x])
            print x,y,point
            flag = 0
            data = []
            for i,key in enumerate(Var_Keys):
                temp = []
                for j,f in enumerate(Files[key]):
                    file = Paths['inPath']+f
                    # get all the data for that variable (multiple files if need be)
                    if flag == 0:
                        d = read_netcdf(file,vars = str(key),coords = [slice(None),y,x],verbose=True)
                    if (j == 0 and flag == 0):
                        if d[key].all() is np.ma.masked:
                            flag = 1
                            #print 'skipping', point, '(no data here)'
                    if (j > 0 and flag == 0):
                        temp = np.append(temp,d[key],axis=0)
                    elif flag==0:
                        temp = d[key]
                if (i == 0 and flag == 0):
                    data = np.empty((len(temp),len(Var_Keys)))
                    data[:,i] = temp[:]
                elif flag==0:
                    data[:,i] = temp[:]
            if flag == 0:
                if output['Binary']:
                    write_binary(data*Binary_Mult,point,Binary_Type,outPrefix,Paths,verbose)
                if output['ASCII']:
                    write_ASCII(data,point,outPrefix,Paths,verbose)
                              
##################################################################################
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
    parser.add_argument("-C","--config", type=str, help="Use Configuration File")
    parser.add_argument("-B","--Binary", help="output in binary format (default)",action="store_true")
    parser.add_argument("-T","--ASCII", help="output in ASCII format",action="store_true")
    args = parser.parse_args()
    if args.config:
        print 'Using configuration file'
        configFile = args.config
        (Files,Coord_Keys,Var_Keys,Points,Full_Grid,output,Binary_Mult,Binary_Type,Paths,outPrefix,verbose) = process_config(configFile)
    else:
        verbose = args.verbose
        if args.allFile: 
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
        Paths = {'inPath': os.getcwd, 'ASCIIoutPath':os.getcwd(),'BinaryoutPath':os.getcwd}
        Binary_Mult = {'precip':40,'tmax':100,'tmin':100,'wind':100}
        Binary_Type = "Hhhh"
    return (Files,Coord_Keys,Var_Keys,Points,Full_Grid,output,Binary_Mult,Binary_Type,Paths,outPrefix,verbose)

##################################################################################  
def process_config(configFile):
    """
    Parse arguments and assign flags for further loading of files.  Configuration
    flag must be raised on command line (-C) to configuration file must be provided.
    Usage:  netcdf2vic.py -C netcdf2vic.cfg
    """
    config = ConfigParser.ConfigParser()
    config.read(configFile)
    config.sections()
    Var_Keys = config.get('Basics','Var_Keys').split(',')
    Coord_Keys = config.get('Basics','Coord_Keys').split(',')
    verbose = config.getboolean('Basics','verbose')
    WC_Start = config.getint('Basics','WC_Start')
    WC_End = config.getint('Basics','WC_End')
    Full_Grid = config.getboolean('Basics','Full_Grid')
    output = {'Binary':config.getboolean('Basics','Binary'),
              'ASCII':config.getboolean('Basics','ASCII')}
    Paths = {'inPath':config.get('Paths','inPath'),
             'ASCIIoutPath':config.get('Paths','ASCIIoutPath'),
             'BinaryoutPath':config.get('Paths','BinaryoutPath')}
    Binary_Mult = np.zeros(len(Var_Keys))
    Binary_Type = ""
    Files = {}
    for i,key in enumerate(Var_Keys):
        Binary_Mult[i] = config.getint(key,'Binary_Mult')
        Binary_Type += config.get(key,'Binary_Type')
        Files[key] = config.get(key,'Files').split(',') 
    if Full_Grid == True:
        Points = [()]
    else:
        Points = zip(config.get('Points','Lons'),config.get('Points','Lats'))
    outPrefix = config.get('Basics','outPrefix')
    # Make list of files from wild cards
    if (WC_Start and WC_End):
        nums = np.arange(WC_Start,WC_End+1)
        for key in Var_Keys:
            temp = []
            for i in nums:
                f = Files[key][0]
                temp.append(f.replace('**',str(i)))
            Files[key] = temp
    return (Files,Coord_Keys,Var_Keys,Points,Full_Grid,output,Binary_Mult,Binary_Type,Paths,outPrefix,verbose)

##################################################################################
## Read netCDF Inputs
## Read data from input netCDF.
##################################################################################
def read_netcdf(ncFile,vars = False,coords = False, verbose = False):
    """
    Read data from input netCDF. Will read all variables if none provided.
    """
    if verbose: print 'Reading input data vars:', vars, ', from file:',ncFile
    f = Dataset(ncFile,'r')
    if vars==[]: vars = f.variables.keys()
    d={}
    if coords:
        if isinstance(vars,str):
            d[vars] = f.variables[vars][coords]
        else:
            for var in vars:
                d[var] = f.variables[var][coords]
    else:
        if isinstance(vars,str):
            d[vars] = f.variables[vars][:]
        else:
            for var in vars:
                d[var] = f.variables[var][:]
    f.close()
    return d
##################################################################################
def get_points_list(ncFile,Coord_Keys,verbose = False):
    coords = read_netcdf(ncFile,vars = Coord_Keys)
    if coords[Coord_Keys[0]].ndim == 1:
        x = coords[Coord_Keys[0]]
        y = coords[Coord_Keys[1]]
        #num = len(x)*len(y)
        xs,ys = np.meshgrid(x,y)
    elif coords[Coord_Keys[0]].ndim == 2:
        xs = coords[Coord_Keys[0]]
        ys = coords[Coord_Keys[1]]

    else:
        raise NameError("Coord_Keys given provided coordinate variable with an inapropriate number of dimensions",coords[Coord_Keys[0]].ndim,'was expecting 1 or 2')
    return xs,ys

##################################################################################
def write_ASCII(array,point,outPrefix,Paths,verbose):
    """
    Write an array to standard VIC ASCII output.
    """
    dt = np.dtype([('precip',np.uint16),('tmax',np.int16),('tmin',np.int16),('wind',np.int16)])
    outFile = Paths['ASCIIoutPath']+outPrefix+('%.3f' % point[0])+'_'+('%.3f' % point[1])
    if verbose: print 'Writing ASCII Data to',outFile
    np.savetxt(outFile,array,fmt='%1.4f')

##################################################################################
def write_binary(array,point,Binary_Type,outPrefix,Paths,verbose):
    """
    Write a given array to standard binary short int format.  
    """
    outFile = Paths['BinaryoutPath']+outPrefix+('%.3f' % point[0])+'_'+('%.3f' % point[1])
    if verbose: print 'Writing Binary Data to',outFile
    f = open(outFile, 'w')

    for row in array:
        data = struct.pack(Binary_Type,*row)
        f.write(data)
    f.close()
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
