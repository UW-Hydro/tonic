#!/usr/local/bin/python

import numpy as np
import struct
from netCDF4 import Dataset
import time as tm
import sys, os
import datetime
import argparse
import ConfigParser

##################################################################################
def main():
    (Files,Coord_Keys,Var_Keys,output,Binary_Mult,Binary_Type,Paths,outPrefix,verbose) = process_command_line()

    mask = read_netcdf(Paths['maskPath'],vars = ['mask'])['mask']
    yi,xi = np.nonzero(mask)
    print 'found', len(yi), 'points in mask file'
        
    for i,file in enumerate(Files):
        if i == 0:
            xlist = []
            ylist = []
            pointlist = []
            Append = False
            print Paths['inPath']+file
            d = read_netcdf(Paths['inPath']+file,verbose=verbose)
            xs = d['xc']
            ys = d['yc']
            posinds = np.nonzero(xs>180)
            xs[posinds] -= 360
            print 'adjusted xs lon minimum'
            for j in xrange(len(yi)):
                y,x = yi[j],xi[j]
                flag = 0
                for key in Var_Keys:
                    if (d[key][:,y,x].all() is np.ma.masked or mask[y,x]==0):
                        flag = 1
                if flag == 0:
                    point = (ys[y,x],xs[y,x])
                    xlist.append(x)
                    ylist.append(y)
                    pointlist.append(point)
                    data = np.empty((d[Var_Keys[0]].shape[0],len(Var_Keys)))
                    for j,key in enumerate(Var_Keys):
                        data[:,j] = d[key][:,y,x]
                        
                    if output['Binary']:
                        write_binary(data*Binary_Mult,point,Binary_Type,outPrefix,Paths,Append)
                    if output['ASCII']:
                        write_ASCII(data,point,outPrefix,Paths,Append)
        else:
            Append = True
            d = []
            data = []
            d = read_netcdf(Paths['inPath']+file, verbose=verbose)
            if verbose:
                print 'on file',str(i), 'point list is length: ',len(xlist)
            
            for k in xrange(len(xlist)):
                x = xlist[k]
                y = ylist[k]
                point = pointlist[k]
                data = np.empty((d[Var_Keys[0]].shape[0],len(Var_Keys)))
                for j,key in enumerate(Var_Keys):
                    data[:,j] = d[key][:,y,x]

                if output['Binary']:
                    write_binary(data*Binary_Mult,point,Binary_Type,outPrefix,Paths,Append)
                if output['ASCII']:
                    write_ASCII(data,point,outPrefix,Paths,Append)
                              
##################################################################################
def process_command_line():
    """
    Parse arguments and assign flags for further loading of variables, for
    information on input arguments, type rout.py -h
    """
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=str, help="Input Configuration File")
    args = parser.parse_args()
    if args.config:
        configFile = args.config
        (Files,Coord_Keys,Var_Keys,output,Binary_Mult,Binary_Type,Paths,outPrefix,verbose) = process_config(configFile)
    else:
        raise IOError('Need to supply configuration file, ./netcdf2vic.py configfile.cfg, type ./netcdf2vic.py -h for more info')

    return (Files,Coord_Keys,Var_Keys,output,Binary_Mult,Binary_Type,Paths,outPrefix,verbose)
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
    output = {'Binary':config.getboolean('Basics','Binary'),
              'ASCII':config.getboolean('Basics','ASCII')}
    Paths = {'inPath':config.get('Paths','inPath'),'maskPath':config.get('Paths','Mask')}
    File_str = config.get('Basics','File')

    if output['ASCII']:
        Paths['ASCIIoutPath'] = config.get('Paths','ASCIIoutPath')

    if output['Binary']:
        Binary_Mult = map(int,config.get('Basics','Mult').split(','))
        Binary_Type = config.get('Basics','Type')
        Paths['BinaryoutPath'] = config.get('Paths','BinaryoutPath')
    else:
        Binary_Mult = []
        Binary_Type = ""

    outPrefix = config.get('Basics','outPrefix')
    # Make list of files from wild cards
    if (WC_Start and WC_End):
        nums = np.arange(WC_Start,WC_End+1)
        temp = []
        for i in nums:
            temp.append(File_str.replace('**',str(i)))
        Files = temp
    
    return (Files,Coord_Keys,Var_Keys,output,Binary_Mult,Binary_Type,Paths,outPrefix,verbose)

##################################################################################
## Read netCDF Inputs
## Read data from input netCDF.
##################################################################################
def read_netcdf(ncFile,vars = [],coords = False, verbose = False):
    """
    Read data from input netCDF. Will read all variables if none provided.
    """
    f = Dataset(ncFile,'r')
    if vars==[]: vars = f.variables.keys()
    if verbose: print 'Reading input data vars:', vars, ', from file:',ncFile
    d={}
    if coords:
        if isinstance(vars,str):
            d[vars] = np.squeeze(f.variables[vars][coords])
        else:
            for var in vars:
                d[var] = np.squeeze(f.variables[var][coords])
    else:
        if isinstance(vars,str):
            d[vars] = np.squeeze(f.variables[vars][:])
        else:
            for var in vars:
                d[var] = np.squeeze(f.variables[var][:])
    f.close()
    return d

##################################################################################
def write_ASCII(array,point,outPrefix,Paths,Append,verbose = False):
    """
    Write an array to standard VIC ASCII output.
    """
    outFile = Paths['ASCIIoutPath']+outPrefix+('%1.3f' % point[0])+'_'+('%1.3f' % point[1])
    if Append:
        f = open(outFile, 'a')
    else:
        f = open(outFile, 'w')
    
    if verbose: print 'Writing ASCII Data to',outFile
    np.savetxt(f,array,fmt='%1.4f')
    f.close()

##################################################################################
def write_binary(array,point,Binary_Type,outPrefix,Paths,Append,verbose = False):
    """
    Write a given array to standard binary short int format.  
    """
    outFile = Paths['BinaryoutPath']+outPrefix+('%.3f' % point[0])+'_'+('%.3f' % point[1])
    if verbose: print 'Writing Binary Data to',outFile
    if Append:
        f = open(outFile, 'a')
    else:
        f = open(outFile, 'w')

    for row in array:
        data = struct.pack(Binary_Type,*row)
        f.write(data)
    f.close()

if __name__ == "__main__":
    main()    
