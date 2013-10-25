#!/usr/local/bin/python

import numpy as np
from netCDF4 import Dataset
import time as tm
import sys, os, shutil, glob, fnmatch
import datetime
import argparse, ConfigParser
from scipy.spatial import cKDTree


def main():
    # Get command line arguments
    targetFile, inputPath, variables, inputPrefix, outputPrefix,output,period = process_command_line()

    # Get target grid information
    gridData, gridAttrs = read_netcdf(targetFile)

    # Get list of files and lat/lon coordinates
    files,lats,lons,prefixes = get_files(inputPath,inputPrefix)

    # Get Grid inds
    yinds,xinds = get_grid_inds(gridData,lats,lons)

    # Look in first file and get the time series info (split by month)
    time, step_names,tInds = get_dates(files[0],output)

    # Make netcdf files
    make_files(tInds,yinds,xinds,step_names,time,targetFile,files,variables,outputPrefix,inputPath,period)    

def make_files(tInds,yinds,xinds,step_names,time,targetFile,files,variables,outputPrefix,inputPath,period):

    if period:
        period_num = [i for i,x in enumerate(step_names) if x == period]
        print period_num
        if period_num == []:
            raise RuntimeError('Period string (', period , ') did not match any of the possible strings:',step_names)
        else:
            step_names = [step_names[period_num[0]]]
            print 'output step name is: ',step_names

    for i in xrange(len(step_names)):
        newFile = inputPath+outputPrefix+step_names[i]+'.nc'
        if period:
            inds = tInds==period_num
        else:
            inds = tInds==i        
        times = time[inds]
        #start new netcdf file
        f,d = start_netcdf(targetFile,newFile,variables,times)
        
        for j, vicFile in enumerate(files):
            y,x = yinds[j],xinds[j]
            data = np.loadtxt(vicFile)
            #Loop over all variables, slicing from data and placing in netcdf file
            for variable in variables:
                col = variables[variable]["column"]
                d[variable][:,y,x] = data[inds,col]
        f.close()
    return

def start_netcdf(targetFile,newFile,variables,times):
    """
    Copy target grid and return file handle.
    Add times dimension and data (this wont change for any of the inputs to come
    """
    shutil.copyfile(targetFile,newFile)
    f = Dataset(newFile,'a')
    # Add time dimension
    time = f.createDimension('time', len(times))
    time = f.createVariable('time','f8',('time',))
    time.units = 'days since 1-1-1'
    time[:]= times

    # Add variables and attributes
    d = {}
    for variable in variables:
        d[variable] = f.createVariable(variable,'f8',('time','nj','ni',))
        d[variable].description = variables[variable]["description"]
        d[variable].units = variables[variable]["units"]
        d[variable].coordinates = variables[variable]["coordinates"]

    f.description = 'VIC output'
    f.history = 'Created: {}\n'.format(tm.ctime(tm.time()))
    f.history += ' '.join(sys.argv) + '\n'
    f.source = sys.argv[0] # prints the name of script used

    print 'started netcdf',newFile

    return f,d

def process_command_line():
    """
    Process command line arguments. Must have target grid (in netcdf format) and inputPath where VIC outputs are located
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("targetFile", type=str, help="Input netCDF target grid")
    parser.add_argument("configFile",type=str,help="Input configuration file")
    parser.add_argument("inputPath", type=str, help="Input path where VIC outputs are located")
    parser.add_argument("inputPrefix",type=str, help="Input prefix out VIC output files")
    parser.add_argument("-OP","--outputPrefix",type=str, help="Output file prefix", default="VIC.offline.")
    parser.add_argument('--output',type=str,help='Type of output file (options are monthly, yearly, or full).  Yearly is the default.', default='yearly')
    parser.add_argument('--period',type=str,help='Limit the output to only one time period', default=False)
    args = parser.parse_args()
    targetFile = args.targetFile
    inputPath = args.inputPath
    inputPrefix = args.inputPrefix
    outputPrefix = args.outputPrefix
    configFile = args.configFile
    output = args.output
    period = args.period

    variables = process_config_file(configFile)

    return targetFile, inputPath, variables, inputPrefix,outputPrefix,output,period

def process_config_file(configFile):
    config = ConfigParser.ConfigParser()
    config.read(configFile)
    vars = config.sections()
    variables = {}
    for var in vars:
        variables[var] = ConfigSectionMap(config,var)

    print 'processed configuration file and found this data:', variables

    return variables

def ConfigSectionMap(config,section):
    dict1 = {}
    options = config.options(section)
    for option in options:
        try:
            dict1[option] = config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1

def get_files(inputPath,inputPrefix):
    """
    Get list of VIC output file:
    Return a list of file names, a list of corresponing prefixes, and arrays of lats and lons
    """
    names = os.listdir(inputPath)
    files = fnmatch.filter(names, inputPrefix+"*")
    lats = np.zeros(len(files))
    lons = np.zeros(len(files))

    prefix = []
    for i in xrange(len(files)):
        temp = files[i].split('_')
        lons[i] = float(temp[2])
        lats[i] = float(temp[1])
        prefix.append(temp[0])

    prefixes = list(set(prefix))
    
    files = glob.glob(inputPath+inputPrefix+"*")

    print 'found',len(files),'files in',inputPath

    return files,lats,lons,prefixes

def get_dates(file,output='yearly'):
    """
    Read the first file in the input directory and create a ordinal based timeseries
    Also find the indicies to split the time series into months and years
    """
    data = np.loadtxt(file)
    sYear,sMonth,sDay,sHour = data[0,0:4].astype(int)
    eYear,eMonth,eDay,eHour = data[-1,0:4].astype(int)
    # make a timeseries of dates between these dates
    start = datetime.datetime(sYear,sMonth,sDay,sHour).toordinal()
    end = float(datetime.datetime(eYear,eMonth,eDay,eHour).toordinal())+eHour/24.
    time = np.linspace(start,end,num=len(data))+1
    
    # split between months
    years,months = data[:,0].astype(int),data[:,1].astype(int)
    tInds = np.zeros(len(years),dtype=int)
    step = 0
    step_names = []

    # if we want monthly output
    if output == 'monthly':
        for y in np.unique(years):
            for m in np.unique(months):
                inds = np.intersect1d(np.nonzero(years==y)[0],np.nonzero(months==m)[0])
                tInds[inds] = step
                step_names.append("%i-%02d" % (y,m))
                step += 1

    # If we want yearly output
    elif output == 'yearly':
        for y in np.unique(years):
            inds = np.nonzero(years==y)[0]
            tInds[inds] = step
            step_names.append("%i" % y)
            step += 1

    # If we want a single output file
    else:
        step_names.append('full')

    print 'made a timeseries of length', len(time), 'and prepared inds for these output files', step_names

    return time,step_names,tInds
        
##################################################################################
## Read netCDF Inputs
## Read data from input netCDF.
##################################################################################
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

def get_grid_inds(targetGrid,lats,lons):
    """
    Find location of lat/lon points in 2d target grid.
    Uses cKdtree nearest neighbor mapping.
    """
    if (lons.min()<0 and targetGrid['xc'].min()>=0):
        posinds = np.nonzero(lons<0)
        lons[posinds] += 360
        print 'adjusted VIC lon minimum (+360 for negative lons)'

    combined = np.dstack(([targetGrid['yc'].ravel(),targetGrid['xc'].ravel()]))[0]
    points = list(np.vstack((lats,lons)).transpose())

    mytree = cKDTree(combined)
    dist, indexes = mytree.query(points,k=1)

    yinds,xinds = np.unravel_index(indexes,(targetGrid['xc'].shape[0],targetGrid['xc'].shape[1]))

    print 'got x and y inds'
    return yinds, xinds

if __name__ == "__main__":
    main()    
