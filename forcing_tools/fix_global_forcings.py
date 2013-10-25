#!/usr/local/bin/python

import numpy as np
from netCDF4 import Dataset
import glob
import os
import sys
import time as tm
from cdo import *
cdo = Cdo()

grid = '/raid/jhamman/RASM_masks/domain.lnd.wr50a_ar9v4.100920.nc'
inpath = '/raid/jhamman/raw_forcings/Adam2006/Global/all/'
outpath = '/raid/jhamman/RASM_met_forcings/Adam2006_Global/'
tempfile = 'temp.nc'

def main():
    os.chdir(inpath)
    files = glob.glob("*nc")
    for i,infile in enumerate(files):
        fix_netcdf(infile,tempfile)
        
        outfile1 = outpath+infile[:-3]+'_RASM_BIL.nc'
        outfile2 = outpath+infile[:-3]+'_RASM_CON.nc'

        print 'trying input:', infile,tempfile, 'output:', outfile1,outfile2
        cdo.remapbil(grid,input = tempfile,output = outfile1)

        cdo.remapbil(grid,input = tempfile,output = outfile2)

        
        print 'done with outfile', outfile2, '(',i,'of',len(files),')'


##################################################################################
## Read netCDF Inputs
## Read data from input netCDF.
##################################################################################
def read_netcdf(ncFile,vars = [],coords = False, verbose = False):
    """
    Read data from input netCDF. Will read all variables if none provided.
    Will also return all variable attributes.
    Both variables (data and attributes) are returned as dictionaries named by variable
    """
    #if verbose: print 'Reading input data vars:', vars, ', from file:',ncFile
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

def fix_netcdf(infile,outfile):
    """
    Write a new netcdf but this time do the coordinate vars correctly
    """
    rootgrp = Dataset(outfile,'w', format='NETCDF3_64BIT')

    data, targetAttrs = read_netcdf(infile,vars=('Prec','Wind','Tmax','Tmin','time','nav_lat','nav_lon'))
    res = 0.5
    # set dimensions
    lon = rootgrp.createDimension('lon',data['Prec'].shape[2])
    lat = rootgrp.createDimension('lat',data['Prec'].shape[1])
    time = rootgrp.createDimension('time',data['Prec'].shape[0])

    # do vars
    
    times = rootgrp.createVariable('time','f8',('time',))
    times[:] = np.arange(data['Prec'].shape[0])*86400
    times.units = targetAttrs['time']['units']
    times.long_name = targetAttrs['time']['long_name']
    
    lat = rootgrp.createVariable('lat','f8',('lat',))
    lat[:] = np.arange(data['nav_lat'].min(),data['nav_lat'].max()+res,res)
    lat.units = 'degrees_north'
    lat.long_name = 'Latitude'
    
    lon = rootgrp.createVariable('lon','f8',('lon',))
    lon[:] = np.arange(data['nav_lon'].min(),data['nav_lon'].max()+res,res)
    lon.units = 'degrees_east'
    lon.long_name = 'Longitude'
    
    Precip = rootgrp.createVariable('Precip','f8',('time','lat','lon',),fill_value=data['Prec'].fill_value)
    Precip[:,:,:] = data['Prec']
    Precip.units = targetAttrs['Prec']['units']
    Precip.long_name = targetAttrs['Prec']['long_name']
    
    Tmax = rootgrp.createVariable('Tmax','f8',('time','lat','lon',),fill_value=data['Tmax'].fill_value)
    Tmax[:,:,:] = data['Tmax']
    Tmax.units = targetAttrs['Tmax']['units']
    Tmax.long_name = targetAttrs['Tmax']['long_name']

    Tmin = rootgrp.createVariable('Tmin','f8',('time','lat','lon',),fill_value=data['Tmin'].fill_value)
    Tmin[:,:,:] = data['Tmin']
    Tmin.units = targetAttrs['Tmin']['units']
    Tmin.long_name = targetAttrs['Tmin']['long_name']    

    Wind = rootgrp.createVariable('Wind','f8',('time','lat','lon',),fill_value=data['Wind'].fill_value)
    Wind[:,:,:] = data['Wind']
    Wind.units = targetAttrs['Wind']['units']
    Wind.long_name = targetAttrs['Wind']['long_name']
    
    rootgrp.description = 'Global 1/2 Degree Gridded Meteorological VIC Forcing Data Set '
    rootgrp.history = 'Created: {}\n'.format(tm.ctime(tm.time()))
    rootgrp.source = sys.argv[0] # prints the name of script used
    rootgrp.institution = "University of Washington Dept. of Civil and Environmental Engineering"
    rootgrp.sources = "UDel (Willmott and Matsuura 2007), CRU (Mitchell et al., 2004), NCEP/NCAR (Kalnay et al. 1996)"
    rootgrp.projection = "Geographic"
    rootgrp.surfSng_convention = "Traditional"

    rootgrp.close()

if __name__ == "__main__":
    main()



    

