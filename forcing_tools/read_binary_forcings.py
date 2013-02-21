#!/usr/local/bin/python

import numpy as np
import struct
from netCDF4 import Dataset
import time as tm
import sys
import datetime

def main():
    latlonf = '/raid/jhamman/raw_forcings/Adam2006/Arctic/latlon.txt'
    datapath = '/raid/jhamman/raw_forcings/Adam2006/Arctic/'
    outpath = '/raid/jhamman/raw_forcings/Adam2006/'
    NODATA = -9999
    gridstep = 0.5
    start_year = 1948
    start_month = 1
    start_day = 1

    # Get input file names and lat/lon locations
    print 'getting file list'
    lats,lons,files = load_file_list(latlonf)
    
    # make grid
    print 'making grids'
    x_grid,y_grid,xs,ys = make_grid(lats,lons,gridstep)

    #read in data
    for i,file in enumerate(files):
        #print 'reading data from file:',file[0]
        # find grid indicies to place data into
        x = find_nearest(xs,lons[i])
        y = find_nearest(ys,lats[i])
        #print 'placing into grid location', y, x
        # get data from file
        if i >0:
            precip[:,y,x],tmax[:,y,x],tmin[:,y,x],wind[:,y,x] = get_data(datapath+files[i][0])
            mask.append((x,y))
        else:
            p,t1,t2,w = get_data(datapath+files[i][0])
            #initialize target grids
            precip = np.empty((len(p),x_grid.shape[0],x_grid.shape[1]),dtype=np.uint16)
            tmin = np.zeros((len(p),x_grid.shape[0],x_grid.shape[1]),dtype=np.int16)
            tmax = np.zeros((len(p),x_grid.shape[0],x_grid.shape[1]),dtype=np.int16)
            wind = np.zeros((len(p),x_grid.shape[0],x_grid.shape[1]),dtype=np.int16)
            precip[:,y,x] = p
            tmin[:,y,x] = t1
            tmax[:,y,x] = t2
            wind[:,y,x] = w
            mask = [(x,y)]

    # Split by year and write to netCDF4
    len_time = precip.shape[0]
    time_series,year_series = time_fun(start_year,start_month,start_day,len_time)
    for year in np.unique(year_series):
        inds = np.nonzero(year==year_series)[0]
        times = time_series[inds]
        outfile = outpath+'Precip_Arctic_'+str(int(year))+'.nc'
        write_netcdf(outfile,xs,ys,times)
        write_var_netcdf(outfile,'precip','f4','mm',NODATA,transform(precip[inds,:,:],mask,NODATA,1/40.))

        outfile = outpath+'Data_Arctic'+str(int(year))+'.nc'
        write_netcdf(outfile,xs,ys,times)
        write_var_netcdf(outfile,'tmax','f4','celsius',NODATA,transform(tmax[inds,:,:],mask,NODATA,1/100.))
        write_var_netcdf(outfile,'tmin','f4','celsius',NODATA,transform(tmin[inds,:,:],mask,NODATA,1/100.))
        write_var_netcdf(outfile,'wind','f4','m/s',NODATA,transform(wind[inds,:,:],mask,NODATA,1/100.))
        print 'done with year', str(int(year))

##################################################################################
def load_file_list(latlonf):
    files = [i.strip().split() for i in open(latlonf).readlines()]
    # There's a better way to do this but I can't think of it now
    lats = np.zeros(len(files))
    lons = np.zeros(len(files))
    for i in xrange(len(files)):
        temp = files[i][0].split('_')
        lons[i] = float(temp[2])
        lats[i] = float(temp[1])
    return lats,lons,files
        
##################################################################################
def make_grid(lats,lons,gridstep):
    xs = np.arange(np.min(lons),np.max(lons)+gridstep,gridstep)
    ys = np.arange(np.min(lats),np.max(lats)+gridstep,gridstep)
    x_grid,y_grid = np.meshgrid(xs,ys)
    return x_grid,y_grid,xs,ys

##################################################################################
##  Find Indicies
##  Given an input lat or lon, the function returns the nearest index location
##################################################################################
def find_nearest(array,value):
    """ Find the index location in (array) with value nearest to (value)"""
    idx = (np.abs(array-value)).argmin()
    return idx

##################################################################################
def get_data(file):
    dt = np.dtype([('precip',np.uint16),('tmax',np.int16),('tmin',np.int16),('wind',np.int16)])
    Data = np.fromfile(file,dtype=dt)
    precip = Data['precip'].ravel()
    tmax = Data['tmax'].ravel()
    tmin = Data['tmin'].ravel()
    wind = Data['wind'].ravel()
    return precip,tmax,tmin,wind

def time_fun(start_year,start_month,start_day,len_time):
    start = datetime.datetime(start_year,start_month,start_day).toordinal()
    end = start + len_time
    time_series = np.arange(start,end)
    end_year = datetime.datetime.fromordinal(end).year
    year_series = np.zeros(len(time_series))
    for i,ord in enumerate(time_series):
        year_series[i] = datetime.datetime.fromordinal(ord).year
    time_series = time_series+1 # shift to days since 1-1-1
    return time_series,year_series

def transform(data,mask,NODATA,scale_factor):
    data_t = np.empty((data.shape))+NODATA
    for i in mask:
        x = i[0]
        y = i[1]
        data_t[:,y,x] = scale_factor * data[:,y,x].astype(np.float32)
    return data_t
    
def write_netcdf(outfile,lons,lats,times):
    """
    Write output to netCDF.  Writes out a netCDF4-64BIT data file containing
    vars defined above and a full set of history and description attributes.
    """
    f = Dataset(outfile,'w', format='NETCDF4')

    # set dimensions
    time = f.createDimension('time', (len(times)))
    lon = f.createDimension('lon', (len(lons)))
    lat = f.createDimension('lat', (len(lats)))

    # initialize variables
    time = f.createVariable('time','f8',('time',))
    lon = f.createVariable('lon','f8',('lon',))
    lat = f.createVariable('lat','f8',('lat',))

    # write attributes for netcdf
    f.description = 'Adam2006 Arctic data set in regular grid'
    f.citation1 = 'Adam, J. C., and D. P. Lettenmaier, 2003: Adjustment of global gridded precipitation for systematic bias. J. Geophys. Res., 108, 1-14'
    f.citation2 = 'Adam, J. C., E. A. Clark, D. P. Lettenmaier, and E. F. Wood, 2006: Correction of global precipitation products for orographic effects. J. Climate,  19, 15-38.'
    f.history = 'Created: {}\n'.format(tm.ctime(tm.time()))
    f.history += ' '.join(sys.argv) + '\n'
    f.source = sys.argv[0] # prints the name of script used
    f.user = 'jhamman@hydro.washington.edu'

    lat.long_name = 'latitude coordinate'
    lat.standard_name = 'latitude'
    lat.units = 'degrees_north'

    lon.long_name = 'longitude coordinate'
    lon.standard_name = 'longitude'
    lon.units = 'degrees_east'

    time.units = 'days since 1-1-1'
    
    # write data to variables initialized above
    time[:]= times
    lon[:] = lons
    lat[:] = lats

    f.close()

def write_var_netcdf(outfile,var_name,dtype,units,fill_value,data):
    """
    Write output to netCDF.  Writes out a netCDF4-64BIT data file containing
    vars defined above and a full set of history and description attributes.
    """
    f = Dataset(outfile,'a', format='NETCDF4')

    # initialize variables
    Var = f.createVariable(var_name,dtype,('time','lat','lon',),fill_value=fill_value,zlib=True)

    Var.units = units

    # write data to variables initialized above
    Var[:,:,:] = data

    f.close()

if __name__ == "__main__":
    main()
