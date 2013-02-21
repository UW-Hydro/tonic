#!/usr/bin/python

import numpy as np
from netCDF4 import Dataset
from shapely.geometry.polygon import *
from shapely.geometry import *
from shapely.prepared import prep
import itertools
import pickle



TARGET_IN = '/raid/jhamman/temp_uh_files/run7/OBSAL_RASM_UH.nc'
SOURCE_IN = '/usr1/jhamman/Dropbox/RASM_Joe/pour_points/POP_pour_points.txt'
outfile = '/usr1/jhamman/Dropbox/RASM_Joe/pour_points/mapping.pkl'

def main():
    # Load inputs
    Input_TARGET = read_netcdf(TARGET_IN,('xc','xc_bnds','yc','yc_bnds'))
    (lons,lats) = load_points(SOURCE_IN)
    # Find Points within
    d = points_within(Input_TARGET['xc'], Input_TARGET['xc_bnds'],Input_TARGET['yc'], Input_TARGET['yc_bnds'],lons,lats)
    # Save dictionary to disk
    if not d:
        print 'dictionary is empty'
    else:
        save_dict(outfile,d)

##################################################################################
##  Find points within grid box
##  Return a list of all grid points inside box
##################################################################################
def points_within(xc, xc_bnds,yc,yc_bnds,lons,lats):
    print 'looking for points...'
    lons = lons+180
    cell=count=0
    points = MultiPoint(zip(lons,lats))
    length = len(points)
    cells = yc_bnds.shape[0]*yc_bnds.shape[1]
    d = {}
    for y_ind in xrange(yc.shape[0]):
        for x_ind in xrange(yc.shape[1]):
            coords = zip(xc_bnds[y_ind,x_ind,:],yc_bnds[y_ind,x_ind,:])
            coords.append(coords[0])
            polygon = Polygon(coords)
            prepared_polygon = prep(polygon)
            hits = filter(prepared_polygon.contains,points)
            if not hits:
                cell = cell+1
            else: 
                ### points.remove(hits) # so we don't keep looking for points already found
                plist = []
                for point in hits:
                    plist.append((point.x-180,point.y))
                d[(xc[y_ind,x_ind],yc[y_ind,x_ind])]=plist 
                count = count+len(hits)
                cell = cell +1
        print 'found target cell for', 100*count/length, '% of points'
        print '(looked in ', 100*cell/cells, '% of cells)'
    if count != length:
        print 'not all points found'
    return d

##################################################################################
##  Read netCDF Inputs
##  Read data from input netCDF.  Input netCDF should include the following vars:
##  ('Basin_ID','Flow_Direction','Flow_Distance','Land_Mask','lon','lat')
##  Velocity and Diffusion are assumed to be constant in this version
##################################################################################
def read_netcdf(nc_str,vars):
    print 'Reading input data vars:', vars
    f = Dataset(nc_str,'r')
    d={}
    for var in vars:
        d[var] = f.variables[var][:]
    f.close()
    return d

#################################################################################
## Read Pour Points
##   format: (lon,lat,basins,area,lon_min,lat_min, lon_max,lat_max)
#################################################################################
def load_points(infile):
    print 'Reading Pour Points from file: ', infile
    (lons,lats,basins,area,lon_min,lat_min, lon_max,lat_max) = np.genfromtxt(infile,delimiter=',',unpack=True)
    lons = lons
    lats = lats
    return (lons,lats)

#################################################################################
## Save dictionary
##   format: (lon,lat,basins,area,lon_min,lat_min, lon_max,lat_max)
#################################################################################
def save_dict(outfile,d):
    print 'Saving dictionary to: ', outfile
    f = open(outfile,'wb')
    pickle.dump(d,f)
    f.close()

main()
