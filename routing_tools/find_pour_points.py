#!/usr/bin/python

def read_netcdf(nc_str,vars):
    print 'Reading input data vars:', vars
    f = Dataset(nc_str,'r')
    d={}
    for var in vars:
        d[var] = f.variables[var][:]
    f.close()
    return d

###################
def Find_Points(areas,basins,lons,lats):
    #Make arrays of same dimensions as input arrays of lat/lon values
    x,y = np.meshgrid(lons, lats)
    #Setup basin in/out arrays
    basin_ids = np.arange(np.min(basins),np.max(basins))
    #Test basins = old...basin_ids = np.array([57087,53728,55400,37250,60966,71339])
    #basin_ids = np.array([46049,35684,49218,51765])
    basin = np.zeros(len(basin_ids),dtype='i')
    max_area =  np.zeros(len(basin_ids),dtype='i')
    x_outlet = np.zeros(len(basin_ids),dtype='f')
    y_outlet = np.zeros(len(basin_ids),dtype='f')
    min_x = np.zeros(len(basin_ids),dtype='f')
    max_x = np.zeros(len(basin_ids),dtype='f')
    min_y = np.zeros(len(basin_ids),dtype='f')
    max_y = np.zeros(len(basin_ids),dtype='f')
    #Loop over every basin id, finding the maximum upstream area [and location]
    #and record the basin#,longitude,latitude,area
    length = len(basin_ids)
    for i,j in enumerate(basin_ids):
        print 'basin:', i, ' of:',length
        basin[i]=np.int(j)
        inds = np.nonzero(basins==j)
        x_basin = x[inds]
        y_basin = y[inds]
        max_area[i] = np.int(max(areas[inds]))
        max_ind = np.argmax(areas[inds])
        x_outlet[i] = x_basin[max_ind]
        y_outlet[i] = y_basin[max_ind]
        min_x[i] = min(x_basin)
        max_x[i] = max(x_basin)
        min_y[i] = min(y_basin)
        max_y[i] = max(y_basin)
    #save the list of pour points as a comma seperated text file
    # This file is directly importable into arcgis for validation purposes
    out_file = 'global_pour_points.txt'
    print 'saving to outfile: ', out_file
    with file(out_file,'w') as outfile:
        outfile.write('OID,longitude,latitude,basin_area,min_lon,min_lat,max_lon,max_lat\n')
        np.savetxt(outfile,(np.array([basin,x_outlet,y_outlet,max_area,min_x,min_y,max_x,max_y]).T),delimiter=',')
    return i

print 'compiling and loading libraries'
import numpy as np
from netCDF4 import Dataset

nc_in = '/usr1/jhamman/Dropbox/RASM_Joe/Wu_netCDF/Wu_routing_inputs.nc'
print 'using input file: ', nc_in

print 'reading inputs...'
Inputs = read_netcdf(nc_in,('Basin_ID','Source_Area','lon','lat'))
print 'done reading inputs...'

print 'Finding pour points'
I = Find_Points(Inputs['Source_Area'],Inputs['Basin_ID'],Inputs['lon'],Inputs['lat'])
print 'Found ', str(I+1), ' Basins'
