#!/usr/local/bin/python

import sys, os
import numpy as np
from netCDF4 import Dataset
from scipy.spatial import cKDTree
import time as tm
from mpl_toolkits.basemap import Basemap, cm
import matplotlib
import matplotlib.pyplot as plt
import params
import argparse

def main():
    """
    If called from the command line, main() will process command line arguments and run make_grid,
    Make grid reads the soil file first (mandatory) and then includes any other paramter files
    Parameter column locations and units/descriptions are from params.py.
    Output is to a netcdf containing all the parameterfiles (default name is params.nc)
    """
    gridFile,soilFile,snowFile,vegFile,veglFile,outFile = process_command_line()
    grids = make_grid(gridFile, soilFile,snowFile = snowFile,vegFile = vegFile,veglFile = veglFile,ncFile = outFile)

    print 'completed grid_parms.main(), output file was:', outFile

def process_command_line():
    """
    Process command line arguments. Must have target grid (in netcdf format) and soil file (in standard vic format)
    Optionally may include snow and vegitation parameter files.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-g","--gridFile", type=str, help="Input netCDF target grid")
    parser.add_argument("-s","--soilFile", type=str, help="Input file containing soil parameter data in standard VIC format")
    parser.add_argument("-e","--snowFile", type=str, help="Input file containing snowband/elevation band data in standard VIC format",default = False)
    parser.add_argument("-v","--vegFile", type=str, help="Input file containing vegitation parameter data in standard VIC format",default = False)
    parser.add_argument("-l","--veglFile", type=str, help="Input file containing vegitation library data in standard VIC format",default = False)
    parser.add_argument("-o","--outFile", type=str, help="Input file containing vegitation parameter data in standard VIC format",default = 'params.nc')
    args = parser.parse_args()
    
    gridFile = args.gridFile
    soilFile = args.soilFile
    snowFile = args.snowFile
    vegFile = args.vegFile
    veglFile = args.veglFile
    outFile = args.outFile

    return gridFile,soilFile,snowFile,vegFile,veglFile,outFile

def make_grid(gridFile, soilFile,snowFile = False,vegFile = False,veglFile = False,ncFile = 'params.nc'):
    """
    Make grid uses routines from params.py to read standard vic format parameter files.
    After the parameter files are read, the files are placed onto the target grid using
    nearest neighbor mapping.  If a land mask is present in the target grid it will be used
    to exclude areas in the ocean.  Finally, if the ncFile = 'any_string.nc', a netcdf file
    be written with the parameter data, if ncFile = False, the dictionary of grids is returned.
    """
    targetGrid,targetAttrs = read_netcdf(gridFile)

    soilDict = params.soil(soilFile)

    if snowFile:
        snowDict = params.snow(snowFile,soilDict)

    if vegFile:
        vegDict = params.veg(vegFile,soilDict,LAIndex = True)

    if veglFile:
        veglDict = params.veg_class(veglFile)

    try:
        mask = np.ma.masked_values(targetGrid['mask'],0).mask
    except:
        mask = False

    gridDict = grid_params(soilDict,targetGrid,snowDict = snowDict, vegDict = vegDict, mask = mask)
    
    if ncFile:
        write_netcdf_main(ncFile)
        write_netcdf_vars(ncFile,targetAttrs,soilGrid = gridDict['soilDict'],snowGrid = gridDict['snowDict'],
                          veglDict = veglDict, vegGrid = gridDict['vegDict'])
    
    return gridDict

def grid_params(soilDict,targetGrid,snowDict = False, vegDict=False, fill_value = -9999,mask = False):
    """
    Reads the coordinate information from the soilDict and targetGrid and maps all input dictionaries
    to the target grid.  Returns a gridDict with the mapped input dictionary data.  
    """

    if (soilDict['lon'].min()<0 and targetGrid['xc'].min()>=0):
        posinds = np.nonzero(targetGrid['xc']>180)
        targetGrid['xc'][posinds] -= 360
        print 'adjusted xc lon minimum'
    try:
        combined = np.dstack(([targetGrid['yc'].ravel(),targetGrid['xc'].ravel()]))[0]
    except:
        targetGrid['xc'],targetGrid['yc'] = np.meshgrid(targetGrid['lon'],targetGrid['lat'])
        combined = np.dstack(([yc.ravel(),xc.ravel()]))[0]
    points = list(np.vstack((soilDict['lat'],soilDict['lon'])).transpose())
   
    mytree = cKDTree(combined)
    dist, indexes = mytree.query(points,k=1)

    inds = np.unravel_index(indexes,(targetGrid['xc'].shape[0],targetGrid['xc'].shape[1]))

    inDicts = {'soilDict':soilDict}
    outDicts = {}
    if snowDict:
        inDicts['snowDict'] = snowDict
    if vegDict:
        inDicts['vegDict'] = vegDict
    for i,name in enumerate(inDicts):

        if name == 'soilDict':
            outDict = targetGrid
        else:
            outDict = {}
            
        Dict = inDicts[name]
        for var in Dict:
            if Dict[var].ndim == 1:
                outDict[var] = np.ma.zeros((targetGrid['xc'].shape[0],targetGrid['xc'].shape[1]))+fill_value
                outDict[var][inds[0],inds[1]] = Dict[var][indexes]
                
                if isinstance(mask,np.ndarray):
                    outDict[var] = np.ma.masked_array(outDict[var],mask = mask)
                else:
                    outDict[var] = np.ma.masked_values(outDict[var],fill_value)
                
       
            elif Dict[var].ndim == 2:
                steps = Dict[var].shape[1]
                outDict[var] = np.ma.zeros((steps,targetGrid['xc'].shape[0],targetGrid['xc'].shape[1]))+fill_value
                for i in xrange(steps):
                    outDict[var][i,inds[0],inds[1]] = Dict[var][indexes,i]
                    
                    if isinstance(mask,np.ndarray):
                        outDict[var][i,:,:] = np.ma.masked_array(outDict[var][i,:,:],mask = mask)
                    else:
                        outDict[var][i,:,:] = np.ma.masked_values(outDict[var][i,:,:],fill_value)

            elif Dict[var].ndim == 3:
                y = Dict[var].shape[1]
                x = Dict[var].shape[2]
                outDict[var] = np.ma.zeros((y,x,targetGrid['xc'].shape[0],targetGrid['xc'].shape[1]))+fill_value
                for yy in xrange(y):
                    for xx in xrange(x):
                        outDict[var][yy,xx,inds[0],inds[1]] = Dict[var][indexes,yy,xx]
                    
                        if isinstance(mask,np.ndarray):
                            outDict[var][yy,xx,:,:] = np.ma.masked_array(outDict[var][yy,xx,:,:],mask = mask)
                        else:
                            outDict[var][yy,xx,:,:] = np.ma.masked_values(outDict[var][yy,xx,:,:],fill_value)

        outDicts[name] = outDict

    return outDicts

def plot_test(data,yc,xc,lats,lons):
    """
    Plot grid test to make sure the target data is landing in the same place as it should.
    data, yc, and xc are MxN grids
    lats and lons are length l vecors where l= M*N
    """
    Projection_Parameters = {'projection': 'npstere', 'boundinglat': 15,\
                             'lon_0': -114, 'lat_ts': 80.5}

    fig = plt.figure()

    m = mm = Basemap(**Projection_Parameters)
    m.drawcoastlines()
    m.drawparallels(np.arange(-80, 81, 20))
    m.drawmeridians(np.arange(-180, 181, 20))
    
    xi,yi = m(xc,yc)
    cm = matplotlib.cm.get_cmap('jet')
    m.drawlsmask(land_color = 'white',ocean_color='0.8')
    cs = m.pcolor(xi,yi,data,cmap=cm)
    cbar = m.colorbar(cs,location='right')
    cbar.set_label('GridCell')

    xx,yy = mm(lons,lats)
    xx.shape
    xx
    mm.plot(xx,yy,'go',markersize=0.5)

    plt.title('RASM Grid Cell Number')

    plt.show()

def plot_grid(data,yc,xc):
    """
    Plot grid test to make sure the target data is landing in the same place as it should.
    data, yc, and xc are MxN grids
    """
    Projection_Parameters = {'projection': 'npstere', 'boundinglat': 50,
                             'lon_0': -114, 'lat_ts': 80.5}

    fig = plt.figure()

    m = Basemap(**Projection_Parameters)
    m.drawcoastlines()
    m.drawparallels(np.arange(-80, 81, 20))
    m.drawmeridians(np.arange(-180, 181, 20))
    
    xi,yi = m(xc,yc)
    cm = matplotlib.cm.get_cmap('jet')
    m.drawlsmask(land_color = 'white',ocean_color='0.8')
    cs = m.pcolor(xi,yi,data,cmap=cm)
    cbar = m.colorbar(cs,location='right')
    cbar.set_label('fraction')

    plt.title('RASM Bare Soil')

    plt.show()
    return

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

##################################################################################
##  Write output to netCDF
##  Writes out a netCDF4 
##################################################################################
def write_netcdf_main(file):
    """
    Write output to netCDF.  Writes out a netCDF4 data file.
    Just establish the file here and write the rootgrp attributes
    """
    rootgrp = Dataset(file,'w', format='NETCDF4')

    # write attributes for netcdf
    rootgrp.description = 'VIC parameter files'
    rootgrp.history = 'Created: {}\n'.format(tm.ctime(tm.time()))
    rootgrp.history += ' '.join(sys.argv) + '\n'
    rootgrp.source = sys.argv[0] # prints the name of script used

    # write data to variables initialized above
    rootgrp.close()

def write_netcdf_vars(file,targetAttrs,soilGrid = False, snowGrid = False, vegGrid = False,veglDict = False):
    """
    Write the gridded parameters to a netcdf4 file
    Will only write paramters that it is given
    Reads attributes from params.py and from targetAtters dictionary read from gridFile
    """
    rootgrp = Dataset(file,'a', format='NETCDF4')
    
    unit = params.units(GLOBAL_LAI=True)
    desc = params.description()
        
    # set dimensions
    if soilGrid['lon'].ndim == 1:
        x = rootgrp.createDimension('lon',len(soilGrid['lon']))
        y = rootgrp.createDimension('lat',len(soilGrid['lat']))
        xcord,ycord = 'lon','lat'
    else:
        x = rootgrp.createDimension('x',soilGrid['xc'].shape[1])
        y = rootgrp.createDimension('y',soilGrid['xc'].shape[0])
        nv4 = rootgrp.createDimension('nv4',4)
        xcord,ycord = 'xc','yc'
    layers = rootgrp.createDimension('Nlayer',None)

    for var in soilGrid:

        try:
            fv = soilGrid[var].fill_value
        except:
            fv = None

        if soilGrid[var].ndim == 1:
            v = rootgrp.createVariable(var,'f8',('y','x',))
            v[:] = soilGrid[var]

        elif soilGrid[var].ndim == 2:
            v = rootgrp.createVariable(var,'f8',(ycord,xcord,),fill_value=fv)
            v[:,:] = soilGrid[var]

        elif soilGrid[var].ndim == 3:
            if (var == 'xv' or var=='yv'):
                v = rootgrp.createVariable(var,'f8',('nv4',ycord,xcord,),fill_value=fv)
            else:
                v = rootgrp.createVariable(var,'f8',('Nlayer',ycord,xcord,),fill_value=fv)
            v[:,:,:] = soilGrid[var]

        else:
            raise IOError('all vars should be 2 or 3 dimensions')


        try:
            v.units = unit.soilParam[var]
            v.description = desc.soilParam[var]
            v.coordinates = 'xc yc'

        except:
            for attr in targetAttrs[var]:
                try:
                    setattr(v,attr,targetAttrs[var][attr])
                except:
                    print 'dont have units or description for',var

    if snowGrid:
        snowbnds = rootgrp.createDimension('Snow_band',None)

        for var in snowGrid:
            try:
                fv = snowGrid[var].fill_value
            except:
                fv = None
            if snowGrid[var].ndim == 2:
                v = rootgrp.createVariable(var,'f8',(ycord,xcord,),fill_value=fv)
                v[:,:] = snowGrid[var]
            elif snowGrid[var].ndim == 3:
                v = rootgrp.createVariable(var,'f8',('Snow_band',ycord,xcord,),fill_value=fv)
                v[:,:,:] = snowGrid[var]
            else:
                raise IOError('all vars should be 2 or 3 dimensions')

            v.units = unit.snowParam[var]
            v.description = desc.snowParam[var]
            v.coordinates = 'xc yc'

    if vegGrid:
        veg_type = rootgrp.createDimension('Veg_class',None)
        root_zone = rootgrp.createDimension('Root_zone',None)
        month = rootgrp.createDimension('Month',12)

        for var in vegGrid:
            try:
                fv = vegGrid[var].fill_value
            except:
                fv = None
            if vegGrid[var].ndim == 2:
                v = rootgrp.createVariable(var,'f8',(ycord,xcord,),fill_value=fv)
                v[:,:] = vegGrid[var]
            elif vegGrid[var].ndim == 3:
                v = rootgrp.createVariable(var,'f8',('Veg_class',ycord,xcord,),fill_value=fv)
                v[:,:,:] = vegGrid[var]
            elif var == 'LAI':
                v = rootgrp.createVariable(var,'f8',('Veg_class','Month',ycord,xcord,),fill_value=fv)
                v[:,:,:,:] = vegGrid[var]
            elif vegGrid[var].ndim == 4:
                v = rootgrp.createVariable(var,'f8',('Veg_class','Root_zone','y','x',),fill_value=fv)
                v[:,:,:,:] = vegGrid[var]
            else:
                raise IOError('only able to handle dimensions <=4')

            v.units = unit.vegParam[var]
            v.description = desc.vegParam[var]
            v.coordinates = 'xc yc'

        if veglDict:
            for var in veglDict:
                if veglDict[var].ndim == 1:
                    v = rootgrp.createVariable(var,'f8',('Veg_class',))
                    v[:] = veglDict[var]
                elif veglDict[var].ndim == 2:
                    v = rootgrp.createVariable(var,'f8',('Veg_class','Month'))
                    v[:,:] = veglDict[var]
                else:
                    raise IOError('veglDict shouldnt have data with more that 2 dimentions')

                v.units = unit.vegLib[var]
                v.description = desc.vegLib[var]

    rootgrp.close()
    return

if __name__ == "__main__":
    main()

