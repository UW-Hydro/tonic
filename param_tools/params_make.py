#!/usr/local/bin/python

import sys, os
import numpy as np
from netCDF4 import Dataset
import params
import grid_params
import argparse


# snowOut = 'RASM_snow_param'
# vegOut = 'RASM_veg_param'
paramNC = 'params.nc'

# # Full Grid
# UL = False
# LR = False
# soilOut = 'RASM_soil_param'
# outfiles = 20

# Greenland 
# UL = (108,170)
# LR = (76,224)
# soilOut = 'Greenland_soil_param'
# outfiles = 1

# Boreal 1 - Canada
# UL = (37,161)
# LR = (19,177)
# soilOut = 'Boreal_1_soil_param'
# outfiles = 1

# Boreal 2 - Russia
# UL = (168,100)
# LR = (155,124)
# soilOut = 'Boreal_2_soil_param'
# outfiles = 1

# BC Coastal Range
UL = (27,134)
LR = (11,157)
soilOut = 'BCcoastal_soil_param'
outfiles = 1

def main():

    subset(paramNC,UL = UL,LR=LR,outfiles = outfiles,soilOut = soilOut,snowOut = False,vegOut = False)
    
# def full(paramNC,outfiles = 1,soilOut = False,snowOut = False,vegOut = False):

#     data,attributes = read_netcdf('params.nc')

#     cells,yinds,xinds = find_gridcells(data['mask'])

#     filesize = np.ceil(cells/outfiles)

#     for i in xrange(outfiles):
#         start = i*filesize
#         end = i*filesize+filesize
#         if end>cells:
#             end = cells
#         if outfiles>1:
#             outFile = soilOut+'_'+str(i)
#         else:
#             outFile = soilOut
#         soil(data,xinds[start:end],yinds[start:end],outFile)

def subset(paramNC,UL=False,LR=False,outfiles = 1,soilOut = False,snowOut = False,vegOut = False):

    data,attributes = read_netcdf('params.nc')

    cells,yinds,xinds = find_gridcells(data['mask'])

    if (UL and LR):
        inds = (yinds<UL[0]) & (yinds>LR[0]) & (xinds<LR[1]) & (xinds>UL[1])
        yinds = yinds[inds]
        xinds = xinds[inds]

    filesize = np.ceil(cells/outfiles)

    for i in xrange(outfiles):
        start = i*filesize
        end = i*filesize+filesize
        if end>cells:
            end = cells
        if outfiles>1:
            outFile = soilOut+'_'+str(i)
        else:
            outFile = soilOut
        soil(data,xinds[start:end],yinds[start:end],outFile)
        print 'finished', outFile
    
def soil(data,xinds,yinds,soilOut):

    c = params.cols(Nlayers=3)
    arraysize = (len(xinds),1+np.max([np.max(c.soilParam[var]) for var in c.soilParam]))
    soilParams = np.empty(arraysize)
    for var in c.soilParam:
        if data[var].ndim == 2:
            soilParams[:,c.soilParam[var]] = np.atleast_2d(data[var][yinds,xinds]).transpose()
        elif data[var].ndim == 3:
            soilParams[:,c.soilParam[var]] = np.atleast_2d(data[var][:,yinds,xinds]).transpose()
    
    np.savetxt(soilOut,soilParams,fmt='%1.6f')

# def snow():

# def veg():

# def vegLib():

def find_gridcells(mask):

    cells = np.sum(mask)

    xinds,yinds = np.nonzero(mask>0)
    
    return cells,xinds,yinds
    

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

    

if __name__ == "__main__":
    main()
