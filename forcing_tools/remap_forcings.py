#!/usr/local/bin/python

import os
import glob
import tarfile
import numpy as np
from cdo import *
cdo = Cdo()

input_file_path = '/raid/jhamman/raw_forcings/Adam2006/Global/'
grid = '/raid/jhamman/RASM_masks/domain.lnd.wr50a_ar9v4.100920.nc'
outpath = '/raid/jhamman/RASM_met_forcings/Adam2006_Global/'

def standard_remap():
    print 'running standard remap program'
    #change to input directory
    os.chdir(input_file_path)
    # get list of files
    files = glob.glob("*.nc")
    for file in files:
        outfile = outpath+file[:-3]+'_RASM.nc'
        if file.startswith("Precip"):
            cdo.remapcon(grid,input = file,output = outfile)
            print 'done with', file
        elif file.startswith("shum"):
            cdo.remapcon(grid,input = file,output = outfile)
            print 'done with', file
        elif file.startswith("prcp"):
            cdo.remapcon(grid,input = file,output = outfile)
            print 'done with', file
        elif file.startswith("dswrf"):
            cdo.remapcon(grid,input = file,output = outfile)
            print 'done with', file
        elif file.startswith("dlwrf"):
            cdo.remapcon(grid,input = file,output = outfile)
            print 'done with', file
        else:
            cdo.remapbil(grid,input = file,output = outfile)
            print 'done with', file
    print 'done remapping all files in', input_file_path

def global_remap():
    print 'running remap program for Adam2006 data (global)'
    # get list of files
    for year in np.arange(1948,2007):
        dir = input_file_path+str(year)+'/'
        files = glob.glob(dir+"*")
        for file in files:
            if file.endswith('.tar.gz'):
                tfile = tarfile.open(file)
                member = tfile.getnames()
                tfile.extractall(dir)
                if member[0].endswith('.nc'):
                    infile = dir+member[0]
                    outfile = outpath+member[0][:-3]+'_RASM_BIL.nc'
                    print 'trying input:', infile, 'output:', outfile
                    cdo.remapbil(grid,input = infile,output = outfile)
                    outfile = outpath+member[0][:-3]+'_RASM_CON.nc'
                    cdo.remapcon(grid,input = infile,output = outfile)
                    print 'done with', member 

#standard_remap()
global_remap()
        
