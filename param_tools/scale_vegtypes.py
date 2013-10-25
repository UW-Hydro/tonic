#!/usr/local/bin/python

import numpy as np
from params_make import read_netcdf, veg, find_gridcells

# snowOut = 'RASM_snow_param'
vegOut = 'RASM_veg_param_nobare'
VparamNC = '/usr1/jhamman/RASM/VIC/soil_param/params_new.nc'
IparamNC = '/raid/jhamman/RASM_masks/RASM_landice_mask.nc'

def main():
    nobare(VparamNC,IparamNC,vegOut)

def nobare(VparamNC,IparamNC,vegOut):
    """
    Scale vegitation tiles to take up an space not taken up by ice
    """

    Vdata,Vattributes = read_netcdf(VparamNC)
    Idata,Iattriubtes = read_netcdf(IparamNC)

    cells,yinds,xinds = find_gridcells(Vdata['mask'])
    
    down, up = 0,0
    for i in xrange(len(yinds)):
        y,x = yinds[i], xinds[i]
        Cv = Vdata['Cv'][:,y,x]
        TotCv = np.sum(Cv)
        TotIce = Idata['mask'][y,x]/100.
        bare = 1-TotCv-TotIce
        if bare>0:
            Cv /= TotCv/(1-TotIce)
            up+=1
            print Vdata['cellnum'][y,x]
        elif bare<0:
            Cv /=TotCv/(1-TotIce)
            down+=1
        Vdata['Cv'][:,y,x] = Cv

    print 'downscaled', down, 'cells, upscaled ', up,'cells'
    
    veg(Vdata,xinds,yinds,vegOut,rootzones=2,GLOBAL_LAI=True)

if __name__ == "__main__":
    main()
