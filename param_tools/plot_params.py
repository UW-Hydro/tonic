#!/usr/local/bin/python

import numpy as np
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, cm
import matplotlib
from matplotlib import colors
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import grid_params

paramNC = '/usr1/jhamman/RASM/VIC/soil_param/params.nc'
iceNC = '/raid/jhamman/RASM_masks/RASM_landice_mask.nc'

def main():
    Pdata,Pattributes = read_netcdf(paramNC)
    Idata,Iattributes = read_netcdf(iceNC)

    baresoil = 1-np.sum(Pdata['Cv'],axis=0)
    ice = Idata['mask']/100.
    # compare_ice(Pdata['yc'],Pdata['xc'],baresoil,ice)
    plot_veg_types(Pdata['yc'],Pdata['xc'],Pdata['Cv'],baresoil)

def plot_veg_types(yc,xc,Cv,baresoil):
    Projection_Parameters = {'projection': 'npstere', 'boundinglat': 49,\
                             'lon_0': -114, 'lat_ts': 80.5}

    labels = ['Evergreen Needleleaf','Evergreen Broadleaf','Deciduous Needleleaf','Deciduous Broadleaf','Mixed Cover',
              'Woodland','Wooded Grasslands','Closed Shrublands','Open Shrublands','Grasslands','Crop land', 'Bare Soil/Ice']
    
    fig = plt.figure(figsize=(10,12))
    
    gs1 = gridspec.GridSpec(4,3)

    for loc in xrange(11):
        ax = fig.add_subplot(gs1[loc])
        c = plot_map(ax,yc,xc,Cv[loc],Projection_Parameters,vmin=0,cmap ='Jet')
        ax.set_title(labels[loc])
    ax = fig.add_subplot(gs1[11])
    c = plot_map(ax,yc,xc,baresoil,Projection_Parameters,vmin=0,cmap ='Jet')
    ax.set_title(labels[11])

    sm = plt.cm.ScalarMappable(cmap='Jet', norm=plt.normalize(vmin=0, vmax=1))
    colorbar_ax = fig.add_axes([0.92, 0.1, 0.03, 0.8])
    sm._A = []
    plt.colorbar(sm, cax=colorbar_ax)
    fig.suptitle('Fraction of Vegetation Type', fontsize=20, fontweight='bold')
    fig.text(0.5, 0.93, 'Regional Arctic Climate Model' , ha='center',fontsize=18)

    plt.show()

def compare_ice(yc,xc,baresoil,ice):
    
    ice = np.ma.masked_array(ice,mask = baresoil.mask)
    diff = ice-baresoil

    Projection_Parameters = {'projection': 'npstere', 'boundinglat': 53,\
                             'lon_0': -114, 'lat_ts': 80.5}

    fig = plt.figure()
    
    gs1 = gridspec.GridSpec(2, 2)
    gs2 = gridspec.GridSpec(2,1)

    ax1 =fig.add_subplot(gs1[0:2,0:2])
    c = plot_map(ax1,yc,xc,diff,Projection_Parameters,cbar_loc = 'right')
    ax1.set_title("Ice Fraction - Bare Soil Fraction")

    ax2 = fig.add_subplot(gs2[0])
    c =plot_map(ax2,yc,xc,ice,Projection_Parameters)
    ax2.set_title('Ice Fraction')
    
    ax3 = fig.add_subplot(gs2[1])
    c = plot_map(ax3,yc,xc,baresoil,Projection_Parameters)
    ax3.set_title('Bare Soil Fraction')

    gs1.tight_layout(fig,rect = [0.33, 0, 0.95, 1])
    gs2.tight_layout(fig,rect = [0, 0.1, 0.333, 0.9])

    plt.show()

def plot_map(ax,yc,xc,data,Projection_Parameters,cbar_loc = False,vmin=-1,vmax=1,cmap ='GrBG'):
    m = Basemap(**Projection_Parameters)
    m.drawlsmask(land_color = 'white',ocean_color='0.8')
    m.drawcoastlines()
    m.drawparallels(np.arange(-80, 81, 20))
    m.drawmeridians(np.arange(-180, 181, 20))
    xi,yi = m(xc,yc)
    cm = matplotlib.cm.get_cmap(cmap)
    m.drawlsmask(land_color = 'white',ocean_color='0.8')
    ax = m.pcolor(xi,yi,data,cmap=cm,vmin=vmin,vmax=vmax)
    if cbar_loc:
        cbar = m.colorbar(ax,location=cbar_loc)
    else:
        cbar = None

    return cbar
    
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

