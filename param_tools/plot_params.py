#!/usr/local/bin/python

import numpy as np
from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap

projection = {'projection': 'npstere', 'boundinglat': 49,
              'lon_0': -114, 'lat_ts': 80.5}

paramNC = 'param_file'
iceNC = 'ice_file'


def main():
    Pdata, Pattributes = read_netcdf(paramNC)
    Idata, Iattributes = read_netcdf(iceNC)

    baresoil = 1-np.sum(Pdata['Cv'], axis=0)

    plot_veg_types(Pdata['yc'], Pdata['xc'], Pdata['Cv'], baresoil)


def plot_veg_types(yc, xc, Cv, baresoil):
    Projection_Parameters = projection

    labels = ['Evergreen Needleleaf', 'Evergreen Broadleaf',
              'Deciduous Needleleaf', 'Deciduous Broadleaf', 'Mixed Cover',
              'Woodland', 'Wooded Grasslands', 'Closed Shrublands',
              'Open Shrublands', 'Grasslands', 'Crop land', 'Bare Soil/Ice']

    fig = plt.figure(figsize=(10, 12))

    gs1 = gridspec.GridSpec(4, 3)

    for loc in xrange(11):
        ax = fig.add_subplot(gs1[loc])
        c = plot_map(ax, yc, xc, Cv[loc], Projection_Parameters, vmin=0,
                     cmap='Jet')
        ax.set_title(labels[loc])
    ax = fig.add_subplot(gs1[11])
    c = plot_map(ax, yc, xc, baresoil, Projection_Parameters, vmin=0,
                 cmap='Jet')
    ax.set_title(labels[11])

    sm = plt.cm.ScalarMappable(cmap='Jet', norm=plt.normalize(vmin=0, vmax=1))
    colorbar_ax = fig.add_axes([0.92, 0.1, 0.03, 0.8])
    sm._A = []
    plt.colorbar(sm, cax=colorbar_ax)
    fig.suptitle('Fraction of Vegetation Type', fontsize=20, fontweight='bold')
    fig.text(0.5, 0.93, 'Regional Arctic Climate Model',
             ha='center', fontsize=18)

    plt.show()


def compare_ice(yc, xc, baresoil, ice):

    ice = np.ma.masked_array(ice, mask=baresoil.mask)
    diff = ice-baresoil

    Projection_Parameters = projection

    fig = plt.figure()

    gs1 = gridspec.GridSpec(2, 2)
    gs2 = gridspec.GridSpec(2, 1)

    ax1 = fig.add_subplot(gs1[0:2, 0:2])
    c = plot_map(ax1, yc, xc, diff, Projection_Parameters, cbar_loc='right')
    ax1.set_title("Ice Fraction - Bare Soil Fraction")

    ax2 = fig.add_subplot(gs2[0])
    c = plot_map(ax2, yc, xc, ice, Projection_Parameters)
    ax2.set_title('Ice Fraction')

    ax3 = fig.add_subplot(gs2[1])
    c = plot_map(ax3, yc, xc, baresoil, Projection_Parameters)
    ax3.set_title('Bare Soil Fraction')

    gs1.tight_layout(fig, rect=[0.33, 0, 0.95, 1])
    gs2.tight_layout(fig, rect=[0, 0.1, 0.333, 0.9])

    plt.show()


def plot_map(ax, yc, xc, data, Projection_Parameters, cbar_loc=False, vmin=-1,
             vmax=1, cmap='GrBG'):

    m = Basemap(**Projection_Parameters)
    m.drawlsmask(land_color='white', ocean_color='0.8')
    m.drawcoastlines()
    m.drawparallels(np.arange(-80, 81, 20))
    m.drawmeridians(np.arange(-180, 181, 20))
    xi, yi = m(xc, yc)
    cm = matplotlib.cm.get_cmap(cmap)
    m.drawlsmask(land_color='white', ocean_color='0.8')
    ax = m.pcolor(xi, yi, data, cmap=cm, vmin=vmin, vmax=vmax)
    if cbar_loc:
        cbar = m.colorbar(ax, location=cbar_loc)
    else:
        cbar = None

    return cbar


###############################################################################
## Read netCDF Inputs
## Read data from input netCDF.
###############################################################################
def read_netcdf(nc_file, variables=[], coords=False, verbose=False):
    """
    Read data from input netCDF. Will read all variables if none provided.
    Will also return all variable attributes.
    Both variables (data and attributes) are returned as dictionaries named by
    variable
    """
    if verbose:
        print('Reading input data variables: {0} '
              'from file: {1}'.format(variables, nc_file))
    f = Dataset(nc_file, 'r')
    if variables == []:
        variables = f.variables.keys()
    d = {}
    a = {}
    if coords:
        if isinstance(variables, str):
            d[variables] = f.variables[variables][coords]
            a[variables] = f.variables[variables].__dict__
        else:
            for var in variables:
                d[var] = f.variables[var][coords]
                a[var] = f.variables[var].__dict__
    else:
        if isinstance(variables, str):
            d[variables] = f.variables[variables][:]
            a[variables] = f.variables[variables].__dict__
        else:
            for var in variables:
                d[var] = f.variables[var][:]
                a[var] = f.variables[var].__dict__
    f.close()
    return d, a

if __name__ == "__main__":
    main()
