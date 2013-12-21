#!/usr/bin/env python
"""
compare_soil_params.py

Takes two netcdf soil parameter datasets and plots/compares the variables

"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib import cm
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import mpl_toolkits.basemap.cm as GMT
import matplotlib.colors as mcolors
import argparse
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

projection = {'urcrnrlat': 27.511827255753555,
              'urcrnrlon': 16.90845094934209,
              'llcrnrlat': 16.534986367884521,
              'llcrnrlon': 189.2229322311162,
              'projection': 'lcc',
              'rsphere': 6371200.0,
              'lon_0': -114,
              'lat_0': 90}

def main(projection=None, plot_atts_3=None, plot_atts_9=None):

    domain_file, soil_file1, soil_file2, out_path, title1, title2 = process_command_line()

    dom, dom_atts = read_netcdf(domain_file)
    d1, d1a = read_netcdf(soil_file1)
    d2, d2a = read_netcdf(soil_file2)

    if not plot_atts_3:
        plot_atts_3 = {'infilt':{'vmin':0, 'vmax':1, 'amin':-0.5, 'amax':0.5, 'amap':cmap_discretize('cm.RdBu_r')},
                       'Ws':{'vmin':0, 'vmax':100, 'amin':-50, 'amax':50, 'amap':cmap_discretize('cm.RdBu_r')},
                       'Ds':{'vmin':0, 'vmax':1, 'amin':-0.5, 'amax':0.5, 'amap':cmap_discretize('cm.RdBu_r')},
                       'Dsmax':{'vmin':0, 'vmax':1, 'amin':-0.5, 'amax':0.5, 'amap':cmap_discretize('cm.RdBu_r')},
                       'avg_T':{'vmin':-25, 'vmax':25, 'amin':-2, 'amax':2, 'amap':cmap_discretize('cm.RdBu_r')},
                       'c':{'vmin':0, 'vmax':2.5, 'amin':-0.5, 'amax':0.5, 'amap':cmap_discretize('cm.RdBu_r')},
                       'elev':{'vmin':0, 'vmax':2500, 'amin':-200, 'amax':200, 'amap':cmap_discretize('cm.RdBu_r')},
                       'annual_prec':{'vmin':0, 'vmax':2000, 'amin':-500, 'amax':500, 'amap':cmap_discretize('cm.RdBu')},}
                       # 'Nveg':{'vmin':0, 'vmax':10, 'amin':-2, 'amax':2, 'amap':cmap_discretize('cm.RdBu_r')}}

    if not plot_atts_9:
        plot_atts_9 = {'soil_density':{'vmin':0, 'vmax':4000, 'amin':-500, 'amax':500, 'amap':cmap_discretize('cm.RdBu_r')},
                       'bulk_density':{'vmin':0, 'vmax':1800, 'amin':-100, 'amax':100, 'amap':cmap_discretize('cm.RdBu_r')},
                       'Wpwp_FRACT':{'vmin':0, 'vmax':1, 'amin':-0.4, 'amax':0.4, 'amap':cmap_discretize('cm.RdBu_r')},
                       'bubble':{'vmin':0, 'vmax':100, 'amin':-20, 'amax':20, 'amap':cmap_discretize('cm.RdBu_r')},
                       'quartz':{'vmin':0, 'vmax':1, 'amin':-0.25, 'amax':0.25, 'amap':cmap_discretize('cm.RdBu_r')},
                       'resid_moist':{'vmin':0, 'vmax':0.1, 'amin':-0.05, 'amax':0.05, 'amap':cmap_discretize('cm.RdBu')},
                       'Wcr_FRACT':{'vmin':0, 'vmax':1, 'amin':-0.5, 'amax':0.5, 'amap':cmap_discretize('cm.RdBu_r')},
                       'expt':{'vmin':0, 'vmax':75, 'amin':-50, 'amax':50, 'amap':cmap_discretize('cm.RdBu_r')},
                       'depth':{'vmin':0, 'vmax':2.5, 'amin':-2, 'amax':2, 'amap':cmap_discretize('cm.RdBu_r')},
                       'Ksat':{'vmin':0, 'vmax':4000, 'amin':-1000, 'amax':1000, 'amap':cmap_discretize('cm.RdBu_r')},
                       'init_moist':{'vmin':0, 'vmax':200, 'amin':-100, 'amax':100, 'amap':cmap_discretize('cm.RdBu')}}

    # surface plots
    for var in plot_atts_3.keys():
        print('making plot3 for {}'.format(var))
        try:
            units = d1a[var]['units']
        except:
            units = ''
        try:
            f = my_plot3(dom['xc'],
                         dom['yc'],
                         d1[var],
                         d2[var],
                         units = units,
                         mask=(dom['mask']==0),
                         t1=title1,
                         t2=title2,
                         **plot_atts_3[var])

            plt.figtext(.5,0.94, var, fontsize=18, ha='center')

            plt.figtext(.5,0.90, d1a[var]['description'] ,fontsize=12,ha='center')

            fname = os.path.join(out_path, '{}-{}-{}.png'.format(title1, title2, var))
            f.savefig(fname, format='png', dpi=150, bbox_inches='tight', pad_inches=0)
            print('finished {}'.format(fname))
        except:
            print('problem with {}'.format(fname))

	# level plots
    for var in plot_atts_9.keys():
        print('making plot9 for {}'.format(var))
        try:
            units = d1a[var]['units']
        except:
            units = ''
        f = my_plot9(dom['xc'],
                     dom['yc'],
                     d1[var],
                     d2[var],
                     units = units,
                     mask=(dom['mask']==0),
                     t1=title1,
                     t2=title2,
                     **plot_atts_9[var])

        plt.figtext(.5,1.06, var, fontsize=18, ha='center')
        plt.figtext(.5,1.02, d1a[var]['description'], fontsize=12,ha='center')

        fname = os.path.join(out_path, '{}-{}-{}.png'.format(title1, title2, var))
        f.savefig(fname, format='png', dpi=150, bbox_inches='tight', pad_inches=0)
        print('finished {}'.format(fname))

    return

def my_plot3(lons, lats, d1, d2, units=None,
             vmin=None, vmax=None, amin=None,
             amax=None, mask=True, amap=None,
             t1='', t2='', cmap=None):
    """ 3 pannel plot to compare two different datasets"""
    if not vmin:
        vmin=d1.min()
    if not vmax:
        vmax=d1.max()

    if not cmap:
        cmap = cmap_discretize('cm.winter')

    if not amap:
        amap = cmap_discretize('cm.RdBu')

    d1 = np.ma.masked_where(mask, d1)
    d2 = np.ma.masked_where(mask, d2)

    anom=d1-d2

    if not amin:
        amin = -1*np.max(np.abs(anom))
    if not amin:
        amax = np.max(np.abs(anom))

    f, axarr = plt.subplots(1, 3, figsize=(13.5, 4), dpi=150)
    f.tight_layout()

    plt.sca(axarr[0])
    sub_plot_pcolor(lons, lats, d1, vmin=vmin, vmax=vmax, ncolors=9,
                    title='{} (A)'.format(t1), cbar=False, cmap=cmap)
    plt.sca(axarr[1])
    sub_plot_pcolor(lons, lats, d2, vmin=vmin, vmax=vmax, ncolors=9,
                    title='{} (B)'.format(t2), units='', cbar_location='right', cmap=cmap)
    plt.sca(axarr[2])
    sub_plot_pcolor(lons, lats, anom, ncolors=9, cmap=amap,
                    title='Difference (A-B)', units=units, cbar_location='right', vmin=amin, vmax=amax)
    return f

def my_plot9(lons, lats, d1, d2, units=None,
             vmin=None, vmax=None, amin=None,
             amax=None, mask=True, amap=None,
             t1='', t2='', cmap=None):
    """ 9 pannel plot to compare soil layers for two different datasets"""
    if not vmin:
        vmin=d1.min()
    if not vmax:
        vmax=d1.max()

    anom=d1-d2
    if not amin:
        amin = -1*np.max(np.abs(anom))
    if not amin:
        amax = np.max(np.abs(anom))

    if not cmap:
        cmap = cmap_discretize('cm.winter')

    if not amap:
        amap = cmap_discretize('cm.RdBu')

    f, axarr = plt.subplots(3, 3, figsize=(13.5, 9), dpi=150)
    f.tight_layout()

    for layer in xrange(3):
        plt.sca(axarr[layer, 0])
        sub_plot_pcolor(lons, lats, np.ma.masked_where(mask, d1[layer]), vmin=vmin, vmax=vmax,
                        ncolors=9, title='{} (A)'.format(t1), cbar=False, cmap=cmap)
        plt.ylabel('Layer {}'.format(layer))
        plt.sca(axarr[layer, 1])
        sub_plot_pcolor(lons, lats, np.ma.masked_where(mask, d2[layer]), vmin=vmin, vmax=vmax,
                        ncolors=9, title='{} (B)'.format(t2), units='',
                        cbar_location='right', cmap=cmap)
        plt.sca(axarr[layer, 2])

        sub_plot_pcolor(lons, lats, np.ma.masked_where(mask, anom[layer]), ncolors=9, cmap=amap,
                        title='Difference (A-B)', units=units,
                        cbar_location='right', vmin=amin, vmax=amax)
    return f

def sub_plot_pcolor(lons, lats, data, title=None, cmap=cm.jet,
                    vmin=None, vmax=None, cbar=True, cbar_location='bottom',
                    units=None, ncolors=5):

    if vmin==None:
        vmin=data.min()
    if vmax==None:
        vmax=data.max()

    m = Basemap(**projection)
    m.drawlsmask(land_color='grey', lakes=False)
    xi,yi = m(np.squeeze(lons),np.squeeze(lats))
    sp = m.pcolormesh(xi, yi, np.squeeze(data), vmin=vmin,vmax=vmax, cmap=cmap)
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))
    m.drawcoastlines(color='k',linewidth=0.25);
    if title:
        plt.title(title,size=13)
    if cbar:
        cmap = cmap_discretize(cmap, ncolors)
        mappable = cm.ScalarMappable(cmap=cmap)
        mappable.set_array([])
        mappable.set_clim(-0.5, ncolors+0.5)
        colorbar = m.colorbar(mappable, location=cbar_location)
        colorbar.set_ticks(np.linspace(0, ncolors, ncolors))
        colorbar.set_ticklabels(np.linspace(vmin, vmax, ncolors))
        colorbar.set_label(units)
    return sp

def cmap_discretize(cmap, N=9):
    """Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet.
        N: number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    if type(cmap) == str:
        cmap = cm.get_cmap(eval(cmap))
    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki])
                       for i in xrange(N+1) ]
    # Return colormap object.
    return mcolors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)


def read_netcdf(nc_file, variables = [],coords = False, verbose = False):
    """
    Read data from input netCDF. Will read all variables if none provided.
    Will also return all variable attributes.
    Both variables (data and attributes) are returned as dictionaries named by variable
    """
    if verbose:
      print('Reading input data variables:', variables, ', from file:', nc_file)

    f = Dataset(nc_file, 'r')

    if variables==[]: variables = f.variables.keys()

    d={}
    a={}

    if coords:
        if isinstance(variables,str):
            d[variables] = f.variables[variables][coords]
            a[variables] = f.variables[variables].__dict__
        else:
            for var in variables:
                d[var] = f.variables[var][coords]
                a[var] = f.variables[var].__dict__
    else:
        if isinstance(variables,str):
            d[variables] = f.variables[variables][:]
            a[variables] = f.variables[variables].__dict__
        else:
            for var in variables:
                d[var] = f.variables[var][:]
                a[var] = f.variables[var].__dict__
    f.close()
    return d, a

def process_command_line():
    """
    Process command line arguments. Must have target grid (in netcdf format)
    and soil file (in standard vic format)
    Optionally may include snow and vegitation parameter files.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--domain_file",
                        type=str,
                        help="Input netCDF target grid",
                        required=True)

    parser.add_argument("-s1", "--soil_file1",
                        type=str,
                        help="Input file containing soil parameter",
                        required=True)
    parser.add_argument("-s2", "--soil_file2",
                        type=str,
                        help="Input file containing soil parameter",
                        required=True)
    parser.add_argument("-o", "--out_path",
                        type=str,
                        help="outpath for files",
                        default='./')
    parser.add_argument("-t1", "--title1",
                        type=str,
                        help="Input tile for soil_file1",
                        default=None)
    parser.add_argument("-t2", "--title2",
                        type=str,
                        help="Input tile for soil_file2",
                        default=None)
    args = parser.parse_args()

    domain_file = args.domain_file
    soil_file1 = args.soil_file1
    soil_file2 = args.soil_file2
    out_path = args.out_path

    if args.title1:
        title1 = args.title1
    else:
        title1 = os.path.splitext(os.path.split(soil_file1)[1])[0]

    if args.title2:
        title2 = args.title2
    else:
        title2 = os.path.splitext(os.path.split(soil_file2)[1])[0]

    return domain_file, soil_file1, soil_file2, out_path, title1, title2

if __name__ == "__main__":
    main()