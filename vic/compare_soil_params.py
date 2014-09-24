#!/usr/bin/env python
"""
compare_soil_params.py

Takes two netcdf soil parameter datasets and plots/compares the variables

"""
from __future__ import print_function
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from share import read_netcdf
from plot_utils import sub_plot_pcolor, cmap_discretize

description = 'Create plots comparing two sets of VIC soil parameters'
help = 'Create plots comparing two sets of VIC soil parameters'

rcParams.update({'figure.autolayout': True})

projection = {'urcrnrlat': 27.511827255753555,
              'urcrnrlon': 16.90845094934209,
              'llcrnrlat': 16.534986367884521,
              'llcrnrlon': 189.2229322311162,
              'projection': 'lcc',
              'rsphere': 6371200.0,
              'lon_0': -114,
              'lat_0': 90}


def _run(args, ):
    plot_atts_3 = None
    plot_atts_9 = None

    dom, dom_atts = read_netcdf(args.domain_file)
    d1, d1a = read_netcdf(args.soil_file1)
    d2, d2a = read_netcdf(args.soil_file2)
    out_path = args.out_path
    title1 = args.title1
    title2 = args.title2

    if not plot_atts_3:
        plot_atts_3 = {'infilt': {'vmin': 0,
                                  'vmax': 1,
                                  'amin': -0.5,
                                  'amax': 0.5,
                                  'amap': cmap_discretize('cm.RdBu_r')},
                       'Ws': {'vmin': 0,
                              'vmax': 100,
                              'amin': -50,
                              'amax': 50,
                              'amap': cmap_discretize('cm.RdBu_r')},
                       'Ds': {'vmin': 0,
                              'vmax': 1,
                              'amin': -0.5,
                              'amax': 0.5,
                              'amap': cmap_discretize('cm.RdBu_r')},
                       'Dsmax': {'vmin': 0,
                                 'vmax': 1,
                                 'amin': -0.5,
                                 'amax': 0.5,
                                 'amap': cmap_discretize('cm.RdBu_r')},
                       'avg_T': {'vmin': -25,
                                 'vmax': 25,
                                 'amin': -2,
                                 'amax': 2,
                                 'amap': cmap_discretize('cm.RdBu_r')},
                       'c': {'vmin': 0,
                             'vmax': 2.5,
                             'amin': -0.5,
                             'amax': 0.5,
                             'amap': cmap_discretize('cm.RdBu_r')},
                       'elev': {'vmin': 0,
                                'vmax': 2500,
                                'amin': -200,
                                'amax': 200,
                                'amap': cmap_discretize('cm.RdBu_r')},
                       'annual_prec': {'vmin': 0,
                                       'vmax': 2000,
                                       'amin': -500,
                                       'amax': 500,
                                       'amap': cmap_discretize('cm.RdBu')}}
                       # 'Nveg': {'vmin': 0, 'vmax': 10, 'amin': -2, 'amax': 2,
                            # 'amap': cmap_discretize('cm.RdBu_r')}}

    if not plot_atts_9:
        plot_atts_9 = {'soil_density': {'vmin': 0,
                                        'vmax': 4000,
                                        'amin': -500,
                                        'amax': 500,
                                        'amap': cmap_discretize('cm.RdBu_r')},
                       'bulk_density': {'vmin': 0,
                                        'vmax': 1800,
                                        'amin': -100,
                                        'amax': 100,
                                        'amap': cmap_discretize('cm.RdBu_r')},
                       'Wpwp_FRACT': {'vmin': 0,
                                      'vmax': 1,
                                      'amin': -0.4,
                                      'amax': 0.4,
                                      'amap': cmap_discretize('cm.RdBu_r')},
                       'bubble': {'vmin': 0,
                                  'vmax': 100,
                                  'amin': -20,
                                  'amax': 20,
                                  'amap': cmap_discretize('cm.RdBu_r')},
                       'quartz': {'vmin': 0,
                                  'vmax': 1,
                                  'amin': -0.25,
                                  'amax': 0.25,
                                  'amap': cmap_discretize('cm.RdBu_r')},
                       'resid_moist': {'vmin': 0,
                                       'vmax': 0.1,
                                       'amin': -0.05,
                                       'amax': 0.05,
                                       'amap': cmap_discretize('cm.RdBu')},
                       'Wcr_FRACT': {'vmin': 0,
                                     'vmax': 1,
                                     'amin': -0.5,
                                     'amax': 0.5,
                                     'amap': cmap_discretize('cm.RdBu_r')},
                       'expt': {'vmin': 0,
                                'vmax': 75,
                                'amin': -50,
                                'amax': 50,
                                'amap': cmap_discretize('cm.RdBu_r')},
                       'depth': {'vmin': 0,
                                 'vmax': 2.5,
                                 'amin': -2,
                                 'amax': 2,
                                 'amap': cmap_discretize('cm.RdBu_r')},
                       'Ksat': {'vmin': 0,
                                'vmax': 4000,
                                'amin': -1000,
                                'amax': 1000,
                                'amap': cmap_discretize('cm.RdBu_r')},
                       'init_moist': {'vmin': 0,
                                      'vmax': 200,
                                      'amin': -100,
                                      'amax': 100,
                                      'amap': cmap_discretize('cm.RdBu')}}

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
                         units=units,
                         mask=(dom['mask'] == 0),
                         t1=title1,
                         t2=title2,
                         **plot_atts_3[var])

            plt.figtext(.5, 0.94, var, fontsize=18, ha='center')

            plt.figtext(.5, 0.90, d1a[var]['description'], fontsize=12,
                        ha='center')

            fname = os.path.join(out_path,
                                 '{}-{}-{}.png'.format(title1, title2, var))
            f.savefig(fname, format='png', dpi=150, bbox_inches='tight',
                      pad_inches=0)
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
                     units=units,
                     mask=(dom['mask'] == 0),
                     t1=title1,
                     t2=title2,
                     **plot_atts_9[var])

        plt.figtext(.5, 1.06, var, fontsize=18, ha='center')
        plt.figtext(.5, 1.02, d1a[var]['description'], fontsize=12,
                    ha='center')

        fname = os.path.join(out_path,
                             '{}-{}-{}.png'.format(title1, title2, var))
        f.savefig(fname, format='png', dpi=150, bbox_inches='tight',
                  pad_inches=0)
        print('finished {}'.format(fname))

    return


def my_plot3(lons, lats, d1, d2, units=None,
             vmin=None, vmax=None, amin=None,
             amax=None, mask=True, amap=None,
             t1='', t2='', cmap=None):
    """ 3 pannel plot to compare two different datasets"""
    if vmin is None:
        vmin = d1.min()
    if vmax is None:
        vmax = d1.max()

    if not cmap:
        cmap = cmap_discretize('cm.winter')

    if not amap:
        amap = cmap_discretize('cm.RdBu')

    d1 = np.ma.masked_where(mask, d1)
    d2 = np.ma.masked_where(mask, d2)

    anom = d1-d2

    if amin is None:
        amin = -1*np.max(np.abs(anom))
    if amax is None:
        amax = np.max(np.abs(anom))

    f, axarr = plt.subplots(1, 3, figsize=(13.5, 4), dpi=150)
    f.tight_layout()

    plt.sca(axarr[0])
    sub_plot_pcolor(lons, lats, d1, vmin=vmin, vmax=vmax, ncolors=9,
                    title='{} (A)'.format(t1), cmap=cmap, units=units,
                    cbar_location='right')
    plt.sca(axarr[1])
    sub_plot_pcolor(lons, lats, d2, vmin=vmin, vmax=vmax, ncolors=9,
                    title='{} (B)'.format(t2), units=units,
                    cbar_location='right', cmap=cmap)
    plt.sca(axarr[2])
    sub_plot_pcolor(lons, lats, anom, ncolors=9, cmap=amap,
                    title='Difference (A-B)', units=units,
                    cbar_location='right', vmin=amin, vmax=amax)
    return f


# -------------------------------------------------------------------- #
def my_plot9(lons, lats, d1, d2, units=None,
             vmin=None, vmax=None, amin=None,
             amax=None, mask=True, amap=None,
             t1='', t2='', cmap=None):
    """ 9 pannel plot to compare soil layers for two different datasets"""
    if vmin is None:
        vmin = d1.min()
    if vmax is None:
        vmax = d1.max()

    anom = d1-d2
    if amin is None:
        amin = -1*np.max(np.abs(anom))
    if amax is None:
        amax = np.max(np.abs(anom))

    if cmap is None:
        cmap = cmap_discretize('cm.winter')

    if amap is None:
        amap = cmap_discretize('cm.RdBu')

    f, axarr = plt.subplots(3, 3, figsize=(13.5, 9), dpi=150)
    f.tight_layout()

    for layer in xrange(3):
        plt.sca(axarr[layer, 0])
        sub_plot_pcolor(lons, lats, np.ma.masked_where(mask, d1[layer]),
                        vmin=vmin, vmax=vmax, units=units,
                        ncolors=9, title='{} (A)'.format(t1),
                        cmap=cmap, cbar_location='right',)
        plt.ylabel('Layer {}'.format(layer))
        plt.sca(axarr[layer, 1])
        sub_plot_pcolor(lons, lats, np.ma.masked_where(mask, d2[layer]),
                        vmin=vmin, vmax=vmax,
                        ncolors=9, title='{} (B)'.format(t2), units=units,
                        cbar_location='right', cmap=cmap)
        plt.sca(axarr[layer, 2])

        sub_plot_pcolor(lons, lats, np.ma.masked_where(mask, anom[layer]),
                        ncolors=9, cmap=amap,
                        title='Difference (A-B)', units=units,
                        cbar_location='right', vmin=amin, vmax=amax)
    return f
# -------------------------------------------------------------------- #
