#!/usr/bin/env python
"""netcdf2vic.py"""

from __future__ import print_function
import numpy as np
import struct
import os
from tonic.io import read_netcdf, read_config
from tonic.pycompat import pyzip

description = 'Convert netCDF meteorological forcings to VIC sytle format'
help = 'Convert netCDF meteorological forcings to VIC sytle format'


# -------------------------------------------------------------------- #
# top level run function
def _run(args):

    config = read_config(args.config)
    files = config['options']['files']
    var_keys = config['options']['var_keys']
    output = config['options']['output']
    binary_mult = config['options']['binary_mult']
    binary_type = config['options']['binary_type'],
    paths = config['options']['paths']
    out_prefix = config['options']['out_prefix']
    verbose = config['options']['verbose']

    mask = read_netcdf(paths['mask_path'], nc_vars=['mask'])['mask']
    yi, xi = np.nonzero(mask)
    print('found {0} points in mask file.'.format(len(yi)))

    xlist = []
    ylist = []
    pointlist = []
    append = False

    for i, fname in enumerate(files):
        d = read_netcdf(os.path.join(paths['in_path'], fname),
                        verbose=verbose)

        if i == 0:

            # find point locations
            xs = d['xc']
            ys = d['yc']
            posinds = np.nonzero(xs > 180)
            xs[posinds] -= 360
            print('adjusted xs lon minimum')

            for y, x in pyzip(yi, xi):
                active_flag = False
                for key in var_keys:
                    if (d[key][:, y, x].all() is np.ma.masked) \
                            or (mask[y, x] == 0):
                        active_flag = True
                if not active_flag:
                    point = (ys[y, x], xs[y, x])
                    xlist.append(x)
                    ylist.append(y)
                    pointlist.append(point)

        else:
            append = True

        for y, x, point in pyzip(ylist, xlist, pointlist):

            data = np.empty((d[var_keys[0]].shape[0], len(var_keys)))

            for j, key in enumerate(var_keys):
                data[:, j] = d[key][:, y, x]

            if output['Binary']:
                write_binary(data * binary_mult, point, binary_type,
                             out_prefix, paths['BinaryoutPath'], append)
            if output['ASCII']:
                write_ascii(data, point, out_prefix, paths['ASCIIoutPath'],
                            append)
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Write ASCII
def write_ascii(array, point, out_prefix, path, append, verbose=False):
    """
    Write an array to standard VIC ASCII output.
    """
    fname = out_prefix + ('%1.3f' % point[0]) + '_' + ('%1.3f' % point[1])
    out_file = os.path.join(path, fname)
    if append:
        f = open(out_file, 'a')
    else:
        f = open(out_file, 'w')

    if verbose:
        print('Writing ASCII Data to'.format(out_file))
    np.savetxt(f, array, fmt='%12.7g')
    f.close()
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Write Binary
def write_binary(array, point, binary_type, out_prefix, path, append,
                 verbose=False):
    """
    Write a given array to standard binary short int format.
    """
    fname = out_prefix + ('%.3f' % point[0]) + '_' + ('%.3f' % point[1])
    out_file = os.path.join(path, fname)
    if verbose:
        print('Writing Binary Data to'.format(out_file))
    if append:
        f = open(out_file, 'a')
    else:
        f = open(out_file, 'w')

    for row in array:
        data = struct.pack(binary_type, *row)
        f.write(data)
    f.close()
    return
# -------------------------------------------------------------------- #
