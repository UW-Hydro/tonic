#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import struct
from netCDF4 import Dataset
import os
import argparse
import ConfigParser


###############################################################################
def main():
    (files, coord_keys, var_keys, output, binary_mult, binary_type, paths,
     out_prefix, verbose) = process_command_line()

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

            for y, x in zip(yi, xi):
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

        for y, x, point in zip(ylist, xlist, pointlist):

            data = np.empty((d[var_keys[0]].shape[0], len(var_keys)))

            for j, key in enumerate(var_keys):
                data[:, j] = d[key][:, y, x]

            if output['Binary']:
                write_binary(data*binary_mult, point, binary_type,
                             out_prefix, paths['BinaryoutPath'], append)
            if output['ASCII']:
                write_ASCII(data, point, out_prefix, paths['ASCIIoutPath'],
                            append)


###############################################################################
def process_command_line():
    """
    Parse arguments and assign flags for further loading of variables, for
    information on input arguments, type rout.py -h
    """
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=str, help="Input Configuration File")
    args = parser.parse_args()

    (files, coord_keys, var_keys, output, binary_mult, binary_type, paths,
        out_prefix, verbose) = process_config(args.config)

    return (files, coord_keys, var_keys, output, binary_mult, binary_type,
            paths, out_prefix, verbose)


###############################################################################
def process_config(config_file):
    """
    Parse arguments and assign flags for further loading of files,
    Configuration flag must be raised on command line (-C) to configuration
    file must be provided.
    Usage:  netcdf2vic.py -C netcdf2vic.cfg
    """
    config = ConfigParser.ConfigParser()
    config.read(config_file)
    config.sections()
    var_keys = config.get('Basics', 'var_keys').split(',')
    coord_keys = config.get('Basics', 'coord_keys').split(',')
    verbose = config.getboolean('Basics', 'verbose')
    WC_Start = config.getint('Basics', 'WC_Start')
    WC_End = config.getint('Basics', 'WC_End')
    output = {'Binary': config.getboolean('Basics', 'Binary'),
              'ASCII': config.getboolean('Basics', 'ASCII')}
    paths = {'in_path': config.get('Paths', 'in_path'),
             'mask_path': config.get('Paths', 'mask')}
    File_str = config.get('Basics', 'File')

    if output['ASCII']:
        paths['ASCIIoutPath'] = config.get('Paths', 'ASCIIoutPath')

    if output['Binary']:
        binary_mult = map(int, config.get('Basics', 'Mult').split(','))
        binary_type = config.get('Basics', 'Type')
        paths['BinaryoutPath'] = config.get('Paths', 'BinaryoutPath')
    else:
        binary_mult = []
        binary_type = ""

    out_prefix = config.get('Basics', 'out_prefix')
    # Make list of files, from wild cards
    if (WC_Start and WC_End):
        nums = np.arange(WC_Start, WC_End+1)
        files = []
        for i in nums:
            files.append(File_str.replace('**', str(i)))

    return (files, coord_keys, var_keys, output, binary_mult, binary_type,
            paths, out_prefix, verbose)


###############################################################################
## Read netCDF Inputs
## Read data from input netCDF.
###############################################################################
def read_netcdf(nc_file, nc_vars=[], coords=False, verbose=False):
    """
    Read data from input netCDF. Will read all variables if none provided.
    """
    f = Dataset(nc_file, 'r')
    if nc_vars == []:
        nc_vars = f.variables.keys()
    if verbose:
        print('Reading input data nc_vars: '
              '{0} from file: {1}'.format(nc_vars, nc_file))
    d = {}
    if coords:
        if isinstance(nc_vars, str):
            d[nc_vars] = np.squeeze(f.variables[nc_vars][coords])
        else:
            for var in nc_vars:
                d[var] = np.squeeze(f.variables[var][coords])
    else:
        if isinstance(nc_vars, str):
            d[nc_vars] = np.squeeze(f.variables[nc_vars][:])
        else:
            for var in nc_vars:
                d[var] = np.squeeze(f.variables[var][:])
    f.close()
    return d


###############################################################################
def write_ASCII(array, point, out_prefix, path, append, verbose=False):
    """
    Write an array to standard VIC ASCII output.
    """
    fname = out_prefix+('%1.3f' % point[0])+'_'+('%1.3f' % point[1])
    out_file = os.path.join(path, fname)
    if append:
        f = open(out_file, 'a')
    else:
        f = open(out_file, 'w')

    if verbose:
        print('Writing ASCII Data to'.format(out_file))
    np.savetxt(f, array, fmt='%12.7g')
    f.close()


###############################################################################
def write_binary(array, point, binary_type, out_prefix, path, append,
                 verbose=False):
    """
    Write a given array to standard binary short int format.
    """
    fname = out_prefix+('%.3f' % point[0])+'_'+('%.3f' % point[1])
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

if __name__ == "__main__":
    main()
