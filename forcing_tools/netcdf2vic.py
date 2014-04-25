#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import struct
import glob
from netCDF4 import Dataset
import os
import argparse
import ConfigParser
from collections import OrderedDict


class Point(object):
    '''
    Creates a point class for intellegently storing coordinate information and
    writing to disk
    '''

    def __init__(self, lat='', lon='', x='', y='', prefix='', griddecimal=4,
                 outpath='./'):
        '''Defines x and y variables'''
        self.lat = lat
        self.lon = lon
        self.x = x
        self.y = y
        self.filename = os.path.join(outpath,
                                     "{0}_{1:.{prec}f}_"
                                     "{2:.{prec}f}".format(prefix, lon, lat,
                                                           prec=griddecimal))

    def _open_binary(self):
        print('opening binary file: {0}'.format(self.filename))
        self.f = open(self.filename, 'wb')

    def _open_ascii(self):
        print('opening ascii file: {0}'.format(self.filename))
        self.f = open(self.filename, 'w')

    def _append_ascii(self, data):
        np.savetxt(self.f, data, fmt=self.float_format)
        return

    def _append_binary(self, data):
        for row in data:
            d = struct.pack(self.binary_type, *row)
            self.f.write(d)
        return

    def close(self):
        print('closing file: {0}'.format(self.filename))
        try:
            self.f.close()
        except:
            pass

    def __str__(self):
        return "Point({0},{1},{2},{3})".format(self.lat, self.lon,
                                               self.y, self.x)

    def __repr__(self):
        return "Point(lat={0}, lon={1}, \
                      y={2}, x={3}, \
                      filename={4})".format(self.lat,
                                            self.lon,
                                            self.y,
                                            self.x,
                                            self.filename)

# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def main():
    config_dict = process_command_line()

    points, ncfiles = init(config_dict)

    points = run(config_dict, ncfiles, points)

    final(config_dict, ncfiles, points)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def init(config_dict):

    # ---------------------------------------------------------------- #
    # Read domain file
    domain_dict = config_dict['DOMAIN']
    dom_data, dom_vatts, dom_gatts = read_domain(config_dict['DOMAIN'])
    ys, xs = np.nonzero(dom_data[domain_dict['MASK_VAR']])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Make sure OUTPATH is available
    if not os.path.exists(config_dict['OPTIONS']['OUTPATH']):
        os.makedirs(config_dict['OPTIONS']['OUTPATH'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Setup Points
    points = []
    for y, x in zip(ys, xs):
        points.append(Point(y=y,
                            x=x,
                            lat=dom_data[domain_dict['LATITUDE_VAR']][y, x],
                            lon=dom_data[domain_dict['LONGITUDE_VAR']][y, x],
                            prefix=config_dict['OPTIONS']['PREFIX'],
                            griddecimal=config_dict['OPTIONS']['GRIDDECIMAL'],
                            outpath=config_dict['OPTIONS']['OUTPATH']))
    points = points[:100]
    # Add fileformat specific attributes and assign open and write functions
    if config_dict['OPTIONS']['FILEFORMAT'].lower() == 'ascii':
        float_format = "%1.{0}f".format(config_dict['OPTIONS']['FLOAT_PREC'])
        for p in points:
            p._open_ascii()
            p.append = p._append_ascii
            p.float_format = float_format
    elif config_dict['OPTIONS']['FILEFORMAT'].lower() == 'binary':
        raise ValueError('Binary output is not yet complete, use ascii')
        for p in points:
            p._open_binary()
            p.append = p._append_binary
    else:
        raise ValueError('Unsupported filetype: '
                         '{0}'.format(config_dict['OPTIONS']['FILEFORMAT']))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Get a list of netCDF files
    file_dict = config_dict['NCFILES']
    if ('WC_START' in file_dict) and ('WC_END' in file_dict):
        ncfiles = []
        nums = np.arange(file_dict['WC_START'], file_dict['WC_END']+1)
        for i in nums:
            ncfiles.append(os.path.join(file_dict['PATH'],
                                        file_dict['FILENAME'].replace('**',
                                                                      str(i))))
    else:
        ncfiles = glob.glob(os.path.join(file_dict['PATH'],
                                         file_dict['FILENAME']))
    # ---------------------------------------------------------------- #


    return points, ncfiles
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def run(config_dict, ncfiles, points):

    options = config_dict['OPTIONS']
    varlist = list(options['VARIABLES'])
    nvars = len(varlist)
    file_dict = config_dict['NCFILES']

    # ---------------------------------------------------------------- #
    # Loop over ncfiles
    for ncfile in ncfiles:
        print('Reading from netCDF file: {0}'.format(ncfile))
        f = Dataset(ncfile)
        grid_data = {}
        tsize = len(f.dimensions[file_dict['TIME_DIM']])
        for v in varlist:
            grid_data[v] = f.variables[v][:]
        f.close()

        # Loop over points
        point_data = np.empty((tsize, nvars))
        for p in points:
            for i, v in enumerate(varlist):
                point_data[:, i] = grid_data[v][:, p.y, p.x]
            p.append(point_data)
    # ---------------------------------------------------------------- #

    return points
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def final(config_dict, ncfiles, points):

    # ---------------------------------------------------------------- #
    # close point files
    for p in points:
        p.close()
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # report what was done
    print('Finished netcdf2vic.py')
    print('-->Wrote to {0} {1} '
          'files'.format(len(points), config_dict['OPTIONS']['FILEFORMAT']))
    print('-->Read from {0} netCDF files'.format(len(ncfiles)))
    print('-->Output directory is here: '
          '{0}'.format(config_dict['OPTIONS']['OUTPATH']))
    # ---------------------------------------------------------------- #

    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def process_command_line():
    """
    Parse arguments and assign flags for further loading of variables, for
    information on input arguments, type netcdf2vic.py -h
    """
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=str, help="Input Configuration File")
    args = parser.parse_args()

    config_dict = read_config(args.config)

    return config_dict


# -------------------------------------------------------------------- #
# Read the Configuration File
def read_config(config_file):
    """
    Return a dictionary with subdictionaries of all configFile options/values
    """
    config = ConfigParser.SafeConfigParser()
    config.optionxform = str
    config.read(config_file)
    sections = config.sections()
    dict1 = OrderedDict()
    for section in sections:
        options = config.options(section)
        dict2 = OrderedDict()
        for option in options:
            dict2[option] = config_type(config.get(section, option))
        dict1[section] = dict2
    return dict1
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Find the type of the config options
def config_type(value):
    """
    Parse the type of the configuration file option.
    First see the value is a bool, then try float, finally return a string.
    """
    val_list = [x.strip() for x in value.split(',')]
    if len(val_list) == 1:
        value = val_list[0]
        if value in ['true', 'True', 'TRUE', 'T']:
            return True
        elif value in ['false', 'False', 'FALSE', 'F']:
            return False
        elif value in ['none', 'None', 'NONE', '']:
            return None
        elif isint(value):
            return int(value)
        elif isfloat(value):
            return float(value)
        else:
            return os.path.expandvars(value)
    else:
        try:
            return map(float, val_list)
        except:
            pass
        try:
            return map(int, val_list)
        except:
            return val_list
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def isfloat(x):
    """Test of value is a float"""
    try:
        a = float(x)
    except ValueError:
        return False
    else:
        return True
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def isint(x):
    """Test if value is an integer"""
    try:
        a = float(x)
        b = int(a)
    except ValueError:
        return False
    else:
        return a == b
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# A helper function to read a netcdf file
def read_netcdf(nc_file, variables=None, coords=None):
    """
    Read data from input netCDF. Will read all variables if none provided.
    Will also return all variable attributes.
    Both variables (data and attributes) are returned as dictionaries named
    by variable
    """

    f = Dataset(nc_file, 'r')

    if not variables:
        variables = f.variables.keys()
    if not coords:
        coords = slice(None)

    d = {}
    a = {}
    g = {}

    for var in variables:
        d[var] = f.variables[var][coords]
        a[var] = f.variables[var].__dict__

    for attr in f.ncattrs():
        g[attr] = getattr(f, attr)

    f.close()

    return d, a, g
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Read the domain
def read_domain(domain_dict):
    """
    Read the domain file and return all the variables and attributes.
    """
    dom_data, dom_vatts, dom_gatts = read_netcdf(domain_dict['FILE_NAME'])

    dom_lat = domain_dict['LATITUDE_VAR']
    dom_lon = domain_dict['LONGITUDE_VAR']

    # ---------------------------------------------------------------- #
    # Make sure longitude coords are as we expect
    posinds = np.nonzero(dom_data[dom_lon] > 180)
    dom_data[dom_lon][posinds] -= 360
    if len(posinds) > 0:
        print('adjusted lon minimum')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Make sure the longitude / latitude vars are 2d
    dom_data['cord_lons'] = dom_data[dom_lon][:]
    dom_data['cord_lats'] = dom_data[dom_lat][:]

    if dom_data[dom_lon].ndim == 1:
        # ------------------------------------------------------------- #
        # Make 2d coordinate vars
        dom_data[dom_lon], dom_data[dom_lat] = np.meshgrid(dom_data[dom_lon],
                                                           dom_data[dom_lat])
        # ------------------------------------------------------------- #
    # ---------------------------------------------------------------- #

    return dom_data, dom_vatts, dom_gatts
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
if __name__ == "__main__":
    main()
# -------------------------------------------------------------------- #
