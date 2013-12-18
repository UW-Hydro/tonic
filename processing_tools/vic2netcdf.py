#!/usr/bin/env python
"""
Python implementation of vic2nc

This module facilitates the conversion of ascii VIC output files into 3 or 4
dimenstional netcdf files.

References:
 - VIC: http://www.hydro.washington.edu/Lettenmaier/Models/VIC/index.shtml
 - netCDF: http://www.unidata.ucar.edu/software/netcdf/
 - Python netCDF4: https://code.google.com/p/netcdf4-python/
 - NetCDF Climate and Forecast (CF) Metadata Convention: http://cf-pcmdi.llnl.gov/
 - Pandas: http://pandas.pydata.org/
"""
# Imports
from os import path
from glob import glob
from re import findall
from collections import OrderedDict
from bisect import bisect_right, bisect_left
from argparse import ArgumentParser
from getpass import getuser
from datetime import datetime
from pandas import read_csv
from netCDF4 import Dataset, date2num, num2date, default_fillvals
from ConfigParser import SafeConfigParser
from scipy.spatial import cKDTree
from getpass import getuser
import socket
import subprocess
import dateutil.relativedelta as relativedelta
import os, sys
import numpy as np
import time as tm

REFERENCE_STRING = '0001-1-1 0:0:0'
TIMEUNITS = 'days since ' + REFERENCE_STRING     # do not change (MUST BE DAYS)!
TIMESTAMPFORM = '%Y-%m-%d-%H'

# Precision
NC_DOUBLE = 'f8'
NC_FLOAT = 'f4'

# Default configuration
default_config = {'OPTIONS':{'out_file_format': 'NETCDF3_64BIT',
                             'precision': 'single',
                             'calendar': 'standard',
                             'time_segment': 'month',
                             'snow_bands': False,
                             'veg_tiles': False,
                             'soil_layers': False},
                  'DOMAIN':{'longitude_var': 'longitude',
                            'latitude_var': 'latitude',
                            'y_x_dims': ['y', 'x']}}

# -------------------------------------------------------------------- #
# Point object
class Point(object):
    '''Creates a point class for intellegently storing coordinate information'''

    def __init__(self, lat='', lon='', x='', y='', filename=''):
        '''Defines x and y variables'''
        self.lat = lat
        self.lon = lon
        self.x = x
        self.y = y
        self.filename = filename

    def read(self, names=[], usecols=[]):
        print('reading %s' %self.filename)
        self.df = read_csv(self.filename,
                           delimiter='\t',
                           header=None,
                           usecols=usecols)
        self.df.columns = names
        return

    def __str__(self):
        return "Point(%s,%s,%s,%s)" % (self.lat, self.lon, self.y, self.x)

    def __repr__(self):
        return "Point(lat=%s, lon=%s, y=%s, x=%s, filename=%s)" % (self.lat, self.lon, self.y, self.x, self.filename)

# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Point List object
class Plist(list):
    '''List subclass that has a few helper methods for adding and obtaining coordinates'''

    def get_lons(self):
        return np.array([p.lon for p in self])

    def get_lats(self):
        return np.array([p.lat for p in self])

    def add_xs(self, xinds):
        for i in xrange(len(self)):
            self[i].x = xinds[i]
        return

    def add_ys(self, yinds):
        for i in xrange(len(self)):
            self[i].y = yinds[i]
        return
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Segment object
class Segment(object):
    def __init__(self, num, i0, i1, nc_format, filename, big_memory=False):
        '''Class used for holding segment information '''
        self.num = num
        self.i0 = i0
        self.i1 = i1
        self.filename = filename
        self.fields = {}
        self.big_memory = big_memory

        self.nc_write(nc_format)


    def nc_globals(self,
                   title='VIC netCDF file',
                   history='Created: {} by {}'.format(tm.ctime(tm.time()), getuser()),
                   institution='University of Washington',
                   source=sys.argv[0],
                   references='Primary Historical Reference for VIC: Liang, X., D. P. Lettenmaier, E. F. Wood, and S. J. Burges, 1994: A Simple hydrologically Based Model of Land Surface Water and Energy Fluxes for GSMs, J. Geophys. Res., 99(D7), 14,415-14,428.',
                   comment='Output from the Variable Infiltration Capacity (VIC) Macroscale Hydrologic Model',
                   Conventions='CF-1.6',
                   target_grid_file='unknown',
                   username=None,
                   hostname=None,
                   version=None,
                   **kwargs):

        self.f.title = title
        self.f.history = history
        self.f.institution = institution
        self.f.source = source
        self.f.references = references
        self.f.comment = comment
        self.f.Conventions = Conventions
        if hostname:
            self.f.hostname = hostname
        else:
            self.f.hostname = socket.gethostname()
        if username:
            self.f.username = username
        else:
            self.f.username = getuser()
        if version:
            self.f.version = version
        else:
            try:
                self.f.version = subprocess.check_output(["git", "describe"]).rstrip()
            except:
                self.f.version = 'unknown'

        for attribute, value in kwargs.iteritems():
            if hasattr(self.f, attribute):
                print('WARNING: Attribute %s already exists.')
                print('Renaming to g_%s to ovoid overwriting.' %(attribute, attribute))
                attribute='g_'+attribute
            setattr(self.f, attribute, value)
        return

    def __str__(self):
        return "Segment Object(%s)" % (self.filename)

    def __repr__(self):
        return """-------------------------- Segment %s --------------------------
Filename: %s
Start Index: %s
End Index: %s
------------------------------------------------------------------""" %(self.num, self.filename, self.i0, self.i1)

    def nc_time(self, times, calendar):
        """ define time dimension (and write data) """
        time = self.f.createDimension('time', len(times[self.i0:self.i1]))
        time = self.f.createVariable('time','f8',('time',))
        time[:] = times[self.i0:self.i1]
        time.units = TIMEUNITS
        time.calendar = calendar

    def nc_domain(self, domain):
        """ define the coordinate dimension (and write data) """
        # Setup dimensions
        dimensions = []
        for name, ncvar in domain.iteritems():
            # Setup dimensions
            for dim in ncvar.dimensions:
                if dim not in dimensions:
                    dimensions.append(dim)
                    d = self.f.createDimension(dim, getattr(ncvar, dim))
            # Create variable
            if "_FillValue" in ncvar.attributes.keys():
                fill_val = ncvar.attributes['_FillValue']
                del ncvar.attributes['_FillValue']
            else:
                fill_val = None
            self.fields[name] = self.f.createVariable(name, NC_DOUBLE,
                                                      ncvar.dimensions,
                                                      fill_value=fill_val)
            # Apply the data
            self.fields[name][:] = ncvar
            # Add the attributes
            for key, val in ncvar.attributes.iteritems():
                setattr(self.fields[name], key, val)

        return

    def nc_dimensions(self, snow_bands=False, veg_tiles=False, soil_layers=False):
        """ Define 4th dimensions """
        if snow_bands:
            d = self.f.createDimension('snow_bands', snow_bands)
        if veg_tiles:
            d = self.f.createDimension('veg_tiles', veg_tiles)
        if soil_layers:
            d = self.f.createDimension('soil_layers', soil_layers)
        return

    def nc_fields(self, fields, y_x_dims, precision):
        """ define each field """
        coords = ('time',)+tuple(y_x_dims)

        self.var_list = fields

        if precision == 'single':
            prec_global = NC_FLOAT
        elif precision == 'double':
            prec_global = NC_DOUBLE
        else:
            raise ValueError('Unkown value for OPTIONS[precision] field: %s', precision)

        self.three_dim_vars = []
        self.four_dim_vars = []

        if self.big_memory:
            self.data = {}

        for name, field in fields.iteritems():

            if 'dim4' in field.keys():
                if len(field['column'])==len(self.f.dimensions[field['dim4']]):
                    # 4d var
                    coords = ('time',)+tuple([field['dim4']])+tuple(y_x_dims)
                    self.four_dim_vars.append(name)
                elif len(field['column'])!=len(self.f.dimensions[field['dim4']]):
                    raise ValueError('Number of columns for variable %s \
                                     does not match the length (%s) of the \
                                     %s dimension' %(name, len(self.f.dimensions[field['dim4']]), field['dim4']))
            else:
                # standard 3d var
                coords = ('time',)+tuple(y_x_dims)
                self.three_dim_vars.append(name)

            if 'type' in field.keys():
                prec = field['type']
            else:
                prec = prec_global
            fill_val = default_fillvals[prec]

            self.fields[name] = self.f.createVariable(name, prec, coords, fill_value=fill_val)

            if 'units' in field.keys():
                self.fields[name].long_name = name
                self.fields[name].coordinates = 'lon lat'
                for key, val in field.iteritems():
                    setattr(self.fields[name], key, val)
            else:
                raise ValueError('Field %s missing units attribute' %name)

            # setup array for big_memory
            if self.big_memory:
                self.data[name] = np.zeros(self.fields[name].shape, dtype=prec) + fill_val

        return

    def nc_add_data_big_memory(self, point):
        for name in self.three_dim_vars:
            self.data[name][:, point.y, point.x] = point.df[name].values[self.i0:self.i1]
        for name in self.four_dim_vars:
            varshape = self.f.variables[name].shape[1]
            for i in xrange(varshape):
                subname = name + str(i)
                self.data[name][:, i, point.y, point.x] = point.df[subname].values[self.i0:self.i1]

    def nc_add_data_standard(self, data):
        for filename, point in data.iteritems():
            for name in self.three_dim_vars:
                self.f.variables[name][:, point.y, point.x] = point.df[name].values[self.i0:self.i1]
            for name in self.four_dim_vars:
                varshape = self.f.variables[name].shape[1]
                for i in xrange(varshape):
                    subname = name + str(i)
                    self.f.variables[name][:, i, point.y, point.x] = point.df[subname].values[self.i0:self.i1]

    def nc_write_data_big_memory(self):
        """ write completed data arrays to disk """
        for name in self.three_dim_vars:
            self.f.variables[name][:, :, :] = self.data[name]
        for name in self.four_dim_vars:
                self.f.variables[name][:, :, :, :] = self.data[name]

    def nc_write(self, nc_format):
        self.f = Dataset(self.filename, mode="w", clobber=True, format=nc_format)
        print('Opened in write mode: %s' %self.filename)

    def nc_close(self):
        self.f.close()
        print('Closed: %s' %self.filename)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
class NcVar(np.ndarray):
    """ Subclass of numpy array to cary netcdf attributes"""
    def __new__(cls, f, varname):
        obj = np.asarray(f.variables[varname][:]).view(cls)
        # add the new attribute to the created instance
        obj.dimensions = f.variables[varname].dimensions
        obj.attributes = f.variables[varname].__dict__
        for dim in obj.dimensions:
            setattr(obj, dim, len(f.dimensions[dim]))
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
class FakeNcVar(np.ndarray):
    """ Subclass of numpy array to carry netcdf attributes"""
    def __new__(cls, data, dimensions, attributes):
        obj = np.asarray(data).view(cls)
        # add the new attribute to the created instance
        obj.dimensions = dimensions
        obj.attributes = attributes
        shape = data.shape
        for i, dim in enumerate(obj.dimensions):
            setattr(obj, dim, shape[i])
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Top level driver
def main():

    # ---------------------------------------------------------------- #
    # Read command Line
    config_file, big_memory, create_batch, batch_dir = process_command_line()
    # ---------------------------------------------------------------- #

    if create_batch:
        # ------------------------------------------------------------ #
        # Create batch files and exit
        batch(config_file, create_batch, batch_dir)
        # ------------------------------------------------------------ #
    else:
        # ------------------------------------------------------------ #
        # Read Configuration files
        config_dict = read_config(config_file)
        options = config_dict.pop('OPTIONS')
        global_atts = config_dict.pop('GLOBAL_ATTRIBUTES')
        if not options['regular_grid']:
            domain_dict = config_dict.pop('DOMAIN')
        else:
            domain_dict = None
        fields = config_dict

        vic2nc(options, global_atts, domain_dict, fields, big_memory)
        # ------------------------------------------------------------ #
    return
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# VIC2NC program
def vic2nc(options, global_atts, domain_dict, fields, big_memory):
    """ Convert ascii VIC files to netCDF format"""

    # determine run mode
    if big_memory \
        or not options['chunksize'] \
        or options['chunksize'] in ['all', 'All', 'ALL']:
        big_memory = True
    else:
        big_memory = False

    print("\n-------------------------------")
    print("Configuration File Options")
    print("-------------OPTIONS-------------")
    for pair in options.iteritems():
        print("%s: %s" %(pair))
    print('Fields %s' %fields.keys())
    print("-------------DOMAIN--------------")
    if domain_dict:
        for pair in domain_dict.iteritems():
            print("%s: %s" %(pair))
    print("--------GLOBAL_ATTRIBUTES--------")
    for pair in global_atts.iteritems():
        print("%s: %s" %(pair))
    print("--------RUN MODE--------")
    if big_memory:
        print('Run Mode: big_memory')
    else:
        print('Run Mode: standard (chunking)')
    print("---------------------------------\n")
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Make output directory
    if not os.path.exists(options['out_directory']):
        os.makedirs(options['out_directory'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Make pairs (i.e. find inds)
    files = glob(options['input_files'])
    points = get_file_coords(files)
    # ---------------------------------------------------------------- #


    # ---------------------------------------------------------------- #
    # Get target grid information
    if domain_dict:
        domain = read_domain(domain_dict)
        target_grid_file = path.split(domain_dict['filename'])[1]
        global_atts['target_grid_file'] = target_grid_file
    else:
        # must be a regular grid, build from file names
        domain = calc_grid(points.get_lats(), points.get_lons())
        target_grid_file = None
        domain_dict = {'y_x_dims': ['lat', 'lon']}
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Get grid index locations
    points = get_grid_inds(domain, points)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read first file and get a list of dates
    vic_datelist = get_dates(files[0])
    vic_ordtime = date2num(vic_datelist, TIMEUNITS, calendar=options['calendar'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Determine time segmentation
    if options['start_date']:
        start_date = datetime.strptime(options['start_date'], TIMESTAMPFORM)
        if start_date < vic_datelist[0]:
            print("WARNING: Start date in configuration file is before first date in file.")
            start_date = vic_datelist[0]
            print('WARNING: New start date is %s' %start_date)
    else:
        start_date = vic_datelist[0]

    if options['end_date']:
        end_date = datetime.strptime(options['end_date'], TIMESTAMPFORM)
        if end_date > vic_datelist[-1]:
            print("WARNING: End date in configuration file is after last date in file.")
            end_date = vic_datelist[-1]
            print('WARNING: New end date is %s' %end_date)
    else:
        end_date = vic_datelist[-1]

    # Ordinal Time
    start_ord = date2num(start_date, TIMEUNITS, calendar=options['calendar'])
    end_ord = date2num(end_date, TIMEUNITS, calendar=options['calendar'])

    print("netCDF Start Date: %s" %start_date)
    print("netCDF End Date: %s" %end_date)

    segment_dates = []
    if options['time_segment'] == 'day':
        # calendar insensitive
        num_segments = np.ceil(end_ord - start_ord)
        if start_date.hour == 0:
            segment_dates = num2date(np.arange(start_ord, end_ord+1, 1), TIMEUNITS, calendar=options['calendar'])
        else:
            # allow start at time other than 0
            temp = [start_ord].append(np.arange(np.ceil(start_ord), end_ord+1, 1))
            segment_dates = num2date(temp , TIMEUNITS, calendar=options['calendar'])
    elif options['time_segment'] == 'month':
        num_segments = (end_date.year - start_date.year)*12 + end_date.month - start_date.month+ 1
        month = start_date.month
        year = start_date.year
        for i in xrange(num_segments+1):
            segment_dates.append(datetime(year, month, 1))
            month += 1
            if month == 13:
                month = 1
                year += 1
    elif options['time_segment'] == 'year':
        num_segments = end_date.year - start_date.year + 1
        year = start_date.year
        for i in xrange(num_segments+1):
            segment_dates.append(datetime(year, 1, 1))
            year += 1
    elif options['time_segment'] == 'decade':
        num_segments = (end_date.year - start_date.year)/10 + 1
        year = start_date.year
        for i in xrange(num_segments+1):
            segment_dates.append(datetime(year, 1, 1))
            year += 10
    elif options['time_segment'] == 'all':
        num_segments = 1
        segment_dates=[start_date, end_date]
    else:
        raise ValueError('Unknown timesegment options %s', options['time_segment'])
    print("Number of files: %s" %(len(segment_dates)-1))
    assert len(segment_dates) == num_segments+1

    # Make sure the first and last dates are start/end_date
    segment_dates[0] = start_date
    segment_dates[-1] = end_date
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Setup Segments
    segments = np.empty(num_segments, dtype=object)

    for num in xrange(num_segments):
        # Segment time bounds
        t0 = segment_dates[num]
        t1 = segment_dates[num+1]

        # Get segment inds
        i0 = bisect_left(vic_datelist, t0)
        i1 = bisect_right(vic_datelist, t1)

        # Make segment filename (with path)
        if options['time_segment'] == 'day':
            filename = "%s.%s.nc" %(options['out_file_prefix'], t0.strftime('%Y-%m-%d'))
        elif options['time_segment'] == 'month':
            filename = "%s.%s.nc" %(options['out_file_prefix'], t0.strftime('%Y-%m'))
        elif options['time_segment'] == 'year':
            filename = "%s.%s.nc" %(options['out_file_prefix'], t0.strftime('%Y'))
        elif options['time_segment'] == 'all':
            filename = "%s.%s-%s.nc" %(options['out_file_prefix'], t0.strftime('%Y%m%d'), t1.strftime('%Y%m%d'))

        filename = path.join(options['out_directory'], filename)

        # Setup segment and initialize netcdf
        segments[num] = Segment(num, i0, i1, options['out_file_format'],
                                filename, big_memory=big_memory)
        segments[num].nc_globals(**global_atts)
        segments[num].nc_time(vic_ordtime, options['calendar'])
        segments[num].nc_dimensions(snow_bands=options['snow_bands'],
                                    veg_tiles=options['veg_tiles'],
                                    soil_layers=options['soil_layers'])

        segments[num].nc_domain(domain)
        segments[num].nc_fields(fields,
                                domain_dict['y_x_dims'], options['precision'])

        print(repr(segments[num]))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Get column numbers and names (will help speed up reading)
    names = []
    usecols = []
    for name, field in fields.iteritems():
        print name, field
        if type(field['column']) == list:
            for i, col in enumerate(field['column']):
                names.append(name+str(i))
                usecols.append(col)
        else:
            names.append(name)
            usecols.append(field['column'])


    # ---------------------------------------------------------------- #

    if big_memory:
        # run in big memory mode
        while points:
            point = points.pop()
            point.read(names=names, usecols=usecols)

            for num in xrange(num_segments):
                segments[num].nc_add_data_big_memory(point)

        for num in xrange(num_segments):
	    segments[num].nc_write_data_big_memory()

    else:
        # ------------------------------------------------------------ #
        # Chunk the input files
        point_chunks = chunks(points, int(options['chunksize']))
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Open VIC files and put data into netcdfs
        for chunk in point_chunks:
            data_points = {}
            for point in chunk:
                point.read(names=names, usecols=usecols)
                data_points[point.filename] = point

            for segment in segments:
                segment.nc_add_data_standard(data_points)
        # ------------------------------------------------------------ #

    # ---------------------------------------------------------------- #
    # Close the netcdf files
    for segment in segments:
        segment.nc_close()
    # ---------------------------------------------------------------- #
    return
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Generator to chunk points list
def chunks(l, n):
    """ Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i+n]
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Read the Configuration File
def read_config(config_file):
    """
    Return a dictionary with subdictionaries of all configFile options/values
    """
    config = SafeConfigParser()
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

    for name, section in dict1.iteritems():
        if name in default_config.keys():
            for option, key in default_config[name].iteritems():
                if option not in section.keys():
                    dict1[name][option] = key

    return dict1

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
        else:
            try:
                return int(value)
            except:
                pass
            try:
                return float(value)
            except:
                return value
    else:
        try:
            return map(int, val_list)
        except:
            pass
        try:
            return map(float, val_list)
        except:
            return val_list
# -------------------------------------------------------------------- #

def flatten(foo):
    for x in foo:
        if hasattr(x, '__iter__'):
            for y in flatten(x):
                yield y
        else:
            yield x

# -------------------------------------------------------------------- #
def get_file_coords(files):
    """
    Get list of Point objects
    """

    points = Plist()

    for i,filename in enumerate(files):
        # fname = path.split(f)[1][-16:] # just look at last 16 characters
        f = filename[-22:] # just look at last 16 characters
        lat, lon = map(float, findall(r"[-+]?\d*\.\d+|\d+", f))[-2:]
        points.append(Point(lat=lat, lon=lon, filename=filename))

    return points
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
def get_dates(file):
    """
    Read the first file in the input directory and create a ordinal based timeseries
    Also find the indicies to split the time series into months and years
    """
    hours = (0, 1, 2, 3)
    days = (0, 1, 2)
    try:
        data = np.loadtxt(file, usecols=hours, dtype=int)
        datelist = [datetime(*d) for d in data]
    except (ValueError, TypeError):
        data = np.loadtxt(file, usecols=days, dtype=int)
        datelist = [datetime(*d) for d in data]

    # check to make sure we haven't used used daily by mistake
    # (creating a bunch of duplicate times)
    newlist = []
    for i in datelist:
        if i not in newlist:
            newlist.append(i)
        else:
            raise ValueError('Found duplicate datetimes in datelist')

    print('VIC startdate: %s' %datelist[0])
    print('VIC enddate: %s' %datelist[-1])

    return datelist
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
def read_domain(domain_dict):

    print('reading domain file: %s' %domain_dict['filename'])
    f = Dataset(domain_dict['filename'])

    domain = {'lon': NcVar(f, domain_dict['longitude_var']),
              'lat': NcVar(f, domain_dict['latitude_var'])}

    if domain_dict['copy_vars']:
        for varname in domain_dict['copy_vars']:
            domain[varname] = NcVar(f, varname)

    f.close()

    return domain

# -------------------------------------------------------------------- #
def get_grid_inds(domain, points):
    """
    Find location of lat/lon points in 2d target grid.
    Uses cKdtree nearest neighbor mapping.
    """
    lons = points.get_lons()
    lats = points.get_lats()

    if (lons.min()<0 and domain['lon'].min()>=0):
        posinds = np.nonzero(lons<0)
        lons[posinds] += 360
        print('adjusted VIC lon minimum (+360 for negative lons)')

    # Make sure the longitude / latitude vars are 2d
    if domain['lat'].ndim == 1 or domain['lon'].ndim == 1:
        dlons, dlats = np.meshgrid(domain['lon'], domain['lat'])

    combined = np.dstack(([dlats.ravel(), dlons.ravel()]))[0]
    point_list = list(np.vstack((lats, lons)).transpose())

    mytree = cKDTree(combined)
    dist, indexes = mytree.query(point_list, k=1)

    yinds, xinds = np.unravel_index(indexes, dlons.shape)

    points.add_xs(xinds)
    points.add_ys(yinds)

    return points
# -------------------------------------------------------------------- #

def batch(config_file, create_batch, batch_dir):
    """Create a set of batch configuration files"""

    # Read Configuration files
    config_dict = read_config(config_file)
    options = config_dict.pop('OPTIONS')
    global_atts = config_dict.pop('GLOBAL_ATTRIBUTES')
    domain_dict = config_dict.pop('DOMAIN')
    fields = config_dict

    config = SafeConfigParser()
    config.optionxform = str

    # Figure out what to call the new files
    nameprefix = os.path.splitext(os.path.split(config_file)[1])[0]

    if create_batch == 'variables':
        # batch by variables
        # options section
        config.add_section('OPTIONS')
        for option, value in options.iteritems():
            if type(value) == list:
                try:
                    value = ", ".join(value)
                except TypeError:
                    value = ", ".join( repr(e) for e in value)
            elif type(value) != str:
                value = str(value)
            config.set('OPTIONS', option, str(value))

        # global_atts section
        config.add_section('GLOBAL_ATTRIBUTES')
        for option, value in global_atts.iteritems():
            if type(value) == list:
                try:
                    value = ", ".join(value)
                except TypeError:
                    value = ", ".join( repr(e) for e in value)
            elif type(value) != str:
                value = str(value)
            config.set('GLOBAL_ATTRIBUTES', option, str(value))

        # domain dict section
        config.add_section('DOMAIN')
        for option, value in domain_dict.iteritems():
            if type(value) == list:
                try:
                    value = ", ".join(value)
                except TypeError:
                    value = ", ".join( repr(e) for e in value)
            elif type(value) != str:
                value = str(value)

            config.set('DOMAIN', option, value.strip("'"))

        for var, field in fields.iteritems():
            suffix="_%s.cfg" %var
            new_cfg_file = os.path.join(batch_dir, nameprefix+suffix)

            # this var
            config.add_section(var)
            for option, value in field.iteritems():
                if type(value) == list:
                    try:
                        value = ", ".join(value)
                    except TypeError:
                        value = ", ".join( repr(e) for e in value)
                elif type(value) != str:
                    value = str(value)
                config.set(var, option, str(value))

            # write that config
            with open(new_cfg_file, 'wb') as cf:
                config.write(cf)

            # clear the var section
            config.remove_section(var)

    else:
        # start with existing config
        config.read(config_file)

        # by time
        start_date = datetime.strptime(options['start_date'], TIMESTAMPFORM)
        end_date = datetime.strptime(options['end_date'], TIMESTAMPFORM)

        t0 = start_date

        if create_batch == 'years':
            td = relativedelta.relativedelta(years=1)
            t1 = datetime(t0.year, 12, 31, end_date.hour)
        elif create_batch == 'months':
            td = relativedelta.relativedelta(months=1)
        elif create_batch == 'days':
            # days option is only valid for gregorian calendar
            td = relativedelta.relativedelta(days=1)

        hour = relativedelta.relativedelta(hours=-1)

        while t0 < end_date:
            t1 = t0 + td
            if t1> end_date:
                t1 = end_date
            else:
                t1 += hour

            suffix = "_%s-%s.cfg" %(t0.strftime("%Y%m%d%H"), t1.strftime("%Y%m%d%H"))
            new_cfg_file = os.path.join(batch_dir, nameprefix+suffix)

            # Write config replacing start and end dates
            config.set('OPTIONS', 'start_date', t0.strftime(TIMESTAMPFORM))
            config.set('OPTIONS', 'end_date', t1.strftime(TIMESTAMPFORM))

            with open(new_cfg_file, 'wb') as cf:
                config.write(cf)

            t0 += td
    return

# -------------------------------------------------------------------- #
# find x y coordinates
def latlon2yx(plats, plons, glats, glons):
    """find y x coordinates """

    if glons.ndim == 1 or glats.ndim == 1:
        glons, glats = np.meshgrid(glons, glats)

    combined = np.dstack(([glats.ravel(), glons.ravel()]))[0]
    points = list(np.vstack((np.array(plats), np.array(plons))).transpose())

    mytree = cKDTree(combined)
    dist, indexes = mytree.query(points, k=1)
    y, x = np.unravel_index(np.array(indexes), glons.shape)
    return y, x

# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
def calc_grid(lats, lons, decimals=4):
    """ determine shape of regular grid from lons and lats"""

    print('Calculating grid size now...')

    target_grid = {}

    # get unique lats and lons
    lon = np.sort(np.unique(lons.round(decimals=decimals)))
    print('found %s unique lons' %len(lon))

    lat = np.sort(np.unique(lats.round(decimals=decimals)))
    print('found %s unique lats' %len(lat))

    y, x = latlon2yx(lats, lons, lat, lon)

    mask = np.zeros((len(lat), len(lon)))

    mask[y, x] = 1

    # Create fake NcVar Types
    target_grid['lon'] = FakeNcVar(lon, ('lon', ), {'long_name': 'longitude coordinate',
                                                    'units': 'degrees_east'})
    target_grid['lat'] = FakeNcVar(lat, ('lat', ), {'long_name': 'latitude coordinate',
                                                    'units': 'degrees_north'})
    target_grid['mask'] = FakeNcVar(mask, ('lat', 'lon', ), {'long_name': 'domain mask',
                                                             'comment': '0 indicates grid cell is not active'})

    print('Created a target grid based on the lats and lon in the input file names')
    print('Grid Size: {}'.format(mask.shape))

    return target_grid
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
def process_command_line():
    """
    Get the path to the config_file
    """
    # Parse arguments
    parser = ArgumentParser(description='convert VIC ascii output to netCDF format')
    parser.add_argument("config_file", type=str,
                        help="Input configuration file")
    parser.add_argument("--big_memory", action='store_true',
                        help="Operate in high memory mode (all data will be stored in memory until write time)")
    parser.add_argument("--create_batch", type=str, choices=['days', 'months', 'years', 'variables'],
                        default=False, help="Create a batch of config files")
    parser.add_argument("--batch_dir", type=str, default="./",
                        help="Location to put batch config files")
    args = parser.parse_args()

    if not os.path.isfile(args.config_file):
        raise IOError('Configuration File: %s is not a valid file' %(args.config_file))

    if not os.path.isdir(args.batch_dir) and args.create_batch:
        raise IOError('Configuration File: %s is not a valid file' %(args.config_file))

    return args.config_file, args.big_memory, args.create_batch, args.batch_dir
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
if __name__ == "__main__":
    main()
# -------------------------------------------------------------------- #
