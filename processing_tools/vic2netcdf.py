#!/usr/bin/env python
import traceback
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
from multiprocessing.pool import Pool
import multiprocessing
multiprocessing.log_to_stderr()
import sys
import numpy as np
import time as tm
'''
Development plan
1.  read configuration
2.  create points dict
3.  create GLOBAL dict of netcdfs (for each segment)
4.  use Pool to split file Point list by chunksize
5.  for each point,
    - open file
    - loop over segments, writing data to put appropriate file and variables
    - close file
6.  close GLOBAL dict of netcdfs
'''

REFERENCE_STRING = '0001-1-1 0:0:0'
TIMEUNITS = 'days since ' + REFERENCE_STRING     # do not change (MUST BE DAYS)!
TIMESTAMPFORM = '%Y-%m-%d-%H'

# precision
NC_DOUBLE = 'f8'
NC_FLOAT = 'f4'

# Global declaration of data list
data_points = []
segments = {}

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
        print 'reading', self.filename
        self.df = read_csv(self.filename,
                           delimiter='\t',
                           header=None,
                           usecols=usecols)
        self.df.columns = names
        return

    def __str__(self):
        return "Point(%s,%s,%s,%s)" % (self.lat, self.lon, self.y, self.x)

    def __repr__(self):
        return '__repr__'
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
    def __init__(self, num, i0, i1, filename):
        '''Class used for holding segment information '''
        self.num = num
        self.i0 = i0
        self.i1 = i1
        self.filename = filename
        self.fields = {}

        self.nc_write()


    def nc_globals(self,
                   title='VIC netCDF file',
                   history='Created: {} by {}'.format(tm.ctime(tm.time()), getuser()),
                   institution='University of Washington',
                   source=sys.argv[0],
                   references='Primary Historical Reference for VIC: Liang, X., D. P. Lettenmaier, E. F. Wood, and S. J. Burges, 1994: A Simple hydrologically Based Model of Land Surface Water and Energy Fluxes for GSMs, J. Geophys. Res., 99(D7), 14,415-14,428.',
                   comment='Output from the Variable Infiltration Capacity (VIC) Macroscale Hydrologic Model',
                   Conventions='CF-1.6',
                   target_grid_file='unknown',
                   **kwargs):

        self.f.title = title
        self.f.history = history
        self.f.institution = institution
        self.f.source = source
        self.f.references = references
        self.f.comment = comment
        self.f.Conventions = Conventions

        for attribute, value in kwargs.iteritems():
            if hasattr(self.f, attribute):
                print 'WARNING: Attribute %s already exists.'
                print 'Renaming to g_%s to ovoid overwriting.' %(attribute, attribute)
                attribute='g_'+attribute
            setattr(self.f, attribute, value)
        return

    def __str__(self):
        return "Segment(%s,%s)" % (self.i0, self.i1)

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

        for name, field in fields.iteritems():

            if 'type' in field.keys():
                prec = field['type']
            else:
                prec = prec_global
            fill_val = default_fillvals[prec]

            self.fields[name] = self.f.createVariable(name,
                                                      prec,
                                                      coords,
                                                      fill_value=fill_val)
            self.fields[name][:, :] = fill_val

            if 'units' in field.keys():
                self.fields[name].long_name = name
                for key, val in field.iteritems():
                    setattr(self.fields[name], key, val)
            else:
                raise ValueError('Field %s missing units attribute', name)
        return

    def nc_add_data(self, data):
        print 'adding data to %s' %self.filename
        for point in data:
            for name in self.var_list:
                self.f.variables[name][:, point.y, point.x]
        return

    def nc_write(self):
        self.f = Dataset(self.filename, mode="w", clobber=True, format='NETCDF4_CLASSIC')
        print 'Opened in write mode: %s' %self.filename
        return

    def nc_append(self):
        self.f = Dataset(self.filename, mode='a')
        print 'Opened in append mode: %s' %self.filename
        return

    def nc_close(self):
        self.f.close()
        print 'Closed: %s' %self.filename
        return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
class NcVar(np.ndarray):
    """ Subclass of numpy array to carry netcdf attributes"""
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

# Shortcut to multiprocessing's logger
def error(msg, *args):
    return multiprocessing.get_logger().error(msg, *args)

class LogExceptions(object):
    def __init__(self, callable):
        self.__callable = callable
        return

    def __call__(self, *args, **kwargs):
        try:
            result = self.__callable(*args, **kwargs)

        except Exception as e:
            # Here we add some debugging help. If multiprocessing's
            # debugging is on, it will arrange to log the traceback
            error(traceback.format_exc())
            # Re-raise the original exception so the Pool worker can
            # clean up
            raise

        # It was fine, give a normal answer
        return result
    pass

class LoggingPool(Pool):
    def apply_async(self, func, args=(), kwds={}, callback=None):
        return Pool.apply_async(self, LogExceptions(func), args, kwds, callback)


# -------------------------------------------------------------------- #
# Top level driver
def main():

    # ---------------------------------------------------------------- #
    # Read command Line
    config_file, numofproc = process_command_line()
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read Configuration files
    config_dict = read_config(config_file)
    options = config_dict.pop('OPTIONS')
    global_atts = config_dict.pop('GLOBAL_ATTRIBUTES')
    domain_dict = config_dict.pop('DOMAIN')
    fields = config_dict

    print "\n-------------------------------"
    print "Configuration File Options"
    print "-------------OPTIONS-------------"
    for pair in options.iteritems():
        print "%s: %s" %(pair)
    print 'Fields %s' %fields.keys()
    print "-------------DOMAIN--------------"
    for pair in domain_dict.iteritems():
        print "%s: %s" %(pair)
    print "--------GLOBAL_ATTRIBUTES--------"
    for pair in global_atts.iteritems():
        print "%s: %s" %(pair)
    print "---------------------------------\n"

    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Get target grid information
    domain = read_domain(domain_dict)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # make pairs (i.e. find inds)
    files = glob(options['input_files'])
    points = get_file_coords(files)
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
            print "WARNING: Start date in configuration file is before first date in file."
            start_date = vic_datelist[0]
            print 'WARNING: New start date is %s' %start_date
    else:
        start_date = vic_datelist[0]

    if options['end_date']:
        end_date = datetime.strptime(options['end_date'], TIMESTAMPFORM)
        if end_date > vic_datelist[-1]:
            print "WARNING: End date in configuration file is after last date in file."
            end_date = vic_datelist[-1]
            print 'WARNING: New end date is %s' %end_date

    # Ordinal Time
    start_ord = date2num(start_date, TIMEUNITS, calendar=options['calendar'])
    end_ord = date2num(end_date, TIMEUNITS, calendar=options['calendar'])

    print "netCDF Start Date:", start_date
    print "netCDF End Date:", end_date

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
    print "Number of files: ", len(segment_dates)-1
    assert len(segment_dates) == num_segments+1

    # Make sure the first and last dates are start/end_date
    segment_dates[0] = start_date
    segment_dates[-1] = end_date

    target_grid_file = path.split(domain_dict['filename'])[1]

    if numofproc > 1:
        pool0 = LoggingPool(processes=numofproc)

        for num in xrange(len(segment_dates)-1):
            # Segment time bounds
            t0 = segment_dates[num]
            t1 = segment_dates[num+1]

            # Get segment inds
            i0 = bisect_right(vic_datelist, t0)
            i1 = bisect_left(vic_datelist, t1)

            # Make segment filename (with path)
            filename = "%s.%s-%s.nc" %(options['out_file_prefix'], t0.strftime('%Y%m%d'), t1.strftime('%Y%m%d'))
            filename = path.join(options['out_directory'], filename, )

            # Setup segment and initialize netcdf
            pool0.apply_async(setup_segment,
                             args=(num, i0, i1, filename, global_atts,
                                   vic_ordtime, options, domain_dict, fields),
                             callback=store_segment)
        pool0.close()
        pool0.join()
    else:
        for num in xrange(len(segment_dates)-1):
            # Segment time bounds
            t0 = segment_dates[num]
            t1 = segment_dates[num+1]

            # Get segment inds
            i0 = bisect_right(vic_datelist, t0)
            i1 = bisect_left(vic_datelist, t1)

            # Make segment filename (with path)
            filename = "%s.%s-%s.nc" %(options['out_file_prefix'], t0.strftime('%Y%m%d'), t1.strftime('%Y%m%d'))
            filename = path.join(options['out_directory'], filename, )

            # Setup segment and initialize netcdf
            segments[num] = Segment(num, i0, i1, filename)
            segments[num].nc_globals(target_grid_file=target_grid_file,
                                     **global_atts)
            segments[num].nc_time(vic_ordtime, options['calendar'])
            segments[num].nc_domain(domain)
            segments[num].nc_fields(fields, domain_dict['y_x_dims'], options['precision'])

    for num, segment in segments.iteritems():
        print "\n-------------------------- Segment %s --------------------------" %num
        print "Filename: %s" %segment.filename
        print "Start Index: %s" %segment.i0
        print "End Index: %s" %segment.i1
        print "------------------------------------------------------------------"
    # ---------------------------------------------------------------- #

    # After profiling, writing seems to be taking about 10 times as long as reading
    # We cant pass around file descriptors so instead, we'll use multiprocessing to
    # read big chunks of files and then write them to each segment.
    # It also turns out that writing is slow enough that it makes sense to sacrafice
    # opening and closing the nc files so that we can use the other processors.

    # ---------------------------------------------------------------- #
    # Get column numbers and names (will help speed up reading)
    names = []
    usecols = []
    for name, field in fields.iteritems():
        names.append(name)
        usecols.append(field['column'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Get chunksize
    point_chunks = chunks(points, int(options['chunksize']))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Open VIC files and put in netcdfs
    global data_points
    if numofproc > 1:
        for chunk in point_chunks:
            # read chunk of points
            pool0 = LoggingPool(processes=numofproc)
            for point in chunk:
                pool0.apply_async(read_point,
                                 args=(point, names, usecols),
                                 callback=store_data)
            pool0.close()
            pool0.join()
            print 'done reading points for this chunk'
            # Use a new pool for writing
            pool2 = LoggingPool(processes=numofproc)
            print 'setup pool2 for writing / appending data'
            for num, segment in segments.iteritems():
                pool2.apply_async(add_data,
                                 args=(data_points, segment),
                                 callback=store_segment)
            pool2.close()
            pool2.join()
            data_points = []
    else:
        for chunk in point_chunks:
            for point in chunk:
                point.read(names=names, usecols=usecols)
                data_points.append(point)

            for num, segment in segments.iteritems():
                segment.nc_add_data(data_points)
            data_points = []
        # Close the netcdf files
        for num, segment in segments.iteritems():
            segment.nc_close()
    # ---------------------------------------------------------------- #
    return

# -------------------------------------------------------------------- #
# Helper functions to support multiprocessing
# wrapper function for seting up netcdf segments
def setup_segment(num, i0, i1, filename, global_atts,
                  vic_ordtime, options, domain_dict, fields):
    domain = read_domain(domain_dict)
    target_grid_file = path.split(domain_dict['filename'])[1]
    seg = Segment(num, i0, i1, filename)
    seg.nc_globals(target_grid_file=target_grid_file,
                             **global_atts)
    seg.nc_time(vic_ordtime, options['calendar'])
    seg.nc_domain(domain)
    seg.nc_fields(fields, domain_dict['y_x_dims'], options['precision'])
    seg.nc_close()
    return seg

def add_data(data_points, seg):
    seg.nc_append()
    seg.nc_add_data(data_points)
    seg.nc_close()
    return seg

# Callback function for setup_segment
def store_segment(seg):
    global segments
    segments[seg.num] = seg

# point.read wrapper function to ovoid pickling error
def read_point(point, names, usecols):
    point.read(names=names, usecols=usecols)
    return point

# Callback function for read_point
def store_data(point):
    global data_points
    data_points.append(point)

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


# -------------------------------------------------------------------- #
def process_command_line():
    """
    Get the path to the config_file
    """
    # Parse arguments
    parser = ArgumentParser(description='convert VIC ascii output to netCDF format')
    parser.add_argument("config_file", type=str,
                        help="Input configuration file")
    parser.add_argument("-np", "--numofproc", type=int,
                        help="Number of processors used to run job", default=1)

    args = parser.parse_args()

    return args.config_file, args.numofproc
# -------------------------------------------------------------------- #


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
    data = np.loadtxt(file, usecols=(0, 1, 2, 3), dtype=int)
    datelist = [datetime(*d) for d in data]

    print 'VIC startdate:', datelist[0]
    print 'VIC enddate:', datelist[-1]

    return datelist
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
def read_domain(domain_dict):

    print 'reading domain file: %s' %domain_dict['filename']
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
        print 'adjusted VIC lon minimum (+360 for negative lons)'

    combined = np.dstack(([domain['lat'].ravel(), domain['lon'].ravel()]))[0]
    point_list = list(np.vstack((lats, lons)).transpose())

    mytree = cKDTree(combined)
    dist, indexes = mytree.query(point_list, k=1)

    yinds, xinds = np.unravel_index(indexes, domain['lon'].shape)

    points.add_xs(xinds)
    points.add_ys(yinds)

    return points
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
if __name__ == "__main__":
    main()
# -------------------------------------------------------------------- #
