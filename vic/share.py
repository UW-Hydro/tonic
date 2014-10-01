""""""
from netCDF4 import Dataset
from collections import OrderedDict
import numpy as np
try:
    from configparser import SafeConfigParser
except:
    from ConfigParser import SafeConfigParser
from scipy.spatial import cKDTree


MMPERMETER = 1000.


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
        if obj is None:
            return
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
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
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
def read_netcdf(nc_file, variables=[], coords=False, verbose=True):
    """
    Read data from input netCDF. Will read all variables if none provided.
    Will also return all variable attributes.
    Both variables (data and attributes) are returned as dictionaries named by
    variable.
    """

    f = Dataset(nc_file, 'r')

    if variables == []:
        variables = f.variables.keys()

    if verbose:
        print('Reading input data variables: '
              ' {0} from file: {1}'.format(variables, nc_file))

    d = OrderedDict()
    a = OrderedDict()

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
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
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
    print('found {0} unique lons'.format(len(lon)))

    lat = np.sort(np.unique(lats.round(decimals=decimals)))
    print('found {0} unique lats'.format(len(lat)))

    y, x = latlon2yx(lats, lons, lat, lon)

    mask = np.zeros((len(lat), len(lon)), dtype=int)

    mask[y, x] = 1

    # Create fake NcVar Types
    target_grid['lon'] = FakeNcVar(lon, ('lon', ),
                                   {'long_name': 'longitude coordinate',
                                    'units': 'degrees_east'})
    target_grid['lat'] = FakeNcVar(lat, ('lat', ),
                                   {'long_name': 'latitude coordinate',
                                    'units': 'degrees_north'})
    target_grid['mask'] = FakeNcVar(mask, ('lat', 'lon', ),
                                    {'long_name': 'domain mask',
                                     'comment': '0 indicates grid cell is not \
                                     active'})

    print('Created a target grid based on the lats and lons in the '
          'input file names')
    print('Grid Size: {}'.format(mask.shape))

    return target_grid
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def get_grid_inds(domain, points):
    """
    Find location of lat/lon points in 2d target grid.
    Uses cKdtree nearest neighbor mapping.
    """
    lons = points.get_lons()
    lats = points.get_lats()

    if (lons.min() < 0) and (domain['lon'].min() >= 0):
        posinds = np.nonzero(lons < 0)
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
