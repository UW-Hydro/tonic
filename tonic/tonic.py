"""tonic"""
import numpy as np
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
