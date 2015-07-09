from osgeo import osr
import mpl_toolkits.basemap.pyproj as pyproj

import numpy as np
import pandas as pd

from scipy import stats
import quantities as pq
from scipy.spatial import cKDTree

# Default projection from Alan Hamlet 2015-06-29
default_proj = pyproj.Proj('+proj=lcc +lat_1=47.5 +lat_2=48.73333333333333 '
                           '+lat_0=47'' +lon_0=-120.8333333333333 '
                           '+x_0=609601.2192024385 +y_0=0 +datum=NAD83 '
                           '+units=us-ft')


def esriprj2standards(shapeprj_path, kind='proj4'):
    '''Return string representation of projection in shapeprj_path.

    Parameters
    ----------
    shapeprj_path : str
        Filepath to ESRI style projection file.
    kind : str {'prj', 'wkt', 'proj4', 'epsg'}
        Projection standard output type.

    Returns
    ----------
    projection_params : str
        Projection parameters.
    '''
    prj_file = open(shapeprj_path, 'r')
    prj_txt = prj_file.read()
    srs = osr.SpatialReference()
    srs.ImportFromESRI([prj_txt])
    if kind.lower() == 'prj':
        return prj_txt
    if kind.lower() == 'wkt':
        return srs.ExportToWkt()
    elif kind.lower() == 'proj4':
        return srs.ExportToProj4()
    elif kind.lower() == 'epsg':
        srs.AutoIdentifyEPSG()
        return srs.GetAuthorityCode(None)
    else:
        raise ValueError('Kind is unknown')


def read_flo2d_depth_file(fname, names=('x', 'y', 'depth')):
    '''Parser helper function to read ASCII output from FLO-2D

    Parameters
    ----------
    fname : str
        Filepath to FLO-2D depth file.
    names : tuple of str, optional
        Header names

    Returns
    ----------
    df : Pandas.DataFrame
        DataFrame of FLO-2D output.
    '''
    return pd.read_table(fname, sep='\s+',
                         names=names,
                         index_col=0)


def make_coordinates(arr):
    '''Calculate 1D coordinate vector from array of coordinates'''
    min_arr = arr.min()
    max_arr = arr.max()
    d_arr = stats.mode(np.diff(sorted(arr.unique())))[0][0]
    return np.arange(min_arr, max_arr + d_arr, d_arr)


def grid_flo2d_depth(fname):
    '''Read and grid FLO-2D depth file output

    Parameters
    ----------
    fname : str
        Filepath to FLO-2D depth file.

    Returns
    ----------
    data: 2d numpy array
        Gridded depth data
    gys: 2d numpy array
        Gridded y coordinates
    gxs: 2d numpy array
        Gridded x coordinates
    '''
    df = read_flo2d_depth_file(fname)

    xs = make_coordinates(df.x)
    ys = make_coordinates(df.y)
    gys, gxs = np.meshgrid(ys, xs)

    combined = np.dstack(([gys.ravel(), gxs.ravel()]))[0]
    points = list(np.vstack((df.y.values, df.x.values)).transpose())

    mytree = cKDTree(combined)
    _, indexes = mytree.query(points, k=1)
    y, x = np.unravel_index(indexes, gxs.shape)

    data = np.zeros_like(gys)
    data[y, x] = df.depth

    return data, gys, gxs


def flo2d_coords_to_geographic(x, y, proj=default_proj):
    '''Convert local projection coordinates to geographic (lat/lon) coordinates

    Parameters
    ----------
    x : array_like
        x coordinates, units=ft
    y : array_like
        y coordinates, units=ft
    proj : pyproj.Proj
        Projection object

    Returns
    ----------
    lons: array_like
        longitude values corresponding to x
    lats: array_like
        latitude values corresponding to y
    '''
    x_ft = x * pq.ft
    x_m = x_ft.rescale(pq.m)

    y_ft = y * pq.ft
    y_m = y_ft.rescale(pq.m)

    return proj(x_m, y_m, inverse=True)
