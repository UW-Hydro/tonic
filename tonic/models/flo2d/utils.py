from osgeo import osr
import mpl_toolkits.basemap.pyproj as pyproj

import numpy as np
import pandas as pd

from scipy import stats
import quantities as pq
from scipy.spatial import cKDTree


default_proj = pyproj.Proj('+proj=lcc +lat_1=47.5 +lat_2=48.73333333333333 '
                           '+lat_0=47'' +lon_0=-120.8333333333333 '
                           '+x_0=609601.2192024385 +y_0=0 +datum=NAD83 '
                           '+units=us-ft')


def esriprj2standards(shapeprj_path, kind='proj4'):
    prj_file = open(shapeprj_path, 'r')
    prj_txt = prj_file.read()
    srs = osr.SpatialReference()
    srs.ImportFromESRI([prj_txt])
    if kind.lower() == 'prj':
        return prj_txt
    if kind.lower() == 'wtk':
        return srs.ExportToWkt()
    elif kind.lower() == 'proj4':
        return srs.ExportToProj4()
    elif kind.lower() == 'epsg':
        srs.AutoIdentifyEPSG()
        return srs.GetAuthorityCode(None)
    else:
        raise ValueError('Kind is unknown')


def read_flo2d_depth_file(fname):
    return pd.read_table(fname, sep='\s+',
                         names=['x', 'y', 'depth'],
                         index_col=0)


def make_coordinates(arr):
    min_arr = arr.min()
    max_arr = arr.max()
    d_arr = stats.mode(np.diff(sorted(arr.unique())))[0][0]
    return np.arange(min_arr, max_arr + d_arr, d_arr)


def grid_flo2d_depth(fname):
    df = read_flo2d_depth_file(fname)

    xs = make_coordinates(df.x)
    ys = make_coordinates(df.y)
    gys, gxs = np.meshgrid(ys, xs)

    combined = np.dstack(([gys.ravel(), gxs.ravel()]))[0]
    points = list(np.vstack((np.array(df.y.values),
                             np.array(df.x.values))).transpose())

    mytree = cKDTree(combined)
    dist, indexes = mytree.query(points, k=1)
    y, x = np.unravel_index(indexes, gxs.shape)

    data = np.zeros_like(gys)
    data[y, x] = df.depth

    return data, gys, gxs


def flo2d_coords_to_geographic(x, y, proj=default_proj):
    x_ft = x * pq.ft
    x_m = x_ft.rescale(pq.m)

    y_ft = y * pq.ft
    y_m = y_ft.rescale(pq.m)

    return proj(x_m, y_m, inverse=True)
