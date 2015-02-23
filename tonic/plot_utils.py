"""plot_utils"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.basemap import Basemap
import numpy as np

projections = {}
projections['wr50a'] = {'urcrnrlat': 27.511827255753555,
                        'urcrnrlon': 16.90845094934209,
                        'llcrnrlat': 16.534986367884521,
                        'llcrnrlon': 189.2229322311162,
                        'projection': 'lcc',
                        'rsphere': 6371200.0,
                        'lon_0': -114,
                        'lat_0': 90}
default_projection = projections['wr50a']

default_map = None  # set again at end of module


class Bmap(object):
    """Wrapper container for Basemap object and map indices"""
    def __init__(self, projection=default_projection):
        self.projection = projection
        self.m = Basemap(**self.projection)

    def set_map_inds(self, lons, lats):
        """Set the map indices
        lons : numpy.ndarray
            Array of grid cell longitudes.
        lons : numpy.ndarray
            Array of grid cell latitudes.
        """
        self.xi, self.yi = self.m(np.squeeze(lons), np.squeeze(lats))
        self.inds_set = True


def make_bmap(projection=default_projection, lons=None, lats=None):
    """
    Make a bmap object storing the projection and basemap indices for a series
    of plots.

    Parameters
    ----------
    projection : dict
        Projection keywords to be passed to `Basemap`.
    lons : numpy.ndarray
        Array of grid cell longitudes.
    lons : numpy.ndarray
        Array of grid cell latitudes.

    Returns
    ----------
    bmap : Bmap
        Bmap object.
    """
    bmap = Bmap(projection=projection)

    if lons is not None and lats is not None:
        bmap.set_map_inds(lons, lats)
    return bmap


# -------------------------------------------------------------------- #
def cmap_discretize(cmap, n_colors=10):
    """Return discretized colormap.

    Parameters
    ----------
    cmap : str or colormap object
        Colormap to discretize.
    n_colors : int
        Number of discrete colors to divide `cmap` into.

    Returns
    ----------
    disc_cmap : LinearSegmentedColormap
        Discretized colormap.
    """
    try:
        cmap = cm.get_cmap(cmap)
    except:
        cmap = cm.get_cmap(eval(cmap))
    colors_i = np.concatenate((np.linspace(0, 1., n_colors), (0., 0., 0., 0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., n_colors + 1)
    cdict = {}
    for ki, key in enumerate(('red', 'green', 'blue')):
        cdict[key] = [(indices[i], colors_rgba[i - 1, ki], colors_rgba[i, ki])
                      for i in range(n_colors + 1)]

    return mpl.colors.LinearSegmentedColormap(cmap.name + "_%d" % n_colors,
                                              cdict, 1024)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #

def sub_plot_pcolor(data,
                    title=None,
                    cmap='Spectral_r',
                    vmin=None,
                    vmax=None,
                    cbar=True,
                    cbar_location='bottom',
                    units=None,
                    cbar_extend='neither',
                    map_obj=default_map,
                    ax=None):
    """Plot data into a subplot using pcolormesh.

    Parameters
    ----------
    data : 2d numpy array
        Array of data to plot using pcolormesh.
    title : str, optional
        Title to add to subplot.
    cmap : str or colormap object, optional
        Colormap to use in pcolormesh.
    vmin : scalar, optional
        Minimum value for color range.
    vmax : scalar, optional
        Maximum value for color range.
    cbar_location : str, optional
        Location of colorbar.  Value is passed to `Basemap.colorbar()`
    cbar_extend : str, optional
        Extend colorbar ends.  Value is passed to `Basemap.colorbar()`
    units : str
        Colorbar label.
    map_obj : Bmap, optional
        `Bmap` object containing `Basemap` object as well as plot coordinates.
        Set using `make_bmap()`.
    ax : axis object, optional
        Axis object to use for subplot.  If None, `ax=plt.gca()`.
    """

    if vmin is None:
        vmin = data.min()
    if vmax is None:
        vmax = data.max()

    if ax is None:
        ax = plt.gca()

    map_obj.m.pcolormesh(map_obj.xi, map_obj.yi, np.squeeze(data),
                         vmin=vmin, vmax=vmax, cmap=cmap, ax=ax)
    map_obj.m.drawparallels(np.arange(-80., 81., 20.))
    map_obj.m.drawmeridians(np.arange(-180., 181., 20.))
    map_obj.m.drawcoastlines(color='k', linewidth=0.25)

    if title is not None:
        plt.title(title, size=13)
    if cbar is not None:
        cbar = map_obj.m.colorbar(location=cbar_location, extend=cbar_extend)
    cbar.set_label(units)

    return
# -------------------------------------------------------------------- #

default_map = make_bmap()
