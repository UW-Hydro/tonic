"""plot_utils"""
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as mcolors
import numpy as np


# -------------------------------------------------------------------- #
def cmap_discretize(cmap, N=9):
    """Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet.
        N: number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    if type(cmap) == str:
        cmap = cm.get_cmap(eval(cmap))
    colors_i = np.concatenate((np.linspace(0, 1., N), (0., 0., 0., 0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki, key in enumerate(('red', 'green', 'blue')):
        cdict[key] = [(indices[i], colors_rgba[i-1, ki], colors_rgba[i, ki])
                      for i in xrange(N+1)]
    # Return colormap object.
    return mcolors.LinearSegmentedColormap(cmap.name + "_%d" % N, cdict, 1024)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def sub_plot_pcolor(lons, lats, data, title=None, cmap=cm.jet,
                    vmin=None, vmax=None, cbar=True, cbar_location='bottom',
                    units=None, ncolors=5, projection=None):

    if vmin is None:
        vmin = data.min()
    if vmax is None:
        vmax = data.max()

    m = Basemap(**projection)
    m.drawlsmask(land_color='grey', lakes=False)
    xi, yi = m(np.squeeze(lons), np.squeeze(lats))
    sp = m.pcolormesh(xi, yi, np.squeeze(data), vmin=vmin, vmax=vmax,
                      cmap=cmap)
    m.drawparallels(np.arange(-80., 81., 20.))
    m.drawmeridians(np.arange(-180., 181., 20.))
    m.drawcoastlines(color='k', linewidth=0.25)
    if title:
        plt.title(title, size=13)
    if cbar:
        cmap = cmap_discretize(cmap, ncolors)
        mappable = cm.ScalarMappable(cmap=cmap)
        mappable.set_array([])
        mappable.set_clim(-0.5, ncolors+0.5)
        colorbar = m.colorbar(mappable, location=cbar_location)
        colorbar.set_ticks(np.linspace(0, ncolors, ncolors))
        colorbar.set_ticklabels(np.linspace(vmin, vmax, ncolors))
        colorbar.set_label(units)
    return sp
# -------------------------------------------------------------------- #
