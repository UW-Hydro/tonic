#!/usr/bin/env python
"""
grid_parmas.pyc

A toolkit for converting classic vic parameters to netcdf format
"""

from __future__ import print_function
import sys
import numpy as np
from netCDF4 import Dataset, default_fillvals
from scipy.spatial import cKDTree
from scipy import stats
import time as tm
import argparse
import socket
from getpass import getuser
from collections import OrderedDict

# -------------------------------------------------------------------- #
# precision
PRECISION = 1.0e-30
NC_DOUBLE = 'f8'
NC_FLOAT = 'f4'
NC_INT = 'i4'
NC_CHAR = 'S1'
MAX_NC_CHARS = 256

# fill values
FILLVALUE_F = default_fillvals[NC_DOUBLE]
FILLVALUE_I = default_fillvals[NC_INT]

XVAR = 'lon'
YVAR = 'lat'

# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
class cols(object):
    def __init__(self, nlayers=3, snow_bands=5):
        self.soil_param = OrderedDict([('run_cell', np.array([0])),
                                       ('gridcell', np.array([1])),
                                       ('lats', np.array([2])),
                                       ('lons', np.array([3])),
                                       ('infilt', np.array([4])),
                                       ('Ds', np.array([5])),
                                       ('Dsmax',  np.array([6])),
                                       ('Ws', np.array([7])),
                                       ('c', np.array([8])),
                                       ('expt', np.arange(9, nlayers+9)),
                                       ('Ksat', np.arange(nlayers+9,
                                                          2*nlayers+9)),
                                       ('phi_s', np.arange(2*nlayers+9,
                                                           3*nlayers+9)),
                                       ('init_moist', np.arange(3*nlayers+9,
                                                                4*nlayers+9)),
                                       ('elev', np.array([4*nlayers+9])),
                                       ('depth', np.arange(4*nlayers+10,
                                                           5*nlayers+10)),
                                       ('avg_T', np.array([5*nlayers+10])),
                                       ('dp', np.array([5*nlayers+11])),
                                       ('bubble', np.arange(5*nlayers+12,
                                                            6*nlayers+12)),
                                       ('quartz', np.arange(6*nlayers+12,
                                                            7*nlayers+12)),
                                       ('bulk_density', np.arange(7*nlayers+12,
                                                                  8*nlayers+12)),
                                       ('soil_density', np.arange(8*nlayers+12,
                                                                  9*nlayers+12)),
                                       ('off_gmt', np.array([9*nlayers+12])),
                                       ('Wcr_FRACT', np.arange(9*nlayers+13,
                                                                10*nlayers+13)),
                                       ('Wpwp_FRACT', np.arange(10*nlayers+13,
                                                                11*nlayers+13)),
                                       ('rough', np.array([11*nlayers+13])),
                                       ('snow_rough', np.array([11*nlayers+14])),
                                       ('annual_prec',np.array([11*nlayers+15])),
                                       ('resid_moist', np.arange(11*nlayers+16,
                                                                 12*nlayers+16)),
                                       ('fs_active',  np.array([12*nlayers+16]))])

        self.snow_param = OrderedDict([('cellnum', np.array([0])),
                                       ('AreaFract', np.arange(1,
                                                               snow_bands+1)),
                                       ('elevation', np.arange(snow_bands+1,
                                                               2*snow_bands+1)),
                                       ('Pfactor', np.arange(2*snow_bands+1,
                                                             3*snow_bands+1))])

        self.veglib = OrderedDict([('Veg_class', np.array([0])),
                                   ('lib_overstory', np.array([1])),
                                   ('lib_rarc', np.array([2])),
                                   ('lib_rmin', np.array([3])),
                                   ('lib_LAI', np.arange(4, 16)),
                                   ('lib_albedo', np.arange(16, 28)),
                                   ('lib_rough', np.arange(28, 40)),
                                   ('lib_displacement', np.arange(40, 52)),
                                   ('lib_wind_h', np.array(52)),
                                   ('lib_RGL', np.array(53)),
                                   ('lib_rad_atten', np.array(54)),
                                   ('lib_wind_atten', np.array(55)),
                                   ('lib_trunk_ratio', np.array(56))])
                        #'comment',  [57]}
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
class format(object):
    def __init__(self, nlayers=3, snow_bands=5):
        self.soil_param = {'run_cell': '%1i',
                           'gridcell': '%1i',
                           'lats': '%1.4f',
                           'lons': '%1.4f',
                           'infilt': '%1.6f',
                           'Ds': '%1.6f',
                           'Dsmax': '%1.6f',
                           'Ws': '%1.6f',
                           'c': '%1.6f',
                           'expt': '%1.6f',
                           'Ksat': '%1.6f',
                           'phi_s': '%1.6f',
                           'init_moist': '%1.6f',
                           'elev': '%1.6f',
                           'depth': '%1.6f',
                           'avg_T': '%1.6f',
                           'dp': '%1.6f',
                           'bubble': '%1.6f',
                           'quartz': '%1.6f',
                           'bulk_density': '%1.6f',
                           'soil_density': '%1.6f',
                           'off_gmt': '%1.6f',
                           'Wcr_FRACT': '%1.6f',
                           'Wpwp_FRACT': '%1.6f',
                           'rough': '%1.6f',
                           'snow_rough': '%1.6f',
                           'annual_prec': '%1.6f',
                           'resid_moist': '%1.6f',
                           'fs_active': '%1i'}

        self.snow_param = {'cellnum': '%1i',
                           'AreaFract': '%1.6f',
                           'elevation': '%1.6f',
                           'Pfactor': '%1.6f'}

        self.veglib = {'Veg_class': '%1i',
                       'lib_overstory': '%1.6f',
                       'lib_rarc': '%1.6f',
                       'lib_rmin': '%1.6f',
                       'lib_LAI': '%1.6f',
                       'lib_albedo': '%1.6f',
                       'lib_rough': '%1.6f',
                       'lib_displacement': '%1.6f',
                       'lib_wind_h': '%1.6f',
                       'lib_RGL': '%1.6f',
                       'lib_rad_atten': '%1.6f',
                       'lib_wind_atten': '%1.6f',
                       'lib_trunk_ratio': '%1.6f'}
                        #'comment': [57]}
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
class description(object):
    def __init__(self, ORGANIC_FRACT=False, SPATIAL_FROST=False,
                 SPATIAL_SNOW=False, EXCESS_ICE=False,
                 JULY_TAVG_SUPPLIED=False, BLOWING_SNOW=False,
                 GLOBAL_LAI=False):
        self.soil_param = {'run_cell': '1 = Run Grid Cell, 0 = Do Not Run',
                           'gridcell': 'Grid cell number',
                           'lats': 'Latitude of grid cell',
                           'lons': 'Longitude of grid cell',
                           'infilt': 'Variable infiltration curve parameter '
                                     '(binfilt)',
                           'Ds': 'Fraction of Dsmax where non-linear baseflow '
                                 'begins',
                           'Dsmax': 'Maximum velocity of baseflow',
                           'Ws': 'Fraction of maximum soil moisture where '
                                 'non-linear baseflow occurs',
                           'c': 'Exponent used in baseflow curve, normally set'
                                ' to 2',
                           'expt': 'Exponent n (=3+2/lambda) in Campbells eqn'
                                   ' for hydraulic conductivity, HBH 5.6 '
                                   '(where lambda = soil pore size '
                                   'distribution parameter).  Values should be'
                                   ' > 3.0.',
                           'Ksat': 'Saturated hydrologic conductivity',
                           'phi_s': 'Soil moisture diffusion parameter',
                           'init_moist': 'Initial layer moisture content',
                           'elev': 'Average elevation of grid cell',
                           'depth': 'Thickness of each soil moisture layer',
                           'avg_T': 'Average soil temperature, used as the '
                                    'bottom boundary for soil heat flux '
                                    'solutions',
                           'dp': 'Soil thermal damping depth (depth at which '
                                 'soil temperature remains constant through '
                                 'the year, ~4 m)',
                           'bubble': 'Bubbling pressure of soil. Values should'
                                     ' be > 0.0',
                           'quartz': 'Quartz content of soil',
                           'bulk_density': 'Bulk density of soil layer',
                           'soil_density': 'Soil particle density, normally '
                                           '2685 kg/m3',
                           'off_gmt': 'Time zone offset from GMT. This '
                                      'parameter determines how VIC interprets'
                                      ' sub-daily time steps relative to the '
                                      'model start date and time.',
                           'Wcr_FRACT': 'Fractional soil moisture content at '
                                        'the critical point (~70%% of field '
                                        'capacity) (fraction of maximum '
                                        'moisture)',
                           'Wpwp_FRACT': 'Fractional soil moisture content at '
                                         'the wilting point (fraction of '
                                         'maximum moisture)',
                           'rough': 'Surface roughness of bare soil',
                           'snow_rough': 'Surface roughness of snowpack',
                           'annual_prec': 'Average annual precipitation.',
                           'resid_moist': 'Soil moisture layer residual '
                                          'moisture.',
                           'fs_active': 'If set to 1, then frozen soil '
                                        'algorithm is activated for the grid '
                                        'cell. A 0 indicates that frozen soils'
                                        ' are not computed even if soil '
                                        'temperatures fall below 0C.'}

        self.snow_param = {'cellnum': 'Grid cell number (should match numbers '
                                      'assigned in soil parameter file)',
                           'AreaFract': 'Fraction of grid cell covered by each'
                                        ' elevation band. Sum of the fractions'
                                        ' must equal 1.',
                           'elevation': 'Mean (or median) elevation of '
                                        'elevation band. This is used to '
                                        'compute the change in air temperature'
                                        ' from the grid cell mean elevation.',
                           'Pfactor': 'Fraction of cell precipitation that'
                                      'falls on each elevation band. Total '
                                      'must equal 1. To ignore effects of '
                                      'elevation on precipitation, set these '
                                      'fractions equal to the area fractions.'}

        self.veglib = {'Veg_class': 'Vegetation class identification number '
                                    '(reference index for library table)',
                       'lib_overstory': 'Flag to indicate whether or not the '
                                        'current vegetation type has an '
                                        'overstory (TRUE for overstory present'
                                        ' [e.g. trees], FALSE for overstory '
                                        'not present [e.g. grass])',
                       'lib_rarc': 'Architectural resistance of vegetation '
                                   'type (~2 s/m)',
                       'lib_rmin': 'Minimum stomatal resistance of vegetation '
                                   'type (~100 s/m)',
                       'lib_LAI': 'Leaf-area index of vegetation type',
                       'lib_albedo': 'Shortwave albedo for vegetation type',
                       'lib_rough': 'Vegetation roughness length (typically '
                                    '0.123 * vegetation height)',
                       'lib_displacement': 'Vegetation displacement height '
                                           '(typically 0.67 * vegetation '
                                           'height)',
                       'lib_wind_h': 'Height at which wind speed is measured.',
                       'lib_RGL': 'Minimum incoming shortwave radiation at '
                                  'which there will be transpiration. For '
                                  'trees this is about 30 W/m^2, for crops '
                                  'about 100 W/m^2.',
                       'lib_rad_atten': 'Radiation attenuation factor. '
                                        'Normally set to 0.5, though may need '
                                        'to be adjusted for high latitudes.',
                       'lib_wind_atten': 'Wind speed attenuation through the '
                                         'overstory. The default value has '
                                         'been 0.5.',
                       'lib_trunk_ratio': 'Ratio of total tree height that is '
                                          'trunk (no branches). The default '
                                          'value has been 0.2.',
                       'lib_comment': 'Comment block for vegetation type. '
                                      'Model skips end of line so spaces are '
                                      'valid entrys.'}

        self.veg_param = {'gridcell': 'Grid cell number',
                          'Nveg': 'Number of vegetation tiles in the grid '
                                  'cell',
                          'veg_class': 'Vegetation class identification number'
                                      ' (reference index to vegetation '
                                      'library)',
                          'Cv': 'Fraction of grid cell covered by vegetation '
                                'tile',
                          'root_depth': 'Root zone thickness (sum of depths is'
                                        ' total depth of root penetration)',
                          'root_fract': 'Fraction of root in the current root '
                                        'zone.',
                          'LAI': 'Leaf Area Index, one per month'}

        # if ORGANIC_FRACT:
        #     self.soil_param['organic'] = 'Fraction of soil layer that is organic'
        #     self.soil_param['bul_dens_org'] = 'Bulk density of organic portion of soil'
        #     self.soil_param['soil_dens_org'] = 'Soil particle density of organic portion of soil, normally 1300 kg/m3 '

        # if SPATIAL_FROST:
        #     self.soil_param['frost_slope'] = 'Slope of uniform distribution of soil temperature'

        # if SPATIAL_SNOW:
        #     self.soil_param['max_snow_distrib_slope'] = 'Maximum slope of the snow depth distribution.'

        # if EXCESS_ICE:
        #     self.soil_param['initial_ice_content'] = 'Initial volumetric ice content in soil '

        # if JULY_TAVG_SUPPLIED:
        #     self.soil_param['July_Tavg'] = 'Average July air temperature, used for treeline computations'

        # if BLOWING_SNOW:
        #     self.veg_param['sigma_slope'] = 'Standard deviation of terrain slopes within vegetation tile'
        #     self.veg_param['lag_one'] = 'Standard deviation of terrain slopes within vegetation tile'
        #     self.veg_param['fetch'] = 'Standard deviation of terrain slopes within vegetation tile'

        # if GLOBAL_LAI:
        #     print 'global lai included'
        #     self.veg_param['LAI'] = 'Leaf Area Index, one per month'
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
class units(object):
    def __init__(self, ORGANIC_FRACT=False, SPATIAL_FROST=False,
                 SPATIAL_SNOW=False, EXCESS_ICE=False,
                 JULY_TAVG_SUPPLIED=False, BLOWING_SNOW=False,
                 GLOBAL_LAI=False):

        self.soil_param = {'run_cell': 'N/A',
                           'gridcell': 'N/A',
                           'lats': 'degrees',
                           'lons': 'degrees',
                           'infilt': 'mm/day',
                           'Ds': 'fraction',
                           'Dsmax': 'mm/day',
                           'Ws': 'fraction',
                           'c': 'N/A',
                           'expt': 'N/A',
                           'Ksat': 'mm/day',
                           'phi_s': 'mm/mm',
                           'init_moist': 'mm',
                           'elev': 'm',
                           'depth': 'm',
                           'avg_T': 'C',
                           'dp': 'm',
                           'bubble': 'cm',
                           'quartz': 'fraction',
                           'bulk_density': 'kg/m3',
                           'soil_density': 'kg/m3',
                           'off_gmt': 'hours',
                           'Wcr_FRACT': 'fraction',
                           'Wpwp_FRACT': 'fraction',
                           'rough': 'm',
                           'snow_rough': 'm',
                           'annual_prec': 'mm',
                           'resid_moist': 'fraction',
                           'fs_active': 'binary'}

        self.snow_param = {'cellnum': 'N/A',
                           'AreaFract': 'fraction',
                           'elevation': 'm',
                           'Pfactor': 'fraction'}

        self.veglib = {'Veg_class': 'N/A',
                       'lib_overstory': 'N/A',
                       'lib_rarc': 's/m',
                       'lib_rmin': 's/m',
                       'lib_LAI': 'N/A',
                       'lib_albedo': 'fraction',
                       'lib_rough': 'm',
                       'lib_displacement': 'm',
                       'lib_wind_h': 'm',
                       'lib_RGL': 'W/m^2.',
                       'lib_rad_atten': 'fraction',
                       'lib_wind_atten': 'fraction',
                       'lib_trunk_ratio': 'fraction',
                       'lib_comment': 'N/A'}

        self.veg_param = {'gridcell': 'N/A',
                          'Nveg': 'N/A',
                          'veg_class': 'N/A',
                          'Cv': 'fraction',
                          'root_depth': 'm',
                          'root_fract': 'fraction',
                          'LAI': 'N/A'}
        # if ORGANIC_FRACT:
        #     self.soil_param['organic'] = 'fraction'
        #     self.soil_param['bul_dens_org'] = 'kg/m3'
        #     self.soil_param['soil_dens_org'] = 'kg/m3 '

        # if SPATIAL_FROST:
        #     self.soil_param['frost_slope'] = 'C'

        # if SPATIAL_SNOW:
        #     self.soil_param['max_snow_distrib_slope'] = 'm'

        # if EXCESS_ICE:
        #     self.soil_param['initial_ice_content'] = 'N/A'

        # if JULY_TAVG_SUPPLIED:
        #     self.soil_param['July_Tavg'] = 'C'

        # if BLOWING_SNOW:
        #     self.veg_param['sigma_slope'] = 'N/A'
        #     self.veg_param['lag_one'] = 'N/A'
        #     self.veg_param['fetch'] = 'm'

        #     if GLOBAL_LAI:
        #         self.veg_param['LAI'] = 'N/A'
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def main():
    """
    If called from the command line, main() will process command line arguments
    and run make_grid, Make grid reads the soil file first (mandatory) and then
    includes any other paramter files.  Parameter column locations and
    units/descriptions are from params.py.  Output is to a netcdf containing
    all the parameterfiles (default name is params.nc)
    """
    grid_file, soil_file, snow_file, veg_file, \
        vegl_file, out_file = process_command_line()
    grids = make_grid(grid_file, soil_file,
                      snow_file=snow_file,
                      veg_file=veg_file,
                      vegl_file=vegl_file,
                      nc_file=out_file)

    print('completed grid_parms.main(), output file was: {0}'.format(out_file))
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
#
def process_command_line():
    """
    Process command line arguments. Must have target grid (in netcdf format)
    and soil file (in standard vic format)
    Optionally may include snow and vegitation parameter files.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--grid_file",
                        type=str,
                        help="Input netCDF target grid",
                        default=None)
    parser.add_argument("-s", "--soil_file",
                        type=str,
                        help="Input file containing soil parameter data in "
                             "standard VIC format",
                        required=True)
    parser.add_argument("-e", "--snow_file",
                        type=str,
                        help="Input file containing snowband/elevation band "
                             "data in standard VIC format",
                        default=None)
    parser.add_argument("-v", "--veg_file",
                        type=str,
                        help="Input file containing vegitation parameter data "
                             "in standard VIC format",
                        default=None)
    parser.add_argument("-l", "--vegl_file",
                        type=str,
                        help="Input file containing vegitation library data in"
                             " standard VIC format",
                        default=None)
    parser.add_argument("-o", "--out_file",
                        type=str,
                        help="Input file containing vegitation parameter data "
                             "in standard VIC format",
                        default='params.nc')

    args = parser.parse_args()

    grid_file = args.grid_file
    soil_file = args.soil_file
    snow_file = args.snow_file
    veg_file = args.veg_file
    vegl_file = args.vegl_file
    out_file = args.out_file

    return grid_file, soil_file, snow_file, veg_file, vegl_file, out_file
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def make_grid(grid_file, soil_file, snow_file, veg_file, vegl_file,
              nc_file='params.nc'):
    """
    Make grid uses routines from params.py to read standard vic format
    parameter files.  After the parameter files are read, the files are placed
    onto the target grid using nearest neighbor mapping.  If a land mask is
    present in the target grid it will be used to exclude areas in the ocean.
    Finally, if the nc_file = 'any_string.nc', a netcdf file be written with
    the parameter data, if nc_file = False, the dictionary of grids is
    returned.
    """
    print('making grided parameters now...')

    soil_dict = soil(soil_file)

    if snow_file:
        snow_dict = snow(snow_file, soil_dict)
    else:
        snow_dict = False

    if veg_file:
        veg_dict = veg(veg_file, soil_dict, LAIindex=True)
    else:
        veg_dict = False

    if vegl_file:
        veglib_dict = veg_class(vegl_file)
    else:
        veglib_dict = False

    if grid_file:
        target_grid, target_attrs = read_netcdf(grid_file)
    else:
        target_grid, target_attrs = calc_grid(soil_dict['lats'],
                                              soil_dict['lons'])

    grid_dict = grid_params(soil_dict, target_grid, snow_dict=snow_dict,
                            veg_dict=veg_dict)

    if nc_file:
        write_netcdf(nc_file, target_attrs,
                     target_grid=target_grid,
                     soil_grid=grid_dict['soil_dict'],
                     snow_grid=grid_dict['snow_dict'],
                     veglib_dict=veglib_dict,
                     veg_grid=grid_dict['veg_dict'])

    return grid_dict
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def calc_grid(lats, lons, decimals=4):
    """ determine shape of regular grid from lons and lats"""

    print('Calculating grid size now...')

    target_grid = {}

    # get unique lats and lons
    ulons = np.sort(np.unique(lons.round(decimals=decimals)))
    print('found {0} unique lons'.format(len(ulons)))
    lon_step, lon_count = stats.mode(np.diff(ulons))

    ulats = np.sort(np.unique(lats.round(decimals=decimals)))
    print('found {0} unique lats'.format(len(ulats)))
    lat_step, lat_count = stats.mode(np.diff(ulats))

    # check that counts and steps make sense
    if lon_step != lat_step:
        print('WARNING: lon_step ({0}) and lat_step ({1}) do not '
              'match'.format(lon_step, lat_step))
    if lat_count/len(ulats) < 0.95:
        print('WARNING: lat_count of mode is less than 95% ({0}%) of'
              ' len(lats)'.format(lat_count/len(ulats)))
    if lon_count/len(ulons) < 0.95:
        print('WARNING: lon_count of mode is less than 95% ({0}%) of'
              ' len(lons)'.format(lon_count/len(ulons)))

    if lats.min() < -55 and lats.max() > 70:
        # assume global grid
        print('assuming grid is meant to be global...')
        target_grid[XVAR] = np.linspace(-180+lon_step[0]/2,
                                        180-lon_step[0]/2,
                                        360/lon_step[0])
        target_grid[YVAR] = np.linspace(-90+lat_step[0]/2,
                                        90-lat_step[0]/2,
                                        180/lat_step[0])
    else:
        target_grid[XVAR] = np.arange(lons.min(),
                                      lons.max()+lon_step[0],
                                      lon_step[0])
        target_grid[YVAR] = np.arange(lats.min(),
                                      lats.max()+lat_step[0],
                                      lat_step[0])

    y, x = latlon2yx(lats, lons, target_grid[YVAR],
                     target_grid[XVAR])

    mask = np.zeros((len(target_grid[YVAR]),
                     len(target_grid[XVAR])), dtype=int)

    mask[y, x] = 1

    target_grid['mask'] = mask

    target_attrs = {YVAR: {'long_name': 'latitude coordinate',
                           'units': 'degrees_north'},
                    XVAR: {'long_name': 'longitude coordinate',
                           'units': 'degrees_east'},
                    'mask': {'long_name': 'domain mask',
                             'comment': '0 indicates grid cell is not active'}
                    }

    print('Created a target grid based on the lats '
          'and lon in the soil parameter file')
    print('Grid Size: {0}'.format(mask.shape))

    return target_grid, target_attrs
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# find x y coordinates
def latlon2yx(plats, plons, glats, glons):
    """find y x coordinates """
    # use astronomical conventions for longitude
    # (i.e. negative longitudes to the east of 0)
    if (glons.max() > 180):
        posinds = np.nonzero(glons > 180)
        glons[posinds] -= 360
        print('adjusted grid lon minimum ')

    if (plons.max() > 180):
        posinds = np.nonzero(plons > 180)
        plons[posinds] -= 360
        print('adjusted points lon minimum')

    if glons.ndim == 1 or glats.ndim == 1:
        print('creating 2d coordinate arrays')
        glats, glons = np.meshgrid(glats, glons, indexing='ij')

    combined = np.dstack(([glats.ravel(), glons.ravel()]))[0]
    points = list(np.vstack((np.array(plats), np.array(plons))).transpose())

    mytree = cKDTree(combined)
    dist, indexes = mytree.query(points, k=1)
    y, x = np.unravel_index(indexes, glons.shape)

    return y, x
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def grid_params(soil_dict, target_grid, snow_dict, veg_dict):
    """
    Reads the coordinate information from the soil_dict and target_grid and
    maps all input dictionaries to the target grid.  Returns a grid_dict with
    the mapped input dictionary data.
    """
    print('gridding params now...')

    yi, xi = latlon2yx(soil_dict['lats'], soil_dict['lons'],
                       target_grid[YVAR], target_grid[XVAR])

    in_dicts = {'soil_dict': soil_dict}
    out_dicts = OrderedDict()
    if snow_dict:
        in_dicts['snow_dict'] = snow_dict
    else:
        out_dicts['snow_dict'] = False
    if veg_dict:
        in_dicts['veg_dict'] = veg_dict
    else:
        out_dicts['veg_dict'] = False

    # get "unmasked" mask
    mask = target_grid['mask'].data

    ysize, xsize = target_grid['mask'].shape

    ymask, xmask = np.nonzero(mask != 1)
    print(ymask, xmask)

    for name, mydict in in_dicts.iteritems():
        out_dict = OrderedDict()

        for var in mydict:
            if mydict[var].dtype in [np.int, np.int64, np.int32]:
                fill_val = FILLVALUE_I
                dtype = np.int
            else:
                fill_val = FILLVALUE_F
                dtype = np.float

            if mydict[var].ndim == 1:
                out_dict[var] = np.ma.zeros((ysize, xsize),
                                            dtype=dtype) + fill_val
                out_dict[var][yi, xi] = mydict[var]
                out_dict[var][ymask, xmask] = fill_val

            elif mydict[var].ndim == 2:
                steps = mydict[var].shape[1]
                out_dict[var] = np.ma.zeros((steps, ysize, xsize),
                                            dtype=dtype) + fill_val
                for i in xrange(steps):
                    out_dict[var][i, yi, xi] = mydict[var][:, i]
                out_dict[var][:, ymask, xmask] = fill_val

            elif mydict[var].ndim == 3:
                j = mydict[var].shape[1]
                k = mydict[var].shape[2]
                out_dict[var] = np.ma.zeros((j, k, ysize, xsize),
                                            dtype=dtype) + fill_val
                for jj in xrange(j):
                    for kk in xrange(k):
                        out_dict[var][jj, kk, yi, xi] = mydict[var][:, jj, kk]
                for y, x in zip(ymask, xmask):
                    out_dict[var][:, :, y, x] = fill_val

            out_dict[var] = np.ma.masked_values(out_dict[var], fill_val)

        out_dicts[name] = out_dict

    return out_dicts
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
#  Write output to netCDF
def write_netcdf(myfile, target_attrs, target_grid,
                 soil_grid=None, snow_grid=None, veg_grid=None,
                 veglib_dict=None):
    """
    Write the gridded parameters to a netcdf4 file
    Will only write paramters that it is given
    Reads attributes from params.py and from targetAtters dictionary read from
    grid_file
    """
    f = Dataset(myfile, 'w', format='NETCDF4')

    # write attributes for netcdf
    f.description = 'VIC parameter file'
    f.history = 'Created: {0}\n'.format(tm.ctime(tm.time()))
    f.history += ' '.join(sys.argv) + '\n'
    f.source = sys.argv[0]  # prints the name of script used
    f.username = getuser()
    f.host = socket.gethostname()

    unit = units(GLOBAL_LAI=True)
    desc = description()

    # target grid
    # coordinates
    if target_grid[XVAR].ndim == 1:
        x = f.createDimension('lon', len(target_grid[XVAR]))
        y = f.createDimension('lat', len(target_grid[YVAR]))
        dims2 = ('lat', 'lon', )
        coordinates = None

        v = f.createVariable('lat', NC_DOUBLE, ('lat',))
        v[:] = target_grid[YVAR]
        v.units = 'degrees_north'
        v.long_name = "latitude of grid cell center"

        v = f.createVariable('lon', NC_DOUBLE, ('lon',))
        v[:] = target_grid[XVAR]
        v.units = 'degrees_east'
        v.long_name = "longitude of grid cell center"

    else:
        x = f.createDimension('nj', target_grid[XVAR].shape[0])
        y = f.createDimension('ni', target_grid[YVAR].shape[1])
        dims2 = ('nj', 'ni', )
        coordinates = "{0} {1}".format(XVAR, YVAR)

        v = f.createVariable(YVAR, NC_DOUBLE, dims2)
        v[:, :] = target_grid[YVAR]
        v.units = 'degrees_north'

        v = f.createVariable(XVAR, NC_DOUBLE, dims2)
        v[:, :] = target_grid[XVAR]
        v.units = 'degrees_east'

        # corners
        if ('xv' in target_grid) and ('yv' in target_grid):
            nv4 = f.createDimension('nv4', 4)
            dims_corner = ('nv4', 'nj', 'ni', )

            v = f.createVariable('xv', NC_DOUBLE, dims_corner)
            v[:, :, :] = target_grid['xv']

            v = f.createVariable('yv', NC_DOUBLE, dims_corner)
            v[:, :, :] = target_grid['yv']

    # mask
    v = f.createVariable('mask', NC_DOUBLE, dims2)
    v[:, :] = target_grid['mask']
    v.long_name = 'land mask'
    if coordinates:
        v.coordinates = coordinates

    # set attributes
    for var in target_grid:
        for name, attr in target_attrs[var].iteritems():
            try:
                setattr(v, name, attr)
            except:
                print('dont have units or description for {0}'.format(var))

    # Layers
    layers = f.createDimension('nlayer', soil_grid['soil_density'].shape[0])
    layer_dims = ('nlayer', ) + dims2

    # soil grid
    for var, data in soil_grid.iteritems():
        print('writing var: {0}'.format(var))

        if data.ndim == 1:
            v = f.createVariable(var, NC_DOUBLE, ('nlayer', ),
                                 fill_value=FILLVALUE_F)
            v[:] = data

        elif data.ndim == 2:
            v = f.createVariable(var, NC_DOUBLE, dims2, fill_value=FILLVALUE_F)
            v[:, :] = data

        elif data.ndim == 3:
            v = f.createVariable(var, NC_DOUBLE, layer_dims,
                                 fill_value=FILLVALUE_F)
            v[:, :, :] = data
        else:
            raise IOError('all soil vars should be 2 or 3 dimensions')

        # add attributes
        v.units = unit.soil_param[var]
        v.description = desc.soil_param[var]
        v.long_name = var
        if coordinates:
            v.coordinates = coordinates

    if snow_grid:
        try:
            del snow_grid['gridcell']
        except:
            pass

        snowbnds = f.createDimension('snow_band',
                                     snow_grid['AreaFract'].shape[0])
        snow_dims = ('snow_band', ) + dims2

        for var, data in snow_grid.iteritems():
            print('writing var: {0}'.format(var))

            if data.ndim == 2:
                v = f.createVariable(var, NC_DOUBLE, dims2,
                                     fill_value=FILLVALUE_F)
                v[:, :] = data
            elif data.ndim == 3:
                v = f.createVariable(var, NC_DOUBLE, snow_dims,
                                     fill_value=FILLVALUE_F)
                v[:, :, :] = data
            else:
                raise IOError('all snow vars should be 2 or 3 dimensions')

            v.units = unit.snow_param[var]
            v.description = desc.snow_param[var]
            if coordinates:
                v.coordinates = coordinates

    if veg_grid:
        try:
            del veg_grid['gridcell']
        except:
            pass

        veg_type = f.createDimension('veg_class', veg_grid['Cv'].shape[0])
        root_zone = f.createDimension('root_zone', 3)
        month = f.createDimension('month', 12)

        v = f.createVariable('month', NC_INT, ('month', ))
        v[:] = np.arange(1, 13)
        v.long_name = 'month of year'

        for var, data in veg_grid.iteritems():
            print('writing var: {0}'.format(var))

            if veg_grid[var].ndim == 2:
                v = f.createVariable(var, NC_DOUBLE, dims2,
                                     fill_value=FILLVALUE_F)
                v[:, :] = data

            elif veg_grid[var].ndim == 3:
                mycoords = ('veg_class', ) + dims2
                v = f.createVariable(var, NC_DOUBLE, mycoords,
                                     fill_value=FILLVALUE_F)
                v[:, :, :] = data

            elif var == 'LAI':
                mycoords = ('veg_class', 'month', ) + dims2
                v = f.createVariable(var, NC_DOUBLE, mycoords,
                                     fill_value=FILLVALUE_F)
                v[:, :, :, :] = data

            elif veg_grid[var].ndim == 4:
                mycoords = ('veg_class', 'root_zone', ) + dims2
                v = f.createVariable(var, NC_DOUBLE, mycoords,
                                     fill_value=FILLVALUE_F)
                v[:, :, :, :] = data

            else:
                raise ValueError('only able to handle dimensions <=4')

            v.units = unit.veg_param[var]
            v.description = desc.veg_param[var]
            v.long_name = var
            if coordinates:
                v.coordinates = coordinates

        if veglib_dict:
            print('writing var: {0}'.format(var))

            for var, data in veglib_dict.iteritems():
                if data.ndim == 1:
                    v = f.createVariable(var, NC_DOUBLE, ('veg_class', ),
                                         fill_value=FILLVALUE_F)
                    v[:] = data
                elif data.ndim == 2:
                    v = f.createVariable(var, NC_DOUBLE,
                                         ('veg_class', 'month', ),
                                         fill_value=FILLVALUE_F)
                    v[:, :] = data
                else:
                    raise IOError('veglib_dict shouldnt have data with more \
                                   that 2 dimentions')

                v.units = unit.veglib[var]
                v.description = desc.veglib[var]
                v.long_name = var

    f.close()

    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def soil(in_file, nlayers=3):
    """
    Load the entire soil file into a dictionary of numpy arrays.
    Also reorders data to match gridcell order of soil file.
    """
    print('reading {0}'.format(in_file))
    data = np.loadtxt(in_file)

    c = cols(nlayers=nlayers)

    soil_dict = OrderedDict()
    for var, columns in c.soil_param.iteritems():
        soil_dict[var] = np.squeeze(data[:, columns])

    return soil_dict
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def snow(snow_file, soil_dict, snow_bands=5):
    """
    Load the entire snow file into a dictionary of numpy arrays.
    Also reorders data to match gridcell order of soil file.
    """

    print('reading {0}'.format(snow_file))

    data = np.loadtxt(snow_file)

    c = cols(snow_bands=snow_bands)

    snow_dict = OrderedDict()
    for var in c.snow_param:
        snow_dict[var] = data[:, c.snow_param[var]]

    target = soil_dict['gridcell'].argsort()
    indexes = target[np.searchsorted(soil_dict['gridcell'][target],
                                     snow_dict['cellnum'])]

    for var in snow_dict:
        snow_dict[var] = np.squeeze(snow_dict[var][indexes])

    return snow_dict
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def veg(veg_file, soil_dict, maxRoots=3, vegClasses=11,
        cells=False, BLOWING_SNOW=False, LAIindex=False):
    """
    Read the vegetation file from vegFile.  Assumes max length for rootzones
    and vegclasses.  Also reorders data to match gridcell order of soil file.
    """

    print('reading {0}'.format(veg_file))

    with open(veg_file) as f:
        lines = f.readlines()

    if not cells:
        cells = len(lines)

    gridcel = np.zeros(cells)
    Nveg = np.zeros(cells)
    Cv = np.zeros((cells, vegClasses))
    root_depth = np.zeros((cells, vegClasses, maxRoots))
    root_fract = np.zeros((cells, vegClasses, maxRoots))
    if BLOWING_SNOW:
        sigma_slope = np.zeros((cells, vegClasses))
        lag_one = np.zeros((cells, vegClasses))
        fetch = np.zeros((cells, vegClasses))
    if LAIindex:
        lfactor = 2
        LAI = np.zeros((cells, vegClasses, 12))
    else:
        lfactor = 1

    row = 0
    cell = 0

    while row < len(lines):
        line = lines[row].strip('\n').split(' ')
        gridcel[cell], Nveg[cell] = np.array(line).astype(int)
        numrows = Nveg[cell]*lfactor + row
        row += 1

        while row < numrows:
            lines[row] = lines[row].strip()
            line = lines[row].strip('\n').split(' ')
            temp = np.array(line).astype(float)
            vind = int(temp[0])-1
            Cv[cell, vind] = temp[1]

            if not BLOWING_SNOW:
                rind = (len(temp)-2)/2

            else:
                rind = (len(temp)-5)/2
                sigma_slope[cell, vind, :rind] = temp[-3]
                lag_one[cell, vind, :rind] = temp[-2]
                fetch[cell, vind, :rind] = temp[-1]

            root_depth[cell, vind, :rind] = temp[2::2]
            root_fract[cell, vind, :rind] = temp[3::2]
            row += 1
            if LAIindex:
                lines[row] = lines[row].strip()
                line = lines[row].strip('\n').split(' ')
                LAI[cell, vind, :] = np.array(line).astype(float)
                row += 1
        cell += 1
    veg_dict = OrderedDict()
    veg_dict['gridcell'] = gridcel[:cell]
    veg_dict['Nveg'] = Nveg[:cell]
    veg_dict['Cv'] = Cv[:cell, :]
    veg_dict['root_depth'] = root_depth[:cell, :]
    veg_dict['root_fract'] = root_fract[:cell, :]

    if BLOWING_SNOW:
        veg_dict['sigma_slope'] = sigma_slope[:cell, :]
        veg_dict['lag_one'] = lag_one[:cell, :]
        veg_dict['fetch'] = fetch[:cell, :]

    if LAIindex:
        veg_dict['LAI'] = LAI[:cell, :, :]

    inds = []
    for sn in soil_dict['gridcell']:
        inds.append(np.nonzero(veg_dict['gridcell'] == sn))

    for var in veg_dict:
        veg_dict[var] = np.squeeze(veg_dict[var][inds])
    return veg_dict


# -------------------------------------------------------------------- #
def veg_class(veg_file, maxcols=57, skiprows=1):
    """
    Load the entire vegetation library file into a dictionary of numpy arrays.
    Also reorders data to match gridcell order of soil file.
    """

    print('reading {0}'.format(veg_file))

    usecols = np.arange(maxcols)

    data = np.loadtxt(veg_file, usecols=usecols, skiprows=skiprows)

    c = cols()

    veglib_dict = OrderedDict()
    for var in c.veglib:
        veglib_dict[var] = np.squeeze(data[:, c.veglib[var]])

    return veglib_dict
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
if __name__ == "__main__":
    main()
# -------------------------------------------------------------------- #
