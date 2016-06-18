#!/usr/bin/env python
"""
grid_params.py

A toolkit for converting classic vic parameters to netcdf format
"""

from __future__ import print_function
import sys
import numpy as np
from netCDF4 import Dataset, default_fillvals
from scipy.spatial import cKDTree
from scipy import stats
import time as tm
import socket
from getpass import getuser
from collections import OrderedDict
from warnings import warn
from tonic.io import read_netcdf
from tonic.pycompat import pyrange, pyzip
import re


# -------------------------------------------------------------------- #
description = 'Converter for VIC ASCII style parameters to gridded netCDF'
help = 'Converter for VIC ASCII style parameters to gridded netCDF'
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# precision
PRECISION = 1.0e-30
NC_DOUBLE = 'f8'
NC_FLOAT = 'f4'
NC_INT = 'i4'
NC_CHAR = 'S1'
MAX_NC_CHARS = 256

# Months per year
MONTHS_PER_YEAR = 12

# Default values of veg params over bare soil
bare_vegparam = {'overstory': 0,
                 'rarc': 100,
                 'rmin': 0,
                 'LAI': 0,
                 'fcanopy': 0,
                 'albedo': 0.2,
                 'veg_rough': 0,
                 'displacement': 0,
                 'wind_h': 2,
                 'RGL': 0,
                 'rad_atten': 0,
                 'wind_atten': 0,
                 'trunk_ratio': 0,
                 'comment': 'Bare Soil',
                 'root_depth': 0,
                 'root_fract': 0,
                 'sigma_slope': 0.08,
                 'lag_one': 0.8,
                 'fetch': 1000.0,
                 'Ctype': 0,
                 'MaxCarboxRate': 0,
                 'MaxE_or_CO2Spec': 0,
                 'CO2Specificity': 0,
                 'LUE': 0,
                 'Nscale': 0,
                 'Wnpp_inhib': 0,
                 'NPPfactor_sat': 0}

# fill values
FILLVALUE_F = default_fillvals[NC_DOUBLE]
FILLVALUE_I = default_fillvals[NC_INT]

XVAR = 'xc'
YVAR = 'yc'

# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
class Cols(object):
    def __init__(self, nlayers=3, snow_bands=5, organic_fract=False,
                 spatial_frost=False, spatial_snow=False,
                 july_tavg_supplied=False, veglib_fcan=False,
                 veglib_photo=False):

        # Soil Parameters
        self.soil_param = OrderedDict([('run_cell', np.array([0])),
                                       ('gridcell', np.array([1])),
                                       ('lats', np.array([2])),
                                       ('lons', np.array([3])),
                                       ('infilt', np.array([4])),
                                       ('Ds', np.array([5])),
                                       ('Dsmax', np.array([6])),
                                       ('Ws', np.array([7])),
                                       ('c', np.array([8]))])

        i = 9
        for var in ['expt', 'Ksat', 'phi_s', 'init_moist']:
            self.soil_param[var] = np.arange(i, i + nlayers)
            i += nlayers

        self.soil_param['elev'] = np.array([i])
        i += 1
        self.soil_param['depth'] = np.arange(i, i + nlayers)
        i += nlayers
        self.soil_param['avg_T'] = np.array([i])
        i += 1
        self.soil_param['dp'] = np.array([i])
        i += 1

        varnames = ['bubble', 'quartz', 'bulk_density', 'soil_density']
        if organic_fract:
            varnames.append(['organic', 'bulk_dens_org', 'soil_dens_org'])
        for var in varnames:
            self.soil_param[var] = np.arange(i, i + nlayers)
            i += nlayers

        self.soil_param['off_gmt'] = np.array([i])
        i += 1

        for var in ['Wcr_FRACT', 'Wpwp_FRACT']:
            self.soil_param[var] = np.arange(i, i + nlayers)
            i += nlayers

        for var in ['rough', 'snow_rough', 'annual_prec']:
            self.soil_param[var] = np.array([i])
            i += 1

        self.soil_param['resid_moist'] = np.arange(i, i + nlayers)
        i += nlayers

        self.soil_param['fs_active'] = np.array([i])
        i += 1

        if spatial_frost:
            self.soil_param['frost_slope'] = np.array([i])
            i += 1

        if spatial_snow:
            self.soil_param['max_snow_distrib_slope'] = np.array([i])
            i += 1

        if july_tavg_supplied:
            self.soil_param['July_Tavg'] = np.array([i])
            i += 1

        # Snow Parameters
        self.snow_param = OrderedDict([('cellnum', np.array([0]))])
        i = 1
        for var in ['AreaFract', 'elevation', 'Pfactor']:
            self.snow_param[var] = np.arange(i, i + snow_bands)
            i += snow_bands

        # Veg Library
        self.veglib = OrderedDict([('Veg_class', np.array([0])),
                                   ('lib_overstory', np.array([1])),
                                   ('lib_rarc', np.array([2])),
                                   ('lib_rmin', np.array([3])),
                                   ('lib_LAI', np.arange(4, 16))])

        varnames = ['lib_albedo', 'lib_veg_rough', 'lib_displacement']
        if veglib_fcan:
            varnames = ['lib_fcanopy'] + varnames
        i = 16
        for var in varnames:
            self.veglib[var] = np.arange(i, i + MONTHS_PER_YEAR)
            i += MONTHS_PER_YEAR

        for var in ['lib_wind_h', 'lib_RGL', 'lib_rad_atten',
                    'lib_wind_atten', 'lib_trunk_ratio']:
            self.veglib[var] = np.array([i])
            i += 1

        if veglib_photo:
            varnames = ['lib_Ctype', 'lib_MaxCarboxRate',
                        'lib_MaxE_or_CO2Spec', 'lib_LUE',
                        'lib_Nscale', 'lib_Wnpp_inhib', 'lib_NPPfactor_sat']
            for var in varnames:
                self.veglib[var] = np.array([i])
                i += 1

        self.veglib['lib_comment'] = np.array([i])


# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
class Format(object):
    def __init__(self, nlayers=3, snow_bands=5, organic_fract=False,
                 spatial_frost=False, spatial_snow=False,
                 july_tavg_supplied=False, veglib_fcan=False,
                 veglib_photo=False, blowing_snow=False,
                 vegparam_lai=False, vegparam_fcan=False,
                 vegparam_albedo=False, lakes=False):

        # Soil Params
        self.soil_param = {'run_cell': '%1i',
                           'gridcell': '%1i',
                           'lats': '%12.7g',
                           'lons': '%12.7g',
                           'infilt': '%12.7g',
                           'Ds': '%12.7g',
                           'Dsmax': '%12.7g',
                           'Ws': '%12.7g',
                           'c': '%12.7g',
                           'expt': '%12.7g',
                           'Ksat': '%12.7g',
                           'phi_s': '%12.7g',
                           'init_moist': '%12.7g',
                           'elev': '%12.7g',
                           'depth': '%12.7g',
                           'avg_T': '%12.7g',
                           'dp': '%12.7g',
                           'bubble': '%12.7g',
                           'quartz': '%12.7g',
                           'bulk_density': '%12.7g',
                           'soil_density': '%12.7g',
                           'off_gmt': '%12.7g',
                           'Wcr_FRACT': '%12.7g',
                           'Wpwp_FRACT': '%12.7g',
                           'rough': '%12.7g',
                           'snow_rough': '%12.7g',
                           'annual_prec': '%12.7g',
                           'resid_moist': '%12.7g',
                           'fs_active': '%1i'}
        if organic_fract:
            self.soil_param['organic'] = '%12.7g'
            self.soil_param['bulk_dens_org'] = '%12.7g'
            self.soil_param['soil_dens_org'] = '%12.7g'
        if spatial_frost:
            self.soil_param['frost_slope'] = '%12.7g'
        if spatial_snow:
            self.soil_param['max_snow_distrib_slope'] = '%12.7g'
        if july_tavg_supplied:
            self.soil_param['July_Tavg'] = '%12.7g'

        # Snow Band Params
        self.snow_param = {'cellnum': '%1i',
                           'AreaFract': '%12.7g',
                           'elevation': '%12.7g',
                           'Pfactor': '%12.7g'}

        # Veg Library
        self.veglib = {'Veg_class': '%1i',
                       'lib_overstory': '%1i',
                       'lib_rarc': '%12.7g',
                       'lib_rmin': '%12.7g',
                       'lib_LAI': '%12.7g',
                       'lib_albedo': '%12.7g',
                       'lib_veg_rough': '%12.7g',
                       'lib_displacement': '%12.7g',
                       'lib_wind_h': '%12.7g',
                       'lib_RGL': '%12.7g',
                       'lib_rad_atten': '%12.7g',
                       'lib_wind_atten': '%12.7g',
                       'lib_trunk_ratio': '%12.7g',
                       'lib_comment': '%s'}
        if veglib_fcan:
            self.veglib['lib_fcanopy'] = '%12.7g'
        if veglib_photo:
            self.veglib['lib_Ctype'] = '%1i'
            self.veglib['lib_MaxCarboxRate'] = '%12.7g'
            self.veglib['lib_MaxE_or_CO2Spec'] = '%12.7g'
            self.veglib['lib_LUE'] = '%12.7g'
            self.veglib['lib_Nscale'] = '%1i'
            self.veglib['lib_Wnpp_inhib'] = '%12.7g'
            self.veglib['lib_NPPfactor_sat'] = '%12.7g'

        # Veg Params
        self.veg_param = {'gridcell': '%1i',
                          'Nveg': '%1i',
                          'veg_class': '%1i',
                          'Cv': '%12.7g',
                          'root_depth': '%12.7g',
                          'root_fract': '%12.7g'}
        if blowing_snow:
            self.veg_param['sigma_slope'] = '%12.7g'
            self.veg_param['lag_one'] = '%12.7g'
            self.veg_param['fetch'] = '%12.7g'
        if vegparam_lai:
            self.veg_param['LAI'] = '%12.7g'
        if vegparam_fcan:
            self.veg_param['fcanopy'] = '%12.7g'
        if vegparam_albedo:
            self.veg_param['albedo'] = '%12.7g'

        # Lake Params
        if lakes:
            self.lake_param = {'gridcell': '%1i',
                               'lake_idx': '%1i',
                               'numnod': '%1i',
                               'mindepth': '%12.7g',
                               'wfrac': '%12.7g',
                               'depth_in': '%12.7g',
                               'rpercent': '%12.7g',
                               'basin_depth': '%12.7g',
                               'basin_area': '%12.7g'}

# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
class Desc(object):
    def __init__(self, organic_fract=False, spatial_frost=False,
                 spatial_snow=False, july_tavg_supplied=False,
                 veglib_fcan=False, veglib_photo=False,
                 blowing_snow=False, vegparam_lai=False,
                 vegparam_fcan=False, vegparam_albedo=False,
                 lakes=False):

        # Soil Params
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
                                        'cell. A 0 indicates that frozen '
                                        'soils are not computed even if soil '
                                        'temperatures fall below 0C.'}
        if organic_fract:
            self.soil_param['organic'] = 'Organic content of soil'
            self.soil_param['bulk_dens_org'] = 'Bulk density of organic portion of soil'
            self.soil_param['soil_dens_org'] = 'Soil density of organic portion of soil'
        if spatial_frost:
            self.soil_param['frost_slope'] = 'Slope of uniform distribution of soil temperature'
        if spatial_snow:
            self.soil_param['max_snow_distrib_slope'] = 'Maximum slope of the snow depth distribution'
        if july_tavg_supplied:
            self.soil_param['July_Tavg'] = 'Average July air temperature'

        # Snow Band Params
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

        # Veg Library
        self.veglib = {'Veg_class': 'Vegetation class identification number '
                                    '(reference index for library table)',
                       'lib_overstory': 'Flag to indicate whether or not the '
                                        'current vegetation type has an '
                                        'overstory (1 for overstory present'
                                        ' [e.g. trees], 0 for overstory '
                                        'not present [e.g. grass])',
                       'lib_rarc': 'Architectural resistance of vegetation '
                                   'type (~2 s/m)',
                       'lib_rmin': 'Minimum stomatal resistance of vegetation '
                                   'type (~100 s/m)',
                       'lib_LAI': 'Leaf-area index of vegetation type',
                       'lib_albedo': 'Shortwave albedo for vegetation type',
                       'lib_veg_rough': 'Vegetation roughness length ('
                                        'typically 0.123 * vegetation height)',
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
        if veglib_fcan:
            self.veglib['lib_fcanopy'] = 'Canopy cover fraction, one per month'
        if veglib_photo:
            self.veglib['lib_Ctype'] = 'Photosynthetic pathway (0 = C3; 1 = C4)'
            self.veglib['lib_MaxCarboxRate'] = 'Maximum carboxylation rate at 25 C'
            self.veglib['lib_MaxE_or_CO2Spec'] = 'Maximum electron transport rate at 25 C (C3) or CO2 specificity (C4)'
            self.veglib['lib_LUE'] = 'Light use efficiency'
            self.veglib['lib_Nscale'] = '1 = this class employs nitrogen scaling factors; 0 = no nitrogen scaling factors'
            self.veglib['lib_Wnpp_inhib'] = 'Fraction of maximum moisture storage in top soil layer above which photosynthesis begins to be inhibited by wet conditions'
            self.veglib['lib_NPPfactor_sat'] = 'NPP inhibition factor under saturated conditions (when moisture = 100% of maximum)'

        # Veg Params
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
                                        'zone'}
        if blowing_snow:
            self.veg_param['sigma_slope'] = 'Std. deviation of terrain slope within veg tile'
            self.veg_param['lag_one'] = 'Lag one gradient autocorrelation of terrain slope'
            self.veg_param['fetch'] = 'Average fetch length within veg tile'
        if vegparam_lai:
            self.veg_param['LAI'] = 'Leaf Area Index, one per month'
        if vegparam_fcan:
            self.veg_param['fcanopy'] = 'Canopy cover fraction, one per month'
        if vegparam_albedo:
            self.veg_param['albedo'] = 'Albedo, one per month'

        # Lake Params
        if lakes:
            self.lake_param = {'gridcell': 'Grid cell number',
                               'lake_idx': 'index of veg tile that contains '
                               'the lake/wetland',
                               'numnod': 'Maxium number of lake layers in the '
                                         'grid cell',
                               'mindepth': 'Minimum lake water depth for '
                                           'channel runoff to occur',
                               'wfrac': 'Channel outlet width (expressed as a '
                                        'fraction of lake perimeter',
                               'depth_in': 'Initial lake depth',
                               'rpercent': 'Fraction of runoff from other veg '
                                           'tiles that flows into the lake',
                               'basin_depth': 'Elevation (above lake bottom) '
                                              'of points on lake depth-area '
                                              'curve',
                               'basin_area': 'Surface area (expressed as '
                                             'fraction of grid cell area) of '
                                             'points on lake depth-area curve'}
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
class Units(object):
    def __init__(self, organic_fract=False, spatial_frost=False,
                 spatial_snow=False, july_tavg_supplied=False,
                 veglib_fcan=False, veglib_photo=False,
                 blowing_snow=False, vegparam_lai=False,
                 vegparam_fcan=False, vegparam_albedo=False,
                 lakes=False):

        # Soil Params
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
        if organic_fract:
            self.soil_param['organic'] = 'fraction'
            self.soil_param['bulk_dens_org'] = 'kg/m3'
            self.soil_param['soil_dens_org'] = 'kg/m3'
        if spatial_frost:
            self.soil_param['frost_slope'] = 'C'
        if spatial_snow:
            self.soil_param['max_snow_distrib_slope'] = 'm'
        if july_tavg_supplied:
            self.soil_param['July_Tavg'] = 'C'

        # Snow Band Params
        self.snow_param = {'cellnum': 'N/A',
                           'AreaFract': 'fraction',
                           'elevation': 'm',
                           'Pfactor': 'fraction'}

        # Veg Library
        self.veglib = {'Veg_class': 'N/A',
                       'lib_overstory': 'N/A',
                       'lib_rarc': 's/m',
                       'lib_rmin': 's/m',
                       'lib_LAI': 'm2/m2',
                       'lib_albedo': 'fraction',
                       'lib_veg_rough': 'm',
                       'lib_displacement': 'm',
                       'lib_wind_h': 'm',
                       'lib_RGL': 'W/m^2.',
                       'lib_rad_atten': 'fraction',
                       'lib_wind_atten': 'fraction',
                       'lib_trunk_ratio': 'fraction',
                       'lib_comment': 'N/A'}
        if veglib_fcan:
            self.veglib['lib_fcanopy'] = 'fraction'
        if veglib_photo:
            self.veglib['lib_Ctype'] = '0 or 1'
            self.veglib['lib_MaxCarboxRate'] = 'mol CO2/m2s'
            self.veglib['lib_MaxE_or_CO2Spec'] = 'mol CO2/m2s'
            self.veglib['lib_LUE'] = 'mol CO2/mol photons'
            self.veglib['lib_Nscale'] = '0 or 1'
            self.veglib['lib_Wnpp_inhib'] = 'fraction'
            self.veglib['lib_NPPfactor_sat'] = 'fraction'

        # Veg Params
        self.veg_param = {'gridcell': 'N/A',
                          'Nveg': 'N/A',
                          'veg_class': 'N/A',
                          'Cv': 'fraction',
                          'root_depth': 'm',
                          'root_fract': 'fraction'}
        if blowing_snow:
            self.veg_param['sigma_slope'] = 'm'
            self.veg_param['lag_one'] = 'fraction'
            self.veg_param['fetch'] = 'm'
        if vegparam_lai:
            self.veg_param['LAI'] = 'm2/m2'
        if vegparam_fcan:
            self.veg_param['fcanopy'] = 'fraction'
        if vegparam_albedo:
            self.veg_param['albedo'] = 'fraction'

        # Lake Params
        if lakes:
            self.lake_param = {'gridcell': 'N/A',
                               'lake_idx': 'N/A',
                               'numnod': 'N/A',
                               'mindepth': 'm',
                               'wfrac': 'fraction',
                               'depth_in': 'm',
                               'rpercent': 'fraction',
                               'basin_depth': 'm',
                               'basin_area': 'fraction'}

# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def _run(args):
    """
    """
    nc_file = make_grid(grid_file=args.grid_file,
                        soil_file=args.soil_file,
                        snow_file=args.snow_file,
                        vegl_file=args.vegl_file,
                        veg_file=args.veg_file,
                        lake_file=args.lake_file,
                        nc_file=args.out_file,
                        version_in=args.VIC_version_in,
                        grid_decimal=args.grid_decimal,
                        nlayers=args.nlayers,
                        snow_bands=args.snow_bands,
                        veg_classes=args.veg_classes,
                        max_roots=args.max_roots,
                        max_numnod=args.max_numnod,
                        cells=args.cells,
                        organic_fract=args.organic_fract,
                        spatial_frost=args.spatial_frost,
                        spatial_snow=args.spatial_snow,
                        july_tavg_supplied=args.july_tavg_supplied,
                        veglib_fcan=args.veglib_fcan,
                        veglib_photo=args.veglib_photo,
                        blowing_snow=args.blowing_snow,
                        vegparam_lai=args.vegparam_lai,
                        vegparam_fcan=args.vegparam_fcan,
                        vegparam_albedo=args.vegparam_albedo,
                        lai_src=args.lai_src,
                        fcan_src=args.fcan_src,
                        alb_src=args.alb_src,
                        lake_profile=args.lake_profile)

    print('completed grid_params.main(), output file was: {0}'.format(nc_file))
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def make_grid(grid_file, soil_file, snow_file, vegl_file, veg_file,
              lake_file, nc_file='params.nc', version_in='4.2',
              grid_decimal=4, nlayers=3, snow_bands=5, veg_classes=11,
              max_roots=3, max_numnod=10, cells=None,
              organic_fract=False, spatial_frost=False,
              spatial_snow=False, july_tavg_supplied=False,
              veglib_fcan=False, veglib_photo=False,
              blowing_snow=False, vegparam_lai=False,
              vegparam_fcan=False, vegparam_albedo=False,
              lai_src='FROM_VEGLIB', fcan_src='FROM_DEFAULT',
              alb_src='FROM_VEGLIB', lake_profile=False):
    """
    Make grid uses routines from params.py to read standard vic format
    parameter files.  After the parameter files are read, the files are placed
    onto the target grid using nearest neighbor mapping.  If a land mask is
    present in the target grid it will be used to exclude areas in the ocean.
    Finally, if the nc_file = 'any_string.nc', a netcdf file be written with
    the parameter data, if nc_file = False, the dictionary of grids is
    returned.
    """
    print('making gridded parameters now...')

    soil_dict = soil(soil_file, c=Cols(nlayers=nlayers,
                     organic_fract=organic_fract,
                     spatial_frost=spatial_frost,
                     spatial_snow=spatial_snow,
                     july_tavg_supplied=july_tavg_supplied))

    if cells is None:
        cells = len(soil_dict['gridcell'])

    if snow_file:
        snow_dict = snow(snow_file, soil_dict, c=Cols(snow_bands=snow_bands))
    else:
        snow_dict = False

    if vegl_file:
        veglib_dict, lib_bare_idx = veg_class(vegl_file,
                                              veglib_photo=veglib_photo,
                                              c=Cols(veglib_fcan=veglib_fcan,
                                                     veglib_photo=veglib_photo))
        veg_classes = len(veglib_dict['Veg_class'])
    else:
        veglib_dict = False
        lib_bare_idx = None

    if veg_file:
        veg_dict = veg(veg_file, soil_dict, veg_classes,
                       max_roots, cells, blowing_snow,
                       vegparam_lai, vegparam_fcan,
                       vegparam_albedo, lai_src,
                       fcan_src, alb_src)
    else:
        veg_dict = False

    if lake_file:
        lake_dict = lake(lake_file, soil_dict, max_numnod,
                         cells, lake_profile)
    else:
        lake_dict = False

    if grid_file:
        target_grid, target_attrs = read_netcdf(grid_file)
    else:
        target_grid, target_attrs = calc_grid(soil_dict['lats'],
                                              soil_dict['lons'],
                                              grid_decimal)

    grid_dict = grid_params(soil_dict, target_grid, snow_dict,
                            veglib_dict, veg_dict, lake_dict,
                            version_in, veglib_fcan,
                            veglib_photo, lib_bare_idx,
                            blowing_snow, vegparam_lai,
                            vegparam_fcan, vegparam_albedo,
                            lai_src, fcan_src, alb_src)

    if nc_file:
        write_netcdf(nc_file, target_attrs, target_grid,
                     grid_dict['soil_dict'], grid_dict['snow_dict'],
                     grid_dict['veg_dict'], grid_dict['lake_dict'],
                     version_in, organic_fract, spatial_frost,
                     spatial_snow, july_tavg_supplied,
                     veglib_fcan, veglib_photo, blowing_snow,
                     vegparam_lai, vegparam_fcan, vegparam_albedo,
                     lai_src, fcan_src, alb_src)
        return nc_file
    else:
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
        warn('lon_step ({0}) and lat_step ({1}) do not '
             'match'.format(lon_step, lat_step))
    if lat_count / len(ulats) < 0.95:
        warn('lat_count of mode is less than 95% ({0}%) of'
             ' len(lats)'.format(lat_count / len(ulats)))
    if lon_count / len(ulons) < 0.95:
        warn('lon_count of mode is less than 95% ({0}%) of'
             ' len(lons)'.format(lon_count / len(ulons)))

    if lats.min() < -55 and lats.max() > 70:
        # assume global grid
        print('assuming grid is meant to be global...')
        target_grid[XVAR] = np.linspace(-180 + lon_step[0] / 2,
                                        180 - lon_step[0] / 2,
                                        360 / lon_step[0])
        target_grid[YVAR] = np.linspace(-90 + lat_step[0] / 2,
                                        90 - lat_step[0] / 2,
                                        180 / lat_step[0])
    else:
        target_grid[XVAR] = np.arange(lons.min(),
                                      lons.max() + lon_step[0],
                                      lon_step[0])
        target_grid[YVAR] = np.arange(lats.min(),
                                      lats.max() + lat_step[0],
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
def grid_params(soil_dict, target_grid, snow_dict, veglib_dict, veg_dict,
                lake_dict=None, version_in='4.2', veglib_fcan=False,
                veglib_photo=False, lib_bare_idx=None, blowing_snow=False,
                vegparam_lai=False, vegparam_fcan=False,
                vegparam_albedo=False, lai_src='FROM_VEGLIB',
                fcan_src='FROM_DEFAULT', alb_src='FROM_VEGLIB'):
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
    if lake_dict:
        in_dicts['lake_dict'] = lake_dict
    else:
        out_dicts['lake_dict'] = False

    # get "unmasked" mask
    mask = target_grid['mask']

    ysize, xsize = target_grid['mask'].shape

    ymask, xmask = np.nonzero(mask != 1)

    print('{0} masked values'.format(len(ymask)))

    for name, mydict in in_dicts.items():
        out_dict = OrderedDict()

        for var in mydict:
            if mydict[var].dtype in [np.int, np.int64, np.int32]:
                fill_val = FILLVALUE_I
                dtype = np.int
            elif isinstance(mydict[var].dtype, np.str):
                fill_val = ''
                dtype = np.str
            else:
                fill_val = FILLVALUE_F
                dtype = np.float

            if mydict[var].ndim == 1:
                out_dict[var] = np.ma.zeros((ysize, xsize),
                                            dtype=dtype)
                out_dict[var][yi, xi] = mydict[var]
                out_dict[var][ymask, xmask] = fill_val

            elif mydict[var].ndim == 2:
                steps = mydict[var].shape[1]
                out_dict[var] = np.ma.zeros((steps, ysize, xsize),
                                            dtype=dtype)
                for i in pyrange(steps):
                    out_dict[var][i, yi, xi] = mydict[var][:, i]
                out_dict[var][:, ymask, xmask] = fill_val

            elif mydict[var].ndim == 3:
                j = mydict[var].shape[1]
                k = mydict[var].shape[2]
                out_dict[var] = np.ma.zeros((j, k, ysize, xsize),
                                            dtype=dtype)
                for jj in pyrange(j):
                    for kk in pyrange(k):
                        out_dict[var][jj, kk, yi, xi] = mydict[var][:, jj, kk]
                for y, x in pyzip(ymask, xmask):
                    out_dict[var][:, :, y, x] = fill_val

            out_dict[var] = np.ma.masked_values(out_dict[var], fill_val)

        out_dicts[name] = out_dict

    # Merge information from veglib and veg and transfer to veg_dict
    if veglib_dict and veg_dict:

        # Check for any missing Cv areas that will become bare soil
        var = 'Cv'
        bare = 1 - out_dicts['veg_dict'][var].sum(axis=0)
        bare[bare < 0.0] = 0.0

        # Determine the final number of veg classes, accounting for
        # potential new bare soil class, and determine that class's idx
        var = 'Cv'
        if lib_bare_idx is not None:
            # Bare soil class already exists at lib_bare_idx
            extra_class = 0
        else:
            # Bare soil class needs to be created
            extra_class = 1
            lib_bare_idx = out_dicts['veg_dict'][var].shape[0]
        nveg_classes = out_dicts['veg_dict'][var].shape[0] + extra_class

        # Transfer Cv info, accounting for additional bare soil areas
        var = 'Cv'
        shape = (nveg_classes, ) + out_dicts['veg_dict'][var].shape[1:]
        new = np.zeros(shape)
        if extra_class:
            new[:-1, :, :] = out_dicts['veg_dict'][var]
        else:
            new[:, :, :] = out_dicts['veg_dict'][var]
        new[lib_bare_idx, :, :] += bare
        # Ensure that Cvs sum to 1.0
        new /= new.sum(axis=0)
        new[:, ymask, xmask] = FILLVALUE_F
        out_dicts['veg_dict'][var] = new

        # Distribute the vegparam variables (geographically-varying)
        #   double root_depth(veg_class, root_zone, nj, ni) ;
        #   double root_fract(veg_class, root_zone, nj, ni) ;
        #   double LAI(veg_class, month, nj, ni) ;
        #   double fcan(veg_class, month, nj, ni) ;
        #   double albedo(veg_class, month, nj, ni) ;
        varnames = ['root_depth', 'root_fract']
        if vegparam_lai and lai_src == 'FROM_VEGPARAM':
            varnames.append('LAI')
        if vegparam_fcan and fcan_src == 'FROM_VEGPARAM':
            varnames.append('fcanopy')
        if vegparam_albedo and alb_src == 'FROM_VEGPARAM':
            varnames.append('albedo')
        for var in varnames:
            shape = (nveg_classes, ) + out_dicts['veg_dict'][var].shape[1:]
            new = np.full(shape, FILLVALUE_F)
            if extra_class:
                new[:-1, :, :] = out_dicts['veg_dict'][var]
                new[-1, :, :] = bare_vegparam[var]
            else:
                new[:, :, :] = out_dicts['veg_dict'][var]
            out_dicts['veg_dict'][var] = np.ma.masked_values(new, FILLVALUE_F)

        if blowing_snow:
            for var in ['sigma_slope', 'lag_one', 'fetch']:
                shape = (nveg_classes, ) + out_dicts['veg_dict'][var].shape[1:]
                new = np.full(shape, FILLVALUE_F)
                if extra_class:
                    new[:-1, :, :] = out_dicts['veg_dict'][var]
                    new[-1, :, :] = bare_vegparam[var]
                else:
                    new[:, :, :] = out_dicts['veg_dict'][var]
                out_dicts['veg_dict'][var] = np.ma.masked_values(new, FILLVALUE_F)

        # Distribute the veglib variables
        # 1st - the 1d vars
        #   double lib_overstory(veg_class) ;  --> (veg_class, nj, ni)
        varnames = ['overstory', 'rarc', 'rmin', 'wind_h', 'RGL', 'rad_atten',
                    'rad_atten', 'wind_atten', 'trunk_ratio']
        if veglib_photo:
            varnames.append('Ctype')
            varnames.append('MaxCarboxRate')
            varnames.append('MaxE_or_CO2Spec')
            varnames.append('LUE')
            varnames.append('Nscale')
            varnames.append('Wnpp_inhib')
            varnames.append('NPPfactor_sat')
        for var in varnames:
            lib_var = 'lib_{0}'.format(var)
            if var in ['Ctype', 'Nscale']:
                fill_val = FILLVALUE_I
                new = np.full((nveg_classes, ysize, xsize), fill_val,
                              dtype=np.int)
            else:
                fill_val = FILLVALUE_F
                new = np.full((nveg_classes, ysize, xsize), fill_val)
            if extra_class:
                new[:-1, yi, xi] = veglib_dict[lib_var][:, np.newaxis]
                new[-1, yi, xi] = bare_vegparam[var]
            else:
                new[:, yi, xi] = veglib_dict[lib_var][:, np.newaxis]
            new[:, ymask, xmask] = fill_val
            out_dicts['veg_dict'][var] = np.ma.masked_values(new, fill_val)

        # 2nd - the 2d vars
        varnames = ['veg_rough', 'displacement']
        if alb_src == 'FROM_VEGLIB':
            varnames = ['albedo'] + varnames
        if veglib_fcan and fcan_src == 'FROM_VEGLIB':
            varnames = ['fcanopy'] + varnames
        if lai_src == 'FROM_VEGLIB':
            varnames = ['LAI'] + varnames
        for var in varnames:
            lib_var = 'lib_{0}'.format(var)
            shape = (nveg_classes, veglib_dict[lib_var].shape[1],
                     ysize, xsize)
            new = np.full(shape, FILLVALUE_F)
            if extra_class:
                new[:-1, :, yi, xi] = veglib_dict[lib_var][:, :, np.newaxis]
                new[-1, :, yi, xi] = bare_vegparam[var]
            else:
                new[:, :, yi, xi] = veglib_dict[lib_var][:, :, np.newaxis]
            for y, x in pyzip(ymask, xmask):
                new[:, :, y, x] = fill_val
            out_dicts['veg_dict'][var] = np.ma.masked_values(new, FILLVALUE_F)

        # Finally, transfer veglib class descriptions (don't distribute)
        # This deviates from dimensions of other grid vars
        var = 'comment'
        lib_var = 'lib_{0}'.format(var)
        new = np.empty(nveg_classes, 'O')
        if extra_class:
            new[:-1] = veglib_dict['lib_comment']
            new[-1] = 'Bare Soil'
        else:
            new[:] = veglib_dict['lib_comment']
        out_dicts['veg_dict']['comment'] = new

    return out_dicts
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
#  Write output to netCDF
def write_netcdf(myfile, target_attrs, target_grid,
                 soil_grid=None, snow_grid=None, veg_grid=None,
                 lake_grid=None, version_in='4.2',
                 organic_fract=False, spatial_frost=False,
                 spatial_snow=False, july_tavg_supplied=False,
                 veglib_fcan=False, veglib_photo=False, blowing_snow=False,
                 vegparam_lai=False, vegparam_fcan=False,
                 vegparam_albedo=False, lai_src='FROM_VEGLIB',
                 fcan_src='FROM_DEFAULT', alb_src='FROM_VEGLIB'):
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

    if lake_grid:
        lakes = True
    else:
        lakes = False

    unit = Units(organic_fract=organic_fract, spatial_frost=spatial_frost,
                 spatial_snow=spatial_snow,
                 july_tavg_supplied=july_tavg_supplied,
                 veglib_fcan=veglib_fcan, veglib_photo=veglib_photo,
                 blowing_snow=blowing_snow,
                 vegparam_lai=vegparam_lai, vegparam_fcan=vegparam_fcan,
                 vegparam_albedo=vegparam_albedo, lakes=lakes)
    desc = Desc(organic_fract=organic_fract, spatial_frost=spatial_frost,
                spatial_snow=spatial_snow,
                july_tavg_supplied=july_tavg_supplied,
                veglib_fcan=veglib_fcan, veglib_photo=veglib_photo,
                blowing_snow=blowing_snow,
                vegparam_lai=vegparam_lai, vegparam_fcan=vegparam_fcan,
                vegparam_albedo=vegparam_albedo, lakes=lakes)

    # target grid
    # coordinates
    if target_grid[XVAR].ndim == 1:
        f.createDimension('lon', len(target_grid[XVAR]))
        f.createDimension('lat', len(target_grid[YVAR]))
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
        f.createDimension('nj', target_grid[XVAR].shape[0])
        f.createDimension('ni', target_grid[YVAR].shape[1])
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
            f.createDimension('nv4', 4)
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
        for name, attr in target_attrs[var].items():
            try:
                setattr(v, name, attr)
            except:
                print('dont have units or description for {0}'.format(var))

    # Layers
    f.createDimension('nlayer', soil_grid['soil_density'].shape[0])
    layer_dims = ('nlayer', ) + dims2

    v = f.createVariable('layer', NC_INT, ('nlayer', ))
    v[:] = np.arange(1, soil_grid['soil_density'].shape[0] + 1)
    v.long_name = 'soil layer'

    # soil grid
    for var, data in soil_grid.items():
        print('writing var: {0}'.format(var))

        if data.ndim == 1:
            v = f.createVariable(var, NC_DOUBLE, ('nlayer', ),
                                 fill_value=FILLVALUE_F)
            v[:] = data

        elif data.ndim == 2:
            if var in ['gridcell', 'run_cell', 'fs_active']:
                v = f.createVariable(var, NC_INT, dims2,
                                     fill_value=FILLVALUE_I)
            else:
                v = f.createVariable(var, NC_DOUBLE, dims2,
                                     fill_value=FILLVALUE_F)
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

        f.createDimension('snow_band', snow_grid['AreaFract'].shape[0])
        snow_dims = ('snow_band', ) + dims2

        v = f.createVariable('snow_band', NC_INT, ('snow_band', ))
        v[:] = np.arange(1, snow_grid['AreaFract'].shape[0] + 1)
        v.long_name = 'snow band'

        for var, data in snow_grid.items():
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

        f.createDimension('veg_class', veg_grid['Cv'].shape[0])
        f.createDimension('root_zone', veg_grid['root_depth'].shape[1])
        f.createDimension('month', MONTHS_PER_YEAR)

        v = f.createVariable('veg_class', NC_INT, ('veg_class', ))
        v[:] = np.arange(1, veg_grid['Cv'].shape[0] + 1)
        v.long_name = 'Vegetation Class'

        v = f.createVariable('veg_descr', np.dtype(str), ('veg_class', ))
        v[:] = veg_grid['comment']
        v.long_name = 'Vegetation Class Description'

        v = f.createVariable('root_zone', NC_INT, ('root_zone', ))
        v[:] = np.arange(1, veg_grid['root_depth'].shape[1] + 1)
        v.long_name = 'root zone'

        v = f.createVariable('month', NC_INT, ('month', ))
        v[:] = np.arange(1, 13)
        v.long_name = 'month of year'

        for var, data in veg_grid.items():
            if var != 'comment':
                print('writing var: {0} {1}'.format(var, data.shape))

                if veg_grid[var].ndim == 2:
                    if var  == 'Nveg':
                        v = f.createVariable(var, NC_INT, dims2,
                                             fill_value=FILLVALUE_I)
                    else:
                        v = f.createVariable(var, NC_DOUBLE, dims2,
                                             fill_value=FILLVALUE_F)
                    v[:, :] = data

                elif veg_grid[var].ndim == 3:
                    mycoords = ('veg_class', ) + dims2
                    if var in ['overstory', 'Ctype', 'Nscale']:
                        v = f.createVariable(var, NC_INT, mycoords,
                                             fill_value=FILLVALUE_I)
                    else:
                        v = f.createVariable(var, NC_DOUBLE, mycoords,
                                             fill_value=FILLVALUE_F)
                    v[:, :, :] = data

                elif var in ['LAI', 'fcanopy', 'albedo', 'veg_rough',
                             'displacement']:
                    mycoords = ('veg_class', 'month') + dims2
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

                v.long_name = var
                try:
                    v.units = unit.veg_param[var]
                    v.description = desc.veg_param[var]
                except KeyError:
                    lib_var = 'lib_{0}'.format(var)
                    v.units = unit.veglib[lib_var]
                    v.description = desc.veglib[lib_var]

                if coordinates:
                    v.coordinates = coordinates

    if lake_grid:
        if 'gridcell' in lake_grid:
            del lake_grid['gridcell']

        f.createDimension('lake_node', lake_grid['basin_depth'].shape[0])

        v = f.createVariable('lake_node', NC_INT, ('lake_node', ))
        v[:] = np.arange(1, lake_grid['basin_depth'].shape[0] + 1)
        v.long_name = 'lake basin node'

        for var, data in lake_grid.items():
            print('writing var: {0} {1}'.format(var, data.shape))

            if lake_grid[var].ndim == 2:
                if var in ['lake_idx', 'numnod']:
                    v = f.createVariable(var, NC_INT, dims2,
                                         fill_value=FILLVALUE_I)
                else:
                    v = f.createVariable(var, NC_DOUBLE, dims2,
                                         fill_value=FILLVALUE_F)
                v[:, :] = data

            elif lake_grid[var].ndim == 3:
                mycoords = ('lake_node', ) + dims2
                v = f.createVariable(var, NC_DOUBLE, mycoords,
                                     fill_value=FILLVALUE_F)
                v[:, :, :] = data

            else:
                raise ValueError('only able to handle dimensions <=3')

            v.long_name = var
            v.units = unit.lake_param[var]
            v.description = desc.lake_param[var]

            if coordinates:
                v.coordinates = coordinates

    f.close()

    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def soil(in_file, c=Cols(nlayers=3, organic_fract=False,
                         spatial_frost=False, spatial_snow=False,
                         july_tavg_supplied=False)):
    """
    Load the entire soil file into a dictionary of numpy arrays.
    Also reorders data to match gridcell order of soil file.
    """
    print('reading {0}'.format(in_file))
    data = np.loadtxt(in_file)

    soil_dict = OrderedDict()
    for var, columns in c.soil_param.items():
        if var in ['gridcell', 'run_cell', 'fs_active']:
            soil_dict[var] = np.squeeze(data[:, columns]).astype(int)
        else:
            soil_dict[var] = np.squeeze(data[:, columns])

    unique_grid_cells, inds = np.unique(soil_dict['gridcell'], return_index=True)
    if len(unique_grid_cells) != len(soil_dict['gridcell']):
        warn('found duplicate grid cells in soil file')
        for key, val in soil_dict.items():
            soil_dict[key] = val[inds]

    return soil_dict
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def snow(snow_file, soil_dict, c=Cols(snow_bands=5)):
    """
    Load the entire snow file into a dictionary of numpy arrays.
    Also reorders data to match gridcell order of soil file.
    """

    print('reading {0}'.format(snow_file))

    data = np.loadtxt(snow_file)

    snow_dict = OrderedDict()
    for var in c.snow_param:
        snow_dict[var] = data[:, c.snow_param[var]]

    # Make gridcell order match that of soil_dict
    cell_nums = snow_dict['cellnum'].astype(np.int)
    indexes = np.zeros_like(soil_dict['gridcell'], dtype=np.int)
    for i, sn in enumerate(soil_dict['gridcell'].astype(np.int)):
        try:
            indexes[i] = np.nonzero(cell_nums == sn)[0]
            last = indexes[i]
        except IndexError:
            warn('grid cell %d not found, using last known snowband' % sn)
            indexes[i] = last

    for var in snow_dict:
        snow_dict[var] = np.squeeze(snow_dict[var][indexes])

    return snow_dict
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def veg_class(vegl_file, veglib_photo=False,
              c=Cols(veglib_fcan=False, veglib_photo=False)):
    """
    Load the entire vegetation library file into a dictionary of lists.
    """

    print('reading {0}'.format(vegl_file))

    data = []
    sep = ' '
    lib_bare_idx = None
    row = 0
    with open(vegl_file, 'r') as f:
        for line in f:
            words = line.split()
            if row == 0:
                col_desc = len(words) - 1
            elif line.startswith('#'):
                continue  # skip additional lines with comments
            else:
                data.append([])
                for col in np.arange(0, col_desc):
                    data[row - 1][:] = words[:col_desc]
                    data[row - 1].append(sep.join(words[col_desc:]))
                if veglib_photo:
                    if data[row - 1][c.veglib['lib_Ctype']] == 'C3':
                        data[row - 1][c.veglib['lib_Ctype']] = 0
                    elif data[row - 1][c.veglib['lib_Ctype']] == 'C4':
                        data[row - 1][c.veglib['lib_Ctype']] = 1
                if re.match('(bare|barren|unvegetated)',
                            sep.join(words[col_desc:]), re.I):
                    lib_bare_idx = row - 1
            row += 1

    veglib_dict = OrderedDict()
    data = np.array(data)
    for var in c.veglib:
        veglib_dict[var] = np.squeeze(data[..., c.veglib[var]])
    return veglib_dict, lib_bare_idx
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def veg(veg_file, soil_dict, veg_classes=11, max_roots=3,
        cells=None, blowing_snow=False, vegparam_lai=False,
        vegparam_fcan=False, vegparam_albedo=False,
        lai_src='FROM_VEGLIB', fcan_src='FROM_DEFAULT',
        alb_src='FROM_VEGLIB'):
    """
    Read the vegetation file from vegFile.  Assumes max length for rootzones
    and vegclasses.  Also reorders data to match gridcell order of soil file.
    """

    print('reading {0}'.format(veg_file))

    with open(veg_file) as f:
        lines = f.readlines()

    if cells is None:
        cells = len(lines)

    gridcel = np.zeros(cells, dtype=np.int)
    nveg = np.zeros(cells, dtype=np.int)
    cv = np.zeros((cells, veg_classes))
    root_depth = np.zeros((cells, veg_classes, max_roots))
    root_fract = np.zeros((cells, veg_classes, max_roots))
    if blowing_snow:
        sigma_slope = np.zeros((cells, veg_classes))
        lag_one = np.zeros((cells, veg_classes))
        fetch = np.zeros((cells, veg_classes))
    lfactor = 1
    if vegparam_lai:
        lfactor += 1
        lai = np.zeros((cells, veg_classes, MONTHS_PER_YEAR))
    if vegparam_fcan:
        lfactor += 1
        fcan = np.zeros((cells, veg_classes, MONTHS_PER_YEAR))
    if vegparam_albedo:
        lfactor += 1
        albedo = np.zeros((cells, veg_classes, MONTHS_PER_YEAR))

    row = 0
    cell = 0
    while row < len(lines):
        line = lines[row].strip('\n').split(' ')
        gridcel[cell], nveg[cell] = np.array(line).astype(int)
        numrows = nveg[cell] * lfactor + row + 1
        row += 1

        while row < numrows:
            lines[row] = lines[row].strip()
            line = lines[row].strip('\n').split(' ')
            temp = np.array(line).astype(float)
            vind = int(temp[0]) - 1
            cv[cell, vind] = temp[1]

            tmp = 1 + max_roots * 2
            root_depth[cell, vind, :] = temp[2:tmp:2]
            root_fract[cell, vind, :] = temp[3:1+tmp:2]
            tmp += 1

            if blowing_snow:
                sigma_slope[cell, vind] = temp[tmp]
                lag_one[cell, vind] = temp[1+tmp]
                fetch[cell, vind] = temp[2+tmp]
                tmp += 2

            row += 1

            if vegparam_lai:
                lines[row] = lines[row].strip()
                line = lines[row].strip('\n').split(' ')
                lai[cell, vind, :] = np.array(line, dtype=np.float)
                row += 1
            if vegparam_fcan:
                lines[row] = lines[row].strip()
                line = lines[row].strip('\n').split(' ')
                fcan[cell, vind, :] = np.array(line, dtype=np.float)
                row += 1
            if vegparam_albedo:
                lines[row] = lines[row].strip()
                line = lines[row].strip('\n').split(' ')
                albedo[cell, vind, :] = np.array(line, dtype=np.float)
                row += 1
        cell += 1
    veg_dict = OrderedDict()
    veg_dict['gridcell'] = gridcel[:cell]
    veg_dict['Nveg'] = nveg[:cell]
    veg_dict['Cv'] = cv[:cell, :]
    veg_dict['root_depth'] = root_depth[:cell, :, :]
    veg_dict['root_fract'] = root_fract[:cell, :, :]

    if blowing_snow:
        veg_dict['sigma_slope'] = sigma_slope[:cell, :]
        veg_dict['lag_one'] = lag_one[:cell, :]
        veg_dict['fetch'] = fetch[:cell, :]

    if vegparam_lai and lai_src == 'FROM_VEGPARAM':
        veg_dict['LAI'] = lai[:cell, :, :]

    if vegparam_fcan and fcan_src == 'FROM_VEGPARAM':
        veg_dict['fcanopy'] = fcan[:cell, :, :]

    if vegparam_albedo and alb_src == 'FROM_VEGPARAM':
        veg_dict['albedo'] = albedo[:cell, :, :]

    # Make gridcell order match that of soil_dict
    inds = []
    for sn in soil_dict['gridcell']:
        inds.append(np.nonzero(veg_dict['gridcell'] == sn))

    new_veg_dict = OrderedDict()
    for var in veg_dict:
        new_veg_dict[var] = np.squeeze([veg_dict[var][i] for i in inds])

    return new_veg_dict
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def lake(lake_file, soil_dict, max_numnod=10,
         cells=None, lake_profile=False):
    """
    Read the lake file from lakeFile.  Assumes max length for depth-area
    relationship.  Also reorders data to match gridcell order of soil file.
    """

    print('reading {0}'.format(lake_file))

    with open(lake_file) as f:
        lines = f.readlines()

    if cells is None:
        cells = len(lines)

    gridcel = np.zeros(cells, dtype=np.int)
    lake_idx = np.zeros(cells, dtype=np.int)
    numnod = np.zeros(cells, dtype=np.int)
    mindepth = np.zeros(cells)
    wfrac = np.zeros(cells)
    depth_in = np.zeros(cells)
    rpercent = np.zeros(cells)
    if lake_profile:
        basin_depth = np.zeros((cells, max_numnod))
        basin_area = np.zeros((cells, max_numnod))
    else:
        basin_depth = np.zeros(cells)
        basin_area = np.zeros(cells)

    row = 0
    cell = 0
    while row < len(lines):
        line = lines[row].strip('\n').split(' ')
        gridcel[cell], lake_idx[cell], numnod[cell] = np.array(line[0:3],
                                                               dtype=np.int)
        temp = np.array(line, dtype=np.float)
        mindepth[cell] = temp[3]
        wfrac[cell] = temp[4]
        depth_in[cell] = temp[5]
        rpercent[cell] = temp[6]
        numrows = row + 1
        if lake_idx[cell] >= 0:
            numrows += 1
        row += 1

        while row < numrows:
            lines[row] = lines[row].strip()
            line = lines[row].strip('\n').split(' ')
            temp = np.array(line, dtype=np.float)

            if lake_profile:
                rind = len(temp) / 2
                basin_depth[cell, :rind] = temp[0::2]
                basin_area[cell, :rind] = temp[1::2]
            else:
                basin_depth[cell] = temp[0]
                basin_area[cell] = temp[1]

            row += 1
        cell += 1

    lake_dict = OrderedDict()
    lake_dict['gridcell'] = gridcel[:cell]
    lake_dict['lake_idx'] = lake_idx[:cell]
    lake_dict['numnod'] = numnod[:cell]
    lake_dict['mindepth'] = mindepth[:cell]
    lake_dict['wfrac'] = wfrac[:cell]
    lake_dict['depth_in'] = depth_in[:cell]
    lake_dict['rpercent'] = rpercent[:cell]
    if lake_profile:
        lake_dict['basin_depth'] = basin_depth[:cell, :]
        lake_dict['basin_area'] = basin_area[:cell, :]
    else:
        lake_dict['basin_depth'] = basin_depth[:cell]
        lake_dict['basin_area'] = basin_area[:cell]

    # Make gridcell order match that of soil_dict
    inds = []
    for sn in soil_dict['gridcell']:
        inds.append(np.nonzero(lake_dict['gridcell'] == sn))

    new_lake_dict = OrderedDict()
    for var in lake_dict:
        new_lake_dict[var] = np.squeeze([lake_dict[var][i] for i in inds])

    return new_lake_dict
# -------------------------------------------------------------------- #
