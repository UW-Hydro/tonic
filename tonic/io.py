#!/usr/bin/env python
"""Input/Output functions"""
from netCDF4 import Dataset
import configobj
from .pycompat import OrderedDict, SafeConfigParser


# -------------------------------------------------------------------- #
def read_config(config_file, default_config=None):
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

    if default_config is not None:
        for name, section in dict1.items():
            if name in default_config.keys():
                for option, key in default_config[name].items():
                    if option not in section.keys():
                        dict1[name][option] = key

    return dict1
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def read_configobj(config_file, default_config=None):
    """
    Return a dictionary with nested dictionaries
    """

    config = configobj.ConfigObj(config_file)

    config = type_configobj(config)

    return config
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def type_configobj(d):
    """recursively loop through dictionary calling config_type"""
    for k, v in d.items():
        if isinstance(v, dict):
            type_configobj(v)
        else:
            d[k] = config_type(v)
    return d
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
            return list(map(int, val_list))
        except:
            pass
        try:
            return list(map(float, val_list))
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
