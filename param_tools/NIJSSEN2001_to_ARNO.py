#!/usr/bin/env python

"""
NIJSSEN2001_to_ARNO.py

A simple script that converts VIC NIJSSEN2001 soil parameters to
standard VIC ARNO baseflow parameters.
"""

import argparse
from netCDF4 import Dataset
from ncparam2ascii import read_netcdf


# -------------------------------------------------------------------- #
def main():

    nc_params = process_command_line()

    data, atts = read_netcdf(nc_params)

    data = convert(data)

    write(nc_params, data)

    return
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
#
def write(nc_params, data):
    """Update netCDF file with new parameters"""

    f = Dataset(nc_params, 'r+')

    note1 = 'ARNO baseflow parameter'
    variables = ['Ds', 'Dsmax', 'Ws', 'c']
    for var in variables:
        print('Writing updated param: {}'.format(var))
        f.variables[var][:] = data[var]
        f.variables[var].note = note1

    note2 = 'Converted NIJSSEN2001 baseflow params to ARNO baseflow params'
    try:
        f.history += note2
    except:
        f.history = note2

    f.close()

    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def convert(data):
    """Convert baseflow parameters to ARNO style"""

    max_moist = calc_max_moist(data['depth'],
                               data['bulk_density'],
                               data['soil_density'])

    Ds, Dsmax, Ws, c = calc_params(data['Ds'], data['Dsmax'],
                                   data['Ws'], data['c'],
                                   max_moist[-1])

    data['Ds'], data['Dsmax'], data['Ws'], data['c'] = Ds, Dsmax, Ws, c

    return data
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def calc_max_moist(depth, bulk_density, soil_density):
    """ calculate the maximum soil moisture of each layer """
    porosity = 1.0 - bulk_density / soil_density
    max_moist = depth * porosity * 1000.
    return max_moist
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def calc_params(d1, d2, d3, d4, max_moist):
    """
    Convert parameters

    VIC Code:
        if(options.BASEFLOW == NIJSSEN2001) {
            layer = options.Nlayer-1;
            temp.Dsmax = temp.Dsmax *
                pow((double)(1./(temp.max_moist[layer]-temp.Ws)), -temp.c) +
                temp.Ds * temp.max_moist[layer];
            temp.Ds = temp.Ds * temp.Ws / temp.Dsmax;
            temp.Ws = temp.Ws/temp.max_moist[layer];


    d? - values corresponding to the 4 NIJSSEN2001 baseflow parameters
    max_moist - maximum moisture of the bottom soil layer
    """

    Dsmax = d2 * pow(1./(max_moist - d3), -d4) + d1 * max_moist
    Ds = d1 * d3 / Dsmax
    Ws = d3 / max_moist
    c = d4

    return Ds, Dsmax, Ws, c
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def process_command_line():
    """
    Process command line arguments.
    """
    parser = argparse.ArgumentParser(description='Simple script to convert '
                                                 'between NIJSSEN2001 and ARNO'
                                                 ' baseflow parameters')
    parser.add_argument("nc_params",
                        type=str,
                        help="Input netCDF VIC parameter file")

    args = parser.parse_args()

    return args.nc_params
# -------------------------------------------------------------------- #

if __name__ == "__main__":
    main()
