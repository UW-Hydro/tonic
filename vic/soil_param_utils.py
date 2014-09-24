#!/usr/bin/env python

"""
soil_param_utils.py

A simple script that converts VIC NIJSSEN2001 soil parameters to
standard VIC ARNO baseflow parameters.
"""

from __future__ import print_function
from .share import MMPERMETER


# -------------------------------------------------------------------- #
def calc_max_moist(depth, bulk_density, soil_density):
    """ calculate the maximum soil moisture of each layer """
    porosity = 1.0 - bulk_density / soil_density
    max_moist = depth * porosity * MMPERMETER
    return max_moist
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def nijssen2001_to_arno(d1, d2, d3, d4, max_moist):
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
def anro_to_nijssen2001(Ds, Dsmax, Ws, c, max_moist):
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
    d4 = c
    d3 = max_moist * Ws
    d1 = Ds * Dsmax / d3
    d2 = (Dsmax - d1 * max_moist) / pow(max_moist - d3, d4)

    return d1, d2, d3, d4
# -------------------------------------------------------------------- #
