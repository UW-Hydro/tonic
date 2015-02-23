#!/usr/bin/env python
"""functions to support testing of VIC output"""

from __future__ import print_function
# import pandas as pd


# -------------------------------------------------------------------- #
class VICTestError(BaseException):
    pass
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def check_for_nans(df):
    """check if dataframe has nans in it"""
    if df.isnull().any().any():
        raise VICTestError('VIC output has nans in it!')
    else:
        return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def check_completed(df, start, end):
    if not df.index[0] == start:
        raise VICTestError('Start dates ({0} and {1}) do not match'.format(df.index[0], start))
    if not df.index[-1] == end:
        raise VICTestError('End dates ({0} and {1}) do not match'.format(df.index[-1], end))
    return

# -------------------------------------------------------------------- #


# # -------------------------------------------------------------------- #
# def check_water_balance(df):
#     """
#     Calculate the water balance error from VIC output.

#     Parameters
#     ----------
#     df : Pandas DataFrame with the following variables:
#         WATER_ERROR     :   Water balance error

#             or

#         PREC            :   Precipitation for current record
#         EVAP            :   Evaporation for current record
#         RUNOFF          :   Runoff for current record
#         BASEFLOW        :   Baseflow for current record
#         DELINTERCEPT    :   Change in intercepted storage
#         DELSWE          :   Change in snow water equivalent
#         DELSURFSTOR     :   Change in surface water storage
#         DELSOILMOIST    :   Change in soil moisture content
#         WDEW            :   Canopy interception of liquid water
#         SWE             :   Snow water equivalent
#         SURFSTOR        :   Surface storage
#         SOIL_MOIST      :   Soil Moisture

#     Returns
#     -------
#     df : Pandas Dataframe.  As original including WATER_ERROR time series of water balance errors.
#     total_error : scalar (mm)

# inflow  = out_data[OUT_PREC].data[0] + out_data[OUT_LAKE_CHAN_IN].data[0]; // mm over grid cell
#   outflow = out_data[OUT_EVAP].data[0] + out_data[OUT_RUNOFF].data[0] + out_data[OUT_BASEFLOW].data[0]; // mm over grid cell
#   storage = 0.;
#   for(index=0;index<options.Nlayer;index++)
#     if(options.MOISTFRACT)
#       storage += (out_data[OUT_SOIL_LIQ].data[index] + out_data[OUT_SOIL_ICE].data[index]) 
#     * depth[index] * 1000;
#     else
#       storage += out_data[OUT_SOIL_LIQ].data[index] + out_data[OUT_SOIL_ICE].data[index];
#   storage += out_data[OUT_SWE].data[0] + out_data[OUT_SNOW_CANOPY].data[0] + out_data[OUT_WDEW].data[0] + out_data[OUT_SURFSTOR].data[0];
#   out_data[OUT_WATER_ERROR].data[0] = calc_water_balance_error(rec,inflow,outflow,storage);


# what to do for storage?

#     """

#     if 'WATER_ERROR' not in df:
#         storage = 
#         df['WATER_ERROR'] = (df['PREC'] + df['LAKE_CHAN_IN']) - \
#                             (df['EVAP'] + df['RUNOFF'] + df['BASEFLOW']) - \
#                             (storage - last_storage)

#     return df, total_error
# # -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def check_energy_balance():
    """
    Calculate the energy balance error from VIC output.

    Parameters
    ----------
    df : Pandas DataFrame with the following variables:
        NET_SHORT       :   Net shortwave radiation at the surface
        R_NET           :   Net radiation at the surface, includes long and shortwave radiation
        LATENT          :   Latent heat from the surface
        SENSIBLE        :   Sensible heat flux from the surface
        GRND_FLUX       :   Ground heat flux plus heat storage in the top soil layer
        DELTAH          :   Rate of change in heat storage
        FUSION          :   Net energy used to melt/freeze soil moisture
        IN_LONG         :   Incoming longwave at ground surface (under vegetation)
        SHORTWAVE       :   Incoming shortwave at ground surface (above vegetation)
        LONGWAVE        :   Incoming longwave at ground surface (above vegetation)
        LATENT          :   Latent heat from the surface


this...
        ENERGY_ERROR
        NET_SHORT
        NET_LONG
        LATENT
        LATENT_SUB
        SENSIBLE
        ADV_SENS
        GRND_FLUX
        DELTAH
        FUSION
        ADVECTION
        DELTACC
        SNOW_FLUX
        RFRZ_ENERGY

    Returns
    -------
    df : Pandas DataFrame.  As original including ENERGY_ERROR time series of energy balance errors.
    total_error : scalar (mm)
    """

    if 'ENERGY_ERROR' not in df:
        df['ENERGY_ERROR'] = (df['NET_SHORT'] + df['NET_LONG']) - \
                             (df['LATENT'] + df['LATENT_SUB'])  - \
                             (df['SENSIBLE'] + df['ADV_SENS'])  - \
                             (df['GRND_FLUX'] + df['DELTAH'] + df['FUSION']) + \
                             (df['ADVECTION'] - df['DELTACC'] + 
                              df['SNOW_FLUX'] + df['RFRZ_ENERGY'])

    total_error = df['ENERGY_ERROR'].sum()

    return df, total_error
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def check_forcings_integrity():
    return
# -------------------------------------------------------------------- #
