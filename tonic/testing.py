#!/usr/bin/env python
"""functions to support testing of VIC output"""

from __future__ import print_function
# import pandas as pd


# -------------------------------------------------------------------- #
class VICTestError(Exception):
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
        raise VICTestError(
            'Start dates ({0} and {1}) do not match'.format(df.index[0], start))
    if not df.index[-1] == end:
        raise VICTestError(
            'End dates ({0} and {1}) do not match'.format(df.index[-1], end))
    return

# -------------------------------------------------------------------- #
