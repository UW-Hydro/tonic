"""
vic2netcdf_unit_test.py

Set to run with pytest

Usage: py.test (from RVIC or test directory)
"""
import numpy as np
from tonic.tonic import calc_grid


def test_calc_grid_standard():
    lons = np.arange(0, 5, 0.25)
    lats = np.arange(0, 10, 0.5)
    target_grid = calc_grid(lons, lats, decimals=4)
    assert type(target_grid) == dict
