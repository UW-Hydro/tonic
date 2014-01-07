#!/usr/local/env python
"""
vic2netcdf_unit_test.py

Set to run with pytest

Usage: py.test (from RVIC or test directory)
"""
import pytest
import numpy as np
import sys
sys.path.append("../")

# -------------------------------------------------------------------- #
# Unit tests for vic2netcdf.py


@pytest.fixture
def lons():
	pass

def test_calc_grid_standard():
    from processing_tools.vic2netcdf import calc_grid
    lons = np.arange(0, 5, 0.5)
    lats = np.arange(0, 10, 0.5)
    target_grid = calc_grid(lons, lats, decimals=4)

# -------------------------------------------------------------------- #


