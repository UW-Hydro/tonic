"""Set to run with pytest

Usage: py.test
"""
import pytest

import numpy as np
from tonic.tonic import calc_grid


@pytest.fixture(scope="function")
def lons():
    return np.arange(0, 5, 0.25)

    
@pytest.fixture(scope="function")
def lats():
    return np.arange(0, 10, 0.5)


def test_calc_grid_standard(lons, lats):
    shape = np.meshgrid(lons,lats)[0].shape
    target_grid = calc_grid(lons, lats, decimals=4)
    assert type(target_grid) == dict
    assert target_grid['mask'].shape == shape
