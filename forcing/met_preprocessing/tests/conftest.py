import pytest
import pandas as pd
import xarray as xr
import numpy as np
import yaml

TEST_PARAM_MAP_FILE = "tests/data/test_param_map.yaml"


@pytest.fixture(scope="module")
def param_map():
    with open(TEST_PARAM_MAP_FILE) as file:
        param_map = yaml.safe_load(file)
    return param_map


# TODO: Add integration test for unchanged temp (have NetCDF)
# TODO: Similar for leap year
# TODO: Trusted short script for vpd / calc


@pytest.fixture(scope="module")
def sample_xarray_data():
    lon = [-99.83, -99.32]
    lat = [42.25, 42.21]
    time = pd.date_range("2024-09-06", periods=4)
    temperature = 15 + 8 * np.random.randn(2, 2, 4)
    rain = 500 + 1 * np.random.randn(2, 2, 4)
    vpd = 15 + 8 * np.random.randn(2, 2, 4)
    reference_time = pd.Timestamp("2014-09-05")
    ds = xr.Dataset(
        data_vars=dict(
            Tair=(["lon", "lat", "time"], temperature, {"units": "celsius"}),
            Rainf=(["lon", "lat", "time"], rain, {"units": "mm/month"}),
            vpd=(["lon", "lat", "time"], vpd, {"units": "Pa"}),
        ),
        coords=dict(
            lon=lon,
            lat=lat,
            time=time,
            reference_time=reference_time,
        ),
        attrs=dict(description="Test dataset."),
    )
    return ds
