from met_preprocessor.unit_conv import UnitConversion
from met_preprocessor.met_preprocessing import get_unit_conv_params
import pytest


@pytest.fixture(scope="module")
def unit_conv_params(param_map):
    return get_unit_conv_params(param_map)


@pytest.fixture(scope="module")
def param_conv(unit_conv_params):
    return UnitConversion(unit_conv_params)


def test_unit_conv(sample_xarray_data, param_conv):
    print(param_conv.convert_param(sample_xarray_data["Tair"], "kelvin"))


# TODO: Unit conv for month
