from met_preprocessor.met_preprocessing import order_load_dep, process_dependencies
from met_preprocessor.opt_param import calc_lwdown_swinbank
from met_preprocessor.standard_param import vpd_tair_sh
import pytest


@pytest.fixture()
def test_params(sample_xarray_data):
    return list(sample_xarray_data.keys())


@pytest.fixture()
def test_ord_deps():
    return [
        ("Qair", ["vpd", "Tair"], vpd_tair_sh),
        ("LWDown", ["Tair"], calc_lwdown_swinbank),
    ]


@pytest.fixture()
def test_all_deps(param_map):
    return process_dependencies(param_map)


def test_order_load_dep(test_params, test_all_deps, test_ord_deps):
    assert test_ord_deps == order_load_dep([], test_all_deps, test_params)
