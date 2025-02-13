from metpy.units import units, check_units
from metpy.xarray import preprocess_and_wrap
import metpy.constants as c
import xarray as xr


# REVIEW: Metpy has recently added triple point as well (c.T0)
T0 = units.Quantity(273.16, "kelvin")


# TODO: Check whether magnitude can be replace with to_magnitude and to_base_units with process_units
@preprocess_and_wrap(wrap_like="temperature")
@check_units("[temperature]")
def calc_lwdown_swinbank(temperature):
    """Longwave radiation (W/m^2) calculation using Swinbank's formula
    https://doi.org/10.1002/qj.49708938105
    """
    t = temperature.to("kelvin").m
    return 0.0000094 * 0.0000000567 * (t**6.0) * units("W/m^2")


@preprocess_and_wrap(wrap_like="temperature", broadcast=("temperature", "elevation"))
@check_units("[temperature]", "[length]")
def calc_psurf(temperature, elevation):
    t = temperature.to("kelvin").m
    e = elevation.to("m").m
    return (
        1013.25 * t + (0.0065 * e) ** (c.earth_gravity.m / c.Rd.m / 0.0065)
    ) * units("Pa")


@preprocess_and_wrap(wrap_like="temperature", broadcast=("temperature", "rain"))
@check_units("[temperature]", "[length] / [time]")
def calc_snow(temperature, rain):
    t = temperature.to("kelvin").m
    return xr.where(t < T0.m, rain, 0.0, keep_attrs=True)
