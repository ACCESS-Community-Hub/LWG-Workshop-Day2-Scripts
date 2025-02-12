from metpy.units import units, check_units
from metpy.xarray import preprocess_and_wrap


# TODO: Check whether magnitude can be replace with to_magnitude and to_base_units with process_units
@preprocess_and_wrap(wrap_like="temperature")
@check_units("[temperature]")
def calc_lwdown_swinbank(temperature):
    """Longwave radiation (W/m^2) calculation using Swinbank's formula
    https://doi.org/10.1002/qj.49708938105
    """
    t = temperature.to("kelvin").magnitude
    return 0.0000094 * 0.0000000567 * (t**6.0) * units("W/m^2")
