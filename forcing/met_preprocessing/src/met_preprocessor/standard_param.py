import metpy.calc as mpcalc
from metpy.units import check_units
from metpy.xarray import preprocess_and_wrap


# Specific humidity calculations
@preprocess_and_wrap(wrap_like="tair", broadcast=("vp", "vpd", "tair"))
@check_units("[pressure]", "[pressure]", "[temperature]")
def vp_vpd_tair_sh(vp, vpd, tair):
    svp = vp + vpd
    return _calc_sh(vp, svp, tair)


@preprocess_and_wrap(wrap_like="tair", broadcast=("vp", "tair"))
@check_units("[pressure]", "[temperature]")
def vp_tair_sh(vp, tair):
    svp = mpcalc.saturation_vapor_pressure(tair)
    return _calc_sh(vp, svp, tair)


@preprocess_and_wrap(wrap_like="tair", broadcast=("vpd", "tair"))
@check_units("[pressure]", "[temperature]")
def vpd_tair_sh(vpd, tair):
    svp = mpcalc.saturation_vapor_pressure(tair)
    vp = svp - vpd
    return _calc_sh(vp, svp, tair)


def _calc_sh(vp, svp, tair):
    relative_humidity = vp / svp
    mixing_ratio = mpcalc.mixing_ratio_from_relative_humidity(
        vp, tair, relative_humidity
    )
    specific_humidity = mpcalc.specific_humidity_from_mixing_ratio(mixing_ratio)

    return specific_humidity
