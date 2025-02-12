import pint
from pint import Unit
from xarray import DataArray
from metpy.units import units


def rain_conversion(units: Unit, depth_time: Unit):
    """Conversion from precipitation in terms of depth
    (<length> / <time>) to intensity (kg / m^2 / s^-1)"""
    return (1000 * units.water) * (1 * depth_time).to_base_units()


class UnitConversion:

    def __init__(self, params: list[str]) -> None:
        self.contexts = {}
        self._add_unit_conversions(params)

    def _add_unit_conversions(self, params: list[str]) -> None:
        """Define additional unit conversions other than default ones
        in `metpy`/`pint`."""
        # TODO: Load from file
        units.define("Celsius = degC")
        units.define("HPa = 100 Pa")

        for param in params:
            self._add_new_empty_context(param)

        if "Rainf" in self.contexts:
            self.contexts["Rainf"].add_transformation(
                "[length] / [time]",
                "[mass] / [length] ** 2 / [time]",
                rain_conversion,
            )

        self._add_new_empty_context("month")

    def _add_new_empty_context(self, param: str) -> None:
        self.contexts[param] = pint.Context(param)
        units.add_context(self.contexts[param])

    def _convert_units(self, da: DataArray, new_unit: str) -> DataArray:
        # Prevent unnecessarily quantify-ing (expensive operation)
        if new_unit != da.units:
            da_q = da.metpy.quantify()
            da_q = da_q.metpy.convert_units(new_unit)
            da_q = da_q.metpy.dequantify()
        return da_q

    def _monthly_conversions(self, da: DataArray) -> DataArray:
        """Convert monthly to daily data."""
        print(f"Monthly conversions for {da.name}")
        daily_units = pint.util.to_units_container(units(da.units))
        daily_units = daily_units.rename("month", "day")
        for monthly_days in range(28, 32):
            da_time = da.isel(time=da.time.dt.days_in_month == monthly_days)
            scaling_factor = monthly_days
            self.contexts["month"].redefine(f"month = {scaling_factor} * days")
            if da_time.time.size == 0:
                continue
            time_filter_dict = dict(time=da_time.time)
            with units.context("month"):
                da.loc[time_filter_dict] = self._convert_units(
                    da.loc[time_filter_dict],
                    str(daily_units),
                )
        da.attrs["units"] = str(daily_units)
        return da

    def convert_param(self, da: DataArray, out_units: str) -> DataArray:
        """Convert parameter into necessary units."""
        print(f"Converting param: {da.name}")
        with units.context(da.name):
            if "month" in str(da.units):
                da = self._monthly_conversions(da)
            return self._convert_units(da, out_units)
