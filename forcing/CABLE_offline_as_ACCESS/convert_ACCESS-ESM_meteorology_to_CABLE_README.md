## Purpose

Convert the meteorology data from a given ACCESS-ESM1.5 experiment/variation to a format readable by CABLE. The format chosen at this time is picked to match the GSWP format, which is a single file for each year of data for each variable.

## Usage

Call from the command line via

```convert_ACCESS-ESM1p5_meteorology_to_CABLE.py -e <experiment> -v <variation> -y <yearstart,yearend> -o <output>```,

where:

* ```--experiment(-e) <experiment>``` is the experiment name e.g. historical, piControl. Defaults to ```historical```.
* ```--variation(-v) <variation>``` is the variation code for the run. Defaults to ```r1i1p1f1```.
* ```--years(-y) <startyear,endyear>``` is the range of years to include in the conversion. Defaults to ```None```, which converts all years in the run.
* ```--output(-o) <output>``` is a template to write to. The template is appended by ```_<variable>_<year>.nc``` e.g. year 1980 in the precipitation would be ```<output>_rainf_1980.nc```. Defaults to ```ACCESS-ESM1p5-to-CABLE```.

The script facilitates writing the variables in parallel, with the number of processes being the number of available CPUs. The only package requirement is xarray, which is available through the ```hh5``` ```conda_concept``` environment on Gadi.

## Method

All the variables bar the surface wind speed have a 1 to 1 mapping between the ACCESS-ESM output and the required CABLE forcing. [xarray](https://docs.xarray.dev/en/stable/) is used for the conversion process, and Python's intrinsic [multiprocessing](https://docs.xarray.dev/en/stable/) package is used for the parallelism.

The time units are converted from days, which is the ACCESS-ESM1.5 default, to seconds, which the GSWP format in CABLE requires.

The wind speed is not available at 3 hour frequency, but the wind components are available on the u/v offset grids at 3 hour frequency. The components on the offset grids are bilinearly interpolated (see the [xesmf documentation](https://xesmf.readthedocs.io/en/stable/notebooks/Dataset.html); unfortunately details on the interpolation implementation are lacking) to a grid common with the rest of the meteorological variables then the Euclidean norm is applied to compute the magnitude.

The reference date for the time axis is changed to the first year selected for the given conversion. The time between reference date of 1850-01-01 and the start of the ACCESS-ESM1.5 data at 1960-01-01, when converted to seconds, is actually larger than the maximum signed integer value. The reference time is changed to reduce the magnitude of the values in the time axis.
