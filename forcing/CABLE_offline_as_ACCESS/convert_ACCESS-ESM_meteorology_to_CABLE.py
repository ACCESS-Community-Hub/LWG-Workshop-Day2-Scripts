import argparse
import datetime
import itertools
import multiprocessing
import numpy

import xarray
import xesmf

#-----------------------------------------------------------------------------#

def read_CLI_args():
    """Process the command line arguments and return the args."""

    parser = argparse.ArgumentParser(
            prog='convert_ESM_to_CABLE',
            description='Converts the meteorology from a given experiment ' +\
                    'to a format readable by CABLE. Splits file into ' +\
                    'individual years.'
                    )

    parser.add_argument(
            '-e',
            '--experiment',
            help='CMIP6 experiment to read from.',
            default='historical'
            )

    parser.add_argument(
            '-p',
            '--perturbation',
            help='Perturbation within the experiment to read from',
            default='r1i1p1f1'
            )

    parser.add_argument(
            '-y',
            '--years',
            help='<StartYear,EndYear> to set the range of year to convert. ' +\
                    'Defaults to all available years.',
            default=None
            )

    parser.add_argument(
            '-o',
            '--output',
            help='Prefix to apply to the output files ' +\
                    'i.e. <output>_<variable>_<year>.nc.',
            default='ACCESS-ESM1p5-to-CABLE'
            )

    return parser.parse_args()

#-----------------------------------------------------------------------------#

def map_to_pm180(Data):
    """Map the given dataset to a longitude domain of (-180.0, 180.0]."""
    RollDist = (Data['lon'] > 180.0).argmax().item()

    return Data.roll(lon=RollDist)

#-----------------------------------------------------------------------------#

def convert_variable(
        InDataset,
        InVariable,
        OutTemplate,
        OutVariable,
        Scaling=1,
        YearRange=None
        ):
    """Convert the variable InVariable in InDataset to data stored by year in
    <OutTemplate>_<year>.nc in variable OutVariable. Scale the data by scaling
    and add the supplied attrs to the variable. A YearRange can be supplied as
    a tuple of (YearStart, YearEnd)."""

    print(f"Getting data from {InDataset}")
    # Open as an mfdataset
    InDataset = xarray.open_mfdataset(
            InDataset,
            chunks={'time': 365*8},     # 3 hourly timesteps, a year of data
            engine = 'h5netcdf'
            )

    # Rename the variable to one recogniseble by GSWP routines and set attrs
    InDataset = InDataset.rename(
            {InVariable: OutVariable}
            )
    
    # Determine the time domain of the data, if not provided
    if YearRange is not None:
        StartTimeDomain = YearRange[0]
        EndTimeDomain = YearRange[1]
    else:
        StartTimeDomain = InDataset.time[0].dt.year.item()
        EndTimeDomain = InDataset.time[-1].dt.year.item()

    # The input data is in days since ..., but we convert to seconds since
    # by setting the time units attribute
    # The ACCESS-ESM data uses 1850 as the reference year, and when converting
    # the time since reference to seconds, we actually hit the size limit for
    # signed integers. To temporarily get around this, set the reference year
    # to the start year of the selected slice
    StartTimeAsSec = f'seconds since {StartTimeDomain}-01-01'

    # Set the units encoding
    InDataset.time.encoding['units'] = StartTimeAsSec
    if 'time_bnds' in InDataset.variables:
        InDataset.time_bnds.encoding['units'] = StartTimeAsSec

    # For each year, write a new file with that year's data
    for Year in range(StartTimeDomain, EndTimeDomain+1):
        # Select the data from that year
        InDatasetSlice = InDataset.sel(
                time=slice(
                    datetime.datetime(year=Year, month=1, day=1),
                    datetime.datetime(year=Year+1, month=1, day=1))
                )
        
        # Apply the scaling- even if our scaling is 1, it should concretize
        # the data for writing
        InDatasetSlice[OutVariable] *= Scaling

        # Map the slice to (-180, 180]
        InDatasetSlice = map_to_pm180(InDatasetSlice)

        # Add commentary on this processing
        InDatasetSlice.attrs["history"] += f'{datetime.datetime.now}; ' +\
                'separated into single years and modified variable names ' +\
                'for CABLE compatibility.'

        # Generate the filename to write to
        OutFileName = f'{OutTemplate}_{OutVariable}_{str(Year)}.nc'

        InDatasetSlice.to_netcdf(OutFileName)

#-----------------------------------------------------------------------------#

def convert_wind_components_to_magnitude(
        UWindDataset,
        VWindDataset,
        OutTemplate,
        YearRange=None
        ):
    """Convert the u and v wind components provided by ACCESS-ESM1.5 and
    convert them to a single wind magnitude."""

    # Open the datasets
    UWindDataset = xarray.open_mfdataset(
            UWindDataset,
            chunks={'time': 365*8},     # 3 hourly timesteps, a year of data
            engine = 'h5netcdf'
            )

    VWindDataset = xarray.open_mfdataset(
            VWindDataset,
            chunks={'time': 365*8},     # 3 hourly timesteps, a year of data
            engine = 'h5netcdf'
            )

    # Determine the time domain of the data, if not provided
    if YearRange is not None:
        StartTimeDomain = YearRange[0]
        EndTimeDomain = YearRange[1]
    else:
        StartTimeDomain = UWindDataset.time[0].dt.year.item()
        EndTimeDomain = UWindDataset.time[-1].dt.year.item()

    # The ACCESS-ESM data uses 1850 as the reference year, and when converting
    # the time since reference to seconds, we actually hit the size limit for
    # signed integers. To temporarily get around this, set the reference year
    # to the start year of the selected slice
    StartTimeAsSec = f'seconds since {StartTimeDomain}-01-01'

    # Set the units encoding
    UWindDataset.time.encoding['units'] = StartTimeAsSec
    if 'time_bnds' in UWindDataset.variables:
        UWindDataset.time_bnds.encoding['units'] = StartTimeAsSec

    VWindDataset.time.encoding['units'] = StartTimeAsSec
    if 'time_bnds' in VWindDataset.variables:
        VWindDataset.time_bnds.encoding['units'] = StartTimeAsSec

    # Set up the longitudes/latitudes to interpolate onto
    Longitudes = numpy.linspace(
            0.0,
            360.0,
            len(UWindDataset['lon']),
            endpoint=False
            )

    Latitudes = numpy.linspace(-90.0, 90.0, len(UWindDataset['lat']))
    
    # Create the template dataset to interpolate onto with xESMF
    MagWindTemplate = xarray.Dataset(
            coords=dict(
                lat=('lat', Latitudes),
                lon=('lon', Longitudes)
                ),
            )

    # Build the regridders using xESMF
    UWindRegridder = xesmf.Regridder(
            UWindDataset,
            MagWindTemplate,
            'bilinear',
            periodic=True
            )
    VWindRegridder = xesmf.Regridder(
            VWindDataset,
            MagWindTemplate,
            'bilinear',
            periodic=True
            )

    # Do the operations by year
    for Year in range(StartTimeDomain, EndTimeDomain+1):
        UWindSlice = UWindDataset.sel(
                time=slice(
                    datetime.datetime(year=Year, month=1, day=1),
                    datetime.datetime(year=Year+1, month=1, day=1)
                    ),
                drop=True
                )

        VWindSlice = VWindDataset.sel(
                time=slice(
                    datetime.datetime(year=Year, month=1, day=1),
                    datetime.datetime(year=Year+1, month=1, day=1)
                    ),
                drop=True
                )

        MagWindSlice = numpy.sqrt(
                UWindRegridder(UWindSlice)['uas'] ** 2 +\
                VWindRegridder(VWindSlice)['vas'] ** 2
                )

        # This is now a DataArray- bind it to a Dataset
        # Should it inherit all the metadata from the original U and V wind
        # datasets, or an abridged dataset referring to the original data?
        MagWindSlice = xarray.Dataset(
                data_vars=dict(wind_speed=MagWindSlice),
                coords=MagWindSlice.coords,
                attrs=UWindDataset.attrs
                )

        # Append to the history attribute?
        MagWindSlice.attrs['history'] += ' Wind speed computed by taking the ' +\
                'Euclidean norm of the the eastward and northwind wind ' +\
                'components, each regridded onto a grid common with the ' +\
                'remaining variables using bilinear interpolation. Mapped ' +\
                'from [0.0, 360.0) coordinates to (-180.0, 180.0].'

        # We need to set all the variable attributes manually
        MagWindSlice['wind_speed'].attrs = {
                'standard_name' : 'wind_speed',
                'long_name'     : 'Near-Surface Wind Speed',
                'comment'       : 'Magnitude of the Euclidean norm of the ' +\
                        'east-ward and north-ward wind components at 10m.',
                'cell_methods'  : 'time: point'
                }

        # Map the slice to (-180, 180]
        MagWindSlice = map_to_pm180(MagWindSlice)

        # Generate the filename to write to
        OutFileName = f'{OutTemplate}_windspeed_{str(Year)}.nc'

        MagWindSlice.to_netcdf(OutFileName)
       
#-----------------------------------------------------------------------------#

if __name__ == '__main__':
    args = read_CLI_args()

    # Define the mapping from ACCESS-ESM names to CABLE
    VariablesToConvert = {
            'rainf' : 'pr',
            'snow'  : 'prsn',
            'lw'    : 'rlds',
            'sw'    : 'rsds',
            'ps'    : 'ps',
            'qa'    : 'huss',
            'ta'    : 'tas',
            }
    
    # Split the year range, if given
    if args.years is not None:
        YearRange = tuple(int(Year) for Year in args.years.split(','))
    else:
        YearRange = None

    # A convenience lambda for file building
    InFileName = lambda x: f'/g/data/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/{args.experiment}/{args.perturbation}/3hr/{x}/gn/latest/{x}_3hr_ACCESS-ESM1-5_{args.experiment}_{args.perturbation}_gn_*.nc'

    # Get number of CPUs available
    NCPU = multiprocessing.cpu_count()

    # Call the converter for each variable in parallel
    with multiprocessing.Pool(NCPU) as Proc:
        CABLEVar = VariablesToConvert.keys()
        ACCESSVar = VariablesToConvert.values()
        InDatasets = [InFileName(ACCVar) for ACCVar in ACCESSVar]

        # Need to convert the set of inputs to a single iterable
        ConvertInputs = tuple((InData,
                               ACCVar,
                               args.output,
                               CABVar,
                               1,
                               YearRange)
            for InData, ACCVar, CABVar in zip(InDatasets, ACCESSVar, CABLEVar))

        Proc.starmap(convert_variable, ConvertInputs)

    # Handle the wind
    UWindDataset = InFileName('uas')
    VWindDataset = InFileName('vas')

    convert_wind_components_to_magnitude(
            UWindDataset,
            VWindDataset,
            args.output,
            YearRange=YearRange
            )
