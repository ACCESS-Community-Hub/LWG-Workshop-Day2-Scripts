import argparse
import datetime
import itertools
import multiprocessing
import numpy

import xarray

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

    # Start by setting up the metadata
    # The input data is in days since ..., but we convert to seconds since
    # by setting the time units attribute
    StartTime = InDataset.time.encoding['units']
    StartTimeAsSec = StartTime.replace('days', 'seconds')

    # Set the units encoding
    InDataset.time.encoding['units'] = StartTimeAsSec
    if 'time_bnds' in InDataset.variables:
        InDataset.time_bnds.encoding['units'] = StartTimeAsSec

    # Rename the variable to one recogniseble by GSWP routines and set attrs
    InDataset = InDataset.rename(
            {InVariable: OutVariable}
            )
    
    # Determine the time domain of the data, if not provided
    if YearRange is not None:
        StartTimeDomain = YearRange[0]
        EndTimeDomain = YearRange[1]
    else:
        StartTimeDomain = InDataset.time[0].dt.year
        EndTimeDomain = InDataset.time[-1].dt.year

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
        StartTimeDomain = UWindDataset.time[0].dt.year
        EndTimeDomain = UWindDataset.time[-1].dt.year

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

        # Convert to numpy arrays, then compute magnitude
        MagWindSlice = numpy.sqrt(
                UWindSlice['uas'].to_numpy()**2 +\
                VWindSlice['vas'].to_numpy()**2
                )

        # To make handling of the metadata simpler, we'll just write over the 
        # UWindDataset component
        UWindSlice['uas'] = MagWindSlice

        # Now rename the variable, and make appropriate modifications to the
        # metadata
        MagWindSlice = UWindSlice.rename({'uas': 'wind'})

        # Set master attrs
        MagWindSlice.attrs['variable_id'] = 'wind_speed'

        # Set variable attrs
        MagWindSlice['wind'].attrs['standard_name'] = 'wind_speed'
        MagWindSlice['wind'].attrs['long_name'] = 'near surface wind speed'
        MagWindSlice['wind'].attrs['history'] = 'Computed from eastward ' +\
                'and northword wind components.'

        MagWindSlice['wind'].encoding['units'] = \
            MagWindSlice['wind'].encoding['units'].replace('days', 'seconds')

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

    # for CABLEVar, ACCESSVar in VariablesToConvert.items():
        # # Construct the input dataset location
        # InDataset = InFileName(ACCESSVar)

        # convert_variable(
                # InDataset,
                # ACCESSVar,
                # args.output,
                # CABLEVar,
                # YearRange=YearRange
                # )

    # # Handle the wind
    # UWindDataset = InFileName('uas')
    # VWindDataset = InFileName('vas')

    # convert_wind_components_to_magnitude(
            # UWindDataset,
            # VWindDataset,
            # args.output,
            # YearRange=YearRange
            # )
