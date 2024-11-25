# This script used the convenience tools created in convert_netcdf.py to
# convert the NetCDF file created by xconv with the process described in
# the README to a CABLE ancillary NetCDF file.

import argparse
import netCDF4
import numpy

from convert_netcdf import add_dimensions, add_metadata, convert_variable, \
        add_variable

# Define all the mappers we will use
def convert_to_dominant(ivegFractions):
    """Convert series of vegetations to dominant type. Input data contains
    vegetation type fractions for all grid cells; this is to be converted to
    a single dominant vegetation type per grid cell, by taking the vegetation
    type with the largest fraction."""

    # Squeeze the array to remove length-1 time axis
    ivegData = ivegFractions[:].squeeze()

    # Initial dimensions are t, nveg, lat, lon
    nveg, lat, lon = ivegData.shape

    # Initialise the new array to save to
    # Get cast warnings if I try to use ma.empty_like on ivegData[0, :, :]
    ivegArray = numpy.ma.empty_like(ivegData[0, :, :], dtype = numpy.int32)
    ivegArray.set_fill_value(numpy.int32(-1))

    # Iterate through the non-masked elements
    for (lat, lon), v in numpy.ma.ndenumerate(ivegArray):
        ivegArray[lat, lon] = ivegData[:, lat, lon].argmax() + 1
    
    return ivegArray

def scale_soil_moisture(SoilMoisture):
    """Convert soil moisture from kg/m^3 to m^3/m^3."""

    # Set the layer thicknesses- known from cable_parameters.F90
    LayerThickness = [0.022, 0.058, 0.154, 0.409, 1.085, 6.872]
    # Density of water
    RhoWater = 1000.0

    # CABLE requires monthly soil moisture, with months being the first dim
    SoilMoistArray = spread_across_months(SoilMoisture)

    # Now perform the scaling
    for Layer in range(6):
        SoilMoistArray[:, Layer, :, :] *= LayerThickness[Layer] * RhoWater

    # Set fill value to be the same as CABLE
    SoilMoistArray.set_fill_value(-1.0)

    return SoilMoistArray

def convert_to_int(FloatVariable):
    """Convert dataset from floating point values to integer values."""

    # Simple int conversion
    IntArray = FloatVariable[:].astype(numpy.int32)

    # Given most integers are classifications, we should be able to set
    # the fill value to -1
    IntArray.set_fill_value(-1) 

    return IntArray

def spread_across_months(Variable):
    """Spread a given variable across 12 months along the first axis. Many
    CABLE ancillaries require monthly values (even though in most cases they
    are not used)."""
    Data = numpy.repeat(Variable[:], 12, axis = 0)
    return Data

def compute_snow_depth(InDatasets, InVariables):
    """Convert from 3 snow layers distributed across PFTs to a single snow
    depth. The original snow layers in UMrst['SnowLayer{1,2,3}'] are given per
    PFT, so we do a weighted average by scaling the snow layer depth in a PFT
    by the patch fraction in UMrst['PFTs'] then summing over the layers."""

    # We know that the dimensions of each snow layer are (z, lat, lon)
    # So we can take the last 2 dimensions to form the new snowdepth array
    Template = InDatasets["UMrst"][InVariables["UMrst"]["SnowLayer1"]]
    _, nPFT, nLat, nLon = Template.shape

    # Initialise the empty array
    SnowDepth = numpy.ma.empty_like(
            Template[0, 0, :, :],
            dtype = numpy.float32)

    # Iterate over the land points
    for PFT in range(nPFT):
        SnowDepth = InDatasets["UMrst"][InVariables["UMrst"]["SnowLayer1"]]\
                [0, PFT, :, :]
        SnowDepth += InDatasets["UMrst"][InVariables["UMrst"]["SnowLayer2"]]\
                [0, PFT, :, :]
        SnowDepth += InDatasets["UMrst"][InVariables["UMrst"]["SnowLayer3"]]\
                [0, PFT, :, :]
        SnowDepth *= InDatasets["UMrst"][InVariables["UMrst"]["PFTs"]]\
                [0, PFT, :, :]

    # Add the time axis
    SnowDepth = SnowDepth[numpy.newaxis, ...]

    # Spread over months
    SnowDepth = spread_across_months(SnowDepth)

    return SnowDepth

def add_rad_dimension(InAlbedo):
    """Take the albedo and spread it over the 3 radiation dimensions.
    Apparently this albedo dimension is not actually used by CABLE, but is
    checked for rational values?"""

    # Set it to all 0.2
    Albedo = numpy.ma.ones_like(InAlbedo[:]).__imul__(0.2)

    # Add the rad dimension
    Albedo = Albedo[numpy.newaxis, ...]
    Albedo = numpy.repeat(Albedo, 3, axis=0)

    return Albedo

def compute_LAI(InDatasets, InVariables):
    """Take the LAI from the dominant vegetation type. The LAI is described
    per PFT in the UM file, so take the the LAI from the dominant vegetation
    type in CABLErst['PFT'] to pick the LAI from the UMrst['LAI']."""

    # Use the iveg data to build a template
    iveg = InDatasets["CABLErst"][InVariables["CABLErst"]["PFT"]][:]
    LAI = numpy.ma.empty_like(iveg, dtype = numpy.float32)
    
    for (lat, lon), PFT in numpy.ma.ndenumerate(iveg):
        LAI[lat, lon] = InDatasets["UMrst"][InVariables["UMrst"]["LAI"]]\
                [0, PFT-1, lat, lon]
    
    # Spread over months
    LAI = LAI[numpy.newaxis, ...]
    LAI = spread_across_months(LAI)

    return LAI


if __name__ == "__main__":
    # Read the command line arguments

    parser = argparse.ArgumentParser(
            prog='convert_ACCESS_to_CABLE',
            description='Converts the NetCDF file created by xconv from an ' +\
                    'ACCESS-ESM1.5 fields file to a NetCDF file useable by ' +\
                    'CABLE.'
            )

    parser.add_argument('-i', '--inputfile', type=str)
    parser.add_argument('-o', '--outputfile', type=str)
    parser.add_argument('--areafile', '--areafile', type=str)

    args = parser.parse_args()

    InputFile = args.inputfile
    OutputFile = args.outputfile
    AreaFile = args.areafile

    # Invoke the conversions
    InDataset = netCDF4.Dataset(InputFile)
    AreaDataset = netCDF4.Dataset(AreaFile)
    OutDataset = netCDF4.Dataset(
            OutputFile,
            'w',
            format='NETCDF4'
            )

    DimMap = {
            'longitude': 'longitude',
            'latitude': 'latitude',
            'level6': 'soil',
            }
    # Add dimensions for months, patches and radiation bands
    ExtraDims = {'time': 12, "patch": 1, "rad": 3}

    add_dimensions(InDataset, OutDataset, DimMap, ExtraDims)

    add_metadata(OutDataset, {
        'source': 'Generated from the ACCESS-ESM1.5 restart file at ' +\
                '/g/data/vk83/configurations/inputs/access-esm1p5/modern/pre-industrial/restart/atmosphere/PI-02.astart-01010101. ' +\
                'Mapped to CABLE ancillaries via the process described at ' +\
                '[TODO github location].',
        'author': 'Lachlan Whyborn, lachlan.whyborn@anu.edu.au'
        }
    )

    convert_variable(
            InDataset,
            OutDataset,
            "longitude",
            "longitude",
            mapper = lambda x: numpy.where(x[:] < 180.0, x[:], x[:] - 360.0),
            attribs = {
                "units": "degrees_east",
                "spacing": "even"
                }
            )

    convert_variable(
            InDataset,
            OutDataset,
            "latitude",
            "latitude",
            attribs = {
                "units": "degrees_north",
                "spacing": "even"
                }
            )

    # area
    convert_variable(
            AreaDataset,
            OutDataset,
            "areacella",
            "area",
            dims = ("latitude", "longitude")
            )

    #patchfrac- just use existing data to get the landmask
    convert_variable(
            InDataset,
            OutDataset,
            "lsm",
            "patchfrac",
            mapper = lambda x: numpy.ma.ones_like(x[:]),
            dims = ("patch", "latitude", "longitude"),
            attribs = {
                "comment": "Single PFT type per grid cell, so patchfrac set " +
                        "to 1.0 everywhere.",
                "units": 1,
                }
            )

    # iveg
    convert_variable(
            InDataset,
            OutDataset,
            "field1391",
            "iveg",
            mapper = convert_to_dominant,
            dims = ("latitude", "longitude"),
            attribs = {
                "long_name": "Vegetation type classification",
                "comment": "Computed by taking the dominant vegetation " +
                    "type in each grid cell from the UM fields file.",
                "valid_min": 1,
                "valid_max": 17
                }
            )

    # isoil
    convert_variable(
            OutDataset,
            OutDataset,
            "iveg",
            "isoil",
            mapper = lambda x: numpy.ma.where(x[:] == 17, 9, 2),
            attribs = {
                "comment": "Determined based on iveg. Set to 2 for all PFTs " +
                        "bar 17 which is set to 9.",
                "units": 1,
                "name": "isoil",
                "long_name": "soil type",
                }
            )

    # Soil moisture
    convert_variable(
            InDataset,
            OutDataset,
            "sm",
            "SoilMoist",
            mapper = scale_soil_moisture,
            dims = ("time", "soil", "latitude", "longitude"),
            attribs = {
                "comment": "Rescaled from kg/m^3 to m^3/m^3 by dividing " +
                    "by the soil layer thickness and density of water. " +
                    "Single time point spread across all 12 months.",
                "units": "m^3/m^3"
                }
            )

    # Soil temperature
    convert_variable(
            InDataset,
            OutDataset,
            "soiltemp",
            "SoilTemp",
            mapper = spread_across_months,
            dims = ("time", "soil", "latitude", "longitude")
            )

    # Snow depth- need many to one here
    convert_variable(
            {"UMrst": InDataset},
            OutDataset,
            {"UMrst": {"SnowLayer1": "temp",
                       "SnowLayer2": "temp_1",
                       "SnowLayer3": "temp_2",
                       "PFTs": "field1391"}
            },
            "SnowDepth",
            dims = ("time", "latitude", "longitude"),
            mapper = compute_snow_depth,
            attribs = {
                "comment": "Computed from the UM snow layers by scaling " +
                       "the layer depth by the PFT fraction and summating " +
                       "over the layers. Spread across 12 months.",
                "units": "m"
                }
            )
                            
    # First albedo (?)
    convert_variable(
            InDataset,
            OutDataset,
            "field1395",
            "Albedo",
            dims = ("rad", "latitude", "longitude"),
            mapper = add_rad_dimension,
            attribs = {
                "units": 1,
                "long_name": "Albedo by radiation band",
                "comment": "Set to 0.2 everywhere, since apparently it is " +
                        "not used but checked?"
                }
            )

    # Soil order
    convert_variable(
            InDataset,
            OutDataset,
            "temp_3",
            "SoilOrder",
            mapper = convert_to_int,
            attribs = {
                "comment": "Converted to integer values.",
                "units": 1
                }
            )

    # LAI
    convert_variable(
            {"UMrst": InDataset, "CABLErst": OutDataset},
            OutDataset,
            {"UMrst": {"LAI": "temp_8"},
             "CABLErst": {"PFT": "iveg"}},
            "LAI",
            dims = ("time", "latitude", "longitude"),
            mapper = compute_LAI,
            attribs = {
                "comment": "Computed by taking the LAI attached to the " +
                        "dominant PFT type.",
                "units": 1,
                "name": "LAI",
                "long_name": "leaf area index",
                }
            )

    # Ndep
    convert_variable(
            InDataset,
            OutDataset,
            "temp_4",
            "Ndep",
            mapper = lambda x: x[:].__imul__(365),
            attribs = {
                "comment": "Converted from g/m^2/day to g/m^2/year.",
                "units": "g/m^2/year",
                "long_name": "annual Nitrogen deposition rate"
                }
            )

    # Nfix
    convert_variable(
            InDataset,
            OutDataset,
            "temp_5",
            "Nfix",
            attribs = {
                "units": "g/m^2/year",
                "long_name": "annual Nitrogen fixation rate"
                }
            )

    # Pwea
    convert_variable(
            InDataset,
            OutDataset,
            "temp_7",
            "Pwea",
            attribs = {
                "units": "g/m^2/year",
                "long_name": "annual yield of Phosphorous from weathering"
                }
            )

    # Pdust
    convert_variable(
            InDataset,
            OutDataset,
            "temp_6",
            "Pdust",
            attribs = {
                "units": "g/m^2/year",
                "long_name": "annual yield of Phosphorous from dust deposition"
                }
            )

    # clay
    convert_variable(
            InDataset,
            OutDataset,
            "field1630",
            "clay",
            attribs = {
                "units": 1,
                "long_name": "UM soil texture - clay fraction"
                }
            )

    # silt
    convert_variable(
            InDataset,
            OutDataset,
            "field1631",
            "silt",
            attribs = {
                "units": 1,
                "long_name": "UM soil texture - silt fraction"
                }
            )

    # clay
    convert_variable(
            InDataset,
            OutDataset,
            "field1632",
            "sand",
            attribs = {
                "units": 1,
                "long_name": "UM soil texture - sand fraction"
                }
            )

    #Swilt
    convert_variable(
            InDataset,
            OutDataset,
            "field329",
            "swilt",
            attribs = {
                "units": 1,
                "long_name": "VOL SMC AT WILTING"
                }
            )

    #sfc
    convert_variable(
            InDataset,
            OutDataset,
            "field330",
            "sfc",
            attribs = {
                "units": 1,
                "long_name": "VOL SMC AT CRIT PT"
                }
            )

    #ssat
    convert_variable(
            InDataset,
            OutDataset,
            "field332",
            "ssat",
            attribs = {
                "units": 1,
                "long_name": "VOL SMC AT SATURATION"
                }
            )

    #bch
    convert_variable(
            InDataset,
            OutDataset,
            "field1381",
            "bch",
            attribs = {
                "units": 1,
                "long_name": "VOL SMC AT SATURATION"
                }
            )

    #hyds
    convert_variable(
            InDataset,
            OutDataset,
            "field333",
            "hyds",
            mapper = lambda x: x[:] * 1e-3,
            attribs = {
                "comment": "Converted from kg/m^s/s to m/s by dividing " +
                        "by 1000",
                "units": 1
                }
            )

    #sucs
    convert_variable(
            InDataset,
            OutDataset,
            "field342",
            "sucs",
            attribs = {
                "units": "m",
                "long_name": "SATURATED SOIL WATER SUCTION",
                "missing_value": -1.0
                }
            )

    # rhosoil
    convert_variable(
            OutDataset,
            OutDataset,
            "isoil",
            "rhosoil",
            mapper = lambda x: numpy.ma.where(x[:] == 2, 1600.0, 910.0),
            attribs = {
                "comment": "Set according to isoil- 1600.0 where " + 
                        "isoil = 2 and 910.0 where isoil = 9",
                "units": "kg/m^3",
                "name": "rhosoil",
                "long_name": "density of soil",
                }
            )

    # cnsd
    convert_variable(
            InDataset,
            OutDataset,
            "field336",
            "cnsd",
            attribs = {
                "units": "W/m/K",
                "long_name": "THERMAL CONDUCTIVITY"
                }
            )

    # css
    convert_variable(
            OutDataset,
            OutDataset,
            "isoil",
            "css",
            mapper = lambda x: numpy.ma.where(x[:] == 2, 850.0, 2100.0),
            attribs = {
                "comment": "Set according to isoil- 850.0 where " +
                        "isoil = 2, 2100.0 where isoil = 9",
                "units": "J/kg/K",
                "name": "css",
                "long_name": "soil specific heat capacity",
                }
            )

    # albedo2
    convert_variable(
            InDataset,
            OutDataset,
            "field1395",
            "albedo2",
            attribs = {
                "long_name": "Snow free soil albedo",
                "units": 1
                }
            )
