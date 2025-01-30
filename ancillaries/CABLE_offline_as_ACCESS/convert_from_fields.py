import mule
import xarray
import numpy
import argparse

#----------------------------------------------------------------------------#    

def read_CLI_args():
    """Read the command line arguments."""

    # Set up the argparser with default values
    parser = argparse.ArgumentParser(
            prog="convert_from_fields",
            description="Create a CABLE gridinfo file from ACCESS-ESM1.5 " +\
            "restart file.",
            )

    parser.add_argument(
            "-i",
            "--input",
            help="ACCESS-ESM1.5 restart file to use as input",
            default="/g/data/vk83/configurations/inputs/access-esm1p5/modern/pre-industrial/restart/atmosphere/PI-02.astart-01010101"
            )

    parser.add_argument(
            "-s",
            "--stash",
            help="STASHmaster(s) to attach to the fields file",
            default="/g/data/access/umdir/vn7.3/ctldata/STASHmaster/STASHmaster_A,/g/data/rp23/experiments/2024-03-12_CABLE4-dev/lw5085/CABLE-as-ACCESS/prefix.PRESM_A"
            )

    parser.add_argument(
            "-a",
            "--areafile",
            help="NetCDF file containing areacella variable to use",
            default="/g/data/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/fx/areacella/gn/v20191115/areacella_fx_ACCESS-ESM1-5_historical_r1i1p1f1_gn.nc"
            )

    parser.add_argument(
            "-o",
            "--output",
            help="File to write the CABLE gridinfo to",
            default="ACCESS-ESM1p5-1p875x1p25-gridinfo-CABLE.nc"
            )

    parser.add_argument(
            "--create_landmask",
            action='store_true',
            help="Create a landmask from the given ESM restart file"
            )

    parser.add_argument(
            "--landmask_file",
            help="File to save the land mask to",
            default="ACCESS-ESM1p5-1p875x1p25-landmask.nc"
            )

    args = parser.parse_args()

    return args

#----------------------------------------------------------------------------#    

def define_ncvar(Data, NCDataset, VarName, dims, attribs, RollDist):
    """Add a variable to the given netCDF file, with the given data, name,
    dims and attributes."""

    # Set the fill value, depending on the type of the data
    FillVal = -1.0 if Data.dtype.kind == "f" else -1
    DType = "float32" if Data.dtype.kind == "f" else "int32"
    numpy.ma.set_fill_value(Data, FillVal)

    # Roll the array so that it matches with monotonic longitudes
    RollAxis = dims.index("longitude")
    Data = numpy.roll(Data, RollDist, axis=RollAxis)

    # Create the variable
    NCDataset[VarName] = (dims, Data)

    # Set the encoding
    NCDataset[VarName].encoding = {
            "dtype": DType,
            "_FillValue": FillVal
            }

    # Assign the attributes
    NCDataset[VarName] = NCDataset[VarName].assign_attrs(attribs)

#----------------------------------------------------------------------------#    

def get_area(AreaFile, NCDataset, dims, attribs, RollDist):
    """Get the area attribute from a different netCDF file."""

    # Assuming we have areacella attribute
    Area = AreaFile["areacella"].to_numpy()

    define_ncvar(
            Area,
            NCDataset,
            "area",
            dims["area"],
            attribs["area"],
            RollDist
            )
    
#----------------------------------------------------------------------------#    

def retrieve_landmask(FieldsFile, RollDist):
    """Retrieve the landmask from the fields file, to use with other
    variables."""

    # Get the land mask from the stash
    MaskStash = FieldsFile.stashmaster.by_regex("LAND MASK")

    # Should be a single item
    StashCode = list(MaskStash.values())[0].item

    for Field in FieldsFile.fields:
        if Field.lbuser4 == StashCode:
            # Somehow there are two fields with the same stash code? I thought
            # this was supposed to be impossible
            Landmask = Field.get_data()
            # For some reason, this single level field has duplicate stash
            # codes, but the second entry is full of 0s. Only take the first
            break

    # Convert it to boolean (it seems to be a float initially?), and invert it,
    # since we want to use it where numpy.masked_where
    Landmask = numpy.invert(Landmask.astype(bool))

    return Landmask

#----------------------------------------------------------------------------#    

def compute_iveg_and_dep_vars(
        FieldsFile,
        NCDataset,
        dims,
        attribs,
        Mask,
        RollDist
        ):
    """Determine the dominant vegetation type and use it to get the LAI,
    soil moisture, soil temperature, snow depth and soil type."""

    # Retrieve the iveg fraction
    PFTStash = FieldsFile.stashmaster.by_regex("FRACTIONS OF SURFACE TYPES")

    # Get the stash code from what should be a single item list
    StashCode = list(PFTStash.values())[0].item

    # We"re going to need to manually iterate through indices
    nLat = NCDataset.sizes["latitude"]
    nLon = NCDataset.sizes["longitude"]

    # Yank out the arrays corresponding to the land fractions first and place
    # them in a single array
    PFTArrays = numpy.zeros((17, nLat, nLon), dtype = numpy.float32)

    for Field in FieldsFile.fields:
        if Field.lbuser4 == StashCode:
            PFTArrays[Field.lbuser5-1, :, :] = Field.get_data()

    # Now we want the max index across each page of the array (+1 due to 0
    # based indexing).
    PFTs = numpy.argmax(PFTArrays, axis=0) + 1

    # Now apply the mask, based on the sum of the PFT fractions
    PFTs = numpy.ma.masked_array(PFTs, mask=Mask)

    # Send it to the ncvar creator
    define_ncvar(
            PFTs,
            NCDataset,
            "iveg",
            dims["iveg"],
            attribs["iveg"],
            RollDist
            )

    # Set all active cells with a land fraction of 1.0
    PatchFrac = numpy.ma.masked_array(
            numpy.ones((nLat, nLon), dtype=numpy.float32),
            mask=Mask
            )[numpy.newaxis, :, :]

    define_ncvar(
            PatchFrac,
            NCDataset,
            "patchfrac",
            dims["patchfrac"],
            attribs["patchfrac"],
            RollDist
            )

    # Now use this information to get the attributes which are dependent on
    # the PFT- this is LAI and isoil (with rhosoil and css being dependent on
    # isoil).
    
    # Soil Moisture (all layers)
    # Should retrieve stash entries for layer 1, 2, 3, 4, 5, 6
    SoilMoistureStash = FieldsFile.stashmaster.by_regex("SOIL MOISTURE LAYER")

    # We know the desired size of the system from the dims- initialise it with
    # values of -1.0 for easier masking
    SoilMoisture = numpy.ones(
            tuple(NCDataset.sizes[dim] for dim in dims["SoilMoist"]),
            dtype = numpy.float32
            ) * -1.0
    
    # Thickness of the soil layers
    SoilThickness = [0.022, 0.058, 0.154, 0.409, 1.085, 6.872]

    # Density of water
    RhoWater = 1000.0

    # Iterate through the soil layers
    for (Layer, SoilLayer) in enumerate(SoilMoistureStash.values()):
        StashCode = SoilLayer.item

        # Get the set of soil moisture arrays, so we can retrieve the value
        # from the desired PFT
        SoilMoistureArrays = []
        for Field in FieldsFile.fields:
            if Field.lbuser4 == StashCode:
                SoilMoistureArrays.append(Field.get_data())


        # Now iterate through the valid indices of the PFT array
        for (i, j), PFT in numpy.ma.ndenumerate(PFTs):
            SoilMoisture[:, Layer, i, j] = SoilMoistureArrays[PFT-1][i, j]


        # Scale the layer by the thickness and water density
        SoilMoisture[:, Layer, :, :] /= (SoilThickness[Layer] * RhoWater)

    # Apply the mask- anywhere the value is less than 0
    SoilMoisture = numpy.ma.masked_less(SoilMoisture, 0.0)

    define_ncvar(
            SoilMoisture,
            NCDataset,
            "SoilMoist",
            dims["SoilMoist"],
            attribs["SoilMoist"],
            RollDist
            )

    # Now do soil temperature, in the same manner as soil moisture
    # Should retrieve stash entries for layer 1, 2, 3, 4, 5, 6
    SoilTempStash = FieldsFile.stashmaster.by_regex("SOIL TEMPERATURE LAYER")

    # We know the desired size of the system from the dims- initialise it with
    # values of -1.0 for easier masking
    SoilTemp = numpy.ones(
            tuple(NCDataset.sizes[dim] for dim in dims["SoilTemp"]),
            dtype = numpy.float32
            ) * -1.0
    
    # Iterate through the soil layers
    for (Layer, SoilLayer) in enumerate(SoilTempStash.values()):
        StashCode = SoilLayer.item

        # Get the set of soil moisture arrays, so we can retrieve the value
        # from the desired PFT
        SoilTempArrays = []
        for Field in FieldsFile.fields:
            if Field.lbuser4 == StashCode:
                SoilTempArrays.append(Field.get_data())


        # Now iterate through the valid indices of the PFT array
        for (i, j), PFT in numpy.ma.ndenumerate(PFTs):
            SoilTemp[:, Layer, i, j] = SoilTempArrays[PFT-1][i, j]

    # Apply the mask- anywhere where the value is less than 0 must be a masked
    # point
    SoilTemp = numpy.ma.masked_less(SoilTemp, 0.0)

    define_ncvar(
            SoilTemp,
            NCDataset,
            "SoilTemp",
            dims["SoilTemp"],
            attribs["SoilTemp"],
            RollDist
            )

    # Now do snow depth- ACCESS-ESM1.5 has 3 snow layers, but CABLE only has 1.
    # Retrieve the snow depth stash items
    # WARNING: The STASH file has the snow fields incorrectly labelled, as of
    # 17/12/2024. The fields with Stash code 819,820,821 are labelled as
    # SNOW TEMPERATURE LAYER 1,2,3 but they are actually SNOW DEPTH LAYER 1,2,3
    # Should get layer 1, 2, 3
    SnowDepthStash = FieldsFile.stashmaster.by_regex("SNOW TEMPERATURE LAYER")

    # Initialise the array- don"t initialise to -1.0, as we want to summate
    # over layers
    SnowDepth = numpy.zeros(
            tuple(NCDataset.sizes[dim] for dim in dims["SnowDepth"]),
            dtype = numpy.float32
            )

    for (Layer, SnowLayer) in enumerate(SnowDepthStash.values()):
        StashCode = SnowLayer.item

        # Get the arrays corresponding to each PFT
        SnowLayerArrays = []
        for Field in FieldsFile.fields:
            if Field.lbuser4 == StashCode:
                SnowLayerArrays.append(Field.get_data())

        # Iterate through the valid indices of the PFT array, summating the
        # snow depth
        for (i, j), PFT in numpy.ma.ndenumerate(PFTs):
            SnowDepth[:, i, j] += SnowLayerArrays[PFT-1][i, j]

    # Apply the mask- since we didn"t initialise the array to -1.0, we can"t
    # just use numpy.ma.masked_less. Inherit the mask from the PFT array, and
    # broadcast it to the desired shape by repeating over the time axis.
    SnowDepth = numpy.ma.masked_array(
            SnowDepth,
            mask=numpy.repeat(Mask[numpy.newaxis, :, :], 12, axis=0)
            )

    define_ncvar(
            SnowDepth,
            NCDataset,
            "SnowDepth",
            dims["SnowDepth"],
            attribs["SnowDepth"],
            RollDist
            )

    # LAI
    # Retrieve the LAI stash item
    LAIStash = FieldsFile.stashmaster.by_regex("LAI")

    # Should be a single entry- retrieve the stash code
    StashCode = list(LAIStash.values())[0].item

    # Initialise the array of LAI values
    LAI = numpy.zeros(
            tuple(NCDataset.sizes[dim] for dim in dims["LAI"]),
            dtype = numpy.float32
            )

    # Here we just place the LAI arrays in a list, so we can retrieve the
    # relevant LAI from the desired array
    LAIArrays = []
    for Field in FieldsFile.fields:
        if Field.lbuser4 == StashCode:
            LAIArrays.append(Field.get_data())

    # Get the LAI from the dominant PFT (-1 due to 0 based)
    for (i, j), PFT in numpy.ma.ndenumerate(PFTs):
        LAI[:, i, j] = LAIArrays[PFT-1][i, j]

    # Apply the mask
    LAI = numpy.ma.masked_array(
            LAI,
            mask=numpy.repeat(Mask[numpy.newaxis, :, :], 12, axis=0)
            )

    define_ncvar(
            LAI,
            NCDataset,
            "LAI",
            dims["LAI"],
            attribs["LAI"],
            RollDist
            )

    # Now do isoil- 2 everywhere except where PFT=19, where isoil=9
    iSoil = numpy.ma.where(PFTs == 17, 9, 2)

    define_ncvar(
            iSoil,
            NCDataset,
            "isoil",
            dims["isoil"],
            attribs["isoil"],
            RollDist
            )

    # Variables that depend on isoil
    # rhosoil
    rhoSoil = numpy.ma.where(iSoil == 2, 1600.0, 900.0)

    define_ncvar(
            rhoSoil,
            NCDataset,
            "rhosoil",
            dims["rhosoil"],
            attribs["rhosoil"],
            RollDist
            )

    # css
    css = numpy.ma.where(iSoil == 2, 850.0, 2100.0)

    define_ncvar(
            css,
            NCDataset,
            "css",
            dims["css"],
            attribs["css"],
            RollDist)

    # Give back the mask to pass to other builders
    return Mask

#----------------------------------------------------------------------------#    

def write_albedo(FieldsFile, NCDataset, dims, attribs, Mask, RollDist):
    """Take the Albedo from the UM file as is."""
    
    # Retrieve the Albedo stash item
    AlbedoStash = FieldsFile.stashmaster.by_regex("SNOW-FREE ALBEDO OF SOIL")

    # Should be a single item dictionary- retrieve the value
    StashCode = list(AlbedoStash.values())[0].item

    # Retrieve the array
    for Field in FieldsFile.fields:
        if Field.lbuser4 == StashCode:
            Albedo = Field.get_data()

    # Mask where the values are less than 0
    Albedo2 = numpy.ma.masked_array(Albedo, mask=Mask)

    # Send it to the ncvar creator
    define_ncvar(Albedo2,
                 NCDataset,
                 "albedo2",
                 dims["albedo2"],
                 attribs["albedo2"],
                 RollDist)

    # We also need the albedo variable, although it"s not used
    # Use the template of Albedo2 for the mask
    Albedo = numpy.repeat(
            (numpy.ma.ones_like(Albedo2) * 0.2)[numpy.newaxis, :, :],
            3,
            axis=0
            )

    define_ncvar(Albedo,
                 NCDataset,
                 "Albedo",
                 dims["Albedo"],
                 attribs["Albedo"],
                 RollDist
                 )

#----------------------------------------------------------------------------#    

def convert_soilorder_to_int(
        FieldsFile,
        NCDataset,
        dims,
        attribs,
        Mask,
        RollDist
        ):
    """Convert the soil order from floating point to integer."""

    # Get the stash entry associated
    SoilOrderStash = FieldsFile.stashmaster.by_regex("SOIL ORDER")
    StashCode = list(SoilOrderStash.values())[0].item

    # Get the field associated with the stash code
    for Field in FieldsFile.fields:
        if Field.lbuser4 == StashCode:
            break

    # Convert it to int
    SoilOrder = Field.get_data().astype(numpy.int32)

    # Mask the array
    SoilOrder = numpy.ma.masked_array(SoilOrder, mask=Mask)

    # Pass it to the variable creator
    define_ncvar(
            SoilOrder,
            NCDataset,
            "SoilOrder",
            dims["SoilOrder"],
            attribs["SoilOrder"],
            RollDist
            )

#----------------------------------------------------------------------------#    

def direct_conversions(FieldsFile, NCDataset, dims, attribs, Mask, RollDist):
    """Convert all variables that are a 1:1 conversion, with no level
    dimension, from a UM field to a CABLE field."""

    # Build a dictionary of conversion mappings
    mappings = {
            "NITROGEN FIXATION": "Nfix",
            "PHOSPHORUS FROM WEATHERING": "Pwea",
            "PHOSPHORUS FROM DUST": "Pdust",
            "Dust parent soil clay fraction": "clay",
            "Dust parent soil silt fraction": "silt",
            "Dust parent soil sand fraction": "sand",
            "VOL SMC AT WILTING AFTER TIMESTEP": "swilt",
            "VOL SMC AT CRIT PT AFTER TIMESTEP": "sfc",
            "VOL SMC AT SATURATION AFTER TIMESTEP": "ssat",
            "CLAPP-HORNBERGER \"B\" COEFFICIENT": "bch",
            "SATURATED SOIL WATER SUCTION": "sucs",
            "THERMAL CAPACITY AFTER TIMESTEP": "cnsd"
            }

    # Now iterate through the mappings and add the NCvars for each
    for UMVar, CABLEVar in mappings.items():
        # Get the stash code
        VarStash = FieldsFile.stashmaster.by_regex(UMVar)
        StashCode = list(VarStash.values())[0].item

        # Get the field from the stash code
        for Field in FieldsFile.fields:
            if Field.lbuser4 == StashCode:
                break

        # Mask the values in a reasonable range- UM seems to set missing
        # float values to 1e9?
        Data = numpy.ma.masked_array(Field.get_data(), mask=Mask)

        define_ncvar(
                Data,
                NCDataset,
                CABLEVar,
                dims[CABLEVar],
                attribs[CABLEVar],
                RollDist
                )

#----------------------------------------------------------------------------#    

def scaling_conversions(
        FieldsFile,
        NCDataset,
        dims,
        attribs,
        scalings,
        Mask,
        RollDist
        ):
    """Perform all conversions that require scalar/unit conversions to go from
    UM fields to CABLE fields."""

    # Current variables that fit this description are Ndep and hyds
    mappings = {
            "NITROGEN DEPOSITION": "Ndep",
            "SAT SOIL CONDUCTIVITY AFTER TIMESTEP": "hyds"
            }

    for UMVar, CABLEVar in mappings.items():
        print(f"Converting {UMVar}")
        # Get the stash code
        VarStash = FieldsFile.stashmaster.by_regex(UMVar)
        StashCode = list(VarStash.values())[0].item

        # Get the field from the stash code
        for Field in FieldsFile.fields:
            if Field.lbuser4 == StashCode:
                break

        Data = numpy.ma.masked_array(
                Field.get_data() * scalings[CABLEVar],
                mask=Mask
                )

        define_ncvar(
                Data,
                NCDataset,
                CABLEVar,
                dims[CABLEVar],
                attribs[CABLEVar],
                RollDist
                )

def apply_ice_fixes(NCDataset):
    """Apply the hard-coded ice fixes, since the ACCESS-ESM1.5 restart file
    does not contain the values actually used by the land model. See
    https://github.com/CABLE-LSM/CABLE/pull/497#issuecomment-2589039055"""

    # We want to use the DataArray.where method to keep everything not iveg=17
    # as is, and set everywhere else to a hardcoded value
    IceMask = NCDataset['iveg'] != 17

    # Use the mask to modify select values for each parameter
    NCDataset['isoil'] = NCDataset['isoil'].where(IceMask, 9)
    NCDataset['clay'] = NCDataset['clay'].where(IceMask, 0.3)
    NCDataset['silt'] = NCDataset['silt'].where(IceMask, 0.33)
    NCDataset['sand'] = NCDataset['sand'].where(IceMask, 0.37)
    NCDataset['swilt'] = NCDataset['swilt'].where(IceMask, 0.216)
    NCDataset['sfc'] = NCDataset['sfc'].where(IceMask, 0.301)
    NCDataset['ssat'] = NCDataset['ssat'].where(IceMask, 0.479)
    NCDataset['bch'] = NCDataset['bch'].where(IceMask, 7.1)
    NCDataset['hyds'] = NCDataset['hyds'].where(IceMask, 1e-6)
    NCDataset['sucs'] = NCDataset['sucs'].where(IceMask, -0.153)
    NCDataset['rhosoil'] = NCDataset['rhosoil'].where(IceMask, 910.0)
    NCDataset['css'] = NCDataset['css'].where(IceMask, 2100.0)

    return NCDataset

#----------------------------------------------------------------------------#    

if __name__ == "__main__":

    # Read the command line args
    args = read_CLI_args()

    # Start by opening the desired UM fields file and attaching the stash
    FieldsFile = mule.FieldsFile.from_file(args.input)

    # We currently build the ESM stash as a merge of the base stash and an
    # extension. The args.stash may contain multiple comma separated stashes
    ESMStashes = [mule.STASHmaster.from_file(SM) for SM in
                args.stash.split(',')]
    
    # Merge the stashes into 1
    ESMStash = mule.STASHmaster()
    for Stash in ESMStashes:
        ESMStash.update(Stash)

    # Only take the first section- others are for diagnostic variables
    ESMStash = ESMStash.by_section(0)
    FieldsFile.attach_stashmaster_info(ESMStash)

    AreaFile = xarray.open_dataset(
            args.areafile,
            engine="netcdf4"
            )

    # Define the dimensions of the NetCDF file
    NCDimensions = {
            "time": 12,
            "rad": 3,
            "soil": 6,
            "patch": 1,
            "latitude": 145,
            "longitude": 192
            }

    # Remember to map longitude from 0.0 -> 360.0 to -180.0 -> 180.0
    Longitudes = numpy.linspace(
            0.0,
            360.0,
            NCDimensions["longitude"],
            endpoint=False)

    Longitudes = numpy.where(Longitudes > 180.0, -360 + Longitudes, Longitudes)
    
    # We also want to roll the longitude axis so that it is monotonically
    # ascending- and also store the amount we want to shift, to use on the
    # data arrays
    RollDist = NCDimensions["longitude"] - numpy.argmax(Longitudes < 0.0)
    Longitudes = numpy.roll(Longitudes, RollDist)
    Latitudes = numpy.linspace(-90.0, 90.0, NCDimensions["latitude"])

    # For some reason xarray doesn"t accept coordinates without values in the
    # constuctor (what an awesome package), need to fill the other dimensions
    # with something
    rad = list(range(3))
    time = list(range(12))
    patch = list(range(1))
    soil = list(range(6))

    NCCoords = {
            "longitude": (["longitude"], Longitudes),
            "latitude": (["latitude"], Latitudes),
            "patch": (["patch"], patch),
            "soil": (["soil"], soil),
            "rad": (["rad"], rad),
            "time": (["time"], time)
            }

    # Add top level attributes
    NCAttribs = {
            "source": "ACCESS-ESM1.5 atmosphere restart file from " +\
                    "CMIP6 pre-industrial runs.",
            "comments": "Converted from the UM fields file on Gadi at " +\
                    "/g/data/vk83/configurations/inputs/access-esm1p5/modern/pre-industrial/restart/atmosphere/PI-02.astart-01010101 " +\
                    "with the python script at ." +\
                    "https://github.com/ACCESS-Community-Hub/Land-ancillary-creation/ancillaries/CABLE_offline_as_ACCESS",
            "contact": "lachlan.whyborn@anu.edu.au"
            }

    # Apparently you can"t give dimensions to a Dataset?
    NCDataset = xarray.Dataset(
            coords = NCCoords,
            attrs = NCAttribs
            )

    # Define a single dictionary which contains the dimension information for
    # all the variables
    VarDimensions = {
            "area"      : ("latitude", "longitude"),
            "iveg"      : ("latitude", "longitude"),
            "patchfrac" : ("patch", "latitude", "longitude"),
            "isoil"     : ("latitude", "longitude"),
            "SoilMoist" : ("time", "soil", "latitude", "longitude"),
            "SoilTemp"  : ("time", "soil", "latitude", "longitude"),
            "SnowDepth" : ("time", "latitude", "longitude"),
            "Albedo"    : ("rad", "latitude", "longitude"),
            "LAI"       : ("time", "latitude", "longitude"),
            "SoilOrder" : ("latitude", "longitude"),
            "Ndep"      : ("latitude", "longitude"),
            "Nfix"      : ("latitude", "longitude"),
            "Pwea"      : ("latitude", "longitude"),
            "Pdust"     : ("latitude", "longitude"),
            "clay"      : ("latitude", "longitude"),
            "silt"      : ("latitude", "longitude"),
            "sand"      : ("latitude", "longitude"),
            "swilt"     : ("latitude", "longitude"),
            "sfc"       : ("latitude", "longitude"),
            "ssat"      : ("latitude", "longitude"),
            "bch"       : ("latitude", "longitude"),
            "hyds"      : ("latitude", "longitude"),
            "sucs"      : ("latitude", "longitude"),
            "rhosoil"   : ("latitude", "longitude"),
            "cnsd"      : ("latitude", "longitude"),
            "css"       : ("latitude", "longitude"),
            "albedo2"   : ("latitude", "longitude"),
            "mask"      : ("latitude", "longitude")
            }

    # And do the same with the attributes
    VarAttribs = {
            "longitude" : {
                "units"         : "degrees east"
                },
            "latitude"  : {
                "units"         : "degrees north"
                },
            "area"      : {
                "long_name"     : "Area of grid cell",
                "units"         : "m^2",
                "missing_value" : "-1.0"
                },
            "iveg"      : {
                "long_name"     : "CSIRO classification of vegetation type",
                "units"         : 1,
                "flag_values"   : [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15,
                                   16, 17],
                "flag_meanings" : "Evergreen_Needleleaf Evergreen_Broadleaf "+\
                                  "Deciduous Needleleaf Deciduous Broadleaf "+\
                                  "Shrub C3_Grassland C4_Grassland Tundra "+\
                                  "C3_Cropland C4_Cropland Wetland Barren "+\
                                  "Urban Lakes Ice",
                "missing_value": -1
                },
            "patchfrac" : {
                "long_name"     : "Fraction of grid cell assigned to " +\
                        "vegetation type",
                "units"         : 1,
                "missing_value" : -1.0
                },
            "isoil"     : {
                "long_name"     : "Zobler soil type",
                "units"         : 1,
                "flag_values"   : [2, 9],
                "flag_meanings" : "Unknown Unknown",
                "missing_value" : -1,
                "comments"      : "Set to 2 for all vegetation types bar " +\
                        "17, which is set to 9."
                },
            "SoilMoist" : {
                "long_name"     : "Soil moisture in a layer",
                "units"         : "m^3/m^3",
                "missing_value" : -1.0,
                "comments"      : "Converted from kg/m^2 to m^3/m^3 by " +\
                        "dividing by layer thickness and water density."
                },
            "SoilTemp"  : {
                "long_name"     : "Soil temperature in a layer",
                "units"         : "K",
                "missing_value" : -1.0
                },
            "SnowDepth" : {
                "long_name"     : "Depth of snow",
                "units"         : "m",
                "missing_value" : -1.0
                },
            "Albedo"    : {
                "long_name"     : "Snow free albedo of soil",
                "units"         : 1,
                "missing_value" : -1.0,
                "comments"      : "Set to 0.2 everywhere, as it is not " +\
                        "used in CABLE but it is checked."
                },
            "LAI"       : {
                "long_name"     : "Leaf area index",
                "units"         : 1,
                "missing_value" : -1.0
                },
            "SoilOrder" : {
                "long_name"     : "Soil order",
                "units"         : 1,
                "missing_value" : -1,
                "comments"      : "Converted to integer from ACCESS restart."
                },
            "Ndep"      : {
                "long_name"     : "annual Nitrogen deposition rate",
                "units"         : "g/m^2/year",
                "missing_value" : -1.0,
                "comments"      : "Converted from daily rate to yearly rate"
                },
            "Nfix"      : {
                "long_name"     : "annual Nitrogen fixation rate",
                "units"         : "g/m^2/year",
                "missing_value" : -1.0
                },
            "Pwea"      : {
                "long_name"     : "annual yield of Phosphorus from weathering",
                "units"         : "g/m^2/year",
                "missing_value" : -1.0
                },
            "Pdust"     : {
                "long_name"     : "annual yield of Phosphorus from dust",
                "units"         : "g/m^2/year",
                "missing_value" : -1.0
                },
            "clay"      : {
                "long_name"     : "soil clay fraction",
                "units"         : 1,
                "missing_value" : -1.0
                },
            "silt"      : {
                "long_name"     : "soil silt fraction",
                "units"         : 1,
                "missing_value" : -1.0
                },
            "sand"      : {
                "long_name"     : "soil sand fraction",
                "units"         : 1,
                "missing_value" : -1.0
                },
            "swilt"     : {
                "long_name"     : "Volume soil moisture content at wilting",
                "units"         : "m^3",
                "missing_value" : -1.0
                },
            "sfc"       : {
                "long_name"     : "Volume soil moisture content at crit point",
                "units"         : "m^3",
                "missing_value" : -1.0
                },
            "ssat"      : {
                "long_name"     : "Volume soil moisture content at saturation",
                "units"         : "m^3",
                "missing_value" : -1.0
                },
            "bch"       : {
                "long_name"     : "CLAPP-HORNBERGER coefficient",
                "units"         : 1,
                "missing_value" : -1.0
                },
            "hyds"      : {
                "long_name"     : "Saturated soil conductivity",
                "units"         : "m/s",
                "missing_value" : -1.0,
                "comments"      : "Divided by 1000 for unit conversion."
                },
            "sucs"      : {
                "long_name"     : "Saturated soil water suction",
                "units"         : "m",
                "missing_value" : -1.0
                },
            "rhosoil"   : {
                "long_name"     : "Soil density",
                "units"         : "kg/m^3",
                "missing_value" : -1.0,
                "comments"      : "Set to 1600.0 where isoil=2, " +\
                        "910.0 where isoil=9."
                },
            "cnsd"      : {
                "long_name"     : "Soil thermal conductivity",
                "units"         : "W/m/K",
                "missing_value" : -1.0
                },
            "css"       : {
                "long_name"     : "Soil specific heat capacity",
                "units"         : "J/kg/K",
                "missing_value" : -1.0,
                "comments"      : "Set to 850.0 where isoil=2, " +\
                        "2100.0 where isoil=9."
                },
            "albedo2"   : {
                "long_name"     : "Snow-free albedo of soil",
                "units"         : 1,
                "missing_value" : -1.0
                },
            "mask"      : {
                "long_name"     : "land_binary_mask",
                "comments"      : "1=land, 0=sea"
                }
            }

    scalings = {
            "Ndep"  : 365.0,
            "hyds"  : 1.0 / 1000.0
            }

    # Call the variable builders

    # Need to get area from elsewhere- there is a GRIDBOX AREA field
    # in the UM file, but for some reason it"s an array of 0.0?
    get_area(
            AreaFile,
            NCDataset,
            VarDimensions,
            VarAttribs,
            RollDist
            )

    # Retrieve the mask used to act as a numpy mask
    Mask = retrieve_landmask(FieldsFile, RollDist)

    # If required, create the landmask file
    if args.create_landmask:
        MaskCoords = {
                "longitude": (["longitude"], Longitudes),
                "latitude": (["latitude"], Latitudes)
                }
        MaskAttribs = {
                "source": "ACCESS-ESM1.5 atmosphere restart file from " +\
                        "CMIP6 pre-industrial runs.",
                "comments": "From the \"LAND MASK (No halo)\" field in " +\
                            "/g/data/vk83/configurations/inputs/access-esm1p5/modern/pre-industrial/restart/atmosphere/PI-02.astart-01010101.",
                "contact": "lachlan.whyborn@anu.edu.au"
                }

        MaskDataset = xarray.Dataset(
                data_vars = {
                    "mask": (
                        ("latitude", "longitude"),
                        numpy.roll(Mask.astype(numpy.int8), RollDist, axis=1),
                        MaskAttribs
                        )
                    },
                coords = MaskCoords,
                attrs = MaskAttribs
                )

        MaskDataset.to_netcdf(args.landmask_file)

    # We have a series of variables which are dependent on the vegetation type
    # calculation, so do them all at once
    compute_iveg_and_dep_vars(
            FieldsFile,
            NCDataset,
            VarDimensions,
            VarAttribs,
            Mask,
            RollDist
            )

    write_albedo(
            FieldsFile,
            NCDataset,
            VarDimensions,
            VarAttribs,
            Mask,
            RollDist
            )

    convert_soilorder_to_int(
            FieldsFile,
            NCDataset,
            VarDimensions,
            VarAttribs,
            Mask,
            RollDist
            )

    direct_conversions(
            FieldsFile,
            NCDataset,
            VarDimensions,
            VarAttribs,
            Mask,
            RollDist
            )

    scaling_conversions(
            FieldsFile,
            NCDataset,
            VarDimensions,
            VarAttribs,
            scalings,
            Mask,
            RollDist
            )

    NCDataset = apply_ice_fixes(NCDataset)

    NCDataset.to_netcdf(args.output)
