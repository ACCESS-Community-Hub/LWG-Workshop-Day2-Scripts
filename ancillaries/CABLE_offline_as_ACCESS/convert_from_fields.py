import mule
import netCDF4
import numpy

#----------------------------------------------------------------------------#    

def define_ncvar(Data, NCFile, VarName, dims, attribs):
    """Add a variable to the given netCDF file, with the given data, name,
    dims and attributes."""

    # Set the fill value, depending on the type of the data
    FillVal = -1.0 if Data.dtype.kind == 'f' else -1

    # Create the variable
    NCFile.createVariable(
            VarName,
            Data.dtype,
            dims,
            fill_value = FillVal
            )

    # Assign the data
    NCFile[VarName][:] = Data

    # Assign the attributes
    for (attrib, value) in attribs.items():
        NCFile[VarName].setncattr(attrib, value)

#----------------------------------------------------------------------------#    

def write_albedo(FieldsFile, NCFile, dims, attribs):
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
    Albedo2 = numpy.ma.masked_less(Albedo, 0.0)

    # Send it to the ncvar creator
    define_ncvar(Albedo2,
                 NCFile,
                 "albedo2",
                 dims["albedo2"],
                 attribs["albedo2"])

    # We also need the albedo variable, although it's not used
    # Use the template of Albedo2 for the mask
    Albedo = numpy.repeat(
            (numpy.ma.ones_like(Albedo2) * 0.2)[numpy.newaxis, :, :],
            len(NCFile.dimensions["rad"]),
            axis=0
            )
    define_ncvar(Albedo,
                 NCFile,
                 "Albedo",
                 dims["Albedo"],
                 attribs["Albedo"]
                 )

#----------------------------------------------------------------------------#    

def compute_iveg_and_dependent_vars(FieldsFile, NCFile, dims, attribs):
    """Determine the dominant vegetation type and use it to get the LAI,
    soil moisture, soil temperature, snow depth and soil type."""

    # Retrieve the iveg fraction
    PFTStash = FieldsFile.stashmaster.by_regex("FRACTIONS OF SURFACE TYPES")

    # Get the stash code from what should be a single item list
    StashCode = list(PFTStash.values())[0].item

    # We're going to need to manually iterate through indices
    nLat = len(NCFile.dimensions["latitude"])
    nLon = len(NCFile.dimensions["longitude"])

    # Yank out the arrays corresponding to the land fractions first and place
    # them in a single array
    PFTArrays = numpy.zeros((17, nLat, nLon), dtype = numpy.float32)

    for Field in FieldsFile.fields:
        if Field.lbuser4 == StashCode:
            PFTArrays[Field.lbuser5-1, :, :] = Field.get_data()

    # Build the mask that we will use- in most instances we'll also want to
    # broadcast it over another dimension.
    Mask = numpy.sum(PFTArrays, axis=0) <= 0.0

    # Now we want the max index across each page of the array (+1 due to 0
    # based indexing).
    PFTs = numpy.argmax(PFTArrays, axis=0) + 1

    # Now apply the mask, based on the sum of the PFT fractions
    PFTs = numpy.ma.masked_array(PFTs, mask=Mask)

    # Send it to the ncvar creator
    define_ncvar(
            PFTs,
            NCFile,
            "iveg",
            dims["iveg"],
            attribs["iveg"]
            )

    # Now use this information to get the attributes which are dependent on
    # the PFT- this is LAI and isoil (with rhosoil and css being dependent on
    # isoil).
    
    # Soil Moisture (all layers)
    # Should retrieve stash entries for layer 1, 2, 3, 4, 5, 6
    SoilMoistureStash = FieldsFile.stashmaster.by_regex("CABLE SOIL MOISTURE")

    # We know the desired size of the system from the dims- initialise it with
    # values of -1.0 for easier masking
    SoilMoisture = numpy.ones(
            tuple(len(NCFile.dimensions[dim]) for dim in dims["SoilMoist"]),
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
            NCFile,
            "SoilMoist",
            dims["SoilMoist"],
            attribs["SoilMoist"]
            )

    # Now do soil temperature, in the same manner as soil moisture
    # Should retrieve stash entries for layer 1, 2, 3, 4, 5, 6
    SoilTempStash = FieldsFile.stashmaster.by_regex("CABLE SOIL TEMPERATURE")

    # We know the desired size of the system from the dims- initialise it with
    # values of -1.0 for easier masking
    SoilTemp = numpy.ones(
            tuple(len(NCFile.dimensions[dim]) for dim in dims["SoilTemp"]),
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
            NCFile,
            "SoilTemp",
            dims["SoilTemp"],
            attribs["SoilTemp"],
            )

    # Now do snow depth- ACCESS-ESM1.5 has 3 snow layers, but CABLE only has 1.
    # Retrieve the snow depth stash items
    # Should get layer 1, 2, 3
    SnowDepthStash = FieldsFile.stashmaster.by_regex("SNOW DEPTH LAYER")

    # Initialise the array- don't initialise to -1.0, as we want to summate
    # over layers
    SnowDepth = numpy.zeros(
            tuple(len(NCFile.dimensions[dim]) for dim in dims["SnowDepth"]),
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

    # Apply the mask- since we didn't initialise the array to -1.0, we can't
    # just use numpy.ma.masked_less. Inherit the mask from the PFT array, and
    # broadcast it to the desired shape by repeating over the time axis.
    SnowDepth = numpy.ma.masked_array(
            SnowDepth,
            mask=numpy.repeat(Mask[numpy.newaxis, :, :], 12, axis=0)
            )

    define_ncvar(
            SnowDepth,
            NCFile,
            "SnowDepth",
            dims["SnowDepth"],
            attribs["SnowDepth"]
            )

    # LAI
    # Retrieve the LAI stash item
    LAIStash = FieldsFile.stashmaster.by_regex("CASA LEAF AREA INDEX")

    # Should be a single entry- retrieve the stash code
    StashCode = list(LAIStash.values())[0].item

    # Initialise the array of LAI values
    LAI = numpy.zeros(
            tuple(len(NCFile.dimensions[dim]) for dim in dims["LAI"]),
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

    define_ncvar(LAI, NCFile, "LAI", dims["LAI"], attribs["LAI"])

    # Now do isoil- 2 everywhere except where PFT=19, where isoil=9
    iSoil = numpy.ma.where(PFTs == 17, 9, 2)

    define_ncvar(iSoil, NCFile, "isoil", dims["isoil"], attribs["isoil"])

    # Variables that depend on isoil
    # rhosoil
    rhoSoil = numpy.ma.where(iSoil == 2, 1600.0, 900.0)

    define_ncvar(
            rhoSoil,
            NCFile,
            "rhosoil",
            dims["rhosoil"],
            attribs["rhosoil"]
            )

    # css
    css = numpy.ma.where(iSoil == 2, 850.0, 2100.0)

    define_ncvar(css, NCFile, "css", dims["css"], attribs["css"])

#----------------------------------------------------------------------------#    

def convert_soilorder_to_int(FieldsFile, NCFile, dims, attribs):
    """Convert the soil order from floating point to integer."""

    # Get the stash entry associated
    SoilOrderStash = FieldsFile.stashmaster.by_regex("CASA SOIL ORDER")
    StashCode = list(SoilOrderStash.values())[0].item

    # Get the field associated with the stash code
    for Field in FieldsFile.fields:
        if Field.lbuser4 == StashCode:
            break

    # Convert it to int
    SoilOrder = Field.get_data().astype(numpy.int32)

    # Mask the array
    SoilOrder = numpy.ma.masked_less(SoilOrder, 1)

    # Pass it to the variable creator
    define_ncvar(
            SoilOrder,
            NCFile,
            "SoilOrder",
            dims["SoilOrder"],
            attribs["SoilOrder"]
            )

#----------------------------------------------------------------------------#    

def direct_conversions(FieldsFile, NCFile, dims, attribs):
    """Convert all variables that are a 1:1 conversion, with no level
    dimension, from a UM field to a CABLE field."""

    # Build a dictionary of conversion mappings
    mappings = {
            "GRIDBOX AREAS": "area",
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
        Data = numpy.ma.masked_outside(Field.get_data(), 0.0, 1e6)

        define_ncvar(
                Data,
                NCFile,
                CABLEVar,
                dims[CABLEVar],
                attribs[CABLEVar]
                )

#----------------------------------------------------------------------------#    

def scaling_conversions(FieldsFile, NCFile, dims, attribs, scalings):
    """Perform all conversions that require scalar/unit conversions to go from
    UM fields to CABLE fields."""

    # Current variables that fit this description are Ndep and hyds
    mappings = {
            "CASA NITROGEN DEPOSITION": "Ndep",
            "SAT SOIL CONDUCTIVITY AFTER TIMESTEP": "hyds"
            }

    for UMVar, CABLEVar in mappings.items():
        # Get the stash code
        VarStash = FieldsFile.stashmaster.by_regex(UMVar)
        StashCode = list(VarStash.values())[0].item

        # Get the field from the stash code
        for Field in FieldsFile.fields:
            if Field.lbuser4 == StashCode:
                break

        define_ncvar(
                Field.get_data() * scalings[CABLEVar],
                NCFile,
                CABLEVar,
                dims[CABLEVar],
                attribs[CABLEVar]
                )

#----------------------------------------------------------------------------#    

if __name__ == "__main__":
    # Start by opening the desired UM fields file and attaching the stash
    FieldsFile = mule.FieldsFile.from_file("PI-02.astart-01010101")
    # Only take section 0, other sections are diagnostics
    UMStash = mule.STASHmaster.from_file("STASHmaster_A").by_section(0)
    FieldsFile.attach_stashmaster_info(UMStash)

    # Open the NCfile to write to
    NCFile = netCDF4.Dataset(
            "ACCESS-ESM1p5-gridinfo-1p875x1p25-CABLE.nc",
            "w",
            format="NETCDF4"
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

    for Dimension, DimLength in NCDimensions.items():
        NCFile.createDimension(Dimension, DimLength)

    # Add top level attributes
    NCAttribs = {
            "source": "ACCESS-ESM1.5 atmosphere restart file from " +\
                    "CMIP6 pre-industrial runs.",
            "comments": "Converted from the UM fields file on Gadi at " +\
                    "/g/data/vk83/configurations/inputs/access-esm1p5/modern/pre-industrial/restart/atmosphere/PI-02.astart-01010101 " +\
                    "with the python script at [SCRIPT LOCATION HERE].",
            "contact": "lachlan.whyborn@anu.edu.au"
            }

    for AttribName, AttribContent in NCAttribs.items():
        NCFile.setncattr(AttribName, AttribContent)

    # Start by adding dimension variables
    NCFile.createVariable(
            "longitude",
            'f4',
            ('longitude',)
            )

    NCFile.createVariable(
            "latitude",
            'f4',
            ('latitude',)
            )

    # Remember to map longitude from 0.0 -> 360.0 to -180.0 -> 180.0
    Longitudes = numpy.linspace(0.0, 360.0, 192, endpoint=False)
    Longitudes = numpy.where(Longitudes > 180.0, -360 + Longitudes, Longitudes)
    NCFile["longitude"][:] = Longitudes

    Latitudes = numpy.linspace(-90.0, 90.0, 145)
    NCFile["latitude"][:] = Latitudes

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
            "albedo2"   : ("latitude", "longitude")
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
                }
            }

    scalings = {
            "Ndep"  : 1.0 / 365.0,
            "hyds"  : 1.0 / 1000.0
            }

    # Call the variable builders
    compute_iveg_and_dependent_vars(
            FieldsFile,
            NCFile,
            VarDimensions,
            VarAttribs
            )

    write_albedo(
            FieldsFile,
            NCFile,
            VarDimensions,
            VarAttribs
            )

    convert_soilorder_to_int(
            FieldsFile,
            NCFile,
            VarDimensions,
            VarAttribs
            )

    direct_conversions(
            FieldsFile,
            NCFile,
            VarDimensions,
            VarAttribs
            )

    scaling_conversions(
            FieldsFile,
            NCFile,
            VarDimensions,
            VarAttribs,
            scalings
            )
             

            
