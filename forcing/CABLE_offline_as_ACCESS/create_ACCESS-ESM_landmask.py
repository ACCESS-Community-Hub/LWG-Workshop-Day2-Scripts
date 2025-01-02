# Create a reduced landmask based on provided command line arguments

def read_CLI_args():
    """Read the command line arguments."""

    # Set up the parser with default values
    parser = argparse.ArgumentParser(
            prog="create_landmask",
            description="Create a landmask file from a given reference mask."
            )

    parser.add_argument(
            "-r",
            "--ref_landmask",
            help="Landmask file to use as reference",
            default="/g/data/rp23/experiments/2024-03-12_CABLE4-dev/lw5085/CABLE-as-ACCESS/landmasks/ACCESS-ESM1p5-1p875x1p25-landmask.nc"
            )

    parser.add_argument(
            "-d",
            "--domain",
            help="Domain of the new landmask, given as lon_min,lon_max," +\
                    "lat_min,lat_max",
            default="-180.0,180.0,-90.0,90.0"
            )

    parser.add_argument(
            "-p",
            "--points",
            help="List of points to include in the landmask, given as " +\
                    "lon1,lat1:lon2,lat2:lon3,lat3.... Takes the nearest " +\
                    "land point to each requested point."
            default=""
            )

    args = parser.parse_args()

    return args

#-----------------------------------------------------------------------------#

def create_landmask(RefLandmask, Domain, Points, Outfile):
    """Use the reference landmask and the given domain to create a landmask
    usable with the GSWP format."""

    # Load up the reference landmask
    Landmask = xarray.open_dataset(RefLandmask)

    # Before we restrict the mask to the domain, determine the indices of the
    # land points closest to the given points

    # Destructure the domain
    LonMin, LonMax, LatMin, LatMax = Domain

    # Set everything outside the selected domain to 0
    Landmask["mask"].loc[{"
