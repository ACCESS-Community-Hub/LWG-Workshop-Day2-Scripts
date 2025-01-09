# Create a reduced landmask based on provided command line arguments
import argparse
import numpy
import scipy
import xarray

def read_CLI_args():
    """Read the command line arguments."""

    # Set up the parser with default values
    parser = argparse.ArgumentParser(
            prog='create_landmask',
            description='Create a landmask file from a given reference mask.'
            )

    parser.add_argument(
            '-r',
            '--ref_landmask',
            help='Landmask file to use as reference',
            default='/g/data/rp23/experiments/2024-03-12_CABLE4-dev/lw5085/CABLE-as-ACCESS/landmasks/ACCESS-ESM1p5-1p875x1p25-landmask.nc'
            )

    parser.add_argument(
            '-d',
            '--domain',
            help='Domain of the new landmask, given as lon_min,lon_max,' +\
                    'lat_min,lat_max.' +\
                    'NOTE: use -d=... rather than -d ..., as minus signs ' +\
                    'are indistinguishable from dashes indicating a new arg',
            default='-180.0,180.0,-90.0,90.0'
            )

    parser.add_argument(
            '-p',
            '--points',
            help='List of points to include in the landmask, given as ' +\
                    'lon1,lat1:lon2,lat2:lon3,lat3.... Takes the nearest ' +\
                    'land point to each requested point, within a radius ' +\
                    'of 5 degrees. Overrides the --domain argument. ' +\
                    'NOTE: use -p=... rather than -p ..., as minus signs ' +\
                    'are indistinguishable from dashes indicating a new arg',
            default=''
            )

    parser.add_argument(
            '-o',
            '--outfile',
            help='Name to write the landmask out to.',
            default='ACCESS-ESM1p5-1p875x1p25-landmask-processed.nc'
            )

    args = parser.parse_args()

    return args

#-----------------------------------------------------------------------------#

def restrict_to_domain(Landmask, Domain):
    """Restrict the land points to the given domain."""
    
    # Destructure the domain
    LonMin, LonMax, LatMin, LatMax = Domain

    # Set everything outside the domain to int8(0)
    DomainMask = (
            (Landmask.coords['latitude'] < LatMin) |
            (Landmask.coords['latitude'] > LatMax) |
            (Landmask.coords['longitude'] < LonMin) |
            (Landmask.coords['longitude'] > LonMax)
            )

    Landmask['mask'] = xarray.where(DomainMask, 0, Landmask['mask'])
    
    return Landmask

#-----------------------------------------------------------------------------#

def select_points_from_landmask(Landmask, Points):
    """Create a new landmask with only the land points closest to the selected
    points included."""

    # First, we want to remove all non-land points from the original mask
    OnlyLandmask = Landmask.where(Landmask['mask'] == 1, drop=True)

    # Flatten the DataSet to acquire a single MultiIndex
    OnlyLandmask = OnlyLandmask.stack(coordinates=['latitude', 'longitude'])

    # Calling to.numpy() on the MultiIndex unfortunately returns an array of
    # tuples, rather than an array of dimension (N,2). So create the (N,2)
    # array and fill it with the coordinates
    CoordinateArray = numpy.empty(
            (len(OnlyLandmask.indexes['coordinates']), 2),
            dtype=numpy.float32
            )

    for i, Coord in enumerate(OnlyLandmask.indexes['coordinates']):
        CoordinateArray[i, 0], CoordinateArray[i, 1] = Coord

    # Now we can pass this array to the KDTree builder
    Tree = scipy.spatial.KDTree(CoordinateArray)

    # We can pass the set of Points directly to the tree
    Distances, Indexes = Tree.query(Points)

    # Now we can set the original mask to all zeroes, and then set the selected
    # points to 1
    Landmask['mask'][:] = 0

    # NearestPoints is tuple of arrays of (Distances, Indices) tuples
    for PtIndex, (Distance, Index) in enumerate(zip(Distances, Indexes)):
        # First, check whether the distance is less than a given tolerance
        if Distance > 5.0:
            print(f'Warning: the point {Points[PtIndex]} is more than 5 ' +\
                    'degrees from the nearest land point, so will be ignored.')
        else:
            # We need to go from an index in the MultiIndex to a lon/lat
            Lat, Lon = CoordinateArray[Index]

            # Pass these to the loc accessor and set the value to 1
            Landmask['mask'].loc[{'latitude': Lat, 'longitude': Lon}] = 1
    
    return Landmask

#-----------------------------------------------------------------------------#

def parse_domain(Domain):
    """Process the passed domain argument into a list of numbers denoting the
    domain corners."""
    Domain = [float(corner) for corner in Domain.split(',')]

    # Check validity of processed domain
    assert len(Domain) == 4, f'{Domain} is invalid input for --domain.'

    return Domain

#-----------------------------------------------------------------------------#
def parse_points(Points):
    """Process the passed points arguments into a list of coordinates."""
    Points = [tuple(
        float(Component) for Component in reversed(Coord.split(','))
        )
        for Coord in Points.split(':')]

    # Check validity of processed points
    for Point in Points:
        assert len(Point) == 2, f'{Points} is invalid input for --points.'

    return Points

#-----------------------------------------------------------------------------#

def create_landmask(RefLandmask, Domain, Points, Outfile):
    """Use the reference landmask and the given domain to create a landmask
    usable with the GSWP format."""

    # Load up the reference landmask
    Landmask = xarray.open_dataset(RefLandmask)

    if Points is not None:
        # The user has specified a list of points 
        ProcessedLandmask = select_points_from_landmask(Landmask, Points)
    else:
        # Just specified domain
        ProcessedLandmask = restrict_to_domain(Landmask, Domain)

    # Should be able to just write out the new landmask?
    ProcessedLandmask.to_netcdf(Outfile)

#-----------------------------------------------------------------------------#

if __name__ == '__main__':
    args = read_CLI_args()

    if args.points == '':
        Points = None
        Domain = parse_domain(args.domain)
    else:
        Points = parse_points(args.points)
        Domain = None
    
    create_landmask(args.ref_landmask, Domain, Points, args.outfile)
