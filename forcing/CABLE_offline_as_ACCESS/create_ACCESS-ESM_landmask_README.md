## Purpose

Create a new landmask usable by CABLE from a reference landmask either by restricting the domain or by selecting a series of points.

## Usage

Call from the command line via

```create_ACCESS-ESM_landmask.py -r <reference_landmask> -d=<domain> -p=<points> -o <outfile>```

where:
* ```-r <reference_landmask>``` is the landmask to use as a reference. Must be at the same resolution as the meteorology.
* ```-d=<domain>``` is the domain to capture by the landmask, given as "LonMin,LonMax,LatMin,LatMax". Use the ```d=...``` format rather than ```d ...``` format, as the interpreter cannot distingush between negative signs and dashes (which indicate a new argument) in the latter.
* ```-p=<points>``` is a list of points to include in the mask, given as "Lon1,Lat1:Lon2,Lat2:Lon3,Lat3...". If points are specified, the domain argument is ignored. Any points more than 5 degrees from any land point are ignored. For the same reasons as ```-d=<domain>```, use the ```p=<points>``` format.
* ```-o=<outfile>``` is the file to write the new landmask to.

The script uses the ```xarray``` and ```scipy.spatial``` packages, which are both available through the ```hh5``` ```conda_concept``` environment on Gadi.

## Method

The reference landmask is read by the ```xarray``` package. When the landmask is processed by the ```domain``` argument, all points outside the given domain are set to 0. When using the ```points``` argument, the land points on the globe are reduced to a one-dimensional array, with the nearest points being located by querying a ```scipy.spatial.KDTree``` constructed from the stacked latitude and longitude coordinates.

