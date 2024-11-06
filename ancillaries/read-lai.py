#!/usr/bin/env python
__title__ = "Test pre-processing of sample parameter of JULES ancillary data"
__version__ = "2024-10-10"
__author__ = "Abhaas Goyal"
__email__ = "abhaas.goyal@anu.edu.au"
__institution__ = "ACCESS-NRI"


"""
This is a simple use case to learn how to use ANTS.
Reads in LAI variable and applies AGCD Landmask (australia-wide).

NOTE: LAI file needs to be CF compliant w.r.t. units for `ants.analysis.make_consistent_with_lsm` to work.
For example, used `ncatted` to add units for longitude/latitude
```
ncatted -h -a units,latitude,o,c,degrees_north -a units,longitude,o,c,degrees_east lai_bom_2000_2022_climatology.nc lai_bom_2000_2022_climatology.nc
```

To use `ants` python module on gadi use: 

module purge
module use /g/data/access/ngm/modules
module load ants/0.18
"""

import ants


def main():

    lai_path = "/g/data/rp23/experiments/2024-10-10_LWG_workingbee/ag9761/lai_bom_2000_2022_climatology_v2.nc"
    lsm_path = "/g/data/rp23/experiments/2024-10-10_LWG_workingbee/JULES_ancil_bom/Landmask_Australia_agcd_5km.nc"
    output_path = "/g/data/rp23/experiments/2024-10-10_LWG_workingbee/ag9761/lai.nc"

    # Load LAI map
    source_cube = ants.load_cube(lai_path)

    print("Loaded Source Cube")

    # Load Land-Sea Mask in boolean form
    target_lsm = ants.load_cube(lsm_path, "land_mask")
    target_lsm = target_lsm.copy(target_lsm.data.astype("bool", copy=False))

    print("Loaded Land-Sea Mask")

    # Make LAI consistent with Land Sea Mask
    # REVIEW: Getting stuck with invert_mask=False, but does work as expected
    # `ant` documentation hints at denoting land values as True by default
    # So not understanding why the opposite is required
    # https://code.metoffice.gov.uk/doc/ancil/ants/latest/_modules/ants/analysis.html#make_consistent_with_lsm
    ants.analysis.make_consistent_with_lsm(source_cube, target_lsm, invert_mask=True)

    # Save
    print(f"Saving to {output_path}")
    ants.fileformats.netcdf.cf.save(
        source_cube,
        output_path,
        unlimited_dimensions=["time"],
        fill_value=float("nan"),
        netcdf_format="NETCDF4",
    )
    print("Saved file")


if __name__ == "__main__":
    main()
