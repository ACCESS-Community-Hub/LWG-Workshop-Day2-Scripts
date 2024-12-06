## Purpose

Converts directly from the ACCESS-ESM1.5 atmosphere restart file in UM Fields format to a NetCDF file usable by CABLE offline.

## Usage

Invoked at the command line via

```
python3 convert_from_fields.py --input(-i) <ACCESS_restart_file> --stash(-s) <UM_stash_file> --output(-o) <file_to_write.nc> --areafile(-a) <file_with_gridcell_area>
```

There are defaults for each argument which point to existing files on Gadi:
--input:    /g/data/vk83/configurations/inputs/access-esm1p5/modern/pre-industrial/restart/atmosphere/PI-02.astart-01010101
--stash:    /g/data/rp23/experiments/2024-03-12\_CABLE4-dev/lw5085/CABLE-as-ACCESS/STASHmaster\_A
--areafile: /g/data/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/fx/areacella/gn/v20191115/areacella\_fx\_ACCESS-ESM1-5\_historical\_r1i1p1f1\_gn.nc
--output:   ACCESS-ESM1p5-1p875x1p25-gridinfo-CABLE.nc

The script the ```xarray``` and ```mule``` python packages, which are available through the hh5 ```conda_concept``` module.

## Method

The specifics of the conversion process for each of the variables is described in the table below.

| GridInfo Name | UM Stash Code | Comments |
|:-------------:|:--------------|----------|
| longitude | N/A | ACCESS-ESM1.5 operates on a 0.0 to 360.0 domain, but CABLE expects -180.0 to 180.0. |
| latitude | N/A | |
| area | N/A | Taken from the "areacella" variable in the areafile. The ACCESS restart contains a "GRIDBOX AREA" field, but it is filled with 0.0. |
| iveg | FRACTIONS OF SURFACE TYPES | ACCESS restart provides respective fractions of each vegetation type in a lon x lat x nvegtypes array. Take the dominant vegetation type for the CABLE gridinfo. The generated iveg array is used to generate the landmask for the other derived variables. Given units of 1 and categorical flag descriptions. |
| patchfrac | Land fraction in grid box | Currently set to 1.0 for all land points. |
| isoil | N/A | Set to 2 for all vegetation types bar 17, which is set to 9. |
| SoilMoist | CABLE SOIL MOISTURE | There are 6 lon x lat x PFT arrays that match this stash code, one for each layer. The soil moisture is taken from the dominant PFT in that grid cell, and converted from kg/m^2 to m^3/m^3 by dividing by layer thickness and the density of water. CABLE takes monthly values, so copy the value across each month. |
| SoilTemp | CABLE SOIL TEMPERATURE | As with soil moisture, there are 6 lon x lat x PFT arrays that match this stash code, one for each layer. Soil temperature is taken from the dominant PFT in that grid cell. Copied across each month. |
| SnowDepth | SNOW DEPTH LAYER | There are 3 lon x lat x PFT arrays that match this stash, one for each layer of snow in the ACCESS model. CABLE only takes one soil layer, so summate the snow depth in each layer from the dominant PFT. Copied across each month. |
| Albedo | N/A | Set to 0.2. Not used but bounds checked (?). |
| LAI | CASA LEAF AREA INDEX| Take the LAI from the dominant PFT. CABLE needs monthly values, so copy values across length 12 time axis. Given units of 1. |
| SoilOrder | CASA SOIL ORDER | Convert from floating point values to integer. Given units of 1. |
| Ndep | NITROGEN DEPOSITION | Convert from g/m<sup>2</sup>/day to g/m<sup>2</sup>/year by multiplying by 365, and supplied with units metadata. Given long name "annual Nitrogen deposition rate". |
| Nfix | NITROGEN FIXATION | Appears to already by in g/m<sup>2</sup>/year. Given appropriate units and long name "annual Nitrogen fixation rate". |
| Pwea | PHOSPHORUS FROM WEATHERING | Given units of g/m<sup>2</sup>/year and long name "annual yield of Phosphorous from weathering". |
| Pdust | PHOSPHORUS FROM DUST | Given units of g/m<sup>2</sup>/year and long name "annual yield of Phosphorous from dust". |
| clay | Dust parent soil clay fraction | Given units of 1 and long name "soil clay fraction". |
| silt | Dust parent soil silt fraction | Given units of 1 and long name "soil silt fraction". |
| sand | Dust parent soil sand fraction | Given units of 1 and long name "soil sand fraction". |
| swilt | VOL SMC AT WILTING AFTER TIMESTEP | Given units of 1 and long name "Volume soil moisture content at wilting". |
| sfc | VOL SMC AT CRIT PT AFTER TIMESTEP | Given units of 1 and long name "Volume soil moisture content at crit point". |
| ssat | VOL SMC AT SATURATION AFTER TIMESTEP | Given units of 1 and long name "Volume soil moisture content at saturation". |
| bch | CLAPP-HORNBERGER "B" COEFFICIENT | Given units of 1 and long name "CLAPP-HORNBERGER coefficient". |
| hyds | SAT SOIL CONDUCTIVITY AFTER TIMESTEP | Divide by 1000 to go from kg/m<sup>2</sup>/s to m/s, and given appropriate units. Given long name "Saturated soil conductivity". |
| sucs | SATURATED SOIL WATER SUCTION | Given units of m and long name "Saturated soil water suction". |
| rhosoil | N/A | Set according to isoil. 1600.0 where isoil=2, 910.0 where isoil=9. Given units of kg/m<sup>3</sup> and long name "Soil density". |
| cnsd | THERMAL CONDUCTIVITY AFTER TIMESTEP | Given units of W/m/K and long name "Soil thermal conductivity". |
| css | N/A | N/A | N/A | Set according to isoil. 850.0 where isoil=2, 2100.0 where isoil=9. Given units of J/kg/K and long name "Soil specific heat capacity". |
| albedo2 | SNOW-FREE ALBEDO OF SOIL | Given units of 1 and long name "Snow-free albedo of soil". |
