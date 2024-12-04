author: mathew lipson m.lipson@unsw.edu.au

# Purpose
This script emulates the rose/cycl suites that generate ancillaries 
(e.g. vegetation fraction and land-sea mask) required to run JULES offline or 
with the UM. It is based on RAS suite (u-bu503) but uses ants/1.1.0 per 
Siyuan's offline code.

# Usage
from the command line: `qsub run_ancil_lct.sh`
this uses 1 cpu and takes ~ 2 mins to process this example domain

The domain is defined in `./ants_lct_inputs/grid.nl` using the same format as in the RAS.
By default the lsm mask is defined from the land cover source, however an alternative lsm
can be provided as an ancil_lct.py argument (not yet tested).

Other ANTS configuration inputs are defined in `./ants_lct_<app-name>-app.conf`, for example
the ants x/y decomposition and regridding schemes.

# Method

Uses a variety of methods and ants functionality, including
`ants.analysis.make_consistent_with_lsm`
`ants.regrid.GeneralRegridScheme` with a linear horizontal scheme

Outputs are placed by default into `./outputs_1p1` (configurable).
Following the ancillary generation suites, outputs are saved in both UM (no extension) and .nc formats.
Some files are not required at all by current versions of JULES (e.g. `qrparm.landfrac`).
Some files are only required as an intermediatory step (e.g. `qrparm.veg.frac_cci_pre_c4`, `c4_percent_1d.nc`)

The final files important for JULES are `qrparm.mask` and `qrparm.veg.frac_cci`

# Other notes

## Siyuans offline JULES ancillary suite
svn co https://code.metoffice.gov.uk/svn/roses-u/b/x/0/3/8/nci-gadi

## The land-atmos Regional Ancillary Suite (RAS)
svn co https://code.metoffice.gov.uk/svn/roses-u/b/u/5/0/3/trunk

# The RAS uses ANTS 1.0 branch: 
https://code.metoffice.gov.uk/svn/ancil/ants/branches/dev/claudiosanchez/r13230_RAS_ANTS_1p0

# Siyuans ancil suite uses ANTS 1.1 branch:
https://code.metoffice.gov.uk/svn/ancil/ants/branches/dev/martinbest/r13910_ancil_ants_JULES_vn1.1.0/
