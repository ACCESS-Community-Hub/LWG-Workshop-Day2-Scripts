#!/bin/bash
#PBS -q normal
#PBS -l ncpus=4
#PBS -l mem=16GB
#PBS -l walltime=01:00:00
#PBS -l storage=gdata/hh5+gdata/access+gdata/rp23
#PBS -l wd
#PBS -l jobfs=1GB
#PBS -P rp23

# This script emulates the rose/cycl workflows that generate ancillaries 
# required to run JULES offline or with the UM-JULES Regional Ancillary Suite (RAS). 
# It is based on RAS suite (u-bu503) but uses ants/1.1.0 per 
# Siyuan's offline code. The main difference between the RAS (1.0) and Siyuan's (1.1)
# is the use of --use-new-saver, and the merging of soil roughness ancil with general
# soil parameters ancil

# this example with grid.nl of 490x386 takes approxiamately 35 mins to process on 4 cpus.

set -e
module purge
module use /g/data/access/ngm/modules
module load ants/1.1.0

ANTS_SRC_PATH=/g/data/rp23/experiments/2024-10-10_LWG_workingbee/mjl561/ants_src/ants_1p1
ANCIL_MASTER=/g/data/access/TIDS/UM/ancil/atmos/master
INPUT_PATH=/home/561/mjl561/git/LWG-Workshop-Day2-Scripts/ancillaries/JULES/ants_lct_inputs
ANCIL_TARGET_PATH=/g/data/rp23/experiments/2024-10-10_LWG_workingbee/mjl561/outputs_1p1
ANCIL_PREPROC_PATH=/g/data/access/TIDS/RMED/ANTS/preproc

# check if the output directory exists
if [ ! -d ${ANCIL_TARGET_PATH} ]; then
    mkdir -p ${ANCIL_TARGET_PATH}
fi

# ============================================================================
# config based on RAS (u-bu503) app/ancil_lct/rose-app.conf, walltime ~ 1 min
source=${ANCIL_MASTER}/vegetation/cover/cci/v3/vegetation_fraction.nc
target_grid=${INPUT_PATH}/grid.nl
transformpath=/g/data/access/TIDS/UM/ancil/data/transforms/cci2jules_ra1.json
output_vegfrac=${ANCIL_TARGET_PATH}/qrparm.veg.frac_cci_pre_c4
output_lsm=${ANCIL_TARGET_PATH}
ANTS_CONFIG=${INPUT_PATH}/ancil_lct-app.conf

# ancil_lct.py takes the global categorical land cover product CCI [source]
# and converts it to JULES tiles using a crosswalking table [transformpath]
# it outputs landcover ancillary file [output_vegfrac], but as CCI does not 
# distinguish between c3/c4 grasses, another step is required (see below).
# see: https://code.metoffice.gov.uk/doc/ancil/ants/1.1/bin/ancil_lct.html
# outputs:
#   qrparm.veg.frac_cci_pre_c4.nc (intermediate file)
#   qrparm.mask_sea.nc (unknown use)
#   qrparm.mask.nc (land-sea mask)
python ${ANTS_SRC_PATH}/ancil_lct.py ${source} \
         --target-grid ${target_grid} --transform-path ${transformpath} \
         -o ${output_vegfrac} --landseamask-output ${output_lsm}       \
         --use-new-saver --ants-config ${ANTS_CONFIG}

# ============================================================================
# config based on RAS app/ancil_lct_postproc_c4/rose-app.conf, walltime ~ 1 min
source=${ANCIL_TARGET_PATH}/qrparm.veg.frac_cci_pre_c4.nc
target_lsm=${ANCIL_TARGET_PATH}/qrparm.mask
output=${ANCIL_TARGET_PATH}/qrparm.veg.frac_cci
c4source=${ANCIL_MASTER}/vegetation/cover/cci/v3/c4_percent_1d.nc
ANTS_CONFIG=${INPUT_PATH}/ancil_lct_postproc_c4-app.conf

# ancil_general_regrid.py regrids the c4 fraction [c4source] to the target grid [target_lsm]
# The config [ANTS_CONFIG] defines a linear horizontal regriddin scheme because the source
# data is very low resolution (1 degree) where normal area-weighted regridding results in
# sharp changes to c3/c4 fractions.
# see: https://code.metoffice.gov.uk/doc/ancil/ants/1.1/bin/ancil_general_regrid.html

# Then ancil_lct_postproc_c4.py uses the intermediate regridded c4 fraction source 
# [c4_percent_1d.nc] and includes it in a new 9 tile landcover fraction ancillary [output]
# see: https://code.metoffice.gov.uk/doc/ancil/ants/1.1/bin/ancil_lct_postproc_c4.html
# original 1 degree source for C3/C4 partitioning: ISLSCP II C4 Vegetation Percentage
# http://doi.org/10.3334/ORNLDAAC/932 (this being converted to nc on gadi)
# outputs: 
#   c4_percent_1d.nc (intermediate file)
#   qrparm.veg.frac_cci.nc (land fraction including c3 & c4 partition)

python ${ANTS_SRC_PATH}/ancil_general_regrid.py --ants-config ${ANTS_CONFIG} \
       ${c4source} --target-lsm ${target_lsm} -o ${ANCIL_TARGET_PATH}/c4_percent_1d.nc

python ${ANTS_SRC_PATH}/ancil_lct_postproc_c4.py \
       ${source} --islscpiic4 ${ANCIL_TARGET_PATH}/c4_percent_1d.nc \
       --use-new-saver -o ${output}

# ============================================================================
# config based on RAS (u-bu503) app/ancil_lai/rose-app.conf, walltime ~ 1 min
TRANSFORM_DIR=/g/data/access/TIDS/UM/ancil/data/transforms
output=${ANCIL_TARGET_PATH}/lai
source="${ANCIL_MASTER}/vegetation/lai/modis_4km/v2/lai_preprocessed.nc \
        ${ANCIL_TARGET_PATH}/qrparm.veg.frac_cci.nc"
relative_weights=${TRANSFORM_DIR}/lai_weights_gl9.json
ANTS_CONFIG=${INPUT_PATH}/ancil_lai-app.conf

# ancil_lai.py is a CONTRIB app that derives leaf area index for each vegetation 
# tile based on # the total MODIS LAI split into weighted fractions depending on 
# [relative_weights] which by default is: [5.0, 4.0, 2.0, 4.0, 1.0] on each tile.
# more documentation is available in the file (no online documentation)
# outputs:
#   lai.nc (intermediate file)
python ${ANTS_SRC_PATH}/Apps/Lai/ancil_lai.py ${source} \
       --relative-weights ${relative_weights} -o ${output} \
       --use-new-saver --ants-config ${ANTS_CONFIG}

# ============================================================================
# config based on RAS (u-bu503) app/ancil_canopy_heights/rose-app.conf, walltime ~ 1 min
canopy_height_factors=${TRANSFORM_DIR}/canopy_height_factors_gl9.json
lai=${ANCIL_TARGET_PATH}/lai.nc
output=${ANCIL_TARGET_PATH}/qrparm.veg.func
trees=${ANCIL_MASTER}/vegetation/canopy/simard/v1/Simard_Pinto_3DGlobalVeg_JGR.nc
ANTS_CONFIG=${INPUT_PATH}/ancil_lai-app.conf

# ancil_canopy_heights.py is a CONTRIB app that uses a global source [trees] to 
# represent the maximum canopy height across the year in each month for each PFT
# more documentation is available in the file (no online documentation)
# Standard settings for this app fail for some domains because of problems on
# isolated islands (e.g. Macquarie Island). A fix is to increase the neighbourhood
# search readias from 500 km to, say 500km:
# loop_lim_y = index_nearest_neighbour.ydist2index(trees, 2500)
# see discussion here: 
# https://forum.access-hive.org.au/t/aus2200-vegetation-fraction-ancil-creation-issues/1972/19
# outputs:
#   canopy_heights.nc (intermediate file)
python ${ANTS_SRC_PATH}/Apps/CanopyHeights/ancil_canopy_heights.py ${lai} \
       --canopy-height-factors ${canopy_height_factors} --trees-dataset ${trees} \
       -o ${ANCIL_TARGET_PATH}/canopy_heights.nc --use-new-saver \
       --ants-config ${ANTS_CONFIG}

# ancil_2anc.py combines [lai.nc] and [canopy_heights.nc] into one ancillary
# see: https://code.metoffice.gov.uk/doc/ancil/ants/1.1/bin/ancil_2anc.html
# outputs:
#   qrparm.veg.func.nc (veg canopy properties)
python ${ANTS_SRC_PATH}/ancil_2anc.py \
       ${lai} ${ANCIL_TARGET_PATH}/canopy_heights.nc --use-new-saver -o ${output}

# ============================================================================
# config based on RAS (u-bu503) app/ancil_soils_hydr/rose-app.conf, walltime ~ 19 mins
source=${ANCIL_PREPROC_PATH}/soils_hydrology.nc
vegfrac=${ANCIL_TARGET_PATH}/qrparm.veg.frac_cci
soils_lookup=${ANCIL_PREPROC_PATH}/soils_hydrology_lookup_clapp_hornberger.json
output=${ANCIL_TARGET_PATH}/soil_hydrology
ANTS_CONFIG=${INPUT_PATH}/ancil_soils_hydr-app.conf

# ancil_soils.py is a CONTRIB app that creates a large number of soils hydrology, carbon
# and thermal property parameter ancillaries using a dominant approach. This task takes 
# the longest in this workflow. It is only slightly faster with more cpus.
# more documentation is available in the file (no online documentation)
# outputs: 
#   soil_hydrology.nc (soil properties)
python ${ANTS_SRC_PATH}/Apps/SoilParameters/ancil_soils.py ${source} \
       --lct-ancillary ${vegfrac} --soils-lookup ${soils_lookup} \
       -o ${output} --use-new-saver --ants-config ${ANTS_CONFIG}

# ============================================================================
# config based on RAS (u-bu503) app/ancil_topographic_index/rose-app.conf, walltime ~ 1 min
vegfrac=${ANCIL_TARGET_PATH}/qrparm.veg.frac_cci.nc
output=${ANCIL_TARGET_PATH}/qrparm.hydtop
source=${ANCIL_PREPROC_PATH}/topographic_index.nc
ANTS_CONFIG=${INPUT_PATH}/ancil_topographic_index-app.conf

# ancil_topgoraphic_index.py derives paramaters for TOP model hydrology
# see: https://code.metoffice.gov.uk/doc/ancil/ants/1.1/bin/ancil_topographic_index.html
# outputs: 
#      qrparm.hydtop.nc (for TOP model hydrology)
python ${ANTS_SRC_PATH}/ancil_topographic_index.py \
       ${source} --lct-ancillary ${vegfrac}  --use-new-saver -o ${output} \
       --ants-config ${ANTS_CONFIG}

# ============================================================================
# config based on RAS (u-bu503) app/ancil_snowfree_albedo/rose-app.conf, walltime ~ 1 min
source=${ANCIL_MASTER}/land_clim/GlobAlbedo/v2/qrclim.land.nc
target_lsm=${ANCIL_TARGET_PATH}/qrparm.mask
output=${ANCIL_TARGET_PATH}/qrclim.land
ANTS_CONFIG=${INPUT_PATH}/ancil_snowfree_albedo-app.conf

# ancil_general_regrid.py regrids global albedo [source] to the target grid [target_lsm]
# a Linear horizontal regridding scheme is used by the RAS
# see: https://code.metoffice.gov.uk/doc/ancil/ants/1.1/bin/ancil_general_regrid.html
# outputs:
#      qrclim.land.nc (surface albedo)
python ${ANTS_SRC_PATH}/ancil_general_regrid.py --ants-config ${ANTS_CONFIG} \
       ${source} --target-lsm ${target_lsm} -o ${output}

# ============================================================================
# config based on RAS (u-bu503) app/ancil_soil_albedo/rose-app.conf, walltime ~ 1 minute
source=${ANCIL_MASTER}/soil_albedo/classic/v3/soil_albedo.nc
soil_hydrology=${ANCIL_TARGET_PATH}/soil_hydrology.nc
output=${ANCIL_TARGET_PATH}/qrparm.soil_cci
vegfrac=${ANCIL_TARGET_PATH}/qrparm.veg.frac_cci
ANTS_CONFIG=${INPUT_PATH}/ancil_soil_albedo-app.conf

# ancil_soil_albedo is a CONTRIB app that generates soil albedo from a global [source]
# to the target lct grid [vegfrac]. It also sets soil albedo to 0.75 where there is ice.
# more documentation is available in the file (no online documentation)
# outputs:
#    soil_albedo.nc (intermediate file)
python ${ANTS_SRC_PATH}/Apps/SoilAlbedo/ancil_soil_albedo.py ${source} \
       -o ${ANCIL_TARGET_PATH}/soil_albedo  \
       --use-new-saver --ants-config ${ANTS_CONFIG} \
       --lct-ancillary ${vegfrac}

# append.py is a RAS suite app that appends soil albedo to the soil parameter ancillary
# more documentation is available in the file (no online documentation)
# outputs:
#      qrparm.soil_cci.nc (appending soil_hydrology with soil_albedo)
python ${ANTS_SRC_PATH}/app/ancil_soil_albedo/bin/append.py \
       ${soil_hydrology} ${ANCIL_TARGET_PATH}/soil_albedo.nc -o ${output}

# ============================================================================
# config based on RAS (u-bu503) app/ancil_soil_dust/rose-app.conf, walltime ~ 5 seconds
source=${ANCIL_TARGET_PATH}/qrparm.soil_cci
vegfrac=${ANCIL_TARGET_PATH}/qrparm.veg.frac_cci
output=${ANCIL_TARGET_PATH}/qrparm.soil.dust
ANTS_CONFIG=${INPUT_PATH}/ancil_lai-app.conf

# ancil_soil_dust creates 6 x mass fraction and 3 x volume fraction dust 
# ancillaries based on clay silt and sand fractions.
# more documentation is available in the file (no online documentation)
# outputs: 
#   qrparm.soil.dust.nc
python ${ANTS_SRC_PATH}/Apps/SoilDust/ancil_soil_dust.py ${source} \
       -o ${output} --use-new-saver --ants-config ${ANTS_CONFIG} \
       --lct-ancillary ${vegfrac}

# ============================================================================
# config based on RAS (u-bu503) app/ancil_soil_roughness/rose-app.conf, walltime ~ 3 mins
# NOTE: Siyuan's offline code is different, it appends soil roughness to the soil parameter ancillary
#      therefore must be the last dependency (as here)
ANCIL_CONFIG=${INPUT_PATH}/ancil_lai-app.conf
source=${ANCIL_MASTER}/soil_roughness/prigent12/v1/soil_roughness_preproc_prigent12.nc
target_lai=${ANCIL_MASTER}/vegetation/lai/modis_4km/v2/lai_preprocessed.nc
target_lsm=${ANCIL_TARGET_PATH}/qrparm.mask
output=${ANCIL_TARGET_PATH}/qrparm.soilz0
merged_output=${ANCIL_TARGET_PATH}/qrparm.soil_cci.merged

# ancil_soil_roughness.py is a CONTRIB app that calculates soil roughness based on LAI
# more documentation is available in the file (no online documentation)
# outputs: 
#      qrparm.soilz0.nc
python ${ANTS_SRC_PATH}/Apps/SoilRoughness/ancil_soil_roughness.py \
       --target-lsm ${target_lsm} \
       --leaf-area-index ${target_lai} \
       --output ${output} \
       --use-new-saver ${source}

# # append.py appends the soil parameter ancillary [qrparm.soil_cci.nc] with 
# # soil roughness [output]. This is done in Siyuan's offline ancil suite, but not
# # done in the RAS suite. Therefore the merged and unmerged versions are retained.
# # more documentation is available in the file (no online documentation)
# python ${ANTS_SRC_PATH}/app/ancil_soil_roughness/bin/append.py \
#        ${ANCIL_TARGET_PATH}/qrparm.soil_cci.nc ${output}.nc -o ${merged_output}

# remove intermediate or unused files
rm ${ANCIL_TARGET_PATH}/qrparm.landfrac.nc
rm ${ANCIL_TARGET_PATH}/canopy_heights*
rm ${ANCIL_TARGET_PATH}/lai*
rm ${ANCIL_TARGET_PATH}/c4_percent_1d*
rm ${ANCIL_TARGET_PATH}/qrparm.veg.frac_cci_pre_c4*
rm ${ANCIL_TARGET_PATH}/qrparm.mask_sea*
rm ${ANCIL_TARGET_PATH}/soil_albedo*
rm ${ANCIL_TARGET_PATH}/soil_hydrology*
