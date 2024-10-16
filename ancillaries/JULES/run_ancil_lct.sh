#!/bin/bash
#PBS -q normal
#PBS -l ncpus=1
#PBS -l mem=4GB
#PBS -l walltime=01:00:00
#PBS -l storage=gdata/hh5+gdata/access+gdata/rp23
#PBS -l wd
#PBS -l jobfs=1GB
#PBS -P rp23

module purge
module use /g/data/access/ngm/modules
module load ants/1.1.0  # as used by Siyuan's code, also works with the RAS-based code
# module load ants/0.18  # as used in RAS

ANTS_SRC_PATH=/g/data/rp23/experiments/2024-10-10_LWG_workingbee/mjl561/ants_src/ants_1p1
ANCIL_MASTER=/g/data/access/TIDS/UM/ancil/atmos/master
INPUT_PATH=/g/data/rp23/experiments/2024-10-10_LWG_workingbee/mjl561/ants_lct_inputs
ANCIL_TARGET_PATH=/g/data/rp23/experiments/2024-10-10_LWG_workingbee/mjl561/outputs_1p1

# check if the output directory exists
if [ ! -d ${ANCIL_TARGET_PATH} ]; then
    mkdir -p ${ANCIL_TARGET_PATH}
fi

# based on RAS (u-bu503) app/ancil_lct/rose-app.conf
source=${ANCIL_MASTER}/vegetation/cover/cci/v3/vegetation_fraction.nc
target_grid=${INPUT_PATH}/grid.nl
transformpath=/g/data/access/TIDS/UM/ancil/data/transforms/cci2jules_ra1.json
output_vegfrac=${ANCIL_TARGET_PATH}/qrparm.veg.frac_cci_pre_c4
output_lsm=${ANCIL_TARGET_PATH}
ANTS_CONFIG=${INPUT_PATH}/ancil_lct-app.conf

python ${ANTS_SRC_PATH}/ancil_lct.py ${source} \
         --target-grid ${target_grid} --transform-path ${transformpath} \
         -o ${output_vegfrac} --landseamask-output ${output_lsm}       \
         --ants-config ${ANTS_CONFIG}

# based on RAS app/ancil_lct_postproc_c4/rose-app.conf
ANCIL_CONFIG=${INPUT_PATH}/ancil_lct_postproc_c4-app.conf
source=${ANCIL_TARGET_PATH}/qrparm.veg.frac_cci_pre_c4.nc
target_lsm=${ANCIL_TARGET_PATH}/qrparm.mask
output=${ANCIL_TARGET_PATH}/qrparm.veg.frac_cci
c4source=${ANCIL_MASTER}/vegetation/cover/cci/v3/c4_percent_1d.nc
# c4source=/g/data/dp9/mjl561/au_ancils_source/Oceania_CSIRO_ISLSCP_C4_percent_300m.nc # high res version by Luigi

python ${ANTS_SRC_PATH}/ancil_general_regrid.py --ants-config ${ANCIL_CONFIG} \
       ${c4source} --target-lsm ${target_lsm} -o ${ANCIL_TARGET_PATH}/c4_percent_1d.nc

python ${ANTS_SRC_PATH}/ancil_lct_postproc_c4.py \
       ${source} --islscpiic4 ${ANCIL_TARGET_PATH}/c4_percent_1d.nc -o ${output}
