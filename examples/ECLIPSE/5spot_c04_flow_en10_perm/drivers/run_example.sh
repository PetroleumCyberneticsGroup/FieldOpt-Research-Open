#!/bin/bash
# runscript.r001.c04-5spot-flow.wbhp-opt.trdfo-test.mb.sh
# clear

# +---------------------------------------------------------
#. ../../proc/sh-post/set_color_functions.sh

# +---------------------------------------------------------
# Dirs
MAIN_DIR=~/projects/PCG/FieldOpt-Research/FieldOpt/cmake-build-debug/examples/ECLIPSE
PROJ_DIR=${MAIN_DIR}/5spot_c04_flow_en10_perm
XDIR_DIR=~/projects/PCG/FieldOpt-Research/fieldopt-output
MODS_DIR=${PROJ_DIR}
DRVR_DIR=${MODS_DIR}/drivers
ENS_PATH=${PROJ_DIR}/en_5spot.ens

# +---------------------------------------------------------
# FORE_DIR=~/git/PCG/FieldOpt-Research/FieldOpt
FORE_DIR=~/projects/PCG/FieldOpt-Research/FieldOpt
FORE_EXE=${FORE_DIR}/cmake-build-debug/bin/FieldOpt
# SIMS_EXE=${DRVR_DIR}/bash_e300.sh
SIMS_EXE=${DRVR_DIR}/bash_flow.sh
# SIMS_EXE=${DRVR_DIR}/bash_ecl.sh

# +---------------------------------------------------------
# Case vars
if [[ ${1} == "" ]]; then
  # PROB=r001.c04-5spot-flow.wbhp-opt.trdfo-t01
  PROB=r001.c04-5spot-flow.wbhp-opt.apps-t01
  DRVR_PTH=${DRVR_DIR}/fo-drv.${PROB}.npv.json
else
  PROB=${1}
  DRVR_PTH=${DRVR_DIR}/fo-drv.${PROB}.npv.json
fi

MPIX=""
NPRC=2

MODL=ECL_5SPOT_C04_FLOW__R001
MODL_DIR=${MODS_DIR}/include
DECK_PTH=${PROJ_DIR}/r001/${MODL}.DATA
GRID_PTH=${PROJ_DIR}/r001/${MODL}.EGRID
XDIR_PTH=${XDIR_DIR}/${PROB}

# +---------------------------------------------------------
# Dir checks
if [ ! -d ${PROJ_DIR} ]; then echo "Project dir not found at ${PROJ_DIR}"; fi
if [ ! -d ${XDIR_DIR} ]; then echo "X dir not found at ${XDIR_DIR}"; fi
if [ ! -d ${MODS_DIR} ]; then echo "Mods dir not found at ${MODS_DIR}"; fi
if [ ! -d ${DRVR_DIR} ]; then echo "Driver dir not found at ${DRVR_DIR}"; fi

if [ ! -d ${FORE_DIR} ]; then echo "FoRe dir not found at ${FORE_DIR}"; fi
if [ ! -f ${FORE_EXE} ]; then echo "Fo exec dir not found at ${FORE_EXE}"; fi
if [ ! -f ${SIMS_EXE} ]; then echo "Sim exec dir not found at ${SIMS_EXE}"; fi

if [ ! -f ${DRVR_PTH} ]; then echo "JSON driver file not found at ${DRVR_PTH}"; fi
if [ ! -f ${DECK_PTH} ]; then echo "Deck file not found at ${DECK_PTH}"; fi
if [ ! -f ${GRID_PTH} ]; then echo "Grid file not found at ${GRID_PTH}"; fi

# +---------------------------------------------------------
LINE=`printf '=%.0s' {1..80}`

echo "${LINE}"
echo "PROJ_DIR: ${PROJ_DIR}"
echo "XDIR_DIR: ${XDIR_DIR}"
echo "MODS_DIR: ${MODS_DIR}"
echo "DRVR_DIR: ${DRVR_DIR}"

echo "${LINE}"
echo "FORE_DIR: ${FORE_DIR}"
echo "FORE_EXE: ${FORE_EXE}"
echo "SIMS_EXE: ${SIMS_EXE}"

echo "${LINE}"
echo "DRVR_PTH: ${DRVR_PTH}"
echo "DECK_PTH: ${DECK_PTH}"
echo "GRID_PTH: ${GRID_PTH}"

# +---------------------------------------------------------
rm -rf ${XDIR_PTH}/*
mkdir -p ${XDIR_PTH}

# +---------------------------------------------------------
if [[ ${MPIX} == "" ]]; then

  ${FORE_EXE} ${DRVR_PTH} ${XDIR_PTH} \
  -r serial -g ${GRID_PTH} -s ${DECK_PTH} -f \
  -e ${SIMS_EXE} --sim-aux ${MODL_DIR} --ensemble-path ${ENS_PATH}

# elif [[ ${MPIX} == "mpiexec" ]]; then

# ${MPICMD} -n ${N_PROCS} \
# ${FIELDOPT_BIN} ${DRIVER_PATH} ${OUTPUT_DIR} \
# -r mpisync -t5 -g ${GRID_PATH} -s ${DECK_PATH} \
# -e ${SCR_PATH} --sim-aux ${AUX_DIR}

fi
