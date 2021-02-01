#!/bin/bash
# runscript.r001.c04-5spot-flow.wbhp-opt.trdfo-test.mb.sh
# clear

# +---------------------------------------------------------
. ../../proc/sh-post/set_color_functions.sh

# +---------------------------------------------------------
# Dirs
MAIN_DIR=~/WORK-3/SCRATCH_RUNS-PROJECTS
PROJ_DIR=${MAIN_DIR}/P18-TLS-tr-dfo-method-dev_scratch-runs
XDIR_DIR=${PROJ_DIR}/xdir
MODS_DIR=${PROJ_DIR}/mods
DRVR_DIR=${MODS_DIR}/drivers

# +---------------------------------------------------------
# FORE_DIR=~/git/PCG/FieldOpt-Research/FieldOpt
FORE_DIR=~/git/IOC/FieldOpt-Research/FieldOpt
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

MODL=ECL_5SPOT_C04_FLOW
MODL_DIR=${MODS_DIR}/5spot_c04_flow
DECK_PTH=${MODL_DIR}/${MODL}.DATA
GRID_PTH=${MODL_DIR}/${MODL}.EGRID
XDIR_PTH=${XDIR_DIR}/20191106_${PROB}

# +---------------------------------------------------------
# Dir checks
if [ ! -d ${PROJ_DIR} ]; then pred "Project dir not found at ${PROJ_DIR}"; fi
if [ ! -d ${XDIR_DIR} ]; then pred "X dir not found at ${XDIR_DIR}"; fi
if [ ! -d ${MODS_DIR} ]; then pred "Mods dir not found at ${MODS_DIR}"; fi
if [ ! -d ${DRVR_DIR} ]; then pred "Driver dir not found at ${DRVR_DIR}"; fi

if [ ! -d ${FORE_DIR} ]; then pred "FoRe dir not found at ${FORE_DIR}"; fi
if [ ! -f ${FORE_EXE} ]; then pred "Fo exec dir not found at ${FORE_EXE}"; fi
if [ ! -f ${SIMS_EXE} ]; then pred "Sim exec dir not found at ${SIMS_EXE}"; fi

if [ ! -f ${DRVR_PTH} ]; then pred "JSON driver file not found at ${DRVR_PTH}"; fi
if [ ! -f ${DECK_PTH} ]; then pred "Deck file not found at ${DECK_PTH}"; fi
if [ ! -f ${GRID_PTH} ]; then pred "Grid file not found at ${GRID_PTH}"; fi

# +---------------------------------------------------------
LINE=`printf '=%.0s' {1..80}`

pcyn "${LINE}"
pcyn "PROJ_DIR: ${PROJ_DIR}"
pcyn "XDIR_DIR: ${XDIR_DIR}"
pcyn "MODS_DIR: ${MODS_DIR}"
pcyn "DRVR_DIR: ${DRVR_DIR}"

pcyn "${LINE}"
pcyn "FORE_DIR: ${FORE_DIR}"
pcyn "FORE_EXE: ${FORE_EXE}"
pcyn "SIMS_EXE: ${SIMS_EXE}"

pcyn "${LINE}"
pcyn "DRVR_PTH: ${DRVR_PTH}"
pcyn "DECK_PTH: ${DECK_PTH}"
pcyn "GRID_PTH: ${GRID_PTH}"

# +---------------------------------------------------------
rm -rf ${XDIR_PTH}/*
mkdir -p ${XDIR_PTH}

# +---------------------------------------------------------
if [[ ${MPIX} == "" ]]; then

  ${FORE_EXE} ${DRVR_PTH} ${XDIR_PTH} \
  -r serial -g ${GRID_PTH} -s ${DECK_PTH} -f \
  -e ${SIMS_EXE} --sim-aux ${MODL_DIR}

# elif [[ ${MPIX} == "mpiexec" ]]; then

# ${MPICMD} -n ${N_PROCS} \
# ${FIELDOPT_BIN} ${DRIVER_PATH} ${OUTPUT_DIR} \
# -r mpisync -t5 -g ${GRID_PATH} -s ${DECK_PATH} \
# -e ${SCR_PATH} --sim-aux ${AUX_DIR}

fi