#!/bin/bash
DBG=false

# ----------------------------------------------------------
if [[ $DBG == true ]]; then
  printf '%s\n' "[ starting FO_PRESIM.sh ]"
fi

# ----------------------------------------------------------
if [[ -z $1 ]] && [[ -z $2 ]] && [[ -z $3 ]]; then
  WORK_DIR=\
  "/work/bellout/WORK-3/SCRATCH_RUNS-PROJECTS/"\
  "P10-SLB-run-001-test-icd-opt-case-eb-20190611/xdir/"\
  "fo.r001.c02-olymr37.icd-opt.apps/rx00.BC00/test/R37/F37/W01"
  OPTDIMS_STR="7 7 /"
  #GRD_MODE="NoG"
  GRD_MODE="AdjG"
  MSW_MODE=""
  # OPTMODEL="5spot"
  OPTMODEL="olympr37"
  VAR_TEMP=""
  ACTIONX=false
  FLSTRUCT_TYPE="Flat"
  FLSTRUCT_LVL="0"
  FLSTRUCT_STR=".\/"
  FLSTRUCT_STR=`echo "${FLSTRUCT_STR}" | sed 's#\/#\\\/#g'`

else
  WORK_DIR="${1}"
  OPTDIMS_STR="${2}"
  GRD_MODE="${3}"
  MSW_MODE="${4}"
  OPTMODEL="${5}"
  VAR_TEMP="${6}"
  ACTIONX="${7}"
  FLSTRUCT_TYPE="${8}"
  FLSTRUCT_LVL="${9}"
  FLSTRUCT_STR="${10}"
  FLSTRUCT_STR=`echo "${FLSTRUCT_STR}" | sed 's#\/#\\\/#g'`

fi

# ----------------------------------------------------------
if [[ ${OPTMODEL} == "5spot" ]]; then
  OPTFILE=(`ls ${WORK_DIR}/*OPT.INC`)
elif [[ ${OPTMODEL} == "olympr37" ]]; then
  OPTFILE=(`ls ${WORK_DIR}/OLPZ_BCXX*.INC`)
fi

DATFILE=(`ls ${WORK_DIR}/*.DATA`)
DBGFILE=${WORK_DIR}/FO_PRESIM.DBG

# ----------------------------------------------------------
if [[ ${FLSTRUCT_TYPE} == "Flat" ]] && [[ ${OPTMODEL} == "5spot" ]]; then
  ROOTNME=${DATFILE##*/}
  ROOTNME=${ROOTNME%*.DATA}
  SCHDNME=${ROOTNME}.SCH

elif [[ ${FLSTRUCT_TYPE} == "Branched" ]] && [[ ${OPTMODEL} == "olympr37" ]]; then
  SCHDNME=(`ls ${WORK_DIR}/${FLSTRUCT_STR}*SCH.INC`)
  SCHDNME=${SCHDNME##*/}
fi

if [[ $DBG == true ]]; then
  echo "OPT=${OPTFILE}" > ${DBGFILE}
  echo "DAT=${DATFILE}" >> ${DBGFILE}
  echo "SCH=${SCHDNME}" >> ${DBGFILE}
  echo "DBG=${DBGFILE}" >> ${DBGFILE}

  echo "OPTDIMS_STR=${OPTDIMS_STR}" >> ${DBGFILE}
  echo "GRD_MODE=${GRD_MODE}" >> ${DBGFILE}
  echo "MSW_MODE=${MSW_MODE}" >> ${DBGFILE}
  echo "OPTMODEL=${OPTMODEL}" >> ${DBGFILE}
  echo "VAR_TEMP=${VAR_TEMP}" >> ${DBGFILE}
  echo "ACTIONX=${ACTIONX}" >> ${DBGFILE}

  echo "FLSTRUCT_TYPE=${FLSTRUCT_TYPE}" >> ${DBGFILE}
  echo "FLSTRUCT_LVL=${FLSTRUCT_LVL}" >> ${DBGFILE}
  echo "FLSTRUCT_STR=${FLSTRUCT_STR}" >> ${DBGFILE}

  echo "OPTDIMS_STR: ${OPTDIMS_STR[@]}"
  echo "GRD_MODE: ${GRD_MODE[@]}"
  echo "MSW_MODE: ${MSW_MODE[@]}"
  echo "OPTMODEL: ${OPTMODEL[@]}"
  echo "VAR_TEMP: ${VAR_TEMP[@]}"
  echo "ACTIONX: ${ACTIONX[@]}"

  echo "FLSTRUCT_TYPE: ${FLSTRUCT_TYPE[@]}"
  echo "FLSTRUCT_LVL: ${FLSTRUCT_LVL[@]}"
  echo "FLSTRUCT_STR: ${FLSTRUCT_STR[@]}"
fi

# ----------------------------------------------------------
OPTDIMS_LN=(`grep -n OPTDIMS ${OPTFILE} |\
           awk -F[':'] '{ print (($1 + 1))}'`)

if [[ $DBG == true ]]; then
  echo "OPTDIMS_LN=${OPTDIMS_LN}" >> ${DBGFILE}
  echo "OPTDIMS_LN: ${OPTDIMS_LN[@]}"
fi

# ==========================================================
# RUNSPEC MODES

# ----------------------------------------------------------
EXST_ACTDIMS=(`grep ACTDIMS ${DATFILE} | tr -d "\n"`)

if [[ ${EXST_ACTDIMS} == "" ]] && [[ ${ACTIONX} == true ]]; then
  CIN="BLACKOIL"; \
  COU="BLACKOIL\n\nACTDIMS\n50 200 \/"; \
  sed -i 's/'"$CIN"'/'"$COU"'/g' ${DATFILE}
fi

# ==========================================================
# SCHED MODES

# ----------------------------------------------------------
EXST_FOEDITS=(`grep "fo_edits.INC" ${DATFILE} | tr -d "\n"`)

if [[ $DBG == true ]]; then
  echo "EXST_FOEDITS=${EXST_FOEDITS}" >> ${DBGFILE}
  echo "EXST_FOEDITS: ${EXST_FOEDITS[@]}"
fi

if [[ ${EXST_FOEDITS} == "" ]]; then

  if [[ ${ACTIONX} == true ]]; then

    # find line of SCH file (2 lines after)
    SCHED_LN=(`grep -n ${SCHDNME} ${DATFILE} |\
              awk -F[':'] '{ print (($1 + 2))}'`)

    # introduce ACTIONX keywords right after SCH file
    LINE="-- _________________________________________________________"
    FOEDIT_STR="\n${LINE}\nINCLUDE\n\'fo_edits.INC\'\/\n";

  elif [[ ${ACTIONX} == false ]] && [[ ${MSW_MODE} == "MSW" ]]; then

    # find line of SCH file in DATA file
    SCHED_LN=(`grep -n ${SCHDNME} ${DATFILE} |\
              awk -F[':'] '{ print (($1 - 0))}'`)

    # replace SCH file with fo_edit file
    if [[ ${FLSTRUCT_TYPE} == "Flat" ]]; then
      FOEDIT_STR="\'fo_edits.INC\'\/\n";
    elif [[ ${FLSTRUCT_TYPE} == "Branched" ]]; then
      FOEDIT_STR="\'${FLSTRUCT_STR}fo_edits.INC\'\/\n";
    fi

  fi

  if [[ $DBG == true ]]; then
    echo "SCHED_LN=${SCHED_LN}" >> ${DBGFILE}
    echo "FOEDIT_STR: ${FOEDIT_STR[@]}"
  fi

  # replace line
  sed -i "${SCHED_LN}"'s/.*/'"${FOEDIT_STR}"'/' ${DATFILE}

fi

if [[ $DBG == true ]]; then
  echo "fo_edits.INC replacement done." >> ${DBGFILE}
  echo "fo_edits.INC replacement done."
fi

# ==========================================================
# G MODES
if [[ ${GRD_MODE} == "NoG" ]]; then


  EXST_RESOPT=(`grep -- "--RESOPT" ${DATFILE} | tr -d "\n"`)

  if [[ ${EXST_RESOPT} == "" ]]; then
    CIN="RESOPT"; \
    COU="--RESOPT"; \
    sed -i 's/'"$CIN"'/'"$COU"'/g' ${DATFILE}
  fi

  EXST_OPTIMIZE=(`grep "OPTIMIZE" ${DATFILE} -B2 | tr -d "\n"`)

  if [[ ${EXST_OPTIMIZE} != "ENDOPTIMIZE" ]]; then
    CIN="OPTIMIZE"; \
    COU="END\n\nOPTIMIZE"; \
    sed -i 's/'"$CIN"'/'"$COU"'/g' ${DATFILE}
  fi

  # insert OPTDIMS string
  sed -i "${OPTDIMS_LN}"'s/.*/'"${OPTDIMS_STR}"'/' ${OPTFILE}

elif [[ ${GRD_MODE} == "AdjG" || ${GRD_MODE} == "Optz" ]]; then

  sed -i "${OPTDIMS_LN}"'s/.*/'"${OPTDIMS_STR}"'/' ${OPTFILE}

fi
