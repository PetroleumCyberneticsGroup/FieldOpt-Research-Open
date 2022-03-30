#!/bin/bash
DBG=false
STR_CUT=156

# ----------------------------------------------------------
if [[ $DBG == true ]]; then
  printf '%s\n' "[ starting FO_POSTSIM.sh ]"
fi

# ----------------------------------------------------------
if [[ -z $1 ]] && [[ -z $2 ]] && [[ -z $3 ]]; then
  WORK_DIR=\
  "/work/bellout/WORK-3/SCRATCH_RUNS-PROJECTS/"\
  "P10-SLB-run-001-test-icd-opt-case-eb-20190611/xdir/"\
  "fo.r001.c02-olymr37.icd-opt.apps/rx00.BC00/test/R37/F37/W01"
  GRD_MODE=AdjG
  NUM_GRDLNS=8
else
  WORK_DIR="${1}"
  GRD_MODE="${2}"
  NUM_GRDLNS="${3}"
fi

if [[ $DBG == true ]]; then
  echo "WORK_DIR=${WORK_DIR}" | fold -w ${STR_CUT}
  echo "GRD_MODE=${GRD_MODE}" | fold -w ${STR_CUT}
  echo "NUM_GRDLNS=${NUM_GRDLNS}" | fold -w ${STR_CUT}
fi

NUM_TABLNS=$((${NUM_GRDLNS} + 2))
TLT_STR="ID  DOMAIN   TYPE    GRAD          TIME          NORM"

if [[ $DBG == true ]]; then
  echo "NUM_TABLNS=${NUM_TABLNS}" | fold -w ${STR_CUT}
fi

# ----------------------------------------------------------
. check-prt-lic-b.sh "${WORK_DIR}" false
DBGFILE=${WORK_DIR}/FO_POSTSIM.DBG

# ----------------------------------------------------------
DATFILE=`(ls ${WORK_DIR}/*.DATA)`
ECLGRAD=`(ls ${WORK_DIR}/*.DBG)`

GRDFILE=`(ls ${DATFILE})`
# GRDFILE=${GRDFILE[@]%%.DATA*}
GRDFILE=${GRDFILE%*.DATA}

GRDJSON=${GRDFILE}.JSON
FVLFILE=${GRDFILE}.FVAL
BCKGRAD=${GRDFILE}.GBCK
GRDFILE=${GRDFILE}.GRAD

# ----------------------------------------------------------
if [[ $DBG == true ]]; then
  echo "DAT=${DATFILE#*SCRATCH_RUNS-PROJECTS}" | fold -w ${STR_CUT}
  echo "EGR=${ECLGRAD#*SCRATCH_RUNS-PROJECTS}" | fold -w ${STR_CUT}
  echo "GRD=${GRDFILE#*SCRATCH_RUNS-PROJECTS}" | fold -w ${STR_CUT}
  echo "JSN=${GRDJSON#*SCRATCH_RUNS-PROJECTS}" | fold -w ${STR_CUT}
  echo "DBG=${DBGFILE#*SCRATCH_RUNS-PROJECTS}" | fold -w ${STR_CUT}

  echo "DAT=${DATFILE#*SCRATCH_RUNS-PROJECTS}" >  ${DBGFILE}
  echo "EGR=${ECLGRAD#*SCRATCH_RUNS-PROJECTS}" >> ${DBGFILE}
  echo "GRD=${GRDFILE#*SCRATCH_RUNS-PROJECTS}" >> ${DBGFILE}
  echo "FVL=${FVLFILE#*SCRATCH_RUNS-PROJECTS}" >> ${DBGFILE}
  echo "JSN=${GRDJSON#*SCRATCH_RUNS-PROJECTS}" >> ${DBGFILE}
  echo "DBG=${DBGFILE#*SCRATCH_RUNS-PROJECTS}" >> ${DBGFILE}
fi

printf "%s\n" "FVLFILE: ${FVLFILE}" > ${FVLFILE}
printf "%s\n" "BCKGRAD: ${BCKGRAD}" > ${BCKGRAD}

# ----------------------------------------------------------
if [[ ${GRD_MODE} == "AdjG" ]]; then

  # cp "${ECLGRAD}" "${BCKGRAD}"

  FUNCVAL=(`grep "FUNC VAL" "${ECLGRAD}" | awk '{ print $7 }' | sed 's/BETA//g' `)
  printf "FUNCVAL: % .6e\n" "${FUNCVAL}" >> ${FVLFILE}
  grep "FUNC VAL" "${ECLGRAD}" >> ${BCKGRAD}


  grep -m1 "${TLT_STR}" "${ECLGRAD}" -A${NUM_TABLNS} | \
  tail -n${NUM_GRDLNS} > ${GRDFILE}

  TAB="  "
  TTAB=${TAB}${TAB}
  printf "%s\n" "{" > ${GRDJSON}
  printf "${TAB}%s\n" "\"AdjG\": [" >> ${GRDJSON}

  while read -r LINE; do
    IFS=' ' read -r -a ROW <<< "${LINE}"

    # echo "${ROW[@]}"
    # echo "${ROW[3]}"
    # echo "${ROW[4]}"
    # echo "${ROW[5]}"

    printf "${TTAB}{ \"ID\": %s,\n" ${ROW[0]} >> ${GRDJSON}
    printf "${TTAB}  \"DOMAIN\": \"%s\",\n" ${ROW[1]} >> ${GRDJSON}
    printf "${TTAB}  \"TYPE\": \"%s\",\n" ${ROW[2]} >> ${GRDJSON}
    printf "${TTAB}  \"GRAD\": % .6e,\n" ${ROW[3]} >> ${GRDJSON}
    printf "${TTAB}  \"TIME\": % .6e,\n" ${ROW[4]} >> ${GRDJSON}
    printf "${TTAB}  \"NORM\": % .6e, \n" ${ROW[5]} >> ${GRDJSON}
    printf "${TTAB}  \"FVAL\": % .6e \n" ${FUNCVAL} >> ${GRDJSON}

    if [[ ${NUM_GRDLNS} == ${ROW[0]} ]]; then
      printf "${TTAB}}\n" >> ${GRDJSON}
    else
      printf "${TTAB}},\n" >> ${GRDJSON}
    fi

  done < ${GRDFILE}

  printf "${TAB}%s\n" "]" >> ${GRDJSON}
  printf "%s\n" "}" >> ${GRDJSON}

  # CIN="-."; \
  # COU="-0."; \
  # sed -i 's/'"$CIN"'/'"$COU"'/g' ${GRDJSON}

elif [[ ${GRD_MODE} == "NoG" ]]; then

  printf "%s\n" "${GRD_MODE}" > ${GRDJSON}

fi

# ID  DOMAIN   TYPE    GRAD          TIME          NORM


#  1  PRODX2:2 SCSA  -.395934E+04  0.109500E+02  -.361583E+03
#  2  PRODX2:3 SCSA  0.499825E+04  0.109500E+02  0.456461E+03
#  3  PRODX2:4 SCSA  0.413274E+04  0.109500E+02  0.377420E+03
#  4  PRODX2:5 SCSA  0.121217E+05  0.109500E+02  0.110700E+04
#  5  PRODX2:6 SCSA  -.453512E+02  0.109500E+02  -.414166E+01
#  6  PRODX2:7 SCSA  -.153142E+04  0.109500E+02  -.139856E+03
#  7  PRODX2:8 SCSA  -.536249E+05  0.109500E+02  -.489725E+04
#  8  PRODX2:9 SCSA  -.865249E+05  0.109500E+02  -.790182E+04



#   ID  DOMAIN   TYPE    GRAD          TIME          NORM


#  1  PRODX2   WBHP  -.249614E+04  0.109500E+02  -.227958E+03
#  2  PRODX2   WBHP  -.244644E+04  0.229950E+03  -.223419E+03
#  3  PRODX2   WBHP  -.262300E+04  0.448950E+03  -.239543E+03
#  4  PRODX2   WBHP  -.231894E+04  0.667950E+03  -.211776E+03
#  5  PRODX2   WBHP  -.229919E+04  0.886950E+03  -.209971E+03
#  6  PRODX2   WBHP  -.221349E+04  0.110595E+04  -.202145E+03
#  7  PRODX2   WBHP  -.189771E+04  0.132495E+04  -.173307E+03
#  8  PRODX2   WBHP  -.202795E+04  0.154395E+04  -.185201E+03
#  9  PRODX2   WBHP  -.250451E+04  0.176295E+04  -.228723E+03
# 10  PRODX2   WBHP  -.413506E+04  0.198195E+04  -.377631E+03
# 11  INJD-15  WBHP  -.591788E+04  0.109500E+02  -.540446E+03
# 12  INJD-15  WBHP  -.483444E+04  0.229950E+03  -.441501E+03
# 13  INJD-15  WBHP  -.373138E+04  0.448950E+03  -.340765E+03
# 14  INJD-15  WBHP  -.287487E+04  0.667950E+03  -.262545E+03
# 15  INJD-15  WBHP  -.236900E+04  0.886950E+03  -.216347E+03
# 16  INJD-15  WBHP  -.216151E+04  0.110595E+04  -.197398E+03
# 17  INJD-15  WBHP  -.214263E+04  0.132495E+04  -.195674E+03
# 18  INJD-15  WBHP  -.237729E+04  0.154395E+04  -.217104E+03
# 19  INJD-15  WBHP  -.246978E+04  0.176295E+04  -.225550E+03
# 20  INJD-15  WBHP  -.266788E+04  0.198195E+04  -.243642E+03
# 21  INJD-16  WBHP  0.500279E+04  0.109500E+02  0.456876E+03
# 22  INJD-16  WBHP  0.410778E+04  0.229950E+03  0.375140E+03
# 23  INJD-16  WBHP  0.290535E+04  0.448950E+03  0.265329E+03
# 24  INJD-16  WBHP  0.174916E+04  0.667950E+03  0.159741E+03
# 25  INJD-16  WBHP  0.109230E+04  0.886950E+03  0.997531E+02
# 26  INJD-16  WBHP  0.875021E+03  0.110595E+04  0.799106E+02
# 27  INJD-16  WBHP  0.970909E+03  0.132495E+04  0.886675E+02
# 28  INJD-16  WBHP  0.111049E+04  0.154395E+04  0.101414E+03
# 29  INJD-16  WBHP  0.150171E+04  0.176295E+04  0.137142E+03
# 30  INJD-16  WBHP  0.167414E+04  0.198195E+04  0.152890E+03
# 31  PRODX2:2 SCSA  -.395934E+04  0.109500E+02  -.361583E+03
# 32  PRODX2:3 SCSA  0.499825E+04  0.109500E+02  0.456461E+03
# 33  PRODX2:4 SCSA  0.413274E+04  0.109500E+02  0.377420E+03
# 34  PRODX2:5 SCSA  0.121217E+05  0.109500E+02  0.110700E+04
# 35  PRODX2:6 SCSA  -.453512E+02  0.109500E+02  -.414166E+01
# 36  PRODX2:7 SCSA  -.153142E+04  0.109500E+02  -.139856E+03
# 37  PRODX2:8 SCSA  -.536249E+05  0.109500E+02  -.489725E+04
# 38  PRODX2:9 SCSA  -.865249E+05  0.109500E+02  -.790182E+04

