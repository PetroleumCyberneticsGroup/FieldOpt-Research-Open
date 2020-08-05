#!/bin/bash
# make_5spot_rlzs_perm.sh

NRLZ=10
RSRC=../data

rm -rf ../r0*

# ----------------------------------------------------------
for (( R=1; R<$(( ${NRLZ}+1 )); R++ )); do

  NR=`printf "%03.0f" ${R}`
  RTRG=../r${NR}

  mkdir ${RTRG}
  # rm -rf ../${RTRG}/*

  cp -r ${RSRC}/*.DATA ${RSRC}/*.EGRID ${RTRG}
  cp -r ${RSRC}/*.INIT ${RSRC}/*.UNSMRY ${RTRG}
  cp -r ${RSRC}/fo_edits.INC ${RTRG}

  # --------------------------------------------------------  
  cd ${RTRG}

  DAT_FILE0=`ls *.DATA`
  DAT_FILE1=${DAT_FILE0%__*}__R${NR}.DATA
  mv ${DAT_FILE0} ${DAT_FILE1}

  # --------------------------------------------------------
  PX_STR0="PERMX    1.0"
  PY_STR0="PERMY    1.0"
  PZ_STR0="PERMZ    1.0"

  NP=`echo "1.0 + 1/${NRLZ}*${R}" | bc -l`
  NP=`printf "%2.1f" ${NP}`

  PX_STR1="PERMX    "$NP
  PY_STR1="PERMY    "$NP
  PZ_STR1="PERMZ    "$NP

  sed -i 's/'"${PX_STR0}"'/'"${PX_STR1}"'/g' ${DAT_FILE1}
  sed -i 's/'"${PY_STR0}"'/'"${PY_STR1}"'/g' ${DAT_FILE1}
  sed -i 's/'"${PZ_STR0}"'/'"${PZ_STR1}"'/g' ${DAT_FILE1}

  echo "---------------------------------------------------"
  echo "PERM at ${PWD##*/}:"    
  grep PERM ${DAT_FILE1}

  # --------------------------------------------------------
  EGRID_FILE0=`ls *.EGRID`
  EGRID_FILE1=${EGRID_FILE0%__*}__R${NR}.EGRID
  mv ${EGRID_FILE0} ${EGRID_FILE1}

  INIT_FILE0=`ls *.INIT`
  INIT_FILE1=${INIT_FILE0%__*}__R${NR}.INIT
  mv ${INIT_FILE0} ${INIT_FILE1}

  UNSMRY_FILE0=`ls *.UNSMRY`
  UNSMRY_FILE1=${UNSMRY_FILE0%__*}__R${NR}.UNSMRY
  mv ${UNSMRY_FILE0} ${UNSMRY_FILE1}

done

