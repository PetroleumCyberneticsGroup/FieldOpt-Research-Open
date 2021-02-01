#!/bin/bash
# test_5spot_rlzs_perm.sh

RDIRS=(`ls -d ../r0*`)

# ----------------------------------------------------------
for RDIR in "${RDIRS[@]}"; do

  cd $RDIR
  DATA_FILE=`ls *.DATA`
  CMD_STR="flow ${DATA_FILE}"
  # eval ${CMD_STR}

  # --------------------------------------------------------
  RSLT_FILE=`ls *.PRT`

  # DATA AT REPORT STEP 1
  STR0="Balance  at     182.5  Days"
  STR1="Porv volumes are taken at reference conditions"
  STR2="Currently   in place"
  BEG_STR=`grep "${STR0}" -B23 ${RSLT_FILE} \
  | grep "${STR1}" -A4 | grep "${STR2}"`

  # DATA AT END
  STR0="Balance  at      2555  Days"
  STR1="Porv volumes are taken at reference conditions"
  STR2="Currently   in place"
  END_STR=`grep "${STR0}" -B23 ${RSLT_FILE} \
  | grep "${STR1}" -A4 | grep "${STR2}"`

  # echo "BEG_STR: $BEG_STR"
  # echo "END_STR: $END_STR"

  # --------------------------------------------------------
  TOTALS_END=`echo "${END_STR}" \
  | awk -F[':'] '{ print $3 $4 }'`

  # echo "TOTALS_END=${TOTALS_END}"

  FOPT_END=`echo "${TOTALS_END}" \
  | awk '{ print $3 }'`

  FWPT_END=`echo "${TOTALS_END}" \
  | awk '{ print $4 }'`

  echo "${RDIR}: FOPT_END=${FOPT_END} | FWPT_END=${FWPT_END}"

done 