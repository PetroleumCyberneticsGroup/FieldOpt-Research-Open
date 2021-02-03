#!/bin/bash
# Parameter 1: Work direcory.
# Parameter 2: Path to .DATA file

printf "%s\n"   "[bash_e300_slbs.sh] Running simulation"
printf "%s\n"   "WDIR: $1"
printf "%s\n\n" "DATA FILE NAME: $2"

# Switch to work directory
cd $1

HOSTN=$(hostname | tr "/" "\n")

if [[ "${HOSTN}" == "nextron-SYS-7049GP-TRT" ]]; then
  E300=/xhome/bellout/loc_opt/ecl/2013.2/2013.2/bin/linux_x86_64/e300.exe

elif [[ "${HOSTN}" == "bellout-X1" ]]; then
  E300="/opt/ecl/2013.2/bin/linux_x86_64/e300.exe"

else
  echo "[bash_e300_slbs.sh] Provide path to e300"
  E300=""
fi

if [[ ${E300} != "" ]]; then
  XCMD="${E300} ${2} >/dev/null"
  exec ${XCMD}
fi