#!/bin/bash
# Parameter 1: Work direcory.
# Parameter 2: Path to .DATA file

DBG=false
STR_CUT=156

# ----------------------------------------------------------
if [[ $DBG == true ]]; then
  printf "\n%s\n"   "[bash_e300_slbs.sh] Running simulation"
  printf "%s\n"   "WDIR: $1"  | fold -w ${STR_CUT}
  printf "%s\n\n" "DATA FILE NAME: $2" | fold -w ${STR_CUT}
else
  printf ". "
fi

# Switch to work directory
cd $1

# Execute eclipse with the file path as parameter
# echo "IN FieldOpt-Research-Open/FieldOpt/Simulation/execution_scripts/bash_e300_slbs.sh:"
echo "exec INSERT_THE_FULL_PATH_OF_YOUR_ECLIPSE300_EXECUTABLE_HERE $2 > sim.out"

echo 1 >> nexec.out

if [[ $DBG == true ]]; then
  printf "\n"
fi