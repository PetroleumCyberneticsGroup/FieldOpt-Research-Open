#!/bin/bash
# Parameter 1: Work direcory.
# Parameter 2: Path to .DATA file

FOUT=sim_exe.out
# echo "Running simulation"
echo "Work directory: " $1 > $1/${FOUT}
echo "Data file path: " $2 >> $1/${FOUT}

# Switch to work directory
cd $1

# exec /opt/ecl/2013.2/bin/linux_x86_64/eclipse.exe $2 >/dev/null
exec /opt/ecl/ecl-2013.2/macros/@eclipse -file $2 -local >/dev/null
# exec /opt/ecl/ecl-2013.2/bin/linux_x86_64/eclipse.exe $2 >/dev/null

STR="Errors                 0"
echo $STR
grep "${STR}" ${1}/*.PRT >> $1/${FOUT}