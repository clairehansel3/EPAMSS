#!/usr/bin/env bash
#$ -cwd
#$ -o joblog.$JOB_ID
#$ -l h_rt=24:00:00,h_data=4G
#$ -M clairehansel3@gmail.com
#$ -m bea

echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `

. /u/local/Modules/default/init/modules.sh
export LD_LIBRARY_PATH="/u/home/c/claireha/prefix/lib:$LD_LIBRARY_PATH"
module unload gcc/4.9.3
module load gcc/7.2.0
module load openmpi/3.1.3
module unload intel/13.cs
mpirun -n $NSLOTS /u/home/c/claireha/EPAMSS/epamss $1

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
