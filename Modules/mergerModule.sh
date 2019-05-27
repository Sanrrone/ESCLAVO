#!/bin/bash 
#SBATCH -J mergeMod
#SBATCH -o merge_%a.out
#SBATCH -e merge_%a.err
#SBATCH -t MERGERTIME
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=USERMAIL
#SBATCH --partition=MERGERPARTITION
#SBATCH -N 1
#SBATCH --array=1-GNUMS

if [[ "$@" =~ "--debug" ]]; then
		DEBUG="--debug"
        set -ex
else
        set -e
        DEBUG=""
fi

source CHAINREACTION

CORES=$(grep -c ^processor /proc/cpuinfo)

i=$(sed -n "$SLURM_ARRAY_TASK_ID"p INPUTREADSLIST)
R1=$(echo "$i" |awk '{print $1}')
R2=$(echo "$i" |awk '{print $2}')
jobid=$(echo "$i" |awk '{print $3}')

name=$(echo "$R1" |rev |cut -d "/" -f 1 |rev |sed "s/\.fastq//g")

PEARBIN -f $R1 -r $R2 -o $name -q PEARQ -u 0.1 -y PEARMEMORY -j $CORES PEAREXTRACODE

echo "$jobid" >> MERGE_JOBS.txt

if [ $((jobid)) -eq GNUMS ];then

	completedjobs=$(wc -l MERGE_JOBS.txt |awk '{print $1}')
	if [ $((completedjobs)) -eq GNUMS ];then
		echo "MERGE module ended, calling chainReaction Module"
		chainReaction "mergeMod" $1 $DEBUG
	else
		echo "MERGE module is still working, calling CALLSLURM Module"
		sbatch CALLSLURM "mergeMod" "MERGE_JOBS.txt" GNUMS $1  $DEBUG
	fi
fi

