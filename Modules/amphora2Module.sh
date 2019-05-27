#!/bin/bash 
#SBATCH -J AMPH2P1
#SBATCH -o amphora2_%a.out
#SBATCH -e amphora2_%a.err
#SBATCH -t AMPHORA2TIME
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=USERMAIL
#SBATCH --partition=AMPHORA2PARTITIONS
#SBATCH --nodes=1
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


in=$(sed -n "$SLURM_ARRAY_TASK_ID"p samples.txt)
bac=$(echo "$in" |awk '{print $1}')
jobid=$(echo "$in" |awk '{print $2}')
bacname=$(echo "$bac" |rev |cut -d "/" -f 1 |rev)

mkdir amphora_$bacname
cd amphora_$bacname
	PERLBIN AMPHORA2HOME/Scripts/MarkerScanner.pl -DNA -Bacteria $bac
	PERLBIN AMPHORA2HOME/Scripts/MarkerAlignTrim.pl -WithReference -OutputFormat fasta
	#PERLBIN AMPHORA2HOME/Scripts/Phylotyping.pl -CPUs $CORES > ${bacname}.result
cd ..


echo "$jobid" >> AMPHORA2P1_JOBS.txt

if [ $((jobid)) -eq GNUMS ];then

	completedjobs=$(wc -l AMPHORA2P1_JOBS.txt |awk '{print $1}')
	if [ $((completedjobs)) -eq GNUMS ];then
		chainReaction "AMPH2P1" $1 $DEBUG
	else
		sbatch CALLSLURM "AMPH2P1" "AMPHORA2P1_JOBS.txt" GNUMS $1  $DEBUG
	fi
fi

	