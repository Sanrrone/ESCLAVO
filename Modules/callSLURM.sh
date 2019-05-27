#!/bin/bash 
#SBATCH -J CHECKMOD
#SBATCH -o checkmodule.out
#SBATCH -e checkmodule.err
#SBATCH -t 04:00:00
#SBATCH --partition=defq,short,debug
#SBATCH --nodes=1

if [[ "$@" =~ "--debug" ]]; then
		DEBUG="--debug"
        set -ex
else
        set -e
        DEBUG=""
fi

	title=$1
	jobfile=$2
	jobnums=$3
	configfile=$4

grep -v "^#" $configfile |sed "s/#.*//g" |awk '{if($0!="")print}' > tmpconf

while read parameter
do
	Pname=$(echo "$parameter" |awk -F"=" '{print $1}')	
	case $Pname in
		"ESCLAVOPATH")
			ESCLAVOPATH=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
	esac
done < <(grep "" tmpconf)
rm -f tmpconf

source $ESCLAVOPATH/Modules/chainReaction.sh


	checktime=0
	status="onjob"
	echo "waiting for $software results"
	while [ "$status" != "" ]
	do
	        sleep 60
	        checktime=$(echo $checktime |awk '{print $1 + 61}')
	        if [ $((checktime)) -ge 14000 ];then
	        	sbatch ${BASH_SOURCE[0]} $title $jobfile $jobnums $configfile $DEBUG
	        	exit
	        fi
	        status=$(squeue |grep "$title" |awk 'BEGIN{status=""}{if($0!=""){status="onjob"}}END{print status}')
	done
	
	completedjobs=$(wc -l $jobfile |awk '{print $1}')
	if [ $((completedjobs)) -eq $((jobnums)) ];then
		chainReaction $title $configfile $DEBUG
		exit
	else
		echo "$title finished but didn't complete the task, check slurm log files"
		exit
	fi

