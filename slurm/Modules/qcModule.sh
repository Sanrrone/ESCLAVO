#!/bin/bash 
#SBATCH -J qcModule
#SBATCH -o qc_%a.out
#SBATCH -e qc_%a.err
#SBATCH -t QCTIME
#SBATCH --partition=QCPARTITION
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=USERMAIL
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
pair=$(echo "$i"| awk '{if(NF>1){print "PAIR"}else{print "SINGLE"}}')
R1=$(echo "$i" |awk '{print $1}')
R2=$(echo "$i" |awk '{print $2}')


total=$(wc -l $R1 |awk '{print int($1/4)}')
lengthSeq=$(awk -v total=$total '{if(NR%4==2 && NR<(total/4)) print length($1)}' $R1 | sort -n | uniq -c |tail -n1 |awk '{print $2}')
cutoff=$(echo $lengthSeq |awk '{readl=int($1*0.05); if(readl<=10){print "5"}else{print readl}}')

if [ $cutoff -eq 5 ];then
	tleft=5
else
	tleft=10
fi
name=$(echo "$R1" |rev |cut -d "/" -f 1 |rev |sed "s/.fastq//g")

if [ "$pair" == "PAIR" ];then
	jobid=$(echo "$i" |awk '{print $3}')
	PERLBIN PRINSEQHOME/prinseq-lite.pl -fastq $R1 -fastq2 $R2 -out_good $name.pass -out_bad null -trim_left $tleft -trim_right $cutoff -trim_qual_right 20 -min_len 75 -trim_qual_window 15 -trim_qual_step 5 -ns_max_n 0
else
	jobid=$(echo "$i" |awk '{print $2}')
	PERLBIN PRINSEQHOME/prinseq-lite.pl -fastq $R1 -out_good $name.pass -out_bad null -trim_left $tleft -trim_right $cutoff -trim_qual_right 20 -min_len 50 -trim_qual_window 15 -trim_qual_step 5 -ns_max_n 0
fi

echo "$jobid" >> QC_JOBS.txt

if [ $((jobid)) -eq GNUMS ];then

	completedjobs=$(wc -l QC_JOBS.txt |awk '{print $1}')
	if [ $((completedjobs)) -eq GNUMS ];then
		echo "QC ended, calling chainReaction Module"
		chainReaction "qcModule" $1 $DEBUG
	else
		echo "QC is still working, calling CALLSLURM Module"
		sbatch CALLSLURM "qcModule" "QC_JOBS.txt" GNUMS $1  $DEBUG
	fi
fi

