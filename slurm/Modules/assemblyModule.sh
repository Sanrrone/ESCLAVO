#!/bin/bash 
#SBATCH -J assembly
#SBATCH -o assembly_%a.out
#SBATCH -e assembly_%a.err
#SBATCH -t ASSEMBLYTIME
#SBATCH --partition=ASSEMBLYPARTITION
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=USERMAIL
#SBATCH -N 1
#SBATCH --array=1-$GNUMS

if [[ "$@" =~ "--debug" ]]; then
		DEBUG="--debug"
        set -ex
else
        set -e
        DEBUG=""
fi

source CHAINREACTION

CORES=$(grep -c ^processor /proc/cpuinfo |awk '{print $1*2}')

i=$(sed -n "$SLURM_ARRAY_TASK_ID"p INPUTREADSLIST)
pair=$(echo "$i" |awk '{if(NF>1){print "PAIR"}else{print "SINGLE"}}')
R1=$(echo "$i" |awk '{print $1}')
R2=$(echo "$i" |awk '{print $2}')
name=$(echo "$R1" |rev |cut -d "/" -f 1 |rev |sed "s/\.fastq//g")

if [ "$pair" == "PAIR" ];then
	SPADESBIN -1 $R1 -2 $R2 -t $CORES -m SPADESMEMORY -o $name.assembly SPADESEXTRACODE
else
	SPADESBIN -s $R1 -t $CORES -m SPADESMEMORY -o $name.assembly SPADESEXTRACODE
fi

cd $name.assembly
echo '#!/usr/bin/perl

use strict;
use warnings;
my $minlen = shift or die "Error: `minlen` parameter not provided\n";
{
    local $/=">";
    while(<>) {
        chomp;
        next unless /\w/;
        s/>$//gs;
        my @chunk = split /\n/;
        my $header = shift @chunk;
        my $seqlen = length join "", @chunk;
        print ">$_" if($seqlen >= $minlen);
    }
    local $/="\n";
}' > removeSmalls.pl
perl removeSmalls.pl 500 scaffolds.fasta > scaffolds_filter500.fasta
rm removeSmalls.pl
cd ..

echo "$jobid" >> ASSEMBLY_JOBS.txt

if [ $((jobid)) -eq GNUMS ];then

	completedjobs=$(wc -l ASSEMBLY_JOBS.txt |awk '{print $1}')
	if [ $((completedjobs)) -eq GNUMS ];then
		echo "ASSEMBLY ended, calling chainReaction Module"
		chainReaction "assembly" $1 $DEBUG
	else
		echo "ASSEMBLY module is still working, calling CALLSLURM Module"
		sbatch CALLSLURM "assembly" "ASSEMBLY_JOBS.txt" GNUMS $1  $DEBUG
	fi
fi


