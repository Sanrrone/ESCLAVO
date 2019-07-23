#!/usr/bin/env bash

if [[ "$@" =~ "--debug" ]]; then
		DEBUG="--debug"
        set -ex
else
        set -e
        DEBUG=""
fi

#usage: bash 16s18sits.sh 0-raw fastq.gz
#software requirements
# FASTQC (better download and create alias), MULTIQC, 

fastqFolder=$1
fastqPattern=$2

function statusb {

	fastqFolder = $1
	fastqPattern = $2
	
	if [ -d $fastqFolder ];then
		cd $fastqFolder
		nfiles=$(ls -1 *${fastqPattern} | wc -l | awk '{print $1}' )
		if [ $((nfiles)) -eq 0 ];then
			echo "* No files found using the pattern $fastqPattern in the folder $fastqFolder."
		fi
	else
		echo "* $fastqFolder doesn't exist, check the path name."
		exit
	fi

	fastqc -f fastq $fastqPattern -o . -t $(nproc)
	multiqc .
}

statusb $fastqFolder