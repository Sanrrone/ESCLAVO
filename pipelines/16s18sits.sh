#!/usr/bin/env bash
shopt -s direxpand
shopt -s expand_aliases
if [ -f ~/.bash_profile ]; then source ~/.bash_profile; fi
if [ -f ~/.bashrc ]; then source ~/.bashrc; fi
if [ -f ~/.bash_alias ]; then source ~/.bash_alias; fi

set -e


#usage: bash 16s18sits.sh 0-raw fastq.gz
#software requirements
# FASTQC (better download and create alias), MULTIQC, DADA2

function statusb {
	set -e
	FORCE=$1
	PCONF=$2
	
	if [ -d "$FASTQFOLDER" ];then
		cd $FASTQFOLDER
		nfiles=$(ls -1 *${PATTERN} | wc -l | awk '{print $1}' )
		if [ $((nfiles)) -eq 0 ];then
			echo "ESCLAVO: No files found using the pattern '$PATTERN' in the folder $FASTQFOLDER."
			exit
		fi
	else
		echo "ESCLAVO: $FASTQFOLDER doesn't exist, maybe the path is wrong?. (from: $(pwd))"
		exit
	fi
	if [ $FORCE ];then
		rm -f multiqc_report.html
	fi
	if [ -f multiqc_report.html ];then
		echo "ESCLAVO: multiqc_report.html exist, omitting read status (use --force to run it anyway)"
	else
		rm -f *.html *${PATTERN}.zip rbqc.conf
		rm -rf multiqc_data*

		echo -e "inputFiles\tsize\tstepStatus" > tmp0
		ls -lh *${PATTERN} |awk '{print $NF"\t"$5"\trunning"}' >> tmp1
		cat tmp0 tmp1 >> rbqc.conf && rm -f tmp0 tmp1
		echo "timeElpased" > tmp1
		for fastqfile in $(ls *$PATTERN)
		do
			SECONDS=0
			fastqc -f fastq $fastqfile -o . -t $(nproc)
			duration=$SECONDS
			echo "$(($duration/60/60)):$(($duration/60)):$(($duration % 60))" >> tmp1
		done
		ls -1 *.html | grep -v "multiqc_report.html" | awk 'BEGIN{print "qcFiles"}{print $1}' >> tmp2
		paste rbqc.conf tmp2 tmp1 > tmp && rm -f tmp2 rbqc.conf tmp1 && mv tmp rbqc.conf
		multiqc .

		echo "ESCLAVO: Updating config file: $PCONF"
		sed -i "s/running/done/g" rbqc.conf
		sed -i "s/pPercent.*/pPercent\t25/g" $PCONF
		sed -i "s/lastStep.*/lastStep\tRead status before QC/g" $PCONF
	fi

}

function qc {
	set -e
	echo "ESCLAVO: QC step"
	if [ ! -d 1-qc ];then
		mkdir 1-qc
	fi
	cd 1-qc

	sampleR1=$(ls -1 $FASTQFOLDER/*$PATTERN | head -n1)
	total=$(if [[ "$PATTERN" =~ "gz" ]] || [[ "$PATTERN" =~ "zip" ]]; then zcat $sampleR1 ;else cat $sampleR1 ;fi | wc -l |awk '{print int($1/4)}') #get total fastq sequences
	lengthSeq=$(if [[ "$PATTERN" =~ "gz" ]] || [[ "$PATTERN" =~ "zip" ]]; then zcat $sampleR1 ;else cat $sampleR1 ;fi | awk -v total=$total '{if(NR%4==2 && NR<(total/4)) print length($1)}' | sort -n | uniq -c |tail -n1 |awk '{print $2}')
	nfiles=$(ls -1 $FASTQFOLDER/*${PATTERN} |wc -l |awk '{print $1}')

	echo "timeElpased" > tmp0
	echo $nfiles |awk '{for(i=1;i<=$1;i++)print "0"}' >> tmp0
	echo "inputFiles" > tmp1
	ls -1 $FASTQFOLDER/*${PATTERN} >> tmp1
	echo "stepStatus" > tmp2
	echo $nfiles |awk '{for(i=1;i<=$1;i++)print "running"}' >> tmp2
	paste tmp0 tmp1 tmp2 > qc.conf && rm tmp0 tmp1 tmp2
	
	echo '
	rm(list=ls())
	library("dada2")
	args<-commandArgs()
	path<-args[6]
	fqpattern<-args[7]
	tolerance<-as.numeric(args[8])
	readlength<-as.numeric(args[9])
	projectfolder<-args[10]

	fnFs <- sort(list.files(path, pattern=paste0("1",fqpattern), full.names = TRUE))
	fnRs <- sort(list.files(path, pattern=paste0("2",fqpattern), full.names = TRUE))
	sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
	# Place filtered files in filtered/ subdirectory
	filtFs <- file.path(projectfolder, "1-qc", paste0(sample.names, "_F_filt",fqpattern))
	filtRs <- file.path(projectfolder, "1-qc", paste0(sample.names, "_R_filt",fqpattern))
	names(filtFs) <- sample.names
	names(filtRs) <- sample.names
	readtolerance<-readlength*tolerance
	maxeeformula<- (0.01*readlength)+(0.012589254*readtolerance)
	out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=(readlength-readlength*tolerance),
	              maxN=0, maxEE=maxeeformula, truncQ=2, rm.phix=TRUE, minLen = 80,
	              compress=TRUE, multithread=TRUE)
	write.table(out,"qc_filt.tsv",sep="\t")
	errF <- learnErrors(filtFs, multithread=TRUE)
	errR <- learnErrors(filtRs, multithread=TRUE)

	write.table(errF,"dada2_filt_errF.tsv",sep="\t")
	write.table(errR,"dada2_filt_errR.tsv",sep="\t") ' > dada2_filt.R
	SECONDS=0
	Rscript dada2_filt.R $FASTQFOLDER $PATTERN $TOLERANCE $lengthSeq $PROJECTFOLDER > dada2_filt.log
	duration=$SECONDS
	echo "timeElpased" > tmp0
	echo "$(($duration/60/60)):$(($duration/60)):$(($duration % 60))" |awk -v nfiles=$nfiles '{for(i=1;i<=nfiles;i++){print $1/nfiles}}' >> tmp0
	echo "inputFiles" > tmp1
	ls -1 $FASTQFOLDER/*${PATTERN} >> tmp1
	echo "stepStatus" > tmp2
	echo $nfiles |awk '{for(i=1;i<=$1;i++)print "done"}' >> tmp2
	paste tmp0 tmp1 tmp2 > qc.conf && rm tmp0 tmp1 tmp2
	rm dada2_filt.R

	echo "ESCLAVO: Updating config file: $PCONF"
	sed -i "s/running/done/g" qc.conf
	sed -i "s/pPercent.*/pPercent\t50/g" $PCONF
	sed -i "s/lastStep.*/lastStep\tQC/g" $PCONF
}

function statusa {

	set -e
	echo "ESCLAVO: Read status after QC step"
	FORCE=$1
	PCONF=$2
	
	if [ -d 1-qc ];then
		cd 1-qc
		nfiles=$(ls -1 *${PATTERN} | wc -l | awk '{print $1}' )
		if [ $((nfiles)) -eq 0 ];then
			echo "ESCLAVO: No filtered reads found using the pattern '$PATTERN' in the folder 1-qc."
			exit
		fi
	else
		echo "ESCLAVO: 1-qc doesn't exist, maybe the pipeline bypass the QC step? (from: $(pwd))"
		exit
	fi
	if [ $FORCE ];then
		rm -f multiqc_report.html
	fi
	if [ -f multiqc_report.html ];then
		echo "ESCLAVO: multiqc_report.html exist in QC step, omitting read status (use --force to run it anyway)"
	else
		rm -f *.html *${PATTERN}.zip raqc.conf
		rm -rf multiqc_data*

		echo -e "inputFiles\tsize\tstepStatus" > tmp0
		ls -lh *${PATTERN} |awk '{print $NF"\t"$5"\trunning"}' >> tmp1
		cat tmp0 tmp1 >> raqc.conf && rm -f tmp0 tmp1
		echo "timeElpased" > tmp1
		for fastqfile in $(ls *${PATTERN})
		do
			SECONDS=0
			fastqc -f fastq $fastqfile -o . -t $(nproc)
			duration=$SECONDS
			echo "$(($duration/60/60)):$(($duration/60)):$(($duration % 60))" >> tmp1
		done
		ls -1 *.html | grep -v "multiqc_report.html" | awk 'BEGIN{print "qcFiles"}{print $1}' >> tmp2
		paste raqc.conf tmp2 tmp1 > tmp && rm -f tmp2 raqc.conf tmp1 && mv tmp raqc.conf
		multiqc .

		echo "ESCLAVO: Updating config file: $PCONF"
		sed -i "s/running/done/g" raqc.conf
		sed -i "s/pPercent.*/pPercent\t75/g" $PCONF
		sed -i "s/lastStep.*/lastStep\tRead status after QC/g" $PCONF
	fi

}

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"
case $key in
    -p|--projectfolder)
    PROJECTFOLDER="$2"
    shift # past argument
    shift # past value
    ;;
    -f|--fastqfolder)
    FASTQFOLDER="$2"
    shift # past argument
    shift # past value
    ;;
    -pt|--fastqpattern)
    PATTERN="$2"
    shift # past argument
    shift # past value
    ;;
    -to|--tolerance)
    TOLERANCE="$2"
    shift # past argument
    shift # past value
    ;;
    --force)
    FORCE=true
    shift # past argument
    ;;
    --debug)
    DEBUG=true
    set -ex
    shift # past argument
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
if [ "$TOLERANCE" == "" ];then
	TOLERANCE=0.1
fi
set -- "${POSITIONAL[@]}" # restore positional parameters

if [[ -n $1 ]]; then
    echo "the following argument lack of parameter: $1"
    exit
fi

if [ "$PROJECTFOLDER" == "" ];then
	PROJECTFOLDER="."
fi

cd $PROJECTFOLDER
PNAME=$(pwd | awk -F"/" '{print $NF"_eConf.tsv"}')
PCONF=$(pwd | awk -v pname=$PNAME '{print $0"/"pname}')

if [ ! -f "$PNAME" ];then
	echo "	value" > $PNAME
	echo "pfolder $(pwd)" >> $PNAME
	echo "ffolder $FASTQFOLDER" >> $PNAME
	echo "fqpattern $PATTERN" >> $PNAME
	echo "analysis 16s18sits" >> $PNAME
	echo "aversion 1.0" >> $PNAME
	echo "created $(date --iso-8601)" >> $PNAME
	echo "status open" >> $PNAME
	echo "pPercent 0" >> $PNAME
	echo "lastStep statusb" >> $PNAME
fi


statusb $FORCE $PCONF
cd $PROJECTFOLDER
qc
cd $PROJECTFOLDER
statusa $FORCE $PCONF