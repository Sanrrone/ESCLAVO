function statusb {
	set -e

	echo "ESCLAVO: status before begin"

	FORCE=$1
	
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
		echo "ESCLAVO: multiqc_report.html exist (use --force to run it anyway)"
		exit
	else
		rm -f *.html *${PATTERN}.zip rbqc.conf
		rm -rf multiqc_data*

		echo -e "inputFiles\tsize\tstepStatus" > tmp0
		ls -lh *${PATTERN} |awk '{print $NF"\t"$5"\trunning"}' >> tmp1
		echo "timeElapsed" > tmp2
		nfiles=$(ls -1 $FASTQFOLDER/*${PATTERN} |wc -l |awk '{print $1}')
		echo $nfiles |awk '{for(i=1;i<=$1;i++)print "0:0:0"}' >> tmp2
		
		cat tmp0 tmp1 >> tmp && rm -f tmp0 tmp1
		paste tmp tmp2 > rbqc.conf && rm -f tmp2

		echo "timeElapsed" > tmp2
		for fastqfile in $(ls *$PATTERN)
		do
			SECONDS=0
			fastqc -f fastq $fastqfile -o . -t $(nproc)
			duration=$SECONDS
			echo "$(($duration/60/60)):$(($duration/60)):$(($duration % 60))" >> tmp2
		done
		rm -f rbqc.conf
		echo -e "inputFiles\tsize\tstepStatus" > tmp0
		ls -lh *${PATTERN} |awk '{print $NF"\t"$5"\trunning"}' >> tmp1
		cat tmp0 tmp1 > rbqc.conf && rm -f tmp0 tmp1
		ls -1 *.html | grep -v "multiqc_report.html" | awk 'BEGIN{print "qcFiles"}{print $1}' >> tmp3
		paste rbqc.conf tmp2 tmp3 > tmp && rm -f tmp2 rbqc.conf tmp1 tmp3 && mv tmp rbqc.conf
		multiqc .

		sed -i "s/running/done/g" rbqc.conf
		if [ "$PCONF" != "" ]; then
			echo "ESCLAVO: Updating config file: $PCONF"
			sed -i "s/pPercent.*/pPercent\t25/g" $PCONF
			sed -i "s/status.*/status\topen/g" $PCONF
			sed -i "s/lastStep.*/lastStep\tRead status before QC/g" $PCONF
		fi
	fi

	echo "ESCLAVO: status before end"


}