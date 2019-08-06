function statusa {

	set -e

	echo "ESCLAVO: Read status after QC begin"
	FORCE=$1
	
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

		echo "timeElapsed" > tmp2
		nfiles=$(ls -1 $FASTQFOLDER/*${PATTERN} |wc -l |awk '{print $1}')
		echo $nfiles |awk '{for(i=1;i<=$1;i++)print "0:0:0"}' >> tmp2
		paste raqc.conf tmp2 > tmp && rm raqc.conf tmp2 && mv tmp raqc.conf

		echo "timeElapsed" > tmp2
		for fastqfile in $(ls *${PATTERN})
		do
			SECONDS=0
			fastqc -f fastq $fastqfile -o . -t $(nproc)
			duration=$SECONDS
			echo "$(($duration/60/60)):$(($duration/60)):$(($duration % 60))" >> tmp2
		done
		rm raqc.conf
		echo -e "inputFiles\tsize\tstepStatus" > tmp0
		ls -lh *${PATTERN} |awk '{print $NF"\t"$5"\trunning"}' >> tmp1
		cat tmp0 tmp1 > raqc.conf && rm -f tmp0 tmp1
		ls -1 *.html | grep -v "multiqc_report.html" | awk 'BEGIN{print "qcFiles"}{print $1}' >> tmp3
		paste raqc.conf tmp2 tmp3 > tmp && rm -f tmp2 raqc.conf tmp1 tmp3 && mv tmp raqc.conf
		multiqc .
		
		sed -i "s/running/done/g" raqc.conf

		if [ "$PCONF" != "" ]; then
			echo "ESCLAVO: Updating config file: $PCONF"
			sed -i "s/pPercent.*/pPercent\t75/g" $PCONF
			sed -i "s/lastStep.*/lastStep\tRead status after QC/g" $PCONF
		fi
	fi

	echo "ESCLAVO: status after QC end"

}