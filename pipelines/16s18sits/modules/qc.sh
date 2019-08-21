function qc {
	set -e
	echo "ESCLAVO: QC begin"
	if [ ! -d 1-qc ];then
		mkdir 1-qc
	fi
	
	cd 1-qc

	if [ $FORCE ];then
		rm -rf *
	fi

	sampleR1=$(ls -1 $FASTQFOLDER/*$PATTERN | head -n1)
	total=$(if [[ "$PATTERN" =~ "gz" ]] || [[ "$PATTERN" =~ "zip" ]]; then zcat $sampleR1 ;else cat $sampleR1 ;fi | wc -l |awk '{print int($1/4)}') #get total fastq sequences
	lengthSeq=$(if [[ "$PATTERN" =~ "gz" ]] || [[ "$PATTERN" =~ "zip" ]]; then zcat $sampleR1 ;else cat $sampleR1 ;fi | awk -v total=$total '{if(NR%4==2 && NR<(total/4)) print length($1)}' | sort -n | uniq -c |tail -n1 |awk '{print $2}')
	nfiles=$(ls -1 $FASTQFOLDER/*${PATTERN} |wc -l |awk '{print $1}')

	echo "timeElapsed" > tmp0
	echo $nfiles |awk '{for(i=1;i<=$1;i++)print "0:0:0"}' >> tmp0
	echo "inputFiles" > tmp1
	ls -1 $FASTQFOLDER/*${PATTERN} >> tmp1
	echo "stepStatus" > tmp2
	echo $nfiles |awk '{for(i=1;i<=$1;i++)print "running"}' >> tmp2
	paste tmp0 tmp1 tmp2 > qc.conf && rm tmp0 tmp1 tmp2
	
	echo "
	rm(list=ls())
	if('BiocManager' %in% rownames(installed.packages()) == FALSE) {install.packages('BiocManager')}
	if('dada2' %in% rownames(installed.packages()) == FALSE) {BiocManager::install('dada2')}
	library(dada2)
	options(scipen=999)
	
	path<-'$FASTQFOLDER'
	fqpattern<-'$PATTERN'
	tolerance<-$TOLERANCE
	readlength<-$lengthSeq
	projectfolder<-'$PROJECTFOLDER'

	fnFs <- sort(list.files(path, pattern=paste0('1',fqpattern), full.names = TRUE))
	fnRs <- sort(list.files(path, pattern=paste0('2',fqpattern), full.names = TRUE))
	sample.names <- sapply(strsplit(basename(fnFs), '_'), \`[\`, 1)
	# Place filtered files in filtered/ subdirectory
	filtFs <- file.path(projectfolder, '1-qc', paste0(sample.names, '_F_filt',fqpattern))
	filtRs <- file.path(projectfolder, '1-qc', paste0(sample.names, '_R_filt',fqpattern))
	names(filtFs) <- sample.names
	names(filtRs) <- sample.names
	readtolerance<-readlength*tolerance
	maxeeformula<- (0.01*readlength)+(0.012589254*readtolerance)
	print(paste0('Doing filtering at maxEE ',maxeeformula))
	out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=(readlength-readlength*tolerance),
	              maxN=0, maxEE=maxeeformula, truncQ=2, rm.phix=TRUE, minLen = 80,
	              compress=TRUE, multithread=TRUE)
	print('Done')
	write.table(out,'qc_filt.tsv',sep='\t')
	print('learning from errors')
	errF <- learnErrors(filtFs, multithread=TRUE)
	errR <- learnErrors(filtRs, multithread=TRUE)
	print('Done')
	write.table(errF,'dada2_filt_errF.tsv',sep='\t')
	write.table(errR,'dada2_filt_errR.tsv',sep='\t') 

	print('Dereplication')
	derepFs <- derepFastq(filtFs, verbose=TRUE)
	derepRs <- derepFastq(filtRs, verbose=TRUE)
	print('Done')

	print('Inferring sample composition')
	dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
	dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
	print('done')
	
	print('Finally merge reads')
	mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)	
	seqtab <- makeSequenceTable(mergers)
	print('Done')

	#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(250,256)]) # to select specific seq length
	seqtab.nochim <- removeBimeraDenovo(seqtab, method='consensus', multithread=TRUE, verbose=TRUE)
	print(paste0('sequences kept after removing chimera step: ',sum(seqtab.nochim)/sum(seqtab)))

	#make summary table
	getN <- function(x) sum(getUniques(x))
	if(length(fnFs)==1){
	  dadaFs<-list(dadaFs)
	  dadaRs<-list(dadaRs)
	  mergers<-list(mergers)
	}
	track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
	# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
	colnames(track) <- c('input', 'filtered', 'denoisedF', 'denoisedR', 'merged', 'nonchim')
	rownames(track) <- sample.names

	write.table(track,'qc_summary.tsv', sep='\t')
	save(seqtab.nochim,file = 'seqtab.nochim.RData')

	" > dada2_filt.R
	SECONDS=0
	Rscript --vanilla dada2_filt.R > dada2_filt.log
	duration=$SECONDS
	echo "timeElapsed" > tmp0
	echo $duration | awk -v nfiles=$nfiles -v duration=$duration '{for(i=1;i<=nfiles;i++){print int($1/60/60/nfiles)":"int($1/60/nfiles)":"($1%60)/nfiles}}' >> tmp0
	echo "inputFiles" > tmp1
	ls -1 $FASTQFOLDER/*${PATTERN} >> tmp1
	echo "stepStatus" > tmp2
	echo $nfiles |awk '{for(i=1;i<=$1;i++)print "done"}' >> tmp2
	paste tmp0 tmp1 tmp2 > qc.conf && rm tmp0 tmp1 tmp2
	rm dada2_filt.R

	if [ "$PCONF" != "" ]; then
		echo "ESCLAVO: Updating config file: $PCONF"
		sed -i "s/running/done/g" qc.conf
		sed -i "s/pPercent.*/pPercent\t50/g" $PCONF
		sed -i "s/lastStep.*/lastStep\tQC/g" $PCONF
	fi


	echo "ESCLAVO: QC end"
}
