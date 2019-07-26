
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
	write.table(errR,"dada2_filt_errR.tsv",sep="\t") 
