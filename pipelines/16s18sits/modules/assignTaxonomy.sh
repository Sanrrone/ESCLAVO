function assignTaxonomy {
	##### to make new silva db for dada2 (instead of decipher)
    #path <- "~/Desktop/Silva/Silva.nr_v132"
    #dada2:::makeTaxonomyFasta_Silva(file.path(path, "silva.nr_v132.align"), file.path(path, "silva.nr_v132.tax"), "~/tax/silva_nr_v132_train_set.fa.gz")
    #dada2:::makeSpeciesFasta_Silva("~/Desktop/Silva/SILVA_132_SSURef_tax_silva.fasta.gz", "~/tax/silva_species_assignment_v132.fa.gz")
	
	if [ "$1" == "" ];then
		WORKFOLDER="2-taxInsight"
	else
    	WORKFOLDER=$1
    fi

    if [ ! -d "$WORKFOLDER" ]; then
    	mkdir "$WORKFOLDER"
    fi

    cd 2-taxInsight

    nfiles=$(ls -1 $FASTQFOLDER/*${PATTERN} |wc -l |awk '{print $1}')
	echo "timeElpased" > tmp0
	echo $nfiles |awk '{for(i=1;i<=$1;i++)print "0:0:0"}' >> tmp0
	echo "inputFiles" > tmp1
	ls -1 $FASTQFOLDER/*${PATTERN} >> tmp1
	echo "stepStatus" > tmp2
	echo $nfiles |awk '{for(i=1;i<=$1;i++)print "running"}' >> tmp2
	paste tmp0 tmp1 tmp2 > tc.conf && rm tmp0 tmp1 tmp2


	echo "ESCLAVO: assignTaxonomy begin"
	echo "
	library(dada2)
	library(DECIPHER)
	library(phyloseq)
	library(Biostrings)
	library(ggplot2)
	
	load('$PROJECTFOLDER/1-qc/seqtab.nochim.RData')
	dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
	print('Doing assigment')
	taxa <- assignTaxonomy(seqtab.nochim, '$ESCLAVOHOME/DB/silva_nr_v132_train_set.fa.gz', multithread=TRUE, verbose = T)
	taxa <- addSpecies(taxa, '$ESCLAVOHOME/DB/silva_species_assignment_v132.fa.gz')
	print('Done')
	#taxa.print <- taxa # Removing sequence rownames for display only
	#rownames(taxa.print) <- NULL
	#head(taxa.print)

	abudancedf<-as.data.frame(t(seqtab.nochim))
	abudancedf<-cbind(taxa,abudancedf[rownames(taxa),])
	rownames(abudancedf)<-1:nrow(abudancedf)
	
	write.table(abudancedf,'abundance.tsv',row.names = F,sep = '\t')

	ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(taxa))
	write.table(otu_table(ps),'otu_table.tsv',sep='\t')
	write.table(tax_table(ps),'tax_table.tsv',sep='\t')


	top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
	ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
	ps.top20 <- prune_taxa(top20, ps.top20)
	wformula=4 + length(sample_names(ps))*2.5
	hformula=10
	pdf('sampleTaxComposition.pdf', width=wformula, height=hformula)
	plot_bar(ps.top20, x='Sample', fill='Genus') + geom_bar(position='fill', stat='identity', color='black') + 
	  theme_minimal() + theme(axis.text.x = element_text(angle = 65, hjust = 1))
	plot_bar(ps.top20, x='Sample', fill='Family') + geom_bar(position='fill', stat='identity', color='black') + 
	  theme_minimal() + theme(axis.text.x = element_text(angle = 65, hjust = 1))

	dev.off()

	" > dada2_assign.R

	SECONDS=0
	Rscript --vanilla dada2_assign.R > dada2_assign.log
	duration=$SECONDS
	
	#rm dada2_assign.R
	
	duration=$(echo $duration | awk -v nfiles=$nfiles -v duration=$duration '{print int($1/60/60/nfiles)":"int($1/60/nfiles)":"($1%60)/nfiles}')

	if [ "$PCONF" != "" ]; then
		echo "ESCLAVO: Updating config file: $PCONF"
		sed -i "s/running/done/g" tc.conf
		sed -i "s/0:0:0/$duration/g" tc.conf
		sed -i "s/status.*/status\tdone/g" $PCONF
		sed -i "s/pPercent.*/pPercent\t100/g" $PCONF
		sed -i "s/lastStep.*/lastStep\tTaxonomic counts/g" $PCONF
	fi

	echo "ESCLAVO: assignTaxonomy end"
}