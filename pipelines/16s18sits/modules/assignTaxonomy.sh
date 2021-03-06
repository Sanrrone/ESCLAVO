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
	echo "timeElapsed" > tmp0
	echo $nfiles |awk '{for(i=1;i<=$1;i++)print "0:0:0"}' >> tmp0
	echo "inputFiles" > tmp1
	ls -1 ../1-qc/*${PATTERN} >> tmp1
	echo "stepStatus" > tmp2
	echo $nfiles |awk '{for(i=1;i<=$1;i++)print "running"}' >> tmp2
	paste tmp0 tmp1 tmp2 > tc.conf && rm tmp0 tmp1 tmp2


	echo "ESCLAVO: assignTaxonomy begin"
	echo "
	if('BiocManager' %in% rownames(installed.packages()) == FALSE) {install.packages('BiocManager')}
	if('dada2' %in% rownames(installed.packages()) == FALSE) {BiocManager::install('dada2')}
	if('phyloseq' %in% rownames(installed.packages()) == FALSE) {BiocManager::install('phyloseq')}
	if('Biostrings' %in% rownames(installed.packages()) == FALSE) {BiocManager::install('Biostrings')}
	if('ggplot2' %in% rownames(installed.packages()) == FALSE) {install.packages('ggplot2')}

	library(dada2)
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

	if(file.exists('$FASTQFOLDER/metadata.tsv')){
	  metadata<-read.table('$FASTQFOLDER/metadata.tsv', sep='\t', header=T, stringsAsFactors = F)
	  rownames(metadata)<-metadata\$sample
	  write.table(sample_data(metadata),'sample_table.tsv',sep='\t')
	  sample_data(ps)<-sample_data(metadata)
	  ordu = ordinate(ps, 'PCoA', weighted=TRUE)
	  png('pcoa.png', width = 800,height = 600, units = 'px')
	  print(plot_ordination(ps, ordu, color='treatment') + geom_point(size=4))
	  dev.off()
	}else{
	  metadata<-data.frame()
	}


	top10 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:10]
	ps.top10 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
	ps.top10 <- prune_taxa(top10, ps.top10)
	hformula=10

	if(nrow(metadata)!=0 & 'treatment' %in% colnames(metadata) & 'sample' %in% colnames(metadata)){

	  wformula=4 + length(sample_names(ps))*2.5 + length(unique(metadata$treatment))
	  pdf('sampleTaxComposition.pdf', width=wformula, height=hformula)
	  plot_bar(ps.top10, x='Sample', fill='Genus') + geom_bar(position='fill', stat='identity', color='black') + 
	    theme_minimal() + theme(axis.text.x = element_text(angle = 65, hjust = 1)) + facet_wrap(.~treatment,scales = 'free_x')
	  plot_bar(ps.top10, x='Sample', fill='Family') + geom_bar(position='fill', stat='identity', color='black') + 
	    theme_minimal() + theme(axis.text.x = element_text(angle = 65, hjust = 1)) + facet_wrap(.~treatment,scales = 'free_x')
	}else{
	  wformula=4 + length(sample_names(ps))*2.5
	  
	  pdf('sampleTaxComposition.pdf', width=wformula, height=hformula)
	  plot_bar(ps.top10, x='Sample', fill='Genus') + geom_bar(position='fill', stat='identity', color='black') + 
	    theme_minimal() + theme(axis.text.x = element_text(angle = 65, hjust = 1))
	  plot_bar(ps.top10, x='Sample', fill='Family') + geom_bar(position='fill', stat='identity', color='black') + 
	    theme_minimal() + theme(axis.text.x = element_text(angle = 65, hjust = 1))
	}
	dev.off()

	" > dada2_assign.R

	SECONDS=0
	Rscript --vanilla dada2_assign.R > dada2_assign.log
	duration=$SECONDS
	
	#rm dada2_assign.R
	
	duration=$(echo $duration | awk -v nfiles=$nfiles '{print int($1/60/60/nfiles)":"int($1/60/nfiles)":"($1%60)/nfiles}')

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