
	library(dada2)
	library(DECIPHER)
	library(phyloseq)
	library(Biostrings)
	library(ggplot2)
	
	load('/home/sandro/Programas/ESCLAVO/projects/1-qc/seqtab.nochim.RData')
	dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
	print('Doing assigment')
	taxa <- assignTaxonomy(seqtab.nochim, '/home/sandro/Programas/ESCLAVO/pipelines/16s18sits/DB/silva_nr_v132_train_set.fa.gz', multithread=TRUE, verbose = T)
	taxa <- addSpecies(taxa, '/home/sandro/Programas/ESCLAVO/pipelines/16s18sits/DB/silva_species_assignment_v132.fa.gz')
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


	top10 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:10]
	ps.top10 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
	ps.top10 <- prune_taxa(top10, ps.top10)
	wformula=4 + length(sample_names(ps))*2.5
	hformula=10
	pdf('sampleTaxComposition.pdf', width=wformula, height=hformula)
	plot_bar(ps.top10, x='Sample', fill='Genus') + geom_bar(position='fill', stat='identity', color='black') + 
	  theme_minimal() + theme(axis.text.x = element_text(angle = 65, hjust = 1))
	plot_bar(ps.top10, x='Sample', fill='Family') + geom_bar(position='fill', stat='identity', color='black') + 
	  theme_minimal() + theme(axis.text.x = element_text(angle = 65, hjust = 1))

	dev.off()

	
