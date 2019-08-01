
	library(dada2)
	library(DECIPHER)
	library(phyloseq)
	library(Biostrings)
	library(ggplot2)
	
	load('/home/sandro/Programas/ESCLAVO/projects/1-qc/seqtab.nochim.RData')
	dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
	print('Doing assigment')
	taxa <- assignTaxonomy(seqtab.nochim, '/home/sandro/Descargas/silva_nr_v132_train_set.fa', multithread=TRUE, verbose = T)
	taxa <- addSpecies(taxa, '/home/sandro/Descargas/silva_species_assignment_v132.fa')
	print('Done')
	taxa.print <- taxa # Removing sequence rownames for display only
	rownames(taxa.print) <- NULL
	head(taxa.print)

	ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(taxa))
	ps <- prune_samples(sample_names(ps) != 'Mock', ps) # Remove mock sample

	top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
	ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
	ps.top20 <- prune_taxa(top20, ps.top20)

	pdf('sampleTaxComposition.pdf', width=10, height=10)
	plot_bar(ps.top20, x='Sample', fill='Family') + theme_minimal()
	dev.off()

	
