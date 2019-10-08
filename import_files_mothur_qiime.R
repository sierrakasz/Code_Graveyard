biom <- import_biom("C:/Users/sierr/Documents/Thesis project/mothur biom files/mothur.sixty.biom")
(biom)

metadata=(read.csv("C:/Users/sierr/Documents/Thesis project/Metadata/HPMMMeta_r.csv",header=TRUE))
View(metadata)


tree=read_tree("C:/Users/sierr/Documents/Thesis project/R Script/tree.nwk", errorIfNULL=TRUE)
(tree)

sampdat=sample_data(metadata)
row.names(sampdat)
(sampdat)
sample_names(sampdat)=metadata$SampleID


physeq=merge_phyloseq(biom, sampdat)
rank_names(physeq)
#Qiime
colnames(tax_table(physeq)) <- c("Domain",
                             "Phylum",
                             "Class",
                             "Order",
                             "Family",
                             "Genus",
                             "Species")
#mothur
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus")
(physeq)
taxa_names(physeq)

pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}
tax_table(physeq)
badTaxa <- c("Otu000056")
View(badTaxa)
physeq <- pop_taxa(PhylumLevel, badTaxa)
physeq

tax_table(physeq)
options(max.print=1000000)
