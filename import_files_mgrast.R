#Import OTU table
otufull=read.csv("C:\\Users\\sierr\\Documents\\Thesis project\\MG-Rast biom files\\otu_mgrast_thirty.csv",header=TRUE)
View(otufull)

#Import metadata
metadata=(read.csv("C:\\Users\\sierr\\Documents\\Thesis project\\MG-Rast biom files\\HPMMMeta_r.csv",header=TRUE))
View(metadata)

#Import taxa
taxa=as.matrix(read.csv("C:\\Users\\sierr\\Documents\\Thesis project\\MG-Rast biom files\\taxa_thirty.csv"))#For other column add Other to end of name
View(taxa)

#Formating OTU table
OTU=otu_table(otufull, taxa_are_rows=TRUE)
row.names(OTU)
(OTU)

TAX=tax_table(taxa)
taxa_names(TAX)
TAX

sampdat=sample_data(metadata)
row.names(sampdat)
(sampdat)
sample_names(sampdat)=metadata$SampleID

taxa_names(TAX)=row.names(OTU)

physeq=phyloseq(OTU,TAX, sampdat)
sample_variables(physeq)
physeq
rank_names(physeq)
colnames(tax_table(physeq)) <- c("Domain",
                                 "Phylum",
                                 "Class",
                                 "Order",
                                 "Family",
                                 "Genus"
                                 )

#go find row names of unclassified taxa and remove them 
View(taxa)

pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}


badTaxa <- c("sp273", "sp272", "sp251", "sp157", "sp209", "sp295")
View(badTaxa)
physeq <- pop_taxa(physeq, badTaxa)
physeq


