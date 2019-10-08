#change name depending on rank names
rank_names(physeq)
GPrPhylum=tax_glom(physeq, "Phylum")
GPrPhylum
PhylumLevel = filter_taxa(GPrPhylum, function(x) mean(x) > 0.001, TRUE) 
PhylumLevel

plot_tree(PhylumLevel, color = "Sample_Area", label.tips = "Phylum", nodelabf = nodeplotblank, plot.margin = 0.5, ladderize = TRUE)

GPrFamily=tax_glom(physeq,"Rank6")
FamilyLevel = filter_taxa(GPrFamily, function(x) mean(x) > 0.001, TRUE)
FamilyLevel

plot_tree(FamilyLevel, color = "Sample_Area", label.tips = "Rank6", nodelabf = nodeplotblank, plot.margin = 0.5, ladderize = TRUE)
