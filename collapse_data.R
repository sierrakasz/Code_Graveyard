setwd("C:/Users/sierr/Documents/Thesis project")

library(tidyverse)
library(tidyr)
library(tidyselect)

x <- data.frame(otu_table(physeq_mothur_60))
y <- data.frame(tax_table(physeq_mothur_60))

x <- read.csv("otu_all.csv")
y <- read.csv("tax_all.csv")

z <- merge(x,y)
View(z)
z <- z[-1]

tax_group <- z %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))

View(tax_group) 

tax_all <- tax_group[,1:6]
View(tax_all)

otu_all <- tax_group[,-1:-6]
View(otu_all)

write.csv(tax_all, "tax_all.csv")
write.csv(otu_all, "otu_all.csv")

otufull=read.csv("otu_all.csv",header=TRUE)
otufull <- otufull [,-1]
View(otufull)

taxa=as.matrix(read.csv("tax_all.csv"))
taxa <- taxa[,-1]
View(taxa)

OTU=otu_table(otufull, taxa_are_rows=TRUE)
row.names(OTU)

TAX=tax_table(taxa)
taxa_names(TAX)

taxa_names(TAX)=row.names(OTU)

physeq_all=phyloseq(OTU,TAX)

metadata=(read.csv("C:/Users/sierr/Documents/Thesis project/Metadata/HPMMMeta_r_merge.csv",header=TRUE))
View(metadata)

sampdat=sample_data(metadata)
row.names(sampdat)
(sampdat)
sample_names(sampdat)=metadata$SampleID


merge=merge_phyloseq(physeq_all, sampdat)
