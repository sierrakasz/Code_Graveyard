setwd("C:/Users/sierr/Documents/ThesisProject_Pipelines/")

data <- read.csv("unclassified_reads_heatmap.csv", header = T)

library(ggplot2)

theme_set(theme_classic(base_size = 20))
tiff("uncl_heatmap.TIF", width = 700, height = 300)
p <- ggplot(data, aes(Pipeline, Taxonomy)) + geom_tile(aes(fill = Unclassified_Reads)) +
  scale_fill_gradient(low = "#BE93D6", high = "#4C55B1", name = "Number of Unclassified Reads") + 
  geom_text(aes(label = round(Unclassified_Reads, 1))) +
  scale_y_discrete(labels=c('Phylum_R' = 'Phylum/Rectum', 'Phylum_M' = 'Phylum/Mouth',
                            'Family_R' = 'Family/Rectum', 'Family_M' = 'Family/Mouth')) +
  ylab("Taxonomy/Body Site") 
p
dev.off()
