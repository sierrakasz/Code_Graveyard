#Import either rectum, mouth, or total
sampletype <- read.csv("C:/Users/sierr/Documents/Thesis project/R Script/sample_number.csv")
View(sampletype)


dunnTest(area~pipeline, data=sampletype, method = 'bonferroni')


levels(sampletype$pipeline)

class(sampletype$area)

group_by(sampletype, pipeline) %>%
  summarise(
    count = n(),
    mean = mean(area, na.rm = TRUE),
    sd = sd(area, na.rm = TRUE)
  )

#boxplot
mycompare <- list(c("mothur", "QIIME"), c("mothur", "MG-RAST"), c("QIIME", "MG-RAST"))
p <- ggbarplot(sampletype, x = "pipeline", y = "area", 
               order = c("mothur", "QIIME", "MG-RAST"),
               fill = "Subsample", palette = c("#00AFBB", "#E7B800", "#FC4E07", "#00CC33"),
               ylab = "Number of Samples Post Filtering", xlab = "Pipeline",
               title= "Sample Size Post Filtering",
               position= position_dodge()) + stat_compare_means(label = 'p.signif', method= 'kruskal.test')#, comparisons = mycompare)
facet(p, facet.by = "Sample_Area")

#Is the data normal? Signifcant = non normal
shapiro.test(sampletype$area)

#Decides hetero or homoscedastic t test
#Normal data, variance, signficant = equal variance
bartlett.test(area ~ pipeline, data = sampletype)

#Non normal data variance, signficant = equal variance
fligner.test(area ~ pipeline, data = sampletype)

#Anova
res.aov <- aov(reads ~ pipeline, data = sampletype)
summary(res.aov)

#comparing means between groups 
TukeyHSD(res.aov)

#Kruskal test, non para
kruskal.test(sampletype$area, sampletype$pipeline)

#Pairwise, non para
dunn.test(sampletype$area, sampletype$pipeline)

#Pairwise, non para bigger sample
pairwise.wilcox.test(sampletype$area, sampletype$pipeline, p.adjust.method = "BH")



