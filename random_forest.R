
ForestData=physeq
predictors=t(otu_table(ForestData))
dim(predictors)
response <- as.factor(sample_data(ForestData)$Sample_Area)
rf.data <- data.frame(response, predictors)
MozzieForest <- randomForest(response~., data = rf.data, ntree = 1000)
print(MozzieForest)#returns overall Random Forest results



imp <- importance(MozzieForest)#all the steps that are imp or imp. are building a dataframe that contains info about the taxa used by the Random Forest test to classify treatment 
imp <- data.frame(predictors = rownames(imp), imp)
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
imp.20 <- imp.sort[1:20, ]

ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +
  ggtitle("Important OTUs for Classifying Samples by Area")
imp.20$MeanDecreaseGini
otunames <- imp.20$predictors
r <- rownames(tax_table(ForestData)) %in% otunames
kable(tax_table(ForestData)[r, ])#returns a list of the most important predictors for Random Forest Classification
