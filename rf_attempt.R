library(AmesHousing)
library(randomForest)
library(rsample)
library(ggplot2)
library(tidyr)

ames_split <- initial_split(AmesHousing::make_ames(), prop = .7)
ames_train <- training(ames_split)
ames_test  <- testing(ames_split)

set.seed(123)

# default RF model
m1 <- randomForest(
  formula = Sale_Price ~ .,
  data    = ames_train
)

m1

plot(m1)

# number of trees with lowest MSE
which.min(m1$mse)

# RMSE of this optimal random forest
sqrt(m1$mse[which.min(m1$mse)])

#test set training set
valid_split <- initial_split(ames_train, .8)

# training data
ames_train_v2 <- analysis(valid_split)

# validation data
ames_valid <- assessment(valid_split)
x_test <- ames_valid[setdiff(names(ames_valid), "Sale_Price")]
y_test <- ames_valid$Sale_Price

rf_oob_comp <- randomForest(
  formula = Sale_Price ~ .,
  data    = ames_train_v2,
  xtest   = x_test,
  ytest   = y_test
)

# extract OOB & validation errors
oob <- sqrt(rf_oob_comp$mse)
validation <- sqrt(rf_oob_comp$test$mse)

# compare error rates
tibble::tibble(
  `Out of Bag Error` = oob,
  `Test error` = validation,
  ntrees = 1:rf_oob_comp$ntree
) %>%
  gather(Metric, RMSE, -ntrees) %>%
  ggplot(aes(ntrees, RMSE, color = Metric)) +
  geom_line() +
  scale_y_continuous(labels = scales::dollar) +
  xlab("Number of trees")

#tune the tree
features <- setdiff(names(ames_train), "Sale_Price")

m2 <- tuneRF(
  x          = ames_train[features],
  y          = ames_train$Sale_Price,
  ntreeTry   = 500,
  mtryStart  = 5,
  stepFactor = 1.5,
  improve    = 0.01,
  trace      = FALSE      # to not show real-time progress 
)

#trying it with microbiome data
otu <- as.data.frame(t(otu_table(physeq_q_f)))
otu$SampleID <- rownames(otu)
meta_sa <- metadata %>% select(SampleID, MoD)
otu <- merge(meta_sa, otu, by = 'SampleID')
otu <- otu[,-1]

otu_split <- initial_split(otu, prop = .7)
otu_train <- training(otu_split)
otu_test  <- testing(otu_split)




m1 <- randomForest(
  formula = MoD ~ .,
  data    = otu,
  ntree= 1000, 
  mtry = 5
)

m1

plot(m1)


valid_split <- initial_split(otu, .7)
otu_train <- analysis(valid_split)
otu_valid <- assessment(valid_split)
x_test <- otu_valid[setdiff(names(otu_valid), "MoD")]
y_test <- otu_valid$MoD

rf_oob_comp <- randomForest(
  formula = MoD ~ .,
  data    = otu_train,
  xtest   = x_test,
  ytest   = y_test,
  ntree =1000
)

m1
rf_oob_comp

oob <- sqrt(rf_oob_comp$err.rate)
validation <- sqrt(rf_oob_comp$test$err.rate)

# compare error rates
tibble::tibble(
  `Out of Bag Error` = oob,
  `Test error` = validation,
  ntrees = 1:rf_oob_comp$ntree
) %>%
  gather(Metric, Error, -ntrees) %>%
  ggplot(aes(ntrees, Error, color = Metric)) +
  geom_line() +
  xlab("Number of trees")

#tune the tree
features <- setdiff(names(otu), "Sample_Area")

m2 <- tuneRF(
  x          = otu[features],
  y          = otu$Sample_Area,
  ntreeTry   = 1000,
  mtryStart  = 15,
  stepFactor = 1.5,
  improve    = 0.01,
  trace      = FALSE      # to not show real-time progress 
)


