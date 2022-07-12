############# MACHINE LEARNING WITH CARET #################

# This workflow was performed in R 4.2.0

# Load packages
library(tidyverse)
library(dplyr)
library(plyr)
library(magrittr)
library(caret)
library(MLeval)
library(caTools)
library(doParallel)
library(plotROC)

# Load data
df <- read.csv("/Users/mimiferreiro/Desktop/Lab/2_urothelial/4_results/2_optimal_cutpoint/censored/cens_regression_data_coxdata_sd.csv")
sig_genes <- read.csv("/Users/mimiferreiro/Desktop/Lab/2_urothelial/4_results/1_filtering/1_data/censored/cens_significative_genes_cox_optimalcut_bonferroni_res_sd.csv")[,-1]

#### 1. PREPARE DATA ####
# Transform response variable to DC and NR
# Set response levels
df$response <- factor(df$response, levels = c("CR", "PR", "SD", "PD", "NE"))

# Selecting significative genes and variables
genelist <- sig_genes$Variable

# Creating new variables and make it factor type
df <- df %>%
  mutate(dc = case_when(
    response == "CR" ~ "DC",
    response == "PR" ~ "DC",
    response == "SD" ~ "DC",
    TRUE ~ "NDC"
  ))

df <- df %>%
  mutate(response = case_when(
    response == "CR" ~ "R",
    response == "PR" ~ "R",
    TRUE ~ "NR"
  ))

df$dc %<>% factor()
df$response %<>% factor()

# Check factor levels
levels(df$dc)
levels(df$response)

# Reversing factor levels
df$dc <- fct_rev(df$dc)
levels(df$dc)

# Reorder df
df <- df %>% dplyr::relocate(response, dc, .after = 3)
df <- df[,5:ncol(df)]

# Reorder genes levels
all_genes <- names(df[,2:ncol(df)])
df[,all_genes] <- as.data.frame(unclass(df[,all_genes]),stringsAsFactors=TRUE)
df[all_genes] <- lapply(df[all_genes],relevel, ref="low") 

#### 2. PREPROCESS DATA ####
# Near zero variables
nzv <- nearZeroVar(df)
head(nzv) # there are not

# Split data (DCR)
index_dcr <- caret::createDataPartition(df$dc, p=0.8, list = F)
train_dcr <- df[index_dcr,]
test_dcr <- df[-index_dcr,]

train_nodcr <- train_dcr[,-1]
test_nodcr <- test_dcr[,-1]

#### 3. TRAINING RANDOM FOREST MODEL ####
# Paralellize process
cl <- makePSOCKcluster(3)
registerDoParallel(cl)

# Define hyperparameters
cv.k  <- 10
reps <- 1

hyperparameters <- expand.grid(mtry = c(2:6), # √15 = 3.8
                               min.node.size = c(1:3),
                               splitrule = "gini")

# Set seeds
set.seed(1234)
seeds <- vector(mode = "list", length = (cv.k * reps) + 1)
for (i in 1:(cv.k * reps)) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[(cv.k * reps) + 1]] <- sample.int(1000, 1)

# Crossvalidation train control
train.control <- trainControl(method = "repeatedcv", number = cv.k,
                              repeats = reps, seeds = seeds,
                              returnResamp = "final", verboseIter = TRUE,
                              allowParallel = TRUE, classProbs = TRUE, 
                              summaryFunction = twoClassSummary, savePredictions = TRUE)

# Train Random Forest Model
set.seed(5678)
rf.model <- train(x = train_nodcr[,genelist],
                  y = train_dcr$dc,
                  method = "ranger",
                  tuneGrid = hyperparameters,
                  metric = "ROC",
                  importance= 'impurity',
                  trControl = train.control,
                  num.trees = 1000)

rf.model

# RF final model
rf.model$finalModel

#### 4. PLOTS ####
# Mtry (hyperparameter)
ggplot(rf.model, highlight = TRUE) +
  scale_x_continuous(breaks = 1:30) +
  labs(title = "RF", x = "mtry") +
  guides(color = guide_legend(title = "min.node.size"),
         shape = guide_legend(title = "min.node.size")) +
  theme_bw()

# Density plot
densityplot(rf.model, metric = "ROC")

# Gene importance
importance.rf <- varImp(rf.model,scale = TRUE)
plot(importance.rf, main = 'Gene importance in RF model')

# AUC-ROC
# Select a parameter setting
selectedIndices <- rf.model$pred$mtry == 6 & rf.model$pred$min.node.size == 1

g <- ggplot(rf.model$pred[selectedIndices, ], aes(m=NDC, d=factor(obs, levels = c("NDC", "DC")))) + 
  geom_roc(n.cuts=0) + 
  coord_equal() +
  style_roc()

g + annotate("text", x=0.75, y=0.25, label=paste("AUC =", round((calc_auc(g))$AUC, 4)))

#### 5. PREDICTION ####
# Training metrics and ROC
train_pred <- predict(rf.model, newdata=train_nodcr[,genelist], type = "prob")
train_metrics <- evalm(data.frame(train_pred, train_dcr$dc, Group = "RF (train)"))
train_metrics

# Test metrics and ROC
test_pred <- predict(rf.model, newdata=test_nodcr[,genelist], type = "prob")
test_metrics <- evalm(data.frame(test_pred, test_dcr$dc, Group = "RF (test)"))
test_metrics

# Calculate test/train accuracies
results <- list(RF=rf.model)
prediction <- extractPrediction(results, 
                                testX = test_nodcr[,genelist], 
                                testY = test_dcr$dc)
pred_test <- prediction %>% filter(dataType == "Test")
pred_train <- prediction %>% filter(dataType == "Training")

accuracy_test <- pred_test %>% mutate(success = ifelse(obs == pred, TRUE, FALSE)) %>% summarise(accuracy = mean(success))
accuracy_train <- pred_train %>% mutate(success = ifelse(obs == pred, TRUE, FALSE)) %>% summarise(accuracy = mean(success))

#### 6. TEST IF OUR GENES ARE INFORMATIVES ####
# Paralellize process
cl <- makePSOCKcluster(3)
registerDoParallel(cl)
reps <- 1

# Train Random Forest Model
set.seed(5678)
rf.model <- train(x = train_nodcr,
                  y = train_dcr$dc,
                  method = "ranger",
                  tuneGrid = hyperparameters,
                  metric = "Accuracy",
                  importance= 'impurity',
                  trControl = train.control,
                  num.trees = 1000)

rf.model
rf.model$finalModel

# Training metrics and ROC
train_pred <- predict(rf.model, newdata=train_nodcr, type = "prob")
train_metrics <- evalm(data.frame(train_pred, train_dcr$dc, Group = "RF (train)"))
train_metrics

# Test metrics and ROC
test_pred <- predict(rf.model, newdata=test_nodcr, type = "prob")
test_metrics <- evalm(data.frame(test_pred, test_dcr$dc, Group = "RF (test)"))
test_metrics

# Calculate test/train accuracies
results <- list(RF=rf.model)
prediction <- extractPrediction(results, 
                                testX = test_nodcr, 
                                testY = test_dcr$dc)
pred_test <- prediction %>% filter(dataType == "Test")
pred_train <- prediction %>% filter(dataType == "Training")

accuracy_test <- pred_test %>% mutate(success = ifelse(obs == pred, TRUE, FALSE)) %>% summarise(accuracy = mean(success))
accuracy_train <- pred_train %>% mutate(success = ifelse(obs == pred, TRUE, FALSE)) %>% summarise(accuracy = mean(success))

# Gene importance
importance.rf <- varImp(rf.model,scale = TRUE)
plot(importance.rf, main = 'Gene importance in RF model')
