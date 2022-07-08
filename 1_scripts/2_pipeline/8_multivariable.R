############# MULTIVARIABLE ANALYSIS WITH DIFFERENTS MODELS #################

# This workflow was performed in R 4.2.0

# Load packages
library(tidyverse)
library(survival)
library(glmnet)
library(svglite)
library(Hmisc)
library(corrplot)
library(MASS)
library(ggRandomForests)
library(gtable)
library(grid)

#Change working directory
setwd("/Users/mimiferreiro/Documents/GitHub/tfm_mUC")

# Load functions
source("1_scripts/1_functions/flattenCorrMatrix.R")

# Load data
df <- read.csv("3_results/2_optimal_cutpoint/cens_regression_data_coxdata_sd.csv")
genelist <- read.csv("3_results/1_filtering/1_data/cens_significative_genes_cox_optimalcut_bonferroni_res_sd.csv")[,2]

#### 1. PREPARE DATA ####
df <- df %>% dplyr::select(os, censOS, all_of(genelist))

# Relevel factors
df <- as.data.frame(unclass(df),stringsAsFactors=TRUE)
df[genelist] <- lapply(df[genelist],relevel,ref="low")   

x <- df[,3:length(df)]
x <- model.matrix( ~ ., x)

#### 2. CORRELATION ANALYSIS ####
# A. Calculate correlation matrix by Pearson method
corr_pearson <- rcorr(x, type = c("pearson"))
res <- corr_pearson$r[-1,-1] #results
pval <- corr_pearson$P[-1,-1] #p-values
corr_results <- flattenCorrMatrix(corr_pearson$r, corr_pearson$P)
summary(corr_results)

# B. Plot correlation results
corrplot(res, type="lower", method = "color", order="hclust", insig = "label_sig",
         tl.col = "black", title = "Correlation matrix (Pearson)",
         col = colorRampPalette(c("purple", "white", "orange"))(100),
         is.corr = TRUE, mar=c(1,0,1,0), p.mat =pval, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9)

#write.csv(corr_results, "3_results/9_ml/1_multivariant/1_data/cens_corr_mat_pearson.csv")

#### 3. MULTIVARIANT ANALYSIS ####
## A. CROSS VALIDATION GLMNET ##
cv_fit <-  cv.glmnet(x, # X matrix
                     Surv(df$os,df$censOS), # create survival object from the data
                     alpha = 1, # lasso: alpha = 1; ridge: alpha=0
                     family = "cox", # specify Cox PH model
                     type.measure = "C",
                     seed = 1,
                     maxit = 1000)

cv_fit2 <-  cv.glmnet(x, # X matrix
                      Surv(df$os,df$censOS), # create survival object from the data
                      alpha = 1, # lasso: alpha = 1; ridge: alpha=0
                      family = "cox", # specify Cox PH model
                      seed = 1,
                      maxit = 1000)

# Standarize coeficcients to compare
sds <- apply(x, 2, sd)

mat_coef <- as.matrix(coef(cv_fit, s = "lambda.min"))
std_coef <- mat_coef[-1, 1] * sds [-1]
coef(cv_fit, s = "lambda.1se")
active.k = which(mat_coef != 0)
active.k.vals = mat_coef[active.k]
active.k
active.k.vals

mat_coef <- as.matrix(coef(cv_fit2, s = "lambda.min"))
std_coef <- mat_coef[-1, 1] * sds [-1]
coef(cv_fit2, s = "lambda.1se")
active.k = which(mat_coef != 0)
active.k.vals = mat_coef[active.k]
active.k
active.k.vals

plot(cv_fit)
plot(cv_fit2)

## B. COX PH BACKWARD AND FORWARD ##
formula = Surv(os, censOS) ~ .
coxph_fit <-  coxph(formula=formula, data=df)
coxph_fit

# Backward, forward and stepwise selection using AIC
coxph_fit_b <-  stepAIC(coxph_fit, direction="backward", k=2)
coxph_fit_f <-  stepAIC(coxph_fit, scope=formula(coxph_fit) ,direction="forward", k=2)
coxph_fit_s <-  stepAIC(coxph_fit ,direction="both", k=2)


summary(coxph_fit_b)
summary(coxph_fit_f)
summary(coxph_fit_s)

## C. RANDOM SURVIVAL FOREST ##
pbc_rf <- rfsrc(formula = formula, 
                data = df, 
                nsplit = 0,
                block.size = 1, 
                nodesize = 5,
                importance = TRUE,
                na.action = "na.impute", 
                ntree = 1000)

plot(pbc_rf)

# Get C-index
get.cindex(pbc_rf$yvar[,1], pbc_rf$yvar[,2], pbc_rf$predicted.oob) #0.314