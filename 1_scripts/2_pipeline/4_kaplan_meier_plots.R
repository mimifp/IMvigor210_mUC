############ KAPLAN-MEIER CURVES ############

# This workflow was performed in R 4.2.0

# Load packages 
library(dplyr)
library(readr)
library(survival)
library(survminer)
library(patchwork)
library(svglite)
library(gridExtra)

# Change working directory
setwd("/Users/mimiferreiro/Documents/GitHub/tfm_mUC")

# Set input/output paths
input_file <- "3_results/2_survival_analysis/cens_regression_data_coxdata_sd.csv"
gene_list <- read.csv("3_results/2_survival_analysis/cens_significative_genes_cox_optimalcut_bonferroni_res_sd.csv")[,2]

# Load survival data. 
data <- read_csv(file.path(input_file))[,-4]

#### 1. PREPARE DATA ####
# Generate list of categorical factors (genes or other)
filter_variables <- c("sample_id", "os", "censOS", gene_list)
data <- data %>% select_(.dots = filter_variables)
data <- as.data.frame(unclass(data), stringsAsFactors = TRUE)
data[4:ncol(data)] <- lapply(data[,4:ncol(data)],relevel,ref="low")  
data$sample_id <- as.character(data$sample_id)

#### 2. KM PLOTS ####
# Generate combined KM figure
plots <- list() # generate empty list

for(gene in gene_list) {
  fit <- surv_fit(as.formula(paste("Surv(os, censOS)~", gene)), data=data)
  
  plots[[gene]] <- ggsurvplot(fit, 
                              pval = T,
                              conf.int = T,
                              title = gene, # Uses name of factor as plot title
                              legend.labs = c("Low", "High"),
                              palette = c("#791BC2", "#F28E2B"),
                              risk.table = T,
                              xlim = c(0, max(data$os)),
                              legend.title = "",
                              xlab = "Time (months)",
                              ylab = "Overall survival",
                              tables.theme = theme_cleantable()
  )
  print(plots[gene]$gene)
}