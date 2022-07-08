############ TCGA BLADDER CANCER (BLCA) ANALYSIS ############

# This workflow was performed in R 4.2.0

setwd("/Users/mimiferreiro/Documents/GitHub/tfm_mUC")

# Load functions
source("1_scripts/1_functions/categorized_cox.R")

# Load packages
library(dplyr)
library(survminer)
library(survival)
library(magrittr)
library(RegParallel)
library(tidyverse)
library(gtable)
library(grid)
library(gridExtra)

# Load data
extab <- read.csv("3_results/6_TCGA_BLCA/1_data/cens_fkpms_significative_genes_TCGA_BLCA.csv")[,-1]
fData <- read.csv("3_results/6_TCGA_BLCA/1_data/filter_clinical_data_TCGA_BLCA.csv")[,-1]
gene_list <- colnames(extab[,4:ncol(extab)])

#### 3. SURVIVAL ANALYSIS ####
# Calculate optimal cutpoint
opt_cut <- surv_cutpoint(extab, 
                         time = "os", 
                         event = "censOS", 
                         variables = all_of(gene_list),
                         progressbar = F)

cox_cat <- surv_categorize(opt_cut)
cox_cat <- as_tibble(cox_cat)
cox_cat$sample <- extab$sample
cox_cat <- cox_cat %>% relocate(c(sample,os, censOS), .before = 1)

# Relevel genes variables
cox_cat <- as.data.frame(unclass(cox_cat),stringsAsFactors=TRUE)
cox_cat[,4:length(cox_cat)] <- lapply(cox_cat[,4:length(cox_cat)],relevel,ref="low")
cox_cat$sample <- as.character(cox_cat$sample)

# Preparing response data
extab <- dplyr::left_join(cox_cat, fData, by = "sample")
extab <- extab %>% relocate(c(response, stage), .after = 3)
extab$response[is.na(extab$response)] <- "NE"
extab$stage[is.na(extab$stage)] <- "not reported"
extab$response %<>% as.factor()
extab$stage %<>% as.factor()

levels(extab$stage) <- list(NR  = "not reported", I = "stage i", 
                            II = "stage ii", III = "stage iii", IV = "stage iv")

levels(extab$stage)

#write.csv(extab, "3_results/6_TCGA_BLCA/1_data/cens_categorized_significative_genes_TCGA_BLCA.csv")

# Cox Regression
results <- RegParallel(
  data = extab,
  formula = 'Surv(os, censOS) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE), #skip over columns that are linear combinations of earlier columns
  FUNtype = 'coxph',
  variables = all_of(gene_list),
  blocksize = 5, # optimize parallelism
  cores = 4,
  p.adjust = "bonferroni") # change this to BH or Bonferroni to change method 

results <- results[!is.na(results$P),]
res_filter <- results %>% filter(P <= 0.05)

# Save results
#write.csv(results, "3_results/6_TCGA_BLCA/1_data/cens_no_filter_significative_genes_TCGA_BLCA.csv")
#write.csv(res_filter, "3_results/6_TCGA_BLCA/1_data/cens_filter_p005_significative_genes_TCGA_BLCA.csv")

#### 4. PLOTTING ####
# Load data 
extab <- read.csv("3_results/6_TCGA_BLCA/1_data/cens_categorized_significative_genes_TCGA_BLCA.csv")[,-1]
km_data <- read.csv("3_results/6_TCGA_BLCA/1_data/cens_filter_p005_significative_genes_TCGA_BLCA.csv")[,-1]

## A. KM plots #probar con divisiones de 2 años
genes <- km_data$Variable
data <- extab %>% dplyr::select(sample, os, censOS, all_of(genes))

# As factor
data <- as.data.frame(unclass(data),stringsAsFactors=TRUE)
data[4:ncol(data)] <- lapply(data[,4:ncol(data)],relevel,ref="low") 
data$sample <- as.character(data$sample)

plots <- list() # generate empty list

for(gene in genes) {
  fit <- surv_fit(as.formula(paste("Surv(os, censOS)~", gene)), data=data)
  
  plots[[gene]] <- ggsurvplot(fit, 
                              pval = T,
                              conf.int = T,
                              title = gene, # Uses name of factor as plot title
                              legend.labs = c("Low", "High"),
                              palette = c("#791BC2", "#F28E2B"),
                              risk.table = T,
                              xlim = c(0, max(data$os)),
                              break.x.by = 2,
                              legend.title = "",
                              xlab = "Time (years)",
                              ylab = "Overall survival",
                              tables.theme = theme_cleantable()
  )
  print(plots[gene]$gene)
}


#### 5. STAGE IV ONLY ####
stage <- extab %>% filter(stage == "IV") %>% dplyr::select(-response)
stage <- as.data.frame(unclass(stage),stringsAsFactors=TRUE)
stage[5:ncol(stage)] <- lapply(stage[,5:ncol(stage)],relevel,ref="low") 
stage$sample <- as.character(stage$sample)

# A. Cox Regression
results <- RegParallel(
  data = stage,
  formula = 'Surv(os, censOS) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE), #skip over columns that are linear combinations of earlier columns
  FUNtype = 'coxph',
  variables = all_of(gene_list),
  blocksize = 5, # optimize parallelism
  cores = 4,
  p.adjust = "bonferroni") # change this to BH or Bonferroni to change method 

results <- results[!is.na(results$P),]
res_filter <- results %>% filter(P <= 0.05)

#write.csv(res_filter, "3_results/6_TCGA_BLCA/3_stageIV/cens_results_stage_IV.csv")


## B. KM plots
genes <- res_filter$Variable
plots <- list() # generate empty list

for(gene in genes) {
  fit <- surv_fit(as.formula(paste("Surv(os, censOS)~", gene)), data=stage)
  
  plots[[gene]] <- ggsurvplot(fit, 
                              pval = T,
                              conf.int = T,
                              title = gene, # Uses name of factor as plot title
                              legend.labs = c("Low", "High"),
                              palette = c( "#791BC2", "#F28E2B"),
                              risk.table = T,
                              xlim = c(0, max(stage$os)),
                              break.x.by = 2,
                              legend.title = "",
                              xlab = "Time (years)",
                              ylab = "Overall survival",
                              tables.theme = theme_cleantable()
  )
  print(plots[gene]$gene)
}