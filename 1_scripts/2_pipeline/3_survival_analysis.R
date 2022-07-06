############# SURVIVAL ANALYSIS AND STATISTICS #################

# Load packages
library(tidyverse)
library(dplyr)
library(survminer)
library(RegParallel)
library(survival)
library(ggplot2)
library(gridExtra)
library(viridis)
library(hrbrthemes)
library(magrittr)
library(RColorBrewer)
library(ggpubr)

# Load functions
source("../1_functions/tpm_coxdata.R")
source("../1_functions/cox_regression.R")
source("../1_functions/categorized_cox.R")

#### 1. SURVIVAL ANALYSIS ####
## Preparing clincal data and Cox regression data
# We merge clinical data with matrix obtained after apply each filter. For each 
# gene (columns) in each sample (rows) we have a expression value, an overall 
# survival value (os) and patient status in the moment of data collection 
# (censOS; 1 = dead, 0 = alive).

# Preparing clinical data
clinical <- pData[c("os", "censOS")]
clinical$sample_id <- colnames(exmat)

# For filter 1
coxdata_iqr <- tpm_coxdata(tpm_iqr5_log10)

# For filter 2
coxdata_sd <- tpm_coxdata(tpm_sd_log10)

## Calculate optimal cutpoints and Cox Regression. 
# After that, the p-values obtained are filtered so that we are only left with the values that are < 0.05.

# For filter 1
res_iqr <- cox_regression(coxdata_iqr)

# For filter 2
res_sd <- cox_regression(coxdata_sd)

## Process and save survival analysis results
results <- c("res_iqr", "res_sd")

for (item in results) {
  res <- get(item)
  assign(as.character(item), res, envir= .GlobalEnv)
  res <- res[!is.na(res$P),]
  res_fil <- res %>% filter(P.adjust <= 0.05) 
  out_path <- paste0("../../3_results/2_optimal_cutpoint/cens_significative_genes_cox_optimalcut_bonferroni_",paste(item),".csv")
  # Saving results
  write.csv(res_fil, out_path)
}

## Preparing dichotomized data por further analysis
dicho_dat <- c("coxdata_iqr", "coxdata_sd")

clinical_2 <- pData[c("Best.Confirmed.Overall.Response", "sample_id")]

for (df in dicho_dat){
  data <- get(df)
  assign(as.character(df), data, envir= .GlobalEnv)
  cox_cat <- categorize_cox(data)
  # Merge clinical data with results of calculate optimal cutpoint
  cox_cat$sample_id <- clinical_2$sample_id
  reg_data <- merge(cox_cat, clinical_2, by = "sample_id")
  reg_data <- reg_data %>% relocate(c(sample_id, os, censOS, 'Best.Confirmed.Overall.Response'), .before = 1) %>% 
    dplyr::rename(response = 'Best.Confirmed.Overall.Response')
  #  # Write results
  out_path <- paste0("../../3_results/2_optimal_cutpoint/cens_regression_data_",paste(df),".csv")
  write.csv(reg_data, out_path, row.names = F)
}

#### 2. DISTRIBUTION STATISTICS ####
## Calculate statistics for DC and Response by separate. 
filters <- c("iqr", "sd")

name <- "sd"

for (name in filters){
  data <- read.csv(paste0("../../3_results/2_optimal_cutpoint/cens_regression_data_coxdata_",paste0(name),".csv"))
  sig_genes <- read.csv(paste0("../../3_results/2_optimal_cutpoint/cens_significative_genes_cox_optimalcut_bonferroni_res_",paste(name),".csv"))[-1]
  
  ## A. PREPARE DATA ##
  # Set response levels
  data$response <- factor(data$response, levels = c("CR", "PR", "SD", "PD", "NE"))
  
  # DC = Disease Control ; NDC = No Disease Control
  data <- data %>%
    mutate(dc = case_when(
      response == "CR" ~ "DC",
      response == "PR" ~ "DC",
      response == "SD" ~ "DC",
      TRUE ~ "NDC"
    ))
  
  # R = Response; NR = No Response
  data <- data %>%
    mutate(response = case_when(
      response == "CR" ~ "R",
      response == "PR" ~ "R",
      TRUE ~ "NR"
    ))
  
  data$dc %<>% factor()
  data$response %<>% factor()
  
  # Check factor levels
  levels(data$dc)
  levels(data$response)
  
  # Reversing factor levels
  data$dc <- fct_rev(data$dc)
  levels(data$dc)
  
  data <- data %>% relocate(dc, response, .after = 4)
  
  # Converting genes variables to factor
  data <- as.data.frame(unclass(data),stringsAsFactors=TRUE)
  data[,6:length(data)] <- lapply(data[,6:length(data)],relevel,ref="low")
  data$sample_id <- as.character(data$sample_id)
  
  # Selecting significative genes and variables
  genes <- sig_genes$Variable
  
  # Converting into long format
  datalong <- gather(data, gene_id, level, all_of(genes))
  datalong$gene_id %<>% factor()
  datalong$gene_id %<>% order()
  datalong$level %<>% factor()
  datalong$level %<>% fct_rev()
  
  ## B. DISTRIBUTION STATISTICS ##
  sink(file = paste0("../../3_results/3_distribution_statistics/distribution_statistics_output_",paste(name),".txt"))
  print(paste(name))
  for (gene in genes) {
    data_gene_dc <- datalong %>% 
      filter(gene_id == gene) %>% 
      dplyr::select(level, dc) %>% table
    chisq_results_dc <- chisq.test(data_gene_dc)
    
    data_gene_res <- datalong %>% 
      filter(gene_id == gene) %>% 
      dplyr::select(level, response) %>% table
    chisq_results_res <- chisq.test(data_gene_res)
    
    assign(paste("results_", gene, sep=''), list(chisq_results_dc, chisq_results_res))
    print(paste("results_", gene, sep=''))
    print(get(paste("results_", gene, sep='')))
  }
  sink(file = NULL)
}

#### 3. PLOTS ####
## A. STACKED BARPLOTS
customcolors <- c("#791BC2", "#F28E2B")

plot <- ggplot(datalong,
               aes(x = dc, #change dc for response
                   fill = level)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(~gene_id)+
  scale_fill_manual(values = customcolors) +
  ylab("Percentage") + xlab(paste0("Clinical outcome (DC)")) + #change (DC) for (R)
  labs(fill="Expression level") +
  theme(
    panel.background = element_blank(),
    strip.background = element_blank(),
    text = element_text(size = 12),
    strip.text.x = element_text(face = "bold", size = 6.5),
    legend.text = element_text(size=12)
  )

print(plot)