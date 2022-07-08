############ DECONVOLUTION WITH EPIC ############

# This workflow was performed in R 4.2.0

# Install
#install.packages("devtools")
#devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)

# Change working directory
setwd("/Users/mimiferreiro/Documents/GitHub/tfm_mUC")

# Load packages
library(dplyr)
library(tibble)
library(readr)
library(immunedeconv)
library(tidyverse)
library(magrittr)
library(nortest)
library(gridExtra)


#### 1. DECONVOLUTION WITH EPIC 0.1 ####
exmat <- read.csv("1_data/tpm_sd_log10_filter.csv") #need NO CENSORED data
output_file <- "3_results/7_deconvolution/1_data/epic_results.csv"

# Load expression matrix
if (sum(duplicated(exmat[1])) > 0) {
  stop('Expression matrix contains duplicated genes! Check gene_id or first 
       column')
}

# Genes as rownames (required)
exmat <- column_to_rownames(exmat, var = colnames(exmat[1]))

# Revert log2
exmat <- exp(exmat)-1 

# Deconvolution
epic_results <- deconvolute_epic(exmat, 
                                 tumor = T, # The sample is tumoral?
                                 scale_mrna = F) # We don't have reference RNA

# Process results
epic_results <- data.frame(epic_results)
epic_results <- epic_results %>% rownames_to_column(var = "cell_type")

#write.csv(epic_results, file.path(output_file))

#### 2. PLOT DECONVOLUTION RESULTS ####
# Load results
epic <- read.csv("3_results/7_deconvolution/1_data/epic_results.csv")
gene_list <- read.csv("3_results/1_filtering/1_data/cens_significative_genes_cox_optimalcut_bonferroni_res_iqr.csv")[,2]
extab <- read.csv("3_results/2_optimal_cutpoint/cens_regression_data_coxdata_sd.csv")

extab <- extab %>% dplyr::select(sample_id, all_of(gene_list))

# Prepare epic data
cell_types <- epic[,1]
epic <- as.data.frame(t(epic))
colnames(epic) <- epic[1,]
epic <- epic[-1,]
epic$sample_id <- rownames(epic)

# Prepare datalong
extab_cells <- left_join(extab, epic, by = "sample_id")

datalong <- gather(extab_cells, gene_id, level, all_of(gene_list))
datalong <- gather(datalong, cell_type, cell_fraction,  all_of(cell_types))

datalong$level %<>% as.factor
datalong$gene_id %<>% as.factor
datalong$cell_type %<>% as.factor
datalong$cell_fraction %<>% as.numeric()

# Plots
for (type in cell_types){
  df <- datalong %>% filter(cell_type == type)
  
  plot <- df %>%
    group_by(gene_id) %>%
    ggplot(aes(x=gene_id, y=cell_fraction, fill = level)) +
    geom_boxplot() +
    scale_fill_manual(values=c("#F28E2B", "#791BC2")) +
    scale_y_continuous()+
    ggtitle(type) +
    labs(x = "Gene", y = "Cell fraction", fill="Expression level") +
    theme(text = element_text(size = 12),
          legend.text = element_text(size=12),
          plot.title = element_text(hjust = 0.5)) 
  print(plot)
  # Check normality of variables
  hist(df$cell_fraction, main = type, xlab = "Cell fraction", col = "#F28E2B")
  # Lilliefors Normality Test (Kolmogorov-Smirnov)
  test <- lillie.test(df$cell_fraction)
  print(paste0("Normality test for ", type))
  print(test)
}

#### 3. DISTRIBUTION STATISTICS ####
# Wilcoxon Test
sink(file = paste0("3_results/7_deconvolution/1_data/cens_MWUtest_output.txt"))
for (type in cell_types){
  print(paste("ANALYZING", type, sep = " "))
  for (gene in gene_list){
    print(paste("Mann-Whitney U Wilcoxon test for", gene, sep = " "))
    df <-  datalong %>% filter(cell_type == type, gene_id == gene)
    res <- wilcox.test(cell_fraction ~ level, data = df)
    print(res)
  }
}
sink(file = NULL)

#### 4. SIGNIFICATIVE CELLTYPES PLOTS ####
# Load MWUtest results table
data <- read.csv("3_results/7_deconvolution/1_data/MWUtest_table.csv", sep = ";")

data <- dplyr::rename(data, gene_id = X)
cell_types <- names(data[,-1])

# Prepare datalong
tablelong <- gather(data, cell_type, pvalue, all_of(cell_types))
datalong <- dplyr::left_join(datalong, tablelong, by = c("gene_id", "cell_type"))

datalong$gene_id %<>% as.factor
datalong$cell_type %<>% as.factor
datalong$level %<>% fct_rev()

plot <- list()

for (type in cell_types){
  df <- datalong %>% filter(cell_type == type, pvalue <= 0.05)
  ylim_real = boxplot.stats(df$cell_fraction)$stats[c(1, 5)]
  
  # Plot only if p-value is < 0.1
  if (dim(df)[1] != 0) {
    plot <- df %>%
      group_by(gene_id) %>%
      ggplot(aes(x=gene_id, y=cell_fraction, fill = level)) +
      geom_boxplot() +
      scale_fill_manual(values=c("#791BC2", "#F28E2B")) +
      scale_y_continuous()+
      ggtitle(type) +
      labs(x = "Gene", y = "Cell fraction", fill="Expression level") +
      theme(text = element_text(size = 14),
            legend.text = element_text(size=12),
            plot.title = element_text(hjust = 0.5)) +
      stat_compare_means(method = "wilcox.test", aes(label = paste0("p = ", ..p.format..)), label.y = ylim_real[2]*1.1) +
      coord_cartesian(ylim = ylim_real*1.2)  # Reduce ylim in a 20%
    print(plot)
  }
}

#### 5. CD8T CELLS MARKERS ####
# Load data
tpm <- read.csv("3_results/1_filtering/1_data/exmat_tpm.csv")
sig_genes <- read.csv("3_results/1_filtering/1_data/cens_significative_genes_cox_optimalcut_bonferroni_res_sd.csv")[,2]
dicho_data <- read.csv("3_results/2_optimal_cutpoint/cens_regression_data_coxdata_sd.csv")
markers <- c("gene_4341", "gene_8056", "gene_9979", "gene_22450", "gene_28173")

# Prepare data
tpm <- tpm %>% column_to_rownames(var = "X")
tpm <- log(tpm+1)
tpm <- as.data.frame(t(tpm))

# Filter data
tpm_fil <- tpm %>% dplyr::select(all_of(markers))
tpm_fil <- tibble::rownames_to_column(tpm_fil, "sample_id")

dicho_data <- dicho_data %>% dplyr::select(sample_id, all_of(sig_genes))

# Converting genes variables to factor
dicho_data <- as.data.frame(unclass(dicho_data),stringsAsFactors=TRUE)
dicho_data[,5:length(dicho_data)] <- lapply(dicho_data[,5:length(dicho_data)],relevel,ref="low")
dicho_data$sample_id <- as.character(dicho_data$sample_id)

# Convert into long format
datalong_dicho <- gather(dicho_data, gene_id, level, all_of(sig_genes))
datalong_tpm <- gather(tpm_fil, marker_id, expression_value, all_of(markers))

data <- left_join(datalong_dicho, datalong_tpm, by = "sample_id")
data$level <- as.factor(data$level)
data$level %<>% fct_rev()

# Unanonimize markers
data$marker_id[data$marker_id == "gene_4341"] <- "CTLA4"
data$marker_id[data$marker_id == "gene_8056"] <- "HAVCR2"
data$marker_id[data$marker_id == "gene_9979"] <- "LAG3"
data$marker_id[data$marker_id == "gene_22450"] <- "PDCD1"
data$marker_id[data$marker_id == "gene_28173"] <- "TIGIT"


# Boxplot
for (gene in sig_genes){
  data_sel <- data %>% dplyr::filter(gene_id == gene)
  
  ylim_real = boxplot.stats(data_sel$expression_value)$stats[c(1, 5)]
  plot <- data_sel %>%
    dplyr::group_by(marker_id) %>%
    ggplot(aes(x=marker_id, y=expression_value, fill = level)) +
    geom_boxplot() +
    scale_fill_manual(values=c("#791BC2", "#F28E2B")) +
    scale_y_continuous()+
    ggtitle(gene) +
    labs( y = "log2(TPM+1)", fill="Expression level") +
    theme(text = element_text(size = 14),
          axis.title.x=element_blank(),
          legend.text = element_text(size=12),
          plot.title = element_text(hjust = 0.5)) +
    stat_compare_means(method = "wilcox.test", aes(label = paste0("p = ", ..p.format..)), label.y = ylim_real[2]*1.1) +
    coord_cartesian(ylim = ylim_real*1.2)  # Reduce ylim in a 20%
  print(plot)
}