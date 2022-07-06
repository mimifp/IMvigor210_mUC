############# FILTER AND NORMALIZE DATA #################

# This workflow was performed in R 4.2.0

# Load packages
library(tidyverse)
library(dplyr)
library(genefilter)
library(ggplot2)
library(gridExtra)
library(viridis)
library(hrbrthemes)
library(magrittr)
library(RColorBrewer)
library(ggpubr)

# Change working directory
setwd("../../2_data/1_IMvigor210")

# Load data
exmat <- read.csv("exmat_censored_IMvigor210.csv") #31085 genes
pData <- read.csv("pData_IMvigor210.csv")
fData <- read.csv("fData_IMvigor210.csv")

# Change rownames
exmat <- exmat %>% column_to_rownames(var = "X")
pData <- pData %>% dplyr::rename(sample_id = "X")

#### 1. QC ANALYSIS ####
# A. BOXPLOT
# Explore distribution of counts. 
boxplot_exmat <- stack(as.data.frame(log(exmat)+1))

bp1 <- ggplot(boxplot_exmat) + 
  geom_boxplot(aes(x = ind, y = values, alpha = 1)) +
  theme(
    legend.position="none",
    plot.title = element_text(size=15, hjust = 0.5),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()
  ) +
  ggtitle("Raw data") +
  ylab(expression('Log'[2]~'Read counts')) +
  xlab("Samples")

# B. BARPLOT
# Check if any samples need to be discarded based on the number of genes detected.
pars=list(par(mar=c(10,6,4,2)))
{barplot(colSums(exmat>35),
         main="Number of genes detected in samples",
         las=2,
         ylim= c(0,16000)
)
  abline(h=median(colSums(exmat>35)), 
         col = "red", 
         lwd = 1.5)
  text(410, 14750, "Median", col = "red") #abline is median
}

# C. HISTOGRAM
# Check gene expression above samples (at least 10%)
hist(rowSums(exmat>35), main = "Gene expression", xlab = "")

# Mutate fData
colnames(pData)[1] <- "sample_id"

#### 2. CLEAN UP AND NORMALIZATION ####
## Data clean up
# Check if there are some NA or duplicated values
table(is.na(exmat))
table(duplicated(rownames(exmat)))
table(duplicated(names(exmat)))

## Calculate TPMs
flen <- as.numeric(fData[,4])
tpm <- counts_to_tpm(exmat, flen)

#write.csv(tpm, "../../3_results/1_filtering/1_data/exmat_tpm.csv")

#### 3. FILTER DATA ####
## First of all we applied one filter in whole TPM matrix:
#  - Remove genes with 0 values in at least 95% of samples
tpm <- tpm[rowMeans(tpm == 0) <= 0.95,] # 22517 genes

# After that, we follow two different routes to filter TPMs:
#  1. Drop genes with IQR <= 0.5 and include only genes with log(TPM+1) >= 2 in at least 10% of samples
#  2. Drop genes with s.d <= 1 and inlude only genes with log(TPM+1) >= 2 in at least 10% of samples

## Filter 1. IQR < 0.5 + log(TPM+1) > 2
f1 <- function(x) ( IQR(x) > 0.5 )
selected <- genefilter(tpm, f1)
tpm_iqr5 <- tpm[selected,] #12376 genes

tpm_iqr5_log <- log(tpm_iqr5 + 1)
tpm_iqr5_log10 <- tpm_iqr5_log[apply(tpm_iqr5_log, 1, function(x) sum(x > 2) / length(x)) >= 0.1, ] # 8280 genes

dim(tpm_iqr5)
dim(tpm_iqr5_log10)

# Save data
#write.csv(tpm_iqr5, "../../3_results/1_filtering/1_data/tpm_iqr5_filter_cens.csv")
#write.csv(tpm_iqr5_log10, "../../3_results/1_filtering/1_data/tpm_iqr5_log10_filter_cens.csv")

## Filter 2. S.D < 1 + log(TPM+1) > 2
tpm_sd <- tpm[!apply(tpm, 1, sd, na.rm=TRUE) <= 1,] # 11398 genes

tpm_sd_log <- log(tpm_sd + 1)
tpm_sd_log10 <- tpm_sd_log[apply(tpm_sd_log, 1, function(x) sum(x > 2) / length(x)) >= 0.1, ] # 8282 genes

dim(tpm_sd)
dim(tpm_sd_log10)

# Save data
#write.csv(tpm_sd, "../../3_results/1_filtering/1_data/tpm_sd_filter_cens.csv")
#write.csv(tpm_sd_log10, "../../3_results/1_filtering/1_data/tpm_sd_log10_filter_cens.csv")

#### 4. COMPARISION PLOTS ####
# Compare raw data with TPM normalized data (ggplot2)
boxplot_iqr <- stack(as.data.frame(tpm_iqr5_log10))
boxplot_sd <- stack(as.data.frame(tpm_sd_log10))

bp2 <- ggplot(boxplot_iqr) + 
  geom_boxplot(aes(x = ind, y = values, alpha = 1)) +
  theme(
    legend.position="none",
    plot.title = element_text(size=15, hjust = 0.5),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()
  ) +
  ggtitle("IQR data") +
  ylab(expression('Normalized counts')) +
  xlab("Samples")

bp3 <- ggplot(boxplot_sd) + 
  geom_boxplot(aes(x = ind, y = values, alpha = 1)) +
  theme(
    legend.position="none",
    plot.title = element_text(size=15, hjust = 0.5),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()
  ) +
  ggtitle("SD data") +
  ylab(expression('Normalized counts')) +
  xlab("Samples")

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE),
       widths=c(1,1), heights=c(1,1))
bp1
bp2
bp3