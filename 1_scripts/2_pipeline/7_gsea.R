############ Gen Set Enrichment Analysis (GSEA) with clusterProfiler ############

# This workflow was performed in R 4.2.0

# Install packages
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("clusterProfiler")

# Load packages
library(dplyr)
library(plyr)
library(ggplot2) # plotting
library(gridExtra) # combining ggplots
library(DESeq2) # rna-seq
library(edgeR) # rna-seq
library(tidyverse)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)

# Load data 
exmat <- read.csv("../../2_data/1_IMvigor210/exmat_IMvigor210.csv")
pData <- read.csv("../../3_results/2_optimal_cutpoint/regression_data_coxdata_sd.csv")[,-c(2,3,4)] #need NO CENSORED data
tpm <- read.csv("../../3_results/1_filtering/1_data/tpm_sd_log10_filter.csv") #need NO CENSORED data
sig_genes <- read.csv("../../3_results/1_filtering/1_data/significative_genes_cox_optimalcut_bonferroni_res_sd.csv")[,"Variable"] #need NO CENSORED DATA

# Check rownames
exmat <- exmat %>% remove_rownames %>% column_to_rownames(var="gene_id")
pData <- pData %>% remove_rownames %>% column_to_rownames(var="sample_id")
tpm <- tpm %>% remove_rownames %>% column_to_rownames(var="X")

#### 1. FILTERING DATA ####
# We use our filtered genes with sd method
tpm$gene_id <- rownames(tpm)
exmat$gene_id <- rownames(exmat)

# Keep genes in filter tpm with values of raw count matrix
genes <- intersect(rownames(tpm),rownames(exmat))
cfiltered <- exmat[genes,]
cfiltered <- cfiltered[,-length(cfiltered)]

# Check if columns and rows are equal
cfiltered <- cfiltered[,order(colnames(cfiltered))]
pData <- pData[order(rownames(pData)),]
all.equal(colnames(cfiltered),rownames(pData))

# Relevel pData for high vs. low comparision
pData <- as.data.frame(unclass(pData),stringsAsFactors=TRUE)
pData <- catcolwise(function(x) relevel(x, "low"))(pData)

#### 2. DIFFERENTIAL EXPRESSION ANALYSIS (DEA) ####
for (gene in sig_genes){
  # Filter metadata
  pData_filt <- pData %>% dplyr::select(all_of(gene))
  # Create deSeq2 object
  d <- DESeqDataSetFromMatrix(countData=cfiltered,colData=pData_filt,design=formula(paste("~", gene)))
  # Normalise by deseq method (median-of-ratios). 
  d <- DESeq2::estimateSizeFactors(d,type="ratio", locfunc = stats::median)
  # Gene dispersion
  d <- DESeq2::estimateDispersions(d)
  # Testing
  dg <- nbinomWaldTest(d)
  # Results
  res <- results(dg,name=paste0(paste(gene),"_high_vs_low"),alpha=0.05)
  hist(res$pvalue[res$baseMean>1],main=paste0(gene," Pval distribution"),xlab="P-values")
  # Save results
  df_res <- as.data.frame(res)
  write.csv(df_res, paste0("../../3_results/8_gsea/1_data/dif_exp_analysis_",paste(gene),".csv"))
}  

##### 3. GSEA ####
for (gene in sig_genes){
  df_res <- read.csv(paste0("../../3_results/8_gsea/1_data/dif_exp_analysis_",paste(gene),".csv"))
  df_res <- df_res %>% remove_rownames %>% column_to_rownames(var="X")
  
  # Load our data
  geneList <- df_res$log2FoldChange
  names(geneList) <- rownames(df_res)
  geneList <- sort(geneList, decreasing = TRUE)
  
  gsea_BH <- gseGO(geneList  = geneList,
                   ont       = "BP",
                   OrgDb     = org.Hs.eg.db,
                   keyType   = 'SYMBOL',
                   verbose   = T, # print message
                   minGSSize = 10,
                   maxGSSize = 500,
                   pAdjustMethod = "BH",
                   seed      = T,
                   eps = 0)
  
  # GSEA dotplot
  print(
    dotplot(gsea_BH, showCategory=10, font.size = 8) + ggtitle(paste0("GSEA for ",gene))
  )
  
  # Save results
  gsea_res <- gsea_BH@result
  write.csv(gsea_res, paste0("../../3_results/8_gsea/1_data/gsea_results_",paste(gene),".csv"))
}


#### 4. KEGG pathways ####
for (gene in sig_genes){
  # Load data
  df_res <- read.csv(paste0("../../3_results/8_gsea/1_data/dif_exp_analysis_",paste(gene),".csv"))
  df_res <- df_res %>% remove_rownames %>% column_to_rownames(var="X")
  
  # Convert Symbol to EntrezID
  symbol <- rownames(df_res)
  entrez <- mapIds(org.Hs.eg.db, symbol, 'ENTREZID', 'SYMBOL')
  
  # Check and drop NA and duplicates
  table(is.na(entrez))
  entrez <- entrez[!is.na(entrez)]
  table(duplicated(entrez))
  
  geneList_kegg = df_res$log2FoldChange
  names(geneList_kegg) = entrez
  geneList_kegg = sort(geneList_kegg, decreasing = TRUE)
  
  gsea_kegg = gseKEGG(geneList     = geneList_kegg,
                      organism     = 'hsa',
                      keyType      = "ncbi-geneid",
                      pvalueCutoff = 0.05,
                      verbose      = TRUE,
                      eps = 0)
  dim(gsea_kegg@result)
  
  # Save results
  kegg_res <- gsea_kegg@result
  write.csv(kegg_res, paste0("../../3_results/8_gsea/1_data/kegg_results_",paste(gene),".csv"))
}

# Set genes that have enrichment in KEGG pathways
sig_genes <- c("ERCC4", "MLLT3")

#### 5. PLOTS ####
for (gene in sig_genes){
  # Dotplot of top enriched functions 
  # Load data
  gsea_BH <- read.csv(paste0("../../3_results/8_gsea/1_data/gsea_results_",paste(gene),".csv"))
  gsea_KEGG <- read.csv(paste0("../../3_results/8_gsea/1_data/kegg_results_",paste(gene),".csv"))
  gsea_BH <- gsea_BH %>% remove_rownames %>% column_to_rownames(var="X")
  gsea_KEGG <- gsea_KEGG %>% remove_rownames %>% column_to_rownames(var="X")
  
  
  df_gsea_BH <- as.data.frame(gsea_BH)
  df_gsea_KEGG <- as.data.frame(gsea_KEGG)
  
  go_dot <- df_gsea_BH[order(df_gsea_BH$NES, decreasing = T), ]
  go_dot <- go_dot[c(1:10, (nrow(go_dot)-9):nrow(go_dot)), ]
  go_dot$Cluster <- as.factor("GO")
  go_dot$Up_Down <- ifelse(go_dot$NES > 0, "Enriched", "Depleted")
  
  kegg_dot <- df_gsea_KEGG[order(df_gsea_KEGG$NES, decreasing = T), ]
  kegg_dot$Cluster <- as.factor("KEGG")
  kegg_dot$Up_Down <- ifelse(kegg_dot$NES > 0, "Enriched", "Depleted")
  
  df <- rbind(go_dot, kegg_dot)
  
  # plot
  print(
    ggplot(df, aes(Description, Up_Down)) +
      geom_point(aes(color = NES, size = setSize)) +  
      facet_grid(~Cluster) + 
      coord_flip() +
      scale_color_gradient2(low="darkviolet", high="darkorange", midpoint=0, 
                            limits=c(min(df$NES), max(df$NES))) +
      labs(title = paste0(gene))
  )
}