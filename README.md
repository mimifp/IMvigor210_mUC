# IMvigor210_mUC
This project will consist of carrying out an analysis of gene expression in tumor tissue from patients diagnosed with urothelial carcinoma, with the aim of identifying genes not only with prognostic impact, but also correlated with response to immunotherapy. Those final candidate genes will be validated as prognostic and/or predictive biomarkers in future studies. With the aim of analyzing the weight of each biomarker in the prognosis of the disease, various prognostic models will be optimized. On the other hand, cellular deconvolution methods will be apply to those groups of patients defined by our genes of interest in order to better characterize their tumor microenviroment. The same approach will be follow for GSEA.

This repository contains all the scripts and data necessary to carry out this analysis.
* 1_scripts
  + 1_functions
    + categorized_cox.R: Categorize gene expression as "high" or "low" according to optimal cutoff point.
    + counts_to_tpm.R: Transform raw counts in trancripts per million (TPMs).
    + cox_regression.R: Calculate optimal cutoff point and do survival analysis by Cox PH regression.
    + flattenCorrMatrix.R: Calculate flatten correlation matrix.
    + remove_duplicated_genes.R: Remove duplicated genes leaving leaving the one with the highest expression.
    + tpm_coxdata.R: Prepare tpm matrix for categorization and survival analysis. Add clinical data.
  + 2_pipeline
    + 1_download_data.R: download IMvigor210_mUC expression, phenotype and annotation data.
    + 2_preprocess_data.R: filter and normalize (TPMs) expression data.
    + 3_survival_analysis.R: survival analysis by Cox PH regression and calculate distribution statistics and logistic regression.
    + 4_kaplan_meier_plots.R: plot Kaplan-Meier curves and hazard ratio tables.
    + 5_TCGA_BLCA_analysis.R: survival analysis by Cox PH regression and plot Kaplan-Meier curves for TCGA_BLCA data.
    + 6_deconvolution.R: deconvolution with EPIC and analysis of CD-T cells exhaustion markers.
    + 7_gsea.R: Gen Set Enrichment Anlaysis (GSEA) with clusterProfiler.
    + 8_multivariant.R: multivariant analysis with differents models.
* 2_data: necessary data to run analysis with each dataset.
  + 1_IMvigor210
  + 2_TCGA_BLCA
  + 3_hallmark_geneset
    + exp_mat.txt: expression matrix in [.txt format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#TXT:_Text_file_format_for_expression_dataset_.28.2A.txt.29). First column must be gene name (Name), second an description (Description) and then one column per sample (Sample1, Sample2 ... SampleX)
    + h.all.v7.5.1.symbols.gmt: gene set that summarize and represent specific well-defined biological states or processes and display coherent expression in [.gmt format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29).
    + phenotype_x.cls: matrix in [.cls format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#CLS:_Categorical_.28e.g_tumor_vs_normal.29_class_file_format_.28.2A.cls.29) indicating the samples with high or low phenotype for each of the genes.
* 3_results: data and plots resulting from different parts of analysis
  + 1_filtering
  + 2_survival_analysis
  + 3_distribution_statistics
  + 4_logistic_regression
  + 5_km_plots
  + 6_TCGA_BLCA
  + 7_deconvolution
  + 8_gsea
  + 9_multivariant
