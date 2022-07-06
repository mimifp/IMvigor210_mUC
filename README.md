# IMvigor210_mUC
This project will consist of carrying out an analysis of gene expression in tumor tissue from patients diagnosed with urothelial carcinoma, with the aim of identifying genes not only with prognostic impact, but also correlated with response to immunotherapy. Those final candidate genes will be validated as prognostic and/or predictive biomarkers in future studies. With the aim of analyzing the weight of each biomarker in the prognosis of the disease, various prognostic models will be optimized. On the other hand, cellular deconvolution methods will be apply to those groups of patients defined by our genes of interest in order to better characterize their tumor microenviroment. The same approach will be follow for GSEA.

### 1. CONTENTS
This repository contains all the scripts and data necessary to carry out this analysis:
* 1_scripts
  + 1_functions
    + [categorized_cox.R](1_scripts/1_functions/categorized_cox.R): Categorize gene expression as "high" or "low" according to optimal cutoff point.
    + [counts_to_tpm.R](1_scripts/1_functions/counts_to_tpm.R): Transform raw counts in transcripts per million (TPMs).
    + [cox_regression.R](1_scripts/1_functions/cox_regression.R): Calculate optimal cutoff point and do survival analysis by Cox PH regression.
    + [flattenCorrMatrix.R](1_scripts/1_functions/flattenCorrMatrix.R): Calculate flatten correlation matrix.
    + [remove_duplicated_genes.R](1_scripts/1_functions/remove_duplicated_genes.R): Remove duplicated genes leaving leaving the one with the highest expression.
    + [tpm_coxdata.R](1_scripts/1_functions/tpm_coxdata.R): Prepare tpm matrix for categorization and survival analysis. Add clinical data.
  + 2_pipeline
    + [1_download_data.R](1_scripts/2_pipeline/1_download_data.R): download IMvigor210_mUC expression, phenotype and annotation data.
    + [2_preprocess_data.R](1_scripts/2_pipeline/2_preprocess_data.R): filter and normalize (TPMs) expression data.
    + [3_survival_analysis.R](1_scripts/2_pipeline/3_survival_analysis.R): survival analysis by Cox PH regression and calculate distribution statistics and logistic regression.
    + [4_kaplan_meier_plots.R](1_scripts/2_pipeline/4_kaplan_meier_plots.R): plot Kaplan-Meier curves and hazard ratio tables.
    + [5_TCGA_BLCA_analysis.R](1_scripts/2_pipeline/5_TCGA_BLCA_analysis.R): survival analysis by Cox PH regression and plot Kaplan-Meier curves for TCGA_BLCA data.
    + [6_deconvolution.R](1_scripts/2_pipeline/6_deconvolution.R): deconvolution with [EPIC](https://github.com/GfellerLab/EPIC) and analysis of CD-T cells exhaustion markers.
    + [7_gsea.R](1_scripts/2_pipeline/7_gsea.R): Gen Set Enrichment Anlaysis (GSEA) with [clusterProfiler](https://github.com/YuLab-SMU/clusterProfiler).
    + [8_multivariant.R](1_scripts/2_pipeline/8_multivariant.R): multivariant analysis with differents models.
* 2_data: necessary data to run analysis with each dataset.
  + 1_IMvigor210
    + [exmat_censored_IMvigor210.csv](2_data/1_IMvigor210/exmat_censored_IMvigor210.csv): expression matrix
    + [fData_IMvigor210.csv](2_data/1_IMvigor210/fData_IMvigor210.csv): feature data
    + [pData_IMvigor210.csv](2_data/1_IMvigor210/pData_IMvigor210.csv): phenotype data
  + 2_TCGA_BLCA
    + [TCGA-BLCA.GDC_phenotype.tsv](2_data/2_TCGA_BLCA/TCGA-BLCA.GDC_phenotype.tsv): phenotype data
    + [TCGA-BLCA.survival.tsv](2_data/2_TCGA_BLCA/TCGA-BLCA.survival.tsv): survival time and status data
  + 3_hallmark_geneset
    + [exp_mat.txt](2_data/3_hallmark_geneset/exp_mat.txt): expression matrix in [.txt format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#TXT:_Text_file_format_for_expression_dataset_.28.2A.txt.29). First column must be gene name (Name), second an description (Description) and then one column per sample (Sample1, Sample2 ... SampleX)
    + [h.all.v7.5.1.symbols.gmt](2_data/3_hallmark_geneset/h.all.v7.5.1.symbols.gmt): gene set that summarize and represent specific well-defined biological states or processes and display coherent expression in [.gmt format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29).
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

### 2. HOW TO
Two different versions of R are required to run the entire pipeline, due to incompatibility of the DEseq package with the most recent versions of the software. The required versions can be installed with the [RSwitch](https://rud.is/rswitch/) program. During the development of this project, R v3.3.3 and R v4.2.0 for in macOS Monterrey 12.3 were used.

All datasets (IMvigor210 and TCGA-BLCA) have been anonymized because of the potential patentability of the results. The details of how the anonymization has been done are not included, but the pipeline would work the same with any similar data type.

#### A. Download Imvigor210 Core Biologies data
To download and install [IMvigor210 Core Biologies package v1.0.0](http://research-pub.gene.com/IMvigor210CoreBiologies/packageVersions/), R 3.3.3 version is required.

```
wget http://research-pub.gene.com/IMvigor210CoreBiologies/packageVersions/IMvigor210CoreBiologies_1.0.0.tar.gz
```

Next step is to install and load IMvigor210CoreBiologies package. If an error appears regarding the installation of a package, it must be installed before trying again.

```
install.packages("IMvigor210CoreBiologies_1.0.0.tar.gz", repos=NULL)
library(IMvigor210CoreBiologies)
```

You can obtain matrix expression through `counts(cds)` object. Feature characterize data its saved in `fData` and phenotype (clinical) data its in `pData`. These three datasets are listed on the github as [exmat_censored_IMvigor210.csv](2_data/1_IMvigor210/exmat_censored_IMvigor210.csv) (expression matrix), [fData_IMvigor210.csv](2_data/1_IMvigor210/fData_IMvigor210.csv) (feature data) and [pData_IMvigor210.csv](2_data/1_IMvigor210/pData_IMvigor210.csv) (phenotype data). The code used is in the script [1_download_data.R](1_scripts/2_pipeline/1_download_data.R).

#### B. Preprocess gene expression data
Data quality analysis and normalization to TPMs is performed. In addition, three filters are applied on the data:
* Remove genes with 0 values in at least 95% of samples
* Drop genes with IQR <= 0.5 and include only genes with log(TPM+1) >= 2 in at least 10% of samples
* Drop genes with s.d <= 1 and include only genes with log(TPM+1) >= 2 in at least 10% of samples

The last two filters are influential. In our analysis we use the second one (SD), however, the code necessary to apply the other one (IQR) is attached. The necessary code is in the script [2_preprocess_data.R](1_scripts/2_pipeline/2_preprocess_data.R).
