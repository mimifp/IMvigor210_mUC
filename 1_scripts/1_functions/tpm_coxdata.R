#' Prepare tpm matrix for categorization and survival analysis. Add clinical data
#' (survival time and status)
#' @param tpm_matrix expression matrix with TPM normalization. Genes,
#' survival time (os), survival status (censOS) and sample id (sample_id) as 
#' cols. 
#' @reutrn expression matrix with ordered genes and clinical data

tpm_coxdata <- function(tpm_matrix) {
  tpm_matrix <- as.data.frame(t(tpm_matrix))
  tpm_matrix  <- rownames_to_column(tpm_matrix)
  names(tpm_matrix)[1] <- "sample_id"
  
  coxdata <- merge(tpm_matrix, clinical, by = "sample_id")
  coxdata <- coxdata[ , order(c(names(coxdata)))]
  coxdata <- coxdata %>% relocate(c(sample_id, os, censOS), .before = 1)
  
  return(coxdata)
}