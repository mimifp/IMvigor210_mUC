#' Transform raw counts in trancripts per million (TPMs)
#' @param counts matrix with raw counts for samples as cols and genes as rows
#' @param featureLength matrix with name of gene and their length
#' @return matrix with TPMs 

counts_to_tpm <- function(counts, featureLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- featureLength
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen)
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  tpm <- as.data.frame(tpm)
  return(tpm)
}