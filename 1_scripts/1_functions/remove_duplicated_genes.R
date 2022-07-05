remove_duplicated_genes <- function(exmat){
  # Generate a column with the mean of each row after removing gene_id
  gene_means <- rowMeans(exmat %>% dplyr::select(-gene_id))
  
  # Add the column to exmat and put it in second position, after gene_id
  exmat <- exmat %>%
    dplyr::mutate(gene_means) %>%
    dplyr::select(gene_id, gene_means, dplyr::everything())
  
  # Reorder exmat by gene_means in descending order and perform distinct()
  exmat <- exmat %>%
    dplyr::arrange(dplyr::desc(gene_means)) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE)
  
  # Remove gene_means and order genes alphabetically
  exmat <- exmat %>% dplyr::select(-gene_means)
  exmat <- exmat[order(exmat$gene_id), ]
  
  return(exmat)
}