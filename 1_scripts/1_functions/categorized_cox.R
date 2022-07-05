#' Categorize gene expression as "high" or "low" according to optimal cutoff point
#' @param coxdata matrix of expression with genes as cols and samples as rows
#' @return matrix with high or low value for each gene in each sample

categorize_cox <- function (coxdata){
  genes <- names(coxdata)[4:ncol(coxdata)]
  opt_cut <- surv_cutpoint(coxdata, 
                           time = "os", 
                           event = "censOS", 
                           variables = all_of(genes),
                           progressbar = F)
  
  cox_cat <- surv_categorize(opt_cut)
  cox_cat <- as_tibble(cox_cat)
  cox_cat <- cox_cat %>% relocate(c(os, censOS), .before = 1)
  
  return(cox_cat)
}