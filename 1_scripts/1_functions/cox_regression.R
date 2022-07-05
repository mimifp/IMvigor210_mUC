#' Calculate optimal cutoff point and do survival analysis by Cox PH regression
#' model with bonferroni p-value adjust
#' @param coxdata matrix with expression values for genes (columns) and samples 
#' (rows)
#' @return matrix with genes and their p-values for survival analysis

cox_regression <- function (coxdata) {
  genes <- names(coxdata)[4:ncol(coxdata)]
  opt_cut <- surv_cutpoint(coxdata, 
                           time = "os", 
                           event = "censOS", 
                           variables = all_of(genes),
                           progressbar = F)
  
  cox_cat <- surv_categorize(opt_cut)
  
  # Relevel for high vs. low comparision
  cox_cat <- as.data.frame(unclass(cox_cat), stringsAsFactors = TRUE)
  cox_cat[genes] <- lapply(cox_cat[genes], relevel, ref = "low")   
  
  cox_cat <- as_tibble(cox_cat)
  cox_cat <- cox_cat %>% relocate(c(os, censOS), .before = 1)
  genes <- names(cox_cat)[3:ncol(cox_cat)]
  
  # CoxRegression
  results <- RegParallel(
    data = cox_cat,
    formula = 'Surv(os, censOS) ~ [*]',
    FUN = function(formula, data)
      coxph(formula = formula,
            data = data,
            ties = 'breslow',
            singular.ok = TRUE), #skip over columns that are linear combinations of earlier columns
    FUNtype = 'coxph',
    variables = all_of(genes),
    blocksize = 2000, # optimize parallelism
    cores = 4,
    p.adjust = "bonferroni") # change this to BH or Bonferroni to change method 
  
  return(results)
}
