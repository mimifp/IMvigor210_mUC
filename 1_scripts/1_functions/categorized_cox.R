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