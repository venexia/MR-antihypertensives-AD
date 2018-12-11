# Required libraries: TwoSampleMR, tidyverse, furrr

MR <- function(input,outcome,wait=FALSE) {
  
  # System waiting for parallel processing ====================================
  
  if (wait) {
    wait_time <- runif(1, 0.5, 2) 
    Sys.sleep(wait_time)
  }
  
  # Clear previous results ====================================================
  
  exp <- NULL
  out <- NULL
  dat <- NULL
  res <- NULL
  ptest <- NULL
  output <- NULL
  
  # Format exposure data ======================================================
  
  exp <- suppressWarnings(format_data(input, type="exposure"))
  
  # Clump data if required ====================================================
  
  if (nrow(exp)>1) {

    exp <- clump_data(exp)

  }

  # Extract outcome data ======================================================
  
  out <- suppressWarnings(extract_outcome_data(exp$SNP,
                                               outcome,
                                               proxies = 1,
                                               rsq = 0.8,
                                               align_alleles = 1,
                                               palindromes = 1,
                                               maf_threshold = 0.3))
  
  # Harmonise exposure and outcome data =======================================
  
  if (!is.null(out)) {
    
    dat <- harmonise_data(exposure_dat = exp,
                          outcome_dat = out)
    
    if(!is.null(dat)) {
      
      # Perform MR ================================================================
      
      if (nrow(dat)==1 && dat$mr_keep==TRUE) {
        
        res <- mr(dat, method_list = c("mr_wald_ratio"))
        
      } 
      
      if (nrow(dat)>1) {
        
        res <- mr(dat, method_list = c("mr_ivw","mr_egger_regression"))
        ptest <- mr_pleiotropy_test(dat)
        colnames(ptest) <- c("id.exposure","id.outcome","outcome","exposure",
                             "egger_intercept","egger_se","egger_pval")
        res <- merge(res,ptest)
        
      }
    }
  }
  
  # Save results ==============================================================
  
  if (!is.null(res)) {
    if (nrow(res)>0 & is.null(ptest)) {
      output <- res
      output$snps.mr <- paste(dat$SNP,collapse=";")
      output$beta.exposure <- paste(dat$beta.exposure,collapse=";")
      output$se.exposure <- paste(dat$se.exposure,collapse=";")
      output$pval.exposure <- paste(dat$pval.exposure,collapse=";")
      output$beta.outcome <- paste(dat$beta.outcome,collapse=";")
      output$se.outcome <- paste(dat$se.outcome,collapse=";")
      output$pval.outcome <- paste(dat$pval.outcome,collapse=";")
      output$egger_intercept <- NA
      output$egger_se <- NA
      output$egger_pval <- NA
    } 
    else if (nrow(res)>0 & !is.null(ptest)) {
      output <- res
      output$snps.mr <- paste(dat$SNP,collapse=";")
      output$beta.exposure <- paste(dat$beta.exposure,collapse=";")
      output$se.exposure <- paste(dat$se.exposure,collapse=";")
      output$pval.exposure <- paste(dat$pval.exposure,collapse=";")
      output$beta.outcome <- paste(dat$beta.outcome,collapse=";")
      output$se.outcome <- paste(dat$se.outcome,collapse=";")
      output$pval.outcome <- paste(dat$pval.outcome,collapse=";")
    } 
    else {
      output <- setNames(data.frame(matrix(ncol = 20, nrow = 1)),
                        c("id.exposure","id.outcome","outcome","exposure",
                          "method","nsnp","b","se","pval","snps.mr",
                          "egger_intercept","egger_se","egger_pval","snps.in",
                          "beta.exposure","se.exposure","pval.exposure",
                          "beta.outcome","se.outcome","pval.outcome"))
    } 
  } else {
    output <- setNames(data.frame(matrix(ncol = 20, nrow = 1)),
                      c("id.exposure","id.outcome","outcome","exposure",
                        "method","nsnp","b","se","pval","snps.mr",
                        "egger_intercept","egger_se","egger_pval","snps.in",
                        "beta.exposure","se.exposure","pval.exposure",
                        "beta.outcome","se.outcome","pval.outcome"))
  }
  
  output$snps.in <- paste(input$SNP,collapse=";")
  
  output$lci.mr <- output$b - qnorm(0.975)*output$se
  output$uci.mr <- output$b + qnorm(0.975)*output$se
  
  output <- output[,c("method","nsnp","b","lci.mr","uci.mr","se","pval",
                      "snps.mr","egger_intercept","egger_se","egger_pval","snps.in",
                      "beta.exposure","se.exposure","pval.exposure",
                      "beta.outcome","se.outcome","pval.outcome")]
  
  colnames(output) <- c("method","nsnp","beta.mr","lci.mr","uci.mr","se.mr","pval.mr",
                        "snps.mr","egger_intercept","egger_se","egger_pval","snps.in",
                        "beta.exposure","se.exposure","pval.exposure",
                        "beta.outcome","se.outcome","pval.outcome")
  
  return(output)
  
}
