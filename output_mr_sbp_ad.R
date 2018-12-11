
  

  # Load drug targets ===========================================================
  
  drug_targets <- fread("data/drug_targets.csv",
                        stringsAsFactors = FALSE,
                        data.table = FALSE)
  
  # Load results ================================================================
  
  results <- data.frame(read_csv("data/mr_sbp_ad.csv"))
  results <- results[,c("leveli","method","nsnp","pval.mr","level" ,"or","uci","lci","methodtype")]
  
  # Remove results with no OR ===================================================
  
  results <- results[!is.na(results$or),]
  
  # Annotate gene and add drug info =============================================
  
  results <- merge(results,drug_targets,by.x=c("leveli"),by.y=c("gene"),all.x = TRUE)
  results$drug <- ifelse(is.na(results$drug),results$leveli,results$drug)
  
  # Create OR labels ============================================================
  
  results$orlab <- paste0("OR: ",
                          ifelse(sprintf("%.2f",results$or)<0.0051,
                                 format(results$or,scientific = TRUE,digits=3),
                                 sprintf("%.2f",results$or)),
                          " (95% CI: ", 
                          ifelse(sprintf("%.2f",results$lci)<0.0051,
                                 format(results$lci,scientific = TRUE,digits=3),
                                 sprintf("%.2f",results$lci)),
                          " to ",
                          ifelse(sprintf("%.2f",results$uci)<0.0051,
                                 format(results$uci,scientific = TRUE,digits=3),
                                 sprintf("%.2f",results$uci)),
                          ")")
  
  results$orlab <- ifelse(results$pval.mr>=0.01 & results$pval.mr<=0.99,
                          paste0(results$orlab,"; p = ", sprintf("%.2f",results$pval.mr),"; #SNPs = ",results$nsnp),
                          results$orlab)
  
  results$orlab <- ifelse(results$pval.mr<0.01,
                          paste0(results$orlab,"; p = ", sprintf("%.3f",results$pval.mr),"; #SNPs = ",results$nsnp),
                          results$orlab)
  
  results$orlab <- ifelse(results$pval.mr<0.001,
                          paste0(results$orlab,"; p < 0.001; #SNPs = ",results$nsnp),
                          results$orlab)
  
  results$orlab <- ifelse(results$pval.mr>0.99,
                          paste0(results$orlab,"; p > 0.99; #SNPs = ",results$nsnp),
                          results$orlab)
  
  # Refine gene labels ==========================================================
  
  results$gene <- results$leveli
  results$gene <- ifelse(results$level=="Drug","All targets",results$gene)
  results$genelab <- ifelse(results$level %in% c("Overall","Combined"),
                            paste0(results$orlab),
                            paste0(results$gene,"; ",results$orlab))
  results$genelab <- fct_rev(factor(results$genelab))
  results$genelab <- factor(results$genelab,levels(results$genelab)[c(12:13,34:35,49:70,1:11,14:33,36:48,71:81)])
  
  # Distinguish individual target results =======================================
  
  results$alltargs <- ifelse(results$level=="Target","No","Yes")
  
  # Refine drug labels ==========================================================
  
  results$drug <- ifelse(results$level!="Target",results$leveli,results$drug)
  results$drug <- ifelse(results$drug=="Overall","Systolic blood pressure",results$drug)
  results$drug <- ifelse(results$drug=="Combined","Antihypertensive drugs",results$drug)
  results$drug <- fct_rev(factor(results$drug))
  results$drug <- factor(results$drug,levels(results$drug)[c(3,10,14:11,9:4,2:1)])
  
  # Plot target effects =========================================================
  
  ggplot(results[results$methodtype=="main",], 
         aes(x = genelab,y = or, col = alltargs)) + 
    geom_point(aes(shape = alltargs)) + 
    geom_linerange(aes(ymin = lci, ymax = uci)) +
    geom_hline(yintercept=1, linetype = 2) +
    scale_color_manual(values = c("#999999","black")) +
    scale_shape_manual(values=c(16,15))+
    scale_x_discrete(name = "", position = "top") +
    scale_y_log10(name = "OR and 95% CI for developing Alzheimer's disease for a 10mmHg\ndecrease in systolic blood pressure (presented on log scale)") +
    theme_minimal() +
    coord_flip() +
    facet_grid(drug~., scales = "free", space = "free", switch = "y") +
    theme(panel.grid.major.y = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text=element_text(size=8),
          text=element_text(size=8),
          strip.text.y = element_text(size=8,hjust = 1, angle = 180),
          legend.position="none")
  
  ggsave("output/mr_sbp_ad_target.jpeg", height = 15, width = 13, unit = "cm", dpi = 600, scale = 1.75)
  
  # Plot drug class effects =========================================================
  
  ggplot(results[results$methodtype=="main" & results$level!="Target",], 
         aes(x = orlab,y = or, col = alltargs)) + 
    geom_point(aes(shape = alltargs)) + 
    geom_linerange(aes(ymin = lci, ymax = uci)) +
    geom_hline(yintercept=1, linetype = 2) +
    scale_color_manual(values = c("black")) +
    scale_x_discrete(name = "", position = "top") +
    scale_y_log10(name = "OR and 95% CI for developing Alzheimer's disease for a 10mmHg\ndecrease in systolic blood pressure (presented on log scale)") +
    theme_minimal() +
    coord_flip() +
    facet_grid(drug~., scales = "free", space = "free", switch = "y") +
    theme(panel.grid.major.y = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text=element_text(size=8),
          text=element_text(size=8),
          strip.text.y = element_text(size=8,hjust = 1, angle = 180),
          legend.position="none")
  
  ggsave("output/mr_sbp_ad_class.jpeg", height = 7, width = 12, unit = "cm", dpi = 600, scale = 1.75)
  
  # Map main and sensitivity analyses to same label =============================
  
  sensitivity <-  results[results$level!="Target",]
  
  tmp <- sensitivity[,c("leveli","methodtype","orlab")]
  tmp <- spread(tmp, methodtype, orlab)
  tmp$sensitivity <- ifelse(is.na(tmp$sensitivity),"",tmp$sensitivity)
  sensitivity <- merge(sensitivity,tmp)
  
  # Refine drug labels in sensitivity dataset ===================================
  
  sensitivity$leveli <- ifelse(sensitivity$leveli=="Overall",
                               "Systolic blood pressure",
                               sensitivity$leveli)
  
  sensitivity$leveli <- ifelse(sensitivity$leveli=="Combined",
                               "Antihypertensive drugs",
                               sensitivity$leveli)
  
  sensitivity$druglab <- paste0(sensitivity$leveli," (# SNPs = ",sensitivity$nsnp,
                             ")\n",sensitivity$main,
                             "\n",sensitivity$sensitivity)
  
  sensitivity$druglab <- fct_rev(factor(sensitivity$druglab))
  sensitivity$druglab <- factor(sensitivity$druglab,levels(sensitivity$druglab)[c(1:2,4:9,11:14,10,3)])
  
  # Refine method label =========================================================
  
  sensitivity$methodtype <- factor(sensitivity$methodtype)
  sensitivity$methodtype <- fct_rev(factor(sensitivity$methodtype))
  
  # Plot sensitivity results ====================================================
  
  ggplot(sensitivity[sensitivity$nsnp>=10,], aes(x = druglab,y = or,group = methodtype, shape = methodtype)) + 
    geom_point(stat = "identity", position = position_dodge(width = 0.5)) + 
    geom_linerange(aes(ymin = lci, ymax = uci), position = position_dodge(width = 0.5)) +
    geom_hline(yintercept=1, linetype = 2) +
    scale_shape_manual(name  ="Method",
                       breaks=c("main", "sensitivity"),
                       labels=c("IVW", "MR-Egger"),
                       values=c(17,15),
                       guide = guide_legend(nrow = 1)) +
    scale_x_discrete(name = "") +
    scale_y_log10(name = "OR and 95% CI for developing Alzheimer's disease for a 10mmHg\ndecrease in systolic blood pressure (presented on log scale)",breaks = c(1e-5,1e-3,1e-1,1,1e1,1e3,1e5)) +
    theme_minimal() +
    theme(panel.grid.major.y=element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text=element_text(size=8),
          text=element_text(size=8),
          strip.text.y = element_blank(),
          legend.position="bottom")  +
    coord_flip()
  
  ggsave("output/mr_sbp_ad_egger.jpeg", height = 7, width = 12, unit = "cm", dpi = 600, scale = 1.75)
  
}