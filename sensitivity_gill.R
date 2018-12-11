# Beta-adrenoceptor blockers from Gill et al ==================================

gill <- read.xlsx("data/gill.xlsx",sheetName = "Sheet1",
                  stringsAsFactors=FALSE)

input_gill <- suppressWarnings(extract_outcome_data(gill$SNP,
                                                    c('UKB-a:360'),
                                                    proxies = FALSE))

input_gill <- input_gill[,c("SNP","beta.outcome","se.outcome","samplesize.outcome",
                            "pval.outcome","eaf.outcome","effect_allele.outcome",
                            "other_allele.outcome")]

input_gill <- merge(input_gill,gill[,c("SNP","Drug")],all.x=TRUE,by=c("SNP"))

colnames(input_gill) <- c("SNP","beta","se","samplesize","pval",
                          "eaf","effect_allele","other_allele","drug")

write.csv(input_gill,file="data/input_gill.csv",row.names = FALSE)

mr_gill <- group_MR(df = input_gill[input_gill$SNP!="rs1801253",], outcome.data = "mrbase", outcome = c(297), drug)

mr_gill$index <- NULL

mr_gill$study <- "Gill et al"
mr_gill$leveli <- mr_gill$drug

mr_gill <- mr_gill[,c("method","beta.mr","lci.mr","uci.mr","pval.mr","study",
                      "leveli","nsnp")]

# WALKER ======================================================================

mr_walker <- read.csv("data/mr_sbp_ad.csv")

mr_walker$study <- "Present study"

mr_walker <- mr_walker[mr_walker$leveli %in% c("Beta-adrenoceptor blockers","Calcium channel blockers"),
                       c("method","beta.mr","lci.mr","uci.mr","pval.mr","study","leveli","nsnp")]

# Combine all studies =========================================================

df <- rbind(mr_gill,mr_walker)

# Create OR labels ============================================================

df$or <- 1/exp(df$beta.mr*(10/19.268))
df$uci <- 1/exp(df$lci.mr*(10/19.268))
df$lci <- 1/exp(df$uci.mr*(10/19.268))

df$orlab <- paste0("OR: ",sprintf("%.2f",df$or),
                   " (95% CI: ", sprintf("%.2f",df$lci),
                   " to ", sprintf("%.2f",df$uci),
                   "); p = ", sprintf("%.2f",df$pval.mr))

df$orlab <- ifelse(df$pval.mr<0.01,
                   paste0("OR: ",sprintf("%.2f",df$or),
                          " (95% CI: ", sprintf("%.2f",df$lci),
                          " to ", sprintf("%.2f",df$uci),
                          "); p = ", sprintf("%.3f",df$pval.mr)),
                   df$orlab)

df$orlab <- ifelse(df$pval.mr<0.001,
                   paste0("OR: ",sprintf("%.2f",df$or),
                          " (95% CI: ", sprintf("%.2f",df$lci),
                          " to ", sprintf("%.2f",df$uci),
                          "); p < 0.001"),
                   df$orlab)

df$orlab <- ifelse(df$pval.mr>0.99,
                   paste0("OR: ",sprintf("%.2f",df$or),
                          " (95% CI: ", sprintf("%.2f",df$lci),
                          " to ", sprintf("%.2f",df$uci),
                          "); p > 0.99"),
                   df$orlab)

df$lab <- paste0(df$study,": ",df$leveli)

# Map main and sensitivity analyses to same label =============================

df$methodtype <- ifelse(df$method=="MR Egger","sensitivity","main")
tmp <- df[,c("lab","methodtype","orlab")]
tmp <- spread(tmp, methodtype, orlab)
tmp$sensitivity <- ifelse(is.na(tmp$sensitivity),"",tmp$sensitivity)
df <- merge(df,tmp, by = c("lab"))

# Define study labels =========================================================

df$studylab <- paste0(df$study," (# SNPs = ",df$nsnp,
                      ")\n",df$main,
                      "\n",df$sensitivity)

df$studylab <- factor(df$studylab)
df$studylab <- factor(df$studylab,levels(df$studylab)[c(3,1,4,2)])

# Refine method label =========================================================

df$methodtype <- factor(df$methodtype)
df$methodtype <- fct_rev(factor(df$methodtype))

# Plot sensitivity results ====================================================

ggplot(df, aes(x = studylab,y = or,group = methodtype, shape = methodtype)) + 
  geom_point(stat = "identity", position = position_dodge(width = 0.5)) + 
  geom_linerange(aes(ymin = lci, ymax = uci), position = position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, linetype = 2) +
  scale_shape_manual(name  ="Method",
                     breaks=c("main", "sensitivity"),
                     labels=c("Wald ratio / IVW", "MR Egger"),
                     values=c(17,15),
                     guide = guide_legend(nrow = 1)) +
  scale_x_discrete(name = "") +
  scale_y_log10(name = "OR and 95% CI for developing Alzheimer's disease for a 10mmHg\ndecrease in systolic blood pressure (presented on log scale)",
                limits = c(min(df$lci),max(df$uci)), breaks = c(0.015625,0.03125,0.0625,0.125,0.25,0.5,1,2,4,8,16,32,64)) +
  theme_minimal() +
  theme(panel.grid.major.y=element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=8),
        text=element_text(size=8),
        legend.position="bottom")  +
  coord_flip() +
  facet_wrap(~leveli, ncol = 1, scales = "free_y")

ggsave("output/mr_gill.jpeg", height = 6, width = 12, unit = "cm", dpi = 600, scale = 1.75)

