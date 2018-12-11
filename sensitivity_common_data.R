# LARSSON ===========================================================

larsson <- read.xlsx("data/larsson.xlsx",sheetName = "Sheet1",
                     stringsAsFactors=FALSE)

larsson$snp <- gsub(" ","",larsson$snp)

input_lar <- suppressWarnings(extract_outcome_data(larsson$snp,
                                             c('UKB-a:360'),
                                             proxies = FALSE))

input_lar <- input_lar[,c("SNP","beta.outcome","se.outcome","samplesize.outcome","pval.outcome",
                      "eaf.outcome","effect_allele.outcome","other_allele.outcome")]

colnames(input_lar) <- c("SNP","beta","se","samplesize","pval",
                      "eaf","effect_allele","other_allele")

write.csv(input_lar,file="data/input_larsson.csv",row.names = FALSE)

input_lar$study <- "Larsson et al"

mr_larsson <- group_MR(df = input_lar, outcome.data = "mrbase", outcome = c(297), study)

mr_larsson$index <- NULL

# OSTERGAARD ========================================================

ostergaard <- read.xlsx("data/ostergaard.xlsx",sheetName = "Sheet1",
                        stringsAsFactors=FALSE)

input_ost <- suppressWarnings(extract_outcome_data(ostergaard$snp,
                                                 c('UKB-a:360'),
                                                 proxies = FALSE))

input_ost <- input_ost[,c("SNP","beta.outcome","se.outcome","samplesize.outcome","pval.outcome",
                      "eaf.outcome","effect_allele.outcome","other_allele.outcome")]

colnames(input_ost) <- c("SNP","beta","se","samplesize","pval",
                       "eaf","effect_allele","other_allele")

write.csv(input_ost,file="data/input_ostergaard.csv",row.names = FALSE)

input_ost$study <- "Ostergaard et al"

mr_ostergaard <- group_MR(df = input_ost, outcome.data = "mrbase", outcome = c(297), study)

mr_ostergaard$index <- NULL

# WALKER ======================================================================

mr_walker <- read.csv("data/mr_sbp_ad.csv",stringsAsFactors = FALSE)

mr_walker$study <- NA

mr_walker$study <- ifelse(mr_walker$level=="Overall","Present study: systolic blood pressure",mr_walker$study)

mr_walker$study <- ifelse(mr_walker$level=="Combined","Present study: antihypertensive drugs",mr_walker$study)

mr_walker <- mr_walker[mr_walker$level %in% c("Overall","Combined"),colnames(mr_larsson)]


# Combine all studies =========================================================

df <- rbind(mr_larsson,mr_ostergaard,mr_walker)

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

# Map main and sensitivity analyses to same label =============================

df$methodtype <- ifelse(df$method=="MR Egger","sensitivity","main")
tmp <- df[,c("study","methodtype","orlab")]
tmp <- spread(tmp, methodtype, orlab)
tmp$sensitivity <- ifelse(is.na(tmp$sensitivity),"",tmp$sensitivity)
df <- merge(df,tmp, by = c("study"))

# Define study labels =========================================================

df$studylab <- paste0(df$study," (# SNPs = ",df$nsnp,
                              ")\n",df$main,
                              "\n",df$sensitivity)

df$studylab <- factor(df$studylab)

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
                limits = c(0.5,8), breaks = c(0.5,1,2,4,8)) +
  theme_minimal() +
  theme(panel.grid.major.y=element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=8),
        text=element_text(size=8),
        legend.position="bottom")  +
  coord_flip()

ggsave("output/mr_common_data.jpeg", height = 6, width = 12, unit = "cm", dpi = 600, scale = 1.75)

