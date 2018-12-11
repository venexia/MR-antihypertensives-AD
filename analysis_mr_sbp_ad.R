# List drug targets ===========================================================

drug_targets <- fread("data/drug_targets.csv",
                      stringsAsFactors = FALSE,
                      data.table = FALSE)

# List instruments ============================================================

mr_express <- fread("data/mr_exp_sbp.csv",
                 stringsAsFactors = FALSE,
                 data.table = FALSE)

mr_express <- mr_express[!is.na(mr_express$beta.mr) & sign(mr_express$lci.mr)==sign(mr_express$uci.mr),c("snps.in","gene")]

colnames(mr_express) <- c("SNP","gene")

mr_express <- unique(mr_express)

# Load SBP data ===============================================================

input_sbp <- suppressWarnings(extract_outcome_data(mr_express$SNP,
                                                c('UKB-a:360'),
                                                proxies = FALSE))

input_sbp <- input_sbp[,c("beta.outcome","other_allele.outcome","effect_allele.outcome",
                    "SNP","pval.outcome","se.outcome")]

colnames(input_sbp) <- c("beta.orig","other_allele.orig","effect_allele.orig",
                      "SNP","pvalue","se")

input_sbp <- unique(input_sbp)

# Make all effects positive ---------------------------------------------------

input_sbp$beta <- ifelse(sign(input_sbp$beta.orig)==1,input_sbp$beta.orig,-1*input_sbp$beta.orig)

input_sbp$effect_allele <- ifelse(sign(input_sbp$beta.orig)==1,input_sbp$effect_allele.orig,input_sbp$other_allele.orig)

input_sbp$other_allele <- ifelse(sign(input_sbp$beta.orig)==1,input_sbp$other_allele.orig,input_sbp$effect_allele.orig)

input_sbp[,c("beta.orig","other_allele.orig","effect_allele.orig")] <- list(NULL)

# Merge with drug targets =====================================================

input_sbp <- merge(input_sbp,mr_express,by = c("SNP"))
input_sbp <- merge(input_sbp,drug_targets,by = c("gene"))

# Save ========================================================================

input_sbp <- input_sbp[,c(8,1:2,6:7,5,4,3)]
input_sbp$unit <- "sd"
write.csv(input_sbp,file="data/input_sbp.csv",row.names = FALSE,na = "")

# Target analysis =============================================================

mr_target <- group_MR(df = input_sbp, outcome.data = "mrbase", outcome = c(297), gene)
mr_target$level <- "Target"
names(mr_target)[names(mr_target)=="gene"] <- "leveli"

# Drug class analysis =========================================================

mr_drug <- group_MR(df = input_sbp, outcome.data = "mrbase", outcome = c(297), drug)
mr_drug$level <- "Drug"
names(mr_drug)[names(mr_drug)=="drug"] <- "leveli"

# Combined target analysis ====================================================

input_sbp$grp <- "Combined"
mr_overall <- group_MR(df = input_sbp, outcome.data = "mrbase", outcome = c(297), grp)
mr_overall$level <- "Combined"
names(mr_overall)[names(mr_overall)=="grp"] <- "leveli"

# Overall analysis ============================================================

input_ukb <- extract_instruments(outcomes='UKB-a:360')

input_ukb <- input_ukb[,c("SNP","beta.exposure","se.exposure","samplesize.exposure","pval.exposure",
                          "eaf.exposure","effect_allele.exposure","other_allele.exposure")]

colnames(input_ukb) <- c("SNP","beta","se","samplesize","pval",
                         "eaf","effect_allele","other_allele")

write.csv(input_ukb,file="data/input_ukb.csv",row.names = FALSE)

input_ukb$grp <- "Overall"

mr_ukb <- group_MR(df = input_ukb, outcome.data = "mrbase", outcome = c(297), grp)

mr_ukb$level <- "Overall"

names(mr_ukb)[names(mr_ukb)=="grp"] <- "leveli"

# Make results file ===========================================================

results <- rbind(mr_overall,mr_drug,mr_target,mr_ukb)

# Format results data =========================================================
# Divide by sd of SBP in UK Biobank (19.268) and multiply by 10 to get effect per 10mmHg

results$index <- NULL
results$or <- 1/exp(results$beta.mr*(10/19.268))
results$uci <- 1/exp(results$lci.mr*(10/19.268))
results$lci <- 1/exp(results$uci.mr*(10/19.268))
results$methodtype <- ifelse(results$method=="MR Egger" & !is.na(results$method),"sensitivity","main")

# Save results ================================================================

write.csv(results,file="data/mr_sbp_ad.csv",row.names = FALSE,na = "")