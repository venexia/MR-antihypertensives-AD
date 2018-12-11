# List drug targets ===========================================================

drug_targets <- fread("data/drug_targets.csv",
                      stringsAsFactors = FALSE,
                      data.table = FALSE)

# Identify best SNPs for targets ==============================================

gtex_path <- "data/gtex/GTEx_Analysis_v7_eQTL/"
gtex_data <- ".v7.egenes.txt"
tissues <- list.files(path = gtex_path,pattern = paste0("*",gtex_data))

input_exp <- data.frame(gene = rep(unique(drug_targets$gene),each = length(tissues)),
                        tissue = rep(tissues,times = length(unique(drug_targets$gene))),
                        stringsAsFactors = FALSE)

tmp_input_exp <- NULL

for (i in tissues) {
  
  tmp <- suppressWarnings(fread(paste0(gtex_path,i),
                                stringsAsFactors = FALSE,
                                data.table = FALSE))
  
  tmp <- tmp[tmp$gene_name %in% unique(drug_targets$gene),]
  
  tmp$tissue <- gsub(gtex_data,"",i)
  
  tmp_input_exp <- rbind(tmp_input_exp,tmp)
  
}

# Format data for MR analysis =================================================

tmp_input_exp <- tmp_input_exp[,c("slope","ref","alt","gene_name","tissue","rs_id_dbSNP147_GRCh37p13",
                    "gene_chr","gene_start","pval_nominal","slope_se")]

colnames(tmp_input_exp) <- c("beta.orig","other_allele.orig","effect_allele.orig","gene","tissue",
                      "SNP","chr_id","chr_pos","pvalue","se")

# Merge with full gene-tissue df ==============================================

input_exp$tissue <- gsub(".v7.egenes.txt","",input_exp$tissue)
input_exp <- merge(input_exp,tmp_input_exp,all.x = TRUE,by = c("gene","tissue"))

# Make all effects positive ===================================================

input_exp$beta <- ifelse(sign(input_exp$beta.orig)==1,input_exp$beta.orig,-1*input_exp$beta.orig)

input_exp$effect_allele <- ifelse(sign(input_exp$beta.orig)==1,input_exp$effect_allele.orig,input_exp$other_allele.orig)

input_exp$other_allele <- ifelse(sign(input_exp$beta.orig)==1,input_exp$other_allele.orig,input_exp$effect_allele.orig)

input_exp[,c("beta.orig","other_allele.orig","effect_allele.orig")] <- list(NULL)

# Remove unsuitable SNPs ======================================================

input_exp$keep <- ifelse(!is.na(input_exp$SNP) & input_exp$SNP!="." & nchar(input_exp$effect_allele)==1 & nchar(input_exp$other_allele)==1,TRUE,FALSE)

# Save MR input ===============================================================

write.csv(input_exp,file="data/input_exp.csv",row.names = FALSE,na = "")

# Conduct MR of expression on SBP =============================================

input_exp <- input_exp[input_exp$keep==TRUE,]

mr_express <- NULL

for (i in 1:nrow(input_exp)) {
  print(i)
  tmp_input_exp <- input_exp[i,]
  tmp_mr_express <- MR(tmp_input_exp,c('UKB-a:360'))
  tmp_mr_express$gene <- input_exp[i,]$gene
  tmp_mr_express$tissue <- input_exp[i,]$tissue
  mr_express <- rbind(mr_express,tmp_mr_express)
}

mr_express[,c("beta.exposure","se.exposure","pval.exposure")] <- NULL
input_exp <- input_exp[,c("gene","tissue","chr_id","chr_pos","beta","se","pvalue","effect_allele","other_allele")]
colnames(input_exp) <- c("gene","tissue","chr_id","chr_pos","beta.exposure","se.exposure","pval.exposure","effect_allele","other_allele")
mr_express <- merge(mr_express,input_exp,by=c("gene","tissue"),all.x = TRUE)

# mr_express <- group_MR(df = input_exp, outcome.data = "mrbase", outcome = c('UKB-a:360'), tissue, gene)
# mr_express$index <- NULL

# Save output =================================================================

write.csv(mr_express,file="data/mr_exp_sbp.csv",row.names = FALSE,na = "")