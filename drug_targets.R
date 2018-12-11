bnf <- clean_bnf() 
drugbank <- clean_drugbank()
drug_targets <- merge(bnf,drugbank,all.x = TRUE)

# Supplementary table 1 =======================================================

st1 <- unique(drug_targets)

st1$gene <- ifelse(is.na(st1$gene),"",st1$gene)

st1 <- st1 %>% 
  group_by(drug,substance,drugbank_id) %>%
  summarise_all(funs(paste(., collapse=";")))

write.xlsx(data.frame(st1), file="output/Supplementary_Tables_AntihypertensivesMR.xlsx", 
           sheetName="ST1",append=TRUE, row.names = FALSE, showNA = FALSE)

# Match drug substances with targets ==========================================

drug_targets <- drug_targets[!is.na(drug_targets$gene),]

# Restrict to unique drug-gene pairs ==========================================

drug_targets  <- drug_targets[,c("drug","gene")]
drug_targets  <- unique(drug_targets)

# Save output =================================================================

write.csv(drug_targets,file="data/drug_targets.csv",row.names = FALSE)