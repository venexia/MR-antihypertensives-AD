# LOAD DATA ===================================================================

drug_targets <- fread("data/drug_targets.csv",
                      stringsAsFactors = FALSE,
                      data.table = FALSE)

mr_exp_sbp <- read_csv("data/mr_exp_sbp.csv")

# LABEL GTEX TISSUES ===========================================================

gtexlab <- data.frame(rbind(c("Adipose_Subcutaneous","Adipose - Subcutaneous"),
                 c("Adipose_Visceral_Omentum","Adipose - Visceral (Omentum)"),
                 c("Adrenal_Gland","Adrenal Gland"),
                 c("Artery_Aorta","Artery - Aorta"),
                 c("Artery_Coronary","Artery - Coronary"),
                 c("Artery_Tibial","Artery - Tibial"),
                 c("Brain_Amygdala","Brain - Amygdala"),
                 c("Brain_Anterior_cingulate_cortex_BA24","Brain - Anterior cingulate cortex (BA24)"),
                 c("Brain_Caudate_basal_ganglia","Brain - Caudate (basal ganglia)"),
                 c("Brain_Cerebellar_Hemisphere","Brain - Cerebellar Hemisphere"),
                 c("Brain_Cerebellum","Brain - Cerebellum"),
                 c("Brain_Cortex","Brain - Cortex"),
                 c("Brain_Frontal_Cortex_BA9","Brain - Frontal Cortex (BA9)"),
                 c("Brain_Hippocampus","Brain - Hippocampus"),
                 c("Brain_Hypothalamus","Brain - Hypothalamus"),
                 c("Brain_Nucleus_accumbens_basal_ganglia","Brain - Nucleus accumbens (basal ganglia)"),
                 c("Brain_Putamen_basal_ganglia","Brain - Putamen (basal ganglia)"),
                 c("Brain_Spinal_cord_cervical_c-1","Brain - Spinal cord (cervical c-1)"),
                 c("Brain_Substantia_nigra","Brain - Substantia nigra"),
                 c("Breast_Mammary_Tissue","Breast - Mammary Tissue"),
                 c("Cells_EBV-transformed_lymphocytes","Cells - EBV-transformed lymphocytes"),
                 c("Cells_Transformed_fibroblasts","Cells - Transformed fibroblasts"),
                 c("Colon_Sigmoid","Colon - Sigmoid"),
                 c("Colon_Transverse","Colon - Transverse"),
                 c("Esophagus_Gastroesophageal_Junction","Esophagus - Gastroesophageal Junction"),
                 c("Esophagus_Mucosa","Esophagus - Mucosa"),
                 c("Esophagus_Muscularis","Esophagus - Muscularis"),
                 c("Heart_Atrial_Appendage","Heart - Atrial Appendage"),
                 c("Heart_Left_Ventricle","Heart - Left Ventricle"),
                 c("Liver","Liver"),
                 c("Lung","Lung"),
                 c("Minor_Salivary_Gland","Minor Salivary Gland"),
                 c("Muscle_Skeletal","Muscle - Skeletal"),
                 c("Nerve_Tibial","Nerve - Tibial"),
                 c("Ovary","Ovary"),
                 c("Pancreas","Pancreas"),
                 c("Pituitary","Pituitary"),
                 c("Prostate","Prostate"),
                 c("Skin_Not_Sun_Exposed_Suprapubic","Skin - Not Sun Exposed (Suprapubic)"),
                 c("Skin_Sun_Exposed_Lower_leg","Skin - Sun Exposed (Lower leg)"),
                 c("Small_Intestine_Terminal_Ileum","Small Intestine - Terminal Ileum"),
                 c("Spleen","Spleen"),
                 c("Stomach","Stomach"),
                 c("Testis","Testis"),
                 c("Thyroid","Thyroid"),
                 c("Uterus","Uterus"),
                 c("Vagina","Vagina"),
                 c("Whole_Blood","Whole Blood")))

colnames(gtexlab) <- c("tissue","tissue_di")

mr_exp_sbp <- merge(mr_exp_sbp,gtexlab)

# FORMAT DATA =================================================================

mr_exp_sbp <- merge(mr_exp_sbp,drug_targets)
mr_exp_sbp$ex_null <- ifelse(mr_exp_sbp$lci.mr<0 & mr_exp_sbp$uci.mr>0,NA,1)

# MARK GENES WITHOUT EVIDENCE FOR AN EFFECT ON SBP ============================

mr_exp_sbp$proceed <- ifelse(!is.na(mr_exp_sbp$beta.mr) & sign(mr_exp_sbp$lci.mr)==sign(mr_exp_sbp$uci.mr),1,0)

tmp <- mr_exp_sbp %>% 
  group_by(gene) %>% 
  summarise(include = max(proceed)) %>%
  ungroup()

mr_exp_sbp$gene_di <- ifelse(mr_exp_sbp$gene %in% tmp[tmp$include==0,]$gene,
                             paste0(mr_exp_sbp$gene,"*"),mr_exp_sbp$gene)

ggplot(mr_exp_sbp, aes(y = reorder(gene_di, desc(gene_di)), x = tissue_di)) +
  labs(y = "Target", x = "Tissue", fill = "Beta\n") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  facet_grid(mr_exp_sbp$drug~., scales = "free", space = "free", switch = "both") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        legend.position = "bottom",
        legend.text=element_text(size=8),
        legend.title=element_text(size=8),
        axis.text.x = element_text(angle = 90, hjust = 0),
        strip.text.y = element_text(angle = 180, hjust = 1, size=8),
        strip.background = element_blank()) +
  geom_tile(aes(fill = beta.mr)) +
  geom_point(aes(shape = factor(ex_null)), size = 0.5) +
  scale_shape_discrete(name  = " ",
                       breaks=c("1", "0"),
                       labels=c("95% CI excludes null", "Label")) +
  scale_fill_distiller(palette = "Spectral", direction = -1,
                       guide = "colourbar", limits = c(-0.2,0.2),
                       na.value = NA)
ggsave("output/mr_exp_sbp.jpeg",width = 200, height = 300, unit = "mm", dpi = 600)
