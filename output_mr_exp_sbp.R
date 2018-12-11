# LOAD DATA ===================================================================

drug_targets <- fread("data/drug_targets.csv",
                      stringsAsFactors = FALSE,
                      data.table = FALSE)

mr_exp_sbp <- read_csv("data/mr_exp_sbp.csv")

# FORMAT DATA =================================================================

mr_exp_sbp$or <- 1/exp(mr_exp_sbp$beta.mr*(10/19.268))
mr_exp_sbp$uci <- 1/exp(mr_exp_sbp$lci.mr*(10/19.268))
mr_exp_sbp$lci <- 1/exp(mr_exp_sbp$uci.mr*(10/19.268))
mr_exp_sbp <- merge(mr_exp_sbp,drug_targets)
mr_exp_sbp$ex_null <- ifelse(mr_exp_sbp$lci<1 & mr_exp_sbp$uci>1,NA,1)

# MARK GENES WITHOUT EVIDENCE FOR AN EFFECT ON SBP ============================

mr_exp_sbp$proceed <- ifelse(!is.na(mr_exp_sbp$beta.mr) & sign(mr_exp_sbp$lci.mr)==sign(mr_exp_sbp$uci.mr),1,0)

tmp <- mr_exp_sbp %>% 
  group_by(gene) %>% 
  summarise(include = max(proceed)) %>%
  ungroup()

mr_exp_sbp$gene_di <- ifelse(mr_exp_sbp$gene %in% tmp[tmp$include==0,]$gene,
                          paste0(mr_exp_sbp$gene,"*"),mr_exp_sbp$gene)

# LABEL GTEX TISSUES ===========================================================

gtexlab <- data.frame(tissue = unique(mr_exp_sbp$tissue),
                      tissue_di = c("Adipose - Subcutaneous",
                                    "Adipose - Visceral (Omentum)",
                                    "Adrenal Gland",
                                    "Artery - Aorta",
                                    "Artery - Coronary",
                                    "Artery - Tibial",
                                    "Brain - Amygdala",
                                    "Brain - Anterior cingulate cortex (BA24)",
                                    "Brain - Caudate (basal ganglia)",
                                    "Brain - Cerebellar Hemisphere",
                                    "Brain - Cerebellum",
                                    "Brain - Cortex",
                                    "Brain - Frontal Cortex (BA9)",
                                    "Brain - Hippocampus",
                                    "Brain - Hypothalamus",
                                    "Brain - Nucleus accumbens (basal ganglia)",
                                    "Brain - Putamen (basal ganglia)",
                                    "Brain - Spinal cord (cervical c-1)",
                                    "Brain - Substantia nigra",
                                    "Breast - Mammary Tissue",
                                    "Cells - EBV-transformed lymphocytes",
                                    "Cells - Transformed fibroblasts",
                                    "Colon - Sigmoid",
                                    "Colon - Transverse",
                                    "Esophagus - Gastroesophageal Junction",
                                    "Esophagus - Mucosa",
                                    "Esophagus - Muscularis",
                                    "Heart - Atrial Appendage",
                                    "Heart - Left Ventricle",
                                    "Liver",
                                    "Lung",
                                    "Minor Salivary Gland",
                                    "Muscle - Skeletal",
                                    "Nerve - Tibial",
                                    "Ovary",
                                    "Pancreas",
                                    "Pituitary",
                                    "Prostate",
                                    "Skin - Not Sun Exposed (Suprapubic)",
                                    "Skin - Sun Exposed (Lower leg)",
                                    "Small Intestine - Terminal Ileum",
                                    "Spleen",
                                    "Stomach",
                                    "Testis",
                                    "Thyroid",
                                    "Uterus",
                                    "Vagina",
                                    "Whole Blood"),
                      stringsAsFactors = FALSE)

mr_exp_sbp <- merge(mr_exp_sbp,gtexlab)

ggplot(mr_exp_sbp, aes(y = reorder(gene_di, desc(gene_di)), x = tissue_di)) +
  labs(y = "Target", x = "Tissue", fill = "OR\n") +
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
  geom_tile(aes(fill = or)) +
  geom_point(aes(shape = factor(ex_null)), size = 0.5) +
  scale_shape_discrete(name  = " ",
                       breaks=c("1", "0"),
                       labels=c("95% CI excludes null", "Label")) +
  scale_fill_distiller(palette = "Spectral", direction = -1,
                       guide = "colourbar", limits = c(0.9,1.1),
                       na.value = NA)
ggsave("output/mr_exp_sbp.jpeg",width = 200, height = 300, unit = "mm", dpi = 600)
