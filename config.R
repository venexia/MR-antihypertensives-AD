# Set working directory =======================================================

setwd("")
rm(list=ls())
graphics.off()

# Load libraries ==============================================================

if (!require(pacman)) install.packages("pacman"); library(pacman)
p_load("xlsx")
p_load("data.table")
p_load("furrr")
p_load("dplyr")
p_load("ggplot2")
p_load("tibble")
p_load("tidyr")
p_load("purrr")
p_load("readr")
p_load("forcats")
p_load("biomaRt")
p_load_gh("MRCIEU/TwoSampleMR")
p_load_gh("MRCIEU/MRInstruments")
p_load_gh("slowkow/proxysnps")

# Load functions ==============================================================

source("code/func_MR.R")
source("code/func_group_MR.R")
source("code/clean_drugbank.R")
source("code/clean_bnf.R")

# Setup parrallelization ======================================================

plan(multiprocess)

# Remove supplementary tables file ============================================

if (file.exists("output/Supplementary_Tables_AntihypertensivesMR.xlsx")) {
  file.remove("output/Supplementary_Tables_AntihypertensivesMR.xlsx")
}

# Identify drug targets =======================================================

source("code/drug_targets.R")

# Instrument selection ========================================================

source("code/analysis_mr_exp_sbp.R")

# Main analysis ===============================================================

source("code/analysis_mr_sbp_ad.R")

# Generate supplementary figure 2 =============================================

source("code/output_mr_exp_sbp.R")

# Generate main figure and supplementary figure 3 =============================

source("code/output_mr_sbp_ad.R")

# Generate supplementary tables ===============================================

source("code/supplementary_tables.R")

# Run sensitivity analyses ====================================================

source("code/sensitivity_common_data.R")
source("code/sensitivity_gill.R")
source("code/sensitivity_instrument_overlap.R")
