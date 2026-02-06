# Script name: 03_Genotype_stats.R
# ==============================================================================
# Title: Genotype frequency analysis in psychiatric and control cohorts.

# Author: Sergio Pérez Oliveira

# Description: This script performs genotype frequency comparisons across groups 
#              (SCZ, BD, Controls) for genes associated with neurodegenerative
#              disorders: APOE, HTT, ATXN1 (SCA1), and ATXN2 (SCA2). The analysis includes
#              global frequencies, intermediate alleles (IAs), expanded alleles,
#              and subgroup comparisons by clinical subtype and severity (CD).

# Inputs:
#   - Manually selected environment file with custom functions (.R)
#   - Dataframes: DT (full dataset), BD_CONTROLES, SCZ_CONTROLES

# Outputs:
#   - Frequency tables (proportions, row-wise percentages)
#   - Pairwise Fisher's exact tests with Holm correction
# ==============================================================================

# Load environment ----
Env_path <- file.choose()
source(Env_path)
rm(Env_path)

# Function to run genotype frequency analysis ----
run_genotype_analysis <- function(data, group_col, genotype_col, label, drop_cols = NULL) {
  tab <- TABLE(data, group_col, genotype_col, label)
  print(rowPercents(tab))
  my_function(tab)
  
  if (!is.null(drop_cols)) {
    cols <- colnames(tab)
    idx_drop <- if (is.numeric(drop_cols)) intersect(drop_cols, seq_along(cols)) else match(drop_cols, cols, nomatch = 0)
    tab <- tab[, setdiff(seq_along(cols), idx_drop), drop = FALSE]
  }
  
  if (ncol(tab) == 2 && nrow(tab) >= 2) {
    pairwise_fisher_test(tab, p.adjust.method = "holm", conf.int = TRUE, detailed = TRUE)
  } else {
    message("Pairwise Fisher skipped: requires 2 genotype columns and >=2 groups.")
  }
}

# HTT ANALYSIS =================
run_genotype_analysis(DT, "PATHOLOGY", "HTT_CODE", "HTT: Main group comparison")
run_genotype_analysis(DT, "PATHOLOGY", "HTT_CODE", "HTT: Intermediate alleles", drop_cols = 3)
run_genotype_analysis(DT, "PATHOLOGY", "HTT_CODE", "HTT: Expanded alleles", drop_cols = 2)

run_genotype_analysis(BD_CONTROLS, "PATHOLOGY_TYPE_BINARY	", "HTT_CODE", "HTT vs type of BD", drop_cols = 3)
run_genotype_analysis(BD_CONTROLS, "PATHOLOGY_TYPE_BINARY	", "HTT_CODE", "HTT vs type of BD (expanded)", drop_cols = 2)
run_genotype_analysis(BD_CONTROLS, "CD_BINARY", "HTT_CODE", "HTT vs CD severity in BD", drop_cols = 3)

run_genotype_analysis(SCZ_CONTROLS, "PATHOLOGY_TYPE_BINARY	", "HTT_CODE", "HTT vs type of SCZ", drop_cols = 3)
run_genotype_analysis(SCZ_CONTROLS, "CD_BINARY", "HTT_CODE", "HTT vs CD severity in SCZ", drop_cols = 3)

# ATXN1 ANALYSIS =================
run_genotype_analysis(DT, "PATHOLOGY", "ATXN1_CODE", "ATXN1: IA comparison")
run_genotype_analysis(DT, "PATHOLOGY", "ATXN1_CODE", "ATXN1: Intermediate alleles", drop_cols = 3)

run_genotype_analysis(BD_CONTROLS, "PATHOLOGY_TYPE_BINARY	", "ATXN1_CODE", "ATXN1 vs type of BD", drop_cols = 3)
run_genotype_analysis(BD_CONTROLS, "CD_BINARY", "ATXN1_CODE", "ATXN1 vs CD severity in BD", drop_cols = 3)

run_genotype_analysis(SCZ_CONTROLS, "PATHOLOGY_TYPE_BINARY	", "ATXN1_CODE", "ATXN1 vs type of SCZ", drop_cols = 3)
run_genotype_analysis(SCZ_CONTROLS, "CD_BINARY", "ATXN1_CODE", "ATXN1 vs CD severity in SCZ", drop_cols = 3)

# ATXN2 ANALYSIS =================
run_genotype_analysis(DT, "PATHOLOGY", "ATXN2_CODE", "ATXN2: Main group comparison")
run_genotype_analysis(DT, "PATHOLOGY", "ATXN2_CODE", "ATXN2: Intermediate alleles", drop_cols = 3)
run_genotype_analysis(DT, "PATHOLOGY", "ATXN2_CODE", "ATXN2: Expanded alleles", drop_cols = 2)

run_genotype_analysis(BD_CONTROLS, "PATHOLOGY_TYPE_BINARY	", "ATXN2_CODE", "ATXN2 vs type of BD", drop_cols = 3)
run_genotype_analysis(BD_CONTROLS, "PATHOLOGY_TYPE_BINARY	", "ATXN2_CODE", "ATXN2 vs type of BD (expanded)", drop_cols = 2)
run_genotype_analysis(BD_CONTROLS, "CD_BINARY", "ATXN2_CODE", "ATXN2 vs CD severity in BD", drop_cols = 3)

run_genotype_analysis(SCZ_CONTROLS, "PATHOLOGY_TYPE_BINARY	", "ATXN2_CODE", "ATXN2 vs type of SCZ", drop_cols = 3)
run_genotype_analysis(SCZ_CONTROLS, "CD_BINARY", "ATXN2_CODE", "ATXN2 vs CD severity in SCZ", drop_cols = 3)

#APOE & COGNITIVE DECLINE ----
run_genotype_analysis(DT, "PATHOLOGY", "APOE_E4", "APOE E4: main comparison")
run_genotype_analysis(BD_CONTROLS, "PATHOLOGY_TYPE_BINARY	", "APOE_E4", "APOE E4: main comparison")
run_genotype_analysis(BD, "BD_type", "APOE_E4", "APOE E4: main comparison")
run_genotype_analysis(BD_CONTROLS, "CD_BINARY", "APOE_E4", "APOE E4: main comparison")
run_genotype_analysis(SCZ_CONTROLS, "PATHOLOGY_TYPE_BINARY	", "APOE_E4", "APOE E4: main comparison")
run_genotype_analysis(SCZ_CONTROLS, "CD_BINARY", "APOE_E4", "APOE E4: main comparison")

run_genotype_analysis(BD, "CD_BINARY", "APOE_E4", "CD and APOE")
run_genotype_analysis(SCZ, "CD_BINARY", "APOE_E4", "CD and APOE")

# Session info ----
sessionInfo()