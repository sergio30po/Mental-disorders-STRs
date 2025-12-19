# Script name: 03_Genotype_stats.R
# ==============================================================================
# Title: Genotype frequency analysis in psychiatric and control cohorts.

# Author: Sergio Pérez Oliveira

# Description: This script performs genotype frequency comparisons across groups 
#              (SCZ, BD, Controls) for genes associated with neurodegenerative
#              disorders: HTT, ATXN1 (SCA1), and ATXN2 (SCA2). The analysis includes
#              global frequencies, intermediate alleles (IAs), expanded alleles,
#              and subgroup comparisons by clinical subtype and severity (DCO).

# Inputs:
#   - Manually selected environment file with custom functions (.R)
#   - Dataframes: DT (full dataset), BD_CONTROLES, SCH_CONTROLES

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
  res_misc <- my_function(tab)
  fisher_res <- NULL
  
  if (is.null(drop_cols)) {
    message("Pairwise Fisher skipped: requires 'drop_cols' with 2 categories.")
  } else {
    cols <- colnames(tab)
    idx_drop <- if (is.numeric(drop_cols)) intersect(drop_cols, seq_along(cols)) else match(drop_cols, cols, nomatch = 0)
    kept <- setdiff(seq_along(cols), idx_drop)
    mat <- as.matrix(tab[, kept, drop = FALSE])
    if (ncol(mat) == 2) {
      pairwise_fisher_test(mat, p.adjust.method = "holm", conf.int = TRUE, detailed = TRUE)
    } else {
      message("Pairwise Fisher skipped: reduced table not equal to 2 categories.")
    }
  }
}

# HTT ANALYSIS =================
run_genotype_analysis(DT, "PATHOLOGY", "HTT_CODE", "HTT: Main group comparison")
run_genotype_analysis(DT, "PATHOLOGY", "HTT_CODE", "HTT: Intermediate alleles", drop_cols = 3)
run_genotype_analysis(DT, "PATHOLOGY", "HTT_CODE", "HTT: Expanded alleles", drop_cols = 2)

run_genotype_analysis(BD_CONTROLS, "PATHOLOGY_TYPE_BINARY	", "HTT_CODE", "HTT vs type of BD", drop_cols = 3)
run_genotype_analysis(BD_CONTROLS, "PATHOLOGY_TYPE_BINARY	", "HTT_CODE", "HTT vs type of BD (expanded)", drop_cols = 2)
run_genotype_analysis(BD_CONTROLS, "CD_BINARY", "HTT_CODE", "HTT vs DCO severity in BD", drop_cols = 3)

run_genotype_analysis(SCH_CONTROLS, "PATHOLOGY_TYPE_BINARY	", "HTT_CODE", "HTT vs type of SCH", drop_cols = 3)
run_genotype_analysis(SCH_CONTROLS, "CD_BINARY", "HTT_CODE", "HTT vs DCO severity in SCH", drop_cols = 3)

# ATXN1 (SCA1) ANALYSIS =================
run_genotype_analysis(DT, "PATHOLOGY", "SCA1_CODE", "ATXN1: IA comparison")
run_genotype_analysis(DT, "PATHOLOGY", "SCA1_CODE", "ATXN1: Intermediate alleles", drop_cols = 3)

run_genotype_analysis(BD_CONTROLS, "PATHOLOGY_TYPE_BINARY	", "SCA1_CODE", "SCA1 vs type of BD", drop_cols = 3)
run_genotype_analysis(BD_CONTROLS, "CD_BINARY", "SCA1_CODE", "SCA1 vs DCO severity in BD", drop_cols = 3)

run_genotype_analysis(SCH_CONTROLS, "PATHOLOGY_TYPE_BINARY	", "SCA1_CODE", "SCA1 vs type of SCH", drop_cols = 3)
run_genotype_analysis(SCH_CONTROLS, "CD_BINARY", "SCA1_CODE", "SCA1 vs DCO severity in SCH", drop_cols = 3)

# ATXN2 (SCA2) ANALYSIS =================
run_genotype_analysis(DT, "PATHOLOGY", "SCA2_CODE", "ATXN2: Main group comparison")
run_genotype_analysis(DT, "PATHOLOGY", "SCA2_CODE", "ATXN2: Intermediate alleles", drop_cols = 3)
run_genotype_analysis(DT, "PATHOLOGY", "SCA2_CODE", "ATXN2: Expanded alleles", drop_cols = 2)

run_genotype_analysis(BD_CONTROLS, "PATHOLOGY_TYPE_BINARY	", "SCA2_CODE", "SCA2 vs type of BD", drop_cols = 3)
run_genotype_analysis(BD_CONTROLS, "PATHOLOGY_TYPE_BINARY	", "SCA2_CODE", "SCA2 vs type of BD (expanded)", drop_cols = 2)
run_genotype_analysis(BD_CONTROLS, "CD_BINARY", "SCA2_CODE", "SCA2 vs DCO severity in BD", drop_cols = 3)

run_genotype_analysis(SCH_CONTROLS, "PATHOLOGY_TYPE_BINARY	", "SCA2_CODE", "SCA2 vs type of SCH", drop_cols = 3)
run_genotype_analysis(SCH_CONTROLS, "CD_BINARY", "SCA2_CODE", "SCA2 vs DCO severity in SCH", drop_cols = 3)

# Session info ----
sessionInfo()