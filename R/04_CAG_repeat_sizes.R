# Script name: 04_CAG_repeat_sizes.R
# ==============================================================================
# Title: Comparison of CAG repeat sizes across psychiatric and control groups
# Author: Sergio Pérez Oliveira
# Description: This script performs non-parametric statistical analyses to compare
#              CAG repeat lengths of both alleles in HTT, ATXN1 (SCA1), and ATXN2 (SCA2)
#              across diagnostic groups (SCZ, BD, Controls), as well as in clinically
#              defined subgroups based on subtype and severity (DCO).
#              The analysis includes:
#                - Descriptive statistics (mean and SD per group)
#                - Kruskal-Wallis test for global group differences
#                - Dunn and Wilcoxon pairwise post-hoc tests if significant
#                - Effect size estimation (Wilcoxon-based) for selected comparisons
# Inputs:
#   - Manually selected environment file with custom functions (.R)
#   - Dataframes: DT (full dataset), BP_CONTROLES, SCH_CONTROLES
# Outputs:
#   - Printed summaries of statistics, test results and effect sizes
# ==============================================================================

# Load environment ----
Env_path <- file.choose()
source(Env_path)
rm(Env_path)

run_kruskal_analysis <- function(df, group_col, value_col, result_name = NULL, comparisons_list = NULL) {
  # If result_name is not provided, use the name of value_col
  if (is.null(result_name)) {
    result_name <- value_col
  }
  
  cat("\n--- Analysis for:", result_name, "---\n")
  
  # 1. Mean and standard deviation
  cat("Mean and SD by group:\n")
  mean_sd_result <- mean_sd(df, group_col, value_col, result_name = value_col)
  print(mean_sd_result)
  
  # 2. Kruskal-Wallis test
  cat("\nKruskal-Wallis test:\n")
  kruskal <- kruskal.test(df[[value_col]] ~ df[[group_col]])
  print(kruskal)
  
  # 3. If significant, perform post-hoc tests and effect size
  if (kruskal$p.value < 0.05) {
    cat("\nDunn's test (Holm correction):\n")
    dunn_result <- FSA::dunnTest(df[[value_col]] ~ df[[group_col]], method = "holm")
    print(dunn_result)
    
    cat("\nPairwise Wilcoxon test (Holm correction):\n")
    wilcox_result <- pairwise.wilcox.test(df[[value_col]], df[[group_col]], p.adjust.method = "holm")
    print(wilcox_result)
    
    # 4. Effect size using wilcox_effsize if comparisons provided
    if (!is.null(comparisons_list)) {
      cat("\nEffect size for selected comparisons:\n")
      effsize <- df %>%
        wilcox_effsize(as.formula(paste(value_col, "~", group_col)),
                       paired = FALSE,
                       comparisons = comparisons_list,
                       p.adjust.method = "holm")
      print(effsize)
    } else {
      cat("No specific comparisons provided for effect size.\n")
    }
  } else {
    cat("\nKruskal-Wallis not significant; post-hoc tests and effect sizes are not performed.\n")
  }
}


# HTT ANALYSIS ----

### Long allele
run_kruskal_analysis(DT, "PATHOLOGY", "ALLELE2_HTT", "HTT long allele", list(c("BD","SCH"), c("BD","CONTROL"), c("SCH","CONTROL")))
run_kruskal_analysis(BD_CONTROLS, "CD_BINARY", "ALLELE2_HTT", "HTT long allele BD_CONTROLS - Severity", list(c("CD","NO"), c("CD","CONTROL"), c("NO","CONTROL")))
run_kruskal_analysis(BD_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE2_HTT", "HTT long allele BD_CONTROLS - Pathology type", list(c("TBD-1","Other"), c("TBD-1","CONTROL"), c("Other","CONTROL")))
run_kruskal_analysis(SCH_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE2_HTT", "HTT long allele SCH_CONTROLS - Pathology type", list(c("SCH","Other"), c("Other","CONTROL"), c("SCH","CONTROL")))

### Short allele
run_kruskal_analysis(DT, "PATHOLOGY", "ALLELE1_HTT", "HTT short allele", list(c("BD","SCH"), c("BD","CONTROL"), c("SCH","CONTROL")))
run_kruskal_analysis(BD_CONTROLS, "CD_BINARY", "ALLELE1_HTT", "HTT short allele BD_CONTROLES - Severity", list(c("CD","NO"), c("CD","CONTROL"), c("NO","CONTROL")))
run_kruskal_analysis(BD_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE1_HTT", "HTT short allele BD_CONTROLES - Pathology type", list(c("TBD-1","Other"), c("TBD-1","CONTROL"), c("Other","CONTROL")))
run_kruskal_analysis(SCH_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE1_HTT", "HTT short allele SCH_CONTROLS - Pathology type", list(c("SCH","Other"), c("Other","CONTROL"), c("SCH","CONTROL")))

# ATXN1 (SCA1) ANALYSIS ----

### Long allele
run_kruskal_analysis(DT, "PATHOLOGY", "ALLELE2_SCA1", "SCA1 long allele", list(c("BD","SCH"), c("BD","CONTROL"), c("SCH","CONTROL")))
run_kruskal_analysis(BD_CONTROLS, "CD_BINARY", "ALLELE2_SCA1", "SCA1 long allele BD_CONTROLES - Severity", list(c("CD","NO"), c("CD","CONTROL"), c("NO","CONTROL")))
run_kruskal_analysis(BD_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE2_SCA1", "SCA1 long allele BD_CONTROLES - Pathology type", list(c("TBD-1","Other"), c("TBD-1","CONTROL"), c("Other","CONTROL")))
run_kruskal_analysis(SCH_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE2_SCA1", "SCA1 long allele SCH_CONTROLS - Pathology type", list(c("SCH","Other"), c("Other","CONTROL"), c("SCH","CONTROL")))

### Short allele
run_kruskal_analysis(DT, "PATHOLOGY", "ALLELE1_SCA1", "SCA1 short allele", list(c("BD","SCH"), c("BD","CONTROL"), c("SCH","CONTROL")))
run_kruskal_analysis(BD_CONTROLS, "CD_BINARY", "ALLELE1_SCA1", "SCA1 short allele BD_CONTROLES - Severity", list(c("CD","NO"), c("CD","CONTROL"), c("NO","CONTROL")))
run_kruskal_analysis(BD_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE1_SCA1", "SCA1 short allele BD_CONTROLES - Pathology type", list(c("TBD-1","Other"), c("TBD-1","CONTROL"), c("Other","CONTROL")))
run_kruskal_analysis(SCH_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE1_SCA1", "SCA1 short allele SCH_CONTROLS - Pathology type", list(c("SCH","Other"), c("Other","CONTROL"), c("SCH","CONTROL")))

# ATXN2 (SCA2) ANALYSIS ----

### Long allele
run_kruskal_analysis(DT, "PATHOLOGY", "ALLELE2_SCA2", "SCA2 long allele", list(c("BD","SCH"), c("BD","CONTROL"), c("SCH","CONTROL")))
run_kruskal_analysis(BD_CONTROLS, "CD_BINARY", "ALLELE2_SCA2", "SCA2 long allele BD_CONTROLES - Severity", list(c("CD","NO"), c("CD","CONTROL"), c("NO","CONTROL")))
run_kruskal_analysis(BD_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE2_SCA2", "SCA2 long allele BD_CONTROLES - Pathology type", list(c("BD-I","Other"), c("BD-I","CONTROL"), c("Other","CONTROL")))
run_kruskal_analysis(SCH_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE2_SCA2", "SCA2 long allele SCH_CONTROLS - Pathology type", list(c("SCH","Other"), c("Other","CONTROL"), c("SCH","CONTROL")))

### Short allele
run_kruskal_analysis(DT, "PATHOLOGY", "ALLELE1_SCA2", "SCA2 short allele", list(c("BD","SCH"), c("BD","CONTROL"), c("SCH","CONTROL")))
run_kruskal_analysis(BD_CONTROLS, "CD_BINARY", "ALLELE1_SCA2", "SCA2 short allele BD_CONTROLES - Severity", list(c("CD","NO"), c("CD","CONTROL"), c("NO","CONTROL")))
run_kruskal_analysis(BD_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE1_SCA2", "SCA2 short allele BD_CONTROLES - Pathology type", list(c("TBD-1","Other"), c("TBD-1","CONTROL"), c("Other","CONTROL")))
run_kruskal_analysis(SCH_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE1_SCA2", "SCA2 short allele SCH_CONTROLS - Pathology type", list(c("SCH","Other"), c("Other","CONTROL"), c("SCH","CONTROL")))

# Session info ----
sessionInfo()