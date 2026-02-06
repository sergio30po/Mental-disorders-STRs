# Script name: 04_CAG_repeat_sizes.R
# ==============================================================================
# Title: Comparison of CAG repeat sizes across psychiatric and control groups.

# Author: Sergio Pérez Oliveira

# Description: This script performs Non-parametric statistical analyses to compare
#              CAG repeat lengths of both alleles in HTT, ATXN1 (ATXN1), and ATXN2 (ATXN2)
#              across diagnostic groups (SCZ, BD, Controls), as well as in clinically
#              defined subgroups based on subtype and severity (DCO).
#              The analysis includes:
#                - Descriptive statistics (mean and SD per group)
#                - Kruskal-Wallis test for global group differences
#                - Dunn and Wilcoxon pairwise post-hoc tests if significant
#                - Effect size estimation (Wilcoxon-based) for selected comparisons
#                - CAG repeat sizes and its percentages in the different cohorts.

# Inputs:
#   - Manually selected environment file with custom functions (.R)
#   - Dataframes: DT (full dataset), BP_CONTROLES, SCZ_CONTROLES

# Outputs:
#   - Printed summaries of statistics, test results and effect sizes
# ==============================================================================

# Load environment ----
Env_path <- file.choose()
source(Env_path)
rm(Env_path)

run_kruskal_analysis <- function(df, group_col, value_col,
                                 result_name = NULL,
                                 comparisons_list = NULL,
                                 p_adjust = "holm",
                                 run_pairwise_always = TRUE,
                                 exact_wilcox = FALSE) {
  # Packages needed: FSA, rcompanion, rstatix, dplyr
  # library(FSA); library(rcompanion); library(rstatix); library(dplyr)
  
  if (is.null(result_name)) result_name <- value_col
  
  cat("\n--- Analysis for:", result_name, "---\n")
  
  # Basic checks
  if (!group_col %in% names(df)) stop("group_col not found in df")
  if (!value_col %in% names(df)) stop("value_col not found in df")
  
  # Drop rows with missing group/value to keep n consistent
  dsub <- df[!is.na(df[[group_col]]) & !is.na(df[[value_col]]), , drop = FALSE]
  dsub[[group_col]] <- factor(dsub[[group_col]])
  
  # 1) Descriptives
  cat("Mean and SD by group:\n")
  mean_sd_result <- mean_sd(dsub, group_col, value_col, result_name = value_col)
  print(mean_sd_result)
  
  # 2) Kruskal–Wallis
  cat("\nKruskal–Wallis test:\n")
  kw <- kruskal.test(dsub[[value_col]] ~ dsub[[group_col]])
  print(kw)
  
  # 3) Global effect size: epsilon-squared (ε²)
  cat("\nGlobal effect size (epsilon-squared, ε²):\n")
  eps2 <- rcompanion::epsilonSquared(x = dsub[[value_col]], g = dsub[[group_col]])
  print(eps2)
  
  # Decide whether to run pairwise
  run_pairwise <- run_pairwise_always || (kw$p.value < 0.05)
  
  if (run_pairwise) {
    if (kw$p.value < 0.05) {
      cat("\nPost hoc tests (inferential; KW significant):\n")
    } else {
      cat("\nPost hoc tests (exploratory; KW not significant):\n")
    }
    
    # 4) Dunn post hoc (Holm or chosen adjustment)
    cat("\nDunn's test (p-adjust =", p_adjust, "):\n")
    dunn_result <- FSA::dunnTest(
      dsub[[value_col]] ~ dsub[[group_col]],
      method = p_adjust
    )
    print(dunn_result)
    
    # 5) Pairwise Wilcoxon (same p-adjust)
    cat("\nPairwise Wilcoxon test (p-adjust =", p_adjust, "):\n")
    wilcox_result <- pairwise.wilcox.test(
      x = dsub[[value_col]],
      g = dsub[[group_col]],
      p.adjust.method = p_adjust,
      exact = exact_wilcox
    )
    print(wilcox_result)
    
    # 6) Pairwise effect sizes (rank-biserial) for selected comparisons
    if (!is.null(comparisons_list)) {
      cat("\nPairwise effect sizes (rank-biserial correlation; p-adjust =", p_adjust, "):\n")
      effsize <- dsub %>%
        rstatix::wilcox_effsize(
          as.formula(paste(value_col, "~", group_col)),
          paired = FALSE,
          comparisons = comparisons_list,
          p.adjust.method = p_adjust
        )
      print(effsize)
    } else {
      cat("\nNo comparisons_list provided for effect sizes.\n")
    }
    
  } else {
    cat("\nPairwise tests were not run.\n")
  }
  
  # Return objects invisibly (useful for saving)
  invisible(list(
    descriptives = mean_sd_result,
    kw = kw,
    eps2 = eps2,
    dunn = if (run_pairwise) dunn_result else NULL,
    pairwise_wilcox = if (run_pairwise) wilcox_result else NULL,
    pairwise_effsize = if (run_pairwise && !is.null(comparisons_list)) effsize else NULL
  ))
}



# HTT ANALYSIS ----

### Long allele
run_kruskal_analysis(DT, "PATHOLOGY", "ALLELE2_HTT", "HTT long allele", list(c("BD","SCZ"), c("BD","CONTROL"), c("SCZ","CONTROL")))

run_kruskal_analysis(BD_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE2_HTT", "HTT long allele BD_CONTROLS - Pathology type", list(c("BD-I","Other"), c("BD-I","CONTROL"), c("Other","CONTROL")))
run_kruskal_analysis(BD_CONTROLS, "CD_BINARY", "ALLELE2_HTT", "HTT long allele BD_CONTROLS - Severity", list(c("CD","No-CD"), c("CD","CONTROL"), c("No-CD","CONTROL")))

run_kruskal_analysis(SCZ_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE2_HTT", "HTT long allele SCZ_CONTROLS - Pathology type", list(c("SCZ","Other"), c("Other","CONTROL"), c("SCZ","CONTROL")))
run_kruskal_analysis(SCZ_CONTROLS, "CD_BINARY", "ALLELE2_HTT", "HTT long allele SCZ_CONTROLS - Severity", list(c("CD","No-CD"), c("CD","CONTROL"), c("No-CD","CONTROL")))

### Short allele
run_kruskal_analysis(DT, "PATHOLOGY", "ALLELE1_HTT", "HTT short allele", list(c("BD","SCZ"), c("BD","CONTROL"), c("SCZ","CONTROL")))

run_kruskal_analysis(BD_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE1_HTT", "HTT short allele BD_CONTROLES - Pathology type", list(c("BD-I","Other"), c("BD-I","CONTROL"), c("Other","CONTROL")))
run_kruskal_analysis(BD_CONTROLS, "CD_BINARY", "ALLELE1_HTT", "HTT short allele BD_CONTROLES - Severity", list(c("CD","No-CD"), c("CD","CONTROL"), c("No-CD","CONTROL")))

run_kruskal_analysis(SCZ_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE1_HTT", "HTT short allele SCZ_CONTROLS - Pathology type", list(c("SCZ","Other"), c("Other","CONTROL"), c("SCZ","CONTROL")))
run_kruskal_analysis(SCZ_CONTROLS, "CD_BINARY", "ALLELE1_HTT", "HTT short allele SCZ_CONTROLS - Severity", list(c("CD","No-CD"), c("CD","CONTROL"), c("No-CD","CONTROL")))

# ATXN1 ANALYSIS ----

### Long allele
run_kruskal_analysis(DT, "PATHOLOGY", "ALLELE2_ATXN1", "ATXN1 long allele", list(c("BD","SCZ"), c("BD","CONTROL"), c("SCZ","CONTROL")))

run_kruskal_analysis(BD_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE2_ATXN1", "ATXN1 long allele BD_CONTROLES - Pathology type", list(c("BD-I","Other"), c("BD-I","CONTROL"), c("Other","CONTROL")))
run_kruskal_analysis(BD_CONTROLS, "CD_BINARY", "ALLELE2_ATXN1", "ATXN1 long allele BD_CONTROLES - Severity", list(c("CD","No-CD"), c("CD","CONTROL"), c("No-CD","CONTROL")))

run_kruskal_analysis(SCZ_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE2_ATXN1", "ATXN1 long allele SCZ_CONTROLS - Pathology type", list(c("SCZ","Other"), c("Other","CONTROL"), c("SCZ","CONTROL")))
run_kruskal_analysis(SCZ_CONTROLS, "CD_BINARY", "ALLELE2_ATXN1", "ATXN1 long allele SCZ_CONTROLS - Severity", list(c("CD","No-CD"), c("CD","CONTROL"), c("No-CD","CONTROL")))

### Short allele
run_kruskal_analysis(DT, "PATHOLOGY", "ALLELE1_ATXN1", "ATXN1 short allele", list(c("BD","SCZ"), c("BD","CONTROL"), c("SCZ","CONTROL")))

run_kruskal_analysis(BD_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE1_ATXN1", "ATXN1 short allele BD_CONTROLES - Pathology type", list(c("BD-I","Other"), c("BD-I","CONTROL"), c("Other","CONTROL")))
run_kruskal_analysis(BD_CONTROLS, "CD_BINARY", "ALLELE1_ATXN1", "ATXN1 short allele BD_CONTROLES - Severity", list(c("CD","No-CD"), c("CD","CONTROL"), c("No-CD","CONTROL")))

run_kruskal_analysis(SCZ_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE1_ATXN1", "ATXN1 short allele SCZ_CONTROLS - Pathology type", list(c("SCZ","Other"), c("Other","CONTROL"), c("SCZ","CONTROL")))
run_kruskal_analysis(SCZ_CONTROLS, "CD_BINARY", "ALLELE1_ATXN1", "ATXN1 short allele SCZ_CONTROLS - Severity", list(c("CD","No-CD"), c("CD","CONTROL"), c("No-CD","CONTROL")))

# ATXN2 ANALYSIS ----

### Long allele
run_kruskal_analysis(DT, "PATHOLOGY", "ALLELE2_ATXN2", "ATXN2 long allele", list(c("BD","SCZ"), c("BD","CONTROL"), c("SCZ","CONTROL")))

run_kruskal_analysis(BD_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE2_ATXN2", "ATXN2 long allele BD_CONTROLES - Pathology type", list(c("BD-I","Other"), c("BD-I","CONTROL"), c("Other","CONTROL")))
run_kruskal_analysis(BD_CONTROLS, "CD_BINARY", "ALLELE2_ATXN2", "ATXN2 long allele BD_CONTROLES - Severity", list(c("CD","No-CD"), c("CD","CONTROL"), c("No-CD","CONTROL")))

run_kruskal_analysis(SCZ_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE2_ATXN2", "ATXN2 long allele SCZ_CONTROLS - Pathology type", list(c("SCZ","Other"), c("Other","CONTROL"), c("SCZ","CONTROL")))
run_kruskal_analysis(SCZ_CONTROLS, "CD_BINARY", "ALLELE2_ATXN2", "ATXN2 long allele SCZ_CONTROLS - Severity", list(c("CD","No-CD"), c("CD","CONTROL"), c("No-CD","CONTROL")))

### Short allele
run_kruskal_analysis(DT, "PATHOLOGY", "ALLELE1_ATXN2", "ATXN2 short allele", list(c("BD","SCZ"), c("BD","CONTROL"), c("SCZ","CONTROL")))

run_kruskal_analysis(BD_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE1_ATXN2", "ATXN2 short allele BD_CONTROLES - Pathology type", list(c("BD-I","Other"), c("BD-I","CONTROL"), c("Other","CONTROL")))
run_kruskal_analysis(BD_CONTROLS, "CD_BINARY", "ALLELE1_ATXN2", "ATXN2 short allele BD_CONTROLES - Severity", list(c("CD","No-CD"), c("CD","CONTROL"), c("No-CD","CONTROL")))

run_kruskal_analysis(SCZ_CONTROLS, "PATHOLOGY_TYPE_BINARY", "ALLELE1_ATXN2", "ATXN2 short allele SCZ_CONTROLS - Pathology type", list(c("SCZ","Other"), c("Other","CONTROL"), c("SCZ","CONTROL")))
run_kruskal_analysis(SCZ_CONTROLS, "CD_BINARY", "ALLELE1_ATXN2", "ATXN2 short allele SCZ_CONTROLS - Severity", list(c("CD","No-CD"), c("CD","CONTROL"), c("No-CD","CONTROL")))

#CAG repeats frequencies among cohorts----

allele_percentages_by_group <- function(
    data,
    group_var      = "PATHOLOGY",
    allele_cols    = c("ALLELE1_HTT", "ALLELE2_HTT"),  # <- pass ATXN1/ATXN2 here when needed
    group_to_show  = NULL,                              # e.g. "BD", "SCZ", "CONTROL"
    output         = c("long", "wide"),                 # output format
    digits         = 2,
    drop_na        = TRUE
) {
  output <- match.arg(output)
  
  # Check that all requested columns exist in the dataset
  missing_cols <- setdiff(c(group_var, allele_cols), names(data))
  if (length(missing_cols) > 0) {
    stop("These columns were No-CDt found in 'data': ", paste(missing_cols, collapse = ", "))
  }
  
  # Convert selected allele columns to long format
  allele_long <- data %>%
    select(all_of(group_var), all_of(allele_cols)) %>%
    pivot_longer(
      cols      = all_of(allele_cols),
      names_to  = "ALLELE_COL",
      values_to = "CAG_size"
    )
  
  # Optionally remove rows with missing allele values
  if (drop_na) {
    allele_long <- allele_long %>% filter(!is.na(CAG_size))
  }
  
  # Compute frequency and percentage by group and allele size
  allele_freq <- allele_long %>%
    group_by(across(all_of(group_var)), CAG_size) %>%
    summarise(count = n(), .groups = "drop_last") %>%
    mutate(
      percentage = round(100 * count / sum(count), digits)  # % within each group
    ) %>%
    arrange(across(all_of(group_var)), CAG_size)
  
  # Optionally filter to show only one specific group
  if (!is.null(group_to_show)) {
    allele_freq <- allele_freq %>%
      filter(!!sym(group_var) == group_to_show)
  }
  
  # Return results in wide format if requested
  if (output == "wide") {
    res <- allele_freq %>%
      select(all_of(group_var), CAG_size, percentage) %>%
      pivot_wider(
        names_from  = all_of(group_var),
        values_from = percentage,
        values_fill = 0
      ) %>%
      arrange(CAG_size)
  } else {
    res <- allele_freq
  }
  
  # Print and return the result invisibly
  print(res)
  invisible(res)
}

allele_percentages_by_group(DT,  allele_cols = c("ALLELE1_HTT", "ALLELE2_HTT"),  group_to_show = "CONTROL")
allele_percentages_by_group(DT,  allele_cols = c("ALLELE1_HTT", "ALLELE2_HTT"),  group_to_show = "BD")
allele_percentages_by_group(DT,  allele_cols = c("ALLELE1_HTT", "ALLELE2_HTT"),  group_to_show = "SCZ")

allele_percentages_by_group(DT,  allele_cols = c("ALLELE1_ATXN1", "ALLELE2_ATXN1"),  group_to_show = "CONTROL")
allele_percentages_by_group(DT,  allele_cols = c("ALLELE1_ATXN1", "ALLELE2_ATXN1"),  group_to_show = "BD")
allele_percentages_by_group(DT,  allele_cols = c("ALLELE1_ATXN1", "ALLELE2_ATXN1"),  group_to_show = "SCZ")

allele_percentages_by_group(DT,  allele_cols = c("ALLELE1_ATXN2", "ALLELE2_ATXN2"),  group_to_show = "CONTROL")
allele_percentages_by_group(DT,  allele_cols = c("ALLELE1_ATXN2", "ALLELE2_ATXN2"),  group_to_show = "BD")
allele_percentages_by_group(DT,  allele_cols = c("ALLELE1_ATXN2", "ALLELE2_ATXN2"),  group_to_show = "SCZ")


# Session info ----
sessionInfo()