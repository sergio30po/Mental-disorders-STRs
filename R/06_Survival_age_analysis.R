# Script name: 06_Survival_age_analysis.R
# ==============================================================================
# Title: Clinical age-at-onset and cross-sectional disease-duration analysis by genotype.

# Author: Sergio Pérez Oliveira

# Description: This script evaluates the association between STR-based genotype
#              classifications (HTT_CODE, ATXN1_CODE, ATXN2_CODE) and clinical
#              variables (age at onset and observed disease duration) in
#              schizophrenia (SCZ) and bipolar disorder (BD) patients.
#              Disease duration is treated strictly as a cross-sectional clinical
#              variable (age at last assessment minus age at onset), not as a
#              time-to-event or progression endpoint.

# Inputs:
#   - Manually selected environment file with custom functions (.R)
#   - Manually selected data file with clinical and genotype variables
#   - Dataframes: DT, BD, SCZ, BD_I, BD_CD, BD_noCD

# Outputs:
#   - Descriptive statistics and exploratory tests for age at onset
#   - Cross-sectional disease-duration comparisons with prespecified FDR control
#   - Machine-readable duration testing plan and master results tables
#   - Supplementary Figure 2: IA vs NORMAL distributions for age at onset and observed disease duration
# ==============================================================================

# Load environment ----
Env_path <- file.choose()
source(Env_path)
rm(Env_path)

fig_dir <- "figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Subsetting by BD-I and CD
BD_I <- subset(BD, BD$PATHOLOGY_TYPE_BINARY == "BD-I")
BD_II <- subset(BD, BD$BD_type == "TBP2")
BD_CD <- subset(BD, BD$CD_BINARY == "CD")
BD_NOCD <- subset(BD, BD$CD_BINARY == "No-CD")

# Subsetting by CD in SCZ
SCZ_P <- subset(SCZ, SCZ$PATHOLOGY_TYPE_BINARY == "SCZ")
SCZ_CD <- subset(SCZ, SCZ$CD_BINARY == "CD")
SCZ_NOCD <- subset(SCZ, SCZ$CD_BINARY == "No-CD")

# ==============================================================================
# 1. AGE-AT-ONSET ANALYSIS
# ==============================================================================

#
# Reviewer-driven revision:
#   Age at onset is retained as a clinical outcome, but all genetic analyses are
#   treated as exploratory and multiple testing is controlled explicitly.
#
# Two complementary exploratory families are evaluated:
#   A) categorical IA vs NORMAL comparisons (HTT, ATXN1, ATXN2; BD and SCZ);
#   B) continuous CAG gene-block models (short, long, quadratic and interaction
#      terms; HTT, ATXN1, ATXN2; BD and SCZ), adjusted for SEX, COFFEE, SMOKER
#      and APOE_E4.
#
# Expanded carriers are excluded gene-by-gene from these inferential analyses
# because they are clinically distinct rare observations and are described
# separately in the manuscript.
#
# No stepwise selection is used in the reviewer analysis. For the continuous CAG
# models, inference is based on the GLOBAL nested-model F-test for the complete
# gene block. Individual polynomial terms are not interpreted unless the global
# block is supported.

revision_dir <- file.path("results", "reviewer_revision")
if (!dir.exists(revision_dir)) {
  dir.create(revision_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# 1A. Categorical IA vs NORMAL age-at-onset comparisons
# ------------------------------------------------------------------------------

run_onset_binary <- function(data,
                             group_col,
                             positive = "IA",
                             reference = "NORMAL",
                             predictor,
                             cohort_scope,
                             family_id = "ONSET_IA_EXPLORATORY",
                             adjustment_method = "BH") {

  needed <- c("ONSET_AGE", group_col)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "Missing column(s) in ", cohort_scope, ": ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  d <- data %>%
    dplyr::transmute(
      onset_age = as.numeric(ONSET_AGE),
      group = as.character(.data[[group_col]])
    ) %>%
    dplyr::filter(
      !is.na(onset_age),
      is.finite(onset_age),
      group %in% c(reference, positive)
    )

  x_ref <- d$onset_age[d$group == reference]
  x_pos <- d$onset_age[d$group == positive]

  n_ref <- length(x_ref)
  n_pos <- length(x_pos)

  if (n_ref < 2L || n_pos < 2L) {
    return(
      tibble::tibble(
        family_id = family_id,
        adjustment_method = adjustment_method,
        cohort_scope = cohort_scope,
        predictor = predictor,
        group_variable = group_col,
        positive = positive,
        reference = reference,
        contrast = paste0(positive, " vs ", reference),
        n_positive = n_pos,
        n_reference = n_ref,
        n_total = n_pos + n_ref,
        mean_positive = if (n_pos > 0) mean(x_pos) else NA_real_,
        sd_positive = if (n_pos > 1) stats::sd(x_pos) else NA_real_,
        median_positive = if (n_pos > 0) stats::median(x_pos) else NA_real_,
        iqr_positive = if (n_pos > 0) stats::IQR(x_pos) else NA_real_,
        mean_reference = if (n_ref > 0) mean(x_ref) else NA_real_,
        sd_reference = if (n_ref > 1) stats::sd(x_ref) else NA_real_,
        median_reference = if (n_ref > 0) stats::median(x_ref) else NA_real_,
        iqr_reference = if (n_ref > 0) stats::IQR(x_ref) else NA_real_,
        median_difference = if (n_pos > 0 && n_ref > 0) {
          stats::median(x_pos) - stats::median(x_ref)
        } else {
          NA_real_
        },
        rank_biserial = NA_real_,
        p_raw = NA_real_,
        estimable = FALSE
      )
    )
  }

  wt <- stats::wilcox.test(
    x = x_pos,
    y = x_ref,
    alternative = "two.sided",
    exact = FALSE,
    correct = FALSE
  )

  pair_diff <- outer(x_pos, x_ref, FUN = "-")
  r_rb <- (
    sum(pair_diff > 0, na.rm = TRUE) -
      sum(pair_diff < 0, na.rm = TRUE)
  ) / (n_pos * n_ref)

  tibble::tibble(
    family_id = family_id,
    adjustment_method = adjustment_method,
    cohort_scope = cohort_scope,
    predictor = predictor,
    group_variable = group_col,
    positive = positive,
    reference = reference,
    contrast = paste0(positive, " vs ", reference),
    n_positive = n_pos,
    n_reference = n_ref,
    n_total = n_pos + n_ref,
    mean_positive = mean(x_pos),
    sd_positive = stats::sd(x_pos),
    median_positive = stats::median(x_pos),
    iqr_positive = stats::IQR(x_pos),
    mean_reference = mean(x_ref),
    sd_reference = stats::sd(x_ref),
    median_reference = stats::median(x_ref),
    iqr_reference = stats::IQR(x_ref),
    median_difference = stats::median(x_pos) - stats::median(x_ref),
    rank_biserial = r_rb,
    p_raw = wt$p.value,
    estimable = TRUE
  )
}

# Derive within-participant ordered alleles consistently with
# 04_CAG_repeat_sizes.R and 05_Regression_models.R.
derive_ordered_onset_alleles <- function(df) {
  for (gene in c("HTT", "ATXN1", "ATXN2")) {
    raw1 <- paste0("ALLELE1_", gene)
    raw2 <- paste0("ALLELE2_", gene)
    short_name <- paste0("SHORT_", gene)
    long_name <- paste0("LONG_", gene)

    if (!all(c(raw1, raw2) %in% names(df))) {
      stop(
        "Missing raw allele column(s) for ", gene, ": ",
        paste(setdiff(c(raw1, raw2), names(df)), collapse = ", "),
        call. = FALSE
      )
    }

    a1 <- suppressWarnings(as.numeric(df[[raw1]]))
    a2 <- suppressWarnings(as.numeric(df[[raw2]]))

    df[[short_name]] <- ifelse(
      is.na(a1) | is.na(a2),
      NA_real_,
      pmin(a1, a2)
    )
    df[[long_name]] <- ifelse(
      is.na(a1) | is.na(a2),
      NA_real_,
      pmax(a1, a2)
    )
  }

  df
}

onset_datasets <- list(
  BD = derive_ordered_onset_alleles(BD),
  SCZ = derive_ordered_onset_alleles(SCZ)
)

onset_genes <- tibble::tribble(
  ~predictor, ~group_col,
  "HTT",      "HTT_CODE",
  "ATXN1",    "ATXN1_CODE",
  "ATXN2",    "ATXN2_CODE"
)

onset_ia_plan <- tidyr::crossing(
  cohort_scope = names(onset_datasets),
  onset_genes
) %>%
  dplyr::mutate(
    positive = "IA",
    reference = "NORMAL",
    family_id = "ONSET_IA_EXPLORATORY",
    adjustment_method = "BH",
    test_id = paste(
      "ONSET", predictor, cohort_scope, positive, "vs", reference,
      sep = "__"
    )
  )

readr::write_csv(
  onset_ia_plan,
  file.path(revision_dir, "06_onset_IA_testing_plan.csv")
)

onset_ia_results_raw <- purrr::pmap_dfr(
  onset_ia_plan,
  function(cohort_scope,
           predictor,
           group_col,
           positive,
           reference,
           family_id,
           adjustment_method,
           test_id) {

    out <- run_onset_binary(
      data = onset_datasets[[cohort_scope]],
      group_col = group_col,
      positive = positive,
      reference = reference,
      predictor = predictor,
      cohort_scope = cohort_scope,
      family_id = family_id,
      adjustment_method = adjustment_method
    )

    dplyr::mutate(out, test_id = test_id, .before = 1)
  }
)

onset_ia_results <- onset_ia_results_raw %>%
  dplyr::mutate(
    p_adj = stats::p.adjust(p_raw, method = "BH"),
    significant_raw = !is.na(p_raw) & p_raw < 0.05,
    significant_adjusted = !is.na(p_adj) & p_adj < 0.05
  ) %>%
  dplyr::arrange(p_adj, p_raw)

readr::write_csv(
  onset_ia_results,
  file.path(revision_dir, "06_onset_IA_master_results.csv")
)

# ------------------------------------------------------------------------------
# 1B. Continuous CAG age-at-onset models: global gene-block inference
# ------------------------------------------------------------------------------

onset_cag_plan <- tibble::tribble(
  ~cohort_scope, ~gene,   ~short_var,     ~long_var,      ~code_var,
  "BD",          "HTT",   "SHORT_HTT",    "LONG_HTT",    "HTT_CODE",
  "BD",          "ATXN1", "SHORT_ATXN1",  "LONG_ATXN1",  "ATXN1_CODE",
  "BD",          "ATXN2", "SHORT_ATXN2",  "LONG_ATXN2",  "ATXN2_CODE",
  "SCZ",         "HTT",   "SHORT_HTT",    "LONG_HTT",    "HTT_CODE",
  "SCZ",         "ATXN1", "SHORT_ATXN1",  "LONG_ATXN1",  "ATXN1_CODE",
  "SCZ",         "ATXN2", "SHORT_ATXN2",  "LONG_ATXN2",  "ATXN2_CODE"
) %>%
  dplyr::mutate(
    family_id = "ONSET_CAG_GLOBAL_EXPLORATORY",
    adjustment_method = "BH",
    test_id = paste("ONSET_CAG_BLOCK", gene, cohort_scope, sep = "__")
  )

readr::write_csv(
  onset_cag_plan,
  file.path(revision_dir, "06_onset_CAG_global_testing_plan.csv")
)

fit_onset_gene_block <- function(data,
                                 cohort_scope,
                                 gene,
                                 short_var,
                                 long_var,
                                 code_var,
                                 test_id) {

  covars <- c("SEX", "COFFEE", "SMOKER", "APOE_E4")
  needed <- c(
    "ONSET_AGE", covars,
    short_var, long_var, code_var
  )

  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "Missing column(s) for ", test_id, ": ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  d <- data %>%
    dplyr::filter(as.character(.data[[code_var]]) != "EXPANDED") %>%
    dplyr::select(dplyr::all_of(needed)) %>%
    dplyr::filter(dplyr::if_all(dplyr::everything(), ~ !is.na(.))) %>%
    dplyr::mutate(
      SEX = droplevels(factor(SEX)),
      COFFEE = droplevels(factor(COFFEE)),
      SMOKER = droplevels(factor(SMOKER)),
      APOE_E4 = droplevels(factor(APOE_E4))
    )

  if (nrow(d) < 20L) {
    stop("Too few complete cases for ", test_id, ".", call. = FALSE)
  }

  # Center alleles before polynomial expansion for numerical stability.
  d$short_c <- as.numeric(d[[short_var]]) - stats::median(as.numeric(d[[short_var]]))
  d$long_c <- as.numeric(d[[long_var]]) - stats::median(as.numeric(d[[long_var]]))

  m_null <- stats::lm(
    ONSET_AGE ~ SEX + COFFEE + SMOKER + APOE_E4,
    data = d
  )

  m_full <- stats::lm(
    ONSET_AGE ~ SEX + COFFEE + SMOKER + APOE_E4 +
      short_c + I(short_c^2) +
      long_c + I(long_c^2) +
      short_c:long_c,
    data = d
  )

  an <- stats::anova(m_null, m_full)

  global <- tibble::tibble(
    test_id = test_id,
    family_id = "ONSET_CAG_GLOBAL_EXPLORATORY",
    adjustment_method = "BH",
    cohort_scope = cohort_scope,
    gene = gene,
    n = stats::nobs(m_full),
    df_added = unname(an$Df[2]),
    f_statistic = unname(an$F[2]),
    p_raw = unname(an$`Pr(>F)`[2]),
    AIC_null = stats::AIC(m_null),
    AIC_full = stats::AIC(m_full),
    delta_AIC_full_minus_null = stats::AIC(m_full) - stats::AIC(m_null),
    adjusted_r2_null = summary(m_null)$adj.r.squared,
    adjusted_r2_full = summary(m_full)$adj.r.squared,
    rank_deficient = m_full$rank < length(stats::coef(m_full))
  )

  coefficients <- broom::tidy(m_full, conf.int = TRUE) %>%
    dplyr::mutate(
      test_id = test_id,
      cohort_scope = cohort_scope,
      gene = gene,
      .before = 1
    )

  list(
    global = global,
    coefficients = coefficients
  )
}

onset_cag_fits <- purrr::pmap(
  onset_cag_plan,
  function(cohort_scope,
           gene,
           short_var,
           long_var,
           code_var,
           family_id,
           adjustment_method,
           test_id) {

    fit_onset_gene_block(
      data = onset_datasets[[cohort_scope]],
      cohort_scope = cohort_scope,
      gene = gene,
      short_var = short_var,
      long_var = long_var,
      code_var = code_var,
      test_id = test_id
    )
  }
)

onset_cag_global <- purrr::map_dfr(onset_cag_fits, "global") %>%
  dplyr::mutate(
    p_adj = stats::p.adjust(p_raw, method = "BH"),
    significant_raw = p_raw < 0.05,
    significant_adjusted = p_adj < 0.05
  ) %>%
  dplyr::arrange(p_adj, p_raw)

onset_cag_coefficients <- purrr::map_dfr(onset_cag_fits, "coefficients") %>%
  dplyr::left_join(
    onset_cag_global %>%
      dplyr::select(test_id, global_block_p = p_raw, global_block_p_adj = p_adj),
    by = "test_id"
  )

readr::write_csv(
  onset_cag_global,
  file.path(revision_dir, "06_onset_CAG_global_results.csv")
)

readr::write_csv(
  onset_cag_coefficients,
  file.path(revision_dir, "06_onset_CAG_coefficients.csv")
)

# ------------------------------------------------------------------------------
# Reviewer-facing outputs
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("AGE AT ONSET: IA TESTS SIGNIFICANT AFTER BH-FDR\n")
cat("============================================================\n")

print(
  onset_ia_results %>%
    dplyr::filter(significant_adjusted) %>%
    dplyr::select(
      test_id, predictor, cohort_scope,
      n_positive, n_reference,
      median_positive, median_reference,
      median_difference, rank_biserial,
      p_raw, p_adj
    ),
  n = Inf,
  width = Inf
)

cat("\n============================================================\n")
cat("AGE AT ONSET: GLOBAL CAG GENE-BLOCK TESTS\n")
cat("============================================================\n")

print(
  onset_cag_global %>%
    dplyr::select(
      test_id, gene, cohort_scope, n,
      f_statistic, p_raw, p_adj,
      delta_AIC_full_minus_null,
      adjusted_r2_null, adjusted_r2_full,
      rank_deficient
    ),
  n = Inf,
  width = Inf
)

cat("\n============================================================\n")
cat("AGE AT ONSET: GLOBAL CAG BLOCKS SIGNIFICANT AFTER BH-FDR\n")
cat("============================================================\n")

print(
  onset_cag_global %>%
    dplyr::filter(significant_adjusted) %>%
    dplyr::select(
      test_id, gene, cohort_scope, n,
      f_statistic, p_raw, p_adj,
      delta_AIC_full_minus_null
    ),
  n = Inf,
  width = Inf
)


# 3. CROSS-SECTIONAL DISEASE-DURATION ANALYSIS ================================
#
# Reviewer-driven revision:
#   DURATION is the observed interval from age at onset to the last clinical
#   assessment. There is no clinical event and no censoring variable.
#
# Therefore:
#   - no survival objects are created;
#   - no Cox proportional-hazards models are fitted;
#   - no log-rank tests are performed;
#   - no Kaplan-Meier curves are generated;
#   - no hazard ratios are reported.
#
# DURATION is retained only as a cross-sectional, exploratory clinical variable.
# All STR duration comparisons below are IA vs NORMAL and exclude EXPANDED
# carriers from inferential testing. The expanded cases remain described
# individually elsewhere in the study.
#
# Multiple testing:
#   STR duration comparisons form one prespecified exploratory family and use
#   Benjamini-Hochberg FDR adjustment across the full family.
#   APOE duration comparisons are kept in a separate exploratory family.

revision_dir <- file.path("results", "reviewer_revision")
if (!dir.exists(revision_dir)) {
  dir.create(revision_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# Generic two-group duration comparison
# ------------------------------------------------------------------------------

run_duration_binary <- function(data,
                                group_col,
                                positive,
                                reference,
                                predictor,
                                cohort_scope,
                                family_id,
                                adjustment_method = "BH") {

  needed <- c("DURATION", group_col)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "Missing column(s) in ", cohort_scope, ": ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  d <- data %>%
    dplyr::transmute(
      duration = as.numeric(DURATION),
      group = as.character(.data[[group_col]])
    ) %>%
    dplyr::filter(
      !is.na(duration),
      is.finite(duration),
      group %in% c(reference, positive)
    )

  x_ref <- d$duration[d$group == reference]
  x_pos <- d$duration[d$group == positive]

  n_ref <- length(x_ref)
  n_pos <- length(x_pos)

  if (n_ref < 2L || n_pos < 2L) {
    return(
      tibble::tibble(
        family_id = family_id,
        adjustment_method = adjustment_method,
        cohort_scope = cohort_scope,
        predictor = predictor,
        group_variable = group_col,
        positive = positive,
        reference = reference,
        contrast = paste0(positive, " vs ", reference),
        n_positive = n_pos,
        n_reference = n_ref,
        n_total = n_pos + n_ref,
        mean_positive = if (n_pos > 0) mean(x_pos) else NA_real_,
        sd_positive = if (n_pos > 1) stats::sd(x_pos) else NA_real_,
        median_positive = if (n_pos > 0) stats::median(x_pos) else NA_real_,
        iqr_positive = if (n_pos > 0) stats::IQR(x_pos) else NA_real_,
        mean_reference = if (n_ref > 0) mean(x_ref) else NA_real_,
        sd_reference = if (n_ref > 1) stats::sd(x_ref) else NA_real_,
        median_reference = if (n_ref > 0) stats::median(x_ref) else NA_real_,
        iqr_reference = if (n_ref > 0) stats::IQR(x_ref) else NA_real_,
        median_difference = if (n_pos > 0 && n_ref > 0) {
          stats::median(x_pos) - stats::median(x_ref)
        } else {
          NA_real_
        },
        rank_biserial = NA_real_,
        p_raw = NA_real_,
        estimable = FALSE
      )
    )
  }

  wt <- stats::wilcox.test(
    x = x_pos,
    y = x_ref,
    alternative = "two.sided",
    exact = FALSE,
    correct = FALSE
  )

  # Rank-biserial effect with explicit direction:
  # positive values = longer observed duration in POSITIVE vs REFERENCE.
  pair_diff <- outer(x_pos, x_ref, FUN = "-")
  r_rb <- (
    sum(pair_diff > 0, na.rm = TRUE) -
      sum(pair_diff < 0, na.rm = TRUE)
  ) / (n_pos * n_ref)

  tibble::tibble(
    family_id = family_id,
    adjustment_method = adjustment_method,
    cohort_scope = cohort_scope,
    predictor = predictor,
    group_variable = group_col,
    positive = positive,
    reference = reference,
    contrast = paste0(positive, " vs ", reference),
    n_positive = n_pos,
    n_reference = n_ref,
    n_total = n_pos + n_ref,
    mean_positive = mean(x_pos),
    sd_positive = stats::sd(x_pos),
    median_positive = stats::median(x_pos),
    iqr_positive = stats::IQR(x_pos),
    mean_reference = mean(x_ref),
    sd_reference = stats::sd(x_ref),
    median_reference = stats::median(x_ref),
    iqr_reference = stats::IQR(x_ref),
    median_difference = stats::median(x_pos) - stats::median(x_ref),
    rank_biserial = r_rb,
    p_raw = wt$p.value,
    estimable = TRUE
  )
}

# ------------------------------------------------------------------------------
# STR duration family: IA vs NORMAL
# ------------------------------------------------------------------------------

duration_datasets <- list(
  BD_ALL = BD,
  BD_I = BD_I,
  BD_CD = BD_CD,
  BD_NOCD = BD_NOCD,
  SCZ_ALL = SCZ,
  SCZ_MAIN_SUBTYPE = SCZ_P,
  SCZ_CD = SCZ_CD,
  SCZ_NOCD = SCZ_NOCD
)

duration_genes <- tibble::tribble(
  ~predictor, ~group_col,
  "HTT",      "HTT_CODE",
  "ATXN1",    "ATXN1_CODE",
  "ATXN2",    "ATXN2_CODE"
)

duration_str_plan <- tidyr::crossing(
  cohort_scope = names(duration_datasets),
  duration_genes
) %>%
  dplyr::mutate(
    positive = "IA",
    reference = "NORMAL",
    family_id = "DURATION_STR_EXPLORATORY",
    adjustment_method = "BH",
    test_id = paste(
      "DURATION", predictor, cohort_scope, positive, "vs", reference,
      sep = "__"
    )
  )

readr::write_csv(
  duration_str_plan,
  file.path(revision_dir, "06_duration_STR_testing_plan.csv")
)

duration_str_results_raw <- purrr::pmap_dfr(
  duration_str_plan,
  function(cohort_scope,
           predictor,
           group_col,
           positive,
           reference,
           family_id,
           adjustment_method,
           test_id) {

    out <- run_duration_binary(
      data = duration_datasets[[cohort_scope]],
      group_col = group_col,
      positive = positive,
      reference = reference,
      predictor = predictor,
      cohort_scope = cohort_scope,
      family_id = family_id,
      adjustment_method = adjustment_method
    )

    dplyr::mutate(out, test_id = test_id, .before = 1)
  }
)

duration_str_results <- duration_str_results_raw %>%
  dplyr::mutate(
    p_adj = stats::p.adjust(p_raw, method = "BH"),
    significant_raw = !is.na(p_raw) & p_raw < 0.05,
    significant_adjusted = !is.na(p_adj) & p_adj < 0.05,
    interpretation = "Cross-sectional observed disease duration; not progression/time-to-event."
  ) %>%
  dplyr::arrange(p_adj, p_raw)

readr::write_csv(
  duration_str_results,
  file.path(revision_dir, "06_duration_STR_master_results.csv")
)

# ------------------------------------------------------------------------------
# APOE duration family: E4+ vs E4-
# ------------------------------------------------------------------------------

duration_apoe_datasets <- c(
  duration_datasets,
  list(BD_II = BD_II)
)

duration_apoe_plan <- tibble::tibble(
  cohort_scope = names(duration_apoe_datasets),
  predictor = "APOE_E4",
  group_col = "APOE_E4",
  positive = "E4+",
  reference = "E4-",
  family_id = "DURATION_APOE_EXPLORATORY",
  adjustment_method = "BH"
) %>%
  dplyr::mutate(
    test_id = paste(
      "DURATION", predictor, cohort_scope, positive, "vs", reference,
      sep = "__"
    )
  )

readr::write_csv(
  duration_apoe_plan,
  file.path(revision_dir, "06_duration_APOE_testing_plan.csv")
)

duration_apoe_results_raw <- purrr::pmap_dfr(
  duration_apoe_plan,
  function(cohort_scope,
           predictor,
           group_col,
           positive,
           reference,
           family_id,
           adjustment_method,
           test_id) {

    out <- run_duration_binary(
      data = duration_apoe_datasets[[cohort_scope]],
      group_col = group_col,
      positive = positive,
      reference = reference,
      predictor = predictor,
      cohort_scope = cohort_scope,
      family_id = family_id,
      adjustment_method = adjustment_method
    )

    dplyr::mutate(out, test_id = test_id, .before = 1)
  }
)

duration_apoe_results <- duration_apoe_results_raw %>%
  dplyr::mutate(
    p_adj = stats::p.adjust(p_raw, method = "BH"),
    significant_raw = !is.na(p_raw) & p_raw < 0.05,
    significant_adjusted = !is.na(p_adj) & p_adj < 0.05,
    interpretation = "Cross-sectional observed disease duration; not progression/time-to-event."
  ) %>%
  dplyr::arrange(p_adj, p_raw)

readr::write_csv(
  duration_apoe_results,
  file.path(revision_dir, "06_duration_APOE_master_results.csv")
)

# ------------------------------------------------------------------------------
# Reviewer-facing outputs
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("CROSS-SECTIONAL DURATION: ATXN2 IA IN BD\n")
cat("No Cox / KM / log-rank interpretation is used.\n")
cat("============================================================\n")

print(
  duration_str_results %>%
    dplyr::filter(
      predictor == "ATXN2",
      cohort_scope == "BD_ALL"
    ),
  width = Inf
)

cat("\n============================================================\n")
cat("STR DURATION TESTS SIGNIFICANT AFTER BH-FDR\n")
cat("============================================================\n")

print(
  duration_str_results %>%
    dplyr::filter(significant_adjusted) %>%
    dplyr::select(
      test_id, predictor, cohort_scope,
      n_positive, n_reference,
      median_positive, median_reference,
      median_difference, rank_biserial,
      p_raw, p_adj
    ),
  n = Inf,
  width = Inf
)

cat("\n============================================================\n")
cat("TEN SMALLEST RAW STR DURATION P-VALUES\n")
cat("============================================================\n")

print(
  duration_str_results %>%
    dplyr::select(
      test_id, predictor, cohort_scope,
      n_positive, n_reference,
      median_positive, median_reference,
      median_difference, rank_biserial,
      p_raw, p_adj
    ) %>%
    dplyr::slice_head(n = 10),
  n = Inf,
  width = Inf
)

# ==============================================================================
# 4. REVISED SUPPLEMENTARY FIGURE 2
# Age at onset and observed disease duration according to IA status
# ==============================================================================
#
# Reviewer-driven figure:
#
#   A) Age at onset in BD
#   B) Age at onset in SCZ
#   C) Observed disease duration in BD
#   D) Observed disease duration in SCZ
#
# For each clinical outcome, HTT, ATXN1 and ATXN2 are shown separately.
# Only NORMAL and IA carriers are included in inferential comparisons.
# EXPANDED alleles are excluded, consistently with the reviewer-facing analyses.
#
# Raw and BH-FDR-adjusted p-values are taken directly from:
#   - onset_ia_results
#   - duration_str_results
#
# No survival, Cox proportional-hazards, Kaplan-Meier, log-rank or progression
# interpretation is used.
# ==============================================================================


# ------------------------------------------------------------------------------
# Figure settings
# ------------------------------------------------------------------------------

clinical_genes <- c(
  "HTT",
  "ATXN1",
  "ATXN2"
)

clinical_code_vars <- c(
  "HTT"   = "HTT_CODE",
  "ATXN1" = "ATXN1_CODE",
  "ATXN2" = "ATXN2_CODE"
)


# ------------------------------------------------------------------------------
# P-value formatting
# ------------------------------------------------------------------------------

format_p_supfig <- function(x) {

  dplyr::case_when(
    is.na(x) ~ "NA",
    x < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", x)
  )
}


# ------------------------------------------------------------------------------
# Convert each clinical outcome to long format
# ------------------------------------------------------------------------------

make_ia_clinical_long <- function(
    data,
    outcome_var
) {

  purrr::imap_dfr(
    clinical_code_vars,
    function(
      code_var,
      gene
    ) {

      tibble::tibble(
        gene = gene,
        group = as.character(
          data[[code_var]]
        ),
        value = suppressWarnings(
          as.numeric(
            data[[outcome_var]]
          )
        )
      )
    }
  ) %>%

    dplyr::filter(
      group %in% c(
        "NORMAL",
        "IA"
      ),
      !is.na(value),
      is.finite(value)
    ) %>%

    dplyr::mutate(
      gene = factor(
        gene,
        levels = clinical_genes
      ),
      group = factor(
        group,
        levels = c(
          "NORMAL",
          "IA"
        ),
        labels = c(
          "Normal",
          "IA"
        )
      )
    )
}


# ------------------------------------------------------------------------------
# Plot datasets
# ------------------------------------------------------------------------------

onset_bd_plot <- make_ia_clinical_long(
  BD,
  "ONSET_AGE"
)

onset_scz_plot <- make_ia_clinical_long(
  SCZ,
  "ONSET_AGE"
)

duration_bd_plot <- make_ia_clinical_long(
  BD,
  "DURATION"
)

duration_scz_plot <- make_ia_clinical_long(
  SCZ,
  "DURATION"
)


# ------------------------------------------------------------------------------
# Statistical annotations
# ------------------------------------------------------------------------------

make_annotation_table <- function(
    plot_data,
    results_data,
    scope
) {
  
  y_range <- range(
    plot_data$value,
    na.rm = TRUE
  )
  
  y_span <- diff(
    y_range
  )
  
  if (
    !is.finite(y_span) ||
    y_span == 0
  ) {
    y_span <- 1
  }
  
  y_by_gene <- plot_data %>%
    dplyr::group_by(
      gene
    ) %>%
    dplyr::summarise(
      y_max = max(
        value,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  results_data %>%
    dplyr::filter(
      cohort_scope == scope,
      predictor %in% clinical_genes
    ) %>%
    dplyr::transmute(
      gene = factor(
        predictor,
        levels = clinical_genes
      ),
      
      n_normal = n_reference,
      n_ia = n_positive,
      
      p_raw = p_raw,
      p_adj = p_adj
    ) %>%
    dplyr::left_join(
      y_by_gene,
      by = "gene"
    ) %>%
    dplyr::mutate(
      y = y_max + 0.07 * y_span,
      
      label = paste0(
        "n = ", n_normal, "/", n_ia,
        "\nraw p = ", format_p_supfig(p_raw),
        "\nBH p = ", format_p_supfig(p_adj)
      )
    )
}


ann_onset_bd <- make_annotation_table(
  onset_bd_plot,
  onset_ia_results,
  "BD"
)

ann_onset_scz <- make_annotation_table(
  onset_scz_plot,
  onset_ia_results,
  "SCZ"
)

ann_duration_bd <- make_annotation_table(
  duration_bd_plot,
  duration_str_results,
  "BD_ALL"
)

ann_duration_scz <- make_annotation_table(
  duration_scz_plot,
  duration_str_results,
  "SCZ_ALL"
)


# ------------------------------------------------------------------------------
# Common plotting function
# ------------------------------------------------------------------------------

make_clinical_ia_panel <- function(
    plot_data,
    annotation_data,
    title,
    y_label
) {

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = gene,
      y = value,
      fill = group
    )
  ) +

    ggplot2::geom_boxplot(
      ggplot2::aes(
        group = interaction(
          gene,
          group
        )
      ),
      position = ggplot2::position_dodge(
        width = 0.72
      ),
      width = 0.58,
      outlier.shape = NA,
      linewidth = 0.45,
      alpha = 0.70
    ) +

    ggplot2::geom_point(
      ggplot2::aes(
        group = group
      ),
      position = ggplot2::position_jitterdodge(
        jitter.width = 0.12,
        dodge.width = 0.72,
        seed = 123
      ),
      shape = 21,
      size = 1.25,
      stroke = 0.20,
      alpha = 0.60
    ) +

    ggplot2::geom_text(
      data = annotation_data,
      ggplot2::aes(
        x = gene,
        y = y,
        label = label
      ),
      inherit.aes = FALSE,
      size = 2.6,
      lineheight = 0.90,
      vjust = 0
    ) +

    ggplot2::scale_fill_manual(
      values = c(
        "Normal" = "grey80",
        "IA" = "#4C78A8"
      ),
      name = "Allele status"
    ) +

    ggplot2::scale_x_discrete(
      labels = c(
        "HTT" = "HTT",
        "ATXN1" = "ATXN1",
        "ATXN2" = "ATXN2"
      )
    ) +

    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(
          0.04,
          0.20
        )
      )
    ) +

    ggplot2::labs(
      title = title,
      x = NULL,
      y = y_label
    ) +

    ggplot2::theme_classic(
      base_size = 10
    ) +

    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5,
        size = 10
      ),
      axis.text.x = ggplot2::element_text(
        face = "italic",
        size = 9
      ),
      axis.title.y = ggplot2::element_text(
        size = 9
      ),
      legend.position = "top",
      legend.title = ggplot2::element_text(
        size = 8
      ),
      legend.text = ggplot2::element_text(
        size = 8
      )
    )
}


# ------------------------------------------------------------------------------
# Build the four final panels
# ------------------------------------------------------------------------------

p_onset_bd <- make_clinical_ia_panel(
  plot_data = onset_bd_plot,
  annotation_data = ann_onset_bd,
  title = "Bipolar disorder",
  y_label = "Age at onset (years)"
)

p_onset_scz <- make_clinical_ia_panel(
  plot_data = onset_scz_plot,
  annotation_data = ann_onset_scz,
  title = "Schizophrenia",
  y_label = "Age at onset (years)"
)

p_duration_bd <- make_clinical_ia_panel(
  plot_data = duration_bd_plot,
  annotation_data = ann_duration_bd,
  title = "Bipolar disorder",
  y_label = "Observed disease duration (years)"
)

p_duration_scz <- make_clinical_ia_panel(
  plot_data = duration_scz_plot,
  annotation_data = ann_duration_scz,
  title = "Schizophrenia",
  y_label = "Observed disease duration (years)"
)


# ------------------------------------------------------------------------------
# Final Supplementary Figure 2
# ------------------------------------------------------------------------------

Supplementary_Figure_2 <-
  (
    p_onset_bd |
      p_onset_scz
  ) /
  (
    p_duration_bd |
      p_duration_scz
  ) +

  patchwork::plot_layout(
    guides = "collect"
  ) +

  patchwork::plot_annotation(
    tag_levels = "A"
  ) &

  ggplot2::theme(
    legend.position = "top",
    plot.tag = ggplot2::element_text(
      face = "bold",
      size = 13
    )
  )


print(
  Supplementary_Figure_2
)


# ------------------------------------------------------------------------------
# Export TIFF
# ------------------------------------------------------------------------------

ggplot2::ggsave(
  filename = file.path(
    fig_dir,
    "Supplementary_Figure_2_age_onset_duration.tiff"
  ),
  plot = Supplementary_Figure_2,
  device = "tiff",
  width = 180,
  height = 170,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)


# ------------------------------------------------------------------------------
# Export PNG
# ------------------------------------------------------------------------------

ggplot2::ggsave(
  filename = file.path(
    fig_dir,
    "Supplementary_Figure_2_age_onset_duration.png"
  ),
  plot = Supplementary_Figure_2,
  device = "png",
  width = 180,
  height = 170,
  units = "mm",
  dpi = 600
)


# Session info ----
sessionInfo()
