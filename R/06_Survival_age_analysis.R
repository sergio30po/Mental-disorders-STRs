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
#   - Supplementary Figure 2 without survival/Cox/Kaplan-Meier panels
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

# DESCRIPTIVE ONSET-AGE / OBSERVED-DURATION DEPENDENCY ----
# This correlation is descriptive only. DURATION is not interpreted as progression.
cor.test(BD$ONSET_AGE, BD$DURATION, method = "spearman", exact = FALSE)
cor.test(SCZ$ONSET_AGE, SCZ$DURATION, method = "spearman", exact = FALSE)

# Sup. Fig. 2A: Descriptive onset-age / observed-duration plots (BD and SCZ) ----
cols_outcome <- c("BD" = "#8CBDE6", "SCZ" = "#F5A04D")

make_cor_panel <- function(df, group = c("BD", "SCZ"),
                           x = "ONSET_AGE", y = "DURATION",
                           title = "") {
  group <- match.arg(group)
  col_use <- cols_outcome[group]
  
  d <- df %>%
    dplyr::select(all_of(c(x, y))) %>%
    dplyr::filter(!is.na(.data[[x]]), !is.na(.data[[y]]))
  
  ggplot(d, aes(x = .data[[x]], y = .data[[y]])) +
    # CI in group color (only ribbon)
    geom_smooth(method = "lm", se = TRUE, aes(fill = col_use),
                color = "black", linewidth = 0.9, alpha = 0.35) +
    # Points: filled by group color + black border
    geom_point(shape = 20, fill = "black", color = "black",
               stroke = 0.35, size = 2, alpha = 0.85) +
    # Correlation text (Spearman)
    stat_cor(method = "spearman", label.x.npc = "middle", label.y.npc = "top") +
    labs(
      x = "Age at onset (years)",
      y = "Disease duration (years)",
      title = title
    ) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_identity()
}

pA_BD  <- make_cor_panel(BD,  group = "BD",  title = "Bipolar disorder")
pA_SCZ <- make_cor_panel(SCZ, group = "SCZ", title = "Schizophrenia")

# Optional: same axes for comparability
#xlim_all <- range(c(BD$ONSET_AGE, SCZ$ONSET_AGE), na.rm = TRUE)
#ylim_all <- range(c(BD$DURATION,  SCZ$DURATION),  na.rm = TRUE)

#pA_BD  <- pA_BD  + coord_cartesian(xlim = xlim_all, ylim = ylim_all)
#pA_SCZ <- pA_SCZ + coord_cartesian(xlim = xlim_all, ylim = ylim_all)

panel_A <- ggarrange(pA_BD, pA_SCZ, ncol = 2, align = "hv")
panel_A
ggsave(
  filename = file.path(fig_dir, "Sup_Fig_2A.tiff"),
  plot = panel_A,
  device = "tiff",
  width = 500, height = 160, units = "mm",
  dpi = 600, compression = "lzw"
)
# 1. AGE-AT-ONSET ANALYSIS ======================================================
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

onset_datasets <- list(
  BD = BD,
  SCZ = SCZ
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
  ~cohort_scope, ~gene,   ~short_var,       ~long_var,        ~code_var,
  "BD",          "HTT",   "ALLELE1_HTT",    "ALLELE2_HTT",    "HTT_CODE",
  "BD",          "ATXN1", "ALLELE1_ATXN1",  "ALLELE2_ATXN1",  "ATXN1_CODE",
  "BD",          "ATXN2", "ALLELE1_ATXN2",  "ALLELE2_ATXN2",  "ATXN2_CODE",
  "SCZ",         "HTT",   "ALLELE1_HTT",    "ALLELE2_HTT",    "HTT_CODE",
  "SCZ",         "ATXN1", "ALLELE1_ATXN1",  "ALLELE2_ATXN1",  "ATXN1_CODE",
  "SCZ",         "ATXN2", "ALLELE1_ATXN2",  "ALLELE2_ATXN2",  "ATXN2_CODE"
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

# ------------------------------------------------------------------------------
# Revised Supplementary Figure 2 -----------------------------------------------
# The former survival panels are removed. Until the reviewer-driven
# age-at-onset global tests are inspected, only the descriptive onset-age /
# observed-duration panel is retained here. No model-selected ATXN1 panel is
# generated automatically.

Sup_Fig_2 <- patchwork::wrap_elements(full = panel_A)

print(Sup_Fig_2)

ggsave(
  filename = file.path(fig_dir, "Supplementary Figure 2.tiff"),
  plot = Sup_Fig_2,
  device = "tiff",
  width = 400, height = 160, units = "mm",
  dpi = 600, compression = "lzw"
)

# Session info ----
sessionInfo()
