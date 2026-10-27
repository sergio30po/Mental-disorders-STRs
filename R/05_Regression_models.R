#Script name: 05_Regression_models.R
# ==============================================================================
# Title: Psychiatric Disorder Risk and Age-of-Onset Regression Modeling.

# Author: Sergio Pérez Oliveira

# Description: Reviewer-driven regression analyses for:
#              1) the direct HTT intermediate-allele association in BD-I vs controls,
#                 with explicit case/reference coding and covariate adjustment;
#              2) prespecified continuous CAG gene-block models for BD/SCZ vs controls.
#
#              No stepwise, greedy, or AIC-based variable/gene selection is used for
#              reviewer-facing inference. Individual polynomial coefficients are used
#              only to describe model shape; the primary continuous-model inference is
#              the global likelihood-ratio test for the complete gene block.

# Inputs:
#   - Environment file with custom functions (manually selected)
#   - Data file with preprocessed clinical/genetic data
#   - Dataframes: DT, MENTAL, BD, SCZ

# Outputs:
#   - Direct crude/adjusted HTT-IA BD-I logistic models
#   - Prespecified multinomial continuous-CAG models
#   - Global gene-block likelihood-ratio tests with Holm adjustment
#   - Hierarchical ATXN2 component likelihood-ratio follow-up tests
#   - Reviewer-revision CSV outputs
#   - Figure 2

# ==============================================================================

# Load environment  ----
Env_path <- file.choose()
source(Env_path)
rm(Env_path)

# 0. REVIEWER ANALYSIS: HTT IA IN BD-I VS CONTROLS ============================
# Reviewer 1 requested direct clarification of whether the main HTT intermediate-
# allele association remains after APOE adjustment and explicit case/reference
# coding.
#
# Primary orientation:
#   outcome:   BD-I = 1, CONTROL = 0
#   predictor: HTT IA = 1, HTT NORMAL = 0
#
# The crude and adjusted models are fitted on the SAME complete-case dataset so
# that changes in the HTT estimate reflect covariate adjustment rather than
# different missing-data patterns.
#
# The adjusted model is:
#   BD-I ~ HTT_IA + SEX + AGE + APOE_E4
#
# This is a prespecified covariate-adjusted sensitivity analysis of the same HTT
# BD-I association tested in 03_Genotype_stats.R. It is not treated as a new
# discovery family for multiplicity purposes.

revision_dir <- file.path("results", "reviewer_revision")
if (!dir.exists(revision_dir)) {
  dir.create(revision_dir, recursive = TRUE)
}

required_htt_cols <- c(
  "PATHOLOGY_TYPE_BINARY", "HTT_CODE", "SEX", "AGE", "APOE_E4"
)

missing_htt_cols <- setdiff(required_htt_cols, names(BD_CONTROLS))
if (length(missing_htt_cols) > 0) {
  stop(
    "Missing required column(s) for reviewer HTT model: ",
    paste(missing_htt_cols, collapse = ", "),
    call. = FALSE
  )
}

htt_bdi <- BD_CONTROLS %>%
  dplyr::filter(
    as.character(PATHOLOGY_TYPE_BINARY) %in% c("BD-I", "CONTROL"),
    as.character(HTT_CODE) %in% c("NORMAL", "IA")
  ) %>%
  dplyr::transmute(
    outcome = factor(
      as.character(PATHOLOGY_TYPE_BINARY),
      levels = c("CONTROL", "BD-I")
    ),
    case = as.integer(as.character(PATHOLOGY_TYPE_BINARY) == "BD-I"),
    HTT_IA = as.integer(as.character(HTT_CODE) == "IA"),
    SEX = factor(SEX),
    AGE = as.numeric(AGE),
    APOE_E4 = factor(APOE_E4)
  )

# Explicit factor references
if (!"Male" %in% levels(htt_bdi$SEX)) {
  stop("Expected SEX reference level 'Male' was not found.", call. = FALSE)
}
if (!"E4-" %in% levels(htt_bdi$APOE_E4)) {
  stop("Expected APOE_E4 reference level 'E4-' was not found.", call. = FALSE)
}

htt_bdi <- htt_bdi %>%
  dplyr::mutate(
    SEX = stats::relevel(SEX, ref = "Male"),
    APOE_E4 = stats::relevel(APOE_E4, ref = "E4-")
  )

# Common complete-case dataset for crude and adjusted models
htt_bdi_cc <- htt_bdi %>%
  dplyr::filter(
    !is.na(case),
    !is.na(HTT_IA),
    !is.na(SEX),
    !is.na(AGE),
    !is.na(APOE_E4)
  ) %>%
  droplevels()

if (length(unique(htt_bdi_cc$case)) != 2L) {
  stop("Both CONTROL and BD-I must be present in the complete-case dataset.", call. = FALSE)
}
if (length(unique(htt_bdi_cc$HTT_IA)) != 2L) {
  stop("Both NORMAL and IA HTT categories must be present.", call. = FALSE)
}

# Descriptive counts used by the models
htt_model_counts <- htt_bdi_cc %>%
  dplyr::mutate(
    HTT_status = factor(
      ifelse(HTT_IA == 1L, "IA", "NORMAL"),
      levels = c("NORMAL", "IA")
    )
  ) %>%
  dplyr::count(outcome, HTT_status, .drop = FALSE) %>%
  dplyr::group_by(outcome) %>%
  dplyr::mutate(
    group_n = sum(n),
    percentage = 100 * n / group_n
  ) %>%
  dplyr::ungroup()

readr::write_csv(
  htt_model_counts,
  file.path(revision_dir, "05_HTT_BDI_model_counts.csv")
)

# Models
model_HTT_BDI_crude <- stats::glm(
  case ~ HTT_IA,
  data = htt_bdi_cc,
  family = stats::binomial()
)

model_HTT_BDI_adjusted_null <- stats::glm(
  case ~ SEX + AGE + APOE_E4,
  data = htt_bdi_cc,
  family = stats::binomial()
)

model_HTT_BDI_adjusted <- stats::glm(
  case ~ HTT_IA + SEX + AGE + APOE_E4,
  data = htt_bdi_cc,
  family = stats::binomial()
)

# Helper: Wald OR/CI for one coefficient
# NOTE: use "fit" rather than "model" as the function argument. Inside
# tibble::tibble(), columns are evaluated sequentially; creating a column named
# "model" can otherwise mask an object/argument also named "model".
extract_glm_term_wald <- function(fit, term, model_label) {
  co <- summary(fit)$coefficients

  if (!term %in% rownames(co)) {
    stop("Term '", term, "' not found in model ", model_label, ".", call. = FALSE)
  }

  beta <- unname(co[term, "Estimate"])
  se <- unname(co[term, "Std. Error"])
  p <- unname(co[term, "Pr(>|z|)"])

  model_frame <- fit$model
  n_total_fit <- nrow(model_frame)
  n_case_fit <- sum(model_frame$case == 1L)
  n_reference_fit <- sum(model_frame$case == 0L)

  tibble::tibble(
    model = model_label,
    term = term,
    effect_definition = "OR = odds(BD-I vs CONTROL) for HTT IA vs NORMAL",
    n_total = n_total_fit,
    n_case = n_case_fit,
    n_reference = n_reference_fit,
    beta_log_odds = beta,
    std_error = se,
    odds_ratio = exp(beta),
    conf_low = exp(beta - 1.96 * se),
    conf_high = exp(beta + 1.96 * se),
    p_wald = p,
    AIC = stats::AIC(fit)
  )
}

htt_bdi_model_results <- dplyr::bind_rows(
  extract_glm_term_wald(
    model_HTT_BDI_crude,
    "HTT_IA",
    "Crude: HTT IA"
  ),
  extract_glm_term_wald(
    model_HTT_BDI_adjusted,
    "HTT_IA",
    "Adjusted: HTT IA + SEX + AGE + APOE_E4"
  )
)

# Likelihood-ratio test for the adjusted HTT contribution
lrt_adjusted <- stats::anova(
  model_HTT_BDI_adjusted_null,
  model_HTT_BDI_adjusted,
  test = "LRT"
)

p_lrt_adjusted <- lrt_adjusted$`Pr(>Chi)`[2]

htt_bdi_model_results <- htt_bdi_model_results %>%
  dplyr::mutate(
    case = "BD-I",
    reference = "CONTROL",
    predictor_case = "HTT IA",
    predictor_reference = "HTT NORMAL",
    family_id = "A_IA_DIAGNOSIS_SENSITIVITY",
    multiplicity_note = paste(
      "Prespecified covariate-adjusted sensitivity of the HTT BD-I association;",
      "interpret alongside Family A multiplicity results from 03_Genotype_stats.R."
    ),
    p_lrt_HTT_adjusted = dplyr::if_else(
      grepl("^Adjusted", model),
      p_lrt_adjusted,
      NA_real_
    )
  ) %>%
  dplyr::select(
    family_id, model, case, reference,
    predictor_case, predictor_reference,
    n_total, n_case, n_reference,
    odds_ratio, conf_low, conf_high,
    p_wald, p_lrt_HTT_adjusted,
    beta_log_odds, std_error, AIC,
    effect_definition, multiplicity_note
  )

readr::write_csv(
  htt_bdi_model_results,
  file.path(revision_dir, "05_HTT_BDI_logistic_results.csv")
)

cat("\n============================================================\n")
cat("REVIEWER MODEL: HTT IA, BD-I vs CONTROL\n")
cat("Outcome orientation: BD-I = 1; CONTROL = 0\n")
cat("Predictor orientation: HTT IA = 1; NORMAL = 0\n")
cat("Crude and adjusted models use the same complete-case sample.\n")
cat("============================================================\n")
print(htt_model_counts, n = Inf)
print(htt_bdi_model_results, width = Inf)

# ------------------------------------------------------------------------------
# APOE x HTT interaction feasibility check
# ------------------------------------------------------------------------------

htt_apoe_cells <- htt_bdi_cc %>%
  dplyr::mutate(
    HTT_status = factor(
      ifelse(HTT_IA == 1L, "IA", "NORMAL"),
      levels = c("NORMAL", "IA")
    )
  ) %>%
  dplyr::count(outcome, HTT_status, APOE_E4, .drop = FALSE)

readr::write_csv(
  htt_apoe_cells,
  file.path(revision_dir, "05_HTT_BDI_APOE_cells.csv")
)

# Prespecified stability screen:
# do not fit/interpret an APOE interaction if any diagnosis x HTT x APOE cell
# contains fewer than 5 individuals.
min_apoe_cell <- min(htt_apoe_cells$n)

interaction_decision <- tibble::tibble(
  minimum_cell_count = min_apoe_cell,
  threshold = 5L,
  interaction_fitted = min_apoe_cell >= 5L,
  rule = paste(
    "HTT_IA x APOE_E4 interaction is fitted only if every",
    "diagnosis x HTT x APOE cell has n >= 5."
  )
)

readr::write_csv(
  interaction_decision,
  file.path(revision_dir, "05_HTT_BDI_APOE_interaction_decision.csv")
)

if (min_apoe_cell >= 5L) {

  model_HTT_BDI_APOE_interaction <- stats::glm(
    case ~ HTT_IA * APOE_E4 + SEX + AGE,
    data = htt_bdi_cc,
    family = stats::binomial()
  )

  co <- summary(model_HTT_BDI_APOE_interaction)$coefficients

  htt_apoe_interaction_results <- tibble::tibble(
    term = rownames(co),
    beta_log_odds = unname(co[, "Estimate"]),
    std_error = unname(co[, "Std. Error"]),
    odds_ratio = exp(beta_log_odds),
    conf_low = exp(beta_log_odds - 1.96 * std_error),
    conf_high = exp(beta_log_odds + 1.96 * std_error),
    p_wald = unname(co[, "Pr(>|z|)"])
  )

  readr::write_csv(
    htt_apoe_interaction_results,
    file.path(revision_dir, "05_HTT_BDI_APOE_interaction_results.csv")
  )

  cat("\nAPOE x HTT interaction model fitted (all cells n >= 5).\n")
  print(htt_apoe_interaction_results, width = Inf)

} else {

  cat(
    "\nAPOE x HTT interaction NOT fitted: minimum diagnosis x HTT x APOE cell n = ",
    min_apoe_cell,
    " (<5). This sensitivity analysis would be too sparse for stable interpretation.\n",
    sep = ""
  )
}

# ------------------------------------------------------------------------------
# Age-overlap diagnostic
# ------------------------------------------------------------------------------

htt_age_by_group <- htt_bdi_cc %>%
  dplyr::group_by(outcome) %>%
  dplyr::summarise(
    n = dplyr::n(),
    min_age = min(AGE),
    p25_age = stats::quantile(AGE, 0.25),
    median_age = stats::median(AGE),
    p75_age = stats::quantile(AGE, 0.75),
    max_age = max(AGE),
    .groups = "drop"
  )

common_age_min <- max(htt_age_by_group$min_age)
common_age_max <- min(htt_age_by_group$max_age)

htt_age_overlap_counts <- htt_bdi_cc %>%
  dplyr::mutate(
    in_common_age_range = AGE >= common_age_min & AGE <= common_age_max
  ) %>%
  dplyr::count(outcome, in_common_age_range)

htt_age_overlap_summary <- tibble::tibble(
  common_age_min = common_age_min,
  common_age_max = common_age_max,
  common_age_range_exists = common_age_min <= common_age_max
)

readr::write_csv(
  htt_age_by_group,
  file.path(revision_dir, "05_HTT_BDI_age_by_group.csv")
)
readr::write_csv(
  htt_age_overlap_counts,
  file.path(revision_dir, "05_HTT_BDI_age_overlap_counts.csv")
)
readr::write_csv(
  htt_age_overlap_summary,
  file.path(revision_dir, "05_HTT_BDI_age_overlap_summary.csv")
)

# ------------------------------------------------------------------------------
# APOE-specific sensitivity: does adding APOE materially change the HTT estimate?
# ------------------------------------------------------------------------------
# Same complete-case sample, same SEX + AGE adjustment. The only difference
# between these two models is the inclusion of APOE_E4.

model_HTT_BDI_no_APOE <- stats::glm(
  case ~ HTT_IA + SEX + AGE,
  data = htt_bdi_cc,
  family = stats::binomial()
)

model_HTT_BDI_with_APOE <- stats::glm(
  case ~ HTT_IA + SEX + AGE + APOE_E4,
  data = htt_bdi_cc,
  family = stats::binomial()
)

htt_apoe_specific <- dplyr::bind_rows(
  extract_glm_term_wald(
    model_HTT_BDI_no_APOE,
    "HTT_IA",
    "Adjusted: HTT IA + SEX + AGE"
  ),
  extract_glm_term_wald(
    model_HTT_BDI_with_APOE,
    "HTT_IA",
    "Adjusted: HTT IA + SEX + AGE + APOE_E4"
  )
)

htt_apoe_lrt <- stats::anova(
  model_HTT_BDI_no_APOE,
  model_HTT_BDI_with_APOE,
  test = "LRT"
)

apoe_added_p <- htt_apoe_lrt$`Pr(>Chi)`[2]

htt_apoe_specific <- htt_apoe_specific %>%
  dplyr::mutate(
    apoe_added_lrt_p = apoe_added_p,
    interpretation = paste(
      "Comparison isolates the contribution of adding APOE_E4 to the",
      "SEX + AGE adjusted HTT model."
    )
  )

readr::write_csv(
  htt_apoe_specific,
  file.path(revision_dir, "05_HTT_BDI_APOE_specific_sensitivity.csv")
)

readr::write_csv(
  tibble::as_tibble(htt_apoe_lrt, rownames = "model_step"),
  file.path(revision_dir, "05_HTT_BDI_APOE_added_LRT.csv")
)

cat("\n============================================================\n")
cat("APOE-SPECIFIC SENSITIVITY: HTT IA, BD-I vs CONTROL\n")
cat("============================================================\n")
print(htt_apoe_specific, width = Inf)
print(htt_apoe_lrt)


cat("\n============================================================\n")
cat("AGE OVERLAP DIAGNOSTIC: BD-I vs CONTROL\n")
cat("============================================================\n")
print(htt_age_by_group, width = Inf)
print(htt_age_overlap_summary, width = Inf)
print(htt_age_overlap_counts, n = Inf, width = Inf)

# ------------------------------------------------------------------------------
# Sensitivity: flexible age adjustment
# ------------------------------------------------------------------------------
# Controls and BD-I cases have markedly different age distributions. To reduce
# dependence on a strictly linear AGE effect, refit the adjusted model using a
# natural spline with 3 degrees of freedom. This is a sensitivity analysis of the
# HTT coefficient, not a separate discovery test.

model_HTT_BDI_age_spline <- stats::glm(
  case ~ HTT_IA + SEX + splines::ns(AGE, df = 3) + APOE_E4,
  data = htt_bdi_cc,
  family = stats::binomial()
)

htt_age_spline_result <- extract_glm_term_wald(
  model_HTT_BDI_age_spline,
  "HTT_IA",
  "Sensitivity: HTT IA + SEX + ns(AGE, df=3) + APOE_E4"
) %>%
  dplyr::mutate(
    case = "BD-I",
    reference = "CONTROL",
    predictor_case = "HTT IA",
    predictor_reference = "HTT NORMAL",
    family_id = "A_IA_DIAGNOSIS_SENSITIVITY",
    multiplicity_note = paste(
      "Flexible-age sensitivity of the prespecified HTT BD-I association;",
      "interpret alongside Family A multiplicity results from 03_Genotype_stats.R."
    )
  ) %>%
  dplyr::select(
    family_id, model, case, reference,
    predictor_case, predictor_reference,
    n_total, n_case, n_reference,
    odds_ratio, conf_low, conf_high,
    p_wald, beta_log_odds, std_error, AIC,
    effect_definition, multiplicity_note
  )

readr::write_csv(
  htt_age_spline_result,
  file.path(revision_dir, "05_HTT_BDI_age_spline_sensitivity.csv")
)

cat("\n============================================================\n")
cat("AGE-SPLINE SENSITIVITY: HTT IA, BD-I vs CONTROL\n")
cat("============================================================\n")
print(htt_age_spline_result, width = Inf)


# 1. CONTINUOUS CAG DISEASE-RISK MODELS (PRESPECIFIED MULTINOMIAL) ==========
# Reviewer-facing continuous sensitivity analysis.
#
# Outcome:
#   PATHOLOGY with CONTROL as the reference category.
#
# Covariates:
#   SEX + AGE + APOE_E4
#
# Gene block:
#   short + short^2 + long + long^2 + short:long
#
# Primary inference:
#   For each gene, compare the covariate-only multinomial model with the model
#   containing the complete five-term CAG block using a likelihood-ratio test.
#   The three global gene-block tests form family D_CAG_MULTINOMIAL_GLOBAL and
#   are adjusted using Holm.
#
# Important:
#   No AIC/stepwise/greedy model selection is used.
#   Coefficient-level p-values are retained for model-shape inspection only.

# Data prep ---------------------------------------------------------------------
# Work on a local modeling copy so downstream scripts receive the original DT.
# In particular, filtering expanded/missing genotype codes here must not mutate the
# shared cohort object later used by the HTT meta-analysis.
DT_model <- DT %>%
  dplyr::mutate(
    PATHOLOGY = stats::relevel(factor(PATHOLOGY), ref = "CONTROL"),
    SEX = factor(SEX)
  ) %>%
  dplyr::filter(
    ATXN1_CODE != "EXPANDED",
    ATXN2_CODE != "EXPANDED",
    HTT_CODE != "EXPANDED"
  ) %>%
  dplyr::mutate(
    ATXN1_CODE = droplevels(factor(ATXN1_CODE)),
    ATXN2_CODE = droplevels(factor(ATXN2_CODE)),
    HTT_CODE = droplevels(factor(HTT_CODE))
  )

# ALLELE1/ALLELE2 are raw diploid fields. Derive within-person short/long
# consistently with 04_CAG_repeat_sizes.R.
derive_short_long <- function(df, gene) {
  a1_name <- paste0("ALLELE1_", gene)
  a2_name <- paste0("ALLELE2_", gene)
  short_name <- paste0("SHORT_", gene)
  long_name <- paste0("LONG_", gene)

  a1 <- suppressWarnings(as.numeric(df[[a1_name]]))
  a2 <- suppressWarnings(as.numeric(df[[a2_name]]))

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
  df
}

for (g in c("HTT", "ATXN1", "ATXN2")) {
  DT_model <- derive_short_long(DT_model, g)
}

drop_na_vars <- function(df, vars) {
  dplyr::filter(df, dplyr::if_all(dplyr::all_of(vars), ~ !is.na(.)))
}

# Use one common complete-case dataset for all three gene-block comparisons.
# This avoids changes in N when comparing the global evidence across genes.
vars_allele <- c(
  "PATHOLOGY", "SEX", "AGE", "APOE_E4",
  "SHORT_HTT", "LONG_HTT",
  "SHORT_ATXN1", "LONG_ATXN1",
  "SHORT_ATXN2", "LONG_ATXN2"
)

DT_allele <- drop_na_vars(DT_model, vars_allele)

fit_gene_multinom <- function(df,
                              outcome = "PATHOLOGY",
                              covars = c("SEX", "AGE", "APOE_E4"),
                              short,
                              long) {
  f <- stats::as.formula(
    paste(
      outcome, "~",
      paste(covars, collapse = " + "), "+",
      paste0(
        short, " + I(", short, "^2) + ",
        long, " + I(", long, "^2) + ",
        short, ":", long
      )
    )
  )

  nnet::multinom(
    f,
    data = df,
    trace = FALSE
  )
}

model_multinom_null <- nnet::multinom(
  PATHOLOGY ~ SEX + AGE + APOE_E4,
  data = DT_allele,
  trace = FALSE
)

model_HTT_final <- fit_gene_multinom(
  df = DT_allele,
  short = "SHORT_HTT",
  long = "LONG_HTT"
)

model_ATXN1_final <- fit_gene_multinom(
  df = DT_allele,
  short = "SHORT_ATXN1",
  long = "LONG_ATXN1"
)

model_ATXN2_final <- fit_gene_multinom(
  df = DT_allele,
  short = "SHORT_ATXN2",
  long = "LONG_ATXN2"
)

gene_models <- list(
  HTT = model_HTT_final,
  ATXN1 = model_ATXN1_final,
  ATXN2 = model_ATXN2_final
)

extract_multinom_global <- function(gene, fit, null_fit, n_total) {
  lrt <- lmtest::lrtest(null_fit, fit)
  lrt_df <- as.data.frame(lrt)
  last <- nrow(lrt_df)

  tibble::tibble(
    test_id = paste0("CAG_GLOBAL_MULTINOMIAL__", gene),
    family_id = "D_CAG_MULTINOMIAL_GLOBAL",
    adjustment_method = "holm",
    analysis_scope = "MAIN_DIAGNOSIS",
    outcome = "PATHOLOGY: BD/SCZ vs CONTROL",
    gene = gene,
    predictor_block = "short + short^2 + long + long^2 + short:long",
    covariates = "SEX + AGE + APOE_E4",
    n_total = n_total,
    df_added = lrt_df$Df[last],
    lr_chisq = lrt_df$Chisq[last],
    p_raw = lrt_df$`Pr(>Chisq)`[last],
    AIC_null = stats::AIC(null_fit),
    AIC_full = stats::AIC(fit)
  )
}

multinom_global_results <- purrr::imap_dfr(
  gene_models,
  ~ extract_multinom_global(
    gene = .y,
    fit = .x,
    null_fit = model_multinom_null,
    n_total = nrow(DT_allele)
  )
) %>%
  dplyr::mutate(
    p_adj = stats::p.adjust(p_raw, method = "holm"),
    significant_raw = p_raw < 0.05,
    significant_adjusted = p_adj < 0.05
  ) %>%
  dplyr::arrange(p_adj, p_raw)

multinom_testing_plan <- multinom_global_results %>%
  dplyr::select(
    test_id,
    family_id,
    adjustment_method,
    analysis_scope,
    outcome,
    gene,
    predictor_block,
    covariates
  )

readr::write_csv(
  multinom_testing_plan,
  file.path(
    revision_dir,
    "05_CAG_multinomial_global_testing_plan.csv"
  )
)

readr::write_csv(
  multinom_global_results,
  file.path(
    revision_dir,
    "05_CAG_multinomial_global_results.csv"
  )
)

# -----------------------------------------------------------------------------
# ATXN2 hierarchical component decomposition
# -----------------------------------------------------------------------------
# The global ATXN2 gene-block LRT is the primary inferential result. Because that
# global test is significant after Holm correction, decompose the fitted ATXN2
# block using hierarchical likelihood-ratio tests. These follow-up tests ask which
# parts of the five-term ATXN2 block contribute information beyond the remaining
# terms; they do NOT turn individual polynomial coefficients into simple per-repeat
# effects.
#
# Full ATXN2 block:
#   short + short^2 + long + long^2 + short:long
#
# Follow-up contrasts:
#   1) SHORT contribution:
#      remove short + short^2 + short:long
#      -> tests all information involving the short allele conditional on the
#         long-allele linear/quadratic terms.
#   2) LONG contribution:
#      remove long + long^2 + short:long
#      -> tests all information involving the long allele conditional on the
#         short-allele linear/quadratic terms.
#   3) NON-LINEAR contribution:
#      remove short^2 + long^2 while retaining both linear terms and interaction.
#   4) SHORT^2-only contribution:
#      remove only short^2 while retaining short, long, long^2 and short:long.
#      -> tests short-allele curvature conditional on all remaining terms.
#   5) LONG^2-only contribution:
#      remove only long^2 while retaining short, short^2, long and short:long.
#      -> tests long-allele curvature conditional on all remaining terms.
#   6) SHORT x LONG interaction:
#      remove only short:long while retaining both linear and quadratic terms.
#
# The six component tests form a secondary follow-up family and are Holm-adjusted
# together. They are interpreted only after the primary global ATXN2 test.

model_ATXN2_no_short_block <- nnet::multinom(
  PATHOLOGY ~ SEX + AGE + APOE_E4 +
    LONG_ATXN2 + I(LONG_ATXN2^2),
  data = DT_allele,
  trace = FALSE
)

model_ATXN2_no_long_block <- nnet::multinom(
  PATHOLOGY ~ SEX + AGE + APOE_E4 +
    SHORT_ATXN2 + I(SHORT_ATXN2^2),
  data = DT_allele,
  trace = FALSE
)

model_ATXN2_no_quadratic <- nnet::multinom(
  PATHOLOGY ~ SEX + AGE + APOE_E4 +
    SHORT_ATXN2 + LONG_ATXN2 +
    SHORT_ATXN2:LONG_ATXN2,
  data = DT_allele,
  trace = FALSE
)

model_ATXN2_no_short_quadratic <- nnet::multinom(
  PATHOLOGY ~ SEX + AGE + APOE_E4 +
    SHORT_ATXN2 +
    LONG_ATXN2 + I(LONG_ATXN2^2) +
    SHORT_ATXN2:LONG_ATXN2,
  data = DT_allele,
  trace = FALSE
)

model_ATXN2_no_long_quadratic <- nnet::multinom(
  PATHOLOGY ~ SEX + AGE + APOE_E4 +
    SHORT_ATXN2 + I(SHORT_ATXN2^2) +
    LONG_ATXN2 +
    SHORT_ATXN2:LONG_ATXN2,
  data = DT_allele,
  trace = FALSE
)

model_ATXN2_no_interaction <- nnet::multinom(
  PATHOLOGY ~ SEX + AGE + APOE_E4 +
    SHORT_ATXN2 + I(SHORT_ATXN2^2) +
    LONG_ATXN2 + I(LONG_ATXN2^2),
  data = DT_allele,
  trace = FALSE
)

extract_component_lrt <- function(test_id, component, reduced_fit, full_fit,
                                  terms_removed, interpretation) {
  lrt <- lmtest::lrtest(reduced_fit, full_fit)
  lrt_df <- as.data.frame(lrt)
  last <- nrow(lrt_df)

  tibble::tibble(
    test_id = test_id,
    family_id = "D_ATXN2_COMPONENT_LRT_FOLLOWUP",
    adjustment_method = "holm",
    analysis_scope = "ATXN2_GLOBAL_FOLLOWUP",
    outcome = "PATHOLOGY: BD/SCZ vs CONTROL",
    gene = "ATXN2",
    component = component,
    terms_removed = terms_removed,
    n_total = nrow(DT_allele),
    df_removed = lrt_df$Df[last],
    lr_chisq = lrt_df$Chisq[last],
    p_raw = lrt_df$`Pr(>Chisq)`[last],
    AIC_reduced = stats::AIC(reduced_fit),
    AIC_full = stats::AIC(full_fit),
    interpretation = interpretation
  )
}

atxn2_component_results <- dplyr::bind_rows(
  extract_component_lrt(
    test_id = "ATXN2_COMPONENT__SHORT_BLOCK",
    component = "SHORT block",
    reduced_fit = model_ATXN2_no_short_block,
    full_fit = model_ATXN2_final,
    terms_removed = "short + short^2 + short:long",
    interpretation = paste(
      "Joint contribution of all terms involving the ATXN2 short allele",
      "conditional on the long-allele linear and quadratic terms."
    )
  ),
  extract_component_lrt(
    test_id = "ATXN2_COMPONENT__LONG_BLOCK",
    component = "LONG block",
    reduced_fit = model_ATXN2_no_long_block,
    full_fit = model_ATXN2_final,
    terms_removed = "long + long^2 + short:long",
    interpretation = paste(
      "Joint contribution of all terms involving the ATXN2 long allele",
      "conditional on the short-allele linear and quadratic terms."
    )
  ),
  extract_component_lrt(
    test_id = "ATXN2_COMPONENT__QUADRATIC_TERMS",
    component = "Quadratic terms",
    reduced_fit = model_ATXN2_no_quadratic,
    full_fit = model_ATXN2_final,
    terms_removed = "short^2 + long^2",
    interpretation = paste(
      "Joint evidence for non-linearity from the two quadratic terms",
      "after retaining linear short/long effects and their interaction."
    )
  ),
  extract_component_lrt(
    test_id = "ATXN2_COMPONENT__SHORT_QUADRATIC_ONLY",
    component = "SHORT^2 only",
    reduced_fit = model_ATXN2_no_short_quadratic,
    full_fit = model_ATXN2_final,
    terms_removed = "short^2",
    interpretation = paste(
      "Conditional contribution of ATXN2 short-allele curvature after retaining",
      "the short linear term, the complete long-allele terms and the interaction."
    )
  ),
  extract_component_lrt(
    test_id = "ATXN2_COMPONENT__LONG_QUADRATIC_ONLY",
    component = "LONG^2 only",
    reduced_fit = model_ATXN2_no_long_quadratic,
    full_fit = model_ATXN2_final,
    terms_removed = "long^2",
    interpretation = paste(
      "Conditional contribution of ATXN2 long-allele curvature after retaining",
      "the long linear term, the complete short-allele terms and the interaction."
    )
  ),
  extract_component_lrt(
    test_id = "ATXN2_COMPONENT__INTERACTION",
    component = "SHORT x LONG interaction",
    reduced_fit = model_ATXN2_no_interaction,
    full_fit = model_ATXN2_final,
    terms_removed = "short:long",
    interpretation = paste(
      "Additional contribution of the short-by-long interaction beyond",
      "the linear and quadratic short/long terms."
    )
  )
) %>%
  dplyr::mutate(
    p_adj = stats::p.adjust(p_raw, method = "holm"),
    significant_raw = p_raw < 0.05,
    significant_adjusted = p_adj < 0.05,
    global_ATXN2_p_adj = multinom_global_results$p_adj[
      multinom_global_results$gene == "ATXN2"
    ][1],
    global_ATXN2_significant = multinom_global_results$significant_adjusted[
      multinom_global_results$gene == "ATXN2"
    ][1]
  ) %>%
  dplyr::arrange(p_adj, p_raw)

atxn2_component_testing_plan <- atxn2_component_results %>%
  dplyr::select(
    test_id,
    family_id,
    adjustment_method,
    analysis_scope,
    outcome,
    gene,
    component,
    terms_removed,
    interpretation
  )

readr::write_csv(
  atxn2_component_testing_plan,
  file.path(
    revision_dir,
    "05_ATXN2_component_LRT_testing_plan.csv"
  )
)

readr::write_csv(
  atxn2_component_results,
  file.path(
    revision_dir,
    "05_ATXN2_component_LRT_results.csv"
  )
)

cat("\n============================================================\n")
cat("ATXN2 HIERARCHICAL COMPONENT LIKELIHOOD-RATIO TESTS\n")
cat("Primary global ATXN2 Holm-adjusted p = ",
    format(multinom_global_results$p_adj[
      multinom_global_results$gene == "ATXN2"
    ][1], digits = 5),
    "\n", sep = "")
cat("Six component tests are secondary follow-up tests and are Holm-adjusted together.\n")
cat("============================================================\n")
print(
  atxn2_component_results %>%
    dplyr::select(
      component,
      terms_removed,
      df_removed,
      lr_chisq,
      p_raw,
      p_adj,
      significant_adjusted
    ),
  n = Inf,
  width = Inf
)

# Coefficient-level output is descriptive/model-shape information only.
# nnet/broom versions can encode factor terms either as SEXFemale/APOE_E4E4+
# or as SEX[T.Female]/APOE_E4[T.E4+]. Normalize these labels so exported
# results and downstream figures are stable across package versions.
normalize_multinom_term <- function(x) {
  dplyr::case_when(
    x == "SEXFemale" ~ "SEX[T.Female]",
    x == "APOE_E4E4+" ~ "APOE_E4[T.E4+]",
    TRUE ~ x
  )
}

tidy_multinom_or <- function(fit, gene) {
  broom::tidy(
    fit,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::mutate(
      term = normalize_multinom_term(term),
      gene = gene,
      inference_role = "shape inspection; global gene-block LRT is primary"
    )
}

multinom_coefficients <- purrr::imap_dfr(
  gene_models,
  ~ tidy_multinom_or(.x, .y)
)

readr::write_csv(
  multinom_coefficients,
  file.path(
    revision_dir,
    "05_CAG_multinomial_coefficients.csv"
  )
)

cat("\n============================================================\n")
cat("CONTINUOUS CAG MULTINOMIAL: GLOBAL GENE-BLOCK TESTS\n")
cat("============================================================\n")
print(
  multinom_global_results %>%
    dplyr::select(
      test_id,
      gene,
      n_total,
      df_added,
      lr_chisq,
      p_raw,
      p_adj,
      significant_adjusted
    ),
  n = Inf,
  width = Inf
)

# ============================================================================
# FIGURE 2
# Global CAG gene-block tests + ATXN2 predicted probabilities
# ============================================================================
#
# Panel A:
#   Primary inference from global likelihood-ratio tests comparing each complete
#   CAG gene-block model against the covariate-only multinomial model.
#   Holm adjustment is applied across HTT, ATXN1 and ATXN2.
#
# Panel B:
#   Conditional predicted probabilities across the central 95% of the observed
#   ATXN2 short-allele distribution.
#
# Panel C:
#   Conditional predicted probabilities across the central 95% of the observed
#   ATXN2 long-allele distribution.
#
# IMPORTANT:
#   Individual polynomial coefficients are not interpreted as simple per-repeat
#   effects because linear, quadratic and interaction terms are fitted
#   simultaneously. The global gene-block likelihood-ratio test is the primary
#   inferential result.
# ============================================================================


# ----------------------------------------------------------------------------
# General settings
# ----------------------------------------------------------------------------

fig_dir <- "figures"

if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}


cols_outcome <- c(
  "BD"  = "#8CBDE6",
  "SCZ" = "#F5A04D"
)


# ============================================================================
# PANEL A
# GLOBAL GENE-BLOCK LIKELIHOOD-RATIO TESTS
# ============================================================================


format_p_fig <- function(x) {
  
  dplyr::case_when(
    
    is.na(x) ~ "NA",
    
    x < 0.001 ~
      formatC(
        x,
        format = "e",
        digits = 2
      ),
    
    x < 0.1 ~
      formatC(
        x,
        format = "f",
        digits = 4
      ),
    
    TRUE ~
      formatC(
        x,
        format = "f",
        digits = 3
      )
  )
}


df_global_plot <- multinom_global_results %>%
  
  dplyr::mutate(
    
    gene = factor(
      gene,
      levels = c(
        "HTT",
        "ATXN1",
        "ATXN2"
      )
    ),
    
    significance = factor(
      significant_adjusted,
      levels = c(FALSE, TRUE),
      labels = c(
        "Not significant",
        "Holm-adjusted p < 0.05"
      )
    ),
    
    stat_label = paste0(
      "\u03c7\u00b2 = ",
      sprintf("%.2f", lr_chisq),
      "; p = ",
      format_p_fig(p_raw),
      "; Holm p = ",
      format_p_fig(p_adj)
    )
  )


p_global <- ggplot2::ggplot(
  df_global_plot,
  ggplot2::aes(
    x = lr_chisq,
    y = gene
  )
) +
  
  ggplot2::geom_segment(
    ggplot2::aes(
      x = 0,
      xend = lr_chisq,
      y = gene,
      yend = gene
    ),
    linewidth = 0.7,
    color = "grey75"
  ) +
  
  ggplot2::geom_point(
    ggplot2::aes(
      fill = significance
    ),
    shape = 21,
    size = 4.2,
    stroke = 0.7,
    color = "black"
  ) +
  
  ggplot2::geom_text(
    ggplot2::aes(
      label = stat_label
    ),
    hjust = 0,
    nudge_x = 0.7,
    size = 3.1
  ) +
  
  ggplot2::scale_fill_manual(
    values = c(
      "Not significant" = "white",
      "Holm-adjusted p < 0.05" = "#4C78A8"
    ),
    name = NULL
  ) +
  
  ggplot2::scale_y_discrete(
    labels = c(
      "HTT" = "HTT",
      "ATXN1" = "ATXN1",
      "ATXN2" = "ATXN2"
    )
  ) +
  
  ggplot2::scale_x_continuous(
    limits = c(
      0,
      37
    ),
    breaks = c(
      0,
      10,
      20,
      30
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.01)
    )
  ) +
  
  ggplot2::labs(
    tag = "A",
    x = expression(
      "Likelihood-ratio " *
        chi^2 *
        " (10 df)"
    ),
    y = NULL
  ) +
  
  ggplot2::coord_cartesian(
    clip = "off"
  ) +
  
  ggplot2::theme_classic(
    base_size = 10
  ) +
  
  ggplot2::theme(
    
    legend.position = "top",
    legend.justification = "center",
    
    legend.key.width =
      grid::unit(
        0.6,
        "cm"
      ),
    
    axis.text.y =
      ggplot2::element_text(
        size = 10,
        face = "italic"
      ),
    
    plot.tag =
      ggplot2::element_text(
        face = "bold",
        size = 14
      ),
    
    plot.tag.position =
      c(
        0,
        1
      ),
    
    plot.margin =
      ggplot2::margin(
        5,
        20,
        5,
        5
      )
  )


print(p_global)


# ============================================================================
# PANELS B-C
# ATXN2 CONDITIONAL PREDICTED PROBABILITIES
# ============================================================================


DT_plot <- DT_allele

model_plot <- model_ATXN2_final


# ----------------------------------------------------------------------------
# Sex levels
# ----------------------------------------------------------------------------

sex_levels <- levels(
  DT_plot$SEX
)

if (length(sex_levels) < 2) {
  
  stop(
    "SEX must have at least two levels.",
    call. = FALSE
  )
}


# ----------------------------------------------------------------------------
# Age percentiles
# ----------------------------------------------------------------------------

age_q <- stats::quantile(
  DT_plot$AGE,
  probs = c(
    0.25,
    0.50,
    0.75
  ),
  na.rm = TRUE
)


age_df <- tibble::tibble(
  
  AGE = as.numeric(
    age_q
  ),
  
  age_group = factor(
    c(
      "Age P25",
      "Age P50",
      "Age P75"
    ),
    levels = c(
      "Age P25",
      "Age P50",
      "Age P75"
    )
  )
)


# ----------------------------------------------------------------------------
# Hold non-varying allele at its observed median
# ----------------------------------------------------------------------------

med_short <- stats::median(
  DT_plot$SHORT_ATXN2,
  na.rm = TRUE
)


med_long <- stats::median(
  DT_plot$LONG_ATXN2,
  na.rm = TRUE
)


# ----------------------------------------------------------------------------
# APOE reference category
# ----------------------------------------------------------------------------

if (is.factor(DT_plot$APOE_E4)) {
  
  ref_apoe <- levels(
    DT_plot$APOE_E4
  )[1]
  
} else {
  
  ref_apoe <- "E4-"
}


# ============================================================================
# Restrict plotted predictions to central 95% of observed allele distribution
# ============================================================================


short_limits <- stats::quantile(
  DT_plot$SHORT_ATXN2,
  probs = c(
    0.025,
    0.975
  ),
  na.rm = TRUE
)


long_limits <- stats::quantile(
  DT_plot$LONG_ATXN2,
  probs = c(
    0.025,
    0.975
  ),
  na.rm = TRUE
)


cat(
  "\nFigure 2 ATXN2 short-allele plotting range (P2.5-P97.5): ",
  short_limits[1],
  " - ",
  short_limits[2],
  "\n",
  sep = ""
)


cat(
  "Figure 2 ATXN2 long-allele plotting range (P2.5-P97.5): ",
  long_limits[1],
  " - ",
  long_limits[2],
  "\n",
  sep = ""
)


grid_short <- seq(
  from = as.numeric(
    short_limits[1]
  ),
  to = as.numeric(
    short_limits[2]
  ),
  by = 0.05
)


grid_long <- seq(
  from = as.numeric(
    long_limits[1]
  ),
  to = as.numeric(
    long_limits[2]
  ),
  by = 0.05
)


# ----------------------------------------------------------------------------
# Prediction helper
# ----------------------------------------------------------------------------

predict_multinom_probs <- function(
    model,
    newdata
) {
  
  pr <- as.data.frame(
    stats::predict(
      model,
      newdata = newdata,
      type = "probs"
    )
  )
  
  dplyr::bind_cols(
    newdata,
    pr
  ) %>%
    
    tidyr::pivot_longer(
      cols = dplyr::all_of(
        colnames(pr)
      ),
      names_to = "outcome",
      values_to = "p"
    )
}


# ============================================================================
# PANEL B
# ATXN2 SHORT ALLELE
# ============================================================================


nd_short <- tidyr::expand_grid(
  
  SHORT_ATXN2 = grid_short,
  
  SEX = factor(
    sex_levels,
    levels = sex_levels
  ),
  
  age_df,
  
  APOE_E4 = ref_apoe
  
) %>%
  
  dplyr::mutate(
    LONG_ATXN2 = med_long
  )


pred_short <- predict_multinom_probs(
  model_plot,
  nd_short
) %>%
  
  dplyr::filter(
    outcome %in% c(
      "BD",
      "SCZ"
    )
  ) %>%
  
  dplyr::mutate(
    outcome = factor(
      outcome,
      levels = c(
        "BD",
        "SCZ"
      )
    )
  )


p_short <- ggplot2::ggplot(
  pred_short,
  ggplot2::aes(
    x = SHORT_ATXN2,
    y = p,
    color = outcome,
    linetype = SEX
  )
) +
  
  ggplot2::geom_line(
    linewidth = 0.8
  ) +
  
  ggplot2::facet_wrap(
    ~ age_group,
    nrow = 1
  ) +
  
  ggplot2::scale_color_manual(
    values = cols_outcome,
    name = "Diagnosis"
  ) +
  
  ggplot2::scale_linetype_manual(
    values = c(
      "solid",
      "dashed"
    )[seq_along(sex_levels)],
    name = "Sex"
  ) +
  
  ggplot2::scale_x_continuous(
    limits = as.numeric(
      short_limits
    ),
    breaks = seq(
      ceiling(short_limits[1]),
      floor(short_limits[2]),
      by = 1
    )
  ) +
  
  ggplot2::scale_y_continuous(
    limits = c(
      0,
      1
    ),
    breaks = c(
      0,
      0.25,
      0.50,
      0.75,
      1
    )
  ) +
  
  ggplot2::labs(
    tag = "B",
    x = expression(
      italic("ATXN2") ~
        "short allele (CAG repeats)"
    ),
    y = "Predicted probability"
  ) +
  
  ggplot2::theme_classic(
    base_size = 10
  ) +
  
  ggplot2::theme(
    
    legend.position = "top",
    legend.box = "horizontal",
    
    legend.margin =
      ggplot2::margin(
        0,
        0,
        2,
        0
      ),
    
    legend.key.width =
      grid::unit(
        1,
        "cm"
      ),
    
    strip.background =
      ggplot2::element_rect(
        fill = "grey95",
        color = NA
      ),
    
    strip.text =
      ggplot2::element_text(
        face = "bold",
        size = 9
      ),
    
    plot.tag =
      ggplot2::element_text(
        face = "bold",
        size = 14
      ),
    
    plot.tag.position =
      c(
        0,
        1
      )
  )


print(p_short)


# ============================================================================
# PANEL C
# ATXN2 LONG ALLELE
# ============================================================================


nd_long <- tidyr::expand_grid(
  
  LONG_ATXN2 = grid_long,
  
  SEX = factor(
    sex_levels,
    levels = sex_levels
  ),
  
  age_df,
  
  APOE_E4 = ref_apoe
  
) %>%
  
  dplyr::mutate(
    SHORT_ATXN2 = med_short
  )


pred_long <- predict_multinom_probs(
  model_plot,
  nd_long
) %>%
  
  dplyr::filter(
    outcome %in% c(
      "BD",
      "SCZ"
    )
  ) %>%
  
  dplyr::mutate(
    outcome = factor(
      outcome,
      levels = c(
        "BD",
        "SCZ"
      )
    )
  )


p_long <- ggplot2::ggplot(
  pred_long,
  ggplot2::aes(
    x = LONG_ATXN2,
    y = p,
    color = outcome,
    linetype = SEX
  )
) +
  
  ggplot2::geom_line(
    linewidth = 0.8
  ) +
  
  ggplot2::facet_wrap(
    ~ age_group,
    nrow = 1
  ) +
  
  ggplot2::scale_color_manual(
    values = cols_outcome,
    name = "Diagnosis"
  ) +
  
  ggplot2::scale_linetype_manual(
    values = c(
      "solid",
      "dashed"
    )[seq_along(sex_levels)],
    name = "Sex"
  ) +
  
  ggplot2::scale_x_continuous(
    limits = as.numeric(
      long_limits
    ),
    breaks = seq(
      ceiling(long_limits[1]),
      floor(long_limits[2]),
      by = 1
    )
  ) +
  
  ggplot2::scale_y_continuous(
    limits = c(
      0,
      1
    ),
    breaks = c(
      0,
      0.25,
      0.50,
      0.75,
      1
    )
  ) +
  
  ggplot2::labs(
    tag = "C",
    x = expression(
      italic("ATXN2") ~
        "long allele (CAG repeats)"
    ),
    y = "Predicted probability"
  ) +
  
  ggplot2::theme_classic(
    base_size = 10
  ) +
  
  ggplot2::theme(
    
    # Same encoding as B, so legend does not need to be repeated
    legend.position = "none",
    
    strip.background =
      ggplot2::element_rect(
        fill = "grey95",
        color = NA
      ),
    
    strip.text =
      ggplot2::element_text(
        face = "bold",
        size = 9
      ),
    
    plot.tag =
      ggplot2::element_text(
        face = "bold",
        size = 14
      ),
    
    plot.tag.position =
      c(
        0,
        1
      )
  )


print(p_long)


# ============================================================================
# FINAL FIGURE 2
# Vertical publication layout
# ============================================================================


Figure_2 <-
  
  p_global /
  
  p_short /
  
  p_long +
  
  patchwork::plot_layout(
    heights = c(
      0.75,
      1,
      1
    )
  )


print(Figure_2)


# ============================================================================
# EXPORT INDIVIDUAL PANELS
# ============================================================================


ggplot2::ggsave(
  filename = file.path(
    fig_dir,
    "05_CAG_global_gene_block_tests.tiff"
  ),
  plot = p_global,
  device = "tiff",
  width = 180,
  height = 60,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)


ggplot2::ggsave(
  filename = file.path(
    fig_dir,
    "05_ATXN2_predicted_probabilities_short.tiff"
  ),
  plot = p_short,
  device = "tiff",
  width = 180,
  height = 80,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)


ggplot2::ggsave(
  filename = file.path(
    fig_dir,
    "05_ATXN2_predicted_probabilities_long.tiff"
  ),
  plot = p_long,
  device = "tiff",
  width = 180,
  height = 80,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)


# ============================================================================
# EXPORT FINAL FIGURE 2
# ============================================================================


ggplot2::ggsave(
  filename = file.path(
    fig_dir,
    "Figure_2_multinomial_ATXN2.tiff"
  ),
  plot = Figure_2,
  device = "tiff",
  width = 180,
  height = 230,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)


ggplot2::ggsave(
  filename = file.path(
    fig_dir,
    "Figure_2_multinomial_ATXN2.png"
  ),
  plot = Figure_2,
  device = "png",
  width = 180,
  height = 230,
  units = "mm",
  dpi = 600
)

# 2. LEGACY MODEL-SELECTION ANALYSES REMOVED ================================
# The previous block-AIC/greedy binomial workflow is intentionally omitted
# from the reviewer revision. Clinical subgroup inference is handled in the
# prespecified testing frameworks in 03_Genotype_stats.R and
# 06_Survival_age_analysis.R.

# Session info ----
sessionInfo()