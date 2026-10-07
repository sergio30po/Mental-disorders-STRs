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
DT <- DT %>%
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
  DT <- derive_short_long(DT, g)
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

DT_allele <- drop_na_vars(DT, vars_allele)

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

# Coefficient-level output is descriptive/model-shape information only.
tidy_multinom_or <- function(fit, gene) {
  broom::tidy(
    fit,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::mutate(
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

# FIGURE 2.1 ATXN2 multinomial figure ----
# --- Settings
fig_dir <- "figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Ensure model and tidy data exist
model_plot <- model_ATXN2_final
tt_plot <- broom::tidy(model_plot, exponentiate = TRUE, conf.int = TRUE)

# --- Settings Colors & Limits
cols_outcome <- c("BD" = "#8CBDE6", "SCZ" = "#F5A04D")
x_min <- 1e-3
x_max <- 1e4

# --- Data Preparation
df_fp <- tt_plot %>%
  filter(term != "(Intercept)") %>%
  mutate(
    outcome = factor(str_trim(as.character(y.level)), levels = c("BD", "SCZ")),

    # MODIFICATION: Using atop() to split labels into two lines
    term_label = case_when(
      term == "SEX[T.Female]" ~ "'Female sex'",
      term == "AGE" ~ "'Age'",

      term == "APOE_E4[T.E4+]" ~ "italic(APOE)~epsilon[4]~carrier",

      # Split ATXN2 labels:
      term == "SHORT_ATXN2" ~ "atop(italic(ATXN2), 'short allele (linear)')",
      stringr::str_detect(term, "^I\\(SHORT_ATXN2\\^2\\)") ~ "atop(italic(ATXN2), 'short allele (quadratic)')",

      term == "LONG_ATXN2" ~ "atop(italic(ATXN2), 'long allele (linear)')",
      stringr::str_detect(term, "^I\\(LONG_ATXN2\\^2\\)") ~ "atop(italic(ATXN2), 'long allele (quadratic)')",

      term == "SHORT_ATXN2:LONG_ATXN2" ~ "atop(italic(ATXN2), 'short x long allele (interaction)')",
      TRUE ~ term
    ),

    # Stats processing
    p_plot = pmax(p.value, 1e-300),
    sig = -log10(p_plot),
    est_p = pmin(pmax(estimate,  x_min), x_max),
    lo_p  = pmin(pmax(conf.low,  x_min), x_max),
    hi_p  = pmin(pmax(conf.high, x_min), x_max),
    cut_left  = conf.low  < x_min,
    cut_right = conf.high > x_max
  )

# --- Order terms (Must match the strings in case_when EXACTLY)
order_terms <- c(
  "'Female sex'",
  "'Age'",
  "italic(APOE)~epsilon[4]~carrier",
  "atop(italic(ATXN2), 'short allele (linear)')",
  "atop(italic(ATXN2), 'short allele (quadratic)')",
  "atop(italic(ATXN2), 'long allele (linear)')",
  "atop(italic(ATXN2), 'long allele (quadratic)')",
  "atop(italic(ATXN2), 'short x long allele (interaction)')"
)

df_fp <- df_fp %>%
  mutate(term_label = factor(term_label, levels = rev(order_terms)))

# --- Plotting
pd <- position_dodge(width = 0.55)

g_forest <- ggplot(df_fp, aes(x = est_p, y = term_label)) +
  coord_cartesian(xlim = c(x_min, x_max), clip = "off") +

  geom_vline(xintercept = 1, linetype = "dashed",
             linewidth = 0.6, color = "grey45") +

  # Error bars
  geom_errorbarh(
    aes(xmin = lo_p, xmax = hi_p, color = outcome),
    position = pd, height = 0.18, linewidth = 0.9
  ) +

  # Points
  geom_point(
    aes(fill = outcome, size = sig),
    shape = 21, color = "black", stroke = 0.35,
    position = pd, alpha = 0.95
  ) +

  # CI truncation arrows (Left)
  geom_segment(
    data = df_fp %>% filter(cut_left),
    aes(x = x_min * 1.35, xend = x_min * 1.08,
        y = term_label, yend = term_label, color = outcome),
    inherit.aes = FALSE,
    arrow = arrow(type = "closed", length = unit(0.16, "cm")),
    linewidth = 0.9
  ) +

  # CI truncation arrows (Right)
  geom_segment(
    data = df_fp %>% filter(cut_right),
    aes(x = x_max / 1.35, xend = x_max / 1.08,
        y = term_label, yend = term_label, color = outcome),
    inherit.aes = FALSE,
    arrow = arrow(type = "closed", length = unit(0.16, "cm")),
    linewidth = 0.9
  ) +

  scale_x_log10(
    limits = c(x_min, x_max),
    expand = expansion(mult = c(0.05, 0)),
    name = "Odds ratio (log scale)"
  ) +

  # parse=TRUE renders the atop() logic
  scale_y_discrete(labels = function(x) parse(text = x)) +

  scale_fill_manual(
    values = cols_outcome,
    name = "Diagnosis"
  ) +
  scale_color_manual(
    values = cols_outcome,
    name = "Diagnosis"
  ) +

  scale_size_continuous(
    name = expression(atop(-log[10](italic(p)), "(with 95% CI)")),
    range = c(2.6, 6.8)
  ) +

  guides(
    color = "none",
    fill  = guide_legend(order = 1, override.aes = list(size = 5)),
    size  = guide_legend(order = 2)
  ) +

  theme_classic(base_size = 12) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(lineheight = 0.8),
    legend.position = c(0.99, 0.99),
    legend.justification = c(1, 1),
    legend.box = "horizontal",
    legend.box.just = "top",
    legend.spacing.x = unit(0.4, "cm"),
    legend.background = element_rect(fill = "white", colour = "grey70", linewidth = 0.4),
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9),

    plot.margin = margin(8, 20, 8, 8)
  )

print(g_forest)

ggsave(
  filename = file.path(fig_dir, "ATXN2_forest_plot.tiff"),
  plot = g_forest,
  device = "tiff",
  width = 500, height = 160, units = "mm",
  dpi = 600, compression = "lzw"
)
# FIGURE 2.2/3: Predicted probabilities vs allele size----
# --- Model & Data Setup
DT_plot <- DT_allele
model_plot <- model_ATXN2_final

# 1. Define SEX levels
sex_levels <- levels(DT_plot$SEX)
if (length(sex_levels) < 2) stop("SEX must have at least 2 levels.")

# 2. Define AGE quantiles (P25, P50, P75)
age_q <- quantile(DT_plot$AGE, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
age_df <- tibble(
  AGE = as.numeric(age_q),
  age_group = factor(c("Age P25", "Age P50", "Age P75"),
                     levels = c("Age P25", "Age P50", "Age P75"))
)

# 3. Define allele medians (to hold the non-varying allele constant)
med_short <- median(DT_plot$SHORT_ATXN2, na.rm = TRUE)
med_long  <- median(DT_plot$LONG_ATXN2, na.rm = TRUE)

# 4. Define APOE reference
if(is.factor(DT_plot$APOE_E4)) {
  ref_apoe <- levels(DT_plot$APOE_E4)[1]
} else {
  ref_apoe <- "E4-"
}

# 5. Prediction Grids
grid_short <- seq(min(DT_plot$SHORT_ATXN2, na.rm = TRUE),
                  max(DT_plot$SHORT_ATXN2, na.rm = TRUE), by = 0.05)
grid_long  <- seq(min(DT_plot$LONG_ATXN2, na.rm = TRUE),
                  max(DT_plot$LONG_ATXN2, na.rm = TRUE), by = 0.05)

# --- Helper Function
predict_probs_long <- function(model, newdata) {
  pr <- as.data.frame(predict(model, newdata = newdata, type = "probs"))
  bind_cols(newdata, pr) %>%
    pivot_longer(
      cols = all_of(colnames(pr)),
      names_to = "outcome",
      values_to = "p"
    )
}

# FIGURE 2.2: VARY LONG ALLELE ----

# Create newdata
nd_long <- expand_grid(
  LONG_ATXN2 = grid_long,
  SEX = factor(sex_levels, levels = sex_levels),
  age_df,
  APOE_E4 = ref_apoe
) %>%
  mutate(SHORT_ATXN2 = med_short)

# Predict
pred_long <- predict_probs_long(model_plot, nd_long) %>%
  filter(outcome %in% c("BD", "SCZ")) %>%
  mutate(outcome = factor(outcome, levels = c("BD", "SCZ")))

# Plot
p_long <- ggplot(pred_long, aes(x = LONG_ATXN2, y = p, color = outcome, linetype = SEX)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ age_group, nrow = 1) +

  scale_color_manual(values = cols_outcome, name = "Diagnosis") +
  scale_linetype_manual(values = c("solid", "dashed")[seq_along(sex_levels)], name = "Sex") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +

  labs(
    title = NULL,
    subtitle = NULL,
    x = expression(italic("ATXN2") ~ "long allele (CAG repeats)"),
    y = "Predicted probability"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    legend.key.width = unit(1.5, "cm"),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold")
  )

print(p_long)

ggsave(
  filename = file.path(fig_dir, "Predicted_probs_long_byAge.tiff"),
  plot = p_long,
  device = "tiff",
  width = 10, height = 5, units = "in",
  dpi = 600, compression = "lzw"
)


# FIGURE 2.3: VARY SHORT ALLELE ----

# Create newdata
nd_short <- expand_grid(
  SHORT_ATXN2 = grid_short,
  SEX = factor(sex_levels, levels = sex_levels),
  age_df,
  APOE_E4 = ref_apoe
) %>%
  mutate(LONG_ATXN2 = med_long)

# Predict
pred_short <- predict_probs_long(model_plot, nd_short) %>%
  filter(outcome %in% c("BD", "SCZ")) %>%
  mutate(outcome = factor(outcome, levels = c("BD", "SCZ")))

# Plot
p_short <- ggplot(pred_short, aes(x = SHORT_ATXN2, y = p, color = outcome, linetype = SEX)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ age_group, nrow = 1) +

  scale_color_manual(values = cols_outcome, name = "Diagnosis") +
  scale_linetype_manual(values = c("solid", "dashed")[seq_along(sex_levels)], name = "Sex") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +

  labs(
    title = NULL,
    subtitle = NULL,
    x = expression(italic("ATXN2") ~ "short allele (CAG repeats)"),
    y = "Predicted probability"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    legend.key.width = unit(1.5, "cm"),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold")
  )

print(p_short)

ggsave(
  filename = file.path(fig_dir, "Predicted_probs_short_byAge.tiff"),
  plot = p_short,
  device = "tiff",
  width = 10, height = 5, units = "in",
  dpi = 600, compression = "lzw"
)
# FIGURE 2: Layout Configuration ----
Figure_Composite <-
  wrap_elements(g_forest) /
  (p_short | p_long) +

  plot_layout(heights = c(1.3, 1)) +

  plot_annotation(tag_levels = "A") &

  theme(
    plot.tag = element_text(face = "bold", size = 20),
    plot.tag.position = c(0, 1),
    plot.tag.padding = unit(5, "pt")
  )

# Preview
print(Figure_Composite)

# --- Save to File ---
ggsave(
  filename = file.path(fig_dir, "Figure_Composite_ATXN2.tiff"),
  plot = Figure_Composite,
  device = "tiff",
  width = 450, height = 400, units = "mm",
  dpi = 600, compression = "lzw"
)

# 2. LEGACY MODEL-SELECTION ANALYSES REMOVED ================================
# The previous block-AIC/greedy binomial workflow is intentionally omitted
# from the reviewer revision. Clinical subgroup inference is handled in the
# prespecified testing frameworks in 03_Genotype_stats.R and
# 06_Survival_age_analysis.R.

# Session info ----
sessionInfo()