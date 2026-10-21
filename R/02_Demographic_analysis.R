# Script name: 02_Demographic_analysis.R
# ==============================================================================
# Title: Descriptive analysis of psychiatric cohorts (SCZ, BD, Controls).
#
# Author: Sergio Pérez Oliveira
#
# Purpose:
#   Descriptive demographic/clinical comparisons used in Table 1, with explicit
#   and reproducible comparison/reference orientation.
#
# Diagnostic reference convention:
#   - CONTROL is always the reference whenever controls are included.
#   - For the direct BD/SCZ comparison, BD is the reference and SCZ is the
#     comparison group.
#
# Effect orientation:
#   - Continuous outcomes: rank-biserial > 0 means larger values in the
#     comparison group than in the reference group.
#   - Binary outcomes: OR = odds(event in comparison) / odds(event in reference).
#
# Table 1 binary events:
#   - SEX: Male
#   - SMOKER: Smoking
#   - COFFEE: Coffee
#
# Outputs:
#   - Console descriptive summaries.
#   - results/reviewer_revision/02_Table1_continuous_results.csv
#   - results/reviewer_revision/02_Table1_categorical_results.csv
# ==============================================================================

# Load environment ----
Env_path <- file.choose()
source(Env_path)
rm(Env_path)

revision_dir <- file.path("results", "reviewer_revision")
if (!dir.exists(revision_dir)) {
  dir.create(revision_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# Helpers: explicitly oriented Table 1 comparisons
# ------------------------------------------------------------------------------

run_oriented_wilcox <- function(data,
                                group_col,
                                value_col,
                                comparison,
                                reference,
                                test_id,
                                family_id,
                                adjustment_method = "none") {

  d <- data |>
    dplyr::transmute(
      group = as.character(.data[[group_col]]),
      value = suppressWarnings(as.numeric(.data[[value_col]]))
    ) |>
    dplyr::filter(
      group %in% c(reference, comparison),
      !is.na(value)
    )

  x_cmp <- d$value[d$group == comparison]
  x_ref <- d$value[d$group == reference]

  if (length(x_cmp) == 0L || length(x_ref) == 0L) {
    stop("Missing comparison/reference data for ", test_id, call. = FALSE)
  }

  wt <- stats::wilcox.test(
    x_cmp,
    x_ref,
    exact = FALSE,
    correct = FALSE
  )

  # wilcox.test(x, y) reports the Mann-Whitney U statistic for x
  # on the 0..n_x*n_y scale.
  # r_rb > 0 = larger values in comparison than reference.
  n_cmp <- length(x_cmp)
  n_ref <- length(x_ref)
  u_cmp <- unname(wt$statistic)
  r_rb <- 2 * u_cmp / (n_cmp * n_ref) - 1

  tibble::tibble(
    test_id = test_id,
    family_id = family_id,
    adjustment_method = adjustment_method,
    variable = value_col,
    comparison = comparison,
    reference = reference,
    contrast = paste0(comparison, " vs ", reference),
    effect_definition = paste0(
      "rank-biserial > 0 = larger ", value_col,
      " in comparison vs reference"
    ),
    n_comparison = n_cmp,
    n_reference = n_ref,
    mean_comparison = mean(x_cmp),
    sd_comparison = stats::sd(x_cmp),
    mean_reference = mean(x_ref),
    sd_reference = stats::sd(x_ref),
    rank_biserial = r_rb,
    p_raw = wt$p.value
  )
}

run_oriented_fisher_event <- function(data,
                                      group_col,
                                      outcome_col,
                                      event,
                                      comparison,
                                      reference,
                                      test_id,
                                      family_id,
                                      adjustment_method = "none") {

  d <- data |>
    dplyr::transmute(
      group = as.character(.data[[group_col]]),
      outcome = as.character(.data[[outcome_col]])
    ) |>
    dplyr::filter(
      group %in% c(reference, comparison),
      !is.na(outcome)
    )

  if (!all(c(reference, comparison) %in% unique(d$group))) {
    stop("Missing comparison/reference group for ", test_id, call. = FALSE)
  }

  if (!event %in% unique(d$outcome)) {
    stop("Event '", event, "' not found for ", test_id, call. = FALSE)
  }

  a <- sum(d$group == comparison & d$outcome == event)
  b <- sum(d$group == comparison & d$outcome != event)
  c0 <- sum(d$group == reference & d$outcome == event)
  d0 <- sum(d$group == reference & d$outcome != event)

  tab <- matrix(
    c(a, b, c0, d0),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      group = c(comparison, reference),
      outcome = c(event, "Other")
    )
  )

  ft <- stats::fisher.test(tab, conf.int = TRUE)

  tibble::tibble(
    test_id = test_id,
    family_id = family_id,
    adjustment_method = adjustment_method,
    variable = outcome_col,
    event = event,
    comparison = comparison,
    reference = reference,
    contrast = paste0(comparison, " vs ", reference),
    effect_definition = paste0(
      "OR = odds(", event, " in comparison) / odds(",
      event, " in reference)"
    ),
    n_comparison = a + b,
    n_reference = c0 + d0,
    event_comparison = a,
    event_reference = c0,
    event_pct_comparison = 100 * a / (a + b),
    event_pct_reference = 100 * c0 / (c0 + d0),
    odds_ratio = unname(ft$estimate),
    conf_low = unname(ft$conf.int[1]),
    conf_high = unname(ft$conf.int[2]),
    p_raw = ft$p.value
  )
}

apply_family_adjustment <- function(df) {
  df |>
    dplyr::group_by(family_id) |>
    dplyr::group_modify(
      ~ {
        method <- dplyr::first(.x$adjustment_method)
        .x$p_adj <- if (identical(method, "none")) {
          .x$p_raw
        } else {
          stats::p.adjust(.x$p_raw, method = method)
        }
        .x
      }
    ) |>
    dplyr::ungroup()
}

# ------------------------------------------------------------------------------
# Descriptive summaries
# ------------------------------------------------------------------------------

mean_sd(MENTAL, "PATHOLOGY", "ONSET_AGE", "Age at psychiatric onset")
mean_sd(DT, "PATHOLOGY", "AGE", "Age at evaluation / last visit")
mean_sd(MENTAL, "PATHOLOGY", "DURATION", "Observed disease duration")

# Omnibus current-age test retained as descriptive context.
print(stats::kruskal.test(AGE ~ PATHOLOGY, data = DT))
print(effectsize::epsilon_squared(AGE ~ PATHOLOGY, data = DT))

# ------------------------------------------------------------------------------
# Table 1 continuous comparisons
# ------------------------------------------------------------------------------

continuous_results <- dplyr::bind_rows(
  run_oriented_wilcox(
    MENTAL, "PATHOLOGY", "ONSET_AGE",
    comparison = "SCZ", reference = "BD",
    test_id = "TABLE1__ONSET_AGE__SCZ_vs_BD",
    family_id = "TABLE1_ONSET_AGE",
    adjustment_method = "none"
  ),
  run_oriented_wilcox(
    DT, "PATHOLOGY", "AGE",
    comparison = "BD", reference = "CONTROL",
    test_id = "TABLE1__AGE__BD_vs_CONTROL",
    family_id = "TABLE1_AGE_DIAGNOSIS",
    adjustment_method = "holm"
  ),
  run_oriented_wilcox(
    DT, "PATHOLOGY", "AGE",
    comparison = "SCZ", reference = "CONTROL",
    test_id = "TABLE1__AGE__SCZ_vs_CONTROL",
    family_id = "TABLE1_AGE_DIAGNOSIS",
    adjustment_method = "holm"
  ),
  run_oriented_wilcox(
    DT, "PATHOLOGY", "AGE",
    comparison = "SCZ", reference = "BD",
    test_id = "TABLE1__AGE__SCZ_vs_BD",
    family_id = "TABLE1_AGE_DIAGNOSIS",
    adjustment_method = "holm"
  ),
  run_oriented_wilcox(
    MENTAL, "PATHOLOGY", "DURATION",
    comparison = "SCZ", reference = "BD",
    test_id = "TABLE1__DURATION__SCZ_vs_BD",
    family_id = "TABLE1_DURATION",
    adjustment_method = "none"
  )
) |>
  apply_family_adjustment()

readr::write_csv(
  continuous_results,
  file.path(revision_dir, "02_Table1_continuous_results.csv")
)

# ------------------------------------------------------------------------------
# Table 1 categorical comparisons
# ------------------------------------------------------------------------------

categorical_results <- dplyr::bind_rows(
  # SEX: Male event; three diagnostic pairwise comparisons, Holm-adjusted.
  run_oriented_fisher_event(
    DT, "PATHOLOGY", "SEX", "Male",
    comparison = "BD", reference = "CONTROL",
    test_id = "TABLE1__SEX_MALE__BD_vs_CONTROL",
    family_id = "TABLE1_SEX_DIAGNOSIS",
    adjustment_method = "holm"
  ),
  run_oriented_fisher_event(
    DT, "PATHOLOGY", "SEX", "Male",
    comparison = "SCZ", reference = "CONTROL",
    test_id = "TABLE1__SEX_MALE__SCZ_vs_CONTROL",
    family_id = "TABLE1_SEX_DIAGNOSIS",
    adjustment_method = "holm"
  ),
  run_oriented_fisher_event(
    DT, "PATHOLOGY", "SEX", "Male",
    comparison = "SCZ", reference = "BD",
    test_id = "TABLE1__SEX_MALE__SCZ_vs_BD",
    family_id = "TABLE1_SEX_DIAGNOSIS",
    adjustment_method = "holm"
  ),

  # Lifestyle variables are available only in the psychiatric cohorts.
  run_oriented_fisher_event(
    MENTAL, "PATHOLOGY", "SMOKER", "Smoking",
    comparison = "SCZ", reference = "BD",
    test_id = "TABLE1__SMOKING__SCZ_vs_BD",
    family_id = "TABLE1_SMOKING",
    adjustment_method = "none"
  ),
  run_oriented_fisher_event(
    MENTAL, "PATHOLOGY", "COFFEE", "Coffee",
    comparison = "SCZ", reference = "BD",
    test_id = "TABLE1__COFFEE__SCZ_vs_BD",
    family_id = "TABLE1_COFFEE",
    adjustment_method = "none"
  )
) |>
  apply_family_adjustment()

readr::write_csv(
  categorical_results,
  file.path(revision_dir, "02_Table1_categorical_results.csv")
)

cat("\n============================================================\n")
cat("TABLE 1 CONTINUOUS RESULTS — EXPLICIT ORIENTATION\n")
cat("CONTROL reference when present; BD reference for SCZ vs BD.\n")
cat("============================================================\n")
print(continuous_results, n = Inf, width = Inf)

cat("\n============================================================\n")
cat("TABLE 1 CATEGORICAL RESULTS — EXPLICIT ORIENTATION\n")
cat("OR = odds(event in comparison) / odds(event in reference).\n")
cat("============================================================\n")
print(categorical_results, n = Inf, width = Inf)

# ------------------------------------------------------------------------------
# Cognitive status and subtype descriptive summaries
# ------------------------------------------------------------------------------

cd_tab <- table(BD$CD)
cd_percent <- round(prop.table(cd_tab) * 100, 2)
print(cd_tab)
print(cd_percent)

cd_bin_tab <- table(BD$CD_BINARY)
print(round(prop.table(cd_bin_tab) * 100, 2))
print(cd_bin_tab)

cd_tab <- table(SCZ$CD)
cd_percent <- round(prop.table(cd_tab) * 100, 2)
print(cd_tab)
print(cd_percent)

cd_bin_tab <- table(SCZ$CD_BINARY)
print(round(prop.table(cd_bin_tab) * 100, 2))
print(cd_bin_tab)

subtype_tab <- table(BD$PATHOLOGY_TYPE)
print(round(prop.table(subtype_tab) * 100, 2))
print(subtype_tab)

subtype_bin_tab <- table(BD$PATHOLOGY_TYPE_BINARY)
print(round(prop.table(subtype_bin_tab) * 100, 2))
print(subtype_bin_tab)

subtype_TBP_tab <- table(BD$BD_type)
print(round(prop.table(subtype_TBP_tab) * 100, 2))
print(subtype_TBP_tab)

SCZ_type_tab <- table(SCZ$PATHOLOGY_TYPE)
print(round(prop.table(SCZ_type_tab) * 100, 2))
print(SCZ_type_tab)

SCZ_type_bin_tab <- table(SCZ$PATHOLOGY_TYPE_BINARY)
print(round(prop.table(SCZ_type_bin_tab) * 100, 2))
print(SCZ_type_bin_tab)

# Missing-data audit ------------------------------------------------------------

print(
  DT %>%
    dplyr::group_by(PATHOLOGY) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(is.na(.))))
)

print(
  MENTAL %>%
    dplyr::group_by(PATHOLOGY) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(is.na(.))))
)

# Session info ------------------------------------------------------------------
sessionInfo()
