# Script name: 04_CAG_repeat_sizes.R
# ==============================================================================
# Title: Continuous CAG repeat-size analyses and allele-distribution figures.
#
# Author: Sergio Pérez Oliveira
#
# Reviewer-driven revision:
#   1) Treat continuous repeat length as a sensitivity analysis to transferred
#      categorical IA thresholds.
#   2) Define multiplicity families explicitly instead of applying local Holm
#      corrections independently within each call.
#   3) Report short and long allele distributions separately.
#   4) Add a two-dimensional short × long allele representation for each gene.
#
# Main inferential orientation:
#   - "positive" is the first diagnostic/subgroup level in each contrast.
#   - rank-biserial > 0 means larger CAG repeat values in positive vs reference.
#
# Multiplicity families:
#   A0_CAG_DIAGNOSIS_OMNIBUS:
#       3 genes × 2 allele positions = 6 Kruskal-Wallis tests; Holm.
#   A1_CAG_DIAGNOSIS_PAIRWISE:
#       3 genes × 2 allele positions × 3 main diagnostic contrasts = 18 tests;
#       Holm.
#   B0/B1_CAG_SUBTYPE:
#       BD and SCZ subtype analyses; exploratory BH-FDR.
#   C0/C1_CAG_COGNITIVE:
#       cognitive-status analyses within BD/SCZ; exploratory BH-FDR.
#
# Inputs:
#   - Manually selected environment file with custom functions/data import.
#   - Dataframes created by the environment: DT, BD_CONTROLS, SCZ_CONTROLS.
#
# Outputs:
#   - Testing plans and master result tables in results/reviewer_revision/
#   - Allele-order QC table
#   - Separate short- and long-allele frequency tables
#   - Short-allele histogram, long-allele histogram and short×long 2D plot
#     for HTT, ATXN1 and ATXN2
# ==============================================================================

# Load environment ----
Env_path <- file.choose()
source(Env_path)
rm(Env_path)

revision_dir <- file.path("results", "reviewer_revision")
if (!dir.exists(revision_dir)) {
  dir.create(revision_dir, recursive = TRUE)
}

fig_dir <- "figures"
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

# ==============================================================================
# 0. DEFINITIONS AND QC
# ==============================================================================

# ALLELE1/ALLELE2 are raw diploid allele fields and are NOT assumed to be
# consistently ordered. Reviewer-facing "short" and "long" alleles are therefore
# derived within each participant as min(raw allele 1, raw allele 2) and
# max(raw allele 1, raw allele 2), respectively.
gene_map <- tibble::tribble(
  ~gene,   ~raw1_var,       ~raw2_var,       ~short_var,     ~long_var,
  "HTT",   "ALLELE1_HTT",   "ALLELE2_HTT",   "SHORT_HTT",    "LONG_HTT",
  "ATXN1", "ALLELE1_ATXN1", "ALLELE2_ATXN1", "SHORT_ATXN1",  "LONG_ATXN1",
  "ATXN2", "ALLELE1_ATXN2", "ALLELE2_ATXN2", "SHORT_ATXN2",  "LONG_ATXN2"
)

required_raw <- c(
  "PATHOLOGY",
  gene_map$raw1_var,
  gene_map$raw2_var
)

missing_dt <- setdiff(required_raw, names(DT))
if (length(missing_dt) > 0) {
  stop(
    "DT is missing required raw allele column(s): ",
    paste(missing_dt, collapse = ", "),
    call. = FALSE
  )
}

derive_ordered_alleles <- function(df) {

  for (i in seq_len(nrow(gene_map))) {

    raw1 <- gene_map$raw1_var[[i]]
    raw2 <- gene_map$raw2_var[[i]]
    short_name <- gene_map$short_var[[i]]
    long_name <- gene_map$long_var[[i]]

    if (!all(c(raw1, raw2) %in% names(df))) {
      stop(
        "Cannot derive ordered alleles: missing ",
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

# Use the same within-person ordering rule in every dataset used below.
DT <- derive_ordered_alleles(DT)
BD_CONTROLS <- derive_ordered_alleles(BD_CONTROLS)
SCZ_CONTROLS <- derive_ordered_alleles(SCZ_CONTROLS)

# QC documents how often the raw ALLELE1/ALLELE2 order had to be swapped.
# A non-zero n_raw1_gt_raw2 is not an error; it is precisely why derived
# short/long variables are needed.
allele_order_qc <- purrr::pmap_dfr(
  gene_map,
  function(gene, raw1_var, raw2_var, short_var, long_var) {

    d <- DT %>%
      dplyr::transmute(
        raw1 = as.numeric(.data[[raw1_var]]),
        raw2 = as.numeric(.data[[raw2_var]]),
        short = as.numeric(.data[[short_var]]),
        long = as.numeric(.data[[long_var]])
      ) %>%
      dplyr::filter(!is.na(raw1), !is.na(raw2))

    tibble::tibble(
      gene = gene,
      n_complete = nrow(d),
      n_raw1_gt_raw2 = sum(d$raw1 > d$raw2),
      pct_raw1_gt_raw2 = 100 * n_raw1_gt_raw2 / n_complete,
      n_raw1_eq_raw2 = sum(d$raw1 == d$raw2),
      n_derived_short_gt_long = sum(d$short > d$long),
      min_short = min(d$short),
      max_short = max(d$short),
      min_long = min(d$long),
      max_long = max(d$long)
    )
  }
)

readr::write_csv(
  allele_order_qc,
  file.path(revision_dir, "04_allele_order_QC.csv")
)

if (any(allele_order_qc$n_derived_short_gt_long > 0)) {
  stop(
    "Derived short/long QC failed unexpectedly.",
    call. = FALSE
  )
}

cat("\n============================================================\n")
cat("ALLELE ORDER QC AND WITHIN-PERSON REORDERING\n")
cat("============================================================\n")
print(allele_order_qc, n = Inf, width = Inf)

# ==============================================================================
# 1. STATISTICAL HELPERS
# ==============================================================================

run_kw_test <- function(data,
                        group_col,
                        value_col,
                        allowed_groups,
                        analysis_scope,
                        gene,
                        allele_position,
                        family_id,
                        adjustment_method,
                        test_id) {

  d <- data %>%
    dplyr::transmute(
      group = as.character(.data[[group_col]]),
      value = as.numeric(.data[[value_col]])
    ) %>%
    dplyr::filter(
      !is.na(group),
      !is.na(value),
      is.finite(value),
      group %in% allowed_groups
    )

  d$group <- factor(d$group, levels = allowed_groups)
  d <- droplevels(d)

  k <- nlevels(d$group)
  n <- nrow(d)

  if (k < 2L || n <= k) {
    return(
      tibble::tibble(
        test_id = test_id,
        family_id = family_id,
        adjustment_method = adjustment_method,
        analysis_scope = analysis_scope,
        gene = gene,
        allele_position = allele_position,
        value_variable = value_col,
        group_variable = group_col,
        groups = paste(allowed_groups, collapse = " | "),
        n_total = n,
        n_groups = k,
        statistic = NA_real_,
        df = NA_real_,
        epsilon_squared = NA_real_,
        p_raw = NA_real_,
        estimable = FALSE
      )
    )
  }

  kw <- stats::kruskal.test(value ~ group, data = d)

  # Epsilon-squared for Kruskal-Wallis:
  # epsilon^2_H = (H - k + 1) / (n - k)
  eps2 <- (unname(kw$statistic) - k + 1) / (n - k)
  eps2 <- max(0, eps2)

  tibble::tibble(
    test_id = test_id,
    family_id = family_id,
    adjustment_method = adjustment_method,
    analysis_scope = analysis_scope,
    gene = gene,
    allele_position = allele_position,
    value_variable = value_col,
    group_variable = group_col,
    groups = paste(allowed_groups, collapse = " | "),
    n_total = n,
    n_groups = k,
    statistic = unname(kw$statistic),
    df = unname(kw$parameter),
    epsilon_squared = eps2,
    p_raw = kw$p.value,
    estimable = TRUE
  )
}


run_pairwise_wilcox <- function(data,
                                group_col,
                                value_col,
                                positive,
                                reference,
                                analysis_scope,
                                gene,
                                allele_position,
                                family_id,
                                adjustment_method,
                                test_id) {

  d <- data %>%
    dplyr::transmute(
      group = as.character(.data[[group_col]]),
      value = as.numeric(.data[[value_col]])
    ) %>%
    dplyr::filter(
      !is.na(group),
      !is.na(value),
      is.finite(value),
      group %in% c(positive, reference)
    )

  x_pos <- d$value[d$group == positive]
  x_ref <- d$value[d$group == reference]

  n_pos <- length(x_pos)
  n_ref <- length(x_ref)

  if (n_pos < 2L || n_ref < 2L) {
    return(
      tibble::tibble(
        test_id = test_id,
        family_id = family_id,
        adjustment_method = adjustment_method,
        analysis_scope = analysis_scope,
        gene = gene,
        allele_position = allele_position,
        value_variable = value_col,
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

  # wilcox.test(x, y) reports the Mann-Whitney U statistic on a 0..n1*n2 scale.
  u <- as.numeric(wt$statistic)
  rank_biserial <- 2 * u / (n_pos * n_ref) - 1

  tibble::tibble(
    test_id = test_id,
    family_id = family_id,
    adjustment_method = adjustment_method,
    analysis_scope = analysis_scope,
    gene = gene,
    allele_position = allele_position,
    value_variable = value_col,
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
    rank_biserial = rank_biserial,
    p_raw = wt$p.value,
    estimable = TRUE
  )
}


adjust_family <- function(df) {
  df %>%
    dplyr::group_by(family_id) %>%
    dplyr::mutate(
      p_adj = stats::p.adjust(
        p_raw,
        method = dplyr::first(adjustment_method)
      ),
      significant_raw = !is.na(p_raw) & p_raw < 0.05,
      significant_adjusted = !is.na(p_adj) & p_adj < 0.05
    ) %>%
    dplyr::ungroup()
}

# ==============================================================================
# 2. TESTING PLANS
# ==============================================================================

allele_plan <- gene_map %>%
  tidyr::pivot_longer(
    cols = c(short_var, long_var),
    names_to = "allele_key",
    values_to = "value_col"
  ) %>%
  dplyr::mutate(
    allele_position = dplyr::if_else(
      allele_key == "short_var",
      "SHORT",
      "LONG"
    )
  ) %>%
  dplyr::select(gene, allele_position, value_col)

# ---- A. Main diagnosis --------------------------------------------------------

main_omnibus_plan <- allele_plan %>%
  dplyr::mutate(
    dataset_name = "DT",
    group_col = "PATHOLOGY",
    allowed_groups = purrr::map(
      gene,
      ~ c("BD", "SCZ", "CONTROL")
    ),
    analysis_scope = "MAIN_DIAGNOSIS",
    family_id = "A0_CAG_DIAGNOSIS_OMNIBUS",
    adjustment_method = "holm",
    test_id = paste(
      "CAG_OMNIBUS", gene, allele_position, analysis_scope,
      sep = "__"
    )
  )

main_pairwise_contrasts <- tibble::tribble(
  ~positive, ~reference,
  "BD",      "SCZ",
  "BD",      "CONTROL",
  "SCZ",     "CONTROL"
)

main_pairwise_plan <- tidyr::crossing(
  allele_plan,
  main_pairwise_contrasts
) %>%
  dplyr::mutate(
    dataset_name = "DT",
    group_col = "PATHOLOGY",
    analysis_scope = "MAIN_DIAGNOSIS",
    family_id = "A1_CAG_DIAGNOSIS_PAIRWISE",
    adjustment_method = "holm",
    test_id = paste(
      "CAG", gene, allele_position, analysis_scope,
      positive, "vs", reference,
      sep = "__"
    )
  )

# ---- B. Diagnostic subtypes ---------------------------------------------------

subtype_contexts <- list(
  BD_SUBTYPE = list(
    data = BD_CONTROLS,
    group_col = "PATHOLOGY_TYPE_BINARY",
    groups = c("BD-I", "Other", "CONTROL"),
    contrasts = list(
      c("BD-I", "Other"),
      c("BD-I", "CONTROL"),
      c("Other", "CONTROL")
    )
  ),
  SCZ_SUBTYPE = list(
    data = SCZ_CONTROLS,
    group_col = "PATHOLOGY_TYPE_BINARY",
    groups = c("SCZ", "Other", "CONTROL"),
    contrasts = list(
      c("SCZ", "Other"),
      c("SCZ", "CONTROL"),
      c("Other", "CONTROL")
    )
  )
)

build_omnibus_context_plan <- function(scope, spec, family_id) {
  allele_plan %>%
    dplyr::mutate(
      dataset_name = scope,
      group_col = spec$group_col,
      allowed_groups = purrr::map(gene, ~ spec$groups),
      analysis_scope = scope,
      family_id = family_id,
      adjustment_method = "BH",
      test_id = paste(
        "CAG_OMNIBUS", gene, allele_position, scope,
        sep = "__"
      )
    )
}

build_pairwise_context_plan <- function(scope, spec, family_id) {

  contrasts <- purrr::map_dfr(
    spec$contrasts,
    ~ tibble::tibble(
      positive = .x[[1]],
      reference = .x[[2]]
    )
  )

  tidyr::crossing(
    allele_plan,
    contrasts
  ) %>%
    dplyr::mutate(
      dataset_name = scope,
      group_col = spec$group_col,
      analysis_scope = scope,
      family_id = family_id,
      adjustment_method = "BH",
      test_id = paste(
        "CAG", gene, allele_position, scope,
        positive, "vs", reference,
        sep = "__"
      )
    )
}

subtype_omnibus_plan <- dplyr::bind_rows(
  purrr::imap(
    subtype_contexts,
    ~ build_omnibus_context_plan(
      scope = .y,
      spec = .x,
      family_id = "B0_CAG_SUBTYPE_OMNIBUS"
    )
  )
)

subtype_pairwise_plan <- dplyr::bind_rows(
  purrr::imap(
    subtype_contexts,
    ~ build_pairwise_context_plan(
      scope = .y,
      spec = .x,
      family_id = "B1_CAG_SUBTYPE_PAIRWISE"
    )
  )
)

# ---- C. Cognitive status ------------------------------------------------------

cognitive_contexts <- list(
  BD_COGNITIVE = list(
    data = BD_CONTROLS,
    group_col = "CD_BINARY",
    groups = c("CD", "No-CD", "CONTROL"),
    contrasts = list(
      c("CD", "No-CD"),
      c("CD", "CONTROL"),
      c("No-CD", "CONTROL")
    )
  ),
  SCZ_COGNITIVE = list(
    data = SCZ_CONTROLS,
    group_col = "CD_BINARY",
    groups = c("CD", "No-CD", "CONTROL"),
    contrasts = list(
      c("CD", "No-CD"),
      c("CD", "CONTROL"),
      c("No-CD", "CONTROL")
    )
  )
)

cognitive_omnibus_plan <- dplyr::bind_rows(
  purrr::imap(
    cognitive_contexts,
    ~ build_omnibus_context_plan(
      scope = .y,
      spec = .x,
      family_id = "C0_CAG_COGNITIVE_OMNIBUS"
    )
  )
)

cognitive_pairwise_plan <- dplyr::bind_rows(
  purrr::imap(
    cognitive_contexts,
    ~ build_pairwise_context_plan(
      scope = .y,
      spec = .x,
      family_id = "C1_CAG_COGNITIVE_PAIRWISE"
    )
  )
)

omnibus_plan <- dplyr::bind_rows(
  main_omnibus_plan,
  subtype_omnibus_plan,
  cognitive_omnibus_plan
)

pairwise_plan <- dplyr::bind_rows(
  main_pairwise_plan,
  subtype_pairwise_plan,
  cognitive_pairwise_plan
)

# CSV-friendly versions of testing plans
readr::write_csv(
  omnibus_plan %>%
    dplyr::mutate(
      allowed_groups = purrr::map_chr(
        allowed_groups,
        ~ paste(.x, collapse = " | ")
      )
    ),
  file.path(revision_dir, "04_CAG_omnibus_testing_plan.csv")
)

readr::write_csv(
  pairwise_plan,
  file.path(revision_dir, "04_CAG_pairwise_testing_plan.csv")
)

# ==============================================================================
# 3. RUN TESTING PLANS
# ==============================================================================

get_analysis_data <- function(dataset_name) {
  if (dataset_name == "DT") {
    return(DT)
  }
  if (dataset_name == "BD_SUBTYPE" || dataset_name == "BD_COGNITIVE") {
    return(BD_CONTROLS)
  }
  if (dataset_name == "SCZ_SUBTYPE" || dataset_name == "SCZ_COGNITIVE") {
    return(SCZ_CONTROLS)
  }

  stop("Unknown dataset_name: ", dataset_name, call. = FALSE)
}

omnibus_results_raw <- purrr::pmap_dfr(
  omnibus_plan,
  function(gene,
           allele_position,
           value_col,
           dataset_name,
           group_col,
           allowed_groups,
           analysis_scope,
           family_id,
           adjustment_method,
           test_id) {

    run_kw_test(
      data = get_analysis_data(dataset_name),
      group_col = group_col,
      value_col = value_col,
      allowed_groups = allowed_groups,
      analysis_scope = analysis_scope,
      gene = gene,
      allele_position = allele_position,
      family_id = family_id,
      adjustment_method = adjustment_method,
      test_id = test_id
    )
  }
)

omnibus_results <- omnibus_results_raw %>%
  adjust_family() %>%
  dplyr::arrange(family_id, p_adj, p_raw)

pairwise_results_raw <- purrr::pmap_dfr(
  pairwise_plan,
  function(gene,
           allele_position,
           value_col,
           positive,
           reference,
           dataset_name,
           group_col,
           analysis_scope,
           family_id,
           adjustment_method,
           test_id) {

    run_pairwise_wilcox(
      data = get_analysis_data(dataset_name),
      group_col = group_col,
      value_col = value_col,
      positive = positive,
      reference = reference,
      analysis_scope = analysis_scope,
      gene = gene,
      allele_position = allele_position,
      family_id = family_id,
      adjustment_method = adjustment_method,
      test_id = test_id
    )
  }
)

pairwise_results <- pairwise_results_raw %>%
  adjust_family() %>%
  dplyr::arrange(family_id, p_adj, p_raw)

readr::write_csv(
  omnibus_results,
  file.path(revision_dir, "04_CAG_omnibus_master_results.csv")
)

readr::write_csv(
  pairwise_results,
  file.path(revision_dir, "04_CAG_pairwise_master_results.csv")
)

# ==============================================================================
# 4. REVIEWER-FACING STATISTICAL OUTPUT
# ==============================================================================

cat("\n============================================================\n")
cat("MAIN CONTINUOUS CAG TESTS SIGNIFICANT AFTER HOLM\n")
cat("============================================================\n")

print(
  pairwise_results %>%
    dplyr::filter(
      family_id == "A1_CAG_DIAGNOSIS_PAIRWISE",
      significant_adjusted
    ) %>%
    dplyr::select(
      test_id, gene, allele_position, contrast,
      n_positive, n_reference,
      median_positive, median_reference,
      median_difference, rank_biserial,
      p_raw, p_adj
    ),
  n = Inf,
  width = Inf
)

cat("\n============================================================\n")
cat("MAIN CONTINUOUS CAG: TEN SMALLEST RAW P-VALUES\n")
cat("============================================================\n")

print(
  pairwise_results %>%
    dplyr::filter(
      family_id == "A1_CAG_DIAGNOSIS_PAIRWISE"
    ) %>%
    dplyr::select(
      test_id, gene, allele_position, contrast,
      n_positive, n_reference,
      median_positive, median_reference,
      median_difference, rank_biserial,
      p_raw, p_adj
    ) %>%
    dplyr::arrange(p_raw) %>%
    dplyr::slice_head(n = 10),
  n = Inf,
  width = Inf
)

cat("\n============================================================\n")
cat("HTT CONTINUOUS SENSITIVITY: BD-I vs CONTROL\n")
cat("============================================================\n")

print(
  pairwise_results %>%
    dplyr::filter(
      family_id == "B1_CAG_SUBTYPE_PAIRWISE",
      analysis_scope == "BD_SUBTYPE",
      gene == "HTT",
      positive == "BD-I",
      reference == "CONTROL"
    ) %>%
    dplyr::select(
      test_id, allele_position,
      n_positive, n_reference,
      median_positive, median_reference,
      median_difference, rank_biserial,
      p_raw, p_adj
    ),
  n = Inf,
  width = Inf
)

cat("\n============================================================\n")
cat("ALL EXPLORATORY CONTINUOUS CAG TESTS SIGNIFICANT AFTER BH-FDR\n")
cat("============================================================\n")

print(
  pairwise_results %>%
    dplyr::filter(
      family_id %in% c(
        "B1_CAG_SUBTYPE_PAIRWISE",
        "C1_CAG_COGNITIVE_PAIRWISE"
      ),
      significant_adjusted
    ) %>%
    dplyr::select(
      test_id, family_id, gene, allele_position,
      analysis_scope, contrast,
      n_positive, n_reference,
      median_difference, rank_biserial,
      p_raw, p_adj
    ),
  n = Inf,
  width = Inf
)

# ==============================================================================
# 5. SEPARATE SHORT/LONG FREQUENCY TABLES AND 2D DISTRIBUTIONS
# ==============================================================================

diagnosis_levels <- c("CONTROL", "BD", "SCZ")

DT_plot <- DT %>%
  dplyr::filter(as.character(PATHOLOGY) %in% diagnosis_levels) %>%
  dplyr::mutate(
    PATHOLOGY = factor(
      as.character(PATHOLOGY),
      levels = diagnosis_levels
    )
  )

make_gene_distribution_outputs <- function(gene, short_var, long_var) {

  d <- DT_plot %>%
    dplyr::transmute(
      PATHOLOGY = PATHOLOGY,
      short = as.numeric(.data[[short_var]]),
      long = as.numeric(.data[[long_var]])
    ) %>%
    dplyr::filter(
      !is.na(PATHOLOGY),
      !is.na(short),
      !is.na(long)
    )

  # Separate frequency tables for short and long alleles
  freq_short <- d %>%
    dplyr::count(PATHOLOGY, short, name = "count") %>%
    dplyr::group_by(PATHOLOGY) %>%
    dplyr::mutate(
      percentage = 100 * count / sum(count),
      allele_position = "SHORT"
    ) %>%
    dplyr::ungroup() %>%
    dplyr::rename(CAG_size = short)

  freq_long <- d %>%
    dplyr::count(PATHOLOGY, long, name = "count") %>%
    dplyr::group_by(PATHOLOGY) %>%
    dplyr::mutate(
      percentage = 100 * count / sum(count),
      allele_position = "LONG"
    ) %>%
    dplyr::ungroup() %>%
    dplyr::rename(CAG_size = long)

  freq_1d <- dplyr::bind_rows(freq_short, freq_long) %>%
    dplyr::mutate(gene = gene, .before = 1)

  readr::write_csv(
    freq_1d,
    file.path(
      revision_dir,
      paste0("04_", gene, "_short_long_frequencies.csv")
    )
  )

  # 2D genotype distribution: percentage of participants within diagnosis group
  freq_2d <- d %>%
    dplyr::count(PATHOLOGY, short, long, name = "count") %>%
    dplyr::group_by(PATHOLOGY) %>%
    dplyr::mutate(
      percentage = 100 * count / sum(count)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(gene = gene, .before = 1)

  readr::write_csv(
    freq_2d,
    file.path(
      revision_dir,
      paste0("04_", gene, "_short_by_long_2D_frequencies.csv")
    )
  )

  # Histograms use percentages because diagnostic group sizes differ markedly.
  p_short <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = short)
  ) +
    ggplot2::geom_histogram(
      ggplot2::aes(
        y = 100 * ggplot2::after_stat(count / sum(count))
      ),
      binwidth = 1,
      boundary = 0.5,
      color = "black",
      fill = "grey75"
    ) +
    ggplot2::facet_wrap(
      ~ PATHOLOGY,
      ncol = 1,
      scales = "free_y"
    ) +
    ggplot2::labs(
      x = paste0(gene, " short allele (CAG repeats)"),
      y = "Participants (%)"
    ) +
    ggplot2::theme_classic(base_size = 12)

  p_long <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = long)
  ) +
    ggplot2::geom_histogram(
      ggplot2::aes(
        y = 100 * ggplot2::after_stat(count / sum(count))
      ),
      binwidth = 1,
      boundary = 0.5,
      color = "black",
      fill = "grey75"
    ) +
    ggplot2::facet_wrap(
      ~ PATHOLOGY,
      ncol = 1,
      scales = "free_y"
    ) +
    ggplot2::labs(
      x = paste0(gene, " long allele (CAG repeats)"),
      y = "Participants (%)"
    ) +
    ggplot2::theme_classic(base_size = 12)

  # Integer short×long combinations are summarized before plotting, avoiding
  # overplotting. Point area represents the within-diagnosis percentage.
  p_2d <- ggplot2::ggplot(
    freq_2d,
    ggplot2::aes(
      x = short,
      y = long,
      size = percentage
    )
  ) +
    ggplot2::geom_point(
      shape = 21,
      fill = "grey75",
      color = "black",
      alpha = 0.8
    ) +
    ggplot2::facet_wrap(
      ~ PATHOLOGY,
      ncol = 1
    ) +
    ggplot2::scale_size_area(
      max_size = 8,
      name = "Participants (%)"
    ) +
    ggplot2::labs(
      x = paste0(gene, " short allele (CAG repeats)"),
      y = paste0(gene, " long allele (CAG repeats)")
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "right"
    )

  ggplot2::ggsave(
    filename = file.path(
      fig_dir,
      paste0("04_", gene, "_short_allele_histogram.tiff")
    ),
    plot = p_short,
    device = "tiff",
    width = 180,
    height = 240,
    units = "mm",
    dpi = 600,
    compression = "lzw"
  )

  ggplot2::ggsave(
    filename = file.path(
      fig_dir,
      paste0("04_", gene, "_long_allele_histogram.tiff")
    ),
    plot = p_long,
    device = "tiff",
    width = 180,
    height = 240,
    units = "mm",
    dpi = 600,
    compression = "lzw"
  )

  ggplot2::ggsave(
    filename = file.path(
      fig_dir,
      paste0("04_", gene, "_short_by_long_2D.tiff")
    ),
    plot = p_2d,
    device = "tiff",
    width = 190,
    height = 250,
    units = "mm",
    dpi = 600,
    compression = "lzw"
  )

  invisible(
    list(
      frequency_1d = freq_1d,
      frequency_2d = freq_2d,
      short_plot = p_short,
      long_plot = p_long,
      plot_2d = p_2d
    )
  )
}

distribution_outputs <- purrr::pmap(
  gene_map %>%
    dplyr::select(gene, short_var, long_var),
  make_gene_distribution_outputs
)

names(distribution_outputs) <- gene_map$gene

# Combined CSVs for convenient downstream use
all_1d_frequencies <- purrr::map_dfr(
  distribution_outputs,
  "frequency_1d"
)

all_2d_frequencies <- purrr::map_dfr(
  distribution_outputs,
  "frequency_2d"
)

readr::write_csv(
  all_1d_frequencies,
  file.path(revision_dir, "04_all_genes_short_long_frequencies.csv")
)

readr::write_csv(
  all_2d_frequencies,
  file.path(revision_dir, "04_all_genes_short_by_long_2D_frequencies.csv")
)

cat("\n============================================================\n")
cat("04 COMPLETE\n")
cat("============================================================\n")
cat("Continuous CAG testing tables written to results/reviewer_revision/.\n")
cat("Separate short/long histograms and short x long 2D figures written to figures/.\n")

# Session info ----
sessionInfo()
