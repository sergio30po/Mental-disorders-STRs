# Script name: 04_CAG_repeat_sizes.R
# ==============================================================================
# Title: Continuous CAG repeat-size analyses and allele-distribution figures.
#
# Author: Sergio Pérez Oliveira
#
# Analysis framework:
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
#   - Testing plans and master result tables in the results directory
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
# consistently ordered. "Short" and "long" alleles are therefore
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
  "SCZ",     "BD",
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
# 4. KEY STATISTICAL OUTPUT
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
# ============================================================================== 
# 5. SHORT/LONG ALLELE FREQUENCIES AND JOINT 2D DISTRIBUTIONS
# ============================================================================== 

diagnosis_levels <- c("CONTROL", "BD", "SCZ")

diagnosis_labels <- c(
  "CONTROL" = "Controls",
  "BD" = "BD",
  "SCZ" = "SCZ"
)

diagnosis_colors <- c(
  "CONTROL" = "#BDBDBD",
  "BD" = "#4C78A8",
  "SCZ" = "#F28E2B"
)

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
      !is.na(long),
      is.finite(short),
      is.finite(long)
    )

  # Separate short- and long-allele frequency distributions.
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

  # Joint short x long distribution: percentage of participants within diagnosis.
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

  invisible(
    list(
      data = d,
      frequency_1d = freq_1d,
      frequency_2d = freq_2d
    )
  )
}

distribution_outputs <- purrr::pmap(
  gene_map %>%
    dplyr::select(gene, short_var, long_var),
  make_gene_distribution_outputs
)

names(distribution_outputs) <- gene_map$gene

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

# Frequency QC: every within-diagnosis distribution should sum to 100%.
frequency_qc_1d <- all_1d_frequencies %>%
  dplyr::group_by(gene, PATHOLOGY, allele_position) %>%
  dplyr::summarise(
    total_percentage = sum(percentage, na.rm = TRUE),
    .groups = "drop"
  )

frequency_qc_2d <- all_2d_frequencies %>%
  dplyr::group_by(gene, PATHOLOGY) %>%
  dplyr::summarise(
    total_percentage = sum(percentage, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  frequency_qc_1d,
  file.path(revision_dir, "04_allele_frequency_1D_QC.csv")
)

readr::write_csv(
  frequency_qc_2d,
  file.path(revision_dir, "04_allele_frequency_2D_QC.csv")
)

if (
  any(abs(frequency_qc_1d$total_percentage - 100) > 1e-8) ||
  any(abs(frequency_qc_2d$total_percentage - 100) > 1e-8)
) {
  stop(
    "Allele-frequency QC failed: at least one distribution does not sum to 100%.",
    call. = FALSE
  )
}

# ============================================================================== 
# 6. SUPPLEMENTARY FIGURE 1
# ============================================================================== 

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop(
    "Package 'patchwork' is required to compose Supplementary Figure 1.",
    call. = FALSE
  )
}

# For each gene, the figure contains the short-allele distribution, the
# long-allele distribution, and the joint short x long distribution. Percentages
# are calculated within each diagnostic group so that distributions are directly
# comparable despite unequal group sizes.

make_1d_frequency_panel <- function(gene_name, allele_name) {

  d <- all_1d_frequencies %>%
    dplyr::filter(
      gene == gene_name,
      allele_position == allele_name
    ) %>%
    dplyr::select(
      PATHOLOGY,
      CAG_size,
      percentage
    )

  full_sizes <- seq(
    min(d$CAG_size, na.rm = TRUE),
    max(d$CAG_size, na.rm = TRUE),
    by = 1
  )

  d <- d %>%
    tidyr::complete(
      PATHOLOGY = factor(diagnosis_levels, levels = diagnosis_levels),
      CAG_size = full_sizes,
      fill = list(percentage = 0)
    ) %>%
    dplyr::mutate(
      PATHOLOGY = factor(PATHOLOGY, levels = diagnosis_levels)
    )

  allele_label <- if (allele_name == "SHORT") {
    "Short allele"
  } else {
    "Long allele"
  }

  ggplot2::ggplot(
    d,
    ggplot2::aes(
      x = CAG_size,
      y = percentage,
      fill = PATHOLOGY
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge2(
        width = 0.90,
        preserve = "single",
        padding = 0.08
      ),
      width = 0.84,
      color = "black",
      linewidth = 0.18
    ) +
    ggplot2::scale_fill_manual(
      values = diagnosis_colors,
      breaks = diagnosis_levels,
      labels = diagnosis_labels,
      name = "Diagnosis",
      drop = FALSE,
      guide = ggplot2::guide_legend(order = 1)
    ) +
    ggplot2::scale_x_continuous(
      breaks = function(x) {
        br <- pretty(x, n = 8)
        br[br >= min(x) & br <= max(x)]
      },
      expand = ggplot2::expansion(mult = c(0.01, 0.02))
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.06))
    ) +
    ggplot2::labs(
      title = bquote(italic(.(gene_name)) ~ "-" ~ .(allele_label)),
      x = "CAG repeat length",
      y = "Participants (%)"
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5,
        size = 11
      ),
      axis.title = ggplot2::element_text(size = 10),
      axis.text = ggplot2::element_text(size = 8),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 8),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    )
}

# Use a common point-size scale across all 2D panels so bubble areas are directly
# comparable between genes and diagnostic groups. The same size limits and breaks
# are therefore applied to HTT, ATXN1 and ATXN2.
global_2d_max <- max(
  all_2d_frequencies$percentage,
  na.rm = TRUE
)

size_break_candidates <- c(1, 5, 10, 25, 50, 75)
size_breaks_2d <- size_break_candidates[
  size_break_candidates <= global_2d_max
]

if (length(size_breaks_2d) == 0) {
  size_breaks_2d <- pretty(c(0, global_2d_max), n = 4)
  size_breaks_2d <- size_breaks_2d[size_breaks_2d > 0]
}

make_2d_frequency_panel <- function(gene_name) {

  d <- all_2d_frequencies %>%
    dplyr::filter(gene == gene_name) %>%
    dplyr::mutate(
      PATHOLOGY = factor(PATHOLOGY, levels = diagnosis_levels)
    )

  x_range <- range(d$short, na.rm = TRUE)
  y_range <- range(d$long, na.rm = TRUE)

  ggplot2::ggplot(
    d,
    ggplot2::aes(
      x = short,
      y = long,
      size = percentage,
      fill = PATHOLOGY
    )
  ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      linewidth = 0.3,
      color = "grey55"
    ) +
    ggplot2::geom_point(
      shape = 21,
      color = "black",
      stroke = 0.22,
      alpha = 0.82
    ) +
    ggplot2::facet_wrap(
      ~ PATHOLOGY,
      nrow = 1,
      labeller = ggplot2::as_labeller(diagnosis_labels)
    ) +
    ggplot2::scale_fill_manual(
      values = diagnosis_colors,
      breaks = diagnosis_levels,
      labels = diagnosis_labels,
      name = "Diagnosis",
      drop = FALSE,
      guide = "none"
    ) +
    ggplot2::scale_size_area(
      max_size = 7.5,
      limits = c(0, global_2d_max),
      breaks = size_breaks_2d,
      name = "Participants (%)",
      guide = ggplot2::guide_legend(order = 2)
    ) +
    ggplot2::scale_x_continuous(
      limits = x_range,
      breaks = function(x) {
        br <- pretty(x, n = 5)
        br[br >= min(x) & br <= max(x)]
      },
      expand = ggplot2::expansion(mult = c(0.03, 0.03))
    ) +
    ggplot2::scale_y_continuous(
      limits = y_range,
      breaks = function(x) {
        br <- pretty(x, n = 5)
        br[br >= min(x) & br <= max(x)]
      },
      expand = ggplot2::expansion(mult = c(0.03, 0.03))
    ) +
    # Short and long alleles are measured in the same units (CAG repeats), so
    # equal scaling preserves the geometry of the short = long diagonal.
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = bquote(italic(.(gene_name)) ~ "-" ~ "Short × long"),
      x = "Short allele (CAG repeats)",
      y = "Long allele (CAG repeats)"
    ) +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5,
        size = 11
      ),
      axis.title = ggplot2::element_text(size = 9),
      axis.text = ggplot2::element_text(size = 7),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = 8
      ),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 8),
      panel.spacing = grid::unit(0.55, "lines"),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    )
}

# Build the nine panels.
p_HTT_short <- make_1d_frequency_panel("HTT", "SHORT")
p_HTT_long  <- make_1d_frequency_panel("HTT", "LONG")
p_HTT_2d    <- make_2d_frequency_panel("HTT")

p_ATXN1_short <- make_1d_frequency_panel("ATXN1", "SHORT")
p_ATXN1_long  <- make_1d_frequency_panel("ATXN1", "LONG")
p_ATXN1_2d    <- make_2d_frequency_panel("ATXN1")

p_ATXN2_short <- make_1d_frequency_panel("ATXN2", "SHORT")
p_ATXN2_long  <- make_1d_frequency_panel("ATXN2", "LONG")
p_ATXN2_2d    <- make_2d_frequency_panel("ATXN2")

# Save individual panels for inspection and reuse.
individual_panels <- list(
  "A_HTT_short" = p_HTT_short,
  "B_HTT_long" = p_HTT_long,
  "C_HTT_2D" = p_HTT_2d,
  "D_ATXN1_short" = p_ATXN1_short,
  "E_ATXN1_long" = p_ATXN1_long,
  "F_ATXN1_2D" = p_ATXN1_2d,
  "G_ATXN2_short" = p_ATXN2_short,
  "H_ATXN2_long" = p_ATXN2_long,
  "I_ATXN2_2D" = p_ATXN2_2d
)

purrr::iwalk(
  individual_panels,
  function(plot_object, panel_name) {
    ggplot2::ggsave(
      filename = file.path(
        fig_dir,
        paste0("Supplementary_Figure_1_", panel_name, ".tiff")
      ),
      plot = plot_object,
      device = "tiff",
      width = if (grepl("2D$", panel_name)) 240 else 150,
      height = 120,
      units = "mm",
      dpi = 600,
      compression = "lzw"
    )
  }
)

# Compose the complete supplementary figure. The 2D column is wider because each
# panel contains three diagnosis-specific facets. Diagnosis colors are explained
# once in the global legend from the 1D panels; 2D panels contribute only the
# shared Participants (%) bubble-size legend.
# Add panel letters explicitly to each ggplot so their typography is controlled
# at panel level and cannot be reset during guide collection.
add_panel_tag <- function(plot_object, panel_tag) {
  plot_object +
    ggplot2::labs(tag = panel_tag) +
    ggplot2::theme(
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = 16,
        colour = "black"
      )
    )
}

p_HTT_short_tagged  <- add_panel_tag(p_HTT_short, "A")
p_HTT_long_tagged   <- add_panel_tag(p_HTT_long, "B")
p_HTT_2d_tagged     <- add_panel_tag(p_HTT_2d, "C")
p_ATXN1_short_tagged <- add_panel_tag(p_ATXN1_short, "D")
p_ATXN1_long_tagged  <- add_panel_tag(p_ATXN1_long, "E")
p_ATXN1_2d_tagged    <- add_panel_tag(p_ATXN1_2d, "F")
p_ATXN2_short_tagged <- add_panel_tag(p_ATXN2_short, "G")
p_ATXN2_long_tagged  <- add_panel_tag(p_ATXN2_long, "H")
p_ATXN2_2d_tagged    <- add_panel_tag(p_ATXN2_2d, "I")

Supplementary_Figure_1 <- patchwork::wrap_plots(
  p_HTT_short_tagged,
  p_HTT_long_tagged,
  p_HTT_2d_tagged,
  p_ATXN1_short_tagged,
  p_ATXN1_long_tagged,
  p_ATXN1_2d_tagged,
  p_ATXN2_short_tagged,
  p_ATXN2_long_tagged,
  p_ATXN2_2d_tagged,
  ncol = 3,
  widths = c(1, 1, 1.65),
  guides = "collect"
)

Supplementary_Figure_1 <- Supplementary_Figure_1 &
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.tag = ggplot2::element_text(
      face = "bold",
      size = 16,
      colour = "black"
    )
  )

print(Supplementary_Figure_1)

# High-resolution submission figure.
ggplot2::ggsave(
  filename = file.path(
    fig_dir,
    "Supplementary_Figure_1.tiff"
  ),
  plot = Supplementary_Figure_1,
  device = "tiff",
  width = 380,
  height = 285,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

ggplot2::ggsave(
  filename = file.path(
    fig_dir,
    "Supplementary_Figure_1.png"
  ),
  plot = Supplementary_Figure_1,
  device = "png",
  width = 380,
  height = 285,
  units = "mm",
  dpi = 600,
)

# Vector version for inspection and future editing.
ggplot2::ggsave(
  filename = file.path(
    fig_dir,
    "Supplementary_Figure_1.pdf"
  ),
  plot = Supplementary_Figure_1,
  device = grDevices::cairo_pdf,
  width = 380 / 25.4,
  height = 285 / 25.4,
  units = "in"
)

cat("\n============================================================\n")
cat("04 COMPLETE\n")
cat("============================================================\n")
cat("Continuous CAG testing tables written to: ", revision_dir, "\n", sep = "")
cat("Allele-distribution frequency tables written to: ", revision_dir, "\n", sep = "")
cat("Supplementary Figure 1 written to: ", fig_dir, "\n", sep = "")

# Session info ----
sessionInfo()
