# Script name: 03_Genotype_stats.R
# ==============================================================================
# Title: Genotype frequency analysis in psychiatric and control cohorts.
#
# Author: Sergio Pérez Oliveira
#
# Revision purpose:
#   Reviewer-driven reanalysis of intermediate-allele (IA) frequencies with:
#     - explicit case/reference orientation for every odds ratio;
#     - exact Fisher tests for 2x2 comparisons;
#     - predefined multiple-testing families;
#     - one machine-readable master results table.
#
# IMPORTANT:
#   ORs are always reported as:
#       odds(IA in CASE) / odds(IA in REFERENCE)
#
#   Expanded alleles are excluded from IA-vs-normal comparisons.
#   They remain clinically important but are not mixed with the IA hypothesis family.
# ==============================================================================

# Load environment ----
Env_path <- file.choose()
source(Env_path)
rm(Env_path)

# Output directory ----
revision_dir <- file.path("results", "reviewer_revision")
if (!dir.exists(revision_dir)) {
  dir.create(revision_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# Helper: explicit case-vs-reference Fisher test
# ------------------------------------------------------------------------------

run_oriented_fisher <- function(data,
                                group_col,
                                genotype_col,
                                case,
                                reference,
                                gene,
                                test_id,
                                family_id,
                                adjustment_method,
                                analysis_label) {

  required_cols <- c(group_col, genotype_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "Missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  d <- data |>
    dplyr::transmute(
      group = as.character(.data[[group_col]]),
      genotype = as.character(.data[[genotype_col]])
    ) |>
    dplyr::filter(
      group %in% c(case, reference),
      genotype %in% c("NORMAL", "IA")
    )

  if (!all(c(case, reference) %in% unique(d$group))) {
    stop(
      "Test ", test_id, ": case/reference group missing after filtering.",
      call. = FALSE
    )
  }

  # Counts: rows = case/reference; columns = IA/NORMAL
  a <- sum(d$group == case      & d$genotype == "IA")
  b <- sum(d$group == case      & d$genotype == "NORMAL")
  c <- sum(d$group == reference & d$genotype == "IA")
  d0 <- sum(d$group == reference & d$genotype == "NORMAL")

  tab <- matrix(
    c(a, b, c, d0),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      group = c(case, reference),
      genotype = c("IA", "NORMAL")
    )
  )

  ft <- stats::fisher.test(tab, conf.int = TRUE)

  tibble::tibble(
    test_id = test_id,
    family_id = family_id,
    adjustment_method = adjustment_method,
    analysis = analysis_label,
    gene = gene,
    group_variable = group_col,
    genotype_variable = genotype_col,
    case = case,
    reference = reference,
    contrast = paste0(case, " vs ", reference),
    effect_definition = "OR = odds(IA in case) / odds(IA in reference)",
    n_case = a + b,
    n_reference = c + d0,
    n_total = a + b + c + d0,
    ia_case = a,
    normal_case = b,
    ia_reference = c,
    normal_reference = d0,
    ia_pct_case = 100 * a / (a + b),
    ia_pct_reference = 100 * c / (c + d0),
    odds_ratio = unname(ft$estimate),
    conf_low = unname(ft$conf.int[1]),
    conf_high = unname(ft$conf.int[2]),
    p_raw = ft$p.value
  )
}

# ------------------------------------------------------------------------------
# Multiplicity plan for the revision
#
# Family A:
#   IA frequency comparisons by diagnosis / predefined diagnostic subtype
#   across HTT, ATXN1 and ATXN2.
#   Adjustment: Holm (family-wise error control).
#
# Family C:
#   IA comparisons involving cognitive-status subgroups.
#   These are explicitly exploratory.
#   Adjustment: Benjamini-Hochberg FDR.
#
# This family structure is defined in code before the revision rerun and is not
# changed according to which individual p-values are significant.
# ------------------------------------------------------------------------------

genes <- tibble::tribble(
  ~gene,   ~genotype_col,
  "HTT",   "HTT_CODE",
  "ATXN1", "ATXN1_CODE",
  "ATXN2", "ATXN2_CODE"
)

cross_gene_comparisons <- function(comparisons) {
  tibble::as_tibble(
    merge(genes, comparisons, by = NULL, sort = FALSE)
  )
}

# Family A1: main diagnostic groups
plan_main <- cross_gene_comparisons(
  tibble::tribble(
    ~case, ~reference,
    "BD",  "CONTROL",
    "SCZ", "CONTROL",
    "BD",  "SCZ"
  )
) |>
  dplyr::mutate(
    dataset = "DT",
    group_col = "PATHOLOGY",
    family_id = "A_IA_DIAGNOSIS",
    adjustment_method = "holm",
    analysis_scope = "MAIN_DIAGNOSIS",
    analysis_label = "IA frequency: main diagnostic comparison",
    test_id = paste("IA", gene, analysis_scope, case, "vs", reference, sep = "__")
  )

# Family A2: BD subtype comparisons
plan_bd_subtype <- cross_gene_comparisons(
  tibble::tribble(
    ~case, ~reference,
    "BD-I",  "CONTROL",
    "Other", "CONTROL",
    "BD-I",  "Other"
  )
) |>
  dplyr::mutate(
    dataset = "BD_CONTROLS",
    group_col = "PATHOLOGY_TYPE_BINARY",
    family_id = "A_IA_DIAGNOSIS",
    adjustment_method = "holm",
    analysis_scope = "BD_SUBTYPE",
    analysis_label = "IA frequency: BD subtype comparison",
    test_id = paste("IA", gene, analysis_scope, case, "vs", reference, sep = "__")
  )

# Family A3: SCZ subtype comparisons
plan_scz_subtype <- cross_gene_comparisons(
  tibble::tribble(
    ~case, ~reference,
    "SCZ",   "CONTROL",
    "Other", "CONTROL",
    "SCZ",   "Other"
  )
) |>
  dplyr::mutate(
    dataset = "SCZ_CONTROLS",
    group_col = "PATHOLOGY_TYPE_BINARY",
    family_id = "A_IA_DIAGNOSIS",
    adjustment_method = "holm",
    analysis_scope = "SCZ_SUBTYPE",
    analysis_label = "IA frequency: SCZ subtype comparison",
    test_id = paste("IA", gene, analysis_scope, case, "vs", reference, sep = "__")
  )

# Family C1: cognitive-status comparisons in BD
plan_bd_cd <- cross_gene_comparisons(
  tibble::tribble(
    ~case, ~reference,
    "CD",    "No-CD",
    "CD",    "CONTROL",
    "No-CD", "CONTROL"
  )
) |>
  dplyr::mutate(
    dataset = "BD_CONTROLS",
    group_col = "CD_BINARY",
    family_id = "C_IA_COGNITIVE",
    adjustment_method = "BH",
    analysis_scope = "BD_COGNITIVE",
    analysis_label = "Exploratory IA frequency: cognitive status in BD",
    test_id = paste("IA", gene, analysis_scope, case, "vs", reference, sep = "__")
  )

# Family C2: cognitive-status comparisons in SCZ
plan_scz_cd <- cross_gene_comparisons(
  tibble::tribble(
    ~case, ~reference,
    "CD",    "No-CD",
    "CD",    "CONTROL",
    "No-CD", "CONTROL"
  )
) |>
  dplyr::mutate(
    dataset = "SCZ_CONTROLS",
    group_col = "CD_BINARY",
    family_id = "C_IA_COGNITIVE",
    adjustment_method = "BH",
    analysis_scope = "SCZ_COGNITIVE",
    analysis_label = "Exploratory IA frequency: cognitive status in SCZ",
    test_id = paste("IA", gene, analysis_scope, case, "vs", reference, sep = "__")
  )

test_plan <- dplyr::bind_rows(
  plan_main,
  plan_bd_subtype,
  plan_scz_subtype,
  plan_bd_cd,
  plan_scz_cd
)

# Save the testing plan itself for auditability
readr::write_csv(
  test_plan,
  file.path(revision_dir, "03_IA_testing_plan.csv")
)

# ------------------------------------------------------------------------------
# Execute plan
# ------------------------------------------------------------------------------

dataset_lookup <- list(
  DT = DT,
  BD_CONTROLS = BD_CONTROLS,
  SCZ_CONTROLS = SCZ_CONTROLS
)

ia_results_raw <- purrr::pmap_dfr(
  test_plan,
  function(gene,
           genotype_col,
           case,
           reference,
           dataset,
           group_col,
           family_id,
           adjustment_method,
           analysis_scope,
           analysis_label,
           test_id) {

    run_oriented_fisher(
      data = dataset_lookup[[dataset]],
      group_col = group_col,
      genotype_col = genotype_col,
      case = case,
      reference = reference,
      gene = gene,
      test_id = test_id,
      family_id = family_id,
      adjustment_method = adjustment_method,
      analysis_label = analysis_label
    )
  }
)

# Apply multiplicity adjustment ONCE, after all raw tests in each family exist.
ia_results <- ia_results_raw |>
  dplyr::group_by(family_id) |>
  dplyr::mutate(
    p_adj = stats::p.adjust(
      p_raw,
      method = dplyr::first(adjustment_method)
    ),
    significant_raw = p_raw < 0.05,
    significant_adjusted = p_adj < 0.05
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(family_id, p_adj, p_raw)

# Save master table
readr::write_csv(
  ia_results,
  file.path(revision_dir, "03_IA_frequency_master_results.csv")
)

# ------------------------------------------------------------------------------
# Key reviewer-facing checks
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("KEY RESULT: HTT IA, BD-I vs CONTROL\n")
cat("OR orientation: BD-I (case) vs CONTROL (reference)\n")
cat("============================================================\n")

key_htt_bdi <- ia_results |>
  dplyr::filter(
    gene == "HTT",
    analysis == "IA frequency: BD subtype comparison",
    case == "BD-I",
    reference == "CONTROL"
  )

print(key_htt_bdi)

cat("\n============================================================\n")
cat("ATXN1 IA, BD vs SCZ\n")
cat("============================================================\n")

key_atxn1_bd_scz <- ia_results |>
  dplyr::filter(
    gene == "ATXN1",
    analysis == "IA frequency: main diagnostic comparison",
    case == "BD",
    reference == "SCZ"
  )

print(key_atxn1_bd_scz)

cat("\n============================================================\n")
cat("ALL ADJUSTED-SIGNIFICANT IA TESTS\n")
cat("============================================================\n")

print(
  ia_results |>
    dplyr::filter(significant_adjusted) |>
    dplyr::select(
      test_id, family_id, gene, contrast,
      n_total, odds_ratio, conf_low, conf_high,
      p_raw, p_adj
    )
)

# Session info ----
sessionInfo()
