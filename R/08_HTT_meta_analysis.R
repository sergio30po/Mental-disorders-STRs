# Script name: 08_HTT_meta_analysis.R
# ==============================================================================
# Title: Harmonized meta-analysis of HTT intermediate alleles in bipolar disorder.
#
# Reviewer-driven purpose:
#   External synthesis of the present cohort with Ferrari et al. (2023),
#   harmonizing the HTT intermediate-allele (IA) definition and odds-ratio
#   direction.
#
# IMPORTANT INTERPRETATION:
#   - This meta-analysis addresses ALL bipolar disorder (BD) vs controls because
#     Ferrari et al. report IA carrier counts only for their overall BD cohort.
#   - It is NOT a direct meta-analysis of the present study's BD-I subgroup.
#   - Therefore it cannot be described as an independent replication of the
#     BD-I-specific finding.
#   - It is a reviewer-requested external synthesis and must not be used to
#     override the prespecified within-study multiplicity correction.
#
# Harmonized definition:
#   HTT IA = 27-35 CAG repeats
#   HTT normal = <27 CAG repeats
#   Expanded alleles (>35) are excluded from this carrier-level comparison.
#
# External study:
#   Ferrari C, Capacci E, Bagnoli S, Ingannato A, Sorbi S, Nacmias B.
#   Genes (Basel). 2023;14(9):1681. doi:10.3390/genes14091681
#   Published counts: BD IA 7 / 69; controls IA 6 / 104.
#
# Note on Ferrari et al.:
#   The article reports p < 0.001 for 7/69 vs 6/104. This script deliberately
#   re-estimates the association from the published 2x2 counts rather than
#   propagating that reported p-value.
# ==============================================================================

# Load environment --------------------------------------------------------------
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
# 1. PRESENT-COHORT COUNTS: ALL BD VS CONTROL
# ==============================================================================

required_cols <- c("PATHOLOGY", "HTT_CODE")
missing_cols <- setdiff(required_cols, names(DT))
if (length(missing_cols) > 0) {
  stop(
    "DT is missing required column(s): ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

present_counts_long <- DT %>%
  dplyr::transmute(
    group = as.character(PATHOLOGY),
    HTT_status = as.character(HTT_CODE)
  ) %>%
  dplyr::filter(
    group %in% c("BD", "CONTROL"),
    HTT_status %in% c("IA", "NORMAL")
  ) %>%
  dplyr::count(group, HTT_status, name = "n")

get_count <- function(group_value, status_value) {
  x <- present_counts_long %>%
    dplyr::filter(
      group == group_value,
      HTT_status == status_value
    ) %>%
    dplyr::pull(n)

  if (length(x) == 0L) {
    return(0L)
  }

  as.integer(x[[1]])
}

present_bd_ia <- get_count("BD", "IA")
present_bd_normal <- get_count("BD", "NORMAL")
present_ctrl_ia <- get_count("CONTROL", "IA")
present_ctrl_normal <- get_count("CONTROL", "NORMAL")

if (
  present_bd_ia == 0L ||
  present_bd_normal == 0L ||
  present_ctrl_ia == 0L ||
  present_ctrl_normal == 0L
) {
  stop(
    "At least one required present-cohort 2x2 cell is zero or missing. ",
    "Inspect PATHOLOGY and HTT_CODE coding before meta-analysis.",
    call. = FALSE
  )
}

present_counts <- tibble::tibble(
  study = "Present cohort",
  case_definition = "All BD",
  control_definition = "Psychiatric/neurodegenerative disease-free controls",
  ancestry_context = "European ancestry; Spain",
  ia_definition = "27-35 CAG",
  case_IA = present_bd_ia,
  case_normal = present_bd_normal,
  control_IA = present_ctrl_ia,
  control_normal = present_ctrl_normal
)

# ==============================================================================
# 2. FERRARI ET AL. 2023 PUBLISHED COUNTS
# ==============================================================================

ferrari_counts <- tibble::tibble(
  study = "Ferrari et al. 2023",
  case_definition = "BD type I or II",
  control_definition = "Healthy controls from local Italian reference cohort",
  ancestry_context = "Italian cohort",
  ia_definition = "27-35 CAG",
  case_IA = 7L,
  case_normal = 62L,
  control_IA = 6L,
  control_normal = 98L
)

meta_counts <- dplyr::bind_rows(
  present_counts,
  ferrari_counts
)

readr::write_csv(
  meta_counts,
  file.path(revision_dir, "08_HTT_meta_harmonized_counts.csv")
)

# ==============================================================================
# 3. STUDY-SPECIFIC EFFECTS
# ==============================================================================

estimate_study <- function(study,
                           case_definition,
                           control_definition,
                           ancestry_context,
                           ia_definition,
                           case_IA,
                           case_normal,
                           control_IA,
                           control_normal) {

  # OR > 1 means higher odds of being a BD case among HTT-IA carriers
  # compared with HTT-normal individuals.
  log_or <- log(
    (case_IA * control_normal) /
      (case_normal * control_IA)
  )

  se_log_or <- sqrt(
    1 / case_IA +
      1 / case_normal +
      1 / control_IA +
      1 / control_normal
  )

  fisher <- stats::fisher.test(
    matrix(
      c(
        case_IA, case_normal,
        control_IA, control_normal
      ),
      nrow = 2,
      byrow = TRUE
    ),
    alternative = "two.sided"
  )

  tibble::tibble(
    study = study,
    case_definition = case_definition,
    control_definition = control_definition,
    ancestry_context = ancestry_context,
    ia_definition = ia_definition,
    case_IA = case_IA,
    case_normal = case_normal,
    control_IA = control_IA,
    control_normal = control_normal,
    case_IA_pct = 100 * case_IA / (case_IA + case_normal),
    control_IA_pct = 100 * control_IA / (control_IA + control_normal),
    log_OR = log_or,
    SE_log_OR = se_log_or,
    OR = exp(log_or),
    CI_low = exp(log_or - 1.96 * se_log_or),
    CI_high = exp(log_or + 1.96 * se_log_or),
    p_Wald = 2 * stats::pnorm(-abs(log_or / se_log_or)),
    p_Fisher = fisher$p.value,
    effect_definition = paste(
      "OR = odds(BD vs CONTROL) for HTT IA (27-35) vs NORMAL (<27)"
    )
  )
}

study_effects <- purrr::pmap_dfr(
  meta_counts,
  estimate_study
)

readr::write_csv(
  study_effects,
  file.path(revision_dir, "08_HTT_meta_study_effects.csv")
)

# ==============================================================================
# 4. FIXED-EFFECT AND DERSIMONIAN-LAIRD RANDOM-EFFECTS META-ANALYSIS
# ==============================================================================

y <- study_effects$log_OR
v <- study_effects$SE_log_OR^2
k <- length(y)

if (k < 2L) {
  stop("At least two studies are required.", call. = FALSE)
}

# Fixed effect
w_fixed <- 1 / v
mu_fixed <- sum(w_fixed * y) / sum(w_fixed)
se_fixed <- sqrt(1 / sum(w_fixed))

# Cochran Q
Q <- sum(w_fixed * (y - mu_fixed)^2)
df_Q <- k - 1L
p_Q <- stats::pchisq(Q, df = df_Q, lower.tail = FALSE)

# I^2
I2 <- if (Q > 0) {
  max(0, (Q - df_Q) / Q) * 100
} else {
  0
}

# DerSimonian-Laird tau^2
C_dl <- sum(w_fixed) - sum(w_fixed^2) / sum(w_fixed)
tau2_dl <- if (C_dl > 0) {
  max(0, (Q - df_Q) / C_dl)
} else {
  0
}

# Random effects
w_random <- 1 / (v + tau2_dl)
mu_random <- sum(w_random * y) / sum(w_random)
se_random <- sqrt(1 / sum(w_random))

meta_summary <- tibble::tibble(
  model = c(
    "Fixed-effect inverse-variance",
    "Random-effects DerSimonian-Laird"
  ),
  k = k,
  pooled_log_OR = c(mu_fixed, mu_random),
  pooled_OR = exp(c(mu_fixed, mu_random)),
  CI_low = exp(
    c(
      mu_fixed - 1.96 * se_fixed,
      mu_random - 1.96 * se_random
    )
  ),
  CI_high = exp(
    c(
      mu_fixed + 1.96 * se_fixed,
      mu_random + 1.96 * se_random
    )
  ),
  p_value = c(
    2 * stats::pnorm(-abs(mu_fixed / se_fixed)),
    2 * stats::pnorm(-abs(mu_random / se_random))
  ),
  Q = Q,
  Q_df = df_Q,
  Q_p = p_Q,
  I2_percent = I2,
  tau2 = tau2_dl,
  effect_definition = paste(
    "OR = odds(BD vs CONTROL) for HTT IA (27-35) vs NORMAL (<27)"
  )
)

readr::write_csv(
  meta_summary,
  file.path(revision_dir, "08_HTT_meta_summary.csv")
)

# ==============================================================================
# 5. HARMONIZATION / ELIGIBILITY NOTES
# ==============================================================================

eligibility <- tibble::tribble(
  ~source, ~included_primary_meta, ~reason,
  "Present cohort",
  TRUE,
  paste(
    "Carrier-level counts available; all-BD phenotype used to match Ferrari;",
    "HTT IA definition 27-35 CAG."
  ),
  "Ferrari et al. 2023",
  TRUE,
  paste(
    "Carrier-level BD and healthy-control counts published;",
    "HTT IA definition 27-35 CAG."
  ),
  "Ramos et al. 2015",
  FALSE,
  paste(
    "Published main frequency tables are chromosome-level rather than",
    "carrier-level and are not directly harmonized with the present",
    "carrier-level endpoint."
  ),
  "Present cohort BD-I subgroup",
  FALSE,
  paste(
    "Ferrari et al. do not publish HTT-IA carrier counts stratified by BD-I;",
    "therefore a direct BD-I meta-analysis cannot be performed from published data."
  )
)

readr::write_csv(
  eligibility,
  file.path(revision_dir, "08_HTT_meta_eligibility_notes.csv")
)

# ==============================================================================
# 6. REVIEWER-FACING CONSOLE OUTPUT
# ==============================================================================

cat("\n============================================================\n")
cat("HTT IA META-ANALYSIS: HARMONIZED STUDY COUNTS\n")
cat("============================================================\n")
print(meta_counts, n = Inf, width = Inf)

cat("\n============================================================\n")
cat("HTT IA META-ANALYSIS: STUDY-SPECIFIC EFFECTS\n")
cat("============================================================\n")
print(
  study_effects %>%
    dplyr::select(
      study,
      case_IA, case_normal,
      control_IA, control_normal,
      case_IA_pct, control_IA_pct,
      OR, CI_low, CI_high,
      p_Fisher
    ),
  n = Inf,
  width = Inf
)

cat("\n============================================================\n")
cat("FERRARI 2023 INTERNAL CHECK\n")
cat("============================================================\n")

ferrari_check <- study_effects %>%
  dplyr::filter(study == "Ferrari et al. 2023") %>%
  dplyr::select(
    study,
    case_IA, case_normal,
    control_IA, control_normal,
    OR, CI_low, CI_high,
    p_Fisher
  )

print(ferrari_check, n = Inf, width = Inf)

cat(
  "\nPublished text reports p < 0.001 for 7/69 vs 6/104; ",
  "the Fisher exact p-value recalculated from those published counts is shown above.\n",
  sep = ""
)

cat("\n============================================================\n")
cat("HTT IA META-ANALYSIS: POOLED EFFECT AND HETEROGENEITY\n")
cat("============================================================\n")
print(meta_summary, n = Inf, width = Inf)

cat("\n============================================================\n")
cat("META-ANALYSIS SCOPE NOTE\n")
cat("============================================================\n")
cat(
  "Primary synthesis is ALL BD vs CONTROL, not BD-I vs CONTROL.\n",
  "Ferrari et al. do not report IA counts separately for BD-I, so the ",
  "BD-I-specific result cannot be directly meta-analyzed from published data.\n",
  "The pooled result must not be used to override the present cohort's ",
  "prespecified within-study multiplicity correction.\n",
  sep = ""
)

# ==============================================================================
# 7. FOREST PLOT
# ==============================================================================

forest_studies <- study_effects %>%
  dplyr::transmute(
    label = study,
    OR = OR,
    CI_low = CI_low,
    CI_high = CI_high,
    type = "Study"
  )

forest_pooled <- meta_summary %>%
  dplyr::filter(model == "Random-effects DerSimonian-Laird") %>%
  dplyr::transmute(
    label = "Pooled (random effects)",
    OR = pooled_OR,
    CI_low = CI_low,
    CI_high = CI_high,
    type = "Pooled"
  )

forest_df <- dplyr::bind_rows(
  forest_studies,
  forest_pooled
) %>%
  dplyr::mutate(
    label = factor(
      label,
      levels = rev(label)
    )
  )

p_forest <- ggplot2::ggplot(
  forest_df,
  ggplot2::aes(
    x = OR,
    y = label
  )
) +
  ggplot2::geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(
      xmin = CI_low,
      xmax = CI_high
    ),
    height = 0.15
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      shape = type
    ),
    size = 3
  ) +
  ggplot2::scale_x_log10() +
  ggplot2::labs(
    x = "Odds ratio (log scale): HTT IA vs normal",
    y = NULL,
    shape = NULL
  ) +
  ggplot2::theme_classic(base_size = 12)

ggplot2::ggsave(
  filename = file.path(
    fig_dir,
    "08_HTT_IA_meta_analysis_forest.tiff"
  ),
  plot = p_forest,
  device = "tiff",
  width = 180,
  height = 110,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

cat("\n============================================================\n")
cat("08 COMPLETE\n")
cat("============================================================\n")
cat(
  "Harmonized all-BD HTT-IA meta-analysis completed.\n",
  "Study estimates, pooled estimates, heterogeneity metrics, eligibility notes, ",
  "and forest plot have been exported.\n",
  sep = ""
)

# Session info ------------------------------------------------------------------
sessionInfo()
