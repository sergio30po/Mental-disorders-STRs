# Script name: 06_Survival_age_analysis.R
# ==============================================================================
# Title: Clinical and survival analysis by genotype in psychiatric cohorts.

# Author: Sergio Pérez Oliveira

# Description: This script evaluates the association between STR-based genotype 
#              classifications (HTT_CODE, ATXN1_CODE, ATXN2_CODE) and clinical 
#              variables (age at onset, disease duration) in Schizophrenia (SCZ) 
#              and bipolar disorder (BD) patients. Analyses include subgroup 
#              comparisons within BD (e.g., BD-I, BD-CD) and survival models.

# Inputs:
#   - Manually selected environment file with custom functions (.R)
#   - Manually selected data file with clinical and genotype variables
#   - Dataframes: DT, BD, SCZ, BD_I, BD_CD, BD_noCD

# Outputs:
#   - Descriptive statistics for age at onset and disease duration (mean ± SD)
#   - Non-parametric tests (Wilcoxon, Kruskal-Wallis, Dunn with Holm correction)
#   - Subgroup comparisons within BD and SCZ by genotype category
#   - Survival objects using Surv() from survival package
#   - Cox proportional hazards models (with covariates: SEX, SMOKER, genotype)
#   - Log-rank tests and Kaplan-Meier plots (e.g., for ATXN2_CODE in BD)
#   - Supplementary figure 2
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

#CORRELATION TESTS ----
cor.test(BD$ONSET_AGE,BD$DURATION,method = "spearman")
cor.test(SCZ$ONSET_AGE,SCZ$DURATION,method = "spearman")

# Sup. Fig. 2A: Correlation plots (BD and SCZ) ----
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
# 1. AGE AT ONSET - DISEASE DURATION CORRELATIONS ----

# BD - HTT
mean_sd(BD, "HTT_CODE", "ONSET_AGE", "Age at onset BD")
kruskal.test(ONSET_AGE ~ HTT_CODE, data = BD)
dunn <- dunnTest(BD$ONSET_AGE ~ BD$HTT_CODE, method = "holm")
print(dunn, dunn.test.results = TRUE)

mean_sd(BD_I, "HTT_CODE", "ONSET_AGE", "Age at onset BD-I")
kruskal.test(ONSET_AGE ~ HTT_CODE, data = BD_I)
dunn <- dunnTest(BD_I$ONSET_AGE ~ BD_I$HTT_CODE, method = "holm")
print(dunn, dunn.test.results = TRUE)

mean_sd(BD_CD, "HTT_CODE", "ONSET_AGE", "Age at onset BD-CD")
kruskal.test(ONSET_AGE ~ HTT_CODE, data = BD_CD)

mean_sd(BD_NOCD, "HTT_CODE", "ONSET_AGE", "Age at onset BD-NO-CD")
kruskal.test(ONSET_AGE ~ HTT_CODE, data = BD_NOCD)

# SCZ - HTT
mean_sd(SCZ, "HTT_CODE", "ONSET_AGE", "Age at onset SCZ")
wilcox.test(ONSET_AGE ~ HTT_CODE, data = SCZ)

mean_sd(SCZ_P, "HTT_CODE", "ONSET_AGE", "Age at onset SCZ_P")
kruskal.test(ONSET_AGE ~ HTT_CODE, data = SCZ_P)

mean_sd(SCZ_CD, "HTT_CODE", "ONSET_AGE", "Age at onset SCZ_CD")
kruskal.test(ONSET_AGE ~ HTT_CODE, data = SCZ_CD)

mean_sd(SCZ_NOCD, "HTT_CODE", "ONSET_AGE", "Age at onset SCZ_NOCD")
kruskal.test(ONSET_AGE ~ HTT_CODE, data = SCZ_NOCD)

# BD - ATXN1
mean_sd(BD, "ATXN1_CODE", "ONSET_AGE", "Age at onset BD")
wilcox.test(ONSET_AGE ~ ATXN1_CODE, data = BD)

mean_sd(BD_I, "ATXN1_CODE", "ONSET_AGE", "Age at onset BD-I")
kruskal.test(ONSET_AGE ~ ATXN1_CODE, data = BD_I)

mean_sd(BD_CD, "ATXN1_CODE", "ONSET_AGE", "Age at onset BD-CD")
kruskal.test(ONSET_AGE ~ ATXN1_CODE, data = BD_CD)

mean_sd(BD_NOCD, "ATXN1_CODE", "ONSET_AGE", "Age at onset BD-NO-CD")
kruskal.test(ONSET_AGE ~ ATXN1_CODE, data = BD_NOCD)

# SCZ - ATXN1
mean_sd(SCZ, "ATXN1_CODE", "ONSET_AGE", "Age at onset SCZ")
wilcox.test(ONSET_AGE ~ ATXN1_CODE, data = SCZ)

mean_sd(SCZ_P, "ATXN1_CODE", "ONSET_AGE", "Age at onset SCZ_P")
kruskal.test(ONSET_AGE ~ ATXN1_CODE, data = SCZ_P)

mean_sd(SCZ_CD, "ATXN1_CODE", "ONSET_AGE", "Age at onset SCZ_CD")
kruskal.test(ONSET_AGE ~ ATXN1_CODE, data = SCZ_CD)

mean_sd(SCZ_NOCD, "ATXN1_CODE", "ONSET_AGE", "Age at onset SCZ_NOCD")
kruskal.test(ONSET_AGE ~ ATXN1_CODE, data = SCZ_NOCD)

# BD - ATXN2
mean_sd(BD, "ATXN2_CODE", "ONSET_AGE", "Age at onset BD")
wilcox.test(ONSET_AGE ~ ATXN2_CODE, data = BD)

mean_sd(BD_I, "ATXN2_CODE", "ONSET_AGE", "Age at onset BD-I")
wilcox.test(ONSET_AGE ~ ATXN2_CODE, data = BD_I)

mean_sd(BD_CD, "ATXN2_CODE", "ONSET_AGE", "Age at onset BD-CD")
wilcox.test(ONSET_AGE ~ ATXN2_CODE, data = BD_CD)

mean_sd(BD_NOCD, "ATXN2_CODE", "ONSET_AGE", "Age at onset BD-NO-CD")
kruskal.test(ONSET_AGE ~ ATXN2_CODE, data = BD_NOCD)

# SCZ - ATXN2
mean_sd(SCZ, "ATXN2_CODE", "ONSET_AGE", "Age at onset SCZ")
kruskal.test(ONSET_AGE ~ ATXN2_CODE, data = SCZ)

mean_sd(SCZ_P, "ATXN2_CODE", "ONSET_AGE", "Age at onset SCZ_P")
kruskal.test(ONSET_AGE ~ ATXN2_CODE, data = SCZ_P)

mean_sd(SCZ_CD, "ATXN2_CODE", "ONSET_AGE", "Age at onset SCZ_CD")
kruskal.test(ONSET_AGE ~ ATXN2_CODE, data = SCZ_CD)

mean_sd(SCZ_NOCD, "ATXN2_CODE", "ONSET_AGE", "Age at onset SCZ_NOCD")
kruskal.test(ONSET_AGE ~ ATXN2_CODE, data = SCZ_NOCD)

#BD APOE

mean_sd(BD, "APOE_E4", "ONSET_AGE", "Age at onset APOE E4")
wilcox.test(ONSET_AGE ~ APOE_E4, data = BD)
rank_biserial(ONSET_AGE ~ APOE_E4, data = BD)

mean_sd(BD_I, "APOE_E4", "ONSET_AGE", "Age at onset BD-I")
wilcox.test(ONSET_AGE ~ APOE_E4, data = BD_I)
rank_biserial(ONSET_AGE ~ APOE_E4, data = BD_I)

mean_sd(BD_II, "APOE_E4", "DURATION", "Duration BD_II")
wilcox.test(ONSET_AGE ~ APOE_E4, data = BD_II)
rank_biserial(ONSET_AGE ~ APOE_E4, data = BD_II)

mean_sd(BD_CD, "APOE_E4", "ONSET_AGE", "Age at onset BD-CD")
wilcox.test(ONSET_AGE ~ APOE_E4, data = BD_CD)
rank_biserial(ONSET_AGE ~ APOE_E4, data = BD_CD)

mean_sd(BD_NOCD, "APOE_E4", "ONSET_AGE", "Age at onset BD-NO-CD")
wilcox.test(ONSET_AGE ~ APOE_E4, data = BD_NOCD)
rank_biserial(ONSET_AGE ~ APOE_E4, data = BD_NOCD)

#SCZ APOE
mean_sd(SCZ, "APOE_E4", "ONSET_AGE", "Age at onset SCZ")
wilcox.test(ONSET_AGE ~ APOE_E4, data = SCZ)
rank_biserial(ONSET_AGE ~ APOE_E4, data = SCZ)

mean_sd(SCZ_P, "APOE_E4", "ONSET_AGE", "Age at onset SCZ_P")
wilcox.test(ONSET_AGE ~ APOE_E4, data = SCZ_P)
rank_biserial(ONSET_AGE ~ APOE_E4, data = SCZ_P)

mean_sd(SCZ_CD, "APOE_E4", "ONSET_AGE", "Age at onset SCZ_CD")
wilcox.test(ONSET_AGE ~ APOE_E4, data = SCZ_CD)
rank_biserial(ONSET_AGE ~ APOE_E4, data = SCZ_CD)

mean_sd(SCZ_NOCD, "APOE_E4", "ONSET_AGE", "Age at onset SCZ_NOCD")
wilcox.test(ONSET_AGE ~ APOE_E4, data = SCZ_NOCD)
rank_biserial(ONSET_AGE ~ APOE_E4, data = SCZ_NOCD)
# 2. AGE-OF-ONSET MODELING ----
# Linear  models within BD and within SCZ.
# Genes: HTT, ATXN1, ATXN2
#
# For EACH model:
#   - Full model: covariates + genetic block (linear + quadratic + interaction)
#   - Null model: covariates only (always retained)
#   - Global test of genetic block: nested-model ANOVA (F-test)
#   - AIC + delta AIC
#   - Optional stepwise backward selection (exploratory) restricted to genetics
#
# IMPORTANT:
#   - Do NOT use exponentiate=TRUE for lm() (no ORs here).
#   - Report p-values for genetic terms or for the global F-test (full vs null).

# - 2.1 Full model Age of onset and IA: fit + step only genes -----
fit_onset_step_genes <- function(df,
                                 outcome = "ONSET_AGE",
                                 covars = c("SEX","COFFEE","SMOKER","APOE_E4"),
                                 genes  = c("HTT_CODE","ATXN1_CODE","ATXN2_CODE")) {
  
  vars_needed <- c(outcome, covars, genes)
  
  d <- df %>%
    dplyr::select(all_of(vars_needed)) %>%
    dplyr::filter(if_all(everything(), ~ !is.na(.))) %>%
    dplyr::mutate(
      SEX      = droplevels(factor(SEX)),
      COFFEE   = droplevels(factor(COFFEE)),
      SMOKER   = droplevels(factor(SMOKER)),
      HTT_CODE = droplevels(factor(HTT_CODE)),
      ATXN1_CODE= droplevels(factor(ATXN1_CODE)),
      ATXN2_CODE= droplevels(factor(ATXN2_CODE))
    )
  
  f_lower <- as.formula(paste(outcome, "~", paste(covars, collapse = " + ")))
  f_upper <- as.formula(paste(outcome, "~", paste(c(covars, genes), collapse = " + ")))
  
  m_full  <- lm(f_upper, data = d)
  
  # Backward step restricted to genes: cannot drop covariates
  m_step <- step(
    object    = m_full,
    scope     = list(lower = f_lower, upper = f_upper),
    direction = "backward",
    trace     = 0
  )
  
  list(
    data = d,
    full = m_full,
    step = m_step,
    AIC  = AIC(m_full, m_step)
  )
}

# BD
res_BD  <- fit_onset_step_genes(df = BD)
summary(res_BD$step)
res_BD$AIC

res_BD_CD  <- fit_onset_step_genes(df = BD_CD)
summary(res_BD_CD$step)
res_BD_CD$AIC

# SCZ
res_SCZ <- fit_onset_step_genes(df = SCZ)
summary(res_SCZ$step)
res_SCZ$AIC

res_SCZ_CD  <- fit_onset_step_genes(df = SCZ_CD)
summary(res_SCZ_CD$step)
res_SCZ_CD$AIC

# - 2.2 Full model Age of onset and CAG repeats: fit + step only genes -----

drop_na <- function(df, vars) df %>% filter(if_all(all_of(vars), ~ !is.na(.)))

block_terms <- function(s, l) c(
  s, paste0("I(", s, "^2)"),
  l, paste0("I(", l, "^2)"),
  paste0(s, ":", l)
)

enforce_block_hierarchy <- function(keep, s, l) {
  qs <- paste0("I(", s, "^2)")
  ql <- paste0("I(", l, "^2)")
  it <- paste0(s, ":", l)
  
  keep <- unique(keep)
  if (qs %in% keep && !(s %in% keep)) keep <- c(keep, s)
  if (ql %in% keep && !(l %in% keep)) keep <- c(keep, l)
  if (it %in% keep) keep <- unique(c(keep, s, l))
  unique(keep)
}

fit_onset_multiblock <- function(df, outcome, covars, blocks_named, do_step = TRUE) {
  
  # factors for covars if present
  for (v in intersect(covars, names(df))) df[[v]] <- factor(df[[v]])
  
  # vars needed
  gene_vars <- unlist(lapply(blocks_named, \(b) c(b$short, b$long)))
  d <- drop_na(df, unique(c(outcome, covars, gene_vars)))
  
  # formulas
  block_rhs <- unlist(lapply(blocks_named, \(b) block_terms(b$short, b$long)))
  f_null <- as.formula(paste(outcome, "~", paste(covars, collapse = " + ")))
  f_full <- as.formula(paste(outcome, "~", paste(c(covars, block_rhs), collapse = " + ")))
  
  m_null <- lm(f_null, data = d)
  m_full <- lm(f_full, data = d)
  
  out <- list(
    data = d,
    null = m_null,
    full = m_full,
    an_full_vs_null = anova(m_null, m_full),
    step = NULL,
    an_step_vs_null = NULL
  )
  
  if (!isTRUE(do_step)) return(out)
  
  # step: only allow removing gene terms (covars fixed by scope lower)
  m_step_raw <- step(m_full, scope = list(lower = f_null, upper = f_full),
                     direction = "backward", trace = 0)
  sel <- attr(terms(m_step_raw), "term.labels")
  
  # keep covars always
  gene_all <- unique(block_rhs)
  gen_sel <- sel[sel %in% gene_all]
  
  # hierarchy per block
  gen_sel_h <- gen_sel
  for (b in blocks_named) gen_sel_h <- enforce_block_hierarchy(gen_sel_h, b$short, b$long)
  
  rhs <- c(covars, gen_sel_h)
  f_step <- as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))
  m_step <- lm(f_step, data = d)
  
  out$step <- m_step
  out$an_step_vs_null <- anova(m_null, m_step)
  out
}

# BD
blocks <- list(
  HTT  = list(short="ALLELE1_HTT",  long="ALLELE2_HTT"),
  ATXN1= list(short="ALLELE1_ATXN1", long="ALLELE2_ATXN1"),
  ATXN2= list(short="ALLELE1_ATXN2", long="ALLELE2_ATXN2")
)
covars_onset<-c("SEX","COFFEE","SMOKER","APOE_E4")

res_BD <- fit_onset_multiblock(BD, "ONSET_AGE", covars_onset, blocks, do_step = TRUE)

res_BD$an_full_vs_null
res_BD$an_step_vs_null
summary(res_BD$step)

#ATXN1 final model
f_atxn1 <- as.formula(paste(
  "ONSET_AGE ~", paste(covars_onset, collapse=" + "), "+",
  paste(block_terms("ALLELE1_ATXN1","ALLELE2_ATXN1"), collapse=" + ")
))
m_atxn1_BD <- lm(f_atxn1, data = drop_na(BD, c("ONSET_AGE", covars_onset, "ALLELE1_ATXN1","ALLELE2_ATXN1")))
anova(lm(as.formula(paste("ONSET_AGE ~", paste(covars_onset, collapse=" + "))), data = model.frame(m_atxn1_BD)),
      m_atxn1_BD)
summary(m_atxn1_BD)

# SCZ

res_SCZ <- fit_onset_multiblock(SCZ, "ONSET_AGE", covars_onset, blocks, do_step = TRUE)

res_SCZ$an_full_vs_null
res_SCZ$an_step_vs_null
summary(res_SCZ$step)

#ATXN1 final model
f_atxn1 <- as.formula(paste(
  "ONSET_AGE ~", paste(covars_onset, collapse=" + "), "+",
  paste(block_terms("ALLELE1_ATXN1","ALLELE2_ATXN1"), collapse=" + ")
))
m_atxn1_SCZ <- lm(f_atxn1, data = drop_na(SCZ, c("ONSET_AGE", covars_onset, "ALLELE1_ATXN1","ALLELE2_ATXN1")))
anova(lm(as.formula(paste("ONSET_AGE ~", paste(covars_onset, collapse=" + "))), data = model.frame(m_atxn1_SCZ)),
      m_atxn1_SCZ)
summary(m_atxn1_SCZ)

#Sup. Fig. 2B : Age of onset in BD vs ATXN1 ----
# Two clearly separated BD blues (light vs dark)
col_long <- c(
  "Normal (<33)"         = "#8CBDE6",  # light BD blue
  "Intermediate (33–38)" = "#163A5F"   # dark BD blue
)

# References (fixed values used for predictions)
ref_SEX    <- "Female"
ref_COFFEE <- "Coffee"
ref_SMOKER <- "Smoking"

# ---- Data prep 
BD_atxn1 <- BD %>%
  dplyr::select(
    ONSET_AGE, SEX, COFFEE, SMOKER,
    ALLELE1_ATXN1, ALLELE2_ATXN1
  ) %>%
  dplyr::filter(
    !is.na(ONSET_AGE),
    !is.na(SEX), !is.na(COFFEE), !is.na(SMOKER),
    !is.na(ALLELE1_ATXN1), !is.na(ALLELE2_ATXN1)
  ) %>%
  mutate(
    SEX    = factor(SEX),
    COFFEE = factor(COFFEE),
    SMOKER = factor(SMOKER),
    long_bin = case_when(
      ALLELE2_ATXN1 < 33 ~ "Normal (<33)",
      ALLELE2_ATXN1 >= 33 & ALLELE2_ATXN1 <= 38 ~ "Intermediate (33–38)",
      TRUE ~ NA_character_
    ),
    long_bin = factor(long_bin, levels = c("Normal (<33)", "Intermediate (33–38)"))
  ) %>%
  filter(!is.na(long_bin))

# ---- Model (as in your step-selected structure for ATXN1)
m_atxn1 <- lm(
  ONSET_AGE ~ SEX + COFFEE + SMOKER +
    ALLELE1_ATXN1 + I(ALLELE1_ATXN1^2) +
    ALLELE2_ATXN1 + I(ALLELE2_ATXN1^2) +
    ALLELE1_ATXN1:ALLELE2_ATXN1,
  data = BD_atxn1
)

# ---- Prediction grid over OBSERVED short-allele range
x_min_data <- min(BD_atxn1$ALLELE1_ATXN1, na.rm = TRUE)
x_max_data <- max(BD_atxn1$ALLELE1_ATXN1, na.rm = TRUE)
x_grid <- seq(x_min_data, x_max_data, by = 0.05)

# Representative long-allele value per bin (median within bin)
bin_reps <- BD_atxn1 %>%
  group_by(long_bin) %>%
  summarise(long_rep = median(ALLELE2_ATXN1, na.rm = TRUE), .groups = "drop")

# Baseline for centering (median alleles in the data)
a1_ref <- as.numeric(median(BD_atxn1$ALLELE1_ATXN1, na.rm = TRUE))
a2_ref <- as.numeric(median(BD_atxn1$ALLELE2_ATXN1, na.rm = TRUE))

# Build newdata for prediction
newdat <- expand.grid(
  ALLELE1_ATXN1 = x_grid,
  long_bin = levels(BD_atxn1$long_bin),
  stringsAsFactors = FALSE
) %>%
  left_join(bin_reps, by = "long_bin") %>%
  mutate(
    ALLELE2_ATXN1 = long_rep,
    SEX    = factor(ref_SEX,    levels = levels(BD_atxn1$SEX)),
    COFFEE = factor(ref_COFFEE, levels = levels(BD_atxn1$COFFEE)),
    SMOKER = factor(ref_SMOKER, levels = levels(BD_atxn1$SMOKER))
  )

# Predict + SE
pred <- predict(m_atxn1, newdata = newdat, se.fit = TRUE)
newdat$fit <- as.numeric(pred$fit)
newdat$se  <- as.numeric(pred$se.fit)

# Baseline prediction used to center y-axis
base_dat <- data.frame(
  ALLELE1_ATXN1 = a1_ref,
  ALLELE2_ATXN1 = a2_ref,
  SEX    = factor(ref_SEX,    levels = levels(BD_atxn1$SEX)),
  COFFEE = factor(ref_COFFEE, levels = levels(BD_atxn1$COFFEE)),
  SMOKER = factor(ref_SMOKER, levels = levels(BD_atxn1$SMOKER))
)
base_fit <- as.numeric(predict(m_atxn1, newdata = base_dat))

# Centered effects (delta) + centered CI
newdat <- newdat %>%
  mutate(
    delta = fit - base_fit,
    lo = (fit - 1.96 * se) - base_fit,
    hi = (fit + 1.96 * se) - base_fit
  )

# ---- Plot (single legend: fill drives legend; linetype legend removed)
lt_map <- c("Normal (<33)" = "dashed", "Intermediate (33–38)" = "solid")

p_atxn1 <- ggplot() +
  # CI ribbons (colored by long-allele bin)
  geom_ribbon(
    data = newdat,
    aes(
      x = ALLELE1_ATXN1, ymin = lo, ymax = hi,
      fill = long_bin, group = long_bin
    ),
    alpha = 0.25,
    color = NA
  ) +
  # Model lines (black; linetype by bin)
  geom_line(
    data = newdat,
    aes(
      x = ALLELE1_ATXN1, y = delta,
      linetype = long_bin, group = long_bin
    ),
    color = "black",
    linewidth = 0.9
  ) +
  # Points (colored by bin; black border)
  geom_point(
    data = BD_atxn1,
    aes(
      x = ALLELE1_ATXN1,
      y = ONSET_AGE - base_fit,
      fill = long_bin
    ),
    shape = 21,
    color = "black",
    stroke = 0.25,
    alpha = 0.45,
    size = 1.6,
    position = position_jitter(width = 0.10, height = 0)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.6) +
  labs(
    #subtitle = paste0(
    #  "Y-axis shows the change in model-predicted age at onset relative to an individual with median alleles (short=", a1_ref,
    #  ", long=", a2_ref, "), female, with coffee intake, and somker."
    #),
    x = expression(italic("ATXN1") * " short allele (CAG repeats)"),
    y = "Change in predicted age at onset (years)",
    fill = expression(italic("ATXN1") * " long allele group")
  ) +
  scale_fill_manual(values = col_long) +
  scale_linetype_manual(values = lt_map) +
  # One legend only (fill). Keep linetype mapping but hide its guide.
  guides(
    linetype = "none",
    fill = guide_legend(
      title = expression(italic("ATXN1") * " long allele group"),
      override.aes = list(
        linetype = c("dashed", "solid"),
        color = "black"
      )
    )
  ) +
  coord_cartesian(
    xlim = c(x_min_data, x_max_data),
    ylim = c(-30, 60),
    clip = "on"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = c(0.86, 0.82),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
    legend.key = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(6, 6, 6, 6)
  )

p_atxn1

ggsave(
  filename = file.path(fig_dir, "Sup_Fig_2B.tiff"),
  plot = p_atxn1,
  device = "tiff",
  width = 180, height = 120, units = "mm",
  dpi = 600, compression = "lzw"
)


# 3. DURATION CORRELATIONS ----

# BD - HTT
mean_sd(BD, "HTT_CODE", "DURATION", "Duration BD")
kruskal.test(DURATION ~ HTT_CODE, data = BD)

mean_sd(BD_I, "HTT_CODE", "DURATION", "Duration BD-I")
kruskal.test(DURATION ~ HTT_CODE, data = BD_I)

mean_sd(BD_CD, "HTT_CODE", "DURATION", "Duration BD-CD")
kruskal.test(DURATION ~ HTT_CODE, data = BD_CD)

mean_sd(BD_NOCD, "HTT_CODE", "DURATION", "Duration BD-NO-CD")
kruskal.test(DURATION ~ HTT_CODE, data = BD_NOCD)

# SCZ - HTT
mean_sd(SCZ, "HTT_CODE", "DURATION", "Duration SCZ")
wilcox.test(DURATION ~ HTT_CODE, data = SCZ)

mean_sd(SCZ_P, "HTT_CODE", "DURATION", "Duration SCZ_P")
kruskal.test(DURATION ~ HTT_CODE, data = SCZ_P)

mean_sd(SCZ_CD, "HTT_CODE", "DURATION", "Duration SCZ-CD")
kruskal.test(DURATION ~ HTT_CODE, data = SCZ_CD)

mean_sd(SCZ_NOCD, "HTT_CODE", "DURATION", "Duration SCZ-NO-CD")
kruskal.test(DURATION ~ HTT_CODE, data = SCZ_NOCD)

# BD - ATXN1
mean_sd(BD, "ATXN1_CODE", "DURATION", "Duration BD")
wilcox.test(DURATION ~ ATXN1_CODE, data = BD)

mean_sd(BD_I, "ATXN1_CODE", "DURATION", "Duration BD-I")
kruskal.test(DURATION ~ ATXN1_CODE, data = BD_I)

mean_sd(BD_CD, "ATXN1_CODE", "DURATION", "Duration BD-CD")
kruskal.test(DURATION ~ ATXN1_CODE, data = BD_CD)

mean_sd(BD_NOCD, "ATXN1_CODE", "DURATION", "Duration BD-NO-CD")
kruskal.test(DURATION ~ ATXN1_CODE, data = BD_NOCD)

# SCZ - ATXN1
mean_sd(SCZ, "ATXN1_CODE", "DURATION", "Duration SCZ")
wilcox.test(DURATION ~ ATXN1_CODE, data = SCZ)

mean_sd(SCZ_P, "ATXN1_CODE", "DURATION", "Duration SCZ_P")
kruskal.test(DURATION ~ ATXN1_CODE, data = SCZ_P)

mean_sd(SCZ_CD, "ATXN1_CODE", "DURATION", "Duration SCZ-CD")
kruskal.test(DURATION ~ ATXN1_CODE, data = SCZ_CD)

mean_sd(SCZ_NOCD, "ATXN1_CODE", "DURATION", "Duration SCZ-NO-CD")
kruskal.test(DURATION ~ ATXN1_CODE, data = SCZ_NOCD)

# BD - ATXN2
mean_sd(BD, "ATXN2_CODE", "DURATION", "Duration BD")
wilcox.test(DURATION ~ ATXN2_CODE, data = BD)

mean_sd(BD_I, "ATXN2_CODE", "DURATION", "Duration BD-I")
wilcox.test(DURATION ~ ATXN2_CODE, data = BD_I)

mean_sd(BD_CD, "ATXN2_CODE", "DURATION", "Duration BD-CD")
wilcox.test(DURATION ~ ATXN2_CODE, data = BD_CD)

mean_sd(BD_NOCD, "ATXN2_CODE", "DURATION", "Duration BD-NoCD")
wilcox.test(DURATION ~ ATXN2_CODE, data = BD_NOCD)

# SCZ - ATXN2
mean_sd(SCZ, "ATXN2_CODE", "DURATION", "Duration SCZ")
kruskal.test(DURATION ~ ATXN2_CODE, data = SCZ)

mean_sd(SCZ_P, "ATXN2_CODE", "DURATION", "Duration SCZ_P")
kruskal.test(DURATION ~ ATXN2_CODE, data = SCZ_P)

mean_sd(SCZ_CD, "ATXN2_CODE", "DURATION", "Duration SCZ-CD")
kruskal.test(DURATION ~ ATXN2_CODE, data = SCZ_CD)

mean_sd(SCZ_NOCD, "ATXN2_CODE", "DURATION", "Duration SCZ-NO-CD")

# BD - APOE
mean_sd(BD, "APOE_E4", "DURATION", "Duration BD")
wilcox.test(DURATION ~ APOE_E4, data = BD)
rank_biserial(DURATION ~ APOE_E4, data = BD)

mean_sd(BD_I, "APOE_E4", "DURATION", "Duration BD-I")
wilcox.test(DURATION ~ APOE_E4, data = BD_I)
rank_biserial(DURATION ~ APOE_E4, data = BD_I)

mean_sd(BD_II, "APOE_E4", "DURATION", "Duration BD_II")
wilcox.test(DURATION ~ APOE_E4, data = BD_II)
rank_biserial(DURATION ~ APOE_E4, data = BD_II)

mean_sd(BD_CD, "APOE_E4", "DURATION", "Duration BD-CD")
wilcox.test(DURATION ~ APOE_E4, data = BD_CD)
rank_biserial(DURATION ~ APOE_E4, data = BD_CD)

mean_sd(BD_NOCD, "APOE_E4", "DURATION", "Duration BD-NoCD")
wilcox.test(DURATION ~ APOE_E4, data = BD_NOCD)
rank_biserial(DURATION ~ APOE_E4, data = BD_NOCD)

# SCZ - APOE
mean_sd(SCZ, "APOE_E4", "DURATION", "Duration SCZ")
wilcox.test(DURATION ~ APOE_E4, data = SCZ)
rank_biserial(DURATION ~ APOE_E4, data = SCZ)

mean_sd(SCZ_P, "APOE_E4", "DURATION", "Duration SCZ_P")
wilcox.test(DURATION ~ APOE_E4, data = SCZ_P)
rank_biserial(DURATION ~ APOE_E4, data = SCZ_P)

mean_sd(SCZ_CD, "APOE_E4", "DURATION", "Duration SCZ-CD")
wilcox.test(DURATION ~ APOE_E4, data = SCZ_CD)
rank_biserial(DURATION ~ APOE_E4, data = SCZ_CD)

mean_sd(SCZ_NOCD, "APOE_E4", "DURATION", "Duration SCZ-NO-CD")
wilcox.test(DURATION ~ APOE_E4, data = SCZ_CD)
rank_biserial(DURATION ~ APOE_E4, data = SCZ_CD)

# 4. SURVIVAL CURVES ----

# BD
surv_object <- Surv(time = BD$DURATION, event = rep(1, length(BD$DURATION)))
cox_model <- coxph(surv_object ~ SEX + SMOKER + COFFEE + APOE_E4 + ATXN2_CODE + HTT_CODE + ATXN1_CODE, data = BD)
summary(cox_model)

survdiff(surv_object ~ COFFEE, data = BD, rho = 0)
survdiff(surv_object ~ APOE_E4, data = BD, rho = 0)
survdiff(surv_object ~ SEX, data = BD, rho = 0)
survdiff(surv_object ~ SMOKER, data = BD, rho = 0)
survdiff(surv_object ~ ATXN2_CODE, data = BD, rho = 0)
cox_model <- coxph(surv_object ~ ATXN2_CODE, data = BD)
summary(cox_model)
survdiff(surv_object ~ HTT_CODE, data = BD, rho = 0)
survdiff(surv_object ~ ATXN1_CODE, data = BD, rho = 0)

surv_object <- Surv(time = BD_NOCD$DURATION, event = rep(1, length(BD_NOCD$DURATION)))
survdiff(surv_object ~ APOE_E4, data = BD_NOCD, rho = 0)
survdiff(surv_object ~ ATXN2_CODE, data = BD_NOCD, rho = 0)
survdiff(surv_object ~ HTT_CODE, data = BD_NOCD, rho = 0)
survdiff(surv_object ~ ATXN1_CODE, data = BD_NOCD, rho = 0)

surv_object <- Surv(time = BD_CD$DURATION, event = rep(1, length(BD_CD$DURATION)))
cox_model <- coxph(surv_object ~ SEX + SMOKER + COFFEE + APOE_E4 + ATXN2_CODE + HTT_CODE + ATXN1_CODE, data = BD_CD)
summary(cox_model)
survdiff(surv_object ~ ATXN2_CODE, data = BD_CD, rho = 0)
cox_model <- coxph(surv_object ~ ATXN2_CODE, data = BD_CD)
summary(cox_model)
survdiff(surv_object ~ APOE_E4, data = BD_CD, rho = 0)
survdiff(surv_object ~ HTT_CODE, data = BD_CD, rho = 0)
survdiff(surv_object ~ ATXN1_CODE, data = BD_CD, rho = 0)

# SCZ
surv_object <- Surv(time = SCZ$DURATION, event = rep(1, length(SCZ$DURATION)))
cox_model <- coxph(surv_object ~ SEX + SMOKER + COFFEE + APOE_E4 + ATXN2_CODE + HTT_CODE + ATXN1_CODE, data = SCZ)
summary(cox_model)
survdiff(surv_object ~ APOE_E4, data = SCZ, rho = 0)
survdiff(surv_object ~ SEX, data = SCZ, rho = 0)
survdiff(surv_object ~ SMOKER, data = SCZ, rho = 0)
survdiff(surv_object ~ ATXN2_CODE, data = SCZ, rho = 0)
survdiff(surv_object ~ HTT_CODE, data = SCZ, rho = 0)
survdiff(surv_object ~ ATXN1_CODE, data = SCZ, rho = 0)

surv_object <- Surv(time = SCZ_NOCD$DURATION, event = rep(1, length(SCZ_NOCD$DURATION)))
survdiff(surv_object ~ APOE_E4, data = SCZ_NOCD, rho = 0)
survdiff(surv_object ~ HTT_CODE, data = SCZ_NOCD, rho = 0)
survdiff(surv_object ~ ATXN1_CODE, data = SCZ_NOCD, rho = 0)

surv_object <- Surv(time = SCZ_CD$DURATION, event = rep(1, length(SCZ_CD$DURATION)))
survdiff(surv_object ~ APOE_E4, data = SCZ_CD, rho = 0)
survdiff(surv_object ~ ATXN2_CODE, data = SCZ_CD, rho = 0)
survdiff(surv_object ~ HTT_CODE, data = SCZ_CD, rho = 0)
survdiff(surv_object ~ ATXN1_CODE, data = SCZ_CD, rho = 0)

# Sup. Fig. 2C: Forest plot from Cox model (BD) -----
# 1) Prepare BD dataset
BD_cox <- BD %>%
  filter(
    !is.na(DURATION),
    !is.na(SEX),
    !is.na(SMOKER),
    !is.na(COFFEE),
    !is.na(ATXN2_CODE),
    !is.na(HTT_CODE),
    !is.na(ATXN1_CODE),
    ATXN2_CODE != "EXPANDED",
    HTT_CODE  != "EXPANDED"
  ) %>%
  mutate(
    SEX       = droplevels(factor(SEX)),
    SMOKER    = droplevels(factor(SMOKER)),
    HTT_CODE  = droplevels(factor(HTT_CODE)),
    ATXN1_CODE = droplevels(factor(ATXN1_CODE)),
    ATXN2_CODE = droplevels(factor(ATXN2_CODE))
  )

# 2) Survival object
surv_object <- Surv(time = BD_cox$DURATION, event = rep(1, nrow(BD_cox)))

# 3) Cox model
cox_model <- coxph(
  surv_object ~ SEX + SMOKER + COFFEE + APOE_E4 + ATXN2_CODE + HTT_CODE + ATXN1_CODE,
  data = BD_cox
)

# 4) Tidy results
tt <- broom::tidy(cox_model, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term != "(Intercept)")

# 5) Forest plot data
x_min <- 0.25
x_max <- 4.5

col_nsig <- "#8CBDE6"  # blue
col_sig  <- "#d62728"  # red

df_fp <- tt %>%
  mutate(
    term = stringr::str_trim(term),  # por si hay espacios raros
    term_label = case_when(
      term == "SEX[T.Female]"          ~ "Female sex",
      term == "SMOKER[T.Non-smoking]"  ~ "Non-smoking",
      term == "COFFEE[T.Coffee]"       ~ "Coffee consumption",
      term == "HTT_CODE[T.IA]"         ~ "italic(HTT)~' IA'",
      term == "ATXN1_CODE[T.IA]"        ~ "italic(ATXN1)~' IA'",
      term == "ATXN2_CODE[T.IA]"        ~ "italic(ATXN2)~' IA'",
      TRUE ~ NA_character_
    ),
    sig05 = ifelse(p.value < 0.05, "Significant (p < 0.05)", "Not significant"),
    sig05 = factor(sig05, levels = c("Significant (p < 0.05)", "Not significant")),
    sig = -log10(pmax(p.value, 1e-300)),
    est_p = pmin(pmax(estimate,  x_min), x_max),
    lo_p  = pmin(pmax(conf.low,  x_min), x_max),
    hi_p  = pmin(pmax(conf.high, x_min), x_max),
    cut_left  = conf.low  < x_min,
    cut_right = conf.high > x_max
  ) %>%
  filter(!is.na(term_label))

order_terms <- c(
  "Female sex",
  "Non-smoking",
  "Coffee consumption",
  "italic(HTT)~' IA'",
  "italic(ATXN1)~' IA'",
  "italic(ATXN2)~' IA'"
)

df_fp <- df_fp %>%
  mutate(term_label = factor(term_label, levels = rev(order_terms)))

# 6) Forest plot
g_forest_D <- ggplot(df_fp, aes(x = est_p, y = term_label)) +
  coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.6, color = "grey45") +
  
  geom_errorbarh(
    aes(xmin = lo_p, xmax = hi_p, color = sig05),
    height = 0.18, linewidth = 0.9
  ) +
  geom_point(
    aes(size = sig, fill = sig05),
    shape = 21, color = "black", stroke = 0.35
  ) +
  
  geom_segment(
    data = df_fp %>% filter(cut_left),
    aes(x = x_min * 1.35, xend = x_min * 1.08,
        y = term_label, yend = term_label, color = sig05),
    inherit.aes = FALSE,
    arrow = arrow(type = "closed", length = unit(0.16, "cm")),
    linewidth = 0.9
  ) +
  geom_segment(
    data = df_fp %>% filter(cut_right),
    aes(x = x_max / 1.35, xend = x_max / 1.08,
        y = term_label, yend = term_label, color = sig05),
    inherit.aes = FALSE,
    arrow = arrow(type = "closed", length = unit(0.16, "cm")),
    linewidth = 0.9
  ) +
  
  scale_x_log10(name = "Hazard ratio (log scale)") +
  scale_color_manual(
    name = "Statistical significance",
    values = c(
      "Significant (p < 0.05)" = col_sig,
      "Not significant"        = col_nsig
    )
  ) +
  scale_fill_manual(
    name = "Statistical significance",
    values = c(
      "Significant (p < 0.05)" = col_sig,
      "Not significant"        = col_nsig
    )
  ) +
  scale_size_continuous(
    name = expression(-log[10](p)),
    range = c(2.6, 6.8)
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.box = "vertical",
    axis.title.y = element_blank(),
    plot.margin = margin(8, 16, 8, 8)
  ) +
  scale_y_discrete(labels = function(x) {
    out <- x
    is_math <- grepl("^italic\\(", x)
    out[is_math] <- sapply(x[is_math], function(z) as.expression(parse(text = z)))
    out
  })

g_forest_D

# 7) Save
ggsave(
  filename = file.path(fig_dir, "Sup_Fig_2C.tiff"),
  plot = g_forest_D,
  device = "tiff",
  width = 250, height = 160, units = "mm",
  dpi = 600, compression = "lzw"
)

# Sup. Fig. 2D: K-M BD ATXN2 -----
surv_object <- Surv(time = BD$DURATION, event = rep(1, length(BD$DURATION)))
gen.km <- survfit(surv_object ~ ATXN2_CODE, data = BD, type = "kaplan-meier", error = "tsiatis", conf.type = "log-log", conf.int = 0.95)
bd_palette <- c("#F4A6A6", "#C73A3A")

plot_km <- ggsurvplot(
  fit = gen.km,
  data = BD,
  conf.int = TRUE,
  pval = TRUE,
  risk.table = TRUE,
  palette = bd_palette,
  xlab = "Disease duration (years)",
  ylab = "Event-free probability",
  legend.labs = c("Normal", "IA"),
  surv.median.line = "none",
  risk.table.height = 0.25
)

leg_title <- expression(italic("ATXN2") * " genotype")

plot_km$plot <- plot_km$plot +
  scale_color_manual(
    values = bd_palette,
    labels = c("Normal", "IA"),
    name = leg_title
  ) +
  scale_fill_manual(
    values = bd_palette,
    labels = c("Normal", "IA"),
    name = leg_title
  )

plot_km$table <- plot_km$table +
  labs(y = leg_title)

print(plot_km)

tiff(
  filename = file.path(fig_dir, "Sup_Fig_2D.tiff"),
  width = 250, height = 160, units = "mm",
  res = 600, compression = "lzw"
)
print(plot_km)
dev.off()

# Sup. Fig. 2E: K-M BD-CD ATXN2 -----
surv_object <- Surv(time = BD_CD$DURATION, event = rep(1, length(BD_CD$DURATION)))
gen.km <- survfit(surv_object ~ ATXN2_CODE, data = BD_CD, type = "kaplan-meier", error = "tsiatis", conf.type = "log-log", conf.int = 0.95)
bd_palette <- c("#D65C5C", "#7A1F1F")

plot_km_CD <- ggsurvplot(
  fit = gen.km,
  data = BD_CD,
  conf.int = TRUE,
  pval = TRUE,
  risk.table = TRUE,
  palette = bd_palette,
  xlab = "Disease duration (years)",
  ylab = "Event-free probability",
  legend.labs = c("Normal", "IA"),
  surv.median.line = "none",
  risk.table.height = 0.25
)

leg_title <- expression(italic("ATXN2") * " genotype")

plot_km_CD$plot <- plot_km_CD$plot +
  scale_color_manual(
    values = bd_palette,
    labels = c("Normal", "IA"),
    name = leg_title
  ) +
  scale_fill_manual(
    values = bd_palette,
    labels = c("Normal", "IA"),
    name = leg_title
  )

plot_km_CD$table <- plot_km_CD$table +
  labs(y = leg_title)

print(plot_km_CD)

tiff(
  filename = file.path(fig_dir, "Sup_Fig_2E.tiff"),
  width = 250, height = 160, units = "mm",
  res = 600, compression = "lzw"
)
print(plot_km_CD)
dev.off()

# Build composite figure ----

# KM + risk table as ONE patchwork object
km_as_one_panel <- function(km_obj, curve_h = 3.0, table_h = 1.2) {
  
  km_comp <- km_obj$plot / km_obj$table +
    plot_layout(heights = c(curve_h, table_h))
  
  wrap_elements(full = km_comp)
}

# D (left) and E (right), in that order
km_D <- km_as_one_panel(plot_km,    curve_h = 3.0, table_h = 1.2)
km_E <- km_as_one_panel(plot_km_CD, curve_h = 3.0, table_h = 1.2)

Sup_Fig_2 <-
  panel_A /
  (p_atxn1 | g_forest_D) /
  (km_D | km_E) +
  plot_layout(heights = c(1.15, 1.25, 1.75)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 24),
    plot.tag.position = c(0, 1),
    plot.tag.padding = unit(4, "pt")
  )

Sup_Fig_2

ggsave(
  filename = file.path(fig_dir, "Supplementary Figure 2.tiff"),
  plot = Sup_Fig_2,
  device = "tiff",
  width = 500, height = 650, units = "mm",
  dpi = 600, compression = "lzw"
)

# Session info ----
sessionInfo()
